! <YieldModifications_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2018 met.no
!*
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!
! <YieldModifications_mod.f90 - part of the EMEP MSC-W Chemical transport Model>
!_____________________________________________________________________________!
! Added March 2017
! For placement of miscellaneous coefficients and routines which can be
! called from Solver_ml. Currently allows access to VBS and JPC-based yields
! D. Simpson & R. Bergström March-June 2017
!_____________________________________________________________________________!
!
 module YieldModifications_mod

  use CheckStop_ml,       only : StopAll
  use ChemFields_ml             ! => cell_tinv,  NSPEC_TOT, O3, NO2, etc.
  use ChemSpecs                  ! => NSPEC_TOT, O3, NO2, etc.
  use Config_module,  only : MasterProc, DebugCell, DEBUG, USES &
                            ,YieldModifications ! JPCIsoYield
  use NumberConstants,    only : UNDEF_R, UNDEF_I
  use SmallUtils_ml,      only : find_index, trims

  implicit none
  private

  public  :: doYieldModifications

  logical, public, save :: YieldModificationsInUse = .false.
  logical, private, save :: dbg

  ! VBS work ----------------------------------------------------------

  private :: init_VBSyields
  private :: update_VBSyields
  integer, private, parameter   :: NVBS = 5 ! Number VOC 
  real, private, save :: kro2ho2, kro2no    ! rate-coeffs and fraction
 ! Yield arrays used in CM_Reaction2 are hard coded, with OXY=1, ... see below
  real, public, dimension(0:3),  save :: &
    YCOXY =UNDEF_R,  YNOXY =UNDEF_R & !  OXYL carbon and non-carbon yields
   ,YCALK =UNDEF_R,  YNALK =UNDEF_R & !  ALK  carbon and non-carbon yields
   ,YCOLE =UNDEF_R,  YNOLE =UNDEF_R & !  OLE  carbon and non-carbon yields
   ,YCISOP=UNDEF_R,  YNISOP=UNDEF_R & !  ISOP carbon and non-carbon yields
   ,YCTERP=UNDEF_R,  YNTERP=UNDEF_R   !  TERP carbon and non-carbon yields

   type, private :: vbs_t
     character(len=4)   :: name
     real               :: mw
     real               :: omoc     ! OM/OC ratio or  C/non-C ratios
     real               :: ratio    ! eg  C/non-C ratios
     real               :: fHO2NO2  ! modifier to kHO2NO2 rate
     real, dimension(4) :: highnox, lownox
   end type vbs_t
   type(vbs_t), private, dimension(NVBS), save :: Yemep


  ! Coded up some JPC experiments for now
  ! Allows JPC-FY-acid or JPC-VY-neutral or similar

   private :: init_JPCyields
   private :: update_JPCyields

   real, public, save :: & ! CRUDE, but fix later
      YA0APINOH   &  ! Limit of 17% mass at low Isop
     ,YG_APINOH, YG_APINO3, YG_APINNO3, YA_APINOH, YA_APINO3, YA_APINNO3   &
     ,YG_BPINOH, YG_BPINO3, YG_BPINNO3, YA_BPINOH, YA_BPINO3, YA_BPINNO3   &
     ,YG_MTOH,   YG_MTO3,   YG_MTNO3,   YA_MTOH,   YA_MTO3,   YA_MTNO3  &
     ,YA_ISOPOH, YA_ISOPNO3 
   real, private, save  :: IsoOHyield = UNDEF_R  ! 1 or 4%  for JPC


 contains

  ! --------------------------------------------------------------------------
  !> SUBROUTINE doYieldModifications
  !! Called once on Solver first_call (which sets YieldModificationsInUse, and some
  !! constants yield values), at start of each set of chemical iterations
  !! for each grid-cell (to reset yields for this cell), then after each iteration 
  !! to update the yields based.

  subroutine doYieldModifications(txt)
     character(len=*), intent(in) :: txt

     dbg = ( DEBUG%RUNCHEM .and. DebugCell )

     if ( YieldModifications(1:3) == 'VBS' ) then
       YieldModificationsInUse = .true.

       call init_VBSyields(txt)

       if ( txt == 'lastFastChem' )  then ! VBS only used in slower CM_Reactions2
         call update_VBSyields()
       end if

     else if ( YieldModifications(1:3) == 'JPC' ) then
       YieldModificationsInUse = .true.
       call init_JPCyields(txt)
       if ( YieldModifications(1:6) == 'JPC-VY' .and. txt == 'run' )  then
         call update_JPCyields()  ! OH impact
       end if
     end if

  end subroutine doYieldModifications
  ! --------------------------------------------------------------------------
  !> SUBROUTINE init_JPCyields
  !! Converts VBS SOA yields to EMEP stiochiometry for ASOC (mw12) and non-C
  !! aerosol (MW=1).

  subroutine init_VBSyields(txt)
  !>--------------------------------------------------------------------------
   character(len=*), intent(in) :: txt
   character(len=*), parameter :: dtxt='initVBS:'
   logical, save :: first_call = .true.
   real :: mwC = 12.0, mwH = 1.0,  r1, r2
   integer :: isoa

  !! Yields derived from Tsimpidi et al., ACP, 2010, p529
  !! and Robert's VBS_SOAformation
   type(vbs_t), dimension(5), parameter :: vbs = [ &  ! vbs_t(&
      vbs_t('OXYL', 106.0, 2.1, -999., 0.859,  &     ! Tsimpidi ... ARO2
          [ 0.002, 0.195, 0.3,   0.435 ],  [ 0.075, 0.300, 0.375, 0.525 ])&
     ,vbs_t('C4H10',  58.0, 1.7, -999.,  0.625, &     ! Tsimpidi ... ALK4
          [ 0.000, 0.038, 0.0,   0.0   ],  [ 0.000, 0.075, 0.0  , 0.0   ])&
     ,vbs_t('C3H6' ,  42.0, 1.7, -999., 0.52,  &     ! Tsimpidi ... OLE1
          [ 0.001, 0.005, 0.038, 0.150 ],  [ 0.005, 0.009, 0.060, 0.225 ])&
     ,vbs_t('ISOP' ,  68.0, 2.0, -999., 0.706,  &     ! Tsimpidi ... ISOP
          [ 0.001, 0.023, 0.015, 0.000 ],  [ 0.009, 0.030, 0.015, 0.0   ])&
     ,vbs_t('APIN' , 136.0, 1.7, -999., 0.914,  &     ! Tsimpidi ... TERP
          [ 0.012, 0.122, 0.201, 0.500 ],  [ 0.107, 0.092, 0.359, 0.6   ])]


   if ( first_call ) then  ! 1st call for this cell and timestep

!-----------------------------------------------------------------------------------
! Conversions from 'literature' yields to EMEP, based on mail from Robert,
! 2017-06-02.
! The ratio between NON_C_[A,B]SOA_ugN and [A,B]SOC_ugN is set to get
! the OM/OC ratio that we want for the respective type of SOA, e.g. for
! aromatic SOA we have assumed OM/OC = 2.1, which means the mass of
! NON_C_ASOA should be 1.1 × the mass of ASOC (which means that the
! “molar” ratio should be 13.2, since we set M(NON_C)=1).

! The total NON_C + SOC mass yield in each VBS class is based on the
! yields given by Tsimpidi et al. – in the Aromatics high-NOx case we
! have the mass yields 0.002 VBS_ug1 + 0.195 VBS_ug10 + 0.3 VBS_ug1e2 +
! 0.435 VBS_ug1e3, which means that we have the following molar yields:

! Y(ASOC_ug1)       = 0.002 × M[OXYL] / [OM/OC](ASOA) / M[ASOC_ug1]
!                   = 0.002 × 106 / 2.1 / 12
!                   = 0.008413
!
! Y(NON_C_ASOA_ug1) = ( 0.002 × M[OXYL] - Y(ASOC_ug1) × M[ASOC_ug1] ) /
!                       M[NON_C_ASOA_ug1] 
!                   = 0.002 × M[OXYL] × ( 1 - 1 / [OM/OC](ASOA) ) / 
!                                            M[NON_C_ASOA_ug1] 
!                   = 0.002 × 106 (1-1/2.1)/1
!                   = 0.11105

! Y(ASOC_ug10)        = 0.195 × 106 / 2.1 / 12 = 0.820238
! Y(NON_C_ASOA_ug10)  = 0.195 × 106 (1-1/2.1)/1 = 10.8271
! Y(ASOC_ug1e2)       = 0.3 × 106 / 2.1 / 12 = 1.261905
! Y(NON_C_ASOA_ug1e2) = 0.3 × 106 (1-1/2.1)/1 = 16.6571
! Y(ASOC_ug1e3)       = 0.435 × 106 / 2.1 / 12 = 1.82976
! Y(NON_C_ASOA_ug1e3) = 0.435 × 106 (1-1/2.1)/1 = 24.15286
!
!DS the yield of NON_C can thus be determined by the yield of ASOC x some
! constants. Below, we have:
! non-carbon would be  vbs(isoa)%highnox(:) * vbs(isoa)%mw*r2/mwH  ! non-carbon
! non-carbon = Y(ASOC) *mwC/r1 * r2/mwH = 

      Yemep = vbs   ! initialise names, etc.

      do isoa = 1, NVBS 
         r1 = 1/vbs(isoa)%omoc       ! eg 1/2.1 OC/OM ! CHECK for other VOC
         r2 = (1-r1)                 ! for non-C

         Yemep(isoa)%highnox(:) =  vbs(isoa)%highnox(:) * vbs(isoa)%mw*r1/mwC  ! carbon
         Yemep(isoa)%lownox(:)  =  vbs(isoa)%lownox(:)  * vbs(isoa)%mw*r1/mwC  ! carbon
         Yemep(isoa)%ratio =  r2/mwH * mwC/r1  ! non-carbon to carbon ratio
         
         if ( MasterProc ) then
           write(*,'(a,4f10.4,2x,4f10.2)') Yemep(isoa)%name//&
             ',YIELD HighNox ASOC, NOC:: ', Yemep(isoa)%highnox, &
              Yemep(isoa)%highnox * Yemep(isoa)%ratio
           write(*,'(a,4f10.4,2x,4f10.2)') Yemep(isoa)%name//&
             ',YIELD Low Nox ASOC, NOC:: ', Yemep(isoa)%lownox, &
              Yemep(isoa)%lownox * Yemep(isoa)%ratio
         end if
      end do ! isoa

     first_call = .false.
   end if

  end subroutine init_VBSyields

  subroutine update_VBSyields()
     character(len=*), parameter :: dtxt='updateVBSY:'
     real          ::  kNO, kHO2, fNO
     integer ::isoa

     kro2no  = 2.54e-12*exp(360*cell_tinv)
     kro2ho2 = 2.91e-13*exp(1300*cell_tinv) ! will modify below

     kNO  = kro2no*xnew(NO)
     kHO2 = kro2ho2*xnew(HO2)

    ! Loop over bins.
     isoa=1  ! Oxyl
     fNO  = kNO/( kNO+kHO2*Yemep(isoa)%fHO2NO2 )
     YCOXY(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNOXY(:) = YCOXY(:) * Yemep(isoa)%ratio

     isoa=2  ! Alkanes
     fNO  = kNO/( kNO+kHO2*Yemep(isoa)%fHO2NO2 )
     YCALK(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNALK(:) = YCALK(:) * Yemep(isoa)%ratio

     isoa=3  ! Alkenes
     fNO  = kNO/( kNO+kHO2*Yemep(isoa)%fHO2NO2 )
     YCOLE(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNOLE(:) = YCOLE(:) * Yemep(isoa)%ratio

     isoa=4  ! Isop
     fNO  = kNO/( kNO+kHO2*Yemep(isoa)%fHO2NO2 )
     YCISOP(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNISOP(:) = YCISOP(:) * Yemep(isoa)%ratio

     isoa=5  ! terp
     fNO  = kNO/( kNO+kHO2*Yemep(isoa)%fHO2NO2 )
     YCTERP(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNTERP(:) = YCTERP(:) * Yemep(isoa)%ratio

      if ( dbg ) then
          write(*,'(a,f8.5,4es12.3)') 'YIELD RUN '//dtxt,  fNO, &
             kro2no, kro2ho2, xnew(NO), xnew(HO2)
          write(*,'(a,es10.2,4f10.4,2x,4f10.2)') Yemep(1)%name//&
             ',YIELD EMEP ASOC, NOC:: ', fNO, YCOXY(:), YNOXY(:)
      end if

  end subroutine update_VBSyields

  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  !> SUBROUTINE init_JPCyields
  !!
  subroutine init_JPCyields(txt)
  !>--------------------------------------------------------------------------
   character(len=*), intent(in) :: txt
   character(len=*), parameter :: dtxt='initJPC:'
   logical, save :: first_call = .true.

 ! WARNING hard-coded MWs here to allow compilation with EmChem. Will fix later
 ! For isoprene SOA product OM/OC ratio set to 2.0 (as in Bergstrom et al., 
 ! 2012 - based on Chhabra et al., 2010), which gives SOA MW=136, ISOP MW = 68
 ! For MT SOA we use MW=240

   real, parameter :: MWratioISOP = 68.0/136.0  ! = 0.5
   real, parameter :: MWratioMT = 136.0/240.0   ! = 0.5667


   if ( first_call ) then

     if( MasterProc ) write(*,*) dtxt//'YieldModifications:'//&
           trim(YieldModifications)

     if ( index( YieldModifications, 'acid') > 0 ) then
         IsoOHyield = 0.04                   ! default is acid, high yield ?
     else ! neutral
         IsoOHyield = 0.01 
     end if

   ! Isoprene

      YA_ISOPOH  = IsoOHyield * MWratioISOP

   ! Assume 4% (mass-based) BSOA yield for isoprene + NO3 reaction - based
   ! on lowest yield observed by Ng et al. 2008

      YA_ISOPNO3 = 0.04 * MWratioISOP

   ! Monoterpenes. Yields depend on dC5H8/dMT if JPC-VY

     ! OH: 
      YA0APINOH = 0.17 * MWratioMT !species(APINENE)%molwt/species(APINOH_BSOA_NV)%molwt

     ! O3, fixed yields
      YA_APINO3 = 0.15 * MWratioMT !species(APINENE)%molwt/species(APINO3_BSOA_NV)%molwt
      YG_APINO3 = 1 - YA_APINO3

      YA_BPINO3 = YA_APINO3
      YG_BPINO3 = YG_APINO3
      YA_MTO3   = YA_APINO3
      YG_MTO3   = 1- YA_APINO3

     ! NO3, fixed yields
      YA_BPINNO3 = 0.3  * MWratioMT
      YG_BPINNO3 = 1 - YA_BPINNO3

      YA_MTNO3   =  0.3 * MWratioMT
      YG_MTNO3   =  1- YA_APINNO3

      first_call = .false.

      if ( MasterProc ) write(*,*) dtxt//trim(txt), IsoOHYield, YA0APINOH
    end if

    ! Needed at start of each chem timestep (may have changed if JPC-VY used)
    YA_APINOH = YA0APINOH
    YG_APINOH = 1 - YA_APINOH

    YA_BPINOH = YA_APINOH
    YG_BPINOH = YG_APINOH
    YA_MTOH = YA_APINOH
    YG_MTOH = 1- YA_APINOH

    if ( dbg ) then
       write(*,"(3a,3es12.3)") dtxt//trim(txt), &
           DEBUG%datetxt, "====", YA_MTOH, YA_MTO3,YA_MTNO3
    end if
    
  end subroutine init_JPCyields

  ! --------------------------------------------------------------------------
  !> SUBROUTINE update_JPCyields
  !  Only needed for JPC VY runs, and only affects OH yields

   subroutine update_JPCyields()
     character(len=*), parameter :: dtxt='updateJPCY:'
     real, parameter :: Y1=0.48, YD = 1-Y1
     integer, save :: i_iso=UNDEF_I, i_mt=UNDEF_I ! indices of OHLOSS
     logical, save :: first_call = .true.
     real          ::  diso, dmt ! Note, these are kOH*Iso
     real          ::  ratio, yield

     if ( first_call ) then
        i_iso = find_index('OHLOSS_ISO',species(:)%name)
        i_mt  = find_index('OHLOSS_MT',species(:)%name)
        if(MasterProc) write(*,*) dtxt//' i   ', i_iso, i_mt
        if ( i_iso < 1 .or. i_mt < 1 ) then
           print *, dtxt//"MISSING OHLOSS for JPC: ", i_iso, i_mt
           call StopAll(dtxt//'JPC OHLOSS ERR')
       end if
       first_call = .false.
     end if

     ! nb 1.0e5 OH * 10 ppt C5H8 * rc  -> P(OHLOSS)~2500 

     diso=xnew(i_iso)  ! Remember, = kOH*Iso

     if ( diso < 1.0 ) return
      
     dmt =xnew(i_mt)
     ratio = diso /(1.0+diso + dmt )

     if ( dbg ) then
        write(*,"(4a,4es12.3)") dtxt, DEBUG%datetxt, "====", &
             dtxt, diso,dmt, ratio, xnew(C5H8)
     end if

     if ( ratio < 0.01 )  then    ! exact to 3 sig figs
         yield  = 1 - ratio
     else if ( ratio < 5 )  then    ! exp(-1.53*5) = 2.5e-4 anyway: y=0.4802
         !y =          0.48 + 0.52 * exp(-1.53 * ratio )
         yield  =        Y1   + YD * exp(-1.53 * ratio )
     else 
         yield = Y1
     end if

    YA_APINOH = YA0APINOH * yield
    YG_APINOH = 1 - YA_APINOH

    ! We assume same OH behaviour for other MT:
    YA_BPINOH = YA_APINOH
    YG_BPINOH = 1-YG_BPINOH

    YA_MTOH = YA_APINOH
    YG_MTOH = 1- YA_MTOH

    if ( dbg ) then
      write(*,"(2a,9es12.3)") dtxt//' diso dmt rat y ', DEBUG%datetxt, diso, dmt, ratio, yield
      if ( ratio > 0.01 .and. ratio < 0.5  ) then
            write(*,"(3a,9es12.3)") dtxt//"JPC5YIELD ",dtxt, DEBUG%datetxt, &
              diso, dmt, ratio, yield, YA_APINOH
      end if
    end if

  end subroutine update_JPCyields

 end module YieldModifications_mod
