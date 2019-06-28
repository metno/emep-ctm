! <YieldModifications_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2019 met.no
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
! NOTE!!! AUG2018 - from emep-dev, diff coeffs to esx/src version!
! Added March 2017
! For placement of miscellaneous coefficients and routines which can be
! called from Solver_mod. Currently allows access to VBS and JPC-based yields
! D. Simpson & R. Bergström March-June 2017
!_____________________________________________________________________________!
!
 module YieldModifications_mod

  use CheckStop_mod,   only : StopAll
  use ChemFields_mod             ! => cell_tinv,  NSPEC_TOT, O3, NO2, etc.
  use ChemSpecs_mod              ! => NSPEC_TOT, O3, NO2, etc.
  use Config_module,   only : MasterProc,YieldModifications ! JPCIsoYield
  use Debug_module,    only : DebugCell, DEBUG  !-> DEBUG%SOA
  use NumberConstants, only : UNDEF_R, UNDEF_I
  use SmallUtils_mod,  only : find_index, trims

  implicit none
  private

  public  :: doYieldModifications

  logical, public, save :: YieldModificationsInUse = .false.
  logical, private, save :: dbg

  ! VBS work ----------------------------------------------------------

  private :: init_VBSyields
  private :: update_VBSyields
  integer, private, parameter   :: NVBS = 9 !Number VOC 
  real, private, save :: kro2ho2, kro2no, fNO=UNDEF_R  !  rate-coeffs and fraction
 ! Yield arrays used in CM_Reaction2 are hard coded, with OXY=1, ... see below
  real, public, dimension(-2:3),  save :: & !  Hodzic VBS_new range
    YCOXY =UNDEF_R,  YNOXY =UNDEF_R & !  OXYL carbon and non-carbon yields
   ,YCALK =UNDEF_R,  YNALK =UNDEF_R & !  ALK  carbon and non-carbon yields
   ,YCOLE =UNDEF_R,  YNOLE =UNDEF_R & !  OLE  carbon and non-carbon yields
   ,YCISOP=UNDEF_R,  YNISOP=UNDEF_R & !  ISOP carbon and non-carbon yields
   ,YCTERP=UNDEF_R,  YNTERP=UNDEF_R & !  TERP carbon and non-carbon yields
   ,YCBENZ=UNDEF_R,  YNBENZ=UNDEF_R & !  Benzene carbon and non-carbon yields
   ,YCTOL =UNDEF_R,  YNTOL =UNDEF_R & !  Toluene carbon and non-carbon yields
   ,YCIVOC=UNDEF_R,  YNIVOC=UNDEF_R   !  IVOC carbon and non-carbon yields

   type, private :: vbs_t
     character(len=4)   :: name
     real               :: mw
     real               :: omoc     ! OM/OC ratio or  C/non-C ratios
     real               :: ratio    ! eg  C/non-C ratios
     real               :: fHO2RO2  ! modifier to kHO2RO2 rate
     real               :: kRO2     ! RO2 reaction rate
     real, dimension(-2:3) :: highnox, lownox
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
     logical, save :: my_first_call = .true.

     dbg = ( DEBUG%SOA>0 .and. DebugCell )

     if ( YieldModifications(1:3) == 'VBS' ) then
       if( my_first_call .and. txt=='init' ) then
         YieldModificationsInUse = .true.
         call init_VBSyields()
       end if
       
       call update_VBSyields()

     else if ( YieldModifications(1:3) == 'JPC' ) then
       if( my_first_call .and. txt=='init' ) then
          YieldModificationsInUse = .true.
          call init_JPCyields(txt)
       ! QUERY if update needed on 1st call... do LATER
       else if ( YieldModifications(1:6) == 'JPC-VY' .and. txt == 'run' )  then
         call update_JPCyields()  ! OH impact
       end if
     else
       RETURN ! no yield modifications
     end if

     my_first_call = .false.
  end subroutine doYieldModifications

  ! --------------------------------------------------------------------------
  !> SUBROUTINE init_VBSyields
  !! Converts VBS SOA yields to EMEP stiochiometry for ASOC (mw12) and non-C
  !! aerosol (MW=1).

  subroutine init_VBSyields()
  !>--------------------------------------------------------------------------
   character(len=*), parameter :: dtxt='initVBS:'
   real :: mwC = 12.0, mwH = 1.0,  r1, r2
   integer :: bin, lev, isoa

  !! IMPORTANT. Order must be OXYL, C4H10, C3H6, ISOP, APIN,
  !! C15ivoc, C13ivoc (C13 not used)
  !! Yields derived from Tsimpidi et al., ACP, 2010, p529
  !! and Robert's VBS_SOAformation
   type(vbs_t), dimension(9), parameter :: vbs_T10 = [ &  ! vbs_t(&
      vbs_t('OXYL', 106.0, 2.1, -999., 0.859, 9.2e-14,  &     ! Tsimpidi ... ARO2 -- note unusually low RO2 rate for OXYL
          ! 0.01, 0.1,     1     10    100   1000
          [ 0.0,  0.0, 0.002, 0.195, 0.3,   0.435 ],&  ! high
          [ 0.0,  0.0, 0.075, 0.300, 0.375, 0.525 ])&  ! low
     ,vbs_t('C4H10',  58.0, 1.7, -999., 0.625, 2.5e-13, &     ! Tsimpidi ... ALK4
          [ 0.0,  0.0, 0.000, 0.038, 0.0,   0.0   ],&
          [ 0.0,  0.0, 0.000, 0.075, 0.0  , 0.0   ])&
     ,vbs_t('C3H6' ,  42.0, 1.7, -999., 0.52, 8.8e-13,  &     ! Tsimpidi ... OLE1
          [ 0.0,  0.0, 0.001, 0.005, 0.038, 0.150 ],&
          [ 0.0,  0.0, 0.005, 0.009, 0.060, 0.225 ])&
     ,vbs_t('ISOP' ,  68.0, 2.0, -999., 0.706, 8.0e-13,  &     ! Tsimpidi ... ISOP, RO2-rate for ISOPBO2
          [ 0.0,  0.0, 0.001, 0.023, 0.015, 0.000 ],&
          [ 0.0,  0.0, 0.009, 0.030, 0.015, 0.0   ])&
     ,vbs_t('APIN' , 136.0, 1.7, -999., 0.914, 3.6e-13,  &     ! Tsimpidi ... TERP, RO2-rate from MCM (weighted average of three RO2 species, from a-pinene+OH)
          [ 0.0,  0.0, 0.012, 0.122, 0.201, 0.500 ],&
          [ 0.0,  0.0, 0.107, 0.092, 0.359, 0.6   ]) &
  ! Jathar et al 2010 used POA/POC = 1.4
     ,vbs_t('C15ivoc' , 252.0, 1.4, -999., 0.914, 0.e-13,  &   !JAthar 2010 CHANGE
          [ 0.0,  0.044, 0.071, 0.41,  0.30, 0.0 ],&   ! used same for h/l now
          [ 0.0,  0.044, 0.071, 0.41,  0.30, 0.0 ])&
     ,vbs_t('C13ivoc' , 218.4, 1.4, -999., 0.914, 0.e-13,  &   !JAthar 2010 CHANGE
          [ 0.0,  0.014, 0.059, 0.22,  0.40, 0.0 ],&   ! used same for h/l now
          [ 0.0,  0.014, 0.059, 0.22,  0.40, 0.0 ])&
     ,vbs_t('BENZ',  78., 2.1, -999., 0.77, 8.8e-13, &     ! Hodzic BENZ. ? Mw,fHO2..
          [ 0.031, 0.011, 0.507, 0.019, 0.030, 0.142 ],&  ! high
          [ 0.007, 0.003, 0.270, 0.142, 0.400, 0.120 ])&  ! low 
     ,vbs_t('TOLU',  92., 2.1, -999., 0.82, 8.8e-13, &     ! Hodzic TOL. ? Mw,fHO2..
          [ 0.042, 0.123, 0.263, 0.020, 0.319, 0.329 ],&  ! high
          [ 0.371, 0.028, 0.207, 0.586, 0.063, 0.138 ])]  ! low 
! Not used yet. NOTE - starts at C* 0.1 so would need -1 col
!    ! Alkanes - only have high NOx, from Presto, Ots. No RO2, sso skip fHO2RO2
!     ,vbs_t('C15H2' , 212.4, 1.7, -999., -999.,  &     ! Presto, 2010 Alkanes
!          [ 0.044, 0.071, 0.41 , 0.30  ],  [ 0.044, 0.071, 0.41 , 0.30  ])]
!          !    .1      1    10    100

  !! IMPORTANT. Order must be OXYL, C4H10, C3H6, ISOP, APIN,
  !! C15ivoc, C13ivoc (C13 not used)
  !! Yields derived from Hodzic et al., ACP, 2016 combined with Ma et al 2017
  !! For Ma 2017 to fill-in. Use High-NOx for both since that's all we have
  !! and mainly urban anyway
   type(vbs_t), dimension(9), parameter :: vbs_HM  = [ &  
      vbs_t('OXYL', 106.0, 2.1, -999., 0.859, 9.2e-14,  &     ! Hodzic XYL. ? Mw,fHO2.. -- note unusually low RO2 rate for OXYL
          ! 0.01,    0.1,     1     10    100   1000
          [ 0.015, 0.056, 0.006, 0.026, 0.087, 0.193 ],&  ! high
          [ 0.395, 0.041, 0.203, 0.121, 0.232, 0.145 ])&  ! low Loss>40%! 
     ,vbs_t('C4H10',  58.0, 1.7, -999., 0.625, 2.5e-13, &     ! Ma Alk5
          [ 0.0,  0.0, 0.157, 0.0, 0.0,   0.0   ],&
          [ 0.0,  0.0, 0.157, 0.0, 0.0,   0.0   ])&
     ,vbs_t('C3H6' ,  42.0, 1.7, -999., 0.52, 8.8e-13, &      ! Ma Ole2
          [ 0.0,  0.0, 0.052, 0.0, 0.183, 0.157 ],&
          [ 0.0,  0.0, 0.052, 0.0, 0.183, 0.157 ])&
     ,vbs_t('ISOP' ,  68.0, 2.0, -999., 0.706, 8.0e-13, &     ! RO2-rate for ISOPBO2
          [ 0.001,  0.001, 0.027, 0.021, 0.044, 0.185 ],&
          [ 0.012,  0.013, 0.001, 0.100, 0.078, 0.097 ])&
     ,vbs_t('APIN' , 136.0, 1.7, -999., 0.914, 3.6e-13, &     ! RO2-rate from MCM (weighted average of three RO2 species, from a-pinene+OH)
          [ 0.045,  0.015, 0.142, 0.061, 0.074, 0.165 ],&
          [ 0.093,  0.211, 0.064, 0.102, 0.110, 0.125 ])&
     ,vbs_t('C15ivoc', 189.0, 1.7, -999.,  1.0 , 2.5e-13, &  ! H16 IVOC+MW, f=1; using RO2-rate for DDECO2 in MCM
          [ 0.140,  0.136, 0.069, 0.019, 0.010, 0.012 ],&
          [ 0.315,  0.173, 0.046, 0.010, 0.007, 0.008 ])&
     ,vbs_t('C13ivoc', 189.0, 1.7, -999.,  1.0 , 2.5e-13, &  ! H16 IVOC+MW, f=1; using RO2-rate for DDECO2 in MCM
          [ 0.140,  0.136, 0.069, 0.019, 0.010, 0.012 ],&
          [ 0.315,  0.173, 0.046, 0.010, 0.007, 0.008 ])&
     ,vbs_t('BENZ',  78., 2.1, -999., 0.77, 8.8e-13, &     ! Hodzic BENZ. ? Mw,fHO2..
          [ 0.031, 0.011, 0.507, 0.019, 0.030, 0.142 ],&  ! high
          [ 0.007, 0.003, 0.270, 0.142, 0.400, 0.120 ])&  ! low 
     ,vbs_t('TOLU',  92., 2.1, -999., 0.82, 8.8e-13, &     ! Hodzic TOL. ? Mw,fHO2..
          [ 0.042, 0.123, 0.263, 0.020, 0.319, 0.329 ],&  ! high
          [ 0.371, 0.028, 0.207, 0.586, 0.063, 0.138 ])]  ! low 

 ! No Alkane or alkene in VBS_H16
 ! Could add BENZ, TOL, SESQ

   type(vbs_t), dimension(max(size(vbs_T10),size(vbs_HM))), save :: vbs

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

      select case(YieldModifications) ! initialise names, etc.
         case('VBS-T10')
            vbs   = vbs_T10
            if ( MasterProc ) then
               write(*,*) 'Using VBS-T10 yield modifications'
            endif
         case('VBS-HM')
            vbs   = vbs_HM
            if ( MasterProc ) then
               write(*,*) 'Using VBS-HM yield modifications' 
            endif
         case('VBS-TH')  ! Tsimpidi but Hodzic for IVOC
            vbs   = vbs_T10
            vbs(6) = vbs_HM(6)  ! IVOC
            if ( MasterProc ) then
               write(*,*) 'Using VBS-TH yield modifications' 
            endif
         case default
            call StopAll(dtxt//'Unknown YieldModifications '//YieldModifications)
      end select
      Yemep = vbs   ! copies across names etc., before we modify yields

      do isoa = 1, NVBS 
         r1 = 1/vbs(isoa)%omoc       ! eg 1/2.1 OC/OM ! CHECK for other VOC
         r2 = (1-r1)                 ! for non-C

         Yemep(isoa)%highnox(:) =  vbs(isoa)%highnox(:) * vbs(isoa)%mw*r1/mwC  ! carbon
         Yemep(isoa)%lownox(:)  =  vbs(isoa)%lownox(:)  * vbs(isoa)%mw*r1/mwC  ! carbon
         Yemep(isoa)%ratio =  r2/mwH * mwC/r1  ! non-carbon to carbon ratio
         
         if ( MasterProc ) then
           write(*,'(a,7f8.4,2x,7f8.4)') Yemep(isoa)%name//&
             ',YIELD HighNox ASOC, NOC:: ', Yemep(isoa)%highnox, &
              Yemep(isoa)%highnox * Yemep(isoa)%ratio
           write(*,'(a,7f8.4,2x,7f8.4)') Yemep(isoa)%name//&
             ',YIELD Low Nox ASOC, NOC:: ', Yemep(isoa)%lownox, &
              Yemep(isoa)%lownox * Yemep(isoa)%ratio
         end if
      end do ! isoa

  end subroutine init_VBSyields

  subroutine update_VBSyields()
     character(len=*), parameter :: dtxt='updateVBSY:'
     real          ::  fluxNO, fluxHO2, fluxRO2, fNO
     integer ::isoa

     kro2no  = 2.54e-12*exp(360*cell_tinv)
     kro2ho2 = 2.91e-13*exp(1300*cell_tinv) ! will modify below

     fluxNO  = kro2no*xnew(NO)
     fluxHO2 = kro2ho2*xnew(HO2)
!for later     fluxRO2 = xnew(RO2POOL)
     fluxRO2 = 0. ! initially neglect RO2+RO2-reactions

     fNO = 0.0   ! default, eg on 1st call
    ! Loop over bins.
     isoa=1  ! Xylenes
     if(xnew(NO)>1.0 ) fNO  = fluxNO/(fluxNO+fluxHO2*Yemep(isoa)%fHO2RO2+fluxRO2*Yemep(isoa)%kRO2 )
     YCOXY(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNOXY(:) = YCOXY(:) * Yemep(isoa)%ratio

     isoa=2  ! Alkanes
     if(xnew(NO)>1.0 ) fNO  = fluxNO/(fluxNO+fluxHO2*Yemep(isoa)%fHO2RO2+fluxRO2*Yemep(isoa)%kRO2 )
     YCALK(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNALK(:) = YCALK(:) * Yemep(isoa)%ratio

     isoa=3  ! Alkenes
     if(xnew(NO)>1.0 ) fNO  = fluxNO/(fluxNO+fluxHO2*Yemep(isoa)%fHO2RO2+fluxRO2*Yemep(isoa)%kRO2 )
     YCOLE(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNOLE(:) = YCOLE(:) * Yemep(isoa)%ratio

     isoa=4  ! Isop
     if(xnew(NO)>1.0 ) fNO  = fluxNO/(fluxNO+fluxHO2*Yemep(isoa)%fHO2RO2+fluxRO2*Yemep(isoa)%kRO2 )
     YCISOP(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNISOP(:) = YCISOP(:) * Yemep(isoa)%ratio

     isoa=5  ! terp
     if(xnew(NO)>1.0 ) fNO  = fluxNO/(fluxNO+fluxHO2*Yemep(isoa)%fHO2RO2+fluxRO2*Yemep(isoa)%kRO2 )
     YCTERP(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNTERP(:) = YCTERP(:) * Yemep(isoa)%ratio

     isoa=6  ! IVOC
     if(xnew(NO)>1.0 ) fNO  = fluxNO/(fluxNO+fluxHO2*Yemep(isoa)%fHO2RO2+fluxRO2*Yemep(isoa)%kRO2 )
     YCIVOC(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNIVOC(:) = YCIVOC(:) * Yemep(isoa)%ratio

     isoa=8 ! Benzene
     if(xnew(NO)>1.0 ) fNO  = fluxNO/(fluxNO+fluxHO2*Yemep(isoa)%fHO2RO2+fluxRO2*Yemep(isoa)%kRO2 )
     YCBENZ(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNBENZ(:) = YCBENZ(:) * Yemep(isoa)%ratio

     isoa=9 ! Toluene
     if(xnew(NO)>1.0 ) fNO  = fluxNO/(fluxNO+fluxHO2*Yemep(isoa)%fHO2RO2+fluxRO2*Yemep(isoa)%kRO2 )
     YCTOL(:) = fno * Yemep(isoa)%highnox(:) + (1-fNO) * Yemep(isoa)%lownox(:)
     YNTOL(:) = YCTOL(:) * Yemep(isoa)%ratio



      if ( dbg ) then
          !write(*,'(a,f8.5,4es12.3)') 'YIELD RUN '//dtxt,  fNO, &
          !   kro2no, kro2ho2, xnew(NO), xnew(HO2)
          write(*,'(a,es10.2,7f10.4,2x,7f10.2)') Yemep(1)%name//&
             ',YIELD@1ug/m3:: ', fNO, YCOXY(0), YNOXY(0)  ! 0= C*:1ug/m3
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
!TSTEMX program testr
!TSTEMX use Config_module, only :  YieldModifications
!TSTEMX use Config_module,  only : MasterProc
!TSTEMX use Debug_module,   only : DebugCell
!TSTEMX use YieldModifications_mod
!TSTEMX character(len=10), dimension(2) :: vstypes = [ 'VBS-T10', 'VBS-HM ' ]
!TSTEMX !CART DEBUG%RUNCHEM= .true.; DebugCell= .true.
!TSTEMX integer :: i
!TSTEMX !DebugCell= .true.
!TSTEMX do i = 1, 1 !CART 2
!TSTEMX  !CART YieldModifications =  vbstypes(i)
!TSTEMX  call doYieldModifications('first')
!TSTEMX  !! call doYieldModifications('lastFastChem' ) ! Would need xnew, T, etc.
!TSTEMX end do
!TSTEMX end program testr
