! <Setup_1d_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2015 met.no
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
module Setup_1d_ml
!-----------------------------------------------------------------------!
! DESCRIPTION
! Generates arrays for 1-D column , for input to chemical solver. The output
! fields are stored in the Setup_1dfields_ml module.
!-----------------------------------------------------------------------!

use AeroFunctions,       only: umWetRad,WetRad, pmSurfArea, cMolSpeed, UptakeRate
use AirEmis_ml,          only: airn, airlig   ! airborne NOx emissions
use Biogenics_ml,        only: SoilNOx
use Biogenics_ml,        only: EMIS_BioNat, EmisNat  
use Chemfields_ml,       only: xn_adv,xn_bgn,xn_shl, &
                               NSPEC_COL, NSPEC_BGN, xn_2d_bgn
use ChemFunctions_ml, only : S_RiemerN2O5
use ChemGroups_ml,       only: PM10_GROUP, PMFINE_GROUP, SIA_GROUP, SS_GROUP, DUST_GROUP
use CheckStop_ml,        only:  CheckStop, StopAll
use DerivedFields_ml,    only: d_2d, f_2d
use EmisDef_ml,          only:  gridrcemis, gridrcroadd, KEMISTOP,Emis_4D,N_Emis_4D,Found_Emis_4D
use EmisGet_ml,          only:  nrcemis, iqrc2itot  !DSRC added nrcemis
use Emissions_ml,        only:  SumSplitEmis
use ForestFire_ml,       only: Fire_rcemis, burning
use Functions_ml,        only:  Tpot_2_T
use ChemFields_ml,       only: SurfArea_um2cm3
use ChemSpecs  !,           only:  SO4,C5H8,NO,NO2,SO2,CO,
use ChemRates_rct_ml,    only:  set_rct_rates, rct
use GridValues_ml,       only:  xmd, GridArea_m2, & 
                                 debug_proc, debug_li, debug_lj,&
                                 A_mid,B_mid,gridwidth_m,dA,dB,&
                                 i_fdom, j_fdom
use Io_Progs_ml,         only: datewrite !MASS
use LocalVariables_ml,   only: Grid
use MassBudget_ml,       only: totem    ! sum of emissions
use MetFields_ml,        only: ps
use MetFields_ml,        only: roa, th, q, t2_nwp, cc3dmax, &
                               zen, Idirect, Idiffuse,z_bnd
use ModelConstants_ml,   only:  &
   DEBUG,DEBUG_MASS             &
  ,AERO                         & ! for wet radii and surf area.
  ,dt_advec                     & ! time-step
  ,IOU_INST                     & ! for OUTMISC
  ,MasterProc                   &
  ,PPB, PT                      & ! Pressure at top
  ,USES                         & ! Forest fires so far
  ,USE_SEASALT                  &
  ,USE_LIGHTNING_EMIS, USE_AIRCRAFT_EMIS      &
  ,USE_GLOBAL_SOILNOX, USE_DUST, USE_ROADDUST &
  ,VOLCANO_SR                   & ! Reduce Volcanic Emissions
  ,emis_inputlist               & ! Used in EEMEP
  ,KMAX_MID ,KMAX_BND, KCHEMTOP   ! Start and upper k for 1d fields
use My_Derived_ml,       only: EmisSplit_OUT
use Landuse_ml,          only: water_fraction, ice_landcover
use Par_ml,              only: me,MAXLIMAX,MAXLJMAX, & 
                               gi0,gi1,gj0,gj1,IRUNBEG,JRUNBEG
use PhysicalConstants_ml,only: ATWAIR, AVOG, PI, GRAV
use Radiation_ml,        only: PARfrac, Wm2_uE
use Setup_1dfields_ml,   only: &
   xn_2d                &  ! concentration terms
  ,rcemis, deltaZcm     &  ! emission terms and layer thickness
  ,rh, temp, tinv, itemp,pp      &  !
  ,amk, o2, n2, h2o     &  ! Air concentrations
  ,cN2O5, cHO2, cO3, cHNO3 &  ! mol speeds, m/s
  ,DpgNw,S_m2m3   &  ! for wet diameter and surf area
  ,aero_fom, aero_fss, aero_fdust
use SmallUtils_ml,       only: find_index
use Tabulations_ml,      only: tab_esat_Pa
use TimeDate_ml,         only: current_date, date
use ColumnSource_ml,     only: ColumnRate
use Units_ml,            only: to_number_cm3 ! converts roa [kg/m3] to M [molec/cm3]
!!  to_number_cm3=0.001*AVOG/ATWAIR,& ! from density (roa, kg/m3) to molecules/cm3
 

implicit none
private
!-----------------------------------------------------------------------!

public :: setup_1d   ! Extracts results for i,j column from 3-D fields
public :: setup_rcemis ! Emissions  (formerly "poll")
public :: reset_3d     ! Exports final results for i,j column to 3-D fields
                       ! (and XNCOL outputs if asked for)

! Indices for the species defined in this routine. Only set if found
! Hard-coded for 2 specs just now. Could extend and allocate.
integer, private, parameter :: NROADDUST = 2
integer, private, parameter :: iROADF=1,  iROADC=2
integer, private, save :: inat_RDF,  inat_RDC, inat_Rn222
integer, private, save :: itot_RDF=-999,  itot_RDC=-999, itot_Rn222=-999

! Minimum concentration allowed for advected species. Avoids some random
! variations between debugged/normal runs with concs of e.g. 1.0e-23 ppb or
! less. 1.0e-20 ppb is ca. 2.5e-10 molec/cm3 at sea-level, ca. 2.5e-11 at top
! e.g.: MINCONC = 1.0e-29

!DUST_ROAD_F

contains
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   subroutine setup_1d(i,j)

 !..   extracts data along one vertical column for input to chemical
 !     solver concentrations for chemistry......
 !

    integer, intent(in) :: i,j    ! coordinates of column
    character(len=9)  :: sub='setup_1d:'
    character(len=30)  :: fmt="(a,i3,99g13.4)"  ! default format
    logical :: debug_flag
    logical, save :: first_call = .true.
   ! for surface area calcs: 
    real :: ugtmp, ugsiaPM, ugDustF, ugSSaltF, ugDustC, ugSSaltC
    real :: ugSO4, ugNO3f, ugNO3c, ugRemF, ugRemC, ugpmF, ugpmC, rho
    logical :: is_finepm, is_ssalt, is_dust
    real, dimension(size(AERO%Inddry))  :: Ddry ! Dry diameter
    integer :: iw, ipm ! for wet rad
   ! if rates are wanted for d_2d output, we use these indices:
    integer, dimension(10), save :: d2index, id2rct 
    integer, save :: nd2d
    integer :: itmp

   ! local

    integer           :: k, n, ispec   ! loop variables
    real              :: qsat ! saturation water content

    debug_flag =  ( DEBUG%SETUP_1DCHEM .and. debug_proc .and.  &
      i==debug_li .and. j==debug_lj .and. current_date%seconds == 0 )
    if( first_call.and.debug_proc ) debug_flag = .true.  ! make sure we see 1st i,j combo

    if( debug_flag ) write(*,*) sub//"=DBG=======  ", first_call, me


    do k = KCHEMTOP, KMAX_MID

  !- to_number_cm3 - to scale from  density (roa, kg/m3) to  molecules/cm3
  ! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)

       amk(k) = roa(i,j,k,1) * to_number_cm3  ! molecules air/cm3

       h2o(k) = max( 1.e-5*amk(k), &
                     q(i,j,k,1)*amk(k)*ATWAIR/18.0)

      ! nb. max function for h2o  used as semi-lagrangian scheme used
      ! in LAM50 (and HIRLAM) often gives negative H2O....   :-(

       pp(k) = A_mid(k) + B_mid(k)*ps(i,j,1)

       temp(k) = th(i,j,k,1)* Tpot_2_T( pp(k) )

       itemp(k) = nint( temp(k) -1.E-9) ! the "-1.E-9" is put in order to
               ! avoid possible different roundings on different machines.

       qsat  = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
       rh(k) = min( q(i,j,k,1)/qsat , 1.0)
       rh(k) = max( rh(k) , 0.001)

        ! 1)/ Short-lived species - no need to scale with M

         do n = 1, NSPEC_SHL
               xn_2d(n,k) = max(0.0,xn_shl(n,i,j,k))
         end do ! ispec

        ! 2)/ Advected species
        do n = 1, NSPEC_ADV
              ispec = NSPEC_SHL + n
              xn_2d(ispec,k) = max(0.0,xn_adv(n,i,j,k)*amk(k))
        end do ! ispec

        ! 3)/ Background species ( * CTM2 with units in mix. ratio)
        do n = 1, NSPEC_BGN
              xn_2d_bgn(n,k) = max(0.0,xn_bgn(n,i,j,k)*amk(k))
        end do ! ispec

      ! Surf Area
        if ( USES%SURF_AREA ) then ! GERBER

           !ispec=NO3_c ! CRUDE HARD CODE for now, but NO3 is special

           ugSO4       = 0.0
           ugRemF      = 0.0
           ugRemC      = 0.0
           ugpmF       = 0.0
           ugpmC       = 0.0
           ugNO3f      = 0.0 !  0.27*xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG
           ugNO3c      = 0.0 !  0.27*xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG
           ugsiaPM     = 0.0 


           ugDustF  = 0.0
           ugSSaltF = 0.0
           ugDustC  = 0.0
           ugSSaltC = 0.0

           do ipm = 1, size( PM10_GROUP )
             ispec = PM10_GROUP(ipm)

             ugtmp  = xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG
             is_finepm = ( find_index( ispec, PMFINE_GROUP) > 0 )
             is_ssalt  = ( find_index( ispec, SS_GROUP ) >0)
             is_dust   = ( find_index( ispec, DUST_GROUP )>0)
             if( is_finepm ) then
               ugpmF  = ugpmF   + ugtmp
               if(is_ssalt) ugSSaltF = ugSSaltF  +  ugtmp
               if(is_dust ) ugDustF  = ugDustF   +  ugtmp
               if( find_index( ispec, SIA_GROUP )>0 ) &
                    ugsiaPM  = ugsiaPM + ugtmp
             else
                ugpmC  = ugpmC   + ugtmp
               if(is_ssalt) ugSSaltC = ugSSaltC  +  ugtmp
               if(is_dust ) ugDustC  = ugDustC   +  ugtmp
             end if

             if (  species(ispec)%name == 'SO4' ) then
               ugSO4 = ugSO4  +  ugtmp
             else if ( index( species(ispec)%name, 'NO3_F' )>0) then
               ugNO3f= ugNO3f +  ugtmp
             else if ( index( species(ispec)%name, 'NO3_C' )>0) then
               ugNO3c= ugNO3c +  ugtmp
             end if

             if( debug_flag .and. k==20 ) write(*, *) &
                  sub//"UGSIA:"//trim(species(ispec)%name), &
                     find_index( ispec, SIA_GROUP ), ugtmp, ugpmF, ugpmC
           end do

         ! FRACTIONS used for N2O5 hydrolysis
         ! We use mass fractions, since we anyway don't have MW for OM, dust,...
         !  ugRemF will include OM, EC, PPM, Treat as OM 

          ugRemF = ugpmf - ugSIApm -ugSSaltF -ugDustF 

          aero_fss(k)     = ugSSaltF/ugpmF
          aero_fdust(k)   = ugDustF/ugpmF
          aero_fom(k)     = max(0.0, ugRemF)/ugpmF

          if( DEBUG%SETUP_1DCHEM ) then ! extra checks 
             if( aero_fom(k) > 1.0 .or. ugRemF < -1.0e-9 ) then
                print "(a,i4,99es12.3)", sub//"AERO-F ", k, &
                  aero_fom(k), ugRemF,ugpmF, ugSIApm, ugSSaltF, ugDustF
                call CheckStop(sub//"AERO-F problem " )
              end if
          end if

         ! GERBER equations for wet radius

          do iw = 1, AERO%NSAREA 

           Ddry(iw) =  AERO%DpgN( AERO%Inddry(iw))   ! (m)

           if( AERO%Gb(iw) > 0 ) then
              DpgNw(iw,k)  = 2*WetRad( 0.5*Ddry(iw), rh(k), AERO%Gb(iw) ) 

           if( debug_flag .and. k==20 ) write(*,fmt) sub//"WRAD  ", iw, &
               1.0*AERO%Gb(iw), rh(k), Ddry(iw), DpgNw(iw,k), &
                DpgNw(iw,k)/Ddry(iw), WetRad( 0.5*Ddry(iw), rh(k), AERO%Gb(iw) )

           else ! index -1 indicates use dry (for dust)
              DpgNw(iw,k)  = Ddry(iw)
           end if
          end do

           iw= AERO%SIA_F
           rho=AERO%PMdens(AERO%Inddry(iw))
           S_m2m3(iw,k) = pmSurfArea(ugsiaPM,Dp=Ddry(iw), Dpw=DpgNw(iw,k),  &
                                     rho_kgm3=rho )

           iw= AERO%PM_F ! now use for fine PM
           rho=AERO%PMdens(AERO%Inddry(iw))
           S_m2m3(iw,k) = pmSurfArea(ugpmf,Dp=Ddry(iw), Dpw=DpgNw(iw,k),  &
                                     rho_kgm3=rho )

           iw= AERO%SS_F
           rho=AERO%PMdens(AERO%Inddry(iw))
           S_m2m3(iw,k) = pmSurfArea(ugSSaltF,Dp=Ddry(iw), Dpw=DpgNw(iw,k),  &
                                     rho_kgm3=rho )

           iw= AERO%SS_C
           rho=AERO%PMdens(AERO%Inddry(iw))
           S_m2m3(iw,k) = pmSurfArea(ugSSaltC,Dp=Ddry(iw), Dpw=DpgNw(iw,k),  &
                                     rho_kgm3=rho )

          ! dust - just used dry radius
           iw= AERO%DU_F
           rho=AERO%PMdens(AERO%Inddry(iw))
           S_m2m3(iw,k) = pmSurfArea(ugDustF,Dp=Ddry(iw), Dpw=DpgNw(iw,k),  &
                                     rho_kgm3=rho )

           iw= AERO%DU_C
           rho=AERO%PMdens(AERO%Inddry(iw))
           S_m2m3(iw,k) = pmSurfArea(ugDustC,Dp=Ddry(iw), Dpw=DpgNw(iw,k),  &
                                     rho_kgm3=rho )

!           call CheckStop (  S_m2m3(k) < 0.0 , "NEGS_m2m3" )

           if( debug_flag .and. k==20 )  then
            write(*,fmt)  sub//" SAREAugPM in  ", k,  rh(k), temp(k), &
              ugsiaPM, ugpmf, ugSSaltC, ugDustC
            do iw = 1, AERO%NSAREA
             write(*,fmt) sub//"GERB ugDU (S um2/cm3)  ", iw, Ddry(iw), &
              DpgNw(iw,k)/Ddry(iw), 1.0e6*S_m2m3(iw,k)
            end do
           end if

           S_m2m3(:,k) = min( S_m2m3(:,k), 6.0e-3)  !! Allow max 6000 um2/cm3

           ! For coarse, we simply sum. We ignore some non-SS or dust _C.
           iw= AERO%PM_C
           S_m2m3(iw,k) = S_m2m3(AERO%PM_F,k) + S_m2m3(AERO%SS_C,k) + S_m2m3(AERO%DU_C,k) 

           iw= AERO%ORIG
           S_m2m3(iw,k) = S_RiemerN2O5(k)

           ! m2/m3 -> um2/cm3 = 1.0e6, only for output to netcdf
           if( k == KMAX_MID ) then 
              do iw = 1, AERO%NSAREA
                SurfArea_um2cm3(iw,i,j) = 1.0e6* S_m2m3(iw,k)
              end do
           end if

          end if ! GERBER

      ! End Surf Area
    
   end do ! k

! Check that concentrations are not "contaminated" with NaN
   if ( isnan( xn_2d(NO2,KMAX_MID) ) ) then
      print *, "NANAN ", trim(species(NO2)%name)
      call CheckStop( "Detected non numerical concentrations (NaN)")
   end if

   o2(:) = 0.21 *amk(:)
   n2(:) = amk(:) - o2(:)
!   o2(:) = 0.2095 *amk(:) ! more exact, but prefer o3+n2 to add to 100%
!   n2(:) = 0.7808 *amk(:)
   tinv(:) = 1./temp(:)


   cn2o5(:) = cMolSpeed(temp(:),108.0)
   chno3(:) = cMolSpeed(temp(:), 63.0)
   cho2(:)  = cMolSpeed(temp(:), 33.0)
   co3(:)   = cMolSpeed(temp(:), 48.0)

  ! 5 ) Rates  (!!!!!!!!!! NEEDS TO BE AFTER RH, XN, etc. !!!!!!!!!!)


   call set_rct_rates()



   if ( first_call ) then
     call CheckStop( any(isnan(rct(:,:))), sub//"RCT NAN'd")
     call CheckStop( any(rct(:,:) < 0.0 ), sub//"RCT NEG'd") !dsJUL2015
     nd2d = 0
     do itmp = 1, size(f_2d)
           if ( f_2d(itmp)%subclass == 'rct' ) then
             nd2d =  nd2d  + 1
             call CheckStop(nd2d>size(id2rct),sub//"Need bigger id2rct array")
             d2index(nd2d)= itmp
             id2rct(nd2d) = f_2d(itmp)%index  !index of rate constant (config)
             if(MasterProc) write(*,*) 'RCTFOUND', itmp, nd2d, id2rct(nd2d)
           end if
     end do

     first_call = .false.
   end if ! first_call

   do itmp = 1, nd2d
       d_2d(d2index(itmp),i,j,IOU_INST) =  rct(id2rct(itmp),KMAX_MID)
       if( debug_flag ) then
          write(*,"(a,6i5,es12.3)") sub//"OUTRCT "//&
           trim(f_2d(d2index(itmp))%name), me,i,j,itmp,&
           d2index(itmp),id2rct(itmp), d_2d(d2index(itmp),i,j,IOU_INST)
       end if
   end do


   end subroutine setup_1d
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine setup_rcemis(i,j)
!-------------------------------------------------------------------
!    DESCRIPTION:
!    Extracts emissions in column from gridrcemis, for input to chemistry
!    routines. Results in "rcemis" array
!-------------------------------------------------------------------
  !-- arguments
  integer, intent(in) ::  i,j     ! coordinates of column

  !  local
  integer ::  iqrc,k, itot
  real    :: fac, eland   ! for Pb210  - emissions from land

  integer ::  i_help,j_help,i_l,j_l, i_Emis_4D,n
  logical, save     :: first_call = .true. 
  character(len=13) :: sub="setup_rcemis:"

  if(first_call)then
    inat_RDF = find_index( "DUST_ROAD_F", EMIS_BioNat(:) )
    inat_RDC = find_index( "DUST_ROAD_C", EMIS_BioNat(:) )
    itot_RDF = find_index( "DUST_ROAD_F", species(:)%name    )
    itot_RDC = find_index( "DUST_ROAD_C", species(:)%name    )
    itot_Rn222=find_index( "RN222", species(:)%name    )
    first_call = .false.
  endif 

! initilize ! initilize ! initilize ! initilize
  rcemis(:,:)=0.
! initilize ! initilize ! initilize ! initilize

  forall(k=KEMISTOP:KMAX_MID,iqrc=1:NRCEMIS) &
    rcemis(iqrc2itot(iqrc),k) = gridrcemis(iqrc,k,i,j)

  ! Volcanic emissions (SO2 and ASH),
  ! and Contribution from Emergeny scenarios
  if(VOLCANO_SR)then
    rcemis(:,:)=rcemis(:,:)+ColumnRate(i,j,REDUCE_VOLCANO=0.85)
  else
    rcemis(:,:)=rcemis(:,:)+ColumnRate(i,j)
  endif

  ! lightning and aircraft ... Airial NOx emissions if required:
  if(USE_LIGHTNING_EMIS)then
    do k=KCHEMTOP, KMAX_MID
      rcemis(NO ,k) = rcemis(NO ,k) + 0.95 * airlig(k,i,j)
      rcemis(NO2,k) = rcemis(NO2,k) + 0.05 * airlig(k,i,j)
    enddo
  endif
  if(USE_AIRCRAFT_EMIS) then
    do k=KCHEMTOP, KMAX_MID
      rcemis(NO ,k) = rcemis(NO ,k) + 0.95 * airn(k,i,j)
      rcemis(NO2,k) = rcemis(NO2,k) + 0.05 * airn(k,i,j)
    enddo
  end if ! AIRCRAFT NOX
  if(DEBUG%SETUP_1DCHEM.and.debug_proc.and.i==debug_li.and.j==debug_lj)&
    write(*,"(a,2L2,10es10.3)") &
      sub//"AIRNOX ", USE_LIGHTNING_EMIS, USE_AIRCRAFT_EMIS, &
      airn(KMAX_MID,i,j),airlig(KMAX_MID,i,j)

  ! Add sea salt production
  if(USE_ROADDUST.and.itot_RDF>0) then  ! Hard-code indices for now
    rcemis(itot_RDF,KMAX_MID) = gridrcroadd(1,i,j)
    rcemis(itot_RDC,KMAX_MID) = gridrcroadd(2,i,j)
   endif

  if(USES%FOREST_FIRES) then
    if(burning(i,j))call Fire_rcemis(i,j)
  endif

  !Soil NOx
  if(USE_GLOBAL_SOILNOX)then !NEEDS CHECKING NOV2011
    rcemis(NO,KMAX_MID)=rcemis(NO,KMAX_MID)+SoilNOx(i,j)
  endif

  if(Found_Emis_4D>0)then
     do i_Emis_4D=1,N_Emis_4D
        if(emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D)=='NOTSET')exit
        n=find_index(emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D),species(:)%name)
        if(n>0)then
           fac=1.0/1000000.0/3600.0 !convert from Bq/m3/hour into Bq/cm3/s
           do k=KCHEMTOP, KMAX_MID
              rcemis(n,k)=rcemis(n,k)+Emis_4D(i,j,k,i_Emis_4D)*fac
           enddo
        endif
     enddo
  endif


  do k=KCHEMTOP, KMAX_MID

     deltaZcm(k) = 100*( z_bnd(i,j,k)-z_bnd(i,j,k+1) )

    if(DEBUG%SETUP_1DCHEM.and.debug_proc.and.i==debug_li.and.j==debug_lj) then
       write(*,"(a,i3,9es12.4)") "1DCHEM DZ",  k, deltaZcm(k), &
        100*(dA(k)+dB(k)*ps(i,j,1))/(GRAV*roa(i,j,k,1))! , &
!        100*(dA(k)+dB(k)*ps(i,j,1))/(GRAV*amk(k)*ATWAIR)!, &
!        amk(k)*ATWAIR/roa(i,j,k,1)
!       amk(k) = roa(i,j,k,1) * to_number_cm3  ! molecules air/cm3
    end if
  end do


  if(EmisSplit_OUT)then
    !put all added emissions in EmisSplit_OUT, also natural emissions
    SumSplitEmis(i,j,:)=0.0
    do k=KCHEMTOP, KMAX_MID
      do iqrc=1,nrcemis
        !give unit mg/m2/s dt_advec multiplied in Derived_ml 
        itot=iqrc2itot(iqrc)
        SumSplitEmis(i,j,iqrc) = SumSplitEmis(i,j,iqrc)&
          +rcemis(itot,k)*species(itot)%molwt &
          *(dA(k)+dB(k)*ps(i,j,1))/(GRAV*amk(k)*ATWAIR)
      enddo
    enddo
  endif


  ! Soil Rn222 emissions from non-ice covered land, + water
  ! at rate of 1 atom/cm2/s
   eland = 1.0 - water_fraction(i,j) - ice_landcover(i,j)


  ! z_bnd is in m, not cm, so need to divide by 100.
  if(itot_Rn222>0) &
    rcemis(itot_Rn222,KMAX_MID) = (0.00182*water_fraction(i,j)+eland)&
      /deltaZcm(KMAX_MID) 
   !!   /((z_bnd(i,j,KMAX_BND-1)-z_bnd(i,j,KMAX_BND))*100.)

!ESX     rc_Rnwater(KMAX_MID) = water_fraction(i,j)  / &
!ESX            ((z_bnd(i,j,KMAX_BND-1) - z_bnd(i,j,KMAX_BND))*100.)

endsubroutine setup_rcemis
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine reset_3d(i,j)
  integer, intent(in) :: i,j
  integer :: k, n, ispec, id    ! loop variables
! ! XNCOL testing -- sets d_2d for column data from molec/cm3 concs.
 ! if variables are wanted for d_2d output (via USET), we use these indices:
  character(len=*),parameter :: sub='reset3dxncol:'
  character(len=10) :: specname
  integer, dimension(10), save :: d2index, id2col
  logical, save :: first_call = .true.
  integer, save :: nd2d
!XNCOL  end testing


  do k = KCHEMTOP, KMAX_MID
    ! 1)/ Short-lived species - no need to scale with M
    do n = 1, NSPEC_SHL
      xn_shl(n,i,j,k) = xn_2d(n,k)
    enddo ! ispec

    ! 2)/ Advected species
    do n = 1, NSPEC_ADV
      ispec = NSPEC_SHL + n
      xn_adv(n,i,j,k) = xn_2d(ispec,k)/amk(k)
    enddo ! ispec
  enddo ! k

!XNCOL !======================================================================
!! If column totals are wanted, we can do those here also since xn_2d are
!! in molec/cm3, and we want molec/cm2:

   if ( first_call ) then

     nd2d = 0
     do id = 1, size(f_2d)
           !if(MasterProc) write(*,*) 'USET XNCOL SRCH', id, f_2d(id)%subclass
           if ( f_2d(id)%subclass == 'xncol' ) then
             nd2d =  nd2d  + 1
             call CheckStop( nd2d > size(id2col), &
                 sub//"Need bigger id2col array" )
             specname = f_2d(id)%name(7:)  ! Strip XNCOL_
             ispec = find_index( specname, species(:)%name )
             call CheckStop(ispec < 1, sub//"XNCOL not found"//specname )
             d2index(nd2d)= id
             id2col(nd2d) = ispec
             if(MasterProc) write(*,*) 'USET XNCOL FOUND', id, ispec, &
                 trim(specname),nd2d, id2col(nd2d)
           end if
     end do
     first_call = .false.
   end if

   do id = 1, nd2d
        ispec = id2col(id)
        d_2d(d2index(id),i,j,IOU_INST) = dot_product(xn_2d(ispec,:),deltaZcm(:))
    if(DEBUG%SETUP_1DCHEM.and.debug_proc.and.i==debug_li.and.j==debug_lj) then
!    !if( debug_flag ) then
       write(*,"(a,6i5,a,9es12.3)") sub//"OUTXNCOL "//trim(f_2d(d2index(id))%name), me,&
            i,j,id,d2index(id),ispec, trim(species(ispec)%name) &
         ,xn_2d(ispec,20), deltaZcm(20), d_2d(d2index(id),i,j,IOU_INST)
end if
    end do
!!XNCOL


endsubroutine reset_3d
!---------------------------------------------------------------------------
endmodule Setup_1d_ml
!_____________________________________________________________________________!
