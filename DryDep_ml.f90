! <DryDep_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
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
module DryDep_ml

! Dry deposition scheme uses a mosaic approach.

! Latest documentation and ACP eqn references in code below from
!  Simpson, D., Benedictow, A., Berge, H., Bergstr\"om, R., Emberson, L. D.,
!  Fagerli, H., Flechard, C. R.,  Hayman, G. D., Gauss, M., Jonson, J. E., 
!  Jenkin, M. E., Ny\'{\i}ri, A., Richter, C., Semeena, V. S., Tsyro, S.,
!  Tuovinen, J.-P., Valdebenito, \'{A}., and Wind, P.: 
!  The EMEP MSC-W chemical transport model -- technical description, 
!  Atmos. Chem. Phys., 12, 7825--7865, 2012.
!
! History
! Module started from the drag-coefficient based approach of BJ98:
! Berge, E. and  Jakobsen, H.A., A regional scale multi-layer model
! for the calculation of long-term transport and deposition of air
! pollution in Europe, Tellus B (1998), 50, 105-223.

! ** but ** extensively re-written since. See ....
! Emberson, L.,Simpson, D.,Tuovinen, J.-P.,Ashmore, M.R., Cambridge, H.M.",
!  2000, Towards a model of ozone deposition and stomatal uptake over 
!  Europe, EMEP MSC-W Note 6/2000,
! Simpson, D.,Tuovinen, J.-P.,Emberson, L.D.,Ashmore, M.R.,2001,
!  "Characteristics of an ozone deposition module",WASP:Focus,1,253-262
! Simpson, D.,Tuovinen, J.-P.,Emberson, L.D.,Ashmore, M.R.,2003,
!  "Characteristics of an ozone deposition module II: sensitivity analysis",
!   WASP, 143, 123-137
! Tuovinen, J.-P.,Ashmore, M.R.,Emberson, L.D.,Simpson, D., 2004, "Testing
!   and improving the EMEP ozone deposition module", Atmos.Env.,38,2373-2385
!
! Also, handling of dry/wet and co-dep procedure changed following discussions
! with CEH.

! Autumn/Winter 2017-2018
! FUTURE: BiDir functionality (Dave/Roy Wichink Kruit) using Roy's methods
  
use Aero_Vds_ml,      only: SettlingVelocity, GPF_Vds300, Wesely300
use BiDir_emep
use BiDir_module
use Biogenics_ml,     only: SoilNH3  ! for BiDir
use CheckStop_ml,     only: CheckStop, StopAll
use Chemfields_ml ,   only: cfac, so2nh3_24hr,Grid_snow 
use ChemSpecs                ! several species needed
use Config_module,    only: dt_advec,PT, K2=> KMAX_MID, NPROC, &
                            DEBUG, DEBUG_ECOSYSTEMS, DEBUG_VDS,&
                            USES, AERO, &
                            USE_SOILNOX, &
                            MasterProc, &
                            PPBINV,&
                            KUPPER, NLANDUSEMAX
use DO3SE_ml,         only: do3se
use EcoSystem_ml,     only: EcoSystemFrac, Is_EcoSystem,  &
                            NDEF_ECOSYSTEMS, DEF_ECOSYSTEMS
use GasParticleCoeffs_ml         ! ... Init_GasCoeff, DRx, Rb_Cor, ...
use GridValues_ml ,   only: GRIDWIDTH_M,xmd,xm2, glat,dA,dB, &
     glon,   debug_proc, debug_li, debug_lj, i_fdom, j_fdom   ! for testing

!use Io_Nums_ml,       only: IO_SPOD
use Io_Progs_ml,      only: datewrite
use Landuse_ml,       only: SetLandUse, Land_codes  & 
                           ,NLUMAX &  ! Max. no countries per grid
                           ,LandCover   ! Provides codes, SGS, LAI, etc,
use LandDefs_ml,      only: LandType, LandDefs, STUBBLE
use LocalVariables_ml,only: Grid, L, iL ! Grid and sub-scale Met/Veg data
use LocalVariables_ml,only: NLOCDRYDEP_MAX ! Used to store Vg
use MassBudget_ml,    only: totddep
use MetFields_ml,     only: u_ref, rh2m, sst
use MetFields_ml,     only: tau, sdepth, SoilWater_deep, th,pzpbl
use MicroMet_ml,      only: AerRes, Wind_at_h
use MosaicOutputs_ml,     only: Add_MosaicOutput, MMC_RH
use OwnDataTypes_ml,      only: depmap
use Par_ml,               only: limax,ljmax, me,li0,li1,lj0,lj1
use PhysicalConstants_ml, only: ATWAIR,PI,KARMAN,GRAV,RGAS_KG,CP,AVOG,NMOLE_M3
use Rb_ml,                only: Rb_gas
use Rsurface_ml,          only: Rsurface, Rinc
use Setup_1dfields_ml,    only: xn_2d,amk, Fpart, Fgas
use Sites_ml, only : nlocal_sites, site_x, site_y, site_name, site_gn
use SoilWater_ml,         only: fSW !  =1.0 unless set by Met_ml
use StoFlux_ml,  only:   unit_flux, &! = sto. flux per m2
                        lai_flux,  &! = lai * unit_flux
                        Setup_StoFlux, Calc_StoFlux  ! subs
use SubMet_ml,            only: Sub
use TimeDate_ml,          only: daynumber, current_date

implicit none
private

public  :: DryDep, init_drydep
private :: Init_DepMap

!integer, private, save :: P != IO_SPOD + me
! Maps from adv index to one of calc indices
integer, public, save, dimension(NSPEC_ADV) :: DepAdv2Calc 

logical, private, save :: my_first_call = .true.
character(len=30),private, save :: errmsg = "ok"

! WE NEED A FLUX_CDDEP, FLUX_ADV FOR OZONE;
! (set to one for non-ozone models)

integer, public, parameter :: FLUX_CDDEP  = CDDEP_O3
integer, public, parameter :: FLUX_ADV   = IXADV_O3
integer, public, parameter :: FLUX_TOT   = O3

! WE ALSO NEED NO3_f and NH4_f for deposition.
! (set to one for non-no3/nh4 models)

integer, public, parameter :: pNO3  = NO3_f
integer, public, parameter :: pNH4  = NH4_f

!logical, public, parameter :: COMPENSATION_PT = .false. 

!***************************************************************************
!  Specifies which of the possible species (from DryDepDefs list)
!  are required in the current air pollution model   
!***************************************************************************
! .... Define the mapping between the advected species and
!      the specied for which the calculation needs to be done.
!  We also define the number of species which will be deposited in
! total, NDRYDEP_ADV. This number should be >= NDRYDEP_GASES
! The actual species used and their relation to the CDDEP_ indices
! above will be defined in Init_DepMap

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 include 'CM_DryDep.inc'

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

logical, public, dimension(NDRYDEP_ADV), save :: vg_set 

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine init_drydep

  integer, save ::  old_daynumber = -99
  integer ::nadv,n
  integer :: i, j, lc, ilc, nlc, iEco   ! for EcoSystem stuff
  logical :: debug_flag  ! for EcoSystem stuff
  real    :: coverage    ! for EcoSystem stuff

  if ( my_first_call ) then 

     call Init_DepMap()               ! Maps CDDEP to IXADV
     call Init_GasCoeff()             ! Sets DryDepDefs coeffs.

     call CheckStop( NLOCDRYDEP_MAX < NDRYDEP_CALC, &
        "Need to increase size of NLOCDRYDEP_MAX" )

     nadv = 0
     do n = 1, NDRYDEP_ADV  
         nadv       = max( DDepMap(n)%ind, nadv )  ! Looking for highest IXADV
         vg_set(n)  = ( DDepMap(n)%calc == CDDEP_SET ) ! for set vg
     end do

     my_first_call = .false.
     if( MasterProc  .and. DEBUG%DRYDEP) write(*,*) "INIT_DRYDEP day ", &
           daynumber, old_daynumber

!=============================================================================
! From EcoSystems, but caused some circularity problem
! use EcoSystem_ml, only :: EcoSystemFrac, Is_EcoSystem

   EcoSystemFrac(:,:,:) = 0.0
   do j = 1, ljmax
     do i = 1, limax
       debug_flag = ( DEBUG_ECOSYSTEMS .and. debug_proc .and. &
          i == debug_li .and. j == debug_lj )

       nlc = LandCover(i,j)%ncodes
       LCLOOP: do ilc= 1, nlc
           lc       = LandCover(i,j)%codes(ilc)
           coverage = LandCover(i,j)%fraction(ilc)
           ECOLOOP: do iEco= 1, NDEF_ECOSYSTEMS
              if( Is_EcoSystem(iEco,lc) ) then
                 EcoSystemFrac(iEco,i,j) = EcoSystemFrac(iEco,i,j) + coverage
                if( debug_flag ) then
                     write(6,"(a,2i4,a12,3f10.4)") "ECOSYS AREA ",&
                     ilc, lc, "=> "//trim(DEF_ECOSYSTEMS(iEco)), &
                         coverage, EcoSystemFrac(iEco,i,j)
                end if
             end if
          end do ECOLOOP
        end do LCLOOP
      end do ! i
    end do ! j

!     invEcoFrac(:) = 0.0
!
!     do n = 0, size(DEF_ECOSYSTEMS)-1
!        if ( EcoFrac(n) > 1.0e-39 ) invEcoFrac(n) = 1.0/EcoFrac(n)
!     end do
!=============================================================================
  end if !  my_first_call

  end subroutine init_drydep
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Init_DepMap
   integer :: iadv, i

     do i = 1, NDRYDEP_ADV  ! 22
      iadv = DDepMap(i)%ind
      if(DEBUG%DRYDEP .and. MasterProc) &
         write(6,*) "DEPMAP   ", DDepMap(i)%ind, DDepMap(i)%calc
      call CheckStop( iadv < 1, "ERROR: Negative iadv" )
      DepAdv2Calc(iadv) = DDepMap(i)%calc
    end do
  
   ! We process the various combinations of gas-species and ecosystem:
   ! starting with DryDep, e.g. DDEP_SO2_m2CF
 
     if(MasterProc.and.DEBUG%DRYDEP) write(6,*) "Init_DepMap D2D FINISHED"

  end subroutine Init_DepMap


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine DryDep(i,j)
    integer, intent(in):: i,j
    real, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost

    real, dimension(NDRYDEP_GASES ) :: &
          Rb           & ! Quasi-boundary layer rsis.
         ,Rsur         & ! Surface Resistance (s/m) 
         ,Gsto           ! Stomatal conductance (big-leadf)
         
    real, dimension(NDRYDEP_CALC) :: &
          gradient_fac & ! Ratio of conc. at zref (ca. 50m) and 3m
         ,vg_fac       & ! Loss factor due to dry dep.
         ,Vg_ref       & ! Vg at ref ht.
         ,Vg_3m        & ! Vg at  3m
         ,Vg_ratio     & ! Ratio Vg_ref/Vg_3m = ratio C(3m)/C(ref), over land
         ,sea_ratio     ! Ratio Vg_ref/Vg_3m = ratio C(3m)/C(ref), over sea

    character(len=*), parameter :: dtxt='DryDep:' ! debug label
    logical, save      :: dbg, dbghh, dbgBD
    integer n, iiL, nlu, ncalc, nadv, nFlux  ! help indexes
    integer :: imm, idd, ihh, iss     ! date
    integer :: ntot !index of adv species in xn_2d array

    real :: no2fac  ! Reduces Vg for NO2 in ration (NO2-4ppb)/NO2

    real convfac,  & ! rescaling to different units
         convfac2, & ! rescaling to different units
         lossfrac,  & !  If needed in My_DryDep - not used now.
         dtz         ! scaling factor for veff ( = dt/z, where dt=timestep and 
                     ! z = height of layer)

    integer :: nae

    real, save :: inv_gridarea  ! inverse of grid area, m2

    real ::  Sumcover, Sumland   ! Land-coverage
    logical :: debug_flag        ! set true when i,j match DEBUG_i, DEBUG_j
    real :: Vg_scale

    real, dimension(NSPEC_ADV ,NLANDUSEMAX):: fluxfrac_adv
    integer, dimension(NLUMAX)  :: iL_used, iL_fluxes
!BIDIR SKIP    real :: wet, dry         ! Fractions
    real :: snow_iL          !snow_flag fraction for one landuse
    real :: Vds              ! Aerosol near-surface deposition rate (m/s)
    real :: no3nh4ratio      ! Crude NH4/NO3 for Vds ammonium 

    real :: c_hveg, Ra_diff, surf_ppb  ! for O3 fluxes and Fst where needed
    real :: c_hveg3m, o3_45m  ! TESTS ONLY
! temporary for POD/SPOD
!    logical, parameter :: SPOD_OUT = .false.  ! MAKES HUGE FILES. Not for routine use!
!    logical, save      :: first_spod = .true.
    character(len=20), save :: fname
    integer :: nglob
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Extra outputs sometime used. Important that this 
!! line is kept at the end of the variable definitions and the start of real
!!  code - allows both in .inc file
!! Uncomment and make .inc file as required
!   include 'EXTRA_LU_Setup.inc'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     first calculate the 3m deposition velocity following the same
!     procedure as in the lagmod. second the flux and the accumulated 
!     deposition is calculated.
!
!     effective dry deposition velocity applied to the model concentration
!     at the top of the constant flux layer, zdep 
!     Dry deposion rates are specified in subroutine readpar
!

! FOR DEBUGGING
  imm      =    current_date%month            ! for debugging
  idd      =    current_date%day              ! for debugging
  ihh      =    current_date%hour             ! for debugging
  iss      =    current_date%seconds          ! for debugging

     inv_gridarea = 1.0/(GRIDWIDTH_M*GRIDWIDTH_M) 


   ! - Set up debugging coordinates first. ---------------------------!
   ! If location matches debug i,j value, set debug_flag. Also passed
   ! to Rsurface_ml

    debug_flag= ( debug_proc .and. i == debug_li .and. j == debug_lj) 
    dbg       =  DEBUG%DRYDEP .and. debug_flag 
    dbghh     =  dbg .and. iss == 0 
    dbgBD     =  DEBUG%BIDIR .and. debug_flag .and. iss == 0
    if (dbgBD) print *, "BIDIR TEST ", me, debug_flag


   ! -----------------------------------------------------------------!
   !.and conversion factor,  convfac (( ps-pt)/grav... )  ===> 
   !      pressure in kg m-1 s-2
   ! 
    convfac = (dA(K2) + dB(K2)*Grid%psurf)&!dP
               *xmd(i,j)/(ATWAIR*GRAV*inv_gridarea)

   ! -----------------------------------------------------------------!
   ! conver molecules/cm3 to ppb for surface:
    surf_ppb   = PPBINV /amk(K2)
    if ( DEBUG%AOT .and. debug_flag ) write(*,"(a,es12.4)") "CHAMK", surf_ppb

   ! -----------------------------------------------------------------!

!     !.and factor,  kg_air_ij (( ps-pt)/grav... )  ===> 
!     !      pressure in kg m-1 s
!     !      used for converting from mixing ratio to kg
!
!      kg_air_ij = (ps(i,j,1) - PT)*carea(K2) = dP*dx**2/g
!     ! -----------------------------------------------------------------!
!
    lossfrac = 1.0 !  Ratio of xn before and after deposition


    dtz      = dt_advec/Grid%DeltaZ

    if ( dbghh ) call datewrite(dtxt//"DMET ", daynumber, (/ Grid%zen,  &
         1.0e-5*Grid%psurf, Grid%Idiffuse, Grid%Idirect, &
         Grid%Hd, Grid%LE, Grid%invL, Grid%ustar /) )
      
   !/ Initialise Grid-avg Vg for this grid square:
     !call Init_GridMosaic(i,j)

    Vg_ratio(:) = 0.0
    Sumcover = 0.0
    Sumland  = 0.0
    fluxfrac_adv (:,:) = 0.0
    Grid_snow(i,j)=0.0 

    Sub(0)%Gsur = 0.0
    Sub(0)%Gsto = 0.0
    !Sub(0)%Gns  = 0.0
    Sub(0)%Vg_Ref  = 0.0
    Sub(0)%Vg_3m   = 0.0
 
    !/ SO2/NH3 for Rsur calc
    Grid%so2nh3ratio = xn_2d(SO2,K2) / max(1.0,xn_2d(NH3,K2))

    Grid%so2nh3ratio24hr = so2nh3_24hr(i,j)

   !---------------------------------------------------------
   !> NH4NO3 deposition will need this ratio for the NH4 part

    no3nh4ratio = 1.0
    if( xn_2d(pNH4,K2) > 1.0  ) then
       no3nh4ratio = xn_2d(pNO3,K2) / xn_2d(pNH4,K2)
       no3nh4ratio = min( 1.0,  no3nh4ratio )
    end if

   !---------------------------------------------------------
   ! If we do not have soil NO emissions active, we use instead a 
   ! surrogate for NO2 compensation point approach, assuming 
   ! c.p.=4 ppb (actually use 1.0e11 #/cm3):        
   ! Aiming at dC/dt = -Vg.(C-4)/dz instead of -Vg.C/dz
   ! factor difference is then:   C-4/C in ppb units
   ! Note, xn_2d has no2 in #/cm-3
 
    no2fac = 1.0
    if ( .not. USE_SOILNOX ) then
      no2fac = max( 1.0, xn_2d(NO2,K2) )
      no2fac = max(0.00001,  (no2fac-1.0e11)/no2fac)
    end if
   !---------------------------------------------------------

    if ( dbghh ) then
      write(*,"(a,2i4,4es12.4)") "DRYDEP CONCS SO2,NH3,O3 (ppb) ", i,j, &
        xn_2d(SO2,K2)*surf_ppb, xn_2d(NH3,K2)*surf_ppb, &
        xn_2d(O3,K2)*surf_ppb, no2fac
    end if

    call Setup_StoFlux( daynumber )

! - can set settling velcoty here since not landuse dependent

    do nae = 1, AERO%NSIZE
      AERO%Vs(nae) = SettlingVelocity( Grid%t2, Grid%rho_ref, &
                       AERO%sigma(nae), AERO%DpgV(nae), AERO%PMdens(nae) )
    end do

    if ( dbghh ) call datewrite(dtxt//"DRYDEP VS",AERO%NSIZE,&
                                  (/ Grid%t2, Grid%rho_ref, AERO%Vs /) )

    !/ And start the sub-grid stuff over different landuse (iL)

    nlu = LandCover(i,j)%ncodes
    nFlux = 0           ! Numberof LC needing flux calcs
    LULOOP: do iiL= 1, nlu
        iL      = LandCover(i,j)%codes(iiL)

        iL_used (iiL) = iL    ! for eco dep

        if ( LandType(iL)%flux_wanted ) then
          nFlux = nFlux + 1
          iL_fluxes (nFlux) = iL    ! for eco dep
        end if

        Sub(iL)%f_phen = LandCover(i,j)%fphen(iiL) ! for SPOD_OUT
        Sub(iL)%f_sun  = 0.0 ! for SPOD_OUT
        Sub(iL)%g_sun  = 0.0 ! for SPOD_OUT
        Sub(iL)%g_sto  = 0.0 ! for SPOD_OUT
        Sub(iL)%f_temp = 0.0 ! for SPOD_OUT. Can 
        Sub(iL)%f_vpd  = 0.0 ! for SPOD_OUT

        Sub(iL)%SGS = LandCover(i,j)%SGS(iiL)   !used for AOT CHECK?
        Sub(iL)%EGS = LandCover(i,j)%EGS(iiL)

        L = Sub(iL)    ! ! Assign e.g. Sub(iL)ustar to ustar

             if ( dbghh ) then
                write(6,"(a,3i3,f6.1,2i4,3f7.3,i4,9f8.3)") "DVEG: ", &
                    nlu,iiL, iL, glat(i,j), L%SGS, L%EGS, &
                   L%coverage, L%LAI, L%hveg,daynumber, &
                   Grid%sdepth, fSW(i,j),L%fSW,L%t2C

                write(6,"(a,i4,3f7.2,7es10.2)") "DMET SUB", &
                  iL, Grid%ustar, L%ustar, L%rh,  Grid%invL, &
                  L%invL, L%Ra_ref, L%Ra_3m
             end if


         call Rb_gas(L%is_water, L%ustar, L%z0, DRYDEP_GASES ,Rb)

         call Rsurface(i,j,DRYDEP_GASES ,Gsto,Rsur,errmsg,debug_flag,snow_iL)

         if(dbghh) call datewrite(dtxt//"STOFRAC "//LandDefs(iL)%name, &
                 iL, (/ Gsto(2), Rsur(2), Gsto(2)*Rsur(2)  /) ) ! 2 is for WES_O3

           !Sub(iL)%g_sto = L%g_sto   ! needed elsewhere
           !Sub(iL)%g_sun = L%g_sun
         Sub(iL) = L  !Resets Sub with new L values frmom Rsurface

         Grid_snow(i,j) = Grid_snow(i,j) +  L%coverage * snow_iL  ! QUERY why?


       !/... add to grid-average Vg:

!BIDIR SKIP         wet =   Grid%wetarea  ! QUERY Now used for BIDIR
!BIDIR SKIP         dry =   1.0 - wet     !  "  "

         do n = 1, NDRYDEP_CALC

           ! ================================================
           if ( n > NDRYDEP_GASES )  then    ! particles

              nae = AERO_SIZE(n) ! See GasParticleCoeffs_ml

              if ( LandType(iL)%is_forest  ) then 

                 !/ Use eqn *loosely* derived from Petroff results
                 ! ACP67-69

                  Vds = GPF_Vds300(L%ustar,L%invL, L%SAI )

              else !!!  ! ACP67-68

                !/  Use Wesely et al  for other veg & sea
                 ! Vds = Nemitz2004( 0.4, L%ustar, L%invL )

                 Vds = Wesely300( L%ustar, L%invL )

              end if

             ! We allow fine NH4NO3 particles to deposit x 3, in
             ! unstable conditions. (F_N in ACP68)
             ! Now, MARS/EQSAM etc only provide NO3_f, NH4_f, but some of the
             ! latter is (NH4)xSO4. Now, as all NO3_f is associated with NH4NO3, 
             ! we can scale any NH4_f deposition with the ratio

              if (n==CDDEP_PMfN .and. L%invL<0.0 ) then 
                   Vds = Vds * 3.0 ! for nitrate-like
              else if (n==CDDEP_PMfNH4 .and. L%invL<0.0 ) then 
                   !Vds = Vds* [ ( 1 - no3ratio) +  3 * no3nh4ratio ]
                   Vds = Vds * (1 + 2 * no3nh4ratio) 
              end if

            ! Use non-electrical-analogy version of Venkatram+Pleim (AE,1999)
            ! ACP70

              Vg_ref(n) =  AERO%Vs(nae)/ ( 1.0 - exp( -( L%Ra_ref + 1.0/Vds)* AERO%Vs(nae)))
              Vg_3m (n) =  AERO%Vs(nae)/ ( 1.0 - exp( -( L%Ra_3m  + 1.0/Vds)* AERO%Vs(nae)))


              if ( DEBUG_VDS ) then
                if ( debug_flag .or. (Vg_3m(n)>0.50 .or. Vg_ref(n)>0.50 )) then
                  write(*,"(a,5i3,2i4,2f7.3,f8.2,20f7.2)") "AEROCHECK", &
                    imm, idd, ihh, n, iL, i_fdom(i), j_fdom(j),&
                    L%ustar,  L%invL, pzpbl(i,j), &
                    Grid%ustar, Grid%Hd,  100.0*Vds, &
                    100.0*AERO%Vs(nae), 100.0*Vg_ref(n),  100.0*Vg_3m (n) &
                     , Grid%t2, Grid%rho_ref 
                   write(*,"(a,2i4,3es10.3)") "VDS CHECK ",n, nae, Vds, Vg_ref(n)
                  
                   call CheckStop((Vg_3m(n)>0.50 .or. Vg_ref(n)>0.50 ), "AEROSTOP")
                end if
              end if

            ! ================================================
           else   ! gases ! NB no wet-dry difference needed here

              Vg_ref(n) = 1. / ( L%Ra_ref + Rb(n) + Rsur(n) ) 
              Vg_3m (n) = 1. / ( L%Ra_3m  + Rb(n) + Rsur(n) ) 

             ! specials, NH3 and NO2:

              if ( USES%BIDIR .and. n == CDDEP_NH3 ) then

                call StopAll(dtxt//'NOT IMPLEMENTED YET')

           
           ! Surrogate for NO2 compensation point approach, 
           ! assuming c.p.=4 ppb (ca. 1.0e11 #/cm3):        
           ! Note, xn_2d has no2 in #/cm-3

              else if ( .not. USE_SOILNOX .and.  n == CDDEP_NO2 ) then
  
                 Vg_ref(CDDEP_NO2) = Vg_ref(CDDEP_NO2) * no2fac
                 Vg_3m (CDDEP_NO2) = Vg_3m (CDDEP_NO2) * no2fac

              end if ! specials, NH3, NO2

              !QUERY - do we need Gsur for anything now?!
              ! StoFrac should end up with weighted mean. Start with Vg*f
   
              Sub(iL)%Gsur(n) = 1.0/Rsur(n) ! Note iL, not iiL 
              Sub(iL)%Gsto(n)  = Gsto(n)      ! Note iL, not iiL 

              Sub(0)%Gsur(n)  =  Sub(0)%Gsur(n) + L%coverage / Rsur(n)
              Sub(0)%Gsto(n)   =  Sub(0)%Gsto(n)  + L%coverage * Gsto(n)
              if( dbghh.and.n==2 ) call datewrite("CmpSto", iL, &
                         (/ Sub(iL)%Gsto(n) / Sub(iL)%Gsur(n) /) )
           end if ! Gases?

           Sub(0)%Vg_ref(n) = Sub(0)%Vg_ref(n) + L%coverage * Vg_ref(n)
           Sub(0)%Vg_3m(n)  = Sub(0)%Vg_3m(n)  + L%coverage * Vg_3m(n)
           Sub(iL)%Vg_ref(n) = Vg_ref(n)
           Sub(iL)%Vg_3m(n) = Vg_3m(n)

         end do ! chemical species loop

         Sumcover = Sumcover + L%coverage

        !/-- only grab gradients over land-areas

         if ( L%is_water ) then
            do n = 1, NDRYDEP_CALC
               sea_ratio(n) =  Vg_ref(n)/Vg_3m(n)
            end do
         else
            Sumland = Sumland + L%coverage
            do n = 1, NDRYDEP_CALC
                Vg_ratio(n) =  Vg_ratio(n) + L%coverage * Vg_ref(n)/Vg_3m(n)
            end do
         end if

         if ( dbghh ) then
            call CheckStop(  Sumland > 1.011, dtxt// "SUMLAND>1")
            do n = CDDEP_O3 , CDDEP_O3 !!! 1,NDRYDEP_GASES 
               call datewrite("DEPO3 ", iL, &
                   (/ Vg_ref(n), Sub(iL)%Vg_ref(n) /) )
                   !(/ Mosaic_VgRef(n,iL) , Vg_ref(n), Sub(iL)%Vg_ref(n) /) )
               call datewrite("DEPDVGA", iL, (/ L%coverage, 1.0*n,& 
                 L%LAI,100*L%g_sto, L%Ra_ref, Rb(n), min( 999.0,Rsur(n) ) /) )
               call datewrite("DEPDVGB", iL, (/ L%coverage, 1.0*n,& 
                100*Vg_3m(n), 100*Vg_ref(n), Vg_ratio(n) /) )
            end do
         end if

       !=======================

         c_hveg   = -999. ! Just for printout when flux_wanted false
         c_hveg3m = -999.

         if (  LandType(iL)%flux_wanted ) then

           n = CDDEP_O3
           Ra_diff = L%Ra_ref - L%Ra_3m
           c_hveg3m = xn_2d(FLUX_TOT,K2)  &     ! #/cm3 units
                        * ( 1.0-Ra_diff*Vg_ref(n) )

         ! Flux = Vg_ref*c_ref = Vg_h * c_h = (c_ref-c_h)/Ra(z_ref,z_h)
         ! which gives:  c_h = c_ref * [ 1-Ra(z_ref,z_h)*Vg_ref ]
         ! Resistance Ra from z_ref to top of canopy:

           Ra_diff  = AerRes(max( L%hveg-L%d, STUBBLE) , Grid%z_ref-L%d,&
                       L%ustar,L%invL,KARMAN)

           c_hveg = xn_2d(FLUX_TOT,K2)  &     ! #/cm3 units
                        * ( 1.0-Ra_diff*Vg_ref(n) )

           if ( DEBUG%AOT .and. debug_flag .and. iL==1 ) then
              call datewrite(dtxt//"CHVEG ", iL, &
                (/  xn_2d(FLUX_TOT,K2)*surf_ppb, c_hveg*surf_ppb,&
                    c_hveg3m * surf_ppb, 100*Vg_ref(n), 100*Vg_3m(n), &
                    L%Ra_ref, (L%Ra_ref-L%Ra_3m), Ra_diff, Rb(n),Rsur(n) /) )
           end if
           

          ! Need to be careful with scope. L is within iL loop, whereas Sub
          ! will be kept throughout i,j calculations:

           Sub(iL)%cano3_ppb   = c_hveg * surf_ppb  ! change units
           Sub(iL)%cano3_nmole = c_hveg * NMOLE_M3  ! units of nmole/m3

         end if !

          !! Extra outputs sometime used for Sweden/IVL/SEI/CEH
          !! include 'EXTRA_LU_Outputs.inc'

       !=======================
        end do LULOOP
       !=======================
       !=======================

        ! Convert from Vg*f to f for grid-average:
        !where ( Vg_ref > 1.0e-6 )
        where ( Sub(0)%Vg_ref > 1.0e-6 )
           Sub(0)%StoFrac = Sub(0)%StoFrac/Sub(0)%Vg_ref
        end where

         call Calc_StoFlux(nFlux, iL_fluxes(1:nFlux), debug_flag )



        if ( Sumland > 0.01 ) then
            gradient_fac(:) = Vg_ratio(:) / Sumland
        else
            gradient_fac(:) = sea_ratio(:)
        end if

        if ( dbghh ) then
            call datewrite("DEP VGR snow_flag Vg", (/ Grid%sdepth, & 
                (100.0*Sub(0)%Vg_Ref(n), n = 1, min(4,NDRYDEP_GASES )) , &
                (100.0*Sub(iL)%Vg_Ref(n), n = 1, min(4,NDRYDEP_GASES )) /) )
        end if


!-- loop through all affected advected species to calculate changes in
!   concentration (xn_adv), the conc. ratios (cfac), and deposition 

    do ncalc = 1, NDRYDEP_CALC

        vg_fac (ncalc) = 1.0 - exp ( -Sub(0)%Vg_Ref(ncalc) * dtz ) 

    end do ! n


    GASLOOP2 :  do n = 1, NDRYDEP_ADV 
         nadv    = DDepMap(n)%ind
         ntot    = NSPEC_SHL + DDepMap(n)%ind
         ncalc   = DDepMap(n)%calc


         if ( vg_set(n) ) then

             DepLoss(nadv) =   & ! Use directly set Vg
                 ( 1.0 - exp ( -DDepMap(n)%vg * dtz ) ) * xn_2d( ntot,K2)
             cfac(nadv, i,j) = 1.0   ! Crude, for now.
  
         else
           if ( ntot >= FIRST_SEMIVOL .and. ntot <= LAST_SEMIVOL ) THEN
              ! Assuming dry deposition of particulate part of
              ! semi-volatile components as PMfS and the gaseous part as
              ! specified in GenIn.species.

             DepLoss(nadv) =  &
              Fgas(ntot,K2)*vg_fac( ncalc ) * xn_2d(ntot,K2) + &
              Fpart(ntot,K2)*vg_fac( CDDEP_PMfS ) * xn_2d(ntot,K2)

              cfac(nadv, i,j) = Fgas(ntot,K2)*gradient_fac(ncalc) + &
                   Fpart(ntot,K2)*gradient_fac( CDDEP_PMfS )
            else
               DepLoss(nadv) =   vg_fac( ncalc )  * xn_2d( ntot,K2)
               cfac(nadv, i,j) = gradient_fac( ncalc )
            end if
         end if

         if ( DepLoss(nadv) < 0.0 .or. DepLoss(nadv)>xn_2d(ntot,K2) ) then
            print "(a,2i4,a,es12.4,2f8.4,9es11.4)", "NEGXN ", ntot, ncalc, &
              trim(species(ntot)%name), xn_2d(ntot,K2), &
                 Fgas(ntot,K2), Fpart(ntot,K2), DepLoss(nadv), vg_fac(ncalc)
            call CheckStop("NEGXN DEPLOSS" )
         end if


         if ( ntot == O3 ) then 

           o3_45m = xn_2d(O3,K2)*surf_ppb !store for consistency of SPOD outputs
           Grid%surf_o3_ppb  = o3_45m * gradient_fac( ncalc )
           Grid%surf_o3_ppb1 = ( xn_2d(O3,K2) - Deploss(nadv)) * & ! after loss
                gradient_fac( ncalc )*surf_ppb

           if( dbghh ) then
              lossfrac = ( 1 - DepLoss(nadv)/xn_2d( ntot,K2))
              call datewrite("O3_ppb_ratios ", n, (/ o3_45m, &
                 Grid%surf_o3_ppb, Grid%surf_o3_ppb1, &
                 lossfrac, gradient_fac(ncalc), L%StoFrac(ntot) /) )
           end if

         end if

        xn_2d( ntot,K2) = xn_2d( ntot,K2) - DepLoss(nadv)


        if ( ntot == FLUX_TOT ) then

           ! fraction by which xn is reduced - safety measure:
              if( xn_2d(ntot,K2)  > 1.0e-30 ) then
               lossfrac = ( 1 - DepLoss(nadv)/(DepLoss(nadv)+xn_2d( ntot,K2)))
              end if
              if ( DEBUG%DRYDEP .and. lossfrac < 0.1 ) then
                print *, dtxt//"LOSSFRACING ", nadv, (/ 1.0*iL, &
                  Sub(0)%Vg_Ref(n), DepLoss(nadv), vg_fac(ncalc), lossfrac /)
                call CheckStop( lossfrac < 0.1, "ERROR: LOSSFRAC " )
              end if

              if ( DEBUG%AOT .and. debug_flag ) then !FEB2013 testing
                call datewrite(dtxt//"CHVEGX ", me, &
                  (/ xn_2d(FLUX_TOT,K2)*surf_ppb, c_hveg*surf_ppb,&
                   c_hveg3m * surf_ppb, 100*Vg_ref(n), 100*Vg_3m(n) /) )
              end if
        end if


      !.. ecosystem specific deposition - translate from calc to adv 
      !  and normalise

         IIL_LOOP : do iiL = 1, nlu
            iL      = iL_used(iiL)

            if ( vg_set(n) )  then
               fluxfrac_adv(nadv,iL) = Sub(iL)%coverage  ! Since all vg_set equal
            else
               Vg_scale = Sub(iL)%Vg_Ref(ncalc)/ Sub(0)%Vg_Ref(ncalc)
               fluxfrac_adv(nadv,iL) = Sub(iL)%coverage*Vg_scale
            end if


          !=======================
          ! The fraction going to the stomata = g_sto/g_sur = g_sto * R_sur.
          ! Vg*nmole_o3 is the instantaneous deposition flux of ozone, and
          ! the actual amount "deposited" is obtained from DepLoss(O3) using
          ! fluxfrac as obtained above.

          ! Now, DepLoss is loss of molecules/cm3 over time-period dt_advec
          ! and depth z_bnd over 1m2 of grid. For sto fluxes we need to
          ! find values over 1m2 of vegeation (regardless of how much veg
          ! is in grid, so we don't need cover. Instead:


            if ( dbghh .and. iL == 1 ) then ! SO2, CF

               if ( vg_set(n) )  then
                 write(6,"(a,3i3,3f12.3)") "FLUXSET  ", iiL, iL, nadv, &
                     100*DDepMap(n)%vg, Sub(iL)%coverage, fluxfrac_adv(nadv,iL)
               else
                 write(6,"(a,3i3,f8.5,5f8.3)") "FLUXFRAC ", iiL, iL, nadv, &
                  Sub(iL)%coverage, &
                  100*Sub(0)%Vg_Ref(ncalc), &  ! GRID
                  100*Sub(iL)%Vg_Ref(ncalc), & ! Mosaic VgRef &
                  100*Sub(iL)%coverage*Sub(iL)%Vg_Ref(ncalc), & 
                   fluxfrac_adv(nadv,iL)
               end if
            end if !SO2 CF
         end do IIL_LOOP
          

      !..accumulated dry deposition per grid square and summed over the whole
      !  domain

!         totddep( nadv ) = totddep (nadv) + Deploss(nadv) * convfac


        if ( dbghh) then
          if ( vg_set(n) ) then
              write(*, "(a,2i4,f8.3)") "DEBUG DryDep SET ", &
                   n,nadv, DDepMap(n)%vg
          else
              if( n == 1 ) & ! O3
              call datewrite( "DEBUG DDEPxnd: "// trim(species(ntot)%name), &
                n, (/ real(nadv),real(ncalc), gradient_fac( ncalc),&
                   xn_2d(ntot,K2), vg_fac(ncalc) /) )
          end if
        end if

        if ( DEBUG%AOT .and. debug_flag .and. ntot == FLUX_TOT  ) then
              write(*, "(a,3i3,i5,i3,2f9.4,f7.3)") &
               dtxt//"AOTCHXN ", imm, idd, ihh, current_date%seconds, &
                   iL, xn_2d(FLUX_TOT,K2)*surf_ppb, &
                    (xn_2d( FLUX_TOT,K2) + DepLoss(nadv) )*surf_ppb, &
                     gradient_fac( ncalc)
        end if
       end do GASLOOP2 ! n


      convfac =  convfac/amk(K2)

    !  DryDep Budget terms  (do not include values on outer frame)
      if(.not.(i<li0.or.i>li1.or.j<lj0.or.j>lj1))then
         
         do n = 1, NDRYDEP_ADV
            nadv    = DDepMap(n)%ind
            totddep( nadv ) = totddep (nadv) + DepLoss(nadv)*convfac
         end do
      end if

       convfac2 = convfac * xm2(i,j) * inv_gridarea

      !.. Add DepLoss to budgets if needed:

       call Add_MosaicOutput(debug_flag,i,j,convfac2,&
           DepAdv2Calc, fluxfrac_adv, Deploss ) 



  !    !----------------------------------------------------------------
  !    ! HUGE OUTPUTS. Not for routine use ! 
  !     if (SPOD_OUT ) then   ! Extra outputs for ICP folkks
  !        if ( first_spod .and. nlocal_sites > 0  ) then
  !    P = IO_SPOD + me
  !    write(fname,"(a,i2.2)") "OutputSPOD", me
  !    open(P, file=fname)
  !          write(P,"(a25,a5)",advance="no") adjustl("name"), "code"
  !          write(P,"(3a4)",advance="no")   "SGS", "EGS", "jd"
  !          write(P,"(3a3)",advance="no")   "mm", "dd", "hh"
  !          write(P,"(a7)",advance="no")  "cover"
  !          write(P,"(2a7)",advance="no")  "t2C", "ustar"
  !          write(P,"(3a7)",advance="no") "G-rh2m","L-rh","L-vpd"
  !          write(P,"(2a8)",advance="no")  "PARsun", "PARshd"
  !          write(P,"(2a7)",advance="no") "G45-O3","G3-O3"
  !         !if ( LandType(iL)%flux_wanted ) then
  !          write(P,"(a7)",advance="no") "CO3ppb"
  !          write(P,"(2a10)",advance="no") "CO3nmole","FstO3"
  !          ! Some params have value -999 unless sun shining...
  !          write(P,"(5a7)",advance="no") &
  !                "fphen", "f_L","f_T","f_D", "fenv" 
  !          write(P,"(a7,2a10)",advance="no")  "fsun", "gsto", "gsun"
  !          write(P,*)
  !         !end if
  !          first_spod = .false.
  !        end if ! first spod

  !      if( iss==0 .and. nlocal_sites > 0  )then
  !       do  n = 1, nlocal_sites
  !          if( i == site_x(n) .and. j == site_y(n) ) then 
  !             nglob = site_gn(n)  ! number in global list
!        if( iss==0 .and. debug_proc.and.  i==debug_li .and. j==debug_lj)then
  !     
  !       do iiL = 1, nlu
  !          iL  = iL_used(iiL)
  !          L   = Sub(iL)
  !        !print *, "SPOD_OUT:"//trim(fname), me,nlocal_sites, iL, imm, idd, ihh
  !          write(P,"(a25,a8)",advance="no") adjustl(site_name(nglob)), &
  !             adjustl(LandDefs(iL)%code)
  !          write(P,"(3i4)",advance="no")   L%SGS, L%EGS, daynumber
  !          write(P,"(3i3)",advance="no")   imm, idd, ihh
  !          write(P,"(f7.3)",advance="no")  L%coverage
  !          write(P,"(2f7.2)",advance="no") L%t2C,L%ustar
  !          write(P,"(3f7.3)",advance="no") Grid%rh2m,L%rh,L%vpd
  !          write(P,"(2f8.1)",advance="no") L%PARsun,L%PARshade
  !          write(P,"(2f7.2)",advance="no") o3_45m,Grid%surf_o3_ppb
  !         if ( LandType(iL)%flux_wanted ) then
  !          write(P,"(f7.2)",advance="no") L%cano3_ppb
  !          write(P,"(2es10.3)",advance="no") L%cano3_nmole,L%FstO3
  !          ! Some params have value -999 unless sun shining...
  !          write(P,"(5f7.2)",advance="no") &
  !               L%f_phen,max(-1.,L%f_light),max(-1.,L%f_temp),max(-1.,L%f_vpd), L%f_env 
  !          write(P,"(f7.3,2es10.3)",advance="no") L%f_sun, L%g_sto, L%g_sun
  !         end if
  !          !write(P,"(a)",advance="yes") ";"
  !          write(P,*) 
  !       end do !iL
  !      end if
  !       end do !sites
  !      end if !iss
  !    end if !SPOD_OUT
  !   !SPOD----------------------------------------------------
 end subroutine drydep

end module DryDep_ml
