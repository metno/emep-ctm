! <DryDep_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module DryDep_mod

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

! FUTURE: BiDir functionality (Dave/Roy Wichink Kruit) using Roy's methods
  
use AeroConstants_mod,    only: AERO
use Aero_Vds_mod,         only: SettlingVelocity, GPF_Vds300, Wesely300
use BiDir_emep
use BiDir_module
use Biogenics_mod,        only: SoilNH3  ! for BiDir
use CheckStop_mod,        only: CheckStop, StopAll
use Chemfields_mod ,      only: cfac, so2nh3_24hr,Grid_snow 
use ChemDims_mod,         only: NSPEC_ADV, NSPEC_SHL,NDRYDEP_ADV
use ChemSpecs_mod            ! several species needed
use Config_module,        only: dt_advec,PT, K2=> KMAX_MID, NPROC, &
                              USES, USE_SOILNOX, &
                              MasterProc, &
                              PPBINV, IOU_INST,&
                              KUPPER, NLANDUSEMAX
use Debug_module,         only: DEBUG, DEBUG_ECOSYSTEMS
use DerivedFields_mod,    only: d_2d, f_2d, VGtest_out_ix
use DO3SE_mod,            only: do3se
use EcoSystem_mod,        only: EcoSystemFrac, Is_EcoSystem,  &
                                 NDEF_ECOSYSTEMS, DEF_ECOSYSTEMS
use GasParticleCoeffs_mod         ! ... Init_GasCoeff, DRx, Rb_Cor, ...
use GridValues_mod ,      only: GRIDWIDTH_M,xmd,xm2, glat,dA,dB, glon, &
                             debug_proc, debug_li, debug_lj, i_fdom, j_fdom
use Io_Progs_mod,         only: datewrite
use Landuse_mod,          only: SetLandUse, Land_codes  & 
                             ,NLUMAX &  ! Max. no countries per grid
                             ,LandCover   ! Provides codes, SGS, LAI, etc,
use LandDefs_mod,         only: LandType, LandDefs, STUBBLE
use LocalVariables_mod,   only: Grid, L, iL & ! Grid and sub-scale Met/Veg data
                                ,NLOCDRYDEP_MAX ! Used to store Vg
use MassBudget_mod,       only: totddep
use MetFields_mod,        only: u_ref, rh2m, sst, tau, sdepth, &
                             SoilWater_deep, th,pzpbl
use MicroMet_mod,         only: AerRes, Wind_at_h
use MosaicOutputs_mod,    only: Add_MosaicOutput, MMC_RH
use Par_mod,              only: limax,ljmax, me,li0,li1,lj0,lj1
use PhysicalConstants_mod, only: ATWAIR,PI,KARMAN,GRAV,RGAS_KG,CP,AVOG,NMOLE_M3
use Rb_mod,               only: Rb_gas
use Rsurface_mod,         only: Rsurface, Rinc
use ZchemData_mod,        only: xn_2d,M, Fpart, Fgas
use Sites_mod,            only: nlocal_sites, site_x, site_y, &
                                  site_name, site_gn
use SmallUtils_mod,       only:  find_index
use SoilWater_mod,        only: fSW !  =1.0 unless set by Met_mod
use StoFlux_mod,          only: unit_flux, &! = sto. flux per m2
                                 lai_flux,  &! = lai * unit_flux
                                 Setup_StoFlux, Calc_StoFlux  ! subs
use SubMet_mod,           only: Sub
use TimeDate_mod,         only: daynumber, current_date

implicit none
private

public  :: DryDep, init_drydep

!integer, private, save :: P != IO_SPOD + me
! Maps from adv index to one of calc indices

logical, private, save :: my_first_call = .true.
character(len=30),private, save :: errmsg = "ok"

! WE NEED A FLUX_CDDEP, FLUX_ADV FOR OZONE (set to 1 for non-ozone models)

integer, public, parameter :: FLUX_ADV   = IXADV_O3
integer, public, parameter :: FLUX_TOT   = O3

! WE ALSO NEED NO3_f and NH4_f for deposition.
! (set to one for non-no3/nh4 models)

integer, public, parameter :: pNO3  = NO3_f
integer, public, parameter :: pNH4  = NH4_f

!logical, public, parameter :: COMPENSATION_PT = .false. 

logical, public, dimension(NDRYDEP_ADV), save :: vg_set 

!  A type container for big-leaf (bulk) resistances, used by esx
!  A little confusing still, bt Vg_ref etc are part of Sub(iL)
!   SKIP  ,rb_leaf      & ! Quasi-boundary layer rsis.
type, public :: BL_t
  real :: &
    Rsur    & ! Surface Resistance (s/m) 
   ,RsurX   & ! for testing Surface Resistance (s/m) 
   ,Vd      & ! Dep. vel. from zRef to surface
   ,Vd0     & ! Dep. vel. lowest z to surface
   ,VdX     & ! testing with less Rinc influence
   ,Rg      & ! R to ground, e.g. soil
   ,rext    & !  = rextO for O3, otherwise scaled SO2-O3
   ,Rb      & !
   ,Rns       !
  real ::  Gsto  !
  real ::      & ! for Aerosols (but simplest here anyway)
    Vds      & ! Vds term for fine, bulk (m/s)
   ,Vs         ! settling velocity (m/s). Not really BL, but same dims
end type BL_t
type(BL_t), public, save, allocatable, dimension(:) :: BL

real, allocatable, dimension(:), private :: &
    gradient_fac & ! Ratio of conc. at zref (ca. 50m) and 3m
   ,vg_fac       & ! Loss factor due to dry dep.
   ,Vg_ref       & ! Vg at ref ht.
   ,Vg_eff       & ! effective Vg at ref ht.
   ,Vg_3m        & ! Vg at  3m
   ,Vg_ratio     & ! Ratio Vg_ref/Vg_3m = ratio C(3m)/C(ref), over land
   ,sea_ratio    & ! Ratio Vg_ref/Vg_3m = ratio C(3m)/C(ref), over sea
   ,Gsto         & ! Stomatal conductance (big-leaf), only for gases
   ,eff_fac        ! Factor for computing effective resistance

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine init_DryDep()

  integer, save ::  old_daynumber = -99
  integer :: i, j, lc, ilc, nlc, is, iEco   ! for EcoSystem stuff
  logical :: debug_flag  ! for EcoSystem stuff
  real    :: coverage    ! for EcoSystem stuff
  character(len=*), parameter :: dtxt='iniDDep:'

  if ( my_first_call ) then 

     call GetDepMapping()    ! creates DDspec, DDmapping
     call InitGasCoeffs()    ! allocate and set DDspec  
     call InitParticleCoeffs()
     if(MasterProc) write(*,*) dtxt//" GET DEP", nddep

     allocate(BL(nddep))
     allocate(gradient_fac(nddep), vg_fac(nddep), Vg_ref(nddep), &
         Vg_eff(nddep), Vg_3m(nddep), Vg_ratio(nddep), sea_ratio(nddep), Gsto(nddep) ,eff_fac(nddep))

     call CheckStop( NLOCDRYDEP_MAX < nddep, &
        "Need to increase size of NLOCDRYDEP_MAX" )

     !Need to re-implement one day
     !A2018 nadv = 0
     !A2018 do n = 1, nddep  ! A2018  NDRYDEP_ADV  
         !nadv       = max( DDepMap(n)%ind, nadv )  ! Looking for highest IXADV
     !A2018     vg_set(n)  = ( DDepMap(n)%calc == CDDEP_SET ) ! for set vg
     !A2018 end do

     my_first_call = .false.
     if( MasterProc  .and. DEBUG%DRYDEP) write(*,*) "INIT_DRYDEP day ", &
           daynumber, old_daynumber

!=============================================================================
! From EcoSystems, but caused some circularity problem
! use EcoSystem_mod, only :: EcoSystemFrac, Is_EcoSystem

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

  end subroutine init_DryDep

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine DryDep(i,j)
    integer, intent(in):: i,j
    real, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost
    logical, save      :: dbg, dbghh, dbgBD
    integer icmp, ispec, iiL, nlu, nadv, nFlux  ! help indexes
    integer :: imm, idd, ihh, iss     ! date
    integer :: ntot !index of adv species in xn_2d array

    real :: no2fac  ! Reduces Vg for NO2 in ratio (NO2-4ppb)/NO2

    real convfac,  & ! rescaling to different units
         convfac2, & ! rescaling to different units
         lossfrac,  & !  If needed in My_DryDep - not used now.
         dtz         ! scaling factor for veff ( = dt/z, where dt=timestep and 
                     ! z = height of layer)

    integer :: nae

    real, save :: inv_gridarea  ! inverse of grid area, m2

    real ::  Sumcover, Sumland   ! Land-coverage
    logical :: debug_flag        ! set true when i,j match DEBUG_i, DEBUG_j
    character(len=*), parameter :: dtxt='DryDep:' ! debug label
    real :: Vg_scale

    real, dimension(NSPEC_ADV ,NLANDUSEMAX):: fluxfrac_adv
    integer, dimension(NLUMAX)  :: iL_used, iL_fluxes
     !BIDIR SKIP    real :: wet, dry         ! Fractions
    real :: fsnow            ! snow fraction for one landuse
    real :: Vds              ! Aerosol near-surface deposition rate (m/s)
    real :: no3nh4ratio      ! Crude NH4/NO3 for Vds ammonium 

    real :: c_hveg, Ra_diff, surf_ppb  ! for O3 fluxes and Fst where needed
    real :: c_hveg3m, o3_45m  ! TESTS ONLY
    character(len=20), save :: fname
    integer :: nglob
    logical :: first_ddep = .true.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Extra outputs sometime used. Important that this line is kept at the
    !! end of the variable definitions and the start of real code - allows
    !! both in .inc file
    !! Uncomment and make .inc file as required
    ! temporary for POD/SPOD
    !    logical, parameter :: SPOD_OUT = .false. !MAKES HUGE FILES!
    !    logical, save      :: first_spod = .true.
    !   include 'EXTRA_LU_Setup.inc'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! first calculate the 3m deposition velocity following the same
    ! procedure as in the lagmod. second the flux and the accumulated 
    ! deposition is calculated.
    !
    ! effective dry deposition velocity applied to the model concentration
    ! at the top of the constant flux layer, zdep 
    ! Dry deposion rates are specified in subroutine readpar
    !



   ! - Set up debugging stuff first. ---------------------------!
   ! If location matches debug i,j value, set debug_flag. Also passed
   ! to Rsurface_mod

    imm      =    current_date%month            ! for debugging
    idd      =    current_date%day              ! for debugging
    ihh      =    current_date%hour             ! for debugging
    iss      =    current_date%seconds          ! for debugging

    debug_flag= ( debug_proc .and. i == debug_li .and. j == debug_lj) 
    dbg       =  DEBUG%DRYDEP .and. debug_flag 
    dbghh     =  dbg .and. iss == 0 
    dbgBD     =  DEBUG%BIDIR .and. debug_flag .and. iss == 0
    if (dbgBD) print *, "BIDIR TEST ", me, debug_flag



     inv_gridarea = 1.0/(GRIDWIDTH_M*GRIDWIDTH_M) 

   ! -----------------------------------------------------------------!
   !.and conversion factor,  convfac (( ps-pt)/grav... )  ===> 
   !      pressure in kg m-1 s-2
   ! 
    convfac = (dA(K2) + dB(K2)*Grid%psurf)&!dP
               *xmd(i,j)/(ATWAIR*GRAV*inv_gridarea)

   ! -----------------------------------------------------------------!
   ! convert molecules/cm3 to ppb for surface:
    surf_ppb   = PPBINV /M(K2)
    if (DEBUG%AOT.and.debug_flag)write(*,"(a,es12.4)")dtxt//"PPB",surf_ppb

     inv_gridarea = 1.0/(GRIDWIDTH_M*GRIDWIDTH_M) 

   ! -----------------------------------------------------------------!

   !.and factor,  kg_air_ij (( ps-pt)/grav... )  ===> 
   !      pressure in kg m-1 s
   !      used for converting from mixing ratio to kg

   !   kg_air_ij = (ps(i,j,1) - PT)*carea(K2) = dP*dx**2/g
   ! -----------------------------------------------------------------!

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
      write(*,"(a,2i4,L2,9es12.4)") dtxt//" CONCS SO2,NH3,O3,NO2 (ppb),f ", &
       i,j, USE_SOILNOX,  xn_2d(SO2,K2)*surf_ppb, xn_2d(NH3,K2)*surf_ppb, &
          xn_2d(O3,K2)*surf_ppb, xn_2d(NO2,K2)*surf_ppb, no2fac
    end if

    call GasCoeffs(Grid%t2)             ! Sets DDdefs coeffs.

    call ParticleCoeffs(Grid%t2,Grid%rho_s,debug_flag=debug_flag)

    call Setup_StoFlux( daynumber )

! - can set settling velcoty here since not landuse dependent

    do icmp = 1, nddep

      if ( DDspec(icmp)%is_gas ) CYCLE

      BL(icmp)%Vs = SettlingVelocity( Grid%t2, Grid%rho_ref, &
             DDspec(icmp)%sigma, DDspec(icmp)%DpgV, DDspec(icmp)%rho_p )

      if ( dbghh ) then
        call datewrite(dtxt//" VS"//DDspec(icmp)%name, icmp, &
           [ DDspec(icmp)%DpgV, DDspec(icmp)%sigma,  DDspec(icmp)%rho_p, &
              BL(icmp)%Vs ] )
      end if
    end do


!PW_____________________________________________________________________
!Make effective resistance using formula:
! R_effective_i = (Ra_X_i+Rb_i+Rsur_i) * sum_j( coverage_j*(Ra_ref_j+Rb_j+Rsur_j)/(Ra_X_j+Rb_j+Rsur_j) )
!need to make sum_j() in advance

    nlu = LandCover(i,j)%ncodes
    nFlux = 0           ! Numberof LC needing flux calcs
    eff_fac = 0.0
    LULOOP_PRE: do iiL= 1, nlu
      iL      = LandCover(i,j)%codes(iiL)

      iL_used (iiL) = iL    ! for eco dep

      if ( LandType(iL)%flux_wanted ) then
        nFlux = nFlux + 1
        iL_fluxes (nFlux) = iL    ! for eco dep
      end if

      Sub(iL)%f_phen = LandCover(i,j)%fphen(iiL) ! for POD_OUT
      Sub(iL)%f_sun  = 0.0 ! for SPOD_OUT
      Sub(iL)%g_sun  = 0.0 ! for SPOD_OUT
      Sub(iL)%g_sto  = 0.0 ! for SPOD_OUT
      Sub(iL)%f_temp = 0.0 ! for SPOD_OUT. Can 
      Sub(iL)%f_vpd  = 0.0 ! for SPOD_OUT

      Sub(iL)%SGS = LandCover(i,j)%SGS(iiL)   !used for AOT CHECK?
      Sub(iL)%EGS = LandCover(i,j)%EGS(iiL)

      L = Sub(iL)    ! ! Assign e.g. Sub(iL)ustar to ustar


      call Rb_gas(L%is_water, L%ustar, L%z0, BL(:)%Rb)

      call Rsurface(i,j,BL(:)%Gsto,BL(:)%Rsur,errmsg,debug_flag,fsnow)

      do icmp = 1, nddep    !DSQUERY - Check aerosol usage!

        ! L%Ra_ref + BL(icmp)%Rb + BL(icmp)%Rsur is res. from mid cell and down
        ! L%Ra_X   + BL(icmp)%Rb + BL(icmp)%Rsur is res. from 50 meter and down

        eff_fac(icmp) = eff_fac(icmp) &
            + L%coverage*(L%Ra_ref + BL(icmp)%Rb + BL(icmp)%Rsur)&
                         /(L%Ra_X + BL(icmp)%Rb + BL(icmp)%Rsur)

      end do ! dep species loop
      !=======================

    end do LULOOP_PRE

! Vg_eff(icmp) = 1. / ( (L%Ra_X + BL(icmp)%Rb + BL(icmp)%Rsur) * eff_fac(icmp)) 
!PW_____________________________________________________________________
        

    !/ And start the sub-grid stuff over different landuse (iL)

    nlu = LandCover(i,j)%ncodes
    nFlux = 0           ! Numberof LC needing flux calcs
    LULOOP: do iiL= 1, nlu

      iL            = LandCover(i,j)%codes(iiL)
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
         write(6,"(a,3i3,f6.1,2i4,3f7.3,i4,9f8.3)") dtxt//"DVEG: ", &
             nlu,iiL, iL, glat(i,j), L%SGS, L%EGS, &
            L%coverage, L%LAI, L%hveg,daynumber, &
            Grid%sdepth, fSW(i,j),L%fSW,L%t2C

         write(6,"(a,i4,3f7.2,7es10.2)") dtxt//"DMET SUB", &
           iL, Grid%ustar, L%ustar, L%rh,  Grid%invL, &
             L%invL, L%Ra_ref, L%Ra_3m
      end if


      call Rb_gas(L%is_water, L%ustar, L%z0, BL(:)%Rb)

      call Rsurface(i,j,BL(:)%Gsto,BL(:)%Rsur,errmsg,debug_flag,fsnow)

      if(dbghh) call datewrite(dtxt//"STOFRAC "//LandDefs(iL)%name, &
              iL, [  BL(idcmpO3)%Gsto, BL(idcmpO3)%Rsur, &
                    BL(idcmpO3)%Gsto*BL(idcmpO3)%Rsur ] )

        !Sub(iL)%g_sto = L%g_sto   ! needed elsewhere
        !Sub(iL)%g_sun = L%g_sun
      Sub(iL) = L  !Resets Sub with new L values frmom Rsurface

      Grid_snow(i,j) = Grid_snow(i,j) +  L%coverage * fsnow


       !/... add to grid-average Vg:

!BIDIR SKIP         wet =   Grid%wetarea  ! QUERY Now used for BIDIR
!BIDIR SKIP         dry =   1.0 - wet     !  "  "

      CMPLOOP: do icmp = 1, nddep

        ! ================================================
        if ( DDspec(icmp)%is_gas  ) then

           Vg_ref(icmp) = 1. / ( L%Ra_ref + BL(icmp)%Rb + BL(icmp)%Rsur ) 
           Vg_eff(icmp) = 1. / ( (L%Ra_X + BL(icmp)%Rb + BL(icmp)%Rsur) * eff_fac(icmp)) 
           Vg_3m(icmp)  = 1. / ( L%Ra_3m  + BL(icmp)%Rb + BL(icmp)%Rsur ) 

           if( dbg .and. first_ddep .and. icmp==idcmpNO2 ) then
             associate ( b=>BL(icmp) )
               write(*,'(a,2i4,9es10.3)')'DBGXNO2 :',icmp, iL,L%Ra_ref, b%Rb,&
                      b%Rsur, b%Gsto, Vg_ref(icmp), no2fac, 4.0e-11*xn_2d(NO2,K2)
             end associate
           end if

          ! specials, NH3 and NO2:

           if ( USES%BIDIR .and. icmp == idcmpNH3 ) then

              call StopAll(dtxt//'NOT IMPLEMENTED YET')
           
           ! Surrogate for NO2 compensation point approach, 
           ! assuming c.p.=4 ppb (ca. 1.0e11 #/cm3):        
           ! Note, xn_2d has no2 in #/cm-3

            else if ( .not. USE_SOILNOX .and.  icmp == idcmpNO2 ) then

              if( dbg .and. first_ddep .and. icmp==idcmpNO2 ) &
                  write(*,*) 'DBGXNO2 no2fac TRIGGERED', no2fac, 4.0e-11*xn_2d(NO2,K2)
  
              Vg_eff(icmp) = Vg_eff(icmp) * no2fac
              Vg_ref(icmp) = Vg_ref(icmp) * no2fac
              Vg_3m(icmp)  = Vg_3m(icmp)  * no2fac

            end if ! specials, NH3, NO2

            !QUERY - do we need Gsur for anything now?!
            ! StoFrac should end up with weighted mean. Start with Vg*f
   
            Sub(iL)%Gsur(icmp) = 1.0/BL(icmp)%Rsur ! Note iL, not iiL 
            Sub(iL)%Gsto(icmp)  = BL(icmp)%Gsto    ! Note iL, not iiL 

            Sub(0)%Gsur(icmp)  =  Sub(0)%Gsur(icmp) + L%coverage / BL(icmp)%Rsur
            Sub(0)%Gsto(icmp)   =  Sub(0)%Gsto(icmp)  + L%coverage * BL(icmp)%Gsto
            if( dbghh.and.icmp==2 ) call datewrite("CmpSto", iL, &
                       (/ Sub(iL)%Gsto(icmp) / Sub(iL)%Gsur(icmp) /) )


        else ! particles
         ! ================================================

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

           if (icmp==idcmpPMfNO3 .and. L%invL<0.0 ) then 
             Vds = Vds * 3.0 ! for nitrate-like
           else if (icmp==idcmpPMfNH4 .and. L%invL<0.0 ) then 
             !Vds = Vds* [ ( 1 - no3ratio) +  3 * no3nh4ratio ]
             Vds = Vds * (1 + 2 * no3nh4ratio) 
           end if

         ! Use non-electrical-analogy version of Venkatram+Pleim (AE,1999)
         ! ACP70

           Vg_ref(icmp) = BL(icmp)%Vs/ &
                    ( 1.0 - exp( -( L%Ra_ref + 1.0/Vds)* BL(icmp)%Vs))
           Vg_3m (icmp) = BL(icmp)%Vs/ &
                    ( 1.0 - exp( -( L%Ra_3m  + 1.0/Vds)* BL(icmp)%Vs))

           Vg_eff(icmp) = Vg_ref(icmp) !PW

           if ( DEBUG%VDS ) then
             if ( debug_flag .or. &
                 (Vg_3m(icmp)>0.50 .or. Vg_ref(icmp)>0.50 )) then
               write(*,"(a,5i3,2i4,2f7.3,f8.2,20f7.2)") &
                 dtxt//"AEROCHECK:"//trim(DDspec(icmp)%name), &
                 imm, idd, ihh, icmp, iL, i_fdom(i), j_fdom(j),&
                 L%ustar,  L%invL, pzpbl(i,j), &
                 Grid%ustar, Grid%Hd,  100*Vds, &
                 100*BL(icmp)%Vs, 100*Vg_ref(icmp), 100*Vg_3m (icmp) &
                  , Grid%t2, Grid%rho_ref 
               write(*,"(a,i4,3es10.3)") dtxt//"VDS ", icmp, Vds, Vg_ref(icmp)
               
               call CheckStop((Vg_3m(icmp)>0.50 .or. Vg_ref(icmp)>0.50 ),&
                                  dtxt//"AEROSTOP")
             end if
           end if
        end if ! Gases?

        Sub(0)%Vg_eff(icmp) = Sub(0)%Vg_eff(icmp) + L%coverage * Vg_eff(icmp)
        Sub(0)%Vg_ref(icmp) = Sub(0)%Vg_ref(icmp) + L%coverage * Vg_ref(icmp)
        Sub(0)%Vg_3m(icmp)  = Sub(0)%Vg_3m(icmp)  + L%coverage * Vg_3m(icmp)
        Sub(iL)%Vg_ref(icmp) = Vg_ref(icmp)
        Sub(iL)%Vg_3m(icmp) = Vg_3m(icmp)

      end do CMPLOOP ! chemical species loop

      Sumcover = Sumcover + L%coverage

     !/-- only grab gradients over land-areas

      if ( L%is_water ) then
         do icmp = 1, nddep
            if(USES%EFFECTIVE_RESISTANCE)then
               sea_ratio(icmp) =  Vg_eff(icmp)/Vg_3m(icmp)
            else
               sea_ratio(icmp) =  Vg_ref(icmp)/Vg_3m(icmp)
            endif
         end do
      else
         Sumland = Sumland + L%coverage
         do icmp = 1, nddep
            if(USES%EFFECTIVE_RESISTANCE)then
               Vg_ratio(icmp) =  Vg_ratio(icmp) &
                                + L%coverage * Vg_eff(icmp)/Vg_3m(icmp)
            else
               Vg_ratio(icmp) =  Vg_ratio(icmp)&
                                 + L%coverage * Vg_ref(icmp)/Vg_3m(icmp)
            endif
         end do
      end if

      if ( dbghh ) then
         call CheckStop(  Sumland > 1.011, dtxt// "SUMLAND>1")
         do icmp = idcmpO3 , idcmpO3 !!! 1,NDRYDEP_GASES 
            call datewrite(dtxt//"DEPO3 ", iL, &
                (/ Vg_ref(icmp), Sub(iL)%Vg_ref(icmp) /) )

            call datewrite(dtxt//"VGA"//DDspec(icmp)%name, iL, [ L%coverage,&
              1.0*icmp, L%LAI,100*L%g_sto, L%Ra_ref, &
              DDspec(icmp)%Dx, BL(icmp)%Rb, min( 999.0,BL(icmp)%Rsur ) ] )

            call datewrite(dtxt//"VGB", iL, (/ L%coverage, 1.0*icmp,& 
             100*Vg_3m(icmp), 100*Vg_ref(icmp), Vg_ratio(icmp) /) )
         end do
      end if

    !=======================

      c_hveg   = -999. ! Just for printout when flux_wanted false
      c_hveg3m = -999.

      if (  LandType(iL)%flux_wanted ) then

       !PW NOT SURE HOW TO TREAT THIS using effective resistance (which is
       !    species dependent)

        !n = CDDEP_O3
        icmp = idcmpO3
        Ra_diff = L%Ra_ref - L%Ra_3m
        c_hveg3m = xn_2d(FLUX_TOT,K2)  &     ! #/cm3 units
                     * ( 1.0-Ra_diff*Vg_ref(icmp) )

      ! Flux = Vg_ref*c_ref = Vg_h * c_h = (c_ref-c_h)/Ra(z_ref,z_h)
      ! which gives:  c_h = c_ref * [ 1-Ra(z_ref,z_h)*Vg_ref ]
      ! Resistance Ra from z_ref to top of canopy:

        Ra_diff  = AerRes(max( L%hveg-L%d, STUBBLE) , Grid%z_ref-L%d,&
                    L%ustar,L%invL,KARMAN)

        c_hveg = xn_2d(FLUX_TOT,K2)  &     ! #/cm3 units
                     * ( 1.0-Ra_diff*Vg_ref(icmp) )

        if ( DEBUG%AOT .and. debug_flag .and. iL==1 ) then
           call datewrite(dtxt//"CHVEG ", iL, &
             (/  xn_2d(FLUX_TOT,K2)*surf_ppb, c_hveg*surf_ppb,&
                 c_hveg3m * surf_ppb, 100*Vg_ref(icmp), 100*Vg_3m(icmp), &
                 L%Ra_ref, (L%Ra_ref-L%Ra_3m), Ra_diff, &
                 BL(icmp)%Rb,BL(icmp)%Rsur /) )
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
        call datewrite(dtxt//" VGR fsnow Vg", (/ Grid%sdepth, & 
            (100.0*Sub(0)%Vg_Ref(icmp), icmp = 1, min(4,nddep )) , &
            (100.0*Sub(iL)%Vg_Ref(icmp), icmp = 1, min(4,nddep )) /) )
    end if


!-- loop through all affected advected species to calculate changes in
!   concentration (xn_adv), the conc. ratios (cfac), and deposition 

    do icmp = 1, nddep ! NDRYDEP_CALC

       if(USES%EFFECTIVE_RESISTANCE)then
          vg_fac (icmp) = 1.0 - exp ( -Sub(0)%Vg_Ref(icmp) * dtz ) 
       else
          vg_fac (icmp) = 1.0 - exp ( -Sub(0)%Vg_eff(icmp) * dtz ) 
       endif

    end do ! icmp

    if(VGtest_out_ix>0) d_2d(VGtest_out_ix,i,j,IOU_INST) = &
          Sub(0)%Vg_eff(1)/(1.E-20+Sub(0)%Vg_Ref(1))


  ! ===================================================================
   DDEPLOOP: do icmp = 1, nddep

     if ( vg_set(icmp) ) then
         call StopAll(dtxt//'NOT CODED')
         !A2018 DepLoss(nadv) =   & ! Use directly set Vg
         !A2018( 1.0 - exp ( -DDepMap(icmp)%vg * dtz ) ) * xn_2d( ntot,K2)
         cfac(nadv, i,j) = 1.0   ! Crude, for now.
     end if
 
     DCMPLOOP: do ispec = 1, size(DDmapping(icmp)%advspecs)  ! Real species now

        nadv = DDmapping(icmp)%advspecs(ispec)  ! Real species now
        ntot = NSPEC_SHL + nadv
          
        if ( ntot >= FIRST_SEMIVOL .and. ntot <= LAST_SEMIVOL ) THEN
           ! Assuming dry deposition of particulate part of
           ! semi-volatile components as PMf and the gaseous part as
           ! specified in GenIn.species.

          DepLoss(nadv) =  &
           Fgas(ntot,K2)*vg_fac( icmp ) * xn_2d(ntot,K2) + &
           Fpart(ntot,K2)*vg_fac( idcmpPMf ) * xn_2d(ntot,K2)

           cfac(nadv, i,j) = Fgas(ntot,K2)*gradient_fac(icmp) + &
                Fpart(ntot,K2)*gradient_fac( idcmpPMf )
        else
            DepLoss(nadv) =   vg_fac( icmp )  * xn_2d( ntot,K2)
            cfac(nadv, i,j) = gradient_fac( icmp )
        end if !SEMIVOL
        if( dbg .and. first_ddep ) then
           lossfrac = ( 1 - DepLoss(nadv)/(1+xn_2d( ntot,K2))) ! 1 avoids NaN
           write(*,'(a20,i4,4es10.3)') "DBGX "//trim(species(ntot)%name)&
             //' as:'// trim(DDspec(icmp)%name), icmp,  xn_2d(ntot,K2), &
             vg_fac(icmp), lossfrac, gradient_fac(icmp)
        end if

        if ( DepLoss(nadv) < 0.0 .or. DepLoss(nadv)>xn_2d(ntot,K2) ) then
          print "(a,2i4,a,es12.4,2f8.4,9es11.4)", dtxt//"NEGXN ", ntot, icmp,&
           trim(species(ntot)%name), xn_2d(ntot,K2), &
              Fgas(ntot,K2), Fpart(ntot,K2), DepLoss(nadv), vg_fac(icmp)
         call CheckStop("NEGXN DEPLOSS" )
        end if

        if ( ntot == O3 ) then 

          o3_45m = xn_2d(O3,K2)*surf_ppb !store for consistency of SPOD outputs
          Grid%surf_o3_ppb  = o3_45m * gradient_fac( icmp )
          Grid%surf_o3_ppb1 = ( xn_2d(O3,K2) - Deploss(nadv)) * & ! after loss
               gradient_fac( icmp )*surf_ppb

          if( dbghh ) then
             lossfrac = ( 1 - DepLoss(nadv)/xn_2d( ntot,K2))
             call datewrite("O3_ppb_ratios ", icmp, (/ o3_45m, &
                Grid%surf_o3_ppb, Grid%surf_o3_ppb1, &
                lossfrac, gradient_fac(icmp), L%StoFrac(ntot) /) )
          end if
        end if ! ntot==O3

        xn_2d( ntot,K2) = xn_2d( ntot,K2) - DepLoss(nadv)


        if ( ntot == FLUX_TOT ) then

          ! fraction by which xn is reduced - safety measure:
             if( xn_2d(ntot,K2)  > 1.0e-30 ) then
              lossfrac = ( 1 - DepLoss(nadv)/(DepLoss(nadv)+xn_2d( ntot,K2)))
             end if
             if ( DEBUG%DRYDEP .and. lossfrac < 0.1 ) then
               print *, dtxt//"LOSSFRACING ", nadv, (/ 1.0*iL, &
                 Sub(0)%Vg_Ref(icmp), DepLoss(nadv), vg_fac(icmp), lossfrac /)
               call CheckStop( lossfrac < 0.1, "ERROR: LOSSFRAC " )
             end if
  
             if ( DEBUG%AOT .and. debug_flag ) then !FEB2013 testing
               call datewrite(dtxt//"CHVEGX ", me, &
                 (/ xn_2d(FLUX_TOT,K2)*surf_ppb, c_hveg*surf_ppb,&
                  c_hveg3m * surf_ppb, 100*Vg_ref(icmp), 100*Vg_3m(icmp) /) )
             end if
        end if ! not FLUXTOT


        !.. ecosystem specific deposition - translate from calc to adv 
        !  and normalise
  
        IIL_LOOP : do iiL = 1, nlu
           iL      = iL_used(iiL)
  
           !if ( vg_set(icmp) )  then
           !      fluxfrac_adv(nadv,iL) = Sub(iL)%coverage  ! Since all vg_set equal
           !else
              Vg_scale = Sub(iL)%Vg_Ref(icmp)/ Sub(0)%Vg_Ref(icmp)
              fluxfrac_adv(nadv,iL) = Sub(iL)%coverage*Vg_scale
           !end if
  
  
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
  
              !if ( vg_set(icmp) )  then
              !  write(6,"(a,3i3,3f12.3)") "FLUXSET  ", iiL, iL, nadv, &
              !      100*DDepMap(icmp)%vg, Sub(iL)%coverage, fluxfrac_adv(nadv,iL)
              !else
              write(6,"(a,3i3,f8.5,5f8.3)") "FLUXFRAC "//species(ntot)%name,&
                 iiL, iL, ntot, &
                 Sub(iL)%coverage, 100*Sub(0)%Vg_Ref(icmp), &  ! GRID
                    100*Sub(iL)%Vg_Ref(icmp), & ! Mosaic VgRef &
                    100*Sub(iL)%coverage*Sub(iL)%Vg_Ref(icmp), & 
                     fluxfrac_adv(nadv,iL)
              !end if
           end if !SO2 CF
        end do   IIL_LOOP
          

      !..accumulated dry deposition per grid square and summed over the whole
      !  domain

!         totddep( nadv ) = totddep (nadv) + Deploss(nadv) * convfac


        if ( dbghh) then
          !if ( vg_set(icmp) ) then
          !  write(*, "(a,2i4,f8.3)") "DEBUG DryDep SET ", &
          !       n,nadv, DDepMap(icmp)%vg
          !else
          if( ntot == O3 ) & ! O3
            call datewrite( "DEBUG DDEPxnd: "// trim(species(ntot)%name), &
              icmp, (/ real(nadv),real(icmp), gradient_fac( icmp),&
                 xn_2d(ntot,K2), vg_fac(icmp) /) )
          !end if
        end if

        if ( DEBUG%AOT .and. debug_flag .and. ntot == FLUX_TOT  ) then
            write(*, "(a,3i3,i5,i3,2f9.4,f7.3)") &
             dtxt//"AOTCHXN ", imm, idd, ihh, current_date%seconds, &
                 iL, xn_2d(FLUX_TOT,K2)*surf_ppb, &
                  (xn_2d( FLUX_TOT,K2) + DepLoss(nadv) )*surf_ppb, &
                   gradient_fac( icmp)
        end if
     end do DCMPLOOP ! is
   end do DDEPLOOP ! n
  ! ===================================================================


   convfac =  convfac/M(K2)

    !  DryDep Budget terms  (do not include values on outer frame)
   if(.not.(i<li0.or.i>li1.or.j<lj0.or.j>lj1))then
      
     do icmp = 1, nddep
       do ispec = 1, size(DDmapping(icmp)%advspecs)  ! Real species now
          nadv = DDmapping(icmp)%advspecs(ispec)  ! Real species now
           totddep( nadv ) = totddep (nadv) + DepLoss(nadv)*convfac
       end do ! is
     end do ! icmp
   end if

   convfac2 = convfac * xm2(i,j) * inv_gridarea

  !.. Add DepLoss to budgets if needed:

   call Add_MosaicOutput(debug_flag,i,j,convfac2,&
            itot2DDspec, fluxfrac_adv, Deploss ) 

   ! SPOD outputs were put here

   ! 
   if( dbg ) first_ddep = .false.

 end subroutine drydep

end module DryDep_mod
