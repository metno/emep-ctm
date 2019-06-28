! <MetFields_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module MetFields_mod

  use Config_module,  only : USES,USE_WRF_MET_NAMES,NPROC, PBL, meteo, startdate, TopoFile
  use MPI_Groups_mod     , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_REAL8, MPI_INTEGER, MPI_LOGICAL, &
                                 MPI_COMM_CALC, MPI_COMM_WORLD, MPI_COMM_SUB, MPISTATUS, &
                                 IERROR, ME_MPI, NPROC_MPI, largeLIMAX,largeLJMAX, share, share_logical
  use Par_mod            , only : me
  use TimeDate_ExtraUtil_mod,only: date2string
  implicit none
  private


!----------------- basic met fields ----------------------------------!
!  Here we declare the meteorological fields used in the model        !
!
! Horizonal alignments....
! Placement of q(i,j), u(i,j), v(i,j) 
!
!    --------------- --------------- --------------- ---------------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u03   q13      u13    q23      u23     q33     u33    q43      u43 ... u(LIMAX,3)
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v12------- -----v22------- ------v32------ -----v42-------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u02   q12      u12    q22      u22     q32     u32    q42      u42
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v11------- -----v21------- ------v31------ -----v41-------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u01   q11      u11    q21      u21     q31     u31    q41      u41
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v10------- -----v20------- ------v30------ -----v40-------
!
  !---------------------------------------------------------------------!
  !
  !
  ! Vertical levels: z_mid,  z_bnd, sigma_mid, sigma_bnd
  !=============================================================================
  !*   "mid" and "bnd" are used as suffixes on z and sigma as shown in
  !*   the sketch below. "bnd" is the boundary between two layers and
  !*   "mid" the midddle of the layer. The numbering of layers starts
  !*   from 1 at the surface.
  !*
  !*
  !*
  !* ---------------------------
  !*
  !*
  !* - - - - - - - - - - - -     KMAX_MID -1
  !*
  !*
  !* --------------------------  KMAX_BDN-1       (z_bnd)   (sigma_bnd)
  !*
  !*
  !* - - - - - - - - -           KMAX_MID(old kmax2) = 20    (z_mid)   (sigma_mid)   (old z2)
  !*
  !* ------------------------    KMAX_BND = 21    (z_bnd)                 (old z1)
  !* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  surface \\\\\\\\\\\\\\\\
  !*
  ! RESULT FROM ONE TEST
  !Tested zm3d (=zm) which is sometimes used against
  !z_mid. Seems almost identical. Diff in exner functions maybe?
  !                   z_bnd        z_mid       zm3d 
  !     DEBUG_Z  2  1  16276.2359  15244.9023  15244.7131
  !     DEBUG_Z  2  2  14319.8682  13536.3354  13531.2981
  !     DEBUG_Z  2  3  12815.7707  12185.9090  12177.9740
  !     DEBUG_Z  2  4  11598.2202  11011.4342  11004.0286
  !     DEBUG_Z  2  5  10461.1309   9792.2354   9788.6193
  !     DEBUG_Z  2  6   9168.4701   8431.0622   8430.8546
  !     DEBUG_Z  2  7   7747.1534   7010.9878   7013.5860
  !     DEBUG_Z  2  8   6324.9278   5696.2857   5700.0226
  !     DEBUG_Z  2  9   5103.2253   4595.6436   4598.2893
  !     DEBUG_Z  2 10   4110.6562   3690.7686   3692.3391
  !     DEBUG_Z  2 11   3286.0770   2933.6493   2934.3594
  !     DEBUG_Z  2 12   2591.8972   2296.1632   2296.2730
  !     DEBUG_Z  2 13   2007.7924   1761.4943   1761.6936
  !     DEBUG_Z  2 14   1520.4434   1316.9528   1317.4565
  !     DEBUG_Z  2 15   1116.9477    950.9995    951.4678
  !     DEBUG_Z  2 16    787.3321    655.5900    655.8842
  !     DEBUG_Z  2 17    525.1633    424.5627    424.7443
  !     DEBUG_Z  2 18    324.8232    253.9725    254.0498
  !     DEBUG_Z  2 19    183.6019    137.4135    137.4284
  !     DEBUG_Z  2 20     91.3017     45.6234     45.6146


  !   Vertical level geopotential heights:

  real,public, save, allocatable,&
       dimension(:,:,:) :: z_bnd ! height of full layers. Updated each timestep
  real,public, save,allocatable, &
       dimension(:,:,:) :: z_mid ! height of half layers. Updated each METSTEP

  !   Two sets of Met. fields are read in, and a linear interpolation is made
  !   between these two points in time. NMET  == 2 (two points in time)
  !   note u_xmj, v_xmi are not "real" m/s wind speeds
  !   - they are actually divided by the mapping factor in the perpendicular direction).
  !
  real,target,public, save,allocatable, dimension(:,:,:,:) :: u_xmj 
  real,target,public, save,allocatable, dimension(:,:,:,:) :: v_xmi

  real,target,allocatable,public, dimension(:,:,:,:) :: q

  real,target,public, save,allocatable, dimension(:,:,:,:) :: &
       th      &  ! Potential teperature  ( deg. k )
!       ,q      &  ! Specific humidity
       ,roa    &  ! kg/m3
       ,cw_met    ! cloudwater from meteo file (not always defined!)
  real,target,public, save,allocatable, dimension(:,:,:,:) :: &
        EtaKz    &! vertical diffusivity in Eta coords
       ,SigmaKz  &! vertical diffusivity in sigma coords
       ,Etadot     &! vertical velocity, Eta coords, Pa/s
       ,Kz_met    ! vertical diffusivity in sigma coordinates from meteorology

  real,target,public, save,allocatable, dimension(:,:,:,:) :: rain! used only by wrf met
  real,target,public, save,allocatable, dimension(:,:) :: irainc! used only by wrf met
  real,target,public, save,allocatable, dimension(:,:) :: irainnc! used only by wrf met


  ! since pr,cc3d,cc3dmax,cnvuf,cnvdf used only for 1 time layer - define without NMET
  real,target,public, save,allocatable, dimension(:,:,:) :: &
        pr      & ! Precipitation in mm/s "passing" the layer
       ,cc3d    & ! 3-d cloud cover (cc3d),
       ,cc3dmax & ! and maximum for layers above a given layer
       ,lwc     & !liquid water content
  ! QUERY - should xksig be MID, not BND? Is it needed at all?
       ,Kz_m2s     ! estimated Kz, in intermediate sigma levels, m2/s

  real,target,public, save,allocatable, dimension(:,:,:) :: &
        cnvuf   & ! convective_updraft_flux (kg/s/m2)
       ,cnvdf    ! convective_downdraft_flux (kg/s/m2)


 ! We don't need to calculate u,v for RiB, Kz for all layer in future maybe
 ! Still, for safety  we let this extent to K=1 for now

  real,target,public, save,allocatable, dimension(:,:,:) :: &
        u_mid   & ! wind u-compnent, m/s (real, not projected)
       ,v_mid     ! wind v-compnent, m/s
  

 real,target,public,save,allocatable, dimension(:,:,:) :: &
       tau        ! surf. stress  N/m^2

! Surface fields, interpolated:
 real,target,public, save,allocatable, dimension(:,:,:) :: &
        ps        &! Surface pressure Pa
       ,t2_nwp    & ! Temp 2 m   deg. K
       ,pbl_nwp   & ! Planetary boundary layer height (m)
       ,fh        & ! surf.flux.sens.heat W/m^2
       ,fl        & ! latent heat flux W/m^2
!       ,tau       & ! surf. stress  N/m^2
  ! These fields only available for EMEP/PARLAM from 2002 on
       ,rh2m            & !  RH at 2m
       ,SoilWater_uppr  & !  Shallow  (Upper 7.2cm in PARLAM)
       ,SoilWater_deep  & !  Deep (Next 6x7cm in PARLAM), converted to relative value 
       ,sdepth          & !  Snowdepth, m
       ,ice_nwp         & ! QUERY why real?
       ,sst       &  ! SST Sea Surface Temprature- ONLY from 2002 in PARLAM
       ,ws_10m    ! wind speed 10m
 

 real,target,public, save,allocatable, dimension(:,:) :: &
     u_ref             & ! wind speed m/s at 45m (real, not projected)
    ,rho_surf          & ! Surface density
    ,surface_precip    & ! Surface precip mm/hr
    ,convective_precip & ! Convective precip mm/hr
    ,Tpot2m            & ! Potential temp at 2m
    ,ustar_nwp         & ! friction velocity m/s ustar^2 = tau/roa
    ,pzpbl             & ! stores H(ABL) for averaging and plotting purposes, m
    ,pwp               & ! Permanent Wilting Point
    ,fc                & ! Field Capacity
    ,invL_nwp          & ! inverse of the Monin-Obuhkov length
    ,model_surf_elevation! height above sea level of model surface. NB: can differ from physical value

!  temporary placement of solar radiation variations QUERY?
 
  real,target, public,allocatable, dimension(:,:), save:: &
       zen          &  ! Zenith angle (degrees)
      ,coszen       &  ! cos of zenith angle
      ,PARdbh       &  ! PAR, direct beam on horizontal surface, W/m2 !WN17
      ,PARdif       &  ! PAR, diffuse, W/m2 !WN17
      ,fCloud       &  ! cloud atten. factor (0-1), for Weiss&Norman approach !WN17
      ,Idiffuse     &  ! diffuse solar radiation (W/m^2)
      ,Idirect         ! total direct solar radiation (W/m^2)

  integer, parameter ::Nspecial2d = 0
  real,target, public,allocatable, dimension(:,:,:), save:: &
       special2d
  integer, parameter ::Nspecial3d = 0
  real,target, public,allocatable, dimension(:,:,:,:), save:: &
       special3d
  integer, public, save   :: ix_special2d(Nspecial2d),ix_special3d(Nspecial3d)

  real,target,public, save,allocatable, dimension(:,:) :: &   !st-dust
       clay_frac  &  ! clay fraction (%) in the soil
      ,sand_frac     ! sand fraction (%) in the soil
  real,target,public, save,allocatable, dimension(:,:) :: &
       surface_precip_old !precip from previous step for making differences (wrf)

  ! Different NWP outputs for soil water are possible. We can currently
  ! cope with two:
  character(len=10), public, save  :: SoilWaterSource="IFS"  ! IFS or PARLAM

  real,target,public, save, allocatable,dimension(:,:) :: &
    fSW     ! fSW= f(relative extractable water) =  (sw-swmin)/(swFC-swmin)

  real,target, public, dimension(:,:), save,allocatable  ::&
         xwf  ! extension of water fraction, save after 1st call

  integer, parameter, public :: NEXTEND = 2 ! no. box to side of (i,j) 

  integer, public, save   :: Nhh &         ! number of field stored per 24 hours
       ,nhour_first  ! time of the first meteo stored
! Logical flags, used to determine if some met fields are present in the
! input or not:
  logical, target, public, save :: &
     foundustar= .false.     & ! Used for MM5-type, where u_xmj* but not tau
    ,foundcc3d = .false.     & ! false if no cc3d in metdata
    ,foundSST= .false.       & ! false if no SeaSurfaceT in metdata
    ,foundSoilWater_uppr= .false.  & ! false if no SW-shallow
    ,foundSoilWater_deep= .false.  & ! false if no SW-deep
    ,foundrh2m= .false.   & ! false if no relative_humidity_2m in metdata
    ,foundtau= .false.   & ! false if no surface_stress in metdata
    ,foundsdepth= .false.    & ! false if no snow_flag depth in metdata
    ,foundHmix= .false.& ! false if no PBL height found
    ,foundice= .false.       & ! false if no ice_nwp coverage (%) in metdata
    ,foundKz_met= .false.    & ! false if no Kz from meteorology
    ,foundconv= .false.      & ! false if convection not found or not used
  ! Introduced for FUTURE NH3, but also sea-salt
    ,foundws10_met= .false.   & ! false if no u10 from meteorology
    ,foundu10_met= .false.   & ! false if no u10 from meteorology
    ,foundv10_met= .false.   & ! false if no v10 from meteorology
    ,foundprecip= .false.    & ! false if no precipitationfrom meteorology
    ,foundcloudwater= .false.& !false if no cloudwater found
    ,foundSMI1= .true.& ! false if no Soil Moisture Index level 1 (shallow)
    ,foundSMI3= .true.& ! false if no Soil Moisture Index level 3 (deep)
    ,foundrain= .false.& ! false if no rain found or used
    ,foundirainnc= .false. &! false if no irainnc found or used
    ,foundirainc= .false. &! false if no irainc found or used
    ,foundinvL= .false. &! false if topography file found. 
    ,foundtopo= .false. ! false if topography file found. 
     !NB: the value of model_surf_elevation will be set to default value based on surface pressure anyway

! specific indices of met
  integer, public, save   :: ix_u_xmj,ix_v_xmi, ix_q, ix_th, ix_cc3d, ix_pr, &
      ix_cw_met, ix_cnvuf, ix_cnvdf, ix_Kz_met, ix_roa, ix_SigmaKz, ix_EtaKz,&
      ix_Etadot, ix_cc3dmax, ix_lwc, ix_Kz_m2s, ix_u_mid, ix_v_mid, ix_ps, &
      ix_t2_nwp, ix_rh2m, ix_fh, ix_fl, ix_tau, ix_ustar_nwp, ix_sst, &
      ix_SoilWater_uppr, ix_SoilWater_deep, ix_sdepth, ix_ice_nwp, ix_ws_10m,&
      ix_surface_precip, ix_uw, ix_ue, ix_vs, ix_vn, ix_convective_precip, &
      ix_rain,ix_irainc,ix_irainnc, ix_elev, ix_invL, ix_pblnwp

  type,  public :: metfield
     character(len = 100) :: name = 'empty' !name as defined in external meteo file
     character(len = 100) :: unit = 'notset' !unit required by the model
     character(len = 100) :: validity = 'notset' !special conditions set by the external meteo file
     integer :: dim = 3 !number of dimension (2 for 2D, 3 for 3D)
     integer :: frequency =3  ! How many hours between two fields
     logical :: time_interpolate = .true. ! Interpolate in time  
     logical :: read_meteo = .false. ! The field will be looked for in the external meteo file
     logical :: needed= .true. ! The field must be present in the external meteo file
     logical, pointer :: found => null()  ! The field has been found in the external meteo file
!note that it is not allowed in fortran to define a target in a derived type
     real, pointer :: field(:,:,:,:) => null() !actual values for the fields; must be pointed to
     integer :: zsize = 1 ! field, size of third index
     integer :: msize = 1 ! field, size of fourth index
     real, pointer, dimension(:,:,:)::field_shared
     logical, pointer :: ready ! The field must be present in the external meteo file
     logical, pointer :: copied ! The field must be present in the external meteo file
  end type metfield
  logical, public,save, target::ready=.false.,copied=.false.

  integer, public, parameter   :: NmetfieldsMax=100 !maxnumber of metfields
  type(metfield),  public :: met(NmetfieldsMax)  !To put the metfields that need systematic treatment
  type(metfield),  public :: derivmet(20)  !DSA15 To put the metfields derived from NWP, eg for output
  logical, target :: metfieldfound(NmetfieldsMax)=.false. !default for met(ix)%found 
  integer, public, save   :: Nmetfields! number of fields defined in met
  integer, public, save   :: N3Dmetfields! number of 3D fields defined in met
  real,target, public,save,allocatable, dimension(:,:,:) :: uw,ue
  real,target, public,save,allocatable, dimension(:,:,:) :: vs,vn

  logical, public :: WRF_MET_CORRECTIONS = .false.
  logical, public :: MET_SHORT = .true.!metfields are stored as "short" (integer*2 and scaling)
  logical, public :: MET_C_GRID = .false.!true if u and v wind fields are in a C-staggered, larger grid.
  logical, public :: MET_REVERSE_K = .false.!set true if met fields are stored with lowest k at surface
  logical, public :: found_wrf_bucket = .false.
  real, public    :: wrf_bucket = 0.0 !constant used to define precipitation


  integer, public :: Nshared_2d
  integer, public :: Nshared_3d

  public :: Alloc_MetFields !allocate arrays

contains

subroutine Alloc_MetFields(LIMAX,LJMAX,KMAX_MID,KMAX_BND,NMET)
!allocate MetFields arrays arrays
  implicit none
  
  integer, intent(in) ::LIMAX,LJMAX,KMAX_MID,KMAX_BND,NMET
  integer ::ix,i,j,n,data_shape(3),xsize

  do ix=1,NmetfieldsMax
     met(ix)%found => metfieldfound(ix)!default target
     if(.not. associated(met(ix)%ready))met(ix)%ready=>ready
     if(.not. associated(met(ix)%copied))met(ix)%copied=>copied
  end do

  ix=1
  met(ix)%name             = 'u_wind'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(u_xmj(0:LIMAX,1:LJMAX,KMAX_MID,NMET))
  u_xmj=0.0
  met(ix)%field(0:LIMAX,1:LJMAX,1:KMAX_MID,1:NMET)  => u_xmj
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_u_xmj=ix

  ix=ix+1
  met(ix)%name             = 'v_wind'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(v_xmi(1:LIMAX,0:LJMAX,KMAX_MID,NMET))
  v_xmi=0.0
  met(ix)%field(1:LIMAX,0:LJMAX,1:KMAX_MID,1:NMET)  => v_xmi
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_v_xmi=ix

  ix=ix+1
  met(ix)%name             = 'specific_humidity' ! kg/kg
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(q(LIMAX,LJMAX,KMAX_MID,NMET))
  q=0.0
  met(ix)%field => q 
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_q=ix

  ix=ix+1
  met(ix)%name             = 'potential_temperature'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(th(LIMAX,LJMAX,KMAX_MID,NMET))
  th=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:NMET)  => th
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_th=ix

  ix=ix+1
  met(ix)%name             = '3D_cloudcover'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundcc3d 
  allocate(cc3d(LIMAX,LJMAX,KMAX_MID))
  cc3d=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:1)  => cc3d
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_cc3d=ix

  ix=ix+1
  met(ix)%name             = 'precipitation'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundprecip
  allocate(pr(LIMAX,LJMAX,KMAX_MID))
  pr=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:1)  => pr
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_pr=ix

  ix=ix+1
  met(ix)%name             = 'cloudwater'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            => foundcloudwater
  allocate(cw_met(LIMAX,LJMAX,KMAX_MID,NMET))
  cw_met=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:NMET)  => cw_met
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_cw_met=ix

  ix=ix+1
  met(ix)%name             = 'convective_updraft_flux'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = USES%CONVECTION
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(cnvuf(LIMAX,LJMAX,KMAX_BND))
  cnvuf=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_BND,1:1)  => cnvuf
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = 1
  ix_cnvuf=ix

  ix=ix+1
  met(ix)%name             = 'convective_downdraft_flux'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = USES%CONVECTION
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(cnvdf(LIMAX,LJMAX,KMAX_BND))
  cnvdf=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_BND,1:1)  => cnvdf
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = 1
  ix_cnvdf=ix

  ix=ix+1
  met(ix)%name             = 'eddy_diffusion_coefficient'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundKz_met
  allocate(Kz_met(LIMAX,LJMAX,KMAX_BND,NMET))
  Kz_met=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_BND,1:NMET)  => Kz_met
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET
  ix_Kz_met=ix

  ix=ix+1
  met(ix)%name             = 'air_density'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(roa(LIMAX,LJMAX,KMAX_MID,NMET))
  roa=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:NMET)  => roa
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_roa=ix

  ix=ix+1
  met(ix)%name             = 'Kz_sigmacoordinates'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(SigmaKz(LIMAX,LJMAX,KMAX_BND,NMET))
  SigmaKz=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_BND,1:NMET)  => SigmaKz
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET
  ix_SigmaKz=ix

  ix=ix+1
  met(ix)%name             = 'Kz_Etacoordinates'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(EtaKz(LIMAX,LJMAX,KMAX_BND,NMET))
  EtaKz=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_BND,1:NMET)  => EtaKz
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET
  ix_EtaKz=ix

  ix=ix+1
  met(ix)%name             = 'etadot'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(Etadot(LIMAX,LJMAX,KMAX_BND,NMET))
  Etadot=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_BND,1:NMET)  => Etadot
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET
  ix_Etadot=ix

  ix=ix+1
  met(ix)%name             = 'max_cloudcover'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(cc3dmax(LIMAX,LJMAX,KMAX_MID))
  cc3dmax=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:1)  => cc3dmax
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_cc3dmax=ix

  ix=ix+1
  met(ix)%name             = 'cloud_liquid_water'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(lwc(LIMAX,LJMAX,KMAX_MID))
  lwc=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:1)  => lwc
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_lwc=ix

  ix=ix+1
  met(ix)%name             = 'Kz'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(Kz_m2s(LIMAX,LJMAX,KMAX_MID))
  Kz_m2s=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:1)  => Kz_m2s
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_Kz_m2s=ix

  ix=ix+1
  met(ix)%name             = 'u_wind_3D'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(u_mid(LIMAX,LJMAX,KMAX_MID))
  u_mid=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:1)  => u_mid
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_u_mid=ix

  ix=ix+1
  met(ix)%name             = 'v_wind_3D'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(v_mid(LIMAX,LJMAX,KMAX_MID))
  v_mid=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:1)  => v_mid
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_v_mid=ix

  N3Dmetfields=ix
!  write(*,*)'number of 3D metfields: ',N3Dmetfields
!________________________________________________________________
! 2D fields

  ix=ix+1
  met(ix)%name             = 'surface_pressure'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(ps(LIMAX,LJMAX,NMET))
  ps=1.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => ps
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_ps=ix

!NWPHMIX
  ix=ix+1
  met(ix)%name             = 'PBLH' ! GLOBAL05
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  ! If we ask for this, we need it.
  if ( PBL%HmixMethod == 'NWP' )  then
    met(ix)%needed = .true.
    met(ix)%needed = .false. ! DS testing
    met(ix)%found            => foundHmix
    pbl_nwp=0.0
  else
    met(ix)%needed           = .false.
    met(ix)%found            = .false.
  end if
  allocate(pbl_nwp(LIMAX,LJMAX,NMET))
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => pbl_nwp
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_pblnwp=ix
!NWPHMIX


  ix=ix+1
  met(ix)%name             = 'temperature_2m'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(t2_nwp(LIMAX,LJMAX,NMET))
  t2_nwp=1.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => t2_nwp
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_t2_nwp=ix

  ix=ix+1
  met(ix)%name             = 'relative_humidity_2m'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundrh2m
  allocate(rh2m(LIMAX,LJMAX,NMET))
  rh2m=1.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => rh2m
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_rh2m=ix

  ix=ix+1
  met(ix)%name             = 'surface_flux_sensible_heat'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(fh(LIMAX,LJMAX,NMET))
  fh=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => fh
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_fh=ix

  ix=ix+1
  met(ix)%name             = 'surface_flux_latent_heat'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(fl(LIMAX,LJMAX,NMET))
  fl=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => fl
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_fl=ix

  ix=ix+1
  met(ix)%name             = 'surface_stress'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundtau
  allocate(tau(LIMAX,LJMAX,NMET))
  tau=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => tau
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_tau=ix

  ix=ix+1
  met(ix)%name             = 'ustar_nwp'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(ustar_nwp(LIMAX,LJMAX))
  ustar_nwp=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:1)  => ustar_nwp
  met(ix)%zsize = 1
  met(ix)%msize = 1
  ix_ustar_nwp=ix

  ix=ix+1
  met(ix)%name             = 'sea_surface_temperature'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundSST
  allocate(sst(LIMAX,LJMAX,NMET))
  sst=1.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => sst
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_sst=ix

  ix=ix+1
  met(ix)%name             = 'SMI1'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = USES%SOILWATER
  met(ix)%read_meteo       = USES%SOILWATER
  met(ix)%needed           = .false.
  met(ix)%found            => foundSoilWater_uppr
  allocate(SoilWater_uppr(LIMAX,LJMAX,NMET))
  SoilWater_uppr=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => SoilWater_uppr
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_SoilWater_uppr=ix

  ix=ix+1
  met(ix)%name             = 'SMI3'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = USES%SOILWATER
  met(ix)%read_meteo       = USES%SOILWATER
  met(ix)%needed           = .false.
  met(ix)%found            =>  foundSoilWater_deep
  allocate(SoilWater_deep(LIMAX,LJMAX,NMET))
  SoilWater_deep=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => SoilWater_deep
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_SoilWater_deep=ix

  ix=ix+1
  met(ix)%name             = 'snow_depth'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundsdepth
  allocate(sdepth(LIMAX,LJMAX,NMET))
  sdepth=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => sdepth
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_sdepth=ix

  ix=ix+1
  met(ix)%name             = 'fraction_of_ice'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundice
  allocate(ice_nwp(LIMAX,LJMAX,NMET))
  ice_nwp=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => ice_nwp
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_ice_nwp=ix

  ix=ix+1
  met(ix)%name             = 'u10'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundws10_met
  allocate(ws_10m(LIMAX,LJMAX,NMET))
  ws_10m=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:NMET)  => ws_10m
  met(ix)%zsize = 1
  met(ix)%msize = NMET
  ix_ws_10m=ix

  ix=ix+1
  met(ix)%name             = 'invL_nwp'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            => foundinvL
  allocate(invL_nwp(LIMAX,LJMAX))
  invL_nwp=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:1)  => invL_nwp
  met(ix)%zsize = 1
  met(ix)%msize = 1
  ix_invL=ix

  ix=ix+1
  met(ix)%name             = 'large_scale_precipitations'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  allocate(surface_precip(LIMAX,LJMAX))
  surface_precip=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:1)  => surface_precip
  met(ix)%zsize = 1
  met(ix)%msize = 1
  ix_surface_precip=ix

  ix=ix+1
  met(ix)%name             = 'convective_precipitations'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(convective_precip(LIMAX,LJMAX))
  convective_precip=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:1)  => convective_precip
  met(ix)%zsize = 1
  met(ix)%msize = 1
  ix_convective_precip=ix

  ix=ix+1
  met(ix)%name             = 'topography'!NB: as used in the met model; can differ from physical in mountainous areas.
  met(ix)%dim              = 2
  met(ix)%frequency        = 100000!constant
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.!read once only the first time
  met(ix)%needed           = .false.
  met(ix)%found            => foundtopo
  allocate(model_surf_elevation(LIMAX,LJMAX))
  model_surf_elevation=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:1)  => model_surf_elevation
  met(ix)%zsize = 1
  met(ix)%msize = 1
  ix_elev=ix

  ix=ix+1
  met(ix)%name             = 'neigbors_wind-uw'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(uw(LJMAX,KMAX_MID,NMET))
  uw=0.0
  met(ix)%field(1:1,1:LJMAX,1:KMAX_MID,1:NMET)  => uw
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_uw=ix

  ix=ix+1
  met(ix)%name             = 'neigbors_wind-ue'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(ue(LJMAX,KMAX_MID,NMET))
  ue=0.0
  met(ix)%field(1:1,1:LJMAX,1:KMAX_MID,1:NMET)  => ue
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_ue=ix

  ix=ix+1
  met(ix)%name             = 'neigbors_wind-vs'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(vs(LIMAX,KMAX_MID,NMET))
  vs=0.0
  met(ix)%field(1:LIMAX,1:1,1:KMAX_MID,1:NMET)  => vs
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_vs=ix

  ix=ix+1
  met(ix)%name             = 'neigbors_wind-vn'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .true.
  met(ix)%found            = .false.
  allocate(vn(LIMAX,KMAX_MID,NMET))
  vn=0.0
  met(ix)%field(1:LIMAX,1:1,1:KMAX_MID,1:NMET)  => vn
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET
  ix_vn=ix


!can be used to output any 2d field, using 'MET2D'
  do n = 1, Nspecial2d
  ix=ix+1
  write(met(ix)%name,fmt='(A,I0)')'special2d',n
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  if(n==1)then
     allocate(special2d(LIMAX,LJMAX,Nspecial2d))
     special2d=0.0
  endif
!  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:1)  => special2d(1:LIMAX,1:LJMAX,n)
!Since the syntax above is not allowed, we move the adress of the pointer, instead of the target
  met(ix)%field(1:LIMAX,1:LJMAX,1-(n-1):1-(n-1),1:1)  => special2d
  met(ix)%zsize = 1
  met(ix)%msize = 1
  ix_special2d(n)=ix
  enddo

!can be used to output any 2d field, using 'MET2D'
  do n = 1, Nspecial3d
  ix=ix+1
  write(met(ix)%name,fmt='(A,I0)')'special3d',n
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            = .false.
  if(n==1)then
     allocate(special3d(LIMAX,LJMAX,KMAX_MID,Nspecial3d))
     special3d=0.0
  endif
!  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:1)  => special3d(1:LIMAX,1:LJMAX,1:KMAX_MID,n)
!Since the syntax above is not allowed, we move the adress of the pointer, instead of the target
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1-(n-1):1-(n-1))  => special3d
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_special3d(n)=ix
  enddo



if(USE_WRF_MET_NAMES)then
   WRF_MET_CORRECTIONS = .true.
   MET_C_GRID = .true.
   MET_SHORT = .false. !metfields are stored as "float"
   MET_REVERSE_K = .true.!reverse k coordinates when reading
!names used in WRF metfiles
!3D
   met(ix_u_xmj)%name             = 'U' 
   met(ix_v_xmi)%name             = 'V'
   met(ix_q)%name                 = 'QVAPOR'
   met(ix_th)%name                = 'T'
   met(ix_cc3d)%name              = 'CLDFRA'
   met(ix_cw_met)%name            = 'QCLOUD'

!Use QRAIN to make 3D precip profiles from 2D. NB not accumulated, so cannot be used directly
  ix=ix+1
  met(ix)%name             = 'QRAIN'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.!to get new and old value stored
  met(ix)%read_meteo       = .true.
  met(ix)%needed           = .false.
  met(ix)%found            => foundrain
  allocate(rain(LIMAX,LJMAX,KMAX_MID,NMET))
  rain=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:KMAX_MID,1:NMET)  => rain
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1
  ix_rain=ix

!2D
!Use I_RAINNC and I_RAINC to make 2D precip 
  ix=ix+1
  met(ix)%name             = 'I_RAINNC'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.!to get new and old value stored
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            => foundirainnc
  allocate(irainnc(LIMAX,LJMAX))
  irainnc=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:1)  => irainnc
  met(ix)%zsize = 1
  met(ix)%msize = 1
  ix_irainnc=ix

  ix=ix+1
  met(ix)%name             = 'I_RAINC'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%needed           = .false.
  met(ix)%found            => foundirainc
  allocate(irainc(LIMAX,LJMAX))
  irainc=0.0
  met(ix)%field(1:LIMAX,1:LJMAX,1:1,1:1)  => irainc
  met(ix)%zsize = 1
  met(ix)%msize = 1
  ix_irainc=ix

!2D redefine names for wrf
   met(ix_surface_precip)%name    = 'RAINNC'
   met(ix_convective_precip)%name = 'RAINC'
   met(ix_ps)%name                = 'PSFC'
   met(ix_t2_nwp)%name            = 'T2'
   met(ix_fh)%name                = 'HFX'
   met(ix_fl)%name                = 'LH'
   met(ix_ustar_nwp)%name         = 'UST'
   met(ix_sst)%name               = 'SST'
   met(ix_ws_10m)%name            = 'U10'
   met(ix_SoilWater_uppr)%name    = 'SMI1'!take first level. Do not change name! (name set in Getmeteofield)
   met(ix_SoilWater_deep)%name    = 'SMI3'!take third level. Do not change name! (name set in Getmeteofield)
   met(ix_sdepth)%name            = 'SNOWH'!snowdepth in m
   met(ix_ice_nwp)%name           = 'SEAICE'!flag 0 or 1
   met(ix_rh2m)%name              = 'Q2' ! 2 meter relative humidity

!meteo model topography (assumed constant in time)
   met(ix_elev)%name             = 'HGT'
   TopoFile = date2string(meteo,startdate,mode='YMDH')

!... addmore
end if

  Nmetfields=ix
  if(Nmetfields>NmetfieldsMax)then
     write(*,*)"Increase NmetfieldsMax! "
     stop
  end if

    allocate(u_ref(LIMAX,LJMAX))
    allocate(rho_surf(LIMAX,LJMAX))
    allocate(Tpot2m(LIMAX,LJMAX))
    allocate(pzpbl(LIMAX,LJMAX))
    allocate(pwp(LIMAX,LJMAX))
    allocate(fc(LIMAX,LJMAX))
    allocate(xwf(LIMAX+2*NEXTEND,LJMAX+2*NEXTEND)) 
    allocate(fSW(LIMAX,LJMAX))
    fSW = 1.0
    allocate(zen(LIMAX, LJMAX))
    allocate(coszen(LIMAX, LJMAX))
    coszen=0.0
    allocate(Idiffuse(LIMAX, LJMAX))
    allocate(Idirect(LIMAX, LJMAX))
    allocate(PARdbh(LIMAX, LJMAX)) !WN17
    allocate(PARdif(LIMAX, LJMAX)) !WN17
    allocate(fCloud(LIMAX, LJMAX)) !WN17
    allocate(clay_frac(LIMAX, LJMAX))
    allocate(sand_frac(LIMAX, LJMAX))
    allocate(surface_precip_old(LIMAX,LJMAX))
    surface_precip_old=0.0

    if(NPROC/=NPROC_MPI)then
!allocate shared memory within large subdomains

    CALL MPI_BARRIER(MPI_COMM_SUB, IERROR)

    i=1
    j=1    
    do ix = 1, Nmetfields
       data_shape=(/1,1,1/)
       xsize=1
!       write(*,*)me_mpi,'share ',met(ix)%name
       call share_logical(met(ix)%ready,data_shape,xsize,MPI_COMM_SUB)
       call share_logical(met(ix)%copied,data_shape,xsize,MPI_COMM_SUB)
       if(met(ix)%dim==2)then
          data_shape=(/largeLIMAX,largeLJMAX,1/)
          xsize=largeLIMAX*largeLJMAX
          i=i+1
          call share(met(ix)%field_shared,data_shape,xsize,MPI_COMM_SUB)
       end if
       if(met(ix)%dim==3)then
          j=j+1
          data_shape=(/largeLIMAX,largeLJMAX,KMAX_MID/)
          xsize=largeLIMAX*largeLJMAX*KMAX_MID
          call share(met(ix)%field_shared,data_shape,xsize,MPI_COMM_SUB)
       end if
       CALL MPI_BARRIER(MPI_COMM_SUB, IERROR)
    end do
    Nshared_2d=i
    Nshared_3d=j
    end if
  end subroutine Alloc_MetFields

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


end module MetFields_mod
