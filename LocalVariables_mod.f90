! <LocalVariables_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
! <LocalVariables_mod.f90 - A component of the EMEP MSC-W Euleria Chemical transport Model>
!*****************************************************************************!
module LocalVariables_mod
! -----------------------------------------------------------------------
! Near-surface meteorology and other variables for local area,
! e.g. for a measurement site or for a specific landuse within a grid square
! -----------------------------------------------------------------------

! To simplify some other dependencies we removed this:
! use GasParticleCoeffs_mod,         only: NDRYDEP_CALC
! And use NLOCDRYDEP_MAX below (will test in DryDep )

implicit none
private

integer, public, parameter :: NLOCDRYDEP_MAX = 99   ! 
real,    private, parameter ::  NOT_SET = -999.     ! Fake value to check for
integer, private, parameter :: INOT_SET = -999      ! variables being set

integer, public, save :: iL = INOT_SET              ! Landuse index

! -----------------------------------------------------------------------
! 1) Grid data which should be okay for local
! -----------------------------------------------------------------------
type, public :: GridDat
  real    :: latitude     ! deg N
  real    :: longitude    ! deg E
  integer :: i            ! index
  integer :: j            ! index
  logical :: is_wet       ! true if precip > 0
  logical :: is_frozen    ! true if T< -2  ! BiDir CRUDE
  logical :: is_mainlysea !  Usually > 50% sea/water
  logical :: is_allsea    ! Only sea in grid-square
  real    :: precip       ! Precip at surface
  real    :: wetarea      ! Fraction of grid which is wet
  real    :: cloud        ! Cloud-cover (fraction)
  logical ::  snowice     ! true is sdepth > 0 or ice>0
  real    :: sdepth       ! snowdepth (m)
  real    :: ice_nwp      ! ice_nwp (%)
  real    :: psurf        ! Surface pressure (Pa)
  real    :: z_ref        ! Used top of SL, = min(0.1 zi, z_mid)
  real    :: z_mid        ! Height of grid centre (m)
  real    :: DeltaZ       ! Depth of grid centre (m)
  real    :: qw_ref       ! Specific humidity
  real    :: rho_ref      ! Air density (kg/m3)
! the following are likely used in Sub below also
  real    :: t2C          ! Surface (2m) temperature in degrees C
  real    :: t2           ! Surface (2m) temperature in degrees K
  real    :: sst          ! Sea surface temperature in degrees K
  real    :: rh2m         ! Relative humidity, fraction (0-1)
  real    :: rho_s        !  Air density (kg/m3) at surface, here 2m
  real    :: vpd          ! Vapour pressure deficit  (kPa) ! CHECK UNITS
  real    :: SWP          ! SWP  ! CHECK UNITS
  real    :: fSW          ! fSW  - function for soil-water, 0--1
  real    :: ustar        ! friction velocity, m/s
  real    :: wstar        ! convective velocity scale, m/s
  real    :: invL         ! 1/L, where L is Obukhiov length (1/m)
  real    :: Hd           !  Sensible Heat flux, *away from* surface
  real    :: LE           !  Latent Heat flux, *away from* surface
  real    :: theta_ref    ! Pot. temp at grid center
  real    :: Ra_ref       !
  real    :: u_ref        ! wind speed at ref. height
  real    :: Ra_2m        !
  real    :: Ra_3m        !
  real    :: so2nh3ratio  !  for CEH deposition scheme
  real    :: surf_o3_ppb  !  for  EU AOTs, 3m O3
  real    :: surf_o3_ppb1 !   .... after deploss
! CoDep
  real    :: so2nh3ratio24hr  !  for CEH SO2 deposition scheme
  real    ::              & ! some on this set might not be need
     solar     = NOT_SET  & ! => irradiance (W/m^2)
    ,Idirectn  = NOT_SET  &   ! => irradiance (W/m^2), normal to beam
    ,Idiffuse  = NOT_SET  & ! => diffuse solar radiation (W/m^2)
    ,Idirect   = NOT_SET  & ! => total direct solar radiation (W/m^2)
    ,zen       = NOT_SET  & !   Zenith angle (degrees)
    ,coszen    = NOT_SET    ! = cos(zen)(= sinB, where B is elevation angle)
  integer :: izen = INOT_SET  ! int(zen)
  real, dimension(NLOCDRYDEP_MAX) :: &  ! for species subject to dry depostion
     Vg_ref = 0.0   & ! Grid average of Vg at ref ht. (effective Vg for cell)
    ,StoFrac = 0.0  & ! Fraction of flux (Vg) going through stomata.
    ,Vg_3m    & ! and at 3m
    ,Gsur,Gsto, Gns
end type GridDat

type(GridDat), public, save :: Grid

! -----------------------------------------------------------------------
! 2) Near-surface & Sub-grid Veg/landcover data
! -----------------------------------------------------------------------
type, public :: SubDat
  integer ::              &
   iL  = INOT_SET         & ! Landcover index
  ,SGS = INOT_SET         & ! Start, growing seasons (day num)
  ,EGS = INOT_SET           ! End, growing seasons (day num)
  logical :: &
    is_forest, is_water , is_veg, is_ice, is_crop
  real ::                 &
     t2C       = NOT_SET  & ! Surface (2m) temperature in degrees C
    ,t2        = NOT_SET  & ! Surface (2m) temperature in degrees K
    ,sst       = NOT_SET  & ! Sea surface temperature in degrees K
    ,rh        = NOT_SET  & ! Relative humidity, fraction (0-1)
    ,rho_s     = NOT_SET  & ! Air density (kg/m3) at surface, here 2m
    ,vpd       = NOT_SET  & ! Vapour pressure deficit  (kPa) ! CHECK UNITS
    ,EvapTransp= NOT_SET  & ! Evapotranspiration             ! CHECK UNITS
    ,SWP       = NOT_SET  & ! SWP  ! CHECK UNITS
    ,fSW       = NOT_SET  & ! function for fSWP or fSMD or...
    ,ustar     = NOT_SET  & ! friction velocity, m/s
    ,wstar     = NOT_SET  & ! convective velocity scale, m/s
    ,invL      = NOT_SET  & ! 1/L, where L is Obukhiov length (1/m)
    ,Hd        = NOT_SET  & !  Sensible Heat flux, *away from* surface
    ,LE        = NOT_SET    !  Latent Heat flux, *away from* surface
  real ::                 &
     Ra_ref    = NOT_SET  & !pw: name not appropriate for current use! is used as "Ra_mid"
    ,Ra_X      = NOT_SET  & !pw: temporary name
    ,Ra_2m     = NOT_SET  &
    ,Ra_3m     = NOT_SET  &
    ,RgsO      = NOT_SET  & ! ground-surface resistances - set in DO3SE
    ,RgsS      = NOT_SET  & ! ground-surface resistances - set in DO3SE
    ,coverage  = NOT_SET  & ! Area covered (fraction)
    ,LAI       = NOT_SET  & ! Leaf area index (m2/m2)
    ,SAI       = NOT_SET  & ! Surface area index (m2/m2)
    ,hveg      = NOT_SET  & ! Height of veg.      (m)
    ,d         = NOT_SET  & ! displacement height (m)
    ,z_refd    = NOT_SET  & ! z_ref - d (m)
    ,z0        = NOT_SET  & ! roughness length    (m)
! Canopy-Associated Radiation
    ,PARsun    = NOT_SET  & ! photosynthetic active radn. for sun-leaves, W/m2
    ,PARshade  = NOT_SET  & !  " " for shade leaves, W/m2
    ,LAIsunfrac= NOT_SET    ! fraction of LAI in sun
! outputs from Rsurface will include:
  real ::                 &
     g_sto     = NOT_SET  & ! stomatal conductance (m/s)
    ,g_sun     = NOT_SET  & ! g_sto for sunlit upper-canopy (flag) leaves
    ,f_sun     = NOT_SET  & ! f_env for SPOD?
    ,f_shade   = NOT_SET  & ! f_env for SPOD?
    ,f_env     = NOT_SET  & ! f_env for SPOD?
    ,f_vpd     = NOT_SET  & ! f_env for SPOD?
    ,f_light     = NOT_SET  & ! f_env for SPOD?
    ,f_temp      = NOT_SET  & ! f_env for SPOD?
    ,f_phen      = NOT_SET  & ! f_env for SPOD?
    ,f_min       = NOT_SET  & ! f_env for SPOD?
! and enable concentrations at canopy height:
    ,cano3_ppb  = 0.0     & ! Use 0.0 to make d_2d behave better
    ,cano3_nmole= 0.0     & ! Use 0.0 to make d_2d behave better
    ,FstO3      = 0.0       ! leaf O3 flux, nmole/m2/s
  real, dimension(NLOCDRYDEP_MAX) :: & ! for species subject to dry depostion
     Vg_ref   &  ! Grid average of Vg at ref ht. (effective Vg for cell)
    ,Vg_eff   &  ! Grid average of Vg effective at ref ht. (effective Vg for cell)
    ,Vg_3m    &  ! and at 3m
    ,StoFrac = 0.0  & ! Fraction of flux (Vg) going through stomata.
    ,Gsur, Gsto, Gns
end type SubDat

! MOVED TO Mosaic type(SubDat), public, dimension(0:NLANDUSEMAX), save :: Sub
type(SubDat), public, save :: L         ! For just one land-class
type(SubDat), public, save :: ResetSub  ! Keeps NOT_SET values
end module LocalVariables_mod
