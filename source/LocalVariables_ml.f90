! <LocalVariables_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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
!ALBEDO?
!SNOW?
!ALSO CHECK/CHANGE lu==1 test in fphen, DO3SE_ml
!NNLAND vs NLAND???
!! PARsun, shade is set by Rsurace

module LocalVariables_ml
 use ModelConstants_ml, only : NLANDUSE
 implicit none
 private

 ! public :: set_local   ! Copies Sub(lu)%xxx to xxx

! Near-surface meteorology and other variables for local area,
! e.g. for a measurement site or for a specific landuse within a grid square

    real, private, parameter :: NOT_SET = -999.999  ! Fake value to check for
    integer, private, parameter :: INOT_SET = -999  ! Fake value to check for
                                                    ! variables being set

    integer, public, save :: iL = INOT_SET   ! Landuse index

 ! 1) Grid data which should be okay for local  ======================

  type, public :: GridDat     
    real    :: latitude    ! deg N
    real    :: longitude   ! deg E
    integer :: i           ! index
    integer :: j           ! index
    logical :: is_wet      !  true if precip > 0
    logical :: is_NWPsea      ! NWP model defines this square as sea
    real    :: precip      ! Precip at surface
    real    :: wetarea     ! Fraction of grid which is wet
    real    :: cloud       ! Cloud-cover (fraction)
    integer :: snow        ! 1=snow present, 0 = no snow
    real    :: psurf       ! Surface pressure (Pa)
    real    :: z_ref       !  Height of grid centre (m)
    real    :: DeltaZ      !  Depth of grid centre (m)
    real    :: qw_ref      !  Specific humidity
    real    :: rho_ref     !  Air density (kg/m3)
  ! the following are likely used in Sub below also
    real :: t2C                ! Surface (2m) temperature in degrees C
    real :: t2                 ! Surface (2m) temperature in degrees K
    real :: rh                 ! Relative humidity, fraction (0-1)
    real :: rho_s              !  Air density (kg/m3) at surface, here 2m
    real :: vpd                ! Vapour pressure deficit  (kPa) ! CHECK UNITS
    real :: SWP                ! SWP  ! CHECK UNITS
    real :: ustar              ! friction velocity, m/s
    real :: wstar              ! convective velocity scale, m/s
    real :: invL               ! 1/L, where L is Obukhiov length (1/m)
    real :: Hd                !  Sensible Heat flux, *away from* surface
    real :: LE                !  Latent Heat flux, *away from* surface
    real :: theta_ref         ! Pot. temp at grid center
    real :: Ra_ref            !
    real :: u_ref             ! wind speed at ref. height
    real :: Ra_2m             !
    real :: Ra_3m             !
    real :: so2nh3ratio       !  for CEH deposition scheme
    real :: &    !Not quite sure how many of these we need. 
      solar     = NOT_SET   &  ! => irradiance (W/m^2)
     ,Idirectn  = NOT_SET   &  ! => irradiance (W/m^2), normal to beam
     ,Idiffuse  = NOT_SET   &  ! => diffuse solar radiation (W/m^2)
     ,Idirect   = NOT_SET   &  ! => total direct solar radiation (W/m^2)
     ,zen       = NOT_SET   &  !   Zenith angle (degrees)
     ,coszen    = NOT_SET      ! = cos(zen)
                                ! (= sinB, where B is elevation angle)
     integer :: izen = INOT_SET ! int(zen)
!
!  real, dimension(NDRYDEP_TOT) :: &
!       Vg_ref   &! Grid average of Vg at ref ht. (effective Vg for cell)
!      ,Vg_3m     ! and at 3m
  end type GridDat

  type(GridDat), public, save :: Grid

 ! 2) Near-surface Data  -  Sub-grid ====================
 ! +  Sub-grid Veg/landcover data  ====================================
       

  type, public :: SubDat     
   !*
    integer :: &
        iL = INOT_SET     & ! Landcover index
       ,SGS  =  INOT_SET  & ! Start, growing seasons (day num)
       ,EGS  =  INOT_SET    ! End, growing seasons (day num)
   !*
    logical :: &
        is_forest, is_water 
   !*
    real :: &
     t2C       = NOT_SET &! Surface (2m) temperature in degrees C
    ,t2        = NOT_SET &! Surface (2m) temperature in degrees K
    ,rh        = NOT_SET &! Relative humidity, fraction (0-1)
    ,rho_s     = NOT_SET &! Air density (kg/m3) at surface, here 2m
    ,vpd       = NOT_SET &! Vapour pressure deficit  (kPa) ! CHECK UNITS
    ,SWP       = NOT_SET &! SWP  ! CHECK UNITS
    ,ustar     = NOT_SET &! friction velocity, m/s
    ,wstar     = NOT_SET &! convective velocity scale, m/s
    ,invL      = NOT_SET &! 1/L, where L is Obukhiov length (1/m)
    ,Hd        = NOT_SET &!  Sensible Heat flux, *away from* surface
    ,LE        = NOT_SET &!  Latent Heat flux, *away from* surface
    ,Ra_ref    = NOT_SET &!
    ,Ra_2m     = NOT_SET &!
    ,Ra_3m     = NOT_SET &!
    ,RgsO      = NOT_SET &! ground-surface resistances - set in DO3SE
    ,RgsS      = NOT_SET &! ground-surface resistances - set in DO3SE
   !
    ,coverage  = NOT_SET &! Area covered (fraction)
    ,LAI       = NOT_SET &! Leaf area index (m2/m2)
    ,SAI       = NOT_SET &! Surface area index (m2/m2)
    ,hveg      = NOT_SET &! Height of veg.      (m)
    ,d         = NOT_SET &! displacement height (m)
    ,z_refd    = NOT_SET &! z_ref - d (m)
    ,z0        = NOT_SET &! roughness length    (m)
  !
  ! Canopy-Associated Radiation
    ,PARsun    = NOT_SET &! photosynthetic active radn. for sun-leaves
    ,PARshade  = NOT_SET &!  " " for shade leaves
    ,LAIsunfrac= NOT_SET &! fraction of LAI in sun
  ! outputs from Rsurface will include:
    ,g_sto     = NOT_SET &! stomatal conductance (m/s)
    ,g_sun     = NOT_SET  ! g_sto for sunlit upper-canopy (flag) leaves

    !,ObsRad    = NOT_SET &! Used for box-model, for observed values
    !,snow      = NOT_SET &!  Usually from Grid
    !,wetarea   = NOT_SET &!  Usually from Grid
    !,psurf     = NOT_SET &!  Surface Pressure (Pa), Usually from Grid
    !,soil      = NOT_SET  ! Not used yet.
    !,b_inc     = NOT_SET &! in-canopy factor (Erisman-type)
 !
  end type SubDat

  type(SubDat), public, dimension(NLANDUSE), save :: Sub
  type(SubDat), public, save :: L    ! For just one land-class



!contains
!  subroutine set_local(lu)
!    integer, intent(in) :: lu
!
!       iL  = lu   ! Is this sensible/needed?!
!       t2C = Sub(lu)%t2C
!       LAI = Sub(lu)%LAI
!       SAI = Sub(lu)%SAI !! SAIadded????
!       rh  = Sub(lu)%rh  !! SAIadded????
!       ustar  = Sub(lu)%ustar  !! SAIadded????
!       hveg   = Sub(lu)%hveg  !! SAIadded????
!       PARsun     = Sub(lu)%PARsun 
!       PARshade   = Sub(lu)%PARshade 
!       LAIsunfrac = Sub(lu)%LAIsunfrac 
!
!       snow         = Grid%snow
!       wetarea      = Grid%wetarea
!       so2nh3ratio  = Grid%so2nh3ratio  !! SAIadded????
!
!  end subroutine set_local
!
end module LocalVariables_ml
