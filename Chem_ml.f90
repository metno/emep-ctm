! <Chem_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2011 met.no
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module Chemfields_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
use Par_ml               , only: MAXLIMAX,MAXLJMAX   ! => x, y dimensions
use ModelConstants_ml    , only: KMAX_MID, KCHEMTOP     ! =>  z dimension
use ChemSpecs_adv_ml,  only: NSPEC_ADV         ! => No. species 
use ChemSpecs_shl_ml,  only: NSPEC_SHL         ! => No. species 
!see belowuse ChemSpecs_bgn_ml,  only: NSPEC_BGN         ! => No. species 
implicit none
private

!-------- this snipppet was from older GenSpec_bgn_ml. ------
  ! PRETTY MUCH FAKED FOR NOW. CAN BE DELETED SOON IN HOPE!
  !+ Defines indices and NSPEC for bgn : Background species

   ! Species which can be specified simply for each column, e.g.
   ! as function of local meteorology or zenith angle
   !   o2, m,  and for MADE-like, oh, ch3coo2

   integer, public, parameter ::  NSPEC_BGN = 0 ! No. 3D bgn species
   integer, public, parameter ::  NSPEC_COL = 0 ! total no. prescribed specs

    !/ define xn_2d_bgn here.
     real, public, save, dimension(1,KCHEMTOP:KMAX_MID) :: xn_2d_bgn

!-------- end of this snipppet from older GenSpec_bgn_ml. ------

    !----------------- basic chemical fields ----------------------------------!
    !  Here we declare and initialise to zero the chemical fields used in the  !
    !  model, as well as cfac (converts from 50m to 1m/3m output)         ! 
    !---------------------------------------------------------------------!

  real, save, allocatable, public :: &
     xn_adv(:,:,:,:)  &
    ,xn_shl(:,:,:,:)  &
    ,xn_bgn(:,:,:,:) &
    ,PM25_water(:,:,:) &  !3D PM water
    ,PM25_water_rh50(:,:)    !gravimetric PM water

  real, public, save, allocatable:: Fgas3d (:,:,:,:)  ! for SOA

  real, public, save, allocatable, dimension(:,:) :: AOD

  real, save, allocatable, public :: &
     cfac   (:,:,:)   

  real, save, allocatable, public :: &
     so2nh3_24hr(:,:)!hf CoDep

  real, save, allocatable, public :: &
     Grid_snow(:,:) !snow_flag fraction in grid

  public ::alloc_ChemFields

contains

  subroutine alloc_ChemFields

    implicit none

    allocate(xn_adv(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID))
    xn_adv=0.0
    allocate(xn_shl(NSPEC_SHL,MAXLIMAX,MAXLJMAX,KMAX_MID))
    xn_shl=0.0
    allocate(xn_bgn(NSPEC_BGN,MAXLIMAX,MAXLJMAX,KMAX_MID))
    xn_bgn=0.0
    allocate(PM25_water(MAXLIMAX,MAXLJMAX,KMAX_MID))
    PM25_water=0.0
    allocate(PM25_water_rh50(MAXLIMAX,MAXLJMAX))
    PM25_water_rh50=0.0
    allocate(AOD(MAXLIMAX,MAXLJMAX))
    allocate(cfac(NSPEC_ADV,MAXLIMAX,MAXLJMAX))
    cfac=1.0
    allocate(so2nh3_24hr(MAXLIMAX,MAXLJMAX))
    so2nh3_24hr=0.0
    allocate(Grid_snow(MAXLIMAX,MAXLJMAX))
    Grid_snow=0.0


  end subroutine alloc_ChemFields


!_____________________________________________________________________________
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
                     end module Chemfields_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!_____________________________________________________________________________
