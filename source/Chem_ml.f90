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

  real, save, public :: &
     xn_adv (NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0  &
    ,xn_shl (NSPEC_SHL,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0  &
    ,xn_bgn (NSPEC_BGN,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0  &
    ,PM25_water (MAXLIMAX,MAXLJMAX,KMAX_MID)         = 0.0  &  !3D PM water
    ,PM25_water_rh50 (MAXLIMAX,MAXLJMAX)             = 0.0     !gravimetric PM water

  real, public, dimension(MAXLIMAX,MAXLJMAX) :: AOD

  real, save, public :: &
     cfac   (NSPEC_ADV,MAXLIMAX,MAXLJMAX) = 1.0    

  real, save, public :: &
     so2nh3_24hr(MAXLIMAX,MAXLJMAX) = 0.0 !hf CoDep

  real, save, public :: &
     Grid_snow(MAXLIMAX,MAXLJMAX) = 0.0 !snow_flag fraction in grid


!_____________________________________________________________________________
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
                     end module Chemfields_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!_____________________________________________________________________________
