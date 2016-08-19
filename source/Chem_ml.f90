! <Chem_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module Chemfields_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
use Par_ml               , only: MAXLIMAX,MAXLJMAX   ! => x, y dimensions
use ModelConstants_ml    , only: KMAX_MID     ! =>  z dimension
use GenSpec_adv_ml,  only: NSPEC_ADV         ! => No. species 
use GenSpec_shl_ml,  only: NSPEC_SHL         ! => No. species 
use GenSpec_bgn_ml,  only: NSPEC_BGN         ! => No. species 
implicit none
private

    !----------------- basic chemical fields ----------------------------------!
    !  Here we declare and initialise to zero the chemical fields used in the  !
    !  model, as well as cfac (converts from 50m to 1m/3m output)         ! 
    !---------------------------------------------------------------------!

  real, save, public :: &
     xn_adv (NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0     &
    ,xn_shl (NSPEC_SHL,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0     &
    ,xn_bgn (NSPEC_BGN,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0     &
    ,PM_water (MAXLIMAX,MAXLJMAX,KMAX_MID)           = 0.0        !water

  real, save, public :: &
     cfac   (NSPEC_ADV,MAXLIMAX,MAXLJMAX) = 1.0    

!_____________________________________________________________________________
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
                     end module Chemfields_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!_____________________________________________________________________________
