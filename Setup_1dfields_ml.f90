! <Setup_1dfields_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-201409 met.no
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
!_____________________________________________________________________________!
 module Setup_1dfields_ml

  ! Arrays of meteorology and concentration for 1-D column , for input to
  ! chemical solver ........
  ! The k-dimension spans the ground (KMAX_MID) to the K-values
  ! specified by KCHEMTOP - here 2
  !
  ! - new aray added to keep o2, m, and for MADE oh, etc

  implicit none
  private


  !/ variables to keep track of which call

   logical, public, save :: first_call = .true.
   integer, public, save :: ncalls     =      0

  !/-- the chemistry is calculated for arrays of size:

   integer, public, save  :: CHEMSIZE  !

 !FIELDS ALLOCATED IN Chem_ml.f90

  real, public, allocatable, dimension(:,:), save :: &
                   xn_2d            ! Concentrations [molecules/cm3]

! For semivolatiles we track the farction as gas and particle- used for SOA
! We use NSPEC_TOT to allow us to write Fpart for FFUEL and WOOD also -
! these may be semivol one day.
   !real, public, dimension(FIRST_SOA:LAST_SOA,KCHEMTOP:KMAX_MID), save :: &
   real, public, allocatable, dimension(:,:), save :: &
                   Fgas       &! Fraction as gas-phase
                  ,Fpart      ! Fraction as gas-phase

 !Emissions in column. We assume that these only involve advected species
   real, public, allocatable, dimension(:,:), save ::&
         rcemis     ! emissions
  ! We define a column array for isoprene and terpene for use in
  ! the chemical solver. All values except for k=KMAX_MID will
  ! remain zero however

   real, public, allocatable, dimension(:), save :: &
          rh                  & ! RH (fraction, 0-1)
         ,amk                 & ! M - atmospheric conc.
         ,o2, n2              & ! oxygen, nitrogen
         ,h2o                 & ! water
         ,temp                & ! temperature
         ,tinv                & ! inverse temp
         ,pp                     !pressure

   integer, public, allocatable, dimension(:), save :: &
          itemp                  ! int of temperature

 end module Setup_1dfields_ml
