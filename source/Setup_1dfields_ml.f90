! <Setup_1dfields_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________!
 module Setup_1dfields_ml

  ! Arrays of meteorology and concentration for 1-D column , for input to 
  ! chemical solver ........
  ! The k-dimension spans the ground (KMAX_MID) to the K-values
  ! specified by KCHEMTOP - here 2
  !
  ! - new aray added to keep o2, m, and for MADE oh, etc

  use ModelConstants_ml,     only :  KMAX_MID, KCHEMTOP, KUPPER, NBVOC
  use EmisDef_ml,            only :  NSS, NDU !SeaS, Dust
  use ChemSpecs_tot_ml,      only :  NSPEC_TOT, FIRST_SEMIVOL, LAST_SEMIVOL
  use ChemSpecs_shl_ml,      only :  NSPEC_SHL
  use Chemfields_ml,         only :  NSPEC_COL
  implicit none
  private


  !/ variables to keep track of which call

   logical, public, save :: first_call = .true.
   integer, public, save :: ncalls     =      0

  !/-- the chemistry is calculated for arrays of size: 

   integer, public, parameter  :: CHEMSIZE = KMAX_MID-KCHEMTOP+1 !  
 
   real, public, dimension(NSPEC_TOT,KCHEMTOP:KMAX_MID), save :: &
                   xn_2d            ! Concentrations [molecules/cm3]  

! For semivolatiles we track the farction as gas and particle- used for SOA
! We use NSPEC_TOT to allow us to write Fpart for FFUEL and WOOD also -
! these may be semivol one day.
   !real, public, dimension(FIRST_SOA:LAST_SOA,KCHEMTOP:KMAX_MID), save :: &
   real, public, dimension(NSPEC_TOT,KCHEMTOP:KMAX_MID), save :: &
                   Fgas  = 1.0     &! Fraction as gas-phase
                  ,Fpart = 0.0      ! Fraction as gas-phase

 !Emissions in column. We assume that these only involve advected species
   real, public, dimension(NSPEC_SHL+1:NSPEC_TOT,KCHEMTOP:KMAX_MID), save ::&
         rcemis   !emissions

  ! We define a column array for isoprene and terpene for use in
  ! the chemical solver. All values except for k=KMAX_MID will
  ! remain zero however

  ! Emission arrays:
   real, public, dimension(NBVOC,KCHEMTOP:KMAX_MID), save :: rcbio = 0.0 ! BVOC
  !FUTURE real, public, dimension(KCHEMTOP:KMAX_MID), save   :: rcnh3 
   real, public, dimension(KCHEMTOP:KMAX_MID), save   :: rc_Rn222 = 0.0  ! 210Pb
   real, public, dimension(NSS,KCHEMTOP:KMAX_MID), save :: rcss = 0.0  ! Sea salt
   real, public, dimension(NDU,KCHEMTOP:KMAX_MID), save :: rcwbd = 0.0 ! windblown dust

   real, public, dimension(KCHEMTOP:KMAX_MID), save :: &
          rh                  & ! RH (fraction, 0-1)
         ,amk                 & ! M - atmospheric conc.
         ,o2, n2              & ! oxygen, nitrogen
         ,h2o                 & ! water
         ,temp                & ! temperature
         ,tinv                & ! inverse temp
         ,pp                     !pressure

   integer, public, dimension(KCHEMTOP:KMAX_MID), save :: &
          itemp                  ! int of temperature

 end module Setup_1dfields_ml
