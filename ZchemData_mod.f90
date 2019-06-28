! <ZchemData_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
! <ZchemData_mod.f90 - part of the EMEP MSC-W Chemical transport Model>
!_____________________________________________________________________________!
 module ZchemData_mod

  ! Arrays of meteorology and concentration for 1-D column , for input to
  ! chemical solver ........
  ! The k-dimension spans the ground (KMAX_MID) to the K-values
  ! specified by KCHEMTOP - here 2
  !
  ! - new aray added to keep o2, m, and for MADE oh, etc

  use AeroConstants_mod, only: AERO ! for %NSAREA = No types surface area
  use AllocInits,        only: AllocInit
  use ChemDims_mod,      only: NCHEMRATES, nspec => NSPEC_TOT, &
      NPHOTOLRATES
  use Config_module,     only: KCHEMTOP, KMAX_MID
  implicit none
  private

  !/ variables to keep track of which call

   logical, public, save :: first_call = .true.
   integer, public, save :: ncalls     =      0

  !/-- the chemistry is calculated for arrays of size:

   integer, public, save  :: CHEMSIZE  !

 !Emissions in column. We assume that these only involve advected species
  real, public, save, allocatable, dimension(:,:) :: & ! dims: nspec x nz
      rct               & !< chemical rate coeffs.  molec/cm3/s
     ,rcemis            & !< Emissions rate coeff.  molec/cm3/s
     ,rcbio             & !< Emissions rate coeff.  molec/cm3/s (BVOC, soil-NO, etc.)
     ,rcphot              !< Photolysis rates

  real, public, allocatable, dimension(:,:), save :: &
       xn_2d            ! Concentrations [molecules/cm3]

! For semivolatiles we track the farction as gas and particle- used for SOA
! We use NSPEC_TOT to allow us to write Fpart for FFUEL and WOOD also -
! these may be semivol one day.

   real, public, allocatable, dimension(:,:), save :: &
       Fgas       &! Fraction as gas-phase
      ,Fpart      ! Fraction as gas-phase

  ! We define a column array for isoprene and terpene for use in
  ! the chemical solver. All values except for k=KMAX_MID will
  ! remain zero however

   real, public, allocatable, dimension(:), save :: &
       rh                  & ! RH (fraction, 0-1)
      ,M                   & ! M - atmospheric conc. (was amk)
      ,o2, n2              & ! oxygen, nitrogen
      ,h2o                 & ! water
      ,temp                & ! temperature
      ,tinv                & ! inverse temp
      ,cN2O5               & ! mol speed, N2O5
      ,cHNO3               & ! mol speed, HNO3 
      ,cHO2                & ! mol speed, HO2  
      ,cNO2                & ! mol speed, NO2    ! kHet tests
      ,cNO3                & ! mol speed, NO2    ! kHet tests
      ,cO3                 & ! mol speed, O3   
      ,gamN2O5               ! jAero gamma values for output

   real, public, allocatable, dimension(:,:), save :: &
       DpgNw  & ! wet diameter,           dim:NSAREA,k
      ,S_m2m3   ! surface area, m2/m3     dim:NSAREA,k

   real, public, allocatable, dimension(:), save :: &
       deltaZcm             & ! layer thickness, cm
      ,aero_fom, aero_fss, aero_fdust, aero_fbc &! fractions
      ,pp                     !pressure
!         ,ugdryPM             & ! for wet radius from Gerber, etc.

   integer, public, allocatable, dimension(:), save :: &
       itemp                  ! int of temperature

 
 end module ZchemData_mod
