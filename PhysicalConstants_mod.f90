! <PhysicalConstants_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module PhysicalConstants_mod
!----------------------------------------------------------------------------
!  Defines Physical constants 
!----------------------------------------------------------------------------
implicit none
!F private
!-- contains no subroutine:
real , public, parameter ::         &
  AVOG   = 6.023e23                 & ! Avogadros number
, ATWAIR = 28.964                   & ! mol wt of air, g/mol (Jones, 1992)
, RGAS_ATML = 0.08205               & ! Molar Gas constant (atm M-1 K-1)
, RGAS_KG   = 287.0                 & ! Molar Gas constant (J K-1 kg-1)
, RGAS_J    = 8.3144                  ! Molar Gas constant (J mol-1 K-1)

                                        ! NB. ( J = N m2 = kg m2 s-2 )
                                        !       M = mol l-1
real, public, parameter  ::    &
     GRAV    = 9.807           &   ! Gravity, m s-2
  ,  EARTH_RADIUS = 6.37e6     &   ! 
  ,  CP      = 1004.0          &   ! Specific heat at const. pressure
  ,  KAPPA   = RGAS_KG/CP      &   
  ,  KARMAN  = 0.41            &   ! Von Karman  (=0.35 elsehwere in code!)
  ,  PI      = 3.141592653589793238462643383279 & ! www.verbose.net/Pi.html
  ,  DEG2RAD = PI/180.0        &   ! COnverts degrees to radians
  ,  RAD2DEG = 180.0/PI        &   ! COnverts radians to degrees
  ,  ROWATER = 1000.0          &   ! pw density of water kg m-3
  ,  BOLTZMANN = 1.380e-23     &   ! Boltzmann'c constant[J/deg/molec]
  ,  FREEPATH  = 6.5e-8        &   ! Mean Free Path of air [m]
  ,  VISCO     = 1.46e-5           ! Air viscosity [m2/s]   (was NU)


! Converts from mol/cm3 to nmole/m3
real, public, parameter :: NMOLE_M3 = 1.0e6*1.0e9/AVOG  

! Some definitions for daylight, in terms of zenith angle and cos(zen):
! (calculated from criteria that DAY_COSZEN > 1.0e-10 as daytime)

real, public, parameter  ::  &
     DAY_ZEN   = 89.9999999942704 & !
    ,DAY_COSZEN = 1.0e-10

!=================== DEP CODE ================================Y:0

! CHARNOCK is used to calculate the roughness length for the 
! landuse category water

real, public, parameter  :: &
     PRANDTL = 0.71,            &   ! Prandtl number (see Garratt, 1992)
     Sc_H20  = 0.6,             &   ! Schmidt number for water
  CHARNOCK = 0.0144  !  From Garratt for k=0.41
  !CHARNOCK = 0.032   ! Charnock's alpha:
                     ! see Nordeng (1986), p.31, 
                     ! Nordeng(1991), JGR, 96, no. C4, pp. 7167-7174.
                     ! In the second of these publications, Nordeng uses
                     ! "m" to denote Charnock's alpha whilst in the first
                     ! he specifies the value 0.032.

! Standard temperature :
real, public, parameter :: T0 = 273.15   ! zero degrees Celsius in Kelvin 
!=============================================================Y:0
endmodule PhysicalConstants_mod
