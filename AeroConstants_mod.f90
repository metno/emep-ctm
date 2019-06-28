! <AeroConstants_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module AeroConstants_mod

  ! BoxChem/ESX/EMEP need some specific calculations of aerosol
  ! surface area, and we define 7 aerosol types
  ! We end up with variables here to avoid circularity in
  ! Makefile dependencies

  ! To allow emepctm-like aerosol reactions so we can refer to eg AERO%PM_F:
   
  implicit none

  private 

   integer, parameter, public :: NSAREA_DEF = 6 ! skip SIA_F - not needed!

   type, public :: aero_t
     ! EMEP only
     character(len=15) :: EQUILIB  ='MARS'  ! or 'EQSAM' !aerosol thermodynamics 
     character(len=15) :: EQUILIB_WATER  ='MARS'  ! or 'EQSAM' !aerosol thermodynamics for PM water
     logical          :: DYNAMICS = .false.
     integer          :: NSIZE    = 7
     integer :: PM_F=1,SS_F=2,DU_F=3,SS_C=4,DU_C=5,PM=6  ! Will be set in GasParticleCoeffs_mod
   end type aero_t
   type(aero_t), public, save :: AERO = aero_t()

end module AeroConstants_mod
