! <AeroConstants_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.6>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2025 met.no
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

   integer, parameter, public :: NSAREA_DEF = 9 ! skip SIA_F - not needed!

   type, public :: aero_t
     ! EMEP only
     character(len=15) :: EQUILIB  ='ISORROPIA'        ! 'ISORROPIA', 'EQSAM' or 'MARS' !aerosol thermodynamics 
     character(len=15) :: EQUILIB_WATER  = 'ISORROPIA' ! 'ISORROPIA', 'MARS' or 'EQSAM' !aerosol thermodynamics for PM water
     logical           :: DYNAMICS = .false.
     logical           :: INTERNALMIXED = .true.  ! sea salt assumption, only used by ISORROPIA and EQSAM
     logical           :: CATIONS = .true.        ! dust cat assumption, now only used by ISORROPIA
     real              :: RH_UPLIM_AERO = 0.98    ! RH upper limit used in thermodynamic equilibrium calls
     real              :: RH_LOLIM_AERO = 0.15    ! RH lower limit used in thermodynamic equilibrium calls to avoid div. by zero
     logical           :: ORGANIC_WATER = .true.  ! add organic matter water uptake to PM25 aerosol water
     logical           :: ThermoH2OSurfArea = .true.  ! calculate aerosol surf. area based on thermodynamics water uptake
     real              :: OM_KAPPA = 0.087         ! OM kappa hygroscopicity factor, 0.15 default from ISORROPIA II
     real              :: OM_RHO = 1400           ! aerosol density kg/m3; based on observations (Kakavas, 2023) & florou et al., 2014  
     integer           :: NSIZE = 7               ! can be removed?
     integer :: PM_F=1,SS_F=2,DU_F=3,SS_C=4,DU_C=5,SS_F_LS=6,SS_C_LS=7,PM_F_EQUI=8,PM=9  ! Will be set in GasParticleCoeffs_mod
     logical :: JUN21AERO = .false.   ! Flag to trigger ST's 2021 EQSAM and Aero tests
   end type aero_t
   type(aero_t), public, save :: AERO = aero_t()

end module AeroConstants_mod
