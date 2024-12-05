! <BiDir_module.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.5>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2024 met.no
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
module BiDir_module
  ! DUMMY 
  implicit none
  private

  character(len=*), parameter :: BiDir_module_status='TOBEDONE'

  ! Main:
  !public  :: BiDirXconcs
  !public  :: BiDirFluxes

  !public  :: BiDirResistances
  !public :: BiDirXwaterEuro
  !public :: BiDirXwaterOrig

  type, public :: BiDir_t

    logical :: EuroXwater = .false.
    logical :: OrigXwater = .false.
    logical :: skipForestDisp   = .false. 
    character(len=20) :: Method = 'NOTSET'
    ! allow for long file names
    character(len=500):: InputFile = 'NOTSET'
    character(len=500) :: InputDir  = 'NOTSET'
  end type BiDir_t
  type(BiDir_t), public, save :: BiDir= BiDir_t()

end module BiDir_module
