! <BiDir_emep.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module BiDir_emep
  ! DUMMY
  ! Will act as interface between emep-ctm and BiDir_module

  use Config_module, only : USES, MasterProc
  implicit none
  private

  character(len=*), parameter, public :: BiDir_emep_Status='TOBEDONE'

  public :: Init_BiDir  ! FUTURE

contains
  subroutine Init_BiDir()
     logical, save :: first_call = .true.
     character(len=*), parameter :: dtxt='IniBD:'
     if ( USES%BIDIR .and. MasterProc .and.  first_call) then
        write(*,*) dtxt//' FUTURE INIT'
     end if
  end subroutine Init_BiDir
end module BiDir_emep
