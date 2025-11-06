! <PBAP_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.6>
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
!> <PBAP_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!  **************************************************************************!

module PBAP_mod
  !/-- Module to deal with primary biological aerosol particles (PBAPs).
  !    DUMMY! Will be updated and released in due course
  !
  !    Gunnar Felix Lange 2024
  !---------------------------------------------------------------------------

  use CheckStop_mod,      only: CheckStop, StopAll
  use Config_module, only : USES

  implicit none

  public ::  init_PBAPs,set_PBAPs

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine init_PBAPs()
    if (USES%FUNGAL_SPORES) then
      call StopAll('ERROR: You have set USE%FUNGAL_SPORES = T, but this is not currently implemented.')
    end if
    if (USES%BACTERIA) then
      call StopAll('ERROR: You have set USE%BACTERIA = T, but this is not currently implemented.')
    end if
    if (USES%MARINE_OA) then
      call StopAll('ERROR: You have set USE%MARINE_OA = T, but this is not currently implemented.')
    end if
  end subroutine init_PBAPs


  subroutine set_PBAPs(i,j)
    integer, intent(in) ::  i,j
    call StopAll('ERROR: You are trying to call Primary Biological Aerosol Particles (PBAPs),but they are not currently implemented.')

  end subroutine set_PBAPs

end module PBAP_mod
