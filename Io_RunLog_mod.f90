! <Io_RunLog_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.5>
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
module Io_RunLog_mod
!_____________________________________________________________________________
! -- routines to write out to both screen and RunLog file.
! -- keep simple to avoid circularity
!_____________________________________________________________________________
use Io_Nums_mod,             only: IO_LOG
implicit none


 public :: PrintLog                   !writes message to both RunLog and unit 6

 character(len=200), public :: logtxt  !long text string used to convert mixes
                                      ! of str, float to str

contains
!-------------------------------------------------------------------------
subroutine PrintLog(txt,OutputProc,ioOption)
  character(len=*), intent(in) :: txt
  logical, intent(in), optional :: OutputProc  !typically MasterProc, me==0
  integer, intent(in), optional :: ioOption    !use for other files
  logical :: ok2print
  integer :: io
  ok2print = .true.
  if ( present(OutputProc) ) ok2print = OutputProc
  if ( ok2print) then
    io = IO_LOG
    if ( present(ioOption) ) io = ioOption
    write(*,fmt='(A)')   trim(txt)
    write(io,fmt='(A)')  trim(txt)
  end if
end subroutine PrintLog
!-------------------------------------------------------------------------

end module Io_RunLog_mod
