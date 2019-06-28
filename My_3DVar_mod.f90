! <My_3DVar_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module DA_mod
implicit none
logical, parameter ::     &
  DEBUG_DA_1STEP=.false.    ! run only 1 DA step (no adv/chem)
endmodule DA_mod
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module DA_3DVar_mod
use CheckStop_mod,     only: CheckStop
use Config_module,only: ANALYSIS
implicit none
character(len=*), parameter  ::  &
  mname = 'DA_3DVar_mod', &
  errmsg= 'No 3DVar available. Need to recompile, e.g. make MACC-3DVar'
integer, parameter :: NTIMING_3DVAR=0, T_3DVAR=0
contains
!-----------------------------------------------------------------------
! Empty calls, for "standrd" model compilation
!-----------------------------------------------------------------------
subroutine DA_3DVar_Init(status)
! --- in/out ----------------------------
integer, intent(out)         ::  status
! --- const -----------------------------
character(len=*), parameter  ::  rname = mname//'/DA_3DVar_Init'
! --- begin -----------------------------
  write(*,"(A,': ',A)")rname,errmsg
  status = 1
endsubroutine DA_3DVar_Init
!-----------------------------------------------------------------------
subroutine DA_3DVar_Done(status)
! --- in/out ----------------------------
integer, intent(out)         ::  status
! --- const -----------------------------
character(len=*), parameter  ::  rname = mname//'/DA_3DVar_Done'
! --- begin -----------------------------
  write(*,"(A,': ',A)")rname,errmsg
  status = 1
endsubroutine DA_3DVar_Done
!-----------------------------------------------------------------------
subroutine main_3dvar(status)
! --- in/out ----------------------------
integer, intent(out)         ::  status
! --- const -----------------------------
character(len=*), parameter  ::  rname = mname//'/main_3dvar'
! --- begin -----------------------------
  write(*,"(A,': ',A)")rname,errmsg
  status = 1
endsubroutine main_3dvar
endmodule DA_3DVar_mod
