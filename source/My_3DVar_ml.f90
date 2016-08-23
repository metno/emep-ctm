! <My_3DVar_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_5(2809)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-201409 met.no
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
module DA_3DVar_ml
use CheckStop_ml,     only: CheckStop
use ModelConstants_ml,only: ANALYSIS
implicit none
contains
subroutine main_3dvar()
!-----------------------------------------------------------------------
! Empty call to 3dvar, for "standrd" model compilation
!-----------------------------------------------------------------------
implicit none
logical, save :: first_call=.true.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  if(.not.first_call)return
  call CheckStop(ANALYSIS,&
    "No 3DVar available. Need to recompile, e.g. make MACC-3DVar")
  first_call=.false.
endsubroutine main_3dvar
endmodule DA_3DVar_ml
