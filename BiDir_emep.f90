! <BiDir_emep.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.5>
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
module BiDir_emep
  ! DUMMY
  ! Will act as interface between emep-ctm and BiDir_module

  implicit none
  private

  public :: BiDir_ijInit
  public :: BiDir_ijRGs
  public :: BiDir_ijFluxes
  public :: BiDir_ijFinish 
  public :: BiDir_Derived

contains
 subroutine BiDir_ijInit(i,j,NH3_ix)
     integer, intent(in) :: i,j,NH3_ix
 end subroutine BiDir_ijInit
 subroutine BiDir_ijRGs(ncall,iL,Rsur,Gsto)
     integer, intent(in):: ncall, iL
     real, intent(inout) :: Rsur
     real, intent(in)  :: Gsto
 end subroutine BiDir_ijRGs
 subroutine BiDir_ijFluxes(i,j,iL,Vg_ref,Vg_eff,Rb,Rsur,Gsto)
     integer, intent(in) :: i,j,iL
     real, intent(in) :: Vg_ref,Vg_eff,Rb,Rsur,Gsto
 end subroutine BiDir_ijFluxes
 subroutine BiDir_ijFinish(i,j,gradfac,sumLand,DepLoss)
    integer, intent(in) :: i,j
    real, intent(inout) :: gradfac
    real, intent(in)    :: sumLand, DepLoss
 end subroutine BiDir_ijFinish
 subroutine BiDir_Derived(txt,n,limax,ljmax,nerr)
    character(len=*), intent(in) :: txt
    integer, intent(in) :: n,limax,ljmax
    integer, intent(inout) :: nerr
 end subroutine BiDir_Derived

end module BiDir_emep
