! <DerivedFields_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
! <DerivedFields_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************!

module DerivedFields_mod
  use OwnDataTypes_mod, only: Deriv,TXTLEN_DERIV
  implicit none
  private

  integer, public, parameter ::  &
       MAXDEF_DERIV2D = 488 & ! Max. No. 2D derived fields to be defined
      ,MAXDEF_DERIV3D = 230    ! Max. No. 3D derived fields to be defined


  ! We put definitions of **all** possible variables in def_2d, def_3d
  ! and copy the needed ones into f_xx. The data will go into d_2d, d_3d

    type(Deriv),public, dimension(MAXDEF_DERIV2D), save :: def_2d
    type(Deriv),public, dimension(MAXDEF_DERIV3D), save :: def_3d

    type(Deriv),public, allocatable, dimension(:), save :: f_2d
    type(Deriv),public, allocatable, dimension(:), save :: f_3d


  ! Fields for storing derived-style outputs. Will be allocated
  ! in Derived_mod.

  ! e.g. d_2d( num_deriv2d,LIMAX, LJMAX, LENOUT2D)
  ! &    d_3d( num_deriv3d,LIMAX, LJMAX, KMAX_MID, LENOUT3D )
   real, save, public, allocatable, dimension(:,:,:,:) :: d_2d
   real, save, public, allocatable, dimension(:,:,:,:,:) :: d_3d
   integer, save, public :: VGtest_out_ix = -1


end module DerivedFields_mod

