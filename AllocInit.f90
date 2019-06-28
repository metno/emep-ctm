! <AllocInit.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
!> MODULE AllocInits allocates arrays, and initialises to a given value, 
!! Used to prevent common errors when allocations are made without being
!! initialised. (gfortran doesn't catch this.)

module AllocInits
  use CheckStop_mod, only : CheckStop
  implicit none
  private

  public :: AllocInit

  interface AllocInit
    module procedure alloc_real_1d_init_scalar
    module procedure alloc_real_1d_init_array
    module procedure alloc_real_2d_init_scalar
    module procedure alloc_real_3d_init_scalar
    module procedure alloc_integer_1d_init_scalar
    module procedure alloc_integer_1d_init_array
  end interface AllocInit

contains

  subroutine alloc_real_1d_init_scalar(a, init, n1, txt)
    real, allocatable, dimension(:), intent(inout) :: a
    real, intent(in) :: init
    integer, intent(in) :: n1
    character(len=*), intent(in) :: txt

    integer :: istat

    if (allocated(a) .and. any(shape(a) /= [n1])) deallocate(a)
    if (.not. allocated(a)) then
      allocate(a(n1), stat=istat)
      call CheckStop(istat /= 0, "ERROR alloc_real_1d_init_scalar:"//trim(txt))
    end if
    a = init
  end subroutine alloc_real_1d_init_scalar

  subroutine alloc_real_1d_init_array(a, init, n1, txt)
    real, allocatable, dimension(:), intent(inout) :: a
    real, dimension(:), intent(in) :: init
    integer, intent(in) :: n1
    character(len=*), intent(in) :: txt

    integer :: istat

    if (allocated(a) .and. any(shape(a) /= [n1])) deallocate(a)
    if (.not. allocated(a)) then
      allocate(a(n1), stat=istat)
      call CheckStop(istat /= 0, "ERROR alloc_real_1d_init_array:"//trim(txt))
    end if
    a = init
  end subroutine alloc_real_1d_init_array

  subroutine alloc_real_2d_init_scalar(a, init, n1, n2, txt)
    real, allocatable, dimension(:,:), intent(inout) :: a
    real, intent(in) :: init
    integer, intent(in) :: n1, n2
    character(len=*), intent(in) :: txt

    integer :: istat

    if (allocated(a) .and. any(shape(a) /= [n1, n2])) deallocate(a)
    if (.not. allocated(a)) then
      allocate(a(n1, n2), stat=istat)
      call CheckStop(istat /= 0, "ERROR alloc_real_2d_init_scalar:"//trim(txt))
    end if
    a = init
  end subroutine alloc_real_2d_init_scalar

  subroutine alloc_real_3d_init_scalar(a, init, n1, n2, n3, txt)
    real, allocatable, dimension(:,:,:), intent(inout) :: a
    real, intent(in) :: init
    integer, intent(in) :: n1, n2, n3
    character(len=*), intent(in) :: txt

    integer :: istat

    if (allocated(a) .and. any(shape(a) /= [n1, n2, n3])) deallocate(a)
    if (.not. allocated(a)) then
      allocate(a(n1, n2, n3), stat=istat)
      call CheckStop(istat /= 0, "ERROR alloc_real_3d_init_scalar:"//trim(txt))
    end if
    a = init
  end subroutine alloc_real_3d_init_scalar

  subroutine alloc_integer_1d_init_scalar(a, init, n1, txt)
    integer, allocatable, dimension(:), intent(inout) :: a
    integer, intent(in) :: init
    integer, intent(in) :: n1
    character(len=*), intent(in) :: txt

    integer :: istat

    if (allocated(a) .and. any(shape(a) /= [n1])) deallocate(a)
    if (.not. allocated(a)) then
      allocate(a(n1), stat=istat)
      call CheckStop(istat /= 0, "ERROR alloc_integer_1d_init_scalar:"//trim(txt))
    end if
    a = init
  end subroutine alloc_integer_1d_init_scalar

  subroutine alloc_integer_1d_init_array(a, init, n1, txt)
    integer, allocatable, dimension(:), intent(inout) :: a
    integer, dimension(:), intent(in) :: init
    integer, intent(in) :: n1
    character(len=*), intent(in) :: txt

    integer :: istat

    if (allocated(a) .and. any(shape(a) /= [n1])) deallocate(a)
    if (.not. allocated(a)) then
      allocate(a(n1), stat=istat)
      call CheckStop(istat /= 0, "ERROR alloc_integer_1d_init_scalar:"//trim(txt))
    end if
    a = init
  end subroutine alloc_integer_1d_init_array

end module AllocInits
