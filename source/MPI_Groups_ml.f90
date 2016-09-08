! <MPI_Groups_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_10(3282)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2016 met.no
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
module MPI_Groups_ml

use mpi,only: MPI_REAL8,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,&
              MPI_INTEGER,MPI_LOGICAL, MPI_BYTE,MPI_CHARACTER, &
              MPI_SUM,MPI_LOR,MPI_LAND,MPI_MAX,MPI_MIN, &
              MPI_COMM_WORLD,MPI_IN_PLACE,&
              MPI_ADDRESS_KIND,MPI_INFO_NULL,MPI_STATUS_SIZE
implicit none

integer, public, parameter ::  MasterPE = 0 ! root/master processor
integer, public,save :: MPI_COMM_CALC, MPI_COMM_IO, MPI_COMM_SUB
integer :: MPIInfo, MPISTATUS(MPI_STATUS_SIZE)
integer, public :: request_ps_w, request_ps_e, request_xn_w, request_xn_e
integer, public :: request_ps_s, request_ps_n, request_xn_s, request_xn_n
integer, public :: request_s, request_n, request_w, request_e
integer, public :: IERROR

!dummy
integer, public :: MPI_groups_split, MPI_COMM_TYPE_SHARED
integer, public, save:: ME_MPI,ME_IO,ME_CALC,ME_SUB,largeLIMAX,largeLJMAX
integer, public, save:: NPROCX_IO,NPROCY_IO,NPROC_IO,NPROCX_SUB,NPROCY_SUB,NPROC_SUB
integer, public, save:: NPROC_MPI
integer, public, save:: LargeSub_Ix

public :: MPI_world_init

contains
subroutine MPI_world_init(NPROC,ME)
  integer,intent(out) ::NPROC,ME

  CALL MPI_INIT(IERROR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, IERROR)
  ME_MPI=ME
  ME_CALC=ME
  ME_SUB=ME
  ME_IO=-1
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, IERROR)
  NPROC_MPI=NPROC
  MPI_COMM_CALC=MPI_COMM_WORLD
  MPI_COMM_SUB=MPI_COMM_WORLD
  if(ME==0)write(*,"(A,I5,A)")' Found ',NPROC,' MPI processes available'

endsubroutine MPI_world_init
subroutine share(shared_data,data_shape,xsize,MPI_COMM_SHARED)

!share the array shared_data
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER!fortran 2003 extensions 
  implicit none
  TYPE(C_PTR) :: baseptr!,baseptr2
! TYPE(MPI_Win) :: win
! TYPE(MPI_Comm), intent(in) :: MPI_COMM_SHARED
  integer :: win
  integer, intent(in) :: MPI_COMM_SHARED
  real , dimension(:,:,:),pointer, intent(inout) :: shared_data
  INTEGER(KIND=MPI_ADDRESS_KIND) :: MPI_XSIZE
  INTEGER, intent(inout) :: XSIZE
  INTEGER, intent(in) :: data_shape(3)
  INTEGER DISP_UNIT, IERROR,me_shared
  integer :: i,data_size
  real :: ONE

  DISP_UNIT=INT(sizeof(ONE),KIND(DISP_UNIT))!8, number of bytes for real

  CALL MPI_COMM_RANK(MPI_COMM_SHARED, ME_shared, IERROR)
  
  nullify(shared_data)
  MPI_XSIZE=XSIZE*DISP_UNIT
  data_size=1
  do i=1,size(data_shape)
    data_size=data_size*data_shape(i)
  enddo
  if(data_size/=XSIZE)&
    write(*,*)'WARNING: incompatible dimensions in MPI_groups_ml ',&
      data_size,XSIZE,data_shape

  if(ME_shared/=0) MPI_XSIZE = 0
  
! CALL MPI_WIN_ALLOCATE_SHARED(MPI_XSIZE, DISP_UNIT, MPI_INFO_NULL, MPI_COMM_SHARED, BASEPTR2, WIN,IERROR)

  call MPI_Win_fence(0, win, ierror)
! CALL MPI_Win_shared_query(win, 0, MPI_xsize, disp_unit, baseptr,IERROR)
  CALL C_F_POINTER(baseptr, shared_data, data_shape)
  call MPI_Win_fence(0, win, ierror)

!test if it works
! shared_data(1,1,1)=0
! shared_data(2,1,1)=0
  select case(me_mpi)
  case(1000)
    shared_data(1,1,1)=111.
  case(0)
    shared_data(1,1,1)=11.
!   if(me_io==0.and.me_sub==0)shared_data(2,1,1)=0
!   if(me_io==1.and.me_sub==0)shared_data(2,2,1)=4
!   if(me_io==1.and.me_sub==0)shared_data(1,1,1)=0.
  case(10)
    shared_data(2,1,1)=22.
  case default
    shared_data(3,1,1)=me_mpi
  endselect
  CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)
  call MPI_Win_fence(0, win, ierror)
! write(*,"(A,5i7,12F11.2)")'data in share ',ME_MPI,me_calc,me_io,me_sub,&
!  me_shared,shared_data(1:3,1,1),1.0*mpi_xsize!,shared_data(1:2,2,1)
! if(me_io>=0)write(*,*)' COMM',MPI_COMM_SHARED,MPI_COMM_WORLD

endsubroutine share
subroutine share_logical(shared_data,data_shape,xsize,MPI_COMM_SHARED)

!share the array shared_data
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER!fortran 2003 extensions 
  implicit none
  TYPE(C_PTR) :: baseptr
! TYPE(MPI_Win) :: win
! TYPE(MPI_Comm), intent(in) :: MPI_COMM_SHARED
  integer :: win
  integer, intent(in) :: MPI_COMM_SHARED
  logical , pointer, intent(inout) :: shared_data
  INTEGER(KIND=MPI_ADDRESS_KIND) :: MPI_XSIZE
  INTEGER, intent(inout) :: XSIZE
  INTEGER, intent(in) :: data_shape(1)
  INTEGER DISP_UNIT, IERROR,me_shared
  integer :: i,data_size
  logical :: mybool

  DISP_UNIT=1!number of bytes for logical
  CALL MPI_COMM_RANK(MPI_COMM_SHARED, ME_shared, IERROR)

  nullify(shared_data)
  MPI_XSIZE=XSIZE*sizeof(mybool)
  if(me_mpi==0)write(*,*)'size of logical ',sizeof(mybool)
  data_size=1
  do i=1,size(data_shape)
    data_size=data_size*data_shape(i)
  enddo
  if(data_size/=XSIZE)&
    write(*,*)'WARNING: incompatible dimensions in MPI_groups_ml ',&
      data_size,XSIZE,data_shape

! CALL MPI_WIN_ALLOCATE_SHARED(MPI_XSIZE, DISP_UNIT, MPI_INFO_NULL, MPI_COMM_SHARED, BASEPTR, WIN, IERROR)
  call MPI_Win_fence(0, win, ierror)
! CALL MPI_Win_shared_query(win, 0, MPI_xsize, disp_unit, baseptr, IERROR)
  CALL C_F_POINTER(baseptr, shared_data)
  call MPI_Win_fence(0, win, ierror)

!test if it works
  shared_data=.true.
  call MPI_Win_fence(0, win, ierror)
  if(me_io==1.and.me_sub==0)shared_data=.true.
  if(me_io==0.and.me_sub==0)shared_data=.false.
  call MPI_Win_fence(0, win, ierror)
! write(*,*)'logical data in share ',ME_MPI,me_shared,shared_data
endsubroutine share_logical
endmodule MPI_Groups_ml
