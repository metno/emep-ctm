! <CheckStop_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2018 met.no
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
module CheckStop_ml
! Provides routines to check for errors and if necessary close
! down the code neatly (** all processors **). 

! The generic routine CheckStopAll is defined, so that the code may be 
! stopped if:
!   (a)  errmsg   /= ok
!   (b)  int      /= 0               (e.g. iostat index after read)
!   (c)  int1     /= int2
!   (d)  string1  /= string2
!   (e)  logical  expression = true  (e.g. lu < 0 for landuse index)
!   (f)  rangeR   real outside [range(0)..range(1)]
!   (g)  rangeI   int  outside [range(0)..range(1)]

use netcdf, only: NF90_NOERR,NF90_STRERROR
use MPI_Groups_ml  , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_INTEGER&
                                     ,MPI_COMM_CALC, IERROR

implicit none

public  :: StopAll, CheckStop, CheckNC
private :: CheckStop_ok, CheckStop_okinfo, CheckStop_int1, CheckStop_int2, &
           CheckStop_str2, CheckStop_TF, CheckStop_rangeI, CheckStop_rangeR

interface CheckStop
   module procedure CheckStop_ok
   module procedure CheckStop_okinfo
   module procedure CheckStop_int1
   module procedure CheckStop_int2
   module procedure CheckStop_str2
   module procedure CheckStop_TF
   module procedure CheckStop_rangeI,CheckStop_rangeR
end interface CheckStop

contains

subroutine StopAll(errmsg)
  character(len=*), intent(in) :: errmsg
  ! Stops all processors.
  ! MPI_COMM_CALC indicates all processors, in other programs you could have
  ! different groups of processes.
  ! INFO is the error message from MPI

  if(errmsg/="ok") then
    write(*,*) "STOP-ALL ERROR: ", trim(errmsg)
    call MPI_ABORT(MPI_COMM_CALC,9,IERROR)
  end if
end subroutine StopAll

!---- Variations on CheckStop:
subroutine CheckStop_ok(errmsg)                 ! Test if errmsg /= "ok"
  character(len=*), intent(in) :: errmsg

  if(errmsg/="ok") then
   !write(*,*) "CheckStop_ok Called with:  errmsg ", errmsg
    call StopAll(errmsg)
  end if
end subroutine CheckStop_ok

subroutine CheckStop_okinfo(errmsg,infomsg)     ! Test if errmsg /= "ok"
  character(len=*), intent(in) :: errmsg
  character(len=*), intent(in) :: infomsg

  if(errmsg/="ok") then
   !write(*,*) "CheckStop_ok Called with:  errmsg ", errmsg
    write(*,*) "                          infomsg ", infomsg
    call StopAll(errmsg)
  end if
end subroutine CheckStop_okinfo

subroutine CheckStop_int1(int1,infomsg)         ! Test if int1 /= 0
  integer, intent(in)          :: int1
  character(len=*), intent(in) :: infomsg

  if(int1/=0) then
    write(*,*) "CheckStopl_int1 Called with:    int1 ", int1
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  end if
end subroutine CheckStop_int1

subroutine CheckStop_int2(int1,int2, infomsg)   ! Test if int1 /= int2
  integer, intent(in)          :: int1, int2
  character(len=*), intent(in) :: infomsg

  if(int1/=int2) then
    write(*,*) "CheckStopl_int2 Called with: int1 ", int1, " int2 ", int2
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  end if
end subroutine CheckStop_int2

subroutine CheckStop_str2(str1,str2, infomsg)   ! Test if str1 /= str2
  character(len=*), intent(in) :: str1, str2, infomsg

  if(trim(str1)/=trim(str2)) then
    write(*,*) "CheckStopl_str2 Called with: str1 ", str1, " str2 ", str2
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  end if
end subroutine CheckStop_str2

subroutine CheckStop_TF(is_error, infomsg)   ! Test expression, e.g. lu<0
  logical, intent(in)          :: is_error  
  character(len=*), intent(in) :: infomsg

  if(is_error) then
   !write(*,*) "CheckStopl_TF   Called with: logical ", is_error
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  end if
end subroutine CheckStop_TF

subroutine CheckStop_rangeR(var,vrange,infomsg)  ! test .not.(vrange(0)<=var<=vrange(1))
  real, intent(in) :: var,vrange(0:1)
  character(len=*), intent(in) :: infomsg
  character(len=*), parameter :: &
    errfmt="(A,'=',ES12.3,' is out of range ',ES12.3,'..',F6.2)"

  if(var<vrange(0).or.var>vrange(1))then
    write(*,errfmt) "CheckStopl_range: variable",var,vrange
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  end if
end subroutine CheckStop_rangeR

subroutine CheckStop_rangeI(var,vrange,infomsg)  ! test .not.(vrange(0)<=var<=vrange(1))
  integer, intent(in) :: var,vrange(0:1)
  character(len=*), intent(in) :: infomsg
  character(len=*), parameter :: &
    errfmt="(A,'=',I0,' is out of range ',I0,'..',I0)"

  if(var<vrange(0).or.var>vrange(1))then
    write(*,errfmt) "CheckStopl_range: variable",var,vrange
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  end if
end subroutine CheckStop_rangeI

subroutine CheckNC(status,errmsg)
  implicit none
  integer, intent ( in) :: status
  character(len=*), intent(in), optional :: errmsg

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    if(present(errmsg)) print *, "ERRMSG: ", trim(errmsg)
    call StopAll("Error in netcdf routine")
  end if
end subroutine CheckNC

endmodule CheckStop_ml

