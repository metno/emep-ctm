! <Timing_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module My_Timing_mod
!----------------------------------------------------------------------------
!+
! Timing code, Variables and text array for CPU-timing 
!
! This module may be used to collect information on either system time or
! CPU time. Calling Code_timer from the external routines is the standard
! interface. If system time is required then modify Code_timer below to
! use system_clock, and declare the time variables (tim_before, tim_after, etc.)
! as integers. If CPU time is required modify Code_timer to call 
! CPU_TIME and declare the time variables as real.
!
! Code commented out or marked with !SYS is intended for system_clock
! Code commented out or marked with !CPU is intended for cpu_time
!----------------------------------------------------------------------------
use mpi,only: MPI_DOUBLE_PRECISION,MPI_MAX,MPI_MIN, MPI_COMM_WORLD,MPI_IN_PLACE
implicit none

public :: Init_timing
public :: Add_2timing   ! Calls Code_timer, adds times and descriptions to arrays
public :: Output_timing ! Outputs

integer, public, parameter  :: NTIMING_UNIMOD=39
integer, public, save       :: NTIMING=NTIMING_UNIMOD
real, public, dimension(:), allocatable, save :: &
  mytimm,   &           ! stores CPU-s
  lastptim              ! for final CPU-s
character(len=30), dimension(:), allocatable, public, save :: &
  timing                ! description

!/--- MAKE CHANGE HERE TO SWAP FROM SYSTEM_CLOCK TO SYSTEM_TIME

logical, parameter, private ::  IS_CPU_TIME = .true.   

!SYS   integer, public, save :: &   !SYS
real,    public, save :: &   !CPU
  tim_before,tim_before0,tim_after,tim_after0,tim_before1,tim_before2

contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Init_timing(ntim)
  integer, intent(in), optional :: ntim    ! set NTIMING
!SYS      integer :: iclktck,iclksec,ierr
!SYS  call system_clock(iclktck,iclksec)   ! SYS
!SYS  rclksec = 1./float(iclksec)          ! SYS

  if(present(ntim))NTIMING=ntim
  if(allocated(mytimm))   deallocate(mytimm)
  if(allocated(lastptim)) deallocate(lastptim)
  if(allocated(timing))   deallocate(timing)
  allocate(mytimm(NTIMING),lastptim(NTIMING),timing(NTIMING))

  mytimm(:) = 0.0    !CPU and SYS
  lastptim(:)=0.0    !CPU and SYS
  timing(:)  = ""
end subroutine Init_timing
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Add_2timing(n,after,before,txt)
!+ calculates CPU time and resets "before" to "after"
  integer, intent(in) :: n                  ! No (1..NTIMING)
!SYS  integer, intent(inout) :: before   ! SYS
!SYS  integer, intent(out)   :: after    ! SYS
  real,    intent(inout) :: before   ! CPU
  real,    intent(out)   :: after    ! CPU
  character(len=*), intent(in), optional :: txt

  if((n>NTIMING).or.(n<1)) return
  call Code_Timer(after)
!SYS  mytimm(n) = mytimm(n) + (after-before)*rclksec  ! SYS
  mytimm(n) = mytimm(n) +  after-before           ! CPU

  if(present(txt)) timing(n) =  txt    ! Descriptive text if wanted
  if(after<before) mytimm(n) = -999    ! WARNING CODE
  before = after
end subroutine Add_2timing
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Output_timing(io, me,np,nx,ny)
  integer, intent(in) :: io         !  i/o number
  integer, intent(in) :: me         ! number of this processor
  integer, intent(in) :: np     ! number of processors
  integer, intent(in) :: nx, ny     !  dimensions of grid
  integer :: n, IERROR

  lastptim(:) = mytimm(:)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,lastptim,NTIMING,MPI_DOUBLE_PRECISION, &
       MPI_MAX,MPI_COMM_WORLD,IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,mytimm,NTIMING,MPI_DOUBLE_PRECISION, &
       MPI_MIN,MPI_COMM_WORLD,IERROR)
  if(me==0)then
     open(io,file='Timing.out')
     write(io,"(a18,I8)") "Number of grids = ",nx*ny
     write( 6,"(a18,I8)") "Number of grids = ",nx*ny
     write(io,"(a18,2i5)")"Number of CPUs =  ", np
     write( 6,"(a18,2i5)")"Number of CPUs =  ", np
     
     write(6, fmt="(39x,a)")' min time(s)  max time(s)'
     write(io,fmt="(39x,a)")' min time(s)  max time(s)'
     do n=1,NTIMING
        if((timing(n)=="").and.(mytimm(n)==0.0)) cycle
        write(6, fmt="(a3,i3,1x,a30,2f12.4)")'tim',n,timing(n),mytimm(n),lastptim(n)
        write(io,fmt="(a3,i3,1x,a30,2f12.4)")'tim',n,timing(n),mytimm(n),lastptim(n)
     end do
     close(io)
  endif
end subroutine Output_timing
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Code_timer(call_time)
!SYS integer, intent(inout) :: call_time      !SYS
!SYS call system_clock(call_time)             !SYS
  include 'mpif.h'
  real, intent(inout) :: call_time          !CPU
! call cpu_time(call_time)                  !CPU
  call_time=MPI_WTIME()

end subroutine Code_timer
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
endmodule My_Timing_mod


