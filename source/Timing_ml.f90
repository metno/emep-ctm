! <Timing_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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
module My_Timing_ml
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
implicit none

public :: Init_timing
public :: Add_2timing   ! Calls Code_timer, adds times and descriptions to arrays
public :: Output_timing ! Outputs

integer, public, parameter                :: NTIMING=39+8
real, public, dimension(NTIMING), save    :: mytimm       ! stores CPU-s
real, public, dimension(NTIMING), save    :: lastptim     ! for final CPU-s
character(len=30), public, &
                 dimension(NTIMING), save :: timing = ""  ! description
real, private, save                       :: rclksec      ! rate-of-clock


!/--- MAKE CHANGE HERE TO SWAP FROM SYSTEM_CLOCK TO SYSTEM_TIME

logical, parameter, private ::  IS_CPU_TIME = .true.   

!SYS   integer, public, save :: &   !SYS
real,    public, save :: &   !CPU
  tim_before,tim_before0,tim_after,tim_after0,tim_before1

contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Init_timing()
!SYS      integer :: iclktck,iclksec,ierr
!SYS  call system_clock(iclktck,iclksec)   ! SYS
!SYS  rclksec = 1./float(iclksec)          ! SYS

  mytimm(:) = 0.0    !CPU and SYS
endsubroutine Init_timing
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Add_2timing(n,after,before,txt)
!+ calculates CPU time and resets "before" to "after"
  integer, intent(in) :: n                  ! No (1..NTIMING)
!SYS  integer, intent(inout) :: before   ! SYS
!SYS  integer, intent(out)   :: after    ! SYS
  real,    intent(inout) :: before   ! CPU
  real,    intent(out)   :: after    ! CPU
  character(len=*), intent(in), optional :: txt

  call Code_Timer(after)
!SYS  mytimm(n) = mytimm(n) + (after-before)*rclksec  ! SYS
  mytimm(n) = mytimm(n) +  after-before           ! CPU

  if(present(txt)) timing(n) =  txt    ! Descriptive text if wanted
  if(after<before) mytimm(n) = -999    ! WARNING CODE
  before = after
endsubroutine Add_2timing
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Output_timing(io, me,np,nt,nx,ny)
  integer, intent(in) :: io         !  i/o number
  integer, intent(in) :: me         ! number of this processor
  integer, intent(in) :: np, nt     ! number of processors, time-steps
  integer, intent(in) :: nx, ny     !  dimensions of grid
  integer :: n

  open(io,file='Timing.out')
  write(io,"(a40,I7,2i5)") "Timing for No. grids, procs, time-steps",nx*ny,np,nt
  write( 6,"(a40,I7,2i5)") "Timing for No. grids, procs, time-steps",nx*ny,np,nt

  do n=1,NTIMING
    if((timing(n)=="").and.(mytimm(n)==0.0)) cycle
    write(6, fmt="(a3,i3,1x,a30,2f12.4)")'tim',n,timing(n),mytimm(n),lastptim(n)
    write(io,fmt="(a3,i3,1x,a30,2f12.4)")'tim',n,timing(n),mytimm(n),lastptim(n)
  enddo
  close(io)
endsubroutine Output_timing
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Code_timer(call_time)
!SYS integer, intent(inout) :: call_time      !SYS
!SYS call system_clock(call_time)             !SYS
  include 'mpif.h'
  real, intent(inout) :: call_time          !CPU
! call cpu_time(call_time)                  !CPU
  call_time=MPI_WTIME()

endsubroutine Code_timer
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
endmodule My_Timing_ml


