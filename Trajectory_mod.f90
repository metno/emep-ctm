! <Trajectory_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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
module Trajectory_mod

  !----------------------------------------------------------------!
  ! The purpose of this module is to provide output along a        !
  ! trajectory. The trajectory, defined in i,j (should be          !
  ! LAt LONG), height in m and time in UTC                         !
  ! has to be given as in input file. The trajectory could be an   !
  ! an air parcel trajectory or a flight track for comparison with !
  ! aircraft measurements (MOZAIC etc)                             !
  ! WARNING!!  This module has not been used for a very long time  !
  ! and should be updated before use.                              !
  !----------------------------------------------------------------!

  use Chemfields_mod  ,    only : xn_adv
  use ChemSpecs_mod,        only: species_adv
  use GridValues_mod ,     only : glon, glat
  use Io_mod,              only : IO_AIRCR
  use MetFields_mod,       only : z_bnd,z_mid
  use Config_module , only : dt_advec,PPBINV,KMAX_BND,NPROC, METSTEP
  use MPI_Groups_mod , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_REAL8, MPI_INTEGER, MPI_LOGICAL, &
                             MPI_MIN, MPI_MAX, MPI_SUM, &
                             MPI_COMM_CALC, MPI_COMM_WORLD, MPISTATUS, IERROR, ME_MPI, NPROC_MPI
  use Par_mod   ,          only : gi0,gi1,gj0,gj1,IRUNBEG,JRUNBEG,me
  use TimeDate_mod,        only : current_date
  use SmallUtils_mod,      only : find_index

  implicit none
  private
  
  public trajectory_init
  public trajectory_in
  public trajectory_out

  integer, private, save :: iimax, iii
  integer, private, save ::  fapos(2,999)
  real, private, save ::  kfalc(999), rhour(999)
  integer, public, parameter :: &
       NFLIGHT_MAX =    10         &   ! Max. no sondes allowed
       ,FREQ_FLIGHT =    12        &   ! Interval (hrs) between outputs (not used?)
       ,NADV_FLIGHT_MAX = 10   ! Max No.  advected species that can be requested
  integer, public :: &
       NADV_FLIGHT = 0    ! No.  advected species requested in  trajectory_init

integer, public, dimension(NADV_FLIGHT_MAX) :: FLIGHT_ADV

contains

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine trajectory_init
    implicit none
    integer ::index
    iii = 1
    rhour(1) = 1.
    rhour(2) = 0.
    index=find_index("O3",species_adv(:)%name)
    if(index>0)then
       FLIGHT_ADV=index
       NADV_FLIGHT = NADV_FLIGHT + 1
    endif
  end subroutine trajectory_init

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine trajectory_in

    character*20 falc
    logical tra_exist
    integer ii

    if(current_date%seconds /= 0 .or. (mod(current_date%hour,METSTEP)/=0) )return

  !  Here month and day where trajectory is expected is hardcoded. 
    if(current_date%month == 16) then
       if(current_date%day == 1 .or. current_date%day == 16 ) then
          if(me == 0)then
             write(falc,fmt='(''tra9606'',i2.2,''.pos'')') &
                   current_date%day
             inquire(file=falc,exist=tra_exist)
             write(6,*)'trajectory exists?',tra_exist
             if(.not.tra_exist)goto 912
             open(IO_AIRCR,file=falc,status='unknown')
             ii = 0
             do while (.true.)
                ii = ii + 1
                read(IO_AIRCR,*,end=701) rhour(ii), fapos(1,ii), &
                     fapos(2,ii), kfalc(ii)
             end do
701          continue
             iimax = ii
             rhour(iimax+1) = 0.
             write(6,*) 'falcon positions  ',iimax
             write(6,*) (rhour(ii),                 &
                  fapos(1,ii), fapos(2,ii),kfalc(ii),ii=1,5)
             close(IO_AIRCR)
             open(IO_AIRCR,file='aircraft.dat',position='append')
             write(IO_AIRCR,*) 'month and day ',current_date%month&
                  ,current_date%day
             close(IO_AIRCR)
          end if
          iii = 1

!    read on node 0
912       continue
!    now distribute
          CALL MPI_BCAST( iimax ,4*1,MPI_BYTE, 0,MPI_COMM_CALC,IERROR) 
          CALL MPI_BCAST( rhour ,8*iimax+1,MPI_BYTE, 0,MPI_COMM_CALC,IERROR) 
          CALL MPI_BCAST( kfalc ,8*iimax,MPI_BYTE, 0,MPI_COMM_CALC,IERROR) 
          CALL MPI_BCAST( fapos ,4*2*iimax,MPI_BYTE, 0,MPI_COMM_CALC,IERROR) 
!    all distributed
       end if
    end if

    return
  end subroutine trajectory_in

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine trajectory_out

    real ttt, dtmil
    integer ii,jj,k,jjj
    integer :: i

!    trajectory positions
!    if(me == 0) write(6,*) 'for tidsjekk',current_date%hour    &
!        	,dt_advec,(rhour(ii),ii=1,5)

    dtmil = dt_advec/60./60.
    ttt = current_date%hour+current_date%seconds/3600.
    if (rhour(2) > rhour(1)     &
         .and. ttt+dtmil > rhour(1)) then
       do jjj = 1,10
          if(me == 0) write(6,*) 'inne i tidsjekk',    &
               iii,jjj,ttt,rhour(iii),rhour(iii+1),dt_advec,dtmil
          
          if (ttt > rhour(iii)     &
               .and. ttt < rhour(iii+1)) then
!    we have to synchronise the processors, since for next jjj(iii)
!    the aircraft can be on another processor !!!!

             CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)

             if(me == 0) write(6,*) 'inne i tidsjekk2'    &
                  ,fapos(1,iii),fapos(1,iii), ttt
             if(gi0+IRUNBEG-1 <= fapos(1,iii) .and.      &
                  gi1+IRUNBEG-1 >= fapos(1,iii) .and.   &
                  gj0+JRUNBEG-1 <= fapos(2,iii) .and.   &
                  gj1+JRUNBEG-1 >= fapos(2,iii)) then
                write(6,*) 'inne i tidsjekk3',me,kfalc(iii)
                ii = fapos(1,iii) - gi0-IRUNBEG+2
                jj = fapos(2,iii) - gj0-JRUNBEG+2
                do k = 1,KMAX_BND-1
                   if(z_bnd(ii,jj,k) > kfalc(iii) .and.     &
                        z_bnd(ii,jj,k+1) < kfalc(iii)) then
                      write(6,*) 'inne i tidsjekk4',me,kfalc(iii)
                      open(IO_AIRCR,file='aircraft.dat'     &
                           ,position='append')
                      write(IO_AIRCR,*) ttt                 &
                           ,( xn_adv( FLIGHT_ADV(i),ii,jj,k)*PPBINV,& 
                           i=1, NADV_FLIGHT),k,z_mid(ii,jj,k),&
                           glat(ii,jj),glon(ii,jj)
                      close(IO_AIRCR)
                   end if
                end do
             end if
             iii = iii + 1
          end if
          ttt = ttt + dtmil*0.1
       end do
    end if
     
    return
  end subroutine trajectory_out

   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module Trajectory_mod


