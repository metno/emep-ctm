! <local2global.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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
    subroutine local2global(locarr,gloarr,msnr)
!
!    gather a 2D 'real' array from the processors at the host me=0
!
    use ModelConstants_ml, only : NPROC ! Actual total number of processors
    use Par_ml , only : &
             MAXLIMAX &    ! Maximum number of local points in longitude
             ,MAXLJMAX&    ! Maximum number of local points in latitude
             ,GIMAX     &    ! Number of global points in longitude
             ,GJMAX     &    ! Number of global points in latitude
             ,tgi0     &    ! start points for all processors in longitude
             ,tgj0     &    ! start points for all processors in latitude
             ,tlimax     &    ! number of points for all processors in longitude
             ,tljmax     &    ! number of points for all processors in latitude
             ,me         ! Address of processor, numbering starts at 0 in south-west corner of ground level
!
    implicit none

        INCLUDE 'mpif.h'
        INTEGER STATUS(MPI_STATUS_SIZE),INFO
!
!    input
!
    integer msnr            ! message number
        real locarr(MAXLIMAX,MAXLJMAX)    ! Local array
!
!    output
!
        real gloarr(GIMAX,GJMAX)    ! Global array
!
!    local
    integer i,j,d
!
    if (me .ne. 0) then
!
!    send to host
!
           CALL MPI_SEND( locarr, 8*MAXLIMAX*MAXLJMAX, MPI_BYTE,&
                0, msnr, MPI_COMM_WORLD, INFO) 
!
    else ! me = 0
!
!    copy first local array
!
      do j = 1, tljmax(0)
        do i = 1, tlimax(0)
          gloarr(i,j) = locarr(i,j)
        enddo
      enddo
!
!    now get from the others
!
      do d = 1, NPROC-1
             CALL MPI_RECV(locarr, 8*MAXLIMAX*MAXLJMAX, MPI_BYTE, &
             d, msnr, MPI_COMM_WORLD, STATUS, INFO) 
        do j = 1, tljmax(d)
          do i = 1, tlimax(d)
        gloarr(tgi0(d)-1+i,tgj0(d)-1+j)=locarr(i,j)
          enddo
        enddo
      enddo
!
    endif    ! me = ?
!
    end
