! <global2local.f90 - A component of the EMEP MSC-W Unified Eulerian
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
    subroutine global2local(gloarr,locarr,msnr&
                 ,dim0,dimi,dimj,diml,ibeg,jbeg)
!
!    distribute a 'real' array gloarr among the processors to get locarr
!    the array may have maximum 4 dimensions (i.e. snapemis)
!    , where the dimensions to be distributed, are dimi,dimj
!    the input array gloarr may be already restricted or not
!
    use Config_module, only : NPROC  ! Actual total number of processors
    use MPI_Groups_mod    , only : MPI_BYTE, MPI_COMM_CALC, MPISTATUS,IERROR
    use Par_mod , only : &
             MAXLIMAX&    ! Maximum number of local points in longitude&
             ,MAXLJMAX&    ! Maximum number of local points in latitude&
             ,tgi0    &    ! start points for all processors in longitude&
             ,tgj0    &    ! start points for all processors in latitude&
             ,tlimax    &    ! number of points for all processors in longitude&
             ,tljmax    &    ! number of points for all processors in latitude&
             ,me        ! Address of processor, numbering starts at 0 in south-west corner of ground level
!
    implicit none
!

!    input
    integer msnr        ! message number
    integer dim0&        ! first dimension, possibly = 1 (= NSECTORS for snapemis)&
             ,dimi&        ! dimension in longitude (= GIMAX or IIFULLDOM)&
             ,dimj&        ! dimension in latitude  (= GJMAX or JJFULLDOM)&
             ,diml&        ! 4th dimension, possibly = 1 (= NCMAX for snapemis)&
             ,ibeg&        ! start point of the array in longitude, = 1 for dimi = GIMAX or = IRUNBEG for dimi = IIFULLDOM&
             ,jbeg        ! start point of the array in latitude, = 1 for dimj = GJMAX or = JRUNBEG for dimj = JJFULLDOM
        real gloarr(dim0,dimi,dimj,diml)        ! Global array
!
!    output
        real locarr(dim0,MAXLIMAX,MAXLJMAX,diml)    ! Local array
!
!    local
    integer i,j,d,n0,nl
!
    if (me .ne. 0) then
!
!    receive from host
!
           CALL MPI_RECV( locarr, 8*dim0*MAXLIMAX*MAXLJMAX*diml, &
                MPI_BYTE,  0, msnr, MPI_COMM_CALC, MPISTATUS, IERROR) 

!
    else ! me = 0
!
!    first send to the others
!
      do d = 1, NPROC-1
        do nl = 1,diml
          do j = 1, tljmax(d)
        do i = 1, tlimax(d)
          do n0 = 1,dim0
            locarr(n0,i,j,nl) = gloarr(n0,tgi0(d)+ibeg-2+i&
                     ,tgj0(d)+jbeg-2+j,nl)
          end do
        end do
          end do
        end do
            CALL MPI_SEND(locarr,8*dim0*MAXLIMAX*MAXLJMAX*diml, MPI_BYTE, &
                 d, msnr, MPI_COMM_CALC, IERROR) 
      end do
!
!    now assign processor 0 itself
!
      do nl = 1,diml
        do j = 1, tljmax(0)
          do i = 1, tlimax(0)
        do n0 = 1,dim0
          locarr(n0,i,j,nl) = gloarr(n0,i+ibeg-1,j+jbeg-1,nl)
        end do
          end do
        end do
      end do
!
    end if    ! me=?
!
    return
    end
!
!
    subroutine global2local_int(gloarr,locarr,msnr&
                 ,dimi,dimj,diml,ibeg,jbeg)
!
!    distribute an 'integer' array gloarr among the processors to get locarr
!    the array may have maximum 3 dimensions (i.e. landcode)
!    , where the dimensions to be distributed, are dimi,dimj
!    the input array gloarr may be already restricted or not
!
    use Config_module, only : NPROC  ! Actual total number of processors
    use MPI_Groups_mod
    use Par_mod , only : &
             MAXLIMAX&    ! Maximum number of local points in longitude&
             ,MAXLJMAX&    ! Maximum number of local points in latitude&
             ,tgi0    &    ! start points for all processors in longitude&
             ,tgj0    &    ! start points for all processors in latitude&
             ,tlimax    &    ! number of points for all processors in longitude&
             ,tljmax    &    ! number of points for all processors in latitude&
             ,me        ! Address of processor, numbering starts at 0 in south-west corner of ground level
!
    implicit none
!
!    input
    integer msnr        ! message number
    integer dimi    &    ! dimension in longitude (= GIMAX or IIFULLDOM)&
             ,dimj&        ! dimension in latitude  (= GJMAX or JJFULLDOM)&
             ,diml    &    ! 3rd dimension, possibly = 1 (= NCMAX for landcode)&
             ,ibeg    &    ! start point of the array in longitude, = 1 for dimi = GIMAX or = IRUNBEG for dimi = IIFULLDOM&
             ,jbeg        ! start point of the array in latitude, = 1 for dimj = GJMAX or = JRUNBEG for dimj = JJFULLDOM
        integer gloarr(dimi,dimj,diml)        ! Global array
!
!    output
        integer locarr(MAXLIMAX,MAXLJMAX,diml)    ! Local array
!
!    local
    integer i,j,d,nl
!
    if (me .ne. 0) then
!
!    receive from host
!
        CALL MPI_RECV(locarr, 4*MAXLIMAX*MAXLJMAX*diml, MPI_BYTE,  0, &
        msnr, MPI_COMM_CALC,MPISTATUS, IERROR) 
!
    else ! me = 0
!
!    first send to the others
!
      do d = 1, NPROC-1
        do nl = 1,diml
          do j = 1, tljmax(d)
        do i = 1, tlimax(d)
          locarr(i,j,nl)=gloarr(tgi0(d)+ibeg-2+i&
                     ,tgj0(d)+jbeg-2+j,nl)
        end do
          end do
        end do
          CALL MPI_SEND( locarr, 4*MAXLIMAX*MAXLJMAX*diml, &
              MPI_BYTE, d, msnr, MPI_COMM_CALC, IERROR) 
      end do
!
!    now assign processor 0 itself
!
      do nl = 1,diml
        do j = 1, tljmax(0)
          do i = 1, tlimax(0)
        locarr(i,j,nl) = gloarr(i+ibeg-1,j+jbeg-1,nl)
          end do
        end do
      end do
!
    end if    ! me = ?
!
    return
    end
!
    subroutine global2local_short(gloarr,locarr,msnr&
                 ,dimi,dimj,diml,ibeg,jbeg)
!
!    distribute a 'short integer' (integer *2) array gloarr among the processors to get locarr
!    the array may have maximum 3 dimensions (i.e. landcode)
!    , where the dimensions to be distributed, are dimi,dimj
!    the input array gloarr may be already restricted or not
!

    use Config_module, only : NPROC  ! Actual total number of processors
    use MPI_Groups_mod    , only : MPI_BYTE, MPISTATUS, MPI_COMM_CALC, IERROR
    use Par_mod , only : &
             MAXLIMAX&    ! Maximum number of local points in longitude&
             ,MAXLJMAX&    ! Maximum number of local points in latitude&
             ,tgi0    &    ! start points for all processors in longitude&
             ,tgj0    &    ! start points for all processors in latitude&
             ,tlimax    &    ! number of points for all processors in longitude&
             ,tljmax    &    ! number of points for all processors in latitude&
             ,me        ! Address of processor, numbering starts at 0 in south-west corner of ground level
!
    implicit none
!
!    input
    integer msnr        ! message number
    integer dimi    &    ! dimension in longitude (= GIMAX or IIFULLDOM)&
             ,dimj    &    ! dimension in latitude  (= GJMAX or JJFULLDOM)&
             ,diml    &    ! 3rd dimension, possibly = 1 (= NCMAX for landcode)&
             ,ibeg    &    ! start point of the array in longitude, = 1 for dimi = GIMAX or = IRUNBEG for dimi = IIFULLDOM&
             ,jbeg        ! start point of the array in latitude, = 1 for dimj = GJMAX or = JRUNBEG for dimj = JJFULLDOM
        integer*2 gloarr(dimi,dimj,diml)        ! Global array
!
!    output
        integer*2 locarr(MAXLIMAX,MAXLJMAX,diml)    ! Local array
!
!    local
    integer i,j,d,nl
!
    if (me .ne. 0) then
!
!    receive from host
!

       CALL MPI_RECV(locarr, 2*MAXLIMAX*MAXLJMAX*diml, MPI_BYTE, 0,msnr,&
          MPI_COMM_CALC, MPISTATUS, IERROR)
!
    else ! me = 0
!
!    first send to the others
!
      do d = 1, NPROC-1
        do nl = 1,diml
          do j = 1, tljmax(d)
        do i = 1, tlimax(d)
          locarr(i,j,nl)=gloarr(tgi0(d)+ibeg-2+i&
                     ,tgj0(d)+jbeg-2+j,nl)
        end do
          end do
        end do

            CALL MPI_SEND(locarr, MAXLIMAX*MAXLJMAX*diml*2, MPI_BYTE, &
                 d, msnr,MPI_COMM_CALC, IERROR)

      end do
!
!    now assign processor 0 itself
!
      do nl = 1,diml
        do j = 1, tljmax(0)
          do i = 1, tlimax(0)
        locarr(i,j,nl) = gloarr(i+ibeg-1,j+jbeg-1,nl)
          end do
        end do
      end do
!
    end if    ! me = ?
!
    return
    end
