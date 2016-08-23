! <Par_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
! -*- f90 -*-
                         module Par_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Define parameters and variables specific to each node and
!  the parallel data decomposition.  
!
!----------------------------------------------------------------------------
!  Erik Berge, DNMI    Roar Skaalin, SINTEF Industrial Mathematics
!  Modified to use ModelConstants_ml for domain, July 2007, ds
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
!     Parameter statements for definition of integration area
!     and maximum number of gridpoints per processor
!RESTRI
!  we try to run on a smaller domain with the same input
!  for this reason we now define additional parameters:
!  numbers of points in the 'larger' array iilardom,jjlardom
!  and coordinates of the origin of the smaller domain with
!  respect to the larger domain ismbeg,jsmbeg
!  one can run the large domain by setting  
!  gimax, gjmax to iilardom,jjlardom
!  and  ismbeg,jsmbeg   to 1
!  also we have now to distinguish mfsize for input and output: 
!  mfsizeinp,mfsizeout!!!

use CheckStop_ml,      only : CheckStop
use Io_Nums_ml,        only:  IO_LOG
use ModelConstants_ml, only : RUNDOMAIN, IIFULLDOM, JJFULLDOM, &
            MasterProc, &  ! Set true for me=0 processor
            NPROCX, NPROCY, NPROC
implicit none
private


  integer, public, save ::  &
    IRUNBEG  &!
  , JRUNBEG  &
  , GIMAX  &!Number of rundomain points in x direction
  , GJMAX  &!Number of rundomain points in y direction
  , MAXLIMAX  &! Maximum number of subdomain points in x
  , MAXLJMAX  !&! Maximum number of subdomain points in y

!     Parameter statements for the parameters used to access the tabell 
!     of neighbor processors (neighbor)
!
 integer , public, parameter :: &
     NORTH = 1   &   ! Neighbor to the north
   , SOUTH = 2   &   ! Neighbor to the south
   , EAST  = 3   &   ! Neighbor to the east
   , WEST  = 4   &   ! Neighbor to the west
   , NOPROC = -1     ! Value in neighbor when there is no neighbor in the
                          ! actual direction
!
!
!     Variables for actual number of local points, to be computed
!
 integer , public, save :: &
     limax       & ! Actual number of local points in longitude
   , ljmax         ! Actual number of local points in latitude

!
!     Variables for global address of the start and end points 
!     on each processor
!
 integer, public, save :: & 
     gi0        &  ! Global address of longitude start point (in rundomain)
   , gi1        &  ! Global address of longitude end point (in rundomain)
   , gj0        &  ! Global address of latitute start point (in rundomain)
   , gj1           ! Global address of latitude end point (in rundomain)

!
!     Variables used as loop indexes on each processor. The values are
!     derived from limax and ljmax and the position of the domain belonging
!     to this processor. See PARINIT for assignment of values.
!
 integer, public, save :: &
     li0  & ! First local index in longitude when outer boundary is excluded
   , li1  & ! Last local index in longitude when outer boundary is excluded
   , lj0  & ! First local index in latitude when outer boundary is excluded
   , lj1    ! Last local index in latitude when outer boundary is excluded
!
!
!     Variables for address of this processor
!
 integer , public, save :: &
   me   &  ! Address of processer, numbering starts at 0 in south-west corner of ground level
 , mex  &  ! Longitude address of processor, numbering  starts at 0 on the westmost boundary
 , mey     ! Latitude address of processor, numbering starts at 0 on the southmost boundary

!
!     Variable for the table of local neighbor
!

 integer, public, save, dimension(4) ::  neighbor

!
!     Tables of actual number of points and start and end points for
!     all processors
!
 integer, public, save, allocatable, dimension(:) ::  &
      tlimax, tgi0, tgi1, tljmax, tgj0, tgj1
!
!
 logical, private, parameter :: DEBUG_PAR = .false.

!------------------------------------------------------------------------------
! Define parameters used in the communication
!------------------------------------------------------------------------------
!
!     Code for broadcasting information to all nodes
!
    integer, public, parameter ::  NX_BROADCAST = -1
!
!     The different messages used in the bott version of airpol
!
    character (len=230) :: txt
    integer, public, parameter :: &
               MSG_INIT0 = 10     &
              ,MSG_INIT1 = 11     &
              ,MSG_INIT2 = 12     &
              ,MSG_INIT3 = 13     &
              ,MSG_INIT4 = 14     &
              ,MSG_INIT5 = 15     &
              ,MSG_INIT6 = 16     &
              ,MSG_INIT7 = 17     &
              ,MSG_INIT8 = 18     &
              ,MSG_INIT9 = 19     &
              ,MSG_NORTH1 = 21     &
              ,MSG_NORTH2 = 22     &
              ,MSG_EAST1 = 31     &
              ,MSG_EAST2 = 32     &
              ,MSG_SOUTH1 = 41     &
              ,MSG_SOUTH2 = 42     &
              ,MSG_WEST1 = 51     &
              ,MSG_WEST2 = 52     &
              ,MSG_TOPO1 = 61     &
              ,MSG_TOPO2 = 62  
    integer, public, parameter :: & 
               MSG_MAIN1 = 71     &
              ,MSG_MAIN2 = 72     &
              ,MSG_READ1 = 81     &
              ,MSG_READ2 = 82     &
              ,MSG_READ3 = 83     &
              ,MSG_READ4 = 84     &
              ,MSG_READ5 = 85     &
              ,MSG_READ6 = 86     & ! hb NH3emis   
              ,MSG_READ7 = 87     &
              ,MSG_READ8 = 88     & ! HDD
              ,MSG_READ9 = 89     & ! POLLEN              
              ,MSG_FIELD1 = 91     &
              ,MSG_MET1 = 101     &
              ,MSG_MET2 = 102     &
              ,MSG_MET3 = 103     &
              ,MSG_MET4 = 104     &
              ,MSG_MET5 = 105     &
              ,MSG_MET6 = 106     &
              ,MSG_PARI = 107
!-- end of eulnx.inc

 public :: parinit
 public :: Topology

 contains

    subroutine parinit(min_grids,Pole_singular)   
!
!defines size and position of subdomains
!

    implicit none
    integer, intent(in) :: min_grids ,Pole_singular 
    integer ime, imex, imey, rest,i

    !Set size of grid actually used    
    IRUNBEG = RUNDOMAIN(1)  
    JRUNBEG = RUNDOMAIN(3)  
    GIMAX = RUNDOMAIN(2)-RUNDOMAIN(1)+1 !Number of global points in longitude
    GJMAX = RUNDOMAIN(4)-RUNDOMAIN(3)+1 !Number of global points in longitude

    !Determine NPROCX, NPROCY
    NPROCY=1
    if(Pole_singular==0)then
       !try values until it fits
       NPROCX=nint(sqrt(1.0*NPROC))
       do i=1,NPROC-NPROCX
          NPROCY=NPROC/NPROCX
          if(NPROCX*NPROCY==NPROC)then
             !we found some values that divide NPROC
             exit
          endif
          NPROCX=NPROCX+1
          call CheckStop(NPROCX>NPROC,'bug in NPROCX algorithm')             
       enddo
    elseif(Pole_singular==2)then
       !2 poles. divide Y into 2 for max efficiency (load balance)
       NPROCY=2
       NPROCX=NPROC/NPROCY
       if(NPROCX*NPROCY/=NPROC)then
          NPROCY=1
          NPROCX=NPROC/NPROCY
       endif          
    else
       !1 pole. divide only in X direction for  max efficiency (load balance)
       NPROCY=1
       NPROCX=NPROC/NPROCY
    endif
 
 
    if(GJMAX/NPROCY<min_grids.or.GIMAX/NPROCX<min_grids)then
       if(MasterProc) write(*,*)'change number of processors, or rundomain ',min_grids
       call CheckStop(GJMAX/NPROCY<min_grids,'subdomains are too small in Y direction')
       call CheckStop(GIMAX/NPROCX<min_grids,'subdomains are too small in X direction')
    endif

56  format(A,I3,A,I3,A)
66  format(A,I3,A,I5,A)
    if(MasterProc)then
       write(*,56)' Using ',NPROCX*NPROCY ,' processors out of ',NPROC !may be different in future versions
       write(IO_LOG,56)' Using ',NPROCX*NPROCY ,' processors out of ',NPROC !may be different in future versions
       write(*,66)' Divided rundomain into ',NPROCX ,' X',NPROCY ,' subdomains'
       write(IO_LOG,66)' Divided rundomain into ',NPROCX ,' X',NPROCY ,' subdomains'
    endif


    MAXLIMAX = (GIMAX+NPROCX-1)/NPROCX ! Maximum number of local points in lon
    MAXLJMAX = (GJMAX+NPROCY-1)/NPROCY !&! Maximum number of local points in lat
!
!
!     Find the x-, y-, and z-addresses of the domain assigned to the
!     processor
!
!     Check if the subdomain is large enough
!
    mey = me/NPROCX
    mex = me - mey*NPROCX

    allocate(tlimax(0:NPROC-1))
    allocate(tgi0(0:NPROC-1))
    allocate(tgi1(0:NPROC-1))
    allocate(tljmax(0:NPROC-1))
    allocate(tgj0(0:NPROC-1))
    allocate(tgj1(0:NPROC-1))

!
!
!     Find the number of grid points in each direction for this processor.
!     We first try to divide the total number equally among the 
!     processors. Then the rest is distributed one by one to first processor 
!     in each direction. Here we also set the global address of the start
!     and end point in each direction.
!
!     x-direction (longitude)
!
    limax = GIMAX/NPROCX
    rest = GIMAX - limax*NPROCX
    gi0 = mex*limax + 1
    if(rest>0)then
      if(mex.eq.NPROCX-1)then
        limax = limax+1
        gi0 = gi0+rest-1
      elseif (mex < rest-1) then
        limax = limax + 1
        gi0 = gi0 + mex
      else
        gi0 = gi0 + rest-1
      endif
    endif
    gi1 = gi0 + limax - 1
!
!     y-direction (latitude)
!
    ljmax = GJMAX/NPROCY
    rest = GJMAX - ljmax*NPROCY
    gj0 = mey*ljmax + 1
    if(rest>0)then
      if(mey.eq.NPROCY-1)then
        ljmax = ljmax + 1
        gj0 = gj0 + rest-1
      elseif (mey < rest-1) then
        ljmax = ljmax + 1
        gj0 = gj0 + mey
      else
        gj0 = gj0 + rest-1
      endif
    endif
    gj1 = gj0 + ljmax - 1

          if ( DEBUG_PAR ) then
             write(*,"(a12,10i6)") "DEBUG_PAR ", me, IRUNBEG, JRUNBEG, &
                  GIMAX, GJMAX, gi0, gi1, limax, ljmax
          end if
!
!
!     Initialize the tables containing number of gridpoints and
!     addresses of start and endpoint in all directions, for all
!     processors. This is a repetition of the computations above,
!     but now for all processors.
!
    do ime = 0, NPROC-1
!
      imey = ime/NPROCX
      imex = ime - imey*NPROCX
!
!     x-direction (longitude)
!
      tlimax(ime) = GIMAX/NPROCX
      rest = GIMAX - tlimax(ime)*NPROCX
      tgi0(ime) = imex*tlimax(ime) + 1
      if(rest>0)then
        if (imex .eq. NPROCX-1) then
          tlimax(ime) = tlimax(ime) + 1
          tgi0(ime) = tgi0(ime) + rest-1
        elseif (imex < rest-1) then
          tlimax(ime) = tlimax(ime) + 1
          tgi0(ime) = tgi0(ime) + imex
        else
          tgi0(ime) = tgi0(ime) + rest-1
        endif
      endif
      tgi1(ime) = tgi0(ime) + tlimax(ime) - 1
!
!     y-direction (latitude)
!
      tljmax(ime) = GJMAX/NPROCY
      rest = GJMAX - tljmax(ime)*NPROCY
      tgj0(ime) = imey*tljmax(ime) + 1
      if(rest > 0)then
        if (imey .eq. NPROCY-1) then
          tljmax(ime) = tljmax(ime) + 1
          tgj0(ime) = tgj0(ime) + rest-1
        elseif (imey < rest-1) then
          tljmax(ime) = tljmax(ime) + 1
          tgj0(ime) = tgj0(ime) + imey
        else
          tgj0(ime) = tgj0(ime) + rest-1
        endif
      endif
      tgj1(ime) = tgj0(ime) + tljmax(ime) - 1
!
    enddo

!      The size of the grid cannot be too small.       

        call CheckStop( limax < min_grids,   &
            "Subdomain too small! Limax must be at least min_grids")
        call CheckStop( ljmax < min_grids,   &
            "Subdomain too small! Ljmax must be at least min_grids")


    end subroutine parinit


    subroutine Topology(cyclicgrid,poles)   

! Defines the neighbors and boundaries of (sub)domain
! Boundaries are defined as having coordinates 
! between 1 and li0 or between li1 and limax or
! between 1 and lj0 or between lj1 and ljmax

    implicit none
        integer, intent(in) :: cyclicgrid  ! rv2_4_1 1 if cyclic grid
       integer, intent(in) :: poles(2)  !  poles(1)=1 if North pole,
                                        !  poles(2)=1 if South pole
!
!
!     Find the x-, y-, and z-addresses of the domain assigned to the
!     processor
!
!
    mey = me/NPROCX
    mex = me - mey*NPROCX
!
!
!   Find the neighbors of this processor.
!   Allow cyclic map in i direction.
!   Do not define north and south poles as outer Boundaries
!
        lj0 = 1
        lj1 = ljmax
    if (mey > 0) then
      neighbor(SOUTH) = me-NPROCX
    else
      neighbor(SOUTH) = NOPROC
      if(poles(2)==0)lj0 = 2
    endif
    if (mey < NPROCY-1) then
      neighbor(NORTH) = me+NPROCX
    else
      neighbor(NORTH) = NOPROC
      if(poles(1)==0)lj1 = ljmax - 1
    endif
    if (mex > 0) then
      neighbor(WEST) = me-1
      li0 = 1
    else
      neighbor(WEST) = NOPROC
      li0 = 2
          if(Cyclicgrid==1)then
             neighbor(WEST) =  me+NPROCX-1
             li0 = 1
          endif
    endif
    if (mex < NPROCX-1) then
      neighbor(EAST) = me+1
      li1 = limax
    else
      neighbor(EAST) = NOPROC
      li1 = limax - 1
          if(Cyclicgrid==1)then
             neighbor(EAST) = me-NPROCX+1
             li1 = limax
          endif
    endif

    end subroutine topology

end module Par_ml
