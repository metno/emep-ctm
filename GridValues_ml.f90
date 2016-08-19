! <GridValues_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!***************************************************************************! 
!* 
!*  Copyright (C) 2007-2011 met.no
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
!***************************************************************************! 

                          Module GridValues_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Define parameters and variables associated with 3-D grid and its
!  geography.
!
! History: 
! March - changed folllwing Steffen's optimisation/correction of sigma_mid.
! January 2001 : Created by ds from old defconstants, made enough changes
! to get this into F90, and to make x,y inputs to the position subroutine,
! but the basic equations are untouched.
! October 2001 hf added call to ReadField (which now does global2local)
! Nov. 2001 - tidied up a bit (ds). Use statements moved to top of module
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 use Functions_ml,    only : great_circle_distance
                            
 use ModelConstants_ml,    only : KMAX_BND, KMAX_MID  &! vertical extent
      ,DEBUG_i, DEBUG_j  &    ! full-domain coordinate of debug-site
      ,NPROC, IIFULLDOM,JJFULLDOM, PT, Pref
 use Par_ml, only : &
        MAXLIMAX,MAXLJMAX & ! max. possible i, j in this domain
      ,limax,ljmax        & ! actual max.   i, j in this domain
      ,li0,li1,lj0,lj1    & ! for debugging TAB
      ,IRUNBEG,JRUNBEG    & ! start of user-specified domain
      ,gi0,gj0     & ! full-dom coordinates of domain lower l.h. corner
      ,me            ! local processor
 use PhysicalConstants_ml, only : GRAV, PI ! gravity, pi
 implicit none
 private

 !-- contains subroutine:

 Public :: DefGrid ! => GRIDWIDTH_M, map-factor stuff, calls other routines
 Public :: ij2lbm  ! polar stereo grid to longitude latitude
 Public :: lb2ijm  ! longitude latitude to grid in polar stereo
 Public :: ij2ijm  ! polar grid1 to polar grid2
 Public :: lb2ij   ! longitude latitude to (i,j) in any grid projection

 Public :: GlobalPosition
 private :: Position ! => lat(glat), long (glon)


  !** 1) Public (saved) Variables from module:
 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO
 real MPIbuff

  real, public, save :: &
        xp, yp  &   ! Coordinates of North pole (from infield)
      , fi      &   ! projections rotation angle around y axis (from infield)
      , AN      &   ! Distance on the map from pole to equator (No. of cells)
      ,GRIDWIDTH_M &! width of grid at 60N, in meters (old "h")(from infield)
      ,ref_latitude ! latitude at which projection is true (degrees)

!Rotated_Spherical grid prarameters
  real, public, save :: grid_north_pole_latitude,grid_north_pole_longitude

  !/ Variables to define full-domain (fdom) coordinates of local i,j values. 

  integer, public, save, dimension(0:MAXLIMAX+1) :: i_fdom !fdom coordinates
  integer, public, save, dimension(0:MAXLJMAX+1) :: j_fdom !of local i,j

 ! and reverse:
  integer, public, save, dimension(IIFULLDOM) :: i_local  !local coordinates
  integer, public, save, dimension(JJFULLDOM) :: j_local  !of full-domain i,j

!Parameters for Vertical Hybrid coordinates:
  real, public, save,  dimension(KMAX_BND) ::  &
           A_bnd !Unit Pa.  first constant, defined at layer boundary
                 ! (i.e. half levels in EC nomenclature)
  real, public, save,  dimension(KMAX_BND) ::  &
           B_bnd !Unit 1.  second constant, defined at layer boundary
                 ! (i.e. half levels in EC nomenclature)
  real, public, save,  dimension(KMAX_MID) ::  &
           A_mid !Unit Pa.  first constant, defined at middle of layer
                 ! (i.e. full levels in EC nomenclature)
  real, public, save,  dimension(KMAX_MID) ::  &
           B_mid !Unit 1.  second constant, defined at middle of layer
                 ! (i.e. full levels in EC nomenclature)
  real, public, save,  dimension(KMAX_MID) ::  &
           dA !Unit Pa.  A_bnd(k+1)-A_bnd(k) 
  real, public, save,  dimension(KMAX_MID) ::  &
           dB !Unit 1.  B_bnd(k+1)-B_bnd(k) 
! P = A + B*PS
! eta = A/Pref + B

  real, public, save,  dimension(KMAX_BND) ::  &
           sigma_bnd ! sigma, layer boundary 

  real, public, save,  dimension(KMAX_MID) ::  &
           sigma_mid   ! sigma layer midpoint

  real, public, save,  dimension(KMAX_MID) ::  carea    ! for budgets?

  real, public, save,  dimension(MAXLIMAX,MAXLJMAX) :: &
            glon       &         !longitude of gridcell centers
           ,glat                 !latitude  of gridcell centers
  real, public, save,  dimension(0:MAXLIMAX,0:MAXLJMAX) :: &
            gl_stagg   &         !longitude of gridcell corners 
           ,gb_stagg             !latitude  of gridcell corners
!NB: gl_stagg, gb_stagg are here defined as the average of the four
!    surrounding gl gb.
!    These differ slightly from the staggered points in the (i,j) grid. 

  real, public, save,  dimension(IIFULLDOM,JJFULLDOM) :: &
            glat_fdom, &         !latitude of gridcell centers
            glon_fdom            !longitude of gridcell centers


  real, public, save :: gbacmax,gbacmin,glacmax,glacmin

! EMEP grid definitions (old and official)
  real, public, parameter :: xp_EMEP_official=8.&
                            ,yp_EMEP_official=110.&
                            ,fi_EMEP=-32.&
                            ,GRIDWIDTH_M_EMEP=50000.&
                            ,an_EMEP = 237.7316364 &! = 6.370e6*(1.0+0.5*sqrt(3.0))/50000.
                            ,xp_EMEP_old =  43.0&
                            ,yp_EMEP_old = 121.0


  !/** Map factor stuff:

  real, public, save, dimension(0:MAXLIMAX+1,0:MAXLJMAX+1) ::  &
                xm_i     & ! map-factor in i direction, between cell j and j+1
               ,xm_j     & ! map-factor in j direction, between cell i and i+1
               ,xm2      & ! xm*xm: area factor in the middle of a cell (i,j)
               ,xmd        ! 1/xm2  
 
  real, public, save, dimension(0:MAXLJMAX+1,0:MAXLIMAX+1) ::  &
                xm2ji  &
               ,xmdji

  !/** Grid Area

  real, public, save, dimension(MAXLIMAX,MAXLJMAX) :: GridArea_m2

  integer, public, save :: &
              debug_li, debug_lj       ! Local Coordinates of debug-site
  logical, public, save :: debug_proc  ! Processor with debug-site

  integer, public, save :: METEOfelt=0 ! 1 if uses "old" (not CDF) meteo input


  !/** internal parameters

  logical, private, parameter ::  DEBUG_GRID = .false.  ! for debugging
  character (len=100),public::projection
  integer, public, parameter :: MIN_ADVGRIDS = 5 !minimum size of a subdomain
  integer, public :: Poles(2) !Poles(1)=1 if North pole is found, Poles(2)=1:SP
  integer, public :: Pole_included !Pole_included=1 if the grid include at least one pole

contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine DefGrid()
  !-------------------------------------------------------------------! 
  !     defines map parameters and fields for the model.
  !-------------------------------------------------------------------! 

    integer ::  i, j, k, n
    real    ::  an2, x, y, x_j, y_i
    real    ::  rpol2,rpol2_i,rpol2_j ! square of (distance from pole to i,j
                        ! divided by AN )
    real, dimension(0:MAXLIMAX+1,0:MAXLJMAX+1) ::  &
                xm      ! map-factor 

  ! Earth radius = 6.37e6 m, gives gridwidth:

!   GRIDWIDTH_M = 6.370e6*(1.0+0.5*sqrt(3.0))/AN    ! = 50000.0 m


! NB! HIRLAM uses Earth radius = 6.371e6 m : 
! AN = No. grids from pole to equator
! AN = 6.371e6*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M = 237.768957  

  !/ Define full-domain coordinates of local i,j values. We need to account for
  !  the fact that each parallel domain has its starting cordinate
  !  gi0, gj0, and the user may specify a set of lower-left starting
  !  coordinates for running the model, IRUNBEG, JRUNBEG
  !       i_fdom(i)  = i + gi0 + IRUNBEG - 2
  !       j_fdom(j)  = j + gj0 + JRUNBEG - 2

    i_fdom = (/ (n + gi0 + IRUNBEG - 2, n=0,MAXLIMAX+1) /) 
    j_fdom = (/ (n + gj0 + JRUNBEG - 2, n=0,MAXLJMAX+1) /) 

  ! And the reverse, noting that we even define for area
  ! outside local domain

    i_local = (/ (n - gi0 - IRUNBEG + 2, n=1, IIFULLDOM) /)
    j_local = (/ (n - gj0 - JRUNBEG + 2, n=1, JJFULLDOM) /)


    ! -------------- Find debug coords  and processor ------------------

    do i = 1, MAXLIMAX
       do j = 1, MAXLJMAX
          if( i_fdom(i) == DEBUG_i .and. j_fdom(j) == DEBUG_j ) then
             debug_li = i
             debug_lj = j
             debug_proc = .true.
          end if
       end do
    end do
    if( debug_proc ) write(*,*) "GridValues debug:", me, debug_li, debug_lj
    !------------------------------------------------------------------


  !  map factor, and map factor squared.

 if( METEOfelt==1)then

!mapping factor xm and ref_latitude have not been read from the meteo file

    ref_latitude=60.
    AN = 6.370e6*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M ! = 237.7316364 for GRIDWIDTH_M=50 km
    an2 = AN*AN
 
    do j=0,MAXLJMAX+1         ! ds - changed from ljmax+1
          y = j_fdom(j) - yp  ! ds - changed from gj0+JRUNBEG-2+j
          y = y*y
          y_i = j_fdom(j)+0.5 - yp  !in the staggered grid
          y_i = y_i*y_i
          do i=0,MAXLIMAX+1     
              x = i_fdom(i) - xp
              x_j = i_fdom(i)+0.5 - xp !in the staggered grid

              rpol2 = (x*x + y)/an2
              rpol2_i = (x*x + y_i)/an2
              rpol2_j = (x_j*x_j + y)/an2

              xm(i,j) = 0.5*(1.0+sin(PI/3.0))*(1.0 + rpol2)
              xm_i(i,j) = 0.5*(1.0+sin(PI/3.0))*(1.0 + rpol2_i)
              xm_j(i,j) = 0.5*(1.0+sin(PI/3.0))*(1.0 + rpol2_j)

              xm2(i,j) = xm(i,j)*xm(i,j)

          end do
      end do 
   endif
   
   AN = 6.370e6*(1.0+sin( ref_latitude*PI/180.))/GRIDWIDTH_M ! = 237.7316364 for GRIDWIDTH_M=50 km and ref_latitude=60
   do j=0,MAXLJMAX+1
      do i=0,MAXLIMAX+1
         xmd(i,j) = 1.0/xm2(i,j)
         xm2ji(j,i) = xm2(i,j)
         xmdji(j,i) = xmd(i,j)
      enddo
   enddo
   do j=1,MAXLJMAX
      do i=1,MAXLIMAX
          GridArea_m2(i,j) = GRIDWIDTH_M*GRIDWIDTH_M*xmd(i,j)
      enddo
   enddo
   
! definition of the half-sigma levels (boundaries between layers)
! from the full levels. 

   sigma_bnd(KMAX_BND) = 1.
   do k = KMAX_MID,2,-1
      sigma_bnd(k) = 2.*sigma_mid(k) - sigma_bnd(k+1)
   enddo
   sigma_bnd(1) = 0.
          
!TEMPORARY: definition of A and B. Will be read from metfile in the future
   do k = 1,KMAX_BND
      A_bnd(k)=PT * (1-sigma_bnd(k)) 
      B_bnd(k)=sigma_bnd(k) 
   enddo
   do k = 1,KMAX_MID
      dA(k)=A_bnd(k+1)-A_bnd(k)
      dB(k)=B_bnd(k+1)-B_bnd(k)
      A_mid(k)=(A_bnd(k+1)+A_bnd(k))/2.0      
      B_mid(k)=(B_bnd(k+1)+B_bnd(k))/2.0
   enddo


   !
   !     some conversion coefficients needed for budget calculations
   !
    do  k=1,KMAX_MID
        carea(k) = (sigma_bnd(k+1) - sigma_bnd(k))/GRAV*GRIDWIDTH_M*GRIDWIDTH_M
        !write(6,*)'carea,sigma_bnd,h',carea(k),sigma_bnd(k),GRIDWIDTH_M
        !su  cflux(k) = (sigma_bnd(k+1) - sigma_bnd(k))/G*h
    end do

   ! set latitude, longitude
   ! projection='Stereographic'
       call Position()

    if ( DEBUG_GRID ) then
       if ( me == 0 ) then
         write(*,800) "GRIDTAB","me","ISM","JSM","gi0","gj0",&
              "li0","li1","lix","MXI",&
              "lj0","lj1","ljx"," MXJ"," ig1"," igX","jg1","jgX"
         write(*,802) "GRIDLL ","me", "mingl"," maxgl"," mingb"," maxgb",&
              " glat(1,1)"," glat(MAX..)"
       end if

       write(*,804) "GRIDTAB",me,IRUNBEG,JRUNBEG,gi0,gj0,li0,li1,&
                    limax,MAXLIMAX,lj0,lj1,ljmax, MAXLJMAX, i_fdom(1),&
                    i_fdom(MAXLIMAX+1),j_fdom(1), j_fdom(MAXLJMAX+1)

       write(*,806) "GRIDLL ",me, minval(glon), maxval(glon), minval(glat), &
               maxval(glat), glat(1,1), glat(MAXLIMAX,MAXLJMAX)
    end if
 800 format(a10,20a4)
 802 format(a10,a4,10a12)
 804 format(a10,20i4)
 806 format(a10,i4,10f12.4)

  end subroutine DefGrid
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine Position()
  !-------------------------------------------------------------------! 
  !      calculates l(lat),b(long) (geographical coord.) 
  !      in every grid point defined by the i_fdom, j_fdom arrays. 
  !
  !      input:  xp,yp:   coord. of the polar point.
  !              AN:      number of grid-distances from pole to equator.
  !              fi:      rotational angle for the x,y grid.
  !              imax,jmax:   number of points in  x- og y- direction
  !              glmin:   gives min.value of geographical lenght
  !                       =>  glmin <= l <= glmin+360.  
  !                           (example glmin = -180. or 0.)
  !                       if "geopos","georek" is used
  !                       then glmin must be the lenght i(1,1) in the
  !                       geographical grid (gl1 to "geopos")
  !      output: gl(ii,jj): latitude glmin <= l <= glmin+360. 
  !              gb(ii,jj): longitude  -90. <= b <= +90. 
  !-------------------------------------------------------------------! 
  !   - evaluate gl, gb over whole domain given by MAXLIMAX, MAXLJMAX
  !       to safeguard against possible use of non-defined gl,bb squares. 
  !   - note, we could use rpol2(i,j) to save some computations here, 
  !       but for now we leave it. This stuff is only done once anyway

    real    :: glmin, glmax, om, om2, dy, dy2,rp,rb, rl, dx, dr
    integer :: i, j, info

      !su    xp,yp read in infield                
      !su    xp = 43.
      !su    yp = 121.

    glmin = -180.0
    glmax = glmin + 360.0

    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees
    om2   = om * 2.0


  if(trim(projection)=='Stereographic') then     
    
    do j = 1, MAXLJMAX          !  - changed from ljmax
       dy  = yp - j_fdom(j)     !  - changed from gj0+JRUNBEG-2+j
       dy2 = dy*dy
       do i = 1, MAXLIMAX       !  - changed from limax
         dx = i_fdom(i) - xp    !  - changed
         rp = sqrt(dx*dx+dy2)           ! => distance to pole
         rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
         rl = 0.0
         if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         glon(i,j)=rl                   ! longitude
         glat(i,j)=rb                   ! latitude

       end do ! i
    end do ! j
    endif

    ! test to find full-domain min and max lat/long values

    gbacmax = maxval(glat(:,:))
    gbacmin = minval(glat(:,:))
    glacmax = maxval(glon(:,:))
    glacmin = minval(glon(:,:))
    MPIbuff= gbacmax 
    CALL MPI_ALLREDUCE(MPIbuff, gbacmax, 1,MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, INFO) 
    MPIbuff= gbacmin
    CALL MPI_ALLREDUCE(MPIbuff, gbacmin  , 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO) 
    MPIbuff= glacmax 
    CALL MPI_ALLREDUCE(MPIbuff, glacmax, 1,MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, INFO) 
    MPIbuff= glacmin
    CALL MPI_ALLREDUCE(MPIbuff, glacmin  , 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO) 
    
    if(me==0) write(unit=6,fmt="(a,4f9.2)") &
                    " GridValues: max/min for lat,lon ", &
                     gbacmax,gbacmin,glacmax,glacmin

    if ( DEBUG_GRID ) then
       do j = 1, MAXLJMAX
         do i = 1, MAXLIMAX
           if ( i_fdom(i) == DEBUG_i .and. j_fdom(j) == DEBUG_j ) then
             write(*,"(a15,a30,5i4,2f8.2,f7.3)") "DEBUGPosition: ",  &
                        " me,i,j,i_fdom,j_fdom,glon,glat,rp: ", &
              me, i,j, i_fdom(i), j_fdom(j), glon(i,j), glat(i,j),rp
           end if
         end do
       end do
    end if

  end subroutine Position
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine GlobalPosition

    integer i,j
    real :: dr,om,om2,rb,rl,rp,dx,dy,dy2,glmax,glmin
    integer :: im,jm,i0,j0

  if(trim(projection)=='Stereographic') then     

    glmin = -180.0
    glmax = glmin + 360.0

    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees
    om2   = om * 2.0
    AN = 6.370e6*(1.0+sin( ref_latitude*PI/180.))/GRIDWIDTH_M ! = 237.7316364 for GRIDWIDTH_M=50 km and ref_latitude=60

    do j = 1, JJFULLDOM
       dy  = yp - j  
       dy2 = dy*dy
       do i = 1, IIFULLDOM
         dx = i - xp    
         rp = sqrt(dx*dx+dy2)           ! => distance to pole
         rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
         rl = 0.0
         if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         glon_fdom(i,j)=rl              ! longitude
         glat_fdom(i,j)=rb              ! latitude

       end do ! i
    end do ! j

    i0=0
    im=MAXLIMAX
    j0=0
    jm=MAXLJMAX
    if(gi0+MAXLIMAX+1+IRUNBEG-2>IIFULLDOM)im=MAXLIMAX-1!outside fulldomain
    if(gi0+0+IRUNBEG-2<1)i0=1!outside fulldomain
    if(gj0+MAXLJMAX+1+JRUNBEG-2>JJFULLDOM)jm=MAXLJMAX-1!outside fulldomain
    if(gj0+0+JRUNBEG-2<1)j0=1!outside fulldomain
    do j=j0,jm
       do i=i0,im
          gl_stagg(i,j)=0.25*(glon_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               glon_fdom(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               glon_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)+&
               glon_fdom(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2))
          gb_stagg(i,j)=0.25*(glat_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               glat_fdom(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               glat_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)+&
               glat_fdom(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2))
       enddo
    enddo
    do j=0,j0
       do i=i0,im
          gl_stagg(i,j)=2*gl_stagg(i,j+1)-gl_stagg(i,j+2)
          gb_stagg(i,j)=2*gb_stagg(i,j+1)-gb_stagg(i,j+2)
       enddo
    enddo
    do j=jm,MAXLJMAX
       do i=i0,im
          gl_stagg(i,j)=2*gl_stagg(i,j-1)-gl_stagg(i,j-2)
          gb_stagg(i,j)=2*gb_stagg(i,j-1)-gb_stagg(i,j-2)
       enddo
    enddo
    do j=0,MAXLJMAX
       do i=0,i0
          gl_stagg(i,j)=2*gl_stagg(i+1,j)-gl_stagg(i+2,j)
          gb_stagg(i,j)=2*gb_stagg(i+1,j)-gb_stagg(i+2,j)
       enddo
    enddo
    do j=0,MAXLJMAX
       do i=im,MAXLIMAX
          gl_stagg(i,j)=2*gl_stagg(i-1,j)-gl_stagg(i-2,j)
          gb_stagg(i,j)=2*gb_stagg(i-1,j)-gb_stagg(i-2,j)
       enddo
    enddo

 
    endif

end subroutine GlobalPosition


  subroutine lb2ijm(imax,jmax,glon,glat,xr2,yr2,fi2,an2,xp2,yp2)
    !-------------------------------------------------------------------! 
    !   calculates coordinates xr2, yr2 (real values) from gl(lat),gb(long) 
    !
    !   input: xp2,yp2:   coord. of the polar point in grid2
    !          an2:   number of grid-distances from pole to equator in grid2.
    !          fi2:      rotational angle for the grid2 (at i2=0).
    !          i1max,j1max: number of points (grid1) in  x- og y- direction
    !
    !
    !   output: i2(i1,j1): i coordinates in grid2 
    !           j2(i1,j1): j coordinates in grid2 
    !-------------------------------------------------------------------! 


    integer :: imax,jmax,i1, j1
    real    :: fi2,an2,xp2,yp2
    real    :: glon(imax,jmax),glat(imax,jmax)
    real    :: xr2(imax,jmax),yr2(imax,jmax)

    real, parameter :: PI=3.14159265358979323
    real    :: PId4,dr,dr2


    PId4    = PI/4.      
    dr2    = PI/180.0/2.      ! degrees to radians /2
    dr    = PI/180.0      ! degrees to radians 

    do j1 = 1, jmax
       do i1 = 1, imax

          xr2(i1,j1)=xp2+an2*tan(PId4-glat(i1,j1)*dr2) &
                            *sin(dr*(glon(i1,j1)-fi2))
          yr2(i1,j1)=yp2-an2*tan(PId4-glat(i1,j1)*dr2) &
                            *cos(dr*(glon(i1,j1)-fi2))

       end do ! i
    end do ! j

  end subroutine lb2ijm

 subroutine lb2ij(gl2,gb2,xr2,yr2,fi2,an2,xp2,yp2)

!Note: this routine is not supposed to be CPU optimized
 !-------------------------------------------------------------------! 
 !      calculates coordinates xr2, yr2 (real values) from gl(lat),gb(long) 
 !
 !      input:  xp2,yp2:   coord. of the polar point in grid2
 !              an2:   number of grid-distances from pole to equator in grid2.
 !              fi2:      rotational angle for the grid2 (at i2=0).
 !              i1max,j1max: number of points (grid1) in  x- og y- direction
 !
 !
 !      output: i2(i1,j1): i coordinates in grid2 
 !              j2(i1,j1): j coordinates in grid2 
 !-------------------------------------------------------------------! 

    real, intent(in)    :: gl2,gb2 
    real, intent(out)   :: xr2,yr2
    real, intent(in), optional    :: fi2,an2,xp2,yp2

    real  :: fi_loc,an_loc,xp_loc,yp_loc
    real, parameter :: PI=3.14159265358979323
    real    :: PId4,dr,dr2,dist,dist2,dist3
    integer ::i,j,ip1,jp1, ir2, jr2


  if(projection=='Stereographic')then
     PId4  =PI/4.      
     dr2   =PI/180.0/2.   ! degrees to radians /2
     dr    =PI/180.0      ! degrees to radians 
     fi_loc=fi
     an_loc=an
     xp_loc=xp
     yp_loc=yp

    if(present(fi2))fi_loc=fi2
    if(present(an2))an_loc=an2
    if(present(xp2))xp_loc=xp2
    if(present(yp2))yp_loc=yp2

    xr2=xp_loc+an_loc*tan(PId4-gb2*dr2)*sin(dr*(gl2-fi_loc))
    yr2=yp_loc-an_loc*tan(PId4-gb2*dr2)*cos(dr*(gl2-fi_loc))
  else  if(projection=='lon lat')then! lon-lat grid
     xr2=(gl2-glon_fdom(1,1))/(glon_fdom(2,1)-glon_fdom(1,1))+1
     if(xr2<0.5)xr2=xr2+360.0/(glon_fdom(2,1)-glon_fdom(1,1))
     yr2=(gb2-glat_fdom(1,1))/(glat_fdom(1,2)-glat_fdom(1,1))+1
  else!general projection, Use only info from glon_fdom and glat_fdom
     !first find closest by testing all gridcells. 
     dist=10.0!max distance is PI
     do j=1,JJFULLDOM
        do i=1,IIFULLDOM
           if(dist>great_circle_distance(gl2,gb2,glon_fdom(i,j) &
               ,glat_fdom(i,j)))then
              dist=great_circle_distance(gl2,gb2,glon_fdom(i,j) &
               ,glat_fdom(i,j))
              xr2=i
              yr2=j
           endif
        enddo
     enddo
!find the real part of i and j by comparing distances to neighbouring cells
!
!     C
!    /|\
!   / | \
!  /  |  \
! A---D---B
!
!A=(i,j) ,B=(i+1,j), C=(gl2,gb2)
!dist=AC, dist2=BC, dist3=AB
!AD=(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3)
!
     ir2 = nint(xr2)
     jr2 = nint(yr2)
     ip1=ir2+1
     if(ip1>IIFULLDOM)ip1=ip1-2
     dist2=great_circle_distance(gl2,gb2,glon_fdom(ip1,jr2),glat_fdom(ip1,jr2))
     dist3=great_circle_distance( glon_fdom(ir2,jr2), &
                                  glat_fdom(ir2,jr2), &
                                  glon_fdom(ip1,jr2), &
                                  glat_fdom(ip1,jr2))

     xr2=xr2+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)


     jp1=jr2+1
     if(jp1>JJFULLDOM)jp1=jp1-2

     dist2=great_circle_distance(gl2,gb2,glon_fdom(ir2,jp1),glat_fdom(ir2,jp1))
!GFORTRAN CHANGE
     dist3=great_circle_distance( glon_fdom(ir2,jr2), &
                                  glat_fdom(ir2,jr2), &
                                  glon_fdom(ir2,jp1), & 
                                  glat_fdom(ir2,jp1) )

     yr2=yr2+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)

  endif

    return
  end subroutine lb2ij

  subroutine ij2lbm(imax,jmax,glon,glat,fi,an,xp,yp)
  !-------------------------------------------------------------------! 
  !      calculates l(lat),b(long) (geographical coord.) 
  !      in every grid point. 
  !
  !      input:  xp,yp:   coord. of the polar point.
  !              an:      number of grid-distances from pole to equator.
  !              fi:      rotational angle for the x,y grid (at i=0).
  !              imax,jmax:   number of points in  x- og y- direction
  !              glmin:   gives min.value of geographical lenght
  !                       =>  glmin <= l <= glmin+360.  
  !                           (example glmin = -180. or 0.)
  !                       if "geopos","georek" is used
  !                       then glmin must be the lenght i(1,1) in the
  !                       geographical grid (gl1 to "geopos")
  !      output: gl(ii,jj): longitude glmin <= l <= glmin+360. 
  !              gb(ii,jj): latitude  -90. <= b <= +90. 
  !-------------------------------------------------------------------! 

    integer :: i, j, imax, jmax
    real    :: glon(imax,jmax),glat(imax,jmax)
    real    :: fi, an, xp, yp 
    real    :: om, om2, glmin, glmax,dy, dy2,rp,rb, rl, dx, dr
    real, parameter :: PI=3.14159265358979323

    glmin = -180.0

    glmax = glmin + 360.0
    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees (om=Norwegian omvendt?)
    om2   = om * 2.0

    do j = 1, jmax          
       dy  = yp - j            
       dy2 = dy*dy
       do i = 1, imax       

         dx = i - xp    ! ds - changed
         rp = sqrt(dx*dx+dy2)           ! => distance to pole
         rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
         rl = 0.0
         if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         glon(i,j)=rl                   ! longitude
         glat(i,j)=rb                   ! latitude
       end do ! i
    end do ! j

  end subroutine ij2lbm

  subroutine ij2ijm(in_field,imaxin,jmaxin,out_field,imaxout,jmaxout, &
                   fiin,anin,xpin,ypin,fiout,anout,xpout,ypout)

!   Converts data (in_field) stored in coordinates (fiin,anin,xpin,ypin) 
!   into data (out_field) in coordinates (fiout,anout,xpout,ypout)
!   pw august 2002


        integer, intent(in) :: imaxin,jmaxin,imaxout,jmaxout
        real, intent(in) :: fiin,anin,xpin,ypin,fiout,anout,xpout,ypout
        real, intent(in) :: in_field(imaxin,jmaxin)! Field to be transformed
        real, intent(out) :: out_field(imaxout,jmaxout)! Field to be transformed

        real, allocatable,dimension(:,:) :: x,y,glat,glon
        integer alloc_err,i,j,i2,j2
        logical :: interpolate
        real :: f11,f12,f21,f22

        interpolate = .true.
!        interpolate = .false.

       allocate(x(imaxout,jmaxout), stat=alloc_err)
       allocate(y(imaxout,jmaxout), stat=alloc_err)
       allocate(glat(imaxout,jmaxout), stat=alloc_err)
       allocate(glon(imaxout,jmaxout), stat=alloc_err)
       if ( alloc_err /= 0 ) WRITE(*,*) 'MPI_ABORT: ', "ij2ij alloc failed" 
       if ( alloc_err /= 0 ) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 

! find longitude, latitude of wanted area
    call ij2lbm(imaxout,jmaxout,glon,glat,fiout,anout,xpout,ypout)

! find corresponding coordinates (i,j) in in_field coordinates 
    call lb2ijm(imaxout,jmaxout,glon,glat,x,y,fiin,anin,xpin,ypin)


    ! check if the corners of the domain are inside the area covered by the 
    ! in_grid: (In principle we should test for all i,j , but test the corners
    ! should be good enough in practice) 

    if(int(x(1,1)) < 1 .or. int(x(1,1))+1 > imaxin .or. &
       int(x(imaxout,1)) < 1 .or. int(x(imaxout,1))+1 > imaxin .or. &
       int(x(1,jmaxout)) < 1 .or. int(x(1,jmaxout))+1 > imaxin .or. &
       int(x(imaxout,jmaxout)) < 1 .or. &
       int(x(imaxout,jmaxout))+1 > imaxin .or. &
       int(y(1,1)) < 1 .or. int(y(1,1))+1 > jmaxin .or. &
       int(y(imaxout,1)) < 1 .or. int(y(imaxout,1))+1 > jmaxin .or. &
       int(y(1,jmaxout)) < 1 .or. int(y(1,jmaxout))+1 > jmaxin .or. &
       int(y(imaxout,jmaxout)) < 1 .or. &
       int(y(imaxout,jmaxout))+1 > jmaxin ) then
          write(*,*)'Did not find all the necessary data in in_field'
          write(*,*)'values needed: '
          write(*,*)x(1,1),y(1,1)
          write(*,*)x(imaxout,1),y(imaxout,1)
          write(*,*)x(1,jmaxout),y(1,jmaxout)
          write(*,*)x(imaxout,jmaxout),y(imaxout,jmaxout)
          write(*,*)'max values found: ',imaxin ,jmaxin
          write(*,*) 'MPI_ABORT: ', "ij2ij: area to small" 
          call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
    endif


!  interpolate fields if required

    if(interpolate)then
    do j = 1, jmaxout
       do i = 1,imaxout
          i2=int(x(i,j))
          j2=int(y(i,j))
          f11=(1.-(x(i,j)-i2))*(1.-(y(i,j)-j2))
          f12=(1.-(x(i,j)-i2))*((y(i,j)-j2))
          f21=((x(i,j)-i2))*(1.-(y(i,j)-j2))
          f22=((x(i,j)-i2))*((y(i,j)-j2))

          out_field(i,j) =  &
               f11 * in_field(i2,j2) +  & 
               f12 * in_field(i2,j2+1) +  & 
               f21 * in_field(i2+1,j2) +  & 
               f22 * in_field(i2+1,j2+1)

       enddo
    enddo
    else

    do j = 1, jmaxout
       do i = 1,imaxout
          out_field(i,j) =in_field(nint(x(i,j)),nint(y(i,j)))
       enddo
    enddo

    endif

       deallocate(x,stat=alloc_err)
       deallocate(y,stat=alloc_err)
       deallocate(glat,stat=alloc_err)
       deallocate(glon,stat=alloc_err)
       if ( alloc_err /= 0 )   WRITE(*,*) 'MPI_ABORT: ', "ij2ijde-alloc_err" 
         if ( alloc_err /= 0 ) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
    
     end subroutine ij2ijm

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module GridValues_ml
!==============================================================================
