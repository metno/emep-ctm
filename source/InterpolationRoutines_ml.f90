! <InterpolationRoutines_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2010-2011 met.no
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
module InterpolationRoutines_ml
!____________________________________________________________________
! Miscellaneous collection of interpolation routines
! bilinear-interpolation routines,  nearest-four neighbours
!____________________________________________________________________
!
!** includes
!   bilin_interpolate - generic, elemental - guessed bilinera method
!
!   Depends on: none - self-contained.
!   Language: F
!   History:
!   ds - 2000-Jan. 2001 + great_circle_distance, nearest-4 from pw, 2010
!____________________________________________________________________
  implicit none
  private

  public :: inside_1234 !test wether a point is inside the quadrilateral 1234 

  public :: great_circle_distance!distance between two points following the surface on a unit sphere

  !/-- interpolation stuff
  public  :: bilin_interpolate                         !  "Generic" subroutine
  private :: bilin_interp_elem
  private :: bilin_interp_array

  public :: Nearest4interp
  public :: grid2grid_coeff

  real, public, dimension(0:1,0:1) :: wt    ! weighting factors, array version

  interface bilin_interpolate
     module procedure bilin_interp_array
     module procedure bilin_interp_elem
  end interface

  real, private, parameter  ::    &
       PI      = 3.141592653589793238462643383279 & ! www.verbose.net/Pi.html
    ,  DEG2RAD = PI/180.0            ! COnverts degrees to radians

  !========================================
  contains
  !========================================

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !___________________________________________________________________________
  !+ subroutines which can be used in 2-D interpolation
  !  - includes "generic" subroutine bilin_interpolate
  !___________________________________________________________________________

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine bilin_interp_array(xp,yp,ixp,iyp)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  real,    intent(in)  :: xp, yp     ! coordinates of point P (see fig.)
  integer, intent(out) :: ixp, iyp   ! integer coords P

  !/ Output:
  ! real, intent(out), dimension(0:1,0:1) :: wt    ! weights (see below)
  !-----------------------------------------------------------------------------
  ! This subroutine uses a bilinear interpolation method which suuplies the 
  ! weighting factors needed to estimate the value of a field at a point P 
  ! (input coords xp, yp) from the values at the nearest 4 grid points. 
  !
  ! This routine assumes that P is given in the coordinates of the field 
  ! which is being interpolated. If we define ixp = int(xp),iyp=int(yp),
  ! dx = xp - ixp, dy = yp - iyp,  we obtain a system: 
  !
  !        y'
  !        ^
  !        |
  !        0,1--------------------------1,1
  !        |                             |
  !        |                             |
  !        |                             |
  !        p1               *P(dx,dy)    p2
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        |                             |
  !        0,0 -------------------------1,0----------> x'
  !
  ! This subroutine outputs the weight to be given to the four corners
  ! using the array wt(0:1,0:1). 
  !
  ! For the bilinear interpolation we first calculate the weights associated
  ! with points p1,p2 along the y-axis, then interpolate these to P along the 
  ! x-axis
  !
  !  C(0,p1)  = (1-dy) * C(0,0)  + dy * C(0,1)
  !  C(1,p2)  = (1-dy) * C(1,0)  + dy * C(1,1)
  !  C(dx,dy) = (1-dx) * C(0,p1) + dx * C(1,p2)
  !           = (1-dx) * (1-dy) * C(0,0) +(1-dx) * dy * C(0,1) 
  !            +  dx   * (1-dy) * C(1,0) +   dx  * dy * C(1,1)
  !  i.e. Cp  
  !           = (1-dx-dy+dx.dy) * C(0,0)
  !            +(dy  -dx.dy)    * C(0,1)
  !            +(dx  -dx.dy)    * C(1,0)
  !            +(dx.dy)         * C(1,1)
  ! The "wt" array consists of the 4 coefficients of the C terms
  !
  ! Notes:
  !  - robust against P lying on either or both axis - no special cases are 
  !    needed.
  !  - assumes that field values exist at all corners. This is fine as long
  !    as we are using the method to interpolate from global fields.
  !-----------------------------------------------------------------------------
    real :: dx, dy, dxdy      ! local variables

    ixp = int(xp)
    iyp = int(yp)

    dx  = xp - ixp
    dy  = yp - iyp
    dxdy =dx * dy

    wt(0,0) = 1.0 - dx - dy + dxdy
    wt(0,1) = dy - dxdy
    wt(1,0) = dx - dxdy
    wt(1,1) = dxdy

  end subroutine bilin_interp_array

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  elemental subroutine bilin_interp_elem(xp,yp,ixp,iyp,wt_00,wt_01,wt_10,wt_11)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  real,    intent(in)  :: xp, yp     ! coordinates of point P (see fig.)
  integer, intent(out) :: ixp, iyp   ! integer coords P
  real, intent(out)    :: wt_00, wt_01, wt_10, wt_11  ! weights, see below

  !-----------------------------------------------------------------------------
  ! method as for subroutine bilin_interp_array, but now we return scalar
  ! arguments so that the routine can be elemental. Not quite so elegant
  ! maybe, but elemental is nice.
  !  Now we have wt_00 = wt(0,0), wt_01 = wt(0,1), etc.
  ! Note the potential for error if the arguments are not called in the correct
  ! order!
  !-----------------------------------------------------------------------------
    real :: dx, dy, dxdy      ! local variables

    ixp = int(xp)
    iyp = int(yp)

    dx   = xp - ixp
    dy   = yp - iyp
    dxdy = dx * dy

    wt_00   = 1.0 - dx - dy + dxdy
    wt_01   = dy - dxdy
    wt_10   = dx - dxdy
    wt_11   = dxdy

  end subroutine bilin_interp_elem
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function  inside_1234(x1,x2,x3,x4,y1,y2,y3,y4,x,y) result(inside)
    !test wether the point (x,y) is inside the quadrilateral 1234 
    !Note: 1234 = 1432, but not= 1324

    logical ::inside
    real, intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,x,y

    inside=.false.

    if(((x1-x2)*(y-y2)-(x-x2)*(y1-y2))*((x3-x4)*(y-y4)-(x-x4)*(y3-y4))>0 &
         .and. &
       ((x2-x3)*(y-y3)-(x-x3)*(y2-y3))*((x4-x1)*(y-y1)-(x-x1)*(y4-y1))>0)&
      then
       inside=.true.
    endif

  end function inside_1234

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function great_circle_distance(fi1,lambda1,fi2,lambda2) result(dist)

    !compute the great circle distance between to points given in 
    !spherical coordinates. Sphere has radius 1.
    real, intent(in) ::fi1,lambda1,fi2,lambda2 !NB: all in DEGREES here
    real :: dist

    dist=2*asin(sqrt(sin(DEG2RAD*0.5*(lambda1-lambda2))**2+&
         cos(DEG2RAD*lambda1)*cos(DEG2RAD*lambda2)*&
           sin(DEG2RAD*0.5*(fi1-fi2))**2))

  end function great_circle_distance

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine Nearest4interp(glon, glat, values_grid, &
      dlon,dlat,values_data, NXD,NYD,&
      NXG,NYG,limax,ljmax, debug, Undef)

   !makes interpolation coefficients from one grid to the other, 
   !using only latitude and longitudes of individual gridcells

    integer, intent(in) :: NXG, NYG
    integer, intent(in) :: NXD,NYD  ! dimension of data grid
    integer, intent(in) :: limax,ljmax  ! max dims needed (li[j]max<=NI [j]
    real, intent(in),  dimension(NXG,NYG) :: &
        glon,glat ! lat/long of target grid
    real, intent(in), dimension(NXD,NYD) :: &
        dlon,dlat ! lat/long of data to be interpolated
    real, intent(in),  dimension(NXD,NYD) :: values_data
    real, intent(out), dimension(NXG,NYG) :: values_grid
    logical, intent(in) :: debug
    real, optional, intent(in) :: Undef

  real, dimension(NXG,NYG) :: Weight1,Weight2,Weight3,Weight4
  integer, dimension(NXG,NYG,4) :: IIij,JJij
  real, dimension(4) :: Weight 
  real :: Undefined  , sumWeights
  integer :: i, j, ii, jj, k

  Undefined = -999.0e19  ! Should improve later
  if ( present(Undef) ) then
          Undefined = Undef
  end if 

!Get interpolation coefficients.
!Coefficients could be saved and reused if called several times.

     if(debug)then
         print "(a,i3,a,i3)", "Interpolate from data:", NXD, " x ", NYD
         print *, " "
         print "(9x,9f10.3)", ( dlon(i,1),i=1, NXD )  !j=1 fake
         print "(12x,a)" , "-------------------------------------------------"
         do j = NYD, 1, -1
           print "(f9.1,9f10.3)", dlat(1,j), ( values_data(i,j), i = 1, NXD)
         end do
         print "(12x,a)" , "-------------------------------------------------"
     end if

     call grid2grid_coeff( &
       glon,glat, & ! 
       IIij,JJij,    & ! Gives coordinates of 4 nearest pts
       Weight1,Weight2,Weight3,Weight4, & ! and weights
       dlon,dlat,NXD,NYD, NXG, NYG, NXG, NYG, debug, &
        1, 1) !1,1 is just a crude coord, while checking

      do i=1,limax
        do j=1,ljmax
           Weight(1) = Weight1(i,j)  
           Weight(2) = Weight2(i,j)  
           Weight(3) = Weight3(i,j)  
           Weight(4) = Weight4(i,j)  

           values_grid(i,j)= 0.0
           sumWeights      = 0.0

           do k = 1, 4
              ii = IIij(i,j,k)
              jj = JJij(i,j,k)
              if ( values_data(ii,jj) > Undefined ) then
                values_grid(i,j)= values_grid(i,j)+Weight(k)*values_data(ii,jj)
                sumWeights = sumWeights + Weight(k)
              end if

           end do

           if ( sumWeights > 1.0e-9 ) then
                  values_grid(i,j)= values_grid(i,j)/sumWeights
           else
                  values_grid(i,j)= Undef
           end if

        enddo
     enddo

    if(debug)then
         print *, " "
         print *, "To model grid:", NXG, " x ", NYG
         print *, " "
         print "(9x,9f10.3)", ( glon(i,1),i=1, NXG )  !j=1 fake
         print "(12x,a)" , "--------------------------------------------------"
         do j =  NYG, 1, -1
           print "(f9.1,9f10.3)", glat(1,j), ( values_grid(i,j), i = 1, NXG)
         end do
         print "(12x,a)" , "--------------------------------------------------"
    end if

  end subroutine Nearest4interp

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine grid2grid_coeff(glon, glat &
       ,IIij,JJij,Weight1,Weight2,Weight3,Weight4&
       ,dlon,dlat,NXD,NYD,NX,NY,limax,ljmax &
       , debug, debug_li, debug_lj)

!makes interpolation coefficients from one grid to the other, 
!using only latitude and longitudes of individual gridcells

    integer, intent(in) :: NX, NY
    integer , dimension(NX,NY,4), intent(out)::IIij, JJij
    real ,intent(out), dimension(NX,NY) :: Weight1,Weight2,Weight3,Weight4
    real, intent(in), dimension(NX,NY) :: &
        glon,glat ! lat/long of target grid
    integer,intent(in) :: NXD,NYD ! dimension of data grid
    real, intent(in), dimension(NXD,NYD) :: &
        dlon,dlat ! lat/long of data to be interpolated
    integer,intent(in) :: limax,ljmax ! max dims needed (limax<=NX, ljmax<=NY)
    logical, intent(in) :: debug
    integer, intent(in) :: debug_li, debug_lj

    real :: dist(0:4)
    integer :: i,j,II,JJ

    !find interpolation constants
    !note that i,j are local
    !find the four closest points
    do j=1,ljmax
       do i=1,limax
          dist=1.0e30
          do JJ=1,NYD
             do II=1,NXD
                !distance between (i,j) and (II,JJ)
                dist(0)=great_circle_distance(dlon(II,JJ),dlat(II,JJ),&
                    glon(i,j),glat(i,j))

                if(dist(0)<dist(1))then
                   dist(4)=dist(3)
                   dist(3)=dist(2)
                   dist(2)=dist(1)
                   dist(1)=dist(0)
                   IIij(i,j,4)=IIij(i,j,3)
                   JJij(i,j,4)=JJij(i,j,3)
                   IIij(i,j,3)=IIij(i,j,2)
                   JJij(i,j,3)=JJij(i,j,2)
                   IIij(i,j,2)=IIij(i,j,1)
                   JJij(i,j,2)=JJij(i,j,1)
                   IIij(i,j,1)=II
                   JJij(i,j,1)=JJ
                elseif(dist(0)<dist(2))then
                   dist(4)=dist(3)
                   dist(3)=dist(2)
                   dist(2)=dist(0)
                   IIij(i,j,4)=IIij(i,j,3)
                   JJij(i,j,4)=JJij(i,j,3)
                   IIij(i,j,3)=IIij(i,j,2)
                   JJij(i,j,3)=JJij(i,j,2)
                   IIij(i,j,2)=II
                   JJij(i,j,2)=JJ
                elseif(dist(0)<dist(3))then
                   dist(4)=dist(3)
                   dist(3)=dist(0)
                   IIij(i,j,4)=IIij(i,j,3)
                   JJij(i,j,4)=JJij(i,j,3)
                   IIij(i,j,3)=II
                   JJij(i,j,3)=JJ
                elseif(dist(0)<dist(4))then
                   dist(4)=dist(0)
                   IIij(i,j,4)=II
                   JJij(i,j,4)=JJ
                endif

                !if(debug.and.i==debug_li.and.j==debug_lj) then
                !    write(*,"(a,2i4,f10.3,2i4,4f9.3,4es12.3)") "DEBUG-g2g", II, JJ, dist(0),&
                !      IIij(i,j,1),JJij(i,j,1),dlon(II,JJ),dlat(II,JJ),&
                !         glon(i,j),glat(i,j), dist(1), dist(2), dist(3), dist(4)
                !    !write(*,*) "DEBUG-star", II, JJ, dist(0),&
                !    !  IIij(i,j,1),JJij(i,j,1),dlon(II,JJ),dlat(II,JJ)
                !end if
             enddo ! II
          enddo ! JJ

          dist(0)=(dist(1)+dist(2)+dist(3)+dist(4))
          Weight1(i,j)=1.0-3.0*dist(1)/dist(0)
          dist(0)=(dist(2)+dist(3)+dist(4))
          Weight2(i,j)=(1.0-Weight1(i,j))*(1.0-2.0*dist(2)/dist(0))
          dist(0)=(dist(3)+dist(4))
          Weight3(i,j)=(1.0-Weight1(i,j)-Weight2(i,j))*(1.0-dist(3)/dist(0))
          Weight4(i,j)=1.0-Weight1(i,j)-Weight2(i,j)-Weight3(i,j)
          if(debug.and.i==debug_li.and.j==debug_lj) then
             write(*,"(a,4es12.3)") "DEBUG-g2gFinal0", dist(0)
             write(*,"(a,4es12.3)") "DEBUG-g2gFinal1", dist(1), Weight1(i,j)
             write(*,"(a,4es12.3)") "DEBUG-g2gFinal2", dist(2), Weight2(i,j)
             write(*,"(a,4es12.3)") "DEBUG-g2gFinal3", dist(3), Weight3(i,j)
             write(*,"(a,4es12.3)") "DEBUG-g2gFinal4", dist(4), Weight4(i,j)
          end if

       enddo
    enddo

  end subroutine grid2grid_coeff

end module InterpolationRoutines_ml

!UNCOMMENT THE FOLLOWING TO TEST-RUN THIS MODULE AS STAND-ALONE
!program test_int
!  use InterpolationRoutines_ml, only : Nearest4interp
!  implicit none
!
!   integer, parameter :: NXG=6, NYG=3 ! fake model grid
!   integer, parameter :: NXD=4, NYD=4 ! fake data grid
!   real, dimension(NXG,NYG) :: g, glat, glon
!   real, dimension(NXD,NYD) :: d, dlat, dlon
!   integer :: ix, iy
!
!   ! Let model grid have 
!     do ix = 1, NXG
!     do iy = 1, NYG
!         glat(ix,iy) = 20.0+iy*10
!         glon(ix,iy) = -20.0+ix*10
!     end do
!     end do
!
!     do ix = 1, NXD
!     do iy = 1, NYD
!         dlat(ix,iy) =  10.0+iy*1 + 0.5 ! ix*3
!         dlon(ix,iy) = -10.0+ix*5 + 0.5 ! iy*3
!         d   (ix,iy) = -10.0+ix*5 + 0.5 ! iy*3 ! same as lon
!     end do
!     end do
!
!     print *, "TEST 1, all data present"
!
!     call Nearest4interp(glon, glat, g, &
!         dlon,dlat,d, NXD,NYD,NXG,NYG,NXG,NYG,debug=.true.)
!
!     print *, "TEST 2 ======================================================"
!     print *, "TEST 2, some undefined data"
!     where(d>6.0)
!             d = -999.0
!     end where
!
!     call Nearest4interp(glon, glat, g, &
!         dlon,dlat,d, NXD,NYD,NXG,NYG,NXG,NYG,debug=.true.,Undef=-999.0)
!
!end program test_int
