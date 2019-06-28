! <InterpolationRoutines_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module InterpolationRoutines_mod
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

!great_circle_distance already defined in Functions_mod, but keep here for use as standalone
  private :: great_circle_distance!distance between two points following the surface on a unit sphere

  !/-- interpolation stuff
  public  :: bilin_interpolate                         !  "Generic" subroutine
  private :: bilin_interp_elem
  private :: bilin_interp_array

  public :: Nearest4interp
  public :: grid2grid_coeff,point2grid_coeff
  public :: Averageconserved_interpolate

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
    !Note: 1234 = 1432, but not= 1324 .1234 is taken anticlockwise.

    logical ::inside
    real, intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,x,y

    inside=.false.

    if(((x1-x2)*(y-y2)-(x-x2)*(y1-y2))*((x3-x4)*(y-y4)-(x-x4)*(y3-y4))>0 &
         .and. &
       ((x2-x3)*(y-y3)-(x-x3)*(y2-y3))*((x4-x1)*(y-y1)-(x-x1)*(y4-y1))>0)&
      then
       inside=.true.
    end if

  end function inside_1234

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function great_circle_distance(fi1,lambda1,fi2,lambda2) result(dist)
  ! compute the great circle distance between to points given in 
  ! spherical coordinates. Sphere has radius 1.
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

! makes interpolation coefficients from one grid to the other, 
! using only latitude and longitudes of individual gridcells

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

  real, dimension(4,NXG,NYG), target :: Weight
  integer, dimension(4,NXG,NYG) :: IIij,JJij
  real, dimension(:), pointer :: ww 
  real :: Undefined  , sumWeights
  integer :: i, j, ii, jj, k

  Undefined = -999.0e19  ! Should improve later
  if(present(Undef)) Undefined = Undef

! Get interpolation coefficients.
! Coefficients could be saved and reused if called several times.

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
    glon,glat,          &
    IIij,JJij,Weight,   & ! Returns coordinates of 4 nearest pts and weights
    dlon,dlat,NXD,NYD,NXG,NYG,NXG,NYG,&
    debug, 1, 1) !1,1 is just a crude coord, while checking

  do i=1,limax
    do j=1,ljmax
      ww => Weight(:,i,j)  
      values_grid(i,j)= 0.0
      sumWeights      = 0.0

      do k = 1, 4
        ii = IIij(k,i,j)
        jj = JJij(k,i,j)
        if ( values_data(ii,jj) > Undefined ) then
          values_grid(i,j)=values_grid(i,j)+ww(k)*values_data(ii,jj)
          sumWeights      =sumWeights      +ww(k)
        end if
      end do

      if(sumWeights>1.0e-9) then
        values_grid(i,j)= values_grid(i,j)/sumWeights
      else
        values_grid(i,j)= Undef
      end if

      end do
    end do

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

subroutine grid2grid_coeff(glon,glat,IIij,JJij,Weight,&
     dlon,dlat,NXD,NYD,NX,NY,limax,ljmax,debug,debug_li,debug_lj)

!makes interpolation coefficients from one grid to the other, 
!using only latitude and longitudes of individual gridcells

  integer, intent(in) :: NX, NY
  integer , dimension(4,NX,NY), intent(out)::IIij, JJij
  real ,intent(out), dimension(4,NX,NY) :: Weight
  real, intent(in), dimension(NX,NY) :: &
      glon,glat ! lat/long of target grid
  integer,intent(in) :: NXD,NYD ! dimension of data grid
  real, intent(in), dimension(NXD,NYD) :: &
      dlon,dlat ! lat/long of data to be interpolated
  integer,intent(in) :: limax,ljmax ! max dims needed (limax<=NX, ljmax<=NY)
  logical, intent(in) :: debug
  integer, intent(in) :: debug_li, debug_lj
  

  integer :: i,j

  !find interpolation constants
  !note that i,j are local
  !find the four closest points
  do j=1,ljmax
    do i=1,limax
      call point2grid_coeff(glon(i,j),glat(i,j),&
             IIij(:,i,j),JJij(:,i,j),Weight(:,i,j),&
             dlon,dlat,NXD,NYD,all((/debug,i==debug_li,j==debug_lj/)))
    end do
  end do
end subroutine grid2grid_coeff
subroutine point2grid_coeff(glon,glat,IIij,JJij,Weight,dlon,dlat,NXD,NYD,debug)
  real, intent(in)    :: glon,glat ! lat/long of target grid
  integer, intent(in) :: NXD,NYD ! dimension of data grid
  real, dimension(NXD,NYD), intent(in) :: &
                         dlon,dlat ! lat/long of data to be interpolated
  integer, dimension(4), intent(out) :: IIij, JJij
  real,    dimension(4), intent(out) :: Weight
  logical, intent(in) :: debug

  real :: dist(4),DD
  integer :: II,JJ,n
  dist=1.0e30
  do JJ=1,NYD
    do II=1,NXD
      !distance between data:dlon/dlat(II,JJ) and target:glon/glat
      DD=great_circle_distance(dlon(II,JJ),dlat(II,JJ),glon,glat)
      if(DD>=dist(4))cycle
      n=MINVAL([1,2,3,4],MASK=(DD<dist))
      dist(n:4)=EOSHIFT(dist(n:4),-1,BOUNDARY=DD)
      IIij(n:4)=EOSHIFT(IIij(n:4),-1,BOUNDARY=II)
      JJij(n:4)=EOSHIFT(JJij(n:4),-1,BOUNDARY=JJ)
!     if(debug) write(*,"(a,2i4,f10.3,2i4,4f9.3,4es12.3)") "DEBUG-g2g", &
!       II,JJ,DD,IIij(1),JJij(1),dlon(II,JJ),dlat(II,JJ),glon,glat,dist(:)
    end do ! II
  end do   ! JJ

  Weight(1)=1.0-3.0*dist(1)/sum(dist(1:4))
  Weight(2)=(1.0-Weight(1))*(1.0-2.0*dist(2)/sum(dist(2:4)))
  Weight(3)=(1.0-Weight(1)-Weight(2))*(1.0-dist(3)/sum(dist(3:4)))
  Weight(4)=1.0-Weight(1)-Weight(2)-Weight(3)
  if(debug) write(*,"(a,I1,2es12.3)") &
     "DEBUG-g2gFinal",0,sum(dist),sum(Weight),&
    ("DEBUG-g2gFinal",n,dist(n),Weight(n),n=1,4)
end subroutine point2grid_coeff

   subroutine Averageconserved_interpolate(Start,Endval,Average,Nvalues,i,x)
     !this routine interpolates a function, and evaluate it at i 
     !f(n), n=1...Nvalues, x=f(i)

     !The properties imposed to f are:
     !1)The average of f is equal Average, i.e. sum(f(n),n=1,Nvalues)=Average
     !2)The function varies linearly between 0.5 and Nvalues/2 and between (Nvalues+1)/2 and Nvalues+0.5
     !  Nvalues/2 is here considered as integer arithmetics, i.e. 
     !  Nvalues/2 = (Nvalues+1)/2 for Nvalues even, and Nvalues/2 = (Nvalues+1)/2 -1 for Nvalues odd.
     !3) f(Nvalues/2)=f((Nvalues+1)/2)  (still integer arithmetics)
     !4) Start=f(0.5) is the value the function would have at i = 0.5 (input)
     !5) Endval=f(Nvalues+0.5) is the value the function would have at i = Nvalues+0.5 (input)
  
     !Nvalues is the number of values n can have starting from 1(input)
     !i is the actual values of n (input)
     !x=f(i) ; x is the calculated value at n=i (output)

     real, intent(in)::Start,Endval,Average
     integer, intent(in)::Nvalues,i
     real, intent(out)::x

     real::Middle,frac
!A) Calculate the value of the function in the middle (at (N+1)/2) Middle=f((N+1)/2):

!It is not the same value for Nvalues odd or even 
     if(2*(Nvalues/2)==Nvalues)then
        Middle=2.0*Average-(Start+Endval)*0.5
     else
        Middle=(2.0*Nvalues*Average-(Nvalues-1)*(Start+Endval)*0.5)/(Nvalues+1)
     end if

!B) Evaluate the function at i

     !Check that i is within allowed boundaries:
     if(i<1.or.i>Nvalues)stop

     if(i<=Nvalues/2)then
        !linear interpolation between Start and Middle
        frac=(i-0.5)/(Nvalues/2)
        x = Start*(1.0-frac)+Middle*frac
     elseif(i>(Nvalues+1)/2)then
        !linear interpolation between Middle and End
        frac=(Nvalues+0.5-i)/(Nvalues/2)
        x = Endval*(1.0-frac)+Middle*frac        
     elseif(i==(Nvalues+1)/2)then
        !in the middle and odd number of values 
        x=Middle
     else
        !should not be possible
        stop
     end if

   end subroutine Averageconserved_interpolate


   subroutine tabulated_interpolate(x1,x2,x,fx1,fx2,fx,tab_f,Ntab)
     !interpolate a function f(x) given as tabulated values
     !f(x) defined for 0<= x <=1
     !very robust:  min(fx2,fx1) <= fx <=max(fx2,fx1)
     !but will force a "convex" interpolation
     integer, intent(in)  ::Ntab
     real,intent(in)  ::x1,x2,x !x is the value at which f(x) must be evaluated
     real,intent(in)  ::fx1,fx2
     real,intent(out) ::fx !the interpolated value of f(x)
     real,intent(in)  ::tab_f(0:Ntab) !tabulated values of f(x). x=i/Ntab, i=0,1,...Ntab
     
     integer:: ix1,ix2,ix
     real ::weight1,weight2,y,y1,y2

     ix1=nint(x1*Ntab)
     ix2=nint(x2*Ntab)
     ix=nint(x*Ntab)
     y=tab_f(ix)
     y1=tab_f(ix1)
     y2=tab_f(ix2)
     
     weight1 = (abs(y-y2)+0.5e-8)/(abs(y-y2)+abs(y-y1)+1.0e-8)
     weight2 = 1.0 - weight1 != (abs(y-y1)+0.5e-8)/(abs(y-y2)+abs(y-y1)+1.0e-8)
     
     fx = weight1*fx1 + weight2*fx2
     
   end subroutine tabulated_interpolate
      
end module InterpolationRoutines_mod

!UNCOMMENT THE FOLLOWING TO TEST-RUN THIS MODULE AS STAND-ALONE
!program test_int
!  use InterpolationRoutines_mod, only : Nearest4interp
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
