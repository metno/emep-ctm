! <Functions_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Functions_mod
!-------------------------------------------------------------------
! Miscellaneous collection of "standard" (or guessed ) functions
! Including Troe, sine and cosine curves,
! and Standard Atmosphere p -> H conversion
!-------------------------------------------------------------------
use PhysicalConstants_mod, only: KAPPA, PI, DEG2RAD
implicit none
private

public ::  Daily_cosine   ! Generates daily values of a variable
                          ! specified as a cosine curve over the year.
public ::  Daily_sine     ! Generates daily values of a variable
                          ! specified as a sine curve over the year.
public ::  Daily_halfsine ! Similar, but only half-sine curve (0..pi)
                          ! used. (E.g. for H2O2 in ACID versions)

public :: StandardAtmos_kPa_2_km   ! US Standard Atmosphere conversion
public :: StandardAtmos_km_2_kPa   ! US Standard Atmosphere conversion

public :: great_circle_distance!distance between two points following the surface on a unit sphere

public :: heaviside ! The heaviside function, 0 for x<0 and 1 for x>0 (x==0?)

!/- Exner subroutines: ------------------------------------------------------

public :: Exner_nd        ! (p/P0)**KAPPA
public :: Tpot_2_T        ! Same as Exner_nd - but easier to remember
public :: T_2_Tpot        ! Inverse as Exner_nd
public :: Exner_tab       ! Tabulation. Must be called first

public :: ERFfunc    ! Error functions

!/- Interpolation constants

real, private, parameter  ::  &
     PINC=1000.0              &
    ,P0  =1.0e5               &    ! Standard pressure
    ,PBAS=-PINC

real, save, private, dimension(131) ::  tab_exf  ! Tabulated Exner
!-------------------------------------------------------------------
contains
!-------------------------------------------------------------------
function Daily_cosine(mean, amp, dmax, ndays) result (daily)
!+
!   Specifies cosine curve for a variable over a year
  real,    intent(in)  :: mean, amp    ! Annual mean and amplitude of sine
  integer, intent(in)  :: dmax         ! Day where maximum occurs
  integer, intent(in)  :: ndays        ! No. days per year   (365/366)

  real, dimension(ndays) :: daily
  integer    :: d
  real, save :: twopi                  ! Could use PhysiclConstants_mod
  twopi = 8.0 * atan(1.0)              ! but I prefer to keep Functions_mod
                                      ! standalone

  do d = 1, ndays
    daily(d) = mean + amp * cos ( twopi * (d - dmax)/ ndays )
  end do
end function Daily_cosine
!-------------------------------------------------------------------
function Daily_sine(mean, amp, dmax, ndays) result (daily)
!+
!   Specifies sine curve for a variable over a year
!   25/9/2002, ds, dmax redifined to be true dmax. Before it was
!   80 and actually the day when the mean ocrrurred (spotted by hf)
  real,    intent(in)  :: mean, amp    ! Annual mean and amplitude of sine
  integer, intent(in)  :: dmax         ! Day where maximum occurs
  integer, intent(in)  :: ndays        ! No. days per year   (365/366)

  real, dimension(ndays) :: daily
  integer    :: d
  real, save :: shift                  ! Shifts sine curve to give max
                                      ! when d = dmax
  real, save :: twopi                  ! Could use PhysiclConstants_mod
  twopi = 8.0 * atan(1.0)              ! but I prefer to keep Functions_mod
                                      ! standalone
  shift = ndays/4.0

  do d = 1, ndays
    daily(d) = mean + amp * sin ( twopi * (d + shift - dmax)/ ndays )
  end do
end function Daily_sine
!-------------------------------------------------------------------
function Daily_halfsine(base, amp, ndays) result (daily)
!+
!   Specifies half-sine curve for a variable over a year, with
!   values 1.0 at start and end, and max in mid-summer.
  real,    intent(in)  :: base, amp    ! Annual base and amplitude of sine
  integer, intent(in)  :: ndays        ! No. days per year   (365/366)

  real, dimension(ndays) :: daily
  integer    :: d
  real, save :: pi                     ! Could use PhysiclConstants_mod
  pi = 4.0 * atan(1.0)                 ! but I prefer to keep Functions_mod
                                        ! standalone

  do d = 1, ndays
    daily(d) = base + amp * sin ( pi * (ndays - d )/ ndays )
  end do
end function Daily_halfsine
!-------------------------------------------------------------------
elemental function StandardAtmos_km_2_kPa(h_km) result (p_kPa)
!-------------------------------------------------------------------
  implicit none
!+ Converts height (km)  to pressure (kPa) for a US standard Atmosphere
!  Valid up to 20 km
!
! pw 07/4/2010
  real :: p_kPa
  real   , intent(in)          :: h_km
! real :: t    ! Temperature (K)

  if( h_km < 11.0 ) then   ! = p_kPa > 22.632
    ! t = 288.15/(p_kPa/101.325)**(-1.0/5.255876)
    !- use the power function replacament, m**n == exp(n*log m)
    p_kPa = 101.325*exp(-5.255876*log(288.15/(288.15-6.5*h_km)))
  else
    p_kPa =  22.632*exp(-0.1576884*(h_km - 11.0)  )
  end if
end function StandardAtmos_km_2_kPa
!-------------------------------------------------------------------
elemental function StandardAtmos_kPa_2_km(p_kPa) result (h_km)
!-------------------------------------------------------------------
  implicit none
!+ Converts pressure (kPa)  to height (km) for a US standard Atmosphere
!  Valid up to 20 km
!
! ds 27/7/2003
  real, intent(in) :: p_kPa
  real             :: h_km
  real :: t    ! Temperature (K)

  if( p_kPa > 22.632 ) then   ! = 11 km height
    ! t = 288.15/(p_kPa/101.325)**(-1.0/5.255876)
    !- use the power function replacament, m**n == exp(n*log m)

    t = 288.15/exp(-1.0/5.255876*log(p_kPa/101.325))
    h_km = (288.15-t)/6.5
  else
    h_km = 11.0 + log( p_kPa/22.632)/(-0.1576884)
  end if
end function StandardAtmos_kPa_2_km
!=======================================================================
!+
!  Exner functions
!
!  Defined here as (p/P0)**KAPPA

! Where KAPPA = R/CP = 0.286
! P0 = 1.0e5 Pa

! CAREFUL:  The term Exner function  can also be used for CP * (p/P0)**KAPPA
! Hence notation Exner_nd  - non dimensional version
!
! Tabulate  :
  ! defines the exner-function for every 1000 pa from zero to 1.3e+5 pa
  ! in a table for efficient interpolation (same procedure as used in
  ! the nwp-model, see mb1e.f)
!
! Exner_nd returns the non-dimesnional excner function (p/p0)**R/CP
!
! Added 7/4/2005, Dave, based upon tpi code from tiphys
! Test prog at end along with results.
!----------------------------------------------------------------------------
subroutine Exner_tab()
  real    :: p
  integer :: i

  do i = 1,131
    p = PBAS + i*PINC
    ! tpi(i) = CP*(p/1.0e+5)**KAPPA ! With CP!!!!
    tab_exf(i) = (p/1.0e+5)**KAPPA  ! Without CP
  end do
end subroutine Exner_tab
!-------------------------------------------------------------------
elemental function Exner_nd(p) result(exf)
  real, intent(in) :: p    ! Pressure, p
  real :: exf, x1
  integer :: ix1

  x1 = (p-PBAS)/PINC
  ix1 = floor(x1)
  exf =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))
end function Exner_nd
!-------------------------------------------------------------------
elemental function Tpot_2_T(p) result(fTpot)
! Identical to Exner_nd
! Usage:   T = Tpot * Tpot_2_T(p)

  real, intent(in) :: p    ! Pressure, p
  real :: fTpot, x1
  integer :: ix1

  x1 = (p-PBAS)/PINC
  ix1 = int( x1 )
  fTpot =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))
end function Tpot_2_T
!-------------------------------------------------------------------
elemental function T_2_Tpot(p) result(fT)
! Iinvese of Exner_nd
! Usage:   Tpot = T * T_2_Tpot(p)

  real, intent(in) :: p    ! Pressure, p
  real :: fT, exf, x1
  integer :: ix1

  x1 = (p-PBAS)/PINC
  ix1 = int( x1 )
  exf =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))
  fT = 1.0/exf
end function T_2_Tpot
!-------------------------------------------------------------------
real function ERFfunc(x)
  implicit none
  ! This subprogram computes approximate values for erf(x)
  ! (see comments heading calerf).
  ! Author/Date: W. J. Cody, January 8, 1985

  real, intent (in) :: x
  integer :: jint
  real     :: result
  jint=0

  call calerf(x,result,jint)

  ERFfunc=result
end function ERFfunc
!--------------------------------------------------------------------
subroutine calerf(arg,result,jint)
!--------------------------------------------------------------------
  ! This packet evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
  ! for a real argument  x.  It contains three function type
  ! subprograms: erf, erfc, and erfcx (or derf, derfc, and derfcx),
  ! and one subroutine type subprogram, calerf.  The calling
  ! statements for the primary entries are:

  ! y=erf(x)     (or   y=derf(x)),
  ! y=erfc(x)    (or   y=derfc(x)),
  ! and
  ! y=erfcx(x)   (or   y=derfcx(x)).

  ! The routine  calerf  is intended for internal packet use only,
  ! all computations within the packet being concentrated in this
  ! routine.  The function subprograms invoke  calerf  with the
  ! statement
  ! call calerf(arg,result,jint)
  ! where the parameter usage is as follows

  ! Function                     Parameters for calerf
  ! Call              Arg                  Result          Jint
  !
  ! erf(arg)      any real argument         erf(arg)          0
  ! erfc(arg)     abs(arg)  <  xbig        erfc(arg)          1
  ! erfcx(arg)    xneg  <  arg  <  xmax   erfcx(arg)          2

  ! The main computation evaluates near-minimax approximations:
  ! from "Rational Chebyshev Approximations for the Error Function"
  ! by W. J. Cody, Math. Comp., 1969, pp. 631-638.  This
  ! transportable program uses rational functions that theoretically
  ! approximate  erf(x)  and  erfc(x)  to at least 18 significant
  ! decimal digits.  The accuracy achieved depends on the arithmetic
  ! system, the compiler, the intrinsic functions, and proper
  ! selection of the machine-dependent constants.

  ! Explanation of machine-dependent constants:
  ! xmin   = The smallest positive floating-point number.
  ! xinf   = The largest positive finite floating-point number.
  ! xneg   = The largest negative argument acceptable to erfcx;
  ! the negative of the solution to the equation
  ! 2*exp(x*x) = xinf.
  ! xsmall = Argument below which erf(x) may be represented by
  ! 2*x/sqrt(pi)  and above which  x*x  will not underflow.
  ! A conservative value is the largest machine number x
  ! such that   1.0 + x = 1.0   to machine precision.
  ! xbig   = Largest argument acceptable to erfc;  solution to
  ! the equation:  w(x)* (1-0.5/x**2) = xmin,  where
  ! w(x) = exp(-x*x)/[x*sqrt(pi)].
  ! xhuge  = Argument above which  1.0 - 1/(2*x*x) = 1.0  to
  ! machine precision.  a conservative value is
  ! 1/[2*sqrt(xsmall)]
  ! xmax   = Largest acceptable argument to erfcx; the minimum
  ! of xinf and 1/[sqrt(pi)*xmin].

  ! Approximate values for some important machines are:
  ! xmin       xinf        xneg     xsmall
  ! CDC 7600      (s.p.)  3.13e-294   1.26e+322   -27.220  7.11e-15
  ! Cray-1        (s.p.)  4.58e-2467  5.45e+2465  -75.345  7.11e-15
  ! IEEE (IBM/XT,
  ! Sun, etc.)  (s.p.)  1.18e-38    3.40e+38     -9.382  5.96e-8
  ! IEEE (IBM/XT,
  ! Sun, etc.)  (d.p.)  2.23d-308   1.79d+308   -26.628  1.11d-16
  ! IBM 195       (d.p.)  5.40d-79    7.23e+75    -13.190  1.39d-17
  ! Univac 1108   (d.p.)  2.78d-309   8.98d+307   -26.615  1.73d-18
  ! Vax d-format  (d.p.)  2.94d-39    1.70d+38     -9.345  1.39d-17
  ! Vax g-format  (d.p.)  5.56d-309   8.98d+307   -26.615  1.11d-16

  ! xbig       xhuge       xmax
  ! CDC 7600      (s.p.)  25.922      8.39e+6     1.80x+293
  ! Cray-1        (s.p.)  75.326      8.39e+6     5.45e+2465
  ! IEEE (IBM/XT,
  ! Sun, etc.)  (s.p.)   9.194      2.90e+3     4.79e+37
  ! IEEE (IBM/XT,
  ! Sun, etc.)  (d.p.)  26.543      6.71d+7     2.53d+307
  ! IBM 195       (d.p.)  13.306      1.90d+8     7.23e+75
  ! Univac 1108   (d.p.)  26.582      5.37d+8     8.98d+307
  ! Vax d-format  (d.p.)   9.269      1.90d+8     1.70d+38
  ! Vax g-format  (d.p.)  26.569      6.71d+7     8.98d+307

  ! Error returns:
  ! The program returns  erfc = 0      for  arg  >=  xbig;
  ! erfcx = xinf  for  arg  <  xneg;
  ! and
  ! erfcx = 0     for  arg  >=  xmax.

  ! Intrinsic functions required are:
  ! abs, aint, exp

  ! Author: W. J. Cody
  ! Mathematics And Computer Science Division
  ! Argonne National Laboratory
  ! Argonne, IL 60439
  ! Latest modification: March 19, 1990
  implicit none
  integer :: i,jint
  real    :: result, x, &
             arg,del,xden,xnum, y,ysq

  ! Mathematical constants
  real :: four = 4.,one = 1.,half = 0.5,two = 2.,zero = 0., &
       sqrpi = 5.6418958354775628695e-1,thresh=0.46875, &
       sixten=16.0

  ! Machine-dependent constants
  real :: xinf=3.40e+38,xneg=-9.382e0,xsmall=5.96e-8, &
       xbig=9.194, xhuge=2.90e3,xmax=4.79e37

  ! Coefficients for approximation to  erf  in first interval
  real, dimension(5) :: a =(/3.16112374387056560e00,1.13864154151050156e02, &
       3.77485237685302021e02,3.20937758913846947e03, &
       1.85777706184603153e-1/)
  real, dimension(4) :: b =(/2.36012909523441209e01,2.44024637934444173e02, &
       1.28261652607737228e03,2.84423683343917062e03/)

  ! Coefficients for approximation to  erfc  in second interval
  real, dimension(9) ::  c = &
     (/5.64188496988670089e-1, 8.88314979438837594e0, &
       6.61191906371416295e01, 2.98635138197400131e02, &
       8.81952221241769090e02, 1.71204761263407058e03, &
       2.05107837782607147e03, 1.23033935479799725e03, &
       2.15311535474403846e-8/)
  real, dimension(8) :: d = &
     (/1.57449261107098347e01,1.17693950891312499e02, &
       5.37181101862009858e02,1.62138957456669019e03, &
       3.29079923573345963e03,4.36261909014324716e03, &
       3.43936767414372164e03,1.23033935480374942e03/)

  ! Coefficients for approximation to  erfc  in third interval
  real, dimension(6) :: p =  &
     (/3.05326634961232344e-1, 3.60344899949804439e-1, &
       1.25781726111229246e-1, 1.60837851487422766e-2, &
       6.58749161529837803e-4, 1.63153871373020978e-2/)
   real, dimension(5) :: q =  &
     (/2.56852019228982242e0 ,1.87295284992346047e0 , &
       5.27905102951428412e-1,6.05183413124413191e-2, &
       2.33520497626869185e-3/)

  ! Main Code
  x=arg
  y=abs(x)
  if (y <= thresh) then
     ! Evaluate  erf  for  |x| <= 0.46875
     ysq=zero
     if (y > xsmall) ysq=y*y
     xnum=a(5)*ysq
     xden=ysq
     do i=1,3
        xnum=(xnum+a(i))*ysq
        xden=(xden+b(i))*ysq
     end do
     result=x*(xnum+a(4))/(xden+b(4))
     if (jint /= 0) result=one-result
     if (jint == 2) result=exp(ysq)*result
     go to 800
     ! Evaluate  erfc  for 0.46875 <= |x| <= 4.0
  else if (y <= four) then
     xnum=c(9)*y
     xden=y
     do i=1,7
        xnum=(xnum+c(i))*y
        xden=(xden+d(i))*y
     end do
     result=(xnum+c(8))/(xden+d(8))
     if (jint /= 2) then
        ysq=aint(y*sixten)/sixten
        del=(y-ysq)*(y+ysq)
        result=exp(-ysq*ysq)*exp(-del)*result
     end if
     ! Evaluate  erfc  for |x| > 4.0
  else
     result=zero
     if (y >= xbig) then
        if ((jint /= 2).or.(y >= xmax)) go to 300
        if (y >= xhuge) then
           result=sqrpi/y
           go to 300
        end if
     end if
     ysq=one/(y*y)
     xnum=p(6)*ysq
     xden=ysq
     do i=1,4
        xnum=(xnum+p(i))*ysq
        xden=(xden+q(i))*ysq
     end do
     result=ysq*(xnum+p(5))/(xden+q(5))
     result=(sqrpi-result)/y
     if (jint /= 2) then
        ysq=aint(y*sixten)/sixten
        del=(y-ysq)*(y+ysq)
        result=exp(-ysq*ysq)*exp(-del)*result
     end if
  end if
  ! Fix up for negative argument, erf, etc.
300 if (jint == 0) then
     result=(half-result)+half
     if (x < zero) result=-result
  else if (jint == 1) then
     if (x < zero) result=two-result
  else
     if (x < zero) then
        if (x < xneg) then
           result=xinf
        else
           ysq=aint(x*sixten)/sixten
           del=(x-ysq)*(x+ysq)
           y=exp(ysq*ysq)*exp(del)
           result=(y+y)-result
        end if
     end if
  end if
800 return
end subroutine calerf
!-------------------------------------------------------------------
PURE function great_circle_distance(fi1,lambda1,fi2,lambda2) result(dist)
!compute the great circle distance between to points given in
!spherical coordinates. Sphere has radius 1.
  real, intent(in) ::fi1,lambda1,fi2,lambda2 !NB: in DEGREES here
  real :: dist

  !ds sind not allowed in gfortran, so replaced. Also, 360 removed
  !dist=2*asin(sqrt(sind(0.5*(lambda1-lambda2+360.0))**2+&
  !     cosd(lambda1+360.0)*cosd(lambda2+360.0)*sind(0.5*(fi1-fi2+360.0))**2))

  dist=2*asin(sqrt(sin(DEG2RAD*0.5*(lambda1-lambda2))**2+&
       cos(DEG2RAD*lambda1)*cos(DEG2RAD*lambda2)*&
         sin(DEG2RAD*0.5*(fi1-fi2))**2))
end function great_circle_distance
!-----------------------------------------------------------------------
! The heaviside function, 0 for x<0 and 1 for x>0 (x==0?)
! For x=0, one could have 0.5, but numerically this is too tricky to code
! and with double precision a very rare event.
function heaviside(x)
 real, intent(in) :: x
 real             :: heaviside

 if(x<0) then
   heaviside = 0.0
 else
   heaviside = 1.0
 end if
end function heaviside
!-----------------------------------------------------------------------
!program Test_exn
!  use Exner_mod
!  use PhysicalConstants_mod, only : KAPPA
!  implicit none
!
!  real :: p, exf1, exf2
!  integer :: i
!
!  call Exner_tab()
!
!   do i = 1, 20
!     p = 0.05 * i * 1.0e5
!     exf1 = Exner_nd(p)
!     exf2 = (p*1.0e-5)**KAPPA
!     print "(f8.3,4f12.5)", 1.0e-2*p, exf1, exf2, Tpot_2_T(p), T_2_Tpot(p)
!  end do
!end program Test_exn

! Results:
! p(mb)        exf1         exf2      Tpot_2_T    T_2_Tpot
!  50.000     0.42471     0.42471     0.42471     2.35455
! 100.000     0.51778     0.51778     0.51778     1.93133
! 150.000     0.58141     0.58141     0.58141     1.71997
! 200.000     0.63124     0.63124     0.63124     1.58418
! 250.000     0.67282     0.67282     0.67282     1.48629
! 300.000     0.70881     0.70881     0.70881     1.41081
! 350.000     0.74075     0.74075     0.74075     1.34999
! 400.000     0.76957     0.76957     0.76957     1.29943
! 450.000     0.79592     0.79592     0.79592     1.25641
! 500.000     0.82025     0.82025     0.82025     1.21913
! 550.000     0.84291     0.84291     0.84291     1.18637
! 600.000     0.86414     0.86414     0.86414     1.15722
! 650.000     0.88414     0.88414     0.88414     1.13105
! 700.000     0.90307     0.90307     0.90307     1.10734
! 750.000     0.92105     0.92105     0.92105     1.08571
! 800.000     0.93820     0.93820     0.93820     1.06587
! 850.000     0.95461     0.95461     0.95461     1.04755
! 900.000     0.97033     0.97033     0.97033     1.03058
! 950.000     0.98544     0.98544     0.98544     1.01477
!1000.000     1.00000     1.00000     1.00000     1.00000
!+------------------------------------------------------------------
endmodule Functions_mod
