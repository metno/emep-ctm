! <MARS_Aero_water_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
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
      module MARS_Aero_water_mod

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! DESCRIPTION
! 
! Purpose   : Calculates aerosols' liquid water content
!
! Subroutine: Awater
!             Input:  - relative humidity: relh
!                     - number of micromoles/(m^3 of air) for sulfate,
!                       ammonium, and nitrate: mso4, mnh4, mno3
!             Output: - water amount in micrograms/(m^3 of air): wh2o
!
! Author    : Dr. Francis S. Binkowski, 4/8/96
!           : modified slightly for EMEP model
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit none
      private

! subroutines:
      public :: Awater

! functions
      private :: poly4, poly6


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! coefficients for polynomials (function poly4) to be defined
! at start of routine:

! Define saved variables:
      real, private,parameter, dimension(4) :: & !(for x = 0, 1, 1.5 and 2)
        C0 =  (/ 0.798079,  -1.574367,  2.536686,   -1.735297 /),&
        C1 =  (/ 0.9995178, -0.7952896, 0.99683673, -1.143874 /),&
        C15=  (/ 1.697092,  -4.045936,  5.833688, -3.463783 /),&
        C2 =  (/ 2.085067,  -6.024139,  8.967967, -5.002934 /)

      real, private,parameter, dimension(6) ::  &
        KNO3 = (/ 0.2906,    6.83665, -26.9093,  46.6983, -38.803, 11.8837/),&
        KSO4 = (/ 2.27515, -11.147,    36.3369, -64.2134, 56.8341 ,-20.0953/)

     ! Set molecular weights:

      real, private, parameter :: &
        MWSO4   = 96.0 &
       ,MWNH4   = 18.0 &
       ,MWNO3   = 62.0 &
       ,MW2     = MWSO4 + 2.0 * MWNH4  &
       ,MWANO3  = MWNO3 + MWNH4


      contains

      subroutine Awater(relh,mso4,mnh4,mno3,wh2o)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! This routine uses polynomials rather than tables, and uses empirical
! polynomials for the mass fraction of solute (mfs) as a function of
! water activity where:
!
!          mfs = ms / ( ms + mw)
!          ms - the mass of solute
!          mw - the mass of water
!
! Define   y = mw / ms
!
! then     mfs = 1 / (1 + y)
!
!     y can then be obtained from the values of mfs as
!
!          y = (1 - mfs) / mfs
!
! The aerosol is assumed to be in a metastable state if the relative
! humidity (rh) is is below the rh of deliquescence, but above the
! rh of crystallization.
!
! The Zdanovskii-Stokes-Robinson relation ('ZSR') is used for sulfates
! with x (molar ratio of ammonium to sulfate) in the range 0 <= x <= 2,
! subdivided into four sections:
!
!          section 1: 0 <= x < 1
!          section 2: 1 <= x < 1.5
!          section 3: 1.5 <= x < 2.0
!          section 4: 2 <= x
!
! In sections 1 through 3 only the sulfates can affect the amount
! of water on the particles. In section 4 we have fully neutralized
! sulfate, and extra ammonium which allows more nitrate to be present.
! Thus, the ammount of water is calculated using ZSR for ammonium
! sulfate and ammonium nitrate. Crystallization is assumed to occur
! in sections 2, 3, and 4 (see detailed discussion below).
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Description of variables used in the routine:
!
!
! mso4, mnh4, and mno3 : number of micromoles/(cubic meter of air) for
!                        sulfate, ammonium, and nitrate, respectively
! relh : relative humidity (%)
! wh2o : returned water amount in micrograms /(cubic meter of air)
! x    : molar ratio of ammonium to sulfate
! y0, y1, y15, y2 : water contents in mass of water/mass of solute
!                   for pure aqueous solutions with x equal to 0, 1,
!                   1.5, and 2, respectively
! y3 : the value of the mass ratio of water to solute for a pure
!      ammonium nitrate solution.
!
!
! C1, C15, C2: (defined at start of routine!)
! The polynomials use data for relh as a function of mfs from Tang
! and Munkelwitz, JGR. 99: 18801-18808, 1994.
! The polynomials were fit to Tang's values of water activity as a
! function of mfs.
! C1, C15, C2 are the coefficients of polynomials fit to Tang and
! Munkelwitz data giving mfs as a function of water activity.
!
! C0: (defined at start of routine!)
! fit to data from
! Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
! Giaque et al. J. Am. Chem. Soc., 82: 62-70, 1960
! Zeleznik J. Phys. Chem. ref. data, 20: 157-1200
!
! KNO3, KSO4: (defined at start of routine!)
! The polynomials for ammonium nitrate and ammonium sulfate are from:
! Chan et al.1992, atmospheric environment (26a): 1661-1673.
!
!
! tso4, tnh4, tno3: mole concentrations used in the calculations
!   ( tso4 = max(mso4,0.),  tnh4 = max(mnh4,0.),  tno3 = max(mno3,0.) )
!
! aw  : relative humidity used in the calculations (.01<aw<.95 required)
! awc : relative humidity along the crystallization curve
! u   : help variable for rh used for interpolation (u=40%)
! x   : molar ratio (x = tnh4 / tso4 for non-zero sulfate, and =10 else)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!-- in          
      real, intent(in)  :: relh, mso4, mnh4, mno3

!-- out
      real, intent(out) :: wh2o


!-- local
      real :: tso4, tnh4, tno3
      real :: x, awc, aw, u                          &
            , mfs0, mfs1, mfs15, mfs2                &
            , mfsso4, mfsno3                         &
            , y, y0, y1, y15, y2, y3                 &
            , y40, y140, y1540, yc

! Check range of per cent relative humidity:
! aw - water activity = fractional relative humidity
      aw = relh
! Set aw between 0.01 and 0.95:
      aw = max(0.01,aw)
      aw = min(aw,0.95)


      tso4 = max(mso4,0.)
      tnh4 = max(mnh4,0.)
      tno3 = max(mno3,0.)

! Initialise molar ratio:
      x = 0.

! If there is non-zero sulfate calculate the molar ratio:
      if (tso4 > 0.) then
        x = tnh4 / tso4
      else
! ... otherwise check for non-zero nitrate and ammonium
        if (tno3 > 0. .and. tnh4 > 0.) x = 10.
      end if

! Begin screen on x for calculating wh2o:
      if ( x < 1. ) then

        mfs0 = poly4(C0,aw)
        mfs1 = poly4(C1,aw)
        y0 = ( 1. - mfs0 ) / mfs0
        y1 = ( 1. - mfs1 ) / mfs1
        y  = ( 1. - x ) * y0 + x * y1

      else if ( x < 1.5) then

        if ( aw >= 0.40 ) then

          mfs1  = poly4(C1,aw)
          mfs15 = poly4(C15,aw)
          y1    = (1. - mfs1 ) / mfs1
          y15   = (1. - mfs15) / mfs15
          y     = 2. * ( y1 * (1.5 - x) + y15 *( x - 1.) )

        else

! Setup for crystalization:

! Crystallization is done as follows:
!
!   for 1.5 <= x : crystallization is assumed to occur at rh = 0.4
!   for x <= 1.0 : crystallization is assumed to occur at an rh < 0.01
!
!   and since the code does not allow rh < 0.01, crystallization is
!   assumed not to occur in this range.
!
!   for 1.0 <= x <= 1.5 : the crystallization curve is a straight line
!                         from a value of y15 at rh = 0.4 to a value of
!                         zero at y1. From point b to point a in the
!                         diagram the algorithm does a double inter-
!                         polation to calculate the amount of water.
!
!               y1(0.40)               y15(0.40)
!                +                     + point b
!
!
!
!                +---------------------+
!                x=1                   x=1.5
!              point a
!

          awc = 0.80 * (x - 1.0) ! rh along the crystallization curve.

          y = 0.0
          u=0.40

          if ( aw >= awc ) then

! Interpolate using crystalization curve:

            mfs1  = poly4(C1,u)
            mfs15 = poly4(C15,u)
            y140  = (1.0 - mfs1 ) / mfs1
            y1540 = (1.0 - mfs15) / mfs15
            y40 = 2.0 * ( y140 * (1.5 - x) + y1540 *( x - 1.0) )
            yc = 2.0 * y1540 * (x -1.0) ! y along crystallization curve
            y = y40 - (y40 - yc) * (u - aw) / (u - awc)

          end if ! end of "if ( aw >= awc ) then"

        end if ! end of "if ( aw >= 0.40 ) then"

      else if ( x < 1.9999) then

         y= 0.0

         if (aw >= 0.40) then
           mfs15 = poly4(C15,aw)
           mfs2  = poly4(C2,aw)
           y15   = (1.0 - mfs15) / mfs15
           y2    = (1.0 - mfs2) / mfs2
           y     = 2.0 * (y15 * (2.0 - x) + y2 * (x - 1.5) )
         end if ! end of check for crystallization

      else

! i.e. x >= 1.9999
!
! Regime where ammonium sulfate and ammonium nitrate are in solution!
!
! Following cf&s for both ammonium sulfate and ammonium nitrate
! Check for crystallization here. Their data indicate a 40% value
! is appropriate.

        y2 = 0.0
        y3 = 0.0

        if (aw >= 0.40) then
          mfsso4 = poly6(KSO4,aw)
          mfsno3 = poly6(KNO3,aw)
          y2     = (1.0 - mfsso4) / mfsso4
          y3     = (1.0 - mfsno3) / mfsno3
        end if

      end if ! end of 'if ( x < 1. ) then'

! Now set up output of wh2o
!
! wh2o units are micrograms(liquid water) / (cubic meter of air)

      if ( x < 1.9999) then

        wh2o =  y * (tso4 * MWSO4 + tnh4 * MWNH4 )

      else

! This is the case when all the sulfate is ammonium sulfate
! and the excess ammonium forms ammonum nitrate

        wh2o =   y2 * tso4 * MW2 + y3 * tno3 * MWANO3

      end if


      end subroutine Awater

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real function poly4 (a,x)

! Calculates the polynomial based on 4 coefficients a(1:4):

!-- arguments
      real, dimension(4), intent(in)  ::  a
      real, intent(in)  ::  x

      poly4 = a(1) + x * ( a(2) + x * ( a(3) + x * ( a(4) ) ) )
      
      end function poly4

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real function poly6(a,x)

! Calculates the polynomial based on 6 coefficients a(1:6):

!-- arguments
      real, dimension(6), intent(in)  ::  a
      real, intent(in)  ::  x

      poly6 = a(1) + x * ( a(2) + x * ( a(3) + x * ( a(4) +  &
              x * ( a(5) + x * (a(6)  ) ) ) ) )

      end  function poly6 

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end module MARS_Aero_water_mod

