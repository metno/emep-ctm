! <PlumeRise_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
! This module contains different plume rise formula.


module PlumeRise_mod
  implicit none

  public :: Plume_ASME      ! Seems to work quite well, for large sources
  public :: Plume_HOllandS  ! Not sure if ok
  public :: Plume_HOllandQ  ! Not sure if ok
  public :: Plume_PreggerFriedrich   ! Seems to work quite well
  public :: Plume_NILU

  private :: penetr     ! NILU
  private :: building   ! NILU
  private :: uz_f       ! NILU
  private :: phi_m_i_f  ! NILU
  private :: phi_m_f    ! NILU

  logical, private, parameter :: debug = .false.
contains

 ! Plume rise functions for buoyant plumes, as coded in Seingfeld+Pandis

  function Plume_ASME( hs, F, u, Ta, dtdz ) result(dh)
     real, intent(in):: hs, F,u, Ta, dtdz  ! Ta in deg.C
     real :: dh

      if ( dtdz < 0.005 ) then
       dh = 7.4*(F*hs*hs)**0.333 / u  ! neutral and unstable
      else
       ! S1 = g. dt/dz / Ta, F/S1 = F*Ta/(g.dtdz)
       ! Seinfeld had 29. must have meant 2.9?
       ! 2.9 repdroduces Massimo's rs, see ASME.py
       !dh = 29* ( F*(Ta+273.15)/(9.81*dtdz*u ))**0.333 ! Ignore pressure term
       dh = 2.9* ( F*(Ta+273.15)/(9.81*dtdz*u ))**0.333 ! Ignore pressure term
      end if

  end function Plume_ASME

  function Plume_HollandS( hs,d, u, Vs, Ta, Ts, dtdz ) result(dh)
     real, intent(in):: hs, d,u, Vs, Ta, Ts, dtdz  ! Ta in deg.C
     real :: dh

      dh = Vs*d/u * ( 1.5 + 2.68e-3 * 1013.0 * (Ts-Ta)/(Ts+273.15) * d)

      if ( dtdz > 0.5 ) then
        dh = 0.8 * dh
      end if

  end function Plume_HollandS

  function Plume_HollandQ( hs, d, u, Vs, Ta, Ts, dtdz ) result(dh)
     real, intent(in):: hs, d,u, Vs, Ta, Ts, dtdz  ! Ta in deg.C
     real :: dh
     real :: QH
      QH = 41.868 ! * MW??
     !Pregger+Friedrich
     !MW = 1.36e-3 * Flow * (Ts-Ta)

      QH = 1.36e-3  * Vs * 3.142*d*d/4.0  * (Ts-Ta)
      QH = 41.868 * QH   ! MW -> cal/s

     ! Formulae from Thomas et al.,
      dh = (1.5 * Vs*d + 4.0e-5 * QH ) /u

      if ( dtdz > 0.5 ) then
        dh = 0.8 * dh
      end if

  end function Plume_HollandQ

  function Plume_PreggerFriedrich( M, u, Vs, d ) result(dh)
  !    Pregger, T. and R. Friedrich, Effective pollutant emission
  !      heights for atmospheric transport modelling based on 
  !      real-world information.
  !      Environmental Pollution, 157, 442-560, 2009
  !    input
  !    M   -  heat flux  (MW)
  !    u   -  wind speed at stack height (m/s)
  !    Vs  -  stack effluent velocity (m/s)
  !    d   -  stack diameter (m)  
     real, intent(in):: M,u , Vs, d ! Ta in deg.C
     real :: dh, dhm

      if ( M > 6.0 ) then 
         dh = 102 * M**0.6 / u
         !print *, "V.hot ", M, Vs,u, dh
      else if ( M > 1.4 ) then
         dh = 78.4 * M**0.75 / u
         !print *, "Q.hot ", M, Vs,u, dh
      else
         dh  = (0.35*Vs*d + 84*sqrt(M) ) / u
         dhm =  3.0 * Vs * d / u
         dh=max(dh, dhm)
         !print *, "Cold ", M, Vs,u, dh, dhm
      end if

  end function Plume_PreggerFriedrich
!  function PlumeHt( hs, d, Ts, Ta, Vs, F, u ) result(dh)
!     real, intent(in):: hs, d, Ts,Ta, Vs, F,u
!     real :: dh
!       d = max(1.0,d)
!
!     if ( Ts < 0 ) then
!          !h = hs
!     else ! MV Eqn 4-2, 
!         dh = 7.4*( 0.25 * 9.81 * d*d * Vs * ( Ts-Ta)*hs*hs/Ts)**0.3333 / u
!         dh = Vs*d/u * ( 1.5 + 2.68e-3 * 1013.0 * (Ts-Ta)/Ts * d)
!       ! modified from Pregger+Friedrich
!         dh =102 *  ( 1.36e-3*F*(Ts-Ta) )**0.6  ! PF Eqn 4
!
!         h = hs + dh
!
!     end if
!  end function PlumeHt
   
 
! ------------------ START NILU ------------------------------
 function plume_nilu(hs, ust, w,d,tg,ta,z2,hmix,v,cl,bh,bw) result (hnew)
    !----------------------------------------------------------------------
    !
    !****  calculates plume rise and returns
    !      injection height after penetration
    !
    !      plume-rise method (briggs,1969,1971,1975)'!
    !      nilu 30-5-92 trond boehler
    !      nilu 2012     sam-erik walker
    !      modified for use in unimod by matthias steffen karl (24.10.2012)
    !      copyright - nilu
    !
    !      purpose
    !      -------
    !
    !
    !      interface
    !      ---------
    !      input - see intent(in)
    !
    !
    !
    !       output
    !     
    !       hnew - final height after penetration (m)
    !
    !      method
    !      ------
    !
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      none
    !
    !------------------------------------------------------------------

    implicit none

      real, intent(in)        ::  &
          hs             &! Stack height (m)
         ,ust            &! Friction velocity
         ,w              &! stack emission velocity (m/s)
         ,d              &! Stack diameter (m)
         ,tg             &! stack gas temperature (k)
         ,ta             &! air temperature at stack height (k)
         ,z2             &! height of wind speed at 45m (m)
         ,v              &! wind speed at 45m (m/s)  (hard-coded???)
         ,hmix           &! mixing height (m)
         , cl             ! monin-obukhov length (m)

      real, intent(in), optional :: bh !  height of building
      real, intent(in), optional :: bw !  width

    ! output
      real    ::    hnew  ! FInal height after penetration

      real ::  & 
          xf      &! final distance
        , bl      &! minimum of building height and width
        , diffht  &! test height difference
        , dh      &! plume rise
        , dhb     &! plume rise due to buoyancy
        , dhb1    &! plume rise due to buoyancy
        , dhb2    &! plume rise due to buoyancy
        , dhm     &! plume rise due to momentum
        , dhp     &! plume rise interpolation height
        , f       &! plume rise variable
        , hf      &! plume height accounting for stack downwash
        , hp      &! plume rise height
        , qhmw    &! plume rise variable
        , rs      &! plume rise variable
        , s       &! plume rise variable
        , testhfl &! test height
        , tsd     &! test stack downwash
        , tz0     &! test surface height
        , zp       ! test height

! local variables
      real      :: testhp
      real      :: up,us,uz_min
      
      integer   :: id     &! calculated stack downwash indicator
                  ,idh    &! building cavity zone trap indicator
                 ,iwash   &! stack downwash indicator
                 ,niter    ! number of iterations

      real :: hfl    &! emission height (final value)
             ,sd     &! stack downwash height
             ,z0      ! surface roughness


! initialization start

       hfl   = 0.
       sd    = 0.
       iwash = 0

! surface roughness (m)
! Wieringa (1992) : "rough": 0.25
!    Wieringa, J., 1992, Updating the Davenport roughness classification, 
!       J. Wind Eng. Ind.Aerodynam., 41-44, pp 357-368.   
       z0    = 0.25
       id    = 0
       !ust   = 1.0
       uz_min  = 0.0   ! minimum wind speed (m/s)

! mongstad ccs test example
! building       bh  = 25.0
! building       bw  = 50.0
!       hs  = 60.0
       idh = 0
       !tg  = 313.0
       !d   = 7.14
       !w   = 10.0
       !hmix= 150.0
       !ta  = 273.0

! initialization end

! minimum for testing building influence
     if ( present(bh) .and. present(bw) ) then
        bl = min(bh,bw)
     else
        bl = 0.0
     end if
     

! for emep the 45m wind speed is used as wind speed at stack

      us = v

! test parameter for stack induced downwash

      tsd = w/us

! momentum rise for all stability classes

      dhm = 3.*d*w/us

! interpolation between momentum rise and stack induced downwash

      if (tsd  >  2.) then

! full momentum

          sd = 0.
          hf = hs + dhm
          if(debug) print *, "tsd>2", hf

      else if (tsd  <  1.) then

! full stack downwash

          iwash = 1
          sd = 2.*(w/us - 1.5)*d
          if (id  ==  1) sd = 0.
          hf = hs + sd
          if(debug) print *, "tsd<1", hf

      else

! partial momentum and downwash, perform interpolation

          dhp = (7.*tsd - 8.)*d
          if (dhp  <  0.) then
              iwash = 1
              if (id  ==  1) dhp = 0.
          end if
          hf = hs + dhp
          if(debug) print *, "tsd: ", hf

      end if

! no downwash below ground

      if (hf  <  0) hf = 0.

! evaluation of building influence

     if ( present(bh) .and. present(bw) ) then
      call building(bh,bw,hs,v,hf,hp,dhm,idh,iwash)
     else   !no building effects
      idh = 1
      hp = hs
     end if

! define new test height

      testhp = hp

      niter = 0

! do loop replaces goto      
      do 100 niter=1,100

! wind speed at plume height
! later could be replaced by emep wind speed at respective model layer

      up  = uz_f(testhp,uz_min,z2,v,ust,1/cl)
      
          if(debug) print *, "testhp: ", niter,  testhp, up

      if (tg  >  ta) then

! bouyant plumes

          f = 9.81*w*(d/2.)**2*(tg - ta)/tg
!          qhmw = 0.11*f
          if(debug) print *, "testhp: ", niter,  testhp, up, tg, ta, f, cl

! plume-rise method (briggs,1969,1971,1975)

          if (cl  <  0. .or. cl  >  200.) then

! unstable and neutral condition

! bouyancy rise

              if (f  <  55.) then
                  xf  = 0.049*f**(5./8.)
                  dhb = 21.425*f**(3./4.)/up
              else
                  xf  = 0.119*f**(2./5.)
                  dhb = 38.71*f**(3./5.)/up
              end if
          if(debug) print *, "neutrl:  ", niter,f,  xf, up, dhb

          else

! stable condition

              rs = 0.02
              if (cl .le. 40.) rs = 0.035

! calculate stability parameter

              s = 9.81*rs/ta

! distance to final rise

              xf = 0.0020715*up*s**(-1./2.)

! buoyancy rise, calm conditions

              dhb1 = 4.*f**(1./4.)*s**(-3./8.)

! buoyancy rise, wind

              dhb2 = 2.6*(f/(up*s))**(1./3.)

! choose the lower of the two values

              dhb = min(dhb1,dhb2)
          if(debug) print *, "stab:  ", niter, rs, s, xf, up, dhb1, dhb2

          end if

      else

! cold releases

          dhb = 0.

          if (cl  <  0. .or. cl  >  200.) then

! momentum rise, stable for cold releases

              rs = 0.02

! calculate stability parameter

              s = 9.81*rs/ta

! calculate momentum rise

              dhm = 1.5*((w*w*d*d*ta)/(4.*tg*up))**(1./3.)*s**(-1./6.)

          end if

      end if

! no momentum rise if stack downwash

      if (iwash  ==  1) dhm = 0.

! choose the highest value of momentum and buoyancy rise

      dh = max(dhb,dhm)

! if momentum rise then zero distance to final height

      if (dhm  >  dhb) xf = 0.

! final distance in meters

      xf = 1000.*xf

! final plume rise

      if (idh  ==  3) then

! trapped in the cavity-zone

          hfl = 0.5*bl
          xf  = 0.
          exit

      else

          hfl = hp + dh

! iteration on plume rise
          diffht  = abs(testhp - hfl)
          testhfl = 0.05*hfl
          testhp  = hfl
          if (diffht <= testhfl) exit
      
      end if

!      if (idh .ne. 3 .and. diffht  >  testhfl .and.   &
!         niter .le. 100) goto 100

  100  continue    

! calculate influence of penetration

      call penetr(hs,hmix,hfl,idh,hnew)

! minimum height for windspeed calculation

      tz0 = 10.*z0
      if (hnew .le. tz0) hnew = 10.0


  end function plume_nilu



  subroutine penetr(hs,hmix,hfl,idh,hnew)

!  ----------------------------------------------------------------
! the subroutine calculates plume penetration
!
!  input
! hs   - height of stack
! hmix - mixing height
! hfl  - final plume rise height
! idh  - stack downwash indicator
!
!  output
! hnew - new final plume rise height
!
!  local
! dm   - difference value 
! dp   - difference value
! hpen - penetration height
! rmp  - ratio value
! ps   - degree of penetration
!
!  ----------------------------------------------------------------
   implicit none

         integer, intent(in)  :: idh
         real, intent(in)     :: hs,hmix,hfl
         real, intent(out)    :: hnew

! local variables
         real                 :: dm,dp,hpen,rmp,ps 



      hnew = 10.0     

! if the stack is above the mixing height then total penetration

      if (hs  >  hmix) then

! plume height set equal to the stack height

          ps = 1.0
          hnew = hs
          return

      end if

      if (idh  ==  3) then

! if trapped in the cavity sone then no penetration

          ps = 0.0
          hnew = hfl
          return

      end if

      dm = hmix - hs
      dp = hfl  - hs

      if (dp <= 0.0) then

! downwash or fixed plume height

          rmp = hmix/hfl

      else

          rmp = dm/dp

      end if
!print *, "rmp  ", hmix, hfl, rmp  

      if (rmp >= 1.5) then

! no penetration

          ps   = 0.0
          hpen = hs + 0.62*dm
!print *, "no      p", hnew, hpen, dm

      elseif (rmp <= 0.5) then

! total penetration

          ps   = 1.0
          hpen = hs + dm
!print *, "total   p", hnew, hpen

      else

! partial penetration

          ps   = 1.5 - rmp
          hpen = hs + (0.62 + 0.38*ps)*dm
!print *, "partial p", hnew, hpen

      end if

! choose the lowest of effective plume rise and rise due to penetration

      hnew = hfl
      if (hpen  <  hnew) hnew = hpen
!print *, "Subr penetr: hpen hnew", hpen, hnew
!stop 

  end subroutine penetr



  subroutine building(bh,bw,hs,v,hf,hp,dhm,idh,iwash)

!  ----------------------------------------------------------------
! the subroutine calculates modified stack height due to the
! influence of nearby buildings.
!
!  input
! bh    - building height
! bw    - building width
! hs    - height of stack
! v     - windspeed
! hf    - plume height accounting for stack downwash
! dhm   - d <  height due to momentum rise
! iwash - stack downwash indicator
!
!  output
! hp    - modified stack height
! idh   - building cavity zone trap indicator
!
!  local
! bl - minimum of building height and width
! bt - building influence test variable
! he - height value
! th - building influence test variable
!
!  ----------------------------------------------------------------
   implicit none

         integer, intent(in)   :: iwash
         real, intent(in)      :: bh,bw,hs,v,hf,dhm
         real, intent(out)     :: hp
         integer, intent(out)  :: idh


! local variables
         real                  :: bl,bt,he,th


      bl = min(bh,bw)

! no building effects for v < 1.5 m/s

      if (v  <  1.5) then

          idh = 1
          hp  = hs

! if stack downwash then stack height is modified

          if (iwash  ==  1) hp = hf
          return

      end if

      bt = bh + 1.5*bl

      if (bt  ==  0. .or. hf  >  bt) then

! no building effects

          idh = 1
          hp  = hs

! if stack downwash then stack height is modified

          if (iwash  ==  1) hp = hf

      else

! building effects

          if (hf  >  bh) then
              he = 2.*hf - (bh + 1.5*bl)
          else
              he = hf - 1.5*bl
          end if

          th = 0.5*bl

          if (he  >  th) then

              idh = 2

              if (iwash  ==  0) then

                  hp = he - dhm

              else

! stack downwash

                  hp = he

              end if

          else

! trapped in the cavity zone, ground area source = hb*bb

              idh = 3
              hp  = 0.0

          end if

      end if

      return

  end subroutine building


  real function uz_f(z,uz_min,zref,uref,ustar,lmo_inv)
!  ----------------------------------------------------------------  
! The function calculates the wind speed at height z above ground
!
! z       - Actual height above ground (m)
! uz_min  - Minimum allowed wind speed value (m/s)
! zref    - Reference height above ground (m)
! uref    - Wind speed at reference height (m/s)
! ustar   - Monin-Obukhov friction velocity (m/s)
! lmo_inv - Inverse Monin-Obukhov length (1/m)
!
!      nilu 2012 - sam-erik walker
!      copyright nilu
!
!  ----------------------------------------------------------------

   implicit none
   
! Arguments
      real, intent(in) :: z
      real, intent(in) :: uz_min
      real, intent(in) :: zref
      real, intent(in) :: uref
      real, intent(in) :: ustar
      real, intent(in) :: lmo_inv
      real, parameter  :: kappa = 0.4


! Calculate the wind speed at height z above ground

      uz_f = uref + (ustar/kappa)*phi_m_i_f(zref,z,lmo_inv)

! Ensure wind speed is not below minimum allowed value

      uz_f = max(uz_f,uz_min)

  end function uz_f
  

  real function phi_m_i_f(za,zb,lmo_inv)
!  ----------------------------------------------------------------  
!
! The subroutine calculates the integral of the function phi_m between 
! the two heights za and zb
!
!      nilu 2012 - sam-erik walker
!      copyright nilu
!
!  ----------------------------------------------------------------
   implicit none
   
! Arguments
      real, intent(in) :: za
      real, intent(in) :: zb
      real, intent(in) :: lmo_inv

! za      - Height above ground (m)
! zb      - Height above ground (m)
! lmo_inv - Inverse Monin-Obukhov length (1/m)

      real, parameter :: beta_m = 5.3

! Local variables

      real :: am
      real :: bm
      real :: cm
      real :: dm
      real :: lmo
      real :: xl
      real :: xu
      real :: zl
      real :: zu

! xl - Phi_m function value at zl
! xu - Phi_m function value at zu
! zl - Lower height above ground (m)
! zu - Upper height above ground (m)

! Initialize
      phi_m_i_f = 0.0

! Calculate ordered heights

      if (za <= zb) then
        zl = za
        zu = zb
      else
        zl = zb
        zu = za
      end if

      if (lmo_inv < 0) then

! Unstable case

        xl = phi_m_f(zl,lmo_inv)
        xu = phi_m_f(zu,lmo_inv)

        phi_m_i_f = log(zu/zl) - 2.0*log((1.0 + xu)/(1.0 + xl)) -     &
                  log((1.0 + xu**2)/(1.0 + xl**2)) +                  &
                  2.0*(atan(xu) - atan(xl)) 

      else if (lmo_inv > 0) then

! Stable case

        am = 1.0
        bm = beta_m
        cm = beta_m
        dm = 1.0

        lmo = 1.0/lmo_inv

        if (zl <= zu .and. zu <= lmo) then
          phi_m_i_f = am*log(zu/zl) + bm*((zu - zl)/lmo)
        else if (zl < lmo .and. lmo < zu) then
          phi_m_i_f = am*log(lmo/zl) + bm*(1.0 - zl/lmo) +    &
                     cm*log(zu/lmo) + dm*(zu/lmo - 1.0)
        else if (lmo <= zl .and. zl <= zu) then
          phi_m_i_f = cm*log(zu/zl) + dm*((zu - zl)/lmo)
        end if

!       xl = phi_m_f(zl,lmo_inv)
!       xu = phi_m_f(zu,lmo_inv)

!       phi_m_i_f = log(zu/zl) + (xu - xl)

      else

! Neutral case

        phi_m_i_f = log(zu/zl)

      end if

! Change sign of integral if za > zb

      if (za > zb) phi_m_i_f = -phi_m_i_f

      end function phi_m_i_f


      real function phi_m_f(z,lmo_inv)
!  ----------------------------------------------------------------  
!
! The subroutine calculates the Monin-Obukhov universal stability 
! function for wind speed at height z above ground
!
!      nilu 2012 - sam-erik walker
!      copyright nilu
!
!  ----------------------------------------------------------------
   implicit none
   
! Arguments

      real, intent(in) :: z
      real, intent(in) :: lmo_inv

! z       - Height above ground (m)
! lmo_inv - Inverse Monin-Obukhov length (1/m)

      real, parameter :: alfa_m = -19.0
      real, parameter :: beta_m = 5.3
      real, parameter :: mete_miss = -9900.0

! Local variables

      real :: zeta

! zeta - Stability parameter

! Negative height is not allowed

      if (z < 0.0) then
        phi_m_f = mete_miss
        return
      end if

      if (lmo_inv < 0) then

! Unstable case

        zeta = z*lmo_inv

        if (-2.0 <= zeta .and. zeta <= 0.0) then

! Use formulation from (Hoegstroem, 1996) (in validity range)

          phi_m_f = (1.0 + alfa_m*zeta)**(-0.25)

        else

! Use formulation from (Hoegstroem, 1996) (extended)

          phi_m_f = (1.0 + alfa_m*zeta)**(-0.25)

        end if

      else if (lmo_inv > 0) then

! Stable case

        zeta = z*lmo_inv

        if (0.0 <= zeta .and. zeta <= 0.5) then

! Use formulation from (Hoegstroem, 1996) (in validity range)

          phi_m_f = 1.0 + beta_m*zeta

        else if (0.0 <= zeta .and. zeta <= 1.0) then

! Use formulation from (Hoegstroem, 1996) (extended)

          phi_m_f = 1.0 + beta_m*zeta

        else

! Use formulation from (Holtslag, 1990) (very stable conditions)

          phi_m_f = beta_m + zeta

        end if

      else

! Neutral case

        phi_m_f = 1.0

      end if

    end function phi_m_f


! ------------------ end nilu ------------------------------

 
end module PlumeRise_mod

