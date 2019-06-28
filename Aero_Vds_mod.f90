! <Aero_Vds_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
!=============================================================================
  module Aero_Vds_mod
!==============================================================================
  use PhysicalConstants_mod, only : FREEPATH, VISCO, BOLTZMANN, PI, GRAV, ROWATER
  !use Config_module,    only : DEBUG, MasterProc

  ! DESCRIPTION
  ! Calculates laminar sub-layer resistance (rb) and gravitational settling
  ! velocity (vs) for particles
  ! In DryDep_mod: Vd= Vs/(1.0 - exp ( -(Ra + Rb)*Vs ),   where
  ! Vs - gravitational settling velocity,
  ! Ra - aerodynamic resistance, Rb=1/Vds - viscous sub-layer resistance,
  ! Surface resistance is assumed zero for particles
  !---------------------------------------------------------------------------

 implicit none
  private

!  public  :: Aero_Vds
     ! A variety of equations are presented. Two types of
     ! stability corrections are used in the literature,
     ! those based upon 1/L only, and those based upon
     ! zi/L, where zi is the PBL height.
     !
     ! Using zi form gives max stab-fac of ca. 10
     ! (range max 2.8 to 10.7 for zi = 200 - 2500)
     ! Using 300 form gives max stab-fac of ca. 6

     ! For EMEP, use:
     ! WT versions : zi/L and limit 1/L to -0.04, cf Pryor

  public  :: SettlingVelocity     ! uses sigma as arg
  public  :: SettlingVelocityLn2  ! uses ln(sigma)**2 as arg
  public  :: PetroffFit
  public  :: GPF_VdsZi  ! General, Gallager-Petroff fits (loosely!)
  public  :: GPF_Vds300 !  "  "
  public  :: WeselyWT
  public  :: Wesely300
  public  :: Gallagher1997
  public  :: Gallagher2002
  public  :: GallagherWT
  public  :: Nemitz2004
  public  :: RuijgrokDrySO4
  public  :: RuijgrokWetSO4
  public  :: Wesely1985

contains

   !------------------------------------------------------------------------
   ! We have two SettlingVelocity fuctions, one using pre-calculated
   ! long(sigma)**2 and one using sigma.

   elemental function SettlingVelocityLn2(tsK,rho,lnsig2,Dg) result(Vs)
    ! gravitational settling for 2 modes
    ! Equations Axx referred to here are from Appendix A,
    !  Binkowski+Shankar, JGR, 1995
    !  Note confusing notation in B+S, dp was used for diff. coeff.
    !  here we use Di

     real, intent(in) :: tsK, rho   ! temp, air density
     real, intent(in) :: lnsig2, Dg ! geometric diameter 
     real  :: Vs

     real, parameter :: one2three = 1.0/3.0
     real    :: knut, Di,   & ! Knudsen number, Diffusion coefficient
                Di_help, vs_help
     !vind, vsmo, slip, stoke, schmidt ! slip correction, Stokes and Schmidt numbers

        knut = 2*FREEPATH/Dg   ! Knut's number

        Di_help =BOLTZMANN*tsK/(3*PI *VISCO *rho *Dg)       ! A30, Dpg
        vs_help= Dg*Dg * rho * GRAV / (18*VISCO*rho)        ! A32

       !... Diffusion coefficient for poly-disperse 
        Di = Di_help*(exp(-2.5*lnsig2)+1.246*knut*exp(-4.*lnsig2))      ! A29, Dpk 
       !... Settling velocity for poly-disperse 
        Vs = vs_help*(exp(8*lnsig2)+1.246*knut*exp(3.5*lnsig2)) ! A31, k=3

   end function SettlingVelocityLn2
   !------------------------------------------------------------------------
   elemental function SettlingVelocity(tsK,roa,sigma,diam,PMdens) result(Vs)
    ! gravitational settling velocity
    ! Equations Axx referred to here are from Appendix A,
    !  Binkowski+Shankar, JGR, 1995
    !  Note confusing notation in B+S, dp was used for diff. coeff.
    !  here we use Di

     real, intent(in) :: tsK, roa ! temp, air density, Rel.hum
     real, intent(in) :: sigma,diam,PMdens

     real             :: Vs ! (NSIZE)

     real    :: lnsig2, dg, & 
                knut, Di,   & ! Knudsen number, Diffusion coefficient
                Di_help, vs_help
 !-----------------------------------------------------------------------------------


        lnsig2 = log(sigma)**2

       !... mass median diameter -> geometric diameter

        dg = exp (log(diam) - 3.* lnsig2 )

        knut = 2.0*FREEPATH/dg   ! Knut's number

        Di_help =BOLTZMANN*tsK/(3*PI *VISCO *roa *dg)                    ! A30, Dpg
        vs_help= dg*dg * PMdens * GRAV / (18.0* VISCO*roa)        ! A32

       !... Diffusion coefficient for poly-disperse
        Di = Di_help*(exp(-2.5*lnsig2)+1.246*knut*exp(-4.*lnsig2)) ! A29, dpk
       !... Settling velocity for poly-disperse
        Vs = vs_help*(exp(8.0*lnsig2)+1.246*knut*exp(3.5*lnsig2))  ! A31, k=3

! Can't have output from elemental
!        if (DEBUG_VDS.and.MasterProc ) &
!            write(6,'(a,i3,es12.3,f10.3,5es12.3,3f9.2,f9.3)') &
!             "** Settling Vd ** ", roa, tsK, &
!             dg,knut,Di_help,vs_help,Di, lnsig2, &
!             1.0e6*diam, PMdens, sigma, Vs*100.0
!OLDER TEXT/QUERY
  ! Restrict settling velocity to 2cm/s. Seems
  ! very high otherwise,  e.g. see Fig. 4, Petroff et al., 2008 (Part I), where
  ! observed Vg for forests is usually < 2cm/s.


   end function SettlingVelocity

! -------------------------------------------------------------------
   function PetroffFit(ustar,invL,SAI) result(Vds)
     ! "Simple" fitting of Petroff et al. 2008
     ! Fig12b suggests  Vd = 0.3 cm/s for u* = 0.45, LAI=22
     ! ->  Vds = 0.007 * u*
     ! Fig.15 suggests that vds is approx prop to LAI for dp~0.4um
     ! We use SAI to keep some winter dep in decid forests
     ! To keep Vds above grassland, we use max(3.0, SAI)
     real, intent(in) :: ustar, invL,SAI
     real :: Vds

        Vds   = 0.007 * ustar * 0.1*max(SAI, 3.0)

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-300.0 * max(-0.04,invL))**0.6667)
       end if
   end function PetroffFit
   !------------------------------------------------------------------------
   ! GallagherPetrof fits
   ! Two functions here, for different stability methods
     ! "Simple" fitting of Gallagher et al. (1997) and Petroff et al., which 
     ! roughly captures the differences between Speulderbos-type and typical
     ! forests, because of LAI.
     ! Gallagher et al. had   Vds/u* = 0.0135 * Dp * stab function
     ! which gives 0.3 cm/s for neutral conditions, Dp=0.5
     !
     ! Fig12b suggests  Vd ~ 0.3 cm/s for u* = 0.45, LAI=22
     ! ->  Vds = 0.008 * u*
     ! Fig.15 suggests that vds is approx prop to LAI for dp~0.5um
     ! We use SAI to keep some winter dep in decid forests
     ! As Petroff started with a total LAI of 22, which is ca.
     ! 1-sided LAI=10, SAI=11, so we scale with SAI/11 = 0.09
     !
     ! We also limit the lowest Vds/u* to be 0.002, consistent with
     ! Wesely.

   function GPF_VdsZi(ustar,invL,SAI,zi) result(Vds)
     real, intent(in) :: ustar, invL,SAI, zi
     real :: Vds

        Vds   = ustar * max( 0.002, 0.008 * 0.1 * SAI )

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-0.3*zi*max(-0.04, invL))**0.6667)
       end if
   end function GPF_VdsZi
   !------------------------------------------------------------------------
   function GPF_Vds300(ustar,invL,SAI) result(Vds)
     real, intent(in) :: ustar, invL,SAI
     real :: Vds

        Vds   = ustar * max( 0.002, 0.008 * 0.1 * SAI )

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-300.0 * max(-0.04,invL))**0.6667)
       end if
   end function GPF_Vds300
   !------------------------------------------------------------------------
   function Wesely1985(ustar,invL, zi) result(Vds)
     real, intent(in) :: ustar, invL, zi
     real :: Vds

        Vds   = 0.002 * ustar  ! from grass

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-0.3*zi*max(-0.04, invL))**0.6667)
         !Alt: Vds = Vds * 0.0009 ( -zi*invL)**0.6667
       end if
   end function Wesely1985
   !------------------------------------------------------------------------
   function WeselyWT(ustar,invL, zi) result(Vds)
     ! Same as Wesely1985, but with 1/L limit
     real, intent(in) :: ustar, invL, zi
     real :: Vds

        Vds   = 0.002 * ustar  ! from grass

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-0.3*zi*max(-0.04,invL))**0.6667)
       end if
   end function WeselyWT
   !-- ----------------------------------------------------------------------
   function Wesely300(ustar,invL) result(Vds)
     ! Same as Wesely1985, but with 1/L limit and 300/L
     ! version of stability fac
     real, intent(in) :: ustar, invL
     real :: Vds

        Vds   = 0.002 * ustar  ! from grass

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-300.0 * max(-0.04,invL))**0.6667)
       end if
   end function Wesely300
   !-- ----------------------------------------------------------------------
   function RuijgrokDrySO4(ustar,RH,u_hveg) result(Vds)
     real, intent(in) :: u_hveg, RH, ustar ! RH in %
     real :: Vds

        Vds = 0.05*ustar**0.28 / u_hveg

        if ( RH > 80.0  ) then
           Vds = Vds * ( 1+ 0.18*exp( (RH-80)/20.0 ) )
        end if
   end function RuijgrokDrySO4
   !------------------------------------------------------------------------
   function RuijgrokWetSO4(u_hveg,RH,ustar) result(Vds)
     real, intent(in) :: u_hveg, RH, ustar ! RH in %
     real :: Vds

        Vds = 0.08*ustar**0.45 / u_hveg

        if ( RH > 80.0  ) then
           Vds = Vds * ( 1+ 0.37*exp( (RH-80)/20.0 ) )
        end if
   end function RuijgrokWetSO4
   !------------------------------------------------------------------------
   function Nemitz2004(dp,ustar,invL) result(Vds)
     real, intent(in) :: dp, ustar, invL
     real :: Vds

        Vds = 0.001*ustar

        if ( invL < 0.0 ) then
           Vds = Vds *( 1+( -(960*dp-88.0)*invL )**0.6667)
        end if
   end function Nemitz2004
   !------------------------------------------------------------------------
   function Gallagher1997(dp,ustar,invL) result(Vds)
     real, intent(in) :: dp, ustar, invL
     real :: Vds

        Vds = 0.0135 * ustar * dp

        if ( invL < 0.0 ) then
           Vds = Vds * (1.0+(-300*invL)**0.6667 )
        end if
   end function Gallagher1997
   !------------------------------------------------------------------------
   function Gallagher2002(ustar,invL,z0) result(Vds)
     real, intent(in) :: ustar, invL, z0
     real :: Vds
     real :: k1, k2

     !if( log(z0) > 0.0 ) then ! z0 > ~0.04 m
     !if( log10(z0) > 0.0 ) then ! z0 > ~0.04 m
       !k1 = 0.001222 * log(z0) + 0.003906
       k1 = 0.001222 * log10(z0) + 0.003906  ! Eqn (13)

     ! This equation has negative solutions. We set
     ! a small min value,  consistent wth 0.2 m/s from G02, Fig
     ! and u* = 0.5 m/s.
       k1 = min( 0.0004, k1)
     !else
     !  k1 = 0.001222
     !end if
     k2 = 0.0009

        if ( invL < 0.0 ) then
           Vds = ustar * &
          ( k1 + k2* (-300*invL)**0.6667 )
        else
           Vds = ustar * k1
        end if
   end function Gallagher2002
   !------------------------------------------------------------------------
   function GallagherWT(dp,ustar,invL, zi) result(Vds)
     ! Same as Gallagher1997, but with Wesely's zi form of
     ! stability correction, and 1/L limit
     real, intent(in) :: dp, ustar, invL,  zi
     real :: Vds

        Vds = 0.0135 * ustar * dp

        if ( invL < 0.0 ) then
          Vds = Vds * (1.0+(-0.3*zi*max(-0.04,invL))**0.6667)
        end if
   end function GallagherWT
   !-- ----------------------------------------------------------------------

! =================================================================

end module Aero_Vds_mod
