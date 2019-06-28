! <AeroFunctions.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
! <AeroFunctions.f90 - A component of the EMEP MSC-W Chemical transport Model>
! A collection of functions as used in the EMEP model or in some cases
! just for testing.
!*****************************************************************************! 
module AeroFunctions_mod
!____________________________________________________________________
!
 use PhysicalConstants_mod, only : AVOG, RGAS_J, PI
  implicit none
  private

 ! From VsLogNormalCalcs.
  public :: DpgV2DpgN ! from mass/volume Dp to number Dp
  public :: DpgN2DpgV  ! from number Dp to mass/volume Dp
  public :: LogNormFracBelow  ! Frac Mass below Dp

  public :: cMolSpeed

  public ::  kaero !aerosol production rate
  public :: SurfArea_Mono
  public :: SurfArea_Poly
  public :: pmSurfArea

  ! Gerber functions used to get wet radius
  public :: pmH2O_gerber
  public :: WetRad
  public :: umWetRad
  public :: cmWetRad
  public :: WetRadS

  integer, public, parameter :: GbRural=1, GbSeaSalt=2, GbUrban=3, GbAmmSO4=4

  public :: UptakeRate

  public :: GammaN2O5
  public :: GammaN2O5_so4
  public :: GammaN2O5_om  
  !public :: GammaN2O5_ome !testing, v. low value for OM
  public :: GammaN2O5_EJSS
  public :: GammaN2O5_DavisAN

  public :: self_test
  public :: self_test_fracs



  !========================================
  contains
  !========================================

  !---------------------------------------------------------------------
  !> FUNCTION DpgV2DpgN ! from mass/volume Dp median to number median Dp
  ! NB:  Senfeld & Pandis use \overline{D_{pgV}} for volume median and
  ! \overline{D_{pg}} for number median

  function DpgV2DpgN(DpgV,sigma_g) result (DpgN)
      real, intent(in) :: DpgV,sigma_g
      real :: DpgN
      DpgN = exp( log(DpgV) - 3 * ( log(sigma_g) )**2 ) ! mass to N, S7.52
   end function DpgV2DpgN

  !---------------------------------------------------------------------
  !> FUNCTION DpgN2DpgV  ! from number Dp to mass/volume median Dp

   function DpgN2DpgV(DpgN,sigma_g) result (DpgV)
      real, intent(in) :: DpgN,sigma_g
      real :: DpgV
      DpgV = exp( log(DpgN) + 3 * ( log(sigma_g) )**2 ) ! N to mass
   end function DpgN2DpgV

  !---------------------------------------------------------------------
  !> FUNCTION LogNormFracBelow
  !! Calculates fraction of aerosol mass below given diameter Dp, for a given
  !! MMD of Dpg and sigma sig. Uses eqn 7.46 for Fn from Seinfeld + Pandis

 function LogNormFracBelow(Dpg,sig,Dp) result(Fn)
   real, intent(in) :: Dpg, sig, Dp
   real :: Fn
   real, parameter :: sqrt_two = sqrt(2.0)

      Fn =  0.5 + 0.5 * erf( (log(Dp)-log(Dpg)) / ( sqrt_two * log(sig) ))

  end function LogNormFracBelow

  !---------------------------------------------------------------------
  !> FUNCTION cMolSpeed    molecular speed calcs

  elemental function cMolSpeed(tK,mw) result(c) 
     real, intent(in) :: tK            !< Temp, K
     real, intent(in) :: mw            !< mol.wt, g/mole
     real             :: c             ! mean mol. speed, m/s
     ! 3 RT/(0.001 mw) 

     c = sqrt(3.0e3 * RGAS_J * tK / mw)

  end function cMolSpeed

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
!  elemental

 function WetRadS(rdry,fRH,pmtype) result (rwet)
   real, intent(in) :: rdry                 !< UNITS, any, but same as rwet
   real, intent(in) :: fRH                  !< humidity, 0< fRH < 1
   integer, intent(in), optional :: pmtype  !< Type of aerosol
   real             :: rwet                 !< 

   ! Gerber, TAble 2, simplified for ext,. coeffs.
   real, parameter, dimension(4) :: &
     !          R    U     SS2    SS3    ! NAM types, TAble 2, Gerber
      C7 =  (/ 1.17, 1.28, 1.83, 1.97 /) &
     ,C8 =  (/ 1.87, 2.41, 5.13, 5.83 /)
     ! BUG C7 =  (/ 1.17, 1.28, 1.97, 1.83 /) &
     ! BUG ,C8 =  (/ 1.87, 2.41, 5.83, 5.13 /)
   real, parameter :: THIRD = 1.0/3.0
   real :: f
   integer :: ind  

   ind = 1 ! default = rural
   if ( present(pmtype) ) ind = pmtype

   f = ( (C7(ind)-fRH)/(C8(ind)*(1-fRH)) )**THIRD
   rwet = rdry * f
!   print *, "gerber simplified ", f, rdry, rwet
  end function wetradS 

  !---------------------------------------------------------------------

  elemental function WetRad(rdry,fRH,pmtype) result (rwet)
   real, intent(in) :: rdry                 !< m !NOTE UNITS DO MATTER!
   real, intent(in) :: fRH                  !< humidity, 0< fRH < 1
   integer, intent(in), optional :: pmtype  !< Type of aerosol
   real             :: rwet                 !< m

   ! Zhang, 2001,  Table 1, constants for Gerber, Gong typr calc
   real, parameter, dimension(4,4) :: K = reshape(  &
      (/  0.2789,3.115,5.415e-11,-1.399   &  ! rural, default
        , 0.7674,3.079,2.573e-11,-1.424   &  ! sea-salt
        , 0.3926,3.101,4.190e-11,-1.404   &  ! urban
        , 0.4809,3.082,3.110e-11,-1.428 /)&  ! (NH4)2SO4
        ,(/4,4/) )
   real, parameter :: THIRD = 1.0/3.0, cm2m = 1.0e-2
   real :: rd, mrh
   integer :: ind

   ind = 1 ! default = rural
   if ( present(pmtype) ) ind = pmtype

   mrh = max(0.01,fRH)
   rd = rdry * 1.0e2  ! m -> cm
   rwet = cm2m * ( K(1,ind)*rd**K(2,ind) / &
     (K(3,ind) *rd**K(4,ind) - log10(mrh))+rd**3) ** THIRD
  end function WetRad 
  !---------------------------------------------------------------------

  elemental function umWetRad(rdry,fRH,pmtype) result (rwet)
   real, intent(in) :: rdry                 !< um !NOTE UNITS DO MATTER!
   real, intent(in) :: fRH                  !< humidity, 0< fRH < 1
   integer, intent(in), optional :: pmtype  !< Type of aerosol
   real             :: rwet                 !< um

   ! Zhang, 2001,  Table 1, constants for Gerber, Gong typr calc
   real, parameter, dimension(4,4) :: K = reshape(  &
      (/  0.2789,3.115,5.415e-11,-1.399   &  ! rural, default
        , 0.7674,3.079,2.573e-11,-1.424   &  ! sea-salt
        , 0.3926,3.101,4.190e-11,-1.404   &  ! urban
        , 0.4809,3.082,3.110e-11,-1.428 /)&  ! (NH4)2SO4
        ,(/4,4/) )
   real, parameter :: THIRD = 1.0/3.0
   integer :: ind  

   ind = 1 ! default = rural
   if ( present(pmtype) ) ind = pmtype

   rwet = ( K(1,ind)*rdry**K(2,ind) / &
     (K(3,ind) *rdry**K(4,ind) - log10(fRH))+rdry**3) ** THIRD

   ! to stop huge wet particles...
   rwet = min( rwet, 10*rdry)
  end function umWetRad 
 !---------------------------------------------------------------------

  elemental function cmWetRad(rdry,fRH,pmtype) result (rwet)
   real, intent(in) :: rdry                 !< m !NOTE UNITS DO MATTER!
   real, intent(in) :: fRH                  !< humidity, 0< fRH < 1
   integer, intent(in), optional :: pmtype  !< Type of aerosol
   real             :: rwet                 !< m

   ! Zhang, 2001,  Table 1, constants for Gerber, Gong typr calc
   real, parameter, dimension(4,4) :: K = reshape(  &
      (/  0.2789,3.115,5.415e-11,-1.399   &  ! rural, default
        , 0.7674,3.079,2.573e-11,-1.424   &  ! sea-salt
        , 0.3926,3.101,4.190e-11,-1.404   &  ! urban
        , 0.4809,3.082,3.110e-11,-1.428 /)&  ! (NH4)2SO4
        ,(/4,4/) )
   real, parameter :: THIRD = 1.0/3.0, um2cm = 1.0e-4, cm2um = 1.0e4
   real :: rd, mrh
   integer :: ind  

   ind = 1 ! default = rural
   if ( present(pmtype) ) ind = pmtype

   mrh = max(0.01,fRH)
   rd = rdry * um2cm ! um -> cm
   rwet = cm2um * ( K(1,ind)*rd**K(2,ind) / &
     (K(3,ind) *rd**K(4,ind) - log10(mrh))+rd**3) ** THIRD
  end function cmWetRad 
  !---------------------------------------------------------------------

!OLD:
! N2O5 -> nitrate calculation
!===========================================================================
! N2O5 -> nitrate calculation. Some constants for
! calculation of volume fraction of sulphate aerosol, and rate of uptake
! 
!
! The first order reaction coefficient K (corrected for gas phase diffusion, 
! Schwartz, 1986) is given by
!
! K= A* alpha* v/4
!    alpha=sticking coeff. for N2O5 =0.02
!    v=mean molecular speed for N2O5
!    A=aerosol surfac
!
! The surface area of the aerosols can be calculated as
! 
! A = V * surface/volume of aerosols
!     V=volume fraction of sulphate (cm3 aerosol/cm3 air)
!     (similar for nitrate and ammonium):
!
! The surface/volume ratio is calculated using Whitby particle distribution
! with number mean radius r=0.5*DpgN  and standars deviation (Sigma, eg 2). 
!OLD Then surface/volume=3/r *  exp( -5/2 *(lnSigma)^2)=26.54 
!OLD 3* exp( -5/2 *(lnSigma)^2)=0.90236
!OLD (monodisperse aerosols; 4*pi*r^2/(4/3 pi*r^3)= 3/r =88.2)
!
! According to Riemer et al, 2003, we weight the reaction probability
! according to the composition of the aerosol
!
! alpha(N2O5)=f*alpha1 +(1-f)alpha2
!   alpha1=0.02
!   alpha2=0.002
!   f= Mso4/(Mso4+Mno3), M=aerosol mass concentration
 
! N2O5 -> aerosol based upon  based on Riemer 2003 and
! (May 2011) updated based upon results shown in Riemer et al., 2009.
! We do not attempt to model OC, but simply reduce the rate by
! a factor of two to loosely account for this effect. 
! J08 - changed from use of more accurate xnew to xn_2d, since
! surface area won't change so much, and anyway the uncertainties
! are large. (and xn_2d leads to fewer dependencies)

  !---------------------------------------------------------------------
  !> FUNCTION SurfArea_Mono
  !! Mainly as sanity check for comparisons - simple spheres of diameter Dp

  function SurfArea_Mono(xn,MolWt,rho_kgm3,Dp) result(S) 
     real, intent(in) :: xn            !< molecules/cm3
     real, intent(in), optional ::&
          MolWt      &    !< mean mol.wt, g/mole
         ,rho_kgm3   &    !< density, kg/m3 
         ,Dp              !< particle diameter, wich???, here um
     real :: S            !< Surface area, cm2 per cc air
     real :: mass, rho, Rad, mw
     real :: m1, v1, s1, pn     ! Mass, volume, surface 1 particle, number of particles

     rho = 1.6                  !< g/cm3   default
     Rad = 0.034*1.0e-4         !< 0.0341 um default, in cm
     mw  = 100.0                !< g/mole, default

     if ( present(rho_kgm3) ) rho = 0.001 * rho_kgm3   !ca 1.6 g/cm3
     if ( present(Dp)       ) Rad = 0.5 * Dp
     if ( present(MolWt)    ) mw  = MolWt

     v1  = 4.0/3*PI * Rad**3    ! Vol 1 particle, cm3, in 1cc air
     m1  = v1*rho               ! Mass 1 particle, g, in 1cc air
     mass = xn*mw/AVOG ! g/cm3
     pn =  mass/m1              ! Number particles , in 1cc air
     s1   = 4.0 * PI * Rad*Rad   ! cm2 1 particle
     S    = pn * s1              ! cm2 in 1cc air
    !print "(A,99g12.3)", "SM", xn, rho, Rad, mw, pn*v1, mass, S
  end function SurfArea_Mono

  !---------------------------------------------------------------------
  !> FUNCTION SurfArea_Poly 
  !! IS NOT BE USED; so ignore !!

  function SurfArea_Poly(xn,MolWt,rho_kgm3,Dp,Dpw,sigma,sigmaFac) result(S) 
     real, intent(in) :: xn            !< molecules/cm3
     real, intent(in), optional ::&
          MolWt      &    !< mean mol.wt, g/mole
         ,rho_kgm3   &    !< density, kg/m3 
         ,Dp         &    !< particle diameter, m
         ,Dpw        &    !< wet particle diameter, m !QUERY!
         ,sigma      &    !< Sigma of mode
         ,sigmaFac        !< pre-calculated Sigma factor (saves CPU)

     real :: S            !< Surface area, m2 per m3 air
     real :: mass, rho, rdry, rwet, mw, sigFac, vol, sig
     real :: rhod, fwetvol

     rho = 1600.0                                  !< kg/m3 default
     rdry= 0.034*1.0e-6                            !< 0.0341 um default, in m
     mw  = 100.0*1.0e-3                            !< kg/mole, default

     if ( present(rho_kgm3) ) rho = rho_kgm3       ! kg/m3
     if ( present(Dp)       ) rdry= 0.5 * Dp       ! m
     if ( present(MolWt)    ) mw  = MolWt*1.0e-3   ! kg/mole
     if ( present(SigmaFac) ) then
        sigFac  = sigmaFac
     else
        sig = 1.8                                  !< default
        if ( present(sigma) ) sig = sigma
        sigFac = exp( -2.5*log(sig)**2 )
     end if

     if ( present(Dpw) )  then !  we have water
        !call CheckStop( Dpw < Dp, "Dpw<Dp1")
        rwet = 0.5*Dpw
        fwetvol = (rwet/rdry)**3  ! Ratio of wet to dry vols
    ! Moist density =  [ (fwetvol - 1) * 1000.0 + 1600.0 ] / fwetvol
         rhod=rho ! for print
         rho = ( ( fwetvol-1.0 ) * 1000.0  + rho ) / fwetvol
        print *, "Wetting ",  rhod, rho, fwetvol
     else
        rwet= rdry                 ! just for printout
        fwetvol = 1.0
     end if

     mass = xn*mw/AVOG  &              ! kg dry PM          per cc air
           +xn*(fwetvol-1)*0.018/AVOG  ! kg H2O             per cc air
     vol  = mass*1.0e6/rho             ! m3 PM              per m3 air

     ! log normal, S/V = 3/R * exp(-2.5*(log(sigma))**2)
     !S = 3.0 * vol /Rad * sigfac   ! =>  m2 polydisperse

     ! the volume calculated so far is from xn, and knows nothing about wet or
     ! dry

        S = 3.0 * vol /rwet * sigfac   ! =>  m2 polydisperse


    !print "(A,99g12.3)", "S2", xn, rho, rdry, mw, vol, mass, S
  end function SurfArea_Poly

  !---------------------------------------------------------------------
  !> FUNCTION pmSurfArea
  !! Calculates the surface area for a given mass of dry PM, together
  !! with sigma for log-normal
  !! If Dpw provided, calculates surface area of wet aerosol.
  !! NOTE - unlike emepctm's use of VOLFACs, here we use simple mass (ug/m3) and
  !! density. (rho was anyway species independent, so use of specific mol wts.
  !! for SO4, NO3, SEASALT, etc seems incosistent.)

  elemental function pmSurfArea(dry_ug,Dp,Dpw,sigma,sigmaFac,rho_kgm3) result(S) 
     real, intent(in) :: dry_ug        !< mass PM, dry, ug/m3
     real, intent(in), optional ::&
          Dp         &    !< particle diameter, m
         ,Dpw        &    !< wet particle diameter, m !QUERY!
         ,sigma      &    !< Sigma of mode
         ,sigmaFac   &    !< pre-calculated Sigma factor (saves CPU)
         ,rho_kgm3        !< density, kg/m3 

     real :: S            !< Surface area, m2 per m3 air
     real :: dryvol, rho, rdry, rwet, sigFac, totvol
     real :: rhod, fwetvol

     rho = 1600.0                                  !< kg/m3 default
     rdry= 0.034e-6                            !< 0.0341 um default, in m

! RDRY AND MASS AND RHO Should be self-consistent. No! Number!

     if ( present(rho_kgm3) ) rho = rho_kgm3       ! kg/m3
     if ( present(Dp)       ) rdry= 0.5 * Dp       ! m
     if ( present(SigmaFac) ) then
        sigFac  = sigmaFac
     else
        if ( present(sigma) )then
           sigFac = exp( -2.5*log(sigma)**2 )
        else              
           sigFac = 0.421585401578311!=exp( -2.5*log(1.8)**2 )
        endif
     end if

     rhod =rho ! for print
     if ( present(Dpw) )  then ! Dave's 1st approx.
        rwet = 0.5*Dpw
        fwetvol = (rwet/rdry)**3  ! Ratio of wet to dry vols
        rho  = ( ( fwetvol-1.0 ) * 1000.0  + rho ) / fwetvol
        ! print *, "Wetting ",  rhod, rho, fwetvol
     else
        rwet= rdry                 ! just for printout
        fwetvol = 1.0
     end if

   ! Spell out each tiny step.....

     dryvol = (dry_ug*1.0e-9)/rhod   ! m3 dry PM per m3 air
     totvol = fwetvol * dryvol       ! m3 dry+wet PM per m3 air

     ! log normal, S/V = 3/R * exp(-2.5*(log(sigma))**2)
     !S = 3.0 * vol /Rad * sigfac   ! =>  m2 polydisperse

     ! the volume calculated so far is from xn, and knows nothing about wet
     ! or dry

      S = 3.0 * totvol /rwet * sigfac   ! =>  m2 polydisperse


    !print "(A,99g12.3)", "S2", xn, rho, rdry, mw, vol, mass, S
  end function pmSurfArea

  !---------------------------------------------------------------------
  !> FUNCTION pmSurfArea
  !! Calculates the surface area for a given mass of dry PM, together
  !! with sigma for log-normal
  !! If Dpw provided, calculates surface area of wet aerosol.
  !! NOTE - unlike emepctm's earlier use of VOLFACs, here we use simple
  !! mass (ug/m3) and density. (rho was anyway species independent, so use of
  !! specific mol wts.  for SO4, NO3, SEASALT, etc seems incosistent.)

  elemental function pmH2O_gerber(dry_ug,rho_kgm3,Dp,Dpw,sigma,sigmaFac) result(ugH2O) 
     real, intent(in) :: dry_ug        !< mass PM, dry, ug/m3
     real, intent(in), optional ::&
          rho_kgm3   &    !< density, kg/m3 
         ,Dp         &    !< particle diameter, m
         ,Dpw        &    !< wet particle diameter, m !QUERY!
         ,sigma      &    !< Sigma of mode
         ,sigmaFac        !< pre-calculated Sigma factor (saves CPU)

     real :: ugH2O        !<  ug/m3 water associated with aerosol
     real :: dryvol, rho, rdry, rwet, sigFac, totvol, sig
     real :: rhod, fwetvol

     rho = 1600.0                                  !< kg/m3 default
     rdry= 0.034*1.0e-6                            !< 0.0341 um default, in m

! RDRY AND MASS AND RHO Should be self-consistent. No! Number!

     if ( present(rho_kgm3) ) rho = rho_kgm3       ! kg/m3
     if ( present(Dp)       ) rdry= 0.5 * Dp       ! m
     if ( present(SigmaFac) ) then
        sigFac  = sigmaFac
     else
        sig = 1.8                                  !< default
        if ( present(sigma) ) sig = sigma
        sigFac = exp( -2.5*log(sig)**2 )
     end if

     rhod =rho ! for print

     if ( present(Dpw) )  then ! Dave's 1st approx.
        rwet = 0.5*Dpw
        fwetvol = (rwet/rdry)**3  ! Ratio of wet to dry vols
        !rho  = ( ( fwetvol-1.0 ) * 1000.0  + rho ) / fwetvol
        ! print *, "Wetting ",  rhod, rho, fwetvol
     else
        !rwet= rdry                 ! just for printout
        fwetvol = 1.0
     end if

   ! Spell out each tiny step.....

     dryvol = (dry_ug*1.0e-9)/rhod   ! m3 dry PM per m3 air
     totvol = fwetvol * dryvol       ! m3 dry+wet PM per m3 air

     ! log normal, S/V = 3/R * exp(-2.5*(log(sigma))**2)
     !S = 3.0 * vol /Rad * sigfac   ! =>  m2 polydisperse

     ! the volume calculated so far is from xn, and knows nothing about wet
     ! or dry

      ugH2O   = (totvol-dryvol) * 1000.0 * 1.0e9  ! kg/m3
      !S = 3.0 * totvol /rwet * sigfac   ! =>  m2 polydisperse


    !print "(A,99g12.3)", "S2", xn, rho, rdry, mw, vol, mass, S
  end function pmH2O_gerber

  !----------------------------------------------------------------------
  ! Gamma for seasalt + N2O5, Evans & Jacob, 2005, Chang 2011

  elemental function GammaN2O5_EJSS(rh) result(rate) 
     real, intent(in) :: rh             !< rel. hum., 0-1
     real :: rate

     if ( rh <0.62 ) then 

          rate =  0.005
     else 
          rate = 0.03
     end if

  end function GammaN2O5_EJSS
  !---------------------------------------------------------------------
  ! Gamma for NH4NO3, Davis eqn (6)

  elemental function GammaN2O5_DavisAN(frh) result(rate) 
     real, intent(in) :: frh             !< rel. hum., 0-1
     real :: l3, rate
     real, parameter :: b30 = -8.10774, b31 = 0.04902, gmax=0.0154
     real, parameter :: b31p = 100*b31    ! Davis eqn uses RH in %

     l3 = b30 + b31p * frh 
     rate = min( gmax, 1.0/(1+exp(-l3)) )

  end function GammaN2O5_DavisAN
  !---------------------------------------------------------------------
  ! As in ACP 2012, but now elemental:

  elemental function kaero(rh) result(rate) 
     real, intent(in) :: rh  ! fractional RH
    ! Former rate for HNO3 -> NO3_c, not now used
     real :: rate
     
      if ( rh  > 0.9)  then
         rate = 1.0e-4
      else
         rate = 5.0e-6
      end if

  end function kaero
  !---------------------------------------------------------------------
  ! Gamma functions as mixture of SIA, OM, and other 

  elemental function GammaN2O5(t,frh,fso4sia,fom,fss,fdust,fbc) result(gam) 
    real, intent(in) :: t, frh   ! t(K), frh(0-1)
   ! fraction of organic matter, sea-salt, dust in aerosol
    real, intent(in) :: fso4sia, fom, fss, fdust,fbc
    real :: fsia
    real :: gam
    fsia = 1 - ( fom + fss + fdust + fbc )

     gam =  fsia *  (  fso4sia *  GammaN2O5_so4(t,frh) +&
                      (1- fso4sia) * GammaN2O5_DavisAN(frh) )
    if( fom   > 1.0e-6 ) gam = gam + fom   * GammaN2O5_om(frh)
    if( fss   > 1.0e-6 ) gam = gam + fss   * GammaN2O5_EJSS(frh)
    if( fdust > 1.0e-6 ) gam = gam + fdust * 0.01 ! EJ 2005
    if( fbc   > 1.0e-6 ) gam = gam + fbc   * 0.005! EJ 2005,kHet

  end function GammaN2O5
  !---------------------------------------------------------------------
  ! Gamma functions for N2O5 from Evans and Jacob, 2005 as 
  elemental function GammaN2O5_om(frh) result(gam) 
    real, intent(in) :: frh   ! t(K), frh(0-1)
    real :: gam

    if( frh> 0.57 ) then
       gam = 0.03
    else
       gam = 5.2e-2*frh
    end if

  end function GammaN2O5_om
  !---------------------------------------------------------------------
  ! TESTING:
  ! OM generally seems to reduce gamma, and is presumably needed partly
  ! to explain Brown's low measurements. Here we make the crude assumption
  ! that OM has the same gamma as EJ's BC. Much lower than GammaN2O5_om
  ! above.
  !elemental function GammaN2O5_ome() result(gam) 
  !  real :: gam
  !  gam = 0.005
  !end function GammaN2O5_ome
  !---------------------------------------------------------------------
  ! Gamma for sulphate from Evans & Jacob 2005
  elemental function GammaN2O5_so4(t,frh) result(gam) 
    real, intent(in) :: t, frh   ! t(K), frh(0-1)
    real :: a, b, rh
    real :: gam

    if( t> 282) then
      b = 4.0e-2*(t-294)
    else
      b = -0.48
    end if

    rh = 100*frh
    a = 2.79e-4 + 1.3e-4 * rh - 3.43e-6 * rh**2 + 7.52e-8 * rh**3

    !TYPO IN PAPER: gam = a * 10.0**b, should be MINUS b
    gam = a * 10.0**(-b)

  end function GammaN2O5_so4
  !---------------------------------------------------------------------
  elemental function UptakeRate(molSpeed,gam,S,rad) result (k)
   real, intent(in) :: molSpeed             !< mean molec. vel, m/s
   real, intent(in) :: gam                  !< gamma value
   real, intent(in) :: S                    !< Aerosol surface area, m2 per m3 air
   real, intent(in), optional :: rad        !< aerosol radius, m
   real, parameter :: Dg = 0.1 * 1.0e-4     ! 0.1 cm2/s -> m2/s
   real :: k
   if( present(rad) )  then
      k = S / ( rad/Dg  + 4/(molSpeed * gam) )
   else ! ok for > 100 nm
      k = molSpeed * gam * S /  4
   end if
   !print "(a,f8.2,4es10.2)", "Uptake terms ", S*toum2cm3, &
   !    rad, rad/Dg, 4/(molSpeed * gam), k
   
  end function UptakeRate 

  !---------------------------------------------------------------------
  ! This routine is just for offline testing.

  subroutine self_test
   real, parameter :: cm2cm3toum2cm3 = 1.0e8, m2m3toum2cm3 = 1.0e12*1.0e-6
   real, parameter :: nm=1.0e9,  um=1.0e6, um2 = m2m3toum2cm3 ! shorthands
   real :: cn2o5, Smono, Spm, Sdry, Swet, rdry, rwet, kd, kw1,kw2,kw3
   real :: DpgN, DpgV, frh, t
   integer :: ind, iRH, iTK, i
   real, dimension(10) :: ugPM, S_m2m3, Kn2o5
   integer, parameter :: io1=20,io2=22 ! TMP stallo gfortran doesn't handle newunit :-(

   cn2o5 = cMolSpeed( 298.0, 108.0)
   print *, "Speed N2O5 ", cn2o5

   ! emepctm has MMD=0.33e-6 for PMf, sigma 1.8
   DpgN = DpgV2DpgN(DpgV=0.33e-6,sigma_g=1.8)
   DpgV = DpgN2DpgV(DpgN,sigma_g=1.8) ! test reverse
   print "(2(a,f8.4))", "EMEP DpgV 0.33um -> DpgN ",DpgN*um," ->DpgV ",DpgV*um

   ! Test Gerber's eqns (NB Units MUST be m, DpgN must be number)
   rdry=0.5*DpgN
   print *, "Gerber's full scheme: rwet(rdry for RpgNd ", rdry
   print *, " 98%rh  ->Rw(nm)", wetrad(rdry, fRH=0.98, pmtype=1)/rdry
   print *, " 99.9%rh->Rw(nm)", wetrad(rdry, fRH=0.999, pmtype=1)/rdry

   !print *, "Gerber simp", wetradS(rdry=0.5*DpgN, fRH=0.98, pmtype=1) ! m

 ! Start with methods based on emepctm-like xn calcs
 ! 1 ppb SO4 ~ 4 ug/m3
   Smono = SurfArea_Mono(2.55e10,Dp=DpgV)
   Sdry  = SurfArea_Poly(2.55e10,Dp=DpgV)

 !cf Riemer had 200-600 um2/cm3
   print *, "Surf Mono 1 ppb SO4 (from cm2) ", cm2cm3toum2cm3*Smono 
   print *, "Surf Poly 1 p_Polypb SO4 (from m3) ",  m2m3toum2cm3*Sdry

   !print *, "RH98:",   wetrad( rdry=1.0e-6, fRH=98.0, pmtype=2)/rdry  !ssalt
   rdry = 1.0e-6 ! UNITS=m !! cf Fan & Toon examples QUERY!
   print *, "Gerber's full scheme: RpgNd, rdry=1um, -> rwet/rdry"   
   print *, "SS RH0.80:", wetrad( rdry, fRH=0.80, pmtype=2)/rdry !ssalt
   print *, "SS RH0.98:", wetrad( rdry, fRH=0.98, pmtype=2)/rdry !ssalt
   print *, "R  RH0.98:", wetrad( rdry, fRH=0.98, pmtype=1)/rdry !rural
   print *, "Rx2RH0.98:", wetrad( 2*rdry, fRH=0.98, pmtype=1)/(2*rdry) !rural
   print *, "R  RH0.999:", wetrad( rdry, fRH=0.999, pmtype=1)/rdry !rural
   print *, "Rx2RH0.999:", wetrad( 2*rdry, fRH=0.999, pmtype=1)/(2*rdry) !rural
   rdry = 0.5 * DpgN  ! back to EMEP aero
   print *, "Try EMEP RpgN, um:", rdry*1.0e6
   print *, "R  RH0.999:", wetrad( rdry, fRH=0.999, pmtype=1)/rdry !rural
   print *, "Rx2RH0.999:", wetrad( 2*rdry, fRH=0.999, pmtype=1)/(2*rdry) !rural
   print *, "=> Conclusion: rd/rw not very sensitive to rd"

   rdry = 0.5 * DpgN  ! back to EMEP aero
   rwet = wetrad( rdry, fRH=0.98, pmtype=1)  !seasalt
   print "(a,3f12.5)" , "Rd,w(um),w/d, EMEP RH0.98:", rdry*nm, rwet*nm, rwet/rdry
   rwet = wetrad( rdry, fRH=0.999, pmtype=1)  !seasalt
   print "(a,3f12.5)" , "Rd,w(um),w/d, EMEP RH0.999:", rdry*nm, rwet*nm, rwet/rdry
   rwet = wetrad( rdry, fRH=0.9999, pmtype=1)  !seasalt
   print "(a,3f12.5)" , "Rd,w(um),w/d, EMEP RH0.9999:", rdry*nm, rwet*nm, rwet/rdry
   rwet = wetrad( rdry, fRH=1.0, pmtype=1)  !seasalt
   print "(a,3f12.5)" , "Rd,w(um),w/d, EMEP RH1.0:", rdry*nm, rwet*nm, rwet/rdry

   rdry = 0.5 * DpgN  ! back to EMEP aero
    ! QUERY V or N forrwet/rdry and Gerber
   do ind = 2,2   !4  !  1=rural, 2=sea-salt
   print "(a,4f7.2)", " -------------- " !, 1.0e6*rdry, 1.0e6*rwet, Sdry*m2m3toum2cm3, Swet*m2m3toum2cm3
   do iRH = 80, 99, 2

     rwet   = wetrad( rdry, fRH=0.01*iRH, pmtype=ind)      ! m
!print *, "RH rwet/rdry", iRH, rdry, rwet/rdry, 2*rwet/DpgN
     Swet   = SurfArea_Poly(2.55e10,Dp=DpgN,Dpw=2*rwet)  ! xn-based approach
!print *, "Swet= ", Swet
     Spm    = pmSurfArea(dry_ug=4.0,Dp=DpgN,Dpw=2*rwet)    ! mass based
!print *, "Spm = ", Spm 

     kd     = UptakeRate(cn2o5,gam=0.01,S=Sdry,rad=rdry)    ! s-1 
     kw1    = UptakeRate(cn2o5,gam=0.01,S=Sdry,rad=rwet)    ! Using wet radius, dry S
     kw2    = UptakeRate(cn2o5,gam=0.01,S=Swet,rad=rwet)    ! Using wet radius, wet S 
     kw3    = UptakeRate(cn2o5,gam=0.01,S=Spm,rad=rwet )    ! Using wet radius, wet S 
     print "(2i3,5f8.1,a,9es11.2)", ind, iRH, rdry*nm, rwet*nm, &
         Sdry*um2, Swet*um2,Spm*um2, " ks= ", kd, kw1,kw2,kw3
   end do ! iRH
   end do ! ind

   !TMP open(newunit=io1,file="AeroSurf.txt")
   !TMP open(newunit=io2,file="AeroRate.txt")
   open(io1,file="AeroSurf.txt")
   open(io2,file="AeroRate.txt")
   ugPM = (/(1.0*i*i, i=1,size(ugPM)) /)
   do iRH = 100, 0, -10
     frh = min(99.9, 0.01*iRH)
     rwet = wetrad( rdry, 0.01*iRH )
     S_m2m3  = pmSurfArea(ugPM,Dp=2*rdry,Dpw=2*rwet)
     Kn2o5 = UptakeRate(cn2o5,gam=0.01,S=S_m2m3,rad=rwet)    ! Using wet radius, wet S 
! print *, 'High RH', frh, rwet/rdry, S_m2m3(:)*um2
!if(iRH > 97 ) print *, 'High RH', frh, rwet/rdry, S_m2m3(:)*um2
     write(io1,"(f5.2,f9.3,20g11.2)") frh, rwet/rdry, S_m2m3(:)*um2
     write(io2,"(f5.2,f9.3,20es11.2)") frh, rwet/rdry, Kn2o5(:)
     write(*,"(a,f5.2,f9.3,20g11.2)") 'RH-S:', frh, rwet/rdry, S_m2m3(:)*um2
     write(*,"(a,f5.2,f9.3,20es11.2)")'RH-K:', frh, rwet/rdry, Kn2o5(:)
   end do
   print *, "*** See AeroSurf.txt, AeroRate.txt *** "

   rdry = 0.05853 ! um
   frh  = 0.8943
   rwet = umWetRad( rdry, frh, pmtype=1 )
   Spm    = pmSurfArea(dry_ug=0.3820,Dp=2.0e-6*rdry,Dpw=2.0e-6*rwet)    ! mass based
   print *, "EMEP COMP ", 2.0e-6*rdry, 2.0e-6*rwet, rwet/rdry, Spm*1.0e6

   t= 298.0
   do iRH = 20, 100, 20
     frh = 0.01 * iRH
     print "(a,i4,9f12.3)", "Davis ", iRH, GammaN2O5_DavisAN(frh), &
      GammaN2O5_EJSS(frh), GammaN2O5_so4(t,frh), GammaN2O5_so4(t+10,frh), &
      GammaN2O5(t,frh,fso4sia=0.3,fom=0.33,fss=0.1,fdust=0.01,fbc=0.0),&
      GammaN2O5(t,frh,fso4sia=0.3,fom=0.0,fss=0.0,fdust=0.0,fbc=0.0)
   end do
   
   do iTK = 270, 290
      t=iTK
      print "(a,i5,3es12.3)", "GTK ", iTK, &
        GammaN2O5_so4(t,0.3), GammaN2O5_so4(t,0.8), GammaN2O5_so4(t,0.99)
   end do

  end subroutine self_test

  !---------------------------------------------------------------------

  subroutine self_test_fracs
    real :: Dpg, sig, Dp, FnR, FnA, fA

   print *, "Give   Dpg, sigma,   Dp "
   do
     read(*,*) Dpg, sig, Dp

     FnR =  LogNormFracBelow(Dpg,sig,Dp)
   ! Aeodynamic ?
     fA = sqrt(1.6)    ! 1.6g/cm3 - check Seinfeld ch. 8 
     FnA =  LogNormFracBelow(Dpg*fA,sig,Dp)
     print "(2f4.1,2f8.2)", Dp, Dpg,  FnR, FnA
   end do

  end subroutine self_test_fracs
end module AeroFunctions_mod
!TSTEMX program tstr
!TSTEMX use AeroFunctions_mod, only : self_test, self_test_fracs
!TSTEMX implicit none
!TSTEMX call self_test()
!TSTEMX !call self_test_fracs()
!TSTEMX end program tstr
