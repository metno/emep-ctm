! <ChemFunctions_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2015 met.no
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
! <ChemFunctions_ml.f90 - the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
module ChemFunctions_ml
!____________________________________________________________________
! Miscellaneous collection of "standard" (or guessed ) functions
! Including Troe, sine and cosine curves, 
! bilinear-interpolation routines, 
! and Standard Atmosphere p -> H conversion
!
! Where possible, reference to the EMEP documentation paper, Simpson
! et al., ACP, 2012,  are given, indicated by ACP:
! 
!____________________________________________________________________
!
!** includes
!   troe - standrad chemical function
!____________________________________________________________________
 use AeroFunctions,         only : UptakeRate, GammaN2O5_EJSS, GammaN2O5
 use CheckStop_ml,          only : StopAll
 use ChemSpecs,             only : SO4, NO3_f, NH4_f, NO3_c
 use LocalVariables_ml,     only : Grid   ! => izen, is_mainlysea
 use ModelConstants_ml,     only : K1  => KCHEMTOP, K2 => KMAX_MID, USES, AERO
 use PhysicalConstants_ml,  only : AVOG, RGAS_J, DAY_ZEN
 use Setup_1dfields_ml,     only : itemp, tinv, rh, x=> xn_2d, amk, &
     aero_fom,aero_fss,aero_fdust, &
     cN2O5, temp, DpgNw, S_m2m3 ! for surface area
 use ChemSpecs,             only : SO4, NO3_f, NH4_f, NO3_c
  implicit none
  private

  public :: troe
  public :: troeInLog  ! When log(Fc) provided
  public :: IUPAC_troe ! Using the approximate expression for F from 
                       !  Atkinson et al., 2006 (ACP6, 3625)
  !April2015 to AeroFunctions
  public ::  xkaero 
  public ::  kaero2    ! for testing
  public ::  RiemerN2O5
  public ::  S_RiemerN2O5 !TES
  public ::  HydrolysisN2O5
  public ::  ec_ageing_rate
  public ::  kmt3      ! For 3-body reactions, from Robert OCt 2009


! weighting factor for N2O5 hydrolysis
! Some help factors (VOLFAC)  pre-defined here. 0.068e-6 is
! number median radius, assumed for fine aerosol
! 1.2648 is the term 3* exp( -2.5 * (log(sig=1.8))**2 ) used below
! We also assume generic aerosol median number radius of 0.068um

  real, parameter, public :: VOLFACSO4 = 96.0/(AVOG) * 1.2648  *0.02/0.068e-6 
  real, parameter, public :: VOLFACNO3 = 62.0/(AVOG) * 1.2648  *0.02/0.068e-6 
  real, parameter, public :: VOLFACNH4 = 18.0/(AVOG) * 1.2648  *0.02/0.068e-6 


  !========================================
  contains
  !========================================

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  function kmt3(a1,c1,a3,c3,a4,c4,M) result (rckmt3)
     real, intent(in)  :: a1,c1,a3,c3,a4,c4
     real, dimension(size(tinv)), intent(in) :: m
     real, dimension(size(tinv)) :: rckmt3
     real, dimension(size(tinv)) :: k1, k3, k4
       k1 = a1 * EXP(C1*TINV)
       k3 = a3 * EXP(C3*TINV)
       k4 = a4 * EXP(C4*TINV)
       rckmt3 = k1 + (k3*amk)/(1.0+(k3*m)/k4)
       !rckmt3 = k1 + (k3*M)/(1.0+(k3*M)/k4)
  end function kmt3

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  elemental function troe(k0,kinf,Fc,M) result (rctroe)

  !+ Calculates Troe expression
  ! -----------------------------------------------------------
  ! ds note - this isn't checked or optimised yet. Taken from
  ! Seinfeld+Pandis, 1998, pp 283, eqn. 5.98. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.

     real, intent(in)  :: k0,kinf,Fc,M
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacament, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf    ! used also below
     x    = log10(y)
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

!    could have Fc already as log(Fc) to save CPU, but for now
!    keep as proper Fc. Slower but less confusing

!     rctroe = k0M / ( 1.0 + k0M/kinf) * exp(x*log(Fc))
     rctroe = k0M / ( 1.0 + y) * exp(x*log(Fc))

  end function troe
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  elemental function troeInLog(k0,kinf,LogFc,M) result (rctroe)

  !+ Calculates Troe expression
  ! -----------------------------------------------------------
  ! note - this isn't optimised yet. Taken from
  ! Seinfeld+Pandis, 1998, pp 283, eqn. 5.98. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.

     real, intent(in)  :: k0,kinf,LogFc,M
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacament, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf    ! used also below
     x    = log10(y)
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

!    give Fc already as log(Fc)

!     rctroe = k0M / ( 1.0 + k0M/kinf) * exp(x*log(Fc))
     rctroe = k0M / ( 1.0 + y) * exp(x*LogFc)

  end function troeInLog

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  elemental function IUPAC_troe(k0,kinf,Fc,M,N) result (rctroe)

  !+ Calculates Troe expression 
  ! -----------------------------------------------------------
  ! note - this isn't optimised yet. Taken from
  ! Atkinson et al. ACP 2006, 6, 3625-4055. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.
  ! NOTE that in the IUPAC nomenclature k0 already contains [M] so 
  !  the k0(IUPAC)=k0*M here
  !   N=[0.75-1.27*log10(Fc)]

     real, intent(in)  :: k0,kinf,Fc,M,N
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacement, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf    ! used also below
     x    = log10(y)/N
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

     rctroe = k0M / ( 1.0 + y) * exp(x*log(Fc))

  end function IUPAC_troe
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


! N2O5 -> nitrate calculation
!===========================================================================
! N2O5 -> nitrate calculation. Some constants for
! calculation of volume fraction of sulphate aerosol, and rate of uptake
! Mass of sulfate relative to sulfate+nitrate according to  Riemer N, 
! Vogel H, Vogel B, Schell B, Ackermann I, Kessler C, Hass H
! JGR 108 (D4): FEB 27 2003 
! 
!
! The first order reaction coefficient K (corrected for gas phase diffusion, 
! Schwartz, 1986) is given by
!
! K= S* alpha* v/4                               ACP:44
!    alpha=sticking coeff. for N2O5 =0.02
!    v=mean molecular speed for N2O5
!    S=aerosol surfac
!
! The surface area of the aerosols can be calculated as
! 
! S = V * surface/volume of aerosols
!     V=volume fraction of sulphate (cm3 aerosol/cm3 air)
!     (similar for nitrate and ammonium):
!
!     e.g. simplest form (not used) would be:
!     V = (so4 in moleculescm-3) x atw sulphate
!         ---------------------------------------------------------
!        AVOG X specific density of aerosols (assumed 2g/cm3*rh correction)
!
!    Or, shorter, V = C x M0/(AVOG*rho)
!
!    where C is conc. e.g. sulphate (molecule/cm3), M0 is molwt. 
!    We do not want to include  concentrations  or rho yet, so:
!
!     Let VOL =  M0/AVOG
!   
! E12:47
! The surface/volume ratio is calculated using Whitby particle distribution
! with number mean radius rgn=0.068  and standard deviation (Sigma)=2. 
! Then surface/volume=3/r *  exp( -5/2 *(lnSigma)^2)=26.54 
! 3* exp( -5/2 *(lnSigma)^2)=1.2648  for  sigma=1.8
! (monodisperse aerosols; 4*pi*r^2/(4/3 pi*r^3)= 3/r =88.2)
!
! Then 
!      A = VOL * C * 1.24648 /(0.068e-6*rho) 
! and
!      K = VOL * C * 1.24648 /(0.068e-6*rho) * alpha* v/4
! Set
!      VOLFAC= VOL*1.24648/0.068e-6 *alpha    
! Then
!      K = VOLFAC *C *v/(4*rho)
!
! rcmisc k=v/(4*rho) 
!
!      K = VOLFAC *rcmisc() *C
!
! According to Riemer et al, 2003, we weight the reaction probability
! according to the composition of the aerosol
!
! alpha(N2O5)=f*alpha1 +(1-f)alpha2                           ACP:45
!   alpha1=0.02
!   alpha2=0.002
!   f= Mso4/(Mso4+Mno3), M=aerosol mass concentration         ACP:46
 
! N2O5 -> aerosol based upon  based on Riemer 2003 and
! In testing, we had also tried a simple acounting for 
! results shown in Riemer et al., 2009.
! We did not attempt to model OC, but simply reduce the rate by
! a factor of two to loosely account for this effect. 
! June08 - changed from use of more accurate xnew to xn_2d, since
! surface area won't change so much, and anyway the uncertainties
! are large. (and xn_2d leads to fewer dependencies)

  function RiemerN2O5() result(rate) 
     real, dimension(K1:K2) :: rate
     real    :: rc
     real    :: f   ! Was f_Riemer
     real, parameter :: EPSIL = 1.0  ! One mol/cm3 to stop div by zero
     integer :: k
     real :: xNO3  ! As the partitioning between fine and coarse is so difficult
                   ! we include both in the nitrate used here.

     do k = K1, K2
       if ( rh(k)  > 0.4) then
          xNO3 = x(NO3_f,k) + x(NO3_c,k) 

         !mean molec speed of N2O5 (MW 108), m/s
         ! with density corrected for rh (moderate approx.)
          rc = sqrt(3.0 * RGAS_J * itemp(k) / 0.108) & ! mol.speed (m/s)
             /(4*(2.5 - rh(k)*1.25))                   ! density

          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0* xNO3  + EPSIL )


          rate(k) =  (0.9*f + 0.1) * rc *  &
                !TEST   0.5 * & ! v. loosely based on Reimer 2009 
             ( VOLFACSO4 * x(SO4,k) + VOLFACNO3 * xNO3  &
              + VOLFACNH4 * x(NH4_f,k) )    !SIA aerosol surface
        else
          rate(k) = 0.0
        endif
    end do ! k

  end function RiemerN2O5
  !---------------------------------------------------------------------

  ! crude, but rate = xxxx . S , dvs S = rate / xxxx
  function S_RiemerN2O5(k) result(S) 
     integer, intent(in) :: k
     real :: S
     real    :: c, rate, gam, rho
     real    :: f   ! Was f_Riemer
     real, parameter :: EPSIL = 1.0  ! One mol/cm3 to stop div by zero
     real :: xNO3  ! As the partitioning between fine and coarse is so difficult
                   ! we include both in the nitrate used here.

          xNO3 = x(NO3_f,k) + x(NO3_c,k) 

         !mean molec speed of N2O5 (MW 108), m/s
         ! with density corrected for rh (moderate approx.)
          c = sqrt(3.0 * RGAS_J * itemp(k) / 0.108) ! mol.speed (m/s)
          rho= (2.5 - rh(k)*1.25)                   ! density, g/cm3

          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0* xNO3  + EPSIL )


          rate=  (0.9*f + 0.1) * c /(4.0*rho) *  &
                !TEST   0.5 * & ! v. loosely based on Reimer 2009 
             ( VOLFACSO4 * x(SO4,k) + VOLFACNO3 * xNO3  &
              + VOLFACNH4 * x(NH4_f,k) )    !SIA aerosol surface
          !rate = 1/4 . gamma. c . S
          !rate in s-1
          gam = ( 0.9*f + 0.1 )*0.02
          S = 4.0 * rate/(gam * c ) ! will give S in m2

  end function S_RiemerN2O5
  !---------------------------------------------------------------------

  function HydrolysisN2O5(ormethod) result(rate) 
   character(len=*), intent(in) , optional:: ormethod ! overrides default method if wanted
   character(len=30) :: method
   real, dimension(K1:K2) :: rate
   real    :: rc
   real    :: f   ! Was f_Riemer
   real    :: gam, S, Rwet  ! for newer methods
   real, save :: g1 = 0.02, g2=0.002 ! gammas for 100% SO4, 100% NO3, default
   real, parameter :: EPSIL = 1.0  ! One mol/cm3 to stop div by zero
   integer :: k
   real :: xNO3  ! As the partitioning between fine and coarse is so difficult
                 ! we include both in the nitrate used here.
   logical, save :: first_call = .true.


   method = USES%n2o5HydrolysisMethod
   if ( present(ormethod) )then
     method = ormethod
   end if

   select case ( method )
    case ( "ORIGRIEMER","OrigRiemer") 

     do k = K1, K2
       if ( rh(k)  > 0.4) then
          xNO3 = x(NO3_f,k) + x(NO3_c,k) 

         !mean molec speed of N2O5 (MW 108), m/s
         ! with density corrected for rh (moderate approx.)
          rc = sqrt(3.0 * RGAS_J * itemp(k) / 0.108) & ! mol.speed (m/s)
             /(4*(2.5 - rh(k)*1.25))                   ! density

          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0* xNO3  + EPSIL )


          rate(k) =  (0.9*f + 0.1) * rc *  &
                !TEST   0.5 * & ! v. loosely based on Reimer 2009 
             ( VOLFACSO4 * x(SO4,k) + VOLFACNO3 * xNO3  &
              + VOLFACNH4 * x(NH4_f,k) )    !SIA aerosol surface
        else
          rate(k) = 0.0
        endif
      end do ! k
  !---------------------------------------
   case ( "Smix", "SmixTen" )

     do k = K1, K2

       if ( rh(k)  > 0.4) then ! QUERY???

            xNO3 = x(NO3_f,k) + 0.27 * x(NO3_c,k)  ! fracPM25, crude...
            f = 96*x(SO4,k)/( 96*x(SO4,k) + 62* xNO3  + EPSIL )

            S = S_m2m3(AERO%PM_F,k) !NOW all fine PM
            gam = GammaN2O5(temp(k),rh(k),&
                   f,aero_fom(k),aero_fss(k),aero_fdust(k))

            if( method == "SmixTen") gam = 0.1 * gam ! cf Brown et al, 2009!

            rate(k) = UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM
       else
            rate(k) = 0.0
       end if
    end do

      
  !---------------------------------------
   case ( "RiemerSIA", "RiemerSIAc3", "RiemerPMF", "mixedPMF" )

     do k = K1, K2

       if ( rh(k)  > 0.4) then ! QUERY???

         !Unfortunate hard-coding. Will fix in later stages
         if( method == "RiemerSIA" ) then
            S = S_m2m3(AERO%SIA_F,k)
            Rwet = 0.5*DpgNw(AERO%SIA_F,k)
         else if( method == "RiemerSIAc3" ) then
            S = S_m2m3(AERO%SIA_F,k)
            Rwet = 0.5*DpgNw(AERO%SIA_F,k)
            if( first_call ) g2=g1/3.0   ! Chang notes that the factor  of ten reduction was too high

         ! use whole aerosol area, but Riemer nitrate (factor 3 though):
         else ! if( USES%n2o5HydrolysisMethod == "RiemerPMF" ) then
              !.or.  USES%n2o5HydrolysisMethod == "mixedPMF" ) then

            S = S_m2m3(AERO%PM_F,k)
            !Rwet = 0.5*DpgNw(AERO%PM_F,k)
! Chang notes that the factor  of ten reduction was too high, and in PMF we also
! have EC, OM, etc.
            if( first_call ) g2=g1/3.0   
         end if

          xNO3 = x(NO3_f,k) + x(NO3_c,k) 

         !mean molec speed of N2O5 (MW 108), m/s
          !c=cMolSpeed(temp(k),108.0)

          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0* xNO3  + EPSIL )
          gam = g1 * f + g2 * (1-f)

          !rate(k) = UptakeRate(c,gam,S,Rwet) !1=fine SIA ! +OM
          rate(k) = UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM

          if( method == "mixedPMF" ) then

            ! 2) Add  fine sea-salt
            S = S_m2m3(AERO%SS_F,k)
            !Rwet = 0.5*DpgNw(AERO%SS_F,k)
            gam = GammaN2O5_EJSS(rh(k))
            rate(k) = rate(k) + UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM

            ! 3) Add  fine dust
            S = S_m2m3(AERO%DU_F,k)
            !Rwet = 0.5*DpgNw(AERO%DU_F,k)
            gam = 0.01 ! Evans & Jacob, 2005
            rate(k) = rate(k) + UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM

          end if


      else
         rate(k) = 0.0
      endif
    end do ! k
    case ( "Gamma:0.002")  ! Inspired by Brown et al. 2009
     do k = K1, K2

       if ( rh(k)  > 0.4) then ! QUERY???

          gam = 0.002
          S = S_m2m3(AERO%PM_F,k) !fine SIA +OM + ...
          rate(k) = UptakeRate(cN2O5(k),gam,S) 
      else
         rate(k) = 0.0
       end if

    end do ! k

     case default
       call StopAll("Unknown N2O5 hydrolysis"//method )
       
   end select
   first_call = .false.
  end function HydrolysisN2O5
  !---------------------------------------------------------------------
  function xkaero() result(rate) 
    ! Former rate for HNO3 -> NO3_c, not now used
     real, dimension(K1:K2) :: rate
     integer :: k
     
    do k = K1, K2
      if ( rh(k)  > 0.9) then
         rate(k) = 1.0e-4
      else
         rate(k) = 5.0e-6
      end if
    end do !k


  end function xkaero
  !---------------------------------------------------------------------
  function kaero2() result(rate) 
    ! New rate for HNO3 -> NO3_c, used only over sea squares
    ! as very crude simulation of sea-salt HNO3 interactions
    ! near surface (layer 16 ca. 600m).
     real, dimension(K1:K2) :: rate
     integer :: k
     
    if ( Grid%is_mainlysea) then
      rate(K1:15) = 0.0
      do k = 16, K2
        if ( rh(k)  > 0.9) then
           rate(k) = 1.0e-4
        else
           rate(k) = 5.0e-6
        end if
      end do !k
    else ! over land
      rate(K1:K2) = 0.0
    end if
  end function kaero2
 !---------------------------------------------------------------------
  function ec_ageing_rate() result(rate) 
 
   !.. Sets ageing rates for fresh EC [1/s] loosely based on Riemer etal. ACP(2004)
   !   See also Tsyro et al, JGR, 112, D23S19, 2007
   !   ---------------------------------  

     real, dimension(K1:K2) :: rate
 
    if ( Grid%izen <= DAY_ZEN ) then  ! daytime

       rate (K2-2 : K2)   = 3.5e-5  !  half-lifetime ~ 8h
       rate (K1   : K2-3) = 1.4e-4  !                ~ 2h
      else
       rate (K1 : K2 )    = 9.2e-6  !                ~ 30h
    endif

  end function ec_ageing_rate

end module ChemFunctions_ml
