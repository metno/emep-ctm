! <ChemFunctions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
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
!*****************************************************************************! 
module ChemFunctions_ml
!____________________________________________________________________
! Miscellaneous collection of "standard" (or guessed ) functions
! Including Troe, sine and cosine curves, 
! bilinear-interpolation routines, 
! and Standard Atmosphere p -> H conversion
!____________________________________________________________________
!
!** includes
!   troe - standrad chemical function
!____________________________________________________________________
 use LocalVariables_ml,     only : Grid   ! => izen, is_NWPsea
 use ModelConstants_ml,     only : K1  => KCHEMTOP, K2 => KMAX_MID
 use PhysicalConstants_ml,  only : AVOG, RGAS_J, DAY_ZEN
 use Setup_1dfields_ml,     only : itemp, tinv, rh, x=> xn_2d, amk
 use ChemSpecs_tot_ml,        only : SO4, NO3_f, NH4_f, NO3_c
  implicit none
  private

  public :: troe
  public :: troeInLog  ! When log(Fc) provided
  public :: IUPAC_troe ! Using the approximate expression for F from Atkinson et al., 2006 (ACP6, 3625)
  public ::  kaero
  public ::  kaero2    ! for testing
  public ::  RiemerN2O5
  public ::  ec_ageing_rate
  public ::  kmt3      ! For 3-body reactions, from Robert OCt 2009


! weighting factor for N2O5 hydrolysis
! Mass of sulfate relative to sulfate+nitrate
! according to  Riemer N, Vogel H, Vogel B, 
! Schell B, Ackermann I, Kessler C, Hass H
! JGR 108 (D4): FEB 27 2003 

  real, parameter, public :: VOLFACSO4 = 96.0/(AVOG) * 0.90236 *0.02/0.034e-6 
  real, parameter, public :: VOLFACNO3 = 62.0/(AVOG) * 0.90236 *0.02/0.034e-6 
  real, parameter, public :: VOLFACNH4 = 18.0/(AVOG) * 0.90236 *0.02/0.034e-6 


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
  ! ds note - this isn't checked or optimised yet. Taken from
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
  ! rb note - this isn't checked or optimised yet. Taken from
  ! Atkinson et al. ACP 2006, 6, 3625-4055. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.
  ! NOTE that in the IUPAC nomenclature k0 already contains [M] so the k0(IUPAC)=k0*M here
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
!     e.g.
!     V = (so4 in moleculescm-3) x atw sulphate
!         ---------------------------------------------------------
!        AVOG X specific density of aerosols (assumed 2g/cm3*rh correction)
!
!    Or, shorter, V = S x M0/(AVOG*rho)
!
!    where S is conc. e.g. sulphate (molecule/cm3), M0 is molwt. 
!
!
!    We do not want to include  concentrations  or rho yet, so:
!
!     Let VOL =  M0/AVOG
!   
! The surface/volume ratio is calculated using Whitby particle distribution
! with number mean radius 0.034  and standars deviation (Sigma)=2. 
! Then surface/volume=3/r *  exp( -5/2 *(lnSigma)^2)=26.54 
! 3* exp( -5/2 *(lnSigma)^2)=0.90236
! (monodisperse aerosols; 4*pi*r^2/(4/3 pi*r^3)= 3/r =88.2)
!
! Then 
!      A = VOL * S * 0.90236 /(0.034e-6*rho) 
! and
!      K = VOL * S * 0.90236 /(0.034e-6*rho)    * alpha* v/4
! Set
!      VOLFAC= VOL*0.90236/0.034e-6 *alpha    
! Then
!      K = VOLFAC *S *v/(4*rho)
!
! rcmisc k=v/(4*rho) 
!
!      K = VOLFAC *rcmisc() *S
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

          rc = sqrt(3.0 * RGAS_J * itemp(k) / 0.108) & ! mean mol. speed,m/s
             /(4*(2.5 - rh(k)*1.25)) !density, corrected for rh (moderate approx.)

          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0* xNO3  + EPSIL )


          rate(k) =  (0.9*f + 0.1) * rc *  &
             0.5 * & ! very loosely based on OC effects from Reimer 2009 
             ( VOLFACSO4 * x(SO4,k) + VOLFACNO3 * xNO3  &
              + VOLFACNH4 * x(NH4_f,k) )    !SIA aerosol surface
        else
          rate(k) = 0.0
        endif
    end do ! k

  end function RiemerN2O5
  !---------------------------------------------------------------------
  function kaero() result(rate) 
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


  end function kaero
  !---------------------------------------------------------------------
  function kaero2() result(rate) 
    ! New rate for HNO3 -> NO3_c, used only over sea squares
    ! as very crude simulation of sea-salt HNO3 interactions
    ! near surface (layer 16 ca. 600m).
     real, dimension(K1:K2) :: rate
     integer :: k
     
    if ( Grid%is_NWPsea) then
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
 
!.. Sets ageing rates for fresh EC based on Riemer et al.; ACP (2004). 

     real, dimension(K1:K2) :: rate
 
    if ( Grid%izen <= DAY_ZEN ) then  ! daytime

       rate (K2-2 : K2)   = 3.5e-5  ! t= 2h
       rate (K1   : K2-3) = 1.4e-4  ! t= 8h
      else
       rate (K1 : K2 )    = 9.2e-6  ! t= 30h
    endif

  end function ec_ageing_rate

end module ChemFunctions_ml
