! <ChemFunctions_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.5>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2024 met.no
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
! <ChemFunctions_mod.f90 - the EMEP MSC-W Chemical transport Model>
!*****************************************************************************!
module ChemFunctions_mod
!____________________________________________________________________
! Miscellaneous collection of "standard" functions for chemical
! calculations.  Includes Troe, sine and cosine curves,  and some
! from KPP system.
!
! Where possible, reference to the EMEP documentation paper, Simpson
! et al., ACP, 2012,  are given, indicated by ACP:
!
!____________________________________________________________________
!
!** includes
!   troe - standard chemical function
!____________________________________________________________________
 use AeroConstants_mod,     only: AERO
 use AeroFunctions_mod,     only: UptakeRate, GammaN2O5_EJSS, GammaN2O5
 use CheckStop_mod,         only: CheckStop, StopAll
 use ChemSpecs_mod,         only : species, species_adv
 use Config_module,  only : MasterProc, SO4_ix, NH4_f_ix, NO3_f_ix, NO3_c_ix,&
         OH_ix, O3_ix  ! For ECage, Huang
 use Config_module,     only : K1  => KCHEMTOP, K2 => KMAX_MID, USES
 use Debug_module,      only : DebugCell, DEBUG
 use Io_Progs_mod,           only : datewrite
 use LocalVariables_mod,     only : Grid   ! => izen, is_mainlysea
 use PhysicalConstants_mod,  only : AVOG, RGAS_J, DAY_ZEN
 use SmallUtils_mod,     only : find_index
 use ZchemData_mod,     only : itemp, tinv, rh, x=> xn_2d, M, &
     h2o, & ! for HuangOXD
     aero_fom,aero_fss,aero_fdust, aero_fbc,  &
     gamN2O5, cN2O5, temp, DpgNw, S_m2m3 ! for gammas & surface area
  implicit none
  private

  public :: troe
  public :: troeInLog  ! When log(Fc) provided
  public :: IUPAC_troe ! Using the approximate expression for F from
                       !  Atkinson et al., 2006 (ACP6, 3625)
  public ::  xkaero
  public ::  kaero2    ! for testing
  public ::  RiemerN2O5
  public ::  S_RiemerN2O5 !TES
  public ::  HydrolysisN2O5,HydrolysisN2O5k
  public ::  ec_ageing_rate
  public ::  kmt3      ! For 3-body reactions, from Robert OCt 2009
  public :: Chem2Index_adv, Chem2Index


! weighting factor for N2O5 hydrolysis. OLD SCHEME! NOT USED
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
  ! KMT3 uses air concenrtation (M) and inverse Temp (tinv) from Zmet
  !
  function kmt3(a1,c1,a3,c3,a4,c4) result (rckmt3)
     real, intent(in)  :: a1,c1,a3,c3,a4,c4
     real, dimension(size(M)) :: rckmt3
     real, dimension(size(M)) ::  k1, k3, k4
       k1 = a1 * EXP(C1*tinv)
       k3 = a3 * EXP(C3*tinv)
       k4 = a4 * EXP(C4*tinv)
       rckmt3 = k1 + (k3*M)/(1.0+(k3*M)/k4)
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

!DEPRECATED
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

     integer :: SO4, NO3_f, NH4_f, NO3_c

     SO4 = SO4_ix
     call CheckStop( SO4_ix<1, "SO4 not defined" )
     NH4_f = NH4_f_ix
     call CheckStop( NH4_f_ix<1, "NH4_f not defined" )
     NO3_f = NO3_f_ix
     call CheckStop( NO3_f_ix<1, "NO3_f not defined" )
     NO3_c = NO3_c_ix
     call CheckStop( NO3_c_ix<1, "NO3_c not defined" )

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
        end if
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
     integer :: SO4, NH4_f, NO3_f, NO3_c

     SO4 = SO4_ix
     call CheckStop( SO4_ix<1, "SO4 not defined" )
     NH4_f = NH4_f_ix
     call CheckStop( NH4_f_ix<1, "NH4_f not defined" )
     NO3_f = NO3_f_ix
     call CheckStop( NO3_f_ix<1, "NO3_f not defined" )
     NO3_c = NO3_c_ix
     call CheckStop( NO3_c_ix<1, "NO3_c not defined" )

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
   character(len=30), save :: method
   real, dimension(K1:K2) :: rate
   real    :: rc
   real    :: f   ! Was f_Riemer
   real    :: gam, gamSS, S,  S_ss, S_du  ! for newer methods
   real, save :: g1 = 0.02, g2=0.002 ! gammas for 100% SO4, 100% NO3, default
  ! fixed-value gammas can be specified with e.g. Gamma:0.02. We derive
  ! the numerical value, gFix, from this string
   real, save :: gFix= -999.         ! fixed-value, from Gamma:xxxx values
   character(len=20) :: gtxt         ! for Gamma:xxxx values
   real, parameter :: EPSIL = 1.0  ! One mol/cm3 to stop div by zero
   integer :: k
   real :: xNO3  ! As the partitioning between fine and coarse is so difficult
                 ! we include both in the nitrate used here.
   logical, save :: first_call = .true.
   character(len=*), parameter :: dtxt = 'HydrolN2O5:'
   integer :: SO4, NH4_f, NO3_f, NO3_c

   SO4 = SO4_ix
   call CheckStop( SO4_ix<1, "SO4 not defined" )
   NH4_f = NH4_f_ix
   call CheckStop( NH4_f_ix<1, "NH4_f not defined" )
   NO3_f = NO3_f_ix
   call CheckStop( NO3_f_ix<1, "NO3_f not defined" )
   NO3_c = NO3_c_ix
   call CheckStop( NO3_c_ix<1, "NO3_c not defined" )


   if( first_call ) then
     method = USES%n2o5HydrolysisMethod
     if ( present(ormethod) ) method = ormethod  ! WHEN is this used?
     if( method(1:6)=="Gamma:"  ) then
       gtxt=method(7:)
       read(gtxt,*) gFix
       method='gFixed'
      end if
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
        end if
      end do ! k
  !---------------------------------------
   case ( "Smix", "SmixTen" )

!if ( DEBUG%RUNCHEM .and. DebugCell ) then
!  write(*,*) dtxt//trim(method), rh(K2), S_m2m3(AERO%PM_F,K2) , S_m2m3(AERO%DU_C,K2)
!end if
     do k = K1, K2

       if ( rh(k)  > 0.4) then ! QUERY???

            xNO3 = x(NO3_f,k) + 0.27 * x(NO3_c,k)  ! fracPM25, crude...
            f = 96*x(SO4,k)/( 96*x(SO4,k) + 62* xNO3  + EPSIL )

            S = S_m2m3(AERO%PM_F,k) !NOW all fine PM
            gam = GammaN2O5(temp(k),rh(k),&
                   f,aero_fom(k),aero_fss(k),aero_fdust(k),aero_fbc(k))


            rate(k) = UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM

            !Add coarse model ! was SmixC
                 S_ss = S_m2m3(AERO%SS_C,k)
                 gamSS=GammaN2O5_EJSS(rh(k))
                 S_du = S_m2m3(AERO%DU_C,k)
                ! gamDU=0.01 ! for dust
               ! same as UptakeRate(cN2O5,gam,S), but easier to code here:
                 rate(k) = rate(k) + cN2O5(k)*(gamSS*S_ss+0.01*S_du)/4
                 ! ToDo update gam for export. Currently at fine-mod only
            !Coarse end
            if( method == "SmixTen") then
              gam = 0.1 * gam ! cf Brown et al, 2009!
              rate(k) = 0.1 * rate(k)
            end if
       else
            gam = 0.0 ! just for export
            rate(k) = 0.0
       end if
       gamN2O5(k) = gam ! just for export
    end do


  !---------------------------------------
   case ( "RiemerSIA", "RiemerSIAc3", "RiemerPMF", "mixedPMF" )

     do k = K1, K2

       if ( rh(k)  > 0.4) then ! QUERY???

         !Unfortunate hard-coding. Will fix in later stages
         if( method == "RiemerSIA" ) then
            !M24  S = S_m2m3(AERO%SIA_F,k)
            !M24 Rwet = 0.5*DpgNw(AERO%SIA_F,k)
            call StopAll(dtxt//'Deprecated:'//method)
         else if( method == "RiemerSIAc3" ) then
            !M24 S = S_m2m3(AERO%SIA_F,k)
            !M24 Rwet = 0.5*DpgNw(AERO%SIA_F,k)
            call StopAll(dtxt//'Deprecated:'//method)
            if( first_call ) g2=g1/3.0   ! Chang notes that the factor  of ten
                                         ! reduction was too high

         ! use whole aerosol area, but Riemer nitrate (factor 3 though):
         else ! if( USES%n2o5HydrolysisMethod == "RiemerPMF" ) then
              !.or.  USES%n2o5HydrolysisMethod == "mixedPMF" ) then

            S = S_m2m3(AERO%PM_F,k)
            !Rwet = 0.5*DpgNw(AERO%PM_F,k)
          ! Chang notes that the factor  of ten reduction was too high, and
          ! in PMF we also have EC, OM, etc.
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
      end if
    end do ! k
    !case ( "Gamma:0.002", "Gamma:0.05", "Gamma:0.005")  ! Inspired by Brown et al. 2009
    case ( "gFixed")  !  Fixed gammas
     do k = K1, K2

       if ( rh(k)  > 0.4) then ! QUERY???

          gam = gFix ! Found above

          S = S_m2m3(AERO%PM_F,k) !fine SIA +OM + ...
          rate(k) = UptakeRate(cN2O5(k),gam,S)
      else
         rate(k) = 0.0
         gam     = 0.0 ! just for export
      end if
      gamN2O5(k) = gam ! just for export

    end do ! k

     case default
       call StopAll("Unknown N2O5 hydrolysis"//method )

   end select
   first_call = .false.
  end function HydrolysisN2O5
  !---------------------------------------------------------------------

  function HydrolysisN2O5k(k) result(rate)
   integer, intent(in) :: k
   character(len=30), save :: method
   real    :: rc, rate
   real    :: f   ! Was f_Riemer
   real    :: gam, gamSS, S,  S_ss, S_du ! for newer methods
   real, save :: g1 = 0.02, g2=0.002 ! gammas for 100% SO4, 100% NO3, default
  ! fixed-value gammas can be specified with e.g. Gamma:0.02. We derive
  ! the numerical value, gFix, from this string
   real, save :: gFix= -999.         ! fixed-value, from Gamma:xxxx values
   character(len=20) :: gtxt         ! for Gamma:xxxx values
   real, parameter :: EPSIL = 1.0  ! One mol/cm3 to stop div by zero
   real :: xNO3  ! As the partitioning between fine and coarse is so difficult
                 ! we include both in the nitrate used here.
   logical, save :: first_call = .true.
   character(len=*), parameter :: dtxt = 'HydrolN2O5:'
   integer :: SO4, NH4_f, NO3_f, NO3_c

   SO4 = SO4_ix
   call CheckStop( SO4_ix<1, "SO4 not defined" )
   NH4_f = NH4_f_ix
   call CheckStop( NH4_f_ix<1, "NH4_f not defined" )
   NO3_f = NO3_f_ix
   call CheckStop( NO3_f_ix<1, "NO3_f not defined" )
   NO3_c = NO3_c_ix
   call CheckStop( NO3_c_ix<1, "NO3_c not defined" )


   if( first_call ) then
     method = USES%n2o5HydrolysisMethod
     if( method(1:6)=="Gamma:"  ) then
       gtxt=method(7:)
       read(gtxt,*) gFix
       method='gFixed'
      end if
   end if

   select case ( method )
    case ( "ORIGRIEMER","OrigRiemer")

       if ( rh(k)  > 0.4) then
          xNO3 = x(NO3_f,k) + x(NO3_c,k)

         !mean molec speed of N2O5 (MW 108), m/s
         ! with density corrected for rh (moderate approx.)
          rc = sqrt(3.0 * RGAS_J * itemp(k) / 0.108) & ! mol.speed (m/s)
             /(4*(2.5 - rh(k)*1.25))                   ! density

          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0* xNO3  + EPSIL )


          rate =  (0.9*f + 0.1) * rc *  &
                !TEST   0.5 * & ! v. loosely based on Reimer 2009
             ( VOLFACSO4 * x(SO4,k) + VOLFACNO3 * xNO3  &
              + VOLFACNH4 * x(NH4_f,k) )    !SIA aerosol surface
        else
          rate = 0.0
        end if

  !---------------------------------------
   case ( "Smix", "SmixTen" )

!if ( DEBUG%RUNCHEM .and. DebugCell ) then
!  write(*,*) dtxt//trim(method), rh(K2), S_m2m3(AERO%PM_F,K2) , S_m2m3(AERO%DU_C,K2)
!end if

       if ( rh(k)  > 0.4) then ! QUERY???

            xNO3 = x(NO3_f,k) + 0.27 * x(NO3_c,k)  ! fracPM25, crude...
            f = 96*x(SO4,k)/( 96*x(SO4,k) + 62* xNO3  + EPSIL )

            S = S_m2m3(AERO%PM_F,k) !NOW all fine PM
            gam = GammaN2O5(temp(k),rh(k),&
                   f,aero_fom(k),aero_fss(k),aero_fdust(k),aero_fbc(k))


            rate = UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM

            !Add coarse model ! was SmixC
                 S_ss = S_m2m3(AERO%SS_C,k)
                 gamSS=GammaN2O5_EJSS(rh(k))
                 S_du = S_m2m3(AERO%DU_C,k)
                ! gamDU=0.01 ! for dust
               ! same as UptakeRate(cN2O5,gam,S), but easier to code here:
                 rate = rate + cN2O5(k)*(gamSS*S_ss+0.01*S_du)/4
                 ! ToDo update gam for export. Currently at fine-mod only
            !Coarse end
            if( method == "SmixTen") then
              gam = 0.1 * gam ! cf Brown et al, 2009!
              rate = 0.1 * rate
            end if
       else
            gam = 0.0 ! just for export
            rate = 0.0
       end if
       gamN2O5(k) = gam ! just for export


  !---------------------------------------
   case ( "RiemerSIA", "RiemerSIAc3", "RiemerPMF", "mixedPMF" )

       if ( rh(k)  > 0.4) then ! QUERY???

         !Unfortunate hard-coding. Will fix in later stages
         if( method == "RiemerSIA" ) then
            !M24  S = S_m2m3(AERO%SIA_F,k)
            !M24 Rwet = 0.5*DpgNw(AERO%SIA_F,k)
            call StopAll(dtxt//'Deprecated:'//method)
         else if( method == "RiemerSIAc3" ) then
            !M24 S = S_m2m3(AERO%SIA_F,k)
            !M24 Rwet = 0.5*DpgNw(AERO%SIA_F,k)
            call StopAll(dtxt//'Deprecated:'//method)
            if( first_call ) g2=g1/3.0   ! Chang notes that the factor  of ten
                                         ! reduction was too high

         ! use whole aerosol area, but Riemer nitrate (factor 3 though):
         else ! if( USES%n2o5HydrolysisMethod == "RiemerPMF" ) then
              !.or.  USES%n2o5HydrolysisMethod == "mixedPMF" ) then

            S = S_m2m3(AERO%PM_F,k)
            !Rwet = 0.5*DpgNw(AERO%PM_F,k)
          ! Chang notes that the factor  of ten reduction was too high, and
          ! in PMF we also have EC, OM, etc.
            if( first_call ) g2=g1/3.0
         end if

          xNO3 = x(NO3_f,k) + x(NO3_c,k)

         !mean molec speed of N2O5 (MW 108), m/s
          !c=cMolSpeed(temp(k),108.0)

          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0* xNO3  + EPSIL )
          gam = g1 * f + g2 * (1-f)

          !rate(k) = UptakeRate(c,gam,S,Rwet) !1=fine SIA ! +OM
          rate = UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM

          if( method == "mixedPMF" ) then

            ! 2) Add  fine sea-salt
            S = S_m2m3(AERO%SS_F,k)
            !Rwet = 0.5*DpgNw(AERO%SS_F,k)
            gam = GammaN2O5_EJSS(rh(k))
            rate = rate + UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM

            ! 3) Add  fine dust
            S = S_m2m3(AERO%DU_F,k)
            !Rwet = 0.5*DpgNw(AERO%DU_F,k)
            gam = 0.01 ! Evans & Jacob, 2005
            rate = rate + UptakeRate(cN2O5(k),gam,S) !1=fine SIA ! +OM

          end if


      else
         rate = 0.0
      end if
    !case ( "Gamma:0.002", "Gamma:0.05", "Gamma:0.005")  ! Inspired by Brown et al. 2009
    case ( "gFixed")  !  Fixed gammas

       if ( rh(k)  > 0.4) then ! QUERY???

          gam = gFix ! Found above

          S = S_m2m3(AERO%PM_F,k) !fine SIA +OM + ...
          rate = UptakeRate(cN2O5(k),gam,S)
      else
         rate = 0.0
         gam     = 0.0 ! just for export
      end if
      gamN2O5(k) = gam ! just for export

     case default
       call StopAll("Unknown N2O5 hydrolysis"//method )

   end select
   first_call = .false.
  end function HydrolysisN2O5k
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
   !   Croft et al., ACP, 2005, Huang et al., ACP, 2013
   !   ---------------------------------

     real, dimension(K1:K2) :: rate
     real, parameter :: h=1.0/3600
     real, parameter :: KINF=0.015, KO3=2.0e-13, KH2O=2.1e-17 ! Huang
     real :: kACPg, kOXDg, kCCg, kOCCg  ! debug rates, ground-lev
     real :: kACPm, kOXDm, kCCm, kOCCm  ! debug rates, k=10
     real :: tauFix  ! fixed lifetime, hours

     if ( USES%ECageMethod == 'ACP2012' ) then

       if ( Grid%izen <= DAY_ZEN ) then  ! daytime

         rate (K2-2 : K2)   = 3.5e-5  !  half-lifetime ~ 8h
         rate (K1   : K2-3) = 1.4e-4  !                ~ 2h
       else
         rate (K1 : K2 )    = 9.2e-6  !                ~ 30h
       end if

     else if ( USES%ECageMethod(1:3) == 'Tau' ) then
       read(USES%ECageMethod(4:),*) tauFix  ! e.g. 24h, dvs k=  1.0/(3600*tauFix)
       rate = 1.0/(24*tauFix)
     else if ( USES%ECageMethod == 'HuangOXD' ) then

       rate(K1:K2) =  KINF*KO3*x(O3_ix,K1:K2) / &  !kOXD
               (1 + KO3*x(O3_ix,K1:K2) +KH2O * h2o(K1:K2))
       rate(K1:K2) =  rate(K1:K2) * USES%ECageFac
     else if ( USES%ECageMethod == 'HuangCC' ) then !kCC same as Liu
       rate(K1:K2) =  4.6e-12*x(OH_ix,K1:K2) + 5.8e-7 ! Liu et al.,
       rate(K1:K2) =  rate(K1:K2) * USES%ECageFac
     else if ( USES%ECageMethod == 'HuangOCC' ) then ! kCC+kOXD
       rate(K1:K2) =  4.6e-12*x(OH_ix,K1:K2) + 5.8e-7 + &
                      KINF*KO3*x(O3_ix,K1:K2) / &  !kOXD
               (1 + KO3*x(O3_ix,K1:K2) +KH2O * h2o(K1:K2))
     else
       call StopAll('Unrecognosed ECageMethod'//USES%ECageMethod)
    end if

    if (DEBUG%ECAGE .and. DebugCell ) then
       kACPg=9.2e-6
       kACPm=9.2e-6
       if( Grid%izen <= DAY_ZEN ) kACPg = 3.5e-5
       if( Grid%izen <= DAY_ZEN ) kACPm = 1.4e-4
       kOXDg=KINF*KO3*x(O3_ix,K2) / (1 + KO3*x(O3_ix,K2) +KH2O * h2o(K2))
       kOXDm=KINF*KO3*x(O3_ix,10) / (1 + KO3*x(O3_ix,10) +KH2O * h2o(10))
       kCCg=4.6e-12*x(OH_ix,K2) + 5.8e-7
       kCCm=4.6e-12*x(OH_ix,10) + 5.8e-7
       kOCCg=kCCg+kOXDg
       kOCCm=kCCm+kOXDm

       call datewrite('DBGEC', [Grid%iZen], [ x(OH_ix,K2), x(OH_ix,10), & 
         x(O3_ix,K2), h2o(K2), & 
         h/kACPg, h/kACPm, h/kOXDg, h/kOXDm, h/kCCg, h/kCCm, h/kOCCg, h/kOCCm ],&
         afmt="TXTDATE,a,i4,4es10.2,8f8.2)")
    end if ! DebugCell

  end function ec_ageing_rate

  subroutine Chem2Index_adv(species_names,species_indices,Nfound)
    !given an array of chemicals species by name ("O3", "MACRO2" etc.)
    !returns an arrays of species indices for the advected species found (IXADV_O3 etc.)
    implicit none
    character(len=*), dimension(:), intent(in) ::   species_names
    integer, dimension(:), intent(inout)::species_indices
    integer, intent(out)::Nfound
    integer :: i,index
    Nfound = 0
    do i = 1, size(species_names)
       index=find_index(trim(species_names(i)),species_adv(:)%name)
       if(index>0)then
          Nfound = Nfound + 1
          call CheckStop(Nfound>size(species_indices), "Chem2Index: species array too small")
          species_indices(Nfound) = index
       else
          if(MasterProc.and.trim(species_names(i))/='NOTSET')&
               write(*,*)'Chem2Index: '//trim(species_names(i))//' not found'
       endif
    enddo
  end subroutine Chem2Index_adv

  subroutine Chem2Index(species_names,species_indices,Nfound)
    !given an array of chemicals species by name ("OH", "MACRO2" etc.)
    !returns an arrays of species indices found.
    !Corresponding to indice in "species" array, or short lived (IXSHL_OH...)
    implicit none
    character(len=*), dimension(:), intent(in) ::   species_names
    integer, dimension(:), intent(inout)::species_indices
    integer, intent(out)::Nfound
    integer :: i,index
    Nfound = 0
    do i = 1, size(species_names)
       index=find_index(trim(species_names(i)),species(:)%name)
       if(index>0)then
          Nfound = Nfound + 1
          call CheckStop(Nfound>size(species_indices), "Chem2Index: species array too small")
          species_indices(Nfound) = index
       else
          if(MasterProc.and.trim(species_names(i))/='NOTSET')&
               write(*,*)'Chem2Index: '//trim(species_names(i))//' not found'
       endif
    enddo
  end subroutine Chem2Index


end module ChemFunctions_mod
