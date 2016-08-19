! <Tabulations_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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
module Tabulations_ml
 !+
 ! Tabulates miscellaneous functions and reaction rates which depend on
 ! temperature: 
 !
 !    exner function (tpi), sat. vapour pressure, rh deliq.
 !    Mozurk. parameters.
 !
 ! For gas-phase chemistry, only the "simple" temperature dependant rates
 ! are tabulated here.
 !
 !----------------------------------------------------------------------------
  use PhysicalConstants_ml, only : RGAS_J, CP, T0, KAPPA
  use ModelConstants_ml,    only : CHEMTMIN, CHEMTMAX  ! temperature range
  !u3 use GenRates_rct_ml, only : &  ! From GenOut_ml 
  !u3                NRCT, &         ! No. temperature dependant coefficients
  !u3                rcvt, &         ! Temperature dependant coefficients
  !u3                set_rct_rates   ! Gives RCT as function of temp t
  implicit none
  private


 !/- subroutines:

 public  :: tabulate           ! Sets up most tables, and calls tab_rct_rates
 !u3 private :: tab_rct_rates      ! Subroutine to tabulate rct to rcit  rates


 !/- Outputs:

  real, public, parameter  ::    &
        PINC=1000.0              &  
     ,  PBAS=-PINC

  real, save, public, dimension(131) ::  tpi   ! Exner function of pressure?

 real, save, public, dimension(CHEMTMIN:CHEMTMAX) :: &
                   tab_esat_Pa !&  ! saturated vapour pressure (Pa)
!                  ,tab_rhdel   &  ! RH of deliquescence for ammonium nitrate
!                  ,tab_Kp_amni &  ! Equil. constant, nh3 + hno3  <--> nh4no3
!                  ,tab_MozP1   &  ! Mozurkewich P1 value for Kaq
!                  ,tab_MozP2   &  ! Mozurkewich P2 value for Kaq
!                  ,tab_MozP3   &  ! Mozurkewich P3 value for Kaq
!                  ,tab_vav_n2o5   ! avg. molecular speed N2O5

!u3  !/ Output gas-phase chemical rates:
!u3 
!u3   real, save, public, &
!u3          dimension(NRCT,CHEMTMIN:CHEMTMAX) :: rcit  ! rate-coefficients


 contains

 subroutine tabulate()
 !

    real, dimension(CHEMTMIN:CHEMTMAX) :: temp
    real    :: p, a  
    integer :: i

    !  Exner function
    !-------------------------------------------------------------------
    ! define the exner-function for every 1000 pa from zero to 1.3e+5 pa
    ! in a table for efficient interpolation (same procedure as used in
    ! the nwp-model, see mb1e.f)
    !

	do i = 1,131
	  p = PBAS + i*PINC
	  tpi(i) = CP*(p/1.0e+5)**KAPPA
	enddo

    ! Temperature-dependant rates
    !-------------------------------------------------------------------

    temp = (/ (real(i),i=CHEMTMIN,CHEMTMAX) /)   ! temp = 148..333

    !u3 call tab_rct_rates()  ! = > gives gas-phase rates

    ! Tabulation of other rates:
    !-------------------------------------------------------------------
    !  Saturation vapour pressure
    !  Clausius-Clapyron, as given by Jakobsen, eqn 2.55 (now in Pa):

    !   T0  = 273.15
    !   tab_esat_Pa(:) = 611.2*exp( 6816.0*(1.0/T0 + 1.0/temp(:))  &
    !                               + 5.1309 * (T0/temp(:)      )

    !   where T0 is standrard temperature 273.15
    !   Ref: Bolton's formula - really only valid for -35C < T < 35C 
    !   Units: Pa
    !   (MADE/MACHO notes  : was miscrcit(ICES,it)
    !   Corrected  to Pa, 30/10/01
    !
    !tab_esat_Pa(:) = 611.2*exp(17.67*(temp(:)-273.15)/ (temp(:) - 29.65))

    ! From Bohren+Albrecht, Atmospheric Thermodynamics, 1998
    ! ln e/es = 6808(1/T0-1/T) - 5.09 ln(T/T0)
    !=> e = es * exp( 6808(1/T0-1/T)) * (T/T0)**5.09
    !
    !  tab_ba(:) = 611.2* ( exp( 6808.0*(1.0/T0 - 1.0/temp(:)) )  &
    !                               * ( T0/temp(:) )**5.09  )

    ! But, for now  I chose the formula used by HIRLAM; as provided
    ! by Anna, May 2002
    ! with  eps*zxl/Rd in Clausius-Clapyron

       a = 0.622*2.5e6/287.0   ! = 5418

       tab_esat_Pa(:) = 611.2* exp( a *(1.0/T0 - 1.0/temp(:)) ) 

    !
    ! An alternative would be to use 
    !-------------------------------------------------------------------
    !   relative humidity of deliquescence for ammonium nitrate
    !   Ref:  Mozurkewich (1993)  - Journal???
    !   Units : fraction 0-1
    !   (MADE/MACHO notes  : was miscrcit(ICRHD,it)

!    tab_rhdel(:) = exp( 618.3/temp(:) - 2.551 )

    !-------------------------------------------------------------------
    !    Equilibrium constant (Kp):  NH3 + HNO3  <-------> NH4NO3   
    !    Ref: Mozurkewich (1993)
    !    Units : (molecule/cm3)^2 for Kp
    !   (MADE/MACHO notes  : was miscrcit(ICRS,it)
    !
    !      lnKp = 118.87 - 24084.0/T - 6.025* ln(T)
    !

!    tab_Kp_amni(:) = exp( 118.87 - 24084.0/temp(:)-6.025*log(temp(:)) )

    !-------------------------------------------------------------------
    !    temp. dependant constrants for calcolating dissos. rate 
    !    for  the formation of ammonium nitrate  
    !    Ref: Mozurkewich (1993)
    !   (MADE/MACHO notes  : was miscrcit(ICXK1,it)..miscrcit(ICXK_3,it)
    !    n.b. EMEP report 2/98 had 2446 in P3, but 24.46 is correct

!    tab_MozP1(:) = exp( -135.94 +  8763.0/temp(:) + 19.12*log( temp(:) ) )
!    tab_MozP2(:) = exp( -122.65 +  9969.0/temp(:) + 16.22*log( temp(:) ) )
!    tab_MozP3(:) = exp( -182.61 + 13875.0/temp(:) + 24.46*log( temp(:) ) )

    !-------------------------------------------------------------------
    ! vav_n2o5 is the mean molecular speed for n2o5
    !             - calculated as v = sqrt(3RT/atw)
    ! Units:  m/s
    !
    !  (MADE/MACHO: was miscrcit(IC42H,it), with units of cm/s, calculated
    !   as sqrt(3.0 * 8.314e7 * temp(:) / 108.0)
    ! RGAS_J = 8.314 J mol-1 K-1  = 8.314 kg m2 s-2 mol-1 K-1
    ! nb. Seinfeld+Pandis (1998), p.453,  use v = sqrt( 8RT/(pi*atw) )
    !
    ! Note: atwn2o5 = 108 g = 0.108 kg

!    tab_vav_n2o5(:) = sqrt(3.0 * RGAS_J * temp(:) / 0.108)  ! m/s !


  !-------------------------------------------------------------------
  end subroutine tabulate

 !u3   !----------------------------------------------------------------------
 !u3   subroutine tab_rct_rates()
 !u3   !+1) Temperature-dependant rates (rct). Only needs to be called once
 !u3   !    at beginning of simulations to set up table
 !u3 
 !u3       integer :: it          !   Local loop variable
 !u3       real    ::  tinv          !   temperature in K
 !u3 
 !u3       do it = CHEMTMIN, CHEMTMAX
 !u3         tinv = 1./real(it)
 !u3         call set_rct_rates(tinv)
 !u3         rcit(:,it) = rcvt(:)
 !u3       end do
 !u3 
 !u3   end subroutine tab_rct_rates
 !u3   !----------------------------------------------------------------------

end module Tabulations_ml
