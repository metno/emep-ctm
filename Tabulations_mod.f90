! <Tabulations_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Tabulations_mod
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
  use PhysicalConstants_mod, only : RGAS_J, CP, T0, KAPPA
  use Config_module,    only : CHEMTMIN, CHEMTMAX  ! temperature range
  implicit none
  private


 !/- subroutines:

 public  :: tabulate           ! Sets up most tables, and calls tab_rct_rates


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
    end do

    ! Temperature-dependant rates
    !-------------------------------------------------------------------

    temp = (/ (real(i),i=CHEMTMIN,CHEMTMAX) /)   ! temp = 148..333


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

  !-------------------------------------------------------------------
  end subroutine tabulate

end module Tabulations_mod
