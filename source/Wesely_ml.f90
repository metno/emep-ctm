! <Wesely_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Wesely_ml
!..............................................................................
! specifies data for deposition modelling using procedures recommended by
! Wesely, 1989, Atmos. Environ., 23, No.6, pp. 1293-1304
!  
!..............................................................................


! includes Wesely_tab2 for 14 gases
! specifies Henry's coefficients, reactivities for gases
!
use PhysicalConstants_ml, only : PRANDTL, Sc_H20
implicit none
private



!-------------------------------------------------------------------------
!     Table2: (variable, igas)
!       Variable:
!                 1 = DH2O/Dx                  ! ratio of diffusivities
!                 2 = H*       M atm^-1     ! effective Henry coeff.
!                 3 = pe                       !
!                 4 = k        (M s)**(-1)     ! 
!                 5 = f0                       ! 
!
public ::  Init_GasCoeff


integer, public, parameter :: NWESELY = 14   ! no. of gases in Wesely tables

  real, public, parameter,                 & ! Wesely Table 2
  dimension(5,NWESELY)  :: Wesely_tab2 =   &
  reshape (                             &
 (/                                     &
!    D     H*      pe     k     f0       
   1.9, 1.0e5,    -5.0, 9999.0, 0.0,    &! 1 = SO2      Sulphur dioxide
   1.6, 1.0e-2,   28.0,  6.0e8, 1.0,    &! 2 = O3       Ozone
   1.6, 1.0e-2, 9999.0,  2.0e6, 0.1,    &! 3 = NO2      Nitrogen dioxide
   1.3, 2.0e-3, 9999.0, 1.0e-2, 0.0,    &! 4 = NO       Nitric oxide
   1.9, 1.0e14,    7.0, 1.0e-2, 0.0,    &! 5 = HNO3     Nitric acid vapour
   1.4, 1.0e5,   23.0,     7.0, 1.0,    &! 6 = H2O2     Hydrogen peroxide
   1.6, 1.5e1,   -1.0,  9999.0, 0.0,    &! 7 = (ALD)    Acetaldehyde
   1.3, 6.0e3,   -3.0,  9999.0, 0.0,    &! 8 = HCHO     Formaldehyde
   1.6, 2.4e2, 9999.0,     2.0, 0.1,    &! 9 = (OP)     Methyl hydroperoxide
   2.0, 5.4e2, 9999.0,   6.0e2, 0.1,    &! 10 = PAA     Peroxyacetic acid
   1.6, 4.0e6,   -8.0,  9999.0, 0.0,    &! 11 = (ORA)   Formic acid
 ! followed CEH recommendation and set H* NH3 equal to sulphur
 ! (actually, CEH would have set it much higher than SO2!)
 !orig: 2.0e4, 9999.0,  9999.0, 0.0,    &! 12 = NH3     Ammonia
   1.0, 1.0e5, 9999.0,  9999.0, 0.0,    &! 12 = NH3     Ammonia
   2.6, 3.6e0, 9999.0,   3.0e3, 0.1,    &! 13 = PAN     Peroxyacetyl nitrate
   1.6, 1.0e5,    6.0,  4.0e-4, 0.1     &! 14 = HNO2    Nitrous acid
  /), &
  (/5,NWESELY/) )


!/ Ratio of diffusivites compared to ozone..

real, public, dimension(NWESELY), save :: DRx      ! Ratio D(O3)/D(x)

!/ and for the calculation of Rb we need:

real, public, dimension(NWESELY), save :: Rb_cor   ! two-thirds power of the 
                                                   ! Schmidt to Prandtl numbers

integer, public, parameter :: &
       WES_SO2 = 1,  WES_O3 = 2,  WES_NO2 = 3,  WES_NO = 4,  WES_HNO3 = 5,  &
       WES_H2O2= 6,  WES_ALD= 7,  WES_HCHO= 8,  WES_OP = 9,  WES_PAA = 10,  &
       WES_ORA = 11, WES_NH3= 12, WES_PAN = 13, WES_HNO2 = 14



contains
!==========================================================

subroutine Init_GasCoeff()

  !==========================================================
  !Description: 
  !calculates:
  ! 1) DRx - ratio of diffusivities of ozone to gas requried
  ! 2) Rb_corr -  the two-thirds power of the Schmidt to Prandtl 
  !number ratio values for all 14 gases listed in Wesely_tab2

  !==========================================================
  ! -> Calculated Rb_cor

  !Declaration of local variables

  integer :: icmp, iallwes
  real    :: Schmidt !..  number
    

  GASLOOP: do icmp = 1, NWESELY
     DRx   (icmp) = Wesely_tab2(1,WES_O3)/Wesely_tab2(1,icmp)
     Schmidt      = Sc_H20* Wesely_tab2(1,icmp)
     Rb_cor(icmp) = (Schmidt/PRANDTL)**(2.0/3.0)
  end do GASLOOP

  end subroutine Init_GasCoeff
end module Wesely_ml
