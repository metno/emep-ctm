! <Wesely_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
  WES_SO2 =  1, WES_O3 =  2, WES_NO2 =  3, WES_NO =  4, WES_HNO3 = 5, &
  WES_H2O2=  6, WES_ALD=  7, WES_HCHO=  8, WES_OP =  9, WES_PAA = 10, &
  WES_ORA = 11, WES_NH3= 12, WES_PAN = 13, WES_HNO2=14

!/** Variables used in deposition calculations

! DDEP_xx gives the index that will be used in the EMEP model
! WES_xx gives the index of the Wesely gas to which this corresponds
 
! Here we define the minimum set of species which has different
! deposition velocities. We calculate Vg for these, and then
! can use the rates for other similar species. (e.g. AMSU can use
! the Vg for SO4.  Must set NDRYDEP_CALC species

!/** IMPORTANT: the variables below must match up in the sense that, for
! example, if DDEP_NH3=4 then the 4th element of DRYDEP must be WES_NH3.

integer, public, parameter :: NDRYDEP_GASES = 11  ! gases
integer, public, parameter :: NDRYDEP_AER = 6   ! aerosols
integer, public, parameter :: NDRYDEP_CALC = NDRYDEP_GASES + NDRYDEP_AER

integer, public, parameter :: &
  CDDEP_HNO3 =  1, CDDEP_O3  =  2, CDDEP_SO2 = 3, &
  CDDEP_NH3  =  4, CDDEP_NO2 =  5, CDDEP_PAN = 6, &
  CDDEP_H2O2 =  7, CDDEP_ALD =  8, CDDEP_HCHO= 9, &
  CDDEP_ROOH = 10, CDDEP_HNO2= 11   !, CDDEP_PMf = 12, CDDEP_PMc = 13
integer, public, parameter :: CDDEP_RCHO = CDDEP_ALD ! Convenience
!OP renamed to ROOH, FIN to PMf, COA to PMc
! specials for aerosols. we have 2 fine, 1 coarse and 1 'giant'type
integer, public, parameter :: &
  CDDEP_PMfS= 12, CDDEP_PMfN= 13, CDDEP_PMc  = 14, &
  CDDEP_SSc = 15, CDDEP_DUc = 16, CDDEP_POLLd= 17
integer, public, parameter :: &
  CDDEP_ASH1=CDDEP_PMfS,CDDEP_ASH2=CDDEP_PMfS,CDDEP_ASH3=CDDEP_PMfS,&
  CDDEP_ASH4=CDDEP_PMfS,CDDEP_ASH5=CDDEP_PMc ,CDDEP_ASH6=CDDEP_PMc, &
  CDDEP_ASH7=CDDEP_PMc

integer, dimension(CDDEP_PMfS:CDDEP_POLLd), public, parameter :: &
  AERO_SIZE = (/ 1, 1, 2, 3, 4, 5/) !1=fine,2=coarse,3=coarse sea salt, 4=dust, 5 = pollen

integer, public, parameter :: CDDEP_SET = -99

integer, public, parameter, dimension(NDRYDEP_GASES) :: &
  DRYDEP_GASES = (/ WES_HNO3, WES_O3,  WES_SO2, &
                    WES_NH3,  WES_NO2, WES_PAN, &
                    WES_H2O2, WES_ALD, WES_HCHO, WES_OP, WES_HNO2 /)

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

  integer :: icmp
  real    :: Schmidt !..  number
    

  GASLOOP: do icmp = 1, NWESELY
    DRx   (icmp) = Wesely_tab2(1,WES_O3)/Wesely_tab2(1,icmp)
    Schmidt      = Sc_H20* Wesely_tab2(1,icmp)
    Rb_cor(icmp) = (Schmidt/PRANDTL)**(2.0/3.0)
  enddo GASLOOP

  end subroutine Init_GasCoeff
end module Wesely_ml
