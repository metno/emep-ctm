! <GasParticleCoeffs_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2018 met.no
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
module GasParticleCoeffs_ml
!..............................................................................
! specifies data for deposition modelling. Initial parameters for Henry's
! law and reactivit scaling from:
!   Wesely, 1989, Atmos. Environ., 23, No.6, pp. 1293-1304
! extended/modified for EMEP usage.
!..............................................................................


! includes DryDepDefs for 14 gases
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


integer, public, parameter :: NDRYDEP_DEF = 64   ! no. of gases in tables below

  real, public, parameter,                 & ! Extension of Wesely Table
  dimension(5,NDRYDEP_DEF)  :: DryDepDefs =   &
  reshape (                             &
 (/                                     &
!    D     H*      pe     k     f0       
   1.9, 1.0e5,    -5.0, 9999.0, 0.0,    &! 1 = SO2      Sulphur dioxide
!-----------------------------------------------------------------------------
!ORIG   1.6, 1.0e-2,   28.0,  6.0e8, 1.0,    &! 2 = O3       Ozone
!ORIG   1.6, 1.0e-2, 9999.0,  2.0e6, 0.1,    &! 3 = NO2      Nitrogen dioxide
!RB Ozone revised D-value based on Massman 1998 - the diffusion coefficient
!    in air has a large uncertainty - no experimental determination seems to
!    be available.
!-----------------------------------------------------------------------------
     1.51, 1.0e-2,   28.0,  6.0e8, 1.0,    &! 2 = O3       Ozone  RB update
!RB Nitrogen dioxide - reactivity increased - diffusion coeff. uncertain!
!    D could be higher - based on Tang et al. 2014 D is estimated to be
!     ca 1.76 (range 1.3 - 2.7!)
     1.6, 1.0e-2, 9999.0,  2.0e6, 0.5,    &! 3 = NO2     RB update
!-----------------------------------------------------------------------------
!
  1.3, 2.0e-3, 9999.0, 1.0e-2, 0.0,    &! 4 = NO       Nitric oxide
  1.9, 1.0e14,    7.0, 1.0e-2, 0.0,    &! 5 = HNO3     Nitric acid vapour
!-----------------------------------------------------------------------------
!ORIG 1.4, 1.0e5,   23.0,     7.0, 1.0,    &! 6 = H2O2     Hydrogen peroxide
!ORIG 1.6, 1.5e1,   -1.0,  9999.0, 0.0,    &! 7 = (ALD)    Acetaldehyde
!ORIG 1.3, 6.0e3,   -3.0,  9999.0, 0.0,    &! 8 = HCHO     Formaldehyde
!ORIG 1.6, 2.4e2, 9999.0,     2.0, 0.1,    &! 9 = (OP)     Methyl hydroperoxide
!ORIG 2.0, 5.4e2, 9999.0,   6.0e2, 0.1,    &! 10 = PAA     Peroxyacetic acid
!ORIG 1.6, 4.0e6,   -8.0,  9999.0, 0.0,    &! 11 = (ORA)   Formic acid
!-----------------------------------------------------------------------------
  1.36, 1.0e5,   23.0,     7.0, 1.0,     &! 6 = H2O2   Hydrogen peroxide
   2.1, 1.3e1,   -1.0,  9999.0, 0.05,    &! 7 = ALD    Acetaldehyde - RB update
   1.4, 3.2e3,   -3.0,  9999.0, 0.2,     &! 8 = HCHO   Formaldehyde - RB update
!RB MEOOH  Methyl hydroperoxide - maybe reactivity should be higher!
   1.9, 3.0e2, 9999.0,     2.0, 0.2,    &! 9 = MEOOH  Methyl hydroperoxide
   2.4, 8.3e2, 9999.0,   6.0e2, 0.2,    &! 10 = PAA   Peroxyacetic acid - RB
! HCOOH - NOTE - may need updating - solubility is even higher at pH=7 but
!  surface resistance should perhaps not be 1/100 of that for SO2...
   1.6, 1.6e7,   -8.0,  9999.0, 0.0,    &! 11 = HCOOH   Formic acid  RB 
!-----------------------------------------------------------------------------
 ! followed CEH recommendation and set H* NH3 equal to sulphur
 ! (actually, CEH would have set it much higher than SO2!)
 !orig: 2.0e4, 9999.0,  9999.0, 0.0,    &! 12 = NH3     Ammonia
!ORIG 1.0, 1.0e5, 9999.0,  9999.0, 0.0,    &! 12 = NH3     Ammonia
!ORIG 2.6, 3.6e0, 9999.0,   3.0e3, 0.1,    &! 13 = PAN     Peroxyacetyl nitrate
!ORIG 1.6, 1.0e5,    6.0,  4.0e-4, 0.1,    &! 14 = HNO2    Nitrous acid
  1.1, 1.0e5, 9999.0,  9999.0, 0.0,    &! 12 = NH3   Ammonia RB
  2.8, 3.0e0, 9999.0,   3.0e3, 0.5,    &! 13 = PAN   Peroxyacetyl nitrate. RB
!RB HNO2 - uncertain about the f0 - Zhang et al. assume much higher reactivity
  1.6, 2.6e5,    6.0,  4.0e-4, 0.5,    &! 14 = HNO2  Nitrous acid RB
!-----------------------------------------------------------------------------
 ! Now Robert's extension to lots of organics:
   2.1, 6.0e4, 9999.0,  9999.0, 1.0,    &! 15 = HO2NO2  Pernitric acid
   2.9, 2.5e2, 9999.0,  9999.0, 1.0,    &! 16 = ANHY    Maleic anhydride (2,5-furandione)
   3.5, 1.3e4, 9999.0,  9999.0, 0.5,    &! 17 = CO2C3PAN  CH3C(O)CH2C(O)ONO3
   4.0, 1.5e8, 9999.0,  9999.0, 0.3,    &! 18 = VHISOLNO3  Very high solubility (estimated H* > ca 8.8e6 M/atm) multifunctional organic nitrates
   4.3, 5.0e6, 9999.0,  9999.0, 0.3,    &! 19 = HISOLNO3  Fairly high solubility (estimated H* ca 5e6 - 7e6 M/atm) multifunctional organic nitrates
   4.7, 5.5e4, 9999.0,  9999.0, 0.3,    &! 20 = C10H17NO4 moderately soluble C10-nitrates with an OH-group
   3.7, 5.0e4, 9999.0,  9999.0, 0.3,    &! 21 = MDNO3OH medium size moderately soluble organic nitrates with an OH-group
   2.9, 4.0e4, 9999.0,  9999.0, 0.3,    &! 22 = SMNO3OH small moderately soluble organic nitrates with an OH or OOH-group
   3.4, 2.7e4, 9999.0,  9999.0, 0.3,    &! 23 = MNO3OOH small (C3-C4) moderately soluble organic nitrates with a hydro peroxide group
   4.8, 2.2e4, 9999.0,  9999.0, 0.3,    &! 24 = C10NO3OOH moderately soluble C10-organic nitrates with a hydro peroxide group
   3.4, 1.0e4, 9999.0,  9999.0, 0.3,    &! 25 = MDSOLNO3 rather low soluble (H* ca 0.7 - 1.7e4 M/atm) organic nitrates (mixed group)
    4., 6.0e3, 9999.0,  9999.0, 0.3,    &! 26 = LOSOLNO3 low soluble (H* ca 4 - 7e3 M/atm) organic nitrates (mixed group)
   2.3, 2.0e0, 9999.0,  9999.0, 0.3,    &! 27 = CH3NO3 methyl nitrate (and ethyl nitrate)
   3.2, 1.0e0, 9999.0,  9999.0, 0.3,    &! 28 = VLSOLNO3 very low solubility (H* < ca 1e3 M/atm) organic nitrates (mixed group)
   3.8, 3.5e8, 9999.0,  9999.0, 0.2,    &! 29 = VHISOLOOH Very high solubility (estimated H* > ca 1.2e7 M/atm) multifunctional organic hydroperoxides
   2.6, 3.2e6, 9999.0,  9999.0, 0.2,    &! 30 = HCOCO3H 
   4.3, 1.6e6, 9999.0,  9999.0, 0.2,    &! 31 = LHISOLOOH Large (C7-C10) High solubility (estimated H* ca 1.6e6 M/atm) multifunctional organic hydroperoxides
   2.9, 1.2e6, 9999.0,  9999.0, 0.2,    &! 32 = SHISOLOOH Small (C2-C5) High solubility (estimated H* ca 1 - 1.4e6 M/atm) multifunctional organic hydroperoxides
   3.1, 7.2e5, 9999.0,  9999.0, 0.2,    &! 33 = RN12OOH 
   4.5, 4.4e5, 9999.0,  9999.0, 0.2,    &! 34 = PERPINONIC
   4.2, 1.8e5, 9999.0,  9999.0, 0.2,    &! 35 = NOPINAOOH
   3.3, 1.1e5, 9999.0,  9999.0, 0.2,    &! 36 = MDSOLOOH C4/C5 medium solubility (estimated H* ca 1.e5 M/atm) multifunctional organic hydroperoxides
   4.3, 9.0e4, 9999.0,  9999.0, 0.2,    &! 37 = C96OOH
   2.8, 3.1e4, 9999.0,  9999.0, 0.2,    &! 38 = HYPERACET
   4.8, 5.2e3, 9999.0,  9999.0, 0.2,    &! 39 = C10PAN2
   2.6, 4.6e3, 9999.0,  9999.0, 0.2,    &! 40 = HOCH2CO3H
   3.4, 3.0e0, 9999.0,  9999.0, 0.2,    &! 41 = MPAN
   2.7, 8.3e1, 9999.0,  9999.0, 0.2,    &! 42 = C3H7OOH
   4.4, 9.0e3, 9999.0,  9999.0, 0.05,   &! 43 = PINONALDEHYDE
   2.6, 8.0e3, 9999.0,  9999.0, 0.05,   &! 44 = ACETOL
   3.1, 1.5e3, 9999.0,  9999.0, 0.05,   &! 45 = MACROH
   2.7, 2.0e1, 9999.0,  9999.0, 0.05,   &! 46 = MEK
   3.5, 4.0e7, 9999.0,  9999.0, 0.0,    &! 47 = HISOLF0 species with estimated H* >= ca 4e7 M/atm and f0=0.0
   4.6, 1.3e7, 9999.0,  9999.0, 0.0,    &! 48 = PINONIC pinonic acid
   3.2, 5.5e6, 9999.0,  9999.0, 0.0,    &! 49 = CO23C4CHO
   3.2, 1.1e6, 9999.0,  9999.0, 0.0,    &! 50 = CARB13 - more or less guessing since no MCM equivalent to the CARB13 in the CRI scheme identified
   2.0, 7.0e5, 9999.0,  9999.0, 0.0,    &! 51 = CH3CO2H
   3.7, 3.9e5, 9999.0,  9999.0, 0.0,    &! 52 = HCC7CO
   2.1, 3.0e5, 9999.0,  9999.0, 0.0,    &! 53 = GLYOX
   3.0, 2.3e5, 9999.0,  9999.0, 0.0,    &! 54 = DICARB - mixture of C4 and C5 dicarbonyls + UCARB12 (which is not a dicarbonyl but with similar estimated D and H*
   3.45, 5.3e4, 9999.0,  9999.0, 0.0,   &! 55 = MCARB moderately soluble (estimated H* ca 4.8 - 6.1E4 M/atm) carbonyls and dicarbonyls
   2.2, 4.1e4, 9999.0,  9999.0, 0.0,    &! 56 = HOCH2CHO - glycolaldehyde
   3.1, 3.4e4, 9999.0,  9999.0, 0.0,    &! 57 = CARB12 moderately soluble carbonyls (mixed) with estimated H* ca 3.0 - 3.8E4 M/atm
   2.5, 2.4e4, 9999.0,  9999.0, 0.0,    &! 58 = MGLYOX
   3.2, 2.8e3, 9999.0,  9999.0, 0.0,    &! 59 = PHENOL
   2.4, 1.0e14, 9999.0, 9999.0, 0.0,    &! 60 = N2O5
   3.9, 1.3e7, 9999.0,  9999.0, 0.0,    &! 61 = LVASOA - to model Hodzics 0.01 anthropogenic VSOA bin
   3.1, 1.3e5, 9999.0,  9999.0, 0.0,    &! 62 = SVASOA - to model Hodzics 10, 100 and 1000ug/m3 Anthropogenic VSOA bins
   4.6, 6.3e8, 9999.0,  9999.0, 0.0,    &! 63 = LSVBSOA - to model Hodzics 0.01, 0.1, 1 and 10 ug/m3 Biogenic VSOA bins
   4.6, 3.2e7, 9999.0,  9999.0, 0.0     &! 64 = SVBSOA - to model Hodzics 100ug/m3 Biogenic VSOA bin
!END RB 
  /), &
  (/5,NDRYDEP_DEF/) )


!/ Ratio of diffusivites compared to ozone..

real, public, dimension(NDRYDEP_DEF), save :: DRx      ! Ratio D(O3)/D(x)

!/ and for the calculation of Rb we need:

real, public, dimension(NDRYDEP_DEF), save :: Rb_cor   ! two-thirds power of the 
                                                   ! Schmidt to Prandtl numbers

!RB integer, public, parameter :: &
!RB   WES_SO2 =  1, WES_O3 =  2, WES_NO2 =  3, WES_NO =  4, WES_HNO3 = 5, &
!RB   WES_H2O2=  6, WES_ALD=  7, WES_HCHO=  8, WES_OP =  9, WES_PAA = 10, &
!RB   WES_ORA = 11, WES_NH3= 12, WES_PAN = 13, WES_HNO2=14

integer, public, parameter :: &
  WES_SO2 =  1, WES_O3 =  2, WES_NO2 =  3, WES_NO =  4, WES_HNO3 = 5, &
  WES_H2O2=  6, WES_ALD=  7, WES_HCHO=  8, WES_MEOOH =  9, WES_PAA = 10, &
  WES_HCOOH = 11, WES_NH3= 12, WES_PAN = 13, WES_HNO2=14, WES_HO2NO2 = 15, &
  WES_ANHY = 16, WES_CO2C3PAN = 17, WES_VHISOLNO3 = 18, WES_HISOLNO3 = 19, &
  WES_C10H17NO4 = 20, WES_MDNO3OH = 21, WES_SMNO3OH = 22, WES_MNO3OOH = 23, &
  WES_C10NO3OOH = 24, WES_MDSOLNO3 = 25, WES_LOSOLNO3 = 26, WES_CH3NO3 = 27, &
  WES_VLSOLNO3 = 28, WES_VHISOLOOH = 29, WES_HCOCO3H = 30, WES_LHISOLOOH = 31, &
  WES_SHISOLOOH = 32, WES_RN12OOH = 33, WES_PERPINONIC = 34, &
  WES_NOPINAOOH = 35, WES_MDSOLOOH = 36,&
  WES_C96OOH = 37, WES_HYPERACET = 38, WES_C10PAN2 = 39, WES_HOCH2CO3H = 40, &
  WES_MPAN = 41, WES_C3H7OOH = 42, WES_PINONALDEHYDE = 43, WES_ACETOL = 44, &
  WES_MACROH = 45, WES_MEK = 46,  WES_HISOLF0 = 47, WES_PINONIC = 48, &
  WES_CO23C4CHO = 49, WES_CARB13 = 50, WES_CH3CO2H = 51, WES_HCC7CO = 52, &
  WES_GLYOX = 53, WES_DICARB = 54, WES_MCARB = 55, WES_HOCH2CHO = 56, &
  WES_CARB12 = 57, WES_MGLYOX = 58, WES_PHENOL = 59, WES_N2O5 = 60, &
  WES_LVASOA = 61, WES_SVASOA = 62,  WES_LSVBSOA = 63, WES_SVBSOA = 64


!*** Variables used in deposition calculations

! DDEP_xx gives the index that will be used in the EMEP model
! WES_xx gives the index of the DryDepDefs gas to which this corresponds
 
! Here we define the minimum set of species which has different
! deposition velocities. We calculate Vg for these, and then
! can use the rates for other similar species. (e.g. AMSU can use
! the Vg for SO4.  Must set NDRYDEP_CALC species

!*** IMPORTANT: the variables below must match up in the sense that, for
! example, if DDEP_NH3=4 then the 4th element of DRYDEP must be WES_NH3.

!RBinteger, public, parameter :: NDRYDEP_GASES = 11  ! gases
integer, public, parameter :: NDRYDEP_GASES = 63  ! gases

integer, public, parameter :: &
  CDDEP_HNO3 =  1, CDDEP_O3  =  2, CDDEP_SO2 = 3, &
  CDDEP_NH3  =  4, CDDEP_NO2 =  5, CDDEP_PAN = 6, &
  CDDEP_H2O2 =  7, CDDEP_ALD =  8, CDDEP_HCHO= 9, &
  CDDEP_ROOH = 10, CDDEP_HNO2= 11   !, CDDEP_PMf = 12, CDDEP_PMc = 13

! RB:
integer, public, parameter :: &
!RB  CDDEP_HNO3 =  1, CDDEP_O3  =  2, CDDEP_SO2 = 3, &
!RB  CDDEP_NH3  =  4, CDDEP_NO2 =  5, CDDEP_PAN = 6, &
!RB  CDDEP_H2O2 =  7, CDDEP_ALD =  8, CDDEP_HCHO= 9, &
!RB  CDDEP_MEOOH = 10, CDDEP_HNO2= 11,
  CDDEP_MEOOH = 10, &  !DSRD SAME AS ROOH?
  CDDEP_PAA= 12,&
  CDDEP_HCOOH= 13, CDDEP_HO2NO2=14, CDDEP_ANHY=15,&
  CDDEP_CO2C3PAN = 16, CDDEP_VHISOLNO3 = 17, CDDEP_HISOLNO3 = 18,&
  CDDEP_C10H17NO4 = 19, CDDEP_MDNO3OH = 20, &
  CDDEP_SMNO3OH = 21, CDDEP_MNO3OOH = 22, CDDEP_C10NO3OOH = 23, &
  CDDEP_MDSOLNO3 = 24, CDDEP_LOSOLNO3 = 25, &
  CDDEP_CH3NO3 = 26, CDDEP_VLSOLNO3 = 27, CDDEP_VHISOLOOH = 28,&
  CDDEP_HCOCO3H = 29, &
  CDDEP_LHISOLOOH = 30, CDDEP_SHISOLOOH = 31, CDDEP_RN12OOH = 32, CDDEP_PERPINONIC = 33, &
  CDDEP_NOPINAOOH = 34, CDDEP_MDSOLOOH = 35, CDDEP_C96OOH = 36, CDDEP_HYPERACET = 37, &
  CDDEP_C10PAN2 = 38, CDDEP_HOCH2CO3H = 39, CDDEP_MPAN = 40, CDDEP_C3H7OOH = 41, &
  CDDEP_PINONALDEHYDE = 42, CDDEP_ACETOL = 43, CDDEP_MACROH = 44, CDDEP_MEK = 45, &
  CDDEP_HISOLF0 = 46, CDDEP_PINONIC = 47, CDDEP_CO23C4CHO = 48, CDDEP_CARB13 = 49, &
  CDDEP_CH3CO2H = 50, CDDEP_HCC7CO = 51, CDDEP_GLYOX = 52, CDDEP_DICARB = 53, &
  CDDEP_MCARB = 54, CDDEP_HOCH2CHO = 55, CDDEP_CARB12 = 56, CDDEP_MGLYOX = 57, &
  CDDEP_PHENOL = 58, CDDEP_N2O5 = 59, CDDEP_LVASOA = 60, CDDEP_SVASOA = 61, &
  CDDEP_LSVBSOA = 62, CDDEP_SVBSOA = 63 !, CDDEP_PMf = 52, CDDEP_PMc = 53

integer, public, parameter :: CDDEP_RCHO = CDDEP_ALD ! Convenience

!OP renamed to MEOOH, FIN to PMf, COA to PMc
! specials for aerosols. we have 2 fine, 1 coarse and 1 'giant'type
!integer, public, parameter :: &
!  CDDEP_PMfS= 12, CDDEP_PMfN= 13, CDDEP_PMc  = 14, &
!  CDDEP_SSc = 15, CDDEP_DUc = 16, CDDEP_POLLd= 17
!integer, public, parameter :: CDDEP_PMfNH4 = 18  ! TEST_2014
!integer, public, parameter :: CDDEP_LASTPM = 18  ! Safety. Catches changes
integer, private, parameter :: NG = NDRYDEP_GASES
integer, public, parameter :: &
  CDDEP_PMfS= NG+1, CDDEP_PMfN= NG+2, CDDEP_PMc  = NG+3, &
  CDDEP_SSc = NG+4, CDDEP_DUc = NG+5
!RB , CDDEP_POLLd= NG+6
!RB integer, public, parameter :: CDDEP_PMfNH4 = NG+7  ! TEST_2014

!OP renamed to ROOH, FIN to PMf, COA to PMc
! specials for aerosols. we have 2 fine, 1 coarse and 1 'giant'type
integer, public, parameter :: &
!RB   CDDEP_PMfS= 12, CDDEP_PMfN= 13, CDDEP_PMc  = 14, &
!RB    CDDEP_SSc = 15, CDDEP_DUc = 16, &
!RN  CDDEP_BIRCH=17, CDDEP_OLIVE=18, CDDEP_GRASS=19 ! Pollen types
  CDDEP_BIRCH=NG+6, CDDEP_OLIVE=NG+7, CDDEP_GRASS=NG+8 ! Pollen types
integer, public, parameter :: CDDEP_PMfNH4 = NG+9  ! TEST_2014
integer, public, parameter :: CDDEP_LASTPM = NG+9  ! Safety. Catches changes

integer, dimension(CDDEP_PMfS:CDDEP_LASTPM), public, parameter :: &
! 1=fine,2=coarse,3=coarse sea salt, 4=dust, 5/6/7 = birch/olive/grass pollen
  AERO_SIZE = (/ 1, 1, 2, 3, 4, 5, 6, 7, 1 /)

integer, public, parameter :: NDRYDEP_AER = 9    ! aerosols with CDDEP_PMfNH4 
integer, public, parameter :: NDRYDEP_CALC = NDRYDEP_GASES + NDRYDEP_AER

integer, public, parameter :: &
  CDDEP_ASH1=CDDEP_PMfS,CDDEP_ASH2=CDDEP_PMfS,CDDEP_ASH3=CDDEP_PMfS,&
  CDDEP_ASH4=CDDEP_PMfS,CDDEP_ASH5=CDDEP_PMc ,CDDEP_ASH6=CDDEP_PMc, &
  CDDEP_ASH7=CDDEP_PMc

integer, public, parameter :: CDDEP_SET = -99

integer, public, parameter, dimension(NDRYDEP_GASES) :: &
!RB  DRYDEP_GASES = (/ WES_HNO3, WES_O3,  WES_SO2, &
!RB                    WES_NH3,  WES_NO2, WES_PAN, &
!RB                    WES_H2O2, WES_ALD, WES_HCHO, WES_OP, WES_HNO2 /)
 DRYDEP_GASES = (/ WES_HNO3, WES_O3,  WES_SO2, &
                    WES_NH3,  WES_NO2, WES_PAN, &
                    WES_H2O2, WES_ALD, WES_HCHO, WES_MEOOH, &
                    WES_HNO2, WES_PAA, WES_HCOOH, WES_HO2NO2, WES_ANHY, &
                    WES_CO2C3PAN, WES_VHISOLNO3, WES_HISOLNO3, &
                    WES_C10H17NO4, WES_MDNO3OH, WES_SMNO3OH, WES_MNO3OOH, &
                    WES_C10NO3OOH, WES_MDSOLNO3, WES_LOSOLNO3, WES_CH3NO3, &
                    WES_VLSOLNO3, WES_VHISOLOOH, WES_HCOCO3H, WES_LHISOLOOH, &
                    WES_SHISOLOOH, WES_RN12OOH, WES_PERPINONIC, WES_NOPINAOOH, &
                    WES_MDSOLOOH, WES_C96OOH, WES_HYPERACET, WES_C10PAN2, WES_HOCH2CO3H, &
                    WES_MPAN, WES_C3H7OOH, WES_PINONALDEHYDE, WES_ACETOL, &
                    WES_MACROH, WES_MEK, WES_HISOLF0, WES_PINONIC, &
                    WES_CO23C4CHO, WES_CARB13, WES_CH3CO2H, WES_HCC7CO, &
                    WES_GLYOX, WES_DICARB, WES_MCARB, WES_HOCH2CHO, WES_CARB12,&
                    WES_MGLYOX, WES_PHENOL, WES_N2O5, WES_LVASOA, &
                    WES_SVASOA, WES_LSVBSOA, WES_SVBSOA   /)


contains

!==========================================================
subroutine Init_GasCoeff()
!==========================================================
!Description:
!calculates:
! 1) DRx - ratio of diffusivities of ozone to gas requried
! 2) Rb_corr -  the two-thirds power of the Schmidt to Prandtl
!number ratio values for all 14 gases listed in DryDepDefs

!==========================================================
! -> Calculated Rb_cor

  !Declaration of local variables

  integer :: icmp
  real    :: Schmidt !..  number
    

  GASLOOP: do icmp = 1, NDRYDEP_DEF
    DRx   (icmp) = DryDepDefs(1,WES_O3)/DryDepDefs(1,icmp)
    Schmidt      = Sc_H20* DryDepDefs(1,icmp)
    Rb_cor(icmp) = (Schmidt/PRANDTL)**(2.0/3.0)
  end do GASLOOP

  end subroutine Init_GasCoeff
end module GasParticleCoeffs_ml
