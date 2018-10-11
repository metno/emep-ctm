! <CM_ChemRates_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
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
!>_________________________________________________________<

module ChemRates_rct_ml
!-----------------------------------------------------------

 
  use ChemFunctions_ml       ! => kaero, RiemerN2O5

  use AeroFunctions     ! => UpdakeRate, cMolSpeed
  use Setup_1dfields_ml ! => tinv, h2o, m, Fgas
  use Setup_1dfields_ml, m=> amk
  use ChemSpecs_tot_ml  ! => PINALD, .... for FgasJ08
  use Config_module, only: KMAX_MID,KCHEMTOP,DebugCell,DEBUG,AERO
implicit none
private

  !+ Tabulates Rate-coefficients - temperature dependant

    public :: set_rct_rates

    integer, parameter, public :: NRCT = 111   !! No. coefficients

    real, allocatable, save, public, dimension(:,:) :: rct

  contains
  !------------------------------------
  subroutine set_rct_rates() 
     logical,save::first_call=.true.
     real, dimension(KCHEMTOP:KMAX_MID) :: log300divt, logtdiv300
       log300divt(:) = log(300.0*tinv(:))
       logtdiv300(:) = log(temp(:)/300.0)
     if(first_call)then
       allocate(rct(NRCT,KCHEMTOP:KMAX_MID))
       rct=0.0
     endif
       rct(1,:) = (6.0e-34*O2+5.6E-34*N2)*O2*exp(-2.6*LOGTDIV300) 
       rct(2,:) = 1.8e-11*N2*exp(107.0*TINV) 
       rct(3,:) = 3.2e-11*O2*exp(67.0*TINV) 
       rct(4,:) = 2.14e-10*H2O 
       rct(5,:) = 1.4e-12*exp(-1310.0*TINV) 
       rct(6,:) = 1.4e-13*exp(-2470.0*TINV) 
       rct(7,:) = 1.7e-12*exp(-940.0*TINV) 
       rct(8,:) = 2.03e-16*exp(-4.57*LOG300DIVT)*exp(693.0*TINV) 
       rct(9,:) = 1.8e-11*exp(110.0*TINV) 
       rct(10,:) = 3.3e-39*exp(530/TEMP)*O2 
       rct(11,:) = 3.6e-12*exp(270.0*TINV) 
       rct(12,:) = 4.5e-14*exp(-1260.0*TINV) 
       rct(13,:) = 4.8e-11*exp(250.0*TINV) 
       rct(14,:) = 2.9e-12*exp(-160.0*TINV) 
       rct(15,:) = 7.7e-12*exp(-2100.0*TINV) 
       rct(16,:) = KMT3(2.4e-14,460.0,6.5E-34,1335.0,2.7E-17,2199.0,M) 
       rct(17,:) = (1.4e-31*M*(TEMP/300)**(-3.1)*4.0E-12)*10**(LOG10(0.4)/(1+(LOG10(1.4E-31*M*(TEMP/300)**(-3.1)/4.0E-12)/0.75-1.27*(LOG10(0.4)))**2))/(1.4E-31*M*(TEMP/300)**(-3.1)+4.0E-12) 
       rct(18,:) = (4.10e-05*M*exp(-10650/TEMP)*6.0E+15*exp(-11170/TEMP))*10**(LOG10(0.4)/(1+(LOG10(4.10E-05*M*exp(-10650/TEMP)/6.0E+15*exp(-11170/TEMP))/0.75-1.27*(LOG10(0.4)))**2))/(4.10E-05*M*exp(-10650/TEMP)+6.0E+15*exp(-11170/TEMP)) 
       rct(19,:) = 3.2e-13*exp(690/TEMP) 
       rct(20,:) = (1.0+1.4e-21*H2O*exp(2200.0*TINV))*2.2E-13*exp(600.0*TINV) 
       rct(21,:) = (1.0+1.4e-21*H2O*exp(2200.0*TINV))*1.9E-33*exp(980.0*TINV)*M 
       rct(22,:) = 2.5e-12*exp(260.0*TINV) 
       rct(23,:) = 1.85e-20*exp(2.82*LOG(TEMP))*exp(-987.0*TINV) 
       rct(24,:) = 1.44e-13+M*3.43E-33 
       rct(25,:) = 2.3e-12*exp(360.0*TINV) 
       rct(26,:) = 7.4e-13*exp(-520.0*TINV) 
       rct(27,:) = 1.03e-13*exp(365.0*TINV)-7.4E-13*exp(-520.0*TINV) 
       rct(28,:) = 6.38e-18*(TEMP**2)*exp(144.0*TINV) 
       rct(29,:) = 3.8e-13*exp(780.0*TINV) 
       rct(30,:) = 5.3e-12*exp(190.0*TINV) 
       rct(31,:) = 1.25e-17*(TEMP**2)*exp(615.0*TINV) 
       rct(32,:) = 2e-12*exp(-2440.0*TINV) 
       rct(33,:) = 6.9e-12*exp(-1000.0*TINV) 
       rct(34,:) = 2.55e-12*exp(380.0*TINV) 
       rct(35,:) = 3.8e-13*exp(900.0*TINV) 
       rct(36,:) = 1.9e-12*exp(190.0*TINV) 
       rct(37,:) = 4.4e-12*exp(365.0*TINV) 
       rct(38,:) = 7.5e-12*exp(290.0*TINV) 
       rct(39,:) = 2e-12*exp(500.0*TINV) 
       rct(40,:) = 2.9e-12*exp(500.0*TINV) 
       rct(41,:) = 5.2e-13*exp(980.0*TINV) 
       rct(42,:) = 6.7e-18*(TEMP**2)*exp(511.0*TINV) 
       rct(43,:) = 2.03e-17*(TEMP**2)*exp(78.0*TINV) 
       rct(44,:) = 2.54e-12*exp(360.0*TINV) 
       rct(45,:) = 1.81875e-13*exp(1300.0*TINV) 
       rct(46,:) = 2.53e-18*(TEMP**2)*exp(503.0*TINV) 
       rct(47,:) = 9.1e-15*exp(-2580.0*TINV) 
       rct(48,:) = 5.5e-15*exp(-1880.0*TINV) 
       rct(49,:) = 1.5132e-13*exp(1300.0*TINV) 
       rct(50,:) = 2.49969e-13*exp(1300.0*TINV) 
       rct(51,:) = 2.05446e-13*exp(1300.0*TINV) 
       rct(52,:) = 6.6e-18*(TEMP**2)*exp(820.0*TINV) 
       rct(53,:) = 1.9e-12*exp(575.0*TINV) 
       rct(54,:) = 2.7e-11*exp(390*TINV) 
       rct(55,:) = 3.4299e-15*exp(-1995*TINV) 
       rct(56,:) = 3.15e-12*exp(-450*TINV) 
       rct(57,:) = 2.286e-12*exp(360.0*TINV) 
       rct(58,:) = 2.54e-13*exp(360.0*TINV) 
       rct(59,:) = 1.3e-12*exp(610.0*TINV) 
       rct(60,:) = 4e-12*exp(380.0*TINV) 
       rct(61,:) = 2.13e-16*exp(-1520.0*TINV) 
       rct(62,:) = 3.5e-16*exp(-2100.0*TINV) 
       rct(63,:) = 1.27e-12*exp(360.0*TINV) 
       rct(64,:) = 1.6e-12*exp(305*TINV) 
       rct(65,:) = IUPAC_TROE(3.28e-28*exp(6.87*LOG300DIVT)  &
         ,1.125E-11*exp(1.105*LOG300DIVT)  &
         ,0.3  &
         ,M  &
         ,0.75-1.27*LOG10(0.3))*0.107 
       rct(66,:) = IUPAC_TROE(1.1e-5*exp(-10100.0*TINV)  &
         ,1.9E17*exp(-14100.0*TINV)  &
         ,0.3  &
         ,M  &
         ,0.75-1.27*LOG10(0.3)) 
       rct(67,:) = 97760000*exp(-7261*TINV) 
       rct(68,:) = 1450000000000.*exp(-10688*TINV) 
       rct(69,:) = 0.065**2 
       rct(70,:) = HYDROLYSISN2O5() 
       rct(71,:) = IUPAC_TROE(1.0e-31*exp(1.6*LOG300DIVT)  &
         ,5.0E-11*exp(+0.3*LOG300DIVT)  &
         ,0.85  &
         ,M  &
         ,0.75-1.27*LOG10(0.85)) 
       rct(72,:) = IUPAC_TROE(3.6e-30*exp(4.1*LOG300DIVT)  &
         ,1.9E-12*exp(-0.2*LOG300DIVT)  &
         ,0.35  &
         ,M  &
         ,0.75-1.27*LOG10(0.35)) 
       rct(73,:) = IUPAC_TROE(1.3e-3*exp(3.5*LOG300DIVT)*exp(-11000.0*TINV)  &
         ,9.70E14*exp(-0.1*LOG300DIVT)*exp(-11080.0*TINV)  &
         ,0.35  &
         ,M  &
         ,0.75-1.27*LOG10(0.35)) 
       rct(74,:) = IUPAC_TROE(3.2e-30*exp(4.5*LOG300DIVT)  &
         ,3.0E-11  &
         ,0.41  &
         ,M  &
         ,0.75-1.27*LOG10(0.41)) 
       rct(75,:) = IUPAC_TROE(3.28e-28*exp(6.87*LOG300DIVT)  &
         ,1.125E-11*exp(1.105*LOG300DIVT)  &
         ,0.3  &
         ,M  &
         ,0.75-1.27*LOG10(0.3)) 
       rct(76,:) = IUPAC_TROE(8.6e-29*exp(3.1*LOG300DIVT)  &
         ,9.0E-12*exp(0.85*LOG300DIVT)  &
         ,0.48  &
         ,M  &
         ,0.75-1.27*LOG10(0.48)) 
       rct(77,:) = IUPAC_TROE(8.0e-27*exp(3.5*LOG300DIVT)  &
         ,3.0E-11*300.0*TINV  &
         ,0.5  &
         ,M  &
         ,0.75-1.27*LOG10(0.5)) 
       rct(78,:) = IUPAC_TROE(7.4e-31*exp(2.4*LOG300DIVT)  &
         ,3.3E-11*exp(0.3*LOG300DIVT)  &
         ,0.81  &
         ,M  &
         ,0.75-1.27*LOG10(0.81)) 
       rct(79,:) = UPTAKERATE(CNO3,GAM=0.001,S=S_M2M3(AERO%PM,:)) 
       rct(80,:) = UPTAKERATE(CHNO3(:),GAM=0.1,S=S_M2M3(AERO%DU_C,:)) 
       rct(81,:) = UPTAKERATE(CHNO3(:),GAM=0.01,S=S_M2M3(AERO%SS_C,:)) 
       rct(82,:) = UPTAKERATE(CHO2,GAM=0.2,S=S_M2M3(AERO%PM,:)) 
       rct(83,:) = 1.2e-11*exp(440*TINV) 
       rct(84,:) = 2.38e-11*exp(357*TINV) 
       rct(85,:) = 3.948e-11*exp(440*TINV) 
       rct(86,:) = 7.5e-13*exp(700*TINV) 
       rct(87,:) = 3.8e-12*exp(200*TINV) 
       rct(88,:) = 8.05e-16*exp(-640*TINV) 
       rct(89,:) = 1.35e-15*exp(-1270*TINV) 
       rct(90,:) = 2.6887e-15*exp(-640*TINV) 
       rct(91,:) = 1.2e-12*exp(490*TINV) 
       rct(92,:) = 2.6988e-12*exp(490*TINV) 
       rct(93,:) = 1.2e-12*exp(490.0*TINV) 
       rct(94,:) = 4e-12*FGAS(ASOC_UG1,:) 
       rct(95,:) = 4e-12*FGAS(ASOC_UG10,:) 
       rct(96,:) = 4e-12*FGAS(ASOC_UG1E2,:) 
       rct(97,:) = 4e-12*FGAS(ASOC_UG1E3,:) 
       rct(98,:) = 4e-12*FGAS(NON_C_ASOA_UG1,:) 
       rct(99,:) = 4e-12*FGAS(NON_C_ASOA_UG10,:) 
       rct(100,:) = 4e-12*FGAS(NON_C_ASOA_UG1E2,:) 
       rct(101,:) = 4e-12*FGAS(NON_C_ASOA_UG1E3,:) 
       rct(102,:) = 4e-12*FGAS(BSOC_UG1,:) 
       rct(103,:) = 4e-12*FGAS(BSOC_UG10,:) 
       rct(104,:) = 4e-12*FGAS(BSOC_UG1E2,:) 
       rct(105,:) = 4e-12*FGAS(BSOC_UG1E3,:) 
       rct(106,:) = 4e-12*FGAS(NON_C_BSOA_UG1,:) 
       rct(107,:) = 4e-12*FGAS(NON_C_BSOA_UG10,:) 
       rct(108,:) = 4e-12*FGAS(NON_C_BSOA_UG1E2,:) 
       rct(109,:) = 4e-12*FGAS(NON_C_BSOA_UG1E3,:) 
       rct(110,:) = EC_AGEING_RATE() 
       rct(111,:) = EC_AGEING_RATE() 
       first_call=.false.
  end subroutine set_rct_rates
end module  ChemRates_rct_ml
