!>_________________________________________________________<

  module  ChemRates_rcmisc_ml
!-----------------------------------------------------------

  
  use ChemFunctions_ml       ! => kaero, RiemerN2O5
  use Setup_1dfields_ml      ! => tinv, h2o, m, Fgas
  use Setup_1dfields_ml, m=> amk
  use ChemSpecs_tot_ml         ! => PINALD, .... for FgasJ08
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP,DebugCell,DEBUG_RUNCHEM
  implicit none
  private

  !+ Tabulates Rate-coefficients - complex dependancies 

    public :: set_rcmisc_rates

    integer, parameter, public :: NRCMISC = 30   !! No. coefficients

    real, save, public, dimension(NRCMISC,KCHEMTOP:KMAX_MID) :: rcmisc 

  contains
  !------------------------------------
  subroutine set_rcmisc_rates() 
!OLD real, dimension(KCHEMTOP:KMAX_MID) :: lt300
     real, dimension(KCHEMTOP:KMAX_MID) :: log300divt, logtdiv300
!     real, dimension(KCHEMTOP:KMAX_MID) :: BranchingNO
!     real, dimension(KCHEMTOP:KMAX_MID) :: BranchingHO2
!OLD       lt300(:) = log(300.0*tinv(:))
       log300divt(:) = log(300.0*tinv(:))
       logtdiv300(:) = log(temp(:)/300.0)
!       BranchingNO(:) = 1.0e-11*xn_2d(NO,:)/ &
                !( 1.0e-11*xn_2d(NO,:) + 4.2e-12*exp(180*TINV(:))*xn_2d(HO2,:) )
!       BranchingHO2(:) = 1.0 - BranchingNO(:)

       rcmisc(1,:) = (6.0e-34*o2+5.6e-34*n2)*o2*exp(-2.6*logtdiv300) 
       rcmisc(2,:) = 1.8e-11*n2*exp(107.0*tinv) 
       rcmisc(3,:) = 3.2e-11*o2*exp(67.0*tinv) 
       rcmisc(4,:) = 2.2e-10*h2o 
       rcmisc(5,:) = 2.03e-16*exp(-4.57*log300divt)*exp(693.0*tinv) 
       rcmisc(6,:) = kmt3(2.4e-14,460.0,6.5e-34,1335.0,2.7e-17,2199.0,m) 
       rcmisc(7,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*2.2e-13*exp(600.0*tinv) 
       rcmisc(8,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*1.9e-33*exp(980.0*tinv)*m 
       rcmisc(9,:) = 1.85e-20*exp(2.82*log(temp))*exp(-987.0*tinv) 
       rcmisc(10,:) = 1.44e-13+m*3.43e-33 
       rcmisc(11,:) = 6.38e-18*(temp**2)*exp(144.0*tinv) 
       rcmisc(12,:) = 1.25e-17*(temp**2)*exp(615.0*tinv) 
       rcmisc(13,:) = 6.7e-18*(temp**2)*exp(511.0*tinv) 
       rcmisc(14,:) = 2.03e-17*(temp**2)*exp(78.0*tinv) 
       rcmisc(15,:) = 2.53e-18*(temp**2)*exp(503.0*tinv) 
       rcmisc(16,:) = 6.6e-18*(temp**2)*exp(820.0*tinv) 
       rcmisc(17,:) = riemern2o5() 
       rcmisc(18,:) = 1e-12*h2o 
       rcmisc(19,:) = kaero() 
       rcmisc(20,:) = iupac_troe(1.0e-31*exp(1.6*log300divt)  &
         ,3.0e-11*exp(-0.3*log300divt)  &
         ,0.85  &
         ,m  &
         ,0.75-1.27*log10(0.85)) 
       rcmisc(21,:) = iupac_troe(3.6e-30*exp(4.1*log300divt)  &
         ,1.9e-12*exp(-0.2*log300divt)  &
         ,0.35  &
         ,m  &
         ,0.75-1.27*log10(0.35)) 
       rcmisc(22,:) = iupac_troe(1.3e-3*exp(3.5*log300divt)*exp(-11000.0*tinv)  &
         ,9.70e14*exp(-0.1*log300divt)*exp(-11080.0*tinv)  &
         ,0.35  &
         ,m  &
         ,0.75-1.27*log10(0.35)) 
       rcmisc(23,:) = iupac_troe(3.3e-30*exp(3.0*log300divt)  &
         ,4.1e-11  &
         ,0.40  &
         ,m  &
         ,0.75-1.27*log10(0.4)) 
       rcmisc(24,:) = iupac_troe(2.7e-28*exp(7.1*log300divt)  &
         ,1.2e-11*exp(0.9*log300divt)  &
         ,0.3  &
         ,m  &
         ,0.75-1.27*log10(0.3)) 
       rcmisc(25,:) = iupac_troe(4.9e-3*exp(-12100.0*tinv)  &
         ,5.4e16*exp(-13830.0*tinv)  &
         ,0.3  &
         ,m  &
         ,0.75-1.27*log10(0.3)) 
       rcmisc(26,:) = iupac_troe(8.6e-29*exp(3.1*log300divt)  &
         ,9.0e-12*exp(0.85*log300divt)  &
         ,0.48  &
         ,m  &
         ,0.75-1.27*log10(0.48)) 
       rcmisc(27,:) = iupac_troe(8.0e-27*exp(3.5*log300divt)  &
         ,3.0e-11*300.0*tinv  &
         ,0.5  &
         ,m  &
         ,0.75-1.27*log10(0.5)) 
       rcmisc(28,:) = iupac_troe(2.7e-28*exp(7.1*log300divt)  &
         ,1.2e-11*exp(0.9*log300divt)  &
         ,0.3  &
         ,m  &
         ,0.75-1.27*log10(0.3)) 
       rcmisc(29,:) = iupac_troe(4.9e-3*exp(-12100.0*tinv)  &
         ,5.4e16*exp(-13830.0*tinv)  &
         ,0.3  &
         ,m  &
         ,0.75-1.27*log10(0.3)) 
       rcmisc(30,:) = iupac_troe(7.4e-31*exp(2.4*log300divt)  &
         ,3.3e-11*exp(0.3*log300divt)  &
         ,exp(-temp/1420.0)  &
         ,m  &
         ,0.75+3.884e-4*temp) 

  end subroutine set_rcmisc_rates
end module  ChemRates_rcmisc_ml
!>_________________________________________________________<

  module  ChemRates_rct_ml
!-----------------------------------------------------------

  
  use Setup_1dfields_ml      ! => tinv, h2o, m, Fgas
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP
  implicit none
  private

  !+ Tabulates Rate-coefficients - temperature dependant 

    public :: set_rct_rates

    integer, parameter, public :: NRCT = 44   !! No. coefficients

    real, save, public, dimension(NRCT,KCHEMTOP:KMAX_MID) :: rct 

  contains
  !------------------------------------
  subroutine set_rct_rates() 
       rct(1,:) = 1.4e-12*exp(-1310.0*TINV) 
       rct(2,:) = 1.4e-13*exp(-2470.0*TINV) 
       rct(3,:) = 1.7e-12*exp(-940.0*TINV) 
       rct(4,:) = 1.8e-11*exp(110.0*TINV) 
       rct(5,:) = 3.6e-12*exp(270.0*TINV) 
       rct(6,:) = 4.5e-14*exp(-1260.0*TINV) 
       rct(7,:) = 4.8e-11*exp(250.0*TINV) 
       rct(8,:) = 2.9e-12*exp(-160.0*TINV) 
       rct(9,:) = 7.7e-12*exp(-2100.0*TINV) 
       rct(10,:) = 2.5e-12*exp(-260.0*TINV) 
       rct(11,:) = 2.3e-12*exp(360.0*TINV) 
       rct(12,:) = 7.4e-13*exp(-520.0*TINV) 
       rct(13,:) = 1.03e-13*exp(365.0*TINV)-7.4E-13*exp(-520.0*TINV) 
       rct(14,:) = 3.8e-13*exp(780.0*TINV) 
       rct(15,:) = 5.3e-12*exp(190.0*TINV) 
       rct(16,:) = 2e-12*exp(-2440.0*TINV) 
       rct(17,:) = 6.9e-12*exp(-1000.0*TINV) 
       rct(18,:) = 2.55e-12*exp(380.0*TINV) 
       rct(19,:) = 3.8e-13*exp(900.0*TINV) 
       rct(20,:) = 1.9e-12*exp(190.0*TINV) 
       rct(21,:) = 4.4e-12*exp(365.0*TINV) 
       rct(22,:) = 7.5e-12*exp(290.0*TINV) 
       rct(23,:) = 2e-12*exp(500.0*TINV) 
       rct(24,:) = 2.9e-12*exp(500.0*TINV) 
       rct(25,:) = 5.2e-13*exp(980.0*TINV) 
       rct(26,:) = 2.54e-12*exp(360.0*TINV) 
       rct(27,:) = 1.81875e-13*exp(1300.0*TINV) 
       rct(28,:) = 9.1e-15*exp(-2580.0*TINV) 
       rct(29,:) = 5.5e-15*exp(-1880.0*TINV) 
       rct(30,:) = 1.5132e-13*exp(1300.0*TINV) 
       rct(31,:) = 2.49969e-13*exp(1300.0*TINV) 
       rct(32,:) = 2.05446e-13*exp(1300.0*TINV) 
       rct(33,:) = 1.9e-12*exp(575.0*TINV) 
       rct(34,:) = 1.03e-14*exp(-1995.0*TINV) 
       rct(35,:) = 2.7e-11*exp(390.0*TINV) 
       rct(36,:) = 2.6e-12*exp(610.0*TINV) 
       rct(37,:) = 1.36e-15*exp(-2112.0*TINV) 
       rct(38,:) = 8e-12*exp(380.0*TINV) 
       rct(39,:) = 7.6e-12*exp(180.0*TINV) 
       rct(40,:) = 1.6e-12*exp(305.0*TINV) 
       rct(41,:) = 8.7e-12*exp(290.0*TINV) 
       rct(42,:) = 8.5e-16*exp(-1520.0*TINV) 
       rct(43,:) = 3.15e-12*exp(-450.0*TINV) 
       rct(44,:) = 4.3e-13*exp(1040.0*TINV) 

  end subroutine set_rct_rates
end module  ChemRates_rct_ml
