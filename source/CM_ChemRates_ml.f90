!>_________________________________________________________<

module ChemRates_rct_ml
!-----------------------------------------------------------

 
  use ChemFunctions_ml       ! => kaero, RiemerN2O5

  use AeroFunctions     ! => UpdakeRate, cMolSpeed
  use Setup_1dfields_ml ! => tinv, h2o, m, Fgas
  use Setup_1dfields_ml, m=> amk
  use ChemSpecs_tot_ml  ! => PINALD, .... for FgasJ08
  use ModelConstants_ml, only: KMAX_MID,KCHEMTOP,DebugCell,DEBUG,AERO
implicit none
private

  !+ Tabulates Rate-coefficients - temperature dependant

    public :: set_rct_rates

    integer, parameter, public :: NRCT = 9   !! No. coefficients

    real, allocatable, save, public, dimension(:,:) :: rct

  contains
  !------------------------------------
  subroutine set_rct_rates() 
     logical,save::first_call=.true.
     if(first_call)then
       allocate(rct(NRCT,KCHEMTOP:KMAX_MID))
       rct=0.0
     endif
       rct(1,:) = 6.0e-34*M*O2*(TEMP/300.0)**2.6 
       rct(2,:) = 2.2e-10*H2O 
       rct(3,:) = 1.4e-12*exp(-1310.0*TINV) 
       rct(4,:) = 2.03e-17*(TEMP**2)*exp(78.0*TINV) 
       rct(5,:) = 4.4e-12*exp(365.0*TINV) 
       rct(6,:) = 3.6e-12*exp(270.0*TINV) 
       rct(7,:) = 2.54e-12*exp(360.0*TINV) 
       rct(8,:) = 7.5e-12*exp(290.0*TINV) 
       rct(9,:) = 1.95e+16*exp(-13543.0*TINV) 
       first_call=.false.
  end subroutine set_rct_rates
end module  ChemRates_rct_ml
