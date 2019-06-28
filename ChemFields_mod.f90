! <ChemFields_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module ChemFields_mod
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
use AeroConstants_mod,  only: NSAREA_DEF       ! =>  z dimension
use AllocInits,     only: AllocInit
use ChemDims_mod,   only: NSPEC_ADV, NSPEC_SHL, NSPEC_TOT, & ! => No. species 
                          NCHEMRATES, NPHOTOLRATES 
use ChemSpecs_mod,  only: FIRST_SEMIVOL, LAST_SEMIVOL    ! -999 unless SOA used
use Config_module,  only: KMAX_MID, KCHEMTOP &           ! =>  z dimension
                         ,MasterProc &
                         ,NATBIO                         ! for Nrcbio
use DefPhotolysis_mod, only: NRCPHOTextended ! cludge
use NumberConstants,only: UNDEF_R
use Par_mod,        only: LIMAX,LJMAX                    ! => x, y dimensions
use ZchemData_mod    ! rct, h2o, ..
implicit none
private

!-------- this snipppet was from older GenSpec_bgn_mod. ------
  ! PRETTY MUCH FAKED FOR NOW. CAN BE DELETED SOON IN HOPE!
  !+ Defines indices and NSPEC for bgn : Background species

   ! Species which can be specified simply for each column, e.g.
   ! as function of local meteorology or zenith angle
   !   o2, m,  and for MADE-like, oh, ch3coo2

   integer, public, parameter ::  NSPEC_BGN = 0 ! No. 3D bgn species
   integer, public, parameter ::  NSPEC_COL = 0 ! total no. prescribed specs

    !/ define xn_2d_bgn here.
     real, public, save, allocatable, dimension(:,:) :: xn_2d_bgn

!-------- end of this snipppet from older GenSpec_bgn_mod. ------

    !----------------- basic chemical fields ----------------------------------!
    !  Here we declare and initialise to zero the chemical fields used in the  !
    !  model, as well as cfac (converts from 50m to 1m/3m output)         ! 
    !---------------------------------------------------------------------!


    !March 2017: moved from Solver
  real, save, public, dimension(NSPEC_TOT):: &
     x, xold ,xnew  ! Work arrays [molecules/cm3]

  real, save, allocatable, public :: &
     xn_adv(:,:,:,:)  &
    ,xn_shl(:,:,:,:)  &
    ,xn_bgn(:,:,:,:) &
    ,PM25_water(:,:,:) &  !3D PM water
    ,PM25_water_rh50(:,:) &   !gravimetric PM water
    ,Gerber_water(:,:,:)   !3D PM water from GERBER

  real, save, allocatable, public :: &
    SurfArea_um2cm3(:,:,:)  !2D  aerosol surface area, um2/cm3  (n,i,j)

  real, public, save, allocatable:: Fgas3d (:,:,:,:)  ! for SOA

  real, public, save, allocatable :: AOD(:,:,:,:),Extin_coeff(:,:,:,:,:)

  real, save, allocatable, public :: &
     cfac   (:,:,:)   

  real, save, allocatable, public :: &
     so2nh3_24hr(:,:)!hf CoDep

  real, save, allocatable, public :: &
     Grid_snow(:,:) !snow_flag fraction in grid

  real, save, public :: cell_tinv  ! 1/temp,  tmp location

  public ::alloc_ChemFields

contains

  subroutine alloc_ChemFields()

    allocate(xn_adv(NSPEC_ADV,LIMAX,LJMAX,KMAX_MID))
    xn_adv=0.0
    allocate(xn_shl(NSPEC_SHL,LIMAX,LJMAX,KMAX_MID))
    xn_shl=0.0
    allocate(xn_bgn(NSPEC_BGN,LIMAX,LJMAX,KMAX_MID))
    xn_bgn=0.0
    allocate(PM25_water(LIMAX,LJMAX,KMAX_MID))
    PM25_water=0.0
    allocate(PM25_water_rh50(LIMAX,LJMAX))
    PM25_water_rh50=0.0
    allocate(cfac(NSPEC_ADV,LIMAX,LJMAX))
    cfac=1.0
    allocate(so2nh3_24hr(LIMAX,LJMAX))
    so2nh3_24hr=0.0
    allocate(Grid_snow(LIMAX,LJMAX))
    Grid_snow=0.0
    allocate(xn_2d_bgn(1,KCHEMTOP:KMAX_MID))

    allocate(xn_2d(NSPEC_TOT,KCHEMTOP:KMAX_MID))
    xn_2d = 0.0

    allocate(rct(NCHEMRATES,KCHEMTOP:KMAX_MID))
    rct = 0.0

    allocate(rcphot(NRCPHOTextended,KCHEMTOP:KMAX_MID))
    rcphot = 0.0

    allocate(rcbio(NATBIO%Nrcbio,KCHEMTOP:KMAX_MID))
    rcbio = 0.0

    allocate(Fgas(NSPEC_TOT,KCHEMTOP:KMAX_MID),Fpart(NSPEC_TOT,KCHEMTOP:KMAX_MID))
    Fgas  = 1.0! Fraction as gas-phase
    Fpart = 0.0

  ! Fgas3D is only defined for the semivolatile VOC/SOA stuff
  ! We need to assume something on 1st time-step though:

    if(FIRST_SEMIVOL>0)then !FSOA
      allocate(Fgas3d(FIRST_SEMIVOL:LAST_SEMIVOL,LIMAX,LJMAX,KCHEMTOP:KMAX_MID))
      Fgas3d = 1.0
    end if

    allocate(rcemis(NSPEC_SHL+1:NSPEC_TOT,KCHEMTOP:KMAX_MID))
    allocate(deltaZcm(KCHEMTOP:KMAX_MID))
    rcemis = 0.0
    allocate(rh(KCHEMTOP:KMAX_MID),M(KCHEMTOP:KMAX_MID),o2(KCHEMTOP:KMAX_MID))
    allocate(n2(KCHEMTOP:KMAX_MID),h2o(KCHEMTOP:KMAX_MID),temp(KCHEMTOP:KMAX_MID))
    allocate(tinv(KCHEMTOP:KMAX_MID),pp(KCHEMTOP:KMAX_MID))
    allocate(itemp(KCHEMTOP:KMAX_MID))
    CHEMSIZE = KMAX_MID-KCHEMTOP+1

   ! Surface area and water

    allocate(surfarea_um2cm3(NSAREA_DEF,LIMAX,LJMAX))
    SurfArea_um2cm3=0.0
    allocate(Gerber_water(LIMAX,LJMAX,KMAX_MID))
    Gerber_water=0.0

   ! wet DpgN and defaults from dry values
    allocate(DpgNw(NSAREA_DEF, KCHEMTOP:KMAX_MID))

    allocate(S_m2m3(NSAREA_DEF, KCHEMTOP:KMAX_MID)) ! GERBER
    S_m2m3=0.0

    ! Mol speeds
    allocate(cn2o5(KCHEMTOP:KMAX_MID),chno3(KCHEMTOP:KMAX_MID),&
              cho2(KCHEMTOP:KMAX_MID),co3(KCHEMTOP:KMAX_MID))
    cn2o5=UNDEF_R
    chno3=UNDEF_R
    cho2= UNDEF_R
    co3=  UNDEF_R
    allocate(aero_fom(KCHEMTOP:KMAX_MID),aero_fdust(KCHEMTOP:KMAX_MID),&
             aero_fbc(KCHEMTOP:KMAX_MID),aero_fss  (KCHEMTOP:KMAX_MID))
    aero_fom    = UNDEF_R
    aero_fbc    = UNDEF_R
    aero_fss    = UNDEF_R
    aero_fdust  = UNDEF_R

    allocate(gamN2O5(KCHEMTOP:KMAX_MID)) ! kHet  for output
    allocate(cNO2(KCHEMTOP:KMAX_MID),cNO3(KCHEMTOP:KMAX_MID)) ! kHet test
  

  end subroutine alloc_ChemFields


!_____________________________________________________________________________
endmodule ChemFields_mod
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
