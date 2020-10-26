! <Pollen_const_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.36>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2020 met.no
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
module Pollen_const_mod
!-----------------------------------------------------------------------!
! Birch pollen emission calculation based on
! M. Sofiev et al. 2006, doi:10.1007/s00484-006-0027-x
!
! Pollen emission based upon meteorology paparameters, and heatsum.
! Pollen particles are assumed of 22 um diameter and 800 kg/m3 density.
!-----------------------------------------------------------------------!
use PhysicalConstants_mod, only: PI,ATWAIR,AVOG
use Config_module,         only: USES,MasterProc
use Debug_module,          only: DEBUG
use CheckStop_mod,         only: CheckStop
use ChemDims_mod,          only: NSPEC_ADV,NSPEC_SHL
use ChemSpecs_mod,         only: species
use ChemGroups_mod,        only: chemgroups
use SmallUtils_mod,        only: find_index
implicit none
public

real, parameter  :: &
  T_cutoff_birch = 273.2+3.5, & ! Cut-off temperature [K] birch
  T_cutoff_olive = 273.2,     & ! Cut-off temperature [K] olive
  dH_d_birch     =  50,       & ! Flowering period [degree days] birch
  dH_d_olive     = 275,       & ! Flowering period [degree days] olive
  dH_birch=dH_d_birch*24*3600,& ! Flowering period [degree seconds] birch
! dH_olive=dH_d_olive*24*3600,& ! Flowering period [degree seconds] olive
  PREC_MIN = 0.0,             & ! Min cut-off precipitation [mm/h]
  PREC_MAX = 0.5,             & ! Max cut-off precipitation [mm/h]
  N_TOT_birch = 1.0e9,        & ! Available pollen [grains/m2] birch
  N_TOT_olive = 3.9e8,        & ! Available pollen [grains/m2] olive
  N_TOT_rweed = 1.7e7,        & ! Available pollen [grains/m2] ragweed
  N_TOT_grass = 2.0e7,        & ! Available pollen [grains/m2] grass
  RH_LOW   = 0.50,            & ! Min cut-off relative humidity [1]
  RH_HIGH  = 0.80,            & ! Max cut-off relative humidity [1]
  PROB_IN_birch  = 0.2,       & ! Probability for flowering to start
  PROB_OUT_birch = 0.2,       & ! Probability for flowering to end
  PROB_IN_olive  = 0.1,       & ! Probability for flowering to start
  PROB_OUT_olive = 0.1,       & ! Probability for flowering to end
                                ! (could be assumed to be larger than PROB_IN)
  uncert_day_grass = 7,       &
  uncert_tot_grass = 0.2,     & ! end uncertainty for linear releases
  D_POLL_birch = 22.0,        & ! Pollen grain diameter [um] birch
  D_POLL_olive = 28.0,        & ! Pollen grain diameter [um] olive
  D_POLL_rweed = 18.0,        & ! Pollen grain diameter [um] grass
  D_POLL_grass = 32.0,        & ! Pollen grain diameter [um] grass
  POLL_DENS    = 800e3          ! Pollen density [g/m3]

real, parameter ::            &
  temp_min_rweed   = 274.05,  &  !  0.9C (loTemp)
  temp_opt_rweed   = 304.85,  &  ! 31.7C (optTemp)
  temp_max_rweed   = 313.15,  &  ! 40.0C (hiTemp)
  photoperiod_rweed=  14.5,   &  ! date%hour
  HS_startday_rweed=  79.0,   &  ! 20 March, spring equinox
  startThr_HS_rweed=  25.0,   &  ! Deen et al 1998 (StartHSThr)
  uncert_HS_rweed  =   0.1,   &  ! 10.0% (fUncertainty_HS_relative_start)
  TempThr_rweed    = 273.15,  &  !   0C (Deen et al 1998)
  DayTempThr_rweed = 280.65,  &  ! 7.5C (Deen et al 1998)
  uncert_day_rweed =  30.0,   &  ! 30 days (fUncertainty_CD_days_start)
  EndCDThr_rweed   = 265.0       ! 97.5% == 2sigma; 22 sept (autumn equinox)

! pollen arrays indexing, order must match with POLLEN_GROUP: birch,olive,rweed,grass
character(len=*), parameter :: &
  BIRCH = "POLLEN_BIRCH",&
  OLIVE = "POLLEN_OLIVE",&
  RWEED = "POLLEN_RWEED",&
  GRASS = "POLLEN_GRASS",&
  POLLEN_GROUP(4)=[BIRCH,OLIVE,RWEED,GRASS]
integer, parameter :: &
  iBIRCH=1,iOLIVE=2,iRWEED=3,iGRASS=4,POLLEN_NUM=size(POLLEN_GROUP)
real, parameter  :: &
  N_TOT(POLLEN_NUM)=[N_TOT_birch,N_TOT_olive,N_TOT_rweed,N_TOT_grass],&
  T_CUTOFF(iBIRCH:iOLIVE)=[T_cutoff_birch,T_cutoff_olive],&
  PROB_IN(iBIRCH:iOLIVE)=[PROB_IN_birch,PROB_IN_olive],&
  PROB_OUT(iBIRCH:iOLIVE)=[PROB_OUT_birch,PROB_OUT_olive]

real, parameter  :: &
  D_POLL(POLLEN_NUM)=[D_POLL_birch,D_POLL_olive,D_POLL_rweed,D_POLL_grass], & ! pollen diameter
  grain_wt(POLLEN_NUM) = POLL_DENS*PI*(D_POLL*1e-6)**3/6.0       ! 1 grain weight [g]
! weight 1 grain [ug], 1 mol of grains (AVOG*grain_wt) [Tonne=1e3 kg]
! BIRCH: 4.460e-3, 2686e6
! OLIVE: 9.195e-3, 5538e6
! RWEED: 2.443e-3
! GRASS: 13.73e-3, 8267e6

private :: N_TOT_birch,N_TOT_olive,N_TOT_rweed,N_TOT_grass,&
           T_cutoff_birch,T_cutoff_olive,&
           PROB_IN_birch,PROB_IN_olive,PROB_OUT_birch,PROB_OUT_olive,&
           D_POLL_birch,D_POLL_olive,D_POLL_rweed,D_POLL_grass

contains
subroutine pollen_check(igrp,uconv_adv)
  integer, intent(inout), optional :: igrp
  real, dimension(NSPEC_ADV), intent(inout), optional :: uconv_adv
  integer :: poll,g
  logical,save :: first_call=.true.
  poll=find_index("POLLEN",chemgroups(:)%name)
  if(present(igrp))igrp=poll
  if(.not.first_call)return
  first_call=.false.
  call CheckStop(USES%POLLEN.and.(poll<1),&
    "USES%POLLEN on model compiled without pollen")
  call CheckStop(DEBUG%POLLEN.and..not.USES%POLLEN,&
    "DEBUG%POLLEN on run without USES%POLLEN")
  if(.not.USES%POLLEN)return
  call CheckStop(size(chemgroups(poll)%specs),POLLEN_NUM,&
    "pollen_check: Inconsistent POLLEN group size")
  call CheckStop(any(species(chemgroups(poll)%specs)%name/=POLLEN_GROUP),&
    "pollen_check: Inconsistent POLLEN group species")
  do g=1,POLLEN_NUM
    select case(g)
    case(iBIRCH)
      call CheckStop(POLLEN_GROUP(g),BIRCH,&
        "pollen_check: Inconsistent POLLEN group order, "//POLLEN_GROUP(g))
      call CheckStop(N_TOT(g)/=N_TOT_birch,&
        "pollen_check: Inconsistent POLLEN group total, "//POLLEN_GROUP(g))
      call CheckStop(D_POLL(g)/=D_POLL_birch,&
        "pollen_check: Inconsistent POLLEN group diameter, "//POLLEN_GROUP(g))
      call CheckStop(T_CUTOFF(g)/=T_cutoff_birch,&
        "pollen_check: Inconsistent POLLEN group T_cutoff, "//POLLEN_GROUP(g))
      call CheckStop(PROB_IN(g)/=PROB_IN_birch,&
        "pollen_check: Inconsistent POLLEN group PROB_IN, "//POLLEN_GROUP(g))
      call CheckStop(PROB_OUT(g)/=PROB_OUT_birch,&
        "pollen_check: Inconsistent POLLEN group PROB_OUT, "//POLLEN_GROUP(g))
    case(iOLIVE)
      call CheckStop(POLLEN_GROUP(g),OLIVE,&
        "pollen_check: Inconsistent POLLEN group order, "//POLLEN_GROUP(g))
      call CheckStop(N_TOT(g)/=N_TOT_olive,&
        "pollen_check: Inconsistent POLLEN group total, "//POLLEN_GROUP(g))
      call CheckStop(D_POLL(g)/=D_POLL_olive,&
        "pollen_check: Inconsistent POLLEN group diameter, "//POLLEN_GROUP(g))
      call CheckStop(T_CUTOFF(g)/=T_cutoff_olive,&
        "pollen_check: Inconsistent POLLEN group T_cutoff, "//POLLEN_GROUP(g))
      call CheckStop(PROB_IN(g)/=PROB_IN_olive,&
        "pollen_check: Inconsistent POLLEN group PROB_IN, "//POLLEN_GROUP(g))
      call CheckStop(PROB_OUT(g)/=PROB_OUT_olive,&
        "pollen_check: Inconsistent POLLEN group PROB_OUT, "//POLLEN_GROUP(g))
    case(iRWEED)
      call CheckStop(POLLEN_GROUP(g),RWEED,&
        "pollen_check: Inconsistent POLLEN group order, "//POLLEN_GROUP(g))
      call CheckStop(N_TOT(g)/=N_TOT_rweed,&
        "pollen_check: Inconsistent POLLEN group total, "//POLLEN_GROUP(g))
      call CheckStop(D_POLL(g)/=D_POLL_rweed,&
        "pollen_check: Inconsistent POLLEN group diameter, "//POLLEN_GROUP(g))
    case(iGRASS)
      call CheckStop(POLLEN_GROUP(g),GRASS,&
        "pollen_check: Inconsistent POLLEN group order, "//POLLEN_GROUP(g))
      call CheckStop(N_TOT(g)/=N_TOT_grass,&
        "pollen_check: Inconsistent POLLEN group total, "//POLLEN_GROUP(g))
      call CheckStop(D_POLL(g)/=D_POLL_grass,&
        "pollen_check: Inconsistent POLLEN group diameter, "//POLLEN_GROUP(g))
    end select
  end do
  if(present(uconv_adv))then
    uconv_adv(chemgroups(poll)%specs-NSPEC_SHL)=&
      uconv_adv(chemgroups(poll)%specs-NSPEC_SHL)/grain_wt
    if(DEBUG%POLLEN.and.MasterProc) &
      write(*,*)'POLLEN uconv_adv',uconv_adv(chemgroups(poll)%specs-NSPEC_SHL)
  end if
end subroutine pollen_check
end module Pollen_const_mod
