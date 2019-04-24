! <My_Outputs_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.32>
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
module  My_Outputs_mod
! -----------------------------------------------------------------------
! Allows user to specify which species are output to various
! ascii and binary output files.
!
! Sites  - surface sites,     to sites.out
! Sondes - vertical profiles, to sondes.out
! Hourly - ascii output of selected species, selcted domain
! -----------------------------------------------------------------------

use CheckStop_mod,      only: CheckStop
use ChemDims_mod,       only: NSPEC_ADV, NSPEC_SHL
use ChemSpecs_mod
use ChemGroups_mod,     only: chemgroups
use Debug_module,       only: DEBUG   ! -> DEBUG%COLSRC and POLLEN
use DerivedFields_mod,  only: f_2d,f_3d          ! D2D/D3D houtly output type
use Config_module, only: PPBINV, PPTINV, MasterProc, KMAX_MID,&
                             MY_OUTPUTS, USES, &
                             SELECT_LEVELS_HOURLY!, FREQ_HOURLY
use PhysicalConstants_mod, only: ATWAIR
use OwnDataTypes_mod,   only: Asc2D
use Par_mod,            only: GIMAX,GJMAX,IRUNBEG,JRUNBEG,me
use Pollen_const_mod,   only: pollen_check
use SmallUtils_mod,     only: find_index
use TimeDate_mod,       only: date
use Units_mod,          only: Init_Units,&
                             to_molec_cm3,to_molec_cm2,to_mgSIA,to_ugSIA,&
                             to_ug_ADV,to_ug_C,to_ug_N,to_ug_S,to_number_m3

implicit none

logical, public, parameter :: out_binary = .false.
logical, public, parameter :: Ascii3D_WANTED = .false.

! Site outputs   (used in Sites_mod)
!==============================================================
! Specify the species to be output to the sites.out file
! For met params we have no simple index, so we use characters.
! These must be defined in Sites_mod.f90.

integer, private :: isite              ! To assign arrays, if needed
integer, public, parameter :: &
   NSITES_MAX =        99     & ! Max. no surface sites allowed
  ,FREQ_SITE  =         1     & ! Interval (hrs) between outputs
  ,NADV_SITE  = NSPEC_ADV     & ! No. advected species (1 up to NSPEC_ADV)
  ,NSHL_SITE  = NSPEC_SHL     & ! No. short-lived species
  ,NXTRA_SITE_MISC =    2     & ! No. Misc. met. params  ( e.g. T2, d_2d)
  ,NXTRA_SITE_D2D  =  9+8       ! No. Misc. met. params  ( e.g. T2, d_2d)

integer, public, parameter, dimension(NADV_SITE) :: &
  SITE_ADV = [(isite, isite=1,NADV_SITE)]  ! Everything

integer, public, parameter, dimension(NSHL_SITE) :: &
  SITE_SHL = [(isite, isite=1,NSHL_SITE)]  ! All short-lived species

! Extra parameters - need to be coded in Sites_mod also. So far
! we can choose from hmix, T2, or th (pot. temp.) or d_2d fields.
!  d_2d fields can be accessed from Derived_mod by setting common index
!  "D2D" in SITE_XTRA and the actual field name (as defined in Derived_mod)
!  in SITE_XTRA_CODE (e.g. "D2_PM25 " or "D2_SIA") :

!** IMPORTANT!! Make sure the correspondence between selected for output
!** fields in SITE_XTRA and their names in SITE_XTRA_CODE

character(len=18), public, parameter, dimension(NXTRA_SITE_MISC) :: &
  SITE_XTRA_MISC=[character(len=18):: "th","T2"]

!These variables must have been set in My_Derived for them to be used.
character(len=24), public, parameter, dimension(NXTRA_SITE_D2D) :: &
  SITE_XTRA_D2D=[character(len=24):: &
    "HMIX","PSURF","ws_10m","rh2m",&
    "Emis_mgm2_BioNatC5H8","Emis_mgm2_BioNatBIOTERP",&
    "Emis_mgm2_BioNatNO","Emis_mgm2_nox",&
    'WDEP_PREC',&!'Idirect','Idiffuse','SNratio','SMI_deep',&
    'met2d_uref','met2d_u10','met2d_rh2m', &
    'SMI_deep','met2d_SMI_d','SMI_uppr','met2d_SMI_s',&
!   "SoilWater_deep","EVAP_CF","EVAP_DF",&
!   "EVAP_BF","EVAP_NF","WDEP_PREC",&
!   "RH_GR","CanopyO3_GR","VPD_GR","FstO3_GR",&
!   "RH_IAM_DF","CanopyO3_IAM_DF","VPD_IAM_DF","FstO3_IAM_DF",&
!   "COLUMN_CO_k20","COLUMN_C2H6_k20","COLUMN_HCHO_k20","COLUMN_CH4_k20",
    "COLUMN_NO2_k20"]

!**** Aircraft outputs   (used in Polinat_mod)
!==============================================================
!   Specify the species to be output by Polinat for aircraft flight tracks

integer, public, parameter :: &
  NFLIGHT_MAX =    10         &   ! Max. no sondes allowed
 ,FREQ_FLIGHT =    12         &   ! Interval (hrs) between outputs
 ,NADV_FLIGHT =    1              ! No.  advected species

integer, public, parameter, dimension(NADV_FLIGHT) :: &
  FLIGHT_ADV =  (/ IXADV_O3 /)


!**** Sonde outputs   (used in Sites_mod)
!==============================================================
!     Specify the species to be output to the sondes.out file
!  We typically deal with fewer species for sonde output than
!  surface sites, so we use a different method to specify.
! For met params we have no simple index, so we use characters.
! These must be defined in Sites_mod.f90.

integer, public, parameter :: &
   NSONDES_MAX =    99        &   ! Max. no sondes allowed
  ,NLEVELS_SONDE =  20        &   ! No. k-levels (9 => 0--2500 m)
  ,FREQ_SONDE  =     1        &   ! Interval (hrs) between outputs
  ,NADV_SONDE  =    9         &   ! No.  advected species
  ,NSHL_SONDE  =    1         &   ! No. short-lived species
  ,NXTRA_SONDE =    4             ! No. Misc. met. params

integer, public, parameter, dimension(NADV_SONDE) :: &
  SONDE_ADV = [ IXADV_O3, IXADV_NO2, IXADV_NO, IXADV_PAN,  &
                IXADV_NO3_c, IXADV_NO3_f, IXADV_SO4, IXADV_NH4_f, IXADV_NH3 ]

integer, public, parameter, dimension(NSHL_SONDE) :: &
  SONDE_SHL = [ IXSHL_OH ] ! , IXSHL_OD, IXSHL_OP ]
character(len=10), public, parameter, dimension(NXTRA_SONDE) :: &
  SONDE_XTRA= [character(len=10):: &
   "NOy","z_mid","p_mid","th"]!,"Kz_m2s"]

!   can access d_3d fields through index here, by
!   setting "D3D" above and say D3_XKSIG12 here:

!====================================================================
!**** Hourly outputs   (from hourly_out routine) to print out
!     concentrations  or even met. parameters every hour
!     (or multiple: HOURLY_FREQ) for specified sub-grid.
!     Note: as to met. parameters, only temp2m Th arespecified
!           so far- others need change in hourly_out.f also).

!-------------------------------------------------------------------
!  Possibility of multi-layer output. Specify NLEVELS_HOURLY here
!  and in hr_out defs use either:
!
!  ADVppbv to get surface concentrations (only relevant for layer k=20
!  while gives meaningless  number for higher levels.
!
!  Or BCVppbv to get grid-centre concentrations (relevant for all layers)
!----------------------------------------------------------------

integer, public, save :: &
  nhourly_out=0,    & ! No. output variables
  nlevels_hourly=0, & ! No. output levels
  nmax6_hourly = 0    ! No of 6 hourly max outputs

! Output selected model levels
! If SELECT_LEVELS_HOURLY (ModelConstants_config NML)
! levels_hourly contains which levels to print out
!   20<==>uppermost model level (m01)
!   01<==>lowermost model level (m20)
!   00<==>surface approx. from lowermost model level
integer, public, dimension(:), allocatable :: levels_hourly  ! Set below

type(Asc2D), public, dimension(:), allocatable :: hr_out  ! Set below

!*** wanted binary dates... specify days for which full binary
!    output is wanted. Replaces the hard-coding which was in wrtchem:
integer, public, parameter :: NBDATES = 3
type(date), public, save, dimension(NBDATES) :: wanted_dates_inst

!================================================================

public :: set_output_defs

contains

subroutine set_output_defs
   implicit none

  character(len=144) :: errmsg   ! Local error message
  integer            :: i,j,ash,nuc_conc,nuc_wdep,nuc_ddep,& ! Loop & group indexes
                        rn222,pm25,pm10,nmvoc,gpoll,igrp,idx
  character(len=256) :: name     ! Volcano (vent) name

  real, parameter :: atwC=12.0
  real, parameter :: m_s = 100.0 ! From cm/s to m/s

  ! introduce some integers to make specification of domain simpler
  ! and less error-prone. Numbers can be changed as desired.

  !==============================================================
  ! Conversion to ug/m3
  !   xn_adv(ixadv,ix,iy,k)*roa(ix,iy,k,1)*to_ug_ADV(ixadv)
  ! Conversion to ugX/m3
  !   xn_adv(ixadv,ix,iy,k)*roa(ix,iy,k,1)*to_ug_X(ixadv)
  ! Use "ADVugXX" for ug output (ug/m3, ugC/m3, ugN/m3, ugS/m3)
  !   For ug/m3  output use in combination with to_ug_ADV(ixadv).
  !   For ugX/m3 output use in combination with to_ug_X(ixadv).
  !==============================================================
  call Init_Units()

  !*** Hourly outputs
  !    Note that the hourly output uses **lots** of disc space, so specify
  !    as few as you need and with as small format as possible (cf max value).

  ! ** REMEMBER : ADV species are mixing ratios !!
  ! ** REMEMBER : SHL species are in molecules/cm3, not mixing ratio !!
  ! ** REMEMBER : No spaces in name, except at end !!

  if(MasterProc) write(*,*) "TESTHH INSIDE set_output_defs",MY_OUTPUTS

  select case(MY_OUTPUTS)
  case("EMERGENCY")
    nlevels_hourly = 1+18
    nhourly_out=4    !PM*,AOD & Z
    ash=find_index("ASH",chemgroups(:)%name)
    nuc_conc=0  ! find_index("NUCRACT",chemgroups(:)%name)
    nuc_wdep=0  ! find_index("DDEP_NUCRACT",chemgroups(:)%name)
    nuc_ddep=0  ! find_index("WDEP_NUCRACT",chemgroups(:)%name)
    call CheckStop((ash+nuc_conc+nuc_ddep+nuc_wdep)<1,&
      "set_output_defs: Unknown group 'ASH'/'NUC'")
    if(ash>0)then
      name="none"
      do i=1,size(chemgroups(ash)%specs)
        if(species(chemgroups(ash)%specs(i))%name(1:9)==name)cycle
        name=species(chemgroups(ash)%specs(i))%name(1:9)
        nhourly_out=nhourly_out+2+1  !
        nmax6_hourly=nmax6_hourly+1  !
        if(MasterProc.and.DEBUG%COLSRC)&
          write(*,*)'EMERGENCY: Volcanic Ash, Vent=',name
      enddo
    endif
    if(nuc_conc>0)then
      name="none"
      do i=1,size(chemgroups(nuc_conc)%specs)
        name=species(chemgroups(nuc_conc)%specs(i))%name
        nhourly_out=nhourly_out+1
        if(MasterProc.and.DEBUG%COLSRC)&
          write(*,*)'EMERGENCY: Nuclear accident/explosion, NPP/NUC=',name
      enddo
    endif
    if (.false.) then
    if(nuc_ddep>0)then
      name="none"
      do i=1,size(chemgroups(nuc_ddep)%specs)
        name="DDEP_"//species(chemgroups(nuc_ddep)%specs(i))%name
        nhourly_out=nhourly_out+1
        if(MasterProc.and.DEBUG%COLSRC)&
          write(*,*)'EMERGENCY: Nuclear accident/explosion, NPP/NUC=',name
      enddo
    endif
    if(nuc_wdep>0)then
      name="none"
      do i=1,size(chemgroups(nuc_wdep)%specs)
        name="WDEP_"//species(chemgroups(nuc_wdep)%specs(i))%name
        nhourly_out=nhourly_out+1
        if(MasterProc.and.DEBUG%COLSRC)&
          write(*,*)'EMERGENCY: Nuclear accident/explosion, NPP/NUC=',name
      enddo
    endif
    endif
  case("INVERSION") ! use hourly Derived output in the near future
    nhourly_out=19
    nlevels_hourly = 1
  case("3DPROFILES")
    nhourly_out=2
    nlevels_hourly = 2  ! nb zero is one of levels in this system
    SELECT_LEVELS_HOURLY = .true.
  case("IMPACT2C")
    nhourly_out=4  ! Dave's starting set
    nlevels_hourly = 2  ! nb zero is one of levels, so we have 3m and 45m
  case("TFMM")
    nhourly_out=28
    nlevels_hourly = 1  ! nb zero is *not* one of levels
  case("average_vs_instantaneous")
    nhourly_out = 6
    nlevels_hourly = 1
  case("TRENDS")
    nhourly_out=2
    if(USES%AOD     )nhourly_out=nhourly_out+1
    nlevels_hourly = 1  ! nb zero is *not* one of levels
  case("TRENDS@12UTC")
    nhourly_out=3
    nlevels_hourly = KMAX_MID  ! nb zero is *not* one of levels
!- moved to Hourly/Derived
!-case default
!-  nhourly_out=1
!-  nlevels_hourly = 1  ! nb zero is *not* one of levels
  case default
    nhourly_out=0
  endselect

  if(nhourly_out==0)then
    nlevels_hourly=0
    return
  endif
  if(allocated(hr_out))       deallocate(hr_out)
  if(allocated(levels_hourly))deallocate(levels_hourly)
  allocate(hr_out(nhourly_out),levels_hourly(nlevels_hourly))
  hr_out(:)=Asc2D("none","none",-99,-99,"none",-99.9,-99.9)

  select case(MY_OUTPUTS)
  case("EMERGENCY")
    levels_hourly = (/(i-1,i=1,nlevels_hourly)/)

    pm25=find_index("PMFINE",chemgroups(:)%name) !NB There is no "PM25" group
    pm10=find_index("PM10"  ,chemgroups(:)%name)
!**               name     type     ofmt
!**               ispec    nk sellev? unit conv  max
    j=4;hr_out(1:j) = (/&
      Asc2D("pm25_3km"  ,"BCVugXXgroup",pm25,NLEVELS_HOURLY,"ug/m3",1.0,-999.9),&
      Asc2D("pm10_3km"  ,"BCVugXXgroup",pm10,NLEVELS_HOURLY,"ug/m3",1.0,-999.9),&
      Asc2D("pm_h2o_3km","PMwater"     ,00  ,NLEVELS_HOURLY,"ug/m3",1.0,-999.9),&
!     Asc2D("AOD_550nm" ,"D2D",find_index("AOD_550nm",f_2d(:)%name),&
!            1," ",1.0    ,-9999.9),&
      Asc2D("z"         ,"Z_MID" ,00,NLEVELS_HOURLY,"km",1e-3,-9999.9)/)
    if(ash>0)then
      name="none"
      do i=1,size(chemgroups(ash)%specs)
        if(species(chemgroups(ash)%specs(i))%name(1:9)==name)cycle
        name=species(chemgroups(ash)%specs(i))%name(1:9)
        igrp=find_index(name,chemgroups(:)%name)
        call CheckStop(igrp<1,"set_output_defs: Unknown group '"//name//"'")
        j=j+3;hr_out(j-2:j)=(/&
          Asc2D(trim(name)       ,"BCVugXXgroup",igrp,&
                NLEVELS_HOURLY,"ug/m3",1.0,-999.9),&
          Asc2D(trim(name)//"_col","COLUMNgroup" ,igrp,&
                1,"ug/m2",1.0,-999.9),&
          Asc2D(trim(name)//"_6hour_max","MAX6Hgroup" ,igrp,&
                4,"ug/m3",1.0,-999.9)/)
      enddo
    endif
    if(nuc_conc>0)then
      name="none"
      do i=1,size(chemgroups(nuc_conc)%specs)
        name=species(chemgroups(nuc_conc)%specs(i))%name
        !igrp=find_index(name,chemgroups(:)%name) ! position of NPP/NUC -group
        !call CheckStop(igrp<1,"set_output_defs: Unknown conc group '"//name//"'")
        j=j+1;
        idx= chemgroups(nuc_conc)%specs(i) - NSPEC_SHL ! offset between xn_adv and species
        hr_out(j)= Asc2D(trim(name),"BCVugXX",idx,1,"uBq h/m3",&
                PPBINV/ATWAIR*species(chemgroups(nuc_conc)%specs(i))%molwt,-999.9)
      enddo
    endif
    if (.false.) then
    if(nuc_wdep>0)then
      name="none"
      do i=1,size(chemgroups(nuc_wdep)%specs)
        name=species(chemgroups(nuc_wdep)%specs(i))%name
        igrp=find_index(name,chemgroups(:)%name) ! position of NPP/NUC -group
        call CheckStop(igrp<1,"set_output_defs: Unknown wdep group '"//name//"'")
        j=j+1;
!       hr_out(j)=(/&
!         Asc2D(name        ,"BCVugXXgroup",igrp,&
!               1,"uBqh/m2",1.0,-999.9)&
!                 /)
        hr_out(j)=Asc2D(trim(name)//"_WDEP","D2D_accum",&
              find_index("WDEP_"//trim(name),f_2d(:)%name),1,"",1.0,-999.9)

      enddo
    endif
    if(nuc_ddep>0)then
      name="none"
      do i=1,size(chemgroups(nuc_ddep)%specs)
        name=species(chemgroups(nuc_ddep)%specs(i))%name
        igrp=find_index(name,chemgroups(:)%name) ! position of NPP/NUC -group
        call CheckStop(igrp<1,"set_output_defs: Unknown ddep group '"//name//"'")
        j=j+1;
!       hr_out(j:j)=(/&
!         Asc2D(name        ,"BCVugXXgroup",igrp,&
!               1,"uBqh/m2",1.0,-999.9)&
!               /)
        hr_out(j)=Asc2D(trim(name)//"_DDEP","D2D_accum",&
              find_index("DDEP_"//name,f_2d(:)%name),1,"",1.0,-999.9)
      enddo
    endif
    endif ! if false
  case("INVERSION")
    levels_hourly = [0]
    do j=1,19
      write(name,"(A,I2.2)")"ASH_L",j
      igrp=find_index(trim(name),chemgroups(:)%name)
      call CheckStop(igrp<1,"set_output_defs: Unknown ASH group '"//name//"'")
      if(MasterProc)write(*,*) "AshInv: ",trim(name),igrp
      hr_out(j)=Asc2D(trim(name)//"_col","COLUMNgroup" ,igrp,&
           1,"ug/m2",1e0/3600.,-999.9)  ! 1kg/s --> 1kg/h
    enddo
  case("IMPACT2C") ! Dave's starting set. Uses Out3D to get 3m and 45m concs.
    levels_hourly = (/ (i, i= 0,nlevels_hourly-1) /)  ! -1 will give surfac
    hr_out(:)= (/ &
        Asc2D("o3_3dppb"    ,"Out3D",O3   ,&
             nlevels_hourly,"ppbv", PPBINV,600.0*2.0) &
       ,Asc2D("no_3dppb"   ,"Out3D",&
             NO  ,nlevels_hourly,"ppbv",PPBINV ,600.0*1.91) &
       ,Asc2D("no2_3dppb"   ,"Out3D",&
             NO2  ,nlevels_hourly,"ppbv",PPBINV ,600.0*1.91) &
       ,Asc2D("ho2_3dppt"   ,"Out3D",&  !NOTE ppt
             HO2  ,nlevels_hourly,"pptv",PPBINV ,600.0*1.91) &
       ,Asc2D("HMIX","D2D_inst",find_index("HMIX",f_2d(:)%name), &
            1,"m",1.0,10000.0) &
     /)
  case("3DPROFILES")
   ! nb Out3D uses totals, e.g. O3, not IXADV_O3
   ! Number of definitions must match nhourly_out set above
    levels_hourly = (/ (i, i= 0,nlevels_hourly-1) /)  ! -1 will give surfac
    hr_out(:)= (/ &
        Asc2D("o3_3dppb"    ,"Out3D",O3   ,&
             nlevels_hourly,"ppbv",PPBINV,600.0*2.00) &
       ,Asc2D("no2_3dppb"   ,"Out3D",NO2  ,&
             nlevels_hourly,"ppbv",PPBINV,600.0*1.91) &
!      ,Asc2D("o3_3dug"   ,"Out3D",O3     ,&
!            nlevels_hourly,"ug",to_ug_ADV(IXADV_O3),600.0*2.0) &
    /)

    if(MasterProc ) write(*,*) "TESTHH 3D O3 SET", nlevels_hourly
  case("TFMM")

!**               name     type     ofmt
!**               ispec    ix1 ix2 iy1 iy2 nk sellev? unit conv  max

!!!! It seems easiest to just use many variables as given in the d_2d arrays. Thus
!!!! we search for the name as given there using "find_index" below.
!!!! (As a test I tried both pmfine two ways, one as D2D and the other
!!!!  as ADVugXXXgroup. The results were identical.)

   hr_out(:) = (/  &
      Asc2D("o3_3m", "ADVppbv", IXADV_O3    ,1,"ppbv" ,PPBINV,600.0) &
     ,Asc2D("no3_f"  ,"ADVugXX",IXADV_NO3_F ,&
            1,"ug/m3",to_ug_ADV(IXADV_NO3_F),-999.9) &
     ,Asc2D("no3_c"  ,"ADVugXX",IXADV_NO3_C ,&
            1,"ug/m3",to_ug_ADV(IXADV_NO3_C),-999.9) &
     ,Asc2D("nh4_f"  ,"ADVugXX",IXADV_NH4_F ,&
            1,"ug/m3",to_ug_ADV(IXADV_NH4_F),-999.9) &
     ,Asc2D("so4_f"  ,"ADVugXX",IXADV_SO4   ,&
            1,"ug/m3",to_ug_ADV(IXADV_SO4)  ,-999.9) &
   ! Organics
     ,Asc2D("OM25_3m" ,"D2D_inst",find_index("SURF_ug_OM25_p",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("EC25_3m" ,"D2D_inst",find_index("SURF_ug_ECFINE",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("OM_c_3m" ,"D2D_inst",find_index("SURF_ug_OMCOARSE",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("EC_c_3m" ,"D2D_inst",find_index("SURF_ug_ECCOARSE",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("NatDust_f","D2D_inst",find_index("SURF_ug_Dust_NAT_f",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("NatDust_c","D2D_inst",find_index("SURF_ug_Dust_NAT_c",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("RoadDust_f","D2D_inst",find_index("SURF_ug_Dust_ROAD_f",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("RoadDust_c","D2D_inst",find_index("SURF_ug_Dust_ROAD_c",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("SeaSalt_f","D2D_inst",find_index("SURF_ug_SeaSalt_f",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("SeaSalt_c","D2D_inst",find_index("SURF_ug_SeaSalt_c",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
    ! Sums
     ,Asc2D("PM25_3m" ,"D2D_inst",find_index("SURF_ug_PM25_rh50",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("PMFINE"  ,"D2D_inst",find_index("SURF_ug_PMFINE",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("pm_h2o_3m","D2D_inst",find_index("SURF_PM25water",f_2d(:)%name), &
            1,"ug",1.0,-999.9) &
     ,Asc2D("PM25dry_3m","D2D_inst",find_index("SURF_ug_PM25",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("PM10_3m" ,"D2D_inst",find_index("SURF_ug_PM10_rh50",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("PM10dry_3m","D2D_inst",find_index("SURF_ug_PM10",f_2d(:)%name), &
            1, "ug/m3",1.0,-999.9) &
     ,Asc2D("no_3m" ,"ADVppbv", IXADV_NO,&
            1,"ppbv",PPBINV ,600.0*1.91) &
     ,Asc2D("no2_3m","ADVppbv", IXADV_NO2,&
            1,"ppbv",PPBINV ,600.0*1.91) &
     ,Asc2D("T2_C",   "T2_C"  , 00, &
            1, "degC",1.0   ,100.0) &
     ,Asc2D("ws_10m", "ws_10m", 00, &
            1, "m/s",1.0   ,100.0) &
     ,Asc2D("HMIX","D2D_inst",find_index("HMIX",f_2d(:)%name), &
            1, "m",1.0,10000.0) &
     ,Asc2D("USTAR_NWP","D2D_inst",find_index("USTAR_NWP",f_2d(:)%name), &
            1, "m/s",1.0,-999.9) &
     ,Asc2D("Kz_m2s","D2D_inst",find_index("Kz_m2s",f_2d(:)%name), &
            1, "m2/s",1.0,-999.9) &
     /)
  case("average_vs_instantaneous")
  ! Example of different hourly output types - hourly average and hourly instantaneous value
  ! variables of type "ADVppbv" will be outputted instantaneous and
  ! variables of type "D2D_mean" will be hourly averages.
    hr_out(:) = (/  &
      Asc2D("SURF_ppb_O3_inst","ADVppbv", IXADV_O3,1,"ppbv",PPBINV,600.0) ,&
      Asc2D("SURF_ppb_O3_mean","D2D_mean",find_index("SURF_ppb_O3",f_2d(:)%name), &
           1,"ppbv",1.0   ,600.0), &
      Asc2D("SURF_ug_NO3_F_inst","ADVugXX",IXADV_NO3_F,&
          1,"ug/m3",to_ug_ADV(IXADV_NO3_F),-999.9), &
      Asc2D("SURF_ug_NO3_F_mean","D2D_mean",find_index("SURF_ug_NO3_F",f_2d(:)%name),&
          1,"ug/m3",1.0,-999.9), &
      Asc2D("SURF_ug_PMFINE_inst"  ,"ADVugXX",IXADV_NO3_F  ,&
          1,"ug/m3",to_ug_ADV(IXADV_NO3_F)  ,-999.9), &
      Asc2D("SURF_ug_PMFINE_mean","D2D_mean",find_index("SURF_ug_PMFINE",f_2d(:)%name),&
          1,"ug/m3",1.0,-999.9) &
    /)
  case("TRENDS")
    j=2;hr_out(1:j) = (/  &
      Asc2D("o3_3m"  ,"ADVppbv",IXADV_O3 ,1,"ppbv",PPBINV,600.0), &
      Asc2D("no2_col","COLUMN" ,IXADV_NO2,1,"molec/cm2",to_molec_cm2,-999.9) &
!!    Asc2D("no2_col","COLUMN" ,IXADV_NO2,1,"ug",to_ug_ADV(IXADV_NO2),-999.9) &
    /)
    if(USES%AOD)then
      j=j+1;hr_out(j) = &
      Asc2D("AOD_550nm","D2D_inst",find_index("AOD_550nm",f_2d(:)%name),&
            1," ",1.0,-999.9)
    endif
  case("TRENDS@12UTC")
    levels_hourly = (/(i,i=1,nlevels_hourly)/)
    j=3;hr_out(1:j) = (/  &
!!    Asc2D("so2"    ,"BCVppbv",IXADV_SO2, &
!!          nlevels_hourly,"ppbv",PPBINV,-999.9),&
      Asc2D("so2"    ,"BCVugXX",IXADV_SO2,&
            NLEVELS_HOURLY,"ug/m3",to_ug_ADV(IXADV_SO2),-999.9),&
      Asc2D("so2_3m" ,"ADVppbv",IXADV_SO2,1,"ppbv",PPBINV,-999.9),&
      Asc2D("so2_col","COLUMN" ,IXADV_SO2,1,"molec/cm2",to_molec_cm2,-999.9) &
!!    Asc2D("so2_col","COLUMN" ,IXADV_SO2,1,"ug",to_ug_ADV(IXADV_SO2),-999.9) &
!!    Asc2D("ws_10m", "ws_10m", 00,1,"m/s",1.0,100.0) &
    /)

  case default
    hr_out(:) = (/  &
      Asc2D("o3_3m", "ADVppbv", IXADV_O3,1, "ppbv",PPBINV,600.0) &
    /)
  endselect

  !*** Consistency check:
  ! Was the array set?           R: %name/=none
  ! Was the D2D/Group found?     R: %spec>=0
  ! Is the variable name unique? R: find_index((i)%name,(:)%name)==i
  do i=1,nhourly_out
    if(MasterProc)write(*,*) "TESTHH O3 ATEND", i, nlevels_hourly
    idx = find_index(hr_out(i)%name,hr_out(:)%name)
    if((hr_out(i)%name/="none").and.(hr_out(i)%spec>=0).and.&
       (idx==i)) cycle
    if(MasterProc)then
      write(errmsg,"(A,4(1X,A,'=',I0),2(1X,A,':',A))")&
        "set_output_defs: Failed consistency check",&
        "Hourly",i, "Nhourly",NHOURLY_OUT,&
        "HourlyIndex",idx,"Spec",hr_out(i)%spec,&
        "Name",trim(hr_out(i)%name),"Type",trim(hr_out(i)%type)
      call CheckStop(trim(errmsg))
      write(*,*)trim(errmsg)
      write(*,*)"Warning: reduced ourly output"
    endif
    nhourly_out=i-1
    exit
  enddo

  !*** Wanted dates for instantaneous values output:
  !    specify months,days,hours for which full output is wanted.
  wanted_dates_inst(1) = date(-1,1,1,0,0)
  wanted_dates_inst(2) = date(-1,1,1,3,0)
  wanted_dates_inst(3) = date(-1,1,1,6,0)

endsubroutine set_output_defs
endmodule My_Outputs_mod
