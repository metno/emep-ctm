! <My_Outputs_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module  My_Outputs_ml
! -----------------------------------------------------------------------
! Allows user to specify which species are output to various
! ascii and binary output files.
!
! Sites  - surface sites,     to sites.out
! Sondes - vertical profiles, to sondes.out
! Hourly - ascii output of selected species, selcted domain
! -----------------------------------------------------------------------

use CheckStop_ml,      only: CheckStop
use ChemSpecs_tot_ml
use ChemSpecs_adv_ml
use ChemSpecs_shl_ml
use ChemChemicals_ml,  only: species,species_adv

use ChemGroups_ml,     only: chemgroups
use DerivedFields_ml,  only: f_2d               ! D2D houtly output type
use ModelConstants_ml, only: PPBINV, PPTINV, ATWAIR, atwS, atwN, MasterProc, &
                             EXP_NAME, FORECAST, USE_EMERGENCY,DEBUG_EMERGENCY
use OwnDataTypes_ml,   only: Asc2D
use Par_ml,            only: GIMAX,GJMAX,IRUNBEG,JRUNBEG,me
use SmallUtils_ml,     only: find_index
use TimeDate_ml,       only: date
use Units_ml,          only: Init_Units,&
                             to_molec_cm3,to_molec_cm2,to_mgSIA,to_ugSIA,&
                             to_ug_ADV,to_ug_C,to_ug_N,to_ug_S

implicit none

logical, public, parameter :: out_binary = .false.
logical, public, parameter :: Ascii3D_WANTED = .false.

! Site outputs   (used in Sites_ml)
!==============================================================
! Specify the species to be output to the sites.out file
! For met params we have no simple index, so we use characters.
! These must be defined in Sites_ml.f90.

integer, private :: isite              ! To assign arrays, if needed
integer, public, parameter :: &
   NSITES_MAX =       99      & ! Max. no surface sites allowed
  ,FREQ_SITE  =        1      & ! Interval (hrs) between outputs
  ,NADV_SITE  = NSPEC_ADV  & ! No. advected species (1 up to NSPEC_ADV)
  ,NSHL_SITE  = NSPEC_SHL  & ! No. short-lived species
  ,NXTRA_SITE_MISC =    2     & ! No. Misc. met. params  ( e.g. T2, d_2d)
  ,NXTRA_SITE_D2D  =    9       ! No. Misc. met. params  ( e.g. T2, d_2d)

integer, public, parameter, dimension(NADV_SITE) :: &
  SITE_ADV =  (/ (isite, isite=1,NADV_SITE) /)  ! Everything

integer, public, parameter, dimension(NSHL_SITE) :: &
  SITE_SHL =  (/ (isite, isite=1,NSHL_SITE) /)  ! All short-lived species

! Extra parameters - need to be coded in Sites_ml also. So far
! we can choose from hmix, T2, or th (pot. temp.) or d_2d fields.
!  d_2d fields can be accessed from Derived_ml by setting common index
!  "D2D" in SITE_XTRA and the actual field name (as defined in Derived_ml)
!  in SITE_XTRA_CODE (e.g. "D2_PM25 " or "D2_SIA") :

!** IMPORTANT!! Make sure the correspondence between selected for output
!** fields in SITE_XTRA and their names in SITE_XTRA_CODE

character(len=18), public, parameter, dimension(NXTRA_SITE_MISC) :: &
  SITE_XTRA_MISC=(/"th   ","T2   "/)

!These variables must have been set in My_Derived for them to be used.
character(len=24), public, parameter, dimension(NXTRA_SITE_D2D) :: &
  SITE_XTRA_D2D= (/ &
    "HMIX                   ",&
    "PSURF                  ", &
    "ws_10m                 ", &
    "rh2m                   ", &
    "Emis_mgm2_BioNatC5H8   ", &
    "Emis_mgm2_BioNatAPINENE", &
    "Emis_mgm2_BioNatNO     ",&
    "Emis_mgm2_nox          ",&
!   "SoilWater_deep ","EVAP_CF        ","EVAP_DF        ", &
!   "EVAP_BF        ","EVAP_NF        ","WDEP_PREC      ", &
!   "RH_GR          ","CanopyO3_GR    ","VPD_GR         ","FstO3_GR       ", &
!   "RH_IAM_DF      ","CanopyO3_IAM_DF","VPD_IAM_DF     ","FstO3_IAM_DF   ", &
!   "COLUMN_CO_k20  ","COLUMN_C2H6_k20","COLUMN_HCHO_k20","COLUMN_CH4_k20 ",
    "COLUMN_NO2_k20         " /)

!/*** Aircraft outputs   (used in Polinat_ml)
!==============================================================
!   Specify the species to be output by Polinat for aircraft flight tracks

integer, public, parameter :: &
  NFLIGHT_MAX =    10         &   ! Max. no sondes allowed
 ,FREQ_FLIGHT =    12         &   ! Interval (hrs) between outputs
 ,NADV_FLIGHT =    1              ! No.  advected species

integer, public, parameter, dimension(NADV_FLIGHT) :: &
  FLIGHT_ADV =  (/ IXADV_O3 /)


!/*** Sonde outputs   (used in Sites_ml)
!==============================================================
!     Specify the species to be output to the sondes.out file
!  We typically deal with fewer species for sonde output than
!  surface sites, so we use a different method to specify.
! For met params we have no simple index, so we use characters.
! These must be defined in Sites_ml.f90.

integer, public, parameter :: &
   NSONDES_MAX =    99        &   ! Max. no sondes allowed
  ,NLEVELS_SONDE =  20        &   ! No. k-levels (9 => 0--2500 m)
  ,FREQ_SONDE  =     1        &   ! Interval (hrs) between outputs
  ,NADV_SONDE  =    9         &   ! No.  advected species
  ,NSHL_SONDE  =    3         &   ! No. short-lived species
  ,NXTRA_SONDE =    4             ! No. Misc. met. params

integer, public, parameter, dimension(NADV_SONDE) :: &
   SONDE_ADV =  (/ IXADV_O3, IXADV_NO2, IXADV_NO, IXADV_PAN,  &
                IXADV_NO3_c, IXADV_NO3_f, IXADV_SO4, IXADV_NH4_f, IXADV_NH3/)

integer, public, parameter, dimension(NSHL_SONDE) :: &
  SONDE_SHL =  (/ IXSHL_OH, IXSHL_OD, IXSHL_OP /)
character(len=10), public, parameter, dimension(NXTRA_SONDE) :: &
  SONDE_XTRA=  (/ "NOy   ", "z_mid ", "p_mid ", "th    " /) !, "Kz_m2s" /)


!   can access d_3d fields through index here, by
!   setting "D3D" above and say D3_XKSIG12 here:

!====================================================================
!/*** Hourly outputs   (from hourly_out routine) to print out
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

!TESTHH integer, public            :: NHOURLY_OUT =  9 ! No. outputs
!TESTHH integer, public, parameter :: NLEVELS_HOURLY = 4 ! No. outputs
integer, public, save      :: nhourly_out=0    ! No. outputs
integer, public, save      :: nlevels_hourly=0 ! No. outputs
integer, public, parameter :: FREQ_HOURLY = 1  ! 1 hours between outputs

! Output selected model levels
logical, public, parameter ::  &
  SELECT_LEVELS_HOURLY = .false..or.FORECAST.or.(EXP_NAME=="3DPROFILES")
! Decide which levels to print out
! 20<==>uppermost model level (m01)
! 01<==>lowermost model level (m20)
! 00<==>surface approx. from lowermost model level
! 00 and 01 can be both printed out,
! but it might create loads of missing values...
!TESTHH integer, public, parameter, dimension(NLEVELS_HOURLY) :: &
!TESTHH   LEVELS_HOURLY = (/0,4,6,10/)
integer, public, dimension(:), allocatable :: levels_hourly  ! Set below

type(Asc2D), public, dimension(:), allocatable :: hr_out  ! Set below

!/** wanted binary dates... specify days for which full binary
!    output is wanted. Replaces the hard-coding which was in wrtchem:
integer, public, parameter :: NBDATES = 3
type(date), public, save, dimension(NBDATES) :: wanted_dates_inst

!================================================================

public :: set_output_defs

contains

subroutine set_output_defs
   implicit none

  character(len=144) :: errmsg   ! Local error message
  integer            :: i,j,ash,rn222,pm25,pm10,ivent ! Loop & ash group indexes
  character(len=9)   :: vent     ! Volcano (vent) name

  real, parameter :: atwC=12.0
  real, parameter :: m_s = 100.0 ! From cm/s to m/s

  ! introduce some integers to make specification of domain simpler
  ! and less error-prone. Numbers can be changed as desired.

 !integer, save :: ix1 = 36, ix2 = 167, iy1=12, iy2 =  122  ! EMEP
 ! integer, save :: ix1 = 65, ix2 = 167, iy1=12, iy2 =  122  ! restricted EMEP
 integer, save :: ix1, ix2, iy1, iy2

  ! WARNING: If the specification of the subdomain is different for
  !            different components (ix1=125 for ozone and ix1=98 for
  !            NH4 for example) , the variables i_EMEP, j_EMEP
  !            latitude and longitude in NetCDF output will be
  !            wrong.


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

  !/** Hourly outputs
  !    Note that the hourly output uses **lots** of disc space, so specify
  !    as few as you need and with as small format as possible (cf max value).

  ! ** REMEMBER : ADV species are mixing ratios !!
  ! ** REMEMBER : SHL species are in molecules/cm3, not mixing ratio !!
  ! ** REMEMBER : No spaces in name, except at end !!

!default
  ix1=IRUNBEG
  ix2=IRUNBEG+GIMAX-1
  iy1=JRUNBEG
  iy2=JRUNBEG+GJMAX-1 

  if(MasterProc) write(*,*) "TESTHH INSIDE set_output_defs",EXP_NAME

  select case(EXP_NAME)
  case("EMERGENCY")
    nlevels_hourly = 1+18
    nhourly_out=4+1    !PM*,AOD (&Z)
    ash=find_index("ASH",chemgroups(:)%name)
    call CheckStop(ash<1,"set_output_defs: Unknown group 'ASH'")
    vent="none"
    do i=1,size(chemgroups(ash)%ptr)
      if(species(chemgroups(ash)%ptr(i))%name(1:9)==vent)cycle
      vent=species(chemgroups(ash)%ptr(i))%name(1:9)
      nhourly_out=nhourly_out+2
      if(MasterProc.and.DEBUG_EMERGENCY)&
        write(*,*)'EMERGENCY: Volcanic Ash, Vent=',vent
    enddo
  case("FORECAST")
    nhourly_out=11
    nlevels_hourly = 4
  case("EVA2010")
    nhourly_out=4
    nlevels_hourly = 1
  case("3DPROFILES")
    nhourly_out=2
    nlevels_hourly = 10  ! nb zero is one of levels in this system
  case("TFMM")
    nhourly_out=28
    nlevels_hourly = 1  ! nb zero is *not* one of levels
  case default
    nhourly_out=1
    nlevels_hourly = 1  ! nb zero is *not* one of levels
  endselect

  if(allocated(hr_out))       deallocate(hr_out)
  if(allocated(levels_hourly))deallocate(levels_hourly)
  allocate(hr_out(nhourly_out),levels_hourly(nlevels_hourly))
  hr_out(:)=Asc2D("none","none",-99,-99,-99,-99,-99,-99,"none",-99.9,-99.9)

  select case(EXP_NAME)
  case("EMERGENCY")
!   ix1=IRUNBEG;ix2=IRUNBEG+GIMAX-1
!   iy1=JRUNBEG;iy2=JRUNBEG+GJMAX-1
    levels_hourly = (/(i-1,i=1,nlevels_hourly)/)

    pm25 =find_index("PMFINE",chemgroups(:)%name) !NB There is no "PM25" group
    pm10 =find_index("PM10"  ,chemgroups(:)%name)
!**               name     type     ofmt
!**               ispec    ix1 ix2 iy1 iy2 nk sellev? unit conv  max
    j=5;hr_out(:j) = (/&
      Asc2D("pm25_3km"  ,"BCVugXXgroup",pm25,&
             ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",1.0,-999.9),&
      Asc2D("pm10_3km"  ,"BCVugXXgroup",pm10,&
             ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",1.0,-999.9),&
      Asc2D("pm_h2o_3km","PMwater",00         ,&
             ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",1.0,-999.9),&
      Asc2D("AOD_550nm" ,"AOD"   ,00         ,&
             ix1,ix2,iy1,iy2,1," ",1.0    ,-9999.9),&
      Asc2D("z"         ,"Z_MID" ,00         ,&
             ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"km",1e-3,-9999.9)/)
    vent="none"
    do i=1,size(chemgroups(ash)%ptr)
      if(species(chemgroups(ash)%ptr(i))%name(1:9)==vent)cycle
      vent=species(chemgroups(ash)%ptr(i))%name(1:9)
      ivent=find_index(vent,chemgroups(:)%name)
      call CheckStop(ivent<1,"set_output_defs: Unknown group '"//vent//"'")
      j=j+2;hr_out(j-1:j)=(/&
        Asc2D(vent        ,"BCVugXXgroup",ivent,&
              ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",1.0,-999.9),&
        Asc2D(vent//"_col","COLUMNgroup" ,ivent,&
              ix1,ix2,iy1,iy2,1,"ug",1.0,-999.9)/)
    enddo
  case("FORECAST")
!   ix1=IRUNBEG;ix2=IRUNBEG+GIMAX-1
!   iy1=JRUNBEG;iy2=JRUNBEG+GJMAX-1
    levels_hourly = (/0,4,6,10/)
    rn222=find_index("RN222",species_adv(:)%name)
    pm25 =find_index("PMFINE",chemgroups(:)%name) !NB There is no "PM25" group
    pm10 =find_index("PM10"  ,chemgroups(:)%name)  
!**               name     type     ofmt
!**               ispec    ix1 ix2 iy1 iy2 nk sellev? unit conv  max
    hr_out(:) = (/&
     Asc2D("o3_3km"    ,"BCVugXX",IXADV_O3   ,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",to_ug_ADV(IXADV_O3) ,600.0*2.0),&
     Asc2D("no_3km"    ,"BCVugXX",IXADV_NO   ,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",to_ug_ADV(IXADV_NO) ,-999.9),&
     Asc2D("no2_3km"   ,"BCVugXX",IXADV_NO2  ,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",to_ug_ADV(IXADV_NO2),600.0*1.91),&
     Asc2D("so2_3km"   ,"BCVugXX",IXADV_SO2  ,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",to_ug_ADV(IXADV_SO2),-999.9),&
     Asc2D("co_3km"    ,"BCVugXX",IXADV_CO   ,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",to_ug_ADV(IXADV_CO) ,-999.9),&
     Asc2D("Rn222_3km" ,"BCVugXX",rn222,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",to_ug_ADV(rn222)    ,-999.9),&
     Asc2D("pm25_3km"  ,"BCVugXXgroup",pm25,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",1.0                 ,-999.9),&
     Asc2D("pm10_3km"  ,"BCVugXXgroup",pm10,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",1.0                 ,-999.9),&
     Asc2D("pm_h2o_3km","PMwater",00         ,&
           ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",1.0                 ,-999.9),&
!! Partial/Full COLUMN/COLUMgroup calculations:kk$
!!   hr_out%nk indecate the number of levels in the column,
!!     1<%nk<KMAX_MID  ==>  Partial column: %nk lowermost levels
!!       otherwise     ==>  Full column: all model levels
     Asc2D("no2_col"   ,"COLUMN",IXADV_NO2  ,&
           ix1,ix2,iy1,iy2,1,"ug",to_ug_ADV(IXADV_NO2)             ,-999.9),&
!!         ix1,ix2,iy1,iy2,1,"1e15molec/cm2",to_molec_cm2*1e-15    ,-999.9),
     Asc2D("AOD_550nm" ,"AOD"   ,00         ,&
           ix1,ix2,iy1,iy2,1," ",1.0                               ,-999.9)/)
  case("EVA2010")
!   ix1=IRUNBEG;ix2=IRUNBEG+GIMAX-1
!   iy1=JRUNBEG;iy2=JRUNBEG+GJMAX-1
    levels_hourly = (/0/)
    pm25 =find_index("SURF_ug_PM25X_rh50",f_2d(:)%name)
    pm10 =find_index("SURF_ug_PM10_rh50" ,f_2d(:)%name)
!**         name     type     ofmt    ispec    
!**         ix1 ix2 iy1 iy2 nk sellev? unit conv  max
    hr_out = (/&
      Asc2D("o3"  ,"BCVugXX",IXADV_O3   ,&
            ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",to_ug_ADV(IXADV_O3) ,600.0*2.0),&
      Asc2D("no2" ,"BCVugXX",IXADV_NO2  ,&
            ix1,ix2,iy1,iy2,NLEVELS_HOURLY,"ug",to_ug_ADV(IXADV_NO2),600.0*1.91),&
      Asc2D("pm25","D2D",pm25,&
            ix1,ix2,iy1,iy2,1             ,"ug",1.0                 ,-999.9),&
      Asc2D("pm10","D2D"    ,pm10,&
            ix1,ix2,iy1,iy2,1             ,"ug",1.0                 ,-999.9)/)
!  if(MasterProc)then
!  write(*,*)(i,hr_out(i),i=1,nhourly_out)
!  call CheckStop("AMVB")
!  endif
  case("3DPROFILES")
   ! nb Out3D uses totals, e.g. O3, not IXADV_O3
   ! Number of definitions must match nhourly_out set above
    levels_hourly = (/ (i, i= 0,nlevels_hourly-1) /)  ! -1 will give surfac
    hr_out= (/ &
        Asc2D("o3_3dppb"    ,"Out3D",O3   ,&
             ix1,ix2,iy1,iy2,nlevels_hourly,"ppbv", PPBINV,600.0*2.0) &
       ,Asc2D("no2_3dppb"   ,"Out3D",&
             NO2  ,ix1,ix2,iy1,iy2,nlevels_hourly,"ppbv",PPBINV ,600.0*1.91) &
!       ,Asc2D("o3_3dug"   ,"Out3D",&
!         O3, ix1,ix2,iy1,iy2,nlevels_hourly,"ug",to_ug_ADV(IXADV_O3) ,600.0*2.0) &
    /)

    if(MasterProc ) write(*,*) "TESTHH 3D O3 SET", nlevels_hourly
  case("TFMM")

!**               name     type     ofmt
!**               ispec    ix1 ix2 iy1 iy2 nk sellev? unit conv  max

!!!! It seems easiest to just use many variables as given in the d_2d arrays. Thus
!!!! we search for the name as given there using "find_index" below.
!!!! (As a test I tried both pmfine two ways, one as D2D and the other
!!!!  as ADVugXXXgroup. The results were identical.)

   hr_out = (/  &
      Asc2D("o3_3m", "ADVppbv", IXADV_o3, &
            ix1,ix2,iy1,iy2,1, "ppbv",PPBINV,600.0) &
     ,Asc2D("no3_f"  ,"ADVugXX",IXADV_NO3_F  ,&
            ix1,ix2,iy1,iy2,1,"ug/m3",to_ug_ADV(IXADV_NO3_F)  ,-999.9) &
     ,Asc2D("no3_c"  ,"ADVugXX",IXADV_NO3_C  ,&
            ix1,ix2,iy1,iy2,1,"ug/m3",to_ug_ADV(IXADV_NO3_C)  ,-999.9) &
     ,Asc2D("nh4_f"  ,"ADVugXX",IXADV_NH4_F  ,&
            ix1,ix2,iy1,iy2,1,"ug/m3",to_ug_ADV(IXADV_NH4_F)  ,-999.9) &
     ,Asc2D("so4_f"  ,"ADVugXX",IXADV_SO4  ,&
            ix1,ix2,iy1,iy2,1,"ug/m3",to_ug_ADV(IXADV_SO4)  ,-999.9) &
   ! Organics
     ,Asc2D("OM25_3m" ,"D2D", find_index("SURF_ug_PART_OM_F",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("EC25_3m" ,"D2D", find_index("SURF_ug_ECFINE",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("OM_c_3m" ,"D2D", find_index("SURF_ug_OMCOARSE",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("EC_c_3m" ,"D2D", find_index("SURF_ug_ECCOARSE",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("NatDust_f","D2D",find_index("SURF_ug_DUST_NAT_F",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("NatDust_c","D2D",find_index("SURF_ug_DUST_NAT_C",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("RoadDust_f","D2D",find_index("SURF_ug_DUST_ROAD_F",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("RoadDust_c","D2D",find_index("SURF_ug_DUST_ROAD_C",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("SeaSalt_f","D2D",find_index("SURF_ug_SEASALT_F",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("SeaSalt_c","D2D",find_index("SURF_ug_SEASALT_C",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
    ! Sums
     ,Asc2D("PM25_3m" ,"D2D", find_index("SURF_ug_PM25_rh50",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("PMFINE"  ,"D2D", find_index("SURF_ug_PMFINE",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("pm_h2o_3m","D2D",find_index("SURF_PM25water",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1,"ug",1.0,-999.9) &
     ,Asc2D("PM25dry_3m","D2D",find_index("SURF_ug_PM25",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("PM10_3m" ,"D2D", find_index("SURF_ug_PM10_rh50",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("PM10dry_3m","D2D",find_index("SURF_ug_PM10",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "ug/m3",1.0,-999.9) &
     ,Asc2D("no_3m" ,"ADVppbv", IXADV_NO,&
            ix1,ix2,iy1,iy2,1,"ppbv",PPBINV ,600.0*1.91) &
     ,Asc2D("no2_3m","ADVppbv", IXADV_NO2,&
            ix1,ix2,iy1,iy2,1,"ppbv",PPBINV ,600.0*1.91) &
     ,Asc2D("T2_C",   "T2_C"  , 00, &
            ix1,ix2,iy1,iy2,1, "degC",1.0   ,100.0) &
     ,Asc2D("ws_10m", "ws_10m", 00, &
            ix1,ix2,iy1,iy2,1, "m/s",1.0   ,100.0) &
     ,Asc2D("HMIX","D2D", find_index("HMIX",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "m",1.0,10000.0) &
     ,Asc2D("USTAR_NWP","D2D", find_index("USTAR_NWP",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "m/s",1.0,-999.9) &
     ,Asc2D("Kz_m2s","D2D",find_index("Kz_m2s",f_2d(:)%name), &
            ix1,ix2,iy1,iy2,1, "m2/s",1.0,-999.9) &
     /)
  case default
    hr_out = (/  &
      Asc2D("o3_3m", "ADVppbv", IXADV_o3, &
                   ix1,ix2,iy1,iy2,1, "ppbv",PPBINV,600.0) &
    /)
  endselect

  !/** Consistency check:
  ! Was the array set?       R: %name/=none
  ! Was the D2D/Group found? R: %spec>0
  do i=1,nhourly_out
    if(MasterProc) write(*,*) "TESTHH O3 ATEND", i, nlevels_hourly
    if(hr_out(i)%name/="none".and.hr_out(i)%spec>=0) cycle
    write(errmsg,"(A,2(1X,A,'=',I0),2(1X,A,':',A))")&
      "set_output_defs: Failed consistency check",&
      "Hourly",i, "Nhourly",NHOURLY_OUT,&
      "Name",trim(hr_out(i)%name),"Type",trim(hr_out(i)%type)
    call CheckStop(errmsg)
  enddo

  !/** Wanted dates for instantaneous values output:
  !    specify months,days,hours for which full output is wanted.
  wanted_dates_inst(1) = date(-1,1,1,0,0)
  wanted_dates_inst(2) = date(-1,1,1,3,0)
  wanted_dates_inst(3) = date(-1,1,1,6,0)

end subroutine set_output_defs

end module My_Outputs_ml
