! <My_Derived_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

!==============================================================================
module My_Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module specifies the "derived" fields, such as accumulated
  ! precipitation
  ! or sulphate, daily, monthly or yearly averages, depositions. These fields
  ! are all typically output as netCDF fields.
  !
  ! This module provides the user-defined setups which are used in Derived_ml.
  ! Derived fields are identified by a "class", such as "ADV" of "VOC", and
  ! the Derived_ml should perform any integrations for this.
  !
  ! Several often-used routines (e.g. for AOTs, acc. sulphate, are defined
  ! in the Derived_ml.f90, but users can define their own here, since
  ! we do not use "use only" in Derived_ml.
  !
  !   Only text strings used here to define wanted data
  !   All data field characteristics should be defined in Derived_ml, e.g.
  !   in f_2d arrays.
  !   Derived fields such as d_2d only exist in Derived_ml, so are
  !   accessed here through subroutine calls - using just the (i,j) part
  !   of the bigger d_2d arrays
  !---------------------------------------------------------------------------

use AOTx_ml, only : O3cl, VEGO3_OUTPUTS, VEGO3_DEFS
use CheckStop_ml,  only: CheckStop, StopAll
use Chemfields_ml, only : xn_adv, xn_shl, cfac
use ChemSpecs_adv_ml        ! Use IXADV_ indices...
use ChemSpecs_shl_ml        ! Use IXSHL_ indices...
use ChemSpecs_tot_ml        ! eg SO2, SO4, HCHO,  ....
use ChemGroups_ml  ! Allow all groups to ease compilation
                   !,  eg. OXN_GROUP, DDEP_OXNGROUP, BVOC_GROUP
use ChemChemicals_ml, only : species               !  For mol. wts.
use ChemSpecs_adv_ml         ! Use NSPEC_ADV, IXADV_ indices
use EmisDef_ml,     only :  EMIS_FILE
use GridValues_ml, only : debug_li, debug_lj, debug_proc
use Io_Progs_ml,   only: PrintLog
use LandDefs_ml,  only : LandDefs, LandType, Check_LandCoverPresent ! e.g. "CF"
use MetFields_ml,        only : z_bnd, roa
use ModelConstants_ml, only : ATWAIR  &
                        , SOX_INDEX, OXN_INDEX, RDN_INDEX &
                        , MasterProc  &
                        , SOURCE_RECEPTOR  &
                        , USE_SOILNOX &
                        , DEBUG => DEBUG_MY_DERIVED &
                        , M=>IOU_MON, D=>IOU_DAY, H=>IOU_HOUR &
                        , KMAX_MID & ! =>  z dimension
                        , PPBINV  &  !   1.0e9
                        , MFAC       ! converts roa (kg/m3 to M, molec/cm3)
use MosaicOutputs_ml, only : nMosaic, MAX_MOSAIC_OUTPUTS, MosaicOutput, & !
  Init_MosaicMMC,  Add_MosaicMetConcs, &
  Add_NewMosaics, &
  Add_MosaicVEGO3, &
  Add_MosaicDDEP, &
  MMC_USTAR, MMC_INVL, MMC_RH, MMC_CANO3, MMC_VPD, MMC_FST, MMC_GSTO, MMC_EVAP

use OwnDataTypes_ml, only : Deriv, print_deriv_type, TXTLEN_DERIV, &
           TXTLEN_SHORT, typ_ss, typ_s3, typ_s4, typ_s5i
use Par_ml,    only: me, MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     limax, ljmax           ! => used x, y area
use SmallUtils_ml,  only : AddArray, LenArray, NOT_SET_STRING, WriteArray, &
                            find_index
use TimeDate_ml,   only : current_date
implicit none
private

 public  :: Init_My_Deriv
 public  :: My_DerivFunc ! Miscelleaneous functions of xn_adv for output
                         ! (not currently used)


   !/** Depositions are stored in separate arrays for now - to keep size of
   !    derived arrays smaller and to allow possible move to a Deposition
   !    module at a later stage.
   !  Factor 1.0e6 converts from kg/m2/a to mg/m2/a

  !        We normally distinguish source-receptor (SR) stuff from model
  !        evaluation.  The SR runs should use as few as possible outputs
  !        to keep CPU and disc-requirements down. We define first then the
  !        minimum list of outputs for use in SR, then define an extra list
  !        of parameters needed in model evaluation, or even for the base-case
  !        of SR runs.

    logical, parameter, private :: T=.true., F=.false.

  !============ parameters for concentration outputs ========================!

    integer, public, parameter :: MAX_NUM_DERIV2D = 200
    integer, public, parameter :: MAX_NUM_DERIV3D =   5
    character(len=TXTLEN_DERIV), public, save, &
         dimension(MAX_NUM_DERIV2D) :: wanted_deriv2d = NOT_SET_STRING
    character(len=TXTLEN_DERIV), public, save, &
         dimension(MAX_NUM_DERIV3D) ::  wanted_deriv3d = NOT_SET_STRING

    integer, private, save :: mynum_deriv2d
    integer, private, save :: mynum_deriv3d


!Mass-outputs of advected species, will be added to Derived
! time-res: M -> monthly, D-> daily....

  ! some shorthands for this table
   character(len=TXTLEN_SHORT), private, parameter ::&
        D2    = "2d", D3 = "3d", SPEC  = "SPEC", GROUP ="GROUP"

   !REMEMBER - KEEP UPPER CASE FOR ALL GASES
   type(typ_s5i), public, save, dimension(MAX_NUM_DERIV2D) :: OutputFields
   integer, public, save :: nOutputFields = 0

   type(typ_s5i), public, parameter, dimension(71) :: &
      OutputConcs = (/  &
!
! Here we use the 4th text field to give the "class" or "typ". Derived_ml
! will look for this in the select case (typ)
!CAN allow lower case ...
         typ_s5i("HMIX      ", "m",   D2,"HMIX     ","MISC", H)& !hourly tests
        ,typ_s5i("USTAR_NWP ", "m/s", D2,"USTAR_NWP","MISC", H)& !hourly tests
        ,typ_s5i("Kz_m2s    ", "m2/s",D2,"Kz_m2s   ","MISC", H)& !hourly tests
        ,typ_s5i("ws_10m    ", "m",   D2,"ws_10m   ","MISC", H)& !hourly tests
        ,typ_s5i("rh2m      ", "m",   D2,"rh2m     ","MISC", H)& !hourly tests
        ,typ_s5i("T2m       ", "degC",D2,"T2m      ","MISC", D)&
        ,typ_s5i("Snow_m    ", "m",   D2,"SNOW     ","MISC", D)&
        ,typ_s5i("SURF_ppbC_VOC  ", "ppb", D2,"VOC      ","MISC", D)&!?CHECK??
        ,typ_s5i("SMI_deep  ", "-",   D2,"SMI_deep ","MISC", D)&
        ,typ_s5i("SMI_uppr  ", "-",   D2,"SMI_uppr ","MISC", D)& !tno10
        ,typ_s5i("CO        ", "ppb", D2,"AIR_CONCS", SPEC, D)&
       ! ug/m3
        ,typ_s5i("SO2       ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("NH3       ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("HNO3      ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("NO2       ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("NO        ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("SO4       ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("NO3_F     ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("NO3_C     ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("NH4_F     ", "ug ", D2,"AIR_CONCS", SPEC, D)&  !tno20
        ,typ_s5i("SEASALT_F ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("SEASALT_C ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("POLLEN_B  ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("DUST_ROAD_F ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("DUST_ROAD_C ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("DUST_WB_F ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("DUST_WB_C ", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("DUST_SAH_F", "ug ", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("DUST_SAH_C", "ug ", D2,"AIR_CONCS", SPEC, D)&
       ! ppb
        ,typ_s5i("O3        ", "ppb", D2,"AIR_CONCS", SPEC, D)& ! test 3d !tno30
        ,typ_s5i("NO        ", "ppb", D2,"AIR_CONCS", SPEC, D)& !20 also have ugN
        ,typ_s5i("NO2       ", "ppb", D2,"AIR_CONCS", SPEC, D)& ! also have ugN
        ,typ_s5i("NH3       ", "ppb", D2,"AIR_CONCS", SPEC, D)& ! also have ugN
        ,typ_s5i("HNO3      ", "ppb", D2,"AIR_CONCS", SPEC, D)& ! also have ugN
        ,typ_s5i("SO2       ", "ppb", D2,"AIR_CONCS", SPEC, D)& ! also have ugN
        ,typ_s5i("HCHO      ", "ppb", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("C5H8      ", "ppb", D2,"AIR_CONCS", SPEC, D)&
       ! ugC/m3
! GenChem produces a number of groups of species.
! Here we say which ones we want for different units
! ****** UPPER CASE ONLY ************
! Sorry, this is a limitation that GenChem converts all names to
! uppercase:
        ,typ_s5i("OXN       ",  "ugN", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("NOX       ",  "ugN", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("RDN       ",  "ugN", D2,"AIR_CONCS", GROUP, D)&  !tno40
        ,typ_s5i("TNO3      ",  "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("SIA       ",  "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("PMFINE    ",  "ug ", D2,"AIR_CONCS", GROUP, D)& !30
        ,typ_s5i("PM10      ",  "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("PMCO      ",  "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("PPM25     ",  "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("PPM_C     ",  "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("SS        ",  "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("DUST_NAT_F",  "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("DUST_NAT_C",  "ug ", D2,"AIR_CONCS", GROUP, D)& !tno50
        ,typ_s5i("DUST      ",  "ug ", D2,"AIR_CONCS", GROUP, D)&
       ! SOA, PCM_F etc. are special and need appropriate units. Do
       ! not confuse! Only PCM has proper ug units, the others are
       ! carbon-eqiuvalents (PCM is particulate carbonaceous matter
       ! = sum of all EC and OM components.)
       ! ,typ_s5i("AER_ASOA  ", "ugC", D2,"AIR_CONCS", SPEC, D)&  !! ALWAYS as ugC
       ! ,typ_s5i("AER_BSOA  ", "ugC", D2,"AIR_CONCS", SPEC, D)&
         !typ_s5i("DUST      ",  "ug ", D2,"AIR_CONCS", GROUP, D),&   !#35
        ,typ_s5i("PPM25_FIRE", "ug", D2,"AIR_CONCS", GROUP,  D) &
       ! ============================================================
       ! SOA additions (26 entries)
        ,typ_s5i("PART_OM_F  ", "ug ", D2,"AIR_CONCS", SPEC, D)&  !! NEVER as ugC !!
        ,typ_s5i("OMCOARSE   ", "ug ", D2,"AIR_CONCS", GROUP, D)&  !! NEVER as ugC !!
        ,typ_s5i("ECFINE    ", "ug ", D2,"AIR_CONCS", GROUP, D)&
        ,typ_s5i("ECCOARSE  ", "ug ", D2,"AIR_CONCS", GROUP, D)&

        ,typ_s5i("PART_OC10  ", "ug ", D2,"AIR_CONCS", SPEC, M)&  !! NEVER as ugC !!
        ,typ_s5i("PART_OC25  ", "ug ", D2,"AIR_CONCS", SPEC, M)&  !! NEVER as ugC AND note that for nonvolatile type VBS runs (NPNA etc) this lacks the FFUELOC component!!
!X        ,typ_s5i("EC_F      ", "ug ", D2,"AIR_CONCS", GROUP, D)&
       ! SOA, PCM_F etc. are special and need appropriate units. Do
       ! not confuse! Only PCM has proper ug units, the others are
       ! carbon-eqiuvalents (PCM is particulate carbonaceous matter
       ! = sum of all EC and OM components.)
        ,typ_s5i("PART_ASOA_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&  !! ALWAYS as ugC
        ,typ_s5i("PART_BSOA_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&  !tno60
!         ,typ_s5i("PART_FFUELOA25_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&
!         ,typ_s5i("PART_WOODOA25_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&
!        ,typ_s5i("PART_FFIREOA25_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&
!        ,typ_s5i("EC_F_FFUEL_NEW", "ug", D2,"AIR_CONCS", SPEC, D)&
!        ,typ_s5i("EC_F_FFUEL_AGE", "ug", D2,"AIR_CONCS", SPEC, D)&
!        ,typ_s5i("EC_C_FFUEL", "ug", D2,"AIR_CONCS", SPEC, M)&
!        ,typ_s5i("NONVOL_BGNDOC", "ug", D2,"AIR_CONCS", SPEC, D)&
!        ,typ_s5i("NONV_FFUELOC_COARSE", "ug", D2,"AIR_CONCS", SPEC, D)&
        ,typ_s5i("PART_ASOA_OM", "ug", D2,"AIR_CONCS", SPEC, D)& !NEVER as ugC!
        ,typ_s5i("PART_BSOA_OM", "ug", D2,"AIR_CONCS", SPEC, D)& !NEVER as ugC! !tno70
         !zero for NONVOL:
!        ,typ_s5i("PART_FFUELOA25_OM", "ug", D2,"AIR_CONCS", SPEC, D)& !NEVER as ugC!
!        ,typ_s5i("PART_WOODOA25_OM", "ug", D2,"AIR_CONCS", SPEC, D)& !NEVER as ugC!
!        ,typ_s5i("PART_FFIREOA25_OM", "ug", D2,"AIR_CONCS", SPEC, M)& !NEVER as ugC!
         !Sep16 tests
        ,typ_s5i("FFIRE_BC"     , "ug", D2,"AIR_CONCS", SPEC, D)& !NEVER as ugC!
        ,typ_s5i("FFIRE_REMPPM25", "ug", D2,"AIR_CONCS", SPEC, D)& !NEVER as ugC!
        ,typ_s5i("FFIRE_OM"      , "ug", D2,"AIR_CONCS", SPEC, D)& !NEVER as ugC!
!SPECIAL PM25 will be sum of fine + fraction coarse
! PUT AT END OF THIS LIST  !
!---------------------------------------------
        ,typ_s5i("SURF_ug_PM25",  "ug" ,  D2,"PM25     ","MISC", D)&
        ,typ_s5i("SURF_ug_PM25X",  "ug" ,  D2,"PM25X     ","MISC", D)&
        ,typ_s5i("SURF_ug_PM25X_rh50",  "ug" ,  D2,"PM25X_rh50","MISC", D)&
        ,typ_s5i("SURF_ug_PM25_rh50" ,  "ug" ,  D2,"PM25_rh50 ","MISC", D)&
        ,typ_s5i("SURF_ug_PM10_rh50" ,  "ug" ,  D2,"PM10_rh50 ","MISC", D)&  !tno78

!---------------------------------------------
        ,typ_s5i("RN222     ", "ppb", D2,"AIR_CONCS", SPEC, D)&
       ! ============================================================
       /)
!TFMM         typ_s5i("SO2       ", "ugS", D2,"AIR_CONCS", SPEC, D)&
!TFMM        ,typ_s5i("SO4       ", "ugS", D2,"AIR_CONCS", SPEC, D)&
! Could also add from earlier D2_EXTRA array:
      !,"u_ref             " &
!.... down to here
!        ,typ_s5i("RNWATER   ", "ppb", D2,"AIR_CONCS", SPEC, D)&
! Omit for CityZen
!TFMM just keep ppb version below
!TFMM        ,typ_s5i("NO        ", "ugN", D2,"AIR_CONCS", SPEC, D)&
!TFMM        ,typ_s5i("NO2       ", "ugN", D2,"AIR_CONCS", SPEC, D)&
!TFMM        ,typ_s5i("NH3       ", "ugN", D2,"AIR_CONCS", SPEC, D)&
!TFMM        ,typ_s5i("HNO3      ", "ugN", D2,"AIR_CONCS", SPEC, D)&
!TFMM        ,typ_s5i("HONO      ", "ugN", D2,"AIR_CONCS", SPEC, D)&
!TFMM        ,typ_s5i("PAN       ", "ugN", D2,"AIR_CONCS", SPEC, D)&
       ! Remember, species have upper case, so not _f !
!TFMM just keep ug version below
!TFMM    ,typ_s5i("NO3_F     ", "ugN", D2,"AIR_CONCS", SPEC, D)&
!TFMM    ,typ_s5i("NO3_C     ", "ugN", D2,"AIR_CONCS", SPEC, D)&  ! 10 to here
!TFMM    ,typ_s5i("NH4_F     ", "ugN", D2,"AIR_CONCS", SPEC, D)&
!        ,typ_s5i("PART_SOA_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&
       ! ,typ_s5i("PART_OFFUELOA25_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&
       ! ,typ_s5i("PART_OWOODOA25_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&
       ! ,typ_s5i("PART_OFFIREOA25_OC", "ugC", D2,"AIR_CONCS", SPEC, M)&
!none yet        ,typ_s5i("EC_F_WOOD_NEW ", "ug", D2,"AIR_CONCS", SPEC, D)&
!none yet        ,typ_s5i("EC_F_WOOD_AGE ", "ug", D2,"AIR_CONCS", SPEC, D)&
       !DS  ,typ_s5i("EC_C_WOOD ", "ug", D2,"AIR_CONCS", SPEC, M)&
!Needed for forest-fire checks
!        ,typ_s5i("NONVOL_FFUELOC25", "ug", D2,"AIR_CONCS", SPEC, D)&
       !DS  ,typ_s5i("NONVOL_WOODOC25", "ug", D2,"AIR_CONCS", SPEC, D)&
!        ,typ_s5i("NONVOL_FFIREOC25", "ug", D2,"AIR_CONCS", SPEC, D)&
!        ,typ_s5i("PART_SOA_OM", "ug", D2,"AIR_CONCS", SPEC, D)& !NEVER as ugC!
!
      !CityZen Outputs
      !  ,typ_s5i("O3        ", "ug ", D2,"AIR_CONCS", SPEC, D)& ! test 3d
      !  ,typ_s5i("NO2       ", "ug ", D2,"AIR_CONCS", SPEC, D)& ! also have ugN
      !  ,typ_s5i("DUST_NAT_F", "ug ", D2,"AIR_CONCS", SPEC, D)&
      !  ,typ_s5i("DUST_NAT_C", "ug ", D2,"AIR_CONCS", SPEC, D)&
      !  ,typ_s5i("AER_OM_F  ", "ug ", D2,"AIR_CONCS", SPEC, D)&  !! NEVER as ugC !!
      !  ,typ_s5i("AER_OC    ", "ug ", D2,"AIR_CONCS", SPEC, D)&  !! NEVER as ugC !!
      !  ,typ_s5i("EC_F      ", "ug ", D2,"AIR_CONCS", GROUP, D)&
        !,typ_s5i("PART_XO_OFFLOA25_O", "ug", D2,"AIR_CONCS", SPEC, M)& !NEVER as ugC!
        !,typ_s5i("PART_XO_OWDOA25_O", "ug", D2,"AIR_CONCS", SPEC, M)& !NEVER as ugC!
        !,typ_s5i("PART_XO_OFFIOA25_O", "ug", D2,"AIR_CONCS", SPEC, M)& !NEVER as ugC!

! Tropospheric columns
   integer, public, parameter, dimension(1) :: COLUMN_MOLEC_CM2 = &
          (/ NO2 /) ! , CO, CH4, C2H6, HCHO, NO2 /)
   character(len=3), public, parameter, dimension(1) :: COLUMN_LEVELS = &
      (/  "k20" /) ! , "k16", "k12", "k08" /)

    character(len=TXTLEN_DERIV), public, parameter, dimension(4) :: &
  D2_SR = (/ &
       "SURF_MAXO3    " &
      ,"SURF_PM25water" &
      ,"SOMO35        " &
      ,"PSURF         " &  ! Surface  pressure (for cross section):
  /)


    character(len=TXTLEN_DERIV), public, parameter, dimension(1) :: &
  COL_ADD = (/ "AOD" /)


  !============ Extra parameters for model evaluation: ===================!
    !character(len=TXTLEN_DERIV), public, parameter, dimension(13) :: &
    character(len=TXTLEN_DERIV), public, parameter, dimension(5) :: &
  D2_EXTRA = (/ &
       "Area_Grid_km2     " &
      ,"Area_Conif_Frac   " &
      ,"Area_Decid_Frac   " &
      ,"Area_Seminat_Frac " &
      ,"Area_Crops_Frac   " &
!      ,"SoilWater_deep    " & ! See SMI_deep above
!      ,"SoilWater_uppr    " & ! See SMI_uppr above
!      ,"AreaPOLL          " & ! Future usage. Should change name too
  /)


 ! Ecosystem dep output uses receiver land-cover classes (LCs)
 ! which might include several landuse types, e.g. Conif in
!  D2_SO2_m2Conif.


  integer, private, save :: nOutDDep, nOutVEGO3
  integer, private, save :: nOutMET !


   ! Specify some species and land-covers we want to output
   ! depositions for in netcdf files. DDEP_ECOS must match one of
   ! the DEP_RECEIVERS  from EcoSystem_ml.
   !
   ! integer, public, parameter :: NNDRYDEP = 3 ! 7 + 1 !JUST HNO3: size(DDEP_OXNGROUP)
    !TFMM integer, public, parameter, dimension(7+size(DDEP_OXNGROUP)) :: &
    integer, public, parameter, dimension(3) :: &
    !integer, public, parameter, dimension(NNDRYDEP) :: &
      DDEP_SPECS = (/ SOX_INDEX, OXN_INDEX, RDN_INDEX /) !TFMM , & !  /) ! , &
           !SO2,  SO4, NH3, NH4_f, DDEP_OXNGROUP /)
           !SO2,  SO4, NH3, NH4_f, HNO3 /) ! DDEP_OXNGROUP /)

    character(len=TXTLEN_DERIV), public, parameter, dimension(3) :: &
      DDEP_ECOS  = (/ "Grid   " , "Conif  ", "Seminat" /) !&! "Water_D" &
                  !  , "Decid  ", "Crops  " /)

   ! Frequency of dry-dep outputs
    integer, public, parameter ::  DDEP_FREQ = D  ! (D)ay or (M)onth


  ! Have many combinations: species x ecosystems
!  type(Deriv), public, &
!     dimension( size(DDEP_SPECS)*size(DDEP_ECOS) ), save :: OutDDep

   !- specify some species and land-covers we want to output
   ! dep. velocities for in netcdf files. Set in My_DryDep_ml.

!TFMM    type(typ_s5i), public, parameter, dimension(17) :: &
    type(typ_s5i), public, parameter, dimension(1) :: &
         NewMosaic = (/ &
             typ_s5i( "Mosaic", "VG", "O3       ", "Grid","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "O3       ", "CF  ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "O3       ", "SNL ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "HNO3     ", "Grid","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "HNO3     ", "W   ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "HNO3     ", "CF  ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "HNO3     ", "SNL ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "NO3_F    ", "SNL ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "NO3_C    ", "SNL ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "NO3_F    ", "Grid","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "NO3_C    ", "Grid","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "SEASALT_F", "W   ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "SEASALT_C", "W   ","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "SEASALT_F", "Grid","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "VG", "SEASALT_C", "Grid","cms",D ) &
!TFMM            ,typ_s5i( "Mosaic", "Rs", "SO2      ", "Grid","sm",D ) &
!TFMM            ,typ_s5i( "Mosaic", "Rs", "NH3      ", "Grid","sm",D ) &
         /)

! VEGO3 outputs for PODY and AOTX - see AOTnPOD_ml for definitions,
! Any string used here must have been defined in AOTnPOD_ml.
!
    character(len=TXTLEN_DERIV), public, parameter, dimension(13) :: &
     VEGO3_WANTED  =  (/ &
         "POD1_IAM_DF    ",&
         "POD1_IAM_MF    ",&
         "POD1_DF        ",&
         "POD1_CF        ",&
         "POD3_TC        ",&
        ! "POD3_TC30d     ",&
        ! "POD3_TC55d     ",&
         "POD3_IAM_CR    ",&
        ! "POD3_IAM_CR30d ",&
        ! "POD3_IAM_CR55d ",&
        ! "POD6_IAM_CR    ",& ! Not recommended - not robust
        ! "POD6_IAM_CR30d ",&  ! Not recommended - not robust
        ! "POD6_IAM_CR55d ",& ! Not recommended - not robust
         "MMAOT40_TC     ",&
         "MMAOT40_IAM_DF ",&
         "MMAOT40_IAM_MF ",&
         "MMAOT40_IAM_CR ",&
         "EUAOT40_Crops  ", &
         "EUAOT40_Forests", &
         "MMAOT40_IAM_WH " &
    /) !NB -last not found. Could just be skipped, but kept
       !to show behaviour


! For met-data and canopy concs/fluxes ...

!TFMM    character(len=TXTLEN_DERIV), public, parameter, dimension(3) :: &
    character(len=TXTLEN_DERIV), public, parameter, dimension(1) :: &
      MOSAIC_METCONCS = (/ "USTAR" /) ! TFMM "VPD     "  &
                         ! ,"CanopyO3" & !SKIP
         !,"VPD     ", "FstO3   " "EVAP    ", "Gsto    " &
                        !SKIP
   !TFMM                    ,"USTAR   ", "INVL    "  &
   !TFMM                    /)
                          ! "g_sto" needs more work - only set as L%g_sto

    character(len=TXTLEN_DERIV), public, save, dimension(2) :: &
      MET_LCS  = (/ "DF    ", "GR    " /) !, "CF    ", "BF    ", "NF    " /) !,
                                !"IAM_DF", "IAM_MF"/)

      !MET_LCS  = (/ "GR    " , "IAM_CR", "IAM_DF", "IAM_MF"/)
    !character(len=TXTLEN_DERIV), public, parameter, dimension(5) :: &
      !MET_LCS  = (/ "CF", "SNL", "TESTCF", "GR" ,"TC"/)
      ! Can also set dim 4:1 to exclude all - gives zero size MET_LCS

   ! We use Dep_type anyway since we can use the LCC elements
!    type(Dep_type), public, & !Non-stomatal conductance
!     dimension( size(MOSAIC_METCONCS)*size( MET_LCS) ),  save :: OutMET

!----------------------


   type(typ_s3), dimension(7-3+2), public, parameter :: WDEP_WANTED = (/ &
         typ_s3( "PREC     ", "PREC ", "mm  " )  &
        ,typ_s3( "SOX      ", "GROUP", "mgS " )  & ! Will get WDEP_SOX group
        ,typ_s3( "OXN      ", "GROUP", "mgN " )  &
        ,typ_s3( "RDN      ", "GROUP", "mgN " )  &
!TFMM        ,typ_s3( "SS       ", "GROUP", "mg  " )  &
      !
        ,typ_s3( "SO2      ", "SPEC ", "mgS ") &  ! Makes WPEP_SO2
      !  ,typ_s3( "SO4      ", "SPEC ", "mgS ") &
        ,typ_s3( "HNO3     ", "SPEC ", "mgN ") &
      !  ,typ_s3( "NO3_F    ", "SPEC ", "mgN ") &
      !  ,typ_s3( "NO3_C    ", "SPEC ", "mgN ") &
!TFMM        ,typ_s3( "NH4_F    ", "SPEC ", "mgN ") &
!TFMM        ,typ_s3( "NH3      ", "SPEC ", "mgN ") &
      !  ,typ_s3( "SEASALT_F", "SPEC ", "mg  ") &
      ! ,typ_s3( "SEASALT_C", "SPEC ", "mg  ") &
     /)


    ! For some reason having this as a parameter caused problems for
    ! PC-gfortran runs.

    ! other (non-ppb) 3D output, set as zero-size (eg 4:1) for normal runs
     character(len=TXTLEN_DERIV), public, save, dimension(4:1) :: &
     !character(len=TXTLEN_DERIV), public, save, dimension(1) :: &
       D3_OTHER  != (/ "D3_PM25water" /) !**** Under construction *******
     != (/ "D3_m_TH", "D3_m2s_Kz" /)

    integer, private :: i,j,k,n, ivoc, index    ! Local loop variables

   contains

 !=========================================================================
  subroutine Init_My_Deriv()

    integer :: i, itot, nDD, nMET, nVEGO3, n1, n2,istat, nMc
    character(len=TXTLEN_DERIV) :: txt
    character(len=TXTLEN_DERIV), &
    dimension(size(COLUMN_MOLEC_CM2)*size(COLUMN_LEVELS)) ::&
            tmpname ! e.g. DDEP_SO2_m2Conif
    character(len=100) :: errmsg
    character(len=TXTLEN_DERIV), &
       dimension(size( OutputConcs(:)%txt1 ) ) ::&
          tag_name    ! Needed to concatanate some text in AddArray calls
                      ! - older (gcc 4.1?) gfortran's had bug
    character(len=TXTLEN_SHORT) :: outname, outunit, outdim, outtyp, outclass

if(MasterProc ) print *, "TESTHH INSIDE Init_My_Deriv"

    call Init_MosaicMMC(MOSAIC_METCONCS)  ! sets MMC_USTAR etc.


   ! Build up the array wanted_deriv2d with the required field names

     call AddArray( "WDEP_" // WDEP_WANTED(:)%txt1, wanted_deriv2d, NOT_SET_STRING,errmsg)
     call CheckStop( errmsg, errmsg // "WDEP_WANTED too long" )
!TEST     call AddArray( D2_SR,  wanted_deriv2d, NOT_SET_STRING, errmsg)
!TEST     call CheckStop( errmsg, errmsg // "D2_SR too long" )
     call AddArray( COL_ADD,  wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "COL_ADD too long" )

  ! Emission sums - we always add these (good policy!)
   do  i = 1, size(EMIS_FILE)
     tag_name(1) = "Emis_mgm2_" // trim(EMIS_FILE(i))
     call AddArray( tag_name(1:1), wanted_deriv2d, NOT_SET_STRING, errmsg)
   end do
   do  i = 1, size(BVOC_GROUP)
     itot = BVOC_GROUP(i)
     tag_name(1) = "Emis_mgm2_BioNat" // trim(species(itot)%name)
     call AddArray( tag_name(1:1), wanted_deriv2d, NOT_SET_STRING, errmsg)
   end do
  if ( USE_SOILNOX ) then
     tag_name(1) = "Emis_mgm2_BioNatNO"
     call AddArray( tag_name(1:1), wanted_deriv2d, NOT_SET_STRING, errmsg)
   end if



! Do SR last, so we get PM25 after groups have been done
     call AddArray( D2_SR,  wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "D2_SR too long" )

     if ( .not. SOURCE_RECEPTOR ) then !may want extra?
        call AddArray( D2_EXTRA, wanted_deriv2d, NOT_SET_STRING, errmsg)
        call CheckStop( errmsg, errmsg // "D2_EXTRA too long" )
     end if


     ! Column data:
     n = 0
     do n1 = 1, size(COLUMN_MOLEC_CM2)
     do n2 = 1, size(COLUMN_LEVELS)
       n = n + 1
       tmpname(n) = "COLUMN_" // trim( species(COLUMN_MOLEC_CM2(n1))%name ) &
          // "_" // COLUMN_LEVELS(n2)
     end do
     end do
     call AddArray(tmpname, wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "COLUMN too long" )

     ! Didn't work:
     !call AddArray( "COLUMN_" // trim( species(COLUMN_MOLEC_CM2(:))%name ), &
     !  wanted_deriv2d, NOT_SET_STRING)


      !------------- Depositions to ecosystems --------------------------------

      call Add_MosaicDDEP(DDEP_ECOS,DDEP_SPECS,DDEP_FREQ,nDD)
      nOutDDep = nDD

      !------------- VEGO3 stuff ----------------------------------------------
      ! For fluxes or AOTs we start with a formatted name, eg. POD_3.0_CF and
      !untangle it to get threshold Y (=3.0) and landcover type

      allocate(VEGO3_OUTPUTS( size(VEGO3_WANTED) ), stat=istat)
      if( DEBUG .and. istat /= 0 ) then
         print *, "My_Derived ISTAT ERR VEGO3"
      end if

      do n = 1, size(VEGO3_WANTED)
         n1 = find_index(VEGO3_WANTED(n),VEGO3_DEFS(:)%name)
         call CheckStop(  n1>size(VEGO3_DEFS(:)%name) .or. n1<1 , &
                   "VEGO3 not found"//trim(VEGO3_WANTED(n)) )
         VEGO3_OUTPUTS(n) = VEGO3_DEFS(n1)
       if(DEBUG .and. MasterProc)  write(*,*) "VEGO3 NUMS ", n, n1,&
            trim( VEGO3_WANTED(n) )
      end do
      if(MasterProc)call WriteArray(VEGO3_OUTPUTS(:)%name,size(VEGO3_WANTED)," VEGO3 OUTPUTS:")
      call Add_MosaicVEGO3(M, nVEGO3)  ! M=monthly
      nOutVEGO3 = nVEGO3

      !----- some "luxury outputs" -------------------------------------------

 if( .not.SOURCE_RECEPTOR)then

      !------------- Deposition velocities -----------------------------------

      call Add_NewMosaics(NewMosaic, nMc)
       if(DEBUG .and. MasterProc)  write(*, *) "NEWMOSAIC   NUM ", nMc
      !nOutMos = nMos

       if(DEBUG .and. MasterProc)  print *, "VEGO3 FINAL NUM ", nVEGO3

      !------------- Met data for d_2d -------------------------
      ! We find the various combinations of met and ecosystem,
      ! adding them to the derived-type array LCC_Met (e.g. => Met_CF)
      !FEB2011  Daiyl output asked for just now. Change larer

      call Add_MosaicMetConcs(MOSAIC_METCONCS,MET_LCS, D, nMET)
      nOutMET = nMET !not needed?
  end if ! SOURCE_RECEPTOR


      !------------- end LCC data for d_2d -------------------------

     call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics" )
     call AddArray( MosaicOutput(1:nMosaic)%name, &
                        wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "MosaicOutput too long" )

     mynum_deriv2d  = LenArray( wanted_deriv2d, NOT_SET_STRING )


   ! Add the pollutants wanted from OutputConcs:
   ! to both OutputFields and wanted_deriv arrays (TOO MESSY)
   ! Requested species which are not present will trigger warnings
   !type(typ_s5i), public, parameter, dimension(27) :: &
   !   OutputConcs = (/  typ_s5i("SO2", "ugS", D2,"AIR_CONCS", SPEC, M),&
   !                     typ_s5i("SO4", "ugS", D2,"AIR_CONCS", SPEC, M),&

      do n = 1, size( OutputConcs(:)%txt1 )

         outname= trim(OutputConcs(n)%txt1)
         outunit= trim(OutputConcs(n)%txt2)
         outdim = trim(OutputConcs(n)%txt3)
         outtyp = trim(OutputConcs(n)%txt4)
         outclass = trim(OutputConcs(n)%txt5) ! MISC or SPEC or GROUP

         if( outdim == "3d" ) txt = "D3"  ! Will simplify later
         if( outdim == "2d" ) txt = "SURF"  ! Will simplify later

         if( outclass == "MISC" ) then

              tag_name(1)= trim(outname) ! Just use raw name here

              call AddArray( tag_name(1:1) , wanted_deriv2d, &
                     NOT_SET_STRING, errmsg)
              nOutputFields = nOutputFields + 1
              OutputFields(nOutputFields) = OutputConcs(n)

         else if( outtyp == "AIR_CONCS" ) then

             if( outclass == SPEC ) then ! check if available
              n1 = find_index(outname,species(:)%name)
             else if ( outclass == GROUP ) then ! check if available
              n1 = find_index(outname,GROUP_ARRAY(:)%name)
             end if

              if ( n1 < 1 ) then
                if(DEBUG.and.MasterProc) write(*,*) "Xd-2d-SKIP ", n, trim(outname)
                call PrintLog("WARNING: Requested My_Derived OutputField not found: "&
                   // " " //trim(outclass) // ":"   //trim(outname), MasterProc)
                cycle
              end if

!             end if

              tag_name(1) = "SURF_" // trim(outunit) // "_" //  trim(outname)
              call AddArray(  tag_name(1:1) , wanted_deriv2d, &
                     NOT_SET_STRING, errmsg)
              call CheckStop( errmsg, errmsg // trim(outname) // " too long" )
              nOutputFields = nOutputFields + 1
              OutputFields(nOutputFields) = OutputConcs(n)
              if(DEBUG.and.MasterProc) write(*,*) "Xd-2d-DONE ", n, trim(outname)

              if( outdim == "3d" ) then
                  tag_name(1) = "D3_" // trim(outunit) // "_" //  trim(outname)
                  call AddArray(  tag_name(1:1) , wanted_deriv3d, &
                     NOT_SET_STRING, errmsg)
                  call CheckStop( errmsg, errmsg // trim(outname) // " too long" )
                  nOutputFields = nOutputFields + 1
                  OutputFields(nOutputFields) = OutputConcs(n)
                  if(DEBUG.and.MasterProc) write(*,*) "Xd-3d-DONE ", n, trim(tag_name(1))
              end if


         else
            call StopAll("My_Deriv: Not coded yet" // &
              trim(outname) //":"//trim(outtyp) )
         end if

      end do

   ! ditto wanted_deriv3d....

     if ( .not. SOURCE_RECEPTOR ) then

       if ( size(D3_OTHER) > 0 ) then
         call AddArray( D3_OTHER,  wanted_deriv3d, NOT_SET_STRING, errmsg)
         call CheckStop( errmsg, errmsg // "Wanted D3 too long" )
       end if
     end if
! TEST HERE
     mynum_deriv2d  = LenArray( wanted_deriv2d, NOT_SET_STRING )
     mynum_deriv3d  = LenArray( wanted_deriv3d, NOT_SET_STRING )


     if ( MasterProc ) then
        if(  DEBUG ) then
           write(*,*) "Init_My_Deriv, mynum_deriv2d = ", mynum_deriv2d
           write(*,*) "Init_My_Deriv, mynum_deriv3d = ", mynum_deriv3d
           do i = 1, mynum_deriv2d
              write(*,*) "DEBUG DERIV2D ", i, mynum_deriv2d, wanted_deriv2d(i)
           end do
        end if
       call WriteArray(wanted_deriv2d,mynum_deriv2d," Required 2D output ")
       call WriteArray(wanted_deriv3d,mynum_deriv3d," Required 3D output ")

     end if

  end subroutine Init_My_Deriv
 !=========================================================================
  subroutine My_DerivFunc( e_2d, class )!  , density )

    ! We define here here any functions which cannot easily be defined
    ! in the more general Derived_ml.

  real, dimension(:,:), intent(inout) :: e_2d  !  (i,j) 2-d extract of d_2d
  character(len=*), intent(in)    :: class       ! Class of data
  integer, save :: num_warnings = 0  ! crude counter for now

  ! real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density
  ! density = 1 ( or = roa when unit ug)

  select case ( class )

    !      (not currently used)

    case ( "FRNIT" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
             ( xn_adv(IXADV_NO3_f,i,j,KMAX_MID) * cfac(IXADV_NO3_f,i,j)  &
            +  xn_adv(IXADV_NO3_c,i,j,KMAX_MID) * cfac(IXADV_NO3_c,i,j)) &
       /max(1E-80, (xn_adv(IXADV_HNO3,i,j,KMAX_MID) *  cfac(IXADV_HNO3,i,j))&
            +  xn_adv(IXADV_NO3_f,i,j,KMAX_MID) * cfac(IXADV_NO3_f,i,j)    &
            +  xn_adv(IXADV_NO3_c,i,j,KMAX_MID) * cfac(IXADV_NO3_c,i,j))
      end forall

      case  default

          if ( MasterProc .and. num_warnings < 100 ) then
            write(*,*) "My_Deriv:WARNING - REQUEST FOR UNDEFINED OUTPUT:", n, class
            num_warnings = num_warnings + 1
          end if
     end select

  end subroutine My_DerivFunc
 !=========================================================================

end module My_Derived_ml
