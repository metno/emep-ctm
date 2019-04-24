! <Config_module.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.32>
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
module Config_module
!----------------------------------------------------------------------------
! Specifies a number of constants used in the model, and reads namelist
! file to (re-)configure where possible.
! Note that physical constants (e.g. gravity, Cp, etc ( are specified in
! the module PhysicalConstants_mod.f90)
!----------------------------------------------------------------------------
use CheckStop_mod,         only: CheckStop
use ChemSpecs_mod,         only: species, CM_schemes_ChemSpecs
use Debug_module,          only: DEBUG, DebugCell
use Io_Nums_mod,           only: IO_NML, IO_LOG, IO_TMP
use OwnDataTypes_mod,      only: typ_ss, uEMEP_type, Emis_id_type, emis_in,&
                                 EmisFile_id_type, Emis_sourceFile_id_type,&
                                 TXTLEN_NAME, TXTLEN_FILE, Emis_mask_type
use Precision_mod,         only: dp
use SmallUtils_mod,        only: find_index, key2str

implicit none
private

public :: Config_ModelConstants

!=============================================================================
! Experiment name:
!  EMEPSTD      Standard run & output
!  EMEP2010    EMEPSTD with Iceland Volcanic Eruption input
!  TFMM        EMEPSTD, but with INERIS_SNAP & TFMM hourly output
!  FORECAST    Forecast run, MACC-ENS hourly output & BC
!  EVA2010     FORECAST with MACC-EVA2010 hourly output & BC
!  EMERGENCY   FORECAST with ONLY Volcanic Eruption & Nuclear Accident.
!
! We separate the concept of exp_name and the
! variable used to set the type of output in My_outputs_mod.
! The longer term solution puts the outputs into namelists
! but for now we use the MY_OUTPUTS flag. EXP_NAME can
! now be anything descriptive.

!integer, public, parameter :: &
!  TXTLEN_NAME =  64, & !for performance, should be a multiple of 8
!  TXTLEN_FILE = 200    ! large enough for paths from namelists

CHARACTER(LEN=TXTLEN_NAME), public, save :: EXP_NAME="EMEPSTD"
CHARACTER(LEN=TXTLEN_NAME), public, save :: MY_OUTPUTS="EMEPSTD"

! EMEP daily measurements end at 6am, hence we typically adjust
! for that. For global though, zero would be more normal
  integer, save, public :: END_OF_EMEPDAY = 6 ! 

  type, private :: PBL_t
    real :: ZiMIN = 100.0                     ! minimum mixing height
    real :: ZiMAX = 3000.0                    ! maximum mixing height
    character(len=10) :: HmixMethod = "JcRb"  ! Method used for Hmix 
      ! JcRb = Jericevic/Richardson number method
      ! "SbRb"= Seibert !"TIZi" = Original from Trond Iversen tiphysics
  end type PBL_t
  type(PBL_t), public, save :: PBL = PBL_t()

  type, private :: EmBio_t
    character(len=10) :: GlobBvocMethod = '-' ! can be MEGAN
    real :: IsopFac = 1.0                     ! for experiments
    real :: TerpFac = 1.0                     ! for experiments
  ! canopy light factor, 1/1.7=0.59, based on Lamb 1993 (cf MEGAN 0.57)
    real :: CLF     = 1.0                     ! canopy factor, leaf vs branch emissions
  end type EmBio_t
  type(EmBio_t), public, save :: EmBio = EmBio_t()

  ! - allows rcbio in CM_Reactions, but we access elements with
  ! the natbio indices here. These much match the indices used in rcbio
  ! We only use rcbio for isoprene and terpenes so far,  since
  ! soil NO, NH3 emissions etc are dealt with through rcemis.

  type, private :: natbio_t
    integer :: C5H8 = 1
    integer :: TERP = 2
    integer :: Nrcbio = 2  ! No. of rcbio defined in ChemFields/Biogenics_mod
    integer :: NO   = 3    ! used for EmisNat etc
    integer :: NH3  = 4
  end type natbio_t
  type(natbio_t), public, parameter :: NATBIO = natbio_t()

 ! We allow a flexible string which can switch between different
 ! experiments called by e.g. Solver. A but crude, but
 ! it makes sure the experiments are recorded in the config
 ! system

  character(len=100), save, public :: YieldModifications = 'VBS-T10' ! Default for EmChem16mt

  
  type, private :: LandCoverInputs_t
    character(len=TXTLEN_FILE), dimension(2) :: MapFile = 'NOTSET'  ! Usually PS European + global
    character(len=TXTLEN_FILE) :: LandDefs = '-'   !  LAI, h, etc (was Inputs_LandDefs
    character(len=TXTLEN_FILE) :: Do3seDefs = '-'  !  DO3SE inputs
  end type LandCoverInputs_t
  type(LandCoverInputs_t), public, save :: LandCoverInputs=LandCoverInputs_t()


! Namelist controlled:
! Some flags for model setup
!------------ NAMELIST VARIABLES - can be reset by emep_namelist.nml file
! USES system introduced June 2013--------------------------

logical, private, parameter :: F = .false.
type, public :: emep_useconfig
  character(len=10) :: testname = "STD"
  logical :: &                   ! Forest fire options
     FOREST_FIRES     = .true.  &!
    ,SOILWATER        = .false. &!
    ,BVOC             = .true.  &!triggers isoprene and terpene emissions
    ,SEASALT          = .true.  &!
    ,CONVECTION       = .false. &! false works best for Euro runs
    ,AIRCRAFT_EMIS    = .true.  &! Needs global file, see manual 
    ,LIGHTNING_EMIS   = .true.  &! 
    ,ROADDUST         = .false. &! TNO Road Dust routine. So far with simplified "climate-correction" factor
    ,DUST             = .false. &! Experimental
    ,EURO_SOILNOX     = .true.  &! ok, but diff for global + Euro runs
    ,GLOBAL_SOILNOX   = .false. &! Need to design better switch
    ,ASH          = .true.  &! Ash from historical Volcanic Eruption
    ,PreADV       = .false. &! Column Emissions are preadvected when winds are very strong 
    ,NOCHEM       = .false. &! Turns of fchemistry for emergency runs
    ,AOD          = .false. &
    ,POLLEN       = .false. &! EXPERIMENTAL. Only works if start Jan 1
    ,PROGRESS_FILES   = .false. &! write to file.msg for each output in file.nc
    ,SKIP_INCOMPLETE_OUTPUT = .false. & ! skip daily/montly/fullrun output for runs under 1/28/181 days
    ,SURF_AREA        = .true.  &! For improved aerosol uptake
    ,MACEHEADFIX      = .true.  &! Correction to O3 BCs (Mace Head Obs.)
    ,MACEHEAD_AVG     = .false. &! Uses 10-year avg. Good for e.g. RCA runs.
    ,MINCONC          = .false. &! Experimental. To avoid problems with miniscule numbers
    ,ESX              = .false. &! Uses ESX
    ,EMIS             = .false. &! Uses ESX
    ,GRIDDED_EMIS_MONTHLY_FACTOR = .false. & ! .true. triggers ECLIPSE monthly factors
    ,DEGREEDAY_FACTORS = .true.    &! will not be used if not found or global grid
    ,EMISSTACKS       = .false.     &!
    ,PFT_MAPS         = .false.  &! 
    ,EFFECTIVE_RESISTANCE = .false. ! Drydep method designed for shallow layer

 ! Mar 2017. Allow new MEGAN-like VOC
 ! Moved to emep_Config
 ! character(len=10) :: GlobBvocMethod = "-" ! MEGAN

 ! If USES%EMISTACKS, need to set:
  character(len=4) :: PlumeMethod = "none" !MKPS:"ASME","NILU","PVDI"

 ! N2O5 hydrolysis
 ! During 2015 the aersol surface area calculation was much improved, and this
 ! leads to the need for new n2o5 hydrolysis  methods. DO NOT USE EmepReimer,
 ! but one of :; 'SmixTen' , 'Smix', 'Gamma:0.002'

  character(len=20) :: n2o5HydrolysisMethod = 'Smix'

! Selection of method for Whitecap calculation for Seasalt
  character(len=15) :: WHITECAPS  = 'Callaghan'
! In development
   logical :: BIDIR       = .false.  !< FUTURE Bi-directional exchange
   character(len=20)      :: BiDirMethod = 'NOTSET'  ! FUTURE
   character(len=20)      :: MonthlyNH3  = 'NOTSET'  ! can be 'LOTOS'

end type emep_useconfig
type(emep_useconfig), public, save :: USES

logical,  public, save :: &
      FORCE_PFT_MAPS_FALSE = .false. !forces PFT_MAPS  = F, even if global grid

type(emis_in), public, dimension(50) :: emis_inputlist = emis_in()
type(Emis_sourceFile_id_type), public, save:: Emis_sourceFiles(20) !as read from config
type(Emis_mask_type), public, save :: EmisMask(10) !emission mask new format
!MaxNSECTORS to allow reading of SecEmisOutWanted before NSECTORS is defined
integer, public, parameter :: MaxNSECTORS = 100 
logical, public, save :: SecEmisOutWanted(MaxNSECTORS) = .false.
!DSHK - moved from parameter to here:
logical, public, save :: EmisSplit_OUT = .false.

logical, public, save  :: HourlyEmisOut = .false. !to output sector emissions hourly

character(len=40), public, save   :: SECTORS_NAME='SNAP'
character(len=40), public, save   :: USE_SECTORS_NAME='NOTSET'

!Note that we cannot define the settings as logical (T/F), because we need the state "NOTSET" also
character(len=TXTLEN_NAME), public, save :: EUROPEAN_settings = 'NOTSET'! The domain covers Europe
character(len=TXTLEN_NAME), public, save :: GLOBAL_settings = 'NOTSET'!The domain cover other regions

character(len=TXTLEN_FILE), public, save :: &
  EmisDir = '.',  &
  DataDir = '.',  &
  OwnInputDir = '.',  &  ! user-defined location
  GRID = 'EECCA', & ! default grid
  meteo= 'DataDir/GRID/metdata_EC/YYYY/meteoYYYYMMDD.nc', & ! template for meteofile
  DegreeDayFactorsFile = 'MetDir/HDD18-GRID-YYYY.nc'        ! template for DegreeDayFactors.nc

integer, public, save :: startdate(4)=(/0,0,0,0/),enddate(4)=(/0,0,0,24/) ! start and end of the run
integer, public, save :: out_startdate(4)=(/0,0,0,0/) ! start of the output of data
integer, public, save :: spinup_enddate(4)=(/0,0,0,0/) ! end of spinup. Does not average concentration etc before that date

!-----------------------------------------------------------
! Convection factor - reduces convective fluxes (which can be
! too high in some NWPs)
real, public, save :: CONVECTION_FACTOR = 0.33   ! Pragmatic default
!-----------------------------------------------------------
logical, public, save ::             &
  TEGEN_DATA         = .true.        & ! Interpolate global data to make dust if  USE_DUST=.true.
 ,INERIS_SNAP1       = .false.       & !(EXP_NAME=="TFMM"), & ! Switches off decadal trend
 ,INERIS_SNAP2       = .false.       & !(EXP_NAME=="TFMM"), & ! Allows near-zero summer values
 ,USE_AMINEAQ        = .false.       & ! MKPS
 ,ANALYSIS           = .false.       & ! EXPERIMENTAL: 3DVar data assimilation
 ,USE_FASTJ          = .false.       & ! use FastJ_mod for computing rcphot
!
! Output flags
 ,SELECT_LEVELS_HOURLY  = .false.    & ! for 3DPROFILES
 ,ZERO_ORDER_ADVEC  = .false.        & ! force zero order horizontal and vertical advection
 ,JUMPOVER29FEB      = .false.         ! When current date is 29th February, jump to next date.

logical, public, save :: USE_uEMEP = .false.  ! make local fraction of pollutants
type(uEMEP_type), public, save :: uEMEP ! The parameters steering uEMEP
integer, public, save :: NTIMING_uEMEP = 5 !reset to zero if USE_uEMEP = F 

integer, public, save :: &
  FREQ_HOURLY = 1  ! 3Dhourly netcdf special output frequency

! Soil NOx. Choose EURO for better spatial and temp res, but for
! global runs need global monthly. Variable USE_SOILNOX set from
! these below.
!
! Also, is scaling needed for EURO_SOILNOX?
! The Euro soil NO emissions are based upon average Nr-deposition calculated
!  for the 2000s, as given in the AnnualNdep.nc files. For future years a
!  new AnnualNdep.nc could be pre-calculated. A simpler but approximate
!  way is to scale with some other factor, e.g. the ratio of emissions over
!  some area (EMEP, or EU) in year YYYY divided by year 2005 values.
! Remember, soil-NO emissions are *very* uncertain.

  logical, public, save ::             &
!    USE_EURO_SOILNOX      = .true.     & ! ok, but diff for global + Euro runs
!   ,USE_GLOBAL_SOILNOX    = .false.    & ! Need to design better switch
   USE_SOILNOX           = .true.       ! DO NOT ALTER: Set after config
  real, public, save :: EURO_SOILNOX_DEPSCALE = 1.0 !

!NB: *OCEAN*  are internal variables. Cannot be set manually.
  logical, public, save ::  USE_OCEAN_DMS = .false. !set automatically true if found.
  logical, public, save ::  FOUND_OCEAN_DMS = .false. !set automatically true if found
  logical, public, save ::  USE_OCEAN_NH3 = .false. !set automatically true if found

! Methane background:
!  -1 gives defaults in BoundaryConditions_mod
!  -26,-45 or -85 gives RC26, P45, 85 for iyr_trend
  real, public, save :: BGND_CH4 = -1  
! To skip rct value   (jAero work)
  integer, public, save, dimension(10) :: SKIP_RCT  = -1  ! -1 gives defaults
!
  logical, public, save :: USE_WRF_MET_NAMES = .false. !to read directly WRF metdata

!Machine_config variables
 character (len=TXTLEN_FILE), public :: DataPath(20) = 'NOTSET'
!
!------------ END OF NAMELIST VARIABLES ------------------------------------!

! Some flags for model setup will be removed when code is sufficiently tested
! (for convection use foundconv in permanent code)
logical, public, parameter ::         &
  NO_CROPNH3DEP      = .true.,        & ! Stop NH3 deposition for growing crops
  USE_SOILNH3        = .false.,       & ! DUMMY VALUES, DO NOT USE!
  USE_ZREF           = .false.,       & ! testing
  EXTENDEDMASSBUDGET = .false.,       & ! extended massbudget outputs
  LANDIFY_MET        = .false.

logical, public :: &
  USE_EtaCOORDINATES=.true. ! default since October 2014

! Boundary layer profiles. IN-TESTING
character(len=4), parameter, public :: &
  FluxPROFILE = "Iter"
! FluxPROFILE = "Ln95"  ! use Launiainen1995 EXPERIMENTAL. Fails in some areas

!The GEA emission data, which is used for EUCAARI runs on the HIRHAM domains
!have in several sea grid cells non-zero emissions in other sectors than SNAP8
!and there are also NH3 emission over sea areas. The former problem makes
!the code crash if the sea areas are defined  as sea (sea=T), so we treat
!them as land in the EUCAARI/HIRHAM runs (sea=F). This is a problem with GEA
!emission data only, not the HIRHAM domain! When e.g. interpolated EMEP emissions
!are used on the HIRHAM domain, this is not a problem.

logical, public, save :: &
  SEAFIX_GEA_NEEDED = .false. ! only if problems. Read from ModelConstants_config

!=============================================================================
!+ 1) Define first dimensions that might change quite often -  for different
!     run domains

!IN-TESTING (reset in NML if wanted)
!) Emissions. Standard =ascii emislist. CdfFractions possible for INERIS
!  and new cdf emission system in testing. Reset in config_ files
Logical , save, public :: &
  EMIS_OUT    = .false.     ! output emissions in separate files (memory demanding)

integer, public :: &  ! Full domain, set automatically from meteorology
  KMAX_MID, &         ! Number of points (levels) in vertical
  KMAX_BND, &         ! Number of level boundaries (=KMAX_MID+1)
  IIFULLDOM,JJFULLDOM ! Full grid

! Sub-domains, in fulldomain coordinates. Read on ModelConstants_config
integer, public, save, dimension(4) :: &
  RUNDOMAIN =-999,  & ! run sub-domain
  fullrun_DOMAIN,   & ! fullrun (year) output sub-domain
  month_DOMAIN,     & ! montly output sub-domain
  day_DOMAIN,       & ! daily  output sub-domain
  hour_DOMAIN         ! hourly output sub-domain

! 3D Output: all modell levels will be outputed by default
! see Init_My_Deriv and OutputSize_config for details
integer, public, save ::  &
  num_lev3d=0,lev3d(60)=0 ! numbers of levels,list of levels

integer, public, save ::  & ! Actual number of processors in longitude, latitude
  NPROCX, NPROCY, NPROC     ! and total. NPROCY must be 2 for GLOBAL runs.

CHARACTER(LEN=3), public, save :: &
  DOMAIN_DECOM_MODE=''      ! override parinit(Pole_singular) option (Par_mod)

!=============================================================================
!+ 2) Define  debug flags.

! We have one variable, to say if we are on master-processor
! or not: (kept here to avoid too many dependencies for box-model
! codes which don't need Par_mod

logical, public, save ::  MasterProc = .true.
!MOVED logical, public, save ::  DebugCell  = .false.

!=============================================================================
! Some flags for model setup

! Debug flag DEBUG_XXX  applied in subroutine XXX
logical, public, parameter :: PALEO_TEST = .false. 

!=============================================================================
! 3)  Source-receptor runs?
! We don't (generally) want daily outputs for SR runs, so in
! Derived_mod, we set all IOU_DAY false if SOURCE_RECPTOR = .true..

logical, public, save :: SOURCE_RECEPTOR = .false., VOLCANO_SR=.false.

! Compress NetCDF output? (nc4 feature, use -1 for netcdf3 output)
integer, public, save :: NETCDF_DEFLATE_LEVEL=4

!Hourly output in single file or monthly/daily files:
!NB: will not work well by default on Stallo per 14th Feb 2012 because of library bugs!
!Until this is fixed, you must compile with netcdf/4.1.3 and link and run with compiler 12.1.2
character(len=30), public, save :: &    ! ending depeding on date:
! HOURLYFILE_ending="YYYYMM.nc"    ! MM  -> month (01 .. 12)
! HOURLYFILE_ending="YYYYMMDD.nc"  ! DD  -> day of the month (00 .. 31)
! HOURLYFILE_ending="JJJ.nc"   ! JJJ -> the day of the year (001 .. 366)
     HOURLYFILE_ending=".nc"           ! default, just one file
!do not use  HOURLYFILE_ending="_hourExtra.nc"     ! keep the same for the whole run

! NH3 module as set up originally with U10 from met: kept for safety only.
! Will be replaced by sub.grid calculation of wind in future.
! Keep false until code re-implemented
logical, public, parameter :: NH3_U10 = .false.

!=============================================================================
!+ 4)  Define main model dimensions,  things that will
!       generally only change when switching Met-driver
integer, public, parameter ::  &
!TREEX  NLANDUSEMAX  = 19   &   ! Number of land use types in Inputs.Landuse file
  NLANDUSEMAX  = 40    &    ! Max num land use types in Inputs.Landuse file
, KTOP         = 1     &    ! K-value at top of domain
, KWINDTOP     = 5     &    ! Define extent needed for wind-speed array
, NMET         = 2     &    ! No. met fields in memory
, KCHEMTOP     = 2     &    ! chemistry not done for k=1
, KCLOUDTOP    = 8     &    ! limit of clouds (for MADE dj ??)
, KUPPER       = 6          ! limit of clouds (for wet dep.)

integer, public :: METSTEP = 3  ! time-step of met. (h). 3 hours default, but can be reset by metdata
real, public :: Zmix_ref = 50.0 !height at which concentration above different landuse are considered equal 

!> Namelist controlled: which veg do we want flux-outputs for
!! We will put the filename, and params (SGS, EGS, etc) in
!! the _Params array.
character(len=15), public, save, dimension(20) :: FLUX_VEGS=""
character(len=15), public, save, dimension(20) :: FLUX_IGNORE=""   ! e.g. Water, desert..
character(len=15), public, save, dimension(20) :: VEG_2dGS=""
character(len=99), public, save, dimension(10) :: VEG_2dGS_Params=""
integer, public, save :: nFluxVegs = 0 ! reset in Landuse_mod

! To use external maps of plant functional types we need to
! map between EMEP codes and netcdf file codes
!NOT USED YET.
type(typ_ss), public, save, dimension(NLANDUSEMAX) :: &
  PFT_MAPPINGS=typ_ss('-','-')

real, public, save :: &
  dt_advec = -999.9,   & ! time-step for advection (s), grid resolution dependent
  dt_advec_inv  ! =1/dt_advec

! NTDAY:  Number of 2D O3 to be saved each day (for SOMO)
! 24/NTDAY is the time integration step for SOMO
! large value -> large memory use; too small value -> bad approx. for SOMO
! NB must be choosen:  24*3600/dt_advec <= NTDAY >=3 and
! preferably an integer fraction of 24*3600/dt_advec
integer, public, parameter :: NTDAY = 72

!/-- choose temperature range: from 148 K (-125C) ro 333K (+60C).
integer, parameter, public :: CHEMTMIN=148,CHEMTMAX=333

real, public, parameter :: &
  V_RAIN     = 5.           !approximate vertical speed of rain m/s

real, public, parameter :: &
  CW_THRESHOLD = 1.0E-7&!Cloudwater (kg/kg); above threshold allow possibility
                        ! for precipitations. Value could be adjusted.
, RH_THRESHOLD = 0.85  &!Relative humidity (fraction); above threshold allow
                        !possibility for precipitations.Value could be adjusted.
, CW2CC = 1.0E6         !Converts Cloudwater (kg/kg) into CloudCover in %
                        !Value could be adjusted.
!
!  additional parameters
!
integer, public, save   :: nterm, nmax, nstep &
                         , iyr_trend ! Year specified for say BC changes

character(len=120), public, save :: runlabel1&!SHORT Allows explanatory text
                                  , runlabel2 !LONG  Read in from grun.pl

! Typically, we define as mainly sea when > 50% water, and
! likely_coastal when > 20%. See Landuse_mod
real, public, parameter, dimension(2) ::  SEA_LIMIT = (/ 0.2, 0.5 /)

real, public, parameter :: &
  MINCONC  = 1.0e-25   ! Experimental. To avoid problems with miniscule numbers
                       ! so that DEBUG runs have same values as normal. Uses if
                       ! USES%MINCONC = T, set in config

real, public, parameter :: &
  EPSIL  = 1.0e-30         &  ! small number
, PASCAL = 100.0           &  ! Conv. from hPa to Pa
, PPB    = 1.0e-9          &  ! parts per billion (mixing ratio)
, PPBINV = 1.0e+9          &
, PPT    = 1.0e-12         &  ! parts per trillion (mixing ratio)
, PPTINV = 1.0e+12         &
, PT     = 1.0e+4          &  ! Top of model region = 10000 Pa = 100 hPa
, TINY   = 1.0e-9             ! -1.E-9" is sometimes used  in order to avoid
                              ! different roundings on different machines.
real, public :: Pref   = 101325.0  ! Reference pressure in Pa used to define vertical levels


! Define output types.
!   Derived output types: types 1..6 (instantaneous,year,month,day,hour,hour_inst),
!                         refer to output variables defined in Derived_mod.
!   Hourly  output types: types 7..8 (hourly_out inst.,hourly_out_mean),
!                         refer to output variables defined in My_Outputs_mod.
integer, public, parameter ::  &
  IOU_INST=1,IOU_YEAR=2,IOU_MON=3,IOU_DAY=4,IOU_HOUR=5,IOU_HOUR_INST=6, & ! Derived output
  IOU_HOUR_EXTRA=7,IOU_HOUR_EXTRA_MEAN=8, & ! additional hourly output
  IOU_MAX_MAX=8                             ! Max values for of IOU (for array declarations)
  !IOU_MAX_MAX hardcoded in OwnDataTypes: please modify consistently!

character, public, parameter ::  & ! output shorthands, order should match IOU_*
  IOU_KEY(IOU_YEAR:IOU_HOUR_INST)=['Y','M','D','H','I']

character(len=*), public, parameter :: model="EMEP_MSC-W "
character(len=TXTLEN_FILE), public :: fileName_O3_Top = "NOTSET"
! Can use RCP values of CH4 for given iyr_trend.
character(len=TXTLEN_FILE), public :: fileName_CH4_ibcs = "NOTSET" ! eg ch4_rcp45

!DSHK logical, parameter, public :: EmisSplit_OUT = .false.
!DSHK logical, parameter, public :: EmisSplit_OUT = .true.

logical, public, parameter:: MANUAL_GRID=.false.!under developement.

!file names
type, public ::names
character(len=TXTLEN_FILE), pointer:: filename => null()
end type names
integer, public, parameter :: Size_InputFiles = 40
type(names), public, save :: InputFiles(Size_InputFiles)

!To add a new filename:
!1) add a line just here below, XXFile = '/default/Path/Default.name'
!2) add the XXFile in NAMELIST /ModelConstants_config/
!3) add a call associate_File(XXFile) near the end of Config_ModelConstants
!4) In the routine using the file, add the XXFile under  "use Config_module"
!5) replace the name you used in the routine with XXFile
character(len=TXTLEN_FILE), target, save, public :: femisFile = 'DataDir/femis.dat'
character(len=TXTLEN_FILE), target, save, public :: Vertical_levelsFile = 'DataDir/Vertical_levels20_EC.txt'
character(len=TXTLEN_FILE), target, save, public :: EmisHeightsFile = 'DataDir/EmisHeights.txt'
character(len=TXTLEN_FILE), target, save, public :: SoilTypesFile = 'DataDir/SoilTypes_IFS.nc'
character(len=TXTLEN_FILE), target, save, public :: SurfacePressureFile = 'DataDir/SurfacePressure.nc'
character(len=TXTLEN_FILE), target, save, public :: AircraftEmis_FLFile = 'DataDir/AircraftEmis_FL.nc'
character(len=TXTLEN_FILE), target, save, public :: nox_emission_1996_2005File = 'DataDir/nox_emission_1996-2005.nc'
!POLL replaced by name of pollutant in EmisSplit
character(len=TXTLEN_FILE), target, save, public :: MonthlyFacFile = 'DataDir/inputs_emepdefaults_Jun2012/MonthlyFac.POLL'
!POLL replaced by name of pollutant in EmisSplit
character(len=TXTLEN_FILE), target, save, public :: DailyFacFile = 'DataDir/inputs_emepdefaults_Jun2012/DailyFac.POLL'
character(len=TXTLEN_FILE), target, save, public :: HourlyFacFile = 'DataDir/inputs_emepdefaults_Jun2012/HourlyFacs.INERIS'
character(len=TXTLEN_FILE), target, save, public :: HourlyFacSpecialsFile = 'NOTSET'
!POLL replaced by name of pollutant in EmisSplit
character(len=TXTLEN_FILE), target, save, public :: SplitDefaultFile = 'DataDir/ZCM_EmChem16x/emissplit_run/emissplit.defaults.POLL'
!POLL replaced by name of pollutant in EmisSplit
character(len=TXTLEN_FILE), target, save, public :: SplitSpecialsFile = 'DataDir/ZCM_EmChem16x/emissplit_run/emissplit.specials.POLL'
character(len=TXTLEN_FILE), target, save, public :: RoadMapFile = 'DataDir/RoadMap.nc'
character(len=TXTLEN_FILE), target, save, public :: AVG_SMI_2005_2010File = 'DataDir/AVG_SMI_2005_2010.nc'
character(len=TXTLEN_FILE), target, save, public :: Soil_TegenFile = 'DataDir/Soil_Tegen.nc'
character(len=TXTLEN_FILE), target, save, public :: SitesFile = 'DataDir/sitesLL.dat'
character(len=TXTLEN_FILE), target, save, public :: SondesFile = 'DataDir/sondesLL.dat'
character(len=TXTLEN_FILE), target, save, public :: GLOBAL_LAInBVOCFile = 'DataDir/GLOBAL_LAInBVOC.nc'
character(len=TXTLEN_FILE), target, save, public :: EMEP_EuroBVOCFile = 'DataDir/LandInputs_Mar2011/EMEP_EuroBVOC.nc'
!SEASON replace by 'jan', 'apr', 'jul' or 'oct' in readdiss
character(len=TXTLEN_FILE), target, save, public :: jclearFile = 'DataDir/jclear.SEASON'
!SEASON replace by 'jan', 'apr', 'jul' or 'oct' in readdiss
character(len=TXTLEN_FILE), target, save, public :: jcl1kmFile = 'DataDir/jcl1.SEASON'
!SEASON replace by 'jan', 'apr', 'jul' or 'oct' in readdiss
character(len=TXTLEN_FILE), target, save, public :: jcl3kmFile = 'DataDir/jcl3.SEASON'
character(len=TXTLEN_FILE), target, save, public :: NdepFile = 'DataDir/AnnualNdep_PS50x_EECCA2005_2009.nc'
!MM replace by month in lightning()
character(len=TXTLEN_FILE), target, save, public :: lightningFile = 'DataDir/lt21-nox.datMM'
character(len=TXTLEN_FILE), target, save, public :: LoganO3File = 'DataDir/Logan_P.nc'
character(len=TXTLEN_FILE), target, save, public :: DustFile = 'DataDir/Dust2014_month.nc'
character(len=TXTLEN_FILE), target, save, public :: TopoFile = 'DataDir/GRID/topography.nc'
character(len=TXTLEN_FILE), target, save, public :: BiDirInputFile = 'NOTSET' ! FUTURE
character(len=TXTLEN_FILE), target, save, public :: Monthly_patternsFile = 'DataDir/ECLIPSEv5_monthly_patterns.nc'

!----------------------------------------------------------------------------
contains
subroutine Config_ModelConstants(iolog)
  integer, intent(in) :: iolog ! for Log file

  integer :: i, j, ispec, iostat
  logical,save :: first_call = .true.
  character(len=len(meteo)) ::  MetDir='./' ! path from meteo
  character(len=*), parameter ::  dtxt='Config_MC:'

  NAMELIST /ModelConstants_config/ &
    DegreeDayFactorsFile, meteo & !meteo template with full path
   ,END_OF_EMEPDAY &
   ,EXP_NAME &  ! e.g. EMEPSTD, FORECAST, TFMM, TodayTest, ....
   ,USES   & ! just testname so far
   ,PBL    & ! Mar2017 testing
   ,EmBio  & ! Mar2017 testing
   ,YieldModifications &  ! Allows dynamic change of chemical yields
   ,LandCoverInputs  & ! Apr2017 for CLM, etc
   ,DEBUG  & !
   ,MY_OUTPUTS  &  ! e.g. EMEPSTD, TFMM
   ,CONVECTION_FACTOR &
   ,EURO_SOILNOX_DEPSCALE &
   ,USE_uEMEP, uEMEP &
   ,INERIS_SNAP1, INERIS_SNAP2 &   ! Used for TFMM time-factors
   ,SELECT_LEVELS_HOURLY, FREQ_HOURLY  & ! incl. 3DPROFILES
   ,ANALYSIS, SOURCE_RECEPTOR, VOLCANO_SR &
   ,SEAFIX_GEA_NEEDED     & ! only if problems, see text above.
   ,BGND_CH4              & ! Can reset background CH4 values
   ,SKIP_RCT              & ! Can  skip some rct
   ,EMIS_OUT, emis_inputlist, EmisDir&
   ,EmisSplit_OUT         & ! Output of species emissions !DSHK
   ,OwnInputDir           &  !
   ,Emis_sourceFiles      & ! new format
   ,EmisMask              & ! new format
   ,USE_SECTORS_NAME      & ! to force a specific sector (SNAP or GNFR)
   ,SecEmisOutWanted      & ! sector emissions to include in output
   ,HourlyEmisOut         & ! to output emissions hourly
   ,FLUX_VEGS             & ! Allows user to add veg categories for eg IAM ouput
   ,FLUX_IGNORE           & ! Specify which landcovers don't need FLUX
   ,VEG_2dGS              & ! Allows 2d maps of growing seasons
   ,VEG_2dGS_Params       & ! Allows 2d maps of growing seasons
   ,PFT_MAPPINGS          & ! Allows use of external LAI maps
   ,NETCDF_DEFLATE_LEVEL,  RUNDOMAIN, DOMAIN_DECOM_MODE &
   ,JUMPOVER29FEB, HOURLYFILE_ending, USE_WRF_MET_NAMES &
   ,dt_advec & ! can be set to override dt_advec
   ,METSTEP &
   ,ZERO_ORDER_ADVEC &! force zero order horizontal and vertical advection 
   ,EUROPEAN_settings & ! The domain covers Europe -> 
   ,GLOBAL_settings & ! The domain cover other regions too -> Convection
   ,fileName_O3_Top&
   ,fileName_CH4_ibcs&
   ,femisFile&
   ,Vertical_levelsFile&
   ,EmisHeightsFile&
   ,SoilTypesFile&
   ,SurfacePressureFile&
   ,AircraftEmis_FLFile&
   ,nox_emission_1996_2005File&
   ,MonthlyFacFile&
   ,DailyFacFile&
   ,HourlyFacFile&
   ,HourlyFacSpecialsFile&
   ,SplitDefaultFile&
   ,SplitSpecialsFile&
   ,RoadMapFile&
   ,AVG_SMI_2005_2010File&
   ,Soil_TegenFile&
   ,SitesFile&
   ,SondesFile&
   ,GLOBAL_LAInBVOCFile&
   ,EMEP_EuroBVOCFile&
   ,jclearFile&
   ,jcl1kmFile&
   ,jcl3kmFile&
   ,NdepFile&
   ,lightningFile&
   ,BiDirInputFile&
   ,LoganO3File&
   ,DustFile&
   ,TopoFile&
   ,Monthly_patternsFile
  NAMELIST /Machine_config/ DataPath

  NAMELIST /INPUT_PARA/GRID,iyr_trend,runlabel1,runlabel2,&
       startdate,enddate!,meteo

  open(IO_NML,file='config_emep.nml',delim='APOSTROPHE')
  read(IO_NML,NML=ModelConstants_config)
  ! do not close(IO_NML), other modules will be read namelist on this file

  USE_SOILNOX = USES%EURO_SOILNOX .or. USES%GLOBAL_SOILNOx

  ! Convert DEBUG%SPEC to index
  if(first_call)then

    if(MasterProc)&
         write(iolog,*) 'CHEM SCHEMES: ',trim(CM_schemes_ChemSpecs)

    ispec = find_index( DEBUG%SPEC, species(:)%name )
  ! print *, "debug%spec testing", ispec, trim(DEBUG%SPEC)
    call CheckStop(ispec<1,"debug%spec not found"//trim(DEBUG%SPEC))
    DEBUG%ISPEC = ispec
    first_call = .false.

  end if

  if(MasterProc)then
    write(*, * ) dtxt//"NAMELIST START "
  ! write(*, NML=ModelConstants_config)
  ! write(*,* ) "NAMELIST IOLOG IS ", iolog
    write(iolog,*) dtxt//"NAMELIST IS "
    write(iolog, NML=ModelConstants_config)
  end if

  DataPath(1) = '.'!default
  rewind(IO_NML)
  read(IO_NML,NML=Machine_config)

  do i=1,size(DataPath)
    if(DataPath(i)=="NOTSET")then
      if(MasterProc)then
        write(*,*)dtxt//'WARNING: Could not find valid DataDir. Tried:'
        do j=1,i-1
          write(*,*)trim(DataPath(j))
        end do
        stop
      end if
      exit
    end if
!   INQUIRE(...) does not behave consistently across intel/gfortran
    open(IO_TMP,file=trim(DataPath(i)),iostat=iostat,action='read')! does not work without action='read'
    if(iostat==0)then
      DataDir=trim(DataPath(i))
      if(MasterProc)write(*,*)dtxt//'DataDir set to',trim(DataDir)
      close(IO_TMP)
      exit
    end if
  end do


  rewind(IO_NML)
  read(IO_NML,NML=INPUT_PARA)

  meteo = key2str(meteo,'DataDir',DataDir)
  meteo = key2str(meteo,'GRID',GRID)
  MetDir= key2str(meteo,'meteoYYYYMMDD.nc','./')
  DegreeDayFactorsFile=key2str(DegreeDayFactorsFile,'MetDir',MetDir)
  DegreeDayFactorsFile=key2str(DegreeDayFactorsFile,'GRID',GRID)
  DegreeDayFactorsFile=key2str(DegreeDayFactorsFile,'YYYY',startdate(1))
  if(MasterProc)then
    write(*,*)dtxt//'Defined DegreeDayFactorsFile as:'
    write(*,*)trim(DegreeDayFactorsFile)
  end if

  if(trim(fileName_O3_Top)/="NOTSET")then
     fileName_O3_Top = key2str(fileName_O3_Top,'DataDir',DataDir)
     fileName_O3_Top = key2str(fileName_O3_Top,'YYYY',startdate(1))
     if(MasterProc)then
        write(*,*)dtxt//'Reading 3 hourly O3 at top from :'
        write(*,*)trim(fileName_O3_Top)
     end if
  endif

  if(trim(fileName_CH4_ibcs)/="NOTSET")then
     fileName_CH4_ibcs = key2str(fileName_CH4_ibcs,'DataDir',DataDir)
     if(MasterProc)then
        write(*,*)dtxt//'Reading CH4 IBCs from:', iyr_trend, trim(fileName_CH4_ibcs)
     end if
  endif

 ! LandCoverInputs
  !print *, dtxt//'Landcover data:', trim(DataDir)
  do i = 1, size(LandCoverInputs%MapFile(:))
    if ( LandCoverInputs%MapFile(i) /= 'NOTSET' ) then
       LandCoverInputs%MapFile(i)= &
          key2str(LandCoverInputs%MapFile(i),'DataDir',DataDir)
!       print *, dtxt//'Landcover file', i, trim(LandCoverInputs%MapFile(i))
       if(MasterProc)then
          write(*,*)dtxt//'Landcover file', i, trim(LandCoverInputs%MapFile(i))
       end if
    end if
  end do
  LandCoverInputs%LandDefs=&
       key2str(LandCoverInputs%LandDefs,'DataDir',DataDir)
  LandCoverInputs%Do3seDefs=&
       key2str(LandCoverInputs%Do3seDefs,'DataDir',DataDir)
  !print *, dtxt//'Landcover =>', LandCoverInputs


  call associate_File(femisFile)
  call associate_File(Vertical_levelsFile)
  call associate_File(EmisHeightsFile)
  call associate_File(SoilTypesFile)
  call associate_File(SurfacePressureFile)
  call associate_File(AircraftEmis_FLFile)
  call associate_File(nox_emission_1996_2005File)
  call associate_File(MonthlyFacFile)
  call associate_File(DailyFacFile)
  call associate_File(HourlyFacFile)
  call associate_File(HourlyFacSpecialsFile)
  call associate_File(SplitDefaultFile)
  call associate_File(SplitSpecialsFile)
  call associate_File(RoadMapFile)
  call associate_File(AVG_SMI_2005_2010File)
  call associate_File(Soil_TegenFile)
  call associate_File(SitesFile)
  call associate_File(SondesFile)
  call associate_File(GLOBAL_LAInBVOCFile)
  call associate_File(EMEP_EuroBVOCFile)
  call associate_File(jclearFile)
  call associate_File(jcl1kmFile)
  call associate_File(jcl3kmFile)
  call associate_File(NdepFile)
  call associate_File(lightningFile)
  call associate_File(BiDirInputFile)  ! FUTURE INPUT
  call associate_File(LoganO3File)
  call associate_File(DustFile)
  call associate_File(TopoFile)
  call associate_File(Monthly_patternsFile)
  do i = 1, size(Emis_sourceFiles)
     !part of a class cannot be a target (?) must therefore do this separately
     if(Emis_sourceFiles(i)%filename/='NOTSET')then
        Emis_sourceFiles(i)%filename = key2str(Emis_sourceFiles(i)%filename,'DataDir',DataDir)
        Emis_sourceFiles(i)%filename = key2str(Emis_sourceFiles(i)%filename,'GRID',GRID)
        Emis_sourceFiles(i)%filename = &
          key2str(Emis_sourceFiles(i)%filename,'OwnInputDir',OwnInputDir)
     endif
  enddo
  do i = 1,size(InputFiles)
    if(associated(InputFiles(i)%filename))then
     InputFiles(i)%filename =key2str(InputFiles(i)%filename,'DataDir',DataDir)
     InputFiles(i)%filename =key2str(InputFiles(i)%filename,'GRID',GRID)
     InputFiles(i)%filename = &
            key2str(InputFiles(i)%filename,'OwnInputDir',OwnInputDir)
    endif
  enddo

  if(.not. USE_uEMEP)NTIMING_uEMEP = 0

end subroutine Config_ModelConstants

subroutine associate_File(FileName)
  integer, save::ix=0
  character(len=*), target ::FileName
  ix = ix+1
  call CheckStop(ix > size(InputFiles) , "Config_module: Size_InputFiles too small")
  InputFiles(ix)%filename => FileName  
end subroutine associate_File

end module Config_module
!_____________________________________________________________________________
!!TSTEMX program testr
!!TSTEMX use Config_module
!!TSTEXM !FAILS due to mpi call Config_ModelConstants()
!!TSTEMX end program testr
