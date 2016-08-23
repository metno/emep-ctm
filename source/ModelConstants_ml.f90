! <ModelConstants_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-201409 met.no
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
module ModelConstants_ml
 !+
 ! Specifies a number of constants used in the model, and reads namelist
 ! file to (re-)configure where possible. 
 ! Note that physical constants (e.g. gravity, Cp, etc ( are specified in
 ! the module PhysicalConstants_ml.f90)
 !
 !----------------------------------------------------------------------------
use CheckStop_ml,         only : CheckStop
use ChemSpecs,            only : species
use Io_Nums_ml,           only : IO_NML, IO_LOG
use PhysicalConstants_ml, only : AVOG
use Precision_ml,         only : dp
use SmallUtils_ml,        only : find_index

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
! variable used to set the type of output in My_outputs_ml.
! The longer term solution puts the outputs into namelists
! but for now we use the MY_OUTPUTS flag. EXP_NAME can
! now be anything descriptive.
CHARACTER(LEN=30), public, save :: EXP_NAME="EMEPSTD"
CHARACTER(LEN=30), public, save :: MY_OUTPUTS="EMEPSTD"

! FORECAST mode run:
! * Nested IC/BC def in Nest_ml & IFSMOZ_ExternalBICs_ml
! * Special hourly output def in My_Outputs_ml
! * Only dayly and hourly output are required,
!   all other output types to false in Derived_ml.
!NML These will be reset by forecast namelist
!NML logical, public, parameter :: FORECAST=&
!NML (EXP_NAME=="FORECAST").or.(EXP_NAME=="EVA2010").or.(EXP_NAME=="EMERGENCY")

! Namelist controlled:
! Some flags for model setup
!------------ NAMELIST VARIABLES - can be reset by emep_namelist.nml file
! USES system introduced June 2013--------------------------

logical, private, parameter :: F = .false.
type, public :: emep_useconfig
  character(len=10) :: testname = "STD"
  logical :: &                   ! Forest fire options
     FOREST_FIRES     = .true.  &! 
    ,MONTHLY_FF       = .false. &! => monthly emissions (for e.g. climate runs)
    ,MACEHEADFIX      = .true.  &! Correction to O3 BCs (Mace Head Obs.)
    ,MACEHEAD_AVG     = .false. &! Uses 10-year avg. Good for e.g. RCA runs.
    ,MINCONC          = .false. &! Experimental. To avoid problems with miniscule numbers
    ,ESX              = .false. &! Uses ESX
    ,EMIS             = .false. &! Uses ESX
    ,EMISSTACKS      = F        !

 ! If USES%EMISTACKS, need to set:
  character(len=4) :: PlumeMethod = "none" !MKPS:"ASME","NILU","PVDI"
end type emep_useconfig 
type(emep_useconfig), public, save :: USES

type, public :: emep_debug
  logical :: &
     AOT             = .false. &
    ,AQUEOUS         = .false. &
    ,BCS             = .false. & ! BoundaryConditions
    ,COLUMN          = .false. & !  Used in Derived_ml for column integratton
    ,DERIVED         = .false. & ! 
    ,DO3SE           = .false. &
    ,DRYRUN          = .false. & ! Skips fast chemistry to save some CPU
    ,EQUIB           = .false. &   !MARS, EQSAM etc.
    ,FORESTFIRE      = .false. &
    ,GLOBBC          = .false. &
    ,GRIDVALUES      = .false. &
    ,HOURLY_OUTPUTS  = .false. & !  
    ,IOPROG          = .false. &
    ,MAINCODE        = .false. & 
    ,MY_DERIVED      = .false. &
    ,pH              = .false. &
    ,PHYCHEM         = .false. &
    ,RUNCHEM         = .false. & ! DEBUG%RUNCHEM is SPECIAL
       ,MY_WETDEP    = .false. &
    ,SETUP_1DCHEM    = .false. &
    ,SETUP_1DBIO     = .false. &
    ,SOLVER          = .false. &
    ,SOA             = .false.  
   integer, dimension(2) ::   IJ = (/ -999, -999 /)
   character(len=20)     ::   SPEC = 'O3'  ! default. 
   integer               ::   ISPEC = -999 ! Will be set after NML
end type emep_debug
type(emep_debug), public, save :: DEBUG


  !/ Emissions file treatment. Dims more than used.
  type, public :: emis_in
    character(len=50) :: name = "NOTSET" ! e.g. POD1_IAM_DF
    integer :: Nlist=0, Nincl=0, Nexcl=0
    character(len=10), dimension(90) ::  incl = "-"
    character(len=10), dimension(90) ::  excl = "-"
  endtype emis_in
  type(emis_in), public, dimension(5) :: emis_inputlist = emis_in()

!-----------------------------------------------------------
logical, public, save ::             &
  FORECAST              = .false.    &! reset in namelist
 ,USE_SOILWATER         = .false.    &!
 ,USE_DEGREEDAY_FACTORS = .false.    &!
 ,USE_SEASALT           = .true.     & !
 ,USE_CONVECTION        = .false.    & ! false works best for Euro runs,
!
! Might sometimes change for scenario runs (e.g. EnsClim):
 ,USE_AIRCRAFT_EMIS  = .true.        & ! Needs global file, see manual
 ,USE_LIGHTNING_EMIS = .true.        & ! 
!
! More experimental:
 ,DO_SAHARA             = .false.    & ! Turn on/off BG Saharan Dust
 ,USE_ROADDUST          = .false.    & ! TNO Road Dust routine. So far with simplified "climate-correction" factor
 ,USE_DUST              = .false.    &! Experimental
 ,TEGEN_DATA            = .false.     &! Interpolate global data to make dust if  USE_DUST=.true.
 ,INERIS_SNAP1       = .false.       & !(EXP_NAME=="TFMM"), & ! Switches off decadal trend
 ,INERIS_SNAP2       = .false.       & !(EXP_NAME=="TFMM"), & ! Allows near-zero summer values
 ,USE_EMERGENCY      = .false.       & ! Emergency: Volcanic Eruption & Nuclear Accident. Under development.
 ,USE_AOD            = .false.       &
 ,USE_POLLEN         = .false.       &  ! EXPERIMENTAL. Only works if start Jan 1
 ,USE_AMINEAQ        = .false.       &  ! MKPS
 ,ANALYSIS           = .false.       &  ! EXPERIMENTAL: 3DVar data assimilation
!
! Output flags
 ,SELECT_LEVELS_HOURLY  = .false.      ! for FORECAST, 3DPROFILES

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
    USE_EURO_SOILNOX      = .true.     & ! ok, but diff for global + Euro runs
   ,USE_GLOBAL_SOILNOX    = .false.    & ! Need to design better switch
   ,USE_SOILNOX           = .true.       ! DO NOT ALTER: Set after config
  real, public, save :: EURO_SOILNOX_DEPSCALE = 1.0 ! 

! Methane background. 
  real, public, save :: BGND_CH4 = -1  ! -1 gives defaults in BoundaryConditions_ml, 
!
!------------ END OF NAMELIST VARIABLES ------------------------------------!

! Some flags for model setup
! will be removed when code is sufficiently tested
! (for convection use foundconv in permanent code)
logical, public, parameter ::         &
  NO_CROPNH3DEP      = .true.,        & ! Stop NH3 deposition for growing crops
  USE_SOILNH3        = .false.,       & ! DUMMY VALUES, DO NOT USE!
  USE_ZREF           = .false.,       & ! testing
  USE_PFT_MAPS       = .false.,       & ! Future option
  EXTENDEDMASSBUDGET = .false.,       & ! extended massbudget outputs
  LANDIFY_MET        = .false.          ! extended massbudget outputs

logical, public ::  USE_EtaCOORDINATES=.false.!temporay parameter; will be set true and removed after testing



!IN-TESTING (reset in NML if wanted)
!Boundary layer profiles
  character(len=4), parameter, public :: FluxPROFILE = &
     "Iter"   !
!      "Ln95"   ! ! will use Launiainen1995  EXPERIMENTAL. Fails in some areas

! Biogenics. Use 3 even if no terpene chemistry - simplifies
! rest of code.  iso = isoprene, mtp = monoterpenes from pools,
! mtl = monoterpenes with light dependence
!DSA12 integer, public, parameter ::   NSOIL_EMIS = 2 ! NO + NH3
 integer, public, parameter ::   NBVOC = 3
 character(len=4),public, save, dimension(NBVOC) :: &
   BVOC_USED = (/ "Eiso","Emt ","Emtl"/)

!The GEA emission data, which is used for EUCAARI runs on the HIRHAM domains
!have in several sea grid cells non-zero emissions in other sectors than SNAP8
!and there are also NH3 emission over sea areas. The former problem makes
!the code crash if the sea areas are defined  as sea (sea=T), so we treat
!them as land in the EUCAARI/HIRHAM runs (sea=F). This is a problem with GEA
!emission data only, not the HIRHAM domain! When e.g. interpolated EMEP emissions
!are used on the HIRHAM domain, this is not a problem.

logical, public, save :: &
   SEAFIX_GEA_NEEDED = .false. ! only if problems. Reset in nml files. NOT HERE

!=============================================================================
!+ 1) Define first dimensions that might change quite often -  for different
!     run domains
character(len=*), parameter, public :: &
! DomainName = "EMEP-50kmEurope"
  DomainName = "EMEP-50kmEECCA"
! DomainName = "EMEP-1degGLOBAL"
! DomainName = "EMEPCWF-0.25degEurope"
! DomainName = "EMEPCWF-0.20degEurope"
! DomainName = "HIRHAM"

!IN-TESTING (reset in NML if wanted)
!) Emissions. Standard =ascii emislist. CdfFractions possible for INERIS
!  and new cdf emission system in testing. Reset in config_ files
! EMIS_TEST can be merged with EMIS_SOURCE after tests
character(len=20), save, public :: &
   EMIS_SOURCE = "emislist"     &! "emislist" or CdfFractions
  ,EMIS_TEST   = "None"          ! "None" or "CdfSnap" 
Logical , save, public :: &
   EMIS_OUT    = .false.           ! output emissions in separate files (memory demanding)

logical, save, public :: IS_GLOBAL = .false.   !!NML .or.(EXP_NAME=="EMERGENCY")
logical, save, public :: MONTHLY_GRIDEMIS = .false.   !! tmp

integer, public :: KMAX_MID, &  ! Number of points (levels) in vertical
                   KMAX_BND     ! Number of level boundaries (=KMAX_MID+1)
integer, public :: IIFULLDOM,JJFULLDOM!  & SET AUTOMATICALLY BY THE CODE
! IIFULLDOM = 182, JJFULLDOM = 197 ! x,y-Dimensions of full HIRHAM domain
! IIFULLDOM = 170, JJFULLDOM = 133 ! x,y-Dimensions of full EMEP domain
! IIFULLDOM = 132, JJFULLDOM = 159 ! x,y-Dimensions of full EECA domain
! IIFULLDOM = 840, JJFULLDOM = 832 ! x,y-Dimensions of full TNO07 domain     
! IIFULLDOM = 420, JJFULLDOM = 416 ! x,y-Dimensions of full TNO14 domain   
! IIFULLDOM = 210, JJFULLDOM = 208 ! x,y-Dimensions of full TNO28 domain 
! IIFULLDOM = 105, JJFULLDOM = 104 ! x,y-Dimensions of full TNO56 domain
! IIFULLDOM = 360, JJFULLDOM = 180 ! .... full GLOBAL domain
! IIFULLDOM = 201, JJFULLDOM = 161 ! .... full GEMS 0.25 domain
! IIFULLDOM = 301, JJFULLDOM = 221 ! .... full GEMS 0.25 extended domain
! IIFULLDOM = 321, JJFULLDOM = 221 ! .... full MACC 0.20 domain

! The difference between EMEP and EECCA is confusing...
integer, public, parameter :: &
! OFFSET_i=  0, OFFSET_j=  0    ! EMEP or default
  OFFSET_i=-35, OFFSET_j=-11    ! EECCA

integer, public, save, dimension(4) ::   &
!                 x0   x1  y0   y1
  RUNDOMAIN =(/-999,-999,-999,-999/)     ! Set values later
! RUNDOMAIN = (/  1, 182,  1, 197 /)     ! HIRHAM
! RUNDOMAIN = (/  1, 132,  1, 159 /)     ! EECCA = new EMEP domain
!!
!  RUNDOMAIN = (/  1, 100,  1, 100 /)     ! Orig EMEP domain in EECCA (for benchmarks)
! RUNDOMAIN = (/ 40, 210, 12, 184 /)     ! SR TNO28 area
! RUNDOMAIN = (/  1, 210,  1, 208 /)     ! TNO28
! RUNDOMAIN = (/240, 720, 48, 736 /)     ! TNO07 reduced (15W-45E;30N-73N)
! RUNDOMAIN = (/120, 360, 24, 368 /)     ! TNO14 reduced (15W-45E;30N-73N)
! RUNDOMAIN = (/ 60, 180, 12, 184 /)     ! TNO28 reduced (15W-45E;30N-73N)
! RUNDOMAIN = (/ 70, 110, 72, 110 /)     ! TNO28  test
! RUNDOMAIN = (/ 30,  90,  6,  92 /)     ! TNO56 reduced (15W-45E;30N-73N)
! RUNDOMAIN = (/ 60, 180, 12, 184 /)     !  test TNO7 area
!--
! Suggestions, 6th June 2012, for TFMM_RUNS scale-dep -----------------------
! RUNDOMAIN = (/  40, 210, 12, 184 /)     ! TNO28 SR area (25W-60E;30N-73N)
! i.e.  - adds 20 squares west, 30 east to TNO28
! RUNDOMAIN = (/  160, 840, 48, 736 /)    ! TNO07  - add 80, 120
! RUNDOMAIN = (/   80, 420, 24, 368 /)    ! TNO14  - add 40, 60
! RUNDOMAIN = (/   20, 105, 6, 92 /)      ! TNO56  - add 10, 15
!----------------------------------------------------------------------------

! RUNDOMAIN = (/ 36, 167, 12, 122 /)     ! EMEP domain in PARLAM
! RUNDOMAIN = (/  1, 360,  1, 180 /)     ! FULL GLOBAL
! RUNDOMAIN = (/  1, 132,  1, 111 /)     ! EECCA, rep09
! RUNDOMAIN = (/  1, 132,  1, 159 /)     ! EECCA, rep10
! RUNDOMAIN = (/ 20, 167,  1, 122 /)     ! OSPAR/HELCOM domain
! RUNDOMAIN = (/ 18, 169,  1, 124 /)     ! OSPAR/HELCOM domain+borders
! RUNDOMAIN = (/  1, 201,  1, 161 /)     ! EMEP-CWF, GEMS 0.25 domain
! RUNDOMAIN = (/  1, 301, 26, 221 /)     ! EMEP-CWF, GEMS 0.25 extended domain
! RUNDOMAIN = (/  1, 321,  1, 221 /)     ! EMEP-CWF, MACC 0.20 domain
! RUNDOMAIN = (/ 70+OFFSET_i, 90+OFFSET_i, 43+OFFSET_j,  63+OFFSET_j /) ! (UK)
! RUNDOMAIN = (/ 60+OFFSET_i, 86+OFFSET_i, 43+OFFSET_j,  59+OFFSET_j /) ! (UK)
! RUNDOMAIN = (/ 40+OFFSET_i, 96+OFFSET_i, 23+OFFSET_j,  69+OFFSET_j /) ! (UK)
! RUNDOMAIN = (/ 60+OFFSET_i,116+OFFSET_i, 33+OFFSET_j,  69+OFFSET_j /) ! (UK)
! RUNDOMAIN = (/ 85+OFFSET_i,120+OFFSET_i, 55+OFFSET_j,  70+OFFSET_j /) ! (changeable)
! RUNDOMAIN = (/ 85+OFFSET_i,120+OFFSET_i, 15+OFFSET_j,  50+OFFSET_j /) ! (changeable)
! RUNDOMAIN = (/ 75+OFFSET_i,110+OFFSET_i, 45+OFFSET_j,  60+OFFSET_j /) ! (gets Esk)
! RUNDOMAIN = (/ 80+OFFSET_i, 106+OFFSET_i, 33+OFFSET_j,  55+OFFSET_j /) ! (France)
! RUNDOMAIN = (/ 80+OFFSET_i, 106+OFFSET_i, 13+OFFSET_j,  35+OFFSET_j /) ! Southern domain
! RUNDOMAIN = (/ 75+OFFSET_i,110+OFFSET_i, 25+OFFSET_j,  60+OFFSET_j /) ! (gets Esk)

integer, public, save ::  & ! Actual number of processors in longitude, latitude
  NPROCX, NPROCY, NPROC     ! and total. NPROCY must be 2 for GLOBAL runs.

!=============================================================================
!+ 2) Define  debug flags.

! We have one variable, to say if we are on master-processor
! or not: (kept here to avoid too many dependencies for box-model
! codes which don't need Par_ml

logical, public, save ::  MasterProc = .true.
logical, public, save ::  DebugCell  = .false.

! For debugging, we often want to print out for  a specific location
! Set here:

!Apr 2014  ALL MOVED TO CONFIG SYSTEM DEBUG%IJ
! The coordinates given here only apply for the standard EMEP domain
!integer, private, parameter :: &
!  DEBUG_ii= -99, DEBUG_jj= -99 ! none
! DEBUG_ii= 79, DEBUG_jj= 56 ! Eskdalemuir
! DEBUG_ii= 73, DEBUG_jj= 48 ! Mace Head
! DEBUG_ii= 88, DEBUG_jj= 53 ! Sibton
! DEBUG_ii= 88, DEBUG_jj= 53 ! Sibton
! DEBUG_ii= 91, DEBUG_jj= 71 ! Rorvik
! DEBUG_ii= 82, DEBUG_jj= 72 ! Voss, has some snow
! DEBUG_ii=110, DEBUG_jj= 48 ! High Vg!
! DEBUG_ii= 96, DEBUG_jj= 40 ! High VG_SO2_CF!
! DEBUG_ii=111, DEBUG_jj= 54 ! High VG_PMCO_CF!
! DEBUG_ii=101, DEBUG_jj= 51 ! Schauinsland
! DEBUG_ii=103, DEBUG_jj= 50 ! Mid-Europe
! DEBUG_ii= 93, DEBUG_jj= 57 ! Elspeetsche (52d12',5d45') 92.83, 56.64
! DEBUG_ii= 92, DEBUG_jj= 56 ! Cabauw
! DEBUG_ii= 97, DEBUG_jj= 62 ! Waldhof
! DEBUG_ii=116, DEBUG_jj= 63 ! K-Puszta
! DEBUG_ii=102, DEBUG_jj= 48 ! Payerne
! DEBUG_ii= 85, DEBUG_jj= 50 ! Harwell
! DEBUG_ii= 88, DEBUG_jj= 99 ! Harwell TNO TEST
! DEBUG_ii= 93, DEBUG_jj= 47 !  Grignon, France
! DEBUG_ii= 90, DEBUG_jj= 104 !  Wetland, Tundra
! DEBUG_ii= 72-OFFSET_i, DEBUG_jj= 37-OFFSET_j ! biomass burnung, Aug 2003
! DEBUG_ii= 90-OFFSET_i, DEBUG_jj= 27-OFFSET_j ! biomass burnung, Jul 2009
!DEBUG_ii= 58-OFFSET_i, DEBUG_jj= 72-OFFSET_j ! 99% water, SMI problems
!DUST DEBUG_ii= 94-OFFSET_i, DEBUG_jj= 24-OFFSET_j ! 99% water, dust problems
! DEBUG_ii= 85, DEBUG_jj= 35 ! Sea, Bay of Biscay
!DEBUG_ii= 76, DEBUG_jj= 65 ! Sea,  North sea
! DEBUG_ii= 66, DEBUG_jj= 50 ! Sea,  west UK
! DEBUG_ii= 80, DEBUG_jj= 52 ! Irish sea
! DEBUG_ii= 91, DEBUG_jj= 67 ! Tange
! DEBUG_ii=103, DEBUG_jj= 32 ! Prades, SMDge
! DEBUG_ii=128, DEBUG_jj= 13 !  Desert?

!integer, public, parameter :: &
! DEBUG_i= 62, DEBUG_j= 45  ! SEA
! DEBUG_i= 10, DEBUG_j= 140 !NEGSPOD
! DEBUG_i= 70, DEBUG_j= 40 ! Lichtenstein, to test ncc
!  DEBUG_i= DEBUG_II+OFFSET_i, DEBUG_j= DEBUG_JJ+OFFSET_j    ! EMEP/EECCA
!  DEBUG_i= 60, DEBUG_j=  50  ! Bremen
! DEBUG_i= 59, DEBUG_j= 79  ! JCOAST
! DEBUG_i= 48, DEBUG_j= 15  !  BB aug 2006
! DEBUG_i= 9, DEBUG_j= 201                                  ! MACC02
! DEBUG_i= 0, DEBUG_j= 0    ! default

!=============================================================================
! Some flags for model setup

! Debug flag DEBUG_XXX  applied in subroutine XXX
 logical, public, parameter ::    &
   DEBUG_ADV            = .false. &
  ,DEBUG_BIO            = .false. &
  ,DEBUG_BLM            = .false. & ! Produces matrix of differnt Kz and Hmix
  ,DEBUG_DERIVED        = .false. &
    ,DEBUG_COLUMN       = .false. & ! Extra option in Derived
  ,DEBUG_ECOSYSTEMS     = .false. &
  ,DEBUG_EMISSTACKS     = .false. &
  ,DEBUG_Kz             = .false. &
  ,DEBUG_DRYDEP         = .false. &
    ,DEBUG_VDS          = .false. &
    ,DEBUG_MY_DRYDEP    = .false. &
    ,DEBUG_CLOVER       = .false. &
    ,DEBUG_STOFLUX      = .false. &
  ,DEBUG_EMISSIONS      = .false. &
  ,DEBUG_EMISTIMEFACS   = .false. &
  ,DEBUG_GETEMIS        = .false. &
  ,DEBUG_LANDDEFS       = .false. &
  ,DEBUG_LANDUSE        = .false. &
  ,DEBUG_LANDPFTS       = .false. &
  ,DEBUG_LANDIFY        = .false. &
  ,DEBUG_MASS           = .false. &
  ,DEBUG_MET            = .false. &
  ,DEBUG_MOSAICS        = .false. &
  ,DEBUG_NEST           = .false. &
  ,DEBUG_NEST_ICBC      = .false. & ! IFS-MOZART BC
  ,DEBUG_NETCDF         = .false. &
  ,DEBUG_NETCDF_RF      = .false. & ! ReadField_CDF in NetCDF_ml
  ,DEBUG_NH3            = .false. & ! NH3Emis experimental
  ,DEBUG_OUTPUTCHEM     = .false. & ! Output of netcdf results
  ,DEBUG_OUT_HOUR       = .false. & ! Debug Output_hourly.f90
  ,DEBUG_POLLEN         = .false.  &
!MV  ,DEBUG_RUNCHEM        = .false. & ! DEBUG_RUNCHEM is SPECIAL
    ,DEBUG_AEROSOL      = .false. & ! ...needed for intended debugs are to work
    ,DEBUG_DUST           = .false. & ! Skips fast chemistry to save some CPU
    ,DEBUG_ROADDUST     = .false. &
    ,DEBUG_SEASALT      = .false. &
    ,DEBUG_SOA          = .false. &
    ,DEBUG_SUBMET         = .false. &
    ,DEBUG_WETDEP       = .false. &
  ,DEBUG_RSUR           = .false. &
  ,DEBUG_RB             = .false. &
  ,DEBUG_SITES          = .false. &
  ,DEBUG_SOILWATER      = .false. &
  ,DEBUG_SOILNOX        = .false. &
  ,DEBUG_VOLC           = .false. & ! Volcanoes
  ,DEBUG_EMERGENCY      = .false.    ! Emergency: Volcanic Eruption & Nuclear Accident. Under development.

!=============================================================================
! 3)  Source-receptor runs?
! We don't (generally) want daily outputs for SR runs, so in
! Derived_ml, we set all IOU_DAY false if SOURCE_RECPTOR = .true..

logical, public, parameter :: SOURCE_RECEPTOR = .false.

! Compress NetCDF output? (nc4 feature, use -1 for netcdf3 output)
integer, public, save :: NETCDF_DEFLATE_LEVEL=4

!Hourly output in single file or monthly/daily files:
!NB: will not work well by default on Stallo per 14th Feb 2012 because of library bugs!
!Until this is fixed, you must compile with netcdf/4.1.3 and link and run with compiler 12.1.2
character(len=*), public, parameter :: &! ending depeding on date:
! HOURLYFILE_ending="_hour_YYYYMM.nc"   ! MM  -> month (01 .. 12)
! HOURLYFILE_ending="_hour_YYYYMMDD.nc" ! DD  -> day of the month (00 .. 31)
! HOURLYFILE_ending="_hour_YYYYJJJ.nc"  ! JJJ -> the day of the year (001 .. 366)
! HOURLYFILE_ending="+FFF.nc"           ! a new file each forecast hour
  HOURLYFILE_ending="_hour.nc"          ! keep the same for the whole run

! NH3 module as set up originally with U10 from met: kept for safety only.
! Will be replaced by sub.grid calculation of wind in future.
! Keep false until code re-implemented
logical, public, parameter :: NH3_U10 = .false.

!=============================================================================
!+ 4)  Define main model dimensions,  things that will
!       generally only change when switching Met-driver
integer, public, parameter ::  &
!TREEX  NLANDUSEMAX  = 19   &   ! Number of land use types in Inputs.Landuse file
  NLANDUSEMAX  = 30    &    ! Max num land use types in Inputs.Landuse file
, METSTEP      = 3     &    ! time-step of met. (h)
, KTOP         = 1     &    ! K-value at top of domain
, KWINDTOP     = 5     &    ! Define extent needed for wind-speed array
, NMET         = 2     &    ! No. met fields in memory
, KCHEMTOP     = 2     &    ! chemistry not done for k=1
, KCLOUDTOP    = 8     &    ! limit of clouds (for MADE dj ??)
, KUPPER       = 6          ! limit of clouds (for wet dep.)


!Namelist controlled: which veg do we want flux-outputs for
character(len=15), public, save, dimension(20) :: FLUX_VEGS=""
integer, public, save :: nFluxVegs = 0 ! reset in Landuse_ml
    

! EMEP measurements end at 6am, used in  daily averages
integer, public, parameter :: END_OF_EMEPDAY  = 6

real, public, save :: dt_advec     ! time-step for advection (s)
real, public, save :: dt_advec_inv ! =1/dt_advec

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
integer, public, save   :: nterm, nmax, nstep, nprint, nass, nbound &
                         , iyr_trend ! Year specified for say BC changes

integer, public, save , dimension(20)   :: identi   !! ????
integer, public, parameter :: TXTLEN_NAME = 20
character(len=120), public, save :: runlabel1&!SHORT Allows explanatory text
                                  , runlabel2 !LONG  Read in from grun.pl

! Typically, we define as mainly sea when > 50% water, and
! likely_coastal when > 20%. See Landuse_ml
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
, Pref   = 101325.0        &  ! Reference pressure in Pa
, TINY   = 1.0e-9             ! -1.E-9" is sometimes used  in order to avoid
                              ! different roundings on different machines.

real, public, parameter :: &
  ATWAIR = 28.964          & ! Mol. weight of air (Jones, 1992)
, atwS   = 32.             & ! Atomic weight of Sulphur
, atwN   = 14.             & ! Atomic weight of Nitrogen
, atwPM  = 100.

! MFAC replaces earlier use of CHEFAC and ATWAIR - to scale from
! density (roa, kg/m3) to  molecules/cm3
! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)
real, public, parameter :: MFAC = 0.001*AVOG/ATWAIR


! Define output types.
!   Derived output types: First 4 types (instantaneous,year,month,day),
!                         refer to output variables defined in Derived_ml.
!   Hourly  output types: Last 2 types (hourly inst.,hourly mean),
!                         refer to output variables defined in My_Outputs_ml.
!   IOU_YEAR_LASTHH: Auxiliary field for hourly accumulated Derived output
integer, public, parameter ::  &
  IOU_INST=1, IOU_YEAR=2, IOU_MON=3, IOU_DAY=4, & ! Derived output
  IOU_YEAR_LASTHH=5,                          & ! Aux. field
  IOU_HOUR=6, IOU_HOUR_MEAN=7                   & ! Hourly  output
  ,IOU_MAX_MAX=7                                  ! Max values for of IOU (for array declarations)

character(len=*), public, parameter :: model="EMEP_MSC-W"

!----------------------------------------------------------------------------
contains
subroutine Config_ModelConstants(iolog)
    character(len=120)  :: txt
    integer, intent(in) :: iolog ! for Log file
    integer :: i, ispec
    logical :: first_call = .true.

    NAMELIST /ModelConstants_config/ &
      EXP_NAME &  ! e.g. EMEPSTD, FORECAST, TFMM, TodayTest, ....
     ,USES   & ! just testname so far
     ,DEBUG  & !
     ,MY_OUTPUTS  &  ! e.g. EMEPSTD, FORECAST, TFMM 
     ,USE_SOILWATER, USE_DEGREEDAY_FACTORS &
     ,USE_CONVECTION &
     ,USE_AIRCRAFT_EMIS,USE_LIGHTNING_EMIS  &  
     ,DO_SAHARA, USE_ROADDUST, USE_DUST &
     ,USE_EURO_SOILNOX, USE_GLOBAL_SOILNOX, EURO_SOILNOX_DEPSCALE &
     ,USE_SEASALT ,USE_POLLEN &
     ,INERIS_SNAP1, INERIS_SNAP2 &   ! Used for TFMM time-factors
     ,IS_GLOBAL  &     ! Also for EMERGENCY
     ,MONTHLY_GRIDEMIS &     ! was set in Emissions_ml
     ,SELECT_LEVELS_HOURLY &  ! incl. FORECAST, 3DPROFILES
     ,FORECAST, USE_EMERGENCY, ANALYSIS , USE_AOD &
     ,SEAFIX_GEA_NEEDED & ! only if problems, see text above.
     ,BGND_CH4  & ! Can reset background CH4 values 
     ,EMIS_SOURCE, EMIS_TEST, EMIS_OUT & 
     ,emis_inputlist &
     ,FLUX_VEGS  & ! TESTX
     ,NETCDF_DEFLATE_LEVEL,  RUNDOMAIN

    txt = "ok"
    !Can't call check_file due to circularity
    !call check_file('emep_settings.nml', file_exists, needed=.true., errmsg=txt)
    open(IO_NML,file='config_emep.nml',delim='APOSTROPHE')
    read(IO_NML,NML=ModelConstants_config)
    !close(IO_NML)

    USE_SOILNOX = USE_EURO_SOILNOX .or. USE_GLOBAL_SOILNOx

    ! Convert DEBUG%SPEC to index

    if( first_call ) then
        ispec = find_index( DEBUG%SPEC, species(:)%name )
        !print *, "debug%spec testing", ispec, trim(DEBUG%SPEC)
        call CheckStop( ispec < 1, &
              "debug%spec not found"//trim(DEBUG%SPEC))
        DEBUG%ISPEC = ispec
        first_call = .false.
    end if

    if ( MasterProc )  then
     write(*, * ) "NAMELIST IS "
     write(*, NML=ModelConstants_config)
     write(*,* ) "NAMELIST IOLOG IS ", iolog
     write(iolog,*) "NAMELIST IS "
     write(iolog, NML=ModelConstants_config)
    end if

  
end subroutine Config_ModelConstants

end module ModelConstants_ml
!_____________________________________________________________________________
