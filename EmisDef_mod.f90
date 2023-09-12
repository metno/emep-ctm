! <EmisDef_mod.f90 - A component of the EMEP MSC-W Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2023 met.no
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module EmisDef_mod

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
use ChemDims_mod,     only : NEMIS_File
use OwnDataTypes_mod, only : TXTLEN_NAME,TXTLEN_FILE, Emis_id_type, &
                             EmisFile_id_type, Emis_mask_type, &
                             Sector_type

implicit none
private

    !----------------- basic emissions file definitions --------------------!
    !  Here we define the parameters *not* likely to change often           !
    !  between different  model versions - e.g. the size and characteristics!
    !  of emission files and sector splits.                                 !
    !-----------------------------------------------------------------------!

integer i, j

! Sector specific information
integer, parameter, public :: NSECTORS_MAX =  50  ! Max. total number of sector included
integer, save, public :: NSECTORS = 0 ! Number of sectors actually included (<= NSECTORS_MAX)

!Predefine Emission heights distributions
integer, public, parameter :: Emis_Nlevel_MAX=25 ! Max Number of vertical emission levels definable
integer, public :: Emis_Nlevel=8 ! Actual Number of total vertical emission levels (can change according to config settings)
integer, public, parameter :: Emis_Nlevel_pre=8 ! Number of vertical emission levels predefined
integer, public, parameter :: Emis_heights_sec_MAX=25 ! Max Number of vertical emission distributions definable
integer, public :: Emis_heights_sec=8 ! Actual Number of vertical emission distributions (can change according to config settings)
integer, public, parameter :: Emis_heights_sec_pre=8 ! Actual Number of vertical emission distributions predefined

!pressure at top of emission levels
!real, public :: Emis_Plevels_pre(Emis_Nlevel_MAX) = &
!     (/101084.9, 100229.1, 99133.2, 97489.35, 95206.225, 92283.825, 88722.15, &
!     (0.0, i = 1,Emis_Nlevel_MAX-Emis_Nlevel_pre)/)
real, public :: Emis_Zlevels_pre(Emis_Nlevel_MAX) = &
  (/20.0,   50.0,   92.0,  184.0,  324.0,  522.0,  781.0, 1106.0, (0.0, i = 1,Emis_Nlevel_MAX-Emis_Nlevel_pre)/)
real, public :: Emis_Zlevels(Emis_Nlevel_MAX) ! used if set by config_emep.nml

!Fraction released in each vertical level predefined values:
real, public :: Emis_h_pre(Emis_heights_sec_MAX,Emis_Nlevel_MAX) = &
     (reshape((/ &
   0.000,  0.000,  0.000,  0.003,  0.147,  0.400,  0.300,  0.150, (0.0, i=1,Emis_Nlevel_MAX-Emis_Nlevel_pre),  &  ! 1 High
   1.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, (0.0, i=1,Emis_Nlevel_MAX-Emis_Nlevel_pre),  &  ! 2 Surface
  0.0600,  0.067,  0.093,  0.750,  0.030,  0.000,  0.000,  0.000, (0.0, i=1,Emis_Nlevel_MAX-Emis_Nlevel_pre),  &  ! 3 = SNAP3
  0.0500,  0.063,  0.087,  0.700,  0.100,  0.000,  0.000,  0.000, (0.0, i=1,Emis_Nlevel_MAX-Emis_Nlevel_pre),  &  ! 4 = SNAP4
  0.0200,  0.034,  0.046,  0.600,  0.300,  0.000,  0.000,  0.000, (0.0, i=1,Emis_Nlevel_MAX-Emis_Nlevel_pre),  &  ! 5 = SNAP5
   0.000,  0.000,  0.000,  0.410,  0.570,  0.020,  0.000,  0.000, (0.0, i=1,Emis_Nlevel_MAX-Emis_Nlevel_pre),  &  ! 6 = SNAP9
   0.200,  0.300,  0.020,  0.044,  0.066,  0.094,  0.123,  0.153, (0.0, i=1,Emis_Nlevel_MAX-Emis_Nlevel_pre),  &  ! 7 = LTO
   0.200,  0.800,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, (0.0, i=1,Emis_Nlevel_MAX-Emis_Nlevel_pre),  &  ! 8 = Shipping
    (0.0, i=1,Emis_Nlevel_MAX*(Emis_heights_sec_MAX-Emis_heights_sec_pre))/), shape(Emis_h_pre)))

real, save, public :: Emis_h(Emis_heights_sec_MAX,Emis_Nlevel_MAX) ! used if set by config_emep.nml

 ! predefine SNAP sectors from EMEP/CORINAIR (NB: splits indices are supposed to be generated from GNFR setup)
   integer, public, parameter :: &
          NSECTORS_SNAP  = 11    ! Number of sectors defined in GNFR emissions
   type(Sector_type), public :: SNAP_SECTORS(NSECTORS_SNAP) = &
        (/ &
!        general name,  longname, netcdf_name timefac index, height index, split index, description, species
   Sector_type('SNAP', 'SNAP1 ',  'sec01',      1,           1,             1,       'Combustion in energy and transformation industries', 'ALL'),&
   Sector_type('SNAP', 'SNAP2 ',  'sec02',      2,           2,             2,       'Non-industrial combustion', 'ALL'),&
   Sector_type('SNAP', 'SNAP3 ',  'sec03',      3,           3,             3,       'Industrial combustion', 'ALL'),&
   Sector_type('SNAP', 'SNAP4 ',  'sec04',      4,           4,             4,       'Production processes', 'ALL'),&
   Sector_type('SNAP', 'SNAP5 ',  'sec05',      5,           5,             5,       'Extraction and distribution of fossil fuels', 'ALL'),&
   Sector_type('SNAP', 'SNAP6 ',  'sec06',      6,           2,             6,       'Solvent use', 'ALL'),&
   Sector_type('SNAP', 'SNAP7 ',  'sec07',      7,           2,             7,       'Road traffic', 'ALL'),&
   Sector_type('SNAP', 'SNAP8 ',  'sec08',      8,           8,             8,       'Other mobile sources (trains, planes, ships)', 'ALL'),&
   Sector_type('SNAP', 'SNAP9 ',  'sec09',      9,           6,             9,       'Waste treatment and disposal', 'ALL'),&
   Sector_type('SNAP', 'SNAP10',  'sec10',     10,           2,            10,       'Agriculture', 'ALL'),&
   Sector_type('SNAP', 'SNAP11',  'sec11',     11,           2,            11,       'Nature', 'ALL') /)

! predefine GNFR CAMS sectors
   integer, public, parameter :: &
          NSECTORS_GNFR_CAMS  = 19    ! Number of sectors defined in GNFR CAMS emissions
   type(Sector_type), public, parameter:: GNFR_CAMS_SECTORS(NSECTORS_GNFR_CAMS ) = & ! mapping between sector names and what they mean
        (/ &
!           general name,  longname, netcdf_name timefac index, height index, split index, description, species
   Sector_type('GNFR_CAMS', 'GNFR_A',  'sec01',      1,            1,            1,       'Public Power', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_B',  'sec02',      3,            3,            2,       'Industry', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_C',  'sec03',      2,            2,            3,       'OtherStationaryComb', 'ALL'),& !NB: used for domestic/degree-day
   Sector_type('GNFR_CAMS', 'GNFR_D',  'sec04',      5,            5,            4,       'Fugitive', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_E',  'sec05',      6,            2,            5,       'Solvents', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_F',  'sec06',      7,            2,            6,       'RoadTransport', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_G',  'sec07',      8,            8,            7,       'Shipping', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_H',  'sec08',      8,            7,            8,       'Aviation', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_I',  'sec09',      8,            2,            9,       'Offroad', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_J',  'sec10',      9,            6,           10,       'Waste', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_K',  'sec11',     10,            2,           11,       'AgriLivestock', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_L',  'sec12',     10,            2,           12,       'AgriOther', 'ALL'),&
   Sector_type('GNFR_CAMS', 'GNFR_M',  'sec13',      5,            5,           13,       'Other', 'ALL') ,&
   !additional subsectors defined for CAMS
   Sector_type('GNFR_CAMS', 'GNFR_A1',  'sec14',      1,            1,           1,       'PublicPower_Point', 'ALL') ,&
   Sector_type('GNFR_CAMS', 'GNFR_A2',  'sec15',      1,            3,           1,       'PublicPower_Area', 'ALL') ,&
   Sector_type('GNFR_CAMS', 'GNFR_F1',  'sec16',      7,            2,          16,       'RoadTransportExhaustGasoline', 'ALL') ,&
   Sector_type('GNFR_CAMS', 'GNFR_F2',  'sec17',      7,            2,          17,       'RoadTransportExhaustDiesel', 'ALL') ,&
   Sector_type('GNFR_CAMS', 'GNFR_F3',  'sec18',      7,            2,          18,       'RoadTransportExhaustLPGgas', 'ALL') ,&
   Sector_type('GNFR_CAMS', 'GNFR_F4',  'sec19',      7,            2,          19,       'RoadTransportNonExhaustOther', 'ALL')/)

   type(Sector_type), public:: SECTORS(NSECTORS_MAX) ! all sectors used during the run (set by model)

   integer, public, parameter :: NCMAX  =  14  ! Max. No. countries per grid point


   integer, save, public :: & !must be compatible with:
      ! timefac, fac_ehh24x7, fac_edd, fac_emm, fac_min, GridTfac, TFAC_IDX_DOM, TFAC_IDX_TRAF
          N_TFAC  = 0  ! Actual number of timefactor classes defined
   integer, save, public :: & !must be compatible with:
          N_HFAC  = 0  ! Actual Nnumber of height distribution classes defined
   integer, save, public :: & !must be compatible with: emisfrac
          N_SPLITMAX  = 30  ! Max number of speciation classes defined
   integer,  save, public :: & !must be compatible with: emisfrac
          N_SPLIT = 0  ! Actual number of speciation classes defined as defined in SplitDefaultFile
   integer, save, pointer, dimension(:), public :: sec2split_map => null()! mapping of sector to speciation class


!The indices defined here are the one used for timefactors indices
   integer, public, save :: &
          TFAC_IDX_POW,   &   ! Power for winter/summer change
          TFAC_IDX_DOM,   &   ! Domestic/residential, for degree-day Timefactors
          TFAC_IDX_AGR,   &   !
          TFAC_IDX_TRAF       !

! Special sectors that need special timefactors corrections. Could be set by config?
 logical, public, save:: IS_POW(NSECTORS_MAX) = .false.
 logical, public, save:: IS_AGR(NSECTORS_MAX) = .false.
 logical, public, save:: IS_TRAF(NSECTORS_MAX) = .false.
 logical, public, save:: IS_DOM(NSECTORS_MAX) = .false.
 logical, public, save:: IS_IND(NSECTORS_MAX) = .false.
   
   !Dust
!   integer, public, parameter ::  NDU   = 2 &   ! number of dust size modes
!                                 ,QDUFI = 1 &   ! production of fine dust
!                                 ,QDUCO = 2     ! production of coarse dust

   !Road Dust
   integer, public, parameter ::  NROADDUST   = 2 &   ! number of road dust size modes
                                 ,QROADDUST_FI = 1 &   ! production of fine road dust
                                 ,QROADDUST_CO = 2     ! production of coarse road dust
   real, public, parameter    ::  ROADDUST_FINE_FRAC = 0.1 ! PM2.5 fraction of PM10-road dust emission

   !Pollen
!  integer, public, parameter ::  NPOL  = 1 &   ! number of dust size modes
!                                ,QPOL  = 1

   !Volcanos.
   logical, public, parameter :: VOLCANOES_LL  = .true.  ! Read Volcanoes
                                                         ! from VolcanoesLL.dat
                                                         ! and disregard them
                                                         ! from gridSOx

!
real, public, save,  allocatable,dimension(:,:) ::  sumcdfemis ! Only used by MasterProc
real, allocatable, public, save,  dimension(:,:) :: cdfemis
real, allocatable, public, save,  dimension(:,:,:) :: Emis_field
real, allocatable, public, save,  dimension(:,:,:) :: Emis_CO_Profile  ! Forest fire testing
integer,  public, save :: NEmis_id
integer,  public :: NEmisMask = 0 !number of masks defined (new format)
real,  public, allocatable :: EmisMaskValues(:,:,:) ! size will be (LIMAX,LJMAX,NEmisMask)
type(Emis_id_type), public, save:: Emis_id(50)
type(EmisFile_id_type), public, save:: EmisFiles(50) !list of emission files after validation
integer, public, save, allocatable, dimension(:,:):: EmisMaskIntVal
integer, public, parameter :: NEmis_sourcesMAX = 50000 ! max number of sources
type(Emis_id_type), public, save :: Emis_source(NEmis_sourcesMAX) ! list of valid sources found in the emission files
integer,  public, save :: NEmis_sources = 0
integer,  public, save :: NEmis_3Dsources = 0
integer,  public, save :: NEmisFile_sources = 0
integer,  public, save :: ix3Dmap(NEmis_sourcesMAX) = 0
real, allocatable, public, save,  dimension(:,:):: Emis_source_ij !source value for a given ij and index
integer, allocatable, public, save,  dimension(:):: NEmis_source_ij !number of valid source values for a given ij
integer, allocatable, public, save,  dimension(:,:):: Emis_source_ij_ix !source index, for a given ij and index
integer, public, parameter :: NEmis_source_ijMAX = 500 ! max number of non zero sources at any ij
real, allocatable, public, save,  dimension(:,:,:,:):: Emis_source_3D !One 3D map for each source
integer, allocatable, public, save,  dimension(:,:,:):: Emis_country_map !country indices for each gridcell
!type(Emis_id_type), public, save:: Emis_source(10)
integer, allocatable, public, save,  dimension(:,:) :: nGridEmisCodes
integer, allocatable, public, save,  dimension(:,:,:):: GridEmisCodes
real, allocatable, public, save,  dimension(:,:,:,:,:):: GridEmis !yearly sector emissions
logical, public, save,  allocatable,dimension(:,:) :: Emis_mask !old format
logical,  public, save :: Emis_mask_allocate = .false.
real,  public, parameter :: MASK_LIMIT = 1.0E-20
! land-code information in each grid square - needed to know which country
! is emitting.
! nlandcode = No. countries in grid square
! landcode  = Country codes for that grid square
integer, public, save, allocatable, dimension(:,:)   :: nlandcode
integer, public, save, allocatable, dimension(:,:,:) :: landcode
! for road dust emission potentials:
integer, public, save, allocatable, dimension(:,:)   :: road_nlandcode
integer, public, save, allocatable, dimension(:,:,:) :: road_landcode
! Emissions for input to chemistry routines
! KEMISTOP added to avoid hard-coded KMAX_MID-3:
!integer, public, parameter :: KEMISTOP = KMAX_MID - NEMISLAYERS + 1
real, public, allocatable, save, dimension(:,:,:,:) :: &
  gridrcemis        ! varies every hour
real, public, allocatable, save, dimension(:,:,:) :: &
  gridrcroadd,    & ! Road dust emissions
  gridrcroadd0      ! varies every hour

!
! The output emission matrix for the 11-SNAP data is secemis:
!
real, public, allocatable, dimension(:,:,:,:,:), save :: &
  secemis      ! main emission arrays, in kg/m2/s

real, public, allocatable, dimension(:,:,:,:), save :: &
! Not sure if it is really necessary to keep the country info; gives rather messy code but consistent with the rest at least (and can do the seasonal scaling for Nordic countries in the code instead of as preprocessing)
  roaddust_emis_pot ! main road dust emission potential arrays, in kg/m2/s (to be scaled!)

! We store the emissions for output to d_2d files and netcdf in kg/m2/s
real, public, allocatable, dimension(:,:,:), save :: EmisOut!per emitted species
real, public, allocatable, dimension(:,:,:), save :: SplitEmisOut!per splitted species
real, public, allocatable, dimension(:,:,:,:), save :: SecEmisOut !per sector and species

character(len=TXTLEN_NAME), public, save :: mask2name(1000) = 'NOTSET' !name of mask id number

!Ocean variables
type, public :: Ocean
  real,  allocatable, dimension(:,:) :: emis
  real,  allocatable, dimension(:,:) :: map
  real :: sum_month
  real :: sum_year
  integer :: index
end type Ocean

type(Ocean), public, save:: O_NH3, O_DMS

!used for EEMEP
real, allocatable, save, dimension(:,:,:,:), public       ::  Emis_4D !(i,j,k,pollutant)
integer, save, public ::N_Emis_4D=0 !number of pollutants to read
integer, public, save :: Found_Emis_4D = 0

integer, public, save :: KEMISTOP ! not defined yet= KMAX_MID - nemis_kprofile + 1

 ! Names of emis files, generated by GenChem:

     include 'CM_EmisFile.inc'

 ! Road dust emission files (should perhaps be added to GenChem in the future?)

  integer, parameter, public :: NROAD_FILES = 2
  character(len=11), save, dimension(NROAD_FILES), public:: &
       ROAD_FILE =  (/ "HIGHWAYplus", "NONHIGHWAY " /)
  character(len=21), save, public:: &
       ROADDUST_CLIMATE_FILE =  "ROADDUST_CLIMATE_FAC"


 ! FUTURE work
 ! NMR-NH3 project specific variables
 ! NH3 emissions set by meteorology and special activity data
    logical, public, parameter :: NH3EMIS_VAR = .false.
    real, public, save  :: dknh3_agr ! reported nh3emis (IC_NMR)
                                     ! read from gridXXfile

  integer, parameter, public :: MAXFEMISLONLAT = 10!max number of lines with lonlat reductions
  integer,   public          :: N_femis_lonlat    !number of femis lonlat lines defined
  !DS allow femis_lonlat to apply internal or external
  logical, public, dimension(MAXFEMISLONLAT) :: femis_lonlat_internal = .true.


  integer, public, save :: NSecEmisOutWanted = 0 !sum of all sectors not included in this N
  integer, public, allocatable, save :: isec2SecOutWanted(:)

  logical, public, save :: foundYearlySectorEmissions = .false.
  logical, public, save :: foundMonthlySectorEmissions = .false.

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
                     end module EmisDef_mod
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!_____________________________________________________________________________
