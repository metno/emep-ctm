! <OwnDataTypes_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.6>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2025 met.no
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
! <OwnDataTypes_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************!
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!
module OwnDataTypes_mod
use Country_mod,     only : MAXNLAND
use NumberConstants, only : UNDEF_I, UNDEF_R
use TimeDate_mod,    only : date

implicit none
private


public :: print_Deriv_type
public :: print_Asc2D
public :: print_Sector_type
integer, public, parameter :: &
  TXTLEN_DERIV = 44, &
  TXTLEN_SHORT = 28, &
  TXTLEN_IND   =  6, &
  TXTLEN_NAME =  64, & !for performance, should be a multiple of 8
  TXTLEN_FILE = 200    ! large enough for paths from namelists

! Contains some user-defined data-types, and routine associated
! with these. Collecting them here will
! avoid some dependencies, and shorten some My type modules.


  !==================
  !/ generic groups for integers
  type, public :: typ_i2
    integer :: int1
    integer :: int2
  end type typ_i2

  type, public :: typ_i3
    integer :: int1
    integer :: int2
    integer :: int3
  end type typ_i3

  !/ generic group for two (short) strings
  type, public :: typ_ss
    character(len=TXTLEN_SHORT) :: txt1='-' ! e.g. POD1_IAM_DF
    character(len=TXTLEN_SHORT) :: txt2='-' ! e.g. POD1_IAM_DF
  end type typ_ss

  !/ generic group for two (short) strings and float
  !  currently for CMX boundary conditions
  type, public :: typ_ssf
    character(len=TXTLEN_SHORT) :: txt1='-' ! e.g. POD1_IAM_DF
    character(len=TXTLEN_SHORT) :: txt2='-' ! e.g. POD1_IAM_DF
    real :: num  = UNDEF_R
  end type typ_ssf

 !/ generic group for name and pointer to arrays
  type, public :: typ_sp
    character(len=TXTLEN_SHORT) :: name ! e.g. POD1_IAM_DF
    integer, dimension(:), pointer :: specs
  end type typ_sp

  !/ HI: generic group for name and two pointers to one integer and one
  !/  real array
  type, public :: typ_factors
    character(len=TXTLEN_SHORT) :: name ! e.g. POD1_IAM_DF
    integer, dimension(:), pointer :: species ! like ptr in typ_sp
    real, dimension(:), pointer :: factors
  end type typ_factors

  !/ HI: generic group for name and two pointers to one integer and one
  !/  character array
  type, public :: typ_maps
    character(len=TXTLEN_SHORT) :: name ! e.g. POD1_IAM_DF
    integer, dimension(:), pointer :: species ! like ptr in typ_sp
    character(len=TXTLEN_SHORT), dimension(:), pointer :: maps ! other
      ! species' or variables' names to map this one to
  endtype typ_maps

  !/ generic group one (short) string & one integer
  type, public :: typ_si
    character(len=TXTLEN_SHORT) :: name = '-'
    integer :: ind = UNDEF_I
  endtype typ_si

  !/ generic group for one (short) string & one shorter string
  type, public :: typ_s1ind
    character(len=TXTLEN_SHORT) :: name
    character(len=TXTLEN_IND)   :: ind  ! e.g. YMDHI
  end type typ_s1ind

  !/ generic group one (short) string & one real
  type, public :: typ_sr
    character(len=TXTLEN_SHORT) :: name
    real :: rval
  end type typ_sr

!/ generic group for three (short) strings
type, public :: typ_s3
  character(len=TXTLEN_SHORT) :: txt1,txt2,txt3
end type typ_s3

!/ generic group for four (short) strings
type, public :: typ_s4
  character(len=TXTLEN_SHORT) :: txt1,txt2,txt3,txt4 ! e.g. POD1_IAM_DF
end type typ_s4

!/ generic group for five (short) strings & one integer
type, public :: typ_s5i
  character(len=TXTLEN_SHORT) :: txt1,txt2,txt3,txt4, &
                                 txt5 ! e.g. SO2,ugS,2d,AIR_CONC,SPEC
  integer                     :: ind  ! e.g. IOU_DAY
end type typ_s5i
!/ generic group for five (short) strings & one shorter string
type, public :: typ_s5ind
  character(len=TXTLEN_SHORT) :: txt1,txt2,txt3,txt4, &
                                 txt5 ! e.g. SO2,ugS,2d,AIR_CONC,SPEC,
  character(len=TXTLEN_IND)   :: ind  ! e.g. YMDHI
end type typ_s5ind

!==================
!+ Derived output type
type, public:: Deriv
  character(len=TXTLEN_DERIV) :: name     = '-' ! e.g. DDEP_SO2_m2Conif
  character(len=TXTLEN_SHORT) :: class    = '-' ! Type of data, e.g. ADV or Mosaic
  character(len=TXTLEN_SHORT) :: subclass = '-' ! e.g. "VG", "Rns"
  character(len=TXTLEN_SHORT) :: txt      = '-' ! text where needed, e.g. "Conif"
  character(len=TXTLEN_SHORT) :: unit     = '-' ! writen in netCDF output
  integer :: index         =UNDEF_I ! index in concentation array, or other
  integer :: f2d           =UNDEF_I ! index in f_2d arrays
  logical :: dt_scale      =.false. ! used only if we need a factor on dt_advec,
  real    :: scale         =UNDEF_R ! e.g. use 100.0 to get cm/s
  logical :: avg           =.true.  ! True => average data (divide by nav at end),
                                    ! else accumulate over run period
  character(len=TXTLEN_IND)   :: iotype   = '-' ! sets output timing
end type

! Sentinel values (moved to NumberConstants)
! real,    private, parameter :: UNDEF_R = -huge(0.0)
! integer, private, parameter :: UNDEF_I = -huge(0)

!==================
!+ Hourly ASCII/NetCDF output type
type, public:: Asc2D
  character(len=TXTLEN_DERIV):: name = "-"   ! Name (no spaces!)
  character(len=TXTLEN_SHORT):: type = "-"  ! "ADVppbv" or "ADVugm3" or "SHLmcm3"
! character(len=9) :: ofmt      ! Output format (e.g. es12.4)
  integer          :: spec = UNDEF_I   ! Species number in xn_adv or xn_shl array
                                ! or other arrays
  integer          :: nk = UNDEF_I     ! number of vertical levels
  character(len=TXTLEN_SHORT) :: unit   ! Unit used
  real             :: unitconv = UNDEF_R  !  conv. factor
  real             :: max      = UNDEF_R        ! Max allowed value for output
end type

!==================
!+ Defines SOA, NONVOL and VBS params
type, public :: VBST
  integer     :: index    ! just for clarity
  real        :: CiStar   ! ug/m3
 !real        :: Tref     ! Assumed 300
  real        :: DeltaH   ! kJ/mole
end type VBST
!==================


type, public:: O3cl_t
   character(len=TXTLEN_DERIV) :: name = '-'  ! e.g. POD1_IAM_DF
   character(len=TXTLEN_SHORT) :: class = '-' ! POD or AOT
   real    :: Threshold = UNDEF_R     ! Threshold or CL, e.f. AOTx or AFstY
   character(len=TXTLEN_SHORT) :: defn = '-'  !  MM or EU definitions
   character(len=TXTLEN_SHORT) :: txtLC = '-' !  CF, DF, IAM_CF etc.
   logical :: RelSGS = .false.  ! true if accumulation period is relative to
   ! start of growing season (SGS)
   ! can be false it fixed, e.g. April 1st
   integer :: SAccPeriod = UNDEF_I  ! Start of accumulation period, either
   ! rel to SGS or day number, days
   integer :: EAccPeriod = UNDEF_I  ! End ...., days
   character(len=TXTLEN_IND)            :: iotype = '-'! .. 'M'=>IOU_MON, 'D'=>IOU_DAY, ...
end type O3cl_t

!Sector definitions

type, public :: Sector_type
   character(len=TXTLEN_NAME) :: name    = 'NOTSET' ! general name of sector type (GNFR, SNAP ...)
   character(len=TXTLEN_NAME) :: longname = 'NOTSET'! specific name of sector
   character(len=TXTLEN_NAME) :: cdfname = 'NOTSET' ! variable name as used in the netcdf file
   integer                    :: timefac = -1       ! identification number in the timefactor input files
   integer                    :: height  = -1       ! identification number in the Emission height definition
   integer                    :: split   = -1       ! identification number in the splits input files
   character(len=TXTLEN_NAME) :: description = 'NOTSET' ! description in human language
   character(len=TXTLEN_NAME) :: species = 'ALL'    ! which species to include (ALL means loop through all emitted species)
end type Sector_type


type, public :: Emis_id_type
   character(len=TXTLEN_NAME) :: varname = 'NOTSET' !name of variable in netcdf file
   character(len=TXTLEN_NAME) :: species = 'NOTSET' !which species to put emissions into
   character(len=TXTLEN_NAME) :: units = 'NOTSET'!units AFTER netcdf values are multiplied by factor
   character(len=TXTLEN_NAME) :: country_ISO = 'NOTSET' !country name, for example FR for France, as defined in Country_mod
   character(len=TXTLEN_NAME) :: periodicity = 'NOTSET' !how often fresh values must be read from the netcdf file
   character(len=TXTLEN_NAME) :: timevalidity = 'NOTSET' !if the time refers to the start, middle or end of the period
   integer :: countrycode = -1 ! number identifying country in the emission file
   integer :: sector = -1 !sector as defined in this file
   integer :: sector_idx = -1 ! internal index used in SECTORS (set by model)
   integer :: species_ix = -1 ! internal index for species
   integer :: injection_k = -1 !which model k level to put emissions into. Only for individual species
   real    :: factor = -1.0 ! scaling factor. multiply values by this number
   logical :: include_in_local_fractions = .true. !if this is to be accounted in the local fractions (uEMEP)
   logical :: apply_femis = .true. !whether the general femis.dat should be applied to this source
   character(len=TXTLEN_NAME) :: mask_ID = 'NOTSET' ! set to ID of mask, if to be applied
   character(len=TXTLEN_NAME) :: mask_ID_reverse = 'NOTSET' ! set to ID of mask, if to be applied as reversed
   integer :: mask_ix = -1 ! mask index, >0 if set. Internal index, do not set
   integer :: mask_reverse_ix = -1 ! mask index, >0 if set. Internal index, do not set
   integer :: country_ix = 67 !Internal country index. Does not have any meaning outside of code
   integer :: height = 0 !could define own release height. not implemented
   logical :: is3D = .false.
   integer :: istart = -999
   integer :: jstart = -999
   integer :: kstart = -1
   integer :: kend = -1
   logical :: reversek = .true.
   integer :: ix_in = -1!index of the corresponding source in the config defintions (internal use only)
end type Emis_id_type

integer, parameter, public :: NSOURCESMAX = 50
type, public :: Emis_sourceFile_id_type
   character(len=TXTLEN_FILE) :: filename = 'NOTSET'!netcdf filename with path
   character(len=TXTLEN_NAME) :: projection = 'NOTSET' !projection or 'native' if same projection and size as meteo grid
   character(len=TXTLEN_NAME) :: periodicity = 'NOTSET' !how often fresh values must be read from the netcdf file
   character(len=TXTLEN_NAME) :: timevalidity = 'end' !if the time refers to the start, middle or end of the period
   real                       :: grid_resolution = 0.0 !resolution of the emission file
   real :: factor = -1.0 !scaling factor. multiply values for all sources by this number. Comes on top of source factors.
   type(Emis_id_type) :: source(NSOURCESMAX) ! one source defined for each netcdf field to include
!default values for sources:
   character(len=TXTLEN_NAME) :: species = 'NOTSET' !default emep species
   character(len=TXTLEN_NAME) :: units = 'NOTSET'! default units
   character(len=TXTLEN_NAME) :: country_ISO = 'NOTSET' ! default country name
   character(len=TXTLEN_NAME) :: sectorsName = 'NOTSET' !
   integer :: countrycode = -1 ! number identifying country in the emission file
   integer :: sector = -1 !default sector
   logical :: apply_femis = .true. !whether the general femis.dat should be applied to sources from this file
   logical :: include_in_local_fractions = .true. !if this is to be accounted in the local fractions (uEMEP)
   character(len=TXTLEN_NAME) :: mask_ID = 'NOTSET' ! set to ID of mask, if to be applied. Will then be default for all sources in file
   character(len=TXTLEN_NAME) :: mask_ID_reverse = 'NOTSET' ! set to ID of mask, if to be applied as reversed. Will then be default for all sources in file
   character(len=TXTLEN_NAME) :: country_ISO_excl(MAXNLAND) = 'NOTSET' ! Exclude those countries
   character(len=TXTLEN_NAME) :: country_ISO_incl(MAXNLAND) = 'NOTSET' ! If set, include only those countries
end type Emis_sourceFile_id_type

type, public :: Emis_mask_type
   character(len=TXTLEN_FILE) :: filename = 'NOTSET'! netcdf filename with path
   character(len=TXTLEN_NAME) :: cdfname = 'NOTSET' ! name of the mask in the netcdf file
   character(len=TXTLEN_NAME) :: ID = 'NOTSET' ! name that the user set to identify this mask
   character(len=TXTLEN_NAME) :: type = 'CELL-FRACTION' ! Type of mask: 'NUMBER', 'CELL-FRACTION', 'THRESHOLD'
   real                       :: threshold = 1.E-20 !mask is set for where value is above threshold
   real                       :: threshold_max = 1.E60 !mask is not set if value above threshold
   real                       :: fac = 0.0 !multiplicative factor
end type Emis_mask_type

type, public :: EmisFile_id_type
   character(len=TXTLEN_FILE) :: filename = 'NOTSET'!netcdf filename with path
   character(len=TXTLEN_NAME) :: projection = 'NOTSET' !projection or 'native' if same projection and size as meteo grid
   character(len=TXTLEN_NAME) :: periodicity = 'NOTSET' !how often fresh values must be read from the netcdf file
   character(len=TXTLEN_NAME) :: timevalidity = 'NOTSET' !if the time refers to the start, middle or end of the period
   real                       :: grid_resolution = 0.0!resolution of the emission file
   real :: factor = -1.0 !scaling factor. multiply values for all sources by this number. Comes on top of source factors.
   type(date) :: end_of_validity_date = date(0,0,0,0,0)!internal date to know when to fetch new data
   integer :: Nsources = 0 !number of valid sources (i.e variables in the netcdf file)
   integer :: source_start = 0
   integer :: source_end = 0
   character(len=TXTLEN_NAME) :: species = 'NOTSET' !default emep species. NB: not read from attributes
   character(len=TXTLEN_NAME) :: units = 'NOTSET'! default units
   character(len=TXTLEN_NAME) :: country_ISO = 'NOTSET' ! default country name
   character(len=TXTLEN_NAME) :: sectorsName ='NOTSET' !SNAP or GNFR_CAMS or user defined name
   integer :: countrycode = -1 !default countrycode
   integer :: sector = -1 !default sector
   integer :: nsectors = 1 !Number of sectors stored in each variable. Will only change how many sources are read at each variable read
   character(len=TXTLEN_NAME) :: mask_ID = 'NOTSET' ! set to ID of mask, if to be applied. Will then be default for all sources in file .NB: not read from attributes
   character(len=TXTLEN_NAME) :: mask_ID_reverse = 'NOTSET' ! set to ID of mask, if to be applied as reversed. Will then be default for all sources in file .NB: not read from attributes
   integer :: ncFileID = -1 !internal: shows the netcdf file ID if the file is open, or must be < 0.
end type EmisFile_id_type

type, public :: hourly_emis_factor_type
   character(len=TXTLEN_FILE) :: file = 'NOTSET' !filename with path
   character(len=TXTLEN_NAME) :: poll = 'NOTSET' !one of the emitted pollutants, nox, sox, pm25 etc. Case sensitive
   character(len=TXTLEN_NAME) :: cdfname = 'NOTSET' ! name of the variable in the file
end type hourly_emis_factor_type


!/ Emissions file treatment. Dims more than used.
type, public :: emis_in
  character(len=150) :: name = "NOTSET" ! e.g. POD1_IAM_DF
  integer :: Nincl=0, Nexcl=0
  character(len=10), dimension(90) ::  incl = "-"
  character(len=10), dimension(90) ::  excl = "-"
  character(len=40), dimension(20) ::  pollName = "NOTSET"
  character(len=40), dimension(20) ::  pollemepName = "NOTSET"
  character(len=40) ::  periodicity = "NOTSET" !How often new data should be read in
  character(len=40) ::  format = "NOTSET" !set to fraction, if fractions
  character(len=40) ::  type = "sectors" !steers special treatments
  logical ::  use_lonlat_femis = .true. !allows to switch off lonlat femis reductions
                                        !for specific emission files
                                        !Country+sector specific reductions can be dealt
                                        !with with incl/excl, so those are not affected
  logical :: set_mask = .false.  !if T, set mask for each (i,j) where non zero emission is found
  logical :: use_mask = .false.  !if T, do not include emission where mask is set
  character(len=40) ::  sector = "NOTSET" ! e.g. GNFR_CAMS
  real ::  scale = 1.0 ! multiply by scale (not yet implemented)
end type emis_in

!==================
! Local Fractions parameters
integer, public, parameter :: Max_lf_sources = 1000 !max number of sources to track.
integer, public, parameter :: Max_lf_Country_list = 1000 ! max number of countries for each lf pollutant
integer, public, parameter :: Max_lf_Country_groups = 30
integer, public, parameter :: Max_lf_sectors = 50
integer, public, parameter :: Max_lf_res = 50
integer, public, parameter :: Max_lf_spec = 250
integer, public, parameter :: Max_lf_out = 100
type, public :: poll_type
  character(len=TXTLEN_NAME):: name = 'NOTSET'    ! pollutants to include
  integer, dimension(Max_lf_sectors) ::sectors = -1    ! sectors to be included for this pollutant. Zero is sum of all sectors
  integer, dimension(Max_lf_res) ::res = -1    ! resolution of sources to be included for this pollutant.
end type poll_type

type, public :: lf_set_type
  !general
  character(len=TXTLEN_FILE) :: filename_writeatend = 'LF_saveatendDDMMYYYY.nc'
  character(len=TXTLEN_FILE) :: filename_write = 'LF_saveDDMMYYYY.nc'
  character(len=TXTLEN_FILE) :: filename_read = 'LF_saveDDMMYYYY.nc'
  integer :: Nvert = 14 ! vertical extend of the tracking/local window
  integer :: Nvertout = 1 ! number of vertical level to output (non fullchem only)
  logical :: YEAR =.true.! Output frequency
  logical :: MONTH =.false.
  character(len=40)::  MONTH_ENDING = "NOTSET"
  logical :: DAY =.false.
  logical :: HOUR =.false.
  logical :: HOUR_INST =.false.
  integer, dimension(4) :: DOMAIN = -1 ! DOMAIN which will be outputted
  !for fullchem settings
  integer :: dist = -1
  logical :: full_chem =.false.
  logical :: relative =.false. ! compute also grid to grid values 
  integer :: Nfullchem_emis = -1 ! number of emission types to track: 1 {nox+voc+nh3+sox}, 2 {nox,voc}, 4 {nox,voc,nh3,sox} 
  logical :: EmisDer_all =.false. ! reduce voc, sox, nox, nh3 together. Overwritten if Nfullchem_emis is set
  logical :: MDA8 = .false. ! if MDA8 and SOMO35 are to be outputed (if full_chem)
  logical :: restart =.false.
  logical :: Nestsave =.true. !if nesting and LF are used, save also lf values each time Nest is saving 3D
  logical :: saveatend =.false.
end type lf_set_type


type, public :: lf_sources
  character(len=TXTLEN_NAME) :: species = 'NOTSET' !pollutants to include
  character(len=TXTLEN_NAME) :: type = 'relative' !Qualitatively different type of sources: "relative", "country"
  character(len=TXTLEN_NAME) :: name = 'NOTSET' ! name as it appears in output. Only for "relative" type
  integer :: dist = -1 ! window dimension, if defined (window size is 2*dist+1 x 2*dist+1 )
  integer :: res = 1  ! half size of the single source square (square size is 2*res x 2*res )
  integer :: Nvert = 7 ! vertical extend of the tracking/local rwindow
  integer :: sector = 0 ! sector for this source. Zero is sum of all sectors
  integer :: poll = 1 !index of pollutant. One poll for all sources related to that poll (set by model)
  integer :: start = 1 ! first position index in lf_src (set by model)
  integer :: end = 1 ! last position index in lf_src (set by model)
  integer :: iem = 0 ! index of emitted pollutant, emis (set by model)
  integer :: iem_deriv = 0 ! index of emitted pollutant to track, emis (set by model)
  integer :: iem_lf ! index of emitted internal for LF (1 for nox, 2 for voc, ...)
  integer :: Npos = 0 ! number of position indices in lf_src (set by model)
  integer :: nhour = -1 ! number of hours between timestamps, and resets. Not used if <0. Must be <=24
  integer :: time_ix = 0 ! start of hour at which the emissions are set (set by model)
  integer :: Nsplit = 0 ! into how many species the emitted pollutant is split into (set by model)
  integer :: species_ix = -1 !species index, if single pollutant (for example NO or NO2, instead of nox)
  integer :: iqrc = -1 !index for emissplits, if single pollutant (for example NO or NO2, instead of nox)
  integer, dimension(15) :: ix = -1 ! internal index of the  (splitted) species (set by model)
  real, dimension(15) :: mw=0.0  ! molecular weight of the (splitted) species (set by model)
  character(len=TXTLEN_NAME) :: country_ISO = 'NOTSET' !country name, for example FR for France, as defined in Country_mod
  integer :: country_ix = -1 !Internal country index. Does not have any meaning outside of code
  logical :: DryDep = .false. ! if drydep is to be outputed
  logical :: WetDep = .false. ! if wetdep is to be outputed
  logical     :: YEAR = .true.! Output frequency
  logical     :: MONTH = .false.
  logical     :: make_fracsum = .false.
  character(len=40)::  MONTH_ENDING = "NOTSET"
  logical     :: DAY = .false.
  logical     :: HOUR = .false.
  logical     :: HOUR_INST = .false.
  logical     :: is_ASOA = .false.
  logical     :: is_NATURAL = .false.
end type lf_sources

type, public :: lf_out_type
  character(len=TXTLEN_NAME):: name = "NOTSET"
  character(len=TXTLEN_NAME):: species(30) = "NOTSET"
  real                      :: species_fac(30) = 1.0
  integer                   :: ix(30) = -1 ! internal index in loc_frac_drydep
  logical                   :: DryDep
  logical                   :: WetDep
  logical                   :: relative = .false.
end type lf_out_type

integer, parameter, public :: MAX_lf_country_group_size = 50 !max 50 countries in each group
type, public :: lf_country_group_type
   character(len=TXTLEN_NAME) :: name = 'NOTSET' !the overall name of the group (for example 'EU')
   character(len=10), dimension(MAX_lf_country_group_size):: list = 'NOTSET' ! list of countries inside the group
   integer, dimension(MAX_lf_country_group_size):: ix = -1 ! index of the country as defined in Country_ml (set by model)
end type lf_country_group_type

integer, parameter, public :: MAX_lf_sector_group_size = 30 !max 30 sectors in each group
type, public :: lf_sector_group_type
   character(len=TXTLEN_NAME) :: name = 'NOTSET' !the overall name of the sector group (for example 'Low')
   integer :: nsec = 0 ! number of valid sectors defined in this group (set by model)
   integer :: list(MAX_lf_sector_group_size) = -1! list of sectors inside the group
end type lf_sector_group_type

type, public :: lf_country_type
   integer :: mask_val_min = 1
   integer :: mask_val_max = 0
   integer,dimension(Max_lf_Country_list) :: mask_val = -9999999
   character(len=10) :: list(Max_lf_Country_list) = 'NOTSET'
   character(len=TXTLEN_NAME) :: cellmask_name(Max_lf_Country_list) = 'NOTSET'
   type(lf_country_group_type) :: group(Max_lf_Country_groups)
   integer :: sector_list(Max_lf_sectors)=-1
end type lf_country_type


contains
!=========================================================================
subroutine print_Asc2D(w)
  type(Asc2D), intent(in) :: w  ! wanted
  write(*,*) "Prints Asc2D type ========================="
  write(*,"(a,a)")      "Name   :", trim(w%name)
  write(*,"(a,a)")      "type   :", trim(w%type)
  write(*,"(a,i4)")     "spec   :", w%spec
  write(*,"(a,a)")      "unit   :", trim(w%unit)
  write(*,"(a,i4)")     "nk     :", w%nk
  write(*,"(a,es10.3)") "unitconv:",w%unitconv
  write(*,"(a,es10.3)") "max    :", w%max
end subroutine print_Asc2D
!=========================================================================
subroutine print_Deriv_type(w)
  type(Deriv), intent(in) :: w  ! wanted
  write(*,*) "Prints Deriv type ========================="
  write(*,"(a,a)")      "Name   :", trim(w%name)
  write(*,"(a,a)")      "class  :", trim(w%class)
  write(*,"(a,a)")      "subclass:",trim(w%subclass)
  write(*,"(a,a)")      "txt    :", trim(w%txt)
  write(*,"(a,a)")      "units  :", trim(w%unit)
  if(w%index==UNDEF_I)then
    write(*,"(a)")      "index  : UNDEF"
  else
    write(*,"(a,i5)")   "index  :", w%index
  end if
  if(w%f2d==UNDEF_I)then
    write(*,"(a)")      "f2d    : UNDEF"
  else
    write(*,"(a,i3)")   "f2d    :", w%f2d
  end if
  write(*,"(a,es10.3)") "scale  :", w%scale
  write(*,*)            "dt_scale:", w%dt_scale
  write(*,*)            "avg    :", w%avg
end subroutine print_Deriv_type
!=========================================================================
subroutine print_Sector_type(secs,txt)
  type(Sector_type), dimension(:), intent(in) :: secs
  type(Sector_type) :: s
  character(len=*), intent(in) :: txt
  integer :: i
  do i = 1, size(secs)
    s = secs(i)
    if ( s%name=='NOTSET' ) then
      write(*,'(a,i3,a)') 'prSecs'//txt//':',i, 'End print_Sector_type'
      return
    else
      write(*,'(a,i3,3a10,3i3,2x,a)') 'prSecs'//txt//':',i, trim(s%name), &
        s%longname, s%cdfname, s%timefac, s%height, s%split,&
         trim(s%species)//'    :'//trim(s%description)
    end if
  end do
end subroutine print_Sector_type
!=========================================================================
endmodule OwnDataTypes_mod
