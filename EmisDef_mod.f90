! <EmisDef_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module EmisDef_mod

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
use ChemDims_mod,     only : NEMIS_File
use OwnDataTypes_mod, only : TXTLEN_NAME,TXTLEN_FILE, Emis_id_type, &
                             EmisFile_id_type, Emis_mask_type

implicit none
private

    !----------------- basic emissions file definitions --------------------!
    !  Here we define the parameters *not* likely to change often           !
    !  between different  model versions - e.g. the size and characteristics!
    !  of emission files and sector splits.                                 !
    !-----------------------------------------------------------------------!
    ! Note that the names of the emission files are given in My_Emis_mod
    ! and often change from run to run.


    ! SNAP sectors:
    ! ----------------------
    ! SNAP1  =  Combustion in energy and transformation industries,
    !           e.g.  public power stations, 150m nat gas
    ! SNAP2  =  Non-industrial combustion
    ! SNAP3  =  Industrial combustion  !60m nat gas
    ! SNAP4  =  Production processes
    ! SNAP5  =  Extraction and distribution of fossil fuels
    ! SNAP6  =  Solvent use
    ! SNAP7  =  Road traffic
    ! SNAP8  =  Other mobile sources (trains, planes, ships)
    ! SNAP9  =  Waste treatment and disposal
    ! SNAP10 =  Agriculture
    ! SNAP11 =  Nature



  ! First, define here the characteristics of the EMEP/CORINAIR
  ! SNAP sector emissions data-bases:



   integer, public, parameter :: NCMAX  =  14  ! Max. No. countries per grid point

  ! Sector specific information
   integer, save, public :: NSECTORS   ! Number of sectors used in emissions

   integer, save, public :: & !must be compatible with:
      ! timefac, fac_ehh24x7, fac_edd, fac_emm, fac_min, GridTfac, ISNAP_DOM, ISNAP_TRAF
          N_TFAC  = 11  ! Number of timefactor classes defined
   integer, save, pointer, dimension(:), public :: sec2tfac_map => null()! mapping of sector to time factor class
   integer, save, public :: & !must be compatible with:
          N_HFAC  = 11  ! Number of height distribution classes defined
         !DSHK N_HFAC  = 12  ! Number of height distribution classes defined
   integer, save, pointer, dimension(:), public :: sec2hfac_map => null()! mapping of sector to height distribution class
   integer, save, public :: & !must be compatible with: emisfrac
          N_SPLIT  = 11  ! Number of speciation classes defined
   integer, save, pointer, dimension(:), public :: sec2split_map => null()! mapping of sector to speciation class

!SNAP specific definitions
   integer, public, parameter :: &
          NSECTORS_SNAP  = 11    ! Number of sectors defined in SNAP emissions. Do not modify
   integer, save, target, dimension(NSECTORS_SNAP), public :: & ! mapping of sector to time factor class
        SNAP_sec2tfac_map = [1,2,3,4,5,6,7,8,9,10,11] !values must be <= N_TFAC
   integer, save, target, dimension(NSECTORS_SNAP), public :: & ! mapping of sector to height distribution class
        SNAP_sec2hfac_map = [1,2,3,4,5,6,7,8,9,10,11] !values must be <= N_HFAC
   integer, save, target, dimension(NSECTORS_SNAP), public :: & ! mapping of sector to height distribution class
        SNAP_sec2split_map = [1,2,3,4,5,6,7,8,9,10,11] !values must be <= N_SPECIATION
!   integer, save, dimension(NSECTORS_SNAP) ::snap2gnfr=(/1,3,2,4,13,5,6,7,10,11,-1/)
   integer, save, dimension(NSECTORS_SNAP,3), public ::snap2gnfr=reshape([1,3,2,4,13,5,6,7,10,11,-1 &
                                                          ,-1,-1,-1,-1,-1,-1,-1,8,-1,12,-1 &
                                                          ,-1,-1,-1,-1,-1,-1,-1,9,-1,-1,-1 &
                                                          ],shape(snap2gnfr))

!GNFR  specific definitions
   integer, public, parameter :: &
          NSECTORS_GNFR  = 13    ! Number of sectors defined in GNFR emissions
   integer, save, target, dimension(NSECTORS_GNFR), public :: & ! mapping of sector to time factor class
        GNFR_sec2tfac_map = [1,3,2,4,6,7,8,8,8,9,10,10,5] !values must be <= N_TFAC
   integer, save, target, dimension(NSECTORS_GNFR), public :: & ! mapping of sector to height distribution class
        GNFR_sec2hfac_map = [1,3,2,4,6,7,8,8,8,9,10,10,5] !values must be <= N_HFAC
   integer, save, target, dimension(NSECTORS_GNFR), public :: & ! mapping of sector to height distribution class
        GNFR_sec2split_map = [1,3,2,4,6,7,8,8,8,9,10,10,5] !values must be <= N_SPECIATION

   integer, save, dimension(NSECTORS_GNFR), public ::gnfr2snap=(/1,3,2,4,6,7,8,-1,-1,9,10,-1,5/)

!TEST to try own defintions
   integer, public, parameter:: &
          NSECTORS_TEST  = 11 ! Number of sectors defined. Must match the sizes of the maps below
   integer, save, target, dimension(NSECTORS_TEST), public :: &
      ! mapping of sector to time factor class. values must be <= N_TFAC
        TEST_sec2tfac_map  = [1,2,3,4,5,6,7,8,9,10,11] &
      ! mapping of sector to height distribution class. vals must be <= N_HFAC
       ,TEST_sec2hfac_map  = [1,2,3,4,5,6,7,8,9,10,11] & 
      ! mapping of sector to height distribution class ! must be<=N_SPECIATION
       ,TEST_sec2split_map = [1,2,3,4,5,6,7,8,9,10,11] 


!The sectors defined here are always SNAP sectors. Should NOT be changed if other
!categories (for instance GNFR) are used!
   integer, public, parameter :: &
          ANTROP_SECTORS=10, &   ! Non-natural sectors
          ISNAP_DOM  =  2,   &   ! Domestic/residential, for degree-day Timefactors
          ISNAP_AGR  = 10,   &   ! Note that flat emissions do NOT necessarily
          ISNAP_TRAF = 7         ! belong to the same SNAP sector

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
integer,  public, save :: NEmis_id
integer,  public :: NEmisMask = 0 !number of masks defined (new format)
real,  public, allocatable :: EmisMaskValues(:,:,:) ! size will be (LIMAX,LJMAX,NEmisMask)
type(Emis_id_type), public, save:: Emis_id(10)
type(EmisFile_id_type), public, save:: EmisFiles(10) !list of emission files after validation
integer,  public, parameter:: NEmis_sourcesMAX = 50
type(Emis_id_type), public, save:: Emis_source(NEmis_sourcesMAX) !list of valid sources found in the emission files
integer,  public, save :: NEmis_sources = 0
integer,  public, save :: NEmis_3Dsources = 0
integer,  public, save :: NEmisFile_sources = 0
integer,  public, save :: ix3Dmap(1000) = 0
real, allocatable, public, save,  dimension(:,:,:):: Emis_source_2D !One 2D map for each source
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

!should be defined somewhere else?
real, public, allocatable, dimension(:,:,:,:,:,:), save :: &
  loc_frac&    ! Fraction of pollutants that are produced locally
  ,loc_frac_hour_inst&  !Houry local fractions
  ,loc_frac_hour&  !Houry average of local fractions
  ,loc_frac_day&  !Daily average of local fractions
  ,loc_frac_month&  !Monthly average of local fractions
  ,loc_frac_full  !Fullrun average of local fractions
real, public, allocatable, dimension(:,:,:,:), save :: &
   loc_tot_hour_inst&   !all contributions
  ,loc_tot_hour&   !Hourly average of all contributions
  ,loc_tot_day&   !Daily average of all contributions
  ,loc_tot_month&  !Monthly average of all contributions
  ,loc_tot_full  !Fullrun average of all contributions
real, public, allocatable, dimension(:,:,:,:), save :: &
  loc_frac_1d  ! Fraction of pollutants without i or j and extended (0:limax+1 or 0:ljmax+1)
integer, public, parameter:: Nneighbors = 9 !localfractions from 8 neighbors + self

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
