! <EmisDef_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2015 met.no
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

                         module EmisDef_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
implicit none

    !----------------- basic emissions file definitions --------------------!
    !  Here we define the parameters *not* likely to change often           !
    !  between different  model versions - e.g. the size and characteristics!
    !  of emission files and sector splits.                                 !
    !-----------------------------------------------------------------------!
    !  What is "Flat emissions":                                            !
    !  Most emission sources will have a seasonal, weekly, daily cycle. For !
    !  some sources there is no cycle, or the emission cycle is not known.  !
    !  For these source emissions will be constant (or flat) throughout the !
    !  year.                                                                !
    !-----------------------------------------------------------------------!
    ! Note that the names of the emission files are given in My_Emis_ml
    ! and often change from run to run.


    ! Note on SNAP sectors:
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



   integer, public, parameter :: NCMAX  =  11  ! Max. No. countries per grid point
   integer, public, parameter :: FNCMAX =  20  ! Max. No. countries (with
                                               ! flat emissions) per grid

  ! Sector specific information

   integer, public, parameter :: &
          NSECTORS  = 11       ! Number of SNAP-sectors in emissions

! Variables for NMR-NH3 project
! hb NH3emis (ISNAP_AGR, ISNAP_TRAF)
   integer, public, parameter :: &
          ANTROP_SECTORS=10, &   ! Non-natural sectors
          ISNAP_DOM  =  2,   &   ! Domestic/residential, for degree-day Timefactors
          ISNAP_NAT  = 11,   &   ! SNAP index for volcanoe emissions
          ISNAP_SHIP = 8,    &   ! SNAP index for flat emissions, e.g ship
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
   integer, public, parameter ::  NPOL  = 1 &   ! number of dust size modes
                                 ,QPOL  = 1 

   !Volcanos. 
   logical, public, parameter :: VOLCANOES_LL  = .true.  ! Read Volcanoes 
                                                         ! from VolcanoesLL.dat 
                                                         ! and disregard them 
                                                         ! from gridSOx

  integer, public, parameter :: IQ_DMS = 35  ! code for DMS emissions

!
real, public, save,  allocatable,dimension(:,:) ::  sumcdfemis ! Only used fby MasterProc
real, allocatable, public, save,  dimension(:,:) :: cdfemis
integer, allocatable, public, save,  dimension(:,:) :: nGridEmisCodes
integer, allocatable, public, save,  dimension(:,:,:):: GridEmisCodes
real, allocatable, public, save,  dimension(:,:,:,:,:):: GridEmis
! land-code information in each grid square - needed to know which country
! is emitting.                        
! nlandcode = No. countries in grid square
! landcode  = Country codes for that grid square
integer, public, save, allocatable, dimension(:,:)   :: nlandcode
integer, public, save, allocatable, dimension(:,:,:) :: landcode
! for flat emissions, i.e. no vertical extent:
integer, public, save, allocatable, dimension(:,:)   :: flat_nlandcode
integer, public, save, allocatable, dimension(:,:,:) :: flat_landcode
! for road dust emission potentials:
integer, public, save, allocatable, dimension(:,:)   :: road_nlandcode
integer, public, save, allocatable, dimension(:,:,:) :: road_landcode
! Emissions for input to chemistry routines
! KEMISTOP added to avoid hard-coded KMAX_MID-3:
!integer, public, parameter :: KEMISTOP = KMAX_MID - NEMISLAYERS + 1
real, public, allocatable, save, dimension(:,:,:,:) :: &
  gridrcemis,     & ! varies every time-step (as ps changes)
  gridrcemis0       ! varies every hour
real, public, allocatable, save, dimension(:,:,:) :: &
  gridrcroadd,    & ! Road dust emissions
  gridrcroadd0      ! varies every hour

!used for EEMEP 
real, allocatable, save, dimension(:,:,:,:)       ::  Emis_4D !(i,j,k,pollutant)
integer, save ::N_Emis_4D=0 !number of pollutants to read
integer, public, save :: Found_Emis_4D = 0 

integer, public, save :: KEMISTOP ! not defined yet= KMAX_MID - nemis_kprofile + 1


 ! Names of emis files, generated by GenChem:

     include 'CM_EmisFiles.inc'

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
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
                     end module EmisDef_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!_____________________________________________________________________________
