! <Emissions_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_10(3282)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2016 met.no
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
module Emissions_ml
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Calls up emission read/set routines
!  This routine interfaces the stand-alone emission-file reading routines
!  with the 3D model.
!_____________________________________________________________________________

use Biogenics_ml,     only: SoilNOx, AnnualNdep
use CheckStop_ml,     only: CheckStop,StopAll
use ChemSpecs,        only: NSPEC_SHL, NSPEC_TOT,NO2, SO2,species,species_adv
use Chemfields_ml,    only: xn_adv
use Country_ml,       only: MAXNLAND,NLAND,Country,IC_NAT,IC_FI,IC_NO,IC_SE
use Country_ml,       only: EU28,EUMACC2 !CdfSnap
use EmisDef_ml,       only: &
      NEMIS_FILE    & ! No. emission files
     ,EMIS_FILE     & ! Names of species ("sox  ",...)
     ,NCMAX         & ! Max. No. countries per grid
     ,FNCMAX        & ! Max. No. countries (with flat emissions) per grid
     ,ISNAP_DOM     & ! snap index for domestic/resid emis
     ,ISNAP_TRAF    & ! snap index for road-traffic (SNAP7)
     ,ISEC_NAT      & ! index for natural (and flat?) emissions
     ,ISEC_SHIP     & ! index for flat emissions, e.g ship
     ,IQ_DMS        & ! code for DMS emissions
     ,NROAD_FILES   & ! No. road dust emis potential files
     ,ROAD_FILE     & ! Names of road dust emission files
     ,NROADDUST     & ! No. road dust components 
     ,QROADDUST_FI  & ! fine road dust emissions (PM2.5) 
     ,QROADDUST_CO  & ! coarse road dust emis
     ,ROADDUST_FINE_FRAC  & ! fine (PM2.5) fraction of road dust emis
     ,ROADDUST_CLIMATE_FILE &! TEMPORARY! file for road dust climate factors 
     ,nGridEmisCodes,GridEmisCodes,GridEmis,cdfemis&
     ,nlandcode,landcode,flat_nlandcode,flat_landcode&
     ,road_nlandcode,road_landcode&
     ,gridrcemis,gridrcemis0,gridrcroadd,gridrcroadd0&
     ,DMS_natso2_month, DMS_natso2_year,O_NH3, O_DMS&
     ,Emis_4D,N_Emis_4D,Found_Emis_4D & !used for EEMEP 
     ,KEMISTOP&
     ,MAXFEMISLONLAT,N_femis_lonlat,loc_frac &
     ,NSECTORS, N_HFAC, N_TFAC, N_SPLIT     & ! No. emis sectors, height, time and split classes
     ,sec2tfac_map, sec2hfac_map, sec2split_map& !generic mapping of indices
     ,NSECTORS_SNAP, SNAP_sec2tfac_map, SNAP_sec2hfac_map, SNAP_sec2split_map&!SNAP specific mapping
     ,NSECTORS_GNFR, GNFR_sec2tfac_map, GNFR_sec2hfac_map, GNFR_sec2split_map&!GNFR specific mapping
     ,NSECTORS_TEST, TEST_sec2tfac_map, TEST_sec2hfac_map, TEST_sec2split_map!TEST specific mapping

use EmisGet_ml,       only: &
     EmisSplit &
    ,EmisGetCdf &  ! 
    ,EmisGetCdfFrac &
    ,EmisGetASCII &
    ,femis                       &  ! Gets scaling factors -> e_fact
    ,e_fact                      &  ! scaling factors
    ,e_fact_lonlat               &  ! scaling factors
    ,EmisHeights                 &  ! Generates vertical distrib
    ,nrcemis, nrcsplit, emisfrac &  ! speciation routines and array
    ,nemis_kprofile, emis_kprofile &! Vertical emissions profile
    ,iqrc2itot                   &  ! maps from split index to total index
    ,emis_masscorr               &  ! 1/molwt for most species
    ,emis_nsplit                 &  ! No. species per emis file
    ,RoadDustGet                 &  
    ,roaddust_masscorr           &   ! 1/200. Hard coded at the moment, needs proper setting in EmisGet_ml...   
    ,femis_latmin,femis_latmax,femis_lonmin,femis_lonmax
use GridValues_ml,    only: GRIDWIDTH_M    & ! size of grid (m)
                           ,xm2            & ! map factor squared
                           ,debug_proc,debug_li,debug_lj & 
                           ,xmd,dA,dB,i_fdom,j_fdom,glon,glon,glat
use Io_Nums_ml,       only: IO_LOG, IO_DMS, IO_EMIS, IO_TMP
use Io_Progs_ml,      only: ios, open_file, datewrite, PrintLog
use MetFields_ml,     only: u_xmj, v_xmi, roa, ps, z_bnd, surface_precip,EtaKz ! ps in Pa, roa in kg/m3
use MetFields_ml,     only: t2_nwp   ! DS_TEST SOILNO - was zero!
use ModelConstants_ml,only: &
    KMAX_MID, KMAX_BND, PT ,dt_advec, &
    EMIS_SOURCE,   &    ! emislist, CdfFractions
    EMIS_TEST,     &    ! CdfSnap or none
    emis_inputlist, &   !TESTC
    EmisDir,      &    ! template for emission path
    DataDir,      &    ! template for path
    EMIS_OUT,      &    ! output emissions in ASCII or not
!    MONTHLY_GRIDEMIS, &  !NML
    NBVOC,         &    ! > 0 if forest voc wanted
    INERIS_SNAP2 , &    ! INERIS/TFMM HDD20 method
    DEBUG, MYDEBUG => DEBUG_EMISSIONS,  MasterProc, & 
    DEBUG_SOILNOX, DEBUG_EMISTIMEFACS, DEBUG_ROADDUST, &
    USES,  &  ! Gives USES%EMISSTACKS, DEGREEDAY_FACTORS,GRIDDED_EMIS_MONTHLY_FACTOR
    SEAFIX_GEA_NEEDED, & ! see below
    USE_LIGHTNING_EMIS,USE_AIRCRAFT_EMIS,USE_ROADDUST, &
    USE_EURO_SOILNOX, USE_GLOBAL_SOILNOX, EURO_SOILNOX_DEPSCALE,&! one or the other
    USE_OCEAN_NH3,USE_OCEAN_DMS,FOUND_OCEAN_DMS,&
    NPROC, EmisSplit_OUT,USE_uEMEP,uEMEP,SECTORS_NAME
use MPI_Groups_ml  , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_REAL8, MPI_INTEGER&
                            ,MPI_SUM,MPI_COMM_CALC, IERROR
use NetCDF_ml,        only: ReadField_CDF,ReadField_CDF_FL,ReadTimeCDF,IsCDFfractionFormat,&
                            GetCDF_modelgrid,PrintCDF,ReadSectorName
use netcdf
use Par_ml,           only: MAXLIMAX,MAXLJMAX, GIMAX,GJMAX, IRUNBEG,JRUNBEG,&
                            me,limax,ljmax, MSG_READ1,MSG_READ7&
                           ,gi0,gj0,li0,li1,lj0,lj1
use PhysicalConstants_ml,only: GRAV, AVOG, ATWAIR
use PointSource_ml,      only: readstacks !MKPS
use Setup_1dfields_ml,only: rcemis   ! ESX
use SmallUtils_ml,    only: find_index,  key2str
use ReadField_ml,     only: ReadField    ! Reads ascii fields
use TimeDate_ml,      only: nydays, nmdays, date, current_date, &! No. days per 
                            tdif_secs,timestamp,make_timestamp,daynumber,day_of_week ! year, date-type, weekday 
use TimeDate_ExtraUtil_ml,only :nctime2date
use Timefactors_ml,   only: &
     NewDayFactors          & ! subroutines
    ,DegreeDayFactors       & ! degree-days used for SNAP-2
    ,Gridded_SNAP2_Factors, gridfac_HDD & 
    ,fac_min,timefactors   &                  ! subroutine
    ,fac_ehh24x7 ,fac_emm, fac_edd, timefac & ! time-factors
    ,Read_monthly_emis_grid_fac &
    ,GridTfac !array with monthly gridded time factors

implicit none
private

! subroutines:
public :: Emissions         ! Main emissions module 
public :: newmonth
public :: EmisSet           ! Sets emission rates every hour/time-step
public :: EmisOut           ! Outputs emissions in ascii
public :: uemep_emis

! The main code does not need to know about the following 
private :: expandcclist            !  expands e.g. EU28, EUMACC2
private :: consistency_check       ! Safety-checks

!
! The output emission matrix for the 11-SNAP data is snapemis:
!
real, private, allocatable, dimension(:,:,:,:,:), save :: &
  snapemis      ! main emission arrays, in kg/m2/s

real, private, allocatable, dimension(:,:,:,:), save :: &
  snapemis_flat ! main emission arrays, in kg/m2/s  

real, private, allocatable, dimension(:,:,:,:), save :: &
! Not sure if it is really necessary to keep the country info; gives rather messy code but consistent with the rest at least (and can do the seasonal scaling for Nordic countries in the code instead of as preprocessing) 
  roaddust_emis_pot ! main road dust emission potential arrays, in kg/m2/s (to be scaled!)

! We store the emissions for output to d_2d files and netcdf in kg/m2/s
real, public, allocatable, dimension(:,:,:), save :: SumSnapEmis,SumSplitEmis

logical, save, private  :: first_dms_read


! and for budgets (not yet used - not changed dimension)
real, public,  save, dimension(NSPEC_SHL+1:NSPEC_TOT) ::  totemadd
integer, private, save :: iemCO  ! index of CO emissions, for debug

logical :: Cexist,USE_MONTHLY_GRIDEMIS=.false.!internal flag
real ::TimesInDays(120),mpi_out
integer ::NTime_Read=-1,ncFileID,VarID,found, cdfstatus
character(len=125) ::fileName_monthly='NOT_SET'!must be initialized with 'NOT_SET'
character(len=10), private,save ::  incl_monthly(size(emis_inputlist(1)%incl)),&
     excl_monthly(size(emis_inputlist(1)%excl))
integer, private,save :: nin_monthly, nex_monthly, index_monthly

contains
!***********************************************************************
subroutine Emissions(year)
! Calls main emission reading routines
!----------------------------------------------------------------------!
! DESCRIPTION:
!   0) Call set_molwts and set_emisconv_and_iq, followed by
!      consistency check
!   1) Call some setups:
!         Country_Init
!         timefactors: monthly and daily factors, + time zone
!                            -> fac_emm, fac_edd arrays, timezone
!   2) Read in emission correction file femis
!   3) Call emis_split for speciations
!   4) Read the annual emission totals in each grid-square, converts
!      units and corrects using femis data. 
!
!   The output emission matrix for the 11-SNAP data is snapemis:
!
!   real    snapemis (NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILES)
!  
!----------------------------------------------------------------------!
!--arguments
  integer, intent(in)   :: year        ! Year ( 4-digit)

  !-- local variables
  integer :: i, j              ! Loop variables
  real    :: tonne_to_kgm2s    ! Converts tonnes/grid to kg/m2/s
  real    :: ccsum             ! Sum of emissions for one country

  ! arrays for whole EMEP area:
  ! additional arrays on host only for landcode, nlandcode
  ! BIG arrays ... will be needed only on me=0. Make allocatable
  ! to reduce static memory requirements.
  real,    allocatable, dimension(:,:,:)    :: globroad_dust_pot ! Road dust emission potentials
  integer, allocatable, dimension(:,:)      :: road_globnland 
  integer, allocatable, dimension(:,:,:)    :: road_globland 
  real,    allocatable, dimension(:,:)      :: RoadDustEmis_climate_factor ! Climatic factor for scaling road dust emissions (in TNO model based on yearly average soil water)
  integer :: err1, err2, err3, err4, err7, err8, err9 ! Error messages
  integer :: fic ,insec,inland,iemis,iemislist 
  integer :: iic,ic,n         ! country codes 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over pollutants (1..NEMIS_FILE)
  integer :: icc              ! loop variables over  sources
  integer :: nin, nex         !  include, exclude numbers for emis_inputlist
  character(len=*), parameter :: sub='Emissions:'
  character(len=300) :: fname ! txt, File name

  ! Emission sums (after e_fact adjustments):
  real, dimension(NEMIS_FILE)       :: emsum ! Sum emis over all countries
  real, dimension(NLAND,NEMIS_FILE) :: sumemis, sumemis_local ! Sum of emissions per country
  real, dimension(NEMIS_FILE) :: sumEU ! Sum of emissions in EU


  ! Road dust emission potential sums (just for testing the code, the actual emissions are weather dependent!)
  real, dimension(NLAND,NROAD_FILES) :: sumroaddust    ! Sum of emission potentials per country
  real, dimension(NLAND,NROAD_FILES) :: sumroaddust_local    ! Sum of emission potentials per country in subdomain
  real :: fractions(LIMAX,LJMAX,NCMAX),SMI(LIMAX,LJMAX),SMI_roadfactor
  logical ::SMI_defined=.false.
  logical,save :: my_first_call=.true.  ! Used for femis call
  logical :: fileExists            ! to test emission files
  character(len=40) :: varname, fmt,cdf_sector_name
  integer ::allocerr, i_Emis_4D

  if (MasterProc) write(6,*) "Reading emissions for year",  year



  ! Get emission input definitions and other initializations
  !>==============================
  if(my_first_call) then
    sumemis(:,:) =  0.0       ! initialize sums
    ios = 0

    if(EMIS_TEST=="CdfSnap" .or. EMIS_SOURCE=="Mixed") then  

      if(MasterProc)  write(*,*)sub//" Mixed format"
      do iemislist = 1, size( emis_inputlist(:)%name )
        fname = emis_inputlist(iemislist)%name
        if(fname=="NOTSET") cycle
        if(MasterProc)&
          write(*,*)"Emission source number ", iemislist,"from ",sub//trim(fname)

        if(emis_inputlist(iemislist)%type == "sectors".or.&
           emis_inputlist(iemislist)%type == "GNFRsectors")then ! Expand groups, e.g. EUMACC2

          call expandcclist( emis_inputlist(iemislist)%incl , n)
          emis_inputlist(iemislist)%Nincl = n
          if(MasterProc) write(*,*) sub//trim(fname)//" INPUTLIST-INCL", n

          call expandcclist( emis_inputlist(iemislist)%excl , n)
          emis_inputlist(iemislist)%Nexcl = n
          if(MasterProc) write(*,*) sub//trim(fname)//" INPUTLIST-EXCL", n

        endif

        !replace keywords
22      format(5A)

        if(MasterProc)write(*,22)'original emission name ',trim(fname)
        fname = key2str(fname,'EmisDir',EmisDir)
        fname = key2str(fname,'DataDir',DataDir)
        fname = key2str(fname,'YYYY',year)
        emis_inputlist(iemislist)%name=trim(fname)
        if(MasterProc)write(*,22)'filename redefined as: ',&
          trim(emis_inputlist(iemislist)%name)

        nin_monthly = 0
        nex_monthly = 0
        cdf_sector_name='NOTSET'
        call ReadSectorname(fname,cdf_sector_name)
        if(trim(cdf_sector_name)/='NOTSET')then
           SECTORS_NAME=trim(cdf_sector_name)
           if(Masterproc)write(*,*)"Switching sector categories to ",trim(SECTORS_NAME)
           if(Masterproc)write(IO_LOG,*)"Switching sector categories to ",trim(SECTORS_NAME)
        endif
      enddo ! iemislist
    endif

    !>============================
    

    ! 0) set molwts, conversion factors (e.g. tonne NO2 -> tonne N), and
    !    emission indices (IQSO2=.., )
    
    ! init_sectors
    if(SECTORS_NAME=='SNAP')then
       !11 sectors defined in emissions
       NSECTORS = NSECTORS_SNAP
       !map timefactors onto SNAP map
       sec2tfac_map => SNAP_sec2tfac_map
       sec2hfac_map => SNAP_sec2hfac_map
       sec2split_map => SNAP_sec2split_map
    else  if(SECTORS_NAME=='GNFR')then
       !13 sectors defined in emissions
       NSECTORS = NSECTORS_GNFR
       !map timefactors onto GNFR map
       sec2tfac_map => GNFR_sec2tfac_map
       sec2hfac_map => GNFR_sec2hfac_map
       sec2split_map => GNFR_sec2split_map
    else  if(SECTORS_NAME=='TEST')then
       !11 sectors defined in emissions
       NSECTORS = NSECTORS_TEST
       !map timefactors onto TEST map
       sec2tfac_map => TEST_sec2tfac_map
       sec2hfac_map => TEST_sec2hfac_map
       sec2split_map => TEST_sec2split_map
    else
       call StopAll("Sectors not defined")
    endif
    
    allocate(cdfemis(LIMAX,LJMAX))
    allocate(nGridEmisCodes(LIMAX,LJMAX))
    allocate(GridEmisCodes(LIMAX,LJMAX,NCMAX))
    allocate(GridEmis(NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILE),stat=allocerr)
    call CheckStop(allocerr /= 0, &
         "EmisGet:Allocation error for GridEmis")
    GridEmisCodes = -1   
    nGridEmisCodes = 0
    GridEmis = 0.0
    
    allocate(nlandcode(LIMAX,LJMAX),landcode(LIMAX,LJMAX,NCMAX))
    nlandcode=0
    landcode=0
    allocate(flat_nlandcode(LIMAX,LJMAX),flat_landcode(LIMAX,LJMAX,FNCMAX))
    flat_nlandcode=0
    flat_landcode=0
    allocate(road_nlandcode(LIMAX,LJMAX),road_landcode(LIMAX,LJMAX,NCMAX))
    road_nlandcode=0
    road_landcode=0
    allocate(snapemis(NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILE))
    snapemis=0.0
    allocate(snapemis_flat(LIMAX,LJMAX,FNCMAX,NEMIS_FILE))
    snapemis_flat=0.0
    allocate(roaddust_emis_pot(LIMAX,LJMAX,NCMAX,NROAD_FILES))
    roaddust_emis_pot=0.0
    allocate(SumSnapEmis(LIMAX,LJMAX,NEMIS_FILE))
    SumSnapEmis=0.0
    if(USE_uEMEP)then
       allocate(loc_frac(LIMAX,LJMAX,KMAX_MID,1))
       loc_frac=0.0
    endif
    !=========================
    !  call Country_Init() ! In Country_ml, => NLAND, country codes and names, timezone
    
    allocate(e_fact(NSECTORS,NLAND,NEMIS_FILE))!NLAND defined in Country_Init()
    e_fact=1.0
    allocate(e_fact_lonlat(NSECTORS,MAXFEMISLONLAT,NEMIS_FILE))
    e_fact_lonlat=1.0
    if(.not.allocated(timefac))allocate(timefac(NLAND,N_TFAC,NEMIS_FILE))
    if(.not.allocated(fac_ehh24x7))allocate(fac_ehh24x7(N_TFAC,24,7))
    if(.not.allocated(fac_emm))allocate(fac_emm(NLAND,12,N_TFAC,NEMIS_FILE))
    if(.not.allocated(fac_min))allocate(fac_min(NLAND,N_TFAC,NEMIS_FILE))
    if(.not.allocated(fac_edd))allocate(fac_edd(NLAND, 7,N_TFAC,NEMIS_FILE))
        
    call femis()              ! emission factors (femis.dat file)
    if(ios/=0) return
    my_first_call = .false.
  endif


  ! The GEA emission data, which is used for EUCAARI runs on the HIRHAM
  ! domain have in several sea grid cells non-zero emissions in other sectors
  ! than SNAP8 and there are also NH3 emission over sea areas. The former 
  ! problem makes the code crash if the sea areas are defined  as 
  ! sea (sea=T), so we treat them as land in the EUCAARI/HIRHAM runs 
  ! (sea=F). This is a problem with GEA emission data only, not the 
  ! HIRHAM domain! When e.g. interpolated EMEP emissions are used on
  ! the HIRHAM domain, this is not a problem.
  if(SEAFIX_GEA_NEEDED) then ! Special fix for HIRHAM/GEA
    ! Could have hard-coded 30-34, but best to avoid:
    do n=1,5
      select case(n)
      case(1);ic=find_index("BAS",Country(:)%code)!could also use directly ic=IC_BAS
      case(2);ic=find_index("NOS",Country(:)%code)
      case(3);ic=find_index("ATL",Country(:)%code)
      case(4);ic=find_index("MED",Country(:)%code)
      case(5);ic=find_index("BLS",Country(:)%code)
      endselect
      call CheckStop(ic<1,"Country_Init error in HIRHAM/GEA fix")
      Country(ic)%is_sea = .false.
    enddo
  endif ! HIRHAM/GEA fix

  !=========================
  call consistency_check()               ! Below
  !=========================
  ios = 0

  if(USES%DEGREEDAY_FACTORS) call DegreeDayFactors(0)! See if we have gridded SNAP-2

  call EmisHeights()     ! vertical emissions profile
  KEMISTOP = KMAX_MID - nemis_kprofile + 1

  if( USES%EMISSTACKS ) call readstacks(IO_EMIS)

  if(MasterProc) then   !::::::: ALL READ-INS DONE IN HOST PROCESSOR ::::
    write(*,*) "Reading monthly and daily timefactors"
    if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)then
      write(*,*)"Emissions using gridded monhtly timefactors "
      write(IO_LOG,*)"Emissions using gridded monhtly timefactors "       
    endif
    !=========================
    call timefactors(year)               ! => fac_emm, fac_edd
    !=========================
  endif
  !=========================
  call EmisSplit()    ! In EmisGet_ml, => emisfrac
  !=========================
  !Must first call EmisSplit, to get nrcemis defined
  if(EmisSplit_OUT)then
    allocate(SumSplitEmis(LIMAX,LJMAX,nrcemis))
    SumSplitEmis=0.0
  endif
  !=========================
  call CheckStop(ios, "ioserror: EmisSplit")

  ! ####################################
  ! Broadcast  monthly and Daily factors (and hourly factors if needed/wanted)
  CALL MPI_BCAST(fac_emm,8*NLAND*12*N_TFAC*NEMIS_FILE,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 
  CALL MPI_BCAST(fac_edd,8*NLAND*7*N_TFAC*NEMIS_FILE,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 
  CALL MPI_BCAST(fac_ehh24x7,8*N_TFAC*24*7,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 

  !define fac_min for all processors
  forall(iemis=1:NEMIS_FILE,insec=1:N_TFAC,inland=1:NLAND) &
       fac_min(inland,insec,iemis) = minval(fac_emm(inland,:,insec,iemis))
  if(INERIS_SNAP2) & !  INERIS do not use any base-line for SNAP2
       fac_min(:,ISNAP_DOM,:) = 0.

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ! c4b) Set up DMS factors here - to be used in newmonth
  !      Taken from IQ_DMS=35 for SO2 nature (sector 11)
  !      first_dms_read is true until first call to newmonth finished.
  first_dms_read = .true. 

  ! 4) Read emission files 

  ! allocate for MasterProc (me:=0) only:
  err1 = 0
  if(MasterProc) then
    if(USE_ROADDUST)then
      allocate(road_globnland(GIMAX,GJMAX),stat=err7)
      allocate(road_globland(GIMAX,GJMAX,NCMAX),stat=err8)
      allocate(globroad_dust_pot(GIMAX,GJMAX,NCMAX),stat=err9)
      allocate(RoadDustEmis_climate_factor(GIMAX,GJMAX),stat=err1)

      call CheckStop(err7, "Allocation error 7 - globroadland")
      call CheckStop(err8, "Allocation error 8 - globroadland")
      call CheckStop(err9, "Allocation error 9 - globroad_dust_pot")
      call CheckStop(err1, "Allocation error 1 - RoadDustEmis_climate_factor")
    endif ! road dust

    ! Initialise with 0
    sumemis_local(:,:)=0.0
    emsum=0.0

    if(USE_ROADDUST)then
      road_globnland(:,:)=0
      road_globland(:,:,:)=0
      globroad_dust_pot(:,:,:)=0.
      RoadDustEmis_climate_factor(:,:)=1.0 ! default, no scaling
    endif ! road dust
  else
    ! needed for DEBUG=yes compilation options
    if(USE_ROADDUST)then
      allocate(road_globnland(1,1),road_globland(1,1,1),&
           globroad_dust_pot(1,1,1),stat=err9)
      call CheckStop(err9, "Allocation error 9 - dummy roadglob")
    endif ! road dust
  endif


  select case(EMIS_SOURCE)
  case("Mixed") !default, and only one in use

    Found_Emis_4D=0
    do iemislist = 1, size( emis_inputlist(:)%name )

      fname=emis_inputlist(iemislist)%name
      if ( fname == "NOTSET" ) cycle
      38 FORMAT(A,I4,A)
      if(MasterProc)write(*,38)sub//' reading emis_inputlist ',iemislist,trim(fname)

      sumemis=0.0
      sumemis_local(:,:)=0.0

      nin = emis_inputlist(iemislist)%Nincl
      nex = emis_inputlist(iemislist)%Nexcl

      !1) emissions in NetCDF Fractions format 
      if(IsCDFfractionFormat(trim(fname)))then
        !the file is in "fraction" format
        !check first if the file is a monthly file
        NTime_Read=-1
        call ReadTimeCDF(trim(fname),TimesInDays,NTime_Read)

        if(NTime_Read==12)then
          
          USE_MONTHLY_GRIDEMIS=.true.
          emis_inputlist(iemislist)%periodicity = "monthly"
          if(MasterProc)write(*,*)'Assuming monthly emissions in CdfFractions format'
          call CheckStop(trim(fileName_monthly)/='NOT_SET', &
               "Can use only one monthly emissions. Already defined "//trim(fileName_monthly))
          fileName_monthly=trim(fname)!will be read later
          nin_monthly=nin
          nex_monthly=nex
          incl_monthly=emis_inputlist(iemislist)%incl
          excl_monthly=emis_inputlist(iemislist)%excl
          index_monthly=iemislist

          if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)then
            write(*,*)"Uncompatible settings: you use monthly emissions and GRIDDED_EMIS_MONTHLY_FACTOR=T "
            !If you really want this, you can uncomment the stop
            call StopAll("monthly emissions and GRIDDED_EMIS_MONTHLY_FACTOR=T not allowed ")
          endif
        else
          !yearly grid independent netcdf fraction format emissions                                
          do iem = 1, NEMIS_FILE
             if(emis_inputlist(iemislist)%pollName(1)/='NOTSET')then
                if(emis_inputlist(iemislist)%pollName(1)/=trim(EMIS_FILE(iem)))cycle      
             endif
             do isec=1,NSECTORS                      
              write(varname,"(A,I2.2)")trim(EMIS_FILE(iem))//'_sec',isec
              call EmisGetCdfFrac(iem, isec, fname, varname, sumemis_local, &
                   emis_inputlist(iemislist)%incl, nin, emis_inputlist(iemislist)%excl, nex)

            enddo!sectors
          enddo!NEMIS_FILE

          !add together totals from each processor (only me=0 get results)
          sumemis=0.0
          CALL MPI_REDUCE(sumemis_local,sumemis,&
               NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)        

        endif

      elseif(index(emis_inputlist(iemislist)%name,"Emis_4D.nc")>0)then 
        !under development
        Found_Emis_4D=iemislist
        N_Emis_4D = 0
        do i_Emis_4D=1,20
          if(emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)=='NOTSET')exit
          N_Emis_4D = N_Emis_4D +1
          if(MasterProc)write(*,*)'Emis_4D: will read ',emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)
          n=find_index(emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D),species(:)%name)
          if(MasterProc)then
            if(n>0)then
              write(*,*)'Emis_4D: will write to ',n,emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D)        
            else                   
              write(*,*)'Emis_4D: WARNING did not find ',emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D),' among the emep species'
              write(*,*)'Emis_4D: WARNING ',emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D),' is not used'
            endif
          endif
        enddo
        !   else if(IsCDFSnapFormat(trim(emis_inputlist(iemislist)%name)))then !This Does not work because of "POLL"
      elseif(emis_inputlist(iemislist)%type == "sectors" .and. index(emis_inputlist(iemislist)%name,".nc")>1)then 
        !not in "fraction" format. Each land has own set of fields
        !Each pollutant has own file. 
        if(MasterProc)  write(*,*)sub//trim(fname)//" Processing"
        do iem = 1, NEMIS_FILE

           if(emis_inputlist(iemislist)%pollName(1)/='NOTSET')then
              if(emis_inputlist(iemislist)%pollName(1)/=trim(EMIS_FILE(iem)))cycle      
           endif

           fname = key2str(emis_inputlist(iemislist)%name,'POLL',EMIS_FILE(iem))

          if(MasterProc) then
            write(*,*)sub//trim(fname)//" iemProcess",iem,trim(fname)
            write(*,"(a,2i3,a,3i3)") "INPUTLIST:", iem, iemislist, trim(fname), nin, nex,me
            inquire(file=fname,exist=fileExists)
            if(.not.fileExists) write(*,"(a)") 'WARNING EMISFile missing! '//trim(fname)
            write(*,*)sub//trim(fname)//" REPLACE ",iem,trim(fname),&
              key2str(emis_inputlist(iemislist)%name,'POLL',EMIS_FILE(iem))
          endif

          call CheckStop( nin>0 .and. nex > 0, &
              "emis_inputlists cannot have inc and exc")
          if ( nin > 0 ) then
            call EmisGetCdf(iem,fname, sumemis(1,iem), &
                 incl=emis_inputlist(iemislist)%incl(1:nin) )
          elseif (  nex > 0 ) then
            call EmisGetCdf(iem,fname, sumemis(1,iem), &
                 excl=emis_inputlist(iemislist)%excl(1:nex) ) 
          else
            call EmisGetCdf(iem,fname, sumemis(1,iem))
          endif
          if(MasterProc) write(*,*) "PARTEMIS ", iem, trim(fname), sumemis(27,iem) 

        enddo
      elseif(index(emis_inputlist(iemislist)%name,"grid")>0)then
        !ASCII format
        do iem = 1, NEMIS_FILE
          fname=key2str(emis_inputlist(iemislist)%name,'POLL',EMIS_FILE(iem)) ! e.g. POLL -> sox
          if(MasterProc)write(*,fmt='(A)')'Reading ASCII format '//trim(fname)
          call EmisGetASCII(iem, fname, trim(EMIS_FILE(iem)), sumemis_local, &
               emis_inputlist(iemislist)%incl, nin, emis_inputlist(iemislist)%excl, nex)
        enddo

        !add together totals from each processor (only me=0 get results)
        sumemis=0.0
        CALL MPI_REDUCE(sumemis_local,sumemis,&
             NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)        

      elseif(emis_inputlist(iemislist)%type == "OceanNH3")then
        if(MasterProc)write(*,*)' using  OceanNH3'    
        O_NH3%index=find_index("NH3",species(:)%name)
        call CheckStop(O_NH3%index<0,'Index for NH3 not found')
        NTime_Read=-1
        call ReadTimeCDF(trim(fname),TimesInDays,NTime_Read)
        if(NTime_Read==12)then
          emis_inputlist(iemislist)%periodicity = "monthly"
          allocate(O_NH3%emis(LIMAX,LJMAX))
          allocate(O_NH3%map(LIMAX,LJMAX))
          O_NH3%emis=0.0
          O_NH3%map=0.0
          O_NH3%sum_month=0.0
          O_NH3%sum_year=0.0
          if(MasterProc)write(*,*)' found  OceanNH3 monthly'    
        else
          call StopAll("Yearly OceanNH3 not implemented")
        endif
      elseif (emis_inputlist(iemislist)%type == "DMS")then
        if(MasterProc)write(*,*)'using DMS'    
        USE_OCEAN_DMS=.true.
        O_DMS%index=find_index("SO2",species(:)%name)
        call CheckStop(O_DMS%index<0,'Index for SO2 not found')
        NTime_Read=-1
        call ReadTimeCDF(trim(fname),TimesInDays,NTime_Read)
        if(NTime_Read==12)then
          emis_inputlist(iemislist)%periodicity = "monthly"
          allocate(O_DMS%emis(LIMAX,LJMAX))
          allocate(O_DMS%map(LIMAX,LJMAX))
          O_DMS%emis=0.0
          O_DMS%map=0.0
          O_DMS%sum_month=0.0
          O_DMS%sum_year=0.0
          if(MasterProc)write(*,*)' found DMS monthly'    
        else
          call StopAll("Yearly DMS not implemented")
        endif
      else
        if(MasterProc)write(*,*)'WARNING: did not recognize format of '//trim(emis_inputlist(iemislist)%name)
        call StopAll("Emissions file format not recognized ")
      endif


      if(MasterProc.and. emis_inputlist(iemislist)%periodicity == "once") then
        call PrintLog("Total emissions by countries for "//trim(emis_inputlist(iemislist)%name)//" (Gg)")
        write(*     ,"(2a4,3x,30(a12,:))")"  N "," CC ",EMIS_FILE(:)
        write(IO_LOG,"(2a4,3x,30(a12,:))")"  N "," CC ",EMIS_FILE(:)                
        sumEU(:) = 0.0
        fmt="(i4,1x,a4,3x,30(f12.2,:))"
        do ic = 1, NLAND
          ccsum = sum( sumemis(ic,:) )
          icc=Country(ic)%icode
          if ( ccsum > 0.0 )then
            write(*,     fmt) icc, Country(ic)%code, sumemis(ic,:)
            write(IO_LOG,fmt) icc, Country(ic)%code, sumemis(ic,:)
          endif
          if(find_index(Country(ic)%code,EU28(:))>0) sumEU = sumEU + sumemis(ic,:)
        enddo
        if ( sum(sumEU(:))>0.001) then
          write(*     ,fmt) 0, "EU", sumEU(:)
          write(IO_LOG,fmt) 0, "EU", sumEU(:)
        endif
      endif

      !total of emissions from all countries and files into emsum
      do iem = 1, NEMIS_FILE
        emsum(iem)= emsum(iem)+sum(sumemis(:,iem))
      enddo

    enddo

    if(MasterProc)then
      write(*     ,"(a9,3x,30(f12.2,:))")' TOTAL : ',emsum(:)
      write(IO_LOG,"(a9,3x,30(f12.2,:))")' TOTAL : ',emsum(:)
    endif

    !temporary: nlandcode,landcode,snapemis will be completely removed
    nlandcode=nGridEmisCodes
    landcode=GridEmisCodes
    snapemis=GridEmis
    !endif ! EMIS_TEST /Mixed

  case default
    call CheckStop("EMIS_SOURCE not set"//trim(EMIS_SOURCE))
  endselect

  if(USE_ROADDUST) then
    !Use grid-independent Netcdf input files
    call CheckStop(NROAD_FILES>2, "TOO MANY ROADFILES")
    do iem = 1, NROAD_FILES
      !Read data from NetCDF file
      select case(iem)
      case(1);varname='HighwayRoadDustPM10_Jun-Feb'
      case(2);varname='nonHighwayRoadDustPM10_Jun-Feb'
      endselect
      roaddust_emis_pot(:,:,:,iem)=0.0
      call ReadField_CDF('RoadMap.nc',varname,roaddust_emis_pot(1,1,1,iem),&
           nstart=1,interpol='mass_conservative',fractions_out=fractions,&
           CC_out=road_landcode,Ncc_out=road_nlandcode,needed=.true.,&
           debug_flag=.false.,Undef=0.0)
      if(.not.SMI_defined)then
         varname='SMI1'
         call ReadField_CDF('AVG_SMI_2005_2010.nc',varname,SMI,nstart=1,&
              interpol='conservative',needed=.true.,debug_flag=.false.)
         SMI_defined=.true.
      endif

      do i=1,LIMAX
         do j=1,LJMAX
            !Peter: Rough estimate to get something varying between 3.325 (SMI<0.5) and 1.0 (SMI>1)
            SMI_roadfactor=3.325-(min(1.0,max(0.5,SMI(i,j)))-0.5)*2*(3.325-1.0)
            !if(DEBUG_ROADDUST)&
            !  WRITE(*,*)"i,j,RDECF:",i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,SMI_roadfactor
            do iic=road_nlandcode(i,j),1,-1
              roaddust_emis_pot(i,j,iic,iem)=roaddust_emis_pot(i,j,1,iem) &
                   *fractions(i,j,iic)*SMI_roadfactor
            enddo
         enddo
      enddo
      sumroaddust_local(:,iem)=0.0
      do i=1,LIMAX
        do j=1,LJMAX
          do iic=1,road_nlandcode(i,j)
            if(road_landcode(i,j,iic)<=NLAND) &
              ic=find_index(road_landcode(i,j,iic),Country(:)%icode)
            if(Country(ic)%icode/=road_landcode(i,j,iic))then
              write(*,*)"COUNTRY ROAD CODE ERROR: ",road_landcode(i,j,iic),ic,Country(ic)%icode
              call StopAll("COUNTRY CODE ERROR ")
            endif
            if(ic>NLAND)then
              write(*,*)"COUNTRY CODE NOT RECOGNIZED OR UNDEFINED: ",road_landcode(i,j,iic)
              call StopAll("COUNTRY CODE NOT RECOGNIZED ")
            endif
            sumroaddust_local(ic,iem)=sumroaddust_local(ic,iem)&
                 +0.001*roaddust_emis_pot(i,j,iic,iem)
          enddo
        enddo
      enddo
    enddo ! iem = 1, NROAD_FILES-loop
    sumroaddust=0.0
    CALL MPI_REDUCE(sumroaddust_local,sumroaddust,NLAND*NROAD_FILES,MPI_REAL8,&
         MPI_SUM,0,MPI_COMM_CALC,IERROR) 

  endif !USE_ROADDUST

  if(MasterProc) then
    if(USE_ROADDUST)THEN
      call PrintLog("Total road dust emission potentials by countries &
           &(before precipitation and land corrections):")
      write(*     ,"(2a4,11x,30(a12,:))")"  N "," CC ",ROAD_FILE(:)
      write(IO_LOG,"(2a4,11x,30(a12,:))")"  N "," CC ",ROAD_FILE(:)

      do ic = 1, NLAND
        ccsum = sum( sumroaddust(ic,:))
        if(ccsum>0.0) then
          icc=Country(ic)%icode
          write(*     ,"(i4,1x,a4,3x,30(f12.2,:))")icc, Country(ic)%code, sumroaddust(ic,:)
          write(IO_LOG,"(i4,1x,a4,3x,30(f12.2,:))")icc, Country(ic)%code, sumroaddust(ic,:)
        endif
      enddo
    endif ! ROAD DUST
  endif

  ! now all values are read, snapemis is distributed, globnland and 
  ! globland are ready for distribution
  ! print *, "calling glob2local_int for iem", iem, " me ", me
  select case(EMIS_SOURCE)
  case("emislist")
     call StopAll("The emislist option is not available anymore! Use Mixed instead")
  case("CdfFractions")
     ! emissions directly defined into nlandcode,landcode and snapemis
  endselect

  ! Create emislist-type files for both snap emissions and Cdf
  ! Useful for export to other codes, including production of
  ! new emislist for current NWP grid.
  do iem = 1, NEMIS_FILE
    if((EMIS_SOURCE=="Mixed").and.EMIS_OUT) &
      call EmisOut("Snap",iem,nlandcode,landcode,snapemis(:,:,:,:,iem))

    if(EMIS_TEST=="CdfSnap")&
      write(*,"(a,i3,2i4,2es12.3)") "CALLED CDF PRE "//&
        trim(EMIS_FILE(iem)),me,maxval(nlandcode),maxval(nGridEmisCodes), &
        maxval(snapemis(:,:,:,:,iem)),maxval(GridEmis(:,:,:,:,iem))

    if((EMIS_TEST=="CdfSnap").and.EMIS_OUT) &
      call EmisOut("Cdf",iem,nGridEmisCodes,GridEmisCodes,GridEmis(:,:,:,:,iem))
  enddo

  !**  Conversions:
  ! The emission-data file are so far in units of 
  ! tonnes per grid-square. The conversion factor from tonnes per 50*50km2
  ! annual emission values to surface flux (kg/m2/s) is found by division
  ! with (nydays*24*60*60)s and (h*h)m2 and multiply by 1.e+3.
  ! The conversion factor then equals 1.27e-14
  tonne_to_kgm2s  = 1.0e3 / (nydays * 24.0 * 3600.0 * GRIDWIDTH_M * GRIDWIDTH_M)
  if(MYDEBUG.and.MasterProc) then
    write(*,*) "CONV:me, nydays, gridwidth = ",me,nydays,GRIDWIDTH_M
    write(*,*) "No. days in Emissions: ", nydays
    write(*,*) "tonne_to_kgm2s in Emissions: ", tonne_to_kgm2s
    write(*,*) "Emissions sums:"
    do iem = 1, NEMIS_FILE
      write(*,"(a15,f12.2)") EMIS_FILE(iem),emsum(iem)
    enddo
  endif

  iemCO=find_index("co",EMIS_FILE(:)) ! save this index

  if(MYDEBUG.and.debug_proc.and.iemCO>0) &
    write(*,"(a,2es10.3)") "SnapPre:" // trim(EMIS_FILE(iemCO)), &
      sum(snapemis   (:,debug_li,debug_lj,:,iemCO)), &
      sum(snapemis_flat(debug_li,debug_lj,:,iemCO))

  forall (ic=1:NCMAX, j=1:ljmax, i=1:limax, isec=1:NSECTORS,iem=1:NEMIS_FILE)
    snapemis (isec,i,j,ic,iem) = snapemis (isec,i,j,ic,iem) * tonne_to_kgm2s * xm2(i,j)
  endforall

  forall (fic=1:FNCMAX, j=1:ljmax, i=1:limax,iem=1:NEMIS_FILE)
    snapemis_flat(i,j,fic,iem) = snapemis_flat(i,j,fic,iem) * tonne_to_kgm2s * xm2(i,j)
  endforall

  if(MYDEBUG.and.debug_proc.and.iemCO>0) &
    write(*,"(a,2es10.3)") "SnapPos:" // trim(EMIS_FILE(iemCO)), &
      sum(snapemis   (:,debug_li,debug_lj,:,iemCO)), &
      sum(snapemis_flat(debug_li,debug_lj,:,iemCO))

  if(USE_ROADDUST)THEN
    forall (ic=1:NCMAX, j=1:ljmax, i=1:limax, iem=1:NROAD_FILES)
      roaddust_emis_pot(i,j,ic,iem) = &
           roaddust_emis_pot(i,j,ic,iem) * tonne_to_kgm2s * xm2(i,j)
    endforall
  endif !road dust

  err1 = 0
  if(MasterProc) then
    if(USE_ROADDUST)THEN
      deallocate(road_globnland   ,stat=err7)
      deallocate(road_globland    ,stat=err8)
      deallocate(globroad_dust_pot,stat=err9)
      call CheckStop(err7, "De-Allocation error 7 - roadglob")
      call CheckStop(err8, "De-Allocation error 8 - roadglob")
      call CheckStop(err9, "De-Allocation error 9 - roadglob")
    endif
  else
    ! needed for DEBUG=yes compilation options
    if(USE_ROADDUST)THEN
       deallocate(road_globnland,road_globland,globroad_dust_pot,stat=err9)
       call CheckStop(err9, "De-Allocation error 9 - dummy roadglob")
    endif
  endif

  ! now we have nrecmis and can allocate for gridrcemis:
  ! print *, "ALLOCATING GRIDRC", me, NRCEMIS
  allocate(gridrcemis(NRCEMIS,KEMISTOP:KMAX_MID,LIMAX,LJMAX),stat=err1)
  allocate(gridrcemis0(NRCEMIS,KEMISTOP:KMAX_MID,LIMAX,LJMAX),stat=err2)
  call CheckStop(err1, "Allocation error 1 - gridrcemis") 
  call CheckStop(err2, "Allocation error 2 - gridrcemis0")
  if(USE_ROADDUST)THEN
    allocate(gridrcroadd(NROADDUST,LIMAX,LJMAX),stat=err3)
    allocate(gridrcroadd0(NROADDUST,LIMAX,LJMAX),stat=err4)
    call CheckStop(err3, "Allocation error 3 - gridrcroadd")
    call CheckStop(err4, "Allocation error 4 - gridrcroadd0")
  endif
endsubroutine Emissions
!----------------------------------------------------------------------!
!>
!! expandcclist converts e.g. EU28 to indivdual countries
! Only coded for EUMACC2 so far. Should probably use pointers from
! group names.
!----------------------------------------------------------------------!
subroutine expandcclist(xlist, n)
  character(len=*), dimension(:), intent(inout) :: xlist 
  integer, intent(out) ::  n
  character(len=30), dimension(size(xlist)) ::  nlist 
  integer :: i

  nlist(:) = "-"
  n = 1
  CCLIST: do i = 1 , size(xlist)
    !if(MasterProc) print *, "CCNLIST ", me, i, size(xlist), n, xlist(i)
    select case(xlist(i))
    case("EUMACC2")
    !if(MasterProc) print *, "NLIST MACC2 ", me, i, size(EUMACC2), n
      nlist(n:n+size(EUMACC2)-1 ) = (/ EUMACC2 /)
      n=n+size(EUMACC2)
    case("-")
    !if(MasterProc) print *, "NLIST ----- ", me, i, n
      n = n - 1
      exit CCLIST
    case default
    !if(MasterProc) print *, "NLIST DEF - ", me, i, n, xlist(i)
      nlist(n) = xlist(i)
      n=n+1
    endselect
  enddo CCLIST ! i
  xlist(1:n) = nlist(1:n) ! overwrites original
endsubroutine expandcclist
!----------------------------------------------------------------------!
subroutine consistency_check()
!----------------------------------------------------------------------!
!    checks that all the values given so far are consistent
!----------------------------------------------------------------------!
  character(len=30) :: errormsg 
  errormsg = "ok"
  if(size(EMIS_FILE)/=NEMIS_FILE) errormsg = " size EMISNAME wrong "
  call CheckStop(errormsg,"Failed consistency check")
endsubroutine consistency_check
!***********************************************************************
subroutine EmisSet(indate)   !  emission re-set every time-step/hour
!----------------------------------------------------------------------!
! DESCRIPTION:
!   Calculates the emmision-tendencies and the local (instantaneous) dry 
!   deposition in the emission squares.
!   emis set once per hour to allow for day/night variation (and voc 
!   speciation) (based on local time)  for each snap sector.
!   gridrcemis0 calculated every time-step to allow for ps changes.
!   inputs from Emissions in EMISSIONS_ML:
!   country and snap-specific array : 
!          snapemis (NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILES) 
!  
!   Units:
!   snapemis has units of kg/m2/s, SO2 as S, NO2 as N, NH3 as N. 
!   Map factor (xm2) already accounted for. 
!  
!   Data on how many countries contribute per grid square is stored in
!   nlandcode(LIMAX,LJMAX) , with the country-codes given by
!   landcode(LIMAX,LJMAX,NCMAX).
!     
!   Monthly and weekday factors are pre-multiplied and stored in:
!       real timefac(NLAND,NSECTORS,NEMIS_FILES)
!   And day-hour factors in fac_ehh24x7
!----------------------------------------------------------------------!
  implicit none
  type(date), intent(in) :: indate  ! Gives year..seconds
  integer :: i, j, k, f       ! cooridnates, loop variables
  integer :: icc, ncc         ! No. of countries in grid.
  integer :: ficc,fncc        ! No. of countries with
  integer :: iqrc             ! emis indices 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILE
  integer :: itot             ! index in xn()

  ! Save daytime value between calls, initialise to zero
  integer, save, dimension(MAXNLAND) ::  daytime = 0  !  0=night, 1=day
  integer, save, dimension(MAXNLAND) ::  localhour = 1  ! 1-24 local hour in the different countries, ? How to handle Russia, with multiple timezones???
  integer                         ::  hourloc      !  local hour 
  logical                         ::  hourchange   !             
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions

  real ::  ehlpcom,ehlpcom0
  real ::  tfac, dtgrid    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s               ! source term (emis) before splitting
  integer :: iland, iland_timefac  ! country codes, and codes for timefac 
  real ::  sf              ! source term (emis) before splitting  (flat emissions)
  integer :: flat_iland    ! country codes (countries with flat emissions)

  integer, save :: oldday = -1, oldhour = -1
  integer, save :: wday , wday_loc ! wday = day of the week 1-7
  real ::  oldtfac
  logical :: debug_tfac

  ! If timezone=-100, calculate daytime based on longitude rather than timezone
  integer :: daytime_longitude, daytime_iland, hour_longitude, hour_iland,nstart
  integer :: i_Emis_4D
  character(len=125) ::varname 
  TYPE(timestamp)   :: ts1,ts2

  ! Initialize
  ehlpcom0 = GRAV* 0.001*AVOG!0.001 = kg_to_g / m3_to_cm3

  ! Scaling for totemadd:
  dtgrid = dt_advec * GRIDWIDTH_M * GRIDWIDTH_M 

  ! The emis array only needs to be updated every full hour. The 
  ! time-factor calculation needs to know if a local-time causes a shift 
  ! from day to night.  In addition, we reset an overall day's time-factors
  ! at midnight every day. 

  hourchange = (indate%hour/=oldhour).or.(indate%day/=oldday)
  if(hourchange) then
    oldhour = indate%hour
    if(indate%day/=oldday)then
      !==========================
      call NewDayFactors(indate)
      if(USES%DEGREEDAY_FACTORS) call DegreeDayFactors(daynumber) ! => fac_emm, fac_edd
      !==========================
      ! for ROADDUST
      wday=day_of_week(indate%year,indate%month,indate%day)
      if(wday==0)wday=7 ! Sunday -> 7
      oldday = indate%day
    endif

    if(Found_Emis_4D>0)then
      if(.not.allocated(Emis_4D))allocate(Emis_4D(LIMAX,LJMAX,KMAX_MID,N_Emis_4D))
      Emis_4D = 0.0 !default value, must be set to zero when no values are found
      NTime_Read=-1 !read all times    
      call ReadTimeCDF(emis_inputlist(Found_Emis_4D)%Name,TimesInDays,NTime_Read)
      call CheckStop(NTime_Read>size(TimesInDays), "Emissions_ml: increase size of TimesInDays ")
      !if(MasterProc)write(*,*)('found date ',i,TimesInDays(i),i=1,NTime_Read)
      !write(*,*)'compare  ',ts1,ts2
      ts1=make_timestamp(indate)
      do i=1,NTime_Read
        call nctime2date(ts2,TimesInDays(i))   
        if(nint(tdif_secs(ts1,ts2))==0)then
          if(MasterProc)write(*,*)'Emis_4D: found matching date ',i,TimesInDays(i)
          nstart=i
          exit
        endif
      enddo
      if(i>NTime_Read )then
        if(MasterProc)then
          write(*,*)'Emis_4D: WARNING DID NOT FIND ANY MATCHING DATE '
          write(*,*)'Emis_4D: first date found ',TimesInDays(1)
          write(*,*)'Emis_4D: last date found ',TimesInDays(NTime_Read)
          write(*,*)'Emis_4D: difference to last date ',tdif_secs(ts1,ts2)/3600,' hours'
        endif
      else
        do i_Emis_4D=1,N_Emis_4D
          if(emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)=='NOTSET')exit
          varname=emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)
          !if(MasterProc)write(*,*)'Fetching ',trim(varname)
          call GetCDF_modelgrid(varname,emis_inputlist(Found_Emis_4D)%Name,Emis_4D(1,1,1,i_Emis_4D),1,kmax_mid,nstart,1,reverse_k=.true.)
        enddo
      endif
    endif
  endif

  if(DEBUG_EMISTIMEFACS.and.MasterProc) &
    write(*,"(a,2f8.3)") " EmisSet  traffic 24x7", &
      fac_ehh24x7(ISNAP_TRAF,1,4),fac_ehh24x7(ISNAP_TRAF,13,4)
  !..........................................
  !  Look for day-night changes, after local time correction
  !  (daytime(iland) not used if  LONGITUDE_TIME=true)
  do iland = 1, NLAND
    daytime(iland) = 0
    hourloc        = indate%hour + Country(iland)%timezone
    localhour(iland) = hourloc  ! here from 0 to 23
    if(hourloc>=7 .and. hourloc<=18) daytime(iland)=1
  enddo ! iland

  if(hourchange) then 
    totemadd(:)  = 0.
    gridrcemis0(:,:,:,:) = 0.0 
    SumSnapEmis(:,:,:) = 0.0
    if(USE_ROADDUST)gridrcroadd0(:,:,:) = 0.0
    !..........................................
    ! Process each grid:
    do j = 1,ljmax
      do i = 1,limax
        ncc = nlandcode(i,j)            ! No. of countries in grid
        debug_tfac=(DEBUG_EMISTIMEFACS.and.debug_proc.and.i==DEBUG_li.and.j==DEBUG_lj)
        ! find the approximate local time:
        hourloc= mod(nint(indate%hour+24*(1+glon(i,j)/360.0)),24)
        hour_longitude=hourloc
        daytime_longitude=0
        if(hourloc>=7 .and. hourloc<= 18) daytime_longitude=1            
        !*************************************************
        ! First loop over non-flat (one sector) emissions
        !*************************************************
        tmpemis(:)=0.
        do icc = 1, ncc
          !          iland = landcode(i,j,icc)     ! 1=Albania, etc.
          iland=find_index(landcode(i,j,icc),Country(:)%icode) !array index

          !array index of country that should be used as reference for timefactor
          iland_timefac = find_index(Country(iland)%timefac_index,Country(:)%timefac_index)

          if(Country(iland)%timezone==-100)then
            daytime_iland=daytime_longitude
            hour_iland=hour_longitude + 1   ! add 1 to get 1..24 
          else
            daytime_iland=daytime(iland)
            hour_iland=localhour(iland) + 1
          endif
          !if( hour_iland > 24 ) hour_iland = 1 !DSA12
          wday_loc=wday 
          if(hour_iland>24) then
            hour_iland = hour_iland - 24
            wday_loc=wday + 1
            if(wday_loc==0)wday_loc=7 ! Sunday -> 7
            if(wday_loc>7 )wday_loc=1 
          endif
          call CheckStop(hour_iland<1,"ERROR: HOUR Zero in EmisSet")
          if(debug_tfac) then 
          write(*,"(a,i4,2i3,i5,2i4,3x,4i3)") "EmisSet DAYS times ", daynumber, &
            wday, wday_loc, iland, daytime_longitude, daytime_iland,&
            hour_longitude, hour_iland, hourloc, Country(iland)%timezone
          call datewrite("EmisSet DAY 24x7:", &
            (/ icc, iland, wday, wday_loc, hour_iland /), &
            (/ fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc) /) )
          endif
          !  As each emission sector has a different diurnal profile
          !  and possibly speciation, we loop over each sector, adding
          !  the found emission rates to gridrcemis as we go.
          !  ==================================================
          do isec = 1, NSECTORS       ! Loop over snap codes
            ! Calculate emission rates from snapemis, time-factors, 
            ! and if appropriate any speciation fraction (NEMIS_FRAC)
            iqrc = 0   ! index over emisfrac
            do iem = 1, NEMIS_FILE 
              tfac = timefac(iland_timefac,sec2tfac_map(isec),iem) &
                   * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc)


              if(debug_tfac.and.iem==1) &
                write(*,"(a,3i4,f8.3)")"EmisSet DAY TFAC:",isec,sec2tfac_map(isec),hour_iland,tfac

              !it is best to multiply only if USES%GRIDDED_EMIS_MONTHLY_FACTOR
              !in order not to access the array and waste cache if not necessary
              if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)tfac=tfac* GridTfac(i,j,sec2tfac_map(isec),iem)

              !Degree days - only SNAP-2 
              if(USES%DEGREEDAY_FACTORS .and. &
                sec2tfac_map(isec)==ISNAP_DOM .and. Gridded_SNAP2_Factors) then
                oldtfac = tfac
                ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
                ! we make use of a baseload even for SNAP2
                tfac = ( fac_min(iland,sec2tfac_map(isec),iem) & ! constant baseload
                     + ( 1.0-fac_min(iland,sec2tfac_map(isec),iem) )* gridfac_HDD(i,j) ) &
                     * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc)

                if(debug_tfac .and. indate%hour==12 .and. iem==1) &
                  write(*,"(a,3i3,2i4,7f8.3)") "SNAPHDD tfac ",  &
                    isec, sec2tfac_map(isec),iland, daynumber, indate%hour, &
                    timefac(iland_timefac,sec2tfac_map(isec),iem), t2_nwp(i,j,2)-273.15, &
                    fac_min(iland,sec2tfac_map(isec),iem),  gridfac_HDD(i,j), tfac
              endif ! =============== HDD 

              s = tfac * snapemis(isec,i,j,icc,iem)

              ! prelim emis sum kg/m2/s
              SumSnapEmis(i,j,iem) = SumSnapEmis(i,j,iem) + s

              do f = 1,emis_nsplit(iem)
                iqrc = iqrc + 1
                itot = iqrc2itot(iqrc)
                tmpemis(iqrc) = s * emisfrac(iqrc,sec2split_map(isec),iland)
                ! Add up emissions in ktonne 
                totemadd(itot) = totemadd(itot) &
                     + tmpemis(iqrc) * dtgrid * xmd(i,j)
              enddo ! f
            enddo ! iem

            !  Assign to height levels 1-KEMISTOP
            do k=KEMISTOP,KMAX_MID
              do iqrc =1, nrcemis
                gridrcemis0(iqrc,k,i,j) = gridrcemis0(iqrc,k,i,j)   &
                     + tmpemis(iqrc)*ehlpcom0    &
                     *emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) &
                     *emis_masscorr(iqrc)
                !if( debug_tfac.and. iqrc==1 ) then 
                !  write(*,"(a,2i3,2f8.3)") "KPROF ", &
                !    isec, KMAX_BND-k, &
                !    VERTFAC(KMAX_BND-k,sec2hfac_map(isec)),  &
                !    emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))
                !end if
              enddo ! iem
            enddo   ! k
          enddo  ! isec
          !      ==================================================
        enddo ! icc  
        !************************************
        ! Then loop over flat emissions
        !************************************
        tmpemis(:)=0.
        fncc = flat_nlandcode(i,j) ! No. of countries with flat emissions in grid
        do ficc = 1, fncc
          !flat_iland = flat_landcode(i,j,ficc) ! 30=BAS etc.
          flat_iland = find_index(flat_landcode(i,j,ficc),Country(:)%icode) !array index
          if(Country(flat_iland)%is_sea) then  ! saves if statements below
            isec = ISEC_SHIP 
          else
            isec = ISEC_NAT
          endif

          !  As each emission sector has a different diurnal profile
          !  and possibly speciation, we loop over each sector, adding
          !  the found emission rates to gridrcemis as we go.
          !  ==================================================
          !  Calculate emission rates from snapemis, time-factors, 
          !  and if appropriate any speciation  fraction (NEMIS_FRAC)
          iqrc  = 0   ! index over emis
          do iem = 1, NEMIS_FILE 
            sf =  snapemis_flat(i,j,ficc,iem)    
            ! prelim emis sum kg/m2/s
            SumSnapEmis(i,j,iem) = SumSnapEmis(i,j,iem) + sf
            do f = 1,emis_nsplit(iem)
              iqrc = iqrc + 1
              itot = iqrc2itot(iqrc)
              tmpemis(iqrc) = sf * emisfrac(iqrc,sec2split_map(isec),flat_iland)
              ! Add flat emissions in ktonne 
              totemadd(itot) = totemadd(itot) &
                   + tmpemis(iqrc) * dtgrid * xmd(i,j)
            enddo ! f
          enddo ! iem
          ! Assign flat emissions to height levels 1-4. Note, no VERTFAC
          do iqrc =1, nrcemis
            gridrcemis0(iqrc,KMAX_MID,i,j) = gridrcemis0(iqrc,KMAX_MID,i,j) &
                 + tmpemis(iqrc)*ehlpcom0*emis_masscorr(iqrc)
          enddo ! iem
          !      ==================================================
        enddo !ficc 

        if(USE_ROADDUST)then
          ! Limit as in TNO-model (but Lotos/Euros has precip in mm/3h)
          ! In the EMEP case this is in mm/h, so should be equivalent with 2.4mm per day
          NO_PRECIP: if(surface_precip(i,j) < 0.1) then

            ! should use the temporal variation for road dust (SNAP_HOURFAC(HH,7))
            ! and a weekday factor (initially taken from TNO model, could use country
            ! dependent factor in the future)

            ! Temporal variation taken from TNO -> No monthly variation and a single
            ! weekday and diurnal variation (same for all countries)     
            ! -> Need to know day_of_week
            !    Relatively weak variation with day of week so use a simplified approach

            ! if( DEBUG_ROADDUST .and. debug_proc .and. i==DEBUG_li .and. j==DEBUG_lj )THEN
            !    write(*,*)"DEBUG ROADDUST! Dry! ncc=", road_nlandcode(i,j)
            ! endif

            ncc = road_nlandcode(i,j) ! number of countries in grid point
            do icc = 1, ncc    
              !iland = road_landcode(i,j,icc)
              iland = find_index(road_landcode(i,j,icc),Country(:)%icode)
              if(Country(iland)%timezone==-100)then
                hour_iland=hour_longitude+1
              else
                hour_iland=localhour(iland)+1
              endif

              wday_loc = wday ! DS added here also, for fac_ehh24x7
              if( hour_iland > 24 ) then
                hour_iland = 1
                if(wday_loc==0)wday_loc=7 ! Sunday -> 7
                if(wday_loc>7 )wday_loc=1 
              endif

              if(ANY(iland==(/IC_FI,IC_NO,IC_SE/)).and. & ! Nordic countries
                 ANY(indate%month==(/3,4,5/)))then        ! spring road dust
                tfac = fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc)*2.0 ! Doubling in Mar-May (as in TNO model)
              else
                tfac = fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc)
              endif

              do iem = 1, NROAD_FILES
                s = tfac * roaddust_emis_pot(i,j,icc,iem)
                if(DEBUG_ROADDUST.and.debug_proc.and.i==DEBUG_li.and.j==DEBUG_lj)&
                  write(*,*)"DEBUG ROADDUST! iem,tfac,icc,roaddust_emis_pot,s", &
                    iem,tfac,icc,roaddust_emis_pot(i,j,icc,iem),s

                gridrcroadd0(QROADDUST_FI,i,j)=gridrcroadd0(QROADDUST_FI,i,j) &
                     +ROADDUST_FINE_FRAC*s
                gridrcroadd0(QROADDUST_CO,i,j)=gridrcroadd0(QROADDUST_CO,i,j) &
                     +(1.-ROADDUST_FINE_FRAC)*s

                if(all([DEBUG_ROADDUST,debug_proc,i==debug_li,j==debug_lj]))then
                  write(*,*)"gridrcroadfine"  ,gridrcroadd0(QROADDUST_FI,i,j)
                  write(*,*)"gridrcroadcoarse",gridrcroadd0(QROADDUST_CO,i,j)
                endif
              enddo ! nroad files
            enddo   ! icc
            ! should pick the correct emissions (spring or rest of year)
            ! and add the emissions from HIGHWAYplus and NONHIGHWAYS,
            ! using correct fine and coarse fractions.
          else ! precipitation
            gridrcroadd0(:,i,j)=0.
          endif NO_PRECIP
        endif ! ROADDUST
      enddo   ! i
    enddo     ! j
    if(MYDEBUG.and.debug_proc) &    ! emis sum kg/m2/s
      call datewrite("SnapSum, kg/m2/s:"//trim(EMIS_FILE(iemCO)), &
         (/ SumSnapEmis(debug_li,debug_lj,iemCO)  /) )

  endif ! hourchange 

  ! We now scale gridrcemis to get emissions in molecules/cm3/s
  do k= KEMISTOP, KMAX_MID
    do j = 1,ljmax
      do i = 1,limax
        ehlpcom= roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))
        !RB: This should also be done for the road dust emissions
        do iqrc =1, NRCEMIS
          gridrcemis(iqrc,k,i,j) =  gridrcemis0(iqrc,k,i,j)* ehlpcom
        enddo ! iqrc
      enddo   ! i
    enddo     ! j
  enddo       ! k

  if(USE_ROADDUST)THEN
    if(DEBUG_ROADDUST.and.debug_proc) &
      write(*,*)"Before the unit scaling",gridrcroadd(1:2,DEBUG_li,DEBUG_lj)
    do j = 1,ljmax
      do i = 1,limax
        ehlpcom= roa(i,j,KMAX_MID,1)/(ps(i,j,1)-PT)
        do iqrc =1, NROADDUST
          gridrcroadd(iqrc,i,j) =  gridrcroadd0(iqrc,i,j)* ehlpcom * ehlpcom0 &
              * roaddust_masscorr(iqrc)
        enddo ! iqrc
      enddo   ! i
    enddo     ! j
    if(DEBUG_ROADDUST.and.debug_proc) &
      write(*,*)"After the unit scaling",gridrcroadd(1:2,DEBUG_li,DEBUG_lj)
  endif
endsubroutine EmisSet
!***********************************************************************
subroutine newmonth
!----------------------------------------------------------------------!
! DESCRIPTION:
!   Reads in natural DMS emissions at start of each month. Update
!   landcode and nlandcode arrays as needed.
!
!   Reads in snow cover at start of each month. 
!
!   April 2010: read monthly aircraft NOx emissions
!----------------------------------------------------------------------!
  use AirEmis_ml, only : airn
  use ModelConstants_ml, only : KCHEMTOP, KMAX_MID
  use NetCDF_ml, only : ReadField_CDF

  integer i, j,k, iyr, iemislist
  integer n, flat_ncmaxfound         ! Max. no. countries w/flat emissions
  real :: rdemis(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
  character(len=200) :: fname
  real ktonne_to_kgm2s, tonnemonth_to_kgm2s  ! Units conversion
  integer :: IQSO2                   ! Index of sox in  EMIS_FILE
  integer errcode,iland
  integer :: iem,ic,isec, i_gridemis
  real :: conv
  logical , save :: first_call=.true.
  logical :: needed_found

  ! For now, only the global runs use the Monthly files
  integer :: kstart,kend,nstart,Nyears
  real :: buffer(LIMAX,LJMAX),SumSoilNOx,ccsum
  real :: fractions(LIMAX,LJMAX,NCMAX),Reduc(NLAND)
  real, dimension(NEMIS_FILE)       :: emsum ! Sum emis over all countries
  real, dimension(NLAND,NEMIS_FILE) :: sumemis, sumemis_local ! Sum of emissions per country
  character(len=40) :: varname 
  character(len=125) ::fileName
  real :: Mask_ReducFactor
  integer :: NMask_Code,Mask_Code(NLAND), i_femis_lonlat
  real :: lonlat_fac, dms_sum

  if(.not.allocated(airn).and.(USE_LIGHTNING_EMIS.or.USE_AIRCRAFT_EMIS))&
    allocate(airn(KCHEMTOP:KMAX_MID,LIMAX,LJMAX))


  if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)&
    call Read_monthly_emis_grid_fac(current_date%month)

  if(USE_AIRCRAFT_EMIS)then
    airn = 0.0
    kstart=KCHEMTOP
    kend=KMAX_MID

    call ReadField_CDF_FL('AircraftEmis_FL.nc','NOx',airn,&
         current_date%month,kstart,kend,&
         interpol='mass_conservative', needed=.true.,debug_flag=.false.)

    ! convert from kg(NO2)/month into molecules/cm3/s
    ! from kg to molecules: 0.001*AVOG/species(NO2)%molwt
    ! use roa to find dz for consistency with other emissions 
    ! (otherwise could have used z_bnd directly)
    ! dz=dP/(roa*GRAV)  dP=dA(k) + dB(k)*ps(i,j,1)
    ! dV=dz*dx*dy=dz*gridwidth**2/xm**2 *1e6 (1e6 for m3->cm3)
    ! from month to seconds: ndaysmonth*24*3600

    conv=0.001*AVOG/species(NO2)%molwt*GRAV/gridwidth_m**2*1.0e-6&
         /(nmdays(current_date%month)*24*3600)

    do k=KCHEMTOP,KMAX_MID
      do j=1,ljmax
        do i=1,limax
          airn(k,i,j)=airn(k,i,j)*conv*(roa(i,j,k,1))&
               /(dA(k)+dB(k)*ps(i,j,1))*xm2(i,j)
        enddo
      enddo
    enddo
  endif

  if(USE_EURO_SOILNOX)then  ! European Soil NOx emissions
    if(DEBUG_SOILNOX.and.debug_proc) write(*,*)"Emissions DEBUG_SOILNOX START"

    ! read in map of annual N-deposition produced from pre-runs of EMEP model
    ! with script mkcdo.annualNdep
    call ReadField_CDF('annualNdep.nc','Ndep_m2',AnnualNdep,1,&
          interpol='zero_order',needed=.true.,debug_flag=.false.,UnDef=0.0)

    if(DEBUG_SOILNOX.and.debug_proc)&
      write(*,"(a,4es12.3)") "Emissions_ml: SOILNOX AnnualDEBUG ", &
        AnnualNdep(debug_li, debug_lj), maxval(AnnualNdep), minval(AnnualNdep)

    call CheckStop(USE_GLOBAL_SOILNOX, "SOILNOX - cannot use global with Euro")
    ! We then calculate SoulNOx in Biogenics_ml

  elseif(USE_GLOBAL_SOILNOX) then ! Global soil NOx

    SoilNOx(:,:)=0.0      
    buffer(:,:)=0.0      
    nstart=(current_date%year-1996)*12 + current_date%month
    if(nstart>0.and.nstart<=120)then
      !the month is defined
      call ReadField_CDF('nox_emission_1996-2005.nc','NOX_EMISSION',SoilNOx,&
             nstart=nstart,interpol='conservative',known_projection="lon lat",&
             needed=.true.,debug_flag=.false.,UnDef=0.0)
      if(DEBUG_SOILNOX.and.debug_proc) &
        write(*,*) "PROPER YEAR of SOILNO ", current_date%year, nstart
    else
      !the year is not defined; average over all years
      Nyears=10 !10 years defined
      do iyr=1,Nyears 
        nstart=12*(iyr-1) + current_date%month  
        call ReadField_CDF('nox_emission_1996-2005.nc','NOX_EMISSION',buffer,&
              nstart=nstart,interpol='conservative',known_projection="lon lat",&
              needed=.true.,debug_flag=.false.,UnDef=0.0)
        do j=1,ljmax 
          do i=1,limax
            SoilNOx(i,j)=SoilNOx(i,j)+buffer(i,j)
          enddo
        enddo
        if(DEBUG_SOILNOX.and.debug_proc) &
          write(*,"(a,2i6,es10.3,a,2es10.3)") "Averaging SOILNO  inputs", &
            1995+(iyr-1), nstart,SoilNOx(debug_li, debug_lj), &
            "max: ", maxval(buffer), maxval(SoilNOx)
      enddo
      SoilNOx=SoilNOx/Nyears
    endif ! nstart test

     if(DEBUG_SOILNOX.and.debug_proc) then
        write(*,"(a,i3,3es10.3)") "After Global SOILNO ",&
             me,maxval(SoilNOx),SoilNOx(debug_li,debug_lj)
       !write(*,"(a,i3,3es10.3)") "After Global SOILNO ",  me, maxval(SoilNOx), SoilNOx(3, 3)
     endif
  else ! no soil NO
    if(DEBUG_SOILNOX.and.debug_proc) &
      write(*,*) "Emissions DEBUG_SOILNOX - none"
  endif !  SOIL NO

  !for testing, compute total soil NOx emissions within domain
  !convert from g/m2/day into kg/day
  if(USE_GLOBAL_SOILNOX) then 
    SumSoilNOx=0.0
    SoilNOx = max(0.0, SoilNOx)  ! Stops the NEGs!
    do j=1,ljmax
      do i=1,limax      
        SumSoilNOx=SumSoilNOx+0.001*SoilNOx(i,j)*gridwidth_m**2*xmd(i,j)      
      enddo
    enddo
    CALL MPI_ALLREDUCE(SumSoilNOx,mpi_out,1,MPI_DOUBLE_PRECISION, &
         MPI_SUM,MPI_COMM_CALC,IERROR) 
    SumSoilNOx = mpi_out
    if(MasterProc)&
      write(*,*)'GLOBAL SOILNOX emissions this month within domain',&
        SumSoilNOx,' kg per day'

    ! convert from g(N)/m2/day into molecules/cm3/s from g to molecules:
    !  AVOG/14  14=molweight N, use roa to find dz for consistency with other
    !  emissions (otherwise could have used z_bnd directly) dz=dP/(roa*GRAV)
    !  dP=dA(k) + dB(k)*ps(i,j,1) dV=dz*1e6 (1e6 for m3->cm3) from month to
    !  seconds: ndaysmonth*24*3600
    conv=AVOG/14.0*GRAV*1.0e-6/(24*3600)
    k=KMAX_MID!surface
    do j=1,ljmax
      do i=1,limax      
        SoilNOx(i,j)=SoilNOx(i,j)*conv*roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))
      enddo
    enddo
  endif

  ! DMS
  !   Units:
  !   Input files seem to be in ktonne PER YEAR. We convert here to kg/m2/s
  !   The conversion factor from 50*50km2
  !   annual emission values to surface flux (kg/m2/s) is found by division
  !   with (nydays*24*60*60)s and (h*h)m2 and multiply by 1.e+6.
  !   the conversion factor (ktonne_to_kgm2s) then equals 1.27e-8 
  !   NB: a new file is read every month; this means that total emissions 
  !   are NOT the sum of the 12 files emissions (but about 12 times less 
  !   than the sum). 
  !   More precisely: year_emis=sum_months(emis_month*nmdays/nydays)

  ktonne_to_kgm2s = 1.0e6/(nydays*24.*60.*60.*GRIDWIDTH_M*GRIDWIDTH_M)

  if(MasterProc.and.MYDEBUG) then
    write(*,*) 'Enters newmonth, mm, ktonne_to_kgm2s = ', &
         current_date%month,ktonne_to_kgm2s
    write(*,*) ' first_dms_read = ', first_dms_read
  endif
  !...........................................................................
  !        DMS Input - land 35 - SNAP sector 11
  !...........................................................................
  flat_ncmaxfound = 0  ! Max. no. countries(w/flat emissions) per grid
  !  Natural SO2 emissions
  IQSO2=find_index("sox",EMIS_FILE(:))
  if(IQSO2<1) then
    write(*,*) " No SO2 emissions - need to skip DMS also"
    return     ! No need to read DMS fields 
  else    
    ! We have so2 emission so need DMS also
    if(.not.first_call)then
      !write some diagnostic for the past month emissions
      CALL MPI_ALLREDUCE(O_DMS%sum_month, mpi_out, 1,&
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
      O_DMS%sum_month = mpi_out
      if(MasterProc)then
59    format(A,6F14.5)
      write(*,*)'DMS OCEAN emissions '
      write(*,59)'SO2 from ocean DMS cdf file ',O_DMS%sum_month
      write(*,59)'SO2 from natso2.dat ',DMS_natso2_month
      write(*,59)'fraction new/old method',O_DMS%sum_month/DMS_natso2_month
      endif
      O_DMS%sum_year=O_DMS%sum_year+O_DMS%sum_month
      O_DMS%sum_month=0.0
      !natso2 already allreduced
      DMS_natso2_year=DMS_natso2_year+DMS_natso2_month
      DMS_natso2_month=0.0
    endif

    write(fname,'(''natso2'',i2.2,''.dat'')')current_date%month
    if(MasterProc)write(*,*) 'Reading DMS emissions from ',trim(fname)
    needed_found=.false.
    call ReadField(IO_DMS,fname,rdemis,needed_found)

    if(needed_found)then
    errcode = 0
    dms_sum=0.0
    do j=1,ljmax
      do i=1,limax                     
        ! Add DMS for country code IQ_DMS=35  to snap sector 11=Nature.
        ! First time we read we must add DMS to the "countries" 
        ! contributing within the grid square.                      
        ! - for flat emissions:                     
        if(first_dms_read) then 
          flat_nlandcode(i,j) = flat_nlandcode(i,j) + 1 
          n = flat_nlandcode(i,j)
          flat_landcode(i,j,n) = IQ_DMS !=Country(IC_NAT)%icode   ! IQ_DMS country code index 35 
          if(n>flat_ncmaxfound) then
            flat_ncmaxfound = n 
            if (MYDEBUG) write(6,*)'DMS Increased flat_ncmaxfound to ',n 
            call CheckStop( n > FNCMAX, "IncreaseFNCMAX for dms")
          endif
        else  ! We know that DMS lies last in the array, so:
          n = flat_nlandcode(i,j)
          call CheckStop(flat_landcode(i,j,n),IQ_DMS,"Newmonth:DMS not last!")
        endif
        if(.not.USE_OCEAN_DMS)then
        snapemis_flat(i,j,n,IQSO2) = rdemis(i,j) * ktonne_to_kgm2s * xm2(i,j)
        endif
        dms_sum=dms_sum+rdemis(i,j)
      enddo ! i
    enddo   ! j
    CALL MPI_ALLREDUCE(dms_sum, mpi_out, 1,&
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
    dms_sum=mpi_out

    DMS_natso2_month=dms_sum*nmdays(current_date%month)/nydays

    if(MasterProc) then
      write(*,*)'total dms emission from ',trim(fname),dms_sum,DMS_natso2_year+DMS_natso2_month
      write(*,*)' Gg this month ',DMS_natso2_month
    endif
    if(first_dms_read) then
      if(MYDEBUG) &
        write(*,*)'me ',me, ' Increased flat_ncmaxfound to ',flat_ncmaxfound 
      first_dms_read = .false.
    endif
    else!no dms file found
      call PrintLog("WARNING: NO DMS emissions found",MasterProc)
    endif
  endif  ! IQSO2>0

  sumemis=0.0
  do iemislist = 1, size( emis_inputlist(:)%name )

    if(emis_inputlist(iemislist)%name == "NOTSET")cycle
     
    !periodicity set in routine Emissions(year) if 12 records are found
    if(emis_inputlist(iemislist)%periodicity /= "monthly")cycle
    if(MasterProc)write(*,*)'reading monthly emissions for ',trim(emis_inputlist(iemislist)%name) 
     
    if(emis_inputlist(iemislist)%type == "sectors")then !sectors is default
    !Read monthly emission snap sector files

    !reset emissions (except flat emissions!)
    snapemis(:,:,:,:,:) = 0.
    fractions(:,:,:)    = 0.
    cdfemis(:,:)       = 0.
    sumemis_local       = 0.
    emsum               = 0.


    tonnemonth_to_kgm2s = 1.0e3 &
         /(nmdays(current_date%month)*24.*60.*60.*GRIDWIDTH_M*GRIDWIDTH_M)

    do iem = 1,NEMIS_FILE
       if(emis_inputlist(iemislist)%pollName(1)/='NOTSET')then
          if(emis_inputlist(iemislist)%pollName(1)/=trim(EMIS_FILE(iem)))cycle      
       endif
       do  isec = 1,NSECTORS            
        !define mask (can be omitted if not sent to readfield) 
        NMask_Code=0
        !example if code reductions defined from "femis.dat" shoudl be used to define mask countries
        do iland = 1,NLAND
          if (abs(e_fact(isec,iland,iem)-1.0)>0.000001)then
            NMask_Code=NMask_Code+1
            Mask_Code(NMask_Code) = iland!codes from mask file to be reduced
            Mask_ReducFactor=e_fact(isec,iland,iem)!NB: will only take the last defined value!
            !if(MasterProc)write(*,*)'maskcode ',iland,isec,iem,Mask_ReducFactor
          endif
        enddo

        !example if hardcoded definitions are to be used. Here multiply country_code=18 emissions with 0.8
        !Mask_ReducFactor=0.8
        !NMask_Code=1
        !Mask_Code(1)=18

        !NB: "traditional" femis.dat can be used on top of MASK reductions. 
        !    Country codes from Mask_Code will be applied for Mask, and country codes from Reduc 
        !    will be applied to countries defined in emissions file 

        Reduc=e_fact(isec,:,iem) 
        fileName='EmisFracs.nc'
        if(EMIS_SOURCE=="Mixed") fileName=fileName_monthly!will be default in the future
        write(varname,"(A,I2.2)")trim(EMIS_FILE(iem))//'_sec',isec
        call  ReadField_CDF(trim(fileName),varname,cdfemis(1,1),nstart=current_date%month,&
             interpol='mass_conservative',fractions_out=fractions,&
             CC_out=landcode,Ncc_out=nlandcode,Reduc=Reduc,needed=.true., & 
             Mask_fileName='sourceNC.nc',Mask_varname='source_code',Mask_Code=Mask_Code,&
             NMask_Code=NMask_Code,Mask_ReducFactor=Mask_ReducFactor,&
             debug_flag=.false.,Undef=0.0)

        do j=1,ljmax
          do i=1,limax                   
            if(nlandcode(i,j)>NCMAX)then
              write(*,*)"To many emitter countries in one gridcell: ",&
                   me,i,j,nlandcode(i,j)
              call StopAll("To many countries in one gridcell ")
            endif
            lonlat_fac=1.0
            if(N_femis_lonlat>0)then
              do i_femis_lonlat=1,N_femis_lonlat
                if(glat(i,j)>femis_latmin(i_femis_lonlat).and.&
                   glat(i,j)<femis_latmax(i_femis_lonlat).and.&
                   glon(i,j)>femis_lonmin(i_femis_lonlat).and.&
                   glon(i,j)<femis_lonmax(i_femis_lonlat)     )then
                  lonlat_fac=lonlat_fac*e_fact_lonlat(isec,i_femis_lonlat,iem) 
                endif
              enddo
            endif
            do n=1,nlandcode(i,j)
              ic=find_index(landcode(i,j,n),Country(:)%icode)                     
              !exclude or include countries
              !could easily be optimized, by defining a country mask
              if(nin_monthly>0)then
                !1) check that country is in include list
                found=find_index(Country(ic)%code ,incl_monthly(1:nin_monthly),first_only=.true.)
                if(found<=0)cycle!do not include
              endif
              if(nex_monthly>0)then
                !1) check that country is not in exclude list
                found=find_index(Country(ic)%code ,excl_monthly(1:nex_monthly),first_only=.true.)
                if(found>0)cycle!exclude
              endif
              snapemis(isec,i,j,n,iem)=snapemis(isec,i,j,n,iem)&
                   +fractions(i,j,n)*cdfemis(i,j)*lonlat_fac*tonnemonth_to_kgm2s*xm2(i,j)                

              if(Country(ic)%icode/=landcode(i,j,n))then
                write(*,*)"COUNTRY CODE ERROR: ",landcode(i,j,n),ic,Country(ic)%icode
                call StopAll("COUNTRY CODE ERROR ")
              endif
              if(ic>NLAND)then
                write(*,*)"COUNTRY CODE NOT RECOGNIZED OR UNDEFINED: ",landcode(i,j,n)
                call StopAll("COUNTRY CODE NOT RECOGNIZED ")
              endif

              sumemis_local(ic,iem)=sumemis_local(ic,iem)&
                   +0.001*snapemis(isec,i,j,n,iem)*lonlat_fac/(tonnemonth_to_kgm2s*xm2(i,j))!for diagnostics, mass balance
            enddo
          enddo
        enddo
      enddo!sectors

      CALL MPI_REDUCE(sumemis_local(1,iem),sumemis(1,iem),NLAND,MPI_REAL8,&
                      MPI_SUM,0,MPI_COMM_CALC,IERROR) 
      emsum(iem)= sum(sumemis(:,iem))
      if(MasterProc)write(*,"(a,f12.2)")&
        trim(EMIS_FILE(iem))//' monthly TOTAL (tonnes):',emsum(iem)

      if(EMIS_SOURCE=="Mixed".and.emis_inputlist(2)%name /= "NOTSET")then
        !merge with other emissions
        do j=1,ljmax
          do i=1,limax                   
            do i_gridemis=1,nGridEmisCodes(i,j)
              ic=find_index(GridEmisCodes(i,j,i_gridemis),Country(:)%icode)

              !merge other (already read in) emissions into snapemis
              Cexist=.false.
              do n=1,nlandcode(i,j)
                if(GridEmisCodes(i,j,i_gridemis)==landcode(i,j,n))then                           
                  snapemis(isec,i,j,n,iem)=snapemis(isec,i,j,n,iem)&
                       +GridEmis(isec,i,j,i_gridemis,iem)      

                  Cexist=.true.
                  exit
                endif
              enddo
              if(.not.Cexist)then
                !country not included yet. define it now:
                nlandcode(i,j)=nlandcode(i,j)+1
                if(nlandcode(i,j)>NCMAX)then
                  write(*,*)"Too many emitter countries in one gridemiscell: ",&
                       me,i,j,nGridEmisCodes(i,j)
                  call StopAll("To many countries in one gridemiscell ")
                endif
                n=nlandcode(i,j)
                landcode(i,j,n)=GridEmisCodes(i,j,i_gridemis)
                snapemis(:,i,j,n,iem)=snapemis(:,i,j,n,iem)+GridEmis(:,i,j,i_gridemis,iem)
              endif
!             sumemis_local(ic,iem)=sumemis_local(ic,iem)&
!               +0.001*fractions(i,j,n)*cdfemis(i,j)!for diagnostics, mass balance
            enddo
          enddo
        enddo
      endif

    enddo! iem = 1,NEMIS_FILE

  elseif(emis_inputlist(iemislist)%type == 'OceanNH3')then
    if(MasterProc)write(*,*)'reading OceanNH3'
    call ReadField_CDF(emis_inputlist(iemislist)%name,'emiss_ocn',O_NH3%emis,&
         nstart=current_date%month,interpol='conservative',known_projection="lon lat",&
         needed=.true.,debug_flag=.false.,UnDef=0.0)
 
    !diagnostics:
    !call printcdf('OceanNH3',OceanNH3,'kg/m2/s')
    !convert from kg/m2/s into kg/month/subdomain
    O_NH3%sum_month=0.0
    do j=1,ljmax
      do i=1,limax            
        O_NH3%sum_month=O_NH3%sum_month+O_NH3%emis(i,j)&
             *gridwidth_m**2*xmd(i,j)*3600*24*nmdays(current_date%month)             
      enddo
    enddo
    !sum all subdomains
    CALL MPI_ALLREDUCE(O_NH3%sum_month, mpi_out, 1,&
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
    O_NH3%sum_month = mpi_out
    O_NH3%sum_month=O_NH3%sum_month*1e-6 !kg->Gg

    USE_OCEAN_NH3=.true.
    if(MasterProc)write(*,*)'Total monthly NH3 from Oceans (in Gg) ',O_NH3%sum_month
    O_NH3%sum_year=O_NH3%sum_year+O_NH3%sum_month!total for all month

    else if(emis_inputlist(iemislist)%type == 'DMS')then
      if(MasterProc)write(*,*)'reading DMS'
      call ReadField_CDF(emis_inputlist(iemislist)%name,'DMS_sea',O_DMS%emis,&
           nstart=current_date%month,interpol='conservative',known_projection="lon lat",&
           needed=.true.,debug_flag=.false.,UnDef=0.0)

      USE_OCEAN_DMS=.true.
      FOUND_OCEAN_DMS=.true.
      
      !from nanomol/l -> mol/cm3
      O_DMS%emis=O_DMS%emis*1.0e-12

    else !
      call StopAll("Monthly emission type not implemented "//trim(emis_inputlist(iemislist)%type))
    endif

  enddo !iemislist

  if(MasterProc)then
    do ic = 1, NLAND
      ccsum = sum( sumemis(ic,:) )
      !if ( ccsum > 0.0 ) then
      if(ccsum>0.0 .and. MasterProc ) then
        write(*,"(i5,1x,a10,3x,30(f12.2,:))")ic, Country(ic)%code, sumemis(ic,:)
        write(IO_LOG,"(i5,1x,a10,3x,30(f12.2,:))")ic, Country(ic)%code, sumemis(ic,:)
      endif
    enddo
  endif
  first_call=.false.
endsubroutine newmonth
!***********************************************************************
subroutine EmisOut(label, iem,nsources,sources,emis)
!----------------------------------------------------------------------!
! Ascii output of emissions fields (after interpolation and femis.
! To print out emissions, we need to get fields for each country
! in turn, needing info from GridCodes
! e.g. 11-SNAP data is snapemis:
!  real    snapemis (NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILES)
  character(len=*) :: label
  integer, intent(in) :: iem
  real,intent(in), dimension(:,:,:,:):: emis      ! Emission values
  integer, dimension(:,:), intent(in) :: nsources      ! Emission values
  integer, dimension(:,:,:), intent(in) :: sources      ! Emission values
  character(len=100) :: txt
  real :: low, high
  integer :: msg(4), i,j, ii, jj, isec, icc, ncc, iland, iproc
  real ::locemis(LIMAX, LJMAX,NSECTORS)
  real ::lemis(LIMAX, LJMAX)
  txt = trim(label)//"."//trim(EMIS_FILE(iem))
  msg(:) = 0

  if(MYDEBUG) write(*,*)"CALLED "//trim(txt),me,&
    maxval(emis),maxval(nsources),maxval(sources)

!  allocate(locemis(LIMAX, LJMAX,NSECTORS), stat=msg(1) )
!  allocate(lemis(LIMAX, LJMAX), stat=msg(2) )

  call CheckStop(any(msg(:)/=0),"EmisOut alloc error:"//trim(label))

  if(MasterProc)write(*,*)' WARNING - i,j coordinates are not consecutive - dependent on number of processors'
!Each subdomain output its data, one at a time. The others MUST wait.
  do iproc=1,NPROC
    if(me==iproc-1)then
      if(MasterProc)then
        open(IO_TMP,file="EmisOut"//trim(txt))!new file
      else
        open(IO_TMP,file="EmisOut"//trim(txt),access='append')!append
      endif
      EMLAND: do iland = 1, NLAND
        locemis = 0.0
        !    print *,  trim(txt)//" iland ", me, iland, maxval(emis(:,:,:,:))
        
        !/ Collect emissions locally by country code iland
        ncc=-999
        do j = 1, ljmax
          do i = 1, limax
            ncc = nsources(i,j)
            do icc = 1, ncc
              if(sources(i,j,icc)==iland) then
                locemis(i,j,: ) = emis(:, i,j,icc)
                if(MYDEBUG) call CheckStop(any(locemis(i,j,:)< 0.0),"NEG LOCEMIS")
              endif
            enddo
          enddo
        enddo ! j
        if(ncc==-999)cycle!ncountry not in this subdomain
        ! Should never happen, but...
        !call CheckStop( any( lemis < 0.0 ) , "NEG LEMIS")
        
        do i = 1,limax
          do j = 1,ljmax
            ii=gi0+i+IRUNBEG-2
            jj=gj0+j+JRUNBEG-2
            if(sum(locemis(i,j,:))>1.0e-10)then
              high = locemis(i,j,1)
              low =  sum(locemis(i,j,2:NSECTORS))
              write(IO_TMP,"(i3,2i4,2x,13es10.3)") iland, ii,jj, &
                   low, high, (locemis(i,j,isec),isec=1,NSECTORS)
            endif
          enddo
        enddo
      enddo EMLAND
      
      do isec = 1, NSECTORS
         lemis = locemis(:,:,isec)
         if(MYDEBUG.and.debug_proc) write(*,*) trim(txt)//" lemis ",me,iland,isec,maxval(lemis(:,:))
      enddo ! isec
      
    endif
    close(IO_TMP)
    CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)!wait: one should write at a time
  enddo
  
!  deallocate(locemis,lemis)

endsubroutine EmisOut

subroutine uemep_emis(indate)

  implicit none
  type(date), intent(in) :: indate  ! Gives year..seconds
  integer :: i, j, k          ! coordinates, loop variables
  integer :: icc, ncc         ! No. of countries in grid.
  integer :: ficc,fncc        ! No. of countries with
  integer :: iqrc             ! emis indices 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILE
  integer :: itot             ! index in xn()

  ! Save daytime value between calls, initialise to zero
  integer, save, dimension(MAXNLAND) ::  daytime = 0  !  0=night, 1=day
  integer, save, dimension(MAXNLAND) ::  localhour = 1  ! 1-24 local hour in the different countries, ? How to handle Russia, with multiple timezones???
  integer                         ::  hourloc      !  local hour 
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions
  real ::  tfac    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s       ! source term (emis) before splitting
  integer :: iland, iland_timefac  ! country codes, and codes for timefac 
  integer :: daytime_longitude, daytime_iland, hour_longitude, hour_iland,nstart
  integer ::icc_uemep
  integer, save :: wday , wday_loc ! wday = day of the week 1-7
  integer ::ix,iix
  real::dt_uemep, xtot, emis_uemep(KMAX_MID),emis_tot(KMAX_MID)
  logical,save :: first_call=.true.  

  if(first_call)then
    !init uemep
!   uEMEP%emis="pm25"!one of the emission: pm25, sox, nox, voc, pmco, nh3 or co.
!   uEMEP%sector=0!0 = all sectors, or choose one ector

    iem=find_index(uEMEP%emis ,EMIS_FILE(1:NEMIS_FILE))
    call CheckStop( iem<1, "uEMEP did not find corresponding emission file: "//trim(uEMEP%emis) )
    
    uEMEP%Nix=emis_nsplit(iem)
    call CheckStop( uEMEP%Nix>size(uEMEP%ix), "uEMEP: increase size of uEMEP%ix()!" )

    do i=1,uEMEP%Nix
      iqrc=sum(emis_nsplit(1:iem-1)) + i
      itot=iqrc2itot(iqrc)
      uEMEP%ix(i)=itot-NSPEC_SHL
    enddo
    if(MasterProc)then
      write(*,*)'uEMEP sector: ',uEMEP%sector
      write(*,*)'uEMEP emission file: ',uEMEP%emis
      write(*,*)'uEMEP number of species in group: ',uEMEP%Nix
      write(*,"(30A)")'including: ',(trim(species_adv(uEMEP%ix(i))%name),' ', i=1,uEMEP%Nix)
    endif

  endif

  dt_uemep=dt_advec

  wday=day_of_week(indate%year,indate%month,indate%day)
  if(wday==0)wday=7 ! Sunday -> 7
  do iland = 1, NLAND
    daytime(iland) = 0
    hourloc        = indate%hour + Country(iland)%timezone
    localhour(iland) = hourloc  ! here from 0 to 23
    if(hourloc>=7 .and. hourloc<=18) daytime(iland)=1
  enddo ! iland
   
  do j = lj0,lj1
    do i = li0,li1
      ncc = nlandcode(i,j)            ! No. of countries in grid
      fncc = flat_nlandcode(i,j) ! No. of countries with flat emissions in grid
      hourloc= mod(nint(indate%hour+24*(1+glon(i,j)/360.0)),24)
      hour_longitude=hourloc
      daytime_longitude=0
      if(hourloc>=7 .and. hourloc<= 18) daytime_longitude=1            
      !*************************************************
      ! First loop over non-flat (one sector) emissions
      !*************************************************
      tmpemis(:)=0.
      icc_uemep=0
      emis_uemep=0.0
      emis_tot=0.0
      do icc = 1, ncc+fncc
        ficc=icc-ncc
        !          iland = landcode(i,j,icc)     ! 1=Albania, etc.
        if(icc<=ncc)then
          iland=find_index(landcode(i,j,icc),Country(:)%icode) !array index
        else
          iland=find_index(flat_landcode(i,j,ficc),Country(:)%icode) 
        endif
        !array index of country that should be used as reference for timefactor
        iland_timefac = find_index(Country(iland)%timefac_index,Country(:)%timefac_index)

        if(Country(iland)%timezone==-100)then
          daytime_iland=daytime_longitude
          hour_iland=hour_longitude + 1   ! add 1 to get 1..24 
        else
          daytime_iland=daytime(iland)
          hour_iland=localhour(iland) + 1
        endif
        !if( hour_iland > 24 ) hour_iland = 1 !DSA12
        wday_loc=wday 
        if(hour_iland>24) then
          hour_iland = hour_iland - 24
          wday_loc=wday + 1
          if(wday_loc==0)wday_loc=7 ! Sunday -> 7
          if(wday_loc>7 )wday_loc=1 
        endif

        do iem = 1, NEMIS_FILE 
          if(trim(EMIS_File(iem))/=trim(uEMEP%emis))cycle
          do isec = 1, NSECTORS       ! Loop over snap codes
            ! Calculate emission rates from snapemis, time-factors, 
            ! and if appropriate any speciation fraction (NEMIS_FRAC)
            iqrc = 0   ! index over emisfrac
            ! kg/m2/s
            
            if(icc<=ncc)then
              tfac = timefac(iland_timefac,sec2tfac_map(isec),iem) &
                   * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc)

              !Degree days - only SNAP-2 
              if(USES%DEGREEDAY_FACTORS .and. &
                   sec2tfac_map(isec)==ISNAP_DOM .and. Gridded_SNAP2_Factors) then
                 ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
                 ! we make use of a baseload even for SNAP2
                 tfac = ( fac_min(iland,sec2tfac_map(isec),iem) & ! constant baseload
                      + ( 1.0-fac_min(iland,sec2tfac_map(isec),iem) )* gridfac_HDD(i,j) ) &
                      * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc)
              endif ! =============== HDD 
              
              s = tfac * snapemis(isec,i,j,icc,iem)
            else
              s = snapemis_flat(i,j,ficc,iem)                        
            endif

            do k=KEMISTOP,KMAX_MID
              emis_tot(k)=emis_tot(k)+s*emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))*dt_uemep
            enddo

            if(isec==uEMEP%sector .or. uEMEP%sector==0)then
              do k=KEMISTOP,KMAX_MID
                emis_uemep(k)=emis_uemep(k)+s*emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))*dt_uemep
              enddo
            endif

          enddo ! iem

        enddo  ! isec
        !      ==================================================
      enddo ! icc  
      
      do k=KEMISTOP,KMAX_MID
        if(emis_tot(k)<1.E-20)cycle
        !units kg/m2
        !total pollutant
        xtot=0.0
        do iix=1,uEMEP%Nix
          ix=uEMEP%ix(iix)
          xtot=xtot+(xn_adv(ix,i,j,k)*species_adv(ix)%molwt)*(dA(k)+dB(k)*ps(i,j,1))/ATWAIR/GRAV
        enddo
        loc_frac(i,j,k,1)=(loc_frac(i,j,k,1)*xtot+emis_uemep(k))/(xtot+emis_tot(k)+1.e-20)
      enddo

    enddo ! i
  enddo ! j

  first_call=.false. 

end subroutine uemep_emis

endmodule Emissions_ml
