! <Emissions_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Emissions_mod
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Calls up emission read/set routines
!  This routine interfaces the stand-alone emission-file reading routines
!  with the 3D model.
!_____________________________________________________________________________

use AirEmis_mod, only : airn
use Biogenics_mod,     only: SoilNOx, AnnualNdep
use CheckStop_mod,     only: CheckStop,StopAll
use ChemDims_mod,      only: NSPEC_SHL, NSPEC_TOT,&
                             NEMIS_File  ! No. emission files
use ChemSpecs_mod,     only: NO2, SO2,species,species_adv
use Chemfields_mod,    only: xn_adv
use Config_module,only: &
    KMAX_MID, KMAX_BND, PT ,dt_advec, step_main, &
    KCHEMTOP, &         
    emis_inputlist, &   !TESTC
    EmisDir,      &    ! template for emission path
    DataDir,      &    ! template for path
    EMIS_OUT,      &    ! output emissions in ASCII or not
!    MONTHLY_GRIDEMIS, &  !NML
    INERIS_SNAP2 , &    ! INERIS/TFMM HDD20 method
    MasterProc, USES,  &  !
    SEAFIX_GEA_NEEDED, &  !  see below
    EURO_SOILNOX_DEPSCALE,&! one or the other
    USE_OCEAN_NH3,USE_OCEAN_DMS,FOUND_OCEAN_DMS,&
    NPROC, EmisSplit_OUT,USE_uEMEP,uEMEP,SECTORS_NAME,USE_SECTORS_NAME,&
    SecEmisOutWanted,MaxNSECTORS,&
    AircraftEmis_FLFile,nox_emission_1996_2005File,RoadMapFile,&
    AVG_SMI_2005_2010File,NdepFile,&
    startdate, Emis_sourceFiles, EmisMask
use Country_mod,       only: MAXNLAND,NLAND,Country,IC_NAT,IC_FI,IC_NO,IC_SE
use Country_mod,       only: EU28,EUMACC2,IC_DUMMY
use Debug_module,      only: DEBUG, & !DEBUG => DEBUG_EMISSIONS, & 
                                DEBUG_EMISTIMEFACS
use EmisDef_mod,       only: &
      EMIS_FILE     & ! Names of species ("sox  ",...)
     ,NCMAX         & ! Max. No. countries per grid
     ,ISNAP_DOM     & ! snap index for domestic/resid emis
     ,ISNAP_TRAF    & ! snap index for road-traffic (SNAP7)
     ,NROAD_FILES   & ! No. road dust emis potential files
     ,ROAD_FILE     & ! Names of road dust emission files
     ,NROADDUST     & ! No. road dust components 
     ,QROADDUST_FI  & ! fine road dust emissions (PM2.5) 
     ,QROADDUST_CO  & ! coarse road dust emis
     ,ROADDUST_FINE_FRAC  & ! fine (PM2.5) fraction of road dust emis
     ,ROADDUST_CLIMATE_FILE &! TEMPORARY! file for road dust climate factors 
     ,nGridEmisCodes,GridEmisCodes,GridEmis,cdfemis&
     ,secemis,roaddust_emis_pot,SplitEmisOut,EmisOut&
     ,SecEmisOut,NSecEmisOutWanted,isec2SecOutWanted&
     ,nlandcode,landcode&
     ,road_nlandcode,road_landcode&
     ,gridrcemis,gridrcroadd,gridrcroadd0&
     ,O_NH3, O_DMS&
     ,Emis_4D,N_Emis_4D,Found_Emis_4D & !used for EEMEP 
     ,KEMISTOP&
     ,MAXFEMISLONLAT,N_femis_lonlat,loc_frac &
     ,NSECTORS, N_HFAC, N_TFAC, N_SPLIT     & ! No. emis sectors, height, time and split classes
     ,sec2tfac_map, sec2hfac_map, sec2split_map& !generic mapping of indices
     ,Nneighbors & !used for uemep/loc_frac
     ,NSECTORS_SNAP, SNAP_sec2tfac_map, SNAP_sec2hfac_map, SNAP_sec2split_map&!SNAP specific mapping
     ,NSECTORS_GNFR, GNFR_sec2tfac_map, GNFR_sec2hfac_map, GNFR_sec2split_map&!GNFR specific mapping
     ,NSECTORS_TEST, TEST_sec2tfac_map, TEST_sec2hfac_map, TEST_sec2split_map&!TEST specific mapping
     ,gnfr2snap,snap2gnfr&
     ,foundYearlySectorEmissions, foundMonthlySectorEmissions&
     ,Emis_mask, Emis_mask_allocate, MASK_LIMIT & !old format
     , NEmisMask, EmisMaskValues & !new format
     ,Emis_field, Emis_id, NEmis_id &
     ,NEmisFile_sources, EmisFiles,NEmis_sources, Emis_source&
     , Emis_source_2D, Emis_source_3D,ix3Dmap, NEmis_3Dsources

use EmisGet_mod,       only: &
     EmisSplit &
    ,EmisGetCdf &
    ,EmisGetASCII &
    ,femis                       &  ! Gets scaling factors -> e_fact
    ,e_fact                      &  ! scaling factors
    ,e_fact_lonlat               &  ! scaling factors
    ,EmisHeights                 &  ! Generates vertical distrib
    ,nrcemis, nrcsplit, emisfrac &  ! speciation routines and array
    ,nemis_kprofile, emis_kprofile &! Vertical emissions profile
    ,iqrc2itot,itot2iqrc         &  ! maps from split index to total index
    ,iqrc2iem,iemsplit2itot      &  ! maps from split index to emission index
    ,emis_masscorr               &  ! 1/molwt for most species
    ,emis_nsplit                 &  ! No. species per emis file
    ,RoadDustGet                 &  
    ,roaddust_masscorr           &   ! 1/200. Hard coded at the moment, needs proper setting in EmisGet_mod...   
    ,femis_latmin,femis_latmax,femis_lonmin,femis_lonmax,femis_lonlat_ic&
    ,Emis_init_GetCdf, Emis_GetCdf, make_iland_for_time
use GridValues_mod,    only: GRIDWIDTH_M    & ! size of grid (m)
                           ,xm2            & ! map factor squared
                           ,debug_proc,debug_li,debug_lj & 
                           ,xmd,dA,dB,i_fdom,j_fdom,glon,glat
use Io_Nums_mod,       only: IO_LOG, IO_DMS, IO_EMIS, IO_TMP
use Io_Progs_mod,      only: ios, open_file, datewrite, PrintLog
use MetFields_mod,     only: u_xmj, v_xmi, roa, ps, z_bnd, surface_precip,EtaKz ! ps in Pa, roa in kg/m3
use MetFields_mod,     only: t2_nwp   ! DS_TEST SOILNO - was zero!
use MPI_Groups_mod  
use NetCDF_mod,        only: ReadField_CDF,ReadField_CDF_FL,ReadTimeCDF,IsCDFfractionFormat,&
                             GetCDF_modelgrid,PrintCDF,ReadSectorName,check,&
                             create_country_emission_file, output_country_emissions
use netcdf
use OwnDataTypes_mod,  only: TXTLEN_FILE, TXTLEN_NAME,Emis_id_type, EmisFile_id_type, Emis_sourceFile_id_type
use Par_mod,           only: MAXLIMAX,MAXLJMAX, GIMAX,GJMAX, IRUNBEG,JRUNBEG,&
                            me,limax,ljmax, MSG_READ1,MSG_READ7&
                           ,gi0,gj0,li0,li1,lj0,lj1
use PhysicalConstants_mod,only: GRAV, AVOG, ATWAIR
use PointSource_mod,      only: readstacks !MKPS
use ZchemData_mod,only: rcemis   ! ESX
use SmallUtils_mod,    only: find_index,  key2str, trims
use TimeDate_mod,      only: nydays, nmdays, date, current_date,&! No. days per 
                            tdif_secs,timestamp,make_timestamp,daynumber,day_of_week ! year, date-type, weekday 
use TimeDate_ExtraUtil_mod,only :nctime2date, date2string, to_idate,date_is_reached
use Timefactors_mod,   only: &
     NewDayFactors          & ! subroutines
    ,DegreeDayFactors       & ! degree-days used for SNAP-2
    ,Gridded_SNAP2_Factors, gridfac_HDD & 
    ,fac_min,timefactors   &                  ! subroutine
    ,fac_ehh24x7 ,fac_emm, fac_cemm, fac_edd, timefac & ! time-factors
    ,Read_monthly_emis_grid_fac &
    ,GridTfac &!array with monthly gridded time factors
    ,yearly_normalize !renormalize timefactors after reset
implicit none
private

! subroutines:
public :: Init_masks
public :: Init_Emissions    ! defines emission setup and formats
public :: EmisUpdate        ! update emission
public :: Emissions         ! Main emissions module 
public :: newmonth
public :: EmisSet           ! Sets emission rates every hour/time-step
public :: EmisWriteOut           ! Outputs emissions in ascii

! The main code does not need to know about the following 
private :: expandcclist            !  expands e.g. EU28, EUMACC2
private :: consistency_check       ! Safety-checks


logical, save, private  :: first_dms_read


! and for budgets (not yet used - not changed dimension)
real, public,  save, dimension(NSPEC_SHL+1:NSPEC_TOT) ::  totemadd
integer, private, save :: iemCO  ! index of CO emissions, for debug

real ::TimesInDays(120),mpi_out
integer ::NTime_Read=-1,found
logical :: fileExists            ! to test emission files
integer :: nin, nex         !  include, exclude numbers for emis_inputlist
character(len=*), parameter :: sub='Emissions:'
logical, save :: mappingGNFR2SNAP = .false.

contains
!***********************************************************************
  !new (Nov 2018) emission setup and formats
  subroutine Init_Emissions
    !loop through all sources to set parameters for each source.
    !one source is defined as one two dimensional map giving emissions. 
    !one file can have many sources.
    !each source has parameters. Parameters are constant through the run.
    !There are several ways to set parameters. In increasing priority order:
    !1) default
    !2) read from global attributes in file
    !3) read from variable attributes in file
    !4) file parameter set in namelist (config_emep.nml)
    !5) source parameter set in namelist (config_emep.nml)

    !multiplicative factor on top of each other, from the most general:
    !1) from femis.dat (can be switched off for specific files, but not for specific sources)
    !2) global file attribute (can be overwritten by config)
    !3) individual sources (can be overwritten by config)

    !e_fact (femis) is applied through the source%factor in Emis_init_GetCdf

    !Note on units "as SO2":
    !if the sector is defined, units are defined "as SO2" and SO4 must include the 
    !factor 0.6666 (=mw(SO2)/mw(SO4)) to be right. 
    !If sector is defined as zero, no additional factor is required.

    !note on apply_femis and some other parameters of type "logical":
    !these cannot be defined on file and then possibly overwritten by config, 
    !because it is not possible to distiguish between default and set values from config.
    !apply_femis must be either default (true) or set to false by config
    !include_in_local_fractions must be either default (true) or set to false by config
    
    integer, parameter ::maxnames=100
    character(len=TXTLEN_FILE) :: fname, filename, names_in(maxnames)
    integer :: i, ii, n, nn, ix, nemis_old, isource
    integer :: isec, iland, iem, iqrc, itot, f
    integer :: startsource(size(Emis_sourceFiles)), endsource(size(Emis_sourceFiles))
    type(Emis_id_type):: Emis_id_undefined !to get undefined source values
    type(EmisFile_id_type):: Emisfile_undefined !to get undefined file values
    !Note must distinguish between default and not undefined values, to know when config has defined a value, even if it is the default.
    type(Emis_id_type):: Emis_sources_defaults !set values when not specified otherwise 
    type(EmisFile_id_type):: Emisfile_defaults !set values when not specified otherwise 
    integer :: EmisFilesMap(0:size(Emis_sourceFiles)) !index of EmisFile given index of EmisFile_sources
    integer :: max_levels3D
    character(len=*),parameter :: dtxt='Ini_Em:'
    logical :: found
    logical, save  :: first_call=.true., dbg
    
    if ( first_call ) then
      dbg = DEBUG%EMISSIONS .and. MasterProc
      first_call = .false.
    end if

    !1) define lowest level default values 
    Emis_sources_defaults%units = 'mg/m2/h'
    Emis_sources_defaults%country_ISO = 'N/A'
    Emis_sources_defaults%sector = 0
    Emis_sources_defaults%factor = 1.0
    Emis_sources_defaults%include_in_local_fractions = .true.
    Emis_sources_defaults%apply_femis = .true.
    Emis_sources_defaults%injection_k = KMAX_MID
    Emis_sources_defaults%is3D = .false.
    Emis_sources_defaults%istart = 1 
    Emis_sources_defaults%jstart = 1
    Emis_sources_defaults%kstart = 1 
    Emis_sources_defaults%kend = 1
    Emis_sources_defaults%reversek = .false. 
    Emis_sources_defaults%species_ix = -1 

    Emis_source = Emis_sources_defaults!set all initial values to default

    Emisfile_defaults%periodicity = 'time'
    Emisfile_defaults%projection = 'lon lat'
    Emisfile_defaults%factor = 1.0
    Emisfile_defaults%sectorsName = 'SNAPsectors'
    Emisfile_defaults%units = Emis_sources_defaults%units
    EmisFile_defaults%sector = Emis_sources_defaults%sector
    EmisFile_defaults%country_ISO = Emis_sources_defaults%country_ISO
    
    EmisFiles(:) = EmisFile_defaults

    !2) read from global attributes in file
    !3) read from variable attributes in file
    !Emis_sourceFiles is read from config and then not modified
    !EmisFiles collect all valid data and sources
    !Emis_source() contains the metadata of the emissions surces finally used
    !Emis_source_2D() contains the 2D values for each source used
    NEmis_sources = 0 !total number of valid emis (also 3D) variables across all files
    NEmis_3Dsources = 0 !total number of valid 3D emis variables across all files
    max_levels3D = 1
    NEmisFile_sources = 0 !number of valid emission files
    EmisFilesMap = 0 !index of EmisFile given index of EmisFile_sources
    do n = 1, size(Emis_sourceFiles)
       nemis_old = NEmis_sources
       if(Emis_sourceFiles(n)%filename/='NOTSET')then
          !Read all variables and set parameters as needed
          !find which variables names are defined in config
          names_in='NOTSET'
          i=0
          do nn = 1, size(Emis_sourceFiles(n)%source)
             if(Emis_sourceFiles(n)%source(nn)%varname/='NOTSET')then
                i= i + 1
                names_in(i)=trim(Emis_sourceFiles(n)%source(nn)%varname)
             endif
          enddo

          if(MasterProc)write(*,*)dtxt//"Initializing Emissions from ",&
             trim(Emis_sourceFiles(n)%filename)
          call Emis_init_GetCdf(Emis_sourceFiles(n), EmisFiles(NEmisFile_sources+1), names_in, i)
       endif
       startsource(n) = nemis_old + 1
       endsource(n) = NEmis_sources
       if(NEmis_sources>nemis_old)then
          EmisFilesMap(NEmisFile_sources) = n
          EmisFiles(NEmisFile_sources)%Nsources =  NEmis_sources - nemis_old
          EmisFiles(NEmisFile_sources)%source_start =  nemis_old + 1
          EmisFiles(NEmisFile_sources)%source_end =  NEmis_sources
        endif
    enddo
    allocate(Emis_source_2D(LIMAX,LJMAX,NEmis_sources))

    !4)overwrite parameters with settings from config_emep.nml if they are set
    !first overwrite the global attributes: projection and periodicity
    do i = 1, size(Emis_sourceFiles)
       n = EmisFilesMap(i)
       if(n>0)then
          if(Emis_sourceFiles(n)%periodicity /= Emisfile_undefined%periodicity) EmisFiles(i)%periodicity = Emis_sourceFiles(n)%periodicity
          if(Emis_sourceFiles(n)%projection /= Emisfile_undefined%projection) EmisFiles(i)%projection = Emis_sourceFiles(n)%projection
          if(Emis_sourceFiles(n)%grid_resolution /= Emisfile_undefined%grid_resolution) EmisFiles(i)%grid_resolution = Emis_sourceFiles(n)%grid_resolution
          if(Emis_sourceFiles(n)%projection /= 'native')then
          call CheckStop(EmisFiles(i)%grid_resolution <=1.0E-5,'Grid_resolution must be defined for '//trim(Emis_sourceFiles(n)%filename))
          endif
          if(Emis_sourceFiles(n)%mask_ID /= Emisfile_undefined%mask_ID) EmisFiles(i)%mask_ID = Emis_sourceFiles(n)%mask_ID
          if(Emis_sourceFiles(n)%mask_ID_reverse /= Emisfile_undefined%mask_ID_reverse) EmisFiles(i)%mask_ID_reverse = Emis_sourceFiles(n)%mask_ID_reverse
          if(Emis_sourceFiles(n)%species /= Emisfile_undefined%species) EmisFiles(i)%species = Emis_sourceFiles(n)%species
          if(Emis_sourceFiles(n)%units /= Emisfile_undefined%units) EmisFiles(i)%units = Emis_sourceFiles(n)%units 
          if(Emis_sourceFiles(n)%country_ISO /= Emisfile_undefined%country_ISO) EmisFiles(i)%country_ISO = Emis_sourceFiles(n)%country_ISO
          if(Emis_sourceFiles(n)%sector /= Emisfile_undefined%sector) EmisFiles(i)%sector = Emis_sourceFiles(n)%sector
          if(Emis_sourceFiles(n)%sectorsName /= Emisfile_undefined%sectorsName) EmisFiles(i)%sectorsName = Emis_sourceFiles(n)%sectorsName
       endif
    enddo

    !then overwrite the variable attributes
    do i = 1, NEmisFile_sources !loop over files
       n = EmisFilesMap(i)
       found = .false.
       do ii = EmisFiles(i)%source_start, EmisFiles(i)%source_end !loop over sources found in the netcdf file
          !set source default = file parameter if they are set. Defines default for sources in this file, can be also be redefined for individual sources
          if(Emis_sourceFiles(n)%species /= Emisfile_undefined%species) Emis_source(ii)%species = Emis_sourceFiles(n)%species
          if(Emis_sourceFiles(n)%units /= Emisfile_undefined%units) Emis_source(ii)%units = Emis_sourceFiles(n)%units
          if(Emis_sourceFiles(n)%country_ISO /= Emisfile_undefined%country_ISO) Emis_source(ii)%country_ISO = Emis_sourceFiles(n)%country_ISO
          if(Emis_sourceFiles(n)%sector /= Emisfile_undefined%sector) Emis_source(ii)%sector = Emis_sourceFiles(n)%sector
          Emis_source(ii)%periodicity = EmisFiles(i)%periodicity !NB: periodicity cannot be set individually for variables
          Emis_source(ii)%apply_femis = Emis_sourceFiles(n)%apply_femis !cannot be set in the netcdf file
          Emis_source(ii)%include_in_local_fractions = Emis_sourceFiles(n)%include_in_local_fractions !cannot be set in the netcdf file
          Emis_source(ii)%mask_ID = Emis_sourceFiles(n)%mask_ID !cannot be set in the netcdf file
          Emis_source(ii)%mask_ID_reverse = Emis_sourceFiles(n)%mask_ID_reverse !cannot be set in the netcdf file
          
          isource = Emis_source(ii)%ix_in
          if(dbg) write(*,'(a,4i4,1x,a20)') 'writing config attribute on '//trim(Emis_source(ii)%varname),n, ii, isource, NEmis_sources
          if(isource>0)then
             !source defined in config file
             if(trim(Emis_sourceFiles(n)%source(isource)%varname)/=trim(Emis_source(ii)%varname))write(*,*)isource,'ERROR',trim(Emis_sourceFiles(n)%source(isource)%varname),' ',trim(Emis_source(ii)%varname),ii
             
             found = .true.
             if(Emis_sourceFiles(n)%source(isource)%species /= Emis_id_undefined%species) Emis_source(ii)%species = Emis_sourceFiles(n)%source(isource)%species
             if(Emis_sourceFiles(n)%source(isource)%units /= Emis_id_undefined%units) Emis_source(ii)%units = Emis_sourceFiles(n)%source(isource)%units
             if(Emis_sourceFiles(n)%source(isource)%sector /= Emis_id_undefined%sector) Emis_source(ii)%sector = Emis_sourceFiles(n)%source(isource)%sector
             if(Emis_sourceFiles(n)%source(isource)%factor /= Emis_id_undefined%factor) Emis_source(ii)%factor = Emis_sourceFiles(n)%source(isource)%factor
             if(Emis_sourceFiles(n)%source(isource)%country_ISO /= Emis_id_undefined%country_ISO) Emis_source(ii)%country_ISO = Emis_sourceFiles(n)%source(isource)%country_ISO
             if(.not. Emis_sourceFiles(n)%source(isource)%include_in_local_fractions) Emis_source(ii)%include_in_local_fractions = .false.
             if(.not. Emis_sourceFiles(n)%source(isource)%apply_femis) Emis_source(ii)%apply_femis = Emis_sourceFiles(n)%source(isource)%apply_femis
             
             if(Emis_sourceFiles(n)%source(isource)%mask_ID /= Emis_id_undefined%mask_ID) Emis_source(ii)%mask_ID = Emis_sourceFiles(n)%source(isource)%mask_ID
             if(Emis_sourceFiles(n)%source(isource)%mask_ID_reverse /= Emis_id_undefined%mask_ID_reverse) Emis_source(ii)%mask_ID_reverse = Emis_sourceFiles(n)%source(isource)%mask_ID_reverse

             if(Emis_sourceFiles(n)%source(isource)%is3D .neqv. Emis_id_undefined%is3D) Emis_source(ii)%is3D = Emis_sourceFiles(n)%source(isource)%is3D
             if(Emis_sourceFiles(n)%source(isource)%istart /= Emis_id_undefined%istart) Emis_source(ii)%istart = Emis_sourceFiles(n)%source(isource)%istart
             if(Emis_sourceFiles(n)%source(isource)%jstart /= Emis_id_undefined%jstart) Emis_source(ii)%jstart = Emis_sourceFiles(n)%source(isource)%jstart
             if(Emis_sourceFiles(n)%source(isource)%kstart /= Emis_id_undefined%kstart) Emis_source(ii)%kstart = Emis_sourceFiles(n)%source(isource)%kstart
             if(Emis_sourceFiles(n)%source(isource)%kend /= Emis_id_undefined%kend) Emis_source(ii)%kend = Emis_sourceFiles(n)%source(isource)%kend
             if(Emis_sourceFiles(n)%source(isource)%reversek .neqv. Emis_id_undefined%reversek) Emis_source(ii)%reversek = Emis_sourceFiles(n)%source(isource)%reversek
             if(Emis_sourceFiles(n)%source(isource)%injection_k /= Emis_id_undefined%injection_k) Emis_source(ii)%injection_k = Emis_sourceFiles(n)%source(isource)%injection_k
          endif
          ix = find_index(trim(Emis_source(ii)%country_ISO) ,Country(:)%code, first_only=.true.)
          if(ix<0)then
             if(me==0)write(*,*)dtxt//'WARNING: country '//trim(Emis_source(ii)%country_ISO)//' not defined. '
          else
             Emis_source(ii)%country_ix = ix
             if(dbg)write(*,*)dtxt//'country found '//trim(Emis_source(ii)%country_ISO), ix, NEmis_sources
          endif

          !find if it is defined as an individual species
          ix = find_index(Emis_source(ii)%species, species(:)%name )
          if(ix>0)then
             Emis_source(ii)%species_ix = ix
             if(dbg)write(*,'(a,i4,a)')dtxt//' species found '// &
                trim(Emis_source(ii)%country_ISO), ix, ' '//trim(species(ix)%name)
             if(Emis_source(ii)%include_in_local_fractions .and. USE_uEMEP )then
                if(me==0)write(*,*)"WARNING: local fractions will not include single species "//Emis_source(ii)%species
             endif
          else ! ix<=0 
             if(dbg)write(*,'(a,i4,a)')dtxt//' species not found'// &
              trim(Emis_source(ii)%country_ISO),ix,trim(Emis_source(ii)%species)
          endif
          if(EmisFiles(i)%sectorsName == 'SNAPsectors' .and. SECTORS_NAME == 'GNFR' .and. Emis_source(ii)%sector>0)then
             !map to GNFR sectors
             if(me==0)write(*,*)'mapping SNAP ',Emis_source(ii)%sector,' onto GNFR ',snap2gnfr(Emis_source(ii)%sector,1),' for '//trim(Emis_source(ii)%varname)
             Emis_source(ii)%sector = snap2gnfr(Emis_source(ii)%sector,1) !for now we map onto only one sector
          endif
          if(EmisFiles(i)%sectorsName == 'GNFRsectors' .and. SECTORS_NAME == 'SNAP' .and. Emis_source(ii)%sector>0)then
             !map to SNAP sectors
              if(me==0)write(*,*)'mapping GNFR ',Emis_source(ii)%sector,' onto SNAP ',gnfr2snap(Emis_source(ii)%sector),' for '//trim(Emis_source(ii)%varname)
            Emis_source(ii)%sector = gnfr2snap(Emis_source(ii)%sector)
          endif

          max_levels3D=max(max_levels3D, Emis_source(ii)%kend - Emis_source(ii)%kstart + 1)
          if(MasterProc .and. dbg)write(*,*)dtxt//"Final emission source parameters ",Emis_source(ii)
          
       enddo
!       if(.not. found .and. me==0)write(*,*)dtxt//'WARNING: did not find some of the emission sources defined in config in '&
!            //trim(Emis_sourceFiles(n)%filename)

    enddo

    !include reduction factors
    do n = 1, NEmis_sources      
       if(Emis_source(n)%apply_femis)then
          isec = Emis_source(n)%sector
          if(Emis_source(n)%sector>0 .and. Emis_source(n)%sector<=NSECTORS)then
             iland = Emis_source(n)%country_ix
             if(iland<0)iland=IC_DUMMY
             isec = Emis_source(n)%sector
             iem = find_index(Emis_source(n)%species,EMIS_FILE(:))
             if(iem >0 )then
                !apply femis 
                Emis_source(n)%factor = Emis_source(n)%factor * e_fact(isec,iland,iem)
             else
                !see if the species belongs to any of the splitted species
                iqrc = 0
                do iem = 1,NEMIS_FILE
                   do f = 1,emis_nsplit(iem)
                      iqrc = iqrc + 1
                      itot = iqrc2itot(iqrc)
                      if(trim(species(itot)%name)==trim(Emis_source(n)%species))then
                         if(dbg) write(*,'(a,5i4,a,f12.3)')&
                              trim(Emis_source(n)%species)//' included in '//trim(EMIS_FILE(iem)), n, &
                              isec, iland, Emis_source(n)%country_ix, iem, &
                              trim( species(itot)%name ),  e_fact(isec,iland,iem)
                         Emis_source(n)%factor = Emis_source(n)%factor * e_fact(isec,iland,iem)
                         go to 888
                      endif
                   enddo
                enddo ! iem
888             continue
             endif
          endif
          
       endif

       !define  mask indices
       if(Emis_source(n)%mask_ID /= Emis_id_undefined%mask_ID)then
          ix = find_index(Emis_source(n)%mask_ID,EmisMask(:)%ID)
          if(ix > 0)then
             Emis_source(n)%mask_ix = ix
          else
             Emis_source(n)%mask_ix = Emis_id_undefined%mask_ix !in case somebody tries to set in config
             if(MasterProc)write(*,*)'WARNING mask not defined ',trim(Emis_source(n)%mask_ID)
          endif
       endif
       if(Emis_source(n)%mask_ID_reverse /= Emis_id_undefined%mask_ID_reverse)then
          ix = find_index(Emis_source(n)%mask_ID_reverse,EmisMask(:)%ID)
          if(ix > 0)then
             Emis_source(n)%mask_reverse_ix = ix
          else
             Emis_source(n)%mask_reverse_ix = Emis_id_undefined%mask_reverse_ix !in case somebody tries to set in config
             if(MasterProc)write(*,*)'WARNING mask not defined ',trim(Emis_source(n)%mask_ID_reverse)
          endif
       endif

    enddo ! n = 1, NEmis_sources      

    !find and define the 3D emissions
    ix=0
    do n=1, NEmis_sources
       if(Emis_source(n)%is3D)Nemis_3Dsources = Nemis_3Dsources + 1
       ix3Dmap(n)=ix
    enddo
    if(Nemis_3Dsources>0)then
       if(me==0)write(*,*)'found ',Nemis_3Dsources,' 3D sources'
       allocate(Emis_source_3D(LIMAX,LJMAX,max_levels3D,Nemis_3Dsources))
    endif

  end subroutine Init_Emissions

  subroutine Init_masks()
    !sets and define the masks
    real :: mask_cdf(LIMAX,LJMAX), xsum
    logical :: found
    integer :: iEmisMask = 0
    integer :: i,ii,jj,ic

    iEmisMask = 0
    !1) find number of valid masks defined
    do i = 1, size(EmisMask)
       if(EmisMask(i)%filename /= 'NOTSET' .and. EmisMask(i)%cdfname /= 'NOTSET'&
           .and. EmisMask(i)%ID /= 'NOTSET' ) then
          iEmisMask = iEmisMask+1 !assumes the fields are defined, without checking
       endif
    enddo
    NEmisMask = iEmisMask
    allocate(EmisMaskValues(LIMAX,LJMAX,NEmisMask))
   
    !now set the values for the actual masks
    iEmisMask = 0
    do i = 1, size(EmisMask)
       if(EmisMask(i)%filename /= 'NOTSET' .and. EmisMask(i)%cdfname /= 'NOTSET'&
           .and. EmisMask(i)%ID /= 'NOTSET' ) then

          iEmisMask = iEmisMask+1
          call ReadField_CDF(trim(EmisMask(i)%filename),trim(EmisMask(i)%cdfname),mask_cdf,1,&
               interpol='conservative', needed=.false.,found = found, UnDef=-1.0E10, debug_flag=.false.)

          if(found)then
             !set mask value
             ic = 0
             xsum = 0.0
             do jj = 1, LJMAX
             do ii = 1, LIMAX
                xsum = xsum + mask_cdf(ii,jj)
                if(mask_cdf(ii,jj)>EmisMask(i)%threshold)then
                    EmisMaskValues(ii,jj,iEmisMask) = 0.0 !remove everything 
                    ic = ic + 1
                 else
                    EmisMaskValues(ii,jj,iEmisMask) = 1.0 !keep everything
                 endif
             end do
             end do
             if(MasterProc)write(*,*)'defined mask  '//trim(EmisMask(i)%ID)//' based on '//trim(EmisMask(i)%cdfname)
             if(ic>0)write(*,*)me,' masked ',ic,' cells'
          else
             call StopAll("Mask variable not found: "//trim(EmisMask(i)%cdfname)// &
                  ':'//trim(EmisMask(i)%filename))
             EmisMask(i)%ID = 'NOTSET'!cannot be used anymore
          endif
          
          !to keep some compatibility with old format we also set old format masks
          if(any(emis_inputlist(:)%use_mask))then
             if(.not.allocated(Emis_mask))then
                allocate(Emis_mask(LIMAX,LJMAX))
                Emis_mask = .false.
             endif
             if(MasterProc)write(*,*)'defining mask for old formats ',trim(EmisMask(i)%cdfname)
             do jj = 1, LJMAX
                do ii = 1, LIMAX
                   if(EmisMaskValues(ii,jj,iEmisMask)<0.5)Emis_mask(ii,jj) = .true.
                enddo
             enddo
          endif

       endif
   enddo


  end subroutine Init_masks
!***********************************************************************
  subroutine EmisUpdate
    !Update emission arrays, and read new sets as required
    integer :: n, i, j, ix, is, date_limit(5), iem, ic, icc, iqrc
    integer :: itot,isec,iland
    type(date) :: coming_date
    real :: fac, gridyear, ccsum,emsum(NEMIS_FILE)
    character(len=TXTLEN_NAME) :: fmt
    TYPE(timestamp)   :: ts1,ts2
    logical, save ::first_call=.true.
    real, allocatable, dimension(:,:) :: sumemis ! Sum of emissions per country

    ts1=make_timestamp(current_date)
    coming_date = current_date
    coming_date%seconds = coming_date%seconds + 1800!NB: end_of_validity_date is at end of period, for example 1-1-2018 for December 2017
    gridyear = GRIDWIDTH_M * GRIDWIDTH_M * 3600*24*nydays*1.0E-6!kg/m2/s -> kt/year

    if(first_call)then
       ! sum emissions per countries
       allocate(sumemis(NLAND,NEMIS_FILE))
       sumemis=0.0
    endif

    !loop over all sources and see which one need to be reread from files
    do n = 1, NEmisFile_sources     
       if(date_is_reached(to_idate(EmisFiles(n)%end_of_validity_date,5 )))then
          if(me==0 .and. (step_main<10 .or. DEBUG%EMISSIONS))&
               write(*,*)'Emis: update date is reached ',&
               EmisFiles(n)%end_of_validity_date,EmisFiles(n)%periodicity
          !values are no more valid, fetch new one
          if(first_call)sumemis=0.0
          do is = EmisFiles(n)%source_start,EmisFiles(n)%source_end
             if(Emis_source(is)%is3D)then
                ix = ix3Dmap(is)
                Emis_source_3D(1:,1:,1:,ix)=0.0
                call Emis_GetCdf(EmisFiles(n),Emis_source(is),Emis_source_3D(1,1,1,ix),coming_date)
             else
                if(me==0.and. (step_main<10 .or. DEBUG%EMISSIONS))write(*,*)is,&
                     ' getemis '//trim(Emis_source(is)%units)//' '//trim(Emis_source(is)%varname)//' read as '//trim(Emis_source(is)%species)
                Emis_source_2D(1:,1:,is)=0.0
                call Emis_GetCdf(EmisFiles(n),Emis_source(is),Emis_source_2D(1,1,is),coming_date)
             endif
             !reduction factors
             fac = EmisFiles(n)%factor
             fac = fac* Emis_source(is)%factor     
             
             !unit and factor conversions
             !convert into kg/m2/s
             if(Emis_source(is)%units == 'kg/s' .or. Emis_source(is)%units == 'kg/m2/s')then
                fac = fac
             else if(Emis_source(is)%units == 'g/s' .or. Emis_source(is)%units == 'g/m2/s')then
                fac = fac /(1000.0)
             else if(Emis_source(is)%units == 'mg/s' .or. Emis_source(is)%units == 'mg/m2/s')then
                fac = fac /(1000*1000.0)
             else  
                !depends on periodicity
                if(EmisFiles(n)%periodicity == 'yearly')then
                   if(Emis_source(is)%units == 'kt/m2' .or. Emis_source(is)%units == 'kt/m2/year'&
                        .or. Emis_source(is)%units == 'kt' .or. Emis_source(is)%units == 'kt/year')then
                      fac = fac * 1000 *1000/(3600*24*nydays)
                   else if(Emis_source(is)%units == 'tonnes/m2' .or. Emis_source(is)%units == 'tonnes/m2/year'&
                        .or. Emis_source(is)%units == 'tonnes' .or. Emis_source(is)%units == 'tonnes/year')then
                      fac = fac * 1000/(3600*24*nydays)
                   else if(Emis_source(is)%units == 'kg/m2' .or. Emis_source(is)%units == 'kg/m2/year'&
                        .or. Emis_source(is)%units == 'kg' .or. Emis_source(is)%units == 'kg/year')then
                      fac = fac /(3600*24*nydays)
                   else if(Emis_source(is)%units == 'g/m2' .or. Emis_source(is)%units == 'g/m2/year'&
                        .or. Emis_source(is)%units == 'g' .or. Emis_source(is)%units == 'g/year')then
                      fac = fac /(1000.0*3600*24*nydays)
                   else if(Emis_source(is)%units == 'mg/m2' .or. Emis_source(is)%units == 'mg/m2/year'&
                        .or. Emis_source(is)%units == 'mg' .or. Emis_source(is)%units == 'mg/year')then
                      fac = fac /(1.0e6*3600*24*nydays)
                   else
                      call StopAll("B Unit for emissions not recognized: "//trim(Emis_source(is)%units))                 
                   endif
                else if(EmisFiles(n)%periodicity == 'monthly')then
                   if(Emis_source(is)%units == 'kt/m2' .or. Emis_source(is)%units == 'kt/m2/month'&
                        .or. Emis_source(is)%units == 'kt' .or. Emis_source(is)%units == 'kt/month')then
                      fac = fac *1000*1000/(3600*24*nmdays(coming_date%month))
                   else if(Emis_source(is)%units == 'tonnes/m2' .or. Emis_source(is)%units == 'tonnes/m2/month'&
                        .or. Emis_source(is)%units == 'tonnes' .or. Emis_source(is)%units == 'tonnes/month')then
                      fac = fac *1000/(3600*24*nmdays(coming_date%month))
                   else if(Emis_source(is)%units == 'kg/m2' .or. Emis_source(is)%units == 'kg/m2/month'&
                        .or. Emis_source(is)%units == 'kg' .or. Emis_source(is)%units == 'kg/month')then
                      fac = fac /(3600*24*nmdays(coming_date%month))
                   else if(Emis_source(is)%units == 'g/m2' .or. Emis_source(is)%units == 'g/m2/month'&
                           .or. Emis_source(is)%units == 'g' .or. Emis_source(is)%units == 'g/month')then
                      fac = fac /(1000*3600*24*nmdays(coming_date%month))
                   else if(Emis_source(is)%units == 'mg/m2' .or. Emis_source(is)%units == 'mg/m2/month'&
                        .or. Emis_source(is)%units == 'mg' .or. Emis_source(is)%units == 'mg/month')then
                      fac = fac /(1.0e6*3600*24*nmdays(coming_date%month))
                   else
                      call StopAll("C Unit for emissions not recognized: "//trim(Emis_source(is)%units))                 
                   endif
                else
                   !assume hourly
                   if(Emis_source(is)%units == 'mg/m2' .or. Emis_source(is)%units == 'mg/m2/h')then
                      !convert into kg/m2/s
                      fac = fac /(1.0e6*3600.0)
                   else if(Emis_source(is)%units == 'g/m2' .or. Emis_source(is)%units == 'g/m2/h')then
                      fac = fac /(1000.0*3600.0)
                      if(EmisFiles(n)%periodicity /= 'hourly' .and. Emis_source(is)%units == 'g/m2')then
                         call StopAll("Emis_source unit g/m2 only implemented for hourly, monthly or yearly. Found "//trim(EmisFiles(n)%periodicity))  
                      endif
                   else if(Emis_source(is)%units == 'g/s')then
                      fac = fac /(1000.0)
                   else if(Emis_source(is)%units == 'kg/s')then
                      fac = fac
                   else if(Emis_source(is)%units == 'tonnes/s')then
                      fac = fac * 1000.0
                   else
                      call StopAll("Emis_source unit not implemented. Found "//trim(Emis_source(is)%units)//' '//trim(EmisFiles(n)%periodicity))       
                   endif
                   !Note: easy to implement more unit choices. Just add "if" cases here
                endif
             endif
             
             if(Emis_source(is)%units == 'kt' .or. Emis_source(is)%units == 'kt/s' &
                  .or.Emis_source(is)%units == 'kt/month' .or. Emis_source(is)%units == 'kt/year'  &
                  .or.Emis_source(is)%units == 'tonnes' .or. Emis_source(is)%units == 'tonnes/s' &
                  .or.Emis_source(is)%units == 'tonnes/month' .or. Emis_source(is)%units == 'tonnes/year' &
                  .or. Emis_source(is)%units == 'kg' .or. Emis_source(is)%units == 'kg/s' &
                  .or. Emis_source(is)%units == 'kg/month' .or. Emis_source(is)%units == 'kg/year' &
                  .or. Emis_source(is)%units == 'g' .or. Emis_source(is)%units == 'g/s' &
                  .or. Emis_source(is)%units == 'g/month' .or. Emis_source(is)%units == 'g/year' &
                  .or. Emis_source(is)%units == 'mg' .or. Emis_source(is)%units == 'mg/s' &
                  .or. Emis_source(is)%units == 'mg/month' .or. Emis_source(is)%units == 'mg/year' &
                  .or. Emis_source(is)%units == 'g/h' .or. Emis_source(is)%units == 'mg/h')then   
                   !divide by gridarea
                fac = fac / (GRIDWIDTH_M * GRIDWIDTH_M)
                do j = 1,ljmax
                   do i = 1,limax
                      Emis_source_2D(i,j,is) = Emis_source_2D(i,j,is) * xm2(i,j)
                   enddo
                enddo
             endif
             
             do j = 1,ljmax
                do i = 1,limax
                   Emis_source_2D(i,j,is) = Emis_source_2D(i,j,is) * fac
                enddo
             enddo
             
             !now Emis_source_2D should be in kg/m2/s
          
             !apply masks
             !could easily be better CPU optimized if necessary by putting all factors in same i,j loop
             if(Emis_source(is)%mask_ix>0)then
                do j = 1,ljmax
                   do i = 1,limax
                      Emis_source_2D(i,j,is) = Emis_source_2D(i,j,is) * EmisMaskValues(i,j,Emis_source(is)%mask_ix)
                   enddo
                enddo
             endif
             if(Emis_source(is)%mask_reverse_ix>0)then
                do j = 1,ljmax
                   do i = 1,limax
                      Emis_source_2D(i,j,is) = Emis_source_2D(i,j,is) * (1.0-EmisMaskValues(i,j,Emis_source(is)%mask_ix))!take complementary
                   enddo
                enddo
             endif
             if(first_call)then
                ! sum emissions per countries (in ktonnes?)
                itot = Emis_source(is)%species_ix
                isec = Emis_source(is)%sector
                iland = Emis_source(is)%country_ix
                if(itot>0)then
                   iqrc = itot2iqrc(itot)
                   if(isec>0 .and. iqrc>0)then
                      iem = iqrc2iem(iqrc)
                   else
                      iem=-1
                   endif
                else
                   iem=find_index(Emis_source(is)%species,EMIS_FILE(:))
                endif
                if(iem>0)then
                   do j = 1,ljmax
                      do i = 1,limax                         
                         sumemis(iland,iem) = sumemis(iland,iem) + Emis_source_2D(i,j,is) * gridyear * xmd(i,j) !now in kt/year
                      enddo
                   enddo
                endif
             endif

          enddo
          !update date of valitdity
          if(EmisFiles(n)%periodicity == 'yearly')then
             !assumes only one record to read
             EmisFiles(n)%end_of_validity_date = date(current_date%year,1,1,0,0)
             EmisFiles(n)%end_of_validity_date%year = EmisFiles(n)%end_of_validity_date%year + 1
          else if(EmisFiles(n)%periodicity == 'monthly')then
             !assumes 12 records, one for each month
             EmisFiles(n)%end_of_validity_date = date(current_date%year,1,1,0,0)
             EmisFiles(n)%end_of_validity_date%month = current_date%month + 1     
             if(EmisFiles(n)%end_of_validity_date%month>12)then
                 EmisFiles(n)%end_of_validity_date = date(current_date%year+1,1,1,0,0)                 
             endif
          else
             !the correct times must be written in the file and updated in Emis_GetCdf
          endif          

       endif
       if(first_call)then
           CALL MPI_ALLREDUCE(MPI_IN_PLACE,sumemis,&
               NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,MPI_COMM_CALC,IERROR)
           if(me==0)then
              write(*,*)"Emissions per country for "//trim(EmisFiles(n)%filename)//' (Gg/year) '
              write(*     ,"(a14,a5,3x,30(a12,:))")"EMTAB CC Land ","    ",EMIS_FILE(:)
              fmt="(a5,i4,1x,a9,3x,30(f12.2,:))"
              do ic = 1, NLAND
                 ccsum = sum( sumemis(ic,:) )
                 icc=Country(ic)%icode
                 if ( ccsum > 0.0 )then
                    write(*,     fmt) 'EMTAB', icc, Country(ic)%code, sumemis(ic,:)
                 end if
              end do
           end if

           !total of emissions from all countries and files into emsum
           do iem = 1, NEMIS_FILE
              emsum(iem)= emsum(iem)+sum(sumemis(:,iem))
           end do
        endif
    enddo
    if(first_call)then
       fmt="(a5,i4,1x,a9,3x,30(f12.2,:))"
       if(me==0 .and. NEmisFile_sources>0)write(*     ,fmt)'EMTAB', 999,'TOTAL    ',emsum(:)
       deallocate(sumemis)   
       CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)!so that print out comes out nicely
    endif

    first_call=.false.
 
  end subroutine EmisUpdate

!***********************************************************************
  subroutine Emissions(year)
    ! Initialize emission variables, and read yearly emissions
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
    character(len=300) :: fname ! txt, File name
    logical :: fractionformat

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
    character(len=TXTLEN_NAME) :: varname, fmt,cdf_sector_name, species_name
    integer :: allocerr, i_Emis_4D, sec_ix, ih, id, k, ncFileID, varID, iland
    character(len=TXTLEN_FILE) ::fileName
    integer :: emis_inputlist_NEMIS_FILE!number of files for each emis_inputlist(i)
    real :: buffer(LIMAX,LJMAX)
    logical country_owner_map(NLAND,NPROC)

    if (MasterProc) write(6,*) "Reading emissions for year",  year

    ios = 0
    country_owner_map = .false.

    ! initialize emis_inputlist
    !>=========================================================
    do iemislist = 1, size( emis_inputlist(:)%name )
       fname = emis_inputlist(iemislist)%name
       if(fname=="NOTSET") cycle
       if(MasterProc)&
            write(*,*)"Emission source number ", iemislist,"from ",sub//trim(fname)

       if(emis_inputlist(iemislist)%type == "sectors".or.&
            emis_inputlist(iemislist)%type == "GNFRsectors".or.&
            emis_inputlist(iemislist)%type == "SNAPsectors")then ! Expand groups, e.g. EUMACC2

          call expandcclist( emis_inputlist(iemislist)%incl , n)
          emis_inputlist(iemislist)%Nincl = n
          if(MasterProc .and. n>0) then
             write(*,*) sub//trim(fname)//" INPUTLIST-INCL", n
             write(*,*)'including only countries: ', (trim(emis_inputlist(iemislist)%incl(i))//' ',i=1,n)
          endif
          call expandcclist( emis_inputlist(iemislist)%excl , n)
          emis_inputlist(iemislist)%Nexcl = n
          if(MasterProc .and. n>0) then
             write(*,*)'excluding countries: ', (trim(emis_inputlist(iemislist)%excl(i))//' ',i=1,n)
          endif
       end if
       if(emis_inputlist(iemislist)%pollName(1)/='NOTSET')then
          do iem = 1, NEMIS_FILE
             if(all(emis_inputlist(iemislist)%pollName(:)/=trim(EMIS_FILE(iem))))cycle      
             if(Masterproc)write(*,"(A)")'including pollutant '//trim(EMIS_FILE(iem))//' from '//trim(fname)
          enddo
       else
          !include all pollutants
       endif

       !replace keywords
22     format(5A)
       if(MasterProc)write(*,22)'original emission name ',trim(fname)
       fname = key2str(fname,'EmisDir',EmisDir)
       fname = key2str(fname,'DataDir',DataDir)
       fname = key2str(fname,'YYYY',year)
       emis_inputlist(iemislist)%name=trim(fname)
       if(MasterProc)write(*,22)'filename redefined as: ',&
            trim(emis_inputlist(iemislist)%name)

       cdf_sector_name='NOTSET'
       fname = key2str(fname,'POLL',EMIS_FILE(1))
       call ReadSectorname(fname,cdf_sector_name)
       if(trim(cdf_sector_name)/='NOTSET')then
          SECTORS_NAME=trim(cdf_sector_name)
          if(Masterproc .and. USE_SECTORS_NAME =='NOTSET')&
               write(*,*)"Switching sector categories to ",trim(SECTORS_NAME)
          if(Masterproc .and. USE_SECTORS_NAME =='NOTSET')&
               write(IO_LOG,*)"Switching sector categories to ",trim(SECTORS_NAME)
          if(cdf_sector_name == 'GNFR')emis_inputlist(iemislist)%type = "GNFRsectors"
          if(cdf_sector_name == 'SNAP')emis_inputlist(iemislist)%type = "SNAPsectors"
       end if
       if(IsCDFfractionFormat(trim(fname))) emis_inputlist(iemislist)%format='fractions'
       if(emis_inputlist(iemislist)%set_mask.or.emis_inputlist(iemislist)%use_mask)Emis_mask_allocate = .true.
       if(emis_inputlist(iemislist)%periodicity == 'NOTSET')then
          emis_inputlist(iemislist)%periodicity = 'once' !default
          if(index(emis_inputlist(iemislist)%name,".nc")>1)then
             NTime_Read=-1
             call ReadTimeCDF(trim(fname),TimesInDays,NTime_Read)             
             if(NTime_Read == 12)then
                emis_inputlist(iemislist)%periodicity = "monthly"
             endif
          endif
       endif
       if(emis_inputlist(iemislist)%type == "OceanNH3")then
          if(MasterProc)write(*,*)' using  OceanNH3'    
          USE_OCEAN_NH3=.true.
          O_NH3%index=find_index("NH3",species(:)%name)
          call CheckStop(O_NH3%index<0,'NH3 not found. Needed for OceanNH3 emissions')
          allocate(O_NH3%emis(LIMAX,LJMAX))
          allocate(O_NH3%map(LIMAX,LJMAX))
          O_NH3%emis=0.0
          O_NH3%map=0.0
          O_NH3%sum_month=0.0
          O_NH3%sum_year=0.0
       end if
       if (emis_inputlist(iemislist)%type == "DMS")then
          if(MasterProc)write(*,*)'using DMS'    
          USE_OCEAN_DMS=.true.
          O_DMS%index=find_index("SO2",species(:)%name)
          call CheckStop(O_DMS%index<0,'SO2 not found. Needed for DMS emissions')
          allocate(O_DMS%emis(LIMAX,LJMAX))
          allocate(O_DMS%map(LIMAX,LJMAX))
          O_DMS%emis=0.0
          O_DMS%map=0.0
          O_DMS%sum_month=0.0
          O_DMS%sum_year=0.0
       endif
       if(emis_inputlist(iemislist)%type == "Special_ShipEmis")then
          write(*,*)"ERROR: Special_ShipEmis no more supported. Use new format"
          write(*,*)"emissions from "//trim(fname)//" will not be included"
       endif
    end do ! iemislist

    if(Emis_mask_allocate)then
       if(.not.allocated(Emis_mask))then
          allocate(Emis_mask(LIMAX,LJMAX))
          Emis_mask = .false.
       else
          !the masks has been set by init_mask
       endif
    endif

    ! init_sectors
    if(USE_SECTORS_NAME /='NOTSET')then
       SECTORS_NAME = trim(USE_SECTORS_NAME)
       call CheckStop((SECTORS_NAME /= 'GNFR' .and. SECTORS_NAME /= 'SNAP' .and. SECTORS_NAME /= 'TEST'), &
            'Only SNAP and GNFR (and TEST) can be defined as sector names, not '//trim(SECTORS_NAME))
       if(Masterproc)write(*,*)"Forcing sector categories to ",trim(SECTORS_NAME)          
       if(Masterproc .and. SECTORS_NAME == 'TEST')write(*,*)"WARNING: TEST sectors, requires to define sectors consistently yourself"
    endif

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
    end if
    call CheckStop(NSECTORS>MaxNSECTORS, "redefine larger MaxNSECTORS")

    allocate(Emis_field(LIMAX,LJMAX,10))
    NEmis_id = 0

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
    allocate(road_nlandcode(LIMAX,LJMAX),road_landcode(LIMAX,LJMAX,NCMAX))
    road_nlandcode=0
    road_landcode=0
    allocate(secemis(NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILE))
    secemis=0.0
    allocate(roaddust_emis_pot(LIMAX,LJMAX,NCMAX,NROAD_FILES))
    roaddust_emis_pot=0.0
    allocate(EmisOut(LIMAX,LJMAX,NEMIS_FILE))
    EmisOut=0.0

    allocate(e_fact(NSECTORS,NLAND,NEMIS_FILE))!NLAND defined in Country_Init()
    e_fact=1.0
    allocate(e_fact_lonlat(NSECTORS,MAXFEMISLONLAT,NEMIS_FILE))
    e_fact_lonlat=1.0
    if(.not.allocated(timefac))allocate(timefac(NLAND,N_TFAC,NEMIS_FILE))
    if(.not.allocated(fac_ehh24x7))allocate(fac_ehh24x7(N_TFAC,24,7,NLAND))
    if(.not.allocated(fac_emm))allocate(fac_emm(NLAND,12,N_TFAC,NEMIS_FILE))
    if(.not.allocated(fac_min))allocate(fac_min(NLAND,N_TFAC,NEMIS_FILE))
    if(.not.allocated(fac_edd))allocate(fac_edd(NLAND, 7,N_TFAC,NEMIS_FILE))

    allocate(isec2SecOutWanted(0:NSECTORS))
    isec2SecOutWanted = 0 !should never be used
    isec2SecOutWanted(0) = 0
    NSecEmisOutWanted = 0
    do isec = 1, NSECTORS 
       if(SecEmisOutWanted(isec))then
          NSecEmisOutWanted = NSecEmisOutWanted + 1
          isec2SecOutWanted(isec) = NSecEmisOutWanted
       endif
    enddo
    allocate(SecEmisOut(LIMAX,LJMAX,NEMIS_FILE,0:NSecEmisOutWanted))
    SecEmisOut=0.0

    call femis()              ! emission factors (femis.dat file)
    if(ios/=0) return
    my_first_call = .false.

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
       end if
       !=========================
       call timefactors(year)               ! => fac_emm, fac_edd
       !=========================
    end if
    !=========================
    call EmisSplit()    ! In EmisGet_mod, => emisfrac
    !=========================
    !Must first call EmisSplit, to get nrcemis defined
    if(EmisSplit_OUT)then
       allocate(SplitEmisOut(LIMAX,LJMAX,nrcemis))
       SplitEmisOut=0.0
    end if
    !=========================
    call CheckStop(ios, "ioserror: EmisSplit")

    ! ####################################
    ! Broadcast  monthly and Daily factors (and hourly factors if needed/wanted)
    CALL MPI_BCAST(fac_cemm,8*12,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 
    CALL MPI_BCAST(fac_emm,8*NLAND*12*N_TFAC*NEMIS_FILE,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 
    CALL MPI_BCAST(fac_edd,8*NLAND*7*N_TFAC*NEMIS_FILE,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 
    CALL MPI_BCAST(fac_ehh24x7,8*N_TFAC*24*7*NLAND,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 

    !define fac_min for all processors
    forall(iemis=1:NEMIS_FILE,insec=1:N_TFAC,inland=1:NLAND) &
         fac_min(inland,insec,iemis) = minval(fac_emm(inland,:,insec,iemis))
    if(INERIS_SNAP2) & !  INERIS do not use any base-line for SNAP2
         fac_min(:,ISNAP_DOM,:) = 0.

    ! 4) Read emission files 

    ! allocate for MasterProc (me:=0) only:
    err1 = 0
    if(MasterProc) then
       if(USES%ROADDUST)then
          allocate(road_globnland(GIMAX,GJMAX),stat=err7)
          allocate(road_globland(GIMAX,GJMAX,NCMAX),stat=err8)
          allocate(globroad_dust_pot(GIMAX,GJMAX,NCMAX),stat=err9)
          allocate(RoadDustEmis_climate_factor(GIMAX,GJMAX),stat=err1)

          call CheckStop(err7, "Allocation error 7 - globroadland")
          call CheckStop(err8, "Allocation error 8 - globroadland")
          call CheckStop(err9, "Allocation error 9 - globroad_dust_pot")
          call CheckStop(err1, "Allocation error 1 - RoadDustEmis_climate_factor")
       end if ! road dust


       if(USES%ROADDUST)then
          road_globnland(:,:)=0
          road_globland(:,:,:)=0
          globroad_dust_pot(:,:,:)=0.
          RoadDustEmis_climate_factor(:,:)=1.0 ! default, no scaling
       end if ! road dust
    else
       ! needed for DEBUG=yes compilation options
       if(USES%ROADDUST)then
          allocate(road_globnland(1,1),road_globland(1,1,1),&
               globroad_dust_pot(1,1,1),stat=err9)
          call CheckStop(err9, "Allocation error 9 - dummy roadglob")
       end if ! road dust
    end if

    emsum=0.0
    Found_Emis_4D=0
    do iemislist = 1, size( emis_inputlist(:)%name )
       fractionformat = ( emis_inputlist(iemislist)%format=='fractions' )
       fname=emis_inputlist(iemislist)%name
       if ( fname == "NOTSET" ) cycle
       if (emis_inputlist(iemislist)%periodicity /= 'once' ) cycle
38     FORMAT(A,I4,A)
       if(MasterProc)write(*,38)sub//' reading emis_inputlist ',iemislist,trim(fname)

       sumemis=0.0
       sumemis_local(:,:)=0.0

       nin = emis_inputlist(iemislist)%Nincl
       nex = emis_inputlist(iemislist)%Nexcl

       if((emis_inputlist(iemislist)%type == "sectors".or.&
            emis_inputlist(iemislist)%type == "GNFRsectors".or.&
            emis_inputlist(iemislist)%type == "SNAPsectors") .and. index(emis_inputlist(iemislist)%name,".nc")>1)then 

          foundYearlySectorEmissions = .true.
          mappingGNFR2SNAP = .false.
          emis_inputlist_NEMIS_FILE = 1!all pollutants are in same file
          if(index(emis_inputlist(iemislist)%name,"POLL")>0)emis_inputlist_NEMIS_FILE = NEMIS_FILE !one file per pollutant
          do iem = 1, emis_inputlist_NEMIS_FILE
             if(emis_inputlist(iemislist)%pollName(1)/='NOTSET')then
                if(all(emis_inputlist(iemislist)%pollName(:)/=trim(EMIS_FILE(iem))))cycle      
                if(Masterproc)write(*,"(A)")'using PollNames restrictions '
             end if
             fname = key2str(trim(emis_inputlist(iemislist)%name),'POLL',EMIS_FILE(iem))

             call EmisGetCdf(iem, fname, sumemis_local, &
                  GridEmis, GridEmisCodes, nGridEmisCodes, 1,&
                  emis_inputlist(iemislist)%incl, nin, emis_inputlist(iemislist)%excl, nex, &
                  emis_inputlist(iemislist)%use_lonlat_femis,&
                  emis_inputlist(iemislist)%set_mask,emis_inputlist(iemislist)%use_mask,&
                  emis_inputlist(iemislist)%pollName,&
                  fractionformat,emis_inputlist(iemislist)%type)

          end do!NEMIS_FILE

          !add together totals from each processor (only me=0 get results)
          sumemis=0.0
          CALL MPI_REDUCE(sumemis_local,sumemis,&
               NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)        

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
                   write(*,*)'Emis_4D: will write to ',&
                        n,trim(emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D))
                else                   
                   write(*,*)'Emis_4D: WARNING did not find ',&
                        trim(emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D)),&
                        ' among the emep species'
                   write(*,*)'Emis_4D: WARNING ',&
                        trim(emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)),&
                        ' is not used'
                end if
             end if
          end do
          !   else if(IsCDFSnapFormat(trim(emis_inputlist(iemislist)%name)))then !This Does not work because of "POLL"
       elseif(index(emis_inputlist(iemislist)%name,"grid")>0)then
          !ASCII format
          do iem = 1, NEMIS_FILE
             if(emis_inputlist(iemislist)%pollName(1)/='NOTSET')then
                if(all(emis_inputlist(iemislist)%pollName(:)/=trim(EMIS_FILE(iem))))cycle      
                if(Masterproc)write(*,*)'reading '//trim(EMIS_FILE(iem))//' from '//trim(fname)
             endif
             fname=key2str(emis_inputlist(iemislist)%name,'POLL',EMIS_FILE(iem)) ! e.g. POLL -> sox
             if(MasterProc)write(*,fmt='(A)')'Reading ASCII format '//trim(fname)
             call EmisGetASCII(iem, fname, trim(EMIS_FILE(iem)), sumemis_local, &
                  emis_inputlist(iemislist)%incl, nin, emis_inputlist(iemislist)%excl, nex, &
                  emis_inputlist(iemislist)%use_lonlat_femis)
          end do

          !add together totals from each processor (only me=0 get results)
          sumemis=0.0
          CALL MPI_REDUCE(sumemis_local,sumemis,&
               NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)        
       else
          if(MasterProc)write(*,*)'WARNING: did not recognize format of '//trim(emis_inputlist(iemislist)%name)
          call StopAll("Emissions file format not recognized ")
       end if


       if(MasterProc.and. emis_inputlist(iemislist)%periodicity == "once") then
          ! Added EMTAB to make parsing easy. These data are important!
          call PrintLog("#EMTAB Total emissions by countries for "//trim(emis_inputlist(iemislist)%name)//" (Gg)")
          write(*     ,"(a14,a5,3x,30(a12,:))")"EMTAB CC Land ","    ",EMIS_FILE(:)
          write(IO_LOG,"(a14,a5,3x,30(a12,:))")"EMTAB CC Land ","    ",EMIS_FILE(:)                
          sumEU(:) = 0.0
          fmt="(a5,i4,1x,a9,3x,30(f12.2,:))"
          do ic = 1, NLAND
             ccsum = sum( sumemis(ic,:) )
             icc=Country(ic)%icode
             if ( ccsum > 0.0 )then
                write(*,     fmt) 'EMTAB', icc, Country(ic)%code, sumemis(ic,:)
                write(IO_LOG,fmt) 'EMTAB', icc, Country(ic)%code, sumemis(ic,:)
             end if
             if(find_index(Country(ic)%code,EU28(:))>0) sumEU = sumEU + sumemis(ic,:)
          end do
          if ( sum(sumEU(:))>0.001) then
             write(*     ,fmt) 'EMTAB', 998, "EU", sumEU(:)
             write(IO_LOG,fmt) 'EMTAB', 998, "EU", sumEU(:)
          end if
       end if

       !total of emissions from all countries and files into emsum
       do iem = 1, NEMIS_FILE
          emsum(iem)= emsum(iem)+sum(sumemis(:,iem))
       end do

    end do

    if(foundMonthlySectorEmissions .and. .not. foundYearlySectorEmissions)then
       !do not use any additional monthly time factors!
       if(me==0)write(*,*)'WARNING: monthly emissions used. Resetting monthly emission factors!'
       fac_emm=1.0
       call yearly_normalize(year)
    endif


    if(MasterProc)then
       write(*     ,"(a14,3x,30(f12.2,:))")'EMTAB 999 TOTAL',emsum(:)
       write(IO_LOG,"(a14,3x,30(f12.2,:))")' EMTAB 999 TOTAL ',emsum(:)
    end if

    if(USES%ROADDUST) then
       !Use grid-independent Netcdf input files
       call CheckStop(NROAD_FILES>2, "TOO MANY ROADFILES")
       do iem = 1, NROAD_FILES
          !Read data from NetCDF file
          select case(iem)
          case(1);varname='HighwayRoadDustPM10_Jun-Feb'
          case(2);varname='nonHighwayRoadDustPM10_Jun-Feb'
          end select
          roaddust_emis_pot(:,:,:,iem)=0.0
          call ReadField_CDF(RoadMapFile,varname,roaddust_emis_pot(1,1,1,iem),&
               nstart=1,interpol='mass_conservative',fractions_out=fractions,&
               CC_out=road_landcode,Ncc_out=road_nlandcode,needed=.true.,&
               debug_flag=.false.,Undef=0.0)
          if(.not.SMI_defined)then
             varname='SMI1'
             call ReadField_CDF(AVG_SMI_2005_2010File,varname,SMI,nstart=1,&
                  interpol='conservative',needed=.true.,debug_flag=.false.)
             SMI_defined=.true.
          end if

          do i=1,LIMAX
             do j=1,LJMAX
                !Peter: Rough estimate to get something varying between 3.325 (SMI<0.5) and 1.0 (SMI>1)
                SMI_roadfactor=3.325-(min(1.0,max(0.5,SMI(i,j)))-0.5)*2*(3.325-1.0)
                !if(DEBUG%ROADDUST)&
                !  WRITE(*,*)"i,j,RDECF:",i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,SMI_roadfactor
                do iic=road_nlandcode(i,j),1,-1
                   roaddust_emis_pot(i,j,iic,iem)=roaddust_emis_pot(i,j,1,iem) &
                        *fractions(i,j,iic)*SMI_roadfactor
                end do
             end do
          end do
          sumroaddust_local(:,iem)=0.0
          do i=1,LIMAX
             do j=1,LJMAX
                do iic=1,road_nlandcode(i,j)
                   if(road_landcode(i,j,iic)<=NLAND) &
                        ic=find_index(road_landcode(i,j,iic),Country(:)%icode)
                   if(Country(ic)%icode/=road_landcode(i,j,iic))then
                      write(*,*)"COUNTRY ROAD CODE ERROR: ",road_landcode(i,j,iic),ic,Country(ic)%icode
                      call StopAll("COUNTRY CODE ERROR ")
                   end if
                   if(ic>NLAND)then
                      write(*,*)"COUNTRY CODE NOT RECOGNIZED OR UNDEFINED: ",road_landcode(i,j,iic)
                      call StopAll("COUNTRY CODE NOT RECOGNIZED ")
                   end if
                   sumroaddust_local(ic,iem)=sumroaddust_local(ic,iem)&
                        +0.001*roaddust_emis_pot(i,j,iic,iem)
                end do
             end do
          end do
       end do ! iem = 1, NROAD_FILES-loop
       sumroaddust=0.0
       CALL MPI_REDUCE(sumroaddust_local,sumroaddust,NLAND*NROAD_FILES,MPI_REAL8,&
            MPI_SUM,0,MPI_COMM_CALC,IERROR) 

    end if !USES%ROADDUST

    if(MasterProc) then
       if(USES%ROADDUST)THEN
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
             end if
          end do
       end if ! ROAD DUST
    end if

    ! Create emislist-type (ASCII) files for both sec emissions and Cdf
    ! Useful for export to other codes, including production of
    ! new emislist for current NWP grid.
    do iem = 1, NEMIS_FILE
       if(EMIS_OUT) &
            call EmisWriteOut("Sector",iem,nGridEmisCodes,GridEmisCodes,GridEmis(:,:,:,:,iem))
    end do

    !**  Conversions:
    ! The emission-data file are so far in units of 
    ! tonnes per grid-square. The conversion factor from tonnes per gridcell
    ! annual emission values to surface flux (kg/m2/s) is found by division
    ! with (nydays*24*60*60)s and (h*h)m2 and multiply by 1.e+3.
    ! The conversion factor then equals 1.27e-14
    tonne_to_kgm2s  = 1.0e3 / (nydays * 24.0 * 3600.0 * GRIDWIDTH_M * GRIDWIDTH_M)
    if(DEBUG%EMISSIONS .and.MasterProc) then
       write(*,*) "CONV:me, nydays, gridwidth = ",me,nydays,GRIDWIDTH_M
       write(*,*) "No. days in Emissions: ", nydays
       write(*,*) "tonne_to_kgm2s in Emissions: ", tonne_to_kgm2s
       write(*,*) "Emissions sums:"
       do iem = 1, NEMIS_FILE
          write(*,"(a15,f12.2)") EMIS_FILE(iem),emsum(iem)
       end do
    end if

    iemCO=find_index("co",EMIS_FILE(:)) ! save this index

    if(DEBUG%EMISSIONS .and.debug_proc.and.iemCO>0) &
         write(*,"(a,2es10.3)") "SnapPre:" // trim(EMIS_FILE(iemCO)), &
         sum(secemis(:,debug_li,debug_lj,:,iemCO))

    forall (ic=1:NCMAX, j=1:ljmax, i=1:limax, isec=1:NSECTORS,iem=1:NEMIS_FILE)
       GridEmis(isec,i,j,ic,iem) = GridEmis(isec,i,j,ic,iem) * tonne_to_kgm2s * xm2(i,j)
    endforall

    !nGridEmisCodes, GridEmisCodes, GridEmis are yearly and fixed
    !nlandcode,landcode,secemis are updated
    nlandcode=nGridEmisCodes
    landcode=GridEmisCodes
    secemis=GridEmis

    if(DEBUG%EMISSIONS .and.debug_proc.and.iemCO>0) &
         write(*,"(a,2es10.3)") "SnapPos:" // trim(EMIS_FILE(iemCO)), &
         sum(secemis   (:,debug_li,debug_lj,:,iemCO))

    if(USES%ROADDUST)THEN
       forall (ic=1:NCMAX, j=1:ljmax, i=1:limax, iem=1:NROAD_FILES)
          roaddust_emis_pot(i,j,ic,iem) = &
               roaddust_emis_pot(i,j,ic,iem) * tonne_to_kgm2s * xm2(i,j)
       endforall
    end if !road dust

    err1 = 0
    if(MasterProc) then
       if(USES%ROADDUST)THEN
          deallocate(road_globnland   ,stat=err7)
          deallocate(road_globland    ,stat=err8)
          deallocate(globroad_dust_pot,stat=err9)
          call CheckStop(err7, "De-Allocation error 7 - roadglob")
          call CheckStop(err8, "De-Allocation error 8 - roadglob")
          call CheckStop(err9, "De-Allocation error 9 - roadglob")
       end if
    else
       ! needed for DEBUG=yes compilation options
       if(USES%ROADDUST)THEN
          deallocate(road_globnland,road_globland,globroad_dust_pot,stat=err9)
          call CheckStop(err9, "De-Allocation error 9 - dummy roadglob")
       end if
    end if

    ! now we have nrecmis and can allocate for gridrcemis:
    ! print *, "ALLOCATING GRIDRC", me, NRCEMIS
    allocate(gridrcemis(NRCEMIS,KEMISTOP:KMAX_MID,LIMAX,LJMAX),stat=err1)
    call CheckStop(err1, "Allocation error 1 - gridrcemis") 
    if(USES%ROADDUST)THEN
       allocate(gridrcroadd(NROADDUST,LIMAX,LJMAX),stat=err3)
       allocate(gridrcroadd0(NROADDUST,LIMAX,LJMAX),stat=err4)
       call CheckStop(err3, "Allocation error 3 - gridrcroadd")
       call CheckStop(err4, "Allocation error 4 - gridrcroadd0")
    end if

    !output emissions
    fileName = 'EMIS_OUT.nc'
!Not ready
!    call output_country_emissions(filename,GridEmis,GridEmisCodes,nGridEmisCodes,NSECTORS,NCMAX,EMIS_FILE,NEMIS_FILE)
!    CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
!    stop
  end subroutine Emissions
!----------------------------------------------------------------------!
!>
!! expandcclist converts e.g. EU28 to individual countries
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
    end select
  end do CCLIST ! i
  xlist(1:n) = nlist(1:n) ! overwrites original
end subroutine expandcclist
!----------------------------------------------------------------------!
subroutine consistency_check()
!----------------------------------------------------------------------!
!    checks that all the values given so far are consistent
!----------------------------------------------------------------------!
  character(len=30) :: errormsg 
  errormsg = "ok"
  if(size(EMIS_FILE)/=NEMIS_FILE) errormsg = " size EMISNAME wrong "
  call CheckStop(errormsg,"Failed consistency check")
end subroutine consistency_check
!***********************************************************************
subroutine EmisSet(indate)   !  emission re-set every time-step/hour
!----------------------------------------------------------------------!
! DESCRIPTION:
!   Calculates the emission time variations.
!   Time factor recalcualted once per hour to allow for day/night variation (and voc 
!   speciation) based on local time for each sector and country.
!   gridrcemis calculated every time-step to allow for ps changes.
!   inputs from Emissions in EMISSIONS_ML:
!   country and sector-specific array : 
!          secemis (NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILES)
!  
!   Units:
!   secemis has units of kg/m2/s. The values are read from the input file, and multiplied by
!   tonne_to_kgm2s*xm2 (for yearly emissions).  
!   When used, the emissions are converted into molecules/cm3/s every timestep in the routine EmisSet.
!   This gives gridrcemis(NRCEMIS,KEMISTOP:KMAX_MID,LIMAX,LJMAX) 
!   which has no sector and country info, and the NEMIS_FILES are speciated into NRCEMIS species.
!   For conversion from kg to molecules, the molecular masses are assumed to be:
!   SOx as SO2, NOx as NO2 
!   VOC and PM as speciated mw mass
!   (CO and NH3 are not split and have their normal mw mass).
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
  integer :: i, j, k, n, f    ! cooridnates, loop variables
  integer :: icc, ncc         ! No. of countries in grid.
  integer :: ficc,fncc        ! No. of countries with
  integer :: iqrc             ! emis indices 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILE
  integer :: itot             ! index in xn()

  logical                         ::  hourchange   !             
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions

  real ::  ehlpcom,ehlpcom0
  real ::  tfac, dtgrid    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s               ! source term (emis) before splitting
  integer :: iland, iland_timefac, iland_timefac_hour  ! country codes, and codes for timefac 

  integer, save :: oldday = -1, oldhour = -1
  integer, save :: wday , wday_loc ! wday = day of the week 1-7
  real ::  oldtfac, lon
  logical :: debug_tfac

  ! If timezone=-100, calculate time based on longitude rather than timezone
  integer :: hour_iland,nstart
  integer :: i_Emis_4D, iem_treated
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
    end if

    if(Found_Emis_4D>0)then
      if(.not.allocated(Emis_4D))allocate(Emis_4D(LIMAX,LJMAX,KMAX_MID,N_Emis_4D))
      Emis_4D = 0.0 !default value, must be set to zero when no values are found
      NTime_Read=-1 !read all times    
      call ReadTimeCDF(emis_inputlist(Found_Emis_4D)%Name,TimesInDays,NTime_Read)
      call CheckStop(NTime_Read>size(TimesInDays), "Emissions_mod: increase size of TimesInDays ")
      !if(MasterProc)write(*,*)('found date ',i,TimesInDays(i),i=1,NTime_Read)
      !write(*,*)'compare  ',ts1,ts2
      ts1=make_timestamp(indate)
      do i=1,NTime_Read
        call nctime2date(ts2,TimesInDays(i))   
        if(nint(tdif_secs(ts1,ts2))==0)then
          if(MasterProc)write(*,*)'Emis_4D: found matching date ',i,TimesInDays(i)
          nstart=i
          exit
        end if
      end do
      if(i>NTime_Read )then
        if(MasterProc)then
          write(*,*)'Emis_4D: WARNING DID NOT FIND ANY MATCHING DATE '
          write(*,*)'Emis_4D: first date found ',TimesInDays(1)
          write(*,*)'Emis_4D: last date found ',TimesInDays(NTime_Read)
          write(*,*)'Emis_4D: difference to last date ',tdif_secs(ts1,ts2)/3600,' hours'
        end if
      else
        do i_Emis_4D=1,N_Emis_4D
          if(emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)=='NOTSET')exit
          varname=emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)
          !if(MasterProc)write(*,*)'Fetching ',trim(varname)
          call GetCDF_modelgrid(varname,emis_inputlist(Found_Emis_4D)%Name,&
            Emis_4D(1,1,1,i_Emis_4D),1,kmax_mid,nstart,1,reverse_k=.true.)
        end do
      end if
    end if
  end if

  if(DEBUG_EMISTIMEFACS.and.MasterProc) &
    write(*,"(a,2f8.3)") " EmisSet  traffic 24x7", &
      fac_ehh24x7(ISNAP_TRAF,1,4,1),fac_ehh24x7(ISNAP_TRAF,13,4,1)
  !..........................................

  if(hourchange) then 
    totemadd(:)  = 0.
    gridrcemis(:,:,:,:) = 0.0 
    SecEmisOut(:,:,:,:) = 0.0
    if(USES%ROADDUST)gridrcroadd0(:,:,:) = 0.0
    !..........................................
    ! Process each grid:
    do j = 1,ljmax
      do i = 1,limax
        ncc = nlandcode(i,j)            ! No. of countries in grid
        debug_tfac=(DEBUG_EMISTIMEFACS.and.debug_proc.and.i==DEBUG_li.and.j==DEBUG_lj)
        ! find the approximate local time:
        lon = modulo(360+nint(glon(i,j)),360)
        if(lon>180.0)lon=lon-360.0
        
        !*************************************************
        ! First loop over sector emissions
        !*************************************************
        tmpemis(:)=0.
        do icc = 1, ncc
          !          iland = landcode(i,j,icc)     ! 1=Albania, etc.
          iland=find_index(landcode(i,j,icc),Country(:)%icode) !array index
          call make_iland_for_time(debug_tfac, indate, i, j, iland, wday, iland_timefac,hour_iland,wday_loc,iland_timefac_hour)

          !  As each emission sector has a different diurnal profile
          !  and possibly speciation, we loop over each sector, adding
          !  the found emission rates to gridrcemis as we go.
          !  ==================================================
          do isec = 1, NSECTORS       ! Loop over snap codes
            ! Calculate emission rates from secemis, time-factors, 
            ! and if appropriate any speciation fraction (NEMIS_FRAC)
            iqrc = 0   ! index over emisfrac
            do iem = 1, NEMIS_FILE 

              tfac = timefac(iland_timefac,sec2tfac_map(isec),iem) &
                   * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)

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
                     * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)

                if(debug_tfac .and. indate%hour==12 .and. iem==1) &
                  write(*,"(a,3i3,2i4,7f8.3)") "SNAPHDD tfac ",  &
                    isec, sec2tfac_map(isec),iland, daynumber, indate%hour, &
                    timefac(iland_timefac,sec2tfac_map(isec),iem), t2_nwp(i,j,2)-273.15, &
                    fac_min(iland,sec2tfac_map(isec),iem),  gridfac_HDD(i,j), tfac
              end if ! =============== HDD 

              s = tfac * secemis(isec,i,j,icc,iem)
              ! prelim emis sum kg/m2/s
              SecEmisOut(i,j,iem,0) = SecEmisOut(i,j,iem,0) + s !sum of all sectors
              if(SecEmisOutWanted(isec))SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) = &
                   SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) + s

              do f = 1,emis_nsplit(iem)
                iqrc = iqrc + 1
                itot = iqrc2itot(iqrc)
                tmpemis(iqrc) = s * emisfrac(iqrc,sec2split_map(isec),iland)
                ! Add up emissions in ktonne 
                totemadd(itot) = totemadd(itot) &
                     + tmpemis(iqrc) * dtgrid * xmd(i,j)
              end do ! f
            end do ! iem

            !  Assign to height levels 1-KEMISTOP
            do k=KEMISTOP,KMAX_MID
              do iqrc =1, nrcemis
                gridrcemis(iqrc,k,i,j) = gridrcemis(iqrc,k,i,j)   &
                     + tmpemis(iqrc)*ehlpcom0    &
                     *emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) &
                     *emis_masscorr(iqrc)
                !if( debug_tfac.and. iqrc==1 ) then 
                !  write(*,"(a,2i3,2f8.3)") "KPROF ", &
                !    isec, KMAX_BND-k, &
                !    VERTFAC(KMAX_BND-k,sec2hfac_map(isec)),  &
                !    emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))
                !end if
              end do ! iqrc
            end do   ! k
          end do  ! isec
          !      ==================================================
        end do ! icc  

        if(USES%ROADDUST)then
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

            ! if( DEBUG%ROADDUST .and. debug_proc .and. i==DEBUG_li .and. j==DEBUG_lj )THEN
            !    write(*,*)"DEBUG ROADDUST! Dry! ncc=", road_nlandcode(i,j)
            ! end if

            ncc = road_nlandcode(i,j) ! number of countries in grid point
            do icc = 1, ncc    
              !iland = road_landcode(i,j,icc)
              iland = find_index(road_landcode(i,j,icc),Country(:)%icode)
              iland_timefac_hour = find_index(Country(iland)%timefac_index_hourly,Country(:)%icode)
              if(Country(iland)%timezone==-100)then
                hour_iland=mod(nint(indate%hour+24*(lon/360.0)),24) + 1
              else
                hour_iland=indate%hour + Country(iland)%timezone + 1! add 1 to get 1..24 
              end if

              wday_loc = wday ! DS added here also, for fac_ehh24x7
              if( hour_iland > 24 ) then
                hour_iland = 1
                if(wday_loc==0)wday_loc=7 ! Sunday -> 7
                if(wday_loc>7 )wday_loc=1 
              end if

              if(ANY(iland==(/IC_FI,IC_NO,IC_SE/)).and. & ! Nordic countries
                 ANY(indate%month==(/3,4,5/)))then        ! spring road dust
                tfac = fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc,iland_timefac_hour)*2.0 ! Doubling in Mar-May (as in TNO model)
              else
                tfac = fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc,iland_timefac_hour)
              end if

              do iem = 1, NROAD_FILES
                s = tfac * roaddust_emis_pot(i,j,icc,iem)
                if(DEBUG%ROADDUST.and.debug_proc.and.i==DEBUG_li.and.j==DEBUG_lj)&
                  write(*,*)"DEBUG ROADDUST! iem,tfac,icc,roaddust_emis_pot,s", &
                    iem,tfac,icc,roaddust_emis_pot(i,j,icc,iem),s

                gridrcroadd0(QROADDUST_FI,i,j)=gridrcroadd0(QROADDUST_FI,i,j) &
                     +ROADDUST_FINE_FRAC*s
                gridrcroadd0(QROADDUST_CO,i,j)=gridrcroadd0(QROADDUST_CO,i,j) &
                     +(1.-ROADDUST_FINE_FRAC)*s

                if(all([DEBUG%ROADDUST,debug_proc,i==debug_li,j==debug_lj]))then
                  write(*,*)"gridrcroadfine"  ,gridrcroadd0(QROADDUST_FI,i,j)
                  write(*,*)"gridrcroadcoarse",gridrcroadd0(QROADDUST_CO,i,j)
                end if
              end do ! nroad files
            end do   ! icc
            ! should pick the correct emissions (spring or rest of year)
            ! and add the emissions from HIGHWAYplus and NONHIGHWAYS,
            ! using correct fine and coarse fractions.
          else ! precipitation
            gridrcroadd0(:,i,j)=0.
          end if NO_PRECIP
        end if ! ROADDUST
      end do   ! i
    end do     ! j

    do n = 1,NEmis_id
       itot = find_index(Emis_id(n)%species,species(:)%name)
       iqrc = itot2iqrc(itot)
       k=KMAX_MID
       do j = 1,ljmax
          do i = 1,limax
             gridrcemis(iqrc,k,i,j) = gridrcemis(iqrc,k,i,j)   &
                  + Emis_field(i,j,n)*ehlpcom0    &
                  *emis_masscorr(iqrc)
          end do   ! i
       end do     ! j
    enddo

!Add emissions from new format
    do n = 1, NEmis_sources 
       itot = Emis_source(n)%species_ix
       isec = Emis_source(n)%sector
       iland = Emis_source(n)%country_ix
       if(itot>0)then
          !the species is directly defined (no splits)
          iqrc = itot2iqrc(itot)
          if(isec>0)then
             call CheckStop(iqrc<=0,"emitted sector species must belong to one of the splitted species")
             iem = iqrc2iem(iqrc)
             do j = 1,ljmax
                do i = 1,limax
                   
                   if(Emis_source(n)%periodicity == 'yearly' .or. Emis_source(n)%periodicity == 'monthly')then
                      !we need to apply hourly factors
                      call make_iland_for_time(debug_tfac, indate, i, j, iland, wday, iland_timefac,hour_iland,wday_loc,iland_timefac_hour)                      
                      tfac = fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)
                      if(Emis_source(n)%periodicity == 'yearly')then
                         !apply monthly factor on top of hourly factors
                         tfac = tfac * timefac(iland_timefac,sec2tfac_map(isec),iem)                         
                      endif                      
                   else
                      !not monthly or yearly emissions, timefactors must be included in emission values
                      tfac = 1.0
                   endif
                   
                   s = Emis_source_2D(i,j,n) * tfac
                   SecEmisOut(i,j,iem,0) = SecEmisOut(i,j,iem,0) + s !sum of all sectors
                   if(SecEmisOutWanted(isec))SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) = &
                        SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) + s
                      ! Add up emissions in ktonne 
                   totemadd(itot) = totemadd(itot) &
                        + s * dtgrid * xmd(i,j)
                   !  Assign to height levels 1-KEMISTOP
                   do k=KEMISTOP,KMAX_MID
                      gridrcemis(iqrc,k,i,j) = gridrcemis(iqrc,k,i,j)   &
                           + s&
                           *ehlpcom0    &
                           *emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) &
                           *emis_masscorr(iqrc) !NB: assumes mass defined as a split nox, pm25... !
                   end do   ! k
                end do ! i
             enddo
          else
             ! we do not include the emissions as a sector emission
             !directly included in setup_rcemis
          endif
       else
          !the species is defined as a sector emission
          iem=find_index(Emis_source(n)%species,EMIS_FILE(:))
          call CheckStop(iem<0, "did not recognize species "//trim(Emis_source(n)%species))
          call CheckStop(Emis_source(n)%sector<=0," sector must be defined for "//trim(Emis_source(n)%varname))
          do f = 1,emis_nsplit(iem)
             itot = iemsplit2itot(f,iem)
             call CheckStop(itot<0, "did not recognize split "//trim(Emis_source(n)%species))
             iqrc = itot2iqrc(itot)
             do j = 1,ljmax
                do i = 1,limax
                   if(Emis_source(n)%periodicity == 'yearly' .or. Emis_source(n)%periodicity == 'monthly')then
                      !we need to apply hourly factors
                      call make_iland_for_time(debug_tfac, indate, i, j, iland, wday, iland_timefac,hour_iland,wday_loc,iland_timefac_hour)
                      tfac = fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)
                      if(Emis_source(n)%periodicity == 'yearly')then
                         !apply monthly factor on top of hourly factors
                         tfac = tfac * timefac(iland_timefac,sec2tfac_map(isec),iem)                         
                      endif
                   else
                      !not monthly or yearly emissions, timefactors must be included in emission values
                      tfac = 1.0
                   endif
                   
                   if (USES%GRIDDED_EMIS_MONTHLY_FACTOR) tfac=tfac* GridTfac(i,j,sec2tfac_map(isec),iem)

                   !Degree days - only SNAP-2 
                   if(USES%DEGREEDAY_FACTORS .and. &
                        sec2tfac_map(isec)==ISNAP_DOM .and. Gridded_SNAP2_Factors) then
                      oldtfac = tfac
                      ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
                      ! we make use of a baseload even for SNAP2
                      tfac = ( fac_min(iland,sec2tfac_map(isec),iem) & ! constant baseload
                           + ( 1.0-fac_min(iland,sec2tfac_map(isec),iem) )* gridfac_HDD(i,j) ) &
                           * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)
                      
                      if(debug_tfac .and. indate%hour==12 .and. iem==1) &
                           write(*,"(a,3i3,2i4,7f8.3)") "SNAPHDD tfac ",  &
                           isec, sec2tfac_map(isec),iland, daynumber, indate%hour, &
                           timefac(iland_timefac,sec2tfac_map(isec),iem), t2_nwp(i,j,2)-273.15, &
                           fac_min(iland,sec2tfac_map(isec),iem),  gridfac_HDD(i,j), tfac
                   end if ! =============== HDD 
    
                   s = Emis_source_2D(i,j,n) * emisfrac(iqrc,sec2split_map(isec),iland) * tfac
                   
                   SecEmisOut(i,j,iem,0) = SecEmisOut(i,j,iem,0) + s !sum of all sectors
                   if(SecEmisOutWanted(isec))SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) = &
                        SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) + s
                   ! Add up emissions in ktonne 
                   totemadd(itot) = totemadd(itot) + s * dtgrid * xmd(i,j)

                   !  Assign to height levels 1-KEMISTOP
                   do k=KEMISTOP,KMAX_MID
                      gridrcemis(iqrc,k,i,j) = gridrcemis(iqrc,k,i,j)   &
                           + s&
                           *ehlpcom0    &
                           *emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) &
                           *emis_masscorr(iqrc)
                   end do   ! k
                end do ! i
             enddo
          enddo
       endif
    enddo
 end if ! hourchange 
 


  ! We now scale gridrcemis to get emissions in molecules/cm3/s
!MOVED to setup_1d
!  do k= KEMISTOP, KMAX_MID
!    do j = 1,ljmax
!      do i = 1,limax
!        ehlpcom= roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))
!        !RB: This should also be done for the road dust emissions
!        do iqrc =1, NRCEMIS
!          gridrcemis(iqrc,k,i,j) =  gridrcemis0(iqrc,k,i,j)* ehlpcom
!        end do ! iqrc
!      end do   ! i
!    end do     ! j
!  end do       ! k

  if(USES%ROADDUST)THEN
    if(DEBUG%ROADDUST.and.debug_proc) &
      write(*,*)"Before the unit scaling",gridrcroadd(1:2,DEBUG_li,DEBUG_lj)
    do j = 1,ljmax
      do i = 1,limax
        ehlpcom= roa(i,j,KMAX_MID,1)/(ps(i,j,1)-PT)
        do iqrc =1, NROADDUST
          gridrcroadd(iqrc,i,j) =  gridrcroadd0(iqrc,i,j)* ehlpcom * ehlpcom0 &
              * roaddust_masscorr(iqrc)
        end do ! iqrc
      end do   ! i
    end do     ! j
    if(DEBUG%ROADDUST.and.debug_proc) &
      write(*,*)"After the unit scaling",gridrcroadd(1:2,DEBUG_li,DEBUG_lj)
  end if
end subroutine EmisSet
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

  integer :: i, j, k, iyr, iemislist, n, sec_ix
  character(len=200) :: fname
  real ktonne_to_kgm2s, tonnemonth_to_kgm2s  ! Units conversion
  integer :: iland, iem,ic,isec, i_gridemis
  real :: conv
  logical , save :: first_call=.true.

  ! For now, only the global runs use the Monthly files
  integer :: kstart,kend,nstart,Nyears
  real :: buffer(LIMAX,LJMAX),SumSoilNOx,ccsum
  real :: fractions(LIMAX,LJMAX,NCMAX),Reduc(NLAND)
  real, dimension(NEMIS_FILE)       :: emsum ! Sum emis over all countries
  real, dimension(NLAND,NEMIS_FILE) :: sumemis, sumemis_local ! Sum of emissions per country
  real, dimension(NEMIS_FILE) :: sumEU ! Sum of emissions in EU
  character(len=40) :: varname , fmt
  character(len=125) ::fileName
  real :: Mask_ReducFactor
  integer :: NMask_Code,Mask_Code(NLAND), i_femis_lonlat,icc
  real :: lonlat_fac, mw
  logical :: use_lonlat_femis, monthlysectoremisreset, Cexist
  logical :: fractionformat
  integer :: emis_inputlist_NEMIS_FILE!number of files for each emis_inputlist(i)
 ! For AIRCRAFT femis work
  integer, save  :: iFemis, iemNOx

  monthlysectoremisreset =.false.

  if(.not.allocated(airn).and.(USES%LIGHTNING_EMIS.or.USES%AIRCRAFT_EMIS))&
    allocate(airn(KCHEMTOP:KMAX_MID,LIMAX,LJMAX))


  if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)&
    call Read_monthly_emis_grid_fac(current_date%month)

  if(USES%AIRCRAFT_EMIS)then
    airn = 0.0
    kstart=KCHEMTOP
    kend=KMAX_MID

    if ( first_call ) then
        ! Scaling of aircraft. No real emission sector for aircraft (cc 900), so
        ! femis file should have 900 0 fnox etc..
        iFemis = find_index('AIRCRAFT',Country(:)%code )
        iemNOx=find_index("nox",EMIS_FILE(:)) ! save this index
        if(me==0.and.abs(1.0-e_fact(1,iFemis,iemNOx))>1.0E-2)&
             print *,iFemis,iemNOx," AIRCRAFT emission scaled by ",&
             e_fact(1,iFemis,iemNOx)
    end if

    call ReadField_CDF_FL(AircraftEmis_FLFile,'NOx',airn,&
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
             * e_fact(1,iFemis,iemNOx) &  ! sec1, emis1 AIRCRAFT
               /(dA(k)+dB(k)*ps(i,j,1))*xm2(i,j)
        end do
      end do
    end do
  end if

  if(DEBUG%SOILNOX.and.debug_proc) write(*,*)"Emissions DEBUG_SOILNOX ????"
  if(USES%EURO_SOILNOX)then  ! European Soil NOx emissions
    if(DEBUG%SOILNOX.and.debug_proc) write(*,*)"Emissions DEBUG_SOILNOX START"

    ! read in map of annual N-deposition produced from pre-runs of EMEP model
    ! with script mkcdo.annualNdep
    call ReadField_CDF(NdepFile,'Ndep_m2',AnnualNdep,1,&
          interpol='zero_order',needed=.true.,debug_flag=.false.,UnDef=0.0)

    if(DEBUG%SOILNOX.and.debug_proc)&
      write(*,"(a,4es12.3)") "Emissions_mod: SOILNOX AnnualDEBUG ", &
        AnnualNdep(debug_li, debug_lj), maxval(AnnualNdep), minval(AnnualNdep)

    call CheckStop(USES%GLOBAL_SOILNOX, "SOILNOX - cannot use global with Euro")
    ! We then calculate SoulNOx in Biogenics_mod

  elseif(USES%GLOBAL_SOILNOX) then ! Global soil NOx

    SoilNOx(:,:)=0.0      
    buffer(:,:)=0.0      
    nstart=(current_date%year-1996)*12 + current_date%month
    if(nstart>0.and.nstart<=120)then
      !the month is defined
      call ReadField_CDF(nox_emission_1996_2005File,'NOX_EMISSION',SoilNOx,&
             nstart=nstart,interpol='conservative',known_projection="lon lat",&
             needed=.true.,debug_flag=.false.,UnDef=0.0)
      if(DEBUG%SOILNOX.and.debug_proc) &
        write(*,*) "PROPER YEAR of SOILNO ", current_date%year, nstart
    else
      !the year is not defined; average over all years
      Nyears=10 !10 years defined
      do iyr=1,Nyears 
        nstart=12*(iyr-1) + current_date%month  
        call ReadField_CDF(nox_emission_1996_2005File,'NOX_EMISSION',buffer,&
              nstart=nstart,interpol='conservative',known_projection="lon lat",&
              needed=.true.,debug_flag=.false.,UnDef=0.0)
        do j=1,ljmax 
          do i=1,limax
            SoilNOx(i,j)=SoilNOx(i,j)+buffer(i,j)
          end do
        end do
        if(DEBUG%SOILNOX.and.debug_proc) &
          write(*,"(a,2i6,es10.3,a,2es10.3)") "Averaging SOILNO  inputs", &
            1995+(iyr-1), nstart,SoilNOx(debug_li, debug_lj), &
            "max: ", maxval(buffer), maxval(SoilNOx)
      end do
      SoilNOx=SoilNOx/Nyears
    end if ! nstart test

     if(DEBUG%SOILNOX.and.debug_proc) then
        write(*,"(a,i3,3es10.3)") "After Global SOILNO ",&
             me,maxval(SoilNOx),SoilNOx(debug_li,debug_lj)
       !write(*,"(a,i3,3es10.3)") "After Global SOILNO ",  me, maxval(SoilNOx), SoilNOx(3, 3)
     end if
  else ! no soil NO
    if(DEBUG%SOILNOX.and.debug_proc) &
      write(*,*) "Emissions DEBUG_SOILNOX - none"
  end if !  SOIL NO

  !for testing, compute total soil NOx emissions within domain
  !convert from g/m2/day into kg/day
  if(USES%GLOBAL_SOILNOX) then 
    SumSoilNOx=0.0
    SoilNOx = max(0.0, SoilNOx)  ! Stops the NEGs!
    do j=1,ljmax
      do i=1,limax      
        SumSoilNOx=SumSoilNOx+0.001*SoilNOx(i,j)*gridwidth_m**2*xmd(i,j)      
      end do
    end do
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
      end do
    end do
  end if

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

  if(MasterProc.and.DEBUG%EMISSIONS ) then
    write(*,*) 'Enters newmonth, mm, ktonne_to_kgm2s = ', &
         current_date%month,ktonne_to_kgm2s
    write(*,*) ' first_dms_read = ', first_dms_read
  end if

  tonnemonth_to_kgm2s = 1.0e3 &
       /(nmdays(current_date%month)*24.*60.*60.*GRIDWIDTH_M*GRIDWIDTH_M)
 
  do iemislist = 1, size( emis_inputlist(:)%name )

    if(emis_inputlist(iemislist)%name == "NOTSET")cycle
     
    !periodicity set in routine Emissions(year) if 12 records are found
    if(emis_inputlist(iemislist)%periodicity /= "monthly")cycle
    if(MasterProc)write(*,*)'reading monthly emissions for ',trim(emis_inputlist(iemislist)%name)

    use_lonlat_femis = emis_inputlist(iemislist)%use_lonlat_femis

    if((emis_inputlist(iemislist)%type == "sectors" .or.&
         emis_inputlist(iemislist)%type == "GNFRsectors" .or.&
         emis_inputlist(iemislist)%type == "SNAPsectors" )&
         .and. index(emis_inputlist(iemislist)%name,".nc")>1)then 

       !Read monthly emission sector files      
       emis_inputlist_NEMIS_FILE = 1 !all pollutants are in same file
       if(index(emis_inputlist(iemislist)%name,"POLL")>0)emis_inputlist_NEMIS_FILE = NEMIS_FILE !one file per pollutant
       sumemis_local = 0.0 !needed to get total right
       do iem = 1, emis_inputlist_NEMIS_FILE
          sumemis_local(:,iem) = 0.0
          if(emis_inputlist(iemislist)%pollName(1)/='NOTSET')then
             if(all(emis_inputlist(iemislist)%pollName(:)/=trim(EMIS_FILE(iem))))cycle      
             if(Masterproc)write(*,*)'using PollNames restrictions '
          end if
          !reset emissions only once each month
          if(.not.monthlysectoremisreset)then
             secemis = 0.0
             nlandcode=0
             monthlysectoremisreset = .true.
             Emis_field = 0.0
             Emis_id(:)%species = '' !to be removed, when searching only among NEmis_id 
             NEmis_id = 0
          endif
          fractionformat = ( emis_inputlist(iemislist)%format=='fractions' )
          fname = key2str(trim(emis_inputlist(iemislist)%name),'POLL',EMIS_FILE(iem))
          
          call EmisGetCdf(iem, fname, sumemis_local, &
               secemis, landcode, nlandcode, current_date%month,&
               emis_inputlist(iemislist)%incl, nin, emis_inputlist(iemislist)%excl, nex, &
               emis_inputlist(iemislist)%use_lonlat_femis,&
               emis_inputlist(iemislist)%set_mask,emis_inputlist(iemislist)%use_mask,&
               emis_inputlist(iemislist)%pollName,&
               fractionformat,emis_inputlist(iemislist)%type)
          
       end do! iem = 1,NEMIS_FILE
       sumemis=0.0
       CALL MPI_REDUCE(sumemis_local,sumemis,&
            NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)        
       if(MasterProc)then
          call PrintLog("Total emissions by countries for "//trim(emis_inputlist(iemislist)%name)//" (Gg)")
          write(*     ,"(a4,a9,3x,30(a12,:))")" CC ","    ",EMIS_FILE(:)
          write(IO_LOG,"(a4,a9,3x,30(a12,:))")" CC ","    ",EMIS_FILE(:)                
          sumEU(:) = 0.0
          fmt="(i4,1x,a9,3x,30(f12.2,:))"
          do ic = 1, NLAND
             ccsum = sum( sumemis(ic,:) )
             icc=Country(ic)%icode
             if ( ccsum > 0.0 )then
                write(*,     fmt) icc, Country(ic)%code, sumemis(ic,:)
                write(IO_LOG,fmt) icc, Country(ic)%code, sumemis(ic,:)
             end if
             if(find_index(Country(ic)%code,EU28(:))>0) sumEU = sumEU + sumemis(ic,:)
          end do
          if ( sum(sumEU(:))>0.001) then
             write(*     ,fmt) 0, "EU", sumEU(:)
             write(IO_LOG,fmt) 0, "EU", sumEU(:)
          end if
!          do ic = 1, NLAND
!             ccsum = sum( sumemis(ic,:) )
!             if(ccsum>0.0 .and. MasterProc ) then
!                write(*,"(i5,1x,a10,3x,30(f12.2,:))")ic, Country(ic)%code, 1000*sumemis(ic,:)
!                write(IO_LOG,"(i5,1x,a10,3x,30(f12.2,:))")ic, Country(ic)%code, 1000*sumemis(ic,:)
!             end if
!          end do
       end if
       
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
          end do
       end do
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
   end if
  end do !iemislist

  if(monthlysectoremisreset)then
     !unit conversion
     do iem = 1,NEMIS_FILE
        do j=1,ljmax
           do i=1,limax
              do n=1,nlandcode(i,j)
                 do  isec = 1,NSECTORS                    
                    secemis(isec,i,j,n,iem)=secemis(isec,i,j,n,iem)*tonnemonth_to_kgm2s*xm2(i,j)
                 end do
              end do
           end do
        end do
     end do
     do iem = 1,NEmis_id
        do j=1,ljmax
           do i=1,limax
              Emis_field(i,j,iem) = Emis_field(i,j,iem)*tonnemonth_to_kgm2s*xm2(i,j)
           end do
        end do
     end do
     

  endif

  if(foundYearlySectorEmissions .and. monthlysectoremisreset)then
     !merge with yearly emissions (stored in GridEmis)
     if(me==0)write(*,*)'Merging yearly sector emissions with monthly emissions'
     if(me==0)write(*,*)'WARNING: Correcting for monthly timefactor applied on monthly emissions!'
     do iem = 1,NEMIS_FILE
        do j=1,ljmax
           do i=1,limax         
              !monthly emissions should not add an additional monthly time factor
              do n=1,nlandcode(i,j)
                 iland=find_index(landcode(i,j,n),Country(:)%icode) !array index
                 !array index of country that is used as reference for timefactor
                 ic = find_index(Country(iland)%timefac_index,Country(:)%icode)
                 if(ic>0)then
                    !cancel monthly timefactor
                    secemis(:,i,j,n,iem)=secemis(:,i,j,n,iem)/(1.E-10+fac_emm(ic,current_date%month,:,iem))
                 else
                    write(*,*)'warning: Country code from monthly sector emissions not recognized: ',landcode(i,j,n),me,i,j,n
                 end if
              end do

              !merge
              do i_gridemis=1,nGridEmisCodes(i,j)
                 ic=find_index(GridEmisCodes(i,j,i_gridemis),Country(:)%icode)
                 !merge other (already read in) emissions into secemis
                 Cexist=.false.
                 do n=1,nlandcode(i,j)
                    if(GridEmisCodes(i,j,i_gridemis)==landcode(i,j,n))then                           
                       secemis(:,i,j,n,iem)=secemis(:,i,j,n,iem)&
                            +GridEmis(:,i,j,i_gridemis,iem)      
                       Cexist=.true.
                       exit
                    end if
                 end do
                 if(.not.Cexist)then
                    !country not included yet. define it now:
                    nlandcode(i,j)=nlandcode(i,j)+1
                    if(nlandcode(i,j)>NCMAX)then
                       write(*,*)"Too many emitter countries in one gridemiscell: ",&
                            me,i,j,nGridEmisCodes(i,j)
                       call StopAll("To many countries in one gridemiscell ")
                    end if
                    n=nlandcode(i,j)
                    landcode(i,j,n)=GridEmisCodes(i,j,i_gridemis)
                    secemis(:,i,j,n,iem)=GridEmis(:,i,j,i_gridemis,iem)
                 end if
              end do
           end do
        end do
     end do! iem = 1,NEMIS_FILE
  endif
  

  first_call=.false.
end subroutine newmonth
!***********************************************************************
subroutine EmisWriteOut(label, iem,nsources,sources,emis)
!----------------------------------------------------------------------!
! Ascii output of emissions fields (after interpolation and femis.
! To print out emissions, we need to get fields for each country
! in turn, needing info from GridCodes
! e.g. 11-SNAP data is secemis:
!  real    secemis (NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILES)
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

  associate ( dbg => DEBUG%EMISSIONS )
  if(dbg ) write(*,*)"CALLED "//trim(txt),me,&
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
      end if
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
                if(dbg ) call CheckStop(any(locemis(i,j,:)< 0.0),"NEG LOCEMIS")
              end if
            end do
          end do
        end do ! j
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
            end if
          end do
        end do
      end do EMLAND
      
      do isec = 1, NSECTORS
         lemis = locemis(:,:,isec)
         if(dbg.and.debug_proc) write(*,*) trim(txt)//" lemis ",&
             me,iland,isec,maxval(lemis(:,:))
      end do ! isec
      
    end if
    close(IO_TMP)
    CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)!wait: one should write at a time
  end do
  
!  deallocate(locemis,lemis)
  end associate ! dbg

end subroutine EmisWriteOut

end module Emissions_mod
