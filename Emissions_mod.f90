! <Emissions_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.6>
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
module Emissions_mod
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Calls up emission read/set routines
!  This routine interfaces the stand-alone emission-file reading routines
!  with the 3D model.
!_____________________________________________________________________________

use AirEmis_mod, only : airn, TotAircraftEmis
use Biogenics_mod,     only: SoilNOx, SoilNOx3D, AnnualNdep
use CheckStop_mod,     only: CheckStop,StopAll
use ChemDims_mod,      only: NSPEC_SHL, NSPEC_TOT,&
                             NEMIS_File  ! No. emission files
use ChemSpecs_mod,     only: NO2, SO2,species,species_adv
use Chemfields_mod,    only: xn_adv
use Config_module,only: &
    KMAX_MID, KMAX_BND, PT ,dt_advec, step_main, &
    KCHEMTOP, &
    NSECTORS_ADD_MAX, SECTORS_ADD, &
    timeFacs, &  ! e.g. timeFacs%Monthly == 'GRIDDED') &
    emis_inputlist, &
    EmisDir,      &    ! template for emission path
    DataDir,      &    ! template for path
    EMIS_OUT,      &    ! output emissions in ASCII or not
    INERIS_SNAP2 , &    ! INERIS/TFMM HDD20 method
    MasterProc, USES,  &  !
    SEAFIX_GEA_NEEDED, &  !  see below
    DMS,DMS_S_FAC,&
    NPROC, EmisSplit_OUT,&
    SecEmisTotalsWanted,SecEmisOutWanted,MaxNSECTORS,&
    AircraftEmis_FLFile,soilnox_emission_File, RoadMapFile,&
    DMSFile,OceanNH3File,&
    AVG_SMI_2005_2010File,NdepFile,&
    startdate, Emis_sourceFiles, EmisMask, NEmisMaskMAX,&
    hourly_emisfac
use Country_mod,       only: MAXNLAND,NLAND,Country,IC_NAT,IC_FI,IC_NO,IC_SE, IC_HU, IC_AT
use Country_mod,       only: EU28,EUMACC2,IC_DUMMY
use Debug_module,      only: DEBUG ! , & !DEBUG => DEBUG_EMISSIONS, &
                                !DEBUG%EMISTIMEFACS
use EmisDef_mod,       only: &
      EMIS_FILE     & ! Names of species ("sox  ",...)
     ,NCMAX         & ! Max. No. countries per grid
     ,TFAC_IDX_POW  & ! tfac index for public power emis
     ,TFAC_IDX_DOM  & ! tfac index for domestic/resid emis
     ,TFAC_IDX_TRAF & ! tfac index for road-traffic (SNAP7)
     ,TFAC_IDX_AGR  & ! tfac index for agriculture
     ,IS_POW,IS_AGR,IS_TRAF,IS_DOM,IS_IND & ! Also used for special sector managment
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
     ,MAXFEMISLONLAT,N_femis_lonlat &
     ,NSECTORS, N_HFAC, N_TFAC, N_SPLIT     & ! No. emis sectors, height, time and split classes
     ,NSECTORS_GNFR_CAMS, GNFR_CAMS_SECTORS, NSECTORS_SNAP, SNAP_SECTORS, NSECTORS_MAX, SECTORS&
     ,foundYearlySectorEmissions, foundMonthlySectorEmissions&
     ,Emis_mask, Emis_mask_allocate, MASK_LIMIT & !old format
     ,NEmisMask, EmisMaskValues, EmisMaskIntVal,EmisMaskIndex2Name& !new format
     ,Emis_field, Emis_id, NEmis_id &
     ,NEmisFile_sources, EmisFiles,NEmis_sources, Emis_source&
     ,Emis_source_3D,ix3Dmap, NEmis_3Dsources&
     ,NEmis_source_ij, Emis_source_ij,Emis_source_ij_ix,NEmis_source_ijMAX
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
use Io_Progs_mod,      only: ios, open_file, datewrite
use Io_RunLog_mod,     only: PrintLog
use MetFields_mod,     only: u_xmj, v_xmi, roa, ps, z_bnd, surface_precip,EtaKz ! ps in Pa, roa in kg/m3
use MetFields_mod,     only: t2_nwp   ! DS_TEST SOILNO - was zero!
use MPI_Groups_mod
use NetCDF_mod,        only: ReadField_CDF,ReadField_CDF_FL,ReadTimeCDF,IsCDFfractionFormat,&
                             GetCDF_modelgrid,PrintCDF,ReadSectorName,check,&
                             create_country_emission_file, output_country_emissions
use netcdf
use OwnDataTypes_mod,  only: TXTLEN_FILE, TXTLEN_NAME,Emis_id_type, &
                             EmisFile_id_type, Emis_sourceFile_id_type, print_Sector_type
use Par_mod,           only: MAXLIMAX,MAXLJMAX, GIMAX,GJMAX, IRUNBEG,JRUNBEG,&
                            me,limax,ljmax, MSG_READ1,MSG_READ7&
                           ,gi0,gj0,li0,li1,lj0,lj1
use PhysicalConstants_mod,only: GRAV, AVOG, ATWAIR
use PointSource_mod,      only: readstacks !MKPS
use ZchemData_mod,only: rcemis   ! ESX
use SmallUtils_mod,    only: find_index, key2str, trims, basename
use TimeDate_mod,      only: nydays, nmdays, date, current_date,&! No. days per
                            tdif_secs,timestamp,make_timestamp,daynumber,day_of_week ! year, date-type, weekday
use TimeDate_ExtraUtil_mod,only :nctime2date, date2string, to_idate,date_is_reached
use Timefactors_mod,   only: &
     NewDayFactors          & ! subroutines
    ,DegreeDayFactors       & ! degree-days used for SNAP-2
    ,Gridded_SNAP2_Factors, gridfac_HDD &
    ,fac_min,timefactors   &                  ! subroutine
    ,fac_ehh24x7 ,fac_emm, fac_cemm, fac_edd, timefac & ! time-factors
    ,fac_dayofyear & ! time-factors
    ,Read_monthly_emis_grid_fac &
    ,Read_monthly_timezones &
    ,GridTfac &!array with monthly gridded time factors
    ,yearly_normalize !renormalize timefactors after reset
use LocalFractions_mod, only : save_lf_emis,loc_frac, emis_lf_cntry
implicit none
private

! subroutines:
public :: Init_masks
public :: Init_Emissions    ! defines emission setup and formats
public :: EmisUpdate        ! update emission
public :: Emissions         ! Main emissions module
public :: newmonth
public :: EmisSet           ! Sets emission rates every hour
public :: EmisWriteOut           ! Outputs emissions in ascii

! The main code does not need to know about the following
private :: expandcclist            !  expands e.g. EU28, EUMACC2
private :: consistency_check       ! Safety-checks


logical, save, private  :: first_dms_read


! and for budgets (not yet used - not changed dimension)
real, public,  save, dimension(NSPEC_SHL+1:NSPEC_TOT) ::  totemadd
integer, private, save :: iemCO  ! index of CO emissions, for debug

real ::TimesInDays(120),mpi_out
integer ::NTime_Read=-1, found
logical :: fileExists            ! to test emission files
integer :: nin, nex         !  include, exclude numbers for emis_inputlist
character(len=*), parameter :: sub='Emissions:'

contains
!***********************************************************************
  !new (Nov 2018) emission setup and formats + sectors for old format
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

    ! Emis_sourceFiles are read from config
    ! EmisFiles and Emis_source are defined by the model, in Emis_init_GetCdf
    ! Emis_source_ij are the actual emission fields

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
    integer :: i, j, ii, n, nn, ix, nemis_old, isource
    integer :: isec, isec_idx, iland, iem, iqrc, itot, f
    integer :: startsource(size(Emis_sourceFiles)), endsource(size(Emis_sourceFiles))
    type(Emis_id_type):: Emis_id_undefined !to get undefined source values
    type(EmisFile_id_type):: Emisfile_undefined !to get undefined file values
    !Note must distinguish between default and not undefined values, to know when config has defined a value, even if it is the default.
    type(Emis_id_type):: Emis_sources_defaults !set values when not specified otherwise
    type(EmisFile_id_type):: Emisfile_defaults !set values when not specified otherwise
    integer :: EmisFilesMap(0:size(Emis_sourceFiles)) !index of EmisFile given index of Emis_sourceFiles
    integer :: max_levels3D
    character(len=*),parameter :: dtxt='Init_Emissions: '
    logical, save  :: first_call=.true., dbg, dbgX

    if ( first_call ) then
       dbg = DEBUG%EMISSIONS .and. MasterProc
       dbgX = DEBUG%EMISSIONS2 .and. MasterProc
       first_call = .false.
    end if

    !S24: EmisSplit MOVED HERE:
    call EmisSplit()    ! In EmisGet_mod, => emisfrac
    if(MasterProc) write(*,*) dtxt//'TRACK from EmisSplit:', (emis_nsplit(i),i=1,NEMIS_FILE)
    call CheckStop(ios, "ioserror: EmisSplit")

    if(USES%OCEAN_NH3)then
       if(MasterProc)write(*,*)' using  OceanNH3'
       O_NH3%index=find_index("NH3",species(:)%name)
       call CheckStop(O_NH3%index<0,'NH3 not found. Needed for OceanNH3 emissions')
       allocate(O_NH3%emis(LIMAX,LJMAX))
       allocate(O_NH3%map(LIMAX,LJMAX))
       O_NH3%emis=0.0
       O_NH3%map=0.0
       O_NH3%sum_month=0.0
       O_NH3%sum_year=0.0
    end if

    if (USES%OCEAN_DMS)then
       if(MasterProc)write(*,*)'using DMS'
       O_DMS%index=find_index("SO2",species(:)%name)
       call CheckStop(O_DMS%index<0,'SO2 not found. Needed for DMS emissions')
       allocate(O_DMS%emis(LIMAX,LJMAX))
       allocate(O_DMS%map(LIMAX,LJMAX))
       O_DMS%emis=0.0
       O_DMS%map=0.0
       O_DMS%sum_month=0.0
       O_DMS%sum_year=0.0
    endif

    ! new format initializations

    !1) define lowest level default values
    Emis_sources_defaults%units = 'mg/m2/h'
    Emis_sources_defaults%country_ISO = 'N/A'
    Emis_sources_defaults%countrycode = -1
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
    Emisfile_defaults%timevalidity = 'end'
    Emisfile_defaults%projection = 'lon lat'
    Emisfile_defaults%factor = 1.0
    Emisfile_defaults%sectorsName = 'GNFR_CAMS'
    Emisfile_defaults%units = Emis_sources_defaults%units
    EmisFile_defaults%sector = Emis_sources_defaults%sector
    EmisFile_defaults%country_ISO = Emis_sources_defaults%country_ISO
    EmisFile_defaults%countrycode = Emis_sources_defaults%countrycode

    EmisFiles(:) = EmisFile_defaults


    !2) read from global attributes in file
    !3) read from variable attributes in file
    !Emis_sourceFiles is read from config and then not modified
    !EmisFiles collect all valid data and sources
    !Emis_source() contains the metadata of the emissions sources finally used
    !Emis_source_ij contains the non zero values for each source used
    NEmis_sources = 0 !total number of valid emis (also 3D) variables across all files
    NEmis_3Dsources = 0 !total number of valid 3D emis variables across all files
    max_levels3D = 1
    NEmisFile_sources = 0 !number of valid emission files
    EmisFilesMap = 0 !index of EmisFile given index of EmisFile_sources
    do n = 1, size(Emis_sourceFiles)
       if(dbgX) write(*,"(a,3i6,a)") dtxt//'XNEmisHuntP', n, NEmis_sources, size(Emis_sourceFiles), trim(Emis_sourceFiles(n)%filename) 
       nemis_old = NEmis_sources
       if(Emis_sourceFiles(n)%filename/='NOTSET')then
          !Read all variables and set parameters as needed
          !find which variables names are defined in config
          names_in='NOTSET'
          i=0
          if(dbgX) write(*,"(a,3i6,a)") dtxt//'XNEmisHuntQ', size(Emis_sourceFiles(n)%source)
          do nn = 1, size(Emis_sourceFiles(n)%source)
             !if(dbgX) write(*,"(a,i6,a)") dtxt//'XNEmisHuntQQ',nn, trim(Emis_sourceFiles(n)%source(nn)%varname)
             if(Emis_sourceFiles(n)%source(nn)%varname/='NOTSET')then
                i= i + 1
                names_in(i)=trim(Emis_sourceFiles(n)%source(nn)%varname)
                if(dbgX) write(*,"(a,2i6,a)") dtxt//'XNEmisHuntQQi',nn,i, trim(Emis_sourceFiles(n)%source(nn)%varname)
             endif
          enddo

          if(MasterProc)write(*,*)dtxt//"Initializing Emissions from "//&
             trim(Emis_sourceFiles(n)%filename)
          call Emis_init_GetCdf(Emis_sourceFiles(n), EmisFiles(NEmisFile_sources+1), names_in, i)
       else
          !filename/='NOTSET'
          cycle
       endif
       startsource(n) = nemis_old + 1
       endsource(n) = NEmis_sources
       if(dbgX) write(*,"(a,3i6,a)") dtxt//'XNEmisHuntR', NEmis_sources, nemis_old
       if(NEmis_sources>nemis_old)then
          EmisFilesMap(NEmisFile_sources) = n
          EmisFiles(NEmisFile_sources)%Nsources =  NEmis_sources - nemis_old
          EmisFiles(NEmisFile_sources)%source_start =  nemis_old + 1
          EmisFiles(NEmisFile_sources)%source_end =  NEmis_sources
          if(dbgX) write(*,*) dtxt//'XNEmisHuntS', n, NEmis_sources
       else
          if(me==0)write(*,*)dtxt//'WARNING: did not find any valid source in '&
                      //trim(Emis_sourceFiles(n)%filename)
       endif
    enddo
    allocate(Emis_source_ij(LIMAX*LJMAX,NEmis_source_ijMAX))
    allocate(Emis_source_ij_ix(LIMAX*LJMAX,NEmis_source_ijMAX))
    allocate(NEmis_source_ij(LIMAX*LJMAX))
    NEmis_source_ij=0

    !4)overwrite parameters with settings from config_emep.nml if they are set
    !first overwrite the global attributes: projection and periodicity
    !use associate here?
    do i = 1, size(Emis_sourceFiles)
       n = EmisFilesMap(i)
       if(n>0)then

         associate ( emn => Emis_sourceFiles(n),  emUNDEF => Emisfile_undefined, emfi => EmisFiles(i) )

          if(emn%periodicity /= emUNDEF%periodicity)   emfi%periodicity  = emn%periodicity
          if(emn%timevalidity /= emUNDEF%timevalidity) emfi%timevalidity = emn%timevalidity
          if(emn%projection /= emUNDEF%projection)     emfi%projection   = emn%projection
          if(emn%grid_resolution /= emUNDEF%grid_resolution) emfi%grid_resolution = emn%grid_resolution
          if(emn%projection /= 'native')then
            call CheckStop(emfi%grid_resolution <=1.0E-5,&
             'Grid_resolution must be defined for '//trim(emn%filename))
          endif
          if(emn%mask_ID         /= emUNDEF%mask_ID)         emfi%mask_ID = emn%mask_ID
          if(emn%mask_ID_reverse /= emUNDEF%mask_ID_reverse) emfi%mask_ID_reverse = emn%mask_ID_reverse
          if(emn%species         /= emUNDEF%species)         emfi%species = emn%species
          if(emn%units           /= emUNDEF%units)           emfi%units = emn%units
          if(emn%countrycode     /= emUNDEF%countrycode)     emfi%countrycode = emn%countrycode
          if(emn%country_ISO     /= emUNDEF%country_ISO)     emfi%country_ISO = emn%country_ISO
          if(emn%sector          /= emUNDEF%sector)          emfi%sector = emn%sector
          if(emn%sectorsName     /= emUNDEF%sectorsName)     emfi%sectorsName = emn%sectorsName
          if(emn%factor          /= emUNDEF%factor)          emfi%factor = emn%factor
          ! correct for old notations
          if(emfi%sectorsName == 'SNAPsectors') emfi%sectorsName = 'SNAP'
          if(emfi%sectorsName == 'GNFRsectors') emfi%sectorsName = 'GNFR_CAMS'
          if(emfi%sectorsName == 'GNFR') emfi%sectorsName = 'GNFR_CAMS'

        end associate ! emn => Emis_sourceFiles(n), emUNDEF => Emisfile_undefined, emfi => EmisFiles(i)

       endif
    enddo

    !then overwrite the variable attributes
    EMFILE_ILOOP: do i = 1, NEmisFile_sources !loop over files
       n = EmisFilesMap(i)
       found = 0
       EMFILE_IILOOP: do ii = EmisFiles(i)%source_start, EmisFiles(i)%source_end !loop over sources found in the netcdf file
          ! set source default = file parameter if they are set. Defines default
          ! for sources in this file, can be also be redefined for individual sources

          ! - use associate to make this bit more readable ;-)
          associate ( emn => Emis_sourceFiles(n), emundef=>Emisfile_undefined, emsii => Emis_source(ii) )

          if(emn%species     /= emundef%species)     emsii%species = emn%species
          if(emn%units       /= emundef%units)       emsii%units = emn%units
          if(emn%countrycode /= emundef%countrycode) emsii%countrycode = emn%countrycode
          if(emn%country_ISO /= emundef%country_ISO) emsii%country_ISO = emn%country_ISO
          if(emn%sector      /= emundef%sector)      emsii%sector = emn%sector
         !NB: periodicity and timevalidity cannot be set individually for variables
          emsii%periodicity       = EmisFiles(i)%periodicity
          emsii%timevalidity      = EmisFiles(i)%timevalidity
         ! Following cannot be set in netcdf file:
          emsii%apply_femis       = emn%apply_femis
          emsii%include_in_local_fractions = emn%include_in_local_fractions
          emsii%mask_ID           = emn%mask_ID
          emsii%mask_ID_reverse   = emn%mask_ID_reverse

          isource = emsii%ix_in
          if(dbg) write(*,'(a,5i6)') dtxt//'writing config attribute on '//trim(emsii%varname),i, n, ii, isource, NEmis_sources
          if(dbg.and.ii==NEmis_sources) write(*,*)dtxt//'Last ii', emsii
          if(isource>0)then
             !source defined in config file
             if(trim(emn%source(isource)%varname)/=trim(emsii%varname))&
              write(*,*)isource,dtxt//'ERROR', &
                trim(emn%source(isource)%varname),' ',trim(emsii%varname),ii

             found = 1

             associate ( emi => Emis_sourceFiles(n)%source(isource) )

             if(emi%species /= Emis_id_undefined%species) emsii%species = emi%species
             if(emi%units /= Emis_id_undefined%units) emsii%units = emi%units
             if(emi%sector /= Emis_id_undefined%sector) emsii%sector = emi%sector
             if(emi%factor /= Emis_id_undefined%factor) emsii%factor = emi%factor
             if(emi%countrycode /= Emis_id_undefined%countrycode) emsii%countrycode = emi%countrycode
             if(emi%country_ISO /= Emis_id_undefined%country_ISO) emsii%country_ISO = emi%country_ISO
             if(.not. emi%include_in_local_fractions) emsii%include_in_local_fractions = .false.
             if(.not. emi%apply_femis) emsii%apply_femis = emi%apply_femis

             if(emi%mask_ID /= Emis_id_undefined%mask_ID) emsii%mask_ID = emi%mask_ID
             if(emi%mask_ID_reverse /= Emis_id_undefined%mask_ID_reverse) emsii%mask_ID_reverse = emi%mask_ID_reverse

             if(emi%is3D .neqv. Emis_id_undefined%is3D) emsii%is3D   = emi%is3D
             if(emi%istart /= Emis_id_undefined%istart) emsii%istart = emi%istart
             if(emi%jstart /= Emis_id_undefined%jstart) emsii%jstart = emi%jstart
             if(emi%kstart /= Emis_id_undefined%kstart) emsii%kstart = emi%kstart
             if(emi%kend /= Emis_id_undefined%kend)     emsii%kend   = emi%kend
             if(emi%reversek .neqv. Emis_id_undefined%reversek)   emsii%reversek = emi%reversek
             if(emi%injection_k /= Emis_id_undefined%injection_k) emsii%injection_k = emi%injection_k

             end associate ! emn 

          endif

          ix = find_index(trim(emsii%country_ISO) ,Country(:)%code, first_only=.true.)
          if(trim(emsii%country_ISO)=="N/A")then
             ix = find_index(emsii%countrycode ,Country(:)%icode, first_only=.true.)
          end if
          if(ix<0)then
             if(me==0)write(*,*)dtxt//'WARNING: country '//&
                     trim(emsii%country_ISO)//' not defined for '//trim(emsii%varname)
             ix = find_index("N/A" ,Country(:)%code, first_only=.true.)
          else
             if(dbg)write(*,*)dtxt//'country found '//trim(Country(ix)%code), ix, NEmis_sources
          endif
          emsii%country_ix = ix
          emsii%country_ISO = trim(Country(ix)%code)


          ix = find_index(emsii%species, EMIS_FILE(:))
          if ( ix> 0) then
            if(dbg)write(*,*)dtxt//'spec found in EMIS_FILE'//trim(emsii%species), EMIS_FILE(ix)
          else    ! if( ix<=0 )then !find if it is defined as an individual species
             ix = find_index(emsii%species, species(:)%name )
             if(ix<=0)then
                !try case insensitive too
                ix = find_index(emsii%species, species(:)%name, any_case=.true.)
                if(dbg)write(*,*)dtxt//'spec found in species?'//trim(emsii%species), ix
                if(ix>0)then
                   if(me==0)write(*,*)'WARNING: '//trim(emsii%species)//&
                        ' not found, replacing with '//trim(species(ix)%name)
                   emsii%species = trim(species(ix)%name)
                end if
             end if
             if(ix>0)then
                emsii%species_ix = ix
                if(dbg)write(*,'(a,i4,a)')dtxt//' species found '// &
                     trim(emsii%country_ISO), ix, ' '//trim(species(ix)%name)
                if(emsii%include_in_local_fractions .and. USES%LocalFractions )then
                   if(me==0)write(*,*)"WARNING: local fractions will not include single species "//emsii%species
                endif
             else ! ix<=0
                if(dbg)write(*,'(a,i4,a)')dtxt//' species not found'// &
                     trim(emsii%country_ISO),ix,trim(emsii%species)
             endif
          endif


          if (emsii%sector>0) then
             !Add the relevant sector in SECTORS
             !look if it is already defined
             found = 0
             do isec_idx = 1, NSECTORS
                if (SECTORS(isec_idx)%name ==  trim(EmisFiles(i)%sectorsName)) found = 1
                if (found > 0) exit
             end do
             if (found == 0) then
                !add this sector in SECTORS
                do isec_idx = 1, NSECTORS_ADD_MAX
                   if(SECTORS_ADD(isec_idx)%name ==  trim(EmisFiles(i)%sectorsName)) then
                      NSECTORS = NSECTORS + 1
                      call CheckStop(NSECTORS > NSECTORS_MAX, "NSECTORS_MAX too small, please increase value in EmisDef_mod.f90")
                      SECTORS(NSECTORS) = SECTORS_ADD(isec_idx)
                      if(MasterProc) write(*,*) dtxt//'adding sector1st  ',&
                         trim(SECTORS_ADD(isec_idx)%longname),', total ',NSECTORS
                      found = 1
                   end if
                end do
                if (found == 0 .and. trim(EmisFiles(i)%sectorsName) == 'GNFR_CAMS') then ! note that default values are used only if not defined in config
                   do isec_idx = 1, NSECTORS_GNFR_CAMS
                      NSECTORS = NSECTORS + 1
                      call CheckStop(NSECTORS > NSECTORS_MAX, "NSECTORS_MAX too small, please increase value in EmisDef_mod.f90")
                      SECTORS(NSECTORS) = GNFR_CAMS_SECTORS(isec_idx)
                   end do
                end if
                if (found == 0 .and. trim(EmisFiles(i)%sectorsName) == 'SNAP') then ! note that default values are used only if not defined in config
                   do isec_idx = 1, NSECTORS_SNAP
                      NSECTORS = NSECTORS + 1
                      call CheckStop(NSECTORS > NSECTORS_MAX, "NSECTORS_MAX too small, please increase value in EmisDef_mod.f90")
                      SECTORS(NSECTORS) = SNAP_SECTORS(isec_idx)
                   end do
                end if
             end if

             !save the corresponding index in SECTORS
             found = 0
             do isec_idx = 1, NSECTORS
                if (SECTORS(isec_idx)%name ==  trim(EmisFiles(i)%sectorsName)) found = isec_idx
                if(dbg)write(*,'(a,i4)')dtxt//'sec srch:'//trim(SECTORS(isec_idx)%name)//':'//trim(EmisFiles(i)%sectorsName), found
                if (found > 0) exit; ! we want the first entry
             end do
             if (found == 0) call StopAll(trim(emn%filename)//&
              ': sectorsName not recognized! '//trim(EmisFiles(i)%sectorsName))
             emsii%sector_idx = found -1 + emsii%sector !TODO: make more robust (use names, not indices)
             if(dbg)write(*,"(a,2i4,a)")dtxt//'sec found2', emsii%sector,&
                emsii%sector_idx, trim(emsii%varname)//';'//trim(EmisFiles(i)%sectorsName)
          end if

          max_levels3D=max(max_levels3D, emsii%kend - emsii%kstart + 1)
          if (MasterProc .and. dbg) then
             write(*,*)dtxt//"File(i) emission source parameters ", i, EmisFiles(i)%nsectors!,emsii
             if (EmisFiles(i)%nsectors == 1) then
                print *, "FILEI", i, trim(emsii%varname), emsii%sector_idx
                write(*,'(a,a,a)')dtxt//'emsrc include1 '//trim(emsii%varname)//' as '//&
                     trim(emsii%species)//' sector ',&
                     trim(SECTORS(emsii%sector_idx)%longname),&
                     ' country '//trim(emsii%country_ISO)
             else if (emsii%sector == EmisFiles(i)%nsectors) then
                write(*,'(a,I3,a,a)')dtxt//'emsrc include2 '//trim(emsii%varname)//' as '//&
                 trim(emsii%species),EmisFiles(i)%nsectors,&
                 ' sectors ',trim(SECTORS(Emis_source(ii-EmisFiles(i)%nsectors+1)%sector_idx)%longname)//&
                  ' to '// trim(SECTORS(emsii%sector_idx)%longname)//&
                  ' country '//trim(emsii%country_ISO)

             end if
          end if
          end associate ! emn, emsii
       enddo EMFILE_IILOOP
    enddo EMFILE_ILOOP

    do n = 1, NEmis_sources
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
    integer :: i,ii,jj,ic,iEmisMask,EmisMask_ix(NEmisMaskMAX)

    iEmisMask = 0
    !1) find number of valid masks defined
    do i = 1, size(EmisMask)
       if(EmisMask(i)%filename == 'NOTSET') cycle
       call CheckStop(EmisMask(i)%cdfname == 'NOTSET',&
            "EmisMask(i)%cdfname undefined for "//trim(EmisMask(i)%filename))
       call CheckStop(EmisMask(i)%ID == 'NOTSET',&
            "EmisMask(i)%ID undefined for "//trim(EmisMask(i)%filename))
       iEmisMask = iEmisMask+1 !assumes the fields are defined, without checking
       select case(EmisMask(i)%type)
       case('CELL-FRACTION')
          !first check if the ID has already been defined
          found=.false.
          do ii = 1, i-1
             if (EmisMask(ii)%ID==EmisMask(i)%ID) then
                iEmisMask = EmisMask_ix(ii)
                found = .true.
                exit
             end if
          end do
          if(.not. found) iEmisMask = NEmisMask+1
          EmisMask_ix(i) = iEmisMask
          NEmisMask = max(NEmisMask, iEmisMask)
       case('THRESHOLD')
          !first check if the ID has already been defined
          found=.false.
          do ii = 1, i-1
             if (EmisMask(ii)%ID==EmisMask(i)%ID) then
                iEmisMask = EmisMask_ix(ii)
                found = .true.
                exit
             end if
          end do
          if(.not. found) then
             iEmisMask = NEmisMask+1
          else
             call StopAll("Mask type THRESHOLD cannot combine same ID")
          end if
          EmisMask_ix(i) = iEmisMask
          NEmisMask = max(NEmisMask, iEmisMask)
       end select
    enddo
    if(me==0 .and. NEmisMask>NEmisMaskMAX)write(*,*)'found ',NEmisMask,' different IDs. Max is ',NEmisMaskMAX
    call CheckStop(NEmisMask>NEmisMaskMAX, "Max masks exceeded. Increase NEmisMaskMAX")
    allocate(EmisMaskValues(LIMAX,LJMAX,NEmisMask))
    EmisMaskValues = 1.0 !default is to include everything
    allocate(EmisMaskIndex2Name(NEmisMask))

    !now set the values for the actual masks
    iEmisMask = 0
    do i = 1, size(EmisMask)
       if(EmisMask(i)%filename == 'NOTSET' ) cycle
       select case(EmisMask(i)%type)
       case('NUMBER')
          call ReadField_CDF(trim(EmisMask(i)%filename),trim(EmisMask(i)%cdfname),mask_cdf,1,&
               interpol='zero_order', needed=.false., found=found, UnDef=0.0, debug_flag=.false.)
          call CheckStop(.not. found, &
               "Mask variable not found: "//trim(EmisMask(i)%cdfname)//':'//trim(EmisMask(i)%filename))
          if(.not. allocated(EmisMaskIntVal)) &
               allocate(EmisMaskIntVal(LIMAX,LJMAX))
          EmisMaskIntVal(:,:) = nint(mask_cdf(:,:))
          if(MasterProc)write(*,*)'defined mask  '//trim(EmisMask(i)%ID)//' based on '//trim(EmisMask(i)%cdfname)

       case('CELL-FRACTION') !multiplies by mask
          iEmisMask = EmisMask_ix(i)
          mask_cdf = 0.0
          call ReadField_CDF(trim(EmisMask(i)%filename),trim(EmisMask(i)%cdfname),mask_cdf,1,&
               interpol='conservative', needed=.false., found=found, UnDef=0.0, debug_flag=.false.)

          call CheckStop(.not. found, &
               "Mask variable not found: "//trim(EmisMask(i)%cdfname)//':'//trim(EmisMask(i)%filename))
          !NB: masks are additive. Include all parts that have a mask defined.
          !    May be >1 if several masks are included at same gridcell
          !    mask_cdf = 0 -> EmisMaskValues=1 -> do not change emissions
          !    mask_cdf = 1 -> EmisMaskValues=fac=0 -> remove emissions
          EmisMaskValues(:,:,iEmisMask) = EmisMaskValues(:,:,iEmisMask) - mask_cdf(:,:) * (1.0-EmisMask(i)%fac)
          ic = count(mask_cdf(:,:) > 0)
          if(MasterProc)write(*,*)'defined mask  '//trim(EmisMask(i)%ID)//' based on '//trim(EmisMask(i)%cdfname),EmisMask(i)%fac
          ! make table of names
          EmisMaskIndex2Name(iEmisMask) = trim(EmisMask(i)%ID)
          if(ic>0)write(*,*)me,' masked ',ic,' cells'

          call CheckStop(any(EmisMaskValues(:,:,iEmisMask) < 0) .or. any(EmisMaskValues(:,:,iEmisMask) > 1), &
               "Mask variable out of range: "//trim(EmisMask(i)%cdfname)//':'//trim(EmisMask(i)%filename))
       case ('THRESHOLD') !remove (possibly only 1-fac)) where mask is defined
          iEmisMask = EmisMask_ix(i)
          call ReadField_CDF(trim(EmisMask(i)%filename),trim(EmisMask(i)%cdfname),mask_cdf,1,&
               interpol='zero_order', needed=.false., found=found, UnDef=0.0, debug_flag=.false.)
          call CheckStop(.not. found, &
               "Mask variable not found: "//trim(EmisMask(i)%cdfname)//':'//trim(EmisMask(i)%filename))
          !set mask value
          ic = 0
          do jj = 1, LJMAX
             do ii = 1, LIMAX
                if(mask_cdf(ii,jj)>EmisMask(i)%threshold .and. mask_cdf(ii,jj)<EmisMask(i)%threshold_max)then
                   EmisMaskValues(ii,jj,iEmisMask) = EmisMask(i)%fac !remove everything, or fac fraction if defined
                   ic = ic + 1
                else
                   EmisMaskValues(ii,jj,iEmisMask) = 1.0 !keep everything
                end if
             end do
          end do
          if(MasterProc)write(*,*)'defined mask  '//trim(EmisMask(i)%ID)//' based on '//trim(EmisMask(i)%cdfname)
          if(ic>0)write(*,*)me,' masked ',ic,' cells'
       case default
          call StopAll("Mask type not defined: "//trim(EmisMask(i)%type)//':'//trim(EmisMask(i)%filename))
       end select
    end do

   ! to keep some compatibility with old format we also set old format masks
   if(.not. any(emis_inputlist(:)%use_mask)) &
      return
   call CheckStop(NEmisMask <= 0, &
      "No masks found for emis_inputlist(:)%use_mask")
   call CheckStop(any((EmisMask(:NEmisMask)%ID == 'NUMBER' )), &
      "emis_inputlist(:)%use_mask is incompatible with EmisMask(:)%ID='NUMBER'")
   if(.not.allocated(Emis_mask)) &
     allocate(Emis_mask(LIMAX,LJMAX))
   if(MasterProc)write(*,*)'defining mask for old formats'
   forall(jj = 1:LJMAX, ii = 1:LIMAX) &
      Emis_mask(ii,jj) = any(EmisMaskValues(ii,jj,:iEmisMask)<0.999)

 end subroutine Init_masks
!***********************************************************************
subroutine EmisUpdate
   !Update emission arrays, and read new sets as required
   integer :: n, i, j, ij, f, ix, is, is0, date_limit(5), iem, ic, icc, iqrc
   integer :: itot,isec,iland, i_femis_lonlat, ncFileID, maxfound
   type(date) :: coming_date
   real :: fac, gridyear, ccsum
   character(len=TXTLEN_NAME) :: fmt,fmtg
   TYPE(timestamp)   :: ts1,ts2
   logical, save ::first_call = .true.
   real, allocatable, dimension(:) :: emsum ! Sum of emissions per file
   real, allocatable, dimension(:,:) :: sumemis ! Sum of emissions per country
   real, allocatable, dimension(:,:,:) :: sumemis_sec ! Sum of emissions per country and sector
   real, save, allocatable, dimension(:,:) :: xsumemis ! DS testing for EMTBL output. Seems to work for Ncalls = 1, BUT BE CAREFUL
   real, allocatable, dimension(:,:,:) :: cdfemis !to send to  Emis_GetCdf
   logical :: writeoutsums
   logical :: writeout !if something to show and writeoutsums=T
   logical :: read_new_emissions
   logical :: dbgij, dbg, dbgProc, dbgX
   integer, save :: ncalls = 0
   character(len=*), parameter :: dtxt='EmUpdate:'
   character(len=5) :: emtag  ! EMSUM (and future: EMTBL, EMTAB)

   writeoutsums = first_call .or. step_main<10 .or. DEBUG%EMISSIONS
   writeout = .false. !init
   ncalls = ncalls + 1
   dbgX = DEBUG%EMISSIONS2 .and. MasterProc

   ts1=make_timestamp(current_date)
   coming_date = current_date
   coming_date%seconds = coming_date%seconds + 1800!NB: end_of_validity_date is at end of period, for example 1-1-2018 for December 2017
   gridyear = GRIDWIDTH_M * GRIDWIDTH_M * 3600*24*nydays*1.0E-6!kg/m2/s -> kt/year

   read_new_emissions = .false.
   do n = 1, NEmisFile_sources
      if(EmisFiles(n)%periodicity=='monthly') &
         writeoutsums = .true.

      if (.not. date_is_reached(to_idate(EmisFiles(n)%end_of_validity_date, 5))) &
         cycle

      if (MasterProc .and. writeoutsums) &
         write(*,'(5(A,:,X))') sub,'current date has reached past update date',&
            date2string("YYYY-MM-DD hh:mm:ss", EmisFiles(n)%end_of_validity_date),&
            trim(EmisFiles(n)%periodicity)
      read_new_emissions = .true.
   end do
   if (.not. read_new_emissions) then
      if (MasterProc .and. writeoutsums .and. DEBUG%EMISSIONS) &
         write(*,*) sub // " no new emissions to read"
      return
    end if
    if (MasterProc .and. writeoutsums) &
       write(*,*) sub // " read new and old emissions"

   if (writeoutsums) then
      writeout = .true. ! at least something to write

      allocate(emsum(NEMIS_FILE))
      emsum = 0.0

      allocate(sumemis(NLAND,NEMIS_FILE))
      sumemis = 0.0

      if(.not. allocated(xsumemis)) then
         allocate(xsumemis(NLAND,NEMIS_FILE))
         xsumemis = 0.0
      end if

      if(SecEmisTotalsWanted) then
         allocate(sumemis_sec(NLAND,NSECTORS,NEMIS_FILE))
         sumemis_sec = 0.0
      end if
   endif
   maxfound=0
   is0=0
   NEmis_source_ij=0

   !loop over all sources and see which one need to be reread from files
   do n = 1, NEmisFile_sources
      allocate(cdfemis(EmisFiles(n)%nsectors,LIMAX,LJMAX))  !NB: sector is first coordinate
      sumemis = 0.0 !sum for each file
      if(dbgX) write(*,'(a,2i6,L2)') 'XNEmisLoopA', n, EmisFiles(n)%nsectors, SecEmisTotalsWanted
      if(SecEmisTotalsWanted)sumemis_sec = 0.0
      do is = EmisFiles(n)%source_start,EmisFiles(n)%source_end
      !if(debug_proc)write(*,'(a,6i6,a,L2)') 'XNEmisBUGloopC', n, is,is0,&
      !      EmisFiles(n)%source_start,EmisFiles(n)%source_end, EmisFiles(n)%nsectors, &
      !      ' '//basename(EmisFiles(n)%filename)//':'//trim(Emis_source(is)%varname), Emis_source(is)%is3D
         if(Emis_source(is)%is3D)then
            ix = ix3Dmap(is)
            Emis_source_3D(1:,1:,1:,ix)=0.0
           !=========================================
            call Emis_GetCdf(EmisFiles(n),Emis_source(is),Emis_source_3D(1,1,1,ix),coming_date)
           !=========================================
         else
            !if(debug_proc)write(*,'(a,6i6)') 'XNEmisBUGloopD', n, is0, is,&
            !        is-EmisFiles(n)%source_start, EmisFiles(n)%nsectors, &
            !          mod(is-EmisFiles(n)%source_start, EmisFiles(n)%nsectors)
            if (mod(is-EmisFiles(n)%source_start, EmisFiles(n)%nsectors) == 0) then
               is0=is-1 !shift to get at start of cdfemis array
               !we read EmisFiles(n)%nsectors at once
               !if only one sector per variable, EmisFiles(n)%nsectors=1 and, we only rewrite one 2D field
               cdfemis=0.0
               !=========================================
               call Emis_GetCdf(EmisFiles(n),Emis_source(is),cdfemis,coming_date)
               !=========================================
            !   if(debug_proc.and.maxval(cdfemis)>1.0e-10) write(*,'(a,3i6,es12.3,a)') 'XNEmisBUGloopGet', n, is0, is, maxval(cdfemis), &
            !           trims(Emis_source(is)%varname//':' //Emis_source(is)%species)

               if(me==0 .and. DEBUG%EMISSIONS)write(*,'(i6,a,es12.3)') is,&
                  dtxt//' getemis '//trim(Emis_source(is)%units)//' '//&
                     trim(Emis_source(is)%varname)//' read as '//trim(Emis_source(is)%species), maxval(cdfemis)

            end if
         end if
         !reduction factors
         fac = EmisFiles(n)%factor
         fac = fac* Emis_source(is)%factor

         !unit and factor conversions
         !convert into kg/m2/s
         select case(Emis_source(is)%units)
         case ('kg/s', 'kg/m2/s', 'kg m-2 s-1')
            fac = fac
         case ('g/s', 'g/m2/s', 'g m-2 s-1')
            fac = fac /(1000.0)
         case ('mg/s', 'mg/m2/s', 'mg m-2 s-1')
            fac = fac /(1000*1000.0)
         case default
            !if not /s units, depends on periodicity:
            if(EmisFiles(n)%periodicity == 'yearly')then
               select case(Emis_source(is)%units)
               case ('kt/m2', 'kt/m2/year', 'kt', 'kt/year','kt m-2', 'kt m-2 year-1', 'kt year-1')
                     fac = fac * 1000 *1000/(3600*24*nydays)
               case ('tonnes/m2', 'tonnes/m2/year', 'tonnes', 'tonnes/year','tonnes m-2', 'tonnes m-2 year-1', 'tonnes year-1', 'Mg' )
                  fac = fac * 1000/(3600*24*nydays)
               case ('kg/m2', 'kg m-2', 'kg/m2/year', 'kg m-2 year-1', 'kg', 'kg/year', 'kg year-1')
                  fac = fac /(3600*24*nydays)
               case ('g/m2', 'g/m2/year', 'g', 'g/year', 'g m-2', 'g m-2 year-1', 'g year-1')
                  fac = fac /(1000.0*3600*24*nydays)
               case ('mg/m2', 'mg/m2/year', 'mg', 'mg/year', 'mg m-2', 'mg m-2 year-1', 'mg year-1')
                  fac = fac /(1.0e6*3600*24*nydays)
               case default
                  call StopAll("Unit for yearly emissions not recognized: "//trim(Emis_source(is)%units))
               end select
            else if(EmisFiles(n)%periodicity == 'monthly')then
               select case(Emis_source(is)%units)
               case ('kt/m2','kt/m2/month', 'kt', 'kt/month','kt m-2','kt m-2 month-1', 'kt month-1')
                  fac = fac *1000*1000/(3600*24*nmdays(coming_date%month))
               case ('tonnes/m2', 'tonnes/m2/month', 'tonnes', 'tonnes/month','tonnes m-2', 'tonnes m-2 month-1', 'tonnes month-1')
                  fac = fac *1000/(3600*24*nmdays(coming_date%month))
               case ('kg/m2', 'kg/m2/month', 'kg', 'kg/month','kg m-2', 'kg m-2 month-1', 'kg month-1')
                  fac = fac /(3600*24*nmdays(coming_date%month))
               case ('g/m2', 'g/m2/month', 'g', 'g/month','g m-2', 'g m-2 month-1', 'g month-1')
                  fac = fac /(1000*3600*24*nmdays(coming_date%month))
               case ('mg/m2', 'mg/m2/month', 'mg', 'mg/month','mg m-2', 'mg m-2 month-1', 'mg month-1')
                  fac = fac /(1.0e6*3600*24*nmdays(coming_date%month))
               case default
                  call StopAll("Unit for montly emissions not recognized: "//trim(Emis_source(is)%units))
               end select
            else
               !assume hourly
               select case(Emis_source(is)%units)
               case ('mg/m2', 'mg/m2/h', 'mg m-2', 'mg m-2 h-1')
                  !convert into kg/m2/s
                  fac = fac /(1.0e6*3600.0)
               case ('g/m2', 'g/m2/h','g m-2', 'g m-2 h-1')
                  fac = fac /(1000.0*3600.0)
                  if(EmisFiles(n)%periodicity /= 'hourly' .and. Emis_source(is)%units == 'g/m2')then
                     call StopAll("Emis_source unit g/m2 only implemented for hourly, monthly or yearly. Found "//&
                             trim(EmisFiles(n)%periodicity))
                  endif
               case ('g/s','g s-1')
                  fac = fac /(1000.0)
               case ('kg/s','kg s-1')
                  fac = fac
               case ('tonnes/s', 'tonnes s-1')
                  fac = fac * 1000.0
               case default
                  call StopAll("Emis_source unit not implemented. Found "//&
                          trim(Emis_source(is)%units)//' '//trim(EmisFiles(n)%periodicity))
               end select
               !Note: easy to implement more unit choices. Just add appropriate cases here
            endif
         end select

         ! units that are "per gridcell" must be dvided by gridcell areas
         select case(Emis_source(is)%units)
         case ('kt', 'Mg', 'kt/s', 'kt/month', 'kt/year', 'tonnes', 'tonnes/s', 'tonnes/month' &
            , 'tonnes/year', 'kg', 'kg/s', 'kg/month', 'kg/year', 'g', 'g/s', 'g/month' &
            , 'g/year', 'mg', 'mg/s', 'mg/month', 'mg/year', 'g/h', 'mg/h'&
            , 'kt s-1', 'kt month-1', 'kt year-1', 'tonnes s-1', 'tonnes month-1' &
            , 'tonnes year-1', 'kg s-1', 'kg month-1', 'kg year-1', 'g s-1', 'g month-1' &
            , 'g year-1', 'mg s-1', 'mg month-1', 'mg year-1', 'g h-1', 'mg h-1')
            !divide by gridarea
            fac = fac / (GRIDWIDTH_M * GRIDWIDTH_M)
            do j = 1,ljmax
               do i = 1,limax
                  cdfemis(is-is0,i,j) = cdfemis(is-is0,i,j) * xm2(i,j)
               end do
            end do
         end select

         do j = 1,ljmax
            do i = 1,limax
               cdfemis(is-is0,i,j) = cdfemis(is-is0,i,j) * fac
            end do
         end do

         !now cdfemis should be in kg/m2/s
         !apply masks
         !could easily be better CPU optimized if necessary by putting all factors in same i,j loop
         if(Emis_source(is)%mask_ix>0)then
            do j = 1,ljmax
               do i = 1,limax
                  cdfemis(is-is0,i,j) = cdfemis(is-is0,i,j) * EmisMaskValues(i,j,Emis_source(is)%mask_ix)
               end do
            end do
         end if
         if(Emis_source(is)%mask_reverse_ix>0)then
            do j = 1,ljmax
               do i = 1,limax
                  if (abs(EmisMaskValues(i,j,Emis_source(is)%mask_ix)-1.0)>1.0E-10) then
                     !mutiply by fac in the not covered region
                     cdfemis(is-is0,i,j) = cdfemis(is-is0,i,j) * EmisMaskValues(i,j,Emis_source(is)%mask_ix)
                  else
                     !keep unchanged otherwise
                  end if
               end do
            end do
         end if
         if(Emis_source(is)%apply_femis)then
            !apply femis_lonlat
            isec = Emis_source(is)%sector
            if(Emis_source(is)%sector>0 .and. Emis_source(is)%sector<=NSECTORS)then
               iland = Emis_source(is)%country_ix
               if(iland<0)iland=IC_DUMMY
               isec = Emis_source(is)%sector
               iem = find_index(Emis_source(is)%species,EMIS_FILE(:))
               if(iem >0 )then
                  if(N_femis_lonlat>0 .and. isec>0) then
                     do i_femis_lonlat=1,N_femis_lonlat
                        if(femis_lonlat_ic(i_femis_lonlat)==0 .or. &
                           femis_lonlat_ic(i_femis_lonlat)==Country(Emis_source(is)%country_ix)%icode)then
                           do j=1,ljmax
                              do i=1,limax
                                 if(glat(i,j)>femis_latmin(i_femis_lonlat).and.&
                                    glat(i,j)<femis_latmax(i_femis_lonlat).and.&
                                    glon(i,j)>femis_lonmin(i_femis_lonlat).and.&
                                    glon(i,j)<femis_lonmax(i_femis_lonlat))then
                                    cdfemis(is-is0,i,j) = cdfemis(is-is0,i,j) * e_fact_lonlat(isec,i_femis_lonlat,iem)
                                 endif
                              end do
                           end do
                        end if
                     end do
                  end if
               else
                  !see if the species belongs to any of the splitted species
                  iqrc = 0
                  do iem = 1,NEMIS_FILE
                     do f = 1,emis_nsplit(iem)
                        iqrc = iqrc + 1
                        itot = iqrc2itot(iqrc)
                        if(trim(species(itot)%name)==trim(Emis_source(is)%species))then
                           !lonlat format
                           if(N_femis_lonlat>0 .and. isec>0) then
                              do i_femis_lonlat=1,N_femis_lonlat
                                 if(femis_lonlat_ic(i_femis_lonlat)==0 .or. &
                                    femis_lonlat_ic(i_femis_lonlat)==Country(Emis_source(is)%country_ix)%icode)then
                                    do j=1,ljmax
                                       do i=1,limax
                                          if(glat(i,j)>femis_latmin(i_femis_lonlat).and.&
                                             glat(i,j)<femis_latmax(i_femis_lonlat).and.&
                                             glon(i,j)>femis_lonmin(i_femis_lonlat).and.&
                                             glon(i,j)<femis_lonmax(i_femis_lonlat))then
                                             cdfemis(is-is0,i,j) = cdfemis(is-is0,i,j) * e_fact_lonlat(isec,i_femis_lonlat,iem)
                                       endif
                                       end do
                                    end do
                                 end if
                              end do
                           end if
                           go to 888
                        end if
                     end do
                  end do ! iem
888               continue
               end if
            end if
         end if

         !DSTMP writeout=.false.
         if(dbgX)write(*,*) 'XNEmisloopW', is, is0
         if(writeout)then
            ! sum emissions per countries (in ktonnes?)
            itot = Emis_source(is)%species_ix
            isec = Emis_source(is)%sector
            iland = Emis_source(is)%country_ix
            !FOUND BUGif(is>19) print *, 'XNEmisBUGwrite', is, is0, isec, itot, iland
            if(itot>0)then
               iqrc = itot2iqrc(itot)
               if(isec>0 .and. iqrc>0)then
                  iem = iqrc2iem(iqrc)
               else
                  iem=-1
               endif
            else
               iem=find_index(Emis_source(is)%species,EMIS_FILE(:))
            end if
            if(iem>0)then
               do j = 1,ljmax
                  do i = 1,limax
                     sumemis(iland,iem) = sumemis(iland,iem) + cdfemis(is-is0,i,j) * gridyear * xmd(i,j) !now in kt/year
                     xsumemis(iland,iem) = xsumemis(iland,iem) + cdfemis(is-is0,i,j) * gridyear * xmd(i,j) !now in kt/year
                     if(SecEmisTotalsWanted)&
                        sumemis_sec(iland,isec,iem) = sumemis_sec(iland,isec,iem)&
                        + cdfemis(is-is0,i,j) * gridyear * xmd(i,j)
                        !DSBUG was: ??! sumemis_sec(iland,is,iem) = sumemis_sec(iland,is,iem)&
                  end do
               end do
            end if
         end if
         if (is-is0 == EmisFiles(n)%nsectors) then
            ! now we store the emission fields in a more memory effective form,
            ! since most country emissions are zero outside of the country.
            !note that we store only once all the EmisFiles(n)%nsectors are
            ij = 0
            do j=1,ljmax
               do i=1,limax
                  ij=ij+1
                  dbgij= (debug_proc .and. i==debug_li .and. j==debug_lj)
                  do isec=1,EmisFiles(n)%nsectors
                     if (cdfemis(isec,i,j)>1e-20) then
                        if(NEmis_source_ij(ij)>NEmis_source_ijMAX) &
                                write(*,*)dtxt//'>ijmax', me,i,j,NEmis_source_ij(ij)
                        call CheckStop(NEmis_source_ij(ij)>NEmis_source_ijMAX, &
                                dtxt//'NEmis_source_ijMAX too small. Please increase it')
                        NEmis_source_ij(ij)=NEmis_source_ij(ij)+1
                        maxfound= max(maxfound,NEmis_source_ij(ij))
                        Emis_source_ij_ix(ij,NEmis_source_ij(ij))=is0+isec
                        Emis_source_ij(ij,NEmis_source_ij(ij)) =cdfemis(isec,i,j)
                        if (dbgij) write(*,'(a,4i6,es12.3)') dtxt//'BUGdbgIJ',&
                                isec, NEmis_source_ij(ij), is0, isec, cdfemis(isec,i,j) 
                     end if
                  end do ! isec
               end do !i
            end do !j
         end if !is-is0
      end do
      deallocate(cdfemis)

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
         end if
      else
         !the correct times must be written in the file and updated in Emis_GetCdf
      end if
      if(EmisFiles(n)%ncFileID >= 0) then
         call check(nf90_close(EmisFiles(n)%ncFileID))
         EmisFiles(n)%ncFileID = -1 !mark as closed
      end if

      if (writeout) then
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,sumemis,&
            NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,MPI_COMM_CALC,IERROR)
         if(SecEmisTotalsWanted)CALL MPI_ALLREDUCE(MPI_IN_PLACE,sumemis_sec,&
            NLAND*NSECTORS*NEMIS_FILE,MPI_REAL8,MPI_SUM,MPI_COMM_CALC,IERROR)
         if(me==0)then
            !EMTABDS write(*,*)"Emissions per country for "//trim(EmisFiles(n)%filename)//' (Gg/year) '
            if (NEmisFile_sources>1) then !we do not write out separately per file if only one file
               call PrintLog("#EMTBL Total emissions by countries for "//trim(EmisFiles(n)%filename)//' (Gg/year) ')
               write(*     ,"(a14,a5,3x,30(a12,:))")"EMTBL CC Land ","    ",EMIS_FILE(:)
               write(IO_LOG,"(a14,a5,3x,30(a12,:))")"EMTBL CC Land ","    ",EMIS_FILE(:)
               fmt="(a5,i4,1x,a9,3x,30(f12.2,:))"
               fmtg="(a5,i4,1x,a9,3x,30(g12.2,:))"
               do ic = 1, NLAND
                  ccsum = sum( sumemis(ic,:) )
                  icc = Country(ic)%icode
                  if (ccsum > 0.0) then
                     if (ccsum > 1e8 )then
                        write(*,     fmtg) 'EMTBL', icc, Country(ic)%code, sumemis(ic,:)
                        write(IO_LOG,fmtg) 'EMTBL', icc, Country(ic)%code, sumemis(ic,:)
                     else
                        write(*,     fmt) 'EMTBL', icc, Country(ic)%code, sumemis(ic,:)
                        write(IO_LOG,fmt) 'EMTBL', icc, Country(ic)%code, sumemis(ic,:)
                     end if
                  end if
               end do
            end if
            if(SecEmisTotalsWanted)then
               write(*,"(a)") " SECTOR EMISSIONS SUMMARY ==================== "
               write(*,"(a19,2x,30(a12,:))")"  CC Land    Sector",EMIS_FILE(:)
               fmt="(i4,1x,a9,i4,3x,30(f12.2,:))"
               do ic = 1, NLAND
                  do is = 1, NSECTORS
                     ccsum = sum( sumemis_sec(ic,is,:) )
                     icc=Country(ic)%icode
                     if ( ccsum > 0.0 )then
                       if (ccsum > 1e8 )then
                         write(*, fmtg) icc, Country(ic)%code, is, sumemis_sec(ic,is,:)
                       else
                         write(*, fmt) icc, Country(ic)%code, is, sumemis_sec(ic,is,:)
                       end if
                     end if
                  end do
               end do
            end if
         end if
         !total of emissions from all countries and files into emsum
         do iem = 1, NEMIS_FILE
            emsum(iem)= emsum(iem)+sum(sumemis(:,iem))
         end do

      endif
   enddo
   if(writeout)then
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxfound,1,MPI_INTEGER,MPI_MAX,MPI_COMM_CALC,IERROR)
      if(me==0)write(*,*)'max sources found in single gridcell ',maxfound
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,xsumemis,&
            NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,MPI_COMM_CALC,IERROR)
      fmt="(a5,i4,1x,a9,3x,30(f12.2,:))" ! reset
      if(me==0 .and. NEmisFile_sources>0) then
         !call PrintLog('#EMTBL Total emissions by countries, all Emis_source files (Gg/year) ')
         ! Use EMSUM for  final table. Much easier to find in olog and RunLog!
         emtag = 'EMSUM'
         call PrintLog('#'//emtag//' Total emissions by countries, all Emis_source files (Gg/year) ')
         write(*     ,"(a23,a5,3x,30(a12,:))") emtag//"    NCalls CC Land","    ",EMIS_FILE(:)
         write(IO_LOG,"(a23,a5,3x,30(a12,:))") emtag//"    NCalls CC Land ","    ",EMIS_FILE(:)
         fmtg="(a5,i9,i4,1x,a9,3x,30(g12.2,:))"
         fmt="(a5,i9,i4,1x,a9,3x,30(f12.2,:))" !nicer but limited max emissions
         do ic = 1, NLAND
            ccsum = sum( xsumemis(ic,:) )
            icc=Country(ic)%icode
            if ( ccsum > 0.0 )then
               if (ccsum > 1e8 )then
                  write(*,     fmtg) emtag , ncalls, icc, Country(ic)%code, xsumemis(ic,:)
                  write(IO_LOG,fmtg) emtag , ncalls, icc, Country(ic)%code, xsumemis(ic,:)
               else
                  write(*,     fmt) emtag, ncalls, icc, Country(ic)%code, xsumemis(ic,:)
                  write(IO_LOG,fmt) emtag, ncalls, icc, Country(ic)%code, xsumemis(ic,:)
               end if
            end if
         end do
         ccsum = sum( emsum(:))
         if (ccsum > 1e8 )then
            write(*     ,fmtg)emtag,  ncalls,999,'TOTAL    ',emsum(:)
            write(IO_LOG,fmtg)emtag,  ncalls,999,'TOTAL    ',emsum(:)
         else
            write(*     ,fmt)emtag,  ncalls,999,'TOTAL    ',emsum(:)
            write(IO_LOG,fmt)emtag,  ncalls,999,'TOTAL    ',emsum(:)
         end if
      end if
      fmt="(a5,i4,1x,a9,3x,30(f12.2,:))" ! reset
      CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)!so that print out comes out nicely
   endif
   if (allocated(emsum)) deallocate(emsum)
   if (allocated(sumemis)) deallocate(sumemis)
   if (allocated(sumemis_sec)) deallocate(sumemis_sec)

   first_call=.false.

end subroutine EmisUpdate

!***********************************************************************
  subroutine Emissions(year)
    ! Initialize emission variables, and read yearly emissions
    integer, intent(in)   :: year        ! Year ( 4-digit)

    !-- local variables
    integer :: i, j, is          ! Loop variables
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
    real, dimension(NLAND,NEMIS_FILE) :: sumemis ! Sum of emissions per country
    !note: automatic arrays were not accepted for MPI_REDUCE
    real, allocatable, dimension(:,:,:) :: sumemis_sec, sumemis_local ! Sum of emissions per country and sectors
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
    integer largestsplit !used only here to check value
    integer :: f, itot, iqrc
    character(len=*), parameter :: dtxt='YEmis:'

    if (MasterProc) write(6,*) "Initializing emissions for year",  year

    if (MasterProc) write(6,*) "Reading emissions for year",  year

    ios = 0
    country_owner_map = .false.


    ! initialize emis_inputlist
    !>=========================================================
    do iemislist = 1, size( emis_inputlist(:)%name )
       fname = emis_inputlist(iemislist)%name
       if(fname=="NOTSET") cycle
       call StopAll("emis_inputlist style input no more in use. Use Emis_sourceFiles instead ")
       if(MasterProc)&
            write(*,*)"Emission source number ", iemislist,"from ",sub//trim(fname)

       if(emis_inputlist(iemislist)%type == "sectors".or.&
            emis_inputlist(iemislist)%type == "GNFR_CAMSsectors".or.&
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

       ! find which sectors types are being used
       cdf_sector_name='NOTSET'
       fname = key2str(fname,'POLL',EMIS_FILE(1))
       call ReadSectorname(fname,cdf_sector_name)

       if (emis_inputlist(iemislist)%sector /= "NOTSET") then
          ! we will use this even if something else is defined in the file
          if (emis_inputlist(iemislist)%sector == 'GNFR') then
             if(MasterProc)write(*,*)'Redefined GNFR sector as GNFR_CAMS sector'
             emis_inputlist(iemislist)%sector = 'GNFR_CAMS'
          end if
       else if (emis_inputlist(iemislist)%type == "GNFR_CAMSsectors") then
          emis_inputlist(iemislist)%sector = 'GNFR_CAMS'
       else
          ! the sector type must be read from the file
          if (cdf_sector_name /=  'NOTSET') then
             if (cdf_sector_name == 'GNFR') cdf_sector_name = 'GNFR_CAMS' ! GNFR is a subset of GNFR_CAMS
             if (cdf_sector_name == 'GNFR_CAMS') emis_inputlist(iemislist)%type = "GNFR_CAMSsectors"
             if (cdf_sector_name == 'SNAP') emis_inputlist(iemislist)%type = "SNAPsectors"
             emis_inputlist(iemislist)%sector = trim(cdf_sector_name)
          else if (emis_inputlist(iemislist)%type=="sectors") then
             call CheckStop(emis_inputlist(iemislist)%type=="sectors", "Did not find sector name in config or in file"//trim(emis_inputlist(iemislist)%name))
          end if
       end if

       if (MasterProc) then
          if (trim(cdf_sector_name) /= emis_inputlist(iemislist)%sector) then
             write(*,*)'WARNING: using '//trim(emis_inputlist(iemislist)%sector)//' sectors even if '//trim(cdf_sector_name)//' is defined in '//trim(emis_inputlist(iemislist)%name)
          else
             if(emis_inputlist(iemislist)%sector /= "NOTSET") then
                write(*,*)trim(emis_inputlist(iemislist)%sector)//' sectors defined for '//trim(emis_inputlist(iemislist)%name)
             else
                write(*,*)trim(emis_inputlist(iemislist)%name)//' not considered as sectors type emissions '
             end if
          end if
       end if

       !Add the relevant sector in SECTORS
       !look if it is already defined
       found = 0
       if(emis_inputlist(iemislist)%sector == "NOTSET") found = 1 !not a sector type emissions (i.e. DMS etc.)
       do i = 1, NSECTORS
          if(SECTORS(i)%name ==  trim(emis_inputlist(iemislist)%sector)) found = 1
       end do
       if (found == 0) then
          !add this sector in SECTORS
          do i = 1, NSECTORS_ADD_MAX
             if(SECTORS_ADD(i)%name ==  trim(emis_inputlist(iemislist)%sector)) then
                NSECTORS = NSECTORS + 1
                if(MasterProc) write(*,*) dtxt//'adding sector2nd  ',&
                         trim(SECTORS_ADD(i)%longname), i
                call CheckStop(NSECTORS > NSECTORS_MAX, &
                    "NSECTORS_MAX too small, please increase value in EmisDef_mod.f90")
                SECTORS(NSECTORS) = SECTORS_ADD(i)
                found = 1
             end if
          end do
          if (found == 0 .and. (trim(emis_inputlist(iemislist)%sector) == 'GNFR_CAMS' &
               .or. trim(emis_inputlist(iemislist)%sector) == 'GNFR')) then ! note that default values are used only if not defined in config
             do i = 1, NSECTORS_GNFR_CAMS
                NSECTORS = NSECTORS + 1
                call CheckStop(NSECTORS > NSECTORS_MAX, "NSECTORS_MAX too small, please increase value in EmisDef_mod.f90")
                SECTORS(NSECTORS) = GNFR_CAMS_SECTORS(i)
             end do
             if(MasterProc) write(*,*)'including GNFR_CAMS sectors ', NSECTORS
             found = 1
          end if
          if (found == 0 .and. trim(emis_inputlist(iemislist)%sector) == 'SNAP') then ! note that default values are used only if not defined in config
             do i = 1, NSECTORS_SNAP
                NSECTORS = NSECTORS + 1
                call CheckStop(NSECTORS > NSECTORS_MAX, "NSECTORS_MAX too small, please increase value in EmisDef_mod.f90")
                SECTORS(NSECTORS) = SNAP_SECTORS(i)
             end do
             found = 1
          end if
       end if

       call CheckStop(found == 0, 'Sector name not recognized: '//trim(emis_inputlist(iemislist)%sector))

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
          USES%OCEAN_NH3=.true.
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
          USES%OCEAN_DMS=.true.
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
    if (NSECTORS==0) then
       !the code does not work if no sectors are defined. TODO: find why (zero size arrays? no TFAC_IDX_XXX?)
       !first see if some are defined by config
       do i = 1, NSECTORS_ADD_MAX
          if(SECTORS_ADD(i)%name /=  'NOTSET') then
             NSECTORS = NSECTORS + 1
             call CheckStop(NSECTORS > NSECTORS_MAX, "NSECTORS_MAX too small, please increase value in EmisDef_mod.f90")
             SECTORS(NSECTORS) = SECTORS_ADD(i)
             found = 1
          end if
       end do
       if (NSECTORS==0) then
          do i = 1, NSECTORS_GNFR_CAMS
             NSECTORS = NSECTORS + 1
             call CheckStop(NSECTORS > NSECTORS_MAX, "NSECTORS_MAX too small, please increase value in EmisDef_mod.f90")
             SECTORS(NSECTORS) = GNFR_CAMS_SECTORS(i)
          end do
          if(MasterProc) write(*,*)'including default GNFR_CAMS sectors ', NSECTORS
       end if
    end if

    i=0
    N_TFAC = 0
    N_HFAC = 0
    largestsplit = 0
    if (MasterProc) write(*,*)dtxt, NSECTORS,' sectors defined in total'
       if (MasterProc) write(*,*)"   name,   longname, cdfname, timefac,height,split, description"
    do isec = 1, NSECTORS
       91 format(A10,A10,A10,3I7,A)
       if (MasterProc) write(*,91)trim(SECTORS(isec)%name), &
          trim(SECTORS(isec)%longname),trim(SECTORS(isec)%cdfname),&
          SECTORS(isec)%timefac,SECTORS(isec)%height,SECTORS(isec)%split,&
          ' '//trim(SECTORS(isec)%description)
    end do
    ! find indices for special sectors
    do isec = NSECTORS, 1, -1 !reverse order so that TFAC_IDX get value from the lowest (= most "fundamental")
       if (SECTORS(isec)%longname(1:6) == 'GNFR_A' .or. trim(SECTORS(isec)%longname) == 'SNAP1') then
          IS_POW(isec) = .true.
          TFAC_IDX_POW = SECTORS(isec)%timefac
       end if
       if (SECTORS(isec)%longname(1:6) == 'GNFR_C' .or. trim(SECTORS(isec)%longname) == 'SNAP2') then
          IS_DOM(isec) = .true.
          TFAC_IDX_DOM = SECTORS(isec)%timefac
       end if
       if (SECTORS(isec)%longname(1:6) == 'GNFR_F' .or. trim(SECTORS(isec)%longname) == 'SNAP7') then
          IS_TRAF(isec) = .true.
          TFAC_IDX_TRAF = SECTORS(isec)%timefac
       end if
       if (SECTORS(isec)%longname(1:6) == 'GNFR_K' .or. trim(SECTORS(isec)%longname) == 'SNAP10') then
           IS_AGR(isec) = .true.
           TFAC_IDX_AGR = SECTORS(isec)%timefac
       end if
       if (SECTORS(isec)%longname(1:6) == 'GNFR_B' .or. SECTORS(isec)%longname(1:6) == 'GNFR_D' .or.&
            trim(SECTORS(isec)%longname) == 'SNAP3' .or. trim(SECTORS(isec)%longname) == 'SNAP4') then
          IS_IND(isec) = .true.
       end if
    end do

    do isec = 1, NSECTORS
       if (SECTORS(isec)%timefac ==  TFAC_IDX_DOM .or. IS_DOM(isec)) then
          IS_DOM(isec) = .true.
          i=1
          if (MasterProc) write(*,'(a)',advance='no') dtxt//&
             trim(SECTORS(isec)%longname)//" , "//trim(SECTORS(isec)%description)
          if(USES%DEGREEDAY_FACTORS) then
             if (MasterProc) write(*,*)': used with the Degree-Days method'
          else
             if (MasterProc) write(*,*)': recognized as a domestic sector'
          end if
       end if
       if (SECTORS(isec)%timefac ==  TFAC_IDX_AGR .or. IS_AGR(isec)) then
          IS_AGR(isec) = .true.
          if (MasterProc) then
             write(*,'(a)')dtxt//trim(SECTORS(isec)%longname)//" , "//&
                trim(SECTORS(isec)%description) // &
                ': recognized as an agricultur sector'
          end if
       end if
       if (SECTORS(isec)%timefac ==  TFAC_IDX_TRAF .or. IS_TRAF(isec)) then
          IS_TRAF(isec) = .true.
          if (MasterProc) write(*,*)dtxt//trim(SECTORS(isec)%longname)//" , "//&
             trim(SECTORS(isec)%description)// &
             ': recognized as a traffic sector'
       end if
       if (SECTORS(isec)%timefac ==  TFAC_IDX_POW .or. IS_POW(isec)) then
          IS_POW(isec) = .true.
          if (MasterProc) write(*,*)dtxt//trim(SECTORS(isec)%longname)//" , "//&
            trim(SECTORS(isec)%description)// &
            ': recognized as a Public Power sector'
       end if

       N_TFAC  = max(N_TFAC, SECTORS(isec)%timefac)
       N_HFAC  = max(N_HFAC, SECTORS(isec)%height)
       largestsplit = max(largestsplit, SECTORS(isec)%split)
    end do
    call CheckStop(USES%DEGREEDAY_FACTORS .and. i==0,&
            " did not find any sector corresponding to domestic")

    if(Emis_mask_allocate)then
       if(.not.allocated(Emis_mask))then
          allocate(Emis_mask(LIMAX,LJMAX))
          Emis_mask = .false.
       else
          !the masks has been set by init_mask
       endif
    endif

    ! init_sectors

    call CheckStop(NSECTORS>MaxNSECTORS, "redefine larger MaxNSECTORS")

    allocate(Emis_field(LIMAX,LJMAX,10))
    NEmis_id = 0

    allocate(cdfemis(LIMAX,LJMAX))
    allocate(nGridEmisCodes(LIMAX,LJMAX))
    allocate(GridEmisCodes(LIMAX,LJMAX,NCMAX))
    allocate(GridEmis(NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILE),stat=allocerr)
    allocate(sumemis_sec(NLAND,NSECTORS,NEMIS_FILE), sumemis_local(NLAND,NSECTORS,NEMIS_FILE))
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
    if(.not.allocated(fac_ehh24x7))allocate(fac_ehh24x7(NEMIS_FILE,N_TFAC,24,7,NLAND))
    if(.not.allocated(fac_emm))allocate(fac_emm(NLAND,12,N_TFAC,NEMIS_FILE))
    if(.not.allocated(fac_min))allocate(fac_min(NLAND,N_TFAC,NEMIS_FILE))
    if(.not.allocated(fac_edd))allocate(fac_edd(NLAND, 7,N_TFAC,NEMIS_FILE))
    if(timeFacs%Day_of_Year .and. .not.allocated(fac_dayofyear)) &
         allocate(fac_dayofyear(NSECTORS,NLAND,NEMIS_FILE,366))

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
       write(*,*) dtxt//"Reading monthly and daily timefactors"
       if(timeFacs%Monthly == 'GRIDDED') then !_EMIS_MONTHLY_FACTOR)then
          write(*,*)"Emissions using gridded monhtly timefactors "
          write(IO_LOG,*)"Emissions using gridded monhtly timefactors "
       end if
       !=========================
       call timefactors(year)               ! => fac_emm, fac_edd, changes ios
       !QUERY call CheckStop(ios, "ioserror: after timefacs")
       !=========================
    end if
    !=========================
    !S24: MOVED to top call EmisSplit()    ! In EmisGet_mod, => emisfrac
   !if defined own sectors, assumes that we know what we are doing:
    call CheckStop(N_SPLIT == 19 .and. largestsplit<=11 .and. &
         (SECTORS(1)%name=='GNFR_CAMS'.or.SECTORS(1)%name=='SNAP'),&
    " Cannot use GNFR_CAMS splits together with SNAP sector split settings")
    call CheckStop(largestsplit>N_SPLIT," split index not defined ")
    !=========================
    !Must first call EmisSplit, to get nrcemis defined
    if(EmisSplit_OUT)then
       allocate(SplitEmisOut(LIMAX,LJMAX,nrcemis))
       SplitEmisOut=0.0
    end if
    !=========================
    !QUERY WHY HERE? call CheckStop(ios, "ioserror: EmisSplit")

    ! ####################################
    ! Broadcast  monthly and Daily factors (and hourly factors if needed/wanted)
    CALL MPI_BCAST(fac_cemm,8*12,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(fac_emm,8*NLAND*12*N_TFAC*NEMIS_FILE,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(fac_edd,8*NLAND*7*N_TFAC*NEMIS_FILE,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(fac_ehh24x7,8*NEMIS_FILE*N_TFAC*24*7*NLAND,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    if(timeFacs%Day_of_Year ) CALL MPI_BCAST(fac_dayofyear,&
            8*NEMIS_FILE*NSECTORS*366*NLAND,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    !define fac_min for all processors
    forall(iemis=1:NEMIS_FILE,insec=1:N_TFAC,inland=1:NLAND) &
         fac_min(inland,insec,iemis) = minval(fac_emm(inland,:,insec,iemis))
    if(INERIS_SNAP2) then
       !  INERIS do not use any base-line for SNAP2
       do isec = 1, NSECTORS
          if (IS_DOM(isec)) fac_min(:,isec,:) = 0.
       end do
     end if

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
       sumemis_sec=0.0
       sumemis_local=0.0

       nin = emis_inputlist(iemislist)%Nincl
       nex = emis_inputlist(iemislist)%Nexcl

       if((emis_inputlist(iemislist)%type == "sectors".or.&
            emis_inputlist(iemislist)%type == "GNFRsectors".or.&
            emis_inputlist(iemislist)%type == "GNFR_CAMSsectors".or.&
            emis_inputlist(iemislist)%type == "SNAPsectors") .and. &
            index(emis_inputlist(iemislist)%name,".nc")>1)then

          foundYearlySectorEmissions = .true.
          emis_inputlist_NEMIS_FILE = 1!all pollutants are in same file
          if(index(emis_inputlist(iemislist)%name,"POLL")>0) &
                  emis_inputlist_NEMIS_FILE = NEMIS_FILE !one file per pollutant
          do iem = 1, emis_inputlist_NEMIS_FILE
             if(index(emis_inputlist(iemislist)%name,"POLL")>0)then
                !in this case we apply PollName restrictions here, otherwise
                ! it is applied in EmisGetCdf
                if(emis_inputlist(iemislist)%pollName(1)/='NOTSET')then
                   if(all(emis_inputlist(iemislist)%pollName(:)/=trim(EMIS_FILE(iem))))cycle
                   if(Masterproc)write(*,"(A)")dtxt//'using PollNames restrictions '
                end if
             end if
             fname = key2str(trim(emis_inputlist(iemislist)%name),'POLL',EMIS_FILE(iem))

             if(MasterProc) write(*,*) dtxt//'TRACK call EmisGetCdf:'//trim(fname)
             call EmisGetCdf(iem, fname, sumemis_local, &
                  GridEmis, GridEmisCodes, nGridEmisCodes, 1,&
                  emis_inputlist(iemislist)%incl, nin, emis_inputlist(iemislist)%excl, nex, &
                  emis_inputlist(iemislist)%use_lonlat_femis,&
                  emis_inputlist(iemislist)%set_mask,emis_inputlist(iemislist)%use_mask,&
                  emis_inputlist(iemislist)%pollName,&
                  fractionformat,emis_inputlist(iemislist)%sector)
             if ( debug_proc) write(*,*) dtxt//'XNEmisA:'//basename(fname)!//':'&
                     !//trim(emis_inputlist(iemislist)%pollName),maxval(GridEmis)

          end do!NEMIS_FILE

          !add together totals from each processor (only me=0 get results)
          sumemis_sec=0.0
          CALL MPI_REDUCE(sumemis_local,sumemis_sec,&
               NLAND*NSECTORS*NEMIS_FILE,MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)

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
          sumemis_sec=0.0
          CALL MPI_REDUCE(sumemis_local,sumemis_sec,&
               NLAND*NSECTORS*NEMIS_FILE,MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)
       else
          if(MasterProc)write(*,*)'WARNING: did not recognize format of '//trim(emis_inputlist(iemislist)%name)
          call StopAll("Emissions file format not recognized ")
       end if


       if(MasterProc) write(*,*) dtxt//'XNEmis check emislist)', iemislist, &
               trim(emis_inputlist(iemislist)%periodicity)
       if(MasterProc.and. emis_inputlist(iemislist)%periodicity == "once") then
          ! Added EMTAB to make parsing easy. These data are important!
          call PrintLog("#EMTAB Total emissions by countries for "//trim(emis_inputlist(iemislist)%name)//" (Gg)")
          write(*     ,"(a14,a5,3x,30(a12,:))")"EMTAB CC Land ","    ",EMIS_FILE(:)
          write(IO_LOG,"(a14,a5,3x,30(a12,:))")"EMTAB CC Land ","    ",EMIS_FILE(:)
          sumEU(:) = 0.0
          fmt="(a5,i4,1x,a9,3x,30(f12.2,:))"
          do ic = 1, NLAND
             do is = 1,NSECTORS
                sumemis(ic,:) = sumemis(ic,:) + sumemis_sec(ic,is,:)
             end do
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
          if(SecEmisTotalsWanted)then
             write(*,*)
             write(*,*)"Totals per sectors (only when non zero)"
             write(*,"(a19,2x,30(a12,:))")"  CC Land    Sector",EMIS_FILE(:)
             write(IO_LOG,*)
             write(IO_LOG,*)"Totals per sectors"
             write(IO_LOG,"(a19,2x,30(a12,:))")"  CC Land    Sector",EMIS_FILE(:)
             fmt="(i4,1x,a9,i4,3x,30(f12.2,:))"
             do ic = 1, NLAND
                do is = 1, NSECTORS
                   ccsum = sum( sumemis_sec(ic,is,:) )
                   icc=Country(ic)%icode
                   if ( ccsum > 0.0 )then
                      write(*, fmt) icc, Country(ic)%code, is, sumemis_sec(ic,is,:)
                      write(IO_LOG, fmt) icc, Country(ic)%code, is, sumemis_sec(ic,is,:)
                   end if
                end do
             end do
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
       write(*,*) dtxt//"CONV:me, nydays, gridwidth = ",me,nydays,GRIDWIDTH_M
       write(*,*) dtxt//"No. days in Emissions: ", nydays
       write(*,*) dtxt//"tonne_to_kgm2s in Emissions: ", tonne_to_kgm2s
       write(*,*) dtxt//"Emissions sums:"
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

    deallocate(sumemis_sec, sumemis_local)

    !additional new format initializations that must be placed after splits are defined
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
                !lonlat format must be applied after the 2D fields are read
             else
                !see if the species belongs to any of the splitted species
                iqrc = 0
                do iem = 1,NEMIS_FILE
                   do f = 1,emis_nsplit(iem)
                      iqrc = iqrc + 1
                      itot = iqrc2itot(iqrc)
                      if(trim(species(itot)%name)==trim(Emis_source(n)%species))then
                         if(DEBUG%EMISSIONS .and. MasterProc) then
                            write(*,'(a,i6,4i4,a,f12.3)')&
                              trim(Emis_source(n)%species)//' included in '//trim(EMIS_FILE(iem)),&
                              n, isec, iland, Emis_source(n)%country_ix, iem, &
                              trim( species(itot)%name ),  e_fact(isec,iland,iem)
                         end if
                         Emis_source(n)%factor = Emis_source(n)%factor * e_fact(isec,iland,iem)
                          go to 888
                      endif
                   enddo
                enddo ! iem
888             continue
             endif
          endif

       endif


    enddo ! n = 1, NEmis_sources


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
subroutine EmisSet(indate)   !  emission re-set every hour
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
  integer :: i, j, ij, is, k, n, f    ! coordinates, loop variables
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
  logical :: debug_tfac, dbgPoll! dbgPoll set for cc=HU, poll=NH3 now
  integer :: tfac_idx, emish_idx, split_idx, isec_idx

  ! If timezone=-100, calculate time based on longitude rather than timezone
  integer :: hour_iland,nstart
  integer :: i_Emis_4D, iem_treated
  character(len=125) ::varname
  TYPE(timestamp)   :: ts1,ts2
  integer :: iwday,id
  real :: daynorm, roadfac
  character(len=*), parameter:: dtxt='EmisSet:'
  character(len=TXTLEN_FILE) :: filename
  real, allocatable :: Emisfac2D(:,:,:)

  ! Initialize
  ehlpcom0 = GRAV* 0.001*AVOG!0.001 = kg_to_g / m3_to_cm3

  ! Scaling for totemadd:
  dtgrid = 3600 * GRIDWIDTH_M * GRIDWIDTH_M !NB: once every hour, not dt_advec

  ! The emis array only needs to be updated every full hour. The
  ! time-factor calculation needs to know if a local-time causes a shift
  ! from day to night.  In addition, we reset an overall day's time-factors
  ! at midnight every day.

  hourchange = (indate%hour/=oldhour).or.(indate%day/=oldday)
!BUG  hour_iland = -9999
  if(hourchange) then
    oldhour = indate%hour
    if(indate%day/=oldday)then
      !==========================
      call NewDayFactors(indate)
      if(USES%DEGREEDAY_FACTORS) call DegreeDayFactors(daynumber) ! => fac_emm, fac_edd
      !==========================
      wday=day_of_week(indate%year,indate%month,indate%day)
      if(wday==0)wday=7 ! Sunday -> 7
      oldday = indate%day
    end if

    if(Found_Emis_4D>0)then
      if(.not.allocated(Emis_4D))allocate(Emis_4D(LIMAX,LJMAX,KMAX_MID,N_Emis_4D))
      Emis_4D = 0.0 !default value, must be set to zero when no values are found
      NTime_Read=-1 !read all times
      call ReadTimeCDF(emis_inputlist(Found_Emis_4D)%Name,TimesInDays,NTime_Read)
      call CheckStop(NTime_Read>size(TimesInDays), dtxt//"Emissions_mod: increase size of TimesInDays ")
      !if(MasterProc)write(*,*)('found date ',i,TimesInDays(i),i=1,NTime_Read)
      !write(*,*)'compare  ',ts1,ts2
      ts1=make_timestamp(indate)
      do i=1,NTime_Read
        call nctime2date(ts2,TimesInDays(i))
        if(nint(tdif_secs(ts1,ts2))==0)then
          if(MasterProc)write(*,*)dtxt//'Emis_4D: found matching date ',i,TimesInDays(i)
          nstart=i
          exit
        end if
      end do
      if(i>NTime_Read )then
        if(MasterProc)then
          write(*,*)dtxt//'Emis_4D: WARNING DID NOT FIND ANY MATCHING DATE '
          write(*,*)dtxt//'Emis_4D: first date found ',TimesInDays(1)
          write(*,*)dtxt//'Emis_4D: last date found ',TimesInDays(NTime_Read)
          write(*,*)dtxt//'Emis_4D: difference to last date ',tdif_secs(ts1,ts2)/3600,' hours'
        end if
      else
        do i_Emis_4D=1,N_Emis_4D
          if(emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)=='NOTSET')exit
          varname=emis_inputlist(Found_Emis_4D)%pollName(i_Emis_4D)
          !if(MasterProc)write(*,*)'Fetching ',trim(varname)
          call GetCDF_modelgrid(varname,emis_inputlist(Found_Emis_4D)%Name,&
            Emis_4D(1,1,1,i_Emis_4D),1,kmax_mid,nstart,1,reverse_k=.true.)
        end do
      end if !i
    end if !Found_Emis_4D>0
  end if ! hourchange

  if(DEBUG%EMISTIMEFACS.and.MasterProc) &
    write(*,"(a,i4,2f8.3, L2)") dtxt//" traffic? 24x7", N_TFAC, &
      fac_ehh24x7(1,TFAC_IDX_TRAF,1,4,1),fac_ehh24x7(1,TFAC_IDX_TRAF,13,4,1),&
      hourchange
  !..........................................

  if(hourchange) then
    totemadd(:)  = 0.
    gridrcemis(:,:,:,:) = 0.0
    SecEmisOut(:,:,:,:) = 0.0
    if(USES%LocalFractions) emis_lf_cntry(:,:,:,:,:) = 0.0
    if(USES%ROADDUST) gridrcroadd0(:,:,:) = 0.0
    !..........................................

    if (any(hourly_emisfac%file /= 'NOTSET')) then
       allocate(Emisfac2D(LIMAX,LJMAX,NEMIS_File))
       Emisfac2D(1,1,:) = -1.0 ! to indicate that it is not set
       do i = 1, size(hourly_emisfac)
          iem = find_index(hourly_emisfac(i)%poll,EMIS_FILE(:))
          call CheckStop(iem<1 .and. hourly_emisfac(i)%file /= 'NOTSET', &
               trim(hourly_emisfac(i)%poll)//' is not a valid emitted pollutant name')
          if (iem<1) cycle
          n =  current_date%hour + 1 ! record
          filename = date2string(hourly_emisfac(i)%file,current_date,mode='YMDH')
          filename = key2str(filename,'DataDir',DataDir)
          filename = key2str(filename,'DataDir',EmisDir)
          call ReadField_CDF(trim(filename), hourly_emisfac(i)%cdfname, Emisfac2D(1,1,iem), n,&
               interpol='conservative',needed=.true.,debug_flag=.false.,UnDef=1.0)
       end do
    end if

    ! Process each grid:
    if(DEBUG%EMISTIMEFACS.and.debug_proc)write(*,*)'CAMEOiccloop0',NSECTORS, maxval(nlandcode) ! ZERO??
    do j = 1,ljmax
      do i = 1,limax
        ncc = nlandcode(i,j)            ! No. of countries in grid
        debug_tfac=(DEBUG%EMISTIMEFACS.and.debug_proc &
                     .and.i==DEBUG_li.and.j==DEBUG_lj)
        ! find the approximate local time:
        lon = modulo(360+nint(glon(i,j)),360)
        if(lon>180.0)lon=lon-360.0

        !*************************************************
        ! First loop over sector emissions
        !*************************************************
        tmpemis(:)=0.
        if(debug_tfac)write(*,'(a,4i5,2f8.3)')'CAMEOiccloop', ncc,NSECTORS, i_fdom(i),j_fdom(j), glon(i,j), glat(i,j) !CAMEO ncc=0
        do icc = 1, ncc
          !          iland = landcode(i,j,icc)     ! 1=Albania, etc.
          iland=find_index(landcode(i,j,icc),Country(:)%icode) !array index

          call make_iland_for_time(debug_tfac, indate, i, j, iland, wday,&
                 iland_timefac,hour_iland,wday_loc,iland_timefac_hour)

          !if ( hour_iland<1) print *, dtxt//'TZWRONGA', iland, hour_iland
          !if ( hour_iland<1) print *, dtxt//'TZWRONGZ', iland, hour_iland

          !  As each emission sector has a different diurnal profile
          !  and possibly speciation, we loop over each sector, adding
          !  the found emission rates to gridrcemis as we go.
          !  ==================================================
          do isec = 1, NSECTORS       ! Loop over all defined sectors
            ! NB: "isec" is just an indice and not necessarily meaningful
            ! (i.e isec=1 does not always mean "industry" or so, but just the
            ! first sector defined in SECTORS)
            ! Calculate emission rates from secemis, time-factors, emis heights
            ! and if appropriate any speciation fraction (NEMIS_FRAC)
            do iem = 1, NEMIS_FILE
               if (SECTORS(isec)%species /='ALL' .and. SECTORS(isec)%species /= EMIS_FILE(iem)) cycle
               tfac_idx = SECTORS(isec)%timefac
               emish_idx = SECTORS(isec)%height
               split_idx = SECTORS(isec)%split

               if(timeFacs%Day_of_Year ) then
                  tfac = fac_dayofyear(isec, iland_timefac, iem, daynumber)& ! NB: use isec, not tfac_idx
                       * fac_ehh24x7(iem,tfac_idx,hour_iland,wday_loc,iland_timefac_hour)
               else
                  tfac = timefac(iland_timefac,tfac_idx,iem) &
                       * fac_ehh24x7(iem,tfac_idx,hour_iland,wday_loc,iland_timefac_hour)
               endif

               dbgPoll = (debug_tfac.and. EMIS_FILE(iem)=='nh3')
               dbgPoll = (debug_tfac.and. EMIS_FILE(iem)=='sox') ! CAMEO check
               dbgPoll = (debug_tfac.and. EMIS_FILE(iem)=='pm25') ! CAMEO check
               if (dbgPoll) then
                !if (isec==1) write(*,"(a,2i4,2f8.2,i6)")dtxt//"DAY TFAC loc:",&
                write(*,"(a,2i4,2f8.2,i6)")dtxt//"efacs DAY TFAC loc:",&
                     iland,isec, glon(i,j), glat(i,j), Country(iland)%timezone
                write(*,"(a,3i4,f8.3)")dtxt//"efacs DAY TFAC:",isec,tfac_idx,hour_iland,tfac
              end if

              !it is best to multiply only if USES%GRIDDED_EMIS_MONTHLY_FACTOR
              !in order not to access the array and waste cache if not necessary
              if(timeFacs%Monthly == 'GRIDDED') tfac = tfac * GridTfac(i,j,tfac_idx,iem)

              !Degree days - only SNAP-2
              if(USES%DEGREEDAY_FACTORS .and. IS_DOM(isec) .and. Gridded_SNAP2_Factors) then
                 oldtfac = tfac
                ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
                ! we make use of a baseload even for SNAP2
                tfac = ( fac_min(iland_timefac,tfac_idx,iem) & ! constant baseload
                     + ( 1.0-fac_min(iland_timefac,tfac_idx,iem) )* gridfac_HDD(i,j) ) &
                     * fac_ehh24x7(iem,tfac_idx,hour_iland,wday_loc,iland_timefac_hour)

                if(debug_tfac .and. indate%hour==12 .and. iem==1) &
                  write(*,"(a,3i3,2i4,7f8.3)") "SNAPHDD tfac ",  &
                    isec, tfac_idx,iland, daynumber, indate%hour, &
                    timefac(iland_timefac,tfac_idx,iem), t2_nwp(i,j,2)-273.15, &
                    fac_min(iland,tfac_idx,iem),  gridfac_HDD(i,j), tfac
             end if ! =============== HDD

              if (allocated(Emisfac2D)) then
                 if (Emisfac2D(1,1,iem)>-0.5) then
                    if(me==NPROC/2 .and. i==2 .and. j==2 .and. icc == 1 &
                        .and. isec==1) write(*,*)'Warning: reseting all timefactors, and using only gridded hourly factors for '//trim(EMIS_FILE(iem))
                    !the hourly factors are defined for this pollutant
                    tfac = Emisfac2D(i,j,iem)
                 end if
              end if

              s = tfac * secemis(isec,i,j,icc,iem)
              ! prelim emis sum kg/m2/s
              SecEmisOut(i,j,iem,0) = SecEmisOut(i,j,iem,0) + s !sum of all sectors
              if(SecEmisOutWanted(isec))SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) = &
                   SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) + s

              if (debug_tfac ) write(*,*)'CAMEONOTHERE START', iem, split_idx, iland

              do f = 1,emis_nsplit(iem)
                itot = iemsplit2itot(f, iem)
                iqrc = itot2iqrc(itot)
                tmpemis(iqrc) = s * emisfrac(iqrc,split_idx,iland)
                ! Add up emissions in kg
                totemadd(itot) = totemadd(itot) &
                     + tmpemis(iqrc) * dtgrid * xmd(i,j)
                if (debug_tfac ) then
                  write(*,*)'CAMEOTMP', iem, f, itot, iqrc, emisfrac(iqrc,split_idx,iland), tmpemis(iqrc)
                endif
              end do ! f
              if (debug_tfac ) write(*,*)'CAMEOTMP END xxxxxxxxxxxxxxxxx'
              if(USES%LocalFractions) call save_lf_emis(s,i,j,iem,isec,iland)

            end do ! iem

            !  Assign to height levels 1-KEMISTOP
            do k=KEMISTOP,KMAX_MID
              do iqrc =1, nrcemis
                gridrcemis(iqrc,k,i,j) = gridrcemis(iqrc,k,i,j)   &
                     + tmpemis(iqrc)*ehlpcom0    &
                     *emis_kprofile(KMAX_BND-k,emish_idx) &
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

            ncc = road_nlandcode(i,j) ! number of countries in grid point
            do icc = 1, ncc
              !iland = road_landcode(i,j,icc)
              iland = find_index(road_landcode(i,j,icc),Country(:)%icode)
              iland_timefac_hour = find_index(Country(iland)%timefac_index_hourly,Country(:)%icode)
             !print *, 'ROAD ', icc, ncc, iland, hour_iland
             ! Needed here since can be called without hour change
              call make_iland_for_time(debug_tfac, indate, i, j, iland, wday,&
                 iland_timefac,hour_iland,wday_loc,iland_timefac_hour)

              roadfac = 1.0
              if(ANY(iland==(/IC_FI,IC_NO,IC_SE/)).and. & ! Nordic countries
                 ANY(indate%month==(/3,4,5/)))then        ! spring road dust
                 roadfac=2.0
              end if

              if ( hour_iland<1) then
                  print *, dtxt//'HOURCHANGE', hourchange
                  print *, dtxt//'TZWRONGB', iland, hour_iland
                  print '(a,9i5)', 'ROAD TFAC :', iland, iem, Country(iland)%timezone, hour_iland
                  call StopAll(dtxt//'ROADDUST')
              end if

              do iem = 1, NROAD_FILES
                 tfac = fac_ehh24x7(iem, TFAC_IDX_TRAF,hour_iland,wday_loc,iland_timefac_hour)
                 s = tfac * roadfac * roaddust_emis_pot(i,j,icc,iem)
                 if(DEBUG%ROADDUST.and.debug_proc.and.i==DEBUG_li.and.j==DEBUG_lj)&
                  write(*,*)dtxt//"DEBUG ROADDUST! iem,tfac,icc,roaddust_emis_pot,s", &
                    iem,tfac,icc,roaddust_emis_pot(i,j,icc,iem),s

                gridrcroadd0(QROADDUST_FI,i,j)=gridrcroadd0(QROADDUST_FI,i,j) &
                     +ROADDUST_FINE_FRAC*s
                gridrcroadd0(QROADDUST_CO,i,j)=gridrcroadd0(QROADDUST_CO,i,j) &
                     +(1.-ROADDUST_FINE_FRAC)*s

                if(all([DEBUG%ROADDUST,debug_proc,i==debug_li,j==debug_lj]))then
                  write(*,*)dtxt//"gridrcroadfine"  ,gridrcroadd0(QROADDUST_FI,i,j)
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
    do ij = 1,limax*ljmax
       i=mod(ij-1,limax)+1
       j=(ij-1)/limax+1

       debug_tfac=(DEBUG%EMISTIMEFACS.and.debug_proc&
                    .and.i==DEBUG_li.and.j==DEBUG_lj)
       if (debug_tfac) write(*,'(a,i5,2f8.3,i4)') dtxt//'XCAMEO: hereNEW', &
         NEmis_source_ij(ij), glat(i,j), glon(i,j), size(SECTORS)

       do is = 1,NEmis_source_ij(ij)
          n = Emis_source_ij_ix(ij,is)

          call CheckStop(allocated(Emisfac2D), 'hourly gridded emisfactors not implemented yet for emis new formats')
          itot = Emis_source(n)%species_ix
          isec = Emis_source(n)%sector
          isec_idx = Emis_source(n)%sector_idx !index in SECTORS
          tfac_idx = SECTORS(isec_idx)%timefac
          emish_idx = SECTORS(isec_idx)%height
          split_idx = SECTORS(isec_idx)%split

          iland = Emis_source(n)%country_ix ! remember: iland isn't emep cc num!
          if(debug_tfac) write(*,'(2a,99i5)') dtxt//'NewFormStart'// &
            ':'//trim(Emis_source(n)%periodicity), &
               ' '//trim(SECTORS(isec_idx)%longname)//'=>(tfac):'// & 
                    trim(SECTORS(tfac_idx)%longname), &
              is, NEmis_source_ij(ij), iland, itot,isec,isec_idx,tfac_idx

          if(itot>0)then
             !the species is directly defined (no splits)
             iqrc = itot2iqrc(itot)
             if(isec>0)then
                call CheckStop(iqrc<=0,dtxt// &
                     "emitted sector species must be one of the splitted species")
                iem = iqrc2iem(iqrc)

                if(debug_tfac) write(*,'(a,9i5)') dtxt//'NewForm-itot'// &
                  trim(EMIS_FILE(iem))//':'//trim(Emis_source(n)%periodicity), &
                    is, iland, itot,isec,isec_idx

                if(Emis_source(n)%periodicity == 'monthly')then
                   !make normalization factor for daily fac
                   iland_timefac = find_index(Country(iland)%timefac_index,Country(:)%icode)
                   iwday = mod(wday-indate%day+35, 7) + 1 !note both indate%day and wday start at 1
                   daynorm = 0.0
                   do id = 1, nmdays(indate%month)
                      !TODO: CHECK  fac_edd for isec. Why not use timefac?
                      daynorm = daynorm + fac_edd(iland_timefac,iwday,tfac_idx,iem)
                      iwday = iwday + 1
                      if(iwday>7)iwday = 1
                   enddo
                   daynorm = nmdays(indate%month)/daynorm
                endif

                if(Emis_source(n)%periodicity == 'yearly' .or. Emis_source(n)%periodicity == 'monthly')then
                   !we need to apply hourly factors
                   dbgPoll = (debug_tfac.and. EMIS_FILE(iem)=='nh3')

                   call make_iland_for_time(debug_tfac, indate, i, j, iland,&
                           wday, iland_timefac,hour_iland,wday_loc,iland_timefac_hour)
                   tfac = fac_ehh24x7(iem,tfac_idx,hour_iland,wday_loc,iland_timefac_hour)
                   if(debug_tfac) write(*,'(a,4i5,es12.3)') dtxt//'NH3tfacLand', &
                           is, iland,  iland_timefac, tfac

                   if(timeFacs%Day_of_Year ) then
                      tfac = tfac * fac_dayofyear(isec_idx, iland_timefac, iem, daynumber)
                      if(dbgPoll) write(*,*) dtxt//'efacs HERE-DOY', tfac, iland_timefac, daynumber
                   else if(Emis_source(n)%periodicity == 'yearly')then
                      !apply monthly and daily factor on top of hourly factors
                      tfac = tfac * timefac(iland_timefac,tfac_idx,iem)
                      if(dbgPoll) write(*,*) dtxt//'efacs HERE-MD', tfac, iland_timefac, daynumber
                   else if(Emis_source(n)%periodicity == 'monthly')then
                      !apply daily factors, with renormalization to conserve monthly sums
                      tfac = tfac * fac_edd(iland_timefac,wday,tfac_idx,iem) * daynorm
                      if(dbgPoll) write(*,*) dtxt//'efacs HERE-MM', tfac, iland_timefac, daynumber
                   endif
                else
                   !not monthly or yearly emissions, timefactors must be included in emission values
                   tfac = 1.0
                   if(debug_tfac) write(*,'(a,9i5)') dtxt//'NewFormSet1.0',&
                           is, iland,itot,isec,isec_idx
                endif

                s = Emis_source_ij(ij,is) * tfac

                !CAMEOPT
                !debug_tfac = .false.
                !if ( iland==18 .and. SECTORS(isec_idx)%longname=='GNFR_Bb' .and. s>0.0) then
                !   write(*,'(a,3i4,es12.3)') "CAMEOPT", isec_idx, i_fdom(i), j_fdom(j), s
                !end if

                SecEmisOut(i,j,iem,0) = SecEmisOut(i,j,iem,0) + s !sum of all sectors
                if(SecEmisOutWanted(isec))SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) = &
                     SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) + s
                ! Add up emissions in kg
                totemadd(itot) = totemadd(itot) + s * dtgrid * xmd(i,j)

                if(USES%LocalFractions .and. me==0.and.i==1.and.j==1.and.n==1) write(*,*)dtxt//&
                     'WARNING: single emitted spec not implemented for Local Fractions yet'

                !  Assign to height levels 1-KEMISTOP
                do k=KEMISTOP,KMAX_MID
                   gridrcemis(iqrc,k,i,j) = gridrcemis(iqrc,k,i,j)   &
                        + s&
                        *ehlpcom0    &
                        *emis_kprofile(KMAX_BND-k,emish_idx) &
                        *emis_masscorr(iqrc) !NB: assumes mass defined as a split nox, pm25... !
                end do   ! k
                !   end do ! i
                !enddo
             else
                ! we do not include the emissions as a sector emission
                !directly included in setup_rcemis
             endif
          else
             !the species is defined as a sector emission, e.g. sox, nox
             iem=find_index(Emis_source(n)%species,EMIS_FILE(:))

             if(debug_tfac) write(*,'(a,9i5)') dtxt//'SecEmis'//trim(EMIS_FILE(iem))//':'//&
              trim(Emis_source(n)%species)//trim(Emis_source(n)%periodicity ), &
               is, iland, iem
             call CheckStop(iem<0, dtxt//"did not recognize species "//trim(Emis_source(n)%species))
             call CheckStop(Emis_source(n)%sector<=0,dtxt//" sector must be defined for "//trim(Emis_source(n)%varname))

             iland_timefac = find_index(Country(iland)%timefac_index,Country(:)%icode)
             if(Emis_source(n)%periodicity == 'monthly')then
                !make normalization factor for daily fac
                iwday = mod(wday-indate%day+35, 7) + 1
                daynorm = 0.0
                do id = 1, nmdays(indate%month)
                   daynorm = daynorm + fac_edd(iland_timefac,iwday,tfac_idx,iem)
                   iwday = iwday + 1
                   if(iwday>7)iwday = 1
                enddo
                daynorm = nmdays(indate%month)/daynorm
             endif

             do f = 1,emis_nsplit(iem)
                itot = iemsplit2itot(f,iem)
                call CheckStop(itot<0, "did not recognize split "//trim(Emis_source(n)%species))
                iqrc = itot2iqrc(itot)

                if(debug_tfac) write(*,'(a,99i5)') dtxt//'NewFormsplit'//&
                   trim(Emis_source(n)%species)//':'//trim(species(itot)%name), &
                    n, iem, iland, itot, iqrc, emis_nsplit(iem)
                   !trims(EMIS_FILE(iem)//':'//Emis_source(n)%species//':'//species(itot)%name), &

                dbgPoll = (debug_tfac.and. EMIS_FILE(iem)=='nh3' .and. iland_timefac==IC_HU )
                dbgPoll = (debug_tfac.and. EMIS_FILE(iem)=='pm25')   ! XCAMEO

                if(Emis_source(n)%periodicity == 'yearly' .or. Emis_source(n)%periodicity == 'monthly')then
                   !we need to apply hourly factors
                   call make_iland_for_time(debug_tfac, indate, i, j, iland, &
                           wday, iland_timefac,hour_iland,wday_loc,iland_timefac_hour)
                   tfac = fac_ehh24x7(iem,tfac_idx,hour_iland,wday_loc,iland_timefac_hour)

                   if(timeFacs%Day_of_Year ) then
                      tfac = tfac * fac_dayofyear(isec_idx, iland_timefac, iem, daynumber)
                      if(dbgPoll) write(*,'(a,3i4,f12.4)') &
                         dtxt//'efacs HERE-DOYhu', tfac_idx, isec_idx, &
                            daynumber, tfac
                   else if(Emis_source(n)%periodicity == 'yearly')then
                      !apply monthly and daily factor on top of hourly factors
                      tfac = tfac * timefac(iland_timefac,tfac_idx,iem)
                      !if(dbgPoll) write(*,'(a,3i4,f12.4)') dtxt//&
                      if(debug_tfac) write(*,'(a,4i4,f12.4)') dtxt//trim(EMIS_FILE(iem))//&
                         'efacs HERE-CLIMhu', tfac_idx, isec_idx, daynumber, iland_timefac, tfac
                   else if(Emis_source(n)%periodicity == 'monthly')then
                      !apply daily factors, with renormalization to conserve monthly sums
                      tfac = tfac * fac_edd(iland_timefac,wday,tfac_idx,iem) * daynorm
                   endif
                else
                   !not monthly or yearly emissions, timefactors must be included in emission values
                   tfac = 1.0
                   if ( dbgPoll ) write(*,'(a,4i4,f12.4)') dtxt//'efacsUnset HERE', IC_HU, &
                          tfac_idx, isec_idx, daynumber, tfac
                endif

                if(timeFacs%Monthly == 'GRIDDED') tfac=tfac* GridTfac(i,j,tfac_idx,iem)

                !Degree days - only SNAP-2
                if(USES%DEGREEDAY_FACTORS .and. &
                     IS_DOM(isec_idx) .and. Gridded_SNAP2_Factors .and. &
                     (Emis_source(n)%periodicity == 'yearly' .or. &
                      Emis_source(n)%periodicity == 'monthly')) then
                   if(dbgPoll) write(*,*) dtxt//'efacs HERE-PreHDD', tfac,& !CAMEO Gridded_SNAP2????
                           isec,daynumber
                   oldtfac = tfac
                   ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
                   ! we make use of a baseload even for SNAP2
                   tfac = ( fac_min(iland_timefac,tfac_idx,iem) & ! constant baseload
                        + ( 1.0-fac_min(iland_timefac,tfac_idx,iem) )* gridfac_HDD(i,j) ) &
                        * fac_ehh24x7(iem,tfac_idx,hour_iland,wday_loc,iland_timefac_hour)

                   if(debug_tfac .and. indate%hour==12 .and. iem==1) &
                        write(*,"(a,3i3,2i4,7f8.3)") dtxt//"SNAPHDD tfac ",  &
                        isec, tfac_idx,iland, daynumber, indate%hour, &
                        timefac(iland_timefac,tfac_idx,iem), t2_nwp(i,j,2)-273.15, &
                        fac_min(iland_timefac,tfac_idx,iem),  gridfac_HDD(i,j), tfac
                   if(dbgPoll) write(*,*) dtxt//'efacs HERE-HDD', tfac, isec,daynumber
                end if ! =============== HDD

                if(debug_tfac .and. emisfrac(iqrc,split_idx,iland) >0.0) &
                 write(*,'(a,6i5,f12.6)') dtxt//'NewFormSPLIt'//&
                   trims(EMIS_FILE(iem)//':'//species(itot)%name//':'//&
                         trim(SECTORS(isec_idx)%longname)), & 
                    n, iem, iland, itot, iqrc, split_idx, emisfrac(iqrc,split_idx,iland)

                s = Emis_source_ij(ij,is) * tfac
                
                !NB: LF must get unsplitted  emissions
                if(USES%LocalFractions .and. f==1) call save_lf_emis(s,i,j,iem,isec_idx,iland)

                s = s * emisfrac(iqrc,split_idx,iland)

                !if ( iland==2 .and. SECTORS(isec_idx)%longname=='GNFR_Cb' .and. s>0.0) then
                !   write(*,'(a,5i4,es12.3,1x,a)') dtxt//"CAMEOPTA", isec_idx,split_idx,itot, i_fdom(i), j_fdom(j), s &
                !   ,trim(species(itot)%name)
                !else if ( iland==2 .and. SECTORS(isec_idx)%longname=='GNFR_Af' .and. s>0.0) then
                !   write(*,'(a,4i4,es12.3,1x,a)') dtxt//"CAMEOPTB", isec_idx,itot, i_fdom(i), j_fdom(j), s &
                !   ,trim(species(itot)%name)
                !else if ( iland==2 .and. SECTORS(isec_idx)%longname=='GNFR_Af') then
                !   write(*,'(a,3i4,es12.3,a)') dtxt//"CAMEOPTC", isec_idx, i_fdom(i), j_fdom(j), s &
                !   ,trim(species(itot)%name)
                !end if

                SecEmisOut(i,j,iem,0) = SecEmisOut(i,j,iem,0) + s !sum of all sectors
                if(SecEmisOutWanted(isec))SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) = &
                     SecEmisOut(i,j,iem,isec2SecOutWanted(isec)) + s
                ! Add up emissions in kg
                totemadd(itot) = totemadd(itot) + s * dtgrid * xmd(i,j)


                !  Assign to height levels 1-KEMISTOP
                do k=KEMISTOP,KMAX_MID
                   gridrcemis(iqrc,k,i,j) = gridrcemis(iqrc,k,i,j)   &
                        + s&
                        *ehlpcom0    &
                        *emis_kprofile(KMAX_BND-k,emish_idx) &
                        *emis_masscorr(iqrc)
                end do   ! k

             enddo
          endif
       enddo
    enddo
    if(allocated(Emisfac2D))deallocate(Emisfac2D)

 end if ! hourchange

  if(USES%ROADDUST)THEN
    if(DEBUG%ROADDUST.and.debug_proc) &
      write(*,*)dtxt//"Before the unit scaling",gridrcroadd(1:2,DEBUG_li,DEBUG_lj)
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
      write(*,*)dtxt//"After the unit scaling",gridrcroadd(1:2,DEBUG_li,DEBUG_lj)
  end if
end subroutine EmisSet
!***********************************************************************
subroutine newmonth
!----------------------------------------------------------------------!
! DESCRIPTION:
!   Reads in natural DMS emissions at start of each month. Update
!   landcode and nlandcode arrays as needed.
!
!   Reads in snow cover and global soil NOx  at start of each month.
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
  integer :: kstart,kend,nstart,Nyears,mm
  real :: buffer3D(8,LIMAX,LJMAX),buffer(LIMAX,LJMAX),SumSoilNOx,ccsum
  real :: fractions(LIMAX,LJMAX,NCMAX),Reduc(NLAND)
  real, dimension(NEMIS_FILE)       :: emsum ! Sum emis over all countries
  real, dimension(NLAND,NEMIS_FILE) :: sumemis, sumemis_local ! Sum of emissions per country
  real, dimension(NEMIS_FILE) :: sumEU ! Sum of emissions in EU
  character(len=40) :: varname , fmt
  character(len=125) ::fileName
  real :: Mask_ReducFactor
  integer :: NMask_Code,Mask_Code(NLAND), i_femis_lonlat,icc
  real :: lonlat_fac, mw, dbgVal
  logical :: use_lonlat_femis, monthlysectoremisreset, Cexist
  logical :: fractionformat
  integer :: emis_inputlist_NEMIS_FILE!number of files for each emis_inputlist(i)
 ! For AIRCRAFT femis work
  integer, save  :: iFemis, iemNOx
  character(len=*), parameter :: dtxt='Emis:newmonth:'

  monthlysectoremisreset =.false.

  if(.not.allocated(airn).and.(USES%LIGHTNING_EMIS.or.USES%AIRCRAFT_EMIS))&
    allocate(airn(KCHEMTOP:KMAX_MID,LIMAX,LJMAX))


  !F24 if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)&
  if(timeFacs%Monthly == 'GRIDDED') &
    call Read_monthly_emis_grid_fac(current_date%month)

  !Sep2023 always:  if(USES%TIMEZONEMAP)then
  call Read_monthly_timezones(current_date%month)
  !end if

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
             e_fact(TFAC_IDX_POW,iFemis,iemNOx)
    end if

    TotAircraftEmis = 0.0
    if (index(AircraftEmis_FLFile,'CAMS-GLOB-AIR')>0) then
       !assumes the CAMS aircraft emissions (NO in kg/m2/s)
       if(startdate(1)<2000)then
          AircraftEmis_FLFile = key2str(AircraftEmis_FLFile,'YYYY',2000)
       else if(startdate(1)>2020)then
          AircraftEmis_FLFile = key2str(AircraftEmis_FLFile,'YYYY',2020)
       else
          AircraftEmis_FLFile = key2str(AircraftEmis_FLFile,'YYYY',startdate(1))
       end if
       call ReadField_CDF_FL(AircraftEmis_FLFile,'avi',airn,&
            current_date%month,kstart,kend,USES%zero_below3000ft,&
            interpol='conservative', needed=.true.,debug_flag=.false.)
       ! convert from kg(NO)/m2/s into molecules/cm3/s. mw(NO)=30.0
       ! from kg to molecules: 1000*AVOG/30,
       ! 1e-6 for m-3->cm-3
       ! use roa to find dz for consistency with other emissions
       ! (otherwise could have used z_bnd directly)
       ! dz=dP/(roa*GRAV)  dP=dA(k) + dB(k)*ps(i,j,1)

        conv=1000*AVOG/30.0*GRAV*1.0e-6&
             * e_fact(TFAC_IDX_POW,iFemis,iemNOx)  ! sector is assumed to be the same as for "POWER"
        do k=KCHEMTOP,KMAX_MID
          do j=1,ljmax
             do i=1,limax
                TotAircraftEmis = TotAircraftEmis + airn(k,i,j)*xmd(i,j)
                airn(k,i,j)=airn(k,i,j)*conv*(roa(i,j,k,1))&
               /(dA(k)+dB(k)*ps(i,j,1))
             end do
          end do
       end do
       TotAircraftEmis = TotAircraftEmis /30.0*46.0*gridwidth_m**2*(nmdays(current_date%month)*24*3600)
    else
       !assumes NO2 in kg/month
       call ReadField_CDF_FL(AircraftEmis_FLFile,'NOx',airn,&
            current_date%month,kstart,kend,USES%zero_below3000ft,&
            interpol='mass_conservative', needed=.true.,debug_flag=.false.)

       ! convert from kg(NO2)/month into molecules/cm3/s
       ! from kg to molecules: 1000*AVOG/species(NO2)%molwt
       ! use roa to find dz for consistency with other emissions
       ! (otherwise could have used z_bnd directly)
       ! dz=dP/(roa*GRAV)  dP=dA(k) + dB(k)*ps(i,j,1)
       ! dV=dz*dx*dy=dz*gridwidth**2/xm**2 *1e6 (1e6 for m3->cm3)
       ! from month to seconds: ndaysmonth*24*3600

       conv=1000*AVOG/species(NO2)%molwt*GRAV/gridwidth_m**2*1.0e-6&
            /(nmdays(current_date%month)*24*3600)

       do k=KCHEMTOP,KMAX_MID
          do j=1,ljmax
             do i=1,limax
                TotAircraftEmis = TotAircraftEmis + airn(k,i,j)
                airn(k,i,j)=airn(k,i,j)*conv*(roa(i,j,k,1))&
                    * e_fact(TFAC_IDX_POW,iFemis,iemNOx) &  ! sec1, emis1 AIRCRAFT
                     /(dA(k)+dB(k)*ps(i,j,1))*xm2(i,j)
             end do
          end do
       end do
    end if

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,TotAircraftEmis,1,MPI_DOUBLE_PRECISION, &
       MPI_SUM,MPI_COMM_CALC,IERROR)
    if (MasterProc) write(*,*) "Total NOx emissions from aircraft for coming month:"
    if (MasterProc) write(*,*) TotAircraftEmis," kg(NO2)/month"
  end if ! USES%AIRCRAFT_EMIS

  if(DEBUG%SOILNOX.and.debug_proc) write(*,*)"Emissions DEBUG_SOILNOX ????", me
  if(USES%SOILNOx .and. USES%SOILNOX_METHOD=='ACP2012EURO')then  ! European Soil NOx emissions

      ! read in map of annual N-deposition produced from pre-runs of EMEP model
      ! with script mkcdo.annualNdep
      call ReadField_CDF(NdepFile,'Ndep_m2',AnnualNdep,1,&
          interpol='zero_order',needed=.true.,debug_flag=.false.,UnDef=0.0)

      if(DEBUG%SOILNOX.and.debug_proc)&
        write(*,"(a,4es12.3)") dtxt//" SOILNOX AnnualDEBUG ", &
        AnnualNdep(debug_li, debug_lj), maxval(AnnualNdep), minval(AnnualNdep)

  elseif(USES%SOILNOX) then ! Global soil NOx, default from 2021

     !cf MonthlyDiurnalEmisFactor(months, tsteps, lat, lon)
    SoilNOx(:,:)=0.0
    SoilNOx3D(:,:,:)=0.0  ! careful: LIMAX,LJMAX,8
    buffer3D(:,:,:)=0.0   !          8,LIMAX,LJMAX
    buffer(:,:)=0.0
    if(debug_proc) write(*,*)dtxt//' SOILNOx start', me,&
            ' '//trim( USES%SOILNOX_METHOD )//trim(soilnox_emission_File)

    if ( USES%SOILNOX_METHOD == 'Total' .or. USES%SOILNOX_METHOD =='NoFert' ) then


   if (  index(soilnox_emission_File,'v2.4a_GLOBAL05_Clim2000') > 0 .or. & ! EMEP-style
         index(soilnox_emission_File,'v2.4clim') > 0) then                 ! ECCAD-style
     ! New format: tsteps in 4D array. No year info, just month.

      nstart=current_date%month
      i=debug_li
      j=debug_lj

      call ReadField_CDF(soilnox_emission_File,&
           'TotalSoilEmis',buffer3D,&
           nstart=current_date%month, kstart=1, kend=8,& !same variation every year
           interpol='conservative',known_projection="lon lat",&
           needed=.true.,debug_flag=DEBUG%SOILNOX,UnDef=0.0)
      do n=1,8 ! hour= 1.5, 4.5, 7.5 ... 22.5
          SoilNOx3D(:,:,n)=buffer3D(n,:,:)
      end do
      if ( debug_proc) write(*,"(a40,i3,2es12.3)") dtxt//'CLIMSOIL monthTot:',&
         current_date%month, maxval(SoilNOx3D), SoilNOx3D(i,j,1)

      if (USES%SOILNOX_METHOD == 'NoFert') then
        !we must substract Fertilizer emissions
        buffer3D(:,:,:)=0.0   !          8,LIMAX,LJMAX
        call ReadField_CDF(soilnox_emission_File,&
           'FertEmis',buffer3D,&
           nstart=current_date%month, kstart=1, kend=8,& !same variation every year
           interpol='conservative',known_projection="lon lat",&
           needed=.true.,debug_flag=.false.,UnDef=0.0)
        if ( debug_proc) write(*,"(a40,i3,2es12.3)") dtxt//'CLIMSOIL monthBuf:',&
            current_date%month, maxval(buffer3D), buffer3D(1,i,j)
        do n=1,8 ! hour= 1.5, 4.5, 7.5 ... 22.5
          SoilNOx3D(:,:,n)= SoilNOx3D(:,:,n) - buffer3D(n,:,:)
        end do
        if ( debug_proc) write(*,"(a40,i3,2es12.3)") dtxt//'CLIMSOIL monthFer:',&
           current_date%month, maxval(SoilNOx3D), SoilNOx3D(i,j,1)
      end if

      SoilNOx3D = max(0.0, SoilNOx3D)

      !just to get nice BioNat output in .nc files:
      do n=1,8
         SoilNOx(:,:) = SoilNOx(:,:) + SoilNOx3D(:,:,n)/8.0
      end do
      if ( debug_proc) write(*,"(a40,i3,2es12.3)") dtxt//'CLIMSOIL monthOut:',&
         current_date%month, maxval(SoilNOx3D), SoilNOx3D(i,j,1)

   else

      ! safety for < 2000 set here. (more likely to cause errors with
      ! trends)
      !call CheckStop(current_date%year>2018,dtxt//'CAMS81 soil>2018')
      call CheckStop(current_date%year<2000,dtxt//'CAMS81 soil<2000')

      ! CRUDE - will add variable later.
      ! 19 or 20 years defined in files:
      ! 12*19=228 months defined for < v2.4
      ! 12*20=240 months defined for   v2.4
      if ( index(soilnox_emission_File,'v2.4') > 0) then
        nstart=(min(2019,max(2000,current_date%year))-2000)*12 + current_date%month
      else
        nstart=(min(2018,max(2000,current_date%year))-2000)*12 + current_date%month
      end if
      if ( debug_proc) write(*,*) 'YYSOIL file', trim(soilnox_emission_File)// 'XXX', nstart

      call ReadField_CDF(soilnox_emission_File,&
           'TotalSoilEmis',SoilNOx,&
           nstart=nstart,interpol='conservative',known_projection="lon lat",&
           needed=.false.,debug_flag=.false.,UnDef=0.0)

      if (USES%SOILNOX_METHOD == 'NoFert') then
         !we must substract Fertilizer emissions
         call ReadField_CDF(soilnox_emission_File,&
           'FertEmis',buffer,&
           nstart=nstart,interpol='conservative',known_projection="lon lat",&
           needed=.true.,debug_flag=.false.,UnDef=0.0)
         SoilNOx = max(0.0, SoilNOx - buffer)
      end if

      call ReadField_CDF(soilnox_emission_File,&
           'MonthlyDiurnalEmisFactor',buffer3D,&
           nstart=current_date%month, kstart=1, kend=8,& !same variation every year
           interpol='conservative',known_projection="lon lat",&
           needed=.true.,debug_flag=.false.,UnDef=0.0)

      do n=1,8 ! hour= 1.5, 4.5, 7.5 ... 22.5
        do i=1,limax
          do j=1,ljmax
            SoilNOx3D(i,j,n)=buffer3D(n,i,j) * SoilNOx(i,j)
          end do
        end do
        if(debug_proc) write(*,*) dtxt//'3DSOIL ', current_date%month, i, &
            SoilNOx3D(debug_li,debug_lj,n)
      end do
    end if

      if(DEBUG%SOILNOX.and.debug_proc) then
        write(*,*) dtxt//"CAMS81 SOILNO ", current_date%year, current_date%month, nstart, maxval(SoilNOx)
        do n=1,8 ! hour= 1.5, 4.5, 7.5 ... 22.5
          write(*,*) dtxt//'3DSOIL ', current_date%month, n, SoilNOx3D(debug_li,debug_lj,n)
        end do
      end if

    else if (USES%SOILNOX_METHOD == 'Zaehle2011') then
      nstart=(current_date%year-1996)*12 + current_date%month
      if(nstart>0.and.nstart<=120)then
        !the month is defined
        call ReadField_CDF(soilnox_emission_File,'NOX_EMISSION',SoilNOx,&
             nstart=nstart,interpol='conservative',known_projection="lon lat",&
             needed=.true.,debug_flag=.false.,UnDef=0.0)
        if(DEBUG%SOILNOX.and.debug_proc) write(*,*) dtxt// &
          "PROPER YEAR of SOILNO ", current_date%year, nstart
      else
        !the year is not defined; average over all years
        Nyears=10 !10 years defined
        do iyr=1,Nyears
          nstart=12*(iyr-1) + current_date%month
          call ReadField_CDF(soilnox_emission_File,'NOX_EMISSION',buffer,&
              nstart=nstart,interpol='conservative',known_projection="lon lat",&
              needed=.true.,debug_flag=.false.,UnDef=0.0)
          do j=1,ljmax
            do i=1,limax
              SoilNOx(i,j)=SoilNOx(i,j)+buffer(i,j)
            end do
          end do
          if(DEBUG%SOILNOX.and.debug_proc) &
            write(*,"(a,2i6,es10.3,a,2es10.3)") dtxt//&
             "Averaging SOILNO  inputs", 1995+(iyr-1), nstart,&
              SoilNOx(debug_li,debug_lj),"max: ",maxval(buffer),maxval(SoilNOx)
        end do !iyr
        SoilNOx=SoilNOx/Nyears
      end if ! nstart test

    else !  SOILNOX_METHOD Needs to be Fert, NoFert, Zaehle2011, or ACP2012EURO
      call StopAll(dtxt//'WRONG SNOX METHOD:'//USES%SOILNOX_METHOD )
    end if ! CAMS81

     if(DEBUG%SOILNOX.and.debug_proc) then
       write(*,"(a,i3,4es10.3)") dtxt//"After Global SOILNO ",&
             me,maxval(SoilNOx),SoilNOx(debug_li,debug_lj), sum(SoilNOx3D(debug_li,debug_lj,:))/8.0
     end if
  else ! no soil NO
    if(DEBUG%SOILNOX.and.debug_proc) &
      write(*,*) dtxt//"Emissions DEBUG_SOILNOX - none"
  end if !  SOIL NO

  !for testing, compute total soil NOx emissions within domain
  if(USES%SOILNOX .and. USES%SOILNOX_METHOD /= 'ACP2012EURO') then
    SumSoilNOx=0.0
    SoilNOx = max(0.0, SoilNOx)  ! Stops the NEGs!
    ! CAMS uses kg(NO)/m2/s, so we convert to g(N)/m2/day to match earlier Zaehle
    if ( USES%SOILNOX_METHOD /= 'Zaehle2011' ) then ! have kg(NO)/m2/s
      do j=1,ljmax
        do i=1,limax
          SoilNOx(i,j) = SoilNOx(i,j) * 14./30 * 1000.0 * 24*3600 ! from kg(NO) to g(N)
          SoilNOx3D(i,j,:) = SoilNOx3D(i,j,:) * 14./30 * 1000.0 * 24*3600 ! from kg(NO) to g(N)
        end do
      end do
    end if ! CAMS81 test
   !convert from g(N)/m2/day into kg/day
    do j=1,ljmax
      do i=1,limax
        SumSoilNOx=SumSoilNOx+0.001*SoilNOx(i,j)*gridwidth_m**2*xmd(i,j)
      end do
    end do
    CALL MPI_ALLREDUCE(SumSoilNOx,mpi_out,1,MPI_DOUBLE_PRECISION, &
         MPI_SUM,MPI_COMM_CALC,IERROR)

    SumSoilNOx = mpi_out

    mm= current_date%month
    if(MasterProc)&
     write(*,'(a,i3,2(g12.3,a))') dtxt//&
       'GLOBAL SOILNOX emissions within domain, month:', mm, &
        SumSoilNOx, ' kg/day', SumSoilNOx*1.0e-9*nmdays(mm), ' Tg per month'

    ! convert from g(N)/m2/day into molecules/m2/s
    !  from g to molecules: AVOG/14  14=molweight N,
    !  from /day to /seconds : 24*3600
    conv=AVOG/14.0/(24*3600)
    if ( debug_proc) dbgVal = SoilNOx(debug_li,debug_lj)
    k=KMAX_MID!surface
    do j=1,ljmax
      do i=1,limax
        SoilNOx(i,j)=SoilNOx(i,j)*conv
        SoilNOx3D(i,j,:) = SoilNOx3D(i,j,:) * conv
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

  if (USES%OCEAN_DMS) then
     if(MasterProc)write(*,*)'reading OceanDMS'
     call ReadField_CDF(trim(DMSFile),'DMS_sea',O_DMS%emis,&
            nstart=current_date%month,interpol='conservative',known_projection="lon lat",&
            needed=.true.,debug_flag=.false.,UnDef=0.0)
     DMS%FileFound =.true.
     !from nanomol/l -> mol/cm3
     O_DMS%emis=O_DMS%emis*1.0e-12*DMS_S_FAC
  end if

  if (USES%OCEAN_NH3) then
      if(MasterProc)write(*,*)'reading OceanNH3'
       call ReadField_CDF(trim(OceanNH3File),'emiss_ocn',O_NH3%emis,&
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

       if(MasterProc)write(*,*)'Total monthly NH3 from Oceans (in Gg) ',O_NH3%sum_month
       O_NH3%sum_year=O_NH3%sum_year+O_NH3%sum_month!total for all month
  end if

  do iemislist = 1, size( emis_inputlist(:)%name )

    if(emis_inputlist(iemislist)%name == "NOTSET")cycle

    call StopAll("emis_inputlist style input no more in use. Use Emis_sourceFiles instead ")

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
          if ( debug_proc) write(*,*) dtxt//'XEmisB:'//basename(fname),maxval(secemis)!//':'&

       end do! iem = 1,NEMIS_FILE
       sumemis=0.0
       CALL MPI_REDUCE(sumemis_local,sumemis,&
            NLAND*NEMIS_FILE,MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)
       if(MasterProc)then
          call PrintLog("Total emissions by countries for "//&
                  trim(emis_inputlist(iemislist)%name)//" (Gg)")
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
        call StopAll("OceanNH3 no more an emis_inputlist. Use USES%OCEAN_NH3=T instead ")
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

       USES%OCEAN_NH3=.true.
       if(MasterProc)write(*,*)'Total monthly NH3 from Oceans (in Gg) ',O_NH3%sum_month
       O_NH3%sum_year=O_NH3%sum_year+O_NH3%sum_month!total for all month

    else if(emis_inputlist(iemislist)%type == 'DMS')then
       call StopAll("DMS no more an emis_inputlist. Use USES%OCEAN_DMS=T instead ")

       if(MasterProc)write(*,*)'reading DMS'
       call ReadField_CDF(emis_inputlist(iemislist)%name,'DMS_sea',O_DMS%emis,&
            nstart=current_date%month,interpol='conservative',known_projection="lon lat",&
            needed=.true.,debug_flag=.false.,UnDef=0.0)

       USES%OCEAN_DMS=.true.
       DMS%FileFound =.true.

       !from nanomol/l -> mol/cm3
       O_DMS%emis=O_DMS%emis*1.0e-12*DMS_S_FAC

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
  logical :: dbg 
  txt = trim(label)//"."//trim(EMIS_FILE(iem))
  msg(:) = 0

  dbg = DEBUG%EMISSIONS
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
!  end associate ! dbg

end subroutine EmisWriteOut

end module Emissions_mod
