! <Emissions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2013 met.no
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

                    module Emissions_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Calls up emission read/set routines
!  This routine interfaces the stand-alone emission-file reading routines
!  with the 3D model.
!_____________________________________________________________________________

  use Biogenics_ml,     only: SoilNOx, AnnualNdep, Ndep_trends
  use CheckStop_ml,only : CheckStop,StopAll
  use ChemSpecs_shl_ml, only: NSPEC_SHL
  use ChemSpecs_tot_ml, only: NSPEC_TOT,NO2
  use ChemChemicals_ml, only: species
  use Country_ml,    only : NLAND,Country_Init,Country,IC_NAT,IC_FI,IC_NO,IC_SE
  use Country_ml,    only : EU27, EUMACC2 !CdfSnap
  use EmisDef_ml, only : NSECTORS & ! No. sectors
                     ,NEMIS_FILE & ! No. emission files
                     ,EMIS_FILE   & ! Names of species ("sox  ",...)
                     ,NCMAX       & ! Max. No. countries per grid
                     ,FNCMAX      & ! Max. No. countries (with flat emissions)
                                    ! per grid
                     ,ISNAP_DOM   & ! snap index for domestic/resid emis
                     ,ISNAP_TRAF  & ! snap index for road-traffic (SNAP7)
                     ,ISNAP_SHIP  & ! snap index for ship emissions
                     ,ISNAP_NAT   & ! snap index for nat. (dms) emissions
                     ,IQ_DMS      & ! code for DMS emissions
                     ,NROAD_FILES & ! No. road dust emis potential files
                     ,ROAD_FILE   & ! Names of road dust emission files
                     ,NROADDUST   & ! No. road dust components 
                     ,QROADDUST_FI & ! fine road dust emissions (PM2.5) 
                     ,QROADDUST_CO & ! coarse road dust emis
                     ,ROADDUST_FINE_FRAC & ! fine (PM2.5) fraction of road dust emis
                     , ROADDUST_CLIMATE_FILE ! TEMPORARY! file for road dust climate factors 
  use EmisGet_ml, only : EmisGet, EmisSplit &
        ,EmisGetCdf                  &  ! cdfemis
        ,GridEmis                    &  ! cdfemis
        ,nGridEmisCodes               &  ! cdfemis
        ,GridEmisCodes               &  ! cdfemis
        ,sumcdfemis                    &  ! cdfemis
        ,femis                       &  ! Gets scaling factors -> e_fact
        ,EmisHeights                 &  ! Generates vertical distrib
        ,nrcemis, nrcsplit, emisfrac &  ! speciation routines and array
        ,nemis_kprofile, emis_kprofile &! Vertical emissions profile
        ,iqrc2itot                   &  ! maps from split index to total index
        ,emis_masscorr               &  ! 1/molwt for most species
        ,emis_nsplit                 &  ! No. species per emis file
        ,RoadDustGet                 &  
        ,roaddust_masscorr              ! 1/200. Hard coded at the moment, needs proper setting in EmisGet_ml...   
  use GridValues_ml, only:  GRIDWIDTH_M    & ! size of grid (m)
                           ,xm2            & ! map factor squared
                           ,debug_proc,debug_li,debug_lj & 
                           ,sigma_bnd, xmd, glat, glon,dA,dB,i_fdom,j_fdom
  use Io_Nums_ml,       only : IO_LOG, IO_DMS, IO_EMIS, IO_TMP
  use Io_Progs_ml,      only : ios, open_file, datewrite
  use MetFields_ml,     only : roa, ps, z_bnd, surface_precip ! ps in Pa, roa in kg/m3
  use MetFields_ml,     only : t2_nwp
  use ModelConstants_ml,only : KMAX_MID, KMAX_BND, PT ,dt_advec, &
                              EMIS_SOURCE, & ! emislist, CdfFractions
                              EMIS_TEST,   & ! CdfSnap or none
                              EMIS_OUT,   & ! output emissions in ASCII or not
                              IS_GLOBAL, & 
                              MONTHLY_GRIDEMIS, &  !NML
                              NBVOC,     &      ! > 0 if forest voc wanted
                              INERIS_SNAP2 , &     ! INERIS/TFMM HDD20 method
                              DEBUG => DEBUG_EMISSIONS,  MasterProc, & 
                              DEBUG_SOILNOX , DEBUG_EMISTIMEFACS, & 
                              DEBUG_ROADDUST , &
                              DEBUG_I,DEBUG_J, &
                              USE_DEGREEDAY_FACTORS, & 
                              NPROC, IIFULLDOM,JJFULLDOM , & 
                              SEAFIX_GEA_NEEDED, & ! see below
                              USE_LIGHTNING_EMIS,USE_AIRCRAFT_EMIS,USE_ROADDUST, &
                              USE_EURO_SOILNOX, USE_GLOBAL_SOILNOX   ! one or the other
  use NetCDF_ml, only  : ReadField_CDF
  use Par_ml,     only : MAXLIMAX,MAXLJMAX,me,gi0,gi1,gj0,gj1, &
                             GIMAX, GJMAX, IRUNBEG, JRUNBEG,  &   
                             limax,ljmax, &
                             MSG_READ1,MSG_READ7
  use PhysicalConstants_ml,  only :  GRAV,  AVOG
  use Setup_1dfields_ml,    only : rcemis   ! ESX
  use SmallUtils_ml,        only : find_index

  use ReadField_ml, only : ReadField    ! Reads ascii fields
  use TimeDate_ml,  only : nydays, nmdays, date, current_date, &! No. days per 
                           daynumber,day_of_week                ! year, date-type, weekday 
  use Timefactors_ml, only :   &
               NewDayFactors   &         ! subroutines
              ,DegreeDayFactors      &   ! degree-days used for SNAP-2
              ,Gridded_SNAP2_Factors, gridfac_HDD & 
              ,fac_min,timefactors   &                 ! subroutine
              ,fac_ehh24x7 ,fac_emm, fac_edd, timefac   ! time-factors
  use Volcanos_ml


  implicit none
  private


 ! subroutines:

  public :: Emissions                ! Main emissions module 
  public :: newmonth
  public :: EmisSet                  ! Sets emission rates every hour/time-step
  public :: EmisOut           ! Outputs emissions in ascii


  ! The main code does not need to know about the following 
  private :: consistency_check       ! Safety-checks

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO

  
 ! land-code information in each grid square - needed to know which country
 ! is emitting.                        
 ! nlandcode = No. countries in grid square
 ! landcode  = Country codes for that grid square
   integer, private, save, allocatable, dimension(:,:)       :: nlandcode
   integer, private, save, allocatable, dimension(:,:,:) :: landcode
 ! for flat emissions, i.e. no vertical extent:
   integer, private, save, allocatable, dimension(:,:)       :: flat_nlandcode
   integer, private, save, allocatable, dimension(:,:,:):: flat_landcode
 ! for road dust emission potentials:
   integer, private, save, allocatable, dimension(:,:)       :: road_nlandcode
   integer, private, save, allocatable, dimension(:,:,:) :: road_landcode

 !
 ! The output emission matrix for the 11-SNAP data is snapemis:
 !

  real, private, allocatable, dimension(:,:,:,:,:) &
            , save ::  snapemis      ! main emission arrays, in kg/m2/s

  real, private, allocatable, dimension(:,:,:,:) &
            , save ::  snapemis_flat ! main emission arrays, in kg/m2/s  

  real, private, allocatable, dimension(:,:,:,:) & ! Not sure if it is really necessary to keep the country info; gives rather messy code but consistent with the rest at least (and can do the seasonal scaling for Nordic countries in the code instead of as preprocessing) 
            , save ::  roaddust_emis_pot ! main road dust emission potential arrays, in kg/m2/s (to be scaled!)

 ! We store the emissions for output to d_2d files and netcdf in kg/m2/s

  real, public, allocatable, dimension(:,:,:), save :: SumSnapEmis

  logical, save, private  :: first_dms_read

  ! Emissions for input to chemistry routines

  ! KEMISTOP added to avoid hard-coded KMAX_MID-3:

   !integer, public, parameter :: KEMISTOP = KMAX_MID - NEMISLAYERS + 1
   integer, public, save :: KEMISTOP ! not defined yet= KMAX_MID - nemis_kprofile + 1
   real, public, allocatable, save, dimension(:,:,:,:) :: &
        gridrcemis      & ! varies every time-step (as ps changes)
       ,gridrcemis0       ! varies every hour
   real, public, allocatable, save, dimension(:,:,:) :: &
        gridrcroadd     & ! Road dust emissions
       ,gridrcroadd0      ! varies every hour

  ! and for budgets (not yet used - not changed dimension)
   real, public,  save, dimension(NSPEC_SHL+1:NSPEC_TOT) ::  totemadd

   integer, private, save :: iemCO  ! index of CO emissions, for debug

   logical, parameter ::USE_OLDSCHEME_ROADDUST=.false. !temporary until the new scheme is validated

contains
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine Emissions(year)


    ! Calls main emission reading routines
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !***********************************************************************
    !   DESCRIPTION:
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
    !   real    snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES)
    !  
    !**********************************************************************

    !--arguments
    integer, intent(in)   :: year        ! Year ( 4-digit)

    !-- local variables
    real    :: conv              ! Conversion factor
    integer :: i, j              ! Loop variables
    real    :: tonne_to_kgm2s    ! Converts tonnes/grid to kg/m2/s
    real    :: ccsum             ! Sum of emissions for one country

    ! arrays for whole EMEP area:
    ! additional arrays on host only for landcode, nlandcode
    ! BIG arrays ... will be needed only on me=0. Make allocatable
    ! to reduce static memory requirements.

    real,    allocatable, dimension(:,:,:,:)  :: globemis 
    integer, allocatable, dimension(:,:)      :: globnland 
    integer, allocatable, dimension(:,:,:)    :: globland  
    real,    allocatable, dimension(:,:,:)    :: globemis_flat
    integer, allocatable, dimension(:,:)      :: flat_globnland 
    integer, allocatable, dimension(:,:,:)    :: flat_globland 
    real,    allocatable, dimension(:,:,:)    :: globroad_dust_pot ! Road dust emission potentials
    integer, allocatable, dimension(:,:)      :: road_globnland 
    integer, allocatable, dimension(:,:,:)    :: road_globland 
    real,    allocatable, dimension(:,:)      :: RoadDustEmis_climate_factor ! Climatic factor for scaling road dust emissions (in TNO model based on yearly average soil water)
    integer :: err1, err2, err3, err4, err5, err6, err7, err8, err9 ! Error messages
    integer :: fic ,insec,inland,iemis 
    integer :: iic,ic,n         ! country codes 
    integer :: isec             ! loop variables: emission sectors
    integer :: iem,iem2         ! loop variable over pollutants (1..NEMIS_FILE)
    integer :: icc,ncc          ! loop variables over  sources
    character(len=40) :: fname  ! File name
    character(len=300) :: inputline 
    real    :: tmpclimfactor

    ! Emission sums (after e_fact adjustments):
    real, dimension(NEMIS_FILE)       :: emsum ! Sum emis over all countries
    real, dimension(NLAND,NEMIS_FILE) :: sumemis, sumemis_local ! Sum of emissions per country
    real, dimension(NEMIS_FILE) :: sumEU ! Sum of emissions in EU

    ! Road dust emission potential sums (just for testing the code, the actual emissions are weather dependent!)
    real, dimension(NROAD_FILES)       :: roaddustsum    ! Sum emission potential over all countries
    real, dimension(NLAND,NROAD_FILES) :: sumroaddust    ! Sum of emission potentials per country
    real, dimension(NLAND,NROAD_FILES) :: sumroaddust_local    ! Sum of emission potentials per country in subdomain
    real :: fractions(MAXLIMAX,MAXLJMAX,NCMAX),SMI(MAXLIMAX,MAXLJMAX),SMI_roadfactor
    logical ::SMI_defined=.false.
    logical :: my_first_call=.true.  ! Used for femis call
    real, allocatable ::emis_tot(:,:)
    character(len=40) :: varname 

    if (MasterProc) write(6,*) "Reading emissions for year",  year

    ! 0) set molwts, conversion factors (e.g. tonne NO2 -> tonne N), and
    !    emission indices (IQSO2=.., )

    allocate(nlandcode(MAXLIMAX,MAXLJMAX),landcode(MAXLIMAX,MAXLJMAX,NCMAX))
    nlandcode=0
    landcode=0
    allocate(flat_nlandcode(MAXLIMAX,MAXLJMAX),flat_landcode(MAXLIMAX,MAXLJMAX,FNCMAX))
    flat_nlandcode=0
    flat_landcode=0
    allocate(road_nlandcode(MAXLIMAX,MAXLJMAX),road_landcode(MAXLIMAX,MAXLJMAX,NCMAX))
    road_nlandcode=0
    road_landcode=0
    allocate(snapemis(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILE))
    snapemis=0.0
    allocate(snapemis_flat(MAXLIMAX,MAXLJMAX,FNCMAX,NEMIS_FILE))
    snapemis_flat=0.0
    allocate(roaddust_emis_pot(MAXLIMAX,MAXLJMAX,NCMAX,NROAD_FILES))
    roaddust_emis_pot=0.0
    allocate(SumSnapEmis(MAXLIMAX,MAXLJMAX,NEMIS_FILE))
    SumSnapEmis=0.0

    !=========================
    call Country_Init()    ! In Country_ml, => NLAND, country codes and 
    !                   names, timezone
     ! The GEA emission data, which is used for EUCAARI runs on the HIRHAM
     ! domain have in several sea grid cells non-zero emissions in other sectors
     ! than SNAP8 and there are also NH3 emission over sea areas. The former 
     ! problem makes the code crash if the sea areas are defined  as 
     ! sea (sea=T), so we treat them as land in the EUCAARI/HIRHAM runs 
     ! (sea=F). This is a problem with GEA emission data only, not the 
     ! HIRHAM domain! When e.g. interpolated EMEP emissions are used on
     !  the HIRHAM domain, this is not a problem.

     if ( SEAFIX_GEA_NEEDED ) then ! Special fix for HIRHAM/GEA

      ! Could have hard-coded 30-34, but best to avoid:
       ic= find_index("BAS", Country(:)%code )
       Country( ic )%is_sea = .false.
       ic= find_index("NOS", Country(:)%code )
       Country( ic )%is_sea = .false.
       ic= find_index("ATL", Country(:)%code )
       Country( ic )%is_sea = .false.
       ic= find_index("MED", Country(:)%code )
       Country( ic )%is_sea = .false.
       ic= find_index("BLS", Country(:)%code )
       Country( ic )%is_sea = .false.

     end if ! HIRHAM/GEA fix

    !=========================

    call consistency_check()               ! Below
    !=========================

    ios = 0

    if ( USE_DEGREEDAY_FACTORS ) &
         call DegreeDayFactors(0)         ! See if we have gridded SNAP-2

    call EmisHeights()     ! vertical emissions profile
    KEMISTOP = KMAX_MID - nemis_kprofile + 1

    if( MasterProc) then   !::::::: ALL READ-INS DONE IN HOST PROCESSOR ::::

       write(6,*) "Reading monthly and daily timefactors"
       !=========================
       call timefactors(year)               ! => fac_emm, fac_edd
       !=========================

    endif

    !=========================
    call EmisSplit()    ! In EmisGet_ml, => emisfrac
    !=========================
    call CheckStop(ios, "ioserror: EmisSplit")


    ! ####################################
    ! Broadcast  monthly and Daily factors (and hourly factors if needed/wanted)
    CALL MPI_BCAST( fac_emm ,8*NLAND*12*NSECTORS*NEMIS_FILE,MPI_BYTE,  0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST( fac_edd ,8*NLAND*7*NSECTORS*NEMIS_FILE,MPI_BYTE,   0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST( fac_ehh24x7 ,8*NSECTORS*24*7,MPI_BYTE,  0,MPI_COMM_WORLD,INFO) 

    !define fac_min for all processors
    do iemis = 1, NEMIS_FILE
       do insec = 1, NSECTORS
          do inland = 1, NLAND
             fac_min(inland,insec,iemis) = minval( fac_emm(inland,:,insec,iemis) )
          enddo
       enddo
    enddo
    if(INERIS_SNAP2 )THEN !  INERIS do not use any base-line for SNAP2
       fac_min(:,ISNAP_DOM,:) = 0.
    end if

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! c4b) Set up DMS factors here - to be used in newmonth
    !      Taken from IQ_DMS=35 for SO2 nature (sector 11)
    !      first_dms_read is true until first call to newmonth finished.


    first_dms_read = .true. 

    ! 4) Read emission files 

    ! allocate for me=0 only:

    err1 = 0
    if ( MasterProc ) then

       if (DEBUG) write(unit=6,fmt=*) "TTT me ", me , "pre-allocate" 
       allocate(globnland(GIMAX,GJMAX),stat=err1)
       allocate(globland(GIMAX,GJMAX,NCMAX),stat=err2)
       allocate(globemis(NSECTORS,GIMAX,GJMAX,NCMAX),stat=err3)

       allocate(flat_globnland(GIMAX,GJMAX),stat=err4)
       allocate(flat_globland(GIMAX,GJMAX,FNCMAX),stat=err5)
       allocate(globemis_flat(GIMAX,GJMAX,FNCMAX),stat=err6)

       call CheckStop(err1, "Allocation error 1 - globland")
       call CheckStop(err2, "Allocation error 2 - globland")
       call CheckStop(err3, "Allocation error 3 - globland")
       call CheckStop(err4, "Allocation error 4 - globland")
       call CheckStop(err5, "Allocation error 5 - globland")
       call CheckStop(err6, "Allocation error 6 - globland")

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
       globnland(:,:) = 0      
       flat_globnland(:,:)=0  
       globland(:,:,:) = 0    
       globemis(:,:,:,:) = 0  
       flat_globland(:,:,:)=0 
       globemis_flat(:,:,:) =0

       if(USE_ROADDUST)then
          road_globnland(:,:)=0
          road_globland(:,:,:)=0
          globroad_dust_pot(:,:,:)=0.
          RoadDustEmis_climate_factor(:,:)=1.0 ! default, no scaling
       endif ! road dust

    end if

    ! Get emission scaling factors
    !>============================

    if ( my_first_call ) then
         sumemis(:,:) =  0.0       ! initialize sums
         ios = 0
         call femis()              ! emission factors (femis.dat file)
         if ( ios /= 0 )return
         my_first_call = .false.
    endif
    !>============================


    do iem = 1, NEMIS_FILE

       !TESTING NEW SYSTEM, no effect on results so far
       ! NB SLOW!
       if( EMIS_TEST == "CdfSnap" ) then
       !if ( iem < 3) then

          ! *************************
           fname="GriddedSnapEmis_" // trim(EMIS_FILE(iem)) // ".nc"
           call EmisGetCdf(iem,fname,incl=EUMACC2) 
          ! *************************
           fname="GlobalSnapEmis_" // trim(EMIS_FILE(iem)) // ".nc"
           call EmisGetCdf(iem,fname,excl=EUMACC2) 
          ! *************************
 
 
           if( debug_proc ) then
             ncc =  nGridEmisCodes(debug_li, debug_lj)
             do isec = 1, 11
               write(*,"(a,2i4,4i3,2i4,9f12.4)")"CDF out  "// &
                 trim(EMIS_FILE(iem)), i_fdom(debug_li), j_fdom(debug_lj), &
                  (GridEmisCodes(debug_li, debug_lj,icc),icc=1,4), &
                  ncc, isec, (GridEmis(isec,debug_li, debug_lj,i,iem), i=1,4) 
             end do
           end if ! debug
         !end if ! iem
       end if ! EMIS_TEST

       ! now again test for me=0
       if(EMIS_SOURCE == "emislist" )then
          if ( MasterProc ) then
             
             ! read in global emissions for one pollutant
             ! *****************
             call EmisGet(iem,EMIS_FILE(iem),IRUNBEG,JRUNBEG,GIMAX,GJMAX, &
               globemis,globnland,globland,sumemis,&
               globemis_flat,flat_globnland,flat_globland)
             ! *****************
             
             
             emsum(iem) = sum( globemis(:,:,:,:) ) + &
                  sum( globemis_flat(:,:,:) )    ! hf
          endif  ! MasterProc

          call CheckStop(ios, "ios error: EmisGet")
          
          ! Send data to processors 
          ! as  e.g. snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,iem)
          ! send to nodes
          
          call global2local(globemis,snapemis(1,1,1,1,iem),MSG_READ1,   &
               NSECTORS,GIMAX,GJMAX,NCMAX,1,1)
          
          call global2local(globemis_flat,snapemis_flat(1,1,1,iem),MSG_READ1,   &
               1,GIMAX,GJMAX,FNCMAX,1,1)

    !   real    snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES)
    !           GridEmis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILE)


          if ( EMIS_TEST == "CdfSnap" .and. debug_proc ) then
             i=debug_li; j=debug_lj
             do isec = 1, 11
               if( snapemis(isec,i,j,1,iem)+GridEmis(isec,i,j,1,iem) > 1.0e-10) &
                  write(*,"(a,i3,2es10.3)") "CDFCOMP", isec, &
                    snapemis(isec,i,j,1,iem) ,GridEmis(isec,i,j,1,iem)
             end do
          end if
       else if ( EMIS_SOURCE == "CdfFractions" ) then
          !use grid independent netcdf emission file
          !experimental so far. Needs a lot of reorganization
          if(.not.allocated(emis_tot))then
             allocate(emis_tot(MAXLIMAX,MAXLJMAX))          
          endif
77        format(A,I2.2)
          sumemis_local(:,iem)=0.0
          do   isec=1,NSECTORS
             
             if(iem==1)varname='SOx_sec'
             if(iem==2)varname='NOx_sec'
             if(iem==3)varname='CO_sec'
             if(iem==4)varname='NMVOC_sec'
             if(iem==5)varname='NH3_sec'
             if(iem==6)varname='PM25_sec'
             if(iem==7)varname='PMco_sec'
             write(varname,77)trim(varname),isec
             call ReadField_CDF('/global/work/mifapw/emep/Data/Emis_TNO7.nc',varname,emis_tot(1,1),&
                  nstart=1,interpol='mass_conservative', &
                  fractions_out=fractions,CC_out=landcode,Ncc_out=nlandcode,needed=.true.,debug_flag=.true.,Undef=0.0)

             if(debug_proc) write(*,*) "CDFTNO ", me, iem, isec, trim(varname)
             do j=1,ljmax
                do i=1,limax                   
                   if(nlandcode(i,j)>NCMAX)then
                      write(*,*)"To many emitter countries in one gridcell: ", me,i,j,nlandcode(i,j)
                      call StopAll("To many countries in one gridcell ")
                   endif
                   do n=1,nlandcode(i,j)
                      snapemis(isec,i,j,n,iem)=snapemis(isec,i,j,n,iem)+fractions(i,j,n)*emis_tot(i,j)                      
                      ic=1
                      if(landcode(i,j,n)<=NLAND)ic=landcode(i,j,n)!most country_index are equal to country_code
!TESTE
               if ( debug_proc .and. i==debug_li .and. j==debug_lj ) then !  .and. iem == iemCO ) then
                   write(*,"(a,2i3,4es10.3)") "FracCdf:" // trim(EMIS_FILE(iem)), &
                     isec, n, emis_tot(i,j), fractions(i,j,n), fractions(i,j,n)*emis_tot(i,j)
                end if
!TESTE

                      if(Country(ic)%index/=landcode(i,j,n))then
                         !if not equal, find which index correspond to country_code
                         do ic=1,NLAND
                            if((Country(ic)%index==landcode(i,j,n)))exit
                         enddo
                         if(ic>NLAND)then
                            write(*,*)"COUNTRY CODE NOT RECOGNIZED OR UNDEFINED: ", landcode(i,j,n)
                            call StopAll("COUNTRY CODE NOT RECOGNIZED ")
                         endif
                      endif                      
                      sumemis_local(ic,iem)= sumemis_local(ic,iem)+0.001*snapemis(isec,i,j,n,iem)!for diagnostics, mass balance
                   enddo
                enddo
             enddo


          enddo!sectors
          CALL MPI_REDUCE(sumemis_local(1,iem),sumemis(1,iem),NLAND,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,INFO) 
          emsum(iem)= sum(sumemis(:,iem))

          if(EMIS_OUT)call EmisOut ("Frac", iem, nlandcode, landcode, snapemis(:,:,:,:,iem) ) ! cdfemis
       else
          call CheckStop("EMIS_SOURCE not set"//trim(EMIS_SOURCE))
       endif

          
       end do ! iem = 1, NEMIS_FILE-loop

    if(USE_ROADDUST) then

       if(USE_OLDSCHEME_ROADDUST)then
          !use scheme with ASCII and grid dependent input data

          ! First Temporary/Test handling of climate factors. Read from file. Should be enough to do this on MasterProc

          if ( MasterProc )then
             !      if(.true.)then
             if (DEBUG_ROADDUST) write(unit=6,fmt=*) "Setting road dust climate scaling factors from", &
                  trim(roaddust_climate_file)
             fname = roaddust_climate_file
             call open_file(IO_EMIS,"r",fname,needed=.true.)
             call CheckStop(ios,"RoadDust climate file: ios error in emission file")

             read(unit=IO_EMIS,fmt="(a200)",iostat=ios) inputline 
             if( inputline(1:1) .ne. "#" ) then ! Is a  comment
                write(*,*)'ERROR in road dust climate factor file!'
                write(*,*)'First line should be a comment line, starting with #'
             else 
                IF(DEBUG_ROADDUST) write(*,*)'I read the comment line:',inputline
             endif

             READCLIMATEFACTOR: do   ! ************* Loop over emislist files *******************

                read(unit=IO_EMIS,fmt=*,iostat=ios) i,j, tmpclimfactor

                if ( ios <  0 ) exit READCLIMATEFACTOR  ! End of file
                call CheckStop(ios > 0,"RoadDust climate file: ios error2 in climate factor file")

                i = i-IRUNBEG+1     ! for RESTRICTED domain
                j = j-JRUNBEG+1     ! for RESTRICTED domain

                if ( i  <=  0 .or. i  >  GIMAX .or.   &
                     j  <=  0 .or. j  >  GJMAX   )&
                     cycle READCLIMATEFACTOR

                RoadDustEmis_climate_factor(i,j) = tmpclimfactor


                if( DEBUG_ROADDUST .and. i==DEBUG_i .and. j==DEBUG_j ) write(*,*) &
                     "DEBUG RoadDust climate factor (read from file)", RoadDustEmis_climate_factor(i,j)

             enddo READCLIMATEFACTOR
             !      else
             !         write(unit=6,fmt=*) "Test run set road dust climate factor to 1 everywhere!"
             !         RoadDustEmis_climate_factor(:,:) = 1.0
             !      endif

          endif !MasterProc

          do iem = 1, NROAD_FILES
             ! now again test for me=0
             if ( MasterProc ) then

                ! read in road dust emission potentials from one file
                ! There will be two road dust files; one for highways(plus some extras) and one for non-highways
                ! Each file contains spring (Mar-May) and rest-of-the-year (June-February) emission potentials
                ! However, this will be changed so that only "rest-of-the-year" data are read and the spring scaling
                ! for the Nordic countries is handled in the EmisSet routine!
                ! Emission potentials are for PM10 so should be split into PM-fine (10%) and PM-coarse (90%)
                ! There should also be some modifications to take into account temporal variations due to different
                ! traffic intensities and differences due to surface wetness (handled in EmisSet) 
                ! and a climatological factor (handled here).

                ! *****************
                call RoadDustGet(iem,ROAD_FILE(iem),IRUNBEG,JRUNBEG,GIMAX,GJMAX, &
                     sumroaddust,&
                     globroad_dust_pot,road_globnland,road_globland)
                ! *****************

                roaddustsum(iem) = sum( globroad_dust_pot(:,:,:) )    ! 

                ! ToDo-2012-0913

                do i=1,GIMAX
                   do j=1,GJMAX
                      if(DEBUG_ROADDUST)then
                         WRITE(*,*)"i,j,RDECF:",i,j,RoadDustEmis_climate_factor(i,j)
                      endif
                      globroad_dust_pot(i,j,:)=RoadDustEmis_climate_factor(i,j)*globroad_dust_pot(i,j,:)
                   enddo
                enddo


             endif  ! MasterProc

             call CheckStop(ios, "ios error: RoadDustGet")

             ! Send data to processors
             ! as  e.g. snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,iem)
             ! send to nodes

             call global2local(globroad_dust_pot,roaddust_emis_pot(1,1,1,iem),MSG_READ1,   &
                  1,GIMAX,GJMAX,NCMAX,1,1)

          end do ! iem = 1, NROAD_FILES-loop

          
          call global2local_int(road_globnland,road_nlandcode,326,&
               GIMAX,GJMAX,1,1,1)  !extra array
          call global2local_int(road_globland,road_landcode,326,&
               GIMAX,GJMAX,NCMAX,1,1)


       else 
          !Use grid-independent Netcdf input files

          do iem = 1, NROAD_FILES
             !Read data from NetCDF file
             if(iem==1)varname='HighwayRoadDustPM10_Jun-Feb'
             if(iem==2)varname='nonHighwayRoadDustPM10_Jun-Feb'
             call CheckStop(iem>2, "TOO MANY ROADFILES")

             roaddust_emis_pot(:,:,:,iem)=0.0
             call ReadField_CDF('RoadMap.nc',varname,roaddust_emis_pot(1,1,1,iem),&
                  nstart=1,interpol='mass_conservative', &
                  fractions_out=fractions,CC_out=road_landcode,Ncc_out=road_nlandcode,needed=.true.,debug_flag=.true.,Undef=0.0)
             if(.not.SMI_defined)then
                varname='SMI1'
                call ReadField_CDF('AVG_SMI_2005_2010.nc',varname,SMI,nstart=1,&
                     interpol='conservative',needed=.true.,debug_flag=.true.)
                SMI_defined=.true.
             endif

             do i=1,LIMAX
                do j=1,LJMAX
                   SMI_roadfactor=3.325-(min(1.0,max(0.5,SMI(i,j)))-0.5) *2*(3.325-1.0)!Peter: Rough estimate to get something varying between 3.325 (SMI<0.5) and 1.0 (SMI>1)
                      if(DEBUG_ROADDUST)then
                       !  WRITE(*,*)"i,j,RDECF:",i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,SMI_roadfactor
                      endif
                   do iic=road_nlandcode(i,j),1,-1
                      roaddust_emis_pot(i,j,iic,iem)=roaddust_emis_pot(i,j,1,iem)*fractions(i,j,iic)*SMI_roadfactor
                   enddo
                enddo
             enddo
             sumroaddust_local(:,iem)=0.0
             do i=1,LIMAX
                do j=1,LJMAX
                   do iic=1,road_nlandcode(i,j)
                      ic=1
                      if(road_landcode(i,j,iic)<=NLAND)ic=road_landcode(i,j,iic)!most country_index are equal to country_code
                      if(Country(ic)%index/=road_landcode(i,j,iic))then
                         !if not equal, find which index correspond to country_code
                         do ic=1,NLAND
                            if((Country(ic)%index==road_landcode(i,j,iic)))exit
                         enddo
                         if(ic>NLAND)then
                            write(*,*)"COUNTRY CODE NOT RECOGNIZED OR UNDEFINED: ", road_landcode(i,j,iic)
                            call StopAll("COUNTRY CODE NOT RECOGNIZED ")

                         endif
                      endif

                      sumroaddust_local(ic,iem)=sumroaddust_local(ic,iem)+0.001*roaddust_emis_pot(i,j,iic,iem)

                   enddo
                enddo
             enddo

          end do ! iem = 1, NROAD_FILES-loop

          sumroaddust=0.0
          CALL MPI_REDUCE(sumroaddust_local,sumroaddust,NLAND*NROAD_FILES,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,INFO) 

       endif

    end if !USE_ROADDUST

    if ( MasterProc ) then
       write(unit=6,fmt=*) "Total emissions by countries:"
       write(unit=IO_LOG,fmt=*) "Total emissions by countries:"
       write(unit=6,fmt="(2a4,3x,30a12)")  "  N "," CC ",(EMIS_FILE(iem),iem=1,NEMIS_FILE)
       write(unit=IO_LOG,fmt="(2a4,3x,30a12)") "  N "," CC ",(EMIS_FILE(iem),iem=1,NEMIS_FILE)

       sumEU(:) = 0.0
       do ic = 1, NLAND
          ccsum = sum( sumemis(ic,:) )
          !if ( ccsum > 0.0 ) then
          if ( ccsum > 0.0 .or. sum( sumcdfemis(ic,:)) > 0.0 ) then
            if ( EMIS_TEST == "None" ) then
             write(unit=6,fmt="(i3,1x,a4,3x,30f12.2)") &
                  ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS_FILE)
            else
             
             write(unit=6,fmt="(a,i3,1x,a4,3x,30f12.2)") &
                  "ORIG:",ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS_FILE)
             write(unit=6,fmt="(a,i3,1x,a4,3x,30f12.2)") & !cdfemis
                  "CDFS:",ic, Country(ic)%code, (sumcdfemis(ic,i),i=1,NEMIS_FILE)
            end if
             write(unit=IO_LOG,fmt="(i3,1x,a4,3x,30f12.2)")& 
                  ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS_FILE)
             if( find_index( Country(ic)%code , EU27(:) ) > 0 ) then
                 sumEU = sumEU + sumemis(ic,:)
             end if
          end if
       end do
       write(unit=6     ,fmt="(i3,1x,a4,3x,30f12.2)") 0, "EU", (sumEU(i),i=1,NEMIS_FILE)
       write(unit=IO_LOG,fmt="(i3,1x,a4,3x,30f12.2)") 0, "EU", (sumEU(i),i=1,NEMIS_FILE)
      !GP_work, TMP solution with hard-coded emissions
       iem= find_index( "nox", EMIS_FILE(:) )
       iem2=find_index( "nh3", EMIS_FILE(:) )
       Ndep_trends = ( sumEU(iem) + sumEU(iem2) ) / 15167.48
       write(unit=6,fmt="(a,2i3,f12.2,f8.4)")  "Ndep_trends", &
            iem, iem2, sumEU(iem) + sumEU(iem2), Ndep_trends 
       write(unit=IO_LOG,fmt="(a,f8.4)")  "Ndep_trends", Ndep_trends 
   

       if(USE_ROADDUST)THEN

          write(unit=6,fmt=*) "Total road dust emission potentials by countries (before precipitation and land corrections):"
          write(unit=IO_LOG,fmt=*) "Total road dust emission potentials by countries (before precipitation and land corrections):"
          write(unit=6,fmt="(2a4,11x,30a12)")  "  N "," CC ",(ROAD_FILE(iem),iem=1,NROAD_FILES)
          write(unit=IO_LOG,fmt="(2a4,11x,30a12)") "  N "," CC ",(ROAD_FILE(iem),iem=1,NROAD_FILES)

          do ic = 1, NLAND
             ccsum = sum( sumroaddust(ic,:) )
             if ( ccsum > 0.0 ) then
                write(unit=6,fmt="(i3,1x,a4,3x,30f12.2)") &
                     ic, Country(ic)%code, (sumroaddust(ic,i),i=1,NROAD_FILES)
                write(unit=IO_LOG,fmt="(i3,1x,a4,3x,30f12.2)")& 
                     ic, Country(ic)%code, (sumroaddust(ic,i),i=1,NROAD_FILES)
             end if
          end do
       endif ! ROAD DUST
    end if

    ! now all values are read, snapemis is distributed, globnland and 
    ! globland are ready for distribution
    ! print *, "calling glob2local_int for iem", iem, " me ", me

    if( EMIS_SOURCE == "emislist" ) then
      call global2local_int(globnland,nlandcode,326, GIMAX,GJMAX,1,1,1)
      call global2local_int(globland, landcode ,326, GIMAX,GJMAX,NCMAX,1,1)

      call global2local_int(flat_globnland,flat_nlandcode,326,&
         GIMAX,GJMAX,1,1,1)  !extra array
      call global2local_int(flat_globland,flat_landcode,326,&
         GIMAX,GJMAX,FNCMAX,1,1)
    else
       !EMIS_SOURCE = CdfFractions
       !emissions directly defined into 
       !nlandcode,landcode and snapemis
    endif

! Create emislist-type files for both snap emissions and Cdf
! Useful for export to other codes, including production of
! new emislist for current NWP grid.
    do iem = 1, NEMIS_FILE

      if( EMIS_SOURCE == "emislist"  ) then
        if(EMIS_OUT)call EmisOut ("Snap", iem, nlandcode, landcode, snapemis(:,:,:,:,iem) )
      end if

      if( EMIS_TEST  == "CdfSnap"  ) then
        write(6, "(a,i3,2i4,2es12.3)") "CALLED CDF PRE "//&
          trim(EMIS_FILE(iem)), me, &
          maxval( nlandcode ),maxval( nGridEmisCodes ), &
          maxval( snapemis(:,:,:,:,iem)), maxval( GridEmis(:,:,:,:,iem))
         if(EMIS_OUT)call EmisOut ("Cdf", iem, nGridEmisCodes, &
                         GridEmisCodes, GridEmis(:,:,:,:,iem) )
      end if
    end do

    !-------------------------------------------------------------------------
    ! Broadcast volcanoe info derived in EmisGet 

    CALL MPI_BCAST(nvolc,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST(i_volc,4*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST(j_volc,4*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST(emis_volc,8*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
 
    ! And N-dep ttrends
    CALL MPI_BCAST(Ndep_trends,8,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 

    !    Conversions 
    !
    !    The emission-data file are so far in units of 
    !    tonnes per grid-square. The conversion factor from tonnes per 50*50km2
    !    annual emission values to surface flux (kg/m2/s) is found by division
    !    with (nydays*24*60*60)s and (h*h)m2 and multiply by 1.e+3.
    !    The conversion factor then equals 1.27e-14

    tonne_to_kgm2s  = 1.0e3 / (nydays * 24.0 * 3600.0 * &
         GRIDWIDTH_M * GRIDWIDTH_M)

    if ( DEBUG .and. MasterProc ) then
       write(unit=6,fmt=*) "CONV:me, nydays, gridwidth = ",me,nydays,GRIDWIDTH_M
       write(unit=6,fmt=*) "No. days in Emissions: ", nydays
       write(unit=6,fmt=*) "tonne_to_kgm2s in Emissions: ", tonne_to_kgm2s
       write(unit=6,fmt=*) "Emissions sums:"
       do iem = 1, NEMIS_FILE
          write(unit=6,fmt="(a15,f12.2)") EMIS_FILE(iem),emsum(iem)
       end do
    endif


    do iem = 1, NEMIS_FILE
       conv = tonne_to_kgm2s
       if ( trim(EMIS_FILE(iem)) == "co" )  iemCO = iem  ! save this index

       if ( DEBUG .and.  debug_proc .and. iem == iemCO ) then
          write(*,"(a,2es10.3)") "SnapPre:" // trim(EMIS_FILE(iem)), &
               sum( snapemis (:,debug_li,debug_lj,:,iem) ) &
               ,sum( snapemis_flat (debug_li,debug_lj,:,iem) )
       end if

       forall (ic=1:NCMAX, j=1:ljmax, i=1:limax, isec=1:NSECTORS)
          snapemis (isec,i,j,ic,iem) = &
               snapemis (isec,i,j,ic,iem) * conv * xm2(i,j)
       end forall

       forall (fic=1:FNCMAX, j=1:ljmax, i=1:limax)
          snapemis_flat(i,j,fic,iem) = &
               snapemis_flat(i,j,fic,iem) * conv * xm2(i,j)
       end forall


       if ( DEBUG .and.  debug_proc .and. iem == iemCO ) then
          write(*,"(a,2es10.3)") "SnapPos:" // trim(EMIS_FILE(iem)), &
               sum( snapemis (:,debug_li,debug_lj,:,iem) ) &
               ,sum( snapemis_flat (debug_li,debug_lj,:,iem) )
       end if
    enddo !iem

    if(USE_ROADDUST)THEN
       do iem = 1, NROAD_FILES
          conv = tonne_to_kgm2s

          forall (ic=1:NCMAX, j=1:ljmax, i=1:limax)
             roaddust_emis_pot(i,j,ic,iem) = &
                  roaddust_emis_pot(i,j,ic,iem) * conv * xm2(i,j)
          end forall
       enddo ! iem
    endif !road dust


    !    if ( VOLCANOES ) then

    ! Read  Volcanos.dat or VolcanoesLL.dat to get volcano height 
    ! and magnitude in the case of VolcanoesLL.dat
    call VolcGet(height_volc)

    !    endif ! VOLCANOES

    err1 = 0
    if ( MasterProc ) then
       deallocate(globnland,stat=err1)
       deallocate(globland ,stat=err2)
       deallocate(globemis ,stat=err3)

       deallocate(flat_globnland,stat=err4)
       deallocate(flat_globland,stat=err5)
       deallocate(globemis_flat,stat=err6)

       if(USE_ROADDUST)THEN
          deallocate(road_globnland,stat=err7)
          deallocate(road_globland,stat=err8)
          deallocate(globroad_dust_pot,stat=err9)
       endif

       call CheckStop(err1, "De-Allocation error 1 - globland") 
       call CheckStop(err2, "De-Allocation error 2 - globland")
       call CheckStop(err3, "De-Allocation error 3 - globland")
       call CheckStop(err4, "De-Allocation error 4 - globland")
       call CheckStop(err5, "De-Allocation error 5 - globland")
       call CheckStop(err6, "De-Allocation error 6 - globland")
       if(USE_ROADDUST)THEN
          call CheckStop(err7, "De-Allocation error 7 - roadglob")
          call CheckStop(err8, "De-Allocation error 8 - roadglob")
          call CheckStop(err9, "De-Allocation error 9 - roadglob")
       endif
    end if

    ! now we have nrecmis and can allocate for gridrcemis:
    ! print *, "ALLOCATING GRIDRC", me, NRCEMIS
    allocate(gridrcemis(NRCEMIS,KEMISTOP:KMAX_MID,MAXLIMAX,MAXLJMAX),stat=err1)
    allocate(gridrcemis0(NRCEMIS,KEMISTOP:KMAX_MID,MAXLIMAX,MAXLJMAX),stat=err2)
    call CheckStop(err1, "Allocation error 1 - gridrcemis") 
    call CheckStop(err2, "Allocation error 2 - gridrcemis0")
    if(USE_ROADDUST)THEN
       allocate(gridrcroadd(NROADDUST,MAXLIMAX,MAXLJMAX),stat=err3)
       allocate(gridrcroadd0(NROADDUST,MAXLIMAX,MAXLJMAX),stat=err4)
       call CheckStop(err3, "Allocation error 3 - gridrcroadd")
       call CheckStop(err4, "Allocation error 4 - gridrcroadd0")
    endif

  end subroutine Emissions

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine consistency_check()
  !------------------------------------------------------------------!
  !    checks that all the values given so far are consistent        !
  !------------------------------------------------------------------!

  character(len=30) :: errormsg
  
  errormsg = "ok"
  if ( size(EMIS_FILE) /= NEMIS_FILE    ) errormsg = " size EMISNAME wrong "

  call CheckStop(errormsg,"Failed consistency check")

 end subroutine consistency_check
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  !*****************************************************************

  subroutine EmisSet(indate)   !  emission re-set every time-step/hour

  !
  !***********************************************************************
  !   DESCRIPTION:
  !   Calculates the emmision-tendencies and the local (instantaneous) dry 
  !   deposition in the emission squares.
  !   emis set once per hour to allow for day/night variation (and voc 
  !   speciation) (based on local time)  for each snap sector.
  !   gridrcemis0 calculated every time-step to allow for ps changes.
  !   inputs from Emissions in EMISSIONS_ML:
  !   country and snap-specific array : 
  !          snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES) 
  !  
  !   Units:
  !   snapemis has units of kg/m2/s, SO2 as S, NO2 as N, NH3 as N. 
  !   Map factor (xm2) already accounted for. 
  !  
  !   Data on how many countries contribute per grid square is stored in
  !   nlandcode(MAXLIMAX,MAXLJMAX) , with the country-codes given by
  !   landcode(MAXLIMAX,MAXLJMAX,NCMAX).
  !     
  !   Monthly and weekday factors are pre-multiplied and stored in:
  !       real timefac(NLAND,NSECTORS,NEMIS_FILES)
  !   And day-hour factors in fac_ehh24x7
  !
  !*************************************************************************

  implicit none
  type(date), intent(in) :: indate           ! Gives year..seconds
  integer :: i, j, k, f                   ! cooridnates, loop variables
  integer :: icc, ncc                        ! No. of countries in grid.

  integer :: ficc,fncc                       ! No. of countries with
  integer :: iqrc                            ! emis indices 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILE
  integer :: itot             ! index in xn()

  ! Save daytime value between calls, initialise to zero
  integer, save, dimension(NLAND) ::  daytime = 0  !  0=night, 1=day
  integer, save, dimension(NLAND) ::  localhour = 1  ! 1-24 local hour in the different countries, ? How to handle Russia, with multiple timezones???
  integer                         ::  hourloc      !  local hour 
  logical                         ::  hourchange   !             
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions

  real ::  ehlpcom,ehlpcom0
  real ::  tfac, dtgrid    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s               ! source term (emis) before splitting
  integer :: iland, iland_timefac  ! country codes, and codes for timefac 

  real ::  sf              ! source term (emis) before splitting 
                           ! (for flat emissions)
  integer :: flat_iland    ! country codes (countries with flat emissions)

  integer, save :: oldday = -1, oldhour = -1
  integer, save :: wday , wday_loc ! wday = day of the week 1-7
  real ::  oldtfac
  logical :: debug_tfac, debug_kprof

! If timezone=-100, calculate daytime based on longitude rather than timezone
  integer :: daytime_longitude, daytime_iland, hour_longitude, hour_iland
 
! Initialize
  ehlpcom0 = GRAV* 0.001*AVOG!/ (sigma_bnd(k+1) - sigma_bnd(k))

  ! Scaling for totemadd:
     dtgrid = dt_advec * GRIDWIDTH_M * GRIDWIDTH_M 

   ! The emis array only needs to be updated every full hour. The 
   ! time-factor calculation needs to know if a local-time causes a shift 
   ! from day to night.  In addition, we reset an overall day's time-factors
   ! at midnight every day. 

     hourchange = .false.
     if ( indate%hour /= oldhour .or. indate%day /= oldday ) then

           hourchange = .true.
           oldhour = indate%hour

           if ( indate%day /= oldday  )then

              !==========================
               call NewDayFactors(indate)
               if ( USE_DEGREEDAY_FACTORS) & 
                 call DegreeDayFactors(daynumber) ! => fac_emm, fac_edd
              !==========================

               ! for ROADDUST
               wday=day_of_week(indate%year,indate%month,indate%day)
               if(wday==0)wday=7 ! Sunday -> 7

               oldday = indate%day
           endif
     end if

     if( DEBUG_EMISTIMEFACS .and. MasterProc ) then
       write(*,"(a,2f8.3)") " EmisSet  traffic 24x7", &
           fac_ehh24x7(ISNAP_TRAF,1,4),fac_ehh24x7(ISNAP_TRAF,13,4)
     end if


    !..........................................
    !  Look for day-night changes, after local time correction
    !  (daytime(iland) not used if  LONGITUDE_TIME=true)

    do iland = 1, NLAND

       daytime(iland) = 0
       hourloc        = indate%hour + Country(iland)%timezone

       localhour(iland) = hourloc  ! here from 0 to 23

       if ( hourloc  >=   7 .and.  hourloc <= 18 ) daytime(iland)=1

    end do ! iland


    if ( hourchange ) then 

         totemadd(:)  = 0.
         gridrcemis0(:,:,:,:) = 0.0 
         SumSnapEmis(:,:,:) = 0.0
         if(USE_ROADDUST)gridrcroadd0(:,:,:) = 0.0

        !..........................................
        ! Process each grid:

        do j = 1,ljmax
            do i = 1,limax

               ncc = nlandcode(i,j)            ! No. of countries in grid
               debug_tfac =  ( DEBUG_EMISTIMEFACS .and. debug_proc .and. &
                                 i==DEBUG_li .and. j==DEBUG_lj )

                ! find the approximate local time:
                  hourloc= mod(nint(indate%hour+24*(1+glon(i,j)/360.0)),24)
                  hour_longitude=hourloc
                  daytime_longitude=0
                  if( hourloc>=7.and.hourloc<= 18) daytime_longitude=1
             
              !*************************************************
              ! First loop over non-flat (one sector) emissions
              !*************************************************

              tmpemis(:)=0.
              do icc = 1, ncc
                  iland = landcode(i,j,icc)     ! 1=Albania, etc.
                  iland_timefac = Country(iland)%timefac_index

                if(Country(iland)%timezone==-100)then
                   daytime_iland=daytime_longitude
                   hour_iland=hour_longitude + 1   ! add 1 to get 1..24 
                else
                   daytime_iland=daytime(iland)
                   hour_iland=localhour(iland) + 1
                endif
                !if( hour_iland > 24 ) hour_iland = 1 !DSA12
                wday_loc=wday 
                if( hour_iland > 24 ) then
                   hour_iland = hour_iland - 24
                   wday_loc=wday + 1
                   if(wday_loc==0)wday_loc=7 ! Sunday -> 7
                   if(wday_loc>7 )wday_loc=1 
                end if

                call CheckStop( hour_iland < 1, &
                    "ERROR: HOUR Zero in EmisSet")

                if( debug_tfac ) then 
                    write(*,"(a,i4,2i3,i5,2i4,3x,4i3)") "EmisSet DAYS times ", daynumber, &
                     wday, wday_loc, iland, daytime_longitude, daytime_iland,&
                      hour_longitude, hour_iland, hourloc, Country(iland)%timezone
                    call datewrite("EmisSet DAY 24x7:", &
                        (/ icc, iland, wday, wday_loc, hour_iland /), &
                         (/ fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc) /) )
                end if


                 !  As each emission sector has a different diurnal profile
                 !  and possibly speciation, we loop over each sector, adding
                 !  the found emission rates to gridrcemis as we go.
                 !  ==================================================


                do isec = 1, NSECTORS       ! Loop over snap codes

                   ! Calculate emission rates from snapemis, time-factors, 
                   ! and if appropriate any speciation fraction (NEMIS_FRAC)

                   iqrc = 0   ! index over emisfrac

                   do iem = 1, NEMIS_FILE 

                         tfac = timefac(iland_timefac,isec,iem) * &
                              fac_ehh24x7(isec,hour_iland,wday_loc)

                         if( debug_tfac .and. iem == 1 ) then 
                              write(*,"(a,2i4,f8.3)") "EmisSet DAY TFAC:", &
                                 isec, hour_iland, tfac
                         end if

                      !Degree days - only SNAP-2 
                      if ( USE_DEGREEDAY_FACTORS .and. &
                          isec == ISNAP_DOM .and. Gridded_SNAP2_Factors ) then

                        oldtfac = tfac

                        ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
                        ! we make use of a baseload even for SNAP2

                        tfac = ( fac_min(iland,isec,iem) + & ! constant baseload
                                ( 1.0-fac_min(iland,isec,iem) )* gridfac_HDD(i,j) ) &
                                * fac_ehh24x7(isec,hour_iland,wday_loc)

                        if ( debug_tfac .and. indate%hour == 12 .and. iem==1 )  then  ! 
                           write(*,"(a,2i3,2i4,7f8.3)") "SNAPHDD tfac ",  &
                              isec, iland, daynumber, indate%hour, &
                                 timefac(iland_timefac,isec,iem), t2_nwp(i,j,2)-273.15, &
                                   fac_min(iland,isec,iem),  gridfac_HDD(i,j), tfac
                        end if
                      end if ! =============== HDD 

                      s =  tfac * snapemis(isec,i,j,icc,iem)

                    ! prelim emis sum kg/m2/s
                       SumSnapEmis(i,j,iem) = SumSnapEmis(i,j,iem) + s

                      do f = 1, emis_nsplit( iem )
                           iqrc = iqrc + 1
                           itot = iqrc2itot(iqrc)
                           tmpemis(iqrc) = s * emisfrac(iqrc,isec,iland)
                    ! Add up emissions in ktonne 
                           totemadd(itot) = totemadd(itot) + &
                                     tmpemis(iqrc) * dtgrid * xmd(i,j)
                      end do ! f

                   end do ! iem


                   !  Assign to height levels 1-KEMISTOP

                   !do k=KEMISTOP,KMAX_MID
                   do k=KEMISTOP,KMAX_MID
                      do iqrc =1, nrcemis
                         gridrcemis0(iqrc,k,i,j) =   &
                            gridrcemis0(iqrc,k,i,j) + tmpemis(iqrc)*   &
                            ehlpcom0*emis_kprofile(KMAX_BND-k,isec) &
                            !ehlpcom0(k)*VERTFAC(KMAX_BND-k,isec) &
                            * emis_masscorr(iqrc)
                        !if( debug_tfac.and. iqrc==1 ) then 
                        !  write(*,"(a,2i3,2f8.3)") "KPROF ", &
                        !    isec, KMAX_BND-k, &
                        !    VERTFAC(KMAX_BND-k,isec),  &
                        !    emis_kprofile(KMAX_BND-k,isec)
                        !end if
 
                      end do ! iem
                   end do   ! k

                enddo  ! isec
!      ==================================================

       end do ! icc  
 
       !************************************
       ! Then loop over flat emissions
       !************************************
       tmpemis(:)=0.
       fncc = flat_nlandcode(i,j) ! No. of countries with flat 
                                  ! emissions in grid
       
       do ficc = 1, fncc
          flat_iland = flat_landcode(i,j,ficc) ! 30=BAS etc.

          if ( Country(flat_iland)%is_sea ) then  ! saves if statements below
               isec = ISNAP_SHIP 
          else
               isec = ISNAP_NAT
          end if

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

                do f = 1, emis_nsplit( iem )
                   iqrc = iqrc + 1
                   itot = iqrc2itot(iqrc)
                   tmpemis(iqrc) = sf * emisfrac(iqrc,isec,flat_iland)

             ! Add flat emissions in ktonne 
                   totemadd(itot) = totemadd(itot) + &
                             tmpemis(iqrc) * dtgrid * xmd(i,j)

                end do ! f

            end do ! iem

         ! Assign flat emissions to height levels 1-4
         ! Note, no VERTFAC

             do iqrc =1, nrcemis

                gridrcemis0(iqrc,KMAX_MID,i,j) =   &
                  gridrcemis0(iqrc,KMAX_MID,i,j) + tmpemis(iqrc)*&
                    ehlpcom0 * emis_masscorr(iqrc)
             end do ! iem

!      ==================================================
       end do !ficc 

       if(USE_ROADDUST)then
  NO_PRECIP: if(surface_precip(i,j) < 0.1) then ! Limit as in TNO-model (but Lotos/Euros has precip in mm/3h) In the EMEP case this is in mm/h, so should be equivalent with 2.4mm per day

! should use the temporal variation for road dust (SNAP_HOURFAC(HH,7))
! and a weekday factor (initially taken from TNO model, could use country
! dependent factor in the future)

! Temporal variation taken from TNO -> No monthly variation and a single
! weekday and diurnal variation (same for all countries)     
! -> Need to know day_of_week
!    Relatively weak variation with day of week so use a simplified approach
!    

!            if( DEBUG_ROADDUST .and. debug_proc .and. i==DEBUG_li .and. j==DEBUG_lj )THEN
!               write(*,*)"DEBUG ROADDUST! Dry! ncc=", road_nlandcode(i,j)
!            endif

         ncc = road_nlandcode(i,j) ! number of countries in grid point
         do icc = 1, ncc    
            iland = road_landcode(i,j,icc)
            
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
            end if

            if(((icc.eq.IC_FI).or.(icc.eq.IC_NO).or.(icc.eq.IC_SE)).and. & ! Nordic countries
                 ((indate%month.eq.3).or.(indate%month.eq.4).or.(indate%month.eq.5)))then ! spring road dust
               tfac =  fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc) *2.0 ! Doubling in Mar-May (as in TNO model)
            else
               tfac =  fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc)
            endif

            do iem = 1, NROAD_FILES
                s = tfac * roaddust_emis_pot(i,j,icc,iem)
            if( DEBUG_ROADDUST .and. debug_proc .and. i==DEBUG_li .and. j==DEBUG_lj )THEN
               write(*,*)"DEBUG ROADDUST! iem,tfac,icc,roaddust_emis_pot,s", &
                          iem,tfac,icc,roaddust_emis_pot(i,j,icc,iem),s
            endif
 
                gridrcroadd0(QROADDUST_FI,i,j)=gridrcroadd0(QROADDUST_FI,i,j)+ &
                     ROADDUST_FINE_FRAC*s
                gridrcroadd0(QROADDUST_CO,i,j)=gridrcroadd0(QROADDUST_CO,i,j)+ &
                     (1.-ROADDUST_FINE_FRAC)*s

                if ( DEBUG_ROADDUST .AND. debug_proc .and.i==debug_li .and. j==debug_lj)  then  ! 
                   write(*,*)"gridrcroadfine",gridrcroadd0(QROADDUST_FI,i,j)
                   write(*,*)"gridrcroadcoarse",gridrcroadd0(QROADDUST_CO,i,j)
                end if


             enddo ! nroad files
                
          enddo ! icc

! should pick the correct emissions (spring or rest of year)
! and add the emissions from HIGHWAYplus and NONHIGHWAYS,
! using correct fine and coarse fractions.
       else ! precipitation
          gridrcroadd0(:,i,j)=0.
       endif NO_PRECIP
 
  endif ! ROADDUST

   end do ! i
 end do ! j
        if ( DEBUG .and.  debug_proc ) then    ! emis sum kg/m2/s
          call datewrite("SnapSum, kg/m2/s:" // trim(EMIS_FILE(iemCO)), &
                  (/ SumSnapEmis(debug_li,debug_lj,iemCO)  /) )
        end if

    call Set_Volc !set hourly volcano emission(rcemis_volc0)


  end if ! hourchange 


  ! We now scale gridrcemis to get emissions in molecules/cm3/s

   do k= KEMISTOP, KMAX_MID
     do j = 1,ljmax
        do i = 1,limax

              ehlpcom= roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))
!RB: This should also be done for the road dust emissions
              do iqrc =1, NRCEMIS
                 gridrcemis(iqrc,k,i,j) =  &
                    gridrcemis0(iqrc,k,i,j)* ehlpcom
              enddo  ! iqrc
        end do ! i
      end do ! j
   end do ! k

   if(USE_ROADDUST)THEN
      if( DEBUG_ROADDUST .and. debug_proc)then
         write(*,*)"Before the unit scaling", &
              gridrcroadd(1,DEBUG_li,DEBUG_lj), &
              gridrcroadd(2,DEBUG_li,DEBUG_lj)
      endif
      do j = 1,ljmax
         do i = 1,limax
            ehlpcom= roa(i,j,KMAX_MID,1)/(ps(i,j,1)-PT)
              do iqrc =1, NROADDUST
                 gridrcroadd(iqrc,i,j) =  &
                      gridrcroadd0(iqrc,i,j)* ehlpcom * ehlpcom0 * roaddust_masscorr(iqrc)
              enddo  ! iqrc
           end do ! i
        end do ! j
      if( DEBUG_ROADDUST .and. debug_proc)then
         write(*,*)"After the unit scaling", &
              gridrcroadd(1,DEBUG_li,DEBUG_lj), &
              gridrcroadd(2,DEBUG_li,DEBUG_lj)
      endif
    endif

 ! Scale volc emissions to get emissions in molecules/cm3/s (rcemis_volc)

     call Scale_Volc

  end subroutine EmisSet
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    subroutine newmonth

!.....................................................................
!   DESCRIPTION:
!   Reads in natural DMS emissions at start of each month. Update
!   landcode and nlandcode arrays as needed.
!
!   Reads in snow cover at start of each month. 
!
!   April 2010: read monthly aircraft NOx emissions
!
!...........................................................................
      use AirEmis_ml, only : airn
      use ModelConstants_ml, only : KCHEMTOP, KMAX_MID
      use NetCDF_ml, only : ReadField_CDF

    integer i, j,k, iyr
    integer n, flat_ncmaxfound         ! Max. no. countries w/flat emissions
    real :: rdemis(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
    character(len=20) :: fname
    real ktonne_to_kgm2s, tonnemonth_to_kgm2s  ! Units conversion
    integer :: IQSO2                   ! Index of sox in  EMIS_FILE
    integer errcode
    real,    allocatable, dimension(:,:,:,:)  :: globemis 
    integer :: month,iem,ic,iic,isec, err3,icc
    real :: duml,dumh,tmpsec(NSECTORS),conv
    logical , save :: first_call=.true.
    real, dimension(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILE) &
            ::  snapemis_month ! monthly emissions tonne/month
    logical :: needed_found

   ! For now, only the global runs use the Monthly files

        !NML logical, parameter :: MONTHLY_GRIDEMIS= IS_GLOBAL      
        integer :: kstart,kend,nstart,Nyears
        real :: buffer(MAXLIMAX,MAXLJMAX),SumSoilNOx
        
        if(.not.allocated(airn).and.(USE_LIGHTNING_EMIS.or.USE_AIRCRAFT_EMIS))then
           allocate(airn(KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX))
        endif
        
if( USE_AIRCRAFT_EMIS )then
   airn = 0.0
   kstart=KCHEMTOP
   kend=KMAX_MID

   call ReadField_CDF('AircraftEmis_FL.nc','NOx',airn,&
    nstart=current_date%month,kstart=kstart,kend=kend,&
    interpol='mass_conservative', needed=.true.,debug_flag=.false.)

! convert from kg(NO2)/month into molecules/cm3/s
! from kg to molecules: 0.001*AVOG/species(NO2)%molwt
! use roa to find dz for consistency with other emissions 
! (otherwise could have used z_bnd directly)
! dz=dP/(roa*GRAV)  dP=dA(k) + dB(k)*ps(i,j,1)
! dV=dz*dx*dy=dz*gridwidth**2/xm**2 *1e6 (1e6 for m3->cm3)
! from month to seconds: ndaysmonth*24*3600

conv=0.001*AVOG/species(NO2)%molwt*GRAV/gridwidth_m**2*1.0e-6/(nmdays(current_date%month)*24*3600)

do k=KCHEMTOP,KMAX_MID
do j=1,ljmax
do i=1,limax

airn(k,i,j)=airn(k,i,j)*conv*(roa(i,j,k,1))/(dA(k) + dB(k)*ps(i,j,1))*xm2(i,j)

enddo
enddo
enddo

endif

if(USE_EURO_SOILNOX)then  ! European Soil NOx emissions
  if (DEBUG_SOILNOX .and. debug_proc ) write(*,*) "Emissions DEBUG_SOILNOX START"

  ! read in map of annual N-deposition produced from pre-runs of EMEP model
  ! with script mkcdo.annualNdep
  ! 
   call ReadField_CDF('annualNdep.nc',&
     'Ndep_m2',AnnualNdep,1, interpol='zero_order',needed=.true.,debug_flag=.false.,UnDef=0.0)

   if (DEBUG_SOILNOX .and. debug_proc ) then
      write(*,"(a,4es12.3)") "Emissions_ml: SOILNOX AnnualDEBUG ", &
           AnnualNdep(debug_li, debug_lj), maxval(AnnualNdep), minval(AnnualNdep)
   end if
   call CheckStop(USE_GLOBAL_SOILNOX, "SOILNOX - cannot use global with Euro")
   ! We then calculate SoulNOx in Biogenics_ml

else if ( USE_GLOBAL_SOILNOX ) then ! Global soil NOx

  do j=1,ljmax
   do i=1,limax
      SoilNOx(i,j)=0.0      
      buffer(i,j)=0.0      
   enddo
  enddo

  nstart=(current_date%year-1996)*12 + current_date%month
  if(nstart>0.and.nstart<=120)then
   !the month is defined
   call ReadField_CDF('nox_emission_1996-2005.nc','NOX_EMISSION',SoilNOx,nstart=nstart,&
   interpol='conservative',known_projection="lon lat",needed=.true.,debug_flag=.false.)
   if ( DEBUG_SOILNOX .and.debug_proc ) write(*,*) "PROPER YEAR of SOILNO ", current_date%year, nstart
  else
   !the year is not defined; average over all years
   Nyears=10 !10 years defined
   do iyr=1,Nyears 
      nstart=12*(iyr-1) + current_date%month  
      call ReadField_CDF('nox_emission_1996-2005.nc','NOX_EMISSION',buffer,nstart=nstart,&
      interpol='conservative',known_projection="lon lat",needed=.true.,debug_flag=.false.,UnDef=0.0)
      do j=1,ljmax 
         do i=1,limax
           SoilNOx(i,j)=SoilNOx(i,j)+buffer(i,j)
         end do
      end do
      if ( DEBUG_SOILNOX .and.debug_proc ) then
         write(*,"(a,2i6,es10.3,a,2es10.3)") "Averaging SOILNO  inputs", &
           1995+(iyr-1), nstart,SoilNOx(debug_li, debug_lj), &
            "max: ", maxval(buffer), maxval(SoilNOx)
      end if
   enddo
   SoilNOx=SoilNOx/Nyears
  endif ! nstart test

  if ( DEBUG_SOILNOX .and. debug_proc ) then
    write(*,"(a,i3,3es10.3)") "After Global SOILNO ",  me, maxval(SoilNOx), &
           SoilNOx(debug_li, debug_lj)
    !write(*,"(a,i3,3es10.3)") "After Global SOILNO ",  me, maxval(SoilNOx), SoilNOx(3, 3)
  end if 
else ! no soil NO
  if (DEBUG_SOILNOX .and. debug_proc ) write(*,*) "Emissions DEBUG_SOILNOX - none"
end if !  SOIL NO

!for testing, compute total soil NOx emissions within domain
!convert from g/m2/day into kg/day
if ( USE_GLOBAL_SOILNOX ) then 
  SumSoilNOx=0.0
  SoilNOx = max(0.0, SoilNOx)  ! Stops the NEGs!
  do j=1,ljmax
   do i=1,limax      
      SumSoilNOx=SumSoilNOx+0.001*SoilNOx(i,j)*gridwidth_m**2*xmd(i,j)      
   enddo
  enddo
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, SumSoilNOx , 1, &
     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
  if(MasterProc)write(*,*)'GLOBAL SOILNOX emissions this month within domain',&
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
      SoilNOx(i,j)=SoilNOx(i,j)*conv*(roa(i,j,k,1))/(dA(k) + dB(k)*ps(i,j,1))     
   enddo
  enddo

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

        ktonne_to_kgm2s  = 1.0e6 /        &
        (nydays*24.*60.*60.*GRIDWIDTH_M*GRIDWIDTH_M)

    if ( MasterProc .and. DEBUG) then
      write(6,*) 'Enters newmonth, mm, ktonne_to_kgm2s = ',    &
          current_date%month,ktonne_to_kgm2s
      write(6,*) ' first_dms_read = ', first_dms_read
    end if ! me 
!...........................................................................

!...........................................................................
!        DMS Input - land 35 - SNAP sector 11
!...........................................................................
     flat_ncmaxfound = 0  ! Max. no. countries(w/flat emissions) per grid

!  Natural SO2 emissions

          IQSO2 = 0
          do i = 1, NEMIS_FILE
            if ( trim( EMIS_FILE(i) ) == "sox" ) IQSO2 = i
          end do

          if ( IQSO2 < 1 ) then
              write(*,*) " No SO2 emissions - need to skip DMS also"
              return     ! No need to read DMS fields 

          else    
            ! We have so2 emission so need DMS also

             if ( MasterProc ) then
                
                write(fname,fmt='(''natso2'',i2.2,''.dat'')')     &
                     current_date%month
                write(6,*) 'Reading DMS emissions from ',trim(fname)
             endif
             
             needed_found=.false.
             call ReadField(IO_DMS,fname,rdemis,needed_found)
             if(needed_found)then
                errcode = 0
                do j=1,ljmax
                   do i=1,limax
                      
                      ! Add DMS for country code IQ_DMS=35  to snap sector 11=Nature.
                      ! First time we read we must add DMS to the "countries" 
                      ! contributing within the grid square. 
                      
                      ! - for flat emissions:
                      
                      if ( first_dms_read ) then 
                         flat_nlandcode(i,j) = flat_nlandcode(i,j) + 1 
                         n = flat_nlandcode(i,j)
                          flat_landcode(i,j,n) = IQ_DMS   ! country code 35 
                          if ( n > flat_ncmaxfound ) then
                             flat_ncmaxfound = n 
                             if (DEBUG) write(6,*)'DMS Increased flat_ncmaxfound to ',n 
                             call CheckStop( n > FNCMAX, "IncreaseFNCMAX for dms")
                          endif
                       else  ! We know that DMS lies last in the array, so:
                          n = flat_nlandcode(i,j)
                          call CheckStop(flat_landcode(i,j,n), IQ_DMS, &
                               "Newmonth:DMS not last!")
                       endif
                       
                       snapemis_flat(i,j,n,IQSO2) = rdemis(i,j) * ktonne_to_kgm2s &
                            * xm2(i,j)
                    enddo ! i
                 enddo ! j


                 if ( first_dms_read ) then
                    if (DEBUG) write(6,*)'me ',me, ' Increased flat_ncmaxfound to ' &
                         ,flat_ncmaxfound 
                    first_dms_read = .false.
                 end if
              else!no dms file found
                 if ( MasterProc ) write(6,*) 'WARNING: NO DMS emissions found '
                 if ( MasterProc ) write(unit=IO_LOG,fmt=*) "WARNING: NO DMS emissions found "
              end if

        end if  ! IQSO2>0


   if(MONTHLY_GRIDEMIS)then

   !Read monthly emission files

   if(first_call)then
      do j=1,ljmax
         do i=1,limax
            nlandcode(i,j)=nlandcode(i,j)+1
            icc=nlandcode(i,j)
            landcode(i,j,icc)=67 
         enddo
      enddo
      first_call=.false.
   endif
    month = current_date%month
    tonnemonth_to_kgm2s= 1.0e3 /         &
         (nmdays(month)*24.*60.*60.*GRIDWIDTH_M*GRIDWIDTH_M)

    if ( MasterProc ) then
       allocate(globemis(NSECTORS,GIMAX,GJMAX,NCMAX),stat=err3)
       call CheckStop(err3, "Allocation error err3 - globland")
    end if
    do iem = 1, NEMIS_FILE
!      if (trim(EMIS_NAME(iem)).ne.'nox' .and. trim(EMIS_NAME(iem)).ne.'co'.and.&
!           trim(EMIS_NAME(iem)).ne.'pm25'.and.&
!           trim(EMIS_NAME(iem)).ne.'voc'.and.trim(EMIS_NAME(iem)).ne.'nh3'.and.trim(EMIS_NAME(iem)).ne.'sox')cycle !
!      snapemis (:,:,:,:,iem) = 0.0 !NB all previous (snap)emis set to 0
      if ( MasterProc ) then

         globemis = 0.0

         write(fname,fmt='(''grid'',A,i2.2)')     &
              trim(EMIS_FILE(iem))//'.',month
         write(6,*) 'filename for GLOBAL emission',fname
         call open_file(IO_EMIS,"r",fname,needed=.true.)
         call CheckStop( ios , "ios error: emislist" // fname )

         READEMIS: do   ! Loop over emislist files

            read(unit=IO_EMIS,fmt=*,iostat=ios) iic, i, j, duml, dumh, &
                 (tmpsec(isec),isec=1,NSECTORS)
            if ( ios <  0 ) exit READEMIS            ! End of file
            call CheckStop( ios , "GetEmis ios: error on " // fname ) ! exits if ios>0

            ic=1 ! default country



            i = i-IRUNBEG+1     ! for RESTRICTED domain
            j = j-JRUNBEG+1     ! for RESTRICTED domain

            if ( i  <=  0 .or. i  >  GIMAX .or.   & 
                 j  <=  0 .or. j  >  GJMAX .or.   &
                 ic <=  0 .or. ic >  NLAND .or.   &
                 ic == IC_NAT )                   &   !Excludes DMS
                 cycle READEMIS
            globemis(:,i,j,ic) = globemis(:,i,j,ic) & !Put everything in land "1"
                                 + tmpsec(:)

         end do READEMIS
         !
         close(IO_EMIS)
         ios = 0

     endif  ! MasterProc
      call global2local(globemis,snapemis_month(1,1,1,1,iem),MSG_READ1,   &
           NSECTORS,GIMAX,GJMAX,NCMAX,1,1)
   end do ! iem = 1, NEMIS-loop

   ic=1 
   do iem = 1, NEMIS_FILE
!         write(*,*)'iem=',iem
!      if (trim(EMIS_NAME(iem)).ne.'nox' .and. trim(EMIS_NAME(iem)).ne.'co'.and.&
!           trim(EMIS_NAME(iem)).ne.'pm25'.and.&
!           trim(EMIS_NAME(iem)).ne.'voc'.and.trim(EMIS_NAME(iem)).ne.'nh3'.and.trim(EMIS_NAME(iem)).ne.'sox')cycle !
      conv = tonnemonth_to_kgm2s
      do j=1,ljmax
         do i=1,limax
            icc=nlandcode(i,j) !67
           do isec=1,NSECTORS
              snapemis (isec,i,j,icc,iem) = &
                    snapemis_month (isec,i,j,ic,iem) * conv * xm2(i,j)
            enddo
         enddo
     enddo
   enddo !iem
   if ( MasterProc ) then
      deallocate(globemis)
   end if
   endif

    end subroutine newmonth

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine EmisOut(label, iem,nsources,sources,emis)

 ! Ascii output of emissions fields (after interpolation and femis.
 ! To print out emissions, we need to get fields for each country
 ! in turn, needing info from GridCodes
 ! e.g. 11-SNAP data is snapemis:
 !  real    snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES)

    character(len=*) :: label
    integer, intent(in) :: iem
    real,intent(in), dimension(:,:,:,:):: emis      ! Emission values
    integer, dimension(:,:), intent(in) :: nsources      ! Emission values
    integer, dimension(:,:,:), intent(in) :: sources      ! Emission values
    real, allocatable, dimension(:,:):: lemis,gemis    ! 2-D emission fields
    real, allocatable, dimension(:,:,:):: locemis      ! Emission values
    real, allocatable, dimension(:,:,:):: globemis      ! Emission values
    character(len=100) :: txt
    real :: low, high
    integer :: msg(4), i,j, ii, jj, isec, icc, ncc, iland


    txt = trim(label)//"."//trim(EMIS_FILE(iem))
    msg(:) = 0

    if(DEBUG) write(*,*)  "CALLED "//trim(txt), me, &
          maxval(emis), maxval(nsources), maxval(sources)

    if ( MasterProc ) then
      allocate(globemis(GIMAX, GJMAX,NSECTORS), stat=msg(1) )
      open(IO_TMP,file="EmisOut"//trim(txt))
    end if

    allocate(locemis(MAXLIMAX, MAXLJMAX,NSECTORS), stat=msg(2) )
    allocate(lemis(MAXLIMAX, MAXLJMAX), stat=msg(3) )
    allocate(gemis(GIMAX, GJMAX), stat=msg(4) )

    call CheckStop( any(msg(:) /= 0),"EmisOut alloc error:"//trim(label))

    ! TEST LOCAL2GLOBAL. Found that the local array is changed
    !lemis = 100 + me
    !print *, "LEMIS-IN ", me, lemis(2,2)
    !call local2global( lemis(:,:), gemis(:,:), msg )
    !if( MasterProc ) then
    ! do j = 1, GJMAX, 3
    !   print *, "LEMIS-UT ", j, lemis(2,j), gemis(2,j)
    ! end do
    !end if

EMLAND: do iland = 1, NLAND
     locemis = 0.0
!      print *,  trim(txt)//" iland ", me, iland, maxval(emis(:,:,:,:))

     !/ Collect emissions locally by country code iland
     do j = 1, ljmax
       do i = 1, limax
         ncc = nsources(i,j)
         do icc = 1, ncc
            if ( sources(i,j,icc) == iland ) then
                locemis(i,j,: ) = emis(:, i,j,icc)
                if ( DEBUG ) &
                 call CheckStop( any( locemis(i,j,:) < 0.0), "NEG LOCEMIS")
            end if
         end do
       end do
      end do ! j

     ! Should never happen, but...
      !call CheckStop( any( lemis < 0.0 ) , "NEG LEMIS")
  

      !/ Transmit to MasterProc
      !/ (Need to be careful, local2global changes local arrays. Hence lemis)

      do isec = 1, NSECTORS
        lemis = locemis(:,:,isec)
        if( DEBUG) write(*,*)  trim(txt)//" lemis ", &
            me, iland, maxval(lemis(:,:))
        call local2global( lemis, gemis, msg )
        if( MasterProc ) globemis(:,:,isec) = gemis(:,:) !! for output
      end do ! isec
          
      if( MasterProc ) then ! output emislist-type file

        do i = 1, GIMAX
         do j = 1, GJMAX
           ii=i+IRUNBEG-1
           jj= j+JRUNBEG-1
           if(sum(globemis(i,j,:)) > 1.0e-10) then
              high = globemis(i,j,1)
              low =  sum(globemis(i,j,2:NSECTORS))
              write(IO_TMP,"(i3,2i4,2x,13es10.3)") iland, ii,jj, &
                low, high, (globemis(i,j,isec),isec=1,NSECTORS)
            end if
          end do
        end do
      end if ! masterProc
     end do EMLAND
 
     deallocate(locemis)
     deallocate(gemis)
     deallocate(lemis)
     if(MasterProc) then
       close(IO_TMP)
       deallocate(globemis)
     end if

  end subroutine EmisOut


end module Emissions_ml
