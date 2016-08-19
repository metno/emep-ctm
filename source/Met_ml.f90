! <Met_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2011 met.no
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

module Met_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
!  Subroutines:      Frequency    Called from:
!    MeteoGridRead                   Unimod
!    MetModel_LandUse                Unimod
!    MeteoRead         3h            Unimod    - puts data into nr
!    metvar            3h            Unimod    -      data into nr
!    metint            20 min        PhyChem, ends with call met_derived
!    met_derived       20 min        metvar, metint - gets u_mid, rho_sruf,
!                                       ustar_nwp, inL_nwp
!    BLPhysics(numt)   3h            metvar, after met_derived
! 
! Alt:
!  Unimod  do numt = 1, ....
!              call MeteoRead - puts data into nr
!              call metvar(numt)
!                    call met_derived      uses q(...,1)
!                    call BLPhysics(numt)  
!              call phyche(numt) ... do i = 1, nstep
!                    chemistry stuff
!                    call metint 
!                          call met_derived
!  metvar - This routines postprocess the meteo fields:
!     ! Unit changes, special definitions etc...
! 
!     e.g. converts wind to u_xmj, v_xmi
! 
!     call BLPhysics(numt)
! 
!  metint
! 
!     !     this routine does the forward linear stepping of the meteorological
!     !     fields read or derived every 3 hours.
!        q(..1) = q(.., 1) + ...
! 
! 
!============================================================================= 

  use BLPhysics_ml,         only : &
     KZ_MINIMUM, KZ_MAXIMUM, KZ_SBL_LIMIT,PIELKE   &
    ,HmixMethod, UnstableKzMethod, StableKzMethod, KzMethod  &
    ,USE_MIN_KZ              & ! From old code, is it needed?
    ,MIN_USTAR_LAND          & ! sets u* > 0.1 m/s over land
    ,OB_invL_LIMIT           & ! 
    ,Test_BLM                & ! Tests all Kz, Hmix routines
    ,PBL_ZiMAX, PBL_ZiMIN    & ! max  and min PBL heights
    ,JericevicRiB_Hmix       & ! TESTING
    ,Venkatram_Hmix          & ! TESTING
    ,Zilitinkevich_Hmix      & ! TESTING
    ,SeibertRiB_Hmix_3d      & ! TESTING
    ,BrostWyngaardKz         & ! TESTING
    ,JericevicKz             & ! TESTING
    ,TI_Hmix                 & ! TESTING or orig
    ,PielkeBlackadarKz       &
    ,O_BrienKz               &
    ,NWP_Kz                  & ! Kz from meteo 
    ,Kz_m2s_toSigmaKz        & 
    ,SigmaKz_2_m2s

  use CheckStop_ml,         only : CheckStop,StopAll
  use Functions_ml,         only : Exner_tab, Exner_nd
  use GridValues_ml,        only : xmd, i_fdom, j_fdom, METEOfelt, projection &
       ,glon,glat, glat_fdom, glon_fdom, MIN_ADVGRIDS   &
       ,Poles, Pole_included, xm_i, xm_j, xm2, sigma_bnd,sigma_mid &
       ,xp, yp, fi, GRIDWIDTH_M,ref_latitude     &
       ,debug_proc, debug_li, debug_lj &
       ,grid_north_pole_latitude,grid_north_pole_longitude &
       ,GlobalPosition,DefGrid,gl_stagg,gb_stagg,A_mid,B_mid

  use Landuse_ml, only : water_cover, water_cover_set
  use MetFields_ml 
  use MicroMet_ml, only : PsiH  ! Only if USE_MIN_KZ
  use ModelConstants_ml,    only : PASCAL, PT, CLOUDTHRES, METSTEP  &
       ,KMAX_BND,KMAX_MID,NMET &
       ,IIFULLDOM, JJFULLDOM, NPROC  &
       ,MasterProc, DEBUG_MET,DEBUG_i, DEBUG_j, identi, V_RAIN, nmax  &
       ,DEBUG_BLM, DEBUG_Kz, DEBUG_SOILWATER & 
       ,NH3_U10   & !FUTURE
       ,DomainName & !HIRHAM,EMEP,EECCA etc.
       ,USE_DUST, USE_SOILWATER & 
       ,nstep,USE_CONVECTION & 
       ,CW_THRESHOLD,RH_THRESHOLD
  use Par_ml           ,    only : MAXLIMAX,MAXLJMAX,GIMAX,GJMAX, me  &
       ,limax,ljmax,li0,li1,lj0,lj1  &
       ,neighbor,WEST,EAST,SOUTH,NORTH,NOPROC  &
       ,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2  &
       ,MFSIZEINP, IRUNBEG,JRUNBEG, tgi0, tgj0,gi0,gj0  &
       ,MSG_INIT3,MSG_READ4, tlimax, tljmax, parinit
  use PhysicalConstants_ml, only : KARMAN, KAPPA, RGAS_KG, CP, GRAV    &
       ,ROWATER, PI
  use TimeDate_ml,          only : current_date, date,Init_nmdays,nmdays, &
       add_secs,timestamp,&
       make_timestamp, make_current_date, nydays, startdate, enddate
  use Io_ml ,               only : ios, IO_ROUGH, &
                                   IO_CLAY, IO_SAND, open_file, IO_LOG
  use ReadField_ml,         only : ReadField ! reads ascii fields
  use NetCDF_ml,         only : printCDF ! testoutputs
  use netcdf
  use TimeDate_ExtraUtil_ml,only: nctime2idate,date2string



  implicit none
  private


  INCLUDE 'mpif.h'
  INTEGER MPISTATUS(MPI_STATUS_SIZE),INFO

!  logical, private, save      ::  debug_procloc = .false.
  integer, private, save      ::  debug_iloc, debug_jloc  ! local coords


  integer, save   :: Nhh &         ! number of field stored per 24 hours
       ,nhour_first&  ! time of the first meteo stored
       ,nrec          ! nrec=record in meteofile, for example
  ! (Nhh=8): 1=00:00 2=03:00 ... 8=21:00
  ! if nhour_first=3 then 1=03:00 2=06:00...8=24:00

  character (len = 100)        ::  field_not_found='field_not_found'


  public :: MeteoGridRead
  public :: MeteoRead
  public :: MetModel_LandUse
  public :: metvar
  public :: metint
  public :: BLPhysics
  public :: GetCDF_short
  public :: extendarea  ! returns array which includes neighbours

contains



  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine MeteoRead(numt)

    !    the subroutine reads meteorological fields and parameters (every
    !    METSTEP-hours) from NetCDF fields-files, divide the fields into
    !       domains    and sends subfields to the processors


    implicit none

    integer, intent(in):: numt

    character (len = 100), save  ::  meteoname   ! name of the meteofile
    character (len = 100)        ::  namefield & ! name of the requested field
         ,unit='   ',validity='    '    ! field is either instaneous or averaged
    integer ::   ndim,nyear,nmonth,nday,nhour
    integer ::   nr   ! Fields are interpolate in
                      ! time (NMET = 2): between nr=1 and nr=2

    type(date)      ::  next_inptime             ! hfTD,addhours_to_input
    type(timestamp) ::  ts_now                   ! time in timestamp format

    real :: nsec                                 ! step in seconds

    real :: temp(MAXLIMAX,MAXLJMAX)!temporary metfields
    ! Avergaing of soil water used box from +/- NEXTEND (e.g. -1 to +1)
    integer, parameter :: NEXTEND = 2 ! no. box to side of (i,j) 
    real, dimension(MAXLIMAX+2*NEXTEND,MAXLJMAX+2*NEXTEND), save  ::&
         xwf  ! extension of water fraction, save after 1st call
    real, dimension(MAXLIMAX+2*NEXTEND,MAXLJMAX+2*NEXTEND)  ::&
         xsw   ! extension of soil water
    real :: tmpsw, landfrac, sumland  ! for soil water averaging
    real :: SoilMax  ! Max value soil-water-deep, used to normalise SW
    integer :: i, j, ii,jj,ii2,jj2
    logical :: fexist  , xwf_done

    nr=2 !set to one only when the first time meteo is read


    if(numt == 1)then !first time meteo is read
       nr = 1
       nrec = 0
       sdot_at_mid = .false.
       foundustar = .false.
       foundsdot = .false.
       foundSST  = .false.
       foundSoilWater  = .false.
       foundSoilWater_deep  = .false.
       foundKz_met = .false.  ! Kz from meteo
       foundu10_met = .false. ! from FUTURE NH3emis
       foundv10_met = .false. ! from FUTURE NH3emis
       foundprecip = .false.
       foundcloudwater = .false.

       next_inptime = current_date

       ! If origin of meteodomain does not coincide with origin of large domain,
       ! xp and yp should be shifted here, and coordinates must be shifted when
       ! meteofields are read (not yet implemented)



    else

       nsec=METSTEP*3600.0 !from hr to sec
       ts_now = make_timestamp(current_date)
       call add_secs(ts_now,nsec)
       next_inptime=make_current_date(ts_now)


    endif




    nyear=next_inptime%year
    nmonth=next_inptime%month

    nday=next_inptime%day
    nhour=next_inptime%hour

    if(  current_date%month == 1 .and.         &
         current_date%day   == 1 .and.         &
         current_date%hour  == 0 )         &
         call Init_nmdays( current_date )



    if(MasterProc .and. DEBUG_MET) write(6,*) &
         '*** nyear,nmonth,nday,nhour,numt,nmdays2'    &
         ,next_inptime%year,next_inptime%month,next_inptime%day    &
         ,next_inptime%hour,numt,nmdays(2)


    !Read rec=1 both for h=0 and h=3:00 in case 00:00 from 1st January is missing
    if((numt-1)*METSTEP<=nhour_first)nrec=0
    nrec=nrec+1





    if(nrec>Nhh.or.nrec==1) then              ! define a new meteo input file
56     FORMAT(a5,i4.4,i2.2,i2.2,a3)
       write(meteoname,56)'meteo',nyear,nmonth,nday,'.nc'
       nrec = 1
       if(nday==1.and.nmonth==1)then
          !hour 00:00 from 1st January may be missing;checking first:
          inquire(file=meteoname,exist=fexist)
          if(.not.fexist)then
             if(MasterProc)write(*,*)trim(meteoname),&
                  ' does not exist; using data from 31 December'
             write(meteoname,56)'meteo',nyear-1,12,31,'.nc'
             nrec=Nhh
          endif
       endif
       if(MasterProc)write(*,*)'reading ',trim(meteoname)
       !could open and close file here instead of in Getmeteofield
    endif


    if(MasterProc .and. DEBUG_MET) write(*,*)'nrec,nhour=',nrec,nhour



    !==============    3D fields (surface) (i,j,k) ==========================
    ndim=3

  !note that u_xmj and v_xmi have dimensions 0:MAXLIJMAX instead of 1:MAXLIJMAX
  !u_xmj(i=0) and v_xmi(j=0) are set in metvar

    namefield='u_wind'
    call Getmeteofield(meteoname,namefield,nrec,ndim,     &
         unit,validity,u_xmj(1:MAXLIMAX,1:MAXLJMAX,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))


    namefield='v_wind'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity,v_xmi(1:MAXLIMAX,1:MAXLJMAX,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    namefield='specific_humidity'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, q(:,:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    namefield='sigma_dot'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, sdot(:,:,:,nr))
    if(validity==field_not_found)then
       foundsdot = .false.
       if(MasterProc.and.numt==1)write(*,*)'WARNING: sigma_dot will be derived from horizontal winds '
    else
       foundsdot = .true.
    endif

    namefield='potential_temperature'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, th(:,:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    namefield='3D_cloudcover'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, cc3d(:,:,:))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    if(trim(validity)/='averaged')then
       if(MasterProc.and.numt==1)write(*,*)'WARNING: 3D cloud cover is not averaged'
    endif

    namefield='precipitation'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, pr(:,:,:))
    foundprecip = .true.
    if(validity==field_not_found)then
       foundprecip = .false.       
       !Will construct 3D precipitations from 2D precipitations
       namefield='large_scale_precipitations'
       call Getmeteofield(meteoname,namefield,nrec,2,&
            unit,validity, surface_precip(:,:))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
       namefield='convective_precipitations'
       call Getmeteofield(meteoname,namefield,nrec,2,&
            unit,validity, temp(:,:))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
       surface_precip=surface_precip+temp

       !if available, will use cloudwater to determine the height of release
       namefield='cloudwater'
       foundcloudwater = .true.
       if(nr==2)cw(:,:,:,1)=cw(:,:,:,2)
       call Getmeteofield(meteoname,namefield,nrec,ndim,&
            unit,validity, cw(:,:,:,nr))
       if(validity==field_not_found)foundcloudwater = .false.
       if(MasterProc.and.foundcloudwater)then
          if(numt==1)write(*,*)' WARNING: 3D precipitations not found. Using 2D precipitations and cloudwater to make 3D'
          if(numt==1)write(unit=IO_LOG,fmt="(a)")"3D precipitations:  derived from 2D and cloudwater"
       endif
       !if cloudwater not available, will use RH to determine the height of release
       if(MasterProc.and..not.foundcloudwater)then
          if(numt==1)write(*,*)' WARNING: 3D precipitations not found. Using 2D precipitations and relative humidity to make 3D'
          if(numt==1)write(unit=IO_LOG,fmt="(a)")"3D precipitations:  derived from 2D and humidity"
       endif

    else
       pr=max(0.0,pr)  ! positive precipitation
    endif


    if(USE_CONVECTION)then
       namefield='convective_updraft_flux'
       call Getmeteofield(meteoname,namefield,nrec,ndim,&
            unit,validity, cnvuf(:,:,:))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
       cnvuf=max(0.0,cnvuf)!no negative upward fluxes
       cnvuf(:,:,KMAX_BND)=0.0!no flux through surface
       cnvuf(:,:,1)=0.0!no flux through top
       
       namefield='convective_downdraft_flux'
       call Getmeteofield(meteoname,namefield,nrec,ndim,&
            unit,validity, cnvdf(:,:,:))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
       cnvdf=min(0.0,cnvdf)!no positive downward fluxes
       cnvdf(:,:,KMAX_BND)=0.0!no flux through surface
       cnvdf(:,:,1)=0.0!no flux through top
    endif

! Kz from meteo
    if (NWP_Kz) then
        namefield='eddy_diffusion_coefficient'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             unit,validity, Kz_met(:,:,:,nr))
        if(validity==field_not_found)then
           foundKz_met = .false.
        if(MasterProc.and.numt==1)write(*,*)' WARNING: Kz will be derived in model '
        else
           foundKz_met = .true.
        endif
        Kz_met=max(0.0,Kz_met)  ! only positive Kz
    end if
    if( debug_proc .and. DEBUG_Kz)then
          write(6,*)               &
         '*** After Kz', sum(Kz_met(:,:,:,nr)), minval(Kz_met(:,:,:,nr)), &
               maxval(Kz_met(:,:,:,nr)),maxval(Kz_met(:,:,KMAX_BND,nr)), &
               DEBUG_Kz, NWP_Kz, nr, nrec, ndim, namefield

    endif




    !==============    2D fields (surface) (i,j)   ============================

    ndim=2

    namefield='surface_pressure'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, ps(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))

    namefield='temperature_2m'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, t2_nwp(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))


    namefield='relative_humidity_2m'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, rh2m(:,:,nr))
    if(validity==field_not_found)then
        if(MasterProc.and.numt==1)write(*,*)'WARNING: relative_humidity_2m not found'
        rh2m(:,:,nr) = -999.9  ! ?
    else
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
       rh2m(:,:,nr) = 0.01 * rh2m(:,:,nr)  ! Convert from % to fraction 
    endif

    namefield='surface_flux_sensible_heat'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, fh(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
    if(validity=='averaged')fh(:,:,1)=fh(:,:,nr)

    namefield='surface_flux_latent_heat'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, fl(:,:,nr))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
    if(validity=='averaged')fl(:,:,1)=fl(:,:,nr)

    namefield='surface_stress'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
         unit,validity, tau(:,:,nr))
    if(validity==field_not_found)then
       namefield='ustar_nwp'
       call Getmeteofield(meteoname,namefield,nrec,ndim,&
            unit,validity, ustar_nwp(:,:))
       call CheckStop(validity==field_not_found, "meteo field not found:" // trim(namefield))
       foundustar=.true.
    else
       tau=max(0.0,tau)
       if(validity=='averaged')tau(:,:,1)=tau(:,:,nr)
    endif
    
    namefield='sea_surface_temperature'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        unit,validity, sst(:,:,nr))
    if(validity==field_not_found)then
       if(MasterProc.and.numt==1)write(*,*)' WARNING: sea_surface_temperature not found '
       foundSST = .false.
    else
       foundSST = .true.
    endif

    namefield='soil_water_content'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        unit,validity, SoilWater(:,:,nr))
    if(validity==field_not_found)then
       if(MasterProc.and.numt==1)write(*,*)' WARNING: SoilWater not found '
       foundSoilWater = .false.

    !  The shallow soil water is intended for the experimental dust modelling
    !  and only implementted with HIRLAM-type inputs so far
    else if ( trim(unit) /= "m" ) then  ! PARLAM/HIRLAM has metres of water in top 7.2 cm
       if(MasterProc.and.numt==1)write(*,*)'WARNING: SoilWater-shallow not HIRLAM, skipping'
       foundSoilWater = .false.
    else
       foundSoilWater = .true.
    endif

  if ( USE_SOILWATER ) then  !just deep here
  !========================================
    namefield='deep_soil_water_content'
    if(DomainName == "HIRHAM" ) then
         if(MasterProc.and.numt==1)write(*,*) " Rename soil water in HIRHAM"
         namefield='soil_water_second_layer'
    end if
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        unit,validity, SoilWater_deep(:,:,nr))
    if(validity==field_not_found)then
       if(MasterProc.and.numt==1)write(*,*)' WARNING: ',trim(namefield),' not found '
       foundSoilWater_deep = .false.
    else
       !<<<<<<< process SW <<<<<<<<<<<<<<<<<<<<<<<
       foundSoilWater_deep = .true.
       if ( trim(unit) == "m" ) then  ! PARLAM has metres of water
          SoilWaterSource = "PARLAM"
          SoilMax = 0.02   
       else if(trim(unit)=='m3/m3')then
         !IFS has a fairly complex soil water system, with field capacity of 
         ! up to 0.766 for organic soils. More medium soils have ca. 0.43
         ! Typical values in January are even down to 0.2. We choose 0.5
         ! as our max, and SoilWater_ml will consider drying effects after
         ! down to 0.25
          call StopAll("Soil Water Handling with IFS not yet implemented")
          SoilWaterSource = "IFS"
          SoilMax = 0.5
       else   ! units not defined yet
          if(numt==1)write(*,*)trim(unit)
          call StopAll("Need units for deep soil water")
       endif

     ! Make SoilWater relative 0 < SW < 1:

       SoilWater_deep(:,:,nr) = min( SoilWater_deep(:,:,nr)/SoilMax, 1.0 ) 

       if(MasterProc.and.numt==1) write(*,*)' Met_ml Soilwater: ' // &
             trim(SoilWaterSource), SoilMax

       if ( water_cover_set ) then  ! smooth the SoilWater values:

         ! If NWP thinks this is a sea-square, but we anyway have land,
         ! the soil moisture might be very strange.  We search neighbouring
         ! grids and make a land-weighted mean SW
         ! Skip on 1st numt, since water fraction set a little later. No harm done...

           call extendarea( SoilWater_deep(:,:,nr), xsw, DEBUG_SOILWATER)

           if ( .not. xwf_done ) then ! only need to do this once
            call extendarea( water_cover(:,:),       xwf, DEBUG_SOILWATER)
            xwf_done = .true.
           end if

           if(DEBUG_SOILWATER .and. debug_proc)&
             write(*,*)'Met_ml water xwf_done: ', me, xwf_done
   
            do j = 1, ljmax
             do i = 1, limax
   
              ! Take a 5x5 average of the land-weighted values for SW. Seems
              !  best not to "believe" NWP models too much for this param, and
              !  the variation in a grid is so big anyway. We aim at the broad
              !  effect. (Alternative might be to find max values?)

                tmpsw = 0.0  ! Relative SW
                sumland  = 0.0
                if( water_cover(i,j) < 0.999 ) then !some land
                   do jj = -NEXTEND, NEXTEND
                     do ii = -NEXTEND, NEXTEND
                         ii2=i+ii+NEXTEND  ! coord in extended array
                         jj2=j+jj+NEXTEND
                         if( xsw(ii2,jj2) > 1.0e-10 ) then ! have some SW to work with
                            landfrac    =  1.0 - xwf(ii2,jj2) 
                            sumland = sumland +  landfrac
                            tmpsw   = tmpsw + landfrac * xsw(ii2,jj2)
                            if ( DEBUG_SOILWATER .and.i==debug_li.and.j==debug_lj ) then
                              write(*,"(a,2i4,8f10.4)") "METSWX: ", ii2, jj2,&
                                water_cover(i,j), xwf(ii2,jj2), &
                                SoilWater_deep(i,j,nr),&
                                xsw(ii2,jj2), tmpsw, landfrac, sumland
                            end if ! DEBUG
                         end if ! xsw
                     end do!ii
                   end do!jj
                   if( sumland > 0.01) then
                        SoilWater_deep(i,j,nr) = tmpsw/sumland
                   else
                        SoilWater_deep(i,j,nr) = 1.0 ! same as sea
                   end if
                else
                   SoilWater_deep(i,j,nr) = 1.0 ! same as sea
                end if ! water_cover
       
                !if( sumland > 0.1 ) then
                !   SoilWater_deep(i,j,nr) = tmpsw/sumland
                !   if(DEBUG_SOILWATER) call CheckStop( tmpsw > sumland, "METSW ERROR")
                !else
                !   SoilWater_deep(i,j,nr) = 1.0  ! also over water
                !end if
   
             end do ! i
            end do ! j
            if ( DEBUG_SOILWATER.and.debug_proc ) then
               i =  debug_li
               j =  debug_lj
               write(*,"(a,f7.4,2i4,f12.4)") "DEBUG_METSWF: ", &
                water_cover(i,j), nr, current_date%day, SoilWater_deep(i,j,nr)
            end if

       else ! for sea values we usually have zero or negative. Set to 1.0
            where ( SoilWater_deep(:,:,nr) < 1.0e-3 ) SoilWater_deep(:,:,nr) = 1.0
       endif ! water_cover_set test
       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    endif ! validity test
  end if ! USE_SOILWATER

  if ( DEBUG_SOILWATER.and.debug_proc ) then
       i =  debug_li
       j =  debug_lj
       write(*,"(a,2i4,f12.4)") "DEBUG_METSWF2: ", &
        nr, current_date%day, SoilWater_deep(i,j,nr)
   end if
  !========================================

    namefield='snow_depth'
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        unit,validity, sdepth(:,:,nr))
    if(validity==field_not_found)then
       if(MasterProc.and.numt==1)write(*,*)' WARNING: snow_depth not found '
       foundsdepth = .false.
    else
       foundsdepth = .true.
    endif


    namefield='fraction_of_ice' !is really percentage
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        unit,validity, ice_nwp(:,:,nr))
    if(validity==field_not_found)then
       if(MasterProc.and.numt==1)write(*,*)' WARNING: ice_nwp coverage (%) not found '
       foundice = .false.
    else
       foundice = .true.
    endif

    namefield='u10'!first component of ws_10m
    call Getmeteofield(meteoname,namefield,nrec,ndim,&
        unit,validity, temp(:,:))
    if(validity==field_not_found)then
       foundws10_met = .false.
    else
       namefield='v10' !second component of ws_10m
       call Getmeteofield(meteoname,namefield,nrec,ndim,&
            unit,validity, ws_10m(:,:,nr))
       if(validity==field_not_found)then
          foundws10_met = .false.
       else
          foundws10_met = .true.
          ws_10m(:,:,nr)=sqrt(ws_10m(:,:,nr)**2+temp(:,:)**2)
!          call printCDF('ws_10m',ws_10m(:,:,1),unit)
       endif
    endif


  end subroutine Meteoread

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine MeteoGridRead(cyclicgrid)

    !   the subroutine reads the grid parameters (projection, resolution etc.)
    !   defined by the meteorological fields
    !

    implicit none

    integer,  intent(out)      :: cyclicgrid
    integer                    :: nyear,nmonth,nday,nhour,k

    character (len = 100),save :: meteoname !name of the meteofile


    nyear=startdate(1)
    nmonth=startdate(2)
    nday=startdate(3)
    nhour=0
    current_date = date(nyear, nmonth, nday, nhour, 0 )
    call Init_nmdays( current_date )

    !*********initialize grid parameters*********
56  FORMAT(a5,i4.4,i2.2,i2.2,a3)
    write(meteoname,56)'meteo',nyear,nmonth,nday,'.nc'
    if(DEBUG_MET.and.MasterProc)write(*,*)'looking for ',trim(meteoname)


    call Getgridparams(meteoname,GRIDWIDTH_M,xp,yp,fi,xm_i,xm_j,xm2,&
         ref_latitude,sigma_mid,Nhh,nyear,nmonth,nday,nhour,nhour_first&
         ,cyclicgrid)


    if(MasterProc .and. DEBUG_MET)then
       write(*,*)'sigma_mid:',(sigma_mid(k),k=1,20)
       write(*,*)'grid resolution:',GRIDWIDTH_M
       write(*,*)'xcoordinate of North Pole, xp:',xp
       write(*,*)'ycoordinate of North Pole, yp:',yp
       write(*,*)'longitude rotation of grid, fi:',fi
       write(*,*)'true distances latitude, ref_latitude:',ref_latitude
    endif

    call DefGrid()!defines: i_fdom,j_fdom,i_local, j_local,xmd,xm2ji,xmdji,
                  !         sigma_bnd,carea,gbacmax,gbacmin,glacmax,glacmin

  end subroutine Meteogridread


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine metvar(numt)

    ! This routines postprocess the meteo fields:
    ! Unit changes, special definitions etc...


    implicit none

    integer, intent(in):: numt

    !    local

    real,   dimension(KMAX_MID)          ::  prhelp, exf2
    real,   dimension(KMAX_BND)          ::  exf1

    real,   dimension(MAXLJMAX,KMAX_MID) :: usnd   ! send in x
    real,   dimension(MAXLIMAX,KMAX_MID) :: vsnd   ! and in y direction
    real,   dimension(MAXLJMAX,KMAX_MID) :: urcv   ! rcv in x
    real,   dimension(MAXLIMAX,KMAX_MID) :: vrcv   ! and in y direction

    real   divt, p1, p2
    real   prhelp_sum,divk(KMAX_MID),sumdiv
    real   inv_METSTEP

    integer :: i, j, k, kk, nr,info
    integer request_s,request_n,request_e,request_w
    real ::Ps_extended(0:MAXLIMAX+1,0:MAXLJMAX+1),Pmid,Pu1,Pu2,Pv1,Pv2
    real :: relh1,relh2,temperature,swp,wp

    nr = 2
    if (numt == 1) then

       nr = 1

       !-------------------------------------------------------------------
       !  Initialisations:

       call Exner_tab()

       ! Look for processor containing debug coordinates
       debug_iloc    = -999
       debug_jloc    = -999

       do i = 1, limax
          do j = 1, ljmax
!             if (DEBUG_MET .and. &
              if( i_fdom(i) == DEBUG_I .and. j_fdom(j) == DEBUG_J ) then
                debug_proc = .true.
                debug_iloc    = i
                debug_jloc    = j
             end if
          end do
       end do
       if( debug_proc ) write(*,*) "DEBUG EXNER me", me, Exner_nd(99500.0)
       !-------------------------------------------------------------------

    end if !numt == 1


    divt = 1./(3600.0*METSTEP)


    if (neighbor(EAST) .ne. NOPROC) then
       do k = 1,KMAX_MID
          do j = 1,ljmax
             usnd(j,k) = u_xmj(limax,j,k,nr)
          enddo
       enddo
       CALL MPI_ISEND( usnd, 8*MAXLJMAX*KMAX_MID, MPI_BYTE,  &
            neighbor(EAST), MSG_WEST2, MPI_COMM_WORLD, request_e, INFO)
    endif

    if (neighbor(NORTH) .ne. NOPROC) then
       do k = 1,KMAX_MID
          do i = 1,limax
             vsnd(i,k) = v_xmi(i,ljmax,k,nr)
          enddo
       enddo
       CALL MPI_ISEND( vsnd , 8*MAXLIMAX*KMAX_MID, MPI_BYTE,  &
            neighbor(NORTH), MSG_SOUTH2, MPI_COMM_WORLD, request_n, INFO)
    endif


    !     receive from WEST neighbor if any

    if (neighbor(WEST) .ne. NOPROC) then

       CALL MPI_RECV( urcv, 8*MAXLJMAX*KMAX_MID, MPI_BYTE, &
            neighbor(WEST), MSG_WEST2, MPI_COMM_WORLD, MPISTATUS, INFO)
       do k = 1,KMAX_MID
      do j = 1,ljmax
             u_xmj(0,j,k,nr) = urcv(j,k)
          enddo
       enddo

    else

       do k = 1,KMAX_MID
          do j = 1,ljmax
             u_xmj(0,j,k,nr) = u_xmj(1,j,k,nr)
          enddo
       enddo


    endif

    !     receive from SOUTH neighbor if any

    if (neighbor(SOUTH) .ne. NOPROC) then

       CALL MPI_RECV( vrcv, 8*MAXLIMAX*KMAX_MID, MPI_BYTE,  &
            neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_WORLD, MPISTATUS, INFO)
       do k = 1,KMAX_MID
          do i = 1,limax
             v_xmi(i,0,k,nr) = vrcv(i,k)
          enddo
       enddo

    else

       if(Poles(2)/=1)  then
          do k = 1,KMAX_MID
             do i = 1,limax
                v_xmi(i,0,k,nr) = v_xmi(i,1,k,nr)
             enddo
          enddo
       else
          !"close" the South pole
          do k = 1,KMAX_MID
             do i = 1,limax
                v_xmi(i,0,k,nr) = 0.0
             enddo
          enddo
       endif

    endif

    if (neighbor(NORTH) == NOPROC.and.Poles(1)==1) then
       !"close" the North pole
       do k = 1,KMAX_MID
          do i = 1,limax
             v_xmi(i,ljmax,k,nr) = 0.0
          enddo
       enddo
    endif

    if (neighbor(EAST) .ne. NOPROC) then
       CALL MPI_WAIT(request_e, MPISTATUS, INFO)
    endif

    if (neighbor(NORTH) .ne. NOPROC) then
       CALL MPI_WAIT(request_n, MPISTATUS, INFO)
    endif

    if(.not.foundprecip)then
       !have to construct 3D precipitations from 2D precipitations
       if(foundcloudwater)then
          !if available, use cloudwater to determine the height of release
          if(nr==1)cw(:,:,:,2)=cw(:,:,:,nr)       
          do j=1,ljmax
             do i=1,limax
                pr(i,j,KMAX_MID)= surface_precip(i,j)*10800000.0!guarantees precip at surface
                do k=1,KMAX_MID-1
                   if(cw(i,j,k,2)+cw(i,j,k,1)>CW_THRESHOLD)then
                      !fill the column up to this level with constant precip
                      do kk=k,KMAX_MID-1
                         pr(i,j,kk)= surface_precip(i,j)*METSTEP*3600.0*1000.0!from m/s to mm/METSTEP              
                      enddo
                      exit
                   else
                   pr(i,j,k)=0.0               
                endif
             enddo
          enddo
       enddo

       else
          !will use RH to determine the height of release (less accurate than cloudwater)
          do j=1,ljmax
             do i=1,limax 
                pr(i,j,KMAX_MID)= surface_precip(i,j)*10800000.0!guarantees precip at surface
                do k=1,KMAX_MID-1
                   !convert from potential temperature into absolute temperature
                   temperature = th(i,j,k,nr)*exp(KAPPA*log((A_mid(k) + B_mid(k)*ps(i,j,nr)*100)*1.e-5))!Pa, Ps still in hPa here
                   !saturation water pressure
                   swp=611.2*exp(17.67*(temperature-273.15)/(temperature-29.65))
                   !water pressure
                   wp=q(i,j,k,nr)*(A_mid(k) + B_mid(k)*ps(i,j,nr)*100)/0.622
                   relh2=wp/swp
                   !convert from potential temperature into absolute temperature
                   temperature = th(i,j,k,1)*exp(KAPPA*log((A_mid(k) + B_mid(k)*ps(i,j,1)*100)*1.e-5))!Pa, Ps still in hPa here
                   !saturation water pressure
                   swp=611.2*exp(17.67*(temperature-273.15)/(temperature-29.65))
                   !water pressure
                   wp=q(i,j,k,1)*(A_mid(k) + B_mid(k)*ps(i,j,1)*100)/0.622
                   relh1=wp/swp
                   if(relh1>RH_THRESHOLD.or.relh2>RH_THRESHOLD)then
                      !fill the column up to this level with constant precip
                      do kk=k,KMAX_MID-1
                         pr(i,j,kk)= surface_precip(i,j)*10800000.0!3hours and m->mm              
                      enddo
                      exit
                   else
                      pr(i,j,k)=0.0               
                   endif                                     
                enddo
             enddo
          enddo
       endif
       pr=max(0.0,pr)  ! positive precipitation
    endif


    inv_METSTEP = 1.0/METSTEP

    do j = 1,ljmax
       do i = 1,limax

          !     conversion of pressure from hPa to Pascal.
          ps(i,j,nr) = ps(i,j,nr)*PASCAL

          ! surface precipitation, mm/hr
          !NB: surface_precip is different than the one read directly from the 
          !metfile  (which has different units, and is the sum of the 2D 
          !large_scale_precipitations+convective_precipitations)
          surface_precip(i,j) = pr(i,j,KMAX_MID) * inv_METSTEP 

          rho_surf(i,j)  = ps(i,j,nr)/(RGAS_KG * t2_nwp(i,j,nr) )

          !     For MM5 we get u*, not tau. Since it seems better to
          !     interpolate tau than u*  between time-steps we convert
          if ( foundustar) then
             tau(i,j,nr)    = ustar_nwp(i,j)*ustar_nwp(i,j)* rho_surf(i,j)
          end if


          prhelp_sum = 0.0
          prhelp(1) = max(pr(i,j,1),0.)

          prhelp_sum = prhelp_sum + prhelp(1)

          !  pr is 3 hours accumulated precipitation in mm in each
          !  layer summed from above. This is first converted to precipitation
          !  release in each layer.

          do k = 2,KMAX_MID
             prhelp(k) = pr(i,j,k) - pr(i,j,k-1)
          enddo

          ! accumulated deposition over 3 hour interval
          ! k=KMAX_MID now includes accumulated precipitation over all layers
          ! evaporation has been set to zero as it is not accounted for in the
          ! wet deposition
          ! Add up in WetDeposition, to have the prec used in the model

          pr(i,j,:) = prhelp(:)*divt

          !   interpolation of sigma dot for half layers

          if(foundsdot.and.sdot_at_mid)then
             do k = KMAX_MID,2,-1
                sdot(i,j,k,nr) = sdot(i,j,k-1,nr)            &
                     + (sdot(i,j,k,nr)-sdot(i,j,k-1,nr))   &
                     * (sigma_bnd(k)-sigma_mid(k-1))       &
                     / (sigma_mid(k)-sigma_mid(k-1))
             enddo
          endif

          !    set sdot equal to zero at the top and bottom of atmosphere.
          sdot(i,j,KMAX_BND,nr)=0.0
          sdot(i,j,1,nr)=0.0

          !    conversion from % to fractions (<0,1>) for cloud cover
          !    calculation of cc3dmax (see eulmc.inc) -
          !    maximum of cloud fractions for layers above a given layer
          cc3d(i,j,1) = 0.01 * cc3d(i,j,1)
          cc3dmax(i,j,1) = cc3d(i,j,1)

          lwc(i,j,:)=0.
          do k=2,KMAX_MID
             cc3d(i,j,k) = 0.01 * cc3d(i,j,k)
             cc3dmax(i,j,k) = amax1(cc3dmax(i,j,k-1),cc3d(i,j,k-1))
             lwc(i,j,k)=0.6e-6*cc3d(i,j,k)
          enddo

       enddo
    enddo


    !   derive the meteorological parameters from the basic parameters
    !   read from field-files.

    do j = 1,ljmax
       do i = 1,limax
          p1 = sigma_bnd(KMAX_BND)*(ps(i,j,nr) - PT) + PT

          exf1(KMAX_BND) = CP * Exner_nd(p1)

          z_bnd(i,j,KMAX_BND) = 0.0

          do k = KMAX_MID,1,-1

             !   eddy diffusivity in the surface-layer follows the 
             !   formulation usedin the nwp-model which is based on 
             !   Louis (1979), (see mc7e.f).
             !   exner-function of the half-layers

             p1 = sigma_bnd(k)*(ps(i,j,nr) - PT) + PT
             exf1(k) = CP * Exner_nd( p1 )

             !   exner-function of the full-levels
             p2 = sigma_mid(k)*(ps(i,j,nr) - PT) + PT
             exf2(k) = CP * Exner_nd(p2)

             !     height of the half-layers
             z_bnd(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
                  (exf1(k+1) - exf1(k)))/GRAV


             !     height of the full levels.

             z_mid(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
                  (exf1(k+1) - exf2(k)))/GRAV

             roa(i,j,k,nr) = CP*((ps(i,j,nr) - PT)*sigma_mid(k) + PT)/      &
                  (RGAS_KG*th(i,j,k,nr)*exf2(k))

          enddo  ! k

       enddo
    enddo
    !-----------------------------------------------------------------------

    if( DEBUG_MET .and. debug_proc ) then
     ! Note: example results are given in MetFields
       write(*,*) "DEBUG meIJ" , me, limax, ljmax
       do k = 1, KMAX_MID
          write(6,"(a12,2i3,9f12.4)") "DEBUG_Z",me, k, &
            z_bnd(debug_iloc,debug_jloc,k), z_mid(debug_iloc,debug_jloc,k)!,&
!            sigma_mid(k) !,            zm3d(debug_iloc,debug_jloc,k)
       end do
    end if


    !     Horizontal velocity divided by map-factor.
    !divide by the scaling in the perpendicular direction to get effective 
    !u_xmj and v_xmi
    !(for conformal projections like Polar Stereo, xm_i and xm_j are equal)
    do k = 1,KMAX_MID
       do j = 1,ljmax
          do i = 0,limax
             u_xmj(i,j,k,nr) = u_xmj(i,j,k,nr)/xm_j(i,j)
          enddo
       enddo
       do j = 0,ljmax
          do i = 1,limax
             v_xmi(i,j,k,nr) = v_xmi(i,j,k,nr)/xm_i(i,j)
          enddo
       enddo
    enddo



    if(.not.foundsdot)then
       ! sdot derived from divergence=0 principle
       do j = 1,ljmax
          do i = 1,limax
             Ps_extended(i,j) = Ps(i,j,nr)
          enddo
       enddo
!Get Ps at edges from neighbors
!we reuse usnd, vsnd etc
    if (neighbor(EAST) .ne. NOPROC) then
          do j = 1,ljmax
             usnd(j,1) = ps(limax,j,nr)
          enddo
       CALL MPI_ISEND( usnd, 8*MAXLJMAX, MPI_BYTE,  &
            neighbor(EAST), MSG_WEST2, MPI_COMM_WORLD, request_e, INFO)
    endif
    if (neighbor(NORTH) .ne. NOPROC) then
          do i = 1,limax
             vsnd(i,1) = ps(i,ljmax,nr)
          enddo
       CALL MPI_ISEND( vsnd , 8*MAXLIMAX, MPI_BYTE,  &
            neighbor(NORTH), MSG_SOUTH2, MPI_COMM_WORLD, request_n, INFO)
    endif

    !     receive from WEST neighbor if any
    if (neighbor(WEST) .ne. NOPROC) then
       CALL MPI_RECV( urcv, 8*MAXLJMAX, MPI_BYTE, &
            neighbor(WEST), MSG_WEST2, MPI_COMM_WORLD, MPISTATUS, INFO)
      do j = 1,ljmax
             Ps_extended(0,j) = urcv(j,1)
          enddo
    else
          do j = 1,ljmax
             Ps_extended(0,j) = Ps_extended(1,j)
          enddo
    endif
    !     receive from SOUTH neighbor if any
    if (neighbor(SOUTH) .ne. NOPROC) then
       CALL MPI_RECV( vrcv, 8*MAXLIMAX, MPI_BYTE,  &
            neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_WORLD, MPISTATUS, INFO)
          do i = 1,limax
             Ps_extended(i,0) = vrcv(i,1)
          enddo
    else
          do i = 1,limax
             Ps_extended(i,0) = Ps_extended(i,1)
          enddo
    endif
    if (neighbor(WEST) .ne. NOPROC) then
          do j = 1,ljmax
             usnd(j,2) = ps(1,j,nr)
          enddo
       CALL MPI_ISEND( usnd(1,2), 8*MAXLJMAX, MPI_BYTE,  &
            neighbor(WEST), MSG_WEST2, MPI_COMM_WORLD, request_w, INFO)
    endif
    if (neighbor(SOUTH) .ne. NOPROC) then
          do i = 1,limax
             vsnd(i,2) = ps(i,1,nr)
          enddo
       CALL MPI_ISEND( vsnd(1,2) , 8*MAXLIMAX, MPI_BYTE,  &
            neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_WORLD, request_s, INFO)
    endif



    !     receive from EAST neighbor if any
    if (neighbor(EAST) .ne. NOPROC) then
       CALL MPI_RECV( urcv, 8*MAXLJMAX, MPI_BYTE, &
            neighbor(EAST), MSG_WEST2, MPI_COMM_WORLD, MPISTATUS, INFO)
      do j = 1,ljmax
             Ps_extended(limax+1,j) = urcv(j,1)
          enddo
    else
          do j = 1,ljmax
             Ps_extended(limax+1,j) = Ps_extended(limax,j)
          enddo
    endif
    !     receive from NORTH neighbor if any
    if (neighbor(NORTH) .ne. NOPROC) then
       CALL MPI_RECV( vrcv, 8*MAXLIMAX, MPI_BYTE,  &
            neighbor(NORTH), MSG_SOUTH2, MPI_COMM_WORLD, MPISTATUS, INFO)
          do i = 1,limax
             Ps_extended(i,ljmax+1) = vrcv(i,1)
          enddo
    else
          do i = 1,limax
             Ps_extended(i,ljmax+1) = Ps_extended(i,ljmax)
          enddo
    endif

    if (neighbor(EAST) .ne. NOPROC) then
       CALL MPI_WAIT(request_e, MPISTATUS, INFO)
    endif

    if (neighbor(NORTH) .ne. NOPROC) then
       CALL MPI_WAIT(request_n, MPISTATUS, INFO)
    endif
    if (neighbor(WEST) .ne. NOPROC) then
       CALL MPI_WAIT(request_w, MPISTATUS, INFO)
    endif

    if (neighbor(SOUTH) .ne. NOPROC) then
       CALL MPI_WAIT(request_s, MPISTATUS, INFO)
    endif

!note that u_xmj and v_xmi have already been divided by xm here
       do j = 1,ljmax
          do i = 1,limax
             Pmid=Ps_extended(i,j)-PT
             Pu1=0.5*(Ps_extended(i-1,j)+Ps_extended(i,j))-PT
             Pu2=0.5*(Ps_extended(i+1,j)+Ps_extended(i,j))-PT
             Pv1=0.5*(Ps_extended(i,j-1)+Ps_extended(i,j))-PT
             Pv2=0.5*(Ps_extended(i,j+1)+Ps_extended(i,j))-PT

             sdot(i,j,KMAX_BND,nr)=0.0
             sdot(i,j,1,nr)=0.0
             sumdiv=0.0
             do k=1,KMAX_MID
                divk(k)=((u_xmj(i,j,k,nr)*Pu2-u_xmj(i-1,j,k,nr)*Pu1)         &
                     + (v_xmi(i,j,k,nr)*Pv2-v_xmi(i,j-1,k,nr)*Pv1))          &
                     * xm2(i,j)*(sigma_bnd(k+1)-sigma_bnd(k))  &
                     / GRIDWIDTH_M/Pmid
                sumdiv=sumdiv+divk(k)
             enddo
             sdot(i,j,KMAX_MID,nr)=-(sigma_bnd(KMAX_MID+1)-sigma_bnd(KMAX_MID))&
                                    *sumdiv+divk(KMAX_MID)
             do k=KMAX_MID-1,2,-1
                sdot(i,j,k,nr)=sdot(i,j,k+1,nr)-(sigma_bnd(k+1)-sigma_bnd(k))&
                                                *sumdiv+divk(k)
             enddo
          enddo
       enddo
    endif



    call met_derived(nr) !compute derived meteo fields

    call BLPhysics(numt)

  end subroutine metvar

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



  subroutine metint

    !     this routine does the forward linear stepping of the meteorological
    !     fields read or derived every 3 hours.


    implicit none

    real :: div

    if (nstep.lt.nmax) then

       div = 1./real(nmax-(nstep-1))

       u_xmj(:,:,:,1)    = u_xmj(:,:,:,1)                 &
            + (u_xmj(:,:,:,2) - u_xmj(:,:,:,1))*div
       v_xmi(:,:,:,1)    = v_xmi(:,:,:,1)                 &
            + (v_xmi(:,:,:,2) - v_xmi(:,:,:,1))*div
       sdot(:,:,:,1) = sdot(:,:,:,1)             &
            + (sdot(:,:,:,2) - sdot(:,:,:,1))*div
       th(:,:,:,1)   = th(:,:,:,1)                 &
            + (th(:,:,:,2) - th(:,:,:,1))*div
       q(:,:,:,1)    = q(:,:,:,1)                 &
            + (q(:,:,:,2) - q(:,:,:,1))*div
       SigmaKz(:,:,:,1)  = SigmaKz(:,:,:,1)                 &
            + (SigmaKz(:,:,:,2) - SigmaKz(:,:,:,1))*div
       roa(:,:,:,1)  = roa(:,:,:,1)                 &
            + (roa(:,:,:,2) - roa(:,:,:,1))*div
       ps(:,:,1)     = ps(:,:,1)                 &
            + (ps(:,:,2) - ps(:,:,1))*div
       t2_nwp(:,:,1) = t2_nwp(:,:,1)                 &
            + (t2_nwp(:,:,2) - t2_nwp(:,:,1))*div
       rh2m(:,:,1) = rh2m(:,:,1)  &
            + (rh2m(:,:,2) - rh2m(:,:,1))*div
       SoilWater(:,:,1) = SoilWater(:,:,1)   &
            + (SoilWater(:,:,2) - SoilWater(:,:,1))*div
       SoilWater_deep(:,:,1) = SoilWater_deep(:,:,1)    &
            + (SoilWater_deep(:,:,2) - SoilWater_deep(:,:,1))*div


       fh(:,:,1)     = fh(:,:,1)                 &
            + (fh(:,:,2) - fh(:,:,1))*div
       fl(:,:,1)     = fl(:,:,1)                 &
            + (fl(:,:,2) - fl(:,:,1))*div
       tau(:,:,1)    = tau(:,:,1)                 &
            + (tau(:,:,2) - tau(:,:,1))*div
       sst(:,:,1)    = sst(:,:,1)                 &
            + (sst(:,:,2)   - sst(:,:,1))*div
       sdepth(:,:,1)    = sdepth(:,:,1)                 &
            + (sdepth(:,:,2)   - sdepth(:,:,1))*div
       ice_nwp(:,:,1)    = ice_nwp(:,:,1)                 &
            + (ice_nwp(:,:,2)   - ice_nwp(:,:,1))*div
       if(foundws10_met) ws_10m(:,:,1) = ws_10m(:,:,1) &
            + (ws_10m(:,:,2) - ws_10m(:,:,1))*div

       !  precipitation and cloud cover are no longer interpolated

    else

       !     assign the the meteorological data at time-level 2 to level 1 for
       !     the next 6 hours integration period before leaving the inner loop.

       u_xmj(:,:,:,1)    = u_xmj(:,:,:,2)
       v_xmi(:,:,:,1)    = v_xmi(:,:,:,2)
       sdot(:,:,:,1) = sdot(:,:,:,2)
       th(:,:,:,1)   = th(:,:,:,2)
       q(:,:,:,1)    = q(:,:,:,2)
       SigmaKz(:,:,:,1)  = SigmaKz(:,:,:,2)
       roa(:,:,:,1)  = roa(:,:,:,2)
       !  - note we need pressure first before surface_pressure
       ps(:,:,1)     = ps(:,:,2)
       t2_nwp(:,:,1) = t2_nwp(:,:,2)
       rh2m(:,:,1) = rh2m(:,:,2)
       SoilWater(:,:,1) = SoilWater(:,:,2)
       SoilWater_deep(:,:,1) = SoilWater_deep(:,:,2)
       sdepth(:,:,1) = sdepth(:,:,2)
       ice_nwp(:,:,1) = ice_nwp(:,:,2)

       fh(:,:,1)     = fh(:,:,2)
       tau(:,:,1)    = tau(:,:,2)
       fl(:,:,1)     = fl(:,:,2)

       sst(:,:,1)    = sst(:,:,2)

       if(foundws10_met) ws_10m(:,:,1) = ws_10m(:,:,2)

    endif

    call met_derived(1) !update derived meteo fields

  end subroutine metint

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine met_derived(nt)

    ! This routine calculates fields derived from meteofields.
    ! The interpolation in time is done for the meteofields and the
    ! fields here are derived from the interpolated fields after
    ! each interpolation (i.e. every dt_advec).
    ! CPU costly fields (those with special functions like log )
    ! can be computed in metvar only once every METSTEP and interpolated
    ! in metint.

    !horizontal wind speed (averaged over the four edges)
    !Note that u_xmj and v_xmi are wind velocities divided by xm
    !At present u_ref is defined at KMAX_MID


    implicit none
    integer, intent(in) :: nt  ! set to 1 from metint or nr from matvar
    integer ::i,j, k
    logical :: DEBUG_DERIV = .false.

    do k = 1, KMAX_MID
    do j = 1,ljmax
       do i = 1,limax
           u_mid(i,j,k) = 0.5*( u_xmj(i,j,k,nt)*xm_j(i,j) + &
                                u_xmj(i-1,j,k,nt)*xm_j(i-1,j) )
           v_mid(i,j,k) = 0.5*( v_xmi(i,j-1,k,nt)*xm_i(i,j-1) + &
                                v_xmi(i,j,k,nt)*xm_i(i,j))
       enddo
    enddo
    enddo !k
    do j = 1,ljmax
       do i = 1,limax
          u_ref(i,j)= sqrt( u_mid(i,j,KMAX_MID)**2 + v_mid(i,j,KMAX_MID)**2 )
       enddo
    enddo


    ! Tmp ustar solution. May need re-consideration for MM5 etc., but
    ! basic principal should be that fm is interpolated with time, and
    ! ustar derived from this.

    !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

    forall( i=1:limax, j=1:ljmax )
       rho_surf(i,j)  = ps(i,j,nt)/(RGAS_KG * t2_nwp(i,j,nt) )
    end forall

    if(.not. foundustar)then
       forall( i=1:limax, j=1:ljmax )
          ustar_nwp(i,j)   = sqrt( tau(i,j,nt)/rho_surf(i,j) )
       end forall
    endif

    !ds 25/2/2009.. following Branko's comments, 
    ! we limit u* to a physically plausible value over land
    ! to prevent numerical problems, and to account for enhanced
    ! mixing which is usually found over real terrain

    where ( nwp_sea ) 
       ustar_nwp = max( ustar_nwp, 1.0e-5 )
    elsewhere 
       ustar_nwp = max( ustar_nwp, MIN_USTAR_LAND )
    end where



    forall( i=1:limax, j=1:ljmax )
     invL_nwp(i,j)  = KARMAN * GRAV * fh(i,j,nt) & ! - disliked by gfortran
            / (CP*rho_surf(i,j) * ustar_nwp(i,j)**3 * t2_nwp(i,j,1) )
    end forall

    where ( invL_nwp < -1.0 ) 
     invL_nwp  = -1.0
    else where ( invL_nwp > 1.0 ) 
     invL_nwp  = 1.0
    end where 
     

    if ( DEBUG_DERIV .and. debug_proc ) then
       i = debug_iloc
       j = debug_jloc
       write(*,*) "MET_DERIV DONE ", me, nt, ustar_nwp(i,j), rho_surf(i,j), &
            fh(i,j,nt), invL_nwp(i,j)
    end if
    !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa


  end subroutine met_derived

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




  subroutine MetModel_LandUse(callnum)

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !     This subroutine reads parameterfields from file
    !     reading surface roughness classes from file: landsea_mask.dat
    !
    !     ... fields as used in meteorological model
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    implicit none

    integer, intent(in) :: callnum
    integer ::  i,j

    real, dimension(MAXLIMAX,MAXLJMAX) :: r_class  ! Roughness (real)

    character*20 fname
    logical :: needed_found

    ios = 0

    if ( callnum == 1  ) then

       if ( MasterProc  ) then
          write(fname,fmt='(''landsea_mask.dat'')')
          write(6,*) 'reading land-sea map from ',fname
       end if
       needed_found=.false.
       call ReadField(IO_ROUGH,fname,r_class,needed_found,fill_needed=.true.)

       ! And convert from real to integer field

       nwp_sea(:,:) = .false.

       if(needed_found)then
          do j=1,ljmax
             do i=1,limax
                if ( nint(r_class(i,j)) == 0 ) nwp_sea(i,j) = .true.
             enddo
          enddo
       endif

!.. Clay soil content    !
      ios = 0

      if ( USE_DUST ) then
        if ( MasterProc  ) then
           write(fname,fmt='(''clay_frac.dat'')') 
           write(6,*) 'filename for clay fraction ',fname, IO_CLAY, ios
        end if

        call ReadField(IO_CLAY,fname,clay_frac)
   
 ! Convert from percent to fraction
      
        do j=1,ljmax
           do i=1,limax
              clay_frac(i,j) = 0.01 * clay_frac(i,j)
           enddo
        enddo
!.. Sand soil content
     ios = 0

        if ( me == 0  ) then
           write(fname,fmt='(''sand_frac.dat'')') 
           write(6,*) 'filename for sand fraction ',fname, IO_SAND, ios
        end if

        call ReadField(IO_SAND,fname,sand_frac)
   
 ! Convert from percent to fraction
      
        do j=1,ljmax
           do i=1,limax
              sand_frac(i,j) = 0.01 * sand_frac(i,j)
           enddo
        enddo
      end if ! USE_DUST
 
    end if ! callnum == 1

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  end subroutine MetModel_LandUse
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine BLPhysics(numt)
    !c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c    First written by Trond Iversen,  modified by Hugo Jakobsen, 060994
    !c    Extensive modifications to structure by Dave Simpson, March 2010. 
    !c    Some code moved to BLPhysics_ml, together with additinal options. 
    !c    Also now includes Haldis's use of NWP Kz values.
    !c      ** not optimised, bujt called only at 3 h intervals
    !c
    !c-----------------------------------------------------------------
    !c
    !!        This routine calculates the exner function,
    !!        the  geopotential height, and the vertical exchange coefficient
    !!        in sigma surfaces.
    !!        The height zi of the "well mixed layer" or ABL-height
    !!        is also calculated.
    !c
    !c
    !c    if nroa = 1 also roa is calculated.
    !c    if nfirst=1 for the initial timelevel
    !c
    !c
    !c-----------------------------------------------------------------
    !c    routines called:
    !c
    !c        Several options for Kz, Hmix
    !c        smoosp
    !c
    !c
    !c-----------------------------------------------------------------
    !c
    !c**********************************************************************
    logical, parameter :: TKE_DIFF = .false.  !!! CODE NEEDS TESTING/TIDY UP

    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID)::exnm
    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND)::exns

    real :: p_m, p_s, hs

    real, dimension(KMAX_BND) :: p_bnd !TESTzi
    real, dimension(KMAX_MID) :: Kz_nwp 
    real    :: Kz_min, stab_h
!    logical :: Pielke_flag    ! choice in Blackadar/Pielke equations

    integer i,j,k,numt, nr

    call CheckStop( KZ_SBL_LIMIT < 1.01*KZ_MINIMUM,   &
         "SBLlimit too low! in Met_ml")

    !     Preliminary definitions

    nr = 2
    if (numt == 1) nr = 1

    Kz_m2s(:,:,:)= 0.
    Kz_nwp(:)    = -99.0   ! store for printout. only set if read from NWP

    !c..................................
    !c..exner-functions (j/k kg)
    !c
    do  k=1,KMAX_MID
    do j=1,ljmax
       do i=1,limax

          p_m = PT + sigma_mid(k)*(ps(i,j,nr) - PT)
          p_s = PT + sigma_bnd(k)*(ps(i,j,nr) - PT)

             exnm(i,j,k)= CP * Exner_nd(p_m) !c..exner (j/k kg)
             exns(i,j,k)= CP * Exner_nd(p_s)
          end do
       end do
    end do


! Are the invL and fh comparable??
    if  ( debug_proc .and. DEBUG_Kz) then
      i = debug_iloc
      j = debug_jloc
      write(*,"(a,i4,2f12.5)") "TESTNR th ", nr , th(i,j,20,1), th(i,j,20,nr)
      write(*,"(a,i4,2f12.5,es10.2)") "TESTNR fh ", nr , fh(i,j,1), fh(i,j,nr), invL_nwp(i,j)
      write(*,"(a,i4,2es10.2)") "TESTNR ps ", nr , ps(i,j,1), ps(i,j,nr)
    end if



   !SSSSSSSSSSSSSSSSSS Start choice of Kz and Hmix methods SSSSSSSSSSSSSSSSSS


    if (NWP_Kz .and. foundKz_met ) then  ! read from met data

       !hb, +ds rewrote to reduce number of lines. LAter we should remove Kz_met
       ! and Kz_m2s

       forall(i=1:limax,j=1:ljmax,k=2:KMAX_MID)
             SigmaKz(i,j,k,nr)=Kz_met(i,j,k,nr)/(60*60*3)
       end forall

       call SigmaKz_2_m2s( SigmaKz(:,:,:,nr), roa(:,:,:,nr),ps(:,:,nr), Kz_m2s )

       if( debug_proc ) Kz_nwp(:) = Kz_m2s(debug_iloc,debug_jloc,:) !for printout


       if( debug_proc .and. DEBUG_Kz)then            
          write(6,*) '*** After Set SigmaKz', sum(SigmaKz(:,:,:,nr)), &
             minval(SigmaKz(:,:,:,nr)), maxval(SigmaKz(:,:,:,nr)), &
             DEBUG_Kz, 'NWP_Kz:',NWP_Kz, &
            '*** After convert to z',sum(Kz_m2s(:,:,:)), &
            minval(Kz_m2s(:,:,:)), maxval(Kz_m2s(:,:,:))
          write(6,*) 'DS  After Set SigmaKz KTOP', Kz_met(debug_iloc,debug_jloc,1,nr)
       endif

    else   ! Not NWP Kz. Must calculate

         ! 1/  Get Kz first from PielkeBlackadar methods
         ! Use for all methods except NWP_Kz
         ! Do the physics for each i,j for now. Optimise later

          do j=1,ljmax
             do i=1,limax

            call PielkeBlackadarKz ( &
              u_mid(i,j,:),  v_mid(i,j,:),  &
              z_mid(i,j,:),  z_bnd(i,j,:),  &
              th(i,j,:,nr),  Kz_m2s(i,j,:), &
              PIELKE, &     !Pielke_flag, &
              .false. )
              !( debug_proc .and. i == debug_iloc .and. j == debug_jloc )  )

             enddo
          enddo

        !======================================================================
        ! Hmix choices:

          if ( HmixMethod == "TIZi" ) then
    
              ! 2/ Get Mixing height from "orig" method
              !   "old" exner-function of the full-levels

            do j=1,ljmax
               do i=1,limax

                  p_bnd(:) = sigma_bnd(:)*(ps(i,j,nr) - PT) + PT
                 ! p_mid(:) = sigma_mid(:)*(ps(i,j,nr) - PT) + PT
                 !exf2(:) = CP * Exner_nd(p_mid(:))
    
               call TI_Hmix ( &     ! Original EMEP method
                 Kz_m2s(i,j,:), z_mid(i,j,:),  &
                 z_bnd(i,j,:),  fh(i,j,nr),  &
                 th(i,j,:,nr),  exnm(i,j,:),  &
                 p_bnd(:), pzpbl(i,j), &       
                 .false.)

                 pzpbl(i,j) = max( PBL_ZiMIN, pzpbl(i,j))  ! Keep old fixed height ZiMin here
                 pzpbl(i,j)  = min( PBL_ZiMAX, pzpbl(i,j))

               enddo
            enddo

          else ! Newer non-TI methods
             if ( HmixMethod == "SbRb" ) then

                 call SeibertRiB_Hmix_3d(&
                  u_mid(1:limax,1:ljmax,:),  &
                  v_mid(1:limax,1:ljmax,:),  &
                  z_mid(1:limax,1:ljmax,:),  &
                  th(1:limax,1:ljmax,:,nr),  &
                  pzpbl(1:limax,1:ljmax))

              else if ( HmixMethod == "JcRb" ) then

                 do i=1,limax
                    do j=1,ljmax

                      call JericevicRiB_Hmix(&
                          u_mid(i,j,:), v_mid(i,j,:),  &
                          z_mid(i,j,:), th(i,j,:,nr),  &
                          pzpbl(i,j))
                      end do
                   end do
              else
                 call CheckStop("Need HmixMethod")
              end if ! end of newer methods

           ! Set limits on Zi
           ! mid-call at k=19 is lowest we can resolve, so set as min
              forall(i=1:limax,j=1:ljmax)
                 pzpbl(i,j) = max( z_mid(i,j,KMAX_MID-1), pzpbl(i,j))
                 pzpbl(i,j) = min( PBL_ZiMAX, pzpbl(i,j) )
              end forall

          end if ! Hmix done

        !======================================================================
        ! Kz choices:

           if ( KzMethod == "JG" ) then  ! Jericevic/Grisogono for both Stable/Unstable
             do k = 2, KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                 Kz_m2s(i,j,k) = JericevicKz( z_bnd(i,j,k), pzpbl(i,j), ustar_nwp(i,j) )
             end do
             end do
             end do

          else  ! Specify unstable, stable separately:

            if ( StableKzMethod == "JG" ) then  ! Jericevic/Grisogono for both Stable/Unstable
              do j=1,ljmax
                 do i=1,limax
                   if ( invL_nwp(i,j) >= OB_invL_LIMIT ) then !neutral and unstable
                     do k = 2, KMAX_MID
                       Kz_m2s(i,j,k) = &
                          JericevicKz( z_bnd(i,j,k), pzpbl(i,j), ustar_nwp(i,j) )
                     end do
                   end if
                 end do
              end do


            else if ( StableKzMethod == "BW" ) then

                 do k = 2, KMAX_MID
                  do j=1,ljmax
                     do i=1,limax
                       if ( invL_nwp(i,j) > 1.0e-10 ) then !stable ! leaves gap near zero
                           Kz_m2s(i,j,k) = BrostWyngaardKz(z_bnd(i,j,k),pzpbl(i,j),&
                                           ustar_nwp(i,j),invL_nwp(i,j)) 
                       end if
                 end do
                 end do
                 end do

            else if ( StableKzMethod == "PB" ) then
                 ! no change
            else
                 call CheckStop("Need StableKzMethod")
            end if ! Stable Kz

             if ( UnstableKzMethod == "OB" ) then
                do j=1,ljmax
                   do i=1,limax
                     if ( invL_nwp(i,j) < OB_invL_LIMIT ) then !neutral and unstable
                        call O_BrienKz ( &     ! Original EMEP method
                          pzpbl(i,j),  z_bnd(i,j,:),  &
                          ustar_nwp(i,j),  invL_nwp(i,j),  &
                          Kz_m2s(i,j,:), .false.)
                     end if
                   end do
                 end do
         
              else
                 call CheckStop("Need UnstableKzMethod")
              end if

           end if  ! Specify unstable, stable separately:
    end if


    !..spatial smoothing of new zi: Need fixed minimum here. 100 m is okay

     call smoosp(pzpbl,PBL_ZiMIN,PBL_ZiMAX)

  !************************************************************************!
  ! test some alternative options for Kz and Hmix
   if( DEBUG_BLM .and. debug_proc .and. modulo( current_date%hour, 3)  == 0 &
         .and. current_date%seconds == 0  ) then

      i = debug_iloc
      j = debug_jloc
      p_bnd(:) = sigma_bnd(:)*(ps(i,j,nr) - PT) + PT

  !************************************************************************!
  ! We test all the various options here. Pass in  data as keyword arguments
  ! to avoid possible errors!

      call Test_BLM( mm=current_date%month, dd=current_date%day, &
           hh=current_date%hour, fH=fh(i,j,nr), &
           u=u_mid(i,j,:),v=v_mid(i,j,:), zm=z_mid(i,j,:), &
           zb=z_bnd(i,j,:), exnm=exnm(i,j,:), Kz=Kz_m2s(i,j,:), &
           Kz_nwp=Kz_nwp(:), invL=invL_nwp(i,j), &
           ustar=ustar_nwp(i,j), th=th(i,j,:,nr), pb=p_bnd(:), zi=pzpbl(i,j))
  !************************************************************************!
      ! if ( USE_MIN_Kz) then
         !hs = min( z_bnd(i,j,KMAX_MID), 0.04*pzpbl(i,j))
         hs = z_bnd(i,j,KMAX_MID)

        stab_h = min( PsiH(hs*invL_nwp(i,j)), 0.9 )
         Kz_min = ustar_nwp(i,j)*KARMAN*hs /( 1 - stab_h  )
        write(*,"(a,10f10.3)") "PSIH ", stab_h, fh(i,j,nr), invL_nwp(i,j), &
             PsiH(hs*invL_nwp(i,j)),Kz_min
      ! end if

    end if ! end of debug extra options


    !***************************************************
    if ( .not. (NWP_Kz .and. foundKz_met) ) then  ! convert to Sigma units

      call Kz_m2s_toSigmaKz (Kz_m2s,roa(:,:,:,nr),ps(:,:,nr),SigmaKz(:,:,:,nr))

    end if
    !***************************************************


  end subroutine BLPhysics

  !c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  subroutine smoosp(f,rmin,rmax)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c    written by Trond Iversen,  modified by Hugo Jakobsen, 080994
    !       parallellized and modified by Peter February 2003
    !c
    !c    Called from: BLPhysics.f
    !c
    !c----------------------------------------------------------------------
    !c
    !c    This routine applies the shapiro filter with s=0.5 and s=-0.5
    !c    to the field f usinh h as a work space also the boundaries
    !c    are smoothed. f contains the smoothed field upon return.
    !c

    !c    Definition of the variables:
    !c
    !c
    !c    f    : data to be smoothed
    !c    iif    : =limax
    !c    jjf    : =ljmax
    !c    h1,h2    : = help variable
    !c    rmin    : min allowed
    !c    rmax    : max allowed
    !c
    implicit none

    real, intent(inout) :: f(MAXLIMAX,MAXLJMAX)
    real, intent(in)    :: rmin,rmax

    real, dimension(MAXLIMAX+4,MAXLJMAX+4) :: h1, h2
    real, dimension(MAXLIMAX,2)            :: f_south,f_north
    real, dimension(MAXLJMAX+2*2,2)        :: f_west,f_east
    real s

    integer  thick
    integer iif,jjf,is,i,j,ii,jj,iifl,jjfl

    iif=limax
    jjf=ljmax

    thick=2  !we fetch 2 neighbors at once, so that we don't need to call
    ! readneighbours twice
    iifl=iif+2*thick
    jjfl=jjf+2*thick

    call readneighbors(f,f_south,f_north,f_west,f_east,thick)

    do j=1,jjf
       jj=j+thick
       do i=1,iif
          ii=i+thick
          h1(ii,jj) = f(i,j)
       enddo
    enddo
    do j=1,thick
       do i=1,iif
          ii=i+thick
          h1(ii,j) = f_south(i,j)
       enddo
    enddo

    do j=1,thick
       jj=j+jjf+thick
       do i=1,iif
          ii=i+thick
          h1(ii,jj) = f_north(i,j)
       enddo
    enddo

    do j=1,jjfl
       do i=1,thick
          h1(i,j) = f_west(j,i)
       enddo
    enddo

    do j=1,jjfl
       do i=1,thick
          ii=i+iif+thick
          h1(ii,j) = f_east(j,i)
       enddo
    enddo

    do j=1,jjfl
       h2(1,j) = 0.
       h2(iifl,j) = 0.
    enddo

    do i=1,iifl
       h2(i,1) = 0.
       h2(i,jjfl) = 0.
    enddo
    !! 44 format(I2,30F5.0)

    do  is=2,1,-1

       s=is-1.5  !s=0,5 s=-0.5
       if(is /= 2)h1=h2

       !..the smoothing

       do  j=2,jjfl-1
          do  i=2,iifl-1
             h2(i,j)=(1.-2.*s+s*s)*h1(i,j)&
                  + 0.5*s*(1.-s)*(h1(i+1,j)+h1(i-1,j)+h1(i,j+1)+h1(i,j-1))  &
                  + s*s*(h1(i+1,j+1)+h1(i-1,j-1)+h1(i+1,j-1)+h1(i-1,j+1))/4.
             h2(i,j) = amax1(h2(i,j),rmin)
             h2(i,j) = amin1(h2(i,j),rmax)
          end do
       end do

    end do


    do j=1,jjf
       jj=j+thick
       do i=1,iif
          ii=i+thick
          f(i,j)=h2(ii,jj)
       enddo
    enddo

  end subroutine smoosp
  !  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  subroutine extendarea(f,h,debug_flag)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !     based upon the smoosp routine
    !c    adapted by David Simpson March 2011
    !     - returns extended array array, reading neighbour procs as needed
    !c----------------------------------------------------------------------

    real, intent(inout) :: f(:,:)
    real, intent(inout) :: h(:,:)
    logical, intent(in), optional :: debug_flag
    logical :: mydebug = .false.

    real, dimension(size(f,1),2)            :: f_south,f_north
    real, dimension(size(f,2)+2*2,2)        :: f_west,f_east

    integer :: thick ! = size(h,1) - size(f,1) ! Caller has to make h > f 
    integer :: iif,jjf,i,j,ii,jj,iifl,jjfl
    if ( present(debug_flag)  ) mydebug = debug_flag

    thick = ( size(h,1) - size(f,1) ) ! Caller has to make h > f 
    iif=size(f,1)
    jjf=size(f,2)

    if( modulo(thick,2) /= 0 ) then
       print *, "ERROR extendarea para,s ", me, iif , jjf, thick
       print *, "ERROR extendarea mod ", modulo(thick,2)
       call StopAll("ERROR extendarea thickness not even!")
    end if
    thick = thick / 2


    ! readneighbours twice
    iifl=iif+2*thick
    jjfl=jjf+2*thick
    if(mydebug .and. MasterProc ) write(*,*) "DEBUG extendarea", iif,jjf,thick

    call readneighbors(f,f_south,f_north,f_west,f_east,thick)

    do j=1,jjf
       jj=j+thick
       do i=1,iif
          ii=i+thick
          h(ii,jj) = f(i,j)
       enddo
    enddo
    do j=1,thick
       do i=1,iif
          ii=i+thick
          h(ii,j) = f_south(i,j)
       enddo
    enddo

    do j=1,thick
       jj=j+jjf+thick
       do i=1,iif
          ii=i+thick
          h(ii,jj) = f_north(i,j)
       enddo
    enddo

    do j=1,jjfl
       do i=1,thick
          h(i,j) = f_west(j,i)
       enddo
    enddo

    do j=1,jjfl
       do i=1,thick
          ii=i+iif+thick
          h(ii,j) = f_east(j,i)
       enddo
    enddo

!    do j=1,jjf
!       jj=j+thick
!       do i=1,iif
!          ii=i+thick
!          f(i,j)=h2(ii,jj)
!       enddo
!    enddo
!
  end subroutine extendarea
  !  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  subroutine readneighbors(data,data_south,data_north,data_west,data_east,thick)


    ! Read data at the other side of the boundaries
    !
    ! thick is the number of gridcells in each direction to be transferred
    ! Note that we also fetch data from processors in the "diagonal"
    ! directions
    !
    ! Written by Peter February 2003
    !
    !Note,
    !The data_west(jj,:)=data(1,j) is not a bug: when there is no west
    !neighbour,
    !the data is simply copied from the nearest points: data_west(jj,:) should
    !be =data(-thick+1:0,j), but since this data does not exist, we
    !put it =data(1,j).


    implicit none

    integer, intent(in) :: thick
    real,intent(in), dimension(MAXLIMAX,MAXLJMAX) ::data
    real,intent(out), dimension(MAXLIMAX,thick) ::data_south,data_north
    real,intent(out), dimension(MAXLJMAX+2*thick,thick) ::data_west,data_east
    real, dimension(MAXLIMAX,thick) ::data_south_snd,data_north_snd
    real, dimension(MAXLJMAX+2*thick,thick) ::data_west_snd,data_east_snd

    integer :: msgnr,info,request_s,request_n,request_e,request_w
    integer :: j,tj,jj,jt

    !check that limax and ljmax are large enough
    call CheckStop(limax < thick, "ERROR readneighbors in Met_ml")
    call CheckStop(ljmax < thick, "ERROR readneighbors in Met_ml")


    msgnr=1

    data_south_snd(:,:)=data(:,1:thick)
    data_north_snd(:,:)=data(:,ljmax-thick+1:ljmax)
    if(neighbor(SOUTH) >= 0 )then
       CALL MPI_ISEND( data_south_snd , 8*MAXLIMAX*thick, MPI_BYTE,&
            neighbor(SOUTH), msgnr, MPI_COMM_WORLD, request_s,INFO)
    endif
    if(neighbor(NORTH) >= 0 )then
       CALL MPI_ISEND( data_north_snd , 8*MAXLIMAX*thick, MPI_BYTE,&
            neighbor(NORTH), msgnr+9, MPI_COMM_WORLD, request_n,INFO)
    endif

    if(neighbor(SOUTH) >= 0 )then
       CALL MPI_RECV( data_south, 8*MAXLIMAX*thick, MPI_BYTE,&
            neighbor(SOUTH), msgnr+9, MPI_COMM_WORLD, MPISTATUS, INFO)
    else
       do tj=1,thick
          data_south(:,tj)=data(:,1)
       enddo
    endif
    if(neighbor(NORTH) >= 0 )then
       CALL MPI_RECV( data_north, 8*MAXLIMAX*thick, MPI_BYTE,&
            neighbor(NORTH), msgnr, MPI_COMM_WORLD, MPISTATUS, INFO)
    else
       do tj=1,thick
          data_north(:,tj)=data(:,ljmax)
       enddo
    endif

    jj=0
    do jt=1,thick
       jj=jj+1
       data_west_snd(jj,:)=data_south(1:thick,jt)
       data_east_snd(jj,:)=data_south(limax-thick+1:limax,jt)
    enddo
    do j=1,ljmax
       jj=jj+1
       data_west_snd(jj,:)=data(1:thick,j)
       data_east_snd(jj,:)=data(limax-thick+1:limax,j)
    enddo
    do jt=1,thick
       jj=jj+1
       data_west_snd(jj,:)=data_north(1:thick,jt)
       data_east_snd(jj,:)=data_north(limax-thick+1:limax,jt)
    enddo

    if(neighbor(WEST) >= 0 )then
       CALL MPI_ISEND( data_west_snd , 8*(MAXLJMAX+2*thick)*thick, MPI_BYTE,&
            neighbor(WEST), msgnr+3, MPI_COMM_WORLD, request_w,INFO)
    endif
    if(neighbor(EAST) >= 0 )then
       CALL MPI_ISEND( data_east_snd , 8*(MAXLJMAX+2*thick)*thick, MPI_BYTE,&
            neighbor(EAST), msgnr+7, MPI_COMM_WORLD, request_e,INFO)
    endif



    if(neighbor(WEST) >= 0 )then
       CALL MPI_RECV( data_west, 8*(MAXLJMAX+2*thick)*thick, MPI_BYTE,&
            neighbor(WEST), msgnr+7, MPI_COMM_WORLD, MPISTATUS, INFO)
    else
       jj=0
       do jt=1,thick
          jj=jj+1
          data_west(jj,:)=data_south(1,jt)
       enddo
       do j=1,ljmax
          jj=jj+1
          data_west(jj,:)=data(1,j)
       enddo
       do jt=1,thick
          jj=jj+1
          data_west(jj,:)=data_north(1,jt)
       enddo
    endif
    if(neighbor(EAST) >= 0 )then
       CALL MPI_RECV( data_east, 8*(MAXLJMAX+2*thick)*thick, MPI_BYTE, &
            neighbor(EAST), msgnr+3, MPI_COMM_WORLD, MPISTATUS, INFO)
    else
       jj=0
       do jt=1,thick
          jj=jj+1
          data_east(jj,:)=data_south(limax,jt)
       enddo
       do j=1,ljmax
          jj=jj+1
          data_east(jj,:)=data(limax,j)
       enddo
       do jt=1,thick
          jj=jj+1
          data_east(jj,:)=data_north(limax,jt)
       enddo
    endif

    if(neighbor(SOUTH) >= 0 )then
       CALL MPI_WAIT(request_s, MPISTATUS,INFO)
    endif
    if(neighbor(NORTH) >= 0 )then
       CALL MPI_WAIT(request_n, MPISTATUS,INFO)
    endif
    if(neighbor(WEST) >= 0 )then
       CALL MPI_WAIT(request_w, MPISTATUS, INFO)
    endif
    if(neighbor(EAST) >= 0 )then
       CALL MPI_WAIT(request_e, MPISTATUS,INFO)
    endif

 end subroutine readneighbors



  !************************************************************************!
  subroutine tkediff (nr)                                             !
    !************************************************************************!
    !                                                                        !
    !    This routine computes vertical eddy diffusivities as a function     !
    !    altitude, height of PBL, and a velocity scale, square root of       !
    !    turbulent kinetic energy (TKE). This is a non-local scheme.         !
    !    The TKE at the surface is diagnosed using scales for horizontaland  !
    !    vertical velocities (ustar and wstar) in the surface layer          !
    !    (Alapaty 2004; Holstag et al. 1990 and Mihailovic et al. 2004)      !
    !    PBL ht is calculated using the EMEP formulation                     !
    !                                                                        !
    !    Written by DT Mihailovic (October 2004)                             !
    !    EMEP polishing and comments: JE Jonson and P Wind                   !
    !************************************************************************!

    implicit none

    !     Local constants
    real   , parameter :: SZKM=1600.     &   ! Constant (Blackadar, 1976)
         ,CKZ=0.001      &   ! Constant (Zhang and Athens, 1982)
         ,REFPR=1.0E+05  &   ! Referent pressure
         ,KZ0LT=1.0E-04  &   ! Constant (Alapaty et al., 1997)
         ,RIC=0.10       &   ! Critical Richardson number
                                ! (Holstlag et al., 1993)
         ,ROVG=RGAS_KG/GRAV        ! Used in Calculation of R-number
    integer, parameter :: KLM =KMAX_MID-1



    !     INPUT
    integer, intent(in) :: nr  ! Number of meteorological stored
    ! in arrays (1 or 2)

    !     OUTPUT
    !     skh(i,j,k,nr) array
    !     Values of the Kz coefficients (eddyz (i,j,k)) are transformed nto
    !     sigma system and then they stored in this array which is later used
    !     in ADVECTION module


    !     Local arrays

    integer, dimension(MAXLIMAX,MAXLJMAX)      :: iblht   ! Level of the PBL top
    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND):: eddyz   ! Eddy coefficients
    ! (m2/s)
    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID):: &
         t_virt    &! Potential temperature (K)
         ,e         &! Kinetic energy with respect to height (m2/s2)
         ,dzq       &! Thickness of sigma interface layers (m)
         ,u_mid     &! Wind speed in x-direction (m/s)
         ,v_mid      ! Wind speed in y-direction (m/s)

    real, dimension(MAXLIMAX,MAXLJMAX,KLM):: &
         dza         ! Thickness of half sigma layers (m)

    real, dimension(MAXLIMAX,MAXLJMAX):: &
         pblht ,    &! PBL (Holstag, 1990) (m)
         h_flux,    &! Sensible heat flux  (W/m2)
         ust_r ,    &! Friction velocity (m/s)
         mol   ,    &! Monin-obukhov length (m)
         wstar       ! Convective velocity (m/s)

    real, dimension(KMAX_BND) :: rib         ! Bulk Richardson number

    real, dimension(KMAX_MID) :: &
         rich,       &! Richardson number
         psi_zi       ! Used in the vertical integration

    real, dimension (10) ::      &
         psi_z       & ! Used for calculating
         , zovh          ! TKE

    !     Local variables
    real dtmp, tog, wssq1, wssq2, wssq, tconv, wss, wst, PSI_TKE,    &
         dusq, dvsq, ri, ss, dthdz, busfc, zvh,                   &
         part1, part2, fract1, fract2, apbl, kz0,    &
         u_s, goth

    integer i, j, k, l, kcbl

    write(*,*)&
         'This routine is not ready. for example ust_r and kcbl are not set!'
    stop

    !     Functions for averaging the vertical turbulent kinetic energy
    !      (Alapaty, 2003)
    data psi_z /0.00,2.00,1.85,1.51,1.48,1.52,1.43,1.10,1.20,0.25/
    data zovh  /0.00,0.05,0.10,0.20,0.40,0.60,0.80,1.00,1.10,1.20/

    !     Store the NMW meteorology and variables derived from its

    !     Change the sign
    h_flux(1:limax,1:ljmax)=-fh(1:limax,1:ljmax,nr)

    !     Avoid devision by zero later in the code

    where (ABS(h_flux(1:limax,1:ljmax))<0.0001) h_flux(1:limax,1:ljmax)=0.0001

    !     Check PBL height   ! strange tests! Negative pzpbl check? From 1 to 100m
    !   - odd!
    do i=1,limax
       do j=1,ljmax
          if(ABS(pzpbl(i,j)) < 1.) then
             pzpbl(i,j)=100.
          endif
       enddo
    enddo

    !     Calculate velocity components in the (h) poits (Arakawa notation)
    do k=1,KMAX_MID
       do i=1,limax
          do j=1,ljmax
             !          u_mid(i,j,k)=0.5*(u_xmj(i-1,j  ,k,nr)+u_xmj(i,j,k,nr))
             !          v_mid(i,j,k)=0.5*(v_xmi(i  ,j-1,k,nr)+v_xmi(i,j,k,nr))

             u_mid(i,j,k)=u_xmj(i,j  ,k,nr)
             v_mid(i,j,k)=v_xmi(i  ,j,k,nr)

          enddo
       enddo
    enddo

    !     Avoid small values
    where (ABS(u_mid(1:limax,1:ljmax,1:KMAX_MID))<0.001) &
         u_mid(1:limax,1:ljmax,1:KMAX_MID)=0.001
    where (ABS(v_mid(1:limax,1:ljmax,1:KMAX_MID))<0.001) &
         v_mid(1:limax,1:ljmax,1:KMAX_MID)=0.001

    !     Initialize eddy difussivity arrays
    eddyz(1:limax,1:ljmax,1:KMAX_MID)=0.

    !     Calculate tickness of the full layers
    dzq(1:limax,1:ljmax,1:KMAX_MID) = z_bnd(1:limax,1:ljmax,1:KMAX_MID)  &
         - z_bnd(1:limax,1:ljmax,2:KMAX_BND)

    !     ... and the half sigma layers
    dza(1:limax,1:ljmax,1:KLM) = z_mid(1:limax,1:ljmax,1:KLM)          &
         - z_mid(1:limax,1:ljmax,2:KMAX_MID)

    !     Calculate virtual temperature

    t_virt(1:limax,1:ljmax,1:KMAX_MID) = th(1:limax,1:ljmax,1:KMAX_MID,nr)  &
         * (1.0+0.622*q(1:limax,1:ljmax,1:KMAX_MID,nr))


    !     Calculate Monin-Obuhkov length   (Garratt, 1994)

    do i=1,limax
       do j=1,ljmax
          u_s = ustar_nwp(i,j)
          mol(i,j) = -(ps(i,j,nr)*u_s*u_s*u_s)/                        &
               (KARMAN*GRAV*h_flux(i,j)*KAPPA)
       enddo
    enddo

    !     Calculate the convective velocity (wstar)
    do i=1,limax
       do j=1,ljmax
          wstar(i,j) = GRAV*h_flux(i,j)*pzpbl(i,j)/rho_surf(i,j)    &
               /CP/th(i,j,KMAX_MID,nr)
          if(wstar(i,j) < 0.) then
             wstar(i,j)=-ABS(wstar(i,j))**(0.3333)
          else
             wstar(i,j)=(wstar(i,j))**(0.3333)
          endif
       enddo
    enddo

    !                            ------------------------------------------>
    !     Start with a long loop ------------------------------------------>
    !                            ------------------------------------------>
    DO  i=1,limax
       DO  j=1,ljmax

          rib(1:KMAX_MID) = 0.0        ! Initialize bulk Richardson number

          part1=ust_r(i,j)*ust_r(i,j)*ust_r(i,j)
          wst=AMAX1(wstar(i,j),1.0E-20)
          part2=0.6*wst*wst*wst
          wss=AMAX1(1.0E-4,(part1+part2))
          wss=EXP(0.333333*ALOG(wss))

          if (h_flux(i,j) < 0.0) then
             tconv=0.0                                ! Holstlag et al. (1990)
          else
             tconv=8.5*h_flux(i,j)/rho_surf(i,j)/CP/wss   ! Conversion to
             ! kinematic flux
          endif

          do k=KMAX_MID,1,-1
             dtmp=t_virt(i,j,k)-t_virt(i,j,KMAX_MID)-tconv
             tog=0.5*(t_virt(i,j,k)+t_virt(i,j,KMAX_MID))/GRAV
             wssq1=u_mid(i,j,k)*u_mid(i,j,k)
             wssq2=v_mid(i,j,k)*v_mid(i,j,k)
             wssq=wssq1+wssq2
             wssq=AMAX1(wssq,1.0E-4)
             rib(k)=z_mid(i,j,k)*dtmp/(tog*wssq)
             if(rib(k).ge.RIC) go to 9001
          enddo
9001      continue

          !     Calculate PBL height according to Holtslag et al. (1993)
          pblht(i,j)=0.
          if(k.ne.KMAX_MID) then
             fract1=(RIC-rib(k+1))/(rib(k)-rib(k+1))
             fract2=1.-fract1
             apbl=z_mid(i,j,k)*fract1
             pblht(i,j)=apbl+z_mid(i,j,k+1)*fract2
             if(pblht(i,j) > z_bnd(i,j,k+1)) then
                kcbl=k
             else
                kcbl=k+1
             endif
          endif
          iblht(i,j)=kcbl

          if(pblht(i,j)<z_bnd(i,j,KMAX_MID)) then
             pblht(i,j)=z_bnd(i,j,KMAX_MID)
             iblht(i,j)=KMAX_MID
          endif


          if(pblht(i,j).le.100.) then              !Minimum of PBL height
             pblht(i,j)=100.
          endif

          !     Find the critical Richardson number (Shir and Borestein, 1976)
          do k=2,iblht(i,j)-1
             rich(k)=0.257*dza(i,j,k)**0.175
          enddo

          !     Free troposphere and cloudy case Kz values estimation
          do k = 2,iblht(i,j)-1
             dusq = (u_mid(i,j,k-1)-u_mid(i,j,k))*(u_mid(i,j,k-1)-u_mid(i,j,k))
             dvsq = (v_mid(i,j,k-1)-v_mid(i,j,k))*(v_mid(i,j,k-1)-v_mid(i,j,k))
             ss = (dusq+dvsq)/(dza(i,j,k-1)*dza(i,j,k-1))+1.E-9
             goth = 2.*GRAV/(t_virt(i,j,k-1)+t_virt(i,j,k))
             dthdz = (t_virt(i,j,k-1)-t_virt(i,j,k))/dza(i,j,k-1)
             ri = goth*dthdz/ss

             !     (Duran and Clemp, 1982)

             kz0 = CKZ*dzq(i,j,k)
             if (ri-rich(k) > 0.) then
                eddyz(i,j,k)=kz0
             else
                eddyz(i,j,k)=kz0+SZKM*SQRT(ss)*(rich(k)-ri)/rich(k)
             endif
             eddyz(i,j,k)=AMIN1(eddyz(i,j,k),100.)
          enddo

          !     Eddy diffusivity coefficients for all regimes in the mixed layer

          do  k=iblht(i,j),KMAX_MID
             if (mol(i,j) < 0.0) then                 !Unstable conditions
                ri=(1.0-15.*z_mid(i,j,k)/mol(i,j))**(-0.25)
                ri=ri/KARMAN/z_mid(i,j,k)
                ri=ri*AMAX1(0.0,pblht(i,j)-z_mid(i,j,k))
                dthdz=ri*ust_r(i,j)**3.
                goth=AMAX1(wstar(i,j),0.0)
                dusq=0.4*goth**3.
                ri=(dthdz+dusq)**(2./3.)
                e(i,j,k)=0.5*ri*(2.6)**(2./3.)        !Moeng and Sullivan (1994)
             else
                ri=z_bnd(i,j,k)/pblht(i,j)               !Stable
                ri=z_mid(i,j,k)/pblht(i,j)               !New
                ri=(1.0-ri)
                ri=AMAX1(0.0,ri)
                ri=(1.0-ri)**1.75
                e(i,j,k)=6.*ust_r(i,j)*ust_r(i,j)*ri  !Lenshow(1988)
             endif

             !     Calculate Ksi function using interpolation in the vertical
             !     Alapaty (2001, 2003)

             zvh=z_mid(i,j,k)/pblht(i,j)
             do l=1,9
                if (zvh > zovh(l).and. zvh < zovh(l+1)) then
                   psi_zi(k)=(psi_z(l+1)-psi_z(l))/(zovh(l+1)-zovh(l))
                   psi_zi(k)=psi_zi(k)*(zvh-zovh(l))
                   psi_zi(k)=psi_zi(k)+psi_z(l)
                   psi_zi(k)=psi_zi(k)/2.0               !Normalized the value
                endif
             enddo
          enddo

          !      Calculate integral for Ksi
          psi_tke=0.
          do k=KMAX_MID,iblht(i,j),-1
             psi_tke=psi_tke+psi_zi(k)*dzq(i,j,k)*sqrt(e(i,j,k))
          enddo

          psi_tke=psi_tke/pblht(i,j)



          do k=iblht(i,j),KMAX_MID          !Calculate coefficients
             goth=psi_tke
             goth=goth*KARMAN*z_mid(i,j,k)
             dthdz=z_mid(i,j,k)/pblht(i,j)
             dthdz=1.0-dthdz
             dthdz=AMAX1(1.0E-2,dthdz)
             if(mol(i,j) > 0.0) then                         !Stable
                goth=sqrt(e(i,j,iblht(i,j)))                 ! Mihailovic (2004)
                goth=goth*KARMAN*z_mid(i,j,k)                 ! -----------------
                dthdz=z_mid(i,j,k)/pzpbl(i,j)                 ! -----------------
                dthdz=1.0-dthdz
                dthdz=AMAX1(1.0E-2,dthdz)
                busfc=0.74+4.7*z_mid(i,j,KMAX_MID)/mol(i,j)
                busfc=AMAX1(busfc,1.0)
                dthdz=dthdz**1.50                                  !test (2004)
                eddyz(i,j,k)=goth*dthdz/busfc
             else
                dthdz=dthdz*dthdz
                busfc=1.0
                eddyz(i,j,k)=goth*dthdz/busfc
             endif
          enddo

          !      Checking procedure
          do k=2,iblht(i,j)-1
             if(eddyz(i,j,k).le.0.0) THEN
                eddyz(i,j,k)= KZ0LT
             endif
          enddo

          !      Avoid phisically unrealistic values
          do k=2,KMAX_MID
             IF(eddyz(i,j,k).le.0.1) then
                eddyz(i,j,k)=0.1
             endif
          enddo

          !     To avoid loss of mass/energy through top of the model
          !     put eddyz (I,J,K) to zero at the last  level from top
          eddyz(i,j,KMAX_BND)=0.0

          !     Calculate eddy coefficients at the interfaces
          do k=2,KMAX_MID
             eddyz(i,j,k)=0.5*(eddyz(i,j,k-1)+eddyz(i,j,k)) !!

             !              if(i.eq.10.and.j.eq.10.) then
             !             if (abs(u_xmj(i,j  ,k,nr)-u_mid(i,j,k)).gt.5.) then
             !
             !       print *,"NEW ",i,j,u_xmj(i,j  ,KMAX_MID,nr),u_mid(i,j,KMAX_MID)
             !              endif
          enddo

          !     Transform values of the eddy coeficients into the the sigma coordinate

          do k=2,KMAX_MID
             eddyz(i,j,k)=eddyz(i,j,k)*((sigma_mid(k)-sigma_mid(    k-1))/   &
                  (    z_mid(i,j,k)-z_mid(i,j,k-1)))**2.

          enddo

       ENDDO                  !---------------------------------------->
    ENDDO                  !---------------------------------------->
    !---------------------------------------->

    !     Store diffusivity coefficients into skh(i,j,k,nr) array
    do k=2,KMAX_MID
       do i=1,limax
          do j=1,ljmax
             SigmaKz(i,j,k,nr)=eddyz(i,j,k)
          enddo
       enddo
    enddo

    ! For plotting set pblht  =  pzpbl

    pzpbl(:,:) = pblht(:,:)

  end subroutine tkediff
  !---------------------------------------------------------------


  subroutine Getmeteofield(meteoname,namefield,nrec,&
       ndim,unit,validity,field)
    !
    ! Read the meteofields and distribute to nodes
    !


    implicit none

    real, dimension(*),intent(out)  :: field ! dimensions: (MAXLIMAX,MAXLJMAX)
    ! or     (MAXLIMAX,MAXLJMAX,KMAX)

    character (len = *),intent(in)  ::meteoname,namefield
    character (len = *),intent(out) ::unit,validity
    integer,intent(in)              :: nrec,ndim

    integer*2 :: var_local(MAXLIMAX,MAXLJMAX,KMAX_MID)
    integer*2, allocatable ::var_global(:,:,:)   ! faster if defined with
    ! fixed dimensions for all
    ! nodes?
    real :: scalefactors(2)
    integer :: KMAX,ijk,i,k,j,nfetch

    validity=''

    if(ndim==3)KMAX=KMAX_MID
    if(ndim==2)KMAX=1
    if(MasterProc)then
       allocate(var_global(GIMAX,GJMAX,KMAX))
       nfetch=1
       call GetCDF_short(namefield,meteoname,var_global,GIMAX,IRUNBEG,GJMAX, &
            JRUNBEG,KMAX,nrec,nfetch,scalefactors,unit,validity)
    else
       allocate(var_global(1,1,1)) !just to have the array defined
    endif

    !note: var_global is defined only for me=0
    call global2local_short(var_global,var_local,MSG_READ4,GIMAX,GJMAX,&
         KMAX,1,1)
    CALL MPI_BCAST(scalefactors,8*2,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(validity,50,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(unit,50,MPI_BYTE,0,MPI_COMM_WORLD,INFO)


    deallocate(var_global)


    ijk=0
    do k=1,KMAX ! KMAX is =1 for 2D arrays
       do j=1,MAXLJMAX
          do i=1,MAXLIMAX
             ijk=ijk+1
             field(ijk)=var_local(i,j,k)*scalefactors(1)+scalefactors(2)
          enddo
       enddo
    enddo

  end subroutine Getmeteofield



  subroutine GetCDF_short(varname,fileName,var,GIMAX,IRUNBEG,GJMAX,JRUNBEG &
       ,KMAX,nstart,nfetch,scalefactors,unit,validity)
    !
    ! open and reads CDF file
    !
    ! The nf90 are functions which return 0 if no error occur.
    ! check is only a subroutine which check wether the function returns zero
    !
    !
    implicit none

    character (len=*),intent(in) :: fileName

    character (len = *),intent(in) ::varname
    character (len = *),intent(out) ::unit,validity
    real,intent(out) :: scalefactors(2)
    integer, intent(in) :: nstart,GIMAX,IRUNBEG,GJMAX,JRUNBEG,KMAX
    integer, intent(inout) ::  nfetch
    integer*2, dimension(GIMAX*GJMAX*KMAX*NFETCH),intent(out) :: var
    integer :: varID,ndims
    integer :: ncFileID,status
    real :: scale,offset
    character *100 :: period_read

    validity='                                     ' !initialisation
    period_read='                                     ' !initialisation
    scalefactors(1) = 1.0 !default
    scalefactors(2) = 0.  !default

    ndims=3
    if(KMAX==1)ndims=2
    !open an existing netcdf dataset
    call check(nf90_open(path=trim(fileName),mode=nf90_nowrite,ncid=ncFileID))

    !get varID:
    status = nf90_inq_varid(ncid=ncFileID,name=trim(varname),varID=VarID)
    if(status /= nf90_noerr)then
       validity=field_not_found
       var=0
       call check(nf90_close(ncFileID))
       return
    endif


    !get scale factors

    status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
    if(status == nf90_noerr) scalefactors(1) = scale
    status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
    if(status == nf90_noerr) scalefactors(2) = offset

    !find unit
    unit='                                                         '
    status = nf90_get_att(ncFileID, VarID, "units", unit  )
    if(status /= nf90_noerr)then
       unit='unknown' !default
    endif
    !find validity
    status = nf90_get_att(ncFileID, VarID, "validity", period_read  )
    if(status == nf90_noerr)then
       validity  = trim(period_read)
    else
       status = nf90_get_att(ncFileID, VarID, "period_of_validity", &
            period_read  )
       if(status /= nf90_noerr)then
          validity='instantaneous' !default
       endif
    endif

    ! if(Nfetch<nrecords)then
    !    write(*,*)'Reading record',nstart,' to ',nstart+nfetch-1
    ! endif

    !get variable
    if(ndims==2)then
       call check(nf90_get_var(ncFileID, VarID, var,&
            start=(/IRUNBEG,JRUNBEG,nstart/),count=(/ GIMAX,GJMAX,nfetch /)))
    elseif(ndims==3)then
       call check(nf90_get_var(ncFileID, VarID, var,&
            start=(/IRUNBEG,JRUNBEG,1,nstart/),count=(/GIMAX,GJMAX,KMAX,nfetch /)))
    endif

    call check(nf90_close(ncFileID))

  end subroutine GetCDF_short





  subroutine Getgridparams(meteoname,GRIDWIDTH_M,xp,yp,fi,xm_i,xm_j,xm2,&
       ref_latitude,sigma_mid,Nhh,nyear,nmonth,nday,nhour,nhour_first&
       ,cyclicgrid)
    !
    ! Get grid and time parameters as defined in the meteo file
    ! Do some checks on sizes and dates
    !
    ! This routine is called only once (and is therefore not optimized for speed)
    !

    implicit none

    character (len = *), intent(in) ::meteoname
    integer, intent(in):: nyear,nmonth,nday,nhour
    real, intent(out) :: GRIDWIDTH_M,xp,yp,fi, ref_latitude,&
         xm2(0:MAXLIMAX+1,0:MAXLJMAX+1)&
         ,xm_i(0:MAXLIMAX+1,0:MAXLJMAX+1)&
         ,xm_j(0:MAXLIMAX+1,0:MAXLJMAX+1),sigma_mid(KMAX_MID)
    integer, intent(out):: Nhh,nhour_first,cyclicgrid

    integer :: nseconds(1),n1,i,j,im,jm,i0,j0
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID
    integer :: GIMAX_file,GJMAX_file,KMAX_file,ihh,ndate(4)
    real,dimension(-1:GIMAX+2,-1:GJMAX+2) ::xm_global,xm_global_j,xm_global_i
    integer :: status,iglobal,jglobal,info,South_pole,North_pole,Ibuff(2)
    real :: ndays(1),x1,x2,x3,x4
    character (len = 50) :: timeunit

    if(MasterProc)then
       print *,'Defining grid properties from ',trim(meteoname)
      !open an existing netcdf dataset
       status = nf90_open(path=trim(meteoname),mode=nf90_nowrite,ncid=ncFileID)

       if(status /= nf90_noerr) then
          print *,'not found',trim(meteoname)
          METEOfelt=1
       else
!          print *,'  reading ',trim(meteoname)
          projection=''
         call check(nf90_get_att(ncFileID,nf90_global,"projection",projection))
          if(trim(projection)=='Rotated_Spherical'.or.trim(projection)=='rotated_spherical'&
               .or.trim(projection)=='rotated_pole'.or.trim(projection)=='rotated_latitude_longitude')then
             projection='Rotated_Spherical'
          endif
          write(*,*)'projection: ',trim(projection)

          !get dimensions id
          if(trim(projection)=='Stereographic') then
             call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
             call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
          elseif(trim(projection)==trim('lon lat')) then
             call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
             call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
          else
             !     write(*,*)'GENERAL PROJECTION ',trim(projection)
             call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
             call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
          endif

          call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))
          call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeVarID))

          !get dimensions length
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_file))
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_file))
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_file))
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=Nhh))

          write(*,*)'dimensions meteo grid:',GIMAX_file,GJMAX_file,KMAX_file,Nhh

          if(GIMAX_file/=IIFULLDOM.or.GJMAX_file/=JJFULLDOM)then
             write(*,*)'IIFULLDOM,JJFULLDOM ',IIFULLDOM,JJFULLDOM
             write(*,*)'WARNING: large domain and meteorology file have different sizes'
             write(*,*)'WARNING: THIS CASE IS NOT TESTED. Please change large domain'
          endif

          call CheckStop(GIMAX+IRUNBEG-1 > GIMAX_file, "NetCDF_ml: I outside domain" )
          call CheckStop(GJMAX+JRUNBEG-1 > GJMAX_file, "NetCDF_ml: J outside domain" )
          call CheckStop(KMAX_MID > KMAX_file,  "NetCDF_ml: wrong vertical dimension")
          call CheckStop(24/Nhh, METSTEP,          "NetCDF_ml: METSTEP != meteostep" )

          call CheckStop(nhour/=0 .and. nhour /=3,&
               "Met_ml/GetCDF: must start at nhour=0 or 3")

          call check(nf90_get_att(ncFileID,timeVarID,"units",timeunit))

          ihh=1
          n1=1
          if(trim(timeunit)==trim("days since 1900-1-1 0:0:0"))then
             write(*,*)'Meteo date in days since 1900-1-1 0:0:0'
             call check(nf90_get_var(ncFileID,timeVarID,ndays,&
                  start=(/ihh/),count=(/n1 /)))
             call nctime2idate(ndate,ndays(1))  ! for printout: msg="meteo hour YYYY-MM-DD hh"
          else
             call check(nf90_get_var(ncFileID,timeVarID,nseconds,&
                  start=(/ihh/),count=(/n1 /)))
             call nctime2idate(ndate,nseconds(1)) ! default
          endif
          nhour_first=ndate(4)


          call CheckStop(ndate(1), nyear,  "NetCDF_ml: wrong meteo year" )
          call CheckStop(ndate(2), nmonth, "NetCDF_ml: wrong meteo month" )
          call CheckStop(ndate(3), nday,   "NetCDF_ml: wrong meteo day" )

          do ihh=1,Nhh

             if(trim(timeunit)==trim("days since 1900-1-1 0:0:0"))then
                call check(nf90_get_var(ncFileID, timeVarID, ndays,&
                     start=(/ ihh /),count=(/ n1 /)))
                call nctime2idate(ndate,ndays(1))
             else
                call check(nf90_get_var(ncFileID, timeVarID, nseconds,&
                  start=(/ ihh /),count=(/ n1 /)))
                call nctime2idate(ndate,nseconds(1))
             endif
             call CheckStop( mod((ihh-1)*METSTEP+nhour_first,24), ndate(4),  &
                  date2string("NetCDF_ml: wrong meteo hour YYYY-MM-DD hh",ndate))

          enddo


          !get global attributes
          call check(nf90_get_att(ncFileID,nf90_global,"Grid_resolution",GRIDWIDTH_M))
          if(trim(projection)=='Stereographic')then
             call check(nf90_get_att(ncFileID,nf90_global,"ref_latitude",ref_latitude))
             call check(nf90_get_att(ncFileID, nf90_global, "xcoordinate_NorthPole" &
                  ,xp ))
             call check(nf90_get_att(ncFileID, nf90_global, "ycoordinate_NorthPole" &
                  ,yp ))
             call check(nf90_get_att(ncFileID, nf90_global, "fi",fi ))

             call GlobalPosition
          elseif(trim(projection)==trim('lon lat')) then
             ref_latitude=60.
             xp=0.0
             yp=GJMAX
             fi =0.0
             call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
             call check(nf90_get_var(ncFileID, varID, glon_fdom(1:IIFULLDOM,1) ))
             do i=1,IIFULLDOM
                if(glon_fdom(i,1)>180.0)glon_fdom(i,1)=glon_fdom(i,1)-360.0
                if(glon_fdom(i,1)<-180.0)glon_fdom(i,1)=glon_fdom(i,1)+360.0
             enddo
             do j=1,JJFULLDOM
                glon_fdom(:,j)=glon_fdom(:,1)
             enddo
             call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
             call check(nf90_get_var(ncFileID, varID, glat_fdom(1,1:JJFULLDOM) ))
             do i=1,IIFULLDOM
                glat_fdom(i,:)=glat_fdom(1,:)
             enddo
          else
             ref_latitude=60.
             xp=0.0
             yp=GJMAX
             fi =0.0
             if(trim(projection)=='Rotated_Spherical')then
                call check(nf90_get_att(ncFileID,nf90_global,"grid_north_pole_latitude",grid_north_pole_latitude))
                call check(nf90_get_att(ncFileID,nf90_global,"grid_north_pole_longitude",grid_north_pole_longitude))
             endif
             call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
             call check(nf90_get_var(ncFileID, varID, glon_fdom(1:IIFULLDOM,1:JJFULLDOM) ))

             call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
             call check(nf90_get_var(ncFileID, varID, glat_fdom(1:IIFULLDOM,1:JJFULLDOM) ))
             do j=1,JJFULLDOM
             do i=1,IIFULLDOM
                if(glon_fdom(i,j)>180.0)glon_fdom(i,j)=glon_fdom(i,j)-360.0
                if(glon_fdom(i,j)<-180.0)glon_fdom(i,j)=glon_fdom(i,j)+360.0
             enddo
             enddo

          endif
          !get variables
          status=nf90_inq_varid(ncid=ncFileID, name="map_factor", varID=varID)

          if(status == nf90_noerr)then
             !mapping factor at center of cells is defined
             !make "staggered" map factors
             call check(nf90_get_var(ncFileID, varID, xm_global(1:GIMAX,1:GJMAX) &
                  ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
             do j=1,GJMAX
                do i=1,GIMAX-1
                   xm_global_j(i,j)=0.5*(xm_global(i,j)+xm_global(i+1,j))
                enddo
             enddo
             i=GIMAX
             do j=1,GJMAX
                xm_global_j(i,j)=1.5*xm_global(i,j)-0.5*xm_global(i-1,j)
             enddo
             do j=1,GJMAX-1
                do i=1,GIMAX
                   xm_global_i(i,j)=0.5*(xm_global(i,j)+xm_global(i,j+1))
                enddo
             enddo
             j=GJMAX
             do i=1,GIMAX
                xm_global_i(i,j)=1.5*xm_global(i,j)-0.5*xm_global(i,j-1)
             enddo

          else
             !map factor are already staggered
             status=nf90_inq_varid(ncid=ncFileID, name="map_factor_i", varID=varID)

             !BUGCHECK - moved here... (deleted if loop)
             call CheckStop( status, nf90_noerr, "erro rreading map factor" )

             call check(nf90_get_var(ncFileID, varID, xm_global_i(1:GIMAX,1:GJMAX) &
                  ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
             call check(nf90_inq_varid(ncid=ncFileID, name="map_factor_j", varID=varID))
             call check(nf90_get_var(ncFileID, varID, xm_global_j(1:GIMAX,1:GJMAX) &
                  ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
          endif

          call check(nf90_inq_varid(ncid = ncFileID, name = "k", varID = varID))
          call check(nf90_get_var(ncFileID, varID, sigma_mid ))


       endif !found meteo
    endif !me=0

    CALL MPI_BCAST(METEOfelt ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    if( METEOfelt==1)return !do not use NetCDF meteo input


    CALL MPI_BCAST(Nhh,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(GRIDWIDTH_M,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(ref_latitude,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xp,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(yp,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(fi,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(sigma_mid,8*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xm_global_i(1:GIMAX,1:GJMAX),8*GIMAX*GJMAX,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xm_global_j(1:GIMAX,1:GJMAX),8*GIMAX*GJMAX,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(glat_fdom(1:IIFULLDOM,1:JJFULLDOM),8*IIFULLDOM*JJFULLDOM,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(glon_fdom(1:IIFULLDOM,1:JJFULLDOM),8*IIFULLDOM*JJFULLDOM,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(projection,len(projection),MPI_CHARACTER,0,MPI_COMM_WORLD,INFO) 


    do j=1,MAXLJMAX
       do i=1,MAXLIMAX
          glon(i,j)=glon_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
          glat(i,j)=glat_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
       enddo
    enddo
    i0=0
    im=MAXLIMAX
    j0=0
    jm=MAXLJMAX
    if(gi0+MAXLIMAX+1+IRUNBEG-2>IIFULLDOM)im=MAXLIMAX-1!outside fulldomain
    if(gi0+0+IRUNBEG-2<1)i0=1!outside fulldomain
    if(gj0+MAXLJMAX+1+JRUNBEG-2>JJFULLDOM)jm=MAXLJMAX-1!outside fulldomain
    if(gj0+0+JRUNBEG-2<1)j0=1!outside fulldomain

    do j=j0,jm
       do i=i0,im
          x1=glon_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
          x2=glon_fdom(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)
          x3=glon_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)
          x4=glon_fdom(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2)

!8100=90*90; could use any number much larger than zero and much smaller than 180*180  
          if(x1*x2<-8100.0 .or. x1*x3<-8100.0 .or. x1*x4<-8100.0)then
             !Points are on both sides of the longitude -180=180              
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
             if(x3<0)x3=x3+360.0
             if(x4<0)x4=x4+360.0
          endif
          gl_stagg(i,j)=0.25*(x1+x2+x3+x4)

          gb_stagg(i,j)=0.25*(glat_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               glat_fdom(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               glat_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)+&
               glat_fdom(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2))
       enddo
    enddo
    do j=0,j0
       do i=i0,im
          x1=gl_stagg(i,j+1)
          x2=gl_stagg(i,j+2)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i,j+1)-gb_stagg(i,j+2)
       enddo
    enddo
    do j=jm,MAXLJMAX
       do i=i0,im
          x1=gl_stagg(i,j-1)
          x2=gl_stagg(i,j-2)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i,j-1)-gb_stagg(i,j-2)
       enddo
    enddo
    do j=0,MAXLJMAX
       do i=0,i0
          x1=gl_stagg(i+1,j)
          x2=gl_stagg(i+2,j)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i+1,j)-gb_stagg(i+2,j)
       enddo
    enddo
    do j=0,MAXLJMAX
       do i=im,MAXLIMAX
          x1=gl_stagg(i-1,j)
          x2=gl_stagg(i-2,j)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i-1,j)-gb_stagg(i-2,j)
       enddo
    enddo
!ensure that values are within [-180,+180]]
    do j=0,MAXLJMAX
       do i=0,MAXLIMAX
          if(gl_stagg(i,j)>180.0)gl_stagg(i,j)=gl_stagg(i,j)-360.0
          if(gl_stagg(i,j)<-180.0)gl_stagg(i,j)=gl_stagg(i,j)+360.0
       enddo
    enddo

    !test if the grid is cyclicgrid:
    !The last cell + 1 cell = first cell
    Cyclicgrid=1 !Cyclicgrid
    do j=1,JJFULLDOM
       if(mod(nint(glon_fdom(GIMAX,j)+360+360.0/GIMAX),360)/=&
            mod(nint(glon_fdom(IRUNBEG,j)+360.0),360))then
          Cyclicgrid=0  !not cyclicgrid
       endif
    enddo

    if(MasterProc .and. DEBUG_MET)write(*,*)'CYCLICGRID:',Cyclicgrid

    !complete (extrapolate) along the four lateral sides
    do i=1,GIMAX
       xm_global_j(i,0)=1.0/(2.0/(xm_global_j(i,1))-1.0/(xm_global_j(i,2)))
       xm_global_j(i,-1)=1.0/(2.0/(xm_global_j(i,0))-1.0/(xm_global_j(i,1)))
       xm_global_j(i,GJMAX+1)=1.0/(2.0/(xm_global_j(i,GJMAX))-1.0/(xm_global_j(i,GJMAX-1)))
       xm_global_j(i,GJMAX+2)=1.0/(2.0/(xm_global_j(i,GJMAX+1))-1.0/(xm_global_j(i,GJMAX)))
       xm_global_i(i,0)=1.0/(2.0/(xm_global_i(i,1))-1.0/(xm_global_i(i,2)))
       xm_global_i(i,-1)=1.0/(2.0/(xm_global_i(i,0))-1.0/(xm_global_i(i,1)))
       xm_global_i(i,GJMAX+1)=1.0/(2.0/(xm_global_i(i,GJMAX))-1.0/(xm_global_i(i,GJMAX-1)))
       xm_global_i(i,GJMAX+2)=1.0/(2.0/(xm_global_i(i,GJMAX+1))-1.0/(xm_global_i(i,GJMAX)))
    enddo
    do j=-1,GJMAX+2
       xm_global_j(0,j)=1.0/(2.0/(xm_global_j(1,j))-1.0/(xm_global_j(2,j)))
       xm_global_j(-1,j)=1.0/(2.0/(xm_global_j(0,j))-1.0/(xm_global_j(1,j)))
       xm_global_j(GIMAX+1,j)=1.0/(2.0/(xm_global_j(GIMAX,j))-1.0/(xm_global_j(GIMAX-1,j)))
       xm_global_j(GIMAX+2,j)=1.0/(2.0/(xm_global_j(GIMAX+1,j))-1.0/(xm_global_j(GIMAX,j)))
       xm_global_i(0,j)=1.0/(2.0/(xm_global_i(1,j))-1.0/(xm_global_i(2,j)))
       xm_global_i(-1,j)=1.0/(2.0/(xm_global_i(0,j))-1.0/(xm_global_i(1,j)))
       xm_global_i(GIMAX+1,j)=1.0/(2.0/(xm_global_i(GIMAX,j))-1.0/(xm_global_i(GIMAX-1,j)))
       xm_global_i(GIMAX+2,j)=1.0/(1.0/(2*xm_global_i(GIMAX+1,j))-1.0/(xm_global_i(GIMAX,j)))
    enddo

    j=1
    i=1
    if(abs(1.5*glat_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)-0.5*glat_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2))>89.5)then
       write(*,*)'south pole' !xm is infinity
       xm_global_i(:,0)=1.0E19
       xm_global_i(:,-1)=1.0E19
    endif



    !keep only part of xm relevant to the local domain
    !note that xm has dimensions larger than local domain

    call CheckStop( MAXLIMAX+1 > limax+2, "Error in Met_ml X size definition" )
    call CheckStop( MAXLJMAX+1 > ljmax+2, "Error in Met_ml J size definition" )

    do j=0,MAXLJMAX+1
       do i=0,MAXLIMAX+1
          iglobal=gi0+i-1
          jglobal=gj0+j-1
          xm_i(i,j)=xm_global_i(iglobal,jglobal)
          xm_j(i,j)=xm_global_j(iglobal,jglobal)
          !Note that xm is inverse length: interpolate 1/xm rather than xm
          xm2(i,j) = 4.0*( (xm_global_i(iglobal,jglobal-1)*&
               xm_global_i(iglobal,jglobal))/ &
               (xm_global_i(iglobal,jglobal-1)+&
               xm_global_i(iglobal,jglobal)    )   )   *(&
               xm_global_j(iglobal-1,jglobal)*&
               xm_global_j(iglobal,jglobal)     )/(&
               xm_global_j(iglobal-1,jglobal)+&
               xm_global_j(iglobal,jglobal)     )
       enddo
    enddo

    !pw
    !If some cells are to narrow (Poles in lat lon coordinates),
    !this will give too small time steps in the Advection,
    !because of the constraint that the Courant number should be <1.
    !
    !If Poles are found and lon-lat coordinates are used the Advection scheme
    !will be modified to be able to cope with the singularity

    !Look for poles
    !If the northernmost or southernmost lines are poles, they are not
    !considered as outer boundaries and will not be treated 
    !by "BoundaryConditions_ml".
    !Note that "Poles" is defined in subdomains

    North_pole=1
    do i=1,limax
       if(nint(glat(i,ljmax))<=88)then
          North_pole=0  !not north pole
       endif
    enddo

    South_pole=1
    do i=1,limax
       if(nint(glat(i,1))>=-88)then
          South_pole=0  !not south pole
       endif
    enddo

    Poles=0
    if(North_pole==1)then
       Poles(1)=1
       write(*,*)me,'Found North Pole'
    endif

    if(South_pole==1)then
       Poles(2)=1
       write(*,*)me,'Found South Pole'
    endif

     CALL MPI_ALLREDUCE(Poles,Ibuff,2,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,INFO)
     Pole_included=max(Ibuff(1),Ibuff(2))
     
     if(ME==0.and.Pole_included==1)write(*,*)'The grid includes a pole'
     
  end subroutine Getgridparams





  subroutine check(status)
    implicit none
    integer, intent ( in) :: status

    call CheckStop( status, nf90_noerr, "Error in Met_ml/NetCDF stuff" &
         //  trim( nf90_strerror(status) ) )

  end subroutine check
end module met_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
!  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



