! < ForestFire_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

module ForestFire_ml
 !----------------------------------------------------------------
 ! Uses emissions from GFED 3 (Global Forest Emission database)
 ! http://www.falw.vu/~gwerf/GFED/
 ! Currently programmed for 8-daily data (available for 2001 - 2007)
 !----------------------------------------------------------------
  use CheckStop_ml,      only : CheckStop
  use ChemChemicals_ml,  only : species
  use ChemSpecs_tot_ml,  only : NO, CO
  use Country_ml,        only : IC_BB   ! FFIRE
  use EmisDef_ml,        only : ISNAP_NAT ! Fires are assigned to SNAP-11 usually
  use My_Emis_ml,        only :  &
     NEMIS_FILES &
    ,EMIS_NAME  ! lets us know which pollutants are wanted, e.g. sox, pm25

  use EmisGet_ml,        only : &
         nrcemis, nrcsplit, emisfrac &  ! speciation routines and array
        ,iqrc2itot                   &  !maps from split index to total index
        ,emis_nsplit     ! No. spec per file, e.g. nox has 2, for NO and NO2

  use GridValues_ml,     only : i_fdom, j_fdom, debug_li, debug_lj, debug_proc
  use Io_ml,             only : PrintLog
  use MetFields_ml,      only : z_bnd
  use ModelConstants_ml, only : MasterProc, KMAX_MID, &
                                USE_FOREST_FIRES, DEBUG_FORESTFIRE, &
                                IOU_INST,IOU_HOUR,IOU_HOUR_MEAN, IOU_YEAR
  use NetCDF_ml,         only : ReadField_CDF, Out_netCDF,  Real4 ! Reads, writes 
  use OwnDataTypes_ml,   only : Deriv, TXTLEN_SHORT
  use Par_ml,            only : MAXLIMAX, MAXLJMAX, li0, li1, lj0, lj1, me
  use PhysicalConstants_ml, only : AVOG
  use ReadField_ml,      only : ReadField    ! Reads ascii fields
  use Setup_1dfields_ml, only : rcemis
  use SmallUtils_ml,     only : find_index
  use TimeDate_ml,only : nydays, nmdays, date, current_date   ! No. days per year, date-type 
implicit none

!  Unimod calls just call Fire_Emis(daynumber)
!  and put the day-testing code here. This lets the module decide if new
!  emissions are needed, and keeps all forest-fire logic here
!

  public :: Fire_Emis
  public :: Fire_rcemis
  private :: Export_FireNc

  logical, public, dimension(MAXLIMAX,MAXLJMAX), save ::  burning
  real, private, allocatable, dimension(:,:,:), save :: BiomassBurningEmis

  integer, private, save ::  ieCO  ! index for CO

  logical, private, save, dimension(NEMIS_FILES) ::  fires_found 

  real, private, allocatable, dimension(:), save ::   unitsfac

  !/ We use some integers from the general EMEP emission system:

  integer, private :: &
      iem     &! index for emis file, e.g. sox=1,nox=2
     ,iqrc    &! index of species among speciated, e.g. SO2=1, SO4=2,NO=3 etc.
     ,itot     ! index of species in xn_2d. Use iqrc2itot array to map

   character(len=TXTLEN_SHORT), private :: emep_poll, gfed_poll

   type, private :: BB_Defs
      character(len=TXTLEN_SHORT) :: emep ! e.g.  nox
      character(len=TXTLEN_SHORT) :: gfed ! e.g. NOx
      real                        :: MW    ! mol wt. assumed in emission file
   end type

  ! GFED table ============================================================
    integer, private, parameter :: NDEFINED_EMEP  = 16 ! No pollutants in file

   !/ Defintions of GFED data. If known, we assign the GFED pollutant which
   !  corresponds to each possible EMEP emission file. Simply add EMEP 
   !  lines as required - be consistent with EmisDefs though. (We can 
   !  have more definitions than used in EmisDefs, but not vice.versa.

   ! Assign mol. wts of the GFED data  where known. If mol. wt set to
   ! zero, the code in Fire-rcemis will use the values from the 
   ! ChemSpecs_ml, species()%molwt.
   ! 
   ! If GFED doesn't have emissionss, set a "-" for GFED, then the
   ! desired emission factor (g/kg DW), and then follow the
   ! example in Fire_setups for NH3:
   !

    type(BB_Defs), private, dimension(NDEFINED_EMEP) :: gfed_defs = (/ &
      BB_Defs("sox   ", "SO2   ", 64.0 ), & 
      BB_Defs("co    ", "CO    ", 28.0 ), &  
      BB_Defs("pm25  ", "PM25  ", 0    ), & ! species(PPM25)%molwt ), & 
      BB_Defs("nox   ", "NOx   ", 30.0 ), & ! as NO in GFED, assign 100% in emissplit
      BB_Defs("nh3   ", "-     ",  1.0 ), & ! NH3 not available in GFED. Use 1 g/kg DW
      BB_Defs("pocffl", "-     ",  0.0 ), & ! rb: is this really needed?
      BB_Defs("poccfl", "-     ",  0.0 ), & ! rb: is this really needed?
      BB_Defs("pocfwd", "-     ",  0.0 ), & ! rb: is this really needed?
      BB_Defs("eccwd ", "-     ",  0.0 ), & ! rb: is this really needed?
      BB_Defs("ecfwd ", "-     ",  0.0 ), & ! rb: is this really needed?
      BB_Defs("ecffl ", "-     ",  0.0 ), & ! rb: is this really needed?
      BB_Defs("eccfl ", "-     ",  0.0 ), & ! rb: is this really needed?
      BB_Defs("voc   ", "NMHC  ", 0    ), &  
      BB_Defs("forfbc", "BC    ", 12.0 ), &
      BB_Defs("forfoc", "OC    ", 0    ), &
      BB_Defs("pmco  ", "TPM   ", 0    ) /)   ! nearest. QUERY pm25<<pmco though
  ! =======================================================================

contains
   subroutine Fire_Emis(daynumber)

!.....................................................................
!**    DESCRIPTION:

!    Reads forest-fire emissions. So far set up for GFED 8d, but in
!    principal we can re-code by simply adding alternative
!    subroutines, e.g. to cope with pre-2001 monthly emissions

    integer, intent(in) :: daynumber   ! daynumber (1,..... ,365)


    real    :: rdemis(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
    integer :: nstart, alloc_err
    logical :: my_first_call = .true.   ! DSFF
    integer :: dd_old = -1,  n

!// Input emissions are monthly RETRO [kg/m2/s], GFED [g/m2/8days]
    real, save :: to_kgm2s = 1.0e-3 /(8*24.0*60.0*60.0)

    if (current_date%year<2001) then
        if( my_first_call .and. MasterProc  ) then
           call PrintLog("NO 8d GFED FOREST FIRES BEFORE 2001")
        end if
        my_first_call = .false.
        return
    end if
    if (current_date%year>2007) then
        if( my_first_call .and. MasterProc  ) then
           call PrintLog("NO 8d GFED FOREST FIRES AFTER 2007")
        end if
        my_first_call = .false.
        return
    end if

    if ( DEBUG_FORESTFIRE .and. MasterProc ) then 
        write(*,*) "Into the FIRE days:", current_date%year, &
             daynumber, dd_old, mod ( daynumber, 8 ), my_first_call
    end if

    if (dd_old == daynumber) return   ! Only calculate once per day max


   ! Fire emissions are called at 8 days intervals (1, 9, 17, ....)
   ! 46 values available each year: day 361 is the last one.
   ! Return unless new period

    if ( .not. my_first_call .and. mod ( daynumber, 8 ) /= 1  ) return
    dd_old= daynumber


    nstart=(current_date%year-2001)*46+(daynumber+7)/8

    if(DEBUG_FORESTFIRE .and. MasterProc) &
            write(*,*) "FOREST_FIRE: ", daynumber,nstart

    ! We need to look for forest-fire emissions which are equivalent
    ! to the standard emission files:

    ieCO = -999
    fires_found(:) = .false.

    do iem = 1, NEMIS_FILES

       emep_poll = EMIS_NAME(iem)
       n = find_index(emep_poll, gfed_defs(:)%emep )
       gfed_poll = gfed_defs(n)%gfed

       if(DEBUG_FORESTFIRE .and. MasterProc) then
          write(*,"(a,i3,1x,2a8,2i3,a)") "FIRE SETUP: ", &
            iem, trim(emep_poll), trim(gfed_poll), len_trim(gfed_poll) &
           ,n,  trim(gfed_defs(n)%gfed)
       end if
       
       if ( len_trim(gfed_poll) > 1 ) then

         fires_found(iem) = .true.

         call ReadField_CDF('GLOBAL_ForestFireEmis.nc',gfed_poll,&
              rdemis,nstart,interpol='zero_order',needed=.true.)

         if ( my_first_call ) then ! Assume NEMIS_FILES for now
             allocate(BiomassBurningEmis(NEMIS_FILES,MAXLIMAX,MAXLJMAX),&
                          stat=alloc_err)
             call CheckStop( alloc_err, "BB alloc problem")
  
             call Fire_setup()   ! Gets InvMolwWtFac

             my_first_call = .false.

         end if

        !/ CO is special. Keep the index
         if ( trim(gfed_poll) == "CO" ) ieCO = iem

        ! Assign and convert units: GFED [g/m2/month]->[kg/m2/s]

         BiomassBurningEmis(iem,:,:) = rdemis(:,:) * to_kgm2s 

         call PrintLog("ForestFire_ml :: Assigns " // &
                         trim(gfed_poll) , MasterProc)

       else
          call PrintLog("ForestFire_ml :: No GFED emis for " // &
                         trim(gfed_poll) , MasterProc)
       end if
    end do

   !/ If GFED doesn't have emissions, we create them from CO

    do iem = 1, NEMIS_FILES
       emep_poll = EMIS_NAME(iem)
       n = find_index(emep_poll, gfed_defs(:)%emep )
       gfed_poll = gfed_defs(n)%gfed

       if ( gfed_poll == "-" ) then ! Use CO. unitsfac will convert later

         BiomassBurningEmis(iem,:,:) = BiomassBurningEmis(ieCO,:,:)

         fires_found(iem) = .true. !??? not really used yet
         call PrintLog("ForestFire_ml :: Estimates " // trim(emep_poll), &
                MasterProc)
            
       end if
    end do

   !/ Logical to let Unimod know if there is any emission here to
   !  worry about

    burning(:,:) =  ( BiomassBurningEmis(ieCO,:,:) > 1.0e-19 )


  end subroutine Fire_Emis
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Fire_setup()

   ! Pre-calculate conversion factors to get from BiomassBurning's kg/m2/s 
   ! to molecules/cm3/s. An important array (assigned elswhere) is emisfrac
   ! which assigns species such as NOx or VOC to NO, NO2, C3H6 etc.. These
   ! values must be set in the emissplit.specials. files if different from
   ! the default SNAP-11 speciation.
   ! 
   ! We need to assign the correct mol. wt., sometimes from GFED assumptions,
   ! sometimes from EMEP species.
   !
   ! We also handle the case where GFED doesn't have emissions, but an 
   ! emission factor can be assumed (e.g. NH3).

    integer :: ie, f, n, alloc_err
    real :: efCO = 100.0   ! Emission factor of CO for scaling, g/mg DW

    allocate( unitsfac(nrcemis), stat=alloc_err)
    call CheckStop( alloc_err, "BB MWF alloc problem")

     iqrc = 0   ! index over emisfrac

     do ie = 1, NEMIS_FILES

       emep_poll = EMIS_NAME(ie)
       n = find_index(emep_poll, gfed_defs(:)%emep ) ! row in gfed table
       gfed_poll = gfed_defs(n)%gfed

       do f = 1, emis_nsplit( ie )

           iqrc = iqrc + 1
           itot = iqrc2itot(iqrc)  !index in xn_2d array

           if ( len_trim(gfed_poll) > 1 ) then

              if ( gfed_defs(n)%MW > 0 ) then ! use GFED's MW

                unitsfac(iqrc) = emisfrac(iqrc,ISNAP_NAT,IC_BB) / &
                                     gfed_defs(n)%MW

              else ! use EMEP model's MW 

                unitsfac(iqrc) = emisfrac(iqrc,ISNAP_NAT,IC_BB) / &
                                     species(itot)%molwt

              end if

           else if ( gfed_poll == "-" ) then

         ! Factors to get from CO to emissions of other species, here NH3
         !
         ! GFED assumes CO  emission is ca. efCO = 100 g/kg(DM) for extra-trop forest
         ! Andreae+Merlet 2001 have NH3 emission of ca. 1 g/kg as NH3

                unitsfac(iqrc) = &
                    gfed_defs(n)%MW/efCO  & ! When "-", MW is really emis factor
                  * emisfrac(iqrc,ISNAP_NAT,IC_BB) / &
                    species(itot)%molwt

           else 

              call CheckStop( gfed_poll, "GFED case not found")

           end if
           if(DEBUG_FORESTFIRE .and. MasterProc) then
              write(*,"(a,3i3,1x,3a8,2es10.3)") "ForestFire_ml :: Setup-fac " , &
                iem, f, iqrc,  trim(emep_poll), trim(gfed_poll), &
                   trim(species(itot)%name), &
                   emisfrac(iqrc,ISNAP_NAT,IC_BB), unitsfac(iqrc)
           end if
        end do ! f and iqrc
     end do ! ie


  !// And one final conversion factor.
  !// fires [kg/m2/s] -> [kg/m3/s] -> [molec/cm3/s] (after division by DeltaZ and MW)
  !        1 kg ->  1.0e3 g
  !         /m2 ->  1.0e-6 /cm2
  ! Need MW in g/mole and delta-z in cm

     unitsfac(:) =  unitsfac(:)  * 0.001 * AVOG

  end subroutine Fire_setup
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Fire_rcemis(i,j)

!  Disperses the fire emissions vertically and converts to molecules/cm3/s.

!// Injection height: here over 8 levels. Alternative could be PBL
!   or  equally upto ca. 2*PBL (suggested by Sofiev, GEMS)
!ds QUERY - should the emissions be divided equally by level?
!   - will give a higher mixing ratio for thinner levels

   integer, intent(in) :: i,j

   integer, parameter :: KEMISFIRE = 12
   real, dimension(KEMISFIRE:KMAX_MID) :: invDeltaZfac !  height of layer in m div 9
   integer ::  k, f

   integer, parameter ::  N_LEVELS = KMAX_MID - KEMISFIRE + 1  ! = 9.0 here


   real    :: origrc, bbe
   logical :: debug_flag


     debug_flag = ( DEBUG_FORESTFIRE .and. &
                     debug_proc .and. i == debug_li .and. j == debug_lj ) 

     if ( debug_flag ) then
        write(*,"(a,5i4,es12.3,f9.3)") "Burning ", me, i,j, &
              i_fdom(i), j_fdom(j), BiomassBurningEmis(ieCO,i,j)
     end if

    !/ Here we just divide by the number of levels. Biased towards
    !  different levels since thickness and air content differ. Simple though.

     do k = KEMISFIRE, KMAX_MID
       invDeltaZfac(k) = 1.0/ (z_bnd(i,j,k) - z_bnd(i,j,k+1)) /N_LEVELS
     end do

     iqrc = 0   ! index over emisfrac
     EMLOOP : do iem = 1, NEMIS_FILES

        do f = 1, emis_nsplit( iem )

           iqrc = iqrc + 1

           if ( .not. fires_found(iem) ) cycle EMLOOP

           itot = iqrc2itot(iqrc)  !index in xn_2d array

           bbe = BiomassBurningEmis(iem,i,j) * unitsfac( iqrc ) 


           origrc = rcemis( itot, KMAX_MID ) ! just for printout 

           ! distribute vertically:

           do k = KEMISFIRE, KMAX_MID
                rcemis( itot, k ) = rcemis( itot, k )  + bbe * invDeltaZfac(k)
           end do !k

           if ( debug_flag ) then
             k=KMAX_MID
             write(*,"(a,i3,1x,a8,f7.3,i4,es10.2,4es10.2)") "FIRERC ",&
              iem, trim(species(itot)%name), emisfrac(iqrc,ISNAP_NAT,IC_BB), &
                 k, BiomassBurningEmis(iem,i,j),&
                   unitsfac(iqrc), invDeltaZfac(k), origrc, rcemis( itot, k )
           end if

    !--  Add up emissions in ktonne ......
    !   totemadd(itot) = totemadd(itot) + &
    !   tmpemis(iqrc) * dtgrid * xmd(i,j)
  
        end do ! f
     end do  EMLOOP ! iem

  end subroutine Fire_rcemis
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Export_FireNc()
    type(Deriv) :: def1 ! definition of fields
 
    def1%class='ForestFireEmis' !written
    def1%avg=.false.      !not used
    def1%index=0          !not used
    def1%scale=1.0      !not used
!FEB2011    def1%inst=.true.      !not used
!FEB2011    def1%year=.false.     !not used
!FEB2011    def1%month=.false.    !not used
!FEB2011    def1%day=.false.      !not used
    def1%name='NOx'        !written
    def1%unit='g/m2'       !written
    def1%name='NOx_zero'       
    def1%name='CO_ASCII'       

    call Out_netCDF(IOU_INST,def1,2,1, BiomassBurningEmis(ieCO,:,:),1.0,&
           CDFtype=Real4,fileName_given='FF.nc')
  end subroutine Export_FireNc

end module ForestFire_ml

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!ncdump -h /global/work/mifapw/emep/Data/ForestFire/ForestFireEmis.nc:
!---------------------------------------------------------------
! netcdf ForestFireEmis {
! dimensions:
! 	lon = 360 ;
! 	lat = 180 ;
! 	time = UNLIMITED ; // (322 currently)
! variables:
! 	double lon(lon) ;
! 		lon:standard_name = "longitude" ;
! 		lon:long_name = "longitude" ;
! 		lon:units = "degrees_east" ;
! 	double lat(lat) ;
! 		lat:standard_name = "latitude" ;
! 		lat:long_name = "latitude" ;
! 		lat:units = "degrees_north" ;
! 	int time(time) ;
! 		time:units = "days since 1900-1-1 0:0:0" ;
! 	double map_factor_i(lat, lon) ;
! 		map_factor_i:long_name = "mapping factor in i direction" ;
! 		map_factor_i:units = "" ;
! 	double map_factor_j(lat, lon) ;
! 		map_factor_j:long_name = "mapping factor in j direction" ;
! 		map_factor_j:units = "" ;
! 	float PM25(time, lat, lon) ;
! 		PM25:long_name = "PM25" ;
! 		PM25:units = "g/m2/8days" ;
! 		PM25:numberofrecords = 322 ;
! 		PM25:_FillValue = 9.96921e+36f ;
! 	float BC(time, lat, lon) ;
! 		BC:long_name = "BC" ;
! 		BC:units = "g/m2/8days" ;
! 		BC:numberofrecords = 322 ;
! 		BC:_FillValue = 9.96921e+36f ;
! 	float NMHC(time, lat, lon) ;
! 		NMHC:long_name = "NMHC" ;
! 		NMHC:units = "g/m2/8days" ;
! 		NMHC:numberofrecords = 322 ;
! 		NMHC:_FillValue = 9.96921e+36f ;
! 	float CO(time, lat, lon) ;
! 		CO:long_name = "CO" ;
! 		CO:units = "g/m2/8days" ;
! 		CO:numberofrecords = 322 ;
! 		CO:_FillValue = 9.96921e+36f ;
! 	float OC(time, lat, lon) ;
! 		OC:long_name = "OC" ;
! 		OC:units = "g/m2/8days" ;
! 		OC:numberofrecords = 322 ;
! 		OC:_FillValue = 9.96921e+36f ;
! 	float NOx(time, lat, lon) ;
! 		NOx:long_name = "NOx" ;
! 		NOx:units = "g/m2/8days" ;
! 		NOx:numberofrecords = 322 ;
! 		NOx:_FillValue = 9.96921e+36f ;
! 	float SO2(time, lat, lon) ;
! 		SO2:long_name = "SO2" ;
! 		SO2:units = "g/m2/8days" ;
! 		SO2:numberofrecords = 322 ;
! 		SO2:_FillValue = 9.96921e+36f ;
! 	float TPM(time, lat, lon) ;
! 		TPM:long_name = "TPM" ;
! 		TPM:units = "g/m2/8days" ;
! 		TPM:numberofrecords = 322 ;
! 		TPM:_FillValue = 9.96921e+36f ;
! 
! // global attributes:
! 		:Conventions = "CF-1.0" ;
! 		:projection = "lon lat" ;
! 		:vert_coord = "sigma: k ps: PS ptop: PT" ;
! 		:Grid_resolution = 111177.473352039 ;
! 		:created_date = "20091021" ;
! 		:created_hour = "143950.341" ;
! 		:lastmodified_date = "20091021" ;
! 		:lastmodified_hour = "145458.652" ;
! }
