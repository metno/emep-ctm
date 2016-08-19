! <Emissions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

                    module Emissions_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Calls up emission read/set routines
!  This routine interfaces the stand-alone emission-file reading routines
!  with the 3D model.
!_____________________________________________________________________________

  use CheckStop_ml,only : CheckStop
  use ChemSpecs_shl_ml, only: NSPEC_SHL
  use ChemSpecs_tot_ml, only: NSPEC_TOT,NO2
  use ChemChemicals_ml, only: species
  use Country_ml,    only : NLAND,Country_Init,Country, IC_NAT
  use My_Emis_ml, only : NEMIS_FILES & ! No. emission files
                        ,EMIS_NAME    ! Names of species ("sox  ",...)
  use EmisDef_ml, only : NSECTORS & ! No. sectors
                     ,NEMISLAYERS & ! No. vertical layers for emission
                     ,NCMAX       & ! Max. No. countries per grid
                     ,FNCMAX      & ! Max. No. countries (with flat emissions)
                                    ! per grid
                     ,ISNAP_SHIP  & ! snap index for ship emissions
                     ,ISNAP_NAT   & ! snap index for nat. (dms) emissions
                     ,IQ_DMS      & ! code for DMS emissions
                     ,VERTFAC       ! vertical emission split
  use EmisGet_ml, only : EmisGet, EmisSplit, &
         nrcemis, nrcsplit, emisfrac &  ! speciation routines and array
        ,iqrc2itot                   &  ! maps from split index to total index
        ,emis_masscorr               &  ! 1/molwt for most species
        ,emis_nsplit                    ! No. species per emis file
  use GridValues_ml, only:  GRIDWIDTH_M    & ! size of grid (m)
                           ,xm2            & ! map factor squared
                           ,debug_proc,debug_li,debug_lj & 
                           ,sigma_bnd, xmd, glat, glon,dA,dB
  use Io_Nums_ml,       only : IO_LOG, IO_DMS, IO_EMIS
  use Io_Progs_ml,      only : ios, open_file, datewrite
  use MetFields_ml,     only : roa, ps, z_bnd   ! ps in Pa, roa in kg/m3
  use ModelConstants_ml,only : KMAX_MID, KMAX_BND, PT ,dt_advec, &
                              IS_GLOBAL, & 
                              NBVOC,     &      ! > 0 if forest voc wanted
                              DEBUG => DEBUG_EMISSIONS,  MasterProc, & 
                              DEBUG_SOILNO, & 
                              NPROC, IIFULLDOM,JJFULLDOM , & 
                              USE_AIRCRAFT_EMIS, &
                              USE_SOIL_NOX
  use Par_ml,     only : MAXLIMAX,MAXLJMAX,me,gi0,gi1,gj0,gj1, &
                             GIMAX, GJMAX, IRUNBEG, JRUNBEG,  &   
                             limax,ljmax,li0,lj0,li1,lj1, &
                             MSG_READ1,MSG_READ7
  use PhysicalConstants_ml,  only :  GRAV,  AVOG
  use ReadField_ml, only : ReadField    ! Reads ascii fields
  use TimeDate_ml,  only : nydays, nmdays, date, current_date ! No. days per 
                                                              ! year, date-type 
  use Timefactors_ml, only :   &
               NewDayFactors   &         ! subroutines
              ,timefac, day_factor       ! time-factors
  use Timefactors_ml, only : timefactors   &                 ! subroutine
                             ,fac_emm, fac_edd, day_factor   ! time-factors
  use Volcanos_ml


  implicit none
  private


 ! subroutines:

  public :: Emissions                ! Main emissions module 
  public :: newmonth
  public :: EmisSet                  ! Sets emission rates every hour/time-step

  ! The main code does not need to know about the following 
  private :: consistency_check       ! Safety-checks

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO


 ! land-code information in each grid square - needed to know which country
 ! is emitting.                        
 ! nlandcode = No. countries in grid square
 ! landcode  = Country codes for that grid square
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX)       :: nlandcode
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX,NCMAX) :: landcode
 ! for flat emissions, i.e. no vertical extent:
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX)       :: flat_nlandcode
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX,FNCMAX):: flat_landcode

 !
 ! The output emission matrix for the 11-SNAP data is snapemis:
 !

  real, private, dimension(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES) &
            , save ::  snapemis      ! main emission arrays, in kg/m2/s

  real, private, dimension(MAXLIMAX,MAXLJMAX,FNCMAX,NEMIS_FILES) &
            , save ::  snapemis_flat ! main emission arrays, in kg/m2/s  

 ! We store the emissions for output to d_2d files and netcdf in kg/m2/s

  real, public, dimension(MAXLIMAX,MAXLJMAX,NEMIS_FILES), save :: SumSnapEmis

  logical, save, private  :: first_dms_read

  ! Emissions for input to chemistry routines

  ! KEMISTOP added to avoid hard-coded KMAX_MID-3:

   integer, public, parameter :: KEMISTOP = KMAX_MID - NEMISLAYERS + 1
   real, public, allocatable, save, dimension(:,:,:,:) :: &
        gridrcemis      & ! varies every time-step (as ps changes)
       ,gridrcemis0       ! varies every hour

  ! and for budgets (not yet used - not changed dimension)
   real, public,  save, dimension(NSPEC_SHL+1:NSPEC_TOT) ::  totemadd

   real, public, save, dimension(MAXLIMAX,MAXLJMAX) :: SoilNOx
   integer, private, save :: iemCO  ! index of CO emissions, for debug

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
  integer :: err1, err2, err3, err4, err5, err6 ! Error messages
  integer :: fic 
  integer :: ic               ! country codes 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over pollutants (1..NEMIS_FILES)


  ! Emission sums (after e_fact adjustments):
  real, dimension(NEMIS_FILES)       :: emsum    ! Sum emis over all countries
  real, dimension(NLAND,NEMIS_FILES) :: sumemis  ! Sum of emissions per country

  if (MasterProc) write(6,*) "Reading emissions for year",  year

  ! 0) set molwts, conversion factors (e.g. tonne NO2 -> tonne N), and
  !    emission indices (IQSO2=.., )

  !=========================
  call Country_Init()    ! In Country_ml, => NLAND, country codes and 
                         !                   names, timezone
  !=========================

  call consistency_check()               ! Below
  !=========================

  ios = 0

  if( MasterProc) then   !::::::: ALL READ-INS DONE IN HOST PROCESSOR ::::

     write(6,*) "Reading monthly and daily timefactors"
    !=========================
     call timefactors(year)               ! => fac_emm, fac_edd, day_factor
    !=========================

  endif

 !=========================
   call EmisSplit()    ! In EmisGet_ml, => emisfrac
 !=========================
  call CheckStop(ios, "ioserror: EmisSplit")


  ! ####################################
  ! Broadcast  monthly and Daily factors 
    CALL MPI_BCAST( fac_emm ,8*NLAND*12*NSECTORS*NEMIS_FILES,MPI_BYTE,  0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST( fac_edd ,8*NLAND*7*NSECTORS*NEMIS_FILES,MPI_BYTE,   0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST( day_factor ,8*2*NSECTORS,MPI_BYTE,               0,MPI_COMM_WORLD,INFO) 

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


       ! Initialise with 0
       globnland(:,:) = 0      
       flat_globnland(:,:)=0  
       globland(:,:,:) = 0    
       globemis(:,:,:,:) = 0  
       flat_globland(:,:,:)=0 
       globemis_flat(:,:,:) =0

   end if

   do iem = 1, NEMIS_FILES
      ! now again test for me=0
      if ( MasterProc ) then

           ! read in global emissions for one pollutant
           ! *****************
             call EmisGet(iem,EMIS_NAME(iem),IRUNBEG,JRUNBEG,GIMAX,GJMAX, &
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

    end do ! iem = 1, NEMIS_FILES-loop


    if ( MasterProc ) then
        write(unit=6,fmt=*) "Total emissions by countries:"
        write(unit=IO_LOG,fmt=*) "Total emissions by countries:"
        write(unit=6,fmt="(2a4,11x,30a12)")  "  N "," CC ",(EMIS_NAME(iem),iem=1,NEMIS_FILES)
        write(unit=IO_LOG,fmt="(2a4,11x,30a12)") "  N "," CC ",(EMIS_NAME(iem),iem=1,NEMIS_FILES)

        do ic = 1, NLAND
           ccsum = sum( sumemis(ic,:) )
           if ( ccsum > 0.0 ) then
                    write(unit=6,fmt="(i3,1x,a4,3x,30f12.2)") &
                        ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS_FILES)
                    write(unit=IO_LOG,fmt="(i3,1x,a4,3x,30f12.2)")& 
                        ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS_FILES)
           end if
        end do
    end if

    ! now all values are read, snapemis is distributed, globnland and 
    ! globland are ready for distribution
    ! print *, "calling glob2local_int for iem", iem, " me ", me

      call global2local_int(globnland,nlandcode,326, GIMAX,GJMAX,1,1,1)
      call global2local_int(globland, landcode ,326, GIMAX,GJMAX,NCMAX,1,1)

      call global2local_int(flat_globnland,flat_nlandcode,326,&
                            GIMAX,GJMAX,1,1,1)  !extra array
      call global2local_int(flat_globland,flat_landcode,326,&
                            GIMAX,GJMAX,FNCMAX,1,1)

     ! Broadcast volcanoe info derived in EmisGet 

        CALL MPI_BCAST(nvolc,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
        CALL MPI_BCAST(i_volc,4*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
        CALL MPI_BCAST(j_volc,4*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
        CALL MPI_BCAST(emis_volc,8*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 



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
        do iem = 1, NEMIS_FILES
           write(unit=6,fmt="(a15,f12.2)") EMIS_NAME(iem),emsum(iem)
        end do
    endif
              

    do iem = 1, NEMIS_FILES
       conv = tonne_to_kgm2s
       if ( trim(EMIS_NAME(iem)) == "co" )  iemCO = iem  ! save this index
 
        if ( DEBUG .and.  debug_proc .and. iem == iemCO ) then
          write(*,"(a,2es10.3)") "SnapPre:" // trim(EMIS_NAME(iem)), &
                  sum( snapemis (:,debug_li,debug_lj,:,iem) ) &
                 ,sum( snapemis_flat (debug_li,debug_lj,:,iem) )
        end if

       forall (ic=1:NCMAX, j=lj0:lj1, i=li0:li1, isec=1:NSECTORS)
          snapemis (isec,i,j,ic,iem) = &
                 snapemis (isec,i,j,ic,iem) * conv * xm2(i,j)
       end forall

       forall (fic=1:FNCMAX, j=lj0:lj1, i=li0:li1)
          snapemis_flat(i,j,fic,iem) = &
                 snapemis_flat(i,j,fic,iem) * conv * xm2(i,j)
       end forall

        if ( DEBUG .and.  debug_proc .and. iem == iemCO ) then
          write(*,"(a,2es10.3)") "SnapPos:" // trim(EMIS_NAME(iem)), &
                  sum( snapemis (:,debug_li,debug_lj,:,iem) ) &
                 ,sum( snapemis_flat (debug_li,debug_lj,:,iem) )
        end if
    enddo !iem

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

       call CheckStop(err1, "De-Allocation error 1 - globland") 
       call CheckStop(err2, "De-Allocation error 2 - globland")
       call CheckStop(err3, "De-Allocation error 3 - globland")
       call CheckStop(err4, "De-Allocation error 4 - globland")
       call CheckStop(err5, "De-Allocation error 5 - globland")
       call CheckStop(err6, "De-Allocation error 6 - globland")

    end if

    ! now we have nrecmis and can allocate for gridrcemis:
    ! print *, "ALLOCATING GRIDRC", me, NRCEMIS
   allocate(gridrcemis(NRCEMIS,KEMISTOP:KMAX_MID,MAXLIMAX,MAXLJMAX),stat=err1)
   allocate(gridrcemis0(NRCEMIS,KEMISTOP:KMAX_MID,MAXLIMAX,MAXLJMAX),stat=err2)
   call CheckStop(err1, "Allocation error 1 - gridrcemis") 
   call CheckStop(err2, "Allocation error 2 - gridrcemis0")

  end subroutine Emissions

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine consistency_check()
  !------------------------------------------------------------------!
  !    checks that all the values given so far are consistent        !
  !------------------------------------------------------------------!

  character(len=30) :: errormsg
  
  errormsg = "ok"
  if ( size(EMIS_NAME) /= NEMIS_FILES    ) errormsg = " size EMISNAME wrong "

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
  !   And day-night factors are applied here:
  !       day_factor(11,0:1)                  ! 0=night, 1=day
  !
  !*************************************************************************

  implicit none
  type(date), intent(in) :: indate           ! Gives year..seconds
  integer :: i, j, k, f                   ! cooridnates, loop variables
  integer :: icc, ncc                        ! No. of countries in grid.

  integer :: ficc,fncc                       ! No. of countries with
  integer :: iqrc                            ! emis indices 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILES
  integer :: itot             ! index in xn()

  ! Save daytime value between calls, initialise to zero
  integer, save, dimension(NLAND) ::  daytime = 0  !  0=night, 1=day
  integer                         ::  hourloc      !  local hour 
  logical                         ::  hourchange   !             
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions

  real ::  ehlpcom,ehlpcom0(KEMISTOP:KMAX_MID)
  real ::  tfac, dtgrid    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s               ! source term (emis) before splitting
  integer :: iland, iland_timefac  ! country codes, and codes for timefac 

  real ::  sf              ! source term (emis) before splitting 
                           ! (for flat emissions)
  integer :: flat_iland    ! country codes (countries with flat emissions)

  integer, save :: oldday = -1, oldhour = -1

! If timezone=-100, calculate daytime based on longitude rather than timezone
  integer :: daytime_longitude, daytime_iland
 
! Initialize
    ehlpcom0(:)=0.0

   do k=KEMISTOP,KMAX_MID
      ehlpcom0(k) = GRAV* 0.001*AVOG/ (sigma_bnd(k+1) - sigma_bnd(k))
   enddo

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
              !==========================

               oldday = indate%day
           endif
     end if


    !..........................................
    !  Look for day-night changes, after local time correction
    !  (daytime(iland) not used if  LONGITUDE_TIME=true)

    do iland = 1, NLAND

       daytime(iland) = 0
       hourloc        = indate%hour + Country(iland)%timezone

       if ( hourloc  >=   7 .and.  hourloc <= 18 ) daytime(iland)=1

    end do ! iland


    if ( hourchange ) then 

         totemadd(:)  = 0.
         gridrcemis0(:,:,:,:) = 0.0 
         SumSnapEmis(:,:,:) = 0.0


        !..........................................
        ! Process each grid:

        do j = lj0,lj1
            do i = li0,li1

               ncc = nlandcode(i,j)            ! No. of countries in grid

                ! find the approximate local time:
                  hourloc= mod(nint(indate%hour+24*(1+glon(i,j)/360.0)),24)
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
                else
                   daytime_iland=daytime(iland)
                endif

                 !  As each emission sector has a different diurnal profile
                 !  and possibly speciation, we loop over each sector, adding
                 !  the found emission rates to gridrcemis as we go.
                 !  ==================================================


                do isec = 1, NSECTORS       ! Loop over snap codes

                   ! Calculate emission rates from snapemis, time-factors, 
                   ! and if appropriate any speciation fraction (NEMIS_FRAC)

                   iqrc = 0   ! index over emisfrac

                   do iem = 1, NEMIS_FILES 

                      tfac = timefac(iland_timefac,isec,iem) * &
                                 day_factor(isec,daytime_iland)

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

                   do k=KEMISTOP,KMAX_MID
                      do iqrc =1, nrcemis
                         gridrcemis0(iqrc,k,i,j) =   &
                            gridrcemis0(iqrc,k,i,j) + tmpemis(iqrc)*   &
                            ehlpcom0(k)*VERTFAC(KMAX_BND-k,isec) &
                            * emis_masscorr(iqrc)
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

            do iem = 1, NEMIS_FILES 

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
                    ehlpcom0(KMAX_MID) * emis_masscorr(iqrc)
             end do ! iem

!      ==================================================
       end do !ficc 
   end do ! i
 end do ! j
        if ( DEBUG .and.  debug_proc ) then    ! emis sum kg/m2/s
          call datewrite("SnapSum, kg/m2/s:" // trim(EMIS_NAME(iemCO)), &
                  (/ SumSnapEmis(debug_li,debug_lj,iemCO)  /) )
        end if

    call Set_Volc !set hourly volcano emission(rcemis_volc0)

  end if ! hourchange 


  ! We now scale gridrcemis to get emissions in molecules/cm3/s

   do k= KEMISTOP, KMAX_MID
     do j = lj0,lj1
        do i = li0,li1

              ehlpcom= roa(i,j,k,1)/(ps(i,j,1)-PT)
              do iqrc =1, NRCEMIS
                 gridrcemis(iqrc,k,i,j) =  &
                    gridrcemis0(iqrc,k,i,j)* ehlpcom
              enddo  ! iqrc
        end do ! i
      end do ! j
   end do ! k

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
    character*20 fname
    real ktonne_to_kgm2s, tonnemonth_to_kgm2s  ! Units conversion
    integer :: IQSO2                   ! Index of sox in  EMIS_NAME
    integer errcode
    real,    allocatable, dimension(:,:,:,:)  :: globemis 
    integer :: month,iem,ic,iic,isec, err3,icc
    real :: duml,dumh,tmpsec(NSECTORS),conv
    logical , save :: first_call=.true.
    real, dimension(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES) &
            ::  snapemis_month ! monthly emissions tonne/month

   ! For now, only the global runs use the Monthly files

        logical, parameter :: MONTHLY_GRIDEMIS= IS_GLOBAL      
        integer :: kstart,kend,nstart,Nyears
        real :: buffer(MAXLIMAX,MAXLJMAX),SumSoilNOx,SumSoilNOx_buff

if( USE_AIRCRAFT_EMIS )then
airn = 0.0 !ssp8W
!AIRCRAFT
kstart=KCHEMTOP
kend=KMAX_MID
do k=KEMISTOP,KMAX_MID
do j=1,ljmax
do i=1,limax

airn(k,i,j)=0.0

enddo
enddo
enddo
call ReadField_CDF('AircraftEmis_FL.nc','NOx',airn,nstart=current_date%month,kstart=kstart,kend=kend,interpol='mass_conservative', &
     needed=.true.,debug_flag=.true.)

! convert from kg(NO2)/month into molecules/cm3/s
! from kg to molecules: 0.001*AVOG/species(NO2)%molwt
! use roa to find dz for consistency with other emissions 
! (otherwise could have used z_bnd directly)
! dz=dP/(roa*GRAV)  dP=dA(k) + dB(k)*ps(i,j,1)
! dV=dz*dx*dy=dz*gridwidth**2/xm**2 *1e6 (1e6 for m3->cm3)
! from month to seconds: ndaysmonth*24*3600

conv=0.001*AVOG/species(NO2)%molwt*GRAV/gridwidth_m**2*1.0e-6/(nmdays(current_date%month)*24*3600)
!do k=KEMISTOP,KMAX_MID
do k=KCHEMTOP,KMAX_MID  !ssp8X
do j=1,ljmax
do i=1,limax

airn(k,i,j)=airn(k,i,j)*conv*(roa(i,j,k,1))/(dA(k) + dB(k)*ps(i,j,1))*xm2(i,j)

enddo
enddo
enddo

endif

!Soil NOx emissions
if(USE_SOIL_NOX)then 


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
   interpol='conservative',known_projection="lon lat",needed=.true.,debug_flag=.true.)
  if ( DEBUG_SOILNO.and.debug_proc ) write(*,*) "PROPER YEAR of SOILNO ", current_date%year, nstart
else
   !the year is not defined; average over all years
   Nyears=10 !10 years defined
   do iyr=1,Nyears 
      nstart=12*(iyr-1) + current_date%month  
      call ReadField_CDF('nox_emission_1996-2005.nc','NOX_EMISSION',buffer,nstart=nstart,&
      interpol='conservative',known_projection="lon lat",needed=.true.,debug_flag=.true.,UnDef=0.0)
      do j=1,ljmax 
         do i=1,limax
           SoilNOx(i,j)=SoilNOx(i,j)+buffer(i,j)
         end do
      end do
      if ( DEBUG_SOILNO.and.debug_proc ) then
         write(*,"(a,2i6,es10.3,a,2es10.3)") "Averaging SOILNO  inputs", &
           1995+(i-1), nstart,SoilNOx(debug_li, debug_lj), &
            "max: ", maxval(buffer), maxval(SoilNOx)
      !else if ( DEBUG_SOILNO ) then
      !     write(*,"(a,2i6,a,es10.3)") &
      !  "Averaging SOILNO  inputs", 1995+(i-1), nstart, "max: ", maxval(SoilNOx)
      end if
   enddo
   SoilNOx=SoilNOx/Nyears
endif
if ( DEBUG_SOILNO .and. debug_proc ) then
   write(*,"(a,i3,2es10.3)") "After SOILNO ", me, maxval(SoilNOx), SoilNOx(debug_li, debug_lj)
!else if ( DEBUG_SOILNO ) then
!   write(*,"(a,i3,es10.3, 2f8.2)") "After SOILNO ", me, maxval(SoilNOx), gb(1,1), gl(1,1)
end if

!for testing, compute total soil NOx emissions within domain
!convert from g/m2/day into kg/day
SumSoilNOx_buff=0.0
SumSoilNOx=0.0
do j=1,ljmax
   do i=1,limax      
      SumSoilNOx_buff=SumSoilNOx_buff+0.001*SoilNOx(i,j)*gridwidth_m**2*xmd(i,j)      
   enddo
enddo
CALL MPI_ALLREDUCE(SumSoilNOx_buff, SumSoilNOx , 1, &
     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
if(MasterProc)write(*,*)'Soil NOx emissions this month within domain',SumSoilNOx,' kg per day'

! convert from g(N)/m2/day into molecules/cm3/s
! from g to molecules: AVOG/14  14=molweight N
! use roa to find dz for consistency with other emissions 
! (otherwise could have used z_bnd directly)
! dz=dP/(roa*GRAV)  dP=dA(k) + dB(k)*ps(i,j,1)
! dV=dz*1e6 (1e6 for m3->cm3)
! from month to seconds: ndaysmonth*24*3600

conv=AVOG/14.0*GRAV*1.0e-6/(24*3600)
k=KMAX_MID!surface
do j=1,ljmax
   do i=1,limax      
      SoilNOx(i,j)=SoilNOx(i,j)*conv*(roa(i,j,k,1))/(dA(k) + dB(k)*ps(i,j,1))     
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
          do i = 1, NEMIS_FILES
            if ( trim( EMIS_NAME(i) ) == "sox" ) IQSO2 = i
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

              call ReadField(IO_DMS,fname,rdemis)

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

        end if  ! IQSO2>0


   if(MONTHLY_GRIDEMIS)then

   !Read monthly emission files

   if(first_call)then
      do j=lj0,lj1
         do i=li0,li1
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
    do iem = 1, NEMIS_FILES
!      if (trim(EMIS_NAME(iem)).ne.'nox' .and. trim(EMIS_NAME(iem)).ne.'co'.and.&
!           trim(EMIS_NAME(iem)).ne.'pm25'.and.&
!           trim(EMIS_NAME(iem)).ne.'voc'.and.trim(EMIS_NAME(iem)).ne.'nh3'.and.trim(EMIS_NAME(iem)).ne.'sox')cycle !
!      snapemis (:,:,:,:,iem) = 0.0 !NB all previous (snap)emis set to 0
      if ( MasterProc ) then

         globemis = 0.0

         write(fname,fmt='(''grid'',A,i2.2)')     &
              trim(EMIS_NAME(iem))//'.',month
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
   do iem = 1, NEMIS_FILES
!         write(*,*)'iem=',iem
!      if (trim(EMIS_NAME(iem)).ne.'nox' .and. trim(EMIS_NAME(iem)).ne.'co'.and.&
!           trim(EMIS_NAME(iem)).ne.'pm25'.and.&
!           trim(EMIS_NAME(iem)).ne.'voc'.and.trim(EMIS_NAME(iem)).ne.'nh3'.and.trim(EMIS_NAME(iem)).ne.'sox')cycle !
      conv = tonnemonth_to_kgm2s
      do j=lj0,lj1
         do i=li0,li1
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

end module Emissions_ml
