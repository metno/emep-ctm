! <Timefactors_mod.f90 - A component of the EMEP MSC-W Eulerian
!          Chemical transport Model>
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                    module Timefactors_mod

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!.......................................................................
!  DESCRIPTION:
!  Calculates emission temporal variation.
!
! 2024 updates
!  Simplified system to use timeFacs structure:
!   Monthly, Daily, Hourly... 
!  and where Monthly can be "CAMS_TEMPO_CLIM" or "GRIDDED"
! 2023 updates
!  New options
!  TimeFacBasis can be:
!   CAMS_TEMPO - can read day-of-year variations
!   CAMS_TEMPO_CLIM - reads climatological monthly, day-of-week, and hour-of-day files
!    from CAMS-TEMPO system.
!   DAY_OF_YEAR
!  IMPORTANT CAMS_TEMPO_CLIM needs ADD_SECTORS lines in config_emep.nml to
!  reset time-factor indices to 1-19, instead of SNAP's 1-11.
!  *** Check if CAMS_TEMPO day-of-year method should have done this.
!
! Older:
!  Reads monthly and daily (GENEMIS) factors for all emissions from files
!  Monthly.sox, Daily.sox, Monthly.nox, etc., -> in fac_emm, fac_edd arrays 
!  For every day, calculates emission factor "timefac" per country, emission 
!  sector, emission component
!  
!  Sets the day/night emissions variation in day_factor
!
!  Gridded Monthly factors:
!  If gridded monthly timefactors are choosen, the fac_emm only contains a 
!  normalisation factor. This normalization factor ensures that (before 
!  the proper Gridded Monthly factors correction is applied) each month 
!  contributes equally, i.e. the fac_emm will 
!  correct for differences causes by different number of week-ends, but not 
!  for differences in number of days in the month
!_____________________________________________________________________________

  use CheckStop_mod, only: CheckStop, StopAll
  use ChemDims_mod,  only: NEMIS_File 
  use Country_mod,   only: NLAND,Country
  use EmisDef_mod,   only: EMIS_FILE, TFAC_IDX_DOM, TFAC_IDX_TRAF, TFAC_IDX_AGR,&
                           TFAC_IDX_AGRK, TFAC_IDX_AGRL,&
                           IS_POW,IS_AGR,IS_AGRK,IS_AGRL,IS_TRAF,IS_DOM,IS_IND,&
                           TFAC_IDX_POW, SECTORS, NSECTORS, N_TFAC
  use GridValues_mod, only : i_fdom,j_fdom, debug_proc,debug_li,debug_lj
  use InterpolationRoutines_mod, only : Averageconserved_interpolate
  use Met_mod,       only: Getmeteofield
  use Config_module, only: MasterProc
  use Config_module, only: timeFacs !DS Feb 2024
  use Config_module, only: IIFULLDOM, JJFULLDOM
  use Config_module, only: iyr_trend ,USES  ! for GRIDDED_EMIS_MONTHLY_FACTOR 
  use Config_module, only: INERIS_SNAP1, INERIS_SNAP2, DegreeDayFactorsFile,&
                            GriddedMonthlyFacFile, DailyFacFile,MonthlyFacFile,&
                            DayofYearFacFile,&
                            monthly_timezoneFile, &
                            HourlyFacFile,HourlyFacSpecialsFile,&
                            USES
  use Debug_module,  only:   DEBUG ! => DEBUG%EMISTIMEFACS
  use Functions_mod, only : monthly_convolve ! DS added Feb 2024
  use NetCDF_mod,    only: GetCDF , ReadField_CDF
  use OwnDataTypes_mod, only: TXTLEN_FILE,TXTLEN_NAME
  use Par_mod,       only: MAXLIMAX,MAXLJMAX, limax,ljmax, me, li0, lj0, li1, lj1
  use Par_mod,       only: IRUNBEG, JRUNBEG, MSG_READ8
  use PhysicalConstants_mod, only: PI
  use SmallUtils_mod, only: find_index, key2str, basename, basedir, trims
  use Io_mod,        only:            &
                     open_file,       & ! subroutine
                     check_file,       & ! subroutine
                     ios,  IO_TIMEFACS  ! i/o error number, i/o label
  use Io_RunLog_mod, only: PrintLog
  use TimeDate_mod,  only:            &  ! subroutine, sets:
                     date,           &  ! date-type definition
                     nmdays, nydays, &  ! days per month (12), days per year
                     day_of_week,&      ! weekday
                     day_of_year        ! day count in year

  implicit none
  private

  !-- subroutines:

  public :: NewDayFactors
  public :: timefactors
  public :: DegreeDayFactors
  public :: Read_monthly_emis_grid_fac
  public :: Read_monthly_timezones
  public :: yearly_normalize

  !-- time factor stuff: 

  real, public, save, allocatable,&
     dimension(:,:,:) :: timefac ! overall emission 
                                 ! timefactor, calculated daily

  type, private :: timezone_t
    integer, allocatable, dimension(:,:) :: &
     map, Jan, inc   ! mapped tzs, January values, and increments
  end type
  type(timezone_t), public, save :: timezones

  real, public, save, allocatable, dimension(:,:) :: &
     timezone_map, timezone_Jan, timezone_inc

  real, public, save, allocatable, &
     dimension(:,:,:,:) :: fac_emm  ! Monthly factors
  logical, public, save :: InterpolateMonthEmis = .true. ! Whether interpolate (in time) Monthly factors or not

 ! Hourly for each day ! From EURODELTA/INERIS
 ! dims: NEMIS_FILE,N_TFAC,24,7,NLAND were N_TFAC currently 11
  real, public, save, allocatable,  &     
     dimension(:,:,:,:,:) :: fac_ehh24x7 ! Hour factors for 7 days (iemis,insec,ihh,idd2,indexCC)

 ! We keep track of min value for degree-day work
 !
  real, public, save, allocatable, &
     dimension(:,:,:) :: fac_min ! Min of Monthly factors
 !
  real, public, save, &
     dimension(12) :: fac_cemm  ! Change in monthly factors over the years

  real, public, save,allocatable,  &
     dimension(:,:,:,:) :: fac_edd  ! Daily factors

  ! normalized to one, replaces monthly and Day of week. 
  real, public, save,allocatable,  &
     dimension(:,:,:,:) :: fac_dayofyear  ! Daily factors over one year

  ! Heating-degree day factor for SNAP-2. Independent of country:
  logical, public, save :: Gridded_SNAP2_Factors = .false.
  real, public, allocatable,dimension (:,:), save :: gridfac_HDD
  !real, private, dimension (LIMAX,LJMAX), save :: tmpt2

  ! Used for general file calls and mpi routines below

  character(len=TXTLEN_FILE), private :: fname2   ! input filename - do not change 

  real, allocatable, public, save,  dimension(:,:,:,:):: GridTfac

  character(len=100) :: errmsg
  logical, public, save :: dbgTF = .false. ! shorthand for EMISTIMEFACS
  ! EMEP country nums doesn't always match grid-based country index. Be
  ! explicit
  integer, parameter, private :: &
      dbgICCemep  = 17   & ! Use AT for now. 
     ,dbgIsec    = 1   & ! Use GNFR_A for now. 
     ,dbgIemis   = 1     ! 'sox' usually
  integer, save, private  :: &
      dbgICindex   &! not always same as dbgICCemep
     ,dbgIsecRT    ! Use GNFR_F for now. 

contains


  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine timefactors(year)

   !.......................................................................
   !  DESCRIPTION:
   !  Read in monthly and daily factors, -> fac_emm, fac_edd arrays
   !  The input files are Monthly.sox, Daily.sox, Monthly.nox, etc.
   !  Sets the day/night variation in day_factor
   !
   !  D. Simpson, Peter Wind, Alvaro Valdebenito, 1999-2024
   !.......................................................................

  !-- Input
  integer, intent(in) :: year

  !-- Outputs -  module's fac_emm, fac_edd, day_factor, etc.

  !-- local
  integer ::  emepICC, insec    ! EMEP country number and sector value read from femis etc
  integer ::  indexCC           ! Index of country code in Country array. Can differ from emep number
  integer ::  i, isec, n
  integer ::  idd, idd2, ihh, mm ! Loop and count variables
  integer ::  iemis             ! emission count variables

  !integer :: weekday            ! 1=monday, 2=tuesday etc.
  real    :: xday, sumfac,vmin,vmax, vmean ! used in interpolation, testing
  real    :: tmp24(24)          ! used for hourly factors
  character(len=2000) :: inputline ! NB: 24 real number can be many hundred characters long
  real :: fracchange
  real :: buff(366)
  logical :: found_HourlyFacFile, found
  character(len=*), parameter:: dtxt='tfacs:'
  character(len=TXTLEN_NAME):: secname
  character(len=10) :: code
  integer :: maxidx = 0
  logical :: dbgccsec = .false.
  character(len=80) :: dbgmsg
  
  if (DEBUG%EMISTIMEFACS .and. MasterProc )  then
    write(unit=6,fmt=*) dtxt//" N_TFAC:", N_TFAC, "timeFacs%Monthly:", timeFacs%Monthly
    dbgTF = .true.
  end if

   call CheckStop( nydays < 365, &
      dtxt//"Timefactors: ERR:Call set_nmdays before timefactors?")

!  https://github.com/metno/emep-mscw/issues/179#issuecomment-733757660
!  call CheckStop(  N_TFAC /= 11 , &
!     "Timefactors: ERR:Day-Night dimension wrong!")


   fac_emm(:,:,:,:) = 1.0
   fac_min(:,:,:) = 1.0
   fac_edd(:,:,:,:) = 1.0
   fac_ehh24x7 = 1.0

   if(N_TFAC==0)return
   
   !print *, dtxt//' PREG', me, .not. USES%GRIDDED_EMIS_MONTHLY_FACTOR
   !print *, dtxt//' PRET', me, .not. TimeFacBasis == 'DAY_OF_YEAR'
   !print *, dtxt//' PREX', me, trim(timeFacs%Monthly)

   if( timeFacs%Monthly == 'GRIDDED' ) then
      ! will read monthyl gridded later
      if (MasterProc) write(*,*) dtxt//' GRIDDED monthly'
! #################################
   else if ( timeFacs%Day_of_Year )then      
       !read directly day of the year timefactors

       fac_dayofyear = 1.0
       do iemis = 1, NEMIS_FILE
          fname2 = key2str(DayofYearFacFile,'POLL',trim ( EMIS_FILE(iemis) ))
          if(me==0)write(*,*)'reading '//trim(fname2)
          call open_file(IO_TIMEFACS,"r",fname2,needed=.true.,skip=1)
          
          call CheckStop( ios, dtxt//" Opening error in DayofYearfac")
          
          n = 0
          do
             read(IO_TIMEFACS,fmt=*,iostat=ios) code,secname, &
                  (buff(i),i=1,nydays)         
             if ( ios <  0 ) print *, "IOS eof IEMLOOP"//trim(basename(fname2)), ios
             if ( ios <  0 ) exit   ! End of file
             indexCC=find_index(code,Country(:)%code)
             if(indexCC<1.or.indexCC>NLAND)then
                if(me==0.and.insec==1.and.iemis==1)write(*,*)dtxt//"Day of year fac code not used: "//trim(code)
                cycle
             end if
             insec=find_index(secname,SECTORS(:)%longname)             
             call CheckStop(insec>NSECTORS .or. insec==0, secname//' not defined')
             maxidx = max(insec,maxidx)

             fac_dayofyear(insec,indexCC,iemis,1:nydays)=buff(1:nydays)
             
             n = n + 1
             
             !-- Sum over days of the year
             xday =  sum( fac_dayofyear(insec,indexCC,iemis,1:nydays)) / nydays
             
             if (( xday > 1.01 .or. xday < 0.99) .and. MasterProc) then
                write(*,*)insec,trim(code)//" Warning: Day of year- not normalised. Renormalizing with factor ",xday
             end if
             !We renormalize for ensuring exact norm (so that not too many digits need to be written in file)
             fac_dayofyear(insec,indexCC,iemis,1:nydays) = fac_dayofyear(insec,indexCC,iemis,1:nydays)/xday
             
          end do
          
          close(IO_TIMEFACS)
          if (dbgTF) write(*,fmt=*) dtxt//"Read ", n, &
                  " DOY records from ", trim(fname2)
          
       end do  ! NEMIS_FILE
   else ! not GRIDDED or Day_of_Year

!  #################################
!  1) Read in Monthly factors, and determine min value (for baseload)

! Note: the inverse of fac_emm/fac_cemm is applied after reading the monthly 
!       sector emissions, in order to cancel subsequent application of fac_emm, 
!       but keep the fac_cemm variations

  ! Summer/winter SNAP1 ratios reduced from 1990 to 2010 when using "older"
  ! MonthlyFacs (e.g. GENEMIS were from 1994-era).
  ! Follows data presented in Grennfelt & Hov, Ambio, 2005, see Simpson et al 2012

   fac_cemm(:) = 1.0

   if(MasterProc) then
     write(*,*) dtxt//"MonthlyFacBasis:"//trim(timeFacs%Monthly)
     write(*,*) dtxt//"DailyFacBasis:  "//trim(timeFacs%Daily)
     write(*,*) dtxt//"HourlyFacBasis:  "//trim(timeFacs%Hourly)
     write(*,*) dtxt//"DAYOFYEARTIMEFAC:", timeFacs%Day_of_Year
   end if

   select case(timeFacs%Monthly)  ! MonthlyFacBasis)
   case("ECLIPSE")
      call CheckStop(index(MonthlyFacFile,'may2021')<1, & !CRUDE and TMP!
         dtxt//'Incostistent:'//trim(timeFacs%Monthly)//' vs '//MonthlyFacFile)

   case("CAMS_TEMPO", "CAMS_TEMPO_CLIM")
      call CheckStop(index(MonthlyFacFile,'cams_tempo')<1, &
         dtxt//'Incostistent:'//trim(timeFacs%Monthly)//' vs '//MonthlyFacFile)
         
   case("GENEMIS") ! check eclipse not in FacFile
     call CheckStop(index(MonthlyFacFile,'xJun2012')<1, & !CRUDE and TMP!
         dtxt//'Incostistent:'//trim(timeFacs%Monthly)//' vs '//MonthlyFacFile)
     
     fracchange=0.005*(iyr_trend -1990)
     fracchange=max(0.0,fracchange) !do not change before 1990
     fracchange=min(0.1,fracchange) !stop change after 2010 
     !equal 1.1/0.9=1.22 summer/winter change
     do i = 1, NSECTORS
      if(IS_POW(i)) then
         write(unit=6,fmt=*) dtxt//"Change summer/winter ratio in "//&
              trim(SECTORS(i)%longname)//" by ", fracchange
      end if
     end do
     do mm=1,12
      !Assume max change for august and february
      fac_cemm(mm)  = 1.0 + fracchange * cos ( 2 * PI * (mm - 8)/ 12.0 )
      write(unit=6,fmt="(a,i3,f8.3,a,f8.3)") dtxt//"S-W change in fac_cemm ", mm,fac_cemm(mm)
     end do
   case default
      call StopAll(dtxt//'ERROR MonthlyFac:'//trim(timeFacs%Monthly)//' vs '//MonthlyFacFile)
   end select ! timeFacs%Monthly
   write(*,"(a,f8.4)") dtxt//"Mean fac_cemm ", sum( fac_cemm(:) )/12.0
   
   if( INERIS_SNAP1 ) fac_cemm(:) = 1.0   ! Hardly ever used.

   do iemis = 1, NEMIS_FILE

       fname2 = key2str(MonthlyFacFile,'POLL',trim ( EMIS_FILE(iemis) ))

       call open_file(IO_TIMEFACS,"r",fname2,needed=.true.)

       n = 0
       !print *, dtxt//'MFACC', iemis
       do 
         select case(timeFacs%Monthly)
         case("ECLIPSE","GENEMIS")
            read(IO_TIMEFACS,fmt=*,iostat=ios) emepICC,insec,(buff(mm),mm=1,12)
            if( ios <  0 ) exit     ! End of file
            call CheckStop( ios, dtxt//": Read error in Monthlyfac")

            indexCC=find_index(emepICC,Country(:)%icode)
         case("CAMS_TEMPO","CAMS_TEMPO_CLIM")
            !print *, dtxt//'MFACD', iemis
            read(IO_TIMEFACS,fmt=*,iostat=ios) code,secname,(buff(mm),mm=1,12)
             if ( ios <  0 ) print *, "IOS eof CAMS"//trim(basename(fname2)), ios
            !print *, dtxt//'MFACE:'//trim(code)//trim(secname)
            if( ios <  0 ) exit     ! End of file
            call CheckStop( ios, dtxt//": Read error in Monthlyfac")

            indexCC=find_index(code,Country(:)%code)
            emepICC = Country(indexCC)%icode ! just for print out
            if(emepICC==dbgICCemep) dbgICindex=indexCC
            insec=find_index(secname,SECTORS(:)%longname)             
            dbgccsec = dbgTF .and. iemis==1 .and. insec==dbgIsec .and. emepICC==dbgICCemep
            if(dbgccsec) write(*,*)dtxt//'CHECKIC', indexCC, emepICC, code, &
                trim(Country(indexCC)%code) ! e.g. for IT indexCC=16, emepICC=15
         end select
            
         maxidx = max(insec,maxidx)
         if(insec>N_TFAC) cycle
         if(indexCC<1.or.indexCC>NLAND)then
            if(MasterProc.and.insec==1.and.iemis==1)&
               write(*,*)dtxt//" Monthlyfac code not used ",emepICC
            cycle
         end if
         call CheckStop(insec>NSECTORS .or. insec<=0, trim(secname)//' not defined')

         ! Feb 2024: added possibulity to smooth MonthlyFacs:
         if (timeFacs%MonthlySmoothFac < 100) buff(1:12) = &
             monthly_convolve(buff(1:12),timeFacs%MonthlySmoothFac)

         fac_emm(indexCC,1:12,insec,iemis)=buff(1:12)

         ! Temporary and crude implementation for BIDIR tests:
         if ( EMIS_FILE(iemis) == 'nh3' .and. timeFacs%MonthlyNH3 == 'LOTOS' ) then
            fac_emm(indexCC,1:12,insec,iemis) = &
               [ 0.60, 0.66,1.50,1.36,1.02,1.17,1.19,1.27,0.93,0.89,0.77,0.64]
         end if
         !defined after renormalization and send to al processors:
         ! fac_min(emepICC,insec,iemis) = minval( fac_emm(emepICC,:,insec,iemis) ) ???? OR indexCC??

         !if( dbgTF.and.insec==TFAC_IDX_DOM.and.iemis==1  ) &
         if( dbgccsec ) then
            dbgmsg=trims(dtxt//"emm:"// EMIS_FILE(iemis)//":"//&
                    Country(indexCC)%code//":"//SECTORS(insec)%longname)
            write(*,"(a,3i3,f7.3,a,12f6.2)")  dbgmsg, &
              emepICC,insec,iemis, fac_min(indexCC,insec,iemis),&
               " : ",  ( fac_emm(indexCC,mm,insec,iemis), mm=1,12)
         end if

         n = n + 1
       end do

       close(IO_TIMEFACS)

      ! Apply change in monthly factors for PUBLIC POWER (SNAP 1)
       do indexCC = 1, NLAND
          do i = 1, N_TFAC !must loop only once over each timefac index
             found = .false. 
             do isec=1,NSECTORS
                if(SECTORS(isec)%timefac == i .and. IS_POW(isec)) found = .true.
             end do
             if (.not. found) cycle
             sumfac=0.0
             do mm=1,12
                fac_emm(indexCC,mm,i,iemis)=fac_emm(indexCC,mm,i,iemis)*fac_cemm(mm)
                sumfac=sumfac+fac_emm(indexCC,mm,i,iemis)
             end do
             ! normalize
             do mm=1,12
                fac_emm(indexCC,mm,i,iemis)=fac_emm(indexCC,mm,i,iemis)*12./sumfac
             end do
          end do
       end do

       if (dbgTF) then
         if (iemis==1) write(*,fmt='(2a)') dtxt//&
                  "Reading MONTH records from ", trim(basedir(fname2) ) 
         write(*,fmt='(a,i5,2a)') dtxt//"Read ", n, &
               " MONTH records from ", trim(basename(fname2) ) 
       end if
               !" MONTH records from ", trim(fname2) 
   end do  ! iemis

   end if ! .not. GRIDDED MONTHLY/'DAY_OF_YEAR'


   if ( .not. timeFacs%Day_of_Year ) then ! Uses Monthly, Daily facs

     do iemis = 1, NEMIS_FILE
       fname2 = key2str(DailyFacFile,'POLL',trim ( EMIS_FILE(iemis) ))
       call open_file(IO_TIMEFACS,"r",fname2,needed=.true.)

       call CheckStop( ios, dtxt//" Opening error in Dailyfac")

       n = 0
       DoWLOOP: do
         if ( timeFacs%Daily == "CAMS_TEMPO_CLIM") then ! txt-based indices
           read(IO_TIMEFACS,fmt=*,iostat=ios) code,secname, (buff(i), i=1,7)
           indexCC=find_index(code,Country(:)%code)
           insec=find_index(secname,SECTORS(:)%longname)             
           emepICC = Country(indexCC)%icode ! just for print out
         else
           read(IO_TIMEFACS,fmt=*,iostat=ios) emepICC, insec, &
              (buff(i),i=1,7)         
           indexCC=find_index(emepICC,Country(:)%icode)
         end if
         if ( ios <  0 ) print *, "IOS eof DoW2"//trim(basename(fname2)), ios
         if ( ios <  0 ) exit   ! End of file

         maxidx = max(insec,maxidx)
         if(indexCC<1.or.indexCC>NLAND)then
            if(me==0.and.insec==1.and.iemis==1)write(*,*)&
                    dtxt//"Dailyfac code not used",emepICC
            cycle
         end if
         if(insec>N_TFAC) cycle
         fac_edd(indexCC,1:7,insec,iemis)=buff(1:7)
         call CheckStop( ios, dtxt//" Read error in Dailyfac")

         n = n + 1

         !-- Sum over days 1-7
         xday =  sum( fac_edd(indexCC,1:7,insec,iemis) ) / 7.0

         call CheckStop( xday > 1.001 .or. xday < 0.999, &
                dtxt//" ERROR: Dailyfac - not normalised")

         !if( dbgTF .and. emepICC == dbgICCemep .and. iemis == 1  ) then
         !   write(*,"(a,3i3,f7.3,a,12f6.2)") dtxt//"edd tfac:"// &
         !    trim(EMIS_FILE(iemis))//":"//trim(Country(indexCC)%code)//": "!, &
         !    !emepICC,insec,iemis, " : ",  ( fac_edd(indexCC,idd,insec,iemis), idd=1,7)
       end do DoWLOOP

       close(IO_TIMEFACS)
       if (dbgTF) then
         if (iemis==1) write(*,fmt='(2a)') dtxt//&
                  "Reading WEEK records from ", trim(basedir(fname2) ) 
         write(*,fmt='(a,i5,2a)') dtxt//"Read ", n, &
               " WEEK records from ", trim(basename(fname2) ) 
       end if

     end do  ! NEMIS_FILE

   endif ! testing DAY_OF_YEAR

!  #################################
!  3) Read in hourly (24x7) factors, options set in run script.
   ! INERIS option has 11x24x7 emissions factor
   ! TNO2005 option has 11x24 
   ! EMEP2003 option has very simple day night
!
  
   CLIMTEMPO: if ( timeFacs%Hourly == "CAMS_TEMPO_CLIM") then

   ! ===== START New June 2023 ===========================
   ! Will read country dependent hourly facs

     HourlyFacSpecialsFile = HourlyFacFile  ! Will handle in Specials below.
     found_HourlyFacFile=.false.
     if ( dbgTF ) write(*,*) dtxt//' SKIPS 1st HourlyFac'

   ! ===== END New June 2023 ===========================

   else  ! Reads older HourlyFacs, no country dependence

     fname2 = trim(HourlyFacFile) ! From EURODELTA/INERIS/TNO or EMEP2003
     write(unit=6,fmt=*) dtxt//"Starting HOURLY-FACS from "//trim(fname2)
     call open_file(IO_TIMEFACS,"r",fname2,needed=.false.,iostat=ios)
     found_HourlyFacFile=(ios==0)
     FOUNDHRLY: if(found_HourlyFacFile)then
       n = 0
       HRLYLOOP: do 
         read(IO_TIMEFACS,"(a)",iostat=ios) inputline
         n = n + 1
         !if(dbgTF)write(*,*) "HourlyFacs ", n, trim(inputline)
         if ( ios <  0 ) print *, "IOS eof HRLLOOP"//trim(basename(fname2)), ios
         if ( ios <  0 ) exit     ! End of file
         if( index(inputline,"#")>0 ) then ! Headers
            if(n==1) call PrintLog(trim(inputline))
            cycle
         else
            read(inputline,fmt=*,iostat=ios) idd, insec, &
                 (tmp24(ihh),ihh=1,24)
            if( dbgTF ) write(*,"(a,2i3,24f5.2)") dtxt//"HOURLY=> ",&
                 idd, insec, tmp24(:) !(1), tmp24(13)
         end if
         maxidx = max(insec,maxidx)
         if(insec>N_TFAC) cycle
        
         if(  idd == 0 ) then ! same values every day
           do idd2 = 1, 7
              do ihh=1,24
                 do iemis = 1, NEMIS_FILE
                    fac_ehh24x7(iemis,insec,ihh,idd2,:) = tmp24(ihh)
                 end do
              end do
           end do
         else
           do ihh=1,24
              do iemis = 1, NEMIS_FILE
                 fac_ehh24x7(iemis,insec,ihh,idd,:) = tmp24(ihh)
              end do
           end do
         end if
           
         ! Use sumfac for mean, and normalise within each day/sector
         ! (Sector 10 had a sum of 1.00625)

         if(  idd > 0 ) then !day value defined
              sumfac = sum(fac_ehh24x7(1,insec,:,idd,1))/24.0
              !if(dbgTF .and. MasterProc) write(*,"(a,2i3,3f12.5)") &
              !     'HOURLY-FACS mean min max', idd, insec, sumfac, &
              !     minval(fac_ehh24x7(1,insec,:,idd,1)), &
              !     maxval(fac_ehh24x7(1,insec,:,idd,1))
              
              do iemis = 1, NEMIS_FILE    !fac_ehh24x7(NEMIS_FILE,N_TFAC,24,7,NLAND)
                 fac_ehh24x7(iemis,insec,:,idd,:) = fac_ehh24x7(iemis,insec,:,idd,:) / sumfac
              end do
              
         else
              !same value every day
              sumfac = sum(fac_ehh24x7(1,insec,:,1,1))/24.0
              do idd2 = 1, 7
                 do iemis = 1, NEMIS_FILE
                    fac_ehh24x7(iemis,insec,:,idd2,:) = fac_ehh24x7(iemis,insec,:,idd2,:) / sumfac
                 end do
              end do
         end if
       end do HRLYLOOP

       if (dbgTF) write(unit=6,fmt=*) dtxt//"Read ", n, &
              " EHH24x7 records from ", trim(fname2)
       call CheckStop ( any(fac_ehh24x7 < 0.0 ) , dtxt//"Unfilled efac_ehh24x7")

       close(IO_TIMEFACS)
     end if FOUNDHRLY  ! (found_HourlyFacFile)then
   end if  CLIMTEMPO ! CAMS_TEMPO_CLIM test 
   if ( dbgTF ) write(*,*) dtxt//' CLIMTEMPO default hourly? ', found_HourlyFacFile

!3.1)Additional country and species specific hourly time factors

    if(dbgTF) write(*,'(a)') dtxt//'Hourly SPECIALs? '//trim(HourlyFacSpecialsFile)
    SPECEMLOOP: do iemis = 1, NEMIS_FILE
       fname2 = key2str(HourlyFacSpecialsFile,'POLL',trim (EMIS_FILE(iemis)))
       call open_file(IO_TIMEFACS,"r",fname2,needed=.not.found_HourlyFacFile,iostat=ios)
       !if(dbgTF) write(*,'(a,i4,L2,i,a)') dtxt//'Hourly SPECIAL', iemis,&
       !    found_HourlyFacFile,ios, trim(fname2)
       if ( ios /= 0 ) then
          if(me==0 .and. fname2/= 'NOTSET') write(*,*)dtxt//&
            'Special hourly factors not found (but not needed): ',trim(basename(fname2))
          cycle SPECEMLOOP
       end if
       n = 0
       SPECLOOP: do 
          read(IO_TIMEFACS,"(a)",iostat=ios) inputline
          n = n + 1
          !if(dbgTF)write(*,"(a,i4,a)") dtxt//"HourlyFacsSpecials ", n, trim(inputline(1:50))
          if ( ios <  0 ) exit     ! End of file
          if( index(inputline,"#")>0 ) then ! Headers
             if(n==1) call PrintLog(trim(inputline))
             cycle
          else
            ! ===== START New June 2023 ===========================
            !CRUDE hard-code also for name cams_tempo for now:
            if ( timeFacs%Hourly  == "CAMS_TEMPO_CLIM".or. & 
               index( HourlyFacSpecialsFile,'cams_tempo') > 0 ) then
              read(inputline,fmt=*,iostat=ios) code,secname, idd, &
                 (tmp24(ihh),ihh=1,24)
              indexCC=find_index(code,Country(:)%code)
              emepICC = Country(indexCC)%icode ! just for print out
              insec=find_index(secname,SECTORS(:)%longname)             
            ! ===== END New June 2023 ===========================
             else ! e.g. INERIS, 
               read(inputline,fmt=*,iostat=ios) emepICC, idd, insec, &
                  (tmp24(ihh),ihh=1,24)
             end if ! JUNE 2023
             maxidx = max(insec,maxidx)
             dbgccsec = dbgTF .and. iemis==1 .and. insec==dbgIsec .and. emepICC==dbgICCemep

             if(insec>N_TFAC) then
               if ( dbgTF) write(*,*) dtxt//'Warning insec-high', insec, secname
               cycle
             end if
             indexCC =find_index(emepICC,Country(:)%icode)
             !if( dbgTF .and. emepICC==dbgICCemep  ) write(*,'(a,4i4,2f10.4,a)') &
             !  dtxt//"HOURLY SPECIAL=> ", indexCC , emepICC,idd, insec,&
             !   tmp24(1), tmp24(13), " "//trim(timeFacs%Hourly)
             
             if( indexCC <0 .and. emepICC/=0)then
                if(MasterProc) write(*,*)dtxt// &
                 "Warning: HourlyFacsSpecials, country code not recognized:",&
                   trim(inputline)
                cycle
             endif
          end if

          if(  idd == 0 ) then ! same values every day
             do idd2 = 1, 7
                do ihh=1,24
                   if(emepICC/=0)then
                      fac_ehh24x7(iemis,insec,ihh,idd2,indexCC) = tmp24(ihh)
                   else
                      fac_ehh24x7(iemis,insec,ihh,idd2,:) = tmp24(ihh)
                   end if
                end do
                if( dbgccsec ) write(*,'(a,2i4,2f10.4)') dtxt//'HH0:', idd2, insec, tmp24(1), tmp24(13)
             end do
             if(emepICC/=0)then
                sumfac = sum(fac_ehh24x7(iemis,insec,:,1,indexCC))/24.0
                do idd2 = 1, 7
                   fac_ehh24x7(iemis,insec,:,idd2,indexCC) = fac_ehh24x7(iemis,insec,:,idd2,indexCC) / sumfac
                end do
                if( dbgccsec ) write(*,'(a,2i4,2f10.4)') dtxt//'HH1:', insec, emepICC, sumfac, fac_ehh24x7(iemis,insec,dbgIsec,idd2,indexCC)
             else
                sumfac = sum(fac_ehh24x7(iemis,insec,:,1,1))/24.0
                do idd2 = 1, 7
                   fac_ehh24x7(iemis,insec,:,idd2,:) = fac_ehh24x7(iemis,insec,:,idd2,:) / sumfac
                end do
                if( dbgccsec ) write(*,'(a,2i4,2f10.4)') dtxt//'HH2:', insec, emepICC, sumfac, fac_ehh24x7(iemis,insec,dbgIsec,idd2,indexCC)
             end if
          else  
             do ihh=1,24
                if(emepICC/=0)then
                   fac_ehh24x7(iemis,insec,ihh,idd,indexCC) = tmp24(ihh)
                else
                   fac_ehh24x7(iemis,insec,ihh,idd,:) = tmp24(ihh)
                end if
             end do
             if(emepICC/=0)then
                sumfac = sum(fac_ehh24x7(iemis,insec,:,idd,indexCC))/24.0
                fac_ehh24x7(iemis,insec,:,idd,indexCC) = fac_ehh24x7(iemis,insec,:,idd,indexCC) / sumfac
             else
                sumfac = sum(fac_ehh24x7(iemis,insec,:,idd,1))/24.0
                fac_ehh24x7(iemis,insec,:,idd,:) = fac_ehh24x7(iemis,insec,:,idd,:) / sumfac
             end if
             if( dbgccsec ) write(*,'(a,3i4,2f10.4)') dtxt//'HH3:', insec, & !CAMS_TEMPO uses this
                     emepICC, idd, sumfac, fac_ehh24x7(iemis,insec,dbgIsec,idd,indexCC)
          end if
          
          ! Use sumfac for mean, and normalise within each day/sector
          if(emepICC==0)indexCC =1
          !XCAMEO if(dbgTF .and. MasterProc) write(*,"(a,3i3,3f12.5)") &
          ! Just Fri-Sun facs now for dbg settings
          if( dbgTF.and. emepICC == dbgICCemep.and.iemis==dbgIemis.and.idd>4 ) &
            write(*,"(a50,3i3,4f8.4)") &
            dtxt//'HOURLY-FACS mean hr1 min max,'// trim(EMIS_FILE(iemis))// &
             ':'// SECTORS(insec)%longname, &
             insec, idd, emepICC, sumfac, fac_ehh24x7(iemis,insec,1,idd,indexCC),&
             minval(fac_ehh24x7(iemis,insec,:,idd,indexCC)), &
             maxval(fac_ehh24x7(iemis,insec,:,idd,indexCC))
          
       end do SPECLOOP
       close(IO_TIMEFACS)
    end do SPECEMLOOP ! NEMIS_FILE


    if(.not.found_HourlyFacFile) call CheckStop( any(fac_ehh24x7 < 0.0 ) ,&
         dtxt//"Unfilled efac_ehh24x7")

   if(N_TFAC/=maxidx .and. MasterProc) write(*,*)'Warning: ',N_TFAC,&
    ' timefactor indices defined in SECTOR, but ',maxidx,' found in the timefactor files'

! #######################################################################
! 4) Normalise the monthly-daily factors. This is needed in order to
!    account for leap years (nydays=366) and for the fact that different
!    years have different numbers of e.g. Saturdays/Sundays. 
!    Here we execute the same interpolations which are later done
!    in "NewDayFactors", and scale efac_mm if necessary.

       write(unit=6,fmt="(a,I6,a,I5)")dtxt//" Time factors normalisation: ",&
                                        nydays,' days in ',year
       if ( timeFacs%Day_of_Year ) call yearly_normalize(year)

!#########################################################################
!

    if (dbgTF ) then 
       write(unit=6,fmt=*) dtxt//"Ending subroutine timefactors"
           !fac_emm(indexCC,1:12,insec,iemis)=buff(1:12)
           !fac_edd(indexCC,1:7,insec,iemis)=buff(1:7)
           !fac_ehh24x7(iemis,insec,ihh,idd2,:)
       indexCC=dbgICindex; insec=dbgIsec;  iemis=dbgIemis
       write(*,*) dtxt//" test of time factors, dbg: "//&
             trim(EMIS_FILE(iemis))//' EmepIC:'//trim(Country(dbgICindex)%code),&
             ':'//trim(SECTORS(insec)%longname)
       if ( timeFacs%Day_of_Year ) then
          do mm = 1, 12
          write(*, "(a,i2,i6,f8.3,3f8.4)") dtxt//'-doy', mm, nydays, sumfac,  &
               fac_dayofyear(insec,indexCC,iemis,15+((mm-1)*30)),&
                fac_emm(indexCC,mm,insec,iemis), fac_edd(indexCC,1,insec,iemis)
          end do ! mm
       else
          do mm = 1, 12
             write(*, "(a,i3,i6,f8.3,9f8.4)") dtxt//'-md', mm, nydays, sumfac,  &
              fac_emm(indexCC,mm,insec,iemis), &
              fac_edd(indexCC,1,insec,iemis), fac_edd(indexCC,7,insec,iemis), &
              sum(fac_ehh24x7(iemis,insec,:,1,indexCC))/24.0, &
              sum(fac_ehh24x7(iemis,insec,:,7,indexCC))/24.0
          end do ! mm
          mm=1
          do insec = 1, N_TFAC ! NSECTORS
            write(*,'(a,i3,f10.4,a)') dtxt//'-Jan:', insec, &
               fac_emm(indexCC,mm,insec,iemis), ' '//trim(SECTORS(insec)%longname)
            if(trim(SECTORS(insec)%longname)=='GNFR_F') dbgIsecRT=insec
          end do
       end if
       write(*,"(a,4f8.3)") dtxt//" day factors traffic 24x7", &
           fac_ehh24x7(iemis,dbgIsecRT,1,4,dbgICindex),&  ! 1am, Jan4?
           fac_ehh24x7(iemis,dbgIsecRT,13,4,dbgICindex), & ! 1pm, Jan4?
              minval(fac_ehh24x7), maxval(fac_ehh24x7)
       do indexCC = 1, NLAND
         emepICC = Country(indexCC)%icode ! just for print out
                    !fac_ehh24x7(iemis,insec,ihh,idd2,indexCC)
         vmin = minval(fac_ehh24x7(iemis,dbgIsecRT,:,:,indexCC))

         if (vmin>0.999) cycle ! Skip 1.0000 defaults

         vmax  = maxval(fac_ehh24x7(iemis,dbgIsecRT,:,:,indexCC))
         vmean = sum(fac_ehh24x7(iemis,dbgIsecRT,:,:,indexCC))/(7*24)
         write(*,'(a26,2i4,3f12.4)') trims(dtxt//'RTsrch:'//&
            Country(indexCC)%code//':'//EMIS_FILE(iemis)), &
          indexCC, emepICC, vmin, vmax, vmean
         if( vmean < 0.9999 ) then
            dbgmsg=trims(dtxt//"VMEAN-LOW:"// EMIS_FILE(iemis)//":"//&
                    Country(indexCC)%code//":"//SECTORS(insec)%longname)
           print *, 'ERROR VMEAN:', vmean, trim(dbgmsg)
           call StopAll(dbgmsg)
         end if
       end do
    end if ! dbgTF

 end subroutine timefactors
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
 subroutine Read_monthly_timezones(month)
   character(len=*), parameter :: dtxt='readTZ:'
   integer, intent(in) :: month
   logical, save :: first_call = .true.
   real, allocatable, dimension(:,:), save :: tmpwork
   integer :: i,j

    if(.not.allocated(timezones%map))then
       allocate(tmpwork(LIMAX,LJMAX))        ! for input, floats
       allocate(timezones%map(LIMAX,LJMAX))  ! from netcdf timezones
       allocate(timezones%Jan(LIMAX,LJMAX))  ! for january
       allocate(timezones%inc(LIMAX,LJMAX))  ! increments
    end if

    if (dbgTF) write(*,*) dtxt//"reading tzones"

    call ReadField_CDF(monthly_timezoneFile,'tz',tmpwork,month,&
       known_projection='lon lat',interpol='zero_order',needed=.true.,debug_flag=.false.)
    forall(i=1:LIMAX,j=1:LJMAX) &
       timezones%map(i,j) = nint(tmpwork(i,j))


    if ( first_call ) then
      if ( month==1) then
        timezones%Jan = timezones%map
        if(dbgTF) write(*,'(a,2i7)')'AJTIMEZONES',   maxval(timezones%Jan), minval(timezones%Jan)
      else                                 ! Need to get Jan values:
        call ReadField_CDF(monthly_timezoneFile,'tz',tmpwork,1,&
         known_projection='lon lat',interpol='zero_order',needed=.true.,debug_flag=.false.)
        if(dbgTF) write(*,'(a,2es10.3)')'FBJTIMEZONES',   maxval(tmpwork), minval(tmpwork)
        forall(i=1:LIMAX,j=1:LJMAX) &
            timezones%Jan(i,j) = nint(tmpwork(i,j))
        if(dbgTF) write(*,'(a,2i7)')'BJTIMEZONES',   maxval(timezones%Jan), minval(timezones%Jan)
      end if
      first_call = .false.
    end if
    timezones%inc = timezones%map - timezones%Jan  ! used for e.g.summertime

    if (dbgTF) then
      write(*,*) dtxt//"done reading tzones"
      write(*,'(a,2f7.1)')'FTIMEZONES', maxval(tmpwork), minval(tmpwork)
      write(*,'(a,2i7)')'ITIMEZONES',   maxval(timezones%map), minval(timezones%map)
      write(*,'(a,2i7)')'JTIMEZONES',   maxval(timezones%Jan), minval(timezones%Jan)
    end if

 end subroutine Read_monthly_timezones

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine NewDayFactors(newdate)
  
  ! Calculates the monthly and daily factors for emission temporal variation
  ! for each country, emission, and sector.  Called at midnight every day.
  !
  ! Uses arays: 
  !     fac_emm(NLAND,NM,N_TFAC,NEMIS_FILES)    ! Jan - Dec.
  !     fac_edd(NLAND,7,N_TFAC,NEMIS_FILES)     ! Monday=1, Sunday=7
  !
  ! Outputs:
  !    real timefac(NLAND,N_TFAC,NEMIS_FILES)
  !
  !...........................................................................
  ! nyear(1) - year of simulation 
  !...........................................................................

  type(date), intent(in) :: newdate
  integer :: isec           ! index over emission sectors
  integer :: iemis          ! index over emissions (so2,nox,..)
  integer :: iland          ! index over countries 
  integer :: nmnd, nmnd2, nmnd0   ! this month, next month, preceding month.
  integer :: weekday        ! 1=Monday, 2=Tuesday etc.
  real    :: xday           ! used in interpolation
  integer :: yyyy,dd 
  real :: Start, Endval, Average, x

 !-----------------------------

   yyyy=newdate%year
   nmnd=newdate%month
   dd=newdate%day

   weekday = day_of_week(yyyy,nmnd,dd)
   if ( weekday == 0 ) weekday = 7  ! restores Sunday to 7

!   Parameters for time interpolation

    nmnd2 = nmnd + 1                   ! Next month
    if( nmnd2 > 12 ) nmnd2 = 1         ! December+1 => January
    nmnd0 = nmnd - 1                   ! preceding month
    if( nmnd0 < 1 ) nmnd0 = 12         ! January-1 => December

    xday = real( newdate%day - 1 ) / real( nmdays(nmnd) ) 

!   Calculate monthly and daily factors for emissions 

    do iemis = 1, NEMIS_FILE
      do isec = 1, N_TFAC
         do iland = 1, NLAND
            Average=fac_emm(iland ,nmnd,isec,iemis)
            if(InterpolateMonthEmis)then
               Start= 0.5*(fac_emm(iland ,nmnd0,isec,iemis)+fac_emm(iland ,nmnd,isec,iemis))
               Endval= 0.5*(fac_emm(iland ,nmnd,isec,iemis)+fac_emm(iland ,nmnd2,isec,iemis))
               !limit values, to ensure that x never can be negative
               Start=min(Start,2*Average,2*fac_emm(iland ,nmnd0,isec,iemis))
               Endval=min(Endval,2*Average,2*fac_emm(iland ,nmnd2,isec,iemis))
               call Averageconserved_interpolate(Start,Endval,Average,nmdays(nmnd),dd,x)
            else
               x=Average
            endif
            timefac(iland,isec,iemis) = x *  fac_edd(iland,weekday,isec,iemis) 
 
         end do ! iland  
      end do ! isec   
   end do ! iemis 

 end subroutine NewDayFactors
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 subroutine DegreeDayFactors(daynumber)

!.....................................................................
!**    DESCRIPTION:
!   Generally called with daynumber, and then reads the gridded degree-day
!   based factors for emissions.
!   If called with daynumber = 0, just checks existance of file. If not
!   found, can use default country-based (GENEMIS) factors.

    integer, intent(in) :: daynumber   ! daynumber (1,..... ,365)
   
    integer, save :: dd_old = -1
    integer,dimension(2)  :: ijloc   ! debug only 
    integer :: iii, jjj              ! debug only 
    real :: checkmax
    character(len=80) :: errmsg
    real, dimension(IIFULLDOM,JJFULLDOM) :: var2d_global
    integer :: kmax=1, nfetch=1 ! for HDD

!      Gridded_SNAP2_Factors = .false.
!      return

   !/ See if we have a file to work with....
    if ( daynumber == 0 ) then
       if(USES%WRF_MET_NAMES)then
          !we do not require that HDD file has been calculated
          call check_file(trim(DegreeDayFactorsFile), Gridded_SNAP2_Factors,&
               needed=.false., errmsg=errmsg ) ! do not use if file not found
       else
          !the file should be calculated if it does not yet exist
          call check_file(trim(DegreeDayFactorsFile), Gridded_SNAP2_Factors,&
               needed=.true., errmsg=errmsg )
       endif
       if ( Gridded_SNAP2_Factors ) then
          call PrintLog("Found "//trim(DegreeDayFactorsFile), MasterProc)
       else
          call PrintLog("Not-found: "//trim(DegreeDayFactorsFile)//' '//trim(errmsg), MasterProc)
          USES%DEGREEDAY_FACTORS = .false.
          if(me==0)write(*,*)'Warning: cannot use DEGREEDAY_FACTORS because file not found' 
       end if
       return
    end if

    if (dbgTF .and. MasterProc ) then
      write(*,*) 'MFAC DD?', USES%DEGREEDAY_FACTORS,  Gridded_SNAP2_Factors
    end if
    !===============================================
    if ( .not. Gridded_SNAP2_Factors )  return !
    !===============================================

   !/ We have a file, calculate every day ... .

    if (dd_old == daynumber) return   ! Only calculate once per day max
    dd_old= daynumber

!     write(*,*) "HDD inputs", me, " Day ", daynumber

   ! DegreeDays have the same domain/grid as the met data, so we can use:
    if(MasterProc) call GetCDF('HDD_Facs',trim(DegreeDayFactorsFile), &
          var2d_global,IIFULLDOM,JJFULLDOM,1,daynumber,nfetch)

    if(.not.allocated(gridfac_HDD))then
       allocate(gridfac_HDD(MAXLIMAX,MAXLJMAX))
    end if

    call global2local(var2d_global,gridfac_HDD,MSG_READ8,1,IIFULLDOM,JJFULLDOM,&
         kmax,IRUNBEG,JRUNBEG)
         call CheckStop(errmsg=="field_not_found", "INDegreeDay field not found:")

    if ( dbgTF ) then
       ijloc = maxloc( gridfac_HDD(li0:li1,lj0:lj1))
       iii = ijloc(1)+li0-1
       jjj = ijloc(2)+lj0-1
       checkmax = maxval( gridfac_HDD(li0:li1,lj0:lj1))

       write(*,"(a,2i4,2f10.2,20i4)") "DEBUG GRIDFAC MAx", me, daynumber, &
           checkmax, gridfac_HDD(iii,jjj), & !!! tmpt2(iii,jjj), &
             ijloc(1), ijloc(2), i_fdom(iii), j_fdom(jjj)
  
       if( debug_proc ) then
           write(*,"(a,i4,f12.3)") "GRIDFACDAY ", daynumber, &
             gridfac_HDD(debug_li,debug_lj)
       end if
    end if


    if ( dbgTF .and. debug_proc ) then
       iii = debug_li
       jjj = debug_lj
       write(*,*) "DEBUG GRIDFAC", me, daynumber, iii, jjj, gridfac_HDD(iii, jjj)
    end if


   end subroutine DegreeDayFactors

   subroutine Read_monthly_emis_grid_fac(month)

     implicit none
     integer, intent(in) ::month
     integer ::iemis,isec, i
     character(len=20) ::sector_map(NSECTORS,NEMIS_FILE),name
     character(len=3)  :: src ! e.g. ene, ind, 
     character(len=*), parameter:: dtxt='gridtfacs:'
     logical :: GRIDDED_CAMS_TEMPO = .false.
     real :: meanTfac
     integer :: idbg=3, jdbg=3 ! FAKE, since debug_li not set yet
     if ( dbgTF ) write(*,*) dtxt//': '// trim(GriddedMonthlyFacFile)
     if ( index(GriddedMonthlyFacFile,'CAMS_TEMPO')>0 ) then
       if ( dbgTF ) write(*,*) dtxt//' GRIDDED_CAMS_TEMPO'
       GRIDDED_CAMS_TEMPO = .true.
     else
       call CheckStop(N_TFAC/=11, &
         "Only snap timefactors implemented for ECLIPSEv5 gridded timefactors")
     end if
     if ( dbgTF ) then
       do iemis=1,NEMIS_FILE
         write(*,*) dtxt//'EMIS:'//trim(EMIS_File(iemis))
       end do
     end if

     ! sector_map(sector,emis) = name_in_netcdf_file
     ! CRUDELY HARD-CODED for Jan 2024 CAMS-TEMP_GLOB FILE
     sector_map(:,:)='default'
     do i = 1, NSECTORS
        if(IS_POW(i)) then
           src = 'ene'
           !sector_map(i,:)=src  ! no nh3
           do iemis=1,NEMIS_FILE
            if(trim(EMIS_File(iemis))=='voc')  sector_map(i,iemis)=src//'_nmvoc'
            if(trim(EMIS_File(iemis))=='co')   sector_map(i,iemis)=src//'_co'
            if(trim(EMIS_File(iemis))=='nox')  sector_map(i,iemis)=src//'_nox'
            if(trim(EMIS_File(iemis))=='sox')  sector_map(i,iemis)=src//'_sox'
            if(trim(EMIS_File(iemis))=='pm25') sector_map(i,iemis)=src//'_pm2.5'
            if(trim(EMIS_File(iemis))=='pmco') sector_map(i,iemis)=src//'_pm10'
           end do
        else if(IS_DOM(i))  then
           sector_map(i,:)='res'   ! was 'dom'
        else if(IS_TRAF(i)) then
           src = 'tro' 
           !sector_map(i,:)=src  ! no sox
           do iemis=1,NEMIS_FILE
            if(trim(EMIS_File(iemis))=='voc')  sector_map(i,iemis)=src//'_nmvoc'
            if(trim(EMIS_File(iemis))=='co')   sector_map(i,iemis)=src//'_co'
            if(trim(EMIS_File(iemis))=='nox')  sector_map(i,iemis)=src//'_nox'
            if(trim(EMIS_File(iemis))=='nh3')  sector_map(i,iemis)=src//'_nh3'
            if(trim(EMIS_File(iemis))=='pm25') sector_map(i,iemis)=src//'_pm2.5'
            if(trim(EMIS_File(iemis))=='pmco') sector_map(i,iemis)=src//'_pm10'
            if ( dbgTF ) write(*,*) dtxt//' TRAF:'//trim(sector_map(i,iemis)), i, iemis
           end do
        else if(IS_AGRK(i)) then
           src = 'agl'   ! GNFR K
           do iemis=1,NEMIS_FILE
            if(trim(EMIS_File(iemis))=='nox')  sector_map(i,iemis)=src//'_nox'
            if(trim(EMIS_File(iemis))=='nh3')  sector_map(i,iemis)=src//'_nh3'
           end do
        else if(IS_AGRL(i)) then
           src = 'ags'   ! GNFR L inc AWB!!!
           do iemis=1,NEMIS_FILE
            if(trim(EMIS_File(iemis))=='nox')  sector_map(i,iemis)=src//'_nox'
            if(trim(EMIS_File(iemis))=='nh3')  sector_map(i,iemis)=src//'_nh3'
           end do
        else if(IS_AGR(i)) then
           src = 'agl'   ! FIX K n L n AWB!!!
           !sector_map(i,:)= src
           do iemis=1,NEMIS_FILE
            if(trim(EMIS_File(iemis))=='nox')  sector_map(i,iemis)=src//'_nox'
            if(trim(EMIS_File(iemis))=='nh3')  sector_map(i,iemis)=src//'_nh3'
           end do
!SNAP        sector_map(3,:)='ind'
!SNAP        sector_map(4,:)='ind'
        else if(IS_IND(i)) then
              sector_map(i,:)='ind'
        end if
        if ( dbgTF ) then
           write(*,*) dtxt//' SECS:'//trim(sector_map(i,1)), i
           do iemis=1, NEMIS_FILE
             write(*, *) dtxt//'POW', iemis, sector_map(1,iemis)
           end do
        end if
     end do

     if(.not.allocated(GridTfac))then
        allocate(GridTfac(LIMAX,LJMAX,NSECTORS,NEMIS_FILE))
        GridTfac=dble(nmdays(month))/nydays !default, multiplied by inverse later!!
        if(dbgTF) write(*,*)'IJa', month, nmdays(month), maxval(GridTfac), minval(GridTfac)
     end if

     name='none'
     do isec=1,NSECTORS
        do iemis=1,NEMIS_FILE
           if ( dbgTF )write(*,*) dtxt//'maps', isec, iemis, sector_map(isec,iemis)

           if(sector_map(isec,iemis)=='default')then
              GridTfac(:,:,isec,iemis)=dble(nmdays(month))/nydays!default, multiplied by inverse later!!
              name=sector_map(isec,iemis)
              cycle
           end if
           if(sector_map(isec,iemis)==name.and.iemis>1)then
              !has same values as before, no need to read again
              GridTfac(:,:,isec,iemis)=GridTfac(:,:,isec,iemis-1)
              !if(dbgTF) write(*,*)'IJc', isec,iemis, GridTfac(idbg,jdbg,1,1)
           else
              
              name=sector_map(isec,iemis)
              
              !print *, dtxt//'TESTS'//trim(name), isec, iemis, trim(EMIS_FILE(iemis))
              call ReadField_CDF(trim(GriddedMonthlyFacFile),&
                   name,GridTfac(:,:,isec,iemis),month,interpol='conservative',&
                   known_projection='lon lat',needed=.true.,debug_flag=.false.,&
                   Undef=real(nmdays(month))/nydays )!default, multiplied by inverse later!!

              !Feb 7
              if ( GRIDDED_CAMS_TEMPO ) then  ! has values  ca. 1.0
                GridTfac(:,:,isec,iemis)=GridTfac(:,:,isec,iemis)*dble(nmdays(month))/nydays
                ! now has ca. 0.08, ~ 1/12
                ! if(dbgTF) write(*,*)'IJd', isec,iemis, GridTfac(idbg,jdbg,1,1)
              end if

           end if
        end do
     end do

!normalizations:
! in ECLIPSEv5_monthly_patterns.nc the "default" timefactors are defines as
! nmdays(month)/nydays
! In CAMS_TEMPO_GLOB 4emep the default timefactprs are 1.0
   
!the normalization until here is such that GridTfac gives the relative contribution from each month
!However we want to use it as a multiplicative factor, which gives the total as the sum (integral) of the
!factor over the number of days. Therefore as a multiplicative factor it has to be divided by the number 
!of days in the month. 

      !if  ( .not.GRIDDED_CAMS_TEMPO ) then
                        ! LIMAX,LJMAX,NSECTORS,NEMIS_FILE
      meanTfac= sum(GridTfac(:,:,1,1)/(LIMAX*LJMAX) )  ! test if 1.0 or 1/12 or ..
      if  (  abs(meanTfac-1.0) > 1.0e-2 ) then ! .not.GRIDDED_CAMS_TEMPO. Maybe ECLIPSEv5?
        if ( MasterProc) write(*,*) dtxt//'RE-NORMALISE', month, meanTfac, maxval(GridTfac(:,:,1,1))
        !if(dbgTF) write(*,*)dtxt//'IJe:',isec, iemis, GridTfac(idbg,jdbg,1,1)   ! 1/12
        GridTfac = GridTfac*nydays/nmdays(month)
        if(dbgTF) write(*,*)dtxt//'IJf:',isec, iemis, month,  GridTfac(idbg,jdbg,1,1) ! 1.0
      end if

!The normalization now, gives for instance GridTfac = 1 for constant emissions.

      end subroutine Read_monthly_emis_grid_fac

subroutine yearly_normalize(year)

  implicit none

  integer, intent(in) :: year
  integer :: iemis, n, isec, indexCC, iday, mm, idd, weekday, mm2, mm0
  real :: sumfac, Start, Endval, Average, x
  character(len=*), parameter::dtxt='TFyrNorm:'

! #######################################################################
! 4) Normalise the monthly-daily factors. This is needed in order to
!    account for leap years (nydays=366) and for the fact that different
!    years have different numbers of e.g. Saturdays/Sundays. 
!    Here we execute the same interpolations which are later done
!    in "NewDayFactors", and scale efac_mm if necessary.
  
   !F24 if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)then
   if(timeFacs%monthly == 'GRIDDED')then
      write(*,*)'Normalizing monthly emission time factors'
      fac_emm=1.0
      !enforce a constant integral of daily timefactors over each month
      !note that the subsequent interpolation does not change this integral
      do iemis = 1, NEMIS_FILE
         n = 0
         do isec = 1, N_TFAC
            do indexCC = 1, NLAND
               iday = 0
               do mm = 1, 12     ! Jan - Dec
                  sumfac = 0.0
                  do idd = 1, nmdays(mm)
                     
                     weekday=day_of_week (year,mm,idd)
                     
                     if ( weekday == 0 ) weekday = 7  ! restores Sunday to 7
                     sumfac = sumfac + fac_edd(indexCC,weekday,isec,iemis)   
                     
                  end do ! idd
                  ! redefine monthly factor to enforce this
                  fac_emm(indexCC,mm,isec,iemis)=nmdays(mm)/sumfac

               end do ! mm
               
            end do ! indexCC
       end do ! isec

      end do ! iemis

   end if

! normalize the factors over the year 
   do iemis = 1, NEMIS_FILE
       n = 0
       do isec = 1, N_TFAC
           do indexCC = 1, NLAND
             iday = 0
             sumfac = 0.0

             do mm = 1, 12     ! Jan - Dec
                do idd = 1, nmdays(mm)

                   weekday=day_of_week (year,mm,idd)

                   if ( weekday == 0 ) weekday = 7  ! restores Sunday to 7

                   mm2 = mm + 1 
                   if( mm2  > 12 ) mm2 = 1          ! December+1 => January
                   mm0 = mm - 1 
                   if( mm0 < 1   ) mm0 = 12          ! January-1 => December
                   Start= 0.5*(fac_emm(indexCC,mm0,isec,iemis)+fac_emm(indexCC,mm,isec,iemis))
                   Endval= 0.5*(fac_emm(indexCC,mm,isec,iemis)+fac_emm(indexCC,mm2,isec,iemis))                   
                   Average=fac_emm(indexCC,mm,isec,iemis)
                   !limit values, to ensure that x never can be negative
                   Start=min(Start,2*Average,2*fac_emm(indexCC,mm0,isec,iemis))
                   Endval=min(Endval,2*Average,2*fac_emm(indexCC,mm2,isec,iemis))
                   call Averageconserved_interpolate(Start,Endval,Average,nmdays(mm),idd,x)
                   sumfac = sumfac + x * fac_edd(indexCC,weekday,isec,iemis)   

                end do ! idd
             end do ! mm

             sumfac = real(nydays)/sumfac    

            ! if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)then
             if(timeFacs%monthly == 'GRIDDED')then
                ! should already almost be normalized (almost, because there is still some 
                ! variation left due to the differences occuring when week-ends are in the 
                ! middle or end of the month (linear_interpolation*day_factor is not linear)
               if ( sumfac < 0.999 .or. sumfac > 1.001 ) then
                 write(unit=errmsg,fmt=*) &
                   "GRIDDED Time-factor error! for ",iemis, isec, indexCC," sumfac ",sumfac
                 call CheckStop(errmsg)
              end if
               
             end if
              if ( sumfac < 0.99 .or. sumfac > 1.01 )write(*,*)'sumfac: ',iemis,isec,indexCC,sumfac   

              if ( sumfac < 0.97 .or. sumfac > 1.03 ) then
                 write(unit=errmsg,fmt=*) &
                   "Time-factor error! for ",iemis, isec, indexCC," sumfac ",sumfac
                 call CheckStop(errmsg)
              end if

             ! can be caused by variable number of sundays in a month for instance
             ! Slight adjustment of monthly factors
              do mm = 1, 12
                 fac_emm(indexCC,mm,isec,iemis)  =  &
                       fac_emm(indexCC,mm,isec,iemis) * sumfac
              end do ! mm
              !DSTMP if ( dbgTF .and. abs(sumfac-1.0)>0.001) then
              if ( dbgTF .and. abs(sumfac-1.0)>0.01) then
                write(*,fmt='(a,3i4,f14.6)')  dtxt// &
                 "needed for country, isec, iemis, sumfac = " , &
                  indexCC, isec, iemis, sumfac
              end if ! dbgTF
          end do ! indexCC
       end do ! isec


      end do ! iemis


end subroutine yearly_normalize

end module Timefactors_mod

