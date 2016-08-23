! <Timefactors_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2012 met.no
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

                    module Timefactors_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!.......................................................................
!  DESCRIPTION:
!  Calculates emission temporal variation.
!  Reads monthly and daily (GENEMIS) factors for all emissions from files
!  Monthly.sox, Daily.sox, Monthly.nox, etc., -> in fac_emm, fac_edd arrays 
!  For every day, calculates emission factor "timefac" per country, emission 
!  sector, emission component
!  
!  Sets the day/night emissions variation in day_factor
!
!  D. Simpson,    3/2/99-11 0 H. Fagerli, 2011
!_____________________________________________________________________________

  use CheckStop_ml, only : CheckStop
  use Country_ml,   only : NLAND
  use EmisDef_ml,   only : NSECTORS, NEMIS_FILE, EMIS_FILE, ISNAP_DOM
  use GridValues_ml    , only : i_fdom,j_fdom, debug_proc,debug_li,debug_lj
  use Met_ml,       only : Getmeteofield
  use ModelConstants_ml, only : MasterProc, DEBUG => DEBUG_EMISTIMEFACS
  use ModelConstants_ml, only : IIFULLDOM, JJFULLDOM
  use ModelConstants_ml, only : iyr_trend
  use ModelConstants_ml, only : INERIS_SNAP1, INERIS_SNAP2
  use NetCDF_ml,    only : GetCDF 
  use Par_ml,       only : MAXLIMAX,MAXLJMAX, me, li0, lj0, li1, lj1
  use Par_ml,       only : IRUNBEG, JRUNBEG, MSG_READ8
  use PhysicalConstants_ml, only : PI
  use Io_ml,        only :            &
                     open_file,       & ! subroutine
                     check_file,       & ! subroutine
                     PrintLog,        &
                     ios,  IO_TIMEFACS  ! i/o error number, i/o label
  use TimeDate_ml,  only:            &  ! subroutine, sets:
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

  !-- time factor stuff: 

  real, public, save, &
     dimension(NLAND,NSECTORS,NEMIS_FILE) :: timefac ! overall emission 
                                                      ! timefactor 
                                                      ! calculated daily
  real, public, save,  &
     dimension(NLAND,12,NSECTORS,NEMIS_FILE) :: fac_emm  ! Monthly factors

 ! Hourly for each day ! From EURODELTA/INERIS
  real, public, save,  &
     dimension(NSECTORS,24,7) :: fac_ehh24x7  !  Hour factors for 7 days

 ! We keep track of min value for degree-day work
 !
  real, public, save,  &
     dimension(NLAND,NSECTORS,NEMIS_FILE) :: fac_min ! Min of Monthly factors
 !
  real, public, save,  &
     dimension(12) :: fac_cemm  ! Change in monthly factors over the years

  real, public, save,  &
     dimension(NLAND, 7,NSECTORS,NEMIS_FILE) :: fac_edd  ! Daily factors

  ! Heating-degree day factor for SNAP-2. Independent of country:
  logical, public, save :: Gridded_SNAP2_Factors = .false.
  real, public, allocatable,dimension (:,:), save :: gridfac_HDD
  !real, private, dimension (MAXLIMAX,MAXLJMAX), save :: tmpt2

  ! Used for general file calls and mpi routines below

  character(len=30), private :: fname2   ! input filename - do not change 

contains


  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine timefactors(year)

   !.......................................................................
   !  DESCRIPTION:
   !  Read in monthly and daily factors, -> fac_emm, fac_edd arrays
   !  The input files are Monthly.sox, Daily.sox, Monthly.nox, etc.
   !  Sets the day/night variation in day_factor
   !
   !  D. Simpson,    3/2/99
   !.......................................................................

  !-- Input
  integer, intent(in) :: year

  !-- Outputs -  module's fac_emm, fac_edd, day_factor, etc.

  !-- local
  integer ::  inland, insec     ! Country and sector value read from femis
  integer ::  i, ic, isec, n 
  integer ::  idd, idd2, ihh, iday, mm, mm2 ! Loop and count variables
  integer ::  iemis             ! emission count variables

  integer :: weekday            ! 1=monday, 2=tuesday etc.
  real    :: xday, sumfac       ! used in interpolation, testing
  real    :: tmp24(24)          ! used for hourly factors
  character(len=100) :: errmsg
  character(len=200) :: inputline
  real :: fracchange
  real, dimension(NLAND,NEMIS_FILE):: sumfacc !factor to normalize monthly changes                                                         

  if (DEBUG) write(unit=6,fmt=*) "into timefactors "

   call CheckStop( nydays < 365, &
      "Timefactors: ERR:Call set_nmdays before timefactors?")

   call CheckStop(  NSECTORS /= 11 , &
      "Timefactors: ERR:Day-Night dimension wrong!")


!  #################################
!  1) Read in Monthly factors, and determine min value (for baseload)

   fac_emm(:,:,:,:) = 1.0
   fac_min(:,:,:) = 1.0

  ! Summer/winter SNAP1 ratios reduced from 1990 to 2010:
   fac_cemm(:) = 1.0
   fracchange=0.005*(iyr_trend -1990)
   fracchange=max(0.0,fracchange) !do not change before 1990
   fracchange=min(0.1,fracchange) !stop change after 2010 
                                  !equal 1.1/0.9=1.22 summer/winter change
   write(unit=6,fmt=*) "Change summer/winter ratio in SNAP1 by ", fracchange

   do mm=1,12
      !Assume max change for august and february
      fac_cemm(mm)  = 1.0 + fracchange * cos ( 2 * PI * (mm - 8)/ 12.0 )
      write(unit=6,fmt="(a,i3,f8.3,a,f8.3)") "Change in fac_cemm ", mm,fac_cemm(mm)
   enddo
   write(*,"(a,f8.4)") "Mean fac_cemm ", sum( fac_cemm(:) )/12.0

   if( INERIS_SNAP1 ) fac_cemm(:) = 1.0


   do iemis = 1, NEMIS_FILE

       fname2 = "MonthlyFac." // trim ( EMIS_FILE(iemis) )
       call open_file(IO_TIMEFACS,"r",fname2,needed=.true.)

       call CheckStop( ios, &
        "Timefactors: IOS error in Monthlyfac")

       n = 0
       do 
           read(IO_TIMEFACS,fmt=*,iostat=ios) inland, insec, &
             (fac_emm(inland,mm,insec,iemis),mm=1,12)
           if ( ios <  0 ) exit     ! End of file

           !defined after renormalization and send to al processors:
           ! fac_min(inland,insec,iemis) = minval( fac_emm(inland,:,insec,iemis) )

           if( DEBUG.and.insec==ISNAP_DOM  ) &
              write(*,"(a,3i3,f7.3,a,12f6.1)") "emm tfac ", &
               inland,insec,iemis, fac_min(inland,insec,iemis),&
                 " : ",  ( fac_emm(inland,mm,insec,iemis), mm=1,12)

           call CheckStop( ios, "Timefactors: Read error in Monthlyfac")

           n = n + 1
       enddo

       close(IO_TIMEFACS)

! Apply change in monthly factors for SNAP 1
          sumfacc(:,:)=0.0
          do ic = 1, NLAND
             do mm=1,12
                fac_emm(ic,mm,1,iemis)=fac_emm(ic,mm,1,iemis)*fac_cemm(mm)
                sumfacc(ic,iemis)=sumfacc(ic,iemis)+fac_emm(ic,mm,1,iemis)
             enddo
          enddo
! normalize
          do ic = 1, NLAND
             do mm=1,12
                fac_emm(ic,mm,1,iemis)=fac_emm(ic,mm,1,iemis)*12./sumfacc(ic,iemis)
             enddo
          enddo
       if (DEBUG) write(unit=6,fmt=*) "Read ", n, " records from ", fname2 
   enddo  ! iemis


! #################################
! 2) Read in Daily factors

  fac_edd(:,:,:,:) = 1.0

  do iemis = 1, NEMIS_FILE

       fname2 = "DailyFac." // trim ( EMIS_FILE(iemis) )
       call open_file(IO_TIMEFACS,"r",fname2,needed=.true.)

       call CheckStop( ios, &
        "Timefactors: Opening error in Dailyfac")

       n = 0
       do
         read(IO_TIMEFACS,fmt=*,iostat=ios) inland, insec, &
             (fac_edd(inland,i,insec,iemis),i=1,7)
           if ( ios <  0 ) exit   ! End of file

           call CheckStop( ios, "Timefactors: Read error in Dailyfac")

           n = n + 1

           !-- Sum over days 1-7
           xday =  sum( fac_edd(inland,1:7,insec,iemis) ) / 7.0

           call CheckStop( xday > 1.001 .or. xday < 0.999, &
                "Timefactors: ERROR: Dailyfac - not normalised")

       enddo

       close(IO_TIMEFACS)
       if (DEBUG) write(unit=6,fmt=*) "Read ", n, " records from ", fname2

  enddo  ! NEMIS_FILE

!  #################################
!  3) Read in hourly (24x7) factors, options set in run script.
   ! INERIS option has 11x24x7 emissions factor
   ! TNO2005 option has 11x24 
   ! EMEP2003 option has very simple day night
!
       fname2 = "HOURLY-FACS"  ! From EURODELTA/INERIS/TNO or EMEP2003
       write(unit=6,fmt=*) "Starting HOURLY-FACS"
       call open_file(IO_TIMEFACS,"r",fname2,needed=.true.)

       fac_ehh24x7 = -999.

       n = 0
       do 
          read(IO_TIMEFACS,"(a)",iostat=ios) inputline
          n = n + 1
          if(DEBUG)write(*,*) "HourlyFacs ", n, trim(inputline)
          if ( ios <  0 ) exit     ! End of file
          if( index(inputline,"#")>0 ) then ! Headers
            if(n==1) call PrintLog(trim(inputline))
            cycle
          else
            read(inputline,fmt=*,iostat=ios) idd, insec, &
              (tmp24(ihh),ihh=1,24)
            if( DEBUG ) write(*,*) "HOURLY=> ",idd, insec, tmp24(1), tmp24(13)
          end if

            if(  idd == 0 ) then ! same values very day
              do idd2 = 1, 7
                fac_ehh24x7(insec,:,idd2) = tmp24(:)
              end do
              idd = 1 ! Used later
            else
                fac_ehh24x7(insec,:,idd) = tmp24(:)
            end if

             !(fac_ehh24x7(insec,ihh,idd),ihh=1,24)

           ! Use sumfac for mean, and normalise within each day/sector
           ! (Sector 10 had a sum of 1.00625)
           sumfac = sum(fac_ehh24x7(insec,:,idd))/24.0
           if(DEBUG .and. MasterProc) write(*,"(a,2i3,3f12.5)") &
              'HOURLY-FACS mean min max', idd, insec, sumfac, &
                minval(fac_ehh24x7(insec,:,idd)), &
                maxval(fac_ehh24x7(insec,:,idd))

           fac_ehh24x7(insec,:,idd) = fac_ehh24x7(insec,:,idd) * 1.0/sumfac

       !    if ( ios <  0 ) exit     ! End of file
       end do
       !do insec=1, 11; do idd   =1, 7; do ihh   =1, 24
       !    if( fac_ehh24x7(insec,ihh,idd)  < 0.0 ) then 
       !       print *, "Unfilled ", insec, idd, ihh, fac_ehh24x7(insec,ihh,idd)
       !    end if
       !end do; end do; end do
       !call CheckStop ( any(fac_ehh24x7 < 0.0 ) , "Unfilled efac_ehh24x7")
       if (DEBUG) write(unit=6,fmt=*) "Read ", n, " records from ", fname2
       call CheckStop ( any(fac_ehh24x7 < 0.0 ) , "Unfilled efac_ehh24x7")


! #######################################################################
! 4) Normalise the monthly-daily factors. This is needed in order to
!    account for leap years (nydays=366) and for the fact that different
!    years have different numbers of e.g. Saturdays/Sundays. 
!    Here we execute the same interpolations which are later done
!    in "NewDayFactors", and scale efac_mm if necessary.


  write(unit=6,fmt="(a,I6,a,I5)")" Time factors normalisation: ",nydays,' days in ',year 

  do iemis = 1, NEMIS_FILE
       n = 0
       do isec = 1, NSECTORS
           do ic = 1, NLAND
             iday = 0
             sumfac = 0.0

             do mm = 1, 12     ! Jan - Dec
                do idd = 1, nmdays(mm)

                   weekday=day_of_week (year,mm,idd)

                   if ( weekday == 0 ) weekday = 7  ! restores Sunday to 7

                   mm2 = mm + 1 
                   if( mm2  > 12 ) mm2 = 1          ! December+1 => January

                   xday = real(idd-1) /real(nmdays(mm))

                   sumfac = sumfac +                            &  ! timefac 
                      ( fac_emm(ic,mm,isec,iemis) +             &
                        ( fac_emm(ic,mm2,isec,iemis)            &
                         - fac_emm(ic,mm,isec,iemis) ) * xday ) &
                      * fac_edd(ic,weekday,isec,iemis)   

                end do ! idd
             end do ! mm

             sumfac = real(nydays)/sumfac    


              if ( sumfac < 0.97 .or. sumfac > 1.03 ) then
                 write(unit=errmsg,fmt=*) &
                   "Time-factor error! for ",iemis, isec, ic," sumfac ",sumfac
                 call CheckStop(errmsg)
              end if

             if ( sumfac < 0.999 .or. sumfac > 1.001 ) then
                 n = n+1
             ! Slight adjustment of monthly factors
                  do mm = 1, 12
                    fac_emm(ic,mm,isec,iemis)  =  &
                       fac_emm(ic,mm,isec,iemis) * sumfac
                  end do ! mm
             end if

          end do ! ic
       enddo ! isec

       if ( n ==  0 ) &
            write(unit=6,fmt=*)  &
           "Correction not needed for iemis, sumfac = " ,iemis, sumfac

      enddo ! iemis


!#########################################################################
!
    if (DEBUG) write(unit=6,fmt=*) "End of subroutine timefactors"

    if (DEBUG ) then 
       write( *,*) " test of time factors, UK: "
       do mm = 1, 12
           write(*, "(i2,i6,f8.3,3f8.4)") mm, nydays, sumfac,  &
            fac_emm(27,mm,2,1), fac_edd(27,1,2,1), fac_edd(27,7,2,1)
       end do ! mm
       write(*,"(a,4f8.3)") " day factors traffic 24x7", &
           fac_ehh24x7(7,1,4),fac_ehh24x7(7,13,4), &
              minval(fac_ehh24x7), maxval(fac_ehh24x7)
    end if ! DEBUG

 end subroutine timefactors
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine NewDayFactors(newdate)
  
  ! Calculates the monthly and daily factors for emission temporal variation
  ! for each country, emission, and sector.  Called at midnight every day.
  !
  ! Uses arays: 
  !     fac_emm(NLAND,NM,NSECTORS,NEMIS_FILES)    ! Jan - Dec.
  !     fac_edd(NLAND,7,NSECTORS,NEMIS_FILES)     ! Monday=1, Sunday=7
  !
  ! Outputs:
  !    real timefac(NLAND,NSECTORS,NEMIS_FILES)
  !
  !...........................................................................
  ! nyear(1) - year of simulation 
  !...........................................................................

  type(date), intent(in) :: newdate
  integer :: isec           ! index over emission sectors
  integer :: iemis          ! index over emissions (so2,nox,..)
  integer :: iland          ! index over countries 
  integer :: nmnd, nmnd2    ! this month, next month.
  integer :: weekday        ! 1=Monday, 2=Tuesday etc.
  real    :: xday           ! used in interpolation
  integer :: yyyy,dd 

 !-----------------------------

   yyyy=newdate%year
   nmnd=newdate%month
   dd=newdate%day

   weekday = day_of_week(yyyy,nmnd,dd)
   if ( weekday == 0 ) weekday = 7  ! restores Sunday to 7

!   Parameters for time interpolation

    nmnd2 = nmnd + 1                   ! Next month
    if( nmnd2 > 12 ) nmnd2 = 1         ! December+1 => January

    xday = real( newdate%day - 1 ) / real( nmdays(nmnd) ) 

!   Calculate monthly and daily factors for emissions 

    do iemis = 1, NEMIS_FILE
      do isec = 1, NSECTORS 
         do iland = 1, NLAND

             timefac(iland,isec,iemis) =                           &
                ( fac_emm(iland,nmnd,isec,iemis)  +                &
                   ( fac_emm(iland,nmnd2,isec,iemis) -             &
                      fac_emm(iland,nmnd,isec,iemis ) ) * xday )   &
               *  fac_edd(iland,weekday,isec,iemis) 

         enddo ! iland  
      enddo ! isec   
   enddo ! iemis 

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
    character(len=80) :: errmsg, units, varname
    real, dimension(IIFULLDOM,JJFULLDOM) :: var2d_global
    integer :: kmax=1, nfetch=1 ! for HDD

!      Gridded_SNAP2_Factors = .false.
!      return

   !/ See if we have a file to work with....
    if ( daynumber == 0 ) then
      call check_file("DegreeDayFactors.nc", Gridded_SNAP2_Factors,&
        needed=.false., errmsg=errmsg )
      if ( Gridded_SNAP2_Factors ) then
         call PrintLog("Found DEGREE-day factors", MasterProc)
      else
         call PrintLog("Not-found: DEGREE-day factors", MasterProc)
      end if
      return
    end if

    !===============================================
    if ( .not. Gridded_SNAP2_Factors )  return !
    !===============================================

   !/ We have a file, calculate every day ... .

    if (dd_old == daynumber) return   ! Only calculate once per day max
    dd_old= daynumber

!     write(*,*) "HDD inputs", me, " Day ", daynumber

   ! DegreeDays have the same domain/grid as the met data, so we can use:
    if(MasterProc) call GetCDF('HDD_Facs','DegreeDayFactors.nc', &
          var2d_global,IIFULLDOM,JJFULLDOM,1,daynumber,nfetch)

    if(.not.allocated(gridfac_HDD))then
       allocate(gridfac_HDD(MAXLIMAX,MAXLJMAX))
    endif

    call global2local(var2d_global,gridfac_HDD,MSG_READ8,1,IIFULLDOM,JJFULLDOM,&
         kmax,IRUNBEG,JRUNBEG)
         call CheckStop(errmsg=="field_not_found", "INDegreeDay field not found:")

    if ( DEBUG ) then
       ijloc = maxloc( gridfac_HDD(li0:li1,lj0:lj1))
       iii = ijloc(1)+li0-1
       jjj = ijloc(2)+lj0-1
       checkmax = maxval( gridfac_HDD(li0:li1,lj0:lj1))

       write(*,"(a,2i4,2f10.2,20i4)") "DEBUG GRIDFAC MAx", me, daynumber, &
           checkmax, gridfac_HDD(iii,jjj), & !!! tmpt2(iii,jjj), &
             ijloc(1), ijloc(2), i_fdom(iii), j_fdom(jjj)
  
       if( debug_proc ) then
           !write(*,"(a,3i4,2f12.3)") "GRIDFACDAYGEN ", daynumber, &
           !  i_fdom(iii), j_fdom(jjj),  tmpt2(iii,jjj), gridfac_HDD(iii,jjj)
           write(*,"(a,i4,f12.3)") "GRIDFACDAY ", daynumber, &
             gridfac_HDD(debug_li,debug_lj)
       end if
    end if


    if ( DEBUG .and. debug_proc ) then
       iii = debug_li
       jjj = debug_lj
       !write(*,*) "DEBUG GRIDFAC", me, daynumber, iii, jjj, tmpt2(iii,jjj),  gridfac_HDD(iii, jjj)
       write(*,*) "DEBUG GRIDFAC", me, daynumber, iii, jjj, gridfac_HDD(iii, jjj)
    end if


   end subroutine DegreeDayFactors

end module Timefactors_ml

