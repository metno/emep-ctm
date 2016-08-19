! <Timefactors_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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
!**    DESCRIPTION:
!  Calculates emissions temporal variation.
!  Reads monthly and daily (GENEMIS) factors for all emissions from files
!  Monthly.sox, Daily.sox, Monthly.nox, etc., -> in fac_emm, fac_edd arrays 
!  For every day, calculates emission factor "timefac" per country, emission 
!  sector, emission component
!  
!  Sets the day/night emissions variation in day_factor
!
!  D. Simpson,    3/2/99
!_____________________________________________________________________________
  use CheckStop_ml, only : CheckStop
  use Country_ml,   only : NLAND
  use My_Emis_ml,   only : NEMIS, EMIS_NAME
  use EmisDef_ml,   only : NSECTORS
  use TimeDate_ml,  only:            &  ! subroutine, sets:
                     date,           &  ! date-type definition
                     nmdays, nydays, &  ! days per month (12), days per year
                     day_of_week,&      ! weekday, 0=sun, 1=thuesday...
                     day_of_year        ! day count in year
  use Io_ml,        only :            &
                     open_file,       & ! subroutine
                     ios,  IO_TIMEFACS  ! i/o error number, i/o label

  implicit none
  private

  !-- subroutines:

  public :: NewDayFactors
  public :: timefactors

  !-- time factor stuff: 

  real, public, save, &
     dimension(NLAND,NSECTORS,NEMIS) :: timefac ! overall emission timefactor 
                                                ! calculated daily
  real, public, save,  &
     dimension(NLAND,12,NSECTORS,NEMIS) :: fac_emm  ! Monthly factors
  real, public, save,  &
     dimension(NLAND, 7,NSECTORS,NEMIS) :: fac_edd  ! Daily factors

  real, public, save, dimension(NSECTORS,0:1):: day_factor  ! Day/night factor 

  logical, private, parameter :: DEBUG = .false.

  !/** used for general file calls and mpi routines below **/

  character(len=30), private :: fname2   ! input filename - do not change 

contains


  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine timefactors(year)

   !.......................................................................
   !**    DESCRIPTION:
   !  Read in monthly and daily factors, -> fac_emm, fac_edd arrays
   !  The input files are Monthly.sox, Daily.sox, Monthly.nox, etc.
   !  Sets the day/night variation in day_factor
   !
   !      D. Simpson,    3/2/99
   !.......................................................................

  !--Input
  integer, intent(in) :: year

  !-- Outputs -  module's fac_emm, fac_edd, day_factor, etc.

  !-- local
  integer ::  inland, insec     ! Country and sector value read from femis
  integer ::  i, ic, isec, n, idd, iday, mm, mm2 ! Loop and count variables
  integer ::  iemis              ! emission count variables

  integer :: weekday         ! 1=monday, 2=tuesday etc.
  real    :: xday, sumfac    ! used in interpolation, testing
  character(len=100) :: errmsg


!/** Factor giving nighttime  emission ratio. 
! ** note this is hard-coded with NSECTORS=11. Checked in code

   real, parameter, dimension(NSECTORS) ::  & 
        DAY_NIGHT = (/      & 
                       1.0  &! 1.  Power production
                     , 0.8  &! 2.  Comm/res. combustion
                     , 0.8  &! 3.  Industrial combustion
                     , 1.0  &! 4.  Non-industrial combustion
                     , 1.0  &! 5.  Processes
                     , 0.5  &! 6.  Solvent use
                     , 0.5  &! 7.  Road transport
                     , 0.8  &! 8.  Other transport
                     , 1.0  &! 9.  Waste
                     , 1.0  &! 10. Agriculture
                     , 1.0  &! 11. Nature
                     /)

  if (DEBUG) write(unit=6,fmt=*) "into timefactors.f "

   call CheckStop( nydays < 365, &
      "Timefactors: ERR:Call set_nmdays before timefactors?")

   call CheckStop(  NSECTORS /= 11 , &
      "Timefactors: ERR:Day-Night dimension wrong!")


!  #################################
!  *** 1) Read in Monthly factors

   fac_emm(:,:,:,:) = 1.0

   do iemis = 1, NEMIS

       fname2 = "MonthlyFac." // trim ( EMIS_NAME(iemis) )
       call open_file(IO_TIMEFACS,"r",fname2,needed=.true.)

       call CheckStop( ios, &
        "Timefactors: IOS error in Monthlyfac")

       n = 0
       do 
           read(IO_TIMEFACS,fmt=*,iostat=ios) inland, insec, &
             (fac_emm(inland,mm,insec,iemis),mm=1,12)
           if ( ios <  0 ) exit     ! End of file

           call CheckStop( ios, "Timefactors: Read error in Monthlyfac")

           n = n + 1
       enddo

       close(IO_TIMEFACS)

       if (DEBUG) write(unit=6,fmt=*) "Read ", n, " records from ", fname2 
   enddo  ! iemis


! #################################
!CCC*** 2) Read in Daily factors

  fac_edd(:,:,:,:) = 1.0

  do iemis = 1, NEMIS

       fname2 = "DailyFac." // trim ( EMIS_NAME(iemis) )
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

  enddo  ! NEMIS

! #######################################################################
!cccc  3) Normalise the monthly-daily factors. This is needed in order to
!         account for leap years (nydays=366) and for the fact that different
!         years have different numbers of e.g. saturdays/sundays. 
!         Here we execute the same interpolations which are later done
!         in "NewDayFactors", and scale efac_mm if necessary


  write(unit=6,fmt=*) "Time factor interpolation "
  write(unit=6,fmt=*) "for nmdays(2) = ", nmdays(2), " gives nydays= ", nydays

  do iemis = 1, NEMIS
       n = 0
       do isec = 1, NSECTORS
           do ic = 1, NLAND
             iday = 0
             sumfac = 0.0

             do mm = 1, 12     ! Jan - Dec
                do idd = 1, nmdays(mm)

                   weekday=day_of_week (year,mm,idd)

                   if ( weekday == 0 ) weekday = 7  ! restores sunday to 7

                   mm2 = mm + 1 
                   if( mm2  > 12 ) mm2 = 1      ! December+1 => January

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
! Day/night factors are set from parameter DAY_NIGHT in emisdef_ml
!     daytime = 2 - nightime :

  day_factor(:,0)  =  DAY_NIGHT(:)             ! Night
  day_factor(:,1) = 2.0 - day_factor(:,0)      ! Day

!     #################################

    if (DEBUG) write(unit=6,fmt=*) "End of subroutine timefactors"

    if (DEBUG ) then 
       print *, " test of time factors, UK: "
       do mm = 1, 12
           print "(i2,i6,f8.3,3f8.4)", mm, nydays, sumfac,  &
            fac_emm(27,mm,2,1), fac_edd(27,1,2,1), fac_edd(27,7,2,1)
       end do ! mm
       print *, " day factors traffic are", day_factor(7,0), day_factor(7,1)
    end if ! DEBUG

 end subroutine timefactors
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine NewDayFactors(newdate)
  !
  !  Calculates the monthly and daily factors for emission temporal variation
  !  for each country, emission, and sector.  Called at midnight every day.
  !
  !  Uses arays: 
  !     fac_emm(NLAND,NM,NSECTORS,NEMIS)    ! Jan - Dec.
  !     fac_edd(NLAND,7,NSECTORS,NEMIS)     ! Monday=1, Sunday=7
  !
  ! Outputs:
  !    real timefac(NLAND,NSECTORS,NEMIS)
  !
  !...........................................................................
  !nyear(1) - year of simulation 
  !...........................................................................

  type(date), intent(in) :: newdate
  integer :: isec           ! index over emission sectors
  integer :: iemis          ! index over emissions (so2,nox,..)
  integer :: iland          ! index over countries 
  integer :: nmnd, nmnd2    ! this month, next month.
  integer :: weekday,nday,n ! 1=monday, 2=tuesday etc.
  real    :: xday           ! used in interpolation
  integer :: yyyy,dd 

 !-----------------------------

   yyyy=newdate%year
   nmnd=newdate%month
   dd=newdate%day

   weekday = day_of_week(yyyy,nmnd,dd)
   if ( weekday == 0 ) weekday = 7  ! restores sunday to 7

!   Parameters for time interpolation

    nmnd2 = nmnd + 1                   ! Next month
    if( nmnd2 > 12 ) nmnd2 = 1         ! December+1 => January

    xday = real( newdate%day - 1 ) / real( nmdays(nmnd) ) 

!   Calculate monthly and daily factors for emissions 

    do iemis = 1, NEMIS
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

end module Timefactors_ml

