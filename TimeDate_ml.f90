! <TimeDate_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
MODULE TimeDate_ml
        IMPLICIT NONE

! Originally timedate.f90 from Paul Curtis, found on web, 31/8/04: 
! Removed some Windows-specific or un-needed routines, 
! and converted to F. Some routines given longer names, eg. dow -> day_in_week,
! ndiy -> day_of_year. Change IFIX to INT, FLOAT to REAL, MAX0 to MAX, etc.
!===================Routines =================================================

 !/ Functions ...............
  public :: make_current_date      ! convert timestamp to current_date 
  public :: make_timestamp         ! convert current_date(yyyy,mon,day,hour,
                                   ! secs) to timestamp(jdate,secs)
  public  :: Init_nmdays   ! sets number of days per month, year

  public :: tdif_secs              ! t2-t1 -> dif (s)
  public :: ts_earlier             ! gets first ts from ts1, ts2
  public :: ts_later               ! gets later ts from ts1, ts2
  public :: InInterval             ! -1 == earlier, 0 == in, +1 == later
  public :: julian_date            ! yyyy,mm,dd -> julian
  public :: day_of_week            ! yyyy, mm, dd -> day of week (0=SUN)
  public :: day_of_year            ! yyyy,mm,dd -> Day count in year  
  public :: max_day                ! month,year -> maxd, e.g. 31 for July
  public :: leapyear               ! year -> true, false
  public :: y2dig                  ! year -> 2-digit yy
  public :: y4dig                  ! year -> 4-digit yyyy
  public :: get_ymd                ! jd -> yyyy, mm, dd
  public :: get_hms                ! secs -> hour,minute,second
 !/ Subroutines...............
  public :: dup_timestamp        ! ts2=ts1
  public :: add_secs             ! ts+seconds -> new ts. fixit option
  public :: add_month            ! jdate+month, force_day option
!===================TIMESTAMP TYPES & DEFINES=================================

!=========================================================================
integer, public,                save ::  daynumber    ! Day no. (1st jan=1).
  ! (Could be done locally?)
integer, public,                save ::  nydays       ! No. days per year
integer, public, dimension(12), save ::  nmdays       ! No. days per month

type, public :: date
  integer :: year
  integer :: month
  integer :: day
  integer :: hour
  integer :: seconds
 end type date


 type(date), public, save :: current_date

!==============================================================================


	TYPE, public :: timestamp
		INTEGER    :: jdate
		REAL       :: secs
	END TYPE timestamp

	REAL,private,PARAMETER	:: spd = 86400.0
	REAL,private,PARAMETER	:: sph =  3600.0
	REAL,private,PARAMETER	:: spm =    60.0

	TYPE(timestamp),private,PARAMETER  :: ts_null = timestamp(0, 0.0)
	TYPE(timestamp),public, save       :: ts_now   ! current local time
	TYPE(timestamp),public, save       :: ts_next   ! next inp

	CHARACTER(LEN=3),DIMENSION(12), public :: short_month =  &
		 (/"Jan","Feb","Mar","Apr","May","Jun",    &  
		   "Jul","Aug","Sep","Oct","Nov","Dec"/)
 
	CHARACTER(LEN=10),DIMENSION(12), public :: long_month = &
		 (/"January   ", "February  ",  "March     ",  &  
		   "April     ", "May       ",  "June      ",  &  
		   "July      ", "August    ",  "September ",  &  
		   "October   ", "November  ",  "December  "  /)  
	CHARACTER(LEN=3),DIMENSION(0:6), public :: short_day =  &
                 (/"Sun","Mon","Tue","Wed","Thu","Fri","Sat" /)

CONTAINS


	FUNCTION make_timestamp (cd) RESULT (ts)
		TYPE(timestamp)              :: ts
		TYPE(date),INTENT(IN)        :: cd
		INTEGER                      :: yyyy, mon, dd, hh, ss
                yyyy=cd%year
                mon=cd%month
                dd=cd%day
                hh=cd%hour
                ss=cd%seconds
		ts%jdate = julian_date (yyyy, mon, dd)
		ts%secs  = sph*REAL(hh) + REAL(ss)
	END FUNCTION make_timestamp


	FUNCTION make_current_date (ts) RESULT (cd)
		TYPE(timestamp),INTENT(IN)    :: ts
		TYPE(date)                    :: cd
                INTEGER                       :: yy,mm,dd,hh,min,sc,jd
                REAL                          :: ss
                jd=ts%jdate
                ss=ts%secs
                call get_ymd(jd,yy,mm,dd)
                call get_hms(ss,hh,min,sc)
                cd%year=yy
                cd%month=mm
                cd%day=dd
                cd%hour=hh
                cd%seconds=min*60.0 + sc
              END FUNCTION make_current_date



	SUBROUTINE dup_timestamp (ts1,ts2)
		TYPE(timestamp),INTENT(IN)   :: ts1
		TYPE(timestamp),INTENT(OUT)  :: ts2
		ts2%jdate = ts1%jdate
		ts2%secs  = ts1%secs
	END SUBROUTINE dup_timestamp


	FUNCTION tdif_secs (ts1, ts2) RESULT (dif)
		TYPE(timestamp),INTENT(IN)   :: ts1, ts2
		REAL                         :: dif
		dif = spd*REAL(ts2%jdate - ts1%jdate) + ts2%secs - ts1%secs
	END FUNCTION tdif_secs


	FUNCTION ts_earlier (ts1, ts2) RESULT (ts_first)
		TYPE(timestamp),INTENT(IN)   :: ts1, ts2
		TYPE(timestamp)				 :: ts_first
		IF (tdif_secs (ts1, ts2) > 0.0) THEN
			ts_first = ts1
		ELSE
			ts_first = ts2
		END IF
	END FUNCTION ts_earlier


	FUNCTION ts_later (ts1, ts2) RESULT (ts_last)
		TYPE(timestamp),INTENT(IN)	:: ts1, ts2
		TYPE(timestamp)		 	:: ts_last
		IF (tdif_secs (ts1, ts2) > 0.0) THEN
			ts_last = ts2
		ELSE
			ts_last = ts1
		END IF
	END FUNCTION ts_later

	!	returns: -1 == earlier, 0 == contained, +1 == later
	FUNCTION InInterval (t1, t2, t3) RESULT(In)
		TYPE(timestamp), INTENT(IN)	:: t1, t2, t3
                INTEGER :: In
		In = -1
		IF (tdif_secs (t1, t2) >= 0.0) THEN
			IF (tdif_secs (t2, t3) >= 0.0) THEN
				In = 0
			ELSE
				In = 1
			END IF	 	
		END IF
	END FUNCTION InInterval
  

	SUBROUTINE add_secs (ts, seconds, fixit)
		TYPE(timestamp), INTENT(INOUT) :: ts
		REAL, INTENT(IN)               :: seconds
		LOGICAL, INTENT(IN), OPTIONAL  :: fixit
		INTEGER	  		       :: hour, minute, sec

		ts%secs  = ts%secs + seconds
		IF (seconds >= 0) THEN
			DO     
                                IF  (ts%secs >= spd) THEN
				     ts%jdate = ts%jdate + 1
				     ts%secs  = ts%secs - spd
                                ENDIF
                                if ( ts%secs < spd) exit  
			END DO
		ELSE
			DO     
				ts%jdate = ts%jdate - 1
				ts%secs  = ts%secs + spd
                                if ( ts%secs > 0) exit    
			END DO
		END IF

		!	adjust to nearest half-hour
		IF (PRESENT(fixit)) THEN
			CALL get_hms (ts%secs, hour, minute, sec)
			SELECT CASE (minute)
			CASE(:29)
				minute = 0
			CASE(30:)
				minute = 30
			END SELECT
			ts%secs  = sph*REAL(hour) + spm*REAL(minute)
		END IF

	END SUBROUTINE add_secs 


	SUBROUTINE add_month (jdate, force_day)
		INTEGER,INTENT(INOUT)			:: jdate
		INTEGER, INTENT(IN), OPTIONAL	:: force_day
		INTEGER				:: year, month, day
		CALL get_ymd (jdate, year, month, day)
		IF (PRESENT(force_day)) day = force_day
		month = month + 1
		IF (month > 12) THEN
			year = year + 1
			month = 1
		END IF
		jdate = julian_date (year,month,MIN(MAX(day, 1),max_day(month, year)))
	END SUBROUTINE add_month
  

	FUNCTION julian_date (yyyy, mm, dd) RESULT (julian)
		!          converts calendar date to Julian date
		!          cf Fliegel & Van Flandern, CACM 11(10):657, 1968
		!          example: julian_date(1970,1,1)=2440588
		INTEGER,INTENT(IN)           :: yyyy,mm,dd
		INTEGER                      :: julian
		julian = dd - 32075 + 1461*(yyyy + 4800 + (mm - 14)/12)/4 +  &
				 367*(mm - 2 - ((mm - 14)/12)*12)/12 -       &
				 3*((yyyy + 4900 + (mm - 14)/12)/100)/4
	END FUNCTION julian_date
  

	SUBROUTINE get_ymd (jd, yyyy, mm, dd)
		!          expands a Julian date into a calendar date
		!          cf Fliegel & Van Flandern, CACM 11(10):657, 1968
		INTEGER,INTENT(IN)           :: jd
		INTEGER,INTENT(OUT)          :: yyyy,mm,dd
		INTEGER                      :: l,n
		l		= jd + 68569
		n		= 4*l/146097
		l		= l - (146097*n + 3)/4
		yyyy	= 4000*(l + 1)/1461001
		l		= l - 1461*yyyy/4 + 31
		mm		= 80*l/2447
		dd		= l - 2447*mm/80
		l		= mm/11
		mm		= mm + 2 - 12*l
		yyyy	= 100*(n - 49) + yyyy + l
	END SUBROUTINE get_ymd


  
	FUNCTION day_of_week (yyyy,mm,dd) RESULT (dow)
		!          Day_Of_Week: (0=Sunday,1=Monday...6=Saturday)
		!          cf J.D.Robertson, CACM 15(10):918
                !          renamed dow->day_of_week, keep dow as internal,  DSF
		INTEGER,INTENT(IN)           :: yyyy,mm,dd
                INTEGER :: dow
		dow = MODULO((13*(mm+10-(mm+10)/13*12)-1)/5+dd+77           &
			  +5*(yyyy+(mm-14)/12-(yyyy+(mm-14)/12)/100*100)/4   &
			  +(yyyy+(mm-14)/12)/400-(yyyy+(mm-14)/12)/100*2,7)
	END FUNCTION day_of_week

              
	FUNCTION day_of_year (yyyy,mm,dd) result (ndiy)
		!          day count in year
		!          cf J.D.Robertson, CACM 15(10):918
                !          renamed ndiy->day_of_year, keep ndiy as internal,  DSF
		INTEGER,INTENT(IN)           :: yyyy,mm,dd
                INTEGER :: ndiy   
		ndiy = 3055*(mm+2)/100-(mm+10)/13*2-91               &
			   +(1-(MODULO(yyyy,4)+3)/4+(MODULO(yyyy,100)+99)/100  &
			   -(MODULO(yyyy,400)+399)/400)*(mm+10)/13+dd 
	END FUNCTION day_of_year


	FUNCTION max_day (month,year) RESULT (maxd)
		INTEGER,INTENT(IN)                 :: month,year
		INTEGER                            :: maxd
		INTEGER,DIMENSION(12),PARAMETER    :: daycount =  & 
		(/31,28,31,30,31,30,31,31,30,31,30,31/)
		!          table lookup for most months
		maxd = daycount(month)
		!          correct February in a leap year
		IF (month == 2) THEN
			IF (leapyear(year)) maxd = maxd + 1
		END IF
	END FUNCTION max_day
      

	FUNCTION leapyear (year) result (leap)   
		INTEGER,INTENT(IN)		:: year
                logical                         :: leap   
		IF (day_of_year(year, 12, 31) > 365) THEN
			leap = .TRUE.
		ELSE
			leap = .FALSE.
		END IF
	END FUNCTION leapyear


	FUNCTION y2dig (year) result(y2)
		INTEGER,INTENT(IN)     :: year
                integer                :: y2  
		SELECT CASE (year)
		CASE (1900:1999) 
			y2 = year - 1900
		CASE (2000:2099)
			y2 = year - 2000
		CASE DEFAULT
			y2 = 0
		END SELECT
	END FUNCTION y2dig

   
	FUNCTION y4dig (year) result(y4)
		INTEGER,INTENT(IN)     :: year
                integer                :: y4  
		SELECT CASE (year)
		CASE (:90)
			y4 = year + 2000
		CASE (91:99)
			y4 = year + 1900
		CASE (1990:)
			y4 = year
		END SELECT
	END FUNCTION y4dig


	SUBROUTINE get_hms (secs,hour,minute,second)
		REAL,INTENT(IN)            :: secs
		INTEGER,INTENT(OUT)        :: hour,minute,second
		hour   = INT(secs/sph)
		minute = INT((secs - sph*REAL(hour))/spm)
		second = INT(secs - sph*REAL(hour) - spm*REAL(minute))
	END SUBROUTINE get_hms

	SUBROUTINE Init_nmdays (indate)
		TYPE(date),INTENT(IN)              :: indate
		INTEGER                            :: month,maxd,yy
		INTEGER,DIMENSION(12),PARAMETER    :: daycount =  & 
		(/31,28,31,30,31,30,31,31,30,31,30,31/)
		!          table lookup for most months

                yy = indate%year
                nydays=0
                do month=1,12
                   maxd = daycount(month)
                   !          correct February in a leap year
                   IF (month == 2) THEN
                      IF (leapyear(yy)) maxd = maxd + 1
                   END IF
                   nmdays(month)=maxd
                   nydays=nydays+maxd
                enddo
        END SUBROUTINE Init_nmdays

END MODULE TimeDate_ml
