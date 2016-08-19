! <TimeDate_ExtraUtil_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
MODULE TimeDate_ExtraUtil_ml

use ModelConstants_ml,only: METSTEP, MasterProc, IOU_MON,IOU_DAY,IOU_HOUR_MEAN
use My_Outputs_ml,    only: FREQ_HOURLY
use TimeDate_ml,      only: max_day,tdif_secs,tdif_days,add_secs,add_days,&
                            ts2date=>make_current_date,date2ts=>make_timestamp,&
                            date,timestamp,startdate,enddate
use CheckStop_ml,     only: CheckStop

IMPLICIT NONE
PRIVATE

public :: &
  assign_NTERM,     & ! set NTERM, the number of 3-hourly periods
  date2string,      & ! date (various formats) --> formatted string
  idate2nctime,     & ! idate (int array)--> secs since(int)/days since(real)
  nctime2idate,     & ! idate2nctime inverse
  nctime2string       ! as date2string, but from  secs/days since...

interface date2string
  module procedure detail2str,cd2str,int2str,ts2str
end interface date2string

interface idate2nctime
  module procedure idate_to_secs1970,idate_to_days1900
end interface idate2nctime

interface nctime2idate
  module procedure secs1970_to_idate,days1900_to_idate
end interface nctime2idate

interface nctime2string
  module procedure secs2str,days2str
end interface nctime2string

character(len=*), public, parameter :: &  ! Keywords for ncdate2string
! secs1970_key="secs  19.70"           ,&
! days1900_key="days  19.00"           ,&
  nctime_key  ="12345678.12"

private ::  &
  key2str,          & ! basic keyword substitution
  ikey2str,rkey2str,& ! auxiliary keyword substitution tools
  detail2str,             & ! detailed date input--> formatted string
  cd2str,int2str,ts2str,  & ! date/idate (int array)/timestamp --> formatted string
  idate_to_secs1970,& ! idate (int array)--> secs since 1970-01-01 00:00 UTC (int)
  idate_to_days1900,& ! idate (int array)--> days since 1900-01-01 00:00 UTC (real)
  secs1970_to_idate,& ! idate_to_secs1970 inverse
  days1900_to_idate,& ! idate_to_days1900 inverse
  to_stamp,         & ! extended interface for make_timestamp (TimeDate_ml)
  to_date,          & ! extended interface for make_current_date (TimeDate_ml)
  to_idate,         & ! create int array from timestap or date
  int2ts,int2date,ts2int,date2int,& ! auxiliary dateformat transformation tools
  init_ts             ! init 1900 & 1970 timestamps

interface key2str
  module procedure ikey2str,rkey2str
end interface key2str

interface to_stamp
  module procedure date2ts,int2ts
end interface to_stamp

interface to_date
  module procedure ts2date,int2date
end interface to_date

interface to_idate
  module procedure ts2int,date2int
end interface to_idate

real, private, parameter        :: spd=86400.0,sph=3600.0,spm=60.0
type(timestamp), private, save  :: ts1900,ts1970
logical, private, save          :: first_call=.true.

CONTAINS

function int2date (id) result (cd)
  type(date)                        :: cd
  integer, intent(in), dimension(:) :: id
  select case (size(id))
    case (3);     cd=date(id(1),id(2),id(3),0,0)
    case (4);     cd=date(id(1),id(2),id(3),id(4),0)
    case (5);     cd=date(id(1),id(2),id(3),id(4),id(5))
    case default; cd=date(-1,-1,-1,-1,-1)
    call CheckStop("ERROR in int2date: undetermined date")
  end select
end function int2date

function date2int (cd,n) result (id)
  type(date), intent(in)            :: cd
  integer, intent(in)               :: n
  integer, dimension(n)             :: id
  select case (n)
    case (1);     id=cd%year*1000000+cd%month*10000+cd%day*100+cd%hour ! ad-hoc convention
    case (3);     id=(/cd%year,cd%month,cd%day/)
    case (4);     id=(/cd%year,cd%month,cd%day,cd%hour/)
    case (5);     id=(/cd%year,cd%month,cd%day,cd%hour,cd%seconds/)
    case default; id(:)=-1
    call CheckStop("ERROR in date2int: undetermined date")
  end select
end function date2int

function int2ts (id) result (ts)
  type(timestamp)                   :: ts
  integer, intent(in), dimension(:) :: id
  ts=date2ts(int2date(id))
end function int2ts

function ts2int (ts,n) result (id)
  type(timestamp), intent(in)       :: ts
  integer, intent(in)               :: n
  integer, dimension(n)             :: id
  id=date2int(ts2date(ts),n)
end function ts2int

subroutine init_ts()
  if(.not.first_call)return
  first_call=.false.
  ts1970=to_stamp(date(1970,1,1,0,0))
  ts1900=to_stamp(date(1900,1,1,0,0))
end subroutine init_ts

function ikey2str(iname,key,val) result(fname)
  implicit none
  character(len=*), intent(in) :: iname,key
  character(len=len(iname))    :: fname
  integer, intent(in)          :: val
  character(len=9)             :: ifmt="(I??.??)"
  integer :: ind=0,n=0
  fname=iname
  ind=index(fname,trim(key))
  if(ind==0)return
  n=len_trim(key)
  write(ifmt,"('(I',I0,'.',I0,')')")n,n
  do while (ind>0)
    write(fname(ind:ind+n-1),ifmt)val
    ind=index(fname,trim(key))
  enddo
end function ikey2str

function rkey2str(iname,key,val) result(fname)
  implicit none
  character(len=*), intent(in) :: iname,key
  character(len=len(iname))    :: fname
  real, intent(in)             :: val
  character(len=9)             :: ifmt="(F??.??)"
  integer :: ind=0,n=0,n1=0
  fname=iname
  ind=index(fname,trim(key))
  if(ind==0)return
  n=len_trim(key)
  n1=index(key,".")
  if(n1==0)n1=n
  write(ifmt,"('(F',I0,'.',I0,')')")n1-1,n-n1
  do while (ind>0)
    write(fname(ind:ind+n-1),ifmt)val
    ind=index(fname,trim(key))
  enddo
end function rkey2str

function detail2str(iname,year,month,day,hour,seconds,minute,second,&
                        fstep,ntme,nlev,nlat,nlon,debug) result(fname)
  implicit none
  character(len=*), intent(in) :: iname
  character(len=len(iname))    :: fname
  integer, intent(in),optional :: year,month,day,hour,seconds,minute,second,&
                                  fstep,ntme,nlev,nlat,nlon
  logical, intent(in),optional :: debug
  fname=iname
  if(present(seconds))fname=key2str(fname,'ssss',seconds)
  if(present(second ))fname=key2str(fname,'ss'  ,second)
  if(present(minute ))fname=key2str(fname,'mm'  ,minute)
  if(present(hour   ))fname=key2str(fname,'hh'  ,hour  )
  if(present(day    ))fname=key2str(fname,'DD'  ,day   )
  if(present(month  ))fname=key2str(fname,'MM'  ,month )
  if(present(year   ))fname=key2str(fname,'YYYY',year  )
  if(present(year   ))fname=key2str(fname,'YY'  ,mod(year,100))
  if(present(nlon   ))fname=key2str(fname,'LON' ,nlon  )
  if(present(nlat   ))fname=key2str(fname,'LAT' ,nlat  )
  if(present(nlev   ))fname=key2str(fname,'LL'  ,nlev  )
  if(present(ntme   ))fname=key2str(fname,'TTT' ,ntme  )
  if(present(fstep  ))fname=key2str(fname,'FFF' ,fstep )
  if(present(debug))then
    if(debug) print *,'date2string: ',trim(iname),'-->',trim(fname)
  endif
end function detail2str

function cd2str(iname,cd,debug) result(fname)
  implicit none
  character(len=*), intent(in)   :: iname
  character(len=len(iname))      :: fname
  type(date),intent(in)          :: cd
  logical, intent(in),  optional :: debug
  fname=detail2str(iname,year=cd%year,month=cd%month,day=cd%day,&
                         hour=cd%hour,seconds=cd%seconds,&
                         minute=cd%seconds/60,second=mod(cd%seconds,60),&
                         debug=debug)
end function cd2str

function int2str(iname,id,debug) result(fname)
  implicit none
  character(len=*), intent(in)      :: iname
  character(len=len(iname))         :: fname
  integer, intent(in), dimension(:) :: id
  logical, intent(in),  optional    :: debug
  fname=cd2str(iname,to_date(id),debug=debug)
end function int2str

function ts2str(iname,ts,debug) result(fname)
  implicit none
  character(len=*), intent(in)      :: iname
  character(len=len(iname))         :: fname
  type(timestamp), intent(in)       :: ts
  logical, intent(in),  optional    :: debug
  fname=cd2str(iname,to_date(ts),debug=debug)
end function ts2str

subroutine idate_to_secs1970(idate,nsecs,iotyp)
!calculate how many seconds have passed since the start of the year 1970
  integer, intent(in), dimension(:) :: idate
  integer, intent(out)              :: nsecs
  integer, optional, intent(in)     :: iotyp

  if(first_call)call init_ts()
  nsecs=tdif_secs(ts1970,to_stamp(idate))
  call CheckStop(nsecs<0,"ERROR in idate_to_secs1970: date previous to 1970-01-01")

  if(present(iotyp))then
    select case (iotyp)  !middle of period: NB WORKS ONLY FOR COMPLETE PERIODS
      case (IOU_MON      ); nsecs=nsecs-spd/2*max_day(idate(2),idate(1)) !#days(jan)=#days(dec)
      case (IOU_DAY      ); nsecs=nsecs-spd/2
      case (IOU_HOUR_MEAN); nsecs=nsecs-sph/2*FREQ_HOURLY
    end select
  endif
end subroutine idate_to_secs1970

subroutine idate_to_days1900(idate,ndays,iotyp)
! calculate how many days have passed since the start of the year 1900
  integer, intent(in), dimension(:) :: idate
  real(kind=8), intent(out)         :: ndays
  integer, optional, intent(in)     :: iotyp

  if(first_call)call init_ts()
  ndays=tdif_days(ts1900,to_stamp(idate))
  call CheckStop(ndays<0,"ERROR in idate_to_days1900: date previous to 1900-01-01")

  if(present(iotyp))then
    select case (iotyp)  !middle of period: NB WORKS ONLY FOR COMPLETE PERIODS
      case (IOU_MON      ); ndays=ndays-0.5*max_day(idate(2),idate(1))
      case (IOU_DAY      ); ndays=ndays-0.5
      case (IOU_HOUR_MEAN); ndays=ndays-FREQ_HOURLY/48.0  !1.0/48.0=half hour
    end select
  endif
end subroutine idate_to_days1900

subroutine secs1970_to_idate(idate,nsecs,msg)
!calculate date from seconds that have passed since the start of the year 1970
  integer, intent(out), dimension(:)      :: idate
  integer, intent(in)                     :: nsecs
  character(len=*), intent(in), optional  :: msg
  type(timestamp) :: ts

  if(first_call)call init_ts()
  ts=ts1970
  call add_days(ts,nsecs/spd)
  idate=to_idate(ts,size(idate))

  if(present(msg)) print *,date2string(msg,idate)
end subroutine secs1970_to_idate

subroutine days1900_to_idate(idate,ndays,msg)
!calculate date from seconds that have passed since the start of the year 1900
  integer, intent(out), dimension(:)      :: idate
  real(kind=8), intent(in)                :: ndays
  character(len=*), intent(in), optional  :: msg
  type(timestamp) :: ts

  if(first_call)call init_ts()
  ts=ts1900
  call add_days(ts,ndays)
  idate=to_idate(ts,size(idate))

  if(present(msg)) print *,date2string(msg,idate)
end subroutine days1900_to_idate

function secs2str(iname,nsecs,debug) result(fname)
  implicit none
  character(len=*), intent(in)            :: iname
  character(len=len(iname))               :: fname
  integer, intent(in)                     :: nsecs
  integer, dimension(5)                   :: idate
  logical, intent(in), optional           :: debug
  call nctime2idate(idate,nsecs)
  fname=date2string(key2str(iname,nctime_key,nsecs),idate,debug=debug)
end function secs2str

function days2str(iname,ndays,debug) result(fname)
  implicit none
  character(len=*), intent(in)            :: iname
  character(len=len(iname))               :: fname
  real(kind=8), intent(in)                :: ndays
  integer, dimension(5)                   :: idate
  logical, intent(in),  optional          :: debug
  call nctime2idate(idate,ndays)
  fname=date2string(key2str(iname,nctime_key,ndays),idate,debug=debug)
end function days2str

subroutine assign_NTERM(NTERM)
! calculate NTERM (the number of metdata periods)
! on the basis of start and enddate
  integer,intent(out) :: NTERM
  type(timestamp)     :: ts1, ts2
  real, parameter     :: spMETSTEP=sph*METSTEP  ! seconds in period of metadat

  ! ensure that a valid day of the month,
  ! e.g. Feb 31=>Feb 28/29 depending the year
  startdate(3)=min(startdate(3),max_day(startdate(2),startdate(1)))
  enddate  (3)=min(enddate  (3),max_day(enddate  (2),enddate  (1)))

  startdate(4)=0  ! simulation starts at 00:00 UTC
  enddate  (4)=24 ! simulation ends   at 24:00 UTC

  ts1=to_stamp(startdate)
  ts2=to_stamp(enddate)

  NTERM=1+ceiling(tdif_secs(ts1,ts2)/spMETSTEP) !NTERM=1+#time-step
  if(NTERM<=1)then
    if(MasterProc)then
      write(*,*)'WARNING: enddate before startdate, running only one metstep'
      write(*,*)'Start date: ',startdate
      write(*,*)'End   date: ',enddate
    endif
    NTERM=max(2,NTERM)!run at least one period
  endif
end subroutine assign_NTERM

END MODULE TimeDate_ExtraUtil_ml
