! <TimeDate_ExtraUtil_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2019 met.no
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
MODULE TimeDate_ExtraUtil_mod

use Par_mod,           only: me
use Config_module,only: METSTEP, MasterProc, startdate,enddate , &
                            IOU_MON,IOU_DAY,IOU_HOUR,FREQ_HOURLY
use SmallUtils_mod,    only: key2str,to_upper
use CheckStop_mod,     only: CheckStop
use TimeDate_mod,      only: max_day,tdif_secs,tdif_days,&
                            date_addSecs=>add2current_date,ts_addSecs=>add_secs,&
                            ts2date=>make_current_date,date2ts=>make_timestamp,&
                            timestamp,day_of_year,date,make_timestamp,current_date

IMPLICIT NONE
PRIVATE

public :: &
  date2string,      & ! date (various formats) --> formatted string
  date2file,        & ! date (various formats) --> file name w/option for old file version
  string2date,      & ! formatted string & format (pattern) --> date (CD format)
  date2nctime,      & ! date (various formats)--> secs since(int)/days since(real)
  nctime2date,      & ! date2nctime inverse
  nctime2string,    & ! as date2string, but from  secs/days since...
  compare_date,     & ! compare dates, accepts wildcards
  to_stamp,         & ! extended interface for make_timestamp (TimeDate_mod)
  to_date,          & ! extended interface for make_current_date (TimeDate_mod)
  to_idate,         & ! create int array from timestap or date
  self_test,        & ! test/example ussages for different interface/module procedures
  assign_startandenddate,& !if needed, correct days and hours of startdate and enddate
  date_is_reached!test if the first date is before or equal to the second
interface date2string!(iname,...,mode,debug) result(fname)
!   character(len=*),intent(in) :: iname
!   character(len=len(iname))   :: fname
!   character(len=*),intent(in),optional :: mode
!   logical,intent(in),optional :: debug
  module procedure detail2str!(iname,year,month,day,hour,seconds,minute,second,days,&
!                              fstep,ntme,nlev,nlat,nlon,nproc,mode,debug)
!   integer,intent(in),optional :: year,month,day,hour,seconds,&
!                                  minute,second,days,&
!                                  fstep,ntme,nlev,nlat,nlon,nproc
  module procedure cd2str !(iname,cd,addsecs,mode,debug)
!   type(date),intent(in)       :: cd
!   real,intent(in),optional    :: addsecs
  module procedure ts2str !(iname,ts,addsecs,mode,debug)
!   type(timestamp),intent(in)  :: ts
!   real,intent(in),optional    :: addsecs
  module procedure int2str!(iname,id,addsecs,mode,debug)
!   integer,intent(in)          :: id(:)
!   real,intent(in),optional    :: addsecs
end interface date2string

interface date2file!(iname,...,max_age,age_unit,mode,last,debug) result(fname)
!   character(len=*),intent(in) :: iname,age_unit
!   character(len=len(iname))   :: fname
!   integer,intent(in)          :: max_age
!   character(len=*),intent(in),optional :: mode
!   integer,intent(in),optional :: last
!   logical,intent(in),optional :: debug
  module procedure cd2file !(iname,cd,max_age,age_unit,mode,last,debug)
!   type(date),intent(in)       :: cd
  module procedure ts2file !(iname,ts,max_age,age_unit,mode,last,debug)
!   type(timestamp),intent(in)  :: ts
  module procedure int2file!(iname,id,max_age,age_unit,mode,last,debug)
!   integer, intent(in)         :: id(:)
end interface date2file

interface date2nctime
  module procedure ts_to_secs1970
  module procedure cd_to_secs1970
  module procedure int_to_secs1970
  module procedure ts_to_days1900
  module procedure cd_to_days1900
  module procedure int_to_days1900
end interface date2nctime

interface nctime2date
  module procedure secs1970_to_ts
  module procedure secs1970_to_cd
  module procedure secs1970_to_int
  module procedure days1900_to_ts
  module procedure days1900_to_cd
  module procedure days1900_to_int
end interface nctime2date

interface nctime2string
  module procedure secs2str
  module procedure days2str
end interface nctime2string

interface to_stamp
  module procedure date2ts
  module procedure int2ts
  module procedure str2ts
end interface to_stamp

interface to_date
  module procedure ts2date
  module procedure int2date
  module procedure string2date
end interface to_date

interface to_idate
  module procedure ts2int
  module procedure date2int
  module procedure str2int
end interface to_idate

character(len=*), public, parameter :: &  ! Keywords for ncdate2string
! secs1970_key="secs  19.70",&
! days1900_key="days  19.00",&
! nctime_key  ="12345678.12"
  nctime_key  ='#NCTIMEKEY#',nctime_fmt='(F11.2)'

private ::  &
  detail2str,str2detail, & ! detailed date input<--> formatted string
  cd2str,int2str,ts2str, & ! date/idate (int array)/timestamp --> formatted string
  cd2file,int2file,ts2file,&  ! * --> file name w/option for old file version
  ts_to_secs1970,cd_to_secs1970,int_to_secs1970,& ! * --> secs since 1970-01-01 00:00 UTC (int)
  ts_to_days1900,cd_to_days1900,int_to_days1900,& ! * --> days since 1900-01-01 00:00 UTC (real)
  secs1970_to_ts,secs1970_to_cd,secs1970_to_int,& ! *_to_secs1970 inverse
  days1900_to_ts,days1900_to_cd,days1900_to_int,& ! *_to_days1900 inverse
  int2ts,int2date,ts2int,date2int,str2ts,str2int,& ! auxiliary dateformat transformation tools
  init_ts             ! init 1900 & 1970 timestamps

real, private, parameter        :: spd=86400.0,sph=3600.0,spm=60.0
type(timestamp), private, save  :: ts1900,ts1970
logical, private, save          :: first_call=.true.

CONTAINS

function int2date (id) result (cd)
  type(date)                        :: cd
  integer, intent(in), dimension(:) :: id
  select case (size(id))
    case (1);! ad-hoc convention
    cd=date(id(1)/1000000,mod(id(1)/10000,100),mod(id(1)/100,100),mod(id(1),100),0)
    case (3);     cd=date(id(1),id(2),id(3),0,0)
    case (4);     cd=date(id(1),id(2),id(3),id(4),0)
    case (5);     cd=date(id(1),id(2),id(3),id(4),id(5))
    case (6);     cd=date(id(1),id(2),id(3),id(4),id(5)*60+id(6))
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
    case (3);     id=[cd%year,cd%month,cd%day]
    case (4);     id=[cd%year,cd%month,cd%day,cd%hour]
    case (5);     id=[cd%year,cd%month,cd%day,cd%hour,cd%seconds]
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

subroutine str2detail(str,fmt,year,month,day,hour,seconds,minute,second,days,&
                     fstep,ntme,nlev,nlat,nlon,debug)
  character(len=*), intent(in) :: str,fmt
  integer,intent(out),optional :: year,month,day,hour,seconds,&
                                  minute,second,days,&
                                  fstep,ntme,nlev,nlat,nlon
  logical, intent(in),optional :: debug
  if(present(seconds))seconds=str2key(str,fmt,'ssss')&
                             +str2key(str,fmt,'ss'  )&
                             +str2key(str,fmt,'mm'  )*60
  if(present(second ))second =str2key(str,fmt,'ss'  )
  if(present(minute ))minute =str2key(str,fmt,'mm'  )
  if(present(hour   ))hour   =str2key(str,fmt,'hh'  )
  if(present(day    ))day    =str2key(str,fmt,'DD'  )
  if(present(month  ))month  =str2key(str,fmt,'MM'  )
  if(present(year   ))year   =str2key(str,fmt,'YYYY')
! if(present(year   ))year   =str2key(str,fmt,'YY'  )+1900
  if(present(days   ))days   =str2key(str,fmt,'JJJ' ) ! day of the year
  if(present(nlon   ))nlon   =str2key(str,fmt,'LON' )
  if(present(nlat   ))nlat   =str2key(str,fmt,'LAT' )
  if(present(nlev   ))nlev   =str2key(str,fmt,'LL'  )
  if(present(ntme   ))ntme   =str2key(str,fmt,'TTT' )
! if(present(fstep  ))fstep  =str2key(str,fmt,'FFFF')
  if(present(fstep  ))fstep  =str2key(str,fmt,'FFF' )
! if(present(nproc  ))nproc  =str2key(str,fmt,'PPPP')
! if(present(nproc  ))nproc  =str2key(str,fmt,'PPP' )
  if(present(debug))then
    if(debug) write(*,*)'string2date: ',trim(str),'/',trim(fmt)
  end if
contains
function str2key(str,xfmt,key) result(val)
  character(len=*), intent(in) :: str,xfmt,key
  integer :: val
  integer :: ind=0
  val=0
  ind=index(xfmt,trim(key))
  if(ind>0)read(str(ind:ind+len_trim(key)-1),*)val
end function str2key
end subroutine str2detail

function string2date(str,fmt,debug) result(cd)
  character(len=*), intent(in)   :: str,fmt
  logical, intent(in),  optional :: debug
  type(date)                     :: cd
  call str2detail(str,fmt,year=cd%year,month=cd%month,day=cd%day,&
                  hour=cd%hour,seconds=cd%seconds,debug=debug)
end function string2date

function str2ts(str,fmt,debug) result(ts)
  character(len=*), intent(in)   :: str,fmt
  logical, intent(in),  optional :: debug
  type(timestamp)                :: ts
  ts=date2ts(string2date(str,fmt,debug=debug))
end function str2ts

function str2int(str,fmt,n,debug) result(id)
  character(len=*), intent(in)   :: str,fmt
  integer, intent(in)               :: n
  logical, intent(in),  optional :: debug
  integer, dimension(n)          :: id
  id=date2int(string2date(str,fmt,debug=debug),n)
end function str2int

function detail2str(iname,year,month,day,hour,seconds,minute,second,days,&
                    fstep,ntme,nlev,nlat,nlon,nproc,mode,debug) result(fname)
  character(len=*),intent(in) :: iname
  character(len=len(iname))   :: fname
  integer,intent(in),optional :: year,month,day,hour,seconds,&
                                 minute,second,days,&
                                 fstep,ntme,nlev,nlat,nlon,nproc
  character(len=*),intent(in),optional :: mode
  character(len=len('YMDHMS'))         :: my_mode
  logical,intent(in),optional :: debug

  fname=iname
  my_mode='ALL';if(present(mode))my_mode=mode
  select case(to_upper(my_mode))
  case('ALL')
    if(present(seconds))fname=key2str(fname,'ssss',seconds)
    if(present(second ))fname=key2str(fname,'ss'  ,second)
    if(present(minute ))fname=key2str(fname,'mm'  ,minute)
    if(present(hour   ))fname=key2str(fname,'hh'  ,hour  )
    if(present(day    ))fname=key2str(fname,'DD'  ,day   )
    if(present(month  ))fname=key2str(fname,'MM'  ,month )
    if(present(year   ))fname=key2str(fname,'YYYY',year  )
    if(present(year   ))fname=key2str(fname,'YY'  ,mod(year,100))
    if(present(days   ))fname=key2str(fname,'JJJ' ,days  ) ! day of the year
    if(present(nlon   ))fname=key2str(fname,'LON' ,nlon  )
    if(present(nlat   ))fname=key2str(fname,'LAT' ,nlat  )
    if(present(nlev   ))fname=key2str(fname,'LL'  ,nlev  )
    if(present(ntme   ))fname=key2str(fname,'TTT' ,ntme  )
    if(present(fstep  ))fname=key2str(fname,'FFFF',fstep )
    if(present(fstep  ))fname=key2str(fname,'FFF' ,fstep )
    if(present(nproc  ))fname=key2str(fname,'PPPP',nproc )
    if(present(nproc  ))fname=key2str(fname,'PPP' ,nproc )
  case('YMDHMS')
    if(present(seconds))fname=key2str(fname,'ssss',seconds)
    if(present(second ))fname=key2str(fname,'ss'  ,second)
    if(present(minute ))fname=key2str(fname,'mm'  ,minute)
    if(present(hour   ))fname=key2str(fname,'hh'  ,hour  )
    if(present(day    ))fname=key2str(fname,'DD'  ,day   )
    if(present(month  ))fname=key2str(fname,'MM'  ,month )
    if(present(year   ))fname=key2str(fname,'YYYY',year  )
    if(present(year   ))fname=key2str(fname,'YY'  ,mod(year,100))
  case('YMDH')
    if(present(hour   ))fname=key2str(fname,'hh'  ,hour  )
    if(present(day    ))fname=key2str(fname,'DD'  ,day   )
    if(present(month  ))fname=key2str(fname,'MM'  ,month )
    if(present(year   ))fname=key2str(fname,'YYYY',year  )
    if(present(year   ))fname=key2str(fname,'YY'  ,mod(year,100))
  case default
    call CheckStop("Unsupported date2string(mode='"//trim(my_mode)//"')")
  end select
  if(present(debug))then
    if(debug) write(*,*)'date2string: ',trim(iname),'-->',trim(fname)
  end if
end function detail2str

function cd2str(iname,cd,addsecs,mode,debug) result(fname)
  character(len=*),intent(in) :: iname
  character(len=len(iname))   :: fname
  type(date),intent(in)       :: cd
  real,intent(in),optional    :: addsecs
  character(len=*),intent(in),optional :: mode
  logical,intent(in),optional :: debug

  type(date) :: ccd
  ccd=cd
  if(present(addsecs))call date_addSecs(ccd,addsecs) 
  fname=detail2str(iname,year=ccd%year,month=ccd%month,day=ccd%day,&
                         hour=ccd%hour,seconds=ccd%seconds,&
                         minute=ccd%seconds/60,second=mod(ccd%seconds,60),&
                         days=day_of_year(ccd%year,ccd%month,ccd%day),&
                         fstep=nint(tdif_days(to_stamp(startdate),to_stamp(ccd))*24),&
                         nproc=me,mode=mode,debug=debug)
end function cd2str

function ts2str(iname,ts,addsecs,mode,debug) result(fname)
  character(len=*),intent(in) :: iname
  character(len=len(iname))   :: fname
  type(timestamp),intent(in)  :: ts
  real,intent(in),optional    :: addsecs
  character(len=*),intent(in),optional :: mode
  logical,intent(in),optional :: debug
  type(timestamp)       :: tts
  tts=ts
  if(present(addsecs))call ts_addSecs(tts,addsecs)
  fname=cd2str(iname,to_date(tts),mode=mode,debug=debug)
end function ts2str

function int2str(iname,id,addsecs,mode,debug) result(fname)
  character(len=*),intent(in) :: iname
  character(len=len(iname))   :: fname
  integer, intent(in)         :: id(:)
  real,intent(in),optional    :: addsecs
  character(len=*),intent(in),optional :: mode
  logical,intent(in),optional :: debug
  fname=ts2str(iname,to_stamp(id),addsecs=addsecs,mode=mode,debug=debug)
end function int2str

subroutine ts_to_secs1970(ts,nsecs,iotyp)
!calculate how many seconds have passed since the start of the year 1970
  type(timestamp), intent(in)       :: ts
  integer, intent(out)              :: nsecs
  integer, optional, intent(in)     :: iotyp
  type(date)                        :: cd
  integer, parameter ::   &
    half_day =INT(spd/2), & ! 12h in seconds
    half_hour=INT(sph/2)    ! 30m in seconds

  if(first_call)call init_ts()
  nsecs=INT(tdif_secs(ts1970,ts))
  call CheckStop(nsecs<0,"ERROR in date2nctime: date previous to 1970-01-01")

  if(present(iotyp))then
    select case (iotyp)  !middle of period: NB WORKS ONLY FOR COMPLETE PERIODS
    case(IOU_MON)
      cd=to_date(ts)
      nsecs=nsecs-half_day*max_day(cd%month,cd%year) !#days(jan)=#days(dec)
    case(IOU_DAY)
      nsecs=nsecs-half_day
    case(IOU_HOUR)
      nsecs=nsecs-half_hour*FREQ_HOURLY
    end select
  end if
end subroutine ts_to_secs1970

subroutine cd_to_secs1970(cd,nsecs,iotyp)
!calculate how many seconds have passed since the start of the year 1970
  type(date), intent(in)            :: cd
  integer, intent(out)              :: nsecs
  integer, optional, intent(in)     :: iotyp

  call ts_to_secs1970(to_stamp(cd),nsecs,iotyp=iotyp)
end subroutine cd_to_secs1970

subroutine int_to_secs1970(id,nsecs,iotyp)
!calculate how many seconds have passed since the start of the year 1970
  integer, intent(in), dimension(:) :: id
  integer, intent(out)              :: nsecs
  integer, optional, intent(in)     :: iotyp

  call ts_to_secs1970(to_stamp(id),nsecs,iotyp=iotyp)
end subroutine int_to_secs1970

subroutine ts_to_days1900(ts,ndays,iotyp)
! calculate how many days have passed since the start of the year 1900
  type(timestamp), intent(in)       :: ts
  real(kind=8), intent(out)         :: ndays
  integer, optional, intent(in)     :: iotyp
  type(date)                        :: cd

  if(first_call)call init_ts()
  ndays=tdif_days(ts1900,ts)
  call CheckStop(ndays<0,"ERROR in date2nctime: date previous to 1900-01-01")

  if(present(iotyp))then
    select case (iotyp)  !middle of period: NB WORKS ONLY FOR COMPLETE PERIODS
    case(IOU_MON)
      cd=to_date(ts)
      ndays=ndays-0.5*max_day(cd%month,cd%year) !#days(jan)=#days(dec)
    case(IOU_DAY)
      ndays=ndays-0.5
    case(IOU_HOUR)
      ndays=ndays-FREQ_HOURLY/48.0  !1.0/48.0=half hour
    end select
  end if
end subroutine ts_to_days1900

subroutine cd_to_days1900(cd,ndays,iotyp)
! calculate how many days have passed since the start of the year 1900
  type(date), intent(in)            :: cd
  real(kind=8), intent(out)         :: ndays
  integer, optional, intent(in)     :: iotyp

  call ts_to_days1900(to_stamp(cd),ndays,iotyp=iotyp)
end subroutine cd_to_days1900

subroutine int_to_days1900(id,ndays,iotyp)
! calculate how many days have passed since the start of the year 1900
  integer, intent(in), dimension(:) :: id
  real(kind=8), intent(out)         :: ndays
  integer, optional, intent(in)     :: iotyp

  call ts_to_days1900(to_stamp(id),ndays,iotyp=iotyp)
end subroutine int_to_days1900

subroutine secs1970_to_ts(ts,nsecs,msg)
!calculate date from seconds that have passed since the start of the year 1970
  type(timestamp), intent(out)            :: ts
  integer, intent(in)                     :: nsecs
  character(len=*), intent(in), optional  :: msg

  if(first_call)call init_ts()
  ts=ts1970
  call ts_addSecs(ts,float(nsecs))

  if(present(msg)) write(*,*)date2string(msg,ts)
end subroutine secs1970_to_ts

subroutine secs1970_to_cd(cd,nsecs,msg)
!calculate date from seconds that have passed since the start of the year 1970
  type(date), intent(out)                 :: cd
  integer, intent(in)                     :: nsecs
  character(len=*), intent(in), optional  :: msg
  type(timestamp) :: ts

  call secs1970_to_ts(ts,nsecs,msg=msg)
  cd=to_date(ts)
end subroutine secs1970_to_cd

subroutine secs1970_to_int(id,nsecs,msg)
!calculate date from seconds that have passed since the start of the year 1970
  integer, intent(out), dimension(:)      :: id
  integer, intent(in)                     :: nsecs
  character(len=*), intent(in), optional  :: msg
  type(timestamp) :: ts

  call secs1970_to_ts(ts,nsecs,msg=msg)
  ts%secs=nint(ts%secs)!to avoid 3599.9999 seconds
  id=to_idate(ts,size(id))
end subroutine secs1970_to_int

subroutine days1900_to_ts(ts,ndays,msg)
!calculate date from seconds that have passed since the start of the year 1900
  type(timestamp), intent(out)            :: ts
  real(kind=8), intent(in)                :: ndays
  character(len=*), intent(in), optional  :: msg

  if(first_call)call init_ts()
  ts=ts1900
  call ts_addSecs(ts,ndays*spd)

  if(present(msg)) write(*,*)date2string(msg,ts)
end subroutine days1900_to_ts

subroutine days1900_to_cd(cd,ndays,msg)
!calculate date from seconds that have passed since the start of the year 1900
  type(date), intent(out)                 :: cd
  real(kind=8), intent(in)                :: ndays
  character(len=*), intent(in), optional  :: msg
  type(timestamp) :: ts

  call days1900_to_ts(ts,ndays,msg=msg)
  cd=to_date(ts)
end subroutine days1900_to_cd

subroutine days1900_to_int(id,ndays,msg)
!calculate date from seconds that have passed since the start of the year 1900
  integer, intent(out), dimension(:)      :: id
  real(kind=8), intent(in)                :: ndays
  character(len=*), intent(in), optional  :: msg
  type(timestamp) :: ts

  call days1900_to_ts(ts,ndays,msg=msg)
  if(size(id)<=4)ts%secs=ts%secs+0.1!correct for rounding errors
  id=to_idate(ts,size(id))
end subroutine days1900_to_int

function secs2str(iname,nsecs,debug) result(fname)
  character(len=*), intent(in)            :: iname
  character(len=len(iname))               :: fname
  integer, intent(in)                     :: nsecs
  integer, dimension(5)                   :: idate
  logical, intent(in), optional           :: debug
  call nctime2date(idate,nsecs)
  fname=date2string(key2str(iname,nctime_key,nsecs,nctime_fmt),idate,debug=debug)
end function secs2str

function days2str(iname,ndays,debug) result(fname)
  character(len=*), intent(in)            :: iname
  character(len=len(iname))               :: fname
  real(kind=8), intent(in)                :: ndays
  integer, dimension(5)                   :: idate
  logical, intent(in),  optional          :: debug
  call nctime2date(idate,ndays)
  fname=date2string(key2str(iname,nctime_key,ndays,nctime_fmt),idate,debug=debug)
end function days2str

function compare_date(n,dateA,dateB,wildcard) result(equal)
  integer,   intent(in)           :: n
  type(date),intent(in)           :: dateA,dateB(n)
  integer,   intent(in), optional :: wildcard
  logical :: equal
  integer :: dA(5),dB(5),i
  equal=.false.
  dA=(/dateA%year,dateA%month,dateA%day,dateA%hour,dateA%seconds/)
  do i=1,n
    dB=(/dateB(i)%year,dateB(i)%month,dateB(i)%day,&
         dateB(i)%hour,dateB(i)%seconds/)
    if(present(wildcard))then
      equal=equal.or.all((dA==dB).or.(dA==wildcard).or.(dB==wildcard))
    else
      equal=equal.or.all(dA==dB)
    end if
  end do
end function compare_date

function ts2file(iname,ts,max_age,age_unit,mode,last,debug) result(fname)
  intent(in) :: iname,ts,max_age,age_unit,mode,last,debug
  optional   :: mode,last,debug
  type(timestamp) :: ts
  character(len=*):: iname,age_unit,mode
  character(len=len(iname)):: fname
  integer         :: max_age,last
  logical         :: debug
  integer :: age,ind,i
  real :: nsecs
  logical :: fexist

! nsecs: age_unit in seconds
  select case(age_unit)
  case('s','second','seconds','S','SECOND','SECONDS')
    nsecs=1e0
  case('h','hour','hours','H','HOUR','HOURS')
    nsecs=36e2
  case('d','day','days','D','DAY','DAYS')
    nsecs=864e2
  case default
    call CheckStop("Unsupported string2file(age_unit='"//trim(age_unit)//"')")
  end select

! find the nth '/' from the end of iname
  ind=0
  if(present(last))then
    ind=len_trim(iname)
    do i=1,last
      if(ind==0)exit
      ind=index(iname(:ind),'/',BACK=.true.)
    end do
  end if
! do not pharse the 1st ind chadacters
  if(ind>0)fname=iname(:ind)
  
  do age=0,max_age
    ! pharse from ind+1
    fname(ind+1:)=date2string(iname(ind+1:),ts,addsecs=-age*nsecs,&
                              mode=mode,debug=debug)
    inquire(file=fname,exist=fexist)
    if(fexist)exit
  end do
end function ts2file

function cd2file(iname,cd,max_age,age_unit,mode,last,debug) result(fname)
  intent(in) :: iname,cd,max_age,age_unit,mode,last,debug
  optional   :: mode,last,debug
  type(date)      :: cd
  character(len=*):: iname,age_unit,mode
  character(len=len(iname)):: fname
  integer         :: max_age,last
  logical         :: debug
  fname=ts2file(iname,to_stamp(cd),max_age,age_unit,&
                mode=mode,last=last,debug=debug)
end function cd2file

function int2file(iname,id,max_age,age_unit,mode,last,debug) result(fname)
  intent(in) :: iname,id,max_age,age_unit,mode,last,debug
  optional   :: mode,last,debug
  integer         :: id(:)
  character(len=*):: iname,age_unit,mode
  character(len=len(iname)):: fname
  integer         :: max_age,last
  logical         :: debug
  fname=ts2file(iname,to_stamp(id),max_age,age_unit,&
                mode=mode,last=last,debug=debug)
end function int2file

subroutine self_test()
  character(len=*),parameter :: &
    hfmt="(/I0,') Self-test - ',A,/32('='))", & ! header format
    tfmt="(A,/2X,A,', ',A,'.')",              & ! test format
    dfmt="YYYY-MM-DD hh:mm:ss",               & ! date format
    test_file="YYYYMMDD.test"
  integer, parameter :: &
    IO_TMP=27
  integer :: secs,intdate(1)
  real    :: days

  call init_ts()
  print hfmt,1,"date2string"
  print tfmt,'time stamp input',&
    date2string(dfmt,ts1900),&    ! timestamps calculated by init_ts
    date2string(dfmt,ts1970)
  print tfmt,'date structure input',&
    date2string(dfmt,date(1900,1,1,0,0)),&
    date2string(dfmt,date(1970,1,1,0,0))
  print tfmt,'int array (3..5 elements) input',&
    date2string(dfmt,[1900,1,1,0,0]),&
    date2string(dfmt,[1970,1,1,0,0])
  print tfmt,'integer ([YYYYMMDDhh]) input',&
    date2string(dfmt,[1900010100]),&
    date2string(dfmt,[1970010100])
  print tfmt,'add 10 days',&
    date2string(dfmt,ts1900)//&       ! spd=86400.0 (seconds per day)
      key2str("FFFFFFFFFF seconds","FFFFFFFFFF",10*spd,"(SP,F10.1)"),&
    date2string(dfmt,ts1900,10*spd)   ! +10 days in seconds
  print tfmt,'subtract 12 hour',&
    date2string(dfmt,ts1970)//&       ! sph=3600.0 (seconds per hour)
      key2str("FFFFFFFFFF seconds","FFFFFFFFFF",-12*sph,"(SP,F10.1)"),&
    date2string(dfmt,ts1970,-12*sph)  ! -12 hours in seconds

  print hfmt,2,"date2file"
  ! create test_file, 1 day before date
  open(IO_TMP,file=date2string(test_file,ts1970,-spd),status='replace')
  print tfmt,'serch for '//test_file//' up to 3 before date',&
    date2string(dfmt,ts1970),&
    date2file(test_file,ts1970,3,"days")  ! found file from 1 days before date
  ! delete test file
  close(IO_TMP,status='delete') 
  print tfmt,'serch for again',&
    date2string(dfmt,ts1970),&
    date2file(test_file,ts1970,3,"days")  ! no file found, return 3 days before date

  print hfmt,3,"nctime2string"
  days=1.25           ! real ==> days since 1900
  secs=INT(days*spd)  ! int  ==> secs since 1970
  print tfmt,'days since 1900 (real input)',&
    key2str("F.FFFF days since ","F.FFFF",days)//date2string(dfmt,ts1900),&
    nctime2string(dfmt,days)
  print tfmt,'secs since 1970 (integer input)',&
    key2str("HHHHHH secs since ","HHHHHH",secs)//date2string(dfmt,ts1970),&
    nctime2string(dfmt,secs)     

  print hfmt,4,"nctime2date/date2nctime"
  call nctime2date(intdate,days)
  print tfmt,'nctime2date',&
    key2str("F.FFFF days since ","F.FFFF",days)//date2string(dfmt,ts1900),&
    date2string(dfmt,intdate)              
  call date2nctime(intdate,days)
  print tfmt,'date2nctime',&
    key2str("F.FFFF days since ","F.FFFF",days)//date2string(dfmt,ts1900),&
    date2string(dfmt,intdate)              
  call nctime2date(intdate,secs)
  print tfmt,'nctime2date',&
    key2str("HHHHHH secs since ","HHHHHH",secs)//date2string(dfmt,ts1970),&
    date2string(dfmt,intdate)              
  call date2nctime(intdate,secs)
  print tfmt,'date2nctime',&
    key2str("HHHHHH secs since ","HHHHHH",secs)//date2string(dfmt,ts1970),&
    date2string(dfmt,intdate)              
end subroutine self_test

subroutine assign_startandenddate()
  ! ensure that a valid day of the month,
  ! e.g. Feb 31=>Feb 28/29 depending the year
  startdate(3)=min(startdate(3),max_day(startdate(2),startdate(1)))
  enddate  (3)=min(enddate  (3),max_day(enddate  (2),enddate  (1)))

  if(startdate(1)/=enddate  (1))then
     write(*,*)'WARNING: start and end dates in different years. Model not testet for that!'
  endif

end subroutine assign_startandenddate

function date_is_reached(date_limit) result(is_reached)
  integer, intent(in), dimension(:)      ::date_limit
  logical is_reached
  TYPE(timestamp)   :: ts1,ts2

  ts1=make_timestamp(current_date)

  if(size(date_limit)==4)then
     ts2=make_timestamp(date(date_limit(1),date_limit(2),date_limit(3),date_limit(4),0))
  elseif(size(date_limit)==5)then
     ts2=make_timestamp(date(date_limit(1),date_limit(2),date_limit(3),date_limit(4),date_limit(5)))
  else
     call CheckStop('date format not supported ')
  endif

  !is_reached =  (nint(tdif_secs(ts1,ts2))<=0)  !gives floating point exception sometime??
  is_reached = ((ts2%jdate-ts1%jdate)*3600.0*24.0+(ts2%secs-ts1%secs)<=1.0)

end function date_is_reached

ENDMODULE TimeDate_ExtraUtil_mod

!DSX program tester
!DSX   use TimeDate_ExtraUtil_mod, only : self_test
!DSX   call self_test()
!DSX end program tester
