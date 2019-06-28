! <OutputChem_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module OutputChem_mod

use CheckStop_mod,      only: CheckStop
use Config_module,     only: num_lev3d, MasterProc, runlabel1,&
                             FREQ_HOURLY, END_OF_EMEPDAY, METSTEP, &
                             IOU_INST, IOU_YEAR, IOU_MON, IOU_DAY, IOU_HOUR, &
                             IOU_HOUR_INST, IOU_MAX_MAX, HOURLYFILE_ending, &
                             startdate, enddate, out_startdate, USE_uEMEP, USES&
                             ,SITE_XTRA_D2D
use Debug_module,       only: DEBUG => DEBUG_OUTPUTCHEM
use Derived_mod,        only: LENOUT2D, nav_2d, num_deriv2d  &
                            ,LENOUT3D, nav_3d, num_deriv3d  &
                            ,wanted_iou, ResetDerived
use DerivedFields_mod,  only: f_2d, d_2d, f_3d, d_3d
use GridValues_mod,     only: debug_proc ,debug_li, debug_lj
use Io_mod,             only: IO_WRTCHEM, IO_TMP, datewrite
use NetCDF_mod,         only: CloseNetCDF, Out_netCDF, filename_iou, Init_new_netCDF
use OwnDataTypes_mod,   only: Deriv, print_deriv_type, TXTLEN_FILE
use Par_mod,            only: LIMAX, LJMAX, me
use SmallUtils_mod,     only: find_duplicates
use TimeDate_mod,       only: tdif_secs,date,timestamp,make_timestamp,current_date, max_day &! days in month
                             ,daynumber,add2current_date,date
use TimeDate_ExtraUtil_mod,only: date2string, date_is_reached
use uEMEP_mod,          only: out_uEMEP
use Units_mod,          only: Init_Units

implicit none

!** subroutines:
public :: Wrtchem
public :: Output_fields   ! (iotyp)
public :: Output_f2d      ! (iotyp, dim, nav, def, dat)
public :: Output_f3d      ! (iotyp, dim, nav, def, dat)
public :: set_output_defs

!*** wanted binary dates... specify days for which full binary
!    output is wanted. Replaces the hard-coding which was in wrtchem:
integer, public, parameter :: NBDATES = 3
type(date), public, save, dimension(NBDATES) :: wanted_dates_inst

contains

subroutine Wrtchem(ONLY_HOUR)
!---------------------------------------------------------------------
! DESCRIPTION:
!   Writes out data fields as NetCDF
!
!   END_OF_EMEPDAY = 6am, i.g. EMEP sampling period from 6am to 6am
!   Daily outputs for "EMEP" days (which end between 0 and 6am) are
!   dated by the date of sampling start, i.e. date of the previous day.
!
!   Thus, the first output should occur just as Jan 2nd starts (e.g.
!   at 6am on 2nd Jan); logical Jan_1st helps dealing with this and
!   it also marks end of a year run.
!   (For runs starting in other months, one partial write-out will
!   occur at 6am of the 1st day, but this should be over-written
!   as soon as a full day of data is available).
!----------------------------------------------------------------------
  integer, intent(in), optional :: ONLY_HOUR  ! output hourly fields

  integer :: n
  integer :: nyear,nmonth,nday,nhour,nmonpr
  integer :: mm_out, dd_out
  logical :: Jan_1st, End_of_Run
  logical,save :: first_call = .true.
  integer, save :: daynumber_old=-999
  character(len=TXTLEN_FILE) :: newname
  TYPE(timestamp)   :: ts1,ts2

!---------------------------------------------------------------------

  if(.not. date_is_reached(out_startdate))then
     if(me==0)write(*,*)'not yet outputting'
     return ! we do not output yet
  endif
  nyear  = current_date%year
  nmonth = current_date%month
  nday   = current_date%day
  nhour  = current_date%hour

  dd_out = nday
  mm_out = nmonth
  Jan_1st = all((/nyear,nmonth,nday/)==startdate(1:3))

!this is a bit complicated because it must account for the fact that for instance 3feb24:00 = 4feb00:00 
  ts1=make_timestamp(current_date)
  ts2=make_timestamp(date(enddate(1),enddate(2),enddate(3),enddate(4),0))
  End_of_Run =  (nint(tdif_secs(ts1,ts2))<=0)
  End_of_Run =  date_is_reached(enddate)

  if((current_date%seconds /= 0 ).and. .not. End_of_Run)return
  if(MasterProc .and. DEBUG) write(6,"(a12,i5,5i4)") "DAILY DD_OUT ",   &
       nmonth, mm_out, nday, dd_out, nhour

  !. END_OF_EMEPDAY = 6am - end of EMEP daily sampling period
  !. Daily outputs are dated with the start of sampling period
  if ( END_OF_EMEPDAY  <= 7 ) then
    dd_out = nday - 1     ! only used for daily outputs

    if( MasterProc .and. DEBUG) write(6,"(a12,i5,5i4)")&
      "DAILY SET ",  nmonth, mm_out, nday, dd_out, nhour

    if(dd_out == 0) then
      mm_out = nmonth - 1

      if(nmonth == 1) mm_out = 12

      dd_out = max_day(mm_out, nyear)  !  Last day of month

      if( MasterProc .and. DEBUG) write(6,"(a12,i5,4i4)") "DAILY FIX ", &
                       nmonth, mm_out, nday, dd_out
    end if
  end if      ! for END_OF_EMEPDAY <= 7

  !we substract one second, so that the day number of midnight is for the day 
  !which has passed, not the next day
  newname = trim(runlabel1)//'_hour'//date2string(HOURLYFILE_ending,current_date,-1.0)
  if(wanted_iou(IOU_HOUR).and.fileName_iou(IOU_HOUR)/=newname)then
     call Init_new_netCDF(trim(newname),IOU_HOUR)
  endif
  newname = trim(runlabel1)//'_hourInst'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
  if(wanted_iou(IOU_HOUR_INST).and.fileName_iou(IOU_HOUR_INST)/=trim(newname))then
     call Init_new_netCDF(trim(newname),IOU_HOUR_INST)
  endif


  !== Instantaneous results output ====
  !   Possible actual array output for specified days and hours
  !   is defined in wanted_dates_bi array in My_Outputs
  do n = 1, NBDATES
    if ( wanted_dates_inst(n)%month == nmonth .and. &
         wanted_dates_inst(n)%day   == nday   .and. &
         wanted_dates_inst(n)%hour  == nhour ) then
      call Output_fields(IOU_INST)
    end if
  end do

  !== Hourly output ====
  if(modulo(current_date%hour,FREQ_HOURLY)==0) then
    call Output_fields(IOU_HOUR_INST)
    if(present(ONLY_HOUR))then
      if(ONLY_HOUR==IOU_HOUR_INST)return
    end if
    call Output_fields(IOU_HOUR)
    call ResetDerived(IOU_HOUR) 
    if(present(ONLY_HOUR))then
      if(ONLY_HOUR==IOU_HOUR)return
    end if
  end if

  !== Daily output ====
  if (nhour ==  END_OF_EMEPDAY ) then
    if (.not.first_call .and. .not.Jan_1st ) &   ! Doesn't write out 1 Jan. at start
      call Output_fields(IOU_DAY)
    call ResetDerived(IOU_DAY)            ! For daily averaging, reset also 1 Jan.
  end if

  !== Output at the end of the run
  if ( End_of_Run ) then
     if(nhour/=END_OF_EMEPDAY) call Output_fields(IOU_DAY)! Daily outputs
     call Output_fields(IOU_YEAR)  ! Yearly outputs
  endif

  !/ NEW MONTH
  if ( End_of_Run .or. (nday == 1 .and. nhour == 0)) then
     nmonpr = nmonth
     if(nday == 1 .and. nhour == 0)nmonpr = nmonth-1
     if (nmonpr == 0) nmonpr=12
     
     !== Monthly output ====
     call Output_fields(IOU_MON)
     
     call ResetDerived(IOU_MON)
  end if              ! End of NEW MONTH

  first_call=.false.
end subroutine Wrtchem

subroutine Output_fields(iotyp)
  integer, intent(in) :: iotyp
  logical, dimension(IOU_MAX_MAX),save       :: myfirstcall = .true.
  logical :: Init_Only
  integer :: i
  if(myfirstcall(iotyp))then
     !only predefine the fields. For increased performance 
     Init_Only = .true.
     if(num_deriv2d > 0) call Output_f2d(iotyp,num_deriv2d,nav_2d,f_2d,d_2d,Init_Only)
     if(num_deriv3d > 0) call Output_f3d(iotyp,num_deriv3d,nav_3d,f_3d,d_3d,Init_Only)
     myfirstcall(iotyp) = .false.
     IF(DEBUG.and.MasterProc)write(*,*)'2d and 3D OUTPUT INITIALIZED',iotyp
  end if
  Init_Only = .false.
  IF(DEBUG.and.MasterProc)write(*,*)'2d and 3D OUTPUT WRITING',iotyp
  !*** 2D fields, e.g. surface SO2, SO4, NO2, NO3 etc.; AOT, fluxes
  !--------------------

  if(num_deriv2d > 0) call Output_f2d(iotyp,num_deriv2d,nav_2d,f_2d,d_2d,Init_Only)

  !*** 3D concentration fields, e.g. O3
  !--------------------
  if(num_deriv3d > 0) call Output_f3d(iotyp,num_deriv3d,nav_3d,f_3d,d_3d,Init_Only)

  call CloseNetCDF

  !uemep use own outputting for now, since it has several extra dimensions
  if(USE_uEMEP)then
    call out_uEMEP(iotyp)
  endif

  ! Write text file to mark output is finished
  if(.not.all([USES%PROGRESS_FILES,MasterProc,wanted_iou(iotyp)]))return
  i=index(filename_iou(iotyp),'.nc')-1
  if(i<1)i=len_trim(filename_iou(iotyp))
  open(IO_TMP,file=filename_iou(iotyp)(1:i)//'.msg',position='append')
  write(IO_TMP,*)date2string('FFFF: YYYY-MM-DD hh',current_date)
  close(IO_TMP)
end subroutine Output_fields

subroutine Output_f2d (iotyp, dim, nav, def, dat, Init_Only)
!---------------------------------------------------------------------
! Sends fields to NetCDF output routines
!---------------------------------------------------------------------
  integer,                         intent(in) :: iotyp
  integer,                         intent(in) :: dim ! No. fields
  integer, dimension(dim,LENOUT2D),intent(in) :: nav ! No. items averaged
  type(Deriv), dimension(dim),     intent(in) :: def ! Definition of fields
  real, dimension(dim,LIMAX,LJMAX,LENOUT2D), intent(in) :: dat
  logical,                         intent(in) :: Init_Only! only define fields

  integer :: my_iotyp,icmp ! output type,component index
  real    :: scale      ! Scaling factor
!---------------------------------------------------------------------

  do icmp = 1, dim
    if ( wanted_iou(iotyp,def(icmp)%iotype) ) then
      my_iotyp=iotyp
      if(iotyp==IOU_HOUR_INST) my_iotyp=IOU_INST
      scale  = def(icmp)%scale
      if(my_iotyp/=IOU_INST) scale = scale / max(1,nav(icmp,my_iotyp))

      !if ( MasterProc .and. DEBUG ) then
      if ( DEBUG .and. debug_proc ) then
          write(*,*) "DEBUG Output_f2d ", icmp, iotyp, trim(def(icmp)%name)
          write(*,"(a,i6,2es10.3)") "Output_f2d n,Scales:"// &
            trim(def(icmp)%name), nav(icmp,my_iotyp), def(icmp)%scale, scale
          write(*,"(a,2es10.3)") "Output_f2d  max/min", &
            maxval(dat(icmp,:,:,my_iotyp)), minval(dat(icmp,:,:,my_iotyp))
          if( def(icmp)%name == "Emis_mgm2_co" ) then
            call print_deriv_type(def(icmp))
            call datewrite("SnapEmis-Output_f2d Emis", iotyp, (/ dat(icmp,debug_li,debug_lj,my_iotyp) /) )
          end if
        end if
      call Out_netCDF(iotyp,def(icmp),2,1,dat(icmp,:,:,my_iotyp),scale,&
                      create_var_only=Init_Only)
    end if     ! wanted
  end do       ! component loop

end subroutine Output_f2d

subroutine Output_f3d (iotyp, dim, nav, def, dat, Init_Only)
!---------------------------------------------------------------------
! Sends fields to NetCDF output routines
!---------------------------------------------------------------------
  implicit none
  integer,                         intent(in) :: iotyp
  integer,                         intent(in) :: dim ! No. fields
  integer, dimension(dim,LENOUT3D),intent(in) :: nav ! No. items averaged
  type(Deriv), dimension(dim),     intent(in) :: def ! definition of fields
  real, dimension(dim,LIMAX,LJMAX,num_lev3d,LENOUT3D), intent(in):: dat
  logical,                         intent(in) :: Init_Only! only define fields

  integer :: my_iotyp,icmp ! output type,component index
  real    :: scale      ! Scaling factor
!---------------------------------------------------------------------

  do icmp = 1, dim
    if ( wanted_iou(iotyp,def(icmp)%iotype) ) then
      my_iotyp=iotyp
      if(iotyp==IOU_HOUR_INST) my_iotyp=IOU_INST
      scale = def(icmp)%scale
      if(my_iotyp/=IOU_INST) scale = scale / max(1,nav(icmp,my_iotyp))

      call Out_netCDF(iotyp,def(icmp),3,num_lev3d,dat(icmp,:,:,:,my_iotyp),scale,&
                      create_var_only=Init_Only)
    end if     ! wanted
  end do       ! component loop

end subroutine Output_f3d



subroutine set_output_defs
  implicit none
  character(len=144) :: errmsg   ! Local error message
  character(len=*), parameter :: dtxt = 'OutputChem:set:' ! debug txt

  call Init_Units()
  
  ! Safety checks for common mistakes:
  errmsg = find_duplicates(SITE_XTRA_D2D)
  call CheckStop ( errmsg /= 'no', dtxt//' Duplicate SITE_XTRA_D2D'//errmsg)

  !*** Wanted dates for instantaneous values output:
  !    specify months,days,hours for which full output is wanted.
  wanted_dates_inst(1) = date(-1,1,1,0,0)
  wanted_dates_inst(2) = date(-1,1,1,3,0)
  wanted_dates_inst(3) = date(-1,1,1,6,0)

endsubroutine set_output_defs


endmodule OutputChem_mod
