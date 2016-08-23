! <OutputChem_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module OutputChem_ml

use CheckStop_ml,      only: CheckStop
use Derived_ml,        only: LENOUT2D, nav_2d, num_deriv2d  &
                            ,LENOUT3D, nav_3d, num_deriv3d  &
                            ,iou_min, iou_max, ResetDerived
use DerivedFields_ml,  only: f_2d, d_2d, f_3d, d_3d
use GridValues_ml,     only: debug_proc ,debug_li, debug_lj
use My_Outputs_ml,     only: NBDATES, wanted_dates_inst,            &
                             Ascii3D_WANTED
use Io_ml,             only: IO_WRTCHEM, datewrite
use ModelConstants_ml, only: nprint, END_OF_EMEPDAY, KMAX_MID, MasterProc&
                            ,DEBUG => DEBUG_OUTPUTCHEM &
                            ,IOU_INST, IOU_YEAR, IOU_MON, IOU_DAY, IOU_MAX_MAX
use NetCDF_ml,         only: CloseNetCDF, Out_netCDF
use OwnDataTypes_ml,   only: Deriv, print_deriv_type
use Par_ml,            only: MAXLIMAX,MAXLJMAX,GIMAX,GJMAX,     &
                              IRUNBEG,JRUNBEG
use TimeDate_ml,       only: current_date, max_day  ! days in month
use TimeDate_ExtraUtil_ml,only: date2string


implicit none

!/* subroutines:
public :: Wrtchem
public :: Output_fields   ! (iotyp)
public :: wanted_iou      ! (iotyp, def%iotyp)
public :: Output_f2d      ! (iotyp, dim, nav, def, dat)
public :: Output_f3d      ! (iotyp, dim, nav, def, dat)

contains

subroutine Wrtchem(numt)
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
  integer, intent(in) ::  numt

  real, dimension(MAXLIMAX, MAXLJMAX) :: local_2d  !local 2D array
  real, dimension(GIMAX, GJMAX)       :: glob_2d   !array for whole domain
  integer :: i,j,n,k,msnr1
  integer :: nyear,nmonth,nday,nhour,nmonpr
  integer :: mm_out, dd_out
  logical :: Jan_1st, End_of_Run
  character(len=30) :: outfilename
!---------------------------------------------------------------------
  nyear  = current_date%year
  nmonth = current_date%month
  nday   = current_date%day
  nhour  = current_date%hour

  dd_out = nday
  mm_out = nmonth
  Jan_1st    = ( nmonth == 1 .and. nday == 1 )
  End_of_Run = ( mod(numt,nprint) == 0       )

  if(MasterProc .and. DEBUG) write(6,"(a12,i5,5i4)") "DAILY DD_OUT ",   &
      numt, nmonth, mm_out, nday, dd_out, nhour

  !. END_OF_EMEPDAY = 6am - end of EMEP daily sampling period
  !. Daily outputs are dated with the start of sampling period
  if ( END_OF_EMEPDAY  <= 7 ) then
    dd_out = nday - 1     ! only used for daily outputs

    if( MasterProc .and. DEBUG) write(6,"(a12,i5,5i4)")&
      "DAILY SET ", numt, nmonth, mm_out, nday, dd_out, nhour

    if(dd_out == 0) then
      mm_out = nmonth - 1

      if(nmonth == 1) mm_out = 12

      dd_out = max_day(mm_out, nyear)  !  Last day of month

      if( MasterProc .and. DEBUG) write(6,"(a12,i5,4i4)") "DAILY FIX ", &
                      numt, nmonth, mm_out, nday, dd_out
    endif
  endif      ! for END_OF_EMEPDAY <= 7

  !== Instantaneous results output ====
  !   Possible actual array output for specified days and hours
  !   is defined in wanted_dates_bi array in My_Outputs
  do n = 1, NBDATES
    if ( wanted_dates_inst(n)%month == nmonth .and. &
         wanted_dates_inst(n)%day   == nday   .and. &
         wanted_dates_inst(n)%hour  == nhour ) then
      call Output_fields(IOU_INST)
    endif
  enddo


  !== Daily output ====
  if (nhour ==  END_OF_EMEPDAY ) then
    if (numt > 1 .and. .not.Jan_1st ) &   ! Doesn't write out 1 Jan.
      call Output_fields(IOU_DAY)
    call ResetDerived(IOU_DAY)            ! For daily averaging, reset also 1 Jan.
  endif

  !== Output at the end of the run
  if ( End_of_Run ) then
    if(nhour/=END_OF_EMEPDAY) call Output_fields(IOU_DAY)! Daily outputs
    call Output_fields(IOU_YEAR)  ! Yearly outputs
  endif


  !/ NEW MONTH
  if (nday == 1 .and. nhour == 0) then
    nmonpr = nmonth-1
    if (nmonpr == 0) nmonpr=12

    !== Monthly output ====
    call Output_fields(IOU_MON)

    !== ASCII output of 3D fields (if wanted)
    if(Ascii3D_WANTED.and.num_deriv3d > 0) then
      msnr1 = 2000

      do n = 1, num_deriv3d
        if( MasterProc ) then
          outfilename=date2string(trim(f_3d(n)%name)//".out.MM",month=nmonpr)
          open (IO_WRTCHEM,file=outfilename)
          write(IO_WRTCHEM,fmt="(4i4)") &
            IRUNBEG, GIMAX+IRUNBEG-1, JRUNBEG, GJMAX+JRUNBEG-1 ! domain
        endif

        if (nav_3d(n,IOU_MON) == 0 ) then
          write(IO_WRTCHEM,*) "ERROR in 3D ASCII output: nav=0"
        else
          do k = 1, KMAX_MID
            local_2d(:,:) = d_3d(n,:,:,k,IOU_MON)/nav_3d(n,IOU_MON)
            call local2global(local_2d,glob_2d,msnr1)

            if( MasterProc ) &
              write(IO_WRTCHEM,"(es10.3)") ((glob_2d(i,j),i=1,GIMAX),j=1,GJMAX)

          enddo ! k
        endif   ! nav == 0

        if( MasterProc ) close(IO_WRTCHEM)

      enddo     ! 3D-variables loop num_deriv3d
    endif       ! Ascii3D_WANTED

    call ResetDerived(IOU_MON)
  endif              ! End of NEW MONTH

end subroutine Wrtchem

subroutine Output_fields(iotyp)
  integer, intent(in) :: iotyp
  logical, dimension(IOU_MAX_MAX),save       :: myfirstcall = .true.
  logical             :: Init_Only
  if(myfirstcall(iotyp))then
     !only predefine the fields. For increased performance 
     Init_Only = .true.
     if(num_deriv2d > 0) call Output_f2d(iotyp,num_deriv2d,nav_2d,f_2d,d_2d,Init_Only)
     if(num_deriv3d > 0) call Output_f3d(iotyp,num_deriv3d,nav_3d,f_3d,d_3d,Init_Only)
     myfirstcall(iotyp) = .false.
     IF(DEBUG.and.MasterProc)write(*,*)'2d and 3D OUTPUT INITIALIZED',iotyp
  endif
  Init_Only = .false.
  IF(DEBUG.and.MasterProc)write(*,*)'2d and 3D OUTPUT WRITING',iotyp
  !*** 2D fields, e.g. surface SO2, SO4, NO2, NO3 etc.; AOT, fluxes
  !--------------------
  if(num_deriv2d > 0) call Output_f2d(iotyp,num_deriv2d,nav_2d,f_2d,d_2d,Init_Only)

  !*** 3D concentration fields, e.g. O3
  !--------------------
  if(num_deriv3d > 0) call Output_f3d(iotyp,num_deriv3d,nav_3d,f_3d,d_3d,Init_Only)

  call CloseNetCDF
end subroutine Output_fields

function wanted_iou(iou,iotype) result(wanted)
  integer, intent(in)           :: iou
  integer, intent(in), optional :: iotype
  logical                       :: wanted
  wanted=(iou>=iou_min).and.(iou<=iou_max)
  if(present(iotype))wanted=wanted.and.(iou<=iotype)
end function wanted_iou

subroutine Output_f2d (iotyp, dim, nav, def, dat, Init_Only)
!---------------------------------------------------------------------
! Sends fields to NetCDF output routines
!---------------------------------------------------------------------
  integer,                         intent(in) :: iotyp
  integer,                         intent(in) :: dim ! No. fields
  integer, dimension(dim,LENOUT2D),intent(in) :: nav ! No. items averaged
  type(Deriv), dimension(dim),     intent(in) :: def ! Definition of fields
  real, dimension(dim,MAXLIMAX,MAXLJMAX,LENOUT2D), intent(in) :: dat
  logical,                         intent(in) :: Init_Only! only define fields

  integer :: icmp       ! component index
  real    :: scale      ! Scaling factor
!---------------------------------------------------------------------

  do icmp = 1, dim
    if ( wanted_iou(iotyp,def(icmp)%iotype) ) then
      scale  = def(icmp)%scale
      if (iotyp /= IOU_INST ) scale = scale / max(1,nav(icmp,iotyp))

      !if ( MasterProc .and. DEBUG ) then
      if ( DEBUG .and. debug_proc ) then
          write(*,*) "DEBUG Output_f2d ", icmp, iotyp, trim(def(icmp)%name)
          write(*,"(a,i6,2es10.3)") "Output_f2d n,Scales:"// &
            trim(def(icmp)%name), nav(icmp,iotyp), def(icmp)%scale, scale
          write(*,"(a,2es10.3)") "Output_f2d  max/min", &
            maxval(dat(icmp,:,:,iotyp)), minval(dat(icmp,:,:,iotyp))
          if( def(icmp)%name == "Emis_mgm2_co" ) then
            call print_deriv_type(def(icmp))
            call datewrite("SnapEmis-Output_f2d Emis", iotyp, (/ dat(icmp,debug_li,debug_lj,iotyp) /) )
          endif
        endif

      call Out_netCDF(iotyp,def(icmp),2,1,dat(icmp,:,:,iotyp),scale,create_var_only=Init_Only)
    endif     ! wanted
  enddo       ! component loop

end subroutine Output_f2d

subroutine  Output_f3d (iotyp, dim, nav, def, dat, Init_Only)
!---------------------------------------------------------------------
! Sends fields to NetCDF output routines
!---------------------------------------------------------------------
  implicit none
  integer,                         intent(in) :: iotyp
  integer,                         intent(in) :: dim ! No. fields
  integer, dimension(dim,LENOUT3D),intent(in) :: nav ! No. items averaged
  type(Deriv), dimension(dim),     intent(in) :: def ! definition of fields
  real, dimension(dim,MAXLIMAX,MAXLJMAX,KMAX_MID,LENOUT3D), intent(in):: dat
  logical,                         intent(in) :: Init_Only! only define fields

  integer :: icmp       ! component index
  real    :: scale      ! Scaling factor
!---------------------------------------------------------------------

  do icmp = 1, dim
    !FEB2011. QUERY on INST ??
    if ( wanted_iou(iotyp,def(icmp)%iotype) ) then
      scale = def(icmp)%scale
      if (iotyp /= IOU_INST) scale = scale /max(1,nav(icmp,iotyp))

      call Out_netCDF(iotyp,def(icmp),3,KMAX_MID,dat(icmp,:,:,:,iotyp),scale,create_var_only=Init_Only)
    endif     ! wanted
  enddo       ! component loop

end subroutine Output_f3d

end module OutputChem_ml
