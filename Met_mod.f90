! <Met_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Met_mod
!=============================================================================
!  Subroutines:        Frequency  Called from:
!    MetModel_LandUse             emepctm
!    MeteoRead         3h         emepctm    - puts data into nr
!    metfieldint       dt_advec   emepctm
!    BLPhysics(numt)   3h         MeteoRead, after met_derived
!
! Alt:
!  emepctm  do every dt_advec , ....
!    call MeteoRead - puts data into nr
!          Unit changes, special definitions etc...
!          e.g. converts wind to u_xmj, v_xmi
!          call met_derived(nr=2)      uses q(...,nr=2) (future)
!          call BLPhysics(nr=2)
!          call met_derived(nr=1)      uses q(...,nr=1) (now)
!    call phyche() ...
!          chemistry stuff
!    call metfieldint
!        call met_derived(nr=1)
!
!  metfieldint
!    does the forward linear stepping of the meteorological
!    fields read or derived every 3 hours.
!
!=============================================================================

use BLPhysics_mod,         only: &
   KZ_MINIMUM, KZ_SBL_LIMIT, PIELKE   &
!NWPHMIX  ,HmixMethod, UnstableKzMethod, StableKzMethod, KzMethod  &
  , UnstableKzMethod, StableKzMethod, KzMethod  &
  ,USE_MIN_KZ              & ! From old code, is it needed?
  ,MIN_USTAR_LAND          & ! sets u* > 0.1 m/s over land
  ,OB_invL_LIMIT           & !
  ,Test_BLM                & ! Tests all Kz, Hmix routines
  ,JericevicRiB_Hmix0      & ! Used, now allows shallow SBL
  ,SeibertRiB_Hmix_3d      &
  ,BrostWyngaardKz         &
  ,JericevicKz             &
  ,TI_Hmix                 &
  ,PielkeBlackadarKz       &
  ,O_BrienKz               &
  ,NWP_Kz                  & ! Kz from meteo
  ,Kz_m2s_toSigmaKz        &
  ,Kz_m2s_toEtaKz        &
  ,SigmaKz_2_m2s

use CheckStop_mod,         only: CheckStop
use Config_module,    only: PASCAL, PT, Pref, METSTEP  &
     ,PBL & ! Has ZiMIN, ZiMAX, HmixMethod
     ,KMAX_BND, KMAX_MID, MasterProc, nmax  &
     ,GRID & !HIRHAM,EMEP,EECCA etc.
     ,TEGEN_DATA, USES &
     ,step_main,USE_EtaCOORDINATES,USE_FASTJ &
     ,CONVECTION_FACTOR &
     ,LANDIFY_MET,MANUAL_GRID  &
     ,CW_THRESHOLD,RH_THRESHOLD, CW2CC, JUMPOVER29FEB, meteo, startdate&
     ,SoilTypesFile, Soil_TegenFile, TopoFile, SurfacePressureFile
use Debug_module,       only: DEBUG,DEBUG_BLM, DEBUG_Kz, &
                              DEBUG_SOILWATER, DEBUG_LANDIFY 
use FastJ_mod,          only: setup_phot_fastj,rcphot_3D
use Functions_mod,      only: Exner_tab, Exner_nd
use Functions_mod,      only: T_2_Tpot, StandardAtmos_kPa_2_km 
use GridValues_mod,     only: glat, xm_i, xm_j, xm2         &
       ,Poles, sigma_bnd, sigma_mid, xp, yp, fi, GRIDWIDTH_M  &
       ,debug_proc, debug_li, debug_lj, A_mid, B_mid          &
       ,Eta_bnd,Eta_mid,dA,dB,A_mid,B_mid,A_bnd,B_bnd         &
       ,KMAX_MET,External_Levels_Def,k1_met,k2_met,x_k1_met,rot_angle&
       ,Meteo_Get_KMAXMET,remake_vertical_levels_interpolation_coeff

use Io_mod ,            only: ios, datewrite, PrintLog, IO_LOG
use Landuse_mod,        only: water_fraction, water_frac_set, &
                                likely_coastal, mainly_sea
use MetFields_mod
use MicroMet_mod,       only: PsiH  ! Only if USE_MIN_KZ
use MPI_Groups_mod,     only: MPI_DOUBLE_PRECISION, MPI_BYTE, MPI_LOGICAL,&
                             MPI_COMM_IO, MPI_COMM_CALC, IERROR, ME_IO, ME_CALC,&
                             request_e,request_n,request_s,request_w,LargeSub_Ix,&
                             largeLIMAX,largeLJMAX, MPISTATUS, MPI_MIN
use OwnDataTypes_mod,      only: TXTLEN_FILE
use Par_mod,               only: MAXLIMAX,MAXLJMAX,GIMAX,GJMAX, me  &
     ,limax,ljmax, neighbor,WEST,EAST,SOUTH,NORTH,NOPROC  &
     ,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2  &
     ,IRUNBEG,JRUNBEG, gi0,gj0, MSG_READ4, tlargegi0,tlargegj0
use PhysicalConstants_mod, only: KARMAN, KAPPA, RGAS_KG, CP, GRAV
use TimeDate_mod,          only: current_date, date,nmdays, &
     add_secs,timestamp, make_timestamp, make_current_date
use NetCDF_mod,        only: ReadField_CDF,vertical_interpolate,GetCDF_modelgrid,ReadTimeCDF
use netcdf
use TimeDate_ExtraUtil_mod,only: nctime2date,date2string,date2nctime

implicit none
private

! logical, private, save      :: debug_procloc = .false.
integer, private, save      :: debug_iloc, debug_jloc  ! local coords
integer, save   :: nrec          ! nrec=record to read in meteofile
! Nhh is the number of record in each meteo file (assume constant)
! (Nhh=8): nrec=1->00:00 nrec=2->03:00 ... nrec=8->21:00
! if nhour_first=3 then nrec=1->03:00 nrec=2->06:00...nrec=8->24:00
integer, save   :: nrec_mult = 1     ! How many records to jump over at each metstep

logical, save, private      :: xwf_done = .false. ! extended water-fraction array

character(len=*),parameter  :: field_not_found='field_not_found'
integer(kind=2),allocatable :: var_global(:,:,:)
integer(kind=2),allocatable :: var_local(:,:,:)

! Aid for debugging check routine
character (len = 100), private, save :: call_msg=" Not set"

public :: MeteoRead, MeteoRead_io
public :: MetModel_LandUse
public :: metfieldint
public :: BLPhysics
public :: GetCDF_short
public :: extendarea  ! returns array which includes neighbours
public :: Getmeteofield
public :: landify     ! replaces met variables from mixed sea/land with land

contains

subroutine MeteoRead_io()

  character(len=TXTLEN_FILE), save :: meteoname   ! name of the meteofile
  character(len=100) ::  namefield  ! name of the requested field
  character(len=*), parameter :: dtxt='MetRd_io:'
  integer :: ix, KMAX, istart,jstart,ijk,i,j,k,k1,k2
  integer :: nr
  integer ::   ndim,nyear,nmonth,nday,nhour
  real ::meteo_3D(largeLIMAX,largeLJMAX,KMAX_MET)
  logical,save :: first_call = .true.
  type(date)      ::  next_inptime             ! hfTD,addhours_to_input
  type(timestamp) ::  ts_now                   ! time in timestamp format
  real :: nsec                                 ! step in seconds
  logical :: fexist

  if(current_date%seconds /= 0 .or. (mod(current_date%hour,METSTEP)/=0) )return

  nr=2 !set to one only when the first time meteo is read
  call_msg = dtxt//"Meteoread"
  if(me_IO>=0)then
    if(first_call)then !first time meteo is read
      nr = 1
      nrec = 0
      next_inptime = current_date

      !On first call, check that date from meteo file correspond to dates requested.
      !Also defines nhour_first and Nhh (and METSTEP in case of WRF metdata).
      call Check_Meteo_Date_Type !note that all procs read this
    else
      nsec=METSTEP*3600.0 !from hr to sec
      ts_now = make_timestamp(current_date)
      call add_secs(ts_now,nsec)
      if(JUMPOVER29FEB.and.current_date%month==2.and.current_date%day==29)then
        if(MasterProc)write(*,*)'Jumping over one day for meteo_date!'
        call add_secs(ts_now,24*3600.)
      end if
      next_inptime=make_current_date(ts_now)
    end if
    nyear=next_inptime%year
    nmonth=next_inptime%month
    nday=next_inptime%day
    nhour=next_inptime%hour
    nrec=nrec+nrec_mult

    if(nrec>Nhh.or.nrec==1) then              ! start reading a new meteo input file
      meteoname = date2string(meteo,next_inptime,mode='YMDH')
      nrec = 1
      if(nday==1.and.nmonth==1)then
        !hour 00:00 from 1st January may be missing;checking first:
        inquire(file=meteoname,exist=fexist)
        if(.not.fexist)then
          if(MasterProc)write(*,*)trim(meteoname),&
               ' does not exist; using data from previous day'
          meteoname=date2string(meteo,next_inptime,-24*3600.0,mode='YMDH')
          nrec=Nhh
        end if
      end if
      if(ME_IO==0)write(*,*)'io procs reading ',trim(meteoname)
    end if
    if(nrec_mult/=1 .and. MasterProc)write(*,*)'reading record ',nrec

    do ix=1,Nmetfields
      if(met(ix)%read_meteo)then
        namefield=met(ix)%name
        ndim=met(ix)%dim

        if(ndim==3)KMAX=KMAX_MET
        if(ndim==2)KMAX=1

        istart=gi0
        jstart=gj0
        call GetCDF_modelgrid(namefield,meteoname,meteo_3D,1,KMAX,nrec,1,&
            imax_in=largeLIMAX,jmax_in=largeLJMAX,needed=met(ix)%needed,&
            found=met(ix)%found)
        if(met(ix)%found)then
          if(KMAX==1)then
            ijk=0
            k=1
            do j=1,largeLJMAX
              do i=1,largeLIMAX
                ijk=ijk+1
                met(ix)%field_shared(i,j,k)=meteo_3D(i,j,k)
              end do
            end do
          else
            if(External_Levels_Def)then
              !interpolate vertically if the levels are not identical
              ijk=0
              do k=1,KMAX_MID
                k1=k1_met(k)
                k2=k2_met(k)
                do j=1,largeLJMAX
                  do i=1,largeLIMAX
                    ijk=ijk+1
                    met(ix)%field_shared(i,j,k)=x_k1_met(k)*meteo_3D(i,j,k1)&
                                         +(1.0-x_k1_met(k))*meteo_3D(i,j,k2)
                  end do
                end do
              end do
            else
              !use same vertical coordinates as meteo
              ijk=0
              do k=1,KMAX_MID
                do j=1,largeLJMAX
                  do i=1,largeLIMAX
                    ijk=ijk+1
                    met(ix)%field_shared(i,j,k)=meteo_3D(i,j,k)
                  end do
                end do
              end do
            end if
          end if
          met(ix)%ready=.true.
          met(ix)%copied=.false.
        else
          met(ix)%ready=.false.
          met(ix)%copied=.false.
        end if

        if(me_io==0)then
          if(met(ix)%found)write(*,*)'found ',trim(namefield),' in ',trim(meteoname)
          if(met(ix)%found.and.ndim==2)write(*,*)'typical value 2D = ',&
                trim(namefield),me_io,met(ix)%field_shared(5,5,1)
          if(met(ix)%found.and.ndim==3)write(*,*)'typical value 3D = ',&
                trim(namefield),me_io,met(ix)%field_shared(5,5,KMAX_MID)
          if(.not.met(ix)%found)write(*,*)'did not find ',trim(namefield),&
                ' in ',trim(meteoname)
        end if
      end if
    end do
    first_call = .false.
  end if

end subroutine MeteoRead_io

subroutine MeteoRead()
  !    the subroutine reads meteorological fields and parameters (every
  !    METSTEP-hours) from NetCDF fields-files, divide the fields into
  !       domains    and sends subfields to the processors
  implicit none

  character(len=TXTLEN_FILE), save :: meteoname   ! name of the meteofile
  character(len=100) ::  namefield & ! name of the requested field
       ,unit='   ',validity='    '    ! field is either instaneous or averaged
  integer ::   ndim,nyear,nmonth,nday,nhour
  integer ::   nr   ! Fields are interpolate in
  character(len=*), parameter :: dtxt='MetRd:'

  type(date)      ::  next_inptime             ! hfTD,addhours_to_input
  type(timestamp) ::  ts_now                   ! time in timestamp format

  real :: nsec                                 ! step in seconds

  real :: buff(LIMAX,LJMAX)!temporary metfields
  integer :: i, j, ix, k, kk, nrix, isw, KMAX
  logical :: fexist
  logical,save :: first_call = .true.

  ! Soil water has many names. Some we can deal with:
  ! (and all need to end up as SMI)
  character(len=*),  parameter :: &
    possible_soilwater_uppr(4) = (/&
        "SMI1                " &
       ,"SMI                 " &
       ,"soil_water_content  " &
       ,"soil_wetness_surface" /), &
    possible_soilwater_deep(3) = (/&
        "SMI3                   " &
       ,"SMI                    " &
       ,"deep_soil_water_content" /)
   ! Hmix ditto:
  character(len=*),  parameter :: &
    possible_HmixNames(2) =  [ "pblh", &
      "blh " ] ! GLOBAL05

  logical :: write_now

  real :: relh1,relh2,temperature,swp,wp, x_out
  real, dimension(KMAX_MID)       ::  exf2
  real, dimension(KMAX_BND)       ::  exf1

  real, dimension(LJMAX,KMAX_MID) :: usnd   ! send in x
  real, dimension(LIMAX,KMAX_MID) :: vsnd   ! and in y direction
  real, dimension(LJMAX,KMAX_MID) :: urcv   ! rcv in x
  real, dimension(LIMAX,KMAX_MID) :: vrcv   ! and in y direction

  real :: p1, p2, x, y
  real :: divk(KMAX_MID),sumdiv,dB_sum

  real :: Ps_extended(0:LIMAX+1,0:LJMAX+1),Pmid,Pu1,Pu2,Pv1,Pv2

  real :: minprecip, tmpmax ! debug

  real, dimension(LJMAX,KMAX_MID) :: buf_uw,buf_ue
  real, dimension(LIMAX,KMAX_MID) :: buf_vn,buf_vs
  integer :: INFO,i_large,j_large
  logical, save:: ps_in_hPa = .true.
  logical, save:: precip_accumulated = .false.
  
  if(.not. first_call)then
     if(current_date%seconds /= 0 .or. (mod(current_date%hour,METSTEP)/=0) )return
  endif

  nr=2 !set to one only when the first time meteo is read
  call_msg = dtxt//"Meteoread"

  nrec=nrec+nrec_mult

  if(first_call)then !first time meteo is read
     nr = 1
     nrec = 1
     next_inptime = current_date

     KMAX=max(KMAX_MID,KMAX_MET)!so that allocated arrays are large for both use
     if(MasterProc)then
       allocate(var_global(GIMAX,GJMAX,KMAX))
     else
       allocate(var_global(1,1,1)) !just to have the array defined
     end if
     allocate(var_local(MAXLIMAX,MAXLJMAX,KMAX))

     !On first call, check that date from meteo file correspond to dates requested.
     !Also defines nhour_first and Nhh (and METSTEP and bucket in case of WRF metdata).
     !Also check if data is in format short
     call Check_Meteo_Date_Type !note that all procs read this

     call Exner_tab()!init table

     debug_iloc = debug_li
     debug_jloc = debug_lj

  else
    nsec=METSTEP*3600.0 !from hr to sec
    ts_now = make_timestamp(current_date)
    call add_secs(ts_now,nsec)
    if(JUMPOVER29FEB.and.current_date%month==2.and.current_date%day==29)then
      if(MasterProc)write(*,*)dtxt//'Jumping over one day for meteo_date!'
      call add_secs(ts_now,24*3600.)
    end if
    next_inptime=make_current_date(ts_now)
  end if

  nyear=next_inptime%year
  nmonth=next_inptime%month
  nday=next_inptime%day
  nhour=next_inptime%hour

  if(MasterProc.and.DEBUG%MET) &
    write(6,*) dtxt//'*** nyear,nmonth,nday,nhour,nmdays2'    &
      ,next_inptime%year,next_inptime%month,next_inptime%day    &
      ,next_inptime%hour,nmdays(2)

  !Read rec=1 both for h=0 and h=3:00 in case 00:00 in 1st meteofile

  if(nrec>Nhh.or.nrec==1.or.first_call) then              ! start reading a new meteo input file
    meteoname = date2string(meteo,next_inptime,mode='YMDH')

    !check if the number of vertical levels has changed
    call Meteo_Get_KMAXMET(meteoname,kmax)

    if(kmax/=KMAX_MET)then
       if(me==0)write(*,*)'WARNING: number of vertical levels in meteo file has changed from ',KMAX_MET,' to ', KMAX
       if(me==0)write(*,*)'recalculating vertical interpolation coefficients'
       if(kmax>KMAX_MET .and. kmax>KMAX_MID)then
          !resize arrays
          if(MasterProc)then
             deallocate(var_global)
             allocate(var_global(GIMAX,GJMAX,KMAX))
          end if
          deallocate(var_local)
          allocate(var_local(MAXLIMAX,MAXLJMAX,KMAX))
       endif

       call remake_vertical_levels_interpolation_coeff(meteoname)

    endif

    if(.not.first_call)nrec = 1
    if(nday==1.and.nmonth==1)then
      !hour 00:00 from 1st January may be missing;checking first:
      inquire(file=meteoname,exist=fexist)
      if(.not.fexist)then
        if(MasterProc)write(*,*)dtxt//trim(meteoname),&
            ' does not exist; using data from previous day'
        meteoname=date2string(meteo,next_inptime,-24*3600.0,mode='YMDH')
        nrec=Nhh
      end if
    end if
    if(MasterProc)write(*,*)dtxt//'reading ',trim(meteoname)
    !could open and close file here instead of in Getmeteofield
  end if

  if(MasterProc.and.DEBUG%MET) write(*,*)dtxt//'nrec,nhour=',nrec,nhour

 !inform of what is done with each field the first time:
  write_now=MasterProc.and.(DEBUG%MET.or.first_call)

  if((nrec_mult/=1 .and. MasterProc) .or. write_now) &
      write(*,*)dtxt//'reading record ',nrec

  !==============    Read the meteo fields  ================================================

  do ix=1,Nmetfields
    if(met(ix)%read_meteo)then
      namefield=met(ix)%name
      ndim=met(ix)%dim
      nrix=min(met(ix)%msize,nr)

!     if(met(ix)%ready)then
      if(.false.)then ! DS QUERY
        select case(ndim)
        case(2)
          do j=1,ljmax
            j_large=j+gj0-tlargegj0(LargeSub_Ix)
            do i=1,limax
              i_large=i+gi0-tlargegi0(LargeSub_Ix)
              met(ix)%field(i,j,1,nrix)=met(ix)%field_shared(i_large,j_large,1)
            end do
          end do
        case(3)
          do k=1,KMAX_MID
            do j=1,ljmax
              j_large=j+gj0-tlargegj0(LargeSub_Ix)
              do i=1,limax
                i_large=i+gi0-tlargegi0(LargeSub_Ix)
                met(ix)%field(i,j,k,nrix)=met(ix)%field_shared(i_large,j_large,k)
              end do
            end do
          end do
        end select
        met(ix)%found=.true.
        met(ix_fh)%validity='averaged'!should be softified
        met(ix_fl)%validity='averaged'!should be softified
      else
        call Getmeteofield(meteoname,namefield,nrec,ndim,unit,met(ix)%validity,&
                met(ix)%field(1:LIMAX,1:LJMAX,:,nrix),needed=met(ix)%needed,&
                found=met(ix)%found)
      end if
      if(write_now)then
        if(met(ix)%found)then
          write(*,*)'found ',trim(namefield),' in ',trim(meteoname)
          select case(ndim)
          case(2)
            write(*,*)'example values = ',met(ix)%field(5,5,1,nrix),&
              maxval(met(ix)%field(:,:,1,nrix))
          case(3)
            write(*,*)'example values = ',met(ix)%field(5,5,kmax_mid,nrix),&
              maxval(met(ix)%field(:,:,kmax_mid,nrix))
          end select
        else
         write(*,*)dtxt//'did not find ',trim(namefield),' in ',trim(meteoname)
        end if
        if(me_calc<0)then
          select case(ndim)
          case(2)
            write(*,*)'met compare 2D ',me,met(ix)%field(5,5,1,nrix),&
              met(ix)%field_shared(i_large,j_large,1)
          case(3)
            write(*,*)'met compare 3D ',me,met(ix)%field(5,5,KMAX_MID,nrix),&
              met(ix)%field_shared(i_large,j_large,KMAX_MID)
          end select
        end if
      end if
    end if
  end do

  !==============  now correct and complete the metfields as needed!  ==========================
  if(MANUAL_GRID)then
    !rotate the wind fields
    !non-staggered grid here
    do k=1,KMAX_MID
      do j=1,ljmax
        do i=1,limax
          x=u_xmj(i,j,k,nr)
          y=v_xmi(i,j,k,nr)
          u_xmj(i,j,k,nr) = x*cos(rot_angle(i,j))-y*sin(rot_angle(i,j))
          v_xmi(i,j,k,nr) = x*sin(rot_angle(i,j))+y*cos(rot_angle(i,j))
        end do
      end do
    end do

    !Now must stagger the wind fields
    !must first fetch values from neighbors
    if (neighbor(WEST) .ne. NOPROC) then
      if(neighbor(WEST) .ne. me)then
        buf_uw(:,:) = u_xmj(1,:,:,nr)
        CALL MPI_ISEND(buf_uw, 8*LJMAX*KMAX_MID, MPI_BYTE, &
             neighbor(WEST), MSG_EAST2, MPI_COMM_CALC, request_w, IERROR)
      else
        ! cyclic grid: own neighbor
        ue(:,:,nr) = u_xmj(1,:,:,nr)
      end if
    end if
    if (neighbor(SOUTH) .ne. NOPROC) then
      buf_vs(:,:) = v_xmi(:,1,:,nr)
      CALL MPI_ISEND(buf_vs, 8*LIMAX*KMAX_MID, MPI_BYTE, &
          neighbor(SOUTH), MSG_NORTH2, MPI_COMM_CALC, request_s, IERROR)
    end if

    if (neighbor(EAST) .ne. NOPROC .and. neighbor(EAST) .ne. me) then
      CALL MPI_RECV(ue(1,1,nr), 8*LJMAX*KMAX_MID, MPI_BYTE, &
            neighbor(EAST), MSG_EAST2, MPI_COMM_CALC, MPISTATUS, IERROR)
    end if
    if (neighbor(NORTH) .ne. NOPROC) then
      CALL MPI_RECV(vn(1,1,nr), 8*LIMAX*KMAX_MID, MPI_BYTE, &
            neighbor(NORTH), MSG_NORTH2, MPI_COMM_CALC, MPISTATUS, IERROR)
    end if

    if (neighbor(WEST) .ne. NOPROC .and. neighbor(WEST) .ne. me) then
      CALL MPI_WAIT(request_w, MPISTATUS, IERROR)
    end if
    if (neighbor(SOUTH) .ne. NOPROC) then
      CALL MPI_WAIT(request_s, MPISTATUS, IERROR)
    end if

    do k=1,KMAX_MID
      do j=1,ljmax
        do i=1,limax-1
          u_xmj(i,j,k,nr) = 0.5*(u_xmj(i,j,k,nr)+u_xmj(i+1,j,k,nr))
        end do
      end do
    end do
    do k=1,KMAX_MID
      do j=1,ljmax
       do i=limax,limax
         u_xmj(i,j,k,nr) = 0.5*(u_xmj(i,j,k,nr)+ue(j,k,nr))
       end do
      end do
    end do

    do k=1,KMAX_MID
      do j=1,ljmax-1
        do i=1,limax
          v_xmi(i,j,k,nr) = 0.5*(v_xmi(i,j,k,nr)+v_xmi(i,j+1,k,nr))
        end do
      end do
    end do
    do k=1,KMAX_MID
      do j=ljmax,ljmax
        do i=1,limax
          v_xmi(i,j,k,nr) = 0.5*(v_xmi(i,j,k,nr)+vn(i,k,nr))
        end do
      end do
    end do
  end if ! MANUAL_GRID

  !extend the i or j index to 0
  if (neighbor(EAST) .ne. NOPROC) then
    usnd(:,:) = u_xmj(limax,:,:,nr)
    CALL MPI_ISEND( usnd, 8*LJMAX*KMAX_MID, MPI_BYTE,  &
        neighbor(EAST), MSG_WEST2, MPI_COMM_CALC, request_e, IERROR)
  end if
  if (neighbor(NORTH) .ne. NOPROC) then
    vsnd(:,:) = v_xmi(:,ljmax,:,nr)
    CALL MPI_ISEND( vsnd , 8*LIMAX*KMAX_MID, MPI_BYTE,  &
        neighbor(NORTH), MSG_SOUTH2, MPI_COMM_CALC, request_n, IERROR)
  end if

  if (neighbor(WEST) .ne. NOPROC) then
    CALL MPI_RECV( u_xmj(0,:,:,nr), 8*LJMAX*KMAX_MID, MPI_BYTE, &
        neighbor(WEST), MSG_WEST2, MPI_COMM_CALC, MPISTATUS, IERROR)
  else
    u_xmj(0,:,:,nr) = u_xmj(1,:,:,nr)
  end if
  if (neighbor(SOUTH) .ne. NOPROC) then
    CALL MPI_RECV( v_xmi(:,0,:,nr) , 8*LIMAX*KMAX_MID, MPI_BYTE,  &
        neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_CALC, MPISTATUS, IERROR)
  else
    if(Poles(2)/=1)  then
      v_xmi(:,0,:,nr) = v_xmi(:,1,:,nr)
    else
      !"close" the South pole
      v_xmi(:,0,:,nr) = 0.0
    end if
  end if
  if (neighbor(NORTH) == NOPROC.and.Poles(1)==1) then
    !"close" the North pole
    v_xmi(:,ljmax,:,nr) = 0.0
  end if
  if(neighbor(EAST) .ne. NOPROC) CALL MPI_WAIT(request_e, MPISTATUS, IERROR)
  if(neighbor(NORTH) .ne. NOPROC)CALL MPI_WAIT(request_n, MPISTATUS, IERROR)

  !divide by the scaling in the perpendicular direction to get effective
  !u_xmj and v_xmi
  !(for conformal projections like Polar Stereo, xm_i and xm_j are equal)
  do k = 1,KMAX_MID
    do j = 1,ljmax
      do i = 0,limax
        u_xmj(i,j,k,nr) = u_xmj(i,j,k,nr)/xm_j(i,j)
      end do
    end do
    do j = 0,ljmax
      do i = 1,limax
        v_xmi(i,j,k,nr) = v_xmi(i,j,k,nr)/xm_i(i,j)
      end do
    end do
  end do

  if(WRF_MET_CORRECTIONS)then
    !WRF temperatures are shifted by  300:
    th(:,:,:,nr) = th(:,:,:,nr) + 300.0
  end if

  ! correct surface pressure here, because we will need it before we get to the 2D block
  if(first_call)then
     if(maxval(ps)<2000.0)then
        ps_in_hPa = .true.
        if(write_now)write(*,*)'Asuming surface pressure in hPa'
      else
        if(write_now)write(*,*)'Asuming surface pressure in Pa'
        ps_in_hPa = .false.
     endif
  endif
  ! conversion of pressure from hPa to Pascal.
  if(ps_in_hPa)ps(1:limax,1:ljmax,nr) = ps(1:limax,1:ljmax,nr)*PASCAL


  if(first_call)then !constant->read only once
     call read_surf_elevation(ix_elev)
  endif

  if(foundcc3d)then
    if(WRF_MET_CORRECTIONS)then
      !WRF clouds in fraction, multiply by 100:
      cc3d(:,:,:) = 100*cc3d(:,:,:)
    end if
    if(trim(met(ix_cc3d)%validity)/='averaged'.and.write_now)&
      write(*,*)'WARNING: 3D cloud cover are instantaneous values'
    cc3d(:,:,:) = 0.01*max(0.0,min(100.0,cc3d(:,:,:)))!0-100 % clouds to fraction
  else !if available, will use cloudwater to determine the height of release
     if(write_now)write(*,*)'WARNING: deriving 3D cloud cover (cc3d) from cloud water '
     namefield='cloudwater'
     call Getmeteofield(meteoname,namefield,nrec,3,unit,validity,&
            cc3d(:,:,:),found=foundcloudwater)
     call CheckStop(.not.foundcloudwater,&
            "meteo field not found: 3D_cloudcover and"//trim(namefield))
     cc3d(:,:,:)=0.01*max(0.0,min(100.0,cc3d(:,:,:)*CW2CC))!from kg/kg water to % clouds to fraction
  end if
  !    maximum of cloud fractions for layers above a given layer
  cc3dmax(:,:,1) = cc3d(:,:,1)
  do k=2,KMAX_MID
    cc3dmax(:,:,k) = amax1(cc3dmax(:,:,k-1),cc3d(:,:,k-1))
  end do

  lwc = 0.6e-6*cc3d

  if(.not.foundprecip)then
    !Will construct 3D precipitations from 2D precipitations
    if(write_now)write(*,*)'WARNING: deriving 3D precipitations from 2D precipitations '

    ix = ix_surface_precip
    met(ix)%found = .false.
    if(.not. precip_accumulated)then
      call Getmeteofield(meteoname,met(ix)%name,nrec,met(ix)%dim,unit,&
           met(ix)%validity,met(ix)%field,needed=met(ix)%needed,found=met(ix)%found)
      if(met(ix)%found)then
         if(write_now )write(*,*)'2D precipitations sum of large_scale and convective precipitations'
         ix = ix_convective_precip
         call Getmeteofield(meteoname,met(ix)%name,nrec,met(ix)%dim,unit,&
              met(ix)%validity,met(ix)%field,needed=met(ix)%needed,found=met(ix)%found)
         surface_precip = surface_precip + convective_precip
      else
         !only set once
         if(write_now )write(*,*)trim(met(ix)%name),' not found. assuming accumulated'
         precip_accumulated = .true.
      endif
    endif
    ix = ix_surface_precip
    if(precip_accumulated)then
       !assume accumulated precipitations (AROME)
       ix = ix_surface_precip
       met(ix)%name = 'precipitation_amount_acc' 
       if(me==0)write(*,*)'assuming 2D precipitations accumulated and named '//trim(met(ix)%name)
       call Getmeteofield(meteoname,met(ix)%name,nrec,met(ix)%dim,unit,&
            met(ix)%validity,met(ix)%field,needed=met(ix)%needed,found=met(ix)%found)
       buff=surface_precip !save to save in old below
       if(nr == 1 .and. nhour_first>0) buff=0.0 !the first reading of metdata is not from 00:00
       !must first check that precipitation is increasing. At some dates acc maybe restarted!
       minprecip=minval(surface_precip(1:limax,1:ljmax) - surface_precip_old(1:limax,1:ljmax))
       
       CALL MPI_ALLREDUCE(minprecip, x_out, 1,MPI_DOUBLE_PRECISION, &
            MPI_MIN, MPI_COMM_CALC, IERROR)
       minprecip=x_out
       if(minprecip<-10)then
          if(me==0)write(*,*)'WARNING: found negative precipitations. '&
             ,'set precipitations to zero!',minprecip
          surface_precip = 0.0
       else
          surface_precip = max(0.0, & ! get only the variation. mm ->m/s
             (surface_precip - surface_precip_old))*0.001/(METSTEP*3600)
       end if
       surface_precip_old = buff ! Accumulated precipitation
     end if

    !write(*,*)'precip ',nrec,Nhh,surface_precip(5,5),convective_precip(5,5),surface_precip_old(5,5)
    if(WRF_MET_CORRECTIONS) then
      if(found_wrf_bucket)then
        !wrf "bucket" definition for surface precipitation:
        !surface_precip = I_RAINNC*bucket + RAINNC + I_RAINC*bucket + RAINC
        surface_precip=surface_precip+irainnc*wrf_bucket+irainc*wrf_bucket
      end if
      buff=surface_precip !save to save in old below

      !must first check that precipitation is increasing. At some dates WRF maybe restarted!
      minprecip=minval(surface_precip(1:limax,1:ljmax) - surface_precip_old(1:limax,1:ljmax))

      CALL MPI_ALLREDUCE(minprecip, x_out, 1,MPI_DOUBLE_PRECISION, &
             MPI_MIN, MPI_COMM_CALC, IERROR)
      minprecip=x_out
      if(minprecip<-10)then
        if(me==0)write(*,*)'WARNING: found negative precipitations. ',&
          ' set precipitations to zero!',minprecip
        surface_precip = 0.0
      else
        surface_precip = max(0.0,(surface_precip - surface_precip_old))*0.001/(METSTEP*3600)! get only the variation. mm ->m/s
      end if
      surface_precip_old = buff ! Accumulated rain in WRF
      ice_nwp = ice_nwp*100!flag->%
      !smooth qrain, because it is instantaneous but rain may move
      do k=1,kmax_mid
        call smoosp(rain(1,1,k,nr),0.0,1.0E10)
      end do
    end if

    ! if available, will use cloudwater to determine the height of release
    ! NB: array cw_met only used here
    if(nr==2)cw_met(:,:,:,1)=cw_met(:,:,:,2)!save previous value

    ix = ix_cw_met
    call Getmeteofield(meteoname,met(ix)%name,nrec,met(ix)%dim,unit,&
           met(ix)%validity,met(ix)%field,needed=met(ix)%needed,found=met(ix)%found)

    if(foundrain)then
      if(write_now)&
        write(*,*)'release profile for 3D precipitations derived from QRAIN'
      do j=1,ljmax
        do i=1,limax
          !note that there is much noise in surface precipitations, because it is the difference of two large  number (1E5)
          !values smaller than 0.01 mm are meaningless
          if(surface_precip(i,j)>1E-9 .and. &
             rain(i,j,kmax_mid,1)+rain(i,j,kmax_mid,nr)>1E-12)then
            do k=1,kmax_mid
              pr(i,j,k)=min(surface_precip(i,j),&
                            surface_precip(i,j)*(rain(i,j,k,1)+rain(i,j,k,nr))&
                           /(rain(i,j,kmax_mid,1)+rain(i,j,kmax_mid,nr)))
              pr(i,j,k)=pr(i,j,k)*METSTEP*3600.0*1000.0! and m/s->mm/METSTEP
            end do
          elseif(surface_precip(i,j)>1.0E-8 )then
            ! surface precipitation but no QRAIN found.
            ! Should not occur too often. use humidity method.
            pr(i,j,KMAX_MID)= surface_precip(i,j)*METSTEP*3600.0*1000.0! guarantees precip at surface
            do k=1,KMAX_MID-1
              !convert from potential temperature into absolute temperature
              temperature = th(i,j,k,nr) &
                  *exp(KAPPA*log((A_mid(k)+B_mid(k)*ps(i,j,nr))*1.e-5))!Pa, Ps in Pa here
              !saturation water pressure
              swp=611.2*exp(17.67*(temperature-273.15)/(temperature-29.65))
              !water pressure
              wp=q(i,j,k,nr)*(A_mid(k)+B_mid(k)*ps(i,j,nr))/0.622!Ps in Pa here
              relh2=wp/swp
              !convert from potential temperature into absolute temperature
              !Ps  in Pa here
              temperature = th(i,j,k,1) &
                  *exp(KAPPA*log((A_mid(k)+B_mid(k)*ps(i,j,1))*1.e-5))!Ps in Pa here
              !saturation water pressure
              swp=611.2*exp(17.67*(temperature-273.15)/(temperature-29.65))
              !water pressure
              wp=q(i,j,k,1)*(A_mid(k)+B_mid(k)*ps(i,j,1))/0.622!Ps in Pa here
              relh1=wp/swp
              if(relh1>RH_THRESHOLD.or.relh2>RH_THRESHOLD)then
                !fill the column up to this level with constant precip
                do kk=k,KMAX_MID-1
                  pr(i,j,kk)= surface_precip(i,j)*METSTEP*3600.0*1000.0!m/s->mm/METSTEP
                end do
                exit
              else
                pr(i,j,k)=0.0
              end if
            end do
          else
            do k=1,kmax_mid
              pr(i,j,k)=0.0
            end do
          end if
        end do
      end do

    elseif(foundcloudwater)then

      if(write_now)&
        write(*,*)'release height for 3D precipitations derived from cloudwater'
      if(MasterProc.and.first_call)&
        write(IO_LOG,*)"3D precipitations: derived from 2D and cloudwater"
      if(nr==1)cw_met(:,:,:,2)=cw_met(:,:,:,nr)!so that nr=2 also is defined
      do j=1,ljmax
        do i=1,limax
          pr(i,j,KMAX_MID)= surface_precip(i,j)*METSTEP*3600.0*1000.0!guarantees precip at surface
          do k=1,KMAX_MID-1
            if(cw_met(i,j,k,2)+cw_met(i,j,k,1)>CW_THRESHOLD)then
              !fill the column up to this level with constant precip
              do kk=k,KMAX_MID-1
                pr(i,j,kk)= surface_precip(i,j)*METSTEP*3600.0*1000.0!from m/s to mm/METSTEP
              end do
              exit
            else
              pr(i,j,k)=0.0
            end if
          end do
        end do
      end do

    else
      !will use RH to determine the height of release (less accurate than cloudwater)
      if(write_now)write(*,*)'release height for 3D precipitations derived from humidity'
      if(MasterProc.and.first_call)write(unit=IO_LOG,fmt="(a)")&
           "3D precipitations: derived from 2D and humidity"
      do j=1,ljmax
        do i=1,limax
          pr(i,j,KMAX_MID)= surface_precip(i,j)*METSTEP*3600.0*1000.0!guarantees precip at surface
          do k=1,KMAX_MID-1
            !convert from potential temperature into absolute temperature
            temperature = th(i,j,k,nr)*exp(KAPPA*log((A_mid(k) + B_mid(k)*ps(i,j,nr))*1.e-5))!Pa, Ps in Pa here
            !saturation water pressure
            swp=611.2*exp(17.67*(temperature-273.15)/(temperature-29.65))
            !water pressure
            wp=q(i,j,k,nr)*(A_mid(k) + B_mid(k)*ps(i,j,nr))/0.622!Ps in Pa here
            relh2=wp/swp
            !convert from potential temperature into absolute temperature
            !Ps  in Pa here
            temperature = th(i,j,k,1)* &
                exp(KAPPA*log((A_mid(k) + B_mid(k)*ps(i,j,1))*1.e-5))!Ps in Pa here
            !saturation water pressure
            swp=611.2*exp(17.67*(temperature-273.15)/(temperature-29.65))
            !water pressure
            wp=q(i,j,k,1)*(A_mid(k) + B_mid(k)*ps(i,j,1))/0.622!Ps in Pa here
            relh1=wp/swp
            if(relh1>RH_THRESHOLD.or.relh2>RH_THRESHOLD)then
              !fill the column up to this level with constant precip
              do kk=k,KMAX_MID-1
                pr(i,j,kk)= surface_precip(i,j)*METSTEP*3600.0*1000.0!METSTEP hours and m->mm
              end do
              exit
            else
              pr(i,j,k)=0.0
            end if
          end do
        end do
      end do
    end if
  end if ! .not.foundprecip
  pr=max(0.0,pr)/(3600.0*METSTEP) ! positive precipitation in mm/s

  ! surface precipitation, mm/hr
  ! NB: surface_precip is different than the one read directly from the
  ! metfile  (which has different units, and is the sum of the 2D
  ! large_scale_precipitations+convective_precipitations)
  surface_precip(:,:) = pr(:,:,KMAX_MID) * 3600 !mm/s -> mm/hr

  if(USES%CONVECTION)then
    cnvuf=max(0.0,cnvuf)      !no negative upward fluxes
    cnvuf(:,:,KMAX_BND)=0.0   !no flux through surface
    cnvuf(:,:,1)=0.0          !no flux through top
    cnvdf=min(0.0,cnvdf)      !no positive downward fluxes
    cnvdf(:,:,KMAX_BND)=0.0   !no flux through surface
    cnvdf(:,:,1)=0.0          !no flux through top

    ! Sometimes the NWP calculated fluxes are too high.
    ! can be scaled via config:
    if( abs( CONVECTION_FACTOR - 1 ) > 0.001 ) then
     cnvuf(:,:,:) = cnvuf(:,:,:) * CONVECTION_FACTOR
     cnvdf(:,:,:) = cnvdf(:,:,:) * CONVECTION_FACTOR
    end if
  end if

  ! Kz from meteo
  if(NWP_Kz) then
    if(foundKz_met)then
      Kz_met=max(0.0,Kz_met)  ! only positive Kz
    elseif(MasterProc.and.first_call)then
      write(*,*)' WARNING: Kz will be derived in model '
    end if
  end if

  !==============    2D fields (surface) (i,j)   ============================
  ndim=2

  if(foundrh2m)then
    rh2m(:,:,nr) = 0.01 * rh2m(:,:,nr)  ! Convert from % to fraction
  else
    if(MasterProc.and.first_call)&
      write(*,*)'WARNING: relative_humidity_2m not found'
    rh2m(:,:,nr) = -999.9  ! ?
  end if

  if(.not. WRF_MET_CORRECTIONS .and. foundsdepth)then
     !IFS defines snowdepth in units of water equivalent
     sdepth(:,:,nr) = sdepth(:,:,nr)*5 !crude conversion from water equivalent into meters of snow
  endif

  if(WRF_MET_CORRECTIONS)then
    ! flux defined with opposite signs
    fh(:,:,nr)=-fh(:,:,nr)
    fl(:,:,nr)=-fl(:,:,nr)
    sst=max(273.0,sst)!WRF set land sst to zero
    ! rough conversion wrf->SMI . Can be improved!
    SoilWater_uppr(:,:,nr)=(SoilWater_uppr(:,:,nr)-0.05)*3
    SoilWater_deep(:,:,nr)=(SoilWater_deep(:,:,nr)-0.05)*3
  end if

  if(LANDIFY_MET)then
    call landify(t2_nwp(:,:,nr),"t2nwp")
    call landify(rh2m(:,:,nr),"rh2m")
    call landify(fh(:,:,nr),"fh")
    call landify(fl(:,:,nr),"fl")
    if(foundtau)call landify(tau(:,:,nr),"tau")
  end if

  if(met(ix_fh)%validity=='averaged')fh(:,:,1)=fh(:,:,nr)

  if(met(ix_fl)%validity=='averaged')fl(:,:,1)=fl(:,:,nr)

  if(foundtau)then
    tau=max(0.0,tau)
    if(met(ix_tau)%validity=='averaged')tau(:,:,1)=tau(:,:,nr)
  else
    !     For WRF we get u*, not tau. Since it seems better to
    !     interpolate tau than u*  between time-steps we convert
    if(write_now)write(*,*)' tau derived from ustar_nwp'
    namefield=met(ix_ustar_nwp)%name
    call Getmeteofield(meteoname,namefield,nrec,ndim,unit,validity,&
        ustar_nwp(:,:),needed=.true.,found=foundustar)
    if(LANDIFY_MET) call landify(ustar_nwp(:,:),"ustar")
    !Ps in Pa here
    forall(i=1:limax,j=1:ljmax)
      rho_surf(i,j) = ps(i,j,nr)/(RGAS_KG*t2_nwp(i,j,nr) )
      tau(i,j,nr)   = ustar_nwp(i,j)*ustar_nwp(i,j)*rho_surf(i,j)
    end forall
  end if

  if(.not.foundSST.and.write_now)&
    write(*,*)' WARNING: sea_surface_temperature not found '

  ! Soil water fields. Somewhat tricky.
  ! Ideal is soil moisture index, available from IFS, = (SW-PWP)/(FC-PWP)
  ! Otherwise m3/m3 or m units are converted
  !
  ! Start with shallow

  call CheckStop(USES%DUST.and..not.USES%SOILWATER,"Inconsistent SM, DUST")

  if ( PBL%HmixMethod == "NWP") then
    foundHmix=.false.
    do isw = 1, size(possible_HmixNames)
      namefield=possible_HmixNames(isw)
      if(first_call.and.MasterProc) &
            write(*,*) "Met_mod: HMIX search ",isw,trim(namefield)
      call Getmeteofield(meteoname,namefield,nrec,ndim,unit,validity,&
               pbl_nwp(:,:,nr),found=foundHmix)
      if(foundHmix) then ! found
          exit
      end if
    end do
  end if ! PBL Hmix

  if(USES%SOILWATER) then
    ! Soil water fields. Somewhat tricky.
    ! Ideal is soil moisture index, available from IFS, = (SW-PWP)/(FC-PWP)
    ! Otherwise m3/m3 or m units are converted
    !
    ! Start with shallow
    if(.not.foundSoilWater_uppr) then
      foundSMI1=.false.
      do isw = 1, size(possible_soilwater_uppr)
        namefield=possible_soilwater_uppr(isw)
        if((DEBUG_SOILWATER.or.first_call).and.MasterProc) &
               write(*,*) "Met_mod: soil water search ",isw,trim(namefield)
        call Getmeteofield(meteoname,namefield,nrec,ndim,unit,validity,&
               SoilWater_uppr(:,:,nr),found=foundSoilWater_uppr)
        if(foundSoilWater_uppr) then ! found
          foundSMI1=(index(namefield,"SMI")>0)
          exit
        end if
      end do
      if(foundSMI1.and.MasterProc.and.first_call) &  ! = 1st call
        call PrintLog("Met: found SMI1:"//trim(namefield))
      if(foundSoilWater_uppr.and.trim(unit)=="m") SoilWaterSource="PARLAM"
    end if ! upper

    if(.not.foundSoilWater_deep) then  !just deep here
      foundSMI3=.false.
      do isw = 1, size(possible_soilwater_deep)
        namefield=possible_soilwater_deep(isw)
        if(GRID=="HIRHAM") then
          if(MasterProc.and.first_call)write(*,*) " Rename soil water in HIRHAM"
          namefield='soil_water_second_layer'
        end if
        if(MasterProc.and.first_call) write(*,*) "Met_mod: ', &
          'deep soil water search ", isw, trim(namefield)
        call Getmeteofield(meteoname,namefield,nrec,ndim,unit,validity,&
            SoilWater_deep(:,:,nr),found=foundSoilWater_deep)
        if(foundSoilWater_deep) then ! found
          foundSMI3=(index(namefield,"SMI")>0)
          if(.not.foundSMI3) &  ! = 1st call
            call PrintLog("Met: found SMI3:"//trim(namefield), MasterProc)
          exit
        end if
      end do
      if(foundSoilWater_deep ) then
        if(trim(unit)=="m") then  ! PARLAM has metres of water
          SoilWaterSource = "PARLAM"
        elseif(unit(1:5)=='m3/m3')then
          SoilWaterSource = "IFS"
        end if
      end if !found deep_soil_water_content
    end if !

    if(SoilWaterSource == "IFS")then
      if(first_call)then
        !needed for transforming IFS soil water
        call ReadField_CDF(SoilTypesFile,'pwp',pwp, 1,&
            interpol='conservative',needed=.true.,UnDef=-999.,debug_flag=.false.)
        call ReadField_CDF(SoilTypesFile,'fc',fc, 1,&
            interpol='conservative',needed=.true.,UnDef=-999.,debug_flag=.false.)

        ! landify(x,intxt,xmin,xmax,wfmin,xmask)
        ! We use a global mask for water_fraction < 100%, but set wfmin to 1.0
        ! to allow all grids with some land to be processed
        ! Fc and PWP should be above zero and  below 1, let's use 0.8

        call landify( pwp(:,:), " PWP ", &
            0.0, 0.8, 1.0, water_fraction < 1.0 ) ! mask for where there is land
        call landify( fc(:,:), " FC  ", &
            0.0, 0.8, 1.0, water_fraction < 1.0 ) ! mask for where there is land
        do i = 1, limax
          do j = 1, ljmax
            if( fc(i,j) > pwp(i,j) ) then ! Land values
              tmpmax = -0.99 * pwp(i,j)/(fc(i,j)-pwp(i,j) )
              SoilWater_uppr(i,j,nr) = max( tmpmax, SoilWater_uppr(i,j,nr) )
            else
              SoilWater_uppr(i,j,nr) = -999.  ! NOT NEEDED????
            end if
          end do
        end do
      end if
    end if

    if(foundSMI3.or.foundSoilWater_deep)then
      if ( water_frac_set ) then  ! smooth the SoilWater values:

        ! If NWP thinks this is a sea-square, but we anyway have land,
        ! the soil moisture might be very strange.  We search neighbouring
        ! grids and make a land-weighted mean SW
        ! Skip on 1st numt, since water fraction set a little later. No harm done...
        ! changed landify routine to accept water_fraction as mask. Should
        ! works almost the same as code below did.
        ! Should move later also, after other units converted to SMI
        ! NB  Some grid squares in EECCA have water cover of 99.998
        call landify( SoilWater_deep(:,:,nr), "SMI_DEEP", &
            0.0, 1.0, 1.0, water_fraction < 1.0 )
        ! Allow some negative SMI for upper levels
        call landify( SoilWater_uppr(:,:,nr), "SMI_UPPR", &
            -1.0, 1.0, 1.0, water_fraction < 1.0 )
      else ! water_frac not set yet
        call CheckStop("ERROR, Met_mod: SMD not set!!  here"  )
      end if ! water_frac_set test
    end if !

    ! SMI = (SW-PWP)/((FC-PWP), therefore min SMI value should be -PWP/(FC-PWP)
    ! Let's use 99% of this:
    if ( SoilWaterSource == "IFS") then !MAR2013
      do i = 1, limax
        do j = 1, ljmax
          if( fc(i,j) > pwp(i,j) ) then ! Land values
            tmpmax = -0.99 * pwp(i,j)/(fc(i,j)-pwp(i,j) )
            SoilWater_uppr(i,j,nr) = max( tmpmax, SoilWater_uppr(i,j,nr) )
          else
            SoilWater_uppr(i,j,nr) = -999.  ! NOT NEEDED????
          end if
        end do
      end do
    end if !MAR2013 test

    ! PARLAM SMI calc: hard-coded FC value of 0.02 (cm)
    if( SoilWaterSource == "PARLAM")then
      !SoilMax = 0.02
      SoilWater_deep(:,:,nr) = SoilWater_deep(:,:,nr) / 0.02 !SoilMax
      SoilWater_uppr(:,:,nr) = SoilWater_uppr(:,:,nr) / 0.02 !SoilMax
    end if

    ! We should now have SMI regardless of soil water data source.
    ! We restrict this to be in range 0 --- 1 for deep soil water.
    ! For upper-soil water, we allow some negative, since evaporation
    ! can dry the soil bellow the PWP.
    SoilWater_deep(:,:,nr) = max(0.0, SoilWater_deep(:,:,nr) )
    SoilWater_uppr(:,:,nr) = min(1.0, SoilWater_uppr(:,:,nr))
    SoilWater_deep(:,:,nr) = min(1.0, SoilWater_deep(:,:,nr) )

  end if ! USES%SOILWATER

  !========================================

  if(.not.foundsdepth.and.write_now)&
    write(*,*)dtxt//' WARNING: snow_depth not found '

  if(.not.foundice.and.write_now)&
    write(*,*)dtxt//' WARNING: ice_nwp coverage (%) not found '


  if(foundws10_met)then
    namefield='v10' !second component of ws_10m
    call Getmeteofield(meteoname,namefield,nrec,ndim,unit,validity,&
        buff(:,:),needed=.false.,found=foundws10_met)
    if(.not.foundws10_met)then
      namefield='V10' !second component of ws_10m
      call Getmeteofield(meteoname,namefield,nrec,ndim,unit,validity,&
           buff(:,:),found=foundws10_met)
    end if
    if(foundws10_met)then
      if(write_now)write(*,*)dtxt//' found v component of 10m wind '
      ws_10m(:,:,nr)=sqrt(ws_10m(:,:,nr)**2+buff(:,:)**2)
      if(LANDIFY_MET) call landify(ws_10m(:,:,nr),"WS10")
    end if
  end if


  ! derive the meteorological parameters from the basic parameters
  ! read from field-files.
  do j = 1,ljmax
    do i = 1,limax
      p1 = A_bnd(KMAX_BND)+B_bnd(KMAX_BND)*ps(i,j,nr)

      exf1(KMAX_BND) = CP * Exner_nd(p1)

      z_bnd(i,j,KMAX_BND) = 0.0

      do k = KMAX_MID,1,-1
        ! eddy diffusivity in the surface-layer follows the
        ! formulation usedin the nwp-model which is based on
        ! Louis (1979), (see mc7e.f).
        ! exner-function of the half-layers

        p1 = A_bnd(k)+B_bnd(k)*ps(i,j,nr)
        exf1(k) = CP * Exner_nd( p1 )

        ! exner-function of the full-levels
        p2 = A_mid(k)+B_mid(k)*ps(i,j,nr)
        exf2(k) = CP * Exner_nd(p2)

        ! height of the half-layers
        z_bnd(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
            (exf1(k+1) - exf1(k)))/GRAV

        ! height of the full levels.
        z_mid(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
            (exf1(k+1) - exf2(k)))/GRAV

        ! derive density from pressure and heights, so that they are consistent
        roa(i,j,k,nr) = (dA(k)+dB(k)*ps(i,j,nr))&
                       /(GRAV*(z_bnd(i,j,k)-z_bnd(i,j,k+1)))

        end do  ! k
     end do
  end do

  ! set etadot equal to zero at the top and bottom of atmosphere.
  Etadot(:,:,KMAX_BND,nr)=0.0
  Etadot(:,:,1,nr)=0.0
  if( .not.met(ix_Etadot)%found)then
    if(write_now)write(*,*)'WARNING: Etadot derived from horizontal winds '
    ! sdot derived from divergence=0 principle
    do j = 1,ljmax
      do i = 1,limax
        Ps_extended(i,j) = Ps(i,j,nr)
      end do
    end do
    ! Get Ps at edges from neighbors
    ! we reuse usnd, vsnd etc
    if (neighbor(EAST) .ne. NOPROC) then
      usnd(:,1) = ps(limax,:,nr)
      CALL MPI_ISEND( usnd, 8*LJMAX, MPI_BYTE,  &
           neighbor(EAST), MSG_WEST2, MPI_COMM_CALC, request_e, IERROR)
    end if
    if (neighbor(NORTH) .ne. NOPROC) then
      vsnd(:,1) = ps(:,ljmax,nr)
      CALL MPI_ISEND( vsnd , 8*LIMAX, MPI_BYTE,  &
           neighbor(NORTH), MSG_SOUTH2, MPI_COMM_CALC, request_n, IERROR)
    end if
    ! receive from WEST neighbor if any
    if (neighbor(WEST) .ne. NOPROC) then
      CALL MPI_RECV( urcv, 8*LJMAX, MPI_BYTE, &
           neighbor(WEST), MSG_WEST2, MPI_COMM_CALC, MPISTATUS, IERROR)
      Ps_extended(0,1:ljmax) = urcv(1:ljmax,1)
    else
      Ps_extended(0,1:ljmax) = Ps_extended(1,1:ljmax)
    end if
    ! receive from SOUTH neighbor if any
    if (neighbor(SOUTH) .ne. NOPROC) then
      CALL MPI_RECV( vrcv, 8*LIMAX, MPI_BYTE,  &
           neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_CALC, MPISTATUS, IERROR)
      Ps_extended(1:limax,0) = vrcv(1:limax,1)
    else
      Ps_extended(1:limax,0) = Ps_extended(1:limax,1)
    end if
    if (neighbor(WEST) .ne. NOPROC) then
      usnd(:,2) = ps(1,:,nr)
      CALL MPI_ISEND( usnd(1,2), 8*LJMAX, MPI_BYTE,  &
           neighbor(WEST), MSG_WEST2, MPI_COMM_CALC, request_w, IERROR)
    end if
    if (neighbor(SOUTH) .ne. NOPROC) then
      vsnd(:,2) = ps(:,1,nr)
      CALL MPI_ISEND( vsnd(1,2) , 8*LIMAX, MPI_BYTE,  &
           neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_CALC, request_s, IERROR)
    end if

    ! receive from EAST neighbor if any
    if (neighbor(EAST) .ne. NOPROC) then
      CALL MPI_RECV( urcv, 8*LJMAX, MPI_BYTE, &
           neighbor(EAST), MSG_WEST2, MPI_COMM_CALC, MPISTATUS, IERROR)
      Ps_extended(limax+1,1:ljmax) = urcv(1:ljmax,1)
    else
      Ps_extended(limax+1,1:ljmax) = Ps_extended(limax,1:ljmax)
    end if
    ! receive from NORTH neighbor if any
    if (neighbor(NORTH) .ne. NOPROC) then
      CALL MPI_RECV( vrcv, 8*LIMAX, MPI_BYTE,  &
           neighbor(NORTH), MSG_SOUTH2, MPI_COMM_CALC, MPISTATUS, IERROR)
      Ps_extended(1:limax,ljmax+1) = vrcv(1:limax,1)
    else
      Ps_extended(1:limax,ljmax+1) = Ps_extended(1:limax,ljmax)
    end if

    if (neighbor(EAST) .ne. NOPROC) then
      CALL MPI_WAIT(request_e, MPISTATUS, IERROR)
    end if

    if (neighbor(NORTH) .ne. NOPROC) then
      CALL MPI_WAIT(request_n, MPISTATUS, IERROR)
    end if
    if (neighbor(WEST) .ne. NOPROC) then
      CALL MPI_WAIT(request_w, MPISTATUS, IERROR)
    end if

    if (neighbor(SOUTH) .ne. NOPROC) then
      CALL MPI_WAIT(request_s, MPISTATUS, IERROR)
    end if
! Eta=A/Pref+B; see
! http://www.ecmwf.int/research/ifsdocs/DYNAMICS/Chap2_Discretization3.html#959545

    !(note that u_xmj and v_xmi have already been divided by xm here)
    dB_sum=1.0/(B_bnd(KMAX_MID+1)-B_bnd(1))!normalisation factor for dB (should be one if entire atmosphere is included)

    do j = 1,ljmax
      do i = 1,limax
        Pmid=Ps_extended(i,j)! without "-PT"
        !surface pressure at gridcell boundaries
        Pu1=0.5*(Ps_extended(i-1,j)+Ps_extended(i,j))
        Pu2=0.5*(Ps_extended(i+1,j)+Ps_extended(i,j))
        Pv1=0.5*(Ps_extended(i,j-1)+Ps_extended(i,j))
        Pv2=0.5*(Ps_extended(i,j+1)+Ps_extended(i,j))

        sumdiv=0.0
        do k=1,KMAX_MID
          divk(k)=((u_xmj(i,j,k,nr)  *(dA(k)+dB(k)*Pu2)   &
                   -u_xmj(i-1,j,k,nr)*(dA(k)+dB(k)*Pu1))  &
                  +(v_xmi(i,j,k,nr)  *(dA(k)+dB(k)*Pv2)   &
                   -v_xmi(i,j-1,k,nr)*(dA(k)+dB(k)*Pv1))) &
                  * xm2(i,j)/ GRIDWIDTH_M
          sumdiv=sumdiv+divk(k)
        end do

        Etadot(i,j,KMAX_MID+1,nr)=0.0
        do k=KMAX_MID,1,-1
          Etadot(i,j,k,nr)=Etadot(i,j,k+1,nr)-dB(k)*dB_sum*sumdiv+divk(k)
          Etadot(i,j,k+1,nr)=Etadot(i,j,k+1,nr)*(dA(k)/Pref+dB(k))/(dA(k)+dB(k)*Pmid)
        end do
        Etadot(i,j,1,nr)=0.0! Is zero anyway from relations above
      end do
    end do

    if(MANUAL_GRID)then
      !around the north pole the interpolation of wind fields is not accurate,
      !and we cannot rely on continuity to derive the vertical winds.
      !Brutal: set Etadot to 0!
      do k=1,KMAX_MID,1
        do j = 1,ljmax
          do i = 1,limax
             if(glat(i,j)>88.0)Etadot(i,j,k,nr)=0.0
          end do
        end do
      end do
    end if

  end if

  if(met(ix_Etadot)%found)then
    call CheckStop(.not.USE_EtaCOORDINATES,&
        "Conflict: requested etadot, but does not use eta coordinates")
    !convert from mid values to boundary values
    if(write_now)write(*,*)'interpolating etadot from mid to boundary levels'
    do k = KMAX_MID,2,-1
      do j = 1,ljmax
        do i = 1,limax
          Etadot(i,j,k,nr) = Etadot(i,j,k-1,nr)       &
            + (Etadot(i,j,k,nr)-Etadot(i,j,k-1,nr))   &
            * (Eta_bnd(k)-Eta_mid(k-1))               &
            / (Eta_mid(k)-Eta_mid(k-1))
        end do
      end do
    end do
    Etadot(:,:,1,nr)=0.0!no exchanges above top (should not be useed anyway)
  end if

  call met_derived(nr) !compute derived meteo fields used in BLPhysics

  call BLPhysics()

  call met_derived(1) !compute derived meteo fields for nr=1 "now"

  !windspeed of neighbor subdomains at edges (used for advection)
  !It the windspeed divided by xm which must be used here.

  ! send to WEST neighbor if any
  if (neighbor(WEST) .ne. NOPROC) then
    if(neighbor(WEST) .ne. me)then
      buf_uw(:,:) = u_xmj(1,:,:,nr)
      CALL MPI_ISEND(buf_uw, 8*LJMAX*KMAX_MID, MPI_BYTE, &
           neighbor(WEST), MSG_EAST2, MPI_COMM_CALC, request_w, IERROR)
    else
    ! cyclic grid: own neighbor
      ue(:,:,nr) = u_xmj(1,:,:,nr)
    end if
  end if

  ! send to EAST neighbor if any
  if (neighbor(EAST) .ne. NOPROC) then
    if (neighbor(EAST) .ne. me) then
      buf_ue(:,:) = u_xmj(limax-1,:,:,nr)
      CALL MPI_ISEND(buf_ue, 8*LJMAX*KMAX_MID, MPI_BYTE, &
          neighbor(EAST), MSG_WEST2, MPI_COMM_CALC, request_e, IERROR)
    else
      ! cyclic grid: own neighbor
      uw(:,:,nr) = u_xmj(limax-1,:,:,nr)
    end if
  end if

  ! send to SOUTH neighbor if any
  if (neighbor(SOUTH) .ne. NOPROC) then
    buf_vs(:,:) = v_xmi(:,1,:,nr)
    CALL MPI_ISEND(buf_vs, 8*LIMAX*KMAX_MID, MPI_BYTE, &
          neighbor(SOUTH), MSG_NORTH2, MPI_COMM_CALC, request_s, IERROR)
  end if

  ! send to NORTH neighbor if any
  if (neighbor(NORTH) .ne. NOPROC) then
    buf_vn(:,:) = v_xmi(:,ljmax-1,:,nr)
    CALL MPI_ISEND(buf_vn, 8*LIMAX*KMAX_MID, MPI_BYTE, &
          neighbor(NORTH), MSG_SOUTH2, MPI_COMM_CALC, request_n, IERROR)
  end if

  ! receive from EAST neighbor if any
  if (neighbor(EAST) .ne. NOPROC .and. neighbor(EAST) .ne. me) then
    CALL MPI_RECV(ue(1,1,nr), 8*LJMAX*KMAX_MID, MPI_BYTE, &
        neighbor(EAST), MSG_EAST2, MPI_COMM_CALC, MPISTATUS, IERROR)
  end if

  ! receive from WEST neighbor if any
  if (neighbor(WEST) .ne. NOPROC .and. neighbor(WEST) .ne. me) then
    CALL MPI_RECV(uw(1,1,nr), 8*LJMAX*KMAX_MID, MPI_BYTE, &
         neighbor(WEST), MSG_WEST2, MPI_COMM_CALC, MPISTATUS, IERROR)
  end if

  ! receive from NORTH neighbor if any
  if (neighbor(NORTH) .ne. NOPROC) then
    CALL MPI_RECV(vn(1,1,nr), 8*LIMAX*KMAX_MID, MPI_BYTE, &
         neighbor(NORTH), MSG_NORTH2, MPI_COMM_CALC, MPISTATUS, IERROR)
  end if

  ! receive from SOUTH neighbor if any
  if (neighbor(SOUTH) .ne. NOPROC) then
    CALL MPI_RECV(vs(1,1,nr), 8*LIMAX*KMAX_MID, MPI_BYTE, &
        neighbor(SOUTH), MSG_SOUTH2, MPI_COMM_CALC, MPISTATUS, IERROR)
  end if

  if (neighbor(EAST) .ne. NOPROC .and. neighbor(EAST) .ne. me) then
    CALL MPI_WAIT(request_e, MPISTATUS, IERROR)
  end if
  if (neighbor(WEST) .ne. NOPROC .and. neighbor(WEST) .ne. me) then
    CALL MPI_WAIT(request_w, MPISTATUS, IERROR)
  end if
  if (neighbor(NORTH) .ne. NOPROC) then
    CALL MPI_WAIT(request_n, MPISTATUS, IERROR)
  end if
  if (neighbor(SOUTH) .ne. NOPROC) then
    CALL MPI_WAIT(request_s, MPISTATUS, IERROR)
  end if

  if(USE_FASTJ)then
    !compute photolysis rates from FastJ
    if(nr==2)rcphot_3D(:,:,:,:,1)=rcphot_3D(:,:,:,:,2)
    do j = 1,ljmax
      do i = 1,limax
        call setup_phot_fastj(i,j,INFO,nr)
      end do
    end do
  end if

  if(first_call.and.next_inptime%hour<nhour_first)&
    nrec=0 ! not yet reached first available time in meteo file

  first_call=.false.
end subroutine Meteoread

subroutine metfieldint
  ! Forward linear stepping of the meteorological
  ! fields read or derived every 3 hours.

  implicit none

  real :: div
  integer :: ix

  if (mod(step_main,nmax) > 0) then
    div = 1./real(nmax-(mod(step_main,nmax)-1))
    do ix=1,Nmetfields
      if(met(ix)%time_interpolate)then
        !if(DEBUG%MET) print *, 'METINTERP ',me,trim(met(ix)%name)
        !if(DEBUG%MET) print *, 'METSIZE   ',me,size(met(ix)%field,dim=4)

        if(me==0.and.DEBUG%MET)&
          write(*,*)'interpolating in time ',ix,met(ix)%name

        met(ix)%field(:,:,:,1)    = met(ix)%field(:,:,:,1)                 &
        + (met(ix)%field(:,:,:,2) - met(ix)%field(:,:,:,1))*div

      end if
    end do

  else
    do ix=1,Nmetfields
      if(met(ix)%time_interpolate)then
               met(ix)%field(:,:,:,1)    = met(ix)%field(:,:,:,2)
      end if
    end do
  end if

  call met_derived(1) !update derived meteo fields
end subroutine metfieldint

subroutine met_derived(nt)
  ! Calculates fields derived from meteofields.
  ! The interpolation in time is done for the meteofields and the
  ! fields here are derived from the interpolated fields after
  ! each interpolation (i.e. every dt_advec).
  ! CPU costly fields (those with special functions like log )
  ! can be computed in MeteoRead  only once every METSTEP and interpolated
  ! in metint.

  ! horizontal wind speed (averaged over the four edges)
  ! Note that u_xmj and v_xmi are wind velocities divided by xm
  ! At present u_ref is defined at KMAX_MID

  implicit none
  integer, intent(in) :: nt  ! set to 1 from metint or nr from matvar
  integer ::i,j, k
  logical :: DEBUG_DERIV = .false.

  do k = 1, KMAX_MID
    do j = 1,ljmax
      do i = 1,limax
        u_mid(i,j,k) = 0.5*(u_xmj(i,j,k,nt)  *xm_j(i,j)   &
                           +u_xmj(i-1,j,k,nt)*xm_j(i-1,j))
        v_mid(i,j,k) = 0.5*(v_xmi(i,j-1,k,nt)*xm_i(i,j-1)  &
                           +v_xmi(i,j,k,nt)  *xm_i(i,j))
      end do
    end do
  end do !k
  do j = 1,ljmax
    do i = 1,limax
      u_ref(i,j)= sqrt( u_mid(i,j,KMAX_MID)**2 + v_mid(i,j,KMAX_MID)**2 )
    end do
  end do
  if(LANDIFY_MET) &
    call landify(u_ref(:,:),"u_ref")


  ! Tmp ustar solution. May need re-consideration for MM5 etc., but
  ! basic principal should be that fm is interpolated with time, and
  ! ustar derived from this.

  forall( i=1:limax, j=1:ljmax )
    rho_surf(i,j) = ps(i,j,nt)/(RGAS_KG * t2_nwp(i,j,nt) )
  end forall

  ! update z_bnd
! z_bnd(:,:,KMAX_MID+1)=0.0 !should never change
  do k = KMAX_MID,1,-1
    do j = 1,ljmax
      do i = 1,limax
        z_bnd(i,j,k)=z_bnd(i,j,k+1)+(dA(k)+dB(k)*ps(i,j,nt))/(roa(i,j,k,nt)*GRAV)
      end do
    end do
  end do

  ! always use tau, since ustar_nwp is not interpolated in time (in metfieldint)
  forall( i=1:limax, j=1:ljmax )
    ustar_nwp(i,j)   = sqrt( tau(i,j,nt)/rho_surf(i,j) )
  end forall

  ! we limit u* to a physically plausible value over land
  ! to prevent numerical problems, and to account for enhanced
  ! mixing which is usually found over real terrain

  where ( mainly_sea )
     ustar_nwp = max( ustar_nwp, 1.0e-5 )
  elsewhere
     ustar_nwp = max( ustar_nwp, MIN_USTAR_LAND )
  end where

  forall( i=1:limax, j=1:ljmax )
    invL_nwp(i,j)  = KARMAN * GRAV * fh(i,j,nt) & ! - disliked by gfortran
            / (CP*rho_surf(i,j) * ustar_nwp(i,j)**3 * t2_nwp(i,j,nt) )
  end forall

  where ( invL_nwp(:,:) < -1.0 )
    invL_nwp(:,:)  = -1.0
  else where ( invL_nwp(:,:) > 1.0 )
    invL_nwp(:,:)  = 1.0
  end where


  if ( DEBUG_DERIV .and. debug_proc ) then
    i = debug_iloc
    j = debug_jloc
    write(*,*) "MET_DERIV DONE ", me, nt, ustar_nwp(i,j), rho_surf(i,j), &
        fh(i,j,nt), invL_nwp(i,j)
  end if
end subroutine met_derived

subroutine MetModel_LandUse(callnum)
  ! Reads parameterfields from file
  ! reading surface roughness classes from file: landsea_mask.dat
  !  ... fields as used in meteorological model

  implicit none

  integer, intent(in) :: callnum

  ios = 0

  if ( callnum == 1  ) then
  !.. Clay soil content    !
    ios = 0
    if(USES%DUST)then
      if(TEGEN_DATA)then
        !use global data interpolated to present grid

        if(MasterProc)write(6,*)'Sand and clay fractions from ',trim(Soil_TegenFile)

        call ReadField_CDF(Soil_TegenFile,'clay',clay_frac,1,  &
             interpol='conservative',needed=.true.,debug_flag=.false.)
        call ReadField_CDF(Soil_TegenFile,'sand',sand_frac,1,  &
             interpol='conservative',needed=.true.,debug_flag=.false.)

      elseif(MasterProc)then
        !use grid specific data
        call CheckStop('ASCII DUST NO MORE AVAILABLE! ')
      end if
    end if ! USES%DUST

  end if ! callnum == 1
end subroutine MetModel_LandUse

subroutine BLPhysics()
  !-----------------------------------------------------------------
  ! First written by Trond Iversen, modified by Hugo Jakobsen, 060994
  ! Extensive modifications to structure by Dave Simpson, March 2010.
  ! Some code moved to BLPhysics_mod, together with additinal options.
  ! Also now includes Haldis's use of NWP Kz values.
  !   ** not optimised, bujt called only at 3 h intervals
  !-----------------------------------------------------------------
  ! Calculate the exner function, geopotential height,
  ! and vertical exchange coefficient in sigma surfaces.
  ! The height zi of the "well mixed layer" or ABL-height is also calculated.
  !-----------------------------------------------------------------

  real, dimension(LIMAX,LJMAX,KMAX_MID)::exnm
  real, dimension(LIMAX,LJMAX,KMAX_BND)::exns

  real :: p_m, p_s, hs

  real, dimension(KMAX_BND) :: p_bnd !TESTzi
  real, dimension(KMAX_MID) :: Kz_nwp
  real    :: Kz_min, stab_h

  integer i,j,k,nr
  real :: theta2
  logical :: debug_flag
  logical,save :: first_call = .true.

  call CheckStop(KZ_SBL_LIMIT < 1.01*KZ_MINIMUM, "SBLlimit too low! in Met_mod")

  ! Preliminary definitions
  nr = 2
  if (first_call) nr = 1

  Kz_m2s(:,:,:)= 0.
  Kz_nwp(:)    = -99.0   ! store for printout. only set if read from NWP

  !..................................
  ! exner-functions (j/k kg)
  !..................................
  do k=1,KMAX_MID
    do j=1,ljmax
      do i=1,limax
        p_m = A_mid(k)+B_mid(k)*ps(i,j,nr)
        p_s = A_bnd(k)+B_bnd(k)*ps(i,j,nr)

        exnm(i,j,k)= CP * Exner_nd(p_m) ! exner (j/k kg)
        exns(i,j,k)= CP * Exner_nd(p_s)
      end do
    end do
  end do


  if ( debug_proc .and. DEBUG_Kz) then
    i = debug_iloc
    j = debug_jloc
    write(*,"(a,i4,2f12.5)") "TESTNR th ",nr,th(i,j,20,[1,nr])
    write(*,"(a,i4,2f12.5,es10.2)") "TESTNR fh ",nr,fh(i,j,[1,nr]),invL_nwp(i,j)
    write(*,"(a,i4,2es10.2)") "TESTNR ps ",nr,ps(i,j,[1,nr])
  end if

  !..................................
  ! Start choice of Kz and Hmix methods
  !..................................

  if (NWP_Kz .and. foundKz_met ) then  ! read from met data
    ! LAter we should remove Kz_met and Kz_m2s
    forall(i=1:limax,j=1:ljmax,k=2:KMAX_MID)
      SigmaKz(i,j,k,nr)=Kz_met(i,j,k,nr)/(60*60*METSTEP)
    end forall

    call SigmaKz_2_m2s( SigmaKz(:,:,:,nr), roa(:,:,:,nr),ps(:,:,nr), Kz_m2s )

    if( debug_proc ) Kz_nwp(:) = Kz_m2s(debug_iloc,debug_jloc,:) !for printout

    if( debug_proc .and. DEBUG_Kz)then
      write(6,*) '*** After Set SigmaKz', sum(SigmaKz(:,:,:,nr)), &
         minval(SigmaKz(:,:,:,nr)), maxval(SigmaKz(:,:,:,nr)), &
         DEBUG_Kz, 'NWP_Kz:',NWP_Kz, &
        '*** After convert to z',sum(Kz_m2s(:,:,:)), &
        minval(Kz_m2s(:,:,:)), maxval(Kz_m2s(:,:,:))
      write(6,*) 'After Set SigmaKz KTOP', Kz_met(debug_iloc,debug_jloc,1,nr)
    end if

  else   ! Not NWP Kz. Must calculate

    ! Get Kz first from PielkeBlackadar methods
    ! Use for all methods except NWP_Kz
    ! Do the physics for each i,j for now. Optimise later

    do j=1,ljmax
      do i=1,limax

        debug_flag = ( DEBUG_Kz .and. debug_proc .and. &
             i == debug_iloc .and. j == debug_jloc )

        call PielkeBlackadarKz ( &
          u_mid(i,j,:),  v_mid(i,j,:),  &
          z_mid(i,j,:),  z_bnd(i,j,:),  &
          th(i,j,:,nr),  Kz_m2s(i,j,:), &
          PIELKE, debug_flag )

      end do
    end do

    !======================================================================
    ! Hmix choices:
    if ( PBL%HmixMethod == "TIZi" ) then
      ! Get Mixing height from "orig" method
      ! "old" exner-function of the full-levels

      do j=1,ljmax
        do i=1,limax
          p_bnd(:) = A_bnd(:)+B_bnd(:)*ps(i,j,nr)

          call TI_Hmix ( &     ! Original EMEP method
             Kz_m2s(i,j,:), z_mid(i,j,:),  &
             z_bnd(i,j,:),  fh(i,j,nr),  &
             th(i,j,:,nr),  exnm(i,j,:),  &
             p_bnd(:), pzpbl(i,j), &
             .false.)

           pzpbl(i,j) = max( PBL%ZiMIN, pzpbl(i,j)) ! Keep old fixed height ZiMin here
           pzpbl(i,j) = min( PBL%ZiMAX, pzpbl(i,j))
        end do
      end do

    else ! Newer non-TI methods
      if ( PBL%HmixMethod == "SbRb" ) then
        call SeibertRiB_Hmix_3d(&
                u_mid(1:limax,1:ljmax,:),  &
                v_mid(1:limax,1:ljmax,:),  &
                z_mid(1:limax,1:ljmax,:),  &
                th(1:limax,1:ljmax,:,nr),  &
                pzpbl(1:limax,1:ljmax))

      elseif ( PBL%HmixMethod == "JcRb" ) then
        do i=1,limax
          do j=1,ljmax
            theta2 = t2_nwp(i,j,nr) * T_2_Tpot(ps(i,j,nr))
            call JericevicRiB_Hmix0(&
                u_mid(i,j,:), v_mid(i,j,:),  &
                z_mid(i,j,:), th(i,j,:,nr),  pzpbl(i,j))
          end do
        end do
      elseif ( PBL%HmixMethod == "NWP" ) then ! NWPHMIX
        do i=1,limax
          do j=1,ljmax
             pzpbl(i,j) = pbl_nwp(i,j,1)
          end do
        end do
        if ( DEBUG%MET .and. debug_proc) call datewrite("NWP HMIX: ", &
           [ pbl_nwp(debug_li,debug_lj,1),pbl_nwp(debug_li,debug_lj,2) ] )
      else
        call CheckStop("Need HmixMethod")
      end if ! end of newer methods

     ! Set limits on Zi
      forall(i=1:limax,j=1:ljmax)
        pzpbl(i,j) = max( PBL%ZiMIN, pzpbl(i,j))
        pzpbl(i,j) = min( PBL%ZiMAX, pzpbl(i,j) )
      end forall

    end if ! Hmix done
    !..spatial smoothing of new zi: Need fixed minimum here. 100 or 50 m is okay
    !  First, we make sure coastal areas had "land-like" values.

     if(LANDIFY_MET) &
         call landify(pzpbl,"pzbpl")
     call smoosp(pzpbl,PBL%ZiMIN,PBL%ZiMAX)

    !======================================================================
    ! Kz choices:

    if ( KzMethod == "JG" ) then  ! Jericevic/Grisogono for both Stable/Unstable
      do k = 2, KMAX_MID
        do j=1,ljmax
          do i=1,limax
            Kz_m2s(i,j,k) = JericevicKz(z_bnd(i,j,k),pzpbl(i,j),&
                                        ustar_nwp(i,j),Kz_m2s(i,j,k))
          end do
        end do
      end do

    else ! Specify unstable, stable separately:
      if ( StableKzMethod == "JG" ) then  ! Jericevic/Grisogono for both Stable/Unstable
        do j=1,ljmax
          do i=1,limax
           if ( invL_nwp(i,j) >= OB_invL_LIMIT ) then !neutral and unstable
             do k = 2, KMAX_MID
               if( z_bnd(i,j,k) <  pzpbl(i,j) ) then
                 Kz_m2s(i,j,k) = JericevicKz(z_bnd(i,j,k),pzpbl(i,j),&
                                             ustar_nwp(i,j),Kz_m2s(i,j,k))
               !else: keep Kz from Pielke/BLackadar
               end if
             end do
           end if
          end do
        end do
        if(debug_proc ) then
          i = debug_iloc
          j = debug_jloc
          if(DEBUG_Kz .and. invL_nwp(i,j) >= OB_invL_LIMIT ) then
            do k = 15, KMAX_MID
              print "(a,i3,f7.1,3es11.3)", "DEBUG SKz_m2s",k,&
                pzpbl(i,j), invL_nwp(i,j), ustar_nwp(i,j), Kz_m2s(i,j,k)
             end do
          end if
        end if

      elseif ( StableKzMethod == "BW" ) then
        do k = 2, KMAX_MID
          do j=1,ljmax
            do i=1,limax
              if ( invL_nwp(i,j) > 1.0e-10 ) then !stable ! leaves gap near zero
                Kz_m2s(i,j,k) = BrostWyngaardKz(z_bnd(i,j,k),pzpbl(i,j),&
                                ustar_nwp(i,j),invL_nwp(i,j),Kz_m2s(i,j,k))
              !else: keep Kz from Pielke/BLackadar
              end if
            end do
          end do
        end do

      else if ( StableKzMethod == "PB" ) then
        ! no change (keep Kz from Pielke/BLackadar)
      else
        call CheckStop("Need StableKzMethod")
      end if ! Stable Kz

      if ( UnstableKzMethod == "OB" ) then
        do j=1,ljmax
          do i=1,limax
            if ( invL_nwp(i,j) < OB_invL_LIMIT ) then !neutral and unstable
              call O_BrienKz ( &     ! Original EMEP method
                pzpbl(i,j),  z_bnd(i,j,:),  &
                ustar_nwp(i,j),  invL_nwp(i,j),  &
                Kz_m2s(i,j,:), .false.)
            end if
          end do
        end do
        if(debug_proc) then
          i = debug_iloc
          j = debug_jloc
          if(DEBUG_Kz .and. invL_nwp(i,j) <  OB_invL_LIMIT ) then
            do k = 15, KMAX_MID
              write(*,"(a,f7.1,3es10.3)") "DEBUG UKz_m2s", &
                pzpbl(i,j), invL_nwp(i,j), ustar_nwp(i,j), Kz_m2s(i,j,k)
            end do
          end if
        end if

      else
        call CheckStop("Need UnstableKzMethod")
      end if

    end if  ! Specify unstable, stable separately:
  end if ! NWP_Kz .and. foundKz_met

  ! spatial smoothing of new zi: Need fixed minimum here. 100 or 50 m is okay
  ! First, we make sure coastal areas had "land-like" values.


  !************************************************************************!
  ! test some alternative options for Kz and Hmix
  if( DEBUG_BLM .and. debug_proc .and. modulo( current_date%hour, 3)  == 0 &
         .and. current_date%seconds == 0  ) then

    i = debug_iloc
    j = debug_jloc
    p_bnd(:) = A_bnd(:)+B_bnd(:)*ps(i,j,nr)

  !************************************************************************!
  ! We test all the various options here. Pass in  data as keyword arguments
  ! to avoid possible errors!

    call Test_BLM( mm=current_date%month, dd=current_date%day, &
           hh=current_date%hour, fH=fh(i,j,nr), &
           u=u_mid(i,j,:),v=v_mid(i,j,:), zm=z_mid(i,j,:), &
           zb=z_bnd(i,j,:), exnm=exnm(i,j,:), Kz=Kz_m2s(i,j,:), &
           Kz_nwp=Kz_nwp(:), invL=invL_nwp(i,j), &
           q=q(i,j,:,nr),  & ! TEST Vogel
           ustar=ustar_nwp(i,j), th=th(i,j,:,nr), pb=p_bnd(:), zi=pzpbl(i,j))
  !************************************************************************!
    hs = z_bnd(i,j,KMAX_MID)

    stab_h = min( PsiH(hs*invL_nwp(i,j)), 0.9 )
    Kz_min = ustar_nwp(i,j)*KARMAN*hs /( 1 - stab_h  )
    write(*,"(a,10f10.3)") "PSIH ", stab_h, fh(i,j,nr), invL_nwp(i,j), &
             PsiH(hs*invL_nwp(i,j)),Kz_min

  end if ! end of debug extra options


  !***************************************************
  if ( .not. (NWP_Kz .and. foundKz_met) ) then  ! convert to Sigma units

    call Kz_m2s_toSigmaKz (Kz_m2s(1:limax,1:ljmax,:),roa(1:limax,1:ljmax,:,nr),&
         ps(1:limax,1:ljmax,nr),SigmaKz(1:limax,1:ljmax,:,nr))
    call Kz_m2s_toEtaKz (Kz_m2s(1:limax,1:ljmax,:),roa(1:limax,1:ljmax,:,nr),&
         ps(1:limax,1:ljmax,nr),EtaKz(1:limax,1:ljmax,:,nr),Eta_mid,A_mid,B_mid)

  end if
  !***************************************************

    first_call=.false.
end subroutine BLPhysics

subroutine smoosp(f,rmin,rmax)
  !-----------------------------------------------------------------
  ! Written by Trond Iversen,  modified by Hugo Jakobsen, 080994
  ! parallellized and modified by Peter February 2003
  !----------------------------------------------------------------------
  ! Apply the shapiro filter with s=0.5 and s=-0.5
  ! to the field f using h as a work space also the boundaries
  ! are smoothed. f contains the smoothed field upon return.
  !----------------------------------------------------------------------
  implicit none

  real, intent(inout) :: f(LIMAX,LJMAX) ! data to be smoothed
  real, intent(in)    :: rmin,rmax      ! min/max allowed

  real, dimension(LIMAX+4,LJMAX+4) :: h1,h2 ! help variables
  real, dimension(LIMAX,2)         :: f_south,f_north
  real, dimension(LJMAX+2*2,2)     :: f_west,f_east
  real s

  integer :: thick
  integer :: iif,jjf,is,i,j,ii,jj,iifl,jjfl

  iif=limax
  jjf=ljmax

  thick=2  !we fetch 2 neighbors at once, so that we don't need to call
  ! readneighbours twice
  iifl=iif+2*thick
  jjfl=jjf+2*thick

  call readneighbors(f,f_south,f_north,f_west,f_east,thick)

  do j=1,jjf
    jj=j+thick
    do i=1,iif
      ii=i+thick
      h1(ii,jj) = f(i,j)
    end do
  end do
  do j=1,thick
    do i=1,iif
      ii=i+thick
      h1(ii,j) = f_south(i,j)
    end do
  end do

  do j=1,thick
    jj=j+jjf+thick
    do i=1,iif
      ii=i+thick
      h1(ii,jj) = f_north(i,j)
    end do
  end do

  do j=1,jjfl
    do i=1,thick
      h1(i,j) = f_west(j,i)
    end do
  end do

  do j=1,jjfl
    do i=1,thick
      ii=i+iif+thick
      h1(ii,j) = f_east(j,i)
    end do
  end do

  do j=1,jjfl
    h2(1,j) = 0.
    h2(iifl,j) = 0.
  end do

  do i=1,iifl
    h2(i,1) = 0.
    h2(i,jjfl) = 0.
  end do

  do is=2,1,-1
    s=is-1.5  !s=0,5 s=-0.5
    if(is /= 2)h1=h2

    ! the smoothing
    do  j=2,jjfl-1
      do  i=2,iifl-1
        h2(i,j)=(1.-2.*s+s*s)*h1(i,j)&
            + 0.5*s*(1.-s)*(h1(i+1,j)+h1(i-1,j)+h1(i,j+1)+h1(i,j-1))  &
            + s*s*(h1(i+1,j+1)+h1(i-1,j-1)+h1(i+1,j-1)+h1(i-1,j+1))/4.
        h2(i,j) = amax1(h2(i,j),rmin)
        h2(i,j) = amin1(h2(i,j),rmax)
      end do
    end do
  end do


  do j=1,jjf
    jj=j+thick
    do i=1,iif
      ii=i+thick
      f(i,j)=h2(ii,jj)
    end do
  end do

end subroutine smoosp

subroutine extendarea(f,h,debug_flag)
  !----------------------------------------------------------------------
  ! based upon the smoosp routine
  ! returns extended array array, reading neighbour procs as needed
  !----------------------------------------------------------------------

  real, intent(in) :: f(LIMAX,LJMAX)
  real, intent(inout) :: h(:,:)
  logical, intent(in), optional :: debug_flag
  logical :: mydebug = .false.

  real, dimension(size(f,1),2)     :: f_south,f_north
  real, dimension(size(f,2)+2*2,2) :: f_west,f_east

  integer :: thick ! = size(h,1) - size(f,1) ! Caller has to make h > f
  integer :: iif,jjf,i,j,ii,jj,iifl,jjfl
  if ( present(debug_flag)  ) mydebug = debug_flag

  thick = ( size(h,1) - size(f,1) ) ! Caller has to make h > f ;NB: NOT SAFE!
  iif=limax
  jjf=ljmax

  if( modulo(thick,2) /= 0 ) then
     write(*,*) "ERROR extendarea para,s ", me, iif , jjf, thick
     write(*,*) "ERROR extendarea mod ", modulo(thick,2)
     call CheckStop("ERROR extendarea thickness not even!")
  end if
  thick = thick / 2


  ! readneighbours twice
  iifl=iif+2*thick
  jjfl=jjf+2*thick
  if(mydebug .and. MasterProc ) write(*,*) "DEBUG extendarea", iif,jjf,thick

  call readneighbors(f,f_south,f_north,f_west,f_east,thick)

  do j=1,jjf
    jj=j+thick
    do i=1,iif
      ii=i+thick
      h(ii,jj) = f(i,j)
    end do
  end do
  do j=1,thick
    do i=1,iif
      ii=i+thick
      h(ii,j) = f_south(i,j)
    end do
  end do

  do j=1,thick
    jj=j+jjf+thick
    do i=1,iif
      ii=i+thick
      h(ii,jj) = f_north(i,j)
    end do
  end do

  do j=1,jjfl
    do i=1,thick
      h(i,j) = f_west(j,i)
    end do
  end do

  do j=1,jjfl
    do i=1,thick
      ii=i+iif+thick
      h(ii,j) = f_east(j,i)
    end do
  end do

end subroutine extendarea

subroutine landify(x,intxt,xmin,xmax,wfmin,xmask)
  real, dimension(LIMAX,LJMAX), intent(inout) :: x
  character(len=*), intent(in), optional :: intxt
  real, intent(in), optional :: xmin, xmax  ! Limits of valid data for x
  real, intent(in), optional :: wfmin ! Limits of valid data for water frac
  logical, dimension(LIMAX,LJMAX), intent(in), optional :: xmask

  logical, dimension(LIMAX,LJMAX) :: mask
  real, dimension(LIMAX+2*NEXTEND,LJMAX+2*NEXTEND)  :: xx  ! extended
  character(len=30) :: txt, masktxt
  real :: xwfmin, xxmin, xxmax
  logical :: debug_flag=.false.
  real :: sumland, sumx, landfrac, oldx
  integer :: i,j, ii, jj, ii2, jj2

  txt = "Landify: "
  if ( present(intxt) )  txt = trim(txt) // trim(intxt)
  xwfmin = 0.5 ! Default fraction of water
  if ( present(wfmin) )  xwfmin = wfmin
  xxmin = 1.0e-10  ! Default  min value x
  if ( present(xmin) )  xxmin = xmin
  xxmax = 1.0e30   ! Default max value x
  if ( present(xmax) )  xxmax = xmax

  if(DEBUG_LANDIFY.and.MasterProc) then
      write(*,*) trim(txt) , water_frac_set
      write(*,"(a,2g12.4)") 'Data Limits ', xxmin, xxmax
      write(*,"(a,g12.4)")  'Water Limit ', xwfmin
  end if

  if( .not. water_frac_set  ) then
     if(MasterProc) write(*,*) trim(txt) //  " skips 1st NTERM"
     write(*,*) trim(txt) //  " skips 1st NTERM"
     return   !  on 1st time-step water_frac hasnt yet been set.
  end if

  if ( present(xmask) )  then
    mask = xmask
    masktxt = "Input mask"
  else
    mask = likely_coastal
    masktxt = "Coastal mask"
  end if

   if ( DEBUG_LANDIFY.and. debug_proc ) then
      write(*,"(a,6i4,L2,1x,a)") "DLandify start ", &
       debug_li, debug_lj, 1, limax, 1, ljmax,  xwf_done, trim(masktxt)
   end if

  ! We need the extended water-fraction too, but just once
  if ( .not. xwf_done ) then ! only need to do this once
       if ( DEBUG%MET .and. debug_proc) write(*,*) "Landify xwf"
       call extendarea( water_fraction(:,:), xwf, debug_flag)
       xwf_done = .true.
  end if

 ! Then the data we are working with:

  call extendarea( x(:,:), xx(:,:), debug_flag )

  if ( DEBUG_LANDIFY .and. debug_proc) write(*,*) "Landify now ", &
     xwf_done , likely_coastal(debug_li,debug_lj), mask(debug_li,debug_lj)

  oldx = 0.0
  if( debug_proc ) oldx = x(debug_li, debug_lj)

  do j = 1, ljmax
    do i = 1, limax
      ! Take a 5x5 average of the land-weighted values for SW. Seems
      ! best not to "believe" NWP models too much for this param, and the
      ! variation in a grid is so big anyway. We aim at the broad effect.

      sumland  = 0.0
      sumx     = 0.0
      debug_flag = ( DEBUG_LANDIFY .and. debug_proc .and. &
             i==debug_li .and. j==debug_lj )

      if( mask(i,j) ) then ! likely coastal or water_frac <0.0 for SW
        do jj = -NEXTEND, NEXTEND
          do ii = -NEXTEND, NEXTEND
            ii2=i+ii+NEXTEND  ! coord in extended array !CHECK!
            jj2=j+jj+NEXTEND

            !Had 0.5, 1.0e-10 in original for met-data
            if( xwf(ii2,jj2)<= wfmin .and. &! was 0.5 likely not NWP sea
                xx(ii2,jj2) <= xxmax .and. &! Valid x range
                xx(ii2,jj2) >= xxmin ) then !
              landfrac    =  1.0 - xwf(ii2,jj2)
              sumland = sumland + landfrac
              sumx    = sumx    + landfrac * xx(ii2,jj2)
              if ( debug_flag ) then
                write(*,"(a,4i4,2f7.4,3g11.3,2f7.4)") "DBG"//trim(intxt), i,j, &
                   ii2, jj2, water_fraction(i,j), xwf(ii2,jj2), x(i,j), &
                   xx(ii2,jj2), sumx, landfrac, sumland
              end if ! DEBUG
            end if ! xsw
          end do ! ii
        end do ! jj

        if ( sumland > 0.001 ) then ! replace x with land-weighted values
          x(i,j) = sumx/sumland
          if ( debug_flag ) then
             write(*,"(a,2i4,8g12.3)") "DBGDONE", i,j, &
                 water_fraction(i,j), sumx, sumland, x(i,j)
          end if
        end if ! water_fraction

       end if ! likely_coastal
     end do ! i
  end do ! j

  if ( DEBUG_LANDIFY .and. debug_proc ) then
    call datewrite("LandifyDONE: "//trim(intxt), (/ oldx, x(debug_li,debug_lj) /) )
  end if

end subroutine landify

subroutine readneighbors(data,data_south,data_north,data_west,data_east,thick)
  !----------------------------------------------------------------------
  ! Read data at the other side of the boundaries
  !
  ! thick is the number of gridcells in each direction to be transferred
  ! Note that we also fetch data from processors in the "diagonal"
  ! directions
  !
  ! Written by Peter February 2003
  !
  ! Note,
  ! The data_west(jj,:)=data(1,j) is not a bug: when there is no west
  ! neighbour,
  ! the data is simply copied from the nearest points: data_west(jj,:) should
  ! be =data(-thick+1:0,j), but since this data does not exist, we
  ! put it =data(1,j).
  !----------------------------------------------------------------------


  implicit none

  integer, intent(in) :: thick
  real,intent(in), dimension(LIMAX,LJMAX) ::data
  real,intent(out), dimension(LIMAX,thick) ::data_south,data_north
  real,intent(out), dimension(LJMAX+2*thick,thick) ::data_west,data_east
  real, dimension(LIMAX,thick) ::data_south_snd,data_north_snd
  real, dimension(LJMAX+2*thick,thick) ::data_west_snd,data_east_snd

  integer :: msgnr
  integer :: j,tj,jj,jt

  !check that limax and ljmax are large enough
  call CheckStop(limax < thick, "ERROR readneighbors in Met_mod")
  call CheckStop(ljmax < thick, "ERROR readneighbors in Met_mod")


  msgnr=1

  data_south_snd(:,:)=data(:,1:thick)
  data_north_snd(:,:)=data(:,ljmax-thick+1:ljmax)
  if(neighbor(SOUTH) >= 0 )then
    CALL MPI_ISEND( data_south_snd , 8*LIMAX*thick, MPI_BYTE,&
          neighbor(SOUTH), msgnr, MPI_COMM_CALC, request_s,IERROR)
  end if
  if(neighbor(NORTH) >= 0 )then
    CALL MPI_ISEND( data_north_snd , 8*LIMAX*thick, MPI_BYTE,&
          neighbor(NORTH), msgnr+9, MPI_COMM_CALC, request_n,IERROR)
  end if

  if(neighbor(SOUTH) >= 0 )then
    CALL MPI_RECV( data_south, 8*LIMAX*thick, MPI_BYTE,&
          neighbor(SOUTH), msgnr+9, MPI_COMM_CALC, MPISTATUS, IERROR)
  else
    do tj=1,thick
      data_south(:,tj)=data(:,1)
    end do
  end if
  if(neighbor(NORTH) >= 0 )then
    CALL MPI_RECV( data_north, 8*LIMAX*thick, MPI_BYTE,&
          neighbor(NORTH), msgnr, MPI_COMM_CALC, MPISTATUS, IERROR)
  else
    do tj=1,thick
      data_north(:,tj)=data(:,ljmax)
    end do
  end if

  jj=0
  do jt=1,thick
    jj=jj+1
    data_west_snd(jj,:)=data_south(1:thick,jt)
    data_east_snd(jj,:)=data_south(limax-thick+1:limax,jt)
  end do
  do j=1,ljmax
    jj=jj+1
    data_west_snd(jj,:)=data(1:thick,j)
    data_east_snd(jj,:)=data(limax-thick+1:limax,j)
  end do
  do jt=1,thick
    jj=jj+1
    data_west_snd(jj,:)=data_north(1:thick,jt)
    data_east_snd(jj,:)=data_north(limax-thick+1:limax,jt)
  end do

  if(neighbor(WEST) >= 0 )then
    CALL MPI_ISEND( data_west_snd , 8*(LJMAX+2*thick)*thick, MPI_BYTE,&
          neighbor(WEST), msgnr+3, MPI_COMM_CALC, request_w,IERROR)
  end if
  if(neighbor(EAST) >= 0 )then
    CALL MPI_ISEND( data_east_snd , 8*(LJMAX+2*thick)*thick, MPI_BYTE,&
          neighbor(EAST), msgnr+7, MPI_COMM_CALC, request_e,IERROR)
  end if

  if(neighbor(WEST) >= 0 )then
    CALL MPI_RECV( data_west, 8*(LJMAX+2*thick)*thick, MPI_BYTE,&
          neighbor(WEST), msgnr+7, MPI_COMM_CALC, MPISTATUS, IERROR)
  else
    jj=0
    do jt=1,thick
      jj=jj+1
      data_west(jj,:)=data_south(1,jt)
    end do
    do j=1,ljmax
      jj=jj+1
      data_west(jj,:)=data(1,j)
    end do
    do jt=1,thick
      jj=jj+1
      data_west(jj,:)=data_north(1,jt)
    end do
  end if
  if(neighbor(EAST) >= 0 )then
    CALL MPI_RECV( data_east, 8*(LJMAX+2*thick)*thick, MPI_BYTE, &
          neighbor(EAST), msgnr+3, MPI_COMM_CALC, MPISTATUS, IERROR)
  else
    jj=0
    do jt=1,thick
      jj=jj+1
      data_east(jj,:)=data_south(limax,jt)
    end do
    do j=1,ljmax
      jj=jj+1
      data_east(jj,:)=data(limax,j)
    end do
    do jt=1,thick
      jj=jj+1
      data_east(jj,:)=data_north(limax,jt)
    end do
  end if

  if(neighbor(SOUTH) >= 0 )then
     CALL MPI_WAIT(request_s, MPISTATUS,IERROR)
  end if
  if(neighbor(NORTH) >= 0 )then
     CALL MPI_WAIT(request_n, MPISTATUS,IERROR)
  end if
  if(neighbor(WEST) >= 0 )then
     CALL MPI_WAIT(request_w, MPISTATUS, IERROR)
  end if
  if(neighbor(EAST) >= 0 )then
     CALL MPI_WAIT(request_e, MPISTATUS,IERROR)
  end if

end subroutine readneighbors

subroutine tkediff (nr)
  !----------------------------------------------------------------------
  ! Compute vertical eddy diffusivities as a function
  ! altitude, height of PBL, and a velocity scale, square root of
  ! turbulent kinetic energy (TKE). This is a non-local scheme.
  ! The TKE at the surface is diagnosed using scales for horizontaland
  ! vertical velocities (ustar and wstar) in the surface layer
  ! (Alapaty 2004; Holstag et al. 1990 and Mihailovic et al. 2004)
  ! PBL ht is calculated using the EMEP formulation
  !
  ! Written by DT Mihailovic (October 2004)
  ! EMEP polishing and comments: JE Jonson and P Wind
  !----------------------------------------------------------------------

  implicit none

  !     Local constants
  real, parameter :: &
    SZKM=1600.     &   ! Constant (Blackadar, 1976)
   ,CKZ=0.001      &   ! Constant (Zhang and Athens, 1982)
!  ,REFPR=1.0E+05  &   ! Referent pressure
   ,KZ0LT=1.0E-04  &   ! Constant (Alapaty et al., 1997)
   ,RIC=0.10      !&   ! Critical Richardson number (Holstlag et al., 1993)
!  ,ROVG=RGAS_KG/GRAV  ! Used in Calculation of R-number

  ! INPUT
  integer, intent(in) :: nr  ! Number of meteorological stored in arrays (1 or 2)

  ! OUTPUT
  !  skh(i,j,k,nr) array
  !  Values of the Kz coefficients (eddyz (i,j,k)) are transformed nto
  !  sigma system and then they stored in this array which is later used
  !  in ADVECTION module

  ! Local arrays
  integer, dimension(LIMAX,LJMAX)      :: iblht   ! Level of the PBL top
  real, dimension(LIMAX,LJMAX,KMAX_BND):: eddyz   ! Eddy coefficients (m2/s)
  real, dimension(LIMAX,LJMAX,KMAX_MID):: &
    t_virt,   & ! Potential temperature (K)
    e,        & ! Kinetic energy with respect to height (m2/s2)
    dzq,      & ! Thickness of sigma interface layers (m)
    u_mid,    & ! Wind speed in x-direction (m/s)
    v_mid       ! Wind speed in y-direction (m/s)

  real, dimension(LIMAX,LJMAX,KMAX_MID-1):: &
    dza         ! Thickness of half sigma layers (m)
  real, dimension(LIMAX,LJMAX):: &
    pblht ,   & ! PBL (Holstag, 1990) (m)
    h_flux,   & ! Sensible heat flux  (W/m2)
    ust_r ,   & ! Friction velocity (m/s)
    mol   ,   & ! Monin-obukhov length (m)
    wstar       ! Convective velocity (m/s)

  real, dimension(KMAX_BND) :: &
    rib         ! Bulk Richardson number

  real, dimension(KMAX_MID) :: &
    rich,      & ! Richardson number
    psi_zi       ! Used in the vertical integration

  real, dimension (10) ::      &
    psi_z,      & ! Used for calculating
    zovh          ! TKE

  !     Local variables
  real :: dtmp, tog, wssq1, wssq2, wssq, tconv, wss, wst, PSI_TKE, &
          dusq, dvsq, ri, ss, dthdz, busfc, zvh,                   &
          part1, part2, fract1, fract2, apbl, kz0, u_s, goth

  integer :: i, j, k, l, kcbl

  call CheckStop('This routine is not ready. for example ust_r and kcbl are not set!')

  ! Functions for averaging the vertical turbulent kinetic energy
  !   (Alapaty, 2003)
  data psi_z /0.00,2.00,1.85,1.51,1.48,1.52,1.43,1.10,1.20,0.25/
  data zovh  /0.00,0.05,0.10,0.20,0.40,0.60,0.80,1.00,1.10,1.20/

  ! Store the NMW meteorology and variables derived from its

  ! Change the sign
  h_flux(1:limax,1:ljmax)=-fh(1:limax,1:ljmax,nr)

  ! Avoid devision by zero later in the code
  where (ABS(h_flux(1:limax,1:ljmax))<0.0001) h_flux(1:limax,1:ljmax)=0.0001

  ! Check PBL height
  ! strange tests! Negative pzpbl check? From 1 to 100m, odd!
  do i=1,limax
    do j=1,ljmax
      if(ABS(pzpbl(i,j)) < 1.) then
        pzpbl(i,j)=100.
      end if
    end do
  end do

  ! Calculate velocity components in the (h) poits (Arakawa notation)
  do k=1,KMAX_MID
    do i=1,limax
      do j=1,ljmax
      ! u_mid(i,j,k)=0.5*(u_xmj(i-1,j  ,k,nr)+u_xmj(i,j,k,nr))
      ! v_mid(i,j,k)=0.5*(v_xmi(i  ,j-1,k,nr)+v_xmi(i,j,k,nr))
        u_mid(i,j,k)=u_xmj(i,j,k,nr)
        v_mid(i,j,k)=v_xmi(i,j,k,nr)
      end do
    end do
  end do

  ! Avoid small values
  where (ABS(u_mid(1:limax,1:ljmax,1:KMAX_MID))<0.001) &
    u_mid(1:limax,1:ljmax,1:KMAX_MID)=0.001
  where (ABS(v_mid(1:limax,1:ljmax,1:KMAX_MID))<0.001) &
    v_mid(1:limax,1:ljmax,1:KMAX_MID)=0.001

  ! Initialize eddy difussivity arrays
  eddyz(1:limax,1:ljmax,1:KMAX_MID)=0.

  ! Calculate tickness of the full layers
  dzq(1:limax,1:ljmax,1:KMAX_MID) = z_bnd(1:limax,1:ljmax,1:KMAX_MID)  &
       - z_bnd(1:limax,1:ljmax,2:KMAX_BND)

  ! ... and the half sigma layers
  dza(1:limax,1:ljmax,1:KMAX_MID-1) = z_mid(1:limax,1:ljmax,1:KMAX_MID-1)          &
       - z_mid(1:limax,1:ljmax,2:KMAX_MID)

  ! Calculate virtual temperature
  t_virt(1:limax,1:ljmax,1:KMAX_MID) = th(1:limax,1:ljmax,1:KMAX_MID,nr)  &
       * (1.0+0.622*q(1:limax,1:ljmax,1:KMAX_MID,nr))

  ! Calculate Monin-Obuhkov length   (Garratt, 1994)
  do i=1,limax
    do j=1,ljmax
      u_s = ustar_nwp(i,j)
      mol(i,j) = -(ps(i,j,nr)*u_s*u_s*u_s)/(KARMAN*GRAV*h_flux(i,j)*KAPPA)
    end do
  end do

  ! Calculate the convective velocity (wstar)
  do i=1,limax
    do j=1,ljmax
      wstar(i,j) = GRAV*h_flux(i,j)*pzpbl(i,j)/rho_surf(i,j)    &
           /CP/th(i,j,KMAX_MID,nr)
      if(wstar(i,j) < 0.) then
        wstar(i,j)=-ABS(wstar(i,j))**(0.3333)
      else
        wstar(i,j)=(wstar(i,j))**(0.3333)
      end if
    end do
  end do

  ! Start with a long loop
  DO i=1,limax
    DO j=1,ljmax

      rib(1:KMAX_MID) = 0.0        ! Initialize bulk Richardson number

      part1=ust_r(i,j)*ust_r(i,j)*ust_r(i,j)
      wst=AMAX1(wstar(i,j),1.0E-20)
      part2=0.6*wst*wst*wst
      wss=AMAX1(1.0E-4,(part1+part2))
      wss=EXP(0.333333*ALOG(wss))

      if (h_flux(i,j) < 0.0) then
        tconv=0.0                                ! Holstlag et al. (1990)
      else
        tconv=8.5*h_flux(i,j)/rho_surf(i,j)/CP/wss ! Conversion to kinematic flux
      end if

      do k=KMAX_MID,1,-1
        dtmp=t_virt(i,j,k)-t_virt(i,j,KMAX_MID)-tconv
        tog=0.5*(t_virt(i,j,k)+t_virt(i,j,KMAX_MID))/GRAV
        wssq1=u_mid(i,j,k)*u_mid(i,j,k)
        wssq2=v_mid(i,j,k)*v_mid(i,j,k)
        wssq=wssq1+wssq2
        wssq=AMAX1(wssq,1.0E-4)
        rib(k)=z_mid(i,j,k)*dtmp/(tog*wssq)
        if(rib(k).ge.RIC) go to 9001
      end do
9001  continue

      ! Calculate PBL height according to Holtslag et al. (1993)
      pblht(i,j)=0.
      if(k.ne.KMAX_MID) then
        fract1=(RIC-rib(k+1))/(rib(k)-rib(k+1))
        fract2=1.-fract1
        apbl=z_mid(i,j,k)*fract1
        pblht(i,j)=apbl+z_mid(i,j,k+1)*fract2
        if(pblht(i,j) > z_bnd(i,j,k+1)) then
          kcbl=k
        else
          kcbl=k+1
        end if
      end if
      iblht(i,j)=kcbl

      if(pblht(i,j)<z_bnd(i,j,KMAX_MID)) then
        pblht(i,j)=z_bnd(i,j,KMAX_MID)
        iblht(i,j)=KMAX_MID
      end if

      if(pblht(i,j).le.100.) then              !Minimum of PBL height
        pblht(i,j)=100.
      end if

      ! Find the critical Richardson number (Shir and Borestein, 1976)
      do k=2,iblht(i,j)-1
        rich(k)=0.257*dza(i,j,k)**0.175
      end do

      ! Free troposphere and cloudy case Kz values estimation
      do k = 2,iblht(i,j)-1
        dusq = (u_mid(i,j,k-1)-u_mid(i,j,k))*(u_mid(i,j,k-1)-u_mid(i,j,k))
        dvsq = (v_mid(i,j,k-1)-v_mid(i,j,k))*(v_mid(i,j,k-1)-v_mid(i,j,k))
        ss = (dusq+dvsq)/(dza(i,j,k-1)*dza(i,j,k-1))+1.E-9
        goth = 2.*GRAV/(t_virt(i,j,k-1)+t_virt(i,j,k))
        dthdz = (t_virt(i,j,k-1)-t_virt(i,j,k))/dza(i,j,k-1)
        ri = goth*dthdz/ss

        ! (Duran and Clemp, 1982)
        kz0 = CKZ*dzq(i,j,k)
        if (ri-rich(k) > 0.) then
          eddyz(i,j,k)=kz0
        else
          eddyz(i,j,k)=kz0+SZKM*SQRT(ss)*(rich(k)-ri)/rich(k)
        end if
        eddyz(i,j,k)=AMIN1(eddyz(i,j,k),100.)
      end do

      ! Eddy diffusivity coefficients for all regimes in the mixed layer
      do  k=iblht(i,j),KMAX_MID
        if (mol(i,j) < 0.0) then                 !Unstable conditions
           ri=(1.0-15.*z_mid(i,j,k)/mol(i,j))**(-0.25)
           ri=ri/KARMAN/z_mid(i,j,k)
           ri=ri*AMAX1(0.0,pblht(i,j)-z_mid(i,j,k))
           dthdz=ri*ust_r(i,j)**3.
           goth=AMAX1(wstar(i,j),0.0)
           dusq=0.4*goth**3.
           ri=(dthdz+dusq)**(2./3.)
           e(i,j,k)=0.5*ri*(2.6)**(2./3.)        !Moeng and Sullivan (1994)
        else
           ri=z_bnd(i,j,k)/pblht(i,j)               !Stable
           ri=z_mid(i,j,k)/pblht(i,j)               !New
           ri=(1.0-ri)
           ri=AMAX1(0.0,ri)
           ri=(1.0-ri)**1.75
           e(i,j,k)=6.*ust_r(i,j)*ust_r(i,j)*ri  !Lenshow(1988)
        end if

        ! Calculate Ksi function using interpolation in the vertical
        ! Alapaty (2001, 2003)
        zvh=z_mid(i,j,k)/pblht(i,j)
        do l=1,9
          if (zvh > zovh(l).and. zvh < zovh(l+1)) then
            psi_zi(k)=(psi_z(l+1)-psi_z(l))/(zovh(l+1)-zovh(l))
            psi_zi(k)=psi_zi(k)*(zvh-zovh(l))
            psi_zi(k)=psi_zi(k)+psi_z(l)
            psi_zi(k)=psi_zi(k)/2.0               !Normalized the value
          end if
        end do
      end do

      ! Calculate integral for Ksi
      psi_tke=0.
      do k=KMAX_MID,iblht(i,j),-1
        psi_tke=psi_tke+psi_zi(k)*dzq(i,j,k)*sqrt(e(i,j,k))
      end do
      psi_tke=psi_tke/pblht(i,j)

      do k=iblht(i,j),KMAX_MID          !Calculate coefficients
        goth=psi_tke
        goth=goth*KARMAN*z_mid(i,j,k)
        dthdz=z_mid(i,j,k)/pblht(i,j)
        dthdz=1.0-dthdz
        dthdz=AMAX1(1.0E-2,dthdz)
        if(mol(i,j) > 0.0) then                        ! Stable
          goth=sqrt(e(i,j,iblht(i,j)))                 ! Mihailovic (2004)
          goth=goth*KARMAN*z_mid(i,j,k)
          dthdz=z_mid(i,j,k)/pzpbl(i,j)
          dthdz=1.0-dthdz
          dthdz=AMAX1(1.0E-2,dthdz)
          busfc=0.74+4.7*z_mid(i,j,KMAX_MID)/mol(i,j)
          busfc=AMAX1(busfc,1.0)
          dthdz=dthdz**1.50                            ! test (2004)
          eddyz(i,j,k)=goth*dthdz/busfc
        else
          dthdz=dthdz*dthdz
          busfc=1.0
          eddyz(i,j,k)=goth*dthdz/busfc
        end if
      end do

      ! Checking procedure
      do k=2,iblht(i,j)-1
        if(eddyz(i,j,k).le.0.0) THEN
          eddyz(i,j,k)= KZ0LT
        end if
      end do

      ! Avoid phisically unrealistic values
      do k=2,KMAX_MID
        IF(eddyz(i,j,k).le.0.1) then
          eddyz(i,j,k)=0.1
        end if
      end do

      ! To avoid loss of mass/energy through top of the model
      ! put eddyz (I,J,K) to zero at the last  level from top
      eddyz(i,j,KMAX_BND)=0.0

      ! Calculate eddy coefficients at the interfaces
      do k=2,KMAX_MID
        eddyz(i,j,k)=0.5*(eddyz(i,j,k-1)+eddyz(i,j,k))
      ! if(i.eq.10.and.j.eq.10.) then
      !   if (abs(u_xmj(i,j ,k,nr)-u_mid(i,j,k)).gt.5.) then
      !      print *,"NEW ",i,j,u_xmj(i,j ,KMAX_MID,nr),u_mid(i,j,KMAX_MID)
      ! end if
      end do

      ! Transform values of the eddy coeficients into the the sigma coordinate
      do k=2,KMAX_MID
        eddyz(i,j,k)=eddyz(i,j,k)*((sigma_mid(k)-sigma_mid(    k-1))/   &
              (    z_mid(i,j,k)-z_mid(i,j,k-1)))**2.
      end do
     ENDDO
  ENDDO ! long loop

  !     Store diffusivity coefficients into skh(i,j,k,nr) array
  do k=2,KMAX_MID
    do i=1,limax
      do j=1,ljmax
         SigmaKz(i,j,k,nr)=eddyz(i,j,k)
      end do
    end do
  end do

  ! For plotting set pblht  =  pzpbl
  pzpbl(:,:) = pblht(:,:)
end subroutine tkediff

subroutine Getmeteofield(meteoname,namefield,nrec,&
                         ndim,unit,validity,field,needed,found)
  !----------------------------------------------------------------------
  ! Read the meteofields and distribute to nodes
  !----------------------------------------------------------------------
  implicit none

  real, dimension(*),intent(out)  :: field ! dimensions: (LIMAX,LJMAX)

  character(len=*),intent(in)  :: meteoname,namefield
  character(len=*),intent(out) :: unit,validity
  integer,intent(in)           :: nrec,ndim
  logical,intent(in) ,optional :: needed
  logical,intent(out),optional :: found
  character(len=len(namefield)) :: namefield_met

  real :: scalefactors(2)
  integer :: KMAX,ijk,i,k,j,nfetch,k1,k2,istart,jstart,Nlevel,kstart,kend
  logical :: reverse_k
  real, allocatable,save ::meteo_3D(:,:,:)

  validity=''
  call_msg = "GetMeteofield" // trim(namefield)

  if(ndim==3)KMAX=KMAX_MET
  if(ndim==2)KMAX=1

  if(MANUAL_GRID)then
    Nlevel=37
!   call ReadField_CDF(meteoname,namefield,field, &
!         Nlevel,interpol='conservative',needed=.true.,debug_flag=.false.)
    if(.not.allocated(meteo_3D))allocate(meteo_3D(Nlevel,LIMAX,LJMAX))
      meteo_3D=0.0
      if(ndim==3)then
        if(trim(namefield)=='u_wind')then
!         call ReadField_CDF(meteoname,namefield,meteo_3D,nstart=nrec,kstart=1,kend=Nlevel,interpol='zero_order', &
          call ReadField_CDF(meteoname,namefield,meteo_3D,nstart=nrec,kstart=1,kend=Nlevel,interpol='conservative', &
!               use_lat_name='lat_u', use_lon_name='lon_u', &
                stagg='stagg_u',&
                needed=needed,found=found,unit=unit,debug_flag=.false.)
        elseif(trim(namefield)=='v_wind') then
!         call ReadField_CDF(meteoname,namefield,meteo_3D,nstart=nrec,kstart=1,kend=Nlevel,interpol='zero_order', &
          call ReadField_CDF(meteoname,namefield,meteo_3D,nstart=nrec,kstart=1,kend=Nlevel,interpol='conservative', &
!               use_lat_name='lat_v', use_lon_name='lon_v', &
                stagg='stagg_v',&
                needed=needed,found=found,unit=unit,debug_flag=.false.)
        else
          call ReadField_CDF(meteoname,namefield,meteo_3D,nstart=nrec,kstart=1,kend=Nlevel,interpol='zero_order', &
                  needed=needed,found=found,unit=unit,debug_flag=.false.)
        end if
        validity='not set'
        ! CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)

        ! interpolate vertically
        call vertical_interpolate(meteoname,meteo_3D,Nlevel,field,.false.)

      elseif(ndim==2)then
        call ReadField_CDF(meteoname,namefield,field(1),nstart=nrec,interpol='zero_order', &
               needed=needed,found=found,unit=unit,debug_flag=.false.)
        ! write(*,*)'METVAL ',trim(namefield),me,field(40),nrec
        ! NB: need to fix validity
        validity='not set'
      end if
    elseif(MET_SHORT)then
      if(MasterProc)then
        nfetch=1
        call GetCDF_short(namefield,meteoname,var_global,GIMAX,IRUNBEG,GJMAX, &
               JRUNBEG,KMAX,nrec,nfetch,scalefactors,unit,validity,needed=needed)
      end if

      !note: var_global is defined only for me=0
      call global2local_short(var_global,var_local,MSG_READ4,GIMAX,GJMAX,&
            KMAX,1,1)

      CALL MPI_BCAST(scalefactors,8*2,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
      CALL MPI_BCAST(validity,50,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
      CALL MPI_BCAST(unit,50,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
      !scalefactors=1.0
      !validity=' '
      !unit=' '
      if(present(found))found=(validity/=field_not_found)

      if(KMAX==1)then
        ijk=0
        k=1
        do j=1,LJMAX
          do i=1,LIMAX
            ijk=ijk+1
            field(ijk)=var_local(i,j,k)*scalefactors(1)+scalefactors(2)
          end do
        end do
      else
        if(External_Levels_Def)then
          !interpolate vertically if the levels are not identical
          ijk=0
          do k=1,KMAX_MID
            k1=k1_met(k)
            k2=k2_met(k)
            do j=1,LJMAX
              do i=1,LIMAX
                ijk=ijk+1
                field(ijk)=(x_k1_met(k)*var_local(i,j,k1)&
                          +(1.0-x_k1_met(k))*var_local(i,j,k2))*scalefactors(1)&
                          +scalefactors(2)
              end do
            end do
          end do
        else
          ijk=0
          do k=1,KMAX_MID!=KMAX
            do j=1,LJMAX
              do i=1,LIMAX
                ijk=ijk+1
                field(ijk)=var_local(i,j,k)*scalefactors(1)+scalefactors(2)
              end do
            end do
          end do
        end if
      end if
    else
      !data are read as real
      !could also use ReadField to interpolate into a different grid!

      nfetch=1
      istart=1
      jstart=1

      if(namefield==met(ix_u_xmj)%name .and. MET_C_GRID)istart=2!set origin at (2,1) for u-wind
      if(namefield==met(ix_v_xmi)%name .and. MET_C_GRID)jstart=2!set origin at (1,2) for v-wind

      reverse_k=.false.
      if(MET_REVERSE_K)reverse_k=.true.

      if(.not.allocated(meteo_3D))then
        allocate(meteo_3D(LIMAX,LJMAX,KMAX_MET))
        meteo_3D=0.0
      end if
      kstart=1
      kend=KMAX
      namefield_met=namefield
      if(namefield=='SMI1' .and. WRF_MET_CORRECTIONS)then
        !has to fetch level 1
        kstart=1
        kend=1
        namefield_met='SMOIS'
      end if
      if(namefield=='SMI3' .and. WRF_MET_CORRECTIONS)then
        !has to fetch level 3
        kstart=3
        kend=3
        namefield_met='SMOIS'
      end if

      call GetCDF_modelgrid(namefield_met,meteoname,meteo_3D,kstart,kend,nrec,&
            nfetch,i_start=istart,j_start=jstart,reverse_k=reverse_k,&
            needed=needed,found=found)

      if(KMAX==1)then
        ijk=0
        k=1
        do j=1,LJMAX
          do i=1,LIMAX
            ijk=ijk+1
            field(ijk)=meteo_3D(i,j,k)
          end do
        end do
      else
        if(External_Levels_Def)then
          ! interpolate vertically if the levels are not identical
          ijk=0
          do k=1,KMAX_MID
            k1=k1_met(k)
            k2=k2_met(k)
            do j=1,LJMAX
              do i=1,LIMAX
                ijk=ijk+1
                field(ijk)=x_k1_met(k)*meteo_3D(i,j,k1)+(1.0-x_k1_met(k))*meteo_3D(i,j,k2)
              end do
            end do
          end do
        else
         ! use same vertical coordinates as meteo
          ijk=0
          do k=1,KMAX_MID
            do j=1,LJMAX
              do i=1,LIMAX
                ijk=ijk+1
                field(ijk)=meteo_3D(i,j,k)
              end do
            end do
          end do

        end if
      end if
    end if
  end subroutine Getmeteofield

subroutine GetCDF_short(varname,fileName,var,GIMAX,IRUNBEG,GJMAX,JRUNBEG &
     ,KMAX,nstart,nfetch,scalefactors,unit,validity,needed)
  !----------------------------------------------------------------------
  ! open and reads CDF file
  !
  ! The nf90 are functions which return 0 if no error occur.
  ! check is only a subroutine which check wether the function returns zero
  !----------------------------------------------------------------------
  implicit none

  character(len=*),intent(in) :: varname,fileName
  character(len=*),intent(out) :: unit,validity
  real,intent(out) :: scalefactors(2)
  integer, intent(in) :: nstart,GIMAX,IRUNBEG,GJMAX,JRUNBEG,KMAX
  integer, intent(inout) ::  nfetch
  integer(kind=2), dimension(GIMAX*GJMAX*KMAX*NFETCH),intent(out) :: var
  logical,intent(in),optional :: needed
  integer :: varID,ndims
  integer :: ncFileID,status
  real :: scale,offset
  character(len=100) :: period_read=' '
  character(len=200),save :: filename_save='notsaved'
  integer,save :: ncFileID_save=-99
  logical :: is_needed=.false.

  validity='                                     ' !initialisation
  period_read='                                     ' !initialisation
  scalefactors(1) = 1.0 !default
  scalefactors(2) = 0.  !default
  call_msg = "GetCDF_short:"//trim(fileName)
  is_needed=.false.;if(present(needed))is_needed=needed

  ndims=3
  if(KMAX==1)ndims=2
  !open an existing netcdf dataset
  if(trim(filename_save)==trim(filename))then
    ncFileID=ncFileID_save
  else
    if(ncFileID_save/=-99)then
      call check(nf90_close(ncFileID_save))
      filename_save='notsaved'
    end if
    call check(nf90_open(path=trim(fileName),mode=nf90_nowrite,ncid=ncFileID))
      ncFileID_save=ncFileID
    filename_save=trim(filename)
  end if

  !get varID:
  status = nf90_inq_varid(ncid=ncFileID,name=trim(varname),varID=VarID)
  if(status/=nf90_noerr)then
    call CheckStop(is_needed,"meteo field not found:"//trim(varname))
    validity=field_not_found
    var=0.0
    ! call check(nf90_close(ncFileID))
    return
  end if

  ! get scale factors
  status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
  if(status == nf90_noerr) scalefactors(1) = scale
  status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
  if(status == nf90_noerr) scalefactors(2) = offset

  ! find unit
  unit='                                                         '
  status = nf90_get_att(ncFileID, VarID, "units", unit )
  if(status /= nf90_noerr) unit='unknown' !default

  ! find validity
  status = nf90_get_att(ncFileID, VarID, "validity", period_read  )
  if(status == nf90_noerr)then
    validity  = trim(period_read)
  else
    status = nf90_get_att(ncFileID, VarID, "period_of_validity", period_read )
    if(status /= nf90_noerr) validity='instantaneous' !default
  end if

  ! if(Nfetch<nrecords)then
  !    write(*,*)'Reading record',nstart,' to ',nstart+nfetch-1
  ! end if

  ! get variable
  if(ndims==2)then
    call check(nf90_get_var(ncFileID, VarID, var,&
          start=(/IRUNBEG,JRUNBEG,nstart/),count=(/ GIMAX,GJMAX,nfetch /)))
  elseif(ndims==3)then
    call check(nf90_get_var(ncFileID, VarID, var,&
          start=(/IRUNBEG,JRUNBEG,1,nstart/),count=(/GIMAX,GJMAX,KMAX,nfetch /)))
  end if
 ! var=0.0
 ! call check(nf90_close(ncFileID))

end subroutine GetCDF_short

subroutine check(status)
  implicit none
  integer, intent (in) :: status
  call CheckStop(status,nf90_noerr,"Error in Met_mod/NetCDF: "//&
     trim(call_msg)//" "//trim(nf90_strerror(status)))
end subroutine check

subroutine Check_Meteo_Date_Type
  !----------------------------------------------------------------------
  ! On first call, check that dates from meteo file correspond to dates requested.
  ! Also defines nhour_first, nrec, Nhh and METSTEP
  ! Also check if data is in format/tyep short
  !----------------------------------------------------------------------
  character(len=len(meteo)) :: meteoname
  integer :: nyear,nmonth,nday
  integer :: status,ncFileID,timeDimID,timeVarID,VarID,xtype
  character (len = 19) ::  Times_string
  integer ::ihh,ndate(4),n1,nseconds(1),n
  real :: ndays(1),Xminutes(24)
  logical :: date_in_days,MasterProc_local
  integer :: NTime_Read,string_length
  real :: TimesInDays(1000)
  real(kind=8), parameter :: &
       halfsecond=0.5/(24.0*3600.0)! used to avoid rounding errors
  real(kind=8) :: ndays_indate

  nyear=startdate(1)
  nmonth=startdate(2)
  nday=startdate(3)

  meteoname=date2string(meteo,startdate,mode='YMDH')
  MasterProc_local = MasterProc.or.(ME_IO==0)

  if(MasterProc_local)then
     status=nf90_open(path=trim(meteoname),mode=nf90_nowrite,ncid=ncFileID)
     call CheckStop(status,nf90_noerr,'meteo file not found: '//trim(meteoname))
     
     NTime_Read=-1
     call ReadTimeCDF(meteoname,TimesInDays,NTime_Read)
     if(NTime_Read>0)then
        call date2nctime(startdate,ndays_indate)
        do n=1,NTime_Read
           if(ndays_indate-TimesInDays(n)<halfsecond) goto 886               
        enddo
        n=1
        call CheckStop('did not find correct meteo date in '//trim(meteoname))
886     continue
        nrec = n
        Nhh=NTime_Read
        nhour_first=0 !hour of the first record
        if(n==1)nhour_first=mod(nint(TimesInDays(1)*24),24)
        write(*,*)n,' nhour_first =',nhour_first
        if(NTime_Read>1)METSTEP = nint((TimesInDays(2)-TimesInDays(1))*24.)
        if(METSTEP<1)then
           if(MasterProc)write(*,*)'WARNING: met timestep less than one hour; not tested ',(TimesInDays(2)-TimesInDays(1))*24.
           METSTEP = 1
           nrec_mult = nint(1.0/((TimesInDays(2)-TimesInDays(1))*24.))
           if(MasterProc)write(*,*)'reading only every ',nrec_mult,' record'
        endif

        call nctime2date(ndate,TimesInDays(n)) 
        if(abs(ndays_indate-TimesInDays(n))<halfsecond)then
           write(*,*)'first meteo date read ',date2string('YYYY-MM-DD hh:mm:ss',ndate)
        else
           if(abs(ndays_indate+METSTEP/24.0-TimesInDays(n))>halfsecond)then
              write(*,*)'WARNING: correct hours not found ',ndays_indate+METSTEP/24.0,TimesInDays(n)
              call CheckStop('could not find correct meteo dates in '//trim(meteoname))
           endif
           write(*,*)'WARNING: correct hour not found, using twice ',date2string('YYYY-MM-DD hh:mm:ss',ndate)
        endif
     else
        if(MasterProc)write(*,*)'time variable not found'
        nhour_first=0 !hour of the first record
        Nhh=8
        nrec=1
        if(MasterProc_local)then
           write(*,*)'Did not check times, and assume nhour_first =',nhour_first
           write(*,*)'Assume  Nhh =',Nhh
           write(*,*)'Assume  nrec =',nrec
        end if
     endif
     call CheckStop(mod(24,METSTEP),"Met_mod: METSTEP must be a fraction of 24")
  
     if(WRF_MET_CORRECTIONS)then
      !check if the "bucket" method is used
      status = nf90_get_att(ncFileID,nf90_global,"BUCKET_MM",wrf_bucket)
      if(status == nf90_noerr)then
        if(wrf_bucket>0.0)then
          found_wrf_bucket = .true.
          write(*,*)'assuming constant bucket size: ',wrf_bucket
        else
          write(*,*)'Not using buckets ',wrf_bucket
        end if
      end if
    end if

    !check wether the meteo fields are defined as short
    call check(nf90_inq_varid(ncid = ncFileID, name = met(ix_ps)%name, varID = VarID))
    call check(nf90_Inquire_Variable(ncFileID,VarID,xtype=xtype))
    if(xtype==NF90_SHORT)then
       write(*,*)'assuming meteo variables type short'
    else
       if(xtype==NF90_FLOAT .or. xtype==NF90_DOUBLE)then
          write(*,*)'assuming meteo variables type float or double'
          MET_SHORT = .false.
       endif
    endif
    call check(nf90_close(ncFileID))
  end if
  if(me_calc>=0)then
    CALL MPI_BCAST(nrec_mult,4*1,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(nrec,4*1,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(nhour_first,4*1,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(Nhh,4*1,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(METSTEP,4*1,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(found_wrf_bucket,1,MPI_LOGICAL,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(MET_SHORT,1,MPI_LOGICAL,0,MPI_COMM_CALC,IERROR)
  else
    CALL MPI_BCAST(nrec_mult,4*1,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    CALL MPI_BCAST(nrec,4*1,MPI_BYTE,0,MPI_COMM_IO,IERROR)
    CALL MPI_BCAST(nhour_first,4*1,MPI_BYTE,0,MPI_COMM_IO,IERROR)
    CALL MPI_BCAST(Nhh,4*1,MPI_BYTE,0,MPI_COMM_IO,IERROR)
    CALL MPI_BCAST(METSTEP,4*1,MPI_BYTE,0,MPI_COMM_IO,IERROR)
    CALL MPI_BCAST(found_wrf_bucket,1,MPI_LOGICAL,0,MPI_COMM_IO,IERROR)
    CALL MPI_BCAST(MET_SHORT,1,MPI_LOGICAL,0,MPI_COMM_CALC,IERROR)
  end if
  if(found_wrf_bucket)then
    if(me_calc>=0)CALL MPI_BCAST(wrf_bucket,8,MPI_BYTE,0,MPI_COMM_CALC,IERROR)
    if(me_calc<0)CALL MPI_BCAST(wrf_bucket,8,MPI_BYTE,0,MPI_COMM_IO,IERROR)
    met(ix_irainc)%read_meteo  = found_wrf_bucket
    met(ix_irainc)%needed      = found_wrf_bucket
    met(ix_irainnc)%read_meteo = found_wrf_bucket
    met(ix_irainnc)%needed     = found_wrf_bucket
  end if
end subroutine Check_Meteo_Date_Type

subroutine read_surf_elevation(ix)

  integer :: ix

  call GetCDF_modelgrid(met(ix)%name,TopoFile,met(ix)%field,&
                        1,1,1,1,needed=met(ix)%needed,found=met(ix)%found)

  if(.not.met(ix)%found) then
     if( me==0 )write(*,*)'WARNING: met topography not found. Approximating elevation using standard map'
     call ReadField_CDF(SurfacePressureFile,'surface_pressure_year',&
          met(ix)%field,1,needed=.true.,interpol='zero_order')
     
     met(ix)%field=1000.0*StandardAtmos_kPa_2_km(0.001*met(ix)%field)
  else
     if( me==0 )write(*,*)'Read met topography '//trim(met(ix)%name)//' from '//trim(TopoFile)//', typical value:',met(ix)%field(5,5,1,1)
  endif

end subroutine read_surf_elevation

end module met_mod
