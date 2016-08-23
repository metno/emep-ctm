! <Nest_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2015 met.no
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
module Nest_ml
! This module performs the reading or writing of data for nested runs
!
! The Nesting modes (MODE in Nest_config nml) are:
! 0=donothing , 1=write , 2=read , 3=read and write
! 10=write at end of run, 11=read at start , 12=read at start and write at end (BIC)
!
! To make a nested run:
! 1) run with MODE=1 (MODE in Nest_config nml) to write out 3d BC (name in filename_write defined below)
! 2) copy or link filename_write to filename_read_BC (for example "ln -s EMEP_OUT.nc EMEP_IN.nc")
! 3) run (in a smaller domain) with MODE=2
!
! Set MODE (in Nest_config nml) and istart,jstart,iend,jend (same namelist)
! Choose NHOURSAVE and NHOURREAD
! Also filename_read_BC and filename_read_3D should point to appropriate files
! Be careful to remove old BC files before making new ones.
!
! Grids may have any projection.
! Horizontal interpolation uses a weighted average of the four closest points
! This will work also if points in the present grid are not covered by the external grid.
! Vertical interpolation is done from hybrid coordinates.
!
!To do:
!  It should be possible to save only xn_adv_bnd if the inner grid is known for the outer grid.
!  The routines should be thought together with GlobalBC_ml (can it replace it?)

!----------------------------------------------------------------------------!
! External Boundary (BC) and Initial Conditions (IC)
!   ExternalBICs_ml should handle different for different external sources.
!   Experiment specific information must be set on ExternalBICs namelists.
!   So far coded for FORECAST and EnsClimRCA(?) work.
use ExternalBICs_ml,     only: set_extbic, icbc, ICBC_FMT,&
       EXTERNAL_BIC_SET, EXTERNAL_BC, EXTERNAL_BIC_NAME, TOP_BC, &
       iw, ie, js, jn, kt, &! i West/East bnd; j North/South bnd; k Top
       filename_eta,BC_DAYS
!----------------------------------------------------------------------------!
use CheckStop_ml,           only: CheckStop,check=>CheckNC
use Chemfields_ml,          only: xn_adv    ! emep model concs.
use ChemSpecs,              only: NSPEC_ADV, NSPEC_SHL, species_adv
use Functions_ml,           only: great_circle_distance
use GridValues_ml,          only: A_mid,B_mid, glon,glat, i_fdom,j_fdom
use Io_ml,                  only: open_file,IO_TMP,IO_NML,PrintLog
use InterpolationRoutines_ml,  only : grid2grid_coeff
use ModelConstants_ml,      only: Pref,PPB,PT,KMAX_MID, MasterProc, NPROC,  &
    IOU_INST,IOU_HOUR,IOU_YEAR,IOU_MON,IOU_DAY, RUNDOMAIN,  &
    FORECAST,USE_POLLEN, DEBUG_NEST,DEBUG_ICBC=>DEBUG_NEST_ICBC
use MetFields_ml,           only: roa
use netcdf,                 only: nf90_open,nf90_close,nf90_inq_dimid,&
                                  nf90_inquire_dimension,nf90_inq_varid,&
                                  nf90_inquire_variable,nf90_get_var,nf90_get_att,&
                                  nf90_noerr,nf90_nowrite,nf90_global
use netcdf_ml,              only: GetCDF,Out_netCDF,Init_new_netCDF,&
                                  CDFtype=>Real4,ReadTimeCDF,max_filename_length
use OwnDataTypes_ml,        only: Deriv,TXTLEN_SHORT
use Par_ml,                 only: MAXLIMAX,MAXLJMAX,GIMAX,GJMAX,IRUNBEG,JRUNBEG, &
                                  me, li0,li1,lj0,lj1,limax,ljmax
use Pollen_const_ml,        only: pollen_check
use TimeDate_ml,            only: date,current_date,nmdays
use TimeDate_ExtraUtil_ml,  only: idate2nctime,nctime2idate,&
                                  date2string,date2file,compare_date
use Units_ml,               only: Units_Scale
use SmallUtils_ml,          only: find_index
use ChemGroups_ml,          only: chemgroups
implicit none

INCLUDE 'mpif.h'
INTEGER :: INFO

! Nesting modes:
! produces netcdf dump of concentrations if wanted, or initialises mode runs
! from such a file. Used in Nest_ml:
!   0=donothing; 1=write; 2=read; 3=read and write;
!  10=write at end of run; 11=read at start; 12=read atstart and write at end (BIC)
integer, public, save :: MODE

! Nested input/output on FORECAST mode
integer,private,parameter :: FORECAST_NDUMP_MAX = 4  ! Number of nested output
integer, public, save     :: FORECAST_NDUMP     = 1  ! Read by Unimod.f90
! on FORECAST mode (1: start next forecast; 2-4: NMC statistics)
type(date), public :: outdate(FORECAST_NDUMP_MAX)=date(-1,-1,-1,-1,-1)

!coordinates of subdomain to write, relative to FULL domain (only used in write mode)
integer, public ::istart,jstart,iend,jend ! Set on Nest_config namelist

!/-- subroutines

public  :: readxn
public  :: wrtxn

private

logical, private, save :: mydebug =  .false.
integer, private, save :: NHOURSAVE,NHOURREAD ! write/read frequency
!if(NHOURREAD<NHOURSAVE) the data is interpolated in time

character(len=max_filename_length),public, save ::  &  
  template_read_3D = 'EMEP_IN.nc',&       ! Different paths can be set here
  template_read_BC = 'EMEP_IN.nc',&       ! for each of the IO IC/BC files,
  template_write   = 'EMEP_OUT.nc'        ! on Nest_config namelist, if needed.
character(len=max_filename_length),private, save ::  &  
  filename_read_3D = 'template_read_3D',& ! Overwritten in readxn and wrtxn.
  filename_read_BC = 'template_read_BC',& ! Filenames are updated according to date
  filename_write   = 'template_write'     ! following respective templates
logical,private, save ::  &               ! if IC/BC are in the same model/run
  native_grid_3D = .false.,&              ! grid, the expensive call to
  native_grid_BC = .false.,&              ! grid2grid_coeff in init_nest can be avoided
  omit_zero_write= .false.                ! skip const=0.0 variables

! Limit output, e.g. for NMC statistics (3DVar)
character(len=TXTLEN_SHORT), private, save, dimension(10) :: &
  WRITE_SPC = "", &   ! If these varables remain ""
  WRITE_GRP = ""      ! all advected species will be written out.

real(kind=8), parameter :: halfsecond=0.5/(24.0*3600.0)!used to avoid rounding errors
!BC values at boundaries in present grid
real, save, allocatable, dimension(:,:,:,:) :: &
  xn_adv_bndw, xn_adv_bnde, & ! west and east
  xn_adv_bnds, xn_adv_bndn, & ! north and south
  xn_adv_bndt                 ! top

real, save, allocatable, dimension(:) :: ndays_ext !time stamps in days since 1900. NB: only defined on MasterProc

!dimension of external grid for BC
integer,save :: N_ext_BC  !NB: only defined on MasterProc
integer,save :: KMAX_ext_BC

integer,save :: itime!itime_saved(2),
real(kind=8),save :: rtime_saved(2)

integer,save :: date_nextfile(4)!date corresponding to the next BC file to read
integer,save :: NHOURS_Stride_BC   !number of hours between start of two consecutive records in BC files
integer, public, parameter :: NHOURS_Stride_BC_default=6 !time between records if only one record per file (RCA for example)

type(icbc), private, target, dimension(NSPEC_ADV) :: &
  adv_ic=icbc('none','none',1.0,.false.,.false.,-1)  ! Initial 3D IC/BC: spcname,varname,wanted,found,ixadv
type(icbc), private, pointer, dimension(:) :: &
  adv_bc=>null()                                     ! Time dependent BC: spcname,varname,wanted,found,ixadv

contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine Config_Nest()
  integer :: ios,i
  logical, save :: first_call=.true.
  NAMELIST /Nest_config/ MODE,NHOURSAVE,NHOURREAD, &
    template_read_3D,template_read_BC,template_write,&
    native_grid_3D,native_grid_BC,omit_zero_write,istart,jstart,iend,jend,&
    WRITE_SPC,WRITE_GRP,FORECAST_NDUMP,outdate

  if(.not.first_call)return
  mydebug = DEBUG_NEST.and.MasterProc
! Default Nest mode
  MODE=0        ! do nothing (unless FORECAST mode)
! write/read frequency: Hours between consecutive saves(wrtxn)/reads(readxn)
  NHOURSAVE=3   ! Between wrtxn calls.  Should be fraction of 24
  NHOURREAD=1   ! Between readxn calls. Should be fraction of 24
! Default domain for write modes 1,3. Modes 10,12 write full RUNDOMAIN regardles
  istart=RUNDOMAIN(1)+1;iend=RUNDOMAIN(2)-1
  jstart=RUNDOMAIN(3)+1;jend=RUNDOMAIN(4)-1
  rewind(IO_NML)
  read(IO_NML,NML=Nest_config,iostat=ios)
  call CheckStop(ios,"NML=Nest_config")  
  if(mydebug)then
    write(*,*) "NAMELIST IS "
    write(*,NML=Nest_config)
  endif
! write/read frequency should be fraction of 24
  if(MasterProc)then
    call CheckStop(mod(24,NHOURSAVE),"Config_Nest: NHOURSAVE should be fraction of 24")
    call CheckStop(mod(24,NHOURREAD),"Config_Nest: NHOURREAD should be fraction of 24")
  endif
! Update filenames according to date following templates defined on Nest_config
  call init_icbc(cdate=current_date)
! Ensure sub-domain is not larger than run-domain
  istart=max(istart,RUNDOMAIN(1));iend=min(iend,RUNDOMAIN(2))
  jstart=max(jstart,RUNDOMAIN(3));jend=min(jend,RUNDOMAIN(4))
! Ensure that only FORECAST_NDUMP are taking into account
  if(FORECAST.and.(FORECAST_NDUMP<FORECAST_NDUMP_MAX))&
    outdate(FORECAST_NDUMP+1:FORECAST_NDUMP_MAX)%day=0
  if(MasterProc.and.FORECAST)&
    write (*,"(1X,A,10(1X,A,:,','))")'Forecast nest/dump at:',&
     (date2string("YYYY-MM-DD hh:mm:ss",outdate(i)),i=1,FORECAST_NDUMP)
  first_call=.false.
endsubroutine Config_Nest
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine readxn(indate)
  type(date), intent(in) :: indate           ! Gives year..seconds

  integer :: n,i,j,k,KMAX_BC,bc !nseconds(1),n1,II,JJ
  integer :: ndate(4) !nstart,nfetch,nseconds_indate
  real(kind=8):: ndays_indate

  !    real , dimension(48,48,20) ::data
  real :: W1,W2
  logical, save :: first_call=.true.
  logical :: fexist_3D=.false.,fexist_BC=.false.
  integer, save :: oldmonth=0

  call Config_Nest()
  if(mydebug) write(*,*)'Nest:Read BC, MODE=',MODE
  if(.not.any(MODE==[2,3,11,12]).and..not.FORECAST)return

  KMAX_BC=KMAX_MID
  ndate(1:4)=[indate%year,indate%month,indate%day,indate%hour]
  call idate2nctime(ndate,ndays_indate)
  if(first_call)date_nextfile=ndate

  select case(MODE)
  case(100) ! monthly input file
    if(indate%month==oldmonth)return
    if(MasterProc.and.oldmonth==0) write(*,*)'Nest: Initialzing IC'
    oldmonth=indate%month
    if(MasterProc) write(*,*)'Nest: New month, reset BC'
  case(11,12)
    if(.not.first_call)return
    first_call=.false.
    filename_read_3D=date2string(template_read_3D,ndate,debug=mydebug)
    if(MasterProc) write(*,*)'Nest RESET ALL XN 3D ',trim(filename_read_3D)
    call reset_3D(ndays_indate)
    return
  case default
   !if(MasterProc) print *,'call to READXN',indate%hour,indate%seconds
    if(mod(indate%hour,NHOURREAD)/=0.or.indate%seconds/=0)return
  endselect
  ! never comes to this point if MODE=100, 11 or 12

  if(DEBUG_NEST.and.MasterProc) write(*,*) 'Nest: kt', kt, first_call

! Update filenames according to date following templates defined on Nest_config nml
  if(FORECAST)then
    filename_read_3D=date2string(template_read_3D,ndate,debug=mydebug)
    filename_read_BC=date2file  (template_read_BC,ndate,BC_DAYS,"days",debug=mydebug)  
    inquire(file=filename_read_3D,exist=fexist_3D)
    inquire(file=filename_read_BC,exist=fexist_BC)
  else
    filename_read_3D=date2string(template_read_3D,ndate,debug=mydebug)
    filename_read_BC=date2string(template_read_BC,date_nextfile,debug=mydebug)
    fexist_3D=.true.  ! assume 3D file exists
    fexist_BC=.true.  ! assume BC file exists
  endif

  if(first_call)then
    first_call=.false.
    if(fexist_3D)then
      if(MasterProc)write(*,*)'Nest RESET ALL XN 3D ',trim(filename_read_3D)
      call reset_3D(ndays_indate)
    else
      if(MasterProc)write(*,*)'No Nest IC file found: ',trim(filename_read_3D)
    endif

! the first hour only these values are used, no real interpolation between two records
    if(fexist_3D)then
      if(mydebug) write(*,*)'Nest: READING FIRST BC DATA from ',&
            trim(filename_read_BC), ndays_indate
      call read_newdata_LATERAL(ndays_indate)
      if(mydebug) write(*,"(a,5i4)")'Nest: iw, ie, js, jn, kt ',iw,ie,js,jn,kt
    endif
  endif
  if(.not.fexist_BC)then
    if(MasterProc)write(*,*)'No Nest BC file found: ',trim(filename_read_BC)
    return
  endif

  if(ndays_indate-rtime_saved(2)>halfsecond.or.MODE==100)then
   !look for a new data set
    if(MasterProc) write(*,*)'Nest: READING NEW BC DATA from ',&
          trim(filename_read_BC)
    call read_newdata_LATERAL(ndays_indate)
  endif

!   make weights for time interpolation
  if(MODE==100)then   ! don't interpolate for now
    W1=0.0;  W2=1.0   ! use last read value
  else
    W1=1.0;  W2=0.0   ! default: use first read value
    if(rtime_saved(2)-rtime_saved(1)<halfsecond)then
      W1=0.0;  W2=1.0 ! use last read value
    elseif(ndays_indate-rtime_saved(1)>halfsecond)then
      W2=(ndays_indate-rtime_saved(1))/(rtime_saved(2)-rtime_saved(1))
      W1=1.0-W2       ! interpolate
    endif
  endif
  if(DEBUG_NEST.and.MasterProc) then
    write(*,*) 'Nesting BC 2D: time weights : ',W1,W2
    write(*,*) 'Nesting BC 2D: time stamps : ',rtime_saved(1),rtime_saved(2)
  endif

  do bc=1,size(adv_bc)
    if(.not.(adv_bc(bc)%wanted.and.adv_bc(bc)%found))cycle
    n=adv_bc(bc)%ixadv
    if(DEBUG_ICBC.and.MasterProc) write(*,"(2(A,1X),I0,'-->',I0)") &
      'NestICBC: Nesting component',trim(adv_bc(bc)%varname),bc,n
    forall (i=iw:iw, k=1:KMAX_BC, j=1:ljmax, i>=1) &
      xn_adv(n,i,j,k)=W1*xn_adv_bndw(n,j,k,1)+W2*xn_adv_bndw(n,j,k,2)
    forall (i=ie:ie, k=1:KMAX_BC, j=1:ljmax, i<=limax) &
      xn_adv(n,i,j,k)=W1*xn_adv_bnde(n,j,k,1)+W2*xn_adv_bnde(n,j,k,2)
    forall (j=js:js, k=1:KMAX_BC, i=1:limax, j>=1) &
      xn_adv(n,i,j,k)=W1*xn_adv_bnds(n,i,k,1)+W2*xn_adv_bnds(n,i,k,2)
    forall (j=jn:jn, k=1:KMAX_BC, i=1:limax, j<=ljmax) &
      xn_adv(n,i,j,k)=W1*xn_adv_bndn(n,i,k,1)+W2*xn_adv_bndn(n,i,k,2)
    forall (k=kt:kt, i=1:limax, j=1:ljmax, k>=1) &
      xn_adv(n,i,j,k)=W1*xn_adv_bndt(n,i,j,1)+W2*xn_adv_bndt(n,i,j,2)
  enddo

  call CheckStop(EXTERNAL_BIC_NAME=="RCA",&
    "WORK NEEDED: RCA BICs commented out in Nest_ml - not consistent with all chem schemes")
endsubroutine readxn

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine wrtxn(indate,WriteNow)
  type(date), intent(in) :: indate
  logical, intent(in) :: WriteNow !Do not check indate value
  real,allocatable, dimension(:,:,:) :: data ! Data arrays

  type(Deriv) :: def1 ! definition of fields
  integer :: n,iotyp,ndim,kmax,i,ncfileID
  real :: scale
  logical :: fexist, wanted
  logical, save :: first_call=.true.

  call Config_Nest()
  if(.not.any(MODE==[1,3,10,12]).and..not.FORECAST)return

  select case(MODE)
  case(10,12)
    if(.not.WriteNow)return
    istart=RUNDOMAIN(1)
    jstart=RUNDOMAIN(3)
    iend=RUNDOMAIN(2)
    jend=RUNDOMAIN(4)
  case default
    if(FORECAST)then
      outdate(:)%seconds=0   ! output only at full hours
      if(.not.compare_date(FORECAST_NDUMP,indate,outdate(:FORECAST_NDUMP),&
                           wildcard=-1))return
      if(MasterProc) write(*,*)&
        date2string(" Forecast nest/dump at YYYY-MM-DD hh:mm:ss",indate)
      istart=RUNDOMAIN(1)
      jstart=RUNDOMAIN(3)
      iend=RUNDOMAIN(2)
      jend=RUNDOMAIN(4)
    else
      if(mod(indate%hour,NHOURSAVE)/=0.or.indate%seconds/=0)return
    endif
  endselect

  iotyp=IOU_INST
  ndim=3 !3-dimensional
  kmax=KMAX_MID
  scale=1.0
  def1%class='Advected' ! written
  def1%avg=.false.      ! not used
  def1%index=0          ! not used
  def1%scale=scale      ! not used
  def1%iotype=iotyp     ! not used
  def1%name=''          ! written
  def1%unit='mix_ratio' ! written
 
! Update filenames according to date following templates defined on Nest_config nml
! e.g. set template_write="EMEP_BC_MMYYYY.nc" on namelist for different names each month
  filename_write=date2string(template_write,indate,debug=mydebug)
  if(MasterProc)then
    inquire(file=fileName_write,exist=fexist)
    write(*,*)'Nest:write data ',trim(fileName_write),fexist
  endif
  CALL MPI_BCAST(fexist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,INFO)

! Limit output, e.g. for NMC statistics (3DVar)
  if(first_call)then
    first_call=.false.
    call init_icbc(cdate=indate)
    if(any([WRITE_GRP,WRITE_SPC]/=""))then
      adv_ic(:)%wanted=.false.
      do n=1,size(WRITE_GRP)
        if(WRITE_GRP(n)=="")cycle
        i=find_index(WRITE_GRP(n),chemgroups(:)%name)
        if(i>0)then
          where(chemgroups(i)%ptr>NSPEC_SHL) &
            adv_ic(chemgroups(i)%ptr-NSPEC_SHL)%wanted=.true.
        elseif(MasterProc)then
          write(*,"(A,':',/2(2X,A,1X,'''',A,'''',1X,A,'.'))")&
           "Warning (wrtxn)", &
           "Wanted group",trim(WRITE_GRP(n)),"was not found", &
           "Can not be written to file:",trim(filename_write),""
        endif
      enddo
      do n=1,size(WRITE_SPC)
        if(WRITE_SPC(n)=="")cycle
        i=find_index(WRITE_SPC(n),species_adv(:)%name)
        if(i>0)then
          adv_ic(i)%wanted=.true.
        elseif((DEBUG_NEST.or.DEBUG_ICBC).and.MasterProc)then
          write(*,"(A,':',/2(2X,A,1X,'''',A,'''',1X,A,'.'))")&
           "Warning (wrtxn)", &
           "Wanted specie",trim(WRITE_SPC(n)),"was not found", &
           "Can not be written to file:",trim(filename_write),""
        endif
      enddo
    elseif(FORECAST.and.USE_POLLEN)then
      ! POLLEN group members are written to pollen restart/dump file
      call pollen_check(igrp=i)
      if(i>0)then
        where(chemgroups(i)%ptr>NSPEC_SHL) &
          adv_ic(chemgroups(i)%ptr-NSPEC_SHL)%wanted=.false.
        if((DEBUG_NEST.or.DEBUG_ICBC).and.MasterProc)&
          write(*,"(A,':',/2(2X,A,1X,'''',A,'''',1X,A,'.'))")&
           "Warning (wrtxn)", &
           "Group","POLLEN","is written to pollen restart/dump file", &
           "Will not be written to file:",trim(filename_write),""
      endif
    endif
    do n=1,NSPEC_ADV
      if(.not.adv_ic(n)%wanted)then
        if((DEBUG_NEST.or.DEBUG_ICBC).and.MasterProc)&
          write(*,"(A,':',/2(2X,A,1X,'''',A,'''',1X,A,'.'))")&
            "Nest(wrtxn) DEBUG_ICBC",&
            "Variable",trim(species_adv(n)%name),"is not wanted as IC",&
            "Will not be written to file:",trim(filename_write),""
      elseif(omit_zero_write)then !  further reduce output
        wanted=any(xn_adv(n,:,:,:)/=0.0)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,wanted,1,MPI_LOGICAL,MPI_LOR,&
                           MPI_COMM_WORLD,INFO)
        adv_ic(n)%wanted=wanted
        if(.not.adv_ic(n)%wanted.and.&
          (DEBUG_NEST.or.DEBUG_ICBC).and.MasterProc)&
          write(*,"(A,':',/2(2X,A,1X,'''',A,'''',1X,A,'.'))")&
            "Nest(wrtxn) DEBUG_ICBC",&
            "Variable",trim(species_adv(n)%name),"was found constant=0.0",&
            "Will not be written to file:",trim(filename_write),""
      endif
    enddo
  endif

  allocate(data(MAXLIMAX,MAXLJMAX,KMAX_MID))
  !do first one loop to define the fields, without writing them (for performance purposes)
  ncfileID=-1 ! must be <0 as initial value
  if(.not.fexist)then
    do n=1,NSPEC_ADV
      if(.not.adv_ic(n)%wanted)cycle
      def1%name=species_adv(n)%name   ! written
!!    data=xn_adv(n,:,:,:)
      call Out_netCDF(iotyp,def1,ndim,kmax,data,scale,CDFtype=CDFtype,&
            ist=istart,jst=jstart,ien=iend,jen=jend,create_var_only=.true.,&
            fileName_given=trim(fileName_write),ncFileID_given=ncFileID)
    enddo
  endif

  do n=1,NSPEC_ADV
    if(.not.adv_ic(n)%wanted)cycle
    def1%name=species_adv(n)%name     ! written
    data=xn_adv(n,:,:,:)
    call Out_netCDF(iotyp,def1,ndim,kmax,data,scale,CDFtype=CDFtype,&
          ist=istart,jst=jstart,ien=iend,jen=jend,create_var_only=.false.,&
          fileName_given=trim(fileName_write),ncFileID_given=ncFileID)
  enddo
  if(MasterProc)call check(nf90_close(ncFileID))
  deallocate(data)
endsubroutine wrtxn

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine init_icbc(idate,cdate,ndays,nsecs)
!----------------------------------------------------------------------------!
! Setup IC/BC detailed description.
! ICs are assumed to come from Unimod.
!
! adv_ic            IC detailed description for all adv species
! adv_bc            BC detailed description relevant adv species
! EXTERNAL_BC       External (non Unimod) BC detailed description/setup
! EXTERNAL_BIC_SET  EXTERNAL_BC has been set (adv_bc=>EXTERNAL_BC)
!        otherwise  Assume Unimod BCs        (adv_bc=>adv_ic)
!----------------------------------------------------------------------------!
  integer,   intent(in), optional :: idate(4)
  type(date),intent(in), optional :: cdate
  real(kind=8),intent(in),optional:: ndays
  integer,     intent(in),optional:: nsecs
  logical, save :: first_call=.true.
  integer :: n,dat(4)

  if(.not.first_call)return
  first_call=.false.

! One of the date formats needs to be provided
  call CheckStop(count([present(idate),present(cdate),present(ndays),&
                        present(nsecs)]),1,"init_icbc: wrong date option")

! Update filenames according to date following templates defined on Nest_config nml
  if(present(idate)) dat=idate 
  if(present(cdate)) dat=[cdate%year,cdate%month,cdate%day,cdate%hour]
  if(present(ndays)) call nctime2idate(dat,ndays)
  if(present(nsecs)) call nctime2idate(dat,nsecs)
  call set_extbic(dat)  ! set mapping, EXTERNAL_BC, TOP_BC

  if(.not.EXTERNAL_BIC_SET .and. MODE==0)return !No nesting

  filename_read_3D=date2string(template_read_3D,dat,debug=mydebug)
  filename_read_BC=date2file  (template_read_BC,dat,BC_DAYS,"days",debug=mydebug)  
  filename_write  =date2string(template_write  ,dat,debug=mydebug)

  adv_ic(:)%ixadv=(/(n,n=1,NSPEC_ADV)/)
  adv_ic(:)%spcname=species_adv(:)%name
  adv_ic(:)%varname=species_adv(:)%name
  adv_ic(:)%frac=1.0
  adv_ic(:)%wanted=.true.
  adv_ic(:)%found=find_icbc(filename_read_3D,adv_ic%varname(:))
  if(EXTERNAL_BIC_SET) then
    adv_bc=>EXTERNAL_BC
    adv_bc(:)%found=find_icbc(filename_read_bc,adv_bc%varname(:))
  else
    adv_bc=>adv_ic
  endif
  
  if(MasterProc)then
    do n = 1,size(adv_ic%varname)
      if(.not.adv_ic(n)%found)then
        call PrintLog("WARNING: IC variable '"//trim(adv_ic(n)%varname)//"' not found")
      elseif(DEBUG_NEST.or.DEBUG_ICBC)then 
        write(*,*) "init_icbc filled adv_ic "//trim(adv_ic(n)%varname)
      endif
    enddo
    do n = 1,size(adv_bc%varname)
      if(.not.adv_bc(n)%found)then
        call PrintLog("WARNING: BC variable '"//trim(adv_bc(n)%varname)//"' not found")
      elseif(DEBUG_NEST.or.DEBUG_ICBC)then 
        write(*,*) "init_icbc filled adv_bc "//trim(adv_bc(n)%varname)
      endif
    enddo
  endif
  
  if((DEBUG_NEST.or.DEBUG_ICBC).and.MasterProc)then
    write(*,"(a)") "Nest: DEBUG_ICBC Variables:",&
      trim(filename_read_3D),trim(filename_read_BC)
    write(*,"((1X,A,I3,'->',"//ICBC_FMT//"))") &
      ('Nest: ADV_IC',n,adv_ic(n),n=1,size(adv_ic)),&
      ('Nest: ADV_BC',n,adv_bc(n),n=1,size(adv_bc))
  endif
contains
function find_icbc(filename_read,varname) result(found)
!----------------------------------------------------------------------------!
! Check if variables (varname) are present on file (filename_read)
!----------------------------------------------------------------------------!
  implicit none
  character(len=*), intent(in)               :: filename_read
  character(len=*), dimension(:), intent(in) :: varname
  logical, dimension(size(varname))          :: found
  integer :: ncFileID,varID,status,n

  found(:)=.false.
  if(MasterProc)then
    status=nf90_open(trim(filename_read),nf90_nowrite,ncFileID)
    if(status/=nf90_noerr) then
      print *,'icbc: not found ',trim(filename_read)
    else
      print *,'icbc: reading ',trim(filename_read)
      do n=1,size(varname)
        if(varname(n)=="") cycle
        status=nf90_inq_varid(ncFileID,trim(varname(n)),varID)
        found(n)=(status==nf90_noerr)
      enddo
      call check(nf90_close(ncFileID))
    endif
  endif
  CALL MPI_BCAST(found,size(found),MPI_LOGICAL,0,MPI_COMM_WORLD,INFO)
endfunction find_icbc
endsubroutine init_icbc

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine init_nest(ndays_indate,filename_read,native_grid,IIij,JJij,Weight,&
     k1_ext,k2_ext,weight_k1,weight_k2,N_ext,KMAX_ext,GIMAX_ext,GJMAX_ext)
  logical, parameter :: USE_LAST_HYBRID_LEVELS=.true.
  character(len=*),intent(in) :: filename_read
  logical,intent(in) :: native_grid
  real ,intent(out):: Weight(4,MAXLIMAX,MAXLJMAX)
  integer ,intent(out)::IIij(4,MAXLIMAX,MAXLJMAX),JJij(4,MAXLIMAX,MAXLJMAX)
  integer, intent(out), dimension(*) :: k1_ext,k2_ext
  real, intent(out), dimension(*) :: weight_k1,weight_k2
  integer ,intent(out)::N_ext,KMAX_ext,GIMAX_ext,GJMAX_ext
  real(kind=8) :: ndays_indate
  integer :: ncFileID,timeDimID,varid,status,dimIDs(3) !,timeVarID
  integer :: ndate(4) !nseconds_indate,
  real :: DD,dist(4),P_emep
  integer :: i,j,k,n,k_ext,II,JJ !nseconds(1),n,n1,k
  real, allocatable, dimension(:,:) ::lon_ext,lat_ext
  real, allocatable, dimension(:) ::hyam,hybm,P_ext,temp_ll
  character(len=80) ::projection,word,iDName,jDName
  logical :: reversed_k_BC,time_exists,fexist

  rtime_saved = -99999.9 !initialization

  !Read dimensions (global)
  if(MasterProc)then
    status = nf90_open(trim(filename_read),nf90_nowrite,ncFileID)
    if(status/=nf90_noerr) then
      print *,'init_Nest: not found',trim(filename_read)
      return
    else
      print *,'init_Nest: reading ',trim(filename_read)
    endif

    projection='Unknown'
    status = nf90_get_att(ncFileID,nf90_global,"projection",projection)
    if(status==nf90_noerr) then
      if(projection=='lon_lat')projection='lon lat'
      write(*,*)'Nest: projection: '//trim(projection)
    else
      projection='lon lat'
      write(*,*)'Nest: projection not found for ',&
           trim(filename_read)//', assuming '//trim(projection)
    endif
    !get dimensions id/name/len: include more dimension names, if necessary
    GIMAX_ext=get_dimLen(["i        ","lon      ","longitude"],name=iDName)
    GJMAX_ext=get_dimLen(["j       ","lat     ","latitude" ],name=jDName)
    KMAX_ext =get_dimLen(["k    ","mlev ","lev  ","level"])

    select case(projection)
    case('Stereographic')
      call CheckStop("i",iDName,"Nest: unsuported "//&
        trim(iDName)//" as i-dimension on "//trim(projection)//" projection")
      call CheckStop("j",jDName,"Nest: unsuported "//&
        trim(jDName)//" as j-dimension on "//trim(projection)//" projection")
    case('lon lat')
      call CheckStop("lon",iDName(1:3),"Nest: unsuported "//&
        trim(iDName)//" as i-dimension on "//trim(projection)//" projection")
      call CheckStop("lat",jDName(1:3),"Nest: unsuported "//&
        trim(jDName)//" as j-dimension on "//trim(projection)//" projection")
    case default
     !call CheckStop("Nest: unsuported projection "//trim(projection))
     !write(*,*)'GENERAL PROJECTION ',trim(projection)
      call CheckStop("i",iDName,"Nest: unsuported "//&
        trim(iDName)//" as i-dimension on "//trim(projection)//" projection")
      call CheckStop("j",jDName,"Nest: unsuported "//&
        trim(jDName)//" as j-dimension on "//trim(projection)//" projection")
    endselect

    N_ext=0
    status = nf90_inq_dimid(ncFileID,"time",timeDimID)
    time_exists=(status==nf90_noerr)
    if(time_exists) then
      call check(nf90_inquire_dimension(ncFileID,timedimID,len=N_ext))
    else
      status = nf90_inq_dimid(ncFileID,"Months",dimID=timeDimID)
      if(status==nf90_noerr)then
        call check(nf90_inquire_dimension(ncFileID,timedimID,len=N_ext))
        call CheckStop(N_ext,12,'Nest BC: did not find 12 monthes')
      else
        write(*,*)'Nest: time dimension not found. Assuming only one record '
        N_ext=1
      endif
    endif

    write(*,*)'Nest: dimensions external grid',GIMAX_ext,GJMAX_ext,KMAX_ext,N_ext
    if(.not.allocated(ndays_ext))then
      allocate(ndays_ext(N_ext))
    elseif(size(ndays_ext)<N_ext)then
      if(Masterproc)write(*,*)'Nest: Sizes times ',N_ext
      deallocate(ndays_ext)
      allocate(ndays_ext(N_ext))
    endif

  endif
  CALL MPI_BCAST(GIMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(GJMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(KMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  !CALL MPI_BCAST(N_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) !not needed by others than MasterProc

  allocate(lon_ext(GIMAX_ext,GJMAX_ext),lat_ext(GIMAX_ext,GJMAX_ext))
  allocate(hyam(KMAX_ext+1),hybm(KMAX_ext+1),P_ext(KMAX_ext))

  if(MasterProc)then
    !Read lon lat of the external grid (global)
    if(trim(projection)==trim('lon lat')) then
      call check(nf90_inq_varid(ncFileID,iDName,varID),&
        "Read lon-variable: "//trim(iDName))
      allocate(temp_ll(GIMAX_ext))
      call check(nf90_get_var(ncFileID,varID,temp_ll))
      lon_ext=SPREAD(temp_ll,2,GJMAX_ext)
      deallocate(temp_ll)
      call check(nf90_inq_varid(ncFileID,jDName,varID),&
        "Read lat-variable: "//trim(jDName))
      allocate(temp_ll(GJMAX_ext))
      call check(nf90_get_var(ncFileID,varID,temp_ll))
      lat_ext=SPREAD(temp_ll,1,GIMAX_ext)
      deallocate(temp_ll)
    else
      call check(nf90_inq_varid(ncFileID,"lon",varID),"dim:lon")
      call check(nf90_get_var(ncFileID,varID,lon_ext),"get:lon")

      call check(nf90_inq_varid(ncFileID,"lat",varID),"dim:lat")
      call check(nf90_get_var(ncFileID,varID,lat_ext),"get:lat")
    endif

    if(time_exists)then
      call ReadTimeCDF(filename_read,ndays_ext,N_ext)
    else
      !cannot read time on file. assumes it is correct
      ndays_ext(1)=ndays_indate
    endif
    if(MODE==100)then
      !asuming 12 monthes for BC, and 12 or 1 values for IC
      ndays_ext(1)=0
      do n=2,N_ext
        ndays_ext(n)=ndays_ext(n-1)+nmdays(n-1)
      enddo
    elseif(ndays_ext(1)-ndays_indate>halfsecond)then
      call nctime2idate(ndate,ndays_indate,&
          'WARNING: Nest did not find BIC for date YYYY-MM-DD hh:mm:ss')
      call nctime2idate(ndate,ndays_ext(1),&
          'Nest first date found YYYY-MM-DD hh:mm:ss')
    endif

    if(N_ext>1)then
      NHOURS_Stride_BC = nint((ndays_ext(2)-ndays_ext(1))*24)
    else
      !use manually set stride:
      NHOURS_Stride_BC = NHOURS_Stride_BC_default
    endif
    write(*,*)'Nest: new BC record every ',NHOURS_Stride_BC,' hours'

    !enddo
    !Read pressure for vertical levels
    write(*,*)'Nest: reading vertical levels'

    status = nf90_inq_varid(ncFileID,"hyam",varID)
    if(status==nf90_noerr) then
      write(*,*)'Found hyam type levels (values at level midpoints)'
      call check(nf90_inquire_variable(ncFileID,varID,dimIDs=dimIDs))
      call check(nf90_inquire_dimension(ncFileID,dimIDs(1),len=k))
      call CheckStop(k<KMAX_ext,"Nest BC, wrong hyam/hybm dimension")
      if(USE_LAST_HYBRID_LEVELS)then
        k_ext=1+k-KMAX_ext ! for 1+k-KMAX_ext .. k levels
      else
        k_ext=1            ! for 1 .. KMAX_ext levels
      endif
      if(k/=KMAX_ext)&
        write(*,"(A,4(1X,A,I0))")'Nest BC warning:',&
             'kdim #lev=',KMAX_ext,'and hyam/hybm #lev=',k,&
             '. Using only levels ',k_ext,'..',k
      call check(nf90_get_var(ncFileID,varID,hyam,start=(/k_ext/),count=(/KMAX_ext/)))
      status = nf90_get_att(ncFileID,VarID,"units",word)
      if(status==nf90_noerr)then
        if(word(1:3)=='hPa')then
          write(*,*)'Changing hyam from hPa to Pa'
          hyam=100*hyam
        endif
      endif
      call check(nf90_inq_varid(ncFileID,"hybm",varID))
      call check(nf90_get_var(ncFileID,varID,hybm,start=(/k_ext/),count=(/KMAX_ext/)))
    else

      status = nf90_inq_varid(ncFileID,"hyai",varID)
      if(status==nf90_noerr) then
        write(*,*)'Found hyai type levels (values at level interfaces)'

        call check(nf90_inquire_variable(ncFileID,varID,dimIDs=dimIDs))
        call check(nf90_inquire_dimension(ncFileID,dimIDs(1),len=k))
        call CheckStop(k<KMAX_ext+1,"Nest BC, wrong hyai/hybi dimension")
        if(USE_LAST_HYBRID_LEVELS)then
          k_ext=1+(k-1)-KMAX_ext ! for 1+k-KMAX_ext .. k levels
        else
          k_ext=1            ! for 1 .. KMAX_ext levels
        endif
        if(k/=KMAX_ext+1.and.MasterProc)&
          write(*,"(A,4(1X,A,I0))")'Nest BC warning:',&
            'kdim #lev=',KMAX_ext,'and hyam/hybm #lev=',k,&
            '. Using only levels ',k_ext,'..',k
        call check(nf90_get_var(ncFileID,varID,hyam,start=(/k_ext/),count=(/KMAX_ext+1/)))
        status = nf90_get_att(ncFileID,VarID,"units",word)
        if(status==nf90_noerr)then
          if(word(1:3)=='hPa')then
            write(*,*)'Changing hyai from hPa to Pa'
            hyam=100*hyam
          endif
        endif
        call check(nf90_inq_varid(ncFileID,"hybi",varID))
        call check(nf90_get_var(ncFileID,varID,hybm,start=(/k_ext/),count=(/KMAX_ext+1/)))
        do k=1,KMAX_ext
          hyam(k)=0.5*(hyam(k)+hyam(k+1))
          hybm(k)=0.5*(hybm(k)+hybm(k+1))
        enddo

      else
        inquire(file=filename_eta,exist=fexist)
        status = nf90_inq_varid(ncFileID,"k",varID)
        if(status==nf90_noerr) then
          write(*,*)'Nest: assuming sigma level and PT=',PT,KMAX_ext
          call check(nf90_get_var(ncFileID, varID, hybm,count=(/ KMAX_ext /) ))!NB: here assume = sigma
          do k=1,KMAX_ext
            hyam(k)=PT*(1.0-hybm(k))
          enddo
        elseif(fexist) then
          !read eta levels from ad-hoc text file
          write(*,*)'Nest: Reading vertical level from ',trim(filename_eta)
          call open_file(IO_TMP,"r",trim(filename_eta),needed=.true.)
          do n=1,10000
             read(IO_TMP,*)word
             if(trim(word)=='vct')exit
          enddo
          read(IO_TMP,*)(hyam(k),k=1,KMAX_ext+1)!NB: here = A_bnd, not mid
          read(IO_TMP,*)(hybm(k),k=1,KMAX_ext+1)!NB: here = B_bnd, not mid
          close(IO_TMP)
          !convert to mid levels coefficients
          do k=1,KMAX_ext
            hyam(k)=0.5*(hyam(k)+hyam(k+1))
            hybm(k)=0.5*(hybm(k)+hybm(k+1))
          enddo
        else
          status = nf90_inq_varid(ncFileID,"lev",varID)
          if(status == nf90_noerr) then
            call CheckStop('Pressure levels not yet implemented')
            write(*,*)'Nest: assuming pressure levels and hPa'
            call check(nf90_get_var(ncFileID, varID, hyam,count=(/ KMAX_ext /) ))
            hyam=100.0*hyam ! hPa ->Pa
            hybm=0.0
          else
            call CheckStop('Vertical coordinate Unknown/Not yet implemented')
            !assumes lev=1000*(A+B) (IFS-MOZART?)
            !call check(nf90_inq_varid(ncFileID,"lev",varID))
            !call check(nf90_get_var(ncFileID,varID,hybm))
            !hybm=hybm/1000.0
            !hyam=0.0
          endif
        endif
      endif
    endif

    call check(nf90_close(ncFileID))
  endif !end MasterProc

  CALL MPI_BCAST(lon_ext,8*GIMAX_ext*GJMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(lat_ext,8*GIMAX_ext*GJMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(hyam,8*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(hybm,8*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

  ! find horizontal interpolation constants
  ! note that i,j are local and but IIij,JJij refer to the full nest-file
  if(native_grid)then   ! nest-file is on the model/run grid
    forall(i=1:limax,j=1:ljmax)
      IIij(:,i,j)=i_fdom(i)
      JJij(:,i,j)=j_fdom(j)
      Weight(:,i,j)=[1.0,0.0,0.0,0.0]
    endforall
  else                  ! find the four closest points
    call grid2grid_coeff(glon,glat,IIij,JJij,Weight,lon_ext,lat_ext,&
      GIMAX_ext,GJMAX_ext,MAXLIMAX,MAXLJMAX,limax,ljmax,mydebug,1,1)
  endif
! if(MasterProc)print "(A,4(X,I0,X,I0,X,F6.3,','))",'Nest interpol_const',&
!    (IIij(n,1,1),JJij(n,1,1),Weight(n,1,1),n=1,4)
  deallocate(lon_ext,lat_ext)

  !find vertical interpolation coefficients
  !use pressure as reference
  !we want, if possible, P_ext(k1) and P_ext(k2) to be on each side of P_emep
  !We assume constant surface pressure, both for emep and external grid; should not be so
  !   important as long as they are both terrain following.
  do k_ext=1,KMAX_EXT
    P_ext(k_ext)=hyam(k_ext)+hybm(k_ext)*Pref
    if(mydebug) write(*,fmt="(A,I3,F10.2)")'Nest: P_ext',k_ext,P_ext(k_ext)
  enddo
  reversed_k_BC=(P_ext(1)>P_ext(2))
  ! .true.  --> assumes k_ext=KMAX_EXT is top and k_ext=1 is surface
  ! .false. --> assumes k_ext=1 is top and k_ext=KMAX_EXT is surface

  if(reversed_k_BC)then
    do k=1,KMAX_MID
      P_emep=A_mid(k)+B_mid(k)*Pref !Pa
      if(mydebug) write(*,fmt="(A,I3,F10.2)")'Nest: P_emep',k,P_emep
      !largest available P smaller than P_emep (if possible)
      k1_ext(k)=1 !start at surface, and go up until P_emep
      do k_ext=1,KMAX_EXT
        if(P_ext(k_ext)<P_emep)exit
        k1_ext(k)=k_ext
      enddo
      !smallest available P larger than P_emep (if possible)
      k2_ext(k)=KMAX_EXT !start at top, and go down until P_emep
      if(k2_ext(k)==k1_ext(k))k2_ext(k)=KMAX_EXT-1 !avoid k2=k1
      do k_ext=KMAX_EXT,1,-1
         if(P_ext(k_ext)>P_emep)exit
         if(k_ext/=k1_ext(k))k2_ext(k)=k_ext
      enddo
      weight_k1(k)=(P_emep-P_ext(k2_ext(k)))/(P_ext(k1_ext(k))-P_ext(k2_ext(k)))
      weight_k2(k)=1.0-weight_k1(k)
      if(mydebug)&
        write(*,fmt="(A,I3,2(A,I2,A,F5.2))")'Nest: level',k,&
          ' is the sum of level ',k1_ext(k),' weight ',weight_k1(k),&
          ' and level ',k2_ext(k),' weight ',weight_k2(k)
    enddo

  else
    do k=1,KMAX_MID
      P_emep=A_mid(k)+B_mid(k)*Pref !Pa
      if(mydebug) write(*,fmt="(A,I3,F10.2)")'Nest: P_emep',k,P_emep
      !largest available P smaller than P_emep (if possible)
      k1_ext(k)=KMAX_EXT !start at surface, and go up until P_emep
      do k_ext=KMAX_EXT,1,-1
        if(P_ext(k_ext)<P_emep)exit
        k1_ext(k)=k_ext
      enddo
      !smallest available P larger than P_emep (if possible)
      k2_ext(k)=1 !start at top, and go down until P_emep
      if(k2_ext(k)==k1_ext(k))k2_ext(k)=2 !avoid k2=k1
      do k_ext=1,KMAX_EXT
        if(P_ext(k_ext)>P_emep)exit
        if(k_ext/=k1_ext(k))k2_ext(k)=k_ext
      enddo
      weight_k1(k)=(P_emep-P_ext(k2_ext(k)))/(P_ext(k1_ext(k))-P_ext(k2_ext(k)))
      weight_k2(k)=1.0-weight_k1(k)
      if(mydebug) &
        write(*,fmt="(A,I3,2(A,I2,A,F5.2))")'Nest: level',k,&
          ' is the sum of level ', k1_ext(k),' weight ',weight_k1(k),&
           ' and level ', k2_ext(k),' weight ',weight_k2(k)
      enddo
    endif
  deallocate(P_ext,hyam,hybm)

  if(mydebug) &
    write(*,*)'Nest: finished determination of interpolation parameters'
contains
function get_dimLen(dimName,id,name) result(len)
  character(len=*), dimension(:), intent(in) :: dimName
  integer,          optional,     intent(out):: id
  character(len=*), optional,     intent(out):: name
  integer :: d, dID, len

  do d=1,size(dimName)
    status = nf90_inq_dimid(ncFileID,dimName(d),dID)
    if(status==nf90_noerr)then
      if(present(id))  id=dID
      if(present(name))name=trim(dimName(d))
      call check(nf90_inquire_dimension(ncFileID,dID,len=len),&
        "get_dimLen: "//trim(dimName(d)))
      exit
    endif
  enddo
  call CheckStop(status,nf90_noerr,'Nest: '//&
    trim(dimName(1))//'-dimension not found: '//&
    trim(filename_read)//'. Include new name in init_nest')
endfunction get_dimLen
endsubroutine init_nest

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine read_newdata_LATERAL(ndays_indate)
  real(kind=8), intent(in)::ndays_indate
  real, allocatable, dimension(:,:,:) ::data
  integer :: ncFileID,varid,status
  integer :: ndate(4),n,i,j,k,bc
  real    :: unitscale
  real(kind=8) :: ndays(1),ndays_old
  logical, save :: first_call=.true.

  !4 nearest points from external grid  (horizontal)
  integer, save,allocatable :: IIij(:,:,:),JJij(:,:,:)
  !weights of the 4 nearest points (horizontal)
  real, save,allocatable :: Weight(:,:,:)

  !2 adjacent levels from external grid  (vertical)
  integer, allocatable,save, dimension(:) :: k1_ext,k2_ext
  !weights of the 2 adjacent levels (vertical)
  real, allocatable,save, dimension(:) :: weight_k1,weight_k2

  integer:: KMAX_BC!which lvels are interpolated, = KMAX_MID for now
  integer:: timedimID

  !dimensions of external grid for BC
  integer, save ::GIMAX_ext,GJMAX_ext
  character (len=80) ::units
  real :: scale_factor,add_offset
  logical :: time_exists,divbyroa

  KMAX_BC=KMAX_MID
  if(mydebug)write(*,*)'Nest: read_newdata_LATERAL, first?', first_call
  if(first_call)then
    if(mydebug)write(*,*)'Nest: initializations 2D'
    allocate(IIij(4,MAXLIMAX,MAXLJMAX),JJij(4,MAXLIMAX,MAXLJMAX))
    allocate(Weight(4,MAXLIMAX,MAXLJMAX))
    allocate(k1_ext(KMAX_MID),k2_ext(KMAX_MID))
    allocate(weight_k1(KMAX_MID),weight_k2(KMAX_MID))

    call init_icbc(ndays=ndays_indate)
    if(mydebug)write(*,*)'calling init_nest for '//trim(filename_read_BC)
    call init_nest(ndays_indate,filename_read_BC,native_grid_BC,&
                   IIij,JJij,Weight,k1_ext,k2_ext,weight_k1,weight_k2,&
                   N_ext_BC,KMAX_ext_BC,GIMAX_ext,GJMAX_ext)
    if(MODE==100.and.N_ext_BC/=12.and.MasterProc)then
      write(*,*)'Nest: WARNING: Expected 12 monthes in BC file, found ',N_ext_BC
      call CheckStop('Nest BC: wrong number of monthes')
    endif

    !Define & allocate West/East/South/Nort Boundaries
    iw=li0-1;ie=li1+1   ! i West/East   boundaries
    js=lj0-1;jn=lj1+1   ! j South/North boundaries
    kt=0;if(TOP_BC)kt=1 ! k Top         boundary
    if(mydebug)then
      if(kt==1)then
        write(*,*)'Nest-kt test: Also including the top layer in BC'
      else
        write(*,*)'Nest-kt test: Not resetting the top layer'
      endif
    endif

    if(iw>=1    .and..not.allocated(xn_adv_bndw)) &
      allocate(xn_adv_bndw(NSPEC_ADV,MAXLJMAX,KMAX_MID,2)) ! West
    if(ie<=limax.and..not.allocated(xn_adv_bnde)) &
      allocate(xn_adv_bnde(NSPEC_ADV,MAXLJMAX,KMAX_MID,2)) ! East
    if(js>=1    .and..not.allocated(xn_adv_bnds)) &
      allocate(xn_adv_bnds(NSPEC_ADV,MAXLIMAX,KMAX_MID,2)) ! South
    if(jn<=ljmax.and..not.allocated(xn_adv_bndn)) &
      allocate(xn_adv_bndn(NSPEC_ADV,MAXLIMAX,KMAX_MID,2)) ! North
    if(kt>=1    .and..not.allocated(xn_adv_bndt)) &
      allocate(xn_adv_bndt(NSPEC_ADV,MAXLIMAX,MAXLJMAX,2)) ! Top
    if(DEBUG_ICBC)then
      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
      if(MasterProc) write(*, "(A)") "Nest: DEBUG_ICBC Boundaries:"
      write(*,"(1X,'me=',i3,5(1X,A,I0,'=',L1))")&
        me,'W:i',iw,allocated(xn_adv_bndw),'E:i',ie,allocated(xn_adv_bnde),&
           'S:j',js,allocated(xn_adv_bnds),'N:j',jn,allocated(xn_adv_bndn),&
           'T:k',kt,allocated(xn_adv_bndt)
      if(MasterProc)flush(6)
      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    endif
    rtime_saved(2)=-99.0!just to put a value
    if(mydebug)write(*,*)'Nest: end initializations 2D'

  endif

  rtime_saved(1)=rtime_saved(2)!put old values in 1
  allocate(data(GIMAX_ext,GJMAX_ext,KMAX_ext_BC), stat=status)
  if(MasterProc)then
    call check(nf90_open(trim(fileName_read_BC),nf90_nowrite,ncFileID))
    status = nf90_inq_dimid(ncFileID,"time",timeDimID)
    time_exists=(status==nf90_noerr)
    if(time_exists) then
      call check(nf90_inquire_dimension(ncFileID,timedimID,len=N_ext_BC))
    else
      status = nf90_inq_dimid(ncFileID,"Months",timeDimID)
      if(status==nf90_noerr)then
        call check(nf90_inquire_dimension(ncFileID,timedimID,len=N_ext_BC))
        call CheckStop(N_ext_BC,12,'Nest BC: did not find 12 monthes')
      else
        N_ext_BC=1
     endif
    endif

    if(size(ndays_ext)<N_ext_BC)then
      write(*,*)'Nest: New size times in BC file ',N_ext_BC
      deallocate(ndays_ext)
      allocate(ndays_ext(N_ext_BC))
    endif

    if(MODE==100)then
      !only care of the month
      call nctime2idate(ndate,ndays_indate,'Nest: BC reading record MM')
      n=ndate(2)
    else
      if(time_exists)then
        call ReadTimeCDF(filename_read_BC,ndays_ext,N_ext_BC)
        do n=1,N_ext_BC
          if(ndays_indate-ndays_ext(n)<halfsecond) goto 876
        enddo
        n=N_ext_BC
        write(*,*)'Nest: WARNING: did not find correct date ',n
876     continue
      else
        !cannot read time on file. assume it is corresponds to date_nextfile
        n=1
        call idate2nctime(date_nextfile,ndays_ext(n))
      endif
    endif

    call nctime2idate(ndate,ndays_ext(n),'Nest: Reading date YYYY-MM-DD hh:mm:ss')
    if(mydebug) write(*,*)'Nest: Record ',n,' of ',N_ext_BC
    itime=n
    rtime_saved(2)=ndays_ext(n)
    if(n==N_ext_BC)then
      !next data to be read should be from another file
      if(mydebug)then
        write(*,*)'Nest: Last record reached ',n,N_ext_BC
        call nctime2idate(date_nextfile,ndays_ext(n)+NHOURS_Stride_BC/24.0,&
              'next BC date to read:  YYYY-MM-DD hh:mm:ss')
        write(*,*)'Nest: date_nextfile ',date_nextfile
      else
        call nctime2idate(date_nextfile,ndays_ext(n)+NHOURS_Stride_BC/24.0)
      endif
    endif

  endif

  CALL MPI_BCAST(rtime_saved,8*2,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

  if(.not.first_call) call store_old_bc() !store the old values in 1
  if(allocated(xn_adv_bndw)) xn_adv_bndw(:,:,:,2)=0.0
  if(allocated(xn_adv_bnde)) xn_adv_bnde(:,:,:,2)=0.0
  if(allocated(xn_adv_bnds)) xn_adv_bnds(:,:,:,2)=0.0
  if(allocated(xn_adv_bndn)) xn_adv_bndn(:,:,:,2)=0.0
  if(allocated(xn_adv_bndt)) xn_adv_bndt(:,:,:,2)=0.0

  DO_BC: do bc=1,size(adv_bc)
    if(.not.(adv_bc(bc)%wanted.and.adv_bc(bc)%found))cycle DO_BC
    n=adv_bc(bc)%ixadv
    if(MasterProc)then
      if(DEBUG_NEST.or.DEBUG_ICBC) write(*,"(2(A,1X),I0,'-->',I0)")&
        'Nest: DO_BC',trim(adv_bc(bc)%varname),bc,n
    !Could fetch one level at a time if sizes becomes too big
      call check(nf90_inq_varid(ncFileID,trim(adv_bc(bc)%varname),varID))

      call check(nf90_get_var(ncFileID, varID, data &
            ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext_BC,1 /) ))
      status = nf90_get_att(ncFileID,VarID,"scale_factor",scale_factor)
      if(status==nf90_noerr) data=data*scale_factor
      status = nf90_get_att(ncFileID,VarID,"add_offset",add_offset)
      if(status==nf90_noerr) data=data+add_offset
      status = nf90_get_att(ncFileID,VarID,"units",units)
      if(units=="1")then
        if(index(adv_bc(bc)%varname,"vmr")>0)units="vmr"
        if(index(adv_bc(bc)%varname,"mmr")>0)units="mmr"
      endif
      if(status==nf90_noerr) then
        if(DEBUG_NEST.or.DEBUG_ICBC) write(*,*)&
          'Nest: variable '//trim(adv_bc(bc)%varname)//' has unit '//trim(units)
        unitscale=adv_bc(bc)%frac/Units_Scale(units,n,needroa=divbyroa,&
                                              debug_msg="read_newdata_LATERAL")
      else
        if(DEBUG_NEST.or.DEBUG_ICBC) write(*,*)&
          'Nest: variable '//trim(adv_bc(bc)%varname//' has no unit attribute')
        unitscale=adv_bc(bc)%frac
      endif
      if(unitscale/=1.0) data=data*unitscale
    endif
    CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext_BC,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(divbyroa,1,MPI_LOGICAL,0,MPI_COMM_WORLD,INFO)

   !overwrite Global Boundaries (lateral faces)
    if(divbyroa)then
      if(allocated(xn_adv_bndw)) forall(k=1:KMAX_BC, j=1:ljmax) &
        xn_adv_bndw(n,j,k,2)=xn_adv_bndw(n,j,k,2) &
                            +(WeightData(iw,j,k1_ext(k))*weight_k1(k) &
                             +WeightData(iw,j,k2_ext(k))*weight_k2(k))&
                            /roa(iw,j,k,1)
      if(allocated(xn_adv_bnde)) forall(k=1:KMAX_BC, j=1:ljmax) &
        xn_adv_bnde(n,j,k,2)=xn_adv_bnde(n,j,k,2) &
                            +(WeightData(ie,j,k1_ext(k))*weight_k1(k) &
                             +WeightData(ie,j,k2_ext(k))*weight_k2(k))&
                            /roa(ie,j,k,1)
      if(allocated(xn_adv_bnds)) forall(k=1:KMAX_BC, i=1:limax) &
        xn_adv_bnds(n,i,k,2)=xn_adv_bnds(n,i,k,2) &
                            +(WeightData(i,js,k1_ext(k))*weight_k1(k) &
                             +WeightData(i,js,k2_ext(k))*weight_k2(k))&
                            /roa(i,js,k,1)
      if(allocated(xn_adv_bndn)) forall(k=1:KMAX_BC, i=1:limax) &
        xn_adv_bndn(n,i,k,2)=xn_adv_bndn(n,i,k,2) &
                            +(WeightData(i,jn,k1_ext(k))*weight_k1(k) &
                             +WeightData(i,jn,k2_ext(k))*weight_k2(k))&
                            /roa(i,jn,k,1)
      if(allocated(xn_adv_bndt)) forall(i=1:limax, j=1:ljmax) &
        xn_adv_bndt(n,i,j,2)=xn_adv_bndt(n,i,j,2) &
                            +(WeightData(i,j,k1_ext(kt))*weight_k1(kt) &
                             +WeightData(i,j,k2_ext(kt))*weight_k2(kt))&
                            /roa(i,j,kt,1)
    else
      if(allocated(xn_adv_bndw)) forall(k=1:KMAX_BC, j=1:ljmax) &
        xn_adv_bndw(n,j,k,2)=xn_adv_bndw(n,j,k,2) &
                            +WeightData(iw,j,k1_ext(k))*weight_k1(k)&
                            +WeightData(iw,j,k2_ext(k))*weight_k2(k)
      if(allocated(xn_adv_bnde)) forall(k=1:KMAX_BC, j=1:ljmax) &
        xn_adv_bnde(n,j,k,2)=xn_adv_bnde(n,j,k,2) &
                            +WeightData(ie,j,k1_ext(k))*weight_k1(k)&
                            +WeightData(ie,j,k2_ext(k))*weight_k2(k)
      if(allocated(xn_adv_bnds)) forall(k=1:KMAX_BC, i=1:limax) &
        xn_adv_bnds(n,i,k,2)=xn_adv_bnds(n,i,k,2) &
                            +WeightData(i,js,k1_ext(k))*weight_k1(k)&
                            +WeightData(i,js,k2_ext(k))*weight_k2(k)
      if(allocated(xn_adv_bndn)) forall(k=1:KMAX_BC, i=1:limax) &
        xn_adv_bndn(n,i,k,2)=xn_adv_bndn(n,i,k,2) &
                            +WeightData(i,jn,k1_ext(k))*weight_k1(k)&
                            +WeightData(i,jn,k2_ext(k))*weight_k2(k)
      if(allocated(xn_adv_bndt)) forall(i=1:limax, j=1:ljmax) &
        xn_adv_bndt(n,i,j,2)=xn_adv_bndt(n,i,j,2) &
                            +WeightData(i,j,k1_ext(kt))*weight_k1(kt)&
                            +WeightData(i,j,k2_ext(kt))*weight_k2(kt)
    endif
  enddo DO_BC

  if(first_call)then
    !copy 2 into 1 so that both are well defined
    rtime_saved(1)=rtime_saved(2)!put  time in 1
    call store_old_bc() !store the old values in 1
  endif

  deallocate(data)
  if(MasterProc) call check(nf90_close(ncFileID))
  first_call=.false.
  return
  contains
  PURE function WeightData(i,j,k) result(wsum)
    integer, intent(in)::i,j,k
    real:: wsum
    wsum=dot_product(Weight(:,i,j),&
      (/data(IIij(1,i,j),JJij(1,i,j),k),data(IIij(2,i,j),JJij(2,i,j),k),&
        data(IIij(3,i,j),JJij(3,i,j),k),data(IIij(4,i,j),JJij(4,i,j),k)/))
  endfunction WeightData
  subroutine store_old_bc !store the old values in 1
    if(allocated(xn_adv_bndw)) xn_adv_bndw(:,:,:,1)=xn_adv_bndw(:,:,:,2)
    if(allocated(xn_adv_bnde)) xn_adv_bnde(:,:,:,1)=xn_adv_bnde(:,:,:,2)
    if(allocated(xn_adv_bnds)) xn_adv_bnds(:,:,:,1)=xn_adv_bnds(:,:,:,2)
    if(allocated(xn_adv_bndn)) xn_adv_bndn(:,:,:,1)=xn_adv_bndn(:,:,:,2)
    if(allocated(xn_adv_bndt)) xn_adv_bndt(:,:,:,1)=xn_adv_bndt(:,:,:,2)
  endsubroutine store_old_bc
endsubroutine read_newdata_LATERAL

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine reset_3D(ndays_indate)
  implicit none
  real(kind=8), intent(in)::ndays_indate
  real, allocatable, dimension(:,:,:) ::data
  integer :: ndate(4),n,i,j,k,itime=0,status
  integer :: ncFileID,varid
  real    :: unitscale
  real(kind=8) :: ndays(1)
  logical, save :: first_call=.true.

  !4 nearest points from external grid
  integer, save,allocatable :: IIij(:,:,:),JJij(:,:,:)

  !weights of the 4 nearest points
  real, save,allocatable :: Weight(:,:,:)

  !dimensions of external grid for 3D
  integer, save ::N_ext,KMAX_ext,GIMAX_ext,GJMAX_ext

  !2 adjacent levels from external grid  (vertical)
  integer, allocatable,save, dimension(:) :: k1_ext,k2_ext
  !weights of the 2 adjacent levels (vertical)
  real, allocatable,save, dimension(:) :: weight_k1,weight_k2

  character (len=80) :: units
  real :: scale_factor,add_offset
  logical :: divbyroa

  if(mydebug) write(*,*) 'Nest: initializations 3D', first_call
 
  if(first_call)then
    if(mydebug) write(*,*)'Nest: initializations 3D'
    allocate(IIij(4,MAXLIMAX,MAXLJMAX),JJij(4,MAXLIMAX,MAXLJMAX))
    allocate(Weight(4,MAXLIMAX,MAXLJMAX))
    allocate(k1_ext(KMAX_MID),k2_ext(KMAX_MID))
    allocate(weight_k1(KMAX_MID),weight_k2(KMAX_MID))
    first_call=.false.
    if(mydebug) write(*,*) 'Nest: init-icbc'
    call init_icbc(ndays=ndays_indate)
    if(mydebug) write(*,*)'calling init_nest for 3D '//trim(filename_read_3D)
    call init_nest(ndays_indate,filename_read_3D,native_grid_3D,&
                  IIij,JJij,Weight,k1_ext,k2_ext,weight_k1,weight_k2,&
                  N_ext,KMAX_ext,GIMAX_ext,GJMAX_ext)
    if(MODE==100.and.(N_ext/=12.and.N_ext/=1.and.MasterProc))then
      write(*,*)'Nest: WARNING: Expected 12 or 1 monthes in IC file, found ',N_ext
      call CheckStop('Nest: IC: wrong number of months')
    endif
    if(mydebug) write(*,*)'Nest: end initializations 3D'
  endif
  allocate(data(GIMAX_ext,GJMAX_ext,KMAX_ext), stat=status)
  if(MasterProc)then
    call check(nf90_open(trim(fileName_read_3D),nf90_nowrite,ncFileID))
    if(MODE==100)then
      if(N_ext==1)then
        n=1
      else
        call nctime2idate(ndate,ndays_indate,'Using record MM')
        n=ndate(2)
      endif
    else
      do n=1,N_ext
        if(ndays_ext(n)>=ndays_indate) goto 876
      enddo
      n=N_ext
      write(*,*)'Nest: WARNING: did not find correct date'
876   continue
      call nctime2idate(ndate,ndays_ext(n),'Using date YYYY-MM-DD hh:mm:ss')
    endif
    itime=n
  endif

  if(mydebug)write(*,*)'Nest: overwrite 3D'
 
  DO_SPEC: do n= 1, NSPEC_ADV
    if(.not.(adv_ic(n)%wanted.and.adv_ic(n)%found)) cycle DO_SPEC
    if(MasterProc)then
      if(DEBUG_NEST) print *,'Nest: 3D component ',trim(adv_ic(n)%varname)
      !Could fetch one level at a time if sizes becomes too big
      call check(nf90_inq_varid(ncFileID,trim(adv_ic(n)%varname),varID))
      call check(nf90_get_var(ncFileID, varID, data &
            ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext,1 /) ))
      status = nf90_get_att(ncFileID,VarID,"scale_factor",scale_factor)
      if(status==nf90_noerr) data=data*scale_factor
      status = nf90_get_att(ncFileID,VarID,"add_offset",add_offset)
      if(status==nf90_noerr) data=data+add_offset
      status = nf90_get_att(ncFileID,VarID,"units",units)
      if(units=="1")then
        if(index(adv_ic(n)%varname,"vmr")>0)units="vmr"
        if(index(adv_ic(n)%varname,"mmr")>0)units="mmr"
      endif
      if(status==nf90_noerr) then
        if(DEBUG_NEST) write(*,*)&
          'Nest: variable '//trim(adv_ic(n)%varname)//' has unit '//trim(units)
        unitscale=adv_ic(n)%frac/Units_Scale(units,n,needroa=divbyroa,&
                                             debug_msg="reset_3D")
      else
        if(DEBUG_NEST) write(*,*)&
          'Nest: variable '//trim(adv_ic(n)%varname//' has no unit attribute')
        unitscale=adv_ic(n)%frac
      endif
      if(unitscale/=1.0) data=data*unitscale
    endif
    CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(divbyroa,1,MPI_LOGICAL,0,MPI_COMM_WORLD,INFO)

     !overwrite everything 3D (init)
    if(divbyroa)then
      forall (k=1:KMAX_MID, j=1:ljmax, i=1:limax) &
        xn_adv(n,i,j,k)=(WeightData(i,j,k1_ext(k))*weight_k1(k) &
                        +WeightData(i,j,k2_ext(k))*weight_k2(k))&
                       /roa(i,j,k,1)
    else
      forall (k=1:KMAX_MID, j=1:ljmax, i=1:limax) &
        xn_adv(n,i,j,k)=WeightData(i,j,k1_ext(k))*weight_k1(k)&
                       +WeightData(i,j,k2_ext(k))*weight_k2(k)
    endif

  enddo DO_SPEC

  deallocate(data)
  if(MasterProc) call check(nf90_close(ncFileID))
  contains
  PURE function WeightData(i,j,k) result(wsum)
    integer, intent(in)::i,j,k
    real:: wsum
    wsum=dot_product(Weight(:,i,j),&
      (/data(IIij(1,i,j),JJij(1,i,j),k),data(IIij(2,i,j),JJij(2,i,j),k),&
        data(IIij(3,i,j),JJij(3,i,j),k),data(IIij(4,i,j),JJij(4,i,j),k)/))
  endfunction WeightData
endsubroutine reset_3D

endmodule Nest_ml
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
