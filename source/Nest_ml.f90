! <Nest_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

module Nest_ml
! This module performs the reading or writing of data for nested runs
!
! The Nesting modes (NEST_MODE in ModelConstants_ml) are:
! 0=donothing , 1=write , 2=read , 3=read and write
! 10=write at end of run, 11=read at start , 12=read at start and write at end (BIC)
!
! To make a nested run:
! 1) run with MODE=1 (NEST_MODE in ModelConstants_ml) to write out 3d BC
! 2) run (in a smaller domain) with MODE=2
!
! Set MODE (NEST_MODE in ModelConstants_ml) and istart,jstart,iend,jend
! Choose NHOURSAVE and NHOURREAD
! Also filename_read_BC_template and filename_read_3D should point to appropriate files
!
! Grids may have any projection.
! Horizontal interpolation uses a weighted average of the four closest points
! This will work also if points in the present grid are not covered by the external grid.
! Vertical interpolation is done from hybrid coordinates.
!
!To do:
!  It should be possible to save only xn_adv_bnd if the inner grid is known for the outer grid.
!  The routines should be thought together with GlobalBC_ml (can it replace it?)

use OwnDataTypes_ml,        only: Deriv
use TimeDate_ml,            only: date
use GridValues_ml,          only: glon,glat
use CheckStop_ml,          only : CheckStop,StopAll
use ChemChemicals_ml,       only: species
use ChemSpecs_shl_ml,       only: NSPEC_SHL
use ChemSpecs_adv_ml
use ChemSpecs_tot_ml,       only: NSPEC_TOT
use GridValues_ml,          only: A_mid,B_mid
use Io_ml,                  only: open_file, IO_TMP
use netcdf
use netcdf_ml,              only: GetCDF,Out_netCDF,Init_new_netCDF&
                                  ,Int1,Int2,Int4,Real4,Real8,ReadTimeCDF
use Functions_ml,           only: great_circle_distance
use ModelConstants_ml,      only: Pref,PPB,PT,KMAX_MID, MasterProc, NPROC     &
  , IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY,RUNDOMAIN  &
  , MODE=>NEST_MODE, FORECAST, DEBUG_NEST, DEBUG_ICBC=>DEBUG_NEST_ICBC
use Par_ml,                 only: MAXLIMAX, MAXLJMAX, GIMAX,GJMAX,IRUNBEG,JRUNBEG &
  , me, li0,li1,lj0,lj1,limax,ljmax, tgi0, tgj0, tlimax, tljmax
use Chemfields_ml,          only: xn_adv, xn_shl    ! emep model concs.
use TimeDate_ExtraUtil_ml,  only: idate2nctime,nctime2idate,date2string

implicit none

INCLUDE 'mpif.h'
INTEGER INFO

! Nested input/output on FORECAST mode
integer, public, parameter :: FORECAST_NDUMP = 2  ! Number of nested output
! on FORECAST mode (1: starnt next forecast; 2: NMC statistics)
type(date), public :: outdate(FORECAST_NDUMP)=date(-1,-1,-1,-1,-1)
logical , public,parameter :: TRANSPHORM = .false.  ! Use limited set of BC components
logical , public,parameter :: RCA = .false.        ! Use limited set of IC and BC components
!RCA: remember to manually set some BC to fixed values!
!coordinates of subdomain to write
!coordinates relative to LARGE domain (only used in write mode)
integer ::istart=60,jstart=11,iend=107,jend=58 !ENEA NB: version has changed, these numbers where for small domain!!!

!/-- subroutines

public  :: readxn
public  :: wrtxn

integer, public, parameter :: NHOURSAVE=3 !time between two saves. should be a fraction of 24
integer, public, parameter :: NHOURREAD=1 !time between two reads. should be a fraction of 24
!if(NHOURREAD<NHOURSAVE) the data is interpolated in time

private

!Use TOP BC on forecast mode
logical, parameter :: TOP_BC=.false..or.FORECAST.or.RCA.or.TRANSPHORM
integer,save :: iw, ie, js, jn, kt ! i West/East bnd; j North/South bnd; k Top
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

!YYYY, YY, MM, DD, hh will be replaced by numbers by the program. For details, see detail2str in TimeDate_ExtraUtil_ml.f90
character(len=130),save  :: filename_read_BC_template='/global/work/mifapw/emep/Data/RCA/Boundary_conditions/YYYY/YYYYMMDDhh00_ll.nc'
!character(len=130),save  :: filename_read_BC_template='/global/work/mifajej/TRANSPHORM/Data/EMAC_YYYYMM_TRANSPHORM.nc'

character(len=130),save  :: filename_read_3D='/global/work/mifapw/emep/Data/RCA/Boundary_conditions/YYYY/YYYYMMDDhh00_ll.nc'
!character(len=130),save  :: filename_read_3D='/global/work/mifajej/TRANSPHORM/MACC/spinnup/EMEP_IN.nc'
character(len=130),save  :: filename_eta='/global/work/mifapw/emep/Data/MACC02/Boundary_conditions/mozart_eta.zaxis'
character(len=30),save  :: filename_write='EMEP_OUT.nc'


character(len=130),save  :: filename_read_BC
integer,save :: date_nextfile(4)!date corresponding to the next BC file to read
integer,save :: NHOURS_Stride_BC   !number of hours between start of two consecutive records in BC files
integer, public, parameter :: NHOURS_Stride_BC_default=6 !time between records if only one record per file (RCA for example)


! Nested output dates on FORECAST mode
! IFS-MOZART BC
type, private :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
  character(len=24) :: varname=""
  logical           :: wanted=.false.,found=.false.
end type icbc
type(icbc), dimension(NSPEC_ADV), private :: &
  adv_ic=icbc('',.false.,.false.), &  ! Initial 3D IC/CB
  adv_bc=icbc('',.false.,.false.)     ! Time dependent BC
type, private :: adv_icbc             ! IC/BC Set, included intended ixadv
  integer           :: ixadv=-1
  type(icbc)        :: icbc
end type adv_icbc
type(adv_icbc), dimension(9), private, parameter :: &  ! BC from IFS-MOZART
  FORECAST_BC=(/adv_icbc(IXADV_O3    ,icbc('O3_VMR_inst'    ,.true.,.false.)), &
                adv_icbc(IXADV_NO    ,icbc('NO_VMR_inst'    ,.true.,.false.)), &
                adv_icbc(IXADV_NO2   ,icbc('NO2_VMR_inst'   ,.true.,.false.)), &
                adv_icbc(IXADV_PAN   ,icbc('PAN_VMR_inst'   ,.true.,.false.)), &
                adv_icbc(IXADV_HNO3  ,icbc('HNO3_VMR_inst'  ,.true.,.false.)), &
                adv_icbc(IXADV_CO    ,icbc('CO_VMR_inst'    ,.true.,.false.)), &
                adv_icbc(IXADV_C2H6  ,icbc('C2H6_VMR_inst'  ,.true.,.false.)), &
                adv_icbc(IXADV_HCHO  ,icbc('CH2O_VMR_inst'  ,.true.,.false.)), &
                adv_icbc(IXADV_CH3CHO,icbc('CH3CHO_VMR_inst',.true.,.false.))/)

type(adv_icbc), dimension(9), private, parameter :: &  ! BC from EMAC-TRANSPHORM
  TRANSPHORM_BC=(/adv_icbc(IXADV_SO2   ,icbc('SO2'    ,.true.,.false.)), &
                 adv_icbc(IXADV_HNO3   ,icbc('HNO3'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NH3    ,icbc('NH3'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NO     ,icbc('NO'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NO2    ,icbc('NO2'    ,.true.,.false.)), &
                 adv_icbc(IXADV_O3     ,icbc('O3'    ,.true.,.false.)), &
                 adv_icbc(IXADV_SO4    ,icbc('SO4_i'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NH4_F  ,icbc('NH4_i'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NO3_F  ,icbc('NO3_i'    ,.true.,.false.))/)

type(adv_icbc), dimension(26), private, parameter :: &  ! BC from EnsClim-RCA
  RCA_BC=(/adv_icbc(IXADV_NO   ,icbc('NO'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NO2   ,icbc('NO2'    ,.true.,.false.)), &
                 adv_icbc(IXADV_O3    ,icbc('O3'    ,.true.,.false.)), &
                 adv_icbc(IXADV_CO     ,icbc('CO'    ,.true.,.false.)), &
                 adv_icbc(IXADV_HCHO    ,icbc('HCHO'    ,.true.,.false.)), &
                 adv_icbc(IXADV_PAN     ,icbc('PAN'    ,.true.,.false.)), &
                 adv_icbc(IXADV_HNO3    ,icbc('HNO3'    ,.true.,.false.)), &
                 adv_icbc(IXADV_H2O2  ,icbc('H2O2'    ,.true.,.false.)), &
                 adv_icbc(IXADV_CH4  ,icbc('CH4'    ,.true.,.false.)), &
                 adv_icbc(IXADV_CH3CHO  ,icbc('CH3CHO'    ,.true.,.false.)), &
                 adv_icbc(IXADV_C2H6  ,icbc('C2H6'    ,.true.,.false.)), &
                 adv_icbc(IXADV_C5H8  ,icbc('Isoprene'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NC4H10  ,icbc('nC4H10'    ,.true.,.false.)), &
                 adv_icbc(IXADV_OXYL  ,icbc('oXylene'    ,.true.,.false.)), &
                 adv_icbc(IXADV_SO2  ,icbc('SO2'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NH3  ,icbc('NH3'    ,.true.,.false.)), &
                 adv_icbc(IXADV_N2O5  ,icbc('N2O5'    ,.true.,.false.)), &
                 adv_icbc(IXADV_SO4  ,icbc('SO4'    ,.true.,.false.)), &!SO4+ NH4HSO4 + (NH4)2SO4
                 adv_icbc(IXADV_SO4  ,icbc('NH4HSO4'    ,.true.,.false.)), &!SO4+ NH4HSO4 + (NH4)2SO4
                 adv_icbc(IXADV_SO4  ,icbc('NH42SO4'    ,.true.,.false.)), &!SO4+ NH4HSO4 + (NH4)2SO4
                 adv_icbc(IXADV_NO3_F  ,icbc('NH4NO3'    ,.true.,.false.)), &!NH4NO3+ANIT 
                 adv_icbc(IXADV_NO3_F  ,icbc('ANIT'    ,.true.,.false.)), &!+ANIT 
                 adv_icbc(IXADV_NH4_F  ,icbc('NH4NO3'    ,.true.,.false.)), &!NH4NO3+ NH4HSO4 + 2 * (NH4)2SO4
                 adv_icbc(IXADV_NH4_F  ,icbc('NH4HSO4'    ,.true.,.false.)), &!NH4NO3+ NH4HSO4 + 2 * (NH4)2SO4
                 adv_icbc(IXADV_NH4_F  ,icbc('NH42SO4'    ,.true.,.false.)), &!NH4NO3+ NH4HSO4 + 2 * (NH4)2SO4
                 adv_icbc(IXADV_NH4_F  ,icbc('NH42SO4'    ,.true.,.false.)) &!NH4NO3+ NH4HSO4 + 2 * (NH4)2SO4
  
/)


contains

subroutine readxn(indate)
  implicit none
  type(date), intent(in) :: indate           ! Gives year..seconds
  integer,save  :: first_data=-1

  integer :: n,i,j,k,KMAX_BC !nseconds(1),n1,II,JJ
  integer :: ndate(4) !nstart,nfetch,nseconds_indate
  real(kind=8):: ndays_indate

  !    real , dimension(48,48,20) ::data
  real :: W1,W2
  logical, save :: first_call=.true.
  logical :: fexist=.false.

  if(DEBUG_NEST.and.MasterProc)write(*,*)'Read BC, MODE=',MODE
  if(MODE /= 2.and.MODE /= 3.and. MODE /= 11.and. MODE /= 12.and. .not.FORECAST)return

  KMAX_BC=KMAX_MID

  ndate(1)  = indate%year
  ndate(2)  = indate%month
  ndate(3)  = indate%day
  ndate(4)  = indate%hour
  call idate2nctime(ndate,ndays_indate)
  if(first_call)date_nextfile=ndate

  if(FORECAST)then ! FORECAST mode superseeds nest MODE
    filename_read_3D='EMEP_IN_IC.nc'          !IC file: dump/re-start
    filename_read_BC_template='EMEP_IN_BC_YYYYMMDD.nc'
    filename_read_BC=date2string(trim(filename_read_BC_template),date_nextfile) 
!    filename_read_BC=date2string('EMEP_IN_BC_YYYYMMDD.nc',indate) !BC file: 01,...,24 UTC rec for 1 day
    if(first_call)then
      first_call=.false.
      inquire(file=filename_read_3D,exist=fexist)
      if(.not.fexist)then
        if(MasterProc) print *,'No nest IC file found: ',trim(filename_read_3D)
      else
        if(MasterProc) print *,'RESET ALL XN 3D'
        call reset_3D(ndays_indate)
      endif
    endif
    if(mod(indate%hour,NHOURREAD)/=0.or.indate%seconds/=0) return
    inquire(file=filename_read_BC,exist=fexist)
    if(.not.fexist)then
      if(MasterProc) print *,'No nest BC file found: ',trim(filename_read_BC)
      return
    endif
  elseif(MODE == 11.or.MODE == 12)then
    if(.not. first_call)return
    first_call=.false.
    if(MasterProc)   print *,'RESET ALL XN 3D'
    call reset_3D(ndays_indate)
    return
  else
   !if(MasterProc) print *,'call to READXN',indate%hour,indate%seconds
    if(mod(indate%hour,NHOURREAD)/=0.or.indate%seconds/=0)return
  endif

!never comes to this point if MODE=11 or 12

  if(MasterProc) print *,'NESTING'


  if(first_data==-1)then
     
    filename_read_BC=date2string(trim(filename_read_BC_template),date_nextfile) 
    filename_read_3D=date2string(trim(filename_read_3D),date_nextfile) !used only once
   !filename_read_BC file must be defined before reset_3D, because it is used by init_icbc

    if(.not.FORECAST) call reset_3D(ndays_indate)
    if(MasterProc) print *,'NEST: READING BC DATA from ',trim(filename_read_BC)
    call read_newdata_LATERAL(ndays_indate) 
    !the first hour only these values are used, no real interpolation between two records
  endif
  

  if(ndays_indate-rtime_saved(2)>halfsecond)then
    !look for a new data set
     filename_read_BC=date2string(trim(filename_read_BC_template),date_nextfile) !
    if(MasterProc) print *,'NEST: READING NEW BC DATA from ',trim(filename_read_BC)
    call read_newdata_LATERAL(ndays_indate)
  endif

!   make weights for time interpolation
  W1=1.0;  W2=0.0 ! default
  if(ndays_indate-rtime_saved(1)>halfsecond)then
    !interpolate
    W2=(ndays_indate-rtime_saved(1))/(rtime_saved(2)-rtime_saved(1))
    W1=1.0-W2
!   if(me==1)then
!     call nctime2idate(ndate,ndays_indate,'YYYY-MM-DD hh:mm:ss')
!     call nctime2idate(ndate,rtime_saved(1),'interpolating between YYYY-MM-DD hh:mm:ss')
!     call nctime2idate(ndate,rtime_saved(2),'and                   YYYY-MM-DD hh:mm:ss')
!     print *,'with weights : ',W1,W2
!   endif
  endif
 if(DEBUG_NEST.and.MasterProc) print *,'nesting BC 2D: time weights : ',W1,W2
 if(DEBUG_NEST.and.MasterProc) print *,'nesting BC 2D: time stamps : ',rtime_saved(1),rtime_saved(2)

  do n=1,NSPEC_ADV
     if(adv_bc(n)%wanted.and.adv_bc(n)%found)then
        if(DEBUG_ICBC.and.MasterProc) print *,'nesting component ',trim(adv_bc(n)%varname)
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
     endif
  enddo

if(RCA)then
!some components put to a fixed value
    if(kt==1)then
       !top
       xn_adv(IXADV_H2,:,:,kt)=       5e-7   
       xn_adv(IXADV_C2H4,:,:,kt)=     2e-10  
       xn_adv(IXADV_C3H6,:,:,kt)=     5e-11  
       xn_adv(IXADV_C2H5OH,:,:,kt)=   4e-10  
       xn_adv(IXADV_MEK,:,:,kt)=      2.5e-11
       xn_adv(IXADV_CH3O2H,:,:,kt)=   7.5e-11
       xn_adv(IXADV_MGLYOX,:,:,kt)=   0      
       xn_adv(IXADV_GLYOX,:,:,kt)=    0      
       xn_adv(IXADV_C2H5OOH,:,:,kt)=  1e-12  
   endif
    if(iw>=1)then
!west
       xn_adv(IXADV_H2,iw,:,:)=       5e-7   
       xn_adv(IXADV_C2H4,iw,:,:)=     2e-10  
       xn_adv(IXADV_C3H6,iw,:,:)=     5e-11  
       xn_adv(IXADV_C2H5OH,iw,:,:)=   4e-10  
       xn_adv(IXADV_MEK,iw,:,:)=      2.5e-11
       xn_adv(IXADV_CH3O2H,iw,:,:)=   1e-10  
       xn_adv(IXADV_MGLYOX,iw,:,:)=   2e-12  
       xn_adv(IXADV_GLYOX,iw,:,:)=    6e-12  
       xn_adv(IXADV_C2H5OOH,iw,:,:)=  1e-12  
    endif
    if(ie<=limax)then
!east
       xn_adv(IXADV_H2,ie,:,:)=        5e-7   
       xn_adv(IXADV_C2H4,ie,:,:)=      2e-10  
       xn_adv(IXADV_C3H6,ie,:,:)=      2e-10  
       xn_adv(IXADV_C2H5OH,ie,:,:)=    6e-10  
       xn_adv(IXADV_MEK,ie,:,:)=       5e-11  
       xn_adv(IXADV_CH3O2H,ie,:,:)=    1e-10  
       xn_adv(IXADV_MGLYOX,ie,:,:)=    1.5e-12
       xn_adv(IXADV_GLYOX,ie,:,:)=     1.3e-11
       xn_adv(IXADV_C2H5OOH,ie,:,:)=   1e-12  
    endif
    if(js>=1)then
!south
       xn_adv(IXADV_H2,:,js,:)=        5e-7
       xn_adv(IXADV_C2H4,:,js,:)=      5e-11
       xn_adv(IXADV_C3H6,:,js,:)=      1.6e-11
       xn_adv(IXADV_C2H5OH,:,js,:)=    7e-11
       xn_adv(IXADV_MEK,:,js,:)=       2.5e-11
       xn_adv(IXADV_CH3O2H,:,js,:)=    1e-10
       xn_adv(IXADV_MGLYOX,:,js,:)=    2e-12
       xn_adv(IXADV_GLYOX,:,js,:)=     4e-12
       xn_adv(IXADV_C2H5OOH,:,js,:)=   1e-12
    endif
    if(jn<=ljmax)then
!north
       xn_adv(IXADV_H2,:,jn,:)=        5e-7
       xn_adv(IXADV_C2H4,:,jn,:)=      2e-10
       xn_adv(IXADV_C3H6,:,jn,:)=      2e-10
       xn_adv(IXADV_C2H5OH,:,jn,:)=    4e-10
       xn_adv(IXADV_MEK,:,jn,:)=       2.5e-11
       xn_adv(IXADV_CH3O2H,:,jn,:)=    1e-12
       xn_adv(IXADV_MGLYOX,:,jn,:)=    2e-12
       xn_adv(IXADV_GLYOX,:,jn,:)=     4e-12
       xn_adv(IXADV_C2H5OOH,:,jn,:)=   1e-12
    endif
endif


  first_data=0
  first_call=.false.
  return
end subroutine readxn

subroutine wrtxn(indate,WriteNow)
  implicit none
  type(date), intent(in) :: indate
  logical, intent(in) :: WriteNow !Do not check indate value
  real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: dat ! Data arrays

  type(Deriv) :: def1 ! definition of fields
  integer :: n,iotyp,ndim,kmax
  real :: scale
  logical, save ::first_call=.true.

  if(MODE /= 1.and.MODE /= 3.and.MODE /= 10.and.MODE /= 12.and. .not.FORECAST)return

  if(FORECAST)then ! FORECAST mode superseeds nest MODE
    outdate(:)%seconds=0   ! output only at full hours
    if(.not.any(indate%year   ==outdate%year  .and.   &
                indate%month  ==outdate%month .and.   &
                indate%day    ==outdate%day   .and.   &
                indate%hour   ==outdate%hour  .and.   &
                indate%seconds==outdate%seconds))return
    if(MasterProc) print *,&
      date2string(" Forecast nest/dump at YYYY-MM-DD hh:mm:ss",indate)
    istart=RUNDOMAIN(1)
    jstart=RUNDOMAIN(3)
    iend=RUNDOMAIN(2)
    jend=RUNDOMAIN(4)
  elseif(MODE == 10.or.MODE == 12)then
    if(.not.WriteNow)return
    istart=RUNDOMAIN(1)
    jstart=RUNDOMAIN(3)
    iend=RUNDOMAIN(2)
    jend=RUNDOMAIN(4)
  else
    if(mod(indate%hour,NHOURSAVE)/=0.or.indate%seconds/=0)return
  endif

! fileName_write=date2string("EMEP_BC_MMYYYY.nc",indate)!for different names each month
                                                        !NB: readxn should have same name
  if(MasterProc)print *,'write Nest data ',trim(fileName_write)

  iotyp=IOU_INST
  if(first_call)then
    if(MasterProc)then
      print *,'Writing BC on ',trim(fileName_write)
     !write(command,*)'rm ',trim(fileName_write)
     !call system(command)
    endif
    first_call=.false.
  endif

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

  do n= 1, NSPEC_ADV
 !do n= 1, NSPEC_ADV-4  !ENEA
    def1%name= species(NSPEC_SHL+n)%name       !written
    dat=xn_adv(n,:,:,:)
    call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real4,&
        ist=istart,jst=jstart,ien=iend,jen=jend,fileName_given=fileName_write)
  enddo

  return
end subroutine wrtxn


subroutine check(status)
  use netcdf
  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
        WRITE(*,*) 'MPI_ABORT: ', "errorin NetCDF_ml"
        call  MPI_ABORT(MPI_COMM_WORLD,9,INFO)
  end if
end subroutine check

subroutine init_icbc()
  implicit none
  logical, save :: first_call=.true.
  integer :: n

  if(.not.first_call)return
  first_call=.false.

  if(all(adv_ic%varname==""))then
    adv_ic(:)%varname=species(NSPEC_SHL+1:NSPEC_SHL+NSPEC_ADV)%name
    adv_ic(:)%wanted=.true.
    adv_ic(:)%found=find_icbc(filename_read_3D,adv_ic%varname(:))
  endif
  if(all(adv_bc%varname==""))then
    if(FORECAST)then ! IFS-MOZART BC
      adv_bc(:)%wanted=.false.
      adv_bc(FORECAST_BC%ixadv)=FORECAST_BC%icbc
      adv_bc(:)%found=find_icbc(filename_read_bc,adv_bc%varname(:))
    elseif(TRANSPHORM)then ! TRANSPHORM BC
      adv_bc(:)%wanted=.false.
      adv_bc(TRANSPHORM_BC%ixadv)=TRANSPHORM_BC%icbc
      adv_bc(:)%found=find_icbc(filename_read_bc,adv_bc%varname(:))
    elseif(RCA)then ! RCA BC
      adv_bc(:)%wanted=.false.
      adv_bc(RCA_BC%ixadv)=RCA_BC%icbc
      adv_bc(:)%found=find_icbc(filename_read_bc,adv_bc%varname(:))
      adv_ic(:)=adv_bc(:)
    else
      adv_bc(:)=adv_ic(:)
    endif
  endif

  if((DEBUG_NEST.or.DEBUG_ICBC).and.MasterProc)then
    print "(A)","DEBUG_ICBC Variables:"
    print "(2(1X,A,I3,'=',A24,2L2))",&
      ('ADV_IC',n,adv_ic(n),'ADV_BC',n,adv_bc(n),n=1,NSPEC_ADV)
  endif
  contains
  function find_icbc(filename_read,varname) result(found)
    implicit none
    character(len=*), intent(in)               :: filename_read
    character(len=*), dimension(:), intent(in) :: varname
    logical, dimension(size(varname))          :: found
    integer :: status,ncFileID,varID,n

    found(:)=.false.
    if(MasterProc)then
      status = nf90_open(path=trim(filename_read),mode=nf90_nowrite,ncid=ncFileID)
      if(status /= nf90_noerr) then
        print *,'icbc: not found ',trim(filename_read)
      else
        print *,'icbc: reading ',trim(filename_read)
        do n=1,size(varname)
          if(varname(n)/="") &
            found(n)=(nf90_inq_varid(ncid=ncFileID,name=trim(varname(n)),varID=varID)==nf90_noerr)
        enddo
      endif
    endif
    CALL MPI_BCAST(found,size(found),MPI_LOGICAL,0,MPI_COMM_WORLD,INFO)
  end function find_icbc
end subroutine init_icbc

subroutine init_nest(ndays_indate,filename_read,IIij,JJij,Weight,&
                      k1_ext,k2_ext,weight_k1,weight_k2,&
                      N_ext,KMAX_ext,GIMAX_ext,GJMAX_ext)

  implicit none
  character(len=*),intent(in) :: filename_read
  real ,intent(out):: Weight(MAXLIMAX,MAXLJMAX,4)
  integer ,intent(out)::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)
  integer, intent(out), dimension(KMAX_MID) :: k1_ext,k2_ext
  real, intent(out), dimension(KMAX_MID) :: weight_k1,weight_k2
  integer ,intent(out)::N_ext,KMAX_ext,GIMAX_ext,GJMAX_ext
  real(kind=8) :: ndays_indate
  integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,status !,timeVarID
  integer :: ndate(4) !nseconds_indate,
  real :: dist(0:4),P_emep
  integer :: i,j,k,n,k_ext,II,JJ !nseconds(1),n,n1,k
  real, allocatable, dimension(:,:) ::lon_ext,lat_ext
  real, allocatable, dimension(:) ::hyam,hybm,P_ext
  character(len=80) ::projection,word
  logical :: reversed_k_BC,time_exists

  rtime_saved = -99999.9 !initialization

!Read dimensions (global)
  if(MasterProc)then
    status = nf90_open(path=trim(filename_read),mode=nf90_nowrite,ncid=ncFileID)
    if(status /= nf90_noerr) then
      print *,'init_nest: not found',trim(filename_read)
      return
    else
      print *,'init_nest: reading ',trim(filename_read)
    endif

    projection=''
    status = nf90_get_att(ncFileID,nf90_global,"projection",projection)
    if(status == nf90_noerr) then
       write(*,*)'projection: '
    else
       write(*,*)'projection not found for ',trim(filename_read)//', assuming lon lat'
       projection='lon lat'
    endif
    !get dimensions id
    if(trim(projection)=='Stereographic') then
      call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
      call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
    elseif(trim(projection)==trim('lon lat').or. &
          trim(projection)==trim('lon_lat')) then
      projection='lon lat'
      call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
      call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
    else
     !write(*,*)'GENERAL PROJECTION ',trim(projection)
      call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
      call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
     !WRITE(*,*) 'MPI_ABORT: ', "PROJECTION NOT RECOGNIZED"
     !call  MPI_ABORT(MPI_COMM_WORLD,9,INFO)
    endif

    status = nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID)
    if(status /= nf90_noerr) then
       status = nf90_inq_dimid(ncid = ncFileID, name = "mlev", dimID = kdimID)
       if(status /= nf90_noerr) then
          status = nf90_inq_dimid(ncid = ncFileID, name = "lev", dimID = kdimID)
          if(status /= nf90_noerr) then
             !include more possible names here
             write(*,*)'vertical levels name not found: ',trim(filename_read)
             call StopAll('Include new name in init_nest')
          endif
       endif
    endif


    N_ext=0
    status = nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID)
    if(status == nf90_noerr) then
       time_exists=.true.
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=N_ext))
    else
       write(*,*)'time dimension not found. Assuming only one record '
       time_exists=.false.
       N_ext=1
    endif

    call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_ext))
    call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_ext))
    call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_ext))

    write(*,*)'dimensions external grid',GIMAX_ext,GJMAX_ext,KMAX_ext,N_ext
    if(allocated(ndays_ext))then
       if(size(ndays_ext)<N_ext)then
          if(Masterproc)write(*,*)'Sizes times ',N_ext
          deallocate(ndays_ext)
          allocate(ndays_ext(N_ext))
       endif
    else
       allocate(ndays_ext(N_ext))
    endif

  endif
  CALL MPI_BCAST(GIMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(GJMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(KMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
!  CALL MPI_BCAST(N_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) !not needed by others than MaterProc

  allocate(lon_ext(GIMAX_ext,GJMAX_ext))
  allocate(lat_ext(GIMAX_ext,GJMAX_ext))
  allocate(hyam(KMAX_ext+1))
  allocate(hybm(KMAX_ext+1))
  allocate(P_ext(KMAX_ext))

  if(MasterProc)then
   !Read lon lat of the external grid (global)
    if(trim(projection)==trim('lon lat')) then
      call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
      call check(nf90_get_var(ncFileID, varID, lon_ext(:,1) ))
      do i=1,GJMAX_ext
        lon_ext(:,i)=lon_ext(:,1)
      enddo
      call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
      call check(nf90_get_var(ncFileID, varID, lat_ext(1,:) ))
      do i=1,GIMAX_ext
        lat_ext(i,:)=lat_ext(1,:)
      enddo
    else
      call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
      call check(nf90_get_var(ncFileID, varID, lon_ext ))

      call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
      call check(nf90_get_var(ncFileID, varID, lat_ext ))
    endif
!    call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))

   !do n=1,N_ext
!    call check(nf90_get_var(ncFileID, varID, ndays,start=(/ 1 /),count=(/ 1 /) ))
    if(time_exists)then
       call ReadTimeCDF(filename_read,ndays_ext,N_ext)
    else
       !cannot read time on file. assumes it is correct
       ndays_ext(1)=ndays_indate
    endif

    if(ndays_ext(1)-ndays_indate>halfsecond)then
      call nctime2idate(ndate,ndays_indate,&
        'WARNING: did not find BIC for date YYYY-MM-DD hh:mm:ss')
      call nctime2idate(ndate,ndays_ext(1),&
        'first date found YYYY-MM-DD hh:mm:ss')
    endif
    
    if(N_ext>1)then
       NHOURS_Stride_BC = nint((ndays_ext(2)-ndays_ext(1))*24)
    else
       !use manually set stride:
       NHOURS_Stride_BC = NHOURS_Stride_BC_default
    endif
    write(*,*)'new BC record every ',NHOURS_Stride_BC,' hours'

   !enddo
    !Read pressure for vertical levels
         write(*,*)'reading vertical level'

      status = nf90_inq_varid(ncid = ncFileID, name = "hyam", varID = varID)
      if(status == nf90_noerr) then
         call check(nf90_get_var(ncFileID, varID, hyam,count=(/ KMAX_ext /) ))
         call check(nf90_inq_varid(ncid = ncFileID, name = "hybm", varID = varID))
         call check(nf90_get_var(ncFileID, varID, hybm,count=(/ KMAX_ext /) ))
      else

      status = nf90_inq_varid(ncid = ncFileID, name = "k", varID = varID)
      if(status == nf90_noerr) then
         write(*,*)'assuming sigma level and PT=',PT,KMAX_ext
         call check(nf90_get_var(ncFileID, varID, hybm,count=(/ KMAX_ext /) ))!NB: here assume = sigma
         do k=1,KMAX_ext
            hyam(k)=PT*(1.0-hybm(k))
         enddo
      else

         !read eta levels from ad-hoc text file       
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
         !assumes lev=1000*(A+B) (IFS-MOZART?)
         !call check(nf90_inq_varid(ncid = ncFileID, name = "lev", varID = varID))
         !call check(nf90_get_var(ncFileID, varID, hybm ))
         !hybm=hybm/1000.0
         !hyam=0.0
      endif
      endif

    call check(nf90_close(ncFileID))
  endif !end MasterProc

  CALL MPI_BCAST(lon_ext,8*GIMAX_ext*GJMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(lat_ext,8*GIMAX_ext*GJMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(hyam,8*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(hybm,8*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

 !find horizontal interpolation constants
 !note that i,j are local
 !find the four closest points
  do j=1,ljmax
    do i=1,limax
      dist=1.0E40
      do JJ=1,GJMAX_ext
        do II=1,GIMAX_ext
         !distance between (i,j) and (II,JJ)
          dist(0)=great_circle_distance(lon_ext(II,JJ),lat_ext(II,JJ),glon(i,j),glat(i,j))
          if(dist(0)<dist(1))then
            dist(4)=dist(3)
            dist(3)=dist(2)
            dist(2)=dist(1)
            dist(1)=dist(0)
            IIij(i,j,4)=IIij(i,j,3)
            JJij(i,j,4)=JJij(i,j,3)
            IIij(i,j,3)=IIij(i,j,2)
            JJij(i,j,3)=JJij(i,j,2)
            IIij(i,j,2)=IIij(i,j,1)
            JJij(i,j,2)=JJij(i,j,1)
            IIij(i,j,1)=II
            JJij(i,j,1)=JJ
          elseif(dist(0)<dist(2))then
            dist(4)=dist(3)
            dist(3)=dist(2)
            dist(2)=dist(0)
            IIij(i,j,4)=IIij(i,j,3)
            JJij(i,j,4)=JJij(i,j,3)
            IIij(i,j,3)=IIij(i,j,2)
            JJij(i,j,3)=JJij(i,j,2)
            IIij(i,j,2)=II
            JJij(i,j,2)=JJ
          elseif(dist(0)<dist(3))then
            dist(4)=dist(3)
            dist(3)=dist(0)
            IIij(i,j,4)=IIij(i,j,3)
            JJij(i,j,4)=JJij(i,j,3)
            IIij(i,j,3)=II
            JJij(i,j,3)=JJ
          elseif(dist(0)<dist(4))then
            dist(4)=dist(0)
            IIij(i,j,4)=II
            JJij(i,j,4)=JJ
          endif
        enddo
      enddo

      dist(0)=(dist(1)+dist(2)+dist(3)+dist(4))
      Weight(i,j,1)=1.0-3.0*dist(1)/dist(0)
      dist(0)=(dist(2)+dist(3)+dist(4))
      Weight(i,j,2)=(1.0-Weight(i,j,1))*(1.0-2.0*dist(2)/dist(0))
      dist(0)=(dist(3)+dist(4))
      Weight(i,j,3)=(1.0-Weight(i,j,1)-Weight(i,j,2))*(1.0-dist(3)/dist(0))
      Weight(i,j,4)=1.0-Weight(i,j,1)-Weight(i,j,2)-Weight(i,j,3)
    enddo
  enddo

  deallocate(lon_ext,lat_ext)


  !find vertical interpolation coefficients
  !use pressure as reference
  !we want, if possible, P_ext(k1) and P_ext(k2) to be on each side of P_emep
  !We assume constant surface pressure, both for emep and external grid; should not be so 
  !   important as long as they are both terrain following.
  do k_ext=1,KMAX_EXT
     P_ext(k_ext)=hyam(k_ext)+hybm(k_ext)*Pref
     if(DEBUG_NEST.and.MasterProc) write(*,fmt="(A,I3,F10.2)")'P_ext',k_ext,P_ext(k_ext)
  enddo
  if(P_ext(1)>P_ext(2))then
  ! assumes that k_ext=KMAX_EXT is top and k_ext=1 is surface
     reversed_k_BC=.true.
  else
  ! assumes that k_ext=1 is top and k_ext=KMAX_EXT is surface
     reversed_k_BC=.false.
  endif

if(reversed_k_BC)then
  do k=1,KMAX_MID
     P_emep=A_mid(k)+B_mid(k)*Pref !Pa
     if(DEBUG_NEST.and.MasterProc) write(*,fmt="(A,I3,F10.2)")'P_emep',k,P_emep
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
     if(DEBUG_NEST.and.MasterProc)then
        write(*,fmt="(A,I3,A,I2,A,f4.2,A,I2,A,F4.2)")'level',k,' is the sum of level ',&
             k1_ext(k),' weight ',weight_k1(k),' and level ',k2_ext(k),' weight ',weight_k2(k)
     endif
  enddo

else
  do k=1,KMAX_MID
     P_emep=A_mid(k)+B_mid(k)*Pref !Pa
     if(DEBUG_NEST.and.MasterProc) write(*,fmt="(A,I3,F10.2)")'P_emep',k,P_emep
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
     if(DEBUG_NEST.and.MasterProc)then
        write(*,fmt="(A,I3,A,I2,A,f4.2,A,I2,A,F4.2)")'level',k,' is the sum of level ',&
             k1_ext(k),' weight ',weight_k1(k),' and level ',k2_ext(k),' weight ',weight_k2(k)
     endif
  enddo
endif
  deallocate(P_ext,hyam,hybm)

  if(DEBUG_NEST.and.MasterProc)write(*,*)'Nesting: finished determination of interpolation parameters'

end subroutine init_nest

subroutine read_newdata_LATERAL(ndays_indate)
  implicit none
  real(kind=8), intent(in)::ndays_indate
  real, allocatable, dimension(:,:,:) ::data
  integer :: ncFileID,varid,status
  integer :: ndate(4),n,i,j,k
  real(kind=8) :: ndays(1),ndays_old
  logical, save :: first_call=.true.

  !4 nearest points from external grid  (horizontal)
  integer, save ::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)
  !weights of the 4 nearest points (horizontal)
  real, save :: Weight(MAXLIMAX,MAXLJMAX,4)

  !2 adjacent levels from external grid  (vertical)
  integer, save, dimension(KMAX_MID) :: k1_ext,k2_ext
  !weights of the 2 adjacent levels (vertical)
  real, save, dimension(KMAX_MID) :: weight_k1,weight_k2

  integer:: KMAX_BC!which lvels are interpolated, = KMAX_MID for now
  integer:: timedimID

  !dimensions of external grid for BC
  integer, save ::GIMAX_ext,GJMAX_ext
  character (len=80) ::units
  real :: scale_factor,add_offset
  logical :: time_exists

  KMAX_BC=KMAX_MID
  if(first_call)then
     if(DEBUG_NEST.and.MasterProc)write(*,*)'Nesting: initializations 2D'
    call init_icbc()
    call init_nest(ndays_indate,filename_read_BC,IIij,JJij,Weight,&
                   k1_ext,k2_ext,weight_k1,weight_k2,&
                   N_ext_BC,KMAX_ext_BC,GIMAX_ext,GJMAX_ext)


    !Define & allocate West/East/South/Nort Boundaries
    iw=li0-1;ie=li1+1   ! i West/East   boundaries
    js=lj0-1;jn=lj1+1   ! j South/North boundaries
    kt=0;if(TOP_BC)kt=1 ! k Top         boundary
    if(DEBUG_NEST.and.MasterProc)then
       if(kt==1)then
          write(*,*)'Also including the top layer in BC'
       else
          write(*,*)'Not resetting the top layer'
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
      if(MasterProc) print "(A)","DEBUG_ICBC Boundaries:"
      print "(1X,'me=',i3,5(1X,A,I0,'=',L1))",&
        me,'W:i',iw,allocated(xn_adv_bndw),'E:i',ie,allocated(xn_adv_bnde),&
           'S:j',js,allocated(xn_adv_bnds),'N:j',jn,allocated(xn_adv_bndn),&
           'T:k',kt,allocated(xn_adv_bndt)
      if(MasterProc)flush(6)
      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    endif
    rtime_saved(2)=-99.0!just to put a value
     if(DEBUG_NEST.and.MasterProc)write(*,*)'Nesting: end initializations 2D'
  endif

  rtime_saved(1)=rtime_saved(2)!put old values in 1
  allocate(data(GIMAX_ext,GJMAX_ext,KMAX_ext_BC), stat=status)
  if(MasterProc)then
    call check(nf90_open(path = trim(fileName_read_BC), mode = nf90_nowrite, ncid = ncFileID))
    status = nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID)
    if(status == nf90_noerr) then
       time_exists=.true.
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=N_ext_BC))
    else
       time_exists=.false.
       N_ext_BC=1
    endif

     if(size(ndays_ext)<N_ext_BC)then
        write(*,*)'New size times in BC file ',N_ext_BC
        deallocate(ndays_ext)
        allocate(ndays_ext(N_ext_BC))
     endif

    if(time_exists)then
       call ReadTimeCDF(filename_read_BC,ndays_ext,N_ext_BC)
       !    call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))
       do n=1,N_ext_BC
          !      call check(nf90_get_var(ncFileID, varID, ndays,start=(/ n /),count=(/ 1 /) ))
          if(ndays_indate-ndays_ext(n)<halfsecond) goto 876
       enddo
       n=N_ext_BC
       write(*,*)'WARNING: did not find correct date ',n
876    continue
    else
       !cannot read time on file. assume it is corresponds to date_nextfile
       n=1
       call idate2nctime(date_nextfile,ndays_ext(n))
    endif
    call nctime2idate(ndate,ndays_ext(n),'Reading date YYYY-MM-DD hh:mm:ss')
    if(DEBUG_NEST.and.MasterProc)write(*,*)'Record ',n,' of ',N_ext_BC
    itime=n
    rtime_saved(2)=ndays_ext(n)
    if(n==N_ext_BC)then
       !next data to be read should be from another file
       if(DEBUG_NEST.and.MasterProc)then
          write(*,*)'Last record reached ',n,N_ext_BC
          call nctime2idate(date_nextfile,ndays_ext(n)+NHOURS_Stride_BC/24.0,'next BC date to read:  YYYY-MM-DD hh:mm:ss')
          write(*,*)'date_nextfile ',date_nextfile
       else
          call nctime2idate(date_nextfile,ndays_ext(n)+NHOURS_Stride_BC/24.0)
       endif
    endif

  endif

  CALL MPI_BCAST(rtime_saved,8*2,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

  DO_SPEC_init: do n= 1, NSPEC_ADV  !move BC used to nr=1 and set nr=2 to zero
    if(.not.(adv_bc(n)%wanted.and.adv_bc(n)%found)) cycle DO_SPEC_init
    if(.not.first_call)then
      !store the old values in 1
      if(iw>=1)     forall (k=1:KMAX_BC, j=1:ljmax) &
        xn_adv_bndw(n,j,k,1)=xn_adv_bndw(n,j,k,2)
      if(ie<=limax) forall (k=1:KMAX_BC, j=1:ljmax) &
        xn_adv_bnde(n,j,k,1)=xn_adv_bnde(n,j,k,2)
      if(js>=1)     forall (k=1:KMAX_BC, i=1:limax) &
        xn_adv_bnds(n,i,k,1)=xn_adv_bnds(n,i,k,2)
      if(jn<=ljmax) forall (k=1:KMAX_BC, i=1:limax) &
        xn_adv_bndn(n,i,k,1)=xn_adv_bndn(n,i,k,2)
      if(kt>=1)     forall (i=1:limax, j=1:ljmax) &
        xn_adv_bndt(n,i,j,1)=xn_adv_bndt(n,i,j,2)
    endif
      if(iw>=1)     forall (k=1:KMAX_BC, j=1:ljmax) &
        xn_adv_bndw(n,j,k,2)=0.0
      if(ie<=limax) forall (k=1:KMAX_BC, j=1:ljmax) &
        xn_adv_bnde(n,j,k,2)=0.0
      if(js>=1)     forall (k=1:KMAX_BC, i=1:limax) &
        xn_adv_bnds(n,i,k,2)=0.0
      if(jn<=ljmax) forall (k=1:KMAX_BC, i=1:limax) &
        xn_adv_bndn(n,i,k,2)=0.0
      if(kt>=1)     forall (i=1:limax, j=1:ljmax) &
        xn_adv_bndt(n,i,j,2)=0.0

  enddo DO_SPEC_init

  DO_SPEC: do n= 1, NSPEC_ADV
    if(.not.(adv_bc(n)%wanted.and.adv_bc(n)%found)) cycle DO_SPEC

    if(MasterProc)then
    !Could fetch one level at a time if sizes becomes too big
      call check(nf90_inq_varid(ncid=ncFileID, name=trim(adv_bc(n)%varname), varID=varID))

      call check(nf90_get_var(ncFileID, varID, data &
            ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext_BC,1 /) ))
      status = nf90_get_att(ncFileID,VarID,"scale_factor",scale_factor)
      if(status == nf90_noerr) then
         data=data*scale_factor
      endif
      status = nf90_get_att(ncFileID,VarID,"add_offset",add_offset)
      if(status == nf90_noerr) then
         data=data+add_offset
      endif
      status = nf90_get_att(ncFileID,VarID,"units",units)
      if(status == nf90_noerr) then
         if(DEBUG_NEST)write(*,*)'variable '//trim(adv_bc(n)%varname)//' has unit '//trim(units)
         if(units(1:3) == 'ppb') then
            if(DEBUG_NEST)write(*,*)'which is ppb unit. Scaling by ',PPB
            data=data*PPB
         else
            if(DEBUG_NEST)write(*,*)'which is not recognized as ppb unit. Assuming mixing ratio'
         endif
      else
        if(DEBUG_NEST)write(*,*)'units attribute not found for variable '//trim(adv_bc(n)%varname)
      endif
    endif
    CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext_BC,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

   !overwrite Global Boundaries (lateral faces)
    if(iw>=1)     forall (k=1:KMAX_BC, j=1:ljmax) &
      xn_adv_bndw(n,j,k,2)=xn_adv_bndw(n,j,k,2)+(Weight(iw,j,1)*data(IIij(iw,j,1),JJij(iw,j,1),k1_ext(k)) &
                           +Weight(iw,j,2)*data(IIij(iw,j,2),JJij(iw,j,2),k1_ext(k)) &
                           +Weight(iw,j,3)*data(IIij(iw,j,3),JJij(iw,j,3),k1_ext(k)) &
                           +Weight(iw,j,4)*data(IIij(iw,j,4),JJij(iw,j,4),k1_ext(k)))*weight_k1(k)&
                          +(Weight(iw,j,1)*data(IIij(iw,j,1),JJij(iw,j,1),k2_ext(k)) &
                           +Weight(iw,j,2)*data(IIij(iw,j,2),JJij(iw,j,2),k2_ext(k)) &
                           +Weight(iw,j,3)*data(IIij(iw,j,3),JJij(iw,j,3),k2_ext(k)) &
                           +Weight(iw,j,4)*data(IIij(iw,j,4),JJij(iw,j,4),k2_ext(k)))*weight_k2(k) 
    if(ie<=limax) forall (k=1:KMAX_BC, j=1:ljmax) &
      xn_adv_bnde(n,j,k,2)=xn_adv_bnde(n,j,k,2)+(Weight(ie,j,1)*data(IIij(ie,j,1),JJij(ie,j,1),k1_ext(k)) &             
                           +Weight(ie,j,2)*data(IIij(ie,j,2),JJij(ie,j,2),k1_ext(k)) &             
                           +Weight(ie,j,3)*data(IIij(ie,j,3),JJij(ie,j,3),k1_ext(k)) &             
                           +Weight(ie,j,4)*data(IIij(ie,j,4),JJij(ie,j,4),k1_ext(k)))*weight_k1(k)&
                          +(Weight(ie,j,1)*data(IIij(ie,j,1),JJij(ie,j,1),k2_ext(k)) &             
                           +Weight(ie,j,2)*data(IIij(ie,j,2),JJij(ie,j,2),k2_ext(k)) &             
                           +Weight(ie,j,3)*data(IIij(ie,j,3),JJij(ie,j,3),k2_ext(k)) &             
                           +Weight(ie,j,4)*data(IIij(ie,j,4),JJij(ie,j,4),k2_ext(k)))*weight_k2(k) 
    if(js>=1)     forall (k=1:KMAX_BC, i=1:limax) &
      xn_adv_bnds(n,i,k,2)=xn_adv_bnds(n,i,k,2)+(Weight(i,js,1)*data(IIij(i,js,1),JJij(i,js,1),k1_ext(k)) &             
                           +Weight(i,js,2)*data(IIij(i,js,2),JJij(i,js,2),k1_ext(k)) &             
                           +Weight(i,js,3)*data(IIij(i,js,3),JJij(i,js,3),k1_ext(k)) &             
                           +Weight(i,js,4)*data(IIij(i,js,4),JJij(i,js,4),k1_ext(k)))*weight_k1(k)&
                          +(Weight(i,js,1)*data(IIij(i,js,1),JJij(i,js,1),k2_ext(k)) &             
                           +Weight(i,js,2)*data(IIij(i,js,2),JJij(i,js,2),k2_ext(k)) &             
                           +Weight(i,js,3)*data(IIij(i,js,3),JJij(i,js,3),k2_ext(k)) &             
                           +Weight(i,js,4)*data(IIij(i,js,4),JJij(i,js,4),k2_ext(k)))*weight_k2(k) 
    if(jn<=ljmax) forall (k=1:KMAX_BC, i=1:limax) &
      xn_adv_bndn(n,i,k,2)=xn_adv_bndn(n,i,k,2)+(Weight(i,jn,1)*data(IIij(i,jn,1),JJij(i,jn,1),k1_ext(k)) &             
                           +Weight(i,jn,2)*data(IIij(i,jn,2),JJij(i,jn,2),k1_ext(k)) &             
                           +Weight(i,jn,3)*data(IIij(i,jn,3),JJij(i,jn,3),k1_ext(k)) &             
                           +Weight(i,jn,4)*data(IIij(i,jn,4),JJij(i,jn,4),k1_ext(k)))*weight_k1(k)&
                          +(Weight(i,jn,1)*data(IIij(i,jn,1),JJij(i,jn,1),k2_ext(k)) &             
                           +Weight(i,jn,2)*data(IIij(i,jn,2),JJij(i,jn,2),k2_ext(k)) &             
                           +Weight(i,jn,3)*data(IIij(i,jn,3),JJij(i,jn,3),k2_ext(k)) &             
                           +Weight(i,jn,4)*data(IIij(i,jn,4),JJij(i,jn,4),k2_ext(k)))*weight_k2(k) 
    if(kt>=1)     forall (i=1:limax, j=1:ljmax) &
      xn_adv_bndt(n,i,j,2)=xn_adv_bndt(n,i,j,2)+(Weight(i,j,1)*data(IIij(i,j,1),JJij(i,j,1),k1_ext(kt)) &            
                           +Weight(i,j,2)*data(IIij(i,j,2),JJij(i,j,2),k1_ext(kt)) &            
                           +Weight(i,j,3)*data(IIij(i,j,3),JJij(i,j,3),k1_ext(kt)) &            
                           +Weight(i,j,4)*data(IIij(i,j,4),JJij(i,j,4),k1_ext(kt)))*weight_k1(kt)&
                          +(Weight(i,j,1)*data(IIij(i,j,1),JJij(i,j,1),k2_ext(kt)) &            
                           +Weight(i,j,2)*data(IIij(i,j,2),JJij(i,j,2),k2_ext(kt)) &            
                           +Weight(i,j,3)*data(IIij(i,j,3),JJij(i,j,3),k2_ext(kt)) &            
                           +Weight(i,j,4)*data(IIij(i,j,4),JJij(i,j,4),k2_ext(kt)))*weight_k2(kt)
  enddo DO_SPEC

  if(first_call)then
     !copy 2 into 1 so that both are well defined
     rtime_saved(1)=rtime_saved(2)!put  time in 1
     DO_SPEC_init1: do n= 1, NSPEC_ADV  !copy BC used to nr=1
        if(.not.(adv_bc(n)%wanted.and.adv_bc(n)%found)) cycle DO_SPEC_init1
        !store the old values in 1
        if(iw>=1)     forall (k=1:KMAX_BC, j=1:ljmax) &
             xn_adv_bndw(n,j,k,1)=xn_adv_bndw(n,j,k,2)
        if(ie<=limax) forall (k=1:KMAX_BC, j=1:ljmax) &
             xn_adv_bnde(n,j,k,1)=xn_adv_bnde(n,j,k,2)
        if(js>=1)     forall (k=1:KMAX_BC, i=1:limax) &
             xn_adv_bnds(n,i,k,1)=xn_adv_bnds(n,i,k,2)
        if(jn<=ljmax) forall (k=1:KMAX_BC, i=1:limax) &
             xn_adv_bndn(n,i,k,1)=xn_adv_bndn(n,i,k,2)
        if(kt>=1)     forall (i=1:limax, j=1:ljmax) &
             xn_adv_bndt(n,i,j,1)=xn_adv_bndt(n,i,j,2)
     enddo DO_SPEC_init1
  endif

  deallocate(data)
  if(MasterProc) call check(nf90_close(ncFileID))
  first_call=.false.
  return
end subroutine read_newdata_LATERAL

subroutine reset_3D(ndays_indate)
  implicit none
  real(kind=8), intent(in)::ndays_indate
  real, allocatable, dimension(:,:,:) ::data
  integer :: ndate(4),n,i,j,k,itime=0,status
  integer :: ncFileID,varid
  real(kind=8) :: ndays(1)
  logical, save :: first_call=.true.

  !4 nearest points from external grid
  integer, save ::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)

  !weights of the 4 nearest points
  real, save :: Weight(MAXLIMAX,MAXLJMAX,4)

  !dimensions of external grid for 3D
  integer, save ::N_ext,KMAX_ext,GIMAX_ext,GJMAX_ext

  !2 adjacent levels from external grid  (vertical)
  integer, save, dimension(KMAX_MID) :: k1_ext,k2_ext
  !weights of the 2 adjacent levels (vertical)
  real, save, dimension(KMAX_MID) :: weight_k1,weight_k2
  character (len=80) ::units
  real :: scale_factor,add_offset

  if(first_call)then
     if(DEBUG_NEST.and.MasterProc)write(*,*)'Nesting: initializations 3D'
     first_call=.false.
     call init_icbc()
     call init_nest(ndays_indate,filename_read_3D,IIij,JJij,Weight,&
                    k1_ext,k2_ext,weight_k1,weight_k2,&
                    N_ext,KMAX_ext,GIMAX_ext,GJMAX_ext)
     if(DEBUG_NEST.and.MasterProc)write(*,*)'Nesting: end initializations 3D'
  endif
  allocate(data(GIMAX_ext,GJMAX_ext,KMAX_ext), stat=status)
  if(MasterProc)then
    call check(nf90_open(path = trim(fileName_read_3D), mode = nf90_nowrite, ncid = ncFileID))

!    call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))
    do n=1,N_ext
      !call check(nf90_get_var(ncFileID, varID, ndays,start=(/ n /),count=(/ 1 /) ))
      if(ndays_ext(n)>=ndays_indate) goto 876
    enddo
    n=N_ext
    write(*,*)'WARNING: did not find correct date'
876 continue
    call nctime2idate(ndate,ndays_ext(n),'Using date YYYY-MM-DD hh:mm:ss')
    itime=n
  endif

  if(DEBUG_NEST.and.MasterProc)write(*,*)'Nesting: overwrite 3D'
  DO_SPEC: do n= 1, NSPEC_ADV
     if(.not.(adv_ic(n)%wanted.and.adv_ic(n)%found)) cycle DO_SPEC
     if(MasterProc)then
        !Could fetch one level at a time if sizes becomes too big
        call check(nf90_inq_varid(ncid=ncFileID, name=trim(adv_ic(n)%varname), varID=varID))

        call check(nf90_get_var(ncFileID, varID, data &
             ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext,1 /) ))
        status = nf90_get_att(ncFileID,VarID,"scale_factor",scale_factor)
        if(status == nf90_noerr) then
           data=data*scale_factor
        endif
        status = nf90_get_att(ncFileID,VarID,"add_offset",add_offset)
        if(status == nf90_noerr) then
           data=data+add_offset
        endif
        if(DEBUG_NEST) print *,'nesting 3D component ',trim(adv_ic(n)%varname)
        status = nf90_get_att(ncFileID,VarID,"units",units)
        if(status == nf90_noerr) then
           if(DEBUG_NEST)write(*,*)'variable '//trim(adv_ic(n)%varname)//' has unit '//trim(units)
           if(units(1:3) == 'ppb') then
              if(DEBUG_NEST)write(*,*)'which is ppb unit. Scaling by ',PPB
              data=data*PPB
           else
              if(DEBUG_NEST)write(*,*)'which is not recognized as ppb unit. Assuming mixing ratio'
           endif
        else
           if(DEBUG_NEST)write(*,*)'units attribute not found for variable '//trim(adv_bc(n)%varname)
        endif
     endif
     CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
     
     !overwrite everything 3D (init)
     forall (k=1:KMAX_ext, j=1:ljmax, i=1:limax) &
      xn_adv(n,i,j,k)=(Weight(i,j,1)*data(IIij(i,j,1),JJij(i,j,1),k1_ext(k)) &
                     +Weight(i,j,2)*data(IIij(i,j,2),JJij(i,j,2),k1_ext(k)) &
                     +Weight(i,j,3)*data(IIij(i,j,3),JJij(i,j,3),k1_ext(k)) &
                     +Weight(i,j,4)*data(IIij(i,j,4),JJij(i,j,4),k1_ext(k)))*weight_k1(k)&
                     +(Weight(i,j,1)*data(IIij(i,j,1),JJij(i,j,1),k2_ext(k)) &
                     +Weight(i,j,2)*data(IIij(i,j,2),JJij(i,j,2),k2_ext(k)) &
                     +Weight(i,j,3)*data(IIij(i,j,3),JJij(i,j,3),k2_ext(k)) &
                     +Weight(i,j,4)*data(IIij(i,j,4),JJij(i,j,4),k2_ext(k)))*weight_k2(k) 

  enddo DO_SPEC

  deallocate(data)
  if(MasterProc) call check(nf90_close(ncFileID))
end subroutine reset_3D




end module Nest_ml

