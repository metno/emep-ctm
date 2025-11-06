! <Pollen_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.6>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2025 met.no
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
module Pollen_mod
!-----------------------------------------------------------------------!
! Birch pollen emission calculation based on
! M. Sofiev et al. 2006, doi:10.1007/s00484-006-0027-x
!
! Pollen emission based upon meteorology paparameters, and heatsum.
! Pollen particles are assumed of 22 um diameter and 800 kg/m3 density.
!-----------------------------------------------------------------------!
use Pollen_const_mod
use PhysicalConstants_mod, only: AVOG
use Biogenics_mod,         only: EMIS_BioNat, EmisNat
use CheckStop_mod,         only: CheckStop,CheckNC
use ChemDims_mod,          only: NSPEC_SHL
use ChemSpecs_mod,         only: species_adv
use Chemfields_mod,        only: xn_adv    ! emep model concs.
use Config_module,         only: USES
use Debug_module,          only: DEBUG
use DerivedFields_mod,     only: f_2d,d_2d ! D2D houtly (debug) output
use GasParticleCoeffs_mod, only: DDdefs
use Functions_mod,         only: heaviside
use GridValues_mod ,       only: glon, glat, debug_proc, debug_li, debug_lj
use Io_mod,                only: IO_NML
use Io_RunLog_mod,         only: PrintLog
use Landuse_mod,           only: LandCover
use LocalFractions_mod,    only: lf_rcemis_nat
use LocalVariables_mod,    only: Grid
use MetFields_mod,         only: surface_precip, ws_10m ,rh2m,t2_nwp,&
                                foundws10_met,foundprecip,pr,u_ref,z_bnd,z_mid
use MicroMet_mod,          only: Wind_at_h
use Config_module,         only: KMAX_MID, DataDir, &
                                 METSTEP, MasterProc, IOU_INST, RUNDOMAIN, &
                                 dt=>dt_advec, GRID_NAME=>GRID,&
                                 out_DOMAIN=>NEST_out_DOMAIN,&
                                 MODE_READ=>NEST_MODE_READ,&
                                 NEST_template_read_3D,NEST_template_write,NEST_template_dump
use MPI_Groups_mod,        only: MPI_INTEGER,MPI_LOGICAL,MPI_COMM_CALC,&
                                MasterPE,IERROR
use NetCDF_mod,            only: ReadField_CDF,Out_netCDF,GetCDF_modelgrid,&
                                ReadTimeCDF,CDFtype=>Real4
use netcdf,                only: nf90_close
use Par_mod,               only: limax, ljmax, LIMAX, LJMAX, me
use PhysicalConstants_mod, only: PI
use Radiation_mod,         only: daylength,sunrise
use OwnDataTypes_mod,      only: Deriv,TXTLEN_FILE
use ZchemData_mod,         only: rcemis
use SmallUtils_mod,        only: find_index,key2str
use SubMet_mod,            only: Sub
use TimeDate_mod,          only: current_date,daynumber,date,day_of_year
use TimeDate_ExtraUtil_mod,only: date2string,compare_date,date2nctime
use MPI_Groups_mod,        only: MasterPE,IERROR,MPI_COMM_WORLD,MPI_COMM_CALC,&
                                 MPI_INTEGER,MPI_LOGICAL,MPI_DOUBLE_PRECISION,&
                                 MPI_IN_PLACE,MPI_MIN,MPI_MAX!,MPI_ALLREDUCE
                                ! openMPI has no explicit interface for MPI_ALLREDUCE
!-------------------------------------
implicit none
private
public:: pollen_flux,pollen_dump

!** 1) Public (saved) Variables from module:

real, public,save , allocatable, dimension(:,:,:)::&
  heatsum,      & ! heatsum, needs to be remembered for forecast
  AreaPOLL,     & ! emission of pollen
  R               ! pollen released so far

! pollen arrays indexing, order must match with POLLEN_GROUP:
! birch,olive,grass,mugwort1, mugwort2,mugwort3,mugwort4,mugwort5
integer, save, dimension(POLLEN_NUM) :: &
  inat=-1,iadv=-1,itot=-1

!-------------------------------------
! Variables read from NetCDF Files
! pollen_flux       NetCDF file     NetCDF var
!   pollen_frac     birch_frac_nc   birch_frac
!   birch_corr      birch_corr_nc   scale_factor (year specific version)
!   birch_corr      birch_data_nc   cross
!   pollen_h_c      birch_data_nc   h_c
!   pollen_frac     olive_data_nc   olive_frac
!   pollen_h_c      olive_data_nc   olive_th
!   pollen_frac     rweed_frac_nc   fraction
!   rweed_start     rweed_data_nc   start_th
!   rweed_corr      rweed_data_nc   scale
!   pollen_frac     grass_field_nc  grass_frac
!   grass_start     grass_time_nc   grass_start
!   grass_end       grass_time_nc   grass_end
!   grass_end       grass_time_nc   grass_length+grass_start if grass_end not found
!   pollen_frac     mugwort*_data_nc   mugwort*_frac    *=1,..,5
!   mugwort*_start  mugwort*_data_nc   mugwort*_start
!   mugwort*_end    mugwort*_data_nc   mugwort*_end
! Variables write/read from dump/restart Files
! pollen_dump/pollen_read   pollen ype (SPC)    NetCDF var
!   xn_adv(i,:,:,:)         BIRCH,OLIVE,GRASS,MUGWT   SPC
!   N_TOT(i)-R(:,:,i)       BIRCH,OLIVE,GRASS,MUGWT   SPC//'_rest'
!   heatsum(:,:,i)          BIRCH,OLIVE         SPC//'_heatsum'
character(len=TXTLEN_FILE), save :: &
  birch_frac_nc ='birch_frac.nc',         &
  birch_data_nc ='pollen_data.nc',        &
  birch_corr_nc ='birch_factor_YYYY.nc',  &
  olive_data_nc ='olive_YYYY.nc',         &
  alder_data_nc ='alder_data.nc',         &
  rweed_frac_nc ='ragweed_frac.nc',       &
  rweed_data_nc ='ragweed_YYYY.nc',       &
  grass_field_nc='grass_frac.nc',         &
  grass_time_nc ='grass_time.nc',         &
  grass_mode    ='linear',      & ! 'linear' (old) | 'gamma' (new)
  mugwort1_data_nc='mugwort1_data.nc',         &
  mugwort2_data_nc='mugwort2_data.nc',         &
  mugwort3_data_nc='mugwort3_data.nc',         &
  mugwort4_data_nc='mugwort4_data.nc',         &
  mugwort5_data_nc='mugwort5_data.nc',         &
  mugwort_mode    ='gamma',      & ! 'linear' (old) | 'gamma' (new)
  template_read ='POLLEN_IN.nc',& ! dump/restart input
  template_write='POLLEN_OUT.nc'  ! dump/restart output
logical,save:: &
  overwrite=.true. ! append to existing template_write if false

type(date), parameter :: &
  date_first_birch=date(-1,3,1,0,0),date_last_birch=date(-1,8,1,0,0),&      ! Mar-Jul
  date_first_olive=date(-1,1,1,0,0),date_last_olive=date(-1,7,1,0,0),&      ! Jan-Jun
  date_first_alder=date_first_olive,date_last_alder=date_last_olive,&       ! Jan-Jun
  ! ragweed uses HS_startday_rweed and EndCDThr_rweed
  date_first_rweed=date(-1,6,1,0,0),date_last_rweed=date(-1,11,1,0,0)       ! Jun-Oct
type(date), save :: & ! will be updated when grass_time.nc is read
  date_first_grass=date(-1,3,1,0,0),date_last_grass=date(-1,9,1,0,0),&      ! Mar-Aug
  date_first_mugwort1=date_first_rweed,date_last_mugwort1=date_last_rweed,& ! Jun-Oct
  date_first_mugwort2=date_first_rweed,date_last_mugwort2=date_last_rweed,&
  date_first_mugwort3=date_first_rweed,date_last_mugwort3=date_last_rweed,&
  date_first_mugwort4=date_first_rweed,date_last_mugwort4=date_last_rweed,&
  date_first_mugwort5=date_first_rweed,date_last_mugwort5=date_last_rweed

real, save, allocatable, dimension(:,:) :: &
  rweed_start,        & ! day when daylight length goes below 14.5h (StartCDThr)
  grass_start,grass_end, & ! Stard/End day of grass, read in
  mugwort1_start,mugwort1_end, & ! Start/End day mugwort
  mugwort2_start,mugwort2_end, & ! Start/End day mugwort
  mugwort3_start,mugwort3_end, & ! Start/End day mugwort
  mugwort4_start,mugwort4_end, & ! Start/End day mugwort
  mugwort5_start,mugwort5_end  ! Start/End day mugwort

logical, parameter :: DEBUG_NC=.false.

!-------------------------------------------------------------------------!
contains
!-------------------------------------------------------------------------!
subroutine Config_Pollen()
  integer :: ios,g,n
  logical, save :: first_call = .true.
  NAMELIST /Pollen_config/&
    birch_frac_nc,birch_data_nc,birch_corr_nc,&
    olive_data_nc,alder_data_nc,rweed_frac_nc,rweed_data_nc,&
    grass_field_nc,grass_time_nc,grass_mode,&
    mugwort1_data_nc,mugwort2_data_nc,mugwort3_data_nc,&
    mugwort4_data_nc,mugwort5_data_nc,mugwort_mode,&   
    template_read,template_write,overwrite

  if(.not.first_call)return
  first_call = .false.

  ! check consistency between Pollen_const_mod and species_adv definitions
  call pollen_check()

  rewind(IO_NML)
  read(IO_NML,NML=Pollen_config,iostat=ios)
  call CheckStop(ios,"NML=Pollen_config")
  if(DEBUG%POLLEN.and.MasterProc)then
    write(*,*) "NAMELIST IS "
    write(*,NML=Pollen_config)
  end if
  if(any([template_read == NEST_template_read_3D, &
          template_read == NEST_template_dump,    &
          template_write == NEST_template_write,  &
          template_write == NEST_template_dump])) &
    call CheckStop(MasterProc,&
      "Pollen template_read/write should different than NEST_template_read_3D/write/dump")

  ! expand DataDir, YYYY and GRID keywords
  birch_frac_nc  = key2str(birch_frac_nc ,'DataDir',DataDir)
  birch_data_nc  = key2str(birch_data_nc ,'DataDir',DataDir)
  birch_corr_nc  = key2str(birch_corr_nc ,'DataDir',DataDir)
  birch_corr_nc  = key2str(birch_corr_nc ,'YYYY'   ,current_date%year)
  olive_data_nc  = key2str(olive_data_nc ,'DataDir',DataDir)
  olive_data_nc  = key2str(olive_data_nc ,'YYYY'   ,current_date%year)
  alder_data_nc  = key2str(alder_data_nc ,'DataDir',DataDir)
  alder_data_nc  = key2str(alder_data_nc ,'GRID'   ,GRID_NAME)
  rweed_frac_nc  = key2str(rweed_frac_nc ,'DataDir',DataDir)
  rweed_frac_nc  = key2str(rweed_frac_nc ,'GRID'   ,GRID_NAME)
  rweed_data_nc  = key2str(rweed_data_nc ,'DataDir',DataDir)
  rweed_data_nc  = key2str(rweed_data_nc ,'GRID'   ,GRID_NAME)
  grass_field_nc = key2str(grass_field_nc,'DataDir',DataDir)
  grass_time_nc  = key2str(grass_time_nc ,'DataDir',DataDir)
  mugwort1_data_nc = key2str(mugwort1_data_nc,'DataDir',DataDir)
  mugwort1_data_nc = key2str(mugwort1_data_nc,'GRID'   ,GRID_NAME)
  mugwort2_data_nc = key2str(mugwort2_data_nc,'DataDir',DataDir)
  mugwort2_data_nc = key2str(mugwort2_data_nc,'GRID'   ,GRID_NAME)
  mugwort3_data_nc = key2str(mugwort3_data_nc,'DataDir',DataDir)
  mugwort3_data_nc = key2str(mugwort3_data_nc,'GRID'   ,GRID_NAME)
  mugwort4_data_nc = key2str(mugwort4_data_nc,'DataDir',DataDir)
  mugwort4_data_nc = key2str(mugwort4_data_nc,'GRID'   ,GRID_NAME)
  mugwort5_data_nc = key2str(mugwort5_data_nc,'DataDir',DataDir)
  mugwort5_data_nc = key2str(mugwort5_data_nc,'GRID'   ,GRID_NAME)
  template_read  = key2str(template_read ,'DataDir',DataDir)
  template_write = key2str(template_write,'DataDir',DataDir)

  do g=1,POLLEN_NUM
    inat(g) = find_index(POLLEN_GROUP(g),EMIS_BioNat(:))
    iadv(g) = find_index(POLLEN_GROUP(g),species_adv(:)%name)
    itot(g) = iadv(g)+NSPEC_SHL
    call CheckStop(inat(g)<0,"EMIS_BioNat misses: "//POLLEN_GROUP(g))
    call CheckStop(iadv(g)<0,"species_adv misses: "//POLLEN_GROUP(g))

    ! check gravitational setling params
    select case(POLLEN_GROUP(g))
      case(RWEED,MUGWORT1:MUGWORT5)
        n=find_index('POLL18',DDdefs(:)%name)
      case(BIRCH,ALDER)
        n=find_index('POLL22',DDdefs(:)%name)
      case(OLIVE)
        n=find_index('POLL28',DDdefs(:)%name)
      case(GRASS)
        n=find_index('POLL32',DDdefs(:)%name)
      case default
        call CheckStop("Not implemented "//POLLEN_GROUP(g))
      end select
    call CheckStop(n<0,"DDdefs misses: "//POLLEN_GROUP(g))
    call CheckStop(DDdefs(n)%umDpgV,D_POLL(g),"Wrong DDdefs%umDpgV: "//POLLEN_GROUP(g))
!   call CheckStop(DDdefs(n)%sigma ,0.01     ,"Wrong DDdefs%sigma: "//POLLEN_GROUP(g))
    call CheckStop(DDdefs(n)%rho_p ,POLL_DENS*1e-3,"Wrong DDdefs%rho_p: "//POLLEN_GROUP(g))
  end do
  if(MasterProc)write(*,"(A,10(' adv#',I3,'=',A,1X,es10.3,:))") &
    "Pollen: ",(iadv(g),POLLEN_GROUP(g),grain_wt(g),g=1,POLLEN_NUM)

  allocate(heatsum(LIMAX,LJMAX,POLLEN_NUM-6),&     ! Grass+mugwort* does not need heatsum
           R(LIMAX,LJMAX,POLLEN_NUM))
  heatsum(:,:,:) = 0.00
  R(:,:,:) = 0.0
end subroutine Config_Pollen
!-------------------------------------------------------------------------!
function checkdates(nday,spc,i,j) result(ok)
  integer, intent(in)           :: nday
  character(len=*), intent(in)  :: spc
  integer, intent(in), optional :: i,j ! update GRASS start/end from coords
  logical :: ok
  logical, save :: first_call = .true.
  integer :: g
  integer, parameter :: iPOLLEN=POLLEN_NUM+1
  integer, save :: day(0:1,iPOLLEN) ! day(0,:)=first,day(1,:)=last
  if(first_call)then
    day(0,iBIRCH)=day_of_year(current_date%year,date_first_birch%month,date_first_birch%day)
    day(0,iOLIVE)=day_of_year(current_date%year,date_first_olive%month,date_first_olive%day)
    day(0,iALDER)=day_of_year(current_date%year,date_first_alder%month,date_first_alder%day)
    day(0,iRWEED)=FLOOR(HS_startday_rweed)
    day(0,iGRASS)=day_of_year(current_date%year,date_first_grass%month,date_first_grass%day)&
                 -FLOOR(uncert_day_grass)
    day(0,iMUGW1)=day_of_year(current_date%year,date_first_mugwort1%month,date_first_mugwort1%day)&
                 -FLOOR(uncert_day_mugwort)
    day(0,iMUGW2)=day_of_year(current_date%year,date_first_mugwort2%month,date_first_mugwort2%day)&
                 -FLOOR(uncert_day_mugwort)
    day(0,iMUGW3)=day_of_year(current_date%year,date_first_mugwort3%month,date_first_mugwort3%day)&
                 -FLOOR(uncert_day_mugwort)
    day(0,iMUGW4)=day_of_year(current_date%year,date_first_mugwort4%month,date_first_mugwort4%day)&
                 -FLOOR(uncert_day_mugwort)
    day(0,iMUGW5)=day_of_year(current_date%year,date_first_mugwort5%month,date_first_mugwort5%day)&
                 -FLOOR(uncert_day_mugwort)
    day(0,iPOLLEN)=MINVAL(day(0,:POLLEN_NUM))

    day(1,iBIRCH)=day_of_year(current_date%year,date_last_birch%month,date_last_birch%day)
    day(1,iOLIVE)=day_of_year(current_date%year,date_last_olive%month,date_last_olive%day)
    day(1,iALDER)=day_of_year(current_date%year,date_last_alder%month,date_last_alder%day)
    day(1,iRWEED)=CEILING(EndCDThr_rweed)
    day(1,iGRASS)=day_of_year(current_date%year,date_last_grass%month,date_last_grass%day)&
                 +CEILING(uncert_day_grass)
    day(1,iMUGW1)=day_of_year(current_date%year,date_last_mugwort1%month,date_last_mugwort1%day)&
                 +CEILING(uncert_day_mugwort)
    day(1,iMUGW2)=day_of_year(current_date%year,date_last_mugwort2%month,date_last_mugwort2%day)&
                 +CEILING(uncert_day_mugwort)
    day(1,iMUGW3)=day_of_year(current_date%year,date_last_mugwort3%month,date_last_mugwort3%day)&
                 +CEILING(uncert_day_mugwort)
    day(1,iMUGW4)=day_of_year(current_date%year,date_last_mugwort4%month,date_last_mugwort4%day)&
                 +CEILING(uncert_day_mugwort)
    day(1,iMUGW5)=day_of_year(current_date%year,date_last_mugwort5%month,date_last_mugwort5%day)&
                 +CEILING(uncert_day_mugwort)
    day(1,iPOLLEN)=MAXVAL(day(1,:POLLEN_NUM))
    first_call = .false.
  end if
  if(present(i).and.present(j))then
    day(0,iRWEED)=MIN(day(0,iRWEED),FLOOR(rweed_start(i,j)-uncert_day_rweed))
    day(0,iGRASS)=MIN(day(0,iGRASS),FLOOR(grass_start(i,j)-uncert_day_grass))
    day(1,iGRASS)=MAX(day(1,iGRASS),CEILING(grass_end(i,j)+uncert_day_grass))
    day(0,iMUGW1)=MIN(day(0,iMUGW1),FLOOR(mugwort1_start(i,j)-uncert_day_mugwort))
    day(1,iMUGW1)=MAX(day(1,iMUGW1),CEILING(mugwort1_end(i,j)+uncert_day_mugwort))
    day(0,iMUGW2)=MIN(day(0,iMUGW2),FLOOR(mugwort2_start(i,j)-uncert_day_mugwort))
    day(1,iMUGW2)=MAX(day(1,iMUGW2),CEILING(mugwort2_end(i,j)+uncert_day_mugwort))
    day(0,iMUGW3)=MIN(day(0,iMUGW3),FLOOR(mugwort3_start(i,j)-uncert_day_mugwort))
    day(1,iMUGW3)=MAX(day(1,iMUGW3),CEILING(mugwort3_end(i,j)+uncert_day_mugwort))
    day(0,iMUGW4)=MIN(day(0,iMUGW4),FLOOR(mugwort4_start(i,j)-uncert_day_mugwort))
    day(1,iMUGW4)=MAX(day(1,iMUGW4),CEILING(mugwort4_end(i,j)+uncert_day_mugwort))
    day(0,iMUGW5)=MIN(day(0,iMUGW5),FLOOR(mugwort5_start(i,j)-uncert_day_mugwort))
    day(1,iMUGW5)=MAX(day(1,iMUGW5),CEILING(mugwort5_end(i,j)+uncert_day_mugwort)) 
 end if
  select case(spc)
    case("POLLEN","P","p","pollen");g=iPOLLEN
    case(BIRCH,"B","b","birch");g=iBIRCH
    case(OLIVE,"O","o","olive");g=iOLIVE
    case(ALDER,"A","a","alder");g=iALDER
    case(RWEED,"R","r","rweed");g=iRWEED
    case(GRASS,"G","g","grass");g=iGRASS
    case(MUGWORT1,"M1","m1","mugwort1");g=iMUGW1
    case(MUGWORT2,"M2","m2","mugwort2");g=iMUGW2
    case(MUGWORT3,"M3","m3","mugwort3");g=iMUGW3
    case(MUGWORT4,"M4","m4","mugwort4");g=iMUGW4
    case(MUGWORT5,"M5","m5","mugwort5");g=iMUGW5
    case default; call CheckStop("Unknown pollen type: "//trim(spc))
  end select
  ok=(day(0,g)<=nday).and.(nday<=day(1,g))
end function checkdates
!-------------------------------------------------------------------------!
subroutine pollen_flux(i,j,debug_flag)
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug_flag

  real, parameter :: &
    kgm2h(POLLEN_NUM)=grain_wt(:)*1e-6*3600, & ! EmisNat: grains/m2/s --> kg/m2/h
    Z10  = 10.0,              & ! 10m height
    UnDef=  0.0

  logical, save :: first_call=.true.
  logical :: debug_ij=.false.,found=.false.,pollen_out(POLLEN_NUM)=.false.

  integer :: ii, jj, nlu, ilu, lu, info, g, n
  real :: scale,heatsum_min,heatsum_max,rcpoll,n2m(POLLEN_NUM),u10,prec,relhum

  real, save, allocatable, dimension(:,:,:) :: &
    pollen_frac,   &  ! fraction of pollen (birch/olive/alder/rweed/grass/mugwort*), read in
    pollen_h_c,    &  ! temperature treshold (birch/olive/alder), read in
    pollen_dH         ! flowering period [degree seconds] (olive/alder)
  real, save, allocatable, dimension(:,:) :: &
    birch_corr,     & ! correction field for p.emission, read in
    rweed_corr,     & ! relative productivity (PollenTotal=rweed_corr*1.7e7)
    t2_day            ! daily temperature

  ! Read in the different fields
  if(first_call) then
    if(.not.checkdates(daynumber,"pollen"))return
    call Config_Pollen()
    call pollen_read()

    allocate(AreaPOLL(LIMAX,LJMAX,POLLEN_NUM),t2_day(LIMAX,LJMAX))
    AreaPOLL(:,:,:)=0.0
    t2_day(:,:)=0.0

    allocate(pollen_frac(LIMAX,LJMAX,POLLEN_NUM),pollen_h_c(LIMAX,LJMAX,iBIRCH:iALDER))
    allocate(birch_corr(LIMAX,LJMAX),pollen_dH(LIMAX,LJMAX,iOLIVE:IALDER))
    allocate(rweed_start(LIMAX,LJMAX),rweed_corr(LIMAX,LJMAX))
    allocate(grass_start(LIMAX,LJMAX),grass_end(LIMAX,LJMAX))
    allocate(mugwort1_start(LIMAX,LJMAX),mugwort1_end(LIMAX,LJMAX))
    allocate(mugwort2_start(LIMAX,LJMAX),mugwort2_end(LIMAX,LJMAX))
    allocate(mugwort3_start(LIMAX,LJMAX),mugwort3_end(LIMAX,LJMAX))
    allocate(mugwort4_start(LIMAX,LJMAX),mugwort4_end(LIMAX,LJMAX))
    allocate(mugwort5_start(LIMAX,LJMAX),mugwort5_end(LIMAX,LJMAX))
    call ReadField_CDF(birch_frac_nc,'birch_frac',pollen_frac(:,:,iBIRCH),1, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(birch_data_nc,'h_c',pollen_h_c(:,:,iBIRCH),1, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
! birch: cross_corr for specific year
    call ReadField_CDF(birch_corr_nc,'scale_factor',birch_corr,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
! birch: cross_corr default
    if(.not.found)&
      call ReadField_CDF(birch_data_nc,'cross',birch_corr,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
! olive
    call ReadField_CDF(olive_data_nc,'olive_frac',pollen_frac(:,:,iOLIVE),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(olive_data_nc,'olive_th',pollen_h_c(:,:,iOLIVE),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(olive_data_nc,'olive_len',pollen_dH(:,:,iOLIVE),1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
    if(.not.found) pollen_dH(:,:,iOLIVE) = dH_d_olive
! alder
    call ReadField_CDF(alder_data_nc,'emis_mask_alder',pollen_frac(:,:,iALDER),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(alder_data_nc,'start_heatsum_threshold',pollen_h_c(:,:,iALDER),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(alder_data_nc,'end_heatsum_threshold',pollen_dH(:,:,iALDER),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
! ragweed
   call ReadField_CDF(rweed_frac_nc,'fraction',pollen_frac(:,:,iRWEED),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(rweed_data_nc,'start_th',rweed_start,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(rweed_data_nc,'scale',rweed_corr,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
! grass
    call ReadField_CDF(grass_field_nc,'grass_frac',pollen_frac(:,:,iGRASS),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(grass_time_nc,'grass_start',grass_start,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(grass_time_nc,'grass_end',grass_end,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
    if(.not.found)then
      call ReadField_CDF(grass_time_nc,'grass_length',grass_end,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
      where(grass_end/=UnDef.and.grass_start/=UnDef)&
        grass_end=grass_end+grass_start
    end if
! mugwort1
    call ReadField_CDF(mugwort1_data_nc,'mugwort1_frac',pollen_frac(:,:,iMUGW1),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort1_data_nc,'mugwort1_start',mugwort1_start,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort1_data_nc,'mugwort1_end',mugwort1_end,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
! mugwort2
    call ReadField_CDF(mugwort2_data_nc,'mugwort2_frac',pollen_frac(:,:,iMUGW2),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort2_data_nc,'mugwort2_start',mugwort2_start,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort2_data_nc,'mugwort2_end',mugwort2_end,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
! mugwort3
    call ReadField_CDF(mugwort3_data_nc,'mugwort3_frac',pollen_frac(:,:,iMUGW3),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort3_data_nc,'mugwort3_start',mugwort3_start,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort3_data_nc,'mugwort3_end',mugwort3_end,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
! mugwort4
    call ReadField_CDF(mugwort4_data_nc,'mugwort4_frac',pollen_frac(:,:,iMUGW4),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort4_data_nc,'mugwort4_start',mugwort4_start,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort4_data_nc,'mugwort4_end',mugwort4_end,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
! mugwort5
    call ReadField_CDF(mugwort5_data_nc,'mugwort5_frac',pollen_frac(:,:,iMUGW5),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort5_data_nc,'mugwort5_start',mugwort5_start,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(mugwort5_data_nc,'mugwort5_end',mugwort5_end,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)


    ! negative fractions as UnDef
    where(pollen_frac/=UnDef .and. pollen_frac<0) &
      pollen_frac = UnDef

    ! reduce birch from 60 degrees north:
    forall(ii=1:limax,jj=1:ljmax,glat(ii,jj)>=60.0 .and. pollen_frac(ii,jj,iBIRCH)/=UnDef) &
      pollen_frac(ii,jj,iBIRCH)=pollen_frac(ii,jj,iBIRCH)&
                       *MAX(0.3,1.0-0.005*(glat(ii,jj)-60.0))

    ! olive fraction [%] --> [1/1]
    where(pollen_frac(:,:,iOLIVE)/=UnDef) &
      pollen_frac(:,:,iOLIVE) = pollen_frac(:,:,iOLIVE) / 100.0

    ! alder [dagree hour]
    where(pollen_dH(:,:,iALDER)/=UnDef)&
      pollen_dH(:,:,iALDER) = pollen_dH(:,:,iALDER) - pollen_h_c(:,:,iALDER)

    ! Set D2D/USET output
    call write_uset()

    ! start of pollen forecast
    if(checkdates(daynumber,"pollen",i=i,j=j).and.MasterProc)&
      write(*,*) "POLLEN setup ",date2string("YYYY-MM-DD hh:mm",current_date)
    first_call = .false.
  end if !first_call
  debug_ij=all([DEBUG%POLLEN.or.debug_flag,debug_proc,i==debug_li,j==debug_lj])

  do g=1,POLLEN_NUM
    select case(g)
    case(iBIRCH)
      pollen_out(g)=any([pollen_frac(i,j,g),pollen_h_c(i,j,g),birch_corr(i,j)]==UnDef)
    case(iOLIVE)
      pollen_out(g)=any([pollen_frac(i,j,g),pollen_h_c(i,j,g),pollen_dH(i,j,g)]==UnDef)
    case(iALDER)
      pollen_out(g)=any([pollen_frac(i,j,g),pollen_h_c(i,j,g),pollen_dH(i,j,g)]==UnDef)
    case(iRWEED)
      pollen_out(g)=any([pollen_frac(i,j,g),rweed_start(i,j),rweed_corr(i,j)]==UnDef)
    case(iGRASS)
      pollen_out(g)=any([pollen_frac(i,j,g),grass_start(i,j),grass_end(i,j)]==UnDef)
    case(iMUGW1)
      pollen_out(g)=any([pollen_frac(i,j,g),mugwort1_start(i,j),mugwort1_end(i,j)]==UnDef)
    case(iMUGW2)
      pollen_out(g)=any([pollen_frac(i,j,g),mugwort2_start(i,j),mugwort2_end(i,j)]==UnDef)
    case(iMUGW3)
      pollen_out(g)=any([pollen_frac(i,j,g),mugwort3_start(i,j),mugwort3_end(i,j)]==UnDef)
    case(iMUGW4)
      pollen_out(g)=any([pollen_frac(i,j,g),mugwort4_start(i,j),mugwort4_end(i,j)]==UnDef)
    case(iMUGW5)
      pollen_out(g)=any([pollen_frac(i,j,g),mugwort5_start(i,j),mugwort5_end(i,j)]==UnDef)
    case default
      call CheckStop("Not implemented "//POLLEN_GROUP(g))
    end select
    pollen_out(g)=pollen_out(g).or..not.checkdates(daynumber,POLLEN_GROUP(g))
  end do

  EmisNat(inat(:),i,j)      = 0.0
  rcemis (itot(:),KMAX_MID) = 0.0
  AreaPOLL(i,j,:)           = 0.0
  if(all(pollen_out))then
    call write_uset()
    return
  end if

  ! Heatsum: sums up the temperatures that day for all time-steps
  do g=1,POLLEN_NUM
    if(pollen_out(g)) cycle
    select case(g)
    case(iBIRCH,iOLIVE)
      call heatsum_day(heatsum(i,j,g),t2_nwp(i,j,1),T_cutoff(g))
    case(iALDER)
      call heatsum_hour(heatsum(i,j,g),t2_nwp(i,j,1),T_cutoff(g))
    case(iRWEED)
      call heatsum_rweed(heatsum(i,j,g),t2_nwp(i,j,1),daylength(glat(i,j)))
    case(iGRASS,iMUGW1:iMUGW5)
      ! nothing to do
    case default
      call CheckStop("Not implemented "//POLLEN_GROUP(g))
    end select
  end do

  ! calculate daily mean temperatures
  if(current_date%hour==0 .and. current_date%seconds==0)&
    t2_day(i,j) = 0.0
  t2_day(i,j) = t2_day(i,j) + t2_nwp(i,j,1)*dt/(3600*24)

  !Need to set the meteorological fields if not defined in nwp.
  relhum = 0.0
  u10    = 0.0
  prec   = 0.0

  if(foundws10_met) u10=ws_10m(i,j,1)
  nlu = LandCover(i,j)%ncodes
  do ilu=1,nlu
    lu = LandCover(i,j)%codes(ilu)
    relhum = Sub(lu)%rh
    !wind
    if(.not.foundws10_met .and. Sub(lu)%z0 > 1.0)&
      u10=Wind_at_h(Grid%u_ref, Grid%z_ref, Z10,&
      Sub(lu)%d, Sub(lu)%z0, Sub(lu)%invL)
    if(relhum>=0.0 .and. u10>=0.0) exit
  end do

  ! precipitation
  if(foundprecip) then
    prec=sum(pr(i,j,:))*60  !mm/h
  else
    prec=surface_precip(i,j)
  end if

  do g=1,POLLEN_NUM
    if(pollen_out(g)) cycle
    select case(g)
    case(iBIRCH)  ! Birch specific emission inhibitors
      heatsum_min = (1-prob_in(g))*pollen_h_c(i,j,g)
      heatsum_max = pollen_h_c(i,j,g) + dH_d_birch
      pollen_out(g)=(R(i,j,g)>N_TOT(g)*birch_corr(i,j)) & ! out of pollen
        .or.(relhum>RH_HIGH)                            & ! too humid
        .or.(prec>prec_max)                             & ! too rainy
        .or.(heatsum(i,j,g)<heatsum_min)                & ! too cold season
        .or.(t2_nwp(i,j,1)<T_cutoff(g))                 & ! too cold
        .or.(heatsum(i,j,g)>heatsum_max)                  ! too warm season
    case(iOLIVE)  ! Olive specific emission inhibitors
      heatsum_min = (1-prob_in(g))*pollen_h_c(i,j,g)
      heatsum_max = pollen_h_c(i,j,g) + pollen_dH(i,j,g)
      pollen_out(g)=(R(i,j,g)>N_TOT(g))                 & ! out of pollen
        .or.(relhum>RH_HIGH)                            & ! too humid
        .or.(prec>prec_max)                             & ! too rainy
        .or.(heatsum(i,j,g)<heatsum_min)                & ! too cold season
        .or.(t2_nwp(i,j,1)<T_cutoff(g))                 & ! too cold
        .or.(heatsum(i,j,g)>heatsum_max)                  ! too warm season
    case(iALDER)  ! Alder specific emission inhibitors
      heatsum_min = (1-prob_in(g))*pollen_h_c(i,j,g)
      heatsum_max = pollen_h_c(i,j,g) + pollen_dH(i,j,g)
      pollen_out(g)=(R(i,j,g)>N_TOT(g))                 & ! out of pollen
        .or.(relhum>RH_HIGH)                            & ! too humid
        .or.(prec>prec_max)                             & ! too rainy
        .or.(heatsum(i,j,g)<heatsum_min)                & ! too cold season
        .or.(t2_nwp(i,j,1)<T_cutoff(g))                 & ! too cold
        .or.(heatsum(i,j,g)>heatsum_max)                  ! too warm season
    case(iRWEED)  ! Ragweed specific emission inhibitors
      heatsum_min = (1.0-uncert_HS_rweed)*startThr_HS_rweed
      ! delay past calendar day threshold; change value to preserve distribution function
      if(heatsum(i,j,g)<heatsum_min .and. daynumber>rweed_start(i,j))&
        rweed_start(i,j) = daynumber
      ! meteo flowering thresholds can zero the heatsum
      if(t2_nwp(i,j,1)<TempThr_rweed .or. & ! inst. treshold
        (current_date%hour==23 .and. current_date%seconds==3600-int(dt) .and. &
        t2_day(i,j)<DayTempThr_rweed))then
        heatsum(i,j,g)=0.0
        if(R(i,j,g)>0.0) &                        ! if flowering started 
          R(i,j,g)=N_TOT(g)*rweed_corr(i,j)+1e-6  ! end emission in gridcell
      end if
      pollen_out(g)=(R(i,j,g)>N_TOT(g)*rweed_corr(i,j)) & ! out of pollen
        .or.(heatsum(i,j,g)<heatsum_min)                  ! too cold season
!     call CheckStop(.not.pollen_out(g),RWEED//': start season')
    case(iGRASS)! Grass specific emission inhibitors
      pollen_out(g)=(R(i,j,g)>N_TOT(g))                 & ! out of pollen
        .or.(relhum>RH_HIGH)                            & ! too humid
        .or.(prec>prec_max)                               ! too rainy
    case(iMUGW1:iMUGW5)! Mugwort specific emission inhibitors
      pollen_out(g)=(R(i,j,g)>N_TOT(g))                 & ! out of pollen
        .or.(relhum>RH_HIGH_MUGW)                       & ! too humid
        .or.(prec>prec_max)                               ! too rainy
    case default
      call CheckStop("Not implemented "//POLLEN_GROUP(g))
    end select
  end do

  ! grains to molecular weight: grains/m2/s --> mol/cm3/s
  n2m(:) = 1e-6/Grid%DeltaZ           & ! 1/dZ [1/cm3]
    *grain_wt*AVOG/species_adv(iadv(:))%molwt
!  write(*,*) "n2m_diff ", n2m(:)

!------------------------
! Emission rates: Birch,Olive,Alder,Ragweed,Grass
!------------------------
  do g=1,POLLEN_NUM
    if(pollen_out(g)) cycle
    ! scale factor meteorological conditions
    select case(g)
    case(iBIRCH,iOLIVE,iALDER,iRWEED)
      scale=scale_factor(POLLEN_GROUP(g))
    case(iGRASS);
      scale=scale_factor(GRASS//trim(grass_mode))
    case(iMUGW1);
      scale=scale_factor(MUGWORT1//trim(mugwort_mode))
    case(iMUGW2);
      scale=scale_factor(MUGWORT2//trim(mugwort_mode))
    case(iMUGW3);
      scale=scale_factor(MUGWORT3//trim(mugwort_mode))
    case(iMUGW4);
      scale=scale_factor(MUGWORT4//trim(mugwort_mode))
    case(iMUGW5);
      scale=scale_factor(MUGWORT5//trim(mugwort_mode))
    case default
      call CheckStop("Not implemented "//POLLEN_GROUP(g))
    end select
    R(i,j,g)=R(i,j,g)+N_TOT(g)*scale*dt       ! pollen grains released so far
    rcpoll=N_TOT(g)*scale*pollen_frac(i,j,g)  ! pollen grains production [grains/m2/s]
    EmisNat(inat(g),i,j)      = rcpoll*kgm2h(g) ! [kg/m2/h]
    rcemis (itot(g),KMAX_MID) = rcpoll*n2m(g) ! [mol/cm3/s]
    AreaPOLL(i,j,g)           = rcpoll*3600   ! [grains/m2/h]

    !NB: lf_rcemis_nat requires units of g/cm3/s
    if (USES%LocalFractions) call lf_rcemis_nat(itot(g), rcemis(itot(g),KMAX_MID)*species(itot(g))%molwt, i, j)
    
    if(debug_ij) write(*,'(a,3(1x,I3),3(1x,es10.3))')&
      POLLEN_GROUP(g),me,i,j,rcemis(itot(g),KMAX_MID),R(i,j,g),AreaPOLL(i,j,g)
  end do
  call write_uset()
contains
!------------------------
! scale factor for different meteorological conditions
!------------------------
function scale_factor(spc) result(scale)
  character(len=*), intent(in)  :: spc
  real :: scale, dH_secs, sec_since_sunrise, seasLen
  integer :: g

  scale=f_wind(u10,Grid%wstar)           & ! wind dependence
       *f_cond(relhum,RH_LOW,RH_HIGH)    & ! rh dependence
       *f_cond(prec,PREC_MIN,PREC_MAX)     ! precipitation dependence
  select case(spc)
  case(BIRCH)
    g=iBIRCH
    scale = scale*birch_corr(i,j) &
      *(t2_nwp(i,j,1)-T_cutoff(g))/dH_birch &
      *f_in(heatsum(i,j,g),pollen_h_c(i,j,g),PROB_IN(g)) &! prob. flowering start
      *f_out(R(i,j,g),N_TOT(g),PROB_OUT(g))               ! prob. flowering end
  case(OLIVE)
    g=iOLIVE
    dH_secs = pollen_dH(i,j,g)*86400.0       ! Flowering period [degree seconds]
    scale = scale &
      *(t2_nwp(i,j,1)-T_cutoff(g))/dH_secs &
      *f_in(heatsum(i,j,g),pollen_h_c(i,j,g),PROB_IN(g)) &! prob. flowering start
      *f_out(R(i,j,g),N_TOT(g),PROB_OUT(g))               ! prob. flowering end
  case(ALDER)
    g=iALDER
    dH_secs = pollen_dH(i,j,g)*3600.0        ! Flowering period [degree seconds]
    scale = scale &
      *(t2_nwp(i,j,1)-T_cutoff(g))/dH_secs &
      *f_in(heatsum(i,j,g),pollen_h_c(i,j,g),PROB_IN(g)) &! prob. flowering start
      *f_out(R(i,j,g),N_TOT(g),PROB_OUT(g))               ! prob. flowering end
  case(RWEED)
    g=iRWEED
    sec_since_sunrise=(current_date%hour-sunrise(glat(i,j),glon(i,j)))*3600 &
                      +current_date%seconds
    sec_since_sunrise=modulo(sec_since_sunrise,86400.0)
    scale = rweed_corr(i,j) &
      *regweed_normal_diurnal(daynumber,sec_since_sunrise,&
          rweed_start(i,j),EndCDThr_rweed)          ! StartCDThr,EndCDThr
    scale = scale/dt                                ! emited fration over dt to emission rate
  case(GRASS//'linear')   ! emission mass assuming linear release
    g=iGRASS
    scale = scale &
      *f_fade_in (real(daynumber)/grass_start(i,j),   &
                  uncert_day_grass/grass_start(i,j))  &   ! fade-in
      *f_fade_out(real(daynumber)/grass_end(i,j),     &
                  uncert_day_grass/grass_start(i,j))  &   ! fade-out
      *f_fade_out(R(i,j,g)/N_TOT(g),1.0-uncert_tot_grass) ! total-pollen fade-out
    ! Full-emission rate is total pollen divided by the total duration of the season
    scale = scale/(grass_end(i,j)-daynumber+uncert_day_grass)&
                 /(grass_end(i,j)-grass_start(i,j))/86400.0
  case(GRASS//'gamma')    ! assume the modified "taily" Gamma distribution of the season
    g=iGRASS
    seasLen = grass_end(i,j)-grass_start(i,j)       ! season length in days
    scale = f_gamma_w_tails(                      &
      (real(daynumber)-grass_start(i,j))/seasLen, & ! days since season start/season length
      (dt/86400.0)/seasLen)                         ! timestep in days/season length
    scale = scale/dt                                ! emited fration over dt to emission rate
  case(MUGWORT1//'gamma')    ! assume the modified "taily" Gamma distribution of the season
    g=iMUGW1
    seasLen = mugwort1_end(i,j)-mugwort1_start(i,j)       ! season length in days
    scale = f_gamma_w_tails(                      &
      (real(daynumber)-mugwort1_start(i,j))/seasLen, & ! days since season start/season length
      (dt/86400.0)/seasLen)                         ! timestep in days/season length
    scale = scale/dt                                ! emited fration over dt to emission rate
  case(MUGWORT2//'gamma')    ! assume the modified "taily" Gamma distribution of the season
    g=iMUGW2
    seasLen = mugwort2_end(i,j)-mugwort2_start(i,j)       ! season length in days
    scale = f_gamma_w_tails(                      &
      (real(daynumber)-mugwort2_start(i,j))/seasLen, & ! days since season start/season length
      (dt/86400.0)/seasLen)                         ! timestep in days/season length
    scale = scale/dt                                ! emited fration over dt to emission rate
  case(MUGWORT3//'gamma')    ! assume the modified "taily" Gamma distribution of the season
    g=iMUGW3
    seasLen = mugwort3_end(i,j)-mugwort3_start(i,j)       ! season length in days
    scale = f_gamma_w_tails(                      &
      (real(daynumber)-mugwort3_start(i,j))/seasLen, & ! days since season start/season length
      (dt/86400.0)/seasLen)                         ! timestep in days/season length
    scale = scale/dt                                ! emited fration over dt to emission rate
  case(MUGWORT4//'gamma')    ! assume the modified "taily" Gamma distribution of the season
    g=iMUGW4
    seasLen = mugwort4_end(i,j)-mugwort4_start(i,j)       ! season length in days
    scale = f_gamma_w_tails(                      &
      (real(daynumber)-mugwort4_start(i,j))/seasLen, & ! days since season start/season length
      (dt/86400.0)/seasLen)                         ! timestep in days/season length
    scale = scale/dt                                ! emited fration over dt to emission rate
  case(MUGWORT5//'gamma')    ! assume the modified "taily" Gamma distribution of the season
    g=iMUGW5
    seasLen = mugwort5_end(i,j)-mugwort5_start(i,j)       ! season length in days
    scale = f_gamma_w_tails(                      &
      (real(daynumber)-mugwort5_start(i,j))/seasLen, & ! days since season start/season length
      (dt/86400.0)/seasLen)                         ! timestep in days/season length
    scale = scale/dt                                ! emited fration over dt to emission rate
  case default
    call CheckStop("Unknown pollen type: "//trim(spc))
  end select
end function scale_factor
!------------------------
! Write D2D/USET output at the end of every hour
!------------------------
subroutine write_uset()
  ! indexes for USET/D2D debug output
  integer, save, dimension(POLLEN_NUM) :: &
    n2d_heatsum=-1,n2d_left=-1,n2d_emiss=-1
  integer :: g,n
  if(first_call)then
    do n=1,size(f_2d)
      if(f_2d(n)%class/='USET')cycle
      g=find_index(f_2d(n)%txt,POLLEN_GROUP)
      if(g<0)cycle
      select case(f_2d(n)%subclass)
      case("heatsum")
        call CheckStop(g,[1,size(heatsum,DIM=3)],&
          "USET: '"//trim(f_2d(n)%subclass)//"' out of bounds!")
        f_2d(n)%unit="degreedays"
        if(g==iALDER)&
          f_2d(n)%unit="degreehours"
        f_2d(n)%scale=1.0
        n2d_heatsum(g)=n
      case("pollen_left")
        call CheckStop(g,[1,size(R,DIM=3)],&
          "USET: '"//trim(f_2d(n)%subclass)//"' out of bounds!")
        f_2d(n)%unit="1"
        f_2d(n)%scale=1.0
        n2d_left(g)=n
      case("pollen_emiss")
        call CheckStop(g,[1,size(AreaPOLL,DIM=3)],&
          "USET: '"//trim(f_2d(n)%subclass)//"' out of bounds!")
        f_2d(n)%unit="grains/m2/h"
        f_2d(n)%scale=1.0
        n2d_emiss(g)=n
      case default
        cycle
      end select
    end do
    return
  end if
  do g=1,POLLEN_NUM
    n=n2d_heatsum(g)  ! heatsum
    if(n>0) d_2d(n,i,j,IOU_INST)=heatsum(i,j,g)
    n=n2d_left(g)     ! pollen_left
    if(n>0) d_2d(n,i,j,IOU_INST)=1.0-R(i,j,g)/N_TOT(g)
    n=n2d_emiss(g)    ! pollen_emiss
    if(n>0) d_2d(n,i,j,IOU_INST)=AreaPOLL(i,j,g)
  end do
end subroutine write_uset
end subroutine pollen_flux
subroutine heatsum_day(hsum,t2,T_cutoff)
! The temperature needs to be over the cutoff temperature
  real, intent(inout) :: hsum
  real, intent(in) :: t2,T_cutoff ! degreedays
  real             :: ff
  ff = DIM(t2,T_cutoff) ! same as MAX(t2-T_cutoff,0.0)
  hsum = hsum + ff*dt/(3600*24)      ! seconds to days
end subroutine heatsum_day
subroutine heatsum_hour(hsum,t2,T_cutoff)
  ! The temperature needs to be over the cutoff temperature
    real, intent(inout) :: hsum
    real, intent(in) :: t2,T_cutoff ! degreehours
    real             :: ff
    ff = DIM(t2,T_cutoff) ! same as MAX(t2-T_cutoff,0.0)
    hsum = hsum + ff*dt/3600         ! seconds to hours
  end subroutine heatsum_hour
  subroutine heatsum_rweed(hsum,t2,daylen)
  real, intent(inout) :: hsum    ! degreedays
  real, intent(in) :: t2,daylen  ! deg K,date%hours
  real             :: ff

  ! Too early to do anything
  if(daynumber<HS_startday_rweed) return

  ! Temperature response for biotime accumulation
  if (t2>=temp_min_rweed .and. t2<=temp_opt_rweed)then
    ff = (t2 - temp_min_rweed) / (temp_opt_rweed - temp_min_rweed)
  elseif (t2>temp_opt_rweed .and. t2<=temp_max_rweed)then
    ff = (temp_max_rweed - t2) / (temp_max_rweed - temp_opt_rweed)
  else 
    return
  endif
  
  ! Photoperiod response for biotime accumulation
  if(daylen>photoperiod_rweed)then
    if(hsum>11.5 .and. hsum<16.)then
      ff = ff * exp((daylen-photoperiod_rweed) * log(1.-0.5))
    elseif(hsum>16. .and. hsum<20.5)then
      ff = ff * exp((daylen-photoperiod_rweed) * log(1.-0.6))
    endif
  endif

  hsum = hsum + ff*dt/(3600*24)      ! seconds to days
end subroutine heatsum_rweed
function regweed_normal_diurnal(dayJ,dayT,startDay,endDay) result(ff)
  integer, intent(in) :: &
    dayJ                ! julian date, daynumber
  real, intent(in) :: &
    dayT,             & ! seconds after sunrise
    startDay,endDay     ! season start/end (StartCDThr/EndCDThr), daynumber
  real, parameter :: &
    sqrt2 = sqrt(2.),&
    dayMid1 = 135*60, daySgm1 =  15*60*sqrt2, & ! 1st peak; sec since sunrise,sec*sqrt2
    dayMid2 = 285*60, daySgm2 = 100*60*sqrt2, & ! 2nd peak; sec since sunrise,sec*sqrt2
    dayfrac1= 0.4,    dayfrac2= 1.0-dayFrac1 ! % in the 1st/2nd peak 
  real(kind=8) :: t1, t2, s, ff, seasSgm

  ff=0.0
  if(startDay>endDay)return
  if(dayT<0.0 .or. dayT>36000.0)return  ! no nighttime emission

  ! Bimodal daily release
  t1 = dayT - dayMid1
  t2 = t1 + dt
  ff =      dayFrac1 * 0.5 * (fu_erf(t2/daySgm1) - fu_erf(t1/daySgm1))

  t1 = dayT - dayMid2
  t2 = t1 + dt
  ff = ff + dayFrac2 * 0.5 * (fu_erf(t2/daySgm2) - fu_erf(t1/daySgm2))

  ! Normal seas day by day
  t1 = dayJ - (startDay+endDay)*0.5  ! dayJ - seasMid
  t2 = t1 + 1.0
  s = (endDay-startDay)*0.25 * sqrt2 ! seasSgm * sqrt(2.)
  ff = ff * 0.5 * (fu_erf(t2/s) - fu_erf(t1/s))

contains
function fu_erf(x) result(y)
  ! error function in double precision
  real(kind=8) :: x, y, t, w
  integer :: k, i
  real(kind=8), dimension(0:64) :: &
    a=[  0.00000000005958930743,-0.00000000113739022964, 0.00000001466005199839,-0.00000016350354461960, &
         0.00000164610044809620,-0.00001492559551950604, 0.00012055331122299265,-0.00085483269811296660, &
         0.00522397762482322257,-0.02686617064507733420, 0.11283791670954881569,-0.37612638903183748117, &
         1.12837916709551257377,   & ! !(a(i), i = 0, 12) 
         0.00000000002372510631,-0.00000000045493253732, 0.00000000590362766598,-0.00000006642090827576, &
         0.00000067595634268133,-0.00000621188515924000, 0.00005103883009709690,-0.00037015410692956173, &
         0.00233307631218880978,-0.01254988477182192210, 0.05657061146827041994,-0.21379664776456006580, &
         0.84270079294971486929,   & !(a(i), i = 13, 25)
         0.00000000000949905026,-0.00000000018310229805, 0.00000000239463074000,-0.00000002721444369609, &
         0.00000028045522331686,-0.00000261830022482897, 0.00002195455056768781,-0.00016358986921372656, &
         0.00107052153564110318,-0.00608284718113590151, 0.02986978465246258244,-0.13055593046562267625, &
         0.67493323603965504676,   & ! (a(i), i = 26, 38)
         0.00000000000382722073,-0.00000000007421598602, 0.00000000097930574080,-0.00000001126008898854, &
         0.00000011775134830784,-0.00000111992758382650, 0.00000962023443095201,-0.00007404402135070773, &
         0.00050689993654144881,-0.00307553051439272889, 0.01668977892553165586,-0.08548534594781312114, &
         0.56909076642393639985,   & ! (a(i), i = 39, 51)
         0.00000000000155296588,-0.00000000003032205868, 0.00000000040424830707,-0.00000000471135111493, &
         0.00000005011915876293,-0.00000048722516178974, 0.00000430683284629395,-0.00003445026145385764, &
         0.00024879276133931664,-0.00162940941748079288, 0.00988786373932350462,-0.05962426839442303805, &
         0.49766113250947636708 ], & ! (a(i), i = 52, 64)
    b=[ -0.00000000029734388465, 0.00000000269776334046,-0.00000000640788827665,-0.00000001667820132100, &
        -0.00000021854388148686, 0.00000266246030457984, 0.00001612722157047886,-0.00025616361025506629, &
         0.00015380842432375365, 0.00815533022524927908,-0.01402283663896319337,-0.19746892495383021487, & 
         0.71511720328842845913, & ! (b(i), i = 0, 12)
        -0.00000000001951073787,-0.00000000032302692214, 0.00000000522461866919, 0.00000000342940918551, &
        -0.00000035772874310272, 0.00000019999935792654, 0.00002687044575042908,-0.00011843240273775776, &
        -0.00080991728956032271, 0.00661062970502241174, 0.00909530922354827295,-0.20160072778491013140, &
         0.51169696718727644908, & ! (b(i), i = 13, 25)
         0.00000000003147682272,-0.00000000048465972408, 0.00000000063675740242, 0.00000003377623323271, &
        -0.00000015451139637086,-0.00000203340624738438, 0.00001947204525295057, 0.00002854147231653228, &
        -0.00101565063152200272, 0.00271187003520095655, 0.02328095035422810727,-0.16725021123116877197, &
         0.32490054966649436974, & ! (b(i), i = 26, 38)
         0.00000000002319363370,-0.00000000006303206648,-0.00000000264888267434, 0.00000002050708040581, &
         0.00000011371857327578,-0.00000211211337219663, 0.00000368797328322935, 0.00009823686253424796, &
        -0.00065860243990455368,-0.00075285814895230877, 0.02585434424202960464,-0.11637092784486193258, &
         0.18267336775296612024, & ! (b(i), i = 39, 51)
        -0.00000000000367789363, 0.00000000020876046746,-0.00000000193319027226,-0.00000000435953392472, &
         0.00000018006992266137,-0.00000078441223763969,-0.00000675407647949153, 0.00008428418334440096, &
        -0.00017604388937031815,-0.00239729611435071610, 0.02064129023876022970,-0.06905562880005864105, &
         0.09084526782065478489 ]  ! (b(i), i = 52, 64)
  w = abs(x)
  if(w<2.2) then
    t = w * w
    k = int(t)
    t = t - k
    k = k * 13
    y = ((((((((((((a(k)*t + a(k+1))*t + a(k+2))*t + a(k+3))*t + a(k+4))*t + &
                 a(k+5))*t + a(k+6))*t + a(k+7))*t + a(k+8))*t + a(k+9))*t + &
                 a(k+10))*t + a(k+11))*t + a(k+12))*w
  elseif(w<6.9) then
    k = int(w)
    t = w - k
    k = 13 * (k - 2)
    y = (((((((((((b(k)*t + b(k+1))*t + b(k+2))*t + b(k+3))*t + b(k+4))*t + &
                b(k+5))*t + b(k+6))*t + b(k+7))*t + b(k+8))*t + b(k+9))*t + &
                b(k+10))*t + b(k+11))*t + b(k+12)
    y = y * y
    y = y * y
    y = y * y
    y = 1d0 - y * y
  else
    y = 1d0
  end if
  if(x<0) y = -y
end function fu_erf
end function regweed_normal_diurnal
function f_gamma_w_tails(relTime,relDt) result(ff)
! Returns the pollen prepared for release assuming the modified "taily" Gamma distribution
! of the season. Tails are the reason for many parameters: have to describe the main peak
! via gamma-type distribution, and both elevated tails via add-on corrections.
! formula: rate(x)=exp(-a_1/beta)* sum(scale_i * a_i^power_i), i=1:3
!          where a_i = MAX(x-timesRel_i,0)
  real, intent(in) :: relTime,relDt ! normalised time, normalised timestep
  real             :: ff
! Fitting parameters
  integer, parameter:: nTerms= 3
  real, parameter::&
    beta = 0.155, &! [main season shape,correction 1,correction 2]
    timesRel(nTerms)=[ 0.164,-0.3 ,-0.6 ],&
    scales  (nTerms)=[13.1  ,12.6 , 0.25],&
    powers  (nTerms)=[ 1.16 , 2.8 , 1.3 ]
! Local variable
  real :: a1

  ff = 0.0
  a1 = DIM(relTime,timesRel(1)) ! same as MAX(relTime-timesRel(1),0.0)
  if(a1>beta*10.0)return ! too far from the season peak, as decided by timesRel(1)
  ! Be careful: the rise of the function can be quite steep
  a1=exp(-a1/beta)
! ff=0.5*relDt*(rate(a1,relTime)+rate(a1,relTime+relDt)) ! trapezoidal integration
! WRONG integration to match SILAM code
  ff=0.5*relDt*(rate(a1,relTime)+rate(a1*a1,relTime+relDt))
contains
real function rate(a,x)
  real, intent(in) :: a,x
  rate=a*sum(scales*DIM(x,timesRel)**powers,MASK=(x>timesRel))
end function rate
end function f_gamma_w_tails
function f_wind(u,wstar) result(ff)
! Pollen emission increases with wind speeds up to 5 m/s (u_sat).
! This term cannot be higher than 1.5 (u_max).
  real, intent(in) :: u,wstar
  real             :: ff
  real, parameter  :: u_max=1.5,u_sat=5.0
  ff = u_max - exp(-(u+wstar)/u_sat)
end function f_wind
function f_cond(x,x_min,x_max) result(ff)
! This function is used for both humitidy and rain, as too much
! humidity and rain stops pollen release
  real, intent(in) :: x,x_min,x_max
  real             :: ff
  if(x>x_max) then
    ff=0.0
  elseif(x<x_min)then
    ff=1.0
  else
    ff=(x_max-x)/(x_max-x_min)
  end if
end function f_cond
function f_in(h,h_c,prob) result(ff)
! takes in account the uncertainity that all the trees start
! to release pollen at the same time
  real, intent(in) :: H,H_c,prob
  real             :: ff
  if(h<(1.0-prob)*h_c) then
    ff=0.0
  elseif(h>(1.0+prob)*h_c)then
    ff=1.0
  else
    ff=(h-(1.0-prob)*h_c)/(2.0*prob*h_c)
  end if
end function f_in
function f_out(R,N_tot,prob) result(ff)
! takes in account the uncertainity that all the trees stop
! releasing pollen at the same time
  real, intent(in) :: R,N_tot,prob
  real             :: ff
  if(R<(1.0-prob)*N_tot) then
    ff=1.0
  elseif (R>(1.0+prob)*N_tot) then
    ff=0.0
  else
    ff=1.0-(R-(1.0-prob)*N_tot)/(2.0*prob*N_tot)
  end if
end function f_out
function f_fade_in(value_rel,uncert_rel_) result(ff)
! Computes the linear fade-in function. It is 0 at vaule_rel = 1.-uncert_rel
! Grows to 0.5 at vaule_rel = 1 and grows to 1 at value_rel = 1.+uncert_rel
  real, intent(in) :: value_rel,uncert_rel_
  real             :: uncert_rel,ff

  uncert_rel = MAX(uncert_rel_,1e-5) ! avoid zero uncertainty
  ff = MIN(1.0,MAX(0.0,(value_rel-1.0+uncert_rel)/(2.0*uncert_rel)))
end function f_fade_in
function f_fade_out(value_rel,uncert_rel_) result(ff)
! Computes the linear fade-in function. It is 0 at vaule_rel = 1.-uncert_rel
! Grows to 0.5 at vaule_rel = 1 and grows to 1 at value_rel = 1.+uncert_rel
  real, intent(in) :: value_rel,uncert_rel_
  real             :: uncert_rel,ff

  uncert_rel = MAX(uncert_rel_,1e-5) ! avoid zero uncertainty
  ff = MIN(1.0,MAX(0.0,(1.0+uncert_rel-value_rel)/(2.0*uncert_rel)))
end function f_fade_out
!-------------------------------------------------------------------------!
integer function getRecord(fileName,findDate,fatal)
  character(len=*), intent(in) :: fileName
  type(date), intent(in) :: findDate
  logical, intent(in) :: fatal

  if(MasterProc) getRecord=startRecord()
  call MPI_BCAST(getRecord,1,MPI_INTEGER,MasterPE,MPI_COMM_CALC,IERROR)
  if(getRecord<1)&  ! switch off Pollen_mod when restart file is corrupted
    call MPI_BCAST(USES%POLLEN,1,MPI_LOGICAL,MasterPE,MPI_COMM_CALC,IERROR)
contains
function startRecord() result(nstart)
  integer :: nstart

  logical :: fexist
  integer :: nread
  real(kind=8), parameter :: sec2day=1e0/(24.0*3600.0)
  real, dimension(366), save :: fdays=-1
  real :: ncday(0:1)

  ! ensure file exists
  inquire(file=filename,exist=fexist)
  if(.not.fexist)then
    call CheckStop(fatal,"File not found: "//trim(filename))
    if(MasterProc) write(*,*)&
      "Warning Pollen dump file not found: "//trim(filename)
    call PrintLog("WARNING: Pollen_mod cold start",MasterProc)
    nstart=-1
    return
  end if

  ! ensure file has records
  nread=-1                                ! read all
  fdays(:)=-1.0                           ! times records
  call ReadTimeCDF(filename,fdays,nread)  ! in fname
  if(nread<1)then
    call CheckStop(fatal,"Corrupted file: "//trim(filename))
    if(MasterProc) write(*,*)&
      "Warning Pollen dump file corrupted: "//trim(filename)
    call PrintLog("WARNING: Pollen_mod forced OFF",MasterProc)
    USES%POLLEN=.false.
    nstart=-1
    return
  end if

  ! look for current_date in fdays (time var read from filename)
  call date2nctime(findDate,ncday(1))
  ncday(0)=ncday(1)-dt*sec2day    ! to avoid rounding errors
  ncday(1)=ncday(1)+dt*sec2day    ! look for records in +/- 1dt

  nstart=MINLOC(fdays(:nread),DIM=1,&
    MASK=(fdays(:nread)>=ncday(0)).and.(fdays(:nread)<ncday(1)))

  ! check if we got a match
  if(nstart<1)then  ! ifort compiler needs option: -assume noold_maxminloc
    call CheckStop(fatal,&
      "No records for"//date2string(" YYYY-MM-DD hh:mm ",findDate)//&
      "found in "//trim(filename))
    if(MasterProc) write(*,*)&
      "No records for"//date2string(" YYYY-MM-DD hh:mm ",findDate)//&
      "found in "//trim(filename)
    call PrintLog("WARNING: Pollen_mod cold start",MasterProc)
    nstart=-1
    return
  end if
end function startRecord
end function getRecord
!-------------------------------------------------------------------------!
subroutine pollen_read()
! Read pollen for forecast restart file: heatsum and pollenrest fields
  character(len=*), parameter :: dfmt="('Read: ',A20,' [',L1,'].')"
  character(len=20) :: spc
  logical, save :: first_call = .true.
  logical :: found
  integer :: nstart,g
  character(len=len(template_read)) :: filename
  real,allocatable, dimension(:,:,:) :: data ! Data arrays

  if(MODE_READ=='NONE') return
  if(.not.checkdates(daynumber,"pollen")) return
  if(.not.first_call) return
  first_call=.false.

  call Config_Pollen()
  filename=date2string(template_read,current_date)
  nstart=getRecord(filename,current_date,MODE_READ=='START')
  if(nstart<1) return
  if(MasterProc)&
    write(*,"(3(A,1X),I0)") "Read Pollen dump",trim(filename),"record",nstart

  allocate(data(LIMAX,LJMAX,KMAX_MID))
  do g=1,size(POLLEN_GROUP)
    spc=trim(POLLEN_GROUP(g))
!------------------------
! pollen adv (not written by Nest_mod)
!------------------------
    call GetCDF_modelgrid(trim(spc),filename,data,&
          1,KMAX_MID,nstart,1,needed=MODE_READ=='START',found=found)
    if(DEBUG%POLLEN.and.MasterProc) write(*,dfmt)spc,found
    if(found)xn_adv(iadv(g),:,:,:)=data(:,:,:)
!------------------------
! pollen_rest
!------------------------
    spc=trim(POLLEN_GROUP(g))//'_rest'
    call GetCDF_modelgrid(trim(spc),filename,data(:,:,1),&
        1,1,nstart,1,needed=MODE_READ=='START',found=found)
    if(DEBUG%POLLEN.and.MasterProc) write(*,dfmt)spc,found
    if(found)R(:,:,g)=N_TOT(g)-data(:,:,1)
!------------------------
! heatsum
!------------------------
    if(g>size(heatsum,DIM=3))cycle
    spc=trim(POLLEN_GROUP(g))//'_heatsum'
    call GetCDF_modelgrid(trim(spc),filename,heatsum(:,:,g),&
        1,1,nstart,1,needed=MODE_READ=='START',found=found)
    if(DEBUG%POLLEN.and.MasterProc) write(*,dfmt)spc,found
  end do
  deallocate(data)
end subroutine pollen_read
!-------------------------------------------------------------------------!
subroutine pollen_dump()
! Write pollen for forecast restart file: heatsum and pollenrest fields
  character(len=len(template_write)) :: filename
  character(len=*), parameter :: dfmt="('Write: ',A20,' [',A,'].')"
  character(len=20) :: spc
  logical     :: fexist,create_var
  integer     :: ncfileID,i,g
  type(Deriv) :: def1
  real,allocatable, dimension(:,:,:) :: data ! Data arrays

  if(.not.checkdates(daynumber,"pollen")) return
  call CheckStop(allocated(R).NEQV.allocated(heatsum),&
    "Pollen: rest/heatsum allocated error")
  if(.not.(allocated(R).and.allocated(heatsum))) return

  filename = date2string(template_write,current_date)
  if(MasterProc) write(*,*) "Write Pollen dump ",trim(filename)
  inquire(file=filename,exist=fexist)

  def1%avg=.false.                  ! not used
  def1%index=0                      ! not used
  def1%scale=1.0                    ! not used
  def1%iotype=''                    ! not used

  allocate(data(LIMAX,LJMAX,KMAX_MID))
  ncfileID=-1 ! must be <0 as initial value
  do i=1,2                          ! do first one loop to define the fields,
    create_var=(i==1)               ! without writing them (for performance purposes),
    if(all([create_var,fexist,.not.overwrite]))&  ! if the file does not exists already
      cycle
    overwrite=overwrite.and.create_var

    do g=1,size(POLLEN_GROUP)
      spc=trim(POLLEN_GROUP(g))
!------------------------
! pollen adv (not written by Nest_mod)
!------------------------
      def1%class='Advected'           ! written
      def1%unit='mix_ratio'           ! written
      def1%name=trim(spc)             ! written
      if(DEBUG%POLLEN.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      if(.not.create_var)data=xn_adv(iadv(g),:,:,:)
      call Out_netCDF(IOU_INST,def1,3,KMAX_MID,data,1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,overwrite=overwrite,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
      overwrite=.false.
!------------------------
! pollen_rest
!------------------------
      def1%class='pollen_out'         ! written
      def1%unit='pollengrains'        ! written
      def1%name=trim(spc)//'_rest'    ! written
      if(DEBUG%POLLEN.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      data(:,:,1)=N_TOT(g)-R(:,:,g)
      call Out_netCDF(IOU_INST,def1,2,1,data(:,:,1),1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
!------------------------
! heatsum
!------------------------
      if(g>size(heatsum,DIM=3))cycle
      def1%class='pollen_out'         ! written
      def1%unit='degreedays'          ! written
      if(g==iALDER)&
        def1%unit='degreehours'
      def1%name=trim(spc)//'_heatsum' ! written
      if(DEBUG%POLLEN.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      call Out_netCDF(IOU_INST,def1,2,1,heatsum(:,:,g),1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
    end do
  end do
  if(MasterProc)&
    call CheckNC(nf90_close(ncFileID),"close:"//trim(filename))
  ! ensure time record can be found, fatal error if not
  i=getRecord(filename,current_date,.true.) ! MPI_BCAST inside
  if(DEBUG%POLLEN.and.MasterProc) &
    write(*,"(3(A,1X),I0)") "Found Pollen dump",trim(filename),"record",i
  deallocate(data)
end subroutine pollen_dump
end module Pollen_mod
