! <GridValues_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
Module GridValues_mod
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Define parameters, variables and transformnations associated with grid and
!  projection.
!
! Nomenclature:
! fulldomain is the largest grid, usually where metdata is defined.
! rundomain is a grid where the run is performed, smaller than fulldomain.
! subdomain: the domain covered by one MPI process or processor.
! restricted domain is a grid smaller than rundomain, where data is outputed;
!  (the restricted domains are for instance, fullrun_DOMAIN,month_DOMAIN,
!  day_DOMAIN,hour_DOMAIN).
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

use CheckStop_mod,           only: CheckStop,StopAll,check=>CheckNC
use Functions_mod,           only: great_circle_distance
use Io_Nums_mod,             only: IO_LOG,IO_TMP
use MetFields_mod
use Config_module,      only: &
     KMAX_BND, KMAX_MID, & ! vertical extent
     MasterProc,NPROC,IIFULLDOM,JJFULLDOM,RUNDOMAIN, JUMPOVER29FEB,&
     PT,Pref,NMET,USE_EtaCOORDINATES,MANUAL_GRID,USE_WRF_MET_NAMES,&
     startdate,NPROCX,NPROCY,Vertical_levelsFile,&
     EUROPEAN_settings, GLOBAL_settings,USES,FORCE_PFT_MAPS_FALSE,&
     USE_SOILNOx
use Debug_module, only:  DEBUG   ! -> DEBUG%GRIDVALUES

use MPI_Groups_mod!, only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_LOGICAL, &
                 !         MPI_MIN, MPI_MAX, &
                 !         MPI_COMM_CALC, MPI_COMM_WORLD, MPISTATUS, IERROR, &
                 !         ME_MPI, NPROC_MPI, ME_CALC, largeLIMAX,largeLJMAX
use Par_mod, only : &
     GIMAX,GJMAX,        & ! Size of rundomain
     IRUNBEG,JRUNBEG,    & ! start of rundomain in fulldomain coordinates
     gi0,gj0,            & ! rundomain coordinates of subdomain lower l.h. corner
     gi1,gj1,            & ! rundomain coordinates of subdomain uppet r.h. corner
     limax,ljmax,        & ! max i,j in this subdomain (can differ on other subdomains)
     li0,li1,lj0,lj1,    & ! start and end of i,j excluding outer frame of rundomain.
                           ! li0=1 or 2, li1=limax or limax-1
     me,                 & ! local processor
     neighbor,WEST,EAST,SOUTH,NORTH,NOPROC,  &
     parinit,parinit_groups,  &
     MAXLIMAX,MAXLJMAX,MINLIMAX,MINLJMAX,tljmax,tlimax
use PhysicalConstants_mod,     only: GRAV,PI,EARTH_RADIUS,deg2rad,rad2deg
use TimeDate_mod,              only: current_date,date,Init_nmdays,nmdays
use TimeDate_ExtraUtil_mod,    only: date2string
use InterpolationRoutines_mod, only: inside_1234
use netcdf,                   only: &
  NF90_OPEN,NF90_NOWRITE,NF90_NOERR,NF90_CLOSE,&
  NF90_GET_ATT,NF90_GLOBAL,NF90_INQ_DIMID,NF90_INQUIRE_DIMENSION,&
  NF90_INQ_VARID,NF90_GET_VAR

implicit none
private

!-- contains subroutine:
public :: DefDebugProc ! =>  sets debug_proc, debug_li, debug_lj
public :: ij2lbm  ! polar stereo grid to longitude latitude
public :: lb2ijm  ! longitude latitude to grid in polar stereo
public :: ij2ijm  ! polar grid1 to polar grid2
public :: lb2ij   ! longitude latitude to (i,j) in any grid projection
public :: ij2lb   ! polar stereo grid to longitude latitude
public :: lb_rot2lb !rotated lon lat to lon lat
public :: lb2UTM ! lon lat to Transverse Mercator
public :: UTM2lb ! Transverse Mercator to lon lat

interface lb2ij
  module procedure lb2ij_real,lb2ij_int
end interface
private :: lb2ij_real,lb2ij_int

public :: coord_check   ! normalize longitudes

public :: &
  coord_in_gridbox,  &  ! Are coord (lon/lat) inside gridbox(i,j)?
  coord_in_processor,&  ! Are coord (lon/lat) inside local domain?
  coord_in_domain       ! Are coord (lon/lat) inside "domain" (full, run or sub)?

public :: RestrictDomain ! mask from full domain to rundomain

public :: GridRead
public :: extendarea_N ! returns array which includes neighbours from other subdomains
public :: set_EuropeanAndGlobal_Config
public :: remake_vertical_levels_interpolation_coeff
public :: Meteo_Get_KMAXMET

private :: Alloc_GridFields
private :: GetFullDomainSize
private :: find_poles

!** 1) Public (saved) Variables from module:

! Polar stereographic projection parameters
real, public, save :: &
  xp=0.0, yp=1.0,  & ! Coordinates of North pole
  fi=0.0,          & ! projections rotation angle around y axis
  AN=1.0,          & ! Distance on the map from pole to equator (No. of cells)
  GRIDWIDTH_M=1.0, & ! width of grid at ref_latitude, in meters
  ref_latitude =60.  ! latitude at which projection is true (degrees)

! Rotated_Spherical grid prarameters
real, public, save :: &
  grid_north_pole_latitude,grid_north_pole_longitude,&
  dx_rot,dx_roti,x1_rot,y1_rot

! Lambert conformal projection parameters
real, public, save :: &
  lon0_lambert,&  ! reference longitude, also called phi, at which y=0 if lat=lat0_lambert
  lat0_lambert,&  ! reference latitude, at which x=0
  lat_stand1_lambert,&! standard latitude at which mapping factor=1
  lat_stand2_lambert,&! second standard latitude
  y0_lambert,&    ! reference y coordinate, also called rho0
  k_lambert,&     ! also called n, = sin(lat_stand1_lambert)
  earth_radius_lambert,&! earth_radius used to define x and y in the met file. NOT USED
  F_lambert,&     ! normalization constant = cos(dr*lat0_lambert)*tan(PI/4+dr2*lat0_lambert)**k_lambert/k_lambert
  x1_lambert,&    ! x value at i=1
  y1_lambert      ! y value at j=1

!/ Variables to define full-domain (fdom) coordinates of local i,j values,
!  and reciprocal variables.
integer, public, allocatable, save, dimension(:) :: &
  i_fdom,j_fdom, & ! fdom coordinates of local i,j
  i_local,j_local  ! local coordinates of full-domain i,j

!Parameters for Vertical Hybrid coordinates:
real, public, save,allocatable,  dimension(:) ::  &
  A_bnd,B_bnd,&         ! first [Pa],second [1] constants at layer boundary
                        ! (i.e. half levels in EC nomenclature)
  A_bnd_met,B_bnd_met,& ! first [Pa],second [1] constants at layer boundary
                        ! (i.e. half levels in EC nomenclature)
  A_mid,B_mid,&         ! first [Pa],second [1] constants at middle of layer
                        ! (i.e. full levels in EC nomenclature)
  dA,dB,&               ! A_bnd(k+1)-A_bnd(k) [Pa],B_bnd(k+1)-B_bnd(k) [1]
                        ! P = A + B*PS; eta = A/Pref + B
  dEta_i,&              ! 1/deta = 1/(dA/Pref + dB)
  Eta_bnd,Eta_mid,&     ! boundary,midpoint of eta layer
  sigma_bnd,sigma_mid   ! boundary,midpoint of sigma layer

real, public, save,allocatable,  dimension(:,:) :: &
  glon     ,glat    ,&  !longitude,latitude of gridcell centers
  gl_stagg ,gb_stagg,&  !longitude,latitude of gridcell corners
      !NB: gl_stagg,gb_stagg are here defined as the average of the four
      !    surrounding glat,glon.
      !    These differ slightly from the staggered points in the (i,j) grid.
  rot_angle

real, public, save :: gbacmax,gbacmin,glacmax,glacmin

! EMEP grid definitions (old and official)
real, public, parameter :: &
  xp_EMEP_official=8.,yp_EMEP_official=110.0,fi_EMEP=-32.0,&
  ref_latitude_EMEP=60.0,GRIDWIDTH_M_EMEP=50000.0,&
  an_EMEP=237.7316364, &! = 6.370e6*(1.0+0.5*sqrt(3.0))/50000.
  xp_EMEP_old=43.0,yp_EMEP_old=121.0

!*** Map factor stuff:
real, public, save,allocatable, dimension(:,:) ::  &
  xm_i,     & ! map-factor in i direction, between cell j and j+1
  xm_j,     & ! map-factor in j direction, between cell i and i+1
  xm2,      & ! xm*xm: area factor in the middle of a cell (i,j)
  xmd,      & ! 1/xm2
  xm2ji,xmdji
!vertical "map factors"
real, public, save, allocatable,dimension(:)  ::  dhs1, dhs1i, dhs2i


!*** Grid Area
real, public, save,allocatable, dimension(:,:) :: GridArea_m2

integer, public, save :: &
  debug_li=-99, debug_lj=-99         ! Local Coordinates of debug-site
logical, public, save :: debug_proc  ! Processor with debug-site

character(len=100),public  :: projection
integer, public, parameter :: MIN_ADVGRIDS = 5 !minimum size of a subdomain
integer, public :: Poles(2)       ! Poles(1)=1 if North pole is found, Poles(2)=1:SP
integer, public :: Pole_Singular  ! Pole_included=1 or 2 if the grid include
! at least one pole and has lat lon projection
logical, public :: Grid_Def_exist

integer, allocatable, save, public :: k1_met(:),k2_met(:)
real, allocatable, save, public :: x_k1_met(:)
logical, public, save ::  External_Levels_Def=.false.
integer, public, save :: KMAX_MET !number of vertical levels from the meteo files
real, private :: u(2),v(2)!array for temporary use

contains

subroutine GridRead(meteo,cyclicgrid)
  ! the subroutine reads the grid parameters (projection, resolution etc.)
  ! defined by the meteorological fields
  implicit none

  character(len=*),intent(in):: meteo   ! template for meteofile
  integer,  intent(out)      :: cyclicgrid
  integer                    :: nyear,nmonth,nday,nhour,k,ios
  integer                    :: MIN_GRIDS
  character(len=len(meteo))  :: filename !name of the input file
  logical :: Use_Grid_Def=.false.!Experimental for now

  nyear=startdate(1)
  nmonth=startdate(2)
  nday=startdate(3)
  nhour=startdate(4)
  current_date = date(nyear, nmonth, nday, nhour, 0 )
  call Init_nmdays( current_date, JUMPOVER29FEB)

  !*********initialize grid parameters*********
  if(MANUAL_GRID)then
    ! define the grid parameter manually (explicitely)
    if(MasterProc)write(*,*)'DEFINING GRID MANUALLY!'
    ! must set all parameters... see example in version from 20151019 (or before)
  else
    ! NOT MANUAL GRID

    !check first if grid is defined in a separate file:
    filename='Grid_Def.nc'
    inquire(file=filename,exist=Grid_Def_exist)
    Grid_Def_exist=Grid_Def_exist.and.Use_Grid_Def
    if(Grid_Def_exist)then
     if(MasterProc)write(*,*)'Found Grid_Def! ',trim(filename)
    else
      if(MasterProc.and.Use_Grid_Def)&
        write(*,*)'Did not found Grid_Def ',trim(filename)
      filename=date2string(meteo,startdate,mode='YMDH')
    end if
    if(MasterProc)write(*,*)'reading domain sizes from ',trim(filename)

    call GetFullDomainSize(filename,IIFULLDOM,JJFULLDOM,KMAX_MET,projection)

    KMAX_MID=0!initialize
    open(IO_TMP,file=Vertical_levelsFile,action="read",iostat=ios)
    if(ios==0)then
      ! define own vertical coordinates
      if(MasterProc)&
        write(*,*)'Define vertical levels from ',trim(Vertical_levelsFile)
      read(IO_TMP,*)KMAX_MID
      if(MasterProc)write(*,*)KMAX_MID, 'vertical levels '
      External_Levels_Def=.true.
      ! Must use eta coordinates
      if(MasterProc.and..not.USE_EtaCOORDINATES)&
        write(*,*)'WARNING: using hybrid levels even if not asked to! '
      USE_EtaCOORDINATES=.true.
    else
      if(MasterProc)write(*,*)'WARNING: could not open '//trim(Vertical_levelsFile)
      External_Levels_Def=.false.
      close(IO_TMP)
      KMAX_MID=KMAX_MET
    end if

    KMAX_BND=KMAX_MID+1

    allocate(A_bnd(KMAX_BND),B_bnd(KMAX_BND))
    allocate(A_mid(KMAX_MID),B_mid(KMAX_MID))
    allocate(dA(KMAX_MID),dB(KMAX_MID),dEta_i(KMAX_MID))
    allocate(sigma_bnd(KMAX_BND),sigma_mid(KMAX_MID))
    allocate(Eta_bnd(KMAX_BND),Eta_mid(KMAX_MID))
    allocate(i_local(IIFULLDOM))
    allocate(j_local(JJFULLDOM))

    ! set RUNDOMAIN default values where not defined
    if(RUNDOMAIN(1)<1)RUNDOMAIN(1)=1
    if(RUNDOMAIN(2)<1 .or. RUNDOMAIN(2)>IIFULLDOM) RUNDOMAIN(2)=IIFULLDOM
    if(RUNDOMAIN(3)<1)RUNDOMAIN(3)=1
    if(RUNDOMAIN(4)<1 .or. RUNDOMAIN(4)>JJFULLDOM) RUNDOMAIN(4)=JJFULLDOM
    if(MasterProc)then
55    format(A,I5,A,I5)
      write(*,55)     'FULLDOMAIN has sizes ',IIFULLDOM,' X ',JJFULLDOM
      write(IO_LOG,55)'FULLDOMAIN has sizes ',IIFULLDOM,' X ',JJFULLDOM
      write(*,55)     'RUNDOMAIN  x coordinates from ',RUNDOMAIN(1),' to ',RUNDOMAIN(2)
      write(IO_LOG,55)'RUNDOMAIN  x coordinates from ',RUNDOMAIN(1),' to ',RUNDOMAIN(2)
      write(*,55)     'RUNDOMAIN  y coordinates from ',RUNDOMAIN(3),' to ',RUNDOMAIN(4)
      write(IO_LOG,55)'RUNDOMAIN  y coordinates from ',RUNDOMAIN(3),' to ',RUNDOMAIN(4)
    end if

    call find_poles(filename,Pole_Singular)

    MIN_GRIDS=5
    if(NPROC==NPROC_MPI)then
      ! partition into subdomains
      call parinit(MIN_GRIDS,Pole_Singular)        ! subdomains sizes and position
    else
      ! partition into largesubdomains and subdomains
      call parinit_groups(MIN_GRIDS,Pole_Singular) ! subdomains sizes and position
    end if

    call Alloc_MetFields(LIMAX,LJMAX,KMAX_MID,KMAX_BND,NMET)

    if(ME_CALC>=0)then
      call Alloc_GridFields(LIMAX,LJMAX,KMAX_MID,KMAX_BND)
    else
      call Alloc_GridFields(largeLIMAX,largeLJMAX,KMAX_MID,KMAX_BND)
    end if

    if(ME_CALC>=0)then
      call Getgridparams(LIMAX,LJMAX,filename,cyclicgrid)
      ! defines i_fdom,j_fdom,i_local,j_local,Cyclicgrid,North_pole,Poles
      ! GRIDWIDTH_M, glon, glat, xm_i,xm_j,xm2,xmd,xmdji,xm2ji,gl_stagg,gb_stagg
      ! P0,A_bnd_met,B_bnd_met,A_bnd,B_bnd,A_mid,B_mid,sigma_mid,sigma_bnd
      ! for Stereographic projection:
      !   ref_latitude,fi,xp,yp,AN
      ! for lon lat projection:
      !   no additional parameters
      ! for Rotated_Spherical projection: 
      !  grid_north_pole_latitude,grid_north_pole_longitude,x1_rot,y1_rot,dx_rot,dx_roti
    else
      call Getgridparams(largeLIMAX,largeLJMAX,filename,cyclicgrid)
    end if


    if(ios==0)close(IO_TMP)

  end if
end subroutine GridRead

subroutine GetFullDomainSize(filename,IIFULLDOM,JJFULLDOM,KMAX,projection)
  ! Get input grid sizes

  implicit none

  character (len = *), intent(in) ::filename
  integer, intent(out):: IIFULLDOM,JJFULLDOM,KMAX
  character (len = *), intent(out) ::projection

  integer :: status,ncFileID,idimID,jdimID, kdimID,timeDimID
  integer :: GIMAX_file,GJMAX_file,KMAX_file,wrf_proj_code
  real :: wrf_POLE_LAT=0.0
  character (len = 30) ::MAP_PROJ_CHAR
  integer :: NTime_Read
  real :: TimesInDays(1000)


  if(ME_MPI==0)then
     print *,'Defining grid properties from ',trim(filename)
     ! open an existing netcdf dataset
     status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
     if(status/=nf90_noerr) then
        print *,'not found',trim(filename)
        call StopAll("GridValues: File not found:"//trim(filename))
     end if

     projection=''
     status = nf90_get_att(ncFileID,nf90_global,"projection",projection)
     if(status/=nf90_noerr) then
        ! WRF projection format
        call check(nf90_get_att(ncFileID,nf90_global,"MAP_PROJ",wrf_proj_code))
        if(.not.USE_WRF_MET_NAMES .and. MasterProc)write(*,*)'Assuming WRF metdata'
        USE_WRF_MET_NAMES = .true.
        select case(wrf_proj_code)
        case(6)
           status = nf90_get_att(ncFileID,nf90_global,"POLE_LAT",wrf_POLE_LAT)
           if(status==nf90_noerr) then
              write(*,*)"POLE_LAT", wrf_POLE_LAT
              if(abs(wrf_POLE_LAT-90.0)<0.001)then
                 projection='lon lat'
              else
                 projection='Rotated_Spherical'
              end if
           else
              write(*,*)"POLE_LAT not found"
              projection='lon lat'
           end if
        case(2)
           projection='Stereographic'
        case(1)
           projection='lambert'
           call check(nf90_get_att(ncFileID,nf90_global,"MAP_PROJ_CHAR",MAP_PROJ_CHAR))
           write(*,*)"wrf projection: "//trim(MAP_PROJ_CHAR)
        case default
           call CheckStop("Projection not recognized")
        end select
     end if

     ! put into emep standard
     if(trim(projection)=='Polar Stereographic')projection='Stereographic'

     if(trim(projection)=='Rotated_Spherical'.or.trim(projection)=='rotated_spherical'&
          .or.trim(projection)=='rotated_pole'.or.trim(projection)=='rotated_latitude_longitude')then
        projection='Rotated_Spherical'
     end if

     if(trim(projection)=='lambert_conformal_conic'.or.trim(projection)=='Lambert Conformal')projection='lambert'

     write(*,*)'projection: ',trim(projection)

     ! get dimensions id
     if(trim(projection)=='Stereographic') then
        status = nf90_inq_dimid(ncid=ncFileID, name="i", dimID=idimID)
        if(status/=nf90_noerr)& ! WRF  format
             call check(nf90_inq_dimid(ncid=ncFileID, name="west_east", dimID=idimID))
        status = nf90_inq_dimid(ncid=ncFileID, name="j", dimID=jdimID)
        if(status/=nf90_noerr)& ! WRF  format
             call check(nf90_inq_dimid(ncid=ncFileID, name="south_north", dimID=jdimID))
     elseif(trim(projection)=='lambert') then
        status = nf90_inq_dimid(ncid=ncFileID, name="x", dimID=idimID)
        if(status/=nf90_noerr)& ! WRF format
             call check(nf90_inq_dimid(ncid=ncFileID, name="west_east", dimID=idimID))
        status = nf90_inq_dimid(ncid=ncFileID, name="y", dimID=jdimID)
        if(status/=nf90_noerr)& ! WRF format
             call check(nf90_inq_dimid(ncid=ncFileID, name="south_north", dimID=jdimID))
        write(*,*)'x y dimensions'
     elseif(trim(projection)==trim('lon lat')) then
        status=nf90_inq_dimid(ncid=ncFileID, name="lon", dimID=idimID)
        if(status/=nf90_noerr)& ! WRF format
             call check(nf90_inq_dimid(ncid=ncFileID, name="west_east", dimID=idimID))
        status=nf90_inq_dimid(ncid=ncFileID, name="lat", dimID=jdimID)
        if(status/=nf90_noerr)& ! WRF format
             call check(nf90_inq_dimid(ncid=ncFileID, name="south_north", dimID=jdimID))
     else
        ! write(*,*)'GENERAL PROJECTION ',trim(projection)
        status=nf90_inq_dimid(ncid=ncFileID, name="i", dimID = idimID)
        if(status/=nf90_noerr)& ! WRF format
             call check(nf90_inq_dimid(ncid=ncFileID, name="west_east", dimID=idimID))
        status=nf90_inq_dimid(ncid=ncFileID, name="j", dimID = jdimID)
        if(status/=nf90_noerr)& ! WRF format
             call check(nf90_inq_dimid(ncid=ncFileID, name="south_north", dimID=jdimID))
     end if

     status=nf90_inq_dimid(ncid=ncFileID, name="k", dimID=kdimID)
     if(status/=nf90_noerr)then
        status=nf90_inq_dimid(ncid=ncFileID, name="lev", dimID=kdimID)!hybrid coordinates
        if(status/=nf90_noerr) then
           status=nf90_inq_dimid(ncid=ncFileID, name="hybrid", dimID=kdimID)!hybrid coordinates
           if(status/=nf90_noerr) then ! WRF format
              call check(nf90_inq_dimid(ncid=ncFileID, name="bottom_top", dimID=kdimID))
           end if
        end if
     end if

     !get dimensions length
     call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_file))
     call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_file))
     call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_file))
     write(*,*)'dimensions input grid:',GIMAX_file,GJMAX_file,KMAX_file!,Nhh

     IIFULLDOM=GIMAX_file
     JJFULLDOM=GJMAX_file
     KMAX     =KMAX_file

     call check(nf90_close(ncFileID))
  end if

  CALL MPI_BCAST(USE_WRF_MET_NAMES ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(IIFULLDOM ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(JJFULLDOM ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(KMAX ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(projection ,len(projection),MPI_BYTE,0,MPI_COMM_WORLD,IERROR)

end subroutine GetFullDomainSize

subroutine find_poles(filename,Pole_Singular)
  ! defines if there is a singularity at poles

  implicit none
  character (len = *), intent(in) ::filename
  integer, intent(out):: Pole_Singular
  integer :: status,ncFileID,varid
  real,allocatable :: latitudes(:)
  
  Pole_Singular=0
  if(trim(projection)==trim('lon lat')) then
    ! find wether poles are included (or almost included) in grid
    !
    ! If some cells are to narrow (Poles in lat lon coordinates),
    ! this will give too small time steps in the Advection,
    ! because of the constraint that the Courant number should be <1.
    !
    ! If Poles are found and lon-lat coordinates are used the Advection scheme
    ! will be modified to be able to cope with the singularity
    ! the advection routine will not work efficiently with NPROCY>2 in this case
    if(ME_MPI==0)then
      ! open an existing netcdf dataset
      status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
      if(status/=nf90_noerr) then
        print *,'not found',trim(filename)
        call StopAll("GridValues: File not found")
      end if
      
      allocate(latitudes(JJFULLDOM))
      status=nf90_inq_varid(ncid=ncFileID, name="lat", varID=varID)
      if(status/=nf90_noerr) then
        !WRF format
        call check(nf90_inq_varid(ncid=ncFileID, name="XLAT", varID=varID))
        call check(nf90_get_var(ncFileID, varID,latitudes ,start=(/1,1/), count=(/1,JJFULLDOM/)   ))
      else
        call check(nf90_get_var(ncFileID, varID,latitudes  ))
      endif
      
      if(latitudes(RUNDOMAIN(4))>88.0)then
        write(*,*)'The grid is singular at North Pole'
        Pole_Singular=Pole_Singular+1
      end if
      if(latitudes(RUNDOMAIN(3))<-88.0)then
        write(*,*)'The grid is singular at South Pole'
        Pole_Singular=Pole_Singular+1
      end if
      deallocate(latitudes)
      call check(nf90_close(ncFileID))
    end if
  end if
  
  CALL MPI_BCAST(Pole_Singular ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)

end subroutine find_poles

subroutine Getgridparams(LIMAX,LJMAX,filename,cyclicgrid)
  ! Get grid and time parameters as defined in the meteo or Grid_Def file
  ! Do some checks on sizes and dates
  !
  ! This routine is called only once (and is therefore not optimized for speed)

  ! defines i_fdom,j_fdom,i_local,j_local,Cyclicgrid,North_pole,Poles
  ! GRIDWIDTH_M, glon, glat, xm_i,xm_j,xm2,xmd,xmdji,xm2ji,gl_stagg,gb_stagg
  ! P0,A_bnd_met,B_bnd_met,A_bnd,B_bnd,A_mid,B_mid,sigma_mid,sigma_bnd
  ! for Stereographic projection:
  !   ref_latitude,fi,xp,yp,AN
  ! for lon lat projection:
  !   no additional parameters
  ! for Rotated_Spherical projection: 
  !  grid_north_pole_latitude,grid_north_pole_longitude,x1_rot,y1_rot,dx_rot,dx_roti


  implicit none

  integer, intent(in):: LIMAX,LJMAX
  character (len = *), intent(in) ::filename
  integer, intent(out):: cyclicgrid

  integer :: n,i,j,k,kk
  integer :: ncFileID,idimID,jdimID,varID
  integer :: status,South_pole,North_pole
  real :: x1,x2,x3,x4,P0,x,y,mpi_out,r,t
  logical::found_hybrid=.false.,found_metlevels=.false.
  real :: CEN_LAT, CEN_LON,P_TOP_MET, WRF_DY
  real :: rb,rl,rp,dx,dy,dy2,glmax,glmin,v2(2),glon_fdom1,glat_fdom1,lat
  integer :: iloc_start, iloc_end,jloc_start, jloc_end

  real, dimension(-1:LIMAX+2,-1:LJMAX+2)::xm,xm_i_ext,xm_j_ext
  real, dimension(0:LIMAX+1,0:LJMAX+1)::lon_ext,lat_ext


  !define longitudes in interval [-180,180]
  glmin = -180.0
  glmax = glmin + 360.0

  !   we can already define some arrays:
  
  !/ Define full-domain coordinates of local i,j values. We need to account for
  !  the fact that each parallel domain has its starting cordinate
  !  gi0, gj0, and the user may specify a set of lower-left starting
  !  coordinates for running the model, IRUNBEG, JRUNBEG
  !       i_fdom(i)  = i + gi0 + IRUNBEG - 2
  !       j_fdom(j)  = j + gj0 + JRUNBEG - 2
  i_fdom = (/ (n + gi0 + IRUNBEG - 2, n=0,LIMAX+1) /)
  j_fdom = (/ (n + gj0 + JRUNBEG - 2, n=0,LJMAX+1) /)
  
  ! And the reverse, noting that we even define for area
  ! outside local domain
  i_local = (/ (n - gi0 - IRUNBEG + 2, n=1, IIFULLDOM) /)
  j_local = (/ (n - gj0 - JRUNBEG + 2, n=1, JJFULLDOM) /)
  
  call CheckStop(GIMAX+IRUNBEG-1 > IIFULLDOM, "GridRead: I outside domain" )
  call CheckStop(GJMAX+JRUNBEG-1 > JJFULLDOM, "GridRead: J outside domain" )
  
  !open an existing netcdf dataset
  status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
  if(status/=nf90_noerr) then
    print *,'not found',trim(filename)
    call CheckStop("GridValues: File not found")
  end if
  if(MasterProc)print *,'Defining grid parameters from ',trim(filename)
  
  !get dimensions id
  select case(projection)
  case('Stereographic')
    status = nf90_inq_dimid(ncid=ncFileID, name="i", dimID=idimID)
    if(status/=nf90_noerr)& ! WRF format
      call check(nf90_inq_dimid(ncid=ncFileID, name="west_east", dimID=idimID))
    status = nf90_inq_dimid(ncid=ncFileID, name="j", dimID=jdimID)
    if(status/=nf90_noerr)& ! WRF format
      call check(nf90_inq_dimid(ncid=ncFileID, name="south_north", dimID=jdimID))
  case('lambert')
    status = nf90_inq_dimid(ncid=ncFileID, name="x", dimID=idimID)
    if(status/=nf90_noerr)& ! WRF format
      call check(nf90_inq_dimid(ncid=ncFileID, name="west_east", dimID=idimID))
    status = nf90_inq_dimid(ncid=ncFileID, name="y", dimID=jdimID)
    if(status/=nf90_noerr)& ! WRF format
      call check(nf90_inq_dimid(ncid=ncFileID, name="south_north", dimID=jdimID))
  case('lon lat')
    status=nf90_inq_dimid(ncid=ncFileID, name="lon", dimID=idimID)
    if(status/=nf90_noerr)& ! WRF format
      call check(nf90_inq_dimid(ncid=ncFileID, name="west_east", dimID=idimID))
    status=nf90_inq_dimid(ncid=ncFileID, name="lat", dimID=jdimID)
    if(status/=nf90_noerr)& ! WRF  format
      call check(nf90_inq_dimid(ncid=ncFileID, name="south_north", dimID=jdimID))
  case default
    ! write(*,*)'GENERAL PROJECTION ',trim(projection)
    status=nf90_inq_dimid(ncid=ncFileID, name="i", dimID=idimID)
    if(status/=nf90_noerr)& ! WRF format
      call check(nf90_inq_dimid(ncid=ncFileID, name="west_east", dimID=idimID))
    status=nf90_inq_dimid(ncid=ncFileID, name="j", dimID=jdimID)
    if(status/=nf90_noerr)& ! WRF  format
      call check(nf90_inq_dimid(ncid=ncFileID, name="south_north", dimID=jdimID))
  end select
  
  !get global attributes
  status = nf90_get_att(ncFileID,nf90_global,"Grid_resolution",GRIDWIDTH_M)
  if(status/=nf90_noerr)then
    !WRF  format
    call check(nf90_get_att(ncFileID,nf90_global,"DX",GRIDWIDTH_M))
    status = nf90_get_att(ncFileID,nf90_global,"DY",v(1))
    if(status==nf90_noerr .and. abs(GRIDWIDTH_M-v(1))>0.01) then
    ! if(MasterProc)write(*,*)'Gridcells not square. Will correct y mapping factor'
    endif
  end if
  if(MasterProc)write(*,*)"Grid_resolution",GRIDWIDTH_M
  
  select case(projection)
  case('Stereographic')
    status=nf90_get_att(ncFileID,nf90_global,"ref_latitude",ref_latitude)
    if(status/=nf90_noerr)&  ! WRF format
      call check(nf90_get_att(ncFileID,nf90_global,"TRUELAT1",ref_latitude))
    status = nf90_get_att(ncFileID, nf90_global, "fi",fi )
    if(status/=nf90_noerr)& ! WRF format
      call check(nf90_get_att(ncFileID, nf90_global, "STAND_LON",fi ))
    status = nf90_get_att(ncFileID, nf90_global, "xcoordinate_NorthPole",xp )
    if(status==nf90_noerr) then
      call check(nf90_get_att(ncFileID, nf90_global, "ycoordinate_NorthPole",yp ))
    else
      !WRF format compute from grid center coordinates
      call check(nf90_get_att(ncFileID, nf90_global, "CEN_LAT", CEN_LAT))
      call check(nf90_get_att(ncFileID, nf90_global, "CEN_LON", CEN_LON))
      xp = 0.5 + 0.5*IIFULLDOM &
         - EARTH_RADIUS/GRIDWIDTH_M*(1+sin(ref_latitude*PI/180.))&
          *tan(PI/4-CEN_LAT*PI/180./2)*sin((CEN_LON-fi)*PI/180.)
      yp = 0.5 + 0.5*JJFULLDOM &
         + EARTH_RADIUS/GRIDWIDTH_M*(1+sin(ref_latitude*PI/180.))&
          *tan(PI/4-CEN_LAT*PI/180./2)*cos((CEN_LON-fi)*PI/180.)
      !correct for last digits. Assume that numbers are close enough to an integer
      if(abs(nint(xp)-xp)<0.01)xp=nint(xp)
      if(abs(nint(yp)-yp)<0.01)yp=nint(yp)
      
      if(MasterProc)then
        write(*,*)"M= ",EARTH_RADIUS/GRIDWIDTH_M*(1+sin(ref_latitude*PI/180.))
        write(*,*)"coordinates of North pole ",xp,yp
      end if
    end if
    
    AN = 6.370e6*(1.0+sin(ref_latitude*PI/180.))/GRIDWIDTH_M
     ! = 237.7316364 for GRIDWIDTH_M=50 km and ref_latitude=60
    
    do j = 0, LJMAX+1
      dy  = yp - j_fdom(j)
      dy2 = dy*dy
      do i = 0, LIMAX+1
        dx = i_fdom(i) - xp
        rp = sqrt(dx*dx+dy2)           ! => distance to pole
        rb = 90.0 - 180.0/PI*2* atan(rp/AN)  ! => latitude
        rl = 0.0
        if (rp >  1.0e-10) rl = fi + 180.0/PI*atan2(dx,dy)
        if (rl <  glmin)   rl = rl + 360.0
        if (rl >  glmax)   rl = rl - 360.0
        lon_ext(i,j)=rl              ! longitude
        lat_ext(i,j)=rb              ! latitude
        
      end do ! i
    end do ! j
    
  case('lambert')
    if(USE_WRF_MET_NAMES)then
      call check(nf90_get_att(ncFileID, nf90_global, "TRUELAT1",lat_stand1_lambert ))
      call check(nf90_get_att(ncFileID, nf90_global, "TRUELAT2",lat_stand2_lambert ))
      lat0_lambert = lat_stand1_lambert
      call check(nf90_get_att(ncFileID, nf90_global, "STAND_LON",lon0_lambert ))
      earth_radius_lambert = earth_radius
    else
      status = nf90_get_att(ncFileID,nf90_global,"latitude_of_projection_origin",lat0_lambert)!reference latitude, at which x=0
      if(status/=nf90_noerr) then
        call check(nf90_inq_varid(ncid=ncFileID, name="projection_lambert", varID=varID))
        call check(nf90_get_att(ncFileID,varID,"latitude_of_projection_origin",lat0_lambert))
      endif
      status = nf90_get_att(ncFileID,nf90_global,"longitude_of_central_meridian",lon0_lambert)!reference longitude
      if(status/=nf90_noerr) then
        call check(nf90_inq_varid(ncid=ncFileID, name="projection_lambert", varID=varID))
        call check(nf90_get_att(ncFileID,varID,"longitude_of_central_meridian",lon0_lambert))
      endif
      status = nf90_get_att(ncFileID,nf90_global,"earth_radius",earth_radius_lambert)!
      if(status/=nf90_noerr) then
        call check(nf90_inq_varid(ncid=ncFileID, name="projection_lambert", varID=varID))
        call check(nf90_get_att(ncFileID,varID,"earth_radius",earth_radius_lambert))
      endif
      !status = nf90_get_att(ncFileID,nf90_global,"standard_parallel",(/lat_stand1_lambert,lat_stand2_lambert/))!standard latitude at which mapping factor=1
      !default lat_stand1_lambert=lat_stand2_lambert=lat0_lambert
      lat_stand1_lambert = lat0_lambert
      lat_stand2_lambert = lat0_lambert
    endif
    if(abs(lat_stand1_lambert-lat_stand2_lambert)>1.E-4)then
      k_lambert = log(cos(deg2rad*lat_stand1_lambert)/cos(deg2rad*lat_stand2_lambert))/&
      (log(tan(0.25*PI+0.5*deg2rad*lat_stand2_lambert)/tan(0.25*PI+0.5*deg2rad*lat_stand1_lambert)))
      lat0_lambert = rad2deg*asin(k_lambert)
      if(MasterProc)then
        write(*,*)'first true latitude ',lat_stand1_lambert
        write(*,*)'second true latitude ',lat_stand2_lambert
        write(*,*)'latitude of projection origin calculated to ',lat0_lambert
      end if
    else
      k_lambert = sin(deg2rad*lat0_lambert)! also called n
    endif
    F_lambert = cos(deg2rad*lat_stand1_lambert) &
              * tan(0.25*PI+0.5*deg2rad*lat_stand1_lambert)**k_lambert /k_lambert!normalization constant
    y0_lambert = F_lambert*tan(0.25*PI-0.5*deg2rad*lat0_lambert)**k_lambert!reference y coordinate, also called rho0
    if(USE_WRF_MET_NAMES)then
      call check(nf90_inq_varid(ncid=ncFileID, name="XLONG", varID=varID))
      call check(nf90_get_var(ncFileID,varID,v,count=(/1/)))
      call check(nf90_inq_varid(ncid=ncFileID, name="XLAT", varID=varID))
      call check(nf90_get_var(ncFileID,varID,u,count=(/1/)))
    else
      status = nf90_inq_varid(ncid=ncFileID, name="lon", varID=varID)
      if(status/=nf90_noerr)&
        call check(nf90_inq_varid(ncid=ncFileID, name="longitude", varID=varID))
      call check(nf90_get_var(ncFileID,varID,v,count=(/1/)))
      x1_lambert=v(1)
      status = nf90_inq_varid(ncid=ncFileID, name="lat", varID=varID)
      if(status/=nf90_noerr)&
        call check(nf90_inq_varid(ncid=ncFileID, name="latitude", varID=varID))
      call check(nf90_get_var(ncFileID,varID,u,count=(/1/)))
      y1_lambert=u(1)
    endif
    
    x1_lambert=0.0
    y1_lambert=0.0
    call lb2ij(v(1),u(1),x1_lambert,y1_lambert)
    x1_lambert = (x1_lambert-1)*GRIDWIDTH_M
    y1_lambert = (y1_lambert-1)*GRIDWIDTH_M

    if(MasterProc)then
    ! test that lon lat from meteo is the same as calculated, for the point (i,j)=(1,1)
       call lb2ij(v(1),u(1),v(2),u(2))
       if(abs(u(2)-1.0)+abs(v(2)-1.0)>1.E-4)then
          write(*,*)'ERROR Lambert projection lon lat of (i,j)=(1,1) in file is',v(1),u(1)
          write(*,*)'BUT Lambert projection (i,j) for this lon lat is not (1,1) but ',v(2),u(2)
          call StopAll('ERROR in Lambert projection parameters')
       else
          write(*,*)'Lambert projection checked. (i,j)=(1,1) has lon lat',v(1),u(1)
          write(*,*)"x and y at (i,j)=(1,1)",x1_lambert,y1_lambert
          write(*,*)"y0_lambert,F_lambert ",y0_lambert,F_lambert
       endif
    endif
    
    !make lon lat and mapping factors
    do j = 0, LJMAX+1
      y = (y1_lambert+(j_fdom(j)-1)*GRIDWIDTH_M)/EARTH_RADIUS
      do i = 0, LIMAX+1
        x = (x1_lambert+(i_fdom(i)-1)*GRIDWIDTH_M)/EARTH_RADIUS
        r = sqrt(x*x+(y0_lambert-y)*(y0_lambert-y))
        if(k_lambert<0.0)r = -r
        t = atan(x/(y0_lambert-y))
        lat_ext(i,j) = 2*rad2deg*atan((F_lambert/r)**(1.0/k_lambert))-90.0
        lon_ext(i,j) = lon0_lambert + rad2deg*t/k_lambert
        !does not work for lat = -90.0
        xm(i,j)=k_lambert*F_lambert&
          *tan(PI*0.25-deg2rad*0.5*lat_ext(i,j))**(k_lambert-1)&
          *0.5/(cos(PI*0.25-deg2rad*0.5*lat_ext(i,j))**2)
        xm2(i,j) = xm(i,j)*xm(i,j)
        xmd(i,j) = 1.0/xm2(i,j)
        xm2ji(j,i) = xm2(i,j)
        xmdji(j,i) = xmd(i,j)
      enddo
    enddo
    !staggered map factors
    do j = 0, LJMAX+1
      y = (y1_lambert+(j_fdom(j)-1+0.5)*GRIDWIDTH_M)/EARTH_RADIUS
      do i = 0, LIMAX+1
        x = (x1_lambert+(i_fdom(i)-1)*GRIDWIDTH_M)/EARTH_RADIUS
        r = sqrt(x*x+(y0_lambert-y)*(y0_lambert-y))
        if(k_lambert<0.0)r = -r
        lat = 2*rad2deg*atan((F_lambert/r)**(1.0/k_lambert))-90.0
        xm_i(i,j)=k_lambert*F_lambert&
          *tan(PI*0.25-deg2rad*0.5*lat)**(k_lambert-1)&
          *0.5/(cos(PI*0.25-deg2rad*0.5*lat)**2)
      enddo
    enddo
    do j = 0, LJMAX+1
      y = (y1_lambert+(j_fdom(j)-1)*GRIDWIDTH_M)/EARTH_RADIUS
      do i = 0, LIMAX+1
        x = (x1_lambert+(i_fdom(i)-1+0.5)*GRIDWIDTH_M)/EARTH_RADIUS
        r = sqrt(x*x+(y0_lambert-y)*(y0_lambert-y))
        if(k_lambert<0.0)r = -r
        lat = 2*rad2deg*atan((F_lambert/r)**(1.0/k_lambert))-90.0
        xm_j(i,j)=k_lambert*F_lambert&
          *tan(PI*0.25-deg2rad*0.5*lat)**(k_lambert-1)&
          *0.5/(cos(PI*0.25-deg2rad*0.5*lat)**2)
      enddo
    enddo
  case('lon lat')
    if(.not. USE_WRF_MET_NAMES)then
      !NB: lon and lat are stored as 1 dimensional arrays
      call check(nf90_inq_varid(ncid=ncFileID, name="lon", varID=varID))
      
      call check(nf90_get_var(ncFileID, varID, lon_ext(1:limax,1),&
        start=(/gi0+IRUNBEG-1/),count=(/limax/) ))
      if(LIMAX>limax)&
        lon_ext(LIMAX,1)=lon_ext(limax,1)+(lon_ext(limax,1)-lon_ext(limax-1,1))
      lon_ext(0,1)=2*lon_ext(1,1)-lon_ext(2,1)
      lon_ext(LIMAX+1,1)=2*lon_ext(LIMAX,1)-lon_ext(LIMAX-1,1)
      do j=0,LJMAX+1
        lon_ext(:,j)=lon_ext(:,1)
      end do
      
      call check(nf90_inq_varid(ncid=ncFileID,name="lat",varID=varID))
      call check(nf90_get_var(ncFileID, varID, lat_ext(1,1:ljmax),&
        start=(/gj0+JRUNBEG-1/),count=(/ljmax/) ))
      lat_ext(1,LJMAX)=min(90.0,lat_ext(1,LJMAX))!should never be used anyway
      lat_ext(1,0)=2*lat_ext(1,1)-lat_ext(1,2)
      lat_ext(1,LJMAX+1)=2*lat_ext(1,LJMAX)-lat_ext(1,LJMAX-1)
      do i=0,LIMAX+1
        lat_ext(i,:)=lat_ext(1,:)
      end do
    else
      !WRF  format
      call check(nf90_inq_varid(ncid=ncFileID, name="XLONG", varID=varID))
      call nf90_get_var_extended(ncFileID,varID,lon_ext,0,LIMAX+1,0,LJMAX+1)
      call check(nf90_inq_varid(ncid=ncFileID, name="XLAT", varID=varID))
      call nf90_get_var_extended(ncFileID,varID,lat_ext,0,LIMAX+1,0,LJMAX+1)
    end if
  case('Rotated_Spherical')
    status=nf90_get_att(ncFileID,nf90_global,"grid_north_pole_latitude",grid_north_pole_latitude)
    if(status/=nf90_noerr)& ! WRF format
      call check(nf90_get_att(ncFileID,nf90_global,"POLE_LAT",grid_north_pole_latitude))
    if(MasterProc)write(*,*)"grid_north_pole_latitude",grid_north_pole_latitude
    status=nf90_get_att(ncFileID,nf90_global,"grid_north_pole_longitude",grid_north_pole_longitude)
    if(status/=nf90_noerr) then
      ! WRF format
      call check(nf90_get_att(ncFileID,nf90_global,"POLE_LON",grid_north_pole_longitude))
      !find resolution in degrees from resolution in km. WRF uses Erath Radius 6370 km(?)
      dx_rot=360./(6370000.*2*PI/GRIDWIDTH_M)
      !round to 6 digits
      dx_rot=0.000001*nint(1000000*dx_rot)
    end if
    if(MasterProc)write(*,*)"grid_north_pole_longitude",grid_north_pole_longitude
    status=nf90_inq_varid(ncid=ncFileID, name="i", varID=varID)
    if(status==nf90_noerr) then
      call check(nf90_get_var(ncFileID, varID, v2))!note that i is one dimensional
      x1_rot=v2(1)
      dx_rot=v2(2)-v2(1)
      call check(nf90_inq_varid(ncid=ncFileID, name="j", varID=varID))
      call check(nf90_get_var(ncFileID, varID, v2(1)))!note that j is one dimensional
      y1_rot=v2(1)
      call check(nf90_inq_varid(ncid=ncFileID, name="lon", varID=varID))
      call nf90_get_var_extended(ncFileID,varID,lon_ext,0,LIMAX+1,0,LJMAX+1)
      call check(nf90_inq_varid(ncid=ncFileID, name="lat", varID=varID))
      call nf90_get_var_extended(ncFileID,varID,lat_ext,0,LIMAX+1,0,LJMAX+1)
    else
      ! WRF format
      call check(nf90_inq_varid(ncid=ncFileID, name="XLONG", varID=varID))
      call nf90_get_var_extended(ncFileID,varID,lon_ext,0,LIMAX+1,0,LJMAX+1)
      call check(nf90_get_var(ncFileID, varID, v2,start=(/1,1/),count=(/1,1/)  ))
      glon_fdom1=v2(1)
    ! glon=0.0!to get some value for outside subdomain too (when limax<LIMAX for instance)
    ! call check(nf90_get_var(ncFileID, varID, glon(1:limax,1:ljmax),&
    !      start=(/gi0+IRUNBEG-1,gj0+JRUNBEG-1/),count=(/limax,ljmax/)  ))
      call check(nf90_inq_varid(ncid=ncFileID, name="XLAT", varID=varID))
      call nf90_get_var_extended(ncFileID,varID,lat_ext,0,LIMAX+1,0,LJMAX+1)
      call check(nf90_get_var(ncFileID, varID, v2,start=(/1,1/),count=(/1,1/)  ))
      glat_fdom1=v2(1)
      
      x1_rot=0.
      y1_rot=0.
      dx_roti=1.0/dx_rot
      call lb2ij(glon_fdom1,glat_fdom1,x,y)
      x1_rot=(x-1)*dx_rot
      y1_rot=(y-1)*dx_rot
      do i=1,10
        if(x1_rot>180.0)then
          x1_rot=x1_rot-360.0
        else if(x1_rot<-180.0)then
          x1_rot=x1_rot+360.0
        else
          exit
        end if
      end do
    ! call lb2ij(glon_fdom(1,1),glat_fdom(1,1),x,y)
    ! write(*,*)'after ',glon_fdom(1,1),glat_fdom(1,1),x,y
    ! call lb_rot2lb(x,y,x1_rot,y1_rot,grid_north_pole_longitude,grid_north_pole_latitude)
    ! write(*,*)"spherical lon lat of (i,j)=(1,1)",x,y,glon_fdom(1,1),glat_fdom(1,1)
      if(MasterProc)write(*,*)"rotated lon lat of (i,j)=(1,1)",x1_rot,y1_rot
      if(MasterProc)write(*,*)"resolution",dx_rot
    end if
    dx_roti=1.0/dx_rot
    
  case default
    ! other projection?
    call check(nf90_inq_varid(ncid=ncFileID, name="lon", varID=varID))
    call nf90_get_var_extended(ncFileID,varID,lon_ext,0,LIMAX+1,0,LJMAX+1)
    call check(nf90_inq_varid(ncid=ncFileID, name="lat", varID=varID))
    call nf90_get_var_extended(ncFileID,varID,lat_ext,0,LIMAX+1,0,LJMAX+1)
  end select
  
  glon(1:LIMAX,1:LJMAX)=lon_ext(1:LIMAX,1:LJMAX)             ! longitude
  glat(1:LIMAX,1:LJMAX)=lat_ext(1:LIMAX,1:LJMAX)             ! latitude
  do j=1,LJMAX
    do i=1,LIMAX
      if(glon(i,j)>glmax)glon(i,j)=glon(i,j)-360.0
      if(glon(i,j)<glmin)glon(i,j)=glon(i,j)+360.0
    end do
  end do

  if(projection=='lon lat')then
     !we require the longitudes to increase if i increases
     if(glon(1,1)>glon(limax,1))then
        write(*,*)'WARNING: shifting longitudes by 360 deg for processor ',me,&
             ' lon was from ', glon(1,1),'to ',glon(limax,1)
        do j=1,LJMAX
           do i=1,LIMAX
              if(glon(i,j)<glon(1,j))glon(i,j)=glon(i,j)+360.0
           end do
        end do
     endif
  endif

  ! map factors
  status=nf90_inq_varid(ncid=ncFileID, name="map_factor", varID=varID)
  
  if(status==nf90_noerr)then
    ! mapping factor at center of cells is defined  
    call nf90_get_var_extended(ncFileID,varID,xm,-1,LIMAX+2,-1,LJMAX+2)
    
    !make "staggered" map factors and other derived map factors
    do j=0,LJMAX+1
      do i=0,LIMAX+1
        xm_i(i,j)=0.5*(xm(i,j)+xm(i,j+1))
        xm_j(i,j)=0.5*(xm(i,j)+xm(i+1,j))
        xm2(i,j)=xm(i,j)*xm(i,j)
        xmd(i,j) =1.0/xm2(i,j)
        xm2ji(j,i) = xm2(i,j)
        xmdji(j,i) = xmd(i,j)
      end do
    end do
  elseif(trim(projection)=='lambert') then
    ! map factors have been computed already (above)
  else
    !map factor are already staggered
    status=nf90_inq_varid(ncid=ncFileID, name="map_factor_i", varID=varID)
    iloc_start=-1
    if(iloc_start+IRUNBEG+gi0-2<1)iloc_start=1!first cell (in i direction)
    iloc_end=LIMAX+2
    if(iloc_end+IRUNBEG+gi0-2>IIFULLDOM)iloc_end=IIFULLDOM+2-gi0-IRUNBEG!last cell
    jloc_start=-1
    if(jloc_start+JRUNBEG+gj0-2<1)jloc_start=1!first cell (in j direction)
    jloc_end=LJMAX+2
    if(jloc_end+JRUNBEG+gj0-2>JJFULLDOM)jloc_end=JJFULLDOM+2-gj0-JRUNBEG!last cell
    
    if(status==nf90_noerr)then
      call nf90_get_var_extended(ncFileID,varID,xm_i_ext,-1,LIMAX+2,-1,LJMAX+2)
    else
      !WRF  format
      call check(nf90_inq_varid(ncid=ncFileID, name="MAPFAC_VX", varID=varID))
      call nf90_get_var_extended(ncFileID,varID,xm_i_ext,-1,LIMAX+2,-1,LJMAX+2,&
        jshift_in=1) !NB:shift j by 1 since wrf start at bottom face
    end if
    
    status=nf90_inq_varid(ncid=ncFileID, name="map_factor_j", varID=varID)
    if(status==nf90_noerr)then
      call nf90_get_var_extended(ncFileID,varID,xm_j_ext,-1,LIMAX+2,-1,LJMAX+2)
    else
      !WRF  format
      call check(nf90_inq_varid(ncid=ncFileID, name="MAPFAC_UY", varID=varID))
      call nf90_get_var_extended(ncFileID,varID,xm_j_ext,-1,LIMAX+2,-1,LJMAX+2,&
        ishift_in=1)  !NB:shift i by 1 since wrf start at left face
      status = nf90_get_att(ncFileID,nf90_global,"DY",WRF_DY)
      if(status==nf90_noerr) then
        !WRF uses DY/MAPFAC_UY while emep uses GRIDWIDTH_M/xm_j for the y size
        if(abs(GRIDWIDTH_M/WRF_DY-1.0)>1.E-6)then
          if(MasterProc)write(*,*)"rescaling y mapfactors with = ",GRIDWIDTH_M/WRF_DY
          xm_j_ext=xm_j_ext*GRIDWIDTH_M/WRF_DY
        endif
      else
        if(MasterProc)write(*,*)"not rescaling y mapfactors"
      endif
    end if
    
    !define xm2, xm_i and xm_j now
    !Note that xm is inverse length: interpolate 1/xm rather than xm
    !for lon lat projection we do not want xm=0.0 at Poles
    do j=0,LJMAX+1
      do i=0,LIMAX+1
        xm_i(i,j)=max(1.0E-5,xm_i_ext(i,j))
        xm_j(i,j)=max(1.0E-5,xm_j_ext(i,j))
        xm2(i,j) = 4.0*( (xm_i_ext(i,j-1)*xm_i_ext(i,j))/&
        (xm_i_ext(i,j-1)+xm_i_ext(i,j))  )&
        *( (xm_j_ext(i-1,j)*xm_j_ext(i,j))/&
        (xm_j_ext(i-1,j)+xm_j_ext(i,j))  )
        xm2(i,j)=max(1.E-7,xm2(i,j))
        xmd(i,j) =1.0/xm2(i,j)
        xm2ji(j,i) = xm2(i,j)
        xmdji(j,i) = xmd(i,j)
      end do
    end do
    
  end if
  
  status=nf90_inq_varid(ncid=ncFileID, name="k", varID=varID)
  if(status/=nf90_noerr)then
    !always use hybrid coordinates at output, if hybrid in input
    if(.not.USE_EtaCOORDINATES)then
      write(*,*)'WARNING: using hybrid levels even if not asked to! ',trim(filename)
      USE_EtaCOORDINATES=.true.
    end if
    if(MasterProc)write(*,*)'reading met hybrid levels from ',trim(filename)
    !          call check(nf90_inq_varid(ncid=ncFileID, name="hyam", varID=varID))
    !          call check(nf90_get_var(ncFileID, varID, A_mid ))
    !          A_mid=P0*A_mid!different definition in modell and grid_Def
    !          call check(nf90_inq_varid(ncid=ncFileID, name="hybm", varID=varID))
    !          call check(nf90_get_var(ncFileID, varID,B_mid))
    status=nf90_inq_varid(ncid=ncFileID, name="P0", varID=varID)
    if(status/=nf90_noerr)&
      status=nf90_inq_varid(ncid=ncFileID, name="p0", varID=varID)
    if(status/=nf90_noerr)then
      status=nf90_inq_varid(ncid=ncFileID, name="P00", varID=varID) !WRF case
      if(status/=nf90_noerr)then
          call StopAll('Do not know how to define vertical levels')
      else
        ! WRF format
        ! asuming sigma levels ZNW=(P-P_TOP_MET)/(PS-P_TOP_MET)
        ! P = A+B*PS = P_TOP_MET*(1-ZNW) + ZNW*PS
        ! B = ZNW
        ! A = P_TOP_MET*(1-ZNW)
        call check(nf90_get_var(ncFileID, varID, P0 ))
        if(.not.allocated(A_bnd_met))allocate(A_bnd_met(KMAX_MET+1),B_bnd_met(KMAX_MET+1))
        call check(nf90_inq_varid(ncid=ncFileID, name="P_TOP", varID=varID))
        call check(nf90_get_var(ncFileID, varID, P_TOP_MET ))
        call check(nf90_inq_varid(ncid=ncFileID, name="ZNW", varID=varID))
        call check(nf90_get_var(ncFileID, varID, B_bnd_met ))
        if(MET_REVERSE_K)then
          A_bnd_met=B_bnd_met!use A_bnd_met as temporary buffer
          do k=1,KMAX_MET+1
            B_bnd_met(k)=A_bnd_met(KMAX_MET+2-k)
          end do
        end if
        A_bnd_met=P_TOP_MET*(1.-B_bnd_met)
        found_metlevels=.true.
      end if
      if(MET_REVERSE_K)then
        if(MasterProc)write(*,*)"Reversed vertical levels from met, P at level boundaries:"
      else
        if(MasterProc)write(*,*)"Vertical levels from met, P at level boundaries:"
      end if
      do k=1,KMAX_MET+1
        if(MasterProc)write(*,44)k, A_bnd_met(k)+P0*B_bnd_met(k)
      end do
    else
      call check(nf90_get_var(ncFileID, varID, P0 ))
      if(MasterProc)write(*,*)'P0 = ',P0
      if(.not.allocated(A_bnd_met))allocate(A_bnd_met(KMAX_MET+1),B_bnd_met(KMAX_MET+1))
      status=nf90_inq_varid(ncid=ncFileID, name="hyai", varID=varID)
      if(status/=nf90_noerr)then
        call check(nf90_inq_varid(ncid=ncFileID, name="ap", varID=varID))
        call check(nf90_get_var(ncFileID, varID, A_bnd_met, count=(/KMAX_MET/)))!read mid values!
        call check(nf90_inq_varid(ncid=ncFileID, name="b", varID=varID))
        call check(nf90_get_var(ncFileID, varID, B_bnd_met, count=(/KMAX_MET/) )) !read mid values!
        A_bnd_met(KMAX_MET+1)=0.0
        B_bnd_met(KMAX_MET+1)=1.0
        do k=KMAX_MET,1,-1
          A_bnd_met(k)=A_bnd_met(k+1)-2.0*(A_bnd_met(k+1)-A_bnd_met(k))!from mid to bnd values!
          B_bnd_met(k)=B_bnd_met(k+1)-2.0*(B_bnd_met(k+1)-B_bnd_met(k))!from mid to bnd values!
        end do
        
        if(MasterProc)write(*,*)'Metdata pressure at level boundaries:'
        do k=1,KMAX_MET+1
          if(MasterProc)write(*,44)k, A_bnd_met(k)+P0*B_bnd_met(k)
        end do
        found_metlevels=.true.
        
      else
        call check(nf90_get_var(ncFileID, varID, A_bnd_met ))
        A_bnd_met=P0*A_bnd_met!different definition in model and grid_Def
        call check(nf90_inq_varid(ncid=ncFileID, name="hybi", varID=varID))
        call check(nf90_get_var(ncFileID, varID, B_bnd_met ))
        found_metlevels=.true.
      end if
    end if
    if(External_Levels_Def)then
      !model levels defined from external text file
      if(MasterProc)&
      write(*,*)'reading external hybrid levels from ',trim(Vertical_levelsFile),&
        A_bnd_met(kMAX_met+1),B_bnd_met(kMAX_met+1)
      P0=Pref
      do k=1,KMAX_MID+1
        read(IO_TMP,*)kk,A_bnd(k),B_bnd(k)
        if(kk/=k.and.MasterProc)write(*,*)'WARNING: unexpected format for vertical levels ',k,kk
      end do
      
      if(.not.found_metlevels)then
        ! assume levels from metdata are defined in Vertical_levelsFile
        if(.not.allocated(A_bnd_met))allocate(A_bnd_met(KMAX_MET+1),B_bnd_met(KMAX_MET+1))
        A_bnd_met=A_bnd
        B_bnd_met=B_bnd
      end if
      
    else
      !vertical model levels are the same as in meteo
      A_bnd=A_bnd_met
      B_bnd=B_bnd_met
    end if
    
    do k=1,KMAX_MID
      A_mid(k)=0.5*(A_bnd(k)+A_bnd(k+1))
      B_mid(k)=0.5*(B_bnd(k)+B_bnd(k+1))
    end do
    sigma_mid =B_mid!for Hybrid coordinates sigma_mid=B if A*P0=PT-sigma_mid*PT
    
    if(MasterProc)write(*,*)"Model pressure at level boundaries:"
    do k=1,KMAX_MID+1
44    FORMAT(i4,10F12.2)
      if(MasterProc)write(*,44)k, A_bnd(k)+P0*B_bnd(k)
    end do
    !test if the top is within the height defined in the meteo files
    if(MasterProc.and.External_Levels_Def.and.(A_bnd(1)+P0*B_bnd(1)+0.01<A_bnd_met(1)+P0*B_bnd_met(1)))then
      write(*,*)'Pressure at top of defined levels is ',A_bnd(1)+P0*B_bnd(1)
      write(*,*)'Pressure at top defined in meteo files is ',A_bnd_met(1)+P0*B_bnd_met(1)
      write(*,*)'Pressure at op must be higher (lower altitude) than top defined in meteo '
      call StopAll('Top level too high! Change values in '//trim(Vertical_levelsFile))
    end if
    
    !test if the levels can cope with highest mountains (400 hPa)
    do k=1,KMAX_MID
      if(MasterProc.and.A_bnd(k+1)+40000*B_bnd(k+1)-(A_bnd(k)+40000*B_bnd(k))<0.0)then
        write(*,*)'WARNING: hybrid vertical level definition may cause negative level thickness when pressure below 400 hPa '
        write(*,*)'Pressure at level ',k,' is ',A_bnd(k)+40000*B_bnd(k)
        write(*,*)'Pressure at level ',k+1,' is ',A_bnd(k+1)+40000*B_bnd(k+1),' (should be higher)'
        if(External_Levels_Def)call StopAll('GridValues_mod: possible negative level thickness ')
      end if
    end do
    
    !test if the lowest levels is thick enough (twice height of highest vegetation?) about 550 Pa = about 46m
    !Deposition scheme is not designes for very thin lowest levels
    if(MasterProc.and.A_bnd(KMAX_MID+1)+P0*B_bnd(KMAX_MID+1)-(A_bnd(KMAX_MID)+P0*B_bnd(KMAX_MID))<550.0)then
      write(*,*)'WARNING: lowest level very shallow; ',A_bnd(KMAX_MID+1)+P0*B_bnd(KMAX_MID+1) -&
      (A_bnd(KMAX_MID)+P0*B_bnd(KMAX_MID)),'Pa'
!      call StopAll('Lowest level too thin! Change vertical levels definition in '//trim(Vertical_levelsFile))
    end if
    
    found_hybrid=.true.
  else
    call check(nf90_get_var(ncFileID, varID, sigma_mid ))
  end if
  call check(nf90_close(ncFileID))
  
  !TEMPORARY: definition of sigma. Only A and B will be used in the future
  ! definition of the half-sigma levels (boundaries between layers)
  ! from the full levels.
  sigma_bnd(KMAX_BND) = 1.
  do k = KMAX_MID,2,-1
    sigma_bnd(k) = 2.*sigma_mid(k) - sigma_bnd(k+1)
  end do
  sigma_bnd(1) = 0.
  
  if(.not.(found_hybrid.or.Grid_Def_exist))then
    !define A and B that gives the correspondin sigma
    do k = 1,KMAX_BND
      A_bnd(k)=PT * (1-sigma_bnd(k))
      B_bnd(k)=sigma_bnd(k)
    end do
    do k = 1,KMAX_MID
      A_mid(k)=(A_bnd(k+1)+A_bnd(k))/2.0
      B_mid(k)=(B_bnd(k+1)+B_bnd(k))/2.0
    end do
  end if
  do k = 1,KMAX_MID
    dA(k)=A_bnd(k+1)-A_bnd(k)
    dB(k)=B_bnd(k+1)-B_bnd(k)
    Eta_bnd(k)=A_bnd(k)/Pref+B_bnd(k)
    Eta_mid(k)=A_mid(k)/Pref+B_mid(k)
    dEta_i(k)=1.0/(dA(k)/Pref+dB(k))
  end do
  Eta_bnd(KMAX_MID+1)=A_bnd(KMAX_MID+1)/Pref+B_bnd(KMAX_MID+1)

  if(External_Levels_Def)call make_vertical_levels_interpolation_coeff
  
  do j=0,LJMAX
    do i=0,LIMAX
      x1=lon_ext(i,j)
      x2=lon_ext(i+1,j)
      x3=lon_ext(i,j+1)
      x4=lon_ext(i+1,j+1)
      
      !8100=90*90; could use any number much larger than zero and much smaller than 180*180
      if(x1*x2<-8100.0 .or. x1*x3<-8100.0 .or. x1*x4<-8100.0)then
        !Points are on both sides of the longitude -180=180
        if(x1<0)x1=x1+360.0
        if(x2<0)x2=x2+360.0
        if(x3<0)x3=x3+360.0
        if(x4<0)x4=x4+360.0
      end if
      gl_stagg(i,j)=0.25*(x1+x2+x3+x4)
      
      gb_stagg(i,j)=0.25*(lat_ext(i,j)+&
      lat_ext(i+1,j)+&
      lat_ext(i,j+1)+&
      lat_ext(i+1,j+1))
   end do
  end do
  
  !ensure that lon values are within [-180,+180]]
  do j=0,LJMAX
    do i=0,LIMAX
      if(gl_stagg(i,j)>180.0)gl_stagg(i,j)=gl_stagg(i,j)-360.0
      if(gl_stagg(i,j)<-180.0)gl_stagg(i,j)=gl_stagg(i,j)+360.0
    end do
  end do
  
  !test if the grid is cyclicgrid:
  !The last cell + 1 cell = first cell
  Cyclicgrid=1 !Cyclicgrid
  do j=1,ljmax
    if(mod(nint(10*(360+GIMAX*(glon(2,j)-glon(1,j)))),3600)/=0)then
      Cyclicgrid=0  !not cyclicgrid
    end if
  end do
  
  if(MasterProc .and. DEBUG%GRIDVALUES)write(*,*)'CYCLICGRID:',Cyclicgrid
  
! Look for poles
! If the northernmost or southernmost lines are poles, they are not considered
! as outer boundaries and will not be treated by "BoundaryConditions_mod".
! If the projection is not lat lon (i.e. the poles are not lines, but points),
! the poles are not a problem and Pole=0, even if the grid actually include a pole.
! Note that "Poles" is defined in subdomains
  
  North_pole=1
  do i=1,limax
    if(nint(glat(i,ljmax))<=88)then
      North_pole=0  !not north pole
    end if
  end do
  
  South_pole=1
  do i=1,limax
    if(nint(glat(i,1))>=-88)then
      South_pole=0  !not south pole
    end if
  end do
  
  Poles=0
  if(North_pole==1)then
    Poles(1)=1
    if(MasterProc .or. me==NPROC-1)write(*,*)me,'Found North Pole'
  end if
  
  if(South_pole==1)then
    Poles(2)=1
    if(MasterProc .or. me==NPROC-1)write(*,*)me,'Found South Pole'
  end if
  do j=1,LJMAX
    do i=1,LIMAX
      GridArea_m2(i,j) = GRIDWIDTH_M*GRIDWIDTH_M*xmd(i,j)
    end do
  end do
  
  if(me_calc>=0)then
    gbacmax = maxval(glat(:,:))
    gbacmin = minval(glat(:,:))
    glacmax = maxval(glon(:,:))
    glacmin = minval(glon(:,:))
  else
    gbacmax = -999.9
    gbacmin =  999.9
    glacmax = -999.9
    glacmin =  999.9
  end if
  
  CALL MPI_ALLREDUCE(gbacmax, mpi_out, 1,MPI_DOUBLE_PRECISION, &
  MPI_MAX, MPI_COMM_WORLD, IERROR)
  gbacmax=mpi_out
  CALL MPI_ALLREDUCE(gbacmin, mpi_out, 1, &
  MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, IERROR)
  gbacmin=mpi_out
  CALL MPI_ALLREDUCE(glacmax, mpi_out, 1,MPI_DOUBLE_PRECISION, &
  MPI_MAX, MPI_COMM_WORLD, IERROR)
  glacmax=mpi_out
  CALL MPI_ALLREDUCE(glacmin, mpi_out, 1, &
  MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, IERROR)
  glacmin=mpi_out
  if(MasterProc) write(unit=6,fmt="(a,40f9.2)") &
  " GridValues: max/min for lat,lon ", &
  gbacmax,gbacmin,glacmax,glacmin

end subroutine Getgridparams

subroutine DefDebugProc()
  ! Find debug coords  and processor

  integer :: i, j
  
  debug_proc = .false.
  
  do i = li0, li1
    do j = lj0, lj1
      if( i_fdom(i)==DEBUG%IJ(1) .and. j_fdom(j)==DEBUG%IJ(2) ) then
        debug_li = i
        debug_lj = j
        debug_proc = .true.
      end if
    end do
  end do
  
  if(debug_proc) write(*,*) "GridValues debug_proc found:", &
  me, debug_li, debug_lj
  if(DEBUG%GRIDVALUES) then
    if(MasterProc) write(*,"(a,2a4,a3,4a4,a2,2a4,4a12)") "GridValues debug:", &
    "D_i", "D_j", "me", "li0", "li1", "lj0", "lj1", &
    "dp" , "d_li", "d_lj", "i_fdom(li0)","i_fdom(li1)", &
    "j_fdom(lj0)", "j_fdom(lj1)"
    
    write(*,"(a,2i4,i3,4i4,L2,2i4,4i12)") "GridValues debug:", &
    DEBUG%IJ(1), DEBUG%IJ(2), me, li0, li1, lj0, lj1, &
    debug_proc , debug_li, debug_lj, &
    i_fdom(li0),i_fdom(li1), j_fdom(lj0), j_fdom(lj1)
  end if

end subroutine DefDebugProc

subroutine lb2ijm(imax,jmax,lon,lat,xr2,yr2,fi2,an2,xp2,yp2)
  !-------------------------------------------------------------------!
  !   calculates coordinates xr2, yr2 (real values) from lat and lon
  !
  !   input: glon,glat:   coord. of the polar point in grid1
  !          an2:   number of grid-distances from pole to equator in grid2.
  !          fi2:      rotational angle for the grid2 (at i2=0).
  !          i1max,j1max: number of points (grid1) in  x- og y- direction
  !
  !
  !   output: xr2(i1,j1): i coordinates in grid2 (with decimals)
  !           yr2(i1,j1): j coordinates in grid2 (with decimals)
  !-------------------------------------------------------------------!

  integer, intent(in) :: imax,jmax
  real, intent(in)    :: lon(imax,jmax),lat(imax,jmax)
  real, intent(out)   :: xr2(imax,jmax),yr2(imax,jmax)
  real, intent(in), optional    :: fi2,an2,xp2,yp2
  real  :: fi_loc,an_loc,xp_loc,yp_loc
  real, parameter :: PI=3.14159265358979323
  real    :: PId4,dr,dr2,dist,dist2,dist3
  integer ::i,j,ip1,jp1, ir2, jr2,i1,j1

  if(projection=='Stereographic'.or.(present(fi2).and.present(an2).and.present(xp2).and.present(yp2)))then
    PId4  =PI/4.
    dr2   =PI/180.0/2.   ! degrees to radians /2
    dr    =PI/180.0      ! degrees to radians
    fi_loc=fi
    an_loc=an
    xp_loc=xp
    yp_loc=yp
    
    if(present(fi2))fi_loc=fi2
    if(present(an2))an_loc=an2
    if(present(xp2))xp_loc=xp2
    if(present(yp2))yp_loc=yp2
    do j1 = 1, jmax
      do i1 = 1, imax
        xr2(i1,j1)=xp_loc+an_loc*tan(PId4-lat(i1,j1)*dr2)*sin(dr*(lon(i1,j1)-fi_loc))
        yr2(i1,j1)=yp_loc-an_loc*tan(PId4-lat(i1,j1)*dr2)*cos(dr*(lon(i1,j1)-fi_loc))
      end do
    end do
  else if(projection=='lon lat')then! lon-lat grid
    do j1 = 1, jmax
      do i1 = 1, imax
        xr2(i1,j1)=(lon(i1,j1)-glon(1,1))/(glon(2,1)-glon(1,1))+i_fdom(1)
        if(xr2(i1,j1)<0.5)xr2=xr2+360.0/(glon(2,1)-glon(1,1))
        yr2(i1,j1)=(lat(i1,j1)-glat(1,1))/(glat(1,2)-glat(1,1))+j_fdom(1)
      end do
    end do
  else ! general projection, Use only info from glon_fdom and glat_fdom
    call StopAll('lb2ijm: not implemented yet')
    !NB: glon_fdom is no more defined. Could easily rewrite if necessary
    dist2=0.0
    dist3=0.0
    
    !VERY SLOW, specially for large grids
    do j1 = 1, jmax
      do i1 = 1, imax
        dist=10.0!max distance is PI
        do j=1,JJFULLDOM
          do i=1,IIFULLDOM
         !  if(dist>great_circle_distance(glon(i1,j1),glat(i1,j1),&
         !            glon_fdom(i,j) ,glat_fdom(i,j)))then
         !    dist=great_circle_distance(glon(i1,j1),glat(i1,j1), &
         !            glon_fdom(i,j),glat_fdom(i,j))
         !    xr2(i1,j1)=i
         !    yr2(i1,j1)=j
         !  end if
          end do
        end do
        
        !find the real part of i and j by comparing distances to neighbouring cells
        !
        !     C
        !    /|\
        !   / | \
        !  /  |  \
        ! A---D---B
        !
        !A=(i,j) ,B=(i+1,j), C=(glon,glat)
        !dist=AC, dist2=BC, dist3=AB
        !AD=(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3)
        !
        ir2 = nint(xr2(i1,j1))
        jr2 = nint(yr2(i1,j1))
        ip1=ir2+1
        if(ip1>IIFULLDOM)ip1=ip1-2
        !             dist2=great_circle_distance(glon(i1,j1),glat(i1,j1),glon_fdom(ip1,jr2),glat_fdom(ip1,jr2))
        !             dist3=great_circle_distance( glon_fdom(ir2,jr2), &
        !                  glat_fdom(ir2,jr2), &
        !                  glon_fdom(ip1,jr2), &
        !                  glat_fdom(ip1,jr2))
        
        xr2(i1,j1)=xr2(i1,j1)+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)
        
        
        jp1=jr2+1
        if(jp1>JJFULLDOM)jp1=jp1-2
        
        !             dist2=great_circle_distance(glon(i1,j1),glat(i1,j1),glon_fdom(ir2,jp1),glat_fdom(ir2,jp1))
        !GFORTRAN CHANGE
        !             dist3=great_circle_distance( glon_fdom(ir2,jr2), &
        !                  glat_fdom(ir2,jr2), &
        !                  glon_fdom(ir2,jp1), &
        !                  glat_fdom(ir2,jp1) )
        
        yr2(i1,j1)=yr2(i1,j1)+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)
        
      end do
    end do
    
  end if
end subroutine lb2ijm

subroutine lb2ij_real(gl,gb,xr2,yr2,fi2,an2,xp2,yp2)
  ! Note: this routine is not yet CPU optimized
  !-------------------------------------------------------------------!
  ! calculates coordinates xr2, yr2 (real values) from gl(lon),gb(lat)
  ! NB: xr2, yr2 are given in FULLDOMAIN coordinates
  !
  ! input:  xp2,yp2:   coord. of the polar point in grid2
  !         an2:   number of grid-distances from pole to equator in grid2.
  !         fi2:      rotational angle for the grid2 (at i2=0).
  !         i1max,j1max: number of points (grid1) in  x- og y- direction
  !
  ! output: xr2(i1,j1): i coordinates in grid2
  !         yr2(i1,j1): j coordinates in grid2
  !-------------------------------------------------------------------!
  real, intent(in)    :: gl,gb
  real, intent(out)   :: xr2,yr2
  real, intent(in), optional    :: fi2,an2,xp2,yp2

  real  :: fi_loc,an_loc,xp_loc,yp_loc
  real, parameter :: PI=3.14159265358979323,PId4=PI/4.0,dr=PI/180.0,dri= 180.0/PI,dr2=dr*0.5
  real    :: dist,dist2,dist3,r
  integer ::i,j,ip1,jp1, ir2, jr2
  real ::xscen ,yscen,zsycen,zcycen ,zxmxc,zsxmxc,zcxmxc,zsysph,zsyrot,yrot,zsxrot,zcysph,zcyrot,zcxrot,xrot
!  real,save :: r_save,lat_save=-999.0

  select case (projection)
  case('Stereographic')
    fi_loc=fi
    an_loc=an
    xp_loc=xp
    yp_loc=yp

    if(present(fi2))fi_loc=fi2
    if(present(an2))an_loc=an2
    if(present(xp2))xp_loc=xp2
    if(present(yp2))yp_loc=yp2

    xr2=xp_loc+an_loc*tan(0.25*PI-0.5*deg2rad*gb)*sin(deg2rad*(gl-fi_loc))
    yr2=yp_loc-an_loc*tan(0.25*PI-0.5*deg2rad*gb)*cos(deg2rad*(gl-fi_loc))

  case('lon lat')           ! lon-lat grid
    if((gl-glon(1,1))+i_fdom(1)*(glon(2,1)-glon(1,1))<360.0)then
       xr2=(gl-glon(1,1))/(glon(2,1)-glon(1,1))+i_fdom(1)
    else
       xr2=(gl-360.0-glon(1,1))/(glon(2,1)-glon(1,1))+i_fdom(1)
    end if
    if(xr2<0.5)xr2=xr2+360.0/(glon(2,1)-glon(1,1))
    yr2=(gb-glat(1,1))/(glat(1,2)-glat(1,1))+j_fdom(1)

  case('Rotated_Spherical') ! rotated lon-lat grid
  ! dx_roti=20.0
  ! grid_north_pole_longitude = -170.0
  ! grid_north_pole_latitude = 40.0
    xscen = grid_north_pole_longitude-180.0
    if(xscen<-180.0)xscen = xscen+360.0
    yscen = 90.0-grid_north_pole_latitude
  ! xscen=grid_north_pole_longitude-180.0
  ! yscen=90.0-grid_north_pole_latitude
    zsycen = sin(dr*yscen)
    zcycen = cos(dr*yscen)
    zxmxc  = dr*(gl - xscen)
    zsxmxc = sin(zxmxc)
    zcxmxc = cos(zxmxc)
    zsysph = sin(dr*gb)
    zcysph = cos(dr*gb)
    zsyrot = zcycen*zsysph - zsycen*zcysph*zcxmxc
    zsyrot = amax1(zsyrot,-1.0)
    zsyrot = amin1(zsyrot,+1.0)
    yrot = asin(zsyrot)
    zcyrot = cos(yrot)
    zcxrot = (zcycen*zcysph*zcxmxc + zsycen*zsysph)/zcyrot
    zcxrot = amax1(zcxrot,-1.0)
    zcxrot = amin1(zcxrot,+1.0)
    zsxrot = zcysph*zsxmxc/zcyrot
    xrot = acos(zcxrot)
    if (zsxrot.lt.0.0) xrot = -xrot
    xrot=xrot*dri
    yrot=yrot*dri
    if(xrot<x1_rot)xrot=xrot+360.0
    if(xrot-x1_rot>360.0-dx_rot*0.499999999)xrot=xrot-360.0
    xr2=(xrot-x1_rot)*dx_roti+1
    yr2=(yrot-y1_rot)*dx_roti+1

  case('lambert')           ! lambert projection
  !  if(gb==lat_save)then
  !     r=r_save
  !  else
      ! r depends only on latitude -> reuse save a little, but not worth it?
        r = F_lambert*tan(PId4-dr2*gb)**k_lambert
  !    r_save=r
  !    lat_save=gb
  ! endif
     xr2 = r*sin(dr*k_lambert*(gl-lon0_lambert))
     yr2 = y0_lambert - r*cos(dr*k_lambert*(gl-lon0_lambert))

    !convert from x,y (erath radius=1) to i,j
     xr2=(xr2*EARTH_RADIUS-x1_lambert)/GRIDWIDTH_M + 1
     yr2=(yr2*EARTH_RADIUS-y1_lambert)/GRIDWIDTH_M + 1

  case default ! general projection, Use only info from glon_fdom and glat_fdom
    !first find closest by testing all gridcells.
    call StopAll('lb2ij: conversion broken 27 Oct 2015, Peter')
    !glon_fdom is no more defined. Could easily rewrite if necessary
    dist2=0.0
    dist3=0.0
    dist=10.0!max distance is PI
    do j=1,JJFULLDOM
      do i=1,IIFULLDOM
!       if(dist>great_circle_distance(gl,gb,glon_fdom(i,j) &
!           ,glat_fdom(i,j)))then
!         dist=great_circle_distance(gl,gb,glon_fdom(i,j) &
!              ,glat_fdom(i,j))
!         xr2=i
!         yr2=j
!       end if
      end do
    end do

    !find the real part of i and j by comparing distances to neighbouring cells
    !
    !     C
    !    /|\
    !   / | \
    !  /  |  \
    ! A---D---B
    !
    !A=(i,j) ,B=(i+1,j), C=(gl,gb)
    !dist=AC, dist2=BC, dist3=AB
    !AD=(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3)
    !
    ir2 = nint(xr2)
    jr2 = nint(yr2)
    ip1=ir2+1
    if(ip1>IIFULLDOM)ip1=ip1-2
!   dist2=great_circle_distance(gl,gb,glon_fdom(ip1,jr2),glat_fdom(ip1,jr2))
!   dist3=great_circle_distance(glon_fdom(ir2,jr2),glat_fdom(ir2,jr2), &
!           glon_fdom(ip1,jr2),glat_fdom(ip1,jr2))

    xr2=xr2+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)


    jp1=jr2+1
    if(jp1>JJFULLDOM)jp1=jp1-2

!   dist2=great_circle_distance(gl,gb,glon_fdom(ir2,jp1),glat_fdom(ir2,jp1))
!   dist3=great_circle_distance(glon_fdom(ir2,jr2),glat_fdom(ir2,jr2), &
!           glon_fdom(ir2,jp1),glat_fdom(ir2,jp1) )

    yr2=yr2+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)

  end select
end subroutine lb2ij_real
subroutine lb2ij_int(gl,gb,ix,iy)
  real, intent(in)    :: gl,gb !gl=lon, gb=lat
  integer, intent(out):: ix,iy
  real ::x,y
! stations can easily be defined exactly at gridcell boundaries
! 1.0E-7 is to ensure same rounding for all situations
  call lb2ij_real(gl,gb,x,y)
  ix=nint(x+1.0E-7)
  iy=nint(y+1.0E-7)
end subroutine lb2ij_int

subroutine ij2lbm(imax,jmax,glon,glat,fi,an,xp,yp)
  !-------------------------------------------------------------------!
  !  calculates lon and lat (geographical coord.)
  !  in every grid point for a polarsteraographic projection.
  !
  ! input:  xp,yp:   coord. of the polar point.
  !         an:      number of grid-distances from pole to equator.
  !         fi:      rotational angle for the x,y grid (at i=0).
  !         imax,jmax:   number of points in  x- og y- direction
  !         glmin:   gives min.value of geographical lenght
  !                  =>  glmin <= l <= glmin+360.
  !                      (example glmin = -180. or 0.)
  !                  if "geopos","georek" is used
  !                  then glmin must be the lenght i(1,1) in the
  !                  geographical grid (gl1 to "geopos")
  ! output: gl(ii,jj): longitude glmin <= l <= glmin+360.
  !         gb(ii,jj): latitude  -90. <= b <= +90.
  !-------------------------------------------------------------------!

  integer :: i, j, imax, jmax
  real    :: glon(imax,jmax),glat(imax,jmax)
  real    :: fi, an, xp, yp
  real    :: om, om2, glmin, glmax,dy, dy2,rp,rb, rl, dx, dr
  real, parameter :: PI=3.14159265358979323

  glmin = -180.0
  
  glmax = glmin + 360.0
  dr    = PI/180.0      ! degrees to radians
  om    = 180.0/PI      ! radians to degrees (om=Norwegian omvendt?)
  om2   = om * 2.0
  
  do j = 1, jmax
    dy  = yp - j
    dy2 = dy*dy
    do i = 1, imax
      
      dx = i - xp    ! ds - changed
      rp = sqrt(dx*dx+dy2)           ! => distance to pole
      rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
      rl = 0.0
      if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
      if (rl <  glmin)   rl = rl + 360.0
      if (rl >  glmax)   rl = rl - 360.0
      glon(i,j)=rl                   ! longitude
      glat(i,j)=rb                   ! latitude
    end do ! i
  end do ! j
  
end subroutine ij2lbm

subroutine ij2lb(i,j,lon,lat,fi,an,xp,yp)
  !-------------------------------------------------------------------!
  !  calculates lon and lat (geographical coord.)
  !  from i,j coordinates in polar stereographic projection
  !
  !  input:  i,j
  !          xp,yp:   coord. of the polar point.
  !          an:      number of grid-distances from pole to equator.
  !          fi:      rotational angle for the x,y grid (at i=0).
  !          imax,jmax:   number of points in  x- og y- direction
  !          glmin:   gives min.value of geographical lenght
  !                   =>  glmin <= l <= glmin+360.
  !                       (example glmin = -180. or 0.)
  !                   if "geopos","georek" is used
  !                   then glmin must be the lenght i(1,1) in the
  !                   geographical grid (gl1 to "geopos")
  !  output: lon: longitude glmin <= lon <= glmin+360.
  !          lat: latitude  -90. <= lat <= +90.
  !-------------------------------------------------------------------!

  integer :: i, j
  real    :: lon,lat
  real    :: fi, an, xp, yp
  real    :: om, om2, glmin, glmax,dy, dy2,rp,rb, rl, dx, dr
  real, parameter :: PI=3.14159265358979323
  
  glmin = -180.0
  
  glmax = glmin + 360.0
  dr    = PI/180.0      ! degrees to radians
  om    = 180.0/PI      ! radians to degrees (om=Norwegian omvendt?)
  om2   = om * 2.0
  
  dy  = yp - j
  dy2 = dy*dy
  
  dx = i - xp    ! ds - changed
  rp = sqrt(dx*dx+dy2)           ! => distance to pole
  rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
  rl = 0.0
  if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
  if (rl <  glmin)   rl = rl + 360.0
  if (rl >  glmax)   rl = rl - 360.0
  lon=rl                   ! longitude
  lat=rb                   ! latitude
end subroutine ij2lb

subroutine ij2ijm(in_field,imaxin,jmaxin,out_field,imaxout,jmaxout, &
     fiin,anin,xpin,ypin,fiout,anout,xpout,ypout)
  !   Converts data (in_field) stored in polar stereo coordinates
  !   with parameters "fiin,anin,xpin,ypin,"
  !   into data (out_field) in polar stereo coordinates with parameters
  !   "fiout,anout,xpout,ypout"

  integer, intent(in) :: imaxin,jmaxin,imaxout,jmaxout
  real, intent(in) :: fiin,anin,xpin,ypin,fiout,anout,xpout,ypout
  real, intent(in) :: in_field(imaxin,jmaxin)! Field to be transformed
  real, intent(out) :: out_field(imaxout,jmaxout)! Field to be transformed
  
  real, allocatable,dimension(:,:) :: x,y,glat,glon
  integer alloc_err,i,j,i2,j2
  logical :: interpolate
  real :: f11,f12,f21,f22
  
  interpolate = .true.
  
  allocate(x(imaxout,jmaxout), stat=alloc_err)
  allocate(y(imaxout,jmaxout), stat=alloc_err)
  allocate(glat(imaxout,jmaxout), stat=alloc_err)
  allocate(glon(imaxout,jmaxout), stat=alloc_err)
  if ( alloc_err/=0 ) WRITE(*,*) 'MPI_ABORT: ', "ij2ij alloc failed"
  if ( alloc_err/=0 ) call  MPI_ABORT(MPI_COMM_CALC,9,IERROR)
  
  ! find longitude, latitude of wanted area
  call ij2lbm(imaxout,jmaxout,glon,glat,fiout,anout,xpout,ypout)

  ! find corresponding coordinates (i,j) in in_field coordinates
  call lb2ijm(imaxout,jmaxout,glon,glat,x,y,fiin,anin,xpin,ypin)


  ! check if the corners of the domain are inside the area covered by the
  ! in_grid: (In principle we should test for all i,j , but test the corners
  ! should be good enough in practice)

  if(int(x(1,1)) < 1 .or. int(x(1,1))+1 > imaxin .or. &
     int(x(imaxout,1)) < 1 .or. int(x(imaxout,1))+1 > imaxin .or. &
     int(x(1,jmaxout)) < 1 .or. int(x(1,jmaxout))+1 > imaxin .or. &
     int(x(imaxout,jmaxout)) < 1 .or. &
     int(x(imaxout,jmaxout))+1 > imaxin .or. &
     int(y(1,1)) < 1 .or. int(y(1,1))+1 > jmaxin .or. &
     int(y(imaxout,1)) < 1 .or. int(y(imaxout,1))+1 > jmaxin .or. &
     int(y(1,jmaxout)) < 1 .or. int(y(1,jmaxout))+1 > jmaxin .or. &
     int(y(imaxout,jmaxout)) < 1 .or. &
     int(y(imaxout,jmaxout))+1 > jmaxin ) then
    write(*,*)'Did not find all the necessary data in in_field'
    write(*,*)'values needed: '
    write(*,*)x(1,1),y(1,1)
    write(*,*)x(imaxout,1),y(imaxout,1)
    write(*,*)x(1,jmaxout),y(1,jmaxout)
    write(*,*)x(imaxout,jmaxout),y(imaxout,jmaxout)
    write(*,*)'max values found: ',imaxin ,jmaxin
    write(*,*) 'MPI_ABORT: ', "ij2ij: area to small"
    call MPI_ABORT(MPI_COMM_CALC,9,IERROR)
  end if


  !  interpolate fields if required

  if(interpolate)then
    do j = 1, jmaxout
      do i = 1,imaxout
        i2=int(x(i,j))
        j2=int(y(i,j))
        f11=(1.-(x(i,j)-i2))*(1.-(y(i,j)-j2))
        f12=(1.-(x(i,j)-i2))*((y(i,j)-j2))
        f21=((x(i,j)-i2))*(1.-(y(i,j)-j2))
        f22=((x(i,j)-i2))*((y(i,j)-j2))
        
        out_field(i,j) =  &
        f11 * in_field(i2,j2) +  &
        f12 * in_field(i2,j2+1) +  &
        f21 * in_field(i2+1,j2) +  &
        f22 * in_field(i2+1,j2+1)
        
      end do
    end do
  else
    do j = 1, jmaxout
      do i = 1,imaxout
        out_field(i,j) =in_field(nint(x(i,j)),nint(y(i,j)))
      end do
    end do
  end if
  
  deallocate(x,stat=alloc_err)
  deallocate(y,stat=alloc_err)
  deallocate(glat,stat=alloc_err)
  deallocate(glon,stat=alloc_err)
  if(alloc_err/=0)then
    WRITE(*,*) 'MPI_ABORT: ', "ij2ijde-alloc_err"
    call  MPI_ABORT(MPI_COMM_CALC,9,IERROR)
  end if

end subroutine ij2ijm

subroutine range_check(vname,var,vrange,fatal)
  character(len=*), intent(in) :: vname
  real, intent(in) :: var,vrange(0:1)
  logical :: fatal
  character(len=*), parameter :: &
       errfmt="(A,'=',F6.2,' is out of range ',F6.2,'..',F6.2)"
  character(len=len_trim(vname)+21+6*3) :: errmsg
  if(var<vrange(0).or.var>vrange(1))then
     write(errmsg,errfmt)trim(vname),var,vrange
     if(fatal)then
        call CheckStop("range_check",trim(errmsg))
     else
        write(*,*)"WARNING: ",trim(errmsg)
     end if
  end if
end subroutine range_check
subroutine coord_check(msg,lon,lat,fix)
  !-------------------------------------------------------------------!
  ! lon/lat range check.
  !   Some longitude range errors can be corrected, when fix=.true.
  !   Latitude range errors are always fatal.
  !-------------------------------------------------------------------!
  character(len=*), intent(in) :: msg
  real, intent(inout) :: lon,lat
  logical :: fix
  call range_check(trim(msg)//" lat",lat,(/ -90.0, 90.0/),fatal=.true.)
  call range_check(trim(msg)//" lon",lon,(/-180.0,180.0/),fatal=.not.fix)
  if(fix)then
     lon=modulo(lon+180.0,360.0)-180.0 ! lon/gl_stagg range -180 .. 180
     call range_check(trim(msg)//" lon",lon,(/-180.0,180.0/),fatal=.true.)
  end if
end subroutine coord_check
function coord_in_domain(domain,lon,lat,iloc,jloc,iglob,jglob) result(in)
  !-------------------------------------------------------------------!
  ! Is coord (lon/lat) is inside global domain|local domain|grid cell?
  !-------------------------------------------------------------------!
  character(len=*), intent(in) :: domain
  real, intent(inout) :: lon,lat
  integer, intent(inout),optional :: iloc,jloc
  integer, intent(out)  ,optional :: iglob,jglob
  logical :: in
  integer :: i,j
  call coord_check("coord_in_"//trim(domain),lon,lat,fix=.true.)
  call lb2ij(lon,lat,i,j)
  if(present(iglob))iglob=i
  if(present(jglob))jglob=j
  in=(i>=1).and.(i<=IIFULLDOM).and.(j>=1).and.(j<=JJFULLDOM)
  i=max(1,min(i,IIFULLDOM));i=i_local(i)
  j=max(1,min(j,JJFULLDOM));j=j_local(j)
  select case(domain)
  case("g","G","global","full")
     if(present(iloc))iloc=i
     if(present(jloc))jloc=j
  case("l","L","local","processor")
     if(in) in=(i>=1).and.(i<=limax).and.(j>=1).and.(j<=ljmax)
     if(present(iloc))iloc=i
     if(present(jloc))jloc=j
  case("c","C","cell","gridbox")
     call CheckStop(.not.(present(iloc).and.present(jloc)),&
          "Wrong options for coord_in_"//trim(domain))
     if(in) in=(i==iloc).and.(j==jloc)
  case default
     call CheckStop("Unsupporter coord_in_"//trim(domain))
  end select
end function coord_in_domain
function coord_in_processor(lon,lat,iloc,jloc,iglob,jglob) result(in)
  !-------------------------------------------------------------------!
  ! Is coord (lon/lat) is inside local domain?
  !-------------------------------------------------------------------!
  real, intent(inout) :: lon,lat
  integer, intent(out),optional:: iloc,jloc,iglob,jglob
  logical :: in
  in=coord_in_domain("processor",lon,lat,iloc,jloc,iglob,jglob)
end function coord_in_processor
function coord_in_gridbox(lon,lat,iloc,jloc,iglob,jglob) result(in)
  !-------------------------------------------------------------------!
  ! Is coord (lon/lat) is inside gridbox(iloc,jloc)?
  !-------------------------------------------------------------------!
  real, intent(inout) :: lon,lat
  integer, intent(inout) :: iloc,jloc
  integer, intent(out),optional:: iglob,jglob
  logical :: in
  in=coord_in_domain("gridbox",lon,lat,iloc,jloc,iglob,jglob)
end function coord_in_gridbox

subroutine Alloc_GridFields(LIMAX,LJMAX,KMAX_MID,KMAX_BND)

  integer, intent(in)::LIMAX,LJMAX,KMAX_MID,KMAX_BND

  allocate(i_fdom(0:LIMAX+1))
  allocate(j_fdom(0:LJMAX+1))
  allocate(glon(LIMAX,LJMAX))
  allocate(glat(LIMAX,LJMAX))
  allocate(gl_stagg(0:LIMAX,0:LJMAX))
  allocate(gb_stagg(0:LIMAX,0:LJMAX))
  allocate(xm_i(0:LIMAX+1,0:LJMAX+1))
  allocate(xm_j(0:LIMAX+1,0:LJMAX+1))
  allocate(xm2(0:LIMAX+1,0:LJMAX+1))
  allocate(xmd(0:LIMAX+1,0:LJMAX+1))
  allocate(xm2ji(0:LJMAX+1,0:LIMAX+1))
  allocate(xmdji(0:LJMAX+1,0:LIMAX+1))
  allocate(GridArea_m2(LIMAX,LJMAX))
  allocate(z_bnd(LIMAX,LJMAX,KMAX_BND))
  allocate(z_mid(LIMAX,LJMAX,KMAX_MID))

end subroutine Alloc_GridFields

subroutine make_vertical_levels_interpolation_coeff
  ! make interpolation coefficients to convert the levels defined in meteo
  ! into the levels defined in Vertical_levelsFile

  integer ::k,k_met
  real ::p_met,p_mod,p1,p2
  if(.not. allocated(k1_met))allocate(k1_met(KMAX_MID),k2_met(KMAX_MID),x_k1_met(KMAX_MID))
  if(.not. allocated(A_bnd_met))then
    allocate(A_bnd_met(KMAX_MID+1),B_bnd_met(KMAX_MID+1))
    A_bnd_met=A_bnd
    B_bnd_met=B_bnd
  end if
  
  if(me_mpi==0)then
    ! only me=0 has the values for A_bnd_met and B_bnd_met
    do k=1,KMAX_MID
      P_mod=A_mid(k)+Pref*B_mid(k)
      !find the lowest met level higher than the model level
      !do k_met=1,KMAX_MET
      k_met=KMAX_MET-1
      p_met=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
      do while(p_met>P_mod.and.k_met>1)
        !if(MasterProc) write(*,*)P_mod,p_met
        k_met=k_met-1
        p_met=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
      end do
      k1_met(k)=k_met
      k2_met(k)=k_met+1
      k_met=k1_met(k)
      p1=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
      k_met=k2_met(k)
      p2=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
      x_k1_met(k)=(p_mod-p2)/(p1-p2)
      write(*,77)k, ' interpolated from levels ', k1_met(k),' and ',k2_met(k),P_mod,p1,p2,x_k1_met(k)
77    format(I4,A,I3,A,I3,13f11.3)
      if(x_k1_met(k)<-0.00001 .or. (1.0-x_k1_met(k))<-0.00001)then
        write(*,*)'WARNING: Extrapolation of data. This is NOT recommended for some metfields'
      end if
      
    end do
  end if
  
  CALL MPI_BCAST(k1_met,4*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(k2_met,4*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(x_k1_met,8*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)

end subroutine make_vertical_levels_interpolation_coeff

subroutine Meteo_Get_KMAXMET(filename, KMAX, ncfileID_in)

  character(len=*), intent(in)  :: filename 
  integer, intent(out)  :: KMAX
  integer, intent(in), optional  :: ncfileID_in
  integer :: kdimid, ncFileID, status

  if(present(ncfileID_in))then
     ncfileID = ncfileID_in
  else
     status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
     if(status/=nf90_noerr) then
        print *,'not found',trim(filename)
        call CheckStop("GridValues get KMAX: File not found")
     end if
  endif

  status=nf90_inq_dimid(ncid=ncFileID, name="k", dimID=kdimID)
  if(status/=nf90_noerr)then
     status=nf90_inq_dimid(ncid=ncFileID, name="lev", dimID=kdimID)!hybrid coordinates                                               
     if(status/=nf90_noerr) then
        status=nf90_inq_dimid(ncid=ncFileID, name="hybrid", dimID=kdimID)!hybrid coordinates                                         
        if(status/=nf90_noerr) then ! WRF format                                                                                     
           call check(nf90_inq_dimid(ncid=ncFileID, name="bottom_top", dimID=kdimID))
           end if
     end if
  end if
  call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX))

  if(.not.present(ncfileID_in))call check(nf90_close(ncFileID))  

end subroutine Meteo_Get_KMAXMET

subroutine remake_vertical_levels_interpolation_coeff(filename)
  ! make again interpolation coefficients to convert the levels defined in meteo
  ! into the levels defined in Vertical_levelsFile
  ! this is used if the number of vertical levels has changed from one meteo 
  ! file to the next

  character(len=*), intent(in)  :: filename 
  integer ::k,kk,k_met
  real ::p_met,p_mod,p1,p2
  real :: P_TOP_MET,P0
  integer :: status, ncFileID, kdimID, varID, ios

  status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
  if(status/=nf90_noerr) then
    print *,'not found',trim(filename)
    call CheckStop("GridValues: File not found")
  end if
  if(MasterProc)print *,'Redefining vertical interpolation parameters for ',trim(filename)

  if(.not. allocated(k1_met))allocate(k1_met(KMAX_MID),k2_met(KMAX_MID),x_k1_met(KMAX_MID))
  if(allocated(A_bnd_met))then
    deallocate(A_bnd_met,B_bnd_met)
  end if

  call Meteo_Get_KMAXMET(filename, KMAX_MET, ncfileID)

  status=nf90_inq_varid(ncid=ncFileID, name="k", varID=varID)
  if(status/=nf90_noerr)then
    if(MasterProc)write(*,*)'reading met hybrid levels from ',trim(filename)
    status=nf90_inq_varid(ncid=ncFileID, name="P0", varID=varID)
    if(status/=nf90_noerr)&
      status=nf90_inq_varid(ncid=ncFileID, name="p0", varID=varID)
    if(status/=nf90_noerr)then
      status=nf90_inq_varid(ncid=ncFileID, name="P00", varID=varID) !WRF case
      if(status/=nf90_noerr)then
          call StopAll('Do not know how to define vertical levels')
      else
        ! WRF format
        ! asuming sigma levels ZNW=(P-P_TOP_MET)/(PS-P_TOP_MET)
        ! P = A+B*PS = P_TOP_MET*(1-ZNW) + ZNW*PS
        ! B = ZNW
        ! A = P_TOP_MET*(1-ZNW)
        call check(nf90_get_var(ncFileID, varID, P0 ))
        if(.not.allocated(A_bnd_met))allocate(A_bnd_met(KMAX_MET+1),B_bnd_met(KMAX_MET+1))
        call check(nf90_inq_varid(ncid=ncFileID, name="P_TOP", varID=varID))
        call check(nf90_get_var(ncFileID, varID, P_TOP_MET ))
        call check(nf90_inq_varid(ncid=ncFileID, name="ZNW", varID=varID))
        call check(nf90_get_var(ncFileID, varID, B_bnd_met ))
        if(MET_REVERSE_K)then
          A_bnd_met=B_bnd_met!use A_bnd_met as temporary buffer
          do k=1,KMAX_MET+1
            B_bnd_met(k)=A_bnd_met(KMAX_MET+2-k)
          end do
        end if
        A_bnd_met=P_TOP_MET*(1.-B_bnd_met)
      end if
      if(MET_REVERSE_K)then
        if(MasterProc)write(*,*)"Reversed vertical levels from met, P at level boundaries:"
      else
        if(MasterProc)write(*,*)"Vertical levels from met, P at level boundaries:"
      end if
      do k=1,KMAX_MET+1
        if(MasterProc)write(*,44)k, A_bnd_met(k)+P0*B_bnd_met(k)
      end do
    else
      call check(nf90_get_var(ncFileID, varID, P0 ))
      if(MasterProc)write(*,*)'P0 = ',P0
      if(.not.allocated(A_bnd_met))allocate(A_bnd_met(KMAX_MET+1),B_bnd_met(KMAX_MET+1))
      status=nf90_inq_varid(ncid=ncFileID, name="hyai", varID=varID)
      if(status/=nf90_noerr)then
        call check(nf90_inq_varid(ncid=ncFileID, name="ap", varID=varID))
        call check(nf90_get_var(ncFileID, varID, A_bnd_met, count=(/KMAX_MET/)))!read mid values!
        call check(nf90_inq_varid(ncid=ncFileID, name="b", varID=varID))
        call check(nf90_get_var(ncFileID, varID, B_bnd_met, count=(/KMAX_MET/) )) !read mid values!
        A_bnd_met(KMAX_MET+1)=0.0
        B_bnd_met(KMAX_MET+1)=1.0
        do k=KMAX_MET,1,-1
          A_bnd_met(k)=A_bnd_met(k+1)-2.0*(A_bnd_met(k+1)-A_bnd_met(k))!from mid to bnd values!
          B_bnd_met(k)=B_bnd_met(k+1)-2.0*(B_bnd_met(k+1)-B_bnd_met(k))!from mid to bnd values!
        end do
        
        if(MasterProc)write(*,*)'Metdata pressure at level boundaries:'
        do k=1,KMAX_MET+1
          if(MasterProc)write(*,44)k, A_bnd_met(k)+P0*B_bnd_met(k)
        end do
        
      else
        call check(nf90_get_var(ncFileID, varID, A_bnd_met ))
        A_bnd_met=P0*A_bnd_met!different definition in model and grid_Def
        call check(nf90_inq_varid(ncid=ncFileID, name="hybi", varID=varID))
        call check(nf90_get_var(ncFileID, varID, B_bnd_met ))
      end if
    end if
    if(External_Levels_Def)then
       !model levels defined from external text file
       if(MasterProc)&
            write(*,*)'reading external hybrid levels from ',trim(Vertical_levelsFile),&
            A_bnd_met(kMAX_met+1),B_bnd_met(kMAX_met+1)
       P0=Pref
       open(IO_TMP,file=trim(Vertical_levelsFile),action="read",iostat=ios)
       read(IO_TMP,*)k
       if(k/=KMAX_MID .and. MasterProc)write(*,*)k,kmax_mid,&
            'WARNING: unexpected number of levels for '//trim(Vertical_levelsFile)
       do k=1,KMAX_MID+1
          read(IO_TMP,*)kk,A_bnd(k),B_bnd(k)
          if(kk/=k.and.MasterProc)write(*,*)'WARNING: unexpected format for vertical levels ',k,kk
       end do
       close(IO_TMP)
    else
      !vertical model levels are the same as in meteo
      A_bnd=A_bnd_met
      B_bnd=B_bnd_met
    end if
    
    do k=1,KMAX_MID
      A_mid(k)=0.5*(A_bnd(k)+A_bnd(k+1))
      B_mid(k)=0.5*(B_bnd(k)+B_bnd(k+1))
    end do
    sigma_mid =B_mid!for Hybrid coordinates sigma_mid=B if A*P0=PT-sigma_mid*PT
    
    if(MasterProc)write(*,*)"Model pressure at level boundaries:"
    do k=1,KMAX_MID+1
44    FORMAT(i4,10F12.2)
      if(MasterProc)write(*,44)k, A_bnd(k)+P0*B_bnd(k)
    end do
    !test if the top is within the height defined in the meteo files
    if(MasterProc.and.External_Levels_Def.and.(A_bnd(1)+P0*B_bnd(1)+0.01<A_bnd_met(1)+P0*B_bnd_met(1)))then
      write(*,*)'Pressure at top of defined levels is ',A_bnd(1)+P0*B_bnd(1)
      write(*,*)'Pressure at top defined in meteo files is ',A_bnd_met(1)+P0*B_bnd_met(1)
      write(*,*)'Pressure at op must be higher (lower altitude) than top defined in meteo '
      call StopAll('Top level too high! Change values in '//trim(Vertical_levelsFile))
    end if
    
    !test if the levels can cope with highest mountains (400 hPa)
    do k=1,KMAX_MID
      if(MasterProc.and.A_bnd(k+1)+40000*B_bnd(k+1)-(A_bnd(k)+40000*B_bnd(k))<0.0)then
        write(*,*)'WARNING: hybrid vertical level definition may cause negative level thickness when pressure below 400 hPa '
        write(*,*)'Pressure at level ',k,' is ',A_bnd(k)+40000*B_bnd(k)
        write(*,*)'Pressure at level ',k+1,' is ',A_bnd(k+1)+40000*B_bnd(k+1),' (should be higher)'
        if(External_Levels_Def)call StopAll('GridValues_mod: possible negative level thickness ')
      end if
    end do
    
    !test if the lowest levels is thick enough (twice height of highest vegetation?) about 550 Pa = about 46m
    !Deposition scheme is not designes for very thin lowest levels
    if(MasterProc.and.A_bnd(KMAX_MID+1)+P0*B_bnd(KMAX_MID+1)-(A_bnd(KMAX_MID)+P0*B_bnd(KMAX_MID))<550.0)then
      write(*,*)'WARNING: lowest level very shallow; ',A_bnd(KMAX_MID+1)+P0*B_bnd(KMAX_MID+1) -&
      (A_bnd(KMAX_MID)+P0*B_bnd(KMAX_MID)),'Pa'
!      call StopAll('Lowest level too thin! Change vertical levels definition in '//trim(Vertical_levelsFile))
    end if
    
  end if
  call check(nf90_close(ncFileID))
  
  if(me_mpi==0)then
    ! only me=0 has the values for A_bnd_met and B_bnd_met
    do k=1,KMAX_MID
      P_mod=A_mid(k)+Pref*B_mid(k)
      !find the lowest met level higher than the model level
      !do k_met=1,KMAX_MET
      k_met=KMAX_MET-1
      p_met=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
      do while(p_met>P_mod.and.k_met>1)
        !if(MasterProc) write(*,*)P_mod,p_met
        k_met=k_met-1
        p_met=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
      end do
      k1_met(k)=k_met
      k2_met(k)=k_met+1
      k_met=k1_met(k)
      p1=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
      k_met=k2_met(k)
      p2=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
      x_k1_met(k)=(p_mod-p2)/(p1-p2)
      write(*,77)k, ' interpolated from levels ', k1_met(k),' and ',k2_met(k),P_mod,p1,p2,x_k1_met(k)
77    format(I4,A,I3,A,I3,13f11.3)
      if(x_k1_met(k)<-0.00001 .or. (1.0-x_k1_met(k))<-0.00001)then
        write(*,*)'WARNING: Extrapolation of data. This is NOT recommended for some metfields'
      end if
      
    end do
  end if
  
  CALL MPI_BCAST(k1_met,4*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(k2_met,4*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_BCAST(x_k1_met,8*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)

end subroutine remake_vertical_levels_interpolation_coeff

subroutine lambert2lb(x,y,gl,gb,lon0,y0,k,F)
  real, intent(in) ::x,y,lon0,y0,k,F
  real, intent(out)::gl,gb
  real ::r,t
  r = sqrt(x*x+(y0-y)*(y0-y))
  if(k<0.0)r = -r
  t = atan(x/(y0-y))
  gb = 2*180./PI*atan((F/r)**(1.0/k))-90.0
  gl = lon0 + 180./PI*t/k
end subroutine lambert2lb
subroutine lb2lambert(x,y,gl,gb,lon0,y0,k,F)
  real, intent(in) ::gl,gb,lon0,y0,k,F
  real, intent(out)::x,y
  real ::r
  r = F*tan(0.25*PI-0.5*deg2rad*gb)**k ! depends only on latitude ->reuse (about 100 cycles/operation)
  x = r*sin(deg2rad*k*(gl-lon0))
  y = y0 - r*cos(deg2rad*k*(gl-lon0))
end subroutine lb2lambert

subroutine lb_rot2lb(xsph,ysph,xrot,yrot,grid_north_pole_longitude,grid_north_pole_latitude)
  !  compute spherical coordinates as function of
  !  spherical rotated coordinates
  !
  !  conversion between spherical (xsph,ysph) and spherical rotated
  !  (xrot,yrot) coordinates. (xcen,ycen) is the position of the
  !  rotated equator/greenwich in terms of (longitude,latitude).
  !  all input and output values are given in degrees.
  !
  ! grid_north_pole_longitude: geographical (non-rotated) coordinates
  ! of the "north pole" from the rotated grid (No polar bears there).
  ! (typically out of the grid, since it is singular).
  !
  ! xcen: geographical (non-rotated) coordinates of the (lon=0 lat=0)
  ! point where lonlat are in the rotated grid
  ! (typically in the middle of the grid, since it is "flat")

  implicit none
  real :: xsph, ysph, xrot, yrot,xcen,ycen,zsycen,zcycen
  real :: zsxrot,zcxrot,zsyrot,zcyrot,zsysph,zcysph,zcxmxc,zsxmxc,zxmxc
  real :: grid_north_pole_longitude,grid_north_pole_latitude
  real :: rad2deg,deg2rad
  !
  deg2rad=3.14159265358979323/180.
  rad2deg=1.0/deg2rad
  
  xcen=(180.+grid_north_pole_longitude)*deg2rad
  ycen=(90.-grid_north_pole_latitude)*deg2rad
  
  zsycen = sin(ycen)
  zcycen = cos(ycen)
  
  zsxrot = sin(xrot*deg2rad)
  zcxrot = cos(xrot*deg2rad)
  zsyrot = sin(yrot*deg2rad)
  zcyrot = cos(yrot*deg2rad)
  zsysph = zcycen*zsyrot + zsycen*zcyrot*zcxrot
  zsysph = amax1(zsysph,-1.0)
  zsysph = amin1(zsysph,+1.0)
  ysph = asin(zsysph)
  zcysph = cos(ysph)
  zcxmxc = (zcycen*zcyrot*zcxrot -&
  zsycen*zsyrot)/zcysph
  zcxmxc = amax1(zcxmxc,-1.0)
  zcxmxc = amin1(zcxmxc,+1.0)
  zsxmxc = zcyrot*zsxrot/zcysph
  zxmxc  = acos(zcxmxc)
  if (zsxmxc.lt.0.0) zxmxc = -zxmxc
  xsph = (zxmxc + xcen)*rad2deg
  ysph = ysph*rad2deg
  
end subroutine lb_rot2lb

subroutine lb2UTM(Long, Lat, UTMEasting, UTMNorthing, UTMZone)
  ! converts lat/long to UTM coords.  Equations from USGS Bulletin 1532
  ! East Longitudes are positive, West longitudes are negative.
  ! North latitudes are positive, South latitudes are negative
  ! Lat and Long are in decimal degrees
  
  ! Angle with North: (called "convergence")
  ! angle = arctan(tan(lon)*sin(lat)) lon is longitude relative
  !         to middle of utm zone: lon=gl-Lambda0
  ! works for northern hemisphere at least
  ! u_utm = u_ll*cos(angle)+v_ll*sin(angle)
  ! v_utm =-u_ll*sin(angle)+v_ll*cos(angle)

  implicit none

  real :: UTMNorthing, UTMEasting,  Lat,  Long
  integer :: UTMZone
  real :: a = 6378137.0 !WGS-84
  real :: eccSquared = 0.00669438 !WGS-84
  real :: k0 = 0.9996
  real :: LongOrigin
  real :: eccPrimeSquared
  real :: N, T, C, AA, M
  real :: rad2deg,deg2rad
  
  real :: LongTemp
  real :: LatRad
  real :: LongRad
  real :: LongOriginRad;
  
  rad2deg=180.0/Pi
  deg2rad=Pi/180.
  LatRad = Lat*deg2rad
  !//Make sure the longitude is between -180.00 .. 179.9
  LongTemp = (Long+180)-int((Long+180)/360)*360-180!; // -180.00 .. 179.9;
  LongRad = LongTemp*deg2rad
  UTMZone = int((LongTemp + 180)/6) + 1;
  
  !Southern Norway, zone 32 is extended by 3 degrees to the West
  if( Lat >= 56.0 .and. Lat < 64.0  .and. LongTemp >= 3.0 .and. LongTemp < 12.0 )UTMZone  = 32
  
  !// Special zones for Svalbard
  if( Lat >= 72.0 .and. Lat < 84.0 ) then
    if(      LongTemp >= 0.0  .and. LongTemp <  9.0 )then
      UTMZone = 31
    else if( LongTemp >= 9.0  .and. LongTemp < 21.0 )then
      UTMZone = 33
    else if( LongTemp >= 21.0 .and. LongTemp < 33.0 )then
      UTMZone = 35
    else if( LongTemp >= 33.0 .and. LongTemp < 42.0 )then
      UTMZone = 37
    endif
  endif
  LongOrigin = (UTMZone - 1)*6 - 180 + 3!  //+3 puts origin in middle of zone
  LongOriginRad = LongOrigin * deg2rad
  
  !//compute the UTM Zone from the latitude and longitude
  
  eccPrimeSquared = (eccSquared)/(1-eccSquared)
  
  N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad))
  T = tan(LatRad)*tan(LatRad)
  C = eccPrimeSquared*cos(LatRad)*cos(LatRad)
  AA = cos(LatRad)*(LongRad-LongOriginRad)
  
  M = a*((1 - eccSquared/4 - 3*eccSquared*eccSquared/64 &
    - 5*eccSquared*eccSquared*eccSquared/256)*LatRad &
    - (3*eccSquared/8 + 3*eccSquared*eccSquared/32 &
      + 45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)&
    + (15*eccSquared*eccSquared/256 &
      + 45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad) &
    - (35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad))
  
  UTMEasting = k0*N*(AA+(1-T+C)*AA*AA*AA/6&
  +(5-18*T+T*T+72*C-58*eccPrimeSquared)*AA*AA*AA*AA*AA/120)+ 500000.0
  
  UTMNorthing = (k0*(M+N*tan(LatRad)*(AA*AA/2+(5-T+9*C+4*C*C)*AA*AA*AA*AA/24+&
  (61-58*T+T*T+600*C-330*eccPrimeSquared)*AA*AA*AA*AA*AA*AA/720)))
  if(Lat < 0)UTMNorthing =UTMNorthing + 10000000.0!; //10000000 meter offset for southern hemisphere
  
end subroutine lb2UTM

subroutine UTM2lb(UTMEasting, UTMNorthing, UTMZone,  Long, Lat )
  ! converts UTM coords to lat/long.  Equations from USGS Bulletin 1532
  ! East Longitudes are positive, West longitudes are negative.
  ! North latitudes are positive, South latitudes are negative
  ! Lat and Long are in decimal degrees.

  implicit none
  
  real :: UTMNorthing, UTMEasting,  Lat,  Long
  integer :: UTMZone
  real :: k0 = 0.9996
  real :: a = 6378137.0 !WGS-84
  real :: eccSquared = 0.00669438 !WGS-84
  real :: eccPrimeSquared
  real :: e1,rad2deg
  real :: N1, T1, C1, R1, D, M
  real :: LongOrigin
  real :: mu, phi1, phi1Rad
  real :: x, y, ww
  !char* ZoneLetter;
  integer :: NorthernHemisphere !1 for northern hemispher, 0 for southern
  
  rad2deg=180.0/Pi
  
  x = UTMEasting - 500000.0 !remove 500,000 meter offset for longitude
  y = UTMNorthing
  
  NorthernHemisphere = 1 !point is in northern hemisphere
  
  LongOrigin = (UTMZone - 1)*6 - 180 + 3 !+3 puts origin in middle of zone
  
  e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared))
  
  eccPrimeSquared = (eccSquared)/(1-eccSquared)
  
  M = y / k0
  mu = M/(a*(1-eccSquared/4-3*eccSquared*eccSquared/64-5*eccSquared*eccSquared*eccSquared/256))
  
  phi1Rad = mu + (3*e1/2-27*e1*e1*e1/32)*sin(2*mu)+ (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)+(151*e1*e1*e1/96)*sin(6*mu)
  phi1 = phi1Rad*rad2deg
  
  N1 = a/sqrt(1-eccSquared*sin(phi1Rad)*sin(phi1Rad))
  T1 = tan(phi1Rad)*tan(phi1Rad)
  C1 = eccPrimeSquared*cos(phi1Rad)*cos(phi1Rad)
  ww = 1-eccSquared*sin(phi1Rad)*sin(phi1Rad)
  R1 = a*(1-eccSquared)/(ww*sqrt(ww))
  D = x/(N1*k0)
  
  Lat = phi1Rad - (N1*tan(phi1Rad)/R1)*(D*D/2-(5+3*T1+10*C1-4*C1*C1-9*eccPrimeSquared)*D*D*D*D/24&
  +(61+90*T1+298*C1+45*T1*T1-252*eccPrimeSquared-3*C1*C1)*D*D*D*D*D*D/720)
  Lat = Lat * rad2deg
  
  Long = (D-(1+2*T1+C1)*D*D*D/6+(5-2*C1+28*T1-3*C1*C1+8*eccPrimeSquared+24*T1*T1)&
  *D*D*D*D*D/120)/cos(phi1Rad)
  Long = LongOrigin + Long * rad2deg

end subroutine UTM2lb

subroutine nf90_get_var_extended(ncFileID,varID,Var,i0,i1,j0,j1,ishift_in,jshift_in)
  ! fetch a 2D array over an extended subdomain.
  ! i.e. extended arrays are overlapping, and parts outside the fulldomain
  ! are extrapolated linearly.
  implicit none
  integer, intent(in) :: ncFileID,varID,i0,i1,j0,j1
  real, intent(inout) :: Var(i0:i1,j0:j1) !the extended local array
  integer, optional, intent(in)::ishift_in,jshift_in

  integer ::iloc_start,iloc_end,jloc_start,jloc_end
  integer ::i,j,ishift,jshift

  iloc_start=i0
  iloc_end=i1
  jloc_start=j0
  jloc_end=j1
  ishift=0
  jshift=0
  if(present(ishift_in))ishift=ishift_in
  if(present(jshift_in))jshift=jshift_in

  if(iloc_start+IRUNBEG+gi0-2<1)iloc_start=1!first cell (in i direction)
  if(iloc_end+IRUNBEG+gi0-2>IIFULLDOM)iloc_end=IIFULLDOM+2-gi0-IRUNBEG!last cell
  if(jloc_start+JRUNBEG+gj0-2<1)jloc_start=1!first cell (in j direction)
  if(jloc_end+JRUNBEG+gj0-2>JJFULLDOM)jloc_end=JJFULLDOM+2-gj0-JRUNBEG!last cell
  
  call check(nf90_get_var(ncFileID, varID, Var(iloc_start:iloc_end,jloc_start:jloc_end) &
  ,start=(/iloc_start+IRUNBEG+gi0-2+ishift,jloc_start+JRUNBEG+gj0-2+jshift/),&
  count=(/iloc_end-iloc_start+1, jloc_end-jloc_start+1/)))
  
  !extrapolate if needed at fulldomain boundaries
  do i=iloc_start-1,i0,-1
    do j=jloc_start,jloc_end
      Var(i,j)=2.0*Var(i+1,j)-Var(i+2,j)
    end do
  end do
  do i=iloc_end+1,i1
    do j=jloc_start,jloc_end
      Var(i,j)=2.0*Var(i-1,j)-Var(i-2,j)
    end do
  end do
  do i=i0,i1!now they are all defined
    do j=jloc_start-1,j0,-1
      Var(i,j)=2.0*Var(i,j+1)-Var(i,j+2)
    end do
  end do
  do i=i0,i1!now they are all defined
    do j=jloc_end+1,j1
      Var(i,j)=2.0*Var(i,j-1)-Var(i,j-2)
    end do
  end do

end subroutine nf90_get_var_extended

subroutine RestrictDomain(DOMAIN)
  integer, dimension(4), intent(inout)::  DOMAIN
  ! allow for wildcards,
  ! eg [10,20,-1,-1] => [10,20,RUNDOMAIN(3),RUNDOMAIN(4)]
  where(DOMAIN<0) DOMAIN=RUNDOMAIN

  ! Ensure DOMAIN is contained by RUNDOMAIN
  DOMAIN(1)=max(RUNDOMAIN(1),DOMAIN(1))
  DOMAIN(2)=min(RUNDOMAIN(2),DOMAIN(2))
  DOMAIN(3)=max(RUNDOMAIN(3),DOMAIN(3))
  DOMAIN(4)=min(RUNDOMAIN(4),DOMAIN(4))

  ! consistency check
  if(any([DOMAIN==0,DOMAIN(1)>DOMAIN(2),DOMAIN(3)>DOMAIN(4)]))then
    write(*,"(A,'=[',I0,3(',',I0),']')")'Inconsistent DOMAIN',DOMAIN
    call CheckStop('Inconsistent DOMAIN')
  end if
end subroutine RestrictDomain

subroutine extendarea_N(f,h,thick,Size1,Size2,debug_flag)
  ! returns extended array array, reading neighbour procs as needed
  ! size of h MUST be as declared below

  integer, intent(in) :: thick,Size1,Size2
  real, intent(in) :: f(Size1,LIMAX,LJMAX,Size2)
  real, intent(inout) :: h(Size1,1-thick:LIMAX+thick,1-thick:LJMAX+thick,Size2)
  logical, intent(in), optional :: debug_flag
  logical :: mydebug = .false.

  real, dimension(Size1,LIMAX,thick,Size2)            :: f_south,f_north
  real, dimension(Size1,thick,1-thick:LJMAX+thick,Size2)        :: f_west,f_east

  integer :: iif,jjf,i,j,iifl,jjfl,i1,i2
  if ( present(debug_flag)  ) mydebug = debug_flag
  
  ! readneighbours twice
  iifl=limax+2*thick
  jjfl=ljmax+2*thick
  if(mydebug .and. MasterProc ) write(*,*) "DEBUG extendarea", iif,jjf,thick
  
  call readneighbors_N(f,f_south,f_north,f_west,f_east,thick,Size1,Size2)
  
  do i2=1,Size2
    do j=1,ljmax
      do i=1,limax
        do i1=1,Size1
          h(i1,i,j,i2) = f(i1,i,j,i2)
        enddo
      end do
    end do
  end do
  do i2=1,Size2
    do j=1,thick
      do i=1,limax
        do i1=1,Size1
          h(i1,i,j-thick,i2) = f_south(i1,i,j,i2)
        end do
      end do
    end do
  end do
  
  do i2=1,Size2
    do j=1,thick
      do i=1,limax
        do i1=1,Size1
          h(i1,i,ljmax+j,i2) = f_north(i1,i,j,i2)
        end do
      end do
    end do
  end do
  
  do i2=1,Size2
    do j=1-thick,ljmax+thick
      do i=1,thick
        do i1=1,Size1
          h(i1,i-thick,j,i2) = f_west(i1,i,j,i2)
        end do
      end do
    end do
  end do
  
  do i2=1,Size2
    do j=1-thick,ljmax+thick
      do i=1,thick
        do i1=1,Size1
          h(i1,limax+i,j,i2) = f_east(i1,i,j,i2)
        end do
      end do
    end do
  end do

end subroutine extendarea_N

subroutine readneighbors_N(data,data_south,data_north,data_west,data_east,thick,Size1,Size2)
  ! Read data at the other side of the boundaries
  !
  ! thick is the number of gridcells in each direction to be transferred
  ! Note that we also fetch data from processors in the "diagonal"
  ! directions
  !
  ! Written by Peter January 2017; remote neighbors June 2017

  ! Note,
  ! The data_west(jj,:)=data(1,j) is not a bug: when there is no west
  ! neighbour,
  ! the data is simply copied from the nearest points: data_west(jj,:) should
  ! be =data(-thick+1:0,j), but since this data does not exist, we
  ! put it =data(1,j).

  implicit none
  integer, intent(in) :: thick,Size1,Size2
  real,intent(in), dimension(Size1,LIMAX,LJMAX,Size2) ::data
  real,intent(out), dimension(Size1,LIMAX,thick,Size2) ::data_south,data_north
  real,intent(out), dimension(Size1,thick,1-thick:LJMAX+thick,Size2) ::data_west,data_east
  real, dimension(Size1,LIMAX,min(thick,MAXLJMAX),Size2) ::data_south_snd,data_north_snd
  real, dimension(Size1,min(thick,MAXLIMAX),1-thick:LJMAX+thick,Size2) ::data_west_snd,data_east_snd
  real, dimension(Size1,LIMAX,min(MAXLJMAX,thick),Size2) ::data_sn_rcv
  real, dimension(Size1,min(MAXLIMAX,thick),1-thick:LJMAX+thick,Size2) ::data_we_rcv
  
  integer :: msgnr
  integer :: i,it,j,jj,jt,i1,i2,n
  integer :: mythick,myithick,myjthick,totthick,ineighbor,limaxloc,ljmaxloc
  
  !check that limax and ljmax are large enough. Can only read neighboring subdomain
  !    call CheckStop(limax < thick, "ERROR readneighbors_N in Met_mod")
  !    call CheckStop(ljmax < thick, "ERROR readneighbors_N in Met_mod")
  limaxloc = min(MAXLIMAX,thick)!array size common for all
  ljmaxloc = min(MAXLJMAX,thick)!array size common for all
  
  msgnr=1
  myjthick = min(thick,LJMAX)
  myithick = min(thick,LIMAX)
  !    data_south_snd(:,:,:,:)=data(:,:,1:thick,:)
  !    data_north_snd(:,:,:,:)=data(:,:,max(1,ljmax-thick+1):ljmax,:)
  do i2=1,Size2
    do jt=1,myjthick
      do i=1,limax
        do i1=1,Size1
          data_south_snd(i1,i,jt,i2)=data(i1,i,jt,i2)
        end do
      end do
    end do
  end do
  do i2=1,Size2
    do jt=1,myjthick
      do i=1,limax
        do i1=1,Size1
          data_north_snd(i1,i,jt,i2)=data(i1,i,LJMAX-myjthick+jt,i2)
        end do
      end do
    end do
  end do
  
  if(neighbor(SOUTH) >= 0 )then
    if(thick<MINLJMAX)then
      CALL MPI_ISEND( data_south_snd , 8*LIMAX*thick*Size1*Size2, MPI_BYTE,&
      neighbor(SOUTH), msgnr, MPI_COMM_CALC, request_s,IERROR)
    else
      ineighbor=neighbor(SOUTH)
      totthick=0
      do n=1,NPROCY
        mythick=min(tljmax(ineighbor),thick-totthick)
        CALL MPI_ISEND( data_south_snd , 8*LIMAX*ljmaxloc*Size1*Size2, MPI_BYTE,&
        ineighbor, msgnr, MPI_COMM_CALC, irequest_s(n),IERROR)
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor-NPROCX
        if(ineighbor<0) exit
      enddo
      
    endif
  end if
  if(neighbor(NORTH) >= 0 )then
    if(thick<MINLJMAX)then
      CALL MPI_ISEND( data_north_snd , 8*LIMAX*thick*Size1*Size2, MPI_BYTE,&
      neighbor(NORTH), msgnr+9, MPI_COMM_CALC, request_n,IERROR)
    else
      ineighbor=neighbor(NORTH)
      totthick=0
      do n=1,NPROCY
        mythick=min(tljmax(ineighbor),thick-totthick)
        CALL MPI_ISEND( data_north_snd , 8*LIMAX*ljmaxloc*Size1*Size2, MPI_BYTE,&!NB: send everything anyway!
        ineighbor, msgnr+9, MPI_COMM_CALC, irequest_n(n),IERROR)
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor+NPROCX
        if(ineighbor>=NPROC) exit
      enddo
      
    endif
  end if
  
  if(neighbor(SOUTH) >= 0 )then
    if(thick<MINLJMAX)then
      CALL MPI_RECV( data_south, 8*LIMAX*thick*Size1*Size2, MPI_BYTE,&
      neighbor(SOUTH), msgnr+9, MPI_COMM_CALC, MPISTATUS, IERROR)
    else
      ineighbor=neighbor(SOUTH)
      totthick=0
      do n=1,NPROCY
        mythick=min(tljmax(ineighbor),thick-totthick)
        CALL MPI_RECV( data_sn_rcv, 8*LIMAX*ljmaxloc*Size1*Size2, MPI_BYTE,&
        ineighbor, msgnr+9, MPI_COMM_CALC, MPISTATUS, IERROR)
        do i2=1,Size2
          do jt=1,mythick
            do i=1,limax
              do i1=1,Size1
                data_south(i1,i,thick-totthick-mythick+jt,i2)=data_sn_rcv(i1,i,tljmax(ineighbor)-mythick+jt,i2)
              end do
            end do
          end do
        end do
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor-NPROCX
        if(ineighbor<0)then
          !must fill remainder data, even if out of rundomain
          do i2=1,Size2
            do jt=1,thick-totthick
              do i=1,limax
                do i1=1,Size1
                  data_south(i1,i,jt,i2)=data_south(i1,i,thick-totthick+1,i2)
                end do
              end do
            end do
          end do
          exit
        endif
      enddo
    endif
  else
    do i2=1,Size2
      do jt=1,thick
        do i=1,limax
          do i1=1,Size1
            data_south(i1,i,jt,i2)=data(i1,i,1,i2)
          end do
        end do
      end do
    end do
  end if
  
  if(neighbor(NORTH) >= 0 )then
    if(thick<MINLJMAX)then
      CALL MPI_RECV( data_north, 8*LIMAX*thick*Size1*Size2, MPI_BYTE,&
      neighbor(NORTH), msgnr, MPI_COMM_CALC, MPISTATUS, IERROR)
    else
      ineighbor=neighbor(NORTH)
      totthick=0
      do n=1,NPROCY
        mythick=min(tljmax(ineighbor),thick-totthick)
        CALL MPI_RECV( data_sn_rcv, 8*LIMAX*ljmaxloc*Size1*Size2, MPI_BYTE,&
        ineighbor, msgnr, MPI_COMM_CALC, MPISTATUS, IERROR)
        do i2=1,Size2
          do jt=1,mythick
            do i=1,limax
              do i1=1,Size1
                data_north(i1,i,totthick+jt,i2)=data_sn_rcv(i1,i,jt,i2)
              end do
            end do
          end do
        end do
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor+NPROCX
        if(ineighbor>=NPROC)then
          do i2=1,Size2
            do jt=totthick+1,thick
              do i=1,limax
                do i1=1,Size1
                  data_north(i1,i,jt,i2) = data_north(i1,i,totthick,i2)
                end do
              end do
            end do
          end do
          exit
        endif
      enddo
    endif
  else
    do i2=1,Size2
      do jt=1,thick
        do i=1,limax
          do i1=1,Size1
            data_north(i1,i,jt,i2)=data(i1,i,ljmax,i2)
          end do
        end do
      end do
    end do
  end if
  
  jj=0
  do i2=1,Size2
    do jt=1,thick
      do it=1,myithick
        do i1=1,Size1
          data_west_snd(i1,it,jt-thick,i2)=data_south(i1,it,jt,i2)
          data_east_snd(i1,it,jt-thick,i2)=data_south(i1,limax-myithick+it,jt,i2)
        end do
      end do
    end do
    do j=1,ljmax
      do it=1,myithick
        do i1=1,Size1
          data_west_snd(i1,it,j,i2)=data(i1,it,j,i2)
          data_east_snd(i1,it,j,i2)=data(i1,limax-myithick+it,j,i2)
        end do
      end do
    end do
    do jt=1,thick
      do it=1,myithick
        do i1=1,Size1
          data_west_snd(i1,it,ljmax+jt,i2)=data_north(i1,it,jt,i2)
          data_east_snd(i1,it,ljmax+jt,i2)=data_north(i1,limax-myithick+it,jt,i2)
        end do
      end do
    end do
  end do
  
  if(neighbor(WEST) >= 0 )then
    
    if(thick<MINLIMAX)then
      CALL MPI_ISEND( data_west_snd , 8*(LJMAX+2*thick)*thick*Size1*Size2, MPI_BYTE,&
      neighbor(WEST), msgnr+3, MPI_COMM_CALC, request_w,IERROR)
    else
      ineighbor=neighbor(WEST)
      totthick=0
      do n=1,NPROCX
        mythick=min(tlimax(ineighbor),thick-totthick)
        CALL MPI_ISEND( data_west_snd , 8*(LJMAX+2*thick)*limaxloc*Size1*Size2, MPI_BYTE,&
        ineighbor, msgnr+3, MPI_COMM_CALC, irequest_w(n),IERROR)
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor-1
        if( (ineighbor/NPROCX)/=(me/NPROCX) .or. ineighbor<0) exit
      enddo
    endif
  end if
  if(neighbor(EAST) >= 0 )then
    if(thick<MINLIMAX)then
      CALL MPI_ISEND( data_east_snd , 8*(LJMAX+2*thick)*thick*Size1*Size2, MPI_BYTE,&
      neighbor(EAST), msgnr+7, MPI_COMM_CALC, request_e,IERROR)
    else
      ineighbor=neighbor(EAST)
      totthick=0
      do n=1,NPROCX
        mythick=min(tlimax(ineighbor),thick-totthick)
        CALL MPI_ISEND( data_east_snd , 8*(LJMAX+2*thick)*limaxloc*Size1*Size2, MPI_BYTE,&
        ineighbor, msgnr+7, MPI_COMM_CALC, irequest_e(n),IERROR)
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor+1
        if( (ineighbor/NPROCX)/=(me/NPROCX)) exit
      enddo
    endif
  end if
  
  
  
  if(neighbor(WEST) >= 0 )then
    if(thick<MINLIMAX)then
      CALL MPI_RECV( data_west, 8*(LJMAX+2*thick)*thick*Size1*Size2, MPI_BYTE,&
      neighbor(WEST), msgnr+7, MPI_COMM_CALC, MPISTATUS, IERROR)
    else
      ineighbor=neighbor(WEST)
      totthick=0
      do n=1,NPROCX
        mythick=min(tlimax(ineighbor),thick-totthick)
        CALL MPI_RECV( data_we_rcv, 8*(LJMAX+2*thick)*limaxloc*Size1*Size2, MPI_BYTE,&
        ineighbor, msgnr+7, MPI_COMM_CALC, MPISTATUS, IERROR)
        
        do i2=1,Size2
          do it=1-thick,LJMAX+thick
            do jt=1,mythick
              do i1=1,Size1
                data_west(i1,thick-totthick-mythick+jt,it,i2)=data_we_rcv(i1,tlimax(ineighbor)-mythick+jt,it,i2)
              end do
            end do
          end do
        end do
        totthick=totthick+mythick
        if(totthick>=thick)exit
        ineighbor=ineighbor-1
        if((ineighbor/NPROCX)/=(me/NPROCX).or. ineighbor<0)then
          do i2=1,Size2
            do it=1-thick,LJMAX+thick
              do jt=1,thick-totthick
                do i1=1,Size1
                  data_west(i1,jt,it,i2)=data_west(i1,thick-totthick+1,it,i2)
                end do
              end do
            end do
          end do
          exit
        endif
      enddo
    endif
  else
    
    do i2=1,Size2
      do jt=1,thick
        do it=1,thick
          do i1=1,Size1
            data_west(i1,it,jt-thick,i2)=data_south(i1,1,jt,i2)
          end do
        end do
      end do
      
      do j=1,ljmax
        do it=1,thick
          do i1=1,Size1
            data_west(i1,it,j,i2)=data(i1,1,j,i2)
          end do
        end do
      end do
      do jt=1,thick
        do it=1,thick
          do i1=1,Size1
            data_west(i1,it,ljmax+jt,i2)=data_north(i1,1,jt,i2)
          end do
        end do
      end do
      
    end do
  end if
  if(neighbor(EAST) >= 0 )then
    if(thick<MINLIMAX)then
      CALL MPI_RECV( data_east, 8*(LJMAX+2*thick)*thick*Size1*Size2, MPI_BYTE, &
      neighbor(EAST), msgnr+3, MPI_COMM_CALC, MPISTATUS, IERROR)
    else
      ineighbor=neighbor(EAST)
      totthick=0
      do n=1,NPROCX
        mythick=min(tlimax(ineighbor),thick-totthick)
        CALL MPI_RECV( data_we_rcv, 8*(LJMAX+2*thick)*limaxloc*Size1*Size2, MPI_BYTE,&
        ineighbor, msgnr+3, MPI_COMM_CALC, MPISTATUS, IERROR)
        
        do i2=1,Size2
          do it=1-thick,LJMAX+thick
            do jt=1,mythick
              do i1=1,Size1
                data_east(i1,totthick+jt,it,i2)=data_we_rcv(i1,jt,it,i2)
              end do
            end do
          end do
        end do
        totthick=totthick+mythick
        if(totthick>=thick)exit
        ineighbor=ineighbor+1
        if((ineighbor/NPROCX)/=(me/NPROCX))then
          do i2=1,Size2
            do it=1-thick,LJMAX+thick
              do jt=totthick+1,thick
                do i1=1,Size1
                  data_east(i1,jt,it,i2)=data_east(i1,totthick,it,i2)
                end do
              end do
            end do
          end do
          exit
        endif
      enddo
    endif
  else
    do i2=1,Size2
      do jt=1,thick
        do it=1,thick
          do i1=1,Size1
            data_east(i1,it,jt-thick,i2)=data_south(i1,limax,jt,i2)
          end do
        end do
      end do
      
      do j=1,ljmax
        do it=1,thick
          do i1=1,Size1
            data_east(i1,it,j,i2)=data(i1,limax,j,i2)
          end do
        end do
      end do
      
      do jt=1,thick
        do it=1,thick
          do i1=1,Size1
            data_east(i1,it,ljmax+jt,i2)=data_north(i1,limax,jt,i2)
          end do
        end do
      end do
    end do
  end if
! if(me==9)then
!    write(*,*)'DATA 9 EAST'
!    54 format(200I7)
!    do j=1-thick,ljmax+thick
!    write(*,54)(nint(data_east(1,i,j,1)),i=1,thick)
!    enddo
!    write(*,*)'DATA 9 WEST'
!    do j=1-thick,ljmax+thick
!    write(*,54)(nint(data_west(1,i,j,1)),i=1,thick)
!    enddo
! endif

  if(neighbor(SOUTH) >= 0 )then
    if(thick<MINLJMAX)then
      CALL MPI_WAIT(request_s, MPISTATUS,IERROR)
    else
      ineighbor=neighbor(SOUTH)
      totthick=0
      do n=1,NPROCY
        mythick=min(tljmax(ineighbor),thick-totthick)
        CALL MPI_WAIT(irequest_s(n), MPISTATUS,IERROR)
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor-NPROCX
        if(ineighbor<0) exit
      enddo
    endif
  end if
  if(neighbor(NORTH) >= 0 )then
    if(thick<MINLJMAX)then
      CALL MPI_WAIT(request_n, MPISTATUS,IERROR)
    else
      ineighbor=neighbor(NORTH)
      totthick=0
      do n=1,NPROCY
        mythick=min(tljmax(ineighbor),thick-totthick)
        CALL MPI_WAIT(irequest_n(n), MPISTATUS,IERROR)
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor+NPROCX
        if(ineighbor>=NPROC) exit
      enddo
    endif
  end if
  if(neighbor(WEST) >= 0 )then
    if(thick<MINLIMAX)then
      CALL MPI_WAIT(request_w, MPISTATUS, IERROR)
    else
      ineighbor=neighbor(WEST)
      totthick=0
      do n=1,NPROCX
        mythick=min(tlimax(ineighbor),thick-totthick)
        CALL MPI_WAIT(irequest_w(n), MPISTATUS, IERROR)
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor-1
        if( (ineighbor/NPROCX)/=(me/NPROCX).or. ineighbor<0) exit
      enddo
    end if
  end if
  if(neighbor(EAST) >= 0 )then
    if(thick<MINLIMAX)then
      CALL MPI_WAIT(request_e, MPISTATUS,IERROR)
    else
      ineighbor=neighbor(EAST)
      totthick=0
      do n=1,NPROCX
        mythick=min(tlimax(ineighbor),thick-totthick)
        CALL MPI_WAIT(irequest_e(n), MPISTATUS,IERROR)
        totthick=totthick+mythick
        if(totthick>=thick) exit
        ineighbor=ineighbor+1
        if( (ineighbor/NPROCX)/=(me/NPROCX)) exit
      enddo
    end if
  end if

end subroutine readneighbors_N

subroutine set_EuropeanAndGlobal_Config()
  !test if the grid covers Europe and the rest of the world
  !note that the two are not exclusive!
  !if the file covers more than just Europe, it is defined as GLOBAL
  !
  !If the file is GLOBAL:
  !1) USES%PFT_MAPS = .true. (can be overridden by FORCE_PFT_MAPS_FALSE)
  !2) USES%DEGREEDAY_FACTORS = .false.
  !3) USES%EURO_SOILNOX = .false.
  !4) USES%GLOBAL_SOILNOX = .true.
  !
  !If the file is EUROPEAN:
  !nothing happens for now


  !could not be in Config_module because "me" is not accessible there
  
  implicit none
  real:: x1,x2,x3,x4,y1,y2,y3,y4,lon,lat,ir,jr
  character(len=*), parameter :: dtxt='EurGlobSettings:'

  if(EUROPEAN_settings == 'NOTSET')then
     !No value set in config input, use grid to see if it covers Europe
     !Test approximatively if any European country is included in rundomain
     !in lon lat coordinates test if the middle of the rundomain is within
     ! 40>lon>-32    70>lat>35 

     if(gbacmax>35 .and. glacmin<40 .and. glacmax>-32)then
        
        if(MasterProc)write(*,*) dtxt//'assuming EUROPEAN_settings'
        EUROPEAN_settings = 'YES' 
     else
        
        ! define middle point of middle subdomain
        if(me==NPROC/2)then 
           lon = glon(limax/2,ljmax/2)
           lat = glat(limax/2,ljmax/2)
        endif
        CALL MPI_BCAST(lon,8,MPI_BYTE,NPROC/2,MPI_COMM_CALC,IERROR)
        CALL MPI_BCAST(lat,8,MPI_BYTE,NPROC/2,MPI_COMM_CALC,IERROR)
        
        x1=-32;x2=40;x3=x2;x4=x1;y1=35;y2=y1;y3=70;y4=y3
18      format(A,2F6.1,A)
        if(inside_1234(x1,x2,x3,x4,y1,y2,y3,y4,lon,lat))then
           EUROPEAN_settings = 'YES' 
           if(MasterProc)write(*,18) dtxt//'assuming EUROPEAN_settings: lon,lat ',lon,lat,' within Europe'
        else
           EUROPEAN_settings = 'NO' !default
           if(MasterProc)write(*,18) dtxt//'Not assuming EUROPEAN_settings: lon,lat ',lon,lat,' outside Europe'
        endif

     endif 
  else
    if(MasterProc) write(*,*) dtxt//'settings from config , EUR GLOB ', EUROPEAN_settings, GLOBAL_settings
  endif

  if(GLOBAL_settings == 'NOTSET')then
     !No value set in config input, use grid to see if it covers regions outside Extended Europe
     !We evaluate if the region extend outside the EMEP01 grid East, West or South
     !Note that it should also allow for PS projection which can cover the North Pole, i.e. any longitude
     GLOBAL_settings = 'NO' !default
     !find if lat < 19 are included within the domain
     if(gbacmin<19.0)then
        GLOBAL_settings = 'YES' !default
        if(MasterProc)write(*,*)dtxt//'Assuming GLOBAL_settings because rundomain extends below 19 degrees latitudes'
     else
        !find if the point with lon = -40 and lat = 45 is within the domain
        call lb2ij(-40.0,45.0,ir,jr)
        if(ir>=RUNDOMAIN(1).and.ir<=RUNDOMAIN(2).and.jr>=RUNDOMAIN(3).and.jr<=RUNDOMAIN(4))then
           GLOBAL_settings = 'YES' !default
           if(MasterProc)write(*,*)dtxt//'Assuming GLOBAL_settings because rundomain contains lon=-40 at lat=45'
        else 
           !find if the point with lon = 92 and lat = 45 is within the domain
           call lb2ij(92.0,45.0,ir,jr)
           if(ir>=RUNDOMAIN(1).and.ir<=RUNDOMAIN(2).and.jr>=RUNDOMAIN(3).and.jr<=RUNDOMAIN(4))then
              GLOBAL_settings = 'YES' !default
              if(MasterProc)write(*,*)dtxt//'Assuming GLOBAL_settings because rundomain contains lon=92 at lat=45'
           else
              if(MasterProc)write(*,*)dtxt//'Not assuming GLOBAL_settings'
           endif           
        endif
     endif
  endif

  if(EUROPEAN_settings == 'YES') then
     if(USES%MonthlyNH3  == 'NOTSET' .and. GLOBAL_settings /= 'YES')then
        USES%MonthlyNH3  = 'LOTOS'
        if(MasterProc)write(*,*)dtxt//' MonthlyNH3 set to '//trim(USES%MonthlyNH3)
     endif
  endif

  if(GLOBAL_settings == 'YES') then
     if(FORCE_PFT_MAPS_FALSE)then
        if(MasterProc)write(*,*)dtxt//'WARNING: NOT USING PFT_MAPS in a GLOBAL grid'
        USES%PFT_MAPS = .false. 
     else
        if(USES%PFT_MAPS)then        
           !nothing to change
        else
           if(MasterProc)write(*,*)dtxt//'Using PFT_MAPS because GLOBAL grid'
           USES%PFT_MAPS = .true. 
        endif
     endif
     if(USES%DEGREEDAY_FACTORS)then        
        if(MasterProc)write(*,*)dtxt//'WARNING: not using DEGREEDAY_FACTORS because GLOBAL grid'
        USES%DEGREEDAY_FACTORS = .false.
     endif

     if(.not. USES%GLOBAL_SOILNOX)then
        if(MasterProc)write(*,*)dtxt//'WARNING: setting GLOBAL_SOILNOX because GLOBAL grid'
        USES%GLOBAL_SOILNOX = .true.
     endif

     if(USES%EURO_SOILNOX)then        
        USES%EURO_SOILNOX = .false. 
        if(MasterProc)write(*,*)dtxt//'Not using EURO_SOILNOX because GLOBAL grid'
     endif

     if(USES%MonthlyNH3  == 'LOTOS')then
        if(MasterProc)write(*,*)dtxt//'WARNING: MonthlyNH3=LOTOS is not advised outside Europe'
     endif

  endif
  if(MasterProc)write(*,'(a,3(a,L))')dtxt//'Final settings, EUR = '//trim(EUROPEAN_settings)&
       //', GLOB = '//trim(GLOBAL_settings)&
       //', MonthlyNH3 = '//trim(USES%MonthlyNH3)&
       ,', E-SNOx = ',USES%EURO_SOILNOX&
       ,', G-SNOx = ',USES%GLOBAL_SOILNOX&
       ,', USE_SOILNOx= ',USE_SOILNOX

!  trim(EUROPEAN_settings)//','//&
!       trim(GLOBAL_settings), USES%EURO_SOILNOX, USES%GLOBAL_SOILNOX, USE_SOILNOX

end subroutine set_EuropeanAndGlobal_Config

end module GridValues_mod
