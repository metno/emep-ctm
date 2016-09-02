! <GridValues_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_10(3282)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2016 met.no
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
Module GridValues_ml
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Define parameters, variables and transformnations associated with grid and
!  projection.
!
! Nomenclature:
! fulldomain is the largest grid, usually where metdata is defined.
! rundomain is a grid where the run is performed, smaller than fulldomain.
! restricted domain is a grid smaller than rundomain, where data is outputed;
!  (the restricted domains are for instance, fullrun_DOMAIN,month_DOMAIN,
!  day_DOMAIN,hour_DOMAIN).
! subdomain: the domain covered by one MPI process or processor.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

use CheckStop_ml,           only: CheckStop,StopAll,check=>CheckNC
use Functions_ml,           only: great_circle_distance
use Io_Nums_ml,             only: IO_LOG,IO_TMP
use MetFields_ml 
use ModelConstants_ml,      only: &
     KMAX_BND, KMAX_MID, & ! vertical extent
     DEBUG,              & ! DEBUG%GRIDVALUES
     MasterProc,NPROC,IIFULLDOM,JJFULLDOM,RUNDOMAIN,&
     PT,Pref,NMET,METSTEP,USE_EtaCOORDINATES,MANUAL_GRID,USE_WRF_MET_NAMES,startdate
use MPI_Groups_ml    , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_REAL8, MPI_INTEGER, MPI_LOGICAL, &
                              MPI_MIN, MPI_MAX, MPI_SUM, &
                              MPI_COMM_CALC, MPI_COMM_WORLD, MPISTATUS, IERROR, ME_MPI, NPROC_MPI,&
                              ME_CALC, largeLIMAX,largeLJMAX 
use Par_ml, only : &
     LIMAX,LJMAX,  & ! max. possible i, j in this domain
     limax,ljmax,        & ! actual max.   i, j in this domain
     li0,li1,lj0,lj1,    & ! for debugging TAB
     GIMAX,GJMAX,        & ! Size of rundomain
     IRUNBEG,JRUNBEG,    & ! start of user-specified domain
     gi0,gj0,            & ! full-dom coordinates of domain lower l.h. corner
     gi1,gj1,            & ! full-dom coordinates of domain uppet r.h. corner
     me,                 & ! local processor
     parinit,parinit_groups
use PhysicalConstants_ml,     only: GRAV, PI, EARTH_RADIUS ! gravity, pi
use TimeDate_ml,              only: current_date,date,Init_nmdays,nmdays
use TimeDate_ExtraUtil_ml,    only: date2string
use InterpolationRoutines_ml, only: inside_1234
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

interface lb2ij
  module procedure lb2ij_real,lb2ij_int
end interface 
private :: lb2ij_real,lb2ij_int

public :: coord_check   ! normalize longitudes

public :: &
  coord_in_gridbox,  &  ! Are coord (lon/lat) inside gridbox(i,j)?
  coord_in_processor,&  ! Are coord (lon/lat) inside local domain?
  coord_in_domain       ! Are coord (lon/lat) inside "domain"?

public :: RestrictDomain !mask from full domain to rundomain

public :: GridRead!,Getgridparams
private :: Alloc_GridFields
private :: GetFullDomainSize
private :: find_poles

!** 1) Public (saved) Variables from module:

real, public, save :: &
     xp=0.0, yp=1.0,     & ! Coordinates of North pole (from infield)
     fi=0.0,         & ! projections rotation angle around y axis (from infield)
     AN=1.0,         & ! Distance on the map from pole to equator (No. of cells)
     GRIDWIDTH_M=1.0,& ! width of grid at 60N, in meters (old "h")(from infield)
     ref_latitude =60. ! latitude at which projection is true (degrees)

!Rotated_Spherical grid prarameters
real, public, save :: grid_north_pole_latitude,grid_north_pole_longitude,&
     dx_rot,dx_roti,x1_rot,y1_rot

!/ Variables to define full-domain (fdom) coordinates of local i,j values,
!  and reciprocal variables.
integer, public, allocatable, save, dimension(:) :: &
     i_fdom,j_fdom,&       ! fdom coordinates of local i,j
     i_local,j_local       ! local coordinates of full-domain i,j

!Parameters for Vertical Hybrid coordinates:
real, public, save,allocatable,  dimension(:) ::  &
     A_bnd,B_bnd,&         ! first [Pa],second [1] constants at layer boundary
                              ! (i.e. half levels in EC nomenclature)
     A_bnd_met,B_bnd_met,&         ! first [Pa],second [1] constants at layer boundary
                              ! (i.e. half levels in EC nomenclature)
     A_mid,B_mid,&         ! first [Pa],second [1] constants at middle of layer
                              ! (i.e. full levels in EC nomenclature)
     dA,dB,&               ! A_bnd(k+1)-A_bnd(k) [Pa],B_bnd(k+1)-B_bnd(k) [1]
                              ! P = A + B*PS; eta = A/Pref + B
     dEta_i,&              !1/deta = 1/(dA/Pref + dB)
     Eta_bnd,Eta_mid,&     ! boundary,midpoint of eta layer 
     sigma_bnd,sigma_mid   ! boundary,midpoint of sigma layer

real, public, save,allocatable,  dimension(:,:) :: &
     glon     ,glat    ,&  ! longitude,latitude of gridcell centers
     gl_stagg ,gb_stagg,&  ! longitude,latitude of gridcell corners 
                              !NB: gl_stagg, gb_stagg are here defined as the average of the four
                              !    surrounding gl gb.
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

character(len=230), public  :: filename_vert
integer, allocatable, save, public :: k1_met(:),k2_met(:)
real, allocatable, save, public :: x_k1_met(:)
logical, public, save ::  External_Levels_Def=.false.
integer, public, save :: KMAX_MET !number of vertical levels from the meteo files

contains

  subroutine GridRead(meteo,cyclicgrid)
    !   the subroutine reads the grid parameters (projection, resolution etc.)
    !   defined by the meteorological fields
    !
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
    nhour=0
    current_date = date(nyear, nmonth, nday, nhour, 0 )
    call Init_nmdays( current_date )

    !*********initialize grid parameters*********
    if(MANUAL_GRID)then
       !define the grid parameter manually (explicitely)
       if(me==0)write(*,*)'DEFINING GRID MANUALLY!'

       !must set all parameters... see example in version 3066 (or before)
    else


       !NOT MANUAL GRID


       !check first if grid is defined in a separate file:
       filename='Grid_Def.nc'
       inquire(file=filename,exist=Grid_Def_exist)
       Grid_Def_exist=Grid_Def_exist.and.Use_Grid_Def
       if(Grid_Def_exist)then
          if(MasterProc)write(*,*)'Found Grid_Def! ',trim(filename)
       else
          if(MasterProc.and.Use_Grid_Def)write(*,*)'Did not found Grid_Def ',trim(filename)
          !56 FORMAT(a5,i4.4,i2.2,i2.2,a3)
          !   write(filename,56)'meteo',nyear,nmonth,nday,'.nc'
          filename = date2string(meteo,startdate)
       endif
       if(MasterProc)write(*,*)'reading domain sizes from ',trim(filename)

       call GetFullDomainSize(filename,IIFULLDOM,JJFULLDOM,KMAX_MET,METSTEP,projection)

       KMAX_MID=0!initialize
       filename_vert='Vertical_levels.txt'
       open(IO_TMP,file=filename_vert,action="read",iostat=ios)
       if(ios==0)then
          !define own vertical coordinates
          if(me==0)write(*,*)'Define vertical levels from ',trim(filename_vert)
          read(IO_TMP,*)KMAX_MID
          if(me==0)write(*,*)KMAX_MID, 'vertical levels '
          External_Levels_Def=.true.
          !Must use eta coordinates
          if(.not.USE_EtaCOORDINATES)write(*,*)'WARNING: using hybrid levels even if not asked to! '
          USE_EtaCOORDINATES=.true.   
       else
          External_Levels_Def=.false.
          close(IO_TMP)
          KMAX_MID=KMAX_MET
       endif

       KMAX_BND=KMAX_MID+1

       allocate(A_bnd(KMAX_BND),B_bnd(KMAX_BND))
       allocate(A_mid(KMAX_MID),B_mid(KMAX_MID))
       allocate(dA(KMAX_MID),dB(KMAX_MID),dEta_i(KMAX_MID))
       allocate(sigma_bnd(KMAX_BND),sigma_mid(KMAX_MID))
       allocate(Eta_bnd(KMAX_BND),Eta_mid(KMAX_MID))

       allocate(i_local(IIFULLDOM))
       allocate(j_local(JJFULLDOM))

       !set RUNDOMAIN default values where not defined
       if(RUNDOMAIN(1)<1)RUNDOMAIN(1)=1
       if(RUNDOMAIN(2)<1 .or. RUNDOMAIN(2)>IIFULLDOM) RUNDOMAIN(2)=IIFULLDOM
       if(RUNDOMAIN(3)<1)RUNDOMAIN(3)=1
       if(RUNDOMAIN(4)<1 .or. RUNDOMAIN(4)>JJFULLDOM) RUNDOMAIN(4)=JJFULLDOM
       if(MasterProc)then
55        format(A,I5,A,I5)
          write(*,55)     'FULLDOMAIN has sizes ',IIFULLDOM,' X ',JJFULLDOM
          write(IO_LOG,55)'FULLDOMAIN has sizes ',IIFULLDOM,' X ',JJFULLDOM
          write(*,55)     'RUNDOMAIN  x coordinates from ',RUNDOMAIN(1),' to ',RUNDOMAIN(2)
          write(IO_LOG,55)'RUNDOMAIN  x coordinates from ',RUNDOMAIN(1),' to ',RUNDOMAIN(2)
          write(*,55)     'RUNDOMAIN  y coordinates from ',RUNDOMAIN(3),' to ',RUNDOMAIN(4)
          write(IO_LOG,55)'RUNDOMAIN  y coordinates from ',RUNDOMAIN(3),' to ',RUNDOMAIN(4)
       endif

       call find_poles(filename,Pole_Singular)

       MIN_GRIDS=5
       if(NPROC==NPROC_MPI)then
          !partition into subdomains
          call parinit(MIN_GRIDS,Pole_Singular)     !subdomains sizes and position
       else
          !partition into largesubdomains and subdomains
          call parinit_groups(MIN_GRIDS,Pole_Singular)     !subdomains sizes and position
       endif

       call Alloc_MetFields(LIMAX,LJMAX,KMAX_MID,KMAX_BND,NMET)
       
       if(ME_CALC>=0)then
          call Alloc_GridFields(LIMAX,LJMAX,KMAX_MID,KMAX_BND)
       else
          call Alloc_GridFields(largeLIMAX,largeLJMAX,KMAX_MID,KMAX_BND)
       endif

       if(ME_CALC>=0)then
          call Getgridparams(LIMAX,LJMAX,filename,cyclicgrid)
          !defines i_fdom,j_fdom,i_local,j_local,Cyclicgrid,North_pole,Poles
          !GRIDWIDTH_M, glon, glat, xm_i,xm_j,xm2,xmd,xmdji,xm2ji,gl_stagg,gb_stagg
          !for Stereographic projection: ref_latitude,fi,xp,yp,AN
          !for lon lat projection: no additional parameters
          !for Rotated_Spherical projection: grid_north_pole_latitude,grid_north_pole_longitude,x1_rot,y1_rot,dx_rot,dx_roti
          !P0,A_bnd_met,B_bnd_met,A_bnd,B_bnd,A_mid,B_mid,sigma_mid,sigma_bnd
      else
          call Getgridparams(largeLIMAX,largeLJMAX,filename,cyclicgrid)
     endif


       if(ios==0)close(IO_TMP)

       if(MasterProc .and. DEBUG%GRIDVALUES)then
          write(*,*)'sigma_mid:',(sigma_mid(k),k=1,20)
          write(*,*)'grid resolution:',GRIDWIDTH_M
          write(*,*)'xcoordinate of North Pole, xp:',xp
          write(*,*)'ycoordinate of North Pole, yp:',yp
          write(*,*)'longitude rotation of grid, fi:',fi
          write(*,*)'true distances latitude, ref_latitude:',ref_latitude
       endif

    endif

!
  end subroutine GridRead



  subroutine GetFullDomainSize(filename,IIFULLDOM,JJFULLDOM,KMAX,METSTEP,projection)

    !
    ! Get input grid sizes 
    !

    implicit none

    character (len = *), intent(in) ::filename
    integer, intent(out):: IIFULLDOM,JJFULLDOM,KMAX,METSTEP
    character (len = *), intent(out) ::projection

    integer :: status,ncFileID,idimID,jdimID, kdimID,timeDimID
    integer :: GIMAX_file,GJMAX_file,KMAX_file,wrf_proj_code
    real :: wrf_POLE_LAT=0.0


    if(ME_MPI==0)then
       print *,'Defining grid properties from ',trim(filename)
       !open an existing netcdf dataset
       status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
       if(status /= nf90_noerr) then
          print *,'not found',trim(filename)
          call StopAll("GridValues: File not found")
       endif

       !          print *,'  reading ',trim(filename)
       projection=''
       status = nf90_get_att(ncFileID,nf90_global,"projection",projection)
       if(status /= nf90_noerr) then
          !WRF projection format
          call check(nf90_get_att(ncFileID,nf90_global,"MAP_PROJ",wrf_proj_code))
          if(.not.USE_WRF_MET_NAMES .and. me==0)write(*,*)'Assuming WRF metdata'
          USE_WRF_MET_NAMES = .true.
          if(wrf_proj_code==6)then
             status = nf90_get_att(ncFileID,nf90_global,"POLE_LAT",wrf_POLE_LAT)
             if(status == nf90_noerr) then
                write(*,*)"POLE_LAT", wrf_POLE_LAT
                if(abs(wrf_POLE_LAT-90.0)<0.001)then
                   projection='lon lat'
                else
                   projection='Rotated_Spherical'                                
                endif
             else
                write(*,*)"POLE_LAT not found"
                projection='lon lat'
             endif
          else if(wrf_proj_code==2)then
             projection='Stereographic'     
          else
             call StopAll("Projection not recognized")  
          endif
       endif

       !put into emep standard
       if(trim(projection)=='Polar Stereographic')projection='Stereographic'

       if(trim(projection)=='Rotated_Spherical'.or.trim(projection)=='rotated_spherical'&
            .or.trim(projection)=='rotated_pole'.or.trim(projection)=='rotated_latitude_longitude')then
          projection='Rotated_Spherical'
       endif

       write(*,*)'projection: ',trim(projection)

       !get dimensions id
       if(trim(projection)=='Stereographic') then
          status = nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID)
          if(status /= nf90_noerr) then
             !WRF  format
             call check(nf90_inq_dimid(ncid = ncFileID, name = "west_east", dimID = idimID))
          endif
          status = nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID)
          if(status /= nf90_noerr) then
             !WRF  format
             call check(nf90_inq_dimid(ncid = ncFileID, name = "south_north", dimID = jdimID))
          endif
       elseif(trim(projection)==trim('lon lat')) then
          status=nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID)
          if(status /= nf90_noerr) then
             !WRF format
             call check(nf90_inq_dimid(ncid = ncFileID, name = "west_east", dimID = idimID))
          endif
          status=nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID)
          if(status /= nf90_noerr) then
             !WRF format
             call check(nf90_inq_dimid(ncid = ncFileID, name = "south_north", dimID = jdimID))
          endif
       else
          !     write(*,*)'GENERAL PROJECTION ',trim(projection)
          status=nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID)
          if(status /= nf90_noerr) then
             !WRF  format
             call check(nf90_inq_dimid(ncid = ncFileID, name = "west_east", dimID = idimID))
          endif
          status=nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID)
          if(status /= nf90_noerr) then
             !WRF  format
             call check(nf90_inq_dimid(ncid = ncFileID, name = "south_north", dimID = jdimID))
          endif
       endif

       status=nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID)
       if(status /= nf90_noerr)then
          status=nf90_inq_dimid(ncid = ncFileID, name = "lev", dimID = kdimID)!hybrid coordinates
          if(status /= nf90_noerr) then
             !WRF  format
             call check(nf90_inq_dimid(ncid = ncFileID, name = "bottom_top", dimID = kdimID))      
          endif
       endif

       !get dimensions length
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_file))
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_file))
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_file))

       write(*,*)'dimensions input grid:',GIMAX_file,GJMAX_file,KMAX_file!,Nhh

       IIFULLDOM=GIMAX_file
       JJFULLDOM=GJMAX_file
       KMAX     =KMAX_file


       !find METSTEP (checked also in first meteo read)
       status=nf90_inq_dimid(ncid=ncFileID,name="time",dimID=timedimID)
       if(status/=nf90_noerr)then
          status=nf90_inq_dimid(ncid=ncFileID,name="Time",dimID=timedimID)! WRF
       endif
       if(status/=nf90_noerr)then
          write(*,*)'time variable not found assuming 8 records'
          Nhh=8
       else
          call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=Nhh))
       endif

       METSTEP=24/Nhh
       write(*,*)'METSTEP set to ',METSTEP,' hours'
       call check(nf90_close(ncFileID))       
    endif

    CALL MPI_BCAST(METSTEP ,4,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_BCAST(USE_WRF_MET_NAMES ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_BCAST(IIFULLDOM ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_BCAST(JJFULLDOM ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_BCAST(KMAX ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_BCAST(projection ,len(projection),MPI_BYTE,0,MPI_COMM_WORLD,IERROR)

  end subroutine GetFullDomainSize

  subroutine find_poles(filename,Pole_Singular)
    !defines if there is a singularity at poles

    implicit none
    character (len = *), intent(in) ::filename
    integer, intent(out):: Pole_Singular
    integer :: status,ncFileID,varid
    real,allocatable :: latitudes(:)

    Pole_Singular=0
    if(trim(projection)==trim('lon lat')) then
       if(ME_MPI==0)then
          !find wether poles are included (or almost included) in grid
          !
          !If some cells are to narrow (Poles in lat lon coordinates),
          !this will give too small time steps in the Advection,
          !because of the constraint that the Courant number should be <1.
          !
          !If Poles are found and lon-lat coordinates are used the Advection scheme
          !will be modified to be able to cope with the singularity
          !the advection routine will not work efficiently with NPROCY>2 in this case
          
          !open an existing netcdf dataset
          status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
          if(status /= nf90_noerr) then
             print *,'not found',trim(filename)
             call StopAll("GridValues: File not found")
          endif
          
          allocate(latitudes(JJFULLDOM))
          status=nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID)
          if(status /= nf90_noerr) then
             !WRF format
             call check(nf90_inq_varid(ncid = ncFileID, name = "XLAT", varID = varID))
          endif
          call check(nf90_get_var(ncFileID, varID,latitudes  ))
          if(latitudes(RUNDOMAIN(4))>88.0)then
             write(*,*)'The grid is singular at North Pole'
             Pole_Singular=Pole_Singular+1
          endif
          if(latitudes(RUNDOMAIN(3))<-88.0)then
             write(*,*)'The grid is singular at South Pole'
             Pole_Singular=Pole_Singular+1
          endif
          deallocate(latitudes)
          call check(nf90_close(ncFileID))       
       endif
    endif

    CALL MPI_BCAST(Pole_Singular ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)

  end subroutine find_poles

  subroutine Getgridparams(LIMAX,LJMAX,filename,cyclicgrid)
    !
    ! Get grid and time parameters as defined in the meteo or Grid_Def file
    ! Do some checks on sizes and dates
    !
    ! This routine is called only once (and is therefore not optimized for speed)
    !
        !defines i_fdom,j_fdom,i_local,j_local,Cyclicgrid,North_pole,Poles
       !GRIDWIDTH_M, glon, glat, xm_i,xm_j,xm2,xmd,xmdji,xm2ji,gl_stagg,gb_stagg
       !for Stereographic projection: ref_latitude,fi,xp,yp,AN
       !for lon lat projection: no additional parameters
       !for Rotated_Spherical projection: grid_north_pole_latitude,grid_north_pole_longitude,x1_rot,y1_rot,dx_rot,dx_roti
       !P0,A_bnd_met,B_bnd_met,A_bnd,B_bnd,A_mid,B_mid,sigma_mid,sigma_bnd
   

    implicit none

    integer, intent(in):: LIMAX,LJMAX
    character (len = *), intent(in) ::filename
    integer, intent(out):: cyclicgrid

    integer :: n,i,j,k,kk
    integer :: ncFileID,idimID,jdimID,varID
    integer :: status,South_pole,North_pole
    real :: x1,x2,x3,x4,P0,x,y,mpi_out
    logical::found_hybrid=.false.
    real :: CEN_LAT, CEN_LON,P_TOP_MET
    real :: rb,rl,rp,dx,dy,dy2,glmax,glmin,v2(2),glon_fdom1,glat_fdom1
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
    if(status /= nf90_noerr) then
       print *,'not found',trim(filename)
       call StopAll("GridValues: File not found")
    endif
    if(MasterProc)print *,'Defining grid parameters from ',trim(filename)

    !get dimensions id
    if(trim(projection)=='Stereographic') then
       status = nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_inq_dimid(ncid = ncFileID, name = "west_east", dimID = idimID))
       endif
       status = nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_inq_dimid(ncid = ncFileID, name = "south_north", dimID = jdimID))
       endif
    elseif(trim(projection)==trim('lon lat')) then
       status=nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_inq_dimid(ncid = ncFileID, name = "west_east", dimID = idimID))
       endif
       status=nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_inq_dimid(ncid = ncFileID, name = "south_north", dimID = jdimID))
       endif
     else
       !     write(*,*)'GENERAL PROJECTION ',trim(projection)
       status=nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_inq_dimid(ncid = ncFileID, name = "west_east", dimID = idimID))
       endif
       status=nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_inq_dimid(ncid = ncFileID, name = "south_north", dimID = jdimID))
       endif
    endif

    !get global attributes
    status = nf90_get_att(ncFileID,nf90_global,"Grid_resolution",GRIDWIDTH_M)
    if(status /= nf90_noerr) then
       !WRF  format
       call check(nf90_get_att(ncFileID,nf90_global,"DX",GRIDWIDTH_M))
    endif
    if(MasterProc)write(*,*)"Grid_resolution",GRIDWIDTH_M

    if(trim(projection)=='Stereographic')then
       status = nf90_get_att(ncFileID,nf90_global,"ref_latitude",ref_latitude)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_get_att(ncFileID,nf90_global,"TRUELAT1",ref_latitude))
       endif
       status = nf90_get_att(ncFileID, nf90_global, "fi",fi )
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_get_att(ncFileID, nf90_global, "STAND_LON",fi ))
       endif
       status = nf90_get_att(ncFileID, nf90_global, "xcoordinate_NorthPole",xp )
       if(status == nf90_noerr) then
          call check(nf90_get_att(ncFileID, nf90_global, "ycoordinate_NorthPole",yp ))
       else
          !WRF  format compute from grid center coordinates
          call check(nf90_get_att(ncFileID, nf90_global, "CEN_LAT" &
               , CEN_LAT))
          call check(nf90_get_att(ncFileID, nf90_global, "CEN_LON" &
               , CEN_LON))
          xp = 0.5 + 0.5*IIFULLDOM - EARTH_RADIUS/GRIDWIDTH_M*(1+sin(ref_latitude*PI/180.))&
               *tan(PI/4-CEN_LAT*PI/180./2)*sin((CEN_LON-fi)*PI/180.)
          yp = 0.5 + 0.5*JJFULLDOM + EARTH_RADIUS/GRIDWIDTH_M*(1+sin(ref_latitude*PI/180.))&
               *tan(PI/4-CEN_LAT*PI/180./2)*cos((CEN_LON-fi)*PI/180.)
          !correct for last digits. Assume that numbers are close enough to an integer
          if(abs(nint(xp)-xp)<0.01)xp=nint(xp)
          if(abs(nint(yp)-yp)<0.01)yp=nint(yp)

          if(MasterProc)write(*,*)"M= ",EARTH_RADIUS/GRIDWIDTH_M*(1+sin(ref_latitude*PI/180.))
          if(MasterProc)write(*,*)"coordinates of North pole ",xp,yp
       endif

       AN = 6.370e6*(1.0+sin( ref_latitude*PI/180.))/GRIDWIDTH_M ! = 237.7316364 for GRIDWIDTH_M=50 km and ref_latitude=60

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

    elseif(trim(projection)==trim('lon lat')) then
       if(.not. USE_WRF_MET_NAMES)then
          !NB: lon and lat are stored as 1 dimensional arrays
          call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
          
          call check(nf90_get_var(ncFileID, varID, lon_ext(1:limax,1),start=(/gi0+IRUNBEG-1/),count=(/limax/) ))
          if(LIMAX>limax)lon_ext(LIMAX,1)=lon_ext(limax,1)+(lon_ext(limax,1)-lon_ext(limax-1,1))
          lon_ext(0,1)=2*lon_ext(1,1)-lon_ext(2,1)
          lon_ext(LIMAX+1,1)=2*lon_ext(LIMAX,1)-lon_ext(LIMAX-1,1)
          do j=0,LJMAX+1
             lon_ext(:,j)=lon_ext(:,1)
          enddo
          
          call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
          call check(nf90_get_var(ncFileID, varID, lat_ext(1,1:ljmax),start=(/gj0+JRUNBEG-1/),count=(/ljmax/) ))
          lat_ext(1,LJMAX)=min(90.0,lat_ext(1,LJMAX))!should never be used anyway
          lat_ext(1,0)=2*lat_ext(1,1)-lat_ext(1,2)
          lat_ext(1,LJMAX+1)=2*lat_ext(1,LJMAX)-lat_ext(1,LJMAX-1)
          do i=0,LIMAX+1
             lat_ext(i,:)=lat_ext(1,:)
          enddo
       else
          !WRF  format
          call check(nf90_inq_varid(ncid = ncFileID, name = "XLONG", varID = varID))
          call nf90_get_var_extended(ncFileID,varID,lon_ext,0,LIMAX+1,0,LJMAX+1)
          call check(nf90_inq_varid(ncid = ncFileID, name = "XLAT", varID = varID))
          call nf90_get_var_extended(ncFileID,varID,lat_ext,0,LIMAX+1,0,LJMAX+1)
       endif
    else if(trim(projection)=='Rotated_Spherical')then
       status=nf90_get_att(ncFileID,nf90_global,"grid_north_pole_latitude",grid_north_pole_latitude)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_get_att(ncFileID,nf90_global,"POLE_LAT",grid_north_pole_latitude))
       endif
       if(MasterProc)write(*,*)"grid_north_pole_latitude",grid_north_pole_latitude
       status=nf90_get_att(ncFileID,nf90_global,"grid_north_pole_longitude",grid_north_pole_longitude)
       if(status /= nf90_noerr) then
          !WRF  format
          call check(nf90_get_att(ncFileID,nf90_global,"POLE_LON",grid_north_pole_longitude))
          !find resolution in degrees from resolution in km. WRF uses Erath Radius 6370 km(?)
          dx_rot=360./(6370000.*2*PI/GRIDWIDTH_M)
          !round to 6 digits
          dx_rot=0.000001*nint(1000000*dx_rot)
       endif
       if(MasterProc)write(*,*)"grid_north_pole_longitude",grid_north_pole_longitude
       status=nf90_inq_varid(ncid = ncFileID, name = "i", varID = varID)
       if(status == nf90_noerr) then
          call check(nf90_get_var(ncFileID, varID, v2))!note that i is one dimensional
          x1_rot=v2(1)
          dx_rot=v2(2)-v2(1)
          call check(nf90_inq_varid(ncid = ncFileID, name = "j", varID = varID))
          call check(nf90_get_var(ncFileID, varID, v2(1)))!note that j is one dimensional
          y1_rot=v2(1)
          call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
          call nf90_get_var_extended(ncFileID,varID,lon_ext,0,LIMAX+1,0,LJMAX+1)
          call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
          call nf90_get_var_extended(ncFileID,varID,lat_ext,0,LIMAX+1,0,LJMAX+1)
       else
          !WRF  format
          call check(nf90_inq_varid(ncid = ncFileID, name = "XLONG", varID = varID))
          call nf90_get_var_extended(ncFileID,varID,lon_ext,0,LIMAX+1,0,LJMAX+1)
          call check(nf90_get_var(ncFileID, varID, v2,start=(/1,1/),count=(/1,1/)  ))
          glon_fdom1=v2(1)
          !glon=0.0!to get some value for outside subdomain too (when limax<LIMAX for instance)
          !call check(nf90_get_var(ncFileID, varID, glon(1:limax,1:ljmax),&
          !     start=(/gi0+IRUNBEG-1,gj0+JRUNBEG-1/),count=(/limax,ljmax/)  ))
          call check(nf90_inq_varid(ncid = ncFileID, name = "XLAT", varID = varID))
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
             endif
          enddo
          !                call lb2ij(glon_fdom(1,1),glat_fdom(1,1),x,y)
          !                write(*,*)'after ',glon_fdom(1,1),glat_fdom(1,1),x,y
          !                call lb_rot2lb(x,y,x1_rot,y1_rot,grid_north_pole_longitude,grid_north_pole_latitude)
          !                write(*,*)"spherical lon lat of (i,j)=(1,1)",x,y,glon_fdom(1,1),glat_fdom(1,1)
          if(MasterProc)write(*,*)"rotated lon lat of (i,j)=(1,1)",x1_rot,y1_rot
          if(MasterProc)write(*,*)"resolution",dx_rot
       endif
       dx_roti=1.0/dx_rot

    else
       !other projection?
       call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
       call nf90_get_var_extended(ncFileID,varID,lon_ext,0,LIMAX+1,0,LJMAX+1)
       call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
       call nf90_get_var_extended(ncFileID,varID,lat_ext,0,LIMAX+1,0,LJMAX+1)
    endif

    glon(1:LIMAX,1:LJMAX)=lon_ext(1:LIMAX,1:LJMAX)             ! longitude
    glat(1:LIMAX,1:LJMAX)=lat_ext(1:LIMAX,1:LJMAX)             ! latitude
    do j=1,LJMAX
       do i=1,LIMAX
          if(glon(i,j)>glmax)glon(i,j)=glon(i,j)-360.0
          if(glon(i,j)<glmin)glon(i,j)=glon(i,j)+360.0
       enddo
    enddo

    !map factors
    status=nf90_inq_varid(ncid=ncFileID, name="map_factor", varID=varID)

    if(status == nf90_noerr)then
       !mapping factor at center of cells is defined

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
          enddo
       enddo

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

       if(status == nf90_noerr)then
          call nf90_get_var_extended(ncFileID,varID,xm_i_ext,-1,LIMAX+2,-1,LJMAX+2)
       else
          !WRF  format
          call check(nf90_inq_varid(ncid=ncFileID, name="MAPFAC_VX", varID=varID))             
          call nf90_get_var_extended(ncFileID,varID,xm_i_ext,-1,LIMAX+2,-1,LJMAX+2,jshift_in=1)!NB:shift j by 1 since wrf start at bottom face
       endif

       status=nf90_inq_varid(ncid=ncFileID, name="map_factor_j", varID=varID)
       if(status == nf90_noerr)then
          call nf90_get_var_extended(ncFileID,varID,xm_j_ext,-1,LIMAX+2,-1,LJMAX+2)
       else
          !WRF  format
          call check(nf90_inq_varid(ncid=ncFileID, name="MAPFAC_UY", varID=varID))
          call nf90_get_var_extended(ncFileID,varID,xm_j_ext,-1,LIMAX+2,-1,LJMAX+2,ishift_in=1)!NB:shift i by 1 since wrf start at left face
       endif

       !define xm2, xm_i and xm_j now
       !Note that xm is inverse length: interpolate 1/xm rather than xm
       do j=0,LJMAX+1
          do i=0,LIMAX+1
             xm_i(i,j)=xm_i_ext(i,j)
             xm_j(i,j)=xm_j_ext(i,j)
             xm2(i,j) = 4.0*( (xm_i_ext(i,j-1)*xm_i_ext(i,j))/&
                  (xm_i_ext(i,j-1)+xm_i_ext(i,j))  )&
                  *( (xm_j_ext(i-1,j)*xm_j_ext(i,j))/&
                  (xm_j_ext(i-1,j)+xm_j_ext(i,j))  )
             xmd(i,j) =1.0/xm2(i,j)
             xm2ji(j,i) = xm2(i,j)
             xmdji(j,i) = xmd(i,j)                
          enddo
       enddo

    endif

    status=nf90_inq_varid(ncid = ncFileID, name = "k", varID = varID)
    if(status /= nf90_noerr)then
       !always use hybrid coordinates at output, if hybrid in input
       if(.not.USE_EtaCOORDINATES)then
          write(*,*)'WARNING: using hybrid levels even if not asked to! ',trim(filename)
          USE_EtaCOORDINATES=.true.
       endif
       if(MasterProc)write(*,*)'reading met hybrid levels from ',trim(filename)
       !          call check(nf90_inq_varid(ncid = ncFileID, name = "hyam", varID = varID))                 
       !          call check(nf90_get_var(ncFileID, varID, A_mid ))
       !          A_mid=P0*A_mid!different definition in modell and grid_Def
       !          call check(nf90_inq_varid(ncid = ncFileID, name = "hybm", varID = varID))                 
       !          call check(nf90_get_var(ncFileID, varID,B_mid))
       status=nf90_inq_varid(ncid = ncFileID, name = "P0", varID = varID) 
       if(status /= nf90_noerr)then
          status=nf90_inq_varid(ncid = ncFileID, name = "P00", varID = varID) !WRF case
          if(status /= nf90_noerr)then
             if(External_Levels_Def)then
                write(*,*)'WARNING: did not find P0. Assuming vertical levels from ',trim(filename_vert)
             else
                write(*,*)'Do not know how to define vertical levels '
                call StopAll('Define levels in Vertical_levels.txt')
             endif
          else
             !WRF
             !asuming sigma levels ZNW=(P-P_TOP_MET)/(PS-P_TOP_MET)
             !P = A+B*PS = P_TOP_MET*(1-ZNW) + ZNW*PS  
             !B = ZNW
             !A = P_TOP_MET*(1-ZNW)
             call check(nf90_get_var(ncFileID, varID, P0 ))
             if(.not.allocated(A_bnd_met))allocate(A_bnd_met(KMAX_MET+1),B_bnd_met(KMAX_MET+1))
             call check(nf90_inq_varid(ncid = ncFileID, name = "P_TOP", varID = varID))                 
             call check(nf90_get_var(ncFileID, varID, P_TOP_MET ))
             call check(nf90_inq_varid(ncid = ncFileID, name = "ZNW", varID = varID))       
             call check(nf90_get_var(ncFileID, varID, B_bnd_met ))
             if(MET_REVERSE_K)then
                A_bnd_met=B_bnd_met!use A_bnd_met as temporary buffer
                do k=1,KMAX_MET+1
                   B_bnd_met(k)=A_bnd_met(KMAX_MET+2-k)
                enddo
             endif
             A_bnd_met=P_TOP_MET*(1.-B_bnd_met)
          endif
          if(MET_REVERSE_K)then
             if(MasterProc)write(*,*)"Reversed vertical levels from met, P at levels boundaries:"
          else
             if(MasterProc)write(*,*)"Vertical levels from met, P at levels boundaries:"
          endif
          do k=1,KMAX_MET+1
             if(MasterProc)write(*,44)k, A_bnd_met(k)+P0*B_bnd_met(k)
          enddo
       else
          call check(nf90_get_var(ncFileID, varID, P0 ))
          if(.not.allocated(A_bnd_met))allocate(A_bnd_met(KMAX_MET+1),B_bnd_met(KMAX_MET+1))
          call check(nf90_inq_varid(ncid = ncFileID, name = "hyai", varID = varID))                 
          call check(nf90_get_var(ncFileID, varID, A_bnd_met ))
          A_bnd_met=P0*A_bnd_met!different definition in model and grid_Def
          call check(nf90_inq_varid(ncid = ncFileID, name = "hybi", varID = varID))                 
          call check(nf90_get_var(ncFileID, varID, B_bnd_met ))          
       endif
       if(External_Levels_Def)then
          !model levels defined from external text file
          if(MasterProc)write(*,*)'reading external hybrid levels from ',trim(filename_vert)
          P0=Pref
          do k=1,KMAX_MID+1
             read(IO_TMP,*)kk,A_bnd(k),B_bnd(k)
             if(kk/=k.and.MasterProc)write(*,*)'WARNING: unexpected format for vertical levels ',k,kk
          enddo
          if(status /= nf90_noerr)then
             !assume levels from metdata are defined in filename_vert
             if(.not.allocated(A_bnd_met))allocate(A_bnd_met(KMAX_MET+1),B_bnd_met(KMAX_MET+1))
             A_bnd_met=A_bnd
             B_bnd_met=B_bnd
          endif
       else
          !vertical model levels are the same as in meteo 
          A_bnd=A_bnd_met
          B_bnd=B_bnd_met
       endif

       do k=1,KMAX_MID
          A_mid(k)=0.5*(A_bnd(k)+A_bnd(k+1))
          B_mid(k)=0.5*(B_bnd(k)+B_bnd(k+1))
       enddo
       sigma_mid =B_mid!for Hybrid coordinates sigma_mid=B if A*P0=PT-sigma_mid*PT

       if(me==0)write(*,*)"Hybrid vertical coordinates, P at levels boundaries:"
       do k=1,KMAX_MID+1
44        FORMAT(i4,10F12.2)
          if(me==0)write(*,44)k, A_bnd(k)+P0*B_bnd(k)
       enddo
       !test if the top is within the height defined in the meteo files
       if(me==0.and.External_Levels_Def.and.(A_bnd(1)+P0*B_bnd(1)<A_bnd_met(1)+P0*B_bnd_met(1)))then
          write(*,*)'Pressure at top of defined levels is ',A_bnd(1)+P0*B_bnd(1)
          write(*,*)'Pressure at top defined in meteo files is ',A_bnd_met(1)+P0*B_bnd_met(1)
          write(*,*)'Pressure at op must be higher (lower altitude) than top defined in meteo '
          call StopAll('Top level too high! Change values in Vertical_levels.txt')
       endif

       !test if the levels can cope with highest mountains (400 hPa)
       do k=1,KMAX_MID
          if(me==0.and.A_bnd(k+1)+40000*B_bnd(k+1)-(A_bnd(k)+40000*B_bnd(k))<0.0)then
             write(*,*)'WARNING: hybrid vertical level definition may cause negative level thickness when pressure below 400 hPa '
             write(*,*)'Pressure at level ',k,' is ',A_bnd(k)+40000*B_bnd(k)
             write(*,*)'Pressure at level ',k+1,' is ',A_bnd(k+1)+40000*B_bnd(k+1),' (should be higher)'
             if(External_Levels_Def)call StopAll('GridValues_ml: possible negative level thickness ')
          endif
       enddo

       !test if the lowest levels is thick enough (twice height of highest vegetation?) about 550 Pa = about 46m 
       !Deposition scheme is not designes for very thin lowest levels
       if(me==0.and.A_bnd(KMAX_MID+1)+P0*B_bnd(KMAX_MID+1)-(A_bnd(KMAX_MID)+P0*B_bnd(KMAX_MID))<550.0)then
          write(*,*)'WARNING: lowest level very shallow; ',A_bnd(KMAX_MID+1)+P0*B_bnd(KMAX_MID+1) -&
               (A_bnd(KMAX_MID)+P0*B_bnd(KMAX_MID)),'Pa'
          call StopAll('Lowest level too thin! Change vertical levels definition in Vertical_levels.txt ')
       endif

       found_hybrid=.true.
    else
       call check(nf90_get_var(ncFileID, varID, sigma_mid ))
    endif
    call check(nf90_close(ncFileID))

    !TEMPORARY: definition of sigma. Only A and B will be used in the future
    ! definition of the half-sigma levels (boundaries between layers)
    ! from the full levels. 
    sigma_bnd(KMAX_BND) = 1.
    do k = KMAX_MID,2,-1
       sigma_bnd(k) = 2.*sigma_mid(k) - sigma_bnd(k+1)
    enddo
    sigma_bnd(1) = 0.

    if(.not.(found_hybrid.or.Grid_Def_exist))then
       !define A and B that gives the correspondin sigma
       do k = 1,KMAX_BND
          A_bnd(k)=PT * (1-sigma_bnd(k)) 
          B_bnd(k)=sigma_bnd(k) 
       enddo
       do k = 1,KMAX_MID
          A_mid(k)=(A_bnd(k+1)+A_bnd(k))/2.0      
          B_mid(k)=(B_bnd(k+1)+B_bnd(k))/2.0
       enddo
    endif
    do k = 1,KMAX_MID
       dA(k)=A_bnd(k+1)-A_bnd(k)
       dB(k)=B_bnd(k+1)-B_bnd(k)
       Eta_bnd(k)=A_bnd(k)/Pref+B_bnd(k)
       Eta_mid(k)=A_mid(k)/Pref+B_mid(k)
       dEta_i(k)=1.0/(dA(k)/Pref+dB(k))
    enddo
    Eta_bnd(KMAX_MID+1)=A_bnd(KMAX_MID+1)/Pref+B_bnd(KMAX_MID+1)
    if(me==0)write(*,*)'External_Levels ',External_Levels_Def
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
          endif
          gl_stagg(i,j)=0.25*(x1+x2+x3+x4)

          gb_stagg(i,j)=0.25*(lat_ext(i,j)+&
               lat_ext(i+1,j)+&
               lat_ext(i,j+1)+&
               lat_ext(i+1,j+1))
       enddo
    enddo

    !ensure that lon values are within [-180,+180]]
    do j=0,LJMAX
       do i=0,LIMAX
          if(gl_stagg(i,j)>180.0)gl_stagg(i,j)=gl_stagg(i,j)-360.0
          if(gl_stagg(i,j)<-180.0)gl_stagg(i,j)=gl_stagg(i,j)+360.0
       enddo
    enddo

    !test if the grid is cyclicgrid:
    !The last cell + 1 cell = first cell
    Cyclicgrid=1 !Cyclicgrid
    do j=1,ljmax
       if(mod(nint(10*(360+GIMAX*(glon(2,j)-glon(1,j)))),3600)/=0)then
          Cyclicgrid=0  !not cyclicgrid
       endif
    enddo

    if(MasterProc .and. DEBUG%GRIDVALUES)write(*,*)'CYCLICGRID:',Cyclicgrid

    !Look for poles
    !If the northernmost or southernmost lines are poles, they are not
    !considered as outer boundaries and will not be treated 
    !by "BoundaryConditions_ml".
    !If the projection is not lat lon (i.e. the poles are not lines, but points), the poles are 
    !not a problem and Pole=0, even if the grid actually include a pole.
    !Note that "Poles" is defined in subdomains

    North_pole=1
    do i=1,limax
       if(nint(glat(i,ljmax))<=88)then
          North_pole=0  !not north pole
       endif
    enddo

    South_pole=1
    do i=1,limax
       if(nint(glat(i,1))>=-88)then
          South_pole=0  !not south pole
       endif
    enddo

    Poles=0
    if(North_pole==1)then
       Poles(1)=1
       write(*,*)me,'Found North Pole'
    endif

    if(South_pole==1)then
       Poles(2)=1
       write(*,*)me,'Found South Pole'
    endif
    do j=1,LJMAX
       do i=1,LIMAX
          GridArea_m2(i,j) = GRIDWIDTH_M*GRIDWIDTH_M*xmd(i,j)
       enddo
    enddo
    
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
    endif

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


  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine DefDebugProc()
    !-------------------------------------------------------------------! 
    ! -------------- Find debug coords  and processor ------------------
    !-------------------------------------------------------------------! 

    integer ::  i, j

    debug_proc = .false.

    do i = li0, li1
       do j = lj0, lj1
          if( i_fdom(i) == DEBUG%IJ(1) .and. j_fdom(j) == DEBUG%IJ(2) ) then
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
    endif

  end subroutine DefDebugProc


  subroutine lb2ijm(imax,jmax,lon,lat,xr2,yr2,fi2,an2,xp2,yp2)
    !-------------------------------------------------------------------! 
    !   calculates coordinates xr2, yr2 (real values) from glat(lat),glon(long) 
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

    real, intent(in)    :: lon(imax,jmax),lat(imax,jmax) 
    real, intent(out)   :: xr2(imax,jmax),yr2(imax,jmax)
    real, intent(in), optional    :: fi2,an2,xp2,yp2
    integer, intent(in) :: imax,jmax
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
          enddo
       enddo
    else if(projection=='lon lat')then! lon-lat grid
       do j1 = 1, jmax
          do i1 = 1, imax
             xr2(i1,j1)=(lon(i1,j1)-glon(1,1))/(glon(2,1)-glon(1,1))+i_fdom(1)
             if(xr2(i1,j1)<0.5)xr2=xr2+360.0/(glon(2,1)-glon(1,1))
             yr2(i1,j1)=(lat(i1,j1)-glat(1,1))/(glat(1,2)-glat(1,1))+j_fdom(1)
          enddo
       enddo
    else!general projection, Use only info from glon_fdom and glat_fdom
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
!                   if(dist>great_circle_distance(glon(i1,j1),glat(i1,j1),glon_fdom(i,j) &
!                        ,glat_fdom(i,j)))then
!                      dist=great_circle_distance(glon(i1,j1),glat(i1,j1),glon_fdom(i,j) &
!                           ,glat_fdom(i,j))
!                      xr2(i1,j1)=i
!                      yr2(i1,j1)=j
!                   endif
                enddo
             enddo

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

          enddo
       enddo

    endif
  end subroutine lb2ijm

subroutine lb2ij_real(gl2,gb2,xr2,yr2,fi2,an2,xp2,yp2)
!Note: this routine is not yet CPU optimized
!-------------------------------------------------------------------! 
!      calculates coordinates xr2, yr2 (real values) from gl(lat),gb(long) 
!
!      input:  xp2,yp2:   coord. of the polar point in grid2
!              an2:   number of grid-distances from pole to equator in grid2.
!              fi2:      rotational angle for the grid2 (at i2=0).
!              i1max,j1max: number of points (grid1) in  x- og y- direction
!
!
!      output: xr2(i1,j1): i coordinates in grid2 
!              yr2(i1,j1): j coordinates in grid2 
!-------------------------------------------------------------------! 
  real, intent(in)    :: gl2,gb2 
  real, intent(out)   :: xr2,yr2
  real, intent(in), optional    :: fi2,an2,xp2,yp2

  real  :: fi_loc,an_loc,xp_loc,yp_loc
  real, parameter :: PI=3.14159265358979323,dr=PI/180.0,dri= 180.0/PI
  real    :: PId4,dr2,dist,dist2,dist3
  integer ::i,j,ip1,jp1, ir2, jr2
  real ::xscen ,yscen,zsycen,zcycen ,zxmxc,zsxmxc,zcxmxc,zsysph,zsyrot,yrot,zsxrot,zcysph,zcyrot,zcxrot,xrot

  select case (projection)
  case('Stereographic')
    PId4  =PI/4.      
    dr2   =dr*0.5   ! degrees to radians /2
    fi_loc=fi
    an_loc=an
    xp_loc=xp
    yp_loc=yp

    if(present(fi2))fi_loc=fi2
    if(present(an2))an_loc=an2
    if(present(xp2))xp_loc=xp2
    if(present(yp2))yp_loc=yp2

    xr2=xp_loc+an_loc*tan(PId4-gb2*dr2)*sin(dr*(gl2-fi_loc))
    yr2=yp_loc-an_loc*tan(PId4-gb2*dr2)*cos(dr*(gl2-fi_loc))

  case('lon lat')           ! lon-lat grid    
    if((gl2-glon(1,1))+i_fdom(1)*(glon(2,1)-glon(1,1))<360.0)then
       xr2=(gl2-glon(1,1))/(glon(2,1)-glon(1,1))+i_fdom(1)
    else          
       xr2=(gl2-360.0-glon(1,1))/(glon(2,1)-glon(1,1))+i_fdom(1)
    endif
    if(xr2<0.5)xr2=xr2+360.0/(glon(2,1)-glon(1,1))
    yr2=(gb2-glat(1,1))/(glat(1,2)-glat(1,1))+j_fdom(1)

  case('Rotated_Spherical') ! rotated lon-lat grid
    !       dx_roti=20.0
    !       grid_north_pole_longitude = -170.0
    !       grid_north_pole_latitude = 40.0
    xscen = grid_north_pole_longitude-180.0
    if(xscen<-180.0)xscen = xscen+360.0
    yscen = 90.0-grid_north_pole_latitude
    !       xscen=grid_north_pole_longitude-180.0
    !       yscen=90.0-grid_north_pole_latitude
    zsycen = sin(dr*yscen)
    zcycen = cos(dr*yscen)
    !
    zxmxc  = dr*(gl2 - xscen)
    zsxmxc = sin(zxmxc)
    zcxmxc = cos(zxmxc)
    zsysph = sin(dr*gb2)
    zcysph = cos(dr*gb2)
    zsyrot = zcycen*zsysph - zsycen*zcysph*zcxmxc
    zsyrot = amax1(zsyrot,-1.0)
    zsyrot = amin1(zsyrot,+1.0)
    yrot = asin(zsyrot)
    zcyrot = cos(yrot)
    zcxrot = (zcycen*zcysph*zcxmxc +&
         zsycen*zsysph)/zcyrot
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

  case default ! general projection, Use only info from glon_fdom and glat_fdom
    !first find closest by testing all gridcells. 
    call StopAll('lb2ij: conversion broken 27 Oct 2015, Peter')
    !glon_fdom is no more defined. Could easily rewrite if necessary
    dist2=0.0
    dist3=0.0
    dist=10.0!max distance is PI
    do j=1,JJFULLDOM
      do i=1,IIFULLDOM
!       if(dist>great_circle_distance(gl2,gb2,glon_fdom(i,j) &
!           ,glat_fdom(i,j)))then
!         dist=great_circle_distance(gl2,gb2,glon_fdom(i,j) &
!              ,glat_fdom(i,j))
!         xr2=i
!         yr2=j
!       endif
      enddo
    enddo

    !find the real part of i and j by comparing distances to neighbouring cells
    !
    !     C
    !    /|\
    !   / | \
    !  /  |  \
    ! A---D---B
    !
    !A=(i,j) ,B=(i+1,j), C=(gl2,gb2)
    !dist=AC, dist2=BC, dist3=AB
    !AD=(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3)
    !
    ir2 = nint(xr2)
    jr2 = nint(yr2)
    ip1=ir2+1
    if(ip1>IIFULLDOM)ip1=ip1-2
!      dist2=great_circle_distance(gl2,gb2,glon_fdom(ip1,jr2),glat_fdom(ip1,jr2))
!      dist3=great_circle_distance( glon_fdom(ir2,jr2), &
!           glat_fdom(ir2,jr2), &
!           glon_fdom(ip1,jr2), &
!           glat_fdom(ip1,jr2))

    xr2=xr2+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)


    jp1=jr2+1
    if(jp1>JJFULLDOM)jp1=jp1-2

!      dist2=great_circle_distance(gl2,gb2,glon_fdom(ir2,jp1),glat_fdom(ir2,jp1))
    !GFORTRAN CHANGE
!      dist3=great_circle_distance( glon_fdom(ir2,jr2), &
!           glat_fdom(ir2,jr2), &
!           glon_fdom(ir2,jp1), & 
!           glat_fdom(ir2,jp1) )

    yr2=yr2+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)

  endselect
endsubroutine lb2ij_real
subroutine lb2ij_int(gl2,gb2,ix,iy)
  real, intent(in)    :: gl2,gb2
  integer, intent(out):: ix,iy
  real ::x,y
! stations can easily be defined exactly at gridcell boundaries
! 1.0E-7 is to ensure same rounding for all CPUs
  call lb2ij_real(gl2,gb2,x,y)
  ix=nint(x+1.0E-7)
  iy=nint(y+1.0E-7)
endsubroutine lb2ij_int 

  subroutine ij2lbm(imax,jmax,glon,glat,fi,an,xp,yp)
    !-------------------------------------------------------------------! 
    !      calculates l(lat),b(long) (geographical coord.) 
    !      in every grid point. 
    !
    !      input:  xp,yp:   coord. of the polar point.
    !              an:      number of grid-distances from pole to equator.
    !              fi:      rotational angle for the x,y grid (at i=0).
    !              imax,jmax:   number of points in  x- og y- direction
    !              glmin:   gives min.value of geographical lenght
    !                       =>  glmin <= l <= glmin+360.  
    !                           (example glmin = -180. or 0.)
    !                       if "geopos","georek" is used
    !                       then glmin must be the lenght i(1,1) in the
    !                       geographical grid (gl1 to "geopos")
    !      output: gl(ii,jj): longitude glmin <= l <= glmin+360. 
    !              gb(ii,jj): latitude  -90. <= b <= +90. 
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
    !      calculates l(lat),b(long) (geographical coord.) 
    !      from i,j coordinates in polar stereographic projection 
    !
    !      input:  i,j
    !              xp,yp:   coord. of the polar point.
    !              an:      number of grid-distances from pole to equator.
    !              fi:      rotational angle for the x,y grid (at i=0).
    !              imax,jmax:   number of points in  x- og y- direction
    !              glmin:   gives min.value of geographical lenght
    !                       =>  glmin <= l <= glmin+360.  
    !                           (example glmin = -180. or 0.)
    !                       if "geopos","georek" is used
    !                       then glmin must be the lenght i(1,1) in the
    !                       geographical grid (gl1 to "geopos")
    !      output: lon: longitude glmin <= lon <= glmin+360. 
    !              lat: latitude  -90. <= lat <= +90. 
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

    !    do j = 1, jmax          
    dy  = yp - j            
    dy2 = dy*dy
    !       do i = 1, imax       

    dx = i - xp    ! ds - changed
    rp = sqrt(dx*dx+dy2)           ! => distance to pole
    rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
    rl = 0.0
    if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
    if (rl <  glmin)   rl = rl + 360.0
    if (rl >  glmax)   rl = rl - 360.0
    lon=rl                   ! longitude
    lat=rb                   ! latitude
    !       end do ! i
    !    end do ! j

  end subroutine ij2lb

  subroutine ij2ijm(in_field,imaxin,jmaxin,out_field,imaxout,jmaxout, &
       fiin,anin,xpin,ypin,fiout,anout,xpout,ypout)

    !   Converts data (in_field) stored in coordinates (fiin,anin,xpin,ypin) 
    !   into data (out_field) in coordinates (fiout,anout,xpout,ypout)
    !   pw august 2002


    integer, intent(in) :: imaxin,jmaxin,imaxout,jmaxout
    real, intent(in) :: fiin,anin,xpin,ypin,fiout,anout,xpout,ypout
    real, intent(in) :: in_field(imaxin,jmaxin)! Field to be transformed
    real, intent(out) :: out_field(imaxout,jmaxout)! Field to be transformed

    real, allocatable,dimension(:,:) :: x,y,glat,glon
    integer alloc_err,i,j,i2,j2
    logical :: interpolate
    real :: f11,f12,f21,f22

    interpolate = .true.
    !        interpolate = .false.

    allocate(x(imaxout,jmaxout), stat=alloc_err)
    allocate(y(imaxout,jmaxout), stat=alloc_err)
    allocate(glat(imaxout,jmaxout), stat=alloc_err)
    allocate(glon(imaxout,jmaxout), stat=alloc_err)
    if ( alloc_err /= 0 ) WRITE(*,*) 'MPI_ABORT: ', "ij2ij alloc failed" 
    if ( alloc_err /= 0 ) call  MPI_ABORT(MPI_COMM_CALC,9,IERROR) 

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
       call  MPI_ABORT(MPI_COMM_CALC,9,IERROR) 
    endif


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

          enddo
       enddo
    else

       do j = 1, jmaxout
          do i = 1,imaxout
             out_field(i,j) =in_field(nint(x(i,j)),nint(y(i,j)))
          enddo
       enddo

    endif

    deallocate(x,stat=alloc_err)
    deallocate(y,stat=alloc_err)
    deallocate(glat,stat=alloc_err)
    deallocate(glon,stat=alloc_err)
    if ( alloc_err /= 0 )   WRITE(*,*) 'MPI_ABORT: ', "ij2ijde-alloc_err" 
    if ( alloc_err /= 0 ) call  MPI_ABORT(MPI_COMM_CALC,9,IERROR) 

  end subroutine ij2ijm

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
       endif
    endif
  endsubroutine range_check
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
       !!  lat=mod   (lat      , 90.0)       ! lat/gb_stagg range  -90 .. 90
       lon=modulo(lon+180.0,360.0)-180.0 ! lon/gl_stagg range -180 .. 180
       call range_check(trim(msg)//" lon",lon,(/-180.0,180.0/),fatal=.true.)
    endif
  endsubroutine coord_check
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
    endselect
  endfunction coord_in_domain
  function coord_in_processor(lon,lat,iloc,jloc,iglob,jglob) result(in)
    !-------------------------------------------------------------------!
    ! Is coord (lon/lat) is inside local domain?
    !-------------------------------------------------------------------!
    real, intent(inout) :: lon,lat
    integer, intent(out),optional:: iloc,jloc,iglob,jglob
    logical :: in
    in=coord_in_domain("processor",lon,lat,iloc,jloc,iglob,jglob)
  endfunction coord_in_processor
  function coord_in_gridbox(lon,lat,iloc,jloc,iglob,jglob) result(in)
    !-------------------------------------------------------------------!
    ! Is coord (lon/lat) is inside gridbox(iloc,jloc)?
    !-------------------------------------------------------------------!
    real, intent(inout) :: lon,lat
    integer, intent(inout) :: iloc,jloc
    integer, intent(out),optional:: iglob,jglob
    logical :: in
    in=coord_in_domain("gridbox",lon,lat,iloc,jloc,iglob,jglob)
  endfunction coord_in_gridbox

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
    !make interpolation coefficients to convert the levels defined in meteo 
    !into the levels defined in Vertical_levels.txt
    integer ::k,k_met
    real ::p_met,p_mod,p1,p2
    if(.not. allocated(k1_met))allocate(k1_met(KMAX_MID),k2_met(KMAX_MID),x_k1_met(KMAX_MID))
    if(.not. allocated(A_bnd_met))then
       allocate(A_bnd_met(KMAX_MID+1),B_bnd_met(KMAX_MID+1))
       A_bnd_met=A_bnd
       B_bnd_met=B_bnd
    endif

    if(me_mpi==0)then

       !only me=0 has the values for A_bnd_met and B_bnd_met
       do k=1,KMAX_MID
          P_mod=A_mid(k)+Pref*B_mid(k)
          !find the lowest met level higher than the model level 
          !do k_met=1,KMAX_MET
          k_met=KMAX_MET-1
          p_met=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
          do while(p_met>P_mod.and.k_met>1)
             ! write(*,*)P_mod,p_met
             k_met=k_met-1
             p_met=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
          enddo
          k1_met(k)=k_met
          k2_met(k)=k_met+1
          k_met=k1_met(k)
          p1=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
          k_met=k2_met(k)
          p2=0.5*(A_bnd_met(k_met+1)+A_bnd_met(k_met))+Pref*0.5*(B_bnd_met(k_met+1)+B_bnd_met(k_met))
          x_k1_met(k)=(p_mod-p2)/(p1-p2)          
          write(*,77)k, ' interpolated from levels ', k1_met(k),' and ',k2_met(k),P_mod,p1,p2,x_k1_met(k)
77        format(I4,A,I3,A,I3,13f11.3)
          if(x_k1_met(k)<-0.00001 .or. (1.0-x_k1_met(k))<-0.00001)then
             write(*,*)'WARNING: Extrapolation of data. This is NOT recommended for several metfields'
          endif

       enddo
    endif

    CALL MPI_BCAST(k1_met,4*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_BCAST(k2_met,4*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_BCAST(x_k1_met,8*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,IERROR)

  end subroutine make_vertical_levels_interpolation_coeff

  subroutine lb_rot2lb(xsph,ysph,xrot,yrot,grid_north_pole_longitude,grid_north_pole_latitude)
    !
    !  compute spherical coordinates as function of
    !  spherical rotated coordinates
    !
    !  conversion between spherical (xsph,ysph) and spherical rotated
    !  (xrot,yrot) coordinates. (xcen,ycen) is the position of the
    !  rotated equator/greenwich in terms of (longitude,latitude).
    !  all input and output values are given in degrees.
    !
    ! grid_north_pole_longitude: geographical (non-rotated) coordinates of the "north pole" from the rotated grid (No polar bears there).
    ! (typically out of the grid, since it is singular).
    !
    ! xcen: geographical (non-rotated) coordinates of the (lon=0 lat=0) point where lonlat are in the rotated grid
    ! (typically in the middle of the grid, since it is "flat")
    !
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


  subroutine nf90_get_var_extended(ncFileID,varID,Var,i0,i1,j0,j1,ishift_in,jshift_in)

    !fetch a 2D array over an extended subdomain.
    !i.e. extended arrays are overlapping, and parts outside the fulldomain are extrapolated linearly.
    implicit none
    integer, intent(in) ::ncFileID,varID
    real, intent(inout) ::Var(i0:i1,j0:j1)!the extended local array 
    integer, intent(in)::i0,i1,j0,j1
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
       enddo
    enddo
    do i=iloc_end+1,i1
       do j=jloc_start,jloc_end
          Var(i,j)=2.0*Var(i-1,j)-Var(i-2,j)
       enddo
    enddo
    do i=i0,i1!now they are all defined
       do j=jloc_start-1,j0,-1
          Var(i,j)=2.0*Var(i,j+1)-Var(i,j+2)
       enddo
    enddo
    do i=i0,i1!now they are all defined
       do j=jloc_end+1,j1
          Var(i,j)=2.0*Var(i,j-1)-Var(i,j-2)
       enddo
    enddo

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
  endif
endsubroutine RestrictDomain
endmodule GridValues_ml
!==============================================================================
