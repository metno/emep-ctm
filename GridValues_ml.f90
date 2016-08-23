! <GridValues_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!***************************************************************************! 
!* 
!*  Copyright (C) 2007-2012 met.no
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
!***************************************************************************! 

Module GridValues_ml

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !
  !  Define parameters and variables associated with 3-D grid and its
  !  geography.
  !
  ! History: 
  ! March - changed folllwing Steffen's optimisation/correction of sigma_mid.
  ! January 2001 : Created by ds from old defconstants, made enough changes
  ! to get this into F90, and to make x,y inputs to the position subroutine,
  ! but the basic equations are untouched.
  ! October 2001 hf added call to ReadField (which now does global2local)
  ! Nov. 2001 - tidied up a bit (ds). Use statements moved to top of module
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  use CheckStop_ml,    only : CheckStop,StopAll
  use Functions_ml,    only : great_circle_distance
  use Io_Nums_ml,       only: IO_LOG

  use MetFields_ml 
  use ModelConstants_ml,    only : KMAX_BND, KMAX_MID  &! vertical extent
       ,DEBUG_i, DEBUG_j  &    ! full-domain coordinate of debug-site
       ,DEBUG_GRIDVALUES  &
       ,MasterProc        &
       ,NPROC, IIFULLDOM,JJFULLDOM, RUNDOMAIN&
       ,PT, Pref          &
       ,NMET, METSTEP
  use Par_ml, only : &
       MAXLIMAX,MAXLJMAX & ! max. possible i, j in this domain
       ,limax,ljmax        & ! actual max.   i, j in this domain
       ,li0,li1,lj0,lj1    & ! for debugging TAB
        ,GIMAX,GJMAX       & ! Size of rundomain
       ,IRUNBEG,JRUNBEG    & ! start of user-specified domain
       ,gi0,gj0     & ! full-dom coordinates of domain lower l.h. corner
       ,gi1,gj1     & ! full-dom coordinates of domain uppet r.h. corner
       ,me          & ! local processor
      , parinit
  use PhysicalConstants_ml, only : GRAV, PI ! gravity, pi
  use TimeDate_ml,          only : current_date, date,Init_nmdays,nmdays,startdate
  use TimeDate_ExtraUtil_ml,only : nctime2idate,date2string

  implicit none
  private

  !-- contains subroutine:

  Public :: DefGrid ! => GRIDWIDTH_M, map-factor stuff, calls other routines
  Public :: DefDebugProc ! =>  sets debug_proc, debug_li, debug_lj
  Public :: ij2lbm  ! polar stereo grid to longitude latitude
  Public :: lb2ijm  ! longitude latitude to grid in polar stereo
  Public :: ij2ijm  ! polar grid1 to polar grid2
  Public :: lb2ij   ! longitude latitude to (i,j) in any grid projection
  Public :: ij2lb  ! polar stereo grid to longitude latitude

  Public :: GlobalPosition
  private :: Position ! => lat(glat), long (glon)
  Public :: coord_in_gridbox,  &  ! Are coord (lon/lat) is inside gridbox(i,j)?
       coord_in_processor    ! Are coord (lon/lat) is inside local domain?

  public :: GridRead,Getgridparams
  private :: Alloc_GridFields
  private :: GetFullDomainSize

  !** 1) Public (saved) Variables from module:
  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
  real MPIbuff

  real, public, save :: &
       xp, yp  &   ! Coordinates of North pole (from infield)
       , fi      &   ! projections rotation angle around y axis (from infield)
       , AN      &   ! Distance on the map from pole to equator (No. of cells)
       ,GRIDWIDTH_M &! width of grid at 60N, in meters (old "h")(from infield)
       ,ref_latitude ! latitude at which projection is true (degrees)

  !Rotated_Spherical grid prarameters
  real, public, save :: grid_north_pole_latitude,grid_north_pole_longitude
  real, public, save :: dx_rot,dx_roti,x1_rot,y1_rot

  !/ Variables to define full-domain (fdom) coordinates of local i,j values. 

  integer, public, allocatable, save, dimension(:) :: i_fdom !fdom coordinates
  integer, public, allocatable, save, dimension(:) :: j_fdom !of local i,j

  ! and reverse:
  integer, public, allocatable, save,dimension(:) :: i_local  !local coordinates
  integer, public, allocatable, save,dimension(:) :: j_local  !of full-domain i,j

  !Parameters for Vertical Hybrid coordinates:
  real, public, save,  dimension(KMAX_BND) ::  &
       A_bnd !Unit Pa.  first constant, defined at layer boundary
  ! (i.e. half levels in EC nomenclature)
  real, public, save,  dimension(KMAX_BND) ::  &
       B_bnd !Unit 1.  second constant, defined at layer boundary
  ! (i.e. half levels in EC nomenclature)
  real, public, save,  dimension(KMAX_MID) ::  &
       A_mid !Unit Pa.  first constant, defined at middle of layer
  ! (i.e. full levels in EC nomenclature)
  real, public, save,  dimension(KMAX_MID) ::  &
       B_mid !Unit 1.  second constant, defined at middle of layer
  ! (i.e. full levels in EC nomenclature)
  real, public, save,  dimension(KMAX_MID) ::  &
       dA !Unit Pa.  A_bnd(k+1)-A_bnd(k) 
  real, public, save,  dimension(KMAX_MID) ::  &
       dB !Unit 1.  B_bnd(k+1)-B_bnd(k) 
  ! P = A + B*PS
  ! eta = A/Pref + B

  real, public, save,  dimension(KMAX_BND) ::  &
       sigma_bnd ! sigma, layer boundary 

  real, public, save,  dimension(KMAX_MID) ::  &
       sigma_mid   ! sigma layer midpoint

  real, public, save,  dimension(KMAX_MID) ::  carea    ! for budgets?

  real, public, save,allocatable,  dimension(:,:) :: &
       glon       &         !longitude of gridcell centers
       ,glat                 !latitude  of gridcell centers
  real, public, save,allocatable,  dimension(:,:) :: &
       gl_stagg   &         !longitude of gridcell corners 
       ,gb_stagg             !latitude  of gridcell corners
  !NB: gl_stagg, gb_stagg are here defined as the average of the four
  !    surrounding gl gb.
  !    These differ slightly from the staggered points in the (i,j) grid. 

  real, public, allocatable,  dimension(:,:) :: &
       glat_fdom, &         !latitude of gridcell centers
       glon_fdom            !longitude of gridcell centers


  real, public, save :: gbacmax,gbacmin,glacmax,glacmin

  ! EMEP grid definitions (old and official)
  real, public, parameter :: xp_EMEP_official=8.&
       ,yp_EMEP_official=110.0&
       ,fi_EMEP=-32.0&
       ,ref_latitude_EMEP=60.0&
       ,GRIDWIDTH_M_EMEP=50000.0&
       ,an_EMEP = 237.7316364 &! = 6.370e6*(1.0+0.5*sqrt(3.0))/50000.
       ,xp_EMEP_old =  43.0&
       ,yp_EMEP_old = 121.0


  !/** Map factor stuff:

  real, public, save,allocatable, dimension(:,:) ::  &
       xm_i     & ! map-factor in i direction, between cell j and j+1
       ,xm_j     & ! map-factor in j direction, between cell i and i+1
       ,xm2      & ! xm*xm: area factor in the middle of a cell (i,j)
       ,xmd        ! 1/xm2  

  real, public, save,allocatable, dimension(:,:) ::  &
       xm2ji  &
       ,xmdji

  !/** Grid Area

  real, public, save,allocatable, dimension(:,:) :: GridArea_m2

  integer, public, save :: &
       debug_li=-99, debug_lj=-99         ! Local Coordinates of debug-site
  logical, public, save :: debug_proc  ! Processor with debug-site

  character (len=100),public::projection
  integer, public, parameter :: MIN_ADVGRIDS = 5 !minimum size of a subdomain
  integer, public :: Poles(2) !Poles(1)=1 if North pole is found, Poles(2)=1:SP
  integer, public :: Pole_Singular !Pole_included=1 or 2 if the grid include at least one pole and has lat lon projection


contains


  subroutine GridRead(cyclicgrid)

    !   the subroutine reads the grid parameters (projection, resolution etc.)
    !   defined by the meteorological fields
    !

    implicit none

    integer,  intent(out)      :: cyclicgrid
    integer                    :: nyear,nmonth,nday,nhour,k
    integer                    :: KMAX,MIN_GRIDS
    character (len = 100),save :: filename !name of the input file
    character (len=230) :: txt


    nyear=startdate(1)
    nmonth=startdate(2)
    nday=startdate(3)
    nhour=0
    current_date = date(nyear, nmonth, nday, nhour, 0 )
    call Init_nmdays( current_date )

    !*********initialize grid parameters*********
56  FORMAT(a5,i4.4,i2.2,i2.2,a3)
    write(filename,56)'meteo',nyear,nmonth,nday,'.nc'
 
    if(MasterProc)write(*,*)'reading domain sizes from ',trim(filename)

    call GetFullDomainSize(filename,IIFULLDOM,JJFULLDOM,KMAX,Pole_Singular,projection)

    call CheckStop(KMAX_MID/=KMAX,"vertical cordinates not yet flexible")

    allocate(i_local(IIFULLDOM))
    allocate(j_local(JJFULLDOM))
    allocate(glat_fdom(IIFULLDOM,JJFULLDOM))!should be removed from code
    allocate(glon_fdom(IIFULLDOM,JJFULLDOM))!should be removed from code

     !set RUNDOMAIN default values where not defined
    if(RUNDOMAIN(1)<1)RUNDOMAIN(1)=1
    if(RUNDOMAIN(2)<1)RUNDOMAIN(2)=IIFULLDOM
    if(RUNDOMAIN(3)<1)RUNDOMAIN(3)=1
    if(RUNDOMAIN(4)<1)RUNDOMAIN(4)=JJFULLDOM
    if(MasterProc)then
    55 format(A,I5,A,I5)
       write(*,55) 'FULLDOMAIN has sizes ',IIFULLDOM,' X ',JJFULLDOM
       write(IO_LOG,55)'FULLDOMAIN has sizes ',IIFULLDOM,' X ',JJFULLDOM
       write(*,55)'RUNDOMAIN  x coordinates from ',RUNDOMAIN(1),' to ',RUNDOMAIN(2)
       write(IO_LOG,55)'RUNDOMAIN  x coordinates from ',RUNDOMAIN(1),' to ',RUNDOMAIN(2)
       write(*,55)'RUNDOMAIN  y coordinates from ',RUNDOMAIN(3),' to ',RUNDOMAIN(4)
       write(IO_LOG,55)'RUNDOMAIN  y coordinates from ',RUNDOMAIN(3),' to ',RUNDOMAIN(4)
    endif


    MIN_GRIDS=5
    call parinit(MIN_GRIDS,Pole_Singular)     !subdomains sizes and position

    call Alloc_MetFields(MAXLIMAX,MAXLJMAX,KMAX_MID,KMAX_BND,NMET)

    call Alloc_GridFields(GIMAX,GJMAX,MAXLIMAX,MAXLJMAX,KMAX_MID,KMAX_BND)

    call Getgridparams(filename,GRIDWIDTH_M,xp,yp,fi,&
         ref_latitude,sigma_mid,Nhh,nyear,nmonth,nday,nhour,nhour_first&
         ,cyclicgrid)


    if(MasterProc .and. DEBUG_GRIDVALUES)then
       write(*,*)'sigma_mid:',(sigma_mid(k),k=1,20)
       write(*,*)'grid resolution:',GRIDWIDTH_M
       write(*,*)'xcoordinate of North Pole, xp:',xp
       write(*,*)'ycoordinate of North Pole, yp:',yp
       write(*,*)'longitude rotation of grid, fi:',fi
       write(*,*)'true distances latitude, ref_latitude:',ref_latitude
    endif

    call DefGrid()!defines: i_fdom,j_fdom,i_local, j_local,xmd,xm2ji,xmdji,
    !         sigma_bnd,carea,gbacmax,gbacmin,glacmax,glacmin

  end subroutine GridRead



  subroutine GetFullDomainSize(filename,IIFULLDOM,JJFULLDOM,KMAX,Pole_Singular,projection)

    !
    ! Get input grid sizes 
    !

  use netcdf

    implicit none

    character (len = *), intent(in) ::filename
    integer, intent(out):: IIFULLDOM,JJFULLDOM,KMAX,Pole_Singular
    character (len = *), intent(out) ::projection

    integer :: status,ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID
    integer :: GIMAX_file,GJMAX_file,KMAX_file
    real,allocatable :: latitudes(:)


    if(MasterProc)then
       print *,'Defining grid properties from ',trim(filename)
       !open an existing netcdf dataset
       status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
       if(status /= nf90_noerr) then
          print *,'not found',trim(filename)
          call StopAll("File not found")
       endif

       !          print *,'  reading ',trim(filename)
       projection=''
       call check(nf90_get_att(ncFileID,nf90_global,"projection",projection))
       if(trim(projection)=='Rotated_Spherical'.or.trim(projection)=='rotated_spherical'&
            .or.trim(projection)=='rotated_pole'.or.trim(projection)=='rotated_latitude_longitude')then
          projection='Rotated_Spherical'
       endif
       write(*,*)'projection: ',trim(projection)

       !get dimensions id
       if(trim(projection)=='Stereographic') then
          call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
       elseif(trim(projection)==trim('lon lat')) then
          call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
       else
          !     write(*,*)'GENERAL PROJECTION ',trim(projection)
          call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
       endif

       call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
       call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))
       call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeVarID))

       !get dimensions length
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_file))
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_file))
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_file))
!       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=Nhh))

       write(*,*)'dimensions input grid:',GIMAX_file,GJMAX_file,KMAX_file!,Nhh

       IIFULLDOM=GIMAX_file
       JJFULLDOM=GJMAX_file
       KMAX     =KMAX_file

       Pole_Singular=0
      if(trim(projection)==trim('lon lat')) then
         !find wether poles are included (or almost included) in grid
         !
         !If some cells are to narrow (Poles in lat lon coordinates),
         !this will give too small time steps in the Advection,
         !because of the constraint that the Courant number should be <1.
         !
         !If Poles are found and lon-lat coordinates are used the Advection scheme
         !will be modified to be able to cope with the singularity
         !the advection routine will not work efficiently with NPROCY>2 in this case

          allocate(latitudes(JJFULLDOM))
          call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
          call check(nf90_get_var(ncFileID, varID,latitudes  ))
          if(latitudes(JJFULLDOM)>88.0)then
             write(*,*)'The grid is singular at North Pole'
             Pole_Singular=Pole_Singular+1
          endif
          if(latitudes(1)<-88.0)then
             write(*,*)'The grid is singular at South Pole'
             Pole_Singular=Pole_Singular+1
          endif
          deallocate(latitudes)
      endif
    endif
    CALL MPI_BCAST(IIFULLDOM ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(JJFULLDOM ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(KMAX ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(Pole_Singular ,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(projection ,len(projection),MPI_BYTE,0,MPI_COMM_WORLD,INFO)

  end subroutine GetFullDomainSize


  subroutine Getgridparams(filename,GRIDWIDTH_M,xp,yp,fi,&
       ref_latitude,sigma_mid,Nhh,nyear,nmonth,nday,nhour,nhour_first&
       ,cyclicgrid)
    !
    ! Get grid and time parameters as defined in the meteo file
    ! Do some checks on sizes and dates
    !
    ! This routine is called only once (and is therefore not optimized for speed)
    !
!    use ChemSpecs_adv_ml,  only: NSPEC_ADV         ! => No. species 
!    use ChemSpecs_shl_ml,  only: NSPEC_SHL         ! => No. species 
!    use Chemfields_ml, only:xn_adv, xn_shl, xn_bgn, PM25_water, PM25_water, PM25_water_rh50&
!         ,AOD,cfac,so2nh3_24hr,Grid_snow&
!         ,NSPEC_BGN,NSPEC_COL

  use netcdf

    implicit none

    character (len = *), intent(in) ::filename
    integer, intent(in):: nyear,nmonth,nday,nhour
    real, intent(out) :: GRIDWIDTH_M,xp,yp,fi, ref_latitude,sigma_mid(KMAX_MID)
    integer, intent(out):: Nhh,nhour_first,cyclicgrid

    integer :: nseconds(1),n1,i,j,im,jm,i0,j0
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID
    integer :: GIMAX_file,GJMAX_file,KMAX_file,ihh,ndate(4)
    !    realdimension(-1:GIMAX+2,-1:GJMAX+2) ::xm_global,xm_global_j,xm_global_i
    real,allocatable,dimension(:,:) ::xm_global,xm_global_j,xm_global_i
    integer :: status,iglobal,jglobal,info,South_pole,North_pole,Ibuff(2)
    real :: ndays(1),x1,x2,x3,x4
    character (len = 50) :: timeunit
 
    allocate(xm_global(-1:GIMAX+2,-1:GJMAX+2))
    allocate(xm_global_j(-1:GIMAX+2,-1:GJMAX+2))
    allocate(xm_global_i(-1:GIMAX+2,-1:GJMAX+2))

    if(MasterProc)then
       call CheckStop(GIMAX+IRUNBEG-1 > IIFULLDOM, "GridRead: I outside domain" )
       call CheckStop(GJMAX+JRUNBEG-1 > JJFULLDOM, "GridRead: J outside domain" )

       call CheckStop(nhour/=0 .and. nhour /=3,&
            "ReadGrid: must start at nhour=0 or 3")

      print *,'Defining grid properties from ',trim(filename)
       !open an existing netcdf dataset
       status = nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID)
       if(status /= nf90_noerr) then
          print *,'not found',trim(filename)
          call StopAll("File not found")
       endif

       !          print *,'  reading ',trim(filename)

       !get dimensions id
       if(trim(projection)=='Stereographic') then
          call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
       elseif(trim(projection)==trim('lon lat')) then
          call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
       else
          !     write(*,*)'GENERAL PROJECTION ',trim(projection)
          call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
       endif

       call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
       call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))
       call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeVarID))
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=Nhh))

       call CheckStop(24/Nhh, METSTEP,          "GridRead: METSTEP != meteostep" )

       call check(nf90_get_att(ncFileID,timeVarID,"units",timeunit))

       ihh=1
       n1=1
       if(trim(timeunit(1:19))==trim("days since 1900-1-1"))then
          write(*,*)'Date in days since 1900-1-1 0:0:0'
          call check(nf90_get_var(ncFileID,timeVarID,ndays,&
               start=(/ihh/),count=(/n1 /)))
          call nctime2idate(ndate,ndays(1))  ! for printout: msg="meteo hour YYYY-MM-DD hh"
       else
          call check(nf90_get_var(ncFileID,timeVarID,nseconds,&
               start=(/ihh/),count=(/n1 /)))
          call nctime2idate(ndate,nseconds(1)) ! default
       endif
       nhour_first=ndate(4)
 
       call CheckStop(ndate(1), nyear,  "NetCDF_ml: wrong year" )
       call CheckStop(ndate(2), nmonth, "NetCDF_ml: wrong month" )
       call CheckStop(ndate(3), nday,   "NetCDF_ml: wrong day" )

       do ihh=1,Nhh

          if(trim(timeunit(1:19))==trim("days since 1900-1-1"))then
             call check(nf90_get_var(ncFileID, timeVarID, ndays,&
                  start=(/ ihh /),count=(/ n1 /)))
             call nctime2idate(ndate,ndays(1))
             write(*,*)'ndays ',ndays(1),ndate(3),ndate(4)
          else
             call check(nf90_get_var(ncFileID, timeVarID, nseconds,&
                  start=(/ ihh /),count=(/ n1 /)))
             call nctime2idate(ndate,nseconds(1))
          endif
          write(*,*)ihh,METSTEP,nhour_first, ndate(4)
          call CheckStop( mod((ihh-1)*METSTEP+nhour_first,24), ndate(4),  &
               date2string("NetCDF_ml: wrong hour YYYY-MM-DD hh",ndate))

       enddo


       !get global attributes
       call check(nf90_get_att(ncFileID,nf90_global,"Grid_resolution",GRIDWIDTH_M))
       write(*,*)"Grid_resolution",GRIDWIDTH_M
       if(trim(projection)=='Stereographic')then
          call check(nf90_get_att(ncFileID,nf90_global,"ref_latitude",ref_latitude))
          call check(nf90_get_att(ncFileID, nf90_global, "xcoordinate_NorthPole" &
               ,xp ))
          call check(nf90_get_att(ncFileID, nf90_global, "ycoordinate_NorthPole" &
               ,yp ))
          call check(nf90_get_att(ncFileID, nf90_global, "fi",fi ))

          call GlobalPosition
       elseif(trim(projection)==trim('lon lat')) then
          ref_latitude=60.
          xp=0.0
          yp=GJMAX
          fi =0.0
          call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
          call check(nf90_get_var(ncFileID, varID, glon_fdom(1:IIFULLDOM,1) ))
          do i=1,IIFULLDOM
             if(glon_fdom(i,1)>180.0)glon_fdom(i,1)=glon_fdom(i,1)-360.0
             if(glon_fdom(i,1)<-180.0)glon_fdom(i,1)=glon_fdom(i,1)+360.0
          enddo
          do j=1,JJFULLDOM
             glon_fdom(:,j)=glon_fdom(:,1)
          enddo
          call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
          call check(nf90_get_var(ncFileID, varID, glat_fdom(1,1:JJFULLDOM) ))
          do i=1,IIFULLDOM
             glat_fdom(i,:)=glat_fdom(1,:)
          enddo
       else
          ref_latitude=60.
          xp=0.0
          yp=GJMAX
          fi =0.0
          if(trim(projection)=='Rotated_Spherical')then
             call check(nf90_get_att(ncFileID,nf90_global,"grid_north_pole_latitude",grid_north_pole_latitude))
             write(*,*)"grid_north_pole_latitude",grid_north_pole_latitude
             call check(nf90_get_att(ncFileID,nf90_global,"grid_north_pole_longitude",grid_north_pole_longitude))
             write(*,*)"grid_north_pole_longitude",grid_north_pole_longitude
             call check(nf90_inq_varid(ncid = ncFileID, name = "i", varID = varID))
             call check(nf90_get_var(ncFileID, varID, glon_fdom(1:2,1)))!note that i is one dimensional
             x1_rot=glon_fdom(1,1)
             dx_rot=glon_fdom(2,1)-glon_fdom(1,1)
             call check(nf90_inq_varid(ncid = ncFileID, name = "j", varID = varID))
             call check(nf90_get_var(ncFileID, varID, glon_fdom(1,1)))!note that j is one dimensional
             y1_rot=glon_fdom(1,1)
             write(*,*)"rotated lon lat of (i,j)=(1,1)",x1_rot,y1_rot
             write(*,*)"resolution",dx_rot
          endif
          call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
          call check(nf90_get_var(ncFileID, varID, glon_fdom(1:IIFULLDOM,1:JJFULLDOM) ))

          call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
          call check(nf90_get_var(ncFileID, varID, glat_fdom(1:IIFULLDOM,1:JJFULLDOM) ))
          do j=1,JJFULLDOM
             do i=1,IIFULLDOM
                if(glon_fdom(i,j)>180.0)glon_fdom(i,j)=glon_fdom(i,j)-360.0
                if(glon_fdom(i,j)<-180.0)glon_fdom(i,j)=glon_fdom(i,j)+360.0
             enddo
          enddo

       endif
       !get variables
       status=nf90_inq_varid(ncid=ncFileID, name="map_factor", varID=varID)

       if(status == nf90_noerr)then
          !mapping factor at center of cells is defined
          !make "staggered" map factors
          call check(nf90_get_var(ncFileID, varID, xm_global(1:GIMAX,1:GJMAX) &
               ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
          do j=1,GJMAX
             do i=1,GIMAX-1
                xm_global_j(i,j)=0.5*(xm_global(i,j)+xm_global(i+1,j))
             enddo
          enddo
          i=GIMAX
          do j=1,GJMAX
             xm_global_j(i,j)=1.5*xm_global(i,j)-0.5*xm_global(i-1,j)
          enddo
          do j=1,GJMAX-1
             do i=1,GIMAX
                xm_global_i(i,j)=0.5*(xm_global(i,j)+xm_global(i,j+1))
             enddo
          enddo
          j=GJMAX
          do i=1,GIMAX
             xm_global_i(i,j)=1.5*xm_global(i,j)-0.5*xm_global(i,j-1)
          enddo

       else
          !map factor are already staggered
          status=nf90_inq_varid(ncid=ncFileID, name="map_factor_i", varID=varID)

          !BUGCHECK - moved here... (deleted if loop)
          call CheckStop( status, nf90_noerr, "erro rreading map factor" )

          write(*,*)GIMAX,GJMAX,IRUNBEG,JRUNBEG
          call check(nf90_get_var(ncFileID, varID, xm_global_i(1:GIMAX,1:GJMAX) &
               ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
          call check(nf90_inq_varid(ncid=ncFileID, name="map_factor_j", varID=varID))
          call check(nf90_get_var(ncFileID, varID, xm_global_j(1:GIMAX,1:GJMAX) &
               ,start=(/ IRUNBEG,JRUNBEG /),count=(/ GIMAX,GJMAX /)))
       endif

       call check(nf90_inq_varid(ncid = ncFileID, name = "k", varID = varID))
       call check(nf90_get_var(ncFileID, varID, sigma_mid ))

       call check(nf90_close(ncFileID))

    endif !me=0



    CALL MPI_BCAST(Nhh,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(GRIDWIDTH_M,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(ref_latitude,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xp,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(yp,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(fi,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(sigma_mid,8*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xm_global_i(1:GIMAX,1:GJMAX),8*GIMAX*GJMAX,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(xm_global_j(1:GIMAX,1:GJMAX),8*GIMAX*GJMAX,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(glat_fdom,8*IIFULLDOM*JJFULLDOM,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(glon_fdom,8*IIFULLDOM*JJFULLDOM,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
    CALL MPI_BCAST(projection,len(projection),MPI_CHARACTER,0,MPI_COMM_WORLD,INFO) 

    if(trim(projection)=='Rotated_Spherical')then
       CALL MPI_BCAST(grid_north_pole_longitude,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
       CALL MPI_BCAST(grid_north_pole_latitude,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
       CALL MPI_BCAST(dx_rot,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
       CALL MPI_BCAST(x1_rot,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
       CALL MPI_BCAST(y1_rot,8*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
       dx_roti=1.0/dx_rot
    endif

    do j=1,MAXLJMAX
       do i=1,MAXLIMAX
          glon(i,j)=glon_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
          glat(i,j)=glat_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
       enddo
    enddo
    i0=0
    im=MAXLIMAX
    j0=0
    jm=MAXLJMAX
    if(gi0+MAXLIMAX+1+IRUNBEG-2>IIFULLDOM)im=MAXLIMAX-1!outside fulldomain
    if(gi0+0+IRUNBEG-2<1)i0=1!outside fulldomain
    if(gj0+MAXLJMAX+1+JRUNBEG-2>JJFULLDOM)jm=MAXLJMAX-1!outside fulldomain
    if(gj0+0+JRUNBEG-2<1)j0=1!outside fulldomain

    do j=j0,jm
       do i=i0,im
          x1=glon_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)
          x2=glon_fdom(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)
          x3=glon_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)
          x4=glon_fdom(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2)

          !8100=90*90; could use any number much larger than zero and much smaller than 180*180  
          if(x1*x2<-8100.0 .or. x1*x3<-8100.0 .or. x1*x4<-8100.0)then
             !Points are on both sides of the longitude -180=180              
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
             if(x3<0)x3=x3+360.0
             if(x4<0)x4=x4+360.0
          endif
          gl_stagg(i,j)=0.25*(x1+x2+x3+x4)

          gb_stagg(i,j)=0.25*(glat_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               glat_fdom(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
               glat_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)+&
               glat_fdom(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2))
       enddo
    enddo
    do j=0,j0
       do i=i0,im
          x1=gl_stagg(i,j+1)
          x2=gl_stagg(i,j+2)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i,j+1)-gb_stagg(i,j+2)
       enddo
    enddo
    do j=jm,MAXLJMAX
       do i=i0,im
          x1=gl_stagg(i,j-1)
          x2=gl_stagg(i,j-2)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i,j-1)-gb_stagg(i,j-2)
       enddo
    enddo
    do j=0,MAXLJMAX
       do i=0,i0
          x1=gl_stagg(i+1,j)
          x2=gl_stagg(i+2,j)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i+1,j)-gb_stagg(i+2,j)
       enddo
    enddo
    do j=0,MAXLJMAX
       do i=im,MAXLIMAX
          x1=gl_stagg(i-1,j)
          x2=gl_stagg(i-2,j)
          if(x1*x2<-8100.0 )then
             if(x1<0)x1=x1+360.0
             if(x2<0)x2=x2+360.0
          endif
          gl_stagg(i,j)=2*x1-x2
          gb_stagg(i,j)=2*gb_stagg(i-1,j)-gb_stagg(i-2,j)
       enddo
    enddo
    !ensure that values are within [-180,+180]]
    do j=0,MAXLJMAX
       do i=0,MAXLIMAX
          if(gl_stagg(i,j)>180.0)gl_stagg(i,j)=gl_stagg(i,j)-360.0
          if(gl_stagg(i,j)<-180.0)gl_stagg(i,j)=gl_stagg(i,j)+360.0
       enddo
    enddo

    !test if the grid is cyclicgrid:
    !The last cell + 1 cell = first cell
    Cyclicgrid=1 !Cyclicgrid
    do j=1,JJFULLDOM
       if(mod(nint(glon_fdom(GIMAX,j)+360+360.0/GIMAX),360)/=&
            mod(nint(glon_fdom(IRUNBEG,j)+360.0),360))then
          Cyclicgrid=0  !not cyclicgrid
       endif
    enddo

    if(MasterProc .and. DEBUG_GRIDVALUES)write(*,*)'CYCLICGRID:',Cyclicgrid

    !complete (extrapolate) along the four lateral sides
    do i=1,GIMAX
       xm_global_j(i,0)=1.0/(2.0/(xm_global_j(i,1))-1.0/(xm_global_j(i,2)))
       xm_global_j(i,-1)=1.0/(2.0/(xm_global_j(i,0))-1.0/(xm_global_j(i,1)))
       xm_global_j(i,GJMAX+1)=1.0/(2.0/(xm_global_j(i,GJMAX))-1.0/(xm_global_j(i,GJMAX-1)))
       xm_global_j(i,GJMAX+2)=1.0/(2.0/(xm_global_j(i,GJMAX+1))-1.0/(xm_global_j(i,GJMAX)))
       xm_global_i(i,0)=1.0/(2.0/(xm_global_i(i,1))-1.0/(xm_global_i(i,2)))
       xm_global_i(i,-1)=1.0/(2.0/(xm_global_i(i,0))-1.0/(xm_global_i(i,1)))
       xm_global_i(i,GJMAX+1)=1.0/(2.0/(xm_global_i(i,GJMAX))-1.0/(xm_global_i(i,GJMAX-1)))
       xm_global_i(i,GJMAX+2)=1.0/(2.0/(xm_global_i(i,GJMAX+1))-1.0/(xm_global_i(i,GJMAX)))
    enddo
    do j=-1,GJMAX+2
       xm_global_j(0,j)=1.0/(2.0/(xm_global_j(1,j))-1.0/(xm_global_j(2,j)))
       xm_global_j(-1,j)=1.0/(2.0/(xm_global_j(0,j))-1.0/(xm_global_j(1,j)))
       xm_global_j(GIMAX+1,j)=1.0/(2.0/(xm_global_j(GIMAX,j))-1.0/(xm_global_j(GIMAX-1,j)))
       xm_global_j(GIMAX+2,j)=1.0/(2.0/(xm_global_j(GIMAX+1,j))-1.0/(xm_global_j(GIMAX,j)))
       xm_global_i(0,j)=1.0/(2.0/(xm_global_i(1,j))-1.0/(xm_global_i(2,j)))
       xm_global_i(-1,j)=1.0/(2.0/(xm_global_i(0,j))-1.0/(xm_global_i(1,j)))
       xm_global_i(GIMAX+1,j)=1.0/(2.0/(xm_global_i(GIMAX,j))-1.0/(xm_global_i(GIMAX-1,j)))
       xm_global_i(GIMAX+2,j)=1.0/(1.0/(2*xm_global_i(GIMAX+1,j))-1.0/(xm_global_i(GIMAX,j)))
    enddo

    j=1
    i=1
    if(abs(1.5*glat_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)-0.5*glat_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2))>89.5)then
       write(*,*)'south pole' !xm is infinity
       xm_global_i(:,0)=1.0E19
       xm_global_i(:,-1)=1.0E19
    endif



    !keep only part of xm relevant to the local domain
    !note that xm has dimensions larger than local domain

    call CheckStop( MAXLIMAX+1 > limax+2, "Error in Met_ml X size definition" )
    call CheckStop( MAXLJMAX+1 > ljmax+2, "Error in Met_ml J size definition" )

    do j=0,MAXLJMAX+1
       do i=0,MAXLIMAX+1
          iglobal=gi0+i-1
          jglobal=gj0+j-1
          xm_i(i,j)=xm_global_i(iglobal,jglobal)
          xm_j(i,j)=xm_global_j(iglobal,jglobal)
          !Note that xm is inverse length: interpolate 1/xm rather than xm
          xm2(i,j) = 4.0*( (xm_global_i(iglobal,jglobal-1)*&
               xm_global_i(iglobal,jglobal))/ &
               (xm_global_i(iglobal,jglobal-1)+&
               xm_global_i(iglobal,jglobal)    )   )   *(&
               xm_global_j(iglobal-1,jglobal)*&
               xm_global_j(iglobal,jglobal)     )/(&
               xm_global_j(iglobal-1,jglobal)+&
               xm_global_j(iglobal,jglobal)     )
       enddo
    enddo


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


  end subroutine Getgridparams




  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine DefGrid()
    !-------------------------------------------------------------------! 
    !     defines map parameters and fields for the model.
    !-------------------------------------------------------------------! 

    integer ::  i, j, k, n
    real    ::  an2, x, y, x_j, y_i
    real    ::  rpol2,rpol2_i,rpol2_j ! square of (distance from pole to i,j
    ! divided by AN )
    real, dimension(0:MAXLIMAX+1,0:MAXLJMAX+1) ::  &
         xm      ! map-factor 

    ! Earth radius = 6.37e6 m, gives gridwidth:

    !   GRIDWIDTH_M = 6.370e6*(1.0+0.5*sqrt(3.0))/AN    ! = 50000.0 m


    ! NB! HIRLAM uses Earth radius = 6.371e6 m : 
    ! AN = No. grids from pole to equator
    ! AN = 6.371e6*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M = 237.768957  

    !/ Define full-domain coordinates of local i,j values. We need to account for
    !  the fact that each parallel domain has its starting cordinate
    !  gi0, gj0, and the user may specify a set of lower-left starting
    !  coordinates for running the model, IRUNBEG, JRUNBEG
    !       i_fdom(i)  = i + gi0 + IRUNBEG - 2
    !       j_fdom(j)  = j + gj0 + JRUNBEG - 2

    i_fdom = (/ (n + gi0 + IRUNBEG - 2, n=0,MAXLIMAX+1) /) 
    j_fdom = (/ (n + gj0 + JRUNBEG - 2, n=0,MAXLJMAX+1) /) 

    ! And the reverse, noting that we even define for area
    ! outside local domain

    i_local = (/ (n - gi0 - IRUNBEG + 2, n=1, IIFULLDOM) /)
    j_local = (/ (n - gj0 - JRUNBEG + 2, n=1, JJFULLDOM) /)


    !------------------------------------------------------------------


    AN = 6.370e6*(1.0+sin( ref_latitude*PI/180.))/GRIDWIDTH_M ! = 237.7316364 for GRIDWIDTH_M=50 km and ref_latitude=60
    do j=0,MAXLJMAX+1
       do i=0,MAXLIMAX+1
          xmd(i,j) = 1.0/xm2(i,j)
          xm2ji(j,i) = xm2(i,j)
          xmdji(j,i) = xmd(i,j)
       enddo
    enddo
    do j=1,MAXLJMAX
       do i=1,MAXLIMAX
          GridArea_m2(i,j) = GRIDWIDTH_M*GRIDWIDTH_M*xmd(i,j)
       enddo
    enddo

    ! definition of the half-sigma levels (boundaries between layers)
    ! from the full levels. 

    sigma_bnd(KMAX_BND) = 1.
    do k = KMAX_MID,2,-1
       sigma_bnd(k) = 2.*sigma_mid(k) - sigma_bnd(k+1)
    enddo
    sigma_bnd(1) = 0.

    !TEMPORARY: definition of A and B. Will be read from metfile in the future
    do k = 1,KMAX_BND
       A_bnd(k)=PT * (1-sigma_bnd(k)) 
       B_bnd(k)=sigma_bnd(k) 
    enddo
    do k = 1,KMAX_MID
       dA(k)=A_bnd(k+1)-A_bnd(k)
       dB(k)=B_bnd(k+1)-B_bnd(k)
       A_mid(k)=(A_bnd(k+1)+A_bnd(k))/2.0      
       B_mid(k)=(B_bnd(k+1)+B_bnd(k))/2.0
    enddo


    !
    !     some conversion coefficients needed for budget calculations
    !
    do  k=1,KMAX_MID
       carea(k) = (sigma_bnd(k+1) - sigma_bnd(k))/GRAV*GRIDWIDTH_M*GRIDWIDTH_M
       !write(6,*)'carea,sigma_bnd,h',carea(k),sigma_bnd(k),GRIDWIDTH_M
       !su  cflux(k) = (sigma_bnd(k+1) - sigma_bnd(k))/G*h
    end do

    ! set latitude, longitude
    ! projection='Stereographic'
    call Position()

    if ( DEBUG_GRIDVALUES ) then
       if ( me == 0 ) then
          write(*,800) "GRIDTAB","me","ISM","JSM","gi0","gj0",&
               "li0","li1","lix","MXI",&
               "lj0","lj1","ljx"," MXJ"," ig1"," igX","jg1","jgX"
          write(*,802) "GRIDLL ","me", "mingl"," maxgl"," mingb"," maxgb",&
               " glat(1,1)"," glat(MAX..)"
       end if

       write(*,804) "GRIDTAB",me,IRUNBEG,JRUNBEG,gi0,gj0,li0,li1,&
            limax,MAXLIMAX,lj0,lj1,ljmax, MAXLJMAX, i_fdom(1),&
            i_fdom(MAXLIMAX+1),j_fdom(1), j_fdom(MAXLJMAX+1)

       write(*,806) "GRIDLL ",me, minval(glon), maxval(glon), minval(glat), &
            maxval(glat), glat(1,1), glat(MAXLIMAX,MAXLJMAX)
    end if
800 format(a10,20a4)
802 format(a10,a4,10a12)
804 format(a10,20i4)
806 format(a10,i4,10f12.4)

  end subroutine DefGrid
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine DefDebugProc()
    !-------------------------------------------------------------------! 
    ! -------------- Find debug coords  and processor ------------------
    !-------------------------------------------------------------------! 

    integer ::  i, j

    debug_proc = .false.

    do i = li0, li1
       do j = lj0, lj1
          if( i_fdom(i) == DEBUG_i .and. j_fdom(j) == DEBUG_j ) then
             debug_li = i
             debug_lj = j
             debug_proc = .true.
          end if
       end do
    end do

    if( debug_proc ) write(*,*) "GridValues debug_proc found:", &
         me, debug_li, debug_lj
    if ( DEBUG_GRIDVALUES ) then
       if(me==0) write(*,"(a,2a4,a3,4a4,a2,2a4,4a12)") "GridValues debug:", &
            "D_i", "D_j", "me", "li0", "li1", "lj0", "lj1", &
            "dp" , "d_li", "d_lj", "i_fdom(li0)","i_fdom(li1)", &
            "j_fdom(lj0)", "j_fdom(lj1)"

       write(*,"(a,2i4,i3,4i4,L2,2i4,4i12)") "GridValues debug:", &
            DEBUG_i, DEBUG_j, me, li0, li1, lj0, lj1, &
            debug_proc , debug_li, debug_lj, &
            i_fdom(li0),i_fdom(li1), j_fdom(lj0), j_fdom(lj1)
    end if

  end subroutine DefDebugProc
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine Position()
    !-------------------------------------------------------------------! 
    !      calculates l(lat),b(long) (geographical coord.) 
    !      in every grid point defined by the i_fdom, j_fdom arrays. 
    !
    !      input:  xp,yp:   coord. of the polar point.
    !              AN:      number of grid-distances from pole to equator.
    !              fi:      rotational angle for the x,y grid.
    !              imax,jmax:   number of points in  x- og y- direction
    !              glmin:   gives min.value of geographical lenght
    !                       =>  glmin <= l <= glmin+360.  
    !                           (example glmin = -180. or 0.)
    !                       if "geopos","georek" is used
    !                       then glmin must be the lenght i(1,1) in the
    !                       geographical grid (gl1 to "geopos")
    !      output: gl(ii,jj): latitude glmin <= l <= glmin+360. 
    !              gb(ii,jj): longitude  -90. <= b <= +90. 
    !-------------------------------------------------------------------! 
    !   - evaluate gl, gb over whole domain given by MAXLIMAX, MAXLJMAX
    !       to safeguard against possible use of non-defined gl,bb squares. 
    !   - note, we could use rpol2(i,j) to save some computations here, 
    !       but for now we leave it. This stuff is only done once anyway

    real    :: glmin, glmax, om, om2, dy, dy2,rp,rb, rl, dx, dr
    integer :: i, j, info

    !su    xp,yp read in infield                
    !su    xp = 43.
    !su    yp = 121.

    glmin = -180.0
    glmax = glmin + 360.0

    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees
    om2   = om * 2.0


    if(trim(projection)=='Stereographic') then     

       do j = 1, MAXLJMAX          !  - changed from ljmax
          dy  = yp - j_fdom(j)     !  - changed from gj0+JRUNBEG-2+j
          dy2 = dy*dy
          do i = 1, MAXLIMAX       !  - changed from limax
             dx = i_fdom(i) - xp    !  - changed
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
    endif

    ! test to find full-domain min and max lat/long values

    gbacmax = maxval(glat(:,:))
    gbacmin = minval(glat(:,:))
    glacmax = maxval(glon(:,:))
    glacmin = minval(glon(:,:))
    MPIbuff= gbacmax 
    CALL MPI_ALLREDUCE(MPIbuff, gbacmax, 1,MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, INFO) 
    MPIbuff= gbacmin
    CALL MPI_ALLREDUCE(MPIbuff, gbacmin  , 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO) 
    MPIbuff= glacmax 
    CALL MPI_ALLREDUCE(MPIbuff, glacmax, 1,MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, INFO) 
    MPIbuff= glacmin
    CALL MPI_ALLREDUCE(MPIbuff, glacmin  , 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO) 

    if(me==0) write(unit=6,fmt="(a,4f9.2)") &
         " GridValues: max/min for lat,lon ", &
         gbacmax,gbacmin,glacmax,glacmin

    if ( DEBUG_GRIDVALUES ) then
       do j = 1, MAXLJMAX
          do i = 1, MAXLIMAX
             if ( i_fdom(i) == DEBUG_i .and. j_fdom(j) == DEBUG_j ) then
                write(*,"(a15,a30,5i4,2f8.2,f7.3)") "DEBUGPosition: ",  &
                     " me,i,j,i_fdom,j_fdom,glon,glat,rp: ", &
                     me, i,j, i_fdom(i), j_fdom(j), glon(i,j), glat(i,j),rp
             end if
          end do
       end do
    end if

  end subroutine Position
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine GlobalPosition

    integer i,j
    real :: dr,om,om2,rb,rl,rp,dx,dy,dy2,glmax,glmin
    integer :: im,jm,i0,j0

    if(trim(projection)=='Stereographic') then     

       glmin = -180.0
       glmax = glmin + 360.0

       dr    = PI/180.0      ! degrees to radians
       om    = 180.0/PI      ! radians to degrees
       om2   = om * 2.0
       AN = 6.370e6*(1.0+sin( ref_latitude*PI/180.))/GRIDWIDTH_M ! = 237.7316364 for GRIDWIDTH_M=50 km and ref_latitude=60

       do j = 1, JJFULLDOM
          dy  = yp - j  
          dy2 = dy*dy
          do i = 1, IIFULLDOM
             dx = i - xp    
             rp = sqrt(dx*dx+dy2)           ! => distance to pole
             rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
             rl = 0.0
             if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
             if (rl <  glmin)   rl = rl + 360.0
             if (rl >  glmax)   rl = rl - 360.0
             glon_fdom(i,j)=rl              ! longitude
             glat_fdom(i,j)=rb              ! latitude

          end do ! i
       end do ! j

       i0=0
       im=MAXLIMAX
       j0=0
       jm=MAXLJMAX
       if(gi0+MAXLIMAX+1+IRUNBEG-2>IIFULLDOM)im=MAXLIMAX-1!outside fulldomain
       if(gi0+0+IRUNBEG-2<1)i0=1!outside fulldomain
       if(gj0+MAXLJMAX+1+JRUNBEG-2>JJFULLDOM)jm=MAXLJMAX-1!outside fulldomain
       if(gj0+0+JRUNBEG-2<1)j0=1!outside fulldomain
       do j=j0,jm
          do i=i0,im
             gl_stagg(i,j)=0.25*(glon_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
                  glon_fdom(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
                  glon_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)+&
                  glon_fdom(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2))
             gb_stagg(i,j)=0.25*(glat_fdom(gi0+i+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
                  glat_fdom(gi0+i+1+IRUNBEG-2,gj0+j+JRUNBEG-2)+&
                  glat_fdom(gi0+i+IRUNBEG-2,gj0+j+1+JRUNBEG-2)+&
                  glat_fdom(gi0+i+1+IRUNBEG-2,gj0+j+1+JRUNBEG-2))
          enddo
       enddo
       do j=0,j0
          do i=i0,im
             gl_stagg(i,j)=2*gl_stagg(i,j+1)-gl_stagg(i,j+2)
             gb_stagg(i,j)=2*gb_stagg(i,j+1)-gb_stagg(i,j+2)
          enddo
       enddo
       do j=jm,MAXLJMAX
          do i=i0,im
             gl_stagg(i,j)=2*gl_stagg(i,j-1)-gl_stagg(i,j-2)
             gb_stagg(i,j)=2*gb_stagg(i,j-1)-gb_stagg(i,j-2)
          enddo
       enddo
       do j=0,MAXLJMAX
          do i=0,i0
             gl_stagg(i,j)=2*gl_stagg(i+1,j)-gl_stagg(i+2,j)
             gb_stagg(i,j)=2*gb_stagg(i+1,j)-gb_stagg(i+2,j)
          enddo
       enddo
       do j=0,MAXLJMAX
          do i=im,MAXLIMAX
             gl_stagg(i,j)=2*gl_stagg(i-1,j)-gl_stagg(i-2,j)
             gb_stagg(i,j)=2*gb_stagg(i-1,j)-gb_stagg(i-2,j)
          enddo
       enddo


    endif

  end subroutine GlobalPosition


  subroutine lb2ijm(imax,jmax,glon,glat,xr2,yr2,fi2,an2,xp2,yp2)
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

    real, intent(in)    :: glon(imax,jmax),glat(imax,jmax) 
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
             xr2(i1,j1)=xp_loc+an_loc*tan(PId4-glat(i1,j1)*dr2)*sin(dr*(glon(i1,j1)-fi_loc))
             yr2(i1,j1)=yp_loc-an_loc*tan(PId4-glat(i1,j1)*dr2)*cos(dr*(glon(i1,j1)-fi_loc))
          enddo
       enddo
    else if(projection=='lon lat')then! lon-lat grid
       do j1 = 1, jmax
          do i1 = 1, imax
             xr2(i1,j1)=(glon(i1,j1)-glon_fdom(1,1))/(glon_fdom(2,1)-glon_fdom(1,1))+1
             if(xr2(i1,j1)<0.5)xr2=xr2+360.0/(glon_fdom(2,1)-glon_fdom(1,1))
             yr2(i1,j1)=(glat(i1,j1)-glat_fdom(1,1))/(glat_fdom(1,2)-glat_fdom(1,1))+1
          enddo
       enddo
    else!general projection, Use only info from glon_fdom and glat_fdom
       call StopAll('lb2ijm: conversion not yet tested. (You could try the method below)')

       !VERY SLOW, specially for large grids
       do j1 = 1, jmax
          do i1 = 1, imax
             dist=10.0!max distance is PI
             do j=1,JJFULLDOM
                do i=1,IIFULLDOM
                   if(dist>great_circle_distance(glon(i1,j1),glat(i1,j1),glon_fdom(i,j) &
                        ,glat_fdom(i,j)))then
                      dist=great_circle_distance(glon(i1,j1),glat(i1,j1),glon_fdom(i,j) &
                           ,glat_fdom(i,j))
                      xr2(i1,j1)=i
                      yr2(i1,j1)=j
                   endif
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
             dist2=great_circle_distance(glon(i1,j1),glat(i1,j1),glon_fdom(ip1,jr2),glat_fdom(ip1,jr2))
             dist3=great_circle_distance( glon_fdom(ir2,jr2), &
                  glat_fdom(ir2,jr2), &
                  glon_fdom(ip1,jr2), &
                  glat_fdom(ip1,jr2))

             xr2(i1,j1)=xr2(i1,j1)+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)


             jp1=jr2+1
             if(jp1>JJFULLDOM)jp1=jp1-2

             dist2=great_circle_distance(glon(i1,j1),glat(i1,j1),glon_fdom(ir2,jp1),glat_fdom(ir2,jp1))
             !GFORTRAN CHANGE
             dist3=great_circle_distance( glon_fdom(ir2,jr2), &
                  glat_fdom(ir2,jr2), &
                  glon_fdom(ir2,jp1), & 
                  glat_fdom(ir2,jp1) )

             yr2(i1,j1)=yr2(i1,j1)+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)

          enddo
       enddo

    endif
  end subroutine lb2ijm

  subroutine lb2ij(gl2,gb2,xr2,yr2,fi2,an2,xp2,yp2)

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
    real ::xscen ,yscen,zsycen,zcycen ,zxmxc,zsxmxc,zcxmxc,zsysph,zsyrot,yrot,zsxrot,zcysph,zcyrot,zcxrot,xrot,dx,x1,y1

    if(projection=='Stereographic')then
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

    else  if(projection=='lon lat')then! lon-lat grid

       xr2=(gl2-glon_fdom(1,1))/(glon_fdom(2,1)-glon_fdom(1,1))+1
       if(xr2<0.5)xr2=xr2+360.0/(glon_fdom(2,1)-glon_fdom(1,1))
       yr2=(gb2-glat_fdom(1,1))/(glat_fdom(1,2)-glat_fdom(1,1))+1

    else  if(projection=='Rotated_Spherical')then! rotated lon-lat grid
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
         dx=0.05
         x1=-13.65
         y1=-1.027
         xr2=(xrot-x1_rot)*dx_roti+1
         yr2=(yrot-y1_rot)*dx_roti+1

    else!general projection, Use only info from glon_fdom and glat_fdom
       !first find closest by testing all gridcells. 

       dist=10.0!max distance is PI
       do j=1,JJFULLDOM
          do i=1,IIFULLDOM
             if(dist>great_circle_distance(gl2,gb2,glon_fdom(i,j) &
                  ,glat_fdom(i,j)))then
                dist=great_circle_distance(gl2,gb2,glon_fdom(i,j) &
                     ,glat_fdom(i,j))
                xr2=i
                yr2=j
             endif
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
       dist2=great_circle_distance(gl2,gb2,glon_fdom(ip1,jr2),glat_fdom(ip1,jr2))
       dist3=great_circle_distance( glon_fdom(ir2,jr2), &
            glat_fdom(ir2,jr2), &
            glon_fdom(ip1,jr2), &
            glat_fdom(ip1,jr2))

       xr2=xr2+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)


       jp1=jr2+1
       if(jp1>JJFULLDOM)jp1=jp1-2

       dist2=great_circle_distance(gl2,gb2,glon_fdom(ir2,jp1),glat_fdom(ir2,jp1))
       !GFORTRAN CHANGE
       dist3=great_circle_distance( glon_fdom(ir2,jr2), &
            glat_fdom(ir2,jr2), &
            glon_fdom(ir2,jp1), & 
            glat_fdom(ir2,jp1) )

       yr2=yr2+(dist*dist+dist3*dist3-dist2*dist2)/(2*dist3*dist3)

    endif

    return
  end subroutine lb2ij

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
    if ( alloc_err /= 0 ) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 

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
       call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
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
    if ( alloc_err /= 0 ) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 

  end subroutine ij2ijm

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  function coord_in_gridbox(lon,lat,i,j) result(in)
    !-------------------------------------------------------------------!
    ! Is coord (lon/lat) is inside gridbox(i,j)?
    !-------------------------------------------------------------------!
    real, intent(in) :: lon,lat
    integer, intent(in) :: i,j
    logical :: in
    in=gl_stagg(i-1,j)<=lon.and.lon<gl_stagg(i,j).and. &
         gb_stagg(i,j-1)<=lat.and.lat<gb_stagg(i,j)
  end function coord_in_gridbox
  function coord_in_processor(lon,lat) result(in)
    !-------------------------------------------------------------------!
    ! Is coord (lon/lat) is inside local domain?
    !-------------------------------------------------------------------!
    real, intent(in) :: lon,lat
    logical :: in
    integer :: i,j
    real    :: xr,yr
    call lb2ij(lon,lat,xr,yr)
    i=nint(xr);j=nint(yr)
    in=(gi0<=i).and.(i<=gi1).and.(gj0<=j).and.(j<=gj1)
  end function coord_in_processor
  subroutine check(status)
    use netcdf
    implicit none
    integer, intent ( in) :: status

    call CheckStop( status, nf90_noerr, "Error in GridValues_ml/NetCDF parts:"  &
         //  trim( nf90_strerror(status) ) )

  end subroutine check

  subroutine Alloc_GridFields(GIMAX,GJMAX,MAXLIMAX,MAXLJMAX,KMAX_MID,KMAX_BND)

    integer, intent(in)::GIMAX,GJMAX,MAXLIMAX,MAXLJMAX,KMAX_MID,KMAX_BND

    allocate(i_fdom(0:MAXLIMAX+1))
    allocate(j_fdom(0:MAXLJMAX+1))
    allocate(glon(MAXLIMAX,MAXLJMAX))
    allocate(glat(MAXLIMAX,MAXLJMAX))
    allocate(gl_stagg(0:MAXLIMAX,0:MAXLJMAX))
    allocate(gb_stagg(0:MAXLIMAX,0:MAXLJMAX))
    allocate(xm_i(0:MAXLIMAX+1,0:MAXLJMAX+1))
    allocate(xm_j(0:MAXLIMAX+1,0:MAXLJMAX+1))
    allocate(xm2(0:MAXLIMAX+1,0:MAXLJMAX+1))
    allocate(xmd(0:MAXLIMAX+1,0:MAXLJMAX+1))
    allocate(xm2ji(0:MAXLJMAX+1,0:MAXLIMAX+1))
    allocate(xmdji(0:MAXLJMAX+1,0:MAXLIMAX+1))
    allocate(GridArea_m2(MAXLIMAX,MAXLJMAX))
    allocate(z_bnd(MAXLIMAX,MAXLJMAX,KMAX_BND))
    allocate(z_mid(MAXLIMAX,MAXLJMAX,KMAX_MID))


  end subroutine Alloc_GridFields

end module GridValues_ml
!==============================================================================
