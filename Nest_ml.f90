! <Nest_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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
  !
  !This module performs the reading or writing of data for nested runs
  !
  !To make a nested run:
  !1) run with MODE=1 to write out 3d BC
  !2) run (in a smaller domain) with MODE=2


  !
  !Set MODE and istart,jstart,iend,jend
  !Choose NHOURSAVE and NHOURREAD 
  !


  !Grids may have any projection.
  !Horizontal interpolation uses a weighted average of the four closest points
  !This will work also if points in the present grid are not covered by the external grid.


  !To do:
  !At present the vertical coordinates cannot be interpolated and must be the same in both grid.
  !It should be possible to save only xn_adv_bnd if the inner grid is known for the outer grid.
  !The routines should be thought together with GlobalBC_ml (can it replace it?)


  !Peter May 2006

!hfTD  use Dates_ml,       only : date   
  use TimeDate_ml,       only : date   
  use GridValues_ml,  only : gl,gb
  use GenChemicals_ml , only :species
  use GenSpec_shl_ml , only :NSPEC_SHL
  use GenSpec_adv_ml , only :NSPEC_ADV
  use GenSpec_tot_ml , only :NSPEC_TOT
  use netcdf
  use netcdf_ml,      only : GetCDF,Out_netCDF,Init_new_netCDF,&
       secondssince1970,Int1,Int2,Int4,Real4,Real8
  use ModelConstants_ml,    only : KMAX_MID, NPROC
  use Par_ml   ,      only : MAXLIMAX, MAXLJMAX, GIMAX,GJMAX,IRUNBEG,JRUNBEG &
       , me, li0,li1,lj0,lj1,limax,ljmax, tgi0, tgj0, tlimax, tljmax
  use Chemfields_ml,  only : xn_adv, xn_shl    ! emep model concs.

  use NetCDF_ml, only :WriteCDF

  implicit none

  INCLUDE 'mpif.h'
  INTEGER INFO

  integer,parameter ::MODE=0   !0=donothing , 1=write , 2=read , 3=read and write
  !10=write at end of run, 11=read at start , 12=read at start and write at end (BIC)

  !coordinates of subdomain to write
  !coordinates relative to small domain (only used in write mode)
  integer ::istart=60,jstart=11,iend=107,jend=58 !ENEA

  !/-- subroutines

  public  :: readxn
  public  :: wrtxn


  !  logical, save, public::Nest_BC,Nest_3D

  integer, public, parameter :: NHOURSAVE=3 !time between two saves. should be a fraction of 24
  integer, public, parameter :: NHOURREAD=1 !time between two reads. should be a fraction of 24
  !if(NHOURREAD<NHOURSAVE) the data is interpolated in time

  private

  !BC values at boundaries in present grid
  real,save, dimension(NSPEC_ADV,MAXLIMAX,KMAX_MID,2) :: xn_adv_bnds,xn_adv_bndn !north and south
  real,save, dimension(NSPEC_ADV,MAXLJMAX,KMAX_MID,2) :: xn_adv_bndw,xn_adv_bnde !west and east

  !4 nearest points from external grid
  integer, save ::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)

  !weights of the 4 nearest points
  real, save :: Weight1(MAXLIMAX,MAXLJMAX),Weight2(MAXLIMAX,MAXLJMAX)
  real, save :: Weight3(MAXLIMAX,MAXLJMAX),Weight4(MAXLIMAX,MAXLJMAX)

  !dimension of external grid
  integer,save :: Next,KMAX_ext,GJMAX_ext,GIMAX_ext

  integer,save :: itime_saved(2),itime
  character*30,save  :: filename_read='EMEP_IN.nc' 
  character*30,save  :: filename_write='EMEP_OUT.nc' 

contains

  subroutine readxn(indate)
    type(date), intent(in) :: indate           ! Gives year..seconds
    integer,save  :: first_data=-1


    integer :: nseconds(1),n1,n,i,j,k,II,JJ
    integer :: nstart,nfetch,ndate(4),nseconds_indate

    !    real , dimension(48,48,20) ::data

    real :: W1,W2
    logical, save :: first_call=.true.

    if(MODE /= 2.and.MODE /= 3.and. MODE /= 11.and. MODE /= 12)return

    ndate(1)  = indate%year
    ndate(2)  = indate%month
    ndate(3)  = indate%day
    ndate(4)  = indate%hour
    call secondssince1970(ndate,nseconds_indate)

    if(MODE == 11.or.MODE == 12)then
       if(.not. first_call)return
       first_call=.false.
       if(me==0)   print *,'RESET ALL XN 3D'
       call init_nest(nseconds_indate)
       call reset_3D(nseconds_indate)
       return
    else
!    if(me==0)   print *,'call to READXN',indate%hour,indate%seconds
       if(mod(indate%hour,NHOURREAD)/=0.or.indate%seconds/=0)return
    endif
    if(me==0)   print *,'NESTING'

    if(first_data==-1)then
       call init_nest(nseconds_indate)
       call reset_3D(nseconds_indate)
       call read_newdata_LATERAL(nseconds_indate,1)
       call read_newdata_LATERAL(nseconds_indate,2)
    endif


    if(itime_saved(2)<nseconds_indate)then
       !look for a new data set
       if(me==0)write(*,*)'NEST: READING NEW BC DATA'
       call read_newdata_LATERAL(nseconds_indate,2)
    endif


    !    make weights for time interpolation
    W1=1.0;  W2=0.0!default
    if(itime_saved(1)<nseconds_indate)then
       !interpolate
       W2=(nseconds_indate-itime_saved(1))/(1.0*itime_saved(2)-itime_saved(1))
       W1=1.0-W2
!       call datefromsecondssince1970(ndate,nseconds_indate,1)
!       write(*,*)'interpolating between'
!       call datefromsecondssince1970(ndate,itime_saved(1),1)
!       write(*,*)'and'
!       call datefromsecondssince1970(ndate,itime_saved(2),1)
!       write(*,*)'with weights : ',W1,W2

    endif
!    if(me==0)write(*,*)'weights : ',W1,W2,itime_saved(1),itime_saved(2)

    do n=1,NSPEC_ADV
       do k=1,KMAX_ext
          do i=1,li0-1
             do j=1,ljmax
                xn_adv(n,i,j,k)=W1*xn_adv_bndw(n,j,k,1)+W2*xn_adv_bndw(n,j,k,2)
             enddo
          enddo
          do j=1,lj0-1
             do i=1,limax 
                xn_adv(n,i,j,k)=W1*xn_adv_bnds(n,i,k,1)+W2*xn_adv_bnds(n,i,k,2)
             enddo
          enddo
          do i=li1+1,limax
             do j=1,ljmax
                xn_adv(n,i,j,k)=W1*xn_adv_bnde(n,j,k,1)+W2*xn_adv_bnde(n,j,k,2)
             enddo
          enddo
          do j=lj1+1,ljmax
             do i=1,limax
                xn_adv(n,i,j,k)=W1*xn_adv_bndn(n,i,k,1)+W2*xn_adv_bndn(n,i,k,2)
             enddo
          enddo
       enddo
    enddo


    first_data=0


    return


  end subroutine readxn

  subroutine wrtxn(indate,WriteNow)

    use Derived_ml, only :Deriv 
    use Derived_ml,    only :IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY  
    implicit none

    type(date), intent(in) :: indate      
    logical, intent(in) :: WriteNow !Do not check indate value 
    real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: dat ! Data arrays

    type(Deriv) :: def1 ! definition of fields
    integer :: n,iotyp,ndim,kmax
    real :: scale
    logical, save ::first_call=.true.
    character *40:: command

    if(MODE /= 1.and.MODE /= 3.and.MODE /= 10.and.MODE /= 12)return

    if(MODE == 10.or.MODE == 12)then
       if(.not.WriteNow)return
       istart=1
       jstart=1
       iend=GIMAX
       jend=GJMAX
    else
       if(mod(indate%hour,NHOURSAVE)/=0.or.indate%seconds/=0)return
    endif


    222 FORMAT(A,I2.2,I4.4,A)
!    write(fileName_write,222)'EMEP_BC_',indate%month,indate%year,'.nc'!for different names each month
                                                                      !NB: readxn should have same name
    if(me==0)write(*,*)'write Nest data',trim(fileName_write)



    iotyp=IOU_INST
    if(first_call)then
       if(me==0)then
          write(*,*)'cleaning possible old ',trim(fileName_write)
          write(command,*)'rm ',trim(fileName_write)
          call system(command)       
       endif
       first_call=.false.
    endif

    ndim=3 !3-dimensional
    kmax=KMAX_MID
    scale=1.0
    def1%class='Advected' !written
    def1%avg=.false.      !not used
    def1%index=0          !not used
    def1%scale=scale      !not used
    def1%rho=.false.      !not used
    def1%inst=.true.      !not used
    def1%year=.false.     !not used
    def1%month=.false.    !not used
    def1%day=.false.      !not used
    def1%name=''        !written
    def1%unit='mix_ratio'       !written

    do n= 1, NSPEC_ADV
!    do n= 1, NSPEC_ADV-4  !ENEA
       def1%name= species(NSPEC_SHL+n)%name       !written
       dat=xn_adv(n,:,:,:)
       call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real4,&
            ist=istart,jst=jstart,ien=iend,jen=jend&
            ,fileName_given=fileName_write)
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

  function great_circle_distance(fi1,lambda1,fi2,lambda2) result(dist)

    !compute the great circle distance between to points given in 
    !spherical coordinates. Sphere has radius 1.
    real, intent(in) ::fi1,lambda1,fi2,lambda2 !NB: in DEGREES here
    real :: dist


    dist=2*asin(sqrt(sind(0.5*(lambda1-lambda2+360.0))**2+&
         cosd(lambda1+360.0)*cosd(lambda2+360.0)*sind(0.5*(fi1-fi2+360.0))**2))

  end function great_circle_distance


  subroutine datefromsecondssince1970(ndate,nseconds,printdate)
    !calculate date from seconds that have passed since the start of the year 1970

    !  use Dates_ml, only : nmdays
    implicit none

    integer, intent(out) :: ndate(4)
    integer, intent(in) :: nseconds
    integer,  intent(in) :: printdate

    integer :: n,nday,nmdays(12),nmdays2(13)
    nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/) 

    nmdays2(1:12)=nmdays
    nmdays2(13)=0
    ndate(1)=1969
    n=0
    do while(n<=nseconds)
       n=n+24*3600*365
       ndate(1)=ndate(1)+1
       if(mod(ndate(1),4)==0)n=n+24*3600
    enddo
    n=n-24*3600*365
    if(mod(ndate(1),4)==0)n=n-24*3600
    if(mod(ndate(1),4)==0)nmdays2(2)=29
    ndate(2)=0
    do while(n<=nseconds)
       ndate(2)=ndate(2)+1
       n=n+24*3600*nmdays2(ndate(2))
    enddo
    n=n-24*3600*nmdays2(ndate(2))
    ndate(3)=0
    do while(n<=nseconds)
       ndate(3)=ndate(3)+1
       n=n+24*3600
    enddo
    n=n-24*3600
    ndate(4)=-1
    do while(n<=nseconds)
       ndate(4)=ndate(4)+1
       n=n+3600
    enddo
    n=n-3600
    !    ndate(5)=nseconds-n
    if(printdate>0)then
       write(*,55)'year: ',ndate(1),', month: ',ndate(2),', day: ',&
            ndate(3),', hour: ',ndate(4),', seconds: ',nseconds-n
    endif
55  format(A,I5,A,I4,A,I4,A,I4,A,I10)
  end subroutine datefromsecondssince1970

  subroutine init_nest(nseconds_indate)

    implicit none
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID,status
    integer :: nseconds_indate,ndate(4)
    real :: dist(0:4)
    integer :: nseconds(1),n1,n,i,j,k,II,JJ  
    real, allocatable, dimension(:,:) ::lon_ext,lat_ext
    character*80 ::projection

    itime_saved = -999999 !initialization

    !Read dimensions (global)
    if(me==0)then
       status = nf90_open(path=trim(filename_read),mode=nf90_nowrite,ncid=ncFileID)

       if(status /= nf90_noerr) then     
          print *,'not found',trim(filename_read)
          return
       else
          print *,'  reading ',trim(filename_read)
       endif
       projection=''
       call check(nf90_get_att(ncFileID,nf90_global,"projection",projection))
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
          !       WRITE(*,*) 'MPI_ABORT: ', "PROJECTION NOT RECOGNIZED" 
            !     call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
       endif

       call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
       call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_ext))
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_ext))
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_ext))
       call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=Next))

       write(*,*)'dimensions external grid',GIMAX_ext,GJMAX_ext,KMAX_ext,Next
    endif
      CALL MPI_BCAST(GIMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
      CALL MPI_BCAST(GJMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
      CALL MPI_BCAST(KMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
      CALL MPI_BCAST(Next,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 

    allocate(lon_ext(GIMAX_ext,GJMAX_ext))
    allocate(lat_ext(GIMAX_ext,GJMAX_ext))

    if(me==0)then
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
       call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))
       
       !          do n=1,Next
       call check(nf90_get_var(ncFileID, varID, nseconds,start=(/ 1 /),count=(/ 1 /) ))
       
       if(nseconds(1)>nseconds_indate)then
          write(*,*)'WARNING: did not find BIC for date:'
          call datefromsecondssince1970(ndate,nseconds_indate,1)
          write(*,*)'first date found:'
          call datefromsecondssince1970(ndate,nseconds(1),1)
       endif
       !          enddo
       
       call check(nf90_close(ncFileID))

    endif

      CALL MPI_BCAST(lon_ext,8*GIMAX_ext*GJMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
      CALL MPI_BCAST(lat_ext,8*GIMAX_ext*GJMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 

    !find interpolation constants
    !note that i,j are local 
    !find the four closest points
    do j=1,ljmax
       do i=1,limax
          dist=1.0E40
          do JJ=1,GJMAX_ext
             do II=1,GIMAX_ext
                !distance between (i,j) and (II,JJ)
                dist(0)=great_circle_distance(lon_ext(II,JJ),lat_ext(II,JJ),gl(i,j),gb(i,j))
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
          Weight1(i,j)=1.0-3.0*dist(1)/dist(0)
          dist(0)=(dist(2)+dist(3)+dist(4))
          Weight2(i,j)=(1.0-Weight1(i,j))*(1.0-2.0*dist(2)/dist(0))
          dist(0)=(dist(3)+dist(4))
          Weight3(i,j)=(1.0-Weight1(i,j)-Weight2(i,j))*(1.0-dist(3)/dist(0))
          Weight4(i,j)=1.0-Weight1(i,j)-Weight2(i,j)-Weight3(i,j)

       enddo
    enddo

    deallocate(lon_ext)
    deallocate(lat_ext)

  end subroutine init_nest

  subroutine read_newdata_LATERAL(nseconds_indate,nr)

    implicit none
    real, allocatable, dimension(:,:,:) ::data
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID,status
    integer :: nseconds(1),ndate(4),n1,n,i,j,k,II,JJ,nseconds_indate,nr  
    integer :: nseconds_old

    nseconds_old=itime_saved(2)
    allocate(data(GIMAX_ext,GJMAX_ext,KMAX_ext), stat=status)
    if(me==0)then 
       call check(nf90_open(path = trim(fileName_read), mode = nf90_nowrite, ncid = ncFileID))

       call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))
       do n=1,Next
          call check(nf90_get_var(ncFileID, varID, nseconds,start=(/ n /),count=(/ 1 /) ))
          if(nseconds(1)>=nseconds_indate)then
             write(*,*)'Using date '
             call datefromsecondssince1970(ndate,nseconds(1),1)
             goto 876
          endif
       enddo
       write(*,*)'WARNING: did not find correct date'
       itime=Next
       write(*,*)'Using date '
       call datefromsecondssince1970(ndate,nseconds(1),1)
876    continue
       itime=n
       itime_saved(2)=nseconds(1)
    endif

      CALL MPI_BCAST(itime_saved,4*2,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 

    do n= 1, NSPEC_ADV
       if(nr==2)then
          !store the old vaules in 1
          itime_saved(1)=nseconds_old
          do k=1,KMAX_ext
             do i=1,li0-1
                do j=1,ljmax
                   xn_adv_bndw(n,j,k,1)=xn_adv_bndw(n,j,k,2)
                enddo
             enddo
             do j=1,lj0-1
                do i=1,limax 
                   xn_adv_bnds(n,i,k,1)=xn_adv_bnds(n,i,k,2)
                enddo
             enddo
             do i=li1+1,limax
                do j=1,ljmax
                   xn_adv_bnde(n,j,k,1)=xn_adv_bnde(n,j,k,2)
                enddo
             enddo
             do j=lj1+1,ljmax
                do i=1,limax
                   xn_adv_bndn(n,i,k,1)=xn_adv_bndn(n,i,k,2)
                enddo
             enddo
          enddo
       endif
       if(me==0)then
          !Could fetch one level at a time if sizes becomes too big

          call check(nf90_inq_varid(ncid=ncFileID, name=trim(species(NSPEC_SHL+n)%name), varID=varID))

          call check(nf90_get_var(ncFileID, varID, data &
               ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext,1 /) ))

       endif
         CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 

       !overwrite Global Boundaries (lateral faces)

       do k=1,KMAX_ext
          do i=1,li0-1
             do j=1,ljmax
                xn_adv_bndw(n,j,k,2)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)

             enddo
          enddo
          do j=1,lj0-1
             do i=1,limax 
                xn_adv_bnds(n,i,k,2)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)
             enddo
          enddo
          do i=li1+1,limax
             do j=1,ljmax
                xn_adv_bnde(n,j,k,2)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)
             enddo
          enddo
          do j=lj1+1,ljmax
             do i=1,limax
                xn_adv_bndn(n,i,k,2)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)
             enddo
          enddo
       enddo
    enddo



    deallocate(data)
    if(me==0)then 
       call check(nf90_close(ncFileID))
    endif

  end subroutine read_newdata_LATERAL

  subroutine reset_3D(nseconds_indate)
    implicit none
    real, allocatable, dimension(:,:,:) ::data
    integer :: nseconds(1),ndate(4),n1,n,i,j,k,II,JJ,itime,status  
    integer :: nseconds_indate
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID

    allocate(data(GIMAX_ext,GJMAX_ext,KMAX_ext), stat=status)
    if(me==0)then 
       call check(nf90_open(path = trim(fileName_read), mode = nf90_nowrite, ncid = ncFileID))

       call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))
       do n=1,Next
          call check(nf90_get_var(ncFileID, varID, nseconds,start=(/ n /),count=(/ 1 /) ))
          if(nseconds(1)>=nseconds_indate)then
             write(*,*)'found date '
             call datefromsecondssince1970(ndate,nseconds(1),1)
             goto 876
          endif
       enddo
       write(*,*)'WARNING: did not find correct date'
       itime=Next
876    continue
       itime=n
    endif

    do n= 1, NSPEC_ADV
       if(me==0)then
          !Could fetch one level at a time if sizes becomes too big

          call check(nf90_inq_varid(ncid=ncFileID, name=trim(species(NSPEC_SHL+n)%name), varID=varID))

          call check(nf90_get_var(ncFileID, varID, data &
               ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext,1 /) ))

       endif
         CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 

       ! overwrite everything 3D (init)
       do k=1,KMAX_ext
          do j=1,ljmax
             do i=1,limax
                xn_adv(n,i,j,k)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)
             enddo
          enddo
       enddo

    enddo

    deallocate(data)
    if(me==0)then 
       call check(nf90_close(ncFileID))
    endif

  end subroutine reset_3D

end module Nest_ml

