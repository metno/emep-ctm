! <NetCDF_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module NetCDF_mod
!
! Routines for netCDF output
!
! Written by Peter january 2003->
!
!for details see:
!http://www.unidata.ucar.edu/software/netcdf/
!
!
!To improve: When output is onto the same file, but with different positions
!for the lower left corner, the coordinates i_EMEP j_EMEP and long lat will
!be wrong
!
use Chemfields_mod,     only : xn_shl,xn_adv
use CheckStop_mod,      only : CheckStop,StopAll,check=>CheckNC
use ChemDims_mod,       only : NSPEC_TOT, NSPEC_ADV, NSPEC_SHL
use ChemSpecs_mod,      only : species
use Config_module,       only: KMAX_MID,KMAX_BND, runlabel1, runlabel2&
                             ,MasterProc, NETCDF_DEFLATE_LEVEL &
                             ,NPROC, IIFULLDOM,JJFULLDOM &
                             ,IOU_INST,IOU_YEAR,IOU_MON,IOU_DAY &
                             ,IOU_HOUR,IOU_HOUR_INST &
                             ,PT,Pref,NLANDUSEMAX, model&
                             ,USE_EtaCOORDINATES,RUNDOMAIN&
                             ,fullrun_DOMAIN,month_DOMAIN,day_DOMAIN,hour_DOMAIN&
                             ,SurfacePressureFile &
                             , num_lev3d,lev3d&      ! 3D levels on 3D output
                             , startdate, step_main
use Country_mod,        only : NLAND, Country
use Debug_module,       only : DEBUG_NETCDF, DEBUG_NETCDF_RF
use Functions_mod,       only: StandardAtmos_km_2_kPa
use GridValues_mod,     only : GRIDWIDTH_M,fi,xp,yp,xp_EMEP_official&
                             ,debug_proc, debug_li, debug_lj &
                             ,yp_EMEP_official,fi_EMEP,GRIDWIDTH_M_EMEP&
                             ,grid_north_pole_latitude&
                             ,grid_north_pole_longitude,dx_rot,dx_roti,x1_rot,y1_rot&
                             !,GlobalPosition,glat_fdom,glon_fdom&
                             ,ref_latitude&
                             ,projection, sigma_mid,gb_stagg,gl_stagg,glon&
                             ,sigma_bnd&
                             ,glat,lb2ij,A_bnd,B_bnd&
                             ,lb2ij,lb2ijm,ij2lb&
                             ,ref_latitude_EMEP,xp_EMEP_old,yp_EMEP_old&
                             ,i_local,j_local,i_fdom,j_fdom&
                             ,Eta_bnd,Eta_mid,A_bnd,B_bnd,A_mid,B_mid&
                             ,lon0_lambert,&!reference longitude, also called phi, at which y=0 if lat=lat0_lambert
                             lat0_lambert,&!reference latitude, at which x=0 
                             lat_stand1_lambert,&!standard latitude at which mapping factor=1
                             lat_stand2_lambert,&!second standard latitude
                             x1_lambert,& !x value at i=1
                             y1_lambert  !y value at j=1                             
use InterpolationRoutines_mod,  only : grid2grid_coeff
use MPI_Groups_mod,     only: MPI_LOGICAL, MPI_SUM,MPI_INTEGER, MPI_BYTE,MPISTATUS, &
                               MPI_COMM_IO, MPI_COMM_CALC, IERROR, ME_IO, ME_CALC
use netcdf
use OwnDataTypes_mod,   only : Deriv, TXTLEN_NAME, TXTLEN_FILE
use Par_mod,            only : me,GIMAX,GJMAX,MAXLIMAX, MAXLJMAX, &
                              IRUNBEG,JRUNBEG,limax,ljmax, &
                              gi0,gj0,tgi0,tgi1,tgj0,tgj1,tlimax,tljmax
use PhysicalConstants_mod,  only : PI, EARTH_RADIUS
use TimeDate_mod,       only: nmdays,leapyear ,current_date, date,julian_date
use TimeDate_ExtraUtil_mod,only: date2nctime
use SmallUtils_mod,      only: wordsplit, find_index

implicit none

character(len=TXTLEN_FILE), save :: &
  fileName      = 'NotSet',&
  fileName_iou(IOU_INST:IOU_HOUR_INST)=&
    ['out_inst.nc    ','out_year.nc    ','out_month.nc   ','out_day.nc     ',&
     'out_hour.nc    ','out_hourInst.nc']
character(len=125) :: period_type !TESTHH

integer,parameter ::closedID=-999     !flag for showing that a file is closed
integer      :: ncFileID_new=closedID  !don't save because should always be
                !redefined (in case several routines are using ncFileID_new
                !with different filename_given)
integer,save :: ncFileID_iou(IOU_INST:IOU_HOUR_INST)=closedID
integer,save :: outCDFtag=0
!CDF types for output:
integer, public, parameter  :: Int1=1,Int2=2,Int4=3,Real4=4,Real8=5
character (len=18),parameter::Default_projection_name = 'General_Projection'

public :: Out_netCDF
public :: printCDF   ! minimal caller for Out_netCDF
public :: CloseNetCDF
public :: Init_new_netCDF
public :: GetCDF
public :: GetCDF_modelgrid
public :: WriteCDF !for testing purposes
public :: ReadField_CDF
public :: ReadField_CDF_FL
public :: ReadTimeCDF
public :: vertical_interpolate
public :: Out_CDF_sondes
public :: Create_CDF_sondes
public :: check
public :: IsCDFfractionFormat
public :: ReadSectorName
public :: check_lon_lat
public :: make_gridresolution
public :: create_country_emission_file
public :: output_country_emissions



private :: CreatenetCDFfile
private :: createnewvariable

contains
!_______________________________________________________________________

subroutine Out_CDF_sondes(fileName,SpecName,NSpec,Values,NLevels,g_ps,debug)
  !writes out a given variable with attributes in 1D (vertical coordinates) + station dimension (and time).
  !Note: defining a different variable for each station/component pair, would be very slow (almost double cpu time of run)
  !Note: writing only one variable at a time would be extremely slow

  character(len=*),  intent(in)  :: fileName,SpecName(NSpec)
  integer, intent(in)  :: NLevels,NSpec
  real,intent(in)  :: Values(Nlevels,NSpec,*),g_ps(*)
  logical, intent(in), optional  :: debug
  logical  :: debug_1D
  integer :: ncFileID
  integer :: varID,dimID
  character(len=8)  :: lastmodified_date
  character(len=10) :: lastmodified_hour
  integer :: i,iSpec,nrecords,nstations
  real :: rdays
  real,save,allocatable :: buff(:,:)


  debug_1D=DEBUG_NETCDF
  if(present(debug))debug_1D = debug .or. debug_1D
  
  if(MasterProc)then
    if(debug_1D)write(*,*)'writing in ',trim(fileName)
    !The file must have been created already with "Create_CDF_1D"
    call check(nf90_open(fileName,nf90_share+nf90_write,ncFileID))
  ! netCDF functions calls (optional=arg)
  ! nf90_open(path,cmode,ncid)
  ! nf90_inq_dimid(ncid,name,dimid)
  ! nf90_inquire_dimension(ncid,dimid,name=name,len=len)
  ! nf90_inq_varid(ncid,name,varid)
  ! nf90_put_var(ncid,varid,values,start=start,count=count)
    call Date_And_Time(date=lastmodified_date,time=lastmodified_hour)
    call check(nf90_inq_dimid(ncFileID,"time",dimID),"dim:time")
    call check(nf90_inquire_dimension(ncFileID,dimID,len=nrecords),"len:time")
    call check(nf90_inq_dimid(ncFileID,"station",dimID),"dim:station")
    call check(nf90_inquire_dimension(ncFileID,dimID,len=nstations),"len:station")
    if(debug_1D)write(*,*)'number of stations ',nstations
    nrecords=nrecords+1
    if(debug_1D)write(*,*)'writing on record ',nrecords
    call check(nf90_inq_varid(ncFileID,"time",varID),"inq:time")
    call date2nctime(current_date,rdays)
    call check(nf90_put_var(ncFileID,varID,rdays,start=[nrecords]),"put:time")
    if(.not.allocated(buff).and.NLevels>1)allocate(buff(nstations,Nlevels))
    do iSpec=1,NSpec
      !The variable must have been created already with "Create_CDF_1D"
      call check(nf90_inq_varid(ncFileID,SpecName(iSpec),varID),&
                "inq:"//trim(SpecName(iSpec)))
      if(debug_1D)write(*,"(A,1X,A20,2(1X,I0))")&
        'writing',trim(SpecName(iSpec)),NLevels,nstations
      if(NLevels==1)then
        call check(nf90_put_var(ncFileID,varID,values(1,iSpec,1:nstations),&
            start=[1,nrecords],count=[nstations,1]),"put:"//trim(SpecName(iSpec)))
      else
        do i=1,nstations     
           buff(i,1:NLevels)=Values(1:NLevels,iSpec,i)!NB: indices are switched (transposed of matrix)
        end do
        call check(nf90_put_var(ncFileID,varID,buff,&
            start=[1,1,nrecords],count=[nstations,NLevels,1]),&
            "put:"//trim(SpecName(iSpec)))
    !   call check(nf90_put_var(ncFileID,varID,buff,start=[1,1,nrecords],count=[1,1,1]))
      end if
      if(NLevels>1)then
        call check(nf90_inq_varid(ncFileID,'PS',varID),"inq:PS")
        call check(nf90_put_var(ncFileID,varID,g_ps(1:nstations),&
            start=[1,nrecords],count=[nstations,1]),"put:PS")
      end if
    end do
    call check(nf90_put_att(ncFileID,nf90_global,"lastmodified_date",lastmodified_date))
    call check(nf90_put_att(ncFileID,nf90_global,"lastmodified_hour",lastmodified_hour))
    call check(nf90_close(ncFileID))
  end if !Masterproc
end subroutine Out_CDF_sondes

subroutine Create_CDF_sondes(fileName,NSpec,NSpec_Att,SpecDef,&
        NStations,NMetaData,MetaData,KMAXcdf,CDFtype,debug)
! Create a NetCDF file with given attributes.
! if Nvalues>1 define vertical coordinates.
  integer, intent(in)  :: KMAXcdf,NSpec,NStations,NSpec_Att,NMetaData
  character(len=*),  intent(in) :: fileName, &
    SpecDef(NSpec,0:NSpec_Att),MetaData(0:NStations,NMetaData)
  integer, intent(in), optional :: CDFtype
  logical, intent(in), optional :: debug
  logical  :: debug_1D
  integer :: OUTtype,ncFileID,levDimID,ilevDimID,timeDimID
  integer :: varID,StationDimID,StringDimID
  character (len=*), parameter :: vert_coord='atmosphere_hybrid_sigma_pressure_coordinate'
  character(len=8)  :: lastmodified_date
  character(len=10) :: lastmodified_hour
  integer :: i,k,n,iSpec,iSta
  real :: kcoord(KMAXcdf+1)
  real :: Acdf(KMAXcdf),Bcdf(KMAXcdf),Aicdf(KMAXcdf+1),Bicdf(KMAXcdf+1)
  integer,parameter :: MAX_String_length=36
  character(len=100) :: auxL(4)
  character(len=MAX_String_length) :: metaName,metaType,auxC(NStations)
  integer :: auxI(NStations),ierr
  real :: auxR(NStations)

  debug_1D=DEBUG_NETCDF
  if(present(debug))debug_1D = debug .or. debug_1D
  
  if(MasterProc)then
    !Create file if it does not yet exist
    if(debug_1D)write(*,*)'creating ',trim(fileName)
    call check(nf90_create(fileName,nf90_clobber,ncFileID),"create:"//trim(fileName))
  ! netCDF functions calls (wo optional arg.)
  ! nf90_create(path,cmode,ncid)
  ! nf90_def_dim(ncid,name,len,dimid)
  ! nf90_def_var(ncid,name,xtype,dimids,varid)
  ! nf90_put_att(ncid,varid,name,values)
  ! nf90_inq_varid(ncid,name,varid)
  ! nf90_put_var(ncid,varid,values)
    call check(nf90_def_dim(ncFileID,"time",nf90_unlimited,timeDimID),"dim:time")
    call check(nf90_def_dim(ncFileID,"station",NStations,StationDimID),"dim:station")
    call check(nf90_def_dim(ncFileID,"string_length",MAX_String_length,StringDimID),"dim:slen")

    if(KMAXcdf>1)then
      if(debug_1D)write(*,*)'defining vertical levels ',KMAXcdf
      call check(nf90_put_att(ncFileID,nf90_global,"vert_coord",vert_coord))
      call check(nf90_def_dim(ncFileID, "lev",KMAXcdf  , levDimID),"dim:lev")
      call check(nf90_def_dim(ncFileID,"ilev",KMAXcdf+1,ilevDimID),"dim:ilev")        
      call check(nf90_def_var(ncFileID,"lev",nf90_double,levDimID,varID))
      call check(nf90_put_att(ncFileID,varID,"standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
       call check(nf90_put_att(ncFileID,varID,"long_name","hybrid level at layer midpoints (A/P0+B)"))
      call check(nf90_put_att(ncFileID,varID,"positive","up"))
      call check(nf90_put_att(ncFileID,varID,"formula_terms","ap: hyam b: hybm ps: PS p0: P0"))
        !p(n,k,j,i) = a(k)+ b(k)*ps(n,j,i)
      call check(nf90_def_var(ncFileID,"ilev",nf90_double,ilevDimID,varID))
      call check(nf90_put_att(ncFileID,varID,"standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
      call check(nf90_put_att(ncFileID,varID,"long_name","hybrid level at layer interfaces (A/P0+B)"))
      call check(nf90_put_att(ncFileID,varID,"positive","up"))
      call check(nf90_put_att(ncFileID,varID,"formula_terms","ap: hyai b: hybi ps: PS p0: P0"))

      call check(nf90_def_var(ncFileID,"P0",nf90_double,varID))
      call check(nf90_put_att(ncFileID,varID,"units","hPa"))
      call check(nf90_def_var(ncFileID,"PS",nf90_double,[StationDimID,timeDimID],varID))
      call check(nf90_put_att(ncFileID,varID,"long_name","Surface pressure"))
      call check(nf90_put_att(ncFileID,varID,"units","hPa"))
        
      !The hybrid sigma-pressure coordinate for level k is defined as ap(k)/p0+b(k). 
      call check(nf90_def_var(ncFileID,"hyam",nf90_double,levDimID,varID))
      call check(nf90_put_att(ncFileID,varID,"long_name","hybrid A coefficient at layer midpoints"))
      call check(nf90_put_att(ncFileID,varID,"units","hPa"))
      call check(nf90_def_var(ncFileID,"hybm",nf90_double,levDimID,varID))
      call check(nf90_put_att(ncFileID,varID,"long_name","hybrid B coefficient at layer midpoints"))
        
      call check(nf90_def_var(ncFileID,"hyai",nf90_double,ilevDimID,varID))
      call check(nf90_put_att(ncFileID,varID,"long_name","hybrid A coefficient at layer interfaces"))
      call check(nf90_put_att(ncFileID,varID,"units","hPa"))
      call check(nf90_def_var(ncFileID,"hybi",nf90_double,ilevDimID,varID) )
      call check(nf90_put_att(ncFileID,varID,"long_name","hybrid B coefficient at layer interfaces"))       
    else
      if(debug_1D)write(*,*)'not defining vertical levels ',KMAXcdf
    end if

    call check(nf90_def_var(ncFileID,"time",nf90_double,timeDimID,varID))
    call check(nf90_put_att(ncFileID,varID,"long_name", "time (instantaneous)"))
    call check(nf90_put_att(ncFileID,varID,"units","days since 1900-1-1 0:0:0"))

    call Date_And_Time(date=lastmodified_date,time=lastmodified_hour)
    call check(nf90_put_att(ncFileID,nf90_global,"created_date",lastmodified_date))
    call check(nf90_put_att(ncFileID,nf90_global,"created_hour",lastmodified_hour))
         
    !global attributes
    do n=1,NMetaData
      if(MetaData(0,n)=="")cycle
      if(debug_1D)write(*,*)'Global Attribute ',trim(MetaData(0,n))
      call wordsplit(trim(MetaData(0,n)),3,auxL,k,ierr,strict_separator=':')
      call CheckStop(3,k,&
        "NetCDF_mod: too short metadata definition "//trim(MetaData(0,n)))
      select case(auxL(2))
      case("c","C","s","S") ! string/char attribute
        call check(nf90_put_att(ncFileID,nf90_global,trim(auxL(1)),trim(auxL(3))),&
                   "MetaData="//trim(MetaData(0,n)))
      case("i","I","n","N") ! integer attribute
        read(auxL(3),*)auxI(1)
        call check(nf90_put_att(ncFileID,nf90_global,trim(auxL(1)),auxI(1)),&
                   "MetaData="//trim(MetaData(0,n)))
      case("f","F","d","D") ! float/double attribute
        read(auxL(3),*)auxR(1)
        call check(nf90_put_att(ncFileID,nf90_global,trim(auxL(1)),auxR(1)),&
                   "MetaData="//trim(MetaData(0,n)))
      case default
        call CheckStop("NetCDF_mod: unknown metadata-type "//trim(MetaData(0,n)))
      end select
    end do
     
    OUTtype=Real4  !default value
    if(present(CDFtype))OUTtype=CDFtype
    select case(OUTtype)
      case(Int1 );OUTtype=nf90_byte
      case(Int2 );OUTtype=nf90_short
      case(Int4 );OUTtype=nf90_int
      case(Real4);OUTtype=nf90_float
      case(Real8);OUTtype=nf90_double
      case default;
        call CheckStop("NetCDF_mod:undefined datatype")
    end select

    do iSpec=1,NSpec
      !define the variables
      if(debug_1D)write(*,*)'Defining new variable ',trim(SpecDef(iSpec,0)),OUTtype
      if(KMAXcdf>1)then
        call check(nf90_def_var(ncFileID,SpecDef(iSpec,0),OUTtype,&
              [StationDimID,levDimID,timeDimID],varID))
      else
        call check(nf90_def_var(ncFileID,SpecDef(iSpec,0),OUTtype,&
              [StationDimID,timeDimID],varID))
      end if
      !species attributes:
      do n=1,NSpec_Att
        if(SpecDef(iSpec,n)=="")cycle
        call wordsplit(trim(SpecDef(iSpec,n)),3,auxL,k,ierr,strict_separator=':')
        call CheckStop(3,k,&
          "NetCDF_mod: too short metadata definition "//trim(SpecDef(iSpec,n)))
        select case(auxL(1))
        case("_FillValue")      ! ensure that OUTtype/_FillValue types match
          select case(OUTtype)  ! even under -r8 compilation option
          case(NF90_BYTE)   ! _FillValue=achar(0)
            call check(nf90_put_att(ncFileID,varID,trim(auxL(1)),NF90_FILL_BYTE))
          case(NF90_SHORT)  ! _FillValue=-32767
            call check(nf90_put_att(ncFileID,varID,trim(auxL(1)),NF90_FILL_SHORT))
          case(NF90_INT)    ! _FillValue=-2147483647
            call check(nf90_put_att(ncFileID,varID,trim(auxL(1)),NF90_FILL_INT))
          case(NF90_FLOAT)  ! _FillValue=9.9692099683868690e+36
            call check(nf90_put_att(ncFileID,varID,trim(auxL(1)),NF90_FILL_FLOAT))
          case(NF90_DOUBLE) ! _FillValue=9.9692099683868690e+36
            call check(nf90_put_att(ncFileID,varID,trim(auxL(1)),NF90_FILL_DOUBLE))
          end select
        case default
          auxC(1)=NF90_FILL_CHAR
          auxI(1)=NF90_FILL_INT
          auxR(1)=NF90_FILL_DOUBLE
          select case(auxL(2))
          case("c","C","s","S") ! string/char attribute
            if(auxL(3)/="missing")auxC(1)=trim(auxL(3))
            call check(nf90_put_att(ncFileID,varID,trim(auxL(1)),trim(auxC(1))),&
                       "SpecDef="//trim(SpecDef(iSpec,n)))
          case("i","I","n","N") ! integer attribute
            if(auxL(3)/="missing")read(auxL(3),*)auxI(1)
            call check(nf90_put_att(ncFileID,varID,trim(auxL(1)),auxI(1)),&
                       "SpecDef="//trim(SpecDef(iSpec,n)))
          case("f","F","d","D") ! float/double attribute
            if(auxL(3)/="missing")read(auxL(3),*)auxR(1)
            call check(nf90_put_att(ncFileID,varID,trim(auxL(1)),auxR(1)),&
                       "SpecDef="//trim(SpecDef(iSpec,n)))
          case default
            call CheckStop("NetCDF_mod: unknown metadata-type "//trim(SpecDef(iSpec,n)))
          end select
        end select
      end do
    end do

    if(debug_1D)write(*,*)'Defining station meta-data'
    do n=1,NMetaData
      if(MetaData(1,n)=="")cycle
      if(debug_1D)write(*,*)'Defining station ',trim(MetaData(1,n))
      call wordsplit(trim(MetaData(1,n)),3,auxL,k,ierr,strict_separator=':')
      call CheckStop(3,k,&
        "NetCDF_mod: too short metadata definition "//trim(MetaData(1,n)))
      select case(auxL(2))
      case("c","C","s","S") ! string/char attribute
        call check(nf90_def_var(ncFileID,auxL(1),nf90_char,&
            [StringDimID,StationDimID],varID),"def:"//trim(auxL(1)))
        call check(nf90_put_att(ncFileID,varID,"_FillValue",NF90_FILL_CHAR),&
            "def:"//trim(auxL(1))//"@_FillValue")
      case("i","I","n","N") ! integer attribute
        call check(nf90_def_var(ncFileID,auxL(1),nf90_int,&
            StationDimID,varID),"def:"//trim(auxL(1)))
       call check(nf90_put_att(ncFileID,varID,"_FillValue",NF90_FILL_INT),&
            "def:"//trim(auxL(1))//"@_FillValue")
      case("f","F","d","D") ! float/double attribute
        call check(nf90_def_var(ncFileID,auxL(1),nf90_double,&
            StationDimID,varID),"def:"//trim(auxL(1)))
        call check(nf90_put_att(ncFileID,varID,"_FillValue",NF90_FILL_DOUBLE),&
            "def:"//trim(auxL(1))//"@_FillValue")
      case default
        call CheckStop("NetCDF_mod: unknown metadata-type "//trim(MetaData(1,n)))
      end select
    end do

    call check(nf90_put_att(ncFileID,nf90_global,"lastmodified_date",lastmodified_date))
    call check(nf90_put_att(ncFileID,nf90_global,"lastmodified_hour",lastmodified_hour))
    call check(nf90_enddef(ncFileID),"define done:"//trim(fileName))

    if(debug_1D)write(*,*)'Writing station meta-data'
    do n=1,NMetaData
      if(MetaData(1,n)=="")cycle
      if(debug_1D)write(*,*)'Writing station ',trim(MetaData(1,n))
      call wordsplit(trim(MetaData(1,n)),3,auxL,k,ierr,strict_separator=':')
      call CheckStop(3,k,&
        "NetCDF_mod: too short metadata definition "//trim(MetaData(1,n)))
      metaName=trim(auxL(1))
      metaType=trim(auxL(2))
      auxC(:)=NF90_FILL_CHAR
      auxI(:)=NF90_FILL_INT
      auxR(:)=NF90_FILL_DOUBLE
      do iSta=1,NStations
        if(MetaData(iSta,n)=="")cycle
        call wordsplit(trim(MetaData(iSta,n)),3,auxL,k,ierr,strict_separator=':')
        call CheckStop(3,k,&
          "NetCDF_mod: too short metadata definition "//trim(MetaData(iSta,n)))
        call CheckStop(metaName,auxL(1),&
          "NetCDF_mod: inconsistent metadata-name definition "//trim(MetaData(iSta,n)))
        call CheckStop(metaType,auxL(2),&
          "NetCDF_mod: inconsistent metadata-type definition "//trim(MetaData(iSta,n)))
        if(auxL(3)=="missing")cycle
        select case(metaType)
          case("c","C","s","S");auxC(iSta)=trim(auxL(3))  ! string/char attribute
          case("i","I","n","N");read(auxL(3),*)auxI(iSta) ! integer attribute
          case("f","F","d","D");read(auxL(3),*)auxR(iSta) ! float/double attribute
        end select
      end do
      call check(nf90_inq_varid(ncFileID,metaName,varID),"inq:"//trim(metaName))
      select case(metaType)
      case("c","C","s","S") ! string/char attribute
        call check(nf90_put_var(ncFileID,varID,auxC),"put:"//trim(metaType))
      case("i","I","n","N") ! integer attribute
        call check(nf90_put_var(ncFileID,varID,auxI),"put:"//trim(metaType))
      case("f","F","d","D") ! float/double attribute
        call check(nf90_put_var(ncFileID,varID,auxR),"put:"//trim(metaType))
      case default
        call CheckStop("NetCDF_mod: wrong metadata-type definition "//trim(metaType))
      end select
    end do

    if(KMAXcdf>1)then
      call check(nf90_inq_varid(ncFileID,"P0",varID),"inq:P0")
      call check(nf90_put_var(ncFileID,varID,Pref/100.0 ),"put:P0")

      do k=1,KMAXcdf
       !REVERSE order of k !
        Acdf(k)=A_mid(KMAX_MID-k+1)
        Bcdf(k)=B_mid(KMAX_MID-k+1)
        Aicdf(k)=A_bnd(KMAX_BND-k+1)
        Bicdf(k)=B_bnd(KMAX_BND-k+1)
        if(k==KMAXcdf)then
          Aicdf(k+1)=A_bnd(KMAX_BND-k)
          Bicdf(k+1)=B_bnd(KMAX_BND-k)
        end if
      end do
        
      call check(nf90_inq_varid(ncFileID,"hyam",varID),"inq:hyam")
      call check(nf90_put_var(ncFileID,varID,Acdf(1:KMAXcdf)/1e2),"put:hyam")
      call check(nf90_inq_varid(ncFileID,"hybm",varID),"inq:hybm")
      call check(nf90_put_var(ncFileID,varID,Bcdf(1:KMAXcdf)),"put:hybm")
      call check(nf90_inq_varid(ncFileID,"hyai",varID),"inq:hyai")
      call check(nf90_put_var(ncFileID,varID,Aicdf(1:KMAXcdf+1)/1e2),"put:hyai")
      call check(nf90_inq_varid(ncFileID,"hybi",varID),"inq:hybi")
      call check(nf90_put_var(ncFileID,varID,Bicdf(1:KMAXcdf+1)),"put:hybi")

      do i=1,KMAXcdf
        kcoord(i)=Acdf(i)/Pref+Bcdf(i)
      end do
      call check(nf90_inq_varid(ncFileID,"lev",varID),"inq:lev")
      call check(nf90_put_var(ncFileID,varID,kcoord(1:KMAXcdf))   ,"put:lev")

      do i=1,KMAXcdf+1
        kcoord(i)=Aicdf(i)/Pref+Bicdf(i)
      end do
      call check(nf90_inq_varid(ncFileID,"ilev",varID),"inq:ilev")
      call check(nf90_put_var(ncFileID,varID,kcoord(1:KMAXcdf+1)),"put:ilev")
    end if
    call check(nf90_close(ncFileID))    
  end if !Masterproc
end subroutine Create_CDF_sondes

subroutine Init_new_netCDF(fileName,iotyp)
integer,  intent(in) :: iotyp
character(len=*),  intent(in)  :: fileName

integer :: GIMAXcdf,GJMAXcdf,IBEGcdf,JBEGcdf,KMAXcdf
integer :: i1,i2,j1,j2

call CloseNetCDF !must be called by all procs, to syncronize outCDFtag

!if(MasterProc)then
if(MasterProc.and.DEBUG_NETCDF )&
  write(*,*)'Init_new_netCDF ',trim(fileName),iotyp

!NB IBEGcdf and JBEGcdf are here defined relative to fulldomain
IBEGcdf=GIMAX+IRUNBEG-1; JBEGcdf=GJMAX+JRUNBEG-1  !initialisations
GIMAXcdf=0; GJMAXcdf=0                            !initialisations
KMAXcdf=1
select case (iotyp)
case (IOU_YEAR)
  fileName_iou(iotyp) = trim(fileName)
  period_type = 'fullrun'
  i1=fullrun_DOMAIN(1);i2=fullrun_DOMAIN(2);j1=fullrun_DOMAIN(3);j2=fullrun_DOMAIN(4)
  IBEGcdf=min(IBEGcdf,i1); JBEGcdf=min(JBEGcdf,j1)
  GIMAXcdf=max(GIMAXcdf,i2-i1+1); GJMAXcdf=max(GJMAXcdf,j2-j1+1)
case(IOU_MON)
  fileName_iou(iotyp) = trim(fileName)
  period_type = 'monthly'
  i1=month_DOMAIN(1);i2=month_DOMAIN(2);j1=month_DOMAIN(3);j2=month_DOMAIN(4)
  IBEGcdf=min(IBEGcdf,i1); JBEGcdf=min(JBEGcdf,j1)
  GIMAXcdf=max(GIMAXcdf,i2-i1+1); GJMAXcdf=max(GJMAXcdf,j2-j1+1)
case(IOU_DAY)
  fileName_iou(iotyp) = trim(fileName)
  period_type = 'daily'
  i1=day_DOMAIN(1);i2=day_DOMAIN(2);j1=day_DOMAIN(3);j2=day_DOMAIN(4)
  IBEGcdf=min(IBEGcdf,i1); JBEGcdf=min(JBEGcdf,j1)
  GIMAXcdf=max(GIMAXcdf,i2-i1+1); GJMAXcdf=max(GJMAXcdf,j2-j1+1)
case(IOU_HOUR)
  fileName_iou(iotyp) = trim(fileName)
  period_type = 'hourlyMean'
  i1=hour_DOMAIN(1);i2=hour_DOMAIN(2);j1=hour_DOMAIN(3);j2=hour_DOMAIN(4)
  IBEGcdf=min(IBEGcdf,i1); JBEGcdf=min(JBEGcdf,j1)
  GIMAXcdf=max(GIMAXcdf,i2-i1+1); GJMAXcdf=max(GJMAXcdf,j2-j1+1)
case(IOU_HOUR_INST)
  fileName_iou(iotyp) = trim(fileName)
  period_type = 'hourly'
  i1=hour_DOMAIN(1);i2=hour_DOMAIN(2);j1=hour_DOMAIN(3);j2=hour_DOMAIN(4)
  IBEGcdf=min(IBEGcdf,i1); JBEGcdf=min(JBEGcdf,j1)
  GIMAXcdf=max(GIMAXcdf,i2-i1+1); GJMAXcdf=max(GJMAXcdf,j2-j1+1)
case(IOU_INST)
  fileName_iou(iotyp) = trim(fileName)
  period_type = 'instant'
case default
  period_type = 'unknown'
end select

if(MasterProc.and.DEBUG_NETCDF)&
  write(*,*) "Creating ", trim(fileName),' ',trim(period_type)
call CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,IBEGcdf,JBEGcdf,&
                      num_lev3d,KLEVcdf=lev3d,KLEVcdf_from_top=.true.)

if(MasterProc.and.DEBUG_NETCDF)&
  write(*,*) "Finished Init_new_netCDF", trim(fileName),' ',trim(period_type)
!end if
end subroutine Init_new_netCDF

subroutine CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,IBEGcdf,JBEGcdf,&
     KMAXcdf,KLEVcdf,KLEVcdf_from_top,RequiredProjection)

!IBEGcdf,JBEGcdf relative to fulldomain
  integer, intent(in) :: GIMAXcdf,GJMAXcdf,IBEGcdf,JBEGcdf,KMAXcdf
  character(len=*),  intent(in)  :: fileName
  integer, intent(in), optional :: KLEVcdf(KMAXcdf)
  logical, intent(in), optional :: KLEVcdf_from_top
  character(len=*),optional, intent(in):: RequiredProjection

  character(len=*), parameter :: author_of_run='emepctm group'
  character(len=19) :: projection_params='90.0 -32.0 0.933013' !set later on

  real :: xcoord(GIMAX),ycoord(GJMAX)

  character(len=8)  :: created_date,lastmodified_date
  character(len=10) :: created_hour,lastmodified_hour
  integer :: iDimID,jDimID,levDimID,ilevDimID,timeDimID,VarID,iVarID,jVarID
  integer :: ncFileID,iEMEPVarID,jEMEPVarID,latVarID,longVarID
  integer :: i1,j1,istart,jstart,ishift,jshift,icount,jcount
  integer :: i,j,k,iproc
  logical :: levels_from_top=.true.
  real :: scale_at_projection_origin
  character(len=80) ::UsedProjection
  character (len=*), parameter :: vert_coord='atmosphere_hybrid_sigma_pressure_coordinate'

  real ::Buff2D(MAXLIMAX,MAXLJMAX,2)

  ! fileName: Name of the new created file
  ! nf90_clobber: protect existing datasets
  ! ncFileID: netcdf ID

  !Check that the dimensions are > 0
  if(GIMAXcdf<=0.or.GJMAXcdf<=0.or.KMAXcdf<=0)then
     write(*,*)'WARNING:'
     write(*,*)trim(fileName),&
          ' not created. Requested area too small (or outside domain) '
     write(*,*)'sizes (IMAX,JMAX,IBEG,JBEG,KMAX) ',&
          GIMAXcdf,GJMAXcdf,IBEGcdf,JBEGcdf,KMAXcdf
     return
  end if

  if(present(RequiredProjection))then
    UsedProjection=trim(RequiredProjection)
  else
    UsedProjection=trim(projection)
  end if

  if(MasterProc)write(*,*)'creating ',trim(fileName)
  if(MasterProc.and.DEBUG_NETCDF)write(*,*)'UsedProjection ',trim(UsedProjection)
  if(MasterProc.and.DEBUG_NETCDF)write(*,fmt='(A,8I7)')'with sizes (IMAX,JMAX,IBEG,JBEG,KMAX) ',&
       GIMAXcdf,GJMAXcdf,IBEGcdf,JBEGcdf,KMAXcdf

  if(MasterProc)then
    if(NETCDF_DEFLATE_LEVEL >= 0)then
      if(MasterProc.and.DEBUG_NETCDF)write(*,*)'nf90_create'
      call check(nf90_create(fileName,nf90_hdf5,ncFileID),"create:"//trim(fileName))
      if(MasterProc.and.DEBUG_NETCDF)write(*,*)'nf90_created'
    else
      call check(nf90_create(fileName,nf90_clobber,ncFileID),"create:"//trim(fileName))
    end if

   ! Define the dimensions
    select case(UsedProjection)
    case('Stereographic')
      ! define coordinate dimensions
      call check(nf90_def_dim(ncFileID,"i",GIMAXcdf,iDimID),"dim:i")
      call check(nf90_def_dim(ncFileID,"j",GJMAXcdf,jDimID),"dim:j")

      scale_at_projection_origin=(1.+sin(ref_latitude*PI/180.))/2.
      write(projection_params,fmt='(''90.0 '',F5.1,F9.6)')fi,scale_at_projection_origin
      call check(nf90_put_att(ncFileID,nf90_global,"projection_params",projection_params))

      ! define coordinate variables
      iVarID    =define_var("i"     ,nf90_float,[iDimID])
      jVarID    =define_var("j"     ,nf90_float,[jDimID])
      iEMEPVarID=define_var("i_EMEP",nf90_float,[iDimID])
      jEMEPVarID=define_var("j_EMEP",nf90_float,[jDimID])
      latVarID  =define_var("lat"   ,nf90_float,[iDimID,jDimID])
      longVarID =define_var("lon"   ,nf90_float,[iDimID,jDimID])
    case('lon lat')
      ! define coordinate dimensions
      call check(nf90_def_dim(ncFileID,"lon",GIMAXcdf,iDimID),"dim:lon")
      call check(nf90_def_dim(ncFileID,"lat",GJMAXcdf,jDimID),"dim:lat")

      ! define coordinate variables
      iVarID=define_var("lon",nf90_double,[iDimID])
      jVarID=define_var("lat",nf90_double,[jDimID])
    case('Rotated_Spherical')
      ! define coordinate dimensions
      call check(nf90_def_dim(ncFileID,"i",GIMAXcdf,iDimID),"dim:i_RotS")
      call check(nf90_def_dim(ncFileID,"j",GJMAXcdf,jDimID),"dim:j_RotS")

      ! define coordinate variables
      iVarID    =define_var("i_RotS",nf90_float,[iDimID])
      jVarID    =define_var("j_RotS",nf90_float,[jDimID])
      latVarID  =define_var("lat"   ,nf90_float,[iDimID,jDimID])
      longVarID =define_var("lon"   ,nf90_float,[iDimID,jDimID])
    case('lambert')
      ! define coordinate dimensions
      call check(nf90_def_dim(ncFileID,"i",GIMAXcdf,iDimID))
      call check(nf90_def_dim(ncFileID,"j",GJMAXcdf,jDimID))

      ! define coordinate variables
      iVarID    =define_var("i_lambert",nf90_float,[iDimID])
      jVarID    =define_var("j_lambert",nf90_float,[jDimID])
      latVarID  =define_var("lat"   ,nf90_float,[iDimID,jDimID])
      longVarID =define_var("lon"   ,nf90_float,[iDimID,jDimID])
    case default !general projection
      ! define coordinate dimensions
      call check(nf90_def_dim(ncFileID,"i",GIMAXcdf,iDimID),"dim:i")
      call check(nf90_def_dim(ncFileID,"j",GJMAXcdf,jDimID),"dim:j")

      ! define coordinate variables
      iVarID    =define_var("i"     ,nf90_float,[iDimID])
      jVarID    =define_var("j"     ,nf90_float,[jDimID])
      latVarID  =define_var("lat"   ,nf90_float,[iDimID,jDimID])
      longVarID =define_var("lon"   ,nf90_float,[iDimID,jDimID])
    end select
     
    if(MasterProc.and.DEBUG_NETCDF)write(*,*)'lon lat dims defined'
    call check(nf90_def_dim(ncFileID, "lev",KMAXcdf  , levDimID))
    call check(nf90_def_dim(ncFileID,"ilev",KMAXcdf+1,ilevDimID))
    call check(nf90_put_att(ncFileID, nf90_global,"vert_coord",vert_coord))

    call check(nf90_def_dim(ncFileID,"time",nf90_unlimited,timeDimID))

    call Date_And_Time(date=created_date,time=created_hour)
    if(DEBUG_NETCDF)write(*,"(2A)")&
      'created_date: ',created_date,'created_hour: ',created_hour

    ! Write global attributes
    call check(nf90_put_att(ncFileID,nf90_global,"Conventions", "CF-1.6" ))
   !call check(nf90_put_att(ncFileID,nf90_global,"version", version ))
    call check(nf90_put_att(ncFileID,nf90_global,"model", model))
    call check(nf90_put_att(ncFileID,nf90_global,"author_of_run", author_of_run))
    call check(nf90_put_att(ncFileID,nf90_global,"created_date", created_date))
    call check(nf90_put_att(ncFileID,nf90_global,"created_hour", created_hour))
    lastmodified_date = created_date
    lastmodified_hour = created_hour
    call check(nf90_put_att(ncFileID,nf90_global,"lastmodified_date", lastmodified_date))
    call check(nf90_put_att(ncFileID,nf90_global,"lastmodified_hour", lastmodified_hour))

    call check(nf90_put_att(ncFileID,nf90_global,"projection",UsedProjection))
    call check(nf90_put_att(ncFileID,nf90_global,"period_type", &
         trim(period_type)), "period_type:"//trim(period_type) )
    call check(nf90_put_att(ncFileID,nf90_global,"run_label", trim(runlabel2)))

    ! The hybrid sigma-pressure coordinate for level k is defined as ap(k)/p0+b(k). 
    ! p(n,k,j,i) = ap(k)+ b(k)*ps(n,j,i)
    varID=define_var("lev" ,nf90_double,[levDimID])
    varID=define_var("P0"  ,nf90_double,[0])
    varID=define_var("hyam",nf90_double,[levDimID])
    varID=define_var("hybm",nf90_double,[levDimID])

    varID=define_var("ilev",nf90_double,[ilevDimID])
    varID=define_var("hyai",nf90_double,[ilevDimID])
    varID=define_var("hybi",nf90_double,[ilevDimID])

    varID=define_var("time",nf90_double,[timeDimID])

    !CF-1.0 definitions:
    select case(UsedProjection)
    case('Stereographic')
      varID=define_var("Polar_Stereographic",nf90_int,[0])
    case('lon lat')
    ! no additional projection meta-data needed
    case('Rotated_Spherical')
      varID=define_var("Rotated_Spherical",nf90_int ,[0])
      varID=define_var("rotated_pole"     ,nf90_char,[0]) ! for NCL
    case('lambert')
      varID=define_var("projection_lambert",nf90_int ,[0])
    case default
      varID=define_var("Default_projection_name",nf90_int,[0])
    end select

    ! Leave define mode
    call check(nf90_enddef(ncFileID), "define_done"//trim(fileName) )

    ! Define horizontal distances
    select case(UsedProjection)
    case('Stereographic')
      xcoord(1)=(IBEGcdf-xp)*GRIDWIDTH_M/1000.
      do i=2,GIMAXcdf
        xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
      end do
      ycoord(1)=(JBEGcdf-yp)*GRIDWIDTH_M/1000.
      do j=2,GJMAXcdf
        ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
      end do
      call check(nf90_put_var(ncFileID,iVarID,xcoord(1:GIMAXcdf)))
      call check(nf90_put_var(ncFileID,jVarID,ycoord(1:GJMAXcdf)))

      ! Define horizontal coordinates in the official EMEP grid
      !xp_EMEP_official=8.
      !yp_EMEP_official=110.
      !GRIDWIDTH_M_EMEP=50000.
      !fi_EMEP=-32.
      if(fi==fi_EMEP)then
        ! Implemented only if fi = fi_EMEP = -32 (Otherwise needs a 2-dimensional mapping)
        ! uses (i-xp)*GRIDWIDTH_M = (i_EMEP-xp_EMEP)*GRIDWIDTH_M_EMEP
        do i=1,GIMAXcdf
          xcoord(i)=(i+IBEGcdf-1-xp)*GRIDWIDTH_M/GRIDWIDTH_M_EMEP + xp_EMEP_official
          !print *, i,xcoord(i)
        end do
        do j=1,GJMAXcdf
          ycoord(j)=(j+JBEGcdf-1-yp)*GRIDWIDTH_M/GRIDWIDTH_M_EMEP + yp_EMEP_official
          !print *, j,ycoord(j)
        end do
      else
        xcoord(1:GIMAXcdf)=NF90_FILL_FLOAT
        ycoord(1:GJMAXcdf)=NF90_FILL_FLOAT
      end if
      call check(nf90_put_var(ncFileID,iEMEPVarID,xcoord(1:GIMAXcdf)))
      call check(nf90_put_var(ncFileID,jEMEPVarID,ycoord(1:GJMAXcdf)))

    case('Rotated_Spherical')
      do i=1,GIMAXcdf
        xcoord(i)= (i+IBEGcdf-2)*dx_rot+x1_rot
      end do
      do j=1,GJMAXcdf
        ycoord(j)= (j+JBEGcdf-2)*dx_rot+y1_rot
      end do
      call check(nf90_put_var(ncFileID,iVarID,xcoord(1:GIMAXcdf)))
      call check(nf90_put_var(ncFileID,jVarID,ycoord(1:GJMAXcdf)))

    case('lambert')
      do i=1,GIMAXcdf
        xcoord(i)= (i+IBEGcdf-2)*GRIDWIDTH_M+x1_lambert
      end do
      do j=1,GJMAXcdf
        ycoord(j)= (j+JBEGcdf-2)*GRIDWIDTH_M+y1_lambert
      end do
      call check(nf90_put_var(ncFileID,iVarID,xcoord(1:GIMAXcdf)))
      call check(nf90_put_var(ncFileID,jVarID,ycoord(1:GJMAXcdf)))
    case('lon lat')
      do i=1,GIMAXcdf
        xcoord(i)= glon(1,1)+(i_local(i+IBEGcdf-1)-1)*(glon(2,1)-glon(1,1))
        !force monotone values:
        if(i>1)then
          !must first check that i>1 before testing xcoord(i-1) (to avoid debug errors)
          if(xcoord(i)<xcoord(i-1).and.xcoord(i)<0)xcoord(i)=xcoord(i)+360.0
        end if
      end do
      do j=1,GJMAXcdf
        ycoord(j)= glat(1,1)+(j_local(j+JBEGcdf-1)-1)*(glat(1,2)-glat(1,1))
      end do
      call check(nf90_put_var(ncFileID,iVarID,xcoord(1:GIMAXcdf)))
      call check(nf90_put_var(ncFileID,jVarID,ycoord(1:GJMAXcdf)))
    case default
      xcoord(1)=(IBEGcdf-0.5)*GRIDWIDTH_M/1000.
      do i=2,GIMAXcdf
        xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
       !print *, i,xcoord(i)
      end do
      ycoord(1)=(JBEGcdf-0.5)*GRIDWIDTH_M/1000.
      do j=2,GJMAXcdf
        ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
      end do
      call check(nf90_put_var(ncFileID,iVarID,xcoord(1:GIMAXcdf)))
      call check(nf90_put_var(ncFileID,jVarID,ycoord(1:GJMAXcdf)))
      !        write(*,*)'coord written'
    end select! projections

    ! Define vertical levels
    if(present(KLEVcdf))then     !order is defined in KLEVcdf
      levels_from_top=.false.
      if(present(KLEVcdf_from_top))levels_from_top=KLEVcdf_from_top
      call write_klev(KMAXcdf,KLEVcdf,levels_from_top)
    else
      levels_from_top=(KMAXcdf==KMAX_MID)
      call write_klev(KMAXcdf,[(k,k=1,KMAXcdf)],levels_from_top)
    end if

    if(DEBUG_NETCDF) write(*,*) "end serial part "
  end if! MasterProc

  !NB: this part is run by all processors
  !make lon and lat
  if(UsedProjection=='lon lat') then
    ! lon and lat are 1 dimensional, already done
  else
    !write one block at a time
    Buff2D(1:limax,1:ljmax,1)=glon(1:limax,1:ljmax)
    Buff2D(1:limax,1:ljmax,2)=glat(1:limax,1:ljmax)
    if(MasterProc)then
      do iproc=0,NPROC-1
        if(iproc>0)&
          CALL MPI_RECV(Buff2D,2*8*MAXLIMAX*MAXLJMAX, MPI_BYTE,iproc,&
                        iproc,MPI_COMM_CALC,MPISTATUS,IERROR)
        !NB IBEGcdf, IRUNBEG, etc. relative to fulldomain
        !    tgi0,tlimax, etc. relative to rundomain

        i1=1                  ! buffer, rundomain
        j1=1                  ! buffer, rundomain
        istart=tgi0(iproc)    ! restricted (IBEGcdf, JBEGcdf) domain
        jstart=tgj0(iproc)    ! restricted domain
        icount=tlimax(iproc)  ! both 
        jcount=tljmax(iproc)  ! both 

        ! must shift from rundomain to restricted grid in case it is smaller
        ishift=IBEGcdf-IRUNBEG
        jshift=JBEGcdf-JRUNBEG

        istart=istart-ishift
        jstart=jstart-jshift

        ! not entire buffer may be needed
        ! if the buffer start at the left of the restricted domain,
        ! must shift start of buffer and take less elements
        if(tgi0(iproc)<=ishift)then
          i1=i1+ishift-tgi0(iproc)+1
          istart=istart+(ishift-tgi0(iproc)+1)
          icount=icount-(ishift-tgi0(iproc)+1)
        end if
        if(tgj0(iproc)<=jshift)then
          j1=j1+jshift-tgj0(iproc)+1
          jstart=jstart+(jshift-tgj0(iproc)+1)
          jcount=jcount-(jshift-tgj0(iproc)+1)
        end if
        ! if we end at the right of the buffer end. must shift end of buffer,
        ! i.e. take less element
        if(istart+icount-1>GIMAXcdf)icount=GIMAXcdf-istart+1
        if(jstart+jcount-1>GJMAXcdf)jcount=GJMAXcdf-jstart+1

        if(all([istart<=GIMAXcdf+ishift,jstart<=GJMAXcdf+jshift,&
               icount>0,jcount>0,i1<=tlimax(iproc),j1<=tljmax(iproc)]))then

          call check(nf90_put_var(ncFileID,longVarID,   &
               Buff2D(i1:i1+icount-1,j1:j1+jcount-1,1), &
               start=[istart,jstart],count=[icount,jcount]))
          call check(nf90_put_var(ncFileID, latVarID,   &
               Buff2D(i1:i1+icount-1,j1:j1+jcount-1,2), &
               start=[istart,jstart],count=[icount,jcount]))
        end if
      end do
      if(DEBUG_NETCDF) write(*,*) "NetCDF: lon lat written"
    else
      CALL MPI_SEND(Buff2D,2*8*MAXLIMAX*MAXLJMAX,MPI_BYTE,0,&
                    me,MPI_COMM_CALC,IERROR)
    end if
  end if


  if(MasterProc)then
     call check(nf90_close(ncFileID))
     if(DEBUG_NETCDF)write(*,*)'NetCDF: file created, end of CreatenetCDFfile ',ncFileID
  end if
contains
function define_var(vname,xtype,dimIDs) result(varID)
  character(len=*),  intent(in) :: vname
  integer, intent(in) :: xtype,dimIDs(:)
  integer :: varID

  select case(vname)
  case("i")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name","projection_x_coordinate"))
    call check(nf90_put_att(ncFileID,varID,"coord_axis","x"))
    call check(nf90_put_att(ncFileID,varID,"long_name","EMEP grid x coordinate"))
    call check(nf90_put_att(ncFileID,varID,"units","km"))
  case("j")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name", "projection_y_coordinate"))
    call check(nf90_put_att(ncFileID,varID,"coord_axis","y"))
    call check(nf90_put_att(ncFileID,varID,"long_name","EMEP grid y coordinate"))
    call check(nf90_put_att(ncFileID,varID,"units","km"))
  case("i_EMEP")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"long_name","official EMEP grid coordinate i"))
    call check(nf90_put_att(ncFileID,varID,"units","gridcells"))
  case("j_EMEP")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"long_name","official EMEP grid coordinate j"))
    call check(nf90_put_att(ncFileID,varID,"units","gridcells"))
  case("i_RotS")
    call check(nf90_def_var(ncFileID,"i",xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name","grid_longitude"))
    call check(nf90_put_att(ncFileID,varID,"long_name","Rotated longitude"))
    call check(nf90_put_att(ncFileID,varID,"units","degrees"))
    call check(nf90_put_att(ncFileID,varID,"axis","X"))
  case("j_RotS")
    call check(nf90_def_var(ncFileID,"j",xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name","grid_latitude"))
    call check(nf90_put_att(ncFileID,varID,"long_name","Rotated latitude"))
    call check(nf90_put_att(ncFileID,varID,"units","degrees"))
    call check(nf90_put_att(ncFileID,varID,"axis","Y"))
  case("i_lambert")
    call check(nf90_def_var(ncFileID,"i",xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name","projection_x_coordinate"))
    call check(nf90_put_att(ncFileID,varID,"long_name","x-coordinate in Cartesian system"))
    call check(nf90_put_att(ncFileID,varID,"units","m"))
    call check(nf90_put_att(ncFileID,varID,"axis","X"))
  case("j_lambert")
    call check(nf90_def_var(ncFileID,"j",xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name","projection_y_coordinate"))
    call check(nf90_put_att(ncFileID,varID,"long_name","y-coordinate in Cartesian system"))
    call check(nf90_put_att(ncFileID,varID,"units","m"))
    call check(nf90_put_att(ncFileID,varID,"axis","Y"))
  case("lat")
    call check(nf90_def_var(ncFileID,"lat",xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"long_name","latitude"))
    call check(nf90_put_att(ncFileID,varID,"units","degrees_north"))
    call check(nf90_put_att(ncFileID,varID,"standard_name","latitude"))
  case("lon")
    call check(nf90_def_var(ncFileID,"lon",xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name","longitude"))
    call check(nf90_put_att(ncFileID,varID,"long_name","longitude"))
    call check(nf90_put_att(ncFileID,varID,"units","degrees_east"))
  case("lev") ! The hybrid sigma-pressure coordinate for level k is defined as ap(k)/p0+b(k). 
              ! p(n,k,j,i) = ap(k)+ b(k)*ps(n,j,i)
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
    call check(nf90_put_att(ncFileID,varID,"long_name", "hybrid level at layer midpoints (A/P0+B)"))
    call check(nf90_put_att(ncFileID,varID,"positive", "down"))
    call check(nf90_put_att(ncFileID,varID,"formula_terms","ap: hyam b: hybm ps: PS p0: P0"))
  case("ilev")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
    call check(nf90_put_att(ncFileID,varID,"long_name", "hybrid level at layer interfaces (A/P0+B)"))
    call check(nf90_put_att(ncFileID,varID,"positive", "down"))
    call check(nf90_put_att(ncFileID,varID,"formula_terms","ap: hyai b: hybi ps: PS p0: P0"))
  case("hyam")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"long_name","hybrid A coefficient at layer midpoints"))
    call check(nf90_put_att(ncFileID,varID,"units","hPa"))
  case("hyai")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"long_name","hybrid A coefficient at layer interfaces"))
    call check(nf90_put_att(ncFileID,varID,"units","hPa"))
  case("hybm")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID, "long_name","hybrid B coefficient at layer midpoints"))
  case("hybi")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"long_name","hybrid B coefficient at layer interfaces"))
  case("P0")
    call check(nf90_def_var(ncFileID,vname,xtype,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"units","hPa"))
  case("x_dist")
    call check(nf90_def_var(ncFileID,vname,xtype,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"long_name","displacement in x direction"))
  case("y_dist")
    call check(nf90_def_var(ncFileID,vname,xtype,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"long_name","displacement in y direction"))
  case("time")
    call check(nf90_def_var(ncFileID,vname,xtype,dimIDs,varID),"def:"//trim(vname))
    select case(period_type)
    case('instant','hourly','unknown')
      call check(nf90_put_att(ncFileID,varID,"long_name","time at end of period"))
    case default
      call check(nf90_put_att(ncFileID,varID,"long_name","time at middle of period"))
    end select
    call check(nf90_put_att(ncFileID,varID,"units","days since 1900-1-1 0:0:0"))
  case('Polar_Stereographic')
    call check(nf90_def_var(ncFileID,vname,xtype,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"grid_mapping_name", "polar_stereographic"))
    call check(nf90_put_att(ncFileID,varID,"straight_vertical_longitude_from_pole",Fi))
    call check(nf90_put_att(ncFileID,varID,"latitude_of_projection_origin",90.0))
    call check(nf90_put_att(ncFileID,varID,"scale_factor_at_projection_origin",scale_at_projection_origin))
  case('Rotated_Spherical')
    call check(nf90_def_var(ncFileID,vname,xtype,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"grid_mapping_name","rotated_latitude_longitude"))
    call check(nf90_put_att(ncFileID,varID,"grid_north_pole_latitude",grid_north_pole_latitude))
    call check(nf90_put_att(ncFileID,varID,"grid_north_pole_longitude",grid_north_pole_longitude))
  case('rotated_pole')
    call check(nf90_def_var(ncFileID,vname,xtype,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"grid_mapping_name","rotated_latitude_longitude"))
    call check(nf90_put_att(ncFileID,varID,"grid_north_pole_latitude",grid_north_pole_latitude))
    call check(nf90_put_att(ncFileID,varID,"grid_north_pole_longitude", grid_north_pole_longitude))
  case('projection_lambert')
    call check(nf90_def_var(ncFileID,vname,xtype,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"grid_mapping_name","lambert_conformal_conic"))
    call check(nf90_put_att(ncFileID,varID,"standard_parallel",(/lat_stand1_lambert,lat_stand2_lambert/)))
    call check(nf90_put_att(ncFileID,varID,"longitude_of_central_meridian",lon0_lambert))
    call check(nf90_put_att(ncFileID,varID,"latitude_of_projection_origin",lat0_lambert))
    call check(nf90_put_att(ncFileID,varID,"earth_radius",earth_radius))!NB: reset to the same as in metdata
  case(Default_projection_name)
    call check(nf90_def_var(ncFileID,vname,xtype,varID),"def:"//trim(vname))
    call check(nf90_put_att(ncFileID,varID,"grid_mapping_name",trim(UsedProjection)))
  case default
    call CheckStop("CreatenetCDFfile: unknown metadata variable "//trim(vname))
  end select
end function define_var
subroutine write_klev(kmax,klev,from_top)
  integer, intent(in) :: kmax,klev(kmax)
  logical, intent(in) :: from_top
  real :: Am(kmax),Bm(kmax)
  real :: Ai(kmax+1),Bi(kmax+1)
  integer :: k,varID
  if(from_top)then            ! from model top to surface (as model levels)
    do k=1,kmax
      if(klev(k)==0)then                  ! 0-->surface
        Am(k)=A_bnd(KMAX_BND) ! definition of level ambiguous, since no thickness
        Bm(k)=B_bnd(KMAX_BND)
        Ai(k)=A_bnd(KMAX_BND)
        Bi(k)=B_bnd(KMAX_BND)
      else                                !20-->20;19-->19;...;1-->1
        Am(k)=A_mid(klev(k))
        Bm(k)=B_mid(klev(k))
        Ai(k)=A_bnd(klev(k))
        Bi(k)=B_bnd(klev(k))
        if(DEBUG_NETCDF) write(*,*) "TESTHH netcdf KLEVcdf ",k,klev(k),Am(k) 
      end if
    end do
    Ai(kmax+1)=A_bnd(klev(kmax)+1)        ! top boundary
    Bi(kmax+1)=B_bnd(klev(kmax)+1)
  else                        ! from surface to model top (REVERSE order)
    do k=1,kmax
      if(klev(k)==0)then                  ! 0-->surface
        Am(k)=A_bnd(KMAX_BND) ! definition of level ambiguous, since no thickness
        Bm(k)=B_bnd(KMAX_BND)
        Ai(k)=A_bnd(KMAX_BND)
        Bi(k)=B_bnd(KMAX_BND)
      else                                ! 1-->20;2-->19;...;20-->1
        Am(k)=A_mid(KMAX_MID-klev(k)+1)
        Bm(k)=B_mid(KMAX_MID-klev(k)+1)
        Ai(k)=A_bnd(KMAX_BND-klev(k)+1)
        Bi(k)=B_bnd(KMAX_BND-klev(k)+1)
        if(DEBUG_NETCDF) write(*,*) "TESTHH netcdf KLEVcdf ",k,klev(k),Am(k) 
      end if
    end do
    Ai(kmax+1)=A_bnd(KMAX_BND-klev(kmax)) ! top boundary
    Bi(kmax+1)=B_bnd(KMAX_BND-klev(kmax))
  end if
  call check(nf90_inq_varid(ncFileID,"P0"  ,varID)  ,"inq:P0")
  call check(nf90_put_var(ncFileID,varID,Pref/100.0),"put:P0")
  call check(nf90_inq_varid(ncFileID,"hyam",varID)  ,"inq:hyam")
  call check(nf90_put_var(ncFileID,varID,Am/100.0)  ,"put:hyam")
  call check(nf90_inq_varid(ncFileID,"hyai",varID)  ,"inq:hyai")
  call check(nf90_put_var(ncFileID,varID,Ai/100.0)  ,"put:hyai")
  call check(nf90_inq_varid(ncFileID,"hybm",varID)  ,"inq:hybm")
  call check(nf90_put_var(ncFileID,varID,Bm)        ,"put:hybm")
  call check(nf90_inq_varid(ncFileID,"hybi",varID)  ,"inq:hybi")
  call check(nf90_put_var(ncFileID,varID,Bi)        ,"put:hybi")
  call check(nf90_inq_varid(ncFileID, "lev",varID)  ,"inq:lev")
  call check(nf90_put_var(ncFileID,varID,Am/Pref+Bm),"put:lev")
  call check(nf90_inq_varid(ncFileID,"ilev",varID)  ,"inq:ilev")
  call check(nf90_put_var(ncFileID,varID,Ai/Pref+Bi),"put:ilev")
end subroutine write_klev
end subroutine CreatenetCDFfile
!_______________________________________________________________________
subroutine Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype,dimSizes,dimNames,out_DOMAIN,ik,&
                      fileName_given,overwrite,create_var_only,chunksizes,ncfileID_given)
!The use of fileName_given is probably slower than the implicit filename used by defining iotyp.

! Mandatary arguments:
  integer ,    intent(in) :: ndim,kmax
  type(Deriv), intent(in) :: def1 ! definition of fields
  integer,     intent(in) :: iotyp
  real,        intent(in) :: scale
  real, dimension(*), intent(in) :: dat ! Data arrays
! Optional arguments:
  integer, optional, intent(in) :: &
       dimSizes(ndim),&
    out_DOMAIN(4),ik,& ! Output subdomain. Only level ik is written if defined
    CDFtype            != OUTtype. (Integer*1, Integer*2,Integer*4, real*8 or real*4)
  character(len=*),intent(in), optional :: dimNames(ndim)
  character (len=*),optional, intent(in) :: &
    fileName_given ! filename to which the data must be written
                   !NB if the file fileName_given exist (also from earlier runs) it will be appended
  logical, optional, intent(in) :: &
    overwrite,      &   ! overwrite if file already exists (in case fileName_given)
    create_var_only     ! only create the variable, without writing the data content
  integer, dimension(ndim), intent(in), optional :: &
    chunksizes          ! nc4zip output written in slices, see NETCDF_DEFLATE_LEVEL
  integer, optional, intent(inout) :: ncFileID_given !if present, do not close the file at return
                                                     !if >=0 at input, the file is already open 
                                                     !       and  ncFileID=ncFileID_given
  logical:: create_var_only_local !only create the variable, without writing the data content

  character(len=len(def1%name)) :: varname
  character(len=08) :: lastmodified_date
  character(len=10) :: lastmodified_hour,lastmodified_hour0,created_hour
  integer :: varID,nrecords,ncFileID=closedID,ndate(4)
  integer :: d,alloc_err,ijk,status,i,j,k,i1,i2,j1,j2
  real :: buff(MAXLIMAX*MAXLJMAX*KMAX_MID)
  real(kind=8),   allocatable,dimension(:,:,:) :: R8data3D
  real(kind=4),   allocatable,dimension(:,:,:) :: R4data3D
  integer(kind=4),allocatable,dimension(:,:,:) :: Idata3D
  integer :: OUTtype !local version of CDFtype
  integer :: iotyp_new
  integer :: iDimID,jDimID,kDimID,timeVarID
  integer :: GIMAX_old,GJMAX_old,KMAX_old
  integer :: GIMAXcdf,GJMAXcdf,IBEGcdf,JBEGcdf
  real(kind=8) :: rdays,rdays_time(1),rdaysstart
  logical :: overwrite_local,createfile=.false.
  integer, parameter :: IOU_GIVEN=-IOU_INST
  integer ::domain(4),startvec(10),countvec(10),Nextradim,n,iextradim,iiextradim,nijk,date_start(4)

  domain=RUNDOMAIN!default domain (in fulldoamin coordinates)
!fullrun, Monthly, Daily and hourly domains may be predefined
  select case(iotyp)
  case(IOU_YEAR)
    domain = fullrun_DOMAIN
  case(IOU_MON)
    domain = month_DOMAIN
  case(IOU_DAY)
    domain = day_DOMAIN
  case(IOU_HOUR:IOU_HOUR_INST)
    domain = hour_DOMAIN
  end select
  if(present(out_DOMAIN)) domain = out_DOMAIN
  !convert into rundomain coordinates
  i1=domain(1)-IRUNBEG+1
  i2=domain(2)-IRUNBEG+1
  j1=domain(3)-JRUNBEG+1
  j2=domain(4)-JRUNBEG+1

  create_var_only_local=.false.
  if(present(create_var_only))create_var_only_local=create_var_only
  !Check that that the area is larger than 0
  if((i2-i1)<0.or.(j2-j1)<0.or.kmax<=0)then
     if(MasterProc) write(*,*)'WARNING: requested boundaries inconsistent '
     if(MasterProc) write(*,*) i1,i2,j1,j2,kmax
     return
  end if

 !make variable name
  write(varname,fmt='(A)')trim(def1%name)
  if(DEBUG_NETCDF.and.MasterProc) write(*,*)'Out_NetCDF: START ',trim(varname)

  !to shorten the output we can save only the components explicitely named here
  !if(varname.ne.'D2_NO2'.and.varname.ne.'D2_O3' &
  !                         .and.varname.ne.'D2_PM10')return

  !do not write 3D fields (in order to shorten outputs)
  !if(ndim==3)return

  iotyp_new=iotyp
  createfile=.false.
  if(present(fileName_given))then
    !NB if the file already exist (also from earlier runs) it will be appended
    overwrite_local=.false.
    if(present(overwrite))overwrite_local=overwrite
    if(MasterProc)then
       if(overwrite_local)then
         if(DEBUG_NETCDF) &
           write(*,*)'Out_NetCDF: overwrite file (if exists) ',trim(fileName_given)
         ! do not try to open the file
         status=nf90_noerr
       elseif(present(ncFileID_given))then
         if(DEBUG_NETCDF) &
           write(*,*)'Out_NetCDF: ncFileID_given ',ncFileID_given
         if(ncFileID_given<0)then
           status=nf90_open(trim(fileName_given),nf90_share+nf90_write,ncFileID)
           ncFileID_given=ncFileID         
           if(DEBUG_NETCDF) &
             write(*,*)'Out_NetCDF: opened a file ',ncFileID_given
         else
           !the file must already be open
           ncFileID=ncFileID_given
           status=nf90_noerr
           if(DEBUG_NETCDF) &
             write(*,*)'Out_NetCDF: assuming file already open ',trim(fileName_given)
         end if
       else
         status=nf90_open(trim(fileName_given),nf90_write,ncFileID)
       end if
       if(DEBUG_NETCDF) write(*,*)'Out_NetCDF: fileName_given ' ,&
            trim(fileName_given),overwrite_local,status==nf90_noerr,ncfileID,&
            trim(nf90_strerror(status))
       createfile=overwrite_local.or.(status/=nf90_noerr)
    end if
    CALL MPI_BCAST(createfile ,1,MPI_LOGICAL,0,MPI_COMM_CALC,IERROR)
    
    IBEGcdf=IRUNBEG+i1-1
    JBEGcdf=JRUNBEG+j1-1
    GIMAXcdf=i2-i1+1
    GJMAXcdf=j2-j1+1

    if(createfile) then ! the file does not exist yet or is overwritten
       
      if(MasterProc)write(6,fmt='(A,12I6)') 'creating file: '//trim(fileName_given)//" with sizes ",GIMAXcdf,GJMAXcdf,KMAX
      period_type = 'unknown'
      call CreatenetCDFfile(trim(fileName_given),GIMAXcdf,GJMAXcdf,IBEGcdf,JBEGcdf,KMAX)
      if(present(ncFileID_given))then
        !the file should be opened, but by MasterProc only
        CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)!wait until the file creation is finished     
        if(MasterProc)then
          call check(nf90_open(trim(fileName_given),nf90_share+nf90_write,ncFileID))
        else
          ncFileID=closedID
        end if
      else
        ncFileID=closedID
      end if
   elseif(MasterProc)then
      !test if the defined dimensions are compatible
        if(DEBUG_NETCDF)then
           write(6,*) 'checking dims file: ',trim(fileName_given),ncFileID
            ! find main properties
            call check(nf90_Inquire(ncFileID,n,i,j))
            print *, trim(fileName),ncfileid,' properties: '
            print *, 'Nb of dimensions: ',n
            print *, 'Nb of variables: ',i
         endif
      select case(projection)
      case('lon lat')
         call check(nf90_inq_dimid(ncFileID,"lon",idimID),"dim:lon")
         call check(nf90_inq_dimid(ncFileID,"lat",jdimID),"dim:lat")
      case default
        call check(nf90_inq_dimid(ncFileID,"i"  ,idimID),"dim:i")
        call check(nf90_inq_dimid(ncFileID,"j"  ,jdimID),"dim:j")
      end select

      if(USE_EtaCOORDINATES)then
        call check(nf90_inq_dimid(ncFileID,"lev",kdimID),"dim:lev")
      else
        call check(nf90_inq_dimid(ncFileID,"k"  ,kdimID),"dim:k")
      end if
      ! only i,j coords can be handled for PS so far.
      ! Possible x,y would give wrong dimID. 
      ! Check if all dims are found:
      call CheckStop(any([idimID,jdimID,kdimID]<0),&
           "Out_NetCDF: no dimID found for"//trim(fileName_given))
      
      call check(nf90_inquire_dimension(ncFileID,idimID,len=GIMAX_old),"len:i")
      call check(nf90_inquire_dimension(ncFileID,jdimID,len=GJMAX_old),"len:j")
      call check(nf90_inquire_dimension(ncFileID,kdimID,len=KMAX_old) ,"len:k")
      
      if(any([GIMAX_old,GJMAX_old,KMAX_old]<[GIMAXcdf,GJMAXcdf,KMAX]))then
        write(6,*)'existing file ', trim(fileName_given),' has wrong dimensions'
        write(6,*)GIMAX_old,GIMAXcdf,GJMAX_old,GJMAXcdf,KMAX_old,KMAX
        write(6,*)'WARNING! OLD ', trim(fileName_given),' MUST BE DELETED'
        Call StopAll('outCDF, inconsistent file size ')
        !write(6,*)'creating new file: ',trim(fileName_given)
        !period_type = 'unknown'
        !NB: this cannot be used directly, because it must be called by all processors
        !call CreatenetCDFfile(trim(fileName_given),GIMAXcdf,GJMAXcdf,&
        !     IBEGcdf,JBEGcdf,KMAX)
        !ncFileID=closedID
      end if
    end if
   
    iotyp_new=IOU_GIVEN
    ncFileID_new=ncFileID
  end if

  if(DEBUG_NETCDF.and.MasterProc)then
    if(iotyp_new==IOU_GIVEN)then
      write(*,*)' Out_NetCDF: cases new file ', trim(fileName_given), iotyp
    else
      write(*,*)' Out_NetCDF: cases old file ', trim(fileName), iotyp
    end if
  end if

  select case(iotyp_new)
  case(IOU_GIVEN)
    fileName = trim(fileName_given)
    ncFileID = ncFileID_new
  case(IOU_INST:IOU_HOUR_INST)
    fileName = trim(fileName_iou(iotyp_new))
    ncFileID = ncFileID_iou(iotyp_new)
  case default
    return
  end select

  if(present(ncFileID_given))ncFileID_given=ncFileID!use rather stored ncFileID_XXX

  if(DEBUG_NETCDF.and.MasterProc) &
    write(*,*)'Out_NetCDF, filename ', trim(fileName), iotyp,ncFileID

  if(.not. present(dimSizes))then 
     call CheckStop(ndim,[2,3], "NetCDF_mod: ndim must be 2 or 3")
  endif

  OUTtype=Real4  !default value
  if(present(CDFtype))OUTtype=CDFtype

  if(MasterProc)then
    ndate(1:4) = [current_date%year,current_date%month,&
                  current_date%day ,current_date%hour]

    !test if the file is already open
    if(ncFileID==closedID)then
      !open an existing netcdf dataset
      call check(nf90_open(fileName,nf90_write,ncFileID),"nf90_open"//trim(fileName))
      select case(iotyp_new)      ! needed in case iotyp is defined
      case(IOU_GIVEN)
        ncFileID_new  = ncFileID  ! not really needed
      case(IOU_INST:IOU_HOUR_INST)
        ncFileID_iou(iotyp_new)=ncFileID
      end select
    end if

    !test first if the variable is already defined:
    status=nf90_inq_varid(ncFileID,varname,VarID)
    if(status==nf90_noerr)then
!     print *, 'variable exists: ',varname
      if(DEBUG_NETCDF) write(*,*) 'Out_NetCDF: variable exists: ',varname
    else
      if(DEBUG_NETCDF) write(*,*) 'Out_NetCDF: creating variable: ',varname
      if(create_var_only_local) &
           call check(nf90_set_fill(ncFileID,NF90_NOFILL,ijk),"nofill:"//trim(varname))
      if(.not.present(dimSizes).and.present(chunksizes))&
           call CheckStop(chunksizes(1)/=(i2-i1+1).or.chunksizes(2)/=(j2-j1+1),&
           "NetCDF_mod: chunksizes has wrong dimensions")
      call createnewvariable(ncFileID,varname,ndim,ndate,def1,OUTtype,chunksizes=chunksizes,&
              dimSizes=dimSizes,dimNames=dimNames)
    end if
  end if!MasterProc

  if(create_var_only_local)then
    ! Don't write the data
    ! For performance: need to create all variables before writing data
    if(MasterProc.and.iotyp_new==IOU_GIVEN)then
      if(present(ncFileID_given))then
        if(DEBUG_NETCDF)write(*,*)'keep open ',trim(fileName),ncFileID
        ncFileID_given=ncFileID
      else
        call check(nf90_close(ncFileID),"close:"//trim(fileName))
      end if
    end if
    if(DEBUG_NETCDF.and.MasterProc)write(*,*)'variable ONLY created. Finished'
    CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)!wait until the file creation is finished     
   return
  end if ! create var only

  Nextradim=1
  if(present(dimSizes))then 
     do n=1,ndim-3
        Nextradim=Nextradim*dimSizes(n)
     end do
  endif
  do iextradim = 1,Nextradim
     startvec(1:ndim+1)=1
     countvec(1:ndim+1)=1
     iiextradim = iextradim-1
     do n=1,ndim-3
        startvec(n)=mod(iiextradim,dimSizes(n))+1
        iiextradim=iiextradim/dimSizes(n)
     end do
     if(ndim-2>0)countvec(ndim-2)=i2-i1+1
     if(ndim-1>0)countvec(ndim-1)=j2-j1+1
     countvec(ndim)=kmax
   !buffer the wanted part of data
     ijk=0
     nijk=iextradim-Nextradim
     do k=1,kmax
        do j = 1,tljmax(me)
           do i = 1,tlimax(me)
              ijk=ijk+1
              nijk=nijk+Nextradim
!              if(i==5.and.j==5.and.k==kmax.and.mod(iextradim-1,11)==0.and.dat(nijk)>1.E-12)write(*,*)trim(varname),me,iextradim,dat(nijk)
              buff(ijk)=dat(nijk)*scale
              !if(isnan(buff(ijk)))then
              !   write(*,*)'ERROR ',me,i,j,k,iextradim,ijk
              !   stop
              !endif
              !if(buff(ijk)>1.1)then
              !   write(*,*)'too large ERROR ',trim(varname),me,i,j,k,iextradim,ijk,buff(ijk)
              !   stop
              !endif
           end do
        end do
     end do
     
     !send all data to me=0
     outCDFtag=outCDFtag+1
     
     if(MasterProc)then
        if(iextradim==1)then
           ! allocate a large array (only on one processor)
           select case(OUTtype)
           case(Int1,Int2,Int4)
              allocate(Idata3D (GIMAX,GJMAX,kmax),stat=alloc_err)
           case(Real4)
              allocate(R4data3D(GIMAX,GJMAX,kmax),stat=alloc_err)
           case(Real8)
              allocate(R8data3D(GIMAX,GJMAX,kmax),stat=alloc_err)
           case default
              WRITE(*,*)'WARNING NetCDF:Data type not supported'
              alloc_err=-1  ! trip error stop
           end select
           call CheckStop(alloc_err, "alloc failed in NetCDF_mod")
        endif
        !write own data in global array
        select case(OUTtype)
        case(Int1,Int2,Int4)
           ijk=0
           do k=1,kmax
              do j = tgj0(me),tgj0(me)+tljmax(me)-1
                 do i = tgi0(me),tgi0(me)+tlimax(me)-1
                    ijk=ijk+1
                    Idata3D(i,j,k)=buff(ijk)
                 end do
              end do
           end do
        case(Real4)
           ijk=0
           do k=1,kmax
              do j = tgj0(me),tgj0(me)+tljmax(me)-1
                 do i = tgi0(me),tgi0(me)+tlimax(me)-1
                    ijk=ijk+1
                    R4data3D(i,j,k)=buff(ijk)
                 end do
              end do
           end do
        case(Real8)
           ijk=0
           do k=1,kmax
              do j = tgj0(me),tgj0(me)+tljmax(me)-1
                 do i = tgi0(me),tgi0(me)+tlimax(me)-1
                    ijk=ijk+1
                    R8data3D(i,j,k)=buff(ijk)
                 end do
              end do
           end do
        end select
        
        do d = 1, NPROC-1
           CALL MPI_RECV(buff, 8*tlimax(d)*tljmax(d)*kmax, MPI_BYTE, d, &
                outCDFtag, MPI_COMM_CALC, MPISTATUS, IERROR)
           
           ! copy data to global buffer
           select case(OUTtype)
           case(Int1,Int2,Int4)
              ijk=0
              do k=1,kmax
                 do j = tgj0(d),tgj0(d)+tljmax(d)-1
                    do i = tgi0(d),tgi0(d)+tlimax(d)-1
                       ijk=ijk+1
                       Idata3D(i,j,k)=buff(ijk)
                    end do
                 end do
              end do
           case(Real4)
              ijk=0
              do k=1,kmax
                 do j = tgj0(d),tgj0(d)+tljmax(d)-1
                    do i = tgi0(d),tgi0(d)+tlimax(d)-1
                       ijk=ijk+1
                       R4data3D(i,j,k)=buff(ijk)
                    end do
                 end do
              end do
           case(Real8)
              ijk=0
              do k=1,kmax
                 do j = tgj0(d),tgj0(d)+tljmax(d)-1
                    do i = tgi0(d),tgi0(d)+tlimax(d)-1
                       ijk=ijk+1
                       R8data3D(i,j,k)=buff(ijk)
                    end do
                 end do
              end do
           end select
        end do
     else
        CALL MPI_SEND( buff, 8*tlimax(me)*tljmax(me)*kmax, MPI_BYTE, 0, &
             outCDFtag, MPI_COMM_CALC, IERROR)
     end if
     !return

     if(MasterProc)then
        if(iextradim==1)then
        ndate(1:4) = [current_date%year,current_date%month,&
             current_date%day ,current_date%hour]
        
        ! get variable id
        call check(nf90_inq_varid(ncFileID,varname,VarID))
        ! find the number of records already written
        call check(nf90_get_att(ncFileID,VarID,"numberofrecords",nrecords))
        if(DEBUG_NETCDF) print *,'number of dataset saved: ',nrecords
        
        ! test if new record is needed
        if(present(ik).and.nrecords>0 .and. iextradim==1)then
           ! The new record may already exist
           call date2nctime(current_date,rdays,iotyp)
           call check(nf90_inq_varid(ncFileID,"time",timeVarID))
           call check(nf90_get_var(ncFileID,timeVarID,rdays_time,start=(/nrecords/)))
           ! check if this is a newer time
           if((abs(rdays-rdays_time(1))>0.00001))then!0.00001 is about 1 second
              nrecords=nrecords+1 !start a new record
           end if
        else
           ! increase nrecords, to define position of new data
           nrecords=nrecords+1
        end if
        if(DEBUG_NETCDF) print *,'writing on dataset: ',nrecords
        endif
        startvec(ndim+1)=nrecords
       
        ! append new values
        select case(OUTtype)
        case(Int1,Int2,Int4)  ! type Integer
           if(ndim==3)then
              if(present(ik))then
                 !     print *, 'write: ',i1,i2, j1,j2,ik
                 call check(nf90_put_var(ncFileID,VarID,Idata3D(i1:i2,j1:j2,1),&
                      start=(/1,1,ik,nrecords/)))
              else
                 call check(nf90_put_var(ncFileID,VarID,Idata3D(i1:i2,j1:j2,1:kmax),&
                      start=(/1,1,1,nrecords/)))
              end if
           else if(ndim==2)then
              call check(nf90_put_var(ncFileID, VarID,Idata3D(i1:i2,j1:j2,1),&
                   start=(/1,1,nrecords/)))
           else
              call check(nf90_put_var(ncFileID, VarID,Idata3D(i1:i2,j1:j2,1:kmax),&
                   start=startvec(1:ndim+1),count=countvec(1:ndim+1)))
           end if
           
        case(Real4)           ! type Real4
           if(ndim==3)then
              if(present(ik))then
                 !     print *, 'write: ',i1,i2, j1,j2,ik
                 call check(nf90_put_var(ncFileID,VarID,R4data3D(i1:i2,j1:j2,1),&
                      start=(/1,1,ik,nrecords/)))
              else
                 call check(nf90_put_var(ncFileID,VarID,R4data3D(i1:i2,j1:j2,1:kmax),&
                      start=(/1,1,1,nrecords/)))
              end if
           else if(ndim==2)then
              call check(nf90_put_var(ncFileID,VarID,R4data3D(i1:i2,j1:j2,1),&
                   start=(/1,1,nrecords/)))
           else
              !write(*,*)'writing slice ',iextradim,' of ',Nextradim
!              if(mod(iextradim-1,12)==0)then
!                 write(*,*)trim(varname),me,startvec(2),startvec(3),iextradim,R4data3D(GIMAX/2,GJMAX/2,kmax)
!              endif
    !if(me>=0)write(*,*)iextradim,'outcdf locfrac ',R4data3D(4,37,7),trim(varname),trim(fileName)
              call check(nf90_put_var(ncFileID, VarID,R4data3D(i1:i2,j1:j2,1:kmax),&
                   start=startvec(1:ndim+1),count=countvec(1:ndim+1)))
          end if
           
        case(Real8)           ! type Real8
           if(ndim==3)then
              if(present(ik))then
                 !     print *, 'write: ',i1,i2, j1,j2,ik
                 call check(nf90_put_var(ncFileID,VarID,R8data3D(i1:i2,j1:j2,1),&
                      start=(/1,1,ik,nrecords/)))
              else
                 call check(nf90_put_var(ncFileID,VarID,R8data3D(i1:i2,j1:j2,1:kmax),&
                      start=(/1,1,1,nrecords/)))
              end if
           else if(ndim==2)then
              call check(nf90_put_var(ncFileID,VarID,R8data3D(i1:i2,j1:j2,1),&
                   start=(/1,1,nrecords/)))
           else
              call check(nf90_put_var(ncFileID, VarID,R8data3D(i1:i2,j1:j2,1:kmax),&
                   start=startvec(1:ndim+1),count=countvec(1:ndim+1)))
           end if
           
        end select !type
     end if !MasterProc
  enddo

  if(MasterProc)then
     select case(OUTtype)
     case(Int1,Int2,Int4)  ! type Integer
        deallocate(Idata3D,stat=alloc_err)
        call CheckStop(alloc_err, "dealloc failed in NetCDF_mod")
     case(Real4)           ! type Real4
        deallocate(R4data3D, stat=alloc_err)
        call CheckStop(alloc_err, "dealloc failed in NetCDF_mod")
     case(Real8)           ! type Real8
        deallocate(R8data3D, stat=alloc_err)
        call CheckStop(alloc_err, "dealloc failed in NetCDF_mod")
     end select !type
     
     call check(nf90_get_att(ncFileID,nf90_global,"lastmodified_hour",lastmodified_hour0 ))
     call check(nf90_get_att(ncFileID,nf90_global,"created_hour",created_hour))
     call Date_And_Time(date=lastmodified_date,time=lastmodified_hour)
     
     ! write or change attributes NB: strings must be of same length as originally
     call check(nf90_put_att(ncFileID,VarID,"numberofrecords",nrecords))
     
     ! update dates
     call check(nf90_put_att(ncFileID,nf90_global,"lastmodified_date",lastmodified_date))
     call check(nf90_put_att(ncFileID,nf90_global,"lastmodified_hour",lastmodified_hour))
     call check(nf90_put_att(ncFileID,VarID, "current_date_last",ndate))
     
     ! update time dim
     call check(nf90_inq_varid(ncFileID,"time",VarID))
     if(iotyp==IOU_YEAR)then
        date_start = startdate!start of period        
        call date2nctime(date_start,rdaysstart)!start of period                
        call date2nctime(current_date,rdays)!now
        rdays = (rdaysstart+rdays)/2 !middle
     else if (iotyp==IOU_MON) then
        !careful month can be incomplete in first month or last month or both
        if((mod(12+current_date%month-startdate(2),12)==1 .and.&
             (current_date%day==1 .and. current_date%hour==0)).or.&
             (mod(12+current_date%month-startdate(2),12)==0).and.&
             (current_date%day/=1 .or. current_date%hour/=0))then       
           !first month outputted
           date_start = startdate!start of period        
        else
           !the start is at the start of the current month
           date_start(1) = current_date%year
           date_start(3) = 1
           date_start(4) = 0
           if(current_date%day==1 .and. current_date%hour==0)then
              date_start(2) = current_date%month-1
           else
              date_start(2) = current_date%month           
           endif
        endif
        call date2nctime(date_start,rdaysstart)!start of period                
        call date2nctime(current_date,rdays)!now
        rdays = (rdaysstart+rdays)/2 !middle        
     else
        call date2nctime(current_date,rdays,iotyp) !routine will subtract half hour, day  if necessary
     endif
     call check(nf90_put_var(ncFileID,VarID,rdays,start=(/nrecords/)))
     
     !close file if present(fileName_given)
     if(iotyp_new==IOU_GIVEN)then
        if(present(ncFileID_given))then
           if(DEBUG_NETCDF)write(*,*)'keep open ',trim(fileName),ncFileID
           ncFileID_given=ncFileID
        else
           call check(nf90_close(ncFileID))
        end if
     end if
  end if !MasterProc

  if(DEBUG_NETCDF.and.MasterProc) write(*,*)'Out_NetCDF: FINISHED '
end subroutine Out_netCDF
!_______________________________________________________________________
subroutine  createnewvariable(ncFileID,varname,ndim,ndate,def1,OUTtype,chunksizes,dimSizes,dimNames)
! create new netCDF variable

  implicit none

  type(Deriv),     intent(in) :: def1 ! definition of fields
  character(len=*),intent(in) :: varname
  integer,         intent(in) :: ndim,ncFileID,OUTtype
  integer,dimension(:),intent(in) :: ndate
  integer,dimension(ndim),intent(in), optional :: chunksizes,dimSizes
  character(len=*),intent(in), optional :: dimNames(ndim)
  integer :: iDimID,jDimID,kDimID,timeDimID,nDimID(10)
  integer ::n, i, isize, jsize
  integer :: varID,dimVarID,nrecords,status
  real :: scale
  integer :: OUTtypeCDF !NetCDF code for type
  real,allocatable, dimension(:) :: tmp 

  select case(OUTtype)
    case(Int1 );OUTtypeCDF=nf90_byte
    case(Int2 );OUTtypeCDF=nf90_short
    case(Int4 );OUTtypeCDF=nf90_int
    case(Real4);OUTtypeCDF=nf90_float
    case(Real8);OUTtypeCDF=nf90_double
    case default
    call CheckStop(MasterProc,"NetCDF_mod: undefined datatype")
  end select

  !define mode
  call check(nf90_redef(ncFileID),"file redef:"//trim(varname))

  !get dimensions id
  select case(projection)
  case('lon lat')
    call check(nf90_inq_dimid(ncFileID,"lon",idimID),"dim:lon")
    call check(nf90_inq_dimid(ncFileID,"lat",jdimID),"dim:lat")
  case default
    call check(nf90_inq_dimid(ncFileID,"i"  ,idimID),"dim:i")
    call check(nf90_inq_dimid(ncFileID,"j"  ,jdimID),"dim:j")
  end select

  status=nf90_inq_dimid(ncFileID,"k",kdimID)
  if(status/=nf90_noerr)&
    call check(nf90_inq_dimid(ncFileID,"lev",kdimID),"dim:lev")

  call check(nf90_inq_dimid(ncFileID,"time",timeDimID),"dim:time")

  !define new variable: nf90_def_var(ncid,name,xtype,dimids,varid)
  if(.not.present(dimSizes))then
     select case(ndim)
     case(3)
        call check(nf90_def_var(ncFileID,varname,OUTtypeCDF,&
             [iDimID,jDimID,kDimID,timeDimID],varID),"def3d:"//trim(varname))
     case(2)
        call check(nf90_def_var(ncFileID,varname,OUTtypeCDF,&
             [iDimID,jDimID,timeDimID]       ,varID),"def2d:"//trim(varname))
     case default
        print *, 'createnewvariable: unexpected ndim ',ndim
        call CheckStop(MasterProc,"NetCDF_mod: unexpected ndim")
     end select
  else
     !define dimensions if needed
     do n=1,ndim
        status = nf90_inq_dimid(ncFileID,dimNames(n),ndimID(n))
        if(status/=nf90_noerr)then
           !define new dimension
           !write(*,*)'defining new dimension: '//trim(dimNames(n))//' for '//trim(varname)
           call check(nf90_def_dim(ncFileID,trim(dimNames(n)),dimSizes(n),ndimID(n)),"dim:"//trim(dimNames(n)))
           call check(nf90_def_var(ncFileID,trim(dimNames(n)),nf90_float,ndimID(n),dimvarID),"defvar:"//trim(dimNames(n)))
 
           !call check(nf90_enddef(ncFileID))
           allocate(tmp(dimSizes(n)))
           do i=1,dimSizes(n)
              if(trim(dimNames(n))=='x_dist'.or.trim(dimNames(n))=='y_dist')then
                 tmp(i)=i-(dimSizes(n)+1)/2
              else if(trim(dimNames(n))=='klevel')then
                 tmp(i)=KMAX_MID+i-dimSizes(n)
              else
                 write(*,*)'ERROR: Dimension Name not recognized: '//trim(dimNames(n))
              endif
           enddo
           call check(nf90_put_var(ncFileID,dimVarID,tmp))
           deallocate(tmp)
           !call check(nf90_redef(ncFileID),"file redef:"//trim(varname))
           write(*,*)'defining new dimension: '//trim(dimNames(n)),', size ',dimSizes(n)
        endif
     enddo

     ndimID(ndim+1)=timeDimID
     call check(nf90_def_var(ncFileID,varname,OUTtypeCDF,&
          ndimID(1:ndim+1),varID),"defnd:"//trim(varname))     
  endif
!define variable as to be compressed
  if(NETCDF_DEFLATE_LEVEL >= 0) then
     call check(nf90_def_var_deflate(ncFileid,varID,shuffle=0,deflate=1,&
          deflate_level=NETCDF_DEFLATE_LEVEL),"compress:"//trim(varname))
     if(present(chunksizes))then      ! set chunk-size for 2d slices of 3d output
        ! write(*,*)' chunksizes ',trim(varname),chunksizes
        call check(nf90_def_var_chunking(ncFileID,varID,NF90_CHUNKED,&
             chunksizes(:)),"chunk:"//trim(varname))
     else
        if(ndim==2)then
           !Recommended by Heiko for faster verification script
           call check(nf90_inquire_dimension(ncFileID,idimID,len=isize))
           call check(nf90_inquire_dimension(ncFileID,jdimID,len=jsize))
           call check(nf90_def_var_chunking(ncFileID,varID,NF90_CHUNKED,&
!                (/300,130,1/)),"chunk2D:"//trim(varname))         
                (/min(isize,300),min(jsize,130),1/)),"chunk2D:"//trim(varname))         
        endif
     endif
  end if
  !     FillValue=0.
  scale=1.
  !define attributes of new variable
  call check(nf90_put_att(ncFileID,varID,"long_name",trim(def1%name)))
  call check(nf90_put_att(ncFileID,varID,"coordinates", "lat lon"))
  select case(projection)
  case('Stereographic')
    call check(nf90_put_att(ncFileID, varID, "grid_mapping", "Polar_Stereographic"))
  case('lon lat')
  case('Rotated_Spherical')
    call check(nf90_put_att(ncFileID, varID, "grid_mapping", "Rotated_Spherical"))
  case('lambert')
    call check(nf90_put_att(ncFileID, varID, "grid_mapping", "projection_lambert"))
  case default
    call check(nf90_put_att(ncFileID, varID, "grid_mapping",Default_projection_name ))
  end select

  nrecords=0
  call check(nf90_put_att(ncFileID, varID, "numberofrecords", nrecords))

  call check(nf90_put_att(ncFileID, varID, "units",   def1%unit))
  call check(nf90_put_att(ncFileID, varID, "class",   def1%class))

  select case(OUTtype)
  case(Int1)
    call check(nf90_put_att(ncFileID,varID,"_FillValue",nf90_fill_byte  ))
    call check(nf90_put_att(ncFileID,varID,"scale_factor",scale))
  case(Int2)
    call check(nf90_put_att(ncFileID,varID,"_FillValue",nf90_fill_short ))
    call check(nf90_put_att(ncFileID,varID,"scale_factor",scale))
  case(Int4)
    call check(nf90_put_att(ncFileID,varID,"_FillValue",nf90_fill_int   ))
    call check(nf90_put_att(ncFileID,varID,"scale_factor",scale))
  case(Real4)
    call check(nf90_put_att(ncFileID,varID,"_FillValue",nf90_fill_float ))
  case(Real8)
    call check(nf90_put_att(ncFileID,varID,"_FillValue",nf90_fill_double))
  end select

  call check(nf90_put_att(ncFileID,varID,"current_date_first",ndate))
  call check(nf90_put_att(ncFileID,varID,"current_date_last" ,ndate))

  call check(nf90_enddef(ncFileID))
end subroutine  createnewvariable
!_______________________________________________________________________

subroutine CloseNetCDF
!close open files
!NB the data in a NetCDF file is not "safe" before the file is closed.
!The files are NOT automatically properly closed after end of program,
! and data may be lost if the files are not closed explicitely.
  integer :: i,ncFileID

  if(MasterProc)then
    do i=IOU_INST,IOU_HOUR_INST
      ncFileID=ncFileID_iou(i)
      if(ncFileID/=closedID)then
        call check(nf90_close(ncFileID))
        ncFileID_iou(i)=closedID
      end if
    end do
  end if
end subroutine CloseNetCDF

subroutine GetCDF(varname,fileName,Rvar,varGIMAX,varGJMAX,varKMAX,nstart,nfetch,needed)
! open and reads CDF file
!
! The nf90 are functions which return 0 if no error occur.
! check is only a subroutine which check wether the function returns zero

  use netcdf
  implicit none
  character (len=*),intent(in) :: fileName

  character (len = *),intent(in) ::varname
  integer, intent(in) :: nstart,varGIMAX,varGJMAX,varKMAX
  integer, intent(inout) ::  nfetch
  real, intent(out) :: Rvar(varGIMAX*varGJMAX*varKMAX*nfetch)
  logical, optional,intent(in) :: needed

  logical :: fileneeded
  integer :: status,ndims,alloc_err
  integer :: totsize,xtype,dimids(NF90_MAX_VAR_DIMS),nAtts
  integer :: dims(NF90_MAX_VAR_DIMS),startvec(NF90_MAX_VAR_DIMS)
  integer :: ncFileID,VarID,i
  character*100::name
  real :: scale,offset,scalefactors(2)
  integer, allocatable:: Ivalues(:)

  if(MasterProc.and.DEBUG_NETCDF)print *,'GetCDF  reading ',trim(fileName), ' nstart ', nstart
  !open an existing netcdf dataset
  fileneeded=.true.!default
  if(present(needed)) fileneeded=needed

  if(fileneeded)then
    call check(nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID))
  else
    status=nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID)
    if(status/= nf90_noerr)then
      write(*,*)trim(fileName),' not found (but not needed)'
      nfetch=0
      return
    end if
  end if

  !get global attributes
  !example:
!  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_hour", attribute ))
!  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_date", attribute2 ))
!  print *,'file last modified (yyyymmdd hhmmss.sss) ',attribute2,' ',attribute

  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)

  if(status == nf90_noerr) then
    if(DEBUG_NETCDF)write(*,*) 'variable exists: ',trim(varname)
  else
    print *, 'variable does not exist: ',trim(varname),nf90_strerror(status)
    nfetch=0
    call CheckStop(fileneeded, "NetCDF_mod : variable needed but not found")
    return
  end if

  !get dimensions
  call check(nf90_Inquire_Variable(ncFileID,VarID,name,xtype,ndims,dimids,nAtts))
  dims=0
  totsize=1
  do i=1,ndims
    call check(nf90_inquire_dimension(ncFileID,dimids(i),len=dims(i)))
    totsize=totsize*dims(i)
  ! write(*,*)'size variable ',i,dims(i)
  end do

  if(MasterProc.and.DEBUG_NETCDF)write(*,*)'dimensions ',(dims(i),i=1,ndims)
  if(dims(1)>varGIMAX.or.dims(2)>varGJMAX)then
    write(*,*)'buffer too small',dims(1),varGIMAX,dims(2),varGJMAX
    Call StopAll('GetCDF buffer too small')
  end if

  if(ndims>3.and.dims(3)>varKMAX)then
    if(me==0)write(*,*)'Warning: not reading all levels ',dims(3),varKMAX,trim(varname)
!   Call StopAll('GetCDF not reading all levels')
  end if

  if(nstart+nfetch-1>dims(ndims))then
    write(*,*)'WARNING: did not find all data'
    nfetch=dims(ndims)-nstart+1
    if(nfetch<=0)Call StopAll('GetCDF  nfetch<0')
  end if

  startvec=1
  startvec(ndims)=nstart
  totsize=totsize/dims(ndims)
  if(nfetch<dims(ndims).and.DEBUG_NETCDF)&
    write(*,*)'fetching only',totsize*nfetch,'of ',totsize*dims(ndims),'elements'
  dims(ndims)=nfetch
  totsize=totsize*dims(ndims)

  select case(xtype)
  case(NF90_SHORT,NF90_INT,NF90_BYTE)
    allocate(Ivalues(totsize), stat=alloc_err)
    call check(nf90_get_var(ncFileID, VarID, Ivalues,start=startvec,count=dims))

    scalefactors(1) = 1.0 !default
    scalefactors(2) = 0.  !default
    status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
    if(status == nf90_noerr) scalefactors(1) = scale
    status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
    if(status == nf90_noerr) scalefactors(2) = offset

    do i=1,totsize
      Rvar(i)=Ivalues(i)*scalefactors(1)+scalefactors(2)
    end do

    deallocate(Ivalues)
  case(NF90_FLOAT,NF90_DOUBLE)
    call check(nf90_get_var(ncFileID, VarID, Rvar,start=startvec,count=dims))
    if(DEBUG_NETCDF) &
      write(*,*)'datatype real, read', me, maxval(Rvar), minval(Rvar)
  case default
    write(*,*)'datatype not yet supported'!Char
    Call StopAll('GetCDF  datatype not yet supported')
  end select

  call check(nf90_close(ncFileID))
end subroutine GetCDF

subroutine GetCDF_modelgrid(varname,fileName,Rvar,k_start,k_end,nstart,nfetch,&
                            unit,validity,i_start,j_start,imax_in,jmax_in,reverse_k,needed,found)
! open and reads CDF file 
! The grid MUST be in the same projection and resolution as the model grid
! the field are read directly into the subdomains
! k_start,k_end are the vertical levels to read, k=k_start,k_end
! for instance k_start=KMAX_MID,k_end=KMAX_MID gives only surface
! for instance k_start1,k_end=KMAX_MID gives all levels
! if reverse_k=.true. , the k dimension in the file is reversed
! i_start,j_start can give a new origin instead of (1,1); i_start=2, menas that only i>=2 is read.
! the dimensions of Rvar are assumed:
! Rvar(LIMAX,LJMAX,k_end-k_start+1,nfetch)

  use netcdf
  implicit none
  character(len=*),intent(in) :: fileName,varname
  integer, intent(in) :: nstart,nfetch,k_start,k_end
! real, intent(out) :: Rvar(LIMAX,LJMAX,k_end-k_start+1,nfetch)
  real, intent(out) :: Rvar(*)
  character(len=*), optional, intent(out) ::unit,validity
  integer, optional, intent(in) :: i_start,j_start,imax_in,jmax_in
  logical, optional, intent(in) :: reverse_k,needed
  logical, optional, intent(out):: found

  logical :: fileneeded
  integer :: status,ndims,alloc_err
  integer :: totsize,xtype,dimids(NF90_MAX_VAR_DIMS),nAtts,imax,jmax
  integer :: dims(NF90_MAX_VAR_DIMS),startvec(NF90_MAX_VAR_DIMS)
  integer :: ncFileID,VarID,i,j,k,n,it,i1,j1,i0,j0,ijkn,ijknR,jkn
  character(len=100)::name
  real :: scale,offset
  integer, allocatable:: Ivalues(:)
  real, allocatable:: Rvalues(:)
  logical :: reverse_k_loc

  reverse_k_loc=.false.
  if(present(reverse_k))reverse_k_loc=reverse_k  
  imax=limax
  jmax=ljmax
  if(present(imax_in))imax=imax_in
  if(present(jmax_in))jmax=jmax_in
  
  i0=0;j0=0!origin, i.e. (i=0,j=0) coordinate
  if(present(i_start))i0=i_start-1
  if(present(j_start))j0=j_start-1
  if(i0>imax .or. j0>jmax)then
     if(.not.needed)then
        if(me==0)write(*,*)'WARNING: '//trim(fileName)//'has incompatible grid. Cannot read '//trim(varname)
        found=.false.
        return
     endif
  endif
  call CheckStop(i0>imax, "NetCDF_mod i: subdomain not compatible. cannot handle this")
  call CheckStop(j0>jmax, "NetCDF_mod j: subdomain not compatible. cannot handle this")
  if(MasterProc.and.DEBUG_NETCDF)&
    print *,'GetCDF_modelgrid  reading ',trim(fileName), ' nstart ', nstart

  !open an existing netcdf dataset
  fileneeded=.true. ! default
  if(present(needed))fileneeded=needed
  if(present(found))found=.false.

  if(fileneeded)then
    call check(nf90_open(trim(fileName),nf90_nowrite,ncFileID),&
               trim(fileName)//" not found")
  else
    status=nf90_open(trim(fileName),nf90_nowrite,ncFileID)
    if(status/=nf90_noerr)then
      if(MasterProc) write(*,*)trim(fileName),' not found (but not needed)'
      return
    end if
  end if

  !test if the variable is defined and get varID:
  status=nf90_inq_varid(ncFileID,varname,VarID)
  if(status==nf90_noerr) then
    if(DEBUG_NETCDF)write(*,*) 'variable exists: ',trim(varname)
    if(present(found))found=.true.
  else
    if(MasterProc .and. DEBUG_NETCDF) &
      write(*,*)'variable does not exist: ',trim(varname),trim(nf90_strerror(status))
    call CheckStop(fileneeded, "NetCDF_mod : variable needed but not found "&
         //trim(varname)//' '//trim(fileName))
    call check(nf90_close(ncFileID))
    return
  end if

  !get dimensions
  call check(nf90_Inquire_Variable(ncFileID,VarID,name,xtype,ndims,dimids,nAtts),&
             "inq:"//trim(name))
  dims=0
  do i=1,ndims
    call check(nf90_inquire_dimension(ncFileID,dimids(i),len=dims(i)),"len:dimN")
  end do

  if(DEBUG_NETCDF)write(*,*)'dimensions ',(dims(i),i=1,ndims)
  if(any(dims(1:2)/=[IIFULLDOM,JJFULLDOM]))then
    if(MasterProc.and..not.present(i_start))&
      write(*,"(3A,2I5,A,2I5)")'WARNING : grid ',trim(fileName),&
        ' has different dimensions than model grid: ',&
        dims(1:2),' instead of ',IIFULLDOM,JJFULLDOM
    call CheckStop(.not.present(i_start),'GetCDF_modelgrid: incompatible grid.&
       & This routine does not interpolate. Give i_start, j_start '//trim(fileName))
  end if

  if(ndims>3.and.dims(3)>k_end)then
     if(me==0)write(*,*)'Warning: not reading all levels ',dims(3),k_end,trim(varname)
!    Call StopAll('GetCDF not reading all levels')
  end if

  if(nstart+nfetch-1>dims(ndims))then
    write(*,*)'WARNING: did not find all data'
    call CheckStop(dims(ndims)-nstart+1<=0,'GetCDF_modelgrid  nfetch<0')
  end if

  startvec(:)=1
  startvec(1)=max(1,i_fdom(1)+i0)
  startvec(2)=max(1,j_fdom(1)+j0)
  startvec(3)=k_start
  startvec(ndims)=nstart

  if(startvec(1)+imax-1<dims(1))then
    dims(1)=imax
  else
    dims(1)=dims(1)-startvec(1)+1
  end if
  if(startvec(2)+jmax-1<dims(2))then
    dims(2)=jmax
  else
    dims(2)=dims(2)-startvec(2)+1
  end if
  if(i_fdom(1)==1)dims(1)=min(dims(1),dims(1)+i0)
  if(j_fdom(1)==1)dims(2)=min(dims(2),dims(2)+j0)

  dims(3)=k_end-k_start+1
  if(dims(3)<=0)then
    write(*,*)'WARNING: k_end<k_start. Not fetching anything'
    return
  end if
  dims(ndims)=nfetch!NB: for 2D ndims=3 and we overwrite it!
  83 format(I4,A,3I5,A,15I5)
  if(DEBUG_NETCDF) &
    write(*,83)me,' Reading from ',startvec(1),startvec(2),startvec(3),' to ',&
       startvec(1)+dims(1)-1,startvec(2)+dims(2)-1,startvec(3)+dims(3)-1

  i1=1;j1=1
  if(i_fdom(1)==1)i1=i0+1
  if(j_fdom(1)==1)j1=j0+1
! Rvar=0.0

  if(dims(1)<1 .or. dims(2)<1)then
!     Rvar(1:LIMAX*LJMAX*(k_end-k_start+1)) = 0.0
     call check(nf90_close(ncFileID))
     return
  endif
  totsize=dims(1)*dims(2)*(k_end-k_start+1)*dims(ndims)

  select case(xtype)
  case(NF90_SHORT,NF90_INT,NF90_BYTE)
    ! read scale/offset if present, otherwise assign default value
    if(nf90_get_att(ncFileID,VarID,"scale_factor",scale )/=nf90_noerr)scale=1.0
    if(nf90_get_att(ncFileID,VarID,"add_offset"  ,offset)/=nf90_noerr)offset=0.0

    allocate(Ivalues(totsize),Rvalues(totsize),stat=alloc_err)
    call check(nf90_get_var(ncFileID, VarID, Ivalues,start=startvec,count=dims))
    Rvalues(:)=Ivalues(:)*scale+offset
    deallocate(Ivalues)
  case(NF90_FLOAT,NF90_DOUBLE)
    allocate(Rvalues(totsize), stat=alloc_err)
    call check(nf90_get_var(ncFileID, VarID, Rvalues,start=startvec,count=dims))
!   if(DEBUG_NETCDF) &
!     write(*,*)'datatype real, read', me, maxval(Rvalues), minval(Rvalues)
  case default
    write(*,*)'datatype not yet supported'!Char
    Call CheckStop('GetCDF_modelgrid datatype not yet supported')
  end select

  if(reverse_k_loc)then
    ijkn=0
    do n=1,nfetch
      do k=k_end-k_start+1,1,-1
        do j=1,dims(2)
          jkn=IMAX*(j-1)+imax*jmax*(k-1)+imax*jmax*(k_end-k_start+1)*(n-1)
          do i=1,dims(1)
            ijkn=ijkn+1
            ijknR=i+jkn
            Rvar(ijknR)=Rvalues(ijkn)
          end do
        end do
      end do
    end do
  else
    it=0
    do n=1,max(1,dims(4))
      do k=1,dims(3)
        do j=1,dims(2)
          jkn=IMAX*(j-1)+imax*jmax*(k-1)+imax*jmax*(k_end-k_start+1)*(n-1)
          do i=1,dims(1)
            it=it+1
!           ijknR=ijk+imax*jmax*(k_end-k_start+1)*(n-1)
            ijknR=i+jkn
            Rvar(ijknR)=Rvalues(it)
          end do
        end do
      end do
    end do
  end if

  if(present(unit))then
    status = nf90_get_att(ncFileID, VarID, "units", unit  )
    if(status /= nf90_noerr)then
       unit='unknown' !default
    end if     
  endif
  if(present(validity))then
     status = nf90_get_att(ncFileID, VarID, "validity", validity  )
     if(status == nf90_noerr)then
        !validity  = trim(period_read)
     else
        status = nf90_get_att(ncFileID, VarID, "period_of_validity", &
             validity  )
        if(status /= nf90_noerr)then
           validity='instantaneous' !default
        end if
     end if
  end if

  deallocate(Rvalues)
  call check(nf90_close(ncFileID))
  if(DEBUG_NETCDF)write(*,*)'finished GetCDF_modelgrid', me
end subroutine GetCDF_modelgrid

subroutine WriteCDF(varname,vardate,filename_given,newfile)

!Routine to write data directly into a new file.
!used for testing only

 character (len=*),intent(in)::varname!variable name, or group of variable name
 type(date), intent(in)::vardate!variable name, or group of variable name
 character (len=*),optional, intent(in):: fileName_given!filename to which the data must be written
 logical,optional, intent(in) :: newfile

 real, dimension(LIMAX,LJMAX,KMAX_MID) :: dat ! Data arrays
 character (len=100):: fileName
 real ::scale
 integer :: n,iotyp,ndim,kmax,nseconds
 type(Deriv) :: def1 ! definition of fields

 call date2nctime(vardate,nseconds)
 nseconds=nseconds+vardate%seconds
 write(*,*)nseconds

 iotyp=IOU_INST

 if(present(filename_given))then
    filename=trim(fileName_given)
 else
    filename='EMEP_OUT.nc'
 end if

 if(present(newfile))then
    if(newfile)then
    !make a new file (i.e. delete possible old one)
    if ( me == 0 )then
       write(*,*)'creating',me
       call Init_new_netCDF(fileName,iotyp)
       write(*,*)'created',me
    end if
    end if
 else
    !append if the file exist
 end if

 scale=1.0

 def1%class='Advected' !written
 def1%avg=.false.      !not used
 def1%index=0          !not used
 def1%scale=scale      !not used
! def1%inst=.true.      !not used
! def1%year=.false.     !not used
! def1%month=.false.    !not used
! def1%day=.false.      !not used
 def1%name='O3'        !written
 def1%unit='mix_ratio'       !written


 if(trim(varname)=='ALL')then
 ndim=3 !3-dimensional
 kmax=KMAX_MID

 if(NSPEC_SHL+ NSPEC_ADV /=  NSPEC_TOT.and. MasterProc)then
    write(*,*)'WARNING: NSPEC_SHL+ NSPEC_ADV /=  NSPEC_TOT'
    write(*,*) NSPEC_SHL,NSPEC_ADV, NSPEC_TOT
    write(*,*)'WRITING ONLY SHL and ADV'
    write(*,*)'Check species names'
 end if

! def1%class='Short_lived' !written
! do n=1, NSPEC_SHL
! def1%name= species(n)%name       !written
! dat=xn_shl(n,:,:,:)
! call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real8,fileName_given=fileName)
! end do

 def1%class='Advected' !written
 do n= 1, NSPEC_ADV
 def1%name= species(NSPEC_SHL+n)%name       !written
 dat=xn_adv(n,:,:,:)
 call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real4,&
                 out_DOMAIN=[60,107,11,58],fileName_given=fileName)
                !out_DOMAIN=[ist,ien,jst,jen]
 end do

  elseif(trim(varname)=='LIST')then
     ndim=3 !3-dimensional
     kmax=KMAX_MID


     def1%class='Advected' !written
     do n= 1, NSPEC_ADV
        def1%name= species(NSPEC_SHL+n)%name       !written
        if(trim(def1%name)=='O3'.or.trim(def1%name)=='NO2')then
           dat=xn_adv(n,:,:,:)
           call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real4,&
                           out_DOMAIN=[10,20,10,20],fileName_given=fileName)
                          !out_DOMAIN=[ist,ien,jst,jen]
        end if
     end do

else

    if(MasterProc)write(*,*)'case not implemented'
 end if


end subroutine WriteCDF

subroutine ReadField_CDF(fileName,varname,Rvar,nstart,kstart,kend,interpol, &
     known_projection,  &! can be provided by user, eg. for MEGAN.
     use_lat_name, use_lon_name,stagg, & !for telling the routine that the field is defined in a staggered grid
     Grid_resolution_in, & !Better to put it as attribute, but can also be given here
     fractions_out,CC_out,Ncc_out,Reduc,&! additional output for emissions given with country-codes
     Mask_fileName,Mask_varname,Mask_Code,NMask_Code,Mask_ReducFactor,&
     needed,found,unit,validity,debug_flag,UnDef,ncFileID_given)
  !
  !reads data from netcdf file and interpolates data into model local (subdomain) grid
  !under development!
  !
  !dimensions and indices:
  !Rvar is assumed to have declared dimensions LIMAX,LJMAX in 2D.
  !If 3D, k coordinate in Rvar assumed as first coordinate. Could consider to change this.
  !nstart (optional) is the index of last dimension in the netcdf file, generally time dimension.
  !
  !undefined field values:
  !Some data can be missing/not defined for some gridpoints;
  !If general projection is used (not lon lat or polarstereo), it takes the nearest value.
  !If Undef is present, it is used as value for undefined gridpoints;
  !If it is not present, an error occurs if data is missing.
  !Data can be undefined either because it is outside the domain in the netcdf file,
  !or because it has the value defined in "FillValue" ("FillValue" defined from netcdf file)
  !

  !projections:
  !Lon-lat projection in the netcdf file is implemented with most functionalities.
  !General projection in the netcdf file is still primitive. Limitations are: no 3D,
  !   no Undef, only linear interpolation, cpu expensive.
  !The netcdf file projection is defined by user in "known_projection" or read from
  !netcdf file (in attribute "projection").
  !If the model grid projection is not lon-lat and not stereographic the method is not
  !very CPU efficient in the present version
  !(except conservative inerpolation, where only nearest neighbor is implemented).
  !Vertical interpolation is not implemented, except from "Fligh Levels", but

  !"Flight Levels" are so specific that they are move in an own routine ReadField_CDF_FL

  !interpolation:
  !'zero_order' gives value at closest gridcell. Probably good enough for most applications.
  !Does not smooth out values
  !'conservative' and 'mass_conservative' give smoother fields and are approximatively
  !integral conservative (integral over a region is conserved). The initial gridcells
  !are subdivided into smaller subcells and each subcell is assigned to a cell in the model grid
  !'conservative' can be used for emissions given in kg/m2 (or kg/m2/s) or landuse or most fields.
  !The value in the netcdf file and in model gridcell are of the similar.
  !'mass_conservative' can be used for emissions in kg (or kg/s). If the gricell in the model are
  !twice as small as the gridcell in the netcdf file, the values will also be reduced by a factor 2.

  !Emissions with country-codes: (July 2012, under development)
  !Emissions are given in each gridcell with:
  !1) Total (Rvar)
  !2) Number of country-codes (Ncc_out)
  !3) Country codes (CC_out)
  !4) fraction assigned to each country (fractions_out)
  !Presently only lat-lon projection of input file supported
  !negative data not finished/tested (can give 0 totals, definition of fractions?)

  !MASK:
  !set up to reduce emissions in the native grid. The file MUST be the same size as the emission grid 
  !the routine will however be able to reverse the latitude direction, if necessary.
  !the reduction can applied on top of the reduction from "Reduc"
  !Mask_fileName -> name of the file
  !Mask_varname -> name of the variable used for mask ,
  !Mask_Code(1:NMask_Code) -> values from the mask_value for which to apply 
  !    the reduction (reduction is applied if any of the Mask_Code=Mask_value)
  !NMask_Code -> size of the array Mask_Code with defined values
  !Mask_ReducFactor -> multiplicative factor for reduction where mask applies

  !ncFileID_given , if present gives the ncFileID value AND assumes the file 
  !is already open. It is not closed on exit either. 
  !Much faster if little to do.

  !Technical, future developements:
  !This routine is likely to change a lot in the future: should be divided into simpler routines;
  !more functionalities will be introduced.
  !Should also be usable as standalone.
  !All MPI processes read the same file simultaneously (i.e. in parallel).
  !They read only the chunk of data they need for themselves.

  use netcdf

  implicit none
  character(len=*), parameter :: dtxt = 'ReadF:' ! debug text
  character(len=99) :: dtxt_msg
  character(len = *),intent(in) ::fileName,varname
  real,intent(out) :: Rvar(*)
  integer,intent(in) :: nstart
  character(len = *), optional,intent(in) :: interpol
  character(len = *), optional,intent(in) :: known_projection
  character(len = *), optional,intent(in) :: use_lat_name, use_lon_name
  character(len = *), optional,intent(in) :: stagg
  real, optional, intent(in) :: Grid_resolution_in !(approximative) resolution of the data, steers Ndiv
  logical, optional, intent(in) :: needed
  logical, optional, intent(out) :: found
  character(len=*), optional,intent(out) ::unit,validity
  integer, optional,intent(in) :: kstart!smallest k (vertical level) to read. Default: assume 2D field
  integer, optional,intent(in) :: kend!largest k to read. Default: assume 2D field
  logical, optional, intent(in) :: debug_flag
  real, optional, intent(in) :: UnDef ! Value put into the undefined gridcells
  real , optional, intent(out) ::fractions_out(LIMAX*LJMAX,*) !fraction assigned to each country 
  integer, optional, intent(out)  ::Ncc_out(*), CC_out(LIMAX*LJMAX,*) !Number of country-codes and Country codes
  real, optional, intent(in) :: Reduc(NLAND)
  character(len = *), optional,intent(in) :: Mask_filename,Mask_varname
  integer , optional, intent(in) :: Mask_Code(*),NMask_Code! reduce only where mask take this value. If not defined, multiply by mask
  real, optional, intent(in) :: Mask_ReducFactor!multiply by this. default 0, i.e. reduce by 100% according to Mask
  integer , optional, intent(in) :: ncFileID_given

  real  :: factor

  logical, save :: debug_ij
  logical ::fractions,interpolate_vertical
  integer :: ncFileID,VarID,lonVarID,latVarID,status,xtype,xtype_lon,xtype_lat,ndims,dimids(NF90_MAX_VAR_DIMS),nAtts
  integer :: VarIDCC,VarIDNCC,VarIDfrac,NdimID, ndimslon
  integer :: ncFileID_Mask, VarID_Mask,ndims_Mask,dimids_Mask(NF90_MAX_VAR_DIMS),xtype_Mask,isize,jsize
  integer :: dims(NF90_MAX_VAR_DIMS),NCdims(NF90_MAX_VAR_DIMS),totsize,i,j,k
  integer :: startvec(NF90_MAX_VAR_DIMS),Nstartvec(NF90_MAX_VAR_DIMS)
  integer ::alloc_err
  character*100 ::name,used_lat_name, used_lon_name
  real :: scale,offset,scalefactors(2),dloni,dlati,dlonx,dlatx,dlony,dlaty,dlon,dlat
  integer ::ij,jdiv,idiv,Ndiv_lon,Ndiv,Ndiv2,igjgk,ig,jg,ijk,n,im,jm,ijm,iw
  integer ::imin,imax,jmin,jjmin,jmax,igjg,k2
  integer, allocatable:: Ivalues(:)  ! I counts all data
  integer, allocatable:: Nvalues(:)  !ds counts all values
  real, allocatable:: Rvalues(:),Rlon(:),Rlat(:),lon_mask(:),lat_mask(:)
  integer, allocatable:: Mask_values(:)
  real ::lat,lon,maxlon,minlon,maxlat,minlat,maxlon_var,minlon_var,maxlat_var,minlat_var
  logical ::fileneeded, debug,data3D
  character(len = TXTLEN_NAME) :: interpol_used, data_projection=""
  real :: Grid_resolution_lon,Grid_resolution
  type(Deriv) :: def1 ! definition of fields
  logical ::  OnlyDefinedValues

  real, pointer :: Weight(:,:,:)=>null(), ww(:)=>null()
  integer, allocatable :: IIij(:,:,:),JJij(:,:,:)
  real :: FillValue=0
  real    :: sumWeights
  integer, dimension(4) :: ijn
  integer :: ii, jj,i_ext,j_ext
  real::an_ext,xp_ext,yp_ext,fi_ext,ref_lat_ext,xp_ext_div,yp_ext_div,Grid_resolution_div,an_ext_div
  real ::buffer1(LIMAX, LJMAX),buffer2(LIMAX, LJMAX)
  real, allocatable ::fraction_in(:,:)
  integer, allocatable ::CC(:,:),Ncc(:)
  real :: UnDef_local
  integer :: Nmax,kstart_loc,kend_loc,lon_shift_Mask,startlat_Mask
  logical :: Reverse_lat_direction_Mask,ilast_corrected
  character(len=*),parameter  :: field_not_found='field_not_found'
  real :: latlon_weight

  real ::Rlonmin,Rlonmax,dRlon,dRloni,frac,frac_j,ir,jr
  real ::Rlatmin,Rlatmax,dRlat,dRlati, Resolution_fac
  integer, allocatable ::ifirst(:),ilast(:),jfirst(:),jlast(:)
  real, allocatable :: fracfirstlon(:),fraclastlon(:),fracfirstlat(:),fraclastlat(:)
  logical, save :: debug1   ! -> output on 1st step, or when debug set

  !_______________________________________________________________________________
  !
  !1)           General checks and init
  !_______________________________________________________________________________

  fileneeded=.true.!default

  debug = .false.
  dtxt_msg = dtxt//trim(fileName)//':'//trim(varname)//':' ! Message for errors
  if(present(debug_flag))then
     debug = debug_flag .and. me==0
     if ( debug ) write(*,*) 'ReadCDF start: ',trim(filename),':', trim(varname)
  end if

  debug1 = ( (MasterProc .and. step_main == 1) .or. debug )

  UnDef_local=0.0
  if(present(UnDef))UnDef_local=UnDef
  if ( debug ) write(*,*) 'ReadCDF UnDef defined as: ',UnDef_local

  data3D=.false.
  kstart_loc=0
  kend_loc=0
  if(present(kstart).or.present(kend))then
     call CheckStop((.not. present(kstart).or. .not. present(kend)), &
          "ReadField_CDF : both or none kstart and kend should be present")
     kstart_loc=kstart
     kend_loc=kend
     data3D=.true.
  end if

  if(present(needed))   fileneeded=needed
  if(present(found))found=.true.!default

  if(present(ncFileID_given))then
     ncFileID = ncFileID_given
  else
     !open an existing netcdf dataset
     status=nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID)

     if(status == nf90_noerr) then
        if ( debug ) write(*,*) 'ReadCDF reading ',trim(filename), ' nstart ', nstart
     else
        if(fileneeded)then
           print *, 'file does not exist: ',trim(fileName),nf90_strerror(status)
           call CheckStop(fileneeded, "ReadField_CDF : file needed but not found")
        else
           if(MasterProc) write(*,*)'file does not exist (but not needed): ',&
                trim(fileName),nf90_strerror(status)
           if(present(found))found=.false.
           if(present(validity))validity=field_not_found
           if(present(unit))unit=' '
           k2=1
           if(data3D)k2=kend_loc-kstart_loc+1
           Rvar(1:LIMAX*LJMAX*k2)=UnDef_local!works only for 2D or when kstart is defined!
           return
        end if
     end if
  end if

  interpol_used='zero_order'!default
  if(present(interpol))then
     interpol_used=interpol
  end if
  call CheckStop(interpol_used/='zero_order'.and.&
       interpol_used/='conservative'.and.&
       interpol_used/='mass_conservative',&
       'interpolation method not recognized')


  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = trim(varname), varID = VarID)
  if(status == nf90_noerr) then
    if(debug) write(*,*) 'ReadCDF variable exists: ',trim(varname)
  else
    !     nfetch=0
    if(fileneeded)then
      print *, 'variable does not exist: ',trim(varname),nf90_strerror(status)
      call CheckStop(fileneeded, "ReadField_CDF : variable needed but not found")
!     if(MasterProc)&
!       write(*,*)'variable does not exist: ',trim(varname),' set to ',UnDef_local
!     Rvar(1:LIMAX*LJMAX)=UnDef_local
!     return
    else
        if(MasterProc)write(*,*) 'variable does not exist (but not needed): ',&
             trim(varname),nf90_strerror(status)
        if(.not. present(ncFileID_given))call check(nf90_close(ncFileID))
        if(present(found))found=.false.
        if(present(validity))validity=field_not_found
        if(present(unit))unit=' '
        k2=1
        if(data3D)k2=kend_loc-kstart_loc+1
        Rvar(1:LIMAX*LJMAX*k2)=UnDef_local!works only for 2D or when kstart is defined!
        return
    end if
  end if
  if(present(unit))then
    !find unit
    unit=' '
    status = nf90_get_att(ncFileID, VarID, "units", unit  )
    if(status /= nf90_noerr)then
       unit='unknown' !default
    end if
  end if
  if(present(validity))then
!     status = nf90_get_att(ncFileID, VarID, "validity",  validity)
     if(status /= nf90_noerr)then
        status = nf90_get_att(ncFileID, VarID, "period_of_validity",validity)
        if(status /= nf90_noerr)then
           validity='instantaneous' !default
        end if
     end if
  end if

  fractions=.false.
  if(present(fractions_out))fractions=.true.


  !Check first that variable has data covering the relevant part of the grid:

  !Find chunk of data required (local)
  if(projection=="lon lat")then
     maxlon=maxval(glon)+0.5*(glon(2,1)-glon(1,1))
     minlon=minval(glon)-0.5*(glon(2,1)-glon(1,1))
     maxlat=maxval(glat)+0.5*(glat(1,2)-glat(1,1))
     minlat=minval(glat)-0.5*(glat(1,2)-glat(1,1))
  else
     maxlon=max(maxval(gl_stagg),maxval(glon))
     minlon=min(minval(gl_stagg),minval(glon))
     maxlat=max(maxval(gb_stagg),maxval(glat))
     minlat=min(minval(gb_stagg),minval(glat))
  endif
  !Read the extension of the data in the file (if available)
  status = nf90_get_att(ncFileID, VarID, "minlat", minlat_var  )
  if(status == nf90_noerr)then
     !found minlat, therfore the other (maxlat,minlon,maxlat) expected too
     if ( debug ) write(*,*) 'minlat attribute found: ',minlat_var
     call CheckStop(fractions, &
          "ReadField_CDF: minlat not implemented for fractions")     
     call CheckStop(kstart_loc>0, &
          "ReadField_CDF: minlat not implemented for vertical interpolation")     
     k2=1
     if(data3D)k2=kend_loc-kstart_loc+1
     ijk=LIMAX*LJMAX*k2
     if(minlat_var>maxlat)then
        !the data is outside range. put zero or Undef.
        Rvar(1:ijk)=UnDef_local
        if ( debug ) write(*,*) 'data out of maxlat range ',maxlat
        if(.not. present(ncFileID_given))call check(nf90_close(ncFileID))
        return
     end if
     status = nf90_get_att(ncFileID, VarID, "maxlat", maxlat_var  )
     if(status == nf90_noerr)then
        if ( debug ) write(*,*) 'maxlat attribute found: ',maxlat_var
        if(maxlat_var<minlat)then
           !the data is outside range. put zero or Undef.
           Rvar(1:ijk)=UnDef_local
           if ( debug ) write(*,*) 'data out of minlat range ',minlat
           if(.not. present(ncFileID_given))call check(nf90_close(ncFileID))
           return
        end if
     end if
     status = nf90_get_att(ncFileID, VarID, "minlon", minlon_var  )
     if(status == nf90_noerr)then
        if ( debug ) write(*,*) 'minlon attribute found: ',minlon_var
        if(minlon_var>maxlon)then
           !the data is outside range. put zero or Undef.
           Rvar(1:ijk)=UnDef_local
           if ( debug ) write(*,*) 'data out of minlon range ',minlon
           if(.not. present(ncFileID_given))call check(nf90_close(ncFileID))
           return
        end if
     end if
     status = nf90_get_att(ncFileID, VarID, "maxlon", maxlon_var  )
     if(status == nf90_noerr)then
        if ( debug ) write(*,*) 'maxlon attribute found: ',maxlon_var
        if(maxlon_var<minlon)then
           !the data is outside range. put zero or Undef.
           Rvar(1:ijk)=UnDef_local
           if ( debug ) write(*,*) 'data out of maxlon range ',maxlon
           if(.not. present(ncFileID_given))call check(nf90_close(ncFileID))
           return
        end if
     end if
  else
     !dont expect to find maxlat,minlon or maxlat, therfore don't check
     if ( debug ) write(*,*) 'minlat attribute not found for ',trim(varname)
  end if

  !get dimensions id
  call check(nf90_Inquire_Variable(ncFileID,VarID,name,&
       xtype,ndims,dimids,nAtts),"GetDimsId")

  !only characters cannot be handled
  call CheckStop(xtype==NF90_CHAR,"ReadField_CDF: Datatype not recognised")

  !Find whether Fill values are defined
  status=nf90_get_att(ncFileID, VarID, "_FillValue", FillValue)
  OnlyDefinedValues=.not.present(UnDef)
  if(status == nf90_noerr)then
    OnlyDefinedValues=.false.
    if( debug ) write(*,*)' FillValue (not counted)',FillValue
  else ! ds Jun 2010
    FillValue = -1.0e35  ! Should not arise in normal data
    status=nf90_get_att(ncFileID, VarID, "missing_value", FillValue)
    if(status == nf90_noerr)then
      OnlyDefinedValues=.false.
      if ( debug ) write(*,*)' FillValue found from missing_value (not counted)',FillValue
    else
      select case(xtype)
        case(NF90_BYTE  );FillValue=NF90_FILL_BYTE
        case(NF90_INT   );FillValue=NF90_FILL_INT
        case(NF90_FLOAT );FillValue=NF90_FILL_FLOAT
      ! case(NF90_REAL  );FillValue=NF90_FILL_REAL ! same as FLOAT
        case(NF90_DOUBLE);FillValue=NF90_FILL_DOUBLE
      end select
      if( debug ) write(*,*) 'FillValue not found, using ',FillValue
    end if
  end if

  !get dimensions
  startvec=1
  dims=0
  do i=1,ndims
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(i), &
          len=dims(i)),"GetDims")
     if ( debug ) write(*,*) 'ReadCDF size variable ',i,dims(i)
  end do

  xtype_lon=NF90_FLOAT !default
  xtype_lat=NF90_FLOAT !default
  if( present(known_projection) ) then
     data_projection = trim(known_projection)
     if(trim(known_projection)=="longitude latitude")data_projection = "lon lat"
     if ( debug ) write(*,*) 'data known_projection ',trim(data_projection)
  else
     status = nf90_get_att(ncFileID, nf90_global, &
                           "projection", data_projection )
     if ( debug ) write(*,*) 'data projection from file ',trim(data_projection)
     if(status == nf90_noerr)then
    !remove invisible char(0)!!
        if(trim(data_projection(1:7))=="lon lat")data_projection="lon lat"!remove invisible char(0)!!
        if ( debug ) write(*,*) 'data projection from file ',trim(data_projection)
     else
        status=nf90_inq_varid(ncid = ncFileID, name='lon', varID = lonVarID)
        if(status /= nf90_noerr) then
           status=nf90_inq_varid(ncid = ncFileID, name = 'LON', varID = lonVarID)
           if(status /= nf90_noerr) then
              status=nf90_inq_varid(ncid = ncFileID, name = 'longitude', varID = lonVarID)
              call CheckStop(status /= nf90_noerr,'did not find projection nor longitude variable '//trim(fileName))
           end if
        end if
        call check(nf90_Inquire_Variable(ncid = ncFileID,  varID = lonVarID,xtype=xtype_lon, ndims=ndimslon ))
        if(ndimslon==1)then
           if(MasterProc)write(*,*)'Assuming lon lat projection, because lon is one-dimensional'
           data_projection = "lon lat"
        else
           if(MasterProc)write(*,*)'Warning: did not find projection for '//trim(fileName) 
           data_projection = "Unknown"
        endif
     endif
  end if
!  if(MasterProc)write(*,*)'Interpolating ',trim(varname),' from ',trim(filename),' to present grid'

  if(trim(data_projection)=="lon lat")then 
     ! here we have simple 1-D lat, lon
     allocate(Rlon(dims(1)), stat=alloc_err)
     allocate(Rlat(dims(2)), stat=alloc_err)
     if ( debug ) write(*,"(a,a,i5,i5,a,i5)") 'alloc_flag lon lat ',&
          trim(data_projection), alloc_err, dims(1), "x", dims(2)
  else
     allocate(Rlon(dims(1)*dims(2)), stat=alloc_err)
     allocate(Rlat(dims(1)*dims(2)), stat=alloc_err)
     if ( debug ) write(*,*) 'data allocElse ',trim(data_projection), alloc_err, trim(fileName)
  end if

  call check_lon_lat(ncFileID, used_lon_name, used_lat_name, n, lonVarID, latVarID)
  if(present(use_lat_name)) used_lat_name=trim(use_lat_name)
  if(present(use_lon_name)) used_lon_name=trim(use_lon_name)
  if ( debug ) write(*,*) 'ReadCDF using lon name',trim(used_lon_name)
  if ( debug ) write(*,*) 'ReadCDF using lat name',trim(used_lat_name)

  if(trim(data_projection)=="lon lat")then
     call check(nf90_get_var(ncFileID, lonVarID, Rlon), 'Getting Rlon')
     call check(nf90_get_var(ncFileID, latVarID, Rlat), 'Getting Rlat')
  else
     call check(nf90_get_var(ncFileID, lonVarID, Rlon,start=(/1,1/),count=(/dims(1),dims(2)/)))
     call check(nf90_get_var(ncFileID, latVarID, Rlat,start=(/1,1/),count=(/dims(1),dims(2)/)))
  end if
  if(xtype_lon==NF90_INT.or.xtype_lon==NF90_SHORT.or.xtype_lon==NF90_BYTE)then
     !scale data if it is packed
     scalefactors(1) = 1.0 !default
     scalefactors(2) = 0.  !default
     status = nf90_get_att(ncFileID, lonVarID, "scale_factor", scale  )
     if(status == nf90_noerr) scalefactors(1) = scale
     status = nf90_get_att(ncFileID, lonVarID, "add_offset",  offset )
     if(status == nf90_noerr) scalefactors(2) = offset
     Rlon=Rlon*scalefactors(1)+scalefactors(2)
  end if
  if(xtype_lat==NF90_INT.or.xtype_lat==NF90_SHORT.or.xtype_lat==NF90_BYTE)then
     !scale data if it is packed
     scalefactors(1) = 1.0 !default
     scalefactors(2) = 0.  !default
     status = nf90_get_att(ncFileID, latVarID, "scale_factor", scale  )
     if(status == nf90_noerr) scalefactors(1) = scale
     status = nf90_get_att(ncFileID, latVarID, "add_offset",  offset )
     if(status == nf90_noerr) scalefactors(2) = offset
     Rlat=Rlat*scalefactors(1)+scalefactors(2)
  end if
  if(present(stagg))then
     !use staggered grid
     if(stagg=='stagg_u')then
        if(trim(data_projection)=="lon lat")then
           do i=1,dims(1)-1
              Rlon(i)=(Rlon(i)+Rlon(i+1))*0.5
           end do
           Rlon(dims(1))=Rlon(dims(1))+(Rlon(dims(1))-Rlon(dims(1)-1))*0.5
        else
!not tested
           do j=1,dims(2)
           do i=1,dims(1)-1
              ij=i+(j-1)*dims(1)
              Rlon(ij)=(Rlon(ij)+Rlon(ij+1))*0.5
              Rlat(ij)=(Rlat(ij)+Rlat(ij+1))*0.5
           end do
              ij=dims(1)+(j-1)*dims(1)
           Rlon(ij)=Rlon(ij)+(Rlon(ij)-Rlon(ij-1))*0.5
           Rlat(ij)=Rlat(ij)+(Rlat(ij)-Rlat(ij-1))*0.5
           end do
        end if
     elseif(stagg=='stagg_v')then
        if(trim(data_projection)=="lon lat")then
           do j=1,dims(2)-1
              Rlat(j)=(Rlat(j)+Rlat(j+1))*0.5
           end do
           Rlat(dims(2))=Rlat(dims(2))+(Rlat(dims(2))-Rlat(dims(2)-1))*0.5
        else
!not tested
           do j=1,dims(2)-1
           do i=1,dims(1)
              ij=i+(j-1)*dims(1)
              Rlon(ij)=(Rlon(ij)+Rlon(ij+dims(1)))*0.5
              Rlat(ij)=(Rlat(ij)+Rlat(ij+dims(1)))*0.5
           end do
           end do
           do i=1,dims(1)
              ij=i+(dims(2)-1)*dims(1)
              Rlon(ij)=Rlon(ij)+(Rlon(ij)-Rlon(ij-dims(1)))*0.5
              Rlat(ij)=Rlat(ij)+(Rlat(ij)-Rlat(ij-dims(1)))*0.5
           end do
        end if
     else
        call StopAll("ReadField_CDF: stagg not recognized")
     end if
  end if
!Should define longitude in the range [-180, 180]
  do ij=1,size(Rlon)
     if(Rlon(ij)>180)Rlon(ij)=Rlon(ij)-360.0
     if(Rlon(ij)<-180)Rlon(ij)=Rlon(ij)+360.0
  end do

  if ( debug ) write(*,*) 'ReadCDF lon bounds',minval(Rlon),maxval(Rlon)
  if ( debug ) write(*,*) 'ReadCDF lat bounds',minval(Rlat),maxval(Rlat)

  call CheckStop(fractions.and.trim(data_projection)/="lon lat", &
       "ReadField_CDF: only implemented lon lat projection for fractions")     

  !if Mask_filename is given, open that file
  if(present(Mask_fileName))then
     if ( debug ) write(*,*) 'ReadCDF reading also ',trim(Mask_fileName)
     N=1
     if(present(NMask_Code))N=NMask_Code
     if(N>0)then!otherwise no need to do anything!
     status=nf90_open(path = trim(Mask_fileName), mode = nf90_nowrite, ncid = ncFileID_Mask)
     if(status /= nf90_noerr) then
        write(*,*)'MASK file does not exist: ',trim(Mask_fileName),nf90_strerror(status)
        write(*,*)'Check the call for monthly emis in Emissions_mod.f90'
        call StopAll("ReadField_CDF : Mask file needed but not found")
     endif
    !verify that x, y dimensions have same size
     call check(nf90_inq_varid(ncid = ncFileID_Mask, name = trim(Mask_varname), varID = VarID_Mask),"Var_Mask")
     call check(nf90_Inquire_Variable(ncFileID_Mask,VarID_Mask,name,xtype_Mask,ndims_Mask,dimids_Mask,nAtts),"GetDimsId_Mask")
     call check(nf90_inquire_dimension(ncid=ncFileID_Mask, dimID=dimids_Mask(1),  len=isize))
     call check(nf90_inquire_dimension(ncid=ncFileID_Mask, dimID=dimids_Mask(2),  len=jsize))
     if(isize/=dims(1).or.jsize/=dims(2))then
        write(*,*)'Mask file and Variable files have different dimensions'
        write(*,*)'Sizes '//trim(fileName),dims(1),dims(2)
        write(*,*)'Sizes '//trim(Mask_fileName),isize,jsize
        call StopAll("ReadField_CDF: Incompatible file sizes")
     end if
     if(ndims_Mask/=2)then
        call StopAll("ReadField_CDF: only 2 dimensional mask implemented")
     end if
     else
        if (debug)write(*,*) 'ReadCDF MASK: no need to reduce anything for ', trim(varname)
     end if
  end if

  !_______________________________________________________________________________
  !
  !2)        Coordinates conversion and interpolation
  !_______________________________________________________________________________


  if(trim(data_projection)=="lon lat")then

     !get coordinates
     !we assume first that data is originally in lon lat grid
     !check that there are dimensions called lon and lat

     call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(1), name=name ),name)
     call CheckStop(trim(name)/='lon'.and.trim(name)/='longitude',"longitude not found")
     call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(2), name=name ),name)
     call CheckStop(trim(name)/='lat'.and.trim(name)/='latitude',"latitude not found")

     if(data3D)then
        call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(3), name=name ))
        call CheckStop(name/='k'.and.name/='N'.and.name/='lev'.and.name/='height',&
          "vertical coordinate (k, lev, N or height) not found")
     end if

     !NB: we assume regular grid
     !inverse of resolution
     dRloni=1.0/(Rlon(2)-Rlon(1))
     dRlati=1.0/(Rlat(2)-Rlat(1))
     
     !not likely to happen, but check. Quick fix: change number of procs
     call CheckStop(abs(Rlon(2)-Rlon(1))>180.0,&
          "longitude error: crossing of -180 not implemented. file "//trim(fileName))
     
     
     Grid_resolution = EARTH_RADIUS*abs(Rlat(2)-Rlat(1))*PI/180.0
     if(present(Grid_resolution_in))Grid_resolution = Grid_resolution_in
     Grid_resolution_lon = EARTH_RADIUS*abs(Rlon(2)-Rlon(1))*PI/180.0!NB: varies with latitude
     Resolution_fac = max(1.0,Grid_resolution/GRIDWIDTH_M)
     !the method chosen depends on the relative resolutions
     if(.not.present(interpol).and.Grid_resolution/GRIDWIDTH_M>4)then
        interpol_used='zero_order'!usually good enough, and keeps gradients
     end if
     if ( debug ) write(*,*) 'interpol_used: ',interpol_used

     if(debug) then
        write(*,*) "SET Grid resolution:" // trim(fileName), Grid_resolution
        write(*,"(a,6f8.2,2x,4f8.2)") 'ReadCDF LL values ',&
             Rlon(2),Rlon(1),dRloni, Rlat(2),Rlat(1), dRlati, &
             maxlon, minlon, maxlat, minlat
     end if

     if(mod(nint(100.*(Rlon(dims(1))+(Rlon(2)-Rlon(1))-Rlon(1))),36000)/=0)then
           if(MasterProc.and.debug)write(*,*)'Grid does not cover 360 degrees ', Rlon(1),Rlon(dims(1))
        !the grid does not cover 360 degrees
        if(Rlon(1)>Rlon(dims(1)))then
           if(MasterProc.and.debug)write(*,*)'Longitudes not monotonic ', Rlon(1),Rlon(dims(1))
           imin=1!cover everything available
           imax=dims(1)!cover everything available
        else
          imin=max(1,floor((minlon-Rlon(1))*dRloni-0.001))+1
          imax=min(dims(1),ceiling((maxlon-Rlon(1))*dRloni+1.001))
        endif
     else if(dRloni>0)then
        !floor((minlon-Rlon(1))*dRloni)<=number of gridcells between minlon and Rlon(1)
        !mod(floor((minlon-Rlon(1))*dRloni)+dims(1),dims(1))+1 = get a number in [1,dims(1)]
        imin=mod( floor((minlon-Rlon(1))*dRloni)+dims(1),dims(1))+1!NB lon  -90 = +270
        imax=mod(ceiling((maxlon-Rlon(1))*dRloni)+dims(1),dims(1))+1!NB lon  -90 = +270
        if(imax==1)imax=dims(1)!covered entire circle
        if(minlon-Rlon(1)<0.0)then
           imin=1!cover entire circle
           imax=dims(1)!cover entire circle
        end if
     else
        call CheckStop("Not tested: negativ dRloni")
        imin=mod(floor((maxlon-Rlon(1))*dRloni)+dims(1),dims(1))+1!NB lon  -90 = +270
        imax=mod(ceiling((minlon-Rlon(1))*dRloni)+dims(1),dims(1))+1!NB lon  -90 = +270
        if(imax==1)imax=dims(1)!covered entire circle
     end if

     if(dRlati>0)then
        jmin=max(1,min(dims(2),floor((minlat-Rlat(1))*dRlati)))
        jmax=max(1,min(dims(2),ceiling((maxlat-Rlat(1))*dRlati)+1))
     else!if starting to count from north pole
        jmin=max(1,min(dims(2),floor((maxlat-Rlat(1))*dRlati)))!maxlat is closest to Rlat(1)
        jmax=max(1,min(dims(2),ceiling((minlat-Rlat(1))*dRlati)+1))
     end if

     if(maxlat>85.0.or.minlat<-85.0  .and. &
       ((.not.(projection=='lon lat')) .or. (.not.data_projection=='lon lat')))then
        !close to poles
        imin=1
        imax=dims(1)
     end if

     !latitude is sometime counted from north pole, sometimes from southpole:
     jjmin=jmin
     jmin=min(jmin,jmax)
     jmax=max(jjmin,jmax)
     if(imax<imin)then
        !crossing longitude border !
        !   write(*,*)'WARNING: crossing end of map'
        !take everything...could be memory expensive
        imin=1
        imax=dims(1)
     end if

     if ( debug ) write(*,"(a,4f8.2,6i8)") 'ReadCDF minmax values ',&
          minlon,maxlon,minlat,maxlat,imin,imax,jmin,jmax

     startvec(1)=imin
     startvec(2)=jmin
     if(ndims>2)startvec(ndims)=nstart
     dims=1
     dims(1)=imax-imin+1
     dims(2)=jmax-jmin+1

     if(data3D)then
        if(kstart_loc<0)then
           !take all levels
           kstart_loc=1
           kend_loc=dims(3)
!NB: size of Rvar is not computable (by ifort)
!           call CheckStop(size(Rvar)<LIMAX*LJMAX*(kend_loc-kstart_loc+1), &
!          "ReadField_CDF: size of Rvar is too small")
           interpolate_vertical=.true.
        end if
        startvec(3)=kstart_loc
        dims(3)=kend_loc-kstart_loc+1
        if(ndims>3)startvec(ndims)=nstart
     end if

     totsize=1
     do i=1,ndims
        totsize=totsize*dims(i)
     end do

     allocate(Rvalues(totsize), stat=alloc_err)
     if ( debug ) then
        do i=1, ndims ! NF90_MAX_VAR_DIMS would be 1024
           write(*,"(a,6i8)") 'ReadCDF reading chunk ',i ,startvec(i), dims(i)
        end do
     end if
     call check(nf90_get_var(ncFileID, VarID, Rvalues,start=startvec,count=dims),&
          errmsg=dtxt_msg//"RRvalues")

     if(present(Mask_filename))then
        N=1
        if(present(NMask_Code))N=NMask_Code
        if(N>0)then!otherwise no need to do anything!
      ! assume that mask file has exactely the same size
        allocate(Mask_values(totsize))
        allocate(lon_mask(isize))
        allocate(lat_mask(jsize))
        ! find origin and dimensions
        Reverse_lat_direction_Mask=.false.
        status=nf90_get_var(ncFileID_Mask, dimids_Mask(2), lat_mask)
        if(status==nf90_noerr.and.(lat_mask(2)-lat_mask(1))*(Rlat(2)-Rlat(1))<0)then
           Reverse_lat_direction_Mask=.true.
           if(debug) write(*,*)'ReadCDF mask: reverting latitude direction'
        else
           if(debug) write(*,*)'ReadCDF mask: not reverting latitude direction',&
             (lat_mask(2)-lat_mask(1))*(Rlat(2)-Rlat(1)),status==nf90_noerr
        end if
        lon_shift_Mask=0
        status=nf90_get_var(ncFileID_Mask, dimids_Mask(1), lon_mask)
        if(status==nf90_noerr)then
           lon_shift_Mask=nint((Rlon(1)-lon_mask(1))*dRloni)
           if(lon_shift_Mask/=0)then
              write(*,*)'ReadCDF mask: should shifting longitude by ',lon_shift_Mask,'=',Rlon(1)-lon_mask(1),'degrees'
              call StopAll("Longitude shift for Mask not implemented")
           end if
        end if
        startlat_Mask=startvec(2)    
        if(Reverse_lat_direction_Mask)startlat_Mask=jsize-startvec(2)-dims(2)+2
        call check(nf90_get_var(ncFileID_Mask, VarID_Mask, Mask_values,&
             start=(/startvec(1)+lon_shift_Mask,startlat_Mask/),count=(/dims(1),dims(2)/)),&
          errmsg="Maskvalues")
        !reduce values according to Mask conditions:
       if(present(Mask_Code))then
         factor=0.0
         if(present(Mask_ReducFactor))factor=Mask_ReducFactor
         92 format(A,F5.2,A,500I4)
         if(debug.and.N>0)write(*,92)'multiplying values by ',factor,'where mask has one of values: ',Mask_Code(1:N)
         j=0
           do i=1,totsize
              im=mod(i-1,dims(1))+1
              jm=mod(i-1,dims(1)*dims(2))/dims(1)+1
              ijm=im+(jm-1)*dims(1)
              if(Reverse_lat_direction_Mask)ijm=im+(dims(2)-jm)*dims(1)
              if( any(Mask_Code(1:N)==Mask_values(ijm)))then
                 Rvalues(i)=Rvalues(i)*factor
                 j=j+1
              end if
           end do
          if ( debug )write(*,*)'reduced ',j,' values by ',factor
       else
           factor=1.0
           if(present(Mask_ReducFactor))factor=Mask_ReducFactor
           if ( debug ) write(*,*)'multiplying values by ',factor,'and Mask '
           do i=1,totsize
              Rvalues(i)=Rvalues(i)*Mask_values(i)*factor
           end do
        end if
        if(.not. present(ncFileID_given))call check(nf90_close(ncFileID_Mask))
        deallocate(Mask_values,lon_mask,lat_mask)
        end if
     end if


     !test if this is "fractions" type data
     fractions=.false.
     if(present(fractions_out).or.present(CC_out).or.present(Ncc_out))then
        if ( debug ) write(*,*) 'ReadField_CDF, fraction arrays  '
        if(.not.(present(fractions_out).and.present(CC_out).and.present(Ncc_out)))then
           write(*,*) 'Fraction interpolation missing some arrays of arrays fractions_out CC_out Ncc_out',&
                present(fractions_out),present(CC_out),present(Ncc_out)
        end if
        fractions=.true.
     end if
     if(fractions)then
        if ( debug ) write(*,*) 'fractions method. Reading data '
        Nstartvec=startvec!set 2 first dimensions
        Nstartvec(3)=1
        NCdims=dims!set 2 first dimensions
        !find size of dimension for N (max number of countries per gridcell)
        call check(nf90_inq_dimid(ncid = ncFileID, name = "N", dimID = NdimID))
        call check(nf90_inquire_dimension(ncid=ncFileID,dimID=NdimID,len=Nmax))
        NCdims(3)=Nmax

        allocate(NCC(dims(1)*dims(2)), stat=alloc_err)
        allocate(CC(dims(1)*dims(2),Nmax), stat=alloc_err)     
        allocate(fraction_in(dims(1)*dims(2),Nmax), stat=alloc_err)     

        call check(nf90_inq_varid(ncid = ncFileID, name = 'NCodes', varID = VarIDNCC),&
             errmsg="NCodes not found")

        call check(nf90_get_var(ncFileID, VarIDNCC,NCC ,start=startvec,count=dims),&
             errmsg="Nvalues")

        call check(nf90_inq_varid(ncid = ncFileID, name = 'Codes', varID = VarIDCC),&
             errmsg="Codes not found")
        call check(nf90_get_var(ncFileID, VarIDCC,CC ,start=Nstartvec,count=NCdims),&
             errmsg="CCvalues")

        call check(nf90_inq_varid(ncid = ncFileID, name = 'fractions_'//trim(varname), varID = VarIDfrac),&
             errmsg="fractions not found")
        call check(nf90_get_var(ncFileID, VarIDfrac,fraction_in ,start=Nstartvec,count=NCdims),&
             errmsg="fractions")

!        if( debug )then
!           write(*,*)'More than 2 countries:'
!           do i=1,dims(1)*dims(2)
!              if(NCC(i)>2)write(*,77)me,i,NCC(i),CC(i,1),fraction_in(i,1),CC(i,NCC(i)),fraction_in(i,NCC(i))
!77         format(3I7,2(I5,F6.3))
!           end do
!        end if

        Ncc_out(1:LIMAX*LJMAX)=0
        CC_out(1:LIMAX*LJMAX,1:Nmax)=0
        fractions_out(1:LIMAX*LJMAX,1)=0.0
        fractions_out(1:LIMAX*LJMAX,2:Nmax)=0.0
     end if


     if ( DEBUG_NETCDF_RF ) write(*,*) 'ReadCDF types ', &
          xtype, NF90_INT, NF90_SHORT, NF90_BYTE

     if(xtype==NF90_INT.or.xtype==NF90_SHORT.or.xtype==NF90_BYTE)then
        !scale data if it is packed
        scalefactors(1) = 1.0 !default
        scalefactors(2) = 0.  !default
        status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
        if(status == nf90_noerr) scalefactors(1) = scale
        status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
        if(status == nf90_noerr) scalefactors(2) = offset
        Rvalues=Rvalues*scalefactors(1)+scalefactors(2)
        FillValue=FillValue*scalefactors(1)+scalefactors(2)
        if ( debug ) then
           write(*,*)' Start scaling mpixtype',xtype
           write(*,*)' FillValue scaled to',FillValue
           write(*,*)' Max(RValues)   ',maxval(RValues)
        end if
     else ! Real
        if ( debug )write(*,*)' xtype real ',xtype
     end if

     if(interpol_used=='conservative'.or.interpol_used=='mass_conservative')then
        if ( debug )write(*,*)'interpol case ',interpol_used
        
        if(projection=='lon lat')then
         if ( debug )write(*,*)'interpol case ',trim(interpol_used),' projection lon lat'
          !exact integrals assuming uniform emission over emitter gridcells
           if(.not.allocated(fracfirstlon))then
              allocate(fracfirstlon(dims(1)),fraclastlon(dims(1)),ifirst(dims(1)),ilast(dims(1)))
              allocate(fracfirstlat(dims(2)),fraclastlat(dims(2)),jfirst(dims(2)),jlast(dims(2)))
           endif
           !precompute factors 1-dimensionally
           dRlon = Rlon(2)-Rlon(1)
           if(dRlon<0)dRlon = dRlon+360.0
           dRloni=1.0/dRlon
           dlon=glon(2,1)-glon(1,1)
           if(dlon<0)dlon = dlon+360.0
           dloni=1.0/dlon
           dRlat = Rlat(2)-Rlat(1)
           dRlati=1.0/dRlat
           dlat=glat(1,2)-glat(1,1)
           dlati=1.0/dlat

           do ig=1,dims(1)
              !longitude of edges
              Rlonmin = Rlon(ig+imin-1) - 0.5*dRlon
              Rlonmax = Rlonmin + dRlon
              !find all cells with at least some part in emitter cell ig
              call lb2ij(Rlonmin,0.0,ir,jr)
              ifirst(ig)=floor(ir+0.5+1.E-6)!first i to be treated
              call lb2ij(Rlonmax,0.0,ir,jr)
              ilast(ig)=floor(ir+0.5-1.E-6)!last i to be treated
              
              !end of the world tests
              if(ilast(ig)>= ifirst(ig))then
                 !no problems with monotonicity
              else                 
                 !Problems: we cross the point where lon+eps=lon-360
                 !several cases according to which points are in the subdomain
                  ilast_corrected=.false.
                 if(ifirst(ig)>=rundomain(1) .and. ifirst(ig)<=rundomain(2))then
                    !inside rundomain
                    if(i_local(ifirst(ig))>=1 .and. i_local(ifirst(ig))<=limax)then
                       !inside subdomain. ilast(ig) is wrong, set it at the end of subdomain
                       ilast(ig)=i_fdom(limax)
                       !we are done with correction
                       ilast_corrected=.true.
                    endif
                 endif
                 if((.not. ilast_corrected) .and. ilast(ig)>=rundomain(1) .and. ilast(ig)<=rundomain(2))then
                    !inside rundomain
                    if(i_local(ilast(ig))>=1 .and. i_local(ilast(ig))<=limax)then
                       !inside subdomain. ifirst(ig) is wrong, set it at the start of subdomain
                       ifirst(ig)=i_fdom(1)
                    endif
                 endif
                 if(ilast(ig)< ifirst(ig))then
                    !none of them are in subdomain. nothing to do
                    ilast(ig)=ifirst(ig)-1
                    cycle
                 endif
              endif              

              if((ifirst(ig)>rundomain(2) .or. ilast(ig)<rundomain(1)))then
                 !outside rundomain no need to spend time with this ig
                 ilast(ig)=ifirst(ig)-1
                 cycle
              endif

               ! convert to local coordinates:
              ifirst(ig)=max(ifirst(ig),rundomain(1)) ! first trim to rundomain
              ilast(ig) =min(ilast(ig) ,rundomain(2))
              ifirst(ig)=max(i_local(ifirst(ig)),1)   ! convert to local coordinates
              ilast(ig) =min(i_local(ilast(ig)) ,limax)
              if((ifirst(ig)>limax .or. ilast(ig)<1) )then
                 ! outside local domain. no need to spend time with this ig
                 ilast(ig)=ifirst(ig)-1
                 cycle
              endif

              !make fraction of overlap. Only first or last i can overlap. 1*dloni means 100% of the incoming data is taken.
              !put all incoming longitudes in the range 0-360
!              fracfirstlon(ig) = min(1.0, mod(360.0 + mod(glon(ifirst(ig),1)+360.0,360.0)+ 0.5*dlon - mod(Rlonmin+360.0,360.0),360.0)*dloni)
              fracfirstlon(ig) = mod(360.0 + mod(glon(ifirst(ig),1)+360.0,360.0)&
                                + 0.5*dlon - mod(Rlonmin+360.0,360.0),360.0)*dloni

              if(fracfirstlon(ig)<0.0 .or. fracfirstlon(ig)>10.0*Resolution_fac)then
                 fracfirstlon(ig)=0.0!numerical noise in glon or i=1
              endif
              fracfirstlon(ig)=min(1.0,fracfirstlon(ig))
 
!             fraclastlon(ig)  = min(1.0, mod(360.0 + mod(Rlonmax+360.0,360.0) - (mod(glon(ilast(ig),1)+360.0,360.0) - 0.5*dlon),360.0)*dloni)
              fraclastlon(ig)  =  mod(360.0 + mod(Rlonmax+360.0,360.0) &
                  - (mod(glon(ilast(ig),1)+360.0,360.0) - 0.5*dlon),360.0)*dloni
              if(fraclastlon(ig)<0.0 .or. fraclastlon(ig)>10.0*Resolution_fac)then
                 fraclastlon(ig)=0.0!numerical noise or i=limax
              endif
              fraclastlon(ig)=min(1.0,fraclastlon(ig))

              !NB: when reducing on both sides need to ADD reductions not multiply
              if(ifirst(ig)==ilast(ig))fraclastlon(ig)  = fraclastlon(ig) -(1.0-fracfirstlon(ig))!include reduction from both sides
              if(fraclastlon(ig)<0.0 .or. fraclastlon(ig)>10.0*Resolution_fac)then
                 fraclastlon(ig)=0.0!numerical noise or i=limax
              endif

              if(fracfirstlon(ig)<0.0)write(*,*)'ERROR A in interpolation',me,ig,fracfirstlon(ig)
              if(fraclastlon(ig)<0.0)write(*,*)'ERROR B in interpolation',me,ig,fraclastlon(ig)
!631           format(I4,A,F10.4,A,F10.4,A,I4,A,F10.4,A,I4,A,F10.4,A,F10.4)
!              if(me==0 .and. ig==1.and. trim(varname)=='sox_sec01')write(*,631)ig,' start '//trim(varname),Rlonmin,' end',Rlonmax,' firsti',ifirst(ig),'lon ',glon(ifirst(ig),1),'last i',ilast(ig),'frac ',fracfirstlon(ig),' and ',fraclastlon(ig)
           enddo

           !make factors for j
           do jg=1,dims(2)
              !latitude of edges. NB:Rlat is over fullgrid, while jg is in restricted grid

              Rlatmin = Rlat(jg+jmin-1) - 0.5*abs(dRlat)
              Rlatmax = Rlatmin + abs(dRlat)
              !find all cells with at least some part in emitter cell jg
              call lb2ij(0.0,Rlatmin,ir,jr)
              jfirst(jg)=floor(jr+0.5+1.E-6)!first j to be treated
              call lb2ij(0.0,Rlatmax,ir,jr)
              jlast(jg)=floor(jr+0.5-1.E-6)!last j to be treated
              !write(*,*)trim(varname)//' ',me,jg,jfirst(jg),Rlatmin,Rlatmax,jr,dims(2)
              if(jfirst(jg)>rundomain(4) .or. jlast(jg)<rundomain(3))then
                 !outside rundomain, no need to spend time with this jg
                 jlast(jg)=jfirst(jg)-1
                 cycle
              endif

              jfirst(jg)=max(1,j_local(max(1,jfirst(jg))))
              jlast(jg)=min(ljmax,j_local(min(rundomain(4),jlast(jg))))
              if(jfirst(jg)>ljmax .or. jlast(jg)<1)then
                 !no need to spend time with this jg
                 jlast(jg)=jfirst(jg)-1
                 cycle
              endif

              fracfirstlat(jg) = min(1.0,(glat(1,jfirst(jg))+ 0.5*dlat - Rlatmin)*dlati)
              fraclastlat(jg)  = min(1.0,(Rlatmax - (glat(1,jlast(jg)) - 0.5*dlat))*dlati)
              fracfirstlat(jg) = max(0.0,fracfirstlat(jg))
              fraclastlat(jg)  = max(0.0,fraclastlat(jg))
              if(jfirst(jg)==jlast(jg))fraclastlat(jg)  = fraclastlat(jg) -(1.0-fracfirstlat(jg))!include reduction from both sides

           enddo
           K2=1
           if(data3D)k2=kend_loc-kstart_loc+1
           ijk=LIMAX*LJMAX*k2
           Rvar(1:ijk)=0.0
           allocate(Ivalues(ijk))
           Ivalues = 0 !just to see if Rvar gets at least one value

           idiv=0
           do jg=1,dims(2)
             do j=jfirst(jg),jlast(jg)
              if(j>=1.and.j<=ljmax)then
                 if(interpol_used=='mass_conservative')then
                    !scale for gridcell size differences

                    frac = abs(dRlati*dRloni*dlat*dlon)!Divide by number of model cell in readin cell
                 else
                    !will give average value (emissions in kg/m2 for instance)
                    frac = 1.0
                 endif
                 frac_j = frac

                 if(j==jfirst(jg))frac_j=frac*fracfirstlat(jg)
                 if(j==jlast(jg))frac_j=frac*fraclastlat(jg)!fracfirstlat not used!

                 do ig=1,dims(1)
                    igjg=ig+(jg-1)*dims(1)
                    
                    do i=ifirst(ig),ilast(ig)
                       if(i>=1.and.i<=limax)then

                          frac = frac_j
                          if(i==ifirst(ig))frac=frac_j*fracfirstlon(ig)
                          if(i==ilast(ig))frac=frac_j*fraclastlon(ig)!fracfirstlon not used!

                          ij=i+(j-1)*LIMAX
                          k2=1
                          if(data3D)k2=kend_loc-kstart_loc+1
                          do k=1,k2
                             ijk=k+(ij-1)*k2
                             igjgk=igjg+(k-1)*dims(1)*dims(2)
                             Ivalues(ijk)=Ivalues(ijk)+1

                             if(fractions)then
                                call readfrac(Ncc(igjgk),CC,Rvalues(igjgk),fraction_in,fractions_out,Ncc_out,CC_out,&
                                     Rvar(ijk),dims(1)*dims(2),igjgk,ijk,frac,Reduc)
                                
                             elseif(OnlyDefinedValues.or.Rvalues(igjgk)/=FillValue)then
                                Rvar(ijk)=Rvar(ijk)+Rvalues(igjgk)*frac
                             else
                                !Not defined: don't include this Rvalue
                                Rvar(ijk) = UnDef_local
                             end if
                          end do
                       end if
                    enddo
                 enddo
              endif              
           enddo

        enddo
        do ijk=1,LIMAX*LJMAX*k2
           if(Ivalues(ijk)<=0) Rvar(ijk) = UnDef_local
        enddo
        deallocate(Ivalues)
        deallocate(fracfirstlon,fraclastlon,ifirst,ilast)
        deallocate(fracfirstlat,fraclastlat,jfirst,jlast)
        
     else

        !conserves integral (almost, does not take into account local differences in mapping factor)
        !takes weighted average over gridcells covered by model gridcell

        !divide the coarse grid into pieces significantly smaller than the fine grid
        !Divide each global gridcell into Ndiv x Ndiv pieces
        Ndiv=nint(5*Grid_resolution/GRIDWIDTH_M)
!        if(interpol_used=='conservative')Ndiv=nint(2*Grid_resolution/GRIDWIDTH_M)!Will be smooth anyway
        Ndiv=max(1,Ndiv)
        Ndiv2=Ndiv*Ndiv
        Ndiv_lon=nint(5*Grid_resolution_lon/GRIDWIDTH_M)
        Ndiv_lon=max(1,Ndiv_lon)
        Ndiv2=Ndiv_lon*Ndiv
        !
        if(projection/='Stereographic'.and.projection/='lon lat'.and.projection/='Rotated_Spherical'.and.projection/='lambert')then
           !the method should be revised or used only occasionally
           if(me==0)write(*,*)'WARNING: interpolation method may be CPU demanding'
        end if
        k2=1
        if(data3D)k2=kend_loc-kstart_loc+1
        allocate(Ivalues(LIMAX*LJMAX*k2))
        allocate(Nvalues(LIMAX*LJMAX*k2))
        do ij=1,LIMAX*LJMAX*k2
           Ivalues(ij)=0
           NValues(ij) = 0
           Rvar(ij)=0.0
        end do

        do jg=1,dims(2)
           do jdiv=1,Ndiv
              lat=Rlat(startvec(2)-1+jg)-0.5/dRlati+(jdiv-0.5)/(dRlati*Ndiv)
              do ig=1,dims(1)
                 igjg=ig+(jg-1)*dims(1)
                 do idiv=1,Ndiv_lon
                    lon=Rlon(startvec(1)-1+ig)-0.5/dRloni+(idiv-0.5)/(dRloni*Ndiv_lon)
                    call lb2ij(lon,lat,i,j)
                    i=i-gi0-IRUNBEG+2
                    j=j-gj0-JRUNBEG+2

                    if(i>=1.and.i<=limax.and.j>=1.and.j<=ljmax)then
                       ij=i+(j-1)*LIMAX
                       k2=1
                       if(data3D)k2=kend_loc-kstart_loc+1
                       do k=1,k2
                          ijk=k+(ij-1)*k2
                          Ivalues(ijk)=Ivalues(ijk)+1
                          Nvalues(ijk)=Nvalues(ijk)+1
                          igjgk=igjg+(k-1)*dims(1)*dims(2)
                          latlon_weight=1.0

                          if(fractions)then
                             call readfrac(Ncc(igjgk),CC,Rvalues(igjgk),&
                                  fraction_in,fractions_out,Ncc_out,CC_out,&
                                  Rvar(ijk),dims(1)*dims(2),igjgk,ijk,&
                                  latlon_weight,Reduc)
                          elseif(OnlyDefinedValues.or.Rvalues(igjgk)/=FillValue)then
                             Rvar(ijk)=Rvar(ijk)+Rvalues(igjgk)
                          else
                             !Not defined: don't include this Rvalue
                             Ivalues(ijk)=Ivalues(ijk)-1                             
                          end if
                       end do

                    end if
                 end do
              end do
           end do
        end do

        k2=1
        if(data3D)k2=kend_loc-kstart_loc+1
        do k=1,k2
           do i=1,limax
              do j=1,ljmax
                 ij=i+(j-1)*LIMAX
                 ijk=k+(ij-1)*k2

                 debug_ij = ( DEBUG_NETCDF_RF .and. debug_proc .and. &
                      i== debug_li .and. j== debug_lj )
                 if ( debug_ij ) write(*,"(a,9i5)") 'DEBUG  -- INValues!',&
                      Ivalues(ijk), Nvalues(ijk), me, i,j,k
                 if(Ivalues(ijk)<=0)then
                    if( .not.present(UnDef))then
                       write(*,"(a,a,4i4,6f10.3,2i6,6f10.3)") &
                            'ERROR, NetCDF_mod no values found here ! ', &
                            trim(fileName) // ":" // trim(varname), &
                            i,j,k,me,minlon,maxlon,minlat,maxlat,glon(i,j),glat(i,j), &
                            Ivalues(ijk),Ndiv,Rlon(startvec(1)),Rlon(startvec(1)+dims(1)-1),&
                            Rlat(startvec(2)),Rlat(startvec(2)-1+dims(2))
                       call CheckStop("Interpolation error")
                    else                       
                       Rvar(ijk)=UnDef_local
                    end if
                 else
                    if(interpol_used=='mass_conservative')then
                       !used for example for emissions in kg (or kg/s)
                       Rvar(ijk)=Rvar(ijk)/Ndiv2! Total sum of values from all cells is constant
                       if ( debug_ij ) write(*,"(a,a,3i5,es12.4)") 'DEBUG  -- mass!' , &
                            trim(varname), Ivalues(ijk), Nvalues(ijk), Ndiv2, Rvar(ijk)
                    else
                       !used for example for emissions in kg/m2 (or kg/m2/s)
                       ! integral is approximately conserved
                       Rvar(ijk)=Rvar(ijk)/Ivalues(ijk)
                       if ( debug_ij ) write(*,"(a,a,3i5,es12.4)") &
                            'DEBUG  -- approx!' ,  trim(varname),&
                            Ivalues(ijk), Nvalues(ijk),Ndiv2, Rvar(ijk)

                    end if
                 end if
              end do
           end do
        end do
        deallocate(Ivalues)
        deallocate(Nvalues)
     endif


  elseif(interpol_used=='zero_order')then
        !interpolation 1:
        !nearest gridcell
        ijk=0
        k2=1
        if(data3D)k2=kend_loc-kstart_loc+1
        do k=1,k2
           do j=1,ljmax
              do i=1,limax
                 ij=i+(j-1)*LIMAX
                 ijk=k+(ij-1)*k2
                 ig=nint((glon(i,j)-Rlon(startvec(1)))*dRloni)+1
                 if(ig<0.5)ig=ig+dims(1)
                 if(ig>dims(1))ig=ig-dims(1)
                 ig=max(1,min(dims(1),ig))
                 jg=max(1,min(dims(2),nint((glat(i,j)-Rlat(startvec(2)))*dRlati)+1))
                 igjgk=ig+(jg-1)*dims(1)+(k-1)*dims(1)*dims(2)
                 if(OnlyDefinedValues.or.(Rvalues(igjgk)/=FillValue.and. .not.isnan(Rvalues(igjgk))))then
                    Rvar(ijk)=Rvalues(igjgk)
                 else
                    Rvar(ijk)=UnDef_local
                 end if
              end do
           end do
        end do
     else
        write(*,*)'interpolation method not implemented'
     end if
    !_________________________________________________________________________________________________________
     !_________________________________________________________________________________________________________
  elseif(data_projection=="Stereographic")then
     !we assume that data is originally in Polar Stereographic projection
     if(debug1)write(*,*)'interpolating from ', trim(data_projection),' to ',trim(projection)

     !get coordinates
     !check that there are dimensions called i and j
     call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(1), name=name ),name)
     call CheckStop(trim(name)/='i',"i not found")
     call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(2), name=name ),name)
     call CheckStop(trim(name)/='j',"j not found")

     call CheckStop(data3D,"3D data in Stereographic projection not yet implemented")

     if(present(Grid_resolution_in))then
        Grid_resolution = Grid_resolution_in
     else
        status=nf90_get_att(ncFileID, nf90_global, "Grid_resolution", Grid_resolution )
        if(status /= nf90_noerr)then
           Grid_resolution=GRIDWIDTH_M_EMEP
           if ( debug )write(*,*)'Grid_resolution assumed =',Grid_resolution
        end if
     endif
     !the method chosen depends on the relative resolutions
     if(interpol_used=='conservative'.and.Grid_resolution/GRIDWIDTH_M>2)then
        interpol_used='zero_order'!usually good enough, and keeps gradients
        if ( MasterProc .and. debug) write(*,*) 'Asked for conservative interpolation, but redefined as ',interpol_used
     end if

     status = nf90_get_att(ncFileID, nf90_global, "xcoordinate_NorthPole", xp_ext )
     if(status /= nf90_noerr)then
        xp_ext=xp_EMEP_old
        if ( debug )write(*,*)'xcoordinate_NorthPole assumed =',xp_ext
     end if
     status=nf90_get_att(ncFileID, nf90_global, "ycoordinate_NorthPole", yp_ext )
     if(status /= nf90_noerr)then
        yp_ext=yp_EMEP_old
        if ( debug )write(*,*)'ycoordinate_NorthPole assumed =',yp_ext
     end if
     status=nf90_get_att(ncFileID, nf90_global, "fi", fi_ext )
     if(status /= nf90_noerr)then
        fi_ext=fi_EMEP
        if ( debug )write(*,*)'fi assumed =',fi_ext
     end if
     status=nf90_get_att(ncFileID, nf90_global, "ref_latitude", ref_lat_ext  )
     if(status /= nf90_noerr)then
        ref_lat_ext=ref_latitude_EMEP
        if ( debug )write(*,*)'ref_latitude assumed =',ref_lat_ext
     end if
     an_ext=EARTH_RADIUS*(1.0+sin(ref_lat_ext*PI/180.0))/Grid_resolution

     !read entire grid in a first implementation
     startvec=1
     totsize=1

     do i=1,ndims
        totsize=totsize*dims(i)
     end do
     if ( debug )write(*,*)'totsize ',totsize,ndims
     allocate(Rvalues(totsize), stat=alloc_err)
     call check(nf90_get_var(ncFileID, VarID, Rvalues,start=startvec,count=dims),&
          errmsg="RRvaluesStereo")
     if(xtype==NF90_INT.or.xtype==NF90_SHORT.or.xtype==NF90_BYTE)then
        !scale data if it is packed
        scalefactors(1) = 1.0 !default
        scalefactors(2) = 0.  !default
        status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
        if(status == nf90_noerr) scalefactors(1) = scale
        status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
        if(status == nf90_noerr) scalefactors(2) = offset
        Rvalues=Rvalues*scalefactors(1)+scalefactors(2)
        FillValue=FillValue*scalefactors(1)+scalefactors(2)
        if ( debug ) then
           write(*,*)' Start scaling mpixtype',xtype
           write(*,*)' FillValue scaled to',FillValue
           write(*,*)' Max(RValues)   ',maxval(RValues)
        end if
     else ! Real
        if ( debug ) then
           write(*,*)' xtype real ',xtype
           write(*,*)' FillValue still',FillValue
           write(*,*)' Max(RValues)   ',maxval(RValues)
           write(*,*)' Min(RValues)   ',minval(RValues)
        end if
     end if

     if(interpol_used=='conservative'.or.interpol_used=='mass_conservative')then
        !conserves integral (almost, does not take into account local differences in mapping factor)
        !takes weighted average over gridcells covered by model gridcell

        !divide the external grid into pieces significantly smaller than the fine grid
        !Divide each global gridcell into Ndiv x Ndiv pieces
        Ndiv=nint(3*Grid_resolution/GRIDWIDTH_M)
        Ndiv=max(1,Ndiv)
        Ndiv2=Ndiv*Ndiv
        Grid_resolution_div=Grid_resolution/Ndiv
        xp_ext_div=(xp_ext-0.5)*Ndiv+0.5
        yp_ext_div=(yp_ext-0.5)*Ndiv+0.5
        an_ext_div=an_ext*Ndiv

        if(projection/='Stereographic'.and.projection/='lon lat'.and.projection/='Rotated_Spherical'.and.projection/='lambert')then
           !the method should be revised or used only occasionally
           if(me==0)write(*,*)'WARNING: interpolation from stereographic to '//trim(projection)//'may be CPU demanding:'
        end if
        k2=1
        if(data3D)k2=kend_loc-kstart_loc+1
        allocate(Ivalues(LIMAX*LJMAX*k2))
        allocate(Nvalues(LIMAX*LJMAX*k2))
        do ij=1,LIMAX*LJMAX*k2
           Ivalues(ij)=0
           NValues(ij) = 0
           !           if(present(UnDef))then
           !              Rvar(ij)=UnDef!default value
           !           else
           Rvar(ij)=0.0
           !           end if
        end do

        do jg=1,dims(2)
           do jdiv=1,Ndiv
              j_ext=(jg-1)*Ndiv+jdiv
              do ig=1,dims(1)
                 igjg=ig+(jg-1)*dims(1)
                 do idiv=1,Ndiv
                    i_ext=(ig-1)*Ndiv+idiv
                    call ij2lb(i_ext,j_ext,lon,lat,fi_ext,an_ext_div,xp_ext_div,yp_ext_div)
                    call lb2ij(lon,lat,i,j)!back to model (fulldomain) coordinates

!                    if(abs(lat-57.0)<0.01 .and. abs(lon-1.3)<0.01)write(*,*)'fullij ',lat,lon,me,i,j
!                    if(abs(lon-15)<0.02 .and. abs(lat-63)<0.02)write(*,*)jg,ig,lon,lat,i,j,me
                    !convert from fulldomain to local domain
                    !i,j may be any integer, therefore should not use i_local array
                    i=i-gi0-IRUNBEG+2
                    j=j-gj0-JRUNBEG+2


!83                  format(2I4,31F9.2)
                    !if ( debug .and.me==0) write(*,83)i,j,lon,lat,fi_ext,an_ext_div,xp_ext_div,yp_ext_div,fi,xp,yp,Rvalues(igjg)

                    if(i>=1.and.i<=limax.and.j>=1.and.j<=ljmax)then
                       ij=i+(j-1)*LIMAX
                       k2=1
                       if(data3D)k2=kend_loc-kstart_loc+1
                       do k=1,k2
                          ijk=k+(ij-1)*k2
                          Ivalues(ijk)=Ivalues(ijk)+1
                          Nvalues(ijk)=Nvalues(ijk)+1
                          igjgk=igjg+(k-1)*dims(1)*dims(2)

                          if(OnlyDefinedValues.or.Rvalues(igjgk)/=FillValue)then
                             Rvar(ijk)=Rvar(ijk)+Rvalues(igjgk)
                          else
                             !Not defined: don't include this Rvalue
                             Ivalues(ijk)=Ivalues(ijk)-1

                          end if
                       end do
                    end if
                 end do
              end do
           end do
        end do
        k2=1
        if(data3D)k2=kend_loc-kstart_loc+1
        do k=1,k2
           do i=1,limax
              do j=1,ljmax
                 ij=i+(j-1)*LIMAX
                 ijk=k+(ij-1)*k2

                 debug_ij = ( DEBUG_NETCDF_RF .and. debug_proc .and. &
                      i== debug_li .and. j== debug_lj )
                 if ( debug_ij ) write(*,*) 'DEBUG  -- INValues!', &
                      Ivalues(ijk), Nvalues(ijk)
                 if(Ivalues(ijk)<=0.)then
                    if( .not.present(UnDef))then
                       write(*,"(a,a,4i4,6g10.3,i6)") &
                            'ERROR, NetCDF_mod no values found! B', &
                            trim(fileName) // ":" // trim(varname), &
                            i,j,k,me,maxlon,minlon,maxlat,minlat,glon(i,j),glat(i,j), &
                            Ivalues(ijk)
                       call CheckStop("Interpolation error")
                    else
                       Rvar(ijk)=UnDef
                    end if
                 else
                    if(interpol_used=='mass_conservative')then
                       !used for example for emissions in kg (or kg/s)
                       Rvar(ijk)=Rvar(ijk)/Ndiv2! Total sum of values from all cells is constant
                       if ( debug_ij ) write(*,"(a,a,3i5,es12.4)") 'DEBUG  -- mass!' , &
                            trim(varname), Ivalues(ijk), Nvalues(ijk), Ndiv2, Rvar(ijk)
                    else
                       !used for example for emissions in kg/m2 (or kg/m2/s)
                       ! integral is approximately conserved
                       Rvar(ijk)=Rvar(ijk)/Ivalues(ijk)
                       if ( debug_ij ) write(*,"(a,a,3i5,es12.4)") &
                            'DEBUG  -- approx!' ,  trim(varname),&
                            Ivalues(ijk), Nvalues(ijk),Ndiv2, Rvar(ijk)

                    end if
                 end if
              end do
           end do
        end do

        deallocate(Ivalues)
        deallocate(Nvalues)

     elseif(interpol_used=='zero_order')then
        !interpolation 1:
        !nearest gridcell

        Ndiv=1
        Grid_resolution_div=Grid_resolution/Ndiv
        xp_ext_div=(xp_ext-0.5)*Ndiv+0.5
        yp_ext_div=(yp_ext-0.5)*Ndiv+0.5
        an_ext_div=an_ext*Ndiv
        if(MasterProc.and.debug)write(*,*)'zero_order interpolation ',an_ext_div,xp_ext_div,yp_ext_div,dims(1),dims(2)

        if(projection/='Stereographic'.and.projection/='lon lat'.and.projection/='Rotated_Spherical'.and.projection/='lambert')then
           !the method should be revised or used only occasionally
           if(me==0)write(*,*)'WARNING: interpolation method may be CPU demanding',projection
        end if


        call lb2ijm(LIMAX,LJMAX,glon,glat,buffer1,buffer2,fi_ext,an_ext_div,xp_ext_div,yp_ext_div)
        i_ext=nint(buffer1(1,1))
        j_ext=nint(buffer2(1,1))
        call ij2lb(i_ext,j_ext,lon,lat,fi_ext,an_ext_div,xp_ext_div,yp_ext_div)
        k2=1
        if(data3D)k2=kend_loc-kstart_loc+1
        do j=1,ljmax
           do i=1,limax
              ij=i+(j-1)*LIMAX
              i_ext=nint(buffer1(i,j))
              j_ext=nint(buffer2(i,j))
              if(i_ext>=1.and.i_ext<=dims(1).and.j_ext>=1.and.j_ext<=dims(2))then

                 do k=1,k2
                    ijk=k+(ij-1)*k2

                    igjgk=i_ext+(j_ext-1)*dims(1)+(k-1)*dims(1)*dims(2)

                    if(OnlyDefinedValues.or.(Rvalues(igjgk)/=FillValue.and. .not.isnan(Rvalues(igjgk))))then
                       Rvar(ijk)=Rvalues(igjgk)
                    else
                       if(present(UnDef))then
                          Rvar(ijk)=UnDef!default value
                       else
                          Rvar(ijk)=Rvalues(igjgk)
                       end if
                    end if
                 end do
              else
                 do k=1,k2
                    ijk=k+(ij-1)*k2
                    if(present(UnDef))then
                       Rvar(ijk)=UnDef!default value
                    else
                       !                    if ( debug ) write(*,*)'WARNING: gridcell out of map. Set to ',FillValue
                       call StopAll("ReadField_CDF: values outside grid required")
                    end if
                 end do
              end if
           end do
        end do


     end if

  else ! data_projection /="lon lat" .and. data_projection/="Stereographic"

     if(debug1) write(*,*)'interpolatingL from ', &
        trim(data_projection), ' to ',trim(projection)

     if(interpol_used=='conservative'.or.interpol_used=='mass_conservative')then

!        call CheckStop((data3D),"3D data in general projection not yet implemented for conservative interpolation")        
        if(present(Grid_resolution_in))then
           Grid_resolution = Grid_resolution_in
        else
           status=nf90_get_att(ncFileID, nf90_global, "Grid_resolution", Grid_resolution )
           if(status /= nf90_noerr)then
              if(MasterProc)write(*,*)'warning did not find grid resolution in '//trim(fileName)
              call make_gridresolution(ncFileID,Grid_resolution)
           endif
        endif
        Ndiv=nint(5*Grid_resolution/GRIDWIDTH_M)
        Ndiv=max(1,Ndiv)
        Ndiv2=Ndiv*Ndiv
        if(Ndiv>1.and.MasterProc .and.step_main == 1 )then
           write(*,*)'dividing each gridcell into ',Ndiv2,' pieces'
        end if

        k2=1
        if(data3D)then
           if(kstart_loc<0)then
              !take all levels
              kstart_loc=1
              kend_loc=dims(3)
           end if
           startvec(3)=kstart_loc
           k2=kend_loc-kstart_loc+1
           dims(3)=kend_loc-kstart_loc+1
           if(ndims>3)then
              startvec(ndims)=nstart
              dims(ndims)=1
           end if
        else
           if(ndims>2)then
              startvec(ndims)=nstart
              dims(ndims)=1
           end if
        end if

        allocate(Ivalues(LIMAX*LJMAX*k2))
        do ij=1,LIMAX*LJMAX*k2
           Ivalues(ij)=0
           Rvar(ij)=0.0
        end do
        totsize=1
        do i=1,ndims
           totsize=totsize*dims(i)
        end do

        allocate(Rvalues(totsize), stat=alloc_err)
        call check(nf90_get_var(ncFileID, VarID, Rvalues,start=startvec,count=dims)&
             ,errmsg="Read Rvalues failed")


        if(xtype==NF90_INT.or.xtype==NF90_SHORT.or.xtype==NF90_BYTE)then
           !scale data if it is packed
           scalefactors(1) = 1.0 !default
           scalefactors(2) = 0.  !default
           status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
           if(status == nf90_noerr) scalefactors(1) = scale
           status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
           if(status == nf90_noerr) scalefactors(2) = offset
           Rvalues=Rvalues*scalefactors(1)+scalefactors(2)
           FillValue=FillValue*scalefactors(1)+scalefactors(2)
           if ( debug ) then
              write(*,*)' Start scaling mpixtype',xtype
              write(*,*)' FillValue scaled to',FillValue
              write(*,*)' Max(RValues)   ',maxval(RValues)
           end if
        else ! Real
           if ( debug )write(*,*)' xtype real ',xtype
        end if

!not 100% robust: assumes i increases with longitude, j with latitude. and i almost parallel with longitudes
        do jg=1,dims(2)
           do jdiv=1,Ndiv              
              do ig=1,dims(1)
                 igjg=ig+(jg-1)*dims(1)
                 if(Ndiv>1)then
                    if(ig>1.and.jg>1)then
                       !                          dlat=(Rlat(igjg)-Rlat(ig-1+(jg-2)*dims(1)))/Ndiv
                       !                          dlon=(Rlon(igjg)-Rlon(ig-1+(jg-2)*dims(1)))/Ndiv
                       dlaty=(Rlat(igjg)-Rlat(igjg-dims(1)))/Ndiv
                       dlonx=(Rlon(igjg)-Rlon(igjg-1))/Ndiv
                       dlatx=(Rlat(igjg)-Rlat(igjg-1))/Ndiv
                       dlony=(Rlon(igjg)-Rlon(igjg-dims(1)))/Ndiv
                    else
                       if(ig>1)then
                          !                             dlat=-((Rlat(igjg)-Rlat(ig-1+(jg)*dims(1)))/Ndiv)
                          !                             dlon=((Rlon(igjg)-Rlon(ig-1+(jg)*dims(1)))/Ndiv)
                          dlaty=-(Rlat(igjg)-Rlat(igjg+dims(1)))/Ndiv
                          dlonx=(Rlon(igjg)-Rlon(igjg-1))/Ndiv
                          dlatx=(Rlat(igjg)-Rlat(igjg-1))/Ndiv
                          dlony=-(Rlon(igjg)-Rlon(igjg+dims(1)))/Ndiv
                       else if(jg>1)then
                          !                             dlat=((Rlat(igjg)-Rlat(ig+1+(jg-2)*dims(1)))/Ndiv)
!                             dlon=-((Rlon(igjg)-Rlon(ig+1+(jg-2)*dims(1)))/Ndiv)
                          dlaty=(Rlat(igjg)-Rlat(igjg-dims(1)))/Ndiv
                          dlonx=-(Rlon(igjg)-Rlon(igjg+1))/Ndiv
                          dlatx=-(Rlat(igjg)-Rlat(igjg+1))/Ndiv
                          dlony=(Rlon(igjg)-Rlon(igjg-dims(1)))/Ndiv
                       else
                          !                             dlat=-((Rlat(igjg)-Rlat(ig+1+(jg)*dims(1)))/Ndiv)
                          !                             dlon=-((Rlon(igjg)-Rlon(ig+1+(jg)*dims(1)))/Ndiv)
                          dlaty=-(Rlat(igjg)-Rlat(igjg+dims(1)))/Ndiv
                          dlonx=-(Rlon(igjg)-Rlon(igjg+1))/Ndiv
                          dlatx=-(Rlat(igjg)-Rlat(igjg+1))/Ndiv
                          dlony=-(Rlon(igjg)-Rlon(igjg+dims(1)))/Ndiv
                       end if
                    end if
                 end if
                 do idiv=1,Ndiv
                    if(Ndiv>1)then
                       lon=Rlon(igjg)+dlonx*(idiv-0.5-0.5*Ndiv)+dlony*(jdiv-0.5-0.5*Ndiv)
                       lat=Rlat(igjg)+dlaty*(jdiv-0.5-0.5*Ndiv)+dlatx*(idiv-0.5-0.5*Ndiv)
                    else
                       lon=Rlon(igjg)
                       lat=Rlat(igjg)
                    end if
                    call lb2ij(lon,lat,i,j)!back to model (fulldomain) coordinates
                    i=i-gi0-IRUNBEG+2
                    j=j-gj0-JRUNBEG+2
                    if(i>=1.and.i<=limax.and.j>=1.and.j<=ljmax)then
                       ij=i+(j-1)*LIMAX
                       do k=1,k2
                          ijk=k+(ij-1)*k2
                          Ivalues(ijk)=Ivalues(ijk)+1
                          igjgk=igjg+(k-1)*dims(1)*dims(2)

                          if(OnlyDefinedValues.or.Rvalues(igjg)/=FillValue)then
                             Rvar(ijk)=Rvar(ijk)+Rvalues(igjgk)
                          else
                             !Not defined: don't include this Rvalue
                             Ivalues(ijk)=Ivalues(ijk)-1
                          end if
                       end do
                    end if
                 end do
              end do
           end do
        end do

        do k=1,k2
        do i=1,limax
           do j=1,ljmax
              ij=i+(j-1)*LIMAX
              ijk=k+(ij-1)*k2

              if(Ivalues(ijk)<=0.)then
                 if( .not.present(UnDef))then
                    write(*,"(a,a,4i4,6g10.3,2i6)") &
                         'ERROR, NetCDF_mod no values found! C', &
                         trim(fileName) // ":" // trim(varname), &
                         i,j,k,me,minlon,maxlon,minlat,maxlat,glon(i,j),glat(i,j), &
                         Ivalues(ijk),Ndiv
                    call CheckStop("Interpolation error")
                 else
                    Rvar(ijk)=UnDef
                 end if
              else
                 if(interpol_used=='mass_conservative')then
                    !used for example for emissions in kg (or kg/s)
                    Rvar(ijk)=Rvar(ijk)/Ndiv2
                 else
                    !used for example for emissions in kg/m2 (or kg/m2/s)
                    ! integral is approximately conserved
                    Rvar(ijk)=Rvar(ijk)/Ivalues(ijk)
                    !if(me==12.and.i==10.and.j==2)write(*,*)'cdfVALUE: ',k,Rvar(ijk),Ivalues(ijk),k2,data3D
                 end if
              end if
           end do
        end do
        end do
        deallocate(Ivalues)
     else
       call CheckStop(interpol_used=='mass_conservative', "ReadField_CDF: only linear interpolation implemented")
       if(interpol_used=='zero_order'.and.MasterProc.and.debug)&
         write(*,*)'zero_order interpolation asked, but performing linear interpolation'

!      call CheckStop(data3D, "ReadField_CDF : 3D not yet implemented for general projection")
      call CheckStop(present(UnDef), "Default values filling not implemented")

      allocate(Weight(4,LIMAX,LJMAX),&
                 IIij(4,LIMAX,LJMAX),&
                 JJij(4,LIMAX,LJMAX))

      !Make interpolation coefficients.
      !Coefficients could be saved and reused if called several times.
      if(DEBUG_NETCDF_RF.and.debug_proc.and.i==debug_li.and.j==debug_lj)&
         write(*,"(a)") "DEBUG_RF G2G ", me, debug_proc
      call grid2grid_coeff(glon,glat,IIij,JJij,Weight,&
           Rlon,Rlat,dims(1),dims(2), LIMAX, LJMAX, limax, ljmax,&
           (DEBUG_NETCDF_RF.and.debug_proc), debug_li, debug_lj )

      startvec(1)=minval(IIij(:,:limax,:ljmax))
      startvec(2)=minval(JJij(:,:limax,:ljmax))
      if(ndims>2)startvec(ndims)=nstart
      dims=1
      dims(1)=maxval(IIij(:,:limax,:ljmax))-startvec(1)+1
      dims(2)=maxval(JJij(:,:limax,:ljmax))-startvec(2)+1

      k2=1
      if(data3D)then
         if(kstart_loc<0)then
            !take all levels
            kstart_loc=1
            kend_loc=dims(3)
         end if
         startvec(3)=kstart_loc
         k2=kend_loc-kstart_loc+1
         dims(3)=kend_loc-kstart_loc+1
         if(ndims>3)then
            startvec(ndims)=nstart
            dims(ndims)=1
         end if
      else
         if(ndims>2)then
            startvec(ndims)=nstart
            dims(ndims)=1
         end if
      end if
      
      totsize=1
      do i=1,ndims
        totsize=totsize*dims(i)
      end do

      allocate(Rvalues(totsize), stat=alloc_err)
      if(debug) then
        write(*,"(2a)") 'ReadCDF VarID ', trim(varname)
        do i=1, ndims
          write(*,"(a,6i8)") 'ReadCDF ',i, dims(i),startvec(i)
        end do
        write(*,*)'total size variable (part read only)',totsize
      end if

      call check(nf90_get_var(ncFileID, VarID, Rvalues,start=startvec,count=dims),&
            errmsg="RRvalues")

      if(xtype==NF90_INT.or.xtype==NF90_SHORT.or.xtype==NF90_BYTE)then
        !scale data if it is packed
        scalefactors(1) = 1.0 !default
        scalefactors(2) = 0.  !default
        status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
        if(status == nf90_noerr) scalefactors(1) = scale
        status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
        if(status == nf90_noerr) scalefactors(2) = offset
        Rvalues=Rvalues*scalefactors(1)+scalefactors(2)
      end if

      k=1
      do i=1,limax
        do j=1,ljmax
          ww => Weight(:,i,j)
          ijn(:)=IIij(:,i,j)-startvec(1)+1+(JJij(:,i,j)-startvec(2))*dims(1)

          do k=1,k2

!          ijk=i+(j-1)*LIMAX
              ij=i+(j-1)*LIMAX
              ijk=k+(ij-1)*k2
              Rvar(ijk)    = 0.0
              sumWeights   = 0.0


              do iw = 1, 4
                 ii = IIij(iw,i,j)
                 jj = JJij(iw,i,j)

                 if(Rvalues(ijn(iw))/=FillValue) then
                    Rvar(ijk)  = Rvar(ijk) + ww(iw)*Rvalues(ijn(iw)+(k-1)*dims(1)*dims(2))
                    sumWeights = sumWeights+ ww(iw)
                 end if
              end do !iw

          if(sumWeights>1.0e-9) then
            Rvar(ijk) =  Rvar(ijk)/sumWeights
          else
            Rvar(ijk) = FillValue
          end if

          end do !k
        end do !j
      end do !i

      deallocate(Weight,IIij,JJij)
    end if!conservative
  end if!general projection

!  if(interpolate_vertical)then
!do the interpolation in the vertical direction only
!filename is used only to define the vertical coordinates
!     call vertical_interpolate(filename,Rvar,dims(3),debug)

!  end if


  deallocate(Rvalues)
  deallocate(Rlon)
  deallocate(Rlat)
  if(fractions)then
     deallocate(NCC,CC,fraction_in)
  end if
  if(.not. present(ncFileID_given))call check(nf90_close(ncFileID))


  !  CALL MPI_FINALIZE(IERROR)
  !  CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)

  return

  !code below only used for testing purposes

  CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
  if(debug)write(*,*)'writing results in file. Variable: ',trim(varname)

  !only for tests:
  def1%class='Readtest' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=1.0        !not used
  !  def1%inst=.true.      !not used
  !  def1%year=.false.     !not used
  !  def1%month=.false.    !not used
  !  def1%day=.false.      !not used
  def1%name=trim(varname)!written
  def1%unit='g/m2'       !written

  if(data3D)then
     return
     k2=kend_loc-kstart_loc+1
     n=3

     call Out_netCDF(IOU_INST,def1,n,k2, &
          Rvar,1.0,CDFtype=Real4,fileName_given='ReadField3D.nc')

     CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
     CALL MPI_FINALIZE(IERROR)
     stop
  else
     if(trim(varname)=='nonHighwayRoadDustPM10_Jun-Feb')then
        !      if(.true.)then
        n=2
        k2=1

        call Out_netCDF(IOU_INST,def1,n,k2, &
             rvar,1.0,CDFtype=Real4,fileName_given='ReadField2D.nc',overwrite=.false.)
        CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
        !    CALL MPI_FINALIZE(IERROR)
        !     stop
     end if
  end if

  return

end subroutine ReadField_CDF

!The "flight levels" are so specific, that we prefer to have an ad-hoc routine rather
!than to pollute otehr routine with special cases.
 subroutine ReadField_CDF_FL(fileName,varname,Rvar,nstart,kstart,kend,interpol,needed,debug_flag)

  use netcdf

  implicit none
  character(len = *),intent(in) ::fileName,varname
  real,intent(out) :: Rvar(*)
  integer,intent(in) :: nstart
  integer, intent(in) :: kstart!smallest k (vertical level) to read. Default: assume 2D field
  integer, intent(in) :: kend!largest k to read. Default: assume 2D field
  character(len = *), optional,intent(in) :: interpol
  logical, optional, intent(in) :: needed
  logical, optional, intent(in) :: debug_flag

  integer :: ncFileID,VarID,lonVarID,latVarID,status,ndims,dimids(NF90_MAX_VAR_DIMS),xtype,nAtts
  integer :: dims(NF90_MAX_VAR_DIMS),totsize,i,j,k
  integer :: startvec(NF90_MAX_VAR_DIMS)
  integer ::alloc_err
  real :: dloni,dlati
  integer ::ij,jdiv,idiv,Ndiv,Ndiv2,igjgk,ig,jg,ijk
  integer ::imin,imax,jmin,jmax,igjg,k2
  integer, allocatable:: Ivalues(:)  ! I counts all data
  integer, allocatable:: Nvalues(:)  !ds counts all values
  real, allocatable:: Rvalues(:),Rlon(:),Rlat(:)
  real ::lat,lon,maxlon,minlon,maxlat,minlat
  logical ::fileneeded, debug,data3D
  character(len = TXTLEN_NAME) :: interpol_used, data_projection="",name
  real :: Grid_resolution
  integer, parameter ::NFL=23,NFLmax=50 !number of flight level (could be read from file)
  real :: P_FL(0:NFLmax),P_FL0,Psurf_ref(LIMAX, LJMAX),P_EMEP,dp!

  real :: Pcounted
  logical :: Flight_Levels
  integer :: k_FL,k_FL2

  !_______________________________________________________________________________
  !
  !1)           General checks and init
  !_______________________________________________________________________________

  fileneeded=.true.!default

  debug = .false.
  if(present(debug_flag))then
     debug = debug_flag .and. MasterProc
     if ( debug ) write(*,*) 'ReadCDF start: ',trim(filename),':', trim(varname)
  end if
  if(present(needed))   fileneeded=needed

  !open an existing netcdf dataset
  status=nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID)
  
  if(status == nf90_noerr) then
     if ( debug ) write(*,*) 'ReadCDF reading ',trim(filename), ' nstart ', nstart
  else
     if(fileneeded)then
        print *, 'file does not exist: ',trim(fileName),nf90_strerror(status)
        call CheckStop(fileneeded, "ReadField_CDF : file needed but not found")
     else
        if(MasterProc) write(*,*)'file does not exist (but not needed): ',&
             trim(fileName),nf90_strerror(status)
          Rvar(1:LIMAX*LJMAX*(KMAX_MID-2))=0.0
         return
     end if
  end if

  interpol_used='mass_conservative'!default for FL
  if(present(interpol))then
     interpol_used=interpol
  end if
  call CheckStop(interpol_used/='mass_conservative',&
       'interpolation method for FL not recognized')

  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = trim(varname), varID = VarID)
  if(status == nf90_noerr) then
    if(debug) write(*,*) 'ReadCDF variable exists: ',trim(varname)
  else
    if(fileneeded)then
      print *, 'variable does not exist: ',trim(varname),nf90_strerror(status)
      call CheckStop(fileneeded, "ReadField_CDF : variable needed but not found")

    else
        if(MasterProc)write(*,*) 'variable does not exist (but not needed): ',&
             trim(varname),nf90_strerror(status)
          call check(nf90_close(ncFileID))
          Rvar(1:LIMAX*LJMAX*(KMAX_MID-2))=0.0
        return
    end if
  end if

  data3D=.true.

  !Check first that variable has data covering the relevant part of the grid:

  !Find chunk of data required (local)
  maxlon=max(maxval(gl_stagg),maxval(glon))
  minlon=min(minval(gl_stagg),minval(glon))
  maxlat=max(maxval(gb_stagg),maxval(glat))
  minlat=min(minval(gb_stagg),minval(glat))


  !get dimensions id
  call check(nf90_Inquire_Variable(ncFileID,VarID,name,&
       xtype,ndims,dimids,nAtts),"GetDimsId")

  !get dimensions
  startvec=1
  dims=0
  do i=1,ndims
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(i), &
          len=dims(i)),"GetDims")
     if ( debug ) write(*,*) 'ReadCDF size variable ',i,dims(i)
  end do

  data_projection = "lon lat"

     allocate(Rlon(dims(1)), stat=alloc_err)
     allocate(Rlat(dims(2)), stat=alloc_err)


  call check(nf90_inq_varid(ncid = ncFileID, name='lon', varID = lonVarID))
  call check(nf90_get_var(ncFileID, lonVarID, Rlon), 'Getting Rlon')

  call check(nf90_inq_varid(ncid = ncFileID, name='lat', varID = latVarID))
  call check(nf90_get_var(ncFileID, latVarID, Rlat), 'Getting Rlat')

  !_______________________________________________________________________________
  !
  !2)        Coordinates conversion and interpolation
  !_______________________________________________________________________________


  if(trim(data_projection)=="lon lat")then
     call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(3), name=name ))
     call CheckStop(trim(name)/='FL','expected Flight Levels')
        !special vertical levels for Aircrafts
        !make table for conversion Flight Level -> Pressure
        !Hard coded because non-standard anyway. 610 meters layers
        do k=0,NFLmax
           P_FL(k)=1000*StandardAtmos_km_2_kPa(k*0.610)
        end do
        P_FL0=P_FL(0)
        Flight_Levels=.true.
        call CheckStop(interpol_used/='mass_conservative',&
             "only mass_conservative interpolation implemented for Flight Levels")
        
        !need average surface pressure for the current month
        !montly average is needed, not instantaneous pressure
        call ReadField_CDF(SurfacePressureFile,'surface_pressure',&
             Psurf_ref,current_date%month,needed=.true.,interpol='zero_order',debug_flag=debug_flag)
  end if


     !NB: we assume regular grid
     !inverse of resolution
     dloni=1.0/(Rlon(2)-Rlon(1))
     dlati=1.0/(Rlat(2)-Rlat(1))

     Grid_resolution = EARTH_RADIUS*abs(Rlat(2)-Rlat(1))*PI/180.0

     imin=mod( floor((minlon-Rlon(1))*dloni)+dims(1),dims(1))+1!NB lon  -90 = +270
     imax=mod(ceiling((maxlon-Rlon(1))*dloni)+dims(1),dims(1))+1!NB lon  -90 = +270
     if(imax==1)imax=dims(1)!covered entire circle
     if(minlon-Rlon(1)<0.0)then
        imin=1!cover entire circle
        imax=dims(1)!cover entire circle
     end if

     jmin=max(1,min(dims(2),floor((minlat-Rlat(1))*dlati)))
     jmax=max(1,min(dims(2),ceiling((maxlat-Rlat(1))*dlati)+1))
     
     if(maxlat>85.0.or.minlat<-85.0)then
        !close to poles
        imin=1
        imax=dims(1)
     end if

     if(imax<imin)then
        !crossing longitude border !
        !   write(*,*)'WARNING: crossing end of map'
        !take everything...could be memory expensive
        imin=1
        imax=dims(1)
     end if

     if ( debug ) write(*,"(a,4f8.2,6i8)") 'ReadCDF minmax values ',&
          minlon,maxlon,minlat,maxlat,imin,imax,jmin,jmax

     startvec(1)=imin
     startvec(2)=jmin
     dims=1
     dims(1)=imax-imin+1
     dims(2)=jmax-jmin+1


     startvec(3)=1
     dims(3)=NFL

     startvec(ndims)=nstart!time dimension

     totsize=1
     do i=1,ndims
        totsize=totsize*dims(i)
     end do

     allocate(Rvalues(totsize), stat=alloc_err)
     if ( debug ) then
        do i=1, ndims ! NF90_MAX_VAR_DIMS would be 1024
           write(*,"(a,6i8)") 'ReadCDF reading chunk ',i ,startvec(i), dims(i)
        end do
     end if
     call check(nf90_get_var(ncFileID, VarID, Rvalues,start=startvec,count=dims),&
          errmsg="RRvalues")

     if(interpol_used=='conservative'.or.interpol_used=='mass_conservative')then
        !conserves integral (almost, does not take into account local differences in mapping factor)
        !takes weighted average over gridcells covered by model gridcell

        !divide the coarse grid into pieces significantly smaller than the fine grid
        !Divide each global gridcell into Ndiv x Ndiv pieces
        Ndiv=nint(3*Grid_resolution/GRIDWIDTH_M)
        Ndiv=max(1,Ndiv)
        Ndiv2=Ndiv*Ndiv
        !

        k2=kend-kstart+1
        allocate(Ivalues(LIMAX*LJMAX*k2))
        allocate(Nvalues(LIMAX*LJMAX*k2))
        do ij=1,LIMAX*LJMAX*k2
           Ivalues(ij)=0
           NValues(ij) = 0
           Rvar(ij)=0.0
        end do

        do jg=1,dims(2)
           do jdiv=1,Ndiv
              lat=Rlat(startvec(2)-1+jg)-0.5/dlati+(jdiv-0.5)/(dlati*Ndiv)
              do ig=1,dims(1)
                 igjg=ig+(jg-1)*dims(1)
                 do idiv=1,Ndiv
                    lon=Rlon(startvec(1)-1+ig)-0.5/dloni+(idiv-0.5)/(dloni*Ndiv)
                    call lb2ij(lon,lat,i,j)
                    i=i-gi0-IRUNBEG+2
                    j=j-gj0-JRUNBEG+2

                    if(i>=1.and.i<=limax.and.j>=1.and.j<=ljmax)then
                      ij=i+(j-1)*LIMAX
                       k2=kend-kstart+1
                       if(Flight_Levels)then
                          !Flight_Levels
                          !start filling levels from surface and upwards
                          !add emissions at every emep OR FL level boundary
                          k_FL=1
                          k_FL2=0!last index of entirely included FL layer
                          P_FL(0)=max(A_bnd(KMAX_MID+1)+B_bnd(KMAX_MID+1)*Psurf_ref(i,j), P_FL(0))
                          Pcounted=P_FL(0)!Lowest Pressure accounted for
                          do k=KMAX_MID,KMAX_MID-k2+1,-1
                             ijk=k-(KMAX_MID-k2)+(ij-1)*k2
                             P_EMEP=A_bnd(k)+B_bnd(k)*Psurf_ref(i,j)
                             do while(P_FL(k_FL)>P_EMEP.and.k_FL<NFL)
                                dp=Pcounted-P_FL(k_FL)
                                igjgk=igjg+(k_FL-1)*dims(1)*dims(2)
                                Rvar(ijk)=Rvar(ijk)+Rvalues(igjgk)*dp/(P_FL(k_FL2)-P_FL(k_FL))
                                k_FL2=k_FL
                                k_FL=k_FL+1
                                Pcounted=P_FL(k_FL2)
                             end do
                             Ivalues(ijk)=Ivalues(ijk)+1
                             Nvalues(ijk)=Nvalues(ijk)+1
                             if(k_FL<=NFL)then
                                dp=Pcounted-P_EMEP
                                igjgk=igjg+(k_FL-1)*dims(1)*dims(2)
                                Rvar(ijk)=Rvar(ijk)+Rvalues(igjgk)*dp/(P_FL(k_FL2)-P_FL(k_FL))
                                Pcounted=P_EMEP
                             end if
                          end do

                          P_FL(0)=P_FL0!may have changed above
                       end if !Flight levels

                    end if
                 end do
              end do
           end do
        end do

        k2=kend-kstart+1
        do k=1,k2
           do i=1,limax
              do j=1,ljmax
                 ij=i+(j-1)*LIMAX
                 ijk=k+(ij-1)*k2


                 if(Ivalues(ijk)<=0.)then
                       write(*,"(a,a,4i4,6f10.3,2i6,6f10.3)") &
                            'ERROR, NetCDF_mod no values found!', &
                            trim(fileName) // ":" // trim(varname), &
                            i,j,k,me,minlon,maxlon,minlat,maxlat,glon(i,j),glat(i,j), &
                            Ivalues(ijk),Ndiv,Rlon(startvec(1)),&
                            Rlon(startvec(1)+dims(1)-1),Rlat(startvec(2)),&
                            Rlat(startvec(2)-1+dims(2))
!                       call CheckStop("Interpolation error")
                    !we simply set emission=0
                    Rvar(ijk)=0.0
                 else
                    if(interpol_used=='mass_conservative')then
                       !used for example for emissions in kg (or kg/s)
                       Rvar(ijk)=Rvar(ijk)/Ndiv2! Total sum of values from all cells is constant
                    else
                       call CheckStop("interpol choice not supported")
                    end if
                 end if
              end do
           end do
        end do

        deallocate(Ivalues)
        deallocate(Nvalues)

     end if! projection

  deallocate(Rvalues)
  deallocate(Rlon)
  deallocate(Rlat)

  call check(nf90_close(ncFileID))

end subroutine ReadField_CDF_FL

subroutine printCDF(name, array,unit)
    ! Minimal print out to cdf, for real numbers, 2-d arrays
    character(len=*), intent(in) :: name
    real, dimension(:,:), intent(in) :: array
    character(len=*), intent(in) :: unit

    character(len=60) :: fname
    type(Deriv) :: def1 ! definition of fields

    def1%class='print-cdf' !written
    def1%avg=.false.      !not used
    def1%index=0          !not used
    def1%scale=1.0      !not used
    def1%name=trim(name)   ! written
    def1%unit=trim(unit)

    fname = "PRINTCDF_" // trim(name) // ".nc"

    !Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype,ist,jst,ien,jen,ik,fileName_given)
    write(*,*)me ,' printcdf'
    if(MasterProc) write(*,*) "OUTPUTS printCDF :"//trim(fname),  maxval(array)
    call Out_netCDF(IOU_INST,def1,2,1, array,1.0,&
           CDFtype=Real4,fileName_given=fname,overwrite=.true.)
  end subroutine printCDF


subroutine ReadTimeCDF(filename,TimesInDays,NTime_Read,time_wanted)
  !Read times in file under CF convention and convert into days since 1900-01-01 00:00:00
  !OR in years since 2000 if the times are defined in years!
  !if present(time_wanted), find corresponding record and return it in NTime_Read
  character(len=*) ,intent(in):: filename
  real,intent(out) :: TimesInDays(:)
  integer, intent(inout) :: NTime_Read ! in:records to read, out:records readed
  real,optional,intent(in) :: time_wanted!if present, find first record after this time (within 0.5 second difference)  

  real, allocatable :: times(:)
  integer, allocatable :: int_times(:,:,:)
  integer :: i,j,ntimes,status
  integer :: varID,ncFileID,ndims
  integer :: xtype,dimids(NF90_MAX_VAR_DIMS),nAtts
  integer, parameter::wordarraysize=20
  character(len=TXTLEN_NAME) :: varname,period,since,name,timeunit,wordarray(wordarraysize),calendar
  character(len=TXTLEN_NAME) :: wordarray2(wordarraysize)
  character, allocatable :: Times_string(:,:)
  integer :: string_length
  integer :: yyyy,mo,dd,hh,mi,ss,julian,julian_1900,diff_1900,nwords,errcode,startrecord
  logical:: proleptic_gregorian, find_record

  call check(nf90_open(path=fileName, mode=nf90_nowrite, ncid=ncFileID),&
       errmsg="ReadTimeCDF, file not found: "//trim(fileName))

  find_record = .false.
  if(present(time_wanted)) find_record = .true.!only find the record needed
  if(find_record) NTime_Read = 1

  varname='time'
  status=nf90_inq_varid(ncid=ncFileID, name=varname, varID=VarID)

  if(status==nf90_noerr)then
     if(DEBUG_NETCDF.and.MasterProc)write(*,*)'time variable exists: ',trim(varname)     
     call check(nf90_Inquire_Variable(ncFileID,VarID,name,xtype,ndims,dimids,nAtts))
     if(ndims>1)write(*,*)'WARNING: time has more than 1 dimension!? ',ndims
     call check(nf90_inquire_dimension(ncid=ncFileID,dimID=dimids(1),len=ntimes))
     if(NTime_Read<1)then
        if(DEBUG_NETCDF.and.MasterProc)write(*,*)'reading all time records'
        NTime_Read=ntimes
     end if
     call CheckStop(ntimes<NTime_Read .and. .not.find_record, "to few records in "//trim(fileName))
     call CheckStop(SIZE(TimesInDays)<NTime_Read,"to many records in "//trim(fileName))
     allocate(times(NTime_Read))
     
     if(.not. find_record)ntimes=1

     do startrecord = 1,ntimes !in case find_record, loop over one record at a time

     call check(nf90_get_var(ncFileID, VarID, times, start = (/startrecord/),count=(/NTime_Read/)))
     
     call check(nf90_get_att(ncFileID, VarID, "units", timeunit  ))
     
     if(trim(timeunit(1:16))=='day as %Y%m%d.%f')then
        if(DEBUG_NETCDF.and.MasterProc)write(*,*)'time unit ',trim(timeunit(1:16))
        do i=1,NTime_Read
           yyyy=int(times(i))/10000
           mo=int(times(i)-10000*yyyy)/100
           dd=int(times(i)-10000*yyyy-100*mo)
           hh=24.0*(times(i)-int(times(i)))
           julian=julian_date(yyyy,mo,dd)
           julian_1900=julian_date(1900,1,1)
           diff_1900=julian-julian_1900
           TimesInDays(i)=diff_1900+hh/24.0
           if(DEBUG_NETCDF.and.MasterProc)write(*,*)'time ',yyyy,mo,dd,hh
        end do
     else if (trim(timeunit)=='month')then
        !used for files with only 12 monthly values, without year
        forall(i=1:NTime_Read) &
             TimesInDays(i) = 30*(times(i)-1)+15
     else
        !must be of the form " xxxx since yyyy-mm-dd hh:mm:ss"
        
        !    read(timeunit,fmt="(a,a,a,a)")period,since,date,time
        call wordsplit(trim(timeunit),wordarraysize,wordarray,nwords,errcode,separator='-')
        if(DEBUG_NETCDF.and.MasterProc)&
             write(*,*)"time@units:",(" ",trim(wordarray(i)),i=1,8)
        period=wordarray(1)
        since=wordarray(2)
        call CheckStop(since/='since',"since error "//trim(since))
        
        read(wordarray(3),*)yyyy
        read(wordarray(4),*)mo
        read(wordarray(5),*)dd
        read(wordarray(6),*)hh
        if( period == 'minutes' .or. period =='seconds' )then
           read(wordarray(7),*)mi
           read(wordarray(8),*)ss  !did not work for others?
        else
           mi=0
           ss=0
        endif

        calendar='unknown'
        
        status=nf90_get_att(ncFileID, VarID, "calendar", calendar )
        proleptic_gregorian=(status==nf90_noerr).and.(calendar=='proleptic_gregorian'.or.calendar=='gregorian')
        if(proleptic_gregorian.and.DEBUG_NETCDF.and.MasterProc)&
             write(*,*)'found proleptic_gregorian calendar'
        
        if(yyyy/=0.or.proleptic_gregorian)then
           !read(date,fmt="(I4.4,a1,I2.2,a1,I2.2)")yyyy,s1,mo,s2,dd
           !read(time,fmt="(I2.2,a1,I2.2,a1,I2.2)")hh,s1,mi,s2,ss
           
           if(DEBUG_NETCDF.and.MasterProc)&
                write(*,"(A,I4.4,2('-',I2.2),A,I2.2,2(':',I2.2))")&
                ' refdate ',yyyy,mo,dd,' time ',hh,mi,ss
           ss=ss+60*mi+3600*hh
           julian=julian_date(yyyy,mo,dd)
           julian_1900=julian_date(1900,1,1)
           diff_1900=julian-julian_1900
           !if(MasterProc)write(*,*)'julians ',diff_1900,julian,julian_1900
           select case(period)
           case('years')
              !NB: not completely "standard" and compatible with the "days" setup
              forall(i=1:NTime_Read) &
                   TimesInDays(i)=times(i)-(2000-yyyy)
           case('days')
              forall(i=1:NTime_Read) &
                   TimesInDays(i)=diff_1900+times(i)+ss/(3600.0*24.0)
           case('hours')
              forall(i=1:NTime_Read) &
                   TimesInDays(i)=diff_1900+(times(i)+ss/3600.0)/24.0
           case('minutes')
              read(wordarray(7),*)mi
              forall(i=1:NTime_Read) &
                   TimesInDays(i)=diff_1900+(times(i)*60.0+ss)/(3600.0*24.0)
           case('seconds')
              forall(i=1:NTime_Read) &
                   TimesInDays(i)=diff_1900+(times(i)+ss)/(3600.0*24.0)
           case('month')!used for files with only 12 monthly values, without year
              forall(i=1:NTime_Read) &
                   TimesInDays(i) = 30*(times(i)-1)+15
           case default
              call StopAll("ReadTimeCDF, time unit not recognized: "//trim(period))
           end select
           
        else
           if(DEBUG_NETCDF.and.MasterProc)&
                write(*,*)'assuming days since 0-01-01 00:00 and 365days'
           call CheckStop(period/='days',"Error: only time in days implemented "//trim(period))
           !assume units = "days since 0-01-01 00:00"
           !and calendar = "365_day"
           yyyy=int(times(1)/365)
           
           julian=julian_date(yyyy,1,1)
           julian_1900=julian_date(1900,1,1)
           diff_1900=julian-julian_1900
           forall(i=1:NTime_Read) &
                TimesInDays(i)=diff_1900+times(i)-yyyy*365
           
           !for leap years and dates after 28th February add one day to get Julian days
           if(mod(yyyy,4)==0)then
              forall(i=1:NTime_Read,times(i)-yyyy*365>59.999) & ! later than midnight on 
                   TimesInDays(i)=TimesInDays(i)+1.0              ! Feb 28th (59th day)
              ! if the current date in the model is 29th of february, then this date
              ! is not defined in the 365 days calendar.
              ! We then assume that the 60th day is 29th of february in the netcdf file
              ! and not the 1st of march.
              ! Keep this separately as this may be defined differently in different situations.
              ! This implementation works for the IFS-MOZART BC
              if(current_date%month==2.and.current_date%day==29)then
                 write(*,*)'WARNING: assuming 29th of February for ',trim(filename)
                 forall(i=1:NTime_Read,int(times(i)-yyyy*365)==60) & ! move Mar 1st
                      TimesInDays(i)=TimesInDays(i)-1.0                 ! to Feb 29th
              end if
           end if
        end if
     end if
     if(find_record)then
        if(TimesInDays(1)-time_wanted>-0.5/(24.0*3600.0))then
           !we have found the right record
           NTime_Read = startrecord
           exit
        endif
     endif
     enddo
     deallocate(times)

     call CheckStop(startrecord>ntimes.and.find_record, "did not find correct time in "//trim(fileName))

  else
!     write(*,*)'ReadTimeCDF '//trim(varname)//" not found in "//trim(fileName)
     varname='Times'!wrf format
     status=nf90_inq_varid(ncid=ncFileID, name=varname, varID=VarID)
     if(status==nf90_noerr)then
        call CheckStop(find_record,"find_record not implemented for time format in "//trim(fileName))
        call check(nf90_Inquire_Variable(ncFileID,VarID,name,xtype,ndims,dimids,nAtts))
        if(ndims>2)write(*,*)'WARNING: Times has more than 2 dimension!? ',ndims
        call check(nf90_inquire_dimension(ncid=ncFileID,dimID=dimids(2),len=ntimes))
        if(NTime_Read<1)then
           if(DEBUG_NETCDF)write(*,*)'reading all time records'
           NTime_Read=ntimes
        end if
        call CheckStop(ntimes<NTime_Read, "to few records in "//trim(fileName))
        call CheckStop(SIZE(TimesInDays)<NTime_Read,"to many records in "//trim(fileName))
        call check(nf90_inquire_dimension(ncid=ncFileID,dimID=dimids(1),len=string_length))
        
        allocate(Times_string(string_length,ntimes))
        call check(nf90_get_var(ncFileID, VarID, Times_string(1:string_length,1:NTime_Read),count=(/string_length,NTime_Read/)))
        do i=1,NTime_Read
           do j=1,string_length
              name(j:j)=Times_string(j,i)
           end do
           call wordsplit(name,string_length,wordarray,nwords,errcode,'_')
           if(DEBUG_NETCDF.and.MasterProc)write(*,*)'date ',trim(wordarray(1)),' hour ',trim(wordarray(2)),' minutes ',trim(wordarray(3))
           call wordsplit(wordarray(1),wordarraysize,wordarray2,nwords,errcode,'-')
           if(DEBUG_NETCDF.and.MasterProc)write(*,*)'year ',trim(wordarray2(1)),' month ',trim(wordarray2(2)),' day ',trim(wordarray2(3))
           read(wordarray2(1),*)yyyy
           read(wordarray2(2),*)mo
           read(wordarray2(3),*)dd
           read(wordarray(2),*)hh !colon, ":", is a default separator for wordsplit
           read(wordarray(3),*)mi !colon, ":", is a default separator for wordsplit
           julian=julian_date(yyyy,mo,dd)
           julian_1900=julian_date(1900,1,1)
           diff_1900=julian-julian_1900
           TimesInDays(i)=diff_1900+hh/24.0+mi/24.0/60.0
        end do
     else
        varname='TFLAG'!SMOKE/CMAQ format
        status=nf90_inq_varid(ncid=ncFileID, name=varname, varID=VarID)
        if(status==nf90_noerr)then
           call check(nf90_Inquire_Variable(ncFileID,VarID,name,xtype,ndims,dimids,nAtts))
           if(ndims/=3)write(*,*)'WARNING: TSTEP does not have 3 dimensions!? ',ndims
           call check(nf90_inquire_dimension(ncid=ncFileID,dimID=dimids(3),len=ntimes))
           if(NTime_Read<1)then
              if(DEBUG_NETCDF)write(*,*)'reading all time records'
              NTime_Read=ntimes
           end if
           call CheckStop(ntimes<NTime_Read, "to few records in "//trim(fileName))
           call CheckStop(SIZE(TimesInDays)<NTime_Read,"to many records in "//trim(fileName))
           
           allocate(int_times(2,1,NTime_Read))
           !NB: we assume all variables have same time stamp. Read only for first

           if(.not. find_record)ntimes=1

           do startrecord = 1,ntimes !in case find_record, loop over one record at a time

              call check(nf90_get_var(ncFileID, VarID, int_times,start = (/1,1,startrecord/),count=(/2,1,NTime_Read/)))
              do i=1,NTime_Read
                 yyyy=int_times(1,1,i)/1000
                 mo=1!start at beginning of year
                 dd=mod(int_times(1,1,i),1000)!nb day of year!!
                 hh= int_times(2,1,i)/10000 + &
                      1.0/60*(mod(int_times(2,1,i),10000)/100) + &
                      1.0/3600*mod(int_times(2,1,i),100)
                 julian=julian_date(yyyy,mo,dd)
                 julian_1900=julian_date(1900,1,1)
                 diff_1900=julian-julian_1900
!                 if(me==0)write(*,*)'TFLAG time ',int_times(1,1,i),int_times(2,1,i),dd,hh
                 TimesInDays(i)=diff_1900+hh/24.0
              end do
              if(find_record)then
                 if(TimesInDays(1)-time_wanted>-0.5/(24.0*3600.0))then
                    !we have found the right record        
                    NTime_Read = startrecord
                    exit
                 endif
                 if(startrecord == ntimes .and. me==0)write(*,*)'WARNING: did not find correct emis time. Last time found ',TimesInDays(1),' wanted ',time_wanted
              endif
           enddo
           deallocate(int_times)
        else
           if(DEBUG_NETCDF)write(*,*)'time variable not found: ',trim(varname)
           NTime_Read = 0
        endif
     endif
  endif
  
  call check(nf90_close(ncFileID))

end subroutine ReadTimeCDF

subroutine   vertical_interpolate(filename,Rvar,KMAX_ext,Rvar_emep,debug)
  character(len = *),intent(in) ::fileName
  integer, intent(in)::KMAX_ext
  real,intent(in) :: Rvar(KMAX_ext,LIMAX,LJMAX)
  real,intent(out):: Rvar_emep(LIMAX,LJMAX,KMAX_MID)
  logical, optional, intent(in) :: debug!output only masterproc
  real, allocatable,dimension(:)::hyam_ext,hybm_ext,P_ext,weight_k1
  integer, allocatable,dimension(:)::k1_ext,k2_ext
  integer :: ncFileID,varID,status,i,j,k,k_ext!,ij,ijk
  real :: P_emep,P0
  character(len=100) :: word
  logical ::reversed_k

  allocate(k1_ext(KMAX_MID),k2_ext(KMAX_MID),weight_k1(KMAX_MID))
 !Read pressure for vertical levels
!  if(MasterProc)write(*,*)'vertical_interpolate: reading vertical levels'
!  if(MasterProc)write(*,*)'filename ',trim(filename)
  allocate(P_ext(KMAX_ext),hyam_ext(KMAX_ext+1),hybm_ext(KMAX_ext+1))
  call check(nf90_open(fileName,nf90_nowrite,ncFileID),&
       errmsg="ReadTimeCDF, file not found: "//trim(fileName))
 
  status=nf90_inq_varid(ncFileID,"hyam",varID)
  if(status==nf90_noerr) then
    if(debug) write(*,*)'Found hyam type levels (values at level midpoints)'
    call check(nf90_get_var(ncFileID,varID,hyam_ext(1:KMAX_ext)))
    status=nf90_get_att(ncFileID,VarID,"units",word)
    if(status==nf90_noerr)then
      if(word(1:3)=='hPa')then
        if(MasterProc) write(*,*)'Changing hyam from hPa to Pa in '//trim(filename)
        hyam_ext=100*hyam_ext
      end if
    end if
    call check(nf90_inq_varid(ncFileID,"hybm",varID))
    call check(nf90_get_var(ncFileID,varID,hybm_ext(1:KMAX_ext)))
  else
   status=nf90_inq_varid(ncFileID,"hyai",varID)
   if(status==nf90_noerr) then
      if(debug)write(*,*)'Found hyai type levels (values at level interfaces)'          
      call check(nf90_get_var(ncFileID,varID,hyam_ext))
      status=nf90_get_att(ncFileID,VarID,"units",word)
      if(status==nf90_noerr)then
        if(word(1:3)=='hPa')then
          if(debug)write(*,*)'Changing hyai from hPa to Pa'
          hyam_ext=100*hyam_ext
        end if
      end if
      call check(nf90_inq_varid(ncFileID,"hybi",varID))
      call check(nf90_get_var(ncFileID,varID,hybm_ext))
      do k=1,KMAX_ext
        hyam_ext(k)=0.5*(hyam_ext(k)+hyam_ext(k+1))
        hybm_ext(k)=0.5*(hybm_ext(k)+hybm_ext(k+1))
      end do
    else
      call StopAll('levels not yet implemented')
    end if
  end if
  call check(nf90_inq_varid(ncFileID,"lev",varID))                 
  call check(nf90_get_att(ncFileID,VarID,"formula_terms",word))
  if(word(1:2)/='ap')then
     status=nf90_inq_varid(ncFileID,"P0",varID)                 
     if(status==nf90_noerr) then
        call check(nf90_get_var(ncFileID,varID,P0))
        status=nf90_get_att(ncFileID,VarID,"units",word)
        if(status==nf90_noerr)then
           if(word(1:3)=='hPa')then
              if(debug)write(*,*)'Changing P0 from hPa to Pa'
              P0=100*P0
           end if
        end if
        if(debug)write(*,*)'Multiplying hyam with P0 ',P0
        hyam_ext=hyam_ext*P0
     else
        if(debug)write(*,*)'did not find P0 ',status
     end if
  end if

  !find vertical interpolation coefficients
  !use pressure as reference
  !we want, if possible, P_ext(k1) and P_ext(k2) to be on each side of P_emep
  !We assume constant surface pressure, both for emep and external grid; should not be so
  !   important as long as they are both terrain following.
  do k_ext=1,KMAX_EXT
    P_ext(k_ext)=hyam_ext(k_ext)+hybm_ext(k_ext)*Pref
    if(debug) write(*,fmt="(A,I3,F10.2)")'vert_inter: P_ext',k_ext,P_ext(k_ext)
  end do

  reversed_k=(P_ext(1)>P_ext(2))
  ! .true.  --> assumes k_ext=KMAX_EXT is top and k_ext=1 is surface
  ! .false. --> assumes k_ext=1 is top and k_ext=KMAX_EXT is surface  
  if(reversed_k)then
    do k=1,KMAX_MID
      P_emep=A_mid(k)+B_mid(k)*Pref !Pa
      if(debug) write(*,fmt="(A,I3,F10.2)")'vert_inter: P_emep',k,P_emep
      if(P_emep<P_ext(KMAX_EXT))then
         !we need levels above the highest level available. 
         !we do not want to extrapolate. take the top level value
         weight_k1(k)=1.0
         k1_ext(k)=KMAX_EXT
         k2_ext(k)=KMAX_EXT
      else if(P_emep>P_ext(1))then
         !we need levels below the lowest level available. 
         !we do not want to extrapolate. take the lowest level value
         weight_k1(k)=1.0
         k1_ext(k)=1
         k2_ext(k)=1
      else
         !largest available P smaller than P_emep (if possible)
         k1_ext(k)=1 !start at surface, and go up until P_emep
         do k_ext=1,KMAX_EXT
            if(P_ext(k_ext)<P_emep)exit
            k1_ext(k)=k_ext
         end do
         !smallest available P larger than P_emep (if possible)
         k2_ext(k)=KMAX_EXT !start at top, and go down until P_emep
         if(k2_ext(k)==k1_ext(k))k2_ext(k)=KMAX_EXT-1 !avoid k2=k1
         do k_ext=KMAX_EXT,1,-1
            if(P_ext(k_ext)>P_emep)exit
            if(k_ext/=k1_ext(k))k2_ext(k)=k_ext
         end do
         weight_k1(k)=(P_emep-P_ext(k2_ext(k)))/(P_ext(k1_ext(k))-P_ext(k2_ext(k)))
      end if
      if(debug)&
        write(*,fmt="(A,I4,2(A,I4,A,F5.2))")'vert_inter: level',k,&
          ' is the sum of level ',k1_ext(k),' weight ',weight_k1(k),&
          ' and level ',k2_ext(k),' weight ',1-weight_k1(k)
    end do     
  else
    do k=1,KMAX_MID
      P_emep=A_mid(k)+B_mid(k)*Pref !Pa
      if(debug) write(*,fmt="(A,I3,F10.2)")'vert_inter: P_emep',k,P_emep
      if(P_emep<P_ext(1))then
         !we need levels above the highest level available. 
         !we do not want to extrapolate. take the top level value
         weight_k1(k)=1.0
         k1_ext(k)=1
         k2_ext(k)=1
      else if(P_emep>P_ext(KMAX_EXT))then
         !we need levels below the lowest level available. 
         !we do not want to extrapolate. take the lowest level value
         weight_k1(k)=1.0
         k1_ext(k)=KMAX_EXT
         k2_ext(k)=KMAX_EXT
      else
         !largest available P smaller than P_emep (if possible)
         k1_ext(k)=KMAX_EXT !start at surface, and go up until P_emep
         do k_ext=KMAX_EXT,1,-1
            if(P_ext(k_ext)<P_emep)exit
            k1_ext(k)=k_ext
         end do
         !smallest available P larger than P_emep (if possible)
         k2_ext(k)=1 !start at top, and go down until P_emep
         if(k2_ext(k)==k1_ext(k))k2_ext(k)=2 !avoid k2=k1
         do k_ext=1,KMAX_EXT
            if(P_ext(k_ext)>P_emep)exit
            if(k_ext/=k1_ext(k))k2_ext(k)=k_ext
         end do
         weight_k1(k)=(P_emep-P_ext(k2_ext(k)))/(P_ext(k1_ext(k))-P_ext(k2_ext(k)))
      end if
      if(debug) &
        write(*,fmt="(A,I4,2(A,I4,A,F5.2))")'vertical_interpolate: level',k,&
          ' is the sum of level ', k1_ext(k),' weight ',weight_k1(k),&
          ' and level ', k2_ext(k),' weight ',1.0-weight_k1(k)
    end do
  end if

  forall(i=1:limax,j=1:ljmax,k=1:KMAX_MID)&
    Rvar_emep(i,j,k)=weight_k1(k) *Rvar(k1_ext(k),i,j)&
               +(1.0-weight_k1(k))*Rvar(k2_ext(k),i,j)

  deallocate(P_ext,hyam_ext,hybm_ext,k1_ext,k2_ext,weight_k1)
  call check(nf90_close(ncFileID))
end subroutine vertical_interpolate

 function IsCDFfractionFormat(filename) result(foundFractions)

   character (len=*), intent(in):: filename
   logical ::foundFractions
   integer :: cdfstatus,VarID,ncFileID

   foundFractions=.false.
   cdfstatus=nf90_open(path = trim(filename), mode = nf90_nowrite, ncid = ncFileID)
   if(cdfstatus == nf90_noerr)then
      cdfstatus = nf90_inq_varid(ncid = ncFileID, name = 'NCodes', varID = VarID)
      call check(nf90_close(ncFileID))
      if(cdfstatus == nf90_noerr)foundFractions=.true.
   end if

 end function IsCDFfractionFormat

 subroutine ReadSectorName(filename,cdf_sector_name)

   character (len=*), intent(in):: filename
   character (len=*), intent(inout):: cdf_sector_name
   integer :: cdfstatus,ncFileID
   character (len=len(cdf_sector_name)):: sector_name

   cdfstatus=nf90_open(path = trim(filename), mode = nf90_nowrite, ncid = ncFileID)
   if(cdfstatus == nf90_noerr)then
      cdfstatus=nf90_get_att(ncFileID,nf90_global,'SECTORS_NAME',sector_name)
      call check(nf90_close(ncFileID))
      if(cdfstatus == nf90_noerr)cdf_sector_name=trim(sector_name)
   end if

 end subroutine ReadSectorName

 subroutine readfrac(Ncc,CC,Rvalues,fraction_in,fractions_out,Ncc_out,CC_out,val,dim1dim2,igjgk,ijk,latlon_weight,Reduc)
!accumulates the values val, and update fractions
   implicit none
   real, optional,intent(in) :: Reduc(NLAND)
   integer, intent(in) :: dim1dim2,igjgk,ijk
   real, intent(inout) :: fractions_out(LIMAX*LJMAX,*),val
   real, intent(in) :: Rvalues,fraction_in(dim1dim2,*),latlon_weight
   integer, intent(inout)  ::Ncc_out(*), CC_out(LIMAX*LJMAX,*)
   integer, intent(in)  ::Ncc,CC(dim1dim2,*)
   integer :: N,Ng,N_out,ic
   real :: factor, total

   do Ng=1,Ncc!number of fields at igjg as read
      do N_out=1,Ncc_out(ijk) !number of fields at ij already saved in the model grid
         if(CC(igjgk,Ng)==CC_out(ijk,N_out))goto 737
      end do
      !the country is not yet used for this gridcell. Define it now
      Ncc_out(ijk)=Ncc_out(ijk)+1
      N_out=Ncc_out(ijk)
      CC_out(ijk, N_out)=CC(igjgk,Ng)
      fractions_out(ijk,N_out)=0.0
737   continue
      factor=1.0!default reduction factor
      !if(present(Reduc).and.CC(igjgk,Ng)>0.and.CC(igjgk,Ng)<=NLAND)factor=Reduc(CC(igjgk,Ng))                               
      if(present(Reduc).and.CC(igjgk,Ng)>0)then
         ic=find_index(CC(igjgk,Ng),Country(:)%icode)
         if(ic>NLAND.or.ic<1)then
            write(*,*)"ReadField_cdf: COUNTRY CODE NOT RECOGNIZED OR UNDEFINED: ",&
                 CC(igjgk,Ng)
            call StopAll("COUNTRY CODE NOT RECOGNIZED ")
         end if
         factor=Reduc(ic)
      end if
      !update fractions
      total=val+Rvalues*fraction_in(igjgk,Ng)*factor*latlon_weight
!      if(debug.and.fraction_in(igjgk,Ng)>1.001)then
!         write(*,*)'fractions_in TOO LARGE ',Ng,ig,jg,k,fraction_in(igjgk,Ng)
!         stop
!      end if
      if(abs(total)>1.0E-30)then
         do N=1,Ncc_out(ijk)
            !reduce previously defined fractions
            fractions_out(ijk,N)=fractions_out(ijk,N)*val/total
         end do
         !increase fraction of this country (yes, after having reduced it!)
         fractions_out(ijk,N_out)=fractions_out(ijk,N_out)+Rvalues*fraction_in(igjgk,Ng)*latlon_weight/total*factor
      else
         !should try to keep proportions right in case cancellation of positive an negative; not finished!
         do N=1,Ncc_out(ijk)
            !reduce existing fractions
            fractions_out(ijk,N)=fractions_out(ijk,N)/Ncc_out(ijk)
         end do
         !increase fraction of this country (yes, after having reduced it!)
         fractions_out(ijk,N_out)=fractions_out(ijk,N_out)+Rvalues*fraction_in(igjgk,Ng)*latlon_weight/Ncc_out(ijk)*factor
      end if
      val=total
   end do
 end subroutine readfrac

 subroutine check_lon_lat(ncFileID, lon_name, lat_name, ndims, lonVarID, latVarID)
   !check that longitude and latitude for every gridcell is defined as a 2 dimensional array
   !file is must be open already
   character(len=*), intent(inout)  :: lon_name, lat_name
   integer, intent(in)  :: ncFileID
   integer, intent(out),optional  ::  lonVarID, latVarID
   integer, intent(out)  :: ndims
   integer :: status,VarID,len(4)

   lat_name = 'lat'
   lon_name = 'lon'
   status=nf90_inq_varid(ncid = ncFileID, name=trim(lon_name), varID = VarID)
   if(status /= nf90_noerr) then
      lon_name = 'LON'
      status=nf90_inq_varid(ncid = ncFileID, name = lon_name, varID = VarID)
      if(status /= nf90_noerr) then
         lon_name = 'longitude'
         status=nf90_inq_varid(ncid = ncFileID, name = lon_name, varID = VarID)
         call CheckStop(status /= nf90_noerr,'did not find longitude variable')
      end if
   end if
   if(present(lonVarID))lonVarID=VarID
   call check(nf90_Inquire_Variable(ncid = ncFileID,  varID = VarID,ndims=ndims))
       
   status=nf90_inq_varid(ncid = ncFileID, name=trim(lat_name), varID = VarID)
   if(status /= nf90_noerr) then
      lat_name = 'LAT'
      status=nf90_inq_varid(ncid = ncFileID, name = lat_name, varID = VarID)
      if(status /= nf90_noerr) then
         lat_name = 'latitude'
         status=nf90_inq_varid(ncid = ncFileID, name = lat_name, varID = VarID)
         call CheckStop(status /= nf90_noerr,'did not find latitude variable')
     end if
   end if
   if(present(latVarID))latVarID=VarID

 end subroutine check_lon_lat

 !this routine returns an approximate grid resolution in meters
 !file must be open already, and is not closed
 subroutine make_gridresolution(ncFileID,resolution)
   integer, intent(in) :: ncFileID
   real, intent(out)::resolution
   real :: dlat
   real, allocatable ::Rlat(:,:)
   character(len=30)  :: lon_name, lat_name
   integer :: nDimensions, xtype, lonVarID, latVarID, dimids(NF90_MAX_VAR_DIMS),dims(10),mindim,i,ndims,nAtts

   call check_lon_lat(ncFileID, lon_name, lat_name, nDimensions, lonVarID, latVarID)
   call check(nf90_Inquire_Variable(ncFileID,latVarID,lat_name,ndims,xtype,dimids,nAtts),"EmisGetDimsId")

   if(nDimensions==1)then
      allocate(Rlat(2,1))
      call check(nf90_get_var(ncFileID, latVarID, Rlat,count=(/2/)))
      dlat = abs(Rlat(2,1)-Rlat(1,1))
   else   if(nDimensions==2)then
      call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(1),len=dims(1)),"EmisGetDims2D")
      call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(2),len=dims(2)),"EmisGetDims2D")
      
      allocate(Rlat(dims(1),dims(2)))
      call check(nf90_get_var(ncFileID, latVarID, Rlat))
      dlat=0.0
      mindim = min(dims(1),dims(2))
      do i = 2,mindim
         !check only along both diagonal
         if(abs(Rlat(i-1,i-1)-Rlat(i,i))>dlat)dlat=abs(Rlat(i-1,i-1)-Rlat(i,i))
         if(abs(Rlat(mindim-i+2,i-1)-Rlat(mindim-i+1,i))>dlat)dlat=abs(Rlat(mindim-i+2,i-1)-Rlat(mindim-i+1,i))            
      enddo
   else
      call CheckStop('did not find one or two dimensional latitudes '//trim(lat_name))
   endif
   deallocate(Rlat)
   resolution = EARTH_RADIUS*dlat*PI/180.0          
   if(me==0)write(*,*)'Will set default resolution (it does not need to be exact) ',resolution
 end subroutine make_gridresolution

 subroutine output_country_emissions(filename,GridEmis,GridEmisCodes,nGridEmisCodes,NSECTORS,NCMAX,EMIS_FILE,NEMIS_FILE)
   !output emissions in a format readable by the code
   implicit none
   character(len=*),  intent(in)  :: filename,EMIS_FILE(NEMIS_FILE)
   integer,  intent(in) :: NSECTORS,NCMAX,NEMIS_FILE
   real, intent(in) ::GridEmis(NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILE)
   integer,  intent(in) :: GridEmisCodes(LIMAX,LJMAX,NCMAX),nGridEmisCodes(LIMAX,LJMAX)
   logical :: country_owner_map(NLAND,NPROC), country_owner(NLAND), country_has_emis
   integer :: varID,iDimID,jDimID,timeDimID,status,ncFileID
   integer :: istart,jstart,icount,jcount
   integer :: d,i,j,ij,k, isec, icc, iland, iem
   real :: buff(MAXLIMAX*MAXLJMAX)
   real, allocatable :: emis_full(:,:) 
   character (len = TXTLEN_NAME) ::varname

   call create_country_emission_file(fileName)

   !make a map that show which CPU has non-zero emissions from each country
   country_owner = .false.
   country_owner_map = .false.
   do j = 1,ljmax
      do i = 1,limax
         do icc = 1,nGridEmisCodes(i,j)
            iland=find_index(GridEmisCodes(i,j,icc),Country(:)%icode)
            country_owner(iland) = .true.
         end do
      end do
   end do
   call MPI_GATHER(country_owner,NLAND,MPI_LOGICAL,country_owner_map,NLAND,MPI_LOGICAL,0,MPI_COMM_CALC, IERROR)
   
   if(MasterProc)then
      call check(nf90_open(trim(fileName),nf90_share+nf90_write,ncFileID))

       allocate(emis_full(GIMAX,GJMAX))
      !first create all variables
       do iland = 1, NLAND
         country_has_emis = .false.
         do d = 1, NPROC
            if(country_owner_map(iland,d))country_has_emis = .true.
         enddo
         if(.not. country_has_emis) cycle
         
         do isec = 1, NSECTORS
            do iem = 1, NEMIS_File
               44 FORMAT(A,I0,A)
               write(varname,44)trim(EMIS_FILE(iem))//'_sec',isec,'_'//trim(Country(iland)%code)
               status=nf90_inq_varid(ncFileID,varname,VarID)
               if(status/=nf90_noerr)then
                  write(*,*)me,' creating ',trim(varname),MasterProc
                  call check(nf90_inq_dimid(ncFileID,"i"  ,idimID),"dim:i")
                  call check(nf90_inq_dimid(ncFileID,"j"  ,jdimID),"dim:j")
                  call check(nf90_inq_dimid(ncFileID,"time"  ,timedimID),"dim:time")
                  call check(nf90_def_var(ncFileID,varname,nf90_float,&
                       [iDimID,jDimID,timeDimID]       ,varID),"def2d:"//trim(varname))
                  call check(nf90_def_var_deflate(ncFileid,varID,shuffle=0,deflate=1,&
                       deflate_level=4),"compress:"//trim(varname))
                  call check(nf90_def_var_chunking(ncFileID,varID,NF90_CHUNKED,&
                       (/min(GIMAX,300),min(GJMAX,130),1/)),"chunk2D:"//trim(varname))    
                  
                  call check(nf90_put_att(ncFileID,varID,"units",'kg'))
                  call check(nf90_put_att(ncFileID,varID,"sector",isec))
                  call check(nf90_put_att(ncFileID,varID,"species",trim(EMIS_File(iem))))
               end if
            end do
         end do
      end do
   endif
   
   !gather data from other CPU into master
   do iland = 1, NLAND
      if(MasterProc)then
         country_has_emis = .false.
         do d = 1, NPROC
            if(country_owner_map(iland,d))country_has_emis = .true.
         enddo
         if(.not. country_has_emis) cycle
      else
         if(.not. country_owner(iland)) cycle
      endif
      
      do iem = 1, NEMIS_File
         do isec = 1, NSECTORS
            if(MasterProc)then
               write(varname,44)trim(EMIS_FILE(iem))//'_sec',isec,'_'//trim(Country(iland)%code)
               emis_full = 0.0 !important that default is zero
               if(country_owner(iland))then
                  do j = 1, ljmax
                     do i = 1, limax
                        do icc = 1,nGridEmisCodes(i,j)
                           if(GridEmisCodes(i,j,icc)==iland)then
                              emis_full(i,j) = emis_full(i,j) + GridEmis(isec,i,j,icc,iem)
                           endif
                        end do
                     end do
                  end do
               endif
               do d = 1, NPROC-1
                  if(country_owner_map(iland,d+1))then
                     CALL MPI_RECV(buff, 8*tlimax(d)*tljmax(d), MPI_BYTE, d, &
                          isec+iem*NSECTORS+iland*NEMIS_File*NSECTORS,& !unique tag
                          MPI_COMM_CALC, MPISTATUS, IERROR)
                     ij = 0
                     do j = tgj0(d),tgj0(d)+tljmax(d)-1
                        do i = tgi0(d),tgi0(d)+tlimax(d)-1
                           ij = ij+1
                           emis_full(i,j) = buff(ij)
                        end do
                     end do
                  endif
               end do
               
               call check(nf90_inq_varid(ncFileID,trim(varname),varID),"inq:"//trim(varname))
               call check(nf90_put_var(ncFileID,VarID,emis_full))
               
            else 
               ! not master

               ij = 0
               do j=1, ljmax
                  do i=1, limax
                     ij = ij+1
                     buff(ij)=0.0
                     do icc = 1,nGridEmisCodes(i,j)
                        if(GridEmisCodes(i,j,icc)==iland)then
                           buff(ij)=buff(ij)+GridEmis(isec,i,j,icc,iem)
                        endif
                     enddo
                  end do
               end do
               CALL MPI_SEND( buff, 8*tlimax(me)*tljmax(me), MPI_BYTE, 0, &
                    isec+iem*NSECTORS+iland*NEMIS_File*NSECTORS,& !unique tag
                    MPI_COMM_CALC, IERROR)
            endif
         end do
      end do
   end do
   if(MasterProc)then
      deallocate(emis_full)
      call check(nf90_close(ncFileID))
   end if
  
 end subroutine output_country_emissions

 subroutine create_country_emission_file(fileName)
   implicit none
   character(len=*),  intent(in)  :: filename
   integer :: iDimID,jDimID,timeDimID,ncFileID,VarID, longVarID,latVarID
   real :: Buff2D(MAXLIMAX,MAXLJMAX,2)
   integer :: istart,jstart,icount,jcount
   integer :: i,j,k,iproc

   ! the file does not exist yet or is overwritten
   if(MasterProc) then 
      if(DEBUG_NETCDF)write(*,*)'nf90_create '//trim(fileName)
      call check(nf90_create(trim(fileName),nf90_hdf5,ncFileID),"create:"//trim(fileName))
      
      ! define coordinate dimensions
      call check(nf90_def_dim(ncFileID,"i",GIMAX,iDimID),"dim:i")
      call check(nf90_def_dim(ncFileID,"j",GJMAX,jDimID),"dim:j")
      call check(nf90_def_dim(ncFileID,"time",nf90_unlimited,timeDimID))
      
      ! define coordinate variables
      call check(nf90_def_var(ncFileID,"lat"   ,nf90_float,[iDimID,jDimID],latvarID))
      call check(nf90_def_var(ncFileID,"lon"   ,nf90_float,[iDimID,jDimID],longvarID))
      
      call check(nf90_put_att(ncFileID,nf90_global,"projection",'generic'))
      call check(nf90_put_att(ncFileID,nf90_global,"Grid_resolution",GRIDWIDTH_M))
      call check(nf90_close(ncFileID))
      
   endif
      
   Buff2D(1:limax,1:ljmax,1)=glon(1:limax,1:ljmax)
   Buff2D(1:limax,1:ljmax,2)=glat(1:limax,1:ljmax)
   if(MasterProc)then
      call check(nf90_open(fileName,nf90_share+nf90_write,ncFileID))
      do iproc=0,NPROC-1
         if(iproc>0)&
              CALL MPI_RECV(Buff2D,2*8*MAXLIMAX*MAXLJMAX, MPI_BYTE,iproc,&
              iproc,MPI_COMM_CALC,MPISTATUS,IERROR)
         !NB IBEGcdf, IRUNBEG, etc. relative to fulldomain
         !    tgi0,tlimax, etc. relative to rundomain
         
        istart=tgi0(iproc)    ! restricted (IBEGcdf, JBEGcdf) domain
        jstart=tgj0(iproc)    ! restricted domain
        icount=tlimax(iproc)  ! both 
        jcount=tljmax(iproc)  ! both 

        call check(nf90_put_var(ncFileID,longVarID,   &
             Buff2D(1:icount,1:jcount,1), &
             start=[istart,jstart],count=[icount,jcount]))
        call check(nf90_put_var(ncFileID, latVarID,   &
             Buff2D(1:icount,1:jcount,2), &
             start=[istart,jstart],count=[icount,jcount]))
     end do
     call check(nf90_close(ncFileID))
  else
     CALL MPI_SEND(Buff2D,2*8*MAXLIMAX*MAXLJMAX,MPI_BYTE,0,&
          me,MPI_COMM_CALC,IERROR)
  end if
    
end subroutine create_country_emission_file
endmodule NetCDF_mod
