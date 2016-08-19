! <NetCDF_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

       module NetCDF_ml
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
  use My_Outputs_ml,     only : FREQ_HOURLY, &
                                NHOURLY_OUT, &      ! No. outputs
                                Asc2D, hr_out, &      ! Required outputs
                                SELECT_LEVELS_HOURLY, LEVELS_HOURLY, &
                                NLEVELS_HOURLY 
  use Chemfields_ml,     only : xn_shl,xn_adv
  use CheckStop_ml,      only : CheckStop,StopAll
  use ChemSpecs_shl_ml,  only : NSPEC_SHL
  use ChemSpecs_adv_ml,  only : NSPEC_ADV
  use ChemSpecs_tot_ml,  only : NSPEC_TOT
  use ChemChemicals_ml,  only : species
  use GridValues_ml,     only : GRIDWIDTH_M,fi,xp,yp,xp_EMEP_official&
                               ,debug_proc, debug_li, debug_lj &
                               ,yp_EMEP_official,fi_EMEP,GRIDWIDTH_M_EMEP&
                               ,grid_north_pole_latitude&
                               ,grid_north_pole_longitude,GlobalPosition&
                               ,glat_fdom,glon_fdom,ref_latitude&
                               ,projection, sigma_mid,gb_stagg,gl_stagg,glon&
                               ,sigma_bnd&
                               ,glat,lb2ij,A_bnd,B_bnd
  use InterpolationRoutines_ml,  only : grid2grid_coeff
  use ModelConstants_ml, only : KMAX_MID,KMAX_BND, runlabel1, runlabel2 &
                               ,MasterProc &
                               ,DEBUG_NETCDF, DEBUG_NETCDF_RF &
                               ,NPROC, IIFULLDOM,JJFULLDOM &
                               ,IOU_INST,IOU_HOUR,IOU_HOUR_MEAN, IOU_YEAR &
                               ,IOU_MON, IOU_DAY ,PT,NLANDUSEMAX, model
  use netcdf
  use OwnDataTypes_ml,   only : Deriv
  use Par_ml,            only : me,GIMAX,GJMAX,tgi0,tgj0,tlimax,tljmax, &
                        MAXLIMAX, MAXLJMAX,IRUNBEG,JRUNBEG,limax,ljmax,gi0,gj0
  use PhysicalConstants_ml,  only : PI, EARTH_RADIUS
  use TimeDate_ml,       only: nmdays,leapyear ,current_date, date,julian_date
  use TimeDate_ExtraUtil_ml,only: idate2nctime
  use Functions_ml,       only: StandardAtmos_km_2_kPa
  use SmallUtils_ml,      only: wordsplit

  implicit none


  INCLUDE 'mpif.h'
  INTEGER MPISTATUS(MPI_STATUS_SIZE),INFO

  character (len=125), save :: fileName_inst = 'out_inst.nc'
  character (len=125), save :: fileName_hour = 'out_hour.nc'
  character (len=125), save :: fileName_day = 'out_day.nc'
  character (len=125), save :: fileName_month = 'out_month.nc'
  character (len=125), save :: fileName_year = 'out_year.nc'
  character (len=125) :: fileName ,period_type

  integer,parameter ::closedID=-999     !flag for showing that a file is closed
  integer      :: ncFileID_new=closedID  !don't save because should always be 
                  !redefined (in case several routines are using ncFileID_new 
                  !with different filename_given)
  integer,save :: ncFileID_inst=closedID
  integer,save :: ncFileID_hour=closedID
  integer,save :: ncFileID_day=closedID
  integer,save :: ncFileID_month=closedID
  integer,save :: ncFileID_year=closedID
  integer,save :: outCDFtag=0
  !CDF types for output:
  integer, public, parameter  :: Int1=1,Int2=2,Int4=3,Real4=4,Real8=5 
  character (len=18),parameter::Default_projection_name = 'General_Projection'

  public :: Out_netCDF
  public :: printCDF   ! minimal caller for Out_netCDF
  public :: CloseNetCDF
  public :: Init_new_netCDF
  public :: GetCDF
  public :: WriteCDF !for testing purposes
  public :: Read_Inter_CDF
  public :: ReadField_CDF
  public :: ReadTimeCDF

  private :: CreatenetCDFfile
  private :: createnewvariable
  private :: check

contains
!_______________________________________________________________________


subroutine Init_new_netCDF(fileName,iotyp)
integer,  intent(in) :: iotyp
character(len=*),  intent(in)  :: fileName

integer :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
integer :: ih

call CloseNetCDF !must be called by all procs, to syncronize outCDFtag

if(MasterProc)then
if( DEBUG_NETCDF ) write(*,*)'Init_new_netCDF ',trim(fileName),iotyp

select case (iotyp)
case (IOU_YEAR)
  fileName_year = trim(fileName)
  period_type = 'fullrun'
  if( DEBUG_NETCDF ) write(*,*) "Creating ", trim(fileName),trim(period_type)
  call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)
case(IOU_MON)
  fileName_month = trim(fileName)
  period_type = 'monthly'
  if( DEBUG_NETCDF ) write(*,*) "Creating ", trim(fileName),trim(period_type)
  call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)
case(IOU_DAY)
  fileName_day = trim(fileName)
  period_type = 'daily'
  if( DEBUG_NETCDF ) write(*,*) "Creating ", trim(fileName),trim(period_type)
  call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)
case(IOU_HOUR)
  fileName_hour = trim(fileName)
  period_type = 'hourly'
  if( DEBUG_NETCDF ) write(*,*) "Creating ", trim(fileName),trim(period_type)
  ISMBEGcdf=GIMAX+IRUNBEG-1; JSMBEGcdf=GJMAX+JRUNBEG-1  !initialisations
  GIMAXcdf=0; GJMAXcdf=0                                !initialisations
  KMAXcdf=1
  do ih=1,NHOURLY_OUT
    ISMBEGcdf=min(ISMBEGcdf,hr_out(ih)%ix1)
    JSMBEGcdf=min(JSMBEGcdf,hr_out(ih)%iy1)
    GIMAXcdf=max(GIMAXcdf,hr_out(ih)%ix2-hr_out(ih)%ix1+1)
    GJMAXcdf=max(GJMAXcdf,hr_out(ih)%iy2-hr_out(ih)%iy1+1)
    KMAXcdf =max(KMAXcdf,hr_out(ih)%nk)
  enddo
  GIMAXcdf=min(GIMAXcdf,GIMAX)
  GJMAXcdf=min(GJMAXcdf,GJMAX)
  KMAXcdf =min(KMAXcdf ,NLEVELS_HOURLY) 

! Output selected model levels
  if(SELECT_LEVELS_HOURLY)then
    call CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,&
                          KMAXcdf,KLEVcdf=LEVELS_HOURLY)
  else
  if( DEBUG_NETCDF ) write(*,*) "Creating ", trim(fileName),trim(period_type)
    call CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,&
                          KMAXcdf)
  endif

case(IOU_INST)
  fileName_inst = trim(fileName)
  period_type = 'instant'
  if( DEBUG_NETCDF ) write(*,*) "Creating ", trim(fileName),trim(period_type)
  call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)
case default
  period_type = 'unknown'
  if( DEBUG_NETCDF ) write(*,*) "Creating ", trim(fileName),trim(period_type)
  call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)
end select

if( DEBUG_NETCDF ) write(*,*) "Finished Init_new_netCDF", trim(fileName),&
                  trim(period_type)
endif

end subroutine Init_new_netCDF

! Output selected model levels
subroutine CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,&
                            KMAXcdf,KLEVcdf,RequiredProjection)

integer, intent(in) :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
character(len=*),  intent(in)  :: fileName
character(len=*),optional, intent(in):: requiredprojection
integer, intent(in), optional :: KLEVcdf(KMAXcdf)

character(len=*), parameter :: author_of_run='Unimod group'
character(len=*), parameter :: vert_coord='sigma: k ps: PS ptop: PT'
character(len=19) :: projection_params='90.0 -32.0 0.933013' !set later on

real :: xcoord(GIMAX),ycoord(GJMAX),kcoord(KMAX_MID)

character(len=8)  :: created_date,lastmodified_date
character(len=10) :: created_hour,lastmodified_hour
integer :: iDimID,jDimID,kDimID,timeDimID,VarID,iVarID,jVarID,kVarID,i,j,k
integer :: ncFileID,iEMEPVarID,jEMEPVarID,latVarID,longVarID,PTVarID
real :: scale_at_projection_origin
character(len=80) ::UsedProjection

  ! fileName: Name of the new created file
  ! nf90_clobber: protect existing datasets
  ! ncFileID: netcdf ID

!Check that the dimensions are > 0
  if(GIMAXcdf<=0.or.GJMAXcdf<=0.or.KMAXcdf<=0)then
    write(*,*)'WARNING:'
    write(*,*)trim(fileName),&
              ' not created. Requested area too small (or outside domain) '
    write(*,*)'sizes (IMAX,JMAX,IBEG,JBEG,KMAX) ',&
              GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
    return
  endif

  if(present(RequiredProjection))then
     UsedProjection=trim(RequiredProjection)
  else
     UsedProjection=trim(projection)
  endif

  write(*,*)'creating ',trim(fileName)
  if(DEBUG_NETCDF)write(*,*)'UsedProjection ',trim(UsedProjection)
  if(DEBUG_NETCDF)write(*,fmt='(A,8I7)')'with sizes (IMAX,JMAX,IBEG,JBEG,KMAX) ',&
    GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
  call check(nf90_create(path = trim(fileName), &
        cmode = nf90_clobber, ncid = ncFileID),"create:"//trim(fileName))

  ! Define the dimensions
  if(UsedProjection=='Stereographic')then
    call check(nf90_def_dim(ncid = ncFileID, name = "i", len = GIMAXcdf, dimid = iDimID))
    call check(nf90_def_dim(ncid = ncFileID, name = "j", len = GJMAXcdf, dimid = jDimID))

  elseif(UsedProjection=='lon lat')then
    call check(nf90_def_dim(ncid = ncFileID, name = "lon", len = GIMAXcdf, dimid = iDimID))
    call check(nf90_def_var(ncFileID, "lon", nf90_double, dimids = iDimID, varID = iVarID) )
    call check(nf90_put_att(ncFileID, iVarID, "standard_name", "longitude"))
    call check(nf90_put_att(ncFileID, iVarID, "long_name", "longitude"))
    call check(nf90_put_att(ncFileID, iVarID, "units", "degrees_east"))
    call check(nf90_def_dim(ncid = ncFileID, name = "lat", len = GJMAXcdf, dimid = jDimID))
    call check(nf90_def_var(ncFileID, "lat", nf90_double, dimids = jDimID, varID =jVarID) )
    call check(nf90_put_att(ncFileID, jVarID, "standard_name", "latitude"))
    call check(nf90_put_att(ncFileID, jVarID, "long_name", "latitude"))
    call check(nf90_put_att(ncFileID, jVarID, "units", "degrees_north"))

  else !general projection
    call check(nf90_def_dim(ncid = ncFileID, name = "i", len = GIMAXcdf, dimid = iDimID))
    call check(nf90_def_dim(ncid = ncFileID, name = "j", len = GJMAXcdf, dimid = jDimID))
    call check(nf90_def_var(ncFileID, "i", nf90_float, dimids = iDimID, varID = iVarID) )
    call check(nf90_put_att(ncFileID, iVarID, "standard_name", "projection_x_coordinate"))
    call check(nf90_put_att(ncFileID, iVarID, "coord_axis", "x"))
    call check(nf90_put_att(ncFileID, iVarID, "long_name", "grid x coordinate"))
    call check(nf90_put_att(ncFileID, iVarID, "units", "km"))
    call check(nf90_def_var(ncFileID, "j", nf90_float, dimids = jDimID, varID = jVarID) )
    call check(nf90_put_att(ncFileID, jVarID, "standard_name", "projection_y_coordinate"))
    call check(nf90_put_att(ncFileID, jVarID, "coord_axis", "y"))
    call check(nf90_put_att(ncFileID, jVarID, "long_name", "grid y coordinate"))
    call check(nf90_put_att(ncFileID, jVarID, "units", "km"))

    call check(nf90_def_var(ncFileID, "lat", nf90_float, dimids = (/ iDimID, jDimID/), varID = latVarID) )
    call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"))
    call check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"))
    call check(nf90_put_att(ncFileID, latVarID, "standard_name", "latitude"))

    call check(nf90_def_var(ncFileID, "lon", nf90_float, dimids = (/ iDimID, jDimID/), varID = longVarID) )
    call check(nf90_put_att(ncFileID, longVarID, "long_name", "longitude"))
    call check(nf90_put_att(ncFileID, longVarID, "units", "degrees_east"))
    call check(nf90_put_att(ncFileID, longVarID, "standard_name", "longitude"))
  endif

  call check(nf90_def_dim(ncid = ncFileID, name = "k", len = KMAXcdf, dimid = kDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "time", len = nf90_unlimited, dimid = timeDimID))

  call Date_And_Time(date=created_date,time=created_hour)
  if(DEBUG_NETCDF)write(6,*) 'created_date: ',created_date
  if(DEBUG_NETCDF)write(6,*) 'created_hour: ',created_hour

  ! Write global attributes
  call check(nf90_put_att(ncFileID, nf90_global, "Conventions", "CF-1.0" ))
 !call check(nf90_put_att(ncFileID, nf90_global, "version", version ))
  call check(nf90_put_att(ncFileID, nf90_global, "model", model))
  call check(nf90_put_att(ncFileID, nf90_global, "author_of_run", author_of_run))
  call check(nf90_put_att(ncFileID, nf90_global, "created_date", created_date))
  call check(nf90_put_att(ncFileID, nf90_global, "created_hour", created_hour))
  lastmodified_date = created_date
  lastmodified_hour = created_hour
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_date", lastmodified_date))
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour))

  call check(nf90_put_att(ncFileID, nf90_global, "projection",UsedProjection))

  if(UsedProjection=='Stereographic')then
    scale_at_projection_origin=(1.+sin(ref_latitude*PI/180.))/2.
    write(projection_params,fmt='(''90.0 '',F5.1,F9.6)')fi,scale_at_projection_origin
    call check(nf90_put_att(ncFileID, nf90_global, "projection_params",projection_params))

! define coordinate variables
    call check(nf90_def_var(ncFileID, "i", nf90_float, dimids = iDimID, varID = iVarID) )
    call check(nf90_put_att(ncFileID, iVarID, "standard_name", "projection_x_coordinate"))
    call check(nf90_put_att(ncFileID, iVarID, "coord_axis", "x"))
    call check(nf90_put_att(ncFileID, iVarID, "long_name", "EMEP grid x coordinate"))
    call check(nf90_put_att(ncFileID, iVarID, "units", "km"))

    call check(nf90_def_var(ncFileID, "i_EMEP", nf90_float, dimids = iDimID, varID = iEMEPVarID) )
    call check(nf90_put_att(ncFileID, iEMEPVarID, "long_name", "official EMEP grid coordinate i"))
    call check(nf90_put_att(ncFileID, iEMEPVarID, "units", "gridcells"))

    call check(nf90_def_var(ncFileID, "j", nf90_float, dimids = jDimID, varID = jVarID) )
    call check(nf90_put_att(ncFileID, jVarID, "standard_name", "projection_y_coordinate"))
    call check(nf90_put_att(ncFileID, jVarID, "coord_axis", "y"))
    call check(nf90_put_att(ncFileID, jVarID, "long_name", "EMEP grid y coordinate"))
    call check(nf90_put_att(ncFileID, jVarID, "units", "km"))

    call check(nf90_def_var(ncFileID, "j_EMEP", nf90_float, dimids = jDimID, varID = jEMEPVarID) )
    call check(nf90_put_att(ncFileID, jEMEPVarID, "long_name", "official EMEP grid coordinate j"))
    call check(nf90_put_att(ncFileID, jEMEPVarID, "units", "gridcells"))

    call check(nf90_def_var(ncFileID, "lat", nf90_float, dimids = (/ iDimID, jDimID/), varID = latVarID) )
    call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"))
    call check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"))
    call check(nf90_put_att(ncFileID, latVarID, "standard_name", "latitude"))

    call check(nf90_def_var(ncFileID, "lon", nf90_float, dimids = (/ iDimID, jDimID/), varID = longVarID) )
    call check(nf90_put_att(ncFileID, longVarID, "long_name", "longitude"))
    call check(nf90_put_att(ncFileID, longVarID, "units", "degrees_east"))
    call check(nf90_put_att(ncFileID, longVarID, "standard_name", "longitude"))
  endif

! call check(nf90_put_att(ncFileID, nf90_global, "vert_coord", vert_coord))
  call check(nf90_put_att(ncFileID, nf90_global, "period_type", &
        trim(period_type)), ":period_type"//trim(period_type) )
  call check(nf90_put_att(ncFileID, nf90_global, "run_label", trim(runlabel2)))

  call check(nf90_def_var(ncFileID, "k", nf90_float, dimids = kDimID, varID = kVarID) )
  call check(nf90_put_att(ncFileID, kVarID, "coord_alias", "level"))
!pwsvs for CF-1.0
  call check(nf90_put_att(ncFileID, kVarID, "standard_name", "atmosphere_sigma_coordinate"))
  call check(nf90_put_att(ncFileID, kVarID, "formula_terms", trim(vert_coord)))
  call check(nf90_put_att(ncFileID, kVarID, "units", "sigma_level"))
  call check(nf90_put_att(ncFileID, kVarID, "positive", "down"))
  call check(nf90_def_var(ncFileID, "PT", nf90_float,  varID = PTVarID) )
  call check(nf90_put_att(ncFileID, PTVarID, "units", "Pa"))
  call check(nf90_put_att(ncFileID, PTVarID, "long_name", "Pressure at top"))

  call check(nf90_def_var(ncFileID, "time", nf90_double, dimids = timeDimID, varID = VarID) )
  if(trim(period_type) /= 'instant'.and.trim(period_type) /= 'unknown'.and.&
     trim(period_type) /= 'hourly' .and.trim(period_type) /= 'fullrun')then 
    call check(nf90_put_att(ncFileID, VarID, "long_name", "time at middle of period"))
  else
    call check(nf90_put_att(ncFileID, VarID, "long_name", "time at end of period"))
  endif
  call check(nf90_put_att(ncFileID, VarID, "units", "days since 1900-1-1 0:0:0"))


!CF-1.0 definitions:
  if(UsedProjection=='Stereographic')then
    call check(nf90_def_var(ncid = ncFileID, name = "Polar_Stereographic", xtype = nf90_int, varID=varID ) )
    call check(nf90_put_att(ncFileID, VarID, "grid_mapping_name", "polar_stereographic"))
    call check(nf90_put_att(ncFileID, VarID, "straight_vertical_longitude_from_pole", Fi))
    call check(nf90_put_att(ncFileID, VarID, "latitude_of_projection_origin", 90.0))
    call check(nf90_put_att(ncFileID, VarID, "scale_factor_at_projection_origin", scale_at_projection_origin))
  elseif(UsedProjection=='lon lat')then

  elseif(UsedProjection=='Rotated_Spherical')then
    call check(nf90_def_var(ncid = ncFileID, name = "Rotated_Spherical", xtype = nf90_int, varID=varID ) )
    call check(nf90_put_att(ncFileID, VarID, "grid_mapping_name", "rotated_latitude_longitude"))
    call check(nf90_put_att(ncFileID, VarID, "grid_north_pole_latitude", grid_north_pole_latitude))
    call check(nf90_put_att(ncFileID, VarID, "grid_north_pole_longitude", grid_north_pole_longitude))
  else
    call check(nf90_def_var(ncid = ncFileID, name = Default_projection_name, xtype = nf90_int, varID=varID ) )
    call check(nf90_put_att(ncFileID, VarID, "grid_mapping_name", trim(UsedProjection)))
  endif

  ! Leave define mode
  call check(nf90_enddef(ncFileID), "define_done"//trim(fileName) )

  call check(nf90_open(path = trim(fileName), mode = nf90_write, ncid = ncFileID))

! Define horizontal distances

  if(UsedProjection=='Stereographic')then
    xcoord(1)=(ISMBEGcdf-xp)*GRIDWIDTH_M/1000.
    do i=2,GIMAXcdf
      xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
    enddo
    call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )

    ycoord(1)=(JSMBEGcdf-yp)*GRIDWIDTH_M/1000.
    do j=2,GJMAXcdf
      ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
    enddo
    call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAXcdf)) )

! Define horizontal coordinates in the official EMEP grid
   !xp_EMEP_official=8.
   !yp_EMEP_official=110.
   !GRIDWIDTH_M_EMEP=50000.
   !fi_EMEP=-32.
    if(fi==fi_EMEP)then
      ! Implemented only if fi = fi_EMEP = -32 (Otherwise needs a 2-dimensional mapping)
      ! uses (i-xp)*GRIDWIDTH_M = (i_EMEP-xp_EMEP)*GRIDWIDTH_M_EMEP
      do i=1,GIMAXcdf
        xcoord(i)=(i+ISMBEGcdf-1-xp)*GRIDWIDTH_M/GRIDWIDTH_M_EMEP + xp_EMEP_official
       !print *, i,xcoord(i)
      enddo
      do j=1,GJMAXcdf
        ycoord(j)=(j+JSMBEGcdf-1-yp)*GRIDWIDTH_M/GRIDWIDTH_M_EMEP + yp_EMEP_official
       !print *, j,ycoord(j)
      enddo
    else
      do i=1,GIMAXcdf
        xcoord(i)=NF90_FILL_FLOAT
      enddo
      do j=1,GJMAXcdf
        ycoord(j)=NF90_FILL_FLOAT
      enddo
    endif
    call check(nf90_put_var(ncFileID, iEMEPVarID, xcoord(1:GIMAXcdf)) )
    call check(nf90_put_var(ncFileID, jEMEPVarID, ycoord(1:GJMAXcdf)) )

    if(DEBUG_NETCDF) write(*,*) "NetCDF: Starting long/lat defs"
    !Define longitude and latitude
    call GlobalPosition !because this may not yet be done if old version of meteo is used
    if(ISMBEGcdf+GIMAXcdf-1<=IIFULLDOM .and. JSMBEGcdf+GJMAXcdf-1<=JJFULLDOM)then
      call check(nf90_put_var(ncFileID, latVarID, &
        glat_fdom(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)))
      call check(nf90_put_var(ncFileID, longVarID, &
        glon_fdom(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)))
    endif

  elseif(UsedProjection=='lon lat') then
    do i=1,GIMAXcdf
      xcoord(i)= glon_fdom(i+ISMBEGcdf-1,1)
    enddo
    do j=1,GJMAXcdf
      ycoord(j)= glat_fdom(1,j+JSMBEGcdf-1)
    enddo
    call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )
    call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAXcdf)) )
  else
    xcoord(1)=(ISMBEGcdf-0.5)*GRIDWIDTH_M/1000.
    do i=2,GIMAXcdf
      xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
     !print *, i,xcoord(i)
    enddo
    call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )

    ycoord(1)=(JSMBEGcdf-0.5)*GRIDWIDTH_M/1000.
    do j=2,GJMAXcdf
      ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
    enddo
    call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )
    call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAXcdf)) )
   !write(*,*)'coord written'

   !Define longitude and latitude
   if(ISMBEGcdf+GIMAXcdf-1<=IIFULLDOM .and. JSMBEGcdf+GJMAXcdf-1<=JJFULLDOM)then
     call check(nf90_put_var(ncFileID, latVarID, &
       glat_fdom(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)))
     call check(nf90_put_var(ncFileID, longVarID, &
       glon_fdom(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)))
     endif
  endif
  if(DEBUG_NETCDF) write(*,*) "NetCDF: lon lat written"

  !Define vertical levels
  if(present(KLEVcdf))then     !order is defined in KLEVcdf
    do k=1,KMAXcdf
      if(KLEVcdf(k)==0)then
        kcoord(k)=sigma_bnd(KMAX_BND)              !0-->surface
      else
        kcoord(k)=sigma_mid(KMAX_MID-KLEVcdf(k)+1) !1-->20;2-->19;...;20-->1
      endif
    enddo
  elseif(KMAXcdf==KMAX_MID)then
    do k=1,KMAX_MID
      kcoord(k)=sigma_mid(k)
    enddo
  else
    do k=1,KMAXcdf
      kcoord(k)=sigma_mid(KMAX_MID-k+1) !REVERSE order of k !
    enddo
  endif
  call check(nf90_put_var(ncFileID, kVarID, kcoord(1:KMAXcdf)) )
  call check(nf90_put_var(ncFileID, PTVarID, PT ))

  call check(nf90_close(ncFileID))
  if(DEBUG_NETCDF)write(*,*)'NetCDF: file created, end of CreatenetCDFfile'
end subroutine CreatenetCDFfile

!_______________________________________________________________________

subroutine Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype,ist,jst,ien,jen,ik,fileName_given,overwrite)

  !The use of fileName_given is probably slower than the implicit filename used by defining iotyp.


  integer ,intent(in) :: ndim,kmax
  type(Deriv),     intent(in) :: def1 ! definition of fields
  integer,                         intent(in) :: iotyp
  real    ,intent(in) :: scale
  real, dimension(MAXLIMAX,MAXLJMAX,KMAX), intent(in) :: dat ! Data arrays
  integer, optional, intent(in) :: ist,jst,ien,jen,ik !start and end of saved area. 
                                                      !Only level ik is written if defined
  integer, optional, intent(in) :: CDFtype != OUTtype. (Integer*1, Integer*2,Integer*4, real*8 or real*4)
  character (len=*),optional, intent(in):: fileName_given!filename to which the data must be written
  logical, optional, intent(in) :: overwrite !overwrite if file already exists (in case fileName_given)
  !NB if the file fileName_given exist (also from earlier runs) it will be appended

  character(len=len(def1%name)) :: varname
  character*8 ::lastmodified_date
  character*10 ::lastmodified_hour,lastmodified_hour0,created_hour
  integer :: varID,nrecords,ncFileID=closedID
  integer :: ndate(4)
  integer :: info,d,alloc_err,ijk,status,i,j,k
  integer :: i1,i2,j1,j2
  real :: buff(MAXLIMAX*MAXLJMAX*KMAX_MID)
  real*8 , allocatable,dimension(:,:,:)  :: R8data3D
  real*4 , allocatable,dimension(:,:,:)  :: R4data3D
  integer*4, allocatable,dimension(:,:,:)  :: Idata3D
  integer :: OUTtype !local version of CDFtype
  integer :: iotyp_new
  integer :: iDimID,jDimID,kDimID,timeVarID
  integer :: GIMAX_old,GJMAX_old,KMAX_old
  integer :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf
  real*8  :: rdays,rdays_time(1)
  logical :: overwrite_local

  i1=1;i2=GIMAX;j1=1;j2=GJMAX  !start and end of saved area
  if(present(ist))i1=max(ist-IRUNBEG+1,i1)
  if(present(ien))i2=min(ien-IRUNBEG+1,i2)
  if(present(jst))j1=max(jst-JRUNBEG+1,j1)
  if(present(jen))j2=min(jen-JRUNBEG+1,j2)

  !Check that that the area is larger than 0
  if((i2-i1)<0.or.(j2-j1)<0.or.kmax<=0)return

  !make variable name
  write(varname,fmt='(A)')trim(def1%name)
  if(DEBUG_NETCDF) write(*,*)'Out_NetCDF: START ' , me, trim(varname)

  !to shorten the output we can save only the components explicitely named here
  !if(varname.ne.'D2_NO2'.and.varname.ne.'D2_O3' &
  !                         .and.varname.ne.'D2_PM10')return

  !do not write 3D fields (in order to shorten outputs)
  !if(ndim==3)return

  iotyp_new=0
  if(present(fileName_given))then
     !NB if the file already exist (also from earlier runs) it will be appended
     overwrite_local=.false.
     if(present(overwrite))overwrite_local=overwrite
     if(me==0)then
        if(DEBUG_NETCDF) write(*,*)'Out_NetCDF: fileName_given ' , me, trim(fileName_given)
        !try to open the file
        status=nf90_open(path = trim(fileName_given), mode = nf90_write, ncid = ncFileID)
        ISMBEGcdf=IRUNBEG+i1-1
        JSMBEGcdf=JRUNBEG+j1-1
        GIMAXcdf=i2-i1+1
        GJMAXcdf=j2-j1+1
        if(status /= nf90_noerr .or. overwrite_local) then !the file does not exist yet or is overwritten
           write(6,*) 'creating file: ',trim(fileName_given)
           period_type = 'unknown'
           call CreatenetCDFfile(trim(fileName_given),GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAX)
           ncFileID=closedID
        else !test if the defined dimensions are compatible
           !         write(6,*) 'exists: ',trim(fileName_given)
           if(projection=='lon lat') then
              call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
              call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
           else
              call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
              call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
           endif
           call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
           call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_old))
           call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_old))
           call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_old))

           !         write(6,*)'existing file ', trim(fileName_given),' has dimensions'
           !         write(6,*)GIMAX_old,GJMAX_old,KMAX_old
           if(GIMAX_old<GIMAXcdf .or. GJMAX_old<GJMAXcdf .or.  KMAX_old<KMAX)then
              write(6,*)'existing file ', trim(fileName_given),' has wrong dimensions'
              write(6,*)GIMAX_old,GIMAXcdf,GJMAX_old,GJMAXcdf,KMAX_old,KMAX
              write(6,*)'WARNING! OLD ', trim(fileName_given),' IS DELETED'
              write(6,*) 'creating new file: ',trim(fileName_given)
              period_type = 'unknown'
              call CreatenetCDFfile(trim(fileName_given),GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAX)
              ncFileID=closedID
           endif
        endif
     endif
     iotyp_new=1
     ncFileID_new=ncFileID
  endif

  if(DEBUG_NETCDF) then
    if(iotyp_new==1)then
      write(*,*)' Out_NetCDF: cases new ', trim(fileName_given), iotyp
    else
      write(*,*)' Out_NetCDF: cases old ', trim(fileName), iotyp
    end if
  endif

  if(iotyp_new==1)then
     fileName=trim(fileName_given)
     ncFileID = ncFileID_new
  elseif(iotyp==IOU_YEAR)then
     fileName = fileName_year
     ncFileID = ncFileID_year
  elseif(iotyp==IOU_MON)then
     fileName = fileName_month
     ncFileID = ncFileID_month
  elseif(iotyp==IOU_DAY)then
     fileName = fileName_day
     ncFileID = ncFileID_day
  elseif(iotyp==IOU_HOUR)then
     fileName = fileName_hour
     ncFileID = ncFileID_hour
  elseif(iotyp==IOU_INST)then
     fileName = fileName_inst
     ncFileID = ncFileID_inst
  else
     return
  endif
  if(DEBUG_NETCDF) write(*,*)'Out_NetCDF, filename ', trim(fileName), iotyp

  call CheckStop(ndim /= 2 .and. ndim /= 3, "NetCDF_ml: ndim must be 2 or 3")

  OUTtype=Real4  !default value
  if(present(CDFtype))OUTtype=CDFtype

  !buffer the wanted part of data
  ijk=0
  do k=1,kmax
     do j = 1,tljmax(me)
        do i = 1,tlimax(me)
           ijk=ijk+1
           buff(ijk)=dat(i,j,k)*scale
        enddo
     enddo
  enddo

  !send all data to me=0
  outCDFtag=outCDFtag+1

  if(me.eq.0)then

     !allocate a large array (only on one processor)
     if(OUTtype==Int1 .or. OUTtype==Int2 .or. OUTtype==Int4)then
        allocate(Idata3D(GIMAX,GJMAX,kmax), stat=alloc_err)
        call CheckStop(alloc_err, "alloc failed in NetCDF_ml")
     elseif(OUTtype==Real4)then
        allocate(R4data3D(GIMAX,GJMAX,kmax), stat=alloc_err)
        call CheckStop(alloc_err, "alloc failed in NetCDF_ml")
     elseif(OUTtype==Real8)then
        allocate(R8data3D(GIMAX,GJMAX,kmax), stat=alloc_err)
        call CheckStop(alloc_err, "alloc failed in NetCDF_ml")
     else
        WRITE(*,*)'WARNING NetCDF:Data type not supported'
     endif

     !write own data in global array
     if(OUTtype==Int1 .or. OUTtype==Int2 .or. OUTtype==Int4)then
        ijk=0
        do k=1,kmax
           do j = tgj0(me),tgj0(me)+tljmax(me)-1
              do i = tgi0(me),tgi0(me)+tlimax(me)-1
                 ijk=ijk+1
                 Idata3D(i,j,k)=buff(ijk)
              enddo
           enddo
        enddo
     elseif(OUTtype==Real4)then
        ijk=0
        do k=1,kmax
           do j = tgj0(me),tgj0(me)+tljmax(me)-1
              do i = tgi0(me),tgi0(me)+tlimax(me)-1
                 ijk=ijk+1
                 R4data3D(i,j,k)=buff(ijk)
              enddo
           enddo
        enddo
     else
        ijk=0
        do k=1,kmax
           do j = tgj0(me),tgj0(me)+tljmax(me)-1
              do i = tgi0(me),tgi0(me)+tlimax(me)-1
                 ijk=ijk+1
                 R8data3D(i,j,k)=buff(ijk)
              enddo
           enddo
        enddo
     endif

     do d = 1, NPROC-1
        CALL MPI_RECV(buff, 8*tlimax(d)*tljmax(d)*kmax, MPI_BYTE, d, &
             outCDFtag, MPI_COMM_WORLD, MPISTATUS, INFO)

        !copy data to global buffer
        if(OUTtype==Int1 .or. OUTtype==Int2 .or. OUTtype==Int4)then
           ijk=0
           do k=1,kmax
              do j = tgj0(d),tgj0(d)+tljmax(d)-1
                 do i = tgi0(d),tgi0(d)+tlimax(d)-1
                    ijk=ijk+1
                    Idata3D(i,j,k)=buff(ijk)
                 enddo
              enddo
           enddo
        elseif(OUTtype==Real4)then
           ijk=0
           do k=1,kmax
              do j = tgj0(d),tgj0(d)+tljmax(d)-1
                 do i = tgi0(d),tgi0(d)+tlimax(d)-1
                    ijk=ijk+1
                    R4data3D(i,j,k)=buff(ijk)
                 enddo
              enddo
           enddo
        else
           ijk=0
           do k=1,kmax
              do j = tgj0(d),tgj0(d)+tljmax(d)-1
                 do i = tgi0(d),tgi0(d)+tlimax(d)-1
                    ijk=ijk+1
                    R8data3D(i,j,k)=buff(ijk)
                 enddo
              enddo
           enddo
        endif
     enddo
  else
     CALL MPI_SEND( buff, 8*tlimax(me)*tljmax(me)*kmax, MPI_BYTE, 0, &
          outCDFtag, MPI_COMM_WORLD, INFO)
  endif
  !return

  if(me==0)then

     ndate(1)  = current_date%year
     ndate(2)  = current_date%month
     ndate(3)  = current_date%day
     ndate(4)  = current_date%hour

     !test if the file is already open
     if(ncFileID==closedID)then
        !open an existing netcdf dataset
        call check(nf90_open(path = trim(fileName), mode = nf90_write, &
              ncid = ncFileID), "nf90_open"//trim(fileName) )
       if(iotyp_new==1)then      !needed in case iotyp is defined
           ncFileID_new = ncFileID!not really needed
        elseif(iotyp==IOU_YEAR)then
           ncFileID_year = ncFileID
        elseif(iotyp==IOU_MON)then
           ncFileID_month = ncFileID
        elseif(iotyp==IOU_DAY)then
           ncFileID_day = ncFileID
        elseif(iotyp==IOU_HOUR)then
           ncFileID_hour = ncFileID
        elseif(iotyp==IOU_INST)then
           ncFileID_inst = ncFileID
        endif
     endif

     !test first if the variable is already defined:
     status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)

     if(status == nf90_noerr) then
!             print *, 'variable exists: ',varname
        if (DEBUG_NETCDF) write(6,*) 'Out_NetCDF: variable exists: ',varname
     else
        if (DEBUG_NETCDF) write(6,*) 'Out_NetCDF: creating variable: ',varname
        call  createnewvariable(ncFileID,varname,ndim,ndate,def1,OUTtype)
     endif


     !get variable id
     call check(nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID))


     !find the number of records already written
     call check(nf90_get_att(ncFileID, VarID, "numberofrecords",   nrecords))
     if(DEBUG_NETCDF)  print *,'number of dataset saved: ',nrecords
     !test if new record is needed
     if(present(ik).and.nrecords>0)then
        !The new record may already exist
        call idate2nctime(ndate,rdays,iotyp)
        call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeVarID))
        call check(nf90_get_var(ncFileID, timeVarID, rdays_time,start=(/ nrecords /)))
        !check if this is a newer time
        if((abs(rdays-rdays_time(1))>0.00001))then!0.00001is about 1 second
           nrecords=nrecords+1 !start a new record
        endif
     else
        !increase nrecords, to define position of new data
        nrecords=nrecords+1
     endif
     if(DEBUG_NETCDF)   print *,'writing on dataset: ',nrecords

     !append new values
     if(OUTtype==Int1 .or. OUTtype==Int2 .or. OUTtype==Int4)then
        !type Integer
        if(ndim==3)then
           if(present(ik))then
              !     print *, 'write: ',i1,i2, j1,j2,ik
              call check(nf90_put_var(ncFileID, VarID, &
                   Idata3D(i1:i2, j1:j2, 1), start = (/ 1, 1, ik,nrecords /)) )
           else
              do k=1,kmax
                 call check(nf90_put_var(ncFileID, VarID,&
                      Idata3D(i1:i2, j1:j2, k), start = (/ 1, 1, k,nrecords /)) )
              enddo
           endif
        else
           call check(nf90_put_var(ncFileID, VarID,&
                Idata3D(i1:i2, j1:j2, 1), start = (/ 1, 1, nrecords /)) )
        endif

        deallocate(Idata3D, stat=alloc_err)
        call CheckStop(alloc_err, "dealloc failed in NetCDF_ml")

     elseif(OUTtype==Real4)then
        !type Real4
        if(ndim==3)then
           if(present(ik))then
              !     print *, 'write: ',i1,i2, j1,j2,ik
              call check(nf90_put_var(ncFileID, VarID, &
                   R4data3D(i1:i2, j1:j2, 1), start = (/ 1, 1, ik,nrecords /)) )
           else
              do k=1,kmax
                 call check(nf90_put_var(ncFileID, VarID,&
                      R4data3D(i1:i2, j1:j2, k), start = (/ 1, 1, k,nrecords /)) )
              enddo
           endif
        else
           call check(nf90_put_var(ncFileID, VarID,&
                R4data3D(i1:i2, j1:j2, 1), start = (/ 1, 1, nrecords /)) )
        endif

        deallocate(R4data3D, stat=alloc_err)
        call CheckStop(alloc_err, "dealloc failed in NetCDF_ml")

     else
        !type Real8
        if(ndim==3)then
           if(present(ik))then
              !     print *, 'write: ',i1,i2, j1,j2,ik
              call check(nf90_put_var(ncFileID, VarID, &
                   R8data3D(i1:i2, j1:j2, 1), start = (/ 1, 1, ik,nrecords /)) )
           else
              do k=1,kmax
                 call check(nf90_put_var(ncFileID, VarID,&
                      R8data3D(i1:i2, j1:j2, k), start = (/ 1, 1, k,nrecords /)) )
              enddo
           endif
        else
           call check(nf90_put_var(ncFileID, VarID,&
                R8data3D(i1:i2, j1:j2, 1), start = (/ 1, 1, nrecords /)) )
        endif

        deallocate(R8data3D, stat=alloc_err)
        call CheckStop(alloc_err, "dealloc failed in NetCDF_ml")

     endif !type


     call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour0  ))
     call check(nf90_get_att(ncFileID, nf90_global, "created_hour", created_hour  ))
     call Date_And_Time(date=lastmodified_date,time=lastmodified_hour)

     !write or change attributes NB: strings must be of same length as originally
     call check(nf90_put_att(ncFileID, VarID, "numberofrecords",   nrecords))

     !update dates
     call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_date", lastmodified_date))
     call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour))
     call check(nf90_put_att(ncFileID, VarID, "current_date_last",ndate ))

     !get variable id
     call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = VarID))
     call idate2nctime(ndate,rdays,iotyp)
     call check(nf90_put_var(ncFileID, VarID, rdays, start = (/nrecords/) ) )

     !close file if present(fileName_given)
     if(iotyp_new==1)then
        call check(nf90_close(ncFileID))
     endif
  endif !me=0

  if(DEBUG_NETCDF) write(*,*)'Out_NetCDF: FINISHED '
end subroutine Out_netCDF

!_______________________________________________________________________


subroutine  createnewvariable(ncFileID,varname,ndim,ndate,def1,OUTtype)

  !create new netCDF variable

  implicit none

  type(Deriv),     intent(in) :: def1 ! definition of fields
  character (len = *),intent(in) ::varname
  integer ,intent(in) ::ndim,ncFileID,OUTtype
  integer, dimension(:) ,intent(in) ::  ndate

  integer :: iDimID,jDimID,kDimID,timeDimID
  integer :: varID,nrecords
  real :: scale
  integer :: OUTtypeCDF !NetCDF code for type

  if(OUTtype==Int1)then
     OUTtypeCDF=nf90_byte
  elseif(OUTtype==Int2)then
     OUTtypeCDF=nf90_short
  elseif(OUTtype==Int4)then
     OUTtypeCDF=nf90_int
  elseif(OUTtype==Real4)then
     OUTtypeCDF=nf90_float
  elseif(OUTtype==Real8)then
     OUTtypeCDF=nf90_double
  else
     call CheckStop("NetCDF_ml:undefined datatype")
  endif

     call check(nf90_redef(ncid = ncFileID))

     !get dimensions id
  if(projection=='lon lat') then
     call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
  else
     call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
  endif

  call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
  call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

  !define new variable
  if(ndim==3)then
     call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = OUTtypeCDF,&
          dimids = (/ iDimID, jDimID, kDimID , timeDimID/), varID=varID ) )
  elseif(ndim==2)then
     call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = OUTtypeCDF,&
          dimids = (/ iDimID, jDimID , timeDimID/), varID=varID ) )
  else
     print *, 'createnewvariable: unexpected ndim ',ndim
  endif
  !     FillValue=0.
  scale=1.
  !define attributes of new variable
  call check(nf90_put_att(ncFileID, varID, "long_name",  def1%name ))
  call check(nf90_put_att(ncFileID, varID, "coordinates", "lat lon"))
  if(trim(projection)=='Stereographic')then
     call check(nf90_put_att(ncFileID, varID, "grid_mapping", "Polar_Stereographic"))
  elseif(projection=='lon lat') then

  elseif(projection=='Rotated_Spherical')then
     call check(nf90_put_att(ncFileID, varID, "grid_mapping", "Rotated_Spherical"))
  else
     call check(nf90_put_att(ncFileID, varID, "grid_mapping",Default_projection_name ))
  endif

  nrecords=0
  call check(nf90_put_att(ncFileID, varID, "numberofrecords", nrecords))

  call check(nf90_put_att(ncFileID, varID, "units",   def1%unit))
  call check(nf90_put_att(ncFileID, varID, "class",   def1%class))

  if(OUTtype==Int1)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_byte  ))
     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))
  elseif(OUTtype==Int2)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_short  ))
     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))
  elseif(OUTtype==Int4)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_int   ))
     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))
  elseif(OUTtype==Real4)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_float  ))
  elseif(OUTtype==Real8)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_double  ))
  endif

  call check(nf90_put_att(ncFileID, varID, "current_date_first",ndate ))
  call check(nf90_put_att(ncFileID, varID, "current_date_last",ndate ))

  call check(nf90_enddef(ncid = ncFileID))

end subroutine  createnewvariable
!_______________________________________________________________________

  subroutine check(status,errmsg)
    implicit none
    integer, intent ( in) :: status
    character(len=*), intent(in), optional :: errmsg

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      if( present(errmsg) ) print *, "ERRMSG: ", trim(errmsg)
      call CheckStop("NetCDF_ml : error in netcdf routine")
    end if
  end subroutine check

  subroutine CloseNetCDF
!close open files
!NB the data in a NetCDF file is not "safe" before the file is closed. 
!The files are NOT automatically properly closed after end of program,
! and data may be lost if the files are not closed explicitely.

integer :: ncFileID

outCDFtag=0 !for avoiding too large integers

if(me==0)then

    if(ncFileID_year/=closedID)then
       ncFileID = ncFileID_year
       call check(nf90_close(ncFileID))
       ncFileID_year=closedID
    endif
    if(ncFileID_month/=closedID)then
       ncFileID = ncFileID_month
       call check(nf90_close(ncFileID))
       ncFileID_month=closedID
    endif
    if(ncFileID_day/=closedID)then
       ncFileID = ncFileID_day
       call check(nf90_close(ncFileID))
       ncFileID_day=closedID
    endif
    if(ncFileID_hour/=closedID)then
       ncFileID = ncFileID_hour
       call check(nf90_close(ncFileID))
       ncFileID_hour=closedID
    endif
    if(ncFileID_inst/=closedID)then
       ncFileID = ncFileID_inst
       call check(nf90_close(ncFileID))
       ncFileID_inst=closedID
    endif
endif
end subroutine CloseNetCDF

subroutine GetCDF(varname,fileName,Rvar,varGIMAX,varGJMAX,varKMAX,nstart,nfetch,needed)
  !
  ! open and reads CDF file
  !
  ! The nf90 are functions which return 0 if no error occur.
  ! check is only a subroutine which check wether the function returns zero
  !
  !
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

  if(me==0.and.DEBUG_NETCDF)print *,'GetCDF  reading ',trim(fileName), ' nstart ', nstart
  !open an existing netcdf dataset
  fileneeded=.true.!default
  if(present(needed))then
     fileneeded=needed
  endif

  if(fileneeded)then
     call check(nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID))
  else
     status=nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID)
     if(status/= nf90_noerr)then
        write(*,*)trim(fileName),' not found (but not needed)'
        nfetch=0
        return
     endif
  endif

  !get global attributes
  !example:
!  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_hour", attribute ))
!  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_date", attribute2 ))
!  print *,'file last modified (yyyymmdd hhmmss.sss) ',attribute2,' ',attribute

  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)

  if(status == nf90_noerr) then
     if(DEBUG_NETCDF)print *, 'variable exists: ',trim(varname)
  else
     print *, 'variable does not exist: ',trim(varname),nf90_strerror(status)
     nfetch=0
     call CheckStop(fileneeded, "NetCDF_ml : variable needed but not found")
     return
  endif

  !get dimensions
  call check(nf90_Inquire_Variable(ncFileID,VarID,name,xtype,ndims,dimids,nAtts))
  dims=0
  totsize=1
  do i=1,ndims
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(i),  len=dims(i)))
     totsize=totsize*dims(i)
      !write(*,*)'size variable ',i,dims(i)
  enddo

  write(*,*)'dimensions ',(dims(i),i=1,ndims)
  if(dims(1)>varGIMAX.or.dims(2)>varGJMAX)then
     write(*,*)'buffer too small',dims(1),varGIMAX,dims(2),varGJMAX
     Call StopAll('GetCDF buffer too small')
  endif

  if(ndims>3.and.dims(3)>varKMAX)then
     write(*,*)'Warning: not reading all levels ',dims(3),varKMAX
!     Call StopAll('GetCDF not reading all levels')
  endif

  if(nstart+nfetch-1>dims(ndims))then
     write(*,*)'WARNING: did not find all data'
     nfetch=dims(ndims)-nstart+1
     if(nfetch<=0)Call StopAll('GetCDF  nfetch<0')
  endif

  startvec=1
  startvec(ndims)=nstart
  totsize=totsize/dims(ndims)
  if(nfetch<dims(ndims))write(*,*)'fetching only',totsize*nfetch,'of ', totsize*dims(ndims),'elements'
  dims(ndims)=nfetch
  totsize=totsize*dims(ndims)

  if(xtype==NF90_SHORT.or.xtype==NF90_INT.or.xtype==NF90_BYTE)then
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
    enddo

     deallocate(Ivalues)
  elseif(xtype==NF90_FLOAT .or. xtype==NF90_DOUBLE)then
     call check(nf90_get_var(ncFileID, VarID, Rvar,start=startvec,count=dims))
  else
     write(*,*)'datatype not yet supported'!Char
     Call StopAll('GetCDF  datatype not yet supported')
  endif

  call check(nf90_close(ncFileID))

end subroutine GetCDF

subroutine WriteCDF(varname,vardate,filename_given,newfile)

!Routine to write data directly into a new file.
!used for testing only

 character (len=*),intent(in)::varname!variable name, or group of variable name
 type(date), intent(in)::vardate!variable name, or group of variable name
 character (len=*),optional, intent(in):: fileName_given!filename to which the data must be written
 logical,optional, intent(in) :: newfile

 real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: dat ! Data arrays
 character (len=100):: fileName
 real ::scale
 integer :: n,iotyp,ndim,kmax,icmp,ndate(4),nseconds
 type(Deriv) :: def1 ! definition of fields

 ndate(1)=vardate%year
 ndate(2)=vardate%month
 ndate(3)=vardate%day
 ndate(4)=vardate%hour
 call idate2nctime(ndate,nseconds)
 nseconds=nseconds+vardate%seconds
 write(*,*)nseconds

 iotyp=IOU_INST

 if(present(filename_given))then
    filename=trim(fileName_given)
 else
    filename='EMEP_OUT.nc'
 endif

 if(present(newfile))then
    if(newfile)then
    !make a new file (i.e. delete possible old one)
    if ( me == 0 )then
       write(*,*)'creating',me
       call Init_new_netCDF(fileName,iotyp)
       write(*,*)'created',me
    endif
    endif
 else
    !append if the file exist
 endif

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

 if(NSPEC_SHL+ NSPEC_ADV /=  NSPEC_TOT.and. me==0)then
    write(*,*)'WARNING: NSPEC_SHL+ NSPEC_ADV /=  NSPEC_TOT'
    write(*,*) NSPEC_SHL,NSPEC_ADV, NSPEC_TOT
    write(*,*)'WRITING ONLY SHL and ADV'
    write(*,*)'Check species names'
 endif

! def1%class='Short_lived' !written
! do n=1, NSPEC_SHL
! def1%name= species(n)%name       !written
! dat=xn_shl(n,:,:,:)
! icmp=n
! call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real8,fileName_given=fileName)
! enddo

 def1%class='Advected' !written
 do n= 1, NSPEC_ADV
 def1%name= species(NSPEC_SHL+n)%name       !written
 dat=xn_adv(n,:,:,:)
 call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real4,ist=60,jst=11,ien=107,jen=58,fileName_given=fileName)
 enddo

  elseif(trim(varname)=='LIST')then
     ndim=3 !3-dimensional
     kmax=KMAX_MID


     def1%class='Advected' !written
     do n= 1, NSPEC_ADV
        def1%name= species(NSPEC_SHL+n)%name       !written
        if(trim(def1%name)=='O3'.or.trim(def1%name)=='NO2')then
           dat=xn_adv(n,:,:,:)
           icmp=NSPEC_SHL+n
           call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real4,ist=10,jst=10,ien=20,jen=20,fileName_given=fileName)
        endif
     enddo

else

    if(me==0)write(*,*)'case not implemented'
 endif


end subroutine WriteCDF

subroutine Read_Inter_CDF(fileName,varname,Rvar,varGIMAX,varGJMAX,varKMAX,nstart,nfetch,interpol,needed)
!
!reads data from netcdf file and interpolates data into full-domain model grid
!
!the data in filename must have global coverage and lat-lon projection 
!
!Typical use: Master node reads the data with this routine and distributes the data to all subdomains.
!
!interpol='zero_order' (default) or 'bilinear'
!
!Since this routine is not working with subdomains, it does not scale (in terms of memory or cpu time)
!Consider using ReadField_CDF instead

  use netcdf
implicit none
character(len = *),intent(in) ::fileName,varname
real,intent(out) :: Rvar(*)
integer,intent(in) :: varGIMAX,varGJMAX,varKMAX,nstart
integer,intent(inout) :: nfetch
character(len = *), optional,intent(in) :: interpol
logical, optional, intent(in) :: needed
integer :: ncFileID,VarID,status,xtype,ndims,dimids(NF90_MAX_VAR_DIMS),nAtts
integer :: dims(NF90_MAX_VAR_DIMS),totsize,i,j,k
integer :: startvec(NF90_MAX_VAR_DIMS)
integer ::alloc_err
character*100 ::name
real :: scale,offset,scalefactors(2),di,dj,dloni,dlati
integer ::ig1jg1k,igjg1k,ig1jgk,igjgk,jg1,ig1,ig,jg,ijk,i361

integer, allocatable:: Ivalues(:)  ! counts valid values
real, allocatable:: Rvalues(:),Rlon(:),Rlat(:)
logical ::fileneeded
character(len = 20) :: interpol_used

fileneeded=.true.!default
if(present(needed))then
   fileneeded=needed
endif

!1)Read data
  !open an existing netcdf dataset
  status=nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID)
  if(status == nf90_noerr) then
     print *, 'reading ',trim(varname), ' from ',trim(filename)
  else
     nfetch=0
     if(fileneeded)then
     print *, 'file does not exist: ',trim(varname),nf90_strerror(status)
     call CheckStop(fileneeded, "Read_Inter_CDF : file needed but not found")
     else
     print *, 'file does not exist (but not needed): ',trim(varname),nf90_strerror(status)
        print *, 'file not needed '
     return
     endif
  endif


  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = trim(varname), varID = VarID)
  if(status == nf90_noerr) then
     if(DEBUG_NETCDF)print *, 'variable exists: ',trim(varname)
  else
     nfetch=0
     if(fileneeded)then
        print *, 'variable does not exist: ',trim(varname),nf90_strerror(status)
        call CheckStop(fileneeded, "Read_Inter_CDF : variable needed but not found")
     else
        print *, 'variable does not exist (but not needed): ',trim(varname),nf90_strerror(status)
        return
     endif
  endif

  !get dimensions id
  call check(nf90_Inquire_Variable(ncFileID,VarID,name,xtype,ndims,dimids,nAtts))
  !get dimensions
  totsize=1
  startvec=1
  dims=0
  do i=1,ndims
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(i),  len=dims(i)))
     !write(*,*)'size variable ',i,dims(i)
  enddo
  startvec(ndims)=nstart
  dims(ndims)=nfetch
  do i=1,ndims
     totsize=totsize*dims(i)
  enddo
!  write(*,*)'total size variable ',totsize
  allocate(Rvalues(totsize), stat=alloc_err)

  if(xtype==NF90_SHORT.or.xtype==NF90_INT)then
     allocate(Ivalues(totsize), stat=alloc_err)
     call check(nf90_get_var(ncFileID, VarID, Ivalues,start=startvec,count=dims))

    scalefactors(1) = 1.0 !default
    scalefactors(2) = 0.  !default
    status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
    if(status == nf90_noerr) scalefactors(1) = scale
    status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
    if(status == nf90_noerr) scalefactors(2) = offset

    do i=1,totsize
       Rvalues(i)=Ivalues(i)*scalefactors(1)+scalefactors(2)
    enddo

     deallocate(Ivalues)
  elseif(xtype==NF90_FLOAT .or. xtype==NF90_DOUBLE)then
     call check(nf90_get_var(ncFileID, VarID, Rvalues,start=startvec,count=dims))
  else
     write(*,*)'datatype not yet supported, contact Peter'!Byte or Char
     call StopAll('datatype not yet supported, contact Peter')
  endif

!2) Interpolate to proper grid
!we assume first that data is originally in lon lat grid
!check that there are dimensions called lon and lat

  call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(1), name=name ))
  if(trim(name)/='lon')goto 444
  call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(2), name=name ))
  if(trim(name)/='lat')goto 444

  allocate(Rlon(dims(1)), stat=alloc_err)
  allocate(Rlat(dims(2)), stat=alloc_err)
  status=nf90_inq_varid(ncid = ncFileID, name = 'lon', varID = VarID)
  if(status /= nf90_noerr) then
     status=nf90_inq_varid(ncid = ncFileID, name = 'LON', varID = VarID)
     if(status /= nf90_noerr) then
        write(*,*)'did not find longitude variable'
        call StopAll('did not find longitude variable')
     endif
  endif
  call check(nf90_get_var(ncFileID, VarID, Rlon))
!normalize such that  -180 < Rlon < 180
  do i=1,dims(1)
     if(Rlon(i)<-180.0)Rlon(i)=Rlon(i)+360.0
     if(Rlon(i)>180.0)Rlon(i)=Rlon(i)-360.0
  enddo
  status=nf90_inq_varid(ncid = ncFileID, name = 'lat', varID = VarID)
  if(status /= nf90_noerr) then
     status=nf90_inq_varid(ncid = ncFileID, name = 'LAT', varID = VarID)
     if(status /= nf90_noerr) then
        write(*,*)'did not find latitude variable'
        call StopAll('did not find latitude variable')
     endif
  endif
  call check(nf90_get_var(ncFileID, VarID, Rlat))

!NB: we assume regular grid
!inverse of resolution
  dloni=1.0/abs(Rlon(2)-Rlon(1))
  dlati=1.0/abs(Rlat(2)-Rlat(1))

i361=dims(1)
if(nint(Rlon(dims(1))-Rlon(1)+1.0/dloni)==360)then
   i361=1!cyclic grid
else
   write(*,*)'Read_Inter_CDF: only cyclic grid implemented'
   call StopAll('STOP')
endif

interpol_used='zero_order'!default
if(present(interpol))interpol_used=interpol

if(interpol_used=='zero_order')then
!interpolation 1:
!nearest gridcell
  ijk=0
  do k=1,varKMAX
     do j=1,varGJMAX
        do i=1,varGIMAX
           ijk=ijk+1
           ig=mod(nint((glon_fdom(i,j)-Rlon(1))*dloni),dims(1))+1!NB lon  -90 = +270
           jg=max(1,min(dims(2),nint((glat_fdom(i,j)-Rlat(1))*dlati)+1))
           igjgk=ig+(jg-1)*dims(1)+(k-1)*dims(1)*dims(2)
           Rvar(ijk)=Rvalues(igjgk)
        enddo
     enddo
  enddo
elseif(interpol_used=='bilinear')then
if(DEBUG_NETCDF)write(*,*)'bilinear interpolation',dims(1)
!interpolation 2:
!bilinear
ijk=0

  do k=1,varKMAX
     do j=1,varGJMAX
        do i=1,varGIMAX
           ijk=ijk+1
           ig=mod(floor(abs(glon_fdom(i,j)-Rlon(1))*dloni),dims(1))+1!NB lon  -90 = +270
           jg=max(1,min(dims(2),floor((glat_fdom(i,j)-Rlat(1))*dlati)+1))
           ig1=ig+1
           jg1=min(jg+1,dims(2))

           if(ig1>dims(1))ig1=i361
           jg1=jg+1

!           if(gb_glob(i,j)<Rlat(jg).or.gb_glob(i,j)>Rlat(jg1))then
!              if(ig1>1)then
!                 write(*,*)'error',gb_glob(i,j),Rlat(ig),Rlat(jg1),i,j,jg1
!                 stop
!              endif
!           endif
           di=(glon_fdom(i,j)-Rlon(ig))*dloni
           dj=(glat_fdom(i,j)-Rlat(jg))*dlati
           igjgk=ig+(jg-1)*dims(1)+(k-1)*dims(1)*dims(2)
           ig1jgk=ig1+(jg-1)*dims(1)+(k-1)*dims(1)*dims(2)
           igjg1k=ig+(jg1-1)*dims(1)+(k-1)*dims(1)*dims(2)
           ig1jg1k=ig1+(jg1-1)*dims(1)+(k-1)*dims(1)*dims(2)
           Rvar(ijk)=Rvalues(igjgk)*(1-di)*(1-dj)+Rvalues(igjg1k)*(1-di)*dj+&
                     Rvalues(ig1jgk)*di*(1-dj)+Rvalues(ig1jg1k)*di*dj
        enddo
     enddo
  enddo
else
   write(*,*)'interpolation method not recognized'
   call StopAll('interpolation method not recognized')
endif


deallocate(Rlon)
deallocate(Rlat)
deallocate(Rvalues)

  return
444 continue
  write(*,*)'NOT a longitude-latitude grid!',trim(fileName)
  write(*,*)'case not yet implemented'
  call StopAll('GetCDF case not yet implemented')


end subroutine Read_Inter_CDF


recursive subroutine ReadField_CDF(fileName,varname,Rvar,nstart,kstart,kend,interpol, &
     known_projection,  &! can be provided by user, eg. for MEGAN.
     needed,debug_flag,UnDef)
  !
  !reads data from netcdf file and interpolates data into model local (subdomain) grid
  !under development!
  !
  !dimensions and indices:
  !Rvar is assumed to have declared dimensions MAXLIMAX,MAXLJMAX in 2D.
  !If 3D, k coordinate in Rvar assumed as first coordinate. Could consider to change this.
  !nstart (otional) is the index of last dimension in the netcdf file, generally time dimension.
  !
  !undefined field values:
  !Some data can be missing/not defined for some gridpoints;
  !If Undef is present, it is used as value for these undefined gridpoints;
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
  !very CPU efficient in th epresent version.
  !Vertical interpolation is not implemented, except from "Fligh Levels", but
  !"Flight Levels" are so specific that we will probably move them in an own routine

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

  !Technical, future developements:
  !This routine is likely to change a lot in the future: can be divided into simpler routines; 
  !more functionalities will be introduced.
  !Should also be usable as standalone.
  !All MPI processes read the same file simultaneously (i.e. in parallel). 
  !They read only the chunk of data they need for themselves.

  use netcdf

  implicit none
  character(len = *),intent(in) ::fileName,varname
  real,intent(out) :: Rvar(*)
  integer,intent(in) :: nstart
  character(len = *), optional,intent(in) :: interpol
  character(len = *), optional,intent(in) :: known_projection
  logical, optional, intent(in) :: needed
  integer, optional,intent(in) :: kstart!smallest k (vertical level) to read. Default: assume 2D field
  integer, optional,intent(in) :: kend!largest k to read. Default: assume 2D field
  logical, optional, intent(in) :: debug_flag
  real, optional, intent(in) :: UnDef ! Value put into the undefined gridcells
  logical, save :: debug_ij

  integer :: ncFileID,VarID,lonVarID,latVarID,status,xtype,ndims,dimids(NF90_MAX_VAR_DIMS),nAtts
  integer :: dims(NF90_MAX_VAR_DIMS),totsize,i,j,k
  integer :: startvec(NF90_MAX_VAR_DIMS)
  integer ::alloc_err
  character*100 ::name
  real :: scale,offset,scalefactors(2),dloni,dlati
  integer ::ij,jdiv,idiv,Ndiv,Ndiv2,igjgk,ig,jg,ijk,n
  integer ::imin,imax,jmin,jjmin,jmax,igjg,k2
  integer, allocatable:: Ivalues(:)  ! I counts all data
  integer, allocatable:: Nvalues(:)  !ds counts all values
  real, allocatable:: Rvalues(:),Rlon(:),Rlat(:)
  real ::lat,lon,maxlon,minlon,maxlat,minlat
  logical ::fileneeded, debug,data3D
  character(len = 50) :: interpol_used, data_projection=""
  real :: ir,jr,Grid_resolution
  type(Deriv) :: def1 ! definition of fields
  integer, parameter ::NFL=23,NFLmax=50 !number of flight level (could be read from file)
  real :: P_FL(0:NFLmax),P_FL0,Psurf_ref(MAXLIMAX, MAXLJMAX),P_EMEP,dp!
  logical ::  OnlyDefinedValues

  real, allocatable :: Weight1(:,:),Weight2(:,:),Weight3(:,:),Weight4(:,:)
  integer, allocatable :: IIij(:,:,:),JJij(:,:,:)
  real :: FillValue=0,Pcounted
  logical :: Flight_Levels
  integer :: k_FL,k_FL2
  real, dimension(4) :: Weight
  real               :: sumWeights
  integer, dimension(4) :: ijkn
  integer :: ii, jj

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
     if ( debug ) write(*,*) 'ReadCDF reading ',trim(filename), 'nstart ', nstart
  else
     if(fileneeded)then
        print *, 'file does not exist: ',trim(fileName),nf90_strerror(status)
        call CheckStop(fileneeded, "ReadField_CDF : file needed but not found")
     else
        print *,'file does not exist (but not needed): ',trim(fileName),nf90_strerror(status)
        return
     endif
  endif


  interpol_used='zero_order'!default
  if(present(interpol))then
     interpol_used=interpol
     if ( debug ) write(*,*) 'ReadCDF interp request: ',trim(filename),':', trim(interpol)
  endif
  call CheckStop(interpol_used/='zero_order'.and.&
                 interpol_used/='conservative'.and.&
                 interpol_used/='mass_conservative',&
         'interpolation method not recognized')
  if ( debug ) write(*,*) 'ReadCDF interp set: ',trim(filename),':', trim(interpol)


  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = trim(varname), varID = VarID)
  if(status == nf90_noerr) then
     if ( debug ) write(*,*) 'ReadCDF variable exists: ',trim(varname)
  else
     !     nfetch=0
     if(fileneeded)then
        print *, 'variable does not exist: ',trim(varname),nf90_strerror(status)
        call CheckStop(fileneeded, "ReadField_CDF : variable needed but not found")
     else
        print *, 'variable does not exist (but not needed): ',trim(varname),nf90_strerror(status)
        return
     endif
  endif

  !get dimensions id
  call check(nf90_Inquire_Variable(ncFileID,VarID,name,&
       xtype,ndims,dimids,nAtts),"GetDimsId")

  !only characters cannot be handled
  call CheckStop(xtype==NF90_CHAR,"ReadField_CDF: Datatype not recognised")

  !Find whether Fill values are defined
  if ( debug ) write(*,*) 'PreFillCheck ',FillValue
  status=nf90_get_att(ncFileID, VarID, "_FillValue", FillValue)
  OnlyDefinedValues=.true.
  if(status == nf90_noerr)then
     OnlyDefinedValues=.false.
     if ( debug ) write(*,*)' FillValue (not counted)',FillValue
  else ! ds Jun 2010
     FillValue = -1.0e35  ! Should not arise in normal data
     status=nf90_get_att(ncFileID, VarID, "missing_value", FillValue)
     if(status == nf90_noerr)then
         OnlyDefinedValues=.false.
         if ( debug ) write(*,*)' FillValue found from missing_value (not counted)',FillValue
     end if
  endif
     if ( debug ) write(*,*) 'PostFillCheck ',FillValue

  !get dimensions
  startvec=1
  dims=0
  do i=1,ndims
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(i), &
          len=dims(i)),"GetDims")
     if ( debug ) write(*,*) 'ReadCDF size variable ',i,dims(i)
  enddo

  data3D=.false.
  if(present(kstart).or.present(kend))then
     call CheckStop((.not. present(kend).or. .not. present(kend)), &
          "ReadField_CDF : both or none kstart and kend should be present")
     data3D=.true.
  endif

  if( present(known_projection) ) then
     data_projection = trim(known_projection)
     if ( debug ) write(*,*) 'data known_projection ',trim(data_projection)
  else
    call check(nf90_get_att(ncFileID, nf90_global, "projection", data_projection ),"Proj")
    if ( debug ) write(*,*) 'data projection from file ',trim(data_projection)
  end if
  if(MasterProc)write(*,*)'Interpolating ',trim(varname),' from ',trim(filename),' to present grid'

  if(trim(data_projection)=="lon lat")then ! here we have simple 1-D lat, lon
     allocate(Rlon(dims(1)), stat=alloc_err)
     allocate(Rlat(dims(2)), stat=alloc_err)
     if ( debug ) write(*,"(a,a,i5,i5,a,i5)") 'alloc lon lat ',&
      trim(data_projection), alloc_err, dims(1), "x", dims(2)
  else
     allocate(Rlon(dims(1)*dims(2)), stat=alloc_err)
     allocate(Rlat(dims(1)*dims(2)), stat=alloc_err)
     if ( debug ) write(*,*) 'data allocElse ',trim(data_projection), alloc_err, trim(fileName)
  endif

  status=nf90_inq_varid(ncid = ncFileID, name = 'lon', varID = lonVarID)
  if(status /= nf90_noerr) then
     status=nf90_inq_varid(ncid = ncFileID, name = 'LON', varID = lonVarID)
     call CheckStop(status /= nf90_noerr,'did not find longitude variable')
  endif

  status=nf90_inq_varid(ncid = ncFileID, name = 'lat', varID = latVarID)
  if(status /= nf90_noerr) then
     status=nf90_inq_varid(ncid = ncFileID, name = 'LAT', varID = latVarID)
     call CheckStop(status /= nf90_noerr,'did not find latitude variable')
  endif
  if(trim(data_projection)=="lon lat")then
     call check(nf90_get_var(ncFileID, lonVarID, Rlon), 'Getting Rlon')
     call check(nf90_get_var(ncFileID, latVarID, Rlat), 'Getting Rlat')
  else
     call check(nf90_get_var(ncFileID, lonVarID, Rlon,start=(/1,1/),count=(/dims(1),dims(2)/)))
     call check(nf90_get_var(ncFileID, latVarID, Rlat,start=(/1,1/),count=(/dims(1),dims(2)/)))
  endif

  Flight_Levels=.false.


  !_______________________________________________________________________________
  !
  !2)        Coordinates conversion and interpolation
  !_______________________________________________________________________________


  if(trim(data_projection)=="lon lat")then

     !get coordinates
     !we assume first that data is originally in lon lat grid
     !check that there are dimensions called lon and lat

     call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(1), name=name ),name)
     call CheckStop(trim(name)/='lon',"longitude not found")
     call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(2), name=name ),name)
     call CheckStop(trim(name)/='lat',"latitude not found")

     if(data3D)then
        call check(nf90_inquire_dimension(ncid = ncFileID, dimID = dimids(3), name=name ))
        if(trim(name)=='FL')then
           !special vertical levels for Aircrafts
           !make table for conversion Flight Level -> Pressure
           !Hard coded because non-standard anyway. 610 meters layers
           do k=0,NFLmax
              P_FL(k)=1000*StandardAtmos_km_2_kPa(k*0.610)
           enddo
           P_FL0=P_FL(0)
           Flight_Levels=.true.
           call CheckStop(interpol_used/='mass_conservative',&
           "only mass_conservative interpolation implemented for Flight Levels")

           !need average surface pressure for the current month
           !montly average is needed, not instantaneous pressure
           call ReadField_CDF('SurfacePressure.nc','surface_pressure',&
                Psurf_ref,current_date%month,needed=.true.,interpol='zero_order',debug_flag=debug_flag)
        else
           call CheckStop(trim(name)/='k',"vertical coordinate k not found")
        endif
     endif


     !NB: we assume regular grid
     !inverse of resolution
     dloni=1.0/(Rlon(2)-Rlon(1))
     dlati=1.0/(Rlat(2)-Rlat(1))

     Grid_resolution = EARTH_RADIUS*360.0/dims(1)*PI/180.0

     !the method chosen depends on the relative resolutions
     if(.not.present(interpol).and.Grid_resolution/GRIDWIDTH_M>4)then
        interpol_used='zero_order'!usually good enough, and keeps gradients
    endif
    if ( debug ) write(*,*) 'interpol_used: ',interpol_used

     !Find chunk of data required (local)
     maxlon=max(maxval(gl_stagg),maxval(glon))
     minlon=min(minval(gl_stagg),minval(glon))
     maxlat=maxval(gb_stagg)
     minlat=minval(gb_stagg)

     if(debug) then
         write(*,*) "SET Grid resolution:" // trim(fileName), Grid_resolution
         write(*,"(a,6f8.2,2x,4f8.2)") 'ReadCDF LL stuff ',&
           Rlon(2),Rlon(1),dloni, Rlat(2),Rlat(1), dlati, &
           maxlon, minlon, maxlat, minlat
     end if

     if(dloni>0)then
        !floor((minlon-Rlon(1))*dloni)<=number of gridcells between minlon and Rlon(1)
        !mod(floor((minlon-Rlon(1))*dloni)+dims(1),dims(1))+1 = get a number in [1,dims(1)]
        imin=mod( floor((minlon-Rlon(1))*dloni)+dims(1),dims(1))+1!NB lon  -90 = +270
        imax=mod(ceiling((maxlon-Rlon(1))*dloni)+dims(1),dims(1))+1!NB lon  -90 = +270
        if(imax==1)imax=dims(1)!covered entire circle
        if(minlon-Rlon(1)<0.0)then
           imin=1!cover entire circle
           imax=dims(1)!cover entire circle
        endif
     else
        call CheckStop("Not tested: negativ dloni")
        imin=mod(floor((maxlon-Rlon(1))*dloni)+dims(1),dims(1))+1!NB lon  -90 = +270
        imax=mod(ceiling((minlon-Rlon(1))*dloni)+dims(1),dims(1))+1!NB lon  -90 = +270
        if(imax==1)imax=dims(1)!covered entire circle
     endif

     if(dlati>0)then
        jmin=max(1,min(dims(2),floor((minlat-Rlat(1))*dlati)))
        jmax=max(1,min(dims(2),ceiling((maxlat-Rlat(1))*dlati)+1))
     else!if starting to count from north pole
        jmin=max(1,min(dims(2),floor((maxlat-Rlat(1))*dlati)))!maxlat is closest to Rlat(1)
        jmax=max(1,min(dims(2),ceiling((minlat-Rlat(1))*dlati)+1))        
     endif


     if(maxlat>85.0.or.minlat<-85.0)then
        !close to poles
        imin=1
        imax=dims(1)
     endif

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
     endif

     if ( debug ) write(*,"(a,4f8.2,6i8)") 'ReadCDF minmax ',&
          minlon,maxlon,minlat,maxlat,imin,imax,jmin,jmax


     startvec(1)=imin
     startvec(2)=jmin
     if(ndims>2)startvec(ndims)=nstart
     dims=1
     dims(1)=imax-imin+1
     dims(2)=jmax-jmin+1

     if(data3D)then
        startvec(3)=kstart
        dims(3)=kend-kstart+1
        if(Flight_Levels)dims(3)=NFL
        if(Flight_Levels)startvec(3)=1
        if(ndims>3)startvec(ndims)=nstart
     endif

     totsize=1
     do i=1,ndims
        totsize=totsize*dims(i)
     enddo

     allocate(Rvalues(totsize), stat=alloc_err)
     if ( debug ) then
        write(*,"(a,1i6,a)") 'ReadCDF VarID ', VarID,trim(varname)
        do i=1, ndims ! NF90_MAX_VAR_DIMS would be 1024
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
        end if
     endif

     if(interpol_used=='conservative'.or.interpol_used=='mass_conservative')then
        !conserves integral (almost, does not take into account local differences in mapping factor)
        !takes weighted average over gridcells covered by model gridcell

        !divide the coarse grid into pieces significantly smaller than the fine grid
        !Divide each global gridcell into Ndiv x Ndiv pieces
        Ndiv=5*nint(Grid_resolution/GRIDWIDTH_M)
        Ndiv=max(1,Ndiv)
        Ndiv2=Ndiv*Ndiv
        !
        if(projection/='Stereographic'.and.projection/='lon lat')then
           !the method should be revised or used only occasionally
           if(me==0)write(*,*)'WARNING: interpolation method may be CPU demanding'
        endif
        k2=1
        if(data3D)k2=kend-kstart+1
        allocate(Ivalues(MAXLIMAX*MAXLJMAX*k2))
        allocate(Nvalues(MAXLIMAX*MAXLJMAX*k2))
        do ij=1,MAXLIMAX*MAXLJMAX*k2
           Ivalues(ij)=0
           NValues(ij) = 0
           if(present(UnDef))then
              Rvar(ij)=UnDef!default value
           else
              Rvar(ij)=0.0
           endif
        enddo

        do jg=1,dims(2)
           do jdiv=1,Ndiv
              lat=Rlat(startvec(2)-1+jg)-0.5/dlati+(jdiv-0.5)/(dlati*Ndiv)
              do ig=1,dims(1)
                 igjg=ig+(jg-1)*dims(1)
                 do idiv=1,Ndiv
                    lon=Rlon(startvec(1)-1+ig)-0.5/dloni+(idiv-0.5)/(dloni*Ndiv)
                    call lb2ij(lon,lat,ir,jr)
                    i=nint(ir)-gi0-IRUNBEG+2
                    j=nint(jr)-gj0-JRUNBEG+2
                    if(i>=1.and.i<=limax.and.j>=1.and.j<=ljmax)then
                       ij=i+(j-1)*MAXLIMAX
                       k2=1
                       if(data3D)k2=kend-kstart+1
                       if(.not. Flight_Levels)then
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

                             endif
                          enddo
                       else
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
                             enddo
                             Ivalues(ijk)=Ivalues(ijk)+1
                             Nvalues(ijk)=Nvalues(ijk)+1
                             if(k_FL<=NFL)then
                                dp=Pcounted-P_EMEP
                                igjgk=igjg+(k_FL-1)*dims(1)*dims(2)
                                Rvar(ijk)=Rvar(ijk)+Rvalues(igjgk)*dp/(P_FL(k_FL2)-P_FL(k_FL))
                                Pcounted=P_EMEP
                             endif
                          enddo

                          P_FL(0)=P_FL0!may have changed above
                       endif !Flight levels

                    endif
                 enddo
              enddo
           enddo
        enddo
        k2=1
        if(data3D)k2=kend-kstart+1
        do k=1,k2
           do i=1,limax
              do j=1,ljmax
                 ij=i+(j-1)*MAXLIMAX
                 ijk=k+(ij-1)*k2

                 debug_ij = ( DEBUG_NETCDF_RF .and. debug_proc .and. &
                                 i== debug_li .and. j== debug_lj )
                 if ( debug_ij ) write(*,*) 'DEBUG  -- INValues!' , Ivalues(ijk), Nvalues(ijk)
                 if(Ivalues(ijk)<=0.)then
                    if( .not.present(UnDef))then
                       write(*,*)'ERROR. no values found!', trim(fileName), &
                    i,j,k,me,maxlon,minlon,maxlat,minlat,glon(i,j),glat(i,j), &
                    Ivalues(ijk)
                       call CheckStop("Interpolation error")
                    endif
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

                    endif
                 endif
              enddo
           enddo
        enddo

        deallocate(Ivalues)
        deallocate(Nvalues)

     elseif(interpol_used=='zero_order')then
        !interpolation 1:
        !nearest gridcell
        ijk=0
        k2=1
        if(data3D)k2=kend-kstart+1
        do k=1,k2
           do j=1,ljmax
              do i=1,limax
                 ij=i+(j-1)*MAXLIMAX
                 ijk=k+(ij-1)*k2
                 ig=nint((glon(i,j)-Rlon(startvec(1)))*dloni)+1
                 ig=max(1,min(dims(1),ig))
                 jg=max(1,min(dims(2),nint((glat(i,j)-Rlat(startvec(2)))*dlati)+1))
                 igjgk=ig+(jg-1)*dims(1)+(k-1)*dims(1)*dims(2)
                 if(OnlyDefinedValues.or.Rvalues(igjgk)/=FillValue)then
                    Rvar(ijk)=Rvalues(igjgk)
                 else
                    Rvar(ijk)=UnDef
                 endif
              enddo
           enddo
        enddo

     endif
!_________________________________________________________________________________________________________
!_________________________________________________________________________________________________________
  else ! data_projection)/="lon lat"

     if(MasterProc.and.debug)write(*,*)'interpolating from ', trim(data_projection),' to ',trim(projection)

     call CheckStop(interpol_used=='mass_conservative', "ReadField_CDF: only linear interpolation implemented")
     if(interpol_used=='zero_order'.and.MasterProc.and.debug)&
          write(*,*)'zero_order interpolation asked, but performing linear interpolation'

     call CheckStop(data3D, "ReadField_CDF : 3D not yet implemented for general projection")
     call CheckStop(present(UnDef), "Default values filling not implemented")

     allocate(Weight1(MAXLIMAX,MAXLJMAX))
     allocate(Weight2(MAXLIMAX,MAXLJMAX))
     allocate(Weight3(MAXLIMAX,MAXLJMAX))
     allocate(Weight4(MAXLIMAX,MAXLJMAX))
     allocate(IIij(MAXLIMAX,MAXLJMAX,4))
     allocate(JJij(MAXLIMAX,MAXLJMAX,4))

!Make interpolation coefficients.
!Coefficients could be saved and reused if called several times.
     if( DEBUG_NETCDF_RF .and. debug_proc .and. i==debug_li .and. j==debug_lj ) then
        write(*,"(a)") "DEBUG_RF G2G ", me, debug_proc
     end if
     call grid2grid_coeff(glon,glat,IIij,JJij,Weight1,Weight2,Weight3,Weight4,&
                  Rlon,Rlat,dims(1),dims(2), MAXLIMAX, MAXLJMAX, limax, ljmax,&
                    ( DEBUG_NETCDF_RF .and. debug_proc ), &
                   debug_li, debug_lj )

     startvec(1)=minval(IIij(1:limax,1:ljmax,1:4))
     startvec(2)=minval(JJij(1:limax,1:ljmax,1:4))
     if(ndims>2)startvec(ndims)=nstart
     dims=1
     dims(1)=maxval(IIij(1:limax,1:ljmax,1:4))-startvec(1)+1
     dims(2)=maxval(JJij(1:limax,1:ljmax,1:4))-startvec(2)+1

     totsize=1
     do i=1,ndims
        totsize=totsize*dims(i)
     enddo

     allocate(Rvalues(totsize), stat=alloc_err)
     if ( debug ) then
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
     endif

     k=1
     do i=1,limax
        do j=1,ljmax

           Weight(1) = Weight1(i,j)  
           Weight(2) = Weight2(i,j)  
           Weight(3) = Weight3(i,j)  
           Weight(4) = Weight4(i,j)  

           ijkn(1)=IIij(i,j,1)-startvec(1)+1+(JJij(i,j,1)-startvec(2))*dims(1)
           ijkn(2)=IIij(i,j,2)-startvec(1)+1+(JJij(i,j,2)-startvec(2))*dims(1)
           ijkn(3)=IIij(i,j,3)-startvec(1)+1+(JJij(i,j,3)-startvec(2))*dims(1)
           ijkn(4)=IIij(i,j,4)-startvec(1)+1+(JJij(i,j,4)-startvec(2))*dims(1)

           ijk=i+(j-1)*MAXLIMAX
           Rvar(ijk)    = 0.0
           sumWeights   = 0.0

           do k = 1, 4
              ii = IIij(i,j,k)
              jj = JJij(i,j,k)
              if ( Rvalues(ijkn(k) )  > FillValue ) then
                   Rvar(ijk) =  Rvar(ijk) + Weight(k)*Rvalues(ijkn(k))
                  sumWeights = sumWeights + Weight(k)
              end if
           enddo !k

           if ( sumWeights > 1.0e-9 ) then
                Rvar(ijk) =  Rvar(ijk)/sumWeights
           else
                Rvar(ijk) = FillValue
           end if

        enddo !j
     enddo !i

     deallocate(Weight1)
     deallocate(Weight2)
     deallocate(Weight3)
     deallocate(Weight4)
     deallocate(IIij)
     deallocate(JJij)

  endif

  deallocate(Rvalues)
  deallocate(Rlon)
  deallocate(Rlat)
  call check(nf90_close(ncFileID))


  !  CALL MPI_FINALIZE(INFO)
  !  CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)

  return

!code below only used for testing purposes

     CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
     if(debug)write(*,*)'writing results in file',trim(varname)

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
     k2=kend-kstart+1
     n=3

     call Out_netCDF(IOU_INST,def1,n,k2, &
          Rvar,1.0,CDFtype=Real4,fileName_given='ReadField3D.nc')
!output Flight levels (reverse order of indices)
!     do k=1,2
!       do j=1,ljmax
!           do i=1,limax
!              temp(i,j,k)=0.0
!           enddo
!        enddo
!     enddo
!     do k=KMAX_MID,KMAX_MID-k2+1,-1
!        write(*,*)k,k+(KMAX_MID-k2),kend,kstart
!       do j=1,ljmax
!           do i=1,limax
!              ijk=k-(KMAX_MID-k2)+(i+(j-1)*MAXLIMAX-1)*k2
!              temp(i,j,k)=Rvar(ijk)
!           enddo
!        enddo
!     enddo
!     call Out_netCDF(IOU_INST,def1,n, KMAX_MID,&
!          temp,1.0,CDFtype=Real4,fileName_given='ReadField3D.nc')

     CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    CALL MPI_FINALIZE(INFO)
     stop
  else
      if(trim(varname)=='NOX_EMISSION')then
    n=2
     k2=1
     call Out_netCDF(IOU_INST,def1,n,k2, &
          Rvar,1.0,CDFtype=Real4,fileName_given='ReadField2D.nc')
    CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    CALL MPI_FINALIZE(INFO)
     stop
     endif
   endif

  return

end subroutine ReadField_CDF

subroutine printCDF(name, array,unit)
    ! Minimal print out to cdf, for real numbe, 2-d arrays
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

    if(MasterProc) write(*,*) "TEST printCDF :"//trim(fname),  maxval(array)
    call Out_netCDF(IOU_INST,def1,2,1, array,1.0,&
           CDFtype=Real4,fileName_given=fname,overwrite=.true.)
  end subroutine printCDF


  subroutine ReadTimeCDF(filename,TimesInDays,NTime_Read)
    !Read times in file under CF convention and convert into days since 1900-01-01 00:00:00
    character(len=*) ,intent(in)::filename
    real,intent(out) ::TimesInDays(*)
    integer, intent(in)  :: NTime_Read !number of records to read

    real, allocatable ::times(:)
    integer :: i,ntimes,status
    integer :: varID,ncFileID,ndims
    integer :: xtype,dimids(NF90_MAX_VAR_DIMS),nAtts
    integer, parameter::wordarraysize=20
    character*50 ::varname,period,since,name,timeunit,wordarray(wordarraysize)

    integer :: yyyy,mo,dd,hh,mi,ss,julian,julian_1900,diff_1900,nwords,errcode


    status = nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID)
    call CheckStop(status /= nf90_noerr, "ReadTimeCDF, file not found: "//trim(fileName))

    varname='time'
    status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)
    
    if(status == nf90_noerr) then
       if(DEBUG_NETCDF)print *, 'variable exists: ',trim(varname)
    else
       print *, 'variable does not exist: ',trim(varname),'file: ',trim(fileName),nf90_strerror(status)
       call StopAll("ReadTimeCDF : time not found")
       return
    endif
    
    call check(nf90_Inquire_Variable(ncFileID,VarID,name,xtype,ndims,dimids,nAtts))
    if(ndims>1)write(*,*)'WARNING: time has more than 1 dimension!? ',ndims
    call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(1),  len=ntimes))
    call CheckStop(ntimes<NTime_Read, "to few records in "//trim(fileName))
    allocate(times(NTime_Read))
    call check(nf90_get_var(ncFileID, VarID, times,count=(/NTime_Read/)))

    call check(nf90_get_att(ncFileID, VarID, "units", timeunit  ))

!must be of the form " xxxx since yyyy-mm-dd hh:mm:ss"

!    read(timeunit,fmt="(a,a,a,a)")period,since,date,time
    call wordsplit(trim(timeunit),wordarraysize,wordarray,nwords,errcode,separator='-')
    if(DEBUG_NETCDF.and.me==0)write(*,*)(trim(wordarray(i)),i=1,8)
    period=wordarray(1)
    call CheckStop(trim(period)/='days', "Error: only time in days implemented "//trim(period))
    since=wordarray(2)
    call CheckStop(trim(since)/='since', "since error "//trim(since))

    read(wordarray(3),*)yyyy
    read(wordarray(4),*)mo
    read(wordarray(5),*)dd
    read(wordarray(6),*)hh
    mi=0
    ss=0
!    read(wordarray(7),*)mi
!    read(wordarray(8),*)ss  !did not work ??? 


    if(yyyy/=0)then
!       read(date,fmt="(I4.4,a1,I2.2,a1,I2.2)")yyyy,s1,mo,s2,dd
!       read(time,fmt="(I2.2,a1,I2.2,a1,I2.2)")hh,s1,mi,s2,ss

       if(DEBUG_NETCDF.and.me==0)write(*,fmt="(a,I4.4,a1,I2.2,a1,I2.2,a,I2.2,a1,I2.2,a1,I2.2)")&
            'nest refdate ',yyyy,'-',mo,'-',dd,' time ',hh,':',mi,':',ss
       ss=ss+60*mi+3600*hh
       julian=julian_date(yyyy,mo,dd)
       julian_1900=julian_date(1900,1,1)
       diff_1900=julian-julian_1900
       !    if(me==0)write(*,*)'julians ',diff_1900,julian,julian_1900
       if(trim(period)=='days')then
          do i=1,NTime_Read
             TimesInDays(i)=diff_1900+times(i)+ss/(24.0*3600.0)
          enddo
       elseif(trim(period)=='seconds')then
          do i=1,NTime_Read
             TimesInDays(i)=diff_1900+(times(i)+ss)/(3600.0*24.0)
          enddo
       else
          call StopAll("ReadTimeCDF : time unit not recognized")
       endif

    else

       if(DEBUG_NETCDF.and.me==0)write(*,*)'assuming days since 0-01-01 00:00 and 365days'
       !assume units = "days since 0-01-01 00:00" 
       !and calendar = "365_day"
       yyyy=int(times(1)/365)
       
       julian=julian_date(yyyy,1,1)
       julian_1900=julian_date(1900,1,1)
       diff_1900=julian-julian_1900
       
       do i=1,NTime_Read
          TimesInDays(i)=diff_1900+times(i)-yyyy*365
       enddo
       !for leap years and dates after 28th February add one day to get Julian days
       if(mod(yyyy,4)==0)then
          do i=1,NTime_Read
             if(times(i)-yyyy*365>59.999)then
                !later than midnight the 28th february (28th Feb is 59th day)
                TimesInDays(i)=TimesInDays(i)+1.0
             endif
          enddo
!if the current date in the model is 29th of february, then this date is not defined in the 
!365 days calendar. We then assume that the 60th day is 29th of february in the netcdf file
!and not the 1st of march.
!Keep this separately as this may be defined differently in different situations.
!This implementation works for the IFS-MOZART BC
          if(current_date%month==2.and.current_date%day==29)then
             write(*,*)'WARNING: assuming 29th of February for ',trim(filename)
             do i=1,NTime_Read
                if(int(times(i)-yyyy*365)==60)then
                   !move 1st march to 29th february
                   TimesInDays(i)=TimesInDays(i)-1.0
                endif
             enddo
          endif

       endif
    endif

    deallocate(times)
  end subroutine ReadTimeCDF

end module NetCDF_ml
