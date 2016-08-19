! <NetCDF_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

       module NetCDF_ml
!
! Routines for netCDF output
!
! Written by Peter january 2003
!
!compile with options:
!f90 -L/home/u4/mifahik/netcdf/lib64 -I/home/u4/mifahik/netcdf/include -64 NetCDF_ml.f90 -lnetcdf
!
!view results.nc with:
!xrdb -load /home/u4/mifahik/.app-defaults/Ncview  (once only)
!/home/u4/mifahik/bin/ncview results.nc
!or
!
!/home/u4/mifahik/bin/ncdump results.nc |less
!
!for details see: 
!http://www.unidata.ucar.edu/packages/netcdf/f90/Documentation/f90-html-docs/
!
!
!To improve: When output is onto the same file, but with different positions for the 
!lower left corner, the coordinates i_EMEP j_EMEP and long lat will be wrong
!
  use My_Derived_ml, only : model
  use My_Outputs_ml, only :FREQ_HOURLY 

  use Chemfields_ml,   only : xn_shl,xn_adv
  use CheckStop_ml,    only: CheckStop
  use Derived_ml,    only : Deriv,IOU_INST,IOU_HOUR, IOU_YEAR ,IOU_MON, IOU_DAY
  use GenSpec_shl_ml , only :NSPEC_SHL
  use GenSpec_adv_ml , only :NSPEC_ADV
  use GenSpec_tot_ml , only :NSPEC_TOT
  use GenChemicals_ml, only :species
  use GridValues_ml,   only : GRIDWIDTH_M,fi,xp,yp,xp_EMEP_official&
                                  ,yp_EMEP_official,fi_EMEP,GRIDWIDTH_M_EMEP&
                                  ,GlobalPosition,gb_glob,gl_glob,ref_latitude&
                                  ,projection, sigma_mid
  use ModelConstants_ml, only : KMAX_MID, runlabel1, runlabel2 &
                                ,NPROC, IIFULLDOM,JJFULLDOM &
                                ,PT, current_date  
  use netcdf
  use Par_ml, only : me,GIMAX,GJMAX,tgi0,tgj0,tlimax,tljmax, &
                        MAXLIMAX, MAXLJMAX,IRUNBEG,JRUNBEG
  use PhysicalConstants_ml,  only : PI       
  use TimeDate_ml, only: nmdays,leapyear ,current_date, date


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
  integer      :: ncFileID_new=closedID  !don't save because should always be redefined (in case several routines are using ncFileID_new with different filename_given)
  integer,save :: ncFileID_inst=closedID  
  integer,save :: ncFileID_hour=closedID  
  integer,save :: ncFileID_day=closedID  
  integer,save :: ncFileID_month=closedID
  integer,save :: ncFileID_year=closedID
  integer,save :: outCDFtag=0
  integer, public, parameter :: Int1=1,Int2=2,Int4=3,Real4=4,Real8=5 !CDF typr for output
  character (len=18),  parameter :: Default_projection_name = 'General_Projection'
  logical, parameter :: MY_DEBUG = .false.

  public :: Out_netCDF
  public :: CloseNetCDF
  public :: Init_new_netCDF
  public :: GetCDF
  public :: WriteCDF
  public :: ReadCDF
  public :: secondssince1970

  private :: CreatenetCDFfile
  private :: createnewvariable
  private :: check

contains
!_______________________________________________________________________


subroutine Init_new_netCDF(fileName,iotyp) 

use Par_ml,           only : GIMAX,GJMAX,IRUNBEG,JRUNBEG
use ModelConstants_ml,only : KMAX_MID   
use My_Outputs_ml,    only : NHOURLY_OUT, &      ! No. outputs
                             Asc2D, hr_out      ! Required outputs
!ds 16/12/2003use My_Derived_ml,    only :IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY  
use Derived_ml,    only :IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY  

integer,  intent(in) :: iotyp
  character(len=*),  intent(in)  :: fileName 

integer :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
integer :: ih

!write(*,*)'Init_new_netCDF ',fileName,iotyp
call CloseNetCDF

if(iotyp==IOU_YEAR)then

fileName_year = trim(fileName)
period_type = 'fullrun'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)

elseif(iotyp==IOU_MON)then

fileName_month = trim(fileName)
period_type = 'monthly'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)

elseif(iotyp==IOU_DAY)then

fileName_day = trim(fileName)
period_type = 'daily'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)

elseif(iotyp==IOU_HOUR)then

fileName_hour = trim(fileName)
period_type = 'hourly'
ISMBEGcdf=GIMAX+IRUNBEG-1; JSMBEGcdf=GJMAX+JRUNBEG-1 !initialisations
GIMAXcdf=0; GJMAXcdf=0!initialisations
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
!write(*,*)'sizes CDF ',GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
call CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf)

elseif(iotyp==IOU_INST)then

fileName_inst = trim(fileName)
period_type = 'instant'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)

else
period_type = 'unknown'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,IRUNBEG,JRUNBEG,KMAX_MID)
endif

end subroutine Init_new_netCDF

subroutine CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf,RequiredProjection)
  ! Create the netCDF file

integer, intent(in) :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
character(len=*),  intent(in)  :: fileName 
character (len=*),optional, intent(in):: requiredprojection
character (len=*), parameter :: author_of_run='Unimod group' 
character(len=*), parameter :: vert_coord='sigma: k ps: PS ptop: PT'
character (len=19) :: projection_params='90.0 -32.0 0.933013' !set later on

real :: xcoord(GIMAX),ycoord(GJMAX),kcoord(KMAX_MID)

character*8 ::created_date,lastmodified_date
character*10 ::created_hour,lastmodified_hour
integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,VarID,iVarID,jVarID,kVarID,i,j,k
integer :: iEMEPVarID,jEMEPVarID,latVarID,longVarID,PTVarID
real :: izero,jzero,scale_at_projection_origin
character*80 ::UsedProjection

  ! fileName: Name of the new created file 
  ! nf90_clobber: protect existing datasets
  ! ncFileID: netcdf ID

!Check that the dimensions are > 0
if(GIMAXcdf<=0.or.GJMAXcdf<=0.or.KMAXcdf<=0)then
write(*,*)'WARNING:'
write(*,*)trim(fileName),' not created. Requested area too small (or outside domain) '
write(*,*)'sizes (IMAX,JMAX,IBEG,JBEG,KMAX) ',GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
return
endif
   
if(present(RequiredProjection))then
   UsedProjection=trim(RequiredProjection)
else
   UsedProjection=trim(projection)
endif
  
write(*,*)'create ',trim(fileName)
write(*,*)'UsedProjection ',trim(UsedProjection)
write(*,*)'with sizes (IMAX,JMAX,IBEG,JBEG,KMAX) ',GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
  call check(nf90_create(path = trim(fileName), cmode = nf90_clobber, ncid = ncFileID))

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

  call check(nf90_def_dim(ncid = ncFileID, name = "i", len = GIMAX, dimid = iDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "j", len = GJMAX, dimid = jDimID))
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
     write(6,*) 'created_date: ',created_date
     write(6,*) 'created_hour: ',created_hour

  ! Write global attributes
  call check(nf90_put_att(ncFileID, nf90_global, "Conventions", "CF-1.0" ))
!  call check(nf90_put_att(ncFileID, nf90_global, "version", version ))
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

!  call check(nf90_put_att(ncFileID, nf90_global, "vert_coord", vert_coord))
  call check(nf90_put_att(ncFileID, nf90_global, "period_type", trim(period_type)))
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

  call check(nf90_def_var(ncFileID, "time", nf90_int, dimids = timeDimID, varID = VarID) )
  if(trim(period_type) /= 'instant'.and.trim(period_type) /= 'unknown'.and.trim(period_type) /= 'fullrun')then
  call check(nf90_put_att(ncFileID, VarID, "long_name", "time at middle of period"))
  else
  call check(nf90_put_att(ncFileID, VarID, "long_name", "time at end of period"))
  endif
  call check(nf90_put_att(ncFileID, VarID, "units", "seconds since 1970-1-1 00:00:00.0 +00:00"))
 

!CF-1.0 definitions:
  if(UsedProjection=='Stereographic')then
  call check(nf90_def_var(ncid = ncFileID, name = "Polar_Stereographic", xtype = nf90_int, varID=varID ) )
  call check(nf90_put_att(ncFileID, VarID, "grid_mapping_name", "polar_stereographic"))
  call check(nf90_put_att(ncFileID, VarID, "straight_vertical_longitude_from_pole", Fi))
  call check(nf90_put_att(ncFileID, VarID, "latitude_of_projection_origin", 90.0))
  call check(nf90_put_att(ncFileID, VarID, "scale_factor_at_projection_origin", scale_at_projection_origin))
! call check(nf90_put_att(ncFileID, VarID, "false_easting", ))
  elseif(UsedProjection=='lon lat')then 

  else

  call check(nf90_def_var(ncid = ncFileID, name = Default_projection_name, xtype = nf90_int, varID=varID ) )
  call check(nf90_put_att(ncFileID, VarID, "grid_mapping_name", trim(UsedProjection)))

  endif

  ! Leave define mode
  call check(nf90_enddef(ncFileID))

  call check(nf90_open(path = trim(fileName), mode = nf90_write, ncid = ncFileID))

! Define horizontal distances

  if(UsedProjection=='Stereographic')then

  xcoord(1)=(ISMBEGcdf-xp)*GRIDWIDTH_M/1000.
  do i=2,GIMAXcdf
     xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
!     print *, i,xcoord(i)
  enddo
  call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )

  ycoord(1)=(JSMBEGcdf-yp)*GRIDWIDTH_M/1000.
  do j=2,GJMAXcdf
     ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
  enddo
  call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAXcdf)) )

! Define horizontal coordinates in the official EMEP grid
!  xp_EMEP_official=8.
!  yp_EMEP_official=110.
!  GRIDWIDTH_M_EMEP=50000.
!  fi_EMEP=-32.
  if(fi==fi_EMEP)then
! Implemented only if fi = fi_EMEP = -32 (Otherwise needs a 2-dimensional mapping)
! uses (i-xp)*GRIDWIDTH_M = (i_EMEP-xp_EMEP)*GRIDWIDTH_M_EMEP
  do i=1,GIMAXcdf
     xcoord(i)=(i+ISMBEGcdf-1-xp)*GRIDWIDTH_M/GRIDWIDTH_M_EMEP + xp_EMEP_official
!     print *, i,xcoord(i)
  enddo
  do j=1,GJMAXcdf
     ycoord(j)=(j+JSMBEGcdf-1-yp)*GRIDWIDTH_M/GRIDWIDTH_M_EMEP + yp_EMEP_official
!     print *, j,ycoord(j)
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

!Define longitude and latitude
  call GlobalPosition !because this may not yet be done if old version of meteo is used
  if(ISMBEGcdf+GIMAXcdf-1<=IIFULLDOM .and. JSMBEGcdf+GJMAXcdf-1<=JJFULLDOM)then
  call check(nf90_put_var(ncFileID, latVarID, gb_glob(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1&
       ,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)) )
  call check(nf90_put_var(ncFileID, longVarID, gl_glob(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1&
       ,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)) )
  endif


  elseif(UsedProjection=='lon lat') then
     do i=1,GIMAXcdf
        xcoord(i)= gl_glob(i+ISMBEGcdf-1,1)
     enddo
     do j=1,GJMAXcdf
        ycoord(j)= gb_glob(1,j+JSMBEGcdf-1)
     enddo
     call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )
     call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAXcdf)) )
  else
  xcoord(1)=(ISMBEGcdf-0.5)*GRIDWIDTH_M/1000.
  do i=2,GIMAXcdf
     xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
!     print *, i,xcoord(i)
  enddo
  call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )

  ycoord(1)=(JSMBEGcdf-0.5)*GRIDWIDTH_M/1000.
  do j=2,GJMAXcdf
     ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
  enddo
  call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )
  call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAXcdf)) )
!  write(*,*)'coord written'

!Define longitude and latitude

  if(ISMBEGcdf+GIMAXcdf-1<=IIFULLDOM .and. JSMBEGcdf+GJMAXcdf-1<=JJFULLDOM)then
  call check(nf90_put_var(ncFileID, latVarID, gb_glob(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1&
       ,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)) )
  call check(nf90_put_var(ncFileID, longVarID, gl_glob(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1&
       ,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)) )
  endif

  endif
!  write(*,*)'lon lat written'

!Define vertical levels
  if(KMAXcdf==KMAX_MID)then
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

  write(*,*)'file created'

end subroutine CreatenetCDFfile

!_______________________________________________________________________

subroutine Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype,ist,jst,ien,jen,ik,fileName_given)

  !The use of fileName_given is probably slower than the implicit filename used by defining iotyp.


  integer ,intent(in) :: ndim,kmax
  type(Deriv),     intent(in) :: def1 ! definition of fields
  integer,                         intent(in) :: iotyp
  real    ,intent(in) :: scale 
  !real, dimension(:,:,:,:), intent(in) :: dat ! Data arrays
  real, dimension(MAXLIMAX,MAXLJMAX,KMAX), intent(in) :: dat ! Data arrays
  integer, optional, intent(in) :: ist,jst,ien,jen,ik !start and end of saved area. Only level ik is written if defined
  integer, optional, intent(in) :: CDFtype != OUTtype. output type (Integer*1, Integer*2,Integer*4, real*8 or real*4) 
  character (len=*),optional, intent(in):: fileName_given!filename to which the data must be written
  !NB if the file fileName_given exist (also from earlier runs) it will be appended

  character*18 :: varname
  character*8 ::lastmodified_date
  character*10 ::lastmodified_hour,lastmodified_hour0,created_hour
  integer :: varID,new,nrecords,ncFileID=closedID
  integer :: nyear,nmonth,nday,nhour,ndate(4)
  integer :: info,d,alloc_err,ijk,itag,status,i,j,k,nseconds
  integer :: i1,i2,j1,j2
  !real*4 :: buff 
  real :: buff(MAXLIMAX*MAXLJMAX*KMAX_MID) 
  real*8 , allocatable,dimension(:,:,:)  :: R8data3D
  real*4 , allocatable,dimension(:,:,:)  :: R4data3D
  integer*4, allocatable,dimension(:,:,:)  :: Idata3D
  integer :: OUTtype !local version of CDFtype
  integer :: iotyp_new
  integer :: iDimID,jDimID,kDimID,timeDimID,timeVarID
  integer :: GIMAX_old,GJMAX_old,KMAX_old
  integer :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf
  integer :: is_leap, nseconds_time(1)


  i1=1;i2=GIMAX;j1=1;j2=GJMAX  !start and end of saved area
  if(present(ist))i1=max(ist,i1)
  if(present(ien))i2=min(ien,i2)
  if(present(jst))j1=max(jst,j1)
  if(present(jen))j2=min(jen,j2)

  !Check that that the area is larger than 0
  if((i2-i1)<0.or.(j2-j1)<0.or.kmax<=0)return 

  !make variable name
  write(varname,fmt='(A)')trim(def1%name)

  !to shorten the output we can save only the components explicitely named here
  !if(varname.ne.'D2_NO2'.and.varname.ne.'D2_O3' &
  !                         .and.varname.ne.'D2_PM10')return

  !do not write 3D fields (in order to shorten outputs)
  !if(ndim==3)return

  iotyp_new=0
  if(present(fileName_given))then
     !NB if the file already exist (also from earlier runs) it will be appended
     if(me==0)then
        !try to open the file
        status=nf90_open(path = trim(fileName_given), mode = nf90_write, ncid = ncFileID)
        ISMBEGcdf=IRUNBEG+i1-1
        JSMBEGcdf=JRUNBEG+j1-1
        GIMAXcdf=i2-i1+1
        GJMAXcdf=j2-j1+1
        if(status /= nf90_noerr) then !the file does not exist yet
           write(6,*) 'creating file: ',trim(fileName_given)
           period_type = 'unknown'
           call CreatenetCDFfile(trim(fileName_given),GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAX)
           ncFileID=closedID
        else !test if the defined dimensions are compatible 
           !         write(6,*) 'exists: ',trim(fileName_given)
           if(trim(projection)=='Stereographic')then
              call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
              call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
           elseif(projection=='lon lat') then
              call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
              call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
           else
              call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
              call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
           endif
!           call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
!           call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
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
        call check(nf90_open(path = trim(fileName), mode = nf90_write, ncid = ncFileID))
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
     else
        !ds print *, 'creating variable: ',varname!,nf90_strerror(status)
        if (MY_DEBUG) write(6,*) 'creating variable: ',varname!,nf90_strerror(status)
        call  createnewvariable(ncFileID,varname,ndim,ndate,def1,OUTtype)
     endif


     !get variable id
     call check(nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID))


     !find the number of records already written
     call check(nf90_get_att(ncFileID, VarID, "numberofrecords",   nrecords))
     !  print *,'number of dataset saved: ',nrecords
     !test if new record is needed
     if(present(ik).and.nrecords>0)then
        !The new record may already exist 
        !use time as record reference, (instead of "numberofrecords")
        call secondssince1970(ndate,nseconds,iotyp)
        call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeVarID))
        call check(nf90_get_var(ncFileID, timeVarID, nseconds_time,start=(/ nrecords /)))   
        !check if this is a newer time
        if((nseconds/=nseconds_time(1)))then
           nrecords=nrecords+1 !start a new record
        endif
     else
        !increase nrecords, to define position of new data
        nrecords=nrecords+1
     endif
     !  print *,'writing on dataset: ',nrecords

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
     !    print *, 'date now: ',lastmodified_hour,' date before ',lastmodified_hour0,' date start ', created_hour

     !write or change attributes NB: strings must be of same length as originally

     call check(nf90_put_att(ncFileID, VarID, "numberofrecords",   nrecords))

     !update dates
     call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_date", lastmodified_date))
     call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour))
     call check(nf90_put_att(ncFileID, VarID, "current_date_last",ndate ))

     !get variable id
     call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = VarID))
     call secondssince1970(ndate,nseconds,iotyp)!middle of period: !NB WORKS ONLY FOR COMPLETE PERIODS
     


     call check(nf90_put_var(ncFileID, VarID, nseconds, start = (/nrecords/) ) )

     !close file if present(fileName_given)
     if(iotyp_new==1)then
        call check(nf90_close(ncFileID))
     endif
  endif !me=0

  return
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
  character (len = 50) ::tmpstring

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
  if(trim(projection)=='Stereographic')then
     call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
  elseif(projection=='lon lat') then
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
        call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = OUTtypeCDF,     &
             dimids = (/ iDimID, jDimID, kDimID , timeDimID/), varID=varID ) )
     elseif(ndim==2)then
        call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = OUTtypeCDF,     &
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
!     call check(nf90_put_att(ncFileID, varID, "periodlength",   "yearly"))

!25/10/2005     call check(nf90_put_att(ncFileID, varID, "xfelt_ident",ident ))
     call check(nf90_put_att(ncFileID, varID, "current_date_first",ndate ))
     call check(nf90_put_att(ncFileID, varID, "current_date_last",ndate ))
   
     call check(nf90_enddef(ncid = ncFileID))

end subroutine  createnewvariable
!_______________________________________________________________________

  subroutine check(status)
    implicit none
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      call CheckStop("NetCDF_ml : error in netcdf routine") 
    end if
  end subroutine check  

  subroutine CloseNetCDF
!close open files
!NB the data in a NetCDF file is not "safe" before the file
!is closed. The files are NOT automatically properly 
!closed after end of program, and data may be lost if the files are not
!closed explicitely.

use Par_ml,           only : me

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

  subroutine secondssince1970(ndate,nseconds,iotyp)
    !calculate how many seconds have passed since the start of the year 1970


    integer, intent(in) :: ndate(4)
    integer, intent(out) :: nseconds
    integer, optional, intent(in):: iotyp
    integer :: n,nday,is_leap 

    nday=0
    do n=1,ndate(2)-1
       nday=nday+nmdays(n)
    enddo
    nday=nday+ndate(3)

    nseconds=3600*(ndate(4)+24*(nday-1))

!add seconds from each year since 1970
    do n=1970,ndate(1)-1
       is_leap=0
       if (leapyear(n))is_leap=1
       nseconds=nseconds+24*3600*365+24*3600*is_leap
    enddo

    if(present(iotyp))then
       !middle of period: !NB WORKS ONLY FOR COMPLETE PERIODS
       is_leap=0
       if (leapyear(ndate(1)-1))is_leap=1
       if(iotyp==IOU_YEAR)then
          !take end of run date
          nseconds=nseconds       
!          nseconds=nseconds-43200*365-43200*is_leap
       elseif(iotyp==IOU_MON)then
          nseconds=nseconds-43200*nmdays(max(ndate(2)-1,1))!nmdays(jan)=nmdays(dec)
       elseif(iotyp==IOU_DAY)then
          nseconds=nseconds-43200 !24*3600/2=43200
       elseif(iotyp==IOU_HOUR)then
          nseconds=nseconds-1800*FREQ_HOURLY  !1800=half hour
       elseif(iotyp==IOU_INST)then
          nseconds=nseconds       
       else
          nseconds=nseconds       
       endif
    endif
  end subroutine secondssince1970


subroutine GetCDF(varname,fileName,var,varGIMAX,varGJMAX,varKMAX,nstart,nfetch,needed)
  !
  ! open and reads CDF file
  !
  ! The nf90 are functions which return 0 if no error occur.
  ! check is only a subroutine which check wether the function returns zero
  !
  !

  character (len=*),intent(in) :: fileName 

  character (len = *),intent(in) ::varname
  integer, intent(in) :: nstart,varGIMAX,varGJMAX,varKMAX
  integer, intent(inout) ::  nfetch
  real, dimension(*),intent(out) :: var
  logical, optional,intent(in) :: needed
!  real, dimension(varGIMAX*varGJMAX*varKMAX*NFETCH),intent(out) :: var
!  real, dimension(132,111,Nrec),intent(out) :: var


  logical :: fileneeded
  integer :: GIMAX,GJMAX,KMAX_MID,nrecords,period
  integer :: status,ndims,alloc_err
  integer :: n,KMAX,Nrec,ijn,ijkn,timeVarID
  integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,VarID,iVarID,jVarID,kVarID,i,j,k
  integer :: var_date(9000),ndate(4)
  real , allocatable,dimension(:,:,:,:)  :: values
  real ::depsum
  character*20::attribute,attribute2

!  Nrec=size(var,3)

  print *,'  reading ',trim(fileName)
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
     print *, 'variable exists: ',trim(varname)
  else
     print *, 'variable does not exist: ',trim(varname),nf90_strerror(status)
     nfetch=0
     call CheckStop(fileneeded, "NetCDF_ml : variable needed but not found") 
     return
  endif

  !get dimensions id
  status=nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID)
  if(status /= nf90_noerr) then     
  call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
  endif
  status=nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID)
  if(status /= nf90_noerr) then     
  call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
  endif
!  call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
!  call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
  call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
  call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

  !get dimensions length
  call check(nf90_inquire_dimension(ncid=ncFileID, dimID=idimID,  len=GIMAX))
  call check(nf90_inquire_dimension(ncid=ncFileID, dimID=jdimID,  len=GJMAX))
  call check(nf90_inquire_dimension(ncid=ncFileID, dimID=kdimID,  len=KMAX_MID))
  call check(nf90_inquire_dimension(ncid=ncFileID, dimID=timedimID,  len=nrecords))
  Nrec=nrecords

  print *, 'dimensions ',GIMAX,GJMAX,KMAX_MID,nrecords
  if(GIMAX>varGIMAX.or.GJMAX>varGJMAX)then
     write(*,*)'buffer too small',GIMAX,varGIMAX,GJMAX,varGJMAX
     stop
  endif

  ndims=4
  if(KMAX_MID==1)ndims=3
  !get variable info
!  call check(nf90_inquire_variable(ncFileID, varID, ndims=ndims))
!       print *, 'dimensions ',ndims
  if(KMAX_MID>varKMAX)then
     write(*,*)'Warning: not reading all levels ',KMAX_MID,varKMAX
!     stop
  endif

  if(nstart+nfetch-1>nrecords)then
     write(*,*)'WARNING: did not find all data'
     nfetch=nrecords-nstart+1
     if(nfetch<=0)stop
  endif
!  if(nfetch>Nrec)then
!     write(*,*)'buffer too small. Increase last dimension',nfetch,Nrec
!     stop
!  endif
  if(ndims==3)then
     kmax=1
     !allocate a 2D array 
     allocate(values(GIMAX,GJMAX,nfetch,1), stat=alloc_err)
     if ( alloc_err /= 0 ) then
        print *, 'alloc failed in ReadCDF_ml: ',alloc_err,ndims
        stop
     endif
  elseif(ndims==4)then
     kmax=KMAX_MID
     !allocate a 3D array 
     allocate(values(GIMAX,GJMAX,KMAX_MID,nfetch), stat=alloc_err)
     if ( alloc_err /= 0 ) then
        print *, 'alloc failed in ReadCDF_ml: ',alloc_err,ndims
        stop
     endif

  else
     print *, 'unexpected number of dimensions: ',ndims
     stop
  endif

  !get variable attributes
  !example:
  attribute=''
!  call check(nf90_get_att(ncFileID, VarID, "long_name", attribute))
!       print *,'long_name ',attribute

!  call check(nf90_get_att(ncFileID, VarID, "xfelt_ident", xfelt_ident))

  !get time variable
  call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeVarID))
  call check(nf90_get_var(ncFileID, timeVarID, var_date,start=(/ nstart /),count=(/ nfetch /)))

  !get variable
  if(ndims==3)then
     !      call check(nf90_get_var(ncFileID, VarID, values,start=(/ 1, 1, nstart /),count=(/ 1, 1, nfetch /)))
     call check(nf90_get_var(ncFileID, VarID, values,start=(/ 1, 1, nstart /),count=(/ GIMAX,GJMAX,nfetch /)))
  elseif(ndims==4)then
     call check(nf90_get_var(ncFileID, VarID, values,start=(/1,1,1,nstart/),count=(/GIMAX,GJMAX,KMAX_MID,nfetch /)))
  endif
  if(Nfetch<nrecords)then
     write(*,*)'Reading record',nstart,' to ',nstart+nfetch-1
  endif
!  write(*,*)'date start '
!  call datefromsecondssince1970(ndate,var_date(1))
!  write(*,*)'date end '
!  call datefromsecondssince1970(ndate,var_date(nfetch))
!  period=(var_date(2)-var_date(1))/3600.

  if(ndims==3)then
     do n=1,nfetch
        do j=1,GJMAX
           do i=1,GIMAX
              ijn=i+(j-1)*varGIMAX+(n-1)*varGJMAX*varGIMAX
              var(ijn)=values(i,j,n,1)
           enddo
        enddo
        !     if(n<10)write(*,*)n,values(1,1,n,1)
     enddo
  else
     if(varKMAX==1)then
     do n=1,nfetch
        do k=KMAX_MID,KMAX_MID
        do j=1,GJMAX
           do i=1,GIMAX
              ijkn=i+(j-1)*varGIMAX+(k-1)*varGJMAX*varGIMAX+(n-1)*varGJMAX*varGIMAX*varKMAX
              var(ijkn)=values(i,j,k,n)
           enddo
        enddo
        enddo
     enddo
     else
     do n=1,nfetch
        do k=1,min(KMAX_MID,varKMAX)
        do j=1,GJMAX
           do i=1,GIMAX
              ijkn=i+(j-1)*varGIMAX+(k-1)*varGJMAX*varGIMAX+(n-1)*varGJMAX*varGIMAX*varKMAX
              var(ijkn)=values(i,j,k,n)
           enddo
        enddo
        enddo
     enddo
     endif
!     stop
  endif
  deallocate(values)
  call check(nf90_close(ncFileID))

end subroutine GetCDF

subroutine WriteCDF(varname,vardate,filename_given,newfile)


 character (len=*),intent(in)::varname!variable name, or group of variable name
 type(date), intent(in)::vardate!variable name, or group of variable name
 character (len=*),optional, intent(in):: fileName_given!filename to which the data must be written
 logical,optional, intent(in) :: newfile

 real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: dat ! Data arrays
 character (len=100):: fileName
 real ::scale
 integer :: n,iotyp,ndim,kmax,icmp,dim,ndate(4),nseconds
 type(Deriv) :: def1 ! definition of fields

 ndate(1)=vardate%year
 ndate(2)=vardate%month
 ndate(3)=vardate%day
 ndate(4)=vardate%hour
 call secondssince1970(ndate,nseconds)
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
 def1%rho=.false.      !not used
 def1%inst=.true.      !not used
 def1%year=.false.     !not used
 def1%month=.false.    !not used
 def1%day=.false.      !not used
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

subroutine ReadCDF(varname,vardate,filename_given)


 character (len=*),intent(in)::varname!variable name, or group of variable name
 type(date), intent(in)::vardate!variable name, or group of variable name
 character (len=*),optional, intent(in):: fileName_given!filename to which the data must be written

 real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: dat ! Data arrays
 real, allocatable,dimension(:,:,:) :: buf ! Data arrays
 character (len=100):: fileName
 integer :: n,ndim,kmax,ndate(4),nseconds,nstart,nfetch

 ndate(1)=vardate%year
 ndate(2)=vardate%month
 ndate(3)=vardate%day
 ndate(4)=vardate%hour
 call secondssince1970(ndate,nseconds)
 nseconds=nseconds+vardate%seconds
 write(*,*)nseconds

 filename=trim(fileName_given)

 
 if(trim(varname)=='ALL')then
 ndim=3 !3-dimensional
 kmax=KMAX_MID

 if(NSPEC_SHL+ NSPEC_ADV /=  NSPEC_TOT.and. me==0)then
    write(*,*)'WARNING: NSPEC_SHL+ NSPEC_ADV /=  NSPEC_TOT'
    write(*,*) NSPEC_SHL,NSPEC_ADV, NSPEC_TOT
    write(*,*)'WRITING ONLY SHL and ADV'
    write(*,*)'Check species names'
 endif
if(me==0)then
   allocate(buf(GIMAX,GJMAX,KMAX))
else
   allocate(buf(1,1,1)) !just to have the array defined
endif

 do n=1, NSPEC_SHL
    nstart=1
    nfetch=1
if(me==0)then
   call  GetCDF(species(n)%name,fileName,buf,GIMAX,GJMAX,KMAX_MID,nstart,nfetch)
endif
 call global2local(buf,dat,77,1,GIMAX,GJMAX,KMAX_MID,1,1)

   xn_shl(n,:,:,:)=dat
 enddo

 do n= 1, NSPEC_ADV
    nstart=1
    nfetch=1
if(me==0)then
    call  GetCDF(species(NSPEC_SHL+n)%name,fileName,buf,GIMAX,GJMAX,KMAX_MID,nstart,nfetch)
endif
 call global2local(buf,dat,78,1,GIMAX,GJMAX,KMAX_MID,1,1)
    xn_adv(n,:,:,:)=dat
 enddo

 else

    if(me==0)write(*,*)'case not implemented'
 endif
deallocate(buf)

end subroutine ReadCDF

end module NetCDF_ml
