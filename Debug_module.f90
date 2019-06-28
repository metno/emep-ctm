! <Debug_module.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
!> TEMPORARY MODULE for consistency with ecosx
!  Will move all DEBUG from Config_module here one day
module Debug_module

  ! DebugCell is set in Solver_mod when DEBUG%RUNCHEM is true, i=debug_li,
  ! j=debug_lj,  and k==20. Allows extra debug for one cell

 logical, public, save ::  DebugCell  = .false.

 type, public :: emep_debug
  logical :: &
     AOT             = .false. &
    ,AEROSOL         = .false. & ! ...needed for intended debugs are to work
    ,AQUEOUS         = .false. &
    ,BCS             = .false. & ! BoundaryConditions
    ,BIO             = .false. & ! Biogenic emissions
    ,BIDIR           = .false. & ! FUTURE Bi-directional exchange
    ,COLUMN          = .false. & ! Used in Derived_mod for column integration
    ,COLSRC          = .false. & ! Volcanic emissions and Emergency scenarios
    ,DERIVED         = .false. & !
    ,DRYDEP          = .false. & ! 
    ,DRYRUN          = .false. & ! Skips fast chemistry to save some CPU
    ,EMISSIONS       = .false. & ! 
    ,EQUIB           = .false. &   !MARS, EQSAM etc.
    ,FORESTFIRE      = .false. &
    ,GETEMIS         = .false. &
    ,GLOBBC          = .false. &
    ,GRIDVALUES      = .false. &
    ,HOURLY_OUTPUTS  = .false. & !
    ,IOPROG          = .false. &
    ,LANDDEFS        = .false. &
    ,MAINCODE        = .false. & !< debugs main code (emepctm) driver
    ,MET             = .false. &
    ,MOSAICS         = .false. &
    ,MY_DERIVED      = .false. &
    ,pH              = .false. &
    ,PHYCHEM         = .false. &
    ,POLLEN          = .false. &
    ,ROADDUST        = .false. &
    ,RSUR            = .false. & ! Surface resistance
    ,RUNCHEM         = .false. & ! DEBUG%RUNCHEM is SPECIAL, need for some other debugs
       ,MY_WETDEP    = .false. &
    ,SEASALT         = .false. &
    ,SETUP_1DCHEM    = .false. &
    ,SETUP_1DBIO     = .false. &
    ,SITES           = .false. &
    ,SOILNOX         = .false. &
    ,SOLVER          = .false. &
    ,STOFLUX         = .false. &
    ,VDS             = .false.
  ! integer debug options allow different levels of verbosity
   integer               :: &
      PFT_MAPS  = 0         & !< Future option
     ,LANDUSE   = 0         & !
     ,DO3SE     = 0         & !
     ,SOA       = 0         &
     ,SUBMET    = 2         & ! Use 999 for all land-cover, otherwise LC index
     ,STOP_HH   = -1          ! If positive, code will quite when hh==STOP_HH
  !----------------------------------------------------------
   integer, dimension(2) :: IJ = [-999,-999]  ! index for debugging print out
   character(len=20)     :: SPEC = 'O3'       ! default.
   character(len=20)     :: datetxt = '-'       ! default.
   integer               :: ISPEC = -999      ! Will be set after NML
end type emep_debug
type(emep_debug), public, save :: DEBUG

! Older style, awaiting conversion
logical, public, parameter ::    &
   DEBUG_ADV            = .false. &
  ,DEBUG_BLM            = .false. & ! Produces matrix of differnt Kz and Hmix
  ,DEBUG_DERIVED        = .false. &
  ,DEBUG_ECOSYSTEMS     = .false. &
  ,DEBUG_EMISSTACKS     = .false. &
  ,DEBUG_Kz             = .false. &
  !!,DEBUG_DRYDEP         = .false. &
    ,DEBUG_MY_DRYDEP    = .false. &
    ,DEBUG_CLOVER       = .false. &
  ,DEBUG_EMISTIMEFACS   = .false. &
  ,DEBUG_LANDIFY        = .false. &
  ,DEBUG_MASS           = .false. &
  ,DEBUG_NEST           = .false. &
  ,DEBUG_NEST_ICBC      = .false. & ! IFS-MOZART/C-IFS BC
  ,DEBUG_NETCDF         = .false. &
  ,DEBUG_NETCDF_RF      = .false. & ! ReadField_CDF in NetCDF_mod
  ,DEBUG_NH3            = .false. & ! NH3Emis experimental
  ,DEBUG_OUTPUTCHEM     = .false. & ! Output of netcdf results
  ,DEBUG_OUT_HOUR       = .false. & ! Debug Output_hourly.f90
    ,DEBUG_DUST           = .false. & ! Skips fast chemistry to save some CPU
    ,DEBUG_WETDEP       = .false. &
  ,DEBUG_RB             = .false. &
  ,DEBUG_SOILWATER      = .false. 

end module Debug_module
