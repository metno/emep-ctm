! <EcoSystem_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module EcoSystem_mod

use Config_module ,   only: MasterProc, NLANDUSEMAX, IOU_YEAR, IOU_KEY
use Debug_module,     only: DEBUG=>DEBUG_ECOSYSTEMS
use LandDefs_mod,     only: LandType
use OwnDataTypes_mod, only: Deriv, print_deriv_type, TXTLEN_DERIV, TXTLEN_SHORT
use Par_mod,          only: LIMAX, LJMAX

implicit none
private
public :: Init_EcoSystems

! depositions are calculated to the following landuse classes, where
! e.g. conif may include both temperate and Medit. forests

integer, public, parameter :: FULL_ECOGRID=1
integer, private, parameter :: &
  CONIF=2, DECID=3, CROP=4, SEMINAT=5, FOREST=6, WATER_D=7 ! try to skipW

! We also keep the parameter for FULL_LCGRID=0 here, which is used
! for e.g. Vg values. Do not confuse LC with ECO stuff!

integer, public, parameter :: FULL_LCGRID=0

! Water_D
! *** Note *** Water_D is introduced for some NEU work, with direct 
! deposition to the water surface. This is not to be used for IAM, 
! since CCE want to have deposition to the watershed, which means 
! the grid in practice.

integer, public, parameter :: NDEF_ECOSYSTEMS = 7
character(len=8),public,dimension(NDEF_ECOSYSTEMS),parameter :: &
  DEF_ECOSYSTEMS = [character(len=8):: &
    "Grid","Conif","Decid","Crops","Seminat","Forest","Water_D"]

type(Deriv),public,dimension(NDEF_ECOSYSTEMS)            ,save:: DepEcoSystem
logical,    public,dimension(NDEF_ECOSYSTEMS,NLANDUSEMAX),save:: Is_EcoSystem
real,       public,dimension(:,:,:),allocatable          ,save:: EcoSystemFrac

contains
 !<---------------------------------------------------------------------------
subroutine Init_EcoSystems()
  character(len=TXTLEN_DERIV) :: name
  character(len=TXTLEN_SHORT) :: unit
  integer :: iEco
  logical, parameter :: T = .true., F = .false. ! shorthands only

  allocate(EcoSystemFrac(NDEF_ECOSYSTEMS,LIMAX,LJMAX))
  if(MasterProc) write(*,*) "Defining ecosystems: ",&
    (trim(DEF_ECOSYSTEMS(iEco))," ",iEco = 1, NDEF_ECOSYSTEMS)

  do iEco = 1, NDEF_ECOSYSTEMS
    name = "Area_"//trim(DEF_ECOSYSTEMS(iEco))//"_Frac"
    unit = "Fraction"
    if(iEco==FULL_ECOGRID) then
      name = "Area_"//trim(DEF_ECOSYSTEMS(iEco))//"_km2"
      unit = "km2"
    end if

    ! Deriv(name, class,    subc,  txt,           unit
    ! Deriv index, f2d, dt_scale, scale, avg? Inst Yr Mn Day
    DepEcoSystem(iEco) = Deriv(  &
      trim(name), "EcoFrac", "Area",trim(DEF_ECOSYSTEMS(iEco)) , trim(unit), &
      iEco, -99, F, 1.0, F, IOU_KEY(IOU_YEAR) )

    if(DEBUG .and. MasterProc) &
      call print_deriv_type( DepEcoSystem(iEco) )
  end do

!  Define which landcovers belong to which ecosystem
  Is_EcoSystem(FULL_ECOGRID,:)    =  .true.
  Is_EcoSystem(CONIF,:)   =  LandType(:)%is_conif
  Is_EcoSystem(DECID,:)   =  LandType(:)%is_decid
  Is_EcoSystem(CROP,:)    =  LandType(:)%is_crop
  Is_EcoSystem(SEMINAT,:) =  LandType(:)%is_seminat
  Is_EcoSystem(FOREST,:) =  LandType(:)%is_decid .or. LandType(:)%is_conif
  Is_EcoSystem(WATER_D,:) =  LandType(:)%is_water

  EcoSystemFrac(:,:,:) = 0.0

end subroutine Init_EcoSystems

endmodule EcoSystem_mod
