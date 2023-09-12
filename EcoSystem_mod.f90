! <EcoSystem_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.0>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2023 met.no
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
use Debug_module,     only: DEBUG ! =>DEBUG_ECOSYSTEMS
use LandDefs_mod,     only: LandType, LandDefs
use OwnDataTypes_mod, only: Deriv, print_deriv_type, TXTLEN_DERIV, TXTLEN_SHORT
use Par_mod,          only: LIMAX, LJMAX

implicit none
private
public :: Init_EcoSystems

! depositions are calculated to the following landuse classes, where
! e.g. conif may include both temperate and Medit. forests

integer, public, parameter :: FULL_ECOGRID=1
! These are "real" ecosystems, generally consisting of more than one land-cover
! We will add land-cover specific below for Tålegrenser project
integer, private, parameter :: &
  CONIF=2, DECID=3, CROP=4, SEMINAT=5, FOREST=6, WATER_D=7, NONFOREST=8, LAST_ECO=8 ! try to skipW

! We also keep the parameter for FULL_LCGRID=0 here, which is used
! for e.g. Vg values. Do not confuse LC with ECO stuff!

integer, public, parameter :: FULL_LCGRID=0

! Water_D
! *** Note *** Water_D is introduced for some NEU work, with direct 
! deposition to the water surface. This is not to be used for IAM, 
! since CCE want to have deposition to the watershed, which means 
! the grid in practice.

! Nov2022 - adding Tålegrenser LC. Qing - continue from here - same order as in Inputs_LandDefs !
integer, public, parameter :: NDEF_ECOSYSTEMS = LAST_ECO + 16 ! first 16 LC
character(len=TXTLEN_SHORT),public,dimension(NDEF_ECOSYSTEMS),parameter :: &
  DEF_ECOSYSTEMS = [character(len=TXTLEN_SHORT):: &
    "Grid","Conif","Decid","Crops","Seminat","Forest","Water_D","nonForest", &
     "CF", "DF", "NF", "BF", "TC", "MC", "RC", "SNL", "GR", "MS", "WE", &
     "TU", "DE", "W", "ICE", "U"]  

type(Deriv),public,dimension(NDEF_ECOSYSTEMS)            ,save:: DepEcoSystem
logical,    public,dimension(NDEF_ECOSYSTEMS,NLANDUSEMAX),save:: Is_EcoSystem
real,       public,dimension(:,:,:),allocatable          ,save:: EcoSystemFrac

contains
 !<---------------------------------------------------------------------------
subroutine Init_EcoSystems()
  character(len=TXTLEN_DERIV) :: name
  character(len=TXTLEN_SHORT) :: unit
  integer :: iEco, iLC
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

    if(DEBUG%ECOSYSTEMS .and. MasterProc) &
      call print_deriv_type( DepEcoSystem(iEco) )
  end do

!  Define which landcovers belong to which ecosystem
  Is_EcoSystem(FULL_ECOGRID,:)    =  .true.
  Is_EcoSystem(CONIF,:)   =  LandType(:)%is_conif
  Is_EcoSystem(DECID,:)   =  LandType(:)%is_decid
  Is_EcoSystem(CROP,:)    =  LandType(:)%is_crop
  Is_EcoSystem(SEMINAT,:) =  LandType(:)%is_seminat
  Is_EcoSystem(FOREST,:)  =  LandType(:)%is_decid .or. LandType(:)%is_conif
  Is_EcoSystem(WATER_D,:) =  LandType(:)%is_water
  Is_EcoSystem(NONFOREST,:) =  .not. Is_EcoSystem(FOREST,:)
  do iEco = 1, NDEF_ECOSYSTEMS-LAST_ECO
    Is_EcoSystem(LAST_ECO+iEco,:) = LandDefs(:)%code == DEF_ECOSYSTEMS(LAST_ECO+iEco)
    if(MasterProc) then
      do iLC = 1, 4
        write(*,"(a,2i3,2a4,L2)") 'ADD LC-ECO', iEco, iLC, LandDefs(iLC)%code,&
           DEF_ECOSYSTEMS(LAST_ECO+iEco), Is_EcoSystem(LAST_ECO+iEco,iLC)
       end do
    end if
  end do
  if ( MasterProc ) then
    do iEco = 1, NDEF_ECOSYSTEMS
      write(*,*) 'ECOSYS', iEco, Is_EcoSystem(FOREST,iEco), Is_EcoSystem(DECID,iEco), Is_EcoSystem(NONFOREST,iEco), Is_EcoSystem(9,iEco)
    end do
  end if


  EcoSystemFrac(:,:,:) = 0.0

end subroutine Init_EcoSystems

endmodule EcoSystem_mod
