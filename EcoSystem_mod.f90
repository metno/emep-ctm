! <EcoSystem_mod.f90 - A component of the EMEP MSC-W Chemical Transport Model, version v5.8>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2026 met.no
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

use Config_module ,   only: MasterProc, NLANDUSEMAX, IOU_YEAR,&
                            IOU_KEY, LandCoverInputs, MAX_NUM_DDEP_ECOS
use Debug_module,     only: DEBUG ! =>DEBUG_ECOSYSTEMS
use LandDefs_mod,     only: LandType, LandDefs
use OwnDataTypes_mod, only: Deriv, print_deriv_type, TXTLEN_DERIV, TXTLEN_SHORT, typ_s1ind
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
  NDLF=2, BDLF=3, CROP=4, SEMINAT=5, FOREST=6, WATER_D=7, NONFOREST=8, LAST_ECO=8 ! try to skipW

! We also keep the parameter for FULL_LCGRID=0 here, which is used
! for e.g. Vg values. Do not confuse LC with ECO stuff!

integer, public, parameter :: FULL_LCGRID=0

! Water_D
! *** Note *** Water_D is introduced for some NEU work, with direct 
! deposition to the water surface. This is not to be used for IAM, 
! since CCE want to have deposition to the watershed, which means 
! the grid in practice.

type(Deriv),public,dimension(:),allocatable, save:: DepEcoSystem  ! N ECOSYSTEM_OUTPUTS
logical,    public,dimension(:,:),allocatable, save:: Is_EcoSystem  ! N ECOSYSTEM_OUTPUTS, NLANDUSEMAX
real,       public,dimension(:,:,:),allocatable, save:: EcoSystemFrac


! integer, public, parameter :: MAX_NUM_DDEP_ECOS = 30  ! Outputs for Grid, Conif, etc.
 integer, public, save :: nEcoSysOutputs = 0
 integer :: i
 character(len=TXTLEN_SHORT),public,dimension(MAX_NUM_DDEP_ECOS),save :: &
  ECOSYSTEM_OUTPUTS  = [character(len=TXTLEN_SHORT):: &
    "Grid","NeedleLeaf","BroadLeaf","Crops","Seminat","Forest","Water_D", &
    "nonForest", ('NOTSET', i=9, MAX_NUM_DDEP_ECOS) ]
! Depositions
!  type(typ_s1ind), public, save, dimension(MAX_NUM_DDEP_ECOS) :: &
!  DDEP_ECOS = typ_s1ind("-",'-') ! e.g. "Grid","YMD",


contains
 !<---------------------------------------------------------------------------
subroutine Init_EcoSystems()
  character(len=TXTLEN_DERIV) :: name
  character(len=TXTLEN_SHORT) :: unit
  character(len=100) :: errmsg
  integer :: iEco, iLC
  logical, parameter :: T = .true., F = .false. ! shorthands only
  character(len=*), parameter:: dtxt='IniEcoS:'

  do iEco = 9, size(ECOSYSTEM_OUTPUTS)  ! 1 to 8 are Grid etc, defined above
    ECOSYSTEM_OUTPUTS(iEco) = LandCoverInputs%PFTs(iEco-8)
    if(MasterProc) write(*,*) dtxt//'EcoOUT:', iEco, trim(ECOSYSTEM_OUTPUTS(iEco))
    if (ECOSYSTEM_OUTPUTS(iEco) == 'NOTSET') exit
    nEcoSysOutputs = iEco 
  end do
  if(MasterProc) write(*,*) dtxt//"END", nEcoSysOutputs, ECOSYSTEM_OUTPUTS(nEcoSysOutputs)

  allocate( EcoSystemFrac(nEcoSysOutputs,LIMAX,LJMAX))
  allocate( DepEcoSystem(nEcoSysOutputs) )
  allocate( Is_EcoSystem(nEcoSysOutputs,NLANDUSEMAX)  )

  if(MasterProc) write(*,*) "Defining ecosystems: ",&
    (trim(ECOSYSTEM_OUTPUTS(iEco))," ",iEco = 1, nEcoSysOutputs)

  do iEco = 1, nEcoSysOutputs
    name = "Area_"//trim(ECOSYSTEM_OUTPUTS(iEco))//"_Frac"
    unit = "Fraction"
    if(iEco==FULL_ECOGRID) then
      name = "Area_"//trim(ECOSYSTEM_OUTPUTS(iEco))//"_km2"
      unit = "km2"
    end if

    ! Deriv(name, class,    subc,  txt,           unit
    ! Deriv index, f2d, dt_scale, scale, avg? Inst Yr Mn Day
    DepEcoSystem(iEco) = Deriv(  &
      trim(name), "EcoFrac", "Area",trim(ECOSYSTEM_OUTPUTS(iEco)) , trim(unit), &
      iEco, -99, F, 1.0, F, IOU_KEY(IOU_YEAR) )

    if(DEBUG%ECOSYSTEMS .and. MasterProc) &
      call print_deriv_type( DepEcoSystem(iEco) )
  end do

!  Define which landcovers belong to which ecosystem
  Is_EcoSystem(FULL_ECOGRID,:)    =  .true.
  Is_EcoSystem(NDLF,:)   =  LandType(:)%is_NDLF
  Is_EcoSystem(BDLF,:)   =  LandType(:)%is_BDLF
  Is_EcoSystem(CROP,:)    =  LandType(:)%is_crop
  Is_EcoSystem(SEMINAT,:) =  LandType(:)%is_seminat
  Is_EcoSystem(FOREST,:)  =  LandType(:)%is_BDLF .or. LandType(:)%is_NDLF
  Is_EcoSystem(WATER_D,:) =  LandType(:)%is_water
  Is_EcoSystem(NONFOREST,:) =  .not. Is_EcoSystem(FOREST,:)

  do iEco = 1, nEcoSysOutputs-LAST_ECO
    Is_EcoSystem(LAST_ECO+iEco,:) = &
             (LandDefs(:)%code == ECOSYSTEM_OUTPUTS(LAST_ECO+iEco))
  end do

  if ( MasterProc .and. DEBUG%ECOSYSTEMS) then
    write(*,"(a,a4,1x, a12,4a8)") 'ECOSYS  iEco', 'LC', 'Forest', 'BDLF', 'NDLF', 'NON-For'
    do iEco = 1, nEcoSysOutputs
      write(*,"(a,i3,1x,a12,4L8)") 'ECOSYS', iEco, ECOSYSTEM_OUTPUTS(iEco), &
        Is_EcoSystem(FOREST,iEco), Is_EcoSystem(BDLF,iEco),  &
        Is_EcoSystem(BDLF,iEco), Is_EcoSystem(NONFOREST,iEco) ! , Is_EcoSystem(9,iEco)
    end do
  end if


  EcoSystemFrac(:,:,:) = 0.0

end subroutine Init_EcoSystems

endmodule EcoSystem_mod
