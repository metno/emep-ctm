! <EcoSystem_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2015 met.no
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
module EcoSystem_ml

 use CheckStop_ml,       only : CheckStop, StopAll
 use GridValues_ml,      only : debug_proc, debug_li, debug_lj
 use LandDefs_ml,        only : LandDefs, LandType
 use ModelConstants_ml , only : MasterProc, DEBUG_ECOSYSTEMS &
                                , NLANDUSEMAX, IOU_YEAR
 use OwnDataTypes_ml,    only : Deriv, print_deriv_type &
                                ,TXTLEN_DERIV, TXTLEN_SHORT
 use Par_ml,             only : li0, lj0, li1, lj1, MAXLIMAX, MAXLJMAX
 implicit none
 private

    public :: Init_EcoSystems

  ! depositions are calculated to the following landuse classes, where
  ! e.g. conif may include both temperate and Medit. forests

   integer, public, parameter :: FULL_ECOGRID=1
   integer, private, parameter :: &
    CONIF=2, DECID=3, CROP=4, SEMINAT=5, WATER_D=6 ! try to skip

  ! We also keep the parameter for FULL_LCGRID=0 here, which is used
  ! for e.g. Vg values. Do not confuse LC with ECO stuff!
  ! 

   integer, public, parameter :: FULL_LCGRID=0

  ! Water_D
  ! *** Note *** Water_D is introduced for some NEU work, with direct 
  ! deposition to the water surface. This is not to be used for IAM, 
  ! since CCE want to have deposition to the watershed, which means 
  ! the grid in practice.

   integer, parameter, public :: NDEF_ECOSYSTEMS = 6
   character(len=8), public, dimension(NDEF_ECOSYSTEMS), parameter :: &
    DEF_ECOSYSTEMS = (/ "Grid    " &
                      , "Conif   " &
                      , "Decid   " &
                      , "Crops   " &
                      , "Seminat " &
                      , "Water_D " /)
   type(Deriv), public, dimension(NDEF_ECOSYSTEMS), save :: DepEcoSystem

   logical, public, dimension(NDEF_ECOSYSTEMS,NLANDUSEMAX), &
            save :: Is_EcoSystem

   real, public, dimension(:,:,:), &
            save,allocatable :: EcoSystemFrac

contains
 !<---------------------------------------------------------------------------
  subroutine Init_EcoSystems()

    character(len=TXTLEN_DERIV) :: name
    character(len=TXTLEN_SHORT) :: unit
    integer :: iEco
    logical, parameter :: T = .true., F = .false. ! shorthands only

    allocate(EcoSystemFrac(NDEF_ECOSYSTEMS,MAXLIMAX,MAXLJMAX))

      if( MasterProc ) &
           write(*,*) "Defining ecosystems: ",(trim(DEF_ECOSYSTEMS(iEco))," ",iEco = 1, NDEF_ECOSYSTEMS)
    do iEco = 1, NDEF_ECOSYSTEMS

 
       name = "Area_"//trim(DEF_ECOSYSTEMS(iEco))//"_Frac"
       unit = "Fraction"
       if(iEco==FULL_ECOGRID) then
          name = "Area_"//trim(DEF_ECOSYSTEMS(iEco))//"_km2"
          unit = "km2"
       end if

          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d, dt_scale, scale, avg? Inst Yr Mn Day
        DepEcoSystem(iEco) = Deriv(  &
               trim(name), "EcoFrac", "Area",trim(DEF_ECOSYSTEMS(iEco)) , trim(unit), &
                  iEco, -99, F, 1.0, F, IOU_YEAR   )

        if(DEBUG_ECOSYSTEMS .and. MasterProc) &
             call print_deriv_type( DepEcoSystem(iEco) )

    end do

 !  Define which landcovers belong to which ecosystem

        Is_EcoSystem(FULL_ECOGRID,:)    =  .true.
        Is_EcoSystem(CONIF,:)   =  LandType(:)%is_conif
        Is_EcoSystem(DECID,:)   =  LandType(:)%is_decid
        Is_EcoSystem(CROP,:)    =  LandType(:)%is_crop
        Is_EcoSystem(SEMINAT,:) =  LandType(:)%is_seminat
        Is_EcoSystem(WATER_D,:) =  LandType(:)%is_water

    EcoSystemFrac(:,:,:) = 0.0

  end subroutine Init_EcoSystems

end module EcoSystem_ml
