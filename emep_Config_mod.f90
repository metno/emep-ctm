! <emep_Config_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.15>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2017 met.no
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
module emep_Config_mod
  !---------------------------------------------------------------------
  ! A start of a more general module to store variables which
  ! can be easily changed though the namelist system. Will
  ! take over some of the job of the current ModelConstants;
  ! the latter stopped using just constants years ago ;-)
  !---------------------------------------------------------------------
  implicit none
  private

  type, private :: PBL_t
    real :: ZiMIN = 100.0                     ! minimum mixing height
    real :: ZiMAX = 3000.0                    ! maximum mixing height
    character(len=10) :: HmixMethod = "JcRb"  ! Method used for Hmix
      ! JcRb = Jericevic/Richardson number method
      ! "SbRb"= Seibert !"TIZi" = Original from Trond Iversen tiphysics
  end type PBL_t
  type(PBL_t), public, save :: PBL = PBL_t()

  type, private :: EmBio_t
    character(len=10) :: GlobBvocMethod = '-' ! can be MEGAN
    real :: IsopFac = 1.0                     ! for experiments
    real :: TerpFac = 1.0                     ! for experiments
  ! canopy light factor, 1/1.7=0.59, based on Lamb 1993 (cf MEGAN 0.57)
    real :: CLF     = 1.0                     ! canopy factor, leaf vs branch emissions
  end type EmBio_t
  type(EmBio_t), public, save :: EmBio = EmBio_t()

 ! We allow a flexible string which can switch between different
 ! experiments called by e.g. Solver. A but crude, but
 ! it makes sure the experiments are recorded in the config
 ! system

  character(len=100), save, public :: YieldModifications = 'VBS' ! Default for EmChem16mt

  
  type, private :: LandCoverInputs_t
    character(len=200), dimension(2) :: MapFile = 'NOTSET'  ! Usually PS European + global
    character(len=200) :: LandDefs = '-'   !  LAI, h, etc (was Inputs_LandDefs
    character(len=200) :: Do3seDefs = '-'  !  DO3SE inputs
  end type LandCoverInputs_t
  type(LandCoverInputs_t), public, save :: LandCoverInputs=LandCoverInputs_t()

end module emep_Config_mod
