! <SoilWater_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2010-2011 met.no
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
module SoilWater_mod
 use Config_module,      only: USES
 use Debug_module,       only: DEBUG_SOILWATER
 use GridValues_mod,     only: debug_proc, debug_li, debug_lj, i_fdom, j_fdom,&
                             longitude => glon
 use Io_Progs_mod,       only: PrintLog
 use Landuse_mod,        only: water_fraction
 use LocalVariables_mod, only: Grid
 use Met_mod,            only: extendarea
 use MetFields_mod,      only: SoilWater_deep, SoilWaterSource,fSW &
                               ,foundSoilWater_deep  ! false if no SW-deep
 use Par_mod,            only: limax, ljmax, me
 use TimeDate_mod,       only: current_date, daynumber


    implicit none

    public :: Set_SoilWater


!
!                        :-------------------
!                       /:
!                      / :
!                     /  :
!                    /   :
!                   /    :
!                  /     :
!  ----------------------:-----------------x
!                 0      D                MAM
!
!   MAM = max possible water, whatever units are available
!    (= 0.02 m in PARLAM data, or ca. 0.7 in IFS (organic soils, wetlands
!     have v. high values). 
!   DAM = D/MAM, min 0, max 1.0

    real, private, parameter :: SoilDAM = 0.5  ! assumed fraction of
                            ! MAXM when decline begins, D on figure
                            ! DO NOT SET TO ZERO!
!    real, dimension(366), public, save :: SWP = 0.0  ! daily soil water potential
                                              ! in  MPa


contains

  ! WARNING - THE SOIL MOISTURE WORK IS STILL UNDERWAY, AND IS NOT
  ! FUNCTIONING FOR ALL POSSIBLE METEOROLOGY INPUTS.
  ! If in doubt, set USE_SOILWATER = .false. in Config_module
   subroutine Set_SoilWater()
      integer :: i, j, hourloc
      logical :: my_first_call = .true.
      logical :: mydebug
      real    :: REW       !  Relative soil water

      if( DEBUG_SOILWATER .and. debug_proc ) write(*,*) "DEBUG_SW START: ", &
        current_date%day, current_date%hour, current_date%seconds

      if ( .not. USES%SOILWATER  ) return ! and fSW has been set to 1. at start
      if ( .not. foundSoilWater_deep  ) then
        if( my_first_call ) &
           call PrintLog("WARNING: USES%SOILWATER=true, but no deep SW found")
        my_first_call = .false.
        return ! and fSW has been set to 1. at start
      end if


      ! We reset once per day, but need to loop through the cells to find
      ! the 3am reset point

      if ( current_date%seconds /= 0 .and. .not. my_first_call ) return

        do j = 1, ljmax
           do i = 1, limax

             hourloc= mod(nint(current_date%hour+24*(1+longitude(i,j)/360.0)),24)

             mydebug = ( DEBUG_SOILWATER .and. debug_proc.and. i==debug_li.and.j==debug_lj ) 
             if ( mydebug ) write(*,*) "CHECK_SWF", hourloc, " date ", current_date


             if ( my_first_call ) hourloc = 3 ! fake to get started
             if ( hourloc /= 3  ) cycle  ! Only set one per day, at 3am

             REW      = SoilWater_deep(i,j,1) !!!!/ SoilMAM ! Now done in Met_mod

             if ( REW < SoilDAM ) then
               fSW(i,j) = REW/SoilDAM
             else
               fSW(i,j) = 1.0
             end if

             if ( mydebug ) then
               write(*,"(a,i4,f7.4,i4,2f12.4)") "RESET_SWF: ", &
                 daynumber, SoilWater_deep(i,j,1), hourloc, REW, fSW(i,j)
             end if

        end do
        end do

      if ( DEBUG_SOILWATER .and. debug_proc ) then
         i = debug_li
         j = debug_lj
         hourloc= mod(nint(current_date%hour+24*(1+longitude(i,j)/360.0)),24)
         REW   = SoilWater_deep(i,j,1) !done:/ SoilMAM

         write(*,"(a,f7.4,i4,f7.4,i4,2f12.4)") "DEBUG_SWF: ", &
           water_fraction(i,j), daynumber, SoilWater_deep(i,j,1), hourloc,&
            REW, fSW(i,j)
             
      end if

      my_first_call = .false.

   end subroutine Set_SoilWater
end module SoilWater_mod
