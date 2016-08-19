! <My_MassBudget_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________

module My_MassBudget_ml
!_____________________________________________________________________________
  use GenSpec_adv_ml      !! Can be many species
  use My_Emis_ml,     only :    QRCSO2,  QRCNO,  QRCCO,  QRCNH3, QRCPM25 &
                              , QRCPMCO, QRCNO2

  implicit none
  private

   !-----------------  "my" mass budget terms    ---------------------------!
   !  Here we define a few indices needed to relate species with IXADV_     !
   !  indices to their equivalent emission with QRC_ index.
   !                                                                        !  
   ! Plus, the array MY_MASS_PRINT to say which species to print out        !
   !------------------------------------------------------------------------!

   !-- contains subroutine:

      public :: set_mass_eqvs     ! Called from Emissions_ml

   ! Mass budget equivalency terms

    integer, public, parameter :: N_MASS_EQVS = 7  
    integer, public, save , dimension( N_MASS_EQVS ):: &
          ixadv_eqv  & !  IXADV_ no. of species
           ,qrc_eqv    !  QRC_   no. of equivalent species


   !/** species to print out in MassBudget_ml  (old myprint)
   !  Note - we can any number of species we need here - the dimensions 
   !  are obtained in MassBudget_ml with a size command.

   integer, public, parameter, dimension(16) :: MY_MASS_PRINT = &    !SeaS
     (/  IXADV_O3, IXADV_HNO3, IXADV_PAN, IXADV_NO3, IXADV_N2O5 ,IXADV_NO, &
         IXADV_NO2,  IXADV_SO2, IXADV_SO4, IXADV_NH3, IXADV_aNH4, IXADV_aNO3 &
        ,IXADV_PM25, IXADV_PMco, IXADV_SSfi, IXADV_SSco /)  !SeaS

  contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine set_mass_eqvs()
  !---------------------------------------------------------------------  
  !+  relates say IXADV_SO2 to QRCSO2

     ! Should have dimsnions N_MASS_EQVS

       ixadv_eqv(1) = IXADV_SO2
       ixadv_eqv(2) = IXADV_CO 
       ixadv_eqv(3) = IXADV_NH3
       ixadv_eqv(4) = IXADV_PM25 
       ixadv_eqv(5) = IXADV_PMco
       ixadv_eqv(6) = IXADV_NO2  
       ixadv_eqv(7) = IXADV_NO  

       qrc_eqv(1) = QRCSO2
       qrc_eqv(2) = QRCCO 
       qrc_eqv(3) = QRCNH3 
       qrc_eqv(4) = QRCPM25
       qrc_eqv(5) = QRCPMco
       qrc_eqv(6) = QRCNO2  
       qrc_eqv(7) = QRCNO  

  end subroutine set_mass_eqvs
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module My_MassBudget_ml
