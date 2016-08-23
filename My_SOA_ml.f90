! <SOA_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module OrganicAerosol_ml

   !--------------------------------------------------------------------------
   ! This module is fake - for initial 2011 public-domain ozone model only, pending
   ! decision as to which SOA scheme to release as default.
   ! Contact David.Simpson@met.no for more information if interested in SOA
   ! schemes
   !--------------------------------------------------------------------------
   use ModelConstants_ml,    only : CHEMTMIN, CHEMTMAX, &
                                    K2 => KMAX_MID, K1 => KCHEMTOP
   use PhysicalConstants_ml, only : AVOG
   use Setup_1dfields_ml,    only : itemp, xn => xn_2d
   use ChemChemicals_ml,      only : species   ! for molwts
   use ChemSpecs_tot_ml,  A1 => FIRST_SEMIVOL , A2 => LAST_SEMIVOL
   implicit none

   !/-- subroutines
    public  :: Init_OrganicAerosol
    public  :: OrganicAerosol


   !/-- public

    logical, public, parameter :: ORGANIC_AEROSOLS = .false.
! From Setup_1dfields now
!   real, public, dimension(A1:A2,K1:K2), save :: Fgas  ! Fraction in gas-phase

    character(len=*), public, parameter :: SOA_MODULE_FLAG="NotUsed"

   contains
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !+ Driver routine for Secondary Organic Aerosol  module

   subroutine Init_OrganicAerosol(i,j,debug_flag)
     integer, intent(in) :: i,j
     logical, intent(in) :: debug_flag  ! for debugging purposes only

     ! empty 
     
   end subroutine Init_OrganicAerosol

   subroutine OrganicAerosol(i,j,debug_flag)
     integer, intent(in) :: i,j
     logical, intent(in) :: debug_flag  ! for debugging purposes only

     ! empty 

   end subroutine OrganicAerosol

end module OrganicAerosol_ml
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

