! <SOA_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module OrganicAerosol_ml
! From BOX3D on 9/8/01
   ! Calculates the amount of condesible species in the gas and aerosol
   ! phases. 
   ! Methodology from Andersson-Sk\"old and Simpson, 2000, Secondary Organic
   ! Aerosol Formation in Northern Europe: a Model Study, to be published
   ! in JGR.
   !
   ! Usage: call OrganicAerosol from Runchem, after setup of 1d-fields
   ! finished.  The subroutine initialises itself on the first call
   ! and thereafter modifies two external variables:
   !   xn(SOA,k) : the concentrations of SOA 
   !   Fgas(X,k) : The fraction of X which is gas and not aeorosl
   !
   ! Dave Simpson, August 2001 
   !--------------------------------------------------------------------------
   ! nb - we use all of GenSpec_tot_ml since different model versions
   !  will have different species names. This module is intended to
   !  insensitive to the actual names one day, so this should be
   !  revised .. one day   - ds.
   !--------------------------------------------------------------------------
   use ModelConstants_ml,    only : CHEMTMIN, CHEMTMAX, &
                                    K2 => KMAX_MID, K1 => KCHEMTOP
   use PhysicalConstants_ml, only : AVOG
   use Setup_1dfields_ml,    only : itemp, xn => xn_2d
   use GenChemicals_ml,      only : species   ! for molwts
   use GenSpec_tot_ml,  A1 => FIRST_SOA , A2 => LAST_SOA
   implicit none

   !/-- subroutines
    public  :: OrganicAerosol


   !/-- public

    real, public, dimension(A1:A2,K1:K2), save :: Fgas  ! Fraction in gas-phase

   contains
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !+ Driver routine for Secondary Organic Aerosol  module

   subroutine OrganicAerosol(debug_flag)
   logical, intent(in) :: debug_flag  ! for debugging purposes only

     ! empty 

   end subroutine OrganicAerosol

end module OrganicAerosol_ml
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

