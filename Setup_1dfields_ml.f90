! <Setup_1dfields_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________!
 module Setup_1dfields_ml

  ! Arrays of meteorology and concentration for 1-D column , for input to 
  ! chemical solver ........
  ! The k-dimension spans the ground (KMAX_MID) to the K-values
  ! specified by KCHEMTOP - here 2
  !
  ! - new aray added to keep o2, m, and for MADE oh, etc

  use ModelConstants_ml,     only :  KMAX_MID, KCHEMTOP, KUPPER
  use My_Emis_ml,            only :  NRCEMIS, NSS, NBVOC   !NSS=SeaS
  use GenSpec_tot_ml,        only :  NSPEC_TOT
  use GenSpec_bgn_ml,        only :  NSPEC_COL
  use GenRates_rct_ml,       only :  NRCT
  use GenRates_rcmisc_ml,    only :  NRCMISC
  implicit none
  private


  !/ variables to keep track of which call

   logical, public, save :: first_call = .true.
   integer, public, save :: ncalls     =      0

  !/-- the chemistry is calculated for arrays of size: 

   integer, public, parameter  :: CHEMSIZE = KMAX_MID-KCHEMTOP+1 !  
 
   real, public, dimension(NSPEC_TOT,KCHEMTOP:KMAX_MID), save :: &
                   xn_2d            ! Concentrations [molecules/cm3]  

   real, public, dimension(NRCEMIS,KCHEMTOP:KMAX_MID), save :: rcemis   !emissions
   real, public, dimension(NRCT   ,KCHEMTOP:KMAX_MID), save :: rct    ! T-dependant
   real, public, dimension(NRCMISC,KCHEMTOP:KMAX_MID), save :: rcmisc ! T,M,H2O-dependant
   real, public, dimension(NBVOC ,KCHEMTOP:KMAX_MID), save   :: rcbio  !  Biogenic emissions
   real, public, dimension(KCHEMTOP:KMAX_MID), save   :: rc_Rn222  ! 210Pb emissions, ds Pb210
   real, public, dimension(NSS,KCHEMTOP:KMAX_MID),     save :: rcss   ! Sea salt emissions

   real, public, dimension(KCHEMTOP:KMAX_MID), save :: &
          rh                  & ! RH (fraction, 0-1)
         ,amk                 & ! M - atmospheric conc.
         ,temp                & ! temperature
         ,pp                     !pressure
   integer, public, dimension(KCHEMTOP:KMAX_MID), save :: &
          itemp                  ! int of temperature

 end module Setup_1dfields_ml
