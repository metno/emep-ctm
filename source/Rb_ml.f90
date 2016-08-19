! <Rb_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Rb_ml

use ModelConstants_ml,    only :  DEBUG_RB
use PhysicalConstants_ml, only :  KARMAN
                    
use Wesely_ml,    only  :Wesely_tab2 &  ! Wesely Table 2 for 14 gases
                       ,Rb_cor          !
implicit none
private

public   ::  Rb_gas

    
contains
! =======================================================================

  subroutine Rb_gas(water,ustar,z0,DRYDEP_GAS,Rb) 
! =======================================================================
! Input:

    logical, intent(in) :: water
    real, intent(in)    :: ustar, z0
    integer, dimension(:), intent(in) :: &
         DRYDEP_GAS    ! Array with Wesely indices of gases wanted

! Output:

   real,dimension(:),intent(out) :: Rb   !  Rs  for dry surfaces 

! Working values:
   
    integer :: icmp             ! gaseous species
    integer :: iwes             ! gaseous species, Wesely tables

    real, parameter :: D_H2O = 0.21e-4  ! Diffusivity of H2O, m2/s
    real            :: D_i              ! Diffusivity of gas species, m2/s


! START OF PROGRAMME: 



!.........  Loop over all required gases   ................................

  GASLOOP: do icmp = 1, size( DRYDEP_GAS )
      iwes = DRYDEP_GAS(icmp)

     if   ( water ) then

          D_i = D_H2O / Wesely_tab2(1,iwes)  ! CORR !

          Rb(icmp) = log( z0 * KARMAN * ustar/ D_i )
          Rb(icmp) = Rb(icmp)/(ustar*KARMAN)

         ! CORR - Rb can be very large or even negative from this
         !        formula. We impose limits:

          Rb(icmp) = min( 1000.0, Rb(icmp) )    ! CORR - - gives min 0.001 m/s!
          Rb(icmp) = max(   10.0, Rb(icmp) )    ! CORR - - gives max 0.10 m/s!

      else

          Rb(icmp) = 2.0 * Rb_cor(iwes)/(KARMAN*ustar)
      end if


  end do GASLOOP

   if ( DEBUG_RB ) then
      print *,   "RB DRYDEP_GAS", size(DRYDEP_GAS), DRYDEP_GAS(1)
      print *,   "RB water",  water, "Rb(1) ", Rb(1)
   end if
 end subroutine Rb_gas

!--------------------------------------------------------------------

end module Rb_ml
