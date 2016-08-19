! <Rsurface_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Rsurface_ml

use CheckStop_ml,        only: CheckStop
use CoDep_ml,            only : CoDep_factors, RgsS_dry, RgsS_wet, &
                                 humidity_fac, Rns_NH3
use DO3SE_ml,          only : g_stomatal, do3se

use LocalVariables_ml, only : iL, L, G => Grid
  ! L (local) provides  t2C, rh, LAI, SAI, hveg, ustar, 
  !      PARsun,PARshade,LAIsunfrac, RgsO, RgsS, is_water, is_forest
  ! G (Grid)  provides snow, so2nh3ratio, 

use Radiation_ml, only : CanopyPAR

use Wesely_ml,  only  :Wesely_tab2 &  ! Wesely Table 2 for 14 gases
   ,WES_HNO3, WES_NH3,DRx    ! Indices and Ratio of diffusivities to ozone
implicit none
private

public   :: Rsurface

real, public, save :: Rinc, RigsO, GnsO
logical, private, parameter :: MY_DEBUG = .false.

 
    
contains
! =======================================================================


  subroutine Rsurface(DRYDEP_CALC,Rsur_dry,Rsur_wet,errmsg,debug_flag) 
! =======================================================================
!
!     Description
!       calculates bulk surface resistance (Rsur) for all required gases.
!
!       For O3 the methodology is derived from EMEP MSC_W Note 6/00; the 
!       following pathways apply for the surface resistance for O3:
!
!        -- Rinc-- Rgs      In-canopy + soil/ground cover
!       |
!       |
!        -------- Rext      cuticular+other external surface
!       |
!       |
!        -------- Rsto      Stomatal

!     Hence, we have a surface conductance: 
!
!       Gsur =      LAI +  SAI  +       1
!                   ___    ___     ____________
!                   Rsto   Rext    Rinc + Rgs
!
!     For SO2 and NH3 we use the CEH suggestion of a simple non-stomatal
!     uptake (which is stringly affected by wetness/RH):
!
!        -------- Rns       Non-stomatal
!       |
!        -------- Rsto      Stomatal
!
!     [ Note that the O3 formulation can be written in the same way when we
!     define GnsO = SAI/Rext + 1/(Rinc+Rgs) ]
!
!     Hence, for all gases, we have a surface conductance: 
!
!       Gsur =      LAI * Gsto  +  Gns

!  Wesely's method for other gases was based upon deriving resistances for 
!  ozone and SO2 first (e.g. RgsO, RgsS for Rgs) and then scaling using
!  effective Henry coefficients (H*) and reactivity coefficients (f0) for
!  each gas. However, here we apply scaling to Gns, not individual resistances.
!
! Structure of routine
!
!  1. Calculate:
!        Rlow             low-temperature correction
!        Rinc             in-canopy resistance
!        Rsur(HNO3)  
!        Gsto(O3)         stomatal conductance (if LAI > 0)
!
!       FOR EACH remaining gas (icmp is used as an index, since cmp is assumed 
!                               to  abbreviate "component".):
!  2. Calculate ground surface resistance, Rgs
!              if (LAI<0.1) go to 4  (for snow/ice/water...)
!  3. if (LAI>0.1)  calculate Gext
!  4. Calculate Rsur(icmp)
!       END
!
! =======================================================================

!......................................
! Input:

    integer, dimension(:), intent(in) :: &
         DRYDEP_CALC   ! Array with Wesely indices of gases wanted

! Output:

   real,dimension(:),intent(out) :: Rsur_dry   !  Rs  for dry surfaces 
   real,dimension(:),intent(out) :: Rsur_wet   !  Rs  for wet surfaces
   character(len=*), intent(out) :: errmsg
! Optional
    logical, intent(in), optional :: debug_flag


 ! external resistance for Ozone
  real, parameter :: RextO =  2500.0   ! gives Gext=0.2 cm/s for LAI=5


! Here, "Gext=0.2cm/s" refers to the external conductance, G_ext, where 
! G_ext=LAI/R_ext. In many studies, it has been assumed 
! that G_ext should be low, particularly relative to stomatal conductance g_s.
! Results from a variety of experiments, however, have made the above 
! estimates  Rext0 and RextS plausible.  The above equation for G_ext has been
! designed on the basis of these experimental results. 

! Notice also that given the equations for the canopy resistance R_sur and the 
! deposition velocity V_g, V_g>=LAI/R_ext. The value of G_ext can therefore be
! interpreted as the minimum value for V_g.


! Working values:
   
    integer :: icmp             ! gaseous species
    integer :: iwes             ! gaseous species, Wesely tables
    logical :: canopy         & ! For SAI>0, .e.g grass, forest, also in winter
        ,leafy_canopy           ! For LAI>0, only when green
    real, parameter :: SMALLSAI= 0.05  ! arbitrary value but small enough
    real :: Hstar, f0           ! Wesely tabulated Henry's coeff.'s, reactivity
    real :: Rlow                ! adjustment for low temperatures (Wesely,
                                ! 1989, p.1296, left column) 
!In Local    real :: Rinc       ! In-canopy adjustment
    real :: Rgs_dry, Rgs_wet    !  

    real ::  GnsS_dry, GnsS_wet, Gns_dry, Gns_wet


! START OF PROGRAMME: 
    errmsg = "ok"

    canopy       = ( L%SAI > SMALLSAI ) ! - can include grass
    leafy_canopy = ( L%LAI > SMALLSAI ) ! - can include grass

  !===========================================================================
  !/**  Adjustment for low temperatures (Wesely, 1989, p.1296, left column)

    Rlow = 1000.0*exp(-L%t2C - 4.0)

  !===========================================================================
  !/** get CEH humidity factor and RgsS_dry and RgsS_wet:

    call CoDep_factors(G%so2nh3ratio,L%t2C,L%rh,L%is_forest,debug_flag)


!##############   1. Calculate In-Canopy Resistance, Rinc    ################

  !/** For canopies:
  !/** Calculate stomatal conductance if daytime and LAI > 0

 
   if( leafy_canopy  .and. G%Idirect > 0.001 ) then  ! Daytime 

        call CanopyPAR(L%LAI, G%coszen, G%Idirect, G%Idiffuse, &
                    L%PARsun, L%PARshade, L%LAIsunfrac)


        call g_stomatal(iL)

   else
        L%g_sun = 0.0
        L%g_sto = 0.0

   end if ! leafy canopy and daytime

if ( MY_DEBUG .and. present(debug_flag) ) then
  if ( debug_flag )  then
    write(*,*) "IN RSUR gsto ", leafy_canopy,  G%Idirect,  L%g_sto
  end if
end if


  !/** Calculate Rinc, Gext 
  !       (use multiplication for snow, since snow=0 or 1)

   if(  canopy ) then   

         Rinc = 14.0 * L%SAI * L%hveg  / L%ustar    ! Erisman's b.LAI.h/u*

         RgsS_dry = RgsS_dry  + Rlow  + G%snow * 2000.0
         RgsS_wet = RgsS_wet  + Rlow  + G%snow * 2000.0

        ! for now, use CEH stuff for canopies, keep Ggs for non-canopy

         GnsS_dry = 1.0 /  RgsS_dry       ! For SO2, dry, low NH3 region
         GnsS_wet = 1.0 /  RgsS_wet   ! For SO2, wet, low NH3 region

   else   ! No canopy present

        Rinc = 0.0

        !/ Here we preserve the values from the ukdep_gfac table
        !  giving higher deposition to water, less to deserts

        RgsS_dry = do3se(iL)%RgsS + Rlow  + G%snow * 2000.0
        RgsS_wet = RgsS_dry    ! Hard to know what's best here
      
   end if !  canopy


!####   2. Calculate Surface Resistance, Rsur, for HNO3 and Ground Surface 
!####      Resistance, Rgs, for the remaining Gases of Interest                                

   !/ Ozone values....

     !!xRgsO  = do3se(lu)%RgsO + Rlow  + snow * 2000.0
     !!GnsO   = SAI/RextO + 1.0/( xRgsO + Rinc ) ! (SAI=0 if no canopy)
     RigsO  = Rinc + do3se(iL)%RgsO + Rlow  + G%snow * 2000.0
     GnsO   = L%SAI/RextO + 1.0/ RigsO     ! (SAI=0 if no canopy)


!.........  Loop over all required gases   ................................

  GASLOOP: do icmp = 1, size( DRYDEP_CALC )
      iwes = DRYDEP_CALC(icmp)

     !-------------------------------------------------------------------------

     !  code obtained from Wesely during 1994 personal communication
     !  but changed (ds) to allow Vg(HNO3) to exceed Vg(SO2)

        if ( iwes == WES_HNO3 ) then
            Rsur_dry(icmp)  = max(1.0,Rlow)
            Rsur_wet(icmp)  = Rsur_dry(icmp)
            cycle GASLOOP
        end if

     !-------------------------------------------------------------------------
     ! Calculate the Wesely variables Hstar (solubility) and f0 (reactivity)

        Hstar =Wesely_tab2(2,iwes)    !Extract H*'s 
        f0    =Wesely_tab2(5,iwes)    !Extract f0's
    
     !-------------------------------------------------------------------------

                          

     !   Use SAI to test for snow, ice, water, urban ...

       if ( canopy  ) then   

         ! ###   3. Calculate Cuticle conductance, Gext   ################
         ! ###      and  Ground surface conductance Ggs:

         ! Corrected for other species using Wesely's eqn. 7 approach. 
         ! (We identify leaf surface resistance with Rext/SAI.)
         ! but for conductances, not resistances (pragmatic, I know!)



         ! ##############   4. Calculate Rsur for canopies   ###############


           if ( DRYDEP_CALC(icmp) == WES_NH3 ) then

             Gns_dry = 1.0/Rns_NH3               !/** r_water  from CoDep_ml
             Gns_wet =  Gns_dry

           else  ! Not NH3

               Gns_dry = 1.0e-5*Hstar*GnsS_dry + f0 * GnsO   ! OLD SO2!
               Gns_wet = 1.0e-5*Hstar*GnsS_wet + f0 * GnsO 

             !.. and allow for partially wet surfaces at high RH, even for Gns_dry

               Gns_dry = Gns_dry * (1.0-humidity_fac) + Gns_wet * humidity_fac


           end if  ! NH3 test

           Rsur_dry(icmp) = 1.0/( L%LAI*DRx(iwes) *L%g_sto + Gns_dry  )
           Rsur_wet(icmp) = 1.0/( L%LAI*DRx(iwes) *L%g_sto + Gns_wet  )

      ! write(*,"(a20,2i3,3g12.3)")  "RSURFACE Gs  (i): ", iL, icmp, GnsO, Gns_dry, Gns_wet

       else   ! Non-Canopy modelling:

           Rgs_dry = 1.0/(1.0e-5*Hstar/RgsS_dry + f0/do3se(iL)%RgsO)  ! Eqn. (9)
           Rgs_wet = 1.0/(1.0e-5*Hstar/RgsS_wet + f0/do3se(iL)%RgsO)  ! Eqn. (9)

           Rsur_dry(icmp)   = Rgs_dry
           Rsur_wet(icmp)   = Rgs_wet
      ! write(*,"(a20,2i3,3g12.3)")  "RSURFACE Rgs (i): ", iL, icmp, Rgs_dry, Rgs_wet

       end if  ! end of canopy tests 

      ! write(*,"(a20,2i3,3g12.3)")  "RSURFACE Rsur(i): ", iL, icmp, Rsur_dry(icmp), Rsur_wet(icmp)

  end do GASLOOP


   if ( MY_DEBUG ) then
     if ( present(debug_flag) ) then
       if ( debug_flag ) then 
      write(*,*)  "RSURFACE DRYDEP_CALC", size(DRYDEP_CALC), DRYDEP_CALC(1)
      write(*,*)  "RSURFACE iL, LAI, SAI, LOGIS ", iL, L%LAI, L%SAI, &
                       L%is_forest, L%is_water, canopy, leafy_canopy
      write(*,"(a20,i3,4g12.3)")  "RSURFACE xed Gs", iL, do3se(iL)%RgsO,do3se(iL)%RgsS, Rlow, Rinc
       end if
     end if
   end if
 end subroutine Rsurface

!--------------------------------------------------------------------

end module Rsurface_ml
