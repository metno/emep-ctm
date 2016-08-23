! <Rsurface_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
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
module Rsurface_ml
use LandDefs_ml,       only : LandDefs, LandType
use CheckStop_ml,      only : CheckStop
use CoDep_ml,          only : CoDep_factors, humidity_fac, Rns_NH3, Rns_SO2
use DO3SE_ml,          only : g_stomatal, do3se

use Io_Progs_ml,       only: datewrite
use LocalVariables_ml, only : iL, L, G => Grid
  ! L (local) provides  t2C, rh, LAI, SAI, hveg, ustar, 
  !      PARsun,PARshade,LAIsunfrac, RgsO, RgsS, is_water, is_forest
  ! G (Grid)  provides snow, sdepth so2nh3ratio, 

use ModelConstants_ml, only: DEBUG, NO_CROPNH3DEP
use Radiation_ml, only : CanopyPAR
use TimeDate_ml,  only : current_date
use Wesely_ml,    only : Wesely_tab2 &  ! Wesely Table 2 for 14 gases
   ,WES_HNO3, WES_NH3,DRx,WES_SO2    ! Indices and Ratio of diffusivities to ozone
use MetFields_ml, only : foundsdepth, foundice
use Par_ml,only :me
implicit none
private

public   :: Rsurface
INCLUDE 'mpif.h'
INTEGER STATUS(MPI_STATUS_SIZE),INFO


real, public, save :: Rinc, RigsO, GnsO, RgsS

 
    
contains
! =======================================================================

  subroutine Rsurface(i,j,DRYDEP_CALC,Gsto,Rsur,errmsg,debug_arg,fsnow) 
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
!        lowTcorr         low-temperature correction
!        Rinc             in-canopy resistance
!        Rsur(HNO3)  
!        Gsto(O3)         stomatal conductance (if LAI > 0)
!
!       FOR EACH remaining gas (icmp is used as an index, since cmp is assumed 
!                               to  abbreviate "component".):
!  2. Calculate ground surface resistance, Rgs
!  3. if (LAI>0.1)  calculate Gext
!  4. Calculate Rsur(icmp)
!       END
!
! =======================================================================

!......................................
! Input:
    integer, intent(in) :: i,j
    integer, dimension(:), intent(in) :: &
         DRYDEP_CALC   ! Array with Wesely indices of gases wanted

! Output:

   real,dimension(:),intent(out) :: Rsur   
   real,dimension(:),intent(out) :: Gsto

   character(len=*), intent(out) :: errmsg
! Optional
    logical, intent(in), optional :: debug_arg
    logical :: debug_flag = .false., dbg


 ! external resistance for Ozone
  real, parameter :: RextO =  2500.0   ! gives Gext=0.2 cm/s for LAI=5
  real,dimension(size(Rsur)) :: Gns   

! Here, "Gext=0.2cm/s" refers to the external conductance, G_ext, where 
! G_ext=LAI/R_ext. In many studies, it has been assumed 
! that G_ext should be low, particularly relative to stomatal conductance g_s.
! Results from a variety of experiments, however, have made the above 
! estimates  Rext0 and RextS plausible.  The above equation for G_ext has been
! designed on the basis of these experimental results. 

! Notice also that given the equations for the canopy resistance R_sur and the 
! deposition velocity V_g, V_g>=LAI/R_ext. The value of G_ext can therefore be
! interpreted as the minimum value for V_g.

    character(len=6), parameter :: sub='Rsurf:'

! Working values:
   
    integer :: icmp             ! gaseous species
    integer :: iwes             ! gaseous species, Wesely tables
    logical :: canopy         & ! For SAI>0, .e.g grass, forest, also in winter
        ,leafy_canopy           ! For LAI>0, only when green
    real, parameter :: SMALLSAI= 0.05  ! arbitrary value but small enough
    real :: Hstar, f0           ! Wesely tabulated Henry's coeff.'s, reactivity
    real :: Rgs    !  
    real :: GigsO
    real :: RsnowS, RsnowO !surface resistance for snow_flag, S and O3
    real :: lowTcorr !low temperature correction  
    real :: lowT     !low temperature correction for HNO3  
    real ::  GnsS
    real,  intent(out) :: fsnow ! the output is max(fsnow,fice)
    real :: fice !fraction ice_nwp cover
    real :: Sdmax  !max snowdepth (fsnow =1)


   if ( present(debug_arg) ) debug_flag = debug_arg
   dbg = DEBUG%RSUR .and. debug_flag

! START OF PROGRAMME: 
    errmsg = "ok"
    Sdmax = max( L%hveg/10.0, 0.01) !meters
    fsnow = 2.0 *G%sdepth/Sdmax

!Treat ice in the same way as snow
    fice=0.01*G%ice_nwp !from percent to fraction
    fsnow = max(fsnow,fice) !if snow_flag, ice_nwp probably has snow_flag
                            !but it might be ice without snow..
    fsnow = min(fsnow,1.0)  
    fsnow = max(fsnow,0.0)

    if (L%is_ice) fsnow=1.0 !ice_nwp in Landuse
                            !to ensure it is treated
                            !same way as met input
                            !and not Ggs from table
   

    if ( dbg ) call datewrite(sub//"snow_flag ",-1, &
           (/ L%hveg, Sdmax, G%ice_nwp, G%sdepth, fsnow /) )


    canopy       = ( L%SAI > SMALLSAI ) ! - can include grass
    leafy_canopy = ( L%LAI > SMALLSAI ) ! - can include grass

  !===========================================================================
  !***  Adjustment for low temperatures (Wesely, 1989, p.1296, left column)
  !     (ACP63)

    lowTcorr = exp(0.2*(-1 -L%t2C))!Zhang,2003 & Erisman 1994
    lowTcorr = min(2.0,lowTcorr)   !Zhang,2003 & Erisman 1994
    lowTcorr = max(1.0,lowTcorr)   !Zhang,2003 & Erisman 1994
                                   !effectivelluy means that it will only
                                   !kick in for T<-1

! Rsnow for sulphur and O3, Erisman, 1994 + Zhang, 2003. Also used for ice. 
    RsnowS = 70.0*(2.0 -L%t2C) !Used for snow_flag and ice_nwp
    if (L%t2C < -1.0) RsnowS = 700.0 !700 from Cadle,1985
    RsnowS = min(700.0,RsnowS) !Erisman 1994=500,very low.. Puts to 2000
    RsnowS = max(70.0,RsnowS)  !Erisman 1994. 70 above 1 degree

    RsnowO = 2000.0 !same for snow_flag, ice_nwp, water. Later corrected with lowTcorr
                    !as recommended by Juha-Pekka


  !===========================================================================
  ! Get Rns_SO2

       call CoDep_factors(G%so2nh3ratio24hr,G%so2nh3ratio,&
              L%t2C,L%rh,L%is_forest, debug_flag)


!##############   1. Calculate In-Canopy Resistance, Rinc    ################

  !*** For canopies:
  !*** Calculate stomatal conductance if daytime and LAI > 0 and snowdepth
  !    less than 1m above vegetation (1m chosen arbitrary)

   !g_sto 0 when snow covering canopy
   if ( dbg ) then
      write(*,*) sub//"testcan ", iL, canopy, leafy_canopy
      call datewrite(sub//"testcan ", iL, (/ G%Idirect, L%PARsun, G%sdepth /) )
   end if
   if( leafy_canopy  .and. G%Idirect > 0.001 .and. G%sdepth< (1.0 +Sdmax) ) then  ! Daytime 

        call CanopyPAR(L%LAI, G%coszen, G%Idirect, G%Idiffuse, &
                    L%PARsun, L%PARshade, L%LAIsunfrac)


        call g_stomatal(iL, debug_flag )

        if ( DEBUG%RSUR .and. debug_flag ) then
           call datewrite(sub//"gstoA ", iL, &
               (/ G%Idirect+G%Idiffuse, L%PARsun,  L%g_sto, L%f_env /) )
        end if

   else
        L%g_sun = 0.0
        L%g_sto = 0.0
        L%f_env = 0.0 !JAN2013
        L%PARsun  = 0.0 !FEB2013
        L%PARshade = 0.0 !FEB2013

   end if ! leafy canopy and daytime

   if ( dbg ) write(*,"(a,5i5,i3,L2,3f10.4)") "IN RSUR gsto ", &
              current_date, iL, leafy_canopy, G%Idirect,  L%g_sto, L%f_env


  !*** Calculate Rinc, Gext   (ACPs8.6.1)

     if(  canopy ) then   

         Rinc = 14.0 * L%SAI * L%hveg  / L%ustar    ! Erisman's b.LAI.h/u*

        ! for now, use CEH stuff for canopies,and soils (canopies ouside 
        ! growing season)
        ! keep Ggs for non-canopy

        GnsS = (1.-fsnow)/(Rns_SO2 * lowTcorr) + fsnow/RsnowS 
        RgsS = 1./GnsS

  
     elseif  ( L%is_veg ) then ! vegetation outside growing season

        Rinc = 0.0
        GnsS = (1.-fsnow)/(Rns_SO2 * lowTcorr) + fsnow/RsnowS 
        RgsS = 1./GnsS


     else   ! No canopy or soil present

        Rinc = 0.0

        !/ Here we preserve the values from the ukdep_gfac table
        !  giving higher deposition to water, less to deserts

        GnsS = (1.-fsnow)/(do3se(iL)%RgsS * lowTcorr) + fsnow/RsnowS 
        RgsS = 1./GnsS

     end if !  canopy

        !snow treated similar to Zhang 2003
        !But Zhang wse 2*fsnow for ground surface because Sdmax(snow depth when total coverage is assumed)
        !for soils under vegetation is assumed to stay snow covered longer than 'the leafs'
        !but - we have underlying surfaces only for O3 and for simplicity we treat them equally
        !RECONSIDER THIS ESPECIALLY BASED ON SATELITTES

        !no snow corrections (or low temperature) for Rinc 
        !RgsO 'corrected for snow' and low temp  (JP)

        GigsO=  (1.-fsnow)/do3se(iL)%RgsO   + fsnow/RsnowO
        RigsO = lowTcorr/GigsO +  Rinc


!####   2. Calculate Surface Resistance, Rsur, for HNO3 and Ground Surface 
!####      Resistance, Rgs, for the remaining Gases of Interest                                

   !/ Ozone values....

      !RextO corrected for low temp (JP)

     GnsO   = L%SAI/(RextO * lowTcorr) + 1.0/ RigsO     ! (SAI=0 if no canopy)


!.........  Loop over all required gases   ................................

  GASLOOP: do icmp = 1, size( DRYDEP_CALC )
      iwes = DRYDEP_CALC(icmp)

     !-------------------------------------------------------------------------

     !  code obtained from Wesely during 1994 personal communication
     !  but changed to allow Vg(HNO3) to exceed Vg(SO2)
     !  lowT based on: Rc=10scm-1 for -5, Rc=50scm-1 at -18 in Johanson&Granat, 1986

        if ( iwes == WES_HNO3 ) then
            lowT= -L%t2C *2.0
            Rsur(icmp)  = max(10.0,lowT) !not so affected by snow_flag, e.g. Erisman 1994,table 6,
            Gsto(icmp) = 0.0              !only for code consistency eleswhere
            !Gns(icmp) = 1.0/Rsur(icmp)   !only for code consistency eleswhere
           ! - reimplement Vg limitation for HNO3. 10 cm/s max is enough anyway!
           ! Rsur(icmp)  = max(1.0,lowT) !Cadle,1985

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

               Gns(icmp) = (1.-fsnow)/(Rns_NH3 * lowTcorr) + fsnow/RsnowS 
           else

               Gns(icmp) = 1.0e-5*Hstar*GnsS + f0 * GnsO   ! OLD SO2!
           end if


           Gsto(icmp) = L%LAI*DRx(iwes) *L%g_sto
           Rsur(icmp) = 1.0/( L%LAI*DRx(iwes) *L%g_sto + Gns(icmp)  )

         ! Stop NH3 deposition for growing crops 
         ! Crude reflection of likely emission

           if ( NO_CROPNH3DEP .and. DRYDEP_CALC(icmp) == WES_NH3 ) then

              if ( L%is_crop .and.  L%LAI > 0.1 ) then
                   if ( DEBUG%RSUR .and. debug_flag .and. L%is_crop ) then 
                      write(*,"(a,i4,2i4,L2,f8.2)")  "NO_CROPNH3DEP ", &
                       iL, DRYDEP_CALC(icmp), WES_NH3, L%is_crop, L%LAI
                   end if

                 Rsur(icmp) =  1.0e10  ! BIG number
                 Gsto(icmp) =  0.0     ! just for code consistecny

              end if
           end if

      ! write(*,"(a20,2i3,3g12.3)")  "RSURFACE Gs  (i): ", iL, icmp, GnsO, Gns_dry, Gns_wet

      elseif (L%is_veg) then !vegetation outside growing season

           if ( DRYDEP_CALC(icmp) == WES_NH3 ) then

               Gns(icmp) = (1.-fsnow)/(Rns_NH3 * lowTcorr) + fsnow/RsnowS 
           else 
               Gns(icmp) = 1.0e-5*Hstar*GnsS + f0 * GnsO   ! OLD SO2!
           end if

           Rsur(icmp) = 1.0/Gns(icmp)  
           Gsto(icmp) =  0.0     ! just for code consistecny

      else   ! Non-Canopy modelling:
           Gns(icmp) = 1.0e-5*Hstar*GnsS + f0*GnsO
           Rgs = 1.0/Gns(icmp)  ! Eqn. (9) !hf was  f0/do3se(iL)%RgsO
           Rsur(icmp)   = Rgs
           Gsto(icmp) =  0.0     ! just for code consistecny

      end if  ! end of canopy tests 


      if ( dbg ) write(*,"(a20,2i3,L2,3g12.3)")  &
             "RSURFACE Rsur(i): ", iL, icmp, L%is_crop,  Rsur(icmp)


  end do GASLOOP


  if ( dbg ) then 
      write(*,"(a,2i4)")  "RSURFACE DRYDEP_CALC", &
            size(DRYDEP_CALC), DRYDEP_CALC(1)
      write(*,"(a,i3,2f7.3,5L2)")  "RSURFACE iL, LAI, SAI, LOGIS ", &
            iL, L%LAI, L%SAI, L%is_forest, L%is_water, L%is_veg, &
              canopy, leafy_canopy
      write(*,"(a,i3,4g12.3)")  "RSURFACE xed Gs", iL, &
             do3se(iL)%RgsO,do3se(iL)%RgsS, lowTcorr, Rinc
  end if
 end subroutine Rsurface

!--------------------------------------------------------------------

end module Rsurface_ml
