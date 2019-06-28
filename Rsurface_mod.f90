! <Rsurface_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2019 met.no
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
module Rsurface_mod
use CheckStop_mod,      only: CheckStop, StopAll
use CoDep_mod,          only: CoDep_factors, humidity_fac, Rns_NH3, Rns_SO2
use Config_module,      only: NO_CROPNH3DEP
use Debug_module,       only: DEBUG   ! -> DEBUG%RSUR
use DO3SE_mod,          only: g_stomatal, do3se
use GasParticleCoeffs_mod, only: nddep, DDspec,  &
                              idcmpO3, idcmpHNO3,idcmpNH3,idcmpSO2
use Io_Progs_mod,       only: datewrite
use LandDefs_mod,       only: LandDefs, LandType
! L (local) provides  t2C, rh, LAI, SAI, hveg, ustar, 
!      PARsun,PARshade,  (in W/m2)
!      LAIsunfrac, RgsO, RgsS, is_water, is_forest
! G (Grid)  provides snow, sdepth so2nh3ratio, 
use LocalVariables_mod, only : iL, L, G => Grid

use MetFields_mod, only : foundsdepth, foundice
use MetFields_mod, only : PARdbh, PARdif  ! W/m2
use Par_mod,only :me
use Radiation_mod, only : CanopyPAR
use SmallUtils_mod, only: find_index
use TimeDate_mod,  only : current_date

implicit none
private

public   :: Rsurface
!INCLUDE 'mpif.h'

real, public, save :: Rinc, RigsO, GnsO, RgsS

 
    
contains
! =======================================================================

  subroutine Rsurface(i,j,Gsto,Rsur,errmsg,debug_arg,fsnow) 
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
!     uptake (which is strongly affected by wetness/RH):
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
!        Rsur(NH3)  
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

    character(len=6), parameter :: dtxt='Rsurf:'

! Working values:
   
    integer :: icmp             ! gaseous species, =index in DDspec list
    logical :: canopy         & ! For SAI>0, .e.g grass, forest, also in winter
        ,leafy_canopy           ! For LAI>0, only when green
    real, parameter :: SMALLSAI= 0.05  ! arbitrary value but small enough
    real :: Hstar, f0           ! DDdefs tabulated Henry's coeff.'s, reactivity
    real :: Rgs
    real :: GigsO
    real :: RsnowS, RsnowO !surface resistance for snow_flag, S and O3
    real :: lowTcorr    !low temperature correction  
    real :: lowThno3    !low temperature correction for HNO3  
    real :: GnsS, RnsS  ! for SO2
    real,  intent(out) :: fsnow ! the output is max(fsnow,fice)
    real :: fice !fraction ice_nwp cover
    real :: Sdmax  !max snowdepth (fsnow =1)


   if ( present(debug_arg) ) debug_flag = debug_arg
   dbg = DEBUG%RSUR .and. debug_flag


! START OF PROGRAMME: 
    errmsg = "ok"
    Sdmax = max( L%hveg/10, 0.01) !meters
    fsnow = 2 *G%sdepth/Sdmax

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
   

    if ( dbg ) call datewrite(dtxt//"snow_flag ",-1, &
           (/ L%hveg, Sdmax, G%ice_nwp, G%sdepth, fsnow /) )


    canopy       = ( L%SAI > SMALLSAI ) ! - can include grass
    leafy_canopy = ( L%LAI > SMALLSAI )

  !===========================================================================
  !***  Adjustment for low temperatures (Wesely, 1989, p.1296, left column)
  !     (ACP63)

    lowTcorr = exp(0.2*(-1 -L%t2C))!Zhang,2003 & Erisman 1994
    lowTcorr = min(2.0,lowTcorr)   !Zhang,2003 & Erisman 1994
    lowTcorr = max(1.0,lowTcorr)   !Zhang,2003 & Erisman 1994
                                   !effectively means that it will only
                                   !kick in for T<-1

  ! Rsnow for sulphur and O3, Erisman, 1994 + Zhang, 2003. Also used for ice. 
    RsnowS = 70.0*(2.0 -L%t2C) !Used for snow_flag and ice_nwp
    if (L%t2C < -1.0) RsnowS = 700.0 !700 from Cadle,1985
    RsnowS = min(700.0,RsnowS) !Erisman 1994=500,very low.. Puts to 2000
    RsnowS = max(70.0,RsnowS)  !Erisman 1994. 70 above 1 degree

    RsnowO = 2000.0 !same for snow_flag, ice_nwp, water. Later corrected with
                    ! lowTcorr as recommended by Juha-Pekka

  !===========================================================================
  ! Get Rns_SO2, Rns_NH3 accounting for co-deposition
  ! (FUTURE: Rns_NH3 will be overwritten later in DryDep_mod if Bi-Directional
  !  triggered)

    call CoDep_factors(G%so2nh3ratio24hr,G%so2nh3ratio,&
              L%t2C,L%rh,L%is_forest, debug_flag)

  !##############   1. Calculate In-Canopy Resistance, Rinc    ################

  !*** For canopies:
  !*** Calculate stomatal conductance if daytime and LAI > 0 and snowdepth
  !    less than 1m above vegetation (1m chosen arbitrary)
  !    g_sto 0 when snow covering canopy

   if ( dbg ) then
     call datewrite(dtxt//"testcan ", iL, [ G%Idirect, L%PARsun, G%sdepth ] )
   end if
 
   if( leafy_canopy  .and. G%Idirect > 0.001 .and.       &!: daytime
                          G%sdepth< (1.0 +Sdmax) ) then   !: above snow 

     call CanopyPAR(L%LAI, G%coszen, PARdbh(i,j), PARdif(i,j), &
                    L%PARsun, L%PARshade, L%LAIsunfrac)

     call g_stomatal(iL, debug_flag )

   else
     L%g_sun = 0.0
     L%g_sto = 0.0
     L%f_env = 0.0
     L%PARsun  = 0.0
     L%PARshade = 0.0
   end if ! leafy canopy and daytime

   if ( dbg ) call datewrite(dtxt//" gsto ", iL, &
       [ L%LAI, G%Idirect, G%Idiffuse, L%PARsun,  L%g_sto, L%f_env ] )


  !*** Calculate in-canopy resistance,  Rinc (ACPs8.6.1)
  ! no snow corrections (or low temperature) for Rinc 

   if  ( canopy ) then ! vegetation inside growing season
      Rinc = 14 * L%SAI * L%hveg  / L%ustar    ! Erisman's b.LAI.h/u*
   else
      Rinc = 0.0
   end if


  ! Sulphur, non-stomatal terms , ACP 59
  ! (fsnow already accounts for factor 2 from ACP 59)
  ! If no canopy or soil present we preserve the values from the original
  ! table, giving higher deposition to water, less to deserts for now,
  ! use CEH stuff for canopies,and soils (canopies ouside growing season)
  ! keep Ggs for non-canopy With canopy we allow co-dep effect on Rns_SO2

   RnsS = do3se(iL)%RgsS
   if(  canopy .or. L%is_veg ) RnsS = Rns_SO2

   GnsS = (1.-fsnow)/(RnsS * lowTcorr) + fsnow/RsnowS 
!TMP DEC7 CHECK    GnsS = min( 0.1, GnsS ) ! OCT2017 FIX
   RgsS = 1./GnsS


  !/ Ozone values...., ACP 60
  ! snow treated similar to Zhang 2003
  ! But Zhang use 2*fsnow for ground surface because Sdmax(snow depth
  ! when total coverage is assumed) for soils under vegetation is
  ! assumed to stay snow covered longer than 'the leafs' but - we have
  ! underlying surfaces only for O3 and for simplicity we treat them
  ! equally RECONSIDER THIS ESPECIALLY BASED ON SATELITTES

  !RgsO 'corrected for snow' and low temp  (JP)
  !RextO corrected for low temp (JP)

   GigsO  =  (1.-fsnow)/do3se(iL)%RgsO + fsnow/RsnowO
   RigsO  = lowTcorr/GigsO +  Rinc
   GnsO   = L%SAI/(RextO * lowTcorr) + 1.0/ RigsO     ! (SAI=0 if no canopy)



  !####   2. Calculate Surface Resistance, Rsur, for HNO3 and Ground Surface ##
  !####      Resistance, Rgs, for the remaining Gases of Interest            ##
  !.........  Loop over all required gases   ................................##

  GASLOOP: do icmp = 1, nddep ! size( DRYDEP_CALC )
      Gsto(icmp) = 0.0                     ! change where needed

      if ( .not. DDspec(icmp)%is_gas ) CYCLE

     !-------------------------------------------------------------------------
     ! HNO3
     !  code obtained from Wesely during 1994 personal communication
     !  but changed to allow Vg(HNO3) to exceed Vg(SO2)
     !  lowThno3 based on: Rc=10s/cm for -5, Rc=50s/cm at -18 in Johanson&Granat, 1986
     ! - reimplement Vg limitation for HNO3. 10 cm/s max is enough anyway!
     !  Earlier had Rsur(icmp)  = max(1.0,lowThno3) !Cadle,1985
     !  Also, not so affected by snow_flag, e.g. Erisman 1994,Table 6,

        if ( icmp == idcmpHNO3 ) then
            lowThno3= -L%t2C *2.0
            Rsur(icmp)  = max(10.0,lowThno3) 
            Gns(icmp)   =  1/Rsur(icmp) ! OCT2017 if needed???!
            cycle GASLOOP
        end if

     !-------------------------------------------------------------------------
     ! Calculate the Wesely variables Hstar (solubility) and f0 (reactivity)

        Hstar =DDspec(icmp)%Hstar
        f0    =DDspec(icmp)%f0

     !-------------------------------------------------------------------------
     ! Ammonia is also special
     ! Has just Gsto and Gns. Uses Rns from co-dep or (FUTURE) BiDir
     ! QUERY - no impact of wet soils?

        if ( icmp == idcmpNH3 ) then

          if  (canopy .or. L%is_veg ) then
            Gns(icmp) = (1.-fsnow)/(Rns_NH3 * lowTcorr) + fsnow/RsnowS 
            Gsto(icmp) = L%LAI*DDspec(icmp)%DxDO3 *L%g_sto
          else !OLD
            Gns(icmp) = 1.0e-5*Hstar*GnsS + f0*GnsO
          end if
          Gns(icmp) = min( 0.1, Gns(icmp) ) ! FIX
          if ( Gns(icmp) > 0.1 ) then
            print *, dtxt//"BIGGNS!",canopy,L%is_veg,Rns_NH3,lowTcorr,RsnowS
            call StopAll('BIGGNS-NH3')
          end if

          Rsur(icmp) = 1.0/( Gsto(icmp) + Gns(icmp)  )

        ! No NH3 dep for growing crops - crude reflection of likely emis

          if ( NO_CROPNH3DEP  ) then
             if ( L%is_crop .and.  L%LAI > 0.1 ) then
               if ( dbg .and. L%is_crop ) then 
                write(*,"(a,i4,L2,f8.2)")  "NO_CROPDEP:"// &
                  DDspec(icmp)%name,  iL, L%is_crop, L%LAI
               end if
               Rsur(icmp) =  1.0e10  ! BIG number
             end if ! is_crop
          end if

          cycle GASLOOP
        end if  ! NH3

    
     !-------------------------------------------------------------------------
     ! ###   3. Calculate Cuticle conductance, Gext   ################
     ! ###      and  Ground surface conductance Ggs:

     ! Corrected for other species using Wesely's eqn. 7 approach. 
     ! (We identify leaf surface resistance with Rext/SAI.)
     ! but for conductances, not resistances (pragmatic, I know!)

       Gns(icmp) = 1.0e-5*Hstar*GnsS + f0 * GnsO 


     ! ##############   4. Calculate Rsur for canopies   ###############

      if ( canopy  ) then   
         Gsto(icmp) = L%LAI*DDspec(icmp)%DxDO3 *L%g_sto
      end if

!WHY?  Rgs = 1.0/Gns(icmp)  ! Eqn. (9) !hf was  f0/do3se(iL)%RgsO

      Rsur(icmp) = 1.0/( Gsto(icmp) + Gns(icmp)  )

      if ( dbg ) write(*,"(a,i3,L2,99g10.2)")  &
        dtxt//" Rsur(i):"//trim(DDspec(icmp)%name)//' '//&
           trim(LandDefs(iL)%name), icmp, L%is_crop,&
           1.0e-5*Hstar, GnsS, f0, GnsO, Gsto(icmp),Gns(icmp), Rsur(icmp)

  end do GASLOOP


  if ( dbg ) write(*,"(a,a10,i4,2f7.3,5L2)") &
    dtxt//" nGas iL, LAI, SAI, LOGIS ", trim(LandDefs(iL)%name), nddep, L%LAI,&
       L%SAI, L%is_forest, L%is_water, L%is_veg, canopy, leafy_canopy
 end subroutine Rsurface

!--------------------------------------------------------------------

end module Rsurface_mod
