! <My_Aerosols_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-201409 met.no
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
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                           module My_Aerosols_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!----------------------------------------------------------------------
! Allows to select aerosol types for the model run:
! 1. AERO_DYNAMICS - for running UNI-AERO (presently not included)
! Options for aeroso-gas equilibrium partitioning:
! 2. EQUILIB_EMEP - old EMEP scheme
! 3. EQUILIB_MARS - run MARS equilibrium model
! 4. EQUILIB_EQSAM - run EQSAM equilibrium model
!----------------------------------------------------------------------

 use CheckStop_ml,         only :  StopAll, CheckStop
 use ChemSpecs,            only :  SO4, NH3, HNO3, NO3_f, NH4_f, &
                                     NSPEC_SHL, species
!CMR  use ChemSpecs_tot_ml,     only :  SO4, NH3, HNO3, NO3_f, NH4_f
!CMR  use ChemSpecs_shl_ml,     only :  NSPEC_SHL
!CMR  use ChemChemicals_ml,     only :  species
 use Chemfields_ml,        only :  PM25_water, PM25_water_rh50,  & !PMwater 
                                   cfac
 use EQSAM_v03d_ml,        only :  eqsam_v03d
 use MARS_ml,              only :  rpmares
 use ModelConstants_ml,    only :  KMAX_MID, KCHEMTOP, DEBUG
 use PhysicalConstants_ml, only :  AVOG
 use Setup_1dfields_ml,    only :  xn_2d,       & ! SIA concentration 
                                   temp, rh, pp


 implicit none

   !/-- public           !!  true if wanted
                    
    logical, public, parameter :: AERO_DYNAMICS     = .false.  &  
                                , EQUILIB_EMEP      = .false.  & !old Ammonium stuff
                                , EQUILIB_MARS      = .true.  & !MARS
                                , EQUILIB_EQSAM     = .false.     !EQSAM
                                
 !.. Number of aerosol sizes (1-fine, 2-coarse, 3-'giant' for sea salt )
    integer, public, parameter :: NSIZE = 5
           !   FINE_PM = 1, COAR_NO3 = 2, COAR_SS = 3, COAR DUST = 4,pollen = 5    


contains

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine My_MARS(debug_flag)

 !..................................................................
 ! Pretty old F. Binkowski code from EPA CMAQ-Models3
 ! JGR, 108, D6, 4183
 !..................................................................

 logical, intent(in) :: debug_flag 
 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  

 !.. local
  real    :: so4in, no3in, nh4in, hno3in, nh3in,   &
             aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out,   &
             coef
  integer :: k, errmark
 !-----------------------------------

   coef = 1.e12 / AVOG

   do k = KCHEMTOP, KMAX_MID
  
!//.... molec/cm3 -> ug/m3
      so4in  = xn_2d(SO4,k) * species(SO4)%molwt  *coef
      hno3in = xn_2d(HNO3,k)* species(HNO3)%molwt *coef 
      nh3in  = xn_2d(NH3,k) * species(NH3)%molwt  *coef
      no3in  = xn_2d(NO3_f,k) * species(NO3_f)%molwt  *coef
      nh4in  = xn_2d(NH4_f,k) * species(NH4_f)%molwt  *coef

 !--------------------------------------------------------------------------                
      call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
                    aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
                    ERRMARK,debug_flag) 
 !--------------------------------------------------------------------------

      if( DEBUG%EQUIB) then
        call CheckStop(gNO3out< 0.0, "XMARS: gNO3out")
        call CheckStop(gNH3out< 0.0, "XMARS: gNH3out")
        call CheckStop(aNO3out< 0.0, "XMARS: aNO3out")
        call CheckStop(aNH4out< 0.0, "XMARS: aNH4out")
      end if ! DEBUG%EQUIB

      xn_2d(HNO3,k)  = max (FLOOR, gNO3out / (species(HNO3)%molwt *coef) )
      xn_2d(NH3,k)   = max (FLOOR, gNH3out / (species(NH3)%molwt  *coef) )
      xn_2d(NO3_f,k)  = max (FLOOR, aNO3out / (species(NO3_f)%molwt  *coef) )
      xn_2d(NH4_f,k)  = max (FLOOR, aNH4out / (species(NH4_f)%molwt  *coef) )

   enddo  ! K-levels

 end subroutine My_MARS

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine My_EQSAM(debug_flag)
 logical, intent(in) :: debug_flag 

    !..................................................................
    !EQSAM - Equlibrium Simplified Aerosol Model by Swen Metzger
    !        version v03d is implemented here
    ! Metzger, S., Dentener, F., Pandis, S., and Lelieveld, J. (a): 
    !       Gas/Aerosol Partitioning 1: A computationally efficient model. 
    !       JGR, 107(D16), 10.1029/2001JD001102, 2002.
    !..................................................................


 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  


 !.. local
  real    :: so4in(KCHEMTOP:KMAX_MID),   &
             no3in(KCHEMTOP:KMAX_MID),   &
             nh4in(KCHEMTOP:KMAX_MID),   &
             hno3in(KCHEMTOP:KMAX_MID),  &
             nh3in(KCHEMTOP:KMAX_MID),   &
! The following input is not in use
             NAin(KCHEMTOP:KMAX_MID),    &
             CLin(KCHEMTOP:KMAX_MID),    &

             aSO4out(KCHEMTOP:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), & 
             gNH3out(KCHEMTOP:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID), &
             aNAout(KCHEMTOP:KMAX_MID),  &
             aCLout(KCHEMTOP:KMAX_MID),  &
             gCLout(KCHEMTOP:KMAX_MID),  &
             gSO4out(KCHEMTOP:KMAX_MID)

 !-----------------------------------


  if ( debug_flag  ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(NO3_f,20),xn_2d(NH4_f,20)
  endif

!//.... molec/cm3 -> micromoles/m**3
    so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4,KCHEMTOP:KMAX_MID)*1.e12/AVOG
    hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG
    nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    no3in(KCHEMTOP:KMAX_MID)  = xn_2d(NO3_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(NH4_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG

    NAin(KCHEMTOP:KMAX_MID)  = 0.0
    CLin(KCHEMTOP:KMAX_MID)  = 0.0

 !--------------------------------------------------------------------------                
  
    call eqsam_v03d (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh,temp,pp, &
                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,             &
                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout)
 
 !--------------------------------------------------------------------------

!//.... micromoles/m**3  -> molec/cm3 
!      xn_2d(NO3,KCHEMTOP:KMAX_MID)  = FLOOR !different for ACID/OZONE

      xn_2d(HNO3,KCHEMTOP:KMAX_MID)  = max(FLOOR,gNO3out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )
      xn_2d(NH3,KCHEMTOP:KMAX_MID)   = max(FLOOR,gNH3out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )
      xn_2d(NO3_f,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNO3out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 ) 
      xn_2d(NH4_f,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNH4out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )
      xn_2d(SO4,KCHEMTOP:KMAX_MID)   = max(FLOOR,aSO4out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )

 if ( debug_flag ) then ! Selected debug cell
    write(*,*)'After EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(NO3_f,20),xn_2d(NH4_f,20)
  endif

 end subroutine My_EQSAM


 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



 !water 
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine Aero_water(i,j, ambient, debug_flag)

  !.....................................................................
  ! EQSAM is called before every daily output to calculate aerosol water 
  ! at T=20C and Rh = 50%. This should model the particle water content 
  ! for gravitationally determined PM mass
  ! Tsyro, S. (2005). To what extent can aerosol water explain the 
  ! discrepancy between model calculated and gravimetric PM10 and PM2.5?. 
  ! Atmos. Chem.. Phys., 5, 602, 1-8, 2005.
  !.....................................................................

 implicit none

 integer, intent(in)  :: i, j
 logical, intent(in)  :: ambient
 logical, intent(in)  :: debug_flag
 !.. local
  real    :: so4in(KCHEMTOP:KMAX_MID),   &
             no3in(KCHEMTOP:KMAX_MID),   &
             nh4in(KCHEMTOP:KMAX_MID),   &
             hno3in(KCHEMTOP:KMAX_MID),  &
             nh3in(KCHEMTOP:KMAX_MID),   &
! fix
             NAin(KCHEMTOP:KMAX_MID)  ,  &
             CLin(KCHEMTOP:KMAX_MID) ,   &
             aSO4out(KCHEMTOP:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), & 
             gNH3out(KCHEMTOP:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID), &
             aNAout(KCHEMTOP:KMAX_MID),  &
             aCLout(KCHEMTOP:KMAX_MID),  &
             gCLout(KCHEMTOP:KMAX_MID),  &
             gSO4out(KCHEMTOP:KMAX_MID), &
             rlhum(KCHEMTOP:KMAX_MID),tmpr(KCHEMTOP:KMAX_MID)

  real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  
 !-----------------------------------


  if ( debug_flag ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(NO3_f,20),xn_2d(NH4_f,20)
  endif

!//.... molec/cm3 -> micromoles/m**3
      so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
      no3in(KCHEMTOP:KMAX_MID)  = xn_2d(NO3_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(NH4_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG

!NEEDS WORK: NAin(KCHEMTOP:KMAX_MID)   = xn_2d(SEASALT_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      call StopAll("Sea-salt not implemented in eqsam in Nov 2011 plus versions")
      CLin(:) =  NAin(:)
!      NAin(KCHEMTOP:KMAX_MID)  = 0.
!      CLin(KCHEMTOP:KMAX_MID)  = 0.

      if ( ambient ) then               ! real LWC
                  rlhum(:) = rh(:)
                  tmpr(:)  = temp(:)
      else                              ! for gravimetric mass
                  rlhum(:) = 0.5
                  tmpr(:)  = 293.15
      endif

 !--------------------------------------------------------------------------                
  
   call eqsam_v03d (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin,rlhum,tmpr,pp, &
                    aSO4out, aNO3out, aNH4out, aNAout, aCLout,               &
                    gSO4out, gNH3out, gNO3out, gClout, aH2Oout)
 
 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 
 if (ambient)  then      ! at ambient conditions (3D)
      PM25_water(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )
 else                    ! In gravimetric PM (Rh=50% and t=20C)
      PM25_water_rh50 (i,j)             = max(0., aH2Oout(KMAX_MID) )
 endif

 if ( debug_flag ) then ! Selected debug cell
    write(*,*)'After EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(NO3_f,20),xn_2d(NH4_f,20)
  endif

 end subroutine  Aero_water
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine Aero_water_MARS(i,j, debug_flag)

 !..................................................................
 ! Pretty old F. Binkowski code from EPA CMAQ-Models3
 ! JGR, 108, D6, 4183
 !..................................................................

 integer, intent(in)  :: i, j
 logical, intent(in)  :: debug_flag
 real, parameter      :: FLOOR = 1.0E-30  ! minimum concentration  

 !.. local
  real    :: rlhum(KCHEMTOP:KMAX_MID), tmpr(KCHEMTOP:KMAX_MID)
  real    :: so4in, no3in, nh4in, hno3in, nh3in,   &
             aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out,   &
             coef
  integer :: k, errmark
 !-----------------------------------

   coef = 1.e12 / AVOG

 !.. PM2.5 water at ambient conditions (3D)
   rlhum(:) = rh(:) 
   tmpr(:)  = temp(:)

    do k = KCHEMTOP, KMAX_MID
  
!//.... molec/cm3 -> ug/m3
      so4in  = xn_2d(SO4,k) * species(SO4)%molwt  *coef
      hno3in = xn_2d(HNO3,k)* species(HNO3)%molwt *coef 
      nh3in  = xn_2d(NH3,k) * species(NH3)%molwt  *coef
      no3in  = xn_2d(NO3_f,k) * species(NO3_f)%molwt  *coef
      nh4in  = xn_2d(NH4_f,k) * species(NH4_f)%molwt  *coef


 !--------------------------------------------------------------------------                
      call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
                    aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
                    ERRMARK,debug_flag) 
 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 
      PM25_water(i,j,k) = max (0., aH2Oout )

    enddo  ! k-loop

!.. PM2.5 water at equilibration conditions for gravimetric PM (Rh=50% and t=20C)
                            
    rlhum(:) = 0.5
    tmpr(:)  = 293.15
    k = KMAX_MID
!//.... molec/cm3 -> ug/m3
      so4in  = xn_2d(SO4,k) * species(SO4)%molwt  *coef *cfac(SO4-NSPEC_SHL,i,j) 
      hno3in = xn_2d(HNO3,k)* species(HNO3)%molwt *coef *cfac(HNO3-NSPEC_SHL,i,j)
      nh3in  = xn_2d(NH3,k) * species(NH3)%molwt  *coef *cfac(NH3-NSPEC_SHL,i,j)
      no3in  = xn_2d(NO3_f,k) * species(NO3_f)%molwt *coef *cfac(NO3_f-NSPEC_SHL,i,j)
      nh4in  = xn_2d(NH4_f,k) * species(NH4_f)%molwt *coef *cfac(NH4_f-NSPEC_SHL,i,j)
!--------------------------------------------------------------------------                
      call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
                    aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
                    ERRMARK,debug_flag) 
 !--------------------------------------------------------------------------

      PM25_water_rh50 (i,j) = max (0., aH2Oout )


 end subroutine  Aero_water_MARS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 end module My_Aerosols_ml






