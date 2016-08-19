! <My_Aerosols_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
! 5. ORGANIC_AEROSOLS - for including Secondary Organic Aerosol (not active)
!----------------------------------------------------------------------

   implicit none

   !/-- public           !!  true if wanted
                    
    logical, public, parameter :: AERO_DYNAMICS     = .false.   &  
                                , EQUILIB_EMEP      = .false.  & !old Ammonium stuff
                                , EQUILIB_MARS      = .false. & !MARS
                                , EQUILIB_EQSAM     = .true. & !EQSAM
                                , ORGANIC_AEROSOLS  = .false.   
   ! Number of aerosol sizes (1-fine, 2-coarse)
    integer, public, parameter :: NSIZE = 2    

    logical, public, parameter :: SEASALT = .true. 

contains

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine My_MARS(deb)

 !..................................................................
 ! Pretty old F. Binkowski code from EPA CMAQ-Models3
 ! JGR, 108, D6, 4183
 !..................................................................

 use Setup_1dfields_ml,  only :  xn_2d     ! SIA concentration 
 use GenSpec_tot_ml,     only :  NH3, HNO3, SO4, aNO3, aNH4, NO3
 use Setup_1dfields_ml,  only :  temp, rh
 use ModelConstants_ml,  only :  KMAX_MID, KCHEMTOP   
 use GenChemicals_ml,    only :  species
 use PhysicalConstants_ml, only : AVOG
 use MARS_ml, only: rpmares

 implicit none
 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  

 logical, intent(in)  :: deb


 !.. local
  real    :: so4in, no3in, nh4in, hno3in, nh3in,   &
             aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out,   &
             coef
  integer :: k, ic, bin, spec, errmark
  logical :: debsub
 !-----------------------------------

   coef = 1.e12 / AVOG

   do k = KCHEMTOP, KMAX_MID
  
!//.... molec/cm3 -> ug/m3
      so4in  = xn_2d(SO4,k) * species(SO4)%molwt  *coef
      hno3in = xn_2d(HNO3,k)* species(HNO3)%molwt *coef 
      nh3in  = xn_2d(NH3,k) * species(NH3)%molwt  *coef
      no3in  = xn_2d(aNO3,k) * species(aNO3)%molwt  *coef
      nh4in  = xn_2d(aNH4,k) * species(aNH4)%molwt  *coef

 !--------------------------------------------------------------------------                
      call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
                    aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
                    ERRMARK,debsub) 
 !--------------------------------------------------------------------------

      xn_2d(HNO3,k)  = max (FLOOR, gNO3out / (species(HNO3)%molwt *coef) )
      xn_2d(NH3,k)   = max (FLOOR, gNH3out / (species(NH3)%molwt  *coef) )
      xn_2d(aNO3,k)  = max (FLOOR, aNO3out / (species(aNO3)%molwt  *coef) )
      xn_2d(aNH4,k)  = max (FLOOR, aNH4out / (species(aNH4)%molwt  *coef) )

   enddo  ! K-levels

 end subroutine My_MARS

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine My_EQSAM(debug_cell)

    !..................................................................
    !EQSAM - Equlibrium Simplified Aerosol Model by Swen Metzger
    !        version v03d is implemented here
    ! Metzger, S., Dentener, F., Pandis, S., and Lelieveld, J. (a): 
    !       Gas/Aerosol Partitioning 1: A computationally efficient model. 
    !       JGR, 107(D16), 10.1029/2001JD001102, 2002.
    !..................................................................

 use EQSAM_v03d_ml,      only :  eqsam_v03d
 use Setup_1dfields_ml,  only :  xn_2d     ! SIA concentration 
 use GenSpec_tot_ml,     only :  NH3, HNO3, SO4, aNO3, aNH4,NO3
 use Setup_1dfields_ml,  only :  temp, rh,pp
 use ModelConstants_ml,  only :  KMAX_MID, KCHEMTOP   
 use PhysicalConstants_ml, only : AVOG

 implicit none

 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  

 logical, intent(in)  :: debug_cell


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

  integer :: i,j,k, errmark
  logical :: debsub = .false.
 !-----------------------------------


  if ( debsub .and. debug_cell ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(aNO3,20),xn_2d(aNH4,20)
  endif

!//.... molec/cm3 -> micromoles/m**3
    so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4,KCHEMTOP:KMAX_MID)*1.e12/AVOG
    hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG
    nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    no3in(KCHEMTOP:KMAX_MID)  = xn_2d(aNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(aNH4,KCHEMTOP:KMAX_MID)*1.e12/AVOG

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
      xn_2d(aNO3,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNO3out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 ) 
      xn_2d(aNH4,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNH4out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )
      xn_2d(SO4,KCHEMTOP:KMAX_MID)   = max(FLOOR,aSO4out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )

 if ( debsub .and. debug_cell ) then ! Selected debug cell
    write(*,*)'After EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(aNO3,20),xn_2d(aNH4,20)
  endif

 end subroutine My_EQSAM


 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



 !water 
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine Aero_water(i,j)

  !.....................................................................
  ! EQSAM is called before every daily output to calculate aerosol water 
  ! at T=20C and Rh = 50%. This should model the particle water content 
  ! for gravitationally determined PM mass
  ! Tsyro, S. (2005). To what extent can aerosol water explain the 
  ! discrepancy between model calculated and gravimetric PM10 and PM2.5?. 
  ! Atmos. Chem.. Phys., 5, 602, 1-8, 2005.
  !.....................................................................

 use EQSAM_v03d_ml,      only :  eqsam_v03d
 use Setup_1dfields_ml,  only :  xn_2d      ! SIA concentration 
 use Chemfields_ml,      only :  PM_water               !PMwater  
 use GenSpec_tot_ml,     only :  NH3, HNO3, SO4, aNO3, aNH4
 use Setup_1dfields_ml,  only :  temp, rh,pp
 use ModelConstants_ml,  only :  KMAX_MID, KCHEMTOP   
 use PhysicalConstants_ml, only : AVOG

 implicit none
 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  

 integer, intent(in)  :: i, j

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

             rh50(KCHEMTOP:KMAX_MID),t20(KCHEMTOP:KMAX_MID)

  integer :: k, errmark
  logical :: debsub = .false.
 !-----------------------------------


!  if ( debsub .and. debug_cell ) then ! Selected debug cell
!    write(*,*)'Before EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
!               xn_2d(NH3,20),xn_2d(aNO3,20),xn_2d(aNH4,20)
!  endif

!//.... molec/cm3 -> micromoles/m**3
      so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
      no3in(KCHEMTOP:KMAX_MID)  = xn_2d(aNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(aNH4,KCHEMTOP:KMAX_MID)*1.e12/AVOG

      NAin(KCHEMTOP:KMAX_MID)  = 0.
      CLin(KCHEMTOP:KMAX_MID)  = 0.

      rh50(:) = 0.5
      t20(:)  = 293.

 !--------------------------------------------------------------------------                
  
   call eqsam_v03d (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh50,t20,pp,  &
                    aSO4out, aNO3out, aNH4out, aNAout, aCLout,               &
                    gSO4out, gNH3out, gNO3out, gClout, aH2Oout)
 
 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3 

      PM_water(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )

! if ( debsub .and. debug_cell ) then ! Selected debug cell
!    write(*,*)'After EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
!               xn_2d(NH3,20),xn_2d(aNO3,20),xn_2d(aNH4,20)
!  endif

 end subroutine  Aero_water
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 end module My_Aerosols_ml






