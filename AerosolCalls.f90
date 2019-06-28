! <AerosolCalls.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
!>  AerosolCalls.f90 - A component of the EMEP MSC-W Chemical transport Model
!!****************************************************************************! 
!> Options for aerosol-gas equilibrium partitioning:
!!
!! * EMEP - old EMEP scheme with (NH4)1.5 SO4
!! * MARS - run MARS equilibrium model
!! * EQSAM - run EQSAM equilibrium model
!! * ISORROPIA - run ISORROPIA II equilibrium model (in testing)

module AerosolCalls

 use AeroConstants_mod,     only: AERO
 use Ammonium_mod,          only: Ammonium
 use CheckStop_mod,         only: StopAll, CheckStop
 use ChemDims_mod,          only: NSPEC_SHL
 use ChemSpecs_mod
 use Chemfields_mod,        only: PM25_water, PM25_water_rh50, & !H2O_eqsam, & !PMwater 
                                   cfac
 use Config_module,         only: KMAX_MID, KCHEMTOP, MasterProc
 use Debug_module,          only: DEBUG   ! -> DEBUG%EQUIB
 use EQSAM4clim_ml,        only :  EQSAM4clim
! use EQSAM_v03d_mod,        only: eqsam_v03d
 use MARS_mod,              only: rpmares, rpmares_2900, DO_RPMARES_new
 use PhysicalConstants_mod, only: AVOG
 use ZchemData_mod,         only: xn_2d, temp, rh, pp
 implicit none
 private

 !/-- public           !!  true if wanted

 public :: AerosolEquilib
 public :: emep2MARS, emep2EQSAM, Aero_Water, Aero_Water_rh50, Aero_Water_MARS
 !FUTURE private :: emep2isorropia
                    
!    logical, public, parameter :: AERO_DYNAMICS     = .false.  &  
!                                , EQUILIB_EMEP      = .false.  & !old Ammonium stuff
!                                , EQUILIB_MARS      = .true.  & !MARS
!                                , EQUILIB_EQSAM     = .false.     !EQSAM
                                
!.. Number of aerosol sizes (1-fine, 2-coarse, 3-'giant' for sea salt )
!    integer, public, parameter :: NSIZE = 5
!           !   FINE_PM = 1, COAR_NO3 = 2, COAR_SS = 3, COAR DUST = 4,pollen = 5    


contains

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine AerosolEquilib(debug_flag)
    logical, intent(in) :: debug_flag
!    integer, intent(in)  :: i, j
    logical, save :: my_first_call=.true.
    
    if( my_first_call .and. MasterProc ) then
       write(*,*) 'AerosolEquilib: chosen: ',AERO%EQUILIB 
       write(*,*) 'AerosolEquilib water: chosen: ',AERO%EQUILIB_WATER
    end if
    select case ( AERO%EQUILIB )
      case ( 'EMEP' )
        call ammonium()
      case ( 'MARS' , 'MARS_2900', 'GEOSCHEM')
        call emep2MARS(debug_flag)
      case ( 'EQSAM' )
        call emep2EQSAM(debug_flag)
      case ( 'ISORROPIA' )
        call StopAll('Isorropia problems found. Removed for now')
        !call emep2Isorropia(debug_flag)
      case default
        if( my_first_call .and. MasterProc ) then
          write(*,*) 'WARNING: AerosolEquilib: nothing chosen:'
        end if
    end select
    my_first_call = .false.

  end subroutine AerosolEquilib

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 ! Adapted from List 10, p130, Isoropia manual

 !FUTURE subroutine emep2isorropia(debug_flag)
!FUTURE   logical, intent(in) :: debug_flag
!FUTURE
!FUTURE   real, dimension(8) :: wi = 0.0, wt
!FUTURE   real, dimension(3) :: gas
!FUTURE   real, dimension(15) :: aerliq
!FUTURE   real, dimension(19) :: aersld
!FUTURE   real, parameter, dimension(2) :: CNTRL =  (/ 0, 0 /)
!FUTURE   real, dimension(9) :: other
!FUTURE   !real :: rhi, tempi
!FUTURE   character(len=15) :: scase
!FUTURE
!FUTURE   !EMEP
!FUTURE   real, parameter :: Ncm3_to_molesm3 = 1.0e6/AVOG    ! #/cm3 to moles/m3
!FUTURE   real, parameter :: molesm3_to_Ncm3 = 1.0/Ncm3_to_molesm3
!FUTURE   integer :: k
!FUTURE!BUG   real :: tmpno3
!FUTURE
!FUTURE   ! WI(1)  = max(FLOOR2, xn_2d(Na,k))  / species(Na)%molwt  * Ncm3_to_molesm3
!FUTURE   ! 5=Cl, 6=Ca, 7=K, 8=Mg
!FUTURE
!FUTURE   do k = KMAX_MID, KMAX_MID  ! TESTING KCHEMTOP, KMAX_MID
!FUTURE
!FUTURE     WI(1)  = 0.0 !FINE sum( xn_2d(SS_GROUP,k) ) * Ncm3_to_molesm3
!FUTURE     WI(2)  = xn_2d(SO4,k)             * Ncm3_to_molesm3
!FUTURE     WI(3)  = sum( xn_2d(RDN_GROUP,k) ) * Ncm3_to_molesm3  !NH3, NH4
!FUTURE     !FINE WI(4)  = ( xn_2d(NO3_F,k) + xn_2d(NO3_C,k) + xn_2d(HNO3,k) )&
!FUTURE     WI(4)  = ( xn_2d(NO3_F,k) + xn_2d(HNO3,k) )&
!FUTURE                * Ncm3_to_molesm3
!FUTURE     WI(5)  =0.0 !FINE  WI(1)  ! Cl only from sea-salt. Needs consideration!
!FUTURE
!FUTURE     call isoropia ( wi, rh(k), temp(k), CNTRL,&
!FUTURE                     wt, gas, aerliq, aersld, scase, other)
!FUTURE
!FUTURE    ! gas outputs are in moles/m3(air)
!FUTURE
!FUTURE     xn_2d(NH3,k)  = gas(1) * molesm3_to_Ncm3
!FUTURE     xn_2d(HNO3,k) = gas(2) * molesm3_to_Ncm3
!FUTURE     !xn_2d(HCl,k) = gas(3) * molesm3_to_Ncm3
!FUTURE
!FUTURE    ! aerosol outputs are in moles/m3(air)
!FUTURE    ! 1=H+, 2=Na+, 3=NH4+, 4=Cl-, 5=SO42-, 6=HSO4-, 7=NO3-, 8=Ca2+
!FUTURE    ! 9=K+, 10=Mg2+
!FUTURE     !xn_2d(NH4_F,k) = MOLAL(3)
!FUTURE
!FUTURE    ! Just use those needed:
!FUTURE    ! QUERY: Is NaNO3 always solid? Ans = No!
!FUTURE
!FUTURE      !xn_2d(NO3_c,k ) = aeroHCl * molesm3_to_Ncm3 ! assume all HCl from NaNO3 formation?
!FUTURE      !FINE xn_2d(NO3_f,k ) = tmpno3 - xn_2d(NO3_c,k ) - xn_2d(HNO3,k)
!FUTURE      xn_2d(NO3_f,k ) = tmpno3 - xn_2d(HNO3,k)
!FUTURE
!FUTURE     if( debug_flag ) then 
!FUTURE       write(*, "(a,2f8.3,99g12.3)") "ISORROPIA ", rh(k), temp(k), gas
!FUTURE     end if
!FUTURE     !call StopAll("ISOR")
!FUTURE     
!FUTURE   end do
!FUTURE end subroutine emep2isorropia
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine emep2MARS(debug_flag)

 !..................................................................
 ! Pretty old F. Binkowski code from EPA CMAQ-Models3
 ! JGR, 108, D6, 4183
 !..................................................................

 logical, intent(in) :: debug_flag 
 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  
 real, parameter ::    FLOOR2 = 1.0E-9         ! minimum concentration  

 !.. local
  real    :: so4in, no3in, nh4in, hno3in, nh3in,   &
             aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out,   &
             coef
  integer :: k, errmark
 !-----------------------------------

   coef = 1.e12 / AVOG

   do k = KCHEMTOP, KMAX_MID
  

!//.... molec/cm3 -> ug/m3
! Use FLOOR2 = 1.0e-8 molec/cm3 for input. Too many problems
      so4in  = max(FLOOR2, xn_2d(SO4,k)) * species(SO4)%molwt  *coef
      hno3in = max(FLOOR2, xn_2d(HNO3,k))* species(HNO3)%molwt *coef 
      nh3in  = max(FLOOR2, xn_2d(NH3,k)) * species(NH3)%molwt  *coef
      no3in  = max(FLOOR2, xn_2d(NO3_f,k)) * species(NO3_f)%molwt  *coef
      nh4in  = max(FLOOR2, xn_2d(NH4_f,k)) * species(NH4_f)%molwt  *coef

 !--------------------------------------------------------------------------                
      if(AERO%EQUILIB=='MARS')then 
         call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='MARS_2900')then!svn version 2908, for testing if there are significant differences. will be deleted
         call rpmares_2900 (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='GEOSCHEM')then
         call DO_RPMARES_new (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag)
      end if

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

   end do  ! K-levels

 end subroutine emep2MARS

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine emep2EQSAM(debug_flag)
 logical, intent(in) :: debug_flag 
! integer, intent(in)  :: i, j

    !..................................................................
    ! EQSAM4clim - Equlibrium Simplified Aerosol Model by Swen Metzger
    !             version v10 is implemented here
    ! Metzger, S., B. Steil, M. Abdelkader, K. Klingmüller, L. Xu, 
    ! J.E. Penner, C. Fountoukis, A. Nenes, and J. Lelieveld, 2016; 
    ! Aerosol Water Parameterization: A single parameter framework; 
    ! Atmos. Chem. Phys., 16, 7213–7237, doi:10.5194/acp-16-7213-2016; 
    ! https://www.atmos-chem-phys.net/16/7213/2016/acp-16-7213-2016.html
    !..................................................................


 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  


 !.. local

  real    :: so4in(KCHEMTOP:KMAX_MID),   &
             no3in(KCHEMTOP:KMAX_MID),   &
             nh4in(KCHEMTOP:KMAX_MID),   &
             hno3in(KCHEMTOP:KMAX_MID),  &
             nh3in(KCHEMTOP:KMAX_MID),   &
             aH2Oin(KCHEMTOP:KMAX_MID),  &
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
  end if

!//.... molec/cm3 -> micromoles/m**3
    so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4,KCHEMTOP:KMAX_MID)  *1.e12/AVOG
    hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3,KCHEMTOP:KMAX_MID) *1.e12/AVOG
    nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3,KCHEMTOP:KMAX_MID)  *1.e12/AVOG 
    no3in(KCHEMTOP:KMAX_MID)  = xn_2d(NO3_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(NH4_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG
!    aH2Oin(KCHEMTOP:KMAX_MID) = H2O_eqsam(i,j,KCHEMTOP:KMAX_MID)
    NAin(KCHEMTOP:KMAX_MID)   = xn_2d(SEASALT_F,KCHEMTOP:KMAX_MID)*0.306e12/AVOG
    CLin(KCHEMTOP:KMAX_MID)   = xn_2d(SEASALT_F,KCHEMTOP:KMAX_MID)*0.55e12/AVOG

 !--------------------------------------------------------------------------                
  
!    call eqsam_v03d (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh,temp,pp, &
!                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,             &
!                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout)
    
    call EQSAM4clim (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh, temp,  & !aH2Oin, &
                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,            &
                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout, KCHEMTOP, KMAX_MID)
  
 !--------------------------------------------------------------------------

!//.... micromoles/m**3  -> molec/cm3 
!      xn_2d(NO3,KCHEMTOP:KMAX_MID)  = FLOOR !different for ACID/OZONE

      xn_2d(HNO3,KCHEMTOP:KMAX_MID)  = max(FLOOR,gNO3out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
      xn_2d(NH3,KCHEMTOP:KMAX_MID)   = max(FLOOR,gNH3out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
      xn_2d(NO3_f,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNO3out(KCHEMTOP:KMAX_MID)*AVOG*1.e-12 ) 
      xn_2d(NH4_f,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNH4out(KCHEMTOP:KMAX_MID)*AVOG*1.e-12 )
      xn_2d(SO4,KCHEMTOP:KMAX_MID)   = max(FLOOR,aSO4out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
!      H2O_eqsam(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )

 if ( debug_flag ) then ! Selected debug cell
    write(*,*)'After EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(NO3_f,20),xn_2d(NH4_f,20)
  end if

 end subroutine emep2EQSAM


 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


 !water 

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine Aero_water_rh50(i,j, debug_flag)

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
 logical, intent(in)  :: debug_flag
 !.. local
! integer, parameter    ::  KCHEMTOP = KMAX_MID

!  AerosolCalls.f90(319): error #6592: This symbol must be a defined parameter, an enumerator, or an argument of 
!  an inquiry function that evaluates to a compile-time constant.   [KMAX_MID]
!  integer, parameter    ::  KCHEMTOP = KMAX_MID
!--------------------------------------^

  real    :: so4in(KMAX_MID:KMAX_MID),   &
             no3in(KMAX_MID:KMAX_MID),   &
             nh4in(KMAX_MID:KMAX_MID),   &
             hno3in(KMAX_MID:KMAX_MID),  &
             nh3in(KMAX_MID:KMAX_MID),   &
             aH2Oin(KMAX_MID:KMAX_MID),  &
             NAin(KMAX_MID:KMAX_MID)  ,  &
             CLin(KMAX_MID:KMAX_MID) ,   &
!-- output
             aSO4out(KMAX_MID:KMAX_MID), &
             aNO3out(KMAX_MID:KMAX_MID), &
             aH2Oout(KMAX_MID:KMAX_MID), &
             aNH4out(KMAX_MID:KMAX_MID), & 
             gNH3out(KMAX_MID:KMAX_MID), &
             gNO3out(KMAX_MID:KMAX_MID), &
             aNAout(KMAX_MID:KMAX_MID),  &
             aCLout(KMAX_MID:KMAX_MID),  &
             gCLout(KMAX_MID:KMAX_MID),  &
             gSO4out(KMAX_MID:KMAX_MID), &

             rlhum(KMAX_MID:KMAX_MID),tmpr(KMAX_MID:KMAX_MID)

 !-----------------------------------


  if ( debug_flag ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(NO3_f,20),xn_2d(NH4_f,20)
  end if


!//.... molec/cm3 -> micromoles/m**3
      so4in(KMAX_MID:KMAX_MID)  = xn_2d(SO4,KMAX_MID:KMAX_MID)*1.e12/AVOG
      hno3in(KMAX_MID:KMAX_MID) = xn_2d(HNO3,KMAX_MID:KMAX_MID)*1.e12/AVOG
      nh3in(KMAX_MID:KMAX_MID)  = xn_2d(NH3,KMAX_MID:KMAX_MID)*1.e12/AVOG 
      no3in(KMAX_MID:KMAX_MID)  = xn_2d(NO3_f,KMAX_MID:KMAX_MID)*1.e12/AVOG
      nh4in(KMAX_MID:KMAX_MID)  = xn_2d(NH4_f,KMAX_MID:KMAX_MID)*1.e12/AVOG
      NAin(KMAX_MID:KMAX_MID)   = xn_2d(SEASALT_F,KMAX_MID:KMAX_MID)*0.306e12/AVOG
      CLin(KMAX_MID:KMAX_MID)   = xn_2d(SEASALT_F,KMAX_MID:KMAX_MID)*0.55e12/AVOG

!.. Rh = 50% and T=20C
      rlhum(KMAX_MID:KMAX_MID) = 0.5
      tmpr (KMAX_MID:KMAX_MID) = 293.15


 !--------------------------------------------------------------------------                
   
    call EQSAM4clim (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rlhum, tmpr,  & 
                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,               &
                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout, KMAX_MID, KMAX_MID)     

 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 

      PM25_water_rh50 (i,j)             = max(0., aH2Oout(KMAX_MID) )


 end subroutine  Aero_water_rh50
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 subroutine Aero_water(i,j, debug_flag)


    !..................................................................
    ! EQSAM4clim calculates PM water at ambient conditions
    !..................................................................

 integer, intent(in)  :: i, j
 logical, intent(in)  :: debug_flag

 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  


 !.. local

  real    :: so4in(KCHEMTOP:KMAX_MID),   &
             no3in(KCHEMTOP:KMAX_MID),   &
             nh4in(KCHEMTOP:KMAX_MID),   &
             hno3in(KCHEMTOP:KMAX_MID),  &
             nh3in(KCHEMTOP:KMAX_MID),   &
             aH2Oin(KCHEMTOP:KMAX_MID),  &
             NAin(KCHEMTOP:KMAX_MID),    &
             CLin(KCHEMTOP:KMAX_MID),    &
!.. output
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

 !-----------------------------------


  if ( debug_flag  ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(NO3_f,20),xn_2d(NH4_f,20)
  end if

!//.... molec/cm3 -> micromoles/m**3
    so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4,KCHEMTOP:KMAX_MID)  *1.e12/AVOG
    hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3,KCHEMTOP:KMAX_MID) *1.e12/AVOG
    nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3,KCHEMTOP:KMAX_MID)  *1.e12/AVOG 
    no3in(KCHEMTOP:KMAX_MID)  = xn_2d(NO3_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(NH4_f,KCHEMTOP:KMAX_MID)*1.e12/AVOG
    NAin(KCHEMTOP:KMAX_MID)   = xn_2d(SEASALT_F,KCHEMTOP:KMAX_MID)*0.306e12/AVOG
    CLin(KCHEMTOP:KMAX_MID)   = xn_2d(SEASALT_F,KCHEMTOP:KMAX_MID)*0.55e12/AVOG

    rlhum(KCHEMTOP:KMAX_MID) = rh(:)
    tmpr(KCHEMTOP:KMAX_MID)  = temp(:)

 !--------------------------------------------------------------------------                
 
    call EQSAM4clim (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rlhum, tmpr,  & !aH2Oin, &
                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,            &
                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout, KCHEMTOP, KMAX_MID)     

 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 

      PM25_water(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )

 end subroutine Aero_water


 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine Aero_water_MARS(i,j, debug_flag)

 !..................................................................
 ! Pretty old F. Binkowski code from EPA CMAQ-Models3
 ! JGR, 108, D6, 4183
 !..................................................................

 integer, intent(in)  :: i, j
 logical, intent(in)  :: debug_flag

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
      if(AERO%EQUILIB=='MARS')then 
         call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='MARS_2900')then!svn version 2908, for testing if there are significant differences. will be deleted
         call rpmares_2900 (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='GEOSCHEM')then
         call DO_RPMARES_new (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      end if
      !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 
      PM25_water(i,j,k) = max (0., aH2Oout )

    end do  ! k-loop

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
     if(AERO%EQUILIB=='MARS')then 
         call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='MARS_2900')then!svn version 2908, for testing if there are significant differences. will be deleted
         call rpmares_2900 (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='GEOSCHEM')then
         call DO_RPMARES_new (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      end if
  !--------------------------------------------------------------------------

      PM25_water_rh50 (i,j) = max (0., aH2Oout )


 end subroutine  Aero_water_MARS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module AerosolCalls
