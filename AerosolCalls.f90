! <AerosolCalls.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.0>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2023 met.no
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
 use ChemSpecs_mod,         only: species
 use Chemfields_mod,        only: PM25_water, PM25_water_rh50, & !H2O_eqsam, & !PMwater 
                                   cfac, pH
 use Config_module,         only: KMAX_MID, KCHEMTOP, MasterProc, USES,&
                                  SO4_ix, HNO3_ix, NO3_f_ix, NH3_ix, NH4_f_ix, OM_ix, &
                                  SSf_ix, SSc_ix, Dustwbf_ix, DustSahf_ix
 use Debug_module,          only: DEBUG   ! -> DEBUG%EQUIB
 use EQSAM4clim_ml,        only :  EQSAM4clim
! use EQSAM_v03d_mod,        only: eqsam_v03d
 use MARS_mod,              only: rpmares, rpmares_2900, DO_RPMARES_new
 use PhysicalConstants_mod, only: AVOG
 use SmallUtils_mod,        only: find_index
 use ZchemData_mod,         only: xn_2d, temp, rh, pp
 use Par_mod,               only: me
 implicit none
 private

 !/-- public           !!  true if wanted

 public :: AerosolEquilib
 public :: emep2MARS, emep2EQSAM,  Aero_Water_MARS ! , Aero_Water  , Aero_Water_rh50
 private :: emep2isorropia
                    
!    logical, public, parameter :: AERO_DYNAMICS     = .false.  &  
!                                , EQUILIB_EMEP      = .false.  & !old Ammonium stuff
!                                , EQUILIB_MARS      = .true.  & !MARS
!                                , EQUILIB_EQSAM     = .false.     !EQSAM
                                
!.. Number of aerosol sizes (1-fine, 2-coarse, 3-'giant' for sea salt )
!    integer, public, parameter :: NSIZE = 5
!           !   FINE_PM = 1, COAR_NO3 = 2, COAR_SS = 3, COAR DUST = 4,pollen = 5    


  integer, private, save :: iSeaSalt,cSeaSalt,cNatDust  ! zero if no seasalt compounds
  real, parameter :: MWCL = 35.453, MWNA = 22.9897,  & ! MW Chloride, Sodium
                     MWCA = 40.078, MWMG = 24.305, MWK = 39.098, & ! MW Calcium, Magnegium, Potassium
                     MWH2O = 18.0153                   ! MW Water

contains

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine AerosolEquilib(i,j,debug_flag)
    logical, intent(in) :: debug_flag
    integer, intent(in)  :: i, j
    logical, save :: my_first_call=.true.
    character(len=*),parameter:: dtxt='AeroEqui:'
    
    if( my_first_call ) then
      iSeaSalt = find_index('SeaSalt_f',species(:)%name )
      call CheckStop(USES%SEASALT.and.iSeaSalt<1,dtxt//"iSeaSalt neg")
      if(  MasterProc ) then
        write(*,*) 'AerosolEquilib: chosen: ',AERO%EQUILIB 
        write(*,*) 'AerosolEquilib water: chosen: ',AERO%EQUILIB_WATER
        write(*,*) 'AerosolEquilib seasalt index: ',iSeaSalt
      end if
    end if
    select case ( AERO%EQUILIB )
      case ( 'EMEP' )
        call ammonium()
      case ( 'MARS' , 'MARS_2900', 'GEOSCHEM')
        call emep2MARS(debug_flag)
      case ( 'EQSAM' )
        call emep2EQSAM(i, j, debug_flag)
      case ( 'ISORROPIA' )
        !NOV22 call StopAll('Isorropia problems found. Removed for now')
        call emep2Isorropia(i,j,debug_flag)
      case default
        if( my_first_call .and. MasterProc ) then
          write(*,*) 'WARNING! AerosolEquilib, nothing valid chosen: '//AERO%EQUILIB
        end if
    end select
    my_first_call = .false.

  end subroutine AerosolEquilib

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 ! Adapted from List 10, p130, Isoropia manual

  subroutine emep2isorropia(i,j,debug_flag)
    integer, intent(in) :: i, j
    logical, intent(in) :: debug_flag
  
    real, dimension(8)  :: wi = 0.0, wt = 0.0
    real, dimension(3)  :: wo = 0.0, gas = 0.0
    real, dimension(39) :: aerliq
    real, dimension(19) :: aersld
  
    real, parameter, dimension(1) :: CNTRL =  0 ! 0 - Forward problem. WI = GAS + AEROSOL.
    real, dimension(9) :: other
    character(len=15), save  :: scase
  
    ! EMEP declarations
    real, parameter :: Ncm3_to_molesm3 = 1.0e6/AVOG    ! #/cm3 to moles/m3
    real, parameter :: molesm3_to_Ncm3 = 1.0/Ncm3_to_molesm3
    integer :: k, n
    real, parameter :: CONMIN = 1.0e-30   ! Concentration lower limit [mole/m3]
    real, parameter :: MWCL   = 35.453,   MWNA   = 22.9897 ! MW Chloride, Sodium
    real, parameter :: MWCA   = 40.078,   MWK    = 39.0973 ! MW Calcium, Potassium
    real, parameter :: MWH2O  = 18.0153,  MWMG   = 24.305  ! MW Water, Magnesium
  
    real :: tmpno3, tmpnh3, tmpnhx, tmphno3
    logical, save :: first_isor = .true.
  
  !  INPUT:
  !  1. [WI]
  !     DOUBLE PRECISION array of length [8].
  !     Concentrations, expressed in moles/m3. Depending on the type of
  !     problem solved (specified in CNTRL(1)), WI contains either
  !     GAS+AEROSOL or AEROSOL only concentratios.
  !     WI(1) - sodium    ! WI(2) - sulfate   ! WI(3) - ammonium    ! WI(4) - nitrate
  !     WI(5) - chloride  ! WI(6) - calcium   ! WI(7) - potassium   ! WI(8) - magnesium
  !
  !  2. [WO]
  !     DOUBLE PRECISION array of length [3].
  !     Organic aerosol present in the inorganic phase. Used to calculate water
  !     uptake that interacts with the inorganic species.
  !     WI(1) - organic mass concentration in air, kg/m3 air
  !     WI(2) - organic hygroscopicity parameter kappa (as in Petters and Kreidenweis, 2007)
  !     WI(3) - organic aerosol density, kg/m3
  !
  !  OUTPUT:
  !  1. [WT]
  !     DOUBLE PRECISION array of length [8].
  !     Total concentrations (GAS+AEROSOL) of inorganic species, expressed in moles/m3.
  !     If the foreward probelm is solved (CNTRL(1)=0), array WT is
  !     identical to array WI.
  !     WT(1) - total sodium     ! WT(2) - total sulfate     ! WT(3) - total ammonium
  !     WT(4) - total nitrate    ! WT(5) - total chloride    ! WT(6) - total calcium
  !     WT(7) - total potassium  ! WT(8) - total magnesium
  !
  !  2. [GAS]
  !     DOUBLE PRECISION array of length [03].
  !     Gaseous species concentrations, expressed in moles/m3.
  !     GAS(1) - NH3  ! GAS(2) - HNO3  ! GAS(3) - HCl
  !
  !  3. [AERLIQ]
  !     DOUBLE PRECISION array of length [39].
  !     Liquid aerosol species concentrations, expressed in moles/m3.
  !     AERLIQ(01) - H+(aq)
  !     AERLIQ(02) - Na+(aq)
  !     AERLIQ(03) - NH4+(aq)
  !     AERLIQ(04) - Cl-(aq)
  !     AERLIQ(05) - SO4--(aq)
  !     AERLIQ(06) - HSO4-(aq)
  !     AERLIQ(07) - NO3-(aq)
  !     AERLIQ(08) - Total aerosol water
  !     + many more deliquesced and undissociated forms
  !
  !  4. [AERSLD]
  !     DOUBLE PRECISION array of length [19].
  !     Solid aerosol species concentrations, expressed in moles/m3.
  !     AERSLD(01) - NaNO3(s)
  !     AERSLD(02) - NH4NO3(s)
  !     AERSLD(03) - NaCl(s)
  !     AERSLD(04) - NH4Cl(s)
  !     AERSLD(05) - Na2SO4(s)
  !     AERSLD(06) - (NH4)2SO4(s)
  !     AERSLD(07) - NaHSO4(s)
  !     AERSLD(08) - NH4HSO4(s)
  !     AERSLD(09) - (NH4)4H(SO4)2(s)
  !     AERSLD(10) - CaSO4(s)
  !     AERSLD(11) - Ca(NO3)2(s)
  !     AERSLD(12) - CaCl2(s)
  !     AERSLD(13) - K2SO4(s)
  !     AERSLD(14) - KHSO4(s)
  !     AERSLD(15) - KNO3(s)
  !     AERSLD(16) - KCl(s)
  !     AERSLD(17) - MgSO4(s)
  !     AERSLD(18) - Mg(NO3)2(s)
  !     AERSLD(19) - MgCl2(s)
  
     
    do n = 1,1 ! (fine, coarse) -- optional. Only fine for now; coarse needs tinkering.
      
      do k = KCHEMTOP,KMAX_MID 
  
        ! isorropia only for when T > 250 K and P > 200 hPa (CMAQ and GEOS-Chem; Shannon Capps discussion)
        if (pp(k) > 20000.0 .and. temp(k) > 250.0) then 
  
          if ( AERO%INTERNALMIXED .and. SSf_ix > 0 ) then
            ! if internally mixed is assumed, Na and Cl as mass-fraction of seasalt aerosol.
            WI(1) = max(0., xn_2d(SSf_ix ,k) * Ncm3_to_molesm3 * species(SSf_ix)%molwt * 0.30 / MWNA )
            WI(5) = max(0., xn_2d(SSf_ix ,k) * Ncm3_to_molesm3 * species(SSf_ix)%molwt * 0.55 / MWCL )
          else
            WI(1)  = 0.0
            WI(5)  = 0.0
          endif
  
          WI(2)  = MAX(  xn_2d(SO4_ix  ,k) * Ncm3_to_molesm3, CONMIN )
          WI(3)  = MAX(( xn_2d(NH3_ix  ,k) + xn_2d(NH4_f_ix,k) ) * Ncm3_to_molesm3, CONMIN ) ! NH3, NH4
          WI(4)  = MAX(( xn_2d(NO3_f_ix,k) + xn_2d(HNO3_ix, k) ) * Ncm3_to_molesm3, CONMIN )
  
          if ( AERO%INTERNALMIXED .and. AERO%CATIONS .and. SSf_ix > 0) then
            ! mass percentages of sea salt for below species taken from GEOS-Chem
            WI(6) = max(0., xn_2d(SSf_ix ,k) * Ncm3_to_molesm3 * species(SSf_ix)%molwt * 0.0116 / MWCA )
            WI(7) = max(0., xn_2d(SSf_ix ,k) * Ncm3_to_molesm3 * species(SSf_ix)%molwt * 0.0110 / MWK  )
            WI(8) = max(0., xn_2d(SSf_ix ,k) * Ncm3_to_molesm3 * species(SSf_ix)%molwt * 0.0369 / MWMG )
          else
            WI(6) = 0.0
            WI(7) = 0.0
            WI(8) = 0.0
          end if
  
          ! contributions from dust using values for 'other' crustal species from Karydis et al. 2015 Table 2
          if (AERO%CATIONS .and. Dustwbf_ix > 0) then 
            WI(1) = WI(1) + max(0., xn_2d(Dustwbf_ix ,k) * Ncm3_to_molesm3 * species(Dustwbf_ix)%molwt * 0.012 / MWNA )
            WI(6) = WI(6) + max(0., xn_2d(Dustwbf_ix ,k) * Ncm3_to_molesm3 * species(Dustwbf_ix)%molwt * 0.024 / MWCA )
            WI(7) = WI(7) + max(0., xn_2d(Dustwbf_ix ,k) * Ncm3_to_molesm3 * species(Dustwbf_ix)%molwt * 0.015 / MWK  )
            WI(8) = WI(8) + max(0., xn_2d(Dustwbf_ix ,k) * Ncm3_to_molesm3 * species(Dustwbf_ix)%molwt * 0.009 / MWMG )
          endif
  
          if (AERO%CATIONS .and. DustSahf_ix > 0) then 
            WI(1) = WI(1) + max(0., xn_2d(DustSahf_ix ,k) * Ncm3_to_molesm3 * species(DustSahf_ix)%molwt * 0.012 / MWNA )
            WI(6) = WI(6) + max(0., xn_2d(DustSahf_ix ,k) * Ncm3_to_molesm3 * species(DustSahf_ix)%molwt * 0.024 / MWCA )
            WI(7) = WI(7) + max(0., xn_2d(DustSahf_ix ,k) * Ncm3_to_molesm3 * species(DustSahf_ix)%molwt * 0.015 / MWK  )
            WI(8) = WI(8) + max(0., xn_2d(DustSahf_ix ,k) * Ncm3_to_molesm3 * species(DustSahf_ix)%molwt * 0.009 / MWMG )
          endif
  
          ! for NOV22 testing:
          tmpnh3  = xn_2d(NH3_ix,k)
          tmpnhx  = tmpnh3 + xn_2d(NH4_f_ix,k)
          tmphno3 = xn_2d(HNO3_ix,k)
          tmpno3  = tmphno3 + xn_2d(NO3_f_ix,k)
  
          ! OM25_p is the sum of the particle-phase OM25, and currently has MW 1 for simplicity
          wo(1) = xn_2d(OM_ix,k) * Ncm3_to_molesm3 * species(OM_ix)%molwt * 1e-3 ! kg/m3 organic matter
          wo(2) = 0.15 ! kappa organic aerosol; standard value from stand-alone input file
          wo(3) = 1400 ! aerosol density kg/m3; Based on observations (Kakavas, 2023) & florou et al., 2014
  
          call isoropia ( wi, wo, rh(k),  temp(k), CNTRL, &       
                          wt, gas, aerliq, aersld, scase, other )  
  
          ! pH = -log10([H+]/M) where M = mol dm-3 in the solution. mol/m3 h2o to kg/m3 as 1e-3 * MWH20.
          ! 1 liter water ~ 1 kg, such that pH = -log10( [H+] / ([H2O] * MWH2O * 1e-3) )
          if ( k == KMAX_MID) &
            pH(i,j) = -SAFELOG10( aerliq(1) / (aerliq(8) * MWH2O) ) - 3
  
          ! gas outputs are in moles/m3(air)
          xn_2d(NH3_ix ,k) = max( gas(1), CONMIN ) * molesm3_to_Ncm3
          xn_2d(HNO3_ix,k) = max( gas(2), CONMIN ) * molesm3_to_Ncm3
          !xn_2d(HCl,k) = gas(3) * molesm3_to_Ncm3
  
          ! aerosol outputs are in moles/m3(air)
          xn_2d(NH4_f_ix,k) = max( wt(3) - gas(1), CONMIN ) * molesm3_to_Ncm3
          
          ! aerosol water (ug/m**3) -- 18.01528 MW H2O
          PM25_water(i,j,k) = max( aerliq(8) * MWH2O * 1e6, CONMIN )
  
          ! QUERY: Is NaNO3 always solid? Ans = No!
  
          !xn_2d(NO3_c,k ) = aeroHCl * molesm3_to_Ncm3 ! assume all HCl from NaNO3 formation?
          !FINE xn_2d(NO3_f,k ) = tmpno3 - xn_2d(NO3_c,k ) - xn_2d(HNO3,k)
  
          !tmpno3 = wt(4) * molesm3_to_Ncm3  ! NOV22  wt4=nitrate
          !xn_2d(NO3_f_ix,k ) = tmpno3 - xn_2d(HNO3_ix,k)
          !NOV22 - get some v.small neg., so use max below. Test properly later.
          xn_2d(NO3_f_ix,k ) = max( wt(4) - gas(2), CONMIN )  * molesm3_to_Ncm3 
  
    !      H2O_eqsam(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )
          
          if (k == KMAX_MID ) then
            ! at Rh=50% and T=20C (comparison with gravimentric PM)
            call isoropia ( wi, wo, 0.5, 293.15, CNTRL, &       
                            wt, gas, aerliq, aersld, scase, other) 
  
            ! aerosol water (ug/m**3) -- 18.01528 MW H2O
            PM25_water_rh50(i,j) =  max( 0., aerliq(8) * MWH2O * 1e6 )
          end if 
  
          !if ( xn_2d(NO3_f_ix,k )  < 0.0 .or.  xn_2d(NH4_f_ix,k) < 0.0 ) then
          !   if ( xn_2d(NO3_f_ix,k ) 
          !   print "(a,99e12.3)", 'NEGNO3 pre',tmpnh3, tmpnhx, tmphno3, tmpno3
          !   print "(a,99e12.3)", 'NEGNO3 xn', xn_2d(NH3_ix,k ), xn_2d(NH4_f_ix,k)+xn_2d(NH3_ix,k ), xn_2d(HNO3_ix,k ), xn_2d(NO3_f_ix,k ), xn_2d(SO4_ix,k )
          !   print "(a,99e12.3)", 'NEGNO3 wi', wi(1:4) ! 1 - total sodium 2 - total sulfate 3 - total ammonium 4 - total nitrate
          !   print "(a,99e12.3)", 'NEGNO3 wt', wt(1:4) ! 1 - total sodium 2 - total sulfate 3 - total ammonium 4 - total nitrate
          !   print "(a,99e12.3)", 'NEGNO3 gas', gas !     GAS(1) - NH3 !     GAS(2) - HNO3 !     GAS(3) - HCl
          !   print "(a,99e12.3)", 'NEGNO3 diffs',  xn_2d(NO3_f_ix,k ),  xn_2d(NH4_f_ix,k)
          !   call StopAll('NEGNO3')
          !end if
  
          if( debug_flag ) then 
            write(*, "(a,2f8.3,99g12.3)") "ISORROPIA ", rh(k), temp(k), gas
          end if
          !call StopAll("ISOR")
  
        endif ! > 200 hPa and > 250 K
      end do ! k = KCHEMTOP, KMAX_MID
    end do ! n = 1,2 (fine, coarse)
  end subroutine emep2isorropia
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
  call CheckStop( SO4_ix<1, "emep2MARS:SO4 not defined" )
  call CheckStop( NH4_f_ix<1, "emep2MARS:NH4_f not defined" )
  call CheckStop( HNO3_ix<1, "emep2MARS:HNO3 not defined" )
  call CheckStop( NO3_f_ix<1, "emep2MARS:NO3_f not defined" )
  call CheckStop( NH3_ix<1, "emep2MARS: NH3 not defined" )

   coef = 1.e12 / AVOG

   do k = KCHEMTOP, KMAX_MID
  

!//.... molec/cm3 -> ug/m3
! Use FLOOR2 = 1.0e-8 molec/cm3 for input. Too many problems
      so4in  = max(FLOOR2, xn_2d(SO4_ix,k)) * species(SO4_ix)%molwt  *coef
      hno3in = max(FLOOR2, xn_2d(HNO3_ix,k))* species(HNO3_ix)%molwt *coef 
      nh3in  = max(FLOOR2, xn_2d(NH3_ix,k)) * species(NH3_ix)%molwt  *coef
      no3in  = max(FLOOR2, xn_2d(NO3_f_ix,k)) * species(NO3_f_ix)%molwt  *coef
      nh4in  = max(FLOOR2, xn_2d(NH4_f_ix,k)) * species(NH4_f_ix)%molwt  *coef

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

      xn_2d(HNO3_ix,k)  = max (FLOOR, gNO3out / (species(HNO3_ix)%molwt *coef) )
      xn_2d(NH3_ix,k)   = max (FLOOR, gNH3out / (species(NH3_ix)%molwt  *coef) )
      xn_2d(NO3_f_ix,k)  = max (FLOOR, aNO3out / (species(NO3_f_ix)%molwt  *coef) )
      xn_2d(NH4_f_ix,k)  = max (FLOOR, aNH4out / (species(NH4_f_ix)%molwt  *coef) )

   end do  ! K-levels

 end subroutine emep2MARS

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 subroutine emep2EQSAM(i, j, debug_flag)

    !..................................................................
    ! EQSAM4clim - Equlibrium Simplified Aerosol Model by Swen Metzger
    !             version v10 is implemented here
    ! Metzger, S., B. Steil, M. Abdelkader, K. Klingmüller, L. Xu, 
    ! J.E. Penner, C. Fountoukis, A. Nenes, and J. Lelieveld, 2016; 
    ! Aerosol Water Parameterization: A single parameter framework; 
    ! Atmos. Chem. Phys., 16, 7213–7237, doi:10.5194/acp-16-7213-2016; 
    ! https://www.atmos-chem-phys.net/16/7213/2016/acp-16-7213-2016.html
    ! EQSAM4clim calculates PM water: 1. at ambient conditions and 
    !          2. at Rh=50% and T=20C (comparison with gravimentric PM)
    !..................................................................

 logical, intent(in) :: debug_flag 
 integer, intent(in) :: i, j
 real, parameter :: FLOOR = 1.0E-30         ! minimum concentration  
 real, parameter :: MWCL   = 35.453,   MWNA   = 22.9897 ! MW Chloride, Sodium
 real, parameter :: MWH2O  = 18.0153

 !.. local

  real    :: so4in(KCHEMTOP:KMAX_MID), so4ins(KMAX_MID:KMAX_MID),  &
             no3in(KCHEMTOP:KMAX_MID), no3ins(KMAX_MID:KMAX_MID),  &
             nh4in(KCHEMTOP:KMAX_MID), nh4ins(KMAX_MID:KMAX_MID),  &
             hno3in(KCHEMTOP:KMAX_MID),hno3ins(KMAX_MID:KMAX_MID), &
             nh3in(KCHEMTOP:KMAX_MID), nh3ins(KMAX_MID:KMAX_MID),  &
             aH2Oin(KCHEMTOP:KMAX_MID),aH2Oins(KMAX_MID:KMAX_MID), &
             NAin(KCHEMTOP:KMAX_MID),  Nains(KMAX_MID:KMAX_MID),   &
             CAin(KCHEMTOP:KMAX_MID),  CAins(KMAX_MID:KMAX_MID),   &
             MGin(KCHEMTOP:KMAX_MID),  MGins(KMAX_MID:KMAX_MID),   &
             Kin (KCHEMTOP:KMAX_MID),  Kins (KMAX_MID:KMAX_MID),   &
             CLin(KCHEMTOP:KMAX_MID),  CLins(KMAX_MID:KMAX_MID),   &
!.. output
             aSO4out(KCHEMTOP:KMAX_MID), aSO4outs(KMAX_MID:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), aNO3outs(KMAX_MID:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID), aH2Oouts(KMAX_MID:KMAX_MID), &
             apHout (KCHEMTOP:KMAX_MID), apHouts (KCHEMTOP:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), aNH4outs(KMAX_MID:KMAX_MID), & 
             gNH3out(KCHEMTOP:KMAX_MID), gNH3outs(KMAX_MID:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID), gNO3outs(KMAX_MID:KMAX_MID), &
             aNAout(KCHEMTOP:KMAX_MID),  aNAouts(KMAX_MID:KMAX_MID),  &
             aCLout(KCHEMTOP:KMAX_MID),  aCLouts(KMAX_MID:KMAX_MID),  &
             gCLout(KCHEMTOP:KMAX_MID),  gCLouts(KMAX_MID:KMAX_MID),  &
             gSO4out(KCHEMTOP:KMAX_MID), gSO4outs(KMAX_MID:KMAX_MID), &

             rlhum(KCHEMTOP:KMAX_MID), tmpr(KCHEMTOP:KMAX_MID),  &
             rlhums(KMAX_MID:KMAX_MID),tmprs(KMAX_MID:KMAX_MID)

 !-----------------------------------


  if ( debug_flag  ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4_ix,20),xn_2d(HNO3_ix,20),&
               xn_2d(NH3_ix,20),xn_2d(NO3_f_ix,20),xn_2d(NH4_f_ix,20)
  end if

!//.... molec/cm3 -> micromoles/m**3
    so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4_ix,KCHEMTOP:KMAX_MID)  *1.e12/AVOG
    hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3_ix,KCHEMTOP:KMAX_MID) *1.e12/AVOG
    nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3_ix,KCHEMTOP:KMAX_MID)  *1.e12/AVOG 
    no3in(KCHEMTOP:KMAX_MID)  = xn_2d(NO3_f_ix,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(NH4_f_ix,KCHEMTOP:KMAX_MID)*1.e12/AVOG
  if ( iSeaSalt > 0 .and. AERO%INTERNALMIXED) then
    NAin(KCHEMTOP:KMAX_MID)  = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*1.e12/AVOG* species(iSeaSalt)%molwt * 0.30 / MWNA
    CLin(KCHEMTOP:KMAX_MID)  = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*1.e12/AVOG* species(iSeaSalt)%molwt * 0.55 / MWCL
    CAin(KCHEMTOP:KMAX_MID) = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*1.e12/AVOG* species(iSeaSalt)%molwt * 0.01 / MWCA
    MGin(KCHEMTOP:KMAX_MID) = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*1.e12/AVOG* species(iSeaSalt)%molwt * 0.04 / MWMG
    Kin (KCHEMTOP:KMAX_MID) = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*1.e12/AVOG* species(iSeaSalt)%molwt * 0.01 / MWK
   else
    NAin(KCHEMTOP:KMAX_MID)  = 0.0
    CLin(KCHEMTOP:KMAX_MID)  = 0.0
    CAin(KCHEMTOP:KMAX_MID)  = 0.0
    MGin(KCHEMTOP:KMAX_MID)  = 0.0
    Kin (KCHEMTOP:KMAX_MID)  = 0.0
   endif


 !===  For ambient conditions =====================

    rlhum(KCHEMTOP:KMAX_MID) = min (0.999, rh(:))
    tmpr(KCHEMTOP:KMAX_MID)  = temp(:)

 !--------------------------------------------------------------------------                
 
!    call EQSAM4clim (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rlhum, tmpr,  & !aH2Oin, &
!                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,            &
!                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout, apHout, KCHEMTOP, KMAX_MID)
!ST:1.05.2023
   call EQSAM4clim (so4in,hno3in,no3in,nh3in,nh4in,NAin,CLin,CAin,MGin,Kin,rlhum,tmpr,   &
                     aSO4out, aNO3out, aNH4out, aCLout, gSO4out, gNH3out, gNO3out, gCLout, &
                     aH2Oout, apHout, KCHEMTOP,KMAX_MID)
 
 !.. Previous call with rh, temp - not sure if makes difference ...........................   
 !   call EQSAM4clim (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh, temp,  & !aH2Oin, &
 !                    aSO4out, aNO3out, aNH4out, aNAout, aCLout,            &
 !                    gSO4out, gNH3out, gNO3out, gClout, aH2Oout, apHout, KCHEMTOP, KMAX_MID)
 !--------------------------------------------------------------------------

!//.... micromoles/m**3  -> molec/cm3 
!      xn_2d(NO3,KCHEMTOP:KMAX_MID)  = FLOOR !different for ACID/OZONE

      xn_2d(HNO3_ix,KCHEMTOP:KMAX_MID)  = max(FLOOR,gNO3out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
      xn_2d(NH3_ix,KCHEMTOP:KMAX_MID)   = max(FLOOR,gNH3out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
      xn_2d(NO3_f_ix,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNO3out(KCHEMTOP:KMAX_MID)*AVOG*1.e-12 ) 
      xn_2d(NH4_f_ix,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNH4out(KCHEMTOP:KMAX_MID)*AVOG*1.e-12 )
!      xn_2d(SO4_ix,KCHEMTOP:KMAX_MID)   = max(FLOOR,aSO4out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
!      H2O_eqsam(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )

!//....aerosol water (ug/m**3) 

      PM25_water(i,j,KCHEMTOP:KMAX_MID) = max(FLOOR, aH2Oout(KCHEMTOP:KMAX_MID) )

!//... apHouts [mole/m3] and aerosol water (ug/m**3) 
!      pH (i,j) = -log10 ( max(FLOOR,apHout(KMAX_MID))/PM25_water(i,j,KMAX_MID) ) - 9.0
      pH (i,j) = apHout(KMAX_MID)

 if ( debug_flag ) then ! Selected debug cell
    write(*,*)'After EQSAM',xn_2d(SO4_ix,KMAX_MID),xn_2d(HNO3_ix,KMAX_MID), &
               xn_2d(NH3_ix,KMAX_MID),xn_2d(NO3_f_ix,KMAX_MID),xn_2d(NH4_f_ix,KMAX_MID), &
               PM25_water(i,j,KMAX_MID), pH (i,j)
  end if


 !===  PM water for Rh=50% and T=20C  conditions ================================

!//.... molec/cm3 -> micromoles/m**3
      so4ins  = so4in(KMAX_MID:KMAX_MID)
      hno3ins = hno3in(KMAX_MID:KMAX_MID)
      nh3ins  = nh3in(KMAX_MID:KMAX_MID)
      no3ins  = no3in(KMAX_MID:KMAX_MID)
      nh4ins  = nh4in(KMAX_MID:KMAX_MID)
  if ( iSeaSalt > 0 .and. AERO%INTERNALMIXED) then
      NAins   = NAin(KMAX_MID:KMAX_MID)
      CLins   = CLin(KMAX_MID:KMAX_MID)
      CAins   = CAin(KMAX_MID:KMAX_MID)
      MGins   = MGin(KMAX_MID:KMAX_MID)
      Kins    = Kin(KMAX_MID:KMAX_MID)
  else
      NAins   = 0.0
      CAins   = 0.0
      MGins   = 0.0
      Kins    = 0.0
      CLins   = 0.0
  endif

!//.... Only for the lowest layer KMAX_MID

      rlhums = 0.5
      tmprs  = 293.15

 !--------------------------------------------------------------------------                
   
!    call EQSAM4clim (so4ins, hno3ins,no3ins,nh3ins,nh4ins,NAins,CLins, rlhums, tmprs,  & 
!                     aSO4outs, aNO3outs, aNH4outs, aNAouts, aCLouts,               &
!                     gSO4outs, gNH3outs, gNO3outs, gClouts, aH2Oouts, apHouts, KMAX_MID, KMAX_MID)     

   call EQSAM4clim (so4ins, hno3ins,no3ins,nh3ins,nh4ins,NAins,CLins,CAins,MGins,Kins, rlhums, tmprs,  & 
                     aSO4outs, aNO3outs, aNH4outs, aCLouts,gSO4outs, gNH3outs, gNO3outs, gClouts,  &
                     aH2Oouts, apHouts, KMAX_MID, KMAX_MID)

 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 

      PM25_water_rh50 (i,j)  = max(0., aH2Oouts(KMAX_MID) )

  if ( debug_flag ) then ! Selected debug cell
    write(*,*)'EQSAM PMwater_rh50', PM25_water_rh50 (i,j)
  end if

 end subroutine emep2EQSAM

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
  if(AERO%EQUILIB/='MARS' .and. AERO%EQUILIB/='MARS_2900' .and. AERO%EQUILIB/='GEOSCHEM')then
     PM25_water(i,j,:) = 0.0
     return
  endif
 
   coef = 1.e12 / AVOG

 !.. PM2.5 water at ambient conditions (3D)
   rlhum(:) = rh(:) 
   tmpr(:)  = temp(:)

    do k = KCHEMTOP, KMAX_MID
  
!//.... molec/cm3 -> ug/m3
      so4in  = xn_2d(SO4_ix,k) * species(SO4_ix)%molwt  *coef
      hno3in = xn_2d(HNO3_ix,k)* species(HNO3_ix)%molwt *coef 
      nh3in  = xn_2d(NH3_ix,k) * species(NH3_ix)%molwt  *coef
      no3in  = xn_2d(NO3_f_ix,k) * species(NO3_f_ix)%molwt  *coef
      nh4in  = xn_2d(NH4_f_ix,k) * species(NH4_f_ix)%molwt  *coef


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
      so4in  = xn_2d(SO4_ix,k) * species(SO4_ix)%molwt  *coef *cfac(SO4_ix-NSPEC_SHL,i,j) 
      hno3in = xn_2d(HNO3_ix,k)* species(HNO3_ix)%molwt *coef *cfac(HNO3_ix-NSPEC_SHL,i,j)
      nh3in  = xn_2d(NH3_ix,k) * species(NH3_ix)%molwt  *coef *cfac(NH3_ix-NSPEC_SHL,i,j)
      no3in  = xn_2d(NO3_f_ix,k) * species(NO3_f_ix)%molwt *coef *cfac(NO3_f_ix-NSPEC_SHL,i,j)
      nh4in  = xn_2d(NH4_f_ix,k) * species(NH4_f_ix)%molwt *coef *cfac(NH4_f_ix-NSPEC_SHL,i,j)
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

  FUNCTION SAFELOG10( X ) RESULT ( SAFLOG )
  !  DESCRIPTION: Calculates the LOG (base 10) of a number X for use in pH
  !  calcuations.  Returns a minimum value if X is too small, in order to 
  !  avoid NaN or Infinity problems.
  !
      REAL, INTENT(IN) :: X        ! Argument for LOG10 function
      REAL             :: SAFLOG   ! LOG10 output 
  !  
      IF ( X <= 1e-20 ) THEN
        SAFLOG = -20   ! if X<0, make pH 20. Also avoids issues with negative H+ concentrations.
      ELSE
        SAFLOG = LOG10(X)
      ENDIF
  !  
  END FUNCTION SAFELOG10

end module AerosolCalls
