! <DustProd_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
! <DustProd_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 

 module DustProd_mod

  !-- Under developments and testing ----------------
  !   Calculates wind blown dust production based on works by
  !   Zender et al. (2003). JGR, 108, 4416-4437
  !   Alfaro & Gomes (2001). JGR, 106, 18075-18084 
  !   Gomes et al. (2003). Catena, 52, 257-271.

!!!! IMPORTANT: 
!    SoilWater is available from IFS, but (due to uncertain quality
!    and inconsistency btw soil properties in IFS and sand and clay 
!    content from JRC used here) in this version 
!    foundSoilWater = .false. from Met_mod.f90 (see comments below)   

! REFS: FMB99 = F\'ecan, F., Marticorena, B. and Bergametti,
! G. (1999). Parameterization of the increase of the aeolian erosion
! threshold wind friction velocity to soil moisture for arid and semi-arid
! areas. Ann. Geophysicae, 17,149-157.

 use Biogenics_mod,      only: EmisNat, EMIS_BioNat
 use CheckStop_mod,      only: CheckStop
 use Config_module,      only : KMAX_MID, KMAX_BND, dt_advec, METSTEP, &
                               NPROC, MasterProc, USES
 use Debug_module,       only: DEBUG_DUST
 use Functions_mod,      only: ERFfunc
 use ChemSpecs_mod,      only: species
 use GridValues_mod,     only: glat, glon, i_fdom, j_fdom 
 use GridValues_mod,     only: debug_proc, debug_li, debug_lj
 use Io_mod,             only: PrintLog, datewrite
 use Landuse_mod,        only: LandCover, NLUMAX 
 use Landuse_mod,        only: water_fraction
 use LandDefs_mod,       only:  LandType
 use LocalVariables_mod, only: Grid
 use MetFields_mod,      only: z_bnd, z_mid, u_ref, ustar_nwp, roa,    &
                                  t2_nwp, sdepth, fh, ps, surface_precip, &
                                  rho_surf, &
                                  SoilWater => SoilWater_uppr,  &
                                  foundSoilWater => foundSoilWater_uppr,    &
                                  foundws10_met, ws_10m,                  &
                                  clay_frac, sand_frac,                   & 
                                  pwp, fc, SoilWaterSource
 use MicroMet_mod,       only: Wind_at_h
 use Par_mod,            only: me,LIMAX,LJMAX
 use Par_mod,            only: limax, ljmax ! Debugging 
 use PhysicalConstants_mod, only: GRAV,  AVOG, PI, KARMAN, RGAS_KG, CP
                                  !! ECO_CROP, ECO_SEMINAT, Desert=13, Urban=17
 use SmallUtils_mod,     only: find_index
 use SubMet_mod,         only: Sub
 use TimeDate_mod,       only: daynumber
 use ZchemData_mod,      only: rcemis 
 use ZchemData_mod,      only: rh

!-----------------------------------------------
  implicit none
  private

  public ::  WindDust       

  real, private, save         :: kg2molecDU, frac_fine, frac_coar,  &
                                 help_ustar_th
  real, parameter             :: soil_dens = 2650.0  ! [kg/m3]
  logical, private, save      :: my_first_call = .true.
  logical, private, save      :: dust_found
  integer, private, save      :: ipoll
  integer, allocatable, save  :: dry_period(:,:)
  character(len=20)           :: soil_type
  real, parameter             :: SMALL=1.0e-10

 ! Indices for the species defined in this routine. Only set if found
 ! Hard-coded 2-mode for now. Could generalise later
  integer, private, save :: inat_DUf ,  inat_DUc  
  integer, private, save :: itot_DUf ,  itot_DUc
  integer, private, dimension(2), save :: dust_indices ! indices in EmisNat


  contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   subroutine WindDust (i,j,debug_flag)

   implicit none

   integer, intent(in) :: i,j                 ! coordinates of column
   logical, intent(in) :: debug_flag

!st080218      integer, parameter  :: LU_DESERT = 13    ! REMOVE HARD-CODE
   real   , parameter  :: Ro_water = 1000.0 

   real, parameter:: Dm_soil = 210.0e-6,  & ! [m] MMD of the coarsest soil  (100)
                     Z0s = Dm_soil/30.0     ! [m] Smooth roughness length MaB95 p.16426, 
                                            !     MaB97 p.4392, GMB98 p.6207

   real ::  Mflux = 0.0 
   real ::  cover, z0, vh2o_sat, gr_h2o, v_h2o, ustar_moist_cor    &
          , gwc_thr, dust_lim, soil_dns_dry, ustar_z0_cor          &
          , alfa, ustar_th, uratio, ustar, clay                    &
          , frac_fin, frac_coa, flx_hrz_slt,  flx_vrt_dst

   logical :: arable, dust_prod = .false., debug
   integer :: nlu, ilu, lu

!_______________________________________________________
    if ( USES%DUST .eqv. .false. ) then
        call PrintLog("Skipping soil dust")
        return
    end if

    if ( my_first_call ) then 

        ipoll = find_index("DUST_WB_F", species(:)%name,any_case=.true. )
        dust_found = .true.
        inat_DUf  = find_index( "DUST_WB_F", EMIS_BioNat(:),any_case=.true. )
        if( ipoll < 1 .or. inat_DUf < 1 ) then
            call PrintLog( "WARNING: Dust asked for, but not found"&
                  ,MasterProc)
            dust_found = .false.
        else

           if(MasterProc) write(6,*)'***  Call for init_dust  ***   '

           call init_dust

           kg2molecDU = 1.0e-3 * AVOG / species(ipoll)%molwt

        end if
        if(DEBUG_DUST.and.MasterProc) print *, "DUSTI ", ipoll, dust_found, debug_proc

          my_first_call = .false.

    end if !  my_first_call
!_______________________________________________________


 !++++++++++++++++++++++++++++
 if ( .not. dust_found  .or.  & 
     (glat(i,j)>65.0 .and. glon(i,j)>50.0)) then  ! Avoid dust production in N. Siberia
       if( inat_DUf>0) EmisNat( inat_DUf,i,j) = 0.0
       if( inat_DUc>0) EmisNat( inat_DUc,i,j) = 0.0
       if( itot_DUf>0) rcemis( itot_DUf,KMAX_MID) = 0.0
       if( itot_DUf>0) rcemis( itot_DUc,KMAX_MID) = 0.0
    return 
 end if
 !++++++++++++++++++++++++++++



!/..  landuse types ..........

  EmisNat(dust_indices ,i,j) = 0.0  
  flx_vrt_dst = 0.0   ! vertical dust flux 
  Mflux = 0.0  

   debug = ( DEBUG_DUST .and. debug_flag)
   if( debug ) write(6,"(a,2i5,a,i5,es12.3)")'***   WINDBLOWN DUST', &
       i_fdom(i), j_fdom(j), '>> RAIN >>', dry_period(i,j),surface_precip(i,j)

!/.. Crude assumption: dust generation if Pr < 2.4 mm per day (surpace_prec[mm/h])
 NO_PRECIP:  if (surface_precip(i,j) < 0.1) then
                  dry_period(i,j) = dry_period(i,j) + 1

   if(debug) write(6,'(a30,i5,es12.3)')'>> NO RAIN >>', &
                                   dry_period(i,j),surface_precip(i,j)

!/.. Crude assumption: Dry soil after 48 hours since precipitation
  DRY:  if (dry_period(i,j)*dt_advec/3600.0 >= 48.0) then

     if(debug) &
     write(6,'(a30,f10.4)')'>> DRY period>>', dry_period(i,j)*dt_advec/3600.0

!/.. No dust production when soil is frozen, or covered by snow,
!/.. or wet (crude approximation by surface Rh)

  FROST: if ( Grid%t2C > SMALL    .and.  Grid%sdepth <  SMALL   .and.  & 
                                    rh(KMAX_MID) < 0.85)  then

    if(debug) write(6,'(a25,2f10.2,i4)')   &
          '>> FAVOURABLE for DUST>>', Grid%t2C, rh(KMAX_MID), Grid%sdepth 

   if ( water_fraction(i,j)  > 0.99) then ! skip dust calcs
      if(DEBUG_DUST) call datewrite("DUST: Skip SEA! ", &
        (/ i_fdom(i), j_fdom(j) /), (/ fc(i,j), water_fraction(i,j) /) )
      return
   end if

!==  Land-use loop ===============================================

 nlu = LandCover(i,j)%ncodes
 LUC: do ilu= 1, nlu
        lu      = LandCover(i,j)%codes(ilu)
        cover   = LandCover(i,j)%fraction(ilu)
        arable = .false.

!/.. Dust production from arables (crops outside the growing period)
   if ( LandType(lu)%is_crop) then  
     if (daynumber <  LandCover(i,j)%SGS(ilu) .or.   &
         daynumber >  LandCover(i,j)%EGS(ilu) )      &
             arable = .true.
   end if

!/.. Dust erosion from Crops (Arable) and Desert (Mediterr.Scrubs lu==11 ???)
!/.. Creates problems on e.g. Greenland, as Bare land is included in Desert

!st080218   DUST: if (arable .or. lu == LU_DESERT) then
   DUST: if (arable .or. LandType(lu)%is_desert) then
     if( debug ) write(6,'(a,i5,f10.3)')   &
                         'DUST: -----> Landuse', lu, cover

     flx_hrz_slt = 0.0

!//  Sandblasting coefficient  Alfa [1/m] = Flux_vert/Flux_horiz
!--------------------------------------------------------------
!   Dust production is EXTREMELY sensitive to this parameter, which changes 
!   flux by 3 orders of magnitude in 0.0 < fr_clay < 0.20
!   However, larger clay content decreases sandblasting effect due to coalescing
!   crust formation, which traps loose soil aggregates (Gomes etal. Catena, 2003)  
!   (CHIMERE uses C = alfa*crust_fac = 5e-5 * 4e-3 = 2e-7)
!--------------------------------------------------------------

!//__ Clay content in soil (limited to 0.2 ___
     clay = min (clay_frac(i,j), 0.20) ! [frc]

!== Testing ==== 
!    ln10 = log(10.0) 
!    alfa = 100.0 * exp(ln10 * (13.4 * clay - 6.0))  ! [m-1]     MaB95 p. 16423 (47)
!    if (DEBUG_DUST .and. debug_flag)  &
!    write(6,'(a15,2f10.5,a10,e10.3)') '>>  CLAY  >>',clay_frac(i,j),  &
!                                       clay, '  Alfa: ', alfa

!! >>>>>  For now values to ALFA are assigned below 


!//__ Threshold U* for saltation for particle (D_opt=70um) - dry ground

     ustar_th = help_ustar_th / sqrt(roa(i,j,KMAX_MID,1)) ! [m s-1] 

     if( debug) write(6,'(a,3f15.5)')  &
         'DUST: >> U*/U*th/ro >>', ustar,ustar_th,roa(i,j,KMAX_MID,1)


!//___  Inhibition of saltation by soil moisture (Fecan et al., 1999)

  ustar_moist_cor = 1.0

!//__ Minimal soil moisture at which U*_thresh icreases
     if(SoilWaterSource /= "IFS")then
        gwc_thr = (0.17 + 0.14 * clay)* clay    ! [kg/kg] [m3 m-3] 
        gwc_thr = max ( gwc_thr, 0.1)   ! Lower threshold limit (Vautard, AE, 2005)
     else
        !use a threshold consistent with the one IFS uses
        gwc_thr=pwp(i,j)
     end if
 

  if (foundSoilWater) then    ! Soil Moisture in met data

!__ Saturated volumetric water content (sand-dependent)
     vh2o_sat = 0.489-0.126 * sand_frac(i,j)       ! [m3 m-3] 

!__ Bulk density of dry surface soil [Zender, (8)]
     soil_dns_dry = soil_dens * (1.0 - vh2o_sat)    ! [kg m-3]

!__ Volumetric soil water content [m3/m3] 
!Now we have soil moisture index, SMI = (SW-PWP)/(FC-PWP)
!     v_h2o = SoilWater(i,j,1)
!-- Note, in HIRLAM SoilWater is in [m/0.072m] -> conversion
!   if ( SoilWater_in_meter )   v_h2o = SoilWater(i,j,1)/0.072 
!
! *** BUT *** the following is using IFS pwp. This is likely
! the best we can do for the grid-average, but for dust it might be more
! appropriate to calculate for specific landcover PWP values.
! 
! (Note, v_h2o should not end up negative here, see Met_mod.f90)

  v_h2o = pwp(i,j) + SoilWater(i,j,1) * (fc(i,j)-pwp(i,j) )

  if( v_h2o < SMALL ) then
   print "(a,2i4,9f10.4)"," DUSTY DRY!!",  i_fdom(i), j_fdom(j), &
      v_h2o, pwp(i,j), fc(i,j), SoilWater(i,j,1), water_fraction(i,j)
     !v_h2o = max( 1.0e-12, v_h2o) 
   call CheckStop(v_h2o <= 0.0 ,  "DUSTY DRY" )
  end if
  if( v_h2o > fc(i,j) + 0.00001  ) then
   print "(a,2i4,9f10.4)"," DUSTY WET!!",  i_fdom(i), j_fdom(j), &
    v_h2o, pwp(i,j), fc(i,j), SoilWater(i,j,1), water_fraction(i,j)

    call CheckStop(v_h2o > fc(i,j),  "DUSTY WET" )

  end if

!__ Gravimetric soil water content [kg kg-1]  
        gr_h2o = v_h2o * Ro_water/soil_dns_dry       

!__ Put also gwc_thr (=pwp) in same unit
        if(SoilWaterSource == "IFS")then
           gwc_thr=gwc_thr* Ro_water/soil_dns_dry
        end if

! Soil water correction
     if (gr_h2o > gwc_thr) &
       ustar_moist_cor = sqrt(1.0 + 1.21 *  &
                 (100.*(gr_h2o - gwc_thr))**0.68) ! [frc] FMB99 p.155(15) 

     if( debug ) then
        write(6,'(a,f8.2,2f12.4)') 'DUST: Sand/Water_sat/soil_dens',  &
                                   sand_frac(i,j),vh2o_sat,soil_dns_dry
        write(6,'(a,3f15.5)') 'DUST: SW/VolW/GrW/ ',SoilWater(i,j,1),v_h2o,gr_h2o
        write(6,'(a,3f15.5)') 'DUST: SW COMPS ',SoilWater(i,j,1), fc(i,j), pwp(i,j)
        write(6,'(a,2f10.4)') 'DUST >> U*_moist_corr  >>',gwc_thr, ustar_moist_cor
     end if
 
  else  !.. No SoilWater in met.input; Moisture correction for U*t will be 1.

!... Correction to U*th can be assigned dependent on the typical soil wetness:
! 1.1-1.75 for sand (gr_h2o = 0.1-2%) ; 1.5-2.5 for loam (gr_h2o = 4-10%) and
!                                               for clay (gr_h2o = 9-15%)

     if( debug ) then
       write(6,'(a,f8.2,2f12.4)') 'DUST ++ No SoilWater in meteorology++'
       write(6,'(a,f10.4)')  'DUST: >> U*_moist_corr  >>', ustar_moist_cor                                 
     end if 

  end if

! ===================================
  
!// Limitation of available erodible elements and aerodynamic roughness length
!st080218  if (lu == LU_DESERT)       then      ! ---------  desert -----
  if (LandType(lu)%is_desert) then
        soil_type = 'Saharan desert'
        z0 = 0.5e-4        !TEST
        dust_lim = 0.3 ! rep15   ! 0.5 !- with my versions
        alfa = 2.0e-5 
 !//  European deserts: lat(gb), long (gl)
     if ( (glat(i,j) > 36.0 .and. glon(i,j) < 0.0) .or.   & 
          (glat(i,j) > 34.0 .and. glon(i,j) < 45.0) )    then 
!rep15          (glat(i,j) > 37.0 .and. glon(i,j) < 45.0) )    then 
        soil_type = 'European Arid'
        z0 = 0.5e-4        !TEST
        dust_lim = 0.05
        alfa = 1.3e-5 
!        alfa = 1.5e-5    ! As for TFMM spring 2005
     end if

!// limit emissions in the Spanish desert grid (covered with greenhouses??)
!     if ( (i_glob(i) == 102.0 .and. j_glob(j) == 18.0) )  & 
!         dust_lim = 0.02
!         alfa = 1.3e-5 
!     elseif (lu == 6)  then           ! ------  Mediter. crops --

 !// Crops
  else                                 ! ------  Crops -------  
        soil_type = 'Crops'
        z0  = max (0.1 * LandCover(i,j)%hveg(ilu), 0.001)
        dust_lim = 0.02
        alfa = 1.0e-5  !1.e-5 
  end if

!  else                                 ! ---------  temp/root crops ------
!        z0  = max (0.1 * landuse_hveg(i,j,ilu), 0.001)
!        dust_lim = 0.01
!        alfa = 0.8e-5  !1.e-5 


!//__ Inhibition of threshold velocity by roughness elements
!//   Roughness length momentum for erodible surfaces MaB95 p.16420, GMB98 p.6205

   ustar_z0_cor = 1.0 -          & ! [frc] MaB95 p. 16420, GMB98 p. 6207
                  log(z0 / Z0s) / log( 0.35 * (0.1/Z0s)**0.8 )

   if( debug ) write(6,'(a,es12.3)') '>> ustar_zo_corr  >>', ustar_z0_cor

   ustar_z0_cor = min ( 1.0, max(0.0001,ustar_z0_cor) )

   call CheckStop(ustar_z0_cor <= 0.0 .or. ustar_z0_cor > 1.0,  &
                        "DUST ERROR: ustar_z0_cor out of range")

   ustar_z0_cor = 1.0 / ustar_z0_cor   ! [frc] 

   if( debug ) write(6,'(a,3es12.3)') '>> U*_zo_corr  >>',z0, Z0s, ustar_z0_cor


!//___   Final threshold friction velocity

!TEST       ustar_th =0.25  ! TEST

   ustar_th =                 & ! [m s-1] Threshold friction velocity for saltation
            ustar_th        * & ! [m s-1] Threshold for dry, flat ground
            ustar_moist_cor * & ! [frc] Adjustment for moisture
            ustar_z0_cor        ! [frc] Adjustment for roughness


!//__ Account for wind gustiness under free convection (Beljaars,QJRMS,1994)
!!!!!!!  Under Testing  !!!!!!!
!
!   if (foundws10_met) then
!       u10=ws_10m(i,j,1) 
!   else 
!       u10 = Wind_at_h (Grid%u_ref, Grid%z_ref, Z10, Sub(lu)%d,   &
!                           Sub(lu)%z0,  Sub(lu)%invL)
!   end if
!   
!   u10g_2 = u10*u10 + 1.44 *Grid%wstar*Grid%wstar
!
!   if ( u10g_2 > 0.0 ) then
!       u10_gust = sqrt(u10*u10 + 1.44 *Grid%wstar*Grid%wstar)
!   else
!       u10_gust = u10
!   end if
!
!   ustar = KARMAN/log(10.0/z0) *   &
!           sqrt(u10_gust*u10_gust + 1.44 *Grid%wstar*Grid%wstar)

!.. Gives too low U*; Sub(lu)%ustar=0.1 always-> both too low to generate dust

!//__ U* from NWP model ____

    ustar = Grid%ustar 

   if (DEBUG_DUST .and. debug_flag)  then
      write(6,'(3es12.3)') ustar_th, ustar_moist_cor, ustar_z0_cor
      write(6,'(a35,f8.3,3(a10,f6.3))') 'FINALLY U*_th= ',ustar_th,' U*=',ustar, &
                                   ' U*NWP=',Grid%ustar,' U*sub=',Sub(lu)%ustar
   end if

! >>>>>  Check for saltation to occur [Whi79 p.4648(19), MaB97 p.16422(28)]        

    if (ustar > ustar_th ) then            ! .and. ustar > 0.35) then   
         dust_prod = .true.

    if(DEBUG_DUST .and. debug_flag) write(6,'(a30,f8.2)')  &
         ' Saltation occur U*th/U* => ',   ustar_th/ustar

!----     !!!!!!!  Under Testing  !!!!!!! --------------------------
!//  Increase of friction velocity by surface roughness length 
!    and friction speeds (Owens effect)
!== Saltation roughens the boundary layer, AKA "Owen's effect"
!  GMB98 p. 6206 Fig. 1 shows observed/computed u* dependence on observed U1m
!  GMB98 p. 6209 (12) has u* in cm s-1 and U, Ut in m s-1 (pers. comm.)  
!  D. Gillette, 19990529
!  With everything in MKS, the 0.3 coefficient in GMB98 (12) becomes 0.003 
!  Increase in friction velocity due to saltation varies as square of 
!  difference between reference wind speed and reference threshold speed

!//__ Threshold 10 m wind speed [m s-1] for saltation
!    u_th10 = ustar_th/KARMAN * log(z10/z0) 
!!..or    u_th10 = u10 * ustar_th / ustar ! [m s-1]
!//__  Friction velocity increase from saltation GMB98 p. 6209
!    owens = 0.003* (u10 - u_th10)*(u10 - u_th10)   !  [m s-1]
!      ustar = ustar + owens ! [m s-1] Saltating friction velocity
!    if(DEBUG_DUST .and. debug_flag)  write(6,'(a20,3f10.3)')   &
!             'Owens effect ',ustar-owens, owens , ustar
!----------------------------------------------------------------------

!//__ Calculate U*th / U* ratio
    uratio = ustar_th / ustar 

!//___ Horizontal saltation flux [kg m-1 s-1] 

    flx_hrz_slt = ustar**3 * (1.0-uratio)*(1.0+uratio)*(1.0+uratio) &   
                    * Grid%rho_ref / GRAV 

    if (DEBUG_DUST .and. debug_flag) then
      write(6,*)' '
      write(6,'(a35,es12.3)')  ' Horizontal Flux => ',   flx_hrz_slt
      write(6,'(/a25,4f8.2)') soil_type,glat(i,j),glon(i,j)
      write(6,'(a25,3f10.3,es10.2)') '>>  U*/U*t/Klim,alfa  >>',  &
            ustar,ustar_th, dust_lim, alfa
      write(6,'(a15,f10.3,2es12.3)') 'FLUXES:',uratio, flx_hrz_slt*1000.0,  &
            flx_hrz_slt*dust_lim*alfa *1000.0
    end if

!TEST  to limit the dust production
!     if (lu == LU_DESERT)  then 
!        flx_hrz_slt  = min(10.e-3, flx_hrz_slt* dust_lim )
!     else
!        flx_hrz_slt  = min(2.e-3, flx_hrz_slt* dust_lim  )
!     end if
!        flx_vrt_dst =  flx_vrt_dst + flx_hrz_slt * alfa * cover

!//  Vertical dust flux [kg/m2/s], scaled with area fraction and
!//  added for erodible landuses. 
!//  (dust production is limited to reasonable (estimated) values
!st080218   if (lu == LU_DESERT)  then
   if (LandType(lu)%is_desert) then
 
      flx_vrt_dst  =  flx_vrt_dst + min(1.e-7, flx_hrz_slt * alfa * dust_lim)
   else
      flx_vrt_dst  =  flx_vrt_dst + min(1.e-8, flx_hrz_slt * alfa * dust_lim)
   end if

   flx_vrt_dst =  flx_vrt_dst * cover

!.. Mass flus for test output
   Mflux = Mflux + flx_hrz_slt * alfa * dust_lim

   if( debug )  then
    write(6,'(a35,es12.3/)')  ' Vertical Flux   => ',  Mflux
    write(6,'(a35,es12.3,i4,f8.3)')  'DUST Flux   => ',   flx_vrt_dst, lu, cover
    write(6,'(a15,f10.3,2es12.3)') 'FLUXES:',uratio, flx_hrz_slt*1000., flx_vrt_dst*1000.
   end if

   end if  ! U* > U*_threshold

   end if DUST

   end do LUC ! landuse

  end if FROST
  end if DRY 

 else  ! PREC
    dry_period(i,j) = 0

    if( debug ) write(6,'(a30,i5,es12.3)')   &
        '>> RAIN-RAIN >>', dry_period(i,j),surface_precip(i,j)
 end if NO_PRECIP

!//__  N production [ 1/m2/s]:   d3(mkm->m) * 1e-18
!TEST       Nflux(n) = Mflux(n) *m_to_nDU / dsoil(n)**3 *1.e18

  if (dust_prod) then

!.. Vertical dust mass flux: fine and coarse [kg/m2/s] -> [ng/m3/s]
 
! TEST: stuff below needs to be tested
!   call get_dustfrac(frac_fin, frac_coa, ustar) 

! fine and coarse dust fractions assigned now  loosely based on Alfaro&Gomes(2001)
     frac_fine = 0.05    ! fine fraction 0.10 was found too large 
     frac_coar = 0.20    ! coarse fraction also 0.15-0.23 was tested
!!.. vertical dust flux [kg/m2/s] -> [kg/m3/s]*AVOG/M e-3 -> [molec/cm3/s]  

      rcemis( itot_DUf, KMAX_MID ) = frac_fine * flx_vrt_dst * kg2molecDU /Grid%DeltaZ
      rcemis( itot_DUc, KMAX_MID ) = frac_coar * flx_vrt_dst * kg2molecDU /Grid%DeltaZ

! Need kg/m2/hr for EmisNat
      EmisNat( inat_DUf,i,j) = frac_fine * flx_vrt_dst * 3600.0
      EmisNat( inat_DUc,i,j) = frac_coar * flx_vrt_dst * 3600.0

!//__Dust flux [kg/m2/s] for check
 
     dust_prod = .false.   ! Zero-setting

    if( debug ) write(6,'(//a15,2es12.4,a15,e12.4,2f6.3)') &
       '<< DUST OUT>>', EmisNat( inat_DUf,i,j), EmisNat( inat_DUc,i,j), &
       ' > TOTAL >',  sum( EmisNat( dust_indices, i,j )),frac_fin, frac_coa

  end if  ! dust_prod

  if( debug ) write(6,*) '<< No DUST production TOTAL >>', sum( EmisNat( dust_indices, i,j ))

   end subroutine WindDust

! >=================================================================<


! <=================================================================>

  subroutine init_dust

!_____ Initialization, assignments etc.
!  Calculation of fine/coarse dust fractions - not used yet
 
  implicit none

  integer :: i
  real    :: Re_opt, k_help1, k_help2, help_ust
  real, parameter   :: D_opt = 75.0e-6  ! [um] Optimal particle size for uplift
  integer, parameter  :: Nsoil = 3, Ndust = 4
  real    :: sqrt2, x1,x2, y, tot_soil
  integer :: isoil, idu
  real, dimension(Nsoil,Ndust)  :: help_diff
  real, dimension(Ndust)    :: sum_soil
  real, dimension (Ndust) :: d1 = (/0.1, 1., 2.5, 5. /)
  real, dimension (Ndust) :: d2 = (/1., 2.5, 5., 10. /)
!//__ different soil size distributions are tested
  real, dimension (Nsoil) :: dsoil   = (/ 0.832,  4.82,  19.38 /)
  real, dimension (Nsoil) :: mass_fr = (/ 0.036,  0.957, 0.007 /) 
  real, dimension (Nsoil) :: sig_soil= (/ 2.10,   1.9 ,  1.60  /)
!  real, dimension (3) :: mass_fr = (/ 2.6e-6, 0.781,  0.219/) 
!  real, dimension (3) :: dsoil   = (/ 0.0111, 2.524, 42.10 /)  ! [um] MMD 
!  real, dimension (3) :: sig_soil= (/ 1.89 ,  2.0 ,   2.13 /)  ! [frc] GSD
!  real, dimension (3) :: dsoil   = (/  1.5 ,  6.7,   14.2  /)  ! [um] MMD 
!  real, dimension (3) :: sig_soil= (/  1.7 ,  1.6 ,  1.5   /)  ! [frc] GSD
!---------------------------------------------

  inat_DUf = find_index( "DUST_WB_F", EMIS_BioNat(:),any_case=.true. )
  inat_DUc = find_index( "DUST_WB_C", EMIS_BioNat(:),any_case=.true. )
  itot_DUf = find_index( "DUST_WB_F", species(:)%name,any_case=.true.    )
  itot_DUc = find_index( "DUST_WB_C", species(:)%name,any_case=.true.    )

  if ( itot_DUf < 1 ) then
      write(*,*) "DUST species not found!"
      return
  end if

  dust_indices = (/   inat_DUf,  inat_DUc /)  

  if (DEBUG_DUST .and. MasterProc) then
   write(6,*)
   write(6,*) ' >> DUST init <<',soil_dens,  inat_DUf,  inat_DUc , itot_DUf, itot_DUc
  end if


  allocate(dry_period(LIMAX, LJMAX))
  dry_period = 72


!//__ Reynold's number ( uses D_opt[cm] ) 
    Re_opt=0.38 + 1331. *(100.*D_opt)**1.56  ! [frc] "B" MaB95 p. 16417 (5)

 !.. Given Re*t(D_opt), compute time independent factors contributing to u*t

!//__ [frc] IvW82 p. 115 (6) MaB95 p. 16417 (4) Inter-particle cohesive forces
    k_help1 = 1.0 + 6.e-7 / ( soil_dens *GRAV *(D_opt**2.5) )      !   SQUARED      
    k_help2 = soil_dens * GRAV * D_opt
                             !   SQUARED
   if (DEBUG_DUST .and. MasterProc) write(6,'(a,f6.1,3e12.4)')  &
               'DUST:ROsoil/Re_opt/K1/K2 ',soil_dens,Re_opt, k_help1, k_help2

!//__ U_star_threshold
    if ( Re_opt < 0.03 ) then
      call CheckStop( 'ERROR: Dust: Reynolds < 0.03' )

    else if ( Re_opt < 10.0 ) then
      help_ust = 1.928 * Re_opt**0.0922 - 1.0 ! [frc] IvW82 p. 114 (3), MaB95 p. 16417 (6)
      help_ust = 0.129 * 0.129 /  help_ust ! [frc]                       SQUARED

    else 
      help_ust = 1.0- 0.0858 * exp(-0.0617 *(Re_opt-10.0)) ! [frc] IvW82 p. 114(3), 
                                                           ! MaB95 p. 16417 (7)
      help_ust = 0.12*0.12 * help_ust*help_ust  ! [frc]                  SQUARED

    end if     ! Re_opt < 0.03

!//__ This method minimizes the number of square root computations performed

    help_ustar_th = sqrt (k_help1 * k_help2 * help_ust)

    if (DEBUG_DUST .and. MasterProc) write(6,'(a,2e12.4)') 'DUST: U*, U*t help =',&
           help_ust, help_ustar_th

!// =======================================
!TEST       sand_frac = 0.6   ! sand fraction in soil
!// Saturated volumetric water content (sand-dependent)
!     vh2o_sat = 0.489-0.126 * sand_frac(i,j)       ! [m3 m-3] 
!// Bulk density of dry surface soil [Zender, (8)]
!     soil_dns_dry = soil_dens * (1.0 - vh2o_sat)    ! [kg m-3]

!//__ Calculate mass dust fractions : fine and coarse - JUST BEING TESTED!!!!!

    if (DEBUG_DUST .and. MasterProc) then
      write(6,*)
      write(6,*) 'DUST: >> fractions <<', Nsoil, Ndust
      write(6,'(a,3e12.4)') 'DUST: Sigma =', (sig_soil(i),i=1,Nsoil)
    end if

    sum_soil(:) = 0.0
    tot_soil    = 0.0

    sqrt2 = sqrt (2.0)

    do isoil = 1, Nsoil
       y = log (sig_soil(isoil) ) * sqrt2
    do idu = 1, Ndust 
       x1 = log ( d1(idu) / dsoil(isoil) ) / y
       x2 = log ( d2(idu) / dsoil(isoil) ) / y

       if (DEBUG_DUST .and. MasterProc) write (6,'(a,4e12.4)') 'DUST: TEST 3', &
                               x1,x2,ERFfunc(x1),ERFfunc(x2)

       help_diff(isoil,idu) = 0.5 * ( ERFfunc(x2) - ERFfunc(x1) )  &
                                  * mass_fr(isoil)

       sum_soil(idu) =  sum_soil(idu) + help_diff(isoil,idu)
       tot_soil = tot_soil + help_diff(isoil,idu)

     end do
   end do

   frac_fine =  sum_soil(1) + sum_soil(2) + sum_soil(3) 
   frac_coar =  sum_soil(4)

   if (DEBUG_DUST .and. MasterProc)  then
     do idu = 1,Ndust 
      write (6,'(a25,2f8.4,3(f8.3),2f12.3)') 'DUST: frac in bins:',  &
             d1(idu), d2(idu), (help_diff(isoil,idu), isoil=1,3),         &
             sum_soil(idu),sum_soil(idu)/tot_soil
     end do
     write (6,'(a,2f8.4)') 'DUST: ** FINE / COARSE **',frac_fine, frac_coar
   end if

  end subroutine init_dust
! >=================================================================<

! <=================================================================>
  subroutine get_dustfrac(frac_fine, frac_coarse, ustar)

 ! Calculates fine/coarse dust fractions dependenden on U*,
 ! based Alfaro & Gomes (2001) - JUST BEING TESTED!!!!!

   real :: frac_fine, frac_coarse, ustar, a
 
   if (ustar < 0.35) then
      frac_fine   = 0.02
      frac_coarse = 0.09
   else if (ustar < 0.40) then
      a           = (ustar-0.35)/0.05
      frac_fine   = (1.0-a)*0.02 + a*0.04
      frac_coarse = (1.0-a)*0.09 + a*0.11
   else if (ustar < 0.55) then
      a           = (ustar-0.40)/0.15
      frac_fine   = (1.0-a)*0.04 + a*0.26
      frac_coarse = (1.0-a)*0.11 + a*0.30
   else if (ustar < 0.80) then
      a           = (ustar-0.55)/0.25
      frac_fine   = (1.0-a)*0.26 + a*0.35
      frac_coarse = (1.0-a)*0.30 + a*0.11
   else
      frac_fine   = 0.35
      frac_coarse = 0.11
   end if

 end subroutine get_dustfrac

! >----------------------------------------------------------<

 end module DustProd_mod

