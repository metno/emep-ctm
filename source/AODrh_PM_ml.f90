! <AODrh_PM_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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


 module AODrh_PM_ml

  !-----------------------------------------------------------------------! 
  ! Calculates Aerosol Optical Depth (AOD) for 0.5 um radiation based on 
  ! aerosol mass concentrations and specific extinction cross-sections 
  ! based on Tegen et al. JGR (1997) and Kinne et al., ACP (2005)
  ! (implicit assumption on wetted aerosols)
  !-----------------------------------------------------------------------!
 use Chemfields_ml,        only : AOD, Extin_coeff              
 use ChemChemicals_ml,     only : species  
 use ChemGroups_ml,        only : AOD_GROUP
 use ChemSpecs_tot_ml
 use GridValues_ml,        only : i_fdom, j_fdom
 use MetFields_ml,         only : z_bnd
 use ModelConstants_ml,    only : KMAX_MID, KMAX_BND, KCHEMTOP,   &
                                  MasterProc, NPROC, DEBUG_i, DEBUG_j
 use Par_ml,               only : MAXLIMAX,MAXLJMAX   ! => x, y dimensions
 use PhysicalConstants_ml, only : AVOG
 use Setup_1dfields_ml,    only : xn_2d, rh 

  implicit none
  private
  !-----------------------------------------------------------------------!
 !// subroutines
  public ::   AOD_Ext 

  logical, private, save :: my_first_call = .true.

  contains

! <---------------------------------------------------------->

    subroutine AOD_Ext (i,j,debug)

 !------------------------------------------------
 ! Calculates AOD
 !-------------------------------------------------

 implicit none


   integer, intent(in) :: i,j    ! coordinates of column
   logical, intent(in) :: debug

   integer :: k, n, itot, irh
   real, parameter :: lambda = 0.55e-6
   real, parameter :: rhoSO4=1.6, rhoOC=1.8, rhoEC=1.0,   &
                      rhoDU=2.6,  rhoSS=2.2 
   real, parameter :: Reff_SO4 = 0.156, Reff_OC = 0.087, Reff_EC = 0.039, &
                      Reff_DUf = 0.80,  Reff_DUc = 4.5, Reff_SSf = 0.80,  &
                      Reff_SSc = 5.73
   real, parameter, dimension(7) ::    &
       RelHum = (/ 0.0,  0.5, 0.7, 0.8, 0.9, 0.95, 0.99 /),  &
       GF_SO4 = (/ 1.0,  1.4, 1.5, 1.6, 1.8,  1.9,  2.2 /),  &
       GF_OC  = (/ 1.0,  1.2, 1.4, 1.5, 1.6,  1.8,  2.2 /),  &
       GF_EC  = (/ 1.0,  1.0, 1.0, 1.2, 1.4,  1.5,  1.9 /),  &
       GF_SS  = (/ 1.0,  1.6, 1.8, 2.0, 2.4,  2.9,  4.8 /),  &
!.. 550 nm
       Ex_SO4 = (/ 1.114, 1.545, 1.742, 1.862, 2.036, 2.206, 2.558/),  &  !H2SO4
       Ex_OC  = (/ 0.560, 0.652, 0.701, 0.741, 0.821, 0.921, 1.181/),  &
       Ex_EC  = (/ 0.483, 0.484, 0.471, 0.417, 0.363, 0.343, 0.332/),  &
       Ex_SSf = (/ 2.699, 2.547, 2.544, 2.508, 2.444, 2.362, 2.221/),  &
       Ex_SSc = (/ 2.143, 2.103, 2.090, 2.106, 2.084, 2.070, 2.064/)
  real, dimension(KMAX_MID):: ext_SO4, ext_NO3, ext_NH4, ext_EC,       &
                              ext_OM, ext_SS, ext_DU, kext
  real :: gfSO4, gfOC, gfEC, gfSS, &
          rhoSO4_wet, rhoOC_wet, rhoEC_wet, rhoSS_wet, rhoNO3c_wet,  &
          massGF_SO4, massGF_OC, massGF_EC, massGF_SS,               &
          extSO4, extEC, extOC, extSSf, extSSc,    AOD_old,          &                       
          AOD_SO4, AOD_NO3, AOD_NH4, AOD_EC, AOD_POM, AOD_OM,  AOD_SS, AOD_DU, &
          SpecExt_SO4, SpecExt_OC, SpecExt_EC , SpecExt_SSf, SpecExt_SSc, &
          SpecExt_DUf, SpecExt_DUc, SpecExt_NO3f, SpecExt_NO3c, SpecExt_NH4

!-----------------------------------------------------------------
!   AOD_GROUP = (/ SO4, NO3_F, NH4_F, EC_F_NEW, EC_F_AGE, POC_F, &
!       EXTC  = (/ 8.5, 8.5,   8.5,   7.5,      11.0,     5.7,   &  &.0 for POC!
!                  SEASALT_F, SEASALT_C, DUST_NAT_F, DUST_NAT_C /)
!                  3.0,        0.4       1.0,         0.3,      /)
!__________________________________________________________________

 if(debug) write(*,*) ' #### in AOD module  ###'

  AOD(i,j)     = 0.0
  kext(:)      = 0.0


  do k =  KCHEMTOP, KMAX_MID   !_______________ vertical layer loop

!   kext(k)  = 0.0
 
!.. ===========================================================================
!..  Extinction coefficients: Kext [1/m] = SpecExtCross [m2/g] * mass [g/m3]
!..                           summed up for all components  
!..    xn_2d(spec,k)*1.e15 *species(ispec)%molwt/AVOG   [molec/m3] -> [ng/m3]
!..                                                   [ng/m3 ] * 1e-9 -> [g/m3]
!..=>  xn_2d(ispec,k) * species(ispec)%molwt * 1.e6 / AVOG  [g/m3]
!.. ===========================================================================

!//.. Calculate Rh dependent specific extinction (or mass extinction efficiency)
!..   according to Chin et.al (J. Atm.Sci., 59, 2001)


  RHloop: do irh = 2, 7
    if ( rh(k) <= RelHum (irh) ) then

!.. Find actual growth factor by interpolation
        
      gfSO4 = GF_SO4(irh-1) + ( (GF_SO4(irh) - GF_SO4(irh-1)) /    &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                      ( rh(k) - RelHum(irh-1)   )
      gfOC  = GF_OC(irh-1) + ( (GF_OC(irh) - GF_OC(irh-1)) /       &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                     ( rh(k) - RelHum(irh-1)   )
      gfEC  = GF_EC(irh-1) + ( (GF_EC(irh) - GF_EC(irh-1)) /       &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                     ( rh(k) - RelHum(irh-1)   )
      gfSS  = GF_SS(irh-1) + ( (GF_SS(irh) - GF_SS(irh-1)) /       &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                     ( rh(k) - RelHum(irh-1)   )
!.. Extinction efficiencies
      extSO4 = Ex_SO4(irh-1) + ( (Ex_SO4(irh) - Ex_SO4(irh-1)) /   &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                      ( rh(k) - RelHum(irh-1)   )
      extOC  = Ex_OC(irh-1)  + ( (Ex_OC(irh) - Ex_OC(irh-1)) /     &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                     ( rh(k) - RelHum(irh-1)   )
      extEC  = Ex_EC(irh-1)  + ( (Ex_EC(irh) - Ex_EC(irh-1)) /     &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                     ( rh(k) - RelHum(irh-1)   )
      extSSf = Ex_SSf(irh-1) + ( (Ex_SSf(irh) - Ex_SSf(irh-1)) /   &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                     ( rh(k) - RelHum(irh-1)   )
      extSSc = Ex_SSc(irh-1) + ( (Ex_SSc(irh) - Ex_SSc(irh-1)) /   &
                                (RelHum(irh) - RelHum(irh-1) ) ) * &
                                     ( rh(k) - RelHum(irh-1)   )

  if(debug .and. (k == KMAX_MID) ) write(6,'(a15,i3,5f8.2)') '## Rh >> ', &
  irh, rh(k), RelHum(irh), RelHum(irh-1), GF_SO4(irh), GF_SO4(irh-1)
   if(debug .and. (k == KMAX_MID) ) write(6,'(a15,i3,f8.2,3f10.3)')   &
   '## Ext > ', irh, rh(k),  Ex_SO4(irh-1), Ex_SO4(irh), extSO4

!.. end leave Rh loop
      exit RHloop

    else 
      gfSO4 = 1.0
      gfOC  = 1.0
      gfEC  = 1.0
      gfSS  = 1.0
      extSO4 = Ex_SO4(1)
      extOC  = Ex_OC(1)
      extEC  = Ex_EC(1)
      extSSf = Ex_SSf(1)
      extSSc = Ex_SSc(1)
      cycle 
    
    endif
  enddo RHloop

  if(debug .and. (k == KMAX_MID) ) write(6,'(a15,5f8.2)')  '## GFs =  ', &
     rh(k), gfSO4, gfOC, gfEC, gfSS
  if(debug .and. (k == KMAX_MID) ) write(6,'(a15,6f10.3)') '## ExtEff ', &
     rh(k), extSO4, extOC, extEC, extSSf, extSSc 
!.. Density of wet aerosol
!.   rho_w = Vfr_dry*Rho_dry + (1-Vfr_dry)*Rho_water 
!..  where   Vfr_dry = 1/GF**3 (dry volume fraction)

  rhoSO4_wet = rhoSO4 / gfSO4**3 + (1.0-1.0/gfSO4**3) ! *1.0 [g/cm3]
  rhoOC_wet  = rhoOC  / gfOC**3  + (1.0-1.0/gfOC**3)
  rhoEC_wet  = rhoEC  / gfEC**3  + (1.0-1.0/gfEC**3)
  rhoSS_wet  = rhoSS  / gfSS**3  + (1.0-1.0/gfSS**3)
! Fake
  rhoNO3c_wet  = rhoSO4  / gfSS**3  + (1.0-1.0/gfSS**3)

!.. Ratio mass_wet / mass_dry (Mwet/Mdry)

  massGF_SO4 = gfSO4**3 * rhoSO4_wet/rhoSO4
  massGF_OC  = gfOC**3  * rhoOC_wet/rhoOC
  massGF_EC  = gfEC**3  * rhoEC_wet/rhoEC
  massGF_SS  = gfSS**3  * rhoSS_wet/rhoSS

!.. Specific extinction [m2/g] 
!   beta = 3/4 * ExtCoef/rho_wet/rad_eff * Mwet/Mdry
  
    SpecExt_SO4  = 0.75 * extSO4 * massGF_SO4 / (rhoSO4_wet * Reff_SO4)
    SpecExt_OC   = 0.75 * extOC  * massGF_OC / (rhoOC_wet * Reff_OC)
    SpecExt_EC   = 0.75 * extEC  * massGF_EC / (rhoEC_wet * Reff_EC)
    SpecExt_SSf  = 0.75 * extSSf * massGF_SS / (rhoSS_wet * Reff_SSf)
    SpecExt_SSc  = 0.75 * extSSc * massGF_SS / (rhoSS_wet * Reff_SSc)
    SpecExt_DUf  = 0.75 * species(DUST_WB_F)%ExtC  / (rhoDU * Reff_DUf)
    SpecExt_DUc  = 0.75 * species(DUST_WB_C)%ExtC  / (rhoDU * Reff_DUc)
!... Faking
    SpecExt_NH4  =  SpecExt_SO4                        
    SpecExt_NO3f = SpecExt_SO4
!.. VERY CRUDE: Assume NOc sitting on SSc: Q and GF for SSc are applied
    SpecExt_NO3c = 0.75 * extSSc * massGF_SS / (rhoNO3c_wet * Reff_SSc)

!=========================================================

!.. Specific extinction

    do n = 1, size(AOD_GROUP)
      itot = AOD_GROUP(n)

      kext(k) = kext(k) +   &
                    xn_2d(itot,k) * species(itot)%molwt * species(itot)%ExtC      

!  if(debug .and. (k == KMAX_MID) ) write(6,'(a15,3i4,2es10.3,i5,f8.1)') '$$ AOD ', &
!     n, k, itot, kext(k), xn_2d(itot,k),species(itot)%molwt, species(itot)%ExtC

    enddo

     kext(k) = kext(k) * 1.0e6 / AVOG 

!     if(debug .and. (k == 18 .or. k == KMAX_MID) )  &
!            write(6,'(a17,i4,es15.3)') '> Ext. coeff', k, kext(k)

!.. Aerosol extinction optical depth : integral over all vertical layers
!.. [1/m} * [m]

! ... Old parameterisation... ratained for comparison
      AOD_old = AOD_old + kext(k) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))

!      if(debug .and. (k == 18 .or. k == KMAX_MID) )  & 
!      write(6,'(a25,i4,2es15.4,2f8.1)') '>> Kext AOD for layer', k,  &
!                kext(k), AOD(i,j), z_bnd(i,j,k), z_bnd(i,j,k+1)

!.. Extinction coefficients for individual components

 ext_SO4(k) =  xn_2d(SO4,k) * species(SO4)%molwt * SpecExt_SO4 * 1.0e6 / AVOG

 ext_NO3(k) = ((xn_2d(NO3_F,k) + 0.3*xn_2d(NO3_C,k)) * SpecExt_NO3f +  &
                0.7*xn_2d(NO3_C,k)                   * SpecExt_NO3c )  &
                                            * species(NO3_F)%molwt * 1.0e6 / AVOG

 ext_NH4(k) =  xn_2d(NH4_F,k) * SpecExt_NH4 * species(NH4_F)%molwt * 1.0e6 / AVOG

 ext_EC(k)  = (( xn_2d(EC_F_FFUEL_NEW,k) + xn_2d(EC_F_FFUEL_AGE,k) +           &
                    xn_2d(EC_F_WOOD_NEW,k)  + xn_2d(EC_F_WOOD_AGE,k)  )        &
                  * SpecExt_EC              * species(EC_F_FFUEL_NEW)%molwt +  &
                    xn_2d(FFIRE_BC,k) * SpecExt_EC * species(FFIRE_BC)%molwt ) &
                                                                   * 1.0e6 / AVOG
!                +( xn_2d(EC_C_FFUEL,k) + xn_2d(EC_C_WOOD,k))            &
!                     * SpecExt_ECc * species(EC_C_FFUEL)%molwt

! ext_POM(k) = ( xn_2d(POM_F_FFUEL,k) * species(POM_F_FFUEL)%molwt +    &
!                xn_2d(POM_F_WOOD,k)  * species(POM_F_WOOD)%molwt    )  &
!                                     * SpecExt_OC * 1.0e6 / AVOG
!!             + xn_2d(POM_C_FFUEL,k) * species(POM_C_FFUEL)%molwt * SpecExt_OCc
   
 ext_OM(k)  = ( xn_2d(PART_OM_F,k) * species(PART_OM_F)%molwt * SpecExt_OC   &
              + xn_2d(FFIRE_OM,k)  * species(FFIRE_OM)%molwt  * SpecExt_OC ) &
!             + xn_2d(POM_C_FFUEL,k) * species(POM_C_FFUEL)%molwt * SpecExt_OCc
                                                                  * 1.0e6 / AVOG

 ext_SS(k)  = ( xn_2d(SEASALT_F,k)  * SpecExt_SSf +   &
                xn_2d(SEASALT_C,k)  * SpecExt_SSc   ) &
                                    * species(SEASALT_F)%molwt * 1.0e6 / AVOG

 ext_DU(k)  = ( (xn_2d(REMPPM25,k) + xn_2d(DUST_WB_F,k)+ xn_2d(DUST_SAH_F,k))   &
                 * SpecExt_DUf  &
               +(xn_2d(REMPPM_C,k) + xn_2d(DUST_WB_C,k)+ xn_2d(DUST_SAH_C,k))   &
                 * SpecExt_DUc )    * species(DUST_WB_F)%molwt     * 1.0e6 / AVOG &
               + xn_2d(FFIRE_REMPPM25,k) * SpecExt_DUf                     &
                                   * species(FFIRE_REMPPM25)%molwt * 1.0e6 / AVOG


 Extin_coeff(i,j,k) =  ext_SO4(k) + ext_NO3(k) + ext_NH4(k) + ext_EC(k)   &
                     + ext_OM(k)  + ext_SS(k)  + ext_DU(k)

 if(debug .and. (k==3 .or. k==KMAX_MID) ) write(6,'(a15,i3,8es10.3)') 'EXTINCs for k =', &
    k, ext_SO4(k), ext_NO3(k), ext_NH4(k), ext_EC(k), ext_OM(k),   &
    ext_SS(k), ext_DU(k), Extin_coeff(i,j,k)

!.. Aerosol optical depth for individual components
  AOD_SO4 = AOD_SO4 + ext_SO4(k) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
  AOD_NO3 = AOD_NO3 + ext_NO3(k) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
  AOD_NH4 = AOD_NH4 + ext_NH4(k) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
  AOD_EC  = AOD_EC  + ext_EC(k)  * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
!  AOD_POM = AOD_POM + ext_POM(k) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
  AOD_OM  = AOD_OM  + ext_OM(k)  * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
  AOD_SS  = AOD_SS  + ext_SS(k)  * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
  AOD_DU  = AOD_DU  + ext_DU(k)  * (z_bnd(i,j,k)-z_bnd(i,j,k+1))


  enddo                        !_______________ vertical layer loop


  AOD(i,j) = AOD_SO4 + AOD_NO3 + AOD_NH4 + AOD_EC + AOD_OM  + AOD_SS  + AOD_DU

  if(debug )  write(6,'(a30,2i5,8es15.3)') '>>>  AOD / AODs  <<<',   &
              i_fdom(i), j_fdom(j), AOD(i,j), AOD_SO4, AOD_NO3,  &
              AOD_NH4, AOD_EC, AOD_OM, AOD_SS, AOD_DU


   end subroutine AOD_Ext

 end module AODrh_PM_ml

