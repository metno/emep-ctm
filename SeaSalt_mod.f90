! <SeaSalt_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
! <SeaSalt_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 

                          module SeaSalt_mod

!*****************************************************************************! 
! Calculates production of sea salt based on: 
! Maartinsson et al.(2003) JGR,100,D8      for particles with Ddry<1.25um  
! Monahan et al.(1986) J.Phys.Oceanogr,10  for particles with Ddry=~1.25-5um  
! This appr. corresponds sea salt aerosol at ambient Rh upto 10 um (inclusion 
! of larger particles is also possible). The sea salt aerosols are aggregated
! in fine, coarse and 'giant' size fractions
! Programmed by Svetlana Tsyro
!-----------------------------------------------------------------------------

 use AeroFunctions_mod,     only: WetRad, cmWetRad, GbSeaSalt
 use Biogenics_mod,         only: EMIS_BioNat, EmisNat  
 use CheckStop_mod,         only: StopAll
 use ChemSpecs_mod,         only: species
 use Config_module,         only: KMAX_MID, KMAX_BND, USES, MasterProc
 use Debug_module,          only: DEBUG  ! -> DEBUG%SEASALT
 use GridValues_mod,        only: glat, glon, i_fdom, j_fdom 
 use Io_Progs_mod,          only: PrintLog
 use Landuse_mod,           only: LandCover, water_fraction
 use LocalVariables_mod,    only: Grid
 use MetFields_mod,         only: u_ref, z_bnd, z_mid, sst,  &
                                  u_ref, foundSST, &
                                   foundws10_met,ws_10m
 use MicroMet_mod,          only: Wind_at_h
 use PhysicalConstants_mod, only: CHARNOCK, AVOG ,PI
 use SmallUtils_mod,        only: find_index
 use SubMet_mod,            only: Sub
 use TimeDate_mod,          only: current_date
 use ZchemData_mod,         only: rcemis 

 !-------------------------------------

  implicit none
  private

  public ::  SeaSalt_flux   ! subroutine

  integer, parameter :: SS_MAAR= 7, SS_MONA= 3, &   !Number size ranges for
                                                    !Maartinsson's and Monahan's
                        NFIN= 7, NCOA= 3      , &   !Number fine&coarse bins 
                        SSdens = 2200.0             ! sea salt density [kg/m3]
!TEST                      NFIN= 7, NCOA= 2, NGIG= 1, &   

  real, save, dimension(SS_MAAR) :: dp3, a, b
  real, save, dimension(SS_MAAR+1) :: log_dbin
  real, save, dimension(SS_MONA) :: temp_Monah, radSS, dSS3
  real, save                     :: n_to_mSS
  real, public, allocatable,dimension(:,:,:) :: SS_prod !Sea salt flux

  logical, private, save :: my_first_call = .true.
  logical, private, save :: seasalt_found

 ! Indices for the species defined in this routine. Only set if found
 ! Hard-coded for 2 specs just now. Could extend and allocate.
  integer, private, parameter :: NSS = 2
  integer, private, parameter :: iSSFI=1,  iSSCO=2
  integer, private, save :: inat_SSFI,  inat_SSCO
  integer, private, save :: itot_SSFI,  itot_SSCO

  contains

! <------------------------------------------------------------------------->

   subroutine SeaSalt_flux (i,j, debug_flag)

  !-----------------------------------------------------------------------
  ! Input: Tw        - sea surface temperature - # -
  !        u10       - wind speed at 10m height 
  ! Output: SS_prod - fluxes of fine and coarse sea salt aerosols [molec/cm3/s]
  !-----------------------------------------------------------------------

   integer, intent(in) :: i,j         ! coordinates
   logical, intent(in) :: debug_flag  ! set true for debug i,j

   real, parameter :: Z10 = 10.0  ! 10m height
   integer :: ii, jj, nlu, ilu, lu
   real    :: invdz, n2m, u10, u10_341, Tw, flux_help, total_flux, &
              whitecap
   real, save  ::   moleccm3s_2_kgm2h 
   real    :: ss_flux(SS_MAAR+SS_MONA), d3(SS_MAAR+SS_MONA) 
   real    :: rcss(NSS)
!//---------------------------------------------------
 
  if ( my_first_call ) then 

    ! We might have USES%SEASALT=.true. in ModelConstants, but the
    ! chemical scheme might not have seasalt species. We check.

    inat_SSFI = find_index( "SEASALT_F", EMIS_BioNat(:) , any_case=.true.)
    inat_SSCO = find_index( "SEASALT_C", EMIS_BioNat(:), any_case=.true. )
    itot_SSFI = find_index( "SEASALT_F", species(:)%name  , any_case=.true.  )
    itot_SSCO = find_index( "SEASALT_C", species(:)%name  , any_case=.true.  )

    if(DEBUG%SEASALT .and. MasterProc ) &
        write(*,*) "SSALT INIT", inat_SSFI, itot_SSFI

    if ( inat_SSFI < 1 ) then
       seasalt_found = .false.
       call PrintLog("WARNING: SeaSalt not found in Emis",MasterProc)
       if ( USES%SEASALT ) call StopAll('SeaSalt not found in Emis')
    else if ( itot_SSFI < 1 ) then
       seasalt_found = .false.
       call PrintLog("WARNING: SeaSalt not found in Specs",MasterProc)
       if ( USES%SEASALT ) call StopAll('SeaSalt not found in Emis')
    else
        seasalt_found = .true.
        call init_seasalt()
    end if
    
    ! For EmisNat, need kg/m2/h from molec/cm3/s
    moleccm3s_2_kgm2h =   Grid%DeltaZ * 1.0e6 * 3600.0  &! /cm3/s > /m2/hr
                          /AVOG * 1.0e-6  ! kg  after *MW
    my_first_call = .false.

  end if !  my_first_call
 !....................................

  if ( .not. seasalt_found  ) return 

 !....................................



    if ( .not. Grid%is_mainlysea .or. Grid%snowice ) then ! quick check
       EmisNat( inat_SSFI,i,j) = 0.0
       EmisNat( inat_SSCO,i,j) = 0.0
       rcemis( itot_SSFI,KMAX_MID) = 0.0
       rcemis( itot_SSCO,KMAX_MID) = 0.0
       return
    end if


!// Loop over the land-use types present in the grid

     rcss(:) = 0.0
     nlu = LandCover(i,j)%ncodes
     do ilu= 1, nlu
       lu =  LandCover(i,j)%codes(ilu)

!// only over water
! Obs!  All water is assumed here to be salt water for now 
!      (as fresh water is not distinguished in the input)

       if ( Sub(lu)%is_water ) then

          if(DEBUG%SEASALT .and. debug_flag) then
              write(6,'(a,2i4,f8.4,f12.4,3f8.3)') &
                'SSALT ** Charnock, ustar_nwp, d, Z0, SST ** ',&
                   i_fdom(i), j_fdom(j), & 
                   CHARNOCK,Grid%ustar,Sub(lu)%d,Sub(lu)%z0, sst(i,j,1)
          end if

         !.. Calculate wind velocity over water at Z10=10m 

          if(foundws10_met)then
              u10=ws_10m(i,j,1) 
          else 
              u10 = Wind_at_h (Grid%u_ref, Grid%z_ref, Z10, Sub(lu)%d,   &
                           Sub(lu)%z0,  Sub(lu)%invL)
          end if

         u10 =  max(0.1, u10)  ! make sure u10!=0 because of LOG(u10),
                               ! (use plausible physical limit for ws here)

         u10_341=exp(log(u10) * (3.41))

         !.. Alternatives for whitecap fraction (Monahan etal.1986; 
         !                     Norris eta. 2013; Callaghan etal. 2008)

    select case ( USES%WHITECAPS )
      case ( 'Monahan' )
          whitecap   = u10_341 * 3.84e-6
      case ( 'Norris' )
          whitecap   = 1.03e-5 * (u10-2.62)**3
         if (u10 <= 2.62)   &
          whitecap   = 1.0e-10
      case ( 'Callaghan')
          whitecap   = 1.0e-10
         if (u10 > 3.71 .and. u10 <= 10.18) then    !(11.25)
          whitecap   = 3.18e-5 * (u10 - 3.70)**3 
         elseif (u10 <= 23.09) then
          whitecap   = 4.82e-6 * (u10 + 1.98)**3 
         else
          whitecap   = 4.82e-6 * (23.09 + 1.98)**3 !7.594  !
         end if
    end select

         if(DEBUG%SEASALT .and. debug_flag) &
             write(6,'(a,L2,4f12.4,es14.4)')'SSALT ** U*, Uref, U10, Uh, invL ** ',&
               foundws10_met, Sub(lu)%ustar, Grid%u_ref, u10, &
               Wind_at_h (Grid%u_ref, Grid%z_ref, Z10, Sub(lu)%d,   &
                           Sub(lu)%z0,  Sub(lu)%invL), &
               Sub(lu)%invL

         !.. Sea surface temperature is not always available (e.g. pre-2001 at
         ! MET.NO), so we need an alternative. As emissions are most
         ! sensitive to u* and not T, we ignore differences between Tw and T2
         ! for the default case if SST isn't avialable.

          if ( foundSST ) then
            Tw = sst(i,j,1)
          else
            Tw = Grid%t2
          end if
          Tw = max(Tw, 270.0)! prevents unrealistic sub.zero values
          Tw = min(Tw, 300.0)! prevents unrealistic high values

! ====    Calculate sea salt fluxes in size bins  [part/m2/s] ========
         total_flux = 0.0
!... Fluxes of small aerosols for each size bin (Mårtensson etal,2004)
          do ii = 1, SS_MAAR

               flux_help  = a(ii) * Tw + b(ii)
  
               ss_flux(ii) = flux_help * ( log_dbin(ii+1) - log_dbin(ii) )    &
                                       * whitecap   ! * u10_341 * 3.84e-6 

               d3(ii) = dp3(ii)  ! diameter cubed

               total_flux =  total_flux + ss_flux(ii)

               if(DEBUG%SEASALT .and. debug_flag) write(6,'(a20,i5,es13.4)') &
                  'SSALT Flux Maarten ->  ',ii, ss_flux(ii)
          end do

!... Fluxes of larger aerosols for each size bin (Monahan etal,1986)
          do ii = 1, SS_MONA
             jj = ii + SS_MAAR

               ss_flux(jj) = temp_Monah(ii) * whitecap
!                                           * u10_341

               d3(jj) = dSS3(ii)  ! diameter cubed

               total_flux =  total_flux + ss_flux(ii) 

               if(DEBUG%SEASALT .and. debug_flag) &
                   write(6,'(a20,i5,es13.4)') 'SSALT Flux Monah ->  ',ii, ss_flux(jj)
          end do

         if(DEBUG%SEASALT .and. debug_flag) &
               write(6,'(a20,es13.3)') 'SSALT Total SS flux ->  ',  total_flux


  !ESX n2m = n_to_mSS * invdz *AVOG / species(iseasalt)%molwt *1.0e-15
  ! convert [part/m2/s] to [molec/cm3/s] required for differential equations.
  !1).  n_to_mSS =PI*SSdens/6.0  - partfrom number to mass conversion 
  !                                  (as M = PI/6 * diam^3 * density * N)
  !2).  1.0e-15 = e-18 * e3 (where e-18 is to convert [um3] to [m3] in diam^3;i
  !                                     and e3 is [kg] to [g] in density).
  !3)    (then mass is converted to [molec] with AVOG/MotW)

          invdz  = 1.0e-6 / Grid%DeltaZ       ! 1/dZ [1/cm3]

          n2m = n_to_mSS * invdz *AVOG / species(itot_SSFI)%molwt *1.0e-15

!.. Fine particles emission [molec/cm3/s] need to be scaled to get units kg/m2/s consistent with
! Emissions_mod (snapemis). Scaling factor is 
          do ii = 1, NFIN
               rcss( iSSFI) = rcss( iSSFI)  &  
                 !! ESX SS_prod(QSSFI,i,j) = SS_prod(QSSFI,i,j)   &
                                  + ss_flux(ii) * d3(ii) * n2m   &
                                  * water_fraction(i,j) 
            if(DEBUG%SEASALT .and. debug_flag) &
            write(6,'(a20,i5,2es13.4)') 'SSALT Flux fine ->  ',ii,d3(ii), rcss( iSSFI ) !ESX SS_prod(QSSFI,i,j)
          end do

!..Coarse particles emission [molec/cm3/s]
          do ii = NFIN+1, NFIN+NCOA
               rcss( iSSCO ) = rcss( iSSCO )  &
                 !!ESX SS_prod(QSSCO,i,j) = SS_prod(QSSCO,i,j)   &
                                  + ss_flux(ii) * d3(ii) * n2m   &
                                  * water_fraction(i,j)
            if(DEBUG%SEASALT .and. debug_flag) &
            write(6,'(a20,i5,2es13.4)') 'SSALT Flux coarse ->  ',ii,d3(ii), rcss( iSSCO ) !ESX SS_prod(QSSCO,i,j)
          end do

!... Crude fix for the effect of lower salinity in the Baltic Sea

          if ( (glat(i,j) > 52.0 .and. glat(i,j) < 67.0)     .and.   &  
               (glon(i,j) > 13.0 .and. glon(i,j) < 30.0) )   then 
          
               rcss( iSSFI ) = 0.2 * rcss( iSSFI )
               rcss( iSSCO ) = 0.2 * rcss( iSSCO )
          end if
  
          if(DEBUG%SEASALT .and. debug_flag) write(6,'(a35,2es15.4)')  &
             '>> SSALT production fine/coarse  >>', &
                rcss(  iSSFI ), rcss( iSSCO )
                          
       end if  ! water
     end do  ! LU classes

     EmisNat( inat_SSFI, i,j )      = rcss( iSSFI ) * moleccm3s_2_kgm2h * species( itot_SSFI )%molwt
     EmisNat( inat_SSCO, i,j )      = rcss( iSSCO ) * moleccm3s_2_kgm2h * species( itot_SSCO )%molwt
     rcemis ( itot_SSFI, KMAX_MID ) = rcss( iSSFI )
     rcemis ( itot_SSCO, KMAX_MID ) = rcss( iSSCO )

  end subroutine SeaSalt_flux


!<<---------------------------------------------------------------------------<<

  subroutine init_seasalt

  !------------------------------------------------------------
  ! Assignments and calculations of some help-parameters
  !------------------------------------------------------------

  implicit none

  integer :: i
  real    :: a1, a2
  real, dimension(SS_MONA) :: Rrange, rdry

!//===== Polynomial coeficients from Maartinsson et al. (2004)
  real, parameter, dimension(5) ::    &
        C1 = (/-2.576e35,  5.932e28, -2.867e21, -3.003e13, -2.881e6 /),  &
        C2 = (/-2.452e33,  2.404e27, -8.148e20,  1.183e14, -6.743e6 /),  &
        C3 = (/ 1.085e29, -9.841e23,  3.132e18, -4.165e12,  2.181e6 /),  &

        D1 = (/ 7.188e37, -1.616e31,  6.791e23,  1.829e16,  7.609e8 /),  &
        D2 = (/ 7.368e35, -7.310e29,  2.528e23, -3.787e16,  2.279e9 /),  &
        D3 = (/-2.859e31,  2.601e26, -8.297e20,  1.105e15, -5.800e8 /)

!=== mikrometer in powers
  real, parameter :: MKM  = 1.e-6,  MKM2 = 1.e-12 ,  &  
                     MKM3 = 1.e-18, MKM4 = 1.e-24 
 
!//.. Size bins for Maartinsson's parameterisation (dry diameters):
  real, parameter, dimension(SS_MAAR)::    &
      DP   = (/0.035, 0.075, 0.125, 0.195, 0.335, 0.51, 0.85 /) ! diameters Ddry
  real, parameter, dimension(SS_MAAR+1)::    &
      Dbin = (/0.02, 0.05, 0.1, 0.145, 0.25, 0.419, 0.6, 1.25 /) ! bins' borders Dd
    ! D_80%rh(/0.04, 0.10, 0.2, 0.290, 0.50, 0.838, 1.2, 2.5 /) ! bins' borders D80


!// Size bin borders (for dry R) for Monahan parameterisation
  real, parameter, dimension(SS_MONA+1) ::    &
        RLIM  = (/ 0.625, 1.0, 2.0, 3.0 /)  ! for Rdry
!extended        RLIM  = (/ 0.625, 1.0, 2.0, 3.0, 5.0, 7.0 /)  ! for Rdry
  real, parameter :: K1 = 0.7674, K2 = 3.079, K3 = 2.573e-11, K4 = -1.424
  real, parameter :: third = 1.0/3.0
  real :: lim1, lim2
  real, dimension(SS_MAAR) ::  dp2, dp4  
 !---------------------------------------------------- 

    n_to_mSS = PI * SSdens / 6.0  ! number to mass convertion
 
    log_dbin(:) = log10(Dbin(:))

!.. powers of diameter
     dp2(:) = DP(:)  * DP(:)
     dp3(:) = dp2(:) * DP(:)
     dp4(:) = dp3(:) * DP(:)

!//===== For Maartinsson et al.(2004) parameterisation =======

     a(1:3) = C1(1)*dp4(1:3)*MKM4   + C1(2)*dp3(1:3)  *MKM3        &
            + C1(3)*dp2(1:3)*MKM2   + C1(4)*DP(1:3)   *MKM + C1(5)
     a(4:5) = C2(1)*dp4(4:5)*MKM4   + C2(2)*dp3(4:5)  *MKM3        &
            + C2(3)*dp2(4:5)*MKM2   + C2(4)*DP(4:5)   *MKM + C2(5)
     a(6:7) = C3(1)*dp4(6:7)*MKM4   + C3(2)*dp3(6:7)  *MKM3        &
            + C3(3)*dp2(6:7)*MKM2   + C3(4)*DP(6:7)   *MKM + C3(5)

     b(1:3) = D1(1)*dp4(1:3)*MKM4   + D1(2)*dp3(1:3)  *MKM3        &
            + D1(3)*dp2(1:3)*MKM2   + D1(4)*DP(1:3)   *MKM + D1(5)
     b(4:5) = D2(1)*dp4(4:5)*MKM4   + D2(2)*dp3(4:5)  *MKM3        &
            + D2(3)*dp2(4:5)*MKM2   + D2(4)*DP(4:5)   *MKM + D2(5)
     b(6:7) = D3(1)*dp4(6:7)*MKM4   + D3(2)*dp3(6:7)  *MKM3        &
            + D3(3)*dp2(6:7)*MKM2   + D3(4)*DP(6:7)   *MKM + D3(5)

!//====== For  Monahan et al. (1986) parameterisation  =====

     rdry(1) = 0.8    ! Diameter at 80% ca. 3.1
     rdry(2) = 1.5    !                     6.3
     rdry(3) = 2.5    !                     10.6
 !.. can be extended 
!     rdry(4) = 4.0    !                     17
!     rdry(5) = 6.0    !                     26       

!.. Equilibrium wet radius from Gerber(1985) (Gong&Barrie [1997], JGR)

     do i = 1, SS_MONA
        !Gb radSS(i) = ( K1*rdry(i)**K2 /(K3 *rdry(i)**K4 -     &
        !Gb              log10(0.8))+rdry(i)**3) ** third
        !Gb lim1 = ( K1*RLIM(i+1)**K2 /(K3 *RLIM(i+1)**K4 -     &
        !Gb          log10(0.8))+RLIM(i+1)**3) ** third
        !Gb lim2 = ( K1*RLIM(i)**K2 /(K3 *RLIM(i)**K4 -         &
        !Gb          log10(0.8))+RLIM(i)**3) ** third

        !ds now use Gerber functions
!STbug        radSS(i) = umWetRad(rdry(i), 0.8, GbSeaSalt)
!STbug        lim1     = umWetRad(rlim(i+1), 0.8, GbSeaSalt)
!STbug        lim2     = umWetRad(rlim(i), 0.8, GbSeaSalt)
        radSS(i) = cmWetRad(rdry(i), 0.8, GbSeaSalt)
        lim1     = cmWetRad(rlim(i+1), 0.8, GbSeaSalt)
        lim2     = cmWetRad(rlim(i), 0.8, GbSeaSalt)
        Rrange(i) = lim1 - lim2       ! bin size intervals 

       if( DEBUG%SEASALT ) then
         ! WetRad takes radius in m
        write(*,"(a,i4,9g10.3)") "SSALT WETRAD ", i, radSS(i), lim1, lim2 !,  &
!          WetRad(rdry(i)*MKM, 0.8, GbSeaSalt)/MKM,&
!          WetRad(RLIM(i+1)*MKM, 0.8, GbSeaSalt)/MKM,&
!          WetRad(RLIM(i)*MKM, 0.8, GbSeaSalt)/MKM, &
!          umWetRad(rdry(i), 0.8, GbSeaSalt),&
!          umWetRad(RLIM(i+1), 0.8, GbSeaSalt),&
!          umWetRad(RLIM(i), 0.8, GbSeaSalt)
       end if
     end do

!.. Help parameter
     do i = 1, SS_MONA
          a1 = ( 0.380 - log10(radSS(i)) ) / 0.650
          a2 = 1.19 * exp(-a1*a1)

!st update /3.84e-6     temp_Monah(i) = 1.373 * radSS(i)**(-3) * Rrange(i) *
          temp_Monah(i) = 3.5755e5 * radSS(i)**(-3) * Rrange(i) *      &
                          ( 1.0 + 0.057 * radSS(i)**1.05 )* 10.0**a2
     end do

!// D_dry^3 -  for production of dry SS mass
     dSS3(:) =  ( 2.0 * rdry(:) )**3


   end subroutine init_seasalt
!>>--------------------------------------------------------------------------->>

 end module SeaSalt_mod

