! <Setup_1d_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Setup_1d_mod
!-----------------------------------------------------------------------!
! DESCRIPTION
! Generates arrays for 1-D column , for input to chemical solver. The output
! fields are stored in the ZchemData_mod module.
!-----------------------------------------------------------------------!

use AeroConstants_mod,only: AERO,NSAREA_DEF   ! for aerosol surface area
use AeroFunctions_mod,only: umWetRad,WetRad, pmSurfArea, cMolSpeed, UptakeRate
use AirEmis_mod,      only: airn, airlig   ! airborne NOx emissions
use Biogenics_mod,    only: SoilNOx
use Biogenics_mod,    only: EMIS_BioNat, EmisNat  
use ChemDims_mod,     only: NSPEC_SHL, NSPEC_ADV, NCHEMRATES, NEMIS_File
use ChemFields_mod,   only: SurfArea_um2cm3, xn_adv,xn_bgn,xn_shl, &
                               NSPEC_COL, NSPEC_BGN, xn_2d_bgn
use ChemFunctions_mod,only: S_RiemerN2O5
use ChemGroups_mod,   only:  chemgroups, PM10_GROUP
use ChemRates_mod,    only:  setChemrates ! rct, NRCT
use ChemSpecs_mod  !,           only:  SO4,C5H8,NO,NO2,SO2,CO,
use CheckStop_mod,    only:  CheckStop, StopAll,checkValidArray
use ColumnSource_mod, only: ColumnRate
use Config_module,    only:  &
   SKIP_RCT                     & ! kHet tests
  ,dt_advec                     & ! time-step
  ,IOU_INST                     & ! for OUTMISC
  ,MasterProc                   & 
  ,PPB, PT                      & ! PT-pressure at top
  ,USES                         & ! forest fires, hydrolysis, dergee_days etc.
  ,USE_OCEAN_NH3,USE_OCEAN_DMS,FOUND_OCEAN_DMS&
  ,VOLCANO_SR                   & ! Reduce Volcanic Emissions
  ,emis_inputlist               & ! Used in EEMEP
  ,KMAX_MID ,KMAX_BND, KCHEMTOP & ! Upper layer (k), upper level, and k for 1d fields
  ,EmisSplit_OUT
use Debug_module,     only: DebugCell, DEBUG_MASS, DEBUG !->DEBUG%SETUP_1DCHEM
use DerivedFields_mod,only: d_2d, f_2d
use EmisDef_mod,      only: gridrcemis, gridrcroadd, KEMISTOP,Emis_4D,N_Emis_4D,Found_Emis_4D&
                                , O_NH3, O_DMS, SplitEmisOut, EmisOut&
                                , NEmis_sources, Emis_source, Emis_source_2D
use EmisGet_mod,      only:  nrcemis, iqrc2itot, emis_nsplit, emis_masscorr
use ForestFire_mod,   only: Fire_rcemis, burning
use Functions_mod,    only:  Tpot_2_T
use GasParticleCoeffs_mod, only: DDspec
use GridValues_mod,   only:  xmd, GridArea_m2, & 
                             debug_proc, debug_li, debug_lj,&
                             A_mid,B_mid,gridwidth_m,dA,dB,&
                             i_fdom, j_fdom
use Io_Progs_mod,     only: datewrite !MASS
use Landuse_mod,      only: water_fraction, ice_landcover
use LocalVariables_mod,   only: Grid
use MassBudget_mod,   only: totem    ! sum of emissions
use MetFields_mod,    only: ps,sst
use MetFields_mod,    only: roa, th, q, t2_nwp, cc3dmax, &
                           zen, Idirect, Idiffuse,z_bnd,ws_10m
use Par_mod,          only: me,gi0,gi1,gj0,gj1,IRUNBEG,JRUNBEG
use PhysicalConstants_mod,only: ATWAIR, AVOG, PI, GRAV, T0
use Radiation_mod,    only: PARfrac, Wm2_uE
use SmallUtils_mod,   only: find_index
use Tabulations_mod,  only: tab_esat_Pa
use TimeDate_mod,     only: current_date, date
use Units_mod,        only: to_number_cm3 ! converts roa [kg/m3] to M [molec/cm3]
!!  to_number_cm3=0.001*AVOG/ATWAIR,& ! from density (roa, kg/m3) to molecules/cm3
use ZchemData_mod,    only: &
   xn_2d                &  ! concentration terms (molec/cm3)
  ,rcemis, deltaZcm     &  ! emission terms and lowest layer thickness
  ,rct, rcphot          &  ! chemical reaction rates
  ,rh, temp, tinv, itemp,pp      &  !
  ,M, o2, n2, h2o     &  ! Air concentrations
  ,cN2O5, cHO2, cO3, cHNO3 &  ! mol speeds, m/s
  ,cNO2, cNO3              &  ! mol speeds, m/s, kHetero tests
  ,gamN2O5                 &  ! kHetero test - for printout
  ,DpgNw, S_m2m3           &  ! for wet diameter and surf area
  ,aero_fom, aero_fbc, aero_fss, aero_fdust
 

implicit none
private
!-----------------------------------------------------------------------!

public :: setup_1d     ! Extracts results for i,j column from 3-D fields
public :: setup_rcemis ! Emissions in i,j column
public :: reset_3d     ! Exports final results from i,j column to 3-D fields
                       ! (and XNCOL outputs if asked for)

! Indices for the species defined in this routine. Only set if found
! Hard-coded for 2 specs just now. Could extend and allocate.
integer, private, parameter :: NROADDUST = 2
integer, private, parameter :: iROADF=1,  iROADC=2
integer, private, save :: inat_RDF,  inat_RDC
integer, private, save :: itot_RDF=-999,  itot_RDC=-999, itot_Rn222=-999



contains
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   subroutine setup_1d(i,j)

 !..   extracts data along one vertical column for input to chemical
 !     solver concentrations for chemistry......
 !

    integer, intent(in) :: i,j    ! coordinates of column
    character(len=9)  :: dtxt='setup_1d:'
    character(len=30)  :: fmt="(a32,i3,99es13.4)"  ! default format
    logical :: debug_flag
    logical, save :: first_call = .true.
   ! for surface area calcs: 
   ! had: SIA_F=1,PM_F=2,SS_F=3,DU_F=4,SS_C=5,DU_C=6,PM=7 ORIG=8,NSAREA=NSAREA_DEF
    integer, parameter :: NSAREA_CALC = NSAREA_DEF-1  ! Skip PM=7
    integer, save, dimension(NSAREA_CALC) :: aeroDDspec 
    real :: ugtmp, ugSIApm, ugDustF, ugSSaltF, ugDustC, ugSSaltC
    real, save ::  ugBCf=0.0, ugBCc=0.0 !not always present
    real :: ugSO4, ugNO3f, ugNO3c, ugRemF, ugRemC, ugpmF, ugpmC, rho
    logical :: is_finepm, is_ssalt, is_dust, is_sia
    logical, dimension(size(PM10_GROUP)), save  :: is_BC 
    ! AERO and DDspec still confusing
    real, dimension(NSAREA_DEF)  :: Ddry ! Dry diameter. 
    integer :: iw, ipm ! for wet rad
   ! if rates are wanted for d_2d output, we use these indices:
    integer, dimension(20), save :: d2index, id2rct 
    integer, save :: nd2d
    integer :: itmp

   ! local

    integer           :: k, n, ispec   ! loop variables
    real              :: qsat ! saturation water content
    integer, save :: nSKIP_RCT = 0
    integer, save :: iSIAgroup,iSSgroup, iDUgroup, iPMfgroup
    integer, save :: iBCfgroup,iBCcgroup
    logical, save, dimension(size(PM10_GROUP)) :: & ! arrays
       is_finepm_a,   &  ! 
       is_ssalt_a,    &  ! 
       is_dust_a,     &  ! 
       is_sia_a,      &  ! not needed?
       is_bc_a           ! will add soon

    debug_flag =  ( DEBUG%SETUP_1DCHEM .and. debug_proc .and.  &
      i==debug_li .and. j==debug_lj .and. current_date%seconds == 0 )

    if( first_call ) then

       if( debug_proc ) debug_flag = .true.  ! make sure we see 1st i,j combo

      !- find surrogates for Gerber area calcs, to get Dp,rho,sigma)
      ! Note, if a species is not found (eg seasalt) then iGerber
      ! is negative
      !QUERY: Should we then set USES%SURF_AREA False?

! for surface area calculations: Skips SIA. Find indices:
!    Skip PM=5 as mixed 

       aeroDDspec(AERO%PM_F)= find_index('PMf',DDspec(:)%name)
       aeroDDspec(AERO%SS_F)= find_index('SSf',DDspec(:)%name)
       aeroDDspec(AERO%DU_F)= find_index('DUf',DDspec(:)%name)
       aeroDDspec(AERO%SS_C)= find_index('SSc',DDspec(:)%name)
       aeroDDspec(AERO%DU_C)= find_index('DUc',DDspec(:)%name)

      ! see which groups we have:
       iSIAgroup = find_index('SIA',chemgroups(:)%name)
       iSSgroup  = find_index('SS',chemgroups(:)%name)
       iDUgroup  = find_index('Dust',chemgroups(:)%name,any_case=.true.)
       iPMfgroup = find_index('PMFINE',chemgroups(:)%name)
       iBCfgroup = find_index('ECFINE',chemgroups(:)%name)
       iBCcgroup = find_index('ECCOARSE',chemgroups(:)%name)

       if( MasterProc ) write(*,'(a,9(a,1x,i4))') dtxt//"Groups: ",&
          'SIA', iSIAgroup, 'SS', iSSgroup, 'DU', iDUgroup, &
          'BCf', iBCfgroup, 'BCc', iBCcgroup, 'PMf', iPMfgroup


       is_BC(:) = .false.
       if ( iBCfgroup > 0 ) then
          do ipm = 1, size( PM10_GROUP )
              ispec = PM10_GROUP(ipm)
              is_BC(ipm)  = ( find_index( ispec, chemgroups(iBCfgroup)%specs ) >0 )
              if( iBCcgroup > 0 ) then ! have coarse BC too
                 if( find_index( ispec, chemgroups(iBCcgroup)%specs ) >0) &
                    is_BC(ipm)  = .true.
              end if
              !if( MasterProc) write(*,*) dtxt//"is_BC ",&
              !    iBCcgroup, species(ispec)%name, is_BC(ipm)
          end do
       end if

      ! changes - SS and DU not essential now
       do ipm = 1, size( PM10_GROUP )
          ispec = PM10_GROUP(ipm)
          if ( iPMfgroup>0 ) is_finepm_a(ipm)  = &
                    ( find_index( ispec, chemgroups(iPMfgroup)%specs ) >0)
          if ( iSSgroup>0 ) is_ssalt_a(ipm)  = &
                    ( find_index( ispec, chemgroups(iSSgroup)%specs ) >0)
          if ( iDUgroup>0 ) is_dust_a(ipm)  = &
                    ( find_index( ispec, chemgroups(iDUgroup)%specs ) >0)
          if ( iSIAgroup>0 ) is_sia_a(ipm)  = &
                    ( find_index( ispec, chemgroups(iSIAgroup)%specs ) >0)
          if( MasterProc) then
            write(*,'(2a,5i5,5L2)') dtxt//"is_ PMf, SS, DU BC SIA? ",&
            species(ispec)%name, ipm, iPMfgroup,iSSgroup,iDUgroup,iBCfgroup,&
               is_finepm_a(ipm), is_ssalt_a(ipm), is_dust_a(ipm),is_BC(ipm), is_sia_a(iPM)
          end if
       enddo

       do n = 1, size(SKIP_RCT)
          if ( SKIP_RCT(n) > 0 ) nSKIP_RCT = nSKIP_RCT  + 1
       end do
       if( MasterProc ) write(*,"(a,10i4)") &
          dtxt//"SKIP_RCT:", SKIP_RCT(1:nSKIP_RCT)

    end if ! first_call

    if( debug_flag ) write(*,*) dtxt//"=DBG=======  ", first_call, me


    do k = KCHEMTOP, KMAX_MID

  !- to_number_cm3 - to scale from  density (roa, kg/m3) to  molecules/cm3
  ! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)

       M(k) = roa(i,j,k,1) * to_number_cm3  ! molecules air/cm3

       h2o(k) = max( 1.e-5*M(k), &
                     q(i,j,k,1)*M(k)*ATWAIR/18.0)

      ! nb. max function for h2o used as some NWP numerics can give 
      ! negative negative H2O....   :-(

       pp(k) = A_mid(k) + B_mid(k)*ps(i,j,1)

       temp(k) = th(i,j,k,1)* Tpot_2_T( pp(k) )

       itemp(k) = nint( temp(k) -1.E-9) ! the "-1.E-9" is put in order to
               ! avoid possible different roundings on different machines.

       qsat  = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
       rh(k) = min( q(i,j,k,1)/qsat , 1.0)
       rh(k) = max( rh(k) , 0.001)

        ! 1)/ Short-lived species - no need to convert units

         do n = 1, NSPEC_SHL
               xn_2d(n,k) = max(0.0,xn_shl(n,i,j,k))
         end do ! ispec

        ! 2)/ Advected species
        do n = 1, NSPEC_ADV
              ispec = NSPEC_SHL + n
              xn_2d(ispec,k) = max(0.0,xn_adv(n,i,j,k)*M(k))
        end do ! ispec

        ! 3)/ Background species (with units in mix. ratio)
        do n = 1, NSPEC_BGN
              xn_2d_bgn(n,k) = max(0.0,xn_bgn(n,i,j,k)*M(k))
        end do ! ispec

      ! Surf Area
        if ( USES%SURF_AREA ) then ! GERBER

            S_m2m3(:,k) = 0.0  !! Allow max 6000 um2/cm3

           !ispec=NO3_c ! CRUDE HARD CODE for now, but NO3 is special

           ugSO4       = 0.0
           ugRemF      = 0.0
           ugRemC      = 0.0
           ugpmF       = 0.0
           ugpmC       = 0.0
           ugNO3f      = 0.0 !  0.27*xn_2d(i,k)*species(i)%molwt*1.0e12/AVOG
           ugNO3c      = 0.0 !  0.27*xn_2d(i,k)*species(i)%molwt*1.0e12/AVOG
           ugSIApm     = 0.0 
           ugBCf       = 0.0
           ugBCc       = 0.0

           ugDustF  = 0.0
           ugSSaltF = 0.0
           ugDustC  = 0.0
           ugSSaltC = 0.0

           do ipm = 1, size( PM10_GROUP )
             ispec = PM10_GROUP(ipm)

             ugtmp  = xn_2d(ispec,k)*species(ispec)%molwt*1.0e12/AVOG
             is_finepm = is_finepm_a(ipm)
             is_ssalt  = is_ssalt_a(ipm)
             is_dust   = is_dust_a(ipm)
             is_sia    = is_sia_a(ipm)
             if( is_finepm ) then
               ugpmF  = ugpmF   + ugtmp
               if(is_ssalt) ugSSaltF = ugSSaltF  +  ugtmp
               if(is_dust ) ugDustF  = ugDustF   +  ugtmp
               if( is_sia ) & !BUG iSIAgroup >0 ) &
                    ugSIApm  = ugSIApm + ugtmp
               if(is_BC(ipm)) ugBCf = ugBCf  +  ugtmp
             else
                ugpmC  = ugpmC   + ugtmp
               if(is_ssalt) ugSSaltC = ugSSaltC  +  ugtmp
               if(is_dust ) ugDustC  = ugDustC   +  ugtmp
               if(is_BC(ipm)) ugBCc = ugBCc  +  ugtmp
             end if
!          ugRemF = ugpmf - ugSIApm -ugSSaltF -ugDustF
            !if ( MasterProc .and. k == KMAX_MID) write(*,"(a,i5,3L3,9es10.3)") &
            !if ( debug_flag .and. k == KMAX_MID) write(*,"(a,i5,3L3,9es10.3)") &
            !   'AREACHECK0: '//trim(species(ispec)%name), ispec, &
            !      is_finepm, is_ssalt, is_dust, ugtmp, ugpmF,ugSIApm, ugSSaltF, ugDustF
 

             if (  species(ispec)%name == 'SO4' ) then
               ugSO4 = ugSO4  +  ugtmp
             else if ( index( species(ispec)%name, 'NO3_f' )>0) then
               ugNO3f= ugNO3f +  ugtmp
             else if ( index( species(ispec)%name, 'NO3_c' )>0) then
               ugNO3c= ugNO3c +  ugtmp
             end if

            ! Handle surface area here:
 
             if( debug_flag .and. k==KMAX_MID ) then
               write(*, fmt) dtxt//"UGSIA:"//trim(species(ispec)%name), &
                     k, ugtmp, ugpmF, ugpmC,&
                     ugSO4,ugNO3f,ugNO3c,ugBCf,ugBCc
               !write(*,'(a,4L2,es12.3)') dtxt//trim(species(ispec)%name), &
               !      is_finepm, is_ssalt,is_dust,is_BC(ipm), ugtmp
             end if
           end do ! ipm



         ! FRACTIONS used for N2O5 hydrolysis
         ! We use mass fractions, since we anyway don't have MW for OM, dust,
         !  ugRemF will include OM, EC, PPM, Treat as OM 

          ugRemF = ugpmf - ugSIApm -ugSSaltF -ugDustF

          ugRemF = ugRemF  -ugBCf  ! kHet BC added 24/10/2015

          aero_fss(k)     = ugSSaltF/ugpmF
          aero_fdust(k)   = ugDustF/ugpmF
          aero_fbc(k)     = ugBCf/ugpmF
          aero_fom(k)     = max(0.0, ugRemF)/ugpmF

          ! GERBER equations for wet radius
          ! New approach to Gerber. We use masses above. 
          ! STILL CRUDE :-( NOT SURE THIS IS BEST APPROACH but
          ! hopefully works 'okay', maybe??!
          ! NOTE: There was a bug anyway in Gerber Table, but we recode
          ! now to avoid need for Inddry index

          do iw = 1, NSAREA_CALC ! Skips PM=7 which is sum

            ipm= aeroDDspec(iw)

            if ( ipm < 1 ) then
               if ( debug_flag .and. k == KMAX_MID) &
                   write(*,"(a,i0)") dtxt//'CHECK:iPMNEG:',ipm
               CYCLE ! Component missing !
            end if

            Ddry(iw) =  DDspec(ipm)%DpgN  ! (m)

            if ( DDspec(ipm)%Gb > 0 ) then
              DpgNw(iw,k)  = 2*WetRad( 0.5*Ddry(iw), rh(k), DDspec(ipm)%Gb ) 
            else
              DpgNw(iw,k)  = Ddry(iw) ! for dust, we keep dry
            end if

            if ( iw == AERO%PM_F )  ugtmp = ugpmF
            if ( iw == AERO%SS_F )  ugtmp = ugSSaltF
            if ( iw == AERO%DU_F )  ugtmp = ugDustF
            if ( iw == AERO%SS_C )  ugtmp = ugSSaltC
            if ( iw == AERO%DU_C )  ugtmp = ugDustC

            S_m2m3(iw,k) = pmSurfArea(ugtmp,Dp=Ddry(iw), Dpw=DpgNw(iw,k),  &
                                  rho_kgm3=DDspec(ipm)%rho_p )

            if ( debug_flag .and. k == KMAX_MID) write(*,"(a,3i5,4es10.3)") &
              "AREACHECK: "//trim(DDspec(ipm)%name), iw, ipm, &
             DDspec(ipm)%Gb, &
               Ddry(iw),DpgNw(iw,k),ugtmp,S_m2m3(iw,k)

            S_m2m3(iw,k) = min( S_m2m3(iw,k), 6.0e-3)  !! Allow max 6000 um2/cm3

          end do ! iw

          ! Add all areas (misses some BC_c etc)
          iw=NSAREA_DEF
          S_m2m3(iw,k) = S_m2m3(AERO%PM_F,k) + S_m2m3(AERO%SS_C,k) + S_m2m3(AERO%DU_C,k)
          if ( debug_flag .and. k == KMAX_MID) write(*,"(a,i5,4es10.3)") &
             "AREACHECK: "//"SUM", iw, ugpmf,ugSSaltC,ugDustC,S_m2m3(iw,k)


          if( DEBUG%SETUP_1DCHEM ) then ! extra checks 
             if( aero_fom(k) > 1.0 .or. ugRemF < -1.0e-9 ) then
                print "(a,i4,99es12.3)", dtxt//"AERO-F ", k, &
                  aero_fom(k), ugRemF,ugpmF, ugSSaltF, ugDustF, ugBCf
                call CheckStop(dtxt//"AERO-F problem " )
              end if
          end if

          if ( DebugCell .and. k==KMAX_MID ) then
            write(*,fmt)  dtxt//" SAREAugPM in  ", k,  rh(k), temp(k), &
            ugpmf, ugSSaltC, ugDustC
            do iw = 1, NSAREA_CALC
             write(*,fmt) dtxt//"GERB ugDU (S um2/cm3)  ", iw, Ddry(iw), 1.0e6*S_m2m3(iw,k)
            end do
          end if

           ! m2/m3 -> um2/cm3 = 1.0e6, only for output to netcdf
           if( k == KMAX_MID ) then 
              do iw = 1, NSAREA_DEF  !J2018 bug:CALC
                SurfArea_um2cm3(iw,i,j) = 1.0e6* S_m2m3(iw,k)
              end do
           end if

          end if ! GERBER

      ! End Surf Area
    
   end do ! k

! Check that concentrations are not "contaminated" with NaN
   if ( isnan( xn_2d(NO2,KMAX_MID) ) ) then
      print *, "NANAN ", trim(species(NO2)%name)
      call CheckStop( "Detected non numerical concentrations (NaN)")
   end if

   o2(:) = 0.21 *M(:)
   n2(:) = M(:) - o2(:)
!   o2(:) = 0.2095 *M(:) ! more exact, but prefer o3+n2 to add to 100%
!   n2(:) = 0.7808 *M(:)
   tinv(:) = 1./temp(:)


   cN2O5(:) = cMolSpeed(temp(:),108.0)
   cHNO3(:) = cMolSpeed(temp(:), 63.0)
   cHO2(:)  = cMolSpeed(temp(:), 33.0)
   cO3(:)   = cMolSpeed(temp(:), 48.0)
   cNO3(:)  = cMolSpeed(temp(:), 62.0)
   cNO2(:)  = cMolSpeed(temp(:), 46.0)

  ! 5 ) Rates  (!!!!!!!!!! NEEDS TO BE AFTER RH, XN, etc. !!!!!!!!!!)


   !====================
   call setChemRates()
   !====================

   if( DEBUG%SETUP_1DCHEM ) then ! extra checks 
     call checkValidArray(rcphot(:,KMAX_MID),dtxt//'arrayCheck RCPHOT ')
     call checkValidArray(rct(:,KMAX_MID),dtxt//'arrayCheck RCT ')
     if ( DebugCell .and. current_date%hour == 12 ) then
       print *, "DBG: RCPHOT SIZE ",  size(rcphot,dim=1)
       do n = 1, size(rcphot,dim=1)
         if( rcphot(n,KMAX_MID)<0.0) then
             print '(a,i3,es12.3)', 'RCPHOT:', n, rcphot(n,KMAX_MID)
         end if
       end do
     end if 
   end if

  ! For sensitivity tests
   do  n = 1, nSKIP_RCT ! can be zero
     !if( SKIP_RCT(n) == 720 ) then
     !   rct(72:75,:) = 0.0    ! HNO3 + SS, DU
     !else if ( SKIP_RCT(n) == 770 ) then
     !   rct(77:78,:) = 0.0    ! DU_f+DU_c, O3
     !else 
     !   if ( first_call ) call CheckStop( SKIP_RCT(n) > 100,&
     !                                     dtxt//"SKIP_RCT too big")
        rct(SKIP_RCT(n),:) = 0.0
     !end if
   end do


   if ( first_call ) then
     call CheckStop( any(isnan(rct(:,:))), dtxt//"RCT NAN'd")
     call CheckStop( any(rct(:,:) < 0.0 ), dtxt//"RCT NEG'd")
     nd2d = 0
     do itmp = 1, size(f_2d)
           if ( f_2d(itmp)%subclass == 'rct' ) then
            !check if index of rate constant (config) possible
             if ( f_2d(itmp)%index > NCHEMRATES ) then
                if(MasterProc) write(*,*) 'RCT NOT AVAILABLE!', itmp, &
                  f_2d(itmp)%index, NCHEMRATES
                cycle
             end if
             nd2d =  nd2d  + 1
             call CheckStop(nd2d>size(id2rct),dtxt//"Need bigger id2rct array")
             d2index(nd2d)= itmp
             id2rct(nd2d) = f_2d(itmp)%index  !index of rate constant (config)
             if(MasterProc) write(*,*) 'RCTFOUND', itmp, nd2d, id2rct(nd2d)
           end if
     end do

     first_call = .false.
   end if ! first_call

   do itmp = 1, nd2d
       d_2d(d2index(itmp),i,j,IOU_INST) =  rct(id2rct(itmp),KMAX_MID)
       if( debug_flag ) then
          write(*,"(a,6i5,es12.3)") dtxt//"OUTRCT "//&
           trim(f_2d(d2index(itmp))%name), me,i,j,itmp,&
           d2index(itmp),id2rct(itmp), d_2d(d2index(itmp),i,j,IOU_INST)
       end if
   end do


   end subroutine setup_1d
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine setup_rcemis(i,j)
!-------------------------------------------------------------------
!    DESCRIPTION:
!    Extracts emissions in column from gridrcemis, for input to chemistry
!    routines. Results in "rcemis" array with unts: molecule/cm3/s
!-------------------------------------------------------------------
  !-- arguments
  integer, intent(in) ::  i,j     ! coordinates of column

  !  local
  integer :: iqrc, k, itot, iem ,f
  real    :: Kw,fac, eland   ! for Pb210  - emissions from land

  integer :: i_Emis_4D,n
  logical, save     :: first_call = .true. 
  character(len=13) :: dtxt="setup_rcemis:"
  real :: SC_DMS,SC_DMS_m23,SC_DMS_msqrt,SST_C,invDeltaZfac
  integer,save ::IC_NH3
  real :: ehlpcom0

  ehlpcom0 = GRAV* 0.001*AVOG!0.001 = kg_to_g / m3_to_cm3              


  if(first_call)then
    inat_RDF = find_index( "Dust_ROAD_f", EMIS_BioNat(:) , any_case=.true. )
    inat_RDC = find_index( "Dust_ROAD_c", EMIS_BioNat(:) , any_case=.true. )
    itot_RDF = find_index( "Dust_ROAD_f", species(:)%name  , any_case=.true.  )
    itot_RDC = find_index( "Dust_ROAD_c", species(:)%name  , any_case=.true.  )
    itot_Rn222=find_index( "RN222", species(:)%name    )
    IC_NH3=find_index( "NH3", species(:)%name    )
    first_call = .false.
  end if 

! initilize ! initilize ! initilize ! initilize
  rcemis(:,:)=0.
! initilize ! initilize ! initilize ! initilize

  forall(k=KEMISTOP:KMAX_MID,iqrc=1:NRCEMIS) &
    rcemis(iqrc2itot(iqrc),k) = gridrcemis(iqrc,k,i,j)&
                                  *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))


   ! add emissions from new format which are "pure", i.e. not defined as one of
   ! the splitted or sector species
    do n = 1, NEmis_sources      
       itot = Emis_source(n)%species_ix
       if(itot>0 .and. Emis_source(n)%sector<=0)then
          k = Emis_source(n)%injection_k
          rcemis(itot,k) = rcemis(itot,k)   &
               + Emis_source_2D(i,j,n)*ehlpcom0    &
               /species(itot)%molwt&
               *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))
       endif
    enddo



  ! Volcanic emissions (SO2 and ASH),
  ! and Contribution from Emergeny scenarios
  if(VOLCANO_SR)then
    rcemis(:,:)=rcemis(:,:)+ColumnRate(i,j,REDUCE_VOLCANO=0.85)
  else
    rcemis(:,:)=rcemis(:,:)+ColumnRate(i,j)
  end if

  ! lightning and aircraft ... Aerial NOx emissions if required:
  if(USES%LIGHTNING_EMIS)then
    do k=KCHEMTOP, KMAX_MID
      rcemis(NO ,k) = rcemis(NO ,k) + 0.95 * airlig(k,i,j)
      rcemis(NO2,k) = rcemis(NO2,k) + 0.05 * airlig(k,i,j)
    end do
  end if
  if(USES%AIRCRAFT_EMIS) then
    do k=KCHEMTOP, KMAX_MID
      rcemis(NO ,k) = rcemis(NO ,k) + 0.95 * airn(k,i,j)
      rcemis(NO2,k) = rcemis(NO2,k) + 0.05 * airn(k,i,j)
    end do
  end if ! AIRCRAFT NOX
  if(DEBUG%SETUP_1DCHEM.and.debug_proc.and.i==debug_li.and.j==debug_lj)&
    write(*,"(a,2L2,10es10.3)") &
      dtxt//"AIRNOX ", USES%LIGHTNING_EMIS, USES%AIRCRAFT_EMIS, &
      airn(KMAX_MID,i,j),airlig(KMAX_MID,i,j), &
      airn(KCHEMTOP,i,j),airlig(KCHEMTOP,i,j)

  ! Road dust
  if(USES%ROADDUST.and.itot_RDF>0) then  ! Hard-code indices for now
    rcemis(itot_RDF,KMAX_MID) = gridrcroadd(1,i,j)
    rcemis(itot_RDC,KMAX_MID) = gridrcroadd(2,i,j)
   end if
   
  ! Forest fires
  if(USES%FOREST_FIRES) then
    if(burning(i,j))call Fire_rcemis(i,j)
  end if

  ! Soil NOx
  if(USES%GLOBAL_SOILNOX)then
    rcemis(NO,KMAX_MID)=rcemis(NO,KMAX_MID)+SoilNOx(i,j)
  end if

  ! Emissions from GEIA, was use for Aerocom NO3 experiment
  if(USE_OCEAN_NH3)then
     !keep separated from snapemis/rcemis in order to be able to include more
     ! advanced processes
     k=KMAX_MID
     !convert from kg/m2/s into molecules/cm3/s . 
     !kg->g=1000 , /m3->/cm3=1e-6 , 1000*1e-6=0.001
     
     rcemis(O_NH3%index,k)=rcemis(O_NH3%index,k)+ &
       O_NH3%emis(i,j)*0.001*AVOG/species(IC_NH3)%molwt&
          *(GRAV*roa(i,j,k,1))/(dA(k)+dB(k)*ps(i,j,1))
    !make map in mg/m2
     O_NH3%map(i,j)=O_NH3%emis(i,j)*dt_advec*1.0E6!kg->mg = 1.0E6

  end if

  ! Ocean DMS -> SO2
  if(FOUND_OCEAN_DMS)then!temporarily always compute budgets if file found
     !keep separated from snapemis/rcemis in order to be able to include more
     ! advanced processes
     k=KMAX_MID
     !convert from mol/cm3 into molecules/cm3/s . 
     !Kw in cm/hour after Leonor Tarrason (1995), after Liss and Merlivat (1986)
     !assumes 20 degrees C
    !NB: misprint in gbc1385.pdf!
     SST_C=max(0.0,min(30.0,(SST(i,j,1)-T0))) !the formula uses degrees C
     SC_DMS=2674 -147.12*SST_C + 3.726*SST_C*SST_C - 0.038 * SST_C*SST_C*SST_C
     SC_DMS=SC_DMS/600.0
     SC_DMS_m23=SC_DMS**(-2.0/3.0)

     SC_DMS_msqrt=sqrt(1./SC_DMS)
     if(ws_10m(i,j,1)<=3.6)then
        Kw=0.17*ws_10m(i,j,1)*SC_DMS_m23
     elseif(ws_10m(i,j,1)<=13.0)then
        Kw=0.17*ws_10m(i,j,1)*SC_DMS_m23+2.68*(ws_10m(i,j,1)-3.6)*SC_DMS_msqrt
     else
        Kw=0.17*ws_10m(i,j,1)*SC_DMS_m23+ &
           2.68*(ws_10m(i,j,1)-3.6)*SC_DMS_msqrt+ &
           3.05*(ws_10m(i,j,1)-13)*SC_DMS_msqrt
     end if
     Kw=Kw/3600!cm/hour -> cm/s

!66% of DMS turns into SO2, Leonor Tarrason (1995)
     if(USE_OCEAN_DMS)then
      rcemis(O_DMS%index,k)=rcemis(O_DMS%index,k)+ &
      0.66*O_DMS%emis(i,j)*Kw*0.01*GRAV*roa(i,j,k,1)/ &
                           (dA(k)+dB(k)*ps(i,j,1)) *AVOG 
     end if
     !in g . Multiply by dz(in cm)  * dx*dx (in cm2) * ...
     !... molwgt(SO2) /AVOG . (dz/AVOG just removed from above) 
     !g->Gg = 1.0E-9
     O_DMS%sum_month = O_DMS%sum_month+0.66*O_DMS%emis(i,j)*Kw *1.e4*xmd(i,j)*&
          gridwidth_m*gridwidth_m*dt_advec *64.0*1.0E-9

!make map in mg/m2,    g->mg = 1.0E3  ; /cm2 -> /m2 =1e4
     O_DMS%map(i,j)=O_DMS%emis(i,j)*Kw *1.e4*62.13*1.0E3  

  end if



  if(Found_Emis_4D>0)then
     do i_Emis_4D=1,N_Emis_4D
        if(emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D)=='NOTSET')exit
        n=find_index(emis_inputlist(Found_Emis_4D)%pollemepName(i_Emis_4D),&
                      species(:)%name)
        if(n>0)then
           fac=1.0/1000000.0/3600.0 !convert from Bq/m3/hour into Bq/cm3/s
           do k=KCHEMTOP, KMAX_MID
              rcemis(n,k)=rcemis(n,k)+Emis_4D(i,j,k,i_Emis_4D)*fac
           end do
        end if
     end do
  end if


  do k=KCHEMTOP, KMAX_MID

     deltaZcm(k) = 100*( z_bnd(i,j,k)-z_bnd(i,j,k+1) )

    if(DEBUG%SETUP_1DCHEM.and.debug_proc.and.i==debug_li.and.j==debug_lj) then
       write(*,"(a,i3,9es12.4)") "1DCHEM DZ",  k, deltaZcm(k), &
        100*(dA(k)+dB(k)*ps(i,j,1))/(GRAV*roa(i,j,k,1))! , &
!        100*(dA(k)+dB(k)*ps(i,j,1))/(GRAV*M(k)*ATWAIR)!, &
!        M(k)*ATWAIR/roa(i,j,k,1)
!       M(k) = roa(i,j,k,1) * to_number_cm3  ! molecules air/cm3
    end if
  end do

  if(EmisSplit_OUT)then
    !put all added emissions in EmisSplit_OUT, also natural emissions
    SplitEmisOut(i,j,:)=0.0
    do k=KCHEMTOP, KMAX_MID
      do iqrc=1,nrcemis
        !give unit mg/m2/s dt_advec multiplied in Derived_mod 
        itot=iqrc2itot(iqrc)
        SplitEmisOut(i,j,iqrc) = SplitEmisOut(i,j,iqrc)&
          +rcemis(itot,k)*species(itot)%molwt &
          *(dA(k)+dB(k)*ps(i,j,1))/(GRAV*M(k)*ATWAIR)
      end do
    end do
  end if

  iqrc = 0 
  do iem = 1, NEMIS_FILE   
     EmisOut(i,j,iem) = 0.0
     do f = 1,emis_nsplit(iem)
        iqrc = iqrc + 1
        itot = iqrc2itot(iqrc)
        do k=KCHEMTOP, KMAX_MID
           EmisOut(i,j,iem) = EmisOut(i,j,iem) + rcemis(itot,k)/emis_masscorr(iqrc) &
          *(dA(k)+dB(k)*ps(i,j,1))/(GRAV*M(k)*ATWAIR)
        enddo
     enddo
  enddo


  ! Soil Rn222 emissions from non-ice covered land, + water
  ! at rate of 1 atom/cm2/s
   eland = 1.0 - water_fraction(i,j) - ice_landcover(i,j)


  ! z_bnd is in m, not cm, so need to divide by 100.
  if(itot_Rn222>0) &
    rcemis(itot_Rn222,KMAX_MID) = (0.00182*water_fraction(i,j)+eland)&
      /deltaZcm(KMAX_MID) 

end subroutine setup_rcemis
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine reset_3d(i,j)
  integer, intent(in) :: i,j
  integer :: k, n, ispec, id    ! loop variables
 !! XNCOL testing -- sets d_2d for column data from molec/cm3 concs.
 ! if variables are wanted for d_2d output (via USET), we use these indices:
  character(len=*),parameter :: dtxt='reset3dxncol:'
  character(len=20) :: specname
  integer, dimension(20), save :: d2index, id2col
  logical, save :: first_call = .true.
  integer, save :: nd2d


  do k = KCHEMTOP, KMAX_MID
    ! 1)/ Short-lived species - no need to scale with M
    do n = 1, NSPEC_SHL
      xn_shl(n,i,j,k) = xn_2d(n,k)
    end do ! ispec

    ! 2)/ Advected species
    do n = 1, NSPEC_ADV
      ispec = NSPEC_SHL + n
      xn_adv(n,i,j,k) = xn_2d(ispec,k)/M(k)
    end do ! ispec
  end do ! k

!!======================================================================
!! If column totals are wanted, we can do those here also since xn_2d are
!! in molec/cm3, and we want molec/cm2:

   if ( first_call ) then

     nd2d = 0
     do id = 1, size(f_2d)

           if ( f_2d(id)%subclass == 'xncol' ) then
             nd2d =  nd2d  + 1
             call CheckStop( nd2d > size(id2col), &
                 dtxt//"Need bigger id2col array" )
             specname = trim(f_2d(id)%name(7:))  ! Strip XNCOL_
             ispec = find_index( specname, species(:)%name )
             call CheckStop(ispec < 1, dtxt//"XNCOL not found"//specname )
             d2index(nd2d)= id
             id2col(nd2d) = ispec
             if(MasterProc) write(*,*) 'USET XNCOL FOUND', id, ispec, &
                 trim(specname),nd2d, id2col(nd2d)
           end if
     end do
     first_call = .false.
   end if

   do id = 1, nd2d
      ispec = id2col(id)
      d_2d(d2index(id),i,j,IOU_INST) = dot_product(xn_2d(ispec,:),deltaZcm(:))
      if(DEBUG%SETUP_1DCHEM.and.debug_proc.and. &
          i==debug_li.and.j==debug_lj) then

        write(*,"(a,6i5,a,9es12.3)") dtxt//"OUTXNCOL "//&
          trim(f_2d(d2index(id))%name), me, i,j,id,d2index(id),ispec, &
          trim(species(ispec)%name), xn_2d(ispec,20), deltaZcm(20), &
           d_2d(d2index(id),i,j,IOU_INST)
    end if
   end do


end subroutine reset_3d
!---------------------------------------------------------------------------
endmodule Setup_1d_mod
!_____________________________________________________________________________!
