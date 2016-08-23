! <MosaicOutputs_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

module MosaicOutputs_ml
 use AOTx_ml,          only : Calc_AOTx, Calc_POD
 use AOTx_ml,  only: O3cl, VEGO3_OUTPUTS
 use CheckStop_ml,  only: CheckStop, StopAll
 use ChemChemicals_ml, only : species
 use ChemSpecs_shl_ml, only : NSPEC_SHL
 use ChemSpecs_adv_ml, only : NSPEC_ADV
 use ChemGroups_ml,       only : DDEP_OXNGROUP,DDEP_SOXGROUP, &
                                 DDEP_RDNGROUP, NMAX_DDEP
 use DerivedFields_ml, only : f_2d, d_2d
 use EcoSystem_ml,     only : NDEF_ECOSYSTEMS, DEF_ECOSYSTEMS, &
    EcoSystemFrac, FULL_ECOGRID, FULL_LCGRID, Is_EcoSystem
 use Io_Progs_ml,      only : datewrite
 use LandDefs_ml,      only : LandDefs, LandType, &
         Check_LandCoverPresent ! e.g. "CF"
 use Landuse_ml,       only : LandCover ! for POD
 use LocalVariables_ml,   only : Sub, Grid
 use MetFields_ml
 use ModelConstants_ml, only : MasterProc, DEBUG => DEBUG_MOSAICS,&
   atwS, atwN, &
   NLANDUSEMAX, IOU_INST, &
   SOX_INDEX, OXN_INDEX, RDN_INDEX ! indices for dep groups

 use OwnDataTypes_ml,  only: Deriv, print_deriv_type, &
       TXTLEN_DERIV, TXTLEN_SHORT, typ_s5i
 use SmallUtils_ml, only: find_index
 use TimeDate_ml, only : current_date, effectivdaynumber
 use Wesely_ml, only : NDRYDEP_CALC
 implicit none
 private

 public ::  Init_MosaicMMC
 public ::  Add_MosaicMetConcs
 public ::  Add_MosaicVEGO3
 public ::  Add_MosaicDDEP
 public ::  Add_NewMosaics

 public :: Add_MosaicOutput
 public :: find_MosaicLC

 
 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO

 integer, public, save :: MMC_RH, MMC_CANO3, MMC_VPD, MMC_FST, &
     MMC_USTAR, MMC_INVL, MMC_GSTO, MMC_EVAP
 character(len=30),private, save :: errmsg = "ok"


 ! Mosaic-specific outputs, e.g. VG_CF_HNO3 or Rns_GR_NH3
  integer, public, save :: nMosaic = 0
  integer, public, parameter :: MAX_MOSAIC_OUTPUTS=100
  logical, private, parameter :: T=.true., F=.false.

  type(Deriv), public, &
     dimension( MAX_MOSAIC_OUTPUTS ), save :: MosaicOutput


 contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     subroutine Init_MosaicMMC(MOSAIC_METCONCS)
        character(len=*), dimension(:), intent(in) :: MOSAIC_METCONCS

    ! Set indices of mosaic metconc params for later use. Will be zero if 
    ! not found, but that's okay I hope...
      MMC_RH    = find_index("RH",MOSAIC_METCONCS)
      MMC_CANO3 = find_index("CanopyO3",MOSAIC_METCONCS)
      MMC_VPD   = find_index("VPD",MOSAIC_METCONCS)
      MMC_FST   = find_index("FstO3",MOSAIC_METCONCS)
      MMC_USTAR = find_index("USTAR",MOSAIC_METCONCS)
      MMC_INVL  = find_index("INVL",MOSAIC_METCONCS)
      MMC_GSTO  = find_index("GSTO",MOSAIC_METCONCS)
      MMC_EVAP  = find_index("EVAP",MOSAIC_METCONCS)

     end subroutine Init_MosaicMMC

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     subroutine Add_MosaicMetConcs(MOSAIC_METCONCS,MET_LCS,iotyp, nMET)
        character(len=*), dimension(:), intent(in) :: MOSAIC_METCONCS
        character(len=*), dimension(:), intent(in) :: MET_LCS
        integer, intent(in)  :: iotyp
        integer, intent(out) :: nMET
        integer :: ilab, n, iLC
        character(len=TXTLEN_DERIV) :: name


      !------------- Met data for d_2d -------------------------
      ! We find the various combinations of met and ecosystem,
      ! adding them to the derived-type array Mosaic_Met (e.g. => Met_CF)

      nMET = 0
      do ilab = 1, size(MOSAIC_METCONCS)
        MET_LC: do n = 1, size(MET_LCS)

          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "MET_LCS", n, MET_LCS, (ilab == 1))
          if ( iLC < 0 ) cycle  MET_LC
          if ( MOSAIC_METCONCS(ilab)(1:6) == "Canopy" .or. & 
               MOSAIC_METCONCS(ilab)(1:5) == "FstO3" ) &
                 LandType(iLC)%flux_wanted  = .true.  ! Canopy calc in StoFlux
          !-------------End of Check if LC present in this array ------!

          nMET = nMET + 1
          name = trim ( MOSAIC_METCONCS(ilab) ) // "_"  // trim( MET_LCS(n) )

          nMosaic = nMosaic + 1
          call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, nMET" )
          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d, scale, avg? Inst Yr Mn Day
           MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", "METCONC", MET_LCS(n), MOSAIC_METCONCS(ilab), &
                ilab, -99,F , 1.0,  T,  iotyp )    !FEBN2011 fix later

          if( MOSAIC_METCONCS(ilab)(1:5) == "USTAR" )  then
              MosaicOutput(nMosaic)%unit  =   "m/s"
          else if( MOSAIC_METCONCS(ilab)(1:4) == "INVL" )  then
              MosaicOutput(nMosaic)%unit  =   "m"
          else if( MOSAIC_METCONCS(ilab)(1:8) == "CanopyO3" )  then
              MosaicOutput(nMosaic)%unit  =   "ppb"
          else if( MOSAIC_METCONCS(ilab)(1:5) == "FstO3" )  then
              MosaicOutput(nMosaic)%unit  =   "mmole/m2" ! accumulated
          else if( MOSAIC_METCONCS(ilab)(1:4) == "EVAP" )  then
              MosaicOutput(nMosaic)%avg       =  .false. ! accumulate
              MosaicOutput(nMosaic)%unit      =  "mm"
              MosaicOutput(nMosaic)%dt_scale  =  .true.
          end if

          if(DEBUG .and. MasterProc) call print_deriv_type(MosaicOutput(nMosaic))
        end do MET_LC !n
      end do ! ilab

 end subroutine Add_MosaicMetConcs
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine Add_NewMosaics(Mc,nMc)
        type(typ_s5i), dimension(:), intent(in) :: Mc ! eg VG
        integer, intent(out) :: nMc
        integer :: n, itot, iLC, iadv
        character(len=TXTLEN_DERIV) :: name
        character(len=TXTLEN_SHORT) :: typ, poll, lctxt

      nMc = 0
      MC_LOOP: do n = 1, size( Mc(:)%txt1 )

          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "MMC-VG", n, Mc(:)%txt4 , &
                   write_condition=.true. )
          if ( iLC < 0 ) cycle  MC_LOOP
          !-------------End of Check if LC present in this array ------!
          nMc = nMc + 1
          typ   = Mc(n)%txt2
          poll  = Mc(n)%txt3
          lctxt = Mc(n)%txt4

         ! name = VG_O3_GRID?
          name = trim( Mc(n)%txt2 ) // "_" //  & ! VG
                 trim( Mc(n)%txt3 ) // "_" //  & ! O3, HNO3, ...
                 trim( Mc(n)%txt4 )              ! Grid, SNL,..

          itot = find_index( poll,  species(:)%name )
          iadv = itot - NSPEC_SHL
          if( iadv < 1 ) then
                if(MasterProc) write(*,*) "MOSSPEC not found ", iadv, trim(name)
                cycle MC_LOOP
          end if
          call CheckStop( iadv < 1 .or. iadv > NSPEC_ADV, &
                 " ERR: Mc  _SPECS: Mc_SPECS" )

          nMosaic = nMosaic + 1
          call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, nVg" )

          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC, scale, avg? rho Inst Yr Mn Day atw

      !------------- Deposition velocities for d_2d -------------------------
          if( trim ( typ ) == "VG" ) then
            MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", typ , lctxt,  "cm/s", &
                iadv, -99, F, 100.0, T, Mc(n)%ind ) ! ind gives iotype
          else if( typ == "Rs" )  then
             MosaicOutput(nMosaic) = Deriv( &
               name, "Mosaic", typ, lctxt, "s/m", &
                iadv, -99, F , 1.0,  T, Mc(n)%ind ) ! ind gives iotype

          end if
          if(  DEBUG .and. MasterProc ) write(*,*) "DEBUG nMc ", &
            trim(name) // ":" // &
            trim(Mc(n)%txt2) // ":" // trim(Mc(n)%txt3), iadv, iLC

      end do MC_LOOP
 end subroutine Add_NewMosaics



!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine Add_MosaicVEGO3(iotype,nVEGO3)

        integer, intent(in) :: iotype
        integer, intent(out) :: nVEGO3
        integer :: n, iLC
        character(len=TXTLEN_DERIV) :: name
        character(len=TXTLEN_SHORT) :: units
        real :: scale
        type(O3cl) :: veg
        logical :: dt_scale

      !------------- VEGO3 stuff ----------------------------------------------
      ! For fluxes or AOTs we start with a formatted name, eg. POD_3.0_CF and
      !untangle it to get threshold Y (=3.0) and landcover type

      nVEGO3 = 0
      dt_scale = .false. ! for POD and AOT we need to reset this.
      VEGO3_LC: do n = 1, size(VEGO3_OUTPUTS%name)

         veg = VEGO3_OUTPUTS(n)
         name = veg%name
         !txt = veg%LC

         if( veg%class == "POD" ) then
            units = "mmole/m2"
            scale = 1.0e-6   ! Accumulates nmole/s to mmole (*dt_advec)
            dt_scale = .true.   ! Accumulates nmole/s to mmole (*dt_advec)
         else if( veg%class == "AOT" )  then
            units = "ppb.h"
            scale = 1.0/3600.0 ! AOT in ppb.hour
            dt_scale = .true.
         else if(DEBUG ) then
           call StopAll( "MosaicOuputs: vegclass errror" // veg%class )
         end if


          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "VEGO3_LCS", veg%TXTLC, .true. )
          if ( iLC < 0  ) cycle  VEGO3_LC
          if ( iLC > 0  ) LandType(iLC)%flux_wanted  = .true. 
          !-------------End of Check if LC present in this array ------!
          nVEGO3 = nVEGO3 + 1

          name = veg%name

           nMosaic = nMosaic + 1
           call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, VEGO3" )
          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC, scale dt_scale avg? Inst Yr Mn Day
          ! Use index for veg array. No need to set iadv for VEGO3. Always O3.
           MosaicOutput(nMosaic) = Deriv(  &
              name, veg%class,  veg%defn, veg%TXTLC, units, &
                n, -99, T,  scale,  F,   iotype ) 

      end do VEGO3_LC !n
 end subroutine Add_MosaicVEGO3

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine Add_MosaicDDEP(DDEP_ECOS,DDEP_SPECS,DDEP_FREQ,nDD)

        character(len=*), dimension(:), intent(in) :: DDEP_ECOS
        integer, dimension(:), intent(in) :: DDEP_SPECS  ! eg NH3
        integer, intent(in) :: DDEP_FREQ ! Day, Month, 
        integer, intent(out) :: nDD
        integer :: i, n, ispec, iadv
        character(len=TXTLEN_DERIV) :: name
        character(len=TXTLEN_SHORT) :: units
        real :: atw


     !------------- Dry Depositions for d_2d -------------------------
     ! Add species and ecosystem depositions if wanted:
     ! We find the various combinations of gas-species and ecosystem,
     ! adding them to the derived-type array OutDDep (e.g. => D2_SO4_m2Conif)

      nDD = 0
      do i = 1, size(DDEP_SPECS)
        do n = 1, size(DDEP_ECOS)

          nDD = nDD + 1
          ispec  = DDEP_SPECS(i)  ! Index in ix_tot arrays

          if ( ispec > 0 ) then ! normal species, e.g. DDEP_SO2_m2CF
             name = "DDEP_"  // trim( species(ispec)%name ) // &
                    "_m2" // trim( DDEP_ECOS(n) )

             iadv  = ispec - NSPEC_SHL ! adv for use in Derived

             if ( species(ispec)%sulphurs > 0 ) then
               atw  = species(ispec)%sulphurs * atwS
               units  =  "mgS/m2"

                call CheckStop( species( ispec )%nitrogens > 0 , &
                 " ERROR in DDEP_SPECS: BOTH S AND N!"// &
                   species(DDEP_SPECS(i))%name)

             else if ( species( ispec )%nitrogens > 0 ) then
               atw  = species( ispec )%nitrogens * atwN
               units  =  "mgN/m2"

             else if ( species(ispec)%nitrogens ==  0 .and. &
                      species(ispec)%sulphurs  == 0 ) then
                atw = species(ispec)%molwt 
                write(*,*) "Mosaic Molweight ", trim(species(ispec)%name), atw
                units = "mg/m2"

             else
               call StopAll("ERROR: OutDDep atw failure "// &
                   species( ispec )%name)
             end if
          else ! GROUP
             iadv   = ispec  ! e.g. -1 for SOX
             atw   = atwN   ! first guess
             units = "mgN/m2"
            if ( ispec == SOX_INDEX ) then
                atw   = atwS
                units = "mgS/m2"
                name = "DDEP_SOX_m2"//trim(DDEP_ECOS(n))
             else if ( ispec == OXN_INDEX ) then
                name = "DDEP_OXN_m2"//trim(DDEP_ECOS(n))
             else if ( ispec == RDN_INDEX ) then
                name = "DDEP_RDN_m2"//trim(DDEP_ECOS(n))
             end if
          end if

             nMosaic = nMosaic + 1
             call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, DDEP" )

          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,dt_scale, scale, avg? Inst Yr Mn Day

             MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", "DDEP", DDEP_ECOS(n), units, &
                  iadv,-99, F, 1.0e6 * atw ,  F,  DDEP_FREQ ) 

          if(DEBUG .and. MasterProc) then
            write(6,*) "DDEP setups"
            call print_deriv_type(MosaicOutput(nMosaic))
          end if
        end do ! DDEP_SPECS
     end do ! DDEP_ECOS
 end subroutine Add_MosaicDDEP

!<==========================================================================

 subroutine Add_MosaicOutput(debug_flag,i,j,convfac,DepAdv2Calc,fluxfrac,&
                Deploss)

  !<==========================================================================
     ! Adds deposition losses to ddep arrays
     logical, intent(in) :: debug_flag
     integer, intent(in) :: i,j             ! coordinates
     real,    intent(in) ::  convfac  !!???, lossfrac
     integer, dimension(:), intent(in) :: DepAdv2Calc
     real, dimension(:,:), intent(in) :: fluxfrac  ! dim (NADV, NLANDUSE)
     real, dimension(:), intent(in) :: Deploss

     integer, dimension(NMAX_DDEP) :: ddep_code 
     integer :: n, nadv, nadv2, iLC, iEco
     integer :: imc, f2d, cdep
     real :: output     ! tmp variable
     character(len=TXTLEN_SHORT) :: subclass, class
     logical :: my_first_call = .true.
     logical :: first_call = .true.     ! reset each subroutine call

  ! Variables added for ecosystem dep
     real, dimension(NDEF_ECOSYSTEMS) :: invEcoFrac, EcoFrac
     real :: Fflux, Gs, Gns

     cdep = -99                      ! set on first_vgr_call

  ! Must match areas given above, e.g. DDEP_CONIF -> Conif

  !Ecosystem areas, which were assigned in Init_DryDep:
  !  EcoFrac(CONIF)   = sum( coverage(:), LandType(:)%is_conif )
  !  EcoFrac(FULL_GRID)    = 1.0

     EcoFrac(:)    = EcoSystemFrac(:,i,j)
     invEcoFrac(:) = 0.0

     do n = 1, NDEF_ECOSYSTEMS
        if ( EcoFrac(n) > 1.0e-39 ) invEcoFrac(n) = 1.0/EcoFrac(n)
     end do 

   !  Query - crops, outisde g.s. ????
    if ( first_call ) then  ! need to find indices
       do imc = 1, nMosaic
         MosaicOutput(imc)%f2d  = find_index(MosaicOutput(imc)%name ,f_2d(:)%name)
         if(DEBUG .and. MasterProc) write(*,*) "MOS f2D", imc, &
             trim(MosaicOutput(imc)%name),  MosaicOutput(imc)%f2d
       end do

       if ( DEBUG .and. debug_flag ) then
          write(*,*)  "ECOAREAS ", i,j
             do n = 1,  NDEF_ECOSYSTEMS
                 write(*,"(a,i3,a,f14.4,g12.3)")  "ECOCHECK ", n, &
               DEF_ECOSYSTEMS(n), EcoFrac(n), invEcoFrac(n)
               end do
          write(*,*) "Done ECOCHECK ========================"
       end if
        
     end if
     first_call = .false.


    ! Traditional depositions to whole grid::
!     d_2d(DDEP_SOX,i,j,IOU_INST) = (  &
!          DepLoss(IXADV_SO2) + DepLoss(IXADV_SO4) ) * convfac * atwS
!
!       d_2d(D2_VddCOA,i,j,IOU_INST) =  100.0* Vg3m(CDDEP_COA)

    ! Ecosystem depositions, for grouped or individual species:

     do imc = 1, nMosaic
        class = MosaicOutput(imc)%class
        subclass = MosaicOutput(imc)%subclass
        f2d      = MosaicOutput(imc)%f2d
        nadv     = MosaicOutput(imc)%Index  ! can be negatve for groups
        iLC  = find_MosaicLC(imc)   ! Used for many cases, but replaced by iEco sometimes

       if( class == "AOT" ) subclass = class
       if( class == "POD" ) subclass = class

        output = 0.0  ! We only have instantaneous outputs, so can initialise
                      ! here and set d-2d at end

        !if ( my_first_call ) then ! Some safety tests.
        !   ! Land-cover index can be zero, for FULL_GRID
        !    iLC = find_MosaicLC( imc )
        !    call CheckStop(iLC<FULL_LCGRID, "ILC ERROR: "//MosaicOutput(imc)%name)
        !    call CheckStop(f2d<1, "f2d ERROR:  "//MosaicOutput(imc)%name)
        !end if
        if ( DEBUG .and. debug_flag ) then
           write(6,"(a,a)",advance='no') "Add_Mosaic: "// &
                  trim(MosaicOutput(imc)%name), ", " // trim(subclass)
        end if

        select case ( subclass )
        case ( "DDEP" )

           iEco   = find_index( MosaicOutput(imc)%txt, DEF_ECOSYSTEMS) 
                         ! Eco landcovers can include several
                         ! land-cover classes, see EcoSystem_ml
           if ( nadv > 0 ) then  ! normal advectde species
              nadv2 = 1
              ddep_code(1) = nadv
           else if ( nadv == SOX_INDEX ) then
              nadv2 = size( DDEP_SOXGROUP )
              ddep_code(1:nadv2) = DDEP_SOXGROUP - NSPEC_SHL
           else if ( nadv == OXN_INDEX ) then
              nadv2 = size( DDEP_OXNGROUP )
              ddep_code(1:nadv2) = DDEP_OXNGROUP - NSPEC_SHL
           else  !if ( nadv == RDN_INDEX ) then
              nadv2 = size( DDEP_RDNGROUP )
              ddep_code(1:nadv2) = DDEP_RDNGROUP - NSPEC_SHL
           end if
   
           Fflux = 0.0
           do n = 1, nadv2
              nadv = ddep_code(n)
              Fflux = Fflux + Deploss(nadv) * &
                sum( fluxfrac(nadv,:), Is_EcoSystem(iEco,:) )
           end do ! n

           if ( DEBUG .and. Fflux < 0.0 ) then
             write(6,"(a,3i4,a)") "DDEP Fflux CATASTR ", imc, f2d, iEco, &
                    trim(MosaicOutput(imc)%name)
             call CheckStop("CATASTROPHE: "//MosaicOutput(imc)%name)
           end if

        ! - invEcoFracCF divides the flux per grid by the landarea of each
        ! ecosystem, to give deposition in units of mg/m2 of ecosystem.

        !JAN30  output = Fflux * convfac * MosaicOutput(imc)%atw * invEcoFrac(iEco)
             output = Fflux * convfac * invEcoFrac(iEco)

          if ( DEBUG .and. debug_flag ) then
             write(6,"(3i4,3es12.3)") imc, nadv, iEco, Fflux, output
          end if ! DEBUG_ECO 

        case ( "METCONC" )    ! hard-coded bit n' pieces

          n   =  MosaicOutput(imc)%Index  !ind = ustar or ..
          if( n == MMC_USTAR  ) output = Sub(iLC)%ustar
          if( n == MMC_RH     ) output = Sub(iLC)%rh
          if( n == MMC_INVL   ) output = Sub(iLC)%invL
          if( n == MMC_CANO3  ) output = Sub(iLC)%cano3_ppb
          if( n == MMC_VPD    ) output = Sub(iLC)%vpd
          if( n == MMC_FST    ) output = Sub(iLC)%FstO3
          if( n == MMC_GSTO   ) output = Sub(iLC)%g_sto
          if( n == MMC_EVAP   ) output = Sub(iLC)%EvapTransp

          if ( DEBUG .and. debug_flag .and. &
               n==MMC_CANO3 .and. iLC == 2 ) then !DF
             write(6,"(a,4i5,f10.4)") "MYDDEP CANO3 ", &
                 current_date%month, current_date%day, &
                 current_date%hour, current_date%seconds,   output
          end if ! DEBUG_ECO 

        case ( "POD" )    ! Fluxes, PODY (was AFstY)

          n   =  MosaicOutput(imc)%Index !Index in VEGO3_OUPUTS
          call Calc_POD( n, iLC, output, debug_flag) 

        case ( "AOT" )    ! AOTX

          n   =  MosaicOutput(imc)%Index !Index in VEGO3_OUPUTS
          call Calc_AOTx( n, iLC, output, debug_flag) 

        case ( "VG", "Rs ", "Rns", "Gns" ) ! could we use RG_LABELS? 

             cdep = DepAdv2Calc( nadv ) ! e.g. IXADV_O3 to calc index
             Gs   = Sub(iLC)%Gsur(cdep)
             Gns  = Sub(iLC)%Gns(cdep)

           !It is easy to make mistakes with Vg, so we have som extra checks
           !here

           if ( DEBUG .and. cdep < 1 ) then
             print *, "ERROR: OutVgR name", MosaicOutput(imc)%name
             print *, "ERROR: Negative cdep", cdep, imc, MosaicOutput(imc)%Index
             print *, "ERROR: DEPADV2CALC had size", size(DepAdv2Calc)
             do n = 1, size( DepAdv2Calc)
                print *, "DEPADVLIST ", n, DepAdv2Calc(n)
             end do
             call CheckStop( cdep  < 1 , "ERROR: Negative cdep")
           end if

           if ( subclass == "VG" ) then

             output = Sub(iLC)%Vg_3m(cdep)  ! CHECK iLC
          if ( DEBUG .and. debug_flag ) then !DF
             call datewrite(":: ",iLC,(/ output /) )
             !call datewrite("MYDDEP VGVG "//trim( MosaicOutput(imc)%name),&
             !    iLC,  (/ output /) )
          end if

           else if ( subclass == "Gs" ) then

             output = Gs

           else if ( subclass == "Gns" ) then

             output = Gns

           else if ( subclass == "Rs" ) then

             if( Gs < 1.0e-44 ) then
                 output = -999.0 
             else
                 output = 1.0/Gs
             end if

           else if ( subclass == "Rns" ) then

            if( Gns < 1.0e-44 ) then
               output = -999.0
            else
               output = 1.0/Gns
            end if
           end if ! subclass

           if( DEBUG .and. debug_flag ) then
              write(*,"(2i4,f9.3)") cdep, iLC, output

           end if
        case default
           if ( MasterProc ) then
              print *, "OUTVEG UNDEF name ",   MosaicOutput(imc)%name
              print *, "OUTVEG UNDEF subclass ",   MosaicOutput(imc)%subclass
           call CheckStop("OUTVEG UNDEF" // subclass )
           end if
        end select

        if( DEBUG .and. debug_flag ) then
              write(*,"(a,es12.3)") "ADDED output: ",  output
        end if
        d_2d( f2d,i,j,IOU_INST) = output
     
     end do ! Mosaic

     my_first_call = .false.


  end subroutine  Add_MosaicOutput
  !<==========================================================================
  function  find_MosaicLC(imc) result(iLC)
     ! Searches for index in LandType array, and sets to zero if grid or EU
     ! (simple find:_index would set negati9ve if not found)
     integer, intent(in) :: imc
     integer :: iLC

     if( MosaicOutput(imc)%txt == "Grid") then
         iLC = FULL_LCGRID   ! zero
     else if( MosaicOutput(imc)%txt == "EU") then
         iLC = FULL_LCGRID   ! zero
     else 
         iLC = find_index( MosaicOutput(imc)%txt, LandDefs(:)%code )
     end if
  end function find_MosaicLC

end module MosaicOutputs_ml
