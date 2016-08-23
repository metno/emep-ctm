! <MosaicOutputs_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2015 met.no
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
! <MosaicOutputs_ml.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
module MosaicOutputs_ml
use AOTx_ml,          only: Calc_AOTx, Calc_POD, O3cl_t, VEGO3_OUTPUTS,&
                                nOutputVegO3, OutputVegO3
!use AOTx_ml,          only: Calc_SPOD
use CheckStop_ml,     only: CheckStop
use ChemSpecs,        only: NSPEC_ADV, species_adv
use ChemGroups_ml,    only: chemgroups
use DerivedFields_ml, only: f_2d, d_2d
use EcoSystem_ml,     only: NDEF_ECOSYSTEMS, DEF_ECOSYSTEMS, EcoSystemFrac, &
                            FULL_ECOGRID, FULL_LCGRID, Is_EcoSystem
use Io_Progs_ml,      only: datewrite
use LandDefs_ml,      only: LandDefs, LandType, Check_LandCoverPresent ! e.g. "CF"
use Landuse_ml,       only: LandCover ! for POD
use LocalVariables_ml,only: Grid,SubDat, L
use MetFields_ml
use ModelConstants_ml,only: MasterProc, DEBUG, &
                            NLANDUSEMAX, IOU_INST
use OwnDataTypes_ml,  only: Deriv, print_deriv_type, typ_s5i, typ_si, typ_s3,&
                            TXTLEN_DERIV, TXTLEN_SHORT
use SmallUtils_ml,    only: find_index
use SubMet_ml,        only: Sub
use TimeDate_ml,      only: current_date, effectivdaynumber, print_date
use Units_ml,         only: Units_Scale,Group_Scale,group_umap
use Wesely_ml,        only: NDRYDEP_CALC, CDDEP_O3

implicit none
private

public :: Init_MosaicMMC
public :: Add_MosaicMetConcs
public :: Add_MosaicVEGO3
public :: Add_MosaicDDEP
public :: Add_NewMosaics
public :: Add_MosaicOutput
public :: find_MosaicLC

INCLUDE 'mpif.h'
INTEGER STATUS(MPI_STATUS_SIZE),INFO

integer, public, save :: MMC_RH, MMC_CANO3, MMC_VPD, MMC_FST, &
  MMC_USTAR, MMC_INVL, MMC_GSTO, MMC_EVAP, MMC_LAI
character(len=30),private, save :: errmsg = "ok"

! Mosaic-specific outputs, e.g. VG_CF_HNO3 or Rns_GR_NH3
integer, public, save :: nMosaic = 0
integer, public, parameter :: MAX_MOSAIC_OUTPUTS=150
logical, private, parameter :: T=.true., F=.false.

type(Deriv), public, &
  dimension( MAX_MOSAIC_OUTPUTS ), save :: MosaicOutput

type(group_umap), private, target, &
  dimension( MAX_MOSAIC_OUTPUTS ), save :: dryGroupUnits

logical, private, save :: dbg0  !if DEBUG and MasterProc


contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Init_MosaicMMC(MOSAIC_METCONCS)
! Set indices of mosaic metconc params for later use. Will be zero if 
! not found, but that's okay I hope...
  character(len=*), dimension(:), intent(in) :: MOSAIC_METCONCS

  dbg0 = (DEBUG%MOSAICS .and. MasterProc) 

  MMC_RH    = find_index("RH",MOSAIC_METCONCS)
  MMC_CANO3 = find_index("CanopyO3",MOSAIC_METCONCS)
  MMC_VPD   = find_index("VPD",MOSAIC_METCONCS)
  MMC_FST   = find_index("FstO3",MOSAIC_METCONCS)
  MMC_USTAR = find_index("USTAR",MOSAIC_METCONCS)
  MMC_INVL  = find_index("INVL",MOSAIC_METCONCS)
  MMC_GSTO  = find_index("GSTO",MOSAIC_METCONCS)
  MMC_EVAP  = find_index("EVAP",MOSAIC_METCONCS)
  MMC_LAI   = find_index("LAI",MOSAIC_METCONCS)
end subroutine Init_MosaicMMC
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Add_MosaicMetConcs(MOSAIC_METCONCS,MET_LCS,iotyp, nMET)
  character(len=*), dimension(:), intent(in) :: MOSAIC_METCONCS
  character(len=*), dimension(:), intent(in) :: MET_LCS
  integer, intent(in)  :: iotyp
  integer, intent(out) :: nMET
  integer :: ilab, n, iLC
  character(len=TXTLEN_DERIV) :: name
  character(len=17), save :: dtxt='AddMosaicMetConc:'

  !------------- Met data for d_2d -------------------------
  ! We find the various combinations of met and ecosystem,
  ! adding them to the derived-type array Mosaic_Met (e.g. => Met_CF)
  nMET = 0
  do ilab = 1, size(MOSAIC_METCONCS)
    MET_LC: do n = 1, size(MET_LCS)

      !------------------- Check if LC present in this array ------!
      iLC = Check_LandCoverPresent( "MET_LCS", n, MET_LCS, (ilab == 1))
      if(iLC<0) cycle  MET_LC
      if(MOSAIC_METCONCS(ilab)(1:6)=="Canopy" .or.& 
         MOSAIC_METCONCS(ilab)(1:5)=="FstO3")     &
             LandType(iLC)%flux_wanted  = .true.  ! Canopy calc in StoFlux
      !-------------End of Check if LC present in this array ------!

      nMET = nMET + 1
      name = trim ( MOSAIC_METCONCS(ilab) ) // "_"  // trim( MET_LCS(n) )

      nMosaic = nMosaic + 1
      call CheckStop(NMosaic>=MAX_MOSAIC_OUTPUTS,dtxt//"too many nMosaics, nMET")
      !Deriv(name, class,    sub,   txt,           unit
      !Deriv index, f2d, scale, avg? Inst Yr Mn Day
      MosaicOutput(nMosaic) = Deriv(  &
          name, "Mosaic", "METCONC", MET_LCS(n), MOSAIC_METCONCS(ilab), &
            ilab, -99,F , 1.0,  T,  iotyp )
      
      select case(MOSAIC_METCONCS(ilab))
        case("USTAR"   );MosaicOutput(nMosaic)%unit = "m/s"
        case("LAI"     );MosaicOutput(nMosaic)%unit = "m2/m2"
        case("INVL"    );MosaicOutput(nMosaic)%unit = "m"
        case("CanopyO3");MosaicOutput(nMosaic)%unit = "ppb"
        case("FstO3"   );MosaicOutput(nMosaic)%unit = "mmole/m2" ! accumulated
        case("EVAP"    );MosaicOutput(nMosaic)%unit = "mm"
          MosaicOutput(nMosaic)%avg       =  .false. ! accumulate
          MosaicOutput(nMosaic)%dt_scale  =  .true.
      endselect

      if(dbg0) call print_deriv_type(MosaicOutput(nMosaic))
    enddo MET_LC !n
  enddo ! ilab
 endsubroutine Add_MosaicMetConcs
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine Add_NewMosaics(Mc,nMc)
  type(typ_s5i), dimension(:), intent(in) :: Mc ! eg VG
  integer, intent(out) :: nMc
  integer :: n, itot, iLC, iadv
  character(len=TXTLEN_DERIV) :: name
  character(len=TXTLEN_SHORT) :: typ, poll, lctxt
  character(len=14), save :: dtxt='AddNewMosaics:'

  nMc = 0
  MC_LOOP: do n = 1, size( Mc(:)%txt1 )
    !------------------- Check if LC present in this array ------!
    iLC = Check_LandCoverPresent( "MMC-VG", n, Mc(:)%txt4 , &
             write_condition=.true. )
    if(iLC<0) cycle MC_LOOP
    !-------------End of Check if LC present in this array ------!
    nMc = nMc + 1
    typ   = Mc(n)%txt2                                  ! VG
    poll  = Mc(n)%txt3                                  ! O3, HNO3, ...
    lctxt = Mc(n)%txt4                                  ! Grid, SNL,..
    name = trim(typ)//"_"//trim(poll)//"_"//trim(lctxt) ! VG_O3_GRID?
    
    iadv = find_index(poll,species_adv(:)%name )
    if(iadv<1) then
      if(MasterProc) write(*,*) "MOSSPEC not found ", iadv, trim(name)
      cycle MC_LOOP
    endif
    call CheckStop(iadv<1 .or. iadv>NSPEC_ADV,dtxt//" ERR: Mc  _SPECS: Mc_SPECS")

    nMosaic = nMosaic + 1
    call CheckStop(NMosaic>=MAX_MOSAIC_OUTPUTS,dtxt//"too many nMosaics, nVg" )

    !------------- Deposition velocities for d_2d -------------------------
    !Deriv(name, class,    subc,  txt,           unit
    !Deriv index, f2d,LC, scale, avg? rho Inst Yr Mn Day atw
    select case(typ)
    case("VG")
      MosaicOutput(nMosaic) = Deriv(  &
        name, "Mosaic", typ , lctxt,  "cm/s", &
          iadv, -99, F, 100.0, T, Mc(n)%ind ) ! ind gives iotype
    case("Rs")
       MosaicOutput(nMosaic) = Deriv( &
         name, "Mosaic", typ, lctxt, "s/m", &
          iadv, -99, F , 1.0,  T, Mc(n)%ind ) ! ind gives iotype
    endselect

    if(dbg0) write(*,*) "DEBUG nMc ", &
      trim(name)//":"//trim(Mc(n)%txt2)//":"//trim(Mc(n)%txt3), iadv, iLC
  enddo MC_LOOP
end subroutine Add_NewMosaics
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Add_MosaicVegO3(nVEGO3)
  integer, intent(out) :: nVEGO3
  integer :: n, iLC
  character(len=TXTLEN_DERIV) :: name
  character(len=TXTLEN_SHORT) :: units
  character(len=*), parameter :: dtxt='AddMosaicVEGO3:'
  real :: scale
  type(O3cl_t) :: veg
  logical :: dt_scale
  !------------- VEGO3 stuff ----------------------------------------------
  ! For fluxes or AOTs we start with a formatted name, eg. POD_3.0_CF and
  !untangle it to get threshold Y (=3.0) and landcover type
  nVEGO3 = 0
  dt_scale = .false. ! for POD and AOT we need to reset this.

  VEGO3_LC: do n = 1, nOutputVegO3

     veg = OutputVegO3(n)
     name = veg%name
   
     select case(veg%class)
     case("POD")
       units = "mmole/m2"
       scale = 1.0e-6     ! Accumulates nmole/s to mmole (*dt_advec)
       dt_scale = .true.  ! Accumulates nmole/s to mmole (*dt_advec)
     case("AOT", "SPOD")
       units = "ppb.h"
       scale = 1.0/3600.0 ! AOT in ppb.hour
       dt_scale = .true.
     case default
       call CheckStop(DEBUG%MOSAICS,dtxt//"vegclass errror"//veg%class )
     endselect
   
     !------------------- Check if LC present in this array ------!
     iLC = Check_LandCoverPresent( "VEGO3_LCS", veg%TXTLC, .true. )
     if(iLC<0) cycle  VEGO3_LC
     if(iLC>0) LandType(iLC)%flux_wanted  = .true. 
     !-------------End of Check if LC present in this array ------!
     nVEGO3 = nVEGO3 + 1
     name = veg%name
   
     nMosaic = nMosaic + 1
     call CheckStop(NMosaic>=MAX_MOSAIC_OUTPUTS,dtxt//"too many nMos..EGO3")
      if(dbg0) then
        write(*,*) "Moscaics", nMosaic, trim(name)// "->" //trim(veg%TXTLC)
      end if
     !Deriv(name, class,    subc,  txt,           unit
     !Deriv index, f2d,LC, scale dt_scale avg? Inst Yr Mn Day
     ! Use index for veg array. No need to set iadv for VEGO3. Always O3.
      MosaicOutput(nMosaic) = Deriv( name, veg%class,  veg%defn, veg%TXTLC&
       , units, n, -99, dt_scale,  scale,  F,  veg%iotype ) 
   
  enddo VEGO3_LC !n
end subroutine Add_MosaicVEGO3
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine Add_MosaicDDEP(DDEP_ECOS,DDEP_WANTED,nDD)
  type(typ_si), dimension(:), intent(in) :: DDEP_ECOS    !e.g. (%name,DDEP_FREQ)
  type(typ_s3), dimension(:), intent(in) :: DDEP_WANTED  !e.g. ("NH3","SPECIE","mgN")
  integer, intent(out) :: nDD
  integer :: i, n, iadv,igrp
  character(len=TXTLEN_DERIV) :: name='-', xname, xtyp, xxname
  character(len=TXTLEN_SHORT) :: units
  character(len=14), save :: dtxt='AddMosaicDDEP:' ! for debug output
  real :: unitscale
!------------- Dry Depositions for d_2d -------------------------
! Add species and ecosystem depositions if wanted:
! We find the various combinations of gas-species and ecosystem,
! adding them to the derived-type array OutDDep (e.g. => D2_SO4_m2Conif)
  nDD = 0
!dbg0=.true.
  do i=1,size(DDEP_WANTED)
    xname  = DDEP_WANTED(i)%txt1
    xtyp   = DDEP_WANTED(i)%txt2
    if ( xname == '-') exit
    if(dbg0) write(*,*) dtxt//"DDEP_WANTED,a:"//trim(xname), i, xtyp

    do n=1,size(DDEP_ECOS)

      if(dbg0) write(*,"(a,2i4,2a12)") dtxt//"DDEP_WANTED,b:", n, DDEP_ECOS(n)%ind, &
         trim(DDEP_ECOS(n)%name), trim(xtyp)
      if (DDEP_ECOS(n)%ind < 1) exit
      nDD = nDD + 1
      nMosaic = nMosaic + 1

      call CheckStop(NMosaic>=MAX_MOSAIC_OUTPUTS,&
        dtxt//"too many nMosaics, DDEP "//trim(xname))

      select case(xtyp)
      case("SPEC")    ! normal species, e.g. DDEP_SO2_m2CF
!HACK :-(
xxname = xname
if( xname(1:4) == "STO_" ) then
     xxname = xname(5:)
     if( dbg0 ) print *, "STO_ ", trim(xname) // "=>"//trim(xxname)
end if
        iadv = find_index(xxname,species_adv(:)%name)         ! Index in ix_adv arrays
if( iadv < 1 ) print *, "OOPADVNOW", iadv, trim(xname) // trim(name)
        call CheckStop(iadv<1,dtxt//"Unknown in DDEP_WANTED SPEC: "//trim(xname))
        unitscale = Units_Scale(DDEP_WANTED(i)%txt3,iadv,units)
        if(dbg0) print *, "ADVNO ",n,iadv,trim(xname), unitscale

      case("GROUP")
        igrp = find_index(xname,chemgroups(:)%name) ! array of members: chemgroups%ptr
if( igrp < 1 ) print *, "OOPNOW", igrp, trim(xname) // trim(name)
        call CheckStop(igrp<1,dtxt//"Unknown in DDEP_WANTED GROUP: "//trim(xname))
        unitscale = Units_Scale(DDEP_WANTED(i)%txt3,-1,units)
                      ! Units_Scale(iadv=-1) returns 1.0
                      ! Add_MosaicOutput gets the unit conversion factor from Group_Scale
        iadv = -igrp  ! use negative values for groups (e.g. DDEP_SOX)
        dryGroupUnits(nMosaic) = Group_Scale(igrp,DDEP_WANTED(i)%txt3,&
                                             debug=dbg0)
      case default
        call CheckStop(DEBUG%MOSAICS,dtxt//" unknown MosaicDDEP type "//trim(DDEP_WANTED(i)%txt2))
      endselect
!print *, "OOP===========================POO"

      name = "DDEP_"//trim(xname)//"_m2"//trim(DDEP_ECOS(n)%name)

     !Deriv(name, class,    subc,  txt,           unit
     !Deriv index, f2d,dt_scale, scale, avg? Inst/Yr/Mn/Day
      MosaicOutput(nMosaic) = Deriv( &
        name, "Mosaic", "DDEP", DDEP_ECOS(n)%name, units, &
          iadv,-99, F, unitscale, F, DDEP_ECOS(n)%ind )

      if(dbg0) then
        write(*,*) "DDEP setups", n, nMosaic
        call print_deriv_type(MosaicOutput(nMosaic))
      endif
    enddo ! DDEP_SPECS
  enddo ! DDEP_ECOS
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

  character(len=*), parameter :: dtxt='AddMosc:'
  type(group_umap), pointer :: gmap=>null()  ! group unit mapping  
  integer :: n, nadv, iLC, iEco
  integer :: imc, f2d, cdep
  real :: output     ! tmp variable
  character(len=TXTLEN_SHORT) :: subclass, class, txtdate='-'
  logical :: my_first_call = .true.
  logical :: first_call = .true.     ! reset each subroutine call
  ! for SDEP (tmp)
   integer, save :: nd2d
   logical, save :: first_SDEP = .true.

  ! Variables added for ecosystem dep
  real, dimension(NDEF_ECOSYSTEMS) :: invEcoFrac, EcoFrac
  real :: Fflux, Gs, Gns
  logical :: dbg, dbghh

  cdep = -99                      ! set on first_vgr_call
  dbg   = DEBUG%MOSAICS.and.debug_flag
  dbghh = dbg .and. current_date%seconds == 0
  if(dbg) txtdate = print_date()

  ! Must match areas given above, e.g. DDEP_CONIF -> Conif

  !Ecosystem areas, which were assigned in Init_DryDep:
  !  EcoFrac(CONIF)   = sum( coverage(:), LandType(:)%is_conif )
  !  EcoFrac(FULL_GRID)    = 1.0

  EcoFrac(:)    = EcoSystemFrac(:,i,j)
  invEcoFrac(:) = 0.0
  do n=1,NDEF_ECOSYSTEMS
    if(EcoFrac(n)>1.0e-39) invEcoFrac(n)=1.0/EcoFrac(n)
  enddo 

  !  Query - crops, outisde g.s. ????
  if(first_call) then  ! need to find indices
    do imc = 1, nMosaic
     MosaicOutput(imc)%f2d  = find_index(MosaicOutput(imc)%name ,f_2d(:)%name)
     if(DEBUG%MOSAICS .and. MasterProc) write(*,*) dtxt//" f2D", imc, &
         trim(MosaicOutput(imc)%name),  MosaicOutput(imc)%f2d
    enddo

    if( dbg ) then
      write(*,*)  dtxt//"ECOAREAS ", i,j
      do n=1,NDEF_ECOSYSTEMS
        write(*,"(a,i3,a,f14.4,g12.3)")  dtxt//"ECOCHECK ", n, &
          DEF_ECOSYSTEMS(n), EcoFrac(n), invEcoFrac(n)
      enddo
      write(*,*) dtxt//"Done ECOCHECK ========================"
    endif       
  endif
  first_call = .false.


  ! Traditional depositions to whole grid::
!     d_2d(DDEP_SOX,i,j,IOU_INST) = (  &
!          DepLoss(IXADV_SO2) + DepLoss(IXADV_SO4) ) * convfac * atwS
!
!       d_2d(D2_VddCOA,i,j,IOU_INST) =  100.0* Vg3m(CDDEP_COA)

  ! Ecosystem depositions, for grouped or individual species:

  do imc = 1, nMosaic
    class    = MosaicOutput(imc)%class
    subclass = MosaicOutput(imc)%subclass
    f2d      = MosaicOutput(imc)%f2d
    nadv     = MosaicOutput(imc)%Index  ! can be negatve for groups
    iLC      = find_MosaicLC(imc)   ! Used for many cases, but replaced
                                    ! by iEco sometimes
    if( iLC > 0 ) L  = Sub(iLC) ! Avoid imc=0==Grid. L only used for POD,AOT

    if(class=="AOT") subclass=class
    if(class=="POD") subclass=class
    if(class=="SPOD") subclass=class

    output = 0.0  ! We only have instantaneous outputs, so can initialise
                  ! here and set d-2d at end

    !if ( my_first_call ) then ! Some safety tests.
    !   ! Land-cover index can be zero, for FULL_GRID
    !    iLC = find_MosaicLC( imc )
    !    call CheckStop(iLC<FULL_LCGRID, "ILC ERROR: "//MosaicOutput(imc)%name)
    !    call CheckStop(f2d<1, "f2d ERROR:  "//MosaicOutput(imc)%name)
    !end if
    !TMP if( dbg ) write(*,"(a,a)",advance='no') "Add_Mosaic: "// &
    if( dbg ) write(*,"(a,a)") dtxt//"Add_Mosaic: "// &
            trim(MosaicOutput(imc)%name), ", " // trim(subclass)

    select case(subclass)
    case("DDEP")
      ! Eco landcovers can include several land-cover classes, see EcoSystem_ml
      iEco = find_index(MosaicOutput(imc)%txt,DEF_ECOSYSTEMS) 
      select case(nadv)
      case(1:NSPEC_ADV)                 ! normal advected species
        Fflux = Deploss(nadv) &
               *sum(fluxfrac(nadv,:),Is_EcoSystem(iEco,:))
       if( species_adv(nadv)%name == "O3" .and.  &
            index( MosaicOutput(imc)%name, "STO_")>0 ) then
 !HACK
          Fflux = Fflux * Sub(0)%Gsto(CDDEP_O3)/Sub(0)%Gsur(CDDEP_O3)
          if( dbghh ) then
            print "(a,2i4,3es12.3)", dtxt//"NOWSDEP ", nadv, current_date%hour,  Sub(0)%Gsur(2), Sub(0)%Gsto(2), Fflux
         end if
      end if
 !HACK
      case(-size(chemgroups):-1)        ! gropups
        gmap=>dryGroupUnits(imc)
        Fflux = 0.0
        do n=1,size(gmap%iadv)
          nadv = gmap%iadv(n)
          Fflux = Fflux + Deploss(nadv)*gmap%uconv(n) &
                         *sum(fluxfrac(nadv,:),Is_EcoSystem(iEco,:))
        enddo ! n
      case default
        call CheckStop(dtxt//" unknown DDEP Specie/Group")
      endselect

      if(DEBUG%MOSAICS.and.Fflux<0.0) then
        write(*,"(a,3i4,a)") dtxt//"DDEP Fflux CATASTR ", imc, f2d, iEco, &
              trim(MosaicOutput(imc)%name)
        call CheckStop(dtxt//"CATASTROPHE: "//MosaicOutput(imc)%name)
      endif

      ! - invEcoFracCF divides the flux per grid by the landarea of each
      ! ecosystem, to give deposition in units of mg/m2 of ecosystem.

      output = Fflux * convfac * invEcoFrac(iEco)
      if( dbg ) write(*,"(3i4,3es12.3)") imc, nadv, iEco, Fflux, output

      case("METCONC")     ! hard-coded bit n' pieces
        n =  MosaicOutput(imc)%Index  !ind = ustar or ..
        if(n==MMC_USTAR  ) output = Sub(iLC)%ustar
        if(n==MMC_RH     ) output = Sub(iLC)%rh
        if(n==MMC_INVL   ) output = Sub(iLC)%invL
        if(n==MMC_CANO3  ) output = Sub(iLC)%cano3_ppb
        if(n==MMC_VPD    ) output = Sub(iLC)%vpd
        if(n==MMC_FST    ) output = Sub(iLC)%FstO3
        if(n==MMC_GSTO   ) output = Sub(iLC)%g_sto
        if(n==MMC_EVAP   ) output = Sub(iLC)%EvapTransp
        if(n==MMC_LAI    ) output = Sub(iLC)%LAI

        if(dbg.and.n==MMC_CANO3.and.iLC==2) & !DF
           write(*,"(2a,g12.4)") dtxt//"MYDDEP CANO3 ", &
            trim(txtdate), output

      case("POD")         ! Fluxes, PODY (was AFstY)
        n =  MosaicOutput(imc)%Index !Index in VEGO3_OUPUTS
        call Calc_POD( n, iLC, output, debug_flag) 
        if(dbg) write(*,"(2a,g12.4)") dtxt//"MYPOD ", &
            trim(txtdate), output


      case("AOT")         ! AOTX
        n =  MosaicOutput(imc)%Index !Index in VEGO3_OUPUTS
        if(dbg.and. Sub(iLC)%cano3_ppb> 40.0) &
             write(*,*) dtxt//" preAOT", n,iLC, Sub(iLC)%cano3_ppb
        call Calc_AOTx( n, iLC, output ) 

      case ( "VG", "Rs ", "Rns", "Gns" ) ! could we use RG_LABELS? 
        cdep = DepAdv2Calc( nadv ) ! e.g. IXADV_O3 to calc index
        Gs   = Sub(iLC)%Gsur(cdep)
        Gns  = Sub(iLC)%Gns(cdep)

       !It is easy to make mistakes with Vg, so we have som extra checks here
        if(DEBUG%MOSAICS.and.cdep<1) then
          write(*,*) "ERROR: OutVgR name", MosaicOutput(imc)%name
          write(*,*) "ERROR: Negative cdep", cdep, imc, MosaicOutput(imc)%Index
          write(*,*) "ERROR: DEPADV2CALC had size", size(DepAdv2Calc)
          do n = 1, size( DepAdv2Calc)
            write(*,*) "DEPADVLIST ", n, DepAdv2Calc(n)
          enddo
          call CheckStop(cdep<1, dtxt//"ERROR: Negative cdep")
        endif

        select case(subclass)
        case("VG" )
          output = Sub(iLC)%Vg_3m(cdep)  ! CHECK iLC
          if(dbg) then !DF
            call datewrite(":: ",iLC,(/ output /) )
           !call datewrite("MYDDEP VGVG "//trim( MosaicOutput(imc)%name),&
           !    iLC,  (/ output /) )
          endif
        case("Gs" )
          output = Gs
        !case("Gns") ! Apr 2015. Not likely to need these?
        !  output = Gns
        !case("Rs" )
        !  if(Gs < 1.0e-44)then
        !    output = -999.0 
        !  else
        !    output = 1.0/Gs
        !  endif
        !case("Rns")
        !  if(Gns < 1.0e-44)then
        !    output = -999.0
        !  else
        !    output = 1.0/Gns
        !  endif
        endselect ! subclass

        if(dbg)  write(*,"(2i4,f9.3)") cdep, iLC, output

      case default
        if(MasterProc) then
          write(*,"(/1('OUTVEG UNDEF',2(1x,A)))"),&
            "name"    ,MosaicOutput(imc)%name,&
            "subclass",MosaicOutput(imc)%subclass
          call CheckStop("OUTVEG UNDEF" // subclass )
        endif
      endselect

      if(dbg) write(*,"(a,es12.3)") dtxt//"ADDED output: "// &
         trim(MosaicOutput(imc)%name),  output
      d_2d( f2d,i,j,IOU_INST) = output
    enddo ! Mosaic

     my_first_call = .false.
end subroutine Add_MosaicOutput
!<==========================================================================
function  find_MosaicLC(imc) result(iLC)
! Searches for index in LandType array, and sets to zero if grid or EU
! (simple find:_index would set negati9ve if not found)
  integer, intent(in) :: imc
  integer :: iLC

  select case(MosaicOutput(imc)%txt)
    case("Grid");iLC = FULL_LCGRID   ! zero
    case("EU"  );iLC = FULL_LCGRID   ! zero
    case default;iLC = find_index(MosaicOutput(imc)%txt, LandDefs(:)%code)
  endselect
endfunction find_MosaicLC
!<==========================================================================
endmodule MosaicOutputs_ml
