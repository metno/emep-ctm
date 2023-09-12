! <Derived_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.0>
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
module Derived_mod
!---------------------------------------------------------------------------
! DESCRIPTION
! This module performs the calculations associated with "derived" 2D and 3D,
! such as accumulated precipitation or sulphate, daily, monthly or yearly
! averages, depositions. These fields are all typically output as netCDF
! fields.
!
! This routine defines many possible derived  outputs.
! The names of the derived fields actualy required should have been specified
!  in the user-defined My_Derived_mod.
!
! User-defined routines and treatments are often needed here. Here there is
! added stuff for AOTs, accsu. In
! general such code should be added in such a way that it isn't activated if
! not needed. It then doesn't need to be commented out if not used.
!  
! a line in OutputMisc will be added to def_2d in addDeriv from Derived_mod:
! call AddDeriv(OutputMisc(n),Is3D=Is3D) 
!---------------------------------------------------------------------------

use AeroConstants_mod, only: AERO
use AeroFunctions_mod, only: LogNormFracBelow !=> Frac Mass below Dp
use AOD_PM_mod,        only: AOD_init,aod_grp,wavelength,& ! group and
                                wanted_wlen,wanted_ext3d      ! wavelengths
use AOTx_mod,          only: Calc_GridAOTx
!PW: removed use BiDir_module,     only : Bidir_2d
use BiDir_emep,       only : Bidir_Derived
use Biogenics_mod,     only: EmisNat, NEMIS_BioNat, EMIS_BioNat
use CheckStop_mod,     only: CheckStop, StopAll
use Chemfields_mod,    only: xn_adv, xn_shl, cfac,xn_bgn, AOD,  &
                            SurfArea_um2cm3, &
                            Fgas3d, & ! FSOA
                            Extin_coeff, PM25_water, PM25_water_rh50 &
                          , PMco_water_rh50   !JUN21AERO
use Chemfields_mod ,   only: so2nh3_24hr,Grid_snow, Dobson, pH
use ChemDims_mod,      only: NSPEC_ADV, NSPEC_SHL,NEMIS_File
use ChemGroups_mod          ! SIA_GROUP, PMCO_GROUP -- use tot indices
use ChemSpecs_mod           ! IXADV_ indices etc
use Config_module,     only: &
   KMAX_MID,KMAX_BND  & ! =>  z dimension: layer number,level number
  ,KCHEMTOP           & ! limit of Fgas3d
  ,NPROC              & ! No. processors
  ,dt_advec           &
  ,PPBINV             & ! 1.0e9, for conversion of units
  ,PPTINV             & ! 1.0e12, for conversion of units
  ,PT                 &
  ,NTDAY              & ! Number of 2D O3 to be saved each day (for SOMO)
  ,num_lev3d,lev3d    & ! 3D levels on 3D output
  ! output types corresponding to instantaneous,year,month,day
  ,IOU_INST,IOU_YEAR,IOU_MON,IOU_DAY,IOU_HOUR,IOU_HOUR_INST,IOU_KEY &
  ,MasterProc, SOURCE_RECEPTOR, AOD_WANTED &
  ,USES, lf_src, startdate,enddate,&
  HourlyEmisOut, DailyEmisOut, SecEmisOutWanted, spinup_enddate, &
  OutputMisc, WDEP_WANTED, O3_ix

use Debug_module,      only: DEBUG   ! -> DEBUG%DERIVED and COLSRC
use DerivedFields_mod, only: MAXDEF_DERIV2D, MAXDEF_DERIV3D, &
                            def_2d, def_3d, f_2d, f_3d, d_2d, d_3d, VGtest_out_ix
use EcoSystem_mod,     only: DepEcoSystem, NDEF_ECOSYSTEMS, &
                            EcoSystemFrac,FULL_ECOGRID
use EmisDef_mod,       only: NSECTORS, EMIS_FILE, O_DMS, O_NH3&
                            ,SecEmisOut, EmisOut, SplitEmisOut, &
                            isec2SecOutWanted,SECTORS, Emis_CO_Profile
use EmisGet_mod,       only: nrcemis,iqrc2itot
use Functions_mod,      only: Tpot_2_T    ! Conversion function
use GasParticleCoeffs_mod, only: DDdefs
use GridValues_mod,    only: debug_li, debug_lj, debug_proc, A_mid, B_mid, &
                            dA,dB,xm2, GRIDWIDTH_M, GridArea_m2,xm_i,xm_j,glon,glat
use Io_Progs_mod,      only: datewrite
use MetFields_mod,     only: roa,Kz_m2s,th,zen, ustar_nwp, u_ref, hmix,&
                            met, derivmet,  q, &
                            ws_10m, rh2m, z_bnd, z_mid, u_mid,v_mid,ps, t2_nwp, &
                            cc3dmax, & ! SEI
                            dTleafHd, dTleafRn, & ! TLEAF
                            SoilWater_deep, SoilWater_uppr ,invL_nwp
use MosaicOutputs_mod,     only: nMosaic, MosaicOutput
use My_Derived_mod, only : &
    wanted_deriv2d, wanted_deriv3d, & ! names of wanted derived fields
    Init_My_Deriv, My_DerivFunc,    &
    OutputFields, nOutputFields,    &
    nOutputMisc,        &
    nOutputWdep, D3_OTHER,&
    PS_needed

use NumberConstants,      only: UNDEF_R
use OwnDataTypes_mod,      only: Deriv, print_Deriv_type, &
                                TXTLEN_DERIV,TXTLEN_SHORT,TXTLEN_IND ! type & length of names
use Par_mod,               only: me,                &      ! for print outs
                                limax, ljmax      ! => used x, y area
use PhysicalConstants_mod, only: PI,KAPPA,ATWAIR,GRAV
use OrganicAerosol_mod, only :  ORGANIC_AEROSOLS, Reset_3dOrganicAerosol
use ZchemData_mod,    only: Fpart ! for FSOA work
use SmallUtils_mod,        only: find_index, LenArray, NOT_SET_STRING, trims
use Tabulations_mod,      only : tab_esat_Pa
use TimeDate_mod,          only: day_of_year,daynumber,current_date,&
                                tdif_days
use TimeDate_ExtraUtil_mod,only: to_stamp, date_is_reached
use LocalFractions_mod,    only: lf_av, loc_frac
use Units_mod,             only: Units_Scale,Group_Units,&
                                to_molec_cm3 ! converts roa [kg/m3] to M [molec/cm3]
implicit none
private

public  :: Init_Derived
public  :: ResetDerived   ! Resets values to zero
public  :: DerivedProds   ! Calculates any production terms
public  :: AddDeriv       ! Adds Deriv type to def_2d, def_3d
public  :: AddNewDeriv    ! Creates & Adds Deriv type to def_2d, def_3d
private :: Define_Derived
public  :: wanted_iou     ! (iotyp, def%iotyp)
private :: write_debug
private :: write_debugadv

public  :: Derived        ! Calculations of sums, avgs etc.
private :: group_calc     ! Calculates sum of groups, e.g. pm25 from group array

logical, private, parameter :: T = .true., F = .false. ! shorthands only
integer, public, save :: num_deriv2d, num_deriv3d
integer, private,save :: Nadded2d = 0, Nadded3d=0 ! No. defined derived

! List of wanted IOUs
integer, parameter :: &
  IOU_MIN=lbound(IOU_KEY,DIM=1), &
  IOU_MAX=ubound(IOU_KEY,DIM=1)
logical, public, save :: &
  iou_list(IOU_MIN:IOU_MAX)=.false.

! The 2-d and 3-d fields use the above as a time-dimension. We define
! LENOUTxD according to how fine resolution we want on output. For 2d
! fields we use daily outputs. For the big 3d fields, monthly output
! is sufficient.

integer, public, parameter :: &
  LENOUT2D = IOU_HOUR,& ! Allows INST..DAY for 2d fields
  LENOUT3D = IOU_HOUR   ! Allows INST..DAY for 3d fields

!will be used for:
!e.g. d_2d( num_deriv2d,LIMAX, LJMAX, LENOUT2D)
! &   d_3d( num_deriv3d,LIMAX, LJMAX, num_lev3d, LENOUT3D )


! save O3 every hour during one day to find running max
real, save  , allocatable , public :: &     ! to be used for SOMO35
  D2_O3_DAY( :,:,:)

! Fraction of NO3_c below 2.5 um (v. crude so far)

real, save, private :: fracPM25 = -999.9

! Counters to keep track of averaging
! Initialise to zero in Init.

integer, public, allocatable, dimension(:,:), save :: nav_2d,nav_3d

logical, private, save :: Is3D
logical, private, save :: dbg0   ! = DEBUG%DERIVED .and. MasterProc
logical, private, save :: dbgP   ! = DEBUG%DERIVED .and. debug_proc
logical, private, save :: dbgP0  ! = DEBUG%DERIVED .and. debug_proc .and. first_call
character(len=100), private :: errmsg
! horizontal line for printouts
character(len=*), private, parameter :: HORIZ_LINE =repeat('=',78) !f2003://new_line('a') 

! NB global use of these common variables is dangerous!
integer, private :: i,j,k,l,n, iou, isec   ! Local loop variables

! Avoid hard codded IXADV_SPCS
 integer, private, save :: iadv_O3=-999, &
  iadv_OM25p=-999, igrp_OM25=-999, iadv_PMf=-999,     &
  iadv_NO3_C=-999,iadv_EC_C_WOOD=-999,iadv_EC_C_FFUEL=-999,iadv_POM_C_FFUEL=-999

real, private, save ::                      & ! Avoid hard codded molwt
  ug_NO3_C=-999.0,ug_EC_C_WOOD=-999.0,ug_EC_C_FFUEL=-999.0,ug_POM_C_FFUEL=-999.0

contains

!=========================================================================
subroutine Init_Derived()
  integer :: alloc_err
  integer :: iddefPMc
  character(len=*), parameter :: dtxt='IniDeriv:' !debug label
  dbg0 = (DEBUG%DERIVED .and. MasterProc )

  allocate(D2_O3_DAY( LIMAX, LJMAX, NTDAY))
  D2_O3_DAY = 0.0

  if (USES%LocalFractions .and. (lf_src(1)%HOUR .or. lf_src(1)%HOUR_INST)) HourlyEmisOut = .true.
  if (USES%LocalFractions .and. lf_src(1)%DAY) DailyEmisOut = .true.

  if(dbg0) write(*,*) dtxt//"INIT STUFF"
  call Init_My_Deriv()  !-> wanted_deriv2d, wanted_deriv3d

  ! get lengths of wanted arrays (excludes notset values)
  num_deriv2d = LenArray(wanted_deriv2d,NOT_SET_STRING)
  num_deriv3d = LenArray(wanted_deriv3d,NOT_SET_STRING)

  call CheckStop(num_deriv2d<1,dtxt//"num_deriv2d<1 !!")

  if(num_deriv2d > 0) then
    if(dbg0) write(*,*) dtxt//"Allocate arrays for 2d:", num_deriv2d
    allocate(f_2d(num_deriv2d),stat=alloc_err)
    call CheckStop(alloc_err,dtxt//"Allocation of f_2d")
    allocate(d_2d(num_deriv2d,LIMAX,LJMAX,LENOUT2D),stat=alloc_err)
    call CheckStop(alloc_err,dtxt//"Allocation of d_2d")
    call CheckStop(alloc_err,dtxt//"Allocation of d_3d")
    allocate(nav_2d(num_deriv2d,LENOUT2D),stat=alloc_err)
    call CheckStop(alloc_err,dtxt//"Allocation of nav_2d")
    nav_2d = 0
  end if
  if(num_deriv3d > 0) then
    if(dbg0) write(*,*) dtxt//"Allocate arrays for 3d: ", num_deriv3d
    allocate(f_3d(num_deriv3d),stat=alloc_err)
    call CheckStop(alloc_err,dtxt//"Allocation of f_3d")
    allocate(d_3d(num_deriv3d,LIMAX,LJMAX,num_lev3d,LENOUT3D),&
            stat=alloc_err)
    allocate(nav_3d(num_deriv3d,LENOUT3D),stat=alloc_err)
    call CheckStop(alloc_err,dtxt//"Allocation of nav_3d")
    nav_3d = 0
  end if

  ! Avoid hard codded IXADV_SPCS
  iadv_O3         =find_index('O3'         ,species_adv(:)%name, any_case=.true. )
  iadv_NO3_C      =find_index('NO3_c'      ,species_adv(:)%name, any_case=.true. )
  iadv_EC_C_WOOD  =find_index('EC_C_WOOD'  ,species_adv(:)%name, any_case=.true. )
  iadv_EC_C_FFUEL =find_index('EC_C_FFUEL' ,species_adv(:)%name, any_case=.true. )
  iadv_POM_C_FFUEL=find_index('POM_C_FFUEL',species_adv(:)%name, any_case=.true. )
 ! Need some special tricks for OM25, since species OM25_p and group OM25 are
 ! hard to compare due to cfac:
  iadv_OM25p      =find_index('OM25_p',     species_adv(:)%name, any_case=.true. )
  igrp_OM25       =find_index('OM25',        chemgroups(:)%name, any_case=.true.)
  ! and for deposition gradients of OM we use:
  iadv_PMf        =find_index('SO4',        species_adv(:)%name, any_case=.true. )

  ! units scaling
  ! e.g. ug_NO3_C = 1.0+e9 * MW(NO3)/MW(air)
  if(iadv_NO3_C      >0)call Units_Scale('ug',iadv_NO3_C      ,ug_NO3_C      )
  if(iadv_EC_C_WOOD  >0)call Units_Scale('ug',iadv_EC_C_WOOD  ,ug_EC_C_WOOD  )
  if(iadv_EC_C_FFUEL >0)call Units_Scale('ug',iadv_EC_C_FFUEL ,ug_EC_C_FFUEL )
  if(iadv_POM_C_FFUEL>0)call Units_Scale('ug',iadv_POM_C_FFUEL,ug_POM_C_FFUEL)

  call Define_Derived()

  iddefPMc = find_index('PMc',DDdefs(:)%name, any_case=.true.)
  !associate ( D=> DDdefs(iddefPMc) ) !does not work with gfortran
  fracPM25 = LogNormFracBelow(DDdefs(iddefPMc)%umDpgV, &
       DDdefs(iddefPMc)%sigma, 2.5, 0.001*DDdefs(iddefPMc)%rho_p)
  if(MasterProc) write(*,*) dtxt//"fracPM25 ", DDdefs(iddefPMc)%umDpgV, &
       trim(DDdefs(iddefPMc)%name), fracPM25
  !end associate ! D=> DDdefs(iddefPMc) )

!  select case(nint(DDdefs(iddefPMc)%umDpgV*10))
!    case(25);fracPM25=0.37
!    case(30);fracPM25=0.27
!end if
!    case default
!      call StopAll(dtxt//' cannot set fracPM25')
!  end select
  if ( USES%fPMc_specs(1) == 'NO3' ) fracPM25=0.0  

  if(dbg0) write(*,"(a,i4,2g12.3,i4)") dtxt//' CFAC INIT PMFRACTION Dpgv(um)',&
    iddefPMc, fracPM25, nint(10* DDdefs(iddefPMc)%umDpgV )

end subroutine Init_Derived
!=========================================================================
subroutine AddNewDeriv( name,class,subclass,txt,unit,ind,f2d,&
       dt_scale,scale, avg,iotype,Is3D)
  character(len=*), intent(in) :: name    ! e.g. DDEP_SO2_m2Conif
  character(len=*), intent(in) :: class   ! Type of data, e.g. ADV
  character(len=*), intent(in) :: subclass
  character(len=*), intent(in) :: txt     ! text where needed, e.g. "Conif"
  character(len=*), intent(in) :: unit    ! writen in netCDF output
  integer, intent(in)  :: ind    ! index in concentation array, or other
  integer, intent(in) :: f2d       ! index in f_2d arrays
  logical, intent(in) :: dt_scale  !  where scaling by dt_advec needed,
  real, intent(in)    :: scale     !  e.g. use 100.0 to get cm/s
  logical, intent(in)  :: avg      ! True => average data (divide by
                     ! nav at end),  else accumulate over run period
  character(len=*), intent(in) :: iotype  ! sets daily, monthly, etc.

  logical, intent(in), optional :: Is3D
  type(Deriv) :: inderiv

  if(trim(name)=="HMIX".and.DEBUG%DERIVED .and. MasterProc)&
     write(*,*) "ADDNEWDERIVE", iotype

  inderiv=Deriv(trim(name),trim(class),trim(subclass),&
                trim(txt),trim(unit),ind,f2d,dt_scale, scale,&
                avg,iotype)

  call AddDeriv(inderiv,Is3D=Is3D)
end subroutine AddNewDeriv
!=========================================================================
subroutine AddDeriv(inderiv,Is3D)
  type(Deriv), intent(in) :: inderiv
  logical, intent(in), optional :: Is3D
  logical :: Is3D_local

  dbg0 = (DEBUG%DERIVED .and. MasterProc )
  Is3D_local = .false.
  if(present(Is3D)) Is3D_local = Is3D

  if(Is3D_local) then
    Nadded3d = Nadded3d + 1
    N = Nadded3d
    if(dbg0) write(*,*) "Define 3d deriv ", N, trim(inderiv%name)
    call CheckStop(N>MAXDEF_DERIV3D,"Nadded3d too big! Increase MAXDEF_DERIV3D in DerivedFields_mod")
    def_3d(N) = inderiv
  else
    Nadded2d = Nadded2d + 1
    N = Nadded2d
    if(dbg0)then
      write(*,"(a,i6)") "DEBUG AddDeriv 2d ", N
      call print_Deriv_type(inderiv)
    end if
   !if(dbg0) write(*,*) "DALL", inderiv
    call CheckStop(N>MAXDEF_DERIV2D,"Nadded2d too big! Increase MAXDEF_DERIV2D in DerivedFields_mod")
    def_2d(N) = inderiv
  end if
end subroutine AddDeriv
!=========================================================================
subroutine Define_Derived()
! Set the parameters for the derived parameters, including the codes
! used by MET.NO/xfelt and scaling factors. (The scaling factors may
! be changed later in Derived_mod.
! And, Initialise the fields to zero.

  real    :: unitscale
  logical :: volunit,semivol  ! set true for volume units (e.g. ppb),group with semivol
  !FAILED logical :: outmm, outdd  ! sets time-intervals

  character(len=31) :: dname, class
  character(len=10) :: unittxt
  character(len=TXTLEN_SHORT) :: outname, outunit, outtyp, outdim, subclass
  character(len=*), parameter:: dtxt="DefDerived:"
  character(len=TXTLEN_IND)  :: outind

  integer :: ind, iadv, ishl, idebug, n, igrp, iout, isec_poll

  if(dbg0) write(6,*) " START DEFINE DERIVED "
  !   same mol.wt assumed for PPM25 and PPMCOARSE


!-- Deposition fields. Define all possible fields and their xfelt codes here:

      !code class  avg? ind scale rho Inst Yr Mn Day   name      unit

  Is3D = .false.
  dbgP = ( DEBUG%DERIVED  .and. debug_proc )
if( dbgP ) write(*,*) 'DBGUREF', u_ref(debug_li,debug_lj)


! We process the various combinations of gas-species and ecosystem:
! stuff from My_Derived

  !Deriv(name, class,    subc,  txt,           unit
  !Deriv index, f2d, dt_scale, scale, avg? rho Inst Yr Mn Day atw
  ! for AOT we can use index for the threshold, usually 40
  call AddNewDeriv( "AOT40_Grid", "GRIDAOT","subclass","-", "ppb h", &
          40, -99, T, 1.0/3600.0, F,   'YM'    )
!-------------------------------------------------------------------------------
  !Deriv(name, class,    subc,  txt,           unit
  !Deriv index, f2d, dt_scale, scale, avg? rho Inst Yr Mn Day atw

  call AddNewDeriv( "PS","PSURF",  "SURF","-",   "hPa", &
          -99,  -99,  F,  1.0,  T,  'YMD'//trim(PS_needed) )

  call AddNewDeriv( "T2m","T2m",  "-","-",   "deg. C", &
               -99,  -99, F, 1.0,  T,  'YM' )
!BIDIR
!AddNewDeriv( 		name,	class,	subclass,	txt,	unit,&
!			ind,	f2d,	dt_scale,	scale, 	avg,	iotype,	Is3D)
  call AddNewDeriv( "ugXTOT ","ugXTOT",  "-","-",   "ug/m3", &
               -99,  -99,  F,  1.0,  T,   'YMD' )
  call AddNewDeriv( "ugNH3_3m ","ugNH3_3m",  "-","-",   "ug/m3", &
               -99,  -99,  F,  1.0,  T,   'YMD' )
  call AddNewDeriv( "ugXH3_3m ","ugXH3_3m",  "-","-",   "ug/m3", &
               -99,  -99,  F,  1.0,  T,   'YMD' )
!!Hazelhos, autumn 2019: Added BiDir_NHx_Emissions
!DS see Emis_BioNatNH3
!  call AddNewDeriv( "BiDir_NHx_Emissions ", "USET", "-", &
!      "NH3 emissions BiDir", "ugN/m2/h", &
!               -99,  -99,  F,   1.0, T, 'YMDI' )
!END BIDIR


  if ( USES%TLEAF_FROM_HD )  &
    call AddNewDeriv( "dTleafHd","dTleafHd",  "-","-",   "deg. C", &
               -99,  -99, F, 1.0,  T,  'YMDH' )
  if ( USES%TLEAF_FROM_RN ) & 
    call AddNewDeriv( "dTleafRn","dTleafRn",  "-","-",   "deg. C", &
               -99,  -99, F, 1.0,  T,  'YMDH' )

!CRUDE for SEI Feb 2021:
!  call AddNewDeriv( "u_ref","u_ref",  "SURF", "-",   "m/s", &
!               -99,  -99, F, 1.0,  T,  'YM' )

! OutputFields can contain both 2d and 3d specs.
! Settings for 2D and 3D are independant.

  if(MasterProc)write(*,"(a,/,4a)") HORIZ_LINE, dtxt//": Start OutputFields"

  do ind = 1, nOutputFields
    outname= trim(OutputFields(ind)%txt1)
    outunit= trim(OutputFields(ind)%txt2)   ! eg ugN, which gives unitstxt ugN/m3
    outdim = trim(OutputFields(ind)%txt3)   ! 2d or 3d or e.g. k20
    outtyp = trim(OutputFields(ind)%txt5)   ! SPEC or GROUP or MISC
    outind = OutputFields(ind)%ind    !  H, D, M - fequency of output
    subclass = '-' ! default
    Is3D = .false.

    if ( dbgP)  write(*,*) dtxt//'Field'//trim(outname),trim(outtyp),trim(trim(OutputFields(ind)%txt4))

    if(outtyp=="MISC") then ! Simple species
      iout  = -99 ! find_index( wanted_deriv2d(i), def_2d(:)%name )
      class = trim(OutputFields(ind)%txt4)
      select case(class)
      case ('Z_MID','Z','Z_BND','Zlev','dZ_BND','dZ','EmFFprof')
        iadv = -1
        unitscale=1.0
        unittxt="m"
        Is3D=.true.
      case('SIA25','PM25','PM25X','PM25_rh50','PM25X_rh50','PM10_rh50',&
           'PM25water','PMco_water','PM25_wet','PM10_wet','PM_coarse')
        iadv = -1 ! Units_Scale(iadv=-1) returns 1.0
                  ! group_calc gets the unit conversion factor from Group_Units
        call Units_Scale(outunit,iadv,unitscale,unittxt)
        Is3D=(outdim=='3d')
      ! if(dbgP) write(*,*)"FRACTION UNITSCALE "//trim(class), unitscale
      case('FLYmax6h','FLYmax6h:SPEC')   ! Fly Level, 6 hourly maximum
        iout=find_index(outname, species_adv(:)%name, any_case=.true. )
  !-- Volcanic Emission: Skipp if not found
        if(outname(1:3)=="ASH")then
          if(MasterProc.and.DEBUG%COLSRC)&
            write(*,"(A,':',A,1X,I0,':',A)")'ColumSource',trim(outtyp),iout,trim(outname)
          if(iout<1)cycle
        end if
        call CheckStop(iout<0,dtxt//"OutputFields "//trim(outtyp)//&
                              " not found "//trim(outname))
        call Units_Scale(outunit,iout,unitscale,unittxt)
        outtyp = "FLYmax6h:SPEC"
        subclass = outdim   ! flxxx-yyy: xxx to yyy 100 feet
        outname = "MAX6h_"//trim(outname)//"_"//trim(subclass)
      case('FLYmax6h:GROUP')          ! Fly Level, 6 hourly maximum
        iout=find_index(outname,chemgroups(:)%name, any_case=.true.)
  !-- Volcanic Emission: Skipp if not found
        if(outname(1:3)=="ASH")then
          if(MasterProc.and.DEBUG%COLSRC)&
            write(*,"(A,':',A,1X,I0,':',A)")'ColumSource',trim(class),iout,trim(outname)
          if(iout<1)cycle
        end if
        call CheckStop(iout<0,dtxt//"OutputFields "//trim(outtyp)//&
                              " not found "//trim(outname))
        call Units_Scale(outunit,-1,unitscale,unittxt)
        outtyp = "FLYmax6h:GROUP"
        subclass = outdim   ! flxxx-yyy: xxx to yyy 100 feet
        outname = "MAX6h_"//trim(outname)//"_"//trim(subclass)
      case('COLUMN','COLUMN:ADV','COLUMN:SPEC')
        !COL  'NO2',          'molec/cm2' ,'k20','COLUMN'   ,'MISC' ,4,
        iout=find_index(outname, species_adv(:)%name, any_case=.true. )
        call CheckStop(iout<0,dtxt//"OutputFields "//trim(outtyp)//&
                              " not found "//trim(outname))
        call Units_Scale(outunit,iout,unitscale,unittxt)
        outtyp = "COLUMN:SPEC"
        subclass = outdim   ! k20, k16...
        outname = "COLUMN_" // trim(outname) // "_" // trim(subclass)
      case('COLUMN:SHL')
        !COL   'OH',          'molec/cm2' ,'k20','COLUMN'   ,'MISC' ,4,
        iout=find_index(outname, species_shl(:)%name, any_case=.true. )
        call CheckStop(iout<0,dtxt//"OutputFields "//trim(outtyp)//&
                              " not found "//trim(outname))
        call CheckStop(outunit, 'molec/cm2',dtxt//"OutputFields "//trim(outtyp)//&
                              " unsupported unit "//trim(outunit))
        unitscale = 1e2 ! dZ[m] to dZ[cm]
        unittxt   = trim(outunit)
        outtyp = "COLUMN:SHL"
        subclass = outdim   ! k20, k16...
        outname = "COLUMN_" // trim(outname) // "_" // trim(subclass)
      case('COLUMN:GROUP')
        iout=find_index(outname,chemgroups(:)%name, any_case=.true.)
        call CheckStop(iout<0,dtxt//"OutputFields "//trim(outtyp)//&
                              " not found "//trim(outname))
        call Units_Scale(outunit,-1,unitscale,unittxt)
        outtyp = "COLUMN:GROUP"
        subclass = outdim   ! k20, k16...
        outname = "COLUMN_" // trim(outname) // "_" // trim(subclass)
      case('AOD','AOD:TOTAL','AOD:SPEC','AOD:GROUP',&
           'EXT','EXT:TOTAL','EXT:SPEC','EXT:GROUP')
        if(.not.AOD_WANTED)cycle
        ! treat 'AOD:GROUP'/'EXT:GROUP' as 'AOD'/'EXT'
        if((outname=='AOD'.and.class=='AOD:GROUP').or.&
           (outname=='EXT'.and.class=='EXT:GROUP'))&
          class=trim(outname)
        select case(class)
        case('AOD','EXT')             ! take the full aod_group
          iout=find_index('EXTINC',chemgroups_maps(:)%name, any_case=.true.)
        case('AOD:GROUP','EXT:GROUP') ! cherry pick from aod_group
          iout=find_index(outname,chemgroups(:)%name, any_case=.true.)
        case('AOD:SPEC','EXT:SPEC' )
          iout=find_index(outname,species(:)%name, any_case=.true.)
        case default
          ! should never reach this clause
          call CheckStop(dtxt//" Unsupported output for "//&
            trim(outtyp)//":"//trim(outname)//":"//trim(outdim))
        end select
        call CheckStop(iout<0,dtxt//"OutputFields%class "//trim(class)//&
                              " not found "//trim(outname))
        unitscale = 1.0
        unittxt   = trim(outunit)
        subclass  = outdim   ! 330nm .. 1020nm
        select case(class)
        case('AOD','EXT')
          outname = trim(outname)//"_"//trim(subclass)
        case('AOD:SPEC','AOD:GROUP' )
          outname = "AOD_"//trim(outname)//"_"//trim(subclass)
        case('EXT:SPEC','EXT:GROUP' )
          outname = "EXT_"//trim(outname)//"_"//trim(subclass)
        case default
          ! should never reach this clause
          call CheckStop(dtxt//" Undefined namaming for "//&
            trim(outtyp)//":"//trim(outname)//":"//trim(outdim))
        end select
        Is3D      = (class(1:3)=="EXT")
        call AOD_init("Derived:"//trim(class),wlen=trim(subclass),out3d=Is3D)
      case('VOC')
        write(unit=errmsg,fmt="(2(A,1X,4(A,','),A),:,1X)") &
          dtxt,&
! txt string became too long for errmsg. Not worth expanding errmsg for this old output
!          trim(outname),trim(outunit),'VOC','MISC',trim(outind),&
          "VOC MISC output is no longer supported, try with",&
          'NMVOC',trim(outunit),'AIR_CONCS','GROUP',trim(outind)
        call CheckStop(errmsg)
      case default
         if(outdim=='3d')Is3D=.true.
         unitscale = 1.0
         if(outunit=="ppb") unitscale = PPBINV
         unittxt=trim(outunit)
      end select

      !if(MasterProc)write(*,"(a,/,4a)") HORIZ_LINE,&
      if(MasterProc)write(*,"(4a)") &
        dtxt//":MISC "//trim(outname),outind,trim(class)

      call AddNewDeriv(outname,class,subclass,"-",trim(unittxt),&
                       iout,-99,F,unitscale,T,outind,Is3D=Is3D)

    else ! SPEC and GROUPS of specs.

      select case(outtyp)
      case("SPEC")  ! Simple species
        iadv = find_index(outname, species_adv(:)%name, any_case=.true. )
  !-- Volcanic Emission: Skipp if not found
        if(outname(1:3)=="ASH")then
          if(MasterProc.and.DEBUG%COLSRC)&
            write(*,"(A,':',A,1X,I0,':',A)")'ColumSource',trim(outtyp),iadv,trim(outname)
          if(iadv<1)cycle
        end if
        call CheckStop(iadv<0,dtxt//"OutputFields Species not found "//trim(outname))
        iout = iadv
        call Units_Scale(outunit,iadv,unitscale,unittxt,volunit)
      case("SHL")
        ishl = find_index(outname,species_shl(:)%name, any_case=.true.)
        call CheckStop(ishl<0,dtxt//"OutputFields Short lived Species not found "//trim(outname))
        if(MasterProc) &
          write(*,*)"OutputFields Short lived Species found: "//trim(outname)
        iout = ishl
        unitscale = 1.0
        unittxt = "molec/cm3" ! No PPB possibility here !!!
        volunit = .true.
      case("GROUP") ! groups of species
        igrp = find_index(outname, chemgroups(:)%name, any_case=.true. )
  !-- Volcanic Emission: Skipp if not found
        if(outname(1:3)=="ASH")then
          if(MasterProc.and.DEBUG%COLSRC)&
            write(*,"(A,':',A,1X,I0,':',A)")'ColumSource',trim(outtyp),igrp,trim(outname)
          if(igrp<1)cycle
        end if
        call CheckStop(igrp<0,dtxt//"OutputFields Group not found "//trim(outname))
        iout = igrp
        call Units_Scale(outunit,-1,unitscale,unittxt,volunit,semivol=semivol)
        ! Units_Scale(iadv=-1) returns 1.0
        ! group_calc gets the unit conversion factor from Group_Units
        if( semivol ) subclass = 'FSOA'
        if( dbgP ) write(*,"(2a)") dtxt//'FSOA GRPOM:', &
          trims( outname // ':' // outunit // ':' // subclass )
      case default
        call CheckStop(dtxt//" Unsupported OutputFields%outtyp "//&
          trim(outtyp)//":"//trim(outname)//":"//trim(outdim))
      end select

      class="MASS";if(volunit)class="PPB"   ! Assign PPB to VOL
      select case(outdim)
      case("2d","2D","SURF")
        Is3D = .false.
        class = "SURF_"//trim(class)  //"_"//trim(outtyp)
        dname = "SURF_"//trim(outunit)//"_"//trim(outname)
        call CheckStop(find_index(dname,def_2d(:)%name, any_case=.true.)>0,&
          dtxt//"OutputFields already defined output "//trim(dname))

        !if(dbg0) write(*,"(2a,2i4,4(1x,a),es10.2)") HORIZ_LINE, dtxt//"ADD",&
        if(dbg0) write(*,"(a,2i4,4(1x,a),es10.2)") dtxt//"ADD",&
          ind, iout, trim(dname),";", trim(class), outind,unitscale

      case("Local_Correct")
        Is3D = .false.
        class = "SURF_"//trim(class)  //"_"//trim(outtyp)
        dname = "SURF_LF_"//trim(outunit)//"_"//trim(outname)
        subclass = 'LocFrac_corrected'
        call CheckStop(find_index(dname,def_2d(:)%name, any_case=.true.)>0,&
          dtxt//"OutputFields already defined output "//trim(dname))

        if(dbg0) write(*,"(a,2i4,4(1x,a),es10.2)")"ADD",&
          ind, iout, trim(dname),";", trim(class), outind,unitscale

      case("3d","3D","MLEV")
        Is3D = .true.
        class = "3D_"//trim(class)  //"_"//trim(outtyp)
        dname = "D3_"//trim(outunit)//"_"//trim(outname)
        call CheckStop(find_index(dname,def_3d(:)%name, any_case=.true.)>0,&
          dtxt//"OutputFields already defined output "//trim(dname))

        ! Always print out 3D info. Good to help avoid using 3d unless really needed!
        if( MasterProc ) write(*,"(a,2i4,4(1x,a),es10.2)")"ADD 3D outputs",  &
          ind, iout, trim(dname),";", trim(class), outind,unitscale
      case default
        call CheckStop(dtxt//" Unsupported OutputFields%outdim "//&
          trim(outtyp)//":"//trim(outname)//":"//trim(outdim))
      end select

      call AddNewDeriv(dname,class,subclass,"-",trim(unittxt),&
                       iout,-99,F,unitscale,T,outind,Is3D=Is3D)
    end if
  end do ! OutputFields
  if(MasterProc)write(*,"(a,/,4a)") HORIZ_LINE, dtxt//": End OutputFields"

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  do n = 1, nOutputMisc
    Is3D=(OutputMisc(n)%class=="MET3D").or.(OutputMisc(n)%name(1:2)=='D3')&
         .or.(OutputMisc(n)%subclass(1:2)=='D3')
    if(MasterProc) write(*,"(3(A,1X),L1)") &
      dtxt//'ADDMISC',trim(OutputMisc(n)%name),'Is3D',Is3D
    call AddDeriv(OutputMisc(n),Is3D=Is3D)
  end do

!-------------------------------------------------------------------------------
  do n = 1, nMosaic
    if ( dbg0 ) write(*,*) dtxt//"DEBUG MOSAIC AddDeriv ", n, MosaicOutput(n)
    call AddDeriv( MosaicOutput(n) )
  end do
!-------------------------------------------------------------------------------
! Areas of deposition-related ecosystems. Set externally
  do n = 1, NDEF_ECOSYSTEMS
     if(dbg0) write(*,*) dtxt//"ECODEF ",n, trim( DepEcoSystem(n)%name )
     call AddDeriv( DepEcoSystem(n) )
  end do
!!-------------------------------------------------------------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  do ind = 1, size(WDEP_WANTED(1:nOutputWdep)%txt1)
    dname = "WDEP_"//trim(WDEP_WANTED(ind)%txt1)
    outind = trim(WDEP_WANTED(ind)%txt4)
    select case(WDEP_WANTED(ind)%txt2)
    case("PREC")
      call AddNewDeriv("WDEP_PREC","PREC ","-","-", "mm",  &
                        -1, -99,   F,    1.0,   F,    outind)
    case("SPEC")
      iadv = find_index(WDEP_WANTED(ind)%txt1, species_adv(:)%name, any_case=.true.)
      call CheckStop(iadv<1, "WDEP_WANTED Species not found " // trim(dname) )

      call Units_Scale(WDEP_WANTED(ind)%txt3,iadv,unitscale,unittxt)
      call AddNewDeriv( dname, "WDEP", "-", "-", unittxt , &
              iadv, -99,   F, unitscale,     F,  outind)
    case("GROUP")
      igrp = find_index(dname, chemgroups(:)%name, any_case=.true.)
      call CheckStop(igrp<1, "WDEP_WANTED Group not found " // trim(dname) )

      ! Just get units text here.
      ! Init_WetDep gets the unit conversion factors from Group_Scale.
      call Units_Scale(WDEP_WANTED(ind)%txt3,-1,unitscale,unittxt)
      call AddNewDeriv( dname,  "WDEP ","-","-", unittxt ,  &
              igrp, -99,   F,      1.0,   F,     outind)
    case default
      call CheckStop("Unknown WDEP_WANTED type " // trim(WDEP_WANTED(ind)%txt2) )
    end select
    if(MasterProc) write(*,*)dtxt//"Wet deposition output: ",&
       trim(dname)," ",trim(unittxt)
  end do

!Emissions:
! We use mg/m2 outputs for consistency with depositions
! Would need to multiply by GridArea_m2 later to get ktonne/grid, but not
! done here.
!
! BVOC called every dt_advec, so use dt_scale=1.0e6 to get from kg/m2/s to
!  mg/m2 accumulated (after multiplication by dt_advec)

    ! AddNewDeriv( name,class,subclass,txt,unit,
    !    ind,f2d, dt_scale,scale, avg,iotype,Is3D)

  do  ind = 1, NEMIS_BioNat
    if(EMIS_BioNat(ind)(1:5)=="ASH_L")cycle   ! skip ASH_LxxByy for AshInversion
    dname = "Emis_mgm2_BioNat" // trim(EMIS_BioNat(ind) )
    if(MasterProc) write(*,'(a,i4,a)') dtxt//'NatEmis ', ind, trim(dname)
    call AddNewDeriv( dname, "NatEmis", "-", "-", "mg/m2", &
                 ind , -99, T ,    1.0e6,     F, 'YM' )
  end do

! SNAP emissions called every hour, given in kg/m2/s, but added to
! d_2d every advection step, so get kg/m2.
! Need 1.0e6 to get from kg/m2 to mg/m2 accumulated.
!
! Future option - might make use of Emis_Molwt to get mg(N)/m2
  do  ind = 1, size(EMIS_FILE)
    dname = "Emis_mgm2_" // trim(EMIS_FILE(ind))
    if (HourlyEmisOut .and. DailyEmisOut) then
       call AddNewDeriv( dname, "TotEmis", "-", "-", "mg/m2", &
            ind , -99, T,  1.0e6,  F,  'YMDH' )
    else if (DailyEmisOut) then
       call AddNewDeriv( dname, "TotEmis", "-", "-", "mg/m2", &
            ind , -99, T,  1.0e6,  F,  'YMD' )
    else  if (HourlyEmisOut) then
       call AddNewDeriv( dname, "TotEmis", "-", "-", "mg/m2", &
            ind , -99, T,  1.0e6,  F,  'YMH' )
    else
       call AddNewDeriv( dname, "TotEmis", "-", "-", "mg/m2", &
            ind , -99, T,  1.0e6,  F,  'YM' )
    endif
  end do ! ind

  isec_poll = 0
  do  i = 1, NEMIS_File
     dname = "Sec_Emis_mgm2_"//trim(EMIS_FILE(i))
     isec_poll = i
     if (HourlyEmisOut .and. DailyEmisOut) then
        call AddNewDeriv( dname, "SecEmis", "-", "-", "mg/m2", &
             isec_poll , -99, T,  1.0e6,  F,  'YMDH' )
     else if (DailyEmisOut) then
        call AddNewDeriv( dname, "SecEmis", "-", "-", "mg/m2", &
             isec_poll , -99, T,  1.0e6,  F,  'YMD' )
     else if (HourlyEmisOut) then
        call AddNewDeriv( dname, "SecEmis", "-", "-", "mg/m2", &
             isec_poll , -99, T,  1.0e6,  F,  'YMH' )
     else
        call AddNewDeriv( dname, "SecEmis", "-", "-", "mg/m2", &
             isec_poll , -99, T,  1.0e6,  F,  'YM' )
     endif
  enddo
  do isec=1,NSECTORS
     if(SecEmisOutWanted(isec))then
        do  i = 1, NEMIS_File
           write(dname,"(A)")trim(SECTORS(isec)%longname)//"_Emis_mgm2_"//trim(EMIS_FILE(i))
           isec_poll = (isec)*(NEMIS_File) + i
           if (HourlyEmisOut .and. DailyEmisOut) then
              call AddNewDeriv( dname, "SecEmis", "-", "-", "mg/m2", &
                   isec_poll , -99, T,  1.0e6,  F,  'YMDH' )
           else if (DailyEmisOut) then
              call AddNewDeriv( dname, "SecEmis", "-", "-", "mg/m2", &
                   isec_poll , -99, T,  1.0e6,  F,  'YMD' )
           else if (HourlyEmisOut) then
              call AddNewDeriv( dname, "SecEmis", "-", "-", "mg/m2", &
                   isec_poll , -99, T,  1.0e6,  F,  'YMH' )
           else
              call AddNewDeriv( dname, "SecEmis", "-", "-", "mg/m2", &
                   isec_poll , -99, T,  1.0e6,  F,  'YM' )
           endif
        end do
     endif
  end do

  if(USES%OCEAN_DMS)then
    dname = "Emis_mgm2_DMS"
    call AddNewDeriv( dname, "Emis_mgm2_DMS", "-", "-", "mg/m2", &
                       ind , -99, T,  1.0,  F,  'YM' )
  end if
  if(USES%OCEAN_NH3)then
     dname = "Emis_mgm2_Ocean_NH3"
    call AddNewDeriv( dname, "Emis_mgm2_Ocean_NH3", "-", "-", "mg/m2", &
                       ind , -99, T,  1.0,  F,  'YM' )
  end if

!Splitted total emissions (inclusive Natural)
  do ind=1,nrcemis
    dname = "EmisSplit_mgm2_"//trim(species(iqrc2itot(ind))%name)
    call AddNewDeriv(dname, "EmisSplit_mgm2", "-", "-", "mg/m2", &
                        ind , -99, T, 1.0e6,   F,  'YM' )
  end do

  if(find_index("SURF_PM25water",def_2d(:)%name,any_case=.true.)<1)&
  call AddNewDeriv("SURF_PM25water", "PM25water", "-", "-","ug/m3", &
                       -99 , -99, F, 1.0,   T,  'YMD' )

  if(AERO%JUN21AERO .and. &
     find_index("SURF_PMco_water",def_2d(:)%name,any_case=.true.)<1)&
  call AddNewDeriv("SURF_PMco_water", "PMco_water", "-", "-","ug/m3", &
                       -99 , -99, F, 1.0,   T,  'YMD' )

! As for GRIDAOT, we can use index for the threshold
  call AddNewDeriv( "SOMO35","SOMO",  "SURF","-",   "ppb.day", &
                  35, -99, F, 1.0,   F,   'YM' )
  call AddNewDeriv( "SOMO0 ","SOMO",  "SURF","-",   "ppb.day", &
                  0 , -99, F, 1.0,   F,   'YM' )
  if(iadv_o3>0) &
  call AddNewDeriv( "SURF_MAXO3","MAXADV", "O3","-",   "ppb", &
           iadv_o3, -99, F, PPBINV,   F,   'YMD')

!-- 3-D fields

Is3D = .true.
  do ind = 1, size(D3_OTHER)
    select case ( trim(D3_OTHER(ind)) )
    case ("D3_PM25water")
      if(find_index("D3_PM25water",def_3d(:)%name,any_case=.true.)<1)&
      call AddNewDeriv("D3_PM25water", "PM25water", "-", "-","ug/m3", &
         -99, -99, F, 1.0,   T,  'YM',    Is3D ) !

    case ("D3_m_TH")
      call AddNewDeriv("D3_m_TH","TH", "-","-",   "m", &
         -99, -99, F,  1.0,  F,  'YM',     Is3D )

    case ("D3_m2s_Kz")
      call AddNewDeriv( "D3_Kz","Kz", "-","-",   "-", &
         -99, -99, F,  1.0,  F,  'YM',     Is3D )

    case ("D3_T")
      call AddNewDeriv("D3_T","T", "-","-",   "K", &
         -99, -99, F,  1.0,  T,  'YM',     Is3D )

    case ("D3_Zmid")
      if(find_index("D3_Zmid",def_3d(:)%name,any_case=.true.)<1)&
      call AddNewDeriv("D3_Zmid", "Z_MID", "-", "-", "m", &
                      -99 , -99, F, 1.0,   T, 'YMD',    Is3D  )

    case ("EmFFprof")
      if(find_index("EmFFprof",def_3d(:)%name,any_case=.true.)<1)&
      call AddNewDeriv("EmFFprof", "EmFFprof", "-", "-", "m", &
                      -99 , -99, F, 1.0,   T, 'YMD',    Is3D  )

     case ("D3_Zlev")
      if(find_index("D3_Zlev",def_3d(:)%name,any_case=.true.)<1)&
      call AddNewDeriv("D3_Zlev", "Z_BND", "-", "-", "m", &
           -99 , -99, F, 1.0,   T, 'YMD',    Is3D  )

!    case ("D3_EmF_CO_")
!      call AddNewDeriv("D3_EmF_CO","CO_Fire_Prof", "-","-",   "kg/m2", &
!         -99, -99, F,  1.0,  F,  'YMD',     Is3D )
 


    end select
  end do

  ! Get indices of wanted fields in larger def_xx arrays:
  do i = 1, num_deriv2d
    if(dbg0) write(*,*)dtxt//"CHECK2d",num_deriv2d,i,trim(wanted_deriv2d(i))
    if(MasterProc) call CheckStop(count(f_2d(:i)%name==wanted_deriv2d(i))>0,&
        dtxt//"REQUESTED 2D DERIVED ALREADY DEFINED: "//trim(wanted_deriv2d(i)))
    ind = find_index( wanted_deriv2d(i), def_2d(:)%name,any_case=.false. )
    if(ind<=0)then
       !try with case insensitive
    ind = find_index( wanted_deriv2d(i), def_2d(:)%name,any_case=.true. )
    endif
    if(ind>0)then
      f_2d(i) = def_2d(ind)
      if(dbg0) write(*,"(2(a,i4),3(1x,a),2i4)") "Index f_2d ",i,  &
        " = def:set "//trim(wanted_deriv2d(i)),ind, &
        trim(def_2d(ind)%name),trim(def_2d(ind)%unit),&
        trim(def_2d(ind)%class), f_2d(i)%index,  f_2d(i)%f2d
    elseif(MasterProc)then
      print *,dtxt//"D2IND OOOPS wanted_deriv2d not found: ", wanted_deriv2d(i)
      print *,dtxt//"OOOPS N,N :", num_deriv2d, Nadded2d
      print "(a,i4,a)",(dtxt//"Had def_2d: ",idebug,&
        trim(def_2d(idebug)%name),idebug = 1, Nadded2d)
      call CheckStop(dtxt//"OOPS1 STOPPED" // trim( wanted_deriv2d(i) ) )
    end if
  end do

  do i = 1, num_deriv3d
    if(dbg0) write(*,*) dtxt//"CHECK3d",num_deriv3d,i,trim(wanted_deriv3d(i))
    if(MasterProc)&
     call CheckStop(count(f_3d(:i)%name==wanted_deriv3d(i))>0,&
      dtxt//"REQUESTED 3D DERIVED ALREADY DEFINED: "//trim(wanted_deriv3d(i)))
    ind = find_index( wanted_deriv3d(i), def_3d(:)%name,any_case=.true. )
    if(ind>0)then
      f_3d(i) = def_3d(ind)
      if(dbg0) print "(2(a,i4),3(1x,a))","Index f_3d ",i,  &
        " = def ",ind,trim(def_3d(ind)%name),trim(def_3d(ind)%unit),&
        trim(def_3d(ind)%class)
    elseif(MasterProc)then
      print *,"D3IND OOOPS wanted_deriv3d not found: ", wanted_deriv3d(i)
      print *,"OOOPS N,N :", num_deriv3d, Nadded3d
      print "(a,i4,a)",("Had def_3d: ",idebug,&
        trim(def_3d(idebug)%name),idebug = 1, Nadded3d)
      call CheckStop(dtxt//"OOPS STOPPED" // trim( wanted_deriv3d(i) ) )
    end if
  end do

  !Initialise to zero
  if (num_deriv2d > 0) d_2d(:,:,:,:) = 0.0
  if (num_deriv3d > 0) d_3d(:,:,:,:,:) = 0.0

  ! Determine actual output time ranges for Wanted output
  iou_list(:)=.false.
  do iou=IOU_MIN,IOU_MAX
    do i=1,num_deriv2d
      if(iou_list(iou))exit
      iou_list(iou)=(index(f_2d(i)%iotype,IOU_KEY(iou))>0)
    end do
    do i=1,num_deriv3d
      if(iou_list(iou))exit
      iou_list(iou)=(index(f_3d(i)%iotype,IOU_KEY(iou))>0)
    end do
  end do

  VGtest_out_ix = 0
  if(allocated(f_2d)) &
       VGtest_out_ix = find_index("VgRatio", f_2d(:)%subclass)
  if(VGtest_out_ix>0 .and. me==0)then
     write(*,*)'will output Vg Ratio'
  else
     if(me==0)  write(*,*)'will NOT output Vg Ratio',&
       find_index("VgRatio", f_2d(:)%subclass),&
       find_index("logz0", f_2d(:)%subclass)
  endif


  if(SOURCE_RECEPTOR)&            ! We include daily and monthly also
    iou_list(IOU_DAY+1:)=.false.  ! for SOURCE_RECEPTOR mode which makes
                                  ! it easy for debugging

  if(USES%SKIP_INCOMPLETE_OUTPUT)then             ! reduce output
    ! skip daily/monthly/full-run for runs under 1/28/181 days
    select case(INT(tdif_days(to_stamp(startdate),to_stamp(enddate))))
      case(      0);iou_list(:IOU_HOUR-1)=.false. ! Only hourly outputs
      case(  1: 27);iou_list(: IOU_DAY-1)=.false. ! .. and dayly
      case( 28:180);iou_list(: IOU_MON-1)=.false. ! .. and monthly
      case(181:   );                              ! .. and full-run
    end select
  end if

  if(dbgP) write(*,"(A,': ',10(I2,A2,L2,:,','))")"Wanted IOUs",&
    (iou,IOU_KEY(iou),iou_list(iou),iou=IOU_MIN,IOU_MAX)
end subroutine Define_Derived
!=========================================================================
function wanted_iou(iou,iotype,only_iou) result(wanted)
  integer, intent(in)                 :: iou
  character(len=*),intent(in),optional:: iotype
  integer         ,intent(in),optional:: only_iou
  logical                             :: wanted
  wanted=(iou>=IOU_MIN).and.(iou<=IOU_MAX)  ! in range ov valid IOUs?
  if(wanted)wanted=iou_list(iou)       ! any output requires iou?
  if(wanted.and.present(iotype))then
    wanted=(index(iotype,IOU_KEY(iou))>0)   ! iotype contains IOU_KEY(iou)?
  end if
  if(wanted.and.present(only_iou))then
    wanted=(iou==only_iou)                  ! is only_iou?
  end if
end function wanted_iou
!=========================================================================
subroutine Derived(dt,End_of_Day,ONLY_IOU)
!*** DESCRIPTION
!  Integration and averaging of chemical fields.
!  Includes AOT40, AOT60 if present

  real, intent(in)    :: dt                   ! time-step used in intergrations
  logical, intent(in) :: End_of_Day           ! e.g. 6am for EMEP sites
  integer, intent(in), optional :: ONLY_IOU   ! IOU_INST update only instantenous fields,
                                              ! IOU_HOUR update (mean) hourly and inst. fields.

  character(len=len(f_2d%name)) :: name  !  See defs of f_2d
  character(len=len(f_2d%class)) :: class  !  See defs of f_2d
  character(len=len(f_2d%subclass)) :: subclass  !  See defs of f_2d
  character(len=TXTLEN_SHORT)    :: txt2
  real :: thour                          ! Time of day (GMT)
  real :: timefrac                       ! dt as fraction of hour (3600/dt)
  real :: dayfrac              ! fraction of day elapsed (in middle of dt)
  real :: af, fl0, fl1
  real, save :: km2_grid
  integer :: ntime                        ! 1...NTDAYS
  integer :: klow                         ! lowest extent of column data
  real, dimension(LIMAX,LJMAX) :: density ! roa[kgair m-3] when scale in ug, else 1
  real, dimension(LIMAX,LJMAX) :: tmpwork
  logical, dimension(LIMAX,LJMAX) :: mask2d
  real, dimension(LIMAX,LJMAX,KMAX_MID) :: inv_air_density3D
            ! Inverse of No. air mols/cm3 = 1/M
            ! where M =  roa (kgair m-3) * to_molec_cm3  when ! scale in ug,  else 1
  logical, save :: first_call = .true.
  integer :: igrp, ngrp   ! group methods
  logical :: needroa
  integer, save :: &      ! needed for PM25*,PM10*
    ind2d_sia=-999 ,ind3d_sia=-999,   &
    ind2d_pmfine=-999 ,ind3d_pmfine=-999,   &
    ind2d_pmwater=-999,ind3d_pmwater=-999,  &
    ind2d_pm10=-999   ,ind3d_pm10=-999,     &
    ind2d_pm25=-999   ,ind2d_pmcowater=-999

  integer :: imet_tmp, ind, iadvDep
  real, pointer, dimension(:,:,:) :: met_p => null()

  logical, allocatable, dimension(:)   :: ingrp
  integer :: wlen,ispc,kmax,iem, nerr=0
  integer :: isec_poll,isec,iisec,ii,ipoll,itemp
  real :: default_frac,tot_frac,loc_frac_corr
  character(len=*), parameter :: dtxt='Deriv:'
  real pp, temp, qsat
  real, save, allocatable :: D8M(:,:,:,:), D8Max(:,:,:), hourM(:,:,:),D8_26Max(:,:,:)

  logical, save :: make_MaxD8M_nth = .false.
  integer , save :: i_MaxD8M_26th = 0 !index in f_2d
  integer , save :: i_MaxD8M_1st = 0 !index in f_2d
  integer , save :: i_MaxD8M_2nd = 0 !index in f_2d
  integer , save :: i_MaxD8M_3rd = 0 !index in f_2d
  integer , save :: i_MaxD8M_4th = 0 !index in f_2d
  integer :: i_26th, i4th, i3rd, i2nd, i1st
  integer, save :: count_AvgMDA8_m=0,count_AvgMDA8_y=0
  integer, save :: count_AvgMDA8AprSep_m=0,count_AvgMDA8AprSep_y=0
  real :: w_m,w_y !weights

  if(.not. date_is_reached(spinup_enddate))return ! we do not average during spinup

  timefrac = dt/3600.0
  thour = current_date%hour+current_date%seconds/3600.0

  daynumber=day_of_year(current_date%year,current_date%month,&
                        current_date%day)

  dbgP0 = dbgP .and. first_call 

  ! Just calculate once, and use where needed
  forall(i=1:limax,j=1:ljmax) density(i,j) = roa(i,j,KMAX_MID,1)

  !****** 2-D fields **************************

  if ( ORGANIC_AEROSOLS ) then

     ! Start with OM25_p as this is used  by several other groups, and may have
     ! diverged from OM25_GROUP due to advection etc.
     ! Make sure OM25_p is still equal to the sum of its components
     ! after advection routines:

     call Reset_3dOrganicAerosol(debug_proc) ! -> revised OM25_p
     call CheckStop(iadv_OM25p<1, dtxt//'Need OM25p spec for OrgAer')
     call CheckStop(igrp_OM25<1,  dtxt//'Need OM25 group for OrgAer')
     if( dbgP0) write(*,*) dtxt//"OrgIndices:", iadv_OM25p, igrp_OM25
  end if

  do n = 1, num_deriv2d

    class = f_2d(n)%class
    subclass = f_2d(n)%subclass
    name  = f_2d(n)%name
    ind   = f_2d(n)%index

    if( dbgP0 ) &
       write(*,"(a,i3,9a)") "Derive2d-name-class", n, " C:", trim(class), &
            " U:", trim(f_2d(n)%unit), " N:", trim(name), ":END"


    !*** user-defined time-averaging. Here we have defined TADV and TVOC
    !    so that 8-hour daytime averages will be calculated.
    !    Just comment out if not wanted, or (better!) don't define any
    !    f_2d as TADV or TVOC

    if ( class == "TADV" .or. class == "TVOC" ) then
      if(thour <= 8.0 .or. thour > 16.0 ) cycle  ! Start next species
    end if

    ! hmix average at 00 and 12:
    if ( class == "HMIX00" .or. class == "XKSIG00" ) then
      if(thour /= 0.0 ) cycle  ! Start next species
    end if

    if ( class == "HMIX12" .or. class == "XKSIG12" ) then
      if(thour /= 12.0 ) cycle  ! Start next species
    end if

    !if ( DEBUG .and. MasterProc .and. first_call ) then
    if(MasterProc.and.first_call)&
      write(*,"(a,i4,1x,a,1x,i0,1x,a)") "1st call Derived 2d", n, &
        trim(name), ind, trim(class)

    select case ( class )

    case ( "MET2D")

     ! Meteo fields are available through their names and a pointer, either
     ! from the read-in NWP fields (met%) or the derived met fields
     ! (metderiv%), see MetFields_mod. We thus use the required name and see
     ! if we can find it in either met% or metderiv%

      imet_tmp = find_index(subclass, met(:)%name ) ! subclass has meteo name from MetFields
      if( imet_tmp > 0 ) then
        !Note: must write bounds explicitly for "special2d" to work
         if(met(imet_tmp)%dim==3) then
            met_p => met(imet_tmp)%field(1:limax,1:ljmax,KMAX_MID:KMAX_MID,1) !take lowest level
         else 
            met_p => met(imet_tmp)%field(1:limax,1:ljmax,1:1,1)
         end if
      else
        imet_tmp = find_index(subclass, derivmet(:)%name )
        if ( imet_tmp > 0 )then
          if(met(imet_tmp)%dim==3) then
            met_p => derivmet(imet_tmp)%field(1:limax,1:ljmax,KMAX_MID:KMAX_MID,1) !take lowest level
          else 
            met_p => derivmet(imet_tmp)%field(1:limax,1:ljmax,1:1,1)
          end if
        end if
      end if

      if( imet_tmp > 0 ) then
         kmax=1
         if( MasterProc.and.first_call) write(*,*) "MET2D"//trim(name), &
              imet_tmp, met_p(1,1,kmax),loc(met(imet_tmp)%field(1,1,1,1))
         forall ( i=1:limax, j=1:ljmax )
            d_2d( n, i,j,IOU_INST) = met_p(i,j,kmax)
         end forall

         met_p => null()

      else ! Not found!
        !make derived fields:
        select case ( subclass )
        case ("inv_u10")
           forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = 1.0/(0.2+ws_10m(i,j,1))
           end forall
 
        case default
           if( first_call)  then
              if( MasterProc) write(*,*) "MET2D NOT FOUND"//trim(name)//":"//trim(subclass)
              forall ( i=1:limax, j=1:ljmax )
                 d_2d( n, i,j,IOU_INST) = 0.0 ! UNDEF_R
              end forall
           end if
        end select
     end if
      
    ! The following can be deleted once testing of MET2D is finished...
    case ( "xm_i" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = xm_i(i,j)
      end forall
    case ( "lon" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = glon(i,j)
      end forall
    case ( "xm_j" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = xm_j(i,j)
      end forall
    case ( "lat" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = glat(i,j)
      end forall
    case ( "Kz_m2s" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = Kz_m2s(i,j,KMAX_BND-1)
      end forall
    case ( "ws_10m" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = ws_10m(i,j,1)
    end forall
    case ( "rh2m" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = rh2m(i,j,1)
    end forall

    case ( "SurfAreaPMF_um2cm3" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SurfArea_um2cm3(AERO%PM_F,i,j)
      end forall
      if ( dbgP ) call write_debug(n,ind, "SurfArea_NSDF")
    case ( "SurfAreaPM_um2cm3" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SurfArea_um2cm3(AERO%PM,i,j)
      end forall
      if ( dbgP ) call write_debug(n,ind, "SurfArea_NSDF")
    case ( "SurfAreaSSF_um2cm3" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SurfArea_um2cm3(AERO%SS_F,i,j)
      end forall
      if ( dbgP ) call write_debug(n,ind, "SurfArea_SSF")
    case ( "SurfAreaSSC_um2cm3" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SurfArea_um2cm3(AERO%SS_C,i,j)
      end forall
    case ( "SurfAreaDUF_um2cm3" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SurfArea_um2cm3(AERO%DU_F,i,j)
      end forall
    case ( "SurfAreaDUC_um2cm3" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SurfArea_um2cm3(AERO%DU_C,i,j)
      end forall

    case ( "u_ref" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = u_ref(i,j)
      end forall
    if ( dbgP ) call write_debug(n,ind, "PUREF")

    case ( "SMI_deep" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SoilWater_deep(i,j,1)
      end forall
    if ( dbgP ) call write_debug(n,ind, "SoilWater_DEEP")
    case ( "SMI_uppr" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SoilWater_uppr(i,j,1)
      end forall
    if ( dbgP ) call write_debug(n,ind, "SoilWater_uppr")

    case ( "T2m" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = t2_nwp(i,j,1) - 273.15
      end forall
!BIDIR
     case( "ugXTOT", "ugNH3_3m", "ugXH3_3m") !, "BiDir_NHx_Emissions" )
      if ( USES%BIDIR ) then
       call BiDir_Derived(class,n,limax,ljmax,nerr)
       if (debug_proc) call write_debug(n,ind, class )
       call CheckStop(nerr>0,'BIDIR ERROR'//trim(class) )
      end if
!END BIDIR

    case ( "dTleafHd" )
      if ( USES%TLEAF_FROM_HD ) then
        forall ( i=1:limax, j=1:ljmax )
          d_2d( n, i,j,IOU_INST) = dTleafHd(i,j)
        end forall
        if ( dbgP ) call write_debug(n,ind, "dTleafHd")
      end if
    case ( "dTleafRn" )
      if ( USES%TLEAF_FROM_RN ) then
        forall ( i=1:limax, j=1:ljmax )
          d_2d( n, i,j,IOU_INST) = dTleafRn(i,j)
        end forall
      end if

    case ( "XSNOW" ) ! Was not snow depth, but rather flag
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = Grid_snow(i,j)
      end forall

    case ( "TotDobs" ) ! Dobson populated in CloudJ_mod.f90
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = Dobson(i,j)
      end forall

    case ( "pH_surface" ) ! surface pH from AerosolCalls
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = pH(i,j)
      end forall

    case ( "SNratio" )
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = min(3.0,so2nh3_24hr(i,j))
      end forall

    case ( "PSURF" )

      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = ps(i,j,1)*0.01
        !NOT YET - keep hPa in sites:d_2d( n, i,j,IOU_INST) = ps(i,j,1)
      end forall
      if ( dbgP ) call write_debug(n,ind, "PPSURF")

   case ( "phih" )
       ind = f_2d(n)%index !height 
       do j=1,ljmax
          do i=1,limax
             if (invl_nwp(i,j) < 0) then
                !As in Garratt and Obrien for phih, so with Prandtl number
                d_2d( n, i,j,IOU_INST) = 1.0/sqrt(1-16.*ind*invl_nwp(i,j)) 
             else
                !As in Garratt, Prandtl number is 1 in stable boundary layer
                d_2d( n, i,j,IOU_INST) = 1+5.*ind*invl_nwp(i,j) 
             end if
          end do
       end do
       
    case( "CloudFrac" )
      forall ( i=1:limax, j=1:ljmax )
         d_2d( n, i,j,IOU_INST) = cc3dmax(i,j,KMAX_MID) 
      end forall
    if ( dbgP ) call write_debug(n,ind, "CloudFrac")


    case ( "HMIX", "HMIX00", "HMIX12" )

      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = hmix(i,j,1)
      end forall

      if ( dbgP ) then
       write(*,fmt="(a12,i4,f12.3)") "HMIX" , n , &
               d_2d(n,debug_li,debug_lj,IOU_INST)
      end if

    case ( "SURF_PPB_SPEC" )

!Do not delete
!Not upgraded to new format yet
!      if(subclass=='LocFrac_corrected')then
!         do ipoll=1,uEMEP%Npoll        
!            do i=1,uEMEP%poll(ipoll)%Nix
!               if(ind==uEMEP%poll(ipoll)%ix(i))goto 44
!            enddo
!         enddo
!         if(me==0)write(*,*)'WARNING, no local fractions found for ',trim(class),' index ',ind
!         44 continue
!         if(me==0.and. first_call)then
!            write(*,*)'local fractions found for ',trim(class),&
!              ' SPCex ',ind,' name ',trim(species_adv(ind)%name),&
!              ' locfrac pollutant ',trim(uEMEP%poll(ipoll)%emis)
!            do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!             write(*,*)'local correction for sector',isec,' pollutant ',trim(species_adv(ind)%name)
!         enddo
!
!         endif
!         do j=1,ljmax
!         do i=1,limax
!         default_frac=0.0!local, but any sector that is not explicit
!         tot_frac=0.0!all local (any sector)
!         loc_frac_corr=0.0
!         !isec is sector (number between 1 and 11)
!         !iisec is index over available sectors
!         !isec_poll is an internal uEMEP index that is a combination of sector and pollutant indices 
!         do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!            if(isec/=0)then
!               default_frac = default_frac - loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!               tot_frac =  tot_frac + loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!            endif
!            if(isec==0)default_frac = default_frac + loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!         enddo
!         default_frac=max(0.0,default_frac)!in case "sec=0" not available
!         tot_frac = tot_frac + default_frac
!         do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!            if(isec/=0)then
!               loc_frac_corr=loc_frac_corr+loc_frac(isec_poll,0,0,i,j,KMAX_MID)*2!*LocEmisFac(isec)
!            endif
!         enddo
!         loc_frac_corr=loc_frac_corr+default_frac*2!*LocEmisFac_default
!         loc_frac_corr=loc_frac_corr+(1-tot_frac)! weight one for pollutants from nonlocal sources
!         
!               d_2d( n, i,j,IOU_INST) = xn_adv(ind,i,j,KMAX_MID) &
!                     * cfac(ind,i,j) * loc_frac_corr !NB: cfac also for local fractions (correct?)
!         enddo
!         enddo
!      else
         forall ( i=1:limax, j=1:ljmax )
            d_2d( n, i,j,IOU_INST) = xn_adv(ind,i,j,KMAX_MID) &
                 * cfac(ind,i,j)
         end forall
!      endif
      if ( dbgP ) call write_debugadv(n,ind, 1.0, "PPB OUTS")
      
    case ( "SURF_MASS_SPEC" )  ! Here we need density

       !Do not delete
       !Not upgraded to new format yet
!      if(subclass=='LocFrac_corrected')then
!         do ipoll=1,uEMEP%Npoll        
!            do i=1,uEMEP%poll(ipoll)%Nix
!               if(ind==uEMEP%poll(ipoll)%ix(i))goto 45
!            enddo
!         enddo
!         if(me==0)write(*,*)'WARNING, no local fractions found for ',trim(class),' index ',ind
!         45 continue
!         if(me==0.and. first_call)then
!            write(*,*)'local fractions found for ',trim(class),' index ',ind,&
!              ' name ',trim(species_adv(ind)%name),&
!              ' locfrac pollutant ',trim(uEMEP%poll(ipoll)%emis)
!            do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!             write(*,*)'local correction for sector',isec,' pollutant ',trim(species_adv(ind)%name)
!         enddo
!
!         endif
!         do j=1,ljmax
!         do i=1,limax
!         default_frac=0.0!local, but any sector that is not explicit
!         tot_frac=0.0!all local (any sector)
!         loc_frac_corr=0.0
!         !isec is sector (number between 1 and 11)
!         !iisec is index over available sectors
!         !isec_poll is an internal uEMEP index that is a combination of sector and pollutant indices 
!         do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!            if(isec/=0)then
!               default_frac = default_frac - loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!               tot_frac =  tot_frac + loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!            endif
!            if(isec==0)default_frac = default_frac + loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!         enddo
!         default_frac=max(0.0,default_frac)!in case "sec=0" not available
!         tot_frac = tot_frac + default_frac
!         do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!            if(isec/=0)then
!               loc_frac_corr=loc_frac_corr+loc_frac(isec_poll,0,0,i,j,KMAX_MID)*2!*LocEmisFac(isec)
!            endif
!         enddo
!         loc_frac_corr=loc_frac_corr+default_frac*2!*LocEmisFac_default
!         loc_frac_corr=loc_frac_corr+(1-tot_frac)!No correction for pollutants from other sources
!         
!               d_2d( n, i,j,IOU_INST) = xn_adv(ind,i,j,KMAX_MID) &
!                    * cfac(ind,i,j)  * density(i,j)* loc_frac_corr
!         enddo
!         enddo
!      else
       iadvDep = ind
       if( ind == iadv_OM25p ) iadvDep=iadv_PMf

       forall ( i=1:limax, j=1:ljmax )
           d_2d( n, i,j,IOU_INST) = xn_adv(ind,i,j,KMAX_MID) &
                               * cfac(iadvDep,i,j) * density(i,j)
       end forall
       if( dbgP .and. ind== iadv_OM25p) then
         i=debug_li
         j=debug_lj
         write(*,'(2a,4es14.5)') dtxt//'IOM_SPC ', f_2d(n)%name, &
             d_2d(n,i,j,IOU_INST)*f_2d(n)%scale, &
             cfac(iadvDep,i,j), density(i,j),f_2d(n)%scale
       end if

      if ( dbgP ) call write_debugadv(n,ind, &
                density(debug_li,debug_lj), "SURF_MASS",indDep=iadvDep)
!   case ( "SURF_molec_SHL" )  ! short-lived. Follows pattern of MAXSHL below
!
!         forall ( i=1:limax, j=1:ljmax )
!           d_2d( n, i,j,IOU_INST) = xn_shl(ind,i,j,KMAX_MID)
!         end forall
!      if ( dbgP ) write(*,'(a,f8.2,3es12.3)') &
!          'SHLSHLmcc'//trim( species(ind)%name), thour, &
!           xn_shl(ind,debug_li,debug_lj,KMAX_MID), density(debug_li,debug_lj), to_molec_cm3
!
!      endif
   ! WARNING CLASS PPB just means volume based..
    case ( "SURF_PPB_SHL" )        ! short-lived. Follows pattern of MAXSHL below
      if (  f_2d(n)%unit == "ppb"  ) then  !  NOT ENABLED SO FAR !
         forall ( i=1:limax, j=1:ljmax )
           d_2d( n, i,j,IOU_INST) = &
               xn_shl(ind,i,j,KMAX_MID)  / (density(i,j)*to_molec_cm3)
         end forall
      else
         forall ( i=1:limax, j=1:ljmax )
           d_2d( n, i,j,IOU_INST) = xn_shl(ind,i,j,KMAX_MID)
         end forall
      end if

      if ( dbgP ) write(*,'(a,f8.2,3es12.3)') &
          'SHLSHLppb'//trim( species(ind)%name), thour, &
           xn_shl(ind,debug_li,debug_lj,KMAX_MID), &
            density(debug_li,debug_lj), to_molec_cm3

    case ( "PM25water" )      !water
      forall ( i=1:limax, j=1:ljmax ) &
        d_2d( n, i,j,IOU_INST) = PM25_water_rh50(i,j)
      ind2d_pmwater = n

    case ( "PMco_water" )      !water JUN21AERO
      forall ( i=1:limax, j=1:ljmax ) &
        d_2d( n, i,j,IOU_INST) = PMco_water_rh50(i,j)
      ind2d_pmcowater = n

    case ( "PM25" )      ! Need to add PMFINE + fraction NO3_c
      if(first_call)then
        call CheckStop(f_2d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind2d_pmfine <1,"Missing PMFINE output for "//trim(class))
        call CheckStop(iadv_NO3_C <1,"Unknown specie NO3_C")
     end if

      forall(i=1:limax,j=1:ljmax) &
        d_2d(n,i,j,IOU_INST) = d_2d(ind2d_pmfine,i,j,IOU_INST) + &
                               fracPM25 * &
            ( xn_adv(iadv_NO3_C,i,j,KMAX_MID) * ug_NO3_C &
            ) * cfac(iadv_NO3_C,i,j) * density(i,j)
        ind2d_pm25 = n

    case ( "PM_coarse" )      !
      if(first_call)then
        call CheckStop(f_2d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind2d_pm25 <1,"Missing PM25 output for "//trim(class))
        call CheckStop(ind2d_pm10 <1,"Missing PM10 output for "//trim(class))
     end if

      forall(i=1:limax,j=1:ljmax) &
        d_2d(n,i,j,IOU_INST) = d_2d(ind2d_pm10,i,j,IOU_INST) - &
                               d_2d(ind2d_pm25,i,j,IOU_INST) 

    case ( "SIA25" )   ! Need to subtract some NO3_c from SIA
      if(first_call)then
        call CheckStop(f_2d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind2d_sia <1,"Missing SIA output for "//trim(class))
        call CheckStop(iadv_NO3_C <1,"Unknown specie NO3_C")
      end if

      forall(i=1:limax,j=1:ljmax) & ! SUBTRACT, CAREFUL!
        d_2d(n,i,j,IOU_INST) = d_2d(ind2d_sia,i,j,IOU_INST) - &
                               (1-fracPM25)  * &
            ( xn_adv(iadv_NO3_C,i,j,KMAX_MID) * ug_NO3_C &
            ) * cfac(iadv_NO3_C,i,j) * density(i,j)
      call CheckStop( minval(d_2d(n,:,:,IOU_INST)) < 0.0 , dtxt//'ERROR NEG SIA!')
      if ( dbgP )  then
        write(*,*) "FRACTION SIA25 2d", n, ind2d_sia
        i= debug_li; j=debug_lj
        write(*,"(a,9es12.3)") "Adding SIA25 FRACTIONS:", &
          d_2d([ind2d_sia,n],i,j,IOU_INST), &
          (1-fracPM25) * xn_adv(iadv_NO3_C,i,j,KMAX_MID) * ug_NO3_C &
                   * cfac(iadv_NO3_C,i,j) * density(i,j)
      end if

    case ( "PM25_rh50" )      ! Need to add PMFINE + fraction NO3_c
      if(first_call)then
        call CheckStop(f_2d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind2d_pmfine <1,"Missing PMFINE output for "//trim(class))
        call CheckStop(ind2d_pmwater<1,"Missing PM25water output for "//trim(class))
        call CheckStop(iadv_NO3_C <1,"Unknown specie NO3_C")
      end if

      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = d_2d(ind2d_pmfine ,i,j,IOU_INST) &
                               + d_2d(ind2d_pmwater,i,j,IOU_INST) &
                               + fracPM25 * &
            ( xn_adv(iadv_NO3_C,i,j,KMAX_MID) * ug_NO3_C &
            ) * cfac(iadv_NO3_C,i,j) * density(i,j)
      end forall

      if( dbgP )  then
        write(*,*) "FRACTION PM25 2d", n, ind2d_pmfine, ind2d_pmwater
        i= debug_li; j=debug_lj
        write(*,"(a,4es12.3)") "Adding PM25 FRACTIONS:", &
          d_2d([ind2d_pmwater,ind2d_pmfine,n],i,j,IOU_INST), &
          fracPM25 * xn_adv(iadv_NO3_C,i,j,KMAX_MID) * ug_NO3_C &
                   * cfac(iadv_NO3_C,i,j) * density(i,j)
      end if

    case("PM25X")      ! Need to add PMFINE + fraction NO3_c
      if(first_call)then
        call CheckStop(f_2d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind2d_pmfine <1,"Missing PMFINE output for "//trim(class))
      end if
      if(any([iadv_NO3_C,iadv_EC_C_WOOD,iadv_EC_C_FFUEL,iadv_POM_C_FFUEL]<1))then
        if(first_call.and.MasterProc) write(*,*) &
          "WARNING: Derived - not all "//trim(class)//" species present. Skipping"
        cycle   !! Skip this case
      end if

      ! All this size class has the same cfac.
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = d_2d(ind2d_pmfine,i,j,IOU_INST) + &
                                 fracPM25 * &
            ( xn_adv(iadv_NO3_C      ,i,j,KMAX_MID) * ug_NO3_C       &
            + xn_adv(iadv_EC_C_WOOD  ,i,j,KMAX_MID) * ug_EC_C_WOOD   &
            + xn_adv(iadv_EC_C_FFUEL ,i,j,KMAX_MID) * ug_EC_C_FFUEL  &
            + xn_adv(iadv_POM_C_FFUEL,i,j,KMAX_MID) * ug_POM_C_FFUEL &
            ) * cfac(iadv_NO3_C,i,j) * density(i,j)
      end forall

    case("PM25X_rh50")      ! Need to add PMFINE + fraction NO3_c + water
      if(first_call)then
        call CheckStop(f_2d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind2d_pmfine <1,"Missing PMFINE output for "//trim(class))
        call CheckStop(ind2d_pmwater<1,"Missing PM25water output for "//trim(class))
      end if
      if(any([iadv_NO3_C,iadv_EC_C_WOOD,iadv_EC_C_FFUEL,iadv_POM_C_FFUEL]<1))then
        if(first_call.and.MasterProc) write(*,*) &
          "WARNING: Derived - not all "//trim(class)//" species present. Skipping"
        cycle   !! Skip this case
      end if

      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = d_2d(ind2d_pmfine ,i,j,IOU_INST) &
                               + d_2d(ind2d_pmwater,i,j,IOU_INST) &
                               + fracPM25 * &
            ( xn_adv(iadv_NO3_C      ,i,j,KMAX_MID) * ug_NO3_C       &
            + xn_adv(iadv_EC_C_WOOD  ,i,j,KMAX_MID) * ug_EC_C_WOOD   &
            + xn_adv(iadv_EC_C_FFUEL ,i,j,KMAX_MID) * ug_EC_C_FFUEL  &
            + xn_adv(iadv_POM_C_FFUEL,i,j,KMAX_MID) * ug_POM_C_FFUEL &
            ) * cfac(iadv_NO3_C,i,j) * density(i,j)
      end forall

    case("PM10_rh50")      ! Need to add PMFINE + fraction NO3_c
      if(first_call)then
        call CheckStop(f_2d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind2d_pm10   <1,"Missing PM10 output for "//trim(class))
        call CheckStop(ind2d_pmwater<1,"Missing PM25water output for "//trim(class))
      end if

      forall(i=1:limax,j=1:ljmax) &
        d_2d(n,i,j,IOU_INST) = d_2d(ind2d_pm10   ,i,j,IOU_INST) &
                             + d_2d(ind2d_pmwater,i,j,IOU_INST)
      if ( AERO%JUN21AERO ) then
        forall(i=1:limax,j=1:ljmax) &
         d_2d(n,i,j,IOU_INST) = d_2d(n,i,j,IOU_INST) + &
           d_2d(ind2d_pmcowater,i,j,IOU_INST)   !ST EQSAM
      end if

    case("AOD","AOD:GROUP","AOD:SPEC")  !/ Aerosol Optical Depth (new system)
      if(first_call)call AOD_init("Derived:"//trim(class))
      wlen=find_index(f_2d(n)%subclass,wavelength)! e.g. search "550nm" on array of wavelengths
      if(first_call)then
        call CheckStop(wlen<1,&
          "Unknown AOD wavelength "//trim(f_2d(n)%subclass))
        call CheckStop(.not.wanted_wlen(wlen),&
          "Unwanted AOD wavelength "//trim(f_2d(n)%subclass))
      end if

      ngrp = size(aod_grp)
      allocate(ingrp(ngrp))
      select case(class)
      case("AOD")
        ingrp(:)=.true.        ! take the full aod_grp
      case("AOD:GROUP")
        igrp = f_2d(n)%index
        do i=1,ngrp
          ingrp(i)=any(aod_grp(i)==chemgroups(igrp)%specs(:))
        end do
      case("AOD:SPEC")
        ispc = f_2d(n)%index
        ingrp(:)=(aod_grp(:)==ispc)
      end select
      forall ( i=1:limax, j=1:ljmax )&
        d_2d( n, i,j,IOU_INST) = SUM(AOD(:,i,j,wlen),MASK=ingrp)
      deallocate(ingrp)

    case ( "MAXADV" )
      if (  f_2d(n)%unit == "ppb"  ) then

         d_2d( n, 1:limax,1:ljmax,IOU_DAY) = &
           max( d_2d( n, 1:limax,1:ljmax,IOU_DAY),  &
                xn_adv(ind,1:limax,1:ljmax,KMAX_MID)  &
               * cfac(ind,1:limax,1:ljmax) )
         txt2 = "MAXADV ppb for " // trim( f_2d(n)%name)
       else
         d_2d( n, 1:limax,1:ljmax,IOU_DAY) = &
           max( d_2d( n, 1:limax,1:ljmax,IOU_DAY),  &
                xn_adv(ind,1:limax,1:ljmax,KMAX_MID)  &
               * cfac(ind,1:limax,1:ljmax) * density(1:limax,1:ljmax) )
         txt2 = "MAXADV ug for " // trim( f_2d(n)%name)
       end if

      if ( dbgP ) call write_debugadv(n,ind, &
                               density(debug_li,debug_lj), txt2 )

      !Monthly and yearly ARE averaged over days
      if(End_of_Day)then
        d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)
        nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
        if(    current_date%month >= 4 &
           .or.current_date%month <= 9 )then
        d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY)
        nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
        end if
      end if


    case ( "MAXSHL" )        ! Daily maxima - short-lived
      if (  f_2d(n)%unit /= "ppb"  ) then  ! Mix ratio so far
         forall ( i=1:limax, j=1:ljmax )
           d_2d( n, i,j,IOU_DAY) = max( d_2d( n, i,j,IOU_DAY), &
               xn_shl(ind,i,j,KMAX_MID) )
         end forall
      else
         forall ( i=1:limax, j=1:ljmax )
           d_2d( n, i,j,IOU_DAY) = max( d_2d( n, i,j,IOU_DAY), &
               xn_shl(ind,i,j,KMAX_MID)  / (density(i,j)*to_molec_cm3) )
         end forall
      end if


      if ( dbgP ) then
         write(*, *) "SHL:MAX.,to_molec_cm3 ", n, ind  , to_molec_cm3
         write(*,fmt="(a12,2i4,4es12.3)") "SHL MAX. ", n, ind  &
                , d_2d(n,debug_li,debug_lj,IOU_DAY) &
                ,  xn_shl(ind,debug_li,debug_lj,KMAX_MID)  &
                ,  density(debug_li,debug_lj), to_molec_cm3
      end if

      !Monthly and yearly ARE averaged over days
      if(End_of_Day)then
        d_2d(n,:,:,IOU_MON ) = d_2d(n,:,:,IOU_MON ) + d_2d(n,:,:,IOU_DAY)
        nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
        if(    current_date%month >= 4 &
           .or.current_date%month <= 9 )then
        d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY)
        nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
        end if
      end if

    case( "GRIDAOT" )!  Hardly used these days. The vegetation-specific
                     !  AOTs are handled in the Mosaic class and as
                     !  part of the dry dep calculations.
      if(first_call)&
        call CheckStop(iadv_o3<1,"Unknown specie O3")

      d_2d(n, 1:limax, 1:ljmax, IOU_INST) = Calc_GridAOTx( f_2d(n)%index )

      if( DEBUG%AOT .and. debug_proc ) then
        call datewrite("AOTDEBUG" // trim(f_2d(n)%name), n, &
         (/ zen(debug_li,debug_lj), real(f_2d(n)%index), &
            xn_adv(iadv_O3,debug_li,debug_lj,KMAX_MID)*&
               cfac(iadv_O3,debug_li,debug_lj)*PPBINV, &
            d_2d(n, debug_li, debug_lj, IOU_INST )  /) )
      end if

    case( "SOMO" )
      if(first_call)&
        call CheckStop(iadv_o3<1,"Unknown specie O3")

      !dt/7200: half a dt time step in hours
      !dayfrac "points" to the middle of the integration step
      dayfrac= (thour-(dt/7200.))/24. !must be < 1
      ntime=int(dayfrac*NTDAY )+1 !must be >=1 and <= NTDAY
      if(dayfrac<0)ntime=NTDAY !midnight

      !last value  (not averaged):
      D2_O3_DAY( : , : , ntime) =&
       xn_adv(iadv_o3,:,:,KMAX_MID)*cfac(iadv_o3,:,:)*PPBINV

      if(dayfrac<0)then !only at midnight: write on d_2d
        call somo_calc( n, f_2d(n)%index, dbgP )
        d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)
        if(f_2d(n)%avg) nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
        ! if(current_date%month>=4.and.current_date%month<=9)then
        d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY)
        if(f_2d(n)%avg) nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
        !NB overwritten anyway D2_O3_DAY = 0.
      end if

    case( "MaxD8M_26th", "MaxD8M_1st", "MaxD8M_2nd", "MaxD8M_3rd", "MaxD8M_4th")
      ! do nothing, it is taken care of by "MaxD8M" case
    case( "MaxD8M" , "AvgMDA8", "AvgMDA8AprSep") ! maximum daily eight-hour mean concentration
      if (first_call) then

        if (.not. allocated(D8M)) then ! we want to go through this only once
          iou = 0 !number of independent MaxD8M fields to output. for simplicity we add again the same if both MaxD8M and AvgMDA8 are chosen
          do i = 1, num_deriv2d
            if (f_2d(i)%class=="MaxD8M" .or. f_2d(i)%class=="AvgMDA8" .or. f_2d(i)%class=="AvgMDA8AprSep" ) then
              ! save index of species in f_2d
              f_2d(i)%index = find_index(f_2d(i)%txt,species_adv(:)%name, any_case=.true.)
              call CheckStop(f_2d(i)%index<0,"D8M: species "//trim(f_2d(i)%txt)//"not found")
              !force ug/m3
              call CheckStop(f_2d(i)%unit/="ug/m3",trim(f_2d(i)%class)//" must be in units of ug/m3 ")
              iou = iou + 1
              ! trick to save 2 integers in one:
              f_2d(i)%index = f_2d(i)%index  + iou * 100000
            end if
            if (f_2d(i)%class=="MaxD8M_26th") then
              call CheckStop(i_MaxD8M_26th > 0, "Only one MaxD8M_26th at a time implemented!")
              i_MaxD8M_26th = i
              make_MaxD8M_nth = .true.
            end if
            if (f_2d(i)%class=="MaxD8M_1st") then
              call CheckStop(i_MaxD8M_1st > 0, "Only one MaxD8M_nth at a time implemented!")
              i_MaxD8M_1st = i
              make_MaxD8M_nth = .true.
            end if
            if (f_2d(i)%class=="MaxD8M_2nd") then
              call CheckStop(i_MaxD8M_2nd > 0, "Only one MaxD8M_nth at a time implemented!")
              i_MaxD8M_2nd = i
              make_MaxD8M_nth = .true.
            end if
            if (f_2d(i)%class=="MaxD8M_3rd") then
              call CheckStop(i_MaxD8M_3rd > 0, "Only one MaxD8M_nth at a time implemented!")
              i_MaxD8M_3rd = i
              make_MaxD8M_nth = .true.
            end if
            if (f_2d(i)%class=="MaxD8M_4th") then
              call CheckStop(i_MaxD8M_4th > 0, "Only one MaxD8M_nth at a time implemented!")
              i_MaxD8M_4th = i
              make_MaxD8M_nth = .true.
            end if
          end do
          !we need one array to save the last 8 hourly concentrations, and one to make
          !the average over the last hour
          allocate(D8M(LIMAX,LJMAX,8,iou)) ! running last 8 hour values
          allocate(D8Max(LIMAX,LJMAX,iou)) ! max value of the 8 hour mean since 00:00
          allocate(hourM(LIMAX,LJMAX,iou)) ! hour Mean
          if (make_MaxD8M_nth) then
            !Note that for simplicity, we keep 26 largest , even in the case we only need
            !the first, second, 3rd or fourth
            !NB: for O3 only!
            allocate(D8_26Max(26,LIMAX,LJMAX)) ! 26 highest daily values up to current date 
            D8_26Max = 0.0 !init with low value
          end if
          hourM = 0.0
          D8Max = 0.0 !init with low value
          D8M = 0.0 !init with low value
        end if
      end if
      
      iou = f_2d(n)%index / 100000
      ii = mod(f_2d(n)%index, 100000)
      do j = 1,ljmax
        do i = 1,limax
          hourM(i,j,iou) = hourM(i,j,iou) + xn_adv(ii,i,j,KMAX_MID) * cfac(ii,i,j) * density(i,j) * 1.0e9 * species_adv(ii)%molwt/ATWAIR
        end do
      end do

      if (current_date%seconds == 0 .and. .not. first_call) then
        !one hour has past since last time here        
        !save last hour average
        ii = mod(current_date%hour,8) + 1 !note: this works only because 24 is a multiple of 8!
        timefrac = dt_advec/3600.0 ! inverse of number of timesteps in an hour
        D8M(:,:,ii,iou) =  hourM(:,:,iou) * timefrac
        hourM(:,:,iou) = 0.0
        tmpwork = 0.0
        ! update max value since 01:00
        do j = 1,ljmax
          do i = 1,limax
            af = 0.0
            do l = 1, 8
              af = af + D8M(i,j,l,iou) * 0.125 ! 8 hour average
            end do
            D8Max(i,j,iou) = max(D8Max(i,j,iou) , af) 
          end do
        end do
        !The last period of a day ends at 00:00
        !since the model output at the end of "EMEP day", which is different,
        !we must save only once per day at the right time
        if (current_date%hour == 0) then
          
          if (make_MaxD8M_nth .and. mod(f_2d(n)%index, 100000) == O3_ix-NSPEC_SHL) then           
            !For O3 we keep the 26 highest values and current value            
            do j = 1,ljmax
              do i = 1,limax
                !1) find smallest value stored so far
                i_26th = 1 
                do ii = 2, 26
                  if (D8_26Max(ii,i,j)<D8_26Max(i_26th,i,j)) i_26th=ii
                end do
                if (D8Max(i,j,iou)>D8_26Max(i_26th,i,j)) then
                  !2) Update i_26th, and D8_26Max(i_26th,i,j)
                  D8_26Max(i_26th,i,j) = D8Max(i,j,iou)
                  !3) find new smallest values among 26 largest
                  i_26th = 1
                  do ii = 2, 26
                    if (D8_26Max(ii,i,j)<D8_26Max(i_26th,i,j)) i_26th=ii
                  end do
                  !4) store the required largest value of all times
                  if (i_MaxD8M_26th > 0) then
                    d_2d(i_MaxD8M_26th,i,j,IOU_YEAR) = D8_26Max(i_26th,i,j)
                  endif
                end if
                !5) find new 4 largest values among 26 largest
                i1st=i_26th
                i2nd=i_26th
                i3rd=i_26th
                i4th=i_26th 
                do ii = 1, 26
                   if (D8_26Max(ii,i,j)>D8_26Max(i4th,i,j)) then
                      if(D8_26Max(ii,i,j) > D8_26Max(i1st,i,j)) then
                         i4th=i3rd
                         i3rd=i2nd
                         i2nd=i1st
                         i1st=ii
                      else if(D8_26Max(ii,i,j) > D8_26Max(i2nd,i,j)) then
                         i4th=i3rd
                         i3rd=i2nd
                         i2nd=ii
                      else if(D8_26Max(ii,i,j) > D8_26Max(i3rd,i,j)) then
                         i4th=i3rd
                         i3rd=ii
                      else 
                         i4th=ii
                      endif
                   endif
                end do
                
                if (i_MaxD8M_1st > 0) d_2d(i_MaxD8M_1st,i,j,IOU_YEAR) = D8_26Max(i1st,i,j)
                if (i_MaxD8M_2nd > 0) d_2d(i_MaxD8M_2nd,i,j,IOU_YEAR) = D8_26Max(i2nd,i,j)
                if (i_MaxD8M_3rd > 0) d_2d(i_MaxD8M_3rd,i,j,IOU_YEAR) = D8_26Max(i3rd,i,j)
                if (i_MaxD8M_4th > 0) d_2d(i_MaxD8M_4th,i,j,IOU_YEAR) = D8_26Max(i4th,i,j)
                
              end do
            end do
          end if
          
          if (class=="MaxD8M") then

            do j = 1,ljmax
              do i = 1,limax
                if(LENOUT2D>=IOU_DAY)d_2d(n,i,j,IOU_DAY) = D8Max(i,j,iou)
                !NB: max value over month, not average over daily max:
                if(LENOUT2D>=IOU_MON)d_2d(n,i,j,IOU_MON) = max(d_2d(n,i,j,IOU_MON) , D8Max(i,j,iou))
                !NB: max value over year, not average over daily max:
                d_2d(n,i,j,IOU_YEAR) = max(d_2d(n,i,j,IOU_YEAR), D8Max(i,j,iou))                                          
                D8Max(i,j,iou) = 0.0
              end do
            end do

         else if (class=="AvgMDA8") then
            !NB: at the end of the first day (day 2 hour 00:00), we actually start to write in the next month
            if (current_date%day == 2) count_AvgMDA8_m = 0
            count_AvgMDA8_m = count_AvgMDA8_m + 1
            count_AvgMDA8_y = count_AvgMDA8_y + 1
            w_m = 1.0/count_AvgMDA8_m
            w_y = 1.0/count_AvgMDA8_y
            do j = 1,ljmax
              do i = 1,limax
                if(LENOUT2D>=IOU_DAY)d_2d(n,i,j,IOU_DAY) = D8Max(i,j,iou)
                !average. makw weight according to already counted days
                if(LENOUT2D>=IOU_MON)d_2d(n,i,j,IOU_MON) = (1.0-w_m) * d_2d(n,i,j,IOU_MON) + w_m * D8Max(i,j,iou)
                d_2d(n,i,j,IOU_YEAR) = (1.0-w_y) * d_2d(n,i,j,IOU_YEAR) + w_y * D8Max(i,j,iou)
                D8Max(i,j,iou) = 0.0
              end do
            end do

          else if ( class=="AvgMDA8AprSep" ) then
            !NB: at the end of the first day (day 2 hour 00:00), we actually start to write in the next month
            if (current_date%day == 2) count_AvgMDA8AprSep_m = 0 
            if (current_date%day == 2) count_AvgMDA8AprSep_y = 0 
            count_AvgMDA8AprSep_m = count_AvgMDA8AprSep_m + 1!for monthes we output all anyway!
            if(current_date%month>=4 .and. current_date%month<=9)count_AvgMDA8AprSep_y = count_AvgMDA8AprSep_y + 1
            w_m = 1.0/count_AvgMDA8AprSep_m
            w_y = 1.0/count_AvgMDA8AprSep_y
            do j = 1,ljmax
              do i = 1,limax
                if(LENOUT2D>=IOU_DAY)d_2d(n,i,j,IOU_DAY) = D8Max(i,j,iou)
                !average. makw weight according to already counted days   
                if(LENOUT2D>=IOU_MON)d_2d(n,i,j,IOU_MON) = (1.0-w_m) * d_2d(n,i,j,IOU_MON) + w_m * D8Max(i,j,iou)
                if(current_date%month>=4 .and. current_date%month<=9) d_2d(n,i,j,IOU_YEAR) = (1.0-w_y) * d_2d(n,i,j,IOU_YEAR) + w_y * D8Max(i,j,iou)
                D8Max(i,j,iou) = 0.0
              end do
            end do
          end if
        end if
      end if

    case("PREC","WDEP","DDEP","VG","Rs","Rns","Gns","Mosaic","POD","SPOD","AOT")
    ! Nothing to do - all set in My_DryDep
      ! if(dbgP) write(*,"(2a,i4,a,es12.3)")"PROCESS ",trim(class),&
      !   n, trim(f_2d(n)%name), d_2d(n,debug_li,debug_lj,IOU_INST)

    case('FLYmax6h','FLYmax6h:SPEC')    ! Fly Level, 6 hourly maximum
      ! fl000-200: 0 to 20 kfeet, fl200-350: 20 to 35 kfeet, fl350-500: 35 to 50 kfeet
      read(subclass,"(a2,i3,a1,i3)") txt2, k, txt2, l
      fl0=k*30.48 ! [100 feet] to [m]
      fl1=l*30.48 ! [100 feet] to [m]
      call Units_Scale(f_2d(n)%unit,ind,af,needroa=needroa) ! only want needroa
      if(needroa)then
        tmpwork=maxval(xn_adv(ind,:,:,:)*roa(:,:,:,1),dim=3,&
                       mask=z_mid>=fl0.and.z_mid<=fl1)
      else
        tmpwork=maxval(xn_adv(ind,:,:,:),dim=3,&
                       mask=z_mid>=fl0.and.z_mid<=fl1)
      end if
      forall(i=1:limax,j=1:ljmax)&  ! use IOU_YEAR as a buffer
        d_2d(n,i,j,IOU_YEAR)=max(d_2d(n,i,j,IOU_YEAR),tmpwork(i,j))
    case('FLYmax6h:GROUP')           ! Fly Level, 6 hourly maximum
      ! fl000-200: 0 to 20 kfeet, fl200-350: 20 to 35 kfeet, fl350-500: 35 to 50 kfeet
      read(subclass,"(a2,i3,a1,i3)") txt2, k, txt2, l
      fl0=k*30.48 ! [100 feet] to [m]
      fl1=l*30.48 ! [100 feet] to [m]
      if(dbgP)print *,trim(subclass),fl0,fl1
      do k=1,KMAX_MID
        mask2d(:,:)=(z_mid(:,:,k)>=fl0.and.z_mid(:,:,k)<=fl1)
        if(.not.(any(mask2d)))cycle
        if(dbgP)print *,trim(subclass),k,count(mask2d)
        call group_calc(tmpwork(:,:),roa(:,:,k,1),f_2d(n)%unit,k,ind)
        forall(i=1:limax,j=1:ljmax,mask2d(i,j))&  ! use IOU_YEAR as a buffer
          d_2d(n,i,j,IOU_YEAR)=max(d_2d(n,i,j,IOU_YEAR),tmpwork(i,j))
      end do

    case ("COLUMN",'COLUMN:ADV',"COLUMN:SPEC") ! unit conversion factor stored in f_2d(n)%scale
      klow = KMAX_MID + 1 ! initialize too large
      if (f_2d(n)%subclass == "kmax") then
        klow = KMAX_MID
      else
        read(f_2d(n)%subclass,"(a1,i2)") txt2, klow ! Connvert e.g. k20 to klow=20
      end if
      call CheckStop(klow>KMAX_MID, "column definition too large: "// f_2d(n)%subclass)
      do j = 1, ljmax
        do i = 1, limax
          k = 1
          tmpwork(i,j) =  &
            xn_adv(ind,i,j,k)*roa(i,j,k,1)*(z_bnd(i,j,k)-z_bnd(i,j,k+1))
          do k = 2, klow   !!! KMAX_MID
            tmpwork(i,j) = tmpwork(i,j) + &
              xn_adv(ind,i,j,k)*roa(i,j,k,1)*(z_bnd(i,j,k)-z_bnd(i,j,k+1))

            if(DEBUG%COLUMN.and.dbgP.and.&
              i==debug_li.and.j==debug_lj) &
              write(*,"(a,3i4,a4,f8.3,f8.1,2es12.3)") &
                trim(f_2d(n)%name), n, ind, k, " => ", &
                  roa(i,j,k,1), z_bnd(i,j,k)-z_bnd(i,j,k+1), &
                  xn_adv(ind,i,j,k),tmpwork(i,j)
          end do ! k
          d_2d(n,i,j,IOU_INST) = tmpwork(i,j) ! unit conversion
                   ! is completed elsewere by *f_2d(n)%scale
        end do !i
      end do !j
      if(dbgP) write(*,"(a18,es12.3)") &
        "COLUMN:SPEC d2_2d",d_2d(n,debug_li,debug_lj,IOU_INST)*f_2d(n)%scale
      case ("COLUMN:SHL") ! unit conversion factor stored in f_2d(n)%scale
        klow = KMAX_MID + 1 ! initialize too large
        if (f_2d(n)%subclass == "kmax") then
          klow = KMAX_MID
        else
          read(f_2d(n)%subclass,"(a1,i2)") txt2, klow ! Connvert e.g. k20 to klow=20
        end if
        call CheckStop(klow>KMAX_MID, "column definition too large: "// f_2d(n)%subclass)
        do j = 1, ljmax
          do i = 1, limax
            k = 1
            tmpwork(i,j) =  &
              xn_shl(ind,i,j,k)*(z_bnd(i,j,k)-z_bnd(i,j,k+1))
            do k = 2, klow   !!! KMAX_MID
              tmpwork(i,j) = tmpwork(i,j) + &
                xn_shl(ind,i,j,k)*(z_bnd(i,j,k)-z_bnd(i,j,k+1))

              if(DEBUG%COLUMN.and.dbgP.and.&
                i==debug_li.and.j==debug_lj) &
                write(*,"(a,3i4,a4,f8.3,f8.1,2es12.3)") &
                  trim(f_2d(n)%name), n, ind, k, " => ", &
                    f_2d(n)%scale, z_bnd(i,j,k)-z_bnd(i,j,k+1), &
                    xn_shl(ind,i,j,k),tmpwork(i,j)
            end do ! k
            d_2d(n,i,j,IOU_INST) = tmpwork(i,j) ! unit conversion
                     ! is completed elsewere by *f_2d(n)%scale
          end do !i
        end do !j
        if(dbgP) write(*,"(a18,es12.3)") &
          "COLUMN:SHL d2_2d",d_2d(n,debug_li,debug_lj,IOU_INST)*f_2d(n)%scale
      case("COLUMN:GROUP")
      igrp = f_2d(n)%index
      call CheckStop(igrp<1,"NEG GRP "//trim(f_2d(n)%name))
      call CheckStop(igrp>size(chemgroups(:)%name), &
                            "Outside GRP "//trim(f_2d(n)%name))
      klow = KMAX_MID + 1 ! initialize too large
      if (f_2d(n)%subclass == "kmax") then
        klow = KMAX_MID
      else
        read(f_2d(n)%subclass,"(a1,i2)") txt2, klow ! Connvert e.g. k20 to klow=20
      end if
      call CheckStop(klow>KMAX_MID, "column definition too large: "// f_2d(n)%subclass)
      d_2d(n,:,:,IOU_INST) = 0.0
      do k=1,klow
        call group_calc(tmpwork(:,:),roa(:,:,k,1),f_2d(n)%unit,k,igrp)
        forall(i=1:limax,j=1:ljmax) &
          d_2d(n,i,j,IOU_INST) = d_2d(n,i,j,IOU_INST) &
            + tmpwork(i,j)*(z_bnd(i,j,k)-z_bnd(i,j,k+1)) ! unit conversion in group_calc
      end do
      if(dbgP) write(*,"(a18,es12.3)") &
        "COLUMN:GROUP d2_2d",d_2d(n,debug_li,debug_lj,IOU_INST)*f_2d(n)%scale

    case ( "EcoFrac" ) ! ODD TO HAVE FRAC AND AREA BELOW:"ECOAREA" )

      if( .not. first_call ) cycle ! Only need to do once
      if( f_2d(n)%Index == FULL_ECOGRID ) then
        km2_grid = (GRIDWIDTH_M*GRIDWIDTH_M) * 1.0e-6 ! km2
        forall ( i=1:limax, j=1:ljmax )
            d_2d(n,i,j,IOU_YEAR) =  EcoSystemFrac( f_2d(n)%Index ,i,j)&
                  * KM2_GRID /xm2(i,j)
        end forall
      else
        forall ( i=1:limax, j=1:ljmax )
            d_2d(n,i,j,IOU_YEAR) =  EcoSystemFrac( f_2d(n)%Index ,i,j)
        end forall
      end if
      if( dbgP ) &
        write(*,"(a18,a,i4,a,2es12.3)") "ECOD2D ", &
           " f2d:", f_2d(n)%Index, &
           " Frac", EcoSystemFrac( f_2d(n)%Index, debug_li,debug_lj), &
           !!" Index: ", DepEcoSystem(n)%Index, &
              d_2d( n, debug_li, debug_lj, IOU_YEAR)

    case ( "NatEmis" ) !emissions in kg/m2/s converted??

      forall ( i=1:limax, j=1:ljmax )
          d_2d(n,i,j,IOU_INST) =  EmisNat( f_2d(n)%Index,i,j )
      end forall
      !Not done, keep mg/m2  * GridArea_m2(i,j)
      if ( dbgP ) call write_debug(n,f_2d(n)%Index, "NatEmis")
      if( dbgP ) &
        call datewrite("NatEmis-in-Derived, kg/m2/s, "//trim(f_2d(n)%name), &
          f_2d(n)%Index, (/ EmisNat( f_2d(n)%Index, debug_li,debug_lj), &
                            maxval(EmisNat( f_2d(n)%Index, :,:)) /) )

    case ( "TotEmis" ) !emissions in kg/m2/s converted??

      forall ( i=1:limax, j=1:ljmax )
          d_2d(n,i,j,IOU_INST) =  EmisOut( i,j, f_2d(n)%Index)
      end forall
      !not done, to keep mg/m2 * GridArea_m2(i,j)
      if( dbgP .and. f_2d(n)%Index == 3  ) & ! CO:
        call datewrite("totEmis-in-Derived, still kg/m2/s", n, & !f_2d(n)%Index,&
              (/   EmisOut( debug_li,debug_lj, f_2d(n)%Index ) /) )

    case ( "SecEmis" ) !emissions in mg/m2 per sector

      iem=mod((f_2d(n)%Index-1),NEMIS_File)+1
      isec=(f_2d(n)%Index-1)/(NEMIS_File)
      forall ( i=1:limax, j=1:ljmax )
         d_2d(n,i,j,IOU_INST) =  SecEmisOut( i,j, iem, isec2SecOutWanted(isec))
      end forall

    case ( "Emis_mgm2_DMS" )      ! DMS
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = O_DMS%map(i,j)
      end forall

    case ( "Emis_mgm2_Ocean_NH3" )      ! Ocean_NH3
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = O_NH3%map(i,j)
      end forall

    case ( "EmisSplit_mgm2" )      ! Splitted total emissions (Inclusive natural)
      forall ( i=1:limax, j=1:ljmax )
        d_2d( n, i,j,IOU_INST) = SplitEmisOut(i,j,f_2d(n)%Index)
      end forall

    case ( "EXT" )
    ! Externally set for IOU_INST (in other routines); so no new work
    ! needed except decision to accumalate to yearly or not.
      if(dbgP) write(*,"(a18,i4,a12,a4,es12.3)")"EXT d_2d",&
        n, f_2d(n)%name, " is ", d_2d(n,debug_li,debug_lj,IOU_INST)

    case ( "SURF_MASS_GROUP","SURF_PPB_GROUP" ) !
      igrp = f_2d(n)%index
      call CheckStop(igrp<1,"NEG GRP "//trim(f_2d(n)%name))
      call CheckStop(igrp>size(chemgroups(:)%name), &
                            "Outside GRP "//trim(f_2d(n)%name))
      ngrp = size(chemgroups(igrp)%specs)
      if(first_call.and.f_2d(n)%unit=="ug/m3")then
        select case(chemgroups(igrp)%name)
          case("PMFINE"); ind2d_pmfine = n
          case("SIA");    ind2d_sia = n
          case("PM10");   ind2d_pm10 = n
        end select
        if(dbgP) write(*,"(a,3i4,2a15)") &
          "Deriv: CASEGRP - 2DFOUND "//trim(class), n, igrp, ngrp,&
          trim(chemgroups(igrp)%name), trim(f_2d(n)%name)
      end if

      if(dbg0) then
        write(*,"(a,3i5,3(1x,a))")dtxt//"CASEGRP:"//trim(f_2d(n)%name), &
          n, igrp, ngrp, trim(class), trim(subclass), &
            trim(chemgroups(igrp)%name)
        write(*,"(a,88i4)") "CASEGRP:", chemgroups(igrp)%specs
        write(*,*) "CASEGRPunit ", trim(f_2d(n)%unit)
      end if
     ! GROUP:§
      call group_calc(d_2d(n,:,:,IOU_INST),density,f_2d(n)%unit,0,igrp,&
                      semivol=(f_2d(n)%subclass=='FSOA'))
     !....
      if( dbgP .and. igrp== igrp_OM25) then
         i=debug_li
         j=debug_lj
         write(*,'(2a,4es14.5,a)') dtxt//'IOM_GRP ', f_2d(n)%name, &
             d_2d(n,i,j,IOU_INST), density(i,j),f_2d(n)%scale
      end if

      if(subclass=='LocFrac_corrected')then
        !TO BE DEVELOPED
!         do ii= 1,ngrp
!            do ipoll=1,uEMEP%Npoll        
!               do i=1,uEMEP%poll(ipoll)%Nix
!                  if(chemgroups(igrp)%specs(ii)==uEMEP%poll(ipoll)%ix(i)+NSPEC_SHL)goto 54
!               enddo
!            enddo
!         enddo
!         if(me==0)write(*,*)'WARNING, no local fractions found for ',trim(class),' group ',chemgroups(igrp)%name
!         54 continue
!         if(me==0.and. first_call)then
!            write(*,*)'local fractions found for ',trim(class),&
!              ' group index ',igrp,' name ',trim(chemgroups(igrp)%name),&
!              ' locfrac pollutant ',trim(uEMEP%poll(ipoll)%emis)
!            do iisec=1,uEMEP%poll(ipoll)%Nsectors
!               isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!               isec=uEMEP%poll(ipoll)%sector(iisec)
!               write(*,*)'local correction for sector',isec,', pollutant ',trim(uEMEP%poll(ipoll)%emis)
!            enddo
!         endif
!
!         do j=1,ljmax
!         do i=1,limax
!         default_frac=0.0!local, but any sector that is not explicit
!         tot_frac=0.0!all local (any sector)
!         loc_frac_corr=0.0
!         !isec is sector (number between 1 and 11)
!         !iisec is index over available sectors
!         !isec_poll is an internal uEMEP index that is a combination of sector and pollutant indices 
!         do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!            if(isec/=0)then
!               default_frac = default_frac - loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!               tot_frac =  tot_frac + loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!            endif
!            if(isec==0)default_frac = default_frac + loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!         enddo
!         default_frac=max(0.0,default_frac)!in case "sec=0" not available
!         tot_frac = tot_frac + default_frac
!         if(i==5.and.j==5)then
!         do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!            write(*,*)me,' sector ',isec,' fraction ',loc_frac(isec_poll,0,0,i,j,KMAX_MID)
!         enddo
!            write(*,*)me,' remaining sectors fraction ', default_frac,' tot ',tot_frac
!         endif
!
!         do iisec=1,uEMEP%poll(ipoll)%Nsectors
!            isec_poll=uEMEP%poll(ipoll)%sec_poll_ishift+iisec
!            isec=uEMEP%poll(ipoll)%sector(iisec)
!            if(isec/=0)then
!               loc_frac_corr=loc_frac_corr+loc_frac(isec_poll,0,0,i,j,KMAX_MID)*2!*LocEmisFac(isec)
!            endif
!         enddo
!         loc_frac_corr=loc_frac_corr+default_frac*2!*LocEmisFac_default
!         loc_frac_corr=loc_frac_corr+(1-tot_frac)!weight one for pollutants from other sources
!               d_2d( n, i,j,IOU_INST) = d_2d( n, i,j,IOU_INST) &
!                     * loc_frac_corr
!         enddo
!         enddo
      endif
      if( dbgP )then
        i= debug_li; j=debug_lj
        if(n==ind2d_pmfine ) &
          write(*,"(a,i4,es12.3)") "PMFINE FRACTION:"   ,n,d_2d(n,i,j,IOU_INST)
        if(n==ind2d_pm10   ) &
          write(*,"(a,i4,es12.3)") "PM10 FRACTION:"     ,n,d_2d(n,i,j,IOU_INST)
        write(*,*) "CASErho     ", density(i,j)
      end if

    case("USET")
      if(dbgP) write(*,"(a18,i4,a12,a4,es12.3)")"USET d_2d",&
        n, f_2d(n)%name, " is ", d_2d(n,debug_li,debug_lj,IOU_INST)

    case  default

      if ( dbgP ) then
         if( i == debug_li .and. j == debug_lj ) &
           write(*,"(a,i3,4a)") "My_Deriv Defaults called n=",&
              n, " Type ",trim(class), " Name ", trim( f_2d(n)%name )

           write(*,"(a,i3,i8,i4,a)") &
              "My_Deriv index?, nav? length?, class? ", ind,&
              nav_2d(n,IOU_INST), len(f_2d%class), trim(f_2d(n)%class)
           write(*,*) "My_Deriv index?, avg ", f_2d(n)%avg
       end if

       call My_DerivFunc( d_2d(n,:,:,IOU_INST), class ) ! , density )

    end select

    !*** add to daily, monthly and yearly average, and increment counters
    select case(f_2d(n)%class)
    case('FLYmax6h','FLYmax6h:SPEC','FLYmax6h:GROUP')
    ! Fly Level, 6 hourly maximum:  only need IOU_HOUR,IOU_HOUR_INST
      d_2d(n,:,:,IOU_INST) = d_2d(n,:,:,IOU_YEAR) ! use IOU_YEAR as a buffer
      d_2d(n,:,:,IOU_HOUR) = d_2d(n,:,:,IOU_YEAR)
      if(mod(current_date%hour,6)==0)&  ! reset buffer
        d_2d(n,:,:,IOU_YEAR)=0.0
    case("MAXADV","MAXSHL","SOMO")
    !  MAXADV and MAXSHL and SOMO needn't be summed here.
    !  These d_2d ( MAXADV, MAXSHL, SOMO) are set elsewhere
    case default
      af = 1.0 ! accumulation factor
      if(f_2d(n)%dt_scale) af=dt_advec !need to scale with dt_advec

      ! only accumulate outputs if they are wanted (will be written out)
      do iou=1,LENOUT2D
        if(iou==IOU_INST)cycle
        if(.not.wanted_iou(iou,f_2d(n)%iotype,ONLY_IOU))cycle
        d_2d(n,:,:,iou) = d_2d(n,:,:,iou) + d_2d(n,:,:,IOU_INST)*af
        if(f_2d(n)%avg) nav_2d(n,iou) = nav_2d(n,iou) + 1
      end do
    end select

  end do   ! num_deriv2d

  !****** 3-D fields **************************

  if(dbgP)& ! RUN through indices etc.
    write(*, "(a12,2i4,f12.3)") "3D3D TIME ",  me, num_deriv3d, &
            (current_date%hour+current_date%seconds/3600.0)


  do n = 1, num_deriv3d

    ind   = f_3d(n)%index
    class = f_3d(n)%class
    name  = f_3d(n)%name

    if(f_3d(n)%unit=="ppb") then
      inv_air_density3D(:,:,:) = 1.0
    else
      forall( i=1:limax, j=1:ljmax, k=1:KMAX_MID )&
        inv_air_density3D(i,j,k) = 1.0/( roa(i,j,k,1) * to_molec_cm3 )
    end if

    select case (class)

    case ( "MET3D" )

      imet_tmp = find_index(f_3d(n)%subclass, met(:)%name ) ! subclass has meteo name from MetFields
      if(imet_tmp>0) then
        if(met(imet_tmp)%dim==3)then
          if( MasterProc.and.first_call) write(*,*) "MET3D "//trim(f_3d(n)%name), &
                imet_tmp, met(imet_tmp)%field(2,2,KMAX_MID,1)
          forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
            d_3d(n,i,j,k,IOU_INST)=met(imet_tmp)%field(i,j,lev3d(k),1)
        elseif(MasterProc.and.first_call)then
          write(*,*) "Warning: requested 2D field with MET3D: ",trim(f_3d(n)%name)
        end if
      else
        !make derived fields:
        select case ( f_3d(n)%subclass )
        case ("inv_wind_speed_3D")
           !take 0.2m/s as minimum wind speed. Otherwise the average will be infinite if the wind gets zero any time.
           forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
                d_3d(n,i,j,k,IOU_INST)=1.0/(0.2+sqrt(u_mid(i,j,lev3d(k))**2+v_mid(i,j,lev3d(k))**2))
        case("wind_speed_3D")
           forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
                d_3d(n,i,j,k,IOU_INST)=sqrt(u_mid(i,j,lev3d(k))**2+v_mid(i,j,lev3d(k))**2)
        case("RH_3D") !relative humidity
           do k=1, num_lev3d
           do j=1, ljmax
           do i=1, limax
              pp = A_mid(k) + B_mid(k)*ps(i,i,1)
              itemp= nint(th(i,i,k,1) * Tpot_2_T(pp))
              qsat = 0.622 * tab_esat_Pa( itemp ) / pp
              d_3d(n,i,j,k,IOU_INST)=min(q(i,i,k,1)/qsat,1.0)       
           end do
           end do
           end do
        case default
           if(MasterProc) write(*,*) "MET3D NOT FOUND"//trim(f_3d(n)%name)//":"//trim(f_3d(n)%subclass)
           d_3d(n,:,:,:,IOU_INST)=0.0
        end select
        if( MasterProc.and.first_call)write(*,*) "MET3D "//trim(f_3d(n)%name)//' '//f_3d(n)%subclass, d_3d(n,2,2,num_lev3d,IOU_INST)
      end if

    ! Simple advected species:
    case ( "ADV" )
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=xn_adv(ind,i,j,lev3d(k))

    case ( "BGN" )
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=xn_bgn(ind,i,j,lev3d(k))

    case ( "PM25water" )         !particle water
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=PM25_water(i,j,lev3d(k))
      ind3d_pmwater = n

    case ( "PM25" )      ! Need to add PMFINE + fraction NO3_c
      if(first_call)then
        call CheckStop(f_3d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind3d_pmfine <1,"Missing PMFINE output for "//trim(class))
        call CheckStop(iadv_NO3_C <1,"Unknown specie NO3_C")
      end if

      forall (i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST) = d_3d(ind3d_pmfine,i,j,k,IOU_INST) + &
                                 fracPM25 * &
            ( xn_adv(iadv_NO3_C,i,j,lev3d(k)) * ug_NO3_C &
            ) * roa(i,j,lev3d(k),1)

    case ( "PM25_wet" )         ! Need to add PMFINE + fraction NO3_c
      if(first_call)then
        call CheckStop(f_3d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind3d_pmfine <1,"Missing PMFINE output for "//trim(class))
        call CheckStop(ind3d_pmwater<1,"Missing PM25water output for "//trim(class))
        call CheckStop(iadv_NO3_C <1,"Unknown specie NO3_C")
      end if

      forall (i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST) = d_3d(ind3d_pmfine ,i,j,k,IOU_INST) &
                               + d_3d(ind3d_pmwater,i,j,k,IOU_INST) &
                               + fracPM25 * &
            ( xn_adv(iadv_NO3_C,i,j,lev3d(k)) * ug_NO3_C &
            ) * roa(i,j,lev3d(k),1)

      if( dbgP )  then
        write(*,*) "FRACTION PM25 3d", n, ind3d_pmfine, ind3d_pmwater
        i= debug_li; j=debug_lj; k=1; l=lev3d(k)
        write(*,"(a,4es12.3)") "Adding PM25FRACTIONS:", &
          d_3d([ind3d_pmwater,ind3d_pmfine,n],i,j,k,IOU_INST), &
          ug_NO3_C * xn_adv(iadv_NO3_C,i,j,l) * roa(i,j,l,1)
      end if

    case ( "SIA25" )      ! Need to subtract some NO3_c from SIA
      if(first_call)then
        call CheckStop(f_3d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind3d_sia <1,"Missing 3D-SIA output for "//trim(class))
        call CheckStop(iadv_NO3_C <1,"Unknown specie NO3_C")
      end if

      forall (i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST) = d_3d(ind3d_sia,i,j,k,IOU_INST) - &
                                 (1-fracPM25) * &
            ( xn_adv(iadv_NO3_C,i,j,lev3d(k)) * ug_NO3_C &
            ) * roa(i,j,lev3d(k),1)

    case("PM10_wet")      ! Need to add PMFINE + fraction NO3_c
      if(first_call)then
        call CheckStop(f_3d(n)%unit,"ug/m3","Wrong unit for "//trim(class))
        call CheckStop(ind3d_pm10   <1,"Missing PM10 output for "//trim(class))
        call CheckStop(ind3d_pmwater<1,"Missing PM25water output for "//trim(class))
      end if

      forall (i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST) = d_3d(ind3d_pm10   ,i,j,k,IOU_INST) &
                               + d_3d(ind3d_pmwater,i,j,k,IOU_INST)

    case ("XKSIG00", "XKSIG12" ) !hf hmix Kz_m2s
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=Kz_m2s(i,j,lev3d(k))

    case ("TH" ) ! Pot. temp (needed for cross sections)
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=th(i,j,lev3d(k),1)

    case ("T" ) ! Absolute Temperature
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=th(i,j,lev3d(k),1)&
            *exp(KAPPA*log((A_mid(lev3d(k))+ B_mid(lev3d(k))*ps(i,j,1))*1.e-5))

    case ("pressure_3D" ) !
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST) = &
                     (A_mid(lev3d(k)) + B_mid(lev3d(k))*ps(i,j,1))*1.e-2 ! hPa

    case ( "MAX3DSHL" ) ! Daily maxima - short-lived
      call CheckStop(f_3d(n)%unit=="ppb","Asked for MAX3DSHL ppb ")
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=max(d_3d(n,i,j,k,IOU_INST),&
                             xn_shl(ind,i,j,lev3d(k)) &
                       *inv_air_density3D(i,j,lev3d(k)))

      if(dbgP) write(*,"(a13,i4,f8.3,3es12.3)") "3D3D MAX3DSHL", n, thour, &
        xn_shl(ind,debug_li,debug_lj,KMAX_MID), &
        1.0/inv_air_density3D(debug_li,debug_lj,KMAX_MID), &
        d_3d(n,debug_li,debug_lj,num_lev3d,IOU_INST)

    case ( "MAX3DADV" )
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=max(d_3d(n,i,j,k,IOU_INST),&
                             xn_adv(ind,i,j,lev3d(k)))

      if(dbgP) write(*,"(a12,i4,f8.3,4es12.3)") "SET MAX3DADV", n, thour, &
        xn_adv(ind,debug_li,debug_lj,KMAX_MID), &
        d_3d(n,debug_li,debug_lj,num_lev3d,IOU_INST)

    case ( "SHL" )
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=xn_shl(ind,i,j,lev3d(k))&
                         *inv_air_density3D(i,j,lev3d(k))

    case ( "3D_PPB_SPEC" )
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=xn_adv(ind,i,j,lev3d(k))

      if(dbgP) call write_debugadv(n,ind, 1.0, "3D PPB OUTS",IS3D=.true.)

    case ( "3D_PPB_SHL" )
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=xn_shl(ind,i,j,lev3d(k))

    case ( "3D_MASS_SPEC" )
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=xn_adv(ind,i,j,lev3d(k))*roa(i,j,lev3d(k),1)

      if(dbgP) call write_debugadv(n,ind, 1.0, "3D UG OUTS",IS3D=.true.)

    case ( "3D_MASS_GROUP","3D_PPB_GROUP" ) !
      igrp = f_3d(n)%index
      call CheckStop(igrp<1,"NEG GRP "//trim(f_3d(n)%name))
      call CheckStop(igrp>size(chemgroups(:)%name), &
                            "Outside GRP "//trim(f_3d(n)%name))
      ngrp = size(chemgroups(igrp)%specs)
      if(first_call.and.f_3d(n)%unit=="ug/m3")then
        select case (chemgroups(igrp)%name)
          case("PMFINE"); ind3d_pmfine = n
          case("SIA");    ind3d_sia = n
          case("PM10");   ind3d_pm10 = n
        end select
        if(MasterProc) write(*,"(a,3i4,2a15)") &
          "Deriv: CASEGRP - 3DFOUND "//trim(class), n, igrp, ngrp,&
          trim(chemgroups(igrp)%name), trim(f_3d(n)%name)
      end if

      if(dbg0) then
        write(*,*) "3DCASEGRP ", n, igrp, ngrp, trim(class)
        write(*,*) "3DCASENAM ", trim(f_3d(n)%name)
        write(*,*) "3DCASEGRP:", chemgroups(igrp)%specs
        write(*,*) "3DCASEunit", trim(f_3d(n)%unit)
      end if
      do k=1,num_lev3d
        call group_calc(d_3d(n,:,:,k,IOU_INST),roa(:,:,lev3d(k),1),&
                        f_3d(n)%unit,lev3d(k),igrp)
      end do

      if( dbgP )then
        i= debug_li; j=debug_lj; k=1; l=lev3d(k)
        if(n==ind3d_pmfine ) &
          write(*,"(a,i4,es12.3)") "PMFINE 3d FRACTION:",n,d_3d(n,i,j,k,IOU_INST)
        if(n==ind3d_pm10   ) &
          write(*,"(a,i4,es12.3)") "PM10 3d FRACTION:"  ,n,d_3d(n,i,j,k,IOU_INST)
        write(*,*) "CASErho     ", roa(i,j,l,1)
      end if

    case ( "Kz" )
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=Kz_m2s(i,j,lev3d(k))

    case ("Z_MID","Z")    ! Mid-layer heigh
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=z_mid(i,j,lev3d(k))

    case ("Z_BND","Zlev") ! Mid-layer heigh
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=z_bnd(i,j,lev3d(k))

    case("dZ_BND","dZ")   ! level thickness
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=z_bnd(i,j,lev3d(k))-z_bnd(i,j,lev3d(k)+1)

    case("EXT","EXT:GROUP","EXT:SPEC")  !/ Extinction coefficient (new system)
      if(first_call)call AOD_init("Derived:"//trim(class))
      wlen=find_index(f_3d(n)%subclass,wavelength)! e.g. search "550nm" on array of wavelengths
      if(first_call)then
        call CheckStop(wlen<1,&
          "Unknown EXT wavelength "//trim(f_3d(n)%subclass))
        call CheckStop(.not.(wanted_wlen(wlen).and.wanted_ext3d),&
          "Unwanted EXT wavelength "//trim(f_3d(n)%subclass))
      end if

      ngrp = size(aod_grp)
      allocate(ingrp(ngrp))
      select case(class)
      case("EXT")       ! take the full aod_grp
        ingrp(:) = .true.
      case("EXT:GROUP") ! cherry pick from aod_grp
        igrp = f_3d(n)%index
        do i=1,ngrp
          ingrp(i)=any(aod_grp(i)==chemgroups(igrp)%specs(:))
        end do
      case("EXT:SPEC")
        ispc = f_3d(n)%index
        ingrp(:)=(aod_grp(:)==ispc)
      end select
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=SUM(Extin_coeff(:,i,j,lev3d(k),wlen),MASK=ingrp)
      deallocate(ingrp)

    case("EmFFprof")
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
           d_3d(n,i,j,k,IOU_INST) = Emis_CO_Profile(i,j,lev3d(k))

   case("USET")
      if(dbgP) write(*,"(a18,i4,a12,a4,es12.3)")"USET d_3d",&
        n, f_3d(n)%name, " is ", d_3d(n,debug_li,debug_lj,num_lev3d,IOU_INST)

    case default
      write(*,"(a,2i3,3a)") "*** NOT FOUND",n,ind, trim(f_3d(n)%name),&
               ";Class:", trim(f_3d(n)%class)
      write(unit=errmsg,fmt=*) "Derived 3D class NOT FOUND", n, ind, &
                       trim(f_3d(n)%name)," ",trim(f_3d(n)%class)," ",trim(class)
      call CheckStop( errmsg )
    end select

    !*** add to monthly and yearly average, and increment counters
    select case(f_3d(n)%class)
    case("MAX3DSHL","MAX3DADV")
    ! For the MAX3D possibilities, we store maximum value of the
    !   current day in the IOU_INST variables.
    !   These are then added into IOU_MON **only** at the end of each day.
    ! (NB there is an error made on 1st day used, since only 1st 6 hours
    !  are checked. Still, not much happens on 1st Jan.... ;-)
      if(End_of_Day)then

        ! only accumulate outputs if they are wanted (will be written out)
        do iou=1,LENOUT3D
          if(iou==IOU_INST)cycle
          if(.not.wanted_iou(iou,f_3d(n)%iotype,ONLY_IOU))cycle
          d_3d(n,:,:,:,iou) = d_3d(n,:,:,:,iou) + d_3d(n,:,:,:,IOU_INST)
          if(f_3d(n)%avg) nav_3d(n,iou) = nav_3d(n,iou) + 1
        end do

        if( dbgP ) then
          write(*,fmt="(a20,a9,i4,f8.3,2es12.3)") "END_OF_DAY MAX3D", &
            f_3d(n)%class, n, thour,  &
            d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_MON ),&
            d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST )
          write(*,"(a20,i4,2x,6i6)") "END_OF_DAY NAV ", &
            n, (nav_3d(n,i), i=1,LENOUT3D)
        end if

        d_3d(n,:,:,:,IOU_INST ) = 0.0  !! Reset d_3d

      end if ! End_of_Day
    case default

      af = 1.0 ! accumulation factor
      if(f_3d(n)%dt_scale) af=dt_advec !need to scale with dt_advec

      ! only accumulate outputs if they are wanted (will be written out)
      do iou=1,LENOUT3D
        if(iou==IOU_INST)cycle
        if(.not.wanted_iou(iou,f_3d(n)%iotype,ONLY_IOU))cycle
        d_3d(n,:,:,:,iou) = d_3d(n,:,:,:,iou) + d_3d(n,:,:,:,IOU_INST)*af
        if(f_3d(n)%avg) nav_3d(n,iou) = nav_3d(n,iou) + 1
      end do

    end select
  end do

  !the uemep fields do not fit in the general d_3d arrays. Use ad hoc routine
  if(USES%LocalFractions .and. .not. present(ONLY_IOU))then
    call lf_av(dt,End_of_Day)
  endif

  first_call = .false.
end subroutine Derived
!=========================================================================
subroutine DerivedProds(text,dt)
!*** DESCRIPTION
!  Calculates chemical changes by comparing values before and  after
!  chemistry subroutine. Intended to be a more flexible version of the old
!  PRODO3  calculation

  character(len=*), intent(in) :: text  ! "Before" or "After"
  real,             intent(in) :: dt    ! timestep (s)
  real :: timefrac                      ! dt as fraction of hour (3600/dt)
  integer :: ind



! if ( num_deriv3d < 1 ) print *, "DerivedProds "//text, num_deriv3d
  if ( num_deriv3d < 1 ) return
  if (.not. any( f_3d%class == "PROD" ) ) return

  timefrac = dt/3600.0

!****** 3-D fields **************************
  do n = 1, num_deriv3d
    if(f_3d(n)%class/="PROD")cycle
    ind = f_3d(n)%index
    select case ( text )
    case("Before")    !! Initialise to xn_adv
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=xn_adv(ind,i,j,lev3d(k))
    case("After")     !! Calculate change
      forall(i=1:limax,j=1:ljmax,k=1:num_lev3d) &
        d_3d(n,i,j,k,IOU_INST)=d_3d(n,i,j,k,IOU_INST)&
                              -xn_adv(ind,i,j,lev3d(k))
    end select
  end do
end subroutine DerivedProds
!=========================================================================
subroutine ResetDerived(period)
  integer, intent(in) :: period   ! Either IOU_DAY or IOU_MON

  if(num_deriv2d>0 .and. period<=LENOUT2D) then
    nav_2d  (:,period) = 0
    d_2d(:,:,:,period) = 0.0
  end if

  if(num_deriv3d>0 .and. period<=LENOUT3D) then
    nav_3d    (:,period) = 0
    d_3d(:,:,:,:,period) = 0.0
  end if

end subroutine ResetDerived
!=========================================================================
subroutine group_calc( g2d, density, unit, ik, igrp,semivol)

  !/--  calulates e.g. SIA = SO4 + pNO3_f + pNO3_c + aNH4
  ! (only SIA converted to new group system so far, rv3_5_6 )
  !/--  calulates also PM10  = SIA + PPM2.5 + PPMCOARSE

  real, dimension(:,:), intent(out) :: g2d  ! i,j section of d_2d arrays
  real, intent(in), dimension(LIMAX,LJMAX)  :: density
  character(len=*), intent(in) :: unit
  integer, intent(in) :: ik,igrp
  logical, intent(in), optional :: semivol

  integer, pointer, dimension(:) :: gspec=>null()       ! group array of indexes
  real,    pointer, dimension(:) :: gunit_conv=>null()  ! & unit conv. factors
  logical :: needroa
  integer :: kk, iadv, itot, nspec
 ! needed for OM25_p, OM25 groups
  integer, save:: iadv_PMf= -999,iadv_OM25p= -999,iadv_BGND= -999,iadvDep=-999
  real :: fac, Fgas                       ! FSOA
  logical ::  semivol_wanted, dbgPt       ! FSOA
  logical ::  first_call    = .true.  ! FSOA
  logical ::  first_semivol_call = .true., ParticlePhaseOutputs   ! FSOA
  character(len=*), parameter :: dtxt = 'DrvGrpCalc:'

!OM25: To solve complications with OM, we need:
  if ( first_call ) then
     iadv_PMf = find_index('SO4',species(:)%name ) - NSPEC_SHL
     iadv_OM25p = find_index('OM25_p',species(:)%name, any_case=.true. ) &
               - NSPEC_SHL
     iadv_BGND  = find_index('OM25_BGND',species(:)%name, any_case=.true. ) &
               - NSPEC_SHL
  end if !OM25 end

  kk = ik
  if(ik==0) kk = KMAX_MID

  semivol_wanted=.false.
  if(present(semivol)) semivol_wanted = semivol
  ! NB global use of n is dangerous!
  ParticlePhaseOutputs = ( (index(f_2d(n)%name, 'ug_PM' )>0) .or. &
                           (index(f_2d(n)%name, 'ugC_PM')>0) )

  if( dbgP .and. first_call ) &
    write(*,"(a,L1,3i4,2a16,L2,i4)") dtxt//"SGROUP:",debug_proc,me,ik, kk, &
      trim(chemgroups(igrp)%name), trim(unit), semivol_wanted, iadv_OM25p

  call Group_Units(igrp,unit,gspec,gunit_conv, debug=dbgP,needroa=needroa)

  if(semivol_wanted.and.debug_proc .and. first_semivol_call ) then
    write(*,"(a,L2,2i4,a16,2i4)") dtxt//"GROUP-FSOA",debug_proc,me,ik, &
      trim(chemgroups(igrp)%name), size(gspec), size(gunit_conv)
    !write(*,*) dtxt//"GROUP-FSOA-GSPEC", gspec
    !write(*,*) dtxt//"GROUP-FSOA-GUNIT", gunit_conv
  end if

  do j=1,ljmax
    do i = 1, limax
      g2d(i,j) = 0.0
      dbgPt = ( dbgP .and. i==debug_li .and. j==debug_lj )
      do nspec = 1, size(gspec)
        iadv  = gspec(nspec)
        iadvDep = iadv ! used for cfac, needed since OM25_p is not in CM_DryDep
        if ( iadv == iadv_OM25p ) iadvDep= iadv_PMf 
        if ( iadv == iadv_BGND  ) iadvDep= iadv_PMf 
        itot  = iadv + NSPEC_SHL
        fac = 1.0

        ! With SOA modelling some compounds are semivolatile and others
        ! non-volatile. If in a group XXX which asks for ugPM the latter's
        ! mass is correct. If semivolatile, we need to calculate the PM
        ! fraction and just add this.

        !if(first_call.and.debug_proc) write(*,"(a,3i4,2(1x,a),2i4)") &
        !  dtxt//"FSOA check ", nspec, itot, igrp, trim(chemgroups(igrp)%name),&
        !  trim(species(itot)%name), FIRST_SEMIVOL, LAST_SEMIVOL

        !if(all([semivol_wanted,itot>=FIRST_SEMIVOL,itot<=LAST_SEMIVOL])) then
        if ( dbgPt .and. first_call  ) write(*,'(a,3i4)')dtxt//'IOM_choice '//trim((f_2d(n)%name))//&
             ':'//trim(species(itot)%name), index(f_2d(n)%name, 'ug_PM' ), nspec, size(gspec)

        if(itot>=FIRST_SEMIVOL .and. itot<=LAST_SEMIVOL) then

          ! Fgas3d only defined for KCHEMTOP down
           Fgas = Fgas3d(itot,i,j, max(kk,KCHEMTOP) )

           if ( ParticlePhaseOutputs ) then ! particle phase wanted
             fac = 1 - Fgas
             !iadvDep= iadv_PMf   ! Gives SO4 dep for OMp
             if ( ik == 0 ) fac = (1 - Fgas ) * cfac(iadv_PMf,i,j)
           else  ! keeps fac=1.0, need to consider dry dep of gas vs particle
              if(ik==0) fac = (1 - Fgas ) * cfac(iadv_PMf,i,j) + & !  PM term
                                   Fgas   * cfac(iadv,i,j)         ! gas-term
           end if

           if ( dbgPt  ) write(*,'(a,i4,2f8.4,a)')dtxt//'IOM_mix',itot, &
                fac, Fgas, trim(species(itot)%name )

        else ! Simple gas or particle.

            if(ik==0) then
               if(iadvDep<1) then
                 print *, dtxt//"IADVDEP ", iadv, itot, nspec, &
                     ParticlePhaseOutputs, trim(unit), semivol_wanted
                 call StopAll(dtxt//'IADVDEP problem')
               end if
               fac = fac * cfac(iadvDep,i,j)
            end if
        
        end if


        !if(all([first_semivol_call,debug_proc,chemgroups(igrp)%name=='BSOA']))&
        if ( dbgPt .and. first_call  ) &
          write(*,"(2(1x,a25),L2,2i4,1x, 2es12.3, f12.5)") &
            dtxt//"GRP fac "//trim(chemgroups(igrp)%name), &
             trim(species(itot)%name), semivol_wanted, nspec, itot, &
              cfac(iadvDep,i,j), fac, & !n.b. can't print Fgas3d for most itot
              xn_adv(iadv,i,j,kk)  * gunit_conv(nspec) * fac

        g2d(i,j) = g2d(i,j) + xn_adv(iadv,i,j,kk)  * gunit_conv(nspec) * fac
      end do ! nspec
      if( first_semivol_call .and. semivol_wanted) first_semivol_call = .false.
      first_call = .false.
    end do ! i
  end do ! j


  if(needroa)&
    forall(i=1:limax,j=1:ljmax) &
      g2d(i,j) = g2d(i,j) * density(i,j)
  deallocate(gspec,gunit_conv)
end subroutine group_calc
!=========================================================================
subroutine somo_calc( n, iX, debug_flag )
!/-- Calculates SOMO (8hours) values for input threshold.
!    which is max. value per day of 8-hour running means

  implicit none
  integer, intent(in) :: n           ! index in Derived_mod::d_2d arrays
  integer, intent(in) :: iX          !  threshold, usually 35 ppb
  logical, intent(in) :: debug_flag

  real :: o3                         ! Ozone (ppb) - needed if SOMOs
  real :: sum8h
  integer, parameter :: N8h = (NTDAY*8)/24 !number of periods in 8 hours
  real, parameter :: N8h_inv=1./N8h
  integer :: nh

  do i=1,limax
    do j=1,ljmax

      !find max running 8h sum O3
      sum8h=0.
      do nh=1,N8h
        sum8h = sum8h + D2_O3_DAY( i , j , nh)
      end do
      o3=sum8h
      do nh=N8h+1,NTDAY
        sum8h =sum8h-D2_O3_DAY( i , j , nh-N8h)+D2_O3_DAY( i , j , nh)
        o3=max(o3,sum8h)
        if(n<0)write(*,*)o3 !pw fake for compiler!!
      end do

      !divide by N8h to find 8h mean
      o3=o3*N8h_inv

      if(debug_flag.and.i==debug_li.and.j==debug_lj)&
        write(*,"(a,2i4,f12.3)") "SOMO DEBUG ", n, iX, o3

      o3 = max( o3 - iX , 0.0 )   ! Definition of SOMOs

      ! d_2d values will be accumulated in Derived_mod
      d_2d(n, i,j,IOU_DAY ) = o3

    end do
  end do
end subroutine somo_calc
!=========================================================================
subroutine write_debugadv(n,ind,rho,txt,Is3D,indDep)
  integer, intent(in) :: n, ind
  real, intent(in) :: rho
  character(len=*) :: txt
  logical, intent(in), optional :: Is3D
  integer, intent(in), optional :: indDep ! Used for surrogate dep species
  integer :: indD
  logical :: Is3D_local
  character(len=*), parameter :: dtxt='DrvWrtAdv:'

  Is3D_local = .false.
  if(present(Is3D)) Is3D_local = Is3D
  indD = ind
  if(present(indDep)) indD = indDep
  if(Is3D_local)then
    k=1; l=lev3d(k)
    write(*,"(2a,3i4,2(1x,a),4f11.3)") dtxt//"PROC3D " , trim(txt) , n, &
       ind, indD, trim(f_3d(n)%name),trim(f_3d(n)%unit)  &
      ,d_3d(n,debug_li,debug_lj,k,IOU_INST)*f_3d(n)%scale &
      ,xn_adv(ind,debug_li,debug_lj,l)*f_3d(n)%scale &
      ,rho, cfac(indD,debug_li,debug_lj)
  else
    write(*,"(2a,3i4,2(1x,a),4f11.3)") dtxt//"PROC2D " , trim(txt) , n, &
       ind, indD, trim(f_2d(n)%name), trim(f_2d(n)%unit)  &
      ,d_2d(n,debug_li,debug_lj,IOU_INST)*f_2d(n)%scale &
      ,xn_adv(ind,debug_li,debug_lj,KMAX_MID)*f_2d(n)%scale &
      ,rho, cfac(indD,debug_li,debug_lj)
  end if
end subroutine write_debugadv
!=========================================================================
subroutine write_debug(n,ind,txt)
  integer, intent(in) :: n, ind
  character(len=*) :: txt

  write(*,"(2a,2i4,a,4g12.3)") "DERIV: GEN " , txt , n, ind  &
    ,trim(f_2d(n)%name),d_2d(n,debug_li,debug_lj,IOU_INST)
end subroutine write_debug
!=========================================================================
end module Derived_mod
