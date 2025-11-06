! <LocalFractions_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.6>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2025 met.no
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
! Emissions for primary source tracking can be summed up over splits (but not necessarily all,
! can be new/age or single species) for each source;
! otherwise emissions are first splitted before including into single species.
! Primary sources can be groups of species, otherwise single species.
! All fullchem species must have the same "sources" (Npos)
! Fullchem output must be specified
! Natural emissions have no sector/countries; can be primary or not. Only predefined names (BVOC, SoilNOx, Pollen, )


module LocalFractions_mod
!
! all subroutines for Local Fractions
!
use CheckStop_mod,     only: CheckStop,StopAll
use Chemfields_mod,    only: xn_adv, cfac, Fgas3d, PM25_water_rh50, x, xnew
use ChemDims_mod,      only: NSPEC_ADV, NSPEC_SHL,NSPEC_TOT,NEMIS_File
use ChemFunctions_mod, only: EC_AGEING_RATE, HydrolysisN2O5k
use ChemSpecs_mod,     only: species_adv,species
use ChemGroups_mod,    only: WDEP_OXN_GROUP
use Config_module,     only: MasterProc,NPROC,KMAX_MID,KMAX_BND,KCHEMTOP,USES, &
                             IOU_INST,IOU_HOUR_INST,IOU_HOUR,IOU_YEAR,IOU_MON,IOU_DAY,IOU_MAX_MAX, &
                             dt_advec,dt_advec_inv,RUNDOMAIN,runlabel1,HOURLYFILE_ending, &
                             lf_set,lf_src,lf_species,lf_country,lf_sector_groups,lf_spec_out, &
                             SO2_ix, O3_ix, NO2_ix, SO4_ix, NH4_f_ix, NO3_ix, NO3_f_ix, &
                             NO3_c_ix, NH3_ix, HNO3_ix, C5H8_ix, APINENE_ix, NO_ix, HO2_ix, OH_ix,&
                             HONO_ix,OP_ix,CH3O2_ix,C2H5O2_ix,CH3CO3_ix,C4H9O2_ix,MEKO2_ix,ETRO2_ix,&
                             PRRO2_ix,OXYO2_ix,C5DICARBO2_ix,ISRO2_ix,MACRO2_ix,TERPO2_ix,H2O2_ix,N2O5_ix, &
                             NATBIO, ASOC_ug1e3_ix, non_C_ASOA_ng1e2_ix
use Convection_mod,    only: convection_1d
use Country_mod,       only: MAXNLAND,NLAND,Country

use DefPhotolysis_mod, only: IDHONO,IDNO3,IDNO2
use EmisDef_mod,       only: NSECTORS,SECTORS,EMIS_FILE,NSECTORS_GNFR_CAMS,IS_TRAF,IS_POW, &
                             nlandcode,landcode,NCMAX,&
                             secemis, roaddust_emis_pot,KEMISTOP,&
                             EmisMaskIntVal,EmisMaskValues,EmisMaskIndex2Name,gridrcemis,&
                             NEmisMask, mask2name
use EmisGet_mod,       only: nrcemis, iqrc2itot, emis_nsplit,nemis_kprofile, emis_kprofile,&
                             emis_masscorr,  &  ! 1/molwt for most species
                             make_iland_for_time,itot2iqrc,iqrc2iem, emisfrac
use GridValues_mod,    only: dA,dB,xm2, dhs1i, glat, glon, projection, extendarea_N,i_fdom,j_fdom,&
                             RestrictDomain, GRIDWIDTH_M, xmd
use MetFields_mod,     only: ps,roa,EtaKz
use MPI_Groups_mod
use NetCDF_mod,        only: Real4,Real8,Out_netCDF,LF_ncFileID_iou,closedID,CloseNetCDF,&
                             GetCDF_modelgrid,masked_output
use OwnDataTypes_mod,  only: Deriv, Max_lf_sources, Max_lf_sectors, MAX_lf_country_group_size, &
                             MAX_lf_sector_group_size, Max_lf_spec, TXTLEN_NAME, TXTLEN_FILE, &
                             Max_lf_res, Max_lf_Country_list, Max_lf_sectors, Max_lf_Country_groups, &
                             Max_lf_out
use Par_mod,           only: me,LIMAX,LJMAX,MAXLIMAX,MAXLJMAX,gi0,gj0,li0,li1,lj0,lj1,GIMAX,GJMAX
use PhysicalConstants_mod, only : GRAV, AVOG, ATWAIR
use SmallUtils_mod,    only: find_index, key2str
use TimeDate_mod,      only: date, current_date,day_of_week,enddate
use TimeDate_ExtraUtil_mod,only: date2string
use My_Timing_mod,     only: Add_2timing, Code_timer, NTIMING
use VerticalDiffusion_mod, only: vertdiffn
use ZchemData_mod,only: rct, rcphot, xn_2d, rcemis, M, rcbio, Fgas

!(dx,dy,i,j) shows contribution of pollutants from (i+dx,j+dy) to (i,j)

implicit none
!external advection_mod_mp_vertdiffn_k

private

integer ::IC_STRATOS = 324566 !contribution from O3 stratosphere
integer ::IC_INIT = 324565 !contribution from pollutants present at the start of the run
integer ::IC_BC = 324564 !contribution from Boundary Conditions
integer ::IC_NAT = 324563 !contribution from BVOC and DMS (separately)
integer ::IC_BVOC = 324562 !contribution from BVOC(C5H8 and TERP)
integer ::IC_DMS = 324561 !contribution from DMS
!contributions from Pollen
integer ::IC_BIRCH = 324560
integer ::IC_OLIVE = 324559
integer ::IC_ALDER = 324558 
integer ::IC_RWEED = 324557 
integer ::IC_GRASS = 324556 
integer ::IC_MUGWORT1 = 324555 
integer ::IC_MUGWORT2 = 324554 
integer ::IC_MUGWORT3 = 324553 
integer ::IC_MUGWORT4 = 324552 
integer ::IC_MUGWORT5 = 324551 
logical, parameter :: DEBUG = .false.
logical, parameter :: DEBUGall = .false.

public  :: lf_init
public  :: lf_out
public  :: lf_av
public  :: lf_adv_x
public  :: lf_adv_y
public  :: lf_adv_k
public  :: lf_adv_k_2nd
public  :: lf_diff
public  :: lf_conv
public  :: lf_chem_emis_deriv
public  :: lf_chem
public  :: lf_chem_pre, lf_chem_mid, lf_chem_pos
public  :: lf_sia_pre, lf_sia_pos
public  :: lf_aqu_pre, lf_aqu_pos
public  :: lf_SurfArea_pre, lf_SurfArea_pos
public  :: lf_drydep, lf_wetdep, lf_PODnAOT
public  :: save_lf_emis ! save all emis for use by lf
public  :: lf_rcemis
public  :: lf_rcemis_nat
public  :: lf_saveall
public  :: lf_read
private  :: addsource
private  :: MeanMasksOut  ! makes also integral over each city mask

real, public, allocatable, dimension(:,:,:,:,:,:), save :: &
  loc_frac&    ! Fraction of pollutants that are produced locally, surrounding sources
  ,loc_frac_hour_inst&  !Houry local fractions
  ,loc_frac_hour&  !Houry average of local fractions
  ,loc_frac_day&  !Daily average of local fractions
  ,loc_frac_month&  !Monthly average of local fractions
  ,loc_frac_full  !Fullrun average of local fractions
real, public, allocatable, dimension(:,:,:,:,:), save :: &
     lf_src_acc ! accumulated local fraction over time periods
real, public, allocatable, dimension(:,:,:,:,:), save :: &
     lf_src_tot ! concentrations of pollutants used for Local Fractions
real, public, allocatable, dimension(:,:,:), save :: &
     lf_src_ps ! surface pressure output for Local Fractions

real, public, allocatable, dimension(:,:,:,:,:), save :: emis_lf_cntry
real, public, allocatable, dimension(:,:,:,:), save :: &
   loc_frac_src &   ! Fraction of pollutants that are produced locally, list of defined sources
  ,lf &   ! Fraction of pollutants that are produced locally, for all defined sources
  ,loc_frac_src_full &   ! Fraction of pollutants that are produced locally, list of defined sources
  ,lf_src_full   ! Fraction of pollutants that are produced locally, list of defined sources
real, public, allocatable, dimension(:,:,:), save :: &
  loc_frac_drydep  ! ddepositions per source (not fractions!)
real, public, allocatable, dimension(:,:,:), save :: &
  loc_frac_wetdep  ! wdepositions per source (not fractions!)
real, public, allocatable, dimension(:,:,:), save :: lf_PM25_water
real, public, allocatable, dimension(:,:,:,:), save :: D8M !
real, public, allocatable, dimension(:,:,:), save :: D8Max !
real, public, allocatable, dimension(:,:,:,:), save :: D8Max_av !
real, public, allocatable, dimension(:,:,:), save :: D8Max_6m_av ! 6 month, peak season
real, public, allocatable, dimension(:,:,:), save :: D1Max !
real, public, allocatable, dimension(:,:,:,:), save :: D1Max_av !
real, public, allocatable, dimension(:,:,:,:), save :: D8Max_av_ppb !
real, public, allocatable, dimension(:,:,:), save :: D8Max_6m_av_ppb ! 6 month, peak season
real, public, allocatable, dimension(:,:,:), save :: hourM !
real, public, allocatable, dimension(:,:,:,:), save :: AOT40 !
real, public, allocatable, dimension(:,:,:,:), save :: SOMO35 !
real, public, allocatable, dimension(:,:,:,:), save :: D8M_ppb !
real, public, allocatable, dimension(:,:,:), save :: D8Max_ppb !
real, public, allocatable, dimension(:,:,:), save :: hourM_ppb !
real, public, allocatable, dimension(:,:,:,:), save :: &
  loc_frac_1d  ! Fraction of pollutants without i or j and extended (0:limax+1 or 0:ljmax+1)
real, public, allocatable, dimension(:,:), save :: &
  loc_frac_src_1d  ! Fraction of pollutants without i or j and extended (0:limax+1 or 0:ljmax+1)
real, allocatable, save ::loc_poll_to(:,:,:,:,:)
real, allocatable, public, dimension(:,:), save ::xderiv !dX_ispec/dX_n
real, allocatable, public, dimension(:,:), save ::xderivSOA !dX_ispec/dX_n
real, allocatable, public, dimension(:,:), save ::ederiv !dX_ispec/demis_n
real, allocatable, public, dimension(:,:), save :: lf0_loc, lf0SOA_loc, lf_loc, lfSOA_loc
real, allocatable, public, dimension(:,:), save ::x_lf, xold_lf ,xnew_lf
real, allocatable, public, dimension(:,:,:,:,:), save ::Dchem_lf !may not be worth the cost?
real, allocatable, public, dimension(:,:,:,:,:), save ::xn_shl_lf!may not be worth the cost?
real, allocatable, public, dimension(:,:,:), save ::totrcemis_lf
real, allocatable, public, dimension(:,:), save ::rcemis_lf
real, allocatable, public, dimension(:,:), save ::rcemis_lf_surf
integer, allocatable, public, dimension(:,:), save ::emis2spec_surf
real, allocatable, public, dimension(:), save ::rcemis_lf_primary !for emissions considered as linear
integer, allocatable, public, dimension(:,:), save ::nic
integer, allocatable, public, dimension(:,:,:), save ::ic2iland
integer, allocatable, public, dimension(:,:), save :: lf_sector_map !add some sectors together
integer, allocatable, dimension(:), save :: lf_nsector_map !how many sector to put together
real, allocatable, public, dimension(:), save ::L_lf, P_lf!,rctA_lf,rctB_lf
!real, allocatable, public, dimension(:,:), save ::rctAk_lf,rctBk_lf
real, allocatable, public, dimension(:,:), save ::fgasso2_lf
real, allocatable, public, dimension(:,:,:), save ::AQRCK_lf
real, private, dimension(0:5,0:8),save ::xn_lf !to save concentrations

logical, public, save :: COMPUTE_LOCAL_TRANSPORT=.false.
integer , public, save :: lf_Nvertout = 1!number of vertical levels to save in output
integer, public, save :: NTIMING_lf=9
real, private :: tim_after,tim_before
integer, public, save :: Ndiv_rel=1
integer, public, save :: Nsources=0, Nsources_chem=0, Nsources_nonew=0
integer, public, save :: lf_Nvert=0

integer, public, save :: LF_SRC_TOTSIZE
integer, public, save :: iotyp2ix(IOU_MAX_MAX),ix2iotyp(IOU_MAX_MAX)
integer, public, save :: av_fac(IOU_MAX_MAX)
integer, public, save :: Niou_ix = 0 ! number of time periods to consider (hourly, monthly, full ...)
integer, private, save :: iou_ix_inst = -2 !set if hourly-instantaneous are requested
integer, public, save :: Npoll = 0 !Number of different pollutants to consider
integer, public, save :: iem2ipoll(NEMIS_File,Max_lf_spec) !internal indices of pollutants for that emis file
integer, public, save :: ipoll2iqrc(Max_lf_spec) = -1 !-1 for primary pollutant

integer, public, save :: Ndrydep_lf = 0
integer, public, save :: Nwetdep_lf = 0
logical, public, save :: wetdep_lf(NSPEC_ADV) = .false.

integer, private, save :: iem2Nipoll(NEMIS_File) !number of pollutants for that emis file
logical :: old_format=.false. !temporary, use old format for input and output
integer, private, save :: isrc_O3=-1, isrc_NO=-1, isrc_NO2=-1, isrc_VOC=-1
integer, private, save :: isrc_SO2=-1, isrc_SO4=-1, isrc_NH3=-1, isrc_NH3_nh3=-1
integer, private, save :: isrc_NO3=-1, isrc_NO3_nh3=-1, isrc_HNO3=-1, isrc_HNO3_nh3=-1
integer, private, save :: isrc_NH4_f=-1,isrc_NH4_f_nh3=-1,isrc_NO3_c=-1
integer, private, save :: isrc_SO4_f=-1,isrc_NO3_f=-1,isrc_NO3_f_nh3=-1
!NB: isrc_pm25_new is overwritten many times. TODO: clean up. For fullchem there is only one value it can take, but not otherwise.
integer, private, save :: isrc_pm25_new=-1, isrc_pm25=-1, isrc_pmco=-1
integer, private, save :: isrc_strato=-1, isrc_ini=-1
integer, private, save :: RO2POOL_ix=-1
real, allocatable, private, save :: lf_NH4(:), lf_NH3(:)
real, allocatable, private, save :: lf_NO3(:), lf_HNO3(:)
integer, private, save :: country_ix_list(Max_lf_Country_list)
integer, private, save :: Ncountry_lf = 0
integer, private, save :: Ncountry_group_lf = 0
integer, private, save :: Nsector_groups_lf = 0
integer, private, save :: Ncountrysectors_lf = 0
integer, private, save :: Ncountry_mask_lf=0 !total number of masks defined
integer, private, save :: Ncountry_mask_lf_val=0 !number of masks defined using lf_country%mask_val
integer, private, save :: country_mask_val(Max_lf_Country_list) = -999999 ! values of all defined masks
integer, private, save :: iic2ilf_countrymask(Max_lf_Country_list) = -1
character(len=TXTLEN_NAME), private, save :: iem2names(NEMIS_File,Max_lf_spec) !name of that pollutant
integer, private, save :: isrc_new(Max_lf_sources)
integer, private, save :: Stratos_ix(1000) !1000 must be larger than Nsources
integer, private, save :: nstratos !number of sources to track for Stratos
integer, private, save :: BC_ix(1000) !1000 must be larger than Nsources
integer, private, save :: nbc !number of sources to track for Boundary Conditions
real   , private, save :: P_NO(100),P_NO2(100)
real   , parameter     :: eps1 = 0.999
real   , parameter     :: eps1_sia = 0.95
real   , parameter     :: eps1_sia_m = 1.0+1.0-eps1_sia
integer, public        :: Nsia_deriv = 7
logical, public, save :: lf_fullchem=.false. ! if the full O3 chemistry is needed
integer, public, save :: NSPEC_fullchem_lf=0 ! number of species to include in the "fullchem" derivatives
integer, public, parameter :: N_lf_derivemisMAX = 200 ! max number of emission source to include in CM_Reactions1 derivatives
integer, public, save :: N_lf_derivemis = 0! actual number of emissions to include in CM_Reactions1 derivatives
integer, public, save :: nemis_primary = 0! number of primary emissions to include in this gridcell
integer, public, save :: Nemis_surf = 0! number of non sector surface emissions in this gridcell
integer, private, save :: iem_nox, iem_voc , iem_nh3, iem_sox !index in EMIS_FILE
integer, public, parameter :: iem_lf_nox = 1, iem_lf_voc = 2, iem_lf_nh3 = 3, iem_lf_sox = 4
integer, public, save :: emis2icis(N_lf_derivemisMAX),emis2pos_primary(N_lf_derivemisMAX)
integer, public, save :: emis2is(N_lf_derivemisMAX)
integer, public, save :: emis2isrc_primary(N_lf_derivemisMAX) ! connect the index in rcemis_lf with the isrc (source). Only make sense for primary
integer, public, save :: emis2isrc(N_lf_derivemisMAX),emis2iem(N_lf_derivemisMAX)
integer, public, save :: emis2iic_surf(N_lf_derivemisMAX), emis2nspec_surf(N_lf_derivemisMAX)
integer, public, save :: lfspec2spec(NSPEC_TOT),spec2lfspec(NSPEC_TOT) !mapping between LF species index and the index from CM_Spec (tot)
integer, public, save :: Nlf_species = 0, NSPEC_chem_lf = 0, NSPEC_deriv_lf, N_deriv_SOA_lf = 0, NSOA
integer, public, save :: Nfullchem_emis = 1 !4 if nox, voc, nh3, sox separately or 1 if all together, or 2 if only nox and voc
integer, public, save :: ix_lf_max
integer, public, save :: Npos_lf
logical, public, save :: makeDMS = .false. ! Each natural emission to track has an ad hoc variable , makeXXX
logical, public, save :: makeFungal = .false. ! Each natural emission to track has an ad hoc variable , makeXXX
logical, public, save :: makePollen = .false. ! Each natural emission to track has an ad hoc variable , makeXXX
logical, private, save :: make_PMwater =.false. !NB: water is not a species and will be treated separately. Index isrc=NSOURCES+1, or %start=LF_SRC_TOTSIZE+1
logical, public, save :: makeBVOC = .false.
integer, public, save :: ix_BVOC, ix_DMS
integer, public, save :: ix_BIRCH, ix_OLIVE, ix_ALDER, ix_RWEED, ix_GRASS, ix_MUGWORT1, ix_MUGWORT2, ix_MUGWORT3, ix_MUGWORT4, ix_MUGWORT5
integer, public, save :: nPODnAOT=0, nDryDep=0 !number of outputs asked for. nDryDep includes nPODnAOT 
logical, private, save:: aero_error = .false.
real, parameter :: lf_limit = 1e-5
real, parameter :: lf_limit_sia = 1e-5
integer, parameter :: NAQUEOUS = 5,ICLOHSO2 = 1,ICLRC1 = 2,ICLRC2 = 3,ICLRC3 = 4,ICLHO2H2O2 = 5 !NB: hardcoded!!

contains

  subroutine lf_init
    integer :: n, n0, is, i, j, ii, iii, ic, ir, ix, iix, isrc, n_mask, mask_val_min, mask_val_max
    integer :: found, itot, iqrc, iem, iemis, ipoll, ixnh3, ixnh4, size, IOU_ix, iem_deriv
    integer :: iout, ig, idep, isec, isec_lf
    integer, allocatable :: MaskVal(:)
    logical is_relative,is_country
    character(len=200) :: filename
! pm25_new and pm25 are considered as two different emitted pollutants
    if(DEBUGall .and. me==0)write(*,*)'start init'

  call Code_timer(tim_before)
  ix=0

  !temporary backward compatibility
  lf_set%YEAR=lf_src(1)%YEAR .or. lf_set%YEAR
  lf_set%MONTH=lf_src(1)%MONTH .or. lf_set%MONTH
  if (trim(lf_src(1)%MONTH_ENDING) /= "NOTSET") lf_set%MONTH_ENDING=lf_src(1)%MONTH_ENDING
  lf_set%DAY=lf_src(1)%DAY .or. lf_set%DAY
  lf_set%HOUR=lf_src(1)%HOUR .or. lf_set%HOUR
  lf_set%HOUR_INST=lf_src(1)%HOUR_INST .or. lf_set%HOUR_INST

  lf_Nvert = lf_set%Nvert !Temporary
  if (lf_Nvert>KMAX_MID-1) then
     lf_Nvert = KMAX_MID-1
     if(me==0)then
        write(*,*)'WARNING: cannot track through level 1 (top). Reducing Nvert to ',lf_Nvert
     end if
  end if
  if (lf_set%dist>-1) then
     lf_src(:)%dist = lf_set%dist
  else
     lf_set%dist = lf_src(1)%dist
  end if
 
  if (lf_set%full_chem) lf_fullchem = .true.
  if (lf_set%Nvertout > lf_Nvert .and. me==0 .and. lf_fullchem ) write(*,*)'lf fullchem multi vertical level not implemented '
  if (lf_fullchem ) lf_set%Nvertout = 1
  if (lf_set%Nvertout > lf_Nvert .and. me==0) write(*,*)'WARNING: will only output ',lf_Nvert,' vertical levels'
  lf_set%Nvertout = min(lf_Nvert, lf_set%Nvertout)
  lf_Nvertout = lf_set%Nvertout
  if (me==0) write(*,*)'output ',lf_Nvertout,' vertical levels'

  Nsources = 0
  if (lf_fullchem) then
     !We ASSUME that shipNOx is the last species (highest index) involved in O3 chemistry
     ix=find_index("shipNOx" ,species_adv(:)%name)
     if(ix>0)NSPEC_fullchem_lf = ix
     ix_lf_max = NH3_ix !assumed to be last used index for rcemis

     !RO2POOL is not a proper chemical species and we must take care to avoid it when using derivatives!
     RO2POOL_ix = find_index('RO2POOL' ,species(:)%name)

     !crude check that SOA are as expected
     call CheckStop(non_C_ASOA_ng1e2_ix-ASOC_ug1e3_ix/=9, "LF: unexpected SOA indices. Revise LF reactions")

  else
     do i = 1, Max_lf_sources
        if(lf_src(i)%name == 'DMS') then
           lf_src(i)%species = 'SO2'
           lf_src(i)%is_NATURAL = .true.
           makeDMS= .true.
        end if
        if(lf_src(i)%name == 'FUNGAL_SPORES') then
         lf_src(i)%species = 'FUNGAL_SPORES'
         lf_src(i)%is_NATURAL = .true.
         makeFungal= .true.
        end if
        if(lf_src(i)%name == 'POLLEN_BIRCH') then
           lf_src(i)%species = 'POLLEN_BIRCH'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_OLIVE') then
           lf_src(i)%species = 'POLLEN_OLIVE'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_ALDER') then
           lf_src(i)%species = 'POLLEN_ALDER'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_RWEED') then
           lf_src(i)%species = 'POLLEN_RWEED'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_GRASS') then
           lf_src(i)%species = 'POLLEN_GRASS'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_MUGWORT1') then
           lf_src(i)%species = 'POLLEN_MUGWORT1'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_MUGWORT2') then
           lf_src(i)%species = 'POLLEN_MUGWORT2'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_MUGWORT3') then
           lf_src(i)%species = 'POLLEN_MUGWORT3'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_MUGWORT4') then
           lf_src(i)%species = 'POLLEN_MUGWORT4'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        else if(lf_src(i)%name == 'POLLEN_MUGWORT5') then
           lf_src(i)%species = 'POLLEN_MUGWORT5'
           lf_src(i)%is_NATURAL = .true.
           makePollen= .true.
        end if
        if (lf_src(i)%species == 'NOTSET') exit
        Nsources = Nsources + 1
     enddo
  end if

  do i = Nsources + 1, Max_lf_sources
     if(lf_src(i)%species /= 'NOTSET') then
        if(me==0)write(*,*)'WARNING: lf_src ',i,' ',trim(lf_src(i)%species),' not included because source ',Nsources+1,' is missing'
     end if
  enddo

  !we includes sources defined using lf_species nomenclature
  do i = 1, Max_lf_spec
     if(lf_species(i)%name == 'NOTSET') exit
     do ii = 1, Max_lf_sectors
        if(lf_species(i)%sectors(ii) < 0) exit
        do iii = 1, Max_lf_res
           if(lf_species(i)%res(iii) < 0) exit
           Nsources = Nsources + 1
           lf_src(Nsources)%species = lf_species(i)%name
           lf_src(Nsources)%sector = lf_species(i)%sectors(ii)
           lf_src(Nsources)%res = lf_species(i)%res(iii)
        end do
     end do
  end do

  !we add sources corresponding to full chemistry.
  lfspec2spec = -1!meaning: not a lf species
  spec2lfspec = -1!meaning: not a lf species
  do n=1, NSPEC_SHL
     spec2lfspec(n)=0!meaning: defined, but not a tracked species
  end do
  if (lf_fullchem) then
     if (lf_set%Nfullchem_emis > 0) then
        Nfullchem_emis = lf_set%Nfullchem_emis
        if (Nfullchem_emis == 1) lf_set%EmisDer_all = .true. ! to improve? could be only nox for example
        if (Nfullchem_emis == 3 .or. Nfullchem_emis >4  ) call StopAll("Nfullchem_emis must be 1, 2 or 4")
     else
        !for backward compatibility
        if (lf_set%EmisDer_all) then
           Nfullchem_emis = 1
        else
           Nfullchem_emis = 4
        end if
     end if
     iem = find_index('nox' ,EMIS_FILE(1:NEMIS_FILE))
     call CheckStop(iem<=0, "LF: did not find nox emissions")
     iem_nox = iem
     iem = find_index('voc' ,EMIS_FILE(1:NEMIS_FILE))
     call CheckStop(iem<=0, "LF: did not find voc emissions")
     iem_voc = iem
     iem = find_index('nh3' ,EMIS_FILE(1:NEMIS_FILE))
     call CheckStop(iem<=0, "LF: did not find nh3 emissions")
     iem_nh3 = iem
     iem = find_index('sox' ,EMIS_FILE(1:NEMIS_FILE))
     call CheckStop(iem<=0, "LF: did not find sox emissions")
     iem_sox = iem
     do i = 1, NSPEC_fullchem_lf !make sources for each species to be fully included in chemistry
        call addsource(species_adv(i)%name)
     end do
     
     !not directly used in Reactions1, but have an effect on aqrck
     call addsource('SO4')
     call addsource('NO3_c')
     call addsource('NO3_f') ! source index must be after HNO3 (for lf_sia)!
     call addsource('NH3')
     call addsource('NH4_f') ! source index must be after NH3 (for lf_sia)!


     NSPEC_deriv_lf = Nlf_species !the number of species "Y" in dX/dY for the O3 block
     !additional species Z for which dO3/dZ = 0: NH3, SOA (SO4, NO3_c see above)
     !SOA all ASOC and ASOA
     do i = ASOC_ug1e3_ix-NSPEC_SHL, non_C_ASOA_ng1e2_ix-NSPEC_SHL !make sources for each SOA species
        call addsource(species_adv(i)%name)
     end do

     N_deriv_SOA_lf  = Nlf_species !the number of species "Y" in dX/dY used for SOA
     NSOA = N_deriv_SOA_lf - NSPEC_deriv_lf

     NSPEC_chem_lf = Nlf_species !the number of species "X" in dX/dY

     !end of sources to be included in chemistry
     Nsources_chem = Nsources ! each species can have several sources (voc, nox, nh3, sox)
     
     ! Additional species used for equilibrium chemistry, but not used in chemistry:  none

     ! add relative species. Each "neighbor" count as one country to track
     ! must account for emissions and chemistry
     ! emissions should (directly) only change middle cell lf
     ! chemistry: lf will change in the same way as a country, just increase Npos_lf with number of surrounding cells + self.
     ! no new species to track, but new "countries"!
     !NB: for non-chem situation, a new isrc is defined for each sector, while in the fullchem case, sectors are part of isrc
    
     !for each relative (2*lf_set%dist+1)*Ncountrysectors_lf new sources are defined
     if (lf_set%relative .and. lf_set%dist>-1) then
        Ndiv_rel = 2*lf_set%dist+1
        if(me==0)write(*,*)"WARNING: all outputs (also SURF) are NOT surface corrected, i.e. they are NOT for 3 m height)"
     else
        Ndiv_rel = 0
     end if
     ! add primary PM
     call addsource('pm25')
     call addsource('pmco')

     ! sizes:
     ! number of species that are taken systematically by indices: NSPEC_fullchem_lf
     ! including additional species included, but not derived wrt: NSPEC_chem_lf
     ! short-lived are not advected, but included as species to derive. They can also produce shift of indices:NSPEC_SHL
     ! other species included in equilibrium chemistry (aerosols), but not present in "normal" chemistry: Nlf_species
     ! highest species index used in chemistry (rcemis_lf): ix_lf_max
     if(me==0) write(*,"((X,A,':',I4))") &
        'LF chemistry, number of chem derivatives calculated excl SOA', NSPEC_deriv_lf, &
        'LF chemistry, additional derivatives for SOA', NSOA, &
        'LF chemistry, number of emis derivatives calculated', Nfullchem_emis, &
        'LF chemistry, species included', NSPEC_chem_lf, &
        'LF total species included', Nlf_species + NSPEC_SHL, &
        'LF chemistry, last index included', ix_lf_max
  end if


  !we define additional sources that are defined at different times
  Nsources_nonew = Nsources !for not changing loop bound in the loop
  do i = 1, Nsources_nonew
     if (lf_src(i)%is_NATURAL) then
        ix=find_index(lf_src(i)%species ,species(:)%name, any_case=.true.)
        if (ix<0) then
           call StopAll('LF: NATURAL only defined for single species defined in CM_ChemSpecs_mod.f90')
        else
           lf_src(i)%species = trim(species(ix)%name) !will correct possible case lower/upper differences
           lf_src(i)%species_ix = ix
        end if
     end if
     if (lf_src(i)%nhour>0) then
        isrc = i
        do n = 0, 23, lf_src(i)%nhour
           if (n > 0) then
              !new source, same definition, but different emis time
              if(me==0)write(*,*)'defining new source for time ', n
              Nsources = Nsources + 1
              isrc = Nsources
              lf_src(isrc)%name = lf_src(i)%name
              lf_src(isrc)%species = lf_src(i)%species
              lf_src(isrc)%species_ix = lf_src(i)%species_ix
              lf_src(isrc)%dist = lf_src(i)%dist
              lf_src(isrc)%res = lf_src(i)%res
              lf_src(isrc)%type = lf_src(i)%type
              lf_src(isrc)%is_NATURAL = lf_src(i)%is_NATURAL
              lf_src(isrc)%nhour = lf_src(i)%nhour
           end if
           lf_src(isrc)%time_ix = n
        end do
     end if
  end do

  !for each pm25 we separate into new and age parts
  Nsources_nonew = Nsources

  isrc_new = -1
  do i = 1, Nsources_nonew
     if(lf_src(i)%species == 'pm25')then
        if(MasterProc)write(*,*)'splitting pm25 for source',i,' into pm25 and pm25_new'
        Nsources = Nsources + 1
        call CheckStop(Nsources>Max_lf_sources,"Number of LF sources exceeds Max_lf_sources")
        lf_src(Nsources) = lf_src(i)
        lf_src(Nsources)%species = 'pm25_new'
        isrc_new(i) = Nsources
        isrc_pm25_new = Nsources
     end if
  enddo

  !countries
  ! The lowest indices are for mask indices, the highest indices are for groups
  if (lf_country%mask_val_min <= lf_country%mask_val_max .or. &
       lf_country%mask_val(1) > -999999 .or. &
       lf_country%list(1)/= 'NOTSET' .or. &
       lf_country%group(1)%name/= 'NOTSET'.or. &
       lf_country%cellmask_name(1)/= 'NOTSET') then

     Ncountry_mask_lf = 0
     Ncountry_mask_lf_val = 0
     Ncountry_lf=0
     mask_val_min = lf_country%mask_val_min
     mask_val_max = lf_country%mask_val_max
     do i = 1, Max_lf_Country_list
        if(lf_country%mask_val(i) < -999999) exit
        mask_val_min = min(lf_country%mask_val(i),mask_val_min)
        mask_val_max = max(lf_country%mask_val(i),mask_val_max)
     end do
     if (mask_val_max >= mask_val_min) then
        !make listof mask values that exist in EmisMaskIntVal from any MPI
        allocate(MaskVal(mask_val_min:mask_val_max))
        MaskVal = 0
        do j=1,ljmax
           do i=1,limax
              if (EmisMaskIntVal(i,j)>=mask_val_min .and. EmisMaskIntVal(i,j)<=mask_val_max) then
                 MaskVal(EmisMaskIntVal(i,j))=1
              end if
         enddo
        enddo
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaskVal,mask_val_max-mask_val_min+1,MPI_INTEGER, MPI_SUM,MPI_COMM_CALC,IERROR)
        do j=1,ljmax
           do i=1,limax
              if (EmisMaskIntVal(i,j)>=mask_val_min .and. EmisMaskIntVal(i,j)<=mask_val_max) then
                 MaskVal(EmisMaskIntVal(i,j))=1
              end if
           enddo
        enddo
     end if
     do i = 1, Max_lf_Country_list
        if(lf_country%mask_val(i) < -999999) exit

        if(MaskVal(lf_country%mask_val(i)) == 0 ) cycle !the value is not defined anywhere on the netcdf masks

        Ncountry_mask_lf = Ncountry_mask_lf + 1
        Ncountry_lf = Ncountry_lf + 1
        country_mask_val(Ncountry_lf) = lf_country%mask_val(i)
     end do
     Ncountry_mask_lf_val = Ncountry_lf !only defined with lf_country%mask_val (not min/max)
     if (Ncountry_mask_lf_val>0 .and. MasterProc) then
        write(*,*)'including  ',Ncountry_mask_lf_val,' individually defined mask sources '
     end if

     do i = lf_country%mask_val_min,lf_country%mask_val_max
        if (MaskVal(i)==0) cycle ! is not anywhere on the mask
        found = 0
        do n = 1, Ncountry_mask_lf_val
           if (i == lf_country%mask_val(n)) found = 1
           if (found == 1) exit
        end do
        if (found == 1) cycle !already included
        Ncountry_mask_lf = Ncountry_mask_lf + 1
        Ncountry_lf = Ncountry_lf + 1
        country_mask_val(Ncountry_lf) = i
        if(i>size(mask2name))then
           write(mask2name(i),fmt='(A)')i
        else
           if(mask2name(i)=='NOTSET')write(mask2name(i),fmt='(A)')i
        end if
     end do
     if (Ncountry_mask_lf>0 .and. MasterProc) then
        write(*,*)'including ',Ncountry_mask_lf,'non fractional mask sources:'
        do n = 1, (Ncountry_mask_lf+29)/30
           write(*,fmt="(30(I0,1x))")(country_mask_val(ii),ii=(n-1)*30+1, min(Ncountry_mask_lf,n*30))
        end do
     end if

     if (mask_val_max >= mask_val_min) deallocate(MaskVal)

     !define masks with "fraction of cell" method
     do i = 1, Max_lf_Country_list
        if (lf_country%cellmask_name(i) /= 'NOTSET') then
           !find a defined mask with that name
           ii = find_index(trim(lf_country%cellmask_name(i)), EmisMaskIndex2Name(:), any_case=.true.)
           if (ii>0) then
              Ncountry_lf = Ncountry_lf + 1
              Ncountry_mask_lf = Ncountry_mask_lf + 1
              iic2ilf_countrymask(Ncountry_lf) = ii !this is the internal index use by lf routine for that "country"
              country_mask_val(Ncountry_lf) = i
              mask2name(i) = lf_country%cellmask_name(i)
              if (me==0) write(*,*)Ncountry_lf,ii,'LF will include cell mask for '//trim(lf_country%cellmask_name(i))
           else
              if (me==0) write(*,*)'WARNING: LF did not find defined cell fraction mask for '//trim(lf_country%cellmask_name(i))
           end if
        end if
     end do
     if (Ncountry_mask_lf>0 .and. MasterProc) write(*,*)'including in total',Ncountry_mask_lf,' masks'

     if(lf_country%list(1)/= 'NOTSET' .or. lf_country%group(1)%name/= 'NOTSET')then
        !list of countries/sectors instead of single country
        do i = 1, Max_lf_Country_list
           if(lf_country%list(i) == 'NOTSET') exit
           Ncountry_lf=Ncountry_lf+1
           ix = find_index(trim(lf_country%list(i)) ,Country(:)%code, first_only=.true.)
           if(ix<0 .and. lf_country%list(i) =='STRATOS')then
              ix = IC_STRATOS
           else if(ix<0 .and. lf_country%list(i) =='LBC')then
              ix = IC_BC
           else if(ix<0 .and. lf_country%list(i) =='INIT')then
              ix = IC_INIT
           else if(ix<0 .and. lf_country%list(i) =='DMS')then
              ix = IC_DMS
              makeDMS = .true.
              ix_DMS = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='BVOC')then
              ix = IC_BVOC
              call CheckStop(C5H8_ix<0 .or. APINENE_ix<0,&
                   'country BVOC cannot be computed, because C5H8 or APINENE not found ')
              makeBVOC = .true.
              ix_BVOC = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='BIRCH')then
              ix = IC_BIRCH
              makePollen = .true.
              ix_BIRCH = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='OLIVE')then
              ix = IC_OLIVE
              makePollen = .true.
              ix_OLIVE = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='ALDER')then
              ix = IC_ALDER
              makePollen = .true.
              ix_ALDER = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='RWEED')then
              ix = IC_RWEED
              makePollen = .true.
              ix_RWEED = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='GRASS')then
              ix = IC_GRASS
              makePollen = .true.
              ix_GRASS = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='MUGWORT1')then
              ix = IC_MUGWORT1
              makePollen = .true.
              ix_MUGWORT1 = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='MUGWORT2')then
              ix = IC_MUGWORT2
              makePollen = .true.
              ix_MUGWORT2 = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='MUGWORT3')then
              ix = IC_MUGWORT3
              makePollen = .true.
              ix_MUGWORT3 = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='MUGWORT4')then
              ix = IC_MUGWORT4
              makePollen = .true.
              ix_MUGWORT4 = i + Ncountry_mask_lf
           else if(ix<0 .and. lf_country%list(i) =='MUGWORT5')then
              ix = IC_MUGWORT5
              makePollen = .true.
              ix_MUGWORT5 = i + Ncountry_mask_lf           
           end if
           call CheckStop(ix<0,'country '//trim(lf_country%list(i))//' not defined. ')
           country_ix_list(i + Ncountry_mask_lf) = ix
           if(MasterProc)write(*,*)'include sources from ',trim(lf_country%list(i))
        enddo
     end if

     Ncountry_group_lf=0
     do i = 1, Max_lf_Country_groups
        if(lf_country%group(i)%name == 'NOTSET') exit
        Ncountry_group_lf = Ncountry_group_lf+1
        do ic = 1, MAX_lf_country_group_size
           if(lf_country%group(i)%list(ic) == 'NOTSET') exit
           ix = find_index(trim(lf_country%group(i)%list(ic)) ,Country(:)%code, first_only=.true.)
           call CheckStop(ix<0,'country '//trim(lf_country%group(i)%list(ic))//' not defined. ')
           lf_country%group(i)%ix(ic) = ix
           if(MasterProc)write(*,*)'include sources from '//&
                trim(lf_country%group(i)%list(ic))//' as '//trim(lf_country%group(i)%name)
        enddo
     end do

     Ncountrysectors_lf=0
     do i = 1, Max_lf_sectors
        if(lf_country%sector_list(i) < 0) then
           if(i==1 .and. trim(lf_sector_groups(i)%name) == 'NOTSET')then
              !default to only sector 0
              lf_country%sector_list(i) = 0
           else
              exit
           end if
        end if
        Ncountrysectors_lf = Ncountrysectors_lf + 1
        if(MasterProc)write(*,*)'country sector ',lf_country%sector_list(i)
     end do
  end if
  
  allocate(lf_sector_map(MAX_lf_sector_group_size,0:NSECTORS+Max_lf_sectors),lf_nsector_map(0:NSECTORS+Max_lf_sectors))
  !note: the loop above, is interrputed by an exit, and cannot be used
  do i = 1, NSECTORS
     lf_sector_map(:,i) = i
  end do
  lf_nsector_map(:) = 1
  !sector zero is all sectors
  do i=1, NSECTORS
     lf_sector_map(i,0) = i
  end do
  lf_nsector_map(0) = NSECTORS
  !hardcoded for now:
  if (Ncountrysectors_lf<NSECTORS .and. NSECTORS == NSECTORS_GNFR_CAMS) then
     if(me==0)write(*,*)'LF TRAF sectors 16, 17, 18, 19 mapped into sector 6'
     if(me==0)write(*,*)'LF POW sectors 14, 15 mapped into sector 1',lf_sector_map(1,7)
     !we map all transport sectors (F1...F4) into F
     if(.not. IS_TRAF(6) .or. .not. IS_TRAF(16) .or. .not. IS_TRAF(17) .or.  .not. IS_TRAF(18).or.  .not. IS_TRAF(19))then
        call StopAll('SECTOR HACK does not work for this setup!')
     end if
     lf_sector_map(1,6) = 6
     lf_sector_map(2,6) = 16
     lf_sector_map(3,6) = 17
     lf_sector_map(4,6) = 18
     lf_sector_map(5,6) = 19
     lf_nsector_map(6) = 5
     !we map all transport sectors (A1, A2) into A
     if(.not. IS_POW(1) .or. .not. IS_POW(14) .or. .not. IS_POW(15))then
        call StopAll('SECTOR HACK does not work for this setup!')
     end if
     lf_sector_map(1,1) = 1
     lf_sector_map(2,1) = 14
     lf_sector_map(3,1) = 15
     lf_nsector_map(1) = 3
  end if
  
  !sectors. Define all as groups
  !first predefined groups from config
  Nsector_groups_lf = 0
  do i = 1, Max_lf_sectors
     if (trim(lf_sector_groups(i)%name) == 'NOTSET') exit
     if (lf_sector_groups(i)%list(1) < 0 .or. lf_sector_groups(i)%list(1) > NSECTORS) then
        if (me==0) write(*,*)'WARNING: did not find any valid sector in '//trim(lf_sector_groups(i)%name),lf_sector_groups(i)%list(1)
        exit
     end if
     Nsector_groups_lf = Nsector_groups_lf + 1
     Ncountrysectors_lf = Ncountrysectors_lf + 1
     !redefine groups as pseudo sectors, with index larger than NSECTORS
     isec_lf = NSECTORS + Nsector_groups_lf ! pseudo sector defined as a group of sectors
     lf_country%sector_list(Ncountrysectors_lf) = isec_lf
     lf_sector_groups(i)%nsec = 0
     lf_nsector_map(isec_lf) = 0
     do n = 1, MAX_lf_sector_group_size
        if (lf_sector_groups(i)%list(n) > 0 .and. lf_sector_groups(i)%list(n) <= NSECTORS) then
           lf_sector_groups(i)%nsec = lf_sector_groups(i)%nsec + 1
           lf_nsector_map(isec_lf) = lf_nsector_map(isec_lf) + 1 !number of sectors included this group
           lf_sector_map(n,isec_lf) = lf_sector_groups(i)%list(n)! sector included this group
           if(me==0)write(*,*)n,' LF defined group ',isec_lf,' includes sector ',lf_sector_groups(i)%list(n)
        else
           exit
        end if
     end do
  end do
  
  Npos_lf = (Ncountry_lf + Ncountry_group_lf + Ndiv_rel*Ndiv_rel)*Ncountrysectors_lf
  
  !TODO: add only one sector and deriv for STRATOS, INIT and BVOC
  do isrc = 1, Nsources
     if(lf_fullchem)then
        !each species has same content
        lf_src(isrc)%Npos = Npos_lf
     else
        !country and relative ave different size
        if(lf_src(isrc)%type == 'country')then
           lf_src(isrc)%Npos = (Ncountry_lf + Ncountry_group_lf)*Ncountrysectors_lf
        else if(lf_src(isrc)%type == 'relative')then
           lf_src(isrc)%Npos = Ndiv_rel*Ndiv_rel
        else
           call StopAll("lf_src(isrc)%type not recognized ")
        end if
     end if
  end do
  if(MasterProc)write(*,*)(Ncountry_lf + Ncountry_group_lf )*Ncountrysectors_lf,' countries x sectors for ',Nsources,' sources '
  if(MasterProc .and. Ndiv_rel > 0)write(*,*)(Ndiv_rel*Ndiv_rel)*Ncountrysectors_lf,' local (relative) source gridcells'
  
  
  ipoll=0
  iem2ipoll = -1
  iem2Nipoll = 0
  do isrc = 1, Nsources
     if(lf_src(isrc)%type == 'relative' .and. .not. lf_fullchem)then
        !for now only one Ndiv possible for all sources
        lf_src(isrc)%Npos =  (2*lf_src(isrc)%dist+1)*(2*lf_src(isrc)%dist+1)
        Ndiv_rel = max(Ndiv_rel,2*lf_src(isrc)%dist+1)
        if(me==0)write(*,*)isrc,'relative: ',trim(lf_src(isrc)%species),lf_src(isrc)%Npos
     end if

     if(lf_src(isrc)%country_ISO /= 'NOTSET')then
        lf_src(isrc)%type = 'country'
        ix = find_index(trim(lf_src(isrc)%country_ISO) ,Country(:)%code, first_only=.true.)
        if(ix<0)then
           if(me==0)write(*,*)'LF: WARNING: country '//trim(lf_src(isrc)%country_ISO)//' not defined. '
        end if
        lf_src(isrc)%country_ix = ix
        if(MasterProc)write(*,*)isrc,' country '//trim(lf_src(isrc)%country_ISO)//' '//trim(lf_src(isrc)%species)
     end if

     call RestrictDomain(lf_set%DOMAIN)

     iem=find_index(lf_src(isrc)%species ,EMIS_FILE(1:NEMIS_FILE))

     if(iem<1)then
        if(lf_src(isrc)%species=='pm25_new')then
           !we separate pm25_new from other pm25
           iem = find_index('pm25' ,EMIS_FILE(1:NEMIS_FILE))
           lf_src(isrc)%iem = iem
           ii = 0 ! index that over only the splits included
           do i=1,emis_nsplit(iem)
              iqrc = sum(emis_nsplit(1:iem-1)) + i
              itot = iqrc2itot(iqrc)
              ix = itot-NSPEC_SHL
              if (index(species(itot)%name,'_new') == 0) cycle !not a "new" : exclude
              ii = ii + 1
              lf_src(isrc)%ix(ii) = ix
              lf_src(isrc)%mw(ii) = species_adv(ix)%molwt
              lf_src(isrc)%Nsplit = ii ! will take value of the last ii
           enddo
        else if(lf_src(isrc)%species=='pm25_tot')then
           !we consider pm25 as a whole, without explicitly including the "chemistry", i.e. conversion new->age
           iem = find_index('pm25' ,EMIS_FILE(1:NEMIS_FILE))
           lf_src(isrc)%iem = iem
           ii = 0 ! index that over only the splits included
           do i=1,emis_nsplit(iem)
              iqrc = sum(emis_nsplit(1:iem-1)) + i
              itot = iqrc2itot(iqrc)
              ix = itot-NSPEC_SHL
              ii = ii + 1
              lf_src(isrc)%ix(ii) = ix
              lf_src(isrc)%mw(ii) = species_adv(ix)%molwt
              lf_src(isrc)%Nsplit = ii ! will take value of the last ii
           enddo
         else
           !defined as single species (NO, NO2, O3..)
           lf_src(isrc)%Nsplit = 1
           ix=find_index(lf_src(isrc)%species ,species(:)%name)
           if(ix<0)then
              ix=find_index(lf_src(isrc)%species ,species(:)%name, any_case=.true.) !NB: index among all species also short lived
              if(me==0 .and. ix>0)then
                 write(*,*)'WARNING: '//trim(lf_src(isrc)%species)//' associated to '//trim(species(ix)%name)
!                 lf_src(isrc)%species=trim(species(ix)%name)
              end if
              if(ix<0)write(*,*)isrc,NSOURCES
           end if
           call CheckStop( ix<1, "Local Fractions did not find corresponding pollutant: "//trim(lf_src(isrc)%species) )
           iem=-1
           lf_src(isrc)%species_ix = ix !NB: index among all species
           lf_src(isrc)%ix(1) = ix - NSPEC_SHL !NB: index among advected species
           lf_src(isrc)%mw(1) = species_adv(lf_src(isrc)%ix(1))%molwt
           lf_src(isrc)%iqrc = itot2iqrc(ix) !negative if not among emitted species
           if(lf_src(isrc)%iqrc>0) iem = iqrc2iem(lf_src(isrc)%iqrc)
           lf_src(isrc)%iem = iem
           if (lf_set%full_chem .and. iem>0) then
              !Do not include emissions which are not NOx, VOC, NH3 or SOx emissions
              if (lf_set%EmisDer_all .or. &
              (EMIS_FILE(iem)/='nox' .and. EMIS_FILE(iem)/='voc' .and. EMIS_FILE(iem)/='nh3' .and. EMIS_FILE(iem)/='sox') .or. &
              (Nfullchem_emis == 2 .and. EMIS_FILE(iem)/='nox' .and. EMIS_FILE(iem)/='voc'))&
              lf_src(isrc)%iem = -1
           end if

           if (.not. lf_set%full_chem) then
              !(cannot be used for full_chem, because several
              !isrc are used for same species: one for each derivative)
              !use simplified methods
              if(trim(species(ix)%name)=='O3')isrc_O3=isrc
              if(trim(species(ix)%name)=='NO')isrc_NO=isrc
              if(trim(species(ix)%name)=='NO2')isrc_NO2=isrc
              if(trim(species(ix)%name)=='SO4')isrc_SO4=isrc
              if(trim(species(ix)%name)=='SO2')isrc_SO2=isrc
              if(trim(species(ix)%name)=='NH4_f')isrc_NH4_f=isrc
              if(trim(species(ix)%name)=='NH3')isrc_NH3=isrc
           end if
        end if
     else
        !species defines as primary emitted
        lf_src(isrc)%iem = iem
        if (lf_src(isrc)%species=="pmco") isrc_pmco = isrc
        ii = 0 ! index that over only the splits included
        do i=1,emis_nsplit(iem)
           iqrc=sum(emis_nsplit(1:iem-1)) + i
           itot=iqrc2itot(iqrc)
           ix=itot-NSPEC_SHL
           if(lf_src(isrc)%species == 'pm25') then
              isrc_pm25 = isrc
              !we include only the parts of pm25 which are not "new"
              if(MasterProc .and. i == 1 .and. index(species(itot)%name,'_new')>0)write(*,*)'including only new pm25 into pm25_new'
              if (index(species(itot)%name,'_new')>0) cycle
           end if
           ii = ii + 1
           lf_src(isrc)%Nsplit = ii
           lf_src(isrc)%ix(ii) = ix
           lf_src(isrc)%mw(ii) = species_adv(ix)%molwt

           if(lf_src(isrc)%species=="voc")then
              isrc_VOC = isrc
              ix=find_index('CH3CO3', species(:)%name)
              call CheckStop( ix<1, "Local Fractions did not find CH3CO3 ")
              ix=find_index('HO2', species(:)%name)
              call CheckStop( ix<1, "Local Fractions did not find HO2 ")
           end if
           if(lf_src(isrc)%species=="nox")then
              ix=find_index("NO2",species_adv(:)%name)
              call CheckStop(ix<0,'Index for NO2 not found')
              lf_src(isrc)%mw(ii)=species_adv(ix)%molwt !"as NO2"
           end if
           if(lf_src(isrc)%species=="sox")then
              ix=find_index("SO2",species_adv(:)%name)
              call CheckStop(ix<0,'Index for SO2 not found')
              lf_src(isrc)%mw(ii)=species_adv(ix)%molwt !"as SO2"
           end if
           if(lf_src(isrc)%species=="nox" .and. (lf_src(isrc)%DryDep .or. lf_src(isrc)%WetDep))then
              ix=find_index("NO3",species_adv(:)%name)
              call CheckStop(ix<0,'Index for NO3 not found')
              ix=find_index("HNO3",species_adv(:)%name)
              call CheckStop(ix<0,'Index for HNO3 not found')
           end if

           if(lf_src(isrc)%species=="nh3")then
              lf_src(isrc)%Nsplit = 0
              ixnh4=find_index("NH4_F",species_adv(:)%name , any_case=.true.)
              ixnh3=find_index("NH3",species_adv(:)%name)
              do ix=1,NSPEC_ADV
                 if(ix/=ixnh4.and.ix/=ixnh3)cycle!not reduced nitrogen
                 if(species_adv(ix)%nitrogens>0)then
                    lf_src(isrc)%Nsplit = lf_src(isrc)%Nsplit + 1
                    lf_src(isrc)%ix(lf_src(isrc)%Nsplit) = ix
                    lf_src(isrc)%mw(lf_src(isrc)%Nsplit) = species_adv(ixnh3)%molwt !use NH3 mw also for NH4
                 end if
              enddo
           end if
        end do

     end if

     if(iem>0)then
        !emitted species
        found=0
        do i=1,iem2Nipoll(iem)
           if (iem2names(iem,i)==lf_src(isrc)%species) then
              found=1
              lf_src(isrc)%poll = iem2ipoll(iem,i)
           end if
        end do
        if (found==0) then
           !add a new pollutant for that emis file
           iem2Nipoll(iem)=iem2Nipoll(iem)+1
           ipoll = ipoll + 1
           call CheckStop(ipoll>Max_lf_spec,"Error: increase Max_lf_spec")
           iem2names(iem, iem2Nipoll(iem)) = trim(lf_src(isrc)%species)
           iem2ipoll(iem, iem2Nipoll(iem)) = ipoll
           lf_src(isrc)%poll = ipoll
           Npoll = ipoll
        end if
        if(lf_src(isrc)%iqrc>0)ipoll2iqrc(lf_src(isrc)%poll)=lf_src(isrc)%iqrc!single species
      else
         !single species not emitted (like O3)
         !check if this has been included already
         found=0
         do ii = 1, isrc-1
            if (lf_src(isrc)%species == lf_src(ii)%species) then
               found = 1
               lf_src(isrc)%poll = lf_src(ii)%poll
            end if
         end do
         if (found == 0) then
            ipoll = ipoll + 1
            lf_src(isrc)%poll = ipoll
            Npoll = ipoll
            call CheckStop(ipoll>Max_lf_spec,"Error: increase Max_lf_spec")
         end if
      end if

      if (MasterProc) then
        if(.not. lf_fullchem .or. lf_src(isrc)%iem_lf == iem_lf_nox) then
!           if (lf_src(isrc)%iem>0) then
!              write(*,*)'lf pollutant : ',lf_src(isrc)%species,' ref index ',lf_src(isrc)%poll,' emitted as ',EMIS_FILE(lf_src(isrc)%iem)
!           else
!              write(*,*)'lf pollutant : ',lf_src(isrc)%species,' ref index ',lf_src(isrc)%poll,' not treated as emitted species'
!           end if
        end if
        if (.not. lf_fullchem) then
        write(*,*)'lf number of species in '//trim(lf_src(isrc)%species)//' group: ',lf_src(isrc)%Nsplit
        write(*,"(A,30(A,F6.2))")'including:',('; '//trim(species_adv(lf_src(isrc)%ix(i))%name)//', mw ',lf_src(isrc)%mw(i),i=1,lf_src(isrc)%Nsplit)
        if (lf_src(isrc)%type/='country') write(*,"(A,I4,A,I4)")'sector:',lf_src(isrc)%sector,',  res:',lf_src(isrc)%res
        !write(*,"(A,30I4)")'ix:',(lf_src(isrc)%ix(i),i=1,lf_src(isrc)%Nsplit)
        end if
     end if
  end do


  if (isrc_SO2>0 .and. (isrc_SO4<0)) then
     if(me==0)write(*,*)'WARNING: SO2 tracking, but no SO4 tracking'
!     stop!may be relaxed in future
  end if

  av_fac=0.0

  Ndrydep_lf=0
  Nwetdep_lf=0
  LF_SRC_TOTSIZE = 0
  is_relative = .false.
  is_country = .false.

  nPODnAOT = 0
  nDryDep = 0
  if (lf_fullchem) then
     !activate drydep for species asked for
     do iout = 1, Max_lf_out
        if (lf_spec_out(iout)%name == "NOTSET") exit
        if (index(lf_spec_out(iout)%name,"POD")>0 .or. index(lf_spec_out(iout)%name,"AOT")>0) then
           !treat separately. out is treated as DryDep, but isrc not marked as DryDep
           lf_spec_out(iout)%DryDep = .true.
           nDryDep = nDryDep + 1
           nPODnAOT = nPODnAOT + 1
           cycle
        end if
        if (.not. lf_spec_out(iout)%DryDep) cycle
        nDryDep = nDryDep + 1
        do ig = 1, 30
           if (lf_spec_out(iout)%species(ig) == "NOTSET" ) exit
           found = 0
           do isrc = 1, Nsources
              if(lf_src(isrc)%species == lf_spec_out(iout)%species(ig)) then
                 lf_src(isrc)%drydep = .true.
                 found = 1
                 !note several sources can be defined for same species (one for each emis for example)
               end if
           end do
           call CheckStop(found == 0," DryDep species not found"//trim(lf_spec_out(iout)%species(ig)))
        end do
     end do
     !define the corresponding indices in the loc_frac_drydep array
     !The indices idep must be defined in the order of isrc
     idep = 0
     do isrc=1,Nsources
        if(.not. lf_src(isrc)%DryDep)cycle !activated above for fullchem
        do iout = 1, Max_lf_out
           if (lf_spec_out(iout)%name == "NOTSET") exit
           if (.not.lf_spec_out(iout)%DryDep) cycle
           do ig = 1, 30
              if (lf_spec_out(iout)%species(ig) == "NOTSET" ) exit
              if(trim(lf_src(isrc)%species) == trim(lf_spec_out(iout)%species(ig))) then
                 if(lf_spec_out(iout)%ix(ig)<0)then
                   lf_spec_out(iout)%ix(ig) = idep !NB: start at zero
                   idep=idep+lf_src(isrc)%Npos * lf_src(isrc)%Nsplit * Nfullchem_emis! each idep in lf_spec_out(iout)%ix(ig) is used for several sources
                   !NB: cannot have several lf_spec_out(iout)%ix(ig) for same isrc
                   exit !only position in drydep for first source is stored
                end if
                end if
             end do
        end do
     end do
     !indices for PODnAOT after indices for other DryDep. isrc not marked as DryDep for POD (but lf_spec_out are marked as DryDep!)
     do iout = 1, Max_lf_out
        if(lf_spec_out(iout)%name == 'NOTSET') exit
        if (index(lf_spec_out(iout)%name,"POD")>0 .or. index(lf_spec_out(iout)%name,"AOT")>0) then
           lf_spec_out(iout)%ix(1) = idep !POD or AOT does not use groups
           idep=idep+Nfullchem_emis*Npos_lf
        end if
     end do

     !Activate wetdep
     Nwetdep_lf = 0
     do iout = 1, Max_lf_out
        if (lf_spec_out(iout)%name == "NOTSET") exit
        if (.not. lf_spec_out(iout)%WetDep) cycle
        Nwetdep_lf = Nwetdep_lf + 1
        do ig = 1, 30
           if (lf_spec_out(iout)%species(ig) == "NOTSET" ) exit
           found = 0
           do isrc = 1, Nsources
              if(lf_src(isrc)%species == lf_spec_out(iout)%species(ig)) then
                 lf_src(isrc)%wetdep = .true.
                 found = 1
                 !note several sources can be defined for same species (one for each emis for example)
               end if
           end do
           call CheckStop(found == 0," WetDep species not found"//trim(lf_spec_out(iout)%species(ig)))
        end do
     end do
     idep = 0
     do isrc=1,Nsources
        if(.not. lf_src(isrc)%WetDep)cycle !activated above for fullchem
        do iout = 1, Max_lf_out
           if (lf_spec_out(iout)%name == "NOTSET") exit
           if (.not.lf_spec_out(iout)%WetDep) cycle
           call CheckStop(lf_spec_out(iout)%DryDep,'Cannot use same output for dry and wet dep')
           do ig = 1, 30
              if (lf_spec_out(iout)%species(ig) == "NOTSET" ) exit
              if(trim(lf_src(isrc)%species) == trim(lf_spec_out(iout)%species(ig))) then
                 if(lf_spec_out(iout)%ix(ig)<0)then
                   lf_spec_out(iout)%ix(ig) = idep !NB: start at zero
                   idep=idep+lf_src(isrc)%Npos * lf_src(isrc)%Nsplit * Nfullchem_emis! each idep in lf_spec_out(iout)%ix(ig) is used for several sources
                   !NB: cannot have several lf_spec_out(iout)%ix(ig) for same isrc
                   exit !only position in drydep for first source is stored
                end if
             end if
             end do
        end do
     end do

  end if

  do isrc = 1, Nsources
     if(lf_src(isrc)%drydep) Ndrydep_lf = Ndrydep_lf + lf_src(isrc)%Npos
     if(lf_src(isrc)%drydep .and. me==0) write(*,*)Ndrydep_lf,trim(lf_src(isrc)%species),isrc,Npos_lf
     if(lf_src(isrc)%wetdep) Nwetdep_lf = Nwetdep_lf + lf_src(isrc)%Npos
     if (lf_src(isrc)%WetDep) then
        if(lf_src(isrc)%species=="nox")then !note that "nox" is not a fullchem species
           if(me==0)write(*,*)'WetDep for nox will include ',(trim(species(WDEP_OXN_GROUP(iix))%name)," ", iix=1,size(WDEP_OXN_GROUP))
           do iix=1,size(WDEP_OXN_GROUP)
              wetdep_lf(WDEP_OXN_GROUP(iix)) = .true.
           end do
        else
           do iix=1,lf_src(isrc)%Nsplit
              ix=lf_src(isrc)%ix(iix)
              wetdep_lf(ix) = .true.
           end do
        end if
     end if
     lf_src(isrc)%start = LF_SRC_TOTSIZE + 1
     lf_src(isrc)%end = LF_SRC_TOTSIZE + lf_src(isrc)%Npos
     LF_SRC_TOTSIZE = LF_SRC_TOTSIZE + lf_src(isrc)%Npos
     if(me==0)then
        write(*,*)isrc,' ',trim(lf_src(isrc)%species),lf_src(isrc)%iem_lf," start ",lf_src(isrc)%start," end ",lf_src(isrc)%end
     end if

     if (lf_src(isrc)%type=='country') is_country = .true.
     if (lf_src(isrc)%type=='relative') is_relative = .true.

  end do
  if (.not. is_country .and. Ncountry_lf>0) then
     if(me==0)write(*,*)'LF: no country output required'
     Ncountry_lf=0
     Ncountry_group_lf=0
     Ncountrysectors_lf=0
     Ncountry_mask_lf=0
     Ncountry_mask_lf_val=0
  end if

  if(me==0 )then
     write(*,*)Ndrydep_lf,' dry deposited sources tracked ',nPODnAOT,' POD or AOT' !Ndrydep_lf does not include PODnAOT
     write(*,*)Nwetdep_lf,' wet deposited sources tracked '
  end if

  allocate(lf(LF_SRC_TOTSIZE,LIMAX,LJMAX,KMAX_MID-lf_Nvert+1:KMAX_MID))
  lf=0.0


  nstratos = 0
  nbc = 0
  do isrc = 1, Nsources
     !IC_STRATOS and IC_BC are special. Make a list of corresponding sources indices
     !IC_INIT is special. init lf to 1.
     if(lf_src(isrc)%type=='country' .and. Ncountry_lf+Ncountry_group_lf>0)then
        n0=lf_src(isrc)%start
        do is=1,Ncountrysectors_lf
           do ic=1,Ncountry_lf+Ncountry_group_lf
              if (country_ix_list(ic)==IC_STRATOS .and. lf_src(isrc)%species == 'O3') then
                 nstratos = nstratos + 1
                 if (nstratos > size(Stratos_ix) )then
                    write(*,*)'Increase size of Stratos_ix to ',Nsources, size(Stratos_ix)
                    stop
                 end if
                 Stratos_ix(nstratos)=n0
              end if
              if (country_ix_list(ic)==IC_BC) then
                 if (lf_set%EmisDer_all) then
                    !we use "CAMS" conventions to pickout the BC to change
                    if (lf_src(isrc)%species == 'O3' .or.&
                         lf_src(isrc)%species == 'CO' .or.&
                         lf_src(isrc)%species == 'NO' .or.&
                         lf_src(isrc)%species == 'NO2' .or.&
                         lf_src(isrc)%species == 'PAN' .or.&
                         lf_src(isrc)%species == 'HNO3' .or.&
                         lf_src(isrc)%species == 'HCHO' .or.&
                         lf_src(isrc)%species == 'SO2' .or.&
                         lf_src(isrc)%species == 'CH4' .or.&
                         lf_src(isrc)%species == 'C5H8' .or.&
                         lf_src(isrc)%species == 'C2H6' .or.&
                         lf_src(isrc)%species == 'SO4' .or.&
                         lf_src(isrc)%species == 'SO2' &
                         !lf_src(isrc)%species == 'NO3_f' .or.&
                         !lf_src(isrc)%species == 'NO3_c' .or.&
                         !lf_src(isrc)%species == 'NH4_f'&
                       ) then
                       nbc = nbc + 1
                       if (nbc > size(BC_ix) )then
                          write(*,*)'Increase size of BC_ix to ',Nsources, size(BC_ix)
                          stop
                       end if
                       BC_ix(nbc)=n0
                    end if
                 else
                    !we track different species in different "emission reductions"
                    if((lf_src(isrc)%iem_lf == iem_lf_voc .and. &
                      (lf_src(isrc)%species == 'NC4H10' .or. lf_src(isrc)%species == 'C2H6'))& !NB: NC4H10 mapped from C4H10 (CMX_BoundaryConditions.txt)
                      .or.&
                      (lf_src(isrc)%iem_lf == iem_lf_sox .and. &
                      (lf_src(isrc)%species == 'SO2' .or. lf_src(isrc)%species == 'SO4'))&
                      .or.&
                      (lf_src(isrc)%iem_lf == iem_lf_nh3 .and. &
                      (lf_src(isrc)%species == 'NH4_f'))&
                      .or.&
                      (lf_src(isrc)%iem_lf == iem_lf_nox .and. &
                      (lf_src(isrc)%species == 'NO3_f' .or. lf_src(isrc)%species == 'NO3_c' .or. &
                        lf_src(isrc)%species == 'HNO3' .or. lf_src(isrc)%species == 'NO2' .or. &
                        lf_src(isrc)%species == 'NO' .or. lf_src(isrc)%species == 'PAN') ) &
                        ) then
                       nbc = nbc + 1
                       if (nbc > size(BC_ix) )then
                          write(*,*)'Increase size of BC_ix to ',Nsources, size(BC_ix)
                          stop
                       end if
                       lf(n0,:,:,:) = 1.0 ! Also initial conditions
                       BC_ix(nbc)=n0
                    end if
                 end if
              end if
              if (country_ix_list(ic)==IC_INIT) lf(n0,:,:,:) = 1.0
              n0=n0+1
           end do
        end do
     end if
  end do

  if(lf_set%HOUR)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_HOUR)=Niou_ix
     iotyp2ix(Niou_ix)=IOU_HOUR
  end if
  if(lf_set%HOUR_INST)then
     Niou_ix = Niou_ix + 1
     iou_ix_inst = Niou_ix !should not be accumulated
     iotyp2ix(IOU_HOUR_inst) = Niou_ix
     iotyp2ix(Niou_ix)=IOU_HOUR_inst
  end if
  if(lf_set%DAY)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_DAY)=Niou_ix
     iotyp2ix(Niou_ix)=IOU_DAY
  end if
  if(lf_set%MONTH)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_MON)=Niou_ix
     iotyp2ix(Niou_ix)=IOU_MON
  end if
  if(lf_set%YEAR)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_YEAR)=Niou_ix
     iotyp2ix(Niou_ix)=IOU_YEAR
  end if

  if (lf_fullchem) then
     !ensure that at least something is outputted
     call CheckStop(lf_spec_out(1)%name == "NOTSET", "At least one output must be asked for")
     !check if PMwater is asked for
     found = 0
     do n = 1, Max_lf_out
       if (lf_spec_out(n)%name == "NOTSET") exit
       if (lf_spec_out(n)%name=="PM_WATER") found = 1
       do i = 1, 30
          if (lf_spec_out(n)%species(i) == "PM_WATER") found = 1
          if (lf_spec_out(n)%species(i) == "NOTSET") exit
       end do
    end do
    if (found == 1) then
       make_PMwater = .true.
       isrc = NSOURCES+1
       lf_src(isrc)%species="PM_WATER"
       lf_src(isrc)%Npos = Npos_lf
       lf_src(isrc)%start = LF_SRC_TOTSIZE + 1
       lf_src(isrc)%end = LF_SRC_TOTSIZE + Npos_lf
       lf_src(isrc)%iem_lf=iem_lf_nox
       lf_src(isrc)%iem_deriv = find_index('nox' ,EMIS_FILE(1:NEMIS_FILE))
       lf_src(isrc)%poll = Npoll + 1
       if (Nfullchem_emis>1) then
          isrc = isrc+1
          lf_src(isrc)%species="PM_WATER"
          lf_src(isrc)%Npos = Npos_lf
          lf_src(isrc)%start = LF_SRC_TOTSIZE + Npos_lf + 1
          lf_src(isrc)%end = LF_SRC_TOTSIZE + 2*Npos_lf
          lf_src(isrc)%iem_lf=iem_lf_voc
          lf_src(isrc)%iem_deriv = find_index('voc' ,EMIS_FILE(1:NEMIS_FILE))
          lf_src(isrc)%poll = Npoll + 1          
       end if
       if (Nfullchem_emis==4) then
          isrc = isrc+1
          lf_src(isrc)%species="PM_WATER"
          lf_src(isrc)%Npos = Npos_lf
          lf_src(isrc)%start = LF_SRC_TOTSIZE + 2*Npos_lf + 1
          lf_src(isrc)%end = LF_SRC_TOTSIZE + 3*Npos_lf
          lf_src(isrc)%iem_lf=iem_lf_nh3
          lf_src(isrc)%iem_deriv = find_index('nh3' ,EMIS_FILE(1:NEMIS_FILE))
          lf_src(isrc)%poll = Npoll + 1
          isrc = isrc+1
          lf_src(isrc)%species="PM_WATER"
          lf_src(isrc)%Npos = Npos_lf
          lf_src(isrc)%start = LF_SRC_TOTSIZE + 3*Npos_lf + 1
          lf_src(isrc)%end = LF_SRC_TOTSIZE + 4*Npos_lf
          lf_src(isrc)%iem_lf=iem_lf_sox
          lf_src(isrc)%iem_deriv = find_index('sox' ,EMIS_FILE(1:NEMIS_FILE))
          lf_src(isrc)%poll = Npoll + 1
       end if
    end if
 end if

  if(isrc_NH4_f>0)then
     allocate(lf_NH4(KMAX_MID-lf_Nvert+1:KMAX_MID))
     allocate(lf_NH3(KMAX_MID-lf_Nvert+1:KMAX_MID))
  end if
  if(make_PMwater)then
     !one extra for water
     allocate(lf_src_acc(LF_SRC_TOTSIZE+Npos_lf*Nfullchem_emis,LIMAX,LJMAX,lf_Nvertout,Niou_ix))
   else
     allocate(lf_src_acc(LF_SRC_TOTSIZE,LIMAX,LJMAX,lf_Nvertout,Niou_ix))
  end if
  lf_src_acc = 0.0
  allocate(lf_src_tot(LIMAX,LJMAX,lf_Nvertout,2*Npoll+1,Niou_ix))
  lf_src_tot = 0.0
  allocate(lf_src_ps(LIMAX,LJMAX,Niou_ix))
  lf_src_ps = 0.0
  allocate(loc_frac_src_1d(LF_SRC_TOTSIZE,0:max(LIMAX,LJMAX)+1))
  loc_frac_src_1d=0.0

  allocate(emis_lf_cntry(LIMAX,LJMAX,NCMAX,Nsectors,NEMIS_File))
  emis_lf_cntry=0.0

  if(NdryDep>0 .or. Ndrydep_lf>0)then
     allocate(loc_frac_drydep(LIMAX,LJMAX,Ndrydep_lf + nPODnAOT*Nfullchem_emis*Npos_lf))
     loc_frac_drydep=0.0
  else
     allocate(loc_frac_drydep(1,1,1))
  end if
  if(Nwetdep_lf>0)then
     allocate(loc_frac_wetdep(LIMAX,LJMAX,Nwetdep_lf))
     loc_frac_wetdep=0.0
  else
     allocate(loc_frac_wetdep(1,1,1))
  end if
  if(lf_fullchem)then
     allocate(xnew_lf(NSPEC_deriv_lf+NSOA+N_lf_derivemisMAX,NSPEC_TOT))
     allocate(x_lf(NSPEC_deriv_lf+NSOA+N_lf_derivemisMAX,NSPEC_TOT))
     allocate(xold_lf(NSPEC_deriv_lf+NSOA+N_lf_derivemisMAX,NSPEC_TOT))
     allocate(Dchem_lf(NSPEC_deriv_lf+NSOA,NSPEC_TOT,KMAX_MID-lf_Nvert+1:KMAX_MID,LIMAX,LJMAX))
     allocate(xderiv(NSPEC_deriv_lf,NSPEC_chem_lf))
     allocate(xderivSOA(NSOA,NSOA))
     allocate(ederiv(N_lf_derivemisMAX,NSPEC_chem_lf))
     allocate(lf0_loc(Npos_lf*Nfullchem_emis,NSPEC_deriv_lf))
     allocate(lf0SOA_loc(Npos_lf*Nfullchem_emis,NSOA))
     allocate(lf_loc(Npos_lf*Nfullchem_emis,NSPEC_chem_lf))
     allocate(lfSOA_loc(Npos_lf*Nfullchem_emis,NSOA))
     allocate(rcemis_lf(N_lf_derivemisMAX,ix_lf_max))
     allocate(totrcemis_lf(Ncountrysectors_lf,Ncountry_lf + Ncountry_group_lf,NEMIS_File))
     totrcemis_lf = 0.0
     allocate(rcemis_lf_primary(NCMAX*NSECTORS*2))
     allocate(rcemis_lf_surf(10,20),emis2spec_surf(10,20)) !up to 10 different sources, 20 species each
     allocate(P_lf(NSPEC_deriv_lf+NSOA+N_lf_derivemisMAX))
     allocate(L_lf(NSPEC_deriv_lf+NSOA+N_lf_derivemisMAX))
!     allocate(rctA_lf(NSPEC_deriv_lf+N_lf_derivemisMAX))
!     allocate(rctAk_lf(NSPEC_deriv_lf+N_lf_derivemisMAX,KMAX_MID-lf_Nvert+1:KMAX_MID))
!     allocate(rctB_lf(NSPEC_deriv_lf+N_lf_derivemisMAX))
!     allocate(rctBk_lf(NSPEC_deriv_lf+N_lf_derivemisMAX,KMAX_MID-lf_Nvert+1:KMAX_MID))
     allocate(fgasso2_lf(NSPEC_deriv_lf+N_lf_derivemisMAX,KMAX_MID-lf_Nvert+1:KMAX_MID))
     fgasso2_lf=1.0
     allocate(AQRCK_lf(NSPEC_deriv_lf+N_lf_derivemisMAX,NAQUEOUS,KMAX_MID-lf_Nvert+1:KMAX_MID))
     AQRCK_lf=0.0
     allocate(xn_shl_lf(NSPEC_deriv_lf+NSOA,NSPEC_SHL,KMAX_MID-lf_Nvert+1:KMAX_MID,LIMAX,LJMAX))
     allocate(lf_PM25_water(Npos_lf*Nfullchem_emis,LIMAX,LJMAX))
     lf_PM25_water = 0.0
     xnew_lf = 0.0
     x_lf = 0.0
     xold_lf = 0.0
     Dchem_lf = 0.0
!     rctA_lf = 0.0
!     rctAk_lf = 0.0
!     rctB_lf = 0.0
!     rctBk_lf = 0.0
     xn_shl_lf = 0.0
     rcemis_lf = 0.0 !NB: important
     rcemis_lf_primary = 0.0 !NB: important
     rcemis_lf_surf = 0.0 !NB: important
     emis2nspec_surf = 0
     emis2iic_surf = 0
     allocate(lf_NO3(KMAX_MID-lf_Nvert+1:KMAX_MID))
     allocate(lf_HNO3(KMAX_MID-lf_Nvert+1:KMAX_MID))
  else
      allocate(rcemis_lf(Nsources*NCMAX*NSECTORS,1))
      allocate(rcemis_lf_primary(Nsources))
      rcemis_lf = 0.0 !NB: important
      rcemis_lf_primary = 0.0
  end if
  allocate(nic(LIMAX,LJMAX))
  nic = 0
  allocate(ic2iland(LIMAX,LJMAX,NCMAX))

  if (lf_set%MDA8) then! maximum daily eight-hour mean concentration
    !we need one array to save the last 8 hourly concentrations, and one to make
    !the average over the last hour
    !index 0 is for the pure O3, and the other indices are for the sensibilities
    allocate(D8M(0:Npos_lf*Nfullchem_emis,LIMAX,LJMAX,8)) ! running last 8 hour values
    allocate(D8Max(0:Npos_lf*Nfullchem_emis,LIMAX,LJMAX)) ! max value of the 8 hour mean since 00:00
    allocate(D8Max_av(LIMAX,LJMAX,0:Npos_lf*Nfullchem_emis,Niou_ix)) ! average over time of D8Max
    allocate(D8Max_av_ppb(LIMAX,LJMAX,0:Npos_lf*Nfullchem_emis,Niou_ix)) ! average over time of D8Max in ppb units
    allocate(D8Max_6m_av(LIMAX,LJMAX,0:Npos_lf*Nfullchem_emis)) ! average over time of D8Max
    allocate(D8Max_6m_av_ppb(LIMAX,LJMAX,0:Npos_lf*Nfullchem_emis)) ! average over time of D8Max in ppb units
    allocate(hourM(0:Npos_lf*Nfullchem_emis,LIMAX,LJMAX)) ! hour Mean
    allocate(D1Max(0:Npos_lf*Nfullchem_emis,LIMAX,LJMAX)) ! Daily max hourM
    allocate(D1Max_av(LIMAX,LJMAX,0:Npos_lf*Nfullchem_emis,Niou_ix)) ! average over time of Daily max hourM

    allocate(AOT40(LIMAX,LJMAX,0:Npos_lf*Nfullchem_emis,Niou_ix)) ! 
    allocate(SOMO35(LIMAX,LJMAX,0:Npos_lf*Nfullchem_emis,Niou_ix)) ! accumulated daily max value of the 8 hour mean since 00:00 over 35 ppb
    allocate(D8M_ppb(0:Npos_lf*Nfullchem_emis,LIMAX,LJMAX,8)) ! running last 8 hour values
    allocate(D8Max_ppb(0:Npos_lf*Nfullchem_emis,LIMAX,LJMAX)) ! max value of the 8 hour mean since 00:00
    allocate(hourM_ppb(0:Npos_lf*Nfullchem_emis,LIMAX,LJMAX)) ! hour Mean

    hourM = 0.0
    D8Max = 0.0 !init with low value
    D1Max = 0.0 !init with low value
    D8M = 0.0 !init with low value
    D8Max_av = 0.0 !init necessary
    D8Max_6m_av = 0.0 !init necessary
    D1Max_av = 0.0 !init necessary
    D8Max_av_ppb = 0.0 !init necessary
    D8Max_6m_av_ppb = 0.0 !init necessary
    SOMO35 = 0.0 !init necessary
    hourM_ppb = 0.0
    D8Max_ppb = 0.0 !init with low value
    D8M_ppb = 0.0 !init with low value

    !output must include O3
    found = 0
    do n = 1, Max_lf_out
       if (lf_spec_out(n)%name == "NOTSET") exit
       if (lf_spec_out(n)%name=="O3") found = 1
       do i = 1, 30
          if (lf_spec_out(n)%species(i) == "O3") found = 1
          if (lf_spec_out(n)%species(i) == "NOTSET") exit
       end do
    end do
    if (found == 0) lf_spec_out(n)%name="O3"
 end if



  if (lf_set%restart) then
     !initialize lf with values save on disk
     filename=date2string(lf_set%filename_read,current_date,mode='YMDH')
     call lf_read(filename, needed=.false.)
  end if

!  call Add_2timing(NTIMING-10,tim_after,tim_before,"lf: init") negligible
  if(DEBUGall .and. me==0)write(*,*)'end init'

end subroutine lf_init

subroutine lf_out(iotyp)
  integer, intent(in) :: iotyp
  character(len=200) ::filename, varname
  real :: xtot,scale,invtot,t1,t2
  integer ::i,j,k,kk,n,n1,n1der,dx,dy,ix,iix,isec,iisec,isec_poll, ideriv
  integer ::ipoll,ipoll_cfac,isec_poll1,isrc,iou_ix,iter,iddep,iwdep
  integer ::ndim,kmax,CDFtype,dimSizes(10),chunksizes(10)
  integer ::ndim_tot,dimSizes_tot(10),chunksizes_tot(10)
  character (len=20) ::dimNames(10),dimNames_tot(10)
  type(Deriv) :: def1 ! definition of fields for local fraction
  type(Deriv) :: def2 ! definition of fields for totals
  type(Deriv) :: def3 ! definition of dry and wet dep fields
  type(Deriv) :: defps ! definition of surface pressure fields
  logical ::overwrite, create_var_only
  logical,save :: first_call(10)=.true.
  real,allocatable ::tmp_out(:,:,:,:)!allocate since it may be heavy for the stack TEMPORARY
  real,allocatable ::tmp_out_cntry(:,:,:)!allocate since it may be heavy for the stack TEMPORARY
  real,allocatable ::tmp_out_base(:,:)! base concentrations
  type(date) :: onesecond = date(0,0,0,0,1)
  character(len=TXTLEN_NAME),save :: oldhourlyname = 'NOTSET'
  character(len=TXTLEN_NAME),save :: oldhourlyInstname = 'NOTSET'
  character(len=TXTLEN_NAME),save :: oldmonthlyname = ''
  character(len=TXTLEN_NAME) :: suffix, specname, sourcename, secname, redname, fullname
  real :: fracsum(LIMAX,LJMAX),fac,invfac
  logical :: pollwritten(2*Max_lf_spec+1),is_surf
  integer :: ncFileID, iout, ig, found, iem_lf, iem, idep, nend, is
  character (len=TXTLEN_NAME) ::countryname(Max_lf_Country_list)
  call Code_timer(tim_before)
  if(DEBUGall .and. me==0)write(*,*)'start out'

  select case (iotyp)
  case (IOU_HOUR_INST)
    if (.not. lf_set%HOUR_INST) return
    fileName = trim(runlabel1)//'_LF_hourInst'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
    if(oldhourlyInstname/=fileName)then
      first_call(iotyp) = .true.
      oldhourlyInstname = fileName
    end if
  case (IOU_HOUR)
    if (.not. lf_set%HOUR) return
    fileName = trim(runlabel1)//'_LF_hour'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
    if(oldhourlyname/=fileName)then
      first_call(iotyp) = .true.
      oldhourlyname = fileName
    end if
  case (IOU_DAY)
    if (.not. lf_set%DAY) return
    fileName=trim(runlabel1)//'_LF_day.nc'
  case (IOU_MON)
    if (.not. lf_set%MONTH) return
    if(lf_set%MONTH_ENDING == "NOTSET")then
      fileName=trim(runlabel1)//'_LF_month.nc'
    else
      fileName=trim(runlabel1)//'_LF_month'//date2string(trim(lf_set%MONTH_ENDING),current_date,-1.0)
      if(oldmonthlyname/=fileName)then
        first_call(iotyp) = .true.
        oldmonthlyname = fileName
      end if
    end if
  case (IOU_YEAR)
    if (.not. lf_set%YEAR) return
    if (lf_set%saveatend) then
      if(me==0)write(*,*)'saving ALL'
      call lf_saveall(date2string(lf_set%filename_writeatend,current_date,mode='YMDH'))
    end if
    fileName=trim(runlabel1)//'_LF_full.nc'
  case (IOU_INST)    
    return ! there is no pure inst outout
  case default
    call CheckStop("LF: unknown iotyp")
  end select

  ncFileID=LF_ncFileID_iou(iotyp)

  ndim=5
  ndim_tot=3
  kmax=lf_Nvertout
  scale=1.0
  CDFtype=Real4
  dimSizes=1

  dimSizes(1)=Ndiv_rel
  dimNames(1)='x_dist'
  dimSizes(2)=Ndiv_rel
  dimNames(2)='y_dist'
  
  dimSizes(3)=min(GIMAX,lf_set%DOMAIN(2)-lf_set%DOMAIN(1)+1)
  dimSizes(4)=min(GJMAX,lf_set%DOMAIN(4)-lf_set%DOMAIN(3)+1)

  dimSizes_tot(1)=min(GIMAX,lf_set%DOMAIN(2)-lf_set%DOMAIN(1)+1)
  dimSizes_tot(2)=min(GJMAX,lf_set%DOMAIN(4)-lf_set%DOMAIN(3)+1)

  select case(projection)
  case('Stereographic')
     dimNames(3)='i'
     dimNames(4)='j'
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'
  case('lon lat')
     dimNames(3)='lon'
     dimNames(4)='lat'
     dimNames_tot(1)='lon'
     dimNames_tot(2)='lat'
  case('Rotated_Spherical')
     dimNames(3)='i'
     dimNames(4)='j'
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'
  case('lambert')
     dimNames(3)='i'
     dimNames(4)='j'
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'
  case default
     dimNames(3)='i'
     dimNames(4)='j'
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'
  end select

  dimSizes(5)=kmax
  dimNames(5)='lev'
  dimSizes_tot(3)=kmax
  dimNames_tot(3)='lev'
  def1%class='LF' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=1.0      !not used
  def1%name='notset'
  def1%unit=''
  def2=def1
  def2%unit='ug/m3'
  def3=def1
  def3%unit='mg/m2'
  defps%class='PSURF' !written
  defps%avg=.false.      !not used
  defps%index=0          !not used
  defps%scale=0.01
  defps%name='PS'
  defps%unit='Pa'
  chunksizes=1
  chunksizes(1)=dimSizes(1)
  chunksizes(2)=dimSizes(2)
  chunksizes(3)=min(34, dimSizes(3))! optimal for read and write?
  chunksizes(4)=min(42, dimSizes(4)) ! optimal for read and write?
  chunksizes(5)=1
  chunksizes_tot=1
  chunksizes_tot(1)=min(MAXLIMAX, dimSizes_tot(1))
  chunksizes_tot(2)=min(MAXLJMAX, dimSizes_tot(2))
  chunksizes_tot(3)=dimSizes_tot(3)

  if(lf_fullchem) then
     allocate(tmp_out(Ndiv_rel*Ndiv_rel,LIMAX,LJMAX,Ncountrysectors_lf))
     allocate(tmp_out_cntry(LIMAX,LJMAX,max(Nfullchem_emis*Npos_lf,(Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf+1)))
  else
     allocate(tmp_out(Ndiv_rel*Ndiv_rel,LIMAX,LJMAX,KMAX))
     allocate(tmp_out_cntry(LIMAX,LJMAX,max((Ncountry_lf + Ncountry_group_lf)*Ncountrysectors_lf,Ndiv_rel*Ndiv_rel)))
  end if

  allocate(tmp_out_base(LIMAX,LJMAX))
  tmp_out_base = 0.0

  iou_ix = iotyp2ix(iotyp)
  if (iou_ix == iou_ix_inst) av_fac(iotyp) = 1

  !first loop only create all variables before writing into them (faster for NetCDF)
  do iter=1,2
     if(iter==1 .and. .not. first_call(iotyp))cycle
     overwrite=(iter==1)      ! overwrite if file already exists
     create_var_only=(iter==1)! only create all variables before writing into them

     pollwritten = .false.
     iddep = 0
     iwdep = 0
     ! we always output surface pressure, unless USES%MEAN_MASK_OUTPUT
     scale = defps%scale/av_fac(iotyp) 
     if(.not. USES%MEAN_MASK_OUTPUT) &
       call Out_netCDF(iotyp,defps,2,kmax,lf_src_ps(1,1,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
          fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
     overwrite=.false.

     if (.not. lf_fullchem) then
        do isrc = 1, Nsources
           if (trim(lf_src(isrc)%species) == 'pm25_new') cycle !we do not output pm25_new (it is included in pm25)
           isec=lf_src(isrc)%sector
           ipoll=lf_src(isrc)%poll
           if(lf_src(isrc)%type=='country')then
              is_surf = .true.
              ipoll_cfac = ipoll + Npoll
           else
              is_surf = .false.
              ipoll_cfac = ipoll
           end if
           if(.not. pollwritten(ipoll_cfac))then !one pollutant may be used for several sources
              def2%name=trim(lf_src(isrc)%species)
              if(iter==1 .and. me==0.and.  first_call(iotyp))write(*,*)' poll '//trim(lf_src(isrc)%species),ipoll_cfac
              scale=1.0/av_fac(iotyp)
              if(is_surf)def2%name='SURF_'//trim(def2%name)
              call Out_netCDF(iotyp,def2,ndim_tot,kmax,lf_src_tot(1,1,1,ipoll_cfac,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                   fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
              pollwritten(ipoll_cfac) = .true.
           end if

           if(iter==2)then
              fracsum=0.0
              tmp_out=0.0
              if(lf_src(isrc)%type == 'country')tmp_out_cntry=0.0
              do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
                 kk = k
                 if (lf_Nvertout < KMAX_MID) kk = KMAX_MID-k+1 !index which is 1 for lowest level and lf_Nvertout at highest level outputted
                 do j=1,ljmax
                    do i=1,limax
                       invtot=1.0/(lf_src_tot(i,j,kk,ipoll_cfac,iou_ix)+1.E-20)
                       n1=0
                       if(lf_src(isrc)%type == 'country')then
                          call CheckStop(lf_Nvertout > 1, "Only one vertical level in output implemented for country type")
                          invfac=1.0/av_fac(iotyp) !could also output fractions?
                          do n=lf_src(isrc)%start, lf_src(isrc)%end
                             n1=n1+1
                             tmp_out_cntry(i,j,n1) = lf_src_acc(n,i,j,kk,iou_ix)*invfac 
                             fracsum(i,j)=fracsum(i,j)+lf_src_acc(n,i,j,kk,iou_ix)*invtot ! sum over all n and k and divided by tot
                          enddo
                       else
                          do n=lf_src(isrc)%start, lf_src(isrc)%end
                             n1=n1+1
                             tmp_out(n1,i,j,kk) = lf_src_acc(n,i,j,kk,iou_ix)*invtot 
                             fracsum(i,j)=fracsum(i,j)+lf_src_acc(n,i,j,kk,iou_ix)*invtot ! sum over all n and k
                          enddo
                       end if
                    enddo
                 enddo
              enddo
           end if

           if(lf_src(isrc)%type == 'country')then
              n1=0
              do j=1,Ncountrysectors_lf
                 do i=1,Ncountry_lf+Ncountry_group_lf
                    n1=n1+1
                    !single cell source
                    isec=lf_country%sector_list(j)
                    if(lf_country%sector_list(j)>=0)isec=lf_country%sector_list(j)
                    if(i<=Ncountry_mask_lf)then
                       if (iic2ilf_countrymask(i) > 0) then
                          write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(EmisMaskIndex2Name(iic2ilf_countrymask(i)))
                          if(isec==0) write(def2%name,"(A)")trim(lf_src(isrc)%species)//'_'//trim(EmisMaskIndex2Name(iic2ilf_countrymask(i)))
                       else
                          write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(mask2name(country_mask_val(i)))
                          if(isec==0) write(def2%name,"(A)")trim(lf_src(isrc)%species)//'_'//trim(mask2name(country_mask_val(i)))
                       end if
                    else if(i<=Ncountry_lf)then
                       write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country%list(i-Ncountry_mask_lf))
                       if(isec==0) write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country%list(i-Ncountry_mask_lf))
                    else
                       !country group
                       write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country%group(i-Ncountry_lf)%name)
                       if(isec==0) write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country%group(i-Ncountry_lf)%name)
                    end if
                    if(me==0 .and. iter==2 .and. (iotyp==IOU_MON .or. iotyp==IOU_YEAR))write(*,*)'writing '//trim(def2%name)
                    def2%unit='ug/m3'
                    scale=1.0
                    call Out_netCDF(iotyp,def2,ndim_tot,1,tmp_out_cntry(1,1,n1),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                        fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                    
                    if(lf_src(isrc)%drydep)then
                       write(def3%name,"(A)")'DDEP_'//trim(def2%name)
                       def3%unit='mg/m2'
                       if(isrc==isrc_SO4 .or. isrc==isrc_SO2 .or. lf_src(isrc)%species=="sox")def3%unit='mgS/m2'
                       if(isrc==isrc_NH3 .or. isrc==isrc_NH4_f .or. lf_src(isrc)%species=="nh3")def3%unit='mgN/m2'

                       iddep=iddep+1
                       call Out_netCDF(iotyp,def3,ndim_tot,1,loc_frac_drydep(1,1,iddep),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                            fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                    end if
                    if(lf_src(isrc)%wetdep)then
                       write(def3%name,"(A)")'WDEP_'//trim(def2%name)
                       def3%unit='mg/m2'
                       if(isrc==isrc_SO4 .or. isrc==isrc_SO2 .or. lf_src(isrc)%species=="sox")def3%unit='mgS/m2'
                       if(isrc==isrc_NH3 .or. isrc==isrc_NH4_f .or. lf_src(isrc)%species=="nh3")def3%unit='mgN/m2'
                       if(isrc==isrc_NO2 .or. isrc==isrc_NO .or. lf_src(isrc)%species=="nox")def3%unit='mgN/m2'

                       iwdep=iwdep+1
                       call Out_netCDF(iotyp,def3,ndim_tot,1,loc_frac_wetdep(1,1,iwdep),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                            fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                    end if
                 enddo
              enddo
           else
              !relative
              def1%unit = ''!default
              def1%class = ''!default
              if(lf_src(isrc)%name=='NOTSET')then
                 write(def1%name,"(A,I2.2,A,I0,A,I0)")trim(lf_src(isrc)%species)//'_sec',isec,'_fraction_',lf_src(isrc)%res,'x',lf_src(isrc)%res
                 if(isec==0) write(def1%name,"(A,I0,A,I0)")trim(lf_src(isrc)%species)//'_fraction_',lf_src(isrc)%res,'x',lf_src(isrc)%res
              else
                 def1%name=trim(lf_src(isrc)%name)
              end if

              if(lf_src(isrc)%type == 'relative')write(def1%unit,fmt='(A)')'fraction'
              if(lf_src(isrc)%type == 'relative')write(def1%class,fmt='(A,I0,A,I0)')'source_size_',lf_src(isrc)%res,'x',lf_src(isrc)%res
              if(lf_src(isrc)%nhour>0)write(def1%name,fmt='(A,I0)')trim(def1%name)//'_t',lf_src(isrc)%time_ix
              scale=1.0
              call Out_netCDF(iotyp,def1,ndim,kmax,tmp_out,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=lf_set%DOMAIN,&
                   fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes,ncFileID_given=ncFileID)
           end if


           if(lf_src(isrc)%make_fracsum)then
              if (lf_src(isrc)%name/='NOTSET')then
                 def1%name=trim(lf_src(isrc)%name)//'_fracsum'
              else if (lf_src(isrc)%type == 'country') then
                 write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_fracsum'
                 if(isec==0) write(def1%name,"(A)")trim(lf_src(isrc)%species)//'_fracsum'
              else
                 write(def1%name,"(A,I2.2,A,I0,A,I0)")trim(lf_src(isrc)%species)//'_sec',isec,'_fracsum_',lf_src(isrc)%res,'x',lf_src(isrc)%res
                 if(isec==0) write(def1%name,"(A,I0,A,I0)")trim(lf_src(isrc)%species)//'_fracsum_',lf_src(isrc)%res,'x',lf_src(isrc)%res
              end if

              call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                   fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

           end if
        enddo


     else !FULLCHEM

        is_surf = .true.
        countryname(Nfullchem_emis*Npos_lf+1) = 'Base'
        do iout = 1, Max_lf_out
           if (lf_spec_out(iout)%name == "NOTSET") exit
           fracsum=0.0
           tmp_out_base = 0.0
           do ideriv = 1, Nfullchem_emis
              tmp_out_cntry = 0.0
              tmp_out = 0.0
              do ig = 1, 30
                 found = 1
                 if (lf_spec_out(iout)%species(ig) == "NOTSET" .and. ig>1) exit
                 if ((lf_spec_out(iout)%species(ig) == "pm25" .or. lf_spec_out(iout)%species(ig) == "pmco") .and. ideriv>1) cycle
                 if (lf_spec_out(iout)%species(ig) == "NOTSET" .and. ig==1) then
                    isrc=find_index(trim(lf_spec_out(iout)%name) ,lf_src(:)%species, nth = ideriv)
                 else
                    isrc=find_index(trim(lf_spec_out(iout)%species(ig)) ,lf_src(:)%species, nth = ideriv)
                 end if
                 if((index(lf_spec_out(iout)%name,"POD")>0 .or. index(lf_spec_out(iout)%name,"AOT")>0 ) .and. ig==1) isrc = isrc_O3
                 if(isrc<1)then
                    found = 0
                    if (me==0 .and. ig>1 .and. ideriv==1) write(*,*)'did not find species '//trim(lf_spec_out(iout)%species(ig))
                    if (me==0 .and. ideriv==1) write(*,*)'WARNING: will not write out '//trim(lf_spec_out(iout)%name)
                    exit
                 end if
                 ipoll_cfac = lf_src(isrc)%poll + Npoll
                 if (iter==2) then
                    if (lf_spec_out(iout)%DryDep) then
                       n1=0
                       idep = lf_spec_out(iout)%ix(ig) + Npos_lf*(ideriv-1) !where to find the dep values
                       fac = 1.0/lf_src(isrc)%mw(1) !to make output unit in S or N
                       if(index(lf_spec_out(iout)%name,"POD")>0 .or. index(lf_spec_out(iout)%name,"AOT")>0)then
                          fac = 1.0
                       end if
                       do n=1,Npos_lf
                          n1 = n1 + 1
                          do j=1,ljmax
                             do i=1,limax
                                tmp_out_cntry(i,j,n1) = tmp_out_cntry(i,j,n1) + loc_frac_drydep(i,j,idep+n1)*lf_spec_out(iout)%species_fac(1)*fac
                             end do
                          end do
                       end do

                    else if (lf_spec_out(iout)%WetDep) then
                       n1=0
                       idep = lf_spec_out(iout)%ix(ig) + Npos_lf*(ideriv-1) !where to find the dep values
                       fac = 1.0/lf_src(isrc)%mw(1) !to make output unit in S or N
                       do n=1,Npos_lf
                          n1 = n1 + 1
                          do j=1,ljmax
                             do i=1,limax
                                tmp_out_cntry(i,j,n1) = tmp_out_cntry(i,j,n1) + loc_frac_wetdep(i,j,idep+n1)*lf_spec_out(iout)%species_fac(1)*fac
                             end do
                          end do
                       end do

                    else
                       invfac=lf_spec_out(iout)%species_fac(ig) / av_fac(iotyp)
                       do j=1,ljmax
                          do i=1,limax
                             k = KMAX_MID
                             kk = 1
                             invtot=1.0/(lf_src_tot(i,j,kk,ipoll_cfac,iou_ix)+1.E-20)
                             if(ideriv == 1)tmp_out_base(i,j) = tmp_out_base(i,j) + lf_src_tot(i,j,kk,ipoll_cfac,iou_ix) * invfac
                             n1=0
                             !we must put the "country" and "relative" into different buffers
                             do n=lf_src(isrc)%start, lf_src(isrc)%start + (Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf-1
                                n1 = n1 + 1
                                tmp_out_cntry(i,j,n1) = tmp_out_cntry(i,j,n1) + lf_src_acc(n,i,j,kk,iou_ix)*invfac
                             end do

                             if (lf_spec_out(iout)%relative) then
                                !continue with n for relatives
                                do is=1,Ncountrysectors_lf !each sector makes a new output
                                    do n1=1,Ndiv_rel*Ndiv_rel
                                       tmp_out(n1,i,j,is) = tmp_out(n1,i,j,is) + lf_src_acc(n,i,j,kk,iou_ix)*invfac 
                                       n=n+1
!                                    fracsum(i,j)=fracsum(i,j)+lf_src_acc(n,i,j,kk,iou_ix)*invtot ! sum over all n and k
                                    end do
                                 end do
                              end if
                           end do
                       end do
                    end if
                 end if
               end do !loop over ig
              if (found==0) cycle
              specname = trim(lf_spec_out(iout)%name)
              if (USES%MEAN_MASK_OUTPUT) then
                if (iter==2 .and. ideriv == Nfullchem_emis) then
                  !only write integral over cities
                  n1= Npos_lf + 1
                  do j=1,ljmax
                    do i=1,limax
                      tmp_out_cntry(i,j,n1) = tmp_out_base(i,j) !store together with sources
                    end do
                  end do
                end if
              else
                if (ideriv == 1 .and. .not. lf_spec_out(iout)%DryDep .and. .not. lf_spec_out(iout)%WetDep) then
                  !first write "base" concentrations, not country contributions
                  scale = 1.0
                  def2%name=trim(specname)
                  !if (is_surf .and. .not. lf_fullchem) def2%name='SURF_ug_'//trim(specname)
                  !if (index(lf_spec_out(iout)%name,"ASO")>0) def2%name='SURF_ug_PM_'//trim(specname)
                  if(lf_spec_out(iout)%name=="SURF_ug_PM25_rh50_ppm")def2%name ="SURF_ug_PPM25"
                  if(lf_spec_out(iout)%name=="SURF_ug_PM10_rh50_ppm")def2%name ="SURF_ug_PPM10"
                  if(iter==2 .and. me==0.and.  first_call(iotyp))write(*,*)'writing Base '//trim(def2%name)
                  call Out_netCDF(iotyp,def2,ndim_tot,kmax,tmp_out_base,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                      fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                  pollwritten(ipoll_cfac) = .true.
                end if
              end if
              !now write out sensibilities for each country and sectors
              !(uses last defined isrc, assumes that value is same for all species in group)
              n1=0
              do j=1,Ncountrysectors_lf
                 do i=1,Ncountry_lf+Ncountry_group_lf
                    n1=n1+1

                    redname='' !add emission species to put in netcdf variable name
                    if(lf_set%EmisDer_all) then
                       redname='_PSAVN'
                    else if (lf_spec_out(iout)%name == "pm25" .or. lf_spec_out(iout)%name == "pmco")then
                       redname ='_P'
                    else if (lf_src(isrc)%iem_deriv > 0) then
                       redname='_'//trim(EMIS_FILE(lf_src(isrc)%iem_deriv))
                    end if
                    
                    isec=lf_country%sector_list(j)
                    secname = ''
                    if (isec<=NSECTORS) then
                       if(isec/=0)write(secname,"(A,I2.2)")'_sec',isec
                    else
                       secname='_'//trim(lf_sector_groups(isec-NSECTORS)%name)
                    end if
                    if(i<=Ncountry_mask_lf)then
                      !mask defined region
                      if (iic2ilf_countrymask(i) > 0) then
                        sourcename = '_'//trim(EmisMaskIndex2Name(iic2ilf_countrymask(i)))
                        if(iter==1 .and. iout==1)countryname((ideriv-1)*Npos_lf+n1)=trim(EmisMaskIndex2Name(iic2ilf_countrymask(i)))//trim(secname)
                      else
                        sourcename = '_'//trim(mask2name(country_mask_val(i)))
                        if(iter==1 .and. iout==1)countryname((ideriv-1)*Npos_lf+n1)=trim(mask2name(country_mask_val(i)))//trim(secname)
                      end if
                    else if(i<=Ncountry_lf)then
                      !regular country
                      sourcename = '_'//trim(lf_country%list(i-Ncountry_mask_lf))
                       ix = find_index(trim(lf_country%list(i-Ncountry_mask_lf)) ,Country(:)%code, first_only=.true.)
                       if (ix < 0) then
                         if(iter==1 .and. iout==1)countryname((ideriv-1)*Npos_lf+n1)=trim(lf_country%list(i-Ncountry_mask_lf))
                       else
                         if(iter==1 .and. iout==1)countryname((ideriv-1)*Npos_lf+n1)=trim(Country(ix)%name)//trim(secname)
                       end if
                     else
                       !country group
                       sourcename = '_'//trim(lf_country%group(i-Ncountry_lf)%name)
                       if(iter==1 .and. iout==1)countryname((ideriv-1)*Npos_lf+n1)=trim(lf_country%group(i-Ncountry_lf)%name)//trim(secname)
                    end if
                    if(iotyp==IOU_HOUR_INST .and. iter==1  .and. iout==1)then
                       if(.not.lf_set%EmisDer_all) then
                         countryname((ideriv-1)*Npos_lf+n1)=trim(countryname((ideriv-1)*Npos_lf+n1))//'_'//trim(EMIS_FILE(lf_src(isrc)%iem_deriv))
                       end if
                    end if
                    if( country_ix_list(i)==IC_STRATOS .or. country_ix_list(i)==IC_INIT .or. &
                         country_ix_list(i)==IC_BVOC .or. &
                         country_ix_list(i)==IC_DMS.or. &
                         country_ix_list(i)==IC_NAT )then
                       !do not add "_nox" suffix and do not output voc,nh3,sox "derivatives"
                       if(EMIS_FILE(lf_src(isrc)%iem_deriv) == 'voc' .or.&
                            EMIS_FILE(lf_src(isrc)%iem_deriv) == 'nh3'.or.&
                            EMIS_FILE(lf_src(isrc)%iem_deriv) == 'sox') cycle
                       redname = ''
                    end if
                    scale=1.0
                    def2%unit='ug/m3'
                    if (lf_spec_out(iout)%DryDep .or. lf_spec_out(iout)%WetDep) then
                       def2%unit='mg/m2'
                       if(index(lf_spec_out(iout)%name,"SO")>0)then
                          def2%unit='mgS/m2'
                          scale = 32.0
                       end if
                       if(index(lf_spec_out(iout)%name,"NO")>0 .or. index(lf_spec_out(iout)%name,"OXN")>0)then
                          def2%unit='mgN/m2'
                          scale = 14.0
                       end if
                       if(index(lf_spec_out(iout)%name,"NH")>0 .or. index(lf_spec_out(iout)%name,"RDN")>0)then
                          def2%unit='mgN/m2'
                          scale = 14.0
                       end if

                       if(index(lf_spec_out(iout)%name,"POD")>0 .or. index(lf_spec_out(iout)%name,"AOT")>0)then
                          if (index(lf_spec_out(iout)%name,"POD")>0 ) then
                             def2%unit='mmole/m2'
                             scale = dt_advec*1e-6
                          else if (index(lf_spec_out(iout)%name,"AOT")>0 ) then
                             def2%unit='ppb h'
                             scale = dt_advec/3600.
                          end if
                          !for PODnAOT we do not have a specific isrc, so we have to define redname explicitely
                          if(ideriv==1)redname='_nox'
                          if(ideriv==2)redname='_voc'
                          if(ideriv==3)redname='_nh3'
                          if(ideriv==4)redname='_sox'
                        end if
                     end if

                    if (USES%MEAN_MASK_OUTPUT) then ! "Compressed" CityMasks output
                       if(iter==1 .or. n1/=1)cycle
                       fullname = 'SURF_ug_'//trim(specname)!//trim(secname)//trim(sourcename)//trim(redname)
                       if(ideriv/=Nfullchem_emis) then
                          call MeanMasksOut(iotyp, tmp_out_cntry, fullname, countryname, &
                              Npos_lf, Nfullchem_emis*Npos_lf + 1, (ideriv-1)*Npos_lf+1)
                       else
                          if(DEBUG .and. me==0) &
                             write(*,*)Nfullchem_emis*Npos_lf + 1,'OUT BASE ',trim(fullname)
                          ! last one write Base
                          call MeanMasksOut(iotyp, tmp_out_cntry, fullname, countryname, &
                              Npos_lf+1, Nfullchem_emis*Npos_lf + 1, (ideriv-1)*Npos_lf+1)
                       end if
                    else
                       def2%name =  trim(specname)//trim(secname)//trim(sourcename)//trim(redname)
                       if(lf_spec_out(iout)%name=="SURF_ug_PM25_rh50_ppm" .and. ideriv>1) cycle
                       if(lf_spec_out(iout)%name=="SURF_ug_PM10_rh50_ppm" .and. ideriv>1) cycle
                       if(lf_spec_out(iout)%name=="SURF_ug_PM25_rh50_ppm")def2%name ="SURF_ug_PM25_rh50"//trim(secname)//trim(sourcename)//'_ppm'
                       if(lf_spec_out(iout)%name=="SURF_ug_PM10_rh50_ppm")def2%name ="SURF_ug_PM10_rh50"//trim(secname)//trim(sourcename)//'_ppm'
                       if(me==0 .and. iter==1 .and. (iotyp==IOU_MON .or. iotyp==IOU_YEAR))write(*,*)'writing '//trim(def2%name)             
                   
                       call Out_netCDF(iotyp,def2,ndim_tot,1,tmp_out_cntry(1,1,n1),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                            fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                       if(lf_set%MDA8 .and. (lf_spec_out(iout)%name == 'O3'.or.lf_spec_out(iout)%name == 'SURF_ug_O3').and. .not. lf_spec_out(iout)%DryDep.and. .not. lf_spec_out(iout)%WetDep)then ! NB: assumes O3 is asked for!
                          write(def2%name,"(A)")"AvgMDA1"//trim(secname)//trim(sourcename)//trim(redname)
                          def2%unit='ug/m3'
                          n1der = (lf_src(isrc)%iem_lf-1)*Npos_lf+n1
                          call Out_netCDF(iotyp,def2,ndim_tot,1,D1Max_av(1,1,n1der,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                               fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                          if (iotyp2ix(iou_ix)==IOU_YEAR) then
                             !6 months, peak season only 
                             write(def2%name,"(A)")"AvgMDA8_6month"//trim(secname)//trim(sourcename)//trim(redname)
                             def2%unit='ug/m3'
                             n1der = (lf_src(isrc)%iem_lf-1)*Npos_lf+n1
                             call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_6m_av(1,1,n1der),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                                  fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                             write(def2%name,"(A)")"AvgMDA8_6month_ppb"//trim(secname)//trim(sourcename)//trim(redname)
                             def2%unit='ppb'
                             n1der = (lf_src(isrc)%iem_lf-1)*Npos_lf+n1
                             call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_6m_av_ppb(1,1,n1der),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                                  fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                          end if
                          
                          if(iotyp2ix(iou_ix)==IOU_YEAR) then
                             write(def2%name,"(A)")"AvgMDA8_full"//trim(secname)//trim(sourcename)//trim(redname)
                          else
                             write(def2%name,"(A)")"AvgMDA8"//trim(secname)//trim(sourcename)//trim(redname)
                          end if
                          def2%unit='ug/m3'
                          n1der = (lf_src(isrc)%iem_lf-1)*Npos_lf+n1
                          call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_av(1,1,n1der,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                               fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                          
                          if(iotyp2ix(iou_ix)==IOU_YEAR) then
                             !to stress difference between peak season and full year
                             write(def2%name,"(A)")"AvgMDA8_full_ppb"//trim(secname)//trim(sourcename)//trim(redname)
                          else
                             write(def2%name,"(A)")"AvgMDA8_ppb"//trim(secname)//trim(sourcename)//trim(redname)
                          end if                          
                          def2%unit='ppb'
                          n1der = (lf_src(isrc)%iem_lf-1)*Npos_lf+n1
                          call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_av_ppb(1,1,n1der,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                               fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                          
                          def2%unit='ppbdays'
                          def2%name = "SOMO35"//trim(secname)//trim(sourcename)//trim(redname)
                          n1der = (lf_src(isrc)%iem_lf-1)*Npos_lf+n1
                          call Out_netCDF(iotyp,def2,ndim_tot,1,SOMO35(1,1,n1der,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                               fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                          if(n1der==1 .and. ideriv == 1)then
                                !also save Base MDA8
                             if(iotyp2ix(iou_ix)==IOU_YEAR) then
                                def2%name="AvgMDA8_6month"
                                def2%unit='ug/m3'
                                call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_6m_av(1,1,0),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                                     fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
 
                                def2%name="AvgMDA8_6month_ppb"
                                def2%unit='ppb'
                                call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_6m_av_ppb(1,1,0),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                                     fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                             end if
                             if(iotyp2ix(iou_ix)==IOU_YEAR) then
                                def2%name="AvgMDA8_full"
                             else
                                def2%name="AvgMDA8"
                             end if
                             def2%unit='ug/m3'
                             call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_av(1,1,0,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                                  fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                             
                             if(iotyp2ix(iou_ix)==IOU_YEAR) then
                                def2%name="AvgMDA8_full_ppb"
                             else
                                def2%name="AvgMDA8_ppb"
                             end if                             
                             def2%unit='ppb'
                             call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_av_ppb(1,1,0,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                                  fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                             def2%name="MDA1"
                             def2%unit='ug/m3'
                             call Out_netCDF(iotyp,def2,ndim_tot,1,D1Max_av(1,1,0,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                                  fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                             def2%name="SOMO35"
                             def2%unit='ppbdays'
                             call Out_netCDF(iotyp,def2,ndim_tot,1,SOMO35(1,1,0,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_set%DOMAIN,&
                                  fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                          end if
                       end if
                    end if
                    if(iter == 1 .and. iotyp==IOU_YEAR .and. iout == 1)then
                       if (i==1.and.j==1.and.ideriv==1) then
                          !sum partial results from each MPI process
                          CALL MPI_ALLREDUCE(MPI_IN_PLACE,totrcemis_lf,&
                               Ncountrysectors_lf*(Ncountry_lf + Ncountry_group_lf)*NEMIS_File,&
                               MPI_REAL8, MPI_SUM,MPI_COMM_CALC,IERROR)
                      end if
                      !add one variable per emission source and put total value
                      def2%name="Total_kg"//trim(secname)//trim(sourcename)//trim(redname)
                      if(ideriv==1)iem = find_index('nox' ,EMIS_FILE(1:NEMIS_FILE))
                      if(ideriv==2)iem = find_index('voc' ,EMIS_FILE(1:NEMIS_FILE))
                      if(ideriv==3)iem = find_index('nh3' ,EMIS_FILE(1:NEMIS_FILE))
                      if(ideriv==4)iem = find_index('sox' ,EMIS_FILE(1:NEMIS_FILE))
                      def2%scale = totrcemis_lf(j,i,iem)
                      call Out_netCDF(iotyp,def2,0,1,lf_src_tot,scale,CDFtype,dimSizes_tot,dimNames_tot,&
                           fileName_given=trim(fileName),create_var_only=create_var_only,ncFileID_given=ncFileID)
                    end if

                 enddo
              enddo
              !output relative
              if (lf_spec_out(iout)%relative .and. Ndiv_rel>0) then
                 do is=1,Ncountrysectors_lf
                    isec=lf_country%sector_list(is)
                    secname = ''
                    if (isec<=NSECTORS) then
                       if(isec/=0)write(secname,"(A,I2.2)")'_sec',isec
                    else
                       secname='_'//trim(lf_sector_groups(isec-NSECTORS)%name)
                    end if
                    write(def1%unit,fmt='(A)')'ug/m3'
                    if (index(lf_spec_out(iout)%name,"ASO")>0) def1%unit='ug_PM/m3'
                    write(def1%class,fmt='(A,I0,A,I0)')'source_size_1x1'
                    def1%name =  trim(specname)//'_fraction'//trim(secname)//trim(redname)
                    if (index(lf_spec_out(iout)%name,"ASO")>0)def1%name = trim(specname)//'_PM_fraction'//trim(secname)//trim(redname)
                    scale=1.0
                    if(me==0 .and. iter==1.and. (iotyp==IOU_MON .or. iotyp==IOU_YEAR))write(*,*)'writing '//trim(def1%name)
                    call Out_netCDF(iotyp,def1,ndim,kmax,tmp_out(1,1,1,is),scale,CDFtype,dimSizes,dimNames,out_DOMAIN=lf_set%DOMAIN,&
                         fileName_given=trim(fileName),create_var_only=create_var_only,chunksizes=chunksizes,ncFileID_given=ncFileID)
                 end do
              end if
           end do
        end do
     end if
  enddo
  deallocate(tmp_out)
  deallocate(tmp_out_cntry)
  deallocate(tmp_out_base)

  n1 = 2*Npoll
  if(make_PMwater)n1=n1+1
  do ipoll = 1, n1
     do k = 1, lf_Nvertout
        do j=1,ljmax
           do i=1,limax
              lf_src_tot(i,j,k,ipoll,iou_ix) = 0.0
           enddo
        enddo
     enddo
  enddo
  do j=1,ljmax
     do i=1,limax
        lf_src_ps(i,j,iou_ix) = 0.0
     enddo
  enddo


! reset the cumulative arrays
  n1=Nsources
  if(make_PMwater)n1=n1+Nfullchem_emis
  do isrc = 1, n1
     do k = 1, lf_Nvertout
        do j=1,ljmax
           do i=1,limax
              do n=lf_src(isrc)%start, lf_src(isrc)%end
                 lf_src_acc(n,i,j,k,iou_ix)=0
               enddo
           enddo
        enddo
     enddo
  enddo

  !reset the cumulative counters
  av_fac(iotyp)=0

  LF_ncFileID_iou(iotyp) = ncFileID !to use next time

  first_call(iotyp)=.false.

  if (iotyp==IOU_HOUR_INST) then
     !reset the lf which are time tagged
     do isrc=1,Nsources
        if (lf_src(isrc)%nhour>0) then
           !during lf_src(isrc)%time_ix=current_date%hour, emissions are on
           !we reset at the very start of the hour.
           !current_date%hour is updated after emission and chemistry, and comes here
           !after hour is updated, and before new emis are included for this hour.
           if (mod(current_date%hour,24) == mod(lf_src(isrc)%time_ix,24)) then
              !reset corresponding lf to zero if time is reached.
              do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
                 do j=1,ljmax
                    do i=1,limax
                       do n = lf_src(isrc)%start, lf_src(isrc)%end
                          lf(n,i,j,k) = 0.0
                       enddo
                    enddo
                 enddo
              enddo
           end if
        end if
     end do
  end if

  call Add_2timing(NTIMING-2,tim_after,tim_before,"lf: output")
  if(DEBUGall .and. me==0)write(*,*)'end out'

! CALL MPI_BARRIER(MPI_COMM_CALC, I)

!stop
end subroutine lf_out

subroutine lf_av(dt)
  real, intent(in)    :: dt                   ! time-step used in integrations
  real :: xtot, x, O3_c
  integer ::i,j,k,kk,ii,l,n,nn,n_new,dx,dy,ix,iix,ipoll,ipoll_cfac, iou_ix, isrc
  logical :: pollwritten(2*Max_lf_spec),is_surf
  integer,save :: count_AvgMDA8_m=0,count_AvgMDA8_y=0,count_AvgMDA8_6m=0
  real :: w_m, w_y, w_6m, timefrac, Fgas
  logical, save :: first_call=.true.
  integer :: iadv_PMf

  if(DEBUGall)write(*,*)me,'start av'

  call Code_timer(tim_before)
  if(.not. lf_set%HOUR_INST .and.&
     .not. lf_set%HOUR  .and.&
     .not. lf_set%DAY   .and.&
     .not. lf_set%MONTH .and.&
     .not. lf_set%YEAR) return

  iadv_PMf = find_index('SO4', species_adv(:)%name, any_case=.true. )

  do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
     do j=1,ljmax
        do i=1,limax
           do n=1,LF_SRC_TOTSIZE
              if (abs(lf(n,i,j,k))<lf_limit) lf(n,i,j,k) = 0.0
           end do
        end do
     end do
  end do

  !do the averaging
  do iou_ix = 1, Niou_ix

     do j=1,ljmax
        do i=1,limax
           lf_src_ps(i,j,iou_ix) = lf_src_ps(i,j,iou_ix) + ps(i,j,1)
        end do
     end do
     pollwritten = .false.
     do isrc=1,Nsources_nonew
        ipoll = lf_src(isrc)%poll
        if(lf_src(isrc)%type=='country' .and. (.not.lf_fullchem .or. .not.lf_set%relative))then
           is_surf = .true.
           ipoll_cfac = ipoll + Npoll
        else
           is_surf = .false.
           ipoll_cfac = ipoll
        end if

        do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
          kk = k
          if (lf_Nvertout<KMAX_MID) kk = KMAX_MID - k + 1 !1 for surface and increasing upwards
          do j=1,ljmax
              do i=1,limax
                 xtot=0.0
                 do iix=1,lf_src(isrc)%Nsplit
                    ix = lf_src(isrc)%ix(iix)
                    if(is_surf)then
                       if(lf_src(isrc)%is_ASOA)then
                          !output the particle phase
                          Fgas = Fgas3d(ix+NSPEC_SHL,i,j, KMAX_MID)
                          xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix))/ATWAIR&
                               *roa(i,j,k,1)*1.E9*(1-Fgas)*cfac(iadv_PMf,i,j) !for ug/m3
                          !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                       else
                          !3m height cfac correction
                          xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix))/ATWAIR&
                               *roa(i,j,k,1)*1.E9* cfac(ix,i,j) !for ug/m3
                          !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                       end if
                    else                       
                       if(lf_src(isrc)%is_ASOA)then
                          !output the particle phase
                          Fgas = Fgas3d(ix+NSPEC_SHL,i,j, KMAX_MID)
                          xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix))/ATWAIR&
                               *roa(i,j,k,1)*1.E9*(1-Fgas) !for ug/m3
                          !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                       else
                          xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix))/ATWAIR&
                               *roa(i,j,k,1)*1.E9 !for ug/m3
                        !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                       end if
                    end if
                 end do
                 if(.not. pollwritten(ipoll_cfac))then !one pollutant may be used for several sources
                    if (iou_ix == iou_ix_inst) then
                       !not accumulated
                       lf_src_tot(i,j,kk,ipoll_cfac,iou_ix) = xtot
                    else
                       lf_src_tot(i,j,kk,ipoll_cfac,iou_ix) = lf_src_tot(i,j,kk,ipoll_cfac,iou_ix) + xtot
                    end if
                 end if
                 do n=lf_src(isrc)%start, lf_src(isrc)%end
                    if (iou_ix == iou_ix_inst) then
                       !not accumulated
                       lf_src_acc(n,i,j,kk,iou_ix) = xtot*lf(n,i,j,k)
                    else
                       lf_src_acc(n,i,j,kk,iou_ix)=lf_src_acc(n,i,j,kk,iou_ix)+xtot*lf(n,i,j,k)
                    end if
                    if(DEBUG .and. isnan(lf_src_acc(n,i,j,kk,iou_ix)))then
                       write(*,*)'lf is NaN ',me,isrc,lf(n,i,j,k),xtot,n,i,j,k,trim(lf_src(isrc)%species)
                       stop
                    end if
                 end do
              enddo
           enddo
        enddo
        if(isrc_new(isrc)>0)then
           !add pm25_new to pm25
           isrc_pm25_new = isrc_new(isrc)
           do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
             kk = k
             if (lf_Nvertout<KMAX_MID) kk = KMAX_MID - k + 1 !1 for surface and increasing upwards
             do j=1,ljmax
                 do i=1,limax
                    xtot=0.0
                    do iix=1,lf_src(isrc_pm25_new)%Nsplit
                       ix=lf_src(isrc_pm25_new)%ix(iix)
                       if(is_surf)then
                          !3m height cfac correction
                          xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc_pm25_new)%mw(iix))/ATWAIR&
                               *roa(i,j,k,1)*1.E9* cfac(ix,i,j) !for ug/m3
                          !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                       else
                          xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc_pm25_new)%mw(iix))/ATWAIR&
                               *roa(i,j,k,1)*1.E9 !for ug/m3
                          !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                       end if
                    end do
                    if(.not. pollwritten(ipoll_cfac))then !one pollutant may be used for several sources
                       !It is added to the _old (= not new) component, also in the inst case
                       lf_src_tot(i,j,kk,ipoll_cfac,iou_ix) = lf_src_tot(i,j,kk,ipoll_cfac,iou_ix) + xtot
                    end if
                    n_new = lf_src(isrc_pm25_new)%start
                    do n=lf_src(isrc)%start, lf_src(isrc)%end !NB: loop over isrc for pm25, not new
                       lf_src_acc(n,i,j,kk,iou_ix)=lf_src_acc(n,i,j,kk,iou_ix) + xtot*lf(n_new,i,j,k)
                       n_new = n_new + 1
                    end do
                 end do
              enddo
           enddo
        end if
        pollwritten(ipoll_cfac) = .true.

     end do
     if(make_PMwater)then
        ipoll_cfac = 2*Npoll+1
        do j=1,ljmax
           do i=1,limax
              if (iou_ix == iou_ix_inst) then
                 !not accumulated
                 lf_src_tot(i,j,1,ipoll_cfac,iou_ix) = PM25_water_rh50(i,j)
              else
                 lf_src_tot(i,j,1,ipoll_cfac,iou_ix) = lf_src_tot(i,j,1,ipoll_cfac,iou_ix) + PM25_water_rh50(i,j)
              end if
              nn=0
              if (iou_ix == iou_ix_inst) then
                 !not accumulated
                 do n=LF_SRC_TOTSIZE+1, LF_SRC_TOTSIZE+Npos_lf*Nfullchem_emis
                    nn=nn+1
                    lf_src_acc(n,i,j,1,iou_ix)=lf_PM25_water(nn,i,j)
                 end do
              else
                 do n=LF_SRC_TOTSIZE+1, LF_SRC_TOTSIZE+Npos_lf*Nfullchem_emis
                    nn=nn+1
                    lf_src_acc(n,i,j,1,iou_ix)=lf_src_acc(n,i,j,1,iou_ix)+lf_PM25_water(nn,i,j)
                 end do
              end if
           end do
        end do
     end if
     if (lf_set%MDA8) then! maximum daily eight-hour mean concentration

        if (iou_ix == 1) then !only once for (daily, monthly and yearly)
           do j = 1,ljmax
              do i = 1,limax
                 ix=O3_ix-NSPEC_SHL
                 O3_c = xn_adv(ix,i,j,KMAX_MID) * cfac(ix,i,j) * roa(i,j,KMAX_MID,1) * 1.0e9 * species_adv(ix)%molwt/ATWAIR
                 hourM(0,i,j) = hourM(0,i,j) + O3_c
                 !how much the hourly values will change for a small change of emissions * 100% (summed over hour, normalize later)
                 do n=1, Npos_lf*Nfullchem_emis
                    hourM(n,i,j) = hourM(n,i,j) + O3_c * lf(lf_src(isrc_O3)%start+n-1,i,j,KMAX_MID)
                 end do
              end do
           end do
           do j = 1,ljmax
              do i = 1,limax
                 ix=O3_ix-NSPEC_SHL
                 O3_c = xn_adv(ix,i,j,KMAX_MID) * cfac(ix,i,j) * 1.0e9
                 hourM_ppb(0,i,j) = hourM_ppb(0,i,j) + O3_c
                 !how much the hourly values will change for a small change of emissions * 100% (summed over hour, normalize later)
                 do n=1, Npos_lf*Nfullchem_emis
                    hourM_ppb(n,i,j) = hourM_ppb(n,i,j) + O3_c * lf(lf_src(isrc_O3)%start+n-1,i,j,KMAX_MID)
                 end do
              end do
           end do

           if (current_date%seconds == 0 .and. .not. first_call) then
             !one hour has past since last time here
              !save last hour average
              ii = mod(current_date%hour,8) + 1 !note: this works only because 24 is a multiple of 8!
              timefrac = dt_advec/3600.0 ! inverse of number of timesteps in an hour
              D8M(:,:,:,ii) =  hourM(:,:,:) * timefrac !save all hourly averages
              hourM(:,:,:) = 0.0
              D8M_ppb(:,:,:,ii) =  hourM_ppb(:,:,:) * timefrac !save all hourly averages
              hourM_ppb(:,:,:) = 0.0
              !one hour has past since last time here
              ! update max value since 01:00
              do j = 1,ljmax
                 do i = 1,limax
                    
                    if (D8M(0,i,j,ii) > D1Max(0,i,j)) then
                       !a new max is found. Update for all fractions too
                       do n=0, lf_src(isrc_O3)%Npos*Nfullchem_emis
                          D1Max(n,i,j) = D8M(n,i,j,ii)
                       end do
                    end if
                    
                    x = 0.0
                    do l = 1, 8
                       x = x + D8M(0,i,j,l) * 0.125 ! 8 hour average
                    end do
                    if(x > D8Max(0,i,j)) then
                       !a new max is found. Update for all fractions too
                       do n=0, lf_src(isrc_O3)%Npos*Nfullchem_emis
                          D8Max(n,i,j) = D8M(n,i,j,1) * 0.125 ! 8 hour average, contribution from first hour
                       end do
                       do l = 2, 8
                          do n=0, lf_src(isrc_O3)%Npos*Nfullchem_emis
                             D8Max(n,i,j) = D8Max(n,i,j) + D8M(n,i,j,l) * 0.125 ! 8 hour average, contribution from second to eighth hour
                          end do
                       end do
                    end if
                 end do
              end do
              do j = 1,ljmax
                 do i = 1,limax
                    x = 0.0
                    do l = 1, 8
                       x = x + D8M_ppb(0,i,j,l) * 0.125 ! 8 hour average
                    end do
                    if(x > D8Max_ppb(0,i,j)) then
                       !a new max is found. Update for all fractions too
                       do n=0, lf_src(isrc_O3)%Npos*Nfullchem_emis
                          D8Max_ppb(n,i,j) = D8M_ppb(n,i,j,1) * 0.125 ! 8 hour average, contribution from first hour
                       end do
                       do l = 2, 8
                          do n=0, lf_src(isrc_O3)%Npos*Nfullchem_emis
                             D8Max_ppb(n,i,j) = D8Max_ppb(n,i,j) + D8M_ppb(n,i,j,l) * 0.125 ! 8 hour average, contribution from second to eighth hour
                          end do
                       end do
                    end if
                 end do
              end do
           end if
        end if
        if (current_date%seconds == 0 .and. current_date%hour == 0 .and. .not. first_call) then
           !end of day, save the values
           if (iotyp2ix(iou_ix)==IOU_DAY) then
              D8Max_av(:,:,:,iou_ix)=0.0
              D8Max_av_ppb(:,:,:,iou_ix)=0.0
              D1Max_av(:,:,:,iou_ix)=0.0
              SOMO35(:,:,:,iou_ix)=0.0
           end if
           !NB: at the end of the first day (day 2 hour 00:00), we actually start to write in the next month
           if (current_date%day == 2 .and. iotyp2ix(iou_ix)==IOU_MON) then
              !new month
              count_AvgMDA8_m = 0
              D8Max_av(:,:,:,iou_ix)=0.0
              D8Max_av_ppb(:,:,:,iou_ix)=0.0
              D1Max_av(:,:,:,iou_ix)=0.0
              SOMO35(:,:,:,iou_ix)=0.0
           end if
           
           if (current_date%day == 2 .and. current_date%month == 4 .and. iotyp2ix(iou_ix)==IOU_YEAR) then
              !new yearly max peak season
              count_AvgMDA8_6m = 0
              D8Max_6m_av(:,:,:)=0.0
              D8Max_6m_av_ppb(:,:,:)=0.0
              !SOMO35(:,:,:,iou_ix)=0.0 !for SOMO35 the integral goes over the entire year
              !for MDA1 the integral goes over the entire year
           end if

           if(iotyp2ix(iou_ix)==IOU_MON)count_AvgMDA8_m = count_AvgMDA8_m + 1
           if(iotyp2ix(iou_ix)==IOU_YEAR)count_AvgMDA8_y = count_AvgMDA8_y + 1
           if(iotyp2ix(iou_ix)==IOU_YEAR)count_AvgMDA8_6m = count_AvgMDA8_6m + 1 !not used after September

           w_m = 1.0/count_AvgMDA8_m
           w_y = 1.0/count_AvgMDA8_y
           w_6m = 1.0/count_AvgMDA8_6m
           do j = 1,ljmax
              do i = 1,limax
                 if (iotyp2ix(iou_ix)==IOU_DAY)then
                    do n=0, Npos_lf*Nfullchem_emis
                       D8Max_av(i,j,n,iou_ix) =  D8Max(n,i,j)
                       D1Max_av(i,j,n,iou_ix) =  D1Max(n,i,j)
                       D8Max_av_ppb(i,j,n,iou_ix) =  D8Max_ppb(n,i,j)
                    end do
                    !NB: if and only if D8Max_ppb>35 , all the SOMO35 fractions must be updated
                    if (D8Max_ppb(0,i,j)>35.0) then
                       SOMO35(i,j,0,iou_ix) =  SOMO35(i,j,0,iou_ix) + D8Max_ppb(0,i,j) !integral over days
                       do n=1, Npos_lf*Nfullchem_emis
                          !NB: derivatives have no threshold
                          SOMO35(i,j,n,iou_ix) =  SOMO35(i,j,n,iou_ix) + D8Max_ppb(n,i,j) !integral over days
                       end do
                    end if
                 else if(iotyp2ix(iou_ix)==IOU_MON)then
                    do n=0, Npos_lf*Nfullchem_emis
                       D8Max_av(i,j,n,iou_ix) =  (1.0-w_m) * D8Max_av(i,j,n,iou_ix) + w_m * D8Max(n,i,j)
                       D1Max_av(i,j,n,iou_ix) =  (1.0-w_m) * D1Max(n,i,j) + w_m * D1Max(n,i,j)
                       D8Max_av_ppb(i,j,n,iou_ix) =  (1.0-w_m) * D8Max_av_ppb(i,j,n,iou_ix) + w_m * D8Max_ppb(n,i,j)
                    end do
                    if (D8Max_ppb(0,i,j)>35.0) then
                       SOMO35(i,j,0,iou_ix) =  SOMO35(i,j,0,iou_ix) + D8Max_ppb(0,i,j)-35.0 !integral over days
                       do n=1, Npos_lf*Nfullchem_emis
                          !NB: derivatives have no threshold
                          SOMO35(i,j,n,iou_ix) =  SOMO35(i,j,n,iou_ix) + D8Max_ppb(n,i,j) !integral over days
                       end do
                    end if
                 else if (iotyp2ix(iou_ix)==IOU_YEAR)then
                    do n=0, Npos_lf*Nfullchem_emis
                       D8Max_av(i,j,n,iou_ix) =  (1.0-w_y) * D8Max_av(i,j,n,iou_ix) + w_y * D8Max(n,i,j)
                       D1Max_av(i,j,n,iou_ix) =  (1.0-w_y) * D1Max_av(i,j,n,iou_ix) + w_y * D1Max(n,i,j)
                       D8Max_av_ppb(i,j,n,iou_ix) =  (1.0-w_y) * D8Max_av_ppb(i,j,n,iou_ix) + w_y * D8Max_ppb(n,i,j)
                    end do
                    if(current_date%month>=4 .and. current_date%month<=9)then
                       !we keep only days from April to September for peak season
                       do n=0, Npos_lf*Nfullchem_emis
                          D8Max_6m_av(i,j,n) =  (1.0-w_6m) * D8Max_6m_av(i,j,n) + w_6m * D8Max(n,i,j)
                          D8Max_6m_av_ppb(i,j,n) =  (1.0-w_6m) * D8Max_6m_av_ppb(i,j,n) + w_6m * D8Max_ppb(n,i,j)
                       end do
                    end if
                    
                    !NB: if and only if D8Max_ppb>35 , all the SOMO35 fractions must be updated
                    if (D8Max_ppb(0,i,j)>35.0) then
                       SOMO35(i,j,0,iou_ix) =  SOMO35(i,j,0,iou_ix) + D8Max_ppb(0,i,j)-35.0 !integral over days
                       do n=1, Npos_lf*Nfullchem_emis
                          !NB: derivatives have no threshold
                          SOMO35(i,j,n,iou_ix) =  SOMO35(i,j,n,iou_ix) + D8Max_ppb(n,i,j)!integral over days
                       end do
                    end if
                 end if
              end do
           end do
           if(iou_ix == Niou_ix) D1Max = 0.0 !we are ready to start a new day
           if(iou_ix == Niou_ix) D8Max = 0.0 !we are ready to start a new day
           if(iou_ix == Niou_ix) D8Max_ppb = 0.0 !we are ready to start a new day
        end if
     end if
  end do

  av_fac=av_fac+1
  first_call=.false.

  call Add_2timing(NTIMING-9,tim_after,tim_before,"lf: averaging")
  if(DEBUGall)write(*,*)me,' end av'

end subroutine lf_av

subroutine lf_adv_x(fluxx,i,j,k)
  real, intent(in)::fluxx(NSPEC_ADV,-1:LIMAX+1)
  integer, intent(in)::i,j,k
  real ::x,xn,xx,f_in,inv_tot
  integer ::n,ii,iix,ix,dx,dy,isrc,dp,dm,nstart,nend
  if(DEBUGall .and. me==0)write(*,*)'start advx'

  if(i==li0)then
     !copy small part (could be avoided, but simpler to copy)
     !note that the right hand side of the lf equations must contain unupdated values, therefore values for j-1 must be buffered
     do ii=li0,li1
        do n=1,LF_SRC_TOTSIZE
           loc_frac_src_1d(n,ii) = lf(n,ii,j,k)
        enddo
     enddo
     if (li0 == 2) then
        !track left boundary values
        do n = 1, nbc
           loc_frac_src_1d(BC_ix(n), 1) = 1.0
        enddo
     end if
     if (li1 < LIMAX) then
        !track right boundary values
        do n = 1, nbc
           loc_frac_src_1d(BC_ix(n), LIMAX) = 1.0
        enddo
     end if
  end if

  call Code_timer(tim_before)
  do isrc=1,Nsources
     xn=0.0
     x=0.0
     xx=0.0
     !positive x or xx means incoming, negative means outgoing
     do iix=1,lf_src(isrc)%Nsplit
        ix=lf_src(isrc)%ix(iix)
        xn=xn+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
        x=x-xm2(i,j)*fluxx(ix,i)*lf_src(isrc)%mw(iix)!flux through "East" face (Right)
        xx=xx+xm2(i,j)*fluxx(ix,i-1)*lf_src(isrc)%mw(iix)!flux through "West" face (Left)
     end do
     !NB: here xn already includes the fluxes. Remove them!
     xn=xn-xx-x
     xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux
     f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
     inv_tot=1.0/(xn+f_in+1.e-40)!incoming dilutes

     x =max(0.0,x)*inv_tot!factor due to flux through "East" face (Right)
     xx=max(0.0,xx)*inv_tot!factor due to flux through "West" face (Left)
     xn = xn * inv_tot
     !often either x or xx is zero
     if(lf_src(isrc)%type=='country')then
        nstart = lf_src(isrc)%start
        nend = lf_src(isrc)%start + (Ncountry_lf + Ncountry_group_lf)*Ncountrysectors_lf - 1
        if(x>1.E-20)then
           do n = nstart, nend
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,i+1)*x
           enddo
           if(xx>1.E-20)then
              do n = nstart, nend
                 lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n,i-1)*xx
              enddo
           end if
        else if (xx>1.E-20)then
           do n = nstart, nend
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,i-1)*xx
           enddo
        end if
     end if
     if(lf_src(isrc)%type=='relative' .or. (lf_set%relative .and. lf_fullchem))then
        !The relative position of the source is either the same for i and incoming lf, or differs by 1
        dp = 0 ! for left boundary
        dm = 0 ! for right boundary
        if(mod(i_fdom(i)-1,lf_src(isrc)%res)==0)dp=1
        if(mod(i_fdom(i),lf_src(isrc)%res)==0)dm=1
        nstart = lf_src(isrc)%start
        !the relative for fullchem are placed after the regular countries, as a set of "countries"
        if (lf_fullchem) nstart = nstart + (Ncountry_lf + Ncountry_group_lf)*Ncountrysectors_lf !start for relative
        do while (nstart <= lf_src(isrc)%end - Ndiv_rel*Ndiv_rel + 1) !loop over sectors, since they are under the same isrc
           if(x>1.E-20)then
              n = nstart
              if (dm==0) then
                 do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                       lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,i+1)*x
                       n=n+1
                    enddo
                 enddo
              else
                 do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    lf(n,i,j,k) = lf(n,i,j,k)*xn ! when dx=-lf_src(isrc)%dist and dm=1, there are no local fractions to transport
                    n=n+1
                    do dx=-lf_src(isrc)%dist+1,lf_src(isrc)%dist
                       lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n-dm,i+1)*x
                       n=n+1
                    enddo
                 enddo
              end if
              
              if(xx>1.E-20)then
                 n = nstart
                 do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                       lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+dp,i-1)*xx
                       n=n+1
                    enddo
                    if (dp==0) lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+dp,i-1)*xx
                    n=n+1! when dx=lf_src(isrc)%dist and dp=1 there are no local fractions to transport
                 enddo
              end if
           else if (xx>1.E-20)then
              n = nstart
              do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                 do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                    lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n+dp,i-1)*xx
                    n=n+1
                 enddo
                 lf(n,i,j,k) = lf(n,i,j,k)*xn! when dx=lf_src(isrc)%dist  and dp=1 there are no local fractions to transport
                 if (dp==0) lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+dp,i-1)*xx
                 n=n+1
              enddo
           else
              !nothing to do if no incoming fluxes
           end if
           if (lf_fullchem) then
              nstart = nstart + Ndiv_rel*Ndiv_rel !start for next sector
           else
              nstart = lf_src(isrc)%end + 1 !just do once
           end if
        end do
     end if
  enddo

  call Add_2timing(NTIMING-8,tim_after,tim_before,"lf: adv_x")

  if(DEBUGall .and. me==0)write(*,*)'end advx'
end subroutine lf_adv_x

subroutine lf_adv_y(fluxy,i,j,k)
  real, intent(in)::fluxy(NSPEC_ADV,-1:LJMAX+1)
  integer, intent(in)::i,j,k
  real ::x,xn,xx,f_in,inv_tot
  integer ::n,jj,iix,ix,dx,dy,isrc,dp,dm,nstart,nend
  if(DEBUGall .and. me==0)write(*,*)'start advy'

  if(j==lj0)then
     !copy small part (could be avoided, but simpler to copy)
     !note that the right hand side of the lf equations must contain unupdated values, therefore values for i-1 must be buffered
     do jj=lj0,lj1
        do n=1,LF_SRC_TOTSIZE
           loc_frac_src_1d(n,jj) = lf(n,i,jj,k)
        enddo
     enddo
     if (lj0 == 2) then
        !track lower boundary values
        do n = 1, nbc
           loc_frac_src_1d(BC_ix(n), 1) = 1.0
        enddo
     end if
     if (lj1 < LJMAX) then
        !track right boundary values
        do n = 1, nbc
           loc_frac_src_1d(BC_ix(n), LJMAX) = 1.0
        enddo
     end if
  end if

  call Code_timer(tim_before)
  do isrc=1,Nsources
     xn=0.0
     x=0.0
     xx=0.0
     !positive x or xx means incoming, negative means outgoing
     do iix=1,lf_src(isrc)%Nsplit
        ix=lf_src(isrc)%ix(iix)
        xn=xn+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
        x=x-xm2(i,j)*fluxy(ix,j)*lf_src(isrc)%mw(iix)!flux through "North" face (Up)
        xx=xx+xm2(i,j)*fluxy(ix,j-1)*lf_src(isrc)%mw(iix)!flux through "South" face (Bottom)
     end do
     !NB: here xn already includes the fluxes. Remove them!
     xn=xn-xx-x
     xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux
     f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
     inv_tot=1.0/(xn+f_in+1.e-40)!incoming dilutes

     x =max(0.0,x)*inv_tot!factor due to flux through "East" face (Right)
     xx=max(0.0,xx)*inv_tot!factor due to flux through "West" face (Left)
     xn = xn * inv_tot
     if(lf_src(isrc)%type=='country')then
        nstart = lf_src(isrc)%start
        nend = lf_src(isrc)%start + (Ncountry_lf + Ncountry_group_lf)*Ncountrysectors_lf - 1
        !often either x or xx is zero
        if(x>1.E-20)then
           do n = nstart, nend
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,j+1)*x
           enddo
           if(xx>1.E-20)then
              do n = nstart, nend
                 lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n,j-1)*xx
              enddo
           end if
        else if (xx>1.E-20)then
           do n = nstart, nend
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,j-1)*xx
           enddo
        end if
     end if
     if(lf_src(isrc)%type=='relative' .or. (lf_set%relative .and. lf_fullchem))then
        !The relative position of the source is either the same for j and incoming lf, or differs by 1
        dp = 0 ! for lower boundary
        dm = 0 ! for upper boundary
        if(mod(j_fdom(j)-1,lf_src(isrc)%res)==0)dp=Ndiv_rel
        if(mod(j_fdom(j),lf_src(isrc)%res)==0)dm=Ndiv_rel
        nstart = lf_src(isrc)%start
        !the relative for fullchem are placed after the regular countries
        if (lf_fullchem) nstart = nstart + (Ncountry_lf + Ncountry_group_lf) *Ncountrysectors_lf !start for relative
        do while (nstart <= lf_src(isrc)%end - Ndiv_rel*Ndiv_rel + 1)
           if(x>1.E-20)then
              n = nstart
              if(dm==0)then
                 do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                       lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,j+1)*x
                       n=n+1
                    enddo
                 enddo
              else
                 dy = -lf_src(isrc)%dist
                 do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    lf(n,i,j,k) = lf(n,i,j,k)*xn
                    n=n+1
                 enddo
                 do dy=-lf_src(isrc)%dist+1,lf_src(isrc)%dist
                    do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                       lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n-dm,j+1)*x
                       n=n+1
                    enddo
                 enddo
              end if
              if(xx>1.E-20)then
                 n = nstart
                 if(dp==0)then
                    do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                       do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                          lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n,j-1)*xx
                          n=n+1
                       enddo
                    enddo
                 else
                    do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                       do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                          lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+dp,j-1)*xx
                          n=n+1
                       enddo
                    enddo
                 end if
              end if
           else if (xx>1.E-20)then
              n = nstart
              if(dp==0)then
                 do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                       lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,j-1)*xx
                       n=n+1
                    enddo
                 enddo
              else
                 do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                    do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                       lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n+dp,j-1)*xx
                       n=n+1
                    enddo
                 enddo
                 dy=lf_src(isrc)%dist
                 do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    lf(n,i,j,k) = lf(n,i,j,k)*xn
                    n=n+1
                 enddo
              end if
           else
           !nothing to do if no incoming fluxes
           end if
           if (lf_fullchem) then
              nstart = nstart + Ndiv_rel*Ndiv_rel !jump to next sector
           else
              nstart = lf_src(isrc)%end + 1 !just do once
           end if
        end do
     end if

  enddo
  call Add_2timing(NTIMING-7,tim_after,tim_before,"lf: adv_y")
  if(DEBUGall .and. me==0)write(*,*)'end advy'

end subroutine lf_adv_y

subroutine lf_adv_k(fluxk,i,j)
    real, intent(in)::fluxk(NSPEC_ADV,KMAX_MID)
    integer, intent(in)::i,j
    real ::x,xn,xx,f_in,inv_tot
    integer ::n,k,iix,ix,dx,dy,isrc,is
    real loc_frac_src_km1(Npos_lf,KMAX_MID-lf_Nvert+1:KMAX_MID)
  if(DEBUGall .and. me==0)write(*,*)'start advk'

    call Code_timer(tim_before)
    !need to be careful to always use non-updated values on the RHS
    
    do isrc=1,Nsources
       do k = KMAX_MID-lf_Nvert+2,KMAX_MID
          is = 0
          do n = lf_src(isrc)%start, lf_src(isrc)%end
             is = is + 1
             loc_frac_src_km1(is,k)=lf(n,i,j,k-1) !NB: k is shifted by 1 in loc_frac_src_km1
          enddo
       enddo
       loc_frac_src_km1(1:is,KMAX_MID-lf_Nvert+1) = 0.0 ! everything above is not tracked, zero local fractions coming from above

       !we assume that all O3 above the LF window is from top
       do n = 1, nstratos
          if(Stratos_ix(n) >= lf_src(isrc)%start .and.  Stratos_ix(n) <=lf_src(isrc)%end) then
             is = Stratos_ix(n)- lf_src(isrc)%start + 1 
             loc_frac_src_km1(Stratos_ix(n),KMAX_MID-lf_Nvert+1) = 1.0 ! NB: k is shifted by 1 in loc_frac_src_km1, i.e. this is level k=1 if lf_Nvert = KMAX_MID - 1
          end if
       end do
       if (lf_set%EmisDer_all) then
          !also other species than O3 will be included
          do n = 1, nbc
             if (BC_ix(n) >= lf_src(isrc)%start .and.  BC_ix(n) <=lf_src(isrc)%end) then
                is = BC_ix(n) - lf_src(isrc)%start + 1 
                loc_frac_src_km1(BC_ix(n),KMAX_MID-lf_Nvert+1) = 1.0 ! NB: k is shifted by 1 in loc_frac_src_km1, i.e. this is level k=1 if lf_Nvert = KMAX_MID - 1
             end if
          end do
       end if

       do k = KMAX_MID-lf_Nvert+1,KMAX_MID!k is increasing-> can use k+1 to access non-updated value
          xn=0.0
          x=0.0
          xx=0.0
          !positive x or xx means incoming, negative means outgoing
          do iix=1,lf_src(isrc)%Nsplit
             ix=lf_src(isrc)%ix(iix)
             xn=xn+xn_2d(ix,k)*lf_src(isrc)%mw(iix) !xn_2d are the  unupdated values NB: not same ix convention than in chemistry!
             if(k<KMAX_MID)x=x-dhs1i(k+1)*fluxk(ix,k+1)*lf_src(isrc)%mw(iix)
             xx=xx+dhs1i(k+1)*fluxk(ix,k)*lf_src(isrc)%mw(iix)
          end do
          xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux
          f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
          inv_tot=1.0/(xn+f_in+1.e-40)!incoming dilutes

          x =max(0.0,x)*inv_tot!factor due to flux through bottom face
          xx=max(0.0,xx)*inv_tot!factor due to flux through top face
          if(k==KMAX_MID) x = 0.0 !no fraction coming from surface
          xn = xn * inv_tot
          !often either x or xx is zero
          if(x>1.E-20)then
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                lf(n,i,j,k) = lf(n,i,j,k)*xn +lf(n,i,j,k+1)*x
             enddo
             if(xx>1.E-20)then
                is = 0
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   is = is + 1
                   lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_km1(is,k)*xx !NB: lf already multiplied by xn here
                enddo
             end if
          else if (xx>1.E-20)then
          is = 0
          do n = lf_src(isrc)%start, lf_src(isrc)%end
               is = is + 1
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_km1(is,k)*xx
             enddo
          else
             !nothing to do if no incoming fluxes
          end if
       enddo

    end do

    if(DEBUGall .and. me==0)write(*,*)'end advk'
    call Add_2timing(NTIMING-6,tim_after,tim_before,"lf: adv_k")
  end subroutine lf_adv_k

  subroutine lf_adv_k_2nd(fluxk,i,j,sdot,dt_s,alfnew,alfbegnew,alfendnew)
    !takes into account the terms that arise from the second order Bott scheme:
    !the lf flux contributions are weighted by zzfl1, zzfl2 and zzfl3
    real, intent(in)::fluxk(NSPEC_ADV,KMAX_MID)
    real, intent(in)::alfnew(9,2:KMAX_MID,0:1),alfbegnew(3),alfendnew(3),sdot(0:LIMAX*LJMAX*KMAX_BND-1),dt_s
    integer, intent(in)::i,j
    real ::x,xn,xx,x0,xx0,f_in,inv_tot,fc1,fc2,fc3,xn_post
    integer ::n,k,iix,ix,dx,dy,isrc,n1k,k1,klimlow,klimhig
    real loc_frac_src(Npos_lf,KMAX_MID-lf_Nvert-1:KMAX_MID+2)
    real :: w1,w2,w3 !weights from the used level for making the flux
    !NB: level used for *incoming* fluxes:
    !    fluxes incoming k from k-1: uses level k-2,k-1 and k to make flux(k)
    !    fluxes incoming k from k+1: uses level k,k+1 and k+2 to make flux(k+1)
    ! Eta=A/Pref+B , P = A + B*PS, high altitude means A and B and Eta decrease. Positive Etadot means downwind.
    real fc(KMAX_MID),zzfl1(KMAX_MID-lf_Nvert:KMAX_MID),zzfl2(KMAX_MID-lf_Nvert:KMAX_MID),zzfl3(KMAX_MID-lf_Nvert:KMAX_MID)
    real :: x1,x2,x3,xx1,xx2,xx3,fk1,lfmax
    integer :: is
    
    if(DEBUGall .and. me==0)write(*,*)'start adv_k_2nd'

    if (.not. lf_fullchem) then
       ! we do not want second order vertical lf advection if there is no O3 anyway:
       ! it creates negative values and uses more cpu
       call lf_adv_k(fluxk,i,j)
       return
    end if
    
    call Code_timer(tim_before)

    lfmax=10.0 !limitation for extreme values
    !note about extreme values:
    ! in some cases, Bott will almost entirely empty a cell. If it is 99.999% or 99.9 % may be dependent on a neighboring gridcells concentration.
    ! Since the final concentration is almost zero, the derivative divided by that amount can be large.
    
    !we copy paste from advvk to get zzfl1,zzfl2,zzfl3
    !we must make both the values used for making fluxk(ix,k) and fluxk(ix,k+1)
    
    do k = 1,KMAX_MID-1
       fc(k) = sdot(k*LIMAX*LJMAX)*dt_s
    end do
    fc(KMAX_MID) = -1.
    klimlow = 1
    if(fc(1).ge.0.)klimlow=2
    klimhig = KMAX_MID-1
    if(fc(KMAX_MID-1).lt.0.)klimhig = KMAX_MID-2
    
    !could merge this loop with next one for better performance. Need only two k values for zzfl at a time
    zzfl1 = 0
    zzfl2 = 0
    zzfl3 = 0
    do k = KMAX_MID-lf_Nvert,KMAX_MID
       if(k==1 .or. k==KMAX_MID) then
          !all zzfl default
       else if(k==2) then
          if(fc(1).ge.0.)then
             fc1 = fc(1)
             fc2 = fc1*fc1
             fc3 = fc1*fc2
             zzfl2(k) = alfbegnew(1)*fc1 + alfbegnew(2)*fc2 + alfbegnew(3)*fc3
             zzfl3(k) = alfnew(7,2,0)*fc1 + alfnew(8,2,0)*fc2 + alfnew(9,2,0)*fc3
          end if
       else if(k==KMAX_MID-1 .and. fc(KMAX_MID-1).lt.0.) then
          fc1 = fc(KMAX_MID-1)
          fc2 = fc1*fc1
          fc3 = fc1*fc2
          zzfl1(k) = alfnew(1,KMAX_MID,1)*fc1 + alfnew(2,KMAX_MID,1)*fc2 + alfnew(3,KMAX_MID,1)*fc3
          zzfl2(k) = alfendnew(1)*fc1 + alfendnew(2)*fc2 + alfendnew(3)*fc3
       else
          fc1 = fc(k)
          fc2 = fc1*fc1
          fc3 = fc1*fc2
          n1k = 0
          if(fc1.lt.0)n1k=1
          zzfl1(k) = alfnew(1,k+1,n1k)*fc1         &
               + alfnew(2,k+1,n1k)*fc2         &
               + alfnew(3,k+1,n1k)*fc3
          zzfl2(k) = alfnew(4,k+1,n1k)*fc1         &
               + alfnew(5,k+1,n1k)*fc2         &
               + alfnew(6,k+1,n1k)*fc3
          zzfl3(k) = alfnew(7,k+1,n1k)*fc1         &
               + alfnew(8,k+1,n1k)*fc2         &
               + alfnew(9,k+1,n1k)*fc3
          k1 = k-1+n1k
          
       end if
    end do
    
    do isrc = 1, Nsources       
       !need to be careful to always use non-updated values on the RHS
       do k = KMAX_MID-lf_Nvert-1,KMAX_MID+2
          if(k>=KMAX_MID-lf_Nvert+1 .and. k<=KMAX_MID)then
             is = 0
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                is = is + 1
                loc_frac_src(is,k)=lf(n,i,j,k)
             enddo
          else
             is = 0
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                is = is + 1
                loc_frac_src(is,k)=0.0
             enddo
          end if
       enddo
       
       !for stratosphere tracking: everything above the tracking region has lf=1
       do n = 1, nstratos
          if (Stratos_ix(n) >= lf_src(isrc)%start .and.  Stratos_ix(n) <=lf_src(isrc)%end) then
             is = Stratos_ix(n)- lf_src(isrc)%start + 1 
             loc_frac_src(is,KMAX_MID-lf_Nvert) = 1.0 ! this is level k=1 if lf_Nvert = KMAX_MID - 1
          end if
       end do       
       
       fk1=0.0
       do k = KMAX_MID-lf_Nvert,KMAX_MID-1
          xn=0.0
          x=0.0
          x1=0.0
          x2=0.0
          x3=0.0
          !flux between k and k+1
          !x, x1, x2, x3 >0
          !dhs1i factor: use (k+1) if xn for level k is treated.
          
          if(fc(k).lt.0)then
             !fluxk(:,k+1) is made from xn(k), xn(k+1), xn(k+2)
             !fluxk(:,k+1)<0 into k from k+1
             do iix=1,lf_src(isrc)%Nsplit
                ix=lf_src(isrc)%ix(iix)
                if(k>=KCHEMTOP) xn=xn+xn_2d(ix,k)*lf_src(isrc)%mw(iix) !xn_2d are the  unupdated values NB: not same ix convention than in chemistry!
                x=x-fluxk(ix,k+1)*lf_src(isrc)%mw(iix)
                if(k>=KCHEMTOP) x1=x1+zzfl1(k)*xn_2d(ix,k)*lf_src(isrc)%mw(iix)
                x2=x2+zzfl2(k)*xn_2d(ix,k+1)*lf_src(isrc)%mw(iix)
                if(k<=KMAX_MID-2) x3=x3+zzfl3(k)*xn_2d(ix,k+2)*lf_src(isrc)%mw(iix)
             end do
             !normally x1+x2+x3=x. The exception is when flux limitations are in effect; in such cases we assume all flux is dependent on upwind cell only
             !also if different splits have different fluxes, this may be the case.
             if(1e-3<abs(x1+x2+x3-x)/(abs(x)+1e-20)) then
                x1=0
                x2=x
                x3=0
             end if
             
             !out of k+1 (first time modified)
             is = 0
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                is = is + 1
                lf(n,i,j,k+1) = -(x1*loc_frac_src(is,k)+x2*loc_frac_src(is,k+1)+x3*loc_frac_src(is,k+2))
             enddo
             xn_post=xn+(fk1+(x1+x2+x3))*dhs1i(k+1)
             fk1=-(x1+x2+x3)
             
             if (k>KMAX_MID-lf_Nvert) then
                !into k  (already initialized. Divide by final concentration)
                if (xn_post>1e-20) then
                   is = 0
                   do n = lf_src(isrc)%start, lf_src(isrc)%end
                      is = is + 1
                      lf(n,i,j,k) = (xn*loc_frac_src(is,k)+(lf(n,i,j,k)+x1*loc_frac_src(is,k)+x2*loc_frac_src(is,k+1)+x3*loc_frac_src(is,k+2))*dhs1i(k+1))/xn_post
                      lf(n,i,j,k) =max(-lfmax,min(lfmax,lf(n,i,j,k)))
                   enddo
                else
                   !if there are no concentration, we set the fractions to zero
                   do n = lf_src(isrc)%start, lf_src(isrc)%end
                      lf(n,i,j,k) = 0.0
                   enddo
                end if
             end if
          else
             !fluxk(:,k+1) is made from xn(k-1), xn(k), xn(k+1)
             !fluxk(:,k+1)>0 out of k into k+1
             do iix=1,lf_src(isrc)%Nsplit
                ix=lf_src(isrc)%ix(iix)
                if(k>=KCHEMTOP) xn=xn+xn_2d(ix,k)*lf_src(isrc)%mw(iix) !xn_2d are the  unupdated values NB: not same ix convention than in chemistry!
                x=x+fluxk(ix,k+1)*lf_src(isrc)%mw(iix)
                if(k>=KCHEMTOP+1) x1=x1+zzfl1(k)*xn_2d(ix,k-1)*lf_src(isrc)%mw(iix)
                if(k>=KCHEMTOP) x2=x2+zzfl2(k)*xn_2d(ix,k)*lf_src(isrc)%mw(iix)
                x3=x3+zzfl3(k)*xn_2d(ix,k+1)*lf_src(isrc)%mw(iix)
             end do
             !normally x1+x2+x3=x/dhs1i(k+1). The exception is when flux limitations are in effect; in such cases we assume all flux is dependent on upwind cell only
             if(1e-3<abs(x1+x2+x3-x)/(abs(x)+1e-20)) then
                x1=0
                x2=x
                x3=0
             end if
             
             !into k+1 (first time, initialize)
             is = 0
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                is = is + 1
                lf(n,i,j,k+1) =x1*loc_frac_src(is,k-1)+x2*loc_frac_src(is,k)+x3*loc_frac_src(is,k+1)
             enddo
             
             xn_post=xn+(fk1-(x1+x2+x3))*dhs1i(k+1)
             fk1=(x1+x2+x3)
             if (k>KMAX_MID-lf_Nvert) then
                !out of k (already initialized. Divide by final concentration)
                if (xn_post>1e-20) then
                   is = 0
                   do n = lf_src(isrc)%start, lf_src(isrc)%end
                      is = is + 1
                      lf(n,i,j,k) =(xn*loc_frac_src(is,k)+(lf(n,i,j,k)-(x1*loc_frac_src(is,k-1)+x2*loc_frac_src(is,k)+x3*loc_frac_src(is,k+1)))*dhs1i(k+1))/xn_post
                      lf(n,i,j,k) =max(-lfmax,min(lfmax,lf(n,i,j,k)))
                   enddo
                else
                   !if there are no concentration, we set the fractions to zero
                   do n = lf_src(isrc)%start, lf_src(isrc)%end
                      lf(n,i,j,k) = 0.0
                   enddo
                end if
             end if
          end if
          
          if(k==KMAX_MID-1)then
             !we won't come back to treat KMAX_MID. Need to divid by xn_post(k+1) now
             xn=0.0
             do iix=1,lf_src(isrc)%Nsplit
                ix=lf_src(isrc)%ix(iix)
                xn=xn+xn_2d(ix,k+1)*lf_src(isrc)%mw(iix) !xn_2d are the  unupdated values NB: not same ix convention than in chemistry!
             end do
             xn_post=xn+fk1*dhs1i(k+2)
             if (xn_post>1e-20) then
                is = 0
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   is = is + 1
                   lf(n,i,j,k+1) =(xn*loc_frac_src(is,k+1)+lf(n,i,j,k+1)*dhs1i(k+2))/xn_post
                   lf(n,i,j,k+1) =max(-lfmax,min(lfmax,lf(n,i,j,k+1)))
                enddo
             else
                !if there are no concentration, we set the fractions to zero
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   lf(n,i,j,k+1) = 0.0
                enddo
             end if
          end if
       end do
    end do

    if(DEBUGall .and. me==0)write(*,*)'end adv_k_2nd'
    call Add_2timing(NTIMING-6,tim_after,tim_before,"lf: adv_k")
  end subroutine lf_adv_k_2nd

  subroutine lf_diff(i,j,dt_diff,ndiff)

    implicit none
    integer, intent(in) :: i,j,ndiff
    real, intent(in) :: dt_diff
    real :: xn_k(LF_SRC_TOTSIZE + Npoll,KMAX_MID),x
    integer ::isrc
    integer ::k,n,ix,iix,dx,dy
    !how far diffusion should take place above lf_Nvert.
    ! KUP = 2 gives less than 0.001 differences in locfrac, except sometimes over sea, because
    !ship emission are higher up and need to come down to diminish locfrac
    integer, parameter :: KUP = 2
    if(DEBUGall .and. me==0)write(*,*)'start diff'

    call Code_timer(tim_before)
    xn_k = 0.0
    do k = 1,KMAX_MID
       do isrc=1,Nsources
          x=0.0
          do iix=1,lf_src(isrc)%Nsplit
             ix=lf_src(isrc)%ix(iix)
             !assumes mixing ratios units, but weight by mass
             x=x+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
          end do
          if(k>KMAX_MID-lf_Nvert)then ! lf zero above
             do n=lf_src(isrc)%start, lf_src(isrc)%end
                xn_k(n,k)=x*lf(n,i,j,k)
             enddo
          end if
          xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k) = x
      enddo
    enddo

    call vertdiffn(xn_k,LF_SRC_TOTSIZE+Npoll,1,KMAX_MID-lf_Nvert-KUP,EtaKz(i,j,1,1),dt_advec,ndiff) !DSKz

    do k = KMAX_MID-lf_Nvert+1,KMAX_MID
       do isrc=1,Nsources
          x =  1.0/(xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k)+1.E-30)
          do n=lf_src(isrc)%start, lf_src(isrc)%end
             lf(n,i,j,k) = xn_k(n,k)*x
          enddo
       enddo
    end do
    call Add_2timing(NTIMING-5,tim_after,tim_before,"lf: diffconv")
    if(DEBUGall .and. me==0)write(*,*)'end diff'

end subroutine lf_diff


  subroutine lf_conv(i,j,dt_conv)

    implicit none
    integer, intent(in) :: i,j
    real, intent(in) :: dt_conv
    real :: xn_k(LF_SRC_TOTSIZE + Npoll,KMAX_MID),x,xd
    integer ::isrc
    integer ::k,n,ix,iix,dx,dy

    if(DEBUGall .and. me==0)write(*,*)'start conv'
    call Code_timer(tim_before)
    xn_k = 0.0
    do k = 1,KMAX_MID
       do isrc=1,Nsources
          x=0.0
          do iix=1,lf_src(isrc)%Nsplit
             ix=lf_src(isrc)%ix(iix)
             !assumes mixing ratios units, but weight by mass
             x=x+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
          end do
          if(k>KMAX_MID-lf_Nvert)then ! lf zero above
             do n=lf_src(isrc)%start, lf_src(isrc)%end
                xn_k(n,k)=x*lf(n,i,j,k)
             enddo
          end if
          xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k) = x
      enddo
    enddo

    call convection_1d(xn_k,LF_SRC_TOTSIZE+Npoll,1,i,j,dt_conv)

    do k = KMAX_MID-lf_Nvert+1,KMAX_MID
       do isrc=1,Nsources
          x =  1.0/(xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k)+1.E-30)
          do n=lf_src(isrc)%start, lf_src(isrc)%end
             lf(n,i,j,k) = xn_k(n,k)*x
          enddo
       enddo
    end do
    call Add_2timing(NTIMING-5,tim_after,tim_before,"lf: diffconv")
    if(DEBUGall .and. me==0)write(*,*)'end conv'

end subroutine lf_conv

!make and use chemical and emission derivatives
subroutine lf_chem_emis_deriv(i,j,k,xn,xnew)
  real, intent(in) :: xn(NSPEC_TOT),xnew(NSPEC_TOT)
  integer, intent(in) :: i,j,k
  integer :: n, n0, ispec, nispec, isrc,isrc_emis, ic, is, ics, ideriv0, found
  real :: efac, xtot,totemis,emiss, xd
  integer :: isec, iem, iix, ix, iiix,iemis, ideriv , iem_deriv, isrc_deriv, n_sp

  if(k<KMAX_MID-lf_Nvert+1)return

  if (i<li0 .or.i>li1 .or.j<lj0.or.j>lj1)return !we avoid outer frame
  if(DEBUGall .and. me==0)write(*,*)'start chememis'
  if (.not. lf_fullchem) then
     !case with no chemistry for local fractions
     !Now species may be group of species (pm25...)

     if(k<max(KEMISTOP,KMAX_MID-lf_Nvert+1))then
!        if(nemis_primary>0)write(*,*)'WARNING nemis_primary not zero',nemis_primary, me,i,j,k
        if(N_lf_derivemis>0)write(*,*)'WARNING N_lf_derivemis not zero',N_lf_derivemis,me,i,j,k
!        nemis_primary = 0
        N_lf_derivemis = 0
        return
     end if
     !include emissions that are not created in chemical reactions, but only by emissions
     call Code_timer(tim_before)
     do isrc=1,Nsources
        if (lf_set%full_chem) cycle !included below
        iem = lf_src(isrc)%iem     !pick out only this emissions

        xtot=0.0
        totemis = 0.0 ! all emis to isrc species
        do iix=1,lf_src(isrc)%Nsplit
           iiix=lf_src(isrc)%ix(iix)
           !xtot=xtot+(xn_adv(iiix,i,j,k)*lf_src(isrc)%mw(iix))*(dA(k)+dB(k)*ps(i,j,1))/ATWAIR/GRAV
           xtot=xtot + xn_2d(iiix+NSPEC_SHL,k)*lf_src(isrc)%mw(iix)/dt_advec
           totemis = totemis +  rcemis(iiix+NSPEC_SHL,k)*lf_src(isrc)%mw(iix)
        end do

        if(totemis < 1e-20) cycle !no emissions here

        if(lf_src(isrc)%type=='country' .and. (Ncountry_lf+Ncountry_group_lf>0))then
           !first dilute lf because of emissions
           n0=lf_src(isrc)%start
           do is=1,Ncountrysectors_lf
              do ic=1,Ncountry_lf+Ncountry_group_lf
                 lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot)/(xtot+totemis+1.e-20)
                 n0=n0+1
              end do
           end do
           !add emissions that are tracked
           do iemis = 1, N_lf_derivemis
               if ( emis2icis(iemis) < 0 ) cycle ! relative type
               if ( emis2isrc(iemis) /= isrc ) cycle
               n0 = lf_src(isrc)%start + emis2icis(iemis)
               lf(n0,i,j,k)=lf(n0,i,j,k) + rcemis_lf(iemis,1)/(xtot+totemis+1.e-20)
               rcemis_lf(iemis,1) = 0.0
           end do
        else if(lf_src(isrc)%type=='relative')then
           !first dilute lf because of emissions
           do n0=lf_src(isrc)%start,lf_src(isrc)%end
              lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot)/(xtot+totemis+1.e-20)
           end do
           !add emissions only at centre of window for each sector
           n0 = lf_src(isrc)%start + (lf_src(isrc)%Npos - 1)/2 !"middle" point is dx=0 dy=0
           do iemis = 1, N_lf_derivemis
              if(emis2icis(iemis) >= 0 ) cycle !country type
              if(emis2isrc(iemis) /= isrc ) cycle
              lf(n0,i,j,k)=lf(n0,i,j,k) + rcemis_lf(iemis,1)/(xtot+totemis+1.e-20)
           end do
           !new notations, more compatible with "fullchem"
           if (k==KMAX_MID) then
              do iemis = 1, nemis_primary
                 if(emis2isrc_primary(iemis) /= isrc ) cycle
                 lf(n0,i,j,k)=lf(n0,i,j,k) + rcemis_lf_primary(iemis)/(xtot+totemis+1.e-20)
              end do
           end if
        else
           if(me==0)write(*,*)'LF type not recognized)'
           stop
        end if

     enddo
     do iemis = 1, N_lf_derivemis
        rcemis_lf(iemis,1) = 0.0 !inititalization for next gridcell. TODO: Needed here?
     end do
     N_lf_derivemis = 0
     if (k==KMAX_MID) then
        do iemis = 1, nemis_primary
           rcemis_lf_primary(iemis) = 0.0
        end do
        nemis_primary = 0
     end if
     call Add_2timing(NTIMING-4,tim_after,tim_before,"lf: emissions")
     return
  else
     !case using Jacobian of chemical reactions and emissions are considered part of chemical reactions
  call Code_timer(tim_before)

  !make derivatives
  do ispec = 1, NSPEC_chem_lf!no loop over short lived
     n_sp = lfspec2spec(ispec) !CM_Spec index
     do ideriv = 1, NSPEC_deriv_lf !loop over chemical derivatives without SOA derivatives
        if (xn(n_sp)>1000.0 .and. xnew(n_sp)>1000.0 .and. n_sp/=RO2POOL_ix .and. ideriv/=RO2POOL_ix-NSPEC_SHL.and. xn(lfspec2spec(ideriv))>1000.0) then !units are molecules/cm3, 1000  means almost no molecules
!        if (xn(n_sp)>1000.0 .and. xnew(n_sp)>1000.0 .and. n_sp/=RO2POOL_ix .and. ideriv/=RO2POOL_ix-NSPEC_SHL) then !units are molecules/cm3, 1000  means almost no molecules
           ! derivatives are normalized to concentration of each species, dx/dy *y/x
           !  y and dy must have same units
           !(xnew_lf(ideriv,ispec) - xnew(ispec))/(0.0001*xn(ideriv+NSPEC_SHL)) * xn(ideriv+NSPEC_SHL)/xnew(ispec) =
           xderiv(ideriv,ispec) = (xnew_lf(ideriv,n_sp) - xnew(n_sp))/((eps1-1.0)*(xnew(n_sp)))
           !remove numerical noise
           xd =  (xnew_lf(ideriv,n_sp) - xnew(n_sp))/((eps1-1.0)*xn(lfspec2spec(ideriv)))

           xderiv(ideriv,ispec) = min(xderiv(ideriv,ispec),10.0)
           xderiv(ideriv,ispec) = max(xderiv(ideriv,ispec),-10.0)

        else
           !Not that concentrations reaching or leaving 0.0 are not treated correctly
           !remove numerical noise
           xderiv(ideriv,ispec) = 0.0! should put diago to 1.0 ?
        end if
        !if(ispec>NSPEC_deriv_lf )xderiv(ideriv,ispec) = 0.0
     end do

     do ideriv = 1, N_lf_derivemis !loop over emission derivatives
        if (xnew(n_sp)>1000.0 .and. xnew_lf(ideriv + NSPEC_deriv_lf,n_sp)>1000.0 .and. k>=KEMISTOP.and. n_sp/=RO2POOL_ix ) then !units are molecules/cm3, 1000 means almost no molecules
           ! derivatives are normalized to concentration of each species, dx/dy *y/x
           !  y and dy must have same units
           !(xnew_lf(ideriv,ispec) - xnew(ispec))/(0.0001*xn(ideriv+NSPEC_SHL)) * xn(ideriv+NSPEC_SHL)/xnew(ispec) =
           ederiv(ideriv,ispec) = (xnew_lf(ideriv  + NSPEC_deriv_lf ,n_sp) - xnew(n_sp))/((eps1-1.0)*(xnew(n_sp)))

           !remove numerical noise
           ederiv(ideriv,ispec) = min(ederiv(ideriv,ispec),10.0)
           ederiv(ideriv,ispec) = max(ederiv(ideriv,ispec),-10.0)

        else
           !Note that concentrations reaching 0.0 are not treated correctly
           ederiv(ideriv,ispec) = 0.0
         end if
!       if(ispec>NSPEC_deriv_lf )ederiv(ideriv,ispec) = 0.0
     end do
  end do

  !remove RO2POOL which is not a proper species
  if (RO2POOL_ix>0) then
     xderiv(:,RO2POOL_ix-NSPEC_SHL) = 0.0
     ederiv(1:N_lf_derivemis,RO2POOL_ix-NSPEC_SHL) = 0.0
     xderiv(RO2POOL_ix-NSPEC_SHL,:) = 0.0
  end if

  ! emissions derivatives taken into account for sensibility to changes of emis *during* chemical process

  do isrc = 1, Nsources_chem !do not include NH4_f
     ix = lf_src(isrc)%ix(1) ! index of species in species_adv()
     n0 = lf_src(isrc)%start
     ispec = spec2lfspec(ix+NSPEC_SHL) !spec index in lf arrays
     if(ispec<0)cycle !not a lf_chem_deriv species
     !save the original lf values
     n = 1
     if (lf_src(isrc)%iem_lf == iem_lf_voc) n = Npos_lf + 1 !Same countries, but different nox/voc
     if (lf_src(isrc)%iem_lf == iem_lf_nh3) n = 2*Npos_lf + 1 !Same countries, but nh3
     if (lf_src(isrc)%iem_lf == iem_lf_sox) n = 3*Npos_lf + 1 !Same countries, but sox
     if (lf_set%EmisDer_all) n = 1

     !we must save the original lf, because we do not want to use and modify them at the same time
     if(ispec<=N_deriv_SOA_lf)then !does not include NH3, but includes SOA
        if(ispec<=NSPEC_deriv_lf)then
           do ics = n0,  n0 + Npos_lf - 1 !loop over sector and country contributions
              lf0_loc(n, ispec) = lf(ics,i,j,k)
              lf(ics,i,j,k) = 0.0 ! initialization
              n=n+1
           end do
        else
           do ics = n0,  n0 + Npos_lf - 1 !loop over sector and country contributions
              lf0SOA_loc(n, ispec-NSPEC_deriv_lf) = lf(ics,i,j,k)
              lf(ics,i,j,k) = 0.0 ! initialization
              n=n+1
           end do
        end if
     else
        !NH3
        !species that are not included in the derivatives (xderiv). Make diagonal element (1 if they do not change).
        if (xnew(ix+NSPEC_SHL)>1000.0)then
           xtot=xn_2d(ix+NSPEC_SHL,k)
           do ics = n0,  n0 + Npos_lf - 1
              lf(ics,i,j,k) = (lf(ics,i,j,k)*xtot)/xnew(ix+NSPEC_SHL) !Use as initialization
           end do
        else
           do ics = n0,  n0 + Npos_lf - 1
              lf(ics,i,j,k) = 0.0
           end do
        end if
     end if
     !include emission derivatives.
     do iemis = 1, N_lf_derivemis
        if(emis2iem(iemis) == lf_src(isrc)%iem_deriv .or. lf_set%EmisDer_all) then
           n0 =lf_src(isrc)%start + emis2icis(iemis)
           !contribution from emissions during this timestep.
           !NB: only one iemis per n0 can be included
           lf(n0,i,j,k) = lf(n0,i,j,k) + ederiv(iemis,ispec)
           if (lf_set%relative) then
              is =  emis2is(iemis)
              n0 = lf_src(isrc)%start + (Ncountry_lf + Ncountry_group_lf) * Ncountrysectors_lf !start of relative
              !relative: for sector is add all countries to the middle gridcell
              n0 = n0 + (is-1)*Ndiv_rel*Ndiv_rel + (Ndiv_rel*Ndiv_rel - 1)/2
              lf(n0,i,j,k) = lf(n0,i,j,k) + ederiv(iemis,ispec)
           end if
        end if
     end do
     !for all countries emitting, and each species and sector add to middle of "relative"
  end do

  ! concentrations have changed in the chemical reactions.
  ! The derivatives give the sensitivity to each species
  ! The local fraction are then propagated with the derivatives as weights
  lf_loc = matmul(lf0_loc, xderiv) !MAIN STEP

  !Note that if a species has emissions and no chemistry, xderiv(x,x) = x(t)/(x(t)+E(dt)) /= 1
  ! xderiv(x,x) = ((x*eps1+E) - (x+E))/((eps1-1)*(x+E)) = ((x*eps1 -x) +E -E)/(eps1-1) /(x+E) = x/(x+E)

  !add contribution from chemistry
  do isrc = 1, Nsources_chem! Also SO4, NO3_c
     ix = lf_src(isrc)%ix(1) ! index of species in species_adv()
     n_sp = spec2lfspec(ix+NSPEC_SHL)
     if(n_sp<0)cycle !not a lf_chem_deriv species
     n0 = lf_src(isrc)%start
     n = 1
     if (lf_src(isrc)%iem_lf == iem_lf_voc) n = Npos_lf + 1 !Same countries, but different nox/voc
     if (lf_src(isrc)%iem_lf == iem_lf_nh3) n = 2*Npos_lf + 1 !Same countries, but nh3
     if (lf_src(isrc)%iem_lf == iem_lf_sox) n = 3*Npos_lf + 1 !Same countries, but sox
     if (lf_set%EmisDer_all) n = 1
     !NB: here we get SIA; since d(SIA)/d(O3) /= 0
     do ics = n0,  n0 + Npos_lf - 1 !loop over sector and country contributions
        lf(ics,i,j,k) = lf(ics,i,j,k) + lf_loc(n, n_sp)
        n=n+1
     end do
  end do

  !add dSOA_i/dSOA_j contributions
  do ispec = 1, NSOA !loop over SOA species only
     n_sp = lfspec2spec(NSPEC_deriv_lf + ispec) !CM_Spec index
     do ideriv = 1, NSOA !loop over SOA derivatives only
        if (xn(n_sp)>1000.0 .and. xnew(n_sp)>1000.0 .and. n_sp/=RO2POOL_ix) then !units are molecules/cm3, 1000 means almost no molecules
           ! derivatives are normalized to concentration of each species, dx/dy *y/x
           !  y and dy must have same units
           !(xnew_lf(ideriv,ispec) - xnew(ispec))/(0.0001*xn(ideriv+NSPEC_SHL)) * xn(ideriv+NSPEC_SHL)/xnew(ispec) =
           xderivSOA(ideriv,ispec) = (xnew_lf(ideriv+NSPEC_deriv_lf+N_lf_derivemis ,n_sp) - xnew(n_sp))/((eps1-1.0)*(xnew(n_sp)))
           !remove numerical noise
           xderivSOA(ideriv,ispec) = min(xderivSOA(ideriv,ispec),10.0)
           xderivSOA(ideriv,ispec) = max(xderivSOA(ideriv,ispec),-10.0)
        else
           !Not that concentrations reaching or leaving 0.0 are not treated correctly
           !remove numerical noise
           xderivSOA(ideriv,ispec) = 0.0
        end if
     end do
     !emissions taken only once
  end do

  !lf0SOA_loc(n, ispec) contains original lf values
  lfSOA_loc = matmul(lf0SOA_loc, xderivSOA) !SOA block

  !add contribution from SOA to lf
  do ispec = NSPEC_deriv_lf + 1, NSPEC_deriv_lf + NSOA
     do isrc = (ispec-1) * Nfullchem_emis + 1, (ispec-1) * Nfullchem_emis + Nfullchem_emis
        ix = lf_src(isrc)%ix(1) ! index of species in species_adv()
        n0 = lf_src(isrc)%start
        n = 1
        if (lf_src(isrc)%iem_lf == iem_lf_voc) n = Npos_lf + 1 !Same countries, but different nox/voc
        if (lf_src(isrc)%iem_lf == iem_lf_nh3) n = 2*Npos_lf + 1 !Same countries, but nh3
        if (lf_src(isrc)%iem_lf == iem_lf_sox) n = 3*Npos_lf + 1 !Same countries, but sox
        if (lf_set%EmisDer_all) n = 1
        do ics = n0,  n0 + Npos_lf - 1 !loop over sector and country contributions
           lf(ics,i,j,k) = lf(ics,i,j,k) + lfSOA_loc(n, ispec-NSPEC_deriv_lf)
           n=n+1
        end do
     end do
  end do

  !add "primary" emissions that are tracked directly
  do ispec=1,3
     if (ispec==1)isrc = isrc_pm25
     if (ispec==2)isrc = isrc_pm25_new
     if (ispec==3)isrc = isrc_pmco
     xtot=0.0
     totemis = 0.0 ! all emis to isrc species
     do iix=1,lf_src(isrc)%Nsplit
        iiix=lf_src(isrc)%ix(iix)
        xtot=xtot + xn_2d(iiix+NSPEC_SHL,k)*lf_src(isrc)%mw(iix)/dt_advec
        totemis = totemis +  rcemis(iiix+NSPEC_SHL,k)*lf_src(isrc)%mw(iix)
     end do
     if(totemis < 1e-20) cycle !no emissions here

     !first dilute lf because of emissions
     do n0 = lf_src(isrc)%start, lf_src(isrc)%end
        lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot)/(xtot+totemis+1.e-20)
     end do

     !add emissions that are tracked
     do iemis = 1, nemis_primary
        if ( emis2isrc_primary(iemis) /= isrc)cycle
        n0 = lf_src(isrc)%start + emis2pos_primary(iemis)
        lf(n0,i,j,k)=lf(n0,i,j,k) + rcemis_lf_primary(iemis)/(xtot+totemis+1.e-20)
     end do
  end do

  do n = 1,ix_lf_max
     do iemis = 1, N_lf_derivemis
        rcemis_lf(iemis,n) = 0.0 !initialization for next gridcell
     end do
  end do
  do iemis = 1, nemis_primary
     rcemis_lf_primary(iemis) = 0.0 !initialization for next gridcell
  end do

  !NB: N_lf_derivemis must be reset at end of  lf_chem_pos
  !nemis_primary = 0
  !N_lf_derivemis = 0


  
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")
  if(DEBUGall .and. me==0)write(*,*)'end chememis'
return
  if(i_fdom(i)==106.and.j_fdom(j)==94 .and.k==kmax_mid )then
     do ispec = 1,NSOA
        n_sp = lfspec2spec(NSPEC_deriv_lf + ispec) !CM_Spec index
        write(*,66)species(n_sp)%name//' ',(xderivSOA(n,ispec),' ',n = 1, NSOA)
     end do

  end if
return
  if(i_fdom(i)==-108.and.j_fdom(j)==95 .and.k==kmax_mid )then
     write(*,*)'concentration after chem ',xnew(18)
     write(*,*)'prediction after 0.001 nox emis change',xnew(18)*(1.0+0.001*xderiv(NSPEC_deriv_lf+1,18))
     write(*,*)'lf O3 ',lf(3,i,j,k)
     write(*,*)'prediction after chem and 0.001 nox emis change',xnew(18)*(1+0.001*lf(3,i,j,k))
  end if
  return
  67 format(10x,50(A11))
  if(i_fdom(i)==58.and.j_fdom(j)==36 .and.k==kmax_mid )then
     68 format(10x,50E11.4)
     write(*,68)(xn(n+NSPEC_SHL),n = 1,10)
     write(*,68)(xnew(n+NSPEC_SHL),n = 1,10)
     write(*,67)(trim(species_adv(n)%name)//' ',n = 1,10),'  NOXEMIS'

     do ispec = NSPEC_deriv_lf+1,NSPEC_deriv_lf+NSOA
66      format(A10,50(F10.5,A1))
        write(*,66)species(ispec)%name//' ',(xderiv(n,ispec),' ',n = 1, 10)
     end do
        write(*,66)'sum             ',(sum(xderiv(n,1:NSPEC_SHL+NSPEC_deriv_lf)),' ',n = 1,10)

  end if

end if
end subroutine lf_chem_emis_deriv

subroutine lf_chem(i,j)
  !track through chemical reactions
  integer, intent(in) ::i,j
  real :: VOC,HO2,O3,NO,NO2
  real :: d_O3,d_NO,d_NO2,d_VOC, k1,k2,J_phot,invt,inv
  real :: SO4,SO2, d_SO2, d_SO4
  real :: xn_new, xn_age
  integer :: k, n, n_O3,n_NO,n_NO2,n_VOC,nsteps,nsteps1,nsteps2, isrc
  integer :: n_SO2,n_SO4,  n_EC_new, n_EC, n_age, n_new, iix, ix
  real :: k_OH, k_H2O2, k_O3
  real ::  d_age, ageing_rate(KCHEMTOP:KMAX_MID), lf_temp

  call Code_timer(tim_before)

  ageing_rate = rct(102,:)!EC_AGEING_RATE()
  if(DEBUGall .and. me==0)write(*,*)'start chem'

  if (isrc_pm25 > 0) then
     do isrc = 1, Nsources_nonew
        if (isrc_new(isrc) < 0) cycle
        isrc_pm25_new = isrc_new(isrc)
        !transform some of the pm25_new into age
        do k = KMAX_MID-lf_Nvert+1,KMAX_MID
           xn_new = 0.0
           xn_age = 0.0
           do iix=1,lf_src(isrc_pm25_new)%Nsplit
              ix=lf_src(isrc_pm25_new)%ix(iix)
              xn_new=xn_new+xn_adv(ix,i,j,k)*lf_src(isrc_pm25_new)%mw(iix)
           end do
           d_age = ageing_rate(k)*xn_new *dt_advec
           if (d_age > 1.0E-30) then ! very important. Wrong results without if, even far from the country source.
              do iix=1,lf_src(isrc)%Nsplit
                 ix=lf_src(isrc)%ix(iix)
                 xn_age=xn_age+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
              end do
              inv = 1.0/( xn_age + d_age )
              n_new = lf_src(isrc_pm25_new)%start
              do n_age=lf_src(isrc)%start, lf_src(isrc)%end
                 lf(n_age,i,j,k) = (lf(n_age,i,j,k)*xn_age + d_age*lf(n_new,i,j,k)) * inv
                 n_new = n_new + 1
              enddo
           end if
        enddo
     end do
  end if

  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")
  if(DEBUGall .and. me==0)write(*,*)'end chem'
end subroutine lf_chem

subroutine lf_sia_pre(i,j,k,deriv_iter) !called just before AerosolEquilib
  integer, intent(in) ::i,j,k,deriv_iter
  if (.not.USES%LocalFractions .or. k<KMAX_MID-lf_Nvert+1) return
  if (deriv_iter == 1) aero_error = .false. !initialize for each i,j,k and timestep
  !save concentrations, to see changes
  if(DEBUGall .and. me==0)write(*,*)'start lf_sia_pre'
  call Code_timer(tim_before)
  if( .not.lf_fullchem .and. deriv_iter>1) return;
  if (lf_fullchem) then
 !    lf_NO3(k) = xn_2d(NO3_ix,k)
 !    lf_HNO3(k) = xn_2d(HNO3_ix,k)
  end if
  if ( .not.lf_fullchem) return;
  if (isrc_NH4_f<0 .and. .not.lf_fullchem) return;
  if(.not.lf_fullchem)then
     lf_NH4(k) = xn_2d(NH4_f_ix,k)
     lf_NH3(k) = xn_2d(NH3_ix,k)
     return
  end if

  !include a perturbation to get sensibilities

  if(deriv_iter == 1) then
     !save original values
     xn_lf(1,0) = xn_2d(NH3_ix,k)
     xn_lf(2,0) = xn_2d(NH4_f_ix,k)
     xn_lf(3,0) = xn_2d(NO3_f_ix,k)
     xn_lf(4,0) = xn_2d(HNO3_ix,k)
     xn_lf(5,0) = xn_2d(SO4_ix,k)
    !perturb HNO3
     xn_2d(HNO3_ix,k) = xn_2d(HNO3_ix,k) * eps1_sia
  else if(deriv_iter == 2) then
     !perturb SO4
     xn_2d(SO4_ix,k) = xn_2d(SO4_ix,k) * eps1_sia
  else if(deriv_iter == 3) then
     !perturb NH3
     xn_2d(NH3_ix,k) = xn_2d(NH3_ix,k) * eps1_sia
  else if(deriv_iter == 4) then
     !perturb HNO3 other direction
     xn_2d(HNO3_ix,k) = xn_2d(HNO3_ix,k) * eps1_sia_m
  else if(deriv_iter == 5) then
     !perturb SO4 other direction
     xn_2d(SO4_ix,k) = xn_2d(SO4_ix,k) * eps1_sia_m
  else if(deriv_iter == 6) then
     !perturb NH3 other direction
     xn_2d(NH3_ix,k) = xn_2d(NH3_ix,k) * eps1_sia_m
  else if(deriv_iter == 7) then
     !base case
     !NB: must be the last, so that code can continue with base case results
  end if

  if(DEBUGall .and. me==0)write(*,*)'end lf_sia_pre'
call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

end subroutine lf_sia_pre

subroutine lf_sia_pos(i,j,k,deriv_iter,pmwater,errmark) !called just after AerosolEquilib
  integer, intent(in) ::i,j,k,deriv_iter,errmark
  integer, intent(in) :: pmwater !0-> do not make water. 1-> make water. 2-> make water but do not reset
  real :: d_NH4, d_NH3, NH4, NH3, inv
  integer :: n_NH3, n_NH4
  real :: d_NO3, d_HNO3, NO3, HNO3
  integer :: n_NO3, n_HNO3,n,ix,iix,isrc,d
  real :: xderiv(6,4), fac, xd, xd3
  real, parameter :: mol_limit = 4 !how many molecules in ouptut are allowed to change for one molecule change in input
  real, parameter :: xd_limit = 0.5 !must be <1
  real, parameter :: xn_limit = 0 !mol/cm3(?)
  real, parameter :: deriv_limit = 10000
  logical:: derivok, NH3derivok, HNO3derivok
  integer, save ::ecount=0

  if(.not.USES%LocalFractions .or. k<KMAX_MID-lf_Nvert+1)return
  if( errmark < 0) aero_error = .true. !if any of the scenario goes wrong we abandon
  
  if (isrc_NH4_f<0 .and. .not.lf_fullchem) return
  if( .not.lf_fullchem .and. deriv_iter>1) return
  call Code_timer(tim_before)
  if(DEBUGall .and. me==0)write(*,*)'start lf_sia_pos'
  if(.not.lf_fullchem)then
     NH3 = xn_2d(NH3_ix,k)
     NH4 = xn_2d(NH4_f_ix,k)
     if (isrc_NH4_f>0)then
        d_NH4 = NH4 - lf_NH4(k)
        d_NH3 = NH3 - lf_NH3(k)

        if(d_NH4>0.0 .and. d_NH3<0.0)then
           !NH3 has been transformed into NH4
           n_NH3 = lf_src(isrc_NH3)%start
           inv = 1.0/(NH4+d_NH4)
           do n_NH4=lf_src(isrc_NH4_f)%start, lf_src(isrc_NH4_f)%end
              lf(n_NH4,i,j,k) = (lf(n_NH4,i,j,k)*NH4 + d_NH4*lf(n_NH3,i,j,k)) * inv
              n_NH3 = n_NH3 + 1
           enddo
        else if(d_NH4<0.0 .and. d_NH3>0.0)then
           !NH4 has been transformed into NH3
           n_NH3 = lf_src(isrc_NH3)%start
           inv = 1.0/(NH3+d_NH3)
           do n_NH4=lf_src(isrc_NH4_f)%start, lf_src(isrc_NH4_f)%end
              lf(n_NH3,i,j,k) = (lf(n_NH3,i,j,k)*NH3 + d_NH3*lf(n_NH4,i,j,k)) * inv
              n_NH3 = n_NH3 + 1
           enddo
        else
           !N is not conserved or concentrations are constant
        end if
     end if
  else if (lf_fullchem) then
     !save results (could avoid for deriv_iter=4?)
     if (pmwater>0) then
        xn_lf(0,deriv_iter) = PM25_water_rh50(i,j)
     else
        xn_lf(1,deriv_iter) = xn_2d(NH3_ix,k)
        xn_lf(2,deriv_iter) = xn_2d(NH4_f_ix,k)
        xn_lf(3,deriv_iter) = xn_2d(NO3_f_ix,k)
        xn_lf(4,deriv_iter) = xn_2d(HNO3_ix,k)
        xn_lf(5,deriv_iter) = xn_2d(SO4_ix,k)
     end if
     if(deriv_iter<Nsia_deriv .and. pmwater<2)then
        !put back original values
        !if pmwater
        xn_2d(NH3_ix,k) = xn_lf(1,0)
        xn_2d(NH4_f_ix,k)=xn_lf(2,0)
        xn_2d(NO3_f_ix,k)=xn_lf(3,0)
        xn_2d(HNO3_ix,k)=xn_lf(4,0)
        xn_2d(SO4_ix,k)=xn_lf(5,0)
     end if
     if (deriv_iter >= Nsia_deriv) then
        !make derivative dependencies (and keep xn_2d results)
        !HNO3 ->deriv 1 and 4
        !SO4  ->deriv 2 and 5
        !NH3  ->deriv 3 and 6
        !base ->deriv Nsia_deriv
        !original values ->deriv 0
        if (pmwater > 0) then
           if (k/=kmax_mid) then
              write(*,*)"wrong k",k,kmax_mid
              stop
           end if
           !make water lf
           do n = 1, Npos_lf*Nfullchem_emis
              lf_PM25_water(n,i,j) = 0.0
           end do
           if(xn_lf(0,Nsia_deriv)>0.0001)then
              xderiv(1,1) = min(10.0,max(-10.0,(xn_lf(0,1)- xn_lf(0,Nsia_deriv))/((eps1_sia-1.0)*xn_lf(0,Nsia_deriv))))
              xderiv(2,1) = min(10.0,max(-10.0,(xn_lf(0,2)- xn_lf(0,Nsia_deriv))/((eps1_sia-1.0)*xn_lf(0,Nsia_deriv))))
              xderiv(3,1) = min(10.0,max(-10.0,(xn_lf(0,3)- xn_lf(0,Nsia_deriv))/((eps1_sia-1.0)*xn_lf(0,Nsia_deriv))))
              !filter out dependency when there are almost no pollutants:
              if(xn_lf(1,0)+xn_lf(2,0)<10000)xderiv(3,1) = 0.0
              if(xn_lf(5,0)<10000)xderiv(2,1) = 0.0
              if(xn_lf(3,0)+xn_lf(4,0)<10000)xderiv(1,1) = 0.0
              !loop over X in dX/dY, but also over Y for each emis!
              do isrc=1, Nsources !NB: each species will appear Nfullchem_emis times, but we keep only the first
                 if(lf_src(isrc)%iem_lf>1)cycle !Treat all iem_lf at once -> only one species left
                 fac = 1.0
                 if(lf_src(isrc)%species == 'NH3')then
                    ix=3
                 else if(lf_src(isrc)%species == 'NH4_f')then
                    ix=3 !make from NH3
                    fac =  xn_lf(2,0)/(1+xn_lf(1,0))
                 else if(lf_src(isrc)%species == 'NO3_f')then
                    ix=1 !make from HNO3
                    fac =  xn_lf(3,0)/(5e-5+xn_lf(4,0))
                 else if(lf_src(isrc)%species == 'HNO3')then
                    ix=1
                 else if(lf_src(isrc)%species == 'SO4')then
                    ix=2
                 else
                    cycle
                 end if
                 ! Note: 5species*Nfullchem_emis passes through here for different isrc
                 ! contributions from 5 species are added to same index in lf_PM25_water
                 fac = PM25_water_rh50(i,j) * fac * xderiv(ix,1)
                 do n = 1, Npos_lf*Nfullchem_emis
                    lf_PM25_water(n,i,j) = lf_PM25_water(n,i,j) + fac * lf(lf_src(isrc)%start+n-1,i,j,k)
                 end do

              end do
           end if
        else
           if(.not. aero_error) then
              !make the concentration derivatives and update lf
              
              !save original lf values, and filter out all lf values larger than lf_limit              
              do n = 0, Npos_lf*Nfullchem_emis-1
                 if(abs(lf(lf_src(isrc_NH3)%start+n,i,j,k))<lf_limit_sia)lf(lf_src(isrc_NH3)%start+n,i,j,k)=0
                 lf0_loc(n+1,1) = lf(lf_src(isrc_NH3)%start+n,i,j,k)
              end do
              do n = 0, Npos_lf*Nfullchem_emis-1
                 if(abs(lf(lf_src(isrc_NH4_f)%start+n,i,j,k))<lf_limit_sia)lf(lf_src(isrc_NH4_f)%start+n,i,j,k)=0
                 lf0_loc(n+1,2) = lf(lf_src(isrc_NH4_f)%start+n,i,j,k)
              end do
              do n = 0, Npos_lf*Nfullchem_emis-1
                 if(abs(lf(lf_src(isrc_NO3_f)%start+n,i,j,k))<lf_limit_sia)lf(lf_src(isrc_NO3_f)%start+n,i,j,k)=0
                 lf0_loc(n+1,3) = lf(lf_src(isrc_NO3_f)%start+n,i,j,k)
              end do
              do n = 0, Npos_lf*Nfullchem_emis-1
                 if(abs(lf(lf_src(isrc_HNO3)%start+n,i,j,k))<lf_limit_sia)lf(lf_src(isrc_HNO3)%start+n,i,j,k)=0
                 lf0_loc(n+1,4) = lf(lf_src(isrc_HNO3)%start+n,i,j,k)
              end do
              do n = 0, Npos_lf*Nfullchem_emis-1
                 if(abs(lf(lf_src(isrc_SO4)%start+n,i,j,k))<lf_limit_sia)lf(lf_src(isrc_SO4)%start+n,i,j,k)=0
                 lf0_loc(n+1,5) = lf(lf_src(isrc_SO4)%start+n,i,j,k)
              end do
              derivok = .true. !if regime is supiscious, we keep old lf (derivok = false)
              NH3derivok = .true. !if NH3 derivative is supiscious, we do not take them into account 
              HNO3derivok = .true. !if NH3 derivative is supiscious, we do not take them into account 
              
              !note that the number of N is only conserved within a given derivative scenario. Here we check only base
              xd=abs(xn_lf(1,Nsia_deriv)+xn_lf(2,Nsia_deriv)-xn_lf(1,0)-xn_lf(2,0))/(1e-3+abs(xn_lf(1,0)+xn_lf(2,0)))!NH3+NH4 conserved
              if(xd>1e-4)then
                 derivok = .false.
              end if
              xd=abs(xn_lf(3,Nsia_deriv)+xn_lf(4,Nsia_deriv)-xn_lf(3,0)-xn_lf(4,0))/(1e-3+abs(xn_lf(3,0)+xn_lf(4,0)))!HNO3+NO3 conserved
              if(xd>1e-4)then
                 derivok = .false.
              end if
              do ix=1,4 !loop over X in dX/dY  NB: SO4 assumed not changed in aero
                 if(abs(xn_lf(ix,Nsia_deriv) - xn_lf(ix,0))>xn_limit .and.xn_lf(ix,Nsia_deriv)>xn_limit .and. xn_lf(ix,0)>xn_limit)then !units molec/cm3
                    !HNO3 derivative
                    if(xn_lf(ix,1)>xn_limit .and. xn_lf(4,0)>xn_limit)then
                       !NB: min and max values MUST be allowed large, because if the concentrations
                       !in for example NH3 and NH4 are very different, then a transfer of molecules from one
                       !to the other will make a big difference in concentrations of one of them.
                       xd = (xn_lf(ix,1)- xn_lf(ix,Nsia_deriv))/(eps1_sia-1.0)
                       xderiv(1,ix) = min(deriv_limit,max(-deriv_limit,xd/(1e-3+xn_lf(ix,Nsia_deriv))))!HNO3
                       
                       xd3 = (xn_lf(ix,1+3)- xn_lf(ix,Nsia_deriv))/(eps1_sia_m-1.0)
                       xderiv(1+3,ix) = min(deriv_limit,max(-deriv_limit,xd3/(1e-3+xn_lf(ix,Nsia_deriv))))!HNO3

                       if(abs(xd/xn_lf(4,0))>mol_limit .and. abs(xd3/xn_lf(4,0))>mol_limit) then
                          !they are both unrealistic
                          derivok=.false.
                       else if (abs(xd/xn_lf(4,0))>mol_limit .or. abs(xd3/xn_lf(4,0))>mol_limit) then
                          !choose the one which seems ok
                          if(abs(xd/xn_lf(4,0))>mol_limit) xderiv(1,ix)=xderiv(1+3,ix)
                       else if((1E-16+abs(xderiv(1,ix)-xderiv(1+3,ix)))/(1E-9+abs(xderiv(1,ix))+abs(xderiv(1+3,ix)))<xd_limit)then
                          !Note that if xd_limit<1, we are sure that xderiv(1,ix)*xderiv(1+3,ix)>0
                          xderiv(1,ix)=(xderiv(1,ix)+xderiv(1+3,ix))/2
                        else
                           HNO3derivok=.false.                       
                          !other things tried:
                          !choose the one which change the fewest number of molecules (NB: cannot filter xderiv, because it must be symmetric for NH4/NH3)
                          !if (abs(xd3) < abs(xd)) xderiv(1,ix) = xderiv(1+3,ix)
                          !average
                          !xderiv(1,ix)=(xderiv(1,ix)+xderiv(1+3,ix))/2
                          !limit
                          !xderiv(1,ix) = min(deriv_limit_small,max(-deriv_limit_small,xderiv(1,ix)))
                       end if
                    else
                       HNO3derivok=.false.
                    end if

                    !SO4 derivative
                    if(xn_lf(ix,2)>xn_limit .and. xn_lf(5,0)>xn_limit)then
                       !Note: even if xd is small, the derivative can be too large, in case SO4
                       !concentration is much higher than ix concentration because xd compares the
                       !number of molecules.
                       !a small change in SO4 can at most move all HNO3 to NO3 or all NO3 to HNO3.
                       xd = (xn_lf(ix,2)- xn_lf(ix,Nsia_deriv))/(eps1_sia-1.0)!SO4
                       xderiv(2,ix) = min(deriv_limit,max(-deriv_limit,xd/(1e-3+xn_lf(ix,Nsia_deriv))))
                       
                       xd3 = (xn_lf(ix,2+3)- xn_lf(ix,Nsia_deriv))/(eps1_sia_m-1.0)!SO4
                       xderiv(2+3,ix) = min(deriv_limit,max(-deriv_limit,xd3/(1e-3+xn_lf(ix,Nsia_deriv))))
                       
                       if(abs(xd/xn_lf(5,0))>mol_limit .and. abs(xd3/xn_lf(5,0))>mol_limit) then
                          !they are both unrealistic
                          derivok=.false.
                       else if (abs(xd/xn_lf(5,0))>mol_limit .or. abs(xd3/xn_lf(5,0))>mol_limit) then
                          !choose the one which seems ok
                          if(abs(xd/xn_lf(5,0))>mol_limit) xderiv(2,ix)=xderiv(2+3,ix)
                       else if((1E-16+abs(xderiv(2,ix)-xderiv(2+3,ix)))/(1E-9+abs(xderiv(2,ix))+abs(xderiv(2+3,ix)))<xd_limit)then
                          xderiv(2,ix)=(xderiv(2,ix)+xderiv(2+3,ix))/2
                       else
!                          derivok=.false.
                          xderiv(2,ix) = 0.0 !no SO4 sensibility included now
                       end if
                    else
                       derivok=.false.
                    end if
                    
                    !NH3 derivative
                    if(xn_lf(ix,3)>xn_limit .and. xn_lf(1,0)>xn_limit)then
                       xd = (xn_lf(ix,3)- xn_lf(ix,Nsia_deriv))/(eps1_sia-1.0)!NH3
                       xderiv(3,ix) = min(deriv_limit,max(-deriv_limit,xd/(1e-3+xn_lf(ix,Nsia_deriv))))

                       xd3 = (xn_lf(ix,3+3)-xn_lf(ix,Nsia_deriv))/(eps1_sia_m-1.0)!NH3
                       xderiv(3+3,ix) = min(deriv_limit,max(-deriv_limit,xd3/(1e-3+xn_lf(ix,Nsia_deriv))))
                       
                       if(abs(xd/xn_lf(1,0))>mol_limit .and. abs(xd3/xn_lf(1,0))>mol_limit) then
                          !they are both unrealistic
                          derivok=.false.
                       else if (abs(xd/xn_lf(1,0))>mol_limit .or. abs(xd3/xn_lf(1,0))>mol_limit) then
                          !choose the one which seems ok
                          if(abs(xd/xn_lf(1,0))>mol_limit) xderiv(3,ix)=xderiv(3+3,ix)                          
                       else if((1E-16+abs(xderiv(3,ix)-xderiv(3+3,ix)))/(1E-9+abs(xderiv(3,ix))+abs(xderiv(3+3,ix)))<xd_limit)then
                          xderiv(3,ix)=(xderiv(3,ix)+xderiv(3+3,ix))/2
                       else
                         NH3derivok=.false.
                       end if
                    else
                       NH3derivok=.false.
                    end if

                    !NB: derivatives can be large when close to discontinuity (NH4/SO4=2)
                    !xderiv(1,ix) = min(10.0,max(-10.0,xderiv(1,ix)))
                    !xderiv(2,ix) = min(10.0,max(-10.0,xderiv(2,ix)))
                    !xderiv(3,ix) = min(10.0,max(-10.0,xderiv(3,ix)))

                    !xderiv for NO3_f = xderiv for HNO3, no mass factor when units are in molecules/ or mixing ratio

                    !delta_in(HNO3)=xn_lf(4,0)*eps1-xn_lf(4,0)
                    !give same result as same difference in mol/m3 as delta_in(NO3_f)=xn_lf(4,0)*eps1-xn_lf(4,0)
                    !However has to correct because xderiv is relative to concentrations and they are different for NO3_f and HNO3
                    !lf(lf_src(isrc_HNO3)%start+n,i,j,k) = (xn_lf(4,1)- xn_lf(4,4)) / (xn_lf(4,0)*eps1-xn_lf(4,0)) * (xn_lf(4,0)/xn_lf(4,4)) *lf0_loc(n+1,4) +&
                    !                                      (xn_lf(4,1)- xn_lf(4,4)) / (xn_lf(4,0)*eps1-xn_lf(4,0)) * (xn_lf(3,0)/xn_lf(4,4)) *lf0_loc(n+1,3) +&
                    !                                      (xn_lf(4,2)- xn_lf(4,4)) / (xn_lf(5,0)*eps1-xn_lf(5,0)) * (xn_lf(5,0)/xn_lf(4,4)) *lf0_loc(n+1,SO4)
                    !lf(lf_src(isrc_NO3_f)%start+n,i,j,k)= (xn_lf(3,1)- xn_lf(3,4)) / (xn_lf(4,0)*eps1-xn_lf(4,0)) * (xn_lf(4,0)/xn_lf(3,4)) *lf0_loc(n+1,4) &
                    !                                      (xn_lf(3,1)- xn_lf(3,4)) / (xn_lf(4,0)*eps1-xn_lf(4,0)) * (xn_lf(3,0)/xn_lf(3,4)) *lf0_loc(n+1,3)
                    !                                      (xn_lf(3,2)- xn_lf(3,4)) / (xn_lf(5,0)*eps1-xn_lf(5,0)) * (xn_lf(5,0)/xn_lf(3,4)) *lf0_loc(n+1,SO4)
                    !lf(lf_src(isrc_NH3)%start+n,i,j,k)  = (xn_lf(1,1)- xn_lf(1,4)) / (xn_lf(4,0)*eps1-xn_lf(4,0)) * (xn_lf(4,0)/xn_lf(1,4)) *lf0_loc(n+1,4) &
                    !                                      (xn_lf(1,1)- xn_lf(1,4)) / (xn_lf(4,0)*eps1-xn_lf(4,0)) * (xn_lf(3,0)/xn_lf(1,4)) *lf0_loc(n+1,3)
                    !                                      (xn_lf(1,2)- xn_lf(1,4)) / (xn_lf(5,0)*eps1-xn_lf(5,0)) * (xn_lf(5,0)/xn_lf(1,4)) *lf0_loc(n+1,SO4)
                    !xderiv for NH4_f = xderiv for NH3, no mass factor when units are in molecules/ or mixing ratio
                    !TODO: set in matrix multiplication form
                 else
                    derivok=.false.
                 end if
              end do

              if(derivok)then
                 do isrc=1,Nsources
                   if(lf_src(isrc)%iem_lf>1)cycle !Treat all iem_lf at once
                   if(lf_src(isrc)%species == 'NH3')then
                       ix=1
                       if(xn_lf(1,Nsia_deriv)<1000)then
                          !no NH3 , we set lf to zero
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc)%start+n-1,i,j,k) = 0.0
                          end do
                          cycle
                       else if (xn_lf(2,Nsia_deriv)<1000)then
                          !no NH4 f , we set all NH4_f lf to zero
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc_NH4_f)%start+n-1,i,j,k) = 0.0
                          end do
                         ! and compute NH3 from RDN conservation
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc)%start+n-1,i,j,k) = (lf0_loc(n,1)*xn_lf(1,0) + lf0_loc(n,2)*xn_lf(2,0) -lf(lf_src(isrc_NH4_f)%start+n-1,i,j,k)*xn_lf(2,Nsia_deriv)) /(xn_lf(1,Nsia_deriv))
                          end do
                          cycle
                       end if
                    else if(lf_src(isrc)%species == 'NH4_f')then                      
                       ix=2
                       !NB: assumes NH3 sources have indices smaller than NH4_f
                       !derive NH4 from conservation of RDN:
                       !lf(NH3)*xn(NH3)+NH4)*xn(NH4)=cst
                       !lf(NH4)=(lf0(NH4)*xn0(NH4)+lf0(NH3)*xn0(NH3)-lf(NH3)*xn(NH3))/xn(NH4)
                       if(xn_lf(2,Nsia_deriv)<1000)then
                          !no NH4 f , we set all lf to zero above
                       else
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc)%start+n-1,i,j,k) = (lf0_loc(n,1)*xn_lf(1,0) + lf0_loc(n,2)*xn_lf(2,0) -lf(lf_src(isrc_NH3)%start+n-1,i,j,k)*xn_lf(1,Nsia_deriv)) /(xn_lf(2,Nsia_deriv))
                          end do
                       end if
                       cycle
                    else if(lf_src(isrc)%species == 'NO3_f')then                       
                       ix=3
                       !NB: assumes HNO3 sources have index smaller NO3_f
                       !derive NO3_f from conservation of OXN:                       
                       if(xn_lf(3,Nsia_deriv)<1000)then
                          !no NO3_f , we set all lf to zero in HNO3 case under
                       else
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc)%start+n-1,i,j,k) = (lf0_loc(n,3)*xn_lf(3,0) + lf0_loc(n,4)*xn_lf(4,0) -lf(lf_src(isrc_HNO3)%start+n-1,i,j,k)*xn_lf(4,Nsia_deriv)) /(xn_lf(3,Nsia_deriv))
                          end do
                       end if
                       cycle
                    else if(lf_src(isrc)%species == 'HNO3')then
                       ix=4
                       if(xn_lf(4,Nsia_deriv)<1000)then
                          !no HNO3 , we set all lf to zero
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc)%start+n-1,i,j,k) = 0.0
                          end do
                          cycle
                       else if(xn_lf(3,Nsia_deriv)<1000)then
                          !no NO3_f , we set all lf to zero
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc)%start+n-1,i,j,k) = 0.0
                          end do
                          ! and put all lf into HNO3, using OXN conservation
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc)%start+n-1,i,j,k) = (lf0_loc(n,3)*xn_lf(3,0) + lf0_loc(n,4)*xn_lf(4,0) -lf(lf_src(isrc_NO3_f)%start+n-1,i,j,k)*xn_lf(3,Nsia_deriv)) /(xn_lf(4,Nsia_deriv))
                          end do
                          cycle
                       end if
                    else if(lf_src(isrc)%species == 'SO4')then
                       cycle !SO4 assumed not changed in aero
                       ix=5
                    else
                       cycle
                    end if

                    if(abs(xn_lf(ix,Nsia_deriv) - xn_lf(ix,0))>xn_limit .and.xn_lf(ix,Nsia_deriv)>xn_limit .and. xn_lf(ix,0)>xn_limit)then !units molec/cm3
                       !d(isrc)/dn = sum_i(xderiv(i,isrc)*di/dn) , xderiv(i,isrc)=d(isrc)/di .  i for xderiv and lf0_loc are not using the same mapping
                       !TODO: write as matrix multiplication
                       if(HNO3derivok .and. NH3derivok)then
                          do n = 1, Npos_lf*Nfullchem_emis
                             lf(lf_src(isrc)%start+n-1,i,j,k) = xderiv(1,ix) * lf0_loc(n,4) &
                                  + xderiv(1,ix) * xn_lf(3,0)/(1+xn_lf(4,0)) * lf0_loc(n,3) & !NB: xn_lf(3,0)/(1+xn_lf(4,0)) can be very large/small
                                  + xderiv(2,ix) * lf0_loc(n,5) &
                                  + xderiv(3,ix) * lf0_loc(n,1) &
                                  + xderiv(3,ix) * xn_lf(2,0)/(1+xn_lf(1,0)) * lf0_loc(n,2) !NB: xn_lf(2,0)/(1+xn_lf(1,0)) can be very large/small
                             
                          end do
                       else if (HNO3derivok .and. .not.NH3derivok) then
                          !we do not include NH3 and NH4 derivatives
                          if(ix==1)then
                             !NH3
                             do n = 1, Npos_lf*Nfullchem_emis
                                lf(lf_src(isrc)%start+n-1,i,j,k) = xderiv(1,ix) * lf0_loc(n,4) &
                                     + xderiv(1,ix) * xn_lf(3,0)/(1+xn_lf(4,0)) * lf0_loc(n,3) & !NB: xn_lf(3,0)/(1+xn_lf(4,0)) can be very large/small
                                     + xderiv(2,ix) * lf0_loc(n,5) &
                                     + lf0_loc(n,1)                             
                             end do
                             !NH4 treated with conservation principle above
                          end if
                          
                       else if (.not.HNO3derivok .and. NH3derivok) then
                          !we do not include HNO3 and NO3 derivatives
                          if(ix==4)then
                             !HNO3
                             do n = 1, Npos_lf*Nfullchem_emis
                                lf(lf_src(isrc)%start+n-1,i,j,k) = lf0_loc(n,4) &
                                     + xderiv(2,ix) * lf0_loc(n,5) &
                                     + xderiv(3,ix) * lf0_loc(n,1) &
                                     + xderiv(3,ix) * xn_lf(2,0)/(1+xn_lf(1,0)) * lf0_loc(n,2) !NB: xn_lf(2,0)/(1+xn_lf(1,0)) can be very large/small
                             end do                             
                             !NO3 treated with conservation principle above
                          end if
                          
                       end if
                       
                    else
                       !unchanged concentrations -> unchanged lf (happens if saturated, or other special cases)
                    end if
                    !filter out very small values 
                    do n = 1, Npos_lf*Nfullchem_emis
                       if(abs(lf(lf_src(isrc)%start+n-1,i,j,k))<lf_limit) lf(lf_src(isrc)%start+n-1,i,j,k) = 0.0
                    end do
                 end do
                 if(.false.)then
                 !check RDN conservation:
                 n=1 !nox reduction
!                 if(xn_2d(NH3_ix,k)>1000.and.xn_2d(NH4_f_ix,k)>1000 .and. ecount<10 .and. abs(lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k))+abs(lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k))>1e-4)then
                 if(xn_2d(NH3_ix,k)+xn_2d(NH4_f_ix,k)>2000 .and. ecount<10 )then
                    !                    if(me==251.and.i==6.and.j==2.and.k==12)then
                    xd= abs(lf0_loc(1,1)*xn_lf(1,0)+lf0_loc(1,2)*xn_lf(2,0)-(lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH3_ix,k)+lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH4_f_ix,k)))/&
                            (abs(lf0_loc(1,1)*xn_lf(1,0))+abs(lf0_loc(1,2)*xn_lf(2,0))+&
                            abs(lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH3_ix,k))+abs(lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH4_f_ix,k))+1)
                       if( abs(lf0_loc(1,1)*xn_lf(1,0)+lf0_loc(1,2)*xn_lf(2,0)-(lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH3_ix,k)+lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH4_f_ix,k)))/&
                            (abs(lf0_loc(1,1)*xn_lf(1,0))+abs(lf0_loc(1,2)*xn_lf(2,0))+&
                            abs(lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH3_ix,k))+abs(lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH4_f_ix,k))+1)>1e-3)then
                          write(*,*)me,i_fdom(i),j_fdom(j),k,'xd ',xd,lf0_loc(1,1)*xn_lf(1,0)+lf0_loc(1,2)*xn_lf(2,0)-(lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH3_ix,k)+lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH4_f_ix,k))
                          write(*,*)me,i,j,k,'err ',lf0_loc(1,1)*xn_lf(1,0),lf0_loc(1,2)*xn_lf(2,0),lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH3_ix,k),lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH4_f_ix,k)
                          write(*,*)me,lf_src(isrc_NH3)%start,lf_src(isrc_NH4_f)%start,'err2 ',lf0_loc(1,1)*xn_lf(1,0)+lf0_loc(1,2)*xn_lf(2,0),lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH3_ix,k)+lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k)*xn_2d(NH4_f_ix,k)
                       write(*,*)me,i,j,k,'lf val ',lf0_loc(1,1),lf0_loc(1,2),lf(lf_src(isrc_NH3)%start+n*Npos_lf-1,i,j,k),lf(lf_src(isrc_NH4_f)%start+n*Npos_lf-1,i,j,k)
                       write(*,*)me,i,j,k,'xn val ',xn_lf(1,0),xn_lf(2,0),xn_2d(NH3_ix,k),xn_2d(NH4_f_ix,k)
                       write(*,*)me,i,j,k,'xderiv HNO3 ',xderiv(1,1)*xn_lf(1,Nsia_deriv),xderiv(1,2)*xn_lf(2,Nsia_deriv),xderiv(1,1),xderiv(1,2)
                       write(*,*)me,i,j,k,'xderiv NH3 ',xn_lf(1,3)+xn_lf(2,3)-(xn_lf(1,Nsia_deriv)+xn_lf(2,Nsia_deriv)),(eps1_sia-1.0)*xn_lf(1,0),xderiv(3,1),xderiv(3,2)
                       ecount=ecount+1
                       if(ecount>100)stop
                    end if
                    end if
                 end if
              else
                 !unchanged concentrations -> unchanged lf (happens if unstable solution, saturated, or other special cases)
              end if
           
           end if
        end if
     end if
  end if
 
  if(DEBUGall .and. me==0)write(*,*)'end lf_sia_pos'
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

end subroutine lf_sia_pos

subroutine lf_aqu_pre(k,deriv_iter) !called just before aqrck rate calculations in setup_aqurates
  integer, intent(in) ::k,deriv_iter
  if (.not.USES%LocalFractions .or. k<KMAX_MID-lf_Nvert+1 .or. .not. lf_fullchem) return
  if(deriv_iter < 0)then
     !initialize only (not incloud)     
     AQRCK_lf(:,:,k) = 0.0
     fgasso2_lf(:,k) = 1.0
     return
  end if

  !save concentrations, to see changes
  if(DEBUGall .and. me==0)write(*,*)'start lf_sia_pre'
  call Code_timer(tim_before)
  if ( .not.lf_fullchem) return;

  !include a perturbation to get sensibilities

  if(deriv_iter == 1) then
     !save original values
     xn_lf(1,0) = xn_2d(NH3_ix,k)
     xn_lf(2,0) = xn_2d(NH4_f_ix,k)
     xn_lf(3,0) = xn_2d(NO3_f_ix,k)
     xn_lf(4,0) = xn_2d(HNO3_ix,k)
     xn_lf(5,0) = xn_2d(SO4_ix,k)
    !perturb HNO3 (same effect as changing NO3)
     xn_2d(HNO3_ix,k) = xn_2d(HNO3_ix,k) * eps1
  else if(deriv_iter == 2) then
     !perturb SO4
     xn_2d(SO4_ix,k) = xn_2d(SO4_ix,k) * eps1
  else if(deriv_iter == 3) then
     !perturb NH3
     xn_2d(NH3_ix,k) = xn_2d(NH3_ix,k) * eps1
  else if(deriv_iter == 4) then
     !perturb NH4
     xn_2d(NH4_f_ix,k) = xn_2d(NH4_f_ix,k) * eps1
  else if(deriv_iter == 5) then
     !base case
     !NB: must be the last, so that code can continue with base case results
  end if

  if(DEBUGall .and. me==0)write(*,*)'end lf_sia_pre'
call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

end subroutine lf_aqu_pre

subroutine lf_aqu_pos(k,deriv_iter,AQRCK)
  integer, intent(in) ::k,deriv_iter
  real,  intent(in) :: AQRCK(NAQUEOUS,KCHEMTOP:KMAX_MID)
  integer :: iter,n,i,n_sp
  real :: xd,xd_limit
  xd_limit=4.0
  if (.not.USES%LocalFractions .or. k<KMAX_MID-lf_Nvert+1 .or. .not. lf_fullchem) return
  call Code_timer(tim_before)
     !save results
  xn_lf(1,deriv_iter) = AQRCK(ICLRC1,K)
  xn_lf(2,deriv_iter) = aqrck(ICLRC2,k) 
  xn_lf(3,deriv_iter) = aqrck(ICLRC3,k) 
  xn_lf(4,deriv_iter) = aqrck(ICLHO2H2O2,k)
  xn_lf(5,deriv_iter) = FGAS(SO2_ix, k)

  if(deriv_iter<5)then
     !put back original values
     xn_2d(NH3_ix,k) = xn_lf(1,0)
     xn_2d(NH4_f_ix,k)=xn_lf(2,0)
     xn_2d(NO3_f_ix,k)=xn_lf(3,0)
     xn_2d(HNO3_ix,k)=xn_lf(4,0)
     xn_2d(SO4_ix,k)=xn_lf(5,0)
  end if
  if (deriv_iter >= 5) then
     !make derivative dependencies (and keep xn_2d results)
     !HNO3 ->deriv 1
     !SO4  ->deriv 2
     !NH3  ->deriv 3
     !NH4  ->deriv 4
     !base ->deriv 5
     do i=1,NAQUEOUS
        do n=1,NSPEC_deriv_lf+N_lf_derivemisMAX
           AQRCK_lf(n,i,K)=AQRCK(i,K)
        end do
     end do
     do n=1,NSPEC_deriv_lf+N_lf_derivemisMAX
        fgasso2_lf(n,k) = 1.0
     end do
     do iter = 1,4
        if(iter==1)n_sp = spec2lfspec(HNO3_ix)
        if(iter==2)n_sp = spec2lfspec(SO4_ix)
        if(iter==3)n_sp = spec2lfspec(NH3_ix)
        if(iter==4)n_sp = spec2lfspec(NH4_f_ix)

        xd=(xn_lf(1,iter)-xn_lf(1,5))/((eps1-1.0)*xn_lf(1,5))
        if(abs(xd)>xd_limit)xn_lf(1,iter)=xn_lf(1,5)
        aqrck_lf(n_sp,ICLRC1,K) =xn_lf(1,iter)
        
        xd=(xn_lf(2,iter)-xn_lf(2,5))/((eps1-1.0)*xn_lf(2,5))
        if(abs(xd)>xd_limit)xn_lf(2,iter)=xn_lf(2,5)
        aqrck_lf(n_sp,ICLRC2,K) =xn_lf(2,iter)
        
         xd=(xn_lf(3,iter)-xn_lf(3,5))/((eps1-1.0)*xn_lf(3,5))
        if(abs(xd)>xd_limit)xn_lf(3,iter)=xn_lf(3,5)
        aqrck_lf(n_sp,ICLRC3,K) =xn_lf(3,iter)
        
        xd=(xn_lf(4,iter)-xn_lf(4,5))/((eps1-1.0)*xn_lf(4,5))
        if(abs(xd)>xd_limit)xn_lf(4,iter)=xn_lf(4,5)
        aqrck_lf(n_sp,ICLHO2H2O2,k) = xn_lf(4,iter)
         
        if(xn_lf(5,5)>1e-5)then
           xd=xn_lf(5,iter)/(xn_lf(5,5))
        else
           xd=1.0
        end if
 
        fgasso2_lf(n_sp, k) = xd !xn_lf(5,iter)/xn_lf(5,5)
        
     end do
     
  end if
  call Code_timer(tim_before)
  if(DEBUGall .and. me==0)write(*,*)'end lf_aqu_pos'
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

end subroutine lf_aqu_pos

subroutine lf_SurfArea_pre(k,deriv_iter)
  integer, intent(in) :: k,deriv_iter
  if(.not.USES%LocalFractions .or. k<KMAX_MID-lf_Nvert+1)return
  if( .not.lf_fullchem ) return;
  if(DEBUGall .and. me==0)write(*,*)'start lf_SurfArea_pre'
  call Code_timer(tim_before)

  !include a perturbation to get sensibilities

  if(deriv_iter == 1) then
     !save original values
     xn_lf(1,0) = xn_2d(SO4_ix,k)
     xn_lf(2,0) = xn_2d(NO3_c_ix,k)
     xn_lf(3,0) = xn_2d(NO3_f_ix,k)
    !perturb SO4
     xn_2d(SO4_ix,k) = xn_2d(SO4_ix,k) * eps1
  else if(deriv_iter == 2) then
     !perturb NO3_c
     xn_2d(NO3_c_ix,k) = xn_2d(NO3_c_ix,k) * eps1
  else if(deriv_iter == 3) then
     !perturb NO3_f
     xn_2d(NO3_f_ix,k) = xn_2d(NO3_f_ix,k) * eps1
  else if(deriv_iter == 4) then
     !base case
     !NB: must be the last, so that code can continue with base case results
  end if
call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

   end subroutine lf_SurfArea_pre
subroutine lf_SurfArea_pos(S_m2m3,i,j,k,deriv_iter)
  integer, intent(in) :: i,j,k,deriv_iter
  real, intent(in) :: S_m2m3
  real rate(KCHEMTOP:KMAX_MID)
  if(.not.USES%LocalFractions .or. k<KMAX_MID-lf_Nvert+1)return
  if( .not.lf_fullchem ) return

  !save results
  xn_lf(1,deriv_iter) = S_m2m3 !NB: xn_lf(1,) is just used as an array, it is not a concentration!

  xn_lf(2,deriv_iter) = HYDROLYSISN2O5k(k) !NB: xn_lf(2,) is just used as an array, it is not a concentration!

  if(deriv_iter<4)then
     !put back original values
     xn_2d(SO4_ix,k)  = xn_lf(1,0)
     xn_2d(NO3_c_ix,k)= xn_lf(2,0)
     xn_2d(NO3_f_ix,k)= xn_lf(3,0)
  else
     !make derivative dependencies (and keep xn_2d base results)
     !SO4   -> deriv_iter = 1
     !NO3_c -> deriv_iter = 2
     !NO3_f -> deriv_iter = 3
     !base  -> deriv_iter = 4
!     if (xn_lf(1,4)>1e-10) then
!        rctAk_lf(spec2lfspec(SO4_ix),k) = (xn_lf(1,1) - xn_lf(1,4))/xn_lf(1,4)
!        rctAk_lf(spec2lfspec(NO3_c_ix),k) = (xn_lf(1,2) - xn_lf(1,4))/xn_lf(1,4)
!        rctAk_lf(spec2lfspec(NO3_f_ix),k) = (xn_lf(1,3) - xn_lf(1,4))/xn_lf(1,4)
!     else
!        !no dependency included
!        rctAk_lf(spec2lfspec(SO4_ix),k) = 0.0
!        rctAk_lf(spec2lfspec(NO3_c_ix),k) = 0.0
!        rctAk_lf(spec2lfspec(NO3_f_ix),k) = 0.0
!     end if

!     if (xn_lf(2,4)>1e-10) then
!        rctBk_lf(spec2lfspec(SO4_ix),k) = (xn_lf(2,1) - xn_lf(2,4))/xn_lf(2,4)
!        rctBk_lf(spec2lfspec(NO3_c_ix),k) = (xn_lf(2,2) - xn_lf(2,4))/xn_lf(2,4)
!        rctBk_lf(spec2lfspec(NO3_f_ix),k) = (xn_lf(2,3) - xn_lf(2,4))/xn_lf(2,4)
!        !should add NH4_f ? (or not contributing?)
!     else
!        !no dependency included
!        rctBk_lf(spec2lfspec(SO4_ix),k) = 0.0
!        rctBk_lf(spec2lfspec(NO3_c_ix),k) = 0.0
!        rctBk_lf(spec2lfspec(NO3_f_ix),k) = 0.0
!     end if
!           if(i==5.and.j==5 .and. k>=KMAX_MID-lf_Nvert+1 .and. me==253)then
!             if(xn_lf(2,4)>0.0)write(*,*)'LF ',me,k,xn_lf(2,1),xn_lf(2,4),rctBk_lf(spec2lfspec(SO4_ix),k)
!          end if
  end if
end subroutine lf_SurfArea_pos

subroutine  lf_drydep(i,j,DepLoss, fac)
  integer, intent(in) :: i,j
  real, intent(in) :: fac
  real, intent(in), dimension(NSPEC_ADV) :: DepLoss
  integer :: n,ix,iix,idep, idep0, isrc
  real :: ffac
  integer :: istart,iend
  idep0=0
  idep=0 !NB: order of indices must match the ones defined in lf_init
  call Code_timer(tim_before)

  if(DEBUGall .and. me==0)write(*,*)'start drydep'

  do isrc=1,Nsources
     if(.not. lf_src(isrc)%DryDep)cycle
     do iix=1,lf_src(isrc)%Nsplit
        ix=lf_src(isrc)%ix(iix)
        ffac = fac*1.e6*lf_src(isrc)%mw(iix) !(units ok?)
        istart = lf_src(isrc)%start
        iend = lf_src(isrc)%end
        !TODO: FFAC to revise for more general cases. or use %species_fac?
        if(.not. lf_fullchem)then
           if(isrc==isrc_SO4 .or. isrc==isrc_SO2 .or. lf_src(isrc)%species=="sox") ffac = ffac*32.0/64.0 !SO2->S
           if(isrc==isrc_NH3 .or. isrc==isrc_NH4_f .or. lf_src(isrc)%species=="nh3") ffac = ffac* 14.0/17.0!NH3->N
           if(isrc==isrc_NO .or. isrc==isrc_NO2 .or. lf_src(isrc)%species=="nox") ffac = ffac*14.0/46.0 !NO2->N
        end if
!        if( ix==SO4_ix - NSPEC_SHL ) then
!           !take directly local fractions from SO4 instead of sox
!           istart = lf_src(isrc_SO4)%start
!           iend= lf_src(isrc_SO4)%end
!        end if
!        if( ix==SO2_ix - NSPEC_SHL ) then
!           !take directly local fractions from SO2 instead of sox
!           istart = lf_src(isrc_SO2)%start
!           iend= lf_src(isrc_SO2)%end
!        end if

!        if( ix==NH4_f_ix - NSPEC_SHL ) then
!           istart = lf_src(isrc_NH4)%start
!           iend= lf_src(isrc_NH4)%end
!        end if
!        if( ix==NH3_ix - NSPEC_SHL ) then
!           istart = lf_src(isrc_NH3)%start
!           iend= lf_src(isrc_NH3)%end
!        end if

        idep=idep0
        do n = istart, iend
           idep=idep+1
           loc_frac_drydep(i,j,idep) = loc_frac_drydep(i,j,idep) + lf(n,i,j,KMAX_MID)*DepLoss(ix)*ffac
        enddo
        if(lf_src(isrc)%species=="nox" .and. iix==lf_src(isrc)%Nsplit)then
           !we add also depositions of NO3 and HNO3
           ix=NO3_ix - NSPEC_SHL
           idep=idep0
           do n = istart, iend
              idep=idep+1
              loc_frac_drydep(i,j,idep) = loc_frac_drydep(i,j,idep) + lf(n,i,j,KMAX_MID)*DepLoss(ix)*ffac
           enddo
           ix=HNO3_ix - NSPEC_SHL
           idep=idep0
           do n = istart, iend
              idep=idep+1
              loc_frac_drydep(i,j,idep) = loc_frac_drydep(i,j,idep) + lf(n,i,j,KMAX_MID)*DepLoss(ix)*ffac
           enddo
        end if
     enddo
     idep0 = idep0 + iend - istart + 1
  enddo
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")
  if(DEBUGall .and. me==0)write(*,*)'end drydep'
end subroutine lf_drydep

subroutine  lf_wetdep(iadv, i,j,k_in,loss, fac)
  !NB: loss is in units of molecules , without molecular weight
  integer, intent(in) :: iadv, i,j,k_in
  real, intent(in) :: loss, fac
  integer :: n,ix,iix,idep, idep0, isrc, k
  real :: ffac
  integer :: istart,iend
  idep0=0
  idep=0
  k = max(k_in,KMAX_MID-lf_Nvert+1) ! for scavenging above the lf window, we assume same fraction as highest level available
  call Code_timer(tim_before)

  do isrc=1,Nsources
     if(.not. lf_src(isrc)%WetDep) cycle
     istart = lf_src(isrc)%start
     iend = lf_src(isrc)%end
     do iix=1,lf_src(isrc)%Nsplit
        if(lf_src(isrc)%species=="nox")then
           !nochem nox case, we inlcude wet dep from other species and assume same lf value
           if(any(WDEP_OXN_GROUP==iadv+NSPEC_SHL) .or. iix==1)then
              continue
           else
              cycle
           end if
        else
           ix=lf_src(isrc)%ix(iix)
           if(ix /= iadv) cycle
        end if
        if (lf_fullchem) then
           ffac = fac*1.e6*lf_src(isrc)%mw(iix) !will be divided by mw(iix)(!) and multiplied by N, or S in lf_out
        else
           ffac = fac*1.e6
           if(isrc==isrc_SO2 .or. isrc==isrc_SO4 .or. lf_src(isrc)%species=="sox") ffac = ffac*32.0 !->S
           if(isrc==isrc_NH3 .or. isrc==isrc_NH4_f .or. lf_src(isrc)%species=="nh3") ffac = ffac* 14.0 !->N
           if(lf_src(isrc)%species=="nox") ffac = ffac*14.0 !->N
        end if

        idep=idep0
        do n = istart, iend
           idep=idep+1
           !wet scavenging from all levels are accumulated at surface.
           loc_frac_wetdep(i,j,idep) = loc_frac_wetdep(i,j,idep) + lf(n,i,j,k)*loss*ffac
        enddo
     enddo
     idep0 = idep0 + iend - istart + 1
  enddo
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

end subroutine lf_wetdep

subroutine  lf_PODnAOT(i,j, name, DepLoss, O3)
  integer, intent(in) :: i,j
  real, intent(in) :: DepLoss, O3
  character(len=*), intent(in) :: name
  integer :: iPODnAOT, iout,n,idep
  if (.not.lf_fullchem) return
  iPODnAOT = 0
  do iout = 1, Max_lf_out
     if (lf_spec_out(iout)%name == "NOTSET") exit
     if (index(lf_spec_out(iout)%name,"POD")>0 .or. index(lf_spec_out(iout)%name,"AOT")>0) then
        idep = NdryDep_lf + iPODnAOT * Nfullchem_emis*Npos_lf !NdryDep_lf does not include POD
        iPODnAOT = iPODnAOT + 1
        if(trim(lf_spec_out(iout)%name) == trim(name)) then
           do n = lf_src(isrc_O3)%start, lf_src(isrc_O3)%start + Nfullchem_emis*Npos_lf - 1
              idep = idep + 1
              if(DepLoss >1e-7)then
                 !Derivatives have no threshold!
                 loc_frac_drydep(i,j,idep) = loc_frac_drydep(i,j,idep) + lf(n,i,j,KMAX_MID)*O3
              end if
           end do
        end if
     end if
  end do

  return
end subroutine lf_PODnAOT
subroutine save_lf_emis(s,i,j,iem,isec,iland)
  ! save emission data
  ! save without splits and k
  ! save all (anthropogenic) emissions, even if not used by LF

  ! NB: This routine is called every hour

  real, intent(in) :: s
  integer, intent(in) :: i,j,iem,isec,iland
  integer :: n, ii, iqrc, isrc, k, ipoll,ic,iic,is,ig,emish_idx,split_idx

  if(DEBUGall .and. me==0)write(*,*)'start lf emis'
  call Code_timer(tim_before)
  if (s<1.E-20) return

  !include all countries
  !we store emissions per sector and country for each pollutant, WITHOUT SPLIT and WITHOUT vertical distribution
  do ic=1,nic(i,j)
    if(iland==ic2iland(i,j,ic)) exit
  end do
  if (ic>nic(i,j)) then
    nic(i,j) = nic(i,j) + 1
    ic2iland(i,j,ic) = iland
  end if
  emis_lf_cntry(i,j,ic,isec,iem) = emis_lf_cntry(i,j,ic,isec,iem) + s  !kg/m2/s
  
  call Add_2timing(NTIMING-4,tim_after,tim_before,"lf: emissions")

  if(DEBUGall .and. me==0)write(*,*)'end lf emis'
end subroutine save_lf_emis

subroutine lf_rcemis(i,j,k,eps)
  !makes the emission differences to be used for derivatives wrt emissions
  !We need to make rcemis_lf for each source (emitted species, country, sector) with emissions at this i,j

  ! The number of emitting countries at i,j is nic(i,j) (small)
  ! The emission without splits and k is emis_lf_cntry(i,j,ic,isec,iem) (ic runs over nic(i,j), i.e. corresponds to different countries for different gridcells)
  ! for "fullchem": same rcemis_lf for all advected_species (but different for source species), but lf countries dependent
  ! for "fullchem" rcemis_lf is splitted, otherwise summed with mw weights (and single species, _new_ old_ pm have not all splits)
  ! for "relative": country independent
  integer,intent(in) :: i,j,k
  real, intent(in) :: eps
  integer :: n, nn, n0, iem, iqrc, itot, ic, iic, ig, is, iland, isec, split_idx
  integer :: emish_idx, nemis, found,found_primary,iiix,isrc, nsectors_loop, iisec, isec_lf
  real :: ehlpcom0 = GRAV* 0.001*AVOG !0.001 = kg_to_g / m3_to_cm3
  real :: emiss,maskfac,y,z ! multiply emissions by
  real :: dtgrid !s*m2

  dtgrid = dt_advec * GRIDWIDTH_M * GRIDWIDTH_M !s*m2 . This routine is called every dt_advec
  call Code_timer(tim_before)
  if(DEBUGall .and. me==0)write(*,*)'start lf rcemis'
  if(k<max(KEMISTOP,KMAX_MID-lf_Nvert+1))return
  if (i<li0 .or.i>li1 .or.j<lj0.or.j>lj1)return !we avoid outer frame

  !rcemis_lf = 0.0 !init done at end of lf_chem_emis_deriv
  !rcemis_lf_primary = 0.0 !init done at end of lf_chem_emis_deriv
  !1) For now, we want to take derivative only from sector emissions, i.e. gridrcemis, and not fire, lightning, natural etc.
  nemis = 0
!  nemis_primary = 0
  N_lf_derivemis = 0 !number of distinct sources that have contributions in this gridcell
  if(k < KEMISTOP) return
  do iem = 1, NEMIS_File
    if (iem2Nipoll(iem) <= 0) cycle
    !for fullchem, we only treat nox , voc, nh3 and sox emissions
    if(lf_fullchem .and. EMIS_FILE(iem)/='nox' .and. EMIS_FILE(iem)/='voc' .and. EMIS_FILE(iem)/='nh3' .and. EMIS_FILE(iem)/='sox' .and. EMIS_FILE(iem)/='pm25' .and. EMIS_FILE(iem)/='pmco') cycle
    !we calculate the delta in emission for each country that has emissions
    if (.not.lf_fullchem) then
       do isrc=1,Nsources
          if (lf_src(isrc)%iem /= iem) cycle
          if (lf_src(isrc)%is_NATURAL) cycle !included in lf_rcemis_nat
          if(lf_src(isrc)%type=="relative")then
              !no lf_country involved, all countries are treated together
              if (lf_src(isrc)%nhour>0)then
                 !we add those emissions only to the sources with correct time index
                 if(lf_src(isrc)%time_ix /= lf_src(isrc)%nhour * (mod(current_date%hour,24)/lf_src(isrc)%nhour)) cycle
              end if
              n0 = 0
              !do not loop over countries, only sectors
              isec_lf = lf_src(isrc)%sector              
              found = 0
              do iisec=1,lf_nsector_map(isec_lf)                 
                 !isec is emep sector
                 isec = lf_sector_map(iisec,isec_lf) !isec will be counted as an isec_lf sector
                 do ic = 1,nic(i,j)
                    iland = ic2iland(i,j,ic)
                    if(emis_lf_cntry(i,j,ic,isec,iem)>1.E-20)then
                       if(found == 0) N_lf_derivemis = N_lf_derivemis +1  !should be increased at most once for each isrc
                       found = 1
                       emish_idx = SECTORS(isec)%height
                       split_idx = SECTORS(isec)%split
                       do n = 1, lf_src(isrc)%Nsplit
                          itot=lf_src(isrc)%ix(n)
                          iqrc = itot2iqrc(itot+NSPEC_SHL)
                          if (iqrc<0) cycle
                          rcemis_lf(N_lf_derivemis,1) = rcemis_lf(N_lf_derivemis,1) +&
                                emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)*ehlpcom0&
                               *emisfrac(iqrc,split_idx,iland)*emis_masscorr(iqrc)&
                               *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))  * lf_src(isrc)%mw(n)
                       end do
                    end if
                 end do
              end do

              if (found == 1) then
                 !emissions found here
                 emis2icis(N_lf_derivemis) = -1 ! only to show it is relative
                 emis2isrc(N_lf_derivemis) = isrc
              end if

           else if(lf_src(isrc)%type=="country")then
             do ic=1,nic(i,j)
                 iland = ic2iland(i,j,ic)
                 !only pick one country at a time
                 do iic=1,Ncountry_lf + Ncountry_group_lf !this loop could be avoided if necessary, by defining iland2iic?
                    maskfac = 1.0
                    if(iic<=Ncountry_lf)then
                      if (iic <= Ncountry_mask_lf) then
                        if (iic2ilf_countrymask(iic)>0) then
                           !Cell fraction of grid mask defined EmisMaskIndex2Name
                           maskfac = (1.0-EmisMaskValues(i,j,iic2ilf_countrymask(iic))) !count only what is covered by the mask
                           if (maskfac<1e-6) cycle
                        else if (country_mask_val(iic)>-99999) then
                           if(EmisMaskIntVal(i,j) /= country_mask_val(iic)) cycle !remove all parts which are not covered by mask
                        end if
                     else if(country_ix_list(iic)/=iland)then
                       cycle !only pick out one country
                     end if
                   else
                     !defined as a group of countries
                     found = 0
                     do ig = 1, MAX_lf_country_group_size
                       if (lf_country%group(iic-Ncountry_lf)%list(ig) == Country(iland)%code &
                            .or. lf_country%group(iic-Ncountry_lf)%name == 'ALL') then
                         found = 1
                         exit
                       end if
                       if (lf_country%group(iic-Ncountry_lf)%list(ig) == 'NOTSET') exit
                     end do
                     if (found == 0) cycle
                   end if
                   do is=1,Ncountrysectors_lf
                       isec_lf=lf_country%sector_list(is) !isec_lf are the sectors included in lf. Can have "holes".Can be pseudosector larger than NSECTOR for groups
                       found = 0 !flag to show if nemis already increased
                       do iisec=1, lf_nsector_map(isec_lf)
                          !isec is emep sector
                          isec = lf_sector_map(iisec,isec_lf) !isec will be counted as an isec_lf sector
                          if(emis_lf_cntry(i,j,ic,isec,iem)>1.E-20)then

                             if(found == 0) then
                                !first see if this emis2icis and emis2isrc already exists (for groups and masks type "countries"):
                                do n=1,N_lf_derivemis
                                   if(emis2icis(n)==(is-1)*(Ncountry_lf + Ncountry_group_lf) + (iic-1) .and. emis2isrc(n) == isrc)then
                                      ! add to this instead
                                      nemis = n
                                      found = 1
                                      exit
                                   end if
                                end do
                             end if
                             if(found == 0)  then
                                N_lf_derivemis = N_lf_derivemis + 1 !should be increased at most once if sector=0
                                nemis = N_lf_derivemis
                             end if

                             found = 1
                             emish_idx = SECTORS(isec)%height
                             split_idx = SECTORS(isec)%split
                             do n = 1, lf_src(isrc)%Nsplit
                                itot=lf_src(isrc)%ix(n)
                                iqrc = itot2iqrc(itot+NSPEC_SHL)
                                if (iqrc<0) cycle
                                rcemis_lf(nemis,1) = rcemis_lf(nemis,1) +&
                               maskfac* emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)*ehlpcom0&
                               *emisfrac(iqrc,split_idx,iland)*emis_masscorr(iqrc)&
                               *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))  * lf_src(isrc)%mw(n)
                             end do
                          end if
                       end do
                       if (found == 1) then
                         !emissions found here
                         emis2icis(nemis) = (is-1)*(Ncountry_lf + Ncountry_group_lf) + (iic-1)
                         emis2isrc(nemis) = isrc
                       end if
                    end do
                 end do
              end do

           end if
        end do

     else

        !fullchem case
        !TODO : merge with case not fullchem? difference only lf_src(isrc)%mw(n) and split summation and emis2 iem/isrc?
         do ic=1,nic(i,j)
 
          !iland is the country index of the emission treated
          iland = ic2iland(i,j,ic)
          !iic is the lf country index of the source
          do iic=1,Ncountry_lf + Ncountry_group_lf !this loop could be avoided if necessary, by defining iland2iic?
             maskfac = 1.0
             if(iic<=Ncountry_lf)then
                if (iic <= Ncountry_mask_lf) then
                   if (iic2ilf_countrymask(iic)>0) then
                      !Cell fraction of grid mask defined EmisMaskIndex2Name
                      maskfac = (1.0-EmisMaskValues(i,j,iic2ilf_countrymask(iic))) !count only what is covered by the mask
                      if (maskfac<1e-6) cycle
                   else if (country_mask_val(iic)>-99999) then
                      if(EmisMaskIntVal(i,j) /= country_mask_val(iic)) cycle !remove all parts which are not covered by mask
                   end if
                else if(country_ix_list(iic)/=iland)then
                   cycle !only pick out one country
                end if
             else
              !defined as a group of countries
              found = 0
              do ig = 1, MAX_lf_country_group_size
                if (lf_country%group(iic-Ncountry_lf)%list(ig) == Country(iland)%code &
                     .or. lf_country%group(iic-Ncountry_lf)%name == 'ALL') then
                  found = 1
                  exit
               end if
               if (lf_country%group(iic-Ncountry_lf)%list(ig) == 'NOTSET') exit
              end do
              if (found == 0) cycle
            end if
            do is=1,Ncountrysectors_lf
              isec_lf=lf_country%sector_list(is) !isec_lf are the sectors included in lf. Can have "holes". Can be pseudosector larger than NSECTOR for groups
              found = 0 !flag to show if nemis already increased
              found_primary = 0 !flag to show if nemis_primary already increased
              do iisec=1, lf_nsector_map(isec_lf)
                isec = lf_sector_map(iisec,isec_lf) !isec will be counted as an isec_lf sector
                emish_idx = SECTORS(isec)%height

                !sum all emissions per sector, country, emis_species. Only for one k value!
                if (k == KMAX_MID)then
                   totrcemis_lf(is,iic,iem) = totrcemis_lf(is,iic,iem) + emis_lf_cntry(i,j,ic,isec,iem) * dtgrid * xmd(i,j)
                end if
                
                if(emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)>1.E-20)then

                   if(EMIS_FILE(iem)=='pm25' .or. EMIS_FILE(iem)=='pmco')then
                      !special treated as primary, do not compute derivatives
                      if(EMIS_FILE(iem)=='pm25')then
                         if(found_primary == 0) then
                            nemis_primary = nemis_primary + 2 !pm25 and pm25_new
                            found_primary = 1
                         end if
                         isrc=isrc_pm25 !treated with index "nemis_primary-1"
                         emis2isrc_primary(nemis_primary-1) = isrc
                         emis2pos_primary(nemis_primary-1) = (is-1)*(Ncountry_lf + Ncountry_group_lf) + (iic-1)
                         split_idx = SECTORS(isec)%split
                         do n = 1, lf_src(isrc)%Nsplit
                            itot=lf_src(isrc)%ix(n)
                            iqrc = itot2iqrc(itot+NSPEC_SHL)
                            if (iqrc<0) cycle
                            rcemis_lf_primary(nemis_primary-1) = rcemis_lf_primary(nemis_primary-1) +&
                                 maskfac* emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)*ehlpcom0&
                                 *emisfrac(iqrc,split_idx,iland)*emis_masscorr(iqrc)&
                                 *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1)) * lf_src(isrc)%mw(n)
                         end do

                         isrc=isrc_pm25_new !treated with index "nemis_primary"

                      else
                         isrc=isrc_pmco
                         if(found_primary == 0) then
                            nemis_primary = nemis_primary + 1
                            found_primary = 1
                         end if
                      end if
                      emis2isrc_primary(nemis_primary) = isrc
                      emis2pos_primary(nemis_primary) =  (is-1)*(Ncountry_lf + Ncountry_group_lf) + (iic-1)
                      emish_idx = SECTORS(isec)%height
                      split_idx = SECTORS(isec)%split
                      do n = 1, lf_src(isrc)%Nsplit
                         itot=lf_src(isrc)%ix(n)
                         iqrc = itot2iqrc(itot+NSPEC_SHL)
                         if (iqrc<0) cycle
                         rcemis_lf_primary(nemis_primary) = rcemis_lf_primary(nemis_primary) +&
                              maskfac* emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)*ehlpcom0&
                              *emisfrac(iqrc,split_idx,iland)*emis_masscorr(iqrc)&
                              *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))  * lf_src(isrc)%mw(n)
                      end do
                   else

                      if(found == 0) then
                         !first see if this emis2icis and emis2iem already exists (for groups and masks type "countries"):
                         !because only one derivative for each source is included ("lf(n0,i,j,k) = ederiv(iemis,ix)" without +=)
                         do n=1,N_lf_derivemis
                            !check that emis2isrc(n) is not required for pm25 here
                            if(emis2icis(n)==(is-1)*(Ncountry_lf + Ncountry_group_lf) + (iic-1).and. emis2iem(n) == iem)then
                               ! add to this instead
                               nemis = n
                               found = 1
                               exit
                            end if
                         end do
                      end if

                      if(found == 0) then
                         N_lf_derivemis = N_lf_derivemis + 1 !should be increased at most once if sector=0
                         nemis = N_lf_derivemis
                         if(nemis>N_lf_derivemisMAX)write(*,*)'TOOLARGE NEMIS'
                      end if
                      found = 1
                      emish_idx = SECTORS(isec)%height
                      split_idx = SECTORS(isec)%split
                      do n = 1, emis_nsplit(iem)
                         iqrc = sum(emis_nsplit(1:iem-1)) + n
                         itot = iqrc2itot(iqrc)
                         if(itot<=ix_lf_max .or. .not.lf_fullchem)then
                            !note: eps < 0
                            rcemis_lf(nemis,itot) = rcemis_lf(nemis,itot) +&
                                 maskfac * eps * emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)*ehlpcom0&
                                 *emisfrac(iqrc,split_idx,iland)*emis_masscorr(iqrc)&
                                 *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))
                        end if
                     end do
                   end if
                end if
              end do
              if (found == 1) then
                 !emissions found here
                 if(nemis<=0)then
                    write(*,*)nemis,me,i,j,k,found,nemis_primary,iic,is
                    call StopAll('error')
                 end if
                 emis2icis(nemis) = (is-1)*(Ncountry_lf + Ncountry_group_lf) + (iic-1)
                 emis2is(nemis) = is !is=1 corresponds to first sector
                 emis2iem(nemis) = iem
              end if
           end do
        end do
     end do
     !add NAT emissions

     if(k==KMAX_MID .and. EMIS_FILE(iem)=="nox")then
        ! add surface emissions
        do n=1, Nemis_surf
           N_lf_derivemis = N_lf_derivemis + 1
           nemis = N_lf_derivemis
           do nn = 1, emis2nspec_surf(n)
              rcemis_lf(nemis,emis2spec_surf(n,nn)) = rcemis_lf_surf(n, nn) * eps
              rcemis_lf_surf(n, nn) = 0.0 !reset after use
           end do
           emis2nspec_surf(n) = 0 !reset after use
           iic = emis2iic_surf(n)
           is = 1
           emis2icis(nemis) = (is-1) * (Ncountry_lf + Ncountry_group_lf) + (iic-1)
           emis2iem(nemis) = iem
        end do
        Nemis_surf = 0 !reset after use
     end if
  end if
  end do

  if(N_lf_derivemis>N_lf_derivemisMAX)then
     write(*,*)me,i,j,k,N_lf_derivemis,N_lf_derivemisMAX
     call StopAll("too many sectors*country in one gridcell. Increase N_lf_derivemisMAX")
  end if

  call Add_2timing(NTIMING-4,tim_after,tim_before,"lf: emissions")
  if(DEBUGall .and. me==0)write(*,*)'end lf rcemis'
  end subroutine lf_rcemis

  subroutine lf_rcemis_nat(species_ix, rcemis, i, j, icountry_in)
    !only surface emissions implemented for now
    !so far only called for DMS, in Setup_1d_mod.f90

    integer, intent(in) :: species_ix, i, j !species_ix, index of species as defined in array "species" in CM_ChemSpecs_mod.f90 (not only advected species)
!    integer, optional, intent(in) :: k_in
    integer, optional, intent(in) :: icountry_in
    real, intent(in) :: rcemis !value of emissions to track
    integer :: isrc, n, ipos, stride, k
    integer ::  icountry !internal index for "country" (ix_BVOC or ix_DMS). Not used if not lf_fullchem

    if (rcemis < 1e-20) return
!    if(present(k_in)) then
!       k = k_in
!    else
       k = KMAX_MID
!    end if
     if (k < KMAX_MID-lf_Nvert+1) return

    if(present(icountry_in)) then
       icountry = icountry_in
    else
       icountry = 0
    end if
    if (i<li0 .or.i>li1 .or.j<lj0.or.j>lj1)return !we avoid outer frame
    if (lf_fullchem) then
       call CheckStop(icountry<=0, "country index for natural emissions not defined")
       if(Nemis_surf>0)then
          if(emis2iic_surf(Nemis_surf) == icountry) then
             !if same source, just add here
          else
             !new source
             Nemis_surf = Nemis_surf + 1
          end if
       else
          Nemis_surf = 1
       end if
       emis2nspec_surf(Nemis_surf) = emis2nspec_surf(Nemis_surf) + 1 !allow for split into several species
       emis2iic_surf(Nemis_surf) =  icountry
       emis2spec_surf(Nemis_surf, emis2nspec_surf(Nemis_surf)) = species_ix
       rcemis_lf_surf(Nemis_surf, emis2nspec_surf(Nemis_surf)) = rcemis
    else
       if (rcemis < 1e-20) return
       do isrc=1,Nsources
          if (lf_src(isrc)%species_ix /= species_ix) cycle
          if (.not. lf_src(isrc)%is_NATURAL) then
             call StopAll('Only natural emissions implemented in subroutine lf_rcemis_nat')
          end if
          if (lf_src(isrc)%nhour>0)then
             !we add those emissions only to the sources with correct time index
             if(lf_src(isrc)%time_ix /= lf_src(isrc)%nhour * (mod(current_date%hour,24)/lf_src(isrc)%nhour)) cycle
          end if
          nemis_primary = nemis_primary + 1
          emis2isrc_primary(nemis_primary) = isrc
          rcemis_lf_primary(nemis_primary) = rcemis_lf_primary(nemis_primary) + rcemis
       end do
    end if
  end subroutine lf_rcemis_nat

  subroutine addsource(species_name)
    character(len=*), intent(in) :: species_name
    integer :: i, ix, iem, isrc

    Nlf_species = Nlf_species + 1 !internal index for LF
    do i = 1, Nfullchem_emis
       Nsources = Nsources + 1
       lf_src(Nsources)%type = 'country'
       ix = find_index(trim(species_name), species_adv(:)%name, any_case=.true.)
       if (ix<=0) then
          iem=find_index(trim(species_name) ,EMIS_FILE(1:NEMIS_FILE))
          if(iem>0)then
             isrc = Nsources
             lf_src(isrc)%species = trim(species_name) !primary emitted species
             lf_src(isrc)%iem = iem
             return !NB: no "derivatives", only one source per country-sector
          else
             write(*,*)trim(species_name),' not found '
             stop
          end if
       end if
       lf_src(Nsources)%species = species_adv(ix)%name
       if(me==0 .and. DEBUG)write(*,*)lf_src(Nsources)%species
       if (index(lf_src(Nsources)%species,"ASO")>0) lf_src(Nsources)%is_ASOA = .true.
       lfspec2spec(Nlf_species) = ix + NSPEC_SHL
       spec2lfspec(ix + NSPEC_SHL) = Nlf_species
       lf_set%full_chem = .true.
       lf_src(Nsources)%make_fracsum = lf_src(1)%make_fracsum
       !iem_deriv differs from iem, because it relates not to the source species, but the derivative.
       !iem_deriv can be for voc, while the species is no2
       if (i==1) then
          if (.not.lf_set%EmisDer_all) lf_src(Nsources)%iem_deriv = find_index('nox' ,EMIS_FILE(1:NEMIS_FILE))
          if (.not.lf_set%EmisDer_all) lf_src(Nsources)%iem_lf = iem_lf_nox
          if (lf_set%EmisDer_all) lf_src(Nsources)%iem_lf = 1 ! so that (lf_src(isrc)%iem_lf-1)*Npos_lf = 0
          if (lf_src(Nsources)%species == 'O3') isrc_O3 = Nsources
          if (lf_src(Nsources)%species == 'NO2') isrc_NO2 = Nsources
          if (lf_src(Nsources)%species == 'NH4_f') isrc_NH4_f = Nsources
          if (lf_src(Nsources)%species == 'NH3') isrc_NH3 = Nsources
          if (lf_src(Nsources)%species == 'NO3') isrc_NO3 = Nsources
          if (lf_src(Nsources)%species == 'NO3_f') isrc_NO3_f = Nsources
          if (lf_src(Nsources)%species == 'HNO3') isrc_HNO3 = Nsources
          if (lf_src(Nsources)%species == 'SO4') isrc_SO4 = Nsources
          if (lf_src(Nsources)%species == 'SO2') isrc_SO2 = Nsources
          if (lf_src(Nsources)%species == 'NO3_c') isrc_NO3_c = Nsources
          if (lf_src(Nsources)%species == 'NH3') isrc_NH3 = Nsources
       else if (i==2) then
          lf_src(Nsources)%iem_lf = iem_lf_voc
          lf_src(Nsources)%iem_deriv = find_index('voc' ,EMIS_FILE(1:NEMIS_FILE))
       else if (i==3) then
          lf_src(Nsources)%iem_lf = iem_lf_nh3
          lf_src(Nsources)%iem_deriv = find_index('nh3' ,EMIS_FILE(1:NEMIS_FILE))
       else if (i==4) then
          lf_src(Nsources)%iem_lf = iem_lf_sox
          lf_src(Nsources)%iem_deriv = find_index('sox' ,EMIS_FILE(1:NEMIS_FILE))
       end if

    end do
  end subroutine addsource

  subroutine lf_chem_pre(i,j,k,dt,Nd)
    !The order of the "scenarios" aka derivatives is:
    !species until NO3_f | emissions derivatives | SOA
    !NB:  Dchem_lf, spec2lfspec and lfspec2specdoes not include emissions derivative indices!
    integer,intent(in) :: i,j,k
    integer,intent(inout) :: Nd
    real, intent(in) :: dt
    integer :: n,i_lf, n_sp
    if (k > KMAX_MID-lf_Nvert) then
       !make rcemis_lf and N_lf_derivemis
       call lf_rcemis(i,j,k,eps1-1.0)
    else
       return
    end if
    if (lf_fullchem .and. k > KMAX_MID-lf_Nvert) then
       !compute chemistry with small changes in input concentrations

       Nd = NSPEC_deriv_lf + N_lf_derivemis !shorter. NB: does not include SOA

       !Careful sometimes the only difference between xnew and xnew_lf are from initial values (Dchem_lf and/or xn_shl_lf)
       !Saving the specific lf values for them can lead to instabilities.
       !Therefore if the concentrations are low compared to the deltas, we fix them to same values as base case.
       do n = 1, NSPEC_SHL
          do i_lf = 1, NSPEC_deriv_lf
             if(abs(xn_shl_lf(i_lf,n,k,i,j)- xnew(n))*100<xnew(n))then
                xnew_lf(i_lf,n) = xn_shl_lf(i_lf,n,k,i,j)
                x_lf(i_lf,n)    = xn_shl_lf(i_lf,n,k,i,j) - Dchem_lf(i_lf,n,k,i,j)*dt*1.5
                x_lf(i_lf,n)    = max (x_lf(i_lf,n), 0.0)
             else
                xnew_lf(i_lf,n) = xnew(n)
                x_lf(i_lf,n)    = x(n)
             end if
          end do
          do i_lf = NSPEC_deriv_lf + 1, NSPEC_deriv_lf + N_lf_derivemis + NSOA
             xnew_lf(i_lf,n) = xnew(n)
             x_lf(i_lf,n)    = x(n)
          end do
       end do

       do n = 1, NSPEC_chem_lf
          n_sp = lfspec2spec(n)
          do i_lf = 1, NSPEC_deriv_lf
             xnew_lf(i_lf,n_sp) = xn_2d(n_sp,k)
             x_lf(i_lf,n_sp)    = xn_2d(n_sp,k) - Dchem_lf(i_lf,n_sp,k,i,j)*dt*1.5
             x_lf(i_lf,n_sp)    = max (x_lf(i_lf,n_sp), 0.0)
          end do
       end do
       do n = 1, NSPEC_chem_lf
          n_sp = lfspec2spec(n)
          do i_lf = NSPEC_deriv_lf + N_lf_derivemis + 1, NSPEC_deriv_lf + N_lf_derivemis + NSOA
             xnew_lf(i_lf,n_sp) = xn_2d(n_sp,k)
             x_lf(i_lf,n_sp)    = xn_2d(n_sp,k) - Dchem_lf(i_lf - N_lf_derivemis,n_sp,k,i,j)*dt*1.5 !NB: Dchem_lf has not emisderiv indices
             x_lf(i_lf,n_sp)    = max (x_lf(i_lf,n_sp), 0.0)
          end do
       end do

       !for emissions we use Dchem and not Dchem_lf. Because in some situations
       !the xnew_lf does not change because of emissions, but because of
       !differences in Dchem_lf only -> creates large derivatives.
       do n = 1, NSPEC_SHL
          do i_lf = NSPEC_deriv_lf  + 1, NSPEC_deriv_lf  + N_lf_derivemis
             xnew_lf(i_lf,n) = xnew(n)
             x_lf(i_lf,n)    = x(n)
          end do
       end do
       do n = 1, NSPEC_chem_lf
          n_sp = lfspec2spec(n)
          do i_lf = NSPEC_deriv_lf + 1, NSPEC_deriv_lf  + N_lf_derivemis
             xnew_lf(i_lf,n_sp) = xnew(n_sp)
             x_lf(i_lf,n_sp)    = x(n_sp)
          end do
       end do

       !change slightly one of the concentrations for chemical derivatives
       do i_lf = 1, NSPEC_deriv_lf
          !for n<=NSPEC_deriv_lf we have  lfspec2spec(n) = n+NSPEC_SHL
          n_sp = lfspec2spec(i_lf)
          xnew_lf(i_lf,n_sp) = xn_2d(n_sp,k) * eps1
          x_lf(i_lf,n_sp)    = xn_2d(n_sp,k) * eps1 - Dchem_lf(i_lf,n_sp,k,i,j)*dt*1.5
          x_lf(i_lf,n_sp)    = max (x_lf(i_lf,n_sp), 0.0)
       end do


       !change slightly one of the concentrations for chemical derivatives for SOA
       !xnew_lf and x_lf have derivemis indices included,  lfspec2spec and Dchem_lf have not.
       do i_lf = NSPEC_deriv_lf  + N_lf_derivemis + 1, NSPEC_deriv_lf  + N_lf_derivemis + NSOA
          !for n<=NSPEC_deriv_lf we have  lfspec2spec(n) = n+NSPEC_SHL
          n_sp = lfspec2spec(i_lf - N_lf_derivemis)
          xnew_lf(i_lf,n_sp) = xn_2d(n_sp,k) * eps1
          x_lf(i_lf,n_sp)    = xn_2d(n_sp,k) * eps1 - Dchem_lf(i_lf-N_lf_derivemis,n_sp,k,i,j)*dt*1.5
          x_lf(i_lf,n_sp)    = max (x_lf(i_lf,n_sp), 0.0)
       end do

       !some reaction rates depend on start concentrations. (NB: those must not be updated within chemistry!)
!       rctA_lf(:) = rctAk_lf(:,k)
!       rctB_lf(:) = rctBk_lf(:,k)

    end if
  end subroutine lf_chem_pre

  subroutine lf_chem_mid(k,cc,coeff1,coeff2,CPINIT)
    integer,intent(in) :: k
    real, intent(in) :: cc,coeff1,coeff2,CPINIT
    integer :: n,n_sp,i_lf, Nd
    real :: xextrapol

    Nd = NSPEC_deriv_lf + N_lf_derivemis !shorter

    if (lf_fullchem .and. k > KMAX_MID-lf_Nvert) then
       do n=1,NSPEC_TOT
          n_sp = spec2lfspec(n)
          if (n_sp>=0) then
             do i_lf = 1, Nd + NSOA
                xextrapol = xnew_lf(i_lf,n) + (xnew_lf(i_lf,n)-x_lf(i_lf,n)) *cc
                xold_lf(i_lf,n) = coeff1*xnew_lf(i_lf,n) - coeff2*x_lf(i_lf,n)
                xold_lf(i_lf,n) = max( xold_lf(i_lf,n), 0.0 )
                x_lf(i_lf,n) = xnew_lf(i_lf,n)
                xnew_lf(i_lf,n) = xextrapol
                if(xnew_lf(i_lf,n) < CPINIT )then
                   xnew_lf(i_lf,n) = CPINIT
                end if
             end do
          end if
       end do
    end if
  end subroutine lf_chem_mid

  subroutine lf_chem_pos(i,j,k)
    integer,intent(in) :: i,j,k
    integer :: n,i_lf, n_sp

    if(lf_fullchem .and. k > KMAX_MID-lf_Nvert) then
       !save tendencies for each derivative
       do n = 1, NSPEC_SHL
          do i_lf = 1, NSPEC_deriv_lf
             Dchem_lf(i_lf,n,k,i,j) = (xnew_lf(i_lf,n) - xn_shl_lf(i_lf,n,k,i,j) )*dt_advec_inv
          end do
       end do
       do n = 1, NSPEC_chem_lf
          n_sp = lfspec2spec(n)
          do i_lf = 1, NSPEC_deriv_lf
             Dchem_lf(i_lf,n_sp,k,i,j) = (xnew_lf(i_lf,n_sp) - xn_2d(n_sp,k) )*dt_advec_inv
          end do
       end do
       !SOA . NB: xnew_lf scenarios have also derivemis included in indices
       do n = 1, NSPEC_chem_lf
          n_sp = lfspec2spec(n)
          do i_lf = NSPEC_deriv_lf + 1, NSPEC_deriv_lf + NSOA
             Dchem_lf(i_lf,n_sp,k,i,j) = (xnew_lf(i_lf + N_lf_derivemis ,n_sp) - xn_2d(n_sp,k) )*dt_advec_inv!NB: Dchem_lf has no emisderiv indices
          end do
       end do
       !save shl
       do n = 1, NSPEC_SHL
          do i_lf = 1, NSPEC_deriv_lf
             xn_shl_lf(i_lf,n,k,i,j) = xnew_lf(i_lf,n) !save short lives for each derivatives
          end do
          !only "base" short lived are used for SOA scenarios
       end do
       !"diagonal" are different, since they used xn_2d(n_sp,k) * eps1
       do i_lf = 1, NSPEC_deriv_lf
          n_sp = lfspec2spec(i_lf)
          Dchem_lf(i_lf,n_sp,k,i,j) = (xnew_lf(i_lf,n_sp) - xn_2d(n_sp,k) * eps1)*dt_advec_inv
       end do
       !SOA . NB: xnew_lf scenarios have also derivemis included in indices
       do i_lf = NSPEC_deriv_lf + 1, NSPEC_deriv_lf + NSOA
          n_sp = lfspec2spec(i_lf)
          Dchem_lf(i_lf,n_sp,k,i,j) = (xnew_lf(i_lf + N_lf_derivemis,n_sp) - xn_2d(n_sp,k) * eps1)*dt_advec_inv
       end do
     end if
     !we are done, reset emission counters
     nemis_primary = 0
     N_lf_derivemis = 0


  end subroutine lf_chem_pos

  subroutine lf_saveall(filename)
    !save all values of lf on disk, for future restart
    !NB: order of vertical levels does not match the value of "lev"!
    character(len=*), intent(in) ::filename
    character(len=200) :: varname
    real :: scale
    integer ::i,j,k,kk,n,iter,isrc,n1,jsec,isec,ic
    integer ::ndim,kmax,CDFtype,dimSizes(10),chunksizes(10),ncFileID
    integer ::ndim_tot,dimSizes_tot(10),chunksizes_tot(10)
    character (len=20) ::dimNames(10),dimNames_tot(10)
    type(Deriv) :: def1 ! definition of fields
    logical ::overwrite, create_var_only
    real,allocatable ::tmp_out_cntry(:,:,:)!allocate since it may be heavy for the stack TEMPORARY

    ncFileID=closedID
    call Code_timer(tim_before)

    select case(projection)
     case('lon lat')
       dimNames(3)='lon'
       dimNames(4)='lat'
       dimNames_tot(1)='lon'
       dimNames_tot(2)='lat'
    case default
       dimNames(3)='i'
       dimNames(4)='j'
       dimNames_tot(1)='i'
       dimNames_tot(2)='j'
    end select

    def1%class='LF' !written
    def1%avg=.false.      !not used
    def1%index=0          !not used
    def1%scale=1.0      !not used
    def1%name='notset'
    def1%unit='ug/m3'
    ndim_tot=3
    kmax=lf_Nvert
    scale=1.0
    CDFtype=Real4
    dimSizes=1
    dimSizes_tot(3)=lf_Nvert
    dimNames_tot(3)='lev'
    chunksizes_tot=1
    chunksizes_tot(1)=MAXLIMAX
    chunksizes_tot(2)=MAXLJMAX
    chunksizes_tot(3)=dimSizes_tot(3)
    if(me==0)write(*,*)'LF save name ',trim(filename)
    allocate(tmp_out_cntry(LIMAX,LJMAX,lf_Nvert))
    do iter=1,2
       overwrite=.false. !only used once per file
       if(iter==1)overwrite=.true.!overwrite file if it exists
       create_var_only=.false.
       if(iter==1)create_var_only=.true.!only create all variables before writing into them
       do isrc = 1, Nsources
          if (lf_src(isrc)%type/='country') then
             write(*,*)'lf_saveall: only country type implemented'
             stop
          end if
          n1=lf_src(isrc)%start
          do ic=1,Ncountry_lf+Ncountry_group_lf
             do jsec=1,Ncountrysectors_lf
                !single cell source
                isec=lf_country%sector_list(jsec)
                if(lf_country%sector_list(jsec)>=0)isec=lf_country%sector_list(jsec)
                if(ic<=Ncountry_mask_lf)then
                   write(def1%name,"(A,I2.2,A5,I0)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(mask2name(country_mask_val(ic)))
                   if(isec==0) write(def1%name,"(A,I0)")trim(lf_src(isrc)%species)//'_'//trim(mask2name(country_mask_val(ic)))
                else if(ic<=Ncountry_lf)then
                   write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country%list(ic-Ncountry_mask_lf))
                   if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country%list(ic-Ncountry_mask_lf))
                else
                   !country group
                   write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country%group(ic-Ncountry_lf)%name)
                   if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country%group(ic-Ncountry_lf)%name)
                end if
                if (lf_set%full_chem .and. .not.lf_set%EmisDer_all)then
                   if(trim(lf_src(isrc)%species) /= 'pm25' .and. trim(lf_src(isrc)%species) /= 'pm25_new' .and. trim(lf_src(isrc)%species) /= 'pmco')then 
                      write(def1%name,"(A)")trim(def1%name)//trim(EMIS_FILE(lf_src(isrc)%iem_deriv))
                   else
                      write(def1%name,"(A)")trim(def1%name)//"_P"
                   end if
                end if
                if (lf_set%full_chem .and. lf_set%EmisDer_all)write(def1%name,"(A)")trim(def1%name)

                if(iter==2)then
                   !must transpose array, and reverse k
                   kk=0
                   do k = KMAX_MID,KMAX_MID-lf_Nvert+1,-1
                      kk=kk+1
                      do j=1,ljmax
                         do i=1,limax
                            tmp_out_cntry(i,j,kk) = lf(n1,i,j,k)
                         end do
                      end do
                   end do
                end if
                n1 = n1 + 1
                if(me==0 .and. create_var_only.and.DEBUG)write(*,*)'creating ',trim(def1%name)
                if(me==0 .and. .not. create_var_only.and.DEBUG)write(*,*)'saving ',trim(def1%name)
                call Out_netCDF(IOU_YEAR,def1,ndim_tot,kmax,tmp_out_cntry,scale,CDFtype,dimSizes_tot,dimNames_tot,&
                     fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                overwrite=.false.
             end do
           end do
       end do
    end do
    call CloseNetCDF(ncFileID)
    deallocate(tmp_out_cntry)
    call Add_2timing(NTIMING-2,tim_after,tim_before,"lf: output")
  end subroutine lf_saveall

  subroutine lf_read(filename, needed)
    !read all values of lf from disk
    !NB: order of vertical levels does not match the value of "lev"!
    character(len=*), intent(in)  :: filename
    logical, optional, intent(in) :: needed     

    integer :: i,j,k,kk,isrc,n1,jsec,isec,ic
    integer :: ncFileID
    type(Deriv) :: def1 ! definition of fields
    real, allocatable :: tmp_out_cntry(:,:,:)
    logical :: found

    if (present(needed)) then
      inquire(file=filename, exist=found)
      if (.not. found) then
         call CheckStop(needed, "Missing LF restart file "//trim(filename))
         if(me==0) write(*,"('WARNING: missing LF restart file ',A)") trim(filename)
         return
      end if
   end if

    ncFileID=closedID
    call Code_timer(tim_before)
    if(me==0) write(*,*)'Initializing LF from '//trim(filename)

    allocate(tmp_out_cntry(LIMAX,LJMAX,lf_Nvert))
    do isrc = 1, Nsources
      call CheckStop(lf_src(isrc)%type/='country', 'lf_read: only country type implemented')
       n1=lf_src(isrc)%start
       do ic=1,Ncountry_lf+Ncountry_group_lf
          do jsec=1,Ncountrysectors_lf
             !single cell source
             isec=lf_country%sector_list(jsec)
             if(lf_country%sector_list(jsec)>=0)isec=lf_country%sector_list(jsec)
             if(ic<=Ncountry_mask_lf)then
                write(def1%name,"(A,I2.2,A5,I0)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(mask2name(country_mask_val(ic)))
                if(isec==0) write(def1%name,"(A,I0)")trim(lf_src(isrc)%species)//'_'//trim(mask2name(country_mask_val(ic)))
             else if(ic<=Ncountry_lf)then
                write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country%list(ic-Ncountry_mask_lf))
                if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country%list(ic-Ncountry_mask_lf))
             else
                !country group
                write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country%group(ic-Ncountry_lf)%name)
                if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country%group(ic-Ncountry_lf)%name)
             end if
             if (lf_set%full_chem .and. .not.lf_set%EmisDer_all)then
                   if(trim(lf_src(isrc)%species) /= 'pm25' .and. trim(lf_src(isrc)%species) /= 'pm25_new' .and. trim(lf_src(isrc)%species) /= 'pmco')then 
                   write(def1%name,"(A)")trim(def1%name)//trim(EMIS_FILE(lf_src(isrc)%iem_deriv))
                else
                   write(def1%name,"(A)")trim(def1%name)//"_P"
                end if
             end if
             if (lf_set%full_chem .and. lf_set%EmisDer_all)write(def1%name,"(A)")trim(def1%name)

             if(me==0 .and. DEBUG)write(*,*)'Reading ',trim(def1%name),ncFileID
             call GetCDF_modelgrid(def1%name,fileName,tmp_out_cntry,1,lf_Nvert,1,1,i_start=2-RUNDOMAIN(1), j_start=2-RUNDOMAIN(3), &
                needed=needed, found=found, ncFileID_in=ncFileID)
             if (.not. found) then
                if(me==0) write(*,"('WARNING: missing ',A,' in LF restart file')") trim(def1%name)
                cycle
             end if
             !must transpose array and reverse k
             kk=0
             do k = KMAX_MID,KMAX_MID-lf_Nvert+1,-1
                kk=kk+1
                do j=1,ljmax
                   do i=1,limax
                      lf(n1,i,j,k) = tmp_out_cntry(i,j,kk)
                   end do
                end do
             end do
             n1 = n1 + 1
          end do
       end do
    end do

    call CloseNetCDF(ncFileID)
    deallocate(tmp_out_cntry)
    call Add_2timing(NTIMING-2,tim_after,tim_before,"lf: output")
  end subroutine lf_read


subroutine MeanMasksOut(iotyp, var, varname, runname, Nrun, Nruntot, Runstart, close_now)
  real, intent(in) :: var(LIMAX,LJMAX,Nrun)
  character(len=*) , intent(in) :: varname
  character(len=TXTLEN_NAME), intent(in):: runname(Nruntot)
  integer, intent(in) :: Nrun, Nruntot, Runstart, iotyp
  logical, intent(in), optional :: close_now

  character(len=TXTLEN_FILE) :: filename
  character(len=TXTLEN_FILE), save :: old_file='NOTSET'
  logical, save :: first_call=.true., overwrite(IOU_MAX_MAX)=.true.
  integer :: irun,imask
  real(kind=4), dimension(Nrun,NEmisMask) :: values
  real(kind=4), allocatable, save :: masksum(:)
 
  if(first_call) then
    allocate(masksum(NEmisMask))
    first_call = .false.

    call masked_total(masksum, 1)
    if(me/=0) deallocate(masksum)
  end if

  call masked_total(values, Nrun, field=var)
  if(me/=0) return

  forall(irun=1:Nrun,imask=1:NEmisMask)
    values(irun,imask) = values(irun,imask)/masksum(imask)
  end forall
 
  select case (iotyp)
  case (IOU_HOUR,IOU_HOUR_INST)
    filename=trim(runlabel1)//'_LF_mask'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
    overwrite(iotyp) = (filename /= old_file)
    if(overwrite(iotyp)) old_file = filename
  case (IOU_DAY)
    filename=trim(runlabel1)//'_LF_mask_day.nc'
  case (IOU_MON)
    filename=trim(runlabel1)//'_LF_mask_month.nc'
  case (IOU_YEAR)
    filename=trim(runlabel1)//'_LF_mask_full.nc'
  case default
    call CheckStop("MeanMasksOut: unknown iotyp")
  end select
  call masked_output(iotyp, filename, values, varname, runname, Nrun, Nruntot, Runstart, &
    overwrite=overwrite(iotyp))
  overwrite(iotyp) = .false.

contains
subroutine masked_total(total, Nrun, field)
  integer, intent(in)        :: Nrun
  real(kind=4),  intent(out) :: total(Nrun,NEmisMask)
  real, optional, intent(in) :: field(LIMAX,LJMAX,Nrun)

  if(present(field))then
    forall(irun=1:Nrun,imask=1:NEmisMask) ! EmisMaskValues = 1 outside city, 0 inside
      total(irun,imask) = sum(field(:,:,irun)*(1.0 - EmisMaskValues(:,:,imask)))
    end forall
  else
    forall(imask=1:NEmisMask) ! EmisMaskValues = 1 outside city, 0 inside
      total(1,imask) = sum(1.0 - EmisMaskValues(:,:,imask))
    end forall
  end if

  ! add together totals from each processor (only me=0 get results)
  if(me==0)then
    call MPI_REDUCE(MPI_IN_PLACE,total,NEmisMask*Nrun,MPI_REAL,MPI_SUM,0,MPI_COMM_CALC,IERROR)
  else
    call MPI_REDUCE(total,total,NEmisMask*Nrun,MPI_REAL,MPI_SUM,0,MPI_COMM_CALC,IERROR)
  end if
end subroutine masked_total
end subroutine MeanMasksOut

end module LocalFractions_mod
