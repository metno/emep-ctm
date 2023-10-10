! <LocalFractions_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.0>
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
module LocalFractions_mod
!
! all subroutines for Local Fractions
!
use CheckStop_mod,     only: CheckStop,StopAll
use Chemfields_mod,    only: xn_adv, cfac
use ChemDims_mod,      only: NSPEC_ADV, NSPEC_SHL,NSPEC_TOT,NEMIS_File
use ChemFunctions_mod, only: EC_AGEING_RATE
use ChemSpecs_mod,     only: species_adv,species
use Config_module,     only: NPROC,KMAX_MID, KMAX_BND,KCHEMTOP,USES, lf_src, IOU_HOUR&
                             , IOU_HOUR_INST,IOU_INST,IOU_YEAR,IOU_MON,IOU_DAY&
                             ,IOU_HOUR,IOU_HOUR_INST, IOU_MAX_MAX &
                             ,MasterProc,dt_advec, RUNDOMAIN, runlabel1 &
                             ,HOURLYFILE_ending, lf_species, lf_country,&
                             SO2_ix, O3_ix, NO2_ix, SO4_ix, NH4_f_ix, NO3_ix, NO3_f_ix, &
                             NO3_c_ix, NH3_ix, HNO3_ix, C5H8_ix, NO_ix, HO2_ix, OH_ix,&
                             HONO_ix,OP_ix,CH3O2_ix,C2H5O2_ix,CH3CO3_ix,C4H9O2_ix,MEKO2_ix,ETRO2_ix,&
                             PRRO2_ix,OXYO2_ix,C5DICARBO2_ix,ISRO2_ix,MACRO2_ix,TERPO2_ix,H2O2_ix,N2O5_ix
use Convection_mod,    only: convection_1d
use Country_mod,       only: MAXNLAND,NLAND,Country&
                             ,IC_TMT,IC_TM,IC_TME,IC_ASM,IC_ASE,IC_ARE,IC_ARL,IC_CAS,IC_UZT,IC_UZ&
                             ,IC_UZE,IC_KZT,IC_KZ,IC_KZE,IC_RU,IC_RFE,IC_RUX,IC_RUE,IC_AST

use DefPhotolysis_mod, only: IDHONO,IDNO3,IDNO2
use EmisDef_mod,       only: NSECTORS,SECTORS,EMIS_FILE, &
                             nlandcode,landcode,NCMAX,&
                             TFAC_IDX_DOM,secemis, roaddust_emis_pot,KEMISTOP,&
                             EmisMaskIntVal,gridrcemis,&
                             mask2name
use EmisGet_mod,       only: nrcemis, iqrc2itot, emis_nsplit,nemis_kprofile, emis_kprofile,&
                             emis_masscorr,  &  ! 1/molwt for most species
                             make_iland_for_time,itot2iqrc,iqrc2iem, emisfrac
use GridValues_mod,    only: dA,dB,xm2, dhs1i, glat, glon, projection, extendarea_N,i_fdom,j_fdom,&
                             RestrictDomain
use MetFields_mod,     only: ps,roa,EtaKz
use MPI_Groups_mod
use NetCDF_mod,        only: Real4,Real8,Out_netCDF,LF_ncFileID_iou
use OwnDataTypes_mod,  only: Deriv, Max_lf_sources, Max_lf_sectors, MAX_lf_country_group_size, &
                             Max_lf_spec, TXTLEN_NAME, TXTLEN_FILE, &
                             Max_lf_res, Max_lf_Country_list, Max_lf_sectors, Max_lf_Country_groups
use Par_mod,           only: me,LIMAX,LJMAX,MAXLIMAX,MAXLJMAX,gi0,gj0,li0,li1,lj0,lj1,GIMAX,GJMAX
use PhysicalConstants_mod, only : GRAV, AVOG, ATWAIR
use SmallUtils_mod,    only: find_index
use TimeDate_mod,      only: date, current_date,day_of_week
use TimeDate_ExtraUtil_mod,only: date2string
use My_Timing_mod,     only: Add_2timing, Code_timer, NTIMING
use VerticalDiffusion_mod, only: vertdiffn
use ZchemData_mod,only: rct, rcphot, xn_2d, rcemis, M

!(dx,dy,i,j) shows contribution of pollutants from (i+dx,j+dy) to (i,j)

implicit none
!external advection_mod_mp_vertdiffn_k

private

integer ::IC_AST_EXTRA = 324567,IC_RUT_EXTRA = 324568 !sum of countries which are not defined as countries
integer ::IC_BIC_EXTRA = 324569
integer ::IC_STRATOS = 324566 !contribution from stratosphere
integer ::IC_INIT = 324565 !contribution from pollutants present at the start of the run
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
public  :: lf_aero_pre, lf_aero_pos
public  :: lf_drydep, lf_wetdep
public  :: lf_rcemis
public  :: save_lf_emis

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

real, public, allocatable, dimension(:,:,:,:,:), save :: emis_lf_cntry
real, public, allocatable, dimension(:,:,:,:), save :: &
   loc_frac_src &   ! Fraction of pollutants that are produced locally, list of defined sources
  ,lf &   ! Fraction of pollutants that are produced locally, for all defined sources
  ,loc_frac_src_full &   ! Fraction of pollutants that are produced locally, list of defined sources
  ,lf_src_full &   ! Fraction of pollutants that are produced locally, list of defined sources
  ,loc_tot_hour_inst&   !all contributions
  ,loc_tot_hour&   !Hourly average of all contributions
  ,loc_tot_day&   !Daily average of all contributions
  ,loc_tot_month&  !Monthly average of all contributions
  ,loc_tot_full  !Fullrun average of all contributions
real, public, allocatable, dimension(:,:,:), save :: &
  loc_frac_drydep  ! ddepositions per source (not fractions!)
real, public, allocatable, dimension(:,:,:), save :: &
  loc_frac_wetdep  ! wdepositions per source (not fractions!)
real, public, allocatable, dimension(:,:,:,:), save :: D8M !
real, public, allocatable, dimension(:,:,:), save :: D8Max !
real, public, allocatable, dimension(:,:,:,:), save :: D8Max_av !
real, public, allocatable, dimension(:,:,:), save :: hourM !
real, public, allocatable, dimension(:,:,:,:), save :: &
  loc_frac_1d  ! Fraction of pollutants without i or j and extended (0:limax+1 or 0:ljmax+1)
real, public, allocatable, dimension(:,:), save :: &
  loc_frac_src_1d  ! Fraction of pollutants without i or j and extended (0:limax+1 or 0:ljmax+1)
real, allocatable, save ::loc_poll_to(:,:,:,:,:)
real, allocatable, public, dimension(:,:), save ::xderiv !dX_ispec/dX_n
real, allocatable, public, dimension(:,:), save ::ederiv !dX_ispec/demis_n
real, allocatable, public, dimension(:,:), save :: lf0_loc, lf_loc
real, allocatable, public, dimension(:,:), save ::x_lf, xold_lf ,xnew_lf
real, allocatable, public, dimension(:,:,:,:,:), save ::Dchem_lf !may not be worth the cost?
real, allocatable, public, dimension(:,:,:,:,:), save ::xn_shl_lf!may not be worth the cost?
real, allocatable, public, dimension(:,:), save ::rcemis_lf
integer, allocatable, public, dimension(:,:), save ::nic
integer, allocatable, public, dimension(:,:,:), save ::ic2iland
real, allocatable, public, dimension(:), save ::L_lf, P_lf

logical, public, save :: COMPUTE_LOCAL_TRANSPORT=.false.
integer , public, save :: lf_Nvertout = 1!number of vertical levels to save in output
integer, public, save :: NTIMING_lf=9
real, private :: tim_after,tim_before
integer, public, save :: Ndiv_coarse=1, Ndiv_rel=1, Ndiv2_coarse=1
integer, public, save :: Nsources=0, Nsources_nonew=0
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
integer, private, save :: isrc_O3=-1, isrc_O3_voc=-1, isrc_NO=-1, isrc_NO2=-1, isrc_VOC=-1
integer, private, save :: isrc_SO2=-1, isrc_SO4=-1, isrc_NH4=-1, isrc_NH3=-1
integer, private, save :: isrc_NO3=-1, isrc_HNO3=-1
integer, private, save :: isrc_pm25_new=-1, isrc_pm25=-1
integer, private, save :: isrc_strato=-1, isrc_ini=-1
integer, private, save :: RO2POOL_ix=-1
real, allocatable, private, save :: lf_NH4(:), lf_NH3(:)
real, allocatable, private, save :: lf_NO3(:), lf_HNO3(:)
integer, private, save :: country_ix_list(Max_lf_Country_list)
integer, private, save :: Ncountry_lf=0
integer, private, save :: Ncountry_group_lf=0
integer, private, save :: Ncountrysectors_lf=0
integer, private, save :: Ncountry_mask_lf=0 !total number of masks defined
integer, private, save :: Ncountry_mask_lf_val=0 !number of masks defined using lf_country%mask_val
integer, private, save :: country_mask_val(Max_lf_Country_list) = -999999 ! values of all defined masks
character(len=TXTLEN_NAME), private, save :: iem2names(NEMIS_File,Max_lf_spec) !name of that pollutant
integer, private, save :: isrc_new(Max_lf_sources)
integer, private, save :: Stratos_ix(1000) !1000 must be larger than Nsources
integer, private, save :: nstratos !number of sources to track for Stratos
real   , private, save :: P_NO(100),P_NO2(100)
integer, private, save :: nfullchem=0 ! >0 if the full O3 chemistry is needed
logical, public, save :: lf_fullchem=.false. ! if the full O3 chemistry is needed
integer, public, save :: NSPEC_fullchem_lf=0 ! number of species to include in the "fullchem" derivatives
integer, public, save :: NSPEC_fullchem_inc_lf=0 ! number of species to included in CM_Reactions1
integer, public, parameter :: N_lf_derivemisMAX = 100 ! max number of emission source to include in CM_Reactions1 derivatives
integer, public, save :: N_lf_derivemis ! actual number of emissions to include in CM_Reactions1 derivatives
integer, public, parameter :: iem_lf_nox = 1, iem_lf_voc = 2
integer, public, save :: emis2icis(N_lf_derivemisMAX),emis2isrc(N_lf_derivemisMAX)

contains

  subroutine lf_init
    integer :: n, n0, is, i, j, ii, iii, ic, ix, iix, isrc, isec, n_mask, mask_val_min, mask_val_max
    integer :: found, itot, iqrc, iem, iemis, ipoll, ixnh3, ixnh4, size, IOU_ix, iem_deriv
    integer, allocatable :: MaskVal(:)
! pm25_new and pm25 are considered as two different emitted pollutants
    if(DEBUGall .and. me==0)write(*,*)'start init'

  call Code_timer(tim_before)
  ix=0
  if(USES%uEMEP)then
     call StopAll("USES%uEMEP no longer in use. Use lf_ syntaks")
  else
     !separate values do not work properly yet
     lf_src(:)%dist = lf_src(1)%dist !Temporary
     do i=1,4
        lf_src(:)%DOMAIN(i) = lf_src(1)%DOMAIN(i) !Temporary
     enddo
     lf_src(:)%YEAR=lf_src(1)%YEAR
     lf_src(:)%MONTH=lf_src(1)%MONTH
     lf_src(:)%MONTH_ENDING=lf_src(1)%MONTH_ENDING
     lf_src(:)%DAY=lf_src(1)%DAY
     lf_src(:)%HOUR=lf_src(1)%HOUR
     lf_src(:)%HOUR_INST=lf_src(1)%HOUR_INST
  endif

  lf_Nvert = lf_src(1)%Nvert !Temporary
  if (lf_Nvert>KMAX_MID-1) then
     lf_Nvert = KMAX_MID-1
     if(me==0)then
        write(*,*)'WARNING: cannot track through level 1 (top). Reducing Nvert to ',lf_Nvert
     end if
  end if
  Nsources = 0
  nfullchem = 0
  do i = 1, Max_lf_sources
     if (lf_src(i)%species == 'NOTSET') exit
     if (lf_src(i)%species == 'FULLCHEM') then
        nfullchem = nfullchem + 1
        cycle
     end if
     Nsources = Nsources + 1
  enddo
  if (nfullchem>0) then
     lf_fullchem = .true.

     NSPEC_fullchem_lf = NSPEC_ADV !default include all
     !We ASSUME that shipNOx is the last species (highest index) involved in O3 chemistry
     ix=find_index("shipNOx" ,species_adv(:)%name)
     if(ix>0)NSPEC_fullchem_lf = ix

     !TEMPORARY, should be set to NSPEC_fullchem_lf+NSPEC_SHL
     !We ASSUME that SQT_SOA_NV is the last species (highest index) included in CM_Reactions1.inc
     !ix=find_index("SQT_SOA_NV" ,species(:)%name) !NB: index among all species, also SHL
     NSPEC_fullchem_inc_lf = ix
     NSPEC_fullchem_inc_lf = max(ix, NSPEC_fullchem_lf+NSPEC_SHL)

     if(me==0)write(*,*)'LF chemistry, number of chem derivatives calculated: ',NSPEC_fullchem_lf
     if(me==0)write(*,*)'LF chemistry, max index of species included: ',NSPEC_fullchem_inc_lf

     !RO2POOL is not a proper chemical species and we must take care to avoid it when using derivatives!
     RO2POOL_ix = find_index('RO2POOL' ,species(:)%name)

  end if

  do i = Nsources + nfullchem + 1, Max_lf_sources
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
  do isrc = 1, Max_lf_sources
     if (lf_src(isrc)%species == 'NOTSET') exit
     if (lf_src(isrc)%species == 'FULLCHEM') then
        if(me==0)write(*,*)isrc,' FULLCHEM '
        call CheckStop(lf_src(isrc)%type /= 'country',"LocalFractions: only country type for FULLCHEM. Not implemented "//trim(lf_src(isrc)%type)) ! Only country sources for now
        lf_src(isrc)%Npos = 0
        lf_src(isrc)%Nsplit = 0
        iem = find_index('nox' ,EMIS_FILE(1:NEMIS_FILE))
        call CheckStop(iem<=0, "LF: did not find nox emissions")
        iem = find_index('voc' ,EMIS_FILE(1:NEMIS_FILE))
        call CheckStop(iem<=0, "LF: did not find voc emissions")
        do i = 1, NSPEC_fullchem_lf !first source per advected species, tracking NOx
           Nsources = Nsources + 1
           lf_src(Nsources)%type = lf_src(isrc)%type
           lf_src(Nsources)%species = species_adv(i)%name
           lf_src(Nsources)%full_chem = .true.
           !iem_deriv differs from iem, because it relates not to the source species, but the derivative.
           !iem_deriv can be for voc, while the species is no2
           iem_deriv = find_index('nox' ,EMIS_FILE(1:NEMIS_FILE))
           lf_src(Nsources)%iem_deriv = iem_deriv
           lf_src(Nsources)%iem_lf = iem_lf_nox
           lf_src(Nsources)%make_fracsum = lf_src(isrc)%make_fracsum
           !lf_src(Nsources)%ix(1) = i set later also
           if (lf_src(Nsources)%species =='O3') isrc_O3 = Nsources
           if (lf_src(Nsources)%species =='NO2') isrc_NO2 = Nsources
           if (lf_src(Nsources)%species =='NH4_f') isrc_NH4 = Nsources
           if (lf_src(Nsources)%species =='NH3') isrc_NH3 = Nsources
           if (lf_src(Nsources)%species =='NO3') isrc_NO3 = Nsources
           if (lf_src(Nsources)%species =='HNO3') isrc_HNO3 = Nsources
        end do

        do i = 1, NSPEC_fullchem_lf !second source per advected species, tracking NMVOC
           Nsources = Nsources + 1
           lf_src(Nsources)%type = lf_src(isrc)%type
           lf_src(Nsources)%species = species_adv(i)%name
           lf_src(Nsources)%full_chem = .true.
           iem_deriv = find_index('voc' ,EMIS_FILE(1:NEMIS_FILE))
           lf_src(Nsources)%iem_deriv = iem_deriv
           lf_src(Nsources)%iem_lf = iem_lf_voc
           lf_src(Nsources)%make_fracsum = lf_src(isrc)%make_fracsum
           if (lf_src(Nsources)%species =='O3') isrc_O3_voc = Nsources
           !lf_src(Nsources)%ix(1) = i set later also
        end do

     end if
  end do

  !for each pm25 we should separate into new and age parts
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
     end if
  enddo

  !countries
  ! The lowest indices are for mask indices, the highest indices are for groups
  if(lf_country%mask_val_min <= lf_country%mask_val_max .or. &
       lf_country%mask_val(1) > -999999 .or. &
       lf_country%list(1)/= 'NOTSET' .or. &
       lf_country%group(1)%name/= 'NOTSET')then

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

     if (mask_val_max >= mask_val_min) deallocate(MaskVal)

     if (Ncountry_mask_lf>0 .and. MasterProc) then
        write(*,*)'including in total',Ncountry_mask_lf,'mask sources:'
        do n = 1, (Ncountry_mask_lf+29)/30
           write(*,fmt="(30(I0,1x))")(country_mask_val(ii),ii=(n-1)*30+1, min(Ncountry_mask_lf,n*30))
        end do
     end if

     if(lf_country%list(1)/= 'NOTSET' .or. lf_country%group(1)%name/= 'NOTSET')then
        !list of countries/sectors instead of single country
        do i = 1, Max_lf_Country_list
           if(lf_country%list(i) == 'NOTSET') exit
           Ncountry_lf=Ncountry_lf+1
           ix = find_index(trim(lf_country%list(i)) ,Country(:)%code, first_only=.true.)
           if (ix<0 .and. lf_country%list(i) =='AST') then
              ix = IC_AST_EXTRA
           else if (ix<0 .and. lf_country%list(i) =='RUT') then
              ix = IC_RUT_EXTRA
           else if (ix<0 .and. lf_country%list(i) =='BIC') then
              ix = IC_BIC_EXTRA
           else if(ix<0 .and. lf_country%list(i) =='STRATOS')then
              ix = IC_STRATOS
           else if(ix<0 .and. lf_country%list(i) =='INIT')then
              ix = IC_INIT
           endif
           call CheckStop(ix<0,'country '//trim(lf_country%list(i))//' not defined. ')
           country_ix_list(i + Ncountry_mask_lf) = ix
           if(MasterProc)write(*,*)'include sources from ',trim(lf_country%list(i))
        enddo
     endif

     Ncountry_group_lf=0
     do i = 1, Max_lf_Country_groups
        if(lf_country%group(i)%name == 'NOTSET') exit
        Ncountry_group_lf = Ncountry_group_lf+1
        do ic = 1, MAX_lf_country_group_size
           if(lf_country%group(i)%list(ic) == 'NOTSET') exit
           ix = find_index(trim(lf_country%group(i)%list(ic)) ,Country(:)%code, first_only=.true.)
           if(ix<0 .and. lf_country%list(i) =='AST')then
              ix = IC_AST_EXTRA
           else if(ix<0 .and. lf_country%list(i) =='RUT')then
              ix = IC_RUT_EXTRA
           else if(ix<0 .and. lf_country%list(i) =='BIC')then
              ix = IC_BIC_EXTRA
           endif
           call CheckStop(ix<0,'country '//trim(lf_country%group(i)%list(ic))//' not defined. ')
           lf_country%group(i)%ix(ic) = ix
           if(MasterProc)write(*,*)'include sources from '//&
                trim(lf_country%group(i)%list(ic))//' as '//trim(lf_country%group(i)%name)
        enddo
     end do
     Ncountrysectors_lf=0
     do i = 1, Max_lf_sectors
        if(lf_country%sector_list(i) < 0) exit
        Ncountrysectors_lf=Ncountrysectors_lf+1
        if(MasterProc)write(*,*)'country sector ',lf_country%sector_list(i)
     end do
     do isrc = 1, Nsources
        lf_src(isrc)%Npos = (Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf
     end do
     if(MasterProc)write(*,*)(Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf,' countries x sectors for ',Nsources,' sources'
  end if

  ipoll=0
  iem2ipoll = -1
  iem2Nipoll = 0
  do isrc = 1, Nsources
     !for now only one Ndiv possible for all sources
     if(lf_src(isrc)%type == 'relative')then
        lf_src(isrc)%Npos =  (2*lf_src(isrc)%dist+1)*(2*lf_src(isrc)%dist+1)
        Ndiv_rel = max(Ndiv_rel,2*lf_src(isrc)%dist+1)
     endif

     if(lf_src(isrc)%country_ISO /= 'NOTSET')then
        lf_src(isrc)%type = 'country'
        ix = find_index(trim(lf_src(isrc)%country_ISO) ,Country(:)%code, first_only=.true.)
        if(ix<0)then
           if(me==0)write(*,*)'LF: WARNING: country '//trim(lf_src(isrc)%country_ISO)//' not defined. '
        endif
        lf_src(isrc)%country_ix = ix
        if(MasterProc)write(*,*)isrc,' country '//trim(lf_src(isrc)%country_ISO)//' '//trim(lf_src(isrc)%species)
     endif

     call RestrictDomain(lf_src(isrc)%DOMAIN)

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
              isrc_pm25_new = isrc
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
              endif
           endif
           call CheckStop( ix<1, "Local Fractions did not find corresponding pollutant: "//trim(lf_src(isrc)%species) )
           iem=-1
           lf_src(isrc)%species_ix = ix !NB: index among all species
           lf_src(isrc)%ix(1) = ix - NSPEC_SHL !NB: index among advected species
           lf_src(isrc)%mw(1) = species_adv(lf_src(isrc)%ix(1))%molwt
           lf_src(isrc)%iqrc = itot2iqrc(ix) !negative if not among emitted species
           if(lf_src(isrc)%iqrc>0) iem = iqrc2iem(lf_src(isrc)%iqrc)
           lf_src(isrc)%iem = iem
           if (lf_src(Nsources)%full_chem .and. iem>0) then
              !Do not include emissions which are not NOx or VOC emissions
              if(isrc>Nsources-2*NSPEC_fullchem_lf) then
                 ! sources that track NOx
                 if (EMIS_FILE(iem)/='nox' .and. EMIS_FILE(iem)/='voc') lf_src(isrc)%iem = -1
              end if
           end if

           if (.not. lf_src(Nsources)%full_chem) then
              !use simplified methods
              if(trim(species(ix)%name)=='O3')isrc_O3=isrc
              if(trim(species(ix)%name)=='NO')isrc_NO=isrc
              if(trim(species(ix)%name)=='NO2')isrc_NO2=isrc
              if(trim(species(ix)%name)=='SO4')isrc_SO4=isrc
              if(trim(species(ix)%name)=='SO2')isrc_SO2=isrc
              if(trim(species(ix)%name)=='NH4_f')isrc_NH4=isrc
              if(trim(species(ix)%name)=='NH3')isrc_NH3=isrc
           end if
        end if
     else
        !species defines as primary emitted
        lf_src(isrc)%iem = iem
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
           endif
           if(lf_src(isrc)%species=="nox")then
              ix=find_index("NO2",species_adv(:)%name)
              call CheckStop(ix<0,'Index for NO2 not found')
              lf_src(isrc)%mw(ii)=species_adv(ix)%molwt !"as NO2"
           endif
           if(lf_src(isrc)%species=="sox")then
              ix=find_index("SO2",species_adv(:)%name)
              call CheckStop(ix<0,'Index for SO2 not found')
              lf_src(isrc)%mw(ii)=species_adv(ix)%molwt !"as SO2"
           endif
           if(lf_src(isrc)%species=="nox" .and. (lf_src(isrc)%DryDep .or. lf_src(isrc)%WetDep))then
              ix=find_index("NO3",species_adv(:)%name)
              call CheckStop(ix<0,'Index for NO3 not found')
              ix=find_index("HNO3",species_adv(:)%name)
              call CheckStop(ix<0,'Index for HNO3 not found')
           endif

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
                 endif
              enddo
           end if
        end do

     endif

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
        if (lf_src(isrc)%iem>0) then
           write(*,*)'lf pollutant : ',lf_src(isrc)%species,' ref index ',lf_src(isrc)%poll,' emitted as ',EMIS_FILE(lf_src(isrc)%iem)
        else
           write(*,*)'lf pollutant : ',lf_src(isrc)%species,' ref index ',lf_src(isrc)%poll,' not treated as emitted species'
        end if
        write(*,*)'lf number of species in '//trim(lf_src(isrc)%species)//' group: ',lf_src(isrc)%Nsplit
        write(*,"(A,30(A,F6.2))")'including:',('; '//trim(species_adv(lf_src(isrc)%ix(i))%name)//', mw ',lf_src(isrc)%mw(i),i=1,lf_src(isrc)%Nsplit)
        if (lf_src(isrc)%type/='country') write(*,"(A,I4,A,I4)")'sector:',lf_src(isrc)%sector,',  res:',lf_src(isrc)%res
        !write(*,"(A,30I4)")'ix:',(lf_src(isrc)%ix(i),i=1,lf_src(isrc)%Nsplit)
     end if
  end do


  if (isrc_SO2>0 .and. (isrc_SO4<0)) then
     if(me==0)write(*,*)'WARNING: SO2 tracking requires SO4'
     stop!may be relaxed in future
  end if

  av_fac=0.0

  Ndrydep_lf=0
  Nwetdep_lf=0
  LF_SRC_TOTSIZE = 0
  do isrc = 1, Nsources
     if(lf_src(isrc)%drydep) Ndrydep_lf = Ndrydep_lf + lf_src(isrc)%Npos
     if(lf_src(isrc)%wetdep) Nwetdep_lf = Nwetdep_lf + lf_src(isrc)%Npos
     if (lf_src(isrc)%WetDep) then
        do iix=1,lf_src(isrc)%Nsplit
           ix=lf_src(isrc)%ix(iix)
           wetdep_lf(ix) = .true.
        end do
     end if
     lf_src(isrc)%start = LF_SRC_TOTSIZE + 1
     lf_src(isrc)%end = LF_SRC_TOTSIZE + lf_src(isrc)%Npos
     LF_SRC_TOTSIZE = LF_SRC_TOTSIZE + lf_src(isrc)%Npos
     if(me==0)then
        write(*,*)isrc,' ',trim(lf_src(isrc)%species)," start ",lf_src(isrc)%start," end ",lf_src(isrc)%end,LF_SRC_TOTSIZE
     end if

  end do

  allocate(lf(LF_SRC_TOTSIZE,LIMAX,LJMAX,KMAX_MID-lf_Nvert+1:KMAX_MID))
  lf=0.0


  nstratos = 0
  do isrc = 1, Nsources
     !IC_STRATOS is special. Make a list of corresponding sources indices
     !IC_INIT is special. init lf to 1.
     if(lf_src(isrc)%type=='country' .and. Ncountry_lf+Ncountry_group_lf>0)then
        n0=lf_src(isrc)%start
        do ic=1,Ncountry_lf+Ncountry_group_lf
           do is=1,Ncountrysectors_lf
              if (country_ix_list(ic)==IC_STRATOS) then
                 nstratos = nstratos + 1
                 if (nstratos > size(Stratos_ix) )then
                    write(*,*)'Increase size of Stratos_ix to ',Nsources, size(Stratos_ix)
                    stop
                 end if
                 Stratos_ix(nstratos)=n0
              end if
              if (country_ix_list(ic)==IC_INIT) lf(n0,:,:,:) = 1.0
              n0=n0+1
           end do
        end do
     end if
  end do

  isrc=1!for now all must be the same
  if(lf_src(isrc)%HOUR)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_HOUR)=Niou_ix
     iotyp2ix(Niou_ix)=IOU_HOUR
  endif
  if(lf_src(isrc)%HOUR_INST)then
     Niou_ix = Niou_ix + 1
     iou_ix_inst = Niou_ix !should not be accumulated
     iotyp2ix(IOU_HOUR_inst) = Niou_ix
     iotyp2ix(Niou_ix)=IOU_HOUR_inst
  endif
  if(lf_src(isrc)%DAY)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_DAY)=Niou_ix
     iotyp2ix(Niou_ix)=IOU_DAY
  endif
  if(lf_src(isrc)%MONTH)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_MON)=Niou_ix
     iotyp2ix(Niou_ix)=IOU_MON
  endif
  if(lf_src(isrc)%YEAR)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_YEAR)=Niou_ix
     iotyp2ix(Niou_ix)=IOU_YEAR
  endif

  if(isrc_NH4>0)then
     allocate(lf_NH4(KMAX_MID-lf_Nvert+1:KMAX_MID))
     allocate(lf_NH3(KMAX_MID-lf_Nvert+1:KMAX_MID))
  endif

  allocate(lf_src_acc(LF_SRC_TOTSIZE,LIMAX,LJMAX,KMAX_MID-lf_Nvertout+1:KMAX_MID,Niou_ix))
  lf_src_acc = 0.0
  allocate(lf_src_tot(LIMAX,LJMAX,KMAX_MID-lf_Nvertout+1:KMAX_MID,Npoll,Niou_ix))
  lf_src_tot = 0.0
   allocate(loc_frac_src_1d(LF_SRC_TOTSIZE,0:max(LIMAX,LJMAX)+1))
  loc_frac_src_1d=0.0

  allocate(emis_lf_cntry(LIMAX,LJMAX,NCMAX,Nsectors,NEMIS_File))
  emis_lf_cntry=0.0

  if(Ndrydep_lf>0)then
     allocate(loc_frac_drydep(LIMAX,LJMAX,Ndrydep_lf))
     loc_frac_drydep=0.0
  else
     allocate(loc_frac_drydep(1,1,1))
  endif
  if(Nwetdep_lf>0)then
     allocate(loc_frac_wetdep(LIMAX,LJMAX,Nwetdep_lf))
     loc_frac_wetdep=0.0
  else
     allocate(loc_frac_wetdep(1,1,1))
  endif
  if(lf_fullchem)then
     allocate(xderiv(NSPEC_fullchem_lf,NSPEC_fullchem_lf))
     allocate(ederiv(N_lf_derivemisMAX,NSPEC_fullchem_lf))
     allocate(lf0_loc(2*(Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf,NSPEC_fullchem_lf))
     allocate(lf_loc(2*(Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf,NSPEC_fullchem_lf))
     allocate(xnew_lf(NSPEC_fullchem_lf+N_lf_derivemisMAX,NSPEC_fullchem_inc_lf))
     allocate(x_lf(NSPEC_fullchem_lf+N_lf_derivemisMAX,NSPEC_fullchem_inc_lf))
     allocate(rcemis_lf(N_lf_derivemisMAX,NSPEC_fullchem_inc_lf))
     allocate(P_lf(NSPEC_fullchem_lf+N_lf_derivemisMAX))
     allocate(L_lf(NSPEC_fullchem_lf+N_lf_derivemisMAX))
     allocate(xold_lf(NSPEC_fullchem_lf+N_lf_derivemisMAX,NSPEC_fullchem_inc_lf))
     allocate(Dchem_lf(NSPEC_fullchem_lf,NSPEC_fullchem_inc_lf,KMAX_MID-lf_Nvert+1:KMAX_MID,LIMAX,LJMAX))
     allocate(xn_shl_lf(NSPEC_fullchem_lf,NSPEC_SHL,KMAX_MID-lf_Nvert+1:KMAX_MID,LIMAX,LJMAX))
     xnew_lf = 0.0
     x_lf = 0.0
     xold_lf = 0.0
     Dchem_lf = 0.0
     xn_shl_lf = 0.0
     allocate(lf_NO3(KMAX_MID-lf_Nvert+1:KMAX_MID))
     allocate(lf_HNO3(KMAX_MID-lf_Nvert+1:KMAX_MID))
  else
      allocate(rcemis_lf(Nsources*NCMAX,1))
  endif
  allocate(nic(LIMAX,LJMAX))
  nic = 0
  allocate(ic2iland(LIMAX,LJMAX,NCMAX))

  if (lf_src(1)%MDA8) then! maximum daily eight-hour mean concentration  
    !we need one array to save the last 8 hourly concentrations, and one to make
    !the average over the last hour
    !index 0 is for the pure O3, and the other indices are for the sensibilities
    allocate(D8M(0:lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos,LIMAX,LJMAX,8)) ! running last 8 hour values
    allocate(D8Max(0:lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos,LIMAX,LJMAX)) ! max value of the 8 hour mean since 00:00
    allocate(D8Max_av(LIMAX,LJMAX,0:lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos,Niou_ix)) ! max value of the 8 hour mean since 00:00
    allocate(hourM(0:lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos,LIMAX,LJMAX)) ! hour Mean
    
    hourM = 0.0
    D8Max = 0.0 !init with low value
    D8M = 0.0 !init with low value
    D8Max_av = 0.0 !init necessary
   end if


!  call Add_2timing(NTIMING-10,tim_after,tim_before,"lf: init") negligible
  if(DEBUGall .and. me==0)write(*,*)'end init'

end subroutine lf_init


subroutine lf_out(iotyp)
  integer, intent(in) :: iotyp
  character(len=200) ::filename, varname
  real :: xtot,scale,invtot,t1,t2
  integer ::i,j,k,n,n1,dx,dy,ix,iix,isec,iisec,isec_poll,ipoll,isec_poll1,isrc,iou_ix,iter,iddep,iwdep
  integer ::ndim,kmax,CDFtype,dimSizes(10),chunksizes(10)
  integer ::ndim_tot,dimSizes_tot(10),chunksizes_tot(10)
  character (len=20) ::dimNames(10),dimNames_tot(10)
  type(Deriv) :: def1 ! definition of fields for local fraction
  type(Deriv) :: def2 ! definition of fields for totals
  type(Deriv) :: def3 ! definition of dry and wet dep fields
  logical ::overwrite, create_var_only
  logical,save :: first_call(10)=.true.
  real,allocatable ::tmp_out(:,:,:)!allocate since it may be heavy for the stack TEMPORARY
  real,allocatable ::tmp_out_cntry(:,:,:)!allocate since it may be heavy for the stack TEMPORARY
  type(date) :: onesecond = date(0,0,0,0,1)
  character(len=TXTLEN_FILE),save :: oldhourlyname = 'NOTSET'
  character(len=TXTLEN_FILE),save :: oldhourlyInstname = 'NOTSET'
  character(len=TXTLEN_FILE),save :: oldmonthlyname
  real :: fracsum(LIMAX,LJMAX),invfac
  logical :: pollwritten(Max_lf_spec)
  integer :: ncFileID

  call Code_timer(tim_before)
  if(DEBUGall .and. me==0)write(*,*)'start out'

  if(iotyp==IOU_HOUR_INST .and. lf_src(1)%HOUR_INST)then
     fileName = trim(runlabel1)//'_LF_hourInst'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
     if(oldhourlyInstname/=fileName)then
        first_call(iotyp) = .true.
        oldhourlyInstname = fileName
     endif
  else if(iotyp==IOU_HOUR .and. lf_src(1)%HOUR)then
     fileName = trim(runlabel1)//'_LF_hour'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
     if(oldhourlyname/=fileName)then
        first_call(iotyp) = .true.
        oldhourlyname = fileName
     endif
  else if(iotyp==IOU_DAY .and. lf_src(1)%DAY)then
     fileName=trim(runlabel1)//'_LF_day.nc'
  else if(iotyp==IOU_MON .and. lf_src(1)%MONTH)then
     if(lf_src(1)%MONTH_ENDING /= "NOTSET")then
        fileName=trim(runlabel1)//'_LF_month'//date2string(trim(lf_src(1)%MONTH_ENDING),current_date,-1.0)
        if(oldmonthlyname/=fileName)then
           first_call(iotyp) = .true.
           oldmonthlyname = fileName
        endif
     else
        fileName=trim(runlabel1)//'_LF_month.nc'
     endif
  else if(iotyp==IOU_YEAR .and. lf_src(1)%YEAR)then
     fileName=trim(runlabel1)//'_LF_full.nc'
  else
     return
  endif
  ncFileID=LF_ncFileID_iou(iotyp)

  ndim=5
  ndim_tot=3
  kmax=lf_Nvertout
  scale=1.0
  CDFtype=Real4
  dimSizes=1

  dimSizes(1)=2*lf_src(1)%dist+1
  dimNames(1)='x_dist'
  dimSizes(2)=2*lf_src(1)%dist+1
  dimNames(2)='y_dist'

  isrc=1!temporary
  dimSizes(3)=min(GIMAX,lf_src(isrc)%DOMAIN(2)-lf_src(isrc)%DOMAIN(1)+1)
  dimSizes(4)=min(GJMAX,lf_src(isrc)%DOMAIN(4)-lf_src(isrc)%DOMAIN(3)+1)

  dimSizes_tot(1)=min(GIMAX,lf_src(isrc)%DOMAIN(2)-lf_src(isrc)%DOMAIN(1)+1)
  dimSizes_tot(2)=min(GJMAX,lf_src(isrc)%DOMAIN(4)-lf_src(isrc)%DOMAIN(3)+1)

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
  dimNames(5)='klevel'
  dimSizes_tot(3)=kmax
  dimNames_tot(3)='klevel'
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
  chunksizes=1
  chunksizes(1)=dimSizes(1)
  chunksizes(2)=dimSizes(2)
  chunksizes(3)=34 ! optimal for read and write?
  chunksizes(4)=42 ! optimal for read and write?
  chunksizes(5)=1
  chunksizes_tot=1
  chunksizes_tot(1)=MAXLIMAX
  chunksizes_tot(2)=MAXLJMAX
  chunksizes_tot(3)=dimSizes_tot(3)

  allocate(tmp_out(max(Ndiv2_coarse,Ndiv_rel*Ndiv_rel),LIMAX,LJMAX)) !NB; assumes KMAX=1 TEMPORARY
  allocate(tmp_out_cntry(LIMAX,LJMAX,(Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf))

  iou_ix = iotyp2ix(iotyp)
  if (iou_ix == iou_ix_inst) av_fac(iotyp) = 1

  !first loop only create all variables before writing into them (faster for NetCDF)
  do iter=1,2
     if(iter==1 .and. .not. first_call(iotyp))cycle

     overwrite=.false. !only used once per file
     if(iter==1)overwrite=.true.!only create all variables before writing into them
     create_var_only=.false.
     if(iter==1)create_var_only=.true.!only create all variables before writing into them

     pollwritten = .false.
     iddep = 0
     iwdep = 0
     do isrc = 1, Nsources
        if (lf_src(isrc)%species == 'FULLCHEM') cycle
        if(lf_src(Nsources)%full_chem .and. lf_src(isrc)%species/='O3' .and. lf_src(isrc)%species/='NO'  .and. lf_src(isrc)%species/='NO2'  .and. lf_src(isrc)%species/='NO3'.and. lf_src(isrc)%species/='NH3')cycle
        if (trim(lf_src(isrc)%species) == 'pm25_new') cycle !we do not output pm25_new (it is included in pm25)
        isec=lf_src(isrc)%sector
        ipoll=lf_src(isrc)%poll
        if(.not. pollwritten(ipoll))then !one pollutant may be used for several sources
           def2%name=trim(lf_src(isrc)%species)
           if(iter==1 .and. me==0)write(*,*)' poll '//trim(lf_src(isrc)%species),ipoll
           scale=1.0/av_fac(iotyp)
           call Out_netCDF(iotyp,def2,ndim_tot,kmax,lf_src_tot(1,1,KMAX_MID-lf_Nvertout+1,ipoll,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
           pollwritten(ipoll) = .true.
           overwrite=.false.
        endif

        if(iter==2)then
           fracsum=0.0
           tmp_out=0.0
           if(lf_src(isrc)%type == 'country')tmp_out_cntry=0.0
           do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                    invtot=1.0/(lf_src_tot(i,j,k,ipoll,iou_ix)+1.E-20)
                    n1=0
                    if(lf_src(isrc)%type == 'country')then
                       invfac=1.0/av_fac(iotyp) !could also output fractions?
                       do n=lf_src(isrc)%start, lf_src(isrc)%end
                          n1=n1+1
                          tmp_out_cntry(i,j,n1) = tmp_out_cntry(i,j,n1) + lf_src_acc(n,i,j,k,iou_ix)*invfac ! sum over all k
                          fracsum(i,j)=fracsum(i,j)+lf_src_acc(n,i,j,k,iou_ix)*invtot ! sum over all n and k and divided by tot
                          !if(tmp_out_cntry(i,j,n1)<1.e-18)tmp_out_cntry(i,j,n1)=0.0
!                          if(isnan(tmp_out_cntry(i,j,n1)).or. tmp_out_cntry(i,j,n1)>1.e19)then
!                             write(*,*)'tmp_out_cntry is nan ',tmp_out_cntry(i,j,n1),lf_src_acc(n,i,j,k,iou_ix),invtot,trim(lf_src(isrc)%species)
!                             stop
!                          endif
                       enddo
                    else
                       do n=lf_src(isrc)%start, lf_src(isrc)%end
                          n1=n1+1
                          tmp_out(n1,i,j) = tmp_out(n1,i,j) + lf_src_acc(n,i,j,k,iou_ix)*invtot ! sum over all k
                          fracsum(i,j)=fracsum(i,j)+lf_src_acc(n,i,j,k,iou_ix)*invtot ! sum over all n and k
                       enddo
                    endif
                 enddo
              enddo
           enddo
        endif

        if(lf_src(isrc)%type == 'country')then
           n1=0
           do i=1,Ncountry_lf+Ncountry_group_lf
              do j=1,Ncountrysectors_lf
                 n1=n1+1
                 !single cell source
                 isec=lf_country%sector_list(j)
                 if(lf_country%sector_list(j)>=0)isec=lf_country%sector_list(j)
                 if(i<=Ncountry_mask_lf)then
 !                   write(def2%name,"(A,I2.2,A5,I0)")trim(lf_src(isrc)%species)//'_sec',isec,'_mask',country_mask_val(i)
 !                   if(isec==0) write(def2%name,"(A,I0)")trim(lf_src(isrc)%species)//'_mask',country_mask_val(i)
                    write(def2%name,"(A,I2.2,A5,I0)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(mask2name(country_mask_val(i)))
                    if(isec==0) write(def2%name,"(A,I0)")trim(lf_src(isrc)%species)//'_'//trim(mask2name(country_mask_val(i)))
                 else if(i<=Ncountry_lf)then
                    write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country%list(i-Ncountry_mask_lf))
                    if(isec==0) write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country%list(i-Ncountry_mask_lf))
                 else
                    !country group
                    write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country%group(i-Ncountry_lf)%name)
                    if(isec==0) write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country%group(i-Ncountry_lf)%name)
                 endif
                 if (lf_src(Nsources)%full_chem) then
                    !add emission species to name
                    if( country_ix_list(i)==IC_STRATOS .or. country_ix_list(i)==IC_INIT)then
                       !do not add "_nox" suffix and do not output voc
                       if(EMIS_FILE(lf_src(isrc)%iem_deriv) == 'voc') cycle
                    else
                       write(def2%name,"(A)")trim(def2%name)//trim(EMIS_FILE(lf_src(isrc)%iem_deriv))
                    end if
                 end if
                 if(me==0 .and. iter==2 .and. (iotyp==IOU_MON .or. iotyp==IOU_YEAR))write(*,*)'writing '//trim(def2%name)
                 def2%unit='ug/m3'
                 scale=1.0
                 call Out_netCDF(iotyp,def2,ndim_tot,1,tmp_out_cntry(1,1,n1),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                      fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                 if(lf_src(isrc)%drydep)then
                    write(def3%name,"(A)")'DDEP_'//trim(def2%name)
                    def3%unit='mg/m2'
                    if(isrc==isrc_SO4 .or. isrc==isrc_SO2 .or. lf_src(isrc)%species=="sox")def3%unit='mgS/m2'
                    if(isrc==isrc_NH3 .or. isrc==isrc_NH4 .or. lf_src(isrc)%species=="nh3")def3%unit='mgN/m2'

                    iddep=iddep+1
                    call Out_netCDF(iotyp,def3,ndim_tot,1,loc_frac_drydep(1,1,iddep),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                         fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                 endif
                 if(lf_src(isrc)%wetdep)then
                    write(def3%name,"(A)")'WDEP_'//trim(def2%name)
                    def3%unit='mg/m2'
                    if(isrc==isrc_SO4 .or. isrc==isrc_SO2 .or. lf_src(isrc)%species=="sox")def3%unit='mgS/m2'
                    if(isrc==isrc_NH3 .or. isrc==isrc_NH4 .or. lf_src(isrc)%species=="nh3")def3%unit='mgN/m2'

                    iwdep=iwdep+1
                    call Out_netCDF(iotyp,def3,ndim_tot,1,loc_frac_wetdep(1,1,iwdep),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                         fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                 endif
                 if(lf_src(1)%MDA8 .and. lf_src(isrc)%species=='O3')then
                    write(def2%name,"(A)")"MDA8_"//trim(def2%name)
                    if(isrc==isrc_O3)then
                      call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_av(1,1,n1,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                          fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                    else
                      !voc reduction
                      call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_av(1,1,lf_src(isrc_O3)%Npos+n1,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                          fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)
                    endif

                      if(n1==1 .and. isrc==isrc_O3)then
                         !also save Base MDA8
                         def2%name="MDA8"
                         call Out_netCDF(iotyp,def2,ndim_tot,1,D8Max_av(1,1,0,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                          fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

                      end if
                 end if
              enddo
           enddo
        else
           if(old_format)then
              !for backward compatibility
              write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_local_fraction'
              if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_local_fraction'
           else
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
           endif
           scale=1.0
           call Out_netCDF(iotyp,def1,ndim,kmax,tmp_out,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes,ncFileID_given=ncFileID)
           overwrite=.false.
        endif


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

           call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)

        endif

     enddo
  enddo
  deallocate(tmp_out)

  do ipoll=1,Npoll
     do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
        do j=1,ljmax
           do i=1,limax
              lf_src_tot(i,j,k,ipoll,iou_ix) = 0.0
           enddo
        enddo
     enddo
  enddo


! reset the cumulative arrays
  do isrc = 1, Nsources
     if (lf_src(isrc)%species == 'FULLCHEM') cycle
     do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
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

  call Add_2timing(NTIMING-2,tim_after,tim_before,"lf: output")
  if(DEBUGall .and. me==0)write(*,*)'end out'

! CALL MPI_BARRIER(MPI_COMM_CALC, I)

!stop
end subroutine lf_out

subroutine lf_av(dt,End_of_Day)
  real, intent(in)    :: dt                   ! time-step used in integrations
  logical, intent(in) :: End_of_Day           ! e.g. 6am for EMEP sites
  real :: xtot, x, O3_c
  integer ::i,j,k,ii,l,n,n_new,dx,dy,ix,iix,ipoll,isec_poll1, iou_ix, isrc
  integer ::isec_poll
  logical :: pollwritten(Max_lf_spec)
  integer,save :: count_AvgMDA8_m=0,count_AvgMDA8_y=0
  real :: w_m, w_y, timefrac
  logical, save :: first_call=.true.

  if(DEBUGall .and. me==0)write(*,*)'start av'

  call Code_timer(tim_before)
  if(.not. lf_src(1)%HOUR.and.&
     .not. lf_src(1)%DAY .and.&
     .not. lf_src(1)%MONTH .and.&
     .not. lf_src(1)%YEAR       )return

  !do the averaging
  do iou_ix = 1, Niou_ix
     pollwritten = .false.
     do isrc=1,Nsources_nonew
        if (lf_src(isrc)%species == 'FULLCHEM') cycle
        ipoll = lf_src(isrc)%poll
        do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
           do j=1,ljmax
              do i=1,limax
                 xtot=0.0
                 do iix=1,lf_src(isrc)%Nsplit
                    ix = lf_src(isrc)%ix(iix)
                     if(lf_src(isrc)%type=='country')then
                       !3m height cfac correction
                       xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix))/ATWAIR&
                         *roa(i,j,k,1)*1.E9* cfac(ix,i,j) !for ug/m3
                       !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                     else
                         xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix))/ATWAIR&
                              *roa(i,j,k,1)*1.E9 !for ug/m3
                        !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                     endif
                 end do
                 if(.not. pollwritten(ipoll))then !one pollutant may be used for several sources
                    if (iou_ix == iou_ix_inst) then
                       !not accumulated
                       lf_src_tot(i,j,k,ipoll,iou_ix) = xtot
                    else
                       lf_src_tot(i,j,k,ipoll,iou_ix) = lf_src_tot(i,j,k,ipoll,iou_ix) + xtot
                    endif
                 endif
                 do n=lf_src(isrc)%start, lf_src(isrc)%end
                    if (iou_ix == iou_ix_inst) then
                       !not accumulated
                       lf_src_acc(n,i,j,k,iou_ix) = xtot*lf(n,i,j,k)
                    else
                       lf_src_acc(n,i,j,k,iou_ix)=lf_src_acc(n,i,j,k,iou_ix)+xtot*lf(n,i,j,k)
                    end if
                    if(DEBUG .and. isnan(lf_src_acc(n,i,j,k,iou_ix)))then
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
              do j=1,ljmax
                 do i=1,limax
                    xtot=0.0
                    do iix=1,lf_src(isrc_pm25_new)%Nsplit
                       ix=lf_src(isrc_pm25_new)%ix(iix)
                       if(lf_src(isrc_pm25_new)%type=='country')then
                          !3m height cfac correction
                          xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc_pm25_new)%mw(iix))/ATWAIR&
                               *roa(i,j,k,1)*1.E9* cfac(ix,i,j) !for ug/m3
                          !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                       else
                          xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc_pm25_new)%mw(iix))/ATWAIR&
                               *roa(i,j,k,1)*1.E9 !for ug/m3
                          !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                       endif
                    end do
                    if(.not. pollwritten(ipoll))then !one pollutant may be used for several sources
                       !It is added to the _old (= not new) component, also in the inst case
                       lf_src_tot(i,j,k,ipoll,iou_ix) = lf_src_tot(i,j,k,ipoll,iou_ix) + xtot
                    end if
                    n_new = lf_src(isrc_pm25_new)%start
                    do n=lf_src(isrc)%start, lf_src(isrc)%end !NB: loop over isrc for pm25, not new
                       lf_src_acc(n,i,j,k,iou_ix)=lf_src_acc(n,i,j,k,iou_ix) + xtot*lf(n_new,i,j,k)
                       n_new = n_new + 1
                    end do
                 end do
              enddo
           enddo
        end if
        pollwritten(ipoll) = .true.

     end do
        
     if (lf_src(1)%MDA8) then! maximum daily eight-hour mean concentration
        
        if (iou_ix == 1) then !only once for (daily, monthly and yearly)
           do j = 1,ljmax
              do i = 1,limax
                 ix=O3_ix-NSPEC_SHL
                 O3_c = xn_adv(ix,i,j,KMAX_MID) * cfac(ix,i,j) * roa(i,j,KMAX_MID,1) * 1.0e9 * species_adv(ix)%molwt/ATWAIR
                 hourM(0,i,j) = hourM(0,i,j) + O3_c
                 !how much the hourly values will change for a small change of emissions * 100% (summed over hour, normalize later)
                 do n=1, lf_src(isrc_O3)%Npos
                    hourM(n,i,j) = hourM(n,i,j) + O3_c * lf(lf_src(isrc_O3)%start+n-1,i,j,KMAX_MID)
                 end do
                 do n=1, lf_src(isrc_O3_voc)%Npos
                    hourM(lf_src(isrc_O3)%Npos+n,i,j) = hourM(lf_src(isrc_O3)%Npos+n,i,j) + O3_c * lf(lf_src(isrc_O3_voc)%start+n-1,i,j,KMAX_MID)
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
              !one hour has past since last time here        
              ! update max value since 01:00
              do j = 1,ljmax
                 do i = 1,limax
                    x = 0.0
                    do l = 1, 8
                       x = x + D8M(0,i,j,l) * 0.125 ! 8 hour average
                    end do
                    if(x > D8Max(0,i,j)) then
                       !a new max is found. Update for all fractions too
                       do n=0, lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos
                          D8Max(n,i,j) = D8M(n,i,j,1) * 0.125 ! 8 hour average, contribution from first hour
                       end do
                       do l = 2, 8
                          do n=0, lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos
                             D8Max(n,i,j) = D8Max(n,i,j) + D8M(n,i,j,l) * 0.125 ! 8 hour average, contribution from second to eighth hour
                          end do
                       end do
                    end if
                 end do
              end do
           end if
        end if
        if (current_date%seconds == 0 .and. current_date%hour == 0 .and. .not. first_call) then
           !end of day, save the values
           !NB: at the end of the first day (day 2 hour 00:00), we actually start to write in the next month
           if (current_date%day == 2 .and. iotyp2ix(iou_ix)==IOU_MON) then
              !new month
              count_AvgMDA8_m = 0
              D8Max_av(:,:,:,iou_ix)=0.0
           end if
           if (current_date%day == 2 .and. current_date%month == 4 .and. iotyp2ix(iou_ix)==IOU_YEAR) then
              !new yearly max
              count_AvgMDA8_y = 0
              D8Max_av(:,:,:,iou_ix)=0.0
           end if

           if(iotyp2ix(iou_ix)==IOU_MON)count_AvgMDA8_m = count_AvgMDA8_m + 1
           if(iotyp2ix(iou_ix)==IOU_YEAR)count_AvgMDA8_y = count_AvgMDA8_y + 1
           w_m = 1.0/count_AvgMDA8_m
           w_y = 1.0/count_AvgMDA8_y
           do j = 1,ljmax
              do i = 1,limax
                 if (iotyp2ix(iou_ix)==IOU_DAY)then
                    do n=0, lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos
                       D8Max_av(i,j,n,iou_ix) =  D8Max(n,i,j)
                    end do
                 else if(iotyp2ix(iou_ix)==IOU_MON)then
                    do n=0, lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos
                       D8Max_av(i,j,n,iou_ix) =  (1.0-w_m) * D8Max_av(i,j,n,iou_ix) + w_m * D8Max(n,i,j)
                    end do
                 else if (iotyp2ix(iou_ix)==IOU_YEAR)then
                    !we keep only days from April to September
                    if(current_date%month>=4 .and. current_date%month<=9)then
                       do n=0, lf_src(isrc_O3)%Npos+lf_src(isrc_O3_voc)%Npos
                          D8Max_av(i,j,n,iou_ix) =  (1.0-w_y) * D8Max_av(i,j,n,iou_ix) + w_y * D8Max(n,i,j)
                       end do
                    end if
                 end if
              end do
           end do
           if(iou_ix == Niou_ix) D8Max = 0.0 !we are ready to start a new day
        end if
     end if
  end do

  av_fac=av_fac+1
  first_call=.false.
  
  call Add_2timing(NTIMING-9,tim_after,tim_before,"lf: averaging")
  if(DEBUGall .and. me==0)write(*,*)'end av'

end subroutine lf_av

subroutine lf_adv_x(fluxx,i,j,k)
  real, intent(in)::fluxx(NSPEC_ADV,-1:LIMAX+1)
  integer, intent(in)::i,j,k
  real ::x,xn,xx,f_in,inv_tot
  integer ::n,ii,iix,ix,dx,dy,isrc,dp,dm
  if(DEBUGall .and. me==0)write(*,*)'start advx'

  if(i==li0)then
     !copy small part (could be avoided, but simpler to copy)
     !note that the right hand side of the lf equations must contain unupdated values, therefore values for j-1 must be buffered
     do ii=li0,li1
        do n=1,LF_SRC_TOTSIZE
           loc_frac_src_1d(n,ii) = lf(n,ii,j,k)
        enddo
     enddo
  endif

  call Code_timer(tim_before)
  do isrc=1,Nsources
     if (lf_src(isrc)%species == 'FULLCHEM') cycle
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
     if(lf_src(isrc)%type=='coarse' .or. lf_src(isrc)%type=='country')then
        if(x>1.E-20)then
           do n = lf_src(isrc)%start, lf_src(isrc)%end
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,i+1)*x
           enddo
           if(xx>1.E-20)then
              do n = lf_src(isrc)%start, lf_src(isrc)%end
                 lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n,i-1)*xx
              enddo
           endif
        else if (xx>1.E-20)then
           do n = lf_src(isrc)%start, lf_src(isrc)%end
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,i-1)*xx
           enddo
        endif
     else if(lf_src(isrc)%type=='relative')then
        !The relative position of the source is either the same for i and incoming lf, or differs by 1
        dp = 0 ! for left boundary
        dm = 0 ! for right boundary
        if(mod(i_fdom(i)-1,lf_src(isrc)%res)==0)dp=1
        if(mod(i_fdom(i),lf_src(isrc)%res)==0)dm=1
        if(x>1.E-20)then
           n = lf_src(isrc)%start
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
           endif

           if(xx>1.E-20)then
              n = lf_src(isrc)%start
              do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                 do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                    lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+dp,i-1)*xx
                    n=n+1
                 enddo
                 if (dp==0) lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+dp,i-1)*xx
                 n=n+1! when dx=lf_src(isrc)%dist and dp=1 there are no local fractions to transport
              enddo
           endif
        else if (xx>1.E-20)then
           n = lf_src(isrc)%start
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
        endif
      else
        if(me==0)write(*,*)'LF type not recognized)'
        stop
     endif
  enddo

  call Add_2timing(NTIMING-8,tim_after,tim_before,"lf: adv_x")

  if(DEBUGall .and. me==0)write(*,*)'end advx'
end subroutine lf_adv_x

subroutine lf_adv_y(fluxy,i,j,k)
  real, intent(in)::fluxy(NSPEC_ADV,-1:LJMAX+1)
  integer, intent(in)::i,j,k
  real ::x,xn,xx,f_in,inv_tot
  integer ::n,jj,iix,ix,dx,dy,isrc,dp,dm
  if(DEBUGall .and. me==0)write(*,*)'start advy'

  if(j==lj0)then
     !copy small part (could be avoided, but simpler to copy)
     !note that the right hand side of the lf equations must contain unupdated values, therefore values for i-1 must be buffered
     do jj=lj0,lj1
        do n=1,LF_SRC_TOTSIZE
           loc_frac_src_1d(n,jj) = lf(n,i,jj,k)
        enddo
     enddo
  endif

  call Code_timer(tim_before)
  do isrc=1,Nsources
     if (lf_src(isrc)%species == 'FULLCHEM') cycle
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
     if(lf_src(isrc)%type=='coarse' .or. lf_src(isrc)%type=='country')then
        !often either x or xx is zero
        if(x>1.E-20)then
           do n = lf_src(isrc)%start, lf_src(isrc)%end
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,j+1)*x
           enddo
           if(xx>1.E-20)then
              do n = lf_src(isrc)%start, lf_src(isrc)%end
                 lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n,j-1)*xx
              enddo
           endif
        else if (xx>1.E-20)then
           do n = lf_src(isrc)%start, lf_src(isrc)%end
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,j-1)*xx
           enddo
        endif
     else if(lf_src(isrc)%type=='relative')then
        !The relative position of the source is either the same for j and incoming lf, or differs by 1 (n differs by Ndiv_rel)
        dp = 0 ! for lower boundary
        dm = 0 ! for upper boundary
        if(mod(j_fdom(j)-1,lf_src(isrc)%res)==0)dp=Ndiv_rel
        if(mod(j_fdom(j),lf_src(isrc)%res)==0)dm=Ndiv_rel
        if(x>1.E-20)then
           n = lf_src(isrc)%start
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
           endif
           if(xx>1.E-20)then
              n = lf_src(isrc)%start
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
              endif
           endif
        else if (xx>1.E-20)then
           n = lf_src(isrc)%start
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
           endif
        else
           !nothing to do if no incoming fluxes
        endif
     else
        if(me==0)write(*,*)'LF type not recognized)'
        stop
     endif

  enddo
  call Add_2timing(NTIMING-7,tim_after,tim_before,"lf: adv_y")
  if(DEBUGall .and. me==0)write(*,*)'end advy'

end subroutine lf_adv_y

subroutine lf_adv_k(fluxk,i,j)
    real, intent(in)::fluxk(NSPEC_ADV,KMAX_MID)
    integer, intent(in)::i,j
    real ::x,xn,xx,f_in,inv_tot
    integer ::n,k,iix,ix,dx,dy,isrc
    real loc_frac_src_km1(LF_SRC_TOTSIZE,KMAX_MID-lf_Nvert+1:KMAX_MID)
  if(DEBUGall .and. me==0)write(*,*)'start advk'

    call Code_timer(tim_before)
    !need to be careful to always use non-updated values on the RHS
    do k = KMAX_MID-lf_Nvert+2,KMAX_MID
       do n = 1, LF_SRC_TOTSIZE
          loc_frac_src_km1(n,k)=lf(n,i,j,k-1) !NB: k is shifted by 1 in loc_frac_src_km1
       enddo
    enddo
    loc_frac_src_km1(:,KMAX_MID-lf_Nvert+1) = 0.0 ! everything above is not tracked, zero local fractions coming from above

    do n = 1, nstratos
       loc_frac_src_km1(Stratos_ix(n),KMAX_MID-lf_Nvert+1) = 1.0 ! NB: k is shifted by 1 in loc_frac_src_km1, i.e. this is level k=1 if lf_Nvert = KMAX_MID - 1
    end do

    do k = KMAX_MID-lf_Nvert+1,KMAX_MID!k is increasing-> can use k+1 to access non-updated value
       do isrc=1,Nsources
          if (lf_src(isrc)%species == 'FULLCHEM') cycle
          xn=0.0
          x=0.0
          xx=0.0
          !positive x or xx means incoming, negative means outgoing
          do iix=1,lf_src(isrc)%Nsplit
             ix=lf_src(isrc)%ix(iix)
             xn=xn+xn_2d(ix,k)*lf_src(isrc)%mw(iix) !xn_2d are the  unupdated values
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
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_km1(n,k)*xx !NB: lf already multiplied by xn here
                enddo
             endif
          else if (xx>1.E-20)then
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_km1(n,k)*xx
             enddo
          else
             !nothing to do if no incoming fluxes
          endif
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
    real loc_frac_src(LF_SRC_TOTSIZE,KMAX_MID-lf_Nvert-1:KMAX_MID+2)
    real :: w1,w2,w3 !weights from the used level for making the flux
    !NB: level used for *incoming* fluxes:
    !    fluxes incoming k from k-1: uses level k-2,k-1 and k to make flux(k)
    !    fluxes incoming k from k+1: uses level k,k+1 and k+2 to make flux(k+1)
    ! Eta=A/Pref+B , P = A + B*PS, high altitude means A and B and Eta decrease. Positive Etadot means downwind.
    real fc(KMAX_MID),zzfl1(KMAX_MID-lf_Nvert:KMAX_MID),zzfl2(KMAX_MID-lf_Nvert:KMAX_MID),zzfl3(KMAX_MID-lf_Nvert:KMAX_MID)
    real :: x1,x2,x3,xx1,xx2,xx3,fk1,lfmax

    if(DEBUGall .and. me==0)write(*,*)'start adv_k_2nd'

    call Code_timer(tim_before)

    lfmax=10.0 !limitation for extreme values
    !note about extreme values:
    ! in some cases, Bott will almost entirely empty a cell. If it is 99.999% or 99.9 % may be dependent on a neighboring gridcells concentration.
    ! Since the final concentration is almost zero, the derivative divided by that amount can be large.


    !need to be careful to always use non-updated values on the RHS
    do k = KMAX_MID-lf_Nvert-1,KMAX_MID+2
       if(k>=KMAX_MID-lf_Nvert+1 .and. k<=KMAX_MID)then
          do n = 1, LF_SRC_TOTSIZE
             loc_frac_src(n,k)=lf(n,i,j,k)
          enddo
       else
          do n = 1, LF_SRC_TOTSIZE
             loc_frac_src(n,k)=0.0
          enddo
       end if
    enddo

    !for stratosphere tracking: everything above the tracking region has lf=1
    do n = 1, nstratos
       loc_frac_src(Stratos_ix(n),KMAX_MID-lf_Nvert) = 1.0 ! this is level k=1 if lf_Nvert = KMAX_MID - 1
    end do

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

    do isrc=1,Nsources
       if (lf_src(isrc)%species == 'FULLCHEM') cycle
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
             n1k=1
             !fluxk(:,k+1) is made from xn(k), xn(k+1), xn(k+2)
             !fluxk(:,k+1)<0 into k from k+1
             do iix=1,lf_src(isrc)%Nsplit
                ix=lf_src(isrc)%ix(iix)
                xn=xn+xn_2d(ix,k)*lf_src(isrc)%mw(iix) !xn_2d are the  unupdated values
                x=x-fluxk(ix,k+1)*lf_src(isrc)%mw(iix)
                x1=x1+zzfl1(k)*xn_2d(ix,k)*lf_src(isrc)%mw(iix)
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
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                lf(n,i,j,k+1) = -(x1*loc_frac_src(n,k)+x2*loc_frac_src(n,k+1)+x3*loc_frac_src(n,k+2))
             enddo
             xn_post=xn+(fk1+(x1+x2+x3))*dhs1i(k+1)
             fk1=-(x1+x2+x3)

             if (k>KMAX_MID-lf_Nvert) then
                !into k  (already initialized. Divide by final concentration)
                if (xn_post>1e-20) then
                   do n = lf_src(isrc)%start, lf_src(isrc)%end
                      lf(n,i,j,k) = (xn*loc_frac_src(n,k)+(lf(n,i,j,k)+x1*loc_frac_src(n,k)+x2*loc_frac_src(n,k+1)+x3*loc_frac_src(n,k+2))*dhs1i(k+1))/xn_post
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
             n1k=0
             !fluxk(:,k+1) is made from xn(k-1), xn(k), xn(k+1)
             !fluxk(:,k+1)>0 out of k into k+1
             do iix=1,lf_src(isrc)%Nsplit
                ix=lf_src(isrc)%ix(iix)
                xn=xn+xn_2d(ix,k)*lf_src(isrc)%mw(iix) !xn_2d are the  unupdated values
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
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                lf(n,i,j,k+1) =x1*loc_frac_src(n,k-1)+x2*loc_frac_src(n,k)+x3*loc_frac_src(n,k+1)
             enddo

             xn_post=xn+(fk1-(x1+x2+x3))*dhs1i(k+1)
             fk1=(x1+x2+x3)
             if (k>KMAX_MID-lf_Nvert) then
                !out of k (already initialized. Divide by final concentration)
                if (xn_post>1e-20) then
                   do n = lf_src(isrc)%start, lf_src(isrc)%end
                      lf(n,i,j,k) =(xn*loc_frac_src(n,k)+(lf(n,i,j,k)-(x1*loc_frac_src(n,k-1)+x2*loc_frac_src(n,k)+x3*loc_frac_src(n,k+1)))*dhs1i(k+1))/xn_post
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
                xn=xn+xn_2d(ix,k+1)*lf_src(isrc)%mw(iix) !xn_2d are the  unupdated values
             end do
             xn_post=xn+fk1*dhs1i(k+2)
             if (xn_post>1e-20) then
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   lf(n,i,j,k+1) =(xn*loc_frac_src(n,k+1)+lf(n,i,j,k+1)*dhs1i(k+2))/xn_post
                   lf(n,i,j,k+1) =max(-lfmax,min(lfmax,lf(n,i,j,k+1)))
                enddo
             else
                !if there are no concentration, we set the fractions to zero
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   lf(n,i,j,k+1) = 0.0
                enddo
             end if
          endif

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
    integer ::isec_poll1,isrc
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
          if (lf_src(isrc)%species == 'FULLCHEM') cycle
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
          endif
          xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k) = x
      enddo
    enddo

    call vertdiffn(xn_k,LF_SRC_TOTSIZE+Npoll,1,KMAX_MID-lf_Nvert-KUP,EtaKz(i,j,1,1),dt_advec,ndiff) !DSKz

    do k = KMAX_MID-lf_Nvert+1,KMAX_MID
       do isrc=1,Nsources
          if (lf_src(isrc)%species == 'FULLCHEM') cycle
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
    real :: xn_k(LF_SRC_TOTSIZE + Npoll,KMAX_MID),x
    integer ::isec_poll1,isrc
    integer ::k,n,ix,iix,dx,dy

    if(DEBUGall .and. me==0)write(*,*)'start conv'
    call Code_timer(tim_before)
    xn_k = 0.0
    do k = 1,KMAX_MID
       do isrc=1,Nsources
          if (lf_src(isrc)%species == 'FULLCHEM') cycle
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
          endif
          xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k) = x
      enddo
    enddo

    call convection_1d(xn_k,LF_SRC_TOTSIZE+Npoll,1,i,j,dt_conv)

    do k = KMAX_MID-lf_Nvert+1,KMAX_MID
       do isrc=1,Nsources
          if (lf_src(isrc)%species == 'FULLCHEM') cycle
          x =  1.0/(xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k)+1.E-30)
          do n=lf_src(isrc)%start, lf_src(isrc)%end
             lf(n,i,j,k) = xn_k(n,k)*x
          enddo
       enddo
    end do
    call Add_2timing(NTIMING-5,tim_after,tim_before,"lf: diffconv")
    if(DEBUGall .and. me==0)write(*,*)'end conv'

end subroutine lf_conv

subroutine lf_chem_emis_deriv(i,j,k,xn,xnew,eps1)
  real, intent(in) :: eps1,xn(NSPEC_TOT),xnew(NSPEC_TOT)
  integer, intent(in) :: i,j,k
  integer :: n, n0, ispec, nispec, isrc,isrc_emis, ic, is, ics, ncs, ideriv0, found
  real :: efac, xtot,totemis,emiss
  real :: xd(NSPEC_fullchem_lf,NSPEC_fullchem_lf)
  integer :: isec, iem, iix, ix, iiix,iemis, ideriv , iem_deriv, isrc_deriv


  if (i<li0 .or.i>li1 .or.j<lj0.or.j>lj1)return !we avoid outer frame
  if(DEBUGall .and. me==0)write(*,*)'start chememis'

  if (nfullchem <= 0) then
     !case with no chemistry for local fractions
     !Now species may be group of species (pm25...)

     if(k<max(KEMISTOP,KMAX_MID-lf_Nvert+1))return
     !include emissions that are not created in chemical reactions, but only by emissions
     call Code_timer(tim_before)
     do isrc=1,Nsources
        if (lf_src(isrc)%full_chem) cycle !included below
        if (lf_src(isrc)%species == 'FULLCHEM') cycle
        iem = lf_src(isrc)%iem     !pick out only this emissions
        isec = lf_src(isrc)%sector !pick out only this sector

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
           do ic=1,Ncountry_lf+Ncountry_group_lf
              do is=1,Ncountrysectors_lf
                 lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot)/(xtot+totemis+1.e-20)
                 n0=n0+1
              end do
           end do
           !add emissions that are tracked
           do iemis = 1, N_lf_derivemis
              if(emis2isrc(iemis) /= isrc ) cycle

              n0 = lf_src(isrc)%start + emis2icis(iemis)

              lf(n0,i,j,k)=lf(n0,i,j,k) + rcemis_lf(iemis,1)/(xtot+totemis+1.e-20)

           end do

        else if(lf_src(isrc)%type=='relative')then
           !first dilute lf because of emissions
           do n0=lf_src(isrc)%start,lf_src(isrc)%end
              lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot)/(xtot+totemis+1.e-20)
           end do
           !add emissions only at centre of window
           n0 = lf_src(isrc)%start + (lf_src(isrc)%Npos - 1)/2 !"middle" point is dx=0 dy=0
           do iemis = 1, N_lf_derivemis
              if(emis2isrc(iemis) /= isrc ) cycle
              lf(n0,i,j,k)=lf(n0,i,j,k) + rcemis_lf(iemis,1)/(xtot+totemis+1.e-20)
           end do

        else
           if(me==0)write(*,*)'LF type not recognized)'
           stop
        endif

     enddo
     call Add_2timing(NTIMING-4,tim_after,tim_before,"lf: emissions")

  else
     !case using Jacobian of chemical reactions and emissions are considered part of chemical reactions
  call Code_timer(tim_before)

  !make derivatives
  do ispec = 1, NSPEC_fullchem_lf!no loop over short lived

     do ideriv = 1, NSPEC_fullchem_lf !loop over chemical derivatives
        if (xn(ideriv+NSPEC_SHL)>1000.0 .and. xnew(ispec+NSPEC_SHL)>1000.0 .and. ispec+NSPEC_SHL/=RO2POOL_ix .and. ideriv/=RO2POOL_ix-NSPEC_SHL) then !units are molecules/cm3 means almost no molecules
           ! derivatives are normalized to concentration of each species, dx/dy *y/x
           !  y and dy must have same units
           !(xnew_lf(ideriv,ispec) - xnew(ispec))/(0.0001*xn(ideriv+NSPEC_SHL)) * xn(ideriv+NSPEC_SHL)/xnew(ispec) =
           xderiv(ideriv,ispec) = (xnew_lf(ideriv,ispec+NSPEC_SHL) - xnew(ispec+NSPEC_SHL))/((eps1-1.0)*(xnew(ispec+NSPEC_SHL)))

           !remove numerical noise
           xderiv(ideriv,ispec) = min(xderiv(ideriv,ispec),10.0)
           xderiv(ideriv,ispec) = max(xderiv(ideriv,ispec),-10.0)
        else
           !Not that concentrations reaching or leaving 0.0 are not treated correctly
           !remove numerical noise
           xderiv(ideriv,ispec) = 0.0
        end if
     end do

     do ideriv = 1, N_lf_derivemis !loop over emission derivatives
        if (xnew(ispec+NSPEC_SHL)>1000.0 .and. xnew_lf(ideriv+ NSPEC_fullchem_lf,ispec+NSPEC_SHL)>1000.0 .and. k>=KEMISTOP.and. ispec+NSPEC_SHL/=RO2POOL_ix ) then !units are molecules/cm3 means almost no molecules
           ! derivatives are normalized to concentration of each species, dx/dy *y/x
           !  y and dy must have same units
           !(xnew_lf(ideriv,ispec) - xnew(ispec))/(0.0001*xn(ideriv+NSPEC_SHL)) * xn(ideriv+NSPEC_SHL)/xnew(ispec) =
           ederiv(ideriv,ispec) = (xnew_lf(ideriv + NSPEC_fullchem_lf ,ispec+NSPEC_SHL) - xnew(ispec+NSPEC_SHL))/((eps1-1.0)*(xnew(ispec+NSPEC_SHL)))

           !remove numerical noise
           ederiv(ideriv,ispec) = min(ederiv(ideriv,ispec),10.0)
           ederiv(ideriv,ispec) = max(ederiv(ideriv,ispec),-10.0)

        else
           !Note that concentrations reaching 0.0 are not treated correctly
           ederiv(ideriv,ispec) = 0.0
        end if
     end do
  end do

  !remove RO2POOL which is not a proper species
  if (RO2POOL_ix>0) then
     xderiv(:,RO2POOL_ix-NSPEC_SHL) = 0.0
     ederiv(1:N_lf_derivemis,RO2POOL_ix-NSPEC_SHL) = 0.0
     xderiv(RO2POOL_ix-NSPEC_SHL,:) = 0.0
  end if

  ncs = (Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf

  ! emissions derivatives taken into account for sensibility to changes of emis *during* chemical process

  do isrc = 1, Nsources
     ix = lf_src(isrc)%ix(1) ! index of species in species_adv()
     n0 = lf_src(isrc)%start
     n=1
     if (isrc > NSPEC_fullchem_lf) n = ncs + 1
     !init
     do ics = n0,  n0 + ncs - 1 !loop over sector and country contributions
        lf0_loc(n, ix) = lf(ics,i,j,k)
        lf(ics,i,j,k) = 0.0
        n=n+1
     end do

     !include emission derivatives
     do iemis = 1, N_lf_derivemis
        if(emis2isrc(iemis) == lf_src(isrc)%iem_deriv ) then
           n0 =lf_src(isrc)%start + emis2icis(iemis)
           lf(n0,i,j,k) = ederiv(iemis,ix) !contribution from emissions during this timestep. NB: only one iemis per n0 can be included
        end if
     end do
  end do

  ! concentrations have changed in the chemical reactions.
  ! The derivatives give the sensitivity to each species
  ! The local fraction are then propagated with the derivatives as weights

  lf_loc = matmul(lf0_loc, xderiv) !MAIN STEP

  !Note that if a species has emissions and no chemistry, xderiv(x,x) = x(t)/(x(t)+E(dt)) /= 1

  do isrc = 1, Nsources
     ix = lf_src(isrc)%ix(1) ! index of species in species_adv()
     n0 = lf_src(isrc)%start
     n=1
     if (isrc > NSPEC_fullchem_lf) n = ncs + 1
     do ics = n0,  n0 + ncs - 1 !loop over sector and country contributions
        lf(ics,i,j,k) = lf(ics,i,j,k) + lf_loc(n, ix)
        n=n+1
     end do
  end do

  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")
  if(DEBUGall .and. me==0)write(*,*)'end chememis'
return
  if(i_fdom(i)==108.and.j_fdom(j)==95 .and.k==kmax_mid )then
     write(*,*)'concentration after chem ',xnew(18)
     write(*,*)'prediction after 0.001 nox emis change',xnew(18)*(1.0+0.001*xderiv(NSPEC_fullchem_lf+1,18))
     write(*,*)'lf O3 ',lf(3,i,j,k)
     write(*,*)'prediction after chem and 0.001 nox emis change',xnew(18)*(1+0.001*lf(3,i,j,k))
  end if
  return
  67 format(10x,50(A11))
  if(i_fdom(i)==58.and.j_fdom(j)==36 .and.k==kmax_mid )then
     68 format(10x,50E11.4)
     write(*,68)(xn(n+NSPEC_SHL),n = 1, min(7,NSPEC_fullchem_lf))
     write(*,68)(xnew(n+NSPEC_SHL),n = 1, min(7,NSPEC_fullchem_lf))
     write(*,67)(trim(species_adv(n)%name)//' ',n = 1, min(7,NSPEC_fullchem_lf)),'  NOXEMIS'

     do ispec = 1, NSPEC_SHL+NSPEC_fullchem_lf
66      format(A10,50(F10.5,A1))
        write(*,66)species(ispec)%name//' ',(xderiv(n,ispec),' ',n = 1, min(7,NSPEC_fullchem_lf)),xderiv(NSPEC_fullchem_lf+1,ispec)
     end do
        write(*,66)'sum             ',(sum(xderiv(n,1:NSPEC_SHL+NSPEC_fullchem_lf)),' ',n = 1, min(7,NSPEC_fullchem_lf))

  endif

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

  ageing_rate = EC_AGEING_RATE()
  if(DEBUGall .and. me==0)write(*,*)'start chem'

  if (isrc_pm25 > 0) then
     do isrc = 1, Nsources_nonew
        isrc_pm25_new = isrc_new(isrc)
        if (isrc_pm25_new < 0) cycle
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

  if(isrc_SO2>0)then
     do k = KMAX_MID-lf_Nvert+1,KMAX_MID
        SO4 = xn_2d(SO4_ix,k)
        n_SO4 = lf_src(isrc_SO4)%start
        write(*,*)'LF:SO2 not implemented stop'
        stop
        !SO4 produced by SO2 , without emitted SO4:
        !d_SO4 = max(0.0,Dchem(NSPEC_SHL+ix_SO4,k,i,j)-rcemis(NSPEC_SHL+ix_SO4,k))*dt_advec
        inv = 1.0/(SO4 + 1.0E-20)

        do n_SO2=lf_src(isrc_SO2)%start, lf_src(isrc_SO2)%end
           lf(n_SO4,i,j,k) = (lf(n_SO4,i,j,k)*(SO4-d_SO4)+d_SO4*lf(n_SO2,i,j,k)) * inv
           n_SO4 = n_SO4 + 1
        enddo

    end do

 end if

  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")
  if(DEBUGall .and. me==0)write(*,*)'end chem'
end subroutine lf_chem

subroutine lf_aero_pre(i,j) !called just before AerosolEquilib
  integer, intent(in) ::i,j
  integer :: k
  !save concentrations, to see changes
  call Code_timer(tim_before)
  if (lf_fullchem .and. NO3_ix>0 .and. HNO3_ix>0) then
     do k = KMAX_MID-lf_Nvert+1,KMAX_MID
        lf_NO3(k) = xn_2d(NO3_ix,k)
        lf_HNO3(k) = xn_2d(HNO3_ix,k)
     enddo
  end if
  if (isrc_NH4<0) return;
  do k = KMAX_MID-lf_Nvert+1,KMAX_MID
     lf_NH4(k) = xn_2d(NH4_f_ix,k)
     lf_NH3(k) = xn_2d(NH3_ix,k)
  enddo
call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

end subroutine lf_aero_pre

subroutine lf_aero_pos (i,j) !called just after AerosolEquilib
  integer, intent(in) ::i,j
  real :: d_NH4, d_NH3, NH4, NH3, inv
  integer :: n_NH3, n_NH4, k
  real :: d_NO3, d_HNO3, NO3, HNO3
  integer :: n_NO3, n_HNO3

  if (isrc_NH4<0 .and. .not.lf_fullchem) return;
  call Code_timer(tim_before)
  do k = KMAX_MID-lf_Nvert+1,KMAX_MID

     NH3 = xn_2d(NH3_ix,k)
     NH4 = xn_2d(NH4_f_ix,k)
     if (isrc_NH4>0)then
     d_NH4 = NH4 - lf_NH4(k)
     d_NH3 = NH3 - lf_NH3(k)

     if(d_NH4>0.0 .and. d_NH3<0.0)then
        !NH3 has been transformed into NH4
        n_NH3 = lf_src(isrc_NH3)%start
        inv = 1.0/(NH4+d_NH4)
        do n_NH4=lf_src(isrc_NH4)%start, lf_src(isrc_NH4)%end
           lf(n_NH4,i,j,k) = (lf(n_NH4,i,j,k)*NH4 + d_NH4*lf(n_NH3,i,j,k)) * inv
           n_NH3 = n_NH3 + 1
        enddo
     else if(d_NH4<0.0 .and. d_NH3>0.0)then
        !NH4 has been transformed into NH3
        n_NH3 = lf_src(isrc_NH3)%start
        inv = 1.0/(NH3+d_NH3)
        do n_NH4=lf_src(isrc_NH4)%start, lf_src(isrc_NH4)%end
           lf(n_NH3,i,j,k) = (lf(n_NH3,i,j,k)*NH3 + d_NH3*lf(n_NH4,i,j,k)) * inv
           n_NH3 = n_NH3 + 1
        enddo
     else
        !N is not conserved or concentrations are constant
     endif
     endif
     if (lf_fullchem) then
        NO3 = xn_2d(NO3_ix,k)
        HNO3 = xn_2d(HNO3_ix,k)
        d_NO3 = NO3 - lf_NO3(k)
        d_HNO3 = HNO3 - lf_HNO3(k)

        if(d_NO3>0.0 .and. d_HNO3<0.0)then
           !HNO3 has been transformed into NO3
           n_HNO3 = lf_src(isrc_HNO3)%start
           inv = 1.0/(NO3+d_NO3)
           do n_NO3=lf_src(isrc_NO3)%start, lf_src(isrc_NO3)%end
              lf(n_NO3,i,j,k) = (lf(n_NO3,i,j,k)*NO3 + d_NO3*lf(n_HNO3,i,j,k)) * inv
              n_HNO3 = n_HNO3 + 1
           enddo
        else if(d_NO3<0.0 .and. d_HNO3>0.0)then
           !NO3 has been transformed into HNO3
           n_HNO3 = lf_src(isrc_HNO3)%start
           inv = 1.0/(HNO3+d_HNO3)
           do n_NO3=lf_src(isrc_NO3)%start, lf_src(isrc_NO3)%end
              lf(n_HNO3,i,j,k) = (lf(n_HNO3,i,j,k)*HNO3 + d_HNO3*lf(n_NO3,i,j,k)) * inv
              n_HNO3 = n_HNO3 + 1
           enddo
        else
           !N is not conserved or concentrations are constant
        endif

     end if
  enddo
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

end subroutine lf_aero_pos

subroutine  lf_drydep(i,j,DepLoss, fac)
  integer, intent(in) :: i,j
  real, intent(in) :: fac
  real, intent(in), dimension(NSPEC_ADV) :: DepLoss
  integer :: n,ix,iix,idep, idep0, isrc
  real :: ffac
  integer :: istart,iend
  idep0=0
  idep=0
  call Code_timer(tim_before)

  if(DEBUGall .and. me==0)write(*,*)'start drydep'
  do isrc=1,Nsources
     if(.not. lf_src(isrc)%DryDep)cycle
     if (lf_src(isrc)%species == 'FULLCHEM') cycle
     do iix=1,lf_src(isrc)%Nsplit
        ix=lf_src(isrc)%ix(iix)
        ffac = fac*1.e6*lf_src(isrc)%mw(iix) !(units ok?)
        istart = lf_src(isrc)%start
        iend = lf_src(isrc)%end
        if(isrc==isrc_SO4 .or. isrc==isrc_SO2 .or. lf_src(isrc)%species=="sox")ffac = ffac*32.0/64.0 !SO2->S
        if(isrc==isrc_NH3 .or. isrc==isrc_NH4 .or. lf_src(isrc)%species=="nh3")ffac = ffac* 14.0/17.0!NH3->N
        if(isrc==isrc_NO .or. isrc==isrc_NO2 .or. lf_src(isrc)%species=="nox")ffac = ffac*14.0/46.0 !NO2->N

!        if( ix==SO4_ix - NSPEC_SHL ) then
!           !take directly local fractions from SO4 instead of sox
!           istart = lf_src(isrc_SO4)%start
!           iend= lf_src(isrc_SO4)%end
!        endif
!        if( ix==SO2_ix - NSPEC_SHL ) then
!           !take directly local fractions from SO2 instead of sox
!           istart = lf_src(isrc_SO2)%start
!           iend= lf_src(isrc_SO2)%end
!        endif

!        if( ix==NH4_f_ix - NSPEC_SHL ) then
!           istart = lf_src(isrc_NH4)%start
!           iend= lf_src(isrc_NH4)%end
!        endif
!        if( ix==NH3_ix - NSPEC_SHL ) then
!           istart = lf_src(isrc_NH3)%start
!           iend= lf_src(isrc_NH3)%end
!        endif

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
        endif
     enddo
     idep0 = idep
  enddo
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")
  if(DEBUGall .and. me==0)write(*,*)'end drydep'
end subroutine lf_drydep

subroutine  lf_wetdep(iadv, i,j,k_in,loss, fac)
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
     if (lf_src(isrc)%species == 'FULLCHEM') cycle
     do iix=1,lf_src(isrc)%Nsplit
        ix=lf_src(isrc)%ix(iix)
        if(ix /= iadv) cycle
        ffac = fac*1.e6*lf_src(isrc)%mw(iix)
        istart = lf_src(isrc)%start
        iend = lf_src(isrc)%end
        if(isrc==isrc_SO2 .or. lf_src(isrc)%species=="sox") ffac = ffac*32.0/64.0 !SO2->S
        if(isrc==isrc_SO4) ffac = ffac*32.0/96.0 !SO4->S
        if(isrc==isrc_NH3 .or. lf_src(isrc)%species=="nh3") ffac = ffac* 14.0/17.0 !NH3->N
        if(isrc==isrc_NH4) ffac = ffac* 14.0/18.0 !NH4->N
        if(isrc==isrc_NO .or. isrc==isrc_NO2 .or. lf_src(isrc)%species=="nox") ffac = ffac*14.0/46.0 !NO2->N

!        if( ix==SO4_ix - NSPEC_SHL ) then
!           !take directly local fractions from SO4 instead of sox
!           istart = lf_src(isrc_SO4)%start
!           iend= lf_src(isrc_SO4)%end
!        endif
!        if( ix==SO2_ix - NSPEC_SHL ) then
!           !take directly local fractions from SO2 instead of sox
!           istart = lf_src(isrc_SO2)%start
!           iend= lf_src(isrc_SO2)%end
!        endif

!        if( ix==NH4_f_ix - NSPEC_SHL ) then
!           istart = lf_src(isrc_NH4)%start
!           iend= lf_src(isrc_NH4)%end
!        endif
!        if( ix==NH3_ix - NSPEC_SHL ) then
!           istart = lf_src(isrc_NH3)%start
!           iend= lf_src(isrc_NH3)%end
!        endif

        idep=idep0
        do n = istart, iend
           idep=idep+1
           loc_frac_wetdep(i,j,idep) = loc_frac_wetdep(i,j,idep) + lf(n,i,j,k)*loss*ffac
        enddo
        if(lf_src(isrc)%species=="nox" .and. iix==lf_src(isrc)%Nsplit)then
           !we add also depositions of NO3 and HNO3
           ix=NO3_ix - NSPEC_SHL
           idep=idep0
           do n = istart, iend
              idep=idep+1
              loc_frac_wetdep(i,j,idep) = loc_frac_wetdep(i,j,idep) + lf(n,i,j,k)*loss*ffac
           enddo
           ix=HNO3_ix - NSPEC_SHL
           idep=idep0
           do n = istart, iend
              idep=idep+1
              loc_frac_wetdep(i,j,idep) = loc_frac_wetdep(i,j,idep) + lf(n,i,j,k)*loss*ffac
           enddo
        endif
     enddo
     idep0 = idep0 + lf_src(isrc)%end-lf_src(isrc)%start+1
  enddo
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")

end subroutine lf_wetdep

subroutine save_lf_emis(s,i,j,iem,isec,iland)
  ! save emission data
  ! save without splits and k
  ! save all (anthropogenic) emissions, even if not used by LF

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
  emis_lf_cntry(i,j,ic,isec,iem)=emis_lf_cntry(i,j,ic,isec,iem) + s

  call Add_2timing(NTIMING-4,tim_after,tim_before,"lf: emissions")

  if(DEBUGall .and. me==0)write(*,*)'end lf emis'
end subroutine save_lf_emis

subroutine lf_rcemis(i,j,k,eps)
  !makes the emission differences to be used for derivatives wrt emissions
  !We need to make rcemis_lf for each source (emitted species, country, sector) with emissions at this i,j

  ! The number of emitting countries at i,j is nic(i,j) (small)
  ! The emission without splits and k is emis_lf_cntry(i,j,ic,isec,iem) (ic runs over nic(i,j), i.e. corresponds to different countries in different gridcells)
  ! for "fullchem": same lf_rcemis for all advected_species (but different for nox and voc sources), but lf countries dependent
  ! for "fullchem" lf_rcemis is splitted, otherwise summed with mw weights
  ! for "relative": country independent
  integer,intent(in) :: i,j,k
  real, intent(in) :: eps
  integer :: n, n0, iem, iqrc, itot, ic, iic, ig, is, iland, isec, split_idx
  integer :: emish_idx, nemis, found,iiix,isrc, nsectors_loop, iisec
  real :: ehlpcom0 = GRAV* 0.001*AVOG !0.001 = kg_to_g / m3_to_cm3
  real :: emiss,fac ! multiply emissions by

  call Code_timer(tim_before)
  if(DEBUGall .and. me==0)write(*,*)'start lf rcemis'

  rcemis_lf = 0.0 !default: no difference TODO: avoid setting entire array to zero!
  !1) For now, we want to take derivative only from sector emissions, i.e. gridrcemis, and not fire, lightning, natural etc.
  if(k < KEMISTOP) return
  nemis = 0 !counts the number of source with non zero emissions in this particular gridcell (i,j,k)
  N_lf_derivemis = 0 !final max value of nemis
  fac = 1.0
  if (lf_fullchem) fac = eps
  do iem = 1, NEMIS_File
    if (iem2Nipoll(iem) <= 0) cycle
    !for fullchem, we only treat nox and voc emissions
    if(nfullchem >0 .and. EMIS_FILE(iem)/='nox' .and. EMIS_FILE(iem)/='voc') cycle
    !we calculate the delta in emission for each country that has emissions
    if (.not.lf_fullchem) then
       do isrc=1,Nsources
          if (lf_src(isrc)%iem /= iem) cycle
          if(lf_src(isrc)%type=="relative")then
              !no lf_country involved, all countries are treated together
              n0 = 0
              !do not loop over countries, only sectors
              if(lf_src(isrc)%sector>0) then
                 nsectors_loop = 1
              else
                nsectors_loop = NSECTORS
              end if
              found = 0 !flag to show if nemis already increased
              do iisec=1,nsectors_loop
                 if(lf_src(isrc)%sector>0) then
                    isec = lf_src(isrc)%sector
                 else
                    isec = iisec
                 end if

                 do ic = 1,nic(i,j)
                    iland = ic2iland(i,j,ic)
                    if(emis_lf_cntry(i,j,ic,isec,iem)>1.E-20)then
                       if(found == 0) nemis = nemis + 1 !should be increased at most once for each isrc
                       found = 1
                       emish_idx = SECTORS(isec)%height
                       split_idx = SECTORS(isec)%split
                       do n = 1, lf_src(isrc)%Nsplit
                          itot=lf_src(isrc)%ix(n)
                          iqrc = itot2iqrc(itot+NSPEC_SHL)
                          if (iqrc<0) cycle
                          rcemis_lf(nemis,1) = rcemis_lf(nemis,1) +&
                               fac* emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)*ehlpcom0&
                               *emisfrac(iqrc,split_idx,iland)*emis_masscorr(iqrc)&
                               *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))  * lf_src(isrc)%mw(n)
                       end do
                    end if
                 end do
              end do

              if (found == 1) then
                 !emissions found here
                 !emis2icis(nemis) = n0 !not country/sector specific
                 emis2isrc(nemis) = isrc
              endif

           else if(lf_src(isrc)%type=="country")then
             do ic=1,nic(i,j)
                 iland = ic2iland(i,j,ic)
                 !only pick one country at a time
                 do iic=1,Ncountry_lf + Ncountry_group_lf !this loop could be avoided if necessary, by defining iland2iic?
                   if(iic<=Ncountry_lf)then
                     if (iic <= Ncountry_mask_lf) then
                       if (country_mask_val(iic)>-99999) then
                         if(EmisMaskIntVal(i,j) /= country_mask_val(iic)) cycle !remove all parts which are not covered by mask
                       end if
                     else if(country_ix_list(iic)==IC_TMT.and.(iland==IC_TM.or.iland==IC_TME))then
                     else if(country_ix_list(iic)==IC_AST.and.(iland==IC_ASM.or.iland==IC_ASE.or.iland==IC_ARE.or.iland==IC_ARL.or.iland==IC_CAS))then
                     else if(country_ix_list(iic)==IC_UZT.and.(iland==IC_UZ.or.iland==IC_UZE))then
                     else if(country_ix_list(iic)==IC_KZT.and.(iland==IC_KZ.or.iland==IC_KZE))then
                     else if(country_ix_list(iic)==IC_RUE.and.(iland==IC_RU.or.iland==IC_RFE.or.iland==IC_RUX))then
                     else if(country_ix_list(iic)/=iland)then
                       cycle !only pick out one country
                     endif
                   else
                     !defined as a group of countries
                     found = 0
                     do ig = 1, MAX_lf_country_group_size
                       if(lf_country%group(iic-Ncountry_lf)%list(ig) == 'NOTSET') exit
                       if(lf_country%group(iic-Ncountry_lf)%list(ig) == Country(iland)%code) then
                         found = 1
                         exit
                       end if
                     end do

                     if (found == 0) cycle
                   end if
                   do is=1,Ncountrysectors_lf
                       isec=lf_country%sector_list(is)

                       if(lf_country%sector_list(is)>0) then
                          nsectors_loop = 1
                       else
                          nsectors_loop = NSECTORS
                       end if
                       found = 0 !flag to show if nemis already increased
                       do iisec=1,nsectors_loop
                          if(lf_country%sector_list(is)>0) then
                             isec = lf_country%sector_list(is)
                          else
                             isec = iisec
                          end if

                          if(emis_lf_cntry(i,j,ic,isec,iem)>1.E-20)then
                             if(found == 0) then
                                !first see if this emis2icis and emis2isrc already exists (for groups and masks type "countries"):
                                do n=1,nemis
                                   if(emis2icis(n)==is-1 + (iic-1)*Ncountrysectors_lf .and. emis2isrc(n) == iem)then
                                      ! add to this instead
                                      nemis = n
                                      found = 1
                                      exit
                                   end if
                                end do
                             end if
                             if(found == 0)  then
                                nemis = N_lf_derivemis + 1 !should be increased at most once if sector=0
                                N_lf_derivemis = max(nemis,N_lf_derivemis)
                             end if

                             found = 1
                             emish_idx = SECTORS(isec)%height
                             split_idx = SECTORS(isec)%split
                             do n = 1, lf_src(isrc)%Nsplit
                                itot=lf_src(isrc)%ix(n)
                                iqrc = itot2iqrc(itot+NSPEC_SHL)
                                if (iqrc<0) cycle
                                rcemis_lf(nemis,1) = rcemis_lf(nemis,1) +&
                               fac* emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)*ehlpcom0&
                               *emisfrac(iqrc,split_idx,iland)*emis_masscorr(iqrc)&
                               *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))  * lf_src(isrc)%mw(n)
                             end do
                          end if
                       end do
                       if (found == 1) then
                         !emissions found here
                         emis2icis(nemis) = is-1 + (iic-1)*Ncountrysectors_lf
                         emis2isrc(nemis) = isrc
                       endif
                    end do
                  end do
              end do

           end if
        end do
      else
        !TODO : merge with case not fullchem? difference only lf_src(isrc)%mw(n) and split summation?
         do ic=1,nic(i,j)

          !iland is the country index f the emission treated
          iland = ic2iland(i,j,ic)
          !iic is the lf country index of the source
          do iic=1,Ncountry_lf + Ncountry_group_lf !this loop could be avoided if necessary, by defining iland2iic?
            if(iic<=Ncountry_lf)then
              if (iic <= Ncountry_mask_lf) then
                if (country_mask_val(iic)>-99999) then
                  if(EmisMaskIntVal(i,j) /= country_mask_val(iic)) cycle !remove all parts which are not covered by mask
                end if
              else if(country_ix_list(iic)==IC_TMT.and.(iland==IC_TM.or.iland==IC_TME))then
              else if(country_ix_list(iic)==IC_AST.and.(iland==IC_ASM.or.iland==IC_ASE.or.iland==IC_ARE.or.iland==IC_ARL.or.iland==IC_CAS))then
              else if(country_ix_list(iic)==IC_UZT.and.(iland==IC_UZ.or.iland==IC_UZE))then
              else if(country_ix_list(iic)==IC_KZT.and.(iland==IC_KZ.or.iland==IC_KZE))then
              else if(country_ix_list(iic)==IC_RUE.and.(iland==IC_RU.or.iland==IC_RFE.or.iland==IC_RUX))then
              else if(country_ix_list(iic)/=iland)then
                cycle !only pick out one country
              endif
            else
              !defined as a group of countries
              found = 0
              do ig = 1, MAX_lf_country_group_size
                if(lf_country%group(iic-Ncountry_lf)%list(ig) == 'NOTSET') exit
                if(lf_country%group(iic-Ncountry_lf)%list(ig) == Country(iland)%code) then
                  found = 1
                  exit
                end if
              end do
              if (found == 0) cycle
            end if
            do is=1,Ncountrysectors_lf
              isec=lf_country%sector_list(is)
              if(isec>0) then
                nsectors_loop = 1
              else
                nsectors_loop = NSECTORS
              end if
              found = 0 !flag to show if nemis already increased
              do iisec=1,nsectors_loop
                if(lf_country%sector_list(is)>0) then
                  isec = lf_country%sector_list(is)
                else
                  isec = iisec
                end if
                if(emis_lf_cntry(i,j,ic,isec,iem)>1.E-20)then
                   if(found == 0) then
                      !first see if this emis2icis and emis2isrc already exists (for groups and masks type "countries"):
                      !because only one derivative for each source is included ("lf(n0,i,j,k) = ederiv(iemis,ix)" without +=)
                      do n=1,nemis
                         if(emis2icis(n)==is-1 + (iic-1)*Ncountrysectors_lf .and. emis2isrc(n) == iem)then
                            ! add to this instead
                            nemis = n
                            found = 1
                            exit
                         end if
                      end do
                   end if
                   if(found == 0) then
                      nemis = N_lf_derivemis + 1 !should be increased at most once for each is
                      N_lf_derivemis = max(nemis,N_lf_derivemis)
                   end if
                   found = 1
                  emish_idx = SECTORS(isec)%height
                  split_idx = SECTORS(isec)%split
                  do n = 1, emis_nsplit(iem)
                    iqrc = sum(emis_nsplit(1:iem-1)) + n
                    itot = iqrc2itot(iqrc)

                    if(itot<=NSPEC_fullchem_inc_lf .or. .not.lf_fullchem)then !temporary: avoid SHIPNOx
                      rcemis_lf(nemis,itot) = rcemis_lf(nemis,itot) +&
                          fac* emis_lf_cntry(i,j,ic,isec,iem)*emis_kprofile(KMAX_BND-k,emish_idx)*ehlpcom0&
                          *emisfrac(iqrc,split_idx,iland)*emis_masscorr(iqrc)&
                          *roa(i,j,k,1)/(dA(k)+dB(k)*ps(i,j,1))

                   end if
                  end do
                end if
              end do
              if (found == 1) then
                !emissions found here
                emis2icis(nemis) = is-1 + (iic-1)*Ncountrysectors_lf
                emis2isrc(nemis) = iem
              endif
            end do
          end do
        end do
      end if
    end do

    N_lf_derivemis = nemis !final number of emis derivatives computed for this i,j,k

    if(nemis>N_lf_derivemisMAX)then
      write(*,*)me,i,j,k,nemis,N_lf_derivemisMAX
      call StopAll("too many sectors*country in one gridcell. Increase N_lf_derivemisMAX")
    end if

    call Add_2timing(NTIMING-4,tim_after,tim_before,"lf: emissions")
  if(DEBUGall .and. me==0)write(*,*)'end lf rcemis'
  end subroutine lf_rcemis

end module LocalFractions_mod
