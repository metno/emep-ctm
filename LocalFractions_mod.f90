! <LocalFractions_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.36>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2020 met.no
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
use ChemDims_mod,      only: NSPEC_ADV, NSPEC_SHL,NEMIS_File
use ChemSpecs_mod,     only: species_adv,species
use Config_module,     only: KMAX_MID, KMAX_BND,USES, uEMEP, lf_src, IOU_HOUR&
                             , IOU_HOUR_INST,IOU_INST,IOU_YEAR,IOU_MON,IOU_DAY&
                             ,IOU_HOUR,IOU_HOUR_INST, IOU_MAX_MAX &
                             ,MasterProc,dt_advec, RUNDOMAIN, runlabel1 &
                             ,HOURLYFILE_ending, MAXSRC, Max_Country_list, Max_Country_sectors &
                             ,lf_country_group, Max_Country_groups&
                             ,lf_country_sector_list,lf_country_list
use Country_mod,       only: MAXNLAND,NLAND,Country&
                             ,IC_TMT,IC_TM,IC_TME,IC_ASM,IC_ASE,IC_ARE,IC_ARL,IC_CAS,IC_UZT,IC_UZ&
                             ,IC_UZE,IC_KZT,IC_KZ,IC_KZE,IC_RU,IC_RFE,IC_RUX,IC_RUE,IC_AST

use DefPhotolysis_mod, only: IDNO2
use EmisDef_mod,       only: lf, emis_lf, lf_emis_tot, emis_lf_cntry, loc_frac_src_1d,&
                            lf_src_acc,lf_src_tot,lf_src_full,loc_tot_full, NSECTORS,EMIS_FILE, &
                            loc_frac_drydep, loc_frac_wetdep, &
                            nlandcode,landcode,sec2hfac_map, sec2split_map,&
                            ISNAP_DOM,secemis, roaddust_emis_pot,KEMISTOP,&
                            NEmis_sources, Emis_source_2D, Emis_source
use EmisGet_mod,       only: nrcemis, iqrc2itot, emis_nsplit,nemis_kprofile, emis_kprofile,&
                             make_iland_for_time,itot2iqrc,iqrc2iem, emisfrac
use GridValues_mod,    only: dA,dB,xm2, dhs1i, glat, glon, projection, extendarea_N,i_fdom,j_fdom
use MetFields_mod,     only: ps,roa,EtaKz
use MPI_Groups_mod
use NetCDF_mod,        only: Real4,Out_netCDF,LF_ncFileID_iou
use OwnDataTypes_mod,  only: Deriv, Npoll_lf_max, Nsector_lf_max, MAX_lf_country_group_size, &
                             TXTLEN_NAME, TXTLEN_FILE
use Par_mod,           only: me,LIMAX,LJMAX,MAXLIMAX,MAXLJMAX,gi0,gj0,li0,li1,lj0,lj1,GIMAX,GJMAX
use PhysicalConstants_mod, only : GRAV, ATWAIR 
use SmallUtils_mod,    only: find_index
!use Chemsolver_mod,    only: Dchem
use TimeDate_mod,      only: date, current_date,day_of_week
use TimeDate_ExtraUtil_mod,only: date2string
use My_Timing_mod,     only: Add_2timing, Code_timer, NTIMING
use ZchemData_mod,only: rct, rcphot, xn_2d, rcemis

!(dx,dy,i,j) shows contribution of pollutants from (i+dx,j+dy) to (i,j)

implicit none
!external advection_mod_mp_vertdiffn_k

private

integer ::IC_AST_EXTRA = 324567,IC_RUT_EXTRA = 324568 !sum of countries which are not defined as countries
integer ::IC_BIC_EXTRA = 324569

public  :: lf_init
public  :: lf_out
public  :: lf_av
public  :: lf_adv_x
public  :: lf_adv_y
public  :: lf_adv_k
public  :: lf_diff
public  :: lf_chem
public  :: lf_aero_pre, lf_aero_pos 
public  :: lf_drydep, lf_wetdep
public  :: lf_emis
public  :: add_lf_emis

real, allocatable, save ::loc_poll_to(:,:,:,:,:)

logical, public, save :: COMPUTE_LOCAL_TRANSPORT=.false.
integer , private, save :: lf_Nvertout = 1!number of vertical levels to save in output
integer, public, save :: NTIMING_lf=9
real, private :: tim_after,tim_before
integer, public, save :: Ndiv_coarse=1, Ndiv_rel=1, Ndiv2_coarse=1
integer, public, save :: Nsources=0
integer, public, save :: lf_Nvert=0


integer, public, save :: LF_SRC_TOTSIZE
integer,  public, save :: iotyp2ix(IOU_MAX_MAX)
integer,  public, save :: av_fac(IOU_MAX_MAX)
integer,  public, save :: Niou_ix = 0 ! number of time periods to consider (hourly, monthly, full ...)
integer,  public, save :: Npoll = 0 !Number of different pollutants to consider
integer, private, parameter :: MAXIPOLL=16
integer,  public, save :: iem2ipoll(NEMIS_File,MAXIPOLL) !internal indices of pollutants for that emis file
integer,  public, save :: ipoll2iqrc(MAXIPOLL)=-1 !-1 for primary pollutant
integer, private, save :: iem2Nipoll(NEMIS_File) !number of pollutants for that emis file
logical :: old_format=.false. !temporary, use old format for input and output
integer, private, save :: isrc_O3=-1, isrc_NO=-1, isrc_NO2=-1, isrc_VOC=-1
integer, private, save :: isrc_SO2=-1, isrc_SO4=-1, isrc_NH4=-1, isrc_NH3=-1
integer, private, save :: ix_O3=-1, ix_NO2=-1, ix_NO=-1, ix_CH3CO3=-1, ix_HO2=-1
integer, private, save :: ix_SO4=-1, ix_SO2=-1, ix_H2O2=-1, ix_OH=-1
integer, private, save :: ix_NH4=-1, ix_NH3=-1
integer, private, save :: ix_NO3=-1, ix_HNO3=-1
integer, private, save :: isrc_EC_f_ffuel_new=-1, isrc_EC_f_ffuel_age=-1, isrc_EC_f_wood_new=-1, isrc_EC_f_wood_age=-1
integer, private, save :: isrc_EC_f_new=-1, isrc_EC_f_age=-1
integer, private, save :: ix_EC_f_ffuel_new=-1, ix_EC_f_ffuel_age=-1
integer, private, save :: ix_EC_f_new=-1, ix_EC_f_age=-1
integer, private, save :: ix_EC_f_wood_new=-1, ix_EC_f_wood_age=-1
real, allocatable, private, save :: lf_NH4(:), lf_NH3(:)
integer, private, save :: country_ix_list(Max_Country_list)
integer, private, save :: Ncountry_lf=0
integer, private, save :: Ncountry_group_lf=0
integer, private, save :: Ncountrysectors_lf=0
integer, private, save :: Ndrydep_lf=0
integer, private, save :: Nwetdep_lf=0

contains

  subroutine lf_init
    integer :: i, ic, ix, itot, iqrc, iem, iemis, isec, ipoll, ixnh3, ixnh4, size, IOU_ix, isrc
    integer :: found
    character(len=TXTLEN_NAME) :: iem2names(NEMIS_File,MAXIPOLL) !name of that pollutant

  call Code_timer(tim_before)
  ix=0
  if(USES%uEMEP)then
     !Temporary: we keep compatibilty with lf input
     old_format=.true.
     lf_src(:)%dist = uEMEP%dist !Temporary
     lf_src(:)%Nvert = uEMEP%Nvert !Temporary
     do i=1,4
        lf_src(:)%DOMAIN(i) = uEMEP%DOMAIN(i) !Temporary
        if(lf_src(1)%DOMAIN(i)<0)lf_src(:)%DOMAIN(i) = RUNDOMAIN(i)
        if(me==0)write(*,*)i,' DOMAIN ',lf_src(1)%DOMAIN(i)
     enddo
     lf_src(:)%YEAR=uEMEP%YEAR
     lf_src(:)%MONTH=uEMEP%MONTH
     lf_src(:)%MONTH_ENDING=uEMEP%MONTH_ENDING
     lf_src(:)%DAY=uEMEP%DAY
     lf_src(:)%HOUR=uEMEP%HOUR
     lf_src(:)%HOUR_INST=uEMEP%HOUR_INST
     do isrc=1,Npoll_lf_max
        if(uEMEP%poll(isrc)%emis=='none')then
           call CheckStop(isrc==1,"init_uEMEP: no pollutant specified")
           exit
        else
           do isec=1,Nsector_lf_max
              if(uEMEP%poll(isrc)%sector(isec)<0)then
                 call CheckStop(isec==0,"init_uEMEP: nosector specified for "//uEMEP%poll(isrc)%emis)
                 exit
              else
                 ix=ix+1
                 lf_src(ix)%species = uEMEP%poll(isrc)%emis
                 lf_src(ix)%sector = uEMEP%poll(isrc)%sector(isec)
              endif
           enddo
        endif
     enddo
  else
     !separate value do not work properly yet
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
  Nsources = 0
  do i = 1, MAXSRC
     if(lf_src(i)%species /= 'NONE') Nsources = Nsources + 1
  enddo

  ipoll=0
  iem2ipoll = -1
  iem2Nipoll = 0
  do isrc = 1, Nsources
     !for now only one Ndiv possible for all sources
     if(lf_src(isrc)%type == 'relative')then
        lf_src(isrc)%Npos =  (2*lf_src(isrc)%dist+1)*(2*lf_src(isrc)%dist+1)
        Ndiv_rel = max(Ndiv_rel,2*lf_src(isrc)%dist+1)
     endif
     if(lf_src(isrc)%type == 'coarse')then
        Ndiv_coarse = max(Ndiv_coarse,2*lf_src(isrc)%dist+1)
        lf_src(isrc)%Npos = (2*lf_src(isrc)%dist+1)*(2*lf_src(isrc)%dist+1)
        Ndiv2_coarse = max(Ndiv2_coarse,Ndiv_coarse*Ndiv_coarse)
     endif
     if(lf_src(isrc)%type == 'country')then
        if(lf_country_list(1)/= 'NOTSET' .or. lf_country_group(1)%name/= 'NOTSET')then
           !list of countries/sectors instead of single country
           Ncountry_lf=0
           do i = 1, Max_Country_list
              if(lf_country_list(i) == 'NOTSET') exit
              Ncountry_lf=Ncountry_lf+1
              ix = find_index(trim(lf_country_list(i)) ,Country(:)%code, first_only=.true.)
              if(ix<0 .and. lf_country_list(i) =='AST')then
                 ix = IC_AST_EXTRA 
              else if(ix<0 .and. lf_country_list(i) =='RUT')then
                 ix = IC_RUT_EXTRA 
              else if(ix<0 .and. lf_country_list(i) =='BIC')then
                 ix = IC_BIC_EXTRA 
              endif 
              call CheckStop(ix<0,'country '//trim(lf_country_list(i))//' not defined. ')        
              country_ix_list(i) = ix
              if(MasterProc)write(*,*)'include '//trim(lf_src(isrc)%species)//' from ',trim(lf_country_list(i))
           enddo
           Ncountry_group_lf=0
           do i = 1, Max_Country_groups
              if(lf_country_group(i)%name == 'NOTSET') exit
              Ncountry_group_lf = Ncountry_group_lf+1
              do ic = 1, MAX_lf_country_group_size
                 if(lf_country_group(i)%list(ic) == 'NOTSET') exit
                 
                 ix = find_index(trim(lf_country_group(i)%list(ic)) ,Country(:)%code, first_only=.true.)
                 if(ix<0 .and. lf_country_list(i) =='AST')then
                    ix = IC_AST_EXTRA 
                 else if(ix<0 .and. lf_country_list(i) =='RUT')then
                    ix = IC_RUT_EXTRA 
                 else if(ix<0 .and. lf_country_list(i) =='BIC')then
                    ix = IC_BIC_EXTRA 
                 endif
                 call CheckStop(ix<0,'country '//trim(lf_country_group(i)%list(ic))//' not defined. ')        
                 lf_country_group(i)%ix(ic) = ix
                 if(MasterProc)write(*,*)'include '//trim(lf_src(isrc)%species)//' from '//&
                      trim(lf_country_group(i)%list(ic))//' as '//trim(lf_country_group(i)%name)
              enddo
           enddo
           Ncountrysectors_lf=0
           do i = 1, Max_Country_sectors
              if(lf_country_sector_list(i) < 0) exit
              Ncountrysectors_lf=Ncountrysectors_lf+1
              if(MasterProc)write(*,*)'country sector ',lf_country_sector_list(i)
          enddo
          lf_src(isrc)%Npos = (Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf
          if(MasterProc)write(*,*)lf_src(isrc)%Npos,' countries x sectors'
        else
           lf_src(isrc)%Npos = 1
           ix = find_index(trim(lf_src(isrc)%country_ISO) ,Country(:)%code, first_only=.true.)
           call CheckStop(ix<0,'country '//trim(lf_src(isrc)%country_ISO)//' not defined. ')        
           lf_src(isrc)%country_ix = ix
           lf_src(isrc)%Npos = 1
        endif
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
     
     do i=1,4
        lf_src(isrc)%DOMAIN(i) = max(RUNDOMAIN(i),lf_src(isrc)%DOMAIN(i))
     enddo
     
     iem=find_index(lf_src(isrc)%species ,EMIS_FILE(1:NEMIS_FILE))
     if(iem<1)then
        !defined as single species (NO, NO2, O3..)
        lf_src(isrc)%Nsplit = 1
        ix=find_index(lf_src(isrc)%species ,species(:)%name)
        if(ix<0)then
           ix=find_index(lf_src(isrc)%species ,species(:)%name, any_case=.true.) !NB: index among all species also short lived
           if(me==0 .and. ix>0)then
              write(*,*)'WARNING: '//trim(lf_src(isrc)%species)//' not found, replacing with '//trim(species(ix)%name)
              lf_src(isrc)%species=trim(species(ix)%name)
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
        if(trim(species(ix)%name)=='O3')isrc_O3=isrc
        if(trim(species(ix)%name)=='NO')isrc_NO=isrc
        if(trim(species(ix)%name)=='NO2')isrc_NO2=isrc
        if(trim(species(ix)%name)=='O3')ix_O3=lf_src(isrc)%ix(1)!shortcut
        if(trim(species(ix)%name)=='NO')ix_NO=lf_src(isrc)%ix(1)!shortcut
        if(trim(species(ix)%name)=='NO2')ix_NO2=lf_src(isrc)%ix(1)!shortcut
        if(trim(species(ix)%name)=='SO4')isrc_SO4=isrc
        if(trim(species(ix)%name)=='SO2')isrc_SO2=isrc
        if(trim(species(ix)%name)=='SO4')ix_SO4=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='SO2')ix_SO2=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='NH4_f')isrc_NH4=isrc
        if(trim(species(ix)%name)=='NH3')isrc_NH3=isrc
        if(trim(species(ix)%name)=='NH4_f')ix_NH4=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='NH3')ix_NH3=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='EC_f_new')isrc_EC_f_ffuel_new=isrc
        if(trim(species(ix)%name)=='EC_f_new')ix_EC_f_ffuel_new=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='EC_f_age')isrc_EC_f_ffuel_age=isrc
        if(trim(species(ix)%name)=='EC_f_age')ix_EC_f_ffuel_age=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='EC_f_ffuel_new')isrc_EC_f_ffuel_new=isrc
        if(trim(species(ix)%name)=='EC_f_ffuel_new')ix_EC_f_ffuel_new=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='EC_f_ffuel_age')isrc_EC_f_ffuel_age=isrc
        if(trim(species(ix)%name)=='EC_f_ffuel_age')ix_EC_f_ffuel_age=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='EC_f_wood_new')isrc_EC_f_wood_new=isrc
        if(trim(species(ix)%name)=='EC_f_wood_new')ix_EC_f_wood_new=lf_src(isrc)%ix(1)
        if(trim(species(ix)%name)=='EC_f_wood_age')isrc_EC_f_wood_age=isrc
        if(trim(species(ix)%name)=='EC_f_wood_age')ix_EC_f_wood_age=lf_src(isrc)%ix(1)
     else
        !species defines as primary emitted 
        lf_src(isrc)%iem = iem
        lf_src(isrc)%Nsplit=emis_nsplit(iem)
        do i=1,lf_src(isrc)%Nsplit
           iqrc=sum(emis_nsplit(1:iem-1)) + i
           itot=iqrc2itot(iqrc)
           ix=itot-NSPEC_SHL
           lf_src(isrc)%ix(i)=ix
           lf_src(isrc)%mw(i)=species_adv(ix)%molwt
!           if(lf_src(isrc)%species=="pm25")then
!              ix=find_index('EC_f_ffuel_new', species(:)%name)
!              call CheckStop( ix<1, "Local Fractions did not find  ")
!              ix_EC_f_ffuel_new = ix - NSPEC_SHL!shortcut
!              ix=find_index('EC_f_ffuel_age', species(:)%name)
!              call CheckStop( ix<1, "Local Fractions did not find  ")
!              ix_EC_f_ffuel_age = ix - NSPEC_SHL!shortcut
!              ix=find_index('isrc_EC_f_wood_new', species(:)%name)
!              call CheckStop( ix<1, "Local Fractions did not find  ")
!              ix_isrc_EC_f_wood_new = ix - NSPEC_SHL!shortcut
!              ix=find_index('isrc_EC_f_wood_age', species(:)%name)
!              call CheckStop( ix<1, "Local Fractions did not find  ")
!              ix_isrc_EC_f_wood_age = ix - NSPEC_SHL!shortcut
!           endif
           if(lf_src(isrc)%species=="voc")then
              isrc_VOC = isrc
              ix=find_index('CH3CO3', species(:)%name)
              call CheckStop( ix<1, "Local Fractions did not find CH3CO3 ")
              ix_CH3CO3 = ix - NSPEC_SHL!shortcut
              ix=find_index('HO2', species(:)%name)
              call CheckStop( ix<1, "Local Fractions did not find HO2 ")
              ix_HO2= ix!shortcut !NB: index as short lived
           endif
           if(lf_src(isrc)%species=="nox")then
              ix=find_index("NO2",species_adv(:)%name)
              call CheckStop(ix<0,'Index for NO2 not found')
              lf_src(isrc)%mw(i)=species_adv(ix)%molwt
           endif
           if(lf_src(isrc)%species=="sox")then
              ix=find_index("SO2",species_adv(:)%name)
              call CheckStop(ix<0,'Index for SO2 not found')
              lf_src(isrc)%mw(i)=species_adv(ix)%molwt
           endif
           if(lf_src(isrc)%species=="nox" .and. (lf_src(isrc)%DryDep .or. lf_src(isrc)%WetDep))then
              ix=find_index("NO3",species_adv(:)%name)
              call CheckStop(ix<0,'Index for NO3 not found')
              ix_NO3=ix
              ix=find_index("HNO3",species_adv(:)%name)
              call CheckStop(ix<0,'Index for HNO3 not found')
              ix_HNO3=ix
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
           endif
        end do
        
     endif
     if(ix_SO2>0)then
        !need some more indices
        ix_O3=find_index('O3' ,species(:)%name, any_case=.true.) - NSPEC_SHL !NB: index among advected species
        ix_OH=find_index('OH' ,species(:)%name, any_case=.true.) !NB: index among all species
        ix_H2O2=find_index('H2O2' ,species(:)%name, any_case=.true.) !NB: index among all species
     endif
     if(iem>0)then
        !emitted species
        found=0
        do i=1,iem2Nipoll(iem)
           if(iem2names(iem,i)==lf_src(isrc)%species)then
              found=1
              lf_src(isrc)%poll = iem2ipoll(iem,i)
           endif
        enddo
        if(found==0)then
           !add a new pollutant for that emis file
           iem2Nipoll(iem)=iem2Nipoll(iem)+1
           ipoll = ipoll + 1
           call CheckStop(ipoll>MAXIPOLL,"Error: increase MAXIPOLL in LocalFractions_mod")
           iem2names(iem, iem2Nipoll(iem)) = trim(lf_src(isrc)%species)    
           iem2ipoll(iem, iem2Nipoll(iem)) = ipoll     
           lf_src(isrc)%poll = ipoll
           Npoll = ipoll
        endif
        if(lf_src(isrc)%iqrc>0)ipoll2iqrc(lf_src(isrc)%poll)=lf_src(isrc)%iqrc!single species
      else
        !single species not emitted (like O3)
        ipoll = ipoll + 1
        lf_src(isrc)%poll = ipoll
        Npoll = ipoll
        call CheckStop(ipoll>MAXIPOLL,"Error: increase MAXIPOLL in LocalFractions_mod")
     endif
        
     if(MasterProc)then
        if(lf_src(isrc)%iem>0)then
           write(*,*)'lf pollutant : ',lf_src(isrc)%species,' ref index ',lf_src(isrc)%poll,' emitted as ',EMIS_FILE(lf_src(isrc)%iem)
        else
           write(*,*)'lf pollutant : ',lf_src(isrc)%species,' ref index ',lf_src(isrc)%poll,' not treated as emitted species'
        endif
        write(*,*)'lf number of species in '//trim(lf_src(isrc)%species)//' group: ',lf_src(isrc)%Nsplit
        write(*,"(A,30(A,F6.2))")'including:',('; '//trim(species_adv(lf_src(isrc)%ix(i))%name)//', mw ',lf_src(isrc)%mw(i),i=1,lf_src(isrc)%Nsplit)
        write(*,"(A,30I4)")'sector:',lf_src(isrc)%sector
        !write(*,"(A,30I4)")'ix:',(lf_src(isrc)%ix(i),i=1,lf_src(isrc)%Nsplit)
     end if
  end do
  if(isrc_O3>0 .and. (isrc_NO2<0 .or. isrc_NO<0))then
     if(me==0)write(*,*)'WARNING: O3 tracking requires NO2 and NO'
     stop!may be relaxed in future
  endif
  if(isrc_SO2>0 .and. (isrc_SO4<0))then
     if(me==0)write(*,*)'WARNING: SO2 tracking requires SO4'
     stop!may be relaxed in future
  endif

  av_fac=0.0
  
  Ndrydep_lf=0
  Nwetdep_lf=0
  LF_SRC_TOTSIZE = 0
  do isrc = 1, Nsources
     if(lf_src(isrc)%drydep) Ndrydep_lf = Ndrydep_lf + lf_src(isrc)%Npos
     if(lf_src(isrc)%wetdep) Nwetdep_lf = Nwetdep_lf + lf_src(isrc)%Npos
     lf_src(isrc)%start = LF_SRC_TOTSIZE + 1
     lf_src(isrc)%end = LF_SRC_TOTSIZE + lf_src(isrc)%Npos
     LF_SRC_TOTSIZE = LF_SRC_TOTSIZE + lf_src(isrc)%Npos
     if(me==0)then
        write(*,*)isrc,' ',trim(lf_src(isrc)%species)," start ",lf_src(isrc)%start," end ",lf_src(isrc)%end,LF_SRC_TOTSIZE
     endif
  enddo

  allocate(lf(LF_SRC_TOTSIZE,LIMAX,LJMAX,KMAX_MID-lf_Nvert+1:KMAX_MID))
  lf=0.0

  isrc=1!for now all must be the same
  if(lf_src(isrc)%HOUR)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_HOUR)=Niou_ix
  endif
  if(lf_src(isrc)%HOUR_INST)then
     !Niou_ix = Niou_ix + 1 !should not be accumulated
     iotyp2ix(IOU_HOUR_inst) = -1; !should not be accumulated
  endif  
  if(lf_src(isrc)%DAY)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_DAY)=Niou_ix
  endif
  if(lf_src(isrc)%MONTH)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_MON)=Niou_ix
  endif
  if(lf_src(isrc)%YEAR)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_YEAR)=Niou_ix
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
  allocate(emis_lf(LIMAX,LJMAX,KMAX_MID-lf_Nvert+1:KMAX_MID,Nsources))
  emis_lf = 0.0
  allocate(lf_emis_tot(LIMAX,LJMAX,KMAX_MID-lf_Nvert+1:KMAX_MID,Npoll))
  lf_emis_tot = 0.0
  if(Ncountry_lf*Ncountrysectors_lf>0)then
     allocate(emis_lf_cntry(LIMAX,LJMAX,KMAX_MID-lf_Nvert+1:KMAX_MID,Ncountry_lf+Ncountry_group_lf,Ncountrysectors_lf,Nsources))
     emis_lf_cntry=0.0
  else
     allocate(emis_lf_cntry(1,1,1,1,1,1))!So can be set to zero etc. without compiler complaining
  endif
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

!  call Add_2timing(NTIMING-10,tim_after,tim_before,"lf: init") negligible

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
  real :: fracsum(LIMAX,LJMAX)
  logical :: pollwritten(Npoll_lf_max)
  integer :: ncFileID
  
  call Code_timer(tim_before)

  if(iotyp==IOU_HOUR_INST .and. lf_src(1)%HOUR_INST)then
     fileName = trim(runlabel1)//'_uEMEP_hourInst'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
     if(oldhourlyInstname/=fileName)then
        first_call(iotyp) = .true.
        oldhourlyInstname = fileName
     endif
  else if(iotyp==IOU_HOUR .and. lf_src(1)%HOUR)then
     fileName = trim(runlabel1)//'_uEMEP_hour'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
     if(oldhourlyname/=fileName)then
        first_call(iotyp) = .true.
        oldhourlyname = fileName
     endif
  else if(iotyp==IOU_DAY .and. lf_src(1)%DAY)then
     fileName=trim(runlabel1)//'_uEMEP_day.nc'
  else if(iotyp==IOU_MON .and. lf_src(1)%MONTH)then
     if(lf_src(1)%MONTH_ENDING /= "NOTSET")then
        fileName=trim(runlabel1)//'_uEMEP_month'//date2string(trim(lf_src(1)%MONTH_ENDING),current_date,-1.0)
        if(oldmonthlyname/=fileName)then
           first_call(iotyp) = .true.
           oldmonthlyname = fileName
        endif
     else
        fileName=trim(runlabel1)//'_uEMEP_month.nc'
     endif
  else if(iotyp==IOU_YEAR .and. lf_src(1)%YEAR)then
     fileName=trim(runlabel1)//'_uEMEP_full.nc'
  else
     return
  endif
  ncFileID=LF_ncFileID_iou(iotyp)
  
  ndim=5
  ndim_tot=3
  kmax=lf_Nvertout
  scale=1.0
  CDFtype=Real4
  !  dimSizes(1)=uEMEP%Nsec_poll
  !  dimNames(1)='sector'
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
  def1%class='uEMEP' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=1.0      !not used
  def1%name=trim(varName)
  def1%unit=''
  def2=def1
  def2%unit='ug/m3'
  def3=def1
  def3%unit='mg/m2'
  chunksizes=1
  !chunksizes(1)=dimSizes(1) !slower!!
  !chunksizes(2)=dimSizes(2) !slower!!
  chunksizes(3)=MAXLIMAX
  chunksizes(4)=MAXLJMAX
  chunksizes(5)=dimSizes(5)
  chunksizes_tot=1
  chunksizes_tot(1)=MAXLIMAX
  chunksizes_tot(2)=MAXLJMAX
  chunksizes_tot(3)=dimSizes_tot(3)

  allocate(tmp_out(max(Ndiv2_coarse,Ndiv_rel*Ndiv_rel),LIMAX,LJMAX)) !NB; assumes KMAX=1 TEMPORARY
  allocate(tmp_out_cntry(LIMAX,LJMAX,(Ncountry_lf+Ncountry_group_lf)*Ncountrysectors_lf))
 
  iou_ix = iotyp2ix(iotyp)

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
                       invtot=1.0/av_fac(iotyp) !could also output fractions?
                       do n=lf_src(isrc)%start, lf_src(isrc)%end
                          n1=n1+1
                          tmp_out_cntry(i,j,n1) = tmp_out_cntry(i,j,n1) + lf_src_acc(n,i,j,k,iou_ix)*invtot ! sum over all k
                          !if(tmp_out_cntry(i,j,n1)<1.e-18)tmp_out_cntry(i,j,n1)=0.0
                          if(isnan(tmp_out_cntry(i,j,n1)).or. tmp_out_cntry(i,j,n1)>1.e19)then
                             write(*,*)'tmp_out_cntry is nan ',tmp_out_cntry(i,j,n1),lf_src_acc(n,i,j,k,iou_ix),invtot,trim(lf_src(isrc)%species)
                             stop
                          endif
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
                 isec=lf_country_sector_list(j)
                 if(i<=Ncountry_lf)then
                    write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country_list(i))
                    if(isec==0) write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country_list(i))
                 else
                    !country group
                    write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_'//trim(lf_country_group(i-Ncountry_lf)%name)
                    if(isec==0) write(def2%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_'//trim(lf_country_group(i-Ncountry_lf)%name)                    
                 endif
                 if(me==0 .and. iter==2)write(*,*)'writing '//trim(def2%name)
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
              enddo
           enddo
        else
           if(old_format)then
              !for backward compatibility
              write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_local_fraction'
              if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_local_fraction'
           else
              write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_fraction_'//trim(lf_src(isrc)%type)
              if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_fraction_'//trim(lf_src(isrc)%type)
              write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_sec',isec,'_fraction_'//trim(lf_src(isrc)%type)
              if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%species)//'_fraction_'//trim(lf_src(isrc)%type)
           endif
           scale=1.0
           call Out_netCDF(iotyp,def1,ndim,kmax,tmp_out,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes,ncFileID_given=ncFileID)
           overwrite=.false.
           if(isrc==1)then
              def1%name=trim(lf_src(isrc)%species)//'_fracsum'
              call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=lf_src(isrc)%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot,ncFileID_given=ncFileID)  
           endif           
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

! CALL MPI_BARRIER(MPI_COMM_CALC, I)
 
!stop
end subroutine lf_out

subroutine lf_av(dt,End_of_Day)
  real, intent(in)    :: dt                   ! time-step used in integrations
  logical, intent(in) :: End_of_Day           ! e.g. 6am for EMEP sites
  real :: xtot
  integer ::i,j,k,n,dx,dy,ix,iix,ipoll,isec_poll1, iou_ix, isrc
  integer ::isec_poll
  logical :: pollwritten(Npoll_lf_max)
  
  call Code_timer(tim_before)
  if(.not. lf_src(1)%HOUR.and.&
     .not. lf_src(1)%DAY .and.&
     .not. lf_src(1)%MONTH .and.&
     .not. lf_src(1)%YEAR       )return

  !do the averaging
  do iou_ix = 1, Niou_ix
     pollwritten = .false.
     do isrc=1,Nsources
        ipoll = lf_src(isrc)%poll
        do k = KMAX_MID-lf_Nvertout+1,KMAX_MID
           do j=1,ljmax
              do i=1,limax
                 xtot=0.0
                 do iix=1,lf_src(isrc)%Nsplit
                    ix=lf_src(isrc)%ix(iix)
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
                    lf_src_tot(i,j,k,ipoll,iou_ix) = lf_src_tot(i,j,k,ipoll,iou_ix) + xtot
                 endif
                 do n=lf_src(isrc)%start, lf_src(isrc)%end
                    lf_src_acc(n,i,j,k,iou_ix)=lf_src_acc(n,i,j,k,iou_ix)+xtot*lf(n,i,j,k)
                 end do
              enddo
           enddo
        enddo
        pollwritten(ipoll) = .true.
      enddo
  enddo
  
  av_fac=av_fac+1
  
  call Add_2timing(NTIMING-9,tim_after,tim_before,"lf: averaging")

end subroutine lf_av

subroutine lf_adv_x(fluxx,i,j,k)
  real, intent(in)::fluxx(NSPEC_ADV,-1:LIMAX+1)
  integer, intent(in)::i,j,k
  real ::x,xn,xx,f_in,inv_tot
  integer ::n,ii,iix,ix,dx,dy,isrc

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
        if(x>1.E-20)then
           n = lf_src(isrc)%start
           do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
              lf(n,i,j,k) = lf(n,i,j,k)*xn ! when dx=-lf_src(isrc)%dist there are no local fractions to transport
              n=n+1
              do dx=-lf_src(isrc)%dist+1,lf_src(isrc)%dist
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n-1,i+1)*x
                 n=n+1
              enddo
           enddo

           if(xx>1.E-20)then
              n = lf_src(isrc)%start
              do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                 do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                    lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+1,i-1)*xx   
                    n=n+1
                 enddo                
                 n=n+1! when dx=lf_src(isrc)%dist there are no local fractions to transport
              enddo
           endif
        else if (xx>1.E-20)then
           n = lf_src(isrc)%start
           do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
              do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n+1,i-1)*xx 
                 n=n+1
              enddo
              lf(n,i,j,k) = lf(n,i,j,k)*xn! when dx=lf_src(isrc)%dist there are no local fractions to transport
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

end subroutine lf_adv_x

subroutine lf_adv_y(fluxy,i,j,k)
  real, intent(in)::fluxy(NSPEC_ADV,-1:LJMAX+1)
  integer, intent(in)::i,j,k
  real ::x,xn,xx,f_in,inv_tot
  integer ::n,jj,iix,ix,dx,dy,isrc

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
        if(x>1.E-20)then
           n = lf_src(isrc)%start
           dy = -lf_src(isrc)%dist
           do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
              lf(n,i,j,k) = lf(n,i,j,k)*xn
              n=n+1
           enddo
           do dy=-lf_src(isrc)%dist+1,lf_src(isrc)%dist
              do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n-Ndiv_rel,j+1)*x
                 n=n+1
              enddo
           enddo
           if(xx>1.E-20)then
              n = lf_src(isrc)%start
              do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                 do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+Ndiv_rel,j-1)*xx
                    n=n+1
                 enddo
              enddo
           endif
        else if (xx>1.E-20)then
           n = lf_src(isrc)%start
           do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
              do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n+Ndiv_rel,j-1)*xx 
                 n=n+1
              enddo
           enddo
           dy=lf_src(isrc)%dist
           do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
              lf(n,i,j,k) = lf(n,i,j,k)*xn
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
  call Add_2timing(NTIMING-7,tim_after,tim_before,"lf: adv_y")

end subroutine lf_adv_y

subroutine lf_adv_k(fluxk,i,j)
    real, intent(in)::fluxk(NSPEC_ADV,KMAX_MID)
    integer, intent(in)::i,j
    real ::x,xn,xx,f_in,inv_tot
    integer ::n,k,iix,ix,dx,dy,isrc
    real loc_frac_src_km1(LF_SRC_TOTSIZE,KMAX_MID-lf_Nvert+2:KMAX_MID)

    call Code_timer(tim_before)
    !need to be careful to always use non-updated values on the RHS
    do k = KMAX_MID-lf_Nvert+2,KMAX_MID
       do n = 1, LF_SRC_TOTSIZE
          loc_frac_src_km1(n,k)=lf(n,i,j,k-1) !NB: k is shifted by 1 in loc_frac_src_km1
       enddo
    enddo
    !loc_frac_src_km1(:,KMAX_MID-lf_Nvert+1)=0.0!Assume zero local fractions coming from above

    do k = KMAX_MID-lf_Nvert+1,KMAX_MID!k is increasing-> can use k+1 to access non-updated value
       do isrc=1,Nsources
          xn=0.0
          x=0.0
          xx=0.0
          !positive x or xx means incoming, negative means outgoing
          do iix=1,lf_src(isrc)%Nsplit
             ix=lf_src(isrc)%ix(iix)
             xn=xn+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
             if(k<KMAX_MID)x=x-dhs1i(k+1)*fluxk(ix,k+1)*lf_src(isrc)%mw(iix)
             xx=xx+dhs1i(k+1)*fluxk(ix,k)*lf_src(isrc)%mw(iix)
          end do
          !NB: here xn already includes the fluxes. Remove them!
          xn=xn-xx-x
          xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux 
          f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
          inv_tot=1.0/(xn+f_in+1.e-40)!incoming dilutes

          x =max(0.0,x)*inv_tot!factor due to flux through bottom face
          xx=max(0.0,xx)*inv_tot!factor due to flux through top face
          if(k<=KMAX_MID-lf_Nvert+1) xx = 0.0 !no fraction coming from higher than KMAX_MID-lf_Nvert+1
          if(k==KMAX_MID) x = 0.0 !no fraction coming from surface
          xn = xn * inv_tot
          !often either x or xx is zero
          if(x>1.E-20)then
             do n = lf_src(isrc)%start, lf_src(isrc)%end
                lf(n,i,j,k) = lf(n,i,j,k)*xn +lf(n,i,j,k+1)*x
             enddo
             if(xx>1.E-20)then
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_km1(n,k)*xx
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
    
    call Add_2timing(NTIMING-6,tim_after,tim_before,"lf: adv_k")
  end subroutine lf_adv_k
 
  subroutine lf_diff(i,j,ds3,ds4,ndiff)
    
    implicit none
    interface 
       subroutine vertdiffn(xn_k,NSPEC,Nij,KMIN_in,SigmaKz,ds3,ds4,ndiff)
         real,intent(inout) :: xn_k(NSPEC,0:*)!dummy
         real,intent(in)::  SigmaKz(*)!dummy
         real,intent(in)::  ds3(*),ds4(*)!dummy
         integer,intent(in)::  NSPEC,ndiff,Nij,KMIN_in
       end subroutine vertdiffn
    end interface

    real, intent(in) :: ds3(2:KMAX_MID),ds4(2:KMAX_MID)
    integer, intent(in) :: i,j,ndiff
    real :: xn_k(LF_SRC_TOTSIZE + Npoll,KMAX_MID),x
    integer ::isec_poll1,isrc
    integer ::k,n,ix,iix,dx,dy
    !how far diffusion should take place above lf_Nvert. 
    ! KUP = 2 gives less than 0.001 differences in locfrac, except sometimes over sea, because
    !ship emission are higher up and need to come down to diminish locfrac
    integer, parameter :: KUP = 2

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
          endif
          xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k) = x
      enddo
    enddo

    call vertdiffn(xn_k,LF_SRC_TOTSIZE+Npoll,1,KMAX_MID-lf_Nvert-KUP,EtaKz(i,j,1,1),ds3,ds4,ndiff)
 
    do k = KMAX_MID-lf_Nvert+1,KMAX_MID
     do isrc=1,Nsources
          x =  1.0/(xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k)+1.E-30)
          do n=lf_src(isrc)%start, lf_src(isrc)%end
             lf(n,i,j,k) = xn_k(n,k)*x                               
          enddo
       enddo
    end do
    call Add_2timing(NTIMING-5,tim_after,tim_before,"lf: diffusion")

end subroutine lf_diff

subroutine lf_chem(i,j)
  !track through chemical reactions
  integer, intent(in) ::i,j
  real :: VOC,HO2,O3,NO,NO2,d_O3,d_NO,d_NO2,d_VOC, k1,k2,J_phot,invt,inv
  real :: SO4,SO2, d_SO2, d_SO4
  integer :: k, n, n_O3,n_NO,n_NO2,n_VOC,nsteps,nsteps1,nsteps2
  integer :: n_SO2,n_SO4,  n_EC_new, n_EC
  real :: k_OH, k_H2O2, k_O3
  real ::  d_age
  
  call Code_timer(tim_before)
  
  if(isrc_EC_f_new>0)then
     do k = KMAX_MID-lf_Nvert+1,KMAX_MID
        ! d_age = amount that has been transformed from EC_f_new to EC_f_age
        d_age = rct(80,k)*xn_adv(ix_EC_f_new,i,j,k) *dt_advec
        inv = 1.0/( xn_adv(ix_EC_f_age,i,j,k) + d_age + 1.0E-20)        
        n_EC_new = lf_src(isrc_EC_f_new)%start
        do n_EC=lf_src(isrc_EC_f_age)%start, lf_src(isrc_EC_f_age)%end
           lf(n_EC,i,j,k) = (lf(n_EC,i,j,k)*xn_adv(ix_EC_f_age,i,j,k) + d_age*lf(n_EC_new,i,j,k)) * inv
           n_EC_new = n_EC_new + 1
        enddo
     enddo
  endif

  if(isrc_EC_f_ffuel_new >0)then
     do k = KMAX_MID-lf_Nvert+1,KMAX_MID
        ! d_age = amount that has been transformed from EC_f_ffuel_new to EC_f_ffuel_age
        d_age = rct(96,k)*xn_adv(ix_EC_f_ffuel_new,i,j,k) *dt_advec
        inv = 1.0/( xn_adv(ix_EC_f_ffuel_age,i,j,k) + d_age + 1.0E-20)        
        n_EC_new = lf_src(isrc_EC_f_ffuel_new)%start
        do n_EC=lf_src(isrc_EC_f_ffuel_age)%start, lf_src(isrc_EC_f_ffuel_age)%end
           lf(n_EC,i,j,k) = (lf(n_EC,i,j,k)*xn_adv(ix_EC_f_ffuel_age,i,j,k) + d_age*lf(n_EC_new,i,j,k)) * inv
           n_EC_new = n_EC_new + 1
        enddo
     enddo
  endif
  if(isrc_EC_f_wood_new >0)then
     do k = KMAX_MID-lf_Nvert+1,KMAX_MID
        d_age = rct(96,k)*xn_adv(ix_EC_f_wood_new,i,j,k) *dt_advec
        inv = 1.0/( xn_adv(ix_EC_f_wood_age,i,j,k) + d_age + 1.0E-20)        
        n_EC_new = lf_src(isrc_EC_f_wood_new)%start
        do n_EC=lf_src(isrc_EC_f_wood_age)%start, lf_src(isrc_EC_f_wood_age)%end
           lf(n_EC,i,j,k) = (lf(n_EC,i,j,k)*xn_adv(ix_EC_f_wood_age,i,j,k) + d_age*lf(n_EC_new,i,j,k)) * inv
           n_EC_new = n_EC_new + 1
        enddo
     enddo     
  endif
  if(isrc_SO2>0)then
     do k = KMAX_MID-lf_Nvert+1,KMAX_MID
        SO4 = xn_2d(NSPEC_SHL+ix_SO4,k)
        n_SO4 = lf_src(isrc_SO4)%start
        stop
        !SO4 produced by SO2 , without emitted SO4:
        !d_SO4 = max(0.0,Dchem(NSPEC_SHL+ix_SO4,k,i,j)-rcemis(NSPEC_SHL+ix_SO4,k))*dt_advec
        inv = 1.0/(SO4 + 1.0E-20)
        
        do n_SO2=lf_src(isrc_SO2)%start, lf_src(isrc_SO2)%end
           lf(n_SO4,i,j,k) = (lf(n_SO4,i,j,k)*(SO4-d_SO4)+d_SO4*lf(n_SO2,i,j,k)) * inv
           n_SO4 = n_SO4 + 1
        enddo
        
    enddo

  endif
  if(isrc_O3<=0)return
  !the source index must give three values, stored at isrc_O3, isrc_NO and isrc_NO2
  do k = KMAX_MID-lf_Nvert+1,KMAX_MID
     !xn_adv is proportional to concentrations (~kg/m3) units in advection routines, (in mass mixing ratio otherwise).
     !xn_2d is in units of molecules/cm3, and defined also with short lived
     !units do not matter for local fractions, as long as all units are the same.
!     if(allocated(Dchem))then
!        !if(me==0.and. i==5 .and.j==5)write(*,*)'DCHEM ',Dchem(NSPEC_SHL+ix_SO2,k,i,j)*0.5*dt_advec,xn_2d(NSPEC_SHL+ix_SO2,k)
!        NO = max(0.0,xn_2d(NSPEC_SHL+ix_NO,k)+Dchem(NSPEC_SHL+ix_NO,k,i,j)*0.5*dt_advec)! put value at approximatively average value during chem timestep
!        NO2 = max(0.0,xn_2d(NSPEC_SHL+ix_NO2,k)+Dchem(NSPEC_SHL+ix_NO2,k,i,j)*0.5*dt_advec)! put value at approximatively average value during chem timestep
!     else
        NO = xn_2d(NSPEC_SHL+ix_NO,k)
        NO2 = xn_2d(NSPEC_SHL+ix_NO2,k)
!     endif
     O3 = xn_2d(NSPEC_SHL+ix_O3,k) 
     VOC = xn_2d(NSPEC_SHL+ix_CH3CO3,k)
     HO2 = xn_2d(ix_HO2,k) !NB: ix_HO2 is already short lived index
     k1 = rct(11,k) * dt_advec
     !CH3CO3 + HO2 -> 0.162*O3 + 0.384*CH3CO3H + 0.454 OH + 0.454 CH3O2 + 0.162 CH3COOH ?
     k2 = 0.162* rct(45,k) * HO2 * dt_advec
     d_VOC = VOC*k2/(1 + k2)
     d_VOC=0.0
     n_O3 = lf_src(isrc_O3)%start
     n_NO = lf_src(isrc_NO)%start
     
     J_phot=rcphot(IDNO2,k)*dt_advec
     !NO2 photodecomposition: NO2+hv->NO+O3 (we skip OP)
     !d_NO2 = -J_phot*NO2
     d_NO2 = NO2*J_phot/(1 +  J_phot)
     !if(i==5 .and. j==5 .and. k==20 .and. me>300  .and. me>400  )write(*,*)d_VOC,d_NO2,rct(1,k),rct(10,k)*NO

     !NO+O3->NO2+O2
     !d_NO = k1*O3*NO
     if(NO<O3)then
        d_NO = (k1*O3*NO)/(1+k1*O3)
     else
        d_NO = (k1*O3*NO)/(1+k1*NO)
     endif

     nsteps = 1
     invt=1.0/nsteps
     do n=1,nsteps
!update fractions
        n_O3 = lf_src(isrc_O3)%start
        n_VOC = lf_src(isrc_VOC)%start
        n_NO = lf_src(isrc_NO)%start
        do n_NO2=lf_src(isrc_NO2)%start, lf_src(isrc_NO2)%end
           if(isnan(lf(n_O3,i,j,k)) .or. lf(n_O3,i,j,k)>10.0)then
              write(*,*)'O3 is nan before chem',O3,NO2,NO,VOC, d_NO2 , d_VOC 
              stop
           endif
           if(isnan(lf(n_NO,i,j,k)))then
              write(*,*)'NO is nan before chem',O3,NO2,NO,VOC
               stop
           endif
           if(isnan(lf(n_NO2,i,j,k)) .or. lf(n_NO2,i,j,k)>10.0)then
              write(*,*)'NO2 is nan before chem',O3,NO2,NO,VOC
               stop
           endif
           lf(n_O3,i,j,k) = (lf(n_O3,i,j,k)*O3 + d_NO2*lf(n_NO2,i,j,k) -d_NO*lf(n_NO,i,j,k) + d_VOC*lf(n_VOC,i,j,k) )/(O3 + d_NO2 - d_NO + d_VOC + 1.0)

           if(isnan(lf(n_O3,i,j,k)) .or. lf(n_O3,i,j,k)>10.0)then
              write(*,*)'O3 is nan after chem ',lf(n_O3,i,j,k),lf(n_NO2,i,j,k),O3,NO2,NO,VOC, d_NO2 , d_VOC 
              stop
           endif

           lf(n_NO,i,j,k) = (lf(n_NO,i,j,k)*NO + d_NO2*lf(n_NO2,i,j,k) )/(NO + d_NO2 + 1.0)

           lf(n_NO2,i,j,k) = (lf(n_NO2,i,j,k)*(NO2-d_NO2) + d_NO*lf(n_NO,i,j,k) + d_NO*lf(n_O3,i,j,k) )/(NO2 - d_NO2 + d_NO + d_NO + 1.0)
           if(isnan(lf(n_NO,i,j,k)) .or. lf(n_NO,i,j,k)>10.0)then
              write(*,*)'NO is nan after chem ',lf(n_NO,i,j,k),O3,NO2,NO,VOC
               stop
           endif
           if(isnan(lf(n_NO2,i,j,k)) .or. lf(n_NO2,i,j,k)>10.0)then
              write(*,*)'NO2 is nan after chem ',lf(n_NO2,i,j,k),O3,NO2,NO,VOC,d_NO2 , d_NO
               stop
           endif
 
           n_O3 = n_O3 + 1
           n_VOC = n_VOC + 1
           n_NO = n_NO + 1
        enddo
     enddo
     
  enddo
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")
end subroutine lf_chem

subroutine lf_aero_pre(i,j) !called just before AerosolEquilib
  integer, intent(in) ::i,j
  integer :: k
  !save concentrations, to see changes
  if(isrc_NH4<0)return;
  call Code_timer(tim_before)
  do k = KMAX_MID-lf_Nvert+1,KMAX_MID
     lf_NH4(k) = xn_2d(NSPEC_SHL+ix_NH4,k)
     lf_NH3(k) = xn_2d(NSPEC_SHL+ix_NH3,k)
  enddo
  call Add_2timing(NTIMING-3,tim_after,tim_before,"lf: chemistry")
  
end subroutine lf_aero_pre

subroutine lf_aero_pos (i,j) !called just after AerosolEquilib
  integer, intent(in) ::i,j
  real :: d_NH4, d_NH3, NH4, NH3, inv
  integer :: n_NH3, n_NH4, k

  if(isrc_NH4<0)return;
  call Code_timer(tim_before)
  do k = KMAX_MID-lf_Nvert+1,KMAX_MID

     NH3 = xn_2d(NSPEC_SHL+ix_NH3,k)
     NH4 = xn_2d(NSPEC_SHL+ix_NH4,k)
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
 
  do isrc=1,Nsources
     if(.not. lf_src(isrc)%DryDep)cycle
     do iix=1,lf_src(isrc)%Nsplit
        ix=lf_src(isrc)%ix(iix)
        ffac = fac*1.e6*lf_src(isrc)%mw(iix) !(units ok?)
        istart = lf_src(isrc)%start 
        iend = lf_src(isrc)%end
        if(isrc==isrc_SO4 .or. isrc==isrc_SO2 .or. lf_src(isrc)%species=="sox")ffac = ffac*32.0/64.0 !SO2->S
        if(isrc==isrc_NH3 .or. isrc==isrc_NH4 .or. lf_src(isrc)%species=="nh3")ffac = ffac* 14.0/17.0!NH3->N
        if(isrc==isrc_NO .or. isrc==isrc_NO2 .or. lf_src(isrc)%species=="nox")ffac = ffac*14.0/46.0 !NO2->N
     
        if( ix==ix_SO4 ) then
           !take directly local fractions from SO4 instead of sox
           istart = lf_src(isrc_SO4)%start
           iend= lf_src(isrc_SO4)%end
        endif
        if( ix==ix_SO2 ) then
           !take directly local fractions from SO2 instead of sox
           istart = lf_src(isrc_SO2)%start
           iend= lf_src(isrc_SO2)%end
        endif
        
        if( ix==ix_NH4 ) then
           istart = lf_src(isrc_NH4)%start
           iend= lf_src(isrc_NH4)%end
        endif
        if( ix==ix_NH3 ) then
           istart = lf_src(isrc_NH3)%start
           iend= lf_src(isrc_NH3)%end
        endif
        
        idep=idep0
        do n = istart, iend 
           idep=idep+1
           loc_frac_drydep(i,j,idep) = loc_frac_drydep(i,j,idep) + lf(n,i,j,KMAX_MID)*DepLoss(ix)*ffac
        enddo
        if(lf_src(isrc)%species=="nox" .and. iix==lf_src(isrc)%Nsplit)then
           !we add also depositions of NO3 and HNO3
           ix=ix_NO3
           idep=idep0
           do n = istart, iend 
              idep=idep+1
              loc_frac_drydep(i,j,idep) = loc_frac_drydep(i,j,idep) + lf(n,i,j,KMAX_MID)*DepLoss(ix)*ffac
           enddo
           ix=ix_HNO3
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
     
        if( ix==ix_SO4 ) then
           !take directly local fractions from SO4 instead of sox
           istart = lf_src(isrc_SO4)%start
           iend= lf_src(isrc_SO4)%end
        endif
        if( ix==ix_SO2 ) then
           !take directly local fractions from SO2 instead of sox
           istart = lf_src(isrc_SO2)%start
           iend= lf_src(isrc_SO2)%end
        endif
        
        if( ix==ix_NH4 ) then
           istart = lf_src(isrc_NH4)%start
           iend= lf_src(isrc_NH4)%end
        endif
        if( ix==ix_NH3 ) then
           istart = lf_src(isrc_NH3)%start
           iend= lf_src(isrc_NH3)%end
        endif
        
        idep=idep0
        do n = istart, iend 
           idep=idep+1
           loc_frac_wetdep(i,j,idep) = loc_frac_wetdep(i,j,idep) + lf(n,i,j,k)*loss*ffac
        enddo
        if(lf_src(isrc)%species=="nox" .and. iix==lf_src(isrc)%Nsplit)then
           !we add also depositions of NO3 and HNO3
           ix=ix_NO3
           idep=idep0
           do n = istart, iend 
              idep=idep+1
              loc_frac_wetdep(i,j,idep) = loc_frac_wetdep(i,j,idep) + lf(n,i,j,k)*loss*ffac
           enddo
           ix=ix_HNO3
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

subroutine lf_emis(indate)
!include emission contributions to local fractions

!NB: should replace most of the stuff and use gridrcemis instead!

  implicit none
  type(date), intent(in) :: indate  ! Gives year..seconds
  integer :: i, j, k, n       ! coordinates, loop variables
  integer :: icc, ncc         ! No. of countries in grid.
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILE

  ! Save daytime value between calls, initialise to zero
  integer, save, dimension(MAXNLAND) ::  daytime(1:MAXNLAND) = 0  !  0=night, 1=day
  integer, save, dimension(MAXNLAND) ::  localhour(1:MAXNLAND) = 1  ! 1-24 local hour in the different countries
  integer                         ::  hourloc      !  local hour 
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions
  real ::  s       ! source term (emis) before splitting
  integer ::icc_lf, iqrc, itot
  integer, save :: wday , wday_loc ! wday = day of the week 1-7
  integer ::ix,iy,iix, iiix,dx, dy, isec_poll, iisec_poll, isec_poll1, isrc, ic, is
  real::dt_lf, xtot, x
  real :: lon
  integer :: jmin,jmax,imin,imax,n0
  
  call Code_timer(tim_before)
 
  do j = lj0,lj1
    do i = li0,li1
       ix=(gi0+i-2)/((GIMAX+Ndiv_coarse-1)/Ndiv_coarse)+1 !i coordinate in coarse domain
       iy=(gj0+j-2)/((GJMAX+Ndiv_coarse-1)/Ndiv_coarse)+1 !j coordinate in coarse domain
       do isrc=1,Nsources
         iem = lf_src(isrc)%iem
         isec = lf_src(isrc)%sector
         do k=max(KEMISTOP,KMAX_MID-lf_Nvert+1),KMAX_MID
            if(lf_src(isrc)%iqrc<0)then
               if(lf_emis_tot(i,j,k,lf_src(isrc)%poll)<1.E-20)cycle
            else
               if(emis_lf(i,j,k,isrc)<1.E-20)cycle
            endif
            xtot=0.0
            do iix=1,lf_src(isrc)%Nsplit
               iiix=lf_src(isrc)%ix(iix)
               xtot=xtot+(xn_adv(iiix,i,j,k)*lf_src(isrc)%mw(iix))*(dA(k)+dB(k)*ps(i,j,1))/ATWAIR/GRAV
            end do
            if(lf_src(isrc)%type=='country' .and. (Ncountry_lf>0 .or. Ncountry_group_lf>0))then
               n0=lf_src(isrc)%start
               do ic=1,Ncountry_lf
                  do is=1,Ncountrysectors_lf
                     lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot+emis_lf_cntry(i,j,k,ic,is,isrc))/(xtot+lf_emis_tot(i,j,k,lf_src(isrc)%poll)+1.e-20)
                     n0=n0+1
                 enddo
               enddo
               do ic=Ncountry_lf+1,Ncountry_lf+Ncountry_group_lf
                  do is=1,Ncountrysectors_lf
                     lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot+emis_lf_cntry(i,j,k,ic,is,isrc))/(xtot+lf_emis_tot(i,j,k,lf_src(isrc)%poll)+1.e-20)
                     n0=n0+1
                  enddo
               enddo
               cycle !only one fraction per country
            else if(lf_src(isrc)%type=='relative' .or. lf_src(isrc)%type=='coarse')then
               !Country constraints already included in emis_lf
               if(lf_src(isrc)%type=='relative') n0 = lf_src(isrc)%start + (lf_src(isrc)%Npos - 1)/2 !"middle" point is dx=0 dy=0
               if(lf_src(isrc)%type=='coarse') n0 = lf_src(isrc)%start+ix-1+(iy-1)*Ndiv_coarse 
               lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot+emis_lf(i,j,k,isrc))/(xtot+lf_emis_tot(i,j,k,lf_src(isrc)%poll)+1.e-20)
            else
               if(me==0)write(*,*)'LF type not recognized)'
               stop
            endif
            do n = lf_src(isrc)%start, lf_src(isrc)%end
               if(n==n0)cycle  !counted above
               lf(n,i,j,k)=(lf(n,i,j,k)*xtot)/(xtot+lf_emis_tot(i,j,k,lf_src(isrc)%poll)+1.e-20)!fractions are diluted
            enddo
         enddo
      enddo
   end do ! i
 end do ! j

  call Add_2timing(NTIMING-4,tim_after,tim_before,"lf: emissions")

end subroutine lf_emis

subroutine add_lf_emis(s,i,j,iem,isec,iland)
  real, intent(in) :: s
  integer, intent(in) :: i,j,iem,isec,iland
  integer :: n, isrc, k, ipoll,ic,is,ig
  real :: emis
  integer :: ngroups, ig2ic(Max_Country_groups)

  call Code_timer(tim_before)

  do n=1,iem2Nipoll(iem)
     ipoll = iem2ipoll(iem,n)
     if(ipoll2iqrc(ipoll)>0)then
        !only extract that single pollutant
        do k=max(KEMISTOP,KMAX_MID-lf_Nvert+1),KMAX_MID
           lf_emis_tot(i,j,k,ipoll) = lf_emis_tot(i,j,k,ipoll) + s * emisfrac(ipoll2iqrc(ipoll),sec2split_map(isec),iland)&
                * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) * dt_advec!total over all sectors and countries for each pollutant
        enddo
     else
        do k=max(KEMISTOP,KMAX_MID-lf_Nvert+1),KMAX_MID
           lf_emis_tot(i,j,k,ipoll) = lf_emis_tot(i,j,k,ipoll) + s * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) * dt_advec!total for each pollutant
        enddo
     endif
  enddo

  do isrc = 1, Nsources
     if(lf_src(isrc)%iem /= iem) cycle
     if(Ncountry_lf>0)then
        !has to store more detailed info
        do ic=1,Ncountry_lf
           if(country_ix_list(ic)==IC_TMT.and.(iland==IC_TM.or.iland==IC_TME))then
           else if(country_ix_list(ic)==IC_AST.and.(iland==IC_ASM.or.iland==IC_ASE.or.iland==IC_ARE.or.iland==IC_ARL.or.iland==IC_CAS))then
           else if(country_ix_list(ic)==IC_UZT.and.(iland==IC_UZ.or.iland==IC_UZE))then
           else if(country_ix_list(ic)==IC_KZT.and.(iland==IC_KZ.or.iland==IC_KZE))then
           else if(country_ix_list(ic)==IC_RUE.and.(iland==IC_RU.or.iland==IC_RFE.or.iland==IC_RUX))then
           else if(country_ix_list(ic)/=iland)then
              cycle
           endif
           do is=1,Ncountrysectors_lf
              if(lf_country_sector_list(is)/=isec .and. lf_country_sector_list(is)/=0)cycle
              emis = s * dt_advec
              if(lf_src(isrc)%iqrc>0)emis = emis *emisfrac(lf_src(isrc)%iqrc,sec2split_map(isec),iland)
              do k=max(KEMISTOP,KMAX_MID-lf_Nvert+1),KMAX_MID
                 emis_lf_cntry(i,j,k,ic,is,isrc)=emis_lf_cntry(i,j,k,ic,is,isrc) + emis * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))
              enddo
           enddo
        enddo
     end if
     if(Ncountry_group_lf>0)then
        !has to store more detailed info
        ngroups = 0
        do ic=1,Ncountry_group_lf
           !find all groups that include iland
           if(any(lf_country_group(ic)%ix(:)==iland))then
              ngroups = ngroups + 1
              ig2ic(ngroups) = ic + Ncountry_lf
           endif
        enddo
        do ig = 1, ngroups
           ic = ig2ic(ig) ! index for country_group
           do is=1,Ncountrysectors_lf
              if(lf_country_sector_list(is)/=isec .and. lf_country_sector_list(is)/=0)cycle
              emis = s * dt_advec
              if(lf_src(isrc)%iqrc>0)emis = emis *emisfrac(lf_src(isrc)%iqrc,sec2split_map(isec),iland)
              do k=max(KEMISTOP,KMAX_MID-lf_Nvert+1),KMAX_MID
                 emis_lf_cntry(i,j,k,ic,is,isrc)=emis_lf_cntry(i,j,k,ic,is,isrc) + emis * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))
              enddo
           enddo
        enddo
     endif
     if(lf_src(isrc)%iqrc>0 .and. (Ncountry_lf>0.or.Ncountry_group_lf>0))then
        !sum of all emissions for that species from all countries and sectors
        do k=max(KEMISTOP,KMAX_MID-lf_Nvert+1),KMAX_MID
           emis_lf(i,j,k,isrc) = emis_lf(i,j,k,isrc) + s * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) &
                * emisfrac(lf_src(isrc)%iqrc,sec2split_map(isec),iland) *dt_advec
        enddo
        
     else
        if(lf_src(isrc)%sector /= isec .and. lf_src(isrc)%sector /= 0) cycle
        if(lf_src(isrc)%country_ix>0 .and. lf_src(isrc)%country_ix/=iland) cycle
        if(lf_src(isrc)%iqrc>0)then
           !single pollutant, part of emitted group of pollutant
           ipoll = lf_src(isrc)%poll
           do k=max(KEMISTOP,KMAX_MID-lf_Nvert+1),KMAX_MID
              emis_lf(i,j,k,isrc) = emis_lf(i,j,k,isrc) + s * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) &
                   * emisfrac(lf_src(isrc)%iqrc,sec2split_map(isec),iland) *dt_advec       
           enddo
        else
           do k=max(KEMISTOP,KMAX_MID-lf_Nvert+1),KMAX_MID
              emis_lf(i,j,k,isrc) = emis_lf(i,j,k,isrc) + s * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) * dt_advec
           enddo
        endif
     endif
  enddo
  call Add_2timing(NTIMING-4,tim_after,tim_before,"lf: emissions")
  
end subroutine add_lf_emis

end module LocalFractions_mod
