! <uEMEP_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module uEMEP_mod
!
! all subroutines for uEMEP
!
use CheckStop_mod,     only: CheckStop,StopAll
use Chemfields_mod,    only: xn_adv
use ChemDims_mod,      only: NSPEC_ADV, NSPEC_SHL,NEMIS_File
use ChemSpecs_mod,     only: species_adv,species
use Country_mod,       only: MAXNLAND,NLAND,Country
use EmisDef_mod,       only: loc_frac, loc_frac_1d, loc_frac_hour, loc_tot_hour,&
                            loc_frac_hour_inst, loc_tot_hour_inst,loc_frac_day, &
                            loc_tot_day, loc_frac_month,loc_tot_month,&
                            loc_frac_full,loc_tot_full, NSECTORS,EMIS_FILE, &
                            nlandcode,landcode,sec2tfac_map,sec2hfac_map, &
                            ISNAP_DOM,secemis, roaddust_emis_pot,KEMISTOP,&
                            NEmis_sources, Emis_source_2D, Emis_source
use EmisGet_mod,       only: nrcemis, iqrc2itot, emis_nsplit,nemis_kprofile, emis_kprofile,&
                             make_iland_for_time,itot2iqrc,iqrc2iem
use GridValues_mod,    only: dA,dB,xm2, dhs1i, glat, glon, projection, extendarea_N
use MetFields_mod,     only: ps,roa,EtaKz
use Config_module,     only: KMAX_MID, KMAX_BND,USES, USE_uEMEP, uEMEP, IOU_HOUR&
                             , IOU_HOUR_INST,IOU_INST,IOU_YEAR,IOU_MON,IOU_DAY&
                             ,IOU_HOUR,IOU_HOUR_INST, KMAX_MID &
                             ,MasterProc,dt_advec, RUNDOMAIN, runlabel1 &
                             ,HOURLYFILE_ending
use MPI_Groups_mod
use NetCDF_mod,        only: Real4,Out_netCDF
use OwnDataTypes_mod,  only: Deriv, Npoll_uemep_max, Nsector_uemep_max, TXTLEN_FILE
use Par_mod,           only: me,LIMAX,LJMAX,MAXLIMAX,MAXLJMAX,gi0,gj0,li0,li1,lj0,lj1,GIMAX,GJMAX
use PhysicalConstants_mod, only : GRAV, ATWAIR 
use SmallUtils_mod,    only: find_index
use TimeDate_mod,      only: date, current_date,day_of_week
use TimeDate_ExtraUtil_mod,only: date2string
use Timefactors_mod,   only: &
    DegreeDayFactors       & ! degree-days used for SNAP-2
    ,Gridded_SNAP2_Factors, gridfac_HDD & 
    ,GridTfac &!array with monthly gridded time factors
    ,fac_min,timefactors   &                  ! subroutine
    ,fac_ehh24x7 ,fac_emm, fac_edd, timefac  ! time-factors
use My_Timing_mod,     only: Add_2timing, Code_timer, NTIMING

!(dx,dy,i,j) shows contribution of pollutants from (i+dx,j+dy) to (i,j)

implicit none
!external advection_mod_mp_vertdiffn_k

private

public  :: init_uEMEP
public  :: out_uEMEP
public  :: av_uEMEP
public  :: uemep_adv_x
public  :: uemep_adv_y
public  :: uemep_adv_k
public  :: uemep_diff
public  :: uemep_emis

integer, public, save :: uEMEP_Size1=0 !total size of the first 3 dimensions of loc_frac
!uEMEP_Size1=uEMEP%Nsec_poll*(2*uEMEP%dist+1)**2
integer, public, save :: uEMEP_Sizedxdy=0 !total size of the first 3 dimensions of loc_frac
!uEMEP_Sizedxdy=(2*uEMEP%dist+1)**2

real, private, save ::av_fac_hour,av_fac_day,av_fac_month,av_fac_full
real, allocatable, save ::loc_poll_to(:,:,:,:,:)

logical, public, save :: COMPUTE_LOCAL_TRANSPORT=.false.
integer , private, save :: uEMEPNvertout = 1!number of vertical levels to save in output
integer, public, save :: NTIMING_uEMEP=7
real, private :: tim_after,tim_before

contains
subroutine init_uEMEP
  integer :: i, ix, itot, iqrc, iem, iemis, isec, ipoll, ixnh3, ixnh4

  call Code_timer(tim_before)
  uEMEP%Nsec_poll = 0
  uEMEP%Npoll = 0
  do iemis=1,Npoll_uemep_max
     if(uEMEP%poll(iemis)%emis=='none')then
        call CheckStop(iemis==1,"init_uEMEP: no pollutant specified")
        exit
     else
        uEMEP%Npoll = uEMEP%Npoll + 1
        uEMEP%poll(iemis)%Nsectors = 0
        uEMEP%poll(iemis)%sec_poll_ishift=uEMEP%Nsec_poll
        do isec=1,Nsector_uemep_max
           if(uEMEP%poll(iemis)%sector(isec)<0)then
              call CheckStop(isec==0,"init_uEMEP: nosector specified for "//uEMEP%poll(iemis)%emis)
              exit
           else
              uEMEP%Nsec_poll = uEMEP%Nsec_poll + 1
              uEMEP%poll(iemis)%Nsectors = uEMEP%poll(iemis)%Nsectors +1
           endif
        enddo
     endif
  enddo

!find indices in EMIS_File
  do ipoll=1,uEMEP%Npoll        
     do iem=1,NEMIS_FILE 
        if(trim(EMIS_File(iem))==trim(uEMEP%poll(ipoll)%emis))then
           uEMEP%poll(ipoll)%EMIS_File_ix = iem
           exit
        endif
     enddo
     call CheckStop(iem>NEMIS_FILE,"uemep pollutant not found: "//trim(uEMEP%poll(ipoll)%emis))
  enddo
 
!  uEMEP%dist = 5
!  uEMEP%Nvert =7
  do i=1,4
     if(uEMEP%DOMAIN(i)<0)uEMEP%DOMAIN(i) = RUNDOMAIN(i)
  enddo

  uEMEP_Sizedxdy = (2*uEMEP%dist+1)**2
  uEMEP_Size1 = uEMEP%Nsec_poll*uEMEP_Sizedxdy

  if(MasterProc)then
     write(*,*)'uEMEP pollutants : ',uEMEP%Npoll
     write(*,*)'total uEMEP pollutants and sectors : ',uEMEP%Nsec_poll
  end if
  do ipoll=1,uEMEP%Npoll        
     iem=find_index(uEMEP%poll(ipoll)%emis ,EMIS_FILE(1:NEMIS_FILE))
     call CheckStop( iem<1, "uEMEP did not find corresponding emission file: "//trim(uEMEP%poll(ipoll)%emis) )
     call CheckStop( iem/=uEMEP%poll(ipoll)%EMIS_File_ix, "uEMEP wrong emis file index for: "//trim(uEMEP%poll(ipoll)%emis) )
     uEMEP%poll(ipoll)%Nix=emis_nsplit(iem)
     do i=1,uEMEP%poll(ipoll)%Nix
        iqrc=sum(emis_nsplit(1:iem-1)) + i
        itot=iqrc2itot(iqrc)
        ix=itot-NSPEC_SHL
        uEMEP%poll(ipoll)%ix(i)=ix
        uEMEP%poll(ipoll)%mw(i)=species_adv(ix)%molwt
        if(uEMEP%poll(ipoll)%emis=="nox ")then
           ix=find_index("NO2",species_adv(:)%name)
           call CheckStop(ix<0,'Index for NO2 not found')
           uEMEP%poll(ipoll)%mw(i)=species_adv(ix)%molwt
        endif
        if(uEMEP%poll(ipoll)%emis=="sox ")then
           ix=find_index("SO2",species_adv(:)%name)
           call CheckStop(ix<0,'Index for SO2 not found')
           uEMEP%poll(ipoll)%mw(i)=species_adv(ix)%molwt
        endif

!!$     if(uEMEP%poll(ipoll)%emis=="nox ")then
!!$        uEMEP%poll(ipoll)%Nix=0
!!$        ixnh4=find_index("NH4_F",species_adv(:)%name)
!!$        ixnh3=find_index("NH3",species_adv(:)%name)
!!$        do ix=1,NSPEC_ADV
!!$           if(ix==ixnh4.or.ix==ixnh3)cycle!reduced nitrogen
!!$           if(species_adv(ix)%nitrogens>0)then
!!$              uEMEP%poll(ipoll)%ix(uEMEP%poll(ipoll)%Nix)=ix
!!$              uEMEP%poll(ipoll)%Nix =  uEMEP%poll(ipoll)%Nix + 1
!!$              if(species_adv(ix)%nitrogens==1)uEMEP%poll(ipoll)%mw(uEMEP%poll(ipoll)%Nix)=46
!!$              if(species_adv(ix)%nitrogens==2)uEMEP%poll(ipoll)%mw(uEMEP%poll(ipoll)%Nix)=92
!!$           endif
!!$        enddo
!!$     endif
     if(uEMEP%poll(ipoll)%emis=="nh3 ")then
        uEMEP%poll(ipoll)%Nix=0
        ixnh4=find_index("NH4_F",species_adv(:)%name , any_case=.true.)
        ixnh3=find_index("NH3",species_adv(:)%name)
        do ix=1,NSPEC_ADV
           if(ix/=ixnh4.and.ix/=ixnh3)cycle!not reduced nitrogen
           if(species_adv(ix)%nitrogens>0)then
              uEMEP%poll(ipoll)%Nix =  uEMEP%poll(ipoll)%Nix + 1
              uEMEP%poll(ipoll)%ix(uEMEP%poll(ipoll)%Nix)=ix
              uEMEP%poll(ipoll)%mw(uEMEP%poll(ipoll)%Nix)=species_adv(ixnh3)%molwt!use NH3 mw also for NH4
           endif
        enddo
     endif
     end do
     if(MasterProc)then
        write(*,*)'uEMEP pollutant : ',uEMEP%poll(ipoll)%emis
        write(*,*)'uEMEP number of species in '//trim(uEMEP%poll(ipoll)%emis)//' group: ',uEMEP%poll(ipoll)%Nix
        write(*,"(A,30(A,F6.2))")'including:',('; '//trim(species_adv(uEMEP%poll(ipoll)%ix(i))%name)//', mw ',uEMEP%poll(ipoll)%mw(i),i=1,uEMEP%poll(ipoll)%Nix)
        write(*,"(A,30I4)")'sectors:',(uEMEP%poll(ipoll)%sector(i),i=1,uEMEP%poll(ipoll)%Nsectors)
        write(*,"(A,30I4)")'ix:',(uEMEP%poll(ipoll)%ix(i),i=1,uEMEP%poll(ipoll)%Nix)
     end if
  end do

  COMPUTE_LOCAL_TRANSPORT = uEMEP%COMPUTE_LOCAL_TRANSPORT

  av_fac_hour=0.0
  av_fac_day=0.0
  av_fac_month=0.0
  av_fac_full=0.0
  
  allocate(loc_frac(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
  loc_frac=0.0 !must be initiated to 0 so that outer frame does not contribute.
  if(COMPUTE_LOCAL_TRANSPORT)then
  allocate(loc_poll_to(-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEPNvertout+1:KMAX_MID))
  loc_poll_to=0.0 
  endif
  allocate(loc_frac_1d(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,0:max(LIMAX,LJMAX)+1))        
  if(uEMEP%HOUR)then
     allocate(loc_frac_hour(-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Nsec_poll))
     loc_frac_hour=0.0
     allocate(loc_tot_hour(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Npoll))
     loc_tot_hour=0.0
  endif  
  if(uEMEP%HOUR_INST)then
     allocate(loc_frac_hour_inst(-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Nsec_poll))
     loc_frac_hour_inst=0.0
     allocate(loc_tot_hour_inst(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Npoll))
     loc_tot_hour_inst=0.0
  endif  
  if(uEMEP%DAY)then
     allocate(loc_frac_day(-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Nsec_poll))
     loc_frac_day=0.0
     allocate(loc_tot_day(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Npoll))
     loc_tot_day=0.0
  endif
  if(uEMEP%MONTH)then
     allocate(loc_frac_month(-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Nsec_poll))
     loc_frac_month=0.0
     allocate(loc_tot_month(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Npoll))
     loc_tot_month=0.0
  endif
  if(uEMEP%YEAR)then
     allocate(loc_frac_full(-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Nsec_poll))
     loc_frac_full=0.0
     allocate(loc_tot_full(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,uEMEP%Npoll))
     loc_tot_full=0.0
  else
     !need to be allocated to avoid debugging error
     allocate(loc_frac_full(1,1,1,1,1,1))
     allocate(loc_tot_full(1,1,1,1))
  endif
  
  call Add_2timing(NTIMING-9,tim_after,tim_before,"uEMEP: init")
end subroutine init_uEMEP


subroutine out_uEMEP(iotyp)
  integer, intent(in) :: iotyp
  character(len=200) ::filename, varname
  real :: xtot,scale,invtot,t1,t2
  integer ::i,j,k,dx,dy,ix,iix,isec,iisec,isec_poll,ipoll,isec_poll1
  integer ::ndim,kmax,CDFtype,dimSizes(10),chunksizes(10)
  integer ::ndim_tot,dimSizes_tot(10),chunksizes_tot(10)
  character (len=20) ::dimNames(10),dimNames_tot(10)
  type(Deriv) :: def1 ! definition of fields
  type(Deriv) :: def2 ! definition of fields
  logical ::overwrite
  logical,save :: first_call(10)=.true.
  real,allocatable ::tmp_ext(:,:,:,:,:)!allocate since it may be heavy for the stack
  type(date) :: onesecond = date(0,0,0,0,1)
  character(len=TXTLEN_FILE),save :: oldhourlyname = 'NOTSET'
  character(len=TXTLEN_FILE),save :: oldhourlyInstname = 'NOTSET'
  character(len=TXTLEN_FILE),save :: oldmonthlyname
  real :: fracsum(LIMAX,LJMAX)

  call Code_timer(tim_before)

  if(COMPUTE_LOCAL_TRANSPORT)allocate(tmp_ext(-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,1-uEMEP%dist:LIMAX+uEMEP%dist,1-uEMEP%dist:LJMAX+uEMEP%dist,KMAX_MID-uEMEPNvertout+1:KMAX_MID))
  if(iotyp==IOU_HOUR_INST .and. uEMEP%HOUR_INST)then
     fileName = trim(runlabel1)//'_uEMEP_hourInst'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
     if(oldhourlyInstname/=fileName)then
        first_call(iotyp) = .true.
        oldhourlyInstname = fileName
     endif
  else if(iotyp==IOU_HOUR .and. uEMEP%HOUR)then
     fileName = trim(runlabel1)//'_uEMEP_hour'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
     if(oldhourlyname/=fileName)then
        first_call(iotyp) = .true.
        oldhourlyname = fileName
     endif
  else if(iotyp==IOU_DAY .and. uEMEP%DAY)then
     fileName=trim(runlabel1)//'_uEMEP_day.nc'
  else if(iotyp==IOU_MON .and. uEMEP%MONTH)then
     if(uEMEP%MONTH_ENDING /= "NOTSET")then
        fileName=trim(runlabel1)//'_uEMEP_month'//date2string(trim(uEMEP%MONTH_ENDING),current_date,-1.0)
        if(oldmonthlyname/=fileName)then
           first_call(iotyp) = .true.
           oldmonthlyname = fileName
        endif
     else
        fileName=trim(runlabel1)//'_uEMEP_month.nc'
     endif
  else if(iotyp==IOU_YEAR .and. uEMEP%YEAR)then
     fileName=trim(runlabel1)//'_uEMEP_full.nc'
  else
     return
  endif
  ndim=5
  ndim_tot=3
  kmax=uEMEPNvertout
  scale=1.0
  CDFtype=Real4
  !  dimSizes(1)=uEMEP%Nsec_poll
  !  dimNames(1)='sector'
  dimSizes(1)=2*uEMEP%dist+1
  dimNames(1)='x_dist'
  dimSizes(2)=2*uEMEP%dist+1
  dimNames(2)='y_dist'
  dimSizes(3)=LIMAX
  dimSizes(4)=LJMAX

  dimSizes_tot(1)=LIMAX
  dimSizes_tot(2)=LJMAX

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
  chunksizes=1
  chunksizes(3)=dimSizes(3)
  chunksizes(4)=dimSizes(4)
  chunksizes(5)=dimSizes(5)
  chunksizes_tot=1
  chunksizes_tot(1)=dimSizes_tot(1)
  chunksizes_tot(2)=dimSizes_tot(2)
  chunksizes_tot(3)=dimSizes_tot(3)

  isec_poll1=1
  overwrite=.true.!only first time when new file is created
  do ipoll=1,uEMEP%Npoll
     def2%name=trim(uEMEP%poll(ipoll)%emis)
     if(first_call(iotyp))then
        do iisec=1,uEMEP%poll(ipoll)%Nsectors
           isec=uEMEP%poll(ipoll)%sector(iisec)
           write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_fraction'
           if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_fraction'
           call Out_netCDF(iotyp,def1,ndim,kmax,loc_frac_full,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=overwrite,create_var_only=.true.,chunksizes=chunksizes)
           overwrite=.false.
           if(isec==0 .and. COMPUTE_LOCAL_TRANSPORT)then
           write(def2%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_transport'
           if(isec==0)write(def2%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_transport'
           if(me==0)write(*,*)'making '//def2%name
           call Out_netCDF(iotyp,def2,ndim,kmax,loc_frac_full,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=.false.,create_var_only=.true.,chunksizes=chunksizes)  
           endif
           if(iisec==1)then
              def1%name=trim(uEMEP%poll(ipoll)%emis)//'_fracsum'
              if(me==0)write(*,*)'making '//def1%name//' with dimsizes ',(dimSizes_tot(i),i=1,ndim_tot)
              call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.true.,chunksizes=chunksizes_tot)  
           endif
        enddo
        
        def2%name=trim(uEMEP%poll(ipoll)%emis)
        call Out_netCDF(iotyp,def2,ndim_tot,kmax,loc_tot_full,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
             fileName_given=trim(fileName),overwrite=overwrite,create_var_only=.true.,chunksizes=chunksizes_tot)
  
     endif

     if(iotyp==IOU_HOUR_INST)then
        !compute instantaneous values
        !need to transpose loc_frac. Could be avoided?
        do iisec=1,uEMEP%poll(ipoll)%Nsectors
           isec_poll=isec_poll1+iisec-1
           isec=uEMEP%poll(ipoll)%sector(iisec)
           if(iisec==1)fracsum=0.0
           do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                    xtot=0.0
                    do iix=1,uEMEP%poll(ipoll)%Nix
                       ix=uEMEP%poll(ipoll)%ix(iix)
                       xtot=xtot+(xn_adv(ix,i,j,k)*uEMEP%poll(ipoll)%mw(iix))/ATWAIR&
                            *roa(i,j,k,1)*1.E9 !for ug/m3
                       !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                    end do
                    loc_tot_hour_inst(i,j,k,ipoll)=xtot
                    do dy=-uEMEP%dist,uEMEP%dist
                       do dx=-uEMEP%dist,uEMEP%dist
                          loc_frac_hour_inst(dx,dy,i,j,k,isec_poll)=loc_frac(isec_poll,dx,dy,i,j,k)
                          if(iisec==1.and.k==KMAX_MID)fracsum(i,j)=fracsum(i,j)+loc_frac_hour_inst(dx,dy,i,j,k,isec_poll)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           scale=1.0
           write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_fraction'
           if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_fraction'
           call Out_netCDF(iotyp,def1,ndim,kmax,loc_frac_hour_inst(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.) 

           if(iisec==1)then
              def1%name=trim(uEMEP%poll(ipoll)%emis)//'_fracsum'
              call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)  
           endif

           if(isec==0 .and. COMPUTE_LOCAL_TRANSPORT)then
              !loc_frac_hour_instare fractions -> convert first to pollutant                          
              do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
                 do j=1,ljmax
                    do i=1,limax
                       do dy=-uEMEP%dist,uEMEP%dist
                          do dx=-uEMEP%dist,uEMEP%dist
                             !in  loc_poll_to (dx,dy,i,j) shows from (i,j) to (i+dx,j+dy)
                             loc_frac_hour_inst(dx,dy,i,j,k,isec_poll)=loc_frac_hour_inst(dx,dy,i,j,k,isec_poll)*loc_tot_hour_inst(i,j,k,ipoll)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              call extendarea_N(loc_frac_hour_inst(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),tmp_ext,uEMEP%dist,uEMEP_Sizedxdy,uEMEPNvertout)

           do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                    do dy=-uEMEP%dist,uEMEP%dist
                       do dx=-uEMEP%dist,uEMEP%dist
                          !in  loc_poll_to (dx,dy,i,j) shows from (i,j) to (i+dx,j+dy)
                          loc_poll_to(dx,dy,i,j,k)=tmp_ext(-dx,-dy,i+dx,j+dy,k)!tmp_ext already converted to pollutant                          
                       enddo
                    enddo
                 enddo
              enddo
           enddo

           write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_transport'
           if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_transport'
           scale=1.0
           call Out_netCDF(iotyp,def1,ndim,kmax,loc_poll_to,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)   
       endif
       enddo

           scale=1.0
           def1%name=trim(uEMEP%poll(ipoll)%emis)
           call Out_netCDF(iotyp,def1,ndim_tot,kmax,loc_tot_hour_inst(1,1,KMAX_MID-uEMEPNvertout+1,ipoll),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.) 
     else if(iotyp==IOU_HOUR )then
        do iisec=1,uEMEP%poll(ipoll)%Nsectors
           isec_poll=isec_poll1+iisec-1
           isec=uEMEP%poll(ipoll)%sector(iisec)

           if(iisec==1)fracsum=0.0
           !copy before dividing by loc_tot_hour
           if(COMPUTE_LOCAL_TRANSPORT)call extendarea_N(loc_frac_hour(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),tmp_ext,uEMEP%dist,uEMEP_Sizedxdy,uEMEPNvertout)

           do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                    invtot=1.0/(loc_tot_hour(i,j,k,ipoll)+1.E-20)
                    do dy=-uEMEP%dist,uEMEP%dist
                       do dx=-uEMEP%dist,uEMEP%dist
                          loc_frac_hour(dx,dy,i,j,k,isec_poll)=loc_frac_hour(dx,dy,i,j,k,isec_poll)*invtot
                          if(iisec==1.and.k==KMAX_MID)fracsum(i,j)=fracsum(i,j)+loc_frac_hour(dx,dy,i,j,k,isec_poll)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           scale=1.0
           write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_fraction'
           if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_fraction'
           
           call Out_netCDF(iotyp,def1,ndim,kmax,loc_frac_hour(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.) 

           if(iisec==1)then
              def1%name=trim(uEMEP%poll(ipoll)%emis)//'_fracsum'
              call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)  
           endif

           if(isec==0 .and. COMPUTE_LOCAL_TRANSPORT)then
              do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
                 do j=1,ljmax
                    do i=1,limax
                       do dy=-uEMEP%dist,uEMEP%dist
                          do dx=-uEMEP%dist,uEMEP%dist
                             !in  loc_poll_to (dx,dy,i,j) shows from (i,j) to (i+dx,j+dy)
                             loc_poll_to(dx,dy,i,j,k)=tmp_ext(-dx,-dy,i+dx,j+dy,k)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_transport'
              if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_transport'
              if(abs(av_fac_hour)>1.E-5)then
                 scale=1.0/av_fac_hour
              else
                 scale=0.0
              endif
              if(abs(av_fac_hour)<1.E-5)scale=0.0
              call Out_netCDF(iotyp,def1,ndim,kmax,loc_poll_to,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
           endif
        enddo
        !              loc_tot_hour=loc_tot_hour/av_fac_hour
        if(abs(av_fac_hour)>1.E-5)then
           scale=1.0/av_fac_hour
        else
           scale=0.0
        endif
        def1%name=trim(uEMEP%poll(ipoll)%emis)
        call Out_netCDF(iotyp,def1,ndim_tot,kmax,loc_tot_hour(1,1,KMAX_MID-uEMEPNvertout+1,ipoll),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
             fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.) 

     else  if(iotyp==IOU_DAY)then
        
        do iisec=1,uEMEP%poll(ipoll)%Nsectors
           isec_poll=isec_poll1+iisec-1
           isec=uEMEP%poll(ipoll)%sector(iisec)
           !copy before dividing by loc_tot
           if(COMPUTE_LOCAL_TRANSPORT)call extendarea_N(loc_frac_day(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),tmp_ext,uEMEP%dist,uEMEP_Sizedxdy,uEMEPNvertout)
           if(iisec==1)fracsum=0.0
           do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                    invtot=1.0/(loc_tot_day(i,j,k,ipoll)+1.E-20)
                    do dy=-uEMEP%dist,uEMEP%dist
                       do dx=-uEMEP%dist,uEMEP%dist
                          loc_frac_day(dx,dy,i,j,k,isec_poll)=loc_frac_day(dx,dy,i,j,k,isec_poll)*invtot
                          if(iisec==1.and.k==KMAX_MID)fracsum(i,j)=fracsum(i,j)+loc_frac_day(dx,dy,i,j,k,isec_poll)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           scale=1.0
           write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_fraction'
           if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_fraction'

           call Out_netCDF(iotyp,def1,ndim,kmax,loc_frac_day(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)  


           if(iisec==1)then
              def1%name=trim(uEMEP%poll(ipoll)%emis)//'_fracsum'
              call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)  
           endif

           if(isec==0 .and. COMPUTE_LOCAL_TRANSPORT)then
              do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
                 do j=1,ljmax
                    do i=1,limax
                       do dy=-uEMEP%dist,uEMEP%dist
                          do dx=-uEMEP%dist,uEMEP%dist
                             !in  loc_poll_to (dx,dy,i,j) shows from (i,j) to (i+dx,j+dy)
                             loc_poll_to(dx,dy,i,j,k)=tmp_ext(-dx,-dy,i+dx,j+dy,k)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_transport'
              if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_transport'
              if(abs(av_fac_day)>1.E-5)then
                 scale=1.0/av_fac_day
              else
                 scale=0.0
              endif
              call Out_netCDF(iotyp,def1,ndim,kmax,loc_poll_to,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
           endif
        enddo

        if(abs(av_fac_day)>1.E-5)then
           scale=1.0/av_fac_day
        else
           scale=0.0
        endif
        def2%name=trim(uEMEP%poll(ipoll)%emis)
        call Out_netCDF(iotyp,def2,ndim_tot,kmax,loc_tot_day(1,1,KMAX_MID-uEMEPNvertout+1,ipoll),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
             fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       

     else  if(iotyp==IOU_MON)then

        do iisec=1,uEMEP%poll(ipoll)%Nsectors
           isec_poll=isec_poll1+iisec-1
           isec=uEMEP%poll(ipoll)%sector(iisec)
           !copy before dividing by loc_tot
           if(COMPUTE_LOCAL_TRANSPORT)call extendarea_N(loc_frac_month(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),tmp_ext,uEMEP%dist,uEMEP_Sizedxdy,uEMEPNvertout)
           if(iisec==1)fracsum=0.0
           do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                    invtot=1.0/(loc_tot_month(i,j,k,ipoll)+1.E-20)
                    do dy=-uEMEP%dist,uEMEP%dist
                       do dx=-uEMEP%dist,uEMEP%dist
                          loc_frac_month(dx,dy,i,j,k,isec_poll)=loc_frac_month(dx,dy,i,j,k,isec_poll)*invtot
                          if(iisec==1.and.k==KMAX_MID)fracsum(i,j)=fracsum(i,j)+loc_frac_month(dx,dy,i,j,k,isec_poll)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           scale=1.0
           write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_fraction'
           if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_fraction'

           call Out_netCDF(iotyp,def1,ndim,kmax,loc_frac_month(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)   

           if(iisec==1)then
              def1%name=trim(uEMEP%poll(ipoll)%emis)//'_fracsum'
              call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)  
           endif

           if(isec==0 .and. COMPUTE_LOCAL_TRANSPORT)then
              do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
                 do j=1,ljmax
                    do i=1,limax
                       do dy=-uEMEP%dist,uEMEP%dist
                          do dx=-uEMEP%dist,uEMEP%dist
                             !in  loc_poll_to (dx,dy,i,j) shows from (i,j) to (i+dx,j+dy)
                             loc_poll_to(dx,dy,i,j,k)=tmp_ext(-dx,-dy,i+dx,j+dy,k)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_transport'
              if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_transport'
              if(abs(av_fac_month)>1.E-5)then
                 scale=1.0/av_fac_month
              else
                 scale=0.0
              endif
              call Out_netCDF(iotyp,def1,ndim,kmax,loc_poll_to,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
           endif
        enddo

       if(abs(av_fac_month)>1.E-5)then
           scale=1.0/av_fac_month
        else
           scale=0.0
        endif
        def2%name=trim(uEMEP%poll(ipoll)%emis)
        call Out_netCDF(iotyp,def2,ndim_tot,kmax,loc_tot_month(1,1,KMAX_MID-uEMEPNvertout+1,ipoll),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
             fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       

     else  if(iotyp==IOU_YEAR)then

        do iisec=1,uEMEP%poll(ipoll)%Nsectors
           isec_poll=isec_poll1+iisec-1
           isec=uEMEP%poll(ipoll)%sector(iisec)
           !copy before dividing by loc_tot_full
           if(COMPUTE_LOCAL_TRANSPORT)then
              t1 = MPI_WTIME()
              call extendarea_N(loc_frac_full(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),tmp_ext,uEMEP%dist,uEMEP_Sizedxdy,uEMEPNvertout)
              t2 = MPI_WTIME()
              if(me==0)write(*,*)'transport extendarea ',t2-t1,' seconds'
           endif
           if(iisec==1)fracsum=0.0
           do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                    invtot=1.0/(loc_tot_full(i,j,k,ipoll)+1.E-20)
                    do dy=-uEMEP%dist,uEMEP%dist
                       do dx=-uEMEP%dist,uEMEP%dist
                          loc_frac_full(dx,dy,i,j,k,isec_poll)=loc_frac_full(dx,dy,i,j,k,isec_poll)*invtot
                         if(iisec==1.and.k==KMAX_MID)fracsum(i,j)=fracsum(i,j)+loc_frac_full(dx,dy,i,j,k,isec_poll)
                       enddo
                    enddo
                 enddo
              enddo
           enddo

           scale=1.0
           write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_fraction'
           if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_fraction'
           call Out_netCDF(iotyp,def1,ndim,kmax,loc_frac_full(-uEMEP%dist,-uEMEP%dist,1,1,KMAX_MID-uEMEPNvertout+1,isec_poll),scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       

           if(iisec==1)then
              def1%name=trim(uEMEP%poll(ipoll)%emis)//'_fracsum'
              call Out_netCDF(iotyp,def1,ndim_tot,1,fracsum,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)  
           endif

           if(isec==0 .and. COMPUTE_LOCAL_TRANSPORT)then
              t1 = MPI_WTIME()
              do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
                 !     do k = KMAX_MID,KMAX_MID
                 do j=1,ljmax
                    do i=1,limax
                       !invtot=1.0/(1.E-20+loc_tot_full(i,j,k))
                       do dy=-uEMEP%dist,uEMEP%dist
                          do dx=-uEMEP%dist,uEMEP%dist
                             !in  loc_poll_to (dx,dy,i,j) shows from (i,j) to (i+dx,j+dy)
                             loc_poll_to(dx,dy,i,j,k)=tmp_ext(-dx,-dy,i+dx,j+dy,k)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              t2 = MPI_WTIME()
              if(me==0)write(*,*)'transport transpose ',t2-t1,' seconds'
              write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_sec',isec,'_local_transport'
              if(isec==0)write(def1%name,"(A,I2.2,A)")trim(uEMEP%poll(ipoll)%emis)//'_local_transport'
              if(abs(av_fac_full)>1.E-5)then
                 scale=1.0/av_fac_full
              else
                 scale=0.0
              endif
             call Out_netCDF(iotyp,def1,ndim,kmax,loc_poll_to,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
                   fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
             t1 = MPI_WTIME()
             if(me==0)write(*,*)'transport out ',t1-t2,' seconds'
           endif
         enddo

        if(abs(av_fac_full)>1.E-5)then
           scale=1.0/av_fac_full
        else
           scale=0.0
        endif
        def2%name=trim(uEMEP%poll(ipoll)%emis)
        call Out_netCDF(iotyp,def2,ndim_tot,kmax,loc_tot_full(1,1,KMAX_MID-uEMEPNvertout+1,ipoll),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
             fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)   
     else
        if(me==0)write(*,*)'IOU not recognized'
     endif
     isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
  enddo

  !reset the cumulative counters
  if(iotyp==IOU_HOUR)then
     av_fac_hour=0
     loc_frac_hour=0.0
     loc_tot_hour=0.0
  else  if(iotyp==IOU_DAY)then
     av_fac_day=0.0
     loc_frac_day=0.0
     loc_tot_day=0.0
  else  if(iotyp==IOU_MON)then
     av_fac_month=0.0
     loc_frac_month=0.0
     loc_tot_month=0.0
  else  if(iotyp==IOU_YEAR)then
     av_fac_full=0.0
     loc_frac_full=0.0
     loc_tot_full=0.0
  endif

  first_call(iotyp)=.false.
  if(COMPUTE_LOCAL_TRANSPORT)deallocate(tmp_ext)

  call Add_2timing(NTIMING-2,tim_after,tim_before,"uEMEP: output")

! CALL MPI_BARRIER(MPI_COMM_CALC, I)
 
!stop
end subroutine out_uEMEP

subroutine av_uEMEP(dt,End_of_Day)
  real, intent(in)    :: dt                   ! time-step used in integrations
  logical, intent(in) :: End_of_Day           ! e.g. 6am for EMEP sites
  real :: xtot
  integer ::i,j,k,dx,dy,ix,iix,ipoll,isec_poll1
  integer ::isec_poll
  
  call Code_timer(tim_before)
  if(.not. uEMEP%HOUR.and.&
     .not. uEMEP%DAY .and.&
     .not. uEMEP%MONTH .and.&
     .not. uEMEP%YEAR       )return

  !do the averaging
  do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
     do j=1,ljmax
        do i=1,limax
           isec_poll1=1
           do ipoll=1,uEMEP%Npoll
              xtot=0.0
              do iix=1,uEMEP%poll(ipoll)%Nix
                 ix=uEMEP%poll(ipoll)%ix(iix)
                 xtot=xtot+(xn_adv(ix,i,j,k)*uEMEP%poll(ipoll)%mw(iix))/ATWAIR&
                      *roa(i,j,k,1)*1.E9 !for ug/m3
                 !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
              end do
              if(uEMEP%HOUR)then
                 loc_tot_hour(i,j,k,ipoll)=loc_tot_hour(i,j,k,ipoll)+xtot
                 do dy=-uEMEP%dist,uEMEP%dist
                    do dx=-uEMEP%dist,uEMEP%dist
                       do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                          loc_frac_hour(dx,dy,i,j,k,isec_poll)=loc_frac_hour(dx,dy,i,j,k,isec_poll)+xtot*loc_frac(isec_poll,dx,dy,i,j,k)
                       enddo
                    enddo
                 enddo
              endif
              if(uEMEP%DAY)then
                 loc_tot_day(i,j,k,ipoll)=loc_tot_day(i,j,k,ipoll)+xtot
                 do dy=-uEMEP%dist,uEMEP%dist
                    do dx=-uEMEP%dist,uEMEP%dist
                       do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                          loc_frac_day(dx,dy,i,j,k,isec_poll)=loc_frac_day(dx,dy,i,j,k,isec_poll)+xtot*loc_frac(isec_poll,dx,dy,i,j,k)
                       enddo
                    enddo
                 enddo
              endif
              if(uEMEP%MONTH)then
                 loc_tot_month(i,j,k,ipoll)=loc_tot_month(i,j,k,ipoll)+xtot
                 do dy=-uEMEP%dist,uEMEP%dist
                    do dx=-uEMEP%dist,uEMEP%dist
                       do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                          loc_frac_month(dx,dy,i,j,k,isec_poll)=loc_frac_month(dx,dy,i,j,k,isec_poll)+xtot*loc_frac(isec_poll,dx,dy,i,j,k)
                       enddo
                    enddo
                 enddo
              endif
              if(uEMEP%YEAR)then
                 loc_tot_full(i,j,k,ipoll)=loc_tot_full(i,j,k,ipoll)+xtot
                 do dy=-uEMEP%dist,uEMEP%dist
                    do dx=-uEMEP%dist,uEMEP%dist
                       do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                          loc_frac_full(dx,dy,i,j,k,isec_poll)=loc_frac_full(dx,dy,i,j,k,isec_poll)+xtot*loc_frac(isec_poll,dx,dy,i,j,k)
                       enddo
                    enddo
                 enddo
              endif
              isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
           enddo
        enddo
     enddo
  enddo
  av_fac_hour=av_fac_hour+1.0
  av_fac_day=av_fac_day+1.0
  av_fac_month=av_fac_month+1.0
  av_fac_full=av_fac_full+1.0

  call Add_2timing(NTIMING-8,tim_after,tim_before,"uEMEP: averaging")

end subroutine av_uEMEP

  subroutine uemep_adv_y(fluxy,i,j,k)
    real, intent(in)::fluxy(NSPEC_ADV,-1:LJMAX+1)
    integer, intent(in)::i,j,k
    real ::x,xn,xx,f_in,inv_tot
    integer ::iix,ix,dx,dy,ipoll,isec_poll,isec_poll1

    call Code_timer(tim_before)
    isec_poll1=1
    do ipoll=1,uEMEP%Npoll
       xn=0.0
       x=0.0
       xx=0.0
       !positive x or xx means incoming, negative means outgoing
       do iix=1,uEMEP%poll(ipoll)%Nix
          ix=uEMEP%poll(ipoll)%ix(iix)
          xn=xn+xn_adv(ix,i,j,k)*uEMEP%poll(ipoll)%mw(iix)
          x=x-xm2(i,j)*fluxy(ix,j)*uEMEP%poll(ipoll)%mw(iix)!flux through "North" face (Up)
          xx=xx+xm2(i,j)*fluxy(ix,j-1)*uEMEP%poll(ipoll)%mw(iix)!flux through "South" face (Bottom)
       end do
       !NB: here xn already includes the fluxes. Remove them!
       xn=xn-xx-x
       
       xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux 
       f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
       inv_tot=1.0/(xn+f_in+1.e-20)!incoming dilutes
       
       xx=max(0.0,xx)*inv_tot!factor due to flux through "South" face (Bottom)
       x =max(0.0,x)*inv_tot!factor due to flux through "North" face (Up)

       do dy=-uEMEP%dist,uEMEP%dist
          do dx=-uEMEP%dist,uEMEP%dist
             do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k) *xn *inv_tot
             enddo
             if(x>0.0.and.dy>-uEMEP%dist)then
                do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                   loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)+ loc_frac_1d(isec_poll,dx,dy-1,j+1)*x
                enddo
             endif
             if(xx>0.0.and.dy<uEMEP%dist)then
                do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                   loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)+ loc_frac_1d(isec_poll,dx,dy+1,j-1)*xx
                enddo
             endif
          enddo
       enddo
       isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
    enddo
    call Add_2timing(NTIMING-6,tim_after,tim_before,"uEMEP: adv_y")
  end subroutine uemep_adv_y

  subroutine uemep_adv_x(fluxx,i,j,k)
    real, intent(in)::fluxx(NSPEC_ADV,-1:LIMAX+1)
    integer, intent(in)::i,j,k
    real ::x,xn,xx,f_in,inv_tot
    integer ::iix,ix,dx,dy,ipoll,isec_poll,isec_poll1

    call Code_timer(tim_before)
    isec_poll1=1
    do ipoll=1,uEMEP%Npoll
       xn=0.0
       x=0.0
       xx=0.0
       !positive x or xx means incoming, negative means outgoing
       do iix=1,uEMEP%poll(ipoll)%Nix
          ix=uEMEP%poll(ipoll)%ix(iix)
          xn=xn+xn_adv(ix,i,j,k)*uEMEP%poll(ipoll)%mw(iix)
          x=x-xm2(i,j)*fluxx(ix,i)*uEMEP%poll(ipoll)%mw(iix)!flux through "East" face (Right)
          xx=xx+xm2(i,j)*fluxx(ix,i-1)*uEMEP%poll(ipoll)%mw(iix)!flux through "West" face (Left)
       end do
       !NB: here xn already includes the fluxes. Remove them!
       xn=xn-xx-x

       xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux 
       f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
       inv_tot=1.0/(xn+f_in+1.e-20)!incoming dilutes
       
       x =max(0.0,x)*inv_tot!factor due to flux through "East" face (Right)
       xx=max(0.0,xx)*inv_tot!factor due to flux through "West" face (Left)
    
       do dy=-uEMEP%dist,uEMEP%dist
          do dx=-uEMEP%dist,uEMEP%dist
             do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k) *xn *inv_tot
             enddo
             if(x>0.0.and.dx>-uEMEP%dist)then
                do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                   loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)+ loc_frac_1d(isec_poll,dx-1,dy,i+1)*x
                enddo
             endif
             if(xx>0.0.and.dx<uEMEP%dist)then
                do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                   loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)+ loc_frac_1d(isec_poll,dx+1,dy,i-1)*xx
                enddo
             endif
          enddo
       enddo
       isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
    enddo
    call Add_2timing(NTIMING-7,tim_after,tim_before,"uEMEP: adv_x")
 
  end subroutine uemep_adv_x

  subroutine uemep_adv_k(fluxk,i,j)
    real, intent(in)::fluxk(NSPEC_ADV,KMAX_MID)
    integer, intent(in)::i,j
    real ::x,xn,xx,f_in,inv_tot
    integer ::k,iix,ix,dx,dy,ipoll,isec_poll,isec_poll1
    real loc_frac_km1(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,KMAX_MID-uEMEP%Nvert:KMAX_MID-1)

    call Code_timer(tim_before)
    !need to be careful to always use non-updated values on the RHS
    do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID-1
       loc_frac_km1(:,:,:,k)=loc_frac(:,:,:,i,j,k)
    enddo
    loc_frac_km1(:,:,:,KMAX_MID-uEMEP%Nvert)=0.0!Assume zero local fractions coming from above

    do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID!k is increasing-> can use k+1 to access non-updated value

       isec_poll1=1
       do ipoll=1,uEMEP%Npoll

          xn=0.0
          x=0.0
          xx=0.0
          do iix=1,uEMEP%poll(ipoll)%Nix
             ix=uEMEP%poll(ipoll)%ix(iix)
             xn=xn+xn_adv(ix,i,j,k)*uEMEP%poll(ipoll)%mw(iix)
             if(k<KMAX_MID)x=x-dhs1i(k+1)*fluxk(ix,k+1)*uEMEP%poll(ipoll)%mw(iix)
             xx=xx+dhs1i(k+1)*fluxk(ix,k)*uEMEP%poll(ipoll)%mw(iix)
          end do
          !NB: here xn already includes the fluxes. Remove them!
          xn=xn-xx-x
          
          xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. outgoing flux 
          f_in=max(0.0,x)+max(0.0,xx)!positive part. incoming flux
          inv_tot = 1.0/(xn+f_in+1.e-20)
          
          x =max(0.0,x)*inv_tot!factor due to flux through bottom face
          xx=max(0.0,xx)*inv_tot!factor due to flux through top face
          if(k<KMAX_MID)then
             if(k>KMAX_MID-uEMEP%Nvert+1)then
                do dy=-uEMEP%dist,uEMEP%dist
                   do dx=-uEMEP%dist,uEMEP%dist
                      do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                         loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k) * xn * inv_tot &
                              + loc_frac_km1(isec_poll,dx,dy,k-1) * xx&
                              + loc_frac(isec_poll,dx,dy,i,j,k+1) * x!k is increasing-> can use k+1 to access non-updated value
                      enddo
                   enddo
                enddo
             else
                !k=KMAX_MID-uEMEP%Nvert+1 , assume no local fractions from above
                do dy=-uEMEP%dist,uEMEP%dist
                   do dx=-uEMEP%dist,uEMEP%dist
                      do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                         loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)* xn * inv_tot &
                              + loc_frac(isec_poll,dx,dy,i,j,k+1) * x!k is increasing-> can use k+1 to access non-updated value
                      enddo
                   enddo
                enddo
             endif
          else
             !k=KMAX_MID , no local fractions from below
             do dy=-uEMEP%dist,uEMEP%dist
                do dx=-uEMEP%dist,uEMEP%dist
                   do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                      loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)* xn * inv_tot &
                           + loc_frac_km1(isec_poll,dx,dy,k-1) * xx
                   enddo
                enddo
             enddo
          endif
          isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
       enddo
       
    end do
    call Add_2timing(NTIMING-5,tim_after,tim_before,"uEMEP: adv_k")
  end subroutine uemep_adv_k
 
  subroutine uemep_diff(i,j,ds3,ds4,ndiff)
    
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
!    real :: xn_k(uEMEP_Size1+uEMEP%Nsec_poll,KMAX_MID),x
    real :: xn_k(uEMEP_Size1+uEMEP%Npoll,KMAX_MID),x
    integer ::isec_poll1,ipoll,isec_poll
    integer ::k,n,ix,iix,dx,dy
    !how far diffusion should take place above uEMEP%Nvert. 
    ! KUP = 2 gives less than 0.001 differences in locfrac, except sometimes over sea, because
    !ship emission are higher up and need to come down to diminish locfrac
    integer, parameter :: KUP = 2

    call Code_timer(tim_before)
    xn_k = 0.0
    do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
       isec_poll1=1
       do ipoll=1,uEMEP%Npoll
          x=0.0
          do iix=1,uEMEP%poll(ipoll)%Nix
             ix=uEMEP%poll(ipoll)%ix(iix)
             !assumes mixing ratios units, but weight by mass
             x=x+xn_adv(ix,i,j,k)*uEMEP%poll(ipoll)%mw(iix)
          end do
          n=0
          do dy=-uEMEP%dist,uEMEP%dist
             do dx=-uEMEP%dist,uEMEP%dist
                n=n+1
                do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                   xn_k((n-1)*(uEMEP%Nsec_poll)+isec_poll,k)=x*loc_frac(isec_poll,dx,dy,i,j,k)
                end do
             end do
          end do
          isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
       enddo
    enddo
    do k = 1,KMAX_MID
       isec_poll1=1
       do ipoll=1,uEMEP%Npoll
          x=0.0
          do iix=1,uEMEP%poll(ipoll)%Nix
             ix=uEMEP%poll(ipoll)%ix(iix)
             !assumes mixing ratios units, but weight by mass
             x=x+xn_adv(ix,i,j,k)*uEMEP%poll(ipoll)%mw(iix)
          end do
          xn_k(uEMEP_Size1+ipoll,k)=x
          isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
       enddo
    enddo
    
    call vertdiffn(xn_k,uEMEP_Size1+uEMEP%Npoll,1,KMAX_MID-uEMEP%Nvert-KUP,EtaKz(i,j,1,1),ds3,ds4,ndiff)
    
    do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
       isec_poll1=1
       do ipoll=1,uEMEP%Npoll
          n=0
          do dy=-uEMEP%dist,uEMEP%dist
             do dx=-uEMEP%dist,uEMEP%dist
                n=n+1
                x = 1.0/(xn_k(uEMEP_Size1+ipoll,k)+1.E-30)
                do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                   loc_frac(isec_poll,dx,dy,i,j,k) = xn_k((n-1)*(uEMEP%Nsec_poll)+isec_poll,k)*x
                end do
             end do
          end do
          isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
       end do
    end do
    call Add_2timing(NTIMING-4,tim_after,tim_before,"uEMEP: diffusion")
  end subroutine uemep_diff

subroutine uEMEP_emis(indate)
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
  real ::  tfac    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s       ! source term (emis) before splitting
  integer :: iland, iland_timefac, iland_timefac_hour  ! country codes, and codes for timefac 
  integer :: hour_iland
  integer ::icc_uemep, iqrc, itot
  integer, save :: wday , wday_loc ! wday = day of the week 1-7
  integer ::ix,iix, dx, dy, isec_poll, iisec_poll, isec_poll1, ipoll
  real::dt_uemep, xtot, emis_uemep(KMAX_MID,NEMIS_FILE,NSECTORS),emis_tot(KMAX_MID,NEMIS_FILE)
  real :: lon

  call Code_timer(tim_before)
  dt_uemep=dt_advec

  wday=day_of_week(indate%year,indate%month,indate%day)
  if(wday==0)wday=7 ! Sunday -> 7
  do iland = 1, NLAND
    daytime(iland) = 0
    hourloc        = indate%hour + Country(iland)%timezone
    localhour(iland) = hourloc  ! here from 0 to 23
    if(hourloc>=7 .and. hourloc<=18) daytime(iland)=1
  end do ! iland
   
  do j = lj0,lj1
    do i = li0,li1
      ncc = nlandcode(i,j)            ! No. of countries in grid

      !*************************************************
      ! loop over sector emissions
      !*************************************************
      tmpemis(:)=0.
      icc_uemep=0
      emis_uemep=0.0
      emis_tot=0.0
      do icc = 1, ncc
        !          iland = landcode(i,j,icc)     ! 1=Albania, etc.
        iland=find_index(landcode(i,j,icc),Country(:)%icode) !array index
        call make_iland_for_time(.false., indate, i, j, iland, wday, iland_timefac,hour_iland,wday_loc,iland_timefac_hour)
        
        do iem = 1, NEMIS_FILE 
           do isec = 1, Nsectors     ! Loop over snap codes
            ! Calculate emission rates from secemis, time-factors, 
            ! and if appropriate any speciation fraction (NEMIS_FRAC)
            ! kg/m2/s
            
            tfac = timefac(iland_timefac,sec2tfac_map(isec),iem) &
                   * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)

            if(USES%GRIDDED_EMIS_MONTHLY_FACTOR)tfac=tfac* GridTfac(i,j,sec2tfac_map(isec),iem)
            
            !Degree days - only SNAP-2 
            if(USES%DEGREEDAY_FACTORS .and. &
                 sec2tfac_map(isec)==ISNAP_DOM .and. Gridded_SNAP2_Factors) then
               ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
               ! we make use of a baseload even for SNAP2
               tfac = ( fac_min(iland,sec2tfac_map(isec),iem) & ! constant baseload
                    + ( 1.0-fac_min(iland,sec2tfac_map(isec),iem) )* gridfac_HDD(i,j) ) &
                    * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)
            end if ! =============== HDD 

            s = tfac * secemis(isec,i,j,icc,iem)

            do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
              emis_tot(k,iem)=emis_tot(k,iem)+s*emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))*dt_uemep
            end do

            !if(isec==uEMEP%sector .or. uEMEP%sector==0)then
            do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
               emis_uemep(k,iem,isec)=emis_uemep(k,iem,isec)+s*emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))*dt_uemep
            end do
            !end if
            
          end do ! iem

        end do  ! isec
        !      ==================================================
      end do ! icc  
                     
      !Add emissions from new format
      do n = 1, NEmis_sources      
         if(Emis_source(n)%include_in_local_fractions)then            
            itot = Emis_source(n)%species_ix
            isec = Emis_source(n)%sector
            iland = Emis_source(n)%country_ix
            if(itot>0)then
               !the species is directly defined (no splits)
               iqrc = itot2iqrc(itot)
               if(isec>0)then
                  call CheckStop(itot2iqrc(itot)<=0,"emitted sector species must belong to one of the splitted species")
               endif
               iem = iqrc2iem(iqrc)
            else
               iem=find_index(Emis_source(n)%species,EMIS_FILE(:))
            endif
            if(isec>0)then
               if(Emis_source(n)%periodicity == 'yearly' .or. Emis_source(n)%periodicity == 'monthly')then
                  !we need to apply hourly factors
                  call make_iland_for_time(.false., indate, i, j, iland, wday, iland_timefac,hour_iland,wday_loc,iland_timefac_hour)                      
                  tfac = fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)
                  if(Emis_source(n)%periodicity == 'yearly')then
                     !apply monthly factor on top of hourly factors
                     tfac = tfac * timefac(iland_timefac,sec2tfac_map(isec),iem)                         
                  endif
               else
                  !not monthly or yearly emissions, timefactors must be included in emission values
                  tfac = 1.0
               endif

               if (USES%GRIDDED_EMIS_MONTHLY_FACTOR) tfac=tfac* GridTfac(i,j,sec2tfac_map(isec),iem)
               
               !Degree days - only SNAP-2 
               if(USES%DEGREEDAY_FACTORS .and. &
                    sec2tfac_map(isec)==ISNAP_DOM .and. Gridded_SNAP2_Factors) then
                  ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
                  ! we make use of a baseload even for SNAP2
                  tfac = ( fac_min(iland,sec2tfac_map(isec),iem) & ! constant baseload
                       + ( 1.0-fac_min(iland,sec2tfac_map(isec),iem) )* gridfac_HDD(i,j) ) &
                       * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc,iland_timefac_hour)
               end if ! =============== HDD 
 
               s = Emis_source_2D(i,j,n) * tfac
               
               do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
                  emis_tot(k,iem)=emis_tot(k,iem)+s*emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))*dt_uemep
               end do
               
               do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
                  emis_uemep(k,iem,isec)=emis_uemep(k,iem,isec)+s*emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))*dt_uemep
               end do
            endif
         endif
      enddo

      isec_poll1=1
      do ipoll=1,uEMEP%Npoll              
         iem = uEMEP%poll(ipoll)%EMIS_File_ix
         
         do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
            if(emis_tot(k,iem)<1.E-20)cycle
            !units kg/m2
            !total pollutant
            xtot=0.0
            do iix=1,uEMEP%poll(ipoll)%Nix
               ix=uEMEP%poll(ipoll)%ix(iix)
               xtot=xtot+(xn_adv(ix,i,j,k)*uEMEP%poll(ipoll)%mw(iix))*(dA(k)+dB(k)*ps(i,j,1))/ATWAIR/GRAV
            end do
            dx=0 ; dy=0!local fraction from this i,j
            do iisec_poll=1,uEMEP%poll(ipoll)%Nsectors
               isec_poll=iisec_poll+isec_poll1-1
               isec = uEMEP%poll(ipoll)%sector(iisec_poll)
               if(isec==0)then
                  !sum of all sectors
                  loc_frac(isec_poll,dx,dy,i,j,k)=(loc_frac(isec_poll,dx,dy,i,j,k)*xtot+emis_tot(k,iem))/(xtot+emis_tot(k,iem)+1.e-20)
              else
                  loc_frac(isec_poll,dx,dy,i,j,k)=(loc_frac(isec_poll,dx,dy,i,j,k)*xtot+emis_uemep(k,iem,isec))/(xtot+emis_tot(k,iem)+1.e-20)
               endif
           enddo

            !local fractions from other cells are updated (reduced)            
            do dy=-uEMEP%dist,uEMEP%dist
               do dx=-uEMEP%dist,uEMEP%dist
                  if(dx==0 .and. dy==0)cycle!local fractions from other cells only
                  do isec_poll=isec_poll1,isec_poll1+uEMEP%poll(ipoll)%Nsectors-1
                     loc_frac(isec_poll,dx,dy,i,j,k)=(loc_frac(isec_poll,dx,dy,i,j,k)*xtot)/(xtot+emis_tot(k,iem)+1.e-20)
                  enddo
               enddo
            enddo
            
         end do! k
         isec_poll1=isec_poll1+uEMEP%poll(ipoll)%Nsectors
      enddo!Npoll
    end do ! i
  end do ! j

  call Add_2timing(NTIMING-3,tim_after,tim_before,"uEMEP: emissions")

end subroutine uEMEP_emis


end module uEMEP_mod
