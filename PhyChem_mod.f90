! <PhyChem_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.6>
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
module PhyChem_mod
!
!     physical and chemical routine calls within one advection step
!     driven from here
!
!     Output of hourly data
!
!-----------------------------------------------------------------------------
use Advection_mod,     only: advecdiff_Eta
use Biogenics_mod,     only: Set_ACP2012EuroSoilNOx
use CheckStop_mod,     only: CheckStop
use Chemfields_mod,    only: xn_adv,cfac,xn_shl
use ChemDims_mod,      only: NSPEC_SHL
use ChemSpecs_mod,     only: species
use CoDep_mod,         only: make_so2nh3_24hr
use Config_module,only: MasterProc, KMAX_MID, nmax, step_main,END_OF_EMEPDAY &
                           ,dt_advec       & ! time-step for phyche/advection
                           ,PPBINV, PPTINV  &
                           ,IOU_INST       &
                           ,ANALYSIS       & ! 3D-VAR Analysis
                           ,SOURCE_RECEPTOR&
                           ,USES&
                           ,FREQ_HOURLY    & ! hourly netcdf output frequency
                           ,JUMPOVER29FEB&
                           ,IOU_HOUR, IOU_HOUR_INST, IOU_YEAR&
                           ,fileName_O3_Top,FREQ_SITE, FREQ_SONDE&
                           , O3_ix,SO2_ix,NH3_ix
use DA_mod,            only: DEBUG_DA_1STEP
use DA_3DVar_mod,      only: main_3dvar, T_3DVAR
use Debug_module,      only: DEBUG, DebugCell  ! -> DEBUG%GRIDVALUES
use Derived_mod,       only: DerivedProds, Derived, num_deriv2d
use DerivedFields_mod, only: d_2d, f_2d
use DryDep_mod,        only: init_drydep
use EmisDef_mod,       only: NSECTORS
use Emissions_mod,     only: EmisSet
use Gravset_mod,       only: gravset
!EXP use GravSettling_mod,  only: gravSettling ! EXPERIMENTAl!!
use GridValues_mod,    only: debug_proc,debug_li,debug_lj,&
                            glon,glat,projection,i_local,j_local,i_fdom,j_fdom
use Io_Progs_mod,       only: datewrite

use MetFields_mod,     only: ps,roa,z_bnd,z_mid,cc3dmax, &
                            PARdbh, PARdif, fCloud, & !WN17, PAR in W/m2
                            zen,coszen
use My_Timing_mod,     only: NTIMING, Code_timer, Add_2timing, &
                             tim_before, tim_before0, tim_after
use NetCDF_mod,        only: ReadField_CDF,Real4,ReadTimeCDF
use Nest_mod,          only: readxn, wrtxn
use OutputChem_mod,    only: WrtChem
use OwnDataTypes_mod, only: TXTLEN_FILE
use Par_mod,           only: me, LIMAX, LJMAX
use PhysicalConstants_mod, only : ATWAIR
use Radiation_mod,     only: SolarSetup,       &! sets up radn params
                            ZenithAngle,      &! gets zenith angle
                            WeissNormanPAR,   &! WN17=WeissNorman radn
                            fCloudAtten,      &!
                            CloudAtten         !
use Runchem_mod,       only: runchem   ! Calls setup subs and runs chemistry
use Sites_mod,         only: siteswrt_surf, siteswrt_sondes    ! outputs
  use SmallUtils_mod,  only: key2str
use SoilWater_mod,     only: Set_SoilWater
use TimeDate_mod,      only: date,daynumber,day_of_year, add_secs, &
                            current_date, timestamp,  &
                            make_timestamp, make_current_date
use TimeDate_ExtraUtil_mod,only : date2string, date2nctime, nctime2date, &
                                  date_is_reached
use Timefactors_mod,   only: NewDayFactors
use Trajectory_mod,    only: trajectory_out     ! 'Aircraft'-type  outputs

!-----------------------------------------------------------------------------
implicit none
private

public  :: phyche
private :: debug_concs

contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine phyche()

  logical, save :: End_of_Day = .false.
  integer :: ndays,status,kstart
  real :: thour
  type(timestamp) :: ts_now !date in timestamp format
  logical,save :: first_call = .true.
  character(len=*), parameter :: dtxt='phyche:'
  integer, save :: O3record = -1 !initialize to show not set
  real :: time_wanted, TimesInDays(1)
  integer,save  :: O3_end_date(5) = (/0,0,0,0,0/)
  character(len=TXTLEN_FILE), save :: fileName_O3_Top_current


  !------------------------------------------------------------------
  !     physical and  chemical routines.

  !     Hours since midnight at any time-step
  !    using current_date we have already step_main taken into account
  thour = real(current_date%hour) + current_date%seconds/3600.0

  if ( DEBUG%STOFLUX .and. debug_proc .and. first_call ) then
   write(*,'(a,3i4,2f10.4)')dtxt//"SEIij ", me,debug_li,debug_lj,& 
     glon(debug_li,debug_lj),glat(debug_li,debug_lj)
  end if

  if(DEBUG%PHYCHEM.and.debug_proc ) then
    if(first_call)  write(*, *) "PhyChe First", me, debug_proc

    call debug_concs("PhyChe start ")

    if(current_date%hour==12) then
      ndays = day_of_year(current_date%year,current_date%month, &
           current_date%day)
      write(*,*) 'thour,ndays,step_main,dt', thour,ndays,step_main,dt_advec
    end if
  end if

  if(trim(fileName_O3_Top)/="NOTSET" .and.&
       date_is_reached(O3_end_date))then
     !read in new O3 top data
     if(first_call .or. &
        (current_date%hour<=3.and.current_date%day==1.and.current_date%month==1))then
        !we redefine the file to be read and where to start
        !should account for:
        !   a) the first record (time 0:0 first of January) may be missing
        !   b) we may start a new year after 31 Dec midnight 
        fileName_O3_Top_current = key2str(fileName_O3_Top,'YYYY',current_date%year)
        if(MasterProc)then
           write(*,*)dtxt//'Reading 3 hourly O3 at top from :'
           write(*,*)trim(fileName_O3_Top_current)
        end if
        !read dates from files to find where to start
        call date2nctime(current_date,time_wanted) !converts into real, number of days since 1900
        O3record = 1 ! on input, means find one date only
        !will find first date after or equal to current date
        call ReadTimeCDF(trim(fileName_O3_Top_current),TimesInDays,O3record,time_wanted)
        call nctime2date(O3_end_date, TimesInDays(1) + 3.0/24.0/2) !assumes 3-hourly values, instantaneous
     else
        O3record = O3record + 1 !next record
        O3_end_date(4) = O3_end_date(4) + 3 !We assume 3-hourly values (ok if > 24!)
     end if
     
     if(DEBUG%PHYCHEM .and. MasterProc) then
       write(*,*)'Updating top O3 with record ',&
          O3record,", file ",trim(fileName_O3_Top_current)
       write(*,*)'valid until time ',O3_end_date(1),O3_end_date(2),O3_end_date(3),O3_end_date(4),O3_end_date(5)
     end if

     kstart=3!NB: must be the level corresponding to model top! Hardcoded for now
     call  ReadField_CDF(trim(fileName_O3_Top_current),'O3',xn_adv(O3_ix-NSPEC_SHL,:,:,1),&
          nstart=O3record,kstart=kstart,kend=kstart,&
          interpol='zero_order',debug_flag=.false.)
   
  endif

  call Code_timer(tim_before)
  call readxn(current_date) !Read xn_adv from earlier runs

  call Add_2timing(15,tim_after,tim_before,"nest: Read")

  ! analysis enabled?
  if(ANALYSIS.and.first_call)then
    ! 3D-VAR Analysis for "Zero hour"
    call main_3dvar(status)
    call CheckStop(status,"main_3dvar in PhyChem_mod/PhyChe")
    ! update timing:
    call Add_2timing(T_3DVAR,tim_after,tim_before)
    ! testing ...
    if(DEBUG_DA_1STEP)then
      ! info ..
      if ( MasterProc ) then
        write(*,*) 'ANALYSIS DEBUG_DA_1STEP: only 1st assimilation step'
      end if
      ! update derived output:
      call Derived(dt_advec,End_of_Day)
      ! leave:
      return
    end if ! debug da 1step
  end if  ! analysis at zero hour

  ! For safety we initialise instant. values here to zero.
  ! Usually not needed, but sometimes
  ! ========================
  d_2d(:,:,:,IOU_INST) = 0.0
  ! ========================


  !===================================

  call SolarSetup(current_date%year,current_date%month,current_date%day,thour)

  call ZenithAngle(thour, glat, glon, zen, coszen )

  if(DEBUG%PHYCHEM.and.debug_proc)&
    write(*,"(a,2i4,3f7.3,9es12.3)") "PhyChem ZenRad ", current_date%day, current_date%hour, &
        thour, glon(debug_li,debug_lj),glat(debug_li,debug_lj), &
        zen(debug_li,debug_lj),coszen(debug_li,debug_lj)

!WN17, Gets PAR values, W/m2 here
  fCloud = fCloudAtten(cc3dmax(:,:,KMAX_MID))
  call WeissNormanPAR(ps(:,:,1),coszen,fCloud,PARdbh,PARdif)
  if ( DEBUG%STOFLUX .and. debug_proc ) then
   call datewrite(dtxt//"SEI ", [me], [ &
     thour, zen(debug_li,debug_lj),coszen(debug_li,debug_lj),&
     cc3dmax(debug_li,debug_lj,KMAX_MID), fCloud(debug_li,debug_lj), &
     PARdbh(debug_li,debug_lj)], afmt='(TXTDATE,a,i4,2f9.3,4es12.3)')
  end if

  !================
  ! advecdiff_poles considers the local Courant number along a 1D line
  ! and divides the advection step "locally" in a number of substeps.
  ! Up north in a LatLong domain such as MACC02, mapfactors go up to four,
  ! so using advecdiff_poles pays off, even though none of the poles are
  ! included in the domain.
  ! For efficient parallellisation each subdomain needs to have the same work
  ! load; this can be obtained by setting NPROCY=1 (number of subdomains in
  ! latitude- or y-direction).
  ! Then, all subdomains have exactly the same geometry.

  call Code_timer(tim_before0)

  call advecdiff_Eta

  call Add_2timing(13,tim_after,tim_before0,"phyche: total advecdiff")

  if(USES%ASH) call gravset

  !================

  call Code_timer(tim_before)

  !/ See if we are calculating any before-after chemistry productions:

  !=============================
  if (mod(step_main,nmax) == 0) call DerivedProds("Before",dt_advec)
  !=============================


  !===================================
  call Set_SoilWater()
  !===================================
  if(USES%SOILNOX_METHOD=="ACP2012EURO") &
     call Set_ACP2012EuroSoilNOx()!hourly, deprecated

  !===================================
  call init_drydep()
  !===================================

  call EmisSet(current_date)

  call Add_2timing(12,tim_after,tim_before,"phyche:EmisSet")

  call Code_timer(tim_before0)

  !=========================================================!
  call debug_concs("PhyChe pre-chem ")

  !************ NOW THE HEAVY BIT **************************!

  call Code_timer(tim_before0)
  call runchem()   !  calls setup subs and runs chemistry

  call Add_2timing(23,tim_after,tim_before0,"Total Runchem")
  call debug_concs("PhyChe post-chem ")

  !*********************************************************!
  !EXP if(USES%GRAVSET) call gravSettling ! After chem/emis?? EXPERIMENTAL!
  !========================================================!


  !/ See if we are calculating any before-after chemistry productions:

  !=============================
  if(mod(step_main,nmax) == 0) call DerivedProds("After",dt_advec)
  !=============================

  !=============================
  ! this output needs the 'old' current_date_hour

  call trajectory_out
  !=============================

  !    the following partly relates to end of time step - hourly output
  !    partly not depends on current_date
  !    => add dt_advec to current_date already here


  !====================================
  ts_now = make_timestamp(current_date)

  call add_secs(ts_now,dt_advec)
  current_date = make_current_date(ts_now)

  if(JUMPOVER29FEB.and.current_date%month==2.and.current_date%day==29)then
    if(MasterProc)write(*,*)'Jumping over one day for current_date!'
    if(MasterProc)print "(2(1X,A))",'current date and time before jump:',&
          date2string("YYYY-MM-DD hh:mm:ss",current_date)
    call add_secs(ts_now,24*3600.)
    current_date = make_current_date(ts_now)
    if(MasterProc)print "(2(1X,A))",'current date and time after jump:',&
        date2string("YYYY-MM-DD hh:mm:ss",current_date)
  end if

  call Code_timer(tim_before)
  !====================================
  ! data-assimilation enabled?
  if(ANALYSIS)then
    ! 3D-VAR Analysis for "non-Zero hours"
    call main_3dvar( status )
    call CheckStop(status,"main_3dvar in PhyChem_mod/PhyChe")
    ! end timing:
    call Add_2timing(T_3DVAR,tim_after,tim_before)
  end if
  call wrtxn(current_date,.false.) !Write xn_adv for future nesting
  call Add_2timing(14,tim_after,tim_before,"nest: Write")

  End_of_Day=(current_date%seconds==0).and.(current_date%hour==END_OF_EMEPDAY)

  if(End_of_Day.and.MasterProc)then
    print "(a,a)",' End of EMEP-day ',date2string("(hh:mm:ss)",current_date)
    if(DEBUG%PHYCHEM)write(*,"(a20,2i4,i6)") "END_OF_EMEPDAY ", &
      END_OF_EMEPDAY, current_date%hour,current_date%seconds
  end if

  call debug_concs("PhyChe pre-Derived ")

  call Derived(dt_advec,End_of_Day)

  call Add_2timing(34,tim_after,tim_before,"phyche:Derived")

  ! Hourly Outputs:

  if(current_date%seconds==0) then

    if(.not.SOURCE_RECEPTOR .and. FREQ_SITE>0 .and.&
      modulo(current_date%hour,FREQ_SITE)==0)  &
      call siteswrt_surf(xn_adv,cfac,xn_shl)

    if(.not.SOURCE_RECEPTOR .and. FREQ_SONDE>0 .and. &
      modulo(current_date%hour,FREQ_SONDE)==0) &
      call siteswrt_sondes(xn_adv,xn_shl)

    call Add_2timing(35,tim_after,tim_before,"phyche:sites and hourly out")

  end if

  ! CoDep
  if(modulo(current_date%hour,1)==0 .and. NH3_ix>0 .and. SO2_ix>0) & ! every hour
    call make_so2nh3_24hr(current_date%hour,&
      xn_adv(SO2_ix-NSPEC_SHL,:,:,KMAX_MID),&
      xn_adv(NH3_ix-NSPEC_SHL,:,:,KMAX_MID),&
      cfac(SO2_ix-NSPEC_SHL,:,:),&
      cfac(NH3_ix-NSPEC_SHL,:,:))

  first_call=.false.
end subroutine phyche
!--------------------------------------------------------------------------
subroutine debug_concs(txt)
  character(len=*), intent(in) :: txt
  character(len=6), save :: unit
  real :: c1, c2
  integer, save :: ispec, iadv
  logical :: first_call = .true.

  ! Simple sub to print out eg O3 concentrations for different stages
  ! of calculation. Saves lots of messy lines above.

  if(DEBUG%PHYCHEM.and.debug_proc)then
    if(first_call) then
      ispec = DEBUG%ISPEC
      iadv  = DEBUG%ISPEC - NSPEC_SHL
      unit = 'ppbv'
      if(ispec<=NSPEC_SHL) unit='pptv'
      first_call = .false.
    end if

    if(ispec>NSPEC_SHL)then
      c1=xn_adv(iadv,debug_li,debug_lj,KMAX_MID)*PPBINV
      c2=c1* cfac(iadv,debug_li,debug_lj)
    else
      c1=xn_shl(ispec,debug_li,debug_lj,KMAX_MID)*PPTINV
      c2=-1.0
    end if
    write(*,"(a,2i3,i5,i3,a12,2g12.4,1x,a4)") "debug_concs:"// &
      trim(txt), me, current_date%hour, current_date%seconds, step_main,&
      trim(species(ispec)%name), c1, c2, unit
  end if
end subroutine debug_concs
!--------------------------------------------------------------------------
end module PhyChem_mod
