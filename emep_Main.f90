! <emep_Main.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
program emep_Main
  !-----------------------------------------------------------------------!
  !
  !     This is the main program for the off-line regional scale multilayer
  !     eulerian model at emep/msc-w. the main program contains the outer
  !     time-loop which runs through all time-levels a new meteorological
  !     data-set is read into the model from file. the inner time-loop
  !     runs through the physical time-step.
  !
  !-----------------------------------------------------------------------!

  use My_Timing_mod,     only: lastptim, mytimm, Output_timing, &
       Init_timing, Add_2timing, Code_timer, &
       tim_before, tim_before1, tim_before2, &
       tim_after, tim_after0, NTIMING_UNIMOD,NTIMING
  use Advection_mod,     only: vgrid_Eta, assign_nmax, assign_dtadvec
  use Aqueous_mod,       only: init_aqueous, Init_WetDep   !  Initialises & tabulates
  use AirEmis_mod,       only: lightning
  use BiDir_emep,        only : Init_BiDir  !  FUTURE
  use Biogenics_mod,     only: Init_BVOC, SetDailyBVOC
  use BoundaryConditions_mod, only: BoundaryConditions
  use CheckStop_mod,     only: CheckStop
  use Chemfields_mod,    only: alloc_ChemFields
  use ChemSpecs_mod,     only: define_chemicals
  use ChemGroups_mod,    only: Init_ChemGroups
  use Config_module,only: MasterProc, &   ! true for host processor, me==0
       RUNDOMAIN,  &   ! Model domain
       NPROC,      &   ! No. processors
       METSTEP,    &   ! Hours between met input
       runlabel1,  &   ! explanatory text
       runlabel2,  &   ! explanatory text
       iyr_trend, nmax,step_main , meteo,     &
       IOU_INST,IOU_HOUR,IOU_HOUR_INST, IOU_YEAR,IOU_MON, IOU_DAY, &
       HOURLYFILE_ending, &
       USES, USE_uEMEP,JUMPOVER29FEB,&
       ANALYSIS, & ! forecast in ANALYSIS mode
       Config_Constants, startdate, enddate
  use Country_mod,       only: init_Country
  use DA_3DVar_mod,      only: NTIMING_3DVAR,DA_3DVar_Init, DA_3DVar_Done
  use Debug_module,      only: DEBUG   ! -> DEBUG%MAINCODE
  use DefPhotolysis_mod, only: readdiss
  use Derived_mod,       only: Init_Derived, wanted_iou
  use EcoSystem_mod,     only: Init_EcoSystems
  use Emissions_mod,     only: Emissions, newmonth, Init_masks, Init_emissions,&
                               EmisUpdate
  use ForestFire_mod,    only: Fire_Emis
  use DryDep_mod,        only: init_DryDep ! sets up dry and wet dep
  !use GasParticleCoeffs_mod, only: init_DryDep ! sets up dry and wet dep
  use GridValues_mod,    only: MIN_ADVGRIDS, GRIDWIDTH_M, Poles,&
                              DefDebugProc, GridRead, set_EuropeanAndGlobal_Config
  use Io_mod,            only: IO_MYTIM,IO_RES,IO_LOG,IO_NML,IO_DO3SE
  use Io_Progs_mod,      only: read_line, PrintLog
  use Landuse_mod,       only: InitLandUse, SetLanduse
  use MassBudget_mod,    only: Init_massbudget, massbudget
  use Met_mod,           only: metfieldint, MetModel_LandUse, Meteoread
  use MPI_Groups_mod!,    only: MPI_BYTE, MPISTATUS, MPI_COMM_CALC,MPI_COMM_WORLD, &
                    !          MasterPE,IERROR, MPI_world_init
  use Nest_mod,          only: wrtxn     ! write nested output (IC/BC)
  use NetCDF_mod,        only: Init_new_netCDF
  use OutputChem_mod,    only: WrtChem, wanted_iou, set_output_defs 
  use Par_mod,           only: me, GIMAX, GJMAX, Topology_io, Topology, parinit
  use PhyChem_mod,       only: phyche    ! Calls phys/chem routines each dt_advec
  use Sites_mod,         only: sitesdef  ! to get output sites
  use SmallUtils_mod,    only: key2str
  use Tabulations_mod,   only: tabulate
  use TimeDate_mod,      only: date, current_date, day_of_year, daynumber,&
       tdif_secs,date,timestamp,make_timestamp,Init_nmdays
  use TimeDate_ExtraUtil_mod,only : date2string, assign_startandenddate,&
                                    date_is_reached
  use Trajectory_mod,    only: trajectory_init,trajectory_in
  use uEMEP_mod,         only: init_uEMEP, NTIMING_uEMEP
  !--------------------------------------------------------------------
  !
  !  Variables. There are too many to list here. Still, here are a
  !  few key variables that  might help:
  !     dt_advec       - length of advection and time splitting time-step
  !     GRIDWIDTH_M    - approximate grid-size (multiply by mapfactor for exact)
  !     gb             - latitude (sorry, still Norwegian influenced..)
  !     NPROC          - number of processors used
  !     me             - rank of local processor, me=0 is host (=MasterProc)
  !                      processor where many read/writes are done
  !     ndays          - number of days since 1 january (max 365 or 366)
  !     thour          - utc-time in hours every time-step
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !     +---------------------------------------------------------+
  !     +                                                         +
  !     +                                                         +
  !     +                  start main programme                   +
  !     +                                                         +
  !     +_________________________________________________________+

  !     declarations in main programme.
  !
  implicit none

  integer :: i, oldseason, newseason, status
  integer :: mm_old   ! month and old-month
  integer :: cyclicgrid
  TYPE(timestamp)   :: ts1,ts2
  logical :: End_of_Run=.false.
  real :: tim_before0 !private
  character(len=*), parameter :: dtxt='eMain:'

  associate ( yyyy => current_date%year, mm => current_date%month, &
       dd => current_date%day,  hh => current_date%hour)
    !
    !     initialize the parallel topology
    !
    
  call MPI_world_init(NPROC,ME)

  ! Set a logical from ModelConstants, which can be used for
  !   specifying the master processor for print-outs and such
  MasterProc = ( me == MasterPE )

  call CheckStop(digits(1.0)<50, &
       dtxt//"COMPILED WRONGLY: Need double precision, e.g. f90 -r8")

  if(MasterProc) open(IO_LOG,file='RunLog.out')
!TEST
  call define_chemicals()    ! sets up species details
  call Config_Constants(IO_LOG)

  call assign_startandenddate()
 
  if(MasterProc)then
     call PrintLog(trim(runlabel1))
     call PrintLog(trim(runlabel2))
     call PrintLog(date2string("startdate = YYYYMMDDhh",startdate(1:4)))
     call PrintLog(date2string("enddate   = YYYYMMDDhh",enddate  (1:4)))
    !call PrintLog(key2str("iyr_trend = YYYY","YYYY",iyr_trend))
  end if


  if(ANALYSIS)then              ! init 3D-var module
    call DA_3DVar_Init(status)  ! pass settings
    call CheckStop(status,dtxt//"DA_3DVar_Init in emepctm")
  end if

  !*** Timing ********
  call Init_timing(NTIMING_UNIMOD+NTIMING_3DVAR+NTIMING_uEMEP)
  call Code_Timer(tim_before0)
  tim_before = tim_before0

  call GridRead(meteo,cyclicgrid) ! define:
  ! 1) grid sizes (IIFULLDOM, JJFULLDOM,KMAX_MID),
  ! 2) projection (lon lat or Stereographic etc and Poles),
  ! 3) rundomain size (GIMAX, GJMAX, IRUNBEG, JRUNBEG)
  ! 4) vertical levels defintion and interpolation coefficients
  ! 5) subdomain partition (NPROCX, NPROCY, limax,ljmax)
  ! 6) topology (neighbor, poles)
  ! 7) grid properties arrays (xm, i_local, j_local etc.)

  call Topology(cyclicgrid,Poles)   ! def GlobalBoundaries & subdomain neighbors
  call DefDebugProc()               ! Sets debug_proc, debug_li, debuglj
  call assign_dtadvec(GRIDWIDTH_M)  ! set dt_advec

  call set_EuropeanAndGlobal_Config() !Set config values that depend on domain coverage

  ! daynumber needed  for BCs, so call here to get initial
  daynumber=day_of_year(yyyy,mm,dd)

  !-------------------------------------------------------------------
  !
  !++  parameters and initial fields.
  !

  call alloc_ChemFields     !allocate chemistry arrays
!TEST  call define_chemicals()    ! sets up species details
  call Init_ChemGroups()      ! sets up species details

  call trajectory_init()

  call init_Country() ! In Country_mod, => NLAND, country codes and names, timezone

  call Add_2timing(1,tim_after,tim_before,"Grid init + chem init")

  call SetLandUse(daynumber, mm) !  Reads Inputs.Landuse, Inputs.LandPhen

  call Add_2timing(1,tim_after,tim_before,"landuse read in")

  call MeteoRead()

  call assign_nmax(METSTEP)   ! No. timesteps in inner loop

  call Add_2timing(2,tim_after,tim_before,"Meteo read first record")

  if (MasterProc.and.DEBUG%MAINCODE) print *,"Calling emissions with year",yyyy

  call Init_masks()
  call Emissions(yyyy)! should be set for the enddate year, not start?
  call Init_emissions !new format

  call Add_2timing(3,tim_after,tim_before,"Yearly emissions read in")

  if(USE_uEMEP) call init_uEMEP

  call MetModel_LandUse(1)   !

  call Init_EcoSystems()     ! Defines ecosystem-groups for dep output

  call init_DryDep()        ! sets up dry and wet dep arrays

  call Init_Derived()        ! Derived field defs.

  call Init_BVOC()

  call Init_BiDir()           ! BIDIR FUTURE 

  call tabulate()             ! sets up tab_esat, etc.

  call Init_WetDep()           ! sets up scavenging ratios

  call set_output_defs()     ! Initialises outputs
  call sitesdef()            ! see if any output for specific sites is wanted
  ! (read input files "sites.dat" and "sondes.dat" )

  call vgrid_Eta           !  initialisation of constants used in vertical advection
  if (MasterProc.and.DEBUG%MAINCODE ) print *,"vgrid finish"

  ! open only output netCDF files if needed
  if(MasterProc.and.DEBUG%MAINCODE )&
    print *, "NETCDFINITS: iou", (i,wanted_iou(i),i=IOU_INST,IOU_HOUR_INST)
  ! The fullrun file contains the accumulated or average results
  ! over the full run period, often a year, but even just for
  ! a few timesteps if that is all that is run:

  if(wanted_iou(IOU_INST)) &
    call Init_new_netCDF(trim(runlabel1)//'_inst.nc',IOU_INST)
  if(wanted_iou(IOU_YEAR)) &
    call Init_new_netCDF(trim(runlabel1)//'_fullrun.nc',IOU_YEAR)
  if(wanted_iou(IOU_MON)) &
    call Init_new_netCDF(trim(runlabel1)//'_month.nc',IOU_MON)
  if(wanted_iou(IOU_DAY)) &
    call Init_new_netCDF(trim(runlabel1)//'_day.nc',IOU_DAY)
  if(wanted_iou(IOU_HOUR)) &
    call Init_new_netCDF(trim(runlabel1)//'_hour'//&
    date2string(HOURLYFILE_ending,current_date),IOU_HOUR)
  if(wanted_iou(IOU_HOUR_INST)) &
    call Init_new_netCDF(trim(runlabel1)//'_hourInst'//&
    date2string(HOURLYFILE_ending,current_date),IOU_HOUR_INST)

  call Add_2timing(4,tim_after,tim_before,"Other init")

  tim_before = tim_before0
  call Add_2timing(5,tim_after,tim_before,"Total until time loop")
  call Code_timer(tim_before1)

  mm_old = 0
  oldseason = 0
  step_main=0

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !     performance of physical and chemical calculations,
  !     three-hourly time loop starts here
  !
  do while(.not.End_of_Run)        ! main time-loop , timestep is dt_advec

    step_main=step_main + 1 !main loop with dt_advec between each step

    !FUTURE if (NH3EMIS_VAR) call SetNH3()  ! NH3emis experimental

    select case (mm)
      case(12,1:2);newseason = 1
      case(3:5)   ;newseason = 2
      case(6:8)   ;newseason = 3
      case(9:11)  ;newseason = 4
    end select

    ! daynumber needed for BCs
    daynumber=day_of_year(yyyy,mm,dd)
     
    if(mm==1 .and. dd==1 .and. hh==0)call Init_nmdays(current_date, JUMPOVER29FEB)!new year starts

    call Code_timer(tim_before)
    if(mm_old/=mm) then   ! START OF NEW MONTH !!!!!

      !subroutines/data that must be updated every month
      call readdiss(newseason)

      if(MasterProc.and.DEBUG%MAINCODE) &
        print *,'maaned og sesong', mm,mm_old,newseason,oldseason

      call MetModel_LandUse(2)   ! e.g.  gets snow_flag
      if(MasterProc.and.DEBUG%MAINCODE) write(*,*)"Newmonth start"

      call newmonth

       if(USES%LIGHTNING_EMIS) call lightning()

      call init_aqueous()

      ! Monthly call to BoundaryConditions.
      if(DEBUG%MAINCODE) print *, "Into BCs" , me
      ! We set BCs using the specified iyr_trend
      !   which may or may not equal the meteorology year
      call Code_timer(tim_before2)
      call BoundaryConditions(yyyy,mm)
      call Add_2timing(6,tim_after,tim_before2,"BoundaryConditions")

      if(DEBUG%MAINCODE) print *, "Finished BCs" , me

      !must be called only once, after BC is set
      if(mm_old==0)call Init_massbudget()
      if(DEBUG%MAINCODE) print *, "Finished Initmass" , me

      call Add_2timing(7,tim_after,tim_before,"Total newmonth setup")

    end if

    oldseason = newseason
    mm_old = mm

    if(DEBUG%MAINCODE) print *, "1st Infield" , me

    call SetLandUse(daynumber, mm) !daily
    call Add_2timing(8,tim_after,tim_before,"SetLanduse")

    call Meteoread() ! 3-hourly or hourly

    call Add_2timing(9,tim_after,tim_before,"Meteoread")

    call SetDailyBVOC() !daily

    call EmisUpdate

    if(USES%FOREST_FIRES) call Fire_Emis(daynumber)

    call Add_2timing(10,tim_after,tim_before,"Fires+BVOC")

    call Code_timer(tim_before2)
    call phyche()
    call Add_2timing(11,tim_after,tim_before2,"Total phyche")

    if(MasterProc) print "(2(1X,A))",'current date and time:',&
      date2string("YYYY-MM-DD hh:mm:ss",current_date)

    call Code_timer(tim_before)

    call WrtChem()

    call Add_2timing(37,tim_after,tim_before,"WrtChem")

    call trajectory_in

    call metfieldint

    call Add_2timing(36,tim_after,tim_before,"metfieldint")

    End_of_Run = date_is_reached(enddate)

    if(DEBUG%STOP_HH>=0 .and. DEBUG%STOP_HH==current_date%hour) &
      End_of_Run=.true.

  end do ! time-loop

  call Code_timer(tim_after0)
  call Add_2timing(NTIMING-1,tim_after0,tim_before1,"total within loops")
  call Add_2timing(NTIMING,tim_after0,tim_before0,"total")
  !
  !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  call wrtxn(current_date,.true.)
  call massbudget()

  call Output_timing(IO_MYTIM,me,NPROC,GIMAX,GJMAX)

  ! write 'modelrun.finished' file to flag the end of the run
  if(MasterProc.and.USES%PROGRESS_FILES)then
    open(1,file='modelrun.finished')
    close(1)
  end if

  if(ANALYSIS)then              ! assimilation enabled
    call DA_3DVar_Done(status)  ! done with 3D-var module:
    call CheckStop(status,"DA_3DVar_Done in emepctm")
  end if

  CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
  if(MasterProc)print *,'programme is finished'
  CALL MPI_FINALIZE(IERROR)

end associate   ! yyyy, mm, dd
end program emep_Main

!===========================================================================
