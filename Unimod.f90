! <Unimod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2018 met.no
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
program myeul
  !-----------------------------------------------------------------------!
  !
  !     This is the main program for the off-line regional scale multilayer
  !     eulerian model at emep/msc-w. the main program contains the outer
  !     time-loop which runs through all time-levels a new meteorological
  !     data-set is read into the model from file. the inner time-loop
  !     runs through the physical time-step.
  !
  !-----------------------------------------------------------------------!

  use My_Outputs_ml,    only: set_output_defs
  use My_Timing_ml,     only: lastptim, mytimm, Output_timing, &
       Init_timing, Add_2timing, Code_timer, &
       tim_before, tim_before1, tim_before2, &
       tim_after, tim_after0, NTIMING_UNIMOD,NTIMING
  use Advection_ml,     only: vgrid_Eta, assign_nmax, assign_dtadvec
  use Aqueous_ml,       only: init_aqueous, Init_WetDep   !  Initialises & tabulates
  use AirEmis_ml,       only: lightning
  use BiDir_emep,       only : Init_BiDir  !  FUTURE
  use Biogenics_ml,     only: Init_BVOC, SetDailyBVOC
  use BoundaryConditions_ml, only: BoundaryConditions
  use CheckStop_ml,     only: CheckStop
  use Chemfields_ml,    only: alloc_ChemFields
  use ChemSpecs,        only: define_chemicals
  use ChemGroups_ml,    only: Init_ChemGroups
  use Country_ml,       only: init_Country
  use DA_3DVar_ml,      only: NTIMING_3DVAR,DA_3DVar_Init, DA_3DVar_Done
  use DefPhotolysis_ml, only: readdiss
  use Derived_ml,       only: Init_Derived, wanted_iou
  use EcoSystem_ml,     only: Init_EcoSystems
  use Emissions_ml,     only: Emissions, newmonth
  use ForestFire_ml,    only: Fire_Emis
  use GridValues_ml,    only: MIN_ADVGRIDS, GRIDWIDTH_M, Poles,&
                              DefDebugProc, GridRead
  use Io_ml,            only: IO_MYTIM,IO_RES,IO_LOG,IO_NML,IO_DO3SE
  use Io_Progs_ml,      only: read_line, PrintLog
  use Landuse_ml,       only: InitLandUse, SetLanduse
  use MassBudget_ml,    only: Init_massbudget, massbudget
  use Met_ml,           only: metfieldint, MetModel_LandUse, Meteoread
  use Config_module,only: MasterProc, &   ! set true for host processor, me==MasterPE
       RUNDOMAIN,  &   ! Model domain
       NPROC,      &   ! No. processors
       METSTEP,    &   ! Hours between met input
       runlabel1,  &   ! explanatory text
       runlabel2,  &   ! explanatory text
       iyr_trend, nmax,nstep , meteo,     &
       IOU_INST,IOU_HOUR,IOU_HOUR_INST, IOU_YEAR,IOU_MON, IOU_DAY, &
       USES, USE_uEMEP,JUMPOVER29FEB,&
       FORECAST,ANALYSIS  ! FORECAST/ANALYSIS mode
  use Config_module,only: Config_ModelConstants,DEBUG, startdate,enddate
  use MPI_Groups_ml,    only: MPI_BYTE, MPISTATUS, MPI_COMM_CALC,MPI_COMM_WORLD, &
                              MasterPE,IERROR, MPI_world_init
  use Nest_ml,          only: wrtxn     ! write nested output (IC/BC)
  use NetCDF_ml,        only: Init_new_netCDF
  use OutputChem_ml,    only: WrtChem, wanted_iou
  use Par_ml,           only: me, GIMAX, GJMAX, Topology_io, Topology, parinit
  use PhyChem_ml,       only: phyche    ! Calls phys/chem routines each dt_advec
  use Sites_ml,         only: sitesdef  ! to get output sites
  use SmallUtils_ml,    only: key2str
  use Tabulations_ml,   only: tabulate
  use TimeDate_ml,      only: date, current_date, day_of_year, daynumber,&
       tdif_secs,date,timestamp,make_timestamp,Init_nmdays
  use TimeDate_ExtraUtil_ml,only : date2string, assign_startandenddate
  use Trajectory_ml,    only: trajectory_init,trajectory_in
  use uEMEP_ml,         only: init_uEMEP
  !--------------------------------------------------------------------
  !
  !  Variables. There are too many to list here. Still, here are a
  !  few key variables that  might help:
  !     dt_advec       - length of advection (phyche) time-step
  !     GRIDWIDTH_M    - grid-distance
  !     gb             - latitude (sorry, still Norwegian influenced..)
  !     NPROC          - number of processors used
  !     me             - number of local processor, me=0 is host (=MasterProc)
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
       "COMPILED WRONGLY: Need double precision, e.g. f90 -r8")

  if(MasterProc) open(IO_LOG,file='RunLog.out')
!TEST
  call define_chemicals()    ! sets up species details
  call Config_ModelConstants(IO_LOG)

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
    call CheckStop(status,"DA_3DVar_Init in Unimod")
  end if

  !*** Timing ********
  call Init_timing(NTIMING_UNIMOD+NTIMING_3DVAR)
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

  ! daynumber needed  for BCs, so call here to get initial
  daynumber=day_of_year(yyyy,mm,dd)

  !-------------------------------------------------------------------
  !
  !++  parameters and initial fields.
  !

  call alloc_ChemFields     !allocate chemistry arrays
!TEST  call define_chemicals()    ! sets up species details
  call Init_ChemGroups()    ! sets up species details

  call assign_nmax(METSTEP)   ! No. timesteps in inner loop

  call trajectory_init()

  call init_Country() ! In Country_ml, => NLAND, country codes and names, timezone

  call Add_2timing(1,tim_after,tim_before,"Grid init + chem init")

  call SetLandUse(daynumber, mm) !  Reads Inputs.Landuse, Inputs.LandPhen

  call Add_2timing(1,tim_after,tim_before,"landuse read in")

  call MeteoRead()

  call Add_2timing(2,tim_after,tim_before,"Meteo read first record")

  if (MasterProc.and.DEBUG%MAINCODE) print *,"Calling emissions with year",yyyy

  call Emissions(yyyy)

  call Add_2timing(3,tim_after,tim_before,"Yearly emissions read in")

  if(USE_uEMEP) call init_uEMEP

  call MetModel_LandUse(1)   !

  call Init_EcoSystems()     ! Defines ecosystem-groups for dep output

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
    call Init_new_netCDF(trim(runlabel1)//'_hour.nc',IOU_HOUR)
  if(wanted_iou(IOU_HOUR_INST)) &
    call Init_new_netCDF(trim(runlabel1)//'_hourInst.nc',IOU_HOUR_INST)

  call Add_2timing(4,tim_after,tim_before,"Other init")

  tim_before = tim_before0
  call Add_2timing(5,tim_after,tim_before,"Total until time loop")
  call Code_timer(tim_before1)

  mm_old = 0
  oldseason = 0
  nstep=0

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !     performance of physical and chemical calculations,
  !     three-hourly time loop starts here
  !
  do while(.not.End_of_Run)        ! main time-loop , timestep is dt_advec

    nstep=mod(nstep,nmax)+1 !loops from 1 to nmax, between two meteo read

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
      call Code_timer(tim_before)

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

    if(USES%FOREST_FIRES) call Fire_Emis(daynumber)

    call Add_2timing(10,tim_after,tim_before,"Fires+BVOC")

    if(MasterProc) print "(2(1X,A))",'current date and time:',&
      date2string("YYYY-MM-DD hh:mm:ss",current_date)

    call Code_timer(tim_before2)
    call phyche()
    call Add_2timing(11,tim_after,tim_before2,"Total phyche")

    call Code_timer(tim_before)

    call WrtChem()

    call Add_2timing(36,tim_after,tim_before,"WrtChem")

    call trajectory_in

    call metfieldint

    call Add_2timing(37,tim_after,tim_before,"metfieldint")

    !this is a bit complicated because it must account for the fact that for instance 3feb24:00 = 4feb00:00 
    ts1=make_timestamp(current_date)
    ts2=make_timestamp(date(enddate(1),enddate(2),enddate(3),enddate(4),0))
    End_of_Run =  (nint(tdif_secs(ts1,ts2))<=0)

    if(DEBUG%STOP_HH>=0 .and. DEBUG%STOP_HH==current_date%hour) &
      End_of_Run=.true.

  end do ! time-loop

  call Code_timer(tim_after0)
  call Add_2timing(38,tim_after0,tim_before1,"total within loops")
  call Add_2timing(39,tim_after0,tim_before0,"total")
  !
  !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(.not.FORECAST) call wrtxn(current_date,.true.)
  call massbudget()

  if(MasterProc)then
    print *,'programme is finished'
    ! Gather timing info:
    if(NPROC-1> 0)then
      CALL MPI_RECV(lastptim,NTIMING*8,MPI_BYTE,NPROC-1,765,MPI_COMM_CALC,MPISTATUS,IERROR)
    else
      lastptim(:) = mytimm(:)
    end if
    call Output_timing(IO_MYTIM,me,NPROC,GIMAX,GJMAX)
  elseif(me==NPROC-1) then
    CALL MPI_SEND(mytimm,NTIMING*8,MPI_BYTE,MasterPE,765,MPI_COMM_CALC,IERROR)
  end if

  ! write 'modelrun.finished' file to flag the end of the FORECAST
  if(MasterProc.and.FORECAST)then
    open(1,file='modelrun.finished')
    close(1)
  end if

  if(ANALYSIS)then              ! assimilation enabled
    call DA_3DVar_Done(status)  ! done with 3D-var module:
    call CheckStop(status,"DA_3DVar_Done in Unimod")
  end if

  CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
  CALL MPI_FINALIZE(IERROR)

endassociate   ! yyyy, mm, dd
endprogram

!===========================================================================
