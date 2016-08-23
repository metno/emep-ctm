! <Unimod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_5(2809)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-201409 met.no
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

  use My_Outputs_ml,    only: set_output_defs, NHOURLY_OUT
  use My_Timing_ml,     only: lastptim, mytimm, Output_timing, &
       Init_timing, Add_2timing, Code_timer, &
       tim_before, tim_before0, tim_before1, &
       tim_after, tim_after0
  use Advection_ml,     only: vgrid, assign_nmax, assign_dtadvec
  use Aqueous_ml,       only: init_aqueous, Init_WetDep   !  Initialises & tabulates
  use AirEmis_ml,       only: lightning
  use Biogenics_ml,     only: Init_BVOC, SetDailyBVOC
  use BoundaryConditions_ml, only: BoundaryConditions
  use CheckStop_ml,     only: CheckStop
  use Chemfields_ml,    only: alloc_ChemFields
!CMR  use ChemChemicals_ml, only: define_chemicals
  use ChemSpecs,        only: define_chemicals
  use ChemGroups_ml,    only: Init_ChemGroups
  use DefPhotolysis_ml, only: readdiss
  use Derived_ml,       only: Init_Derived, iou_min, iou_max
  use DerivedFields_ml, only: f_2d, f_3d
  use EcoSystem_ml,     only: Init_EcoSystems
  use Emissions_ml,     only: Emissions, newmonth
  use ForestFire_ml,    only: Fire_Emis
  use GridValues_ml,    only: MIN_ADVGRIDS, GRIDWIDTH_M, Poles, DefDebugProc, GridRead
  use Io_ml,            only: IO_MYTIM,IO_RES,IO_LOG,IO_NML,IO_DO3SE
  use Io_Progs_ml,      only: read_line, PrintLog
  use Landuse_ml,       only: InitLandUse, SetLanduse, Land_codes
  use MassBudget_ml,    only: Init_massbudget, massbudget
  use Met_ml,           only: metfieldint, MetModel_LandUse, Meteoread, meteo
  use ModelConstants_ml,only: MasterProc, &   ! set true for host processor, me==0
       RUNDOMAIN,  &   ! Model domain
       NPROC,      &   ! No. processors
       METSTEP,    &   ! Hours between met input
       runlabel1,  &   ! explanatory text
       runlabel2,  &   ! explanatory text
       nprint,nterm,iyr_trend,    nmax,nstep ,                  &
       IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY, &
       USES, USE_LIGHTNING_EMIS, &
       FORECAST       ! FORECAST mode
  use ModelConstants_ml,only: Config_ModelConstants,DEBUG
  use NetCDF_ml,        only: Init_new_netCDF
  use OutputChem_ml,    only: WrtChem, wanted_iou
  use Par_ml,           only: me, GIMAX, GJMAX, Topology, parinit
  use PhyChem_ml,       only: phyche    ! Calls phys/chem routines each dt_advec
  use Sites_ml,         only: sitesdef  ! to get output sites
  use Tabulations_ml,   only: tabulate
  use TimeDate_ml,      only: date, current_date, day_of_year, daynumber,&
       tdif_secs,date,timestamp,make_timestamp,startdate, enddate,Init_nmdays
  use TimeDate_ExtraUtil_ml,only : date2string, assign_NTERM
  use Trajectory_ml,    only: trajectory_init,trajectory_in
  use Nest_ml,          only: wrtxn     ! write nested output (IC/BC)
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

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO

  integer :: i, oldseason, newseason
  integer :: mm_old   ! month and old-month
  integer :: nproc_mpi,cyclicgrid
  character (len=230) :: errmsg,txt
  TYPE(timestamp)   :: ts1,ts2
  logical :: End_of_Run=.false.

  namelist /INPUT_PARA/iyr_trend,runlabel1,runlabel2,&
       startdate,enddate,meteo

  associate ( yyyy => current_date%year, mm => current_date%month, &
       dd => current_date%day,  hh => current_date%hour)
    !
    !     initialize the parallel topology
    !
  nproc_mpi = NPROC
  CALL MPI_INIT(INFO)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, INFO)
  ! Set a logical from ModelConstants, which can be used for
  !   specifying the master processor for print-outs and such
  MasterProc = ( me == 0 )
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc_mpi, INFO)
  NPROC=nproc_mpi 
55 format(A,I5,A)
  if(MasterProc)write(*,55)' Found ',NPROC,' MPI processes available'

  call CheckStop(digits(1.0)<50, &
       "COMPILED WRONGLY: Need double precision, e.g. f90 -r8")

  if(MasterProc) open(IO_LOG,file='RunLog.out')
!TEST
  call define_chemicals()    ! sets up species details
  call Config_ModelConstants(IO_LOG)

  rewind(IO_NML)
  read(IO_NML,NML=INPUT_PARA)
  startdate(4)=0                ! meteo hour to start/end the run 
  enddate  (4)=0                ! are set in assign_NTERM

  if(MasterProc)then
     call PrintLog(trim(runlabel1))
     call PrintLog(trim(runlabel2))
     call PrintLog(date2string("startdate = YYYYMMDD",startdate(1:3)))
     call PrintLog(date2string("enddate   = YYYYMMDD",enddate  (1:3)))
     ! write(txt,"(a,i4)") "iyr_trend= ", iyr_trend
     ! call PrintLog(trim(txt))
  endif

  !*** Timing ********
  call Init_timing()
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
  call assign_NTERM(NTERM)          ! set NTERM, the number of 3-hourly periods
  call assign_dtadvec(GRIDWIDTH_M)  ! set dt_advec

  ! daynumber needed  for BCs, so call here to get initial
  daynumber=day_of_year(yyyy,mm,dd)

  !     Decide the frequency of print-out
  !
  nprint = nterm

  if (MasterProc) print *,'nterm, nprint:',nterm, nprint

  !-------------------------------------------------------------------
  !
  !++  parameters and initial fields.
  !
  call Add_2timing(1,tim_after,tim_before,"Before define_Chemicals")

  call alloc_ChemFields     !allocate chemistry arrays
!TEST  call define_chemicals()    ! sets up species details
  call Init_ChemGroups()    ! sets up species details

  call assign_nmax(METSTEP)   ! No. timesteps in inner loop

  call trajectory_init()

  call Add_2timing(2,tim_after,tim_before,"After define_Chems, readpar")

  call SetLandUse(daynumber, mm) !  Reads Inputs.Landuse, Inputs.LandPhen

  call MeteoRead()

  call Add_2timing(3,tim_after,tim_before,"After infield")

  if (MasterProc.and.DEBUG%MAINCODE) print *,"Calling emissions with year",yyyy

  call Emissions(yyyy)


  call MetModel_LandUse(1)   !

  call Init_EcoSystems()     ! Defines ecosystem-groups for dep output

  call Init_Derived()        ! Derived field defs.

  call Init_BVOC()

  call tabulate()             ! sets up tab_esat, etc.

  call Init_WetDep()           ! sets up scavenging ratios

  call set_output_defs()     ! Initialises outputs
  call sitesdef()            ! see if any output for specific sites is wanted
  ! (read input files "sites.dat" and "sondes.dat" )

  call vgrid           !  initialisation of constants used in vertical advection
  if (MasterProc.and.DEBUG%MAINCODE ) print *,"vgrid finish"

  ! open only output netCDF files if needed
  if(MasterProc.and.DEBUG%MAINCODE )print *, "NETCDFINITS: minval, maxval", iou_min, iou_max
  ! The fullrun file contains the accumulated or average results
  ! over the full run period, often a year, but even just for
  ! a few timesteps if that is all that is run:

  if (wanted_iou(IOU_YEAR)) &
       call Init_new_netCDF(trim(runlabel1)//'_fullrun.nc',IOU_YEAR)
  if (wanted_iou(IOU_INST)) &
       call Init_new_netCDF(trim(runlabel1)//'_inst.nc',IOU_INST)
  ! if (wanted_iou(IOU_HOUR).or.NHOURLY_OUT>0) &
  !   call Init_new_netCDF(trim(runlabel1)//'_hour.nc',IOU_HOUR)
  if (wanted_iou(IOU_DAY)) &
       call Init_new_netCDF(trim(runlabel1)//'_day.nc',IOU_DAY)
  if (wanted_iou(IOU_MON)) &
       call Init_new_netCDF(trim(runlabel1)//'_month.nc',IOU_MON)

  call Add_2timing(4,tim_after,tim_before,"After tabs, defs, adv_var")

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
  do while (.not. End_of_Run)        ! main time-loop , timestep is dt_advec

     nstep=mod(nstep,nmax)+1 !loops from 1 to nmax, between two meteo read

     !FUTURE if (NH3EMIS_VAR) call SetNH3()  ! NH3emis experimental

     select case (mm)
     case(12,1:2);newseason = 1
     case(3:5)   ;newseason = 2
     case(6:8)   ;newseason = 3
     case(9:11)  ;newseason = 4
     endselect

     ! daynumber needed for BCs
     daynumber=day_of_year(yyyy,mm,dd)
     
     if(mm==1 .and. dd==1 .and. hh==0)call Init_nmdays(current_date)!new year starts

     call Code_timer(tim_before)
     if (mm_old /= mm) then   ! START OF NEW MONTH !!!!!
        call Code_timer(tim_before)

        !subroutines/data that must be updated every month
        call readdiss(newseason)

        if (MasterProc.and.DEBUG%MAINCODE ) print *,'maaned og sesong', &
             mm,mm_old,newseason,oldseason

        call Add_2timing(6,tim_after,tim_before,"readdiss, aircr_nox")

        call MetModel_LandUse(2)   ! e.g.  gets snow_flag
        if ( MasterProc .and. DEBUG%MAINCODE ) write(6,*)"vnewmonth start"

        call newmonth

        call Add_2timing(7,tim_after,tim_before,"newmonth")

        if (USE_LIGHTNING_EMIS) call lightning()

        call init_aqueous()

        call Add_2timing(8,tim_after,tim_before,"init_aqueous")
        ! Monthly call to BoundaryConditions.
        if (DEBUG%MAINCODE ) print *, "Into BCs" , me
        ! We set BCs using the specified iyr_trend
        !   which may or may not equal the meteorology year
        call BoundaryConditions(yyyy,mm)
        if (DEBUG%MAINCODE ) print *, "Finished BCs" , me

        !must be called only once, after BC is set
        if(mm_old==0)call Init_massbudget()
        if (DEBUG%MAINCODE ) print *, "Finished Initmass" , me

     endif

     oldseason = newseason
     mm_old = mm

     call Add_2timing(9,tim_after,tim_before,"BoundaryConditions")

     if (DEBUG%MAINCODE ) print *, "1st Infield" , me

     call SetLandUse(daynumber, mm) !daily
     call Add_2timing(11,tim_after,tim_before,"SetLanduse")

     call Meteoread() ! 3-hourly or hourly

     call Add_2timing(10,tim_after,tim_before,"Meteoread")

     call SetDailyBVOC() !daily

     if (USES%FOREST_FIRES) call Fire_Emis(daynumber)

     call Add_2timing(12,tim_after,tim_before,"Fires+BVOC")


     ! if(MasterProc) print "(a,2I2.2,I4,3x,i2.2,a,i2.2,a,i2.2)",' current date and time: ',&
     !   current_date%day,current_date%month,current_date%year,&
     !   current_date%hour, ':',current_date%seconds/60,':',current_date%seconds-60*(current_date%seconds/60)
     if(MasterProc) print "(2(1X,A))",'current date and time:',&
          date2string("YYYY-MM-DD hh:mm:ss",current_date)

     call Code_timer(tim_before)

     call phyche()
     call Add_2timing(14,tim_after,tim_before,"phyche")
 
     call WrtChem()

     call trajectory_in
     call Add_2timing(37,tim_after,tim_before,"massbud,wrtchem,trajectory_in")


     call metfieldint
     call Add_2timing(36,tim_after,tim_before,"metfieldint")



     !this is a bit complicated because it must account for the fact that for instance 3feb24:00 = 4feb00:00 
     ts1=make_timestamp(current_date)
     ts2=make_timestamp(date(enddate(1),enddate(2),enddate(3),enddate(4),0))
     End_of_Run =  (nint(tdif_secs(ts1,ts2))<=0)

  enddo ! time-loop

  call Code_timer(tim_after0)
  call Add_2timing(38,tim_after0,tim_before1,"total within loops")
  call Add_2timing(39,tim_after0,tim_before0,"total")
  !
  !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  call wrtxn(current_date,.true.)
  call massbudget()

  if(MasterProc)then
     print *,'programme is finished'
     ! Gather timing info:
     if(NPROC-1> 0)then
        CALL MPI_RECV(lastptim,8*39,MPI_BYTE,NPROC-1,765,MPI_COMM_WORLD,STATUS,INFO)
     else
        lastptim(:) = mytimm(:)
     endif
     call Output_timing(IO_MYTIM,me,NPROC,nterm,GIMAX,GJMAX)
  elseif(me==NPROC-1) then
     CALL MPI_SEND(mytimm,8*39,MPI_BYTE,0,765,MPI_COMM_WORLD,INFO)
  endif

  ! write 'modelrun.finished' file to flag the end of the FORECAST
  if (MasterProc.and.FORECAST) then
     open(1,file='modelrun.finished')
     close(1)
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
  CALL MPI_FINALIZE(INFO)

end associate   ! yyyy, mm, dd
end program

!===========================================================================
