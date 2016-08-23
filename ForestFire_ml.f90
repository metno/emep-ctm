! < ForestFire_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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

module ForestFire_ml
 !----------------------------------------------------------------
 ! Uses emissions from either:
 !
 ! 1) FINNv1 daily data 2002 - 2011
 ! REFERENCE:
 ! Wiedinmyer, C., Akagi, S. K., Yokelson, R. J., Emmons, L. K., Al-Saadi,
 ! J. A., Orlando, J. J., and Soja, A. J.: The Fire INventory from NCAR (FINN) 
 ! - a high resolution global model to estimate the emissions from open 
 !  burning, Geosci. Model Dev. Discuss., 3, 2439-2476, 
 !   doi:10.5194/gmdd-3-2439-2010, 2010.
 ! http://www.geosci-model-dev-discuss.net/3/2439/2010/gmdd-3-2439-2010.html
 !
 ! 2)  GFED 3 (Global Forest Emission database)
 ! http://www.falw.vu/~gwerf/GFED/
 ! Currently programmed for 8-daily data (available for 2001 - 2007)
 !
 ! 3) GFASv1 Real-Time Fire Emissions
 ! Daily data. Available since 2003 or 2011, depending on version, from MARS
 !   http://www.gmes-atmosphere.eu/about/project_structure/input_data/d_fire/ProductsInMARS/
 ! REFERENCE:
 ! Kaiser, J.W., Heil, A., Andreae, M.O., Benedetti, A., Chubarova, N.,
 !   Jones, L., Morcrette, J.-J., Razinger, M., Schultz, M. G., Suttie, M.,
 !   and van der Werf, G. R.: Biomass burning emissions estimated with a global
 !   fire assimilation system based on observed fire radiative power,
 !   Biogeosciences, 9, 527-554, doi:10.5194/bg-9-527-2012, 2012.
 !----------------------------------------------------------------
  use CheckStop_ml,      only : CheckStop
!CMR  use ChemChemicals_ml,  only : species
!CMR  use ChemSpecs_tot_ml
  use ChemSpecs           
  use GridValues_ml,     only : i_fdom, j_fdom, debug_li, debug_lj, &
                                 debug_proc,xm2,GRIDWIDTH_M
  use Io_ml,             only : PrintLog, datewrite
  use MetFields_ml,      only : z_bnd
  use ModelConstants_ml, only : MasterProc, KMAX_MID, &
                                USES, & !TESTING 
                                DEBUG, FORECAST, &
                                IOU_INST,IOU_HOUR,IOU_HOUR_MEAN, IOU_YEAR
  use NetCDF_ml,         only : ReadField_CDF, Out_netCDF,Real4 ! Reads, writes 
  use OwnDataTypes_ml,   only : Deriv, TXTLEN_SHORT
  use Par_ml,            only : MAXLIMAX, MAXLJMAX, me,limax,ljmax
  use PhysicalConstants_ml, only : AVOG
  use Setup_1dfields_ml, only : rcemis
  use SmallUtils_ml,     only : find_index
 ! No. days per year, date-type :
  use TimeDate_ml,only : nydays, nmdays, date, current_date, day_of_year
  use TimeDate_ExtraUtil_ml,  only: date2string
  implicit none

!  Unimod calls just call Fire_Emis(daynumber)
!  and put the day-testing code here. This lets the module decide if new
!  emissions are needed, and keeps all forest-fire logic here
!

  public :: Fire_Emis
  public :: Fire_rcemis
  private :: Export_FireNc

  logical, public, allocatable, dimension(:,:), save ::  burning
  real, private, allocatable, dimension(:,:,:), save :: BiomassBurningEmis

  integer, private, save ::  ieCO=-1 ! index for CO

  character(len=TXTLEN_SHORT), private :: FF_poll
  integer :: iemep

  !/ Defintions of BB data. If known, we assign the BB pollutant which
  !  corresponds to each possible EMEP emission file. 

  ! Assign mol. wts of the BB data  where known. If mol. wt set to
  ! zero, the code in Fire-rcemis will use the values from the
  ! ChemSpecs_ml, species()%molwt.
  !
  ! If BB doesn't have emissionss, set a "-" for BB, then the
  ! desired emission factor (g/kg DW), and then follow the
  ! example in Fire_setups for NH3 (GFED):
  ! (above text if old. May need update)
  !

  ! =======================================================================
  !  Mapping to EMEP species
  !
  type, private :: bbtype
    character(len=TXTLEN_SHORT) :: BBname
    real :: unitsfac
    real :: frac
    integer :: emep
  endtype bbtype

  ! Here we include the relevant mapping file, which depends on
  ! the source of ffire data and the chemical mechanism (CM)
  !----------------------------------------------
  !=> NBB_DEFS, NEMEPSPECS, FF_defs(NBB_DEFS)

    include 'BiomassBurningMapping.inc' 

  !----------------------------------------------
  ! matrix to get from forest-fire species to EMEP ones

  integer, private, save :: emep_used(NEMEPSPECS) = 0
  real   , private, save :: sum_emis(NEMEPSPECS) = 0

  ! =======================================================================



contains
subroutine Fire_Emis(daynumber)
!.....................................................................
!**    DESCRIPTION:
!    Reads forest-fire emissions. So far set up for GFED 8d, but in
!    principal we can re-code by simply adding alternative
!    subroutines, e.g. to cope with pre-2001 monthly emissions

  integer, intent(in) :: daynumber   ! daynumber (1,..... ,365)


  real,allocatable :: rdemis(:,:)  ! Emissions read from file
  integer :: i,j,nstart, alloc_err, iBB, n
  logical, save :: my_first_call = .true.   ! DSFF
  logical :: my_first_defs = .true. 
  integer, save  :: dd_old = -1, mm_old=-1 
  real    :: fac, to_kgm2s   

  integer :: ind, ne
  integer :: loc_maxemis(2) ! debug

  character(len=*), parameter :: &
    GFED_PATTERN = 'GFED_ForestFireEmis.nc',&
    FINN_PATTERN = 'FINN_ForestFireEmis_YYYY.nc',&
    GFAS_PATTERN = 'GFAS_ForestFireEmis_YYYY.nc'
  character(len=len(GFAS_PATTERN)) :: fname = ''
  logical :: my_debug=.false.,fexist=.false.
  integer, parameter :: verbose = 1
  integer :: dd1, dd2 ! TESTING
  real,allocatable :: xrdemis(:,:)  !TESTING
  !
  integer :: yyyy, mm  ! instead of associate

  !FAILED with ifort 12.0 :-(
  !ACDATES: associate ( yyyy => current_date%year, mm => current_date%month )
  yyyy = current_date%year
  mm = current_date%month 

  if(my_first_call) &
    call PrintLog("Biomass Mapping: "//trim(BiomassBurningMapping),MasterProc)

  select case(verbose)
    case(:0);my_debug=.false.
    case(1) ;my_debug=DEBUG%FORESTFIRE.and.MasterProc.and.my_first_call
    case(2) ;my_debug=DEBUG%FORESTFIRE.and.MasterProc
    case(3) ;my_debug=DEBUG%FORESTFIRE
    case(4:);my_debug=.true.
  endselect

  nstart = -1 ! reset for GFED
  select case(BiomassBurningMapping(1:4))

  case("GFED") ! 8-day values

    if(DEBUG%FORESTFIRE.and.MasterProc) write(*,*) "FIRE selects GFED"
    select case(yyyy)
    case(2001:2007)
      if(MasterProc)&
        write(*,*) "WARNING! FFIRE GFED USED! May not be working properly check results!"
    case default
      if(my_first_call)&
        call PrintLog("8d GFED Forest Fires: only between 2001--2007",MasterProc)
      call CheckStop("GFED not available. Use other FF data, or set USE_FOREST_FIRES .false. in ModelConstants")
      my_first_call = .false.
      return
    endselect
    if(DEBUG%FORESTFIRE.and.MasterProc) &
      write(*,*) "GFED FIRE days:", yyyy, daynumber, dd_old, mod(daynumber,8), my_first_call

    ! GFED Fire emissions are called at 8 days intervals (1, 9, 17, ....)
    ! 46 values available each year: day 361 is the last one.
    ! Return unless new period

    if(.not.my_first_call.and.mod(daynumber,8)/= 1) return
    nstart=( yyyy -2001)*46+(daynumber+7)/8

  case("FINN")
    !write(*,*) "FIRE selects FINN",me, DEBUG%FORESTFIRE, MasterProc
    if(DEBUG%FORESTFIRE.and.MasterProc) write(*,*) "FIRE selects FINN"

  case("GFAS")
    if(DEBUG%FORESTFIRE.and.MasterProc) write(*,*) "FIRE selects GFAS"

  case default
    call CheckStop("Unknown B.B.Mapping: "//trim(BiomassBurningMapping))
  endselect


  if(DEBUG%FORESTFIRE.and.MasterProc) &
    write(*,*) "Starting FIRE days:", yyyy, &
      daynumber, dd_old, mod(daynumber,8), my_first_call

  if(dd_old==daynumber) return   ! Only calculate once per day max
  dd_old = daynumber

  if(my_first_call)then

    allocate(BiomassBurningEmis(NEMEPSPECS,MAXLIMAX,MAXLJMAX),&
             burning(MAXLIMAX,MAXLJMAX),stat=alloc_err)
    call CheckStop(alloc_err,"ForestFire BiomassBurningEmis alloc problem")
    my_first_call = .false.
    ne = 0     ! number-index of emep species

    do n=1, NBB_DEFS            ! Only unique EMEP SPECS in emep_used
      iemep = FF_defs(n)%emep 
      if(find_index(iemep,emep_used(:))>0) cycle
      
      ne = ne + 1
      emep_used(ne) = iemep

      ! CO is special. Keep the index
      if(species(iemep)%name=="CO") ieCO=ne

      if(MasterProc) write(*,"(a,2i4,a17)") "FFIRE Mapping EMEP ", &
        ne, iemep, trim(species(iemep)%name)
    enddo !n
    call CheckStop(ieCO<0,"No mapping for 'CO' found on "//BiomassBurningMapping)
    call CheckStop(any(emep_used<0),"UNSET FFIRE EMEP "//BiomassBurningMapping)

  endif !my first call

! Check if ForestFire file exists
  select case(BiomassBurningMapping(1:4))
    case("GFED");fname=date2string(GFED_PATTERN,current_date)
    case("FINN");fname=date2string(FINN_PATTERN,current_date)
    case("GFAS");fname=date2string(GFAS_PATTERN,current_date)
  endselect
  inquire(file=fname,exist=fexist)
  if(.not.fexist)then
    if(MasterProc)then
      write(*,*)"ForestFire file not found: "//trim(fname)
      call CheckStop(.not.FORECAST,"Missing ForestFire file")
    endif
    burning(:,:) = .false.
    return
  endif

  if(DEBUG%FORESTFIRE.and.MasterProc) write(*,*) "FOREST_FIRE: ", daynumber,nstart
  BiomassBurningEmis(:,:,:) = 0.0
  allocate(rdemis(MAXLIMAX,MAXLJMAX),stat=alloc_err)
  call CheckStop(alloc_err,"ForestFire rdemis alloc problem")
  if(USES%MONTHLY_FF.and.mm/=mm_old) then
    if( MasterProc ) write(*,*) "Start monthly FF ", mm
    allocate(xrdemis(MAXLIMAX,MAXLJMAX),stat=alloc_err)
  end if
  
  ! We need to look for forest-fire emissions which are equivalent
  ! to the standard emission files:
  do iBB = 1, NBB_DEFS
    FF_poll = FF_defs(iBB)%BBname
    iemep   = FF_defs(iBB)%emep  ! 
    ind     = find_index( iemep, emep_used )  !  Finds 1st emep in BiomassBurning

    if(DEBUG%FORESTFIRE.and.MasterProc) &
      write(*,"( a,3i5, a8,i3)") "FIRE SETUP: ", iBB,iemep,ind, &
        trim(FF_poll), len_trim(FF_poll)

   ! FORECAST mode: if file/variable/timestep not found it should not crash
    rdemis(:,:)=0.0

    select case(BiomassBurningMapping(1:4))
    case("GFED")
      if(my_debug) &
        write(*,*) "FFIRE GFED ", me, iBB, nstart,  trim(FF_poll), trim(fname)
      call ReadField_CDF(fname,FF_poll,rdemis,nstart,interpol='zero_order',&
        needed=.not.FORECAST,UnDef=0.0,debug_flag=DEBUG%FORESTFIRE)
      !unit conversion to GFED [g/m2/8day]->[kg/m2/s]
      to_kgm2s = 1.0e-3 /(8*24.0*3600.0)
      forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*to_kgm2s

    case("FINN")
      if(my_debug) &
        write(*,*) "FFIRE FINN ", me, iBB, daynumber,  trim(FF_poll), trim(fname)
      if( USES%MONTHLY_FF ) then
   !associate ( yyyy => current_date%year, mm => current_date%month )
         dd1=day_of_year(yyyy,mm,1)
         dd2=day_of_year(yyyy,mm+1,1) - 1  ! Ok for Dec also, gets 367-1.
   !end associate ! yyyy, mm

         rdemis = 0.0
         do i =  1, dd1, dd2
            call ReadField_CDF(fname,FF_poll,xrdemis,daynumber,&
               interpol='mass_conservative',&
               needed=.not.FORECAST,UnDef=0.0,debug_flag=DEBUG%FORESTFIRE)
            rdemis = rdemis + xrdemis/real(dd2-dd1+1)  ! Get monthly avg.
         end do
      else
         call ReadField_CDF(fname,FF_poll,rdemis,daynumber,&
              interpol='mass_conservative',&
              needed=.not.FORECAST,UnDef=0.0,debug_flag=DEBUG%FORESTFIRE)
              
      end if ! USES%MONTHLY_FF 

      !else
      !call ReadField_CDF(fname,FF_poll,rdemis,daynumber,interpol='mass_conservative',&
      !  needed=.not.FORECAST,UnDef=0.0,debug_flag=DEBUG%FORESTFIRE)
      !end if
      !MONTHLY FF
      !unit conversion to FINN: Can be negative if REMPPM to be calculated
      !unit conversion to FINN: Can be negative if REMPPM to be calculated
      fac=FF_defs(iBB)%unitsfac * FF_defs(iBB)%frac  ! --> [kg/day]
      fac=fac/(GRIDWIDTH_M*GRIDWIDTH_M*24.0*3600.0) ! [kg/day]->[kg/m2/s]
      forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*fac*xm2(i,j)

    case("GFAS")
      nstart = daynumber
! something more sophisticated is needed for YYYY_ or YYYYMM_ files,
! e.g. use ReadTimeCDF and nctime2idate/idate2nctime to find the right record:
!  nstart=FindTimeCDFRecord(fname,current_date,prec_ss=3600.0*12)
      if(my_debug) &
        write(*,*) "FFIRE GFAS ", me, iBB, n, nstart,  trim(FF_poll), trim(fname)
      call ReadField_CDF(fname,FF_poll,rdemis,nstart,interpol='conservative',&
          needed=.not.FORECAST,debug_flag=DEBUG%FORESTFIRE)
      ! GFAS units are [kg/m2/s]. No further unit conversion is needed.
      ! However, fac can be /=1, e.g. when REMPPM is calculated
      fac=FF_defs(iBB)%unitsfac * FF_defs(iBB)%frac
      if(fac/=1.0) forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*fac
    end select


   ! Assign . units should be [kg/m2/s] here 
    forall(j=1:ljmax,i=1:limax) 
      BiomassBurningEmis(ind,i,j) = BiomassBurningEmis(ind,i,j) + rdemis(i,j) 
    endforall

    if(my_debug)  write(*,"(3a10,i4,f8.3,es12.3)") "FFIRE SUMS:", &
      trim(FF_poll), trim( species(iemep)%name), ind, &
      species(iemep)%molwt, sum( BiomassBurningEmis(ind,:,:) )


    if(my_first_defs) call PrintLog(&
      "ForestFire_ml :: Assigns "//trim(FF_poll) , MasterProc)

    if(DEBUG%FORESTFIRE) sum_emis(ind)=sum_emis(ind)+sum(BiomassBurningEmis(ind,:,:))
  enddo ! BB_DEFS

  my_first_defs  = .false.
  deallocate(rdemis)
  if(USES%MONTHLY_FF.and.mm/=mm_old) then
     deallocate(xrdemis)
     mm_old = mm
  end if

  ! For cases where REMPPM25 s derived as the difference between PM25 and (BC+1.7*OC)
  ! we need some safety:

  BiomassBurningEmis(:,:,:) = max( BiomassBurningEmis(:,:,:), 0.0 )

  ! Logical to let Unimod know if there is any emission here to worry about
  burning(:,:) = ( BiomassBurningEmis(ieCO,:,:) > 1.0e-19 )



  ! Some databases (e.g. FINN, GFED) have both total PM25 and EC, OC. The difference
  ! REMPPM25, is created by the BiomasBurning mapping procedure, but we just
  ! check here
  if(DEBUG%FORESTFIRE.and.debug_proc) then
    n = ieCO
    loc_maxemis = maxloc(BiomassBurningEmis(n,:,: ) )

    associate ( idbg=>loc_maxemis(1), jdbg=>loc_maxemis(2) )

    write(*,"(a,i4,i3,2i4,2i5,es12.3, 2i4)") "SUM_FF CHECK ME: ",  daynumber, me, loc_maxemis, &
         i_fdom(idbg), j_fdom(jdbg), BiomassBurningEmis(n,idbg,jdbg), debug_li,debug_lj

    call datewrite("SUM_FF CHECK CO: ",  &
      (/ daynumber, n, i_fdom( idbg ), j_fdom( jdbg ) /) ,&
      (/  sum_emis(n), maxval(BiomassBurningEmis(n,:,: ) ), &
          BiomassBurningEmis(n,debug_li,debug_lj) /) )
    end associate ! idbg, jdbg
  endif ! debug_proc
  !end associate ACDATES
end subroutine Fire_Emis

!=============================================================================

subroutine Fire_rcemis(i,j)
!  Disperses the fire emissions vertically and converts to molecules/cm3/s.

!// Injection height: here over 8 levels. Alternative could be PBL
!   or  equally upto ca. 2*PBL (suggested by Sofiev, GEMS)
!   QUERY - should the emissions be divided equally by level?
!   - will give a higher mixing ratio for thinner levels

  integer, intent(in) :: i,j

  integer, parameter :: KEMISFIRE = 12
  real, dimension(KEMISFIRE:KMAX_MID) :: invDeltaZfac !  height of layer in m div 9
  integer ::  k, n, iem

  integer ::  N_LEVELS  ! = 9.0 here

  real    :: origrc, bbe, fac
  logical :: debug_flag

  debug_flag = (DEBUG%FORESTFIRE.and.debug_proc .and.&
                i==debug_li.and.j==debug_lj)
  if(debug_flag.and.BiomassBurningEmis(ieCO,i,j) > 1.0e-10)  &
    write(*,"(a,5i4,es12.3,f9.3)") "BurningDEBUG ", me, i,j, &
      i_fdom(i), j_fdom(j), BiomassBurningEmis(ieCO,i,j)

  N_LEVELS = KMAX_MID - KEMISFIRE + 1 

  !// last conversion factors:
  ! The biomassBurning array is kept in kg/m2/s for consistency with other
  ! emissions. We here convert to molecules/cm3/s after spreading
  ! through a vertical distance dz
  !
  ! If we had E in kg/m2/s, we would then take
  !  E*1.0e3  -> g/m2/s
  !  E*0.1    -> g/cm2/s
  !  E*0.1 /MW * Av -> molec/cm2/s
  !  E*0.001 /MW * Av / DZ -> molec/cm3/s where DZ is spread in m
  !  i.e. fmap should be 0.001*Av/MW
  !  (plus account for the fraction of the inventory assigned to EMEP species)


  !/ Here we just divide by the number of levels. Biased towards
  !  different levels since thickness and air content differ. Simple though.

  do k = KEMISFIRE, KMAX_MID
    invDeltaZfac(k) = 1.0/ (z_bnd(i,j,k) - z_bnd(i,j,k+1)) /N_LEVELS
  enddo
 
  do n = 1, NEMEPSPECS 
    iem = emep_used(n)
    origrc = rcemis(iem,KMAX_MID)   ! just for printout
    fac =  0.001 * AVOG /species(iem)%molwt    ! MW scale if needed

    ! distribute vertically:
    do k = KEMISFIRE, KMAX_MID
      rcemis(iem,k) = rcemis(iem,k) + BiomassBurningEmis(n,i,j)*invDeltaZfac(k)*fac
    enddo !k

    if(debug_flag) then
      k=KMAX_MID
      write(*,"(a,2i3,1x,a8,i4,es10.2,4es10.2)") "FIRERC ",&
        n, iem, trim(species(iem)%name), k, BiomassBurningEmis(iem,i,j),&
        invDeltaZfac(k), origrc, rcemis(iem,k)
    endif

!DSBB    !--  Add up emissions in ktonne ......
!DSBB    !   totemadd(iem) = totemadd(iem) + &
!DSBB    !   tmpemis(iqrc) * dtgrid * xmd(i,j)

  enddo ! n
 !       call Export_FireNc() ! Caused problems on last attempt

endsubroutine Fire_rcemis
!=============================================================================
subroutine Export_FireNc()
  type(Deriv) :: def1 ! definition of fields
 
  def1%class='ForestFireEmis' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=1.0        !not used
  def1%name='CO'        !written
  def1%unit='g/m2'      !written
  def1%name='CO_zero'
  def1%name='CO_ASCII'       

  call Out_netCDF(IOU_INST,def1,2,1, BiomassBurningEmis(ieCO,:,:),1.0,&
                  CDFtype=Real4,fileName_given='FF.nc')
endsubroutine Export_FireNc

endmodule ForestFire_ml
!=============================================================================
