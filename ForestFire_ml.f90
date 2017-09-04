! <ForestFire_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.15>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2017 met.no
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
!  Unimod calls just call Fire_Emis(daynumber)
!  and put the day-testing code here. This lets the module decide if new
!  emissions are needed, and keeps all forest-fire logic here
!----------------------------------------------------------------
use CheckStop_ml,         only: CheckStop,CheckNC
use ChemSpecs
use GridValues_ml,        only: i_fdom, j_fdom, debug_li, debug_lj, &
                                debug_proc,xm2,GRIDWIDTH_M
use Io_ml,                only: PrintLog, datewrite, IO_NML
use MetFields_ml,         only: z_bnd
use ModelConstants_ml,    only: MasterProc, KMAX_MID, DEBUG, IOU_INST
use netcdf,               only: nf90_open, nf90_nowrite, nf90_close
use NetCDF_ml,            only: ReadTimeCDF,ReadField_CDF,Out_netCDF,Real4,&
                                closedID
use OwnDataTypes_ml,      only: Deriv, TXTLEN_SHORT
use Par_ml,               only: LIMAX, LJMAX, me,limax,ljmax
use PhysicalConstants_ml, only: AVOG
use Setup_1dfields_ml,    only: rcemis
use SmallUtils_ml,        only: find_index
! No. days per year, date-type:
use TimeDate_ml,          only: current_date,day_of_year,max_day
use TimeDate_ExtraUtil_ml,only: date2string,nctime2string,date2nctime,date2file

implicit none
private
public :: Fire_Emis, Fire_rcemis, burning

logical, allocatable, dimension(:,:), save :: burning
real,  allocatable, dimension(:,:,:), save :: BiomassBurningEmis

integer, save :: ieCO=-1 ! index for CO

character(len=TXTLEN_SHORT) :: FF_poll
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
type :: bbtype
  character(len=TXTLEN_SHORT) :: BBname
  real :: unitsfac
  real :: frac
  integer :: emep
end type bbtype

! Here we include the relevant mapping file, which depends on
! the source of ffire data and the chemical mechanism (CM)
!----------------------------------------------
!=> NBB_DEFS, NEMEPSPECS, FF_defs(NBB_DEFS)

  include 'BiomassBurningMapping.inc' 

!----------------------------------------------
! matrix to get from forest-fire species to EMEP ones

integer, save :: emep_used(NEMEPSPECS) = 0
real,    save :: sum_emis(NEMEPSPECS) = 0

character(len=4), parameter :: BBMAP=BiomassBurningMapping(1:4)
character(len=TXTLEN_SHORT) :: MODE="DAILY_REC"

integer, parameter :: &
  max_string_length=200 ! large enough for paths to be set on Fire_config namelist
character(len=max_string_length), save :: &
  GFED_PATTERN = 'GFED_ForestFireEmis.nc',&
  FINN_PATTERN = 'FINN_ForestFireEmis_v15_YYYY.nc',&
  GFAS_PATTERN = 'GFAS_ForestFireEmis_YYYY.nc'

! interpolation method in ReadField_CDF
character(len=30), save :: bbinterp = '-'

! Notes on interpolation choices: (from NetCDF_ml)
  !'zero_order' gives value at closest gridcell. Probably good enough for most
  ! applications.  Does not smooth out values
  !'conservative' and 'mass_conservative' give smoother fields and
  !are approximatively integral conservative (integral over a region is
  !conserved). The initial gridcells are subdivided into smaller subcells
  !and each subcell is assigned to a cell in the model grid
  !'conservative' can be used for emissions given in kg/m2 (or kg/m2/s)
  !or landuse or most fields.  The value in the netcdf file and in
  !model gridcell are of the similar.  'mass_conservative' can be used
  !for emissions in kg (or kg/s). If the gricell in the model are !twice
  !as small as the gridcell in the netcdf file, the values will also be
  !reduced by a factor 2.

integer, save ::    &
  verbose=1,        & ! debug verbosity 0,..,4
  persistence=1,    & ! persistence in days
  fire_year=-1,     & ! override current year
  nread=-1            ! records in forest fire file
logical, save ::    &
  need_file=.true., & ! stop if don't find file
  need_date=.true., & ! stop if don't find time record
  need_poll=.true., & ! stop if don't find pollutant
  debug_level(-1:5)=.false.
! =======================================================================
contains
subroutine Config_Fire()
  logical, save :: first_call=.true.
  integer :: ios, ne, n
  NAMELIST /Fire_config/MODE,verbose,persistence,fire_year,&
                        need_file,need_date,need_poll,&
                        GFED_PATTERN,FINN_PATTERN,GFAS_PATTERN
                        
  if(.not.first_call)return
  call PrintLog("Biomass Mapping: "//trim(BiomassBurningMapping),MasterProc)

  if(DEBUG%FORESTFIRE.and.MasterProc) write(*,*) "FIRE selects ",BBMAP
  select case(BBMAP)
    case("GFED")
       persistence=8  ! 8-day records
       bbinterp = 'conservative'
    case("FINN")
       persistence=1  ! 1-day records
       bbinterp = 'mass_conservative'
    case("GFAS")
       persistence=3  ! 1-day records, valid for 3 day in FORECAST mode
       bbinterp = 'conservative'
    case default;call CheckStop("Unknown B.B.Mapping")
  end select

  rewind(IO_NML)
  read(IO_NML,NML=Fire_config,iostat=ios)
  call CheckStop(ios,"NML=Fire_config")  
  if(DEBUG%FORESTFIRE.and.MasterProc)then
    write(*,*) "NAMELIST IS "
    write(*,NML=Fire_config)
  end if
  
  ! set vebosity levels
  verbose=min(max(0,verbose),4) ! debug verbosity 0,..,4
  debug_level(:0)=.false.
  debug_level(1)=DEBUG%FORESTFIRE.and.MasterProc.and.first_call
  debug_level(2)=DEBUG%FORESTFIRE.and.MasterProc
  debug_level(3)=DEBUG%FORESTFIRE
  debug_level(4:)=.true.

  allocate(BiomassBurningEmis(NEMEPSPECS,LIMAX,LJMAX),&
           burning(LIMAX,LJMAX),stat=ios)
  call CheckStop(ios,"ForestFire BiomassBurningEmis alloc problem")
  ne = 0     ! number-index of emep species

  do n=1, NBB_DEFS              ! Only unique EMEP SPECS in emep_used
    iemep = FF_defs(n)%emep 
    if(find_index(iemep,emep_used(:))>0) cycle
    
    ne = ne + 1
    emep_used(ne) = iemep

    ! CO is special. Keep the index
    if(species(iemep)%name=="CO") ieCO=ne

    if(MasterProc) write(*,"(a,2i4,a17)") "FFIRE Mapping EMEP ", &
      ne, iemep, trim(species(iemep)%name)
  end do !n
  call CheckStop(ieCO<1,"No mapping for 'CO' found on "//BiomassBurningMapping)
  call CheckStop(any(emep_used<1),"UNSET FFIRE EMEP "//BiomassBurningMapping)

  first_call=.false.
end subroutine Config_Fire

subroutine Fire_Emis(daynumber)
!.....................................................................
!**    DESCRIPTION:
!    Reads forest-fire emissions. So far set up for GFED 8d, but in
!    principal we can re-code by simply adding alternative
!    subroutines, e.g. to cope with pre-2001 monthly emissions

  integer, intent(in) :: daynumber   ! daynumber (1,..... ,365)


  real,allocatable :: rdemis(:,:)  ! Emissions read from file
  integer :: i,j,nstart, alloc_err, iBB, n
  logical, save :: first_call = .true.
  integer, save  :: nn_old=-1
  real    :: fac, to_kgm2s   

  integer :: ind, ncFileID=closedID
  integer :: loc_maxemis(2) ! debug

  character(len=max_string_length) :: fname='new'
  logical :: debug_ff=.false.,debug_nc=.false.
  real, allocatable :: xrdemis(:,:) ! MODE=*_AVG
  integer :: dn1, dn2, ndn          ! MODE=*_AVG
  integer :: yyyy, mm, dd

  if(first_call) call Config_Fire()
  debug_ff=debug_level(verbose)
  debug_nc=debug_level(verbose-1)

  yyyy = current_date%year
  mm   = current_date%month 
  dd   = current_date%day
  if(fire_year>1) yyyy=fire_year
  
  select case(MODE)
    case("DAILY_REC","D","d")
      if(nn_old==daynumber) return  ! Calculate once per day
      nn_old=daynumber
      dn1=dd
      dn2=dd
    case("MONTHLY_AVG","M","m")
      if(nn_old==mm) return         ! Calculate once per month
      nn_old=mm
      dn1=day_of_year(yyyy,mm,01)
      dn2=dn1+max_day(mm,yyyy)-1
    case("YEARLY_AVG","Y","y")
      if(nn_old==yyyy) return       ! Calculate once per year
      nn_old=yyyy
      dn1=day_of_year(yyyy,01,01)
      dn2=day_of_year(yyyy,12,31)
    case default
      call CheckStop("Unknown ForestFire MODE="//trim(MODE))
  end select

  if(dn1<dn2)then
    allocate(xrdemis(LIMAX,LJMAX),stat=alloc_err)
    nstart=daynumber    ! for debug info
  else
    ! newFFrecord: has pollutant|fname|record changed since last call?
    FF_poll=""
    if(.not.newFFrecord([yyyy,mm,dd]))&
      return                        ! Continue if new record||file  
  end if

  if(DEBUG%FORESTFIRE.and.MasterProc) then
    write(*,*) "Starting FIRE MODE=",trim(MODE),&
      date2string(" YYYY-MM-DD",[yyyy,mm,dd]),first_call,debug_ff,debug_nc
    write(*,*) 'FFIRE Interp= ', trim(bbinterp), dn1, dn2
  end if

  BiomassBurningEmis(:,:,:) = 0.0
  allocate(rdemis(LIMAX,LJMAX),stat=alloc_err)
  call CheckStop(alloc_err,"ForestFire rdemis alloc problem")
  
  ! We need to look for forest-fire emissions which are equivalent
  ! to the standard emission files:
  do iBB = 1, NBB_DEFS
    FF_poll = FF_defs(iBB)%BBname
    iemep   = FF_defs(iBB)%emep  ! 
    ind     = find_index( iemep, emep_used ) !Finds 1st emep in BiomassBurning

    if(DEBUG%FORESTFIRE.and.MasterProc) then 
      write(*,"( a,3i5, a8,i3)") "FIRE SETUP: ", iBB,iemep,ind, &
        trim(FF_poll), len_trim(FF_poll)
      if(debug_ff) &
        write(*,*) 'FFIRE ',BBMAP,':',me,iBB,nstart,trim(FF_poll),trim(fname)
    end if

   ! if(.not.need_file|time|poll) continue if file|time|poll is not found
    rdemis(:,:)=0.0


!--------- Aug 2017: methods merged. Keep UnDef=0 for future safety
    if(dn1<dn2)then
        rdemis = 0.0
        ndn=0
        do dd = dn1, dn2
          if(newFFrecord([yyyy,01,dd])) &
          call ReadField_CDF(fname,FF_poll,xrdemis,nstart,interpol=bbinterp,&
             needed=need_poll,UnDef=0.0,debug_flag=debug_nc,&
             ncFileID_given=ncFileID)
          rdemis = rdemis + xrdemis                 ! month total
          ndn    = ndn + 1
        end do
    else
        ndn=1
        call ReadField_CDF(fname,FF_poll,rdemis,nstart,interpol=bbinterp,&
          needed=need_poll,UnDef=0.0,debug_flag=debug_nc,&
          ncFileID_given=ncFileID)
    end if
!-------- Aug 2017

    select case(BBMAP)
    case("GFED")

     !unit conversion, GFED [g/m2/8day]->[kg/m2/s]
      to_kgm2s = 1.0e-3 /(8*24.0*3600.0)
      if(ndn>1) to_kgm2s=to_kgm2s/ndn               ! total-->avg.
      forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*to_kgm2s

    case("FINN")

     ! unit conversion, FINN [mole/day]->[kg/m2/s]
     ! (Can be negative if REMPPM to be calculated)
      fac=FF_defs(iBB)%unitsfac * FF_defs(iBB)%frac ! --> [kg/day]
      fac=fac/(GRIDWIDTH_M*GRIDWIDTH_M*24.0*3600.0) ! [kg/day]->[kg/m2/s]
      if(ndn>1) fac=fac/ndn                         ! total-->avg.
      forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*fac*xm2(i,j)

    case("GFAS")

     ! GFAS units are [kg/m2/s]. No further unit conversion is needed.
     ! However, fac can be /=1, e.g. when REMPPM is calculated
      fac=FF_defs(iBB)%unitsfac * FF_defs(iBB)%frac
      if(ndn>1) fac=fac/ndn                         ! total-->avg.
      if(fac/=1.0) forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*fac
    end select

   ! Assign . units should be [kg/m2/s] here 
    forall(j=1:ljmax,i=1:limax) &
      BiomassBurningEmis(ind,i,j) = BiomassBurningEmis(ind,i,j) + rdemis(i,j) 

    if(debug_ff.and. debug_proc) &
       write(*,"(3a10,2i4,f8.3,es12.3)") "FFIRE SUMS:", &
        trim(FF_poll), trim( species(iemep)%name), me, ind, &
        species(iemep)%molwt, sum( BiomassBurningEmis(ind,:,:) )

    call PrintLog("ForestFire_ml :: Assigns "//trim(FF_poll),&
      first_call.and.MasterProc)

    if(DEBUG%FORESTFIRE) sum_emis(ind)=sum_emis(ind)+&
          sum(BiomassBurningEmis(ind,:,:))
  end do ! BB_DEFS

  ! have to close the file here
  call CheckNC(nf90_close(ncFileID),"close:"//trim(fname))
  ncFileID=closedID

  first_call  = .false.
  deallocate(rdemis)
  if(allocated(xrdemis)) deallocate(xrdemis)

  ! For cases where REMPPM25 s derived as the difference between PM25 and
  !  (BC+1.7*OC) we need some safety:

  BiomassBurningEmis(:,:,:) = max( BiomassBurningEmis(:,:,:), 0.0 )

  ! Logical to tell if there is any emission here to worry about
  burning(:,:) = ( BiomassBurningEmis(ieCO,:,:) > 1.0e-19 )



  ! Some databases (e.g. FINN, GFED) have both total PM25 and EC, OC. The
  ! difference, REMPPM25, is created by the BiomasBurning mapping procedure,
  ! but we just check here
  if(DEBUG%FORESTFIRE.and.debug_proc) then
    n = ieCO
    loc_maxemis = maxloc(BiomassBurningEmis(n,:,: ) )

    associate ( idbg=>loc_maxemis(1), jdbg=>loc_maxemis(2) )

    write(*,"(a,i4,i3,2i4,2i5,es12.3, 2i4)") "SUM_FF CHECK ME: ",  daynumber,&
       me, loc_maxemis, i_fdom(idbg), j_fdom(jdbg),&
         BiomassBurningEmis(n,idbg,jdbg), debug_li,debug_lj

    call datewrite("SUM_FF CHECK CO: ",  &
      (/ daynumber, n, i_fdom( idbg ), j_fdom( jdbg ) /) ,&
      (/  sum_emis(n), maxval(BiomassBurningEmis(n,:,: ) ), &
          BiomassBurningEmis(n,debug_li,debug_lj) /) )
    end associate ! idbg, jdbg
  end if ! debug_proc
  !end associate ACDATES

contains

function newFFrecord(ymd) result(new)
  integer, intent(in) :: ymd(3)
  logical :: new

  character(len=TXTLEN_SHORT), save      :: poll_old=''
  character(len=max_string_length), save :: file_old=''
  integer, save                          :: record_old=-1
  real, dimension(366), save :: fdays=-1
  logical :: fexist=.false.
  real :: ncday(0:1)

  ! Check: New file
  select case(BBMAP)
    case("GFED");fname=date2file(GFED_PATTERN,ymd,persistence-1,"days")
    case("FINN");fname=date2file(FINN_PATTERN,ymd,persistence-1,"days")
    case("GFAS");fname=date2file(GFAS_PATTERN,ymd,persistence-1,"days")
  end select
  if(fname/=file_old)then
    if(DEBUG%FORESTFIRE.and.MasterProc) &
      write(*,*)"ForestFire new file: ",trim(fname)
  ! close old ncFile, if already open
    if(ncFileID/=closedID)&
      call CheckNC(nf90_close(ncFileID),"close:"//trim(file_old))
    ncFileID=closedID
  ! check if new file exists
    inquire(file=fname,exist=fexist)    ! check if fname exixts
    if(.not.fexist)then
      if(MasterProc)then
        write(*,*)"ForestFire file not found: ",trim(fname)
        call CheckStop(need_file,"Missing ForestFire file")
      end if
      burning(:,:) = .false.
      new=.false.
      return
    end if
  ! read all times records in fname, and process them
    nread=-1                            
    fdays(:)=-1.0                       
    call ReadTimeCDF(fname,fdays,nread) 
    record_old=-1                       
  end if

  ! Check: New pollutant
  if(FF_poll/=poll_old)then
    if(DEBUG%FORESTFIRE.and.MasterProc) &
      write(*,*)"ForestFire new pollutant: ",trim(FF_poll)
  end if

  ! Check: New time record
  call date2nctime(ymd,ncday(1))
  ncday(0)=ncday(1)-persistence+1
  nstart=MAXLOC(fdays(:nread),DIM=1,&
    MASK=(fdays(:nread)>=ncday(0)).and.(fdays(:nread)<(ncday(1)+1.0)))
  if(nstart/=record_old)then
    if(DEBUG%FORESTFIRE.and.MasterProc) &
      write(*,*)"ForestFire new record: ",&
        nstart,nctime2string("(YYYY-MM-DD hh:mm)",fdays(nstart))
    if((fdays(nstart)<ncday(0)).or.(fdays(nstart)>=(ncday(1)+1.0)))then
      if(MasterProc)then
        write(*,*)"ForestFire: no records between ",&
          nctime2string("YYYY-MM-DD 00:00",ncday(0))," and ",&
          nctime2string("YYYY-MM-DD 23:59",ncday(1))
        call CheckStop(need_date,"Missing ForestFire records")
      end if
      burning(:,:) = .false.
      new=.false.
      return
    end if
  end if
  
  ! Update if new
  new=(fname/=file_old).or.(nstart/=record_old).or.(FF_poll/=poll_old)
  if(new)then
    if(ncFileID==closedID) & ! open the file only once
      call CheckNC(nf90_open(fname,nf90_nowrite,ncFileID),"open:"//trim(fname))
    file_old=fname
    poll_old=FF_poll
    record_old=nstart
  end if
end function newFFrecord
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
    write(*,"(a,5i4,es12.3,f9.3)") "FIREBurningDEBUG ", me, i,j, &
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
  end do
 
  do n = 1, NEMEPSPECS 
    iem = emep_used(n)
    origrc = rcemis(iem,KMAX_MID)   ! just for printout
    fac =  0.001 * AVOG /species(iem)%molwt    ! MW scale if needed

    ! distribute vertically:
    do k = KEMISFIRE, KMAX_MID
      rcemis(iem,k) = rcemis(iem,k) + BiomassBurningEmis(n,i,j)*invDeltaZfac(k)*fac
    end do !k

    if(debug_flag) then
      k=KMAX_MID
      write(*,"(a,2i3,1x,a8,i4,es10.2,4es10.2)") "FIRERC ",&
        n, iem, trim(species(iem)%name), k, BiomassBurningEmis(iem,i,j),&
        invDeltaZfac(k), origrc, rcemis(iem,k)
    end if

!DSBB    !--  Add up emissions in ktonne ......
!DSBB    !   totemadd(iem) = totemadd(iem) + &
!DSBB    !   tmpemis(iqrc) * dtgrid * xmd(i,j)

  end do ! n
 !       call Export_FireNc() ! Caused problems on last attempt
end subroutine Fire_rcemis
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
end subroutine Export_FireNc

endmodule ForestFire_ml
!=============================================================================
