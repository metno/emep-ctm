! <ForestFire_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.0>
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
module ForestFire_mod
!----------------------------------------------------------------
! Uses emissions from either:
!
! 0) FINNv2.5
!      mod    - 2002-2021
!      modvrs - 2012 - 2021
!    Converted  (by Qing Mu) from 
!MODIS
!wget https://www.acom.ucar.edu/Data/fire/data/finn2/FINNv2.5_mod_GEOSCHEM_2019_c20211213.txt.gz
!MODIS+VIIRS
!wget https://www.acom.ucar.edu/Data/fire/data/finn2/FINNv2.5_modvrs_GEOSCHEM_2019_c20211213.txt.gz
! REFERENCES:
!Wiedinmyer, C.; Kimura, Y.; McDonald-Buller, E. C.; Emmons, L. K.; Buchholz, R.
!R.; Tang, W.; Seto, K.; Joseph, M. B.; Barsanti, K. C.; Carlton, A. G. &
!Yokelson, R. The Fire Inventory from NCAR version 2.5: an updated global fire
!emissions model for climate and chemistry applications EGUsphere, 2023, 2023,
!1-45, https://egusphere.copernicus.org/preprints/egusphere-2023-124/ 

! 1) FINNv1.5 daily data 2002 - 2015
!  http://bai.acom.ucar.edu/Data/fire/
! REFERENCES:
! Wiedinmyer, C., Akagi, S. K., Yokelson, R. J., Emmons, L. K., Al-Saadi,
! J. A., Orlando, J. J., and Soja, A. J.: The Fire INventory from NCAR (FINN) 
! - a high resolution global model to estimate the emissions from open 
!  burning, Geosci. Model Dev. Discuss., 3, 2439-2476, 
!   doi:10.5194/gmdd-3-2439-2010, 2010.
! http://www.geosci-model-dev-discuss.net/3/2439/2010/gmdd-3-2439-2010.html
!
!  
! Wiedinmyer, C.; Yokelson, R. J. & Gullett, B. K. Global Emissions of
! Trace Gases, Particulate Matter, and Hazardous Air Pollutants from
! Open Burning of Domestic Waste Environmental Science & Technology,
! 2014, 48, 9523-9530
!
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
!  emepctm calls just call Fire_Emis(daynumber)
!  and put the day-testing code here. This lets the module decide if new
!  emissions are needed, and keeps all forest-fire logic here
!----------------------------------------------------------------
use CheckStop_mod,         only: CheckStop,CheckNC
use ChemDims_mod ,         only: NSPEC_TOT, NSPEC_SHL
use ChemSpecs_mod
use Config_module,         only: MasterProc, DataDir, KMAX_MID, USES, &
                                IOU_INST,BBMODE,BBverbose,persistence,&
                                fire_year,&
                                cmxBiomassBurning_FINN, &
                                cmxBiomassBurning_GFASv1, &
                                BBneed_file,BBneed_date,BBneed_poll,&
                                BBMAP,GFED_PATTERN,FINN_PATTERN,GFAS_PATTERN
use Debug_module,          only: DEBUG   ! -> DEBUG%FORESTFIRES
use GridValues_mod,        only: i_fdom, j_fdom, debug_li, debug_lj, &
                                debug_proc,xm2,GRIDWIDTH_M, A_bnd,B_bnd
use Io_mod,                only: datewrite, IO_NML, IO_TMP, open_file, ios
use Io_RunLog_mod,         only: PrintLog
use MetFields_mod,         only: z_mid, z_bnd, hmix
use netcdf,                only: nf90_open, nf90_nowrite, nf90_close
use NetCDF_mod,            only: ReadTimeCDF,ReadField_CDF,Out_netCDF,Real4,&
                                closedID
use NumberConstants,       only: UNDEF_I, UNDEF_R
use OwnDataTypes_mod,      only: Deriv, TXTLEN_SHORT, TXTLEN_FILE
use Par_mod,               only: LIMAX, LJMAX, me,limax,ljmax
use PhysicalConstants_mod, only: AVOG
use SmallUtils_mod,        only: find_index, key2str
! No. days per year, date-type:
use TimeDate_mod,          only: current_date,day_of_year,max_day
use TimeDate_ExtraUtil_mod,only: date2string,nctime2string,date2nctime,date2file
use ZchemData_mod,         only: rcemis
use EmisDef_mod,           only: Emis_CO_Profile   !FOR TESTING

implicit none
private

public  :: Fire_Emis
public  :: Fire_rcemis
public  :: burning

private :: Config_Fire
private :: make_mapping   ! chooses FINN or GFAS 
private :: readCMXmapping ! reads CMX_BiomassBurning_xxx.txt mapping files
private :: checkNewFFrecord
private :: Export_FireNc  ! not working?

logical, allocatable, dimension(:,:), save :: burning
real,  allocatable, dimension(:,:,:), save :: BiomassBurningEmis
logical, save :: monthlyEmis = .false.

integer, save :: ieCO=-1 ! index for CO
integer, save :: KEMISFIRE

character(len=TXTLEN_SHORT), save :: FF_poll = 'NOT_SET'
integer :: iemep

!/ Defintions of BB data. If known, we assign the BB pollutant which
!  corresponds to each possible EMEP emission file. 

! Assign mol. wts of the BB data  where known. If mol. wt set to
! zero, the code in Fire-rcemis will use the values from the
! ChemSpecs_mod, species()%molwt.
!
! If BB doesn't have emissionss, set a "-" for BB, then the
! desired emission factor (g/kg DW), and then follow the
! example in Fire_setups for NH3 (GFED):
! (above text if old. May need update)
!

! =======================================================================
!/ - 2019 update for GenChem -  CMX  system 
! list of possible species names from the FINN and GFAS inputs.
! Note: not all species are needed (some names change for different years and
! versions), but only species from these lists are allowed.
! 2022-12-15 added PM10 to cope with FINNv2.5. Still works with 1.5 data
character(len=TXTLEN_SHORT), dimension(24) :: &
  POSSIBLE_FINNv2p5_SPECS  = [ character(len=TXTLEN_SHORT):: &
   'CO', 'NO' ,'NO2' ,'SO2' ,'NH3' , &
   'ACET' ,'ALD2' ,'ALK4' ,'C2H6' ,'C3H8', &
   'CH2O', 'MEK','PRPE' ,'PM25' ,'PM10', 'OC' ,&
   'BC' ,'C2H4' ,'GLYC' ,'HAC' ,'BENZ' ,&
   'TOLU','XYLE' ,'MGLY' ]

character(len=TXTLEN_SHORT), dimension(22) :: &
  POSSIBLE_GFASv1_SPECS  = [ character(len=TXTLEN_SHORT):: &
   'cofire', 'ch4fire', 'h2fire', 'noxfire', 'pm2p5fire', &
   'tpmfire', 'ocfire', 'bcfire', 'so2fire', &
   'ch3ohfire', 'c2h5ohfire', 'c2h4fire', 'c3h6fire', 'c5h8fire', &
   'toluenefire', 'hialkenesfire', 'hialkanesfire', 'ch2ofire','c2h4ofire', &
   'nh3fire', 'c2h6fire', 'c4h10fire' ]

type, private :: cmxmap_t
  character(len=TXTLEN_SHORT) :: bbSpec = "NOTSET"
  real :: unitsfac = 1.0
  real :: frac     = 1.0
  character(len=TXTLEN_SHORT) :: emSpec = "NOTSET"
end type cmxmap_t
type(cmxmap_t), dimension(100), private, save :: cmxmapping = cmxmap_t()

!
type :: bbtype
  character(len=TXTLEN_SHORT) :: BBname = "NOTSET"
  real :: unitsfac = 1.0
  real :: frac    = 1.0
  integer :: emep = -1
end type bbtype

!----------------------------------------------

  character(len=200), public, save :: &
    BiomassBurningMapping
  integer, private, save:: &
    NBB_DEFS    ,& ! No mapping lines below
    NEMEPSPECS   ! No EMEP chemical mech specs used
  type(bbtype), private, allocatable :: FF_defs_BB(:)

!----------------------------------------------
! matrix to get from forest-fire species to EMEP ones

integer, save, allocatable :: emep_used(:)
real,    save, allocatable :: sum_emis(:)



! interpolation method in ReadField_CDF
character(len=30), save :: bbinterp = '-'

! Notes on interpolation choices: (from NetCDF_mod)
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

integer,   save :: nread=-1   ! records in forest fire file
logical, save ::    &
  debug_level(-1:5)=.false.
! =======================================================================
contains
subroutine Config_Fire()
  logical, save :: first_call=.true.
  logical :: dbg0  ! see also debug_level. Not harmonised
  integer :: ios, ne, n, k
  character(len=*), parameter :: dtxt='BB:Config'
  integer ::  N_LEVELS  ! =9 for <rv4.50 code and standard 20 model levels

  if(.not.first_call) return
  dbg0 = (DEBUG%FORESTFIRE.and.MasterProc)
  if(dbg0) write(*,*) dtxt//" selects ",BBMAP,',', BBMODE

  select case(BBMAP)
    case("GFED")
      persistence=8  ! 8-day records
      bbinterp = 'conservative'
    case("FINN")
      persistence=1  ! 1-day records
      bbinterp = 'mass_conservative'
      BiomassBurningMapping = "FINNv1.5"  ! ok for 2.5 also
   case("GFAS")
      persistence=3  ! 1-day records, valid for 3 day in forecast runs
      bbinterp = 'conservative'
      BiomassBurningMapping = "GFASv1"
   case default
      call CheckStop(dtxt//"Unknown Mapping")
  end select

  call make_mapping() !mapping definitions between forestfire species and emep species

  call PrintLog(dtxt//" Mapping: "//trim(BiomassBurningMapping),MasterProc)

  allocate(emep_used(NEMEPSPECS))
  emep_used = 0
  allocate(sum_emis(NEMEPSPECS))
  sum_emis = 0.0

  
  ! set vebosity levels
  BBverbose=min(max(0,BBverbose),4) ! debug verbosity 0,..,4
  debug_level(:0)=.false.
  debug_level(1)=DEBUG%FORESTFIRE.and.MasterProc.and.first_call
  debug_level(2)=DEBUG%FORESTFIRE.and.MasterProc
  debug_level(3)=DEBUG%FORESTFIRE
  debug_level(4:)=.true.

  allocate(BiomassBurningEmis(NEMEPSPECS,LIMAX,LJMAX),&
           burning(LIMAX,LJMAX),stat=ios)
  allocate(Emis_CO_Profile(LIMAX,LJMAX,KMAX_MID))      !TESTING

  call CheckStop(ios,dtxt//" BiomassBurningEmis alloc problem")
  ne = 0     ! number-index of emep species

  do n=1, NBB_DEFS              ! Only unique EMEP SPECS in emep_used
    iemep = FF_defs_BB(n)%emep 
    if(find_index(iemep,emep_used(:))>0) cycle
    
    ne = ne + 1
    emep_used(ne) = iemep

    if( iemep <1) print *, 'ABB', n, FF_defs_BB(n)%emep, iemep

    ! CO is special. Keep the index
    !TRACER if(species(iemep)%name=="CO") ieCO=ne
    ! Now allow emep species to be FFIRE_CO
    if(ieCO<0 .and. species(iemep)%name=="CO" ) ieCO=ne
    if(ieCO<0 .and.FF_defs_BB(n)%BBname=="CO"    ) ieCO=ne

    if(dbg0) write(*,"(i4,1x,a,i4,2x,a)")ne, &
      dtxt//" Mapping "//trim(FF_defs_BB(n)%BBName)//" onto " &
      ,iemep, trim(species(iemep)%name)
  end do !n
  call CheckStop(ieCO<1,&
     dtxt//"No mapping for 'CO' found on "//BiomassBurningMapping)
  if (any(emep_used<1) ) then
    do n = 1, size(emep_used)
       print *, "ABBn ", n, emep_used(n)
    end do
  end if
  call CheckStop(any(emep_used<1),&
     dtxt//"UNSET FFIRE EMEP "//BiomassBurningMapping)

!P800 method is original version, used to rv4.52(??) = simple release height
!defintion which is model levels independent. The highest level is the highest
! level boundary below 800 hPa = ca 2000 m (standard atmosphere)
! KEMISFIRE from this method is overwritten later if
! USES%FFireDispMethod == 'PBL' , in the config_emep.nml file 
  call CheckStop( .not. ( USES%FFireDispMethod == 'P800' .or. &
                          USES%FFireDispMethod == 'PBL' ),    &
        dtxt//"UNSET FFIRE Hmix Method "//trim(USES%FFireDispMethod) )

  do k = KMAX_MID+1,2,-1
     if ( dbg0 ) write(*,*) dtxt//'Kcheck0', k
     if(A_bnd(k)+101325.0*B_bnd(k)< 80000.0)exit
  enddo
  KEMISFIRE = k
  N_LEVELS = KMAX_MID - KEMISFIRE + 1 

  if(MasterProc .and. first_call)write(*,fmt='(A,I3,A)')&
    dtxt//'Orig (P800) Method => ', N_LEVELS, ' lowest levels'

  first_call=.false.
end subroutine Config_Fire

subroutine make_mapping()
  ! define the mapping between names and species defined in the forest fire
  ! input file, and emep species. Must be compatible for different chemistry
  ! mechanisms (CMs), i.e. all the emep species are not necessarily defined.
  
  integer :: emep_used(NSPEC_TOT)
  integer :: ncmx_defs,ncmx_emep
  type(bbtype), dimension(100) :: tmpFF_defs
  character(len=*), parameter :: dtxt='BB:makemap'
  
  integer i,n

  if(BBMAP=='FINN') &
     call readCMXmapping(cmxBiomassBurning_FINN, POSSIBLE_FINNv2p5_SPECS,&
        NBB_DEFS,NEMEPSPECS,tmpFF_defs)

  if(BBMAP=='GFAS') &
     call readCMXmapping(cmxBiomassBurning_GFASv1, POSSIBLE_GFASv1_SPECS,&
        NBB_DEFS,NEMEPSPECS,tmpFF_defs)

  if(.not.allocated(FF_defs_BB))allocate(FF_defs_BB(NBB_DEFS))

  FF_defs_BB = tmpFF_defs(1:NBB_DEFS) !CMX FF_defs


  if(MasterProc) then
     write(*,fmt='(A,I5,A,I5,A)') dtxt//' will read ',NBB_DEFS,&
        ' species, mapped into ',NEMEPSPECS,' emep species'
  end if
  
  do i = 1, NBB_DEFS
     if(masterproc)write(*,*)i,dtxt//' '//trim(FF_defs_BB(i)%BBName)&
           //' maps to '//species(FF_defs_BB(i)%emep)%name
  enddo

end subroutine make_mapping

subroutine Fire_Emis(daynumber)
!.....................................................................
!**    DESCRIPTION:
!    Reads forest-fire emissions. So far set up for GFED 8d, but in
!    principal we can re-code by simply adding alternative
!    subroutines, e.g. to cope with pre-2001 monthly emissions

  integer, intent(in) :: daynumber   ! daynumber (1,..... ,365)


  real :: rdemis(LIMAX,LJMAX)  ! Emissions read from file
  integer :: i,j,nstart, alloc_err, iBB, n
  logical, save :: first_call = .true.
  logical       :: was_first_call
  integer, save  :: nn_old=-1
  real    :: fac, to_kgm2s   

  integer :: ind, ncFileID=closedID
  integer :: loc_maxemis(2) ! debug

  character(len=TXTLEN_FILE), save :: fname='new'
  logical :: debug_me=.false., debug_ff=.false.,debug_nc=.false.,&
              newFFrecord=.false., BBfound
  real :: xrdemis(LIMAX,LJMAX) ! MODE=*_AVG
  integer :: dn1, dn2, ndn          ! MODE=*_AVG
  integer :: yyyy, mm, dd, hh
  character(len=*), parameter :: dtxt='BB:Fire_Emis:'

  ! copy current flag:
  was_first_call = first_call
  ! first?
  if( first_call ) then
    ! allocate arrays etc:
    call Config_Fire()
    ! reset flag:
    first_call = .false.
  end if ! first
  
  debug_me=DEBUG%FORESTFIRE .and. debug_proc
  if ( debug_me ) write(*,*)  'FFIREDEBUG!!', me, i_fdom(debug_li), j_fdom(debug_lj)
  debug_ff=debug_level(BBverbose)
  debug_nc=debug_level(BBverbose-1)

  yyyy = current_date%year
  mm   = current_date%month 
  dd   = current_date%day
  hh   = current_date%hour
  if(fire_year>1) yyyy=fire_year
  if(debug_me)then
    write(*,'(a,7i5)') dtxt//"current date and time -"//trim(BBMODE),&
      yyyy,mm,dd,fire_year, nn_old, mm
    if(allocated(BiomassBurningEmis)) &
      write(*,'(a,es12.3)') dtxt//'sum', sum(BiomassBurningEmis(1,:,:))
  end if
  select case(BBMODE)
    case("HOURLY_REC","H","h")
      if(nn_old==hh) return         ! Calculate once per hour
      nn_old=hh
      dn1=hh
      dn2=hh
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
    case("MONTHLY_CLIM", "C", "c")
      if(nn_old==mm) return         ! Calculate once per month
      nn_old=mm
      dn1=day_of_year(yyyy,mm,01)
      dn2=dn1
      nstart=mm ! Only 12 records per file
    case("YEARLY_AVG","Y","y")
      if(nn_old==yyyy) return       ! Calculate once per year
      nn_old=yyyy
      dn1=day_of_year(yyyy,01,01)
      dn2=day_of_year(yyyy,12,31)
    case default
      call CheckStop("Unknown ForestFire MODE="//trim(BBMODE))
  end select

  if(debug_me) write(*,'(a,5i5)') dtxt// "newFFrec checks ",&
      yyyy,mm,dd, dn1, dn2
  if(dn1<dn2)then
    nstart=daynumber    ! for debug info
    if(debug_me) write(*,'(a,5i5)') dtxt// "FFalloc nstart= ", nstart
  else
    ! newFFrecord: has pollutant|fname|record changed since last call?
    call checkNewFFrecord([yyyy,mm,dd,hh], ncFileID, fname, newFFrecord, nstart)
    if(debug_me) write(*,'(a,L2,5i5)') dtxt// "FFchecked nstart= ",newFFrecord, nstart
    FF_poll=""
    if(.not.newFFrecord) then 
       if(debug_me) write(*,'(a,5i5)') dtxt//" newFFrec not set ",yyyy,mm,dd,hh
       return                        ! Continue if new record||file  
    end if
  end if

  if(debug_me) then
    write(*,'(a,5i5)') dtxt// "newFFrec WAS set ", yyyy,mm,dd, dn1, dn2
    write(*,*) dtxt//'Starting MODE=',trim(BBMODE),&
      date2string(" YYYY-MM-DD",[yyyy,mm,dd]),was_first_call,debug_ff,debug_nc
    write(*,*) dtxt//' Interp= ', trim(bbinterp), dn1, dn2, nstart
  end if

  BiomassBurningEmis(:,:,:) = 0.0
  Emis_CO_Profile(:,:,:) = 0.0
  
  ! We need to look for forest-fire emissions which are equivalent
  ! to the standard emission files:
  do iBB = 1, NBB_DEFS
    FF_poll = FF_defs_BB(iBB)%BBname
    iemep   = FF_defs_BB(iBB)%emep  ! 
    ind     = find_index( iemep, emep_used ) !Finds 1st emep in BiomassBurning

    if(debug_me)then
      write(*,"( a,3i5, a8,i3)") dtxt//" SETUP: ", iBB,iemep,ind, &
        trim(FF_poll), len_trim(FF_poll)
      !DS if(debug_ff) &
        write(*,*) dtxt//'BBMAP ',BBMAP,':',me,iBB,nstart, monthlyEmis, &
           ncFileID, trim(FF_poll),trim(fname)
    end if

   ! if(.not.need_file|time|poll) continue if file|time|poll is not found
    rdemis(:,:)=0.0

    !---------  read data:

    if(dn1<dn2)then
        rdemis = 0.0
        xrdemis = 0.0
        ndn=0
        do dd = dn1, dn2
          call checkNewFFrecord([yyyy,mm,dd,00], &
                ncFileID, fname, newFFrecord, nstart)
          if(newFFrecord) then
            call ReadField_CDF(fname,FF_poll,xrdemis,nstart,interpol=bbinterp,&
             needed=BBneed_poll,found=BBfound,UnDef=0.0,debug_flag=debug_nc,&
             ncFileID_given=ncFileID)
          end if
          rdemis = rdemis + xrdemis                 ! month total
          ndn    = ndn + 1
        end do
    else
        ndn=1
        call ReadField_CDF(fname,FF_poll,rdemis,nstart,interpol=bbinterp,&
          needed=BBneed_poll,found=BBfound,UnDef=0.0,debug_flag=debug_nc,&
          ncFileID_given=ncFileID)
        if(debug_me.and.FF_poll=="CO" ) write(*,"(a,i5,a,es12.3)") &
           dtxt//" CO READ: ", nstart, trim(FF_poll), maxval(rdemis)
    end if
    if(debug_me.and.FF_poll=="CO" ) write(*,"(a,i5,a,es12.3)") &
           dtxt//" FIRE DEBUGCO READ: ", nstart, me, rdemis(debug_li,debug_lj)
    !-------- 

    if ( BBfound ) then
       select case(BBMAP)
       case("GFED")
   
        !unit conversion, GFED [g/m2/8day]->[kg/m2/s]
         to_kgm2s = 1.0e-3 /(8*24.0*3600.0)
         if(ndn>1) to_kgm2s=to_kgm2s/ndn               ! total-->avg.
         forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*to_kgm2s
   
       case("FINN")
   
        ! unit conversion, FINN [mole/day]->[kg/m2/s]
        ! (Can be negative if REMPPM to be calculated)
         fac=FF_defs_BB(iBB)%unitsfac * FF_defs_BB(iBB)%frac ! --> [kg/day]
         fac=fac/(GRIDWIDTH_M*GRIDWIDTH_M*24.0*3600.0) ! [kg/day]->[kg/m2/s]
         if(ndn>1) fac=fac/ndn                         ! total-->avg.
         forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*fac*xm2(i,j)
   
       case("GFAS")
   
        ! GFAS units are [kg/m2/s]. No further unit conversion is needed.
        ! However, fac can be /=1, e.g. when REMPPM is calculated
         fac=FF_defs_BB(iBB)%unitsfac * FF_defs_BB(iBB)%frac
         if(ndn>1) fac=fac/ndn                         ! total-->avg.
         if(fac/=1.0) forall(j=1:ljmax,i=1:limax) rdemis(i,j)=rdemis(i,j)*fac
       end select
   
      ! Assign . units should be [kg/m2/s] here 
       forall(j=1:ljmax,i=1:limax) &
         BiomassBurningEmis(ind,i,j) = BiomassBurningEmis(ind,i,j) + rdemis(i,j) 
   
       if(debug_me) &
          write(*,"(3a10,2i4,f8.3,es12.3)") dtxt//" SUMS:", &
           trim(FF_poll), trim( species(iemep)%name), me, ind, &
           species(iemep)%molwt, sum( BiomassBurningEmis(ind,:,:) )
   
       call PrintLog(dtxt//":: Assigns "//trim(FF_poll),&
         was_first_call.and.MasterProc)
   
       if(debug_me) sum_emis(ind)=sum_emis(ind)+&
             sum(BiomassBurningEmis(ind,:,:))
    else ! BBfound false
       call PrintLog(dtxt//":: Skips   "//trim(FF_poll),&
         was_first_call.and.MasterProc)
    end if ! BBfound
  end do ! BB_DEFS

  ! have to close the file here
  call CheckNC(nf90_close(ncFileID),dtxt//"close:"//trim(fname))
  ncFileID=closedID

  ! For cases where REMPPM25 s derived as the difference between PM25 and
  !  (BC+1.7*OC) we need some safety:

  BiomassBurningEmis(:,:,:) = max( BiomassBurningEmis(:,:,:), 0.0 )

  ! Logical to tell if there is any emission here to worry about
  burning(:,:) = ( BiomassBurningEmis(ieCO,:,:) > 1.0e-19 )



  ! Some databases (e.g. FINN, GFED) have both total PM25 and EC, OC. The
  ! difference, REMPPM25, is created by the BiomasBurning mapping procedure,
  ! but we just check here
  if(debug_me) then
    n = ieCO
    loc_maxemis = maxloc(BiomassBurningEmis(n,:,: ) )

    associate ( idbg=>loc_maxemis(1), jdbg=>loc_maxemis(2) )

    write(*,*) 'FORESTFIRE DEBUGBURN', burning(debug_li,debug_lj), BiomassBurningEmis(ieCO,debug_li,debug_lj) 
    write(*,"(a,i4,i3,2i4,2i5,es12.3, 2i4)") dtxt//"SUM_FF CHECK ME: ", &
       daynumber, me, loc_maxemis, i_fdom(idbg), j_fdom(jdbg),&
         BiomassBurningEmis(n,idbg,jdbg), debug_li,debug_lj

    call datewrite(dtxt//"SUM_FF CHECK CO: ",  &
      (/ daynumber, n, i_fdom( idbg ), j_fdom( jdbg ) /) ,&
      (/  sum_emis(n), maxval(BiomassBurningEmis(n,:,: ) ), &
          BiomassBurningEmis(n,debug_li,debug_lj) /) )
    end associate ! idbg, jdbg
  end if ! debug_me
  !end associate ACDATES

end subroutine Fire_Emis

subroutine checkNewFFrecord(ymdh, ncFileID,fname,new,nstart)
  integer, intent(in) :: ymdh(4)
  integer, intent(inout) :: ncFileID
  character(len=*), intent(inout) :: fname
  logical, intent(inout) :: new
  integer, intent(out) :: nstart

  character(len=TXTLEN_SHORT), save      :: poll_old=''
  character(len=TXTLEN_FILE), save :: file_old=''
  integer, save                          :: record_old=-1
  real, dimension(745), save :: fdays=-1 ! up to a month of hourly recs
  logical :: fexist=.false.
  real :: ncday(0:1)
  character(len=*),parameter:: dtxt='BB:newFFrecord:'
  logical :: debug_me

  debug_me=DEBUG%FORESTFIRE .and. debug_proc

  ! Check: New file
  select case(BBMAP)
!note: without "mode=YMDH", there might be problems with names containing 'ss'
    case("GFED");fname=date2file(GFED_PATTERN,ymdh,persistence-1,"days",mode='YMDH')
    case("FINN");fname=date2file(FINN_PATTERN,ymdh,persistence-1,"days",mode='YMDH')
    case("GFAS");fname=date2file(GFAS_PATTERN,ymdh,persistence-1,"days",mode='YMDH')
  end select
  fname=key2str(fname,'DataDir',DataDir) ! expand DataDir keysword

  if(debug_me) then
    write(*,*)  dtxt//trim(fname), me, ymdh ! TMP
    write(*,*)  dtxt//" Old:", trim(file_old)
    write(*,*)  dtxt//" IDs ", ncFileID, closedID ! TMP
  end if 
  if(fname/=file_old)then
    if(DEBUG%FORESTFIRE.and.MasterProc)then 
       write(*,*)dtxt//" new file:.. ",trim(fname(36:))
       write(*,*)dtxt//' new fpoll FF_poll:',trim(FF_poll)
    endif
  ! close old ncFile, if already open
    if(ncFileID/=closedID)&
      call CheckNC(nf90_close(ncFileID),dtxt//"close:"//trim(file_old))
    ncFileID=closedID
  ! check if new file exists
    inquire(file=fname,exist=fexist)    ! check if fname exixts
    if(.not.fexist)then
      if(MasterProc)then
        !print *, dtxt//" file not found: ",trim(fname(36:))
        call CheckStop(BBneed_file,dtxt//"Missing file:"//trim(fname))
      end if
      burning(:,:) = .false.
      new=.false.
      return
    end if
  ! read all times records in fname, and process them
    nread=-1
    fdays(:)=-1.0                       
    call ReadTimeCDF(fname,fdays,nread) 
    if ( nread == 12 ) monthlyEmis = .true.
    if(debug_me) write( *,'(a,3f8.1,a,2i6,L2)')  dtxt//" fdays ", &
        fdays(1:3), '... n=', count(fdays>0), nread, monthlyEmis
    record_old=-1                       
  end if

  ! Check: New pollutant
  if(FF_poll/=poll_old)then
    if(DEBUG%FORESTFIRE.and.MasterProc) &
      write(*,*)dtxt//" new pollutant: ",trim(FF_poll)
  end if

  ! Check: New time record
  call date2nctime(ymdh,ncday(1))
  ncday(0)=ncday(1)-persistence+1
  nstart=MAXLOC(fdays(:nread),DIM=1,&
    MASK=(fdays(:nread)>=ncday(0)).and.(fdays(:nread)<(ncday(1)+1.0)))
  if ( monthlyEmis ) nstart = ymdh(2) !AUG
  if(nstart/=record_old)then
    if(DEBUG%FORESTFIRE.and.MasterProc) then
      write(*,'(a,2f8.1,3i5,f8.1)') dtxt//" ncday???    ",&
          ncday(0), ncday(1), persistence, nstart,record_old, fdays(nstart)
      write(*,'(a,2i6,a,f8.1,i4)') dtxt//" new record: ", nstart,&
        record_old, nctime2string("(YYYY-MM-DD hh:mm)",fdays(nstart)), &
         fdays(nstart), persistence
    end if

    if ( .not. monthlyEmis ) then
       if((fdays(nstart)<ncday(0)).or.(fdays(nstart)>=(ncday(1)+1.0)))then
         if(MasterProc)then
           write(*,*)dtxt//" no records between ",&
             nctime2string("YYYY-MM-DD 00:00",ncday(0))," and ",&
             nctime2string("YYYY-MM-DD 23:59",ncday(1))
           call CheckStop(BBneed_date,dtxt//"Missing records")
         end if
         burning(:,:) = .false.
         new=.false.
         return
       end if
    end if !  not monthlyEmis
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
end subroutine checkNewFFrecord
!=============================================================================

subroutine Fire_rcemis(i,j)
!  Disperses the fire emissions vertically and converts to molecules/cm3/s.

!// Injection height: here over 8 levels. Alternative could be PBL
!   or  equally up to ca. 2*PBL (suggested by Sofiev, GEMS)
!   QUERY - should the emissions be divided equally by level?
!   - will give a higher mixing ratio for thinner levels

  integer, intent(in) :: i,j

  real, dimension(KMAX_MID) :: invDeltaZfac !  height of layer in m div 9
  integer ::  k, n, iem

  integer ::  N_LEVELS  ! = 9 for standard 20 model levels

  character(len=*), parameter :: dtxt = 'BB:rcemis'
  real    :: origrc, fac, dP, P0
  logical :: debug_flag

  P0 = 101325.0
  debug_flag = (DEBUG%FORESTFIRE.and.debug_proc .and.&
                i==debug_li.and.j==debug_lj)
  if ( debug_flag ) write(*,*)  'FFIREDEBUG-RC!!', me, i_fdom(i), j_fdom(j), debug_flag, BiomassBurningEmis(ieCO,i,j)
  if(debug_flag.and.BiomassBurningEmis(ieCO,i,j) > 1.0e-10)  &
  !if(DEBUG%FORESTFIRE.and.debug_proc .and.BiomassBurningEmis(ieCO,i,j) > 1.0e-10)  &
  !if(DEBUG%FORESTFIRE .and.BiomassBurningEmis(ieCO,i,j) > 1.0e-10)  &
    write(*,"(a,L1,5i4,es12.3,f9.3)") dtxt//"DEBUGRC",debug_flag, me, i,j, &
      i_fdom(i), j_fdom(j), BiomassBurningEmis(ieCO,i,j)

!JEJ 3/22
  if(USES%FFireDispMethod == 'PBL' )then
    do k = KMAX_MID,2,-1
       if ( debug_flag )  write(*,*) dtxt//'VERT:', k, z_bnd(i,j,k-1), hmix(i,j,1)
       if(z_bnd(i,j,k-1) >  hmix(i,j,1)) exit
!!     if(z_mid(i,j,k-1) >  hmix(i,j,1)) exit
    enddo
    KEMISFIRE = k
  end if 

  N_LEVELS = KMAX_MID - KEMISFIRE + 1 

  if ( debug_flag) then ! MasterProc .and. i == 2 and. j==2 ) then
    write(6,"(a,2f8.1,i3,a)") dtxt//trim(USES%FFireDispMethod)//' => New H ',&
      hmix(i,j,1), z_bnd(i,j,KEMISFIRE), N_LEVELS, ' lowest levels'
  end if


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

  !  We distribute the emissions evenly, i.e. proportionnally to the levels thickness.
  !  The "thickness" is measured in a pressure scale (also the partial pressure of pollutants 
  !  decreases with height). For simplicity we do not account for difference in surface 
  !  pressure for the distribution.
  
  do k = KEMISFIRE, KMAX_MID
     invDeltaZfac(k) = 1.0/ (z_bnd(i,j,k) - z_bnd(i,j,k+1)) 
  end do
 
  do n = 1, NEMEPSPECS 
    iem = emep_used(n)
    origrc = rcemis(iem,KMAX_MID)   ! just for printout
    fac =  0.001 * AVOG /species(iem)%molwt    ! MW scale if needed

    ! distribute vertically:
    !dP=total "thickness"
    dP = A_bnd(KMAX_MID+1)+P0*B_bnd(KMAX_MID+1) - (A_bnd(KEMISFIRE)+P0*B_bnd(KEMISFIRE))    
    do k = KEMISFIRE, KMAX_MID
      rcemis(iem,k) = rcemis(iem,k) + BiomassBurningEmis(n,i,j)*invDeltaZfac(k)*fac&
           *(A_bnd(k+1)+P0*B_bnd(k+1) - (A_bnd(k)+P0*B_bnd(k)))/dP!scale with layer thickness

     ! nb :::: ONLY FOR CO!
      Emis_CO_Profile(i,j,k) = BiomassBurningEmis(n,i,j)*invDeltaZfac(k)*fac&
           *(A_bnd(k+1)+P0*B_bnd(k+1) - (A_bnd(k)+P0*B_bnd(k)))/dP
    end do !k

    if(debug_flag) then
      k=KMAX_MID
      write(*,"(a,2i3,1x,a8,i4,es10.2,4es10.2)") dtxt//"FIRERC ",&
        n, iem, trim(species(iem)%name), k, BiomassBurningEmis(n,i,j),&
        invDeltaZfac(k), origrc, rcemis(iem,k)
    end if

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


  ! SUBROUTINE to read CMX_BoundaryCondition file
  !   
  subroutine readCMXmapping(fname,bbspecs,ndefs,nemep,tmpFF_defs) 
     character(len=*), intent(in) :: fname
     character(len=*),dimension(:), intent(in) :: bbspecs
     integer, intent(out) :: ndefs,nemep
     type(bbtype), dimension(:), intent(inout) :: tmpFF_defs
     character(len=200) :: txtinput
     character(len=*), parameter :: dtxt='CMXbb:'
     integer :: ind1, ind2,  npossible=0
     integer, dimension(NSPEC_TOT) ::emep_used = 0
     type(cmxmap_t) :: bbcmx
     ndefs = 0

     call PrintLog(dtxt//' START '//fname,MasterProc)
    ! check file exists. open_file returns ios for iostat
     call open_file(IO_TMP,'r',fname,needed=.true.) ! returns ios
     call CheckStop(ios,dtxt//"open_file error on " // fname )

     do while(.true.)
       read(IO_TMP,fmt='(a100)',iostat=ios) txtinput
       if ( ios /= 0 ) exit   ! likely end of file
       if (txtinput(1:1) == '#' ) then
          call PrintLog(dtxt//txtinput,MasterProc)
          cycle
       end if
       npossible = npossible + 1
       read(txtinput,*) bbcmx
       if(MasterProc) print *, 'CMXbbTEST ', trim(bbcmx%bbSpec)//';'// trim(bbcmx%emSpec)

      ! Find and set indices
       ind1 = find_index(bbcmx%bbSpec,bbspecs(:),any_case=.true.)
       if (ind1<1) print *, 'PANIClin ', trim(txtinput)
       if (ind1<1) print *, 'PANIC ', ind1, trim(bbcmx%bbSpec)
       call CheckStop(ind1<1,dtxt//'ERROR bbSpec not found:'//bbcmx%bbSpec)

       ind2 = find_index(bbcmx%emSpec,species(:)%name,any_case=.true.)
       call CheckStop(ind2<1,dtxt//'ERROR emep spec not found:'//bbcmx%emSpec)
       emep_used(ind2) = 1

      ! bbtype: BBname unitsfac frac emep
       ndefs = ndefs + 1
       tmpFF_defs(ndefs)%BBname    = bbcmx%bbSpec
       tmpFF_defs(ndefs)%unitsfac  = bbcmx%unitsfac
       tmpFF_defs(ndefs)%frac      = bbcmx%frac
       tmpFF_defs(ndefs)%emep      = ind2
!       print *, 'CMXERR', ind2, trim(txtinput)//';'// trim(species(ind2)%name)
       write(txtinput,'(i2, a, i4,a)') ndefs, trim(txtinput)// ' => ', &
            ind2, trim(species(ind2)%name)
       call PrintLog(dtxt//'SET:'//txtinput,MasterProc)
       
     end do
     close(IO_TMP)
     nemep = sum( emep_used )
     write(txtinput,'(a, 3i4)') ' nposs, ndefs,nemep', npossible, ndefs, nemep
     call PrintLog(dtxt//' DONE '//txtinput,MasterProc)
     
  end subroutine readCMXmapping

endmodule ForestFire_mod
!=============================================================================
