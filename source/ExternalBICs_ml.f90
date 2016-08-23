! <ExternalBICs_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_5(2809)>
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
module ExternalBICs_ml
! External Boundary and Initial Conditions set according to
!  ExternalBICs_config namelist
! Nothing in this file will be used if
!  EXTERNAL_BIC_SET = .false. in ExternalBICs_config namelist
use CheckStop_ml,           only: CheckStop
!CMR use ChemChemicals_ml,       only: species_adv
!CMR use ChemSpecs_adv_ml,       only: NSPEC_ADV
use ChemSpecs,              only: NSPEC_ADV, species_adv
use Io_ml,                  only: PrintLog,IO_NML
use ModelConstants_ml,      only: MasterProc, DEBUG=>DEBUG_NEST_ICBC
use SmallUtils_ml,          only: find_index
use TimeDate_ExtraUtil_ml,  only: date2string
implicit none

private
public :: set_extbic

logical, public, save :: &
  USE_EXTERNAL_BIC = .false., & ! use external (non Unimod) BCs
  EXTERNAL_BIC_SET = .false., & ! external BC description/setup has been found
  TOP_BC           = .false.    ! BCs include top level

integer,private, parameter :: &
  BIC_NAME_LEN=16

character(len=BIC_NAME_LEN),public, save :: &
  EXTERNAL_BIC_NAME    = "DUMMY",           &
  EXTERNAL_BIC_VERSION = "use_any"

integer,save, public :: &
  iw=-1, ie=-1, js=-1, jn=-1, kt=-1 ! i West/East bnd; j North/South bnd; k Top

character(len=100),public, save :: &
  filename_eta     = 'EMEP_IN_BC_eta.zaxis'

character(len=*),public, parameter :: &
  ICBC_FMT="(A24,'=',A24,'*',F7.2,2L2,'=',I3)"
type, public :: icbc                ! Inital (IC) & Boundary Conditions (BC)
  character(len=24) :: spcname="none",varname="none" ! Unimod,BC_file names
  real              :: frac=1.0                      ! fraction to unimod variable
  logical           :: wanted=.false.,found=.false.  ! BC is wanted,found in file
  integer           :: ixadv=-1                      ! adv index, set from %spcname
endtype icbc

type, private :: icbc_desc          ! IC/BC description
  character(len=BIC_NAME_LEN) :: name="none",version="none"
  integer                     :: mapsize=-1
endtype icbc_desc

type(icbc), dimension(:), public, pointer :: &
  EXTERNAL_BC=>null() ! external (non Unimod) BCs detailed description/setup
type(icbc), dimension(NSPEC_ADV), private, target, save :: &
  map_bc              ! detailed description/setup from ExternalBICs_bc namelist

character(len=*),private, parameter :: &
  DEBUG_FMT="(A,' DEBUG: ',A,' ''',A,'''.')"

contains
subroutine Config_ExternalBICs()
!----------------------------------------------------------------------------!
! Read basic configuration for external (non Unimod) BCs.
! ICs are assumed to come from Unimod (Nest_ml.init_icbc).
!
! USE_EXTERNAL_BIC  Use of external BCs  
!        otherwise  Assume Unimod BCs (.not.EXTERNAL_BIC_SET)
! EXTERNAL_BIC_NAME description%name to look for on ExternalBICs_bc namelist
!----------------------------------------------------------------------------!
  integer :: ios
  logical, save     :: first_call=.true.
  NAMELIST /ExternalBICs_config/ &
    USE_EXTERNAL_BIC,EXTERNAL_BIC_NAME,EXTERNAL_BIC_VERSION,TOP_BC,filename_eta

  if(.not.first_call) return
  rewind(IO_NML)
  read(IO_NML,NML=ExternalBICs_config,iostat=ios)
  call CheckStop(ios,"NML=ExternalBICs_config")  
  if(DEBUG.and.MasterProc)then
    write(*,*) "NAMELIST IS "
    write(*,NML=ExternalBICs_config)
  endif
endsubroutine Config_ExternalBICs

subroutine set_extbic(idate)
!----------------------------------------------------------------------------!
! Read external (non Unimod) BCs description/setup.
! ICs are assumed to come from Unimod (Nest_ml.init_icbc).
!
! EXTERNAL_BIC_SET  BC description/setup has been found
!        otherwise  Assume Unimod BCs (Nest_ml.init_icbc)
! description       BCs basic info %name,%version,%mapsize
! map_bc            BCs detailed setup from ExternalBICs_bc namelist
! EXTERNAL_BC       pointer to the records on map_bc with data (%mapsize)
!----------------------------------------------------------------------------!
  integer,intent(in) :: idate(4)
  logical, save     :: first_call=.true.
  integer           :: ydmh=0,ios=0,n=0
  type(icbc_desc) :: description
  NAMELIST /ExternalBICs_bc/description,map_bc

  if(.not.first_call) return
  call Config_ExternalBICs()

  if(.not.USE_EXTERNAL_BIC)then
    EXTERNAL_BIC_SET=.false.
    call PrintLog("No external BICs set",MasterProc)
    first_call = .false.
    return
  endif

!--- Set BC type from idate: on first call only
  if(EXTERNAL_BIC_SET) return

!--- Determine %version to look for
  select case (EXTERNAL_BIC_VERSION)
  case("MACC_ENS","FORECAST")
    ydmh=dot_product(idate,nint((/1e6,1e4,1e2,1e0/)))
    select case (ydmh)
    case(:2011113023)            ! untill 2011-11-30 23:00      
      EXTERNAL_BIC_VERSION='IFS_MOZ_f7kn'
    case(2011120100:2014011423)  ! 2011-12-01 00:00 to 2014-01-14 23:00
      EXTERNAL_BIC_VERSION='IFS_MOZ_fkya'
    case(2014011500:)            ! from 2014-01-15 00:00 (available from 2012-09-04)
!   case(2014030100:)            ! from 2014-03-01 00:00 (available from 2012-09-04)
      EXTERNAL_BIC_VERSION='IFS_MOZ_fnyp'
    endselect
  case("MACC_EVA","IFS_MACC_EVA")  ! EVA2010/2011: GRG & AER    
    EXTERNAL_BIC_VERSION='IFS_MACC_EVA'
  case default
    EXTERNAL_BIC_VERSION='use_any'
  endselect

!--- Look for a ExternalBICs_bc with the correct %name and %version
  rewind(IO_NML)
  READ_NML: do
    read(IO_NML,NML=ExternalBICs_bc,iostat=ios)
    if(ios/=0) exit READ_NML
    if(DEBUG.and.MasterProc) write(*,DEBUG_FMT) "set_extbic","read_nml bc",&
      trim(description%name)//"/"//trim(description%version)
    EXTERNAL_BIC_SET=&
      EXTERNAL_BIC_NAME==description%name.and.&
     (EXTERNAL_BIC_VERSION=='use_any'.or.EXTERNAL_BIC_VERSION==description%version).and.&
      description%mapsize>0
    if(EXTERNAL_BIC_SET)then
      EXTERNAL_BC=>map_bc(1:description%mapsize)
      exit READ_NML
    endif
  enddo READ_NML
  if(.not.EXTERNAL_BIC_SET)then
    call PrintLog("No external BICs found",MasterProc)
    USE_EXTERNAL_BIC=.false.
    return
  endif
  if(DEBUG.and.MasterProc) write(*,DEBUG_FMT) "set_extbic", &
    date2string("BCs for YYYY-MM-DD hh type",idate),&
    trim(EXTERNAL_BIC_NAME)//"/"//trim(EXTERNAL_BIC_VERSION)

  do n = 1,size(EXTERNAL_BC%ixadv)
    EXTERNAL_BC(n)%ixadv=find_index(EXTERNAL_BC(n)%spcname,species_adv(:)%name)
    if(EXTERNAL_BC(n)%ixadv<1)then
      EXTERNAL_BC(n)%wanted=.false.
      if(MasterProc) write(*,DEBUG_FMT) "set_extbic","unknow variable",&
        trim(EXTERNAL_BC(n)%spcname)
    endif
  enddo

  if(MasterProc) &
    call PrintLog("External BICs set for "//EXTERNAL_BIC_NAME)
  EXTERNAL_BIC_SET = .true.
  first_call = .false.
endsubroutine set_extbic

endmodule ExternalBICs_ml
