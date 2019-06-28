! <ExternalBICs_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module ExternalBICs_mod
! External Boundary and Initial Conditions set according to
!  ExternalBICs_config namelist
! Nothing in this file will be used if
!  EXTERNAL_BIC_SET = .false. in ExternalBICs_config namelist
use CheckStop_mod,          only: CheckStop
use ChemDims_mod,           only: NSPEC_ADV
use ChemSpecs_mod,          only: species_adv
use Config_module,          only: MasterProc, BC_DAYS,&
     USE_EXTERNAL_BIC,EXTERNAL_BIC_NAME,EXTERNAL_BIC_VERSION,TOP_BC,filename_eta
use Debug_module,           only: DEBUG=>DEBUG_NEST_ICBC
use Io_mod,                 only: PrintLog,IO_NML
use OwnDataTypes_mod,       only: TXTLEN_SHORT
use SmallUtils_mod,         only: find_index
use TimeDate_mod,           only: date
use TimeDate_ExtraUtil_mod, only: date2string
implicit none

private
public :: set_extbic

interface set_extbic
  module procedure set_extbic_id
  module procedure set_extbic_cd
end interface set_extbic

integer,save, public :: &
  iw=-1, ie=-1, js=-1, jn=-1, kt=-1 ! i West/East bnd; j North/South bnd; k Top
!  BC_DAYS=0   ! #days in one BC file, for use old BCs in a FORECAST
              ! 0 means "do not look for old files"

character(len=*),public, parameter :: &
  ICBC_FMT="(A24,'=',A24,'*',F7.2,2L2,'=',I4)"
type, public :: icbc                ! Inital (IC) & Boundary Conditions (BC)
  character(len=24) :: spcname="none",varname="none" ! emepctm,BC_file names
  real              :: frac=1.0                      ! fraction to unimod variable
  logical           :: wanted=.false.,found=.false.  ! BC is wanted,found in file
  integer           :: ixadv=-1                      ! adv index, set from %spcname
end type icbc

type, private :: icbc_desc          ! IC/BC description
  character(len=TXTLEN_SHORT) :: name="none",version="none"
  integer                     :: mapsize=-1
end type icbc_desc

logical, public, save :: &
  EXTERNAL_BIC_SET = .false. ! external BC description/setup has been found
type(icbc), dimension(:), public, pointer :: &
  EXTERNAL_BC=>null() ! external (non emepctm) BCs detailed description/setup
type(icbc), dimension(NSPEC_ADV), private, target, save :: &
  map_bc              ! detailed description/setup from ExternalBICs_bc namelist

character(len=*),private, parameter :: &
  DEBUG_FMT="(A,' DEBUG: ',A,' ''',A,'''.')"

contains

subroutine set_extbic_id(idate)
!----------------------------------------------------------------------------!
! Read external (non emepctm) BCs description/setup.
! ICs are assumed to come from emepctm (Nest_mod.init_icbc).
!
! EXTERNAL_BIC_SET  BC description/setup has been found
!        otherwise  Assume emepctm BCs (Nest_mod.init_icbc)
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

  if(.not.USE_EXTERNAL_BIC)then
    EXTERNAL_BIC_SET=.false.
    call PrintLog("No external BICs set",MasterProc)
    first_call = .false.
    return
  end if

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
    case(2014030100:2014091723)  ! 2014-03-01 00:00 to 2014-09-17 23:00 (avail. 2012-09-04)
      EXTERNAL_BIC_VERSION='IFS_MOZ_fnyp'
    case(2014091800:)            ! from 2014-09-18 00:00
      EXTERNAL_BIC_VERSION='IFS_CMP_g4e2'
    end select
    BC_DAYS=5   ! if BC file is not found, look for 1..5-day old files
  case("IFS_MOZ_f7kn","IFS_MOZ_fkya","IFS_MOZ_fnyp","IFS_CMP_g4e2") 
    BC_DAYS=5   ! explicit MACC_ENS BC mapping version
  case("MACC_EVA","REANALYSIS")  ! GRG & AER    
    select case (idate(1))
    case(:2012)
      EXTERNAL_BIC_VERSION='EVA_EU_AN'
    case(2013:)
      EXTERNAL_BIC_VERSION='EVA_EU_FC'
    end select
    BC_DAYS=1   ! if BC file is not found, look for 1-day old file
  case("EVA_EU_AN","EVA_EU_FC")
    BC_DAYS=1   ! explicit MACC_EVA BC mapping version
  case default
    EXTERNAL_BIC_VERSION='use_any'
    BC_DAYS=0   ! do not look for old BC files
  end select

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
    end if
  end do READ_NML
  if(.not.EXTERNAL_BIC_SET)then
    call PrintLog("No external BICs found",MasterProc)
    USE_EXTERNAL_BIC=.false.
    return
  end if
  if(DEBUG.and.MasterProc) write(*,DEBUG_FMT) "set_extbic", &
    date2string("BCs for YYYY-MM-DD hh type",idate),&
    trim(EXTERNAL_BIC_NAME)//"/"//trim(EXTERNAL_BIC_VERSION)

  do n = 1,size(EXTERNAL_BC%ixadv)
    EXTERNAL_BC(n)%ixadv=find_index(EXTERNAL_BC(n)%spcname,species_adv(:)%name,any_case=.true.)
    if(EXTERNAL_BC(n)%ixadv<1)then
      EXTERNAL_BC(n)%wanted=.false.
      if(MasterProc) write(*,DEBUG_FMT) "set_extbic","unknow variable",&
        trim(EXTERNAL_BC(n)%spcname)
    end if
  end do

  if(MasterProc) &
    call PrintLog("External BICs set for "//EXTERNAL_BIC_NAME)
  EXTERNAL_BIC_SET = .true.
  first_call = .false.
end subroutine set_extbic_id
subroutine set_extbic_cd(cdate)
  type(date) :: cdate
  call set_extbic_id([cdate%year,cdate%month,cdate%day,cdate%hour])
end subroutine set_extbic_cd
endmodule ExternalBICs_mod
