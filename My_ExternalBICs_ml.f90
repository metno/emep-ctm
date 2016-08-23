module My_ExternalBICs_ml
! External Boundary and Initial Conditions
! are set from the value depending on Experimnt Name (EXP_NAME)
use ModelConstants_ml,     only: MasterProc, DEBUG=>DEBUG_NEST_ICBC
use CheckStop_ml,          only: CheckStop
use Io_ml,                 only: PrintLog
use TimeDate_ExtraUtil_ml, only: date2string
implicit none

private
public :: set_extbic

logical, public, parameter :: &
  EXTERNAL_BIC_SET  = .false.,&
  TOP_BC = .false.

character(len=30),public, parameter :: &
  EXTERNAL_BIC_NAME = "DUMMY"

! i West/East bnd; j North/South bnd; k Top
integer,save, public :: iw=-1, ie=-1, js=-1, jn=-1, kt=-1 ! i West/East bnd; j North/South bnd; k Top

! YYYY, YY, MM, DD, hh will be replaced by numbers by the program.
! Search for date2string in set_extbic and uncomment lines, if necessary.
! For details, see detail2str in TimeDate_ExtraUtil_ml.f90
character(len=*),private, parameter :: &
  template_read_3D = 'EMEP_IN_IC.nc'         , &
  template_read_BC = 'EMEP_IN_BC_YYYYMMDD.nc', &
  template_write   = 'EMEP_OUT.nc'
! template_write   = 'EMEP_OUT_YYYYMMDD.nc'
character(len=len(template_read_3D)),public, save :: &
  filename_read_3D = template_read_3D
character(len=len(template_read_BC)),public, save :: &
  filename_read_BC = template_read_BC
character(len=len(template_write)),public, save :: &
  filename_write   = template_write

character(len=*),public, parameter :: &
  filename_eta     = 'EMEP_IN_BC_eta.zaxis'

type, public :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
  integer           :: ixadv=-1
  character(len=24) :: varname="none"
  logical           :: wanted=.false.,found=.false.
endtype icbc

type(icbc), dimension(:), public, pointer :: &
  EXTERNAL_BC=>null()

contains
subroutine set_extbic(idate)
  implicit none
  integer,intent(in) :: idate(4)

  character(len=*), parameter :: &
    DEBUG_FMT="('set_extbic DEBUG: ',A,' ''',A,'''.')"
  logical, save :: first_call=.true.
  integer :: ydmh=0

!--- Set filename from idate: on every call
  if(MasterProc.and.DEBUG) write(*,DEBUG_FMT) &
    "External BICs filenames for",EXTERNAL_BIC_NAME
! filename_read_3D=date2string(template_read_3D,idate,debug=MasterProc.and.DEBUG)
  filename_read_BC=date2string(template_read_BC,idate,debug=MasterProc.and.DEBUG)
! filename_write  =date2string(template_write  ,idate,debug=MasterProc.and.DEBUG)
  
  if(first_call) return
  call PrintLog("No external BICs set",MasterProc)
  first_call = .false.
endsubroutine set_extbic

endmodule My_ExternalBICs_ml
