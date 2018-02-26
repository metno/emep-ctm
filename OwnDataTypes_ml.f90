! <OwnDataTypes_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
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
!> OwnDataTypes_ml.f90 - A component of the EMEP MSC-W Chemical transport Model
!! ***************************************************************************! 

module OwnDataTypes_ml
use NumberConstants, only : UNDEF_I, UNDEF_R
implicit none
private


public :: print_Deriv_type
public :: print_Asc2D
integer, public, parameter :: &
  TXTLEN_DERIV = 34, &
  TXTLEN_SHORT = 28, &
  TXTLEN_IND   =  6

! Contains some user-defined data-types, and routine associated
! with these. Collecting them here will
! avoid some dependencies, and shorten some My type modules.
!
! depmap used in DryDep_ml and My_WetDep_ml
! Deriv used in My_Derived and  Derived_ml
! VBST from SOA_ml

!/-- we define a type to map indices of species to be deposited
!   to the lesser number of species where Vg is calculated

type, public :: depmap
  integer :: ind   ! Index of species in IXADV_ or IX_ arrays
  integer :: calc  ! Index of species in  calculated dep arrays
  real    :: vg    ! if CDDEP_SET, give vg in m/s
end type depmap

  !==================
  !/ generic groups for integers
  type, public :: typ_i2
    integer :: int1
    integer :: int2
  end type typ_i2

  type, public :: typ_i3
    integer :: int1
    integer :: int2
    integer :: int3
  end type typ_i3
  
  !/ generic group for two (short) strings
  type, public :: typ_ss
    character(len=TXTLEN_SHORT) :: txt1='-' ! e.g. POD1_IAM_DF
    character(len=TXTLEN_SHORT) :: txt2='-' ! e.g. POD1_IAM_DF
  end type typ_ss

 !/ generic group for name and pointer to arrays
  type, public :: typ_sp
    character(len=TXTLEN_SHORT) :: name ! e.g. POD1_IAM_DF
    integer, dimension(:), pointer :: specs
  end type typ_sp

  !/ HI: generic group for name and two pointers to integer arrays
  type, public :: typ_maps
    character(len=TXTLEN_SHORT) :: name ! e.g. POD1_IAM_DF
    integer, dimension(:), pointer :: species ! like specs in typ_sp
    integer, dimension(:), pointer :: maps ! other species to map this
                                           !  one to
  end type typ_maps

  !/ generic group one (short) string & one integer
  type, public :: typ_si
    character(len=TXTLEN_SHORT) :: name
    integer :: ind
  end type typ_si
  !/ generic group for one (short) string & one shorter string
  type, public :: typ_s1ind
    character(len=TXTLEN_SHORT) :: name
    character(len=TXTLEN_IND)   :: ind  ! e.g. YMDHI
  end type typ_s1ind

!/ generic group for three (short) strings
type, public :: typ_s3
  character(len=TXTLEN_SHORT) :: txt1,txt2,txt3
end type typ_s3

!/ generic group for four (short) strings
type, public :: typ_s4
  character(len=TXTLEN_SHORT) :: txt1,txt2,txt3,txt4 ! e.g. POD1_IAM_DF
end type typ_s4

!/ generic group for five (short) strings & one integer
type, public :: typ_s5i
  character(len=TXTLEN_SHORT) :: txt1,txt2,txt3,txt4, &
                                 txt5 ! e.g. SO2,ugS,2d,AIR_CONC,SPEC
  integer                     :: ind  ! e.g. IOU_DAY
end type typ_s5i
!/ generic group for five (short) strings & one shorter string
type, public :: typ_s5ind
  character(len=TXTLEN_SHORT) :: txt1,txt2,txt3,txt4, &
                                 txt5 ! e.g. SO2,ugS,2d,AIR_CONC,SPEC,
  character(len=TXTLEN_IND)   :: ind  ! e.g. YMDHI
end type typ_s5ind

!==================
!+ Derived output type
type, public:: Deriv
  character(len=TXTLEN_DERIV) :: name     = '-' ! e.g. DDEP_SO2_m2Conif
  character(len=TXTLEN_SHORT) :: class    = '-' ! Type of data, e.g. ADV or Mosaic
  character(len=TXTLEN_SHORT) :: subclass = '-' ! e.g. "VG", "Rns"
  character(len=TXTLEN_SHORT) :: txt      = '-' ! text where needed, e.g. "Conif"
  character(len=TXTLEN_SHORT) :: unit     = '-' ! writen in netCDF output
  integer :: index         =UNDEF_I ! index in concentation array, or other
  integer :: f2d           =UNDEF_I ! index in f_2d arrays
  logical :: dt_scale      =.false. ! used only if we need a factor on dt_advec,
  real    :: scale         =UNDEF_R ! e.g. use 100.0 to get cm/s
  logical :: avg           =.true.  ! True => average data (divide by nav at end),
                                    ! else accumulate over run period
  character(len=TXTLEN_IND)   :: iotype   = '-' ! sets output timing
end type

! Sentinel values (moved to NumberConstants)
! real,    private, parameter :: UNDEF_R = -huge(0.0)
! integer, private, parameter :: UNDEF_I = -huge(0)

!==================
!+ Hourly ASCII/NetCDF output type
type, public:: Asc2D
  character(len=TXTLEN_DERIV):: name = "-"   ! Name (no spaces!)
  character(len=TXTLEN_SHORT):: type = "-"  ! "ADVppbv" or "ADVugm3" or "SHLmcm3"
! character(len=9) :: ofmt      ! Output format (e.g. es12.4)
  integer          :: spec = UNDEF_I   ! Species number in xn_adv or xn_shl array
                                ! or other arrays
  integer          :: nk = UNDEF_I     ! number of vertical levels
  character(len=TXTLEN_SHORT) :: unit   ! Unit used
  real             :: unitconv = UNDEF_R  !  conv. factor
  real             :: max      = UNDEF_R        ! Max allowed value for output
end type

!==================
!+ Defines SOA, NONVOL and VBS params
type, public :: VBST
  integer     :: index    ! just for clarity
  real        :: CiStar   ! ug/m3
 !real        :: Tref     ! Assumed 300
  real        :: DeltaH   ! kJ/mole
end type VBST
!==================

!==================
! uEMEP parameters
integer, public, parameter :: Npoll_uemep_max=7 !max number of uEMEP pollutant
integer, public, parameter :: Nsector_uemep_max=10 !max number of sectors for each uEMEP pollutant
type, public :: poll_type
  character(len=4):: emis='none'    ! one of EMIS_File: "sox ", "nox ", "co  ", "voc ", "nh3 ", "pm25", "pmco"
  integer, dimension(Nsector_uemep_max) ::sector=-1    ! sectors to be included for this pollutant. Zero is sum of all sectors
  integer :: EMIS_File_ix = 0 !index in EMIS_File (set by model)
  integer :: Nsectors = 0 !set by model
  integer :: sec_poll_ishift = 0 !The start of index for isec_poll loops
  integer     :: Nix=0      ! Number of components to take (set by model)
  integer, dimension(15) :: ix    ! Index of components to take (set by model)
  real, dimension(15)    :: mw=0.0 ! (set by model)
end type poll_type

type, public :: uEMEP_type
  integer     :: Npoll=0    ! Number of pollutants to treat in total
  integer     :: Nsec_poll=1    ! Number of sector and pollutants to treat in total
  integer     :: dist=0    ! max distance of neighbor to include. (will include a square with edge size=2*dist+1)
  integer     :: Nvert=20   ! number of k levels to include
  integer     :: DOMAIN(4) = -999
  type(poll_type) :: poll(Npoll_uemep_max) !pollutants to include
  logical     :: YEAR =.true.! Output frequency
  logical     :: MONTH =.false.
  logical     :: DAY =.false.
  logical     :: HOUR =.false.
  logical     :: HOUR_INST =.false.
  logical     :: COMPUTE_LOCAL_TRANSPORT=.true.
end type uEMEP_type

contains
!=========================================================================
subroutine print_Asc2D(w)
  type(Asc2D), intent(in) :: w  ! wanted
  write(*,*) "Prints Asc2D type ========================="
  write(*,"(a,a)")      "Name   :", trim(w%name)
  write(*,"(a,a)")      "type   :", trim(w%type)
  write(*,"(a,i4)")     "spec   :", w%spec
  write(*,"(a,a)")      "unit   :", trim(w%unit)
  write(*,"(a,i4)")     "nk     :", w%nk
  write(*,"(a,es10.3)") "unitconv:",w%unitconv
  write(*,"(a,es10.3)") "max    :", w%max
end subroutine print_Asc2D
!=========================================================================
subroutine print_Deriv_type(w)
  type(Deriv), intent(in) :: w  ! wanted
  write(*,*) "Prints Deriv type ========================="
  write(*,"(a,a)")      "Name   :", trim(w%name)
  write(*,"(a,a)")      "class  :", trim(w%class)
  write(*,"(a,a)")      "subclass:",trim(w%subclass)
  write(*,"(a,a)")      "txt    :", trim(w%txt)
  write(*,"(a,a)")      "units  :", trim(w%unit)
  if(w%index==UNDEF_I)then
    write(*,"(a)")      "index  : UNDEF"
  else
    write(*,"(a,i5)")   "index  :", w%index
  end if
  if(w%f2d==UNDEF_I)then
    write(*,"(a)")      "f2d    : UNDEF"
  else
    write(*,"(a,i3)")   "f2d    :", w%f2d
  end if
  write(*,"(a,es10.3)") "scale  :", w%scale
  write(*,*)            "dt_scale:", w%dt_scale
  write(*,*)            "avg    :", w%avg
end subroutine print_Deriv_type
!=========================================================================
endmodule OwnDataTypes_ml
