! <Units_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Units_mod
use CheckStop_mod,         only: CheckStop
use ChemGroups_mod,        only: chemgroups
use ChemDims_mod,          only: NSPEC_ADV, NSPEC_SHL
use ChemSpecs_mod,         only: species_adv
use Config_module,         only: PPBINV
use PhysicalConstants_mod, only: AVOG,ATWAIR
use Pollen_const_mod,      only: pollen_check
use OwnDataTypes_mod,      only: TXTLEN_DERIV,TXTLEN_SHORT,Asc2D
use SmallUtils_mod,        only: find_index
implicit none
private

! Subroutines & Functions
public ::                   &
  Init_Units,               & ! initialize conversion arrays
  Units_Scale,              & ! unit factor for single SPC
  Group_Units,              & ! unit factors for a GROUP
  Group_Scale                 ! function version of Group_Units

interface Group_Units
  module procedure Group_Units_Asc2D,Group_Units_detail
end interface Group_Units

real, private, parameter :: &
  atwS  = 32.0,             & ! Atomic weight of Sulphur
  atwN  = 14.0,             & ! Atomic weight of Nitrogen
  atwC  = 12.0,             & ! Atomic weight of Carbon
  ugXm3 = PPBINV/ATWAIR,    & ! will be multiplied by species(?)%molwt
  ugSm3 = ugXm3*atwS,       &                       ! species(?)%sulphurs
  ugNm3 = ugXm3*atwN,       &                       ! species(?)%nitrogens
  ugCm3 = ugXm3*atwC,       &                       ! species(?)%carbons
  mgXm2 = 1e6,              & ! will be multiplied by species(?)%molwt
  mgSm2 = mgXm2*atwS,       &                       ! species(?)%sulphurs
  mgNm2 = mgXm2*atwN,       &                       ! species(?)%nitrogens
  mgCm2 = mgXm2*atwC,       &                       ! species(?)%carbons
  grainsXm3=ugXm3*1e-6,     & ! [g/m3]/grain_wt=[#/m3]
  grainsXm2=mgXm2*1e-3,     & ! [g/m2]/grain_wt=[#/m2]
! Extinction coefficient [1/m] = %ExtC [m2/g] * mass [g/m3]
  extX  = ugXm3*1e-6,       & ! will be multiplied by species(?)%molwt*%ExtC
  s2h   = 1.0/3600.           ! sec to hour conversion factor

real, public, parameter ::  &
  to_ugSIA=ugXm3,           & ! conversion to ug
  to_mgSIA=to_ugSIA*1e3,    & ! conversion to mg
  to_number_cm3=0.001*AVOG/ATWAIR,& ! from density (roa, kg/m3) to molecules/cm3
  to_molec_cm3=to_number_cm3,&    ! kg/m3=1000 g/m3=0.001*Avog/Atw molecules/cm3
  to_molec_cm2=to_molec_cm3*1e2,&
  to_number_m3=to_number_cm3*1e6  ! 1 [#/cm3]=1e6 [#/m3]

! Conversion to ug/m3
!   xn_adv(ixadv,ix,iy,k)*roa(ix,iy,k,1)*to_ug_ADV(ixadv)
! Conversion to ugX/m3
!   xn_adv(ixadv,ix,iy,k)*roa(ix,iy,k,1)*to_ug_X(ixadv)
! Hourly Output: use "ADVugXX" for ug output (ug/m3, ugC/m3, ugN/m3, ugS/m3)
!  - for ug/m3  output use ADVugXX in combination with to_ug_ADV(ixadv);
!  - for ugX/m3 output use ADVugXX in combination with to_ug_X(ixadv).
real, public, dimension(NSPEC_ADV), save  :: &
  to_ug_ADV,  & ! conversion to ug
  to_ug_C,    & ! conversion to ug of C
  to_ug_N,    & ! conversion to ug of N
  to_ug_S       ! conversion to ug of S

logical, parameter :: T=.true., F=.false.
type, public :: umap
  character(len=TXTLEN_SHORT)  :: utxt,units ! short,NetCDF units,output class
  logical :: volunit  ! volume unit (PPB output class)?
  logical :: needroa  ! need to be multiplied by air density (roa)?,
  real, dimension(0:NSPEC_ADV) :: uconv           ! conversion factor
end type umap

type, public :: group_umap
  character(len=TXTLEN_DERIV)  :: name = 'none' ! short name
  integer,pointer,dimension(:) :: iadv =>null() ! advection index
  real,   pointer,dimension(:) :: uconv=>null() ! conversion factor
end type group_umap

type(umap), public, save :: unit_map(24)=(/&
! Air concentration
  umap("mix_ratio","mol/mol",T,F,1.0),&  ! Internal model unit
  umap("mass_ratio","kg/kg" ,T,F,1.0/ATWAIR), &  ! mass mixing ratio
  umap("ppb" ,"ppb"  ,T,F,PPBINV),&
  umap("ppbC","ppbC" ,T,F,PPBINV),&
  umap("ppbN","ppbN" ,T,F,PPBINV),&
  umap("ppbS","ppbS" ,T,F,PPBINV),&
  umap("ppbh","ppb h",T,F,s2h   ),& ! PPBINV already included in AOT calculations
  umap("ug" ,"ug/m3" ,F,T,ugXm3),&  ! ug* units need to be further multiplied
  umap("ugC","ugC/m3",F,T,ugCm3),&  !   by the air density (roa) as part of the
  umap("ugN","ugN/m3",F,T,ugNm3),&  !   unit conversion
  umap("ugS","ugS/m3",F,T,ugSm3),&
! Dry/Wet deposition
  umap("mm" ,"mm"    ,F,F,1.0  ),&
  umap("mg" ,"mg/m2" ,F,F,mgXm2),&
  umap("mgC","mgC/m2",F,F,mgCm2),&
  umap("mgN","mgN/m2",F,F,mgNm2),&
  umap("mgS","mgS/m2",F,F,mgSm2),&
! Pollen concentration/deposition
  umap("Gm3","grains/m3",F,T,grainsXm3),&
  umap("Gm2","grains/m2",F,F,grainsXm2),&
! Exposure to radioactive material
  umap("uBq" ,"uBq/m3"  ,F,T,ugXm3),& ! inst/mean   exposure
  umap("uBqh","uBq h/m3",F,T,ugXm3),& ! accumulated exposure over 1 hour
  umap("mBq" ,"mBq/m2"  ,F,F,mgXm2),& ! deposition
! Aerosol optical properties
! umap("ext" ,"ext550nm",F,T,extX),&! ext* units need to be further multiplied...
! Column output
  umap("ugm2"   ,"ug/m2"        ,F,T,ugXm3),&  ! ug* units need to be further multiplied
  umap("mcm2"   ,"molec/cm2"    ,F,T,to_molec_cm2),&
  umap("e15mcm2","1e15molec/cm2",F,T,to_molec_cm2*1e-15)/)

logical, private, save :: Initialize_Units = .true.

contains

subroutine Init_Units(update)
  logical, optional :: update
  real, dimension(NSPEC_ADV) :: uconv_spec
  integer :: i
  if(present(update))then
    Initialize_Units = Initialize_Units .or. update
  endif
  if(.not.Initialize_Units) return
  Initialize_Units = .false.

! Use "ADVugXX" for ug output (ug/m3, ugC/m3, ugN/m3, ugS/m3) in Hourly Output
!   For ug/m3  output use in combination with to_ug_ADV(ixadv).
!   For ugX/m3 output use in combination with to_ug_X(ixadv).
  to_ug_ADV = ugXm3*species_adv%molwt
  to_ug_C   = ugCm3*species_adv%carbons
  to_ug_N   = ugNm3*species_adv%nitrogens
  to_ug_S   = ugSm3*species_adv%sulphurs

 do i=1,size(unit_map)
   select case (unit_map(i)%utxt)
    case("ug","mg","uBq","uBqh","mBq","ugm2","mass_ratio")
      uconv_spec = species_adv%molwt
    case("ugC","mgC","ppbC")
      uconv_spec = species_adv%carbons
    case("ugN","mgN","ppbN")
      uconv_spec = species_adv%nitrogens
    case("ugS","mgS","ppbS")
      uconv_spec = species_adv%sulphurs
    case("Gm3","Gm2")
      uconv_spec = species_adv%molwt
      call pollen_check(uconv_adv=uconv_spec)
!   case("ext")
!     uconv_spec = species_adv%molwt*species_adv%ExtC
!     uconv_spec = species_adv%molwt*Qm_grp(NSPEC_ADV,[1..NSPEC_ADV]+NSPEC_SHL,rh,...)
    case default
      uconv_spec = 1.0
   end select
   unit_map(i)%uconv(1:)=unit_map(i)%uconv(0)*uconv_spec
 end do
end subroutine Init_Units

subroutine Group_Units_Asc2D(hr_out,gspec,gunit_conv,debug,name,volunit,needroa)
  type(Asc2D), intent(in)                     :: hr_out
  integer, pointer, dimension(:), intent(out) :: gspec      ! group array of indexes
  real,    pointer, dimension(:), intent(out) :: gunit_conv ! group array of unit conv. factors
  logical, intent(in)                         :: debug
  character(len=TXTLEN_DERIV), intent(out),optional :: name  ! For output file, species names
  logical,intent(out),optional :: volunit,needroa
  character(len=TXTLEN_DERIV)  :: dname
  character(len=*), parameter :: dtxt = 'GrpUniAsc:' ! debug text
  integer :: i

  if(Initialize_Units) call Init_Units
  call CheckStop((hr_out%spec<1).or.(hr_out%spec>size(chemgroups)),&
    dtxt//"Group_Units Error: Unknown group id, "//&
    "variable "//trim(hr_out%name)//" type "//trim(hr_out%type))

  dname=trim(chemgroups(hr_out%spec)%name)//"_"//trim(hr_out%unit)
  if(present(name))name = trim(dname)

  if(associated(gspec)) deallocate(gspec)
  allocate(gspec(size(chemgroups(hr_out%spec)%specs)))
  gspec=chemgroups(hr_out%spec)%specs-NSPEC_SHL
  if(debug) write(*,"(A,'=',30(A,':',I0,:,'+'))") &
    dtxt//trim(dname),(trim(species_adv(gspec(i))%name),gspec(i),i=1,size(gspec))

  i=find_index(hr_out%unit,unit_map(:)%utxt)
!!if(i>0)hr_out%unit=unit_map(i)%units
  if(i<1)i=find_index(hr_out%unit,unit_map(:)%units)
  call CheckStop(i<1,"Group_Units Error: Unknown unit "//trim(hr_out%unit))
  if(present(volunit)) volunit = unit_map(i)%volunit
  if(present(needroa)) needroa = unit_map(i)%needroa

  if(associated(gunit_conv)) deallocate(gunit_conv)
  allocate(gunit_conv(size(gspec)))
  gunit_conv(:)=unit_map(i)%uconv(gspec)
end subroutine Group_Units_Asc2D

subroutine Group_Units_detail(igrp,unit,gspec,gunit_conv,debug,volunit,needroa)
  integer, intent(in)                         :: igrp
  character(len=*), intent(in)                :: unit
  integer, pointer, dimension(:), intent(out) :: gspec      ! group array of indexes
  real,    pointer, dimension(:), intent(out) :: gunit_conv ! group array of unit conv. factors
  logical, intent(in)                         :: debug
  logical,intent(out),optional :: volunit,needroa
  type(Asc2D)                                 :: hr_out
  hr_out%spec=igrp
  hr_out%unit=unit//""
  hr_out%type="Group_Units_detail"
  call Group_Units_Asc2D(hr_out,gspec,gunit_conv,debug,&
                         volunit=volunit,needroa=needroa)
end subroutine Group_Units_detail

function Group_Scale(igrp,unit,debug,volunit,needroa) result(gmap)
  integer, intent(in)          :: igrp
  character(len=*), intent(in) :: unit
  logical, intent(in)          :: debug
  type(group_umap)             :: gmap
  logical,intent(out),optional :: volunit,needroa
  type(Asc2D)                  :: hr_out
  hr_out%spec=igrp
  hr_out%unit=unit//""
  hr_out%type="Group_Scale"
  call Group_Units_Asc2D(hr_out,gmap%iadv,gmap%uconv,debug,&
                         name=gmap%name,volunit=volunit,needroa=needroa)
end function Group_Scale

subroutine Units_Scale(txtin,iadv,unitscale,unitstxt,volunit,needroa,semivol,debug_msg)
  character(len=*), intent(in) :: txtin
  integer, intent(in) :: iadv  ! species_adv index, used if > 0
  real, intent(out) :: unitscale
  character(len=*),intent(out),optional :: unitstxt
  logical,         intent(out),optional :: volunit,needroa,semivol
  character(len=*),intent(in) ,optional :: debug_msg
  character(len=len(txtin)) :: txt
  integer :: i

  if(Initialize_Units) call Init_Units()
  txt=ADJUSTL(txtin)                    ! Remove leading spaces
  do i=1,len(txt)                       ! Remove invisible character
    if(ichar(txt(i:i))==0)txt(i:i)=' '  ! char(0)
  end do
  select case (txt)
  case("ugSS","ugSS/m3","ugP","ugP/m3",&
       "mgSS","mgSS/m2","mgP","mgP/m2")
    txt=txt(1:2)
  case("micro g/m3")
    txt="ug"
  case("ug_PM","ug_PM/m3","ugC_PM","ugC_PM/m3",&
     "ugN_PM","ugN_PM/m3","ugS_PM","ugS_PM/m3")
    txt=txt(1:index(txt,'_')-1)
  case("mol/mol","mole mole-1","mixratio","vmr")
    txt="mix_ratio"
  case("kg/kg","kg kg-1","kg kg**-1","massratio","mmr")
    txt="mass_ratio"
  case("ppbv","ppbV")
    txt="ppb"
  end select
  i=find_index(txt,unit_map(:)%utxt)
  if(i<1)i=find_index(txt,unit_map(:)%units)
  call CheckStop(i<1,"Units_Scale Error: Unknown unit "// trim(txtin) )

  if(present(unitstxt))unitstxt = trim(unit_map(i)%units)
  if(present(volunit )) volunit = unit_map(i)%volunit
  if(present(needroa )) needroa = unit_map(i)%needroa
  if(present(semivol))  semivol = (txtin=='ugPM').or.(index(txtin,'_PM')>0)
  select case (iadv)
  case (-1)
! groups (called iadv==-1) do not get a scaling factor at this stage.
! A second call with a valid iadv will provide the full conversion factor.
    unitscale = 1.0
  case (0)
! return the conversion factor without the specie specific part (eg %molwt)
    unitscale = unit_map(i)%uconv(iadv)
  case (1:NSPEC_ADV)
    unitscale = unit_map(i)%uconv(iadv)
    if(present(debug_msg)) &
      call CheckStop(unitscale==0.0,"Units_Scale Error: 0.0 conversion for "//&
      trim(species_adv(iadv)%name)//" in "//trim(unitstxt)//" at "//trim(debug_msg))
  case default
    call CheckStop(iadv,"Units_Scale Error: Unknown iadv.")
  end select

end subroutine Units_Scale

endmodule Units_mod
