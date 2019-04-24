!> MODULE  <KeyValueTypes.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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
module KeyValueTypes
  
! =========================================================
!  Routines for dealing with a fortran equivalent
!  of a key-value pair - a crude attempt to redproduce
!  some of the nice features of perl's hashes or python's dictionary
!
! Language : F-complaint
! History: Created May 2007, Dave
! =========================================================

implicit none

  public :: KeyValue      !  returns value for given key
  private :: KeyValue_txt  !  returns text  value for given key
  private :: KeyValue_flt  !  returns float value for given key
  private :: KeyValue_int  !  returns int   value for given key
  public :: Self_Test


  !-- for Read_Headers we use a key-value pair, inspired by perl's hash arrays

  integer, public, parameter :: LENKEYVAL = 30   ! max length of key or value

  interface KeyValue
    module procedure KeyValue_txt, KeyValue_flt, KeyValue_int
  end interface KeyValue

  type, public :: KeyValReal
    character(len=LENKEYVAL) :: key
    real                     :: flt
  end type KeyValReal

  type, public :: KeyValInt
    character(len=LENKEYVAL) :: key
    integer                  :: int
  end type KeyValInt

  type, public :: KeyVal
    character(len=LENKEYVAL) :: key
    character(len=LENKEYVAL) :: value
  end type KeyVal

  logical, private, parameter :: MY_DEBUG = .false.


contains

  !=======================================================================
  function KeyValue_txt(KV,txt)  result(val)
    type(KeyVal), dimension(:), intent(in) :: KV
     character(len=*), intent(in) :: txt
     character(len=LENKEYVAL) :: val
     integer :: i

     val = ""
     do i = 1, size(KV)
         if( KV(i)%key == trim(txt) )  then
             val = KV(i)%value
             return
         end if
     end do
       
  end function KeyValue_txt
  !=======================================================================
  function KeyValue_flt(KV,txt)  result(flt)
    type(KeyValReal), dimension(:), intent(in) :: KV
     character(len=*), intent(in) :: txt
     real     :: flt
     integer :: i

     flt = -999.0   ! not completely safe, NaN would be better
     do i = 1, size(KV)
         if( KV(i)%key == trim(txt) )  then
             flt = KV(i)%flt
             return
         end if
     end do
       
  end function KeyValue_flt
  !=======================================================================
  function KeyValue_int(KV,txt)  result(int)
    type(KeyValInt), dimension(:), intent(in) :: KV
     character(len=*), intent(in) :: txt
     integer :: int
     integer :: i

     int = -999   ! not completely safe, NaN would be better
     do i = 1, size(KV)
         if( KV(i)%key == trim(txt) )  then
             int = KV(i)%int
             return
         end if
     end do
       
  end function KeyValue_int
  !=======================================================================
  subroutine Self_Test()
     type(KeyVal), dimension(3) :: KeyValues = (/ &
                       KeyVal("Units","ppb"), &
                       KeyVal("Coords","longlat"), &
                       KeyVal("Version","2007may") /)

     print *, "Self_Test, First key: ", KeyValues(1)%key
     print *, "Self_Test, First value: ", KeyValues(1)%value
     print *, "Self_Test, using function ", KeyValue(KeyValues,"Units")
  end subroutine Self_Test
end module KeyValueTypes
