! <KeyValue_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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
module KeyValue_ml
  
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
  public :: Self_Test


  !-- for Read_Headers we use a key-value pair, inspired by perl's hash arrays

  integer, public, parameter :: LENKEYVAL = 30   ! max length of key or value

  type, public :: KeyVal
    character(len=LENKEYVAL) :: key
    character(len=LENKEYVAL) :: value
  end type KeyVal

  logical, private, parameter :: MY_DEBUG = .false.


contains

  !=======================================================================
  function KeyValue(KV,txt)  result(val)
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
       
  end function KeyValue
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
end module KeyValue_ml
