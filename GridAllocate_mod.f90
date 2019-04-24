! <GridAllocate_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*************************************************************************! 
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
!*************************************************************************! 
module GridAllocate_mod
!____________________________________________________________________
!
 ! GridAllocate subroutines:  GridAllocate_ij, GridAllocate_rarray
 !                            + Self_Test
 !
 ! Some gridded arrays have several data values for the same i,j coordinate
 ! associated with dfferent land or emission codes, e.g.
 ! 
 !  e.g. i=2, j=2, code =  4, data = 1.2
 !       i=2, j=2, code = 12, data = 2.3
 !       i=2, j=2, code = 87, data = 0.3
 !
 ! Or, for each i,j, we may have row inputs for lots of possible codes,
 ! e.g. for landuse types with say 17 possibilities
 !
 ! It would be wasteful to assign arrays with code-dimension 87 or 17 to 
 ! collect this data, when only a few data points might be needed for this i,j.
 !
 ! This routine "compresses" the data for grid i,j into two arrays, 
 ! such that:
 !            ngridc(i,j) = 3    ! Number of data points in coord i,j
 !            gridc(i,j,1) = 4
 !            gridc(i,j,2) = 12
 !            gridc(i,j,3) = 87
 !     
 !-- Checks if a country "code" (or landuse type, lu) whose data has just 
 !   been read in has already been found within the given grid square.
 !   If not, the array "ngridc" is incremented by one and the
 !   country (or landuse) index added to "gridc".
 !
 !   These routines are  used for emissions and landuse 
!____________________________________________________________________
!
!** includes
!
!   Depends on: CheckStop - stops code if needed (MPI-enabled)
!   Language: F
!   History:
!   ds - 2000-Jan. 2001, some re-writing + Self_Test added, May 2007
!____________________________________________________________________
  use CheckStop_mod, only : CheckStop
  use Par_mod, only : me ! TESTS
  use GridValues_mod, only : i_fdom, j_fdom
  implicit none
  private

  !/-- subroutines

   public :: GridAllocate
   public :: Self_Test

   private :: GridAllocate_ij, GridAllocate_rarray

   interface GridAllocate
      module procedure  GridAllocate_ij     ! Call for one i,j value 
      module procedure  GridAllocate_rarray ! Call with full real arrays
   end interface GridAllocate


  !========================================
  contains
  !========================================

  subroutine GridAllocate_ij(label,i,j,code,ncmax,ic,&
                             ncmaxfound,gridc,ngridc,debug_flag)

     character(len=*), intent(in) :: label   ! Type of data
     integer, intent(in) :: i,j
     integer, intent(in) :: code      ! Full code (e.g. country code)
     integer, intent(in) :: ncmax     ! Max. no countries (lu) allowed

     integer, intent(out)   :: ic            ! Index in compressed array
     integer, intent(inout) :: ncmaxfound    ! No. countries found so far
     integer, dimension(:,:,:), intent(inout) :: gridc   ! Land-codes
     integer, dimension(:,:),   intent(inout) ::ngridc   ! No. countries
     logical, intent(in), optional :: debug_flag
     logical :: debug = .false.

     integer :: nc, icc  ! local variables
     character(len=100) :: errmsg

     if ( present(debug_flag) ) debug = debug_flag
       nc=ngridc(i,j)    ! nc = no. countries known so far

       do icc = 1,nc
          if( gridc(i,j,icc) == code ) then
              !write(unit=*,fmt="(a8,a20,i3,3i6,4i4)") trim(label),  &
              if( debug ) write(*,*) trim(label) // &
                 " GridAlloc::Already listed ", me, i,j, nc , icc, code
              ic = icc
              return
          end if
       end do

       ic = nc + 1
       if ( ic > ncmax ) then
          do icc  = 1, ic 
              write(unit=*,fmt=*) "GridAlloc_ij:"//trim(label),&
                  icc,ic, code,ncmax
          end do
       
          write(unit=errmsg,fmt=*) "me", me, " i ", i, " j ", j, &
                                   " iglob ", i_fdom(i), j_fdom(j)
          call CheckStop( "GridAlloc ncmax ERROR" // label // errmsg )

       end if


       ngridc(i,j) = ngridc(i,j) + 1    ! country is a new one
       gridc(i,j,ic) = code 

       if( debug ) write(*,"(a,3i3,2x,2i4)") trim(label) // &
                 " GridAlloc::Set ", me, i,j, code, ngridc(i,j)
       if( ic >  ncmaxfound) then
           ncmaxfound = ic
 
       end if

  end subroutine GridAllocate_ij
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine GridAllocate_rarray(label,ncmax,ncmaxfound,&
          fulldata,data,gridc,ngridc)

    !  See comments above for explanation. This routine was designed for 
    !  the full array as input, e.g. landuse(:,:,17)
 
     character(len=*), intent(in) :: label ! Type of data
     integer, intent(in)  :: ncmax         ! Max. no countries (lu) allowed
     integer, intent(out) :: ncmaxfound    ! No. countries found 

     real, dimension(:,:,:), intent(in)     :: fulldata ! Full data-set
     real, dimension(:,:,:), intent(out)    :: data     ! Reduced data-set
     integer, dimension(:,:,:), intent(out) :: gridc    ! Land-codes
     integer, dimension(:,:),   intent(out) :: ngridc   ! No. countries

     integer :: i,j, nc, ic, icc, code     ! local variables
     real    :: dat

     ncmaxfound = 0
     data(:,:,:) = 0.0
     ngridc(:,:) = 0
     gridc(:,:,:) = 0

     do i = 1, size(fulldata, 1)
       do j = 1, size(fulldata, 2)

         GRIDLOOP: do code = 1, size(fulldata, 3) ! e.g. 1 .. 17 for landuse

             dat = fulldata(i,j,code)
             if ( dat == 0.0 ) then
                 cycle GRIDLOOP
             end if
                 
            ! First, check if country is already in list:
             nc= ngridc(i,j)      ! nc = no. countries/landuse known so far
             do icc = 1, nc 
               if( gridc(i,j,icc) == code ) then
                   data(i,j,icc) = data(i,j,icc) + dat
                   cycle GRIDLOOP        ! Yep, go onto to next k
               end if
             end do
  
           ! Nope, must be new. Add to ngridc and gridc:
            ngridc(i,j) = ngridc(i,j) + 1 
            ic = nc + 1

            call CheckStop( ic > ncmax , "GridAlloc ncmax ERROR" // label )

            gridc(i,j,ic) = code
            data(i,j,ic)  = dat

  
            if( ic >  ncmaxfound) then
               ncmaxfound = ic
               write(unit=*,fmt=*) "GridAlloc ", label, &
                        "  increased ncmaxfound:", i,j,ic
               write(unit=*,fmt=*) "GridAlloc gridc: ", &
                           (gridc(i,j,icc),icc=1,ncmaxfound)
               write(unit=*,fmt=*) "GridAlloc Data:  ", &
                           (data(i,j,icc),icc=1,ncmaxfound)
            end if
         end do GRIDLOOP 
       end do ! j
     end do ! i

  end subroutine GridAllocate_rarray
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine Self_Test()
    !+ Tests this module

     integer, parameter :: NCMAX=3 ! Max. no codes allowed
     integer :: i=2,j=1
     integer :: ic, icc, iland     ! Indices
     integer :: ncmaxfound         ! No. countries found so far
     integer, dimension(2,2,NCMAX) :: land   ! Land-codes
     integer, dimension(2,2)   ::nland    ! No. countries
     real, dimension(2,2,10)   :: cover   !  some landuse data, perhaps
     real, dimension(2,2,NCMAX):: ccover  !  compressed array from cover

     print *, "============================================================"
     print *, "Test 1 GridAllocate integers with up to ncmax", ncmax
     print *, "Should print 3 lines: "
     print *, " "

     land = 0
     nland = 0
     do icc = 1, 4  ! e.g. country or vegetation code
        ic = icc
        if( icc == 4 ) ic = 2  ! test of repeated code
        print *, "Starting ic ", ic
        call GridAllocate("TEST-ij",i,j,ic,ncmax,iland,&
                                 ncmaxfound,land,nland)
     end do
     print *, "Self_Test 1 found ncmax = ", ncmax
     print *, " "
     print *, "============================================================"
     print *, "Test 2 Array input"

     land = 0
     nland = 0
     cover = 0.0
     cover(1,1,1)  = 10.0
     cover(1,1,4)  =  3.0
     cover(1,1,10) =  6.0

     call GridAllocate("TEST-array",ncmax, ncmaxfound,&
                cover, ccover, land,nland)
     print *, "Self_Test 2 found ncmax = ", ncmax
     print *, " "
     print *, "============================================================"
     print *, "Test 3 Array input  SHOULD MPI_ABORT WITH nc > ncmax"
     print *, " "

     land = 0
     nland = 0
     cover(1,1,2)  = 10.0  !! Add some more data - too much
     cover(1,1,3)  =  3.0

     call GridAllocate("TEST-array",ncmax, ncmaxfound,&
                        cover, ccover, land,nland)
     print *, "ERROR IN Self_Test 3 SHOULD NOT GET HERE"

  end subroutine Self_Test

end module GridAllocate_mod
