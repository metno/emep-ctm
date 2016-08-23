! <SmallUtils_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2011 met.no
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
module SmallUtils_ml

!_____________________________________________________________________________
! -- small utility provides routines to process text strings,
!    find array indices, write arrays.
!
! Dave Simpson, 1999-2011
! Language: F-complaint, except system calls in Self_Test
! (Can be run with F is test-input file created manually
!  and system calls commented out, as here)
!_____________________________________________________________________________
  implicit none

  ! -- subroutines in this module:

  public :: wordsplit    !  Splits input text into words
  public :: LenArray     ! count No. set strings in array
  public :: AddArray     ! Adds new char array to old
  public :: WriteArray   ! Writes out char array, one element per line
  public :: find_index   ! Finds index of item in list 
  public :: find_indices ! Finds indices of arrays of items in list 
  public :: Self_Test    ! For testing

  private :: find_index_c, find_index_i

  integer, public, parameter :: NOT_FOUND = -999
  character(len=*), public, parameter :: NOT_SET_STRING = "NOT_SET"

  interface find_index
    module procedure find_index_c   ! For character arrays
    module procedure find_index_i   ! For integer arrays
  end interface find_index

contains

!===========================================================================

subroutine wordsplit(text,nword_max,wordarray,nwords,errcode,separator,&
                     strict_separator,empty_words)
!**************************************************************
!   Subroutine takes in a character string and splits it into
!   a word-array, of length nwords
!   Both spaces and commas are treated as seperators
!**************************************************************

 !-- arguments
  character(len=*), intent(in) ::  text       ! to be split
  integer,          intent(in) ::  nword_max  ! Max. no. words expected

  character(len=*), dimension(:), intent(out) :: wordarray
  integer,          intent(out) :: nwords      ! No. words found
  integer,          intent(out) :: errcode      ! error status
  character(len=1), optional, intent(in) ::  &
                            separator,       &  ! additional separators
                            strict_separator    ! only this separator
  logical, optional, intent(in) :: empty_words  ! keep empty strings

  !-- local
  logical   :: wasinword,&  ! remove leading spaces on request
               keep_empty   ! keep empty strings on request
  integer   :: i, is, iw
  character(len=1) ::  c,s(0:3)

  errcode = 0
  wasinword = .false.   !To be safe, with spaces at start of line
  is = 0 ! string index
  iw = 1 ! Word index
  s(:)=(/' ',' ',',',':'/)
  if(present(separator))s(0)=separator
  if(present(strict_separator))s(:)=strict_separator
  wordarray(1) = ""
  keep_empty=.false.
  if(present(empty_words))then
    keep_empty=empty_words
    wasinword=keep_empty
  endif

  do i = 1, len_trim(text)
    c = text(i:i)
    if( all(c/=s) ) then
      is = is + 1
      wordarray(iw)(is:is) = c
      wasinword = .true.
    elseif ( wasinword ) then
      iw = iw + 1
      wordarray(iw) = ""
      wasinword = keep_empty
      is = 0
    endif
  enddo
  nwords = iw
  if (  nwords > nword_max ) then
    errcode = 2
    print *, "ERROR in WORDSPLIT : Problem at ", text
    print *,"Too many words"
  endif

! Remove leading spaces
  if(keep_empty.or.present(strict_separator))then
    do iw=1,nwords
      wordarray(iw)=ADJUSTL(wordarray(iw))
    enddo
  endif
end subroutine wordsplit

!============================================================================
function LenArray(a,notset) result (N)
  !+ Counts number of elements in a which are not equal to notset string
  character(len=*), dimension(:), intent(in) :: a
  character(len=*), intent(in) :: notset
  integer :: N, i

  N=0
  do i = 1, size(a)
    if ( index(a(i),notset) > 0  ) exit
    N=N+1
  enddo
end function LenArray
!============================================================================
subroutine AddArray(new,old,notset,errmsg)
  !+ Adds elements from new array to old array
  character(len=*), dimension(:), intent(in) :: new
  character(len=*), dimension(:), intent(inout) :: old
  character(len=*), intent(in) :: notset
  character(len=*), intent(inout) :: errmsg
  integer :: N, i
  errmsg = "ok"

  N = LenArray(old,notset) ! Find last set element
  do i = 1,  size(new)
    N = N + 1
    if ( N > size(old) ) then
      errmsg = "ERROR: Array Exceeded! "
      return
    endif
    old(N) = new(i)
  enddo
end subroutine AddArray
!============================================================================
subroutine WriteArray(list,NList,txt,io_num)
  character(len=*), dimension(:), intent(in) :: list
  integer, intent(in) :: Nlist
  character(len=*), intent(in) :: txt   ! Some descriptive text
  integer, intent(in), optional :: io_num
  integer :: io, i

  io = 6
  if ( present(io_num) ) io = io_num

  if ( NList > size(list) ) then
    write(unit=*,fmt=*) "WRITEARRAY PROBLEM Nlist, size(List) ", &
              Nlist, size(list), trim(txt)
    return
  endif
  do i = 1, Nlist
    write(unit=io,fmt=*) txt, i, list(i)
  enddo
end subroutine WriteArray
!============================================================================
! A series of find_index routines, for character (c) and integer (i) arrays:
!============================================================================
function find_index_c(wanted, list, debug)  result(Index)
  character(len=*), intent(in) :: wanted
  character(len=*), dimension(:), intent(in) :: list
  logical, intent(in), optional :: debug
!  Output:
  integer ::   Index

  character(len=*), parameter :: &
             debug_fmt="('debug find_index ',I0,':',A,A2,A)"
  logical :: debug_print
  integer :: n_match ! Count for safety
  integer :: n

  n_match  = 0
  Index =  NOT_FOUND
  debug_print=.false.;if(present(debug))debug_print=debug

  do n = 1, size(list)
    if ( wanted == list(n) ) then
      Index = n
      n_match = n_match + 1
      if(debug_print) &
      print debug_fmt,n,trim(list(n)),"==",trim(wanted)
    elseif ( debug_print ) then
      print debug_fmt,n,trim(list(n)),"/=",trim(wanted)
    endif
  enddo

  if ( n_match >  1 ) then !! Too many!
    n_match = -1 * n_match
  endif
end function find_index_c

!============================================================================
function find_index_i(wanted, list, debug)  result(Index)
  integer, intent(in) :: wanted
  integer, dimension(:), intent(in) :: list
  logical, intent(in), optional :: debug
!  Output:
  integer ::   Index       !

  character(len=*), parameter :: &
             debug_fmt="('debug find_index ',I0,':',I0,A2,I0)"
  logical :: debug_print
  integer :: n_match ! Count for safety
  integer :: n

  n_match  = 0
  Index =  NOT_FOUND
  debug_print=.false.;if(present(debug))debug_print=debug

  do n = 1, size(list)
    if ( wanted == list(n)  ) then
      Index = n
      n_match = n_match + 1
      if(debug_print) &
      print debug_fmt,n,list(n),"==",wanted
    elseif ( debug_print ) then
      print debug_fmt,n,list(n),"/=",wanted
    endif
  enddo

  if ( n_match >  1 ) then !! Too many!
    n_match = -1 * n_match
  endif
end function find_index_i

!=======================================================================
 function find_indices(wanted, list, debug)  result(Indices)
  character(len=*), dimension(:), intent(in) :: wanted
  character(len=*), dimension(:), intent(in) :: list
  logical, intent(in), optional :: debug
!  Output:
  integer, dimension(size(wanted)) ::   Indices

  character(len=*), parameter :: &
             debug_fmt="('debug find_indices ',I0,':',A,A2,I0,':',A)"
  logical :: debug_print
  integer :: w, n

  Indices(:) = NOT_FOUND
  debug_print=.false.;if(present(debug))debug_print=debug

  do w = 1, size(wanted)
    do n = 1, size(list)
      if ( wanted(w) == list(n) ) then
        Indices(w) = n
        if(debug_print) &
        print debug_fmt,n,trim(list(n)),"==",w,trim(wanted(w))
      elseif ( debug_print ) then
        print debug_fmt,n,trim(list(n)),"/=",w,trim(wanted(w))
      endif
    enddo
  enddo
end function find_indices

!============================================================================
subroutine Self_test()

  character(len=100) :: text = "Here is a line,split by spaces: note, commas don't work"
  character(len=5), dimension(5) :: headers = (/ "yy", "mm", "dd", "x1", "zz" /)
  character(len=5), dimension(3) :: wanted1 = (/ "yy", "x1", "zz" /)
  character(len=6), dimension(2) :: wanted2 = (/ " yy", "x1 " /)
  character(len=6), dimension(2) :: wanted3 = (/ "zz  ", "yy  " /)
  character(len=16), dimension(6) :: wantedx  = NOT_SET_STRING
  character(len=100) :: errmsg
  integer, parameter :: NWORD_MAX = 99
  character(len=20), dimension(NWORD_MAX) :: words
  integer :: nwords, errcode
  
  print "(/,a)", "1) Self-test - wordsplit ================================="
  call wordsplit(text,NWORD_MAX,words,nwords,errcode)

  print *, "Found ", nwords, "words"
  print *, "Words: ", words(1:nwords)

  print "(a)", "Note - need exact text:"
  print *, "Index of spaces is ", find_index("spaces",words)
  print *, "Index of spaces: is ", find_index("spaces:",words)

  print "(/,a)", "2) Self-test - find_indices ================================="

  print *, wanted1, " Indices => ", find_indices(wanted1,headers)
  print "(a)", "Note - trailing blanks ok, leading blanks cause error:"
  print *, wanted2, " Indices => ", find_indices(wanted2,headers)
  print *, wanted3, " Indices => ", find_indices(wanted3,headers)

  print "(/,a)", "2) Self-test - WriteArray   ================================="
    
  call WriteArray(wanted1,size(wanted1),"Testing wanted1 array")
  print "(a)", "  (Should write headers array (first 4 elements) to fort.77) "
  call WriteArray(headers,4,"Testing headers array",77)

  print "(/,a)", "4) Self-test - AddArray   ================================="
    
  !call AddArray(wanted1,wanted2,errmsg)
  wantedx(1) =  "first  "
  wantedx(2) =  "second  "
  call AddArray(wanted1,wantedx,NOT_SET_STRING,errmsg)
  call WriteArray(wantedx,size(wantedx),"Testing AddArray")
end subroutine Self_test

end module SmallUtils_ml
