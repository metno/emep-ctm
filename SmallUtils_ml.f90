! <SmallUtils_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_10(3282)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2016 met.no
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
!> SmallUtils_ml.f90 - a component of the EMEP MSC-W Chemical transport Model>
!! - provides small utility routines to process test strings and key-valaue pairs
!*****************************************************************************! 
module SmallUtils_ml

!_____________________________________________________________________________
!> @brief small utility provides routines to process text strings,
!!   find array indices, write arrays.
!!
!! @author Dave Simpson+Alvaro, 1999-2016
!! Language: F-complaint
!<____________________________________________________________________________
  implicit none

  ! -- subroutines in this TEST module:

  public :: wordsplit    !> Splits input text into words
  public :: LenArray     !> count No. set strings in array
  public :: AddArray     !> Adds new char array to old
  public :: WriteArray   !! Writes out char array, one element per line
  public :: find_index   !! Finds index of item in list 
  public :: find_indices !< Finds indices of arrays of items in list 
  public :: trims        !> removes all blanks from string
  public :: key2str      ! replace keyword occurence(s) on a string by given value
  private :: skey2str    !
  private :: ikey2str
  private :: rkey2str
  public :: to_upper     !> Converts string to upper case
  public :: Self_Test    !< For testing

  private :: find_index_c, find_index_i

  integer, public, parameter :: NOT_FOUND = -999
  character(len=*), public, parameter :: NOT_SET_STRING = "NOT_SET"

  interface find_index
    module procedure find_index_c   ! For character arrays
    module procedure find_index_i   ! For integer arrays
  end interface find_index

  interface key2str   ! replace keyword occurence(s) on a string by string/integer/real value
    module procedure skey2str   ! string values
    module procedure ikey2str   ! integer values
    module procedure rkey2str   ! real values
  end interface key2str

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
      if(iw>nword_max ) then
         errcode = 2
         print *, "ERROR in WORDSPLIT : Problem at ", text
         print *,"Too many words"
         iw=iw-1
         exit
      endif
    endif
  enddo
  nwords = iw

! Remove leading spaces
  if(keep_empty.or.present(strict_separator))then
    do iw=1,nwords
      wordarray(iw)=ADJUSTL(wordarray(iw))
    enddo
  endif
end subroutine wordsplit

!============================================================================
!> LenArray counts number of elements in input array (a)
!!  which are not equal to notset string
function LenArray(a,notset) result (N)
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
!> AddArray adds elements from new array to old array
subroutine AddArray(new,old,notset,errmsg)
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
       55 format(A,I0,A)
       write(errmsg, 55)"ERROR: Max Array size (",size(old),") exceeded!"
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
!>===========================================================================
!! A series of find_index routines, for character (c) and integer (i) arrays:
!!===========================================================================
function find_index_c(wanted, list, first_only, debug)  result(Index)
  character(len=*), intent(in) :: wanted
  character(len=*), dimension(:), intent(in) :: list
  logical, intent(in), optional :: first_only
  logical, intent(in), optional :: debug
!  Output:
  integer ::   Index

  character(len=*), parameter :: &
             debug_fmt="('debug find_index ',I0, i4,':',A,A2,A)"
  logical :: debug_print, OnlyFirst
  integer :: n_match ! Count for safety
  integer :: n

  n_match  = 0
  Index =  NOT_FOUND
  debug_print=.false.;if(present(debug))debug_print=debug
  OnlyFirst=.false.;if(present(first_only))OnlyFirst=first_only

  do n = 1, size(list)
    if ( wanted == list(n) ) then
      Index = n
      n_match = n_match + 1
      if( OnlyFirst ) return
      if(debug_print) &
      print debug_fmt,n,n_match,trim(list(n)),"==",trim(wanted)
    elseif ( debug_print ) then
      print debug_fmt,n,n_match,trim(list(n)),"/=",trim(wanted)
    endif
  enddo

  if ( n_match >  1 ) then !! Too many!
    n_match = -1 * n_match
      if(debug_print) &
      print *, "debug find_index REVERSE", n_match
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
!=======================================================================
 function trims(str)  result(trimmed)
  character(len=*), intent(in) :: str
  character(len=len(str)) :: trimmed
  character :: c
  integer :: i
  
  trimmed = ''
  do i = 1, len_trim( str )
     c = str(i:i)
     if (  c == ' ' ) cycle
     trimmed = trim(trimmed) // c
  end do

 end function trims
!============================================================================
!> Function posted by SethMMorton at: 
!! http://stackoverflow.com/questions/10759375/how-can-i-write-a-to-upper-or-to-lower-function-in-f90
!> Simpler to understand than use of iachar etc. (see same web side).

Pure Function to_upper (str) Result (string)

!   ==============================
!   Changes a string to upper case
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

End Function to_upper
!============================================================================
! key2str 
!   replace occurence(s) of keyword key on string iname by value val
! module procedures
!   skey2str replace string  values (main implementation)
!   ikey2str replace integer values (uses skey2str)
!   rkey2str replace real    values (uses skey2str/ikey2str)
pure function skey2str(iname,key,val,xfmt) result(fname)
  character(len=*), intent(in):: iname,key
  character(len=*), intent(in):: val
  character(len=*), intent(in), optional :: xfmt
  character(len=len(iname))   :: fname
  character(len=len(val))     :: aux
  integer :: ind,n
  fname=iname
  ind=index(fname,trim(key))
  if(ind==0)return
  if(present(xfmt))then     ! user provided formated
    write(aux,xfmt)trim(val)
  else                      ! keyword lenght same as key
    aux=trim(val)
  endif
  n=len_trim(key)
  do while (ind>0)
    fname=fname(1:ind-1)//trim(aux)//fname(ind+n:len_trim(fname))
    ind=index(fname,trim(key))
  enddo
endfunction skey2str
pure function ikey2str(iname,key,val,xfmt) result(fname)
  character(len=*), intent(in):: iname,key
  integer, intent(in)         :: val
  character(len=*), intent(in), optional :: xfmt
  character(len=len(iname))   :: fname
  character(len=32)           :: sval
  character(len=9)            :: ifmt!="(I??.??)"
  integer :: n
  if(index(iname,trim(key))==0)then
    fname=iname
    return
  endif
  if(present(xfmt))then     ! user supplied format
    write(sval,xfmt)val
    if(index(sval,'*')>0)&    ! problem with user format,
      write(sval,"(I0)")val   !  reformat value
  else                      ! guess format from keyword (same lenght as key)
    if(val<0)then             ! negative numbers would be printed as ****
      fname=iname             ! keep as it is
      return
    endif
    n=len_trim(key)
    write(ifmt,"('(I',I0,'.',I0,')')")n,n
    write(sval,ifmt)val
  endif
  fname=trim(skey2str(iname,key,sval))
endfunction ikey2str
pure function rkey2str(iname,key,val,xfmt) result(fname)
  character(len=*), intent(in):: iname,key
  real, intent(in)            :: val
  character(len=*), intent(in), optional :: xfmt
  character(len=len(iname))   :: fname
  character(len=32)           :: sval
  character(len=9)            :: ifmt!="(F??.??)"
  integer :: n,n1
  if(index(iname,trim(key))==0)then
    fname=iname
    return
  endif
  if(present(xfmt))then     ! user supplied format
    write(sval,xfmt)val
    if(index(sval,'*')>0)&      ! problem with user format,
      write(sval,"(ES15.3)")val !   reformat value
  else                      ! guess format from keyword (same lenght as key)
    if(val<0)then             ! negative numbers would be printed as ****
      fname=iname             ! keep as it is
      return
    endif
    n=len_trim(key)
    n1=index(key,".")
    if(n1>0)then
      write(ifmt,"('(F',I0,'.',I0,')')")n,n-n1
      write(sval,ifmt)val
    else
      write(sval,"(I0)")int(val)
    endif
  endif
  fname=trim(skey2str(iname,key,sval))
endfunction rkey2str
!============================================================================
subroutine Self_test()

  character(len=100) :: text = "Here is a line,split by spaces: note, commas don't work"
  character(len=5), dimension(5) :: headers = (/ "yy", "mm", "dd", "x1", "zz" /)
  character(len=5), dimension(3) :: wanted1 = (/ "yy", "x1", "zz" /)
  character(len=6), dimension(2) :: wanted2 = (/ " yy", "x1 " /)
  character(len=6), dimension(2) :: wanted3 = (/ "zz  ", "yy  " /)
  character(len=16),dimension(6) :: wantedx  = NOT_SET_STRING
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

  print "(/,a)", "4) Self-test - to_upper  ================================="
  print *, "Upper case of AbCd efG is ", trim(to_Upper("AbCd efG"))

  print "(/,a)", "5) Self-test - key2str ==============================="
  print *, key2str('replace string:  "EmisYYYY.txt".','YYYY','2005')
  print *, key2str('replace integer: "EmisYYYY.txt".','YYYY',2005)
  print *, key2str('replace real:    "EmisYYYY.txt".','YYYY',2005.0)
  print *, key2str('string w/spaces::"EmisYYYY.txt".','YYYY','99   ','(A)')
  print *, key2str('1.23 w/2 decimals: "F.FF".' ,'F.FF' ,1.23)
  print *, key2str('1.23 w/3 decimals: "F.FFF".','F.FFF',1.23)
! note that the sting needs to be long enough to hold the formated value
  print *, key2str('1.23 w/es10.3 fmt: "FFF".        ','FFF',1.23,'(es10.3)')
  print *, key2str('1.23 w/f10.2  fmt: "FFF".        ','FFF',1.23,'(f10.2)')
end subroutine Self_test

end module SmallUtils_ml

!TSTEMX program tester
!TSTEMX   use SmallUtils_ml, only : Self_test
!TSTEMX   call Self_test()
!TSTEMX end program tester
