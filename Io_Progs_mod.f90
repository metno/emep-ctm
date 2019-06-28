! <Io_Progs_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
 module Io_Progs_mod  
!_____________________________________________________________________________
! -- routines to check and open files. Also, the routine Read_Headers
!    is designed for EMEP format input files, using specific headers.
!    Read_Headers can also be used to check that all required columns
!    are available, or to read files with columns in arbitray positions.
!
! Language: F-compliant, except system calls in Self_Test
! (Can be run with F is test-input file created manually
!  and system calls commented out, as here)
!  Dave Simpson, 1999-2007
!_____________________________________________________________________________

use CheckStop_mod,           only: CheckStop
use GridValues_mod,          only: i_local, j_local
use Io_Nums_mod,             only: IO_TMP, IO_LOG
use Config_module,           only: GRID, MasterProc, IIFULLDOM, JJFULLDOM
use Debug_module,            only: DEBUG   ! -> DEBUG%IOPROG
use MPI_Groups_mod,          only : MPI_INTEGER, MPI_CHARACTER, MPI_COMM_CALC, IERROR, ME_MPI
use KeyValueTypes,           only: KeyVal, KeyValue, LENKEYVAL
use OwnDataTypes_mod,        only: TXTLEN_NAME
use Par_mod,                 only: me, limax,ljmax
use SmallUtils_mod,          only: wordsplit, WriteArray
use TimeDate_mod,            only: date,current_date
use TimeDate_ExtraUtil_mod,  only: date2string
implicit none

! -- subroutines in this module:

public :: read_line     !  Reads one line of input on host, broadcasts to other
                        !  (done as text for flexibility)
public :: check_file    !  checks that file exists and stops if required
public :: open_file     !  checks that file exists and opens if required
public :: Read_Headers  !  Reads header information from input files
public :: Read2D        !  Reads x,y,z data for simple case
public :: Read2DN       !  Reads x,y,z1...z2 data for simple case
public :: PrintLog      !  writes message to both RunLog and unit 6
public :: datewrite     ! writes date then data - helper sub
private :: datewrite_ia,&   ! int, array vesion
           datewrite_iia,&  !  array of ints and reals version
           datewrite_a      !  array of reals
public :: Self_Test

logical, public :: fexist                      ! true if file exists
integer, public, parameter :: NO_FILE = 777    ! code for non-existing file
integer, public, save :: ios                   ! i/o error status number

integer, private, parameter :: MAXLINELEN = 9000 ! Max length of ascii inputs
integer, private, parameter :: MAXHEADERS = 900  ! Max  No. headers

interface datewrite
  module procedure datewrite_ia,datewrite_iia,datewrite_a
end interface datewrite

contains
!-------------------------------------------------------------------------
subroutine PrintLog(txt,OutputProc,ioOption)
  character(len=*), intent(in) :: txt
  logical, intent(in), optional :: OutputProc  !typically MasterProc, me==0
  integer, intent(in), optional :: ioOption    !use for other files
  logical :: ok2print
  integer :: io
  ok2print = .true.
  if ( present(OutputProc) ) ok2print = OutputProc
  if ( ok2print) then
    io = IO_LOG
    if ( present(ioOption) ) io = ioOption
    write(*,fmt='(A)')  trim(txt)
    write(io,fmt='(A)')  trim(txt)
  end if
end subroutine PrintLog
!-------------------------------------------------------------------------
subroutine read_line(io_in,txt,status,label,printif)
!  Reads one line of input on host (MasterProc), broadcasts to other processors
!  This is done as text for flexibility, with the inten
!
!    Instead of e.g.
!      if ( MasterProc ) read(unit=IO,fmt=*)  i,j, data(:)
!      call MPI_BROADCAST(....)
!    We use
!      call read_line(IO,txtinput)
!      read(unit=txtinput,fmt=*) i,j, data(:)
!
!  Why? To let read_line hide the sending of data across processors
!  in the MPI framework. Above, txtinput is made available to all
!  processors.
!-------------------------------------------------------------------------
  integer, intent(in) :: io_in
  character(len=*), intent(inout) :: txt
  character(len=len(txt)+30) :: errmsg
  integer, intent(out) :: status
  character(len=*), intent(in), optional :: label
  logical, intent(in), optional :: printif   ! Can switch debug printouts
  logical :: ok2print
  character(len=40) :: label2
  label2 = " No-label"
  ok2print = .true.
  if( present(label)  ) label2   = label
  if( present(printif)) ok2print = printif
     

  if ( MasterProc ) then
    txt = ""
    read(unit=io_in,fmt="(a)",iostat=status) txt

    if ( len_trim(txt) > 0.9*MAXLINELEN ) then ! line too long for comfort
      write(unit=errmsg,fmt=*) "ERROR? Increase MAXLINELEN for IO", &
        io_in, len_trim(txt), "txt = "
      call CheckStop ( errmsg // txt )
    end if

    if ( DEBUG%IOPROG ) then ! nb already MasterProc
      if( ok2print ) write(unit=*,fmt="(a,i3,2a,i5,a,a,i4)") &
        "IOREADLINE ", io_in, trim(label2), " Len ", len_trim(txt), &
        "TXT:" //  trim(txt), " Stat ", status
    end if
  end if
   
  call MPI_BCAST( txt, len(txt), MPI_CHARACTER, 0, MPI_COMM_CALC,IERROR)
  call MPI_BCAST( status, 1, MPI_INTEGER, 0, MPI_COMM_CALC,IERROR)
  if ( DEBUG%IOPROG .and. me==1 ) then
    write(unit=errmsg,fmt=*) "proc(me) ", me, " BCAST_LINE:" // trim(txt)
    write(unit=*,fmt=*) trim(errmsg)
  end if
  CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)

end subroutine read_line
!-------------------------------------------------------------------------
subroutine check_file(fname,fexist,needed,errmsg)
! Checks for the existence of a file. If the file is
! specified as "needed", and missing, an error message is
! printed and the run is stopped.
!-------------------------------------------------------------------------
  character (len=*), intent(in)   :: fname    ! file name
  logical,           intent(out)  :: fexist   ! file exists
  logical,           intent(in)   :: needed   ! see below
  character (len=*), intent(inout):: errmsg

  errmsg = "ok"
  inquire(file=fname,exist=fexist)

  if(DEBUG%IOPROG)write(unit=6,fmt='(A)') "check_file::: ", fname
  if ( .not. fexist .and. .not. needed ) then
    write(unit=6,fmt='(A)') "not needed, skipping....." // trim(fname)
    ios = 0
  elseif ( .not. fexist .and. needed ) then
    ios = -1
    if(MasterProc) print *, "ERROR: Missing!!! in check-file:" // trim(fname)
    call CheckStop("Missing!!! in check-file:" // trim(fname))
  else
    if(MasterProc) write(unit=6,fmt='(A)') "IO check_file: Reading ",trim(fname)
  end if
end subroutine check_file
!-------------------------------------------------------------------------
subroutine open_file(io_num,mode,fname,needed,skip,iostat)
! Checks for the existence of a file and opens if present. If the
! file is specified as "needed", and missing, an error message is
! printed and the run is stopped.
!-------------------------------------------------------------------------
  integer,           intent(in) :: io_num   ! i/o number
  character (len=*), intent(in) :: mode     ! "r" for read, "w" for write
  character (len=*), intent(in) :: fname    ! file name
  logical, optional, intent(in) :: needed   ! see below
  integer, optional, intent(in) :: skip     ! No. text lines to be skipped
  integer, optional, intent(out) :: iostat  ! return ios

  integer :: i  ! local loop counter

  ios = 0   ! Start with  assumed ok status

  inquire(file=fname,exist=fexist)
  select case (mode)
  case ("r")
    if ( .not. fexist ) then
      call CheckStop( needed,"OPEN_FILE ::: missing file is :: "//trim(fname))
      ios = NO_FILE
    else
      open(unit=io_num,file=fname,status="old",action="read",iostat=ios)
      if( MasterProc .and. DEBUG%IOPROG) write(unit=6,fmt=*) "File opened: ", fname, ios
      ! *** skip header lines if requested ****
      if ( present( skip ) ) then ! Read (skip) some lines at start of file
        do i = 1, skip
          read(unit=io_num,fmt=*)
        end do
      end if ! skip
    end if
  case ("w")
    if ( .not. fexist ) then  ! Super-fussy coding!
      open(unit=io_num,file=fname,action="write",&
           status="new",iostat=ios)
    else
      open(unit=io_num,file=fname,action="write",&
           status="replace",position="rewind",iostat=ios)
    end if
    write(unit=6,fmt=*) "File created: ", trim(fname)
  case default
    print *, "OPEN FILE: Incorrect mode: ", trim(mode)
    ios = -1
  end select
  if(present(iostat))iostat=ios
end subroutine open_file
!-------------------------------------------------------------------------
subroutine Read_Headers(io_num,io_msg,NHeaders,NKeys,Headers,Keyvalues,&
    CheckValues, required_fields, alternate_fields )  !<= Optional
! Reads the header lines of an EMEP format input file, and extracts
! any key-value pairs, as well as the column headers. See Self_Test
! routine at end for example
!-------------------------------------------------------------------------
  integer, intent(in)                   :: io_num
  character(len=*), intent(inout)       :: io_msg
  integer, intent(out)                  :: NHeaders, NKeys
  character(len=*),dimension(:), intent(out) :: Headers
  type(KeyVal), dimension(:), intent(out) :: KeyValues
  type(KeyVal), dimension(:), intent(in), optional  :: &
    CheckValues !   Sets of key-values which must be present.

  character(len=*),dimension(:), intent(in), optional :: required_fields
  character(len=*),dimension(:), intent(in), optional :: alternate_fields

  character(len=LENKEYVAL),dimension(size(Headers)) :: xHeaders
  character(len=LENKEYVAL)  :: key, value
  character(len=5)  :: marker   ! e.g. !> or !#
  character(len=MAXLINELEN)  :: inputline
  integer :: i, NxHeaders, ncheck

  call CheckStop(present(required_fields),&
    "Read_Headers: option 'required_fields'  not yet implemented")
  call CheckStop(present(alternate_fields),&
    "Read_Headers: option 'alternate_fields' not yet implemented")
  ! Read in file. Lines beginning with "!#" are taken as comments and skipped
  NKeys = 0
  NHeaders = 0
  Headers = ""
  xHeaders = ""
  io_msg = "ok"
  do
    inputline=""
    call read_line(io_num,inputline,ios,"From ReadHeaders")
    if ( DEBUG%IOPROG .and. MasterProc ) &
      write(*,"(a3,3i3,i6,a)") "IN ", io_num, me, ios, &
        len_trim(inputline) ,trim(inputline)
    if ( ios /= 0 ) exit  ! End of file

    if ( inputline(1:2) == ": " ) then ! Key-values
      read(unit=inputline,fmt=*,iostat=ios) marker, key, value
      call CheckStop(ios, "KeyValue Input error" // trim(inputline) )
      NKeys = NKeys + 1
      KeyValues(NKeys)%key = key
      KeyValues(NKeys)%value = value
      if ( MasterProc .and.  DEBUG%IOPROG) &
        write(unit=*,fmt="(a,i3,a,a,a)") "KEYS LINE NKeys=", &
          NKeys, trim(key), " : ", trim(value)
      cycle

    elseif ( index(inputline,"#HEADER") > 0 ) then ! Header lines
      i =  index(inputline,"#")  ! Finds e.g. #Total as well as #HEADER")
      inputline(i:) = " "        ! and gets rid of this stuff.
      call wordsplit(inputline,MAXHEADERS,xHeaders,NxHeaders,ios)
      call CheckStop(ios, "Header wordsplit error" // trim(inputline) )

      do i = 1, NxHeaders
        if ( xHeaders(i)(1:1) /= "#" .and. len_trim(xHeaders(i)) > 0 ) then
          NHeaders = NHeaders + 1
          Headers(i) = xHeaders(i)
        end if
      end do
      do i = NHeaders+1, size(Headers)
        Headers(i) = ""   ! Remove trailing txt
      end do
      if ( DEBUG%IOPROG .and. MasterProc ) then
        write(*,*) "Read_Headers sizes: ", size(xHeaders) , NHeaders
        write(*,*) "New inputline ", trim( inputline )
      end if
      cycle

    elseif ( inputline(1:3) == ":: " ) then ! WILL DO LATER
      cycle ! Maybe keys with multiple values?

    else if ( inputline(1:5) == "#DATA" ) then ! End of headers.
                                               ! Data follows.
      if ( present(CheckValues) ) then
        !Check that the values specified in CheckValues are the same
        !as those found in input file's KeyValues:
        ncheck = size(CheckValues)
        do i = 1, ncheck
          call CheckStop( KeyValue(KeyValues,CheckValues(i)%key),&
            CheckValues(i)%value ,"Comparing Values: "//CheckValues(i)%key)
        end do
      end if

      if ( MasterProc .and. DEBUG%IOPROG ) then
        write(*,*) "DATA LINE" // trim(inputline)
        write(*,*)("HEADER CHECK ", i, Headers(i), i = 1, NHeaders)
      end if
      return

    elseif ( index(inputline,"#SKIP") > 0 ) then ! Comments
      cycle

    elseif ( inputline(1:1) == "#" ) then ! Comments
      if ( MasterProc .and. DEBUG%IOPROG ) &
        write(unit=*,fmt=*) "COMMENT LINE" // trim(inputline)
      cycle

    else
      call CheckStop( NHeaders < 1,&
             "GOT TO END - NO #HEADER or #DATA STATEMENT MAYBE?")
    end if
  end do
  io_msg = "GOT TO END - NO #DATA STATEMENT MAYBE?"
end subroutine Read_Headers
!-------------------------------------------------------------------------
subroutine Read2D(fname,data2d,idata2d)

  character(len=*), intent(in) :: fname
  real, dimension(:,:), intent(out), optional :: data2d
  integer, dimension(:,:), intent(out), optional :: idata2d

  integer :: i_fdom, j_fdom, i,j
  real :: tmp
  character(len=20), dimension(3) :: Headers
  character(len=200) :: txtinput  ! Big enough to contain one full input record
  type(KeyVal), dimension(20)      :: KeyValues ! Info on units, coords, etc.
  character(len=TXTLEN_NAME) :: errmsg

  integer :: NHeaders, NKeys, Nlines, Nused

  Nlines = 0
  Nused = 0

  if (present(idata2d) ) idata2d(:,:) = 0    !***  initialise  ***
  if (present(data2d)  ) data2d (:,:) = 0.0  !***  initialise  ***

  if ( MasterProc ) then
    call open_file(IO_TMP,"r",fname,needed=.true.)
    call CheckStop(ios,"open_file error on " // fname )
  end if


  call Read_Headers(IO_TMP,errmsg,NHeaders,NKeys,Headers,Keyvalues)

  call CheckStop( errmsg , "Read2D: Read_Headers" // fname )
  call CheckStop( Headers(1), "ix" , "HeaderIX not found" // fname)
  call CheckStop( Headers(2), "iy" , "HeaderIY not found" // fname)
  call CheckStop( KeyValue(KeyValues,"Coords"),"ModelCoords" ,"Landuse: Coords??")
  call CheckStop( KeyValue(KeyValues,"GRID"),GRID ,&
          "GRID - matched to ModelConstants")

  ! The first two columns are assumed for now to be ix,iy, hence:

  Headers(1) = Headers(3)
  if ( DEBUG%IOPROG .and. MasterProc ) then
    write(*,*) "Read2D Headers" // fname, NHeaders, Headers(1)
!   call WriteArray(Headers,NHeaders,"Read2D Headers")
  end if

  READLOOP: do
    call read_line(IO_TMP,txtinput,ios,"ReadLine for "//trim(fname))
    if ( ios /= 0 ) exit   ! likely end of file
    read(unit=txtinput,fmt=*,iostat=ios) i_fdom,j_fdom,tmp
    call CheckStop ( ios, "Read2D txt error:" // trim(txtinput) )
    Nlines = Nlines + 1

    if ( i_fdom > IIFULLDOM .or. j_fdom > JJFULLDOM ) then
      if( MasterProc ) write(*,*) "WARNING: Input Data in ",&
          trim(fname)," coords outside fulldomain: ", i_fdom, j_fdom
      cycle READLOOP
    end if

    i = i_local(i_fdom)   ! Convert to local coordinates
    j = j_local(j_fdom)
    if ( i >= 1 .and. i <= limax .and. j >= 1 .and. j <= ljmax  ) then
      Nused = Nused + 1
      if ( DEBUG%IOPROG .and. i_fdom==DEBUG%IJ(1) .and. j_fdom==DEBUG%IJ(2) ) &
        write(*,*) "READ TXTINPUT", me, i_fdom, j_fdom," => ",&
              i,j,tmp, Nlines, Nused
      if (present(idata2d)) then
        idata2d(i,j) = nint(tmp)
      else
        data2d(i,j) = tmp
      end if
    end if ! i,j
  end do READLOOP

  if ( MasterProc ) then
    close(IO_TMP)
    if(DEBUG%IOPROG)write(6,*) fname // " Read2D: me, Nlines, Nused = ", me, Nlines, Nused
  end if
end subroutine Read2D
!-------------------------------------------------------------------------
subroutine Read2DN(fname,Ndata,data2d,CheckValues,HeadersRead)

  character(len=*), intent(in) :: fname
  integer, intent(in) :: Ndata     ! Number of data columns
  real, dimension(:,:,:), intent(out) :: data2d
  type(KeyVal), dimension(:), intent(in), optional  :: &
    CheckValues !   Sets of key-values which must be present.
  logical, intent(in), optional :: HeadersRead

  integer, parameter  :: NCOORDS = 2   ! for ix, iy - "simple"

  integer :: i_fdom, j_fdom, i,j,kk
  real, dimension(Ndata+NCOORDS) :: tmp
  character(len=20), dimension(Ndata+10) :: Headers
  character(len=(Ndata+10)*20) :: txtinput  ! Big enough to contain one full input record
  type(KeyVal), dimension(20)  :: KeyValues ! Info on units, coords, etc.
  character(len=50) :: errmsg

  integer :: NHeaders, NKeys, Nlines, ncheck
  logical :: Start_Needed

  if ( DEBUG%IOPROG .and. MasterProc ) &
    write(*,*) " Starting Read2DN, me ",me, " Ndata ", Ndata

  Nlines = 0

  data2d  (:,:,:) = 0.0     !***  initialise  ***

  Start_Needed = .true.
  if ( present(HeadersRead) ) then   ! Headers have already been read
    Start_Needed  = .false.
    NHeaders = -1       ! not set in this case
  end if

  !======================================================================
  if ( Start_Needed ) then
    if ( MasterProc ) then
      call open_file(IO_TMP,"r",fname,needed=.true.)
      call CheckStop(ios,"ios error on Inputs.landuse")
    end if

    call Read_Headers(IO_TMP,errmsg,NHeaders,NKeys,Headers,Keyvalues)

    call CheckStop( errmsg , "Read2D: Read_Headers" // fname )
    call CheckStop( Headers(1), "ix" , "HeaderIX not found" // fname)
    call CheckStop( Headers(2), "iy" , "HeaderIY not found" // fname)
    call CheckStop( KeyValue(KeyValues,"Coords"),"ModelCoords" ,"Landuse: Coords??")

    if ( present(CheckValues) ) then
      !Check that the values specified in CheckValues are the same as those
      !found in input file's KeyValues:
      ncheck = size(CheckValues)
      do i = 1, ncheck
        call CheckStop( KeyValue(KeyValues,CheckValues(i)%key),&
          CheckValues(i)%value ,"Comparing Values: " // CheckValues(i)%key )
      end do
    end if

    ! The first two columns are assumed for now to be ix,iy, hence:
    Headers(1:Ndata) = Headers(3:Ndata+2)
    NHeaders = NHeaders -2

  end if ! Start_Needed
  !======================================================================
   if ( DEBUG%IOPROG .and. MasterProc ) then
    write(*,*) "Read2DN for ", fname, "Start_Needed ", Start_Needed, " NHeader", NHeaders
    write(*,*)("Read2D Headers" // fname, i, " Len ", len_trim(Headers(i)), &
               " H: ", trim(Headers(i)),i = 1, NHeaders)
    !call WriteArray(Headers,NHeaders,"Read2D Headers")
   end if

   do
    call read_line(IO_TMP,txtinput,ios,"ReadLine for "//fname, &
        printif=(Nlines<5) )
    if ( ios /= 0 ) exit   ! likely end of file
    read(unit=txtinput,fmt=*,iostat=ios) i_fdom,j_fdom,( tmp(kk), kk=1,Ndata)

    call CheckStop ( ios, "Read2D txt error:" // trim(txtinput) )
    Nlines = Nlines + 1

    if ( i_fdom > IIFULLDOM .or. j_fdom > JJFULLDOM ) then
      if( MasterProc ) write(*,*) "WARNING: Input Data in ",&
          trim(fname)," coords outside fulldomain: ", i_fdom, j_fdom
      cycle
    end if

    i = i_local(i_fdom)   ! Convert to local coordinates
    j = j_local(j_fdom)
    if ( i >= 1 .and. i <= limax .and. j >=1 .and. j <= ljmax  ) then
      if ( DEBUG%IOPROG .and. i_fdom==DEBUG%IJ(1) .and. j_fdom == DEBUG%IJ(2) )&
        write(*,*)"READ TXTINPUT", me, i_fdom, j_fdom, " => ", i,j,tmp(1)
      data2d(i,j,1:Ndata) = tmp(1:Ndata)
    end if ! i,j
  end do

  if ( MasterProc ) then
    close(IO_TMP)
    if(DEBUG%IOPROG)write(6,*) fname // " Read2DN: me, Nlines = ", me, Nlines
  end if
end subroutine Read2DN
!-------------------------------------------------------------------------
subroutine datewrite_ia (txt,ii,array,txt_pattern)
  ! to write out date, integer + supplied data array
  character(len=*), intent(in) :: txt
  integer, intent(in) :: ii  ! any old integer. Often needed
  real, dimension(:), intent(in) :: array
  logical, intent(in), optional :: txt_pattern
  logical :: use_pattern=.false.
  use_pattern=.false.;if(present(txt_pattern))use_pattern=txt_pattern
  if(use_pattern)then
    write(*,"(a,1x, i0, 20es11.2)") "dw:" // date2string(txt,current_date), &
      ii, array
  else
    write(*,"(a,3i3,i5,1x, i0, 20es14.5)") "dw:" // trim(txt), &
      current_date%month, current_date%day, current_date%hour, &
      current_date%seconds, ii, array
  end if
end subroutine datewrite_ia
subroutine datewrite_a (txt,array,txt_pattern)
  ! to write out date + supplied data array
  character(len=*), intent(in) :: txt
  real, dimension(:), intent(in) :: array
  logical, intent(in), optional :: txt_pattern
  logical :: use_pattern=.false.
  use_pattern=.false.;if(present(txt_pattern))use_pattern=txt_pattern
  if(use_pattern)then
    write(*,"(a,1x, 20es11.3)") "dw:" // date2string(txt,current_date), &
      array
  else
    write(*,"(a,3i3,i5,1x, 20es11.3)") "dw:" // trim(txt), &
      current_date%month, current_date%day, current_date%hour, &
      current_date%seconds, array
  end if
end subroutine datewrite_a
subroutine datewrite_iia (txt,ii,array,txt_pattern)
  ! to write out date, integer + supplied data array
  character(len=*), intent(in) :: txt
  integer,  dimension(:), intent(in) :: ii  ! arrays of integers, max 5
  real, dimension(:), intent(in) :: array
  logical, intent(in), optional :: txt_pattern
  logical :: use_pattern=.false.
  integer :: Ni
  integer,  dimension(5):: iout ! arrays of integers, max 5
  Ni = size(ii)
  call CheckStop(Ni>5, "Too many integers in datewrite: only coded for 5")
  call CheckStop(maxval(ii)>9999, "Too big integer in datewrite_iia: only coded for i5")
  iout = -1
  iout(1:Ni) = ii
  use_pattern=.false.;if(present(txt_pattern))use_pattern=txt_pattern
  if(use_pattern)then
    write(*,"(a,1x, 5i5, 20es11.2)") "dw:" // date2string(txt,current_date), &
      iout, array
  else
    write(*,"(a,3i3,i5,1x, 5i5, 20es11.2)") "dw:" // trim(txt), &
      current_date%month, current_date%day, current_date%hour, &
      current_date%seconds, iout, array
  end if
end subroutine datewrite_iia
!-------------------------------------------------------------------------
subroutine Self_Test()
! The input files are designed to read nicely in gnumeric and other spread-
! sheets (excel, oocalc), and can be either space of comma separated.
!
! Lines starting with : are for key-value pairs, e.g. : year 2002
! The line following #HEADERS should contain the headings of each column
! IMPORTANT: One line of column headers *must* be provided, and the
! number of headers must match the number of data items.
! (And second lines, e.g. for units, must be commented out)
!
! All lines starting "# " are ignored, but text will show up nicest in
! spread sheets if enlcosed in quotation marks
!-------------------------------------------------------------------------
use Config_module, only: NPROC
  use Par_mod, only: me
  integer :: NHeaders, NKeyValues, i, ios
  character(len=10), dimension(10) :: Headers
  character(len=10) :: msg = "ok"
  character(len=100) :: inputline
  integer :: yy, mm, dd
  real, dimension(2) :: test_data
  type(KeyVal), dimension(10)      :: KeyValues
  integer, parameter :: IO_IN=88

  if ( MasterProc ) then
    print "(/,a)", "Self-test - Io_Progs_mod ========================="
    print *, "PROCESSOR ", me, "CREATES FILE for TEST READS "
    print *, "NPROC ", NPROC
    call open_file(IO_IN,"w","Self_Test_INPUT.csv")
    write(unit=IO_IN,fmt="((a))") &
        "# ""Example of EMEP Input file""", &
        ":  Key1 Value1",  &
        ":  year  2007",  &
        ":  version rv2_9_8" , &
        " mm yy dd v1 v2  #Total #HEADERS", &
        " -  -  - m/s m/s   -    #SKIP ", &
        "#DATA:", &
        " 02,07, 28,1.2 ,2.3, 3.5", &
        " 02, 07, 29,2.4 ,1.2, 3.6", &
        " 02,07, 30,12.2,6.7, 18.9"
    close(IO_IN)
    print *, "PROCESSOR ", me, "OPENS FILE for TEST READS "
    call open_file(IO_IN,"r","Self_Test_INPUT.csv",needed=.true.)
  end if ! MasterProc
  
  print "(/,a)", "Self-test - Read_Headers ========================"
  call Read_Headers(IO_IN,msg,NHeaders,NKeyValues, Headers, KeyValues)

  if ( me == NPROC-1 ) then
    print *, "Checking data on processor me = ", me
    do i = 1, NKeyValues
      print *, "me ", me, "Keys ", i, &
          trim(KeyValues(i)%key), " => ",  trim(KeyValues(i)%value)
    end do
  
    print *, "NHead ", NHeaders
    do i = 1, NHeaders
      print *, "Headers ", i, trim(Headers(i))
    end do
  end if ! me
   
  print "(/,a,/,a,/)", "Self-test - Now read data =========================",&
      " REMINDER - WAS: mm yy dd v1 v2  #Total #HEADERS"
  do
    call read_line(IO_IN,inputline,ios,"ReadLine for SelfTest")
    if ( ios /= 0 ) exit
    if( MasterProc ) print *, "DATA: read_line -> ", trim(inputline)
    read(unit=inputline,fmt=*,iostat=ios) yy,mm,dd,test_data(1:2)
    if ( ios == 0 ) then
      if ( me == NPROC-1 ) &
        print *, "TEST DATA SPLIT INTO: ", yy, mm, dd, &
          test_data(1), test_data(2)
    else
      print *, "Read failed. Maybe wrong dimensions?"
    end if
  end do
end subroutine Self_Test
!-------------------------------------------------------------------------
end module Io_Progs_mod
