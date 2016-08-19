! <Io_Progs_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
 module Io_Progs_ml  
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
!_____________________________________________________________________________

  use CheckStop_ml, only: CheckStop
  use GridValues_ml, only : i_local, j_local
  use Io_Nums_ml,   only: IO_TMP
  use ModelConstants_ml, only : DEBUG_i, DEBUG_j, DomainName
  use KeyValue_ml,  only: KeyVal, KeyValue, LENKEYVAL
  use Par_ml, only: me, li0, li1, lj0, lj1
  use SmallUtils_ml, only : wordsplit, WriteArray
  implicit none

   INCLUDE 'mpif.h' !MPI needed


  ! -- subroutines in this module:

  public :: read_line    !  Reads one line of input on host, broadcasts to other
                         !  (done as text for flexibility)
  public :: check_file   !  checks that file exists and stops if required
  public :: open_file    !  checks that file exists and opens if required
  public :: Read_Headers !  Reads header information from input files
  public :: Read2D       !  Reads x,y,z data for simple case
  public :: Read2DN      !  Reads x,y,z1...z2 data for simple case
  public :: Self_Test

  logical, public :: fexist                      ! true if file exists
  integer, public, parameter :: NO_FILE = 777    ! code for non-existing file
  integer, public, save :: ios                   ! i/o error status number

  integer, private, parameter :: MAXLINELEN = 9000 ! Max length of ascii inputs
  integer, private, parameter :: MAXHEADERS = 900  ! Max  No. headers
  logical, private, parameter :: MY_DEBUG = .false.


contains

  !=======================================================================
  subroutine read_line(io_in,txt,status)
  !=======================================================================
    !  Reads one line of input on host (me==0), broadcasts to other processors
    !  This is done as text for flexibility, with the inten
    !
    !    Instead of e.g.
    !                    if ( me == 0 ) then
    !                        read(unit=IO,fmt=*)  i,j, data(:) on a serial code, or
    !                    end if
    !                    call MPI_BROADCAST(....)
    !
    !    We use          call read_line(IO,txtinput)
    !                    read(unit=txtinput,fmt=*) i,j, data(:)
    !
    !  Why? To let read_line hide the sending of data across processors
    !  in the MPI framework. Above, txtinput is made available to all
    !  processors.

     integer, intent(in) :: io_in
     character(len=*), intent(inout) :: txt
     character(len=len(txt)+30) :: errmsg
     integer, intent(out) :: status
     integer :: INFO

     if ( me == 0 ) then
        txt = ""
        read(unit=io_in,fmt="(a)",iostat=status) txt

        if ( len_trim(txt) > 0.9*MAXLINELEN ) then ! line too long for comfort 
           write(unit=errmsg,fmt=*) "ERROR? Increase MAXLINELEN for IO", &
                io_in, len_trim(txt), "txt = "
           call CheckStop ( errmsg // txt )
        end if

        if ( MY_DEBUG ) then
           write(unit=*,fmt=*) "READTXT" //  trim(txt)
           write(unit=*,fmt=*) "READLEN", len_trim(txt), MAXLINELEN
           write(unit=*,fmt=*) "READ_LINE ", " STATUS ", status , trim(txt)
        end if
     end if

     call MPI_BCAST( txt, len(txt), MPI_CHARACTER, 0, MPI_COMM_WORLD,INFO)
     call MPI_BCAST( status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,INFO)
     if ( MY_DEBUG ) then
        write(unit=errmsg,fmt=*) "me ", me, " BCAST_LINE:" // trim(txt)
        write(unit=*,fmt=*) trim(errmsg)
     end if
     CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)

  end subroutine read_line

  !=======================================================================
  subroutine check_file(fname,fexist,needed,errmsg)
  !=======================================================================
    ! Checks for the existence of a file. If the file is
    ! specified as "needed", and missing, an error message is
    ! printed and the run is stopped.

    character (len=*),   intent(in)  :: fname    ! file name
    logical,             intent(in)  :: needed   ! see below
    logical,             intent(out) :: fexist   ! file exists
    character (len=*),   intent(inout):: errmsg

    errmsg = "ok"
    inquire(file=fname,exist=fexist)

    write(unit=6,fmt=*) "check_file::: ", fname
    if ( .not. fexist .and. .not. needed ) then
       write(unit=6,fmt=*) "not needed, skipping....."
       ios = 0

    else if ( .not. fexist .and. needed ) then
       ios = -1
       print *, "ERROR: Missing!!! in check-file"

    else
       write(unit=6,fmt=*) "ok. File exists"
    end if
  end subroutine check_file

  !=======================================================================
  subroutine open_file(io_num,mode,fname,needed,skip)
  !=======================================================================

  ! Checks for the existence of a file and opens if present. If the
  ! file is specified as "needed", and missing, an error message is
  ! printed and the run is stopped.
  
    integer,             intent(in)  :: io_num   ! i/o number
    character (len=*),   intent(in)  :: mode     ! "r" for read, "w" for write
    character (len=*),   intent(in)  :: fname    ! file name
    logical, optional,   intent(in)  :: needed   ! see below
    integer, optional,   intent(in)  :: skip     ! No. text lines to be skipped

    integer :: i  ! local loop counter

    ios = 0   ! Start with  assumed ok status

     inquire(file=fname,exist=fexist)
     select case (mode)
     case ("r")

        if ( .not. fexist ) then

               call CheckStop( needed, "OPEN_FILE ::: missing file is :: "// fname )
               ios = NO_FILE
         else
           open(unit=io_num,file=fname,status="old",action="read",iostat=ios)
           write(unit=6,fmt=*) "File opened: ", fname, ios
           ! *** skip header lines if requested ****
           if ( present( skip ) ) then ! Read (skip) some lines at start of file
              do i = 1, skip
                  read(unit=io_num,fmt=*)
              end do
           end if ! skip
         end if

     case ("w")
        if ( .not. fexist ) then  ! Super-fussy coding!
          open(unit=io_num,file=fname,status="new",&
                  action="write",iostat=ios)
        else
          open(unit=io_num,file=fname,status="replace",&
                  position="rewind", &
                  action="write",iostat=ios)
        end if
        write(unit=6,fmt=*) "File created: ", fname
     case default
          print *, "OPEN FILE: Incorrect mode: ", mode
          ios = -1
     end select

  end subroutine open_file

  !=======================================================================
  subroutine Read_Headers(io_num,io_msg,NHeaders,NKeys,Headers,Keyvalues,&
      CheckValues, required_fields, alternate_fields )  !<= Optional
  !=======================================================================
    ! Reads the header lines of an EMEP format input file, and extracts
    ! any key-value pairs, as well as the column headers. See Self_Test
    ! routine at end for example
    !

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
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      !------ Read in file. Lines beginning with "!#" are taken as
      !       comments and skipped

       NKeys = 0
       NHeaders = 0
       Headers = ""
       xHeaders = ""
       io_msg = "ok"

       do

         inputline=""
         call read_line(io_num,inputline,ios)
         if (MY_DEBUG .and. me == 0 ) write(*,*) "IN ", io_num, me, ios, &
               len_trim(inputline) ,trim(inputline)
         if ( ios /= 0 ) then  ! End of file
              exit
         end if

         if ( inputline(1:2) == ": " ) then ! Key-values

              read(unit=inputline,fmt=*,iostat=ios) marker, key, value
              call CheckStop(ios, "KeyValue Input error" // trim(inputline) )
              NKeys = NKeys + 1
              KeyValues(NKeys)%key = key
              KeyValues(NKeys)%value = value
              if ( me == 0 ) then
                 if (MY_DEBUG) write(unit=*,fmt=*) "KEYS FULL =", trim(inputline)
                 if (MY_DEBUG) write(unit=*,fmt=*) "KEYS LINE NKeys=", &
                    NKeys, trim(key), " : ", trim(value)
              end if
              cycle

         else if ( index(inputline,"#HEADER") > 0 ) then ! Header lines


              call wordsplit(inputline,MAXHEADERS,xHeaders,NxHeaders,ios)

              call CheckStop(ios, "Header wordsplit error" // trim(inputline) )

              do i = 1, NxHeaders
                if ( xHeaders(i)(1:1) /= "#" .and. &
                     len_trim(xHeaders(i)) > 0 ) then
                   Nheaders = Nheaders + 1
                   Headers(i) = xHeaders(i)
                end if

              end do
              do i = Nheaders+1, size(Headers)
                  Headers(i) = ""   ! Remove trailing txt
              end do

              cycle

         else if ( inputline(1:3) == ":: " ) then ! WILL DO LATER
              cycle ! Maybe keys with multiple values?

         else if ( inputline(1:5) == "#DATA" ) then ! End of headers. 
                                                    ! Data follows.
   
              if ( present(CheckValues) ) then
                 !Check that the values specified in CheckValues are the same 
                 !as those found in input file's KeyValues:
                  ncheck = size(CheckValues)
                  do i = 1, ncheck
                     call CheckStop( KeyValue(KeyValues,CheckValues(i)%key),&
                        CheckValues(i)%value ,&
                              "Comparing Values: " // CheckValues(i)%key )
                   end do
              end if

              if ( me == 0 .and. MY_DEBUG ) then
               write(unit=*,fmt=*) "DATA LINE" // trim(inputline)
              end if

              return

         else if ( index(inputline,"#SKIP") > 0 ) then ! Comments

              cycle

         else if ( inputline(1:1) == "#" ) then ! Comments

              if ( MY_DEBUG )  write(unit=*,fmt=*) &
                      "COMMENTS LINE" // trim(inputline)
              cycle

         end if

       end do

       io_msg = "GOT TO END - NO #DATA STATEMENT MAYBE?"

  end subroutine Read_Headers
  !=======================================================================
  subroutine Read2D(fname,data2d,idata2d)

   character(len=*), intent(in) :: fname
   real, dimension(:,:), intent(out), optional :: data2d
   integer, dimension(:,:), intent(out), optional :: idata2d

   integer :: i_fdom, j_fdom, i,j
   real :: tmp
   character(len=20), dimension(3) :: Headers
   character(len=200) :: txtinput  ! Big enough to contain one full input record
   type(KeyVal), dimension(20)      :: KeyValues ! Info on units, coords, etc.
   character(len=50) :: errmsg

   integer :: NHeaders, NKeys, Nlines
   logical :: debug_flag


   Nlines = 0

   if (present(idata2d) ) idata2d  (:,:) = 0    !/**  initialise  **/
   if (present(data2d)  ) data2d  (:,:) = 0.0     !/**  initialise  **/

   if ( me == 0 ) then
      call open_file(IO_TMP,"r",fname,needed=.true.)
      call CheckStop(ios,"open_file error on " // fname )
   end if


   call Read_Headers(IO_TMP,errmsg,NHeaders,NKeys,Headers,Keyvalues)

   call CheckStop( errmsg , "Read2D: Read_Headers" // fname )
   call CheckStop( Headers(1), "ix" , "HeaderIX not found" // fname)
   call CheckStop( Headers(2), "iy" , "HeaderIY not found" // fname)
   call CheckStop( KeyValue(KeyValues,"Coords"),"ModelCoords" ,"Landuse: Coords??")
   call CheckStop( KeyValue(KeyValues,"Domain"),DomainName ,&
           "Domain Name - matched to ModelConstants")

   ! The first two columns are assumed for now to be ix,iy, hence:

   Headers(1) = Headers(3)
   if ( MY_DEBUG .and. me == 0 ) then
        write(*,*) "Read2D Headers" // fname, Nheaders, Headers(1)
!       call WriteArray(Headers,NHeaders,"Read2D Headers")
   end if

   do
         call read_line(IO_TMP,txtinput,ios)
         if ( ios /= 0 ) exit   ! likely end of file
         read(unit=txtinput,fmt=*,iostat=ios) i_fdom,j_fdom,tmp
         call CheckStop ( ios, "Read2D txt error:" // trim(txtinput) )
         Nlines = Nlines + 1

         i   = i_local(i_fdom)   ! Convert to local coordinates
         j   = j_local(j_fdom)

         if ( i >= li0 .and. i <= li1 .and. j >= lj0 .and. j <= lj1  ) then  

             if ( MY_DEBUG ) debug_flag = ( i_fdom == DEBUG_i  &
                                   .and. j_fdom == DEBUG_j  )

              if ( debug_flag ) then
                write(*,*) "READ TXTINPUT", me, i_fdom, j_fdom, " => ", i,j,tmp
              endif
              if (present(idata2d)) then 
                 idata2d(i,j) = nint(tmp)
              else
                 data2d(i,j) = tmp
              end if

         end if ! i,j
      end do

      if ( me == 0 ) close(IO_TMP)
      if ( me == 0 ) write(6,*) fname // "Read2D: me, Nlines = ", me, Nlines

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
   type(KeyVal), dimension(20)      :: KeyValues ! Info on units, coords, etc.
   character(len=50) :: errmsg

   integer :: NHeaders, NKeys, Nlines, ncheck
   logical :: debug_flag, Start_Needed

   if ( MY_DEBUG .and. me == 0  ) write(*,*) " Starting Read2DN, me ",me

   Nlines = 0

   data2d  (:,:,:) = 0.0     !/**  initialise  **/

   Start_Needed = .true.
   if ( present(HeadersRead) ) then   ! Headers have already been read
       Start_Needed  = .false.
   end if

  !======================================================================
   if ( Start_Needed ) then 
  !======================================================================
   if ( me == 0 ) then
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
   if ( MY_DEBUG .and. me == 0 ) then
        write(*,*) "Read2DN for ", fname, "Start_Needed ", Start_Needed
        do i = 1, NHeaders
          write(*,*) "Read2D Headers" // fname, i, Nheaders, Headers(i)
        end do
       !call WriteArray(Headers,NHeaders,"Read2D Headers")
   end if

   do
         call read_line(IO_TMP,txtinput,ios)
         if ( ios /= 0 ) exit   ! likely end of file
         read(unit=txtinput,fmt=*,iostat=ios) i_fdom,j_fdom,&
                 ( tmp(kk), kk=1,Ndata)

         call CheckStop ( ios, "Read2D txt error:" // trim(txtinput) )
         Nlines = Nlines + 1

         i   = i_local(i_fdom)   ! Convert to local coordinates
         j   = j_local(j_fdom)

         !SAFER? if ( i >= li0 .and. i <= li1 .and. j >= lj0 .and. j <= lj1  ) then  
         if ( i >= 1 .and. i <= li1 .and. j >=1 .and. j <= lj1  ) then  

             if ( MY_DEBUG .and. &
                  i_fdom == DEBUG_i  .and. j_fdom == DEBUG_j  ) &
                write(*,*) "READ TXTINPUT", me, i_fdom, j_fdom, " => ", i,j,tmp(1)

             data2d(i,j,1:Ndata) = tmp(1:Ndata)

         end if ! i,j
      end do

      if ( me == 0 ) close(IO_TMP)
      if ( me == 0 ) write(6,*) fname // "Read2DN: me, Nlines = ", me, Nlines

  end subroutine Read2DN
 
  !-------------------------------------------------------------------------
  subroutine Self_Test()
   use ModelConstants_ml, only: NPROC
   use Par_ml, only: me
    integer :: NHeaders, NKeyValues, i, ios
    character(len=10), dimension(10) :: Headers
    character(len=10) :: msg = "ok"
    character(len=100) :: inputline
    integer :: yy, mm, dd
    real, dimension(2) :: test_data
    type(KeyVal), dimension(10)      :: KeyValues
    integer, parameter :: IO_IN=88
    

!----------------------------------------------------------------------------
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
!
    if ( me == 0 ) then
    
      print "(/,a)", "Self-test - Io_Progs_ml ==========================="

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
     end if ! me = 0
  
     print "(/,a)", "Self-test - Read_Headers =========================="

     call Read_Headers(IO_IN,msg,Nheaders,NKeyValues, Headers, KeyValues)

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
          call read_line(IO_IN,inputline,ios)
          if ( ios /= 0 ) then
            exit
          end if
          if(me==0) then
               print *, "DATA: read_line -> ", trim(inputline)
          end if
          read(unit=inputline,fmt=*,iostat=ios) yy,mm,dd,test_data(1:2)
          if ( ios == 0 ) then
             if ( me == NPROC-1 )  then 
                 print *, "TEST DATA SPLIT INTO: ", yy, mm, dd, &
                  test_data(1), test_data(2)
             end if
          else
                 print *, "Read failed. Maybe wrong dimensions?"
          end if
        end do
  

   end subroutine Self_Test

 end module Io_Progs_ml  
