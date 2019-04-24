! <ReadField_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!__________________________________________________________________________
!
   Module ReadField_mod
!__________________________________________________________________________
!
! -----------------------------------------------------------
! Reads ascii real and integer fields, usually for the whole 
! model area, and calls global2local to distribute these to
! the calling processor. Fields initialised
!
! Initialisation to zero included, so we do not need to input an array which
! covers the whole domain.
!
! Written October 2001, HF
! Cleaned, 3-D possibility added, JEJ and DS, April-May 2007
!------------------------------------------------------------
  use CheckStop_mod,  only : CheckStop
  use Config_module, only : NPROC, IIFULLDOM, JJFULLDOM
  use MPI_Groups_mod  , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_REAL8, MPI_INTEGER, MPI_LOGICAL, &
                              MPI_MIN, MPI_MAX, MPI_SUM, &
                              MPI_COMM_CALC, MPI_COMM_WORLD, MPISTATUS, IERROR, ME_MPI, NPROC_MPI
  use Par_mod,        only : IRUNBEG,JRUNBEG              &
                            ,MAXLIMAX,MAXLJMAX           &
                            ,MSG_READ7,MSG_READ5         &
                            ,me,GIMAX,GJMAX
  use Io_mod,         only : ios, open_file
  implicit none

  integer,private  :: i, j       ! Local variables

  interface ReadField 
     module procedure ReadField_r
     module procedure ReadField_i
     module procedure ReadField_3dr
     module procedure ReadField_3di
  end interface

contains

 
 !>=========================================================================<
  subroutine ReadField_r(IO_INFILE,fname,local_field,needed_found,fill_needed)

  integer,      intent(in) :: IO_INFILE     ! File no.
  character*20, intent(in) :: fname         ! File name
  real, intent(out) :: local_field(MAXLIMAX,MAXLJMAX)! Local field
  logical, optional, intent(inout):: needed_found    ! input:needed, 
                                                     ! output:found
  logical, optional, intent(in) :: fill_needed  ! If field has to fill
                                                ! whole domain
  logical, dimension(IIFULLDOM,JJFULLDOM) :: cell_set != .false.
  logical :: needed
  character*70 :: errmsg 
  real :: tmpin   ! To allow more than one input line per i,j

  !  Initialisation to zero added, so now we do not need to input an array which
  !  covers the whole domain. 

  real :: in_field(IIFULLDOM,JJFULLDOM)  ! Field to be read

  needed=.true.
  if(present(needed_found))needed=needed_found
  in_field(:,:)    = 0.0
  local_field(:,:) = 0.0
  errmsg = "ok"

    if (me==0)then
       cell_set = .false.
       call open_file(IO_INFILE,"r",fname,needed=needed)
       if(.not.needed .and. ios/=0)then
          write(*,*)trim(fname),' not found (but not needed)'
       else

        call CheckStop(ios,"ReadField: ios error in r" // fname )

        READFIELD : do
            read(IO_INFILE,*,iostat=ios) i,j,tmpin
            if ( ios /= 0 ) exit READFIELD
            if (  i < 1 .or. i > IIFULLDOM  .or. &
                  j < 1 .or. j > JJFULLDOM  ) then  
                  errmsg = "error in i,j index in IO_INFILE="  
                  exit READFIELD
            end if
            in_field(i,j) = in_field(i,j) + tmpin
            cell_set(i,j) = .true.
         end do READFIELD

       close(IO_INFILE)
       call CheckStop( errmsg ,"ReadField_r: errmsg in ReadField")
       if( present(fill_needed) ) then
           call CheckStop( any( cell_set .eqv. .false. ) ,&
                   "ERROR: ReadField_r: cell_not_set "//trim(fname))
       end if

       ios=0
       end if
    end if !me==0

    call MPI_BCAST( ios, 1, MPI_INTEGER, 0, MPI_COMM_CALC,IERROR)

    if(ios/=0)then
       if(present(needed_found))needed_found=.false.
    else
       if(present(needed_found))needed_found=.true.
       call global2local(in_field,local_field,MSG_READ7                 &
                       ,1,IIFULLDOM,JJFULLDOM,1,IRUNBEG,JRUNBEG)
    end if
  end subroutine ReadField_r

 !>=========================================================================<

   subroutine ReadField_i(IO_INFILE,fname,local_field,needed_found,fill_needed)

  integer,      intent(in) :: IO_INFILE     ! File no.
  character*20, intent(in) :: fname         ! File name
  integer, intent(out)     :: local_field(MAXLIMAX,MAXLJMAX)
  logical, optional, intent(inout):: needed_found ! input: needed, 
                                                  ! output:found
  logical, optional, intent(in) :: fill_needed  ! If field has to fill 
                                                ! whole domain
  logical, dimension(IIFULLDOM,JJFULLDOM) :: cell_set! = .false.
  logical :: needed
  character*70 :: errmsg 
  integer :: intmp

  integer :: in_field(IIFULLDOM,JJFULLDOM)! Field to be read

  needed=.true.
  if(present(needed_found))needed=needed_found
  in_field(:,:)    = 0.0       ! Initialise - ds, 15/1/2005
  local_field(:,:) = 0.0       ! Initialise - ds, 15/1/2005
  errmsg = "ok"

   if (me==0)then
       call open_file(IO_INFILE,"r",fname,needed=needed)
       if(.not.needed .and. ios/=0)then
          write(*,*)trim(fname),' not found (but not needed)'
       else

        call CheckStop(ios,"ReadField: ios error " // fname )
        cell_set = .false.
        READFIELD : do
           read(IO_INFILE,*,iostat=ios) i,j, intmp
           if ( ios /= 0 ) exit READFIELD
           if (  i < 1 .or. i > IIFULLDOM  .or. &
                 j < 1 .or. j > JJFULLDOM  ) then  
              errmsg = "error in i,j index in IO_INFILE=" // fname
              exit READFIELD
           end if
           in_field(i,j) =  in_field(i,j) + intmp
           cell_set(i,j) = .true.
        end do READFIELD
        close(IO_INFILE)
        call CheckStop( errmsg ,"ReadField: errmsg in ReadField")
        ios=0
        if( present(fill_needed) ) then
           !GFORTRAN call CheckStop( any( cell_set == .false. ) ,&
           call CheckStop( any( cell_set .eqv. .false. ) ,&
                   "ERROR: ReadField_i: cell_not_set "//trim(fname))
        end if
       end if

    end if !me==0

    call MPI_BCAST( ios, 1, MPI_INTEGER, 0, MPI_COMM_CALC,IERROR)

    if(ios/=0)then
       if(present(needed_found))needed_found=.false.
    else
       if(present(needed_found))needed_found=.true.
       call global2local_int(in_field,local_field,MSG_READ5               &
               ,IIFULLDOM,JJFULLDOM,1,IRUNBEG,JRUNBEG)
    end if
 
  end subroutine ReadField_i
 !>=========================================================================<
 subroutine ReadField_3dr(IO_INFILE,fname,DIM3,local_field,opened)

  integer,      intent(in) :: IO_INFILE     ! File no.
  character(len=*), intent(in) :: fname     ! File name
  integer,      intent(in) :: DIM3          ! Size of k,z dimension
  real, intent(out) :: local_field(MAXLIMAX,MAXLJMAX,DIM3)! Local field

  logical, intent(in),optional :: opened

  real :: in_field(IIFULLDOM,JJFULLDOM,DIM3)! Field to be read
  character*50 :: errmsg
  real, dimension(DIM3) :: tmpin

    in_field(:,:,:)    = 0.0       ! Initialise 
    local_field(:,:,:) = 0.0       ! Initialise 
    errmsg = "ok"

    if (me==0)then
       if ( present(opened) .and. opened ) then
          write(*,*) "file ", fname, " opened before ReadField"
       else
          call open_file(IO_INFILE,"r",fname,needed=.true.)
          call CheckStop(ios,"Readfield_r Error opening " // fname )
       end if

      READFIELD : do
             read(IO_INFILE,*,iostat=ios) i,j,tmpin(:)
             if ( ios /= 0 ) exit READFIELD
             if (  i < 1 .or. i > IIFULLDOM  .or. &
                   j < 1 .or. j > JJFULLDOM  ) then
                  errmsg = "error in i,j index in IO_INFILE=" !!! ,fname, i,j
                  exit READFIELD
             end if
             in_field(i,j,:) = in_field(i,j,:) + tmpin(:)
       end do READFIELD
       close(IO_INFILE)
       call CheckStop(errmsg, "ReadField_r: error reading" // fname )

    end if !me==0

    call global2local(in_field,local_field,MSG_READ7                 &
                       ,1,IIFULLDOM,JJFULLDOM,DIM3,IRUNBEG,JRUNBEG)

  end subroutine ReadField_3dr

 !>=========================================================================<

  subroutine ReadField_3di(IO_INFILE,fname,DIM3,local_field)

  integer,      intent(in) :: IO_INFILE     ! File no.
  character*20, intent(in) :: fname         ! File name
  integer,      intent(in) :: DIM3          ! Size of k,z dimension
  integer, intent(out)     :: local_field(MAXLIMAX,MAXLJMAX,DIM3)
  character*70 :: errmsg
  integer, dimension(DIM3) :: intmp

  integer :: in_field(IIFULLDOM,JJFULLDOM,DIM3)! Field to be read

  errmsg = "ok"

   if (me==0)then
      call open_file(IO_INFILE,"r",fname,needed=.true.)
      call CheckStop(ios," ReadField_i: error opening: " // fname )

      READFIELD : do
           read(IO_INFILE,*,iostat=ios) i,j, intmp(:)
           if ( ios /= 0 ) exit READFIELD
           if (  i < 1 .or. i > IIFULLDOM  .or. &
                 j < 1 .or. j > JJFULLDOM  ) then
              errmsg = "error in i,j index in IO_INFILE=" // fname
              exit READFIELD
           end if
          in_field(i,j,:) =  in_field(i,j,:) + intmp(:)
      end do READFIELD
      close(IO_INFILE)
      call CheckStop(errmsg," ReadField_3di: error reading" // fname)

    end if !me==0

    call global2local_int(in_field,local_field,MSG_READ5               &
               ,IIFULLDOM,JJFULLDOM,DIM3,IRUNBEG,JRUNBEG)

  end subroutine ReadField_3di

 !__________________________________________________________________________
end module ReadField_mod
