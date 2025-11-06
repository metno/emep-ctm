! <LandPFT_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.6>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2025 met.no
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
!> <LandPFT_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!! *************************************************************************!
!! Reads LAI maps from sources specified by LandCoverInputs%LAIsrc=
!!   1. LPJ-EMEP  from LPJ-GUESS model - data provided by Guy Schurgers & Almut Arneth (Lund University)
!! or
!!   2. ECOSG-ENORM from ECOCLIMMAP SG with EMEP adjustments
!! with both normalised to LAI factors for EMEP usage (DS)

module LandPFT_mod

use CheckStop_mod,   only: CheckStop, StopAll
!HICKS use Config_module,   only: MasterProc, PFT_MAPPINGS,GLOBAL_LAInBVOCFile
use Config_module,   only: MasterProc, GLOBAL_LAInBVOCFile
use Config_module,   only: LandCoverInputs ! for LAIsrc    = 'LPJ-EMEP' or 'ECOSG-ENORM' or ECOSG-ELAI
use Debug_module,    only:  DEBUG   ! -> DEBUG%PFT_MAPS
use GridValues_mod,  only: debug_proc, debug_li, debug_lj, glon, glat
use NetCDF_mod,      only: ReadField_CDF
use Par_mod,         only: LIMAX, LJMAX, me
use SmallUtils_mod,  only: find_index, NOT_FOUND, WriteArray, trims
use TimeDate_mod,    only: current_date, print_date  ! HICKS

implicit none
private


!/- subroutines:

  public :: MapPFT_Init
  public :: MapPFT_LAI
  public :: MapPFT_BVOC

 real, public, allocatable :: pft_lai(:,:,:)
 real, public, allocatable :: pft_bvoc(:,:,:,:)

 ! PFTs available from smoothed LPJ fields

  integer, public, save :: N_PFTS 
  character(len=15),public, save, dimension(10) :: PFT_CODES 

  !integer, public, parameter :: N_PFTS = 6
  !character(len=5),public, parameter, dimension(N_PFTS) :: PFT_CODES = &
  !      (/ "CF   ", "DF   ", "NF   ", "BF   ", "C3PFT", "C4PFT" /)

   ! Variables available:

    !NORMED character(len=5),public, parameter :: LAI_VAR =  "LAIv_"
    character(len=5),public, parameter :: LAI_VAR =  "LAIv"
    character(len=5),public, parameter, dimension(2) :: BVOC_VAR = &
        (/ "Eiso_" , "Emt_ " /)

   !skip (/ "Normed_LAIv", "LAIv_      ", "Emt_       ", "Eiso_      " /)

contains

 !==========================================================================
 subroutine MapPFT_Init()

    logical :: dbgProc
    character(len=20) :: varname, lai_src

    dbgProc = ( DEBUG%PFT_MAPS > 0 .and. MasterProc )
    lai_src = LandCoverInputs%LAIsrc

    if ( lai_src == 'LPJ-EMEP' ) then
      n_pfts = 6
      PFT_CODES(1:n_pfts) = &
         [ "CF   ", "DF   ", "NF   ", "BF   ", "C3PFT", "C4PFT" ]

    else if ( lai_src(1:7) == 'ECOSG-E' ) then! either ECOSG-ENORM or ECOSG-ELAI
      n_pfts = 6
      !note: gfortran requires that all elements in PFT_CODES have the same number 
      !      of characters
      PFT_CODES(1:n_pfts) = &
         [ "BrDeTr  ", "NeDeTr  ", "Crops_C3", "Crops_C4", "Crops   ", "LowVeg  " ]
    else
      call StopAll('LAIsrc not set!'  // lai_src )
    end if
    if (dbgProc) write(*,*) 'MapPFT_Init PFTs:',PFT_CODES(1:n_pfts)

 end subroutine MapPFT_Init

 subroutine MapPFT_LAI()! month,day)

!.....................................................................
!**    DESCRIPTION:

!    Reads the processed LPJ-based LAIv and BVOC emission potentials.
!    The LPJ data have been merged into 4 EMEP forest classes and two
!    other veg, for either C3 or C4 vegetation.
!    Normed_LAIv is relative LAI, with max value 1.0


    !integer, intent(in) :: month,day

    real    :: lpj(LIMAX,LJMAX)  ! Emissions read from file
    logical,save :: my_first_call = .true.
    logical :: update_needed, dbgProc
    integer :: pft, month, day, iday10
    integer, save :: nrecord = 0, old_nrecord=-99
    character(len=20) :: varname, lai_src
    character(len=*), parameter :: dtxt='MapPFT_LAI:'

!     if ( my_first_call ) then
! LATER:  allocate ( pft_lai(LIMAX,LJMAX,N_PFTS) )
!         my_first_call = .false.
!call Test_LocalBVOC()
!call GetMEGAN_BVOC() This code needs to be  completed still *****
!     end if

    dbgProc = ( DEBUG%PFT_MAPS > 0 .and. MasterProc )
    lai_src = LandCoverInputs%LAIsrc

    update_needed = .false.
    if ( lai_src == 'LPJ-EMEP' ) then

      nrecord = current_date%month

    else if ( lai_src(1:7) == 'ECOSG-E' ) then
     ! have new data on 5, 15 and 25th of each month, so 3 per month
     ! use 0-9,10-19,20...
      month = current_date%month
      iday10   = current_date%day/10
      iday10   = min(2,iday10)  ! Avoids new record for 30th, 31st of month
      nrecord = (month-1)*3 + iday10 + 1
    else
      call StopAll('LAIsrc not set!'  // lai_src )
    end if

    if ( my_first_call ) then
       allocate ( pft_lai(LIMAX,LJMAX,n_pfts) )
       my_first_call = .false.
    end if

    if ( old_nrecord /= nrecord ) then
      update_needed = .true.
      old_nrecord = nrecord
    end if 
    if(dbgProc) write(*,"(a,i5,L2)") 'inpftB '//print_date(), nrecord, update_needed

    ! Get LAI data:

     if ( .not. update_needed ) then
       if( dbgProc ) write(*,*) trim(dtxt//lai_src), print_date()
       return
     end if

     do pft =1, N_PFTS
           varname = trims( "Normed_" // LAI_VAR // PFT_CODES(pft) )

           if ( dbgProc ) write(*,"(a,3i5,L2)") 'pftGET'// trim(varname), &
                    me, nrecord, DEBUG%PFT_MAPS, debug_proc
           if (nrecord > 36 .or. nrecord < 1 ) then
               print "(a,i5,L2)", dtxt//'XXinpftB '//print_date(), nrecord, update_needed, month, day
               call StopAll('NREC')
           end if
           call ReadField_CDF(GLOBAL_LAInBVOCFile,varname,&
              lpj,nrecord,interpol='zero_order',needed=.true.,debug_flag=.true.)

           pft_lai(:,:,pft ) = lpj(:,:)
           if( dbgProc ) then 
             write(*,"(a20,i3,3f8.3)") dtxt//"PFT_DEBUG "//print_date()//&
                 trim(varname), pft, maxval(lpj)
                 !glon(debug_li, debug_lj), glat(debug_li, debug_lj), &
                 !lpj(debug_li, debug_lj)
           end if

     end do ! pft

  end subroutine MapPFT_LAI

 !==========================================================================

 subroutine MapPFT_BVOC(month,nbvoc)

!.....................................................................
!**    DESCRIPTION: (NOT USED!)

!    Reads the processed LPJ-based LAIv and BVOC emission potentials.
!    The LPJ data have been merged into 4 EMEP forest classes and two
!    other veg, for either C3 or C4 vegetation.
!    Normed_LAIv is relative LAI, with max value 1.0


    integer, intent(in) :: month
    integer, intent(in) :: nbvoc ! usually 3 for Eiso,Emt,Emtl

    real    :: lpj(LIMAX,LJMAX)  ! Emissions read from file
    logical :: my_first_call = .true.
    integer ::  pft, ivar
    character(len=20) :: varname

    ! Ebvoc already includes monthly LAI changes - might be wrong?

return ! JAN31TEST
     if ( my_first_call ) then
         allocate ( pft_bvoc(LIMAX,LJMAX,N_PFTS,nbvoc) )
         my_first_call = .false.
     end if

    ! Get BVOC data. Code assumes that we want isoprene first, then
    !    apinene if provided. Multiple terpenes not considered yet,
    !    but we have just Emt anyway.


     do pft =1, N_PFTS
       do ivar =1, nbvoc ! size( BVOC_USED )
           varname = trim(BVOC_VAR(ivar)) // trim(PFT_CODES(pft))

           call ReadField_CDF(GLOBAL_LAInBVOCFile,varname,&
              lpj,month,interpol='zero_order',needed=.true.,debug_flag=.false.)

           pft_bvoc(:,:,pft, ivar ) = lpj(:,:)

       end do ! ivar
     end do ! pft


 end subroutine MapPFT_BVOC

 !==========================================================================

! subroutine Test_LocalBVOC()
!
!!.....................................................................
!!**    DESCRIPTION:
!
!!    Reads the processed LPJ-based LAIv and BVOC emission potentials.
!!    The LPJ data have been merged into 4 EMEP forest classes and two
!!    other veg, for either C3 or C4 vegetation.
!!    Normed_LAIv is relative LAI, with max value 1.0
!
!
!
!    real    :: loc(LIMAX,LJMAX)  ! Emissions read from file
!    logical :: my_first_call = .true.
!    integer ::  n, pft, ivar, iVeg, iEmis
!    character(len=1000) :: varname
!    character(len=2), dimension(4) :: VegName = (/ "CF", "DF", "NF", "BF" /)
!    character(len=4), dimension(3) :: EmisName = (/ "Eiso", "Emt ", "Emtl" /)
!
!
!       varname = "Fake"
!
!       call ReadField_CDF('LOCAL_BVOC.nc',varname,&
!           loc,1,interpol='zero_order',needed=.true.,debug_flag=.true.)
!
!       if( debug_proc ) print *, "LOCAL_BVOC 0,0", loc(1, 1)
!       if( debug_proc ) print *, "LOCAL_BVOC i,j", loc(debug_li, debug_lj)
!       !print *, "LOCAL_BVOC me,i,j", me, maxval(loc)
!
!     do iVeg = 1, size(VegName)
!     do iEmis = 1, size(EmisName)
!        varname = trim(EmisName(iEmis)) // "_" // trim(VegName(iVeg))
!        call ReadField_CDF('LOCAL_BVOC.nc',varname,&
!           loc,1,interpol='zero_order',needed=.true.,debug_flag=.true.)
!       if( debug_proc ) print "(a,a,f12.3)", "LOCAL_BVOC:E ", trim(varname), loc(debug_li, debug_lj)
!     end do
!     end do
!
!  end subroutine Test_LocalBVOC
!
! !=======================================================================
end module LandPFT_mod
