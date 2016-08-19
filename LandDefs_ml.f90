! <LandDefs_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module LandDefs_ml
 use CheckStop_ml, only : CheckStop
 use Io_ml, only : IO_TMP, open_file, ios, Read_Headers, read_line
 use KeyValue_ml, only :  KeyVal
 use ModelConstants_ml, only : NLANDUSE
 use Par_ml, only : me
  implicit none
  private

!=============================================================================
! This module reads  inthe basics landuse data features, e.g. defaults
! for heights, LAI, growinf-season, etc.
! The list given below can be changed, extended or reduced, but then other
! input data files and codimg are needed.

!-----------------------------------------------------------------------------
! Notes: Basis was Emberson et al, EMEP Report 6/2000
!
! flux_wheat is an artificial species with constant LAI, SAI, h throughout year,
! to allow Fst calculations without knowing details of growing season.


  ! 2 ) Phenology part
  !/*** DESCRIPTION**********************************************************
  !/   reads in or sets phenology data used for the default deposition module
  !/   Users with own phenology data can simply provide their own subroutines
  !/   (replacing Init_phenology and Phenology)
  !/*************************************************************************

 public  :: Init_LandDefs         ! Sets table for LAI, SAI, hveg
 public  :: Growing_season 

 real, public, parameter :: STUBBLE  = 0.01 ! Veg. ht. out of season

 !/*****   Data to be read from Phenology_inputs.dat:

  type, public :: land_input
     character(len=15) :: name
     character(len=9) :: code
     character(len=3) :: type   ! Ecocystem type, see headers
     real    ::  hveg_max
     real    ::  Albedo
     integer ::  eNH4         ! Possible source of NHx
     integer ::  SGS50        ! Start of grow season at 50 deg. N
     real    ::  DSGS         ! Increase in SGS per degree N
     integer ::  EGS50        ! End of grow season at 50 deg. N
     real    ::  DEGS         ! Increase in EGS per degree N
     real    ::  LAImin       ! Min value of LAI
     real    ::  LAImax       ! Max value of LAI
     integer ::  SLAIlen      ! Length of LAI growth periods
     integer ::  ELAIlen      ! Length of LAI decline periods
  end type land_input
                                               !##############
  type(land_input), public, dimension(NLANDUSE) :: LandDefs
                                               !##############
  type(land_input), private :: LandInput

  type, public :: land_type
     logical :: is_forest
     logical :: is_conif
     logical :: is_decid
     logical :: is_crop 
     logical :: is_seminat 
     logical :: is_water
     logical :: is_ice
     logical :: is_veg
     logical :: is_bulk  ! Bulk-surface resistance used 
     logical :: is_iam   ! Fake species for IAM outputs 
  end type land_type
                                               !##############
  type(land_type), public,  dimension(NLANDUSE) :: LandType
                                               !##############
     

  logical, private, parameter :: MY_DEBUG = .false.    ! helps with extra printouts


contains
!=======================================================================
    subroutine Growing_season(lu,lat,SGS,EGS)
!=======================================================================

!   calculates the start and end of growing season for land-use
!   class "lu" and latitude "lat".  

    integer, intent(in) :: lu         ! Land-use index
    real,    intent(in) :: lat        ! Latitude 
    integer, intent(out) :: SGS, EGS  ! start and end of growing season

      if ( LandDefs(lu)%DSGS > 0 )  then ! calculate

        SGS = int ( 0.5 +  LandDefs(lu)%SGS50 + LandDefs(lu)%DSGS * (lat-50.0) )
        EGS = int ( 0.5 +  LandDefs(lu)%EGS50 + LandDefs(lu)%DEGS * (lat-50.0) )
      else
        SGS = LandDefs(lu)%SGS50
        EGS = LandDefs(lu)%EGS50
      end if

      EGS = min(EGS, 366 )  ! Keeps EGS to 366 to allow for leap year
                            ! (and ignore diff 365/366 otherwise)

  end subroutine Growing_season

  !=======================================================================
  subroutine Init_LandDefs(wanted_codes)
  !=======================================================================
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      character(len=*), dimension(:) :: wanted_codes   ! From Inputs-Landuse
      character(len=20), dimension(14) :: Headers
      character(len=200) :: txtinput  ! Big enough to contain one input record
      type(KeyVal), dimension(2) :: KeyValues ! Info on units, coords, etc.
      character(len=50) :: errmsg, fname
      integer :: iL, n, NHeaders, NKeys


      ! Read data


      fname = "Inputs_LandDefs.csv"
      if ( me == 0 ) then
         call open_file(IO_TMP,"r",fname,needed=.true.)
         call CheckStop(ios,"open_file error on " // fname )
      end if

      call Read_Headers(IO_TMP,errmsg,NHeaders,NKeys,Headers,Keyvalues)

      call CheckStop( errmsg , "Read LandDefs Headers" )
 

      !------ Read in file. Lines beginning with "!" are taken as
      !       comments and skipped

       n = 0     
       do
            call read_line(IO_TMP,txtinput,ios)
            if ( ios /= 0 ) exit   ! likely end of file
            if ( txtinput(1:1) == "#" ) cycle
            read(unit=txtinput,fmt=*,iostat=ios) LandInput
            call CheckStop ( ios, fname // " txt error:" // trim(txtinput) )
            n = n + 1
           !############################
            LandDefs(n) = LandInput
           !############################

        !/ Set any input negative values to physical ones (some were set as -1)

           LandDefs(n)%hveg_max = max( LandDefs(n)%hveg_max, 0.0)
           LandDefs(n)%LAImax   = max( LandDefs(n)%LAImax,   0.0)


            if ( MY_DEBUG .and. me == 0 ) write(*,*) "LANDPHEN match? ", n, &
                   LandInput%name, LandInput%code, wanted_codes(n)
            call CheckStop(  LandInput%code, wanted_codes(n), "MATCHING CODES in LandDefs")

            LandType(n)%is_water  =  LandInput%code == "W" 
            LandType(n)%is_ice    =  LandInput%code == "ICE" 
            LandType(n)%is_iam    =  LandInput%code(1:4) == "IAM_" 

            LandType(n)%is_forest =  &
                ( LandInput%type == "ECF" .or. LandInput%type == "EDF" )
            LandType(n)%is_conif = ( LandInput%type == "ECF"  )
            LandType(n)%is_decid = ( LandInput%type == "EDF"  )
            LandType(n)%is_crop  = ( LandInput%type == "ECR"  )
            LandType(n)%is_seminat  = ( LandInput%type == "SNL"  )
            LandType(n)%is_bulk   =  LandInput%type == "BLK" 
            LandType(n)%is_veg    =  LandInput%type /= "U" .and. &
                  LandInput%hveg_max > 0.01   ! Excludes water, ice, desert 
       end do
       if ( me == 0 ) close(unit=IO_TMP)

  end subroutine Init_LandDefs

end module LandDefs_ml
