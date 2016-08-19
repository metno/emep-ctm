! <Landuse_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Landuse_ml

use CheckStop_ml,      only: CheckStop
use KeyValue_ml,       only: KeyVal,KeyValue, LENKEYVAL
use DO3SE_ml,   only : fPhenology
use GridAllocate_ml,only: GridAllocate
use GridValues_ml,  only: gb_glob, gb, i_fdom, j_fdom, & ! latitude, coordinates
                          i_local, j_local, &
                         debug_proc, debug_li, debug_lj
use Io_ml,          only: open_file, ios, Read_Headers, Read2DN, IO_TMP
use KeyValue_ml,    only: KeyVal,KeyValue
use LandDefs_ml,    only: Init_LandDefs, LandType, LandDefs, STUBBLE, Growing_Season
use ModelConstants_ml,  only : DEBUG_i, DEBUG_j, NLANDUSE, &
                          NPROC, IIFULLDOM, JJFULLDOM, &
                          DomainName
use Par_ml,         only: li0, lj0, li1, lj1, MAXLIMAX, MAXLJMAX, &
                          limax, ljmax, me
use SmallUtils_ml,  only: find_index, NOT_FOUND, WriteArray
use TimeDate_ml,    only: daynumber, nydays, current_date
implicit none
private


!/- subroutines:

  public :: InitLanduse
  public :: ReadLanduse
  public :: SetLanduse
  private :: Polygon         ! Used for LAI
 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO

 integer, public, parameter :: NLUMAX = 17 ! max no. landuse per grid

! The headers read from Inputs.Landuse define the "master-list" of
! codes for landuse. Each code must be present in the subsequent
! data files for phenology and DO3SE.

 character(len=6), dimension(NLANDUSE), public, save :: Land_codes ! As used

 !=============================================
 type, public :: LandCov
   integer                   :: ncodes     ! Number of codes in grid
   integer,dimension(NLUMAX) :: &
          codes     &! landcover codes
         ,SGS       &! Start of growing season (days)
         ,EGS       &! End of growing season (days)
         ,Astart    &! Start photosyntgetic activity, for DO3SE
         ,Aend       ! 
   real,   dimension(NLUMAX) :: &
          fraction  &!  (coverage) 
         ,LAI       &! Leaf-area-index (m2/m2)
         ,SAI       &! Surface-area-index (m2/m2)   (leaves+bark, etc.)
         ,hveg      &! Max. height of veg.
         ,fphen     &! Potential (age) factor for Jarvis-calc
         ,SumVPD    &! For critical VPD calcs, reset each day
         ,old_gsun   ! also for flux
 end type LandCov
 !=============================================
 type(LandCov), public, save, dimension(MAXLIMAX,MAXLJMAX) :: LandCover
 !=============================================


  integer, public,save,dimension(MAXLIMAX,MAXLJMAX) :: &
          WheatGrowingSeason   ! Growing season (days), IAM_WHEAT =1 for true
 
 ! For some flux work, experimental  XXXXXXXx

 real,public,save,dimension(MAXLIMAX,MAXLJMAX) :: water_fraction, ice_fraction 

 logical, private, parameter :: DEBUG_LU = .false.
 character(len=80), private :: errmsg


contains

 !==========================================================================
  subroutine InitLanduse()

       !=====================================
        call ReadLandUse()     ! => Land_codes,  Percentage cover per grid

        call Init_LandDefs(Land_codes)   ! => LandType, LandDefs
       !=====================================


  end subroutine InitLanduse
 !==========================================================================
  subroutine ReadLanduse()

   integer :: i,j,lu, index_lu, maxlufound
   real, dimension(NLANDUSE) :: tmp
   character(len=20), dimension(NLANDUSE+10) :: Headers
   character(len=(NLANDUSE+10)*20) :: txtinput  ! Big enough to contain one full input record
   type(KeyVal), dimension(10)      :: KeyValues ! Info on units, coords, etc.
   real, dimension(NLANDUSE+1) :: tmpmay
   character(len=50) :: fname
   integer :: iL, n, NHeaders, NKeys, Nlines
   logical :: debug_flag
   real :: sumfrac

  ! Specify the assumed coords and units - Read2DN will check that the data
  ! conform to these.
    type(keyval), dimension(2) :: CheckValues = (/ keyval("Units","PercentGrid"), &
                                                  keyval("Coords","ModelCoords") /)

 ! temporary arrays used.  Will re-write one day....
   real, dimension(MAXLIMAX,MAXLJMAX,NLANDUSE):: landuse_in ! tmp, with all data
   real, dimension(MAXLIMAX,MAXLJMAX,NLUMAX):: landuse_data ! tmp, with all data
   integer, dimension(MAXLIMAX,MAXLJMAX):: landuse_ncodes ! tmp, with all data
   integer, dimension(MAXLIMAX,MAXLJMAX,NLUMAX):: landuse_codes ! tmp, with all data


   if ( DEBUG_LU .and. me == 0 ) write(*,*) "LANDUSE: Starting ReadLandUse, me ",me

   maxlufound = 0   
   Nlines = 0

   landuse_ncodes(:,:)   = 0       !/**  initialise  **/
   landuse_codes(:,:,:)  = 0       !/**  initialise  **/
   landuse_data  (:,:,:) = 0.0     !/**  initialise  **/

!------------------------------------------------------------------------------

      ! Read Header info - this will define landuse classes for model

      fname = "Inputs.Landuse"
      if ( me == 0 ) then
         call open_file(IO_TMP,"r",fname,needed=.true.)
         call CheckStop(ios,"open_file error on " // fname )
      end if


      call Read_Headers(IO_TMP,errmsg,NHeaders,NKeys,Headers,Keyvalues,CheckValues)

      call CheckStop( errmsg , "Read Headers" // fname )
 

      ! The first two columns are assumed for now to be ix,iy, hence:

        NHeaders = NHeaders -2
        call CheckStop( NHeaders /= NLANDUSE, "Inputs.Landuse not consisternt with NLANDUSE")

      ! **** HERE we set the Landuse_codes *****************

        do i = 1, NLANDUSE
           Land_codes(i) = trim ( Headers(i+2) )
        end do

      ! Then data:

      call Read2DN("Inputs.Landuse",NLANDUSE,landuse_in,HeadersRead=.true.)

!------------------------------------------------------------------------------


    if ( DEBUG_LU .and. me == 0 ) then
        write(*,*) "NOW LAND_CODES ARE ", NHeaders
       call WriteArray(Land_codes,NLANDUSE,"Land_Codes")
    end if

    do i = li0, li1
       do j = lj0, lj1
           debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj ) 
           do lu = 1, NLANDUSE
              if ( landuse_in(i,j,lu) > 0.0 ) then

                 call GridAllocate("LANDUSE",i,j,lu,NLUMAX, &
                         index_lu, maxlufound, landuse_codes, landuse_ncodes)
   
                     landuse_data(i,j,index_lu) = &
                       landuse_data(i,j,index_lu) + 0.01 * landuse_in(i,j,lu)
               end if
               if ( DEBUG_LU .and. debug_flag )  &
                       write(*,"(a15,i3,f8.4,a10,i3,f8.4)") "DEBUG Landuse ",&
                          lu, landuse_in(i,j,lu), &
                           "index_lu ", index_lu, landuse_data(i,j,index_lu)
           end do ! lu
           LandCover(i,j)%ncodes  = landuse_ncodes(i,j)
           LandCover(i,j)%codes(:) = landuse_codes(i,j,:)
           LandCover(i,j)%fraction(:)  = landuse_data(i,j,:)

           sumfrac = sum( LandCover(i,j)%fraction(:) )

             if (  sumfrac < 0.99 .or. sumfrac > 1.01 ) then
               write(unit=errmsg,fmt="(a19,3i4,f12.4,8i4)") &
                 "Land SumFrac Error ", me,  &
                    i_fdom(i),j_fdom(j), sumfrac, li0, li1, lj0, lj1, &
                       i_fdom(li0), j_fdom(lj0), i_fdom(li1), j_fdom(lj1)
               call CheckStop(errmsg)
             end if

      end do  !j
   end do  !i

   if (DEBUG_LU) write(6,*) "Landuse_ml: me, Nlines, maxlufound = ", me, Nlines, maxlufound

  end subroutine  ReadLanduse
 
  !=========================================================================
  subroutine  SetLandUse()
    integer :: i,j,ilu,lu, nlu, n ! indices
    logical, save :: my_first_call = .true.
    logical :: debug_flag = .false.
    real :: hveg, lat_factor
    integer :: effectivdaynumber !6 months shift in Southern hemisphere.
    real :: xSAIadd
    logical :: iam_wheat

! Treatment of growing seasons in the southern hemisphere:
!   all the static definitions (SGS,EGS...) refer to northern hemisphere, but the actual 
!   simulation dates are shifted by 6 monthes in the southern hemisphere by using
!   uses effectivdaynumber and mod(current_date%month+5,12)+1 in southern hemis


    if ( DEBUG_LU .and. debug_proc ) write(*,*) "UKDEP SetLandUse, me, day ", me, daynumber, debug_proc
    if ( DEBUG_LU .and. debug_proc ) write(*,*) "DEBUG_LU SetLandUse, me, day ", me, daynumber

    if ( my_first_call ) then
        if ( DEBUG_LU .and. debug_proc ) write(*,*) "UKDEP FIrst Start SetLandUse, me ", me

     ! effectiv daynumber to shift 6 month when in southern hemisphere
        effectivdaynumber=daynumber


      !/ -- Calculate growing seasons where needed and water_fraction
      !          (for Rn emissions)

        water_fraction(:,:) = 0.0         !ds Pb210 
        ice_fraction(:,:)   = 0.0         !ds Pb210 

        do i = li0, li1
          do j = lj0, lj1

             debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj ) 
             do ilu= 1, LandCover(i,j)%ncodes
                lu      = LandCover(i,j)%codes(ilu)
                call CheckStop( lu < 0 .or. lu > NLANDUSE , &
                                "SetLandUse out of range" )

                if ( LandDefs(lu)%SGS50 > 0 ) then  ! need to set growing seasons 

                    call Growing_season( lu,abs(gb(i,j)),&  
                            LandCover(i,j)%SGS(ilu),LandCover(i,j)%EGS(ilu) )
                    if ( DEBUG_LU .and. debug_flag ) write(*,*)"LU_SETGS", lu,  LandCover(i,j)%SGS(ilu),LandCover(i,j)%EGS(ilu)
                else
                   LandCover(i,j)%SGS(ilu) =  LandDefs(lu)%SGS50
                   LandCover(i,j)%EGS(ilu) =  LandDefs(lu)%EGS50
                    if ( DEBUG_LU .and. debug_flag ) write(*,*)"LU_FIXGS", lu,  LandCover(i,j)%SGS(ilu),LandCover(i,j)%EGS(ilu)
                end if


               !/ for landuse classes with bulk-resistances, we only
               !  need to specify height once. Dummy values are assigned
               !  to LAI and gpot:

             if ( LandType(lu)%is_bulk ) then
                 LandCover(i,j)%hveg(ilu) =  LandDefs(lu)%hveg_max
                 LandCover(i,j)%LAI(ilu)  =  0.0          
                 LandCover(i,j)%fphen(ilu) =  0.0          
             end if

             if ( LandType(lu)%is_water ) water_fraction(i,j) = LandCover(i,j)%fraction(ilu)
             if ( LandType(lu)%is_ice   )   ice_fraction(i,j) = LandCover(i,j)%fraction(ilu)

                if ( DEBUG_LU .and. debug_flag ) then
                      write(*,"(a,2i4,2f12.4)") "DEBUG_LU WATER ", ilu, lu, &
                          water_fraction(i,j), ice_fraction(i,j)
                end if

            end do ! ilu
          end do ! j
        end do ! i

        my_first_call = .false.
   !======================================================================
    end if ! my_first_call
   !======================================================================

    
     do i = li0, li1
       do j = lj0, lj1

          effectivdaynumber=daynumber
         ! effectiv daynumber to shift 6 months when in southern hemisphere
          if(gb(i,j)<0.0)effectivdaynumber=mod(daynumber+182,nydays)+1 

          debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj ) 
          if ( DEBUG_LU .and. debug_flag ) then
                 write(*,"(a12,i3,i4)") "LANDUSE N Day? ", LandCover(i,j)%ncodes, daynumber
          end if
          do ilu= 1, LandCover(i,j)%ncodes
             lu      = LandCover(i,j)%codes(ilu)

             if ( LandType(lu)%is_bulk ) cycle    !else Growing veg present:

             LandCover(i,j)%LAI(ilu) = Polygon(effectivdaynumber, &
                                      0.0, LandDefs(lu)%LAImin, LandDefs(lu)%LAImax,&
                                      LandCover(i,j)%SGS(ilu), LandDefs(lu)%SLAIlen, &
                                      LandCover(i,j)%EGS(ilu), LandDefs(lu)%ELAIlen)

             LandCover(i,j)%fphen(ilu) = fPhenology( lu, LandDefs(lu)%code,effectivdaynumber &
                              ,LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu)  &
                              ,debug_flag )



             hveg = LandDefs(lu)%hveg_max   ! defaults
             xSAIadd = 0.0


             iam_wheat = .false.
             if (  LandType(lu)%is_crop ) then

                if ( LandType(lu)%is_iam  ) then ! IAM wheat
                    iam_wheat = .true.
                    if  ( effectivdaynumber >= LandCover(i,j)%SGS(ilu) .and. &
                          effectivdaynumber <= LandCover(i,j)%EGS(ilu)  ) then
                            WheatGrowingSeason(i,j) =  1
                    else
                            WheatGrowingSeason(i,j) =  0
                    end if
                end if

                if ( effectivdaynumber < LandCover(i,j)%SGS(ilu) .or. &
                     effectivdaynumber > LandCover(i,j)%EGS(ilu)  ) then
                   hveg = STUBBLE
                   xSAIadd = 0.0
                else if ( effectivdaynumber < &
                     (LandCover(i,j)%SGS(ilu) + LandDefs(lu)%SLAIlen) ) then
                   hveg=  LandDefs(lu)%hveg_max * &
                     LandCover(i,j)%LAI(ilu) / LandDefs(lu)%LAImax
                   xSAIadd = ( 5.0/3.5 - 1.0) * LandCover(i,j)%LAI(ilu)
                else if ( effectivdaynumber < LandCover(i,j)%EGS(ilu) ) then
                   hveg = LandDefs(lu)%hveg_max             ! not needed?
                   xSAIadd = 1.5   ! Sensescent
                end if
                LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)  + xSAIadd

             !! end if ! crops


           ! Just used reduced LAI for high latitudes for now, because of tests
           ! which suggest that the big-leaf model as coded will overestimate
           ! Gsto if we allow higher LAI in central Europe.

             else if( LandType(lu)%is_forest ) then
               if ( gb(i,j) >= 60.0 ) then
                       lat_factor  = max(0.3, ( 1.0 - 0.05* (gb(i,j)-60.0)) )
                       hveg  = hveg *  lat_factor
                       LandCover(i,j)%LAI(ilu) = LandCover(i,j)%LAI(ilu)  * lat_factor
               end if
               LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)  + 1.0
             else
               LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)   !defaults
             end if

             LandCover(i,j)%hveg(ilu) =  hveg

            if ( DEBUG_LU .and. debug_flag ) then
                   write(*,"(a12,i3,a16,i4,f7.2,2f8.3,4i4)") "LANDPhen ", lu, trim(LandDefs(lu)%name), daynumber, &
                     LandCover(i,j)%hveg(ilu), LandCover(i,j)%LAI(ilu), LandCover(i,j)%fphen(ilu), &
                     LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu)
            end if
                   

         end do ! lu
       end do ! j
    end do ! i
    if ( DEBUG_LU .and. me==0 ) write(*,*)"UKDEP Finishing SetLandUse "
    if(debug_proc .and. DEBUG_LU) write(*,*) "LAST GROWSEASON ", effectivdaynumber, WheatGrowingSeason(debug_li,debug_lj)

  end subroutine  SetLandUse
! =====================================================================

!=======================================================================
function Polygon(jdayin,Ymin,Ystart,Ymax,Sday,LenS,Eday,LenE) &
result (Poly)
!=======================================================================

!     Calculates the value of a parameter Y with a polygon
!     distribution - currently LAI and g_pot

!            _____________       <- Ymax
!           /             \
!          /               \
!         /                 \
!        /                   \
!       |                     |  <- Ystart
!       |                     |
!       |                     |
!  ----------------------------- <- Ymin
!       S  S1            E1   E
!
!1.4 The following has been simplified
    

!   Inputs
    integer, intent(in) :: jdayin      !day of year
!d1.4 integer, intent(in) :: yydays    !no. days in year (365 or 366)
    real, intent(in) ::    Ymin        !minimum value of Y
    real, intent(in) ::    Ystart      !value Y at start of growing season
    real, intent(in) ::    Ymax        !maximum value of Y
    integer, intent(in) ::    Sday        !start day (e.g. of growing season)
    integer, intent(in) ::    LenS        !length of Start period (S..S1 above)
    integer, intent(in) ::    Eday        !end day (e.g. of growing season)
    integer, intent(in) ::    LenE        !length of end period (E..E1 above)

!  Output:
    real ::   Poly  ! value at day jday

! Local
    integer :: jday  ! day of year, after any co-ordinate change
    integer ::    S   !  start day
    integer ::    E   !  end day
    

    jday = jdayin
    E = Eday
    S = Sday

  ! Here we removed a lot of code associated with the leaf-age
  ! version of g_pot. 
       
    if ( jday  <  S .or. jday >  E ) then
       Poly = Ymin
       return
    end if

  !d1.3 - slightly re-written tests:

    if (jday <=  S+LenS  .and. LenS > 0 ) then

        Poly = (Ymax-Ystart) * (jday-S)/LenS  + Ystart 

    else if ( jday >=  E-LenE .and. LenE > 0.0 ) then   !d1.1 test for LenE

        Poly = (Ymax-Ystart) * (E-jday)/LenE + Ystart

    else

        Poly =Ymax
    end if
    

 end function Polygon

 !=======================================================================
end module Landuse_ml
