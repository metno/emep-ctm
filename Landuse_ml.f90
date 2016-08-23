! <Landuse_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!***************************************************************************! 
!* 
!*  Copyright (C) 2007-2012 met.no
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
!***************************************************************************! 
module Landuse_ml

use CheckStop_ml,   only: CheckStop,StopAll
use DO3SE_ml,       only: fPhenology
use GridAllocate_ml,only: GridAllocate
use GridValues_ml,  only: glat_fdom, glat    & ! latitude,
                          , i_fdom, j_fdom   & ! coordinates
                          , i_local, j_local &
                          , debug_proc, debug_li, debug_lj
use Io_ml,          only: open_file, ios, Read_Headers, Read2DN, IO_TMP
use KeyValue_ml,    only: KeyVal,KeyValue, LENKEYVAL
use LandDefs_ml,    only: Init_LandDefs, LandType, LandDefs, &
                          STUBBLE, Growing_Season,&
                          NLanduse_DEF,NLANDUSE_EMEP
use LandPFT_ml,       only: MapPFT_LAI, pft_lai
use MetFields_ml,     only: nwp_sea ,foundnwp_sea
use ModelConstants_ml,only: DEBUG_i, DEBUG_j, NLANDUSEMAX, &
                            USE_PFT_MAPS, DEBUG_LANDPFTS, &
                            DEBUG_LANDUSE, NPROC, IIFULLDOM, JJFULLDOM, &
                            DomainName, MasterProc
use NetCDF_ml,      only: ReadField_CDF,printcdf
use Par_ml,         only: MAXLIMAX, MAXLJMAX, &
                          limax, ljmax, me
use SmallUtils_ml,  only: find_index, NOT_FOUND, WriteArray
use TimeDate_ml,    only: daynumber, effectivdaynumber, nydays, current_date

implicit none
private


!/- subroutines:

  public :: InitLanduse
  public :: ReadLanduse
  public :: SetLanduse
  private :: Polygon         ! Used for LAI
 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO

 integer, public, parameter :: NLUMAX = 19 ! max no. landuse per grid
 !ds logical, public :: LU_cdf

! The headers read from Inputs.Landuse define the "master-list" of
! codes for landuse. Each code must be present in the subsequent
! data files for phenology and DO3SE.

 character(len=6), dimension(NLANDUSEMAX), public, save :: Land_codes ! As used

 !=============================================
 type, public :: LandCov
   integer                   :: ncodes     ! Number of codes in grid
   integer,dimension(NLUMAX) :: &
          codes     &! landcover codes
         ,SGS       &! Start of growing season (days)
         ,EGS       &! End of growing season (days)
         ,Astart    &! Start photosynthetic activity, for DO3SE
         ,Aend       ! 
   real,   dimension(NLUMAX) :: &
          fraction  &! (coverage)
         ,LAI       &! Leaf-area-index (m2/m2)
         ,SAI       &! Surface-area-index (m2/m2) (leaves+bark, etc.)
         ,hveg      &! Max. height of veg.
         ,fphen     &! Potential (age) factor for Jarvis-calc
         ,Eiso      &! Emission potential, isoprene
         ,Emt       &! Emission potential, monoterpenes
         ,SumVPD    &! For critical VPD calcs, reset each day
         ,old_gsun   ! also for flux
 end type LandCov
 !=============================================
 type(LandCov), public, save, allocatable,dimension(:,:) :: LandCover
 !=============================================


  logical, public,save, allocatable,dimension(:,:) :: likely_coastal 

  integer, public,save, allocatable,dimension(:,:) :: &
          WheatGrowingSeason  ! Growing season (days), IAM_WHEAT =1 for true
 
 ! For some flux work, experimental

 real,public,save, allocatable,dimension(:,:) :: water_fraction, ice_landcover 
 logical,public,save :: water_frac_set = .false.

 character(len=80), private :: errmsg


contains

 !==========================================================================
  subroutine InitLanduse()
    logical :: filefound
    real, parameter ::water_fraction_THRESHOLD=0.5
    integer ::i,j,ilu,lu
    logical :: debug_flag = .false.
    !=====================================

    !ALLOCATE ARRAYS
    allocate(LandCover(MAXLIMAX,MAXLJMAX))
    allocate(likely_coastal(MAXLIMAX,MAXLJMAX) )
    likely_coastal  = .false.
    allocate(WheatGrowingSeason(MAXLIMAX,MAXLJMAX))
    allocate(water_fraction(MAXLIMAX,MAXLJMAX), ice_landcover(MAXLIMAX,MAXLJMAX))


    !ReadLandUse_CDF to be used as default when glc2000 data is improved?


    filefound=.false.
    call ReadLandUse(filefound) !=> Land_codes, Percentage cover per grid

    !ReadLandUse_CDF use Max Posch 5km landuse over emep area and glc200 where this dat is not defined.
    if(.not.filefound)call ReadLandUse_CDF(filefound) !=> Land_codes, Percentage cover per grid

    call CheckStop(.not.filefound,"InitLanduse failed!")

    call Init_LandDefs(Land_codes)   ! => LandType, LandDefs


    ! effectiv daynumber to shift 6 month when in southern hemisphere
    effectivdaynumber=daynumber


    !/ -- Calculate growing seasons where needed and water_fraction
    !          (for Rn emissions)

    water_fraction(:,:) = 0.0    !for Pb210 
    ice_landcover(:,:) = 0.0  !for Pb210 
    likely_coastal(:,:) = .false.   ! already done, but just for clarity

    do i = 1, limax
       do j = 1, ljmax

          debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj )
          do ilu= 1, LandCover(i,j)%ncodes
             lu      = LandCover(i,j)%codes(ilu)
             call CheckStop( lu < 0 .or. lu > NLANDUSEMAX , &
                  "SetLandUse out of range" )

             if ( LandDefs(lu)%SGS50 > 0 ) then ! need to set growing seasons 

                call Growing_season( lu,abs(glat(i,j)),&  
                     LandCover(i,j)%SGS(ilu),LandCover(i,j)%EGS(ilu) )
             else
                LandCover(i,j)%SGS(ilu) =  LandDefs(lu)%SGS50
                LandCover(i,j)%EGS(ilu) =  LandDefs(lu)%EGS50
             end if
             if ( DEBUG_LANDUSE .and. debug_flag ) &
                  write(*,"(a,i3,a20,2i4)")"LANDUSE: LU_SETGS", &
                  lu, LandDefs(lu)%name,&
                  LandCover(i,j)%SGS(ilu),LandCover(i,j)%EGS(ilu)


             !/ for landuse classes with bulk-resistances, we only
             !  need to specify height once. Dummy values are assigned
             !  to LAI and gpot:

             if ( LandType(lu)%is_bulk ) then
                LandCover(i,j)%hveg(ilu) =  LandDefs(lu)%hveg_max
                LandCover(i,j)%LAI(ilu)  =  0.0          
                LandCover(i,j)%fphen(ilu) =  0.0          
             end if

             if ( LandType(lu)%is_water ) water_fraction(i,j) = &
                  LandCover(i,j)%fraction(ilu)
             if ( LandType(lu)%is_ice   ) ice_landcover(i,j) = &
                  LandCover(i,j)%fraction(ilu)


          end do ! ilu

!set default values for nwp_sea
          if(water_fraction(i,j)>water_fraction_THRESHOLD) &
               nwp_sea(i,j) = .true.


          ! We don't want to trust some squares with a mixture of sea
          ! and land for micromet purposes, e.g. T2 can be very wrong
          ! We mark these as likely coastal:
          if ( nwp_sea(i,j) )  then
             if (  water_fraction(i,j) < 1.0  ) likely_coastal(i,j) = .true.
          else if ( water_fraction(i,j) > 0.2 ) then
             likely_coastal(i,j) = .true.
          end if !
       end do ! j
    end do ! i

    water_frac_set = .true.  ! just to inform other routines
    
  end subroutine InitLanduse
 !==========================================================================
  subroutine ReadLanduse(filefound)

   logical :: filefound
   integer :: i,j,lu, index_lu, maxlufound
   character(len=20), dimension(NLANDUSEMAX+10) :: Headers
   type(KeyVal), dimension(10)      :: KeyValues ! Info on units, coords, etc.
   character(len=50) :: fname
   integer :: NHeaders, NKeys, Nlines
   logical :: debug_flag
   real :: sumfrac
   
  ! Specify the assumed coords and units - Read2DN will check that the data
  ! conform to these.
    type(keyval), dimension(2) :: CheckValues = &
        (/ keyval("Units","PercentGrid"), &
           keyval("Coords","ModelCoords") /)

 ! temporary arrays used.  Will re-write one day....
   real, dimension(MAXLIMAX,MAXLJMAX,NLANDUSEMAX):: landuse_in ! tmp, with all data
   real, dimension(MAXLIMAX,MAXLJMAX,NLUMAX):: landuse_data ! tmp, with all data
   integer, dimension(MAXLIMAX,MAXLJMAX):: landuse_ncodes ! tmp, with all data
   integer, dimension(MAXLIMAX,MAXLJMAX,NLUMAX):: landuse_codes ! tmp, with all data

   if ( DEBUG_LANDUSE .and. MasterProc ) &
        write(*,*) "LANDUSE: Starting ReadLandUse "

   maxlufound = 0   
   Nlines = 0

   landuse_ncodes(:,:)   = 0     !/**  initialise  **/
   landuse_codes(:,:,:)  = 0     !/**  initialise  **/
   landuse_data  (:,:,:) = 0.0   !/**  initialise  **/

!------------------------------------------------------------------------------

      ! Read Header info - this will define landuse classes for model

      fname = "Inputs.Landuse"
      if ( MasterProc ) then
!         call open_file(IO_TMP,"r",fname,needed=.true.)
!         call CheckStop(ios,"open_file error on " // fname )
         call open_file(IO_TMP,"r",fname,needed=.false.)
      end if
      call MPI_BCAST( ios, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,INFO)
      if(ios==0)then
         filefound=.true.

         call Read_Headers(IO_TMP,errmsg,NHeaders,NKeys,&
                     Headers,Keyvalues,CheckValues)
         
         call CheckStop( errmsg , "Read Headers" // fname )
         
         ! The first two columns are assumed for now to be ix,iy, hence:

         NHeaders = NHeaders -2
         call CheckStop( NHeaders /= NLANDUSE_EMEP, &
              "Inputs.Landuse not consistent with NLANDUSE_EMEP")

          NLanduse_DEF=NHeaders        
        
         ! *** HERE we set the Landuse_codes ***
         do i = 1,  NLanduse_DEF
            Land_codes(i) = trim ( Headers(i+2) )
         end do
         if(MasterProc)write(*,*)NLanduse_DEF,' landuse categories defined from Inputs.Landuse:'
         if(MasterProc)write(*,fmt="(20(A,1x))")(trim(Land_codes(i)),i=1,NLanduse_DEF)
         ! Then data:
         
         call Read2DN("Inputs.Landuse",NLanduse_DEF,landuse_in,&
                 HeadersRead=.true.)
         
         !-------------------------------------------------------------------
         
         if ( DEBUG_LANDUSE .and. MasterProc ) then
            write(*,*) "LANDUSE: LAND_CODES ARE ", NHeaders
            call WriteArray(Land_codes,NLanduse_DEF,"Land_Codes")
         end if
         
      else
         filefound=.false.
         if(MasterProc)Write(*,*)'Inputs.Landuse not found'
         return
         call StopAll('Inputs.Landuse not found') 
      endif

!      call printCDF('LU', landuse_in(:,:,1),'??')

    do i = 1, limax
       do j = 1, ljmax
           debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj ) 
           do lu = 1, NLanduse_DEF
              if ( landuse_in(i,j,lu) > 0.0 ) then

                 call GridAllocate("LANDUSE",i,j,lu,NLUMAX, &
                         index_lu, maxlufound, landuse_codes, landuse_ncodes)
   
                     landuse_data(i,j,index_lu) = &
                       landuse_data(i,j,index_lu) + 0.01 * landuse_in(i,j,lu)
               end if
               if ( DEBUG_LANDUSE .and. debug_flag )  &
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
                    i_fdom(i),j_fdom(j), sumfrac, limax,  ljmax, &
                       i_fdom(1), j_fdom(1), i_fdom(limax), j_fdom(ljmax)
               call CheckStop(errmsg)
             end if

      end do  !j
   end do  !i

   if (DEBUG_LANDUSE) write(6,*) "Landuse_ml: me, Nlines, maxlufound = ", &
                                  me, Nlines, maxlufound

  end subroutine  ReadLanduse
 
  subroutine ReadLanduse_CDF(filefound)
    !Read data in other grid and interpolate to present grid
    !
    !So far only basic version for use in TNO7. Under construction
    !
    implicit none
    logical :: filefound
    integer :: i,j,lu, index_lu, maxlufound
    logical :: debug_flag
    real :: sumfrac


    ! temporary arrays used.  Will re-write one day....
    real, dimension(MAXLIMAX,MAXLJMAX,NLANDUSEMAX):: landuse_in ! tmp, with all data
    real, dimension(MAXLIMAX,MAXLJMAX):: landuse_tmp ! tmp, with all data
    real, dimension(MAXLIMAX,MAXLJMAX,NLUMAX):: landuse_data ! tmp, with all data
    integer, dimension(MAXLIMAX,MAXLJMAX):: landuse_ncodes ! tmp, with all data
    integer, dimension(MAXLIMAX,MAXLJMAX,NLUMAX):: landuse_codes ! tmp, with all data

    if ( DEBUG_LANDUSE .and. MasterProc ) &
         write(*,*) "LANDUSE: Starting ReadLandUse "


    if (MasterProc ) write(*,*) "LANDUSE_CDF: experimental so far"
!    filefound=.false.
!    return

    maxlufound = 0   

    landuse_ncodes(:,:)   = 0     !/**  initialise  **/
    landuse_codes(:,:,:)  = 0     !/**  initialise  **/
    landuse_data  (:,:,:) = 0.0   !/**  initialise  **/

    !hardcoded so far -> to softimize

    NLanduse_DEF=19
    Land_codes(1) = 'CF' 
    Land_codes(2) = 'DF' 
    Land_codes(3) = 'NF' 
    Land_codes(4) = 'BF' 
    Land_codes(5) = 'TC' 
    Land_codes(6) = 'MC' 
    Land_codes(7) = 'RC' 
    Land_codes(8) = 'SNL' 
    Land_codes(9) = 'GR' 
    Land_codes(10) = 'MS' 
    Land_codes(11) = 'WE' 
    Land_codes(12) = 'TU' 
    Land_codes(13) = 'DE' 
    Land_codes(14) = 'W' 
    Land_codes(15) = 'ICE' 
    Land_codes(16) = 'U' 
    Land_codes(17) = 'IAM_CR' 
    Land_codes(18) = 'IAM_DF'
    Land_codes(19) = 'IAM_MF'
    do lu=1,NLanduse_DEF
       !
       if(me==0)write(*,*)'Reading landuse ',trim(Land_codes(lu))
       !   call ReadField_CDF('/global/work/mifapw/emep/Data/LanduseGLC.nc',&!fast but unprecise
       call ReadField_CDF('Landuse_PS_5km.nc',& !SLOW!
            Land_codes(lu),landuse_in(1,1,lu),1,interpol='conservative', &
            needed=.true.,debug_flag=.false.,UnDef=-9.9E19) !NB: Undef must be largenegative, 
!          because it is averagad over many points, and the final result must still be negative
          call ReadField_CDF('LanduseGLC.nc',&
               Land_codes(lu),landuse_tmp,1,interpol='conservative', &
               needed=.true.,debug_flag=.false.)
          do j = 1, ljmax
             do i = 1, limax
                if(landuse_in(i,j,lu)<-0.1)landuse_in(i,j,lu)=landuse_tmp(i,j)
             end do  !j
          end do  !i
    enddo
     ! call printCDF('LU_cdf', landuse_in(:,:,1),'??')

    do i = 1, limax
       do j = 1, ljmax
          do lu = 1, NLanduse_DEF
             if ( landuse_in(i,j,lu) > 0.0 ) then

                call GridAllocate("LANDUSE",i,j,lu,NLUMAX, &
                     index_lu, maxlufound, landuse_codes, landuse_ncodes)
                landuse_data(i,j,index_lu) = &
                     landuse_data(i,j,index_lu) + landuse_in(i,j,lu)!already in fraction unit
             endif
          end do ! lu
          LandCover(i,j)%ncodes  = landuse_ncodes(i,j)
          LandCover(i,j)%codes(:) = landuse_codes(i,j,:)
          LandCover(i,j)%fraction(:)  = landuse_data(i,j,:)
          sumfrac = sum( LandCover(i,j)%fraction(:) )

          if (  sumfrac < 0.99 .or. sumfrac > 1.01 ) then
             write(unit=errmsg,fmt="(a19,3i4,f12.4,8i4)") &
                  "Land SumFrac Error ", me,  &
                  i_fdom(i),j_fdom(j), sumfrac, limax,  ljmax, &
                  i_fdom(1), j_fdom(1), i_fdom(limax), j_fdom(ljmax)
             call CheckStop(errmsg)
          end if

       end do  !j
    end do  !i

    filefound=.true.

  end subroutine ReadLanduse_CDF

  !=========================================================================
  subroutine  SetLandUse()
    integer :: i,j,ilu,lu ! indices
    integer, save :: old_month = -1
    integer, save :: old_daynumber = -1
    logical, save :: my_first_call = .true.
    logical :: debug_flag = .false.
    real :: hveg, lat_factor
    real :: xSAIadd
     integer :: pft

! Treatment of growing seasons in the southern hemisphere:
!   all the static definitions (SGS,EGS...) refer to northern hemisphere, 
!   but the actual simulation dates are shifted by 6 months in the southern
!   hemisphere by using uses effectivdaynumber and 
!   mod(current_date%month+5,12)+1 in southern hemis


    if ( DEBUG_LANDUSE .and. debug_proc ) then
        write(*,*) "LANDUSE: SetLandUse, me, day ", me, daynumber, debug_proc
    end if

   !======================================================================
    if ( my_first_call ) then
       !read in data from file
        my_first_call   = .false.
        
        call InitLanduse()

    end if ! my_first_call
   !======================================================================

    if ( daynumber == old_daynumber ) then
        return
    end if
    old_daynumber = daynumber


   !Landcover data can be set either from simplified LPJ
   !PFTs, or from the "older" DO3SE inputs file

     if ( USE_PFT_MAPS ) then !- Check for LPJ-derived data -
         if ( current_date%month /= old_month ) then 
           call MapPFT_LAI( current_date%month )
         end if
     end if

    
     do i = 1, limax
       do j = 1, ljmax

          effectivdaynumber=daynumber
         ! effectiv daynumber to shift 6 months when in southern hemisphere
          if(glat(i,j)<0.0)effectivdaynumber=mod(daynumber+182,nydays)+1 

          debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj ) 
          if ( DEBUG_LANDUSE .and. debug_flag ) then
                 write(*,"(a12,i3,i4)") "LANDUSE N Day? ", &
                  LandCover(i,j)%ncodes, daynumber
                 write(*,*) "LANDUSE DATE ", current_date
          end if
          do ilu= 1, LandCover(i,j)%ncodes
             lu      = LandCover(i,j)%codes(ilu)
             pft     = LandType(lu)%pft

             if ( LandType(lu)%is_bulk ) then
                LandCover(i,j)%LAI(ilu) = 0.0
                LandCover(i,j)%SAI(ilu) = 0.0
                cycle    
             endif!else Growing veg present:

             LandCover(i,j)%LAI(ilu) = Polygon(effectivdaynumber, &
                    0.0, LandDefs(lu)%LAImin, LandDefs(lu)%LAImax,&
                      LandCover(i,j)%SGS(ilu), LandDefs(lu)%SLAIlen, &
                         LandCover(i,j)%EGS(ilu), LandDefs(lu)%ELAIlen)

             LandCover(i,j)%fphen(ilu) = fPhenology( lu &
                ,effectivdaynumber &
                ,LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu)&
                ,debug_flag )

          if ( DEBUG_LANDPFTS .and. debug_flag.and. USE_PFT_MAPS ) then
                 if ( pft > 0.0 ) then
                   write(*,"(2a,i4,i6,2f8.3)") "LANDUSE PFTS COMP? ", &
                      LandDefs(lu)%name, daynumber, pft,&
                       LandCover(i,j)%LAI(ilu), pft_lai(i,j, pft)
                 end if
          end if


             hveg = LandDefs(lu)%hveg_max   ! defaults
             xSAIadd = 0.0

             if (  LandType(lu)%is_crop ) then

                if ( LandType(lu)%is_iam  ) then ! IAM wheat
                    if  ( effectivdaynumber >= LandCover(i,j)%SGS(ilu) .and. &
                          effectivdaynumber <= LandCover(i,j)%EGS(ilu)  ) then
                            WheatGrowingSeason(i,j) =  1
                    else
                            WheatGrowingSeason(i,j) =  0
                    end if
                end if

               ! Note that IAM crops have SLAIlen=0, so are immediately
               ! given LAI=3.5, SAI=5.

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
                   hveg = LandDefs(lu)%hveg_max  ! not needed?
                   xSAIadd = 1.5  ! Sensescent
                end if
                LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)  + xSAIadd

             ! end if ! crops


           ! Just used reduced LAI for high latitudes for now, because of tests
           ! which suggest that the big-leaf model as coded will overestimate
           ! Gsto if we allow higher LAI in central Europe.

             else if( LandType(lu)%is_forest ) then
               if ( glat(i,j) >= 60.0 ) then
                       lat_factor  = max(0.3, ( 1.0 - 0.05* (glat(i,j)-60.0)) )
                       hveg  = hveg *  lat_factor
                       LandCover(i,j)%LAI(ilu) = LandCover(i,j)%LAI(ilu)  * lat_factor
               end if
               LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)  + 1.0
             else
               LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)   !defaults
             end if

             LandCover(i,j)%hveg(ilu) =  hveg


         end do ! lu
       end do ! j
    end do ! i

! --- print out for debug cell
    if ( DEBUG_LANDUSE.and.debug_proc ) then
       i=debug_li
       j=debug_lj

       do ilu= 1, LandCover(i,j)%ncodes
          lu      = LandCover(i,j)%codes(ilu)
          pft     = LandType(lu)%pft
          if ( LandType(lu)%is_bulk ) cycle    !else Growing veg present:

            write(*,"(a,i3,a16,i4,f7.2,3f8.3,2i4)") "LANDUSE Phen ", lu,&
             trim(LandDefs(lu)%name), daynumber, LandCover(i,j)%hveg(ilu),&
              LandCover(i,j)%SAI(ilu), LandCover(i,j)%LAI(ilu), &
              LandCover(i,j)%fphen(ilu), &
             LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu)
       end do

   end if
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

!   Inputs
    integer, intent(in) :: jdayin     !day of year
    real, intent(in) ::    Ymin       !minimum value of Y
    real, intent(in) ::    Ystart     !value Y at start of growing season
    real, intent(in) ::    Ymax       !maximum value of Y
    integer, intent(in) ::    Sday    !start day (e.g. of growing season)
    integer, intent(in) ::    LenS    !length of Start period (S..S1 above)
    integer, intent(in) ::    Eday    !end day (e.g. of growing season)
    integer, intent(in) ::    LenE    !length of end period (E..E1 above)

!  Output:
    real ::   Poly  ! value at day jday

! Local
    integer :: jday ! day of year, after any co-ordinate change
    integer ::    S ! start day
    integer ::    E ! end day
    
    jday = jdayin
    E = Eday
    S = Sday

  ! Here we removed a lot of code associated with the leaf-age
  ! version of g_pot. 
       
    if ( jday  <  S .or. jday >  E ) then
       Poly = Ymin
       return
    end if


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
