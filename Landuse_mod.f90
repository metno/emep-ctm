! <Landuse_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
!> <Landuse_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!! ************************************************************************! 

module Landuse_mod

use CheckStop_mod,   only: CheckStop,StopAll
use Config_module,   only: NLANDUSEMAX, SEA_LIMIT, USES,  &
                            FLUX_VEGS, FLUX_IGNORE,  nFluxVegs, & 
                            VEG_2dGS, VEG_2dGS_Params, & 
                            NPROC, IIFULLDOM, JJFULLDOM, &
                            MasterProc, LandCoverInputs 
use Debug_module,    only: DEBUG   ! -> DEBUG%LANDUSE
use DO3SE_mod,       only: fPhenology, Init_DO3SE
use GridAllocate_mod,only: GridAllocate
use GridValues_mod,  only:  glat , glon   & ! latitude,
                          , i_fdom, j_fdom   & ! coordinates
                          , i_local, j_local &
                          , debug_proc, debug_li, debug_lj
use Io_mod,          only: open_file, ios, Read_Headers, Read2DN, IO_TMP &
                         ,IO_DO3SE
use KeyValueTypes,    only: KeyVal,KeyValue, LENKEYVAL
use LandDefs_mod,    only: Init_LandDefs, LandType, LandDefs, &
                          STUBBLE, Growing_Season,&
                          NLANDUSE_EMEP
use LandPFT_mod,       only: MapPFT_LAI, pft_lai
use MPI_Groups_mod, only : MPI_INTEGER,MPI_COMM_CALC, IERROR
use Par_mod,         only: LIMAX, LJMAX, &
                          limax, ljmax, me
use SmallUtils_mod,  only: wordsplit, find_index, NOT_FOUND, WriteArray, trims
use TimeDate_mod,    only: effectivdaynumber, nydays, current_date

use netcdf
use NetCDF_mod, only  : ReadField_CDF,check,printcdf

implicit none
private


!/- subroutines:

  public :: InitLanduse
  public :: SetLanduse
  private :: Polygon         ! Used for LAI
  private :: MedLAI          ! Used for LAI, Medit.

 integer, public, parameter :: NLUMAX = 30 ! max no. landuse per grid
 integer, private, save :: NLand_codes = 0 ! no. landuse in input files

! The LC: entries in the netcdf landuse, or the headers read from 
! Inputs.Landuse define the "master-list" of codes for landuse. Each code must
! be present in the subsequent data files for phenology and DO3SE.

 character(len=20), dimension(NLANDUSEMAX), &
           public, save :: Land_codes = " " ! As used

 !=============================================
 type, public :: LandCov
   integer                   :: ncodes     ! Number of codes in grid
   integer,dimension(NLUMAX) :: &
          codes     &! landcover codes
         ,SGS       &! Start of growing season (days)
         ,EGS       &! End of growing season (days)
         ,Anth      &! Anthesis. For wheat
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
 logical, public,save, allocatable,dimension(:,:) :: mainly_sea 

 real,public,save, allocatable,dimension(:,:) :: water_fraction, ice_landcover 
 logical,public,save :: water_frac_set = .false.

 character(len=200), private :: errmsg

contains

 !==========================================================================
  subroutine InitLanduse(daynumber)
    integer, intent(in) :: daynumber
    logical :: filefound
    integer ::i,j,ilu,lu, ipar
    logical :: dbgij
    ! Some config options, for veg with 2-D growing seasons
    ! needs map2d.. as real for ReadCDF
    integer :: n2dGS, n2dGSpars, i2dGS
    real,dimension(:,:,:),allocatable :: map2dGrowingSeasons
    character(len=len(VEG_2dGS_Params(1))) :: fname
    character(len=20) :: varname
    character(len=*), parameter :: dtxt='InitLanduse:'
    !=====================================

    !ALLOCATE ARRAYS
    allocate(LandCover(LIMAX,LJMAX))
    allocate(likely_coastal(LIMAX,LJMAX) )
    allocate(mainly_sea(LIMAX,LJMAX) )
    allocate(water_fraction(LIMAX,LJMAX), ice_landcover(LIMAX,LJMAX))

    ! First, check the number of "extra" (fake) vegetation 
    nFluxVegs = 0
    do ilu = 1, size( FLUX_VEGS )
        if(len_trim(FLUX_VEGS(ilu))>0) nFluxVegs=nFluxVegs+1
    end do

    if(MasterProc) write(*,*) dtxt//" nFluxVegs= ",nFluxVegs

    call ReadLandUse_CDF(filefound) !=> Land_codes, % cover per grid

    ! Quick safety check
    ! we check that the length of the land-codes isn't equal to our declared
    ! length. That is usualyl a sign that the codes are too long and we may
    ! get into truncation worries.
    call CheckStop(maxval( len_trim(Land_codes(:))) >= len(Land_codes(1)),&
          "Land_Codes: increase size of character array" )

    call CheckStop(.not.filefound,"InitLanduse failed!")

    if(MasterProc) then
        write(*,*)  dtxt//" Into Init_LandDefs ", NLand_codes
        write(*,*)  dtxt//" Codes: ", Land_codes
        write(*,*)  dtxt//" LandCoverInputs: ", LandCoverInputs
    end if

    call Init_LandDefs(LandCoverInputs%LandDefs,NLand_codes, Land_codes)  
        ! => LandType, LandDefs

    !------ 2D maps of growing season, if set in config ----------------------- 

     n2dGS    = count( VEG_2dGS(:) /= "" )

    if ( n2dGS > 0 )  then

       n2dGSpars = count( VEG_2dGS_Params(:) /= "" ) - 1 ! First is fname

       allocate(map2dGrowingSeasons(LIMAX,LJMAX,n2dGS*n2dGSpars))

       fname = VEG_2dGS_Params(1)
       if(MasterProc) write(*,*) "CRU GS MAP ", trim(fname), n2dGS, n2dGSpars

       i2dGS = 0 
       do ilu = 1, n2dGS
         do ipar = 2, n2dGSpars+1

           varname = trims (&
               VEG_2dGS(ilu) // "_" // VEG_2dGS_Params(ipar) // ".txt" )

           i2dGS = i2dGS + 1 

           call ReadField_CDF(fname,varname, map2dGrowingSeasons(1,1,i2dGS),&
             1,interpol='zero_order',needed=.true.,debug_flag=.false.)  

           if( debug_proc ) write(*,"(a20,5i5)") '2dGS '//trim(varname),ilu,&
              debug_li,debug_lj, i2dGS, &
               nint( map2dGrowingSeasons(debug_li,debug_lj,i2dGS) )

         end do ! ipar
       end do ! ilu
    end if
    !------  End 2D maps of growing season ------------------------------------- 


    ! effectiv daynumber to shift 6 month when in southern hemisphere
    effectivdaynumber=daynumber


    !/ -- Calculate growing seasons where needed and water_fraction
    !          (for Rn emissions)

    water_fraction(:,:) = 0.0
    ice_landcover(:,:)  = 0.0  !for Pb210 
    likely_coastal(:,:) = .false.
    mainly_sea(:,:)     = .false.

    do i = 1, limax
       do j = 1, ljmax

          dbgij = ( debug_proc .and. i == debug_li .and. j == debug_lj )
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
            ! Hard-coded for now
             LandCover(i,j)%ANTH(ilu) =  -1
             if( LandDefs(lu)%name == "WinterWheat" ) then
                LandCover(i,j)%SGS(ilu) =  nint( map2dGrowingSeasons(i,j,1) )
                LandCover(i,j)%ANTH(ilu) =  nint( map2dGrowingSeasons(i,j,2) )
                LandCover(i,j)%EGS(ilu) =  nint( map2dGrowingSeasons(i,j,3) )
             else if( LandDefs(lu)%name == "SpringWheat" ) then
                LandCover(i,j)%SGS(ilu) =  nint( map2dGrowingSeasons(i,j,4) )
                LandCover(i,j)%ANTH(ilu) =  nint( map2dGrowingSeasons(i,j,5) )
                LandCover(i,j)%EGS(ilu) =  nint( map2dGrowingSeasons(i,j,6) )
!TODO - find better solution!!
                LandCover(i,j)%SGS(ilu) =  LandCover(i,j)%EGS(ilu) - 90
             end if

             if ( DEBUG%LANDUSE>0 .and. dbgij ) &
                  write(*,"(a,i3,a20,3i4)")"LANDUSE: LU_SETGS", &
                  lu, LandDefs(lu)%name,&
                  LandCover(i,j)%SGS(ilu),LandCover(i,j)%ANTH(ilu), &
                  LandCover(i,j)%EGS(ilu)


             !/ for landuse classes with bulk-resistances, we only
             !  need to specify height once. Dummy values are assigned
             !  to LAI and gpot:

             if ( LandType(lu)%is_bulk ) then
                LandCover(i,j)%hveg(ilu) =  LandDefs(lu)%hveg_max
                LandCover(i,j)%LAI(ilu)  =  0.0          
                LandCover(i,j)%fphen(ilu) =  0.0          
             end if

             if ( LandType(lu)%is_water ) water_fraction(i,j) = &
                  water_fraction(i,j) + LandCover(i,j)%fraction(ilu)
             if ( LandType(lu)%is_ice   ) ice_landcover(i,j) = &
                  ice_landcover(i,j) + LandCover(i,j)%fraction(ilu)


          end do ! ilu

         ! Typically, we define as mainly sea when > 50% water, and
         ! likely_coastal when > 20%. SEA_LIMIT stores these numbers
         ! (We don't want to trust some squares with a mixture of sea
         !  and land for micromet purposes, e.g. T2 can be very wrong
         !  We mark these as likely coastal.)
         ! Unfortunately, we cannot yet determine if true sea or water

          if(water_fraction(i,j)>SEA_LIMIT(2) ) mainly_sea(i,j) = .true.
          if(water_fraction(i,j)>SEA_LIMIT(1).and. &
             water_fraction(i,j) < 0.999 ) likely_coastal(i,j) = .true.

          if ( DEBUG%LANDUSE>0 .and. dbgij )  then
             write(*,"(a,2i4,f7.3,2L2)") "SEACOAST ", i_fdom(i), j_fdom(j), &
                water_fraction(i,j), mainly_sea(i,j), likely_coastal(i,j)
          end if
       end do ! j
    end do ! i

    water_frac_set = .true.  ! just to inform other routines
    
  end subroutine InitLanduse
 
  subroutine ReadLanduse_CDF(filefound)
    ! Read data in other grid and interpolate to present grid
    !
    logical :: filefound
    integer :: i,j,lu, ilu, index_lu, maxlufound, iam, ifile, nFiles 
    real :: sumfrac
    integer, save :: ncalls=0

    character(len=40) :: varname, fShort
    character(len=200) :: fName, msg
    integer :: ncFileID, nDimensions,nVariables,nAttributes,timeDimID,varid
    integer :: nwords, err, xtype,ndims, status
    character(len=90) :: ewords(20)    ! LC:CF:EMEP, or /globa.../xxx/yyy
    logical :: fexist=.false.!file exist flag
  
    real, dimension(LIMAX,LJMAX,NLANDUSEMAX):: landuse_in ! tmp, with all data
    real, dimension(LIMAX,LJMAX,NLANDUSEMAX):: landuse_glob  ! CLM crude
    real, dimension(LIMAX,LJMAX):: landuse_tot ! CLM 
    real, dimension(LIMAX,LJMAX):: landuse_tmp ! tmp, with all data
    real    :: dbgsum, sum_veg

    real, dimension(LIMAX,LJMAX,NLUMAX):: landuse_data ! tmp, with all data
    integer, dimension(LIMAX,LJMAX):: landuse_ncodes ! tmp, with all data
    integer, dimension(LIMAX,LJMAX,NLUMAX):: landuse_codes ! tmp, with all data
    integer, dimension(size(FLUX_VEGS)):: iam_index = -1  !
    logical, dimension(NLANDUSEMAX)  :: is_veg 
    character(len=*), parameter :: dtxt='RdLanduseCDF:'
    logical, save :: mydbg=.false.  ! will set for debug_li, debug_lj 
    logical :: dbgij

    nFiles = size(LandCoverInputs%MapFile(:) )

    if( MasterProc )  then
       ncalls = ncalls + 1
       write(*,*) dtxt//" Starting", nFiles, ncalls
       do ifile = 1, nFiles
         write(*,*) 'MapFile ', trim(LandCoverInputs%MapFile(ifile))
       end do
    end if
    if( DEBUG%LANDUSE>0 .and. debug_proc )  mydbg = .true.

    maxlufound = 0   

    landuse_ncodes(:,:)   = 0     !***  initialise  ***
    landuse_codes(:,:,:)  = 0     !***  initialise  ***
    landuse_data  (:,:,:) = 0.0   !***  initialise  ***
    landuse_in  = 0.0              !***  initialise  ***
    landuse_tot   = 0.0              !***  initialise  *** Oct2017
    landuse_glob  = 0.0              !***  initialise  ***

    !Landusefile where landcodes are not predefined, but read from the file.
    ! Typilaclly, we will read a fine-scale map for the inner area (e.g.
    ! 'Landuse_PS_5km_LC.nc') and then a global map to make sure the
    ! full domain is covered, e.g. 'glc2000mCLM.nc'
 
    !1)check that file exists 
    !  (note that every processor opens/reads the same file)

    ilu=0 

    FILELOOP: do ifile = 1, nFiles ! size(LandCoverInputs%MapFile(:))
      fName = LandCoverInputs%MapFile(ifile)
     
      call wordsplit(fName,30,ewords,nwords,err,separator="/")
      fshort= '.../' //  ewords(nwords)    ! e.g. '.../Landuse_PS_5km.nc'
      status=nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncFileID)
      inquire(file=trim(fName),exist=fexist)

      if ( MasterProc .and. fexist) then
        write(*,'(a,i2,1x,a)') dtxt//"LANDUSE: found ", ifile, trim(fShort)
      end if

      call CheckStop (status /= nf90_noerr, &        ! AUG31
         dtxt//"LANDUSE: NOT found "//trim(fName) )  ! AUG31

     !get list of variables
      call check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,&
         timeDimID))
      ! All the inquire functions are inexpensive to use and require no I/O,
      ! since the information they provide is stored in memory when a netCDF
      ! dataset is first opened.   

      VARIDLOOP1: do varid=1,nVariables ! all variables in file

          if ( DEBUG%LANDUSE>0 )  CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)

          call check(nf90_Inquire_Variable(ncFileID,varid,varname,xtype,ndims))
          msg = dtxt//"checking "//trim(fShort)//': '// &
                trim(varname) !!! , index( varname, "LC:")

          ! landcover terms look like, e.g. LC:CF:EMEP
          if( index( varname, "LC:") < 1 ) then
            if ( mydbg ) write(*,*) trim(msg)//": Skipped, not LC"
            cycle ! ONLY LC: (LandCode) wanted
          end if
          call wordsplit(varname,3,ewords,nwords,err,separator=":")
          if ( mydbg ) write(*,*) trim(msg)//": Used from ", &
              trim(varname), ' ', ewords(3)
          if( ewords(3) /= "EMEP" .and.  ewords(3) /= "CLM" ) &
              cycle ! ONLY EMEP and CLM coded for now

          !=========================
          !! Old.  No need for IAM here now (will remove). See ADD_IAM below
          if( ewords(2)=="IAM_VEG" ) cycle
          !=========================

         !CHECK HERE to see if we already have this landcode:...
          lu = find_index( ewords(2), Land_codes(:) ) 
          if (  lu > 0 ) then
            if ( mydbg ) write(*,*) dtxt//"Already have"//ewords(2), lu
          else
            ilu = ilu + 1
            call CheckStop( ilu>NLANDUSEMAX , dtxt//&
            "NLANDUSEMAX smaller than number of landuses defined in file "//&
            trim(fname) )
            lu  = ilu  
            if ( mydbg ) write(*,*) dtxt//"Adding code"//ewords(2), ilu
            Land_codes(ilu) = ewords(2)    ! Landuse code found on file
          end if

          if(debug_proc) write(*,'(a,i3,1x,2f8.2,a)') dtxt//'IFILE', ifile,&
              glat(debug_li,debug_lj), glon(debug_li,debug_lj),trim(varname)

          call ReadField_CDF(trim(fName),varname,& 
               landuse_tmp,1,interpol='conservative', &
               needed=.true.,debug_flag=.false.,UnDef=-9.9E19) 

          if ( ifile == 1 ) then
             where ( landuse_tmp > 0.0 ) !Oct2017
               landuse_in(:,:,lu) = landuse_tmp
               landuse_tot(:,:) = landuse_tot(:,:) + landuse_tmp ! Sum of data from file1
             end where  !Oct2017
          else
               landuse_glob(:,:,lu) = landuse_tmp ! will merge below
          end if
 
          if ( debug_proc ) then
               print "(a,2i4,4es12.3,1x,a)", "F1 ", ifile,lu, &
                    landuse_tmp(debug_li,debug_lj), &
                    landuse_tot(debug_li,debug_lj), maxval(landuse_tmp(:,:)), &
                    landuse_glob(debug_li,debug_lj,lu), &
                    trim(ewords(2))
          end if
      end do VARIDLOOP1
      call check(nf90_close(ncFileID))!fname1
    end do FILELOOP ! DSLC

    NLand_codes=ilu !DSCLM now here


   ! MERGE inner and outer maps (Euro and Glob usually)
    !AUG31 if  ( EuroFileFound .and. GlobFileFound ) then ! we need to merge

    if ( nFiles > 1 ) then  ! we need to merge
      if ( debug_proc ) print '(a,i3,f12.4,5i6)', "F3  START", &
                  NLand_codes, landuse_tot(debug_li,debug_lj), me, &
                    limax, ljmax, debug_li, debug_lj
      do j = 1, ljmax
         do i = 1, limax
            !landuse_tmp can be numerically larger than 1.0 (1E-15 larger). 
            dbgij = ( mydbg .and. i==debug_li.and.j==debug_lj ) 

            if(landuse_tot(i,j)< 0.99999 ) then
              landuse_in(i,j,:)= 0.0  ! Will overwrite all PS stuff
              dbgsum = 0.0

              do ilu = 1, NLand_codes
                landuse_in(i,j,ilu) = min(1.0, landuse_glob(i,j,ilu) )
                dbgsum = dbgsum + landuse_in(i,j,ilu)
                if ( dbgij ) then
                   write(*, "(a,i3,3es15.6,1x,a)") "F4 ", ilu, &
                      landuse_in(debug_li,debug_lj,ilu), &
                      landuse_tot(debug_li,debug_lj), dbgsum,&
                      trim(Land_Codes(ilu))
                end if
              end do

            end if ! land_tot<0.9999
         end do  !j
      end do  !i
    end if  ! Euro, Glob

    if(MasterProc) then
      write(*,*)"LandDefs DONE ", ilu, Land_codes(ilu), &
          maxval( landuse_in ), minval(landuse_in)
      write( *,*) "CDFLAND_CODES: ", NLand_codes, " :"
      write( *,'((5a20))') Land_codes(1:NLand_codes)
    end if

    is_veg(:) = .false.
    do lu = 1, NLand_codes
      if (find_index( Land_codes(lu), FLUX_IGNORE(:) ) < 1 ) then
        if ( mydbg ) write(*,*) dtxt//'Some veg:'//trim(Land_codes(lu)), lu
        is_veg(lu) = .true.
      end if
    end do

  !!!!!!!!!!!! ADD IAM_VEG 
  ! Some "IAM" veg species can be defined for calculations of ozone fluxes.
  ! These are assigned very small land-area, using the mask which the IAM_VEG
  ! species gives.  We divide the area  by the nFluxVegs to keep the total area
  ! small.

  ! Append flux vegs to Land_codes
    do iam = 1, nFluxVegs
       NLand_codes = NLand_codes + 1
       call CheckStop( NLand_codes>NLANDUSEMAX , dtxt//&
            " flux vegs makes NLANDUSEMAX smaller than number of landuses defined")
       Land_codes(NLand_codes) = FLUX_vegs(iam) 
       iam_index(iam) = NLand_codes
       if ( mydbg ) write(*,*) dtxt//'IAM veg:'//trim(FLUX_VEGS(iam)), &
          iam, iam_index(iam)
    end do

    do i = 1, limax
       do j = 1, ljmax
          dbgij = ( mydbg .and. i==debug_li.and.j==debug_lj ) 
          sum_veg = 0.0
          do lu = 1, NLand_codes
             if ( is_veg(lu) .and. landuse_in(i,j,lu)  > 0.0 ) then
                sum_veg = sum_veg + landuse_in(i,j,lu) 
             end if
             if ( dbgij ) write(*,'(a,4i5,a13,L2,es12.3)') dtxt//'IAM add?', &
                  me,ncalls,i,j, trim(Land_codes(lu)), is_veg(lu), sum_veg
          end do
          if (  sum_veg < 1.0e-6 )  CYCLE

          if ( dbgij ) write(*,*) dtxt//'IAM nnn:', nFluxVegs
          IAM_VEG: do iam = 1, nFluxVegs !   size( FLUX_VEGS )
             !if ( len_trim( FLUX_VEGS(iam) ) < 1 ) then
             !   if(MasterProc) write(*,*)"Landuse SKIPS IAM ", iam
             !   cycle IAM_VEG
             !end if
           
             lu = iam_index(iam)
             landuse_in(i,j,lu) = 1.0e-3/real(nFluxVegs) ! Add
             if ( dbgij ) write(*,*) dtxt//'IAM add:', me, lu, trim(Land_codes(lu))
          end do IAM_VEG
       end do  !j
    end do  !i


   !!!!!!!!!!!! Now, convert to more compact arrays for export
    do i = 1, limax
       do j = 1, ljmax
          dbgij = ( mydbg .and. i==debug_li.and.j==debug_lj ) 
          do lu = 1, NLand_codes
             if ( dbgij ) then
                write(*,*) dtxt//'preGridAll', lu,NLand_codes, landuse_in(i,j,lu)
             end if
             if ( landuse_in(i,j,lu) > 0.0 ) then

                call GridAllocate("LANDUSE",i,j,lu,NLUMAX, &
                     index_lu, maxlufound, landuse_codes, landuse_ncodes)
                landuse_data(i,j,index_lu) = &
                     landuse_data(i,j,index_lu) + landuse_in(i,j,lu)!already in fraction unit
                if ( dbgij ) then
                  write(*,*) dtxt//'GridAll', lu, index_lu,& 
                     landuse_data(i,j,index_lu),  landuse_in(i,j,lu)
                end if
             end if
          end do ! lu
          LandCover(i,j)%ncodes  = landuse_ncodes(i,j)
          LandCover(i,j)%codes(:) = landuse_codes(i,j,:)
          LandCover(i,j)%fraction(:)  = landuse_data(i,j,:)
          sumfrac = sum( LandCover(i,j)%fraction(:) )

          if (  sumfrac < 0.99 .or. sumfrac > 1.01 ) then
               write(unit=errmsg,fmt="(a34,5i4,f12.4,6i4,2f7.2)") & !nb len(dtxt)=13
                 dtxt//" SumFrac Error ", me,i,j,  &
                    i_fdom(i),j_fdom(j), sumfrac, limax,  ljmax, &
                       i_fdom(1), j_fdom(1), i_fdom(limax), j_fdom(ljmax), &
                         glat(i,j), glon(i,j)
               write(*,*)'LandCover(i,j)%ncodes ',LandCover(i,j)%ncodes
               write(*,*)'LandCover(i,j)%codes(:) ',LandCover(i,j)%codes(1),LandCover(i,j)%codes(2:LandCover(i,j)%ncodes)
               write(*,*)'landuse_tot(i,j) ',landuse_tot(i,j)
               write(*,*)'landuse_glob(i,j,:) ',landuse_glob(i,j,:)
               print *, trim(errmsg)

               if(abs(sumfrac-1.0)<0.2.and.abs(glat(i,j))>89.0)then
                  write(*,*)'WARNING: ',trim(errmsg),sumfrac,glat(i,j)
               else
                   write(*,*)'lat/lon: ',trim(errmsg),glat(i,j), glon(i,j)
                 call CheckStop(errmsg)
               end if
          end if

       end do  !j
    end do  !i


    filefound=.true.
    if ( DEBUG%LANDUSE>0 ) then
        CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
        write(6,'(a,3i5)') dtxt//" me,  maxlufound, ncalls, cdf = ", &
                                  me, maxlufound, ncalls
    end if

  end subroutine ReadLanduse_CDF

  !=========================================================================
  subroutine  SetLandUse(daynumber, month)
    integer, intent(in) :: daynumber, month
    integer :: i,j,ilu,lu ! indices
    integer, save :: old_month = -1
    integer, save :: old_daynumber = -1
    logical, save :: my_first_call = .true.
    real :: hveg, lat_factor
    real :: xSAIadd
    integer :: pft
    logical, save :: debugProc = .false.
    logical :: dbgij, debug_sgs
    character (len=*), parameter :: dtxt='SetLandUse:'
    character (len=60) :: dnam !mainly for debug

! Treatment of growing seasons in the southern hemisphere:
!   all the static definitions (SGS,EGS...) refer to northern hemisphere, 
!   but the actual simulation dates are shifted by 6 months in the southern
!   hemisphere by using uses effectivdaynumber and 
!   mod(current_date%month+5,12)+1 in southern hemis



   !======================================================================
    if ( my_first_call ) then

        if ( DEBUG%LANDUSE>0 .and. debug_proc ) debugProc = .true.

       !read in data from file
        my_first_call   = .false.
        
        call InitLanduse(daynumber)

       ! The DO3SE params are needed for the call to fPhenology
      
        call Init_DO3SE(IO_DO3SE, LandCoverInputs%Do3seDefs, &
             NLand_codes, Land_codes, errmsg)
        call CheckStop(errmsg, "Reading DO3SE ")

    end if ! my_first_call
   !======================================================================

    if ( daynumber == old_daynumber ) then
        my_first_call = .false.
        return
    end if
    old_daynumber = daynumber

    if(MasterProc) write(*,*) dtxt//" day, pfts? ", daynumber, USES%PFT_MAPS
    if(debugProc ) write(*,"(a,5i5,L2)") dtxt//" debug me i j pft? ", me, &
         debug_li, debug_lj, limax, ljmax, USES%PFT_MAPS


   !Landcover data can be set either from simplified LPJ
   !PFTs, or from the "older" DO3SE inputs file

    if ( USES%PFT_MAPS ) then !- Check for LPJ-derived data -
         if (MasterProc) write(*,*) dtxt//"New PFTMAPS ", month, old_month
         if ( month /= old_month ) then 
           call MapPFT_LAI( month )
           old_month = month
         end if
    end if
    
!PALEO LANDUSE
!    if( PALEO_TEST ) then
!     call SetPaleo(daynumber, month)
!    end if

    do i = 1, limax
      do j = 1, ljmax

        debug_sgs = ( DEBUG%LANDUSE > 1 .and. &
        current_date%hour == 0  .and.  &
        glat(i,j) <   1.0 .and. glat(i,j) > -1.0 .and.  &
        glon(i,j) > -72.0 .and. glon(i,j) < -70.0 )


        effectivdaynumber=daynumber
       ! effectiv daynumber to shift 6 months when in southern hemisphere
        if(glat(i,j)<0.0)effectivdaynumber=mod(daynumber+182,nydays)+1 

        dbgij = ( debugProc .and. i == debug_li .and. j == debug_lj ) 

        if ( dbgij ) then
           write(*,"(a,i3,9i6)") dtxt//" debug DATE ", &
              LandCover(i,j)%ncodes, daynumber, current_date
        end if

        do ilu= 1, LandCover(i,j)%ncodes
          lu      = LandCover(i,j)%codes(ilu)
          pft     = LandType(lu)%pft
          dnam    = dtxt//trim(LandDefs(lu)%name)

          if ( dbgij ) write(*, *) trim(dnam)//" lu pft", lu, pft,&
            LandType(lu)%is_bulk, LandType(lu)%is_forest

          if ( LandType(lu)%is_bulk ) then
            LandCover(i,j)%LAI(ilu) = 0.0
            LandCover(i,j)%SAI(ilu) = 0.0
            cycle    
          end if!else Growing veg present:

          if ( LandDefs(lu)%name == "MED_OAK" .or.  &
                  LandDefs(lu)%name == "MED_PINE"   ) then

            LandCover(i,j)%LAI(ilu) = MedLAI(effectivdaynumber, &
               100, 166, & ! Hard-code from Mapping Manual
                  LandDefs(lu)%LAImin, LandDefs(lu)%LAImax )
            if ( dbgij ) then
               write(*,"(a,3i4,3f8.3)") "MED_TREE "//&
                  trim(LandDefs(lu)%name), effectivdaynumber,&
                  LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu),  &
                  LandDefs(lu)%LAImin, LandDefs(lu)%LAImax, &
                  LandCover(i,j)%LAI(ilu)
            end if

          else
            LandCover(i,j)%LAI(ilu) = Polygon(effectivdaynumber, &
                0.0, LandDefs(lu)%LAImin, LandDefs(lu)%LAImax,&
                  LandCover(i,j)%SGS(ilu), LandDefs(lu)%SLAIlen, &
                     LandCover(i,j)%EGS(ilu), LandDefs(lu)%ELAIlen)
          end if

          LandCover(i,j)%fphen(ilu) = fPhenology( lu ,effectivdaynumber &
                ,LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu) ,dbgij )

          if ( dbgij ) then
             write(*,"(a,3i4,5f8.3)")trim(dnam)//" CHECK_VEG ",&
                effectivdaynumber, &
                LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu),  &
                LandDefs(lu)%LAImin, LandDefs(lu)%LAImax,&
                LandCover(i,j)%LAI(ilu), LandCover(i,j)%fphen(ilu)
             write(*,"(a,L3)")trim(dnam)//' CHECK_PFT', USES%PFT_MAPS
          end if

          if ( USES%PFT_MAPS ) then
            if ( DEBUG%PFT_MAPS > 0 .and. dbgij ) then
              if ( pft > 0 ) then
                write(*,"(a,i4,i6,2f8.3)") trim(dnam)//" PFTS COMP? ", &
                          daynumber, pft, LandCover(i,j)%LAI(ilu), &
                           pft_lai(i,j, pft)*LandDefs(lu)%LAImax
              else
                write(*,"(2a,i4,i6,2f8.3)") trim(dnam)//" PFTS COMP? ", &
                    daynumber, pft, LandCover(i,j)%LAI(ilu), -1.0
              end if
            end if

            if ( pft > 0 ) then !PFT OVERWRITE!
              LandCover(i,j)%LAI(ilu)= pft_lai(i,j, pft)*LandDefs(lu)%LAImax
              if(dbgij) write(*,"(a,2i4,5f8.3)")dtxt//' CHECK_LAI', lu, pft, &
                      pft_lai(i,j, pft),LandDefs(lu)%LAImax
                    LandCover(i,j)%fphen(ilu)= 1.0  ! Skip fphen if using PFT
                    LandCover(i,j)%SGS(ilu)=  -999  ! Marker, since not used
                    LandCover(i,j)%EGS(ilu)=  -999  ! Marker, since not used
              end if
              if (debug_sgs  )write(*,*) "ESGS in PFTMAPS"
            end if


            hveg = LandDefs(lu)%hveg_max   ! defaults
            xSAIadd = 0.0

            if (  LandType(lu)%is_crop ) then

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
            else if( LandType(lu)%is_seminat ) then !A2017 SNL
              if ( glat(i,j) >= 60.0 ) then
                lat_factor  = max(0.3, ( 1.0 - 0.05* (glat(i,j)-60.0)) )
                hveg  = hveg *  lat_factor
                LandCover(i,j)%LAI(ilu) = LandCover(i,j)%LAI(ilu)  * lat_factor
              end if
              LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)  + 0.5 !  A2017SNL
            else
              LandCover(i,j)%SAI(ilu) = LandCover(i,j)%LAI(ilu)   !defaults
            end if

            LandCover(i,j)%hveg(ilu) =  hveg


            if( debug_sgs .or. dbgij  ) then
              write(*, "(a20,i4,2f7.1,3i5,f8.2)") trim(dnam)//":ESGS:",&
                lu, glat(i,j), glon(i,j), daynumber, effectivdaynumber, &
                Landcover(i,j)%SGS(ilu), Landcover(i,j)%LAI(ilu)
              write(*,"(a,3i4,5f8.3)")trim(dnam)//"CHECK_VEGB:",&
                 effectivdaynumber, &
                 LandCover(i,j)%SGS(ilu), LandCover(i,j)%EGS(ilu),  &
                 LandDefs(lu)%LAImin, LandDefs(lu)%LAImax,&
                 LandCover(i,j)%LAI(ilu), LandCover(i,j)%fphen(ilu)
            end if ! debug_sgs
         end do ! lu
       end do ! j
    end do ! i

    my_first_call   = .false.

! --- print out for debug cell
    if ( DEBUG%LANDUSE>0.and.debug_proc ) then
       i=debug_li
       j=debug_lj

       do ilu= 1, LandCover(i,j)%ncodes
          lu      = LandCover(i,j)%codes(ilu)
          pft     = LandType(lu)%pft
          if ( LandType(lu)%is_bulk ) cycle    !else Growing veg present:

            write(*,"(a,i3,a16,i4,f7.2,3f8.3,2i5)") "LANDUSE Phen ", lu,&
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
function MedLAI(jday,LAIs,LAIe,LAImin,LAImax) result (LAI)
!=======================================================================

!     Calculates the value of LAI from the Mapping manual
!     functions for Mediteranean forests

!   Inputs
    integer, intent(in) :: jday     !day of year, after any co-ordinate change
    integer, intent(in) ::  LAIs,LAIe
    real, intent(in) ::     LAImin,LAImax

!  Output:
    real ::   LAI  ! value at day jday


    if( jday <= LAIs ) then
      LAI = 0.35*((real(LAIs-jday))/LAIs) + LAImin
    else if ( jday < (366.0-LAIe) ) then
      LAI = (LAImax-LAImin)*(real(jday-LAIs)/LAIs) + LAImin
    else
      LAI = (LAImax-( LAImin+0.35))*(real(366-jday)/LAIe) + LAImin+0.35
    end if


 end function MedLAI

 !=======================================================================

end module Landuse_mod
