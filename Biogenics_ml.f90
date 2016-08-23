! <Biogenics_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
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
!*****************************************************************************! 
module Biogenics_ml

  !/-- Reads in BVOC emisions factors 
  !
  !     1) From defaults globally
  !
  !     2) from local file if available (e.g. Europe, used by default)
  !
  !   Terminology:
  !
  !    LCC = land cover class, e.g. DF = decid forest
  !
  !    EF = Emission factor at 30 deg C, full sunlight (1000 uE)
  !         ug/g/hr
  !
  !    Em = Emissions, = EF * LAI * fraction of grid
  !         ug/m2(ground)/hr
  !
  !    The code will assign EFs from the local data if available, otherwise
  !    use defaults which have been readfrom Inputs_Landuse, Eiso, Emtp, Emtl. 
  !    Note that we use emissions per m2 of vegetation, not per 
  !    m2 of grid. This lets the model use the landcover from any veg-map
  !    without having to make this consistent with the EF maps - the latter
  !    are regarded as smoothly varying fields which can be interpolated 
  !    by the ReadField_CDF interpolation routines. No need to worry about
  !    conserving these very imperfect numbers accurately ;-)
  !
  !    Dave Simpson, 2010-2012
  !---------------------------------------------------------------------------

  use CheckStop_ml,      only: CheckStop
  use ChemChemicals_ml,  only : species
  use GridValues_ml    , only : i_fdom,j_fdom, debug_proc,debug_li,debug_lj
  use Io_ml            , only : IO_FORES, open_file, ios, PrintLog, datewrite
  use KeyValue_ml,       only : KeyVal,KeyValue
  use LandDefs_ml,       only: LandType, LandDefs
  use LandPFT_ml,        only: MapPFT_LAI, pft_lai
  use Landuse_ml,        only : LandCover
  use LocalVariables_ml, only : Grid  ! -> izen, DeltaZ
  use MetFields_ml,      only : t2_nwp
  use ModelConstants_ml, only : USE_EURO_SOILNOX, USE_GLOBAL_SOILNOx, DEBUG_SOILNOX, USE_SOILNH3
  use ModelConstants_ml, only : NPROC, MasterProc, TINY, &
                           USE_PFT_MAPS, NLANDUSEMAX, IOU_INST, & 
                           KT => KCHEMTOP, KG => KMAX_MID, & 
                           EURO_SOILNOX_DEPSCALE, & 
                           DEBUG_BIO, BVOC_USED, MasterProc
  use NetCDF_ml,        only : ReadField_CDF, printCDF
  use OwnDataTypes_ml,  only : Deriv, TXTLEN_SHORT
  use Par_ml   , only :  MAXLIMAX,MAXLJMAX,MSG_READ1,me, limax, ljmax
  use Par_ml,            only : limax, ljmax, MAXLIMAX, MAXLJMAX, me
  use PhysicalConstants_ml,  only :  AVOG, GRAV
  use Radiation_ml,          only : PARfrac, Wm2_uE
  use Setup_1dfields_ml,     only : rcemis  
  use SmallUtils_ml, only : find_index
  use TimeDate_ml,       only : current_date, daynumber
  implicit none
  private

  !/-- subroutines for BVOC
  public ::  Init_BVOC
  private :: Get_LCinfo
  public ::  GetEuroBVOC
  private :: MergedBVOC
  public ::  setup_bio
  public ::  SetDailyBVOC
  private :: TabulateECF

  !/-- subroutines for soil NO
  public :: Set_SoilNOx

  INCLUDE 'mpif.h'
  include 'CM_EmisBioNat.inc'
  !e.g.
  !  integer, parameter, public ::  NEMIS_BioNat  = 3
  !  character(len=7), save, dimension(NEMIS_BioNat), public:: &
  !    EMIS_BioNat =  (/ "C5H8   " , "APINENE" , "NO     " /)
 
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
  integer, public, parameter :: N_ECF=2, ECF_ISOP=1, ECF_TERP=2
  integer, public, parameter :: BIO_ISOP=1, BIO_MTP=2, &
                                 BIO_MTL=3 ! , BIO_SOILNO=4, BIO_SOILNH3=5
  integer, public, parameter :: BIO_TERP=2 ! Used for final emis, sum of MTP+MTL
  integer, public, save ::  last_bvoc_LC   !max index land-cover with BVOC (min 4)
                                                        
  ! Soil NOx
   real,public, save, allocatable, dimension(:,:) :: &
      AnnualNdep, &  ! N-dep in mgN/m2/
      SoilNOx, SoilNH3


 ! Set true if LCC read from e.g. EMEP_EuroBVOC.nc:
 ! (Currently for 1st four LCC, CF, DF, BF, NF)
  logical, private, dimension(NLANDUSEMAX), save :: HaveLocalEF 

! real, public, save, dimension(MAXLIMAX,MAXLJMAX,size(BVOC_USED)+NSOIL_EMIS) :: &

  ! EmisNat is used for BVOC; soil-NO, also in futur for sea-salt etc.
  ! Main criteria is not provided in gridded data-bases, often land-use
  ! dependent.

  real, public, save, allocatable, dimension(:,:,:) :: &
     EmisNat       !  will be transferred to d_2d emis sums


  !standard emission factors (EFs) per LC
  !Need to dimension later for Emtp, Emtl, last_bvoc_LC
  real, private, save, allocatable, dimension(:,:,:,:) :: &
     bvocEF       !  Gridded std. emissions per PFT

  !standard emission factors per LC for daily LAI
  real, private, save, allocatable, dimension(:,:,:) :: &
     day_embvoc    !  emissions scaled by daily LAI

  logical, private, save, allocatable, dimension(:,:) :: EuroMask

  !/-- Canopy environmental correction factors-----------------------------
  !
  !    - to correct for temperature and light of the canopy 
  !    - from Guenther's papers. (Limit to 0<T<40 deg C.)

  real, public, save, dimension(N_ECF,40) :: canopy_ecf  ! Canopy env. factors

 ! Indices for the species defined in this routine. Only set if found
  integer, private, save :: ispec_C5H8, ispec_APIN, ispec_NO , ispec_NH3
  integer, private, save :: itot_C5H8,  itot_APIN,  itot_NO , itot_NH3

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    subroutine Init_BVOC()

!    Read natural BVOC emission potentials
!-----------------------------------------------------------------------------
!   Emission potentials (EFs) now read a netcdf of ca. 50x50 resolution.
!   This file, EMEP_EuroBVOC.nc uses the ICP-Forests species map as processed
!   by Renate Koeble at JRC (e.g.
!   EFs now a mixure of rates from Simpson et al., 1999, JGR, Vol 104, D7, 
!   8113-8152, and Keenan, T. et al., ACP, 2009, 9, 4053-4076 
!   See Simpson et al., ACP, 2012

    integer :: alloc_err
    
    allocate(AnnualNdep(MAXLIMAX,MAXLJMAX),SoilNOx(MAXLIMAX,MAXLJMAX), SoilNH3(MAXLIMAX,MAXLJMAX))
    allocate(EmisNat(NEMIS_BioNat,MAXLIMAX,MAXLJMAX))
    EmisNat=0.0
    allocate(day_embvoc(MAXLIMAX,MAXLJMAX,size(BVOC_USED)))
    day_embvoc = 0.0
    allocate(EuroMask(MAXLIMAX,MAXLJMAX))
    EuroMask=.false.

      if ( size(BVOC_USED) == 0 ) then
        call PrintLog("No Biogenic Emissions ", MasterProc)
        return
      end if

   !====================================
   ! get indices.  NH3 not yet used.
      ispec_C5H8 = find_index( "C5H8", EMIS_BioNat(:) ) 
      ispec_APIN = find_index( "APINENE", EMIS_BioNat(:) ) 
      ispec_NO   = find_index( "NO", EMIS_BioNat(:) ) 
      ispec_NH3  = find_index( "NH3", EMIS_BioNat(:) ) 
      call CheckStop( ispec_C5H8 < 1 , "BiogencERROR C5H8")
      !call CheckStop( ispec_APIN < 1 , "BiogencERROR APIN")
      if( ispec_APIN < 0 ) call PrintLog("WARNING: No APINENE Emissions")
     
      call CheckStop( USE_EURO_SOILNOX .and. ispec_NO < 1 , "BiogencERROR NO")
      call CheckStop( USE_GLOBAL_SOILNOX .and. ispec_NO < 1 , "BiogencERROR NO")
      if( MasterProc ) write(*,*) "SOILNOX ispec ", ispec_NO

      itot_C5H8 = find_index( "C5H8", species(:)%name    ) 
      itot_APIN = find_index( "APINENE", species(:)%name )
      itot_NO   = find_index( "NO", species(:)%name      )
      itot_NH3  = find_index( "NH3", species(:)%name      )

   !====================================
 
    call TabulateECF()   ! Tabulates temp functions
   !====================================

    call Get_LCinfo() ! Gets landcover info, last_bvoc_LC

    allocate(  bvocEF(MAXLIMAX,MAXLJMAX,last_bvoc_LC,size(BVOC_USED)),&
        stat=alloc_err )
    call CheckStop( alloc_err , "bvocEF alloc failed"  )

    bvocEF(:,:,:,:) = 0.0

   !========= Read in Standard (30 deg C, full sunlight emissions factors = !
   ! Remember factor used in Emissions_ml to get from ug/m2/s
   ! to molecules/cm2/s  (needs more documentation!)

   !====================================

    call GetEuroBVOC()
   !====================================

   !====================================
   ! Merges Local and global/defaults, and scales with land-cover area
   ! Emissions factors shoudl now by ug/m2(grid)/h

    call MergedBVOC() 
   !====================================

   !========================================================================!

     ! old summation. Kept to demonstrate mpi_allreduce
     !output sums. Remember that "shadow" area not included here.
     ! do i = 1,  2 !! size(BVOC_USED) 
     !    bvocsum   = sum ( emforest(li0:li1,lj0:lj1,i) )
     !    CALL MPI_ALLREDUCE(bvocsum,bvocsum1, 1, &
     !      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
     !    if ( MasterProc  ) write(6,"(a20,i4,2es12.4)") &
     !         'Biogenics_ml, ibio, sum1',i, bvocsum, bvocsum1
     ! end do

   end subroutine Init_BVOC
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine Get_LCinfo() ! Gets landcover info, last_bvoc_LC
      ! Checks for default bvoc emissions from each landcover category
      ! (read from Inputs_LandDefs.csv file)
      integer :: iL

        do iL= 1, size(LandType(:)%pft )
            
            if( LandDefs(iL)%Eiso > 0 ) last_bvoc_LC = iL
            if( LandDefs(iL)%Emtp > 0 ) last_bvoc_LC = iL
            if( LandDefs(iL)%Emtl > 0 ) last_bvoc_LC = iL

         end do
         if( MasterProc.and.DEBUG_BIO)&
              write(*,*) " LAST BVOC LC from LandDefs:",last_bvoc_LC

       ! We need at least 4 for CF, DF, NF, BF in Euro file
         last_bvoc_LC =  max(last_bvoc_LC, 4 ) 

    end subroutine Get_LCinfo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   subroutine  GetEuroBVOC()

!.....................................................................
!**    DESCRIPTION:

!    Reads the processed BVOC emission potentials.

!    real    :: loc(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
    integer :: iVeg, iEmis, ibvoc, i,j
    character(len=1000) :: varname
    character(len=2), dimension(4) :: VegName = (/ "CF", "DF", "NF", "BF" /)


       !varname = "Fake"
       !HaveLocalEF(:) = .false.
       !call ReadField_CDF('EMEP_EuroBVOC.nc',varname,&
       !    loc,1,interpol='zero_order',needed=.true.,debug_flag=.true.)
       !if( debug_proc ) write(*,*)  "EMEP_EuroBVOC i,j fake ", &
       !   loc(debug_li, debug_lj)
     do iVeg = 1, size(VegName)
       ibvoc = find_index( VegName(iveg), LandDefs(:)%code )
       HaveLocalEF(ibvoc) = .true.
       do iEmis = 1, size(BVOC_USED)
          varname = trim(BVOC_USED(iEmis)) // "_" // trim(VegName(iVeg))
          
          call ReadField_CDF('EMEP_EuroBVOC.nc',varname,&
             bvocEF(:,:,ibvoc,iEmis),1,interpol='zero_order',needed=.true.,debug_flag=.false.,UnDef=-999.0)

 
         if( debug_proc ) then
              write(*, "(2a,f12.3,3i2)") "EURO-BVOC:E ", &
             trim(varname), bvocEF(debug_li, debug_lj,ibvoc,iEmis), iVeg, ibvoc, iEmis
              write(*, "(2a,2es12.3)") "EURO-BVOC:minmax ", &
             trim(varname), minval(bvocEF(:,:,ibvoc,iEmis)), maxval(bvocEF(:,:,ibvoc,iEmis))
         end if     

       end do

       
      ! Make a mask where we can use the local bvoc. Should be the same from
      ! all EFs, since only non-def areas set to -999, otherwise zero or +
      ! If any values exist, should exist for all entries, hence check.
       iEmis=size(BVOC_USED)
       if( iVeg == 1 )  then
          where(bvocEF(:,:,ibvoc,iEmis)>-1.0)
            EuroMask = .true.
          end where
       else  ! Just check that following maps are consistent
           do i=1,limax
           do j=1,ljmax
             if ( EuroMask(i,j) .and. bvocEF(i,j,ibvoc,iEmis)<0.0 ) then
               write(*,*) "MASK ERROR", me, i_fdom(i), j_fdom(j)
               call CheckStop("EuroMask BVOC ERROR")
             end if
           end do
           end do
       end if
                
     end do

  end subroutine GetEuroBVOC
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine MergedBVOC()

      integer :: i, j, nlu, iL, iiL
      integer :: pft
      real :: biso, bmt    !  Just for printout
      logical :: use_local, debug_flag

      if( MasterProc .and.DEBUG_BIO) write(*,*) "Into MergedBVOC"

      do i = 1, limax !PPP MAXLIMAX
      do j = 1, ljmax !PPP MAXLJMAX

        nlu = LandCover(i,j)%ncodes

        use_local = EuroMask(i,j) 
        debug_flag = ( debug_proc .and. debug_li == i .and. debug_lj == j )

        LULOOP: do iiL= 1, nlu

            iL      = LandCover(i,j)%codes(iiL)
            pft     = LandType(iL)%pft

           if( debug_flag ) then
               write(*,"(a,2i7,2L2,i3)") &
                   "TryMergeBVOC" //trim(LandDefs(iL)%name), iL, pft, &
                     use_local, HaveLocalEF(iL), last_bvoc_LC
           end if

           if( use_local .and. HaveLocalEF(iL) ) then 

                ! Keep EFs from EuroBVOC
                if( debug_flag ) write(*,*) "MergeBVOC: Inside local"

           else if ( iL <= last_bvoc_LC ) then ! otherwise use defaults
               bvocEF(i,j,iL,BIO_ISOP) = LandDefs(iL)%Eiso * LandDefs(iL)%BiomassD 
               bvocEF(i,j,iL,BIO_MTP)  = LandDefs(iL)%Emtp * LandDefs(iL)%BiomassD
               bvocEF(i,j,iL,BIO_MTL)  = LandDefs(iL)%Emtl * LandDefs(iL)%BiomassD
                if( debug_flag ) write(*,"(a,i3,8f8.2)") &
                  "MergeBVOC: Outside local", iL, LandDefs(iL)%BiomassD,&
                   LandDefs(iL)%Eiso, LandDefs(iL)%Emtp, LandDefs(iL)%Emtl
           else
                if( debug_flag ) write(*,*) "MergeBVOC: Outside LCC", iL
           end if


           if( debug_flag ) then

              biso = 0.0
              bmt  = 0.0
              if ( iL <= last_bvoc_LC ) then
                biso   = bvocEF(i, j,iL, BIO_ISOP) 
                bmt    = bvocEF(i, j,iL, BIO_TERP) 
              end if
              write(*,"(a,2i4,2L2,f9.4,9f10.3)") "MergeBVOC", &
                 iL, pft,  use_local, HaveLocalEF(iL),  &
                   LandCover(i,j)%fraction(iiL), biso, bmt, LandDefs(iL)%Eiso, &
                     LandDefs(iL)%Emtp, LandDefs(iL)%Emtl
           end if 
        end do LULOOP
      end do
      end do

   end subroutine MergedBVOC
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine SetDailyBVOC(daynumber)

      ! Scales emission potentials for daily LAI changes

      integer, intent(in) :: daynumber
      integer, save :: last_daynumber = -999, alloc_err
      integer :: i, j, nlu, iL, iiL, ibvoc
      real :: LAIfac  ! Multiplies by land-fraction
      real :: b       !  Just for printout
      logical :: debug
      logical, save :: my_first_call = .true.
      real, allocatable, dimension(:,:) ::  workarray

      if( MasterProc .and.DEBUG_BIO ) write(*,"(a,3i5)") "Into SetDailyBVOC", &
            daynumber, last_daynumber, last_bvoc_LC

      if ( daynumber == last_daynumber ) return
      last_daynumber = daynumber

      if ( DEBUG_BIO .and.  my_first_call  ) then
           allocate(  workarray(MAXLIMAX,MAXLJMAX), stat=alloc_err )
           call CheckStop( alloc_err , "workarray alloc failed"  )
           workarray = 0.0
      end if

      do i = 1, limax !PPP MAXLIMAX
      do j = 1, ljmax !PPP MAXLJMAX

        nlu = LandCover(i,j)%ncodes

        day_embvoc(i,j,:) = 0.0
        debug = ( DEBUG_BIO .and. debug_proc .and.  &
                   debug_li == i .and. debug_lj == j )

        LULOOP: do iiL= 1, nlu

            iL      = LandCover(i,j)%codes(iiL)

            if ( iL >  last_bvoc_LC ) cycle

            !for tundra and wetlands we have zero LAI, so omit
            !LAI scaling. Not an ideal system.... rewrite one day.
             
              if( LandCover(i,j)%LAI(iiL)< 1.0e-5 ) then ! likely wetlands, tundra
                 LAIfac = 1.0
                 if( debug ) write(*,*)"BVOC TUNDRA/WETLANDS",iL,LandCover(i,j)%LAI(iiL)
              else
                LAIfac = LandCover(i,j)%LAI(iiL)/LandDefs(IL)%LAImax
                LAIfac= min(LAIfac, 1.0)
              end if
              LAIfac = LAIfac * LandCover(i,j)%fraction(iiL)
              

              do ibvoc = 1, size(BVOC_USED) 
                day_embvoc(i,j,ibvoc) = day_embvoc(i,j,ibvoc) + &
                   LAIfac * max(1.0e-10,bvocEF(i,j,iL,ibvoc))
                   !done above LandCover(i,j)%fraction(iiL) *  &
              end do

              if ( debug ) then
                 b = 0.0
                 if ( iL <= last_bvoc_LC ) b = bvocEF(i, j,iL, BIO_ISOP)
                 write(*,"(a,a10,2i5,f9.5,2f7.3,8f10.3)") "SetBVOC", &
                  trim(LandDefs(iL)%name), daynumber, iL, &
                   LandCover(i,j)%fraction(iiL), &
                   LandCover(i,j)%LAI(iiL), LandDefs(iL)%LAImax, b,&
                     ( day_embvoc(i, j, ibvoc), ibvoc = 1, size(BVOC_USED) ) 
                  
              end if
              ! When debugging it helps with an LAI map
              if( DEBUG_BIO .and. my_first_call ) &
                 workarray(i,j) = workarray(i,j) + &
                    LandCover(i,j)%LAI(iiL)*LandCover(i,j)%fraction(iiL)
        end do LULOOP
      end do
      end do

      if ( DEBUG_BIO ) then
!         if ( my_first_call  ) then ! print out 1st day
!              call printCDF("BIO-LAI", workarray, "m2/m2" )
!              workarray(:,:) = day_embvoc(:,:,1)
!              call printCDF("BIO-Eiso", workarray, "ug/m2/h" )
!              workarray(:,:) = day_embvoc(:,:,2) + day_embvoc(:,:,3)
!              call printCDF("BIO-Emt", workarray, "ug/m2/h" )
!              deallocate(  workarray )
!         end if
       end if
       my_first_call = .false.

   end subroutine SetDailyBVOC
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    subroutine TabulateECF()
    !/-- Tabulate canopy environmental correction factors
    !-----------------------------------------------------
    ! ag suffix  = Alex Guenther's parameter values

    real    :: agts, agr, agtm, agct1, agct2, agct ,itk
    integer :: it, i

    agts = 303.
    agr = 8.314
    agtm = 314.
    agct1 = 95000.
    agct2 = 230000.

    do it = 1,40
      itk = it + 273.15
      agct = exp(agct1*(itk - agts)/(agr*agts*itk)) / &
                (1. + exp(agct2*(itk - agtm)/(agr*agts*itk)))

      canopy_ecf(ECF_ISOP,it) = agct

      ! Terpenes
      agct = exp( 0.09*(itk-agts) )
      canopy_ecf(ECF_TERP,it) = agct

      !?? for terpene fac = 0.5*fac(iso): as mass terpene = 2xmass isoprene

      if(DEBUG_BIO  .and.  MasterProc ) &
             write(6,"(A12,i4,5g12.3)") 'Biogenic ecfs: ', &
                  it, ( canopy_ecf(i,it), i=1, N_ECF )
    end do
    end subroutine TabulateECF
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine setup_bio(i,j)
  !
  !---- assign isoprene rates  ------------------------------------------------
  !
  !  So far, assigns isoprene using surface (2m) temperature, and for all
  !  zenith angles <90. Should include light dependance at some stage
  !
  !  Output : rcbio added to rcemis - isoprene emissions for 1d column
  !
  !  Called from setup_ml, every  advection step.
  !----------------------------------------------------------------------------

  integer, intent(in) ::  i,j

  integer :: it2m
  real    :: E_ISOP, E_MTP, E_MTL

! To get from ug/m2/h to molec/cm3/s
! ug -> g  1.0e-9; g -> mole / MW; x AVOG
! will need /Grid%DeltaZ
  real, parameter :: & 
        biofac_ISOP   = 1.0e-12*AVOG/68.0 /3600.0  &
       ,biofac_TERP   = 1.0e-12*AVOG/136.0/3600.0  &
       ,biofac_SOILNO = 1.0e-12*AVOG/14.0 /3600.0  &
       ,biofac_SOILNH3= 1.0e-12*AVOG/14.0 /3600.0  

 ! Light effects added for isoprene emissions

  real            :: par   ! Photosynthetically active radiation
  real            :: cL    ! Factor for light effects
  real, parameter :: &
      CL1 = 1.066  , &    ! Guenther et al's params
      ALPHA = 0.0027      ! Guenther et al's params

  if ( size(BVOC_USED) == 0  ) return   ! e.g. for ACID only


  it2m = nint( Grid%t2C - TINY )
  it2m = max(it2m,1)
  it2m = min(it2m,40)

  !ASSUME C5H8 FOR NOW if ( ispec_C5H8 > 0 ) then
    if ( Grid%izen <= 90) then ! Isoprene in daytime only:

     ! Light effects from Guenther G93

      par = (Grid%Idirect + Grid%Idiffuse) * PARfrac * Wm2_uE

      cL = ALPHA * CL1 * par/ sqrt( 1 + ALPHA*ALPHA * par*par)

     ! E in ug/m2/h

      E_ISOP = day_embvoc(i,j,BIO_ISOP)*canopy_ecf(BIO_ISOP,it2m) * cL 

      ! Add light-dependent terpenes to pool-only
      if(BIO_TERP > 0) E_MTL = &
             day_embvoc(i,j,BIO_MTL)*canopy_ecf(ECF_TERP,it2m)*cL

     !  molecules/cm3/s
     ! And we scale EmisNat to get units kg/m2 consistent with
     ! Emissions_ml (snapemis).  ug/m2/h -> kg/m2/s needs 1.0-9/3600.0. 

      E_ISOP = day_embvoc(i,j,BIO_ISOP)*canopy_ecf(BIO_ISOP,it2m) * cL 

      rcemis(itot_C5H8,KG)   = rcemis(itot_C5H8,KG) + E_ISOP * biofac_ISOP/Grid%DeltaZ
      EmisNat(ispec_C5H8,i,j)= E_ISOP * 1.0e-9/3600.0

  else ! night
     EmisNat(ispec_C5H8,i,j) = 0.0
     E_MTL = 0.0
     E_ISOP = 0.0
  endif ! daytime

    if ( ispec_APIN > 0 ) then

     ! add pool-only terpenes rate;
        E_MTP = day_embvoc(i,j,BIO_MTP)*canopy_ecf(ECF_TERP,it2m)
        rcemis(itot_APIN,KG)    = rcemis(itot_APIN,KG) + &
               (E_MTL+E_MTP) * biofac_TERP/Grid%DeltaZ
        EmisNat(ispec_APIN,i,j) = (E_MTL+E_MTP) * 1.0e-9/3600.0
    end if

    if ( USE_EURO_SOILNOX ) then
        rcemis(itot_NO,KG)    = rcemis(itot_NO,KG) + &
             SoilNOx(i,j) * biofac_SOILNO/Grid%DeltaZ
        EmisNat(ispec_NO,i,j) =  SoilNOx(i,j) * 1.0e-9/3600.0
    else if ( USE_GLOBAL_SOILNOX ) then !TEST
        EmisNat(ispec_NO,i,j) =  SoilNOx(i,j)*Grid%DeltaZ/biofac_SOILNO * 1.0e-9/3600.0
    end if

    !EXPERIMENTAL
    if ( USE_SOILNH3 ) then
        rcemis(itot_NH3,KG)    = rcemis(itot_NH3,KG) + &
            SoilNH3(i,j) * biofac_SOILNH3/Grid%DeltaZ
        EmisNat(ispec_NH3,i,j) =  SoilNH3(i,j) * 1.0e-9/3600.0
    end if
     
 
    if ( DEBUG_BIO .and. debug_proc .and. i==debug_li .and. j==debug_lj .and. &
         current_date%seconds == 0 ) then

      call datewrite("DBIO env ", it2m, (/ max(par,0.0), max(cL,0.0), &
            canopy_ecf(BIO_ISOP,it2m),canopy_ecf(BIO_TERP,it2m) /) )
      call datewrite("DBIO EISOP EMTP EMTL ESOIL ", (/  E_ISOP, &
             E_MTP, E_MTL, SoilNOx(i,j) /) ) 
      call datewrite("DBIO EmisNat ", EmisNat(:,i,j) )

     end if


  end subroutine setup_bio

  !----------------------------------------------------------------------------


   subroutine Set_SoilNOx()
      integer :: i, j, nLC, iLC, LC
      logical :: my_first_call = .true.
      real    :: f, ft, fn, ftn
      real    :: enox, enh3  ! emissions, ugN/m2/h
      real :: beta, bmin, bmax, bx, by ! for beta function
      real :: hfac

      if( DEBUG_SOILNOX .and. debug_proc ) then
         write(*,*)"Biogenic_ml DEBUG_SOILNOX EURO: ",&
          current_date%day, current_date%hour, current_date%seconds,&
          USE_EURO_SOILNOX, EURO_SOILNOX_DEPSCALE
      end if

      if ( .not. USE_EURO_SOILNOX  ) return ! and fSW has been set to 1. at start


      ! We reset once per hour

      if ( current_date%seconds /= 0 .and. .not. my_first_call ) return
      hfac = 0.5 ! Lower at night
      if ( current_date%hour > 7 .and. current_date%hour < 20 ) hfac = 1.5


        do j = 1, ljmax
           do i =  1, limax

             nlc = LandCover(i,j)%ncodes

           ! Temperature function from Rolland et al., 2005, eqn. 6

             ft =  exp( (t2_nwp(i,j,1)-273.15-20)*log(2.1) / 10.0 )

           ! Inspired by e.g. Pilegaard et al, Schaufler et al. (2010)
           ! we scale emissions from seminat with N-depositions
           ! We use a factor normalised to 1.0 at 5000 mgN/m2/a

             fn = AnnualNdep(i,j)/5000.0 ! scale for now
             fn = fn * EURO_SOILNOX_DEPSCALE  ! See ModelConstants_ml

             ftn = ft * fn * hfac 

             enox = 0.0
             enh3 = 0.0 

     LCLOOP: do ilc= 1, nLC

                 LC = LandCover(i,j)%codes(ilc)
                 if ( LandType(LC)%is_water ) cycle
                 if ( LandType(LC)%is_ice   ) cycle
                 if ( LandType(LC)%is_iam   ) cycle

               ! Soil NO
               ! for 1 ugN/m2/hr, the temp funct alone would give
               ! ca. 6 mgN/m2/a in Germany, where dep is about 5000 mgN/m2 max ca. 9
               ! Conif Forests in Germany  

                 f  = LandCover(i,j)%fraction(ilc) 
                 beta = 0.0

                 if ( LandType(LC)%is_conif ) then
                    enox = enox + f*ftn*150.0
                    enh3 = enh3 + f*ftn*1500.0 ! Huge?! W+E ca. 600 ngNH3/m2/s -> 1800 ugN/m2/h
                 else if ( LandType(LC)%is_decid ) then
                    enox = enox + f*ftn* 50.0
                    enh3 = enh3 + f*ftn*500.0 !  Just guessing
                 else if ( LandType(LC)%is_seminat ) then
                    enox = enox + f*ftn* 50.0
                    enh3 = enh3 + f * ftn  *20.0 !mg/m2/h approx from US report 1 ng/m2/s

                 else if ( LandType(LC)%is_crop    ) then ! emissions in 1st 70 days

                    bmin = Landcover(i,j)%SGS(iLC) -30 ! !st March LandCover(i,j)%SGS(iLC) - 30 
                    bmax = Landcover(i,j)%SGS(iLC) +30 ! End April  LandCover(i,j)%SGS(iLC) + 40 

                    ! US p.29, Suttn had ca. 20 ng/m2/s  = 60ugN/m2/hfor crops
                    ! throughout growing season
                    if ( daynumber >= Landcover(i,j)%SGS(iLC) .and. &
                         daynumber <= Landcover(i,j)%EGS(iLC) ) then
                         enox = enox + f* 1.0
                         enh3 = enh3 + f * 60.0
                    end if

                    ! CRUDE - just playing for NH3.
                    ! NH3 from fertilizer? Assume e.g. 120 kg/ha over 1 month
                    ! with 10% giving emission, i.e. 10 kg/ha
                    ! 10 kg/ha/month =  ca. 1000 ugN/m2/h

                    ! For NO, numbers based upon papers by e.g. Rolland,
                    ! Butterbach, etc.
                    if ( daynumber >= bmin .and. daynumber <= bmax ) then

                         bx = (daynumber-bmin)/( bmax-bmin)
                         bx = max(bx,0.0)
                         by = 1.0 - bx
                         beta =  ( bx*by *4.0) 
                         enox = enox + f*80.0*ft* beta 
                         enh3 = enh3 + f * 1000.0*ft * beta
                    end if

                    
                 end if
                 if (  DEBUG_SOILNOX .and. debug_proc .and. &
                     i == debug_li .and. j == debug_lj ) then
                   write(*, "(a,4i4,f7.2,9g12.3)") "LOOPING SOIL", daynumber, &
                   iLC, LC, LandCover(i,j)%SGS(iLC), t2_nwp(i,j,1)-273.15, &
                      f, ft, fn, ftn,  beta, enox, enh3
                   if(iLC==1) &
                     call datewrite("HFAC SOIL", (/ 1.0*daynumber,hfac /) )
                 end if
                 enox = max( 0.001, enox ) ! Just to stop negatives while testing
    
               ! Soil NH3
           end do LCLOOP


     ! And we scale EmisNat to get units kg/m2 consistent with
     ! Emissions_ml (snapemis).  ug/m2/h -> kg/m2/s needs 1.0-9/3600.0. 
 
           SoilNOx(i,j) = enox
           SoilNH3(i,j) = enh3
 
         end do
      end do
 !     if ( DEBUG_SOILNOX .and. debug_proc ) then
 !             SoilNOx(:,:) = 1.0 ! Check scaling
 !     end if

      if ( DEBUG_SOILNOX .and. debug_proc ) then
         i = debug_li
         j = debug_lj
         write(*,"(a,4i4)") "RESET_SOILNOX: ",  1, limax, 1, ljmax
         write(*,"(a,2i4,2f12.4,es12.4)") "RESET_SOILNOX: ", &
                 daynumber, current_date%hour, t2_nwp(i,j,1), SoilNOx(i,j), AnnualNdep(i,j)
      end if

      my_first_call = .false.

   end subroutine Set_SoilNOx
end module Biogenics_ml
