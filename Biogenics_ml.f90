! <Biogenics_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2011 met.no
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
  !    Dave Simpson, 2010-2011
  !---------------------------------------------------------------------------

  use CheckStop_ml,      only: CheckStop
  use GridValues_ml    , only : i_fdom,j_fdom, debug_proc,debug_li,debug_lj
  use Io_ml            , only : IO_FORES, open_file, ios, PrintLog, datewrite
  use KeyValue_ml,       only : KeyVal,KeyValue
  use LandDefs_ml,       only: LandType, LandDefs
  use LandPFT_ml,        only: MapPFT_LAI, pft_lai
  use Landuse_ml,        only : LandCover
  use LocalVariables_ml, only : Grid  ! -> izen, DeltaZ
  use ModelConstants_ml, only : NPROC, MasterProc, TINY, &
                           USE_PFT_MAPS, NLANDUSEMAX, IOU_INST, & 
                           KT => KCHEMTOP, KG => KMAX_MID, & 
                           DEBUG_BIO, BVOC_USED, MasterProc
  use NetCDF_ml,        only : ReadField_CDF, printCDF
  use OwnDataTypes_ml,  only : Deriv, TXTLEN_SHORT
  use Par_ml   , only :  MAXLIMAX,MAXLJMAX,MSG_READ1,me, limax, ljmax
  use PhysicalConstants_ml,  only :  AVOG, GRAV
  use Radiation_ml,          only : PARfrac, Wm2_uE
  use Setup_1dfields_ml,     only : rcbio  
  use SmallUtils_ml, only : find_index
  use TimeDate_ml,       only :  current_date
  implicit none
  private

  !/-- subroutines
  public ::  Init_BVOC
  private :: Get_LCinfo
  public ::  GetEuroBVOC
  private :: MergedBVOC
  public ::  setup_bio
  public ::  SetDailyBVOC
  private :: TabulateECF

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
  integer, public, parameter :: N_ECF=2, ECF_ISOP=1, ECF_TERP=2
  integer, public, parameter :: NBIO_DEF=3, BIO_ISOP=1, BIO_MTP=2, BIO_MTL=3
  integer, public, parameter :: BIO_TERP=2 ! Used for final emis, sum of MTP+MTL
  integer, public, save ::  last_bvoc_LC   !max index land-cover with BVOC (min 4)

 ! Set true if LCC read from e.g. EMEP_EuroBVOC.nc:
 ! (Currently for 1st four LCC, CF, DF, BF, NF)
  logical, private, dimension(NLANDUSEMAX), save :: HaveLocalEF 

  real, public, save, dimension(MAXLIMAX,MAXLJMAX,size(BVOC_USED)) :: &
     EmisNat =0.0      !  will be transferred to d_2d emis sums


  !standard emission factors (EFs) per LC
  !Need to dimension later for Emtp, Emtl, last_bvoc_LC
  real, private, save, allocatable, dimension(:,:,:,:) :: &
     bvocEF       !  Gridded std. emissions per PFT

  !standard emission factors per LC for daily LAI
  real, private, save, dimension(MAXLIMAX,MAXLJMAX,size(BVOC_USED)) :: &
     day_embvoc = 0.0   !  emissions scaled by daily LAI

  logical, private, dimension(MAXLIMAX,MAXLJMAX) :: EuroMask

  !/-- Canopy environmental correction factors-----------------------------
  !
  !    - to correct for temperature and light of the canopy 
  !    - from Guenther's papers. (Limit to 0<T<40 deg C.)

  real, public, save, dimension(N_ECF,40) :: canopy_ecf  ! Canopy env. factors
                                                        
  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    subroutine Init_BVOC()

!    Read natural BVOC emission potentials
!-----------------------------------------------------------------------------
!   Emission potentials (EFs) now read a netcdf of ca. 50x50 resolution.
!   This fie, EMEP_EuroBVOC.nc uses the ICP-Forests species map as processed
!   by Renate Koeble at JRC (e.g.
!   EFs now a mixure of rates from Simpson et al., 1999, JGR, Vol 104, D7, 
!   8113-8152, and Keenan, T. et al., ACP, 2009, 9, 4053-4076 

    integer :: alloc_err

      if ( size(BVOC_USED) == 0 ) then
        call PrintLog("No Biogenic Emissions ", MasterProc)
        return
      end if

   !====================================
   ! Initialise chemical emission rates to zero. In rest of  code
   ! only k=KMAX_MID will be set

     rcbio(:,:) = 0.0  
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

    real    :: loc(MAXLIMAX,MAXLJMAX) = 0.0  ! Emissions read from file
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
             loc,1,interpol='zero_order',needed=.true.,debug_flag=.false.)
         if( debug_proc ) write(*, "(2a,f12.3,3i2)") "EURO-BVOC:E ", &
             trim(varname), loc(debug_li, debug_lj), iVeg, ibvoc, iEmis
         if( debug_proc ) write(*, "(2a,2es12.3)") "EURO-BVOC:minmax ", &
             trim(varname), minval(loc), maxval(loc)
         bvocEF(:,:,ibvoc,iEmis) = loc(:,:)
       end do

      ! Make a mask where we can use the local bvoc. Should be the same from
      ! all EFs, since only non-def areas set to -999, otherwise zero or +
      ! If any values exist, should exist for all entries, hence check.

       if( iVeg == 1 )  then
          where(loc>-1.0)
            EuroMask = .true.
          end where
       else  ! Just check that following maps are consistent
           do i=1,MAXLIMAX
           do j=1,MAXLJMAX
             if ( EuroMask(i,j) .and. loc(i,j)<0.0 ) then
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
         if ( my_first_call  ) then ! print out 1st day
              call printCDF("BIO-LAI", workarray, "m2/m2" )
              workarray(:,:) = day_embvoc(:,:,1)
              call printCDF("BIO-Eiso", workarray, "ug/m2/h" )
              workarray(:,:) = day_embvoc(:,:,2) + day_embvoc(:,:,3)
              call printCDF("BIO-Emt", workarray, "ug/m2/h" )
              deallocate(  workarray )
         end if
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
  real :: E_ISOP , E_MTP, E_MTL

! To get from ug/m2/h to molec/cm3/s
! ug -> g  1.0e-9; g -> mole / MW; x AVOG
  real, save :: biofac_ISOP = 1.0e-12*AVOG/64.0/3600.0   ! needs /Grid%DeltaZ
  real, save :: biofac_TERP = 1.0e-12*AVOG/136.0/3600.0  ! needs /Grid%DeltaZ

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

  rcbio(:,KG) = 0.0

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
      rcbio(BIO_ISOP,KG)   = E_ISOP * biofac_ISOP/Grid%DeltaZ
      EmisNat(i,j,BIO_ISOP)= E_ISOP * 1.0e-9/3600.0

  else ! night
     EmisNat(i,j,BIO_ISOP) = 0.0
     E_MTL = 0.0
     E_ISOP = 0.0
  endif ! daytime

    if ( BIO_TERP > 0 ) then

     ! add pool-only terpenes rate;
        E_MTP = day_embvoc(i,j,BIO_MTP)*canopy_ecf(ECF_TERP,it2m)
        rcbio(BIO_TERP,KG)    = (E_MTL+E_MTP) * biofac_TERP/Grid%DeltaZ
        EmisNat(i,j,BIO_TERP) = (E_MTL+E_MTP) * 1.0e-9/3600.0
    end if
 
    if ( DEBUG_BIO .and. debug_proc .and. i==debug_li .and. j==debug_lj .and. &
         current_date%seconds == 0 ) then

      call datewrite("DBIO env ", it2m, (/ max(par,0.0), max(cL,0.0), &
            canopy_ecf(BIO_ISOP,it2m),canopy_ecf(BIO_TERP,it2m) /) )
      call datewrite("DBIO EISOP EMTP EMTL ", (/  E_ISOP, E_MTP, E_MTL /) ) 
      call datewrite("DBIO rc ", (/ rcbio(BIO_ISOP,KG), rcbio(BIO_TERP,KG) /) )
      call datewrite("DBIO EmisNat ", EmisNat(i,j,:) )

     end if


  end subroutine setup_bio

  !----------------------------------------------------------------------------
end module Biogenics_ml
