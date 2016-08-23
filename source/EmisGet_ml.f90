! <EmisGet_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-201409 met.no
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
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                     module EmisGet_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  use CheckStop_ml,      only: CheckStop, StopAll
!CMR  use ChemSpecs_adv_ml,  only: NSPEC_ADV ! max possible number in split 
!CMR  use ChemSpecs_tot_ml,  only: NSPEC_TOT 
!CMR  use ChemChemicals_ml,  only: species
  use ChemSpecs,         only: NSPEC_ADV, NSPEC_TOT, species 
  use Country_ml,        only: NLAND, IC_NAT, IC_VUL, Country, &
                               ! NMR-NH3 specific variables (hb NH3Emis)
                               IC_NMR 
  use EmisDef_ml,        only: NSECTORS, ANTROP_SECTORS, NCMAX, FNCMAX, & 
                               NEMIS_FILE, EMIS_FILE, & 
                               ISNAP_SHIP, ISNAP_NAT, VOLCANOES_LL, &
                               ! NMR-NH3 specific variables (for FUTURE )
                               NH3EMIS_VAR,dknh3_agr,ISNAP_AGR,ISNAP_TRAF, &
                               NROADDUST
  use GridAllocate_ml,   only: GridAllocate
  use GridValues_ml, only: debug_proc,debug_li,debug_lj, i_fdom, j_fdom !cdfemis
  use GridValues_ml, only: glon, glat, A_bnd, B_bnd
  use Io_ml,             only: open_file, NO_FILE, ios, IO_EMIS, &
                               Read_Headers, read_line, PrintLog
  use KeyValueTypes,       only: KeyVal
  use ModelConstants_ml, only: NPROC, TXTLEN_NAME, &
                               DEBUG,  KMAX_MID,KMAX_BND, Pref,&
              SEAFIX_GEA_NEEDED, & ! only if emission problems over sea
                               MasterProc,DEBUG_GETEMIS,DEBUG_ROADDUST,USE_ROADDUST
  use NetCDF_ml, only  : ReadField_CDF  !CDF_SNAP

  use Par_ml,            only: MAXLIMAX, MAXLJMAX, limax, ljmax, me
  use SmallUtils_ml,     only: wordsplit, find_index
  use Volcanos_ml
  use netcdf
  use NetCDF_ml, only  : check

  implicit none
  private

 ! subroutines:

  public  :: EmisGet           ! Collects emissions of each pollutant
  public  :: EmisGetCdf        ! cdfemis
  public  :: EmisSplit         ! => emisfrac, speciation of voc, pm25, etc.
  public  :: EmisHeights       ! => nemis_kprofile, emis_kprofile
                               !     vertical emissions profile
  public  :: RoadDustGet       ! Collects road dust emission potentials
  public :: femis             ! Sets emissions control factors 
  private :: CountEmisSpecs    !


  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
!  logical, private, save :: my_first_call = .true.
  logical, private, save :: my_first_road = .true.

  ! e_fact is the emission control factor (increase/decrease/switch-off)
  ! e_fact is read in from the femis file and applied within EmisGet
  real, private, save, &
         dimension(NSECTORS,NLAND,NEMIS_FILE)  :: e_fact 

  ! emisfrac is used at each time-step of the model run to split
  ! emissions such as VOC, PM into species. 

  integer, public, parameter :: NMAX = NSPEC_ADV 
  integer, public, save :: nrcemis, nrcsplit
  integer, public, dimension(NEMIS_FILE) , save :: emis_nsplit
  real, public,allocatable, dimension(:,:,:), save :: emisfrac
  integer, public,allocatable, dimension(:), save :: iqrc2itot
  integer, public, dimension(NSPEC_TOT), save :: itot2iqrc
  integer, public, dimension(NEMIS_FILE), save :: Emis_MolWt
  real, public,allocatable, dimension(:), save :: emis_masscorr
  real, public,allocatable, dimension(:), save :: roaddust_masscorr

  ! vertical profiles for SNAP emis, read from EmisHeights.txt
  integer, public, save :: nemis_kprofile
  real, public,allocatable, dimension(:,:), save :: emis_kprofile
  real, public,allocatable, dimension(:,:), save :: emis_hprofile

  ! some common variables
  character(len=80), private :: fname             ! File name
  character(len=80), private :: errmsg

 ! Import list of the emitted species we need to find in the 
 ! emissplit files.
  include 'CM_EmisSpecs.inc'
  logical, dimension(NEMIS_SPECS) :: EmisSpecFound = .false.

  !CDF cdfemis tests
   real, public, save,  dimension(NLAND,NEMIS_FILE) ::&
      sumcdfemis ! Only used for MasterProc
   real, allocatable, private, save,  dimension(:,:) :: cdfemis
   integer, allocatable, public, save,  dimension(:,:) :: nGridEmisCodes
   integer, allocatable, public, save,  dimension(:,:,:):: GridEmisCodes
   real, allocatable, public, save,  dimension(:,:,:,:,:):: GridEmis
 contains


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine EmisGetCdf(iem, fname,incl,excl)
   integer, intent(in) :: iem ! index in EMIS_FILE array and GridEmis output
   character(len=*), intent(in)    :: fname
   character(len=*),dimension(:), optional :: &
       incl, excl ! Arrays of cc to inc/exclude
   integer :: i,j, ic, isec, allocerr(6), icode, status
   real, dimension(NLAND) :: sumcdfemis_loc, sumcdfemis_iem
   integer :: icc, ncc
   character(len=40) :: varname
   integer, save :: ncmaxfound = 0 ! Max no. countries found in grid
   integer, save :: ncalls=0
   integer :: ncFileID, nDimensions,nVariables,nAttributes,timeDimID,varid,&
           xtype,ndims  !TESTE testing
   character(len=10) :: ewords(7), code ! Test Emis:UK:snap:7
   integer :: nwords, err
   logical, save :: my_first_call = .true.

!
   if(MasterProc)  print *, "ME INTO EMISGETCDF ", me, trim(fname)&
        ,present(incl), present(excl)  ! optionals
   if( present(incl) ) then
      if(MasterProc) write(*,"(a,i4,99a4)") trim(fname)//":INCL=", size(incl), incl
   end if
   if( present(excl) ) then
      if(MasterProc) write(*,"(a,i4,99a4)") trim(fname)//":EXCL=", size(excl), excl
   end if

   if( my_first_call  )  then
        allocerr = 0
        allocate(cdfemis(MAXLIMAX,MAXLJMAX),stat=allocerr(1))
        allocate(nGridEmisCodes(MAXLIMAX,MAXLJMAX),stat=allocerr(2))
        allocate(GridEmisCodes(MAXLIMAX,MAXLJMAX,NCMAX),stat=allocerr(3))
        allocate(GridEmis(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILE),&
            stat=allocerr(4))
        call CheckStop(any ( allocerr(:) /= 0), &
              "EmisGet:Allocation error for cdfemis")
        GridEmisCodes = -1   
        nGridEmisCodes = 0
        GridEmis = 0.0
        sumcdfemis = 0.0

        my_first_call = .false.

   end if

    ! Initialise sums for this pollutant and processor
    sumcdfemis_loc = 0.0
    sumcdfemis_iem = 0.0

!---------------------------------------------------------------
! find emis file and  main properties

 !HUNT call check(nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncFileID))
 status=nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncFileID)
if( index(fname, "Ship")>0 ) print *, me, " CDFHUNTTOP ", trim(fname)
 if( status /= nf90_noerr ) then
  if( MasterProc ) print *, "EmisGetCdf - couldn't open "//trim(fName)
   return
 end if

 call check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,timeDimID))
 if( MasterProc .or. index(fname, "Ship")>0 ) then
  write( *,*) 'Nb of global attributes: ',nAttributes
  write( *,*) 'EmisGetCdf '//trim(fName),' properties: '
  write( *,*) 'EmisGetCdf Nb of dimensions: ',nDimensions
  write( *,*) 'EmisGetCdf Nb of variables: ',nVariables
 end if


  do varid=1,nVariables
     call check(nf90_Inquire_Variable(ncFileID,varid,varname,xtype,ndims))

    ! Emission terms look like, e.g. Emis:FR:snap:7
     if( index( varname, "Emis:") < 1 ) cycle ! ONLY emissions wanted
     call wordsplit(varname,4,ewords,nwords,err,separator=":")
     if( ewords(3) /= "snap" ) cycle ! ONLY SNAP coded for now

     code =ewords(2)    ! Country ISO
     read(ewords(4),"(I6)") isec

     !/ We can include or exclude countries:
      if ( present(excl) ) then
           !if(MasterProc) print "(a,50a4)", "INCEXCTEST-E "// &
           !    trim(code)!, (trim(excl(i)),i=1,size(excl))
        if ( find_index( code, excl ) >0 ) then
           !if(MasterProc) print *, "INCEXCTEST-XX ", &
           !    trim(code), find_index( code, excl )
           if(MasterProc) write(*,"(a,a,i5)") "EmisGetCdf"//trim(fname)// &
               "exc-exludes:", trim(code), find_index( code, excl )
           cycle
        end if
      end if
      if ( present(incl) ) then
        if ( find_index( code, incl ) <1 ) then
           if(MasterProc) write(*,"(a,a,i5)") "EmisGetCdf"//trim(fname)// &
               "inc-exludes:", trim(code) !, find_index( code, excl )
           cycle
        end if
      end if

      !if(MasterProc) write(*,"(3a,i5)") "EmisGetCdf"//trim(fname),  &
      !         "==>includes:", trim(code), find_index( code, incl ), find_index( code, excl )
      !if(MasterProc) print *, "!!EmisGetCdf"//trim(fname),  &
      !         "includes:", trim(code), find_index( code, excl )

     ic = find_index( code, Country(:)%code )  !from Country_ml
     !if(MasterProc) print "(a)", "INCEXCTEST-A "//trim(fname)//":"//trim(code)

     if ( Country(ic)%code == "N/A" ) then
          if(MasterProc) print *, "CDFCYCLE ", ic, trim(Country(ic)%code)
          cycle    ! see Country_ml
     end if
     call CheckStop( ic < 1 , "CDFEMIS NegIC:"//trim(fname)//":"//trim(code) )

     !if( DEBUG_GETEMIS .and. debug_proc ) write( *,*) 'EmisGetCdf ', trim(fname), varid,trim(varname), &
     !      " ", trim(code), ic, isec 
     if( DEBUG_GETEMIS .and. debug_proc ) write(*,"(2a,i6,a,2i4)") 'EmisGetCdf ', &
           trim(fname), varid," "//trim(varname)// " "// trim(code), ic, isec 


     cdfemis = 0.0 ! safety, shouldn't be needed though
!if( index(fname, "Ship")>0 ) print *, me, " CDFHUNTIN ", trim(fname), trim(varname)
     call ReadField_CDF(fname,varname,cdfemis,1,&
!??              known_projection='longitude latitude', &
               interpol='mass_conservative',&
                needed=.false.,UnDef=0.0,debug_flag=.false.)
!if( index(fname, "Ship")>0 ) print *, me, " CDFHUNTUT ", trim(fname), trim(varname), sum(cdfemis)

     if( maxval(cdfemis ) < 1.0e-10 ) cycle ! Likely no emiss in domain

     if ( DEBUG_GETEMIS .and. debug_proc ) then
     !if ( DEBUG_GETEMIS .and. MasterProc ) then !me has some sea usually !!SHIPHUNT
         ncalls = ncalls + 1
         write(*,"(3a,3i4,i3,9f12.4)")"CDF emis-in ",&
           trim(fname)//":", &
           trim(Country(ic)%name)//":"//trim(EMIS_FILE(iem)), ic, &
            i_fdom(debug_li),j_fdom(debug_lj),&
             isec, e_fact(isec,ic,iem), &
                cdfemis(debug_li, debug_lj), maxval(cdfemis)
      end if

!call MPI_BARRIER(MPI_COMM_WORLD, INFO)

      sumcdfemis_loc(ic) = sumcdfemis_loc(ic) + &
          0.001 *  e_fact(isec,ic,iem) * sum( cdfemis(:,:) )
!print "(a,4i4,3es12.3)", "CDFSUMing ", me, isec, ic, iem, e_fact(isec,ic,iem), sumcdfemis_loc(ic), sum( cdfemis(:,:) )

        do j = 1, ljmax
            do i = 1, limax
              if( cdfemis(i,j) > 1.0e-10 ) then
 !if( i_fdom(i) < 33 .and. j_fdom(j) < 15 .and. ic == 60 ) then
 ! Some bit of Germany found in Africa in MACC2 data. Correct to 
 if( glat(i,j) < 35 .and. j_fdom(j) < 15 .and. ic == 60 ) then
 !if( i_fdom(i) < 33 .and. j_fdom(j) < 15 ) then
        write(*,"(3a,3i4,2f8.3,i3,9f12.4)")"AFRICA:",&
           trim(fname)//":", &
           trim(Country(ic)%name)//":"//trim(EMIS_FILE(iem)), ic, &
            i_fdom(i),j_fdom(j),&
             glon(i,j), glat(i,j), & 
             isec, e_fact(isec,ic,iem), &
                cdfemis(i, j), maxval(cdfemis)
         ic = 224 ! Change to North Africa from IIASA system
 end if !========== AFRICA
                call GridAllocate("GridCode"// trim ( EMIS_FILE(iem) ),&
                  i,j,ic,NCMAX, icode, ncmaxfound,GridEmisCodes,nGridEmisCodes,&
                  debug_flag=.false.)
                 GridEmis(isec,i,j,icode,iem) = &
                   GridEmis(isec,i,j,icode,iem) + cdfemis(i,j) * e_fact(isec,ic,iem)
                !if ( debug_proc .and. i == debug_li .and. j == debug_lj ) then

                 ncc = nGridEmisCodes(i,j)
                if ( DEBUG_GETEMIS .and. &
                !if ( &
                        debug_proc .and. i==debug_li .and. j==debug_lj ) then
                    write(*,"(a,7i4,2es12.3,3x,99i3)")"CDF emis-alloc ",&
                     i_fdom(i),j_fdom(j), isec, ic, icode, ncc, ncmaxfound, &
                     cdfemis(i,j), GridEmis(isec,i,j,icode,iem), &
                     (GridEmisCodes(i,j,icc),icc=1,ncc)
                end if
              end if
            end do !i
        end do !j
      end do ! variables 

     call check(nf90_close(ncFileID))

     CALL MPI_REDUCE(sumcdfemis_loc(:),sumcdfemis_iem(:),NLAND,&
                          MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,INFO)

    if( debug_proc ) write(*,*) "FINISHED CDF "//trim(fname), me, sum(sumcdfemis_loc)
    write(*,*) "FINISHED CDFLOC "//trim(fname), me, sum(sumcdfemis_loc)
    if( MasterProc ) then
       write(*,*) "FINISHED CDFM0 "//trim(fname), me, iem, sum(sumcdfemis_loc)
       sumcdfemis(:,iem) = sumcdfemis(:,iem) + sumcdfemis_iem(:)
       do ic = 1, NLAND
         if ( sumcdfemis(ic,iem) > 1.0e-10 ) & 
            write(*,"(a,i5,f12.3)") "CDFSUM "//trim(fname), ic, sumcdfemis(ic,iem)
       end do
    end if

  end subroutine EmisGetCdf
!-----------------------------------------------------------------------------
  subroutine EmisGet(iemis,emisname,IRUNBEG,JRUNBEG,GIMAX,GJMAX, &
                     globemis,globnland,globland,sumemis,        &
                     globemis_flat,flat_globnland,flat_globland)

!.......................................................................
!  DESCRIPTION:
!  Reads in emissions from one file, specified by iemis. 
!  The arrays read in here are the global arrays (allocatable)
!.......................................................................

  !--arguments
  integer, intent(in)                     :: iemis     ! emis index 
  character(len=*), intent(in)            :: emisname  ! emission name
  integer, intent(in) :: IRUNBEG,JRUNBEG,GIMAX,GJMAX   ! domain limits
  real,    intent(out), dimension(:,:,:,:):: globemis      ! Emission values
  real,    intent(out), dimension(:,:,:)  :: globemis_flat ! Flat emissions
                                                           ! (e.g. shipping)
  integer, intent(inout), dimension(:,:,:)::   &
                                 globland,     & ! Codes of countries-emitters
                                 flat_globland   ! Flat emis.codes (shipping)
  integer, intent(inout), dimension(:,:)  ::   &
                                 globnland,    & ! No. emitions in grid
                                 flat_globnland  ! No. flat emitions in grid
  real,    intent(inout), dimension(:,:)  :: sumemis ! Emission sums per 
                                                     ! country(after e_fact) 

  !--local
  integer :: flat_iland,                &
             i, j, isec, iland,         &  ! loop variables
             iic,ic                        ! country code (read from file)
  real    :: duml,dumh                     ! dummy variables, low/high emis.
  real, dimension(NSECTORS)  :: tmpsec     ! array for reading emission files
  integer, save :: ncmaxfound = 0          ! Max no. countries found in grid
  integer, save :: flat_ncmaxfound = 0     ! Max no. countries found in grid
                                           ! including flat emissions


   !>============================
!rv4_2.1 
!rv4_2.1    if ( my_first_call ) then
!rv4_2.1         sumemis(:,:) =  0.0       ! initialize sums
!rv4_2.1         ios = 0
!rv4_2.1         call femis()              ! emission factors (femis.dat file)
!rv4_2.1         if ( ios /= 0 )return
!rv4_2.1         my_first_call = .false.
!rv4_2.1    endif
  !>============================


      globemis   (:,:,:,:) = 0.0
      globemis_flat(:,:,:) = 0.0

      !if (DEBUG_GETEMIS) write(unit=6,fmt=*) "Called EmisGet with index, name", &
      !     iemis, trim(emisname)
      write( *,*) "Called EmisGet with index, name", &
           iemis, trim(emisname)
      fname = "emislist." // emisname
      call open_file(IO_EMIS,"r",fname,needed=.true.)
      call CheckStop(ios,"EmisGet: ios error in emission file")

      if (trim ( emisname ) == "nh3" ) dknh3_agr=0.0 ! NH3Emis experimental
 

READEMIS: do   ! ************* Loop over emislist files *******************

            read(unit=IO_EMIS,fmt=*,iostat=ios) iic,i,j, duml,dumh,  &
                                    (tmpsec(isec),isec=1,NSECTORS)

            if( DEBUG_GETEMIS .and. i==DEBUG%IJ(1) .and. j==DEBUG%IJ(2) ) write(*,*) &
                "DEBUG GetEmis "//trim(emisname) // ":" , iic, duml,dumh
            !if( j== DEBUG%IJ(2) ) write(*,*) &
            !    "DEBUG GetEmis "//trim(emisname) // ":" , iic, duml,dumh
            if ( ios <  0 ) exit READEMIS            ! End of file
            call CheckStop(ios > 0,"EmisGet: ios error in emission file")

            ! Check if country code in emisfile (iic) is in the country list
            ! from Countries_ml, i.e. corresponds to numbering index ic

            do ic=1,NLAND
               if((Country(ic)%index==iic))&
                    goto 543
            enddo 
            write(unit=errmsg,fmt=*) &
                   "COUNTRY CODE NOT RECOGNIZED OR UNDEFINED ", iic
            call CheckStop(errmsg)
            ic=0
543         continue

            i = i-IRUNBEG+1     ! for RESTRICTED domain
            j = j-JRUNBEG+1     ! for RESTRICTED domain

            if ( i  <=  0 .or. i  >  GIMAX .or.   & 
                 j  <=  0 .or. j  >  GJMAX .or.   &
                 ic <=  0 .or. ic >  NLAND .or.   &
                 ic == IC_NAT              .or.   &  ! Excludes DMS
                 (ic == IC_VUL .and. VOLCANOES_LL) )&! Excludes Volcanoes 
                                                     ! from gridSOx. Read from
                                                     ! VolcanoesLL.dat instead
             cycle READEMIS

             ! Ship emissions

             if ( Country(ic)%is_sea ) then    ! ship emissions

              ! ..........................................................
              ! Generate new land allocation in 50 km grid for FLAT 
              ! EMISSIONS (ships). First, we check if country "ic"
              ! has already  been found within that grid. If not, then ic is
              ! added to flat_landcode and flat_nlandcode increased by one.
   
              ! Test that ship emissions are only in sector ISNAP_SHIP
               do isec=1,(ISNAP_SHIP-1) 
                  if ( MasterProc.and.tmpsec(isec) /= 0) &
                     write(*,"(a,3i4,i3,f12.4)")"SEA"//trim(emisname), &
                        iic,i,j,isec,tmpsec(isec)

                  call CheckStop(tmpsec(isec) /= 0,  &
                   "EmisGet: NOT FLAT EMISSIONS - check SEAFIX_GEA_NEEDED comments in ModelConstants_ml")
               enddo
               do isec=ISNAP_SHIP+1,NSECTORS
                  if ( MasterProc.and.tmpsec(isec) /= 0) &
                     write(*,"(a,3i4,i3,f12.4)")"SEA"//trim(emisname), &
                        iic,i,j,isec,tmpsec(isec)
                  !call CheckStop(tmpsec(isec) /= 0,  &
                  !      "EmisGet: NOT FLAT EMISSIONS")
               enddo
              ! end test

               call GridAllocate("FLat",i,j,ic,FNCMAX, flat_iland, &
                   flat_ncmaxfound,flat_globland,flat_globnland)
              ! ...................................................
              ! Assign e_fact corrected emissions to global FLAT 
              ! emission matrices.
              ! ...................................................

               globemis_flat(i,j,flat_iland) = globemis_flat(i,j,flat_iland) &
                     + e_fact(ISNAP_SHIP,ic,iemis) * tmpsec(ISNAP_SHIP)

              !......................................................
              !        Sum over all sectors, store as Ktonne:
              !......................................................

                 sumemis(ic,iemis) = sumemis(ic,iemis)   &
                      + 0.001 *  e_fact(ISNAP_SHIP,ic,iemis) * tmpsec(ISNAP_SHIP)

                cycle READEMIS
             endif !ship emissions            


              !.......................................................
              !  Volcanos
              !.......................................................

              if ( trim ( emisname ) == "sox" ) then
                if (ic == IC_VUL) then
                  volc_no=volc_no+1
                  if (DEBUG_GETEMIS) write(*,*)'Volcano no. ',volc_no
                  i_volc(volc_no)=i
                  j_volc(volc_no)=j

                  emis_volc(volc_no) = tmpsec(ISNAP_NAT) *   &
                                       e_fact(ISNAP_NAT,IC_VUL,iemis)
                  nvolc=volc_no

                  call CheckStop(nvolc>NMAX_VOLC,"EMISGET, nvolc>NMAX_VULC")

                  write(*,*)'Found ',nvolc,' volcanoes on sox file'

                  sumemis(IC_VUL,iemis) = sumemis(IC_VUL,iemis)        &
                                          + 0.001 * emis_volc(volc_no)
                   cycle READEMIS ! do not want to count volcano "landcode"
                endif ! ic
             endif ! so2

             !..............................................................
             ! end Volcanoes
             !..............................................................


             !  For VOC natural and agricultur emissions (managed forests) 
             !  set to  zero

             if ( trim ( emisname ) == "voc" ) tmpsec(11:11) = 0.0
  
            ! NH3emis (FUTURE/EXPERIMENTAL for NMR-NH3 project)
            ! For NH3 activity data, set 'static emissions' to zero
            ! For northwestern Europe, read in Sector_NH3Emis.txt in run.pl
            ! Traffic emis are zero in the Danish emissions
             if (NH3EMIS_VAR .and.  trim ( emisname ) == "nh3" .and.&
                 ic == IC_NMR) then
                dknh3_agr=dknh3_agr+ tmpsec(ISNAP_AGR)+tmpsec(ISNAP_TRAF)
                tmpsec(ISNAP_AGR:ISNAP_AGR) = 0.0
                tmpsec(ISNAP_TRAF:ISNAP_TRAF) = 0.0
             endif


             ! ..........................................................
             ! generate new land allocation in 50 km grid. First, we check if
             ! country "ic" has already  been found within that grid. If not,
             ! then ic is added to landcode and nlandcode increased by one.

              call GridAllocate("SNAP"// trim ( emisname ),i,j,ic,NCMAX, &
                                 iland,ncmaxfound,globland,globnland)

              globemis(:,i,j,iland) = globemis(:,i,j,iland) &
                        + e_fact(:,ic,iemis) *  tmpsec(:)


             ! Sum over all sectors, store as Ktonne:

              sumemis(ic,iemis) = sumemis(ic,iemis)   &
                                  + 0.001 * sum (e_fact(:,ic,iemis)*tmpsec(:))

!rb: Old version (below) does not work if same grid point occurs several times in emis-file
!    Probably the same problem for ship emissions etc. Not changed now!
!     + 0.001 * sum (globemis (:,i,j,iland))
              
        end do READEMIS 
        
        close(IO_EMIS)
        ios = 0
  end subroutine EmisGet


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine femis()
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!-------------------------------------------------------------------------
! Read emission control factors for the emissions from (optional) file femis. 
! Emission factors are applied EITHER to specified country and/or
! emission sector, e.g. the femis file with input:
!      Code  5  co  sox  nox    nh3   voc
!      27    6  1.0  0.5  1.0   1.0   1.0
! will reduce SO2 emissions from sector 6 by a factor 0.5 for the UK 
! (country 27);
! OR to all countries/sectors: a zero for country or sector means to
! apply factors to all countries and/or sectors. 
! The number following the first text on line 1 (number 5 above) gives
! the number of pollutants treated in the file (ncols below). 
! Note: ncols can be greater than the number of emitted species we actually 
! have, and the species can be specified in any order, as given on the 
! top line. 
!-------------------------------------------------------------------------

  !/** local variables **/
  integer            :: ie, iq, ic, iland1, iland2 & ! loop variables
                       ,inland                     & ! Country read from femis
                       ,isec, isec1 , isec2        & ! loop vars: emis sectors
                       ,ncols, n, oldn               ! No. cols. in "femis" 
  integer, parameter        :: NCOLS_MAX = 20  ! Max. no. cols. in "femis"
  integer, dimension(NEMIS_FILE) :: qc        ! index for sorting femis columns
  real, dimension(NCOLS_MAX):: e_f             ! factors read from femis
  character(len=200) :: txt                    ! For read-in 
  character(len=20), dimension(NCOLS_MAX)::  polltxt ! to read line 1
 !--------------------------------------------------------


  e_fact(:,:,:) = 1.0            !/*** default value = 1 ***/

  associate ( debugm0 => ( DEBUG_GETEMIS .and. MasterProc ) )

  if( debugm0 ) print *, "Enters femis", me
  call open_file(IO_EMIS,"r","femis.dat",needed=.false.)

  if ( ios == NO_FILE ) then
        ios = 0
        write( *,*) "WARNING: NO FEMIS FILE"
        return !/** if no femis file, e_fact=1 as default **/ 
  endif
  call CheckStop( ios < 0 ,"EmisGet:ios error in femis.dat")


  ! Reads in the header line, e.g. name sec sox nox voc. 
  ! Pollutant names wil be checked against those defined in My_Emis_ml 

  read(unit=IO_EMIS,fmt="(a200)") txt
  if(debugm0)write(unit=6,fmt=*) "In femis, header0 is: ",  trim(txt)

  call wordsplit(txt,NCOLS_MAX,polltxt,ncols,ios)
  if(debugm0) then
     write(unit=6,fmt=*) "In femis, header is: ",  txt
     write(unit=6,fmt=*) "In femis, file has ", ncols, " columns (-2)"
  end if
  call CheckStop( ncols > NCOLS_MAX , "EmisGet:femisncols ncols > NCOLS_MAX" )
   if(ios>0)return

  !    we allow the femis file to give factors in any order, and
  !    for pollutants not needed, so we need to work out the indices
  !    for each column. Remember also that ncols includes the 1st
  !    2 columns (country_code and sector), which are not e_factors

  ncols = ncols - 2
  call CheckStop( ncols > NCOLS_MAX , "EmisGet:femisncols ncols > NCOLS_MAX" )
  call CheckStop( ncols < 1 , "EmisGet:femisncols ncols < 1" )

  n = 0
  COLS: do ic=1,ncols
      oldn = n
      EMLOOP: do ie=1, NEMIS_FILE
                if ( polltxt(ic+2) == trim ( EMIS_FILE(ie) ) ) then
                    qc(ie) = ic
                    n = n + 1
                    if(debugm0)write(unit=6,fmt=*) "In femis: ", &
                       polltxt(ic+2), " assigned to ", ie, EMIS_FILE(ie)
                  exit EMLOOP
                end if
      end do EMLOOP ! ie
       if (oldn == n .and.debugm0)   &
           write(unit=6,fmt=*) "femis: ",polltxt(ic+2)," NOT assigned"
  end do COLS   ! ic

if ( n < NEMIS_FILE ) then
  print *, "FEMIS me n NEMIS_FILE", me, n, NEMIS_FILE
  call CheckStop( n < NEMIS_FILE , "EmisGet: too few femis items" )
end if

  
  n = 0

  READFILE: do  ! ************ read lines of femis ***************

      read(unit=IO_EMIS,fmt=*,iostat=ios) inland, isec, (e_f(ic),ic=1,ncols)

      if ( ios <  0 ) exit READFILE                   ! End of file
      call CheckStop( ios > 0 , "EmisGet: read error in femis" )

      n = n + 1
      if(debugm0) then
        write(unit=6,fmt=*) "FEMIS READ", inland, &
          isec, (e_f(ic),ic=1,ncols)
        write(unit=6,fmt="(2a,I3,a,i3,a)") " Emission factors from femis.dat, ",&
          "landcode =", inland, ",  sector code =",isec, &
          " (sector 0 applies to all sectors) :"
        write(unit=6,fmt="(a,14(a,a,F5.2,a))") " ", (trim(polltxt(qc(ie)+2)),&
         " =",e_f(qc(ie)), ",  ", ie=1,NEMIS_FILE-1), &
          (trim(polltxt(qc(ie)+2))," =",e_f(qc(ie))," ", &
            ie=NEMIS_FILE,NEMIS_FILE)
      end if

      if (inland == 0 ) then     ! Apply factors to all countries
          iland1 = 1 
          iland2 = NLAND

      else                       ! Apply factors to country "inland"

! find country number  corresponding to index as written in emisfile
         do iland1=1,NLAND
            if(Country(iland1)%index==inland) goto 544
         enddo

         if(MasterProc) write(*,*)'COUNTRY CODE NOT RECOGNIZED',inland

         iland1 = 0
         iland2 =-1
544      continue
         if(iland1/=0)   iland2 = iland1
      end if

      if (isec == 0 ) then       ! All sectors
          isec1 = 1
          isec2 = NSECTORS
      elseif (isec==100) then    ! Anthropogenic scenario
          isec1 = 1
          isec2 = ANTROP_SECTORS
      else                       ! one sector: isec
          isec1 = isec
          isec2 = isec
      end if


      do ie = 1,NEMIS_FILE

          do iq = iland1, iland2
              do isec = isec1, isec2
                e_fact(isec,iq,ie) = e_fact(isec,iq,ie) * e_f( qc(ie) )
              end do !isec
          end do !iq

          if (debugm0 ) then
              write(unit=6,fmt=*) "IN NEMIS_FILE LOOP WE HAVE : ", ie, &
                                       qc(ie), e_f( qc(ie) )
              write(unit=6,fmt=*) "loops over ", isec1, isec2, iland1, iland2
          end if ! DEBUG_GETEMIS
      end do !ie
          
  enddo READFILE ! Loop over femis

  close(IO_EMIS)

  if( MasterProc) write(unit=6,fmt=*) "In femis, read ", n, "records from femis."
  if(debugm0) then
    ! Extra checks
     write(unit=6,fmt=*) "DEBUG_EMISGET: UK femis gives: "
     write(unit=6,fmt="(6x, 30a10)") (EMIS_FILE(ie), ie=1,NEMIS_FILE)
     do isec = 1, 11
      write(unit=6,fmt="(i6, 30f10.4)") isec, &
          (e_fact(isec,27,ie),ie=1,NEMIS_FILE)
     end do
  end if ! DEBUG
  ios = 0
  end associate ! debugm0
 end subroutine femis

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine EmisHeights()
   integer :: snap, k, allocerr
   real :: tmp(KMAX_MID)  ! values 
   character(len=200) :: txtinput               ! For read-in 
   character(len=20) :: txt1


   integer :: k_up
   real,allocatable:: emis_P_level(:)
   real :: P_emep,frac,sum
   real, parameter:: PT_EMEP=10000.0!Pa = 100 hPa
   integer :: isec,k_ext,k1_ext(KMAX_BND),nemis_hprofile

   !emis_hprofile are read from file. 
   !emis_kprofile are the fraction values converted into model levels

   !REMARK: if you only change nemis_hprofile, but keep exactely the same
   !        emissions, the results will still be sligthly changed.
   !        This is because nemis_hprofile is used to define KEMISTOP, which 
   !        defines the levels where to use 2 or 3 chemical "2steps" iterations.

   !use old format
   call open_file(IO_EMIS,"r","EmisHeights.txt",needed=.true.)


   do
      call read_line(IO_EMIS,txtinput,ios,'EmisHeight')

      if(me==1) write(*,fmt='(A)') "read from EmisHeights.txt : " // trim(txtinput)!, ios
      if ( ios <  0 ) exit     ! End of file
      if( index(txtinput,"#")>0 ) then ! Headers
         call PrintLog(trim(txtinput),MasterProc)
         cycle
      else if( index(txtinput,"Nklevels")>0 ) then !  Number levels
         read(txtinput,fmt=*,iostat=ios)  txt1, nemis_hprofile
         call PrintLog(trim(txtinput),MasterProc)
         allocate(emis_hprofile(nemis_hprofile+1,NSECTORS),stat=allocerr)
         allocate(emis_P_level(0:nemis_hprofile),stat=allocerr)
         emis_P_level=0.0
         call CheckStop(allocerr, "Allocation error for emis_P_level")
         emis_hprofile(:,:) = -999.9 
         emis_hprofile(1:nemis_hprofile+1,:) = 0.0
         cycle
      else if( index(txtinput,"Plevels")>0 ) then ! Pressure levels 
         read(txtinput,fmt=*,iostat=ios)  txt1, (emis_P_level(k),k=1, nemis_hprofile)
         call PrintLog(trim(txtinput),MasterProc)
         call CheckStop(allocerr, "Allocation error for emis_kprofile")
         emis_hprofile(:,:) = -999.9 
         emis_hprofile(1:nemis_hprofile+1,:) = 0.0
         cycle
      else
         read(txtinput,fmt=*,iostat=ios) snap, (tmp(k),k=1, nemis_hprofile)
         if( DEBUG_GETEMIS.and.MasterProc ) write(*,*) "VER=> ",snap, tmp(1), tmp(3)
         emis_hprofile(1:nemis_hprofile,snap) = tmp(1:nemis_hprofile)
      end if
   end do

   call CheckStop(nemis_hprofile < 1,"EmisGet: No EmisHeights set!!")
   call CheckStop( any( emis_hprofile(:,:) < 0 ), "EmisHeight read failure" )

   close(IO_EMIS)

   !Pressure boundaries for emission levels defined in EmisHeights.txt
   !NB in emis_P_level, k increase means higher up, i.e. smaller P (opposite as emep)
   emis_P_level(0)=Pref

   !can hardcode/override the values here. do not write more than nemis_kprofile (7?)
   !examples emep sigma levels:
   !  emis_P_level=0.0
   !  emis_P_level(1)=0.988 * (Pref-PT_EMEP)+PT_EMEP
   !  emis_P_level(2)=0.976 * (Pref-PT_EMEP)+PT_EMEP
   !  emis_P_level(3)=0.958 * (Pref-PT_EMEP)+PT_EMEP
   !  emis_P_level(4)=0.933 * (Pref-PT_EMEP)+PT_EMEP
   !  emis_P_level(5)=0.901 * (Pref-PT_EMEP)+PT_EMEP
   !  emis_P_level(6)=0.862 * (Pref-PT_EMEP)+PT_EMEP
   !  emis_P_level(7)=0.816 * (Pref-PT_EMEP)+PT_EMEP
   !!  emis_P_level(8)=0.763 * (Pref-PT_EMEP)+PT_EMEP
   !!  emis_P_level(9)=0.703 * (Pref-PT_EMEP)+PT_EMEP
   !!  emis_P_level(10)=0.636 * (Pref-PT_EMEP)+PT_EMEP

   if(emis_P_level(1)<1.0)then
      !the levels were not found in file. Assume lowest model levels
      if(me==1)write(*,*)'emission heights: assuming model levels'
      k_up=0
      do k=KMAX_BND-1,KMAX_BND-nemis_hprofile,-1
         k_up=k_up+1
         emis_P_level(KMAX_BND-k)=A_bnd(k)+B_bnd(k)*Pref !not used
      enddo
      nemis_kprofile=nemis_hprofile
      allocate(emis_kprofile(nemis_kprofile,NSECTORS),stat=allocerr)
      emis_kprofile(1:nemis_kprofile,:)=emis_hprofile(1:nemis_hprofile,:)

   else

      if(Masterproc)then
         write(*,*)'emission heights: defined from pressure levels'
         do k=0,nemis_hprofile
            write(*,*)'P emis levels : ',k,emis_P_level(k)
         enddo
      endif
      !stop
      !find highest level used
      nemis_kprofile = 0
      do k=KMAX_BND-1,1,-1
         nemis_kprofile = nemis_kprofile + 1
         if(A_bnd(k)+B_bnd(k)*Pref-0.0001<emis_P_level(nemis_hprofile))exit
      enddo
      if(MasterProc) write(*,*)'Emissions distributed among ',nemis_kprofile,' lowest levels'

      allocate(emis_kprofile(nemis_kprofile,NSECTORS),stat=allocerr)
      emis_kprofile=0.0

      !convert height (given as pressure) distribution into model level distribution
      k1_ext(KMAX_BND)=0
      do isec=1,NSECTORS! could put inside but easier for debugging to put here now
         if(DEBUG_GETEMIS.and.MasterProc) write(*,*)'Sector ',isec
         do k=KMAX_BND-1,max(1,KMAX_BND-nemis_kprofile),-1

!count all contributions to model level k
!i.e. between P_emep(k+1) and P_emep(k)

            P_emep=A_bnd(k)+B_bnd(k)*Pref !Pa
            if(DEBUG_GETEMIS.and.MasterProc) write(*,fmt="(A,I3,F10.2)")'vert_inter: P_emep',k,P_emep
            !largest available P_ext smaller than P_emep (if possible)
            !k1_ext(k) is the external layer just below P_emep(k)
            k1_ext(k)= 0 !start at surface, and go up until P_emep
            do k_ext=1,nemis_kprofile
               if(emis_P_level(k_ext)<P_emep)exit
               k1_ext(k)=k_ext
            enddo
            if(DEBUG_GETEMIS.and.MasterProc) write(*,*)k,k1_ext(k),k1_ext(k+1)

            !sum all contributions starting from last counted level (i.e. P_emep(k+1))

!part just above P_emep(k+1)
            if(emis_P_level(k1_ext(k+1)+1)>P_emep)then
               !part below k1_ext(k+1)+1  above P_emep(k+1)
               frac=((A_bnd(k+1)+B_bnd(k+1)*Pref )-emis_P_level(k1_ext(k+1)+1))/(emis_P_level(k1_ext(k+1))-emis_P_level(k1_ext(k+1)+1))
               emis_kprofile(KMAX_BND-k,isec)=frac*emis_hprofile(k1_ext(k+1)+1,isec)
               if(DEBUG_GETEMIS.and.MasterProc) write(*,fmt="(A,I5,6F10.2)")'adding fraction of level',&
                    k1_ext(k+1)+1,frac,emis_hprofile(k1_ext(k+1)+1,isec),emis_P_level(k1_ext(k+1)+1),&
                    (A_bnd(k+1)+B_bnd(k+1)*Pref ),emis_P_level(k1_ext(k+1)+1),(emis_P_level(k1_ext(k+1))-emis_P_level(k1_ext(k+1)+1))
            else
               !everything between P_emep(k+1) and P_emep(k)
               frac=((A_bnd(k+1)+B_bnd(k+1)*Pref )-P_emep)/(emis_P_level(k1_ext(k+1))-emis_P_level(k1_ext(k+1)+1))
               emis_kprofile(KMAX_BND-k,isec)=frac*emis_hprofile(k1_ext(k+1)+1,isec)
               if(DEBUG_GETEMIS.and.MasterProc) write(*,fmt="(A,I5,6F10.2)")'adding fraction of level between P_emep(k+1) and P_emep(k)',&
                    k1_ext(k+1)+1,frac,emis_hprofile(k1_ext(k+1)+1,isec),emis_P_level(k1_ext(k+1)+1),&
                    (A_bnd(k+1)+B_bnd(k+1)*Pref ),P_emep,(emis_P_level(k1_ext(k+1))-emis_P_level(k1_ext(k+1)+1))
            endif

            !add all full levels in between
            do k_ext=k1_ext(k+1)+2,k1_ext(k)
               emis_kprofile(KMAX_BND-k,isec)=emis_kprofile(KMAX_BND-k,isec)+emis_hprofile(k_ext,isec)
               if(DEBUG_GETEMIS.and.MasterProc) write(*,fmt="(A,I5,6F10.2)")'adding entire level',k_ext,emis_hprofile(k_ext,isec),emis_P_level(k_ext)
            enddo

            !add level just below P_emep(k), if not already counted, above k1_ext(k)  below P_emep
            if(emis_P_level(k1_ext(k+1)+1)>P_emep)then
               frac=(emis_P_level(k1_ext(k))-P_emep)/(emis_P_level(k1_ext(k))-emis_P_level(k1_ext(k)+1))
               emis_kprofile(KMAX_BND-k,isec)=emis_kprofile(KMAX_BND-k,isec)+frac*emis_hprofile(k1_ext(k)+1,isec)
               if(DEBUG_GETEMIS.and.MasterProc) write(*,fmt="(A,I5,6F10.2)")'adding last fraction of level',k1_ext(k)+1,&
                    frac,emis_hprofile(k1_ext(k)+1,isec),emis_P_level(k1_ext(k+1)+1),emis_P_level(k1_ext(k)),P_emep,&
                    (emis_P_level(k1_ext(k))-emis_P_level(k1_ext(k)+1))
            endif

         enddo
      enddo
      
      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)!so that output from different CPU does not get mixed up
      
      if(MasterProc)then
         write(*,*)'Distribution of emission into levels:'
         do isec=1,NSECTORS! could put inside but easier for debugging to put here now
            write(*,fmt="(A,I5,A,20F6.3)")'sector: ',isec,' fractions: ',(emis_kprofile(k,isec),k=1,nemis_kprofile)
         enddo
      endif

   endif

!check normalization of distributions:
   do isec=1,NSECTORS
      sum=0.0
      do k=1,nemis_kprofile
         sum=sum+emis_kprofile(k,isec)
      enddo
      if(abs(sum-1.0)>0.01)then
        if(MasterProc)then
           write(*,*)'WARNING emis height distribution not normalized : ',sum,(emis_kprofile(k,isec),k=1,nemis_kprofile)
           call StopAll( 'emis height distribution not normalized  ')
        endif
      endif
   enddo

 end subroutine EmisHeights
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine EmisSplit()
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!-------------------------------------------------------------------------   
!  DESCRIPTION:
!  Sets speciation for emissions to be splitted as defined in My_Emissions, 
!  e.g. VOC, PM, NOx, per source category for each country from input files. 
!  Done for species where NEMISFRAC > 1
!
! Input files:
!    xxxsplit.defaults   (e.g. vocsplit.defaults, pm25split.defaults)
!    xxxsplit.special
!        where xxx can be voc or pm or whatever has NEMISFRAC > 1
!
! Output: (module public variable)
!    real emisfrac(NRCSPLIT,NSECTORS,NLAND)
!
!------------------------------------------------------------------------- 

  !-- local
  integer ::  ie             ! emission index in EMIS_FILE (1..NEMIS_FILE)
  integer ::  itot           ! Index in IX_ arrays
  integer ::  iqrc           ! index of split compound in emisfrac   

  integer, parameter :: NONREACTIVE = 1   ! No. columns of non-reactive species
                                          ! enforced for read-ins.

  !-- for read-ins, dimension for max possible number of columns: 
  !-- for CRI we have 100s of VOC, hence
  character(len=10000) :: txtinput
  character(len=TXTLEN_NAME), dimension(0:1, NMAX ) :: intext
  character(len=TXTLEN_NAME), dimension(NMAX) :: Headers
! Now we expect a  "key-value", e.g. 46 for NOx as NO2, or
! zero to just use species()%molwt.:
  type(KeyVal), dimension(1)  :: MassValue ! set for e.f. NOx as NO2, SOx as SO2
  integer :: NKeys
  real, dimension(NSPEC_ADV,NSECTORS,NLAND) :: tmp_emisfrac
  real, dimension(NSPEC_ADV) :: tmp_emis_masscorr
  integer, dimension(NSPEC_ADV) :: tmp_iqrc2itot !maps from iqrc 
  real, dimension(NMAX ) :: tmp 
  real     :: sumtmp
  integer  :: nsplit   &         ! No.columns data to be read
             ,iland,isec,i,n,nn, allocerr
  integer  :: idef              ! Set to   0  for defaults, 1 for specials
  integer  :: iland1, iland2    ! loop variables over countries
  logical  :: defaults          ! Set to true for defaults, false for specials
  logical  :: debugm            ! debug flag on master proc 
!-----------------------------------------------

  iqrc = 0               ! Starting index in emisfrac array
  nrcsplit= 0                 !

  debugm = (DEBUG_GETEMIS.and.MasterProc) 

  do ie = 1, NEMIS_FILE 

    IDEF_LOOP: do idef = 0, 1

       defaults = (idef == 0)

       if ( defaults ) then

          fname = trim( "emissplit.defaults." // EMIS_FILE(ie) )
          call open_file(IO_EMIS,"r",fname,needed=.true.)

          call CheckStop( ios, "EmisGet: ioserror:split.defaults " )

       else 
    !** If specials exists, they will overwrite the defaults

          fname = trim( "emissplit.specials." // EMIS_FILE(ie) )
          call open_file(IO_EMIS,"r",fname,needed=.false.)

          if ( ios == NO_FILE ) then  
              ios = 0
              if(MasterProc) &
                 write(*,fmt=*) "emis_split: no specials for:",EMIS_FILE(ie)

              exit IDEF_LOOP
          endif
       end if


       if (debugm) write(*,*) "DEBUG_EMISGET split defaults=", defaults,fname
 
       !/ Read text line and speciation:
       !  the following lines expect one line of a header text with the
       !  species names, followed by lines of the following format:
       !  iland, isec, tmp1, tmp2.... tmpN+1, where the N+1'th column
       !  is  optional, and for non-reactive species. These non-reactives are
       !  not used in  the rest of the program, but are sometimes needed 
       !  (e.g. VOC) !  to check mass-balance.

        call Read_Headers(IO_EMIS,errmsg,nsplit,NKeys,Headers, MassValue)
        read(MassValue(1)%value,fmt=*) Emis_MolWt(ie)
        call CheckStop( errmsg , "Read Headers" // fname )
        call CheckStop( nsplit < 3 , "nsplit problem " // fname )

        nsplit = nsplit - 2

        if ( MasterProc ) then
          if(DEBUG_GETEMIS) then
             write(unit=6,fmt=*) "Will try to split ", nsplit , " times"
             write(unit=6,fmt=*) "Emis_MolWt  = ", Emis_MolWt(ie)
          end if
          write(unit=6,fmt=*) "Splitting ", trim(EMIS_FILE(ie)), &
             " emissions into ",&
               (trim(Headers(i+2)),' ',i=1,nsplit),'using ',trim(fname)
        end if

           do i = 1, nsplit
              intext(idef,i) = Headers(i+2)   ! 1st 2 columns are cc, isec:

             ! Match spec against EMIS_SPECS:

              call CountEmisSpecs( intext(idef,i) )

              if (debugm) write(*,*) "SPLITINFO iem ", i,idef, intext(idef,i)

              itot = find_index(intext(idef,i), species(:)%name )
              if ( defaults ) then
                if ( Headers(i+2) /= "UNREAC" ) then 
                  iqrc = iqrc + 1
                  emis_nsplit(ie) = emis_nsplit(ie) + 1
               
                  if ( MasterProc .and. itot<1 ) then 
                      print *, "EmisSplit FAILED idef ", me, idef, i, nsplit, trim( intext(idef,i) )
                      print *, " Failed Splitting ", trim(EMIS_FILE(ie)), &
                          " emissions into ",&
                            (trim(Headers(n+2)),' ',n=1,nsplit),'using ',trim(fname)
                      print "(a, i3,30a10)", "EmisSplit FAILED headers ", me, (intext(idef,n),n=1,nsplit)
                    call CheckStop( itot<1, &
                       "EmisSplit FAILED "//trim(intext(idef,i)) //&
                       " possible incorrect Chem in run script?" )
                  end if ! FAILURE

                  tmp_iqrc2itot(iqrc) = itot
                  itot2iqrc(itot)     = iqrc

                 ! Now, get factor needed for scaling emissions to molec
                  if ( Emis_MolWt(ie) == 0 ) then
                      tmp_emis_masscorr(iqrc) = 1.0/species(itot)%molwt
                  else
                      tmp_emis_masscorr(iqrc) = 1.0/Emis_MolWt(ie)
                  end if
                end if ! defaults
                if (debugm .and. itot>0 )  then
                   write(6,"(a,i2,i4,a,i4,a,a,f6.1)") &
                   "Mapping idef,iqrc:", idef, iqrc, "->", itot, &
                     trim(species(itot)%name ), " MW:", &
                       1.0/tmp_emis_masscorr(iqrc)
                end if
              end if
           end do
           if (debugm ) write(6,"(a,i4,a,i4)") "Compare ns: used=", &
                emis_nsplit(ie), "including any UNREAC:", nsplit
           
        n = 0

       READ_DATA: do 

           call read_line(IO_EMIS,txtinput,ios)
           if ( ios /=  0 ) exit READ_DATA     ! End of file
           read(unit=txtinput,fmt=*,iostat=ios)  iland, isec, (tmp(i),i=1, nsplit)
           if( MasterProc .and. ios /= 0 ) then
               print *, "ERROR: EmisGet: Failure reading emispslit file"
               print *, "Expecting to split into nsplit=", nsplit
               print *, "reading this record:", trim(txtinput)
              call CheckStop( ios ,"EmisGet: ios error on "//trim(fname) )
           end if

           n = n + 1
           if (debugm ) then
               write(6,"(a,i3,a,3i3,50f8.2)") "Splits: ",  n, trim(fname),&
                  iland, isec, nsplit, tmp(1:nsplit)
           end if

           !/... some checks:
           sumtmp = sum( tmp(1:nsplit) )
           if ( ( sumtmp  >    100.01 .or. sumtmp   <   99.99   )  .or.  &
                ( defaults .and. iland /= 0                     )  .or.  &
                ( defaults .and. isec  /= n                     )        &
                                                                )   then
               print * , "ERROR: emisfrac:" // trim(fname) // " "
               print *, "ERROR: emisfrac:", idef,  iland, n, isec, sumtmp
               write(unit=errmsg,fmt=*) &
            "ERROR: emisfrac:"//trim(fname)//" ", idef, iland, n, isec, sumtmp
               call CheckStop( errmsg )
           end if
           if ( .not. defaults ) then
             ! Check that specials headers match defaults
              if (debugm ) print *,"SPLIT CHECK SPECIES",idef
              do nn=1,emis_nsplit(ie)  
                  if ( MasterProc ) then
                    if (intext(1,nn) /= intext(0,nn) ) then
                      print *, "SPLIT ERROR DEF ", nn, "D:",trim(intext(0,nn))
                      print *, "SPLIT ERROR SPC ", nn, "S:",trim(intext(1,nn))
                    end if
                  end if ! masterproc
                    
                  call CheckStop( intext(1,nn) /= intext(0,nn), &
                    "EmisGet: ERROR intext(1,nn) /= intext(0,nn) ")
              end do !nn
           end if

           if ( defaults .or. iland == 0 ) then
             iland1 = 1
             iland2 = NLAND
           else  ! specials for one country
             iland1 = iland
             iland2 = iland
           end if
               
           do iland = iland1, iland2
             do i = 1, emis_nsplit(ie) !DSRC do i = 1, EMIS_NSPLIT(isp)

                !/** assign and convert from percent to fractions: **/

                iqrc = sum(emis_nsplit(1:ie-1)) + i
                tmp_emisfrac(iqrc,isec,iland) = 0.01 * tmp(i)

                ! just a check
                !if ( DEBUG .and. iland == 27.and.MasterProc ) then 
                if ( DEBUG_GETEMIS .and. iland == 101.and.MasterProc ) then 
                    itot = tmp_iqrc2itot(iqrc)
                    write(*,"(a35,4i3,i4,a,f10.4)") &
                      "DEBUG_EMISGET splitdef UK", isec, ie, i,  &
                       iqrc, itot, trim(species(itot)%name), &
                         tmp_emisfrac(iqrc,isec,iland)
                endif
             enddo ! i
           enddo ! iland

       enddo READ_DATA 
       close(IO_EMIS)

       call CheckStop(  defaults .and. n  /=  NSECTORS, &
                        "ERROR: EmisGet: defaults .and. n  /=  NSECTORS" )

       if (debugm ) write(*,*) "Read ", n, " records from ",fname

    end do IDEF_LOOP 
  end do ! ie
  ios = 0
  
  ! By now we should have found all the species required for the 
  ! chemical scheme:

  do ie = 1, NEMIS_SPECS
    if ( MasterProc .and. ( EmisSpecFound(ie) .eqv. .false.) ) then
       call PrintLog("WARNING: EmisSpec not found in snapemis. Ok if bio, nat, or fire!! " // trim(EMIS_SPECS(ie)) )
       write(*,*) "WARNING: EmisSpec - emissions of this compound were specified",&
&               " in the CM_reactions files, but not found in the ",&
&               " emissplit.defaults files. Make sure that the sets of files",&
&               " are consistent."
       ! Emissions can now be found in ForestFire module. No need to 
       ! stop
       ! call CheckStop ( any( EmisSpecFound .eqv. .false.), &
       !    "EmisSpecFound Error" )
    end if
  end do

  ! Now, we know how many split species we have, nrcsplit, so we allocate
  !  and fill emisfrac:
  nrcemis = sum( emis_nsplit(:) )
  allocate(emisfrac(nrcemis,NSECTORS,NLAND),stat=allocerr)
  call CheckStop(allocerr, "Allocation error for emisfrac")
  allocate(iqrc2itot(nrcemis),stat=allocerr)
  call CheckStop(allocerr, "Allocation error for iqrc2itot")
  allocate(emis_masscorr(nrcemis),stat=allocerr)
  call CheckStop(allocerr, "Allocation error for emis_masscorr")
  emisfrac(:,:,:)     = tmp_emisfrac(1:nrcemis,:,:)
  iqrc2itot(:)        = tmp_iqrc2itot(1:nrcemis)
  emis_masscorr(:)    = tmp_emis_masscorr(1:nrcemis)

!rb: not ideal place for this but used here for a start
! Temporary solution! Need to find the molweight from the
! GenChem input but just to get something running first set
! a hard coded molar mass of 200. 
  if(USE_ROADDUST)THEN
     allocate(roaddust_masscorr(NROADDUST),stat=allocerr)
     call CheckStop(allocerr, "Allocation error for emis_masscorr")
     if(MasterProc) &
          write(*,fmt=*)"NOTE! WARNING! Molar mass assumed to be 200.0 for all road dust components. Emissions will be in ERROR if another value is set in the GenChem input!"
     do ie=1,NROADDUST
        roaddust_masscorr(ie)=1.0/200.
     enddo
  endif
   
 end subroutine EmisSplit
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


 subroutine CountEmisSpecs( inspec )
  character(len=*), intent(in) :: inspec
  integer :: ind

  ind = find_index(inspec, EMIS_SPECS)
  call CheckStop ( ind > NEMIS_SPECS, "CountEmisSpecs Error" )
  if ( ind > 0 ) EmisSpecFound( ind ) = .true.

 end subroutine CountEmisSpecs

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine RoadDustGet(iemis,emisname,IRUNBEG,JRUNBEG,GIMAX,GJMAX, &
                     sumroaddust,        &
                     globroad_dust_pot,road_globnland,road_globland)

!.......................................................................
!  DESCRIPTION:
!  Reads in Road Dust emission potentials from one file, specified by iemis.
!  The arrays read in here are the global arrays (allocatable)
!.......................................................................

  !--arguments
  integer, intent(in)                     :: iemis     ! emis index
  character(len=*), intent(in)            :: emisname  ! emission name
  integer, intent(in) :: IRUNBEG,JRUNBEG,GIMAX,GJMAX   ! domain limits
  real,    intent(out), dimension(:,:,:)  :: globroad_dust_pot ! Road dust emission potentials
  integer, intent(inout), dimension(:,:,:)::   &
                                 road_globland   !Road emis.codes
  integer, intent(inout), dimension(:,:)  ::   &
                                 road_globnland  ! No. flat emitions in grid
  real,    intent(inout), dimension(:,:)  :: sumroaddust ! Emission potential sums per country

  !--local
  integer :: i, j, isec, iland,         &  ! loop variables
             iic,ic                        ! country code (read from file)
  real    :: tmpdust                       ! for reading road dust emission potential file
  integer, save :: ncmaxfound = 0          ! Max no. countries found in grid
  character(len=300) :: inputline

   !>============================

    if ( my_first_road ) then
         if(DEBUG_ROADDUST)WRITE(*,*)"initializing sumroaddust!"
         sumroaddust(:,:) =  0.0       ! initialize sums
         ios = 0
         my_first_road = .false.
    endif

  !>============================

      globroad_dust_pot(:,:,:) = 0.0

      if (DEBUG_ROADDUST) write(unit=6,fmt=*) "Called RoadDustGet with index, name", &
           iemis, trim(emisname)
!      fname = "emislist." // emisname
      fname = emisname
      call open_file(IO_EMIS,"r",fname,needed=.true.)
      call CheckStop(ios,"RoadDustGet: ios error1 in emission file")
 
      read(unit=IO_EMIS,fmt="(a200)",iostat=ios) inputline 
      if( inputline(1:1) .ne. "#" ) then ! Is a  comment
         write(*,*)'ERROR in road dust emission file!'
         write(*,*)'First line should be a comment line, starting with #'
      else
         write(*,*)'I read the comment line:',inputline
      endif     

READEMIS: do   ! ************* Loop over emislist files *******************

            read(unit=IO_EMIS,fmt=*,iostat=ios) iic,i,j, tmpdust

!            write(*,*)'dust to dust',iic,i,j, tmpdust

            if( DEBUG_ROADDUST .and. i==DEBUG%IJ(1) .and. j==DEBUG%IJ(2) ) write(*,*) &
                "DEBUG RoadDustGet "//trim(emisname) // ":" , iic, tmpdust
            if ( ios <  0 ) exit READEMIS            ! End of file
            call CheckStop(ios > 0,"RoadDustGet: ios error2 in emission file")

            ! Check if country code in emisfile (iic) is in the country list
            ! from Countries_ml, i.e. corresponds to numbering index ic

            do ic=1,NLAND
               if((Country(ic)%index==iic))&
                    goto 654
            enddo
            write(unit=errmsg,fmt=*) &
                   "COUNTRY CODE NOT RECOGNIZED OR UNDEFINED ", iic
            call CheckStop(errmsg)
            ic=0
654         continue

            i = i-IRUNBEG+1     ! for RESTRICTED domain
            j = j-JRUNBEG+1     ! for RESTRICTED domain

            if ( i  <=  0 .or. i  >  GIMAX .or.   &
                 j  <=  0 .or. j  >  GJMAX .or.   &
                 ic <=  0 .or. ic >  NLAND  )&
             cycle READEMIS

             ! ..........................................................
             ! generate new land allocation in 50 km grid. First, we check if
             ! country "ic" has already  been found within that grid. If not,
             ! then ic is added to landcode and nlandcode increased by one.

              call GridAllocate("ROAD"// trim ( emisname ),i,j,ic,NCMAX, &
                                 iland,ncmaxfound,road_globland,road_globnland)

              globroad_dust_pot(i,j,iland) = globroad_dust_pot(i,j,iland) &
                        + tmpdust
            if( DEBUG_ROADDUST .and. i==DEBUG%IJ(1) .and. j==DEBUG%IJ(2) ) write(*,*) &
                "DEBUG RoadDustGet iland, globrdp",iland,globroad_dust_pot(i,j,iland)

!   NOTE!!!! A climatological factor is still missing for the road dust!
!              should increase the emissions in dry areas by up to a factor of ca 3.3
!              Will be based on soil water content
!              (Fixed later....)

             ! Sum over all sectors, store as Ktonne:
!              if(tmpdust.lt.0.)write(*,*)'neg dust!', tmpdust
              sumroaddust(ic,iemis) = sumroaddust(ic,iemis)   &
                                  + 0.001 * tmpdust
!              if(sumroaddust(ic,iemis).lt.0.)write(*,*) &
!                 'We are on a road to nowhere', ic,iemis,tmpdust

        end do READEMIS
        
        if(DEBUG_ROADDUST) write(*,*)'done one, got sumroaddust=',sumroaddust(:,iemis)        

        close(IO_EMIS)
        ios = 0
  end subroutine RoadDustGet



end module EmisGet_ml
