! <EmisGet_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                     module EmisGet_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  use CheckStop_ml,      only: CheckStop
  use ChemSpecs_adv_ml,  only: NSPEC_ADV ! max possible number in split 
  use ChemSpecs_tot_ml,  only: NSPEC_TOT 
  use ChemChemicals_ml,  only: species
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
  use Io_ml,             only: open_file, NO_FILE, ios, IO_EMIS, &
                               Read_Headers, read_line, PrintLog
  use KeyValue_ml,       only: KeyVal
  use ModelConstants_ml, only: NPROC, TXTLEN_NAME, DEBUG => DEBUG_GETEMIS, &
                               DEBUG_i, DEBUG_j, &
                               KMAX_MID, &
              SEAFIX_GEA_NEEDED, & ! only if emission problems over sea
                               MasterProc,DEBUG_GETEMIS,DEBUG_ROADDUST,USE_ROADDUST
  use Par_ml,            only: me
  use SmallUtils_ml,     only: wordsplit, find_index
  use Volcanos_ml

  implicit none
  private

 ! subroutines:

  public  :: EmisGet           ! Collects emissions of each pollutant
  public  :: EmisSplit         ! => emisfrac, speciation of voc, pm25, etc.
  public  :: EmisHeights       ! => nemis_kprofile, emis_kprofile
                               !     vertical emissions profile
  public  :: RoadDustGet       ! Collects road dust emission potentials
  private :: femis             ! Sets emissions control factors 
  private :: CountEmisSpecs    !


  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
  logical, private, save :: my_first_call = .true.
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

  ! some common variables
  character(len=40), private :: fname             ! File name
  character(len=80), private :: errmsg

 ! Import list of the emitted species we need to find in the 
 ! emissplit files.
  include 'CM_EmisSpecs.inc'
  logical, dimension(NEMIS_SPECS) :: EmisSpecFound = .false.

 contains


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
 
    if ( my_first_call ) then
         sumemis(:,:) =  0.0       ! initialize sums
         ios = 0
         call femis()              ! emission factors (femis.dat file)
         if ( ios /= 0 )return
         my_first_call = .false.
    endif
  !>============================


      globemis   (:,:,:,:) = 0.0
      globemis_flat(:,:,:) = 0.0

      if (DEBUG) write(unit=6,fmt=*) "Called EmisGet with index, name", &
           iemis, trim(emisname)
      fname = "emislist." // emisname
      call open_file(IO_EMIS,"r",fname,needed=.true.)
      call CheckStop(ios,"EmisGet: ios error in emission file")

      if (trim ( emisname ) == "nh3" ) dknh3_agr=0.0 ! NH3Emis experimental
 

READEMIS: do   ! ************* Loop over emislist files *******************

            read(unit=IO_EMIS,fmt=*,iostat=ios) iic,i,j, duml,dumh,  &
                                    (tmpsec(isec),isec=1,NSECTORS)

            if( DEBUG .and. i==DEBUG_i .and. j==DEBUG_j ) write(*,*) &
                "DEBUG GetEmis "//trim(emisname) // ":" , iic, duml,dumh
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
                  if (DEBUG) write(*,*)'Volcano no. ',volc_no
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
  if(DEBUG_GETEMIS)write(unit=6,fmt=*) "In femis, header0 is: ",  trim(txt)

  call wordsplit(txt,NCOLS_MAX,polltxt,ncols,ios)
  if(DEBUG_GETEMIS) then
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
                    if(DEBUG_GETEMIS)write(unit=6,fmt=*) "In femis: ", &
                       polltxt(ic+2), " assigned to ", ie, EMIS_FILE(ie)
                  exit EMLOOP
                end if
      end do EMLOOP ! ie
       if (oldn == n .and.DEBUG_GETEMIS)   &
           write(unit=6,fmt=*) "femis: ",polltxt(ic+2)," NOT assigned"
  end do COLS   ! ic

  call CheckStop( n < NEMIS_FILE , "EmisGet: too few femis items" )

  
  n = 0

  READFILE: do  ! ************ read lines of femis ***************

      read(unit=IO_EMIS,fmt=*,iostat=ios) inland, isec, (e_f(ic),ic=1,ncols)

      if ( ios <  0 ) exit READFILE                   ! End of file
      call CheckStop( ios > 0 , "EmisGet: read error in femis" )

      n = n + 1
      if(DEBUG_GETEMIS)write(unit=6,fmt=*) "FEMIS READ", inland, &
          isec, (e_f(ic),ic=1,ncols)
      write(unit=6,fmt="(2a,I3,a,i3,a)") " Emission factors from femis.dat, ",&
        "landcode =", inland, ",  sector code =",isec, &
        " (sector 0 applies to all sectors) :"
      write(unit=6,fmt="(a,14(a,a,F5.2,a))") " ", (trim(polltxt(qc(ie)+2)),&
       " =",e_f(qc(ie)), ",  ", ie=1,NEMIS_FILE-1), &
        (trim(polltxt(qc(ie)+2))," =",e_f(qc(ie))," ", &
            ie=NEMIS_FILE,NEMIS_FILE)

      if (inland == 0 ) then     ! Apply factors to all countries
          iland1 = 1 
          iland2 = NLAND

      else                       ! Apply factors to country "inland"

! find country number  corresponding to index as written in emisfile
         do iland1=1,NLAND
            if(Country(iland1)%index==inland) goto 544
         enddo

         if(me==0) write(*,*)'COUNTRY CODE NOT RECOGNIZED',inland

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

          if (DEBUG ) then
              write(unit=6,fmt=*) "IN NEMIS_FILE LOOP WE HAVE : ", ie, &
                                       qc(ie), e_f( qc(ie) )
              write(unit=6,fmt=*) "loops over ", isec1, isec2, iland1, iland2
          end if ! DEBUG
      end do !ie
          
  enddo READFILE ! Loop over femis

  close(IO_EMIS)

  if(DEBUG_GETEMIS)write(unit=6,fmt=*) "In femis, read ", n, "records from femis."
  if ( DEBUG.and.MasterProc ) then    ! Extra checks
     write(unit=6,fmt=*) "DEBUG_EMISGET: UK femis gives: "
     write(unit=6,fmt="(6x, 30a10)") (EMIS_FILE(ie), ie=1,NEMIS_FILE)
     do isec = 1, 11
      write(unit=6,fmt="(i6, 30f10.4)") isec, &
          (e_fact(isec,27,ie),ie=1,NEMIS_FILE)
     end do
  end if ! DEBUG
  ios = 0
 end subroutine femis

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine EmisHeights()
    integer :: snap, k, allocerr
    real :: tmp(KMAX_MID)  ! values 
    character(len=200) :: txtinput               ! For read-in 
    character(len=20) :: txt1
     call open_file(IO_EMIS,"r","EmisHeights.txt",needed=.true.)
   
     do
        call read_line(IO_EMIS,txtinput,ios,'EmisHeight')
        if(me==1) print *, "EMIS HEIGHTS " // trim(txtinput)!, ios
        if ( ios <  0 ) exit     ! End of file
          if( index(txtinput,"#")>0 ) then ! Headers
            call PrintLog(trim(txtinput),MasterProc)
            cycle
          else if( index(txtinput,"Nklevels")>0 ) then !  Number levels
            read(txtinput,fmt=*,iostat=ios)  txt1, nemis_kprofile
            call PrintLog(trim(txtinput),MasterProc)
            allocate(emis_kprofile(nemis_kprofile,NSECTORS),stat=allocerr)
            call CheckStop(allocerr, "Allocation error for emis_kprofile")
            emis_kprofile(:,:) = -999.9 
            cycle
          else
            read(txtinput,fmt=*,iostat=ios) snap, (tmp(k),k=1, nemis_kprofile)
            if( DEBUG ) write(*,*) "VER=> ",snap, tmp(1), tmp(3)
            emis_kprofile(:,snap) = tmp(:)
          end if
      end do

     call CheckStop(nemis_kprofile < 1,"EmisGet: No EmisHeights set!!")
     call CheckStop( any( emis_kprofile(:,:) < 0 ), "EmisHeight read failure" )

     close(IO_EMIS)
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
!-----------------------------------------------

  iqrc = 0               ! Starting index in emisfrac array
  nrcsplit= 0                 !

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

       if (DEBUG.and.MasterProc) write(unit=6,fmt=*) &
             "DEBUG_EMISGET split defaults=", defaults, fname
 
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

              if(MasterProc.and.DEBUG_GETEMIS) write(*,*) "SPLITINFO iem ",&
                    i,idef, intext(idef,i)
              itot = find_index(intext(idef,i), species(:)%name )
              if ( defaults ) then
                if ( Headers(i+2) /= "UNREAC" ) then 
                  iqrc = iqrc + 1
                  emis_nsplit(ie) = emis_nsplit(ie) + 1
               
                  if ( me ==0 .and. itot<1 ) then 
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
                if ( MasterProc .and.DEBUG_GETEMIS.and. itot>0 ) &
                   write(6,"(a,i2,i4,a,i4,a,a,f6.1)") &
                   "Mapping idef,iqrc:", idef, iqrc, "->", itot, &
                     trim(species(itot)%name ), " MW:", &
                       1.0/tmp_emis_masscorr(iqrc)
                 end if
              !end if defaults
           end do
           if ( MasterProc.and.DEBUG_GETEMIS ) &
               write(6,"(a,i4,a,i4)") "Compare ns: used=", &
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
           if ( DEBUG .and. MasterProc ) then
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
              if ( DEBUG .and. MasterProc ) print *,"SPLIT CHECK SPECIES",idef
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
                if ( DEBUG .and. iland == 101.and.MasterProc ) then 
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
       if(MasterProc.and.DEBUG_GETEMIS) &
           write(unit=6,fmt=*) "Read ", n, " records from ",fname

    end do IDEF_LOOP 
  end do ! ie
  ios = 0
  
  ! By now we should have found all the species required for the 
  ! chemical scheme:

  do ie = 1, NEMIS_SPECS
    if ( MasterProc .and. ( EmisSpecFound(ie) .eqv. .false.) ) then
       call PrintLog("WARNING: EmisSpec not found in snapemis. Ok if in forest fire!! " // trim(EMIS_SPECS(ie)) )
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

            if( DEBUG_ROADDUST .and. i==DEBUG_i .and. j==DEBUG_j ) write(*,*) &
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
            if( DEBUG_ROADDUST .and. i==DEBUG_i .and. j==DEBUG_j ) write(*,*) &
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
