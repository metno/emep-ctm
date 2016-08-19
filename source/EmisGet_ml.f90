! <EmisGet_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                     module EmisGet_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  use CheckStop_ml,      only: CheckStop
  use Country_ml,        only: NLAND, IC_NAT, IC_VUL, Country
  use EmisDef_ml,        only: NSECTORS, ANTROP_SECTORS, NCMAX, FNCMAX, & 
                               ISNAP_SHIP, ISNAP_NAT
  use GridAllocate_ml,   only: GridAllocate
  use Io_ml,             only: open_file, NO_FILE, ios, IO_EMIS
  use ModelConstants_ml, only: NPROC
  use My_Emis_ml,        only: NEMIS, NRCSPLIT, EMIS_NAME, SPLIT_NAME,  &
                               NEMIS_PLAIN, NEMIS_SPLIT, EMIS_NSPLIT
  use Par_ml,            only: me
  use SmallUtils_ml,     only: wordsplit
  use Volcanos_ml

  implicit none
  private

 !/* subroutines:

  public  :: EmisGet           ! Collects emissions of each pollutant
  public  :: EmisSplit         ! => emisfrac, Speciation of voc, pm25, etc.
  private :: femis             ! Sets emissions control factors 


  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
  logical, private, save :: my_first_call = .true.

  logical, private, parameter :: DEBUG = .false.


  ! e_fact is the emission control factor (increase/decrease/switch-off)
  ! e_fact is read in from the femis file and applied within EmisGet
  real, private, save, &
         dimension(NSECTORS,NLAND,NEMIS)  :: e_fact 

  !/* emisfrac is used at each time-step of the model run to split
  !   emissions such as VOC; PM into species. 

  real, public, dimension(NRCSPLIT,NSECTORS,NLAND), save :: emisfrac

  !/ some common variables
  character(len=40), private :: fname             ! File name
  character(len=80), private :: errmsg

 contains


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine EmisGet(iemis,emisname,IRUNBEG,JRUNBEG,GIMAX,GJMAX, &
                     globemis,globnland,globland,sumemis,        &
                     globemis_flat,flat_globnland,flat_globland)

!.......................................................................
!**    DESCRIPTION:
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
  integer :: flat_iland, flat_nc,       &
             i, j, isec, iland, k, nc,  &  ! loop variables
             iic,ic                        ! country code (read from file)
  real    :: duml,dumh                     ! dummy variables, low/high emis.
  real, dimension(NSECTORS)  :: tmpsec     ! array for reading emission files
  integer, save :: ncmaxfound = 0          ! Max no. countries found in grid
  integer, save :: flat_ncmaxfound = 0     ! Max no. countries found in grid
                                           ! including flat emissions

   !>============================
 
    if ( my_first_call ) then
         sumemis(:,:) =  0.0       ! Initialise sums
         ios = 0
         call femis()              ! emission factors (femis.dat file).
         if ( ios /= 0 )return
         my_first_call = .false.
    endif
  !>============================


      globemis   (:,:,:,:) = 0.0
      globemis_flat(:,:,:) = 0.0

      if (DEBUG) write(unit=6,fmt=*) "Called EmisGet with index, name", iemis, emisname
      fname = "emislist." // emisname
      call open_file(IO_EMIS,"r",fname,needed=.true.)
      call CheckStop(ios,"EmisGet: ios error in emission file")

READEMIS: do   ! ************* Loop over emislist files *******************

            read(unit=IO_EMIS,fmt=*,iostat=ios) iic,i,j, duml,dumh,  &
                                    (tmpsec(isec),isec=1,NSECTORS)

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
                 ic == IC_NAT )                   &  ! Excludes DMS
             cycle READEMIS

             !/* Ship emissions

             if ( Country(ic)%is_sea ) then    ! ship emissions

              ! ..........................................................
              ! Generate new land allocation in 50 km grid for FLAT 
              ! EMISSIONS (ships). First, we check if country "ic"
              ! has already  been found within that grid. If not, then ic is
              ! added to flat_landcode and flat_nlandcode increased by one.
   
              !/** Test that ship emissions are only in sector ISNAP_SHIP
               do isec=1,(ISNAP_SHIP-1) 
                  call CheckStop(tmpsec(isec) /= 0,  &
                        "EmisGet: NOT FLAT EMISSIONS")
               enddo
               do isec=ISNAP_SHIP+1,NSECTORS
                  call CheckStop(tmpsec(isec) /= 0,  &
                        "EmisGet: NOT FLAT EMISSIONS")
               enddo
              !/** end test

               call GridAllocate("FLat",i,j,ic,FNCMAX, flat_iland, &
                   flat_ncmaxfound,flat_globland,flat_globnland)
              ! ...................................................
              ! Assign e_fact corrected emissions to global FLAT 
              ! emission matrices.
              ! ...................................................

               globemis_flat(i,j,flat_iland) = globemis_flat(i,j,flat_iland) &
                     + e_fact(ISNAP_SHIP,ic,iemis) * tmpsec(ISNAP_SHIP)

              !......................................................
              !..        Sum over all sectors, store as Ktonne:
              !......................................................

                 sumemis(ic,iemis) = sumemis(ic,iemis)   &
                      + 0.001 * globemis_flat(i,j,flat_iland)  

                cycle READEMIS
             endif !ship emissions            


              !.......................................................
              !/**  Volcanos
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
             !/** end Volcanoes
             !..............................................................


             !  For VOC natural and agricultur emissions (managed forests) 
             !  set to  zero

             if ( trim ( emisname ) == "voc" ) tmpsec(11:11) = 0.0

             ! ..........................................................
             ! generate new land allocation in 50 km grid. First, we check if
             ! country "ic" has already  been found within that grid. If not,
             ! then ic is added to landcode and nlandcode increased by one.

              call GridAllocate("SNAP"// trim ( emisname ),i,j,ic,NCMAX, &
                                 iland,ncmaxfound,globland,globnland)

             ! ...................................................
             ! ...................................................

              globemis(:,i,j,iland) = globemis(:,i,j,iland) &
                        + e_fact(:,ic,iemis) *  tmpsec(:)


             !..        Sum over all sectors, store as Ktonne:

              sumemis(ic,iemis) = sumemis(ic,iemis)   &
                                  + 0.001 * sum (globemis (:,i,j,iland))
              
        end do READEMIS 
        !
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
  integer, parameter        :: NCOLS_MAX = 12      ! Max. no. cols. in "femis"
  integer, dimension(NEMIS) :: qc           ! index for sorting femis columns
  real, dimension(NCOLS_MAX):: e_f                 ! factors read from femis
  character(len=80) :: txt                         ! For read-in 
  character(len=5), dimension(NCOLS_MAX)::  polltxt! to read line 1
 !--------------------------------------------------------


  e_fact(:,:,:) = 1.0            !/*** default value = 1 ***/


  call open_file(IO_EMIS,"r","femis.dat",needed=.false.)

  if ( ios == NO_FILE ) then
        ios = 0
        return !/** if no femis file, e_fact=1 as default **/ 
  endif
  call CheckStop( ios < 0 ,"EmisGet:ios error in femis.dat")


  !/** Reads in the header line, e.g. name sec sox nox voc. 
  !    Pollutant names wil be checked against those defined in My_Emis_ml **/

  read(unit=IO_EMIS,fmt="(a80)") txt

  call wordsplit(txt,NCOLS_MAX,polltxt,ncols,ios)
   if(ios>0)return
  write(unit=6,fmt=*) "In femis, header is: ",  txt
  write(unit=6,fmt=*) "In femis, file has ", ncols, " columns (-2)"

  !/** we allow the femis file to give factors in any order, and
  !    for pollutants not needed, so we need to work out the indices
  !    for each column. Remember also that ncols includes the 1st
  !    2 columns (country_code and sector), which are not e_factors

  ncols = ncols - 2
  call CheckStop( ncols > NCOLS_MAX , "EmisGet:femisncols ncols > NCOLS_MAX" )
  call CheckStop( ncols < 1 , "EmisGet:femisncols ncols < 1" )

  n = 0
  COLS: do ic=1,ncols
      oldn = n
      EMLOOP: do ie=1, NEMIS
                if ( polltxt(ic+2) == trim ( EMIS_NAME(ie) ) ) then
                    qc(ie) = ic
                    n = n + 1
                    write(unit=6,fmt=*) "In femis: ", polltxt(ic+2),  &
                              " assigned to ", ie, EMIS_NAME(ie)
                  exit EMLOOP
                end if
      end do EMLOOP ! ie
       if (oldn == n)   &
           write(unit=6,fmt=*) "femis: ",polltxt(ic+2)," NOT assigned"
  end do COLS   ! ic

  call CheckStop( n < NEMIS , "EmisGet: too few femis items" )

  
  n = 0

  READFILE: do  ! ************ read lines of femis ***************

      read(unit=IO_EMIS,fmt=*,iostat=ios) inland, isec, (e_f(ic),ic=1,ncols)

      if ( ios <  0 ) exit READFILE                   ! End of file
      call CheckStop( ios > 0 , "EmisGet: read error in femis" )

      n = n + 1
      write(unit=6,fmt=*) "FEMIS READ", inland, isec, (e_f(ic),ic=1,ncols)

      if (inland == 0 ) then     ! Apply factors to all countries
          iland1 = 1 
          iland2 = NLAND
      else                       !  Apply factors to country "inland"

!.. find country number  corresponding to index as written in emisfile
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


      do ie = 1,NEMIS

          do iq = iland1, iland2
              do isec = isec1, isec2
                e_fact(isec,iq,ie) = e_fact(isec,iq,ie) * e_f( qc(ie) )
              end do !isec
          end do !iq

          if (DEBUG ) then
              write(unit=6,fmt=*) "IN NEMIS LOOP WE HAVE : ", ie, &
                                       qc(ie), e_f( qc(ie) )
              write(unit=6,fmt=*) "loops over ", isec1, isec2, iland1, iland2
          end if ! DEBUG
      end do !ie
          
  enddo READFILE ! Loop over femis

  close(IO_EMIS)

  write(unit=6,fmt=*) "In femis, read ", n, "records from femis."
  if ( DEBUG ) then    ! Extra checks
     write(unit=6,fmt=*) " For UK this gives: "
     write(unit=6,fmt="(6x, 10a12)") (EMIS_NAME(ie), ie=1,NEMIS)
     do isec = 1, 11
      write(unit=6,fmt="(i6, 10f12.6)") isec,(e_fact(isec,27,ie),ie=1,NEMIS)
     end do
  end if ! DEBUG
  ios = 0
 end subroutine femis

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine EmisSplit()
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!-------------------------------------------------------------------------   
!** DESCRIPTION:
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
  integer ::  ie                  ! emission index in EMIS_NAME   (1..NEMIS)
  integer ::  isp                 ! emission index in 1..NEMIS_SPLIT
  integer ::  ifr0, ifr           ! index of split compound in emisfrac   

  integer, parameter :: NONREACTIVE = 1   ! No.columns of non-reactive species
                                          ! enforced for read-ins.

  !-- for read-ins, dimension for max possible number of columns: 
  character(len=12), dimension(0:1, NRCSPLIT + NONREACTIVE ) :: intext
  real             , dimension     (NRCSPLIT + NONREACTIVE ) :: tmp 
  real     :: sumtmp
  integer  :: nsplit  &         ! No.columns data to be read
             ,iland,isec,i,n,nn
  integer  :: idef              ! Set to   0  for defaults, 1 for specials
  integer  :: iland1, iland2    ! loop variables over countries
  logical  :: defaults          ! Set to true for defaults, false for specials
!-----------------------------------------------

  if ( NEMIS_SPLIT == 0 ) return  !/** for safety **/


  ifr0 =    1                 ! Starting index in emisfrac array

  do isp = 1, NEMIS_SPLIT
    ie =   NEMIS_PLAIN + isp   ! Split species, index in 1..NEMIS

    nsplit =    EMIS_NSPLIT(isp) + NONREACTIVE 

    if (isp > 1) ifr0 = ifr0 + EMIS_NSPLIT(isp-1) !start index of next species

    !/ Just in case ....
    call CheckStop( EMIS_NAME(ie) /= SPLIT_NAME(isp) ,   &
                   "EmisGet: Mis-matchSPLIT" )

    IDEF_LOOP: do idef = 0, 1

       defaults = (idef == 0)

    !** Check if *.split.defaults file for the component exists
       if ( defaults ) then

          fname = trim( EMIS_NAME(ie) ) // "split.defaults"
          call open_file(IO_EMIS,"r",fname,needed=.true.)

          call CheckStop( ios, "EmisGet: ioserror:split.defaults " )

       else 
    !** If specials exists, they will overwrite the defaults

          fname = trim( EMIS_NAME(ie) ) // "split.special"
          call open_file(IO_EMIS,"r",fname,needed=.false.)

          if ( ios == NO_FILE ) then  
            write(unit=6,fmt=*) "emis_split: no specials for",EMIS_NAME(ie)
              ios = 0
              exit IDEF_LOOP
          endif
       end if

       if (DEBUG) write(unit=6,fmt=*) "TTT split defaults=", defaults, fname
 
       !/ Read text line and speciation:
       !  the following lines expect one line of a header text with the
       !  species names, followed by lines of the following format:
       !  iland, isec, tmp1, tmp2.... tmpN+1, where the N+1'th column
       !  is for non-reactive species. These non-reactives are not used in  
       !  the rest of the program, but are required to check mass-balance.

       read (unit=IO_EMIS,fmt=*,iostat=ios)   &
                            iland, isec ,(intext(idef,i), i=1, nsplit)
 
       call CheckStop( ios , "EmisGet: Read error on header, emis_split " )

       write(unit=6,fmt="(a25,i3,/,(12a7))") "SPLIT species for idef=", &
                      idef, (intext(idef,i), i=1, nsplit)
       write(unit=6,fmt=*) "Will try to split ", EMIS_NSPLIT(isp) , " times"
           
       n = 0

       READ_DATA: do 
          read (unit=IO_EMIS,fmt=*,iostat=ios)   &
                             iland, isec, (tmp(i),i=1, nsplit)

           if ( ios <  0 ) exit READ_DATA     ! End of file
           call CheckStop( ios > 0 , "EmisGet: Readerror on emis_split " )

           n = n + 1

           !/... some checks:
           sumtmp = sum( tmp(1:nsplit) )
           if ( ( sumtmp  >    100.01 .or. sumtmp   <   99.99   )  .or.  &
                ( defaults .and. iland /= 0                     )  .or.  &
                ( defaults .and. isec  /= n                     )        &
                                                                )   then
               write(unit=errmsg,fmt=*) "ERROR: emisfrac:", idef, &
                                    iland, isec, sumtmp
               call CheckStop( errmsg )
           end if
           if (  .not. defaults ) then
              do nn=1,nsplit
                  call CheckStop( intext(1,nn) /= intext(0,nn), &
                             "EmisGet: ERROR intext(1,nn) /= intext(0,nn) ")
              enddo
           end if

           if ( defaults .or. iland == 0 ) then
             iland1 = 1
             iland2 = NLAND
           else  ! specials for one country
             iland1 = iland
             iland2 = iland
           end if
               
           do iland = iland1, iland2
             do i = 1, EMIS_NSPLIT(isp)
                ifr = ifr0 + i - 1    ! => index in emisfrac array

                !/** assign and convert from percent to fractions: **/

                emisfrac(ifr,isec,iland) = 0.01 * tmp(i)

                ! just a check
                if ( DEBUG .and. iland == 27 ) then 
                    write(*,"(a15,3i3,f10.4)") "TTT splitdef UK", isec, &
                                   ifr0, ifr, emisfrac(ifr,isec,iland)
                endif
             enddo ! i
           enddo ! iland

       enddo READ_DATA 
       close(IO_EMIS)

       call CheckStop(  defaults .and. n  /=  NSECTORS, &
                        "ERROR: EmisGet: defaults .and. n  /=  NSECTORS" )
       write(unit=6,fmt=*) "Read ", n, " records from ",fname

    end do IDEF_LOOP 
  end do ! ie
  ios = 0

 end subroutine EmisSplit
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module EmisGet_ml
