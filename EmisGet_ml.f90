! <EmisGet_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_10(3282)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2016 met.no
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
module EmisGet_ml

use CheckStop_ml,     only: CheckStop, StopAll, check=>CheckNC
use ChemSpecs,        only: NSPEC_ADV, NSPEC_TOT, species 
use Country_ml,       only: NLAND, IC_NAT, IC_VUL, IC_NOA, Country, &
                             ! NMR-NH3 specific variables (hb NH3Emis)
                             IC_NMR 
use EmisDef_ml,       only: NSECTORS, ANTROP_SECTORS, NCMAX, FNCMAX, & 
                            N_HFAC,N_SPLIT, NEMIS_FILE, EMIS_FILE, & 
                            VOLCANOES_LL, &
                          ! NMR-NH3 specific variables (for FUTURE )
                            NH3EMIS_VAR,dknh3_agr,ISNAP_AGR,ISNAP_TRAF, &
                            NROADDUST, &
                            cdfemis,sumcdfemis,nGridEmisCodes,GridEmisCodes,GridEmis&
                            ,gridrcemis,gridrcemis0&
                            ,landcode,nlandcode,MAXFEMISLONLAT,N_femis_lonlat    
  
use GridAllocate_ml,  only: GridAllocate
use GridValues_ml,    only: debug_proc,debug_li,debug_lj,i_fdom,j_fdom,i_local,j_local
use GridValues_ml,    only: glon, glat, A_bnd, B_bnd
use Io_ml,            only: open_file,IO_LOG, NO_FILE, ios, IO_EMIS, &
                             Read_Headers, read_line, PrintLog
use KeyValueTypes,    only: KeyVal
use ModelConstants_ml,only: NPROC, TXTLEN_NAME, &
                             DEBUG,  KMAX_MID,KMAX_BND, Pref,&
                             SEAFIX_GEA_NEEDED, & ! only if emission problems over sea
                             MasterProc,DEBUG_GETEMIS,DEBUG_ROADDUST,USE_ROADDUST&
                             ,IIFULLDOM,JJFULLDOM, SECTORS_NAME
use MPI_Groups_ml   , only : MPI_BYTE, MPI_REAL8, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_INTEGER&
                                     ,MPI_COMM_CALC, IERROR
use NetCDF_ml, only  : ReadField_CDF  !CDF_SNAP

use Par_ml,            only: LIMAX, LJMAX, limax, ljmax, me
use SmallUtils_ml,     only: wordsplit, find_index, key2str
use netcdf,            only: NF90_OPEN,NF90_NOERR,NF90_NOWRITE,&
                             NF90_INQUIRE,NF90_INQUIRE_VARIABLE,NF90_CLOSE

implicit none
private

 ! subroutines:
public  :: EmisGetCdf        ! cdfemis
public  :: EmisGetCdfFrac    ! cdfemis in Fractions format
public  :: EmisGetASCII      !new version of "EmisGet" which does not use fulldomain arrays
public  :: EmisSplit         ! => emisfrac, speciation of voc, pm25, etc.
public  :: EmisHeights       ! => nemis_kprofile, emis_kprofile
                             !     vertical emissions profile
public  :: RoadDustGet       ! Collects road dust emission potentials
public  :: femis             ! Sets emissions control factors 
private :: CountEmisSpecs    !


!logical, private, save :: my_first_call = .true.
logical, private, save :: my_first_road = .true.

! e_fact is the emission control factor (increase/decrease/switch-off)
! e_fact is read in from the femis file and applied within EmisGet
real, public, save, allocatable, &
!         dimension(NSECTORS,NLAND,NEMIS_FILE)  :: e_fact 
   dimension(:,:,:)  :: e_fact 
! e_fact_lonlat is read in from the femis file and applied within EmisGet using lat lon format
real, public, save, allocatable, &
!         dimension(NSECTORS,MAXFEMISLONLAT,NEMIS_FILE)  :: e_fact_lonlat
   dimension(:,:,:)  :: e_fact_lonlat
real, public, save, dimension(MAXFEMISLONLAT) :: femis_latmin,femis_latmax,femis_lonmin,femis_lonmax

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
integer, private ::i_femis_lonlat

contains

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine EmisGetCdfFrac(iem, isec, fname, varname, sumemis_local, &
    incl, nin, excl, nex)

    implicit none
    integer, intent(in) ::iem, isec, nin, nex
    character(len=*),intent(in) :: fname, varname, incl(*),excl(*)
   ! Sum of emissions per country
    real,intent(inout), dimension(NLAND,NEMIS_FILE) :: sumemis_local 


    real :: fractions(LIMAX,LJMAX,NCMAX),Reduc(NLAND)
    character(len=125) ::Mask_fileName,Mask_varname
    real :: Mask_ReducFactor,lonlat_fac
    integer :: NMask_Code,Mask_Code(NLAND)

    integer ::i,j,k,n,ic,i_gridemis,found
    logical :: Cexist

    !yearly grid independent netcdf fraction format emissions                                
    Reduc=e_fact(isec,:,iem)
    call ReadField_CDF(trim(fname),varname,cdfemis(1,1),nstart=1,&
         interpol='mass_conservative',fractions_out=fractions,&
         CC_out=landcode,Ncc_out=nlandcode,Reduc=Reduc,needed=.true.,debug_flag=.false.,&
         Undef=0.0)

    do j=1,ljmax
       do i=1,limax                   
          if(nlandcode(i,j)>NCMAX)then
             write(*,*)"GetCdfFrac: Too many emitter countries in one gridcell: ",&
                  me,i,j,nlandcode(i,j)
             call StopAll("To many countries in one gridcell ")
          endif
          lonlat_fac=1.0
          if(N_femis_lonlat>0)then
             do i_femis_lonlat=1,N_femis_lonlat
                if(glat(i,j)>femis_latmin(i_femis_lonlat).and.&
                     glat(i,j)<femis_latmax(i_femis_lonlat).and.&
                     glon(i,j)>femis_lonmin(i_femis_lonlat).and.&
                     glon(i,j)<femis_lonmax(i_femis_lonlat)     )then
                   lonlat_fac=lonlat_fac*e_fact_lonlat(isec,i_femis_lonlat,iem) 
                endif
             enddo
          endif
          do n=1,nlandcode(i,j)
             ic=find_index(landcode(i,j,n),Country(:)%icode)

             !exclude or include countries
             !could easily be optimized, by defining a country mask
             if(nin>0)then
                !1) check that country is in include list
                found=find_index(Country(ic)%code ,incl(1:nin),first_only=.true.)
                if(found<=0)cycle!do not include
             endif
             if(nex>0)then
                !1) check that country is not in exclude list
                found=find_index(Country(ic)%code ,excl(1:nex),first_only=.true.)
                if(found>0)cycle!exclude
             endif

             if(Country(ic)%icode/=landcode(i,j,n))then
                write(*,*)"COUNTRY CODE ERROR: ",landcode(i,j,n),ic,Country(ic)%icode
                call StopAll("COUNTRY CODE ERROR ")
             endif
             if(ic>NLAND)then
                write(*,*)"COUNTRY CODE NOT RECOGNIZED OR UNDEFINED: ",landcode(i,j,n)
                call StopAll("COUNTRY CODE NOT RECOGNIZED ")
             endif

             !merge into existing emissions           
             Cexist=.false.
             do i_gridemis=1,nGridEmisCodes(i,j)
                if(GridEmisCodes(i,j,i_gridemis)==landcode(i,j,n))then

                   GridEmis(isec,i,j,i_gridemis,iem)=GridEmis(isec,i,j,i_gridemis,iem)&
                        +fractions(i,j,n)*cdfemis(i,j) *lonlat_fac                     

                   Cexist=.true.
                   exit
                endif
             enddo
             if(.not.Cexist)then
                !country not included yet. define it now:
                nGridEmisCodes(i,j)=nGridEmisCodes(i,j)+1
                if(nGridEmisCodes(i,j)>NCMAX)then
                   write(*,*)"Too many emitter countries in one gridemiscell: ",&
                        me,i,j,nGridEmisCodes(i,j)
                   call StopAll("To many countries in one gridemiscell ")
                endif
                i_gridemis=nGridEmisCodes(i,j)
                GridEmisCodes(i,j,i_gridemis)=landcode(i,j,n)
                GridEmis(isec,i,j,i_gridemis,iem)=fractions(i,j,n)*cdfemis(i,j)*lonlat_fac  
             endif
             sumemis_local(ic,iem)=sumemis_local(ic,iem)&
                  +0.001*fractions(i,j,n)*cdfemis(i,j)*lonlat_fac  !for diagnostics, mass balance
          enddo
       enddo
    enddo

  end subroutine EmisGetCdfFrac
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine EmisGetCdf(iem, fname, sumemis, incl,excl)
   integer, intent(in) :: iem ! index in EMIS_FILE array and GridEmis output
   character(len=*), intent(in)    :: fname
   real, intent(inout) ::sumemis(*)
   character(len=*),dimension(:), optional :: &
       incl, excl ! Arrays of cc to inc/exclude
   integer :: i,j, ic, isec, allocerr(6), icode, status
   real, dimension(NLAND) :: sumcdfemis_loc, sumcdfemis_iem
   integer :: icc, ncc
   character(len=40) :: varname, fmt
   integer, save :: ncmaxfound = 0 ! Max no. countries found in grid
   integer, save :: ncalls=0
   integer :: ncFileID, nDimensions,nVariables,nAttributes,timeDimID,varid,&
           xtype,ndims  !TESTE testing
   character(len=10) :: ewords(7), code ! Test Emis:UK:snap:7
   character(len=*), parameter ::  sub = 'EmisGetCdf:'
   integer :: nwords, err, i_gridemis
   logical, save :: my_first_call = .true., dbg0, dbgp
   logical :: Cexist
   real ::lonlat_fac
!
   if( my_first_call  )  then
        allocerr = 0 ! ?NEEDED?
        dbgp =  DEBUG_GETEMIS .and. debug_proc
        dbg0 = (MasterProc.and.DEBUG_GETEMIS)
        my_first_call = .false.
   end if

   if(MasterProc) write(*,*) sub//"START ", trim(fname)

   if( present(incl) ) then
      fmt="(a,i4,99a4)"
      if(dbg0) write(*,fmt) sub//trim(fname)//":INCL=", size(incl), incl
   end if
   if( present(excl) ) then
      fmt="(a,i4,99a4)"
      if(dbg0) write(*,fmt) sub//trim(fname)//":EXCL=", size(excl), excl
   end if


    ! Initialise sums for this pollutant and processor
    sumcdfemis_loc = 0.0

!---------------------------------------------------------------
! find emis file and  main properties

 status=nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncFileID)

 if( status /= nf90_noerr ) then
  if( MasterProc )write(*,*) "ERROR: EmisGetCdf - couldn't open "//trim(fName)
  CALL MPI_FINALIZE(IERROR)
  stop
 end if

 call check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,timeDimID))
 if( MasterProc .or. index(fname, "Ship")>0 ) then
  write( *,*) 'Nb of global attributes: ',nAttributes
  write( *,*) sub//trim(fName),' properties: '
  write( *,*) sub//' Num dimensions: ',nDimensions
  write( *,*) sub//' Num variables: ',nVariables
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
        if ( find_index( code, excl ) >0 ) then
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

     ic = find_index( code, Country(:)%code )  !from Country_ml

     if ( Country(ic)%code == "N/A" ) then
          if(MasterProc) write(*,*) "CDFCYCLE ", ic, trim(Country(ic)%code)
          cycle    ! see Country_ml
     end if
     call CheckStop( ic < 1 , "CDFEMIS NegIC:"//trim(fname)//":"//trim(code) )

     if( dbgp ) write(*,"(2a,i6,a,2i4)") 'EmisGetCdf ', &
           trim(fname), varid," "//trim(varname)// " "// trim(code), ic, isec 


     cdfemis = 0.0 ! safety, shouldn't be needed though

     call ReadField_CDF(fname,varname,cdfemis,1,&
               interpol='mass_conservative',&
                needed=.false.,UnDef=0.0,&
                debug_flag=.false.,ncFileID_given=ncFileID)

     if( maxval(cdfemis ) < 1.0e-10 ) cycle ! Likely no emiss in domain

     if ( dbgp ) then !me==0 has some sea usually !!SHIPHUNT
         ncalls = ncalls + 1
         write(*,"(3a,3i4,i3,9f12.4)")"CDF emis-in ",&
           trim(fname)//":", &
           trim(Country(ic)%name)//":"//trim(EMIS_FILE(iem)), ic, &
            i_fdom(debug_li),j_fdom(debug_lj),&
             isec, e_fact(isec,ic,iem), &
                cdfemis(debug_li, debug_lj), maxval(cdfemis)
      end if



        do j = 1, ljmax
            do i = 1, limax

              if( cdfemis(i,j) > 1.0e-10 ) then
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
                    ic = IC_NOA ! Change to North Africa from IIASA system
                 end if !========== AFRICA

                 lonlat_fac=1.0
                 if(N_femis_lonlat>0)then
                    do i_femis_lonlat=1,N_femis_lonlat
                       if(glat(i,j)>femis_latmin(i_femis_lonlat).and.&
                            glat(i,j)<femis_latmax(i_femis_lonlat).and.&
                            glon(i,j)>femis_lonmin(i_femis_lonlat).and.&
                            glon(i,j)<femis_lonmax(i_femis_lonlat)     )then
                          lonlat_fac=lonlat_fac*e_fact_lonlat(isec,i_femis_lonlat,iem) 
                       endif
                    enddo
                 endif

               sumcdfemis_loc(ic) = sumcdfemis_loc(ic) + &
                    0.001 *  e_fact(isec,ic,iem) * cdfemis(i,j) *lonlat_fac 

                 !merge into existing emissions           
                 Cexist=.false.
                 do i_gridemis=1,nGridEmisCodes(i,j)
                    if(GridEmisCodes(i,j,i_gridemis)==Country(ic)%icode)then
                       
                       GridEmis(isec,i,j,i_gridemis,iem)= &
                             GridEmis(isec,i,j,i_gridemis,iem)&
                            +cdfemis(i,j) * e_fact(isec,ic,iem)   *lonlat_fac                   
                       
                       Cexist=.true.
                       exit
                    endif
                 enddo
                 if(.not.Cexist)then
                    !country not included yet. define it now:
                    nGridEmisCodes(i,j)=nGridEmisCodes(i,j)+1
                    if(nGridEmisCodes(i,j)>NCMAX)then
                       write(*,*) sub//&
                        "Too many emitter countries in one gridemiscell: ",&
                        me,i,j,nGridEmisCodes(i,j)
                       call StopAll("To many countries in one gridemiscell ")
                    endif
                    i_gridemis=nGridEmisCodes(i,j)
                    GridEmisCodes(i,j,i_gridemis)=Country(ic)%icode
                    GridEmis(isec,i,j,i_gridemis,iem)= &
                      cdfemis(i,j) * e_fact(isec,ic,iem)*lonlat_fac 
                 endif
                 
              end if
            end do !i
        end do !j
      end do ! variables 

     call check(nf90_close(ncFileID))

     CALL MPI_REDUCE(sumcdfemis_loc(1:NLAND),sumemis(1:NLAND),NLAND,&
                          MPI_REAL8,MPI_SUM,0,MPI_COMM_CALC,IERROR)
     i=max( 1, len_trim(fname) - 20 ) ! Just to get last bit of filename
     if(MasterProc) write(*,'(a,2x,a,i4,a,f12.3)') sub//'...'//trim(fname(i:)), &
           trim(EMIS_FILE(iem)), iem, ' DONE, UK now:', sumemis(27)

  end subroutine EmisGetCdf




  subroutine EmisGetASCII(iem, fname, emisname, sumemis_local, incl, nin, excl, nex)
    implicit none
    integer, intent(in) ::iem, nin, nex
    character(len=*),intent(in) :: fname, emisname, incl(*),excl(*)
    real,intent(inout), dimension(NLAND,NEMIS_FILE) &
        :: sumemis_local ! Sum of emissions per country
    integer ::i,j,k,n,ic,CC,isec,i_gridemis,found
    character(len=*), parameter ::  sub = 'EmisGetASCII:'
    logical :: Cexist
    real :: tmpsec(NSECTORS),duml,dumh
    real ::lonlat_fac(NSECTORS)
    character(len=len(fname)+10) :: newfname !need some extra characters if the new name gets longer

      call open_file(IO_EMIS, "r", fname, needed=.false., iostat=ios)

      if(ios/=0)then
         22 format(5A)
         if(MasterProc)write(*,22)'did not find ',trim(fname)
         !try to write pollutant with capital letters
         newfname=trim(fname)
         newfname=key2str(newfname,'gridsox','gridSOx')
         newfname=key2str(newfname,'gridnox','gridNOx')
         newfname=key2str(newfname,'gridco','gridCO')
         newfname=key2str(newfname,'gridvoc','gridNMVOC')
         newfname=key2str(newfname,'gridnh3','gridNH3')
         newfname=key2str(newfname,'gridpm25','gridPM25')
         newfname=key2str(newfname,'gridpmco','gridPMco')
         newfname=key2str(newfname,'gridpm10','gridPM10')
         if(MasterProc)write(*,22)'trying: ',trim(newfname)
         call open_file(IO_EMIS,"r",newfname,needed=.true.)
      endif
READEMIS: do   ! ************* Loop over emislist files *******************

            read(unit=IO_EMIS,fmt=*,iostat=ios) CC,i,j, duml,dumh,  &
                                    (tmpsec(isec),isec=1,NSECTORS)
            if ( ios <  0 ) exit READEMIS            ! End of file
            call CheckStop(ios > 0,"EmisGetASCII: ios error in emission file")

            do ic=1,NLAND
               if((Country(ic)%icode==CC))&
                    goto 543
            enddo 
            write(unit=errmsg,fmt=*) &
                   "COUNTRY CODE NOT RECOGNIZED OR UNDEFINED ", CC
            call CheckStop(errmsg)
            ic=0
543         continue
            
            if ( i  <=  0 .or. i  > IIFULLDOM  .or.   & 
                 j  <=  0 .or. j  > JJFULLDOM) cycle READEMIS
            
            i = i_local(i)     ! for SUB domain
            j = j_local(j)     ! for SUB domain

            if ( i  <=  0 .or. i  >  LIMAX .or.   & 
                 j  <=  0 .or. j  >  LJMAX .or.   &
                 ic <=  0 .or. ic >  NLAND .or.   &
                 ic == IC_NAT              .or.   &  ! Excludes DMS
                 (ic == IC_VUL .and. VOLCANOES_LL) )&! Excludes Volcanoes 
                                                     ! from gridSOx. Read from
                                                     ! VolcanoesLL.dat instead
             cycle READEMIS

             if ( trim ( emisname ) == "voc" ) tmpsec(11:11) = 0.0

             !exclude or include countries
             !could easily be optimized, by defining a country mask
             if(nin>0)then
                !1) check that country is in include list
                found=find_index(Country(ic)%code ,incl(1:nin),first_only=.true.)
                if(found<=0)cycle READEMIS!do not include
             endif
             if(nex>0)then
                !1) check that country is not in exclude list
                found=find_index(Country(ic)%code ,excl(1:nex),first_only=.true.)
                if(found>0)cycle READEMIS!exclude
             endif

             lonlat_fac=1.0
             if(N_femis_lonlat>0)then
                do i_femis_lonlat=1,N_femis_lonlat
                   if(glat(i,j)>femis_latmin(i_femis_lonlat).and.&
                        glat(i,j)<femis_latmax(i_femis_lonlat).and.&
                        glon(i,j)>femis_lonmin(i_femis_lonlat).and.&
                        glon(i,j)<femis_lonmax(i_femis_lonlat)     )then
                      lonlat_fac(:)=lonlat_fac(:)*e_fact_lonlat(:,i_femis_lonlat,iem) 
                   endif
                enddo
             endif

              !merge into existing emissions           
             Cexist=.false.
             do i_gridemis=1,nGridEmisCodes(i,j)
                if(GridEmisCodes(i,j,i_gridemis)==CC)then

                   GridEmis(:,i,j,i_gridemis,iem)=&
                         GridEmis(:,i,j,i_gridemis,iem)&
                        +e_fact(:,ic,iem) *  tmpsec(:)   *lonlat_fac(:)              

                   Cexist=.true.
                   exit
                endif
             enddo
             if(.not.Cexist)then
                !country not included yet. define it now:
                nGridEmisCodes(i,j)=nGridEmisCodes(i,j)+1
                if(nGridEmisCodes(i,j)>NCMAX)then
                   write(*,*) sub// &
                      " Too many emitter countries in one gridemiscell: ",&
                        me,i,j,CC,&
                       (GridEmisCodes(i,j,i_gridemis),i_gridemis=1,NCMAX)
                   call StopAll("To many countries in one gridemiscell ")
                endif
                i_gridemis=nGridEmisCodes(i,j)
                GridEmisCodes(i,j,i_gridemis)=CC
                GridEmis(:,i,j,i_gridemis,iem)=e_fact(:,ic,iem) *  tmpsec(:)*lonlat_fac(:)
             endif
            !for diagnostics, mass balance:
             sumemis_local(ic,iem)=sumemis_local(ic,iem)&
                  +0.001*sum(e_fact(:,ic,iem)*tmpsec(:)*lonlat_fac(:))

        end do READEMIS 
        
        close(IO_EMIS)
        ios = 0
  end subroutine EmisGetASCII

!-----------------------------------------------------------------------------

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

  !*** local variables ***
  integer            :: ie, iq, ic, iland1, iland2 & ! loop variables
                       ,inland                     & ! Country read from femis
                       ,isec, isec1 , isec2        & ! loop vars: emis sectors
                       ,nwords,ncols, n, oldn       ! No. cols. in "femis" 
  integer, parameter        :: NCOLS_MAX = 20  ! Max. no. cols. in "femis"
  integer, dimension(NEMIS_FILE) :: qc        ! index for sorting femis columns
  real, dimension(NCOLS_MAX):: e_f             ! factors read from femis
  real, dimension(NCOLS_MAX):: e_f_lonlat      ! factors read from femis in lonlat format
  character(len=200) :: txt                    ! For read-in 
  character(len=30), dimension(NCOLS_MAX)::  txtinwords ! to read lines
 !--------------------------------------------------------


  e_fact(:,:,:) = 1.0            !**** default value = 1 ****
  e_fact_lonlat(:,:,:) = 1.0            !**** default value = 1 ****

  associate ( debugm0 => ( DEBUG_GETEMIS .and. MasterProc ) )

  if( debugm0 ) write(*,*) "Enters femis", me
  call open_file(IO_EMIS,"r","femis.dat",needed=.false.)

  if ( ios == NO_FILE ) then
        ios = 0
        write( *,*) "WARNING: NO FEMIS FILE"
        return !*** if no femis file, e_fact=1 as default *** 
  endif
  call CheckStop( ios < 0 ,"EmisGet:ios error in femis.dat")


  ! Reads in the header line, e.g. name sec sox nox voc. 
  ! Pollutant names wil be checked against those defined in My_Emis_ml 

  read(unit=IO_EMIS,fmt="(a200)") txt
  if(debugm0)write(unit=6,fmt=*) "In femis, header0 is: ",  trim(txt)

  call wordsplit(txt,NCOLS_MAX,txtinwords,ncols,ios)
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
                if ( txtinwords(ic+2) == trim ( EMIS_FILE(ie) ) ) then
                    qc(ie) = ic
                    n = n + 1
                    if(debugm0)write(unit=6,fmt=*) "In femis: ", &
                       txtinwords(ic+2), " assigned to ", ie, EMIS_FILE(ie)
                  exit EMLOOP
                end if
      end do EMLOOP ! ie
       if (oldn == n .and.debugm0)   &
           write(unit=6,fmt=*) "femis: ",txtinwords(ic+2)," NOT assigned"
  end do COLS   ! ic

if ( n < NEMIS_FILE ) then
  print *, "FEMIS me n NEMIS_FILE", me, n, NEMIS_FILE
  call CheckStop( n < NEMIS_FILE , "EmisGet: too few femis items" )
end if

  
  n = 0
  N_femis_lonlat=0

  if(me==0)write(IO_LOG,55)'femis reductions:'
  READFILE: do  ! ************ read lines of femis ***************
       read(unit=IO_EMIS,fmt="(a200)",iostat=ios) txt
       if ( ios <  0 ) exit READFILE                   ! End of file
       call CheckStop( ios > 0 , "EmisGet: read error in femis" )
       call wordsplit(txt,NCOLS_MAX,txtinwords,nwords,ios)
       if ( nwords<3 ) cycle READFILE                   ! End of file
       55 format(A)
       if(me==0)write(IO_LOG,55)trim(txt)
       if(txtinwords(1)=='lonlat')then
!reductions defined with coordinates
       if(nwords<ncols+5)then
          if(me==0)write(*,*)'femis.dat not understood ',nwords,ncols+5,txt
          call CheckStop( nwords<ncols+5 , "EmisGet: read error in femis lonlat" )
       endif
!latmin,latmax,lonmin,lonmax
          N_femis_lonlat=N_femis_lonlat+1
          call CheckStop( N_femis_lonlat>MAXFEMISLONLAT, "EmisGet: increase MAXFEMISLONLAT" )
!         77 format(A,4F,I,20F)
          read(txt,*)txtinwords(1)&
          ,femis_lonmin(N_femis_lonlat),femis_lonmax(N_femis_lonlat),femis_latmin(N_femis_lonlat),femis_latmax(N_femis_lonlat)&
          ,isec,(e_f_lonlat(ic),ic=1,ncols)

!It is rather easy to get coordinates which are equals. 
!In order to avoid random decisions when this happens, we increase slightly the bounds:
          femis_lonmin(N_femis_lonlat)=femis_lonmin(N_femis_lonlat)-1.0E-10
          femis_lonmax(N_femis_lonlat)=femis_lonmax(N_femis_lonlat)+1.0E-10
          femis_latmin(N_femis_lonlat)=femis_latmin(N_femis_lonlat)-1.0E-10
          femis_latmax(N_femis_lonlat)=femis_latmax(N_femis_lonlat)+1.0E-10

!"normalize" longitudes to the interval -180,180
          if(femis_lonmin(N_femis_lonlat)>180.0)femis_lonmin(N_femis_lonlat)=femis_lonmin(N_femis_lonlat)-360.0
          if(femis_lonmin(N_femis_lonlat)<-180.0)femis_lonmin(N_femis_lonlat)=femis_lonmin(N_femis_lonlat)+360.0
          if(femis_lonmax(N_femis_lonlat)>180.0)femis_lonmax(N_femis_lonlat)=femis_lonmax(N_femis_lonlat)-360.0
          if(femis_lonmax(N_femis_lonlat)<-180.0)femis_lonmax(N_femis_lonlat)=femis_lonmax(N_femis_lonlat)+360.0

!There will still be problem when trying to "cross" the 180 degree longitude, therefore we put a test here:
          call CheckStop(femis_lonmin(N_femis_lonlat)>femis_lonmax(N_femis_lonlat),&
               "femislonlat: crossing 180 degrees longitude not allowed")
       else
!      read(unit=IO_EMIS,fmt=*,iostat=ios) inland, isec, (e_f(ic),ic=1,ncols)
!      if ( ios <  0 ) exit READFILE                   ! End of file
!      call CheckStop( ios > 0 , "EmisGet: read error in femis" )

          read(txt,fmt=*,iostat=ios) inland, isec, (e_f(ic),ic=1,ncols)
          
          n = n + 1
          if(debugm0) then
             write(unit=6,fmt=*) "FEMIS READ", inland, &
                  isec, (e_f(ic),ic=1,ncols)
             write(unit=6,fmt="(2a,I3,a,i3,a)") " Emission factors from femis.dat, ",&
                  "landcode =", inland, ",  sector code =",isec, &
                  " (sector 0 applies to all sectors) :"
             write(unit=6,fmt="(a,14(a,a,F5.2,a))") " ", (trim(txtinwords(qc(ie)+2)),&
                  " =",e_f(qc(ie)), ",  ", ie=1,NEMIS_FILE-1), &
                  (trim(txtinwords(qc(ie)+2))," =",e_f(qc(ie))," ", &
                  ie=NEMIS_FILE,NEMIS_FILE)
          end if
          
          if (inland == 0 ) then     ! Apply factors to all countries
             iland1 = 1 
             iland2 = NLAND
             
          else                       ! Apply factors to country "inland"
             
             ! find country number  corresponding to index as written in emisfile
             do iland1=1,NLAND
                if(Country(iland1)%icode==inland) goto 544
             enddo
             
             if(MasterProc) write(*,*)'COUNTRY CODE NOT RECOGNIZED',inland
             
             iland1 = 0
             iland2 =-1
544          continue
             if(iland1/=0)   iland2 = iland1
          end if
       endif

       if (isec == 0 ) then       ! All sectors
          isec1 = 1
          isec2 = NSECTORS
       elseif (isec==100) then    ! Anthropogenic scenario
          !if you need this option, then just uncomment and check that
          !ANTROP_SECTORS is set according to your categories (10 for SNAP)
          call CheckStop(SECTORS_NAME/='SNAP',"Anthropogenic not compatible with non-SNAP")
          isec1 = 1
          isec2 = ANTROP_SECTORS
       else                       ! one sector: isec
          isec1 = isec
          isec2 = isec
       end if
       do ie = 1,NEMIS_FILE
          if(txtinwords(1)=='lonlat')then
             do isec = isec1, isec2
                e_fact_lonlat(isec,N_femis_lonlat,ie) = e_fact_lonlat(isec,N_femis_lonlat,ie) * e_f_lonlat( qc(ie) )
             end do !isec
          else
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
          endif
       end do !ie
          
  enddo READFILE ! Loop over femis

  close(IO_EMIS)

  if( MasterProc) write(unit=6,fmt=*) "femis, read ", n, "country code lines from femis"
  if( MasterProc) write(unit=6,fmt=*) "femis, read ", N_femis_lonlat, "lonlat lines from femis"
  if(debugm0) then
    ! Extra checks
     write(unit=6,fmt=*) "DEBUG_EMISGET: UK femis gives: "
     write(unit=6,fmt="(6x, 30a10)") (EMIS_FILE(ie), ie=1,NEMIS_FILE)
     do isec = 1, NSECTORS
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
         allocate(emis_hprofile(nemis_hprofile+1,N_HFAC),stat=allocerr)
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
         if(snap > N_HFAC)then
            if(me==0)write(*,*)N_HFAC,' sector classes defined, but found ',snap
            call CheckStop(snap > N_HFAC,"EmisGet: sector class out of bounds")
         endif
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
      allocate(emis_kprofile(nemis_kprofile,N_HFAC),stat=allocerr)
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

      allocate(emis_kprofile(nemis_kprofile,N_HFAC),stat=allocerr)
      emis_kprofile=0.0

      !convert height (given as pressure) distribution into model level distribution
      !ext meaning using levels from file (EmisHeights.txt)
      k1_ext(KMAX_BND)=0
      do isec=1,N_HFAC! could put inside but easier for debugging to put here now
         if(DEBUG_GETEMIS.and.MasterProc) write(*,*)'Sector ',isec
         do k=KMAX_BND-1,max(1,KMAX_BND-nemis_kprofile),-1

!count all contributions to model level k
!i.e. between P_emep(k+1) and P_emep(k)

            P_emep=A_bnd(k)+B_bnd(k)*Pref !Pa
            if(DEBUG_GETEMIS.and.MasterProc) write(*,fmt="(A,I3,F10.2)")'vert_inter: P_emep',k,P_emep
            !largest available P_ext smaller than P_emep (if possible)
            !k1_ext(k) is the external layer just below P_emep(k)
            k1_ext(k)= 0 !start at surface, and go up until P_emep
            do k_ext=1,nemis_hprofile
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

            !add level just below P_emep(k), if not already counted, above k1_ext(k)  below P_emep; and must exist
            if(emis_P_level(k1_ext(k+1)+1)>P_emep.and. k1_ext(k)<nemis_hprofile )then
               frac=(emis_P_level(k1_ext(k))-P_emep)/(emis_P_level(k1_ext(k))-emis_P_level(k1_ext(k)+1))
               emis_kprofile(KMAX_BND-k,isec)=emis_kprofile(KMAX_BND-k,isec)+frac*emis_hprofile(k1_ext(k)+1,isec)
               if(DEBUG_GETEMIS.and.MasterProc) write(*,fmt="(A,I5,6F10.2)")'adding last fraction of level',k1_ext(k)+1,&
                    frac,emis_hprofile(k1_ext(k)+1,isec),emis_P_level(k1_ext(k+1)+1),emis_P_level(k1_ext(k)),P_emep,&
                    (emis_P_level(k1_ext(k))-emis_P_level(k1_ext(k)+1))
            endif

         enddo
      enddo
      
      CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)!so that output from different CPU does not get mixed up
      
      if(MasterProc)then
         write(*,*)'Distribution of emission into levels:'
         do isec=1,N_HFAC! could put inside but easier for debugging to put here now
            write(*,fmt="(A,I5,A,20F6.3)")'sector: ',isec,' fractions: ',(emis_kprofile(k,isec),k=1,nemis_kprofile)
         enddo
      endif

   endif

!check normalization of distributions:
   do isec=1,N_HFAC
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
!    real emisfrac(NRCSPLIT,N_SPLIT,NLAND)
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
  real, dimension(NSPEC_ADV,N_SPLIT,NLAND) :: tmp_emisfrac
  real, dimension(NSPEC_ADV) :: tmp_emis_masscorr
  integer, dimension(NSPEC_ADV) :: tmp_iqrc2itot !maps from iqrc 
  real, dimension(NMAX ) :: tmp 
  real     :: sumtmp
  integer  :: nsplit   &         ! No.columns data to be read
             ,iland,isec,i,n,nn, allocerr,iland_icode
  integer  :: idef              ! Set to   0  for defaults, 1 for specials
  integer  :: iland1, iland2    ! loop variables over countries
  logical  :: defaults          ! Set to true for defaults, false for specials
  logical  :: debugm            ! debug flag on master proc 
  character(len=*), parameter :: dtxt = 'EmisGet:'
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
22      format(25A)
          write(unit=6,fmt=22) "Splitting ", trim(EMIS_FILE(ie)), &
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
           read(unit=txtinput,fmt=*,iostat=ios)  iland_icode, isec, (tmp(i),i=1, nsplit)
           if( MasterProc .and. ios /= 0 ) then
               print *, "ERROR: EmisGet: Failure reading emispslit file"
               print *, "Expecting to split into nsplit=", nsplit
               print *, "reading this record:", trim(txtinput)
              call CheckStop( ios ,"EmisGet: ios error on "//trim(fname) )
           end if
           if(iland_icode==0)then
              iland=0!special meaning
           else
              iland=find_index(iland_icode,Country(:)%icode)!find country array index from code. 
              if( iland < 1) then
                 if(MasterProc) then
                    print "(a,2i6,a,1x,a)", dtxt//" ILAND NOT DEF ", iland, &
                     iland_icode, trim(txtinput), trim(EMIS_FILE(ie))
                 end if
                 call StopAll ( dtxt//' UNDEFINED iland in '//&
                 trim(EMIS_FILE(ie))//':'//trim(txtinput) )
                 cycle READ_DATA
              end if
          endif

           n = n + 1
           if (debugm ) then
               write(6,"(a,i3,a,3i3,50f8.2)") "Splits: ",  n, trim(fname),&
                  iland, isec, nsplit, tmp(1:nsplit)
           end if

           !/... some checks:
           sumtmp = sum( tmp(1:nsplit) )
          !GEOS-CHEM HAD 10 mols O3 + 1 molec HNO3, so >> 100
          ! if( iland_icode ==350 ) sumtmp = 100.0
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

                !*** assign and convert from percent to fractions: ***

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

       call CheckStop(  defaults .and. n  /=  N_SPLIT, &
                        "ERROR: EmisGet: defaults .and. n  /=  N_SPLIT" )

       if (debugm ) write(*,*) "Read ", n, " records from ",fname

    end do IDEF_LOOP 
  end do ! ie
  ios = 0
  
  ! By now we should have found all the species required for the 
  ! chemical scheme:

  do ie = 1, NEMIS_SPECS
    if ( DEBUG_GETEMIS .and. MasterProc .and. ( EmisSpecFound(ie) .eqv. .false.) ) then
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
  allocate(emisfrac(nrcemis,N_SPLIT,NLAND),stat=allocerr)
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
               if((Country(ic)%icode==iic))&
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

              call GridAllocate("ROAD"// trim ( emisname ),i,j,Country(ic)%icode,NCMAX, &
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
