! <EmisGet_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module EmisGet_mod

use CheckStop_mod,     only: CheckStop, StopAll, check=>CheckNC
use ChemDims_mod,      only: NSPEC_ADV, NSPEC_TOT, NEMIS_File, NEMIS_Specs
use ChemSpecs_mod,     only: species 
use Config_module,     only: NPROC, MasterProc,USES,step_main,&
                             KMAX_MID,KMAX_BND, Pref,&
                             SEAFIX_GEA_NEEDED, & ! only if emission problems over sea
                             IIFULLDOM,JJFULLDOM, SECTORS_NAME, &
                             SplitSpecialsFile,SplitDefaultFile,EmisHeightsFile,femisFile,&
                             startdate
use Country_mod,       only: NLAND, IC_NAT, IC_VUL, IC_NOA, Country, &
                             ! NMR-NH3 specific variables (hb NH3Emis)
                             IC_NMR,IC_DUMMY 
use Debug_module,      only: DEBUG
use EmisDef_mod,       only: NSECTORS, ANTROP_SECTORS, NCMAX, & 
                            N_HFAC,N_SPLIT, EMIS_FILE, & 
                            VOLCANOES_LL, &
                          ! NMR-NH3 specific variables (for FUTURE )
                            NH3EMIS_VAR,dknh3_agr,ISNAP_AGR,ISNAP_TRAF, &
                            NROADDUST, &
                            gnfr2snap,snap2gnfr&
                            ,cdfemis,sumcdfemis,nGridEmisCodes,GridEmisCodes&
                            ,GridEmis,gridrcemis, Emis_mask, MASK_LIMIT&
                            ,landcode,nlandcode,MAXFEMISLONLAT,N_femis_lonlat &   
                            ,femis_lonlat_internal & 
                            ,Emis_field, NEmis_id, Emis_id, NEmis_sources&
                            ,EmisFiles, NEmisFile_sources, Emis_source &
                            ,NEmis_sourcesMAX
use GridAllocate_mod,  only: GridAllocate
use GridValues_mod,    only: debug_proc,debug_li,debug_lj,i_fdom,j_fdom,i_local
use GridValues_mod,    only: glon, glat, A_bnd, B_bnd,j_local
use Io_mod,            only: open_file, NO_FILE, ios, IO_EMIS, &
                             Read_Headers, read_line, PrintLog
use Io_Progs_mod,      only: datewrite
use KeyValueTypes,     only: KeyVal
use MPI_Groups_mod  , only : MPI_BYTE, MPI_REAL8, MPI_DOUBLE_PRECISION, MPI_SUM&
                             , MPI_INTEGER, MPI_COMM_CALC, IERROR
use NetCDF_mod,      only  : ReadField_CDF, check_lon_lat, ReadTimeCDF, &
                             GetCDF_modelgrid, make_gridresolution

use OwnDataTypes_mod,  only: TXTLEN_NAME, TXTLEN_FILE, Emis_id_type, &
                             EmisFile_id_type, Emis_sourceFile_id_type
use Par_mod,           only: LIMAX, LJMAX, limax, ljmax, me
use SmallUtils_mod,    only: wordsplit, find_index, key2str
use netcdf,            only: NF90_OPEN,NF90_NOERR,NF90_NOWRITE,&
                             NF90_INQUIRE,NF90_INQUIRE_VARIABLE,NF90_CLOSE,&
                             nf90_global,nf90_get_var,nf90_get_att
use PhysicalConstants_mod,  only : PI, EARTH_RADIUS
use TimeDate_mod, only     : date
use TimeDate_ExtraUtil_mod, only : nctime2date,date2nctime, date2string
use Timefactors_mod,   only: fac_ehh24x7 

implicit none
private

 ! subroutines:
public  :: Emis_init_GetCdf  ! define constant parameters
public  :: Emis_GetCdf       ! new formats
public  :: EmisGetCdf        ! cdfemis
public  :: EmisGetASCII      !new version of "EmisGet" which does not use fulldomain arrays
public  :: EmisSplit         ! => emisfrac, speciation of voc, pm25, etc.
public  :: EmisHeights       ! => nemis_kprofile, emis_kprofile
                             !     vertical emissions profile
public  :: RoadDustGet       ! Collects road dust emission potentials
public  :: femis             ! Sets emissions control factors 
public  :: make_iland_for_time !make land indices for timefactors
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
real, public, save, dimension(MAXFEMISLONLAT) :: femis_latmin,femis_latmax,femis_lonmin,femis_lonmax,femis_lonlat_ic

! emisfrac is used at each time-step of the model run to split
! emissions such as VOC, PM into species. 

integer, public, parameter :: NMAX = NSPEC_ADV 
integer, public, save :: nrcemis, nrcsplit
integer, public, dimension(NEMIS_FILE) , save :: emis_nsplit
real, public,allocatable, dimension(:,:,:), save :: emisfrac
integer, public,allocatable, dimension(:), save :: iqrc2itot
integer, public,allocatable, dimension(:), save :: iqrc2iem
integer, public, dimension(NSPEC_TOT), save :: itot2iqrc
integer, public, dimension(NSPEC_TOT,NEMIS_FILE) :: iemsplit2itot !maps from split and iem to itot 
integer, public, dimension(NEMIS_FILE), save :: Emis_MolWt
real, public,allocatable, dimension(:), save :: emis_masscorr
real, public,allocatable, dimension(:), save :: roaddust_masscorr

! vertical profiles for SNAP emis, read from EmisHeightsFile
integer, public, save :: nemis_kprofile
real, public,allocatable, dimension(:,:), save :: emis_kprofile
real, public,allocatable, dimension(:,:), save :: emis_hprofile

! some common variables
character(len=TXTLEN_FILE), private :: fname             ! File name
character(len=80), private :: errmsg

! Import list of the emitted species we need to find in the 
! emissplit files.
include 'CM_EmisSpecs.inc'
logical, dimension(NEMIS_SPECS) :: EmisSpecFound = .false.
integer, private ::i_femis_lonlat

contains
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Emis_GetCdf(EmisFile, Emis_source, Emis_XD, date_wanted)
    type(Emis_id_type),intent(inout) :: Emis_source 
    type(EmisFile_id_type),intent(inout) ::  EmisFile
    real, intent(out), dimension(*) :: Emis_XD
    type(date), intent(in) :: date_wanted
    real :: date_wanted_in_days, TimesInDays(1)
    integer :: record
    logical, save :: dbg= .false., first_call = .true. 
    character(len=*), parameter :: dtxt = 'Emis_GetCdf:'

    if ( first_call ) then
      dbg =  ( MasterProc .and. DEBUG%GETEMIS )
      first_call = .false.
    end if

    fname = date2string(EmisFile%filename,date_wanted,mode='YMDH')

    if(EmisFile%periodicity == 'yearly' .or. EmisFile%periodicity == 'once')then
       !assumes only one record to read
       record = 1
    else if(EmisFile%periodicity == 'monthly')then
        !assumes 12 records, one for each month
       record = date_wanted%month
    else
       !the correct time must be written in the file
       call date2nctime(date_wanted,date_wanted_in_days)
       call ReadTimeCDF(fname,TimesInDays,record,date_wanted_in_days)
       
       call nctime2date(EmisFile%end_of_validity_date, TimesInDays(1))
       if(me==0 .and. (step_main<10 .or. DEBUG%EMISSIONS))&
            write(*,*)record,'ENDOFVAL ', EmisFile%end_of_validity_date
    endif

    if( dbg ) write(*,*) dtxt//'Reading '//trim(fname)
    if(trim(EmisFile%projection) == 'native')then
       if(me==0.and. (step_main==1 .or. DEBUG%EMISSIONS))&
            write(*,*)'reading  new '//trim(Emis_source%varname)//' from native grid, record ',record
       if(Emis_source%is3D)then
          call GetCDF_modelgrid(Emis_source%varname,fname,Emis_XD,&
                            Emis_source%kstart, Emis_source%kend,record,1,&
                            i_start=Emis_source%istart,j_start=Emis_source%jstart,&
                            reverse_k=Emis_source%reversek,needed=.true.)
       else
          call GetCDF_modelgrid(Emis_source%varname,fname,Emis_XD,1,1,record,1,&
                            i_start=Emis_source%istart,j_start=Emis_source%jstart,&
                            needed=.true.)
       endif
    else
       if(me==0 .and. (step_main==1 .or. DEBUG%EMISSIONS))&
            write(*,*)trim(Emis_source%varname)//' reading new emis from '//trim(fname)//', record ',record
       if(Emis_source%units(1:5) == 'kt/m2'  &
            .or. Emis_source%units(1:9) == 'tonnes/m2'  &
            .or. Emis_source%units(1:5) == 'kg/m2' &
            .or. Emis_source%units(1:4) == 'g/m2'  &
            .or. Emis_source%units(1:5) == 'mg/m2')then
          !per area units
          call ReadField_CDF(fname,Emis_source%varname,Emis_XD,record,&
               known_projection=trim(EmisFile%projection),&
               interpol='conservative',&
               Grid_resolution_in = EmisFile%grid_resolution,&
               needed=.true.,UnDef=0.0,&
               debug_flag=.false.)
       else  if(Emis_source%units == 'kt' .or. Emis_source%units == 'kt/s' &
            .or. Emis_source%units == 'kt/month' .or. Emis_source%units == 'kt/year' &
            .or. Emis_source%units == 'tonnes' .or. Emis_source%units == 'tonnes/s' &
            .or. Emis_source%units == 'tonnes/month' .or. Emis_source%units == 'tonnes/year' &
            .or. Emis_source%units == 'kg' .or. Emis_source%units == 'kg/s' &
            .or. Emis_source%units == 'kg/month' .or. Emis_source%units == 'kg/year' &
            .or. Emis_source%units == 'g' .or. Emis_source%units == 'g/s' &
            .or. Emis_source%units == 'g/month' .or. Emis_source%units == 'g/year' &
            .or. Emis_source%units == 'mg' .or. Emis_source%units == 'mg/s' &
            .or. Emis_source%units == 'mg/month' .or. Emis_source%units == 'mg/year' &
            .or. Emis_source%units == 'g/h' .or. Emis_source%units == 'mg/h')then
          !per gridcell unit
          if(me==0 .and. (step_main<10 .or. DEBUG%EMISSIONS))&
               write(*,*)'reading emis '//trim(Emis_source%varname)//' from '//trim(fname)//', proj ',trim(EmisFile%projection),', res ',EmisFile%grid_resolution
         call ReadField_CDF(fname,Emis_source%varname,Emis_XD,record,&
               known_projection=trim(EmisFile%projection),&
               interpol='mass_conservative',&
               Grid_resolution_in = EmisFile%grid_resolution,&
               needed=.true.,UnDef=0.0,&
               debug_flag=.false.)          
       else
          call StopAll("EmisGet: Unit for emissions not recognized: "//trim(Emis_source%units)//' '//trim(Emis_source%varname))
       endif
    endif

  end subroutine Emis_GetCdf
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine EmisGetCdf(iem, fname, sumemis_local, &
       Emis, EmisCodes, nEmisCodes, nstart,&
       incl, nin, excl, nex, use_lonlat_femis, &
       set_mask,use_mask,pollName, fractionformat, type)

    !read in emissions in fraction format and add results into
    !Emis, nEmisCodes and nEmisCodes
    !NB: landcode and nlandcode arrays are local

    implicit none
    integer, intent(in) ::iem, nin, nex,nstart
    character(len=*),intent(in) :: fname, incl(*),excl(*),pollName(*),type
    real,intent(inout) ::Emis(NSECTORS,LIMAX,LJMAX,NCMAX,NEMIS_FILE)
    integer,intent(inout) ::nEmisCodes(LIMAX,LJMAX)
    integer,intent(inout) ::EmisCodes(LIMAX,LJMAX,NCMAX)
    ! Sum of emissions per country
    real,intent(inout), dimension(NLAND,NEMIS_FILE) :: sumemis_local 
    logical, intent(in) :: use_lonlat_femis,set_mask,use_mask
    logical, intent(in) :: fractionformat

    real :: fractions(LIMAX,LJMAX,NCMAX),Reduc(NLAND)
    integer :: landcode(LIMAX,LJMAX,NCMAX),nlandcode(LIMAX,LJMAX)
    real :: lonlat_fac, fac, file_fac

    integer ::i,j,n,ic,i_cc,found, isec, sec_ix
    logical :: Cexist
    integer :: nwords, err, status
    integer :: ncFileID, nDimensions,nVariables,nAttributes,timeDimID,varid,xtype,ndims
    character(len=100) :: code, ewords(10), varname, cdfvarname, specname
    logical :: mappingGNFR2SNAP = .false., foundEmis_id
    integer :: iem_used, ix_Emis

    if(.not. fractionformat)then
       !not fraction format
       !loop over all variables
       status=nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncFileID)

       if( status /= nf90_noerr ) then
          if( MasterProc )write(*,*) "ERROR: EmisGetCdf - couldn't open "//trim(fName)
          CALL MPI_FINALIZE(IERROR)
          stop
       end if
       call check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,timeDimID))
       nlandcode=1
       fractions=1.0
       if(SECTORS_NAME=='GNFR'.and.type == "SNAPsectors")then
          call CheckStop('GNFR sector cannot mix with SNAP sectors for this format')
       else if(SECTORS_NAME=='SNAP'.and.type == "GNFRsectors")then
          call CheckStop('SNAP sector cannot mix with GNFR sectors for this format')
       else
          !only case implemented for not fractions
       endif
    else
       nVariables = NSECTORS*NEMIS_FILE !in fraction format we loop over all species and sectors
    endif

    mappingGNFR2SNAP = .false.
    do varid=1,nVariables
       sec_ix = 1
762    continue !if several GNFR counterparts to one SNAP
       iem_used = iem !default
       file_fac = 1.0 !default
       foundEmis_id = .false.
       if(.not. fractionformat)then
          call check(nf90_Inquire_Variable(ncFileID,varid,cdfvarname,xtype,ndims))
          if(me==0)write(*,*)'reading ',trim(cdfvarname)
          ewords=''
          if( index(cdfvarname, "Emis:") >0 )then
             ! Emission terms look like, e.g. Emis:FR:snap:7
             ! or Emis:FR:snap:7:spec:sox:fac:1.0
             call wordsplit(cdfvarname,8,ewords,nwords,err,separator=":")       
             if( ewords(3) /= "snap" ) cycle ! ONLY SNAP coded for now     
             read(ewords(4),"(I6)") isec       
             code = ewords(2)    ! Country ISO

             if( ewords(5) == "spec" )then !optional
                iem_used = find_index(ewords(6),EMIS_FILE(:))
                if(iem_used<0)then
                   iem_used = find_index(ewords(6),species(:)%name)
                   if(iem_used<0)then
                      if(MasterProc)write(*,*)trim(ewords(6))//' NOT FOUND'
                      cycle
                   else
                      if(MasterProc)write(*,*)'found '//species(iem_used)%name//' for '//trim(ewords(6))
                      !see if case has been found before
                      !should search only among NEmis_id first!
                      ix_Emis = find_index(ewords(6),Emis_id%species(:))
                      if(ix_Emis>0)then
                         if(MasterProc)write(*,*)trim(ewords(6))//' already defined'
                      else
                          if(MasterProc)write(*,*)NEmis_id,' defining new id '//trim(ewords(6))
                          NEmis_id = NEmis_id + 1
                          ix_Emis = NEmis_id
                          Emis_id(ix_Emis)%species = trim(ewords(6))
                      endif
                      iem_used = 1
                      foundEmis_id = .true.
                   endif
                endif
                !could add case for any species
             endif
             if( ewords(7) == "fac" )then !optional, but "spec" must be used if "fac" is set
                file_fac=-999.9
                read(ewords(8),*)file_fac 
                 if(abs(file_fac+999.9)<0.1)then !not good test
                   if(MasterProc)write(*,*)trim(ewords(8))//' NOT UNDERSTOOD AS REAL'
                   cycle
                endif               
             endif

             cdfemis = 0.0 ! safety, shouldn't be needed though
             if(pollName(1)/='NOTSET')then
                if(all(pollName(1:20)/=trim(EMIS_FILE(iem_used))))cycle      
                if(Masterproc)write(*,"(A)")'reading '//trim(EMIS_FILE(iem_used))//' from '//trim(fname)
             end if             
             call ReadField_CDF(fname,cdfvarname,cdfemis,nstart=nstart,&
                  interpol='mass_conservative',&
                  needed=.false.,UnDef=0.0,&
                  debug_flag=.false.,ncFileID_given=ncFileID)
             if( maxval(cdfemis ) < 1.0e-10 ) cycle ! Likely no emiss in domain
          else
             !other allowed cases to be implemented
             cycle
          endif
       else
          !fraction format
          isec=mod((varid-1),NSECTORS)+1
          iem_used=(varid-1)/NSECTORS+1
          if(SECTORS_NAME=='GNFR'.and.type == "SNAPsectors")then
             if(gnfr2snap(isec)<=0)cycle                   
             write(varname,"(A,I2.2)")trim(EMIS_FILE(iem_used))//'_sec',gnfr2snap(isec)
             if(me==0.and.iem_used==1)write(*,*)'WARNING, mapping snap sector ',gnfr2snap(isec),'onto gnfr',isec
          else if(SECTORS_NAME=='SNAP'.and.type == "GNFRsectors")then
             mappingGNFR2SNAP = .true.
             if(snap2gnfr(isec,sec_ix)<=0)cycle
             if(me==0.and.iem_used==1)write(*,*)'WARNING, mapping gnfr sector ',snap2gnfr(isec,sec_ix),'onto snap',isec
             write(varname,"(A,I2.2)")trim(EMIS_FILE(iem_used))//'_sec',snap2gnfr(isec,sec_ix)
          else
             write(varname,"(A,I2.2)")trim(EMIS_FILE(iem_used))//'_sec',isec
          endif
          
          if(pollName(1)/='NOTSET')then
             cdfemis = 0.0 ! safety, shouldn't be needed though
             if(all(pollName(1:20)/=trim(EMIS_FILE(iem_used))))cycle      
             if(Masterproc.and.isec==1)write(*,"(A)")'reading '//trim(EMIS_FILE(iem_used))//' from '//trim(fname)
          end if
          Reduc=e_fact(isec,:,iem_used)          
          call ReadField_CDF(trim(fname),varname,cdfemis(1,1),nstart=nstart,&
               interpol='mass_conservative',fractions_out=fractions,&
               CC_out=landcode,Ncc_out=nlandcode,Reduc=Reduc,needed=.false.,debug_flag=.false.,&
               Undef=0.0)
       endif
       !end reading of data

       !apply mask and femis_lonlat (femis "not lonlat" is already applied)
       do j=1,ljmax
          do i=1,limax                   
             if(foundEmis_id)then
                !for now no reductions applied
                Emis_field(i,j,ix_Emis) = Emis_field(i,j,ix_Emis) +&
                     cdfemis(i,j) * file_fac
             else
                if(use_mask)then
                   if(Emis_mask(i,j))cdfemis(i,j)=0.0
                endif
                if(set_mask)then
                   if(cdfemis(i,j)>MASK_LIMIT)Emis_mask(i,j)=.true.
                endif
                if(cdfemis(i,j)<100.0*MASK_LIMIT)cycle !important, because otherwise get too many countries per gridcell!
                if(nlandcode(i,j)>NCMAX)then
                   write(*,*)"GetCdfFrac: Too many emitter countries in one gridcell: ",&
                     me,i,j,nlandcode(i,j),trim(fname),trim(varname)
                   call StopAll("To many countries in one gridcell ")
                end if
                do n=1,nlandcode(i,j)
                   if(fractionformat)then
                      ic=find_index(landcode(i,j,n),Country(:)%icode)
                   else
                      ic=find_index(code,Country(:)%code)
                      landcode(i,j,n) = Country(ic)%icode
                   endif
                   
                   if(ic>NLAND .or. ic<1)then
                      write(*,*)"COUNTRY CODE NOT RECOGNIZED OR UNDEFINED: ",ic,landcode(i,j,n)
                      call StopAll("COUNTRY CODE NOT RECOGNIZED ")
                   end if

                   !exclude or include countries
                   !could easily be optimized, by defining a country mask
                   if(nin>0)then
                      !1) check that country is in include list
                      found=find_index(Country(ic)%code ,incl(1:nin),first_only=.true.)
                      if(found<=0)cycle!do not include
                   end if
                   if(nex>0)then
                      !1) check that country is not in exclude list
                      found=find_index(Country(ic)%code ,excl(1:nex),first_only=.true.)
                      if(found>0)cycle!exclude
                   end if
                   
                   
                   lonlat_fac=1.0
                   if(N_femis_lonlat>0 .and. use_lonlat_femis)then
                      do i_femis_lonlat=1,N_femis_lonlat
                         if(glat(i,j)>femis_latmin(i_femis_lonlat).and.&
                              glat(i,j)<femis_latmax(i_femis_lonlat).and.&
                              glon(i,j)>femis_lonmin(i_femis_lonlat).and.&
                              glon(i,j)<femis_lonmax(i_femis_lonlat).and.&
                              (femis_lonlat_ic(i_femis_lonlat)==0 .or. &
                              femis_lonlat_ic(i_femis_lonlat)==landcode(i,j,n)) )then
                            ! IN BOX (DSHK testing):
                            if ( femis_lonlat_internal(i_femis_lonlat) ) &
                              lonlat_fac=lonlat_fac*e_fact_lonlat(isec,i_femis_lonlat,iem_used) 
                         else if ( femis_lonlat_internal(i_femis_lonlat ) .eqv. .false. ) then
                              !HK - apply functions outside box
                              lonlat_fac=lonlat_fac*e_fact_lonlat(isec,i_femis_lonlat,iem_used)
                         end if
                      end do
                   end if
                   fac = file_fac*lonlat_fac
                   if(.not. fractionformat) fac = fac*e_fact(isec,ic,iem_used)
                   
                   !merge into existing emissions
                   Cexist=.false.
                   do i_cc=1,nEmisCodes(i,j)
                      if(EmisCodes(i,j,i_cc)==landcode(i,j,n))then
                         
                         Emis(isec,i,j,i_cc,iem_used)=Emis(isec,i,j,i_cc,iem_used)&
                           +fractions(i,j,n)*cdfemis(i,j) * fac                     

                         Cexist=.true.
                         exit
                      end if
                   end do
                   if(.not.Cexist)then
                      !country not included yet. define it now:
                      nEmisCodes(i,j)=nEmisCodes(i,j)+1
                      if(nEmisCodes(i,j)>NCMAX)then
                         write(*,*)"Too many emitter Countries in one emiscell: ",&
                              me,i,j,nEmisCodes(i,j),trim(fname),trim(varname),(EmisCodes(i,j,i_cc),i_cc=1,nEmisCodes(i,j)-1)
                         call StopAll("To many countries in one emiscell ")
                      end if
                      i_cc=nEmisCodes(i,j)
                      EmisCodes(i,j,i_cc)=landcode(i,j,n)
                      Emis(isec,i,j,i_cc,iem_used)=fractions(i,j,n)*cdfemis(i,j)*fac  
                   end if
                   sumemis_local(ic,iem_used)=sumemis_local(ic,iem_used)&
                        +0.001*fractions(i,j,n)*cdfemis(i,j)*fac  !for diagnostics, mass balance
                end do
             endif
          end do
       end do
                
       if(mappingGNFR2SNAP)then
          !some SNAP sectors have several GNFR counterparts. Must take all of them
          if(sec_ix<3)then
             if(snap2gnfr(isec,sec_ix+1)>0)then
                sec_ix = sec_ix + 1
                goto 762 
             endif
          endif
       endif
     enddo
    if(.not.fractionformat) call check(nf90_close(ncFileID))

  end subroutine EmisGetCdf

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Emis_init_GetCdf(EmisFile_in, EmisFile, names_in, nnames)

    !read in emissions from one file and set parameters

    implicit none
    character(len=*),intent(in) :: names_in(*)
    type(Emis_sourceFile_id_type),intent(in) :: EmisFile_in
    type(EmisFile_id_type),intent(inout) ::  EmisFile
    integer,intent(in) :: nnames
    character(len=100) :: projection, default_projection, name, cdfvarname, cdfspecies
    real :: factor, default_factor
    character(len= TXTLEN_FILE) :: fname
    character(len=30)  :: lon_name, lat_name
    integer :: n, nn, i, ix, varid, status, sector, iem, isec, iqrc, itot, f
    integer :: nDimensions, nVariables, nAttributes, xtype, ndims
    real :: x, resolution, default_resolution, Rlat(2)
    integer :: ncFileID, nemis_old, lonVarID, latVarID, dimids(10),dims(10)
    real, allocatable ::Rlat2D(:,:)
    character(len=*), parameter :: dtxt='Em_inicdf:'
    integer :: countrycode
    logical :: apply_femis

    fname=trim(date2string(EmisFile_in%filename,startdate,mode='YMDH'))
    status=nf90_open(path = trim(fname), mode = nf90_nowrite, ncid = ncFileID)
    if( status /= nf90_noerr ) then
       if( MasterProc )write(*,*)"ERROR: EmisGetCdf - couldn't open "//trim(fname)
       CALL MPI_FINALIZE(IERROR)
       stop
    end if


   !-------------
    associate ( debugm0 => ( DEBUG%GETEMIS .and. MasterProc ) )
   !-------------
    if ( debugm0 ) write(*,*) dtxt//'Start File:',trim(fname)
    
    if(EmisFile_in%projection /= 'native')then
       default_projection = 'Unknown'
       status = nf90_get_att(ncFileID, nf90_global,"projection", projection)
       if(status==nf90_noerr)then
          default_projection = trim(projection)
       endif
       
       default_resolution = 0.0
       status = nf90_get_att(ncFileID, nf90_global,"Grid_resolution", resolution)
       if(status==nf90_noerr)then
          default_resolution = resolution
       else
          call make_gridresolution(ncFileID, default_resolution)
       endif
    endif
       
    default_factor = 1.0
    status = nf90_get_att(ncFileID, nf90_global,"factor", factor)
    if(status==nf90_noerr)then
       default_factor = factor
    endif
    
    status = nf90_get_att(ncFileID,nf90_global,"sectorsName", name) !SNAPsectors or GNFRsectors
    if(status==nf90_noerr)EmisFile%sectorsName = trim(name)

!default values for sources
!species cannot be set global attribute, because it is used to recognize valid variables (sources)
!    status = nf90_get_att(ncFileID,nf90_global,"species",cdfspecies)
!    if(status==nf90_noerr)EmisFile%species = trim(cdfspecies)
    status = nf90_get_att(ncFileID,nf90_global,"units", name)
    if(status==nf90_noerr)EmisFile%units = trim(name)
    status = nf90_get_att(ncFileID,nf90_global,"sector", sector)
    if(status==nf90_noerr)EmisFile%sector = sector
    status = nf90_get_att(ncFileID,nf90_global,"country_ISO", name)
    if(status==nf90_noerr)EmisFile%country_ISO = trim(name)
    

    nemis_old = NEmis_sources
    !loop over all variables
    call check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes))
    do varid=1,nVariables
       call check(nf90_Inquire_Variable(ncFileID,varid,cdfvarname,xtype,ndims))
       status = nf90_get_att(ncFileID,varid,"species",cdfspecies)

       nn = 0
       do i = 1,size(EmisFile_in%source)       
          !if ( debugm0 ) write(*,*) dtxt//'source:',trim(EmisFile_in%source(i)%varname)
          if(EmisFile_in%source(i)%varname == cdfvarname)then
             nn = nn + 1
             call CheckStop(NEmis_sources+nn > NEmis_sourcesMAX,"NEmis_sourcesMAX exceeded (A)")
             Emis_source(NEmis_sources+nn)%ix_in=i
             if ( debugm0 ) write(*,*) dtxt//'var add:',trim(cdfvarname)
          endif
       enddo
       if((status==nf90_noerr .and. ndims>=2) .or. nn>0 )then
          !can be that one source must be taken several times (for instance
          ! into different vertical levels)
           do i = 1,max(1,nn)
             !we define a new emission source
             call CheckStop(NEmis_sources+1 > NEmis_sourcesMAX,"NEmis_sourcesMAX exceeded (B)")
             NEmis_sources = NEmis_sources + 1
             Emis_source(NEmis_sources)%varname = trim(cdfvarname)
             Emis_source(NEmis_sources)%species = trim(cdfspecies)
             if ( debugm0 ) write(*,*) dtxt//'source add:',&
              trim(cdfvarname)//'->'// trim(cdfspecies),EmisFile_in%apply_femis
             Emis_source(NEmis_sources)%units = EmisFile%units !default
             status = nf90_get_att(ncFileID,varid,"units", name)
             if(status==nf90_noerr)Emis_source(NEmis_sources)%units = trim(name)
             Emis_source(NEmis_sources)%sector = EmisFile%sector !default
             status = nf90_get_att(ncFileID,varid,"sector", sector)
             if(status==nf90_noerr)Emis_source(NEmis_sources)%sector = sector
             status = nf90_get_att(ncFileID,varid,"factor", x)
             if(status==nf90_noerr)Emis_source(NEmis_sources)%factor = x
             status = nf90_get_att(ncFileID,varid,"country", countrycode)
             if(status==nf90_noerr)then
                ix = find_index(countrycode, Country(:)%icode)
                if(ix<0)then
                   if(me==0)write(*,*)dtxt//'WARNING: country '//trim(name)//&
                     ' not defined. file'//trim(fname)//&
                     ' variable '//trim(cdfvarname)
                else
                   Emis_source(NEmis_sources)%country_ISO = trim(name)
                   Emis_source(NEmis_sources)%country_ix = ix
                   if ( debugm0 ) write(*,*) dtxt//'ISO add:',ix,trim(name)
                endif
             else
                Emis_source(NEmis_sources)%country_ISO = EmisFile%country_ISO !default
                status = nf90_get_att(ncFileID,varid,"country_ISO", name)
                if(status==nf90_noerr)Emis_source(NEmis_sources)%country_ISO = trim(name)
                ix = find_index(Emis_source(NEmis_sources)%country_ISO ,Country(:)%code, first_only=.true.)
                if(ix<0)then
                   if(me==0)write(*,*)dtxt//'WARNING: country_ISO '//trim(name)//&
                        ' not defined. file'//trim(fname)//&
                        ' variable '//trim(cdfvarname)
                else
                   Emis_source(NEmis_sources)%country_ix = ix
                   if ( debugm0 ) write(*,*) dtxt//'country_ISO add: ',ix,trim(name)
                endif
             endif
          enddo
       endif
        
    enddo
    if(nemis_old /= NEmis_sources)then
       !at least one valid source found in the file       
       NEmisFile_sources = NEmisFile_sources + 1
       EmisFile%filename = EmisFile_in%filename            
       EmisFile%projection = default_projection
       EmisFile%grid_resolution = default_resolution
       EmisFile%factor = default_factor       
       if ( debugm0 ) write(*,*) dtxt//'valid File:',&
           trim(EmisFile_in%filename),trim(default_projection), default_resolution
       status = nf90_get_att(ncFileID,nf90_global,"periodicity", name)
       if(status==nf90_noerr)EmisFile%periodicity = trim(name)

    endif
        
    call check(nf90_close(ncFileID))
    if ( debugm0 ) write(*,*) dtxt//'Finished File',trim(fname)
   !-------------
    end associate
   !-------------
    
  end subroutine Emis_init_GetCdf
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine EmisGetASCII(iem, fname, emisname, sumemis_local, incl, nin, excl, nex, &
                           	use_lonlat_femis)
    implicit none
    integer, intent(in) ::iem, nin, nex
    character(len=*),intent(in) :: fname, emisname, incl(*),excl(*)
    real,intent(inout), dimension(NLAND,NEMIS_FILE) &
        :: sumemis_local ! Sum of emissions per country
    logical, intent(in) :: use_lonlat_femis
    integer ::i,j,ic,CC,isec,i_gridemis,found
    character(len=*), parameter ::  sub = 'EmisGetASCII:'
    logical :: Cexist
    real :: tmpsec(NSECTORS),duml,dumh
    real ::lonlat_fac(NSECTORS)
    character(len=len(fname)+10) :: newfname !need some extra characters if the new name gets longer
    character(len=*),parameter :: dtxt = 'Em_getAsc:'

      call open_file(IO_EMIS, "r", fname, needed=.false., iostat=ios)

      if(ios/=0)then
         22 format(5A)
         if(MasterProc)write(*,22)dtxt//'did not find ',trim(fname)
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
      end if
READEMIS: do   ! ************* Loop over emislist files *******************

            read(unit=IO_EMIS,fmt=*,iostat=ios) CC,i,j, duml,dumh,  &
                                    (tmpsec(isec),isec=1,NSECTORS)
            if ( ios <  0 ) exit READEMIS            ! End of file
            call CheckStop(ios > 0,"EmisGetASCII: ios error in emission file")

            do ic=1,NLAND
               if((Country(ic)%icode==CC))&
                    goto 543
            end do 
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
             end if
             if(nex>0)then
                !1) check that country is not in exclude list
                found=find_index(Country(ic)%code ,excl(1:nex),first_only=.true.)
                if(found>0)cycle READEMIS!exclude
             end if

             lonlat_fac=1.0
             if(N_femis_lonlat>0 .and. use_lonlat_femis)then
                do i_femis_lonlat=1,N_femis_lonlat
                   if(glat(i,j)>femis_latmin(i_femis_lonlat).and.&
                        glat(i,j)<femis_latmax(i_femis_lonlat).and.&
                        glon(i,j)>femis_lonmin(i_femis_lonlat).and.&
                        glon(i,j)<femis_lonmax(i_femis_lonlat).and.&
                        (femis_lonlat_ic(i_femis_lonlat)==0 .or. &
                       femis_lonlat_ic(i_femis_lonlat)==CC) )then
                   !DSHK IN BOX
                       if ( femis_lonlat_internal(i_femis_lonlat) ) &
                         lonlat_fac(:)=lonlat_fac(:)*e_fact_lonlat(:,i_femis_lonlat,iem) 
                   !DSHK OUTSIDE BOX
                   else if ( femis_lonlat_internal(i_femis_lonlat) .eqv. .false. ) then
                     lonlat_fac(:)=lonlat_fac(:)*e_fact_lonlat(:,i_femis_lonlat,iem) 
                   end if
                end do
             end if

              !merge into existing emissions           
             Cexist=.false.
             do i_gridemis=1,nGridEmisCodes(i,j)
                if(GridEmisCodes(i,j,i_gridemis)==CC)then

                   GridEmis(:,i,j,i_gridemis,iem)=&
                         GridEmis(:,i,j,i_gridemis,iem)&
                        +e_fact(:,ic,iem) *  tmpsec(:)   *lonlat_fac(:)              

                   Cexist=.true.
                   exit
                end if
             end do
             if(.not.Cexist)then
                !country not included yet. define it now:
                nGridEmisCodes(i,j)=nGridEmisCodes(i,j)+1
                if(nGridEmisCodes(i,j)>NCMAX)then
                   write(*,*) sub// &
                      " Too many emitter countries in one gridemiscell: ",&
                        me,i,j,CC,&
                       (GridEmisCodes(i,j,i_gridemis),i_gridemis=1,NCMAX)
                   call StopAll("To many countries in one gridemiscell ")
                end if
                i_gridemis=nGridEmisCodes(i,j)
                GridEmisCodes(i,j,i_gridemis)=CC
                GridEmis(:,i,j,i_gridemis,iem)=e_fact(:,ic,iem) *  tmpsec(:)*lonlat_fac(:)
             end if
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
  integer            :: ix, ie, iq, ic, iland1, iland2 & ! loop variables
                       ,inland                     & ! Country read from femis
                       ,isec, isec1 , isec2        & ! loop vars: emis sectors
                       ,nwords,ncols, n, oldn       ! No. cols. in "femis" 
  integer, parameter        :: NCOLS_MAX = 20  ! Max. no. cols. in "femis"
  integer, dimension(NEMIS_FILE) :: qc        ! index for sorting femis columns
  real, dimension(NCOLS_MAX):: e_f             ! factors read from femis
  real, dimension(NCOLS_MAX):: e_f_lonlat      ! factors read from femis in lonlat format
  character(len=200) :: txt                    ! For read-in 
  character(len=30), dimension(NCOLS_MAX)::  txtinwords ! to read lines
  character(len=*), parameter :: dtxt = 'femis:'
  character(len=30) :: country_ISO, word30

 !--------------------------------------------------------


  e_fact(:,:,:) = 1.0            !**** default value = 1 ****
  e_fact_lonlat(:,:,:) = 1.0            !**** default value = 1 ****

  associate ( debugm0 => ( DEBUG%GETEMIS .and. MasterProc ) )

  if( debugm0 ) write(*,*) dtxt//" Enters femis", me, trim(femisFile)
  call open_file(IO_EMIS,"r",femisFile,needed=.false.)

  if ( ios == NO_FILE ) then
        ios = 0
        write( *,*) "WARNING: NO FEMIS FILE"
        return !*** if no femis file, e_fact=1 as default *** 
  end if
  call CheckStop( ios < 0 ,"EmisGet:ios error in "//trim(femisFile))


  ! Reads in the header line, e.g. name sec sox nox voc. 
  ! Pollutant names wil be checked against those defined in My_Emis_mod 

  read(unit=IO_EMIS,fmt="(a200)") txt
  !D if(debugm0)write(unit=6,fmt=*) "In femis, header0 is: ",  trim(txt)

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
  print *, "FEMIS me n NEMIS_FILE", me, n, NEMIS_FILE, trim(femisFile)
  call CheckStop( n < NEMIS_FILE , "EmisGet: too few femis items" )
end if

  
  n = 0
  N_femis_lonlat=0

  !if(me==0)write(IO_LOG,55)'femis reductions:'

  if(MasterProc) call PrintLog(dtxt//' reductions:')

  READFILE: do  ! ************ read lines of femis ***************
       read(unit=IO_EMIS,fmt="(a200)",iostat=ios) txt
       if( debugm0 ) write(*,*) dtxt//' input:'//trim(txt), ios
       if ( ios <  0 ) exit READFILE                   ! End of file
       call CheckStop( ios > 0 , "EmisGet: read error in femis" )

       if( txt(1:1) == '#' ) CYCLE ! Comments

       call wordsplit(txt,NCOLS_MAX,txtinwords,nwords,ios)
       if ( nwords<3 ) cycle READFILE                   ! End of file

       if(MasterProc) call PrintLog(txt)
 
       !lonlat box. reductions defined with coordinates
       ! xlonlat applies reductions outside (in testing)
       if(txtinwords(1)=='lonlat' .or.  txtinwords(1)=='xlonlat')then
          if(debugm0) write(*,*) dtxt//' LONLAT'//trim(txtinwords(1) ), &
               nwords, ncols
          if(nwords<ncols+6)then
             if(me==0)write(*,*)trim(femisFile)//' not understood ',nwords,ncols+5,txt
             call CheckStop( nwords<ncols+5 , "EmisGet: read error in femis lonlat" )
          end if
          !latmin,latmax,lonmin,lonmax
          N_femis_lonlat=N_femis_lonlat+1
          femis_lonlat_internal(N_femis_lonlat) = .true.
          call CheckStop( N_femis_lonlat>MAXFEMISLONLAT, "EmisGet: increase MAXFEMISLONLAT" )

         !DSHK apply reductions area external to box:
          if( txtinwords(1) == 'xlonlat' ) then
             femis_lonlat_internal(N_femis_lonlat) = .false.
             if( debugm0 ) write(*,*) dtxt//' exclude outside:'//trim(txt)
          else
             if( debugm0 ) write(*,*) dtxt//' exclude inside:'//trim(txt)
          end if

          inland = 0!default
          if(nwords>ncols+6)then
             !         76 format(A,4F,2I,20F)
             read(txt,*)txtinwords(1),&
                  femis_lonmin(N_femis_lonlat),&
                  femis_lonmax(N_femis_lonlat),&
                  femis_latmin(N_femis_lonlat),&
                  femis_latmax(N_femis_lonlat),&
                  inland, isec, & !country code and sector
                  (e_f_lonlat(ic),ic=1,ncols)!reductions
          else
             !old format without country code
             !         77 format(A,4F,I,20F)
             if(me==0)write(*,*)'WARNING: USING OLD FORMAT for femis.&
                  Please, add country code (zero for all countries)'
             read(txt,*)txtinwords(1),&
                  femis_lonmin(N_femis_lonlat),&
                  femis_lonmax(N_femis_lonlat),&
                  femis_latmin(N_femis_lonlat),&
                  femis_latmax(N_femis_lonlat)&
                  ,isec,(e_f_lonlat(ic),ic=1,ncols)
          endif

          femis_lonlat_ic(N_femis_lonlat) = inland !country code to reduce (or reduce all for 0)             


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

          if(txtinwords(1)=='Country' .or. txtinwords(1)=='country'  .or. txtinwords(1)=='Country_ISO' )then
             read(txt,fmt=*,iostat=ios) word30, country_ISO, isec, (e_f(ic),ic=1,ncols)
             ix = find_index(trim(country_ISO),Country(:)%code)!find country array index from ISO
             if(ix<0)then
                if(MasterProc)write(*,*)'femis: Country ',trim(country_ISO),' not recognized'
                CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
                CALL MPI_FINALIZE(IERROR)
                stop
             else
                if(MasterProc)write(*,*)'femis: reducing Country ',trim(Country(ix)%name)
             endif
             inland =  Country(ix)%icode
          else
             read(txt,fmt=*,iostat=ios) inland, isec, (e_f(ic),ic=1,ncols)
          endif
          n = n + 1
          if(debugm0) then
             write(unit=6,fmt='(a,2i5,99f9.4)') dtxt//"FEMIS READ", inland, &
                  isec, (e_f(ic),ic=1,ncols)
             write(unit=6,fmt="(2a,I3,a,i3,a)") &
                  dtxt//" Emission factors from femis.dat, ",&
                  "landcode =", inland, ",  sector code =",isec, &
                  " (sector 0 applies to all sectors) :"
             write(unit=6,fmt="(a,14(a,a,F8.2,a))") " ", &
                  (trim(txtinwords(qc(ie)+2)),&
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
               if(Country(iland1)%icode==inland  ) then
                 if (debugm0 .and. inland==6) write(*,*)'DKLAND',inland,iland1
                 goto 544
               end if
             end do
             
             if(MasterProc) write(*,*)'COUNTRY CODE NOT RECOGNIZED',inland
             
             iland1 = 0
             iland2 =-1
544          continue
             if(iland1/=0)   iland2 = iland1
          end if
       end if

       if (isec == 0 ) then       ! All sectors
          isec1 = 1
          isec2 = NSECTORS
       elseif (isec==100) then    ! Anthropogenic scenario
          !if you need this option, then just uncomment and check that
          !ANTROP_SECTORS is set according to your categories (10 for SNAP)
          call CheckStop(SECTORS_NAME/='SNAP',&
                "Anthropogenic not compatible with non-SNAP")
          isec1 = 1
          isec2 = ANTROP_SECTORS
       else                       ! one sector: isec
          isec1 = isec
          isec2 = isec
       end if
       if(debugm0) write(*,'(a,4i5)') dtxt//" SECS", inland, isec, isec1, isec2, trim(txtinwords(1))
       do ie = 1,NEMIS_FILE
          if(txtinwords(1)=='lonlat')then
             do isec = isec1, isec2
                e_fact_lonlat(isec,N_femis_lonlat,ie) = &
                e_fact_lonlat(isec,N_femis_lonlat,ie) * e_f_lonlat( qc(ie) )
             end do !isec
          else
             do iq = iland1, iland2
                do isec = isec1, isec2
                   e_fact(isec,iq,ie) = e_fact(isec,iq,ie) * e_f( qc(ie) )
                end do !isec
             end do !iq
             if (debugm0 ) then
                write(unit=6,fmt='(a,2i5,f8.3,a,4i4)') 'NEMIS_FILE LOOP'//&
                  ' HAS : ', ie, qc(ie), e_f( qc(ie) ), &
                  ' loops over  (secs,ccs)', isec1, isec2, iland1, iland2
             end if ! DEBUG%GETEMIS
          end if
       end do !ie
          
  end do READFILE ! Loop over femis

  close(IO_EMIS)

  if( MasterProc) write(unit=6,fmt=*) "femis, read ", n, "i j lines from femis"
  if( MasterProc) write(unit=6,fmt=*) "femis, read ", N_femis_lonlat, &
                                      "lonlat lines from femis"
  if(debugm0) then
    ! Extra checks
     write(unit=6,fmt=*) "DEBUG%GETEMIS: DK femis gives: "
     write(unit=6,fmt="(6x, 30a10)") (EMIS_FILE(ie), ie=1,NEMIS_FILE)
     do isec = 1, NSECTORS
      write(unit=6,fmt="(i6, 30f10.4)") isec, &
          (e_fact(isec,6,ie),ie=1,NEMIS_FILE)
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
   integer :: isec,k_ext,k1_ext(KMAX_BND),nemis_hprofile

   !emis_hprofile are read from file. 
   !emis_kprofile are the fraction values converted into model levels

   !REMARK: if you only change nemis_hprofile, but keep exactely the same
   !        emissions, the results will still be sligthly changed.
   !        This is because nemis_hprofile is used to define KEMISTOP, which 
   !        defines the levels where to use 2 or 3 chemical "2steps" iterations.

   !use old format
   call open_file(IO_EMIS,"r",EmisHeightsFile,needed=.true.)


   do
      call read_line(IO_EMIS,txtinput,ios,'EmisHeight')

      if(me==1) write(*,fmt='(A)') "read from "//trim(EmisHeightsFile)//": " // trim(txtinput)
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
         end if
         if( DEBUG%GETEMIS.and.MasterProc ) write(*,*) "VER=> ",snap, tmp(1), tmp(3)
         emis_hprofile(1:nemis_hprofile,snap) = tmp(1:nemis_hprofile)
      end if
   end do

   call CheckStop(nemis_hprofile < 1,"EmisGet: No EmisHeights set!!")
   call CheckStop( any( emis_hprofile(:,:) < 0 ), "EmisHeight read failure" )

   close(IO_EMIS)

   !Pressure boundaries for emission levels defined in EmisHeightsFile
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
      end do
      nemis_kprofile=nemis_hprofile
      allocate(emis_kprofile(nemis_kprofile,N_HFAC),stat=allocerr)
      emis_kprofile(1:nemis_kprofile,:)=emis_hprofile(1:nemis_hprofile,:)

   else

      if(Masterproc)then
         write(*,*)'emission heights: defined from pressure levels'
         do k=0,nemis_hprofile
            write(*,*)'P emis levels : ',k,emis_P_level(k)
         end do
      end if
      !stop
      !find highest level used
      nemis_kprofile = 0
      do k=KMAX_BND-1,1,-1
         nemis_kprofile = nemis_kprofile + 1
         if(A_bnd(k)+B_bnd(k)*Pref-0.0001<emis_P_level(nemis_hprofile))exit
      end do
      if(MasterProc) write(*,*)'Emissions distributed among ',nemis_kprofile,' lowest levels'

      allocate(emis_kprofile(nemis_kprofile,N_HFAC),stat=allocerr)
      emis_kprofile=0.0

      !convert height (given as pressure) distribution into model level distribution
      !ext meaning using levels from file (EmisHeightsFile)
      k1_ext(KMAX_BND)=0
      do isec=1,N_HFAC! could put inside but easier for debugging to put here now
         if(DEBUG%GETEMIS.and.MasterProc) write(*,*)'Sector ',isec
         do k=KMAX_BND-1,max(1,KMAX_BND-nemis_kprofile),-1

!count all contributions to model level k
!i.e. between P_emep(k+1) and P_emep(k)

            P_emep=A_bnd(k)+B_bnd(k)*Pref !Pa
            if(DEBUG%GETEMIS.and.MasterProc) write(*,fmt="(A,I3,F10.2)")'vert_inter: P_emep',k,P_emep
            !largest available P_ext smaller than P_emep (if possible)
            !k1_ext(k) is the external layer just below P_emep(k)
            k1_ext(k)= 0 !start at surface, and go up until P_emep
            do k_ext=1,nemis_hprofile
               if(emis_P_level(k_ext)<P_emep)exit
               k1_ext(k)=k_ext
            end do
            if(DEBUG%GETEMIS.and.MasterProc) write(*,*)k,k1_ext(k),k1_ext(k+1)

            !sum all contributions starting from last counted level (i.e. P_emep(k+1))

!part just above P_emep(k+1)
            if(emis_P_level(k1_ext(k+1)+1)>P_emep)then
               !part below k1_ext(k+1)+1  above P_emep(k+1)
               frac=((A_bnd(k+1)+B_bnd(k+1)*Pref )-emis_P_level(k1_ext(k+1)+1))&
                    /(emis_P_level(k1_ext(k+1))-emis_P_level(k1_ext(k+1)+1))
               emis_kprofile(KMAX_BND-k,isec)=frac*emis_hprofile(k1_ext(k+1)+1,isec)
               if(DEBUG%GETEMIS.and.MasterProc) &
                write(*,fmt="(A,I5,6F10.2)")'adding fraction of level',&
                  k1_ext(k+1)+1,frac,emis_hprofile(k1_ext(k+1)+1,isec),&
                  emis_P_level(k1_ext(k+1)+1),(A_bnd(k+1)+B_bnd(k+1)*Pref),&
                  emis_P_level(k1_ext(k+1)+1),(emis_P_level(k1_ext(k+1))-emis_P_level(k1_ext(k+1)+1))
            else
               !everything between P_emep(k+1) and P_emep(k)
               frac=((A_bnd(k+1)+B_bnd(k+1)*Pref )-P_emep)/(emis_P_level(k1_ext(k+1))-emis_P_level(k1_ext(k+1)+1))
               emis_kprofile(KMAX_BND-k,isec)=frac*emis_hprofile(k1_ext(k+1)+1,isec)
               if(DEBUG%GETEMIS.and.MasterProc) &
                write(*,"(A,I5,6F10.2)")'adding fraction of level between P_emep(k+1) and P_emep(k)',&
                    k1_ext(k+1)+1,frac,emis_hprofile(k1_ext(k+1)+1,isec),emis_P_level(k1_ext(k+1)+1),&
                    (A_bnd(k+1)+B_bnd(k+1)*Pref ),P_emep,(emis_P_level(k1_ext(k+1))-emis_P_level(k1_ext(k+1)+1))
            end if

            !add all full levels in between
            do k_ext=k1_ext(k+1)+2,k1_ext(k)
               emis_kprofile(KMAX_BND-k,isec)=emis_kprofile(KMAX_BND-k,isec)+emis_hprofile(k_ext,isec)
               if(DEBUG%GETEMIS.and.MasterProc) &
                 write(*,"(A,I5,6F10.2)")'adding entire level',&
                   k_ext,emis_hprofile(k_ext,isec),emis_P_level(k_ext)
            end do

            !add level just below P_emep(k), if not already counted, above k1_ext(k)  below P_emep; and must exist
            if(emis_P_level(k1_ext(k+1)+1)>P_emep.and. k1_ext(k)<nemis_hprofile )then
               frac=(emis_P_level(k1_ext(k))-P_emep)/(emis_P_level(k1_ext(k))-emis_P_level(k1_ext(k)+1))
               emis_kprofile(KMAX_BND-k,isec)=emis_kprofile(KMAX_BND-k,isec)+frac*emis_hprofile(k1_ext(k)+1,isec)
               if(DEBUG%GETEMIS.and.MasterProc) write(*,fmt="(A,I5,6F10.2)")'adding last fraction of level',k1_ext(k)+1,&
                    frac,emis_hprofile(k1_ext(k)+1,isec),emis_P_level(k1_ext(k+1)+1),emis_P_level(k1_ext(k)),P_emep,&
                    (emis_P_level(k1_ext(k))-emis_P_level(k1_ext(k)+1))
            end if

         end do
      end do
      
      CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)!so that output from different CPU does not get mixed up
      
      if(MasterProc)then
         write(*,*)'Distribution of emission into levels:'
         do isec=1,N_HFAC! could put inside but easier for debugging to put here now
            write(*,fmt="(A,I5,A,20F6.3)")'sector: ',isec,' fractions: ',(emis_kprofile(k,isec),k=1,nemis_kprofile)
         end do
      end if

   end if

!check normalization of distributions:
   do isec=1,N_HFAC
      sum=0.0
      do k=1,nemis_kprofile
         sum=sum+emis_kprofile(k,isec)
      end do
      if(abs(sum-1.0)>0.01)then
        if(MasterProc)then
           write(*,*)'WARNING emis height distribution not normalized : ',sum,(emis_kprofile(k,isec),k=1,nemis_kprofile)
           call StopAll( 'emis height distribution not normalized  ')
        end if
      end if
   end do

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
  integer, dimension(NSPEC_ADV) :: tmp_iqrc2iem !maps from iqrc 
  real, dimension(NMAX ) :: tmp 
  real     :: sumtmp
  integer  :: nsplit   &         ! No.columns data to be read
             ,iland,isec,i,n,nn, allocerr,iland_icode
  integer  :: idef              ! Set to   0  for defaults, 1 for specials
  integer  :: iland1, iland2    ! loop variables over countries
  logical  :: defaults          ! Set to true for defaults, false for specials
  logical  :: debugm            ! debug flag on master proc 
  character(len=*), parameter :: dtxt = 'EmisGet:'
  integer  :: itot_RDF
!-----------------------------------------------

  iqrc = 0               ! Starting index in emisfrac array
  nrcsplit= 0                 !
  itot2iqrc = -1 !init. Used to recognized non-split emissions (individual species)
  iemsplit2itot = -1 !init

  debugm = (DEBUG%GETEMIS.and.MasterProc) 

  do ie = 1, NEMIS_FILE 

    IDEF_LOOP: do idef = 0, 1

       defaults = (idef == 0)

       if ( defaults ) then

          fname = key2str(SplitDefaultFile,'POLL', EMIS_FILE(ie) )
          call open_file(IO_EMIS,"r",fname,needed=.true.)

          call CheckStop( ios, "EmisGet: ioserror:split.defaults " )

       else 
    !** If specials exists, they will overwrite the defaults

          fname = key2str(SplitSpecialsFile,'POLL', EMIS_FILE(ie) )
          call open_file(IO_EMIS,"r",fname,needed=.false.)

          if ( ios == NO_FILE ) then  
              ios = 0
              if(MasterProc) &
                 write(*,fmt=*) "emis_split: no specials for:",EMIS_FILE(ie)

              exit IDEF_LOOP
          end if
       end if


       if (debugm) write(*,*) "DEBUG%GETEMIS split defaults=", defaults,fname
 
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
          if(DEBUG%GETEMIS) then
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
               
                  ! This error, itot<1, is quite common when eg emissplit
                  ! files don't match chemistry. Add extensive output
                  if ( itot<1 ) then 
                    print *, "EmisSplit FAILED idef ", me, idef, i, nsplit,&
                      trim( intext(idef,i) )
                    print *, " Failed Splitting ", trim(EMIS_FILE(ie)), &
                      " emissions into ",&
                       (trim(Headers(n+2)),' ',n=1,nsplit),'using ',trim(fname)
                    print "(a, i3,30a10)", "EmisSplit FAILED headers ", &
                        me, (intext(idef,n),n=1,nsplit)
                    call StopAll( &
                       "EmisSplit FAILED "//trim(intext(idef,i)) //&
                       " possible incorrect Chem in run script?" )
                  end if ! FAILURE

                  tmp_iqrc2itot(iqrc) = itot
                  tmp_iqrc2iem(iqrc) = ie
                  itot2iqrc(itot)     = iqrc
                  iemsplit2itot(emis_nsplit(ie),ie) = itot
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
          end if

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
             do i = 1, emis_nsplit(ie)

                !*** assign and convert from percent to fractions: ***

                iqrc = sum(emis_nsplit(1:ie-1)) + i
                tmp_emisfrac(iqrc,isec,iland) = 0.01 * tmp(i)

                ! just a check
                !if ( DEBUG .and. iland == 27.and.MasterProc ) then 
                if ( DEBUG%GETEMIS .and. iland == 101.and.MasterProc ) then 
                    itot = tmp_iqrc2itot(iqrc)
                    write(*,"(a35,4i3,i4,a,f10.4)") &
                      "DEBUG%GETEMIS splitdef UK", isec, ie, i,  &
                       iqrc, itot, trim(species(itot)%name), &
                         tmp_emisfrac(iqrc,isec,iland)
                end if
             end do ! i
           end do ! iland

       end do READ_DATA 
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
    if ( DEBUG%GETEMIS .and. MasterProc .and. ( EmisSpecFound(ie) .eqv. .false.) ) then
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
  allocate(iqrc2iem(nrcemis),stat=allocerr)
  allocate(iqrc2itot(nrcemis),stat=allocerr)
  call CheckStop(allocerr, "Allocation error for iqrc2itot")
  allocate(emis_masscorr(nrcemis),stat=allocerr)
  call CheckStop(allocerr, "Allocation error for emis_masscorr")
  emisfrac(:,:,:)     = tmp_emisfrac(1:nrcemis,:,:)
  iqrc2itot(:)        = tmp_iqrc2itot(1:nrcemis)
  iqrc2iem(:)        = tmp_iqrc2iem(1:nrcemis)
  emis_masscorr(:)    = tmp_emis_masscorr(1:nrcemis)

!rb: not ideal place for this but used here for a start
! Temporary solution! Need to find the molweight from the
! GenChem input but just to get something running first set
! a hard coded molar mass of 200. 
  if(USES%ROADDUST)THEN
     allocate(roaddust_masscorr(NROADDUST),stat=allocerr)
     call CheckStop(allocerr, "Allocation error for emis_masscorr")
     itot_RDF = find_index( "Dust_ROAD_f", species(:)%name ,any_case=.true. )
     call CheckStop(itot_RDF<=0, "Asked for road dust but did not find Dust_ROAD_f")
     do ie=1,NROADDUST
        roaddust_masscorr(ie)=1.0/species(itot_RDF)%molwt
     end do
  end if
   
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
  integer :: i, j, iland,             &  ! loop variables
             iic,ic                      ! country code (read from file)
  real    :: tmpdust                     ! for reading road dust emission potential file
  integer, save :: ncmaxfound = 0        ! Max no. countries found in grid
  character(len=300) :: inputline

   !>============================

    if ( my_first_road ) then
         if(DEBUG%ROADDUST)WRITE(*,*)"initializing sumroaddust!"
         sumroaddust(:,:) =  0.0       ! initialize sums
         ios = 0
         my_first_road = .false.
    end if

  !>============================

      globroad_dust_pot(:,:,:) = 0.0

      if (DEBUG%ROADDUST) write(unit=6,fmt=*) "Called RoadDustGet with index, name", &
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
      end if     

READEMIS: do   ! ************* Loop over emislist files *******************

            read(unit=IO_EMIS,fmt=*,iostat=ios) iic,i,j, tmpdust

!            write(*,*)'dust to dust',iic,i,j, tmpdust

            if( DEBUG%ROADDUST .and. i==DEBUG%IJ(1) .and. j==DEBUG%IJ(2) ) write(*,*) &
                "DEBUG RoadDustGet "//trim(emisname) // ":" , iic, tmpdust
            if ( ios <  0 ) exit READEMIS            ! End of file
            call CheckStop(ios > 0,"RoadDustGet: ios error2 in emission file")

            ! Check if country code in emisfile (iic) is in the country list
            ! from Countries_mod, i.e. corresponds to numbering index ic

            do ic=1,NLAND
               if((Country(ic)%icode==iic))&
                    goto 654
            end do
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
            if( DEBUG%ROADDUST .and. i==DEBUG%IJ(1) .and. j==DEBUG%IJ(2) ) write(*,*) &
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
        
        if(DEBUG%ROADDUST) write(*,*)'done one, got sumroaddust=',sumroaddust(:,iemis)        

        close(IO_EMIS)
        ios = 0
  end subroutine RoadDustGet


subroutine make_iland_for_time(debug_tfac, indate, i, j, iland, wday, iland_timefac,hour_iland,wday_loc,iland_timefac_hour)
  ! make iland_timefac,hour_iland,wday_loc,iland_timefac_hour
  implicit none
  logical, intent(in):: debug_tfac
  type(date), intent(in):: indate
  integer, intent(in):: i,j,iland, wday
  integer, intent(out):: iland_timefac,hour_iland,wday_loc,iland_timefac_hour

  integer :: lon
  iland_timefac = find_index(Country(iland)%timefac_index,Country(:)%icode)
  iland_timefac_hour = find_index(Country(iland)%timefac_index_hourly,Country(:)%icode)
  if(Country(iland)%timezone==-100)then
     ! find the approximate local time:
     lon = modulo(360+nint(glon(i,j)),360)
     if(lon>180.0)lon=lon-360.0
     hour_iland= mod(nint(indate%hour+24*(lon/360.0)),24) + 1   ! add 1 to get 1..24 
  else
     hour_iland = indate%hour + Country(iland)%timezone + 1! add 1 to get 1..24 
  end if
  wday_loc=wday 
  if(hour_iland>24) then
     hour_iland = hour_iland - 24
     wday_loc=wday + 1
     if(wday_loc==0)wday_loc=7 ! Sunday -> 7
     if(wday_loc>7 )wday_loc=1 
  end if
  if(hour_iland<1) then
     hour_iland = hour_iland + 24
     wday_loc=wday - 1
     if(wday_loc<=0)wday_loc=7 ! Sunday -> 7
     if(wday_loc>7 )wday_loc=1 
  end if
  
  if(debug_tfac) then 
     write(*,"(a,2i3,i5,3x,4i3)") "EmisSet DAYS times ", &
          wday, wday_loc, iland,&
          hour_iland, Country(iland)%timezone
     call datewrite("EmisSet DAY 24x7:", &
          (/ i, iland, wday, wday_loc, hour_iland /), &
          (/ fac_ehh24x7(ISNAP_TRAF,hour_iland,wday_loc,iland_timefac_hour) /) )
  end if
    
end subroutine make_iland_for_time


end module EmisGet_mod
