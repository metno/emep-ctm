! <Emergency_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Emergency_ml
!-----------------------------------------------------------------------!
! Emissions for emergency (eEEMP) scenarios
!   Volcanic eruption: ASH and SO2 (optional)
!   NPP accident: Radiaton exposure tracers
! Any scenario can set emissions dates for other species as well.
! Varaible naming follows the Volcanic Eruption scenario, e.g.
!   vent : location of the emergency
!   erup : emission parameters
!-----------------------------------------------------------------------!
use CheckStop_ml,         only: CheckStop
use ChemSpecs,            only: NSPEC_TOT, NSPEC_SHL, species
!CMR use ChemChemicals_ml,     only: species
!CMR use ChemSpecs_shl_ml,     only: NSPEC_SHL
!CMR use ChemSpecs_tot_ml,     only: NSPEC_TOT
use ChemGroups_ml,        only: chemgroups
use GridValues_ml,        only: xm2,sigma_bnd,GridArea_m2,&
                                coord_in_processor,coord_in_gridbox
use Io_ml,                only: open_file,read_line,IO_TMP
use SmallUtils_ml,        only: wordsplit,find_index
use ModelConstants_ml,    only: KCHEMTOP,KMAX_MID,MasterProc, &
                                USE_EMERGENCY,DEBUG=>DEBUG_EMERGENCY,&
                                TXTLEN_NAME,dt_advec,dt_advec_inv
use OwnDataTypes_ml,      only: TXTLEN_SHORT
use MetFields_ml,         only: roa, z_bnd
use Par_ml,               only: me
use PhysicalConstants_ml, only: AVOG
use TimeDate_ml,          only: nydays,&           ! No. days per year
                                startdate,enddate,current_date,&
                                make_timestamp,tdif_secs
use TimeDate_ExtraUtil_ml,only: date2string,string2date

implicit none
private

 !/* subroutines:
public :: EmergencyRate      ! Emission rate

logical, save ::      &
  Emergency_found=.true.     ! Are sources found on this processor/subdomain?

integer, parameter :: &
  NMAX_VENT = 24, & ! Max number of locations on processor/subdomain
  NMAX_ERUP =360    ! Max number of events def per location (~3 months of 6 hr.records)

integer, save ::   & ! No. of ... found on processor/subdomain
  nvent              = -1,  & ! Emergency locations
  nerup(0:NMAX_VENT) = -1     ! Emergency events per location

type :: vent
  character(len=9)              :: id  =''    ! e.g. V1702A02B
  character(len=TXTLEN_NAME)    :: name=''    ! e.g. Eyjafjöll
  real :: lat=-1.0,lon=-1.0,elev=-1.0         ! vent coords and elevation
  character(len=4)              :: etype=''   ! e.g. S0, 10kt
  integer :: grp=-1,iloc=-1,jloc=-1           ! Which (ash)goup,local i,j indes
 !character(len=TXTLEN_SHORT)   :: location,vtype  ! other info
endtype vent
type(vent), save, dimension(NMAX_VENT):: &
  ventdef=vent('UNDEF','UNKNOWN',-999.0,-999.0,-999.0,"??",-99,-99,-99)

character(len=*), parameter :: SDATE_FMT="YYYY-MM-DD hh:mm:ss"
type :: erup
  character(len=9)              :: id  =''    ! e.g. V1702A02B
  character(len=TXTLEN_NAME)    :: name=''    ! e.g. SO2
  real :: base=-1.0,top=-1.0,&                ! Column Base & Height [m asl]
          rate=-1.0                           ! Source strenght: Total release for period [kg/s]
  character(len=len(SDATE_FMT)) :: sbeg=SDATE_FMT,send=SDATE_FMT
  integer :: vent=-1,spc=-1                   ! Which vent,(adv)spc
  logical :: edef=.true.                      ! default setings?
endtype erup
type(erup), save, dimension(0:NMAX_VENT,NMAX_ERUP):: &
  erupdef=erup('UNDEF','UNKNOWN',-999.0,-999.0,-999.0,"??","??",-1,-1,.true.)

character(len=*),parameter :: &
  fventdef="emergency_location.csv",   &   ! see ventdef
  ferupdef="emergency_emission.csv",   &   ! see erupdef
  ERR_VENT_CSV="Emergency VENT def.file "//fventdef,          &
  ERR_ERUP_CSV="Emergency ERUP def.file "//ferupdef,          &
  ERR_VENT_MAX="Emergency NMAX_VENT exceeded in "//fventdef,  &
  ERR_ERUP_MAX="Emergency NMAX_ERUP exceeded in "//ferupdef,  &
  MSG_FMT="('EMERGENCY:',:,1X,A,5(:,1X,I0,':',A),3(:,1X,G10.3,':',A))"

INCLUDE 'mpif.h'
INTEGER STATUS(MPI_STATUS_SIZE),INFO

contains
!-----------------------------------------------------------------------!
! Emergency scenarios:
!  Volcanic Eruption: Ash & SO2.
!  NPP Accident: Radiaton exposure tracers.
!  Any scenario can set emissions for ther species as well.
!-----------------------------------------------------------------------!
function EmergencyRate(i,j) result(emiss)
!-----------------------------------------------------------------------!
! Emissions from Emergency senarios.
!-----------------------------------------------------------------------!
  integer, intent(in) :: i,j
  real, dimension(NSPEC_SHL+1:NSPEC_TOT,KCHEMTOP:KMAX_MID) :: emiss
  logical, save :: first_call=.true.
  character(len=len(SDATE_FMT)) :: & ! Time strings in SDATE_FMT format
    sbeg=SDATE_FMT,&    ! Begin
    snow=SDATE_FMT,&    ! Now (current date)
    send=SDATE_FMT      ! End
  integer :: v,e,itot,k1,k0
  real    :: uconv
!----------------------------!
!
!----------------------------!
  if(first_call) call setRate()
  first_call=.false.
  emiss(:,:)=0.0
  if(.not.Emergency_found)then
    USE_EMERGENCY=.false.
    return
  endif
  snow=date2string(SDATE_FMT,current_date)
!----------------------------!
!
!----------------------------!
  doVENT: do v=1,nvent
!*! if(.not.coord_in_gridbox(ventdef(v)%lon,ventdef(v)%lat,i,j))&
!*!               cycle doVENT    ! Wrong gridbox
    if((i/=ventdef(v)%iloc).or.(j/=ventdef(v)%jloc))&
                  cycle doVENT    ! Wrong gridbox
    if(nerup(v)<1)cycle doVENT    ! Not erupting
    if(DEBUG) &
      write(*,MSG_FMT)snow//' Vent',me,'me',v,trim(ventdef(v)%id),i,"i",j,"j"
    doERUP: do e=1,nerup(v)
      sbeg=date2string(erupdef(v,e)%sbeg,current_date)
      send=date2string(erupdef(v,e)%send,current_date)
      if(snow<sbeg.or.send<snow)& ! Outside time window
        cycle doERUP
      itot=erupdef(v,e)%spc
      k0=getModLev(i,j,erupdef(v,e)%base)
      k1=getModLev(i,j,erupdef(v,e)%top)
!   uconv=1e6/(tdif_secs(make_timestamp(string2date(sbeg,SDATE_FMT)),&
!                        make_timestamp(string2date(send,SDATE_FMT))+1) ! Tg --> 10^6 g/s
      uconv=1e-3                                                        ! Kg/s --> 10^6 g/s
      uconv=uconv/(GridArea_m2(i,j)*DIM(z_bnd(i,j,k1),z_bnd(i,j,k0+1))) ! --> g/s/cm3
      uconv=uconv*AVOG/species(itot)%molwt                              ! --> molecules/s/cm3
      emiss(itot,k1:k0)=emiss(itot,k1:k0)+erupdef(v,e)%rate*uconv
      if(DEBUG) &
        write(*,MSG_FMT)snow//' Erup.',me,'me',e,erupdef(v,e)%sbeg,&
          itot,trim(species(itot)%name),k1,'k1',k0,'k0',&
          emiss(itot,k1),'emiss',erupdef(v,e)%rate,'rate',uconv,'uconv'
    enddo doERUP
  enddo doVENT
!----------------------------!
 contains
!----------------------------!
! Model level for a given height
!----------------------------!
  function getModLev(i,j,height) result(k)
    integer, intent(in) :: i,j
    real, intent(in) :: height
    integer :: k
    k=KMAX_MID
    if(height<=0.0)return
    do while(k>0.and.height>z_bnd(i,j,k))
      k=k-1
    enddo
  endfunction getModLev
!----------------------------!
! Set Volcanic Eruption Param.
!----------------------------!
  subroutine setRate()
    logical, save :: second_call=.true.
    character(len=80) :: txtline       ! Long enough for a full line
    type(vent)        :: dvent
    type(erup)        :: derup
    integer :: stat,l,v,e,g
! Particles classes & default split as London VAAC setup for NAME
! from Witham&al:2011 Table 1
!   Evolution of volcanic ash modelling at the London VAAC April 2010 --- April 2011
!   UK Met Office. Technical Summary (v1.0). May 2011.
!   C. Witham, M. Hort, D. Thomson, S. Leadbetter, B. Devenish and H. Webster.
    real, target :: & !0.0<0.1<0.3<1.0<3.0<10.0<30.0<100.0
      VAAC_7BIN_SPLIT(7)=(/0.0,0.1,0.5,5.0,20.0,70.0,4.4/)*1e-2,&
      VAAC_2BIN_SPLIT(2)=(/            5.6,20.0         /)*1e-2
    real, pointer, dimension(:) :: binsplit => NULL()
  !----------------------------!
  !
  !----------------------------!
    if(.not.USE_EMERGENCY)then
      Emergency_found=.false.
      if(MasterProc.and.DEBUG.and.first_call) &
        write(*,MSG_FMT)'Code is off. No Volcanic Ash.'
      first_call=.false.
      return
    endif
    if(.not.first_call)then
      if(MasterProc.and.DEBUG.and.second_call) &
        write(*,MSG_FMT)'No need for reset volc.def...'
      second_call=.false.
      return
    endif
    first_call=.false.
  !----------------------------!
  ! Read Vent CVS
  !----------------------------!
    if(DEBUG) CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    if(MasterProc)then
      call open_file(IO_TMP,"r",fventdef,needed=.true.,iostat=stat)
      call CheckStop(stat,ERR_VENT_CSV//' not found')
    endif
    nvent=0
    l = 1
    doVENT: do while (l<=NMAX_VENT)
      call read_line(IO_TMP,txtline,stat)
      if(stat/=0) exit doVENT           ! End of file
      txtline=ADJUSTL(txtline)          ! Remove leading spaces
      if(txtline(1:1)=='#')cycle doVENT ! Comment line
      dvent=getVent(txtline)
      if(coord_in_processor(dvent%lon,dvent%lat,iloc=dvent%iloc,jloc=dvent%jloc))then
        nvent=nvent+1
        call CheckStop(nvent>NMAX_VENT,ERR_VENT_MAX//" read")
        ventdef(nvent)=dvent
        if(DEBUG) &
          write(*,MSG_FMT)'Vent',me,'in',nvent,trim(dvent%id),&
            dvent%grp,trim(dvent%name),dvent%iloc,"i",dvent%jloc,"j",&
            dvent%lon,"lon",dvent%lat,"lat"
      elseif(MasterProc)then
        if(DEBUG) &
          write(*,MSG_FMT)'Vent',me,'out',-1,trim(dvent%id),&
            dvent%grp,trim(dvent%name),dvent%iloc,"i",dvent%jloc,"j",&
            dvent%lon,"lon",dvent%lat,"lat"
      endif
      l = l+1
    enddo doVENT
    if(MasterProc) close(IO_TMP)
    Emergency_found=(nvent>0)
  !----------------------------!
  ! Read Eruption CVS
  !----------------------------!
    if(DEBUG) CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    if(MasterProc)then
      call open_file(IO_TMP,"r",ferupdef,needed=.true.,iostat=stat)
      call CheckStop(stat,ERR_ERUP_CSV//' not found')
    endif
    nerup(:)=0
    l = 1
    doERUP: do while (l <= NMAX_ERUP)
      call read_line(IO_TMP,txtline,stat)
      if(stat/=0) exit doERUP             ! End of file
      if(.not.Emergency_found)cycle doERUP ! There is no vents on subdomain
      txtline=ADJUSTL(txtline)            ! Remove leading spaces
      if(txtline(1:1)=='#')cycle doERUP   ! Comment line
      derup=getErup(txtline)
      if(derup%edef)then                                  ! Default
        nerup(0)=nerup(0)+1
        call CheckStop(nerup(0)>NMAX_ERUP,ERR_ERUP_MAX//" read")
        erupdef(0,nerup(derup%vent))=derup
        if(DEBUG) &
          write(*,MSG_FMT)'Erup.Default',me,'in',nerup(0),trim(derup%id)
      elseif(derup%vent>0.and.derup%spc>0)then            ! Specific
        nerup(derup%vent)=nerup(derup%vent)+1
        call CheckStop(nerup(derup%vent)>NMAX_ERUP,ERR_ERUP_MAX//" read")
        erupdef(derup%vent,nerup(derup%vent))=derup
        if(DEBUG) &
          write(*,MSG_FMT)'Erup.Specific',me,'in',nerup(derup%vent),trim(derup%id),&
            derup%spc,trim(derup%name)
      elseif(MasterProc)then
        if(DEBUG) &                                    ! or Unknown Vent/SPC
          write(*,MSG_FMT)'Erup.Specific',me,'out',-1,trim(derup%id),&
            derup%spc,trim(derup%name)
      endif
      l = l+1
    enddo doERUP
    if(MasterProc) close(IO_TMP)
    Emergency_found=any(nerup(1:nvent)>0)
  !----------------------------!
  ! Expand Eruption Defaults
  !----------------------------!
    if(DEBUG) CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    if(nerup(0)<1)then
     !if(DEBUG) write(*,MSG_FMT)'Erup.Default',me,'not found'
      return
    endif
    doVENTe: do v=1,nvent
      if(nerup(v)>0)cycle doVENTe   ! Specific found --> no need for Default
      e=nerup(0)+1       ! A single defaul can have multiple lines, e.g.
      do                 ! each line with a difinition for a different specie
        e=find_index(ventdef(v)%etype,erupdef(0,:e-1)%id)
        if(e<1)       cycle doVENTe   ! No Default found
        if(DEBUG) &
          write(*,MSG_FMT)'Erup.Default',me,'Expand',&
            v,trim(ventdef(v)%id),e,trim(erupdef(0,e)%id)
        derup=erupdef(0,e)
        if(derup%spc<1.and.any(derup%name(1:3)==(/"ASH","NUC","###"/)))then ! Expand variable name
          derup%name=trim(ventdef(v)%id)//trim(derup%name(4:))  ! e.g. ASH_F --> V1702A02B_F
          derup%spc=find_index(derup%name,species(:)%name)      ! Specie (total)
          if(DEBUG)&
            write(*,MSG_FMT)'Erup.Default',me,'Expand',&
              derup%spc,trim(derup%name)
        endif
        if(derup%spc>0)then           ! Expand single SPC
          nerup(v)=nerup(v)+1
          call CheckStop(nerup(v)>NMAX_ERUP,ERR_ERUP_MAX//" expand")
          erupdef(v,nerup(v))=derup
          if(DEBUG) &
            write(*,MSG_FMT)'Erup.Default',me,'Expand SPC',nerup(v),trim(derup%id)
        elseif(ventdef(v)%grp>0)then   ! Expand GROUP of SPC
          select case (size(chemgroups(ventdef(v)%grp)%ptr))
          case(2);binsplit=>VAAC_2BIN_SPLIT
          case(7);binsplit=>VAAC_7BIN_SPLIT
          case default
            call CheckStop(ERR_ERUP_CSV//' can not expand '//trim(ventdef(v)%id))
          endselect
          do g=1,size(chemgroups(ventdef(v)%grp)%ptr)
            derup%spc=chemgroups(ventdef(v)%grp)%ptr(g)  ! Specie (total)
            derup%name=species(derup%spc)%name
            derup%rate=erupdef(0,e)%rate*binsplit(g)
            nerup(v)=nerup(v)+1
            call CheckStop(nerup(v)>NMAX_ERUP,ERR_ERUP_MAX//" expand")
            erupdef(v,nerup(v))=derup
            if(DEBUG) &
            write(*,MSG_FMT)'Erup.Default',me,'Expand GRP',nerup(v),trim(derup%id),&
              derup%spc,trim(derup%name)
          enddo
        else
          if(DEBUG) &
            write(*,MSG_FMT)'Erup.Default',me,'not found'
        endif
      enddo
    enddo doVENTe
    Emergency_found=any(nerup(1:nvent)>0)
  endsubroutine setRate
!----------------------------!
! Extract Vent info from CVS line
!----------------------------!
  function getVent(line) result(def)
    character(len=*)            :: line
    type(vent)                  :: def
    character(len=TXTLEN_SHORT) :: words(10)=''  ! Array of paramaters
    real    :: lat,lon,elev                      ! vent coord and elevation
    integer :: stat,nwords,igrp
    call wordsplit(line,size(words),words,nwords,stat,strict_separator=",",empty_words=.true.)
    call CheckStop(stat,"EMERGENCY: Wrong/Unknown line format "//trim(line))
    call CheckStop(nwords,size(words),"EMERGENCY: Missing data in line "//trim(line))
!#1:NUMBER,2:NAME,3:LOCATION,4:LATITUDE,5:NS,6:LONGITUDE,7:EW,8:ELEV,9:TYPE,10:ERUPTION TYPE
!V1702A02B,Eyjafjöll,Iceland-S,63.63,N,19.62,W,1666,Stratovolcano,S0
    read(words(4),*)lat
    select case (words(5)) ! NS
    case("N","n","degN")            ! degN
    case("S","s","degS");lat=-lat   ! degS
    case default
      call CheckStop("EMERGENCY: Unknown degN/S "//trim(words(5)))
    endselect
    read(words(6),*)lon
    select case (words(7)) ! EW
    case("E","e","degE")            ! degE
    case("W","w","degW");lon=-lon   ! degW
    case default
      call CheckStop("EMERGENCY: Unknown degE/W "//trim(words(7)))
    endselect
    read(words(8),*)elev
    igrp=find_index(words(1),chemgroups(:)%name)
    def=vent(trim(words(1)),trim(words(2)),lat,lon,elev,trim(words(10)),igrp)
  endfunction getVent
!----------------------------!
! Extract Erup. info from CVS line
!----------------------------!
  function getErup(line) result(def)
    character(len=*)            :: line
    type(erup)                  :: def
    character(len=TXTLEN_SHORT) :: words(10)=''  ! Array of paramaters
    logical :: edef=.true.                       ! default setings?
    integer :: stat,nwords,ivent,ispc=0
    real    :: base,top,rate,frac,dhh
    call wordsplit(line,size(words),words,nwords,stat,strict_separator=",",empty_words=.true.)
    call CheckStop(stat,"EMERGENCY: Wrong/Unknown line format "//trim(line))
    call CheckStop(nwords,size(words),"EMERGENCY: Missing data in line "//trim(line))
!#1:TYPE/VOLCANO,2:VARIABLE,3:BASE[km],4:H[km above vent],5:D[h],6:dM/dt[kg/s],7:m63[-],8:START[code/date],9:END[code/date],10:DESCRIPTION
!S0       ,     ,  , 11.000,   3.00, 4e6, 0.40,SR                 ,SR+D,Silicic standard
!V1702A02B,SO2  , 0,  8.000,  24.00,  15,     ,2010-04-14 00:00:00,SE+D,Eyja 20100414 SO2
!V1702A02B,ASH_F, 0,  2.000,  24.00,   0,     ,2010-05-23 00:00:00,SE+D,Eyja 20100523 PM fine
    ivent=find_index(words(1),ventdef(:nvent)%id)             ! Vent Specific
    edef=(ivent<1).and.any(ventdef(:nvent)%etype==words(1))   ! Vent Default
    if(ivent>0.and.words(2)(1:3)=="ASH")&         ! Expand variable name
      words(2)=trim(words(1))//trim(words(2)(4:)) ! e.g. ASH_F --> V1702A02B_F
    ispc=find_index(words(2),species(:)%name)     ! Specie (total)
    read(words(4),*)top    ! [km]
    top=top*1e3            ! [m]
    select case (words(3)) ! base
    case("VENT"," ")       ! From the vent
      base=0.0
      if(ivent>0)base=ventdef(ivent)%elev ! [m]
      top=top+base         ! [m]
    case("SURF","0")       ! From the model surface
      base=0.0
    case default
      read(words(3),*)base  ! [km]
      base=base*1e3         ! [m]
    endselect
    read(words(6),*)rate
    select case (words(7))  ! m63 or effect.fraction
      case(" ")   ;frac=1.0
      case default;read(words(7),*)frac
    endselect
    select case (words(5))  ! dt[h]
    case("1dt")
      dhh=(dt_advec-0.1)/3600.0 ! only one time step
      frac=frac*dt_advec_inv    ! assume rate=total emission
    case default
      read(words(5),*)dhh
    endselect   
    words(8)=getDate(words(8),words(8),words(9),dhh,debug=DEBUG) ! Start [date/code]
    words(9)=getDate(words(9),words(8),words(9),dhh,debug=DEBUG) ! End   [date/code]
    def=erup(trim(words(1)),trim(words(2)),base,top,rate*frac,&
      trim(words(8)),trim(words(9)),max(ivent,0),max(ispc,0),edef)
  endfunction getErup
!----------------------------!
! Time/Date CODE--> YYYY-MM-DD hh:mm:ss
!----------------------------!
  function getDate(code,se,ee,dh,debug) result(str)
    character(len=*), intent(in) :: code,se,ee
    real, intent(in)             :: dh  ! [hours]
    logical, intent(in), optional:: debug
    character(len=TXTLEN_SHORT)  :: str
    logical :: mydebug=.false.
    mydebug=.false.;if(present(debug))mydebug=debug.and..false.
    select case (code)
    case("SR")          ! Start of the simulation (model actually starts on 2nd time step)
      str=date2string(SDATE_FMT,startdate,addsecs=dt_advec,debug=mydebug)
    case("SR+D")        ! Start of the simulation + dh
      str=date2string(SDATE_FMT,startdate,addsecs=dh*36e2+dt_advec,debug=mydebug)
    case("SE+D")        ! Start eruption + dh; no wildcards in SE allowed
      str=date2string(SDATE_FMT,string2date(se,SDATE_FMT,debug=mydebug),addsecs=dh*36e2,debug=mydebug)
    case("EE+D")        ! End eruption   - dh; no wildcards in EE allowed
      str=date2string(SDATE_FMT,string2date(ee,SDATE_FMT,debug=mydebug),addsecs=-dh*36e2,debug=mydebug)
    case("ER-D")        ! End of the simulation - dh
      str=date2string(SDATE_FMT,enddate,addsecs=-dh*36e2,debug=mydebug)
    case("ER")          ! End of the simulation
      str=date2string(SDATE_FMT,enddate,debug=mydebug)
    case default
      str=code
    endselect
  endfunction getDate
endfunction EmergencyRate

endmodule Emergency_ml
