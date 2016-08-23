! <Volcanos_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Volcanos_ml
!-----------------------------------------------------------------------!
! Processes SOx emission heights from volcanoes
! In volcanos_ml, only the height of the volcanos are read from volcanos.dat,
! the emissions themselves comes from gridSOx, so this is to ensure that if
! suddenly e.g. Iceland starts to report volcano emissions and these are in
! gridSOx, the program will discover this.
! Note that nvolc is set according to the emission input from the gridSOx
! files.  It counts the number of grids with the "country code" for volcanoes
! that are in the gridSOx file.
!
! Note - EmisGet and other routines use "restricted" coords, which introduces
! some complications here. Hopefully we can tidy up one day.
!-----------------------------------------------------------------------!

use CheckStop_ml,         only: CheckStop
use ChemChemicals_ml,     only: species
use ChemSpecs_shl_ml,     only: NSPEC_SHL
use ChemSpecs_tot_ml,     only: NSPEC_TOT, SO2
use ChemGroups_ml,        only: chemgroups
use EmisDef_ml,           only: VOLCANOES_LL
use GridValues_ml,        only: GRIDWIDTH_M, xm2, sigma_bnd,  &
                                i_local, j_local, lb2ij, &
                                GridArea_m2,coord_in_processor,coord_in_gridbox
use Io_ml,                only: ios, NO_FILE, open_file,      &
                                IO_VOLC, Read_Headers, read_line,IO_TMP
use SmallUtils_ml,        only: wordsplit,find_index
use ModelConstants_ml,    only: KCHEMTOP, KMAX_BND, KMAX_MID, PT, MasterProc, &
                                DEBUG=>DEBUG_VOLC,&
                                USE_EMERGENCY,DEBUG_EMERGENCY,& ! Emergency: Volcanic Eruption
                                TXTLEN_NAME
use OwnDataTypes_ml,      only: TXTLEN_SHORT
use MetFields_ml,         only: roa, ps, z_bnd
use Par_ml,               only: me, IRUNBEG, JRUNBEG, ljmax, limax, &
                                gi0, gi1, gj0, gj1 ! Test if on correct processor
use PhysicalConstants_ml, only: GRAV, AVOG
use TimeDate_ml,          only: nydays,&           ! No. days per year
                                startdate,enddate,current_date,&
                                make_timestamp,tdif_secs
use TimeDate_ExtraUtil_ml,only: date2string,string2date
use KeyValue_ml,          only: KeyVal

implicit none
private

 !/* subroutines:
public :: VolcGet
public :: Set_Volc
public :: Scale_Volc
public :: EruptionRate ! Emergency: Volcanic Eruption

integer, public, parameter :: &
  NMAX_VOLC = 12  ! Max number of volcanoes
integer, public, save ::      &
  nvolc = 0,                  & ! No. grids with volcano emissions in gridSOx
  volc_no = 0

logical, public, save ::      &
  Volcanoes_found=.false.       ! Are Volcanoes found on this processor/subdomain?

integer, public, save, dimension(NMAX_VOLC):: &
  height_volc,                & ! Height of volcanoes
  i_volc, j_volc                ! Volcano's EMEP coordinates

real, private, save, dimension(NMAX_VOLC)::   &
  rcemis_volc0                  ! Emissions part varying every hour

real, public, save, dimension(NMAX_VOLC) ::   &
  rcemis_volc,                & ! Emissions part varying every time-step
  emis_volc = 0.0               ! Volcanoes' emissions

! Emergency: Volcanic Eruption
integer, public, parameter :: &
  NMAX_VENT = 10, & ! Max number of Vents (erupting) on processor/subdomain
  NMAX_ERUP = 90    ! Max number of Eruption def (~3 months per vent)

integer, private, save ::   & ! No. of ... found on processor/subdomain
  nvent              = -1,  & ! Vents
  nerup(0:NMAX_VENT) = -1     ! Eruption events per Vent

logical, parameter :: &
  DEBUG_EM=DEBUG.or.DEBUG_EMERGENCY

logical, private, save ::      &
  Vent_found=USE_EMERGENCY      ! Any Vent found on this processor/subdomain?

logical, public, save ::      &
  Eruption_found=USE_EMERGENCY  ! Any Eruption found on this processor/subdomain?

type, private :: vent
  character(len=9)              :: id  =''    ! e.g. V1702A02B
  character(len=TXTLEN_NAME)    :: name=''    ! e.g. Eyjafjöll
  real :: lat=-1.0,lon=-1.0,elev=-1.0         ! vent coords and elevation
  character(len=2)              :: etype=''   ! e.g. S0
  integer :: grp=-1                           ! Which (ash)goup ...
 !character(len=TXTLEN_SHORT)   :: location,vtype  ! other info
endtype vent
type(vent), private, save, dimension(NMAX_VENT):: &
  ventdef=vent('UNDEF','UNKNOWN',-999.0,-999.0,-999.0,"??",-99)

character(len=*), private , parameter :: SDATE_FMT="YYYY-MM-DD hh:mm:ss"
type, private :: erup
  character(len=9)              :: id  =''    ! e.g. V1702A02B
  character(len=TXTLEN_NAME)    :: name=''    ! e.g. SO2
  real :: base=-1.0,top=-1.0,&                ! Column Base & Height [m asl]
          rate=-1.0                           ! Source strenght: Total release for period [kg/s]
  character(len=len(SDATE_FMT)) :: sbeg=SDATE_FMT,send=SDATE_FMT
  integer :: vent=-1,spc=-1                   ! Which vent,(adv)spc...
  logical :: edef=.true.                      ! default setings?
endtype erup
type(erup), private, save, dimension(0:NMAX_VENT,NMAX_ERUP):: &
  erupdef=erup('UNDEF','UNKNOWN',-999.0,-999.0,-999.0,"??","??",.true.)

INCLUDE 'mpif.h'
INTEGER STATUS(MPI_STATUS_SIZE),INFO

contains

subroutine VolcGet(height_volc)
!-----------------------------------------------------------------------!
! Reads volcanoes' coordinates (i,j) and height level(height)
! Returns to EmisSet with height of volcanos (height_volc)
! Input file: Volcanoes.dat
!-----------------------------------------------------------------------!
  integer, intent(out), dimension(NMAX_VOLC) :: height_volc
  integer            :: nvolc_read,height,i,j          ! Local variables
  character (len=23) :: fname

  integer            :: nin,i_l,j_l
  real               :: lon,lat,xr,yr,conv,tonne_to_kgm2s

  character(len=70)               :: errmsg ! Message text
  character(len=20), dimension(5) :: Headers
  type(KeyVal), dimension(20)     :: KeyValues ! Info on units, coords, etc.
  integer                         :: NHeaders, NKeys
  character(len=80)               :: txtinput,s  ! Big enough to contain
                                                ! one full input record
  ios=0    ! Start with  assumed ok status

  height_volc(:)=KMAX_MID
  nvolc_read=0

  if(.not.VOLCANOES_LL)then ! Data read from gridSOx and Volcanoes.dat
    fname = "Volcanoes.dat"
    if (MasterProc)then
      call open_file(IO_VOLC,"r",fname,needed=.true.,skip=1)

      if(ios/=NO_FILE)then  !"Volcanoes.dat" does exist
        call CheckStop(ios,"VolcGet: problems with Volcanoes.dat ")

        READVOLC: do
          read(IO_VOLC,*,iostat=ios) i,j,height
          if (DEBUG) print "(A,2I5,F12.3,I0)", &
            'VOLCIJ:found i,j,height,ios',i,j,height,ios
          if ( ios /= 0 ) exit READVOLC

          !/** Set the volcano number to be the same as in emission data (gridSOx)
          do volc_no=1,nvolc
            if ((i_volc(volc_no)==i) .and. (j_volc(volc_no)==j)) then
              height_volc(volc_no)=height
              nvolc_read=nvolc_read+1
              if (DEBUG) print "(A,I0)",&
                'VOLCIJ:Found volcano with height k=',height
            endif
          enddo
        enddo READVOLC

        print *,nvolc_read,' volcanos on volcanos.dat&
              & match volcanos on emislist.sox'
        print *,nvolc,' volcanos found in emislist.sox'

        call CheckStop(nvolc_read < nvolc, "Volc missing in Volcanos.dat")
      endif
      close(IO_VOLC)
    endif

    CALL MPI_BCAST(ios,1,MPI_INTEGER,0,MPI_COMM_WORLD,INFO)

  else ! use "VolcanoesLL.dat" and disregard volcanoes from gridSOx
    nvolc=0
    !read volcanoes from lon lat format file:
    fname = "VolcanoesLL.dat"
    if (MasterProc) call open_file(IO_VOLC,"r",fname,needed=.true.)
    CALL MPI_BCAST(ios,1,MPI_INTEGER,0,MPI_COMM_WORLD,INFO)

    if(ios/=0)then
      if (MasterProc) print *,' No volcanoes found in VolcanoesLL.dat'
      return
    endif

    call Read_Headers(IO_VOLC,errmsg,NHeaders,NKeys,Headers,Keyvalues)
    VOLCLOOP: do nin = 1, NMAX_VOLC+1
      call read_line(IO_VOLC,txtinput,ios)
      if ( ios /= 0 ) exit  ! End of file
      nvolc=nvolc+1
      call CheckStop(nvolc>NMAX_VOLC,&
        "more Volcanoes in VolcanoesLL.dat than NMAX_VOLC")

      read(unit=txtinput,fmt=*) s,lon,lat,height_volc(nvolc),emis_volc(nvolc)
      call lb2ij(lon,lat,xr,yr)
      i_volc(nvolc)=nint(xr)
      j_volc(nvolc)=nint(yr)
      if (DEBUG.and.MasterProc) print *,'VOLC:Found ',trim(s), &
        lon, lat, i_volc(nvolc),j_volc(nvolc)&
        ,i_volc(nvolc)+IRUNBEG-1,j_volc(nvolc)+JRUNBEG-1, &
        height_volc(nvolc), emis_volc(nvolc)

    enddo VOLCLOOP
    if (MasterProc)close(IO_VOLC)

  endif

  tonne_to_kgm2s=1.0e3 / (nydays * 24.0 * 3600.0 * GRIDWIDTH_M * GRIDWIDTH_M)
  conv = tonne_to_kgm2s !DSRC * EmisDef(eindex_vol)%conv
  do volc_no=1,nvolc
    i=i_volc(volc_no)
    j=j_volc(volc_no)
    !Find global<->local coordinates for xm2
    if ((i>=gi0).and.(i<=gi1).and.(j>=gj0).and.(j<=gj1))then !on the correct processor
      Volcanoes_found=.true.
      i_l = i -gi0 +1
      j_l = j -gj0 +1
      if (DEBUG) then
        print *,'i,j for volcano ',i,j, 'EMIS_VOLC: ',emis_volc(volc_no)
        print *,'Volc Local coords are',i_l,j_l,gi0,gj0
      endif
      emis_volc(volc_no) = emis_volc(volc_no)* conv * xm2(i_l,j_l)
    endif
  enddo !volc_no

  !/** broadcast volcano heights
  CALL MPI_BCAST(height_volc,4*NMAX_VOLC,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

end subroutine VolcGet

subroutine Set_Volc
!-----------------------------------------------------------------------!
! Starts converting emission units from kg/m2/s to.... (hourly)
!-----------------------------------------------------------------------!
  !**Local variables
  integer :: i,j,k
  real    :: unit_conv1

  rcemis_volc0(:) = 0.0
  unit_conv1      = 0.0

  !/** Set volcano
  do volc_no=1,nvolc
    k=height_volc(volc_no)
    i=i_volc(volc_no)
    j=j_volc(volc_no)

    if (DEBUG) print '(A,4I6,6I6,4I6)','Volcan: check1 ',  &
      i,j, i_volc(volc_no),j_volc(volc_no),           &
      i_local(i),j_local(j), 1, limax, 1, ljmax,      &
      gi0,gi1,gj0,gj1

    if ((i_local(i)>=1).and.(i_local(i)<=limax).and.  &
        (j_local(j)>=1).and.(j_local(j)<=ljmax)) then

      unit_conv1 = GRAV* 0.001*AVOG/ &
                  (sigma_bnd(KMAX_BND-k+1) - sigma_bnd(KMAX_BND-k))

      rcemis_volc0(volc_no) = emis_volc(volc_no)*unit_conv1/species(SO2)%molwt

      if (DEBUG) print *,'rc_emis_volc0 is ',rcemis_volc0(volc_no)
    endif
  enddo ! volc_no
end subroutine Set_Volc

subroutine Scale_Volc
!-----------------------------------------------------------------------!
! Finishing converting volcano emissions to molecules/cm3/s
! (every advection timestep)
!-----------------------------------------------------------------------!

  integer :: i,j,k,i_l,j_l
  real    :: unit_conv2

  do volc_no=1,nvolc
    k=height_volc(volc_no)
    i=i_volc(volc_no)
    j=j_volc(volc_no)

    if (DEBUG) print '(A,4I6,6I6,4I6)','Volcan: check2 ', &
      i,j, i_volc(volc_no),j_volc(volc_no),           &
      i_local(i),j_local(j), 1, limax,1, ljmax,      &
      gi0,gi1,gj0,gj1

    if ((i_local(i)>=1).and.(i_local(i)<=limax).and.  &
        (j_local(j)>=1).and.(j_local(j)<=ljmax)) then
      i_l = i_local(i) !local i
      j_l = j_local(j) !local j

      if (DEBUG) print '(A,4I8)','Volcan: check 3: ',   &
        i_l, j_l, i_volc(volc_no)-gi0+1, j_volc(volc_no)-gj0+1

      unit_conv2 = roa(i_l,j_l,KMAX_BND-k,1) / (ps(i_l,j_l,1)-PT)

      rcemis_volc(volc_no) = rcemis_volc0(volc_no) * unit_conv2

      if (DEBUG) print *,'rc_emis_volc is ',rcemis_volc(volc_no)
     endif
  enddo ! volc_no
end subroutine Scale_Volc
!-----------------------------------------------------------------------!
! Emergency: Volcanic Eruption. Ash & SO2, other species are possible.
!-----------------------------------------------------------------------!
!----------------------------!
! Get Volcanic Eruption Emiss.
!----------------------------!
function EruptionRate(i,j) result(emiss)
  integer, intent(in) :: i,j
  real, dimension(NSPEC_SHL+1:NSPEC_TOT,KCHEMTOP:KMAX_MID) :: emiss
  character(len=*),parameter :: &
    MSG_FMT="('EMERGENCY:',1X,A,5(:,1X,I0,':',A),3(:,1X,G,':',A))"
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
  if(.not.Eruption_found)return
  snow=date2string(SDATE_FMT,current_date)
!----------------------------!
!
!----------------------------!
  doVENT: do v=1,nvent
    if(.not.coord_in_gridbox(ventdef(v)%lon,ventdef(v)%lat,i,j))&
                  cycle doVENT    ! Wrong gridbox
    if(nerup(v)<1)cycle doVENT    ! Not erupting
    if(DEBUG_EM) &
      write(*,MSG_FMT)snow//' Vent',me,'me',v,ventdef(v)%id,i,"i",j,"j"
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
      if(DEBUG_EM) &
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
  end function getModLev
!----------------------------!
! Set Volcanic Eruption Param.
!----------------------------!
  subroutine setRate()
    character(len=*),parameter :: &
      fventdef="volcanoes.csv",   &
      ferupdef="eruptions.csv",   &
      ERR_VENT_CSV="EruptionSet VENT def.file "//fventdef,          &
      ERR_ERUP_CSV="EruptionSet ERUP def.file "//ferupdef,          &
      ERR_VENT_MAX="EruptionSet NMAX_VENT exceeded in "//fventdef,  &
      ERR_ERUP_MAX="EruptionSet NMAX_ERUP exceeded in "//ferupdef
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
      Eruption_found=.false.
      if(MasterProc.and.DEBUG_EM.and.first_call) &
        write(*,MSG_FMT)'Code is off. No Volcanic Ash.'
      first_call=.false.
      return
    endif
    if(.not.first_call)then
      if(MasterProc.and.DEBUG_EM.and.second_call) &
        write(*,MSG_FMT)'No need for reset volc.def...'
      second_call=.false.
      return
    endif
    first_call=.false.
  !----------------------------!
  ! Read Vent CVS
  !----------------------------!
    if(DEBUG_EM) CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    if(MasterProc)then
      call open_file(IO_TMP,"r",fventdef,needed=.true.,iostat=stat)
      call CheckStop(stat,ERR_VENT_CSV//' not found')
    endif
    nvent=0
    doVENT: do l=1,NMAX_VENT+1
      call read_line(IO_TMP,txtline,stat)
      if(stat/=0) exit doVENT           ! End of file
      txtline=ADJUSTL(txtline)          ! Remove leading spaces
      if(txtline(1:1)=='#')cycle doVENT ! Comment line
      dvent=getVent(txtline)
      if(coord_in_processor(dvent%lon,dvent%lat))then
        nvent=nvent+1
        call CheckStop(nvent>NMAX_VENT,ERR_VENT_MAX//" read")
        ventdef(nvent)=dvent
        if(DEBUG_EM) &
          write(*,MSG_FMT)'Vent',me,'in',nvent,trim(dvent%id),&
            dvent%grp,trim(dvent%name)
!     else
!       if(DEBUG_EM) &
!         write(*,MSG_FMT)'Vent',me,'out',-1,trim(dvent%id),&
!           dvent%grp,trim(dvent%name)
      endif
    enddo doVENT
    if(MasterProc) close(IO_TMP)
    Eruption_found=(nvent>0)
  !----------------------------!
  ! Read Eruption CVS
  !----------------------------!
    if(DEBUG_EM) CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    if(MasterProc)then
      call open_file(IO_TMP,"r",ferupdef,needed=.true.,iostat=stat)
      call CheckStop(stat,ERR_ERUP_CSV//' not found')
    endif
    nerup(:)=0
    doERUP: do l=1,NMAX_ERUP+1
      call read_line(IO_TMP,txtline,stat)
      if(stat/=0) exit doERUP             ! End of file
      if(.not.Eruption_found)cycle doERUP ! There is no vents on subdomain
      txtline=ADJUSTL(txtline)            ! Remove leading spaces
      if(txtline(1:1)=='#')cycle doERUP   ! Comment line
      derup=getErup(txtline)
      if(derup%edef)then                                  ! Default
        nerup(0)=nerup(0)+1
        call CheckStop(nerup(0)>NMAX_ERUP,ERR_ERUP_MAX//" read")
        erupdef(0,nerup(derup%vent))=derup
        if(DEBUG_EM) &
          write(*,MSG_FMT)'Erup.Default',me,'in',nerup(0),derup%id
      elseif(derup%vent>0.and.derup%spc>0)then            ! Specific
        nerup(derup%vent)=nerup(derup%vent)+1
        call CheckStop(nerup(derup%vent)>NMAX_ERUP,ERR_ERUP_MAX//" read")
        erupdef(derup%vent,nerup(derup%vent))=derup
        if(DEBUG_EM) &
          write(*,MSG_FMT)'Erup.Specific',me,'in',nerup(derup%vent),trim(derup%id),&
            derup%spc,trim(derup%name)
!     else                                                ! Vent Outside domain
!       if(DEBUG_EM) &                                    ! or Unknown Vent/SPC
!         write(*,MSG_FMT)'Erup.Specific',me,'out',-1,trim(derup%id),&
!           derup%spc,trim(derup%name)
      endif
    enddo doERUP
    if(MasterProc) close(IO_TMP)
    Eruption_found=any(nerup(1:nvent)>0)
  !----------------------------!
  ! Expand Eruption Defaults
  !----------------------------!
    if(DEBUG_EM) CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    if(nerup(0)<1)then
      if(DEBUG_EM) &
        write(*,MSG_FMT)'Erup.Default',me,'not found'
      return
    endif
    doVENTe: do v=1,nvent
      if(nerup(v)>0)cycle doVENTe   ! Specific found --> no need for Default
      e=find_index(ventdef(v)%etype,erupdef(0,:nerup(0))%id)
      if(e<1)       cycle doVENTe   ! No Default found
      if(DEBUG_EM) &
        write(*,MSG_FMT)'Erup.Default',me,'Expand',&
          v,trim(ventdef(v)%id),e,trim(erupdef(0,e)%id)
      derup=erupdef(0,e)
      if(derup%spc<0.and.derup%name(1:3)=="ASH")then          ! Expand variable name
        derup%name=trim(ventdef(v)%id)//trim(derup%name(4:))  ! e.g. ASH_F --> V1702A02B_F
        derup%spc=find_index(derup%name,species(:)%name)      ! Specie (total)
      endif
      if(derup%spc>0)then           ! Expand single SPC
        nerup(v)=nerup(v)+1
        call CheckStop(nerup(v)>NMAX_ERUP,ERR_ERUP_MAX//" expand")
        erupdef(v,nerup(v))=derup
        if(DEBUG_EM) &
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
          if(DEBUG_EM) &
          write(*,MSG_FMT)'Erup.Default',me,'Expand GRP',nerup(v),trim(derup%id),&
            derup%spc,trim(derup%name)
        enddo
      else
        if(DEBUG_EM) &
          write(*,MSG_FMT)'Erup.Default',me,'not found'
      endif
    enddo doVENTe
    Eruption_found=any(nerup(1:nvent)>0)
  end subroutine setRate
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
  end function getVent
!----------------------------!
! Extract Erup. info from CVS line
!----------------------------!
  function getErup(line) result(def)
    character(len=*)            :: line
    type(erup)                  :: def
    character(len=TXTLEN_SHORT) :: words(10)=''  ! Array of paramaters
    logical :: edef=.true.                       ! default setings?
    integer :: stat,nwords,ivent,ispc=0,ind
    real    :: base,top,rate,m63,dhh
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
    read(words(5),*)dhh
    read(words(6),*)rate
    select case (words(7)) ! m63
    case(" ")   ;m63=1.0
    case default;read(words(7),*)m63
    endselect
    words(8)=getDate(words(8),words(8),words(9),dhh,debug=DEBUG_EM) ! Start [date/code]
    words(9)=getDate(words(9),words(8),words(9),dhh,debug=DEBUG_EM) ! End   [date/code]
    def=erup(trim(words(1)),trim(words(2)),base,top,rate*m63,&
      trim(words(8)),trim(words(9)),max(ivent,0),max(ispc,0),edef)
  end function getErup
!----------------------------!
! Time/Date CODE--> YYYY-MM-DD hh:mm:ss
!----------------------------!
  function getDate(code,se,ee,dh,debug) result(str)
    character(len=*), intent(in) :: code,se,ee
    real, intent(in)             :: dh  ! [hours]
    logical, intent(in), optional:: debug
    character(len=TXTLEN_SHORT)  :: str
    select case (code)
    case("SR")          ! Start of the simulation
      str=date2string(SDATE_FMT,startdate,debug=debug)
    case("SR+D")        ! Start of the simulation + dh
      str=date2string(SDATE_FMT,startdate,addsecs=dh*36e2,debug=debug)
    case("SE+D")        ! Start eruption + dh; no wildcards in SE allowed
      str=date2string(SDATE_FMT,string2date(se,SDATE_FMT,debug=debug),addsecs=dh*36e2,debug=debug)
    case("EE+D")        ! End eruption   - dh; no wildcards in EE allowed
      str=date2string(SDATE_FMT,string2date(ee,SDATE_FMT,debug=debug),addsecs=-dh*36e2,debug=debug)
    case("ER-D")        ! End of the simulation - dh
      str=date2string(SDATE_FMT,enddate,addsecs=-dh*36e2,debug=debug)
    case("ER")          ! End of the simulation
      str=date2string(SDATE_FMT,enddate,debug=debug)
    case default
      str=code
    endselect
  end function getDate
end function EruptionRate

end module Volcanos_ml
