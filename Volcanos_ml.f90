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
use GridValues_ml,        only: GRIDWIDTH_M, xm2, dA, dB,  &
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
use Par_ml,               only: me, IRUNBEG, JRUNBEG, ljmax, limax
use PhysicalConstants_ml, only: GRAV, AVOG
use TimeDate_ml,          only: nydays,&           ! No. days per year
                                startdate,enddate,current_date,&
                                make_timestamp,tdif_secs
use TimeDate_ExtraUtil_ml,only: date2string,string2date
use KeyValue_ml,          only: KeyVal

implicit none
private

 !/* subroutines:
public :: VolcGet,Set_Volc,Scale_Volc

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
  NMAX_VENT = 12, & ! Max number of Vents (erupting) on processor/subdomain
  NMAX_ERUP = 90    ! Max number of Eruption def (~3 months per vent)

integer, private, save ::   & ! No. of ... found on processor/subdomain
  nvent              = -1,  & ! Vents
  nerup(0:NMAX_VENT) = -1     ! Eruption events per Vent

logical, parameter :: &
  DEBUG_EM=DEBUG.or.DEBUG_EMERGENCY

logical, public, save ::      &
  Eruption_found=.true. !!! NML CHECK USE_EMERGENCY  ! Any Eruption found on this processor/subdomain?

type, private :: vent
  character(len=9)              :: id  =''    ! e.g. V1702A02B
  character(len=TXTLEN_NAME)    :: name=''    ! e.g. Eyjafjöll
  real :: lat=-1.0,lon=-1.0,elev=-1.0         ! vent coords and elevation
  character(len=2)              :: etype=''   ! e.g. S0
  integer :: grp=-1,iloc=-1,jloc=-1           ! Which (ash)goup,local i,j indes
 !character(len=TXTLEN_SHORT)   :: location,vtype  ! other info
endtype vent
type(vent), private, save, dimension(NMAX_VENT):: &
  ventdef=vent('UNDEF','UNKNOWN',-999.0,-999.0,-999.0,"??",-99,-99,-99)

character(len=*), private , parameter :: SDATE_FMT="YYYY-MM-DD hh:mm:ss"
type, private :: erup
  character(len=9)              :: id  =''    ! e.g. V1702A02B
  character(len=TXTLEN_NAME)    :: name=''    ! e.g. SO2
  real :: base=-1.0,top=-1.0,&                ! Column Base & Height [m asl]
          rate=-1.0                           ! Source strenght: Total release for period [kg/s]
  character(len=len(SDATE_FMT)) :: sbeg=SDATE_FMT,send=SDATE_FMT
  integer :: vent=-1,spc=-1                   ! Which vent,(adv)spc
  logical :: edef=.true.                      ! default setings?
endtype erup
type(erup), private, save, dimension(0:NMAX_VENT,NMAX_ERUP):: &
  erupdef=erup('UNDEF','UNKNOWN',-999.0,-999.0,-999.0,"??","??",-1,-1,.true.)

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
  integer            :: nin
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
    i=i_local(i_volc(volc_no))
    j=j_local(j_volc(volc_no))
    if(.not.all((/i>=1,i<=limax,j>=1,j<=ljmax/))) cycle 

    Volcanoes_found=.true. ! on the correct processor
    if(DEBUG)then
      print "(2(A,1X,I0,1X),1(A,':',G10.3))",&
        'Volc',volc_no,'found on',me,'EMIS_VOLC',emis_volc(volc_no)
      print "(2(A,':',I0,'x',I0,1X))",&
        'Volc Glob. coords',i_volc(volc_no),j_volc(volc_no),&
        'Volc Local coords',i,j
    endif
    emis_volc(volc_no) = emis_volc(volc_no)* conv * xm2(i,j)
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
    i=i_local(i_volc(volc_no))
    j=j_local(j_volc(volc_no))
    if(DEBUG) print '(A,10I6)','Volcan: check1 ',  &
      me,volc_no,i_volc(volc_no),j_volc(volc_no),i,j,1,limax,1,ljmax
    if(.not.all((/i>=1,i<=limax,j>=1,j<=ljmax/))) cycle 

    unit_conv1 = GRAV* 0.001*AVOG
    rcemis_volc0(volc_no) = emis_volc(volc_no)*unit_conv1/species(SO2)%molwt
    if(DEBUG) print *,'rc_emis_volc0 is ',rcemis_volc0(volc_no)
  enddo ! volc_no
end subroutine Set_Volc

subroutine Scale_Volc
!-----------------------------------------------------------------------!
! Finishing converting volcano emissions to molecules/cm3/s
! (every advection timestep)
!-----------------------------------------------------------------------!
  integer :: i,j,k
  real    :: unit_conv2

  do volc_no=1,nvolc
    k=height_volc(volc_no)
    i=i_local(i_volc(volc_no))
    j=j_local(j_volc(volc_no))
    if(DEBUG) print '(A,10I6)','Volcan: check2 ', &
      me,volc_no,i_volc(volc_no),j_volc(volc_no),i,j,1,limax,1,ljmax
    if(.not.all((/i>=1,i<=limax,j>=1,j<=ljmax/))) cycle 

    if(DEBUG) print '(A,2I8)','Volcan: check 3: ',i,j
    unit_conv2 = roa(i,j,KMAX_BND-k,1) /(dA(KMAX_BND-k)+dB(KMAX_BND-k)*ps(i,j,1))
    rcemis_volc(volc_no) = rcemis_volc0(volc_no) * unit_conv2
    if(DEBUG) print *,'rc_emis_volc is ',rcemis_volc(volc_no)
  enddo ! volc_no
end subroutine Scale_Volc
endmodule Volcanos_ml
