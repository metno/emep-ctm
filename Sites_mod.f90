! <Sites_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Sites_mod

! -----------------------------------------------------------------------
! Contains subroutines to read in list of measurement stations and/or
! radiosonde locations. Locations should be specified in the input
! files "sites.dat" and "sondes.dat"
!
! -----------------------------------------------------------------------
use Units_mod,          only: to_ug_ADV

use CheckStop_mod,      only: CheckStop, StopAll
use ChemDims_mod,       only: NSPEC_SHL, NSPEC_ADV
use ChemFunctions_mod,  only: Chem2Index_adv, Chem2Index
use ChemSpecs_mod
use ChemGroups_mod,     only: OXN_GROUP, PMFINE_GROUP, PMCO_GROUP
use Config_module,      only: NMET,PPBINV,PPTINV, KMAX_MID, MasterProc&
                              ,RUNDOMAIN, IOU_INST, SOURCE_RECEPTOR, meteo&
                              ,SitesFile,SondesFile,KMAX_BND,PT, NPROC,&  ! for sitesout
                              SITE_SHL_names,SONDE_SHL_names,SONDE_ADV_names,&
      NXTRA_SITE_MISC, NXTRA_SITE_D2D, &
      SITE_XTRA_MISC, SITE_XTRA_D2D, &
      FREQ_SITE, &
      NXTRA_SONDE, & 
       SONDE_XTRA, & 
      FREQ_SONDE
use Debug_module,       only: DEBUG   ! -> DEBUG%SITES
use DerivedFields_mod,  only: f_2d, d_2d  ! not used:, d_3d
use Functions_mod,      only: Tpot_2_T    ! Conversion function
use GridValues_mod,     only: lb2ij, i_fdom, j_fdom ,debug_proc &
                              ,i_local, j_local, A_mid, B_mid
use Io_mod,             only: check_file,open_file,ios &
                              , fexist, IO_SITES, IO_SONDES &
                              , Read_Headers,read_line
use KeyValueTypes,      only : KeyVal, KeyValue, LENKEYVAL
use MetFields_mod,      only : t2_nwp, th, pzpbl  &  ! output with concentrations
                              , z_bnd, z_mid, roa, Kz_m2s, q
use MetFields_mod,      only : u_xmj, v_xmi, ps
use MPI_Groups_mod,     only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_REAL8, MPI_INTEGER, MPI_LOGICAL, &
                             MPI_MIN, MPI_MAX, MPI_SUM, &
                             MPI_COMM_CALC, MPI_COMM_WORLD, MPISTATUS, IERROR, ME_MPI, NPROC_MPI
use PhysicalConstants_mod,only: ATWAIR
use NetCDF_mod,         only : Create_CDF_sondes,Out_CDF_sondes,&
                                NF90_FILL_INT,NF90_FILL_DOUBLE
use OwnDataTypes_mod,    only: TXTLEN_NAME
use Par_mod,            only : li0,lj0,li1,lj1 &
                                ,GIMAX,GJMAX,IRUNBEG,JRUNBEG&
                                ,GI0,GI1,GJ0,GJ1,me,LIMAX,LJMAX
use SmallUtils_mod,     only : find_index
use Tabulations_mod,    only : tab_esat_Pa
use TimeDate_mod,       only : current_date
use TimeDate_ExtraUtil_mod,   only : date2string

implicit none
private                     ! stops variables being accessed outside


! subroutines made available

public :: sitesdef          ! Calls Init_sites for sites and sondes
public :: siteswrt_surf     ! Gets site  data ready for siteswrt_out
public :: siteswrt_sondes   ! Gets sonde data ready for siteswrt_out
private :: Init_sites       ! reads locations, species
private :: set_species      ! Sets species/variable names for output
private :: siteswrt_out     ! Collects output from all nodes and prints


integer, public, parameter :: NSONDES_MAX = 99 ! Max. no sondes allowed
integer, public, parameter ::  NSITES_MAX = 99 ! Max. no surface sites allowed

integer, public, save :: nglobal_sites, nlocal_sites
integer,private, save :: nglobal_sondes, nlocal_sondes

! site_gindex stores the global index n asociated
! with each processor and local site

integer, private, save, allocatable,dimension (:,:)  :: site_gindex
integer, private, save, allocatable,dimension (:,:) :: sonde_gindex

integer, public, save, dimension (NSITES_MAX) :: &
        site_x, site_y, site_z      &! local coordinates
       , site_gn                        ! number in global
real, public, save, dimension (NSITES_MAX) :: &
  Sites_lon= -999, Sites_lat= -999

integer, private, save, dimension (NSITES_MAX) :: &
         site_gx, site_gy, site_gz    ! global coordinates
integer, private, save, dimension (NSONDES_MAX) ::  &
         sonde_gx, sonde_gy   &        ! global coordinates
       , sonde_x, sonde_y     &        ! local coordinates
       , sonde_gn                       ! number in global
real, public, save, dimension (NSONDES_MAX) :: &
  Sondes_lon= -999, Sondes_lat= -999, ps_sonde=0.0

integer, private :: NSPC_SITE, NOUT_SITE, NOUT_SONDE, NSPC_SONDE

character(len=TXTLEN_NAME), public, save, dimension(NSITES_MAX) :: site_name
character(len=TXTLEN_NAME), private, save, dimension(NSONDES_MAX):: sonde_name
character(len=20), private, save, allocatable, dimension(:)  :: site_species
character(len=20), private, save, allocatable, dimension(:)  :: sonde_species

character(len=70), private :: errmsg ! Message text
integer, private :: d                 ! processor index
integer, private :: i, n, nloc, ioerr ! general integers
integer, parameter, private :: Spec_Att_Size=40,N_Spec_Att_MAX=5
integer, save, private :: NSPECMAX
character(len=Spec_Att_Size), allocatable, save  :: Spec_Att(:,:)
integer, private :: i_Att !Spec attribute index
integer :: NSpec_Att !number of Spec attributes defined

integer, public, parameter ::  NADV_SITE  = NSPEC_ADV   ! No. advected species (1 up to NSPEC_ADV)
integer ::isite
integer, public, parameter, dimension(NADV_SITE) :: &
  SITE_ADV = [(isite, isite=1,NADV_SITE)]  ! Everything

integer, public :: NSHL_SITE, NADV_SONDE, NSHL_SONDE !number of requested species found
!indices of requested and found species
integer, public, dimension(NSPEC_SHL) :: SITE_SHL
integer, public, dimension(NSPEC_SHL) :: SONDE_SHL
integer, public, dimension(NSPEC_ADV) :: SONDE_ADV 

contains

!==================================================================== >
subroutine sitesdef()
  ! -------------------------------------------------------------------------
  ! reads in sites.dat and sondes.dat (if present), assigns sites to
  ! local domains, and collects lists of sites/species/variables for output.
  ! -------------------------------------------------------------------------

  ! Dummy arrays
  integer, save, dimension (NSONDES_MAX) ::  &
    sonde_gz, sonde_z ! global coordinates
  character(len=*), parameter :: &
    SITE_XTRA(NXTRA_SITE_MISC+NXTRA_SITE_D2D)=(/SITE_XTRA_MISC,SITE_XTRA_D2D/)

  sonde_gz(:) = 0
  sonde_z(:)  = 0

  allocate(site_gindex(0:NPROC-1,NSITES_MAX))
  allocate(sonde_gindex(0:NPROC-1,NSONDES_MAX))
  site_gindex = -1
  sonde_gindex= -1

  !find indices of species
  call Chem2Index(SITE_SHL_names,SITE_SHL,NSHL_SITE)
  call Chem2Index(SONDE_SHL_names,SONDE_SHL,NSHL_SONDE)
  call Chem2Index_adv(SONDE_ADV_names,SONDE_ADV,NADV_SONDE)
  NSPC_SITE  = NADV_SITE + NSHL_SITE + NXTRA_SITE_MISC + NXTRA_SITE_D2D 
  NOUT_SITE  = NSPC_SITE * 1 
  NSPC_SONDE = NADV_SONDE + NSHL_SONDE + NXTRA_SONDE
  NOUT_SONDE = NSPC_SONDE* KMAX_MID
  NSPECMAX=max(NSPC_SITE,NSPC_SONDE)
  if(.not.allocated(site_species))allocate(site_species(NSPC_SITE))
  if(.not.allocated(sonde_species))allocate(sonde_species(NSPC_SONDE))
  if(.not.allocated(Spec_Att))allocate(Spec_Att(NSPECMAX,N_Spec_Att_MAX))

!could use XXX_names, instead of site_species
  site_species(1:NADV_SITE) = species_adv(SITE_ADV(1:NADV_SITE))%name
  site_species(NADV_SITE+1:NADV_SITE+NSHL_SITE) = &
            species(SITE_SHL(1:NSHL_SITE))%name
  site_species(NADV_SITE+NSHL_SITE+1:NSPC_SITE) = &
            SITE_XTRA(1:NXTRA_SITE_MISC+NXTRA_SITE_D2D)
  sonde_species(1:NADV_SONDE) = species_adv(SONDE_ADV(1:NADV_SONDE))%name
  sonde_species(NADV_SONDE+1:NADV_SONDE+NSHL_SONDE) = &
            species(SONDE_SHL(1:NSHL_SONDE))%name
  sonde_species(NADV_SONDE+NSHL_SONDE+1:NSPC_SONDE) = &
            SONDE_XTRA(1:NXTRA_SONDE)

  call Init_sites(SitesFile,IO_SITES,NSITES_MAX, &
        nglobal_sites,nlocal_sites, &
        site_gindex, site_gx, site_gy, site_gz, &
        site_x, site_y, site_z, site_gn, &
        site_name)

  call Init_sites(SondesFile,IO_SONDES,NSONDES_MAX, &
        nglobal_sondes,nlocal_sondes, &
        sonde_gindex, sonde_gx, sonde_gy, sonde_gz, &
        sonde_x, sonde_y, sonde_z, sonde_gn, &
        sonde_name)

!  call set_species(SITE_ADV,SITE_SHL,SITE_XTRA,site_species)
!  call set_species(SONDE_ADV,SONDE_SHL,SONDE_XTRA,sonde_species)

  if ( DEBUG%SITES ) then
     write(6,*) "sitesdef After nlocal ", nlocal_sites, " on me ", me
     do i = 1, nlocal_sites
       write(6,*) "sitesdef After set_species x,y ", &
                     site_x(i), site_y(i),site_z(i), " on me ", me
     end do
  end if ! DEBUG

end subroutine sitesdef
!==================================================================== >
subroutine set_species(adv,shl,xtra,s)
  ! ---------------------------------------------------------------------
  ! Makes a character array "s" containg the names of the species or
  ! meteorological parameters to be output. Called for sites and sondes.
  ! ---------------------------------------------------------------------

  integer, intent(in), dimension(:) :: adv, shl ! Arrays of indices wanted
  character(len=*), intent(in), dimension(:) :: xtra  !Names of extra params
  character(len=*),intent(out), dimension(:) :: s

  integer :: nadv, nshl, n2, nout  ! local sizes
  nadv = size(adv)
  nshl = size(shl)
  n2 = nadv + nshl
  nout = size(s)                   ! Size of array to be returned

  s(1:nadv)    = species( NSPEC_SHL + adv(:) )%name
  s(nadv+1:n2) = species( shl(1:) )%name
  s(n2+1:nout) = xtra(:)
end subroutine set_species
!==================================================================== >
subroutine Init_sites(fname,io_num,NMAX, nglobal,nlocal, &
        s_gindex, s_gx, s_gy, s_gz, s_x, s_y, s_z, s_n, s_name)
  ! ----------------------------------------------------------------------
  ! Reads the file "sites.dat" and "sondes.dat" to get coordinates of
  ! surface measurement stations or locations where vertical profiles
  ! or extra output are required. (These files may be empty, but this is
  ! not recommended - the sites data provide good diagnostics).
  !
  !  define whether a certain output site belongs to the given processor
  !  and assign the local coordinates
  ! ----------------------------------------------------------------------
  ! NB. global below refers to all nodes (full-domain)
  !     local  below refers to the local node

  character(len=*), intent(in) :: fname
  integer,          intent(in) :: io_num
  integer,          intent(in) :: NMAX     ! Max no. sites
  integer, intent(out) :: nglobal, nlocal  ! No. sites
  integer, intent(out), dimension (0:,:) :: s_gindex  ! index, starts at me=0
  integer, intent(out), dimension (:) ::  &
                            s_gx, s_gy, s_gz   & ! global coordinates
                           ,s_x, s_y, s_z      & ! local coordinates
                           ,s_n                  ! number in global
  character(len=*), intent(out), dimension (:) :: s_name

  !-- Local:
  integer,  dimension (NMAX) :: s_n_recv  ! number in global

  integer           :: nin     ! loop index
  integer           :: ix, iy  ! coordinates read in
  integer           :: lev     ! vertical coordinate (20=ground)
  character(len=20) :: s       ! Name of site read in
  character(len=30) :: comment ! comment on site location
  character(len=40) :: errmsg
  real              :: lat,lon
  character(len=*),parameter :: dtxt='SitesInit:'

  character(len=20), dimension(4) :: Headers
  type(KeyVal), dimension(20)     :: KeyValues ! Info on units, coords, etc.
  integer                         :: NHeaders, NKeys
  character(len=80)               :: txtinput  ! Big enough to contain
                                               ! one full input record


  ios = 0                      ! zero indicates no errors
  errmsg = ''
  if ( MasterProc ) then
    call check_file(fname,fexist,needed=.false.,errmsg=errmsg)
    if ( .not.fexist )then
      ios=1
    else
      call open_file(io_num,"r",fname,needed=.true.)
      call CheckStop(ios,"ios error on "//trim(fname))
    end if
  end if

  call MPI_BCAST( ios, 1, MPI_INTEGER, 0, MPI_COMM_CALC,IERROR)
  if(ios/=0)return

  call CheckStop(NMAX,size(s_name), dtxt//"Error : sitesdefNMAX problem")

  call Read_Headers(io_num,errmsg,NHeaders,NKeys,Headers,Keyvalues)

  ! First, see which sites are within full domain:

  n = 0          ! Number of sites found within domain

  SITELOOP: do nin = 1, NMAX

    if (trim(KeyValue(KeyValues,"Coords"))=='LatLong') then
      call read_line(io_num,txtinput,ios)
      if ( ios /= 0 ) exit  ! End of file
      read(unit=txtinput,fmt=*) s, lat, lon, lev
      call lb2ij(lon,lat,ix,iy)
    else
      call read_line(io_num,txtinput,ios)
      lon=-999.0
      lat=-999.0
      if ( ios /= 0 ) exit  ! End of file
      read(unit=txtinput,fmt=*) s,  ix,  iy, lev
    end if

    if (ioerr < 0) then
      write(6,*) dtxt//" end of file after ", nin-1, trim(fname)
      exit SITELOOP
    end if ! ioerr


    if ( ix<RUNDOMAIN(1) .or. ix>RUNDOMAIN(2) .or. &
         iy<RUNDOMAIN(3) .or. iy>RUNDOMAIN(4) ) then
      if(MasterProc) write(6,"(A,': ',A,2(A,1X,I0),A)") &
        " site", trim(s), ', i =',ix, ', j =',iy, ", outside computational domain"
    elseif ( ix==RUNDOMAIN(1) .or. ix==RUNDOMAIN(2) .or. &
             iy==RUNDOMAIN(3) .or. iy==RUNDOMAIN(4) ) then
      if(MasterProc) write(6,"(A,': ',A,2(A,1X,I0),A)") &
        " site", trim(s), ', i =',ix, ', j =',iy, ", on domain boundary"
    else
      comment = " ok - inside domain         "
      n = n + 1
      
      s_gx(n)   = ix
      s_gy(n)   = iy
      s_gz(n)   = lev

      if(trim(fname)==trim(SitesFile))then
         if(lon>-990)Sites_lon(n) = lon
         if(lat>-990)Sites_lat(n) = lat
      end if
      if(trim(fname)==trim(SondesFile))then
         if(lon>-990)Sondes_lon(n) = lon
         if(lat>-990)Sondes_lat(n) = lat
      end if

      s_name(n)  = s !!! remove comments// comment
      if (DEBUG%SITES.and.MasterProc) write(6,"(a,i4,a)") dtxt//" s_name : ",&
            n, trim(s_name(n))
    end if

  end do SITELOOP

  nglobal = n

  ! NSITES/SONDES_MAX must be _greater_ than the number used, for safety

  call CheckStop(n >= NMAX, dtxt//"Error : increase NGLOBAL_SITES_MAX!")

  if(MasterProc) close(unit=io_num)

  nlocal  = 0

  do n = 1, nglobal

    ix = s_gx(n) ! global-domain coords
    iy = s_gy(n)

    if ( i_local(ix)>=li0 .and. i_local(ix)<=li1 .and. &
         j_local(iy)>=lj0 .and. j_local(iy)<=lj1 ) then

      nlocal      = nlocal + 1
      s_x(nlocal) = i_local(ix)
      s_y(nlocal) = j_local(iy)
      s_z(nlocal) = s_gz(n)
      s_n(nlocal) = n

      if (DEBUG%SITES) &
        write(6,"(a,i3,a,2i3,3i4,a,3i4)") dtxt//" Site on me : ", me, &
         " Nos. ", n, nlocal, s_gx(n), s_gy(n) , s_gz(n), " =>  ", &
          s_x(nlocal), s_y(nlocal), s_z(nlocal)
        write(6,"(a,i3,a,2i3,4a)") dtxt// trim(fname), me, &
         " Nos. ", n, nlocal, " ", trim(s_name(n)), " => ", trim(s_name(s_n(nlocal)))

     end if

  end do ! nglobal

  ! inform me=0 of local array indices:
  if(DEBUG%SITES) write(6,*) dtxt//trim(fname), " before gc NLOCAL_SITES", &
                           me, nlocal

  if ( .not.MasterProc ) then
    call MPI_SEND(nlocal, 4*1, MPI_BYTE, 0, 333, MPI_COMM_CALC, IERROR)
    if(nlocal>0) call MPI_SEND(s_n, 4*nlocal, MPI_BYTE, 0, 334, &
                               MPI_COMM_CALC, IERROR)
  else
    if(DEBUG%SITES) write(6,*) dtxt//" for me =0 LOCAL_SITES", me, nlocal
    do n = 1, nlocal
      s_gindex(me,n) = s_n(n)
    end do
    do d = 1, NPROC-1
      call MPI_RECV(nloc, 4*1, MPI_BYTE, d, 333, MPI_COMM_CALC,MPISTATUS, IERROR)
      if(nloc>0) call MPI_RECV(s_n_recv, 4*nloc, MPI_BYTE, d, 334, &
                               MPI_COMM_CALC,MPISTATUS, IERROR)
      if(DEBUG%SITES) write(6,*) dtxt//" recv d ", fname, d,  &
                  " zzzz nloc : ", nloc, " zzzz me0 nlocal", nlocal
      do n = 1, nloc
        s_gindex(d,n) = s_n_recv(n)
        if(DEBUG%SITES) write(6,*) dtxt//" for d =", fname, d, &
          " nloc = ", nloc, " n: ",  n,  " gives nglob ", s_gindex(d,n)
      end do ! n
    end do ! d
  end if ! MasterProc

  if ( DEBUG%SITES ) write(6,*) dtxt//' on me', me, ' = ', nlocal

end subroutine Init_sites
!==================================================================== >
subroutine siteswrt_surf(xn_adv,cfac,xn_shl)
  ! -------------------------------------------------------------------
  ! writes out just simple concentrations for now....
  ! will be improved later to allow choice of output parameter
  ! should look at chemint also - seems similar for somethings
  ! -------------------------------------------------------------------

  ! arguments
  real, dimension(NSPEC_ADV,LIMAX,LJMAX,KMAX_MID), intent(in) :: xn_adv
  real, dimension(NSPEC_ADV,LIMAX,LJMAX), intent(in)          :: cfac
  real, dimension(NSPEC_SHL,LIMAX,LJMAX,KMAX_MID), intent(in) :: xn_shl

  ! Local
  integer :: ix, iy,iz, ispec                  ! Site indices
  integer :: nn                                ! species index
  logical, save :: my_first_call = .true.      ! for debugging
  integer                           :: d2index ! index for d_2d field access
  character(len=len(SITE_XTRA_D2D)) :: d2code  ! parameter code -- # --
  character(len=*),parameter :: dtxt = 'siteswrt_surf:'

  real,dimension(NOUT_SITE,NSITES_MAX) :: out  ! for output, local node

  if ( DEBUG%SITES ) then
    write(6,*) dtxt//"nlocal ", nlocal_sites, " on me ", me
    do i = 1, nlocal_sites
      write(6,*) dtxt//"x,y ",site_x(i),site_y(i),&
                  site_z(i)," me ", me
    end do

    if ( MasterProc ) then
      write(6,*) "======= site_gindex ======== sitesdef ============"
      do n = 1, nglobal_sites
        write(6,*) dtxt//"def ", n, NPROC, (site_gindex(d,n),d=0,4)
        write(6,'(a12,i4,2x,200i4)') dtxt//"def ", n, &
                (site_gindex(d,n),d=0,NPROC-1)
      end do
      write(6,*) "======= site_end    ======== sitesdef ============"
    end if ! MasterProc
  end if ! DEBUG

  ! assign local data to out

  i_Att=0
  NSpec_Att=1 !number of Spec attributes defined
  do i = 1, nlocal_sites
    ix = site_x(i)
    iy = site_y(i)
    iz = site_z(i)
    if( iz == 0 ) iz = KMAX_MID  ! If ZERO'd, skip surface correction

    i_Att=0
    do ispec = 1, NADV_SITE
      !if (iz == KMAX_MID ) then ! corrected to surface
      if (site_z(i) == KMAX_MID ) then ! corrected to surface
        out(ispec,i) = xn_adv( SITE_ADV(ispec) ,ix,iy,KMAX_MID ) * &
                       cfac( SITE_ADV(ispec),ix,iy) * PPBINV
      else                      ! Mountain sites not corrected to surface
        out(ispec,i)  = xn_adv( SITE_ADV(ispec) ,ix,iy,iz ) * PPBINV
      end if
      i_Att=i_Att+1
      Spec_Att(i_Att,1)='units:C:ppb'
    end do


    do ispec = 1, NSHL_SITE
      out(NADV_SITE+ispec,i)  = xn_shl( SITE_SHL(ispec) ,ix,iy,iz )
      i_Att=i_Att+1
      Spec_Att(i_Att,1)='units:C:molecules/cm3'
    end do

    ! XTRA parameters, usually the temmp or pressure
    if (.not. SOURCE_RECEPTOR) then
      nn = NADV_SITE + NSHL_SITE
      do ispec = 1, NXTRA_SITE_MISC
        nn=nn+1
        select case(SITE_XTRA_MISC(ispec))
        case("T2")
          out(nn,i)   = t2_nwp(ix,iy,1) - 273.15
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:K'
        case("th")
          out(nn,i)   = th(ix,iy,iz,1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:K'
!       case("hmix")
!         out(nn,i)   = pzpbl(ix,iy)
        case default
          call CheckStop("Error, Sites_mod/siteswrt_surf: SITE_XTRA_MISC:"&
                               // trim(SITE_XTRA_MISC(ispec)))
        end select
        call CheckStop( abs(out(nn,i))>1.0e99, &
          "ABS(SITES OUT: '"//trim(SITE_XTRA_MISC(ispec))//"') TOO BIG" )
      end do
      do ispec = 1, NXTRA_SITE_D2D
        d2code  = SITE_XTRA_D2D(ispec)
        d2index = find_index(d2code, f_2d(:)%name)
        nn=nn+1

        !call CheckStop(d2index<1,"SITES D2D NOT FOUND"//trim(d2code))
        if(d2index<1) then
          if(MasterProc.and.my_first_call) &
            write(*,*) dtxt//"WARNING: D2D NOT FOUND"//trim(d2code)
          !cycle
          out(nn,i) = NF90_FILL_DOUBLE
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='output:C:undefined'  !to improve?
        else
          out(nn,i) = d_2d(d2index,ix,iy,IOU_INST)*f_2d(d2index)%scale
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:'//trim(f_2d(d2index)%unit)
        end if

        if( DEBUG%SITES ) &
          write(6,"(a,3i4,a15,i4,es12.3)") dtxt//"D2DEBUG ", me, nn, i,&
            " "//trim(d2code), d2index, out(nn,i)
        call CheckStop( abs(out(nn,i))>1.0e99, &
          dtxt//"ABS(SITES OUT: '"//trim(SITE_XTRA_D2D(ispec))//"') TOO BIG" )
      end do
    end if
  end do

  my_first_call = .false.
  ! collect data into gout on me=0 t
  call siteswrt_out("sites",IO_SITES,NOUT_SITE, FREQ_SITE, &
                     nglobal_sites,nlocal_sites, &
                     site_gindex,site_name,site_gx,site_gy,site_gz,&
                     site_species,out,ps_sonde)
end subroutine siteswrt_surf
!==================================================================== >
subroutine siteswrt_sondes(xn_adv,xn_shl)
  ! -------------------------------------------------------------------
  ! Writes vertical concentration  data to files.
  ! IO_SONDES is set in io_mod to be 30
  ! -------------------------------------------------------------------

  real, dimension(:,:,:,:), intent(in) ::  xn_adv
  real, dimension(:,:,:,:), intent(in) ::  xn_shl

  ! Output variables - none

  ! Local variables
  integer :: n, i, k,  ix, iy, nn, ispec   ! Site and chem indices
  integer ::  KTOP_SONDE 
  integer, dimension(KMAX_MID)      :: itemp
  real, dimension(KMAX_MID)              :: pp, temp, qsat, rh, sum_PM, sum_NOy
  real, dimension(NOUT_SONDE,NSONDES_MAX):: out

  out=0.0
  ! Consistency check
  KTOP_SONDE = KMAX_MID - KMAX_MID + 1
  do ispec=1,NXTRA_SONDE
    select case(SONDE_XTRA(ispec))
    case("PM25","PMco","NOy","SH","RH","roa","dz","z_mid","p_mid","Kz_m2s","th","U","V")
      errmsg = "ok"
    case("T")
      call CheckStop(all(SONDE_XTRA(1:ispec)/="RH"),&
        "Error, Sites_mod/siteswrt_sondes SONDE_XTRA: '"//&
        trim(SONDE_XTRA(ispec))//"' needs to be requested after 'RH'")
      errmsg = "ok"
    case default
      call CheckStop("Error, Sites_mod/siteswrt_sondes SONDE_XTRA: "//&
          trim(SONDE_XTRA(ispec)))
    end select
  end do

  i_Att=0
  NSpec_Att=1 !number of Spec attributes defined
  do i = 1, nlocal_sondes
    n  = sonde_gn(i)
    ix = sonde_x(i)
    iy = sonde_y(i)
    nn = 0

    ! collect and print out with ground-level (KMAX_MID) first, hence &
    ! KMAX_MID:KTOP_SONDE:-1 in arrays
    ! first the advected and short-lived species
    i_Att=0
    do ispec = 1, NADV_SONDE    !/ xn_adv in ppb
      out(nn+1:nn+KMAX_MID,i) = PPBINV *  &
          xn_adv( SONDE_ADV(ispec) , ix,iy,KMAX_MID:KTOP_SONDE:-1)
      nn = nn + KMAX_MID
      i_Att=i_Att+1
      Spec_Att(i_Att,1)='units:C:ppb'
    end do

    do ispec = 1, NSHL_SONDE    !/ xn_shl  in molecules/cm3
      out(nn+1:nn+KMAX_MID,i) = xn_shl( SONDE_SHL(ispec) , &
          ix,iy,KMAX_MID:KTOP_SONDE:-1)
      nn = nn + KMAX_MID
      i_Att=i_Att+1
      Spec_Att(i_Att,1)='units:C:molecules/cm3'
    end do

    ! then print out XTRA stuff first,
    ! usually the height or pressure

    do ispec = 1, NXTRA_SONDE
      select case(SONDE_XTRA(ispec))
        case("PM25")  !!  PM data converted to ug m-3
          sum_PM(:) = 0.
          do k = 1, KMAX_MID
            sum_PM(k) = (dot_product(xn_adv(PMFINE_GROUP-NSPEC_SHL,ix,iy,k),&
                                  to_ug_ADV(PMFINE_GROUP-NSPEC_SHL)) &
                          + 0.5 * &  ! 50% of PMcoare in PM2.5, since Dp=2.5
                          dot_product(xn_adv(PMCO_GROUP-NSPEC_SHL,ix,iy,k),&
                                  to_ug_ADV(PMCO_GROUP-NSPEC_SHL)) &
                        ) * roa(ix,iy,k,1)
          end do !k
!bug?          out(nn+1:nn+KMAX_MID,i) = sum_PM(KMAX_MID:1:-1)
          out(nn+1:nn+KMAX_MID,i) = sum_PM(KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:ug/m3'

        case("PMco")  !!  PM data converted to ug m-3
          sum_PM(:) = 0.
          do k = 1, KMAX_MID
            sum_PM(k) = dot_product(xn_adv(PMCO_GROUP-NSPEC_SHL,ix,iy,k),&
                                  to_ug_ADV(PMCO_GROUP-NSPEC_SHL)) &
                      * roa(ix,iy,k,1)
          end do !k
!bug?          out(nn+1:nn+KMAX_MID,i) = sum_PM(KMAX_MID:1:-1)
          out(nn+1:nn+KMAX_MID,i) = sum_PM(KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:ug/m3'

        case("NOy")
          sum_NOy(:) = 0.
          do k = 1, KMAX_MID
            sum_NOy(k) = sum(xn_adv(OXN_GROUP-NSPEC_SHL,ix,iy,k))
          end do
!bug?          out(nn+1:nn+KMAX_MID,i) = PPBINV * sum_NOy(KMAX_MID:1:-1)
          out(nn+1:nn+KMAX_MID,i) = PPBINV * sum_NOy(KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:mix_ratio'

        case("SH")   ! Specific Humidity
          out(nn+1:nn+KMAX_MID,i) =  q(ix,iy,KMAX_MID:KTOP_SONDE:-1,1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:kg kg-1' 

        case("RH")   ! Relative Humidity
          do k = 1,KMAX_MID
            pp(k) = A_mid(k) + B_mid(k)*ps(ix,iy,1)
            temp(k) = th(ix,iy,k,1)* Tpot_2_T(pp(k))
            itemp(k)= nint(temp(k))
            qsat(k) = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
            rh(k) = min(q(ix,iy,k,1)/qsat(k),1.0)
          end do
          out(nn+1:nn+KMAX_MID,i) = rh(KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:fraction'

        case("roa")
          out(nn+1:nn+KMAX_MID,i)=roa(ix,iy,KMAX_MID:KTOP_SONDE:-1,1)  
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:kg/m3'

        case("dz")  ! level thickness
          out(nn+1:nn+KMAX_MID,i)=z_bnd(ix,iy,KMAX_MID+0:KTOP_SONDE+0:-1)& 
                                      -z_bnd(ix,iy,KMAX_MID+1:KTOP_SONDE+1:-1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:m'

        case("z_mid")
          out(nn+1:nn+KMAX_MID,i) =  z_mid(ix,iy,KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:m'

        case("p_mid")
          out(nn+1:nn+KMAX_MID,i) = A_mid(KMAX_MID:KTOP_SONDE:-1) + &
                                    B_mid(KMAX_MID:KTOP_SONDE:-1)*ps(ix,iy,1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:Pa'

        case("Kz_m2s")
          out(nn+1:nn+KMAX_MID,i) =  Kz_m2s(ix,iy,KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:?/m2/s'

        case("th")
          out(nn+1:nn+KMAX_MID,i) =  th(ix,iy,KMAX_MID:KTOP_SONDE:-1,1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:K'

        case("T")
          out(nn+1:nn+KMAX_MID,i) =  temp(KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:K'

        case("U")
          out(nn+1:nn+KMAX_MID,i) = 0.5 &
             *( u_xmj(ix,iy,KMAX_MID:KTOP_SONDE:-1,1) &
             +u_xmj(ix-1,iy,KMAX_MID:KTOP_SONDE:-1,1) )
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:m/s'

        case("V")
          out(nn+1:nn+KMAX_MID,i) = 0.5 &
             *( v_xmi(ix,iy,KMAX_MID:KTOP_SONDE:-1,1) &
             +v_xmi(ix,iy-1,KMAX_MID:KTOP_SONDE:-1,1) )
          i_Att=i_Att+1
          Spec_Att(i_Att,1)='units:C:m/s'

        case("D3D")
          call StopAll("D3D Sites out not defined")
      end select

      nn=nn+KMAX_MID
    end do ! ispec (NXTRA_SONDE)
    
    ps_sonde(i)=ps(ix,iy,1)!surface pressure always needed to define the vertical levels
  end do ! i (nlocal_sondes)

  ! collect data into gout on me=0 t

  call siteswrt_out("sondes",IO_SONDES,NOUT_SONDE, FREQ_SONDE, &
                     nglobal_sondes,nlocal_sondes, &
                     sonde_gindex,sonde_name,sonde_gx,sonde_gy,sonde_gy, &
                     sonde_species,out,ps_sonde)

end subroutine siteswrt_sondes
!==================================================================== >
subroutine siteswrt_out(fname,io_num,nout,f,nglobal,nlocal, &
     s_gindex,s_name,s_gx,s_gy,s_gz,s_species,out,ps_sonde)
  ! -------------------------------------------------------------------
  ! collects data from local nodes and writes out to sites/sondes.dat
  ! -------------------------------------------------------------------

  character(len=*), intent(in) :: fname
  integer, intent(in) :: io_num, nout
  integer, intent(in) :: f               ! Frequency of write-out (hours)
  integer, intent(in) :: nglobal, nlocal
  integer, intent(in), dimension (0:,:) :: s_gindex  ! index, starts at me=0
  character(len=*), intent(in), dimension (:) ::  s_name    ! site/sonde name
  integer, intent(in), dimension (:) :: s_gx, s_gy, s_gz    ! coordinates
  character(len=*), intent(in), dimension (:) ::  s_species ! Variable names
  real,    intent(in), dimension(:,:) :: out    ! outputs, local node
  real,    intent(in), dimension(:) ::  ps_sonde   ! surface pressure local node

  ! Local
  real,dimension(nout,nglobal) :: g_out ! for output, collected
  real,dimension(nglobal) :: g_ps ! for ps, collected
  integer :: nglob, nloc         ! Site indices
  character(len=40)  :: outfile
  character(len=4)   :: suffix
  integer, parameter :: NTYPES = 2      ! No. types, now 2 (sites, sondes)
  integer ::  type=-1                   ! = 1 for sites, 2 for sondes
  integer, save, dimension(NTYPES):: prev_year = (/ -99, -99 /) ! Initialise
  integer :: ii

  integer, parameter :: NattributesMAX=10
  character(len=200),allocatable  :: SpecName(:),SpecDef(:,:),MetaData(:,:)
  character(len=200)              :: fileName
  integer  :: Nlevels,ispec,NSPEC,NStations,NMetaData
  integer ::i_Att_MPI
  logical :: debug_1d=.false.

  select case (fname)
  case("sites") ;type=1
  case("sondes");type=2
  case default
    write(6,*) "non-possible type in siteswrt_out for ", fname
    return
  end select

  write(suffix,fmt="(i4)") prev_year(type)
  fileName = fname // "_" // suffix // ".nc"!Name of the NetCDF file. Will overwrite any preexisting file

  !  if ( MasterProc .and. current_date%month /= prev_month(type)) then
  if(current_date%year/=prev_year(type) & !We do not consider midnight as a new year
    .and.(current_date%month/=1.or.current_date%hour/=0.or.current_date%seconds/=0)&
       ) then
    prev_year(type) = current_date%year
    write(suffix,fmt="(i4)") current_date%year
    fileName = fname // "_" // suffix // ".nc"!Name of the NetCDF file. Will overwrite any preexisting file

    if(MasterProc) then
      if(prev_year(type)>0) close(io_num)  ! Close last-year file
      prev_year(type) = current_date%year

      ! Open new file for write-out
      write(suffix,fmt="(i4)") current_date%year
      outfile = fname // "_" // suffix // ".csv"

      open(file=outfile,unit=io_num,action="write",form='FORMATTED')
      write(io_num,"(i3,2x,a,a, 4i4)") nglobal, fname, " in domain",RUNDOMAIN
      write(io_num,"(i3,a)") f, " Hours between outputs"
      do n = 1, nglobal
        write(io_num,'(a50,3(",",i4))') s_name(n), s_gx(n), s_gy(n),s_gz(n)
      end do ! nglobal

      write(io_num,'(i3,a)') size(s_species), " Variables units: ppb"
      !MV write(io_num,'(a9,<size(s_species)>(",",a))')"site,date",(trim(s_species(i)),i=1,size(s_species))
      !MAY2019 write(io_num,'(9999a)')"site,date", (",", (trim(s_species(i)) ),i=1,size(s_species))
      write(io_num,'(9999a)')"site,date,hh", (",", (trim(s_species(i)) ),i=1,size(s_species))

      !defintions of file for NetCDF output
      select case(fname)
      case("sondes")
        NLevels=KMAX_MID !number of vertical levels (counting from surface)
        NSPEC=NSPC_SONDE  !number of species defined for sondes
      case("sites")
        NLevels=1
        NSPEC=NSPC_SITE !number of species defined for sites
      end select
      NStations = nglobal !number of sondes or sites defined

      allocate(SpecDef(NSPEC,0:NattributesMAX),MetaData(0:NStations,NattributesMAX))
      SpecDef=""
      MetaData=""

      NMetaData=6
      if(IRUNBEG>1)NMetaData=NMetaData+1
      if(JRUNBEG>1)NMetaData=NMetaData+1
      call CheckStop(NattributesMAX<NMetaData,'NattributesMAX too small')
      
      write(MetaData(0,1),"(A,':C:',A)")"File_Type",trim(fname)
      write(MetaData(0,2),"(A,':C:',A)")"meteo_source",trim(meteo)
      write(MetaData(0,3),"(A,':I:',I0)")"Number_of_hours_bewtween_outputs",f
      write(MetaData(0,4),"(A,':I:',I0)")"Number_of_stations_defined",NStations
      write(MetaData(0,5),"(A,':I:',I0)")"Model_domain_x_size",GIMAX
      write(MetaData(0,6),"(A,':I:',I0)")"Model_domain_y_size",GJMAX
      if(IRUNBEG>1)&
        write(MetaData(0,7),"(A,':I:',I0)")"Model_domain_x_shift_origin",IRUNBEG
      if(JRUNBEG>1)&
        write(MetaData(0,8),"(A,':I:',I0)")"Model_domain_y_shift_origin",JRUNBEG

      call CheckStop(NattributesMAX<8,'NattributesMAX too small')
      NMetaData=MAX(NMetaData,1+5) ! metadata: 1 sting + 4..5 real
      call CheckStop(NattributesMAX<NMetaData,'Coords: NattributesMAX too small')

      do n = 1, NStations
        ! check to avoid two stations with same name, 
        ! since this will give problems when writing in the netcdf output
        call CheckStop(COUNT(s_name(1:NStations)==s_name(n))>1,&
          "Two stations named '"//trim(s_name(n))//"'. Remove or rename one of them")

!eg     write(MetaData(n,?),"(A,':C:',A)")"station_type","urban"
        write(MetaData(n,1),"(A,':C:',A)")"station_name",trim(s_name(n))
        write(MetaData(n,2),"(A,':I:',I0)")"model_x_coordinate",s_gx(n)
        write(MetaData(n,3),"(A,':I:',I0)")"model_y_coordinate",s_gy(n)
        write(MetaData(n,4),"(A,':I:',I0)")"model_level",s_gz(n)
        write(MetaData(n,5),"(A,':D:',A)")"longitude","missing"
        write(MetaData(n,6),"(A,':D:',A)")"latitude" ,"missing"
        select case(fname)
        case("sites")
          if(sites_lon(n)>-990)&
            write(MetaData(n,5),"(A,':D:',F10.3)")"longitude",sites_lon(n)
          if(sites_lat(n)>-990)&
            write(MetaData(n,6),"(A,':D:',F10.3)")"latitude" ,sites_lat(n)
        case("sondes")
          MetaData(n,4)=""  ! skip model_level for sondes
          if(sondes_lon(n)>-990)&
            write(MetaData(n,5),"(A,':D:',F10.3)")"longitude",sondes_lon(n)
          if(sondes_lat(n)>-990)&
            write(MetaData(n,6),"(A,':D:',F10.3)")"latitude" ,sondes_lat(n)
        end select
      end do
      
      !take Spec_Attributes from any processor with at least one site/sonde
      if(i_Att>0.and.i_Att/=NSPEC)then
        write(*,*)'MISSING species attribute? ',i_Att,NSPEC
      end if
      do d = 1, NPROC-1
        call MPI_RECV(i_Att_MPI, 4*1, MPI_BYTE, d, 746, MPI_COMM_CALC,MPISTATUS, IERROR)
        if(i_Att_MPI>0)then
          if(i_Att_MPI/=NSPEC)then
             write(*,*)'MISSING species attribute? ',i_Att_MPI,NSPEC
          end if
          call MPI_RECV(Spec_Att,Spec_Att_Size*N_Spec_Att_MAX*NSPECMAX, &
               MPI_BYTE, d, 747, MPI_COMM_CALC,MPISTATUS, IERROR)
        end if
      end do

      call CheckStop(NattributesMAX<NSpec_Att+3,'Coords: NattributesMAX too small')
      do i=1,NSPEC
        SpecDef(i,0)=trim(s_species(i)) ! name of the  species
        SpecDef(i,1)='long_name:C:'//trim(s_species(i))
        SpecDef(i,2)='units:C:?'
        SpecDef(i,3)='_FillValue:F:missing'
        do n=1,NSpec_Att
          SpecDef(i,n+3)=trim(Spec_Att(i,n))  ! redefine long_name|units|_FillValue
        end do
      end do
      if(NStations>0)then
         call Create_CDF_sondes(fileName,&
              NSPEC,NSpec_Att+3,SpecDef(:,0:NSpec_Att+3),&
              NStations,NMetaData,MetaData(:,1:NMetaData),&
              NLevels,debug=debug_1d)
         
         write(*,*)'Created ',trim(fileName)
      else
         write(*,*)'No Stations found! not creating ',trim(fileName)
      end if

      deallocate(SpecDef,MetaData)


    else
      !not MasterProc
      i_Att_MPI=i_Att
      call MPI_SEND(i_Att_MPI, 4*1, MPI_BYTE, 0, 746, MPI_COMM_CALC, IERROR)
      if(i_Att>0)then
        call MPI_SEND(Spec_Att, Spec_Att_Size*N_Spec_Att_MAX*NSPECMAX, MPI_BYTE, 0, 747, MPI_COMM_CALC, IERROR)
      end if
      prev_year(type) = current_date%year
    end if ! MasterProc 
  end if ! current_date%year /= prev_year(type)

  if(.not.MasterProc) then   ! send data to me=0 (MasterProc)
    call MPI_SEND(nlocal, 4*1, MPI_BYTE, 0, 346, MPI_COMM_CALC, IERROR)
    call MPI_SEND(out, 8*nout*nlocal, MPI_BYTE, 0, 347, MPI_COMM_CALC, IERROR)
    if(trim(fname)=="sondes")&
        call MPI_SEND(ps_sonde, 8*nlocal, MPI_BYTE, 0, 348, MPI_COMM_CALC, IERROR)
  else ! MasterProc
    ! first, assign me=0 local data to g_out
    if ( DEBUG%SITES ) print *, "ASSIGNS ME=0 NLOCAL_SITES", me, nlocal

    do n = 1, nlocal
      nglob = s_gindex(0,n)
      g_out(:,nglob) = out(:,n)
      if(trim(fname)=="sondes")g_ps(n) = ps_sonde(n)
    end do ! n

    do d = 1, NPROC-1
      call MPI_RECV(nloc, 4*1, MPI_BYTE, d, 346, MPI_COMM_CALC,MPISTATUS, IERROR)
      call MPI_RECV(out, 8*nout*nloc, MPI_BYTE, d, 347, MPI_COMM_CALC, &
           MPISTATUS, IERROR)
      if(trim(fname)=="sondes")call MPI_RECV(ps_sonde, 8*nloc, MPI_BYTE, d, &
           348, MPI_COMM_CALC, MPISTATUS, IERROR)
      do n = 1, nloc
         nglob = s_gindex(d,n)
         g_out(:,nglob) = out(:,n)
         if(trim(fname)=="sondes")g_ps(nglob) = ps_sonde(n)
      end do ! n
    end do ! d

    ! some computers print out e.g. "2.23-123" instead of "2.23e-123"
    ! when numbes get too small. Here we make a correction for this:
    where( abs(g_out)>0.0 .and. abs(g_out)<1.0e-99 ) g_out = 0.0

    ! Final output
    do n = 1, nglobal

      !! Massimo Vieno change the ouput style make the output csv
      !! Oct 2012 Formatting changed (DS,AMV) for gfortran compliance
      !! and compactness. 
      write (io_num,'(a,9999(:,",",es10.3))') & 
           !MAY2019 trim(s_name(n)) // date2string(", DD/MM/YYYY hh:00",current_date),& 
           trim(s_name(n)) // date2string(", DD/MM/YYYY,hh:00",current_date),& 
           ( g_out(ii,n), ii =1, nout ) 
      ! (The ':' format control item will stop processing once the g_out
      !  is done, avoiding runtime warnings.)

    end do

    if(trim(fname)=="sondes")then
      NLevels = KMAX_MID !number of vertical levels (counting from surface)
      NSPEC=NSPC_SONDE!number of species defined for sondes
    else
      NLevels=1
      NSPEC=NSPC_SITE!number of species defined for sites
    end if
    allocate(SpecName(NSPEC))
    do ispec=1,NSPEC
      SpecName(ispec)=trim(s_species(ispec))!name of the variable for one sites/sonde and species          
    end do ! n
    if(nglobal>0)then
       call Out_CDF_sondes(fileName,SpecName,NSPEC,g_out,NLevels,g_ps,debug=debug_1d)
    end if

    deallocate(SpecName)
  end if ! MasterProc
end subroutine siteswrt_out
!==================================================================== >
endmodule Sites_mod
