! <Sites_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2012 met.no
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
module Sites_ml

! -----------------------------------------------------------------------
! Contains subroutines to read in list of measurement stations and/or
! radiosonde locations. Locations should be specified in the input
! files "sites.dat" and "sondes.dat"
!
! -----------------------------------------------------------------------

use CheckStop_ml,  only : CheckStop, StopAll
use My_Outputs_ml, only : &  ! for sitesout
      NSITES_MAX, &
      NADV_SITE, NSHL_SITE, NXTRA_SITE_MISC, NXTRA_SITE_D2D, &
      SITE_ADV, SITE_SHL, SITE_XTRA_MISC, SITE_XTRA_D2D, &
      FREQ_SITE, NSONDES_MAX, NLEVELS_SONDE, &
      NADV_SONDE, NSHL_SONDE, NXTRA_SONDE, & 
      SONDE_ADV, SONDE_SHL, SONDE_XTRA, & 
      FREQ_SONDE, to_ug_ADV

use DerivedFields_ml,  only : f_2d, d_2d  ! not used:, d_3d
use Functions_ml,      only : Tpot_2_T    ! Conversion function
use GridValues_ml,     only : sigma_bnd, sigma_mid, lb2ij, i_fdom, j_fdom &
                              , i_local, j_local, A_mid, B_mid
use Io_ml,             only : check_file,open_file,ios &
                              , fexist, IO_SITES, IO_SONDES &
                              , Read_Headers,read_line
use ChemSpecs_adv_ml
use ChemSpecs_shl_ml,  only : NSPEC_SHL
use ChemGroups_ml,     only : OXN_GROUP, PMFINE_GROUP, PMCO_GROUP
use ChemChemicals_ml,  only : species               ! for species names
use MetFields_ml,      only : t2_nwp, th, pzpbl  &  ! output with concentrations
                              , z_bnd, z_mid, roa, Kz_m2s, q
use MetFields_ml,      only : u_xmj, v_xmi, ps
use ModelConstants_ml, only : NMET,PPBINV,PPTINV, KMAX_MID, MasterProc &
                              ,KMAX_BND,PT,ATWAIR, NPROC, DEBUG_SITES &
                              ,DomainName, RUNDOMAIN, IOU_INST, SOURCE_RECEPTOR
use Par_ml,            only : li0,lj0,li1,lj1 &
                              ,GIMAX,GJMAX &
                              ,GI0,GI1,GJ0,GJ1,me,MAXLIMAX,MAXLJMAX
use SmallUtils_ml,     only : find_index
use Tabulations_ml,    only : tab_esat_Pa
use TimeDate_ml,       only : current_date
use KeyValue_ml,       only : KeyVal, KeyValue, LENKEYVAL

implicit none
private                     ! stops variables being accessed outside


! subroutines made available

public :: sitesdef          ! Calls Init_sites for sites and sondes
public :: siteswrt_surf     ! Gets site  data ready for siteswrt_out
public :: siteswrt_sondes   ! Gets sonde data ready for siteswrt_out
private :: Init_sites       ! reads locations, species
private :: set_species      ! Sets species/variable names for output
private :: siteswrt_out     ! Collects output from all nodes and prints


! some variables used in following subroutines

INCLUDE 'mpif.h'
INTEGER STATUS(MPI_STATUS_SIZE),INFO
integer, private, save :: nglobal_sites, nlocal_sites
integer, private, save :: nglobal_sondes, nlocal_sondes

! site_gindex stores the global index n asociated
! with each processor and local site

integer, private, save, allocatable,dimension (:,:)  :: site_gindex
integer, private, save, allocatable,dimension (:,:) :: sonde_gindex

integer, private, save, dimension (NSITES_MAX) :: &
         site_gx, site_gy, site_gz   & ! global coordinates
       , site_x, site_y, site_z      & ! local coordinates
       , site_n                        ! number in global
integer, private, save, dimension (NSONDES_MAX) ::  &
         sonde_gx, sonde_gy   &        ! global coordinates
       , sonde_x, sonde_y     &        ! local coordinates
       , sonde_n                       ! number in global

! Values from My_Outputs_ml gives ... =>
integer, private, parameter :: & ! Total No., without counting levels
   NSPC_SITE  = NADV_SITE + NSHL_SITE + NXTRA_SITE_MISC + NXTRA_SITE_D2D &
  ,NSPC_SONDE = NADV_SONDE + NSHL_SONDE + NXTRA_SONDE
integer, public, parameter :: & ! Total No., levels included
   NOUT_SITE  = NSPC_SITE * 1 &
  ,NOUT_SONDE = NSPC_SONDE* NLEVELS_SONDE

character(len=50), private, save, dimension(NSITES_MAX) :: site_name
character(len=50), private, save, dimension(NSONDES_MAX):: sonde_name
character(len=20), private, save, dimension(NSPC_SITE)  :: site_species
character(len=20), private, save, dimension(NSPC_SONDE) :: sonde_species

character(len=70), private :: errmsg ! Message text
integer, private :: d                 ! processor index
integer, private :: i, n, nloc, ioerr ! general integers

! Debugging parameter:
logical, private, parameter :: MY_DEBUG = .false.

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

  call Init_sites("sites",IO_SITES,NSITES_MAX, &
        nglobal_sites,nlocal_sites, &
        site_gindex, site_gx, site_gy, site_gz, &
        site_x, site_y, site_z, site_n, &
        site_name)

  call Init_sites("sondes",IO_SONDES,NSONDES_MAX, &
        nglobal_sondes,nlocal_sondes, &
        sonde_gindex, sonde_gx, sonde_gy, sonde_gz, &
        sonde_x, sonde_y, sonde_z, sonde_n, &
        sonde_name)

  call set_species(SITE_ADV,SITE_SHL,SITE_XTRA,site_species)
  call set_species(SONDE_ADV,SONDE_SHL,SONDE_XTRA,sonde_species)

  if ( MY_DEBUG ) then
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
  s(nadv+1:n2) = species( shl(:) )%name
  s(n2+1:nout) = xtra(:)

end  subroutine set_species
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
  character(len=40) :: infile, errmsg
  real              :: lat,lon,x,y

  character(len=20), dimension(4) :: Headers
  type(KeyVal), dimension(20)     :: KeyValues ! Info on units, coords, etc.
  integer                         :: NHeaders, NKeys
  character(len=80)               :: txtinput  ! Big enough to contain
                                               ! one full input record

  ios = 0                      ! zero indicates no errors
  errmsg = ''
  if ( MasterProc ) then
    infile  = fname // ".dat"
    call check_file(infile,fexist,needed=.false.,errmsg=errmsg)
    if ( .not.fexist )then
      ios=1
    else
      call open_file(io_num,"r",infile,needed=.true.)
      call CheckStop(ios,"ios error on "//trim(infile))
    endif
  endif

  call MPI_BCAST( ios, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,INFO)
  if(ios/=0)return

  call CheckStop(NMAX,size(s_name), &
     "Error in Sites_ml/Init_sites: sitesdefNMAX problem")

  call Read_Headers(io_num,errmsg,NHeaders,NKeys,Headers,Keyvalues)

  ! First, see which sites are within full domain:

  n = 0          ! Number of sites found within domain
  SITELOOP: do nin = 1, NMAX

    if (trim(KeyValue(KeyValues,"Coords"))=='LatLong') then
      call read_line(io_num,txtinput,ios)
      if ( ios /= 0 ) exit  ! End of file
      read(unit=txtinput,fmt=*) s, lat, lon, lev
      call lb2ij(lon,lat,x,y)
      ix=nint(x)
      iy=nint(y)
    else
      call read_line(io_num,txtinput,ios)
      if ( ios /= 0 ) exit  ! End of file
      read(unit=txtinput,fmt=*) s,  ix,  iy, lev
    endif

    if (ioerr < 0) then
      write(6,*) "sitesdef : end of file after ", nin-1, infile
      exit SITELOOP
    endif ! ioerr


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

      s_name(n)  = s !!! remove comments// comment
    endif

  enddo SITELOOP

  nglobal = n

  ! NSITES/SONDES_MAX must be _greater_ than the number used, for safety

  call CheckStop(n >= NMAX, &
      "Error in Sites_ml/Init_sites: increaseNGLOBAL_SITES_MAX!")

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

      if (MY_DEBUG) &
        write(6,*) "sitesdef Site on me : ", me, " No. ", nlocal, &
          s_gx(n), s_gy(n) , s_gz(n), " =>  ", &
          s_x(nlocal), s_y(nlocal), s_z(nlocal)

     endif

  enddo ! nglobal

  ! inform me=0 of local array indices:
  if(MY_DEBUG) write(6,*) "sitesdef ", fname, " before gc NLOCAL_SITES", &
                           me, nlocal

  if ( .not.MasterProc ) then
    call MPI_SEND(nlocal, 4*1, MPI_BYTE, 0, 333, MPI_COMM_WORLD, INFO)
    if(nlocal>0) call MPI_SEND(s_n, 4*nlocal, MPI_BYTE, 0, 334, &
                               MPI_COMM_WORLD, INFO)
  else
    if(DEBUG_SITES) write(6,*) "sitesdef for me =0 LOCAL_SITES", me, nlocal
    do n = 1, nlocal
      s_gindex(me,n) = s_n(n)
    enddo
    do d = 1, NPROC-1
      call MPI_RECV(nloc, 4*1, MPI_BYTE, d, 333, MPI_COMM_WORLD,STATUS, INFO)
      if(nloc>0) call MPI_RECV(s_n_recv, 4*nloc, MPI_BYTE, d, 334, &
                               MPI_COMM_WORLD,STATUS, INFO)
      if(DEBUG_SITES) write(6,*) "sitesdef: recv d ", fname, d,  &
                  " zzzz nloc : ", nloc, " zzzz me0 nlocal", nlocal
      do n = 1, nloc
        s_gindex(d,n) = s_n_recv(n)
        if(DEBUG_SITES) write(6,*) "sitesdef: for d =", fname, d, &
          " nloc = ", nloc, " n: ",  n,  " gives nglob ", s_gindex(d,n)
      enddo ! n
    enddo ! d
  endif ! MasterProc

  if ( DEBUG_SITES ) write(6,*) 'sitesdef on me', me, ' = ', nlocal

end subroutine Init_sites
!==================================================================== >
subroutine siteswrt_surf(xn_adv,cfac,xn_shl)
  ! -------------------------------------------------------------------
  ! writes out just simple concentrations for now....
  ! will be improved later to allow choice of output parameter
  ! should look at chemint also - seems similar for somethings
  ! -------------------------------------------------------------------

  ! arguments
  real, dimension(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID), intent(in) :: xn_adv
  real, dimension(NSPEC_ADV,MAXLIMAX,MAXLJMAX), intent(in)          :: cfac
  real, dimension(NSPEC_SHL,MAXLIMAX,MAXLJMAX,KMAX_MID), intent(in) :: xn_shl

  ! Local
  integer :: ix, iy,iz, ispec                  ! Site indices
  integer :: nn                                ! species index
  logical, save :: my_first_call = .true.      ! for debugging
  integer                           :: d2index ! index for d_2d field access
  character(len=len(SITE_XTRA_D2D)) :: d2code  ! parameter code -- # --

  real,dimension(NOUT_SITE,NSITES_MAX) :: out  ! for output, local node

  if ( DEBUG_SITES ) then
    write(6,*) "sitesdef Into surf  nlocal ", nlocal_sites, " on me ", me
    do i = 1, nlocal_sites
      write(6,*) "sitesdef Into surf  x,y ",site_x(i),site_y(i),&
                  site_z(i)," me ", me
    enddo

    if ( MasterProc ) then
      write(6,*) "======= site_gindex ======== sitesdef ============"
      do n = 1, nglobal_sites
        write(6,'(a12,i4,2x,80(i4,:))') "sitesdef ", n, &
                (site_gindex(d,n),d=0,NPROC-1)
      enddo
      write(6,*) "======= site_end    ======== sitesdef ============"
    endif ! MasterProc
  endif ! DEBUG_SITES

  ! assign local data to out

  do i = 1, nlocal_sites
    ix = site_x(i)
    iy = site_y(i)
    iz = site_z(i)

    do ispec = 1, NADV_SITE
      if (iz == KMAX_MID ) then ! corrected to surface
        out(ispec,i) = xn_adv( SITE_ADV(ispec) ,ix,iy,KMAX_MID ) * &
                       cfac( SITE_ADV(ispec),ix,iy) * PPBINV
      else                      ! Mountain sites not corrected to surface
        out(ispec,i)  = xn_adv( SITE_ADV(ispec) ,ix,iy,iz ) * PPBINV
      endif
    enddo


    do ispec = 1, NSHL_SITE
      out(NADV_SITE+ispec,i)  = xn_shl( SITE_SHL(ispec) ,ix,iy,iz )
    end do

    ! XTRA parameters, usually the temmp or pressure
    if (.not. SOURCE_RECEPTOR) then
      nn = NADV_SITE + NSHL_SITE
      do ispec = 1, NXTRA_SITE_MISC
        nn=nn+1
        select case ( trim(SITE_XTRA_MISC(ispec)) )
          case ( "T2" )
            out(nn,i)   = t2_nwp(ix,iy,1) - 273.15
          case ( "th" )
            out(nn,i)   = th(ix,iy,iz,1)
!         case ( "hmix" )
!           out(nn,i)   = pzpbl(ix,iy)
         case default
           call CheckStop("Error, Sites_ml/siteswrt_surf: SITE_XTRA_MISC:"&
                               // trim(SITE_XTRA_MISC(ispec)))
        end select
        call CheckStop( abs(out(nn,i))>1.0e99, &
          "ABS(SITES OUT: '"//trim(SITE_XTRA_MISC(ispec))//"') TOO BIG" )
      end do
      do ispec = 1, NXTRA_SITE_D2D
        d2code      = SITE_XTRA_D2D(ispec)
        d2index     = find_index(d2code, f_2d(:)%name)
        nn=nn+1

        if ( d2index < 1 ) then
           if( my_first_call) write(*,*) &
                 "WARNING: SITES D2D NOT FOUND"//trim(d2code)
           !cycle
           out(nn,i)   = -999.9
        else
           !call CheckStop( d2index<1, "SITES D2D NOT FOUND"//trim(d2code) )
           out(nn,i)   = d_2d(d2index,ix,iy,IOU_INST)
        end if

        !May25 call CheckStop( d2index<1, "SITES D2D NOT FOUND"//trim(d2code) )
        !May25 out(nn,i)   = d_2d(d2index,ix,iy,IOU_INST)
        if( DEBUG_SITES ) &
          write(6,"(a,3i3,a,i4,es10.3)") "DEBUG_SITES ", me, nn, i,&
            trim(d2code), d2index, out(nn,i)
        call CheckStop( abs(out(nn,i))>1.0e99, &
          "ABS(SITES OUT: '"//trim(SITE_XTRA_D2D(ispec))//"') TOO BIG" )
      enddo
    endif
  enddo

    my_first_call = .false.

  ! collect data into gout on me=0 t

  call siteswrt_out("sites",IO_SITES,NOUT_SITE, FREQ_SITE, &
                     nglobal_sites,nlocal_sites, &
                     site_gindex,site_name,site_gx,site_gy,site_gz,&
                     site_species,out)

end subroutine siteswrt_surf
!==================================================================== >
subroutine siteswrt_sondes(xn_adv,xn_shl)
  ! -------------------------------------------------------------------
  ! Writes vertical concentration  data to files.
  ! IO_SONDES is set in io_ml to be 30
  ! -------------------------------------------------------------------

  real, dimension(:,:,:,:), intent(in) ::  xn_adv
  real, dimension(:,:,:,:), intent(in) ::  xn_shl

  ! Output variables - none

  ! Local variables
  integer :: n, i, k,  ix, iy, nn, ispec   ! Site and chem indices
  integer, parameter ::  KTOP_SONDE = KMAX_MID - NLEVELS_SONDE + 1
  integer, dimension(NLEVELS_SONDE)      :: itemp
  real, dimension(KMAX_MID)              :: pp, temp, qsat, rh, sum_PM, sum_NOy
  real, dimension(NOUT_SONDE,NSONDES_MAX):: out

  ! Consistency check

  do ispec = 1, NXTRA_SONDE
    select case ( SONDE_XTRA(ispec) )
      case ( "PM25 " ,"PMco " , "NOy ", "RH ","z_mid", "p_mid", &
            "Kz_m2s ", "th   ", "U   ", "V    " )
              errmsg = "ok"
      case default
        call CheckStop("Error, Sites_ml/siteswrt_sondes: SONDE_XTRA:"&
                              // SONDE_XTRA(ispec))
    end select
  enddo

  do i = 1, nlocal_sondes
    n  = sonde_n(i)
    ix = sonde_x(i)
    iy = sonde_y(i)
    nn = 0

    ! collect and print out with ground-level (KMAX_MID) first, hence &
    ! KMAX_MID:KTOP_SONDE:-1 in arrays
    ! first the advected and short-lived species

    do ispec = 1, NADV_SONDE    !/ xn_adv in ppb
      out(nn+1:nn+NLEVELS_SONDE,i) = PPBINV *  &
          xn_adv( SONDE_ADV(ispec) , ix,iy,KMAX_MID:KTOP_SONDE:-1)
      nn = nn + NLEVELS_SONDE
    enddo

    do ispec = 1, NSHL_SONDE    !/ xn_shl  in molecules/cm3
      out(nn+1:nn+NLEVELS_SONDE,i) = xn_shl( SONDE_SHL(ispec) , &
          ix,iy,KMAX_MID:KTOP_SONDE:-1)
      nn = nn + NLEVELS_SONDE
    enddo

    ! then print out XTRA stuff first,
    ! usually the height or pressure

    do ispec = 1, NXTRA_SONDE

      select case ( SONDE_XTRA(ispec) )
        case ( "PM25" )    !!  PM data converted to ug m-3
          sum_PM(:) = 0.
          do k = 1, KMAX_MID
            sum_PM(k) = (dot_product(xn_adv(PMFINE_GROUP-NSPEC_SHL,ix,iy,k),&
                                  to_ug_ADV(PMFINE_GROUP-NSPEC_SHL)) &
                          + 0.5 * &  ! 50% of PMcoare in PM2.5, since Dp=2.5
                          dot_product(xn_adv(PMCO_GROUP-NSPEC_SHL,ix,iy,k),&
                                  to_ug_ADV(PMCO_GROUP-NSPEC_SHL)) &
                        ) * roa(ix,iy,k,1)
          enddo !k
          out(nn+1:nn+KMAX_MID,i) = sum_PM(KMAX_MID:1:-1)

        case ( "PMco" ) !!  PM data converted to ug m-3
          sum_PM(:) = 0.
          do k = 1, KMAX_MID
            sum_PM(k) = dot_product(xn_adv(PMCO_GROUP-NSPEC_SHL,ix,iy,k),&
                                  to_ug_ADV(PMCO_GROUP-NSPEC_SHL)) &
                      * roa(ix,iy,k,1)
          enddo !k
          out(nn+1:nn+KMAX_MID,i) = sum_PM(KMAX_MID:1:-1)

        case ( "NOy" )
          sum_NOy(:) = 0.
          do k = 1, KMAX_MID
            sum_NOy(k) = sum(xn_adv(OXN_GROUP-NSPEC_SHL,ix,iy,k))
          enddo
          out(nn+1:nn+KMAX_MID,i) = PPBINV * sum_NOy(KMAX_MID:1:-1)

        case ( "RH   " )
          do k = 1,KMAX_MID
            pp(k) = A_mid(k) + B_mid(k)*ps(ix,iy,1)
            temp(k) = th(ix,iy,k,1)* Tpot_2_T( pp(k) )
            itemp(k) = nint( temp(k) )
            qsat(k)  = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
            rh(k) = min( q(ix,iy,k,1)/qsat(k) , 1.0)
          end do
          out(nn+1:nn+NLEVELS_SONDE,i) =  rh(KMAX_MID:KTOP_SONDE:-1)

        case ( "z_mid" )
          out(nn+1:nn+NLEVELS_SONDE,i) =  z_mid(ix,iy,KMAX_MID:KTOP_SONDE:-1)

        case ( "p_mid" )
          out(nn+1:nn+NLEVELS_SONDE,i) = A_mid(KMAX_MID:KTOP_SONDE:-1) + &
                                    B_mid(KMAX_MID:KTOP_SONDE:-1)*ps(ix,iy,1)

        case ( "Kz_m2s" )
          out(nn+1:nn+NLEVELS_SONDE,i) =  Kz_m2s(ix,iy,KMAX_MID:KTOP_SONDE:-1)

        case ( "th" )
          out(nn+1:nn+NLEVELS_SONDE,i) =  th(ix,iy,KMAX_MID:KTOP_SONDE:-1,1)

        case ( "U" )
          out(nn+1:nn+NLEVELS_SONDE,i) = 0.5 &
             *( u_xmj(ix,iy,KMAX_MID:KTOP_SONDE:-1,1) &
             +u_xmj(ix-1,iy,KMAX_MID:KTOP_SONDE:-1,1) )

        case ( "V" )
          out(nn+1:nn+NLEVELS_SONDE,i) = 0.5 &
             *( v_xmi(ix,iy,KMAX_MID:KTOP_SONDE:-1,1) &
             +v_xmi(ix,iy-1,KMAX_MID:KTOP_SONDE:-1,1) )

        case ( "D3D" )
            call StopAll("D3D Sites out not defined")

      end select

      nn = nn + NLEVELS_SONDE

    enddo ! ispec (NXTRA_SONDE)

  enddo ! i (nlocal_sondes)


  ! collect data into gout on me=0 t

  call siteswrt_out("sondes",IO_SONDES,NOUT_SONDE, FREQ_SONDE, &
                     nglobal_sondes,nlocal_sondes, &
                     sonde_gindex,sonde_name,sonde_gx,sonde_gy,sonde_gy, &
                     sonde_species,out)

end subroutine siteswrt_sondes
!==================================================================== >
subroutine siteswrt_out(fname,io_num,nout,f,nglobal,nlocal, &
                        s_gindex,s_name,s_gx,s_gy,s_gz,s_species,out)
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

  ! Local
  real,dimension(nout,nglobal) :: g_out ! for output, collected
  integer :: nglob, nloc        ! Site indices
  character(len=40)  :: outfile
  character(len=4)   :: suffix
  integer, parameter :: NTYPES = 2      ! No. types, now 2 (sites, sondes)
  integer ::  type=-1                   ! = 1 for sites, 2 for sondes
  integer, save, dimension(NTYPES):: prev_month = (/ -99, -99 /) ! Initialise
  integer, save, dimension(NTYPES):: prev_year = (/ -99, -99 /) ! Initialise

  select case (fname)
    case ("sites" )
      type = 1
    case ("sondes" )
      type = 2
    case default
      write(6,*) "non-possible tpye in siteswrt_out for ", fname
      return
  end select

!  if ( MasterProc .and. current_date%month /= prev_month(type)) then
  if (MasterProc .and. current_date%year /= prev_year(type) ) then

     if ( prev_year(type) > 0 ) close(io_num)  ! Close last-months file
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
     write(io_num,'(a9,<size(s_species)>(",",a))')"site,date",(trim(s_species(i)),i=1,size(s_species))

  endif ! first call



  if ( .not.MasterProc ) then   ! send data to me=0 (MasterProc)

    call MPI_SEND(nlocal, 4*1, MPI_BYTE, 0, 346, MPI_COMM_WORLD, INFO)
    call MPI_SEND(out, 8*nout*nlocal, MPI_BYTE, 0, 347, MPI_COMM_WORLD, INFO)

  else ! MasterProc

    ! first, assign me=0 local data to g_out
    if ( DEBUG_SITES ) print *, "ASSIGNS ME=0 NLOCAL_SITES", me, nlocal

    do n = 1, nlocal
      nglob = s_gindex(0,n)
      g_out(:,nglob) = out(:,n)
    enddo ! n

    do d = 1, NPROC-1
      call MPI_RECV(nloc, 4*1, MPI_BYTE, d, 346, MPI_COMM_WORLD,STATUS, INFO)
      call MPI_RECV(out, 8*nout*nloc, MPI_BYTE, d, 347, MPI_COMM_WORLD, &
                    STATUS, INFO)
      do n = 1, nloc
        nglob = s_gindex(d,n)
        g_out(:,nglob) = out(:,n)
      enddo ! n
    enddo ! d

    ! some computers print out e.g. "2.23-123" instead of "2.23e-123"
    ! when numbes get too small. Here we make a correction for this:
    where( abs(g_out)>0.0 .and. abs(g_out)<1.0e-99 ) g_out = 0.0

    ! Final output
    do n = 1, nglobal

!! Massimo Vieno change the ouput style make the output csv
         write (io_num,'(a<len(trim(s_name(n)))>,",",i2.2,"/",i2.2,"/",i4.4," ",i2.2,":00",<size(s_species)>(",",es10.3))') &
               trim(s_name(n)),&
               current_date%day,current_date%month,current_date%year,current_date%hour,g_out(:,n)
    enddo ! n

  endif ! MasterProc

end subroutine siteswrt_out
!==================================================================== >
end module Sites_ml
