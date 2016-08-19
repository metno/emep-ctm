! <Sites_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Sites_ml

! -----------------------------------------------------------------------
! Contains subroutines to read in list of measurement stations and/or 
! radiosonde locations. Locations should be specified in the input 
! files "sites.dat" and "sondes.dat"
!
! -----------------------------------------------------------------------

use CheckStop_ml,  only : CheckStop
use My_Outputs_ml, only : &  ! for sitesout
      NSITES_MAX, &
      NADV_SITE, NSHL_SITE, NXTRA_SITE, &
      SITE_ADV, SITE_SHL, SITE_XTRA, SITE_XTRA_INDEX, &
      FREQ_SITE, NSONDES_MAX, NLEVELS_SONDE, &
      NADV_SONDE, NSHL_SONDE, NXTRA_SONDE, N_NOy,  &
      SONDE_ADV, SONDE_SHL, SONDE_XTRA, SONDE_XTRA_INDEX, &
      FREQ_SONDE, NOy_SPEC

use Derived_ml,        only : d_2d, d_3d, IOU_INST  ! for deriv system
use Functions_ml,      only : Tpot_2_T              ! Conversion function
use GridValues_ml,     only : sigma_bnd, sigma_mid, lb2ij, i_fdom, j_fdom &
                              , i_local, j_local
use Io_ml,             only : check_file,open_file,ios &
                              , fexist, IO_SITES, IO_SONDES &
                              , Read_Headers,read_line
use GenSpec_adv_ml
use GenSpec_shl_ml,    only : NSPEC_SHL
use GenChemicals_ml,   only : species               ! for species names
use GenSpec_tot_ml,    only : SO4, HCHO, CH3CHO  &  ! for mol. wts.
                              ,aNO3, pNO3, aNH4, PM25, PMCO &
                              ,SSfi, SSco  !SeaS
use Met_ml,            only : t2_nwp, th, pzpbl  &  ! Output with concentrations
                              , z_bnd, z_mid, roa, xksig, u, v, ps, q
use ModelConstants_ml, only : NMET,PPBINV,PPTINV, KMAX_MID &
                              ,KMAX_BND,PT,ATWAIR, NPROC &
                              ,DomainName, RUNDOMAIN
use Par_ml,            only : li0,lj0,li1,lj1 &
                              ,GIMAX,GJMAX &
                              ,GI0,GI1,GJ0,GJ1,me,MAXLIMAX,MAXLJMAX
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

integer, private, save, dimension (0:NPROC-1,NSITES_MAX) :: site_gindex
integer, private, save, dimension (0:NPROC-1,NSONDES_MAX) :: sonde_gindex

integer, private, save, dimension (NSITES_MAX) :: &
         site_gx, site_gy, site_gz   & ! global coordinates
       , site_x, site_y, site_z      & ! local coordinates
       , site_n                        ! number in global
integer, private, save, dimension (NSONDES_MAX) ::  &
         sonde_gx, sonde_gy   &        ! global coordinates
       , sonde_x, sonde_y     &        ! local coordinates
       , sonde_n                       ! number in global

! Values from My_Outputs_ml gives ... =>

integer, public, parameter :: &
     NOUT_SITE = NADV_SITE + NSHL_SITE + NXTRA_SITE & ! Total No.
    ,NOUT_SONDE = NLEVELS_SONDE * ( NADV_SONDE + NSHL_SONDE+ NXTRA_SONDE )

character(len=50), private, save, dimension (NSITES_MAX) :: site_name
character(len=50), private, save, dimension (NSONDES_MAX):: sonde_name
character(len=20), private, save, &
              dimension (NADV_SITE+NSHL_SITE+NXTRA_SITE) :: site_species
character(len=20), private, save, &
           dimension (NADV_SONDE+NSHL_SONDE+NXTRA_SONDE) :: sonde_species

character(len=40), private :: errmsg ! Message text
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
  sonde_gz(:) = 0
  sonde_z(:)  = 0

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
  character(len=*), intent(in), &
                         dimension(:) :: xtra  !Names of extra params
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
                           , s_x, s_y, s_z      & ! local coordinates
                           , s_n                  ! number in global
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
  integer                         :: NHeaders, NKeys, Nlines
  character(len=80)               :: txtinput  ! Big enough to contain
                                               ! one full input record

  ios = 0                      ! zero indicates no errors
  errmsg = "ios error" // infile

  if (me == 0) then
    infile  = fname // ".dat"
    call check_file(infile,fexist,needed=.false.,errmsg=errmsg)
    if ( .not. fexist ) return
    call open_file(io_num,"r",infile,needed=.true.)
    call CheckStop(ios,"ios error on "//trim(infile))
  end if

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
     end if ! ioerr


     if ( ix < RUNDOMAIN(1) .or. ix > RUNDOMAIN(2) .or. & 
          iy < RUNDOMAIN(3) .or. iy > RUNDOMAIN(4) ) then
        if(me==0) write(6,*) "sitesdef: ", s, ix, iy, &
                " outside computational domain"
     else if ( ix == RUNDOMAIN(1) .or. ix == RUNDOMAIN(2) .or. & 
              iy == RUNDOMAIN(3) .or. iy == RUNDOMAIN(4) ) then
        if(me==0) write(6,*) "sitesdef: ", s, ix, iy, &
                " on computational domain"
     else
        comment = " ok - inside domain         "
        n = n + 1
        s_gx(n)   = ix  
        s_gy(n)   = iy  
        s_gz(n)   = lev

        s_name(n)  = s // comment

     endif

  end do SITELOOP 

  nglobal = n 

  ! NSITES/SONDES_MAX must be _greater_ than the number used, for safety

  call CheckStop(n >= NMAX, &
      "Error in Sites_ml/Init_sites: increaseNGLOBAL_SITES_MAX!")

  if(me==0) close(unit=io_num)

  nlocal  = 0

  do n = 1, nglobal

     ix = s_gx(n) ! global-domain coords
     iy = s_gy(n)

     if ( i_local(ix) >= li0 .and. i_local(ix) <= li1 .and. &
          j_local(iy) >= lj0 .and. j_local(iy) <= lj1 ) then

        nlocal      = nlocal + 1
        s_x(nlocal) = i_local(ix)
        s_y(nlocal) = j_local(iy)
        s_z(nlocal) = s_gz(n)
        s_n(nlocal) = n

        if (MY_DEBUG) then
           write(6,*) "sitesdef Site on me : ", me, " No. ", nlocal, &
                s_gx(n), s_gy(n) , s_gz(n), " =>  ", &
                s_x(nlocal), s_y(nlocal), s_z(nlocal)
        end if

     endif

  end do ! nglobal

  ! inform me=0 of local array indices:
  if (MY_DEBUG) write(6,*) "sitesdef ", fname, " before gc NLOCAL_SITES", &
                           me, nlocal

  if ( me /= 0 ) then
     CALL MPI_SEND(nlocal, 4*1, MPI_BYTE, 0, 333, MPI_COMM_WORLD, INFO)
     if (nlocal > 0) &
        CALL MPI_SEND(s_n, 4*nlocal, MPI_BYTE, 0, 334, MPI_COMM_WORLD, INFO)
  else
     if (MY_DEBUG) write(6,*) "sitesdef for me =0 OCAL_SITES", me, nlocal
     do n = 1, nlocal
        s_gindex(me,n) = s_n(n)
     end do
     do d = 1, NPROC-1
        CALL MPI_RECV(nloc, 4*1, MPI_BYTE, d, &
            333, MPI_COMM_WORLD,STATUS, INFO) 
        if (nloc > 0) &
           CALL MPI_RECV(s_n_recv, 4*nloc, MPI_BYTE, d, &
             334, MPI_COMM_WORLD,STATUS, INFO) 
        if (MY_DEBUG) write(6,*) "sitesdef: recv d ", fname, d,  &
                   " zzzz nloc : ", nloc, " zzzz me0 nlocal", nlocal
        do n = 1, nloc
           s_gindex(d,n) = s_n_recv(n)
           if ( MY_DEBUG ) write(6,*) "sitesdef: for d =", fname, d, &
              " nloc = ", nloc, " n: ",  n,  " gives nglob ", s_gindex(d,n)
        end do ! n
     end do ! d
  end if ! me

  if ( MY_DEBUG ) write(6,*) 'sitesdef on me', me, ' = ', nlocal

end subroutine Init_sites


!==================================================================== >
  subroutine siteswrt_surf(xn_adv,cfac,xn_shl)

  ! ---------------------------------------------------------------------
  ! writes out just simple concentrations for now....
  ! will be improved later to allow choice of output parameter
  ! should look at chemint also - seems similar for somethings
  ! ---------------------------------------------------------------------

  ! arguments
  real, dimension(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID), intent(in) :: xn_adv
  real, dimension(NSPEC_ADV,MAXLIMAX,MAXLJMAX), intent(in)          :: cfac
  real, dimension(NSPEC_SHL,MAXLIMAX,MAXLJMAX,KMAX_MID), intent(in) :: xn_shl

  ! Local
  integer :: nglob, nloc, ix, iy,iz, ispec    ! Site indices
  integer :: nn                               ! species index
  logical, save :: my_first_call = .true.     ! for debugging
  integer :: d2index                          ! index for d_2d field access

  real,dimension(NOUT_SITE,NSITES_MAX) :: out ! for output, local node

  if ( MY_DEBUG ) then 
     write(6,*) "sitesdef Into surf  nlocal ", nlocal_sites, " on me ", me 
     do i = 1, nlocal_sites
        write(6,*) "sitesdef Into surf  x,y ",site_x(i),site_y(i),&
                    site_z(i)," me ", me
     end do
 
     if ( me == 0 ) then
        write(6,*)  "======= site_gindex ======== sitesdef ============"
        do n = 1, nglobal_sites 
           write(6,'(a12,i4, 2x, 8i4)') "sitesdef ", n, &
                    (site_gindex(d,n),d=0, NPROC-1)  
        end do
        write(6,*) "======= site_end    ======== sitesdef ============"
     end if ! me = 0
  end if ! MY_DEBUG

  ! assign local data to out

  do i = 1, nlocal_sites
     ix = site_x(i)
     iy = site_y(i)
     iz = site_z(i)

     do ispec = 1, NADV_SITE
        if (iz == KMAX_MID ) then ! corrected to surface
           out(ispec,i) = xn_adv( SITE_ADV(ispec) ,ix,iy,KMAX_MID ) * &
                          cfac( SITE_ADV(ispec),ix,iy) * PPBINV
        else                     ! Mountain sites not corrected to surface
           out(ispec,i)  = xn_adv( SITE_ADV(ispec) ,ix,iy,iz ) * PPBINV
        end if
     end do

     my_first_call = .false.

     do ispec = 1, NSHL_SITE
        out(NADV_SITE+ispec,i)  = xn_shl( SITE_SHL(ispec) ,ix,iy,iz )
     end do

    ! then print out XTRA stuff, usually the temmp
    ! or pressure
    ! SITE_XTRA=(/ "th  ", "hmix", "Vg_ref", "Vg_1m", "Vg_sto", "Flux_ref", "Flux_sto" /)

     do ispec = 1, NXTRA_SITE
        nn = NADV_SITE + NSHL_SITE + ispec
        select case ( SITE_XTRA(ispec) )
          case ( "T2" ) 
            out(nn,i)   = t2_nwp(ix,iy,1) - 273.15 
          case ( "th" ) 
            out(nn,i)   = th(ix,iy,iz,1)
          case ( "hmix" ) 
            out(nn,i)   = pzpbl(ix,iy)

          case ( "D2D" )    
            d2index     = SITE_XTRA_INDEX(ispec)
            out(nn,i)   = d_2d(d2index,ix,iy,IOU_INST)
        end select 
        call CheckStop( abs( out(nn,i) )  > 1.0e99, "ABS(SITES OUT) TOO BIG" )
     end do
  end do

  ! collect data into gout on me=0 t

  call siteswrt_out("sites",IO_SITES,NOUT_SITE,NSITES_MAX, FREQ_SITE, &
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

  real, dimension(:,:,:,:),    intent(in) ::  xn_adv
  real, dimension(:,:,:,:),    intent(in) ::  xn_shl
 
  ! Output variables - none

  ! Local variables
  integer :: n, i, ii, k,  ix, iy, nn, ispec   ! Site and chem indices
  integer :: d3index    !ds index for d_3d field access
  integer, parameter ::  KTOP_SONDE = KMAX_MID - NLEVELS_SONDE + 1
  integer, dimension(NLEVELS_SONDE)       :: itemp
  real, dimension(KMAX_MID)               :: pp, temp, qsat, rh, sum_PM25 &
                                             , sum_PMco, sum_NOy
  real, dimension(NOUT_SONDE,NSONDES_MAX) :: out
 
  ! Consistency check 

  do ispec = 1, NXTRA_SONDE
     select case ( SONDE_XTRA(ispec) )
       case ( "PM25 " ,"PMco " , "NOy ", "RH ","z_mid", "p_mid", &
              "xksig ", "th   ", "U   ", "V    " )
               errmsg = "ok"
       case default
         call CheckStop("Error, Sites_ml/siteswrt_sondes: SONDE_XTRA:"&
                               // SONDE_XTRA(ispec))
     end select
  end do

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
     end do

     do ispec = 1, NSHL_SONDE    !/ xn_shl  in molecules/cm3 
        out(nn+1:nn+NLEVELS_SONDE,i) = xn_shl( SONDE_SHL(ispec) , &
            ix,iy,KMAX_MID:KTOP_SONDE:-1)
        nn = nn + NLEVELS_SONDE
     end do

     ! then print out XTRA stuff first,
     ! usually the height or pressure

     do ispec = 1, NXTRA_SONDE

        select case ( SONDE_XTRA(ispec) )

          case ( "PM25" )    !!  PM data converted to ug m-3 
            sum_PM25(:) = 0.
            do k = 1, KMAX_MID
              sum_PM25(k) = &
               ( xn_adv(IXADV_SO4,ix,iy,k) *species(SO4)%molwt    &
               + xn_adv(IXADV_aNO3,ix,iy,k)*species(aNO3)%molwt   &
               + xn_adv(IXADV_aNH4,ix,iy,k)*species(aNH4)%molwt   &
               + xn_adv(IXADV_PM25,ix,iy,k)*species(PM25)%molwt   & 
               + xn_adv(IXADV_SSfi,ix,iy,k)*species(SSfi)%molwt  ) & !SeaS
               * roa(ix,iy,k,1) /ATWAIR
            end do 
            out(nn+1:nn+KMAX_MID,i) = PPBINV                         &
                                    * sum_PM25(KMAX_MID:1:-1)

          case ( "PMco" ) !!  PM data converted to ug m-3 
            sum_PMco(:) = 0.
            do k = 1, KMAX_MID
              sum_PMco(k) =   &
               ( xn_adv(IXADV_pNO3,ix,iy,k) * species(pNO3)%molwt &
               + xn_adv(IXADV_PMco,ix,iy,k) * species(PMCO)%molwt & 
               + xn_adv(IXADV_SSco,ix,iy,k) * species(SSco)%molwt ) & !SeaS
               * roa(ix,iy,k,1) /ATWAIR
            end do
            out(nn+1:nn+KMAX_MID,i) = PPBINV                        &
                                    * sum_PMco(KMAX_MID:1:-1)

          case ( "NOy" )
            sum_NOy(:) = 0.
            do k = 1, KMAX_MID
              do ii = 1, N_NOy
                sum_NOy(k) = sum_NOy(k) + xn_adv(NOy_SPEC(ii),ix,iy,k)
              end do
            end do 
            out(nn+1:nn+KMAX_MID,i) = PPBINV &
                                    * sum_NOy(KMAX_MID:1:-1)

          case ( "RH   " ) 
            do k = 1,KMAX_MID
              pp(k) = PT + sigma_mid(k)*(ps(ix,iy,1) - PT)
              temp(k) = th(ix,iy,k,1)* Tpot_2_T( pp(k) )
              itemp(k) = nint( temp(k) )
              qsat(k)  = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
              rh(k) = min( q(ix,iy,k,1)/qsat(k) , 1.0) 
            end do
            out(nn+1:nn+NLEVELS_SONDE,i) =  rh(KMAX_MID:KTOP_SONDE:-1)

          case ( "z_mid" ) 
            out(nn+1:nn+NLEVELS_SONDE,i) =  z_mid(ix,iy,KMAX_MID:KTOP_SONDE:-1)

          case ( "p_mid" ) 
            out(nn+1:nn+NLEVELS_SONDE,i) = PT + sigma_mid(KMAX_MID:KTOP_SONDE:-1) &
                                                 *(ps(ix,iy,1) - PT)

          case ( "xksig" ) 
            out(nn+1:nn+NLEVELS_SONDE,i) =  xksig(ix,iy,KMAX_MID:KTOP_SONDE:-1)

          case ( "th" ) 
            out(nn+1:nn+NLEVELS_SONDE,i) =  th(ix,iy,KMAX_MID:KTOP_SONDE:-1,1)

          case ( "U" ) 
            out(nn+1:nn+NLEVELS_SONDE,i) = 0.5*( u(ix,iy,KMAX_MID:KTOP_SONDE:-1,1) &
                                            +u(ix-1,iy,KMAX_MID:KTOP_SONDE:-1,1) )

          case ( "V" ) 
            out(nn+1:nn+NLEVELS_SONDE,i) = 0.5*( v(ix,iy,KMAX_MID:KTOP_SONDE:-1,1) &
                                            +v(ix,iy-1,KMAX_MID:KTOP_SONDE:-1,1) )

          case ( "D3D" )    
            d3index                 = SONDE_XTRA_INDEX(ispec)
            out(nn+1:nn+NLEVELS_SONDE,i)= &
              d_3d(d3index,ix,iy,KMAX_MID:KTOP_SONDE:-1,IOU_INST)
 
        end select

        nn = nn + NLEVELS_SONDE

     end do ! ispec (NXTRA_SONDE)

  end do ! i (nlocal_sondes)


  ! collect data into gout on me=0 t

  call siteswrt_out("sondes",IO_SONDES,NOUT_SONDE,NSONDES_MAX, FREQ_SONDE, &
                     nglobal_sondes,nlocal_sondes, &
                     sonde_gindex,sonde_name,sonde_gx,sonde_gy,sonde_gy, &
                     sonde_species,out)

end subroutine siteswrt_sondes


!==================================================================== >
  subroutine siteswrt_out(fname,io_num,nout,nsites,f,nglobal,nlocal, &
                          s_gindex,s_name,s_gx,s_gy,s_gz,s_species,out)

  ! -------------------------------------------------------------------
  ! collects data from local nodes and writes out to sites/sondes.dat
  ! -------------------------------------------------------------------

  character(len=*), intent(in) :: fname
  integer, intent(in) :: io_num, nout, nsites
  integer, intent(in) :: f               ! Frequency of write-out (hours)
  integer, intent(in) :: nglobal, nlocal
  integer, intent(in), dimension (0:,:) :: s_gindex  ! index, starts at me=0
  character(len=*), intent(in), dimension (:) ::  s_name    ! site/sonde name
  integer, intent(in), dimension (:) :: s_gx, s_gy, s_gz    ! coordinates
  character(len=*), intent(in), dimension (:) ::  s_species ! Variable names
  real,    intent(in), dimension(:,:) :: out    ! outputs, local node

  ! Local
  real,dimension(nout,nglobal) :: g_out ! for output, collected
  integer :: nglob, nloc, ix, iy        ! Site indices
  character(len=40)  :: outfile
  character(len=4)   :: suffix
  integer, parameter :: NTYPES = 2      ! No. types, now 2 (sites, sondes)
  integer ::  type                      ! = 1 for sites, 2 for sondes
  integer, save, dimension(NTYPES):: prev_month = (/ -99, -99 /) ! Initialise

  select case (fname)
  case ("sites" )
     type = 1
  case ("sondes" )
     type = 2
  case default
     write(6,*) "non-possible tpye in siteswrt_out for ", fname
     return
  end select

  if (me == 0 .and. current_date%month /= prev_month(type) ) then

     if ( prev_month(type) > 0 ) close(io_num)  ! Close last-months file

     ! Open new file for write-out

     write(suffix,fmt="(2i2.2)") current_date%month, modulo ( current_date%year, 100 )
     outfile = fname // "." // suffix
     open(file=outfile,unit=io_num,action="write")
     prev_month(type) = current_date%month

     write(io_num,"(i3,2x,a,a, 4i4)") nglobal, fname, " in domain",  &
                                 RUNDOMAIN
     write(io_num,"(i3,a)") f, " Hours between outputs"

     do n = 1, nglobal
        write(io_num,"(a50,3i4)") s_name(n), s_gx(n), s_gy(n),s_gz(n)
     end do ! nglobal

     write(io_num,"(i3,a)") size(s_species), " Variables:"
     do n = 1, size(s_species)
        write(io_num,"(i3,2x,a)") n, s_species(n)
     end do

  endif ! New month

  if ( me /= 0 ) then   ! send data to me=0

     CALL MPI_SEND(nlocal, 4*1, MPI_BYTE, 0, 346, MPI_COMM_WORLD, INFO) 
     CALL MPI_SEND(out, 8*nout*nlocal, MPI_BYTE, 0, 347, MPI_COMM_WORLD, INFO) 

  else ! me = 0

     ! first, assign me=0 local data to g_out
     if ( MY_DEBUG ) print *, "ASSIGNS ME=0 NLOCAL_SITES", me, nlocal

     do n = 1, nlocal
        nglob = s_gindex(0,n)
        g_out(:,nglob) = out(:,n)
     end do ! n

     do d = 1, NPROC-1
        CALL MPI_RECV(nloc, 4*1, MPI_BYTE, d, &
           346, MPI_COMM_WORLD,STATUS, INFO) 
        CALL MPI_RECV(out, 8*nout*nloc, MPI_BYTE, d, 347, MPI_COMM_WORLD, STATUS, INFO) 
        do n = 1, nloc
           nglob = s_gindex(d,n)
           g_out(:,nglob) = out(:,n)
        end do ! n
     end do ! d

     ! some computers print out e.g. "2.23-123" instead of "2.23e-123"
     ! when numbes get too small. Here we make a correction for this:
     where( abs(g_out) > 0.0 .and. abs(g_out) < 1.0e-99 )
        g_out = 0.0
     end where

     ! Final output
     do n = 1, nglobal
        write(unit=io_num, fmt="(a20,i5,3i3,i5)" ) &
              s_name(n), current_date
        write(unit=io_num, fmt="(5es11.3)" ) g_out(:,n)
     end do ! n

  end if ! me

  end subroutine siteswrt_out


!==================================================================== >
end module Sites_ml
