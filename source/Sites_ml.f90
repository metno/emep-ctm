! <Sites_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
use GridValues_ml,     only : lb2ij, i_fdom, j_fdom &
                              , i_local, j_local, A_mid, B_mid
use Io_ml,             only : check_file,open_file,ios &
                              , fexist, IO_SITES, IO_SONDES &
                              , Read_Headers,read_line
!CMR use ChemSpecs_adv_ml
!CMR use ChemSpecs_shl_ml,  only : NSPEC_SHL
!CMR use ChemChemicals_ml,  only : species               ! for species names
use ChemSpecs
use ChemGroups_ml,     only : OXN_GROUP, PMFINE_GROUP, PMCO_GROUP
use Met_ml,            only : meteo
use MetFields_ml,      only : t2_nwp, th, pzpbl  &  ! output with concentrations
                              , z_bnd, z_mid, roa, Kz_m2s, q
use MetFields_ml,      only : u_xmj, v_xmi, ps
use ModelConstants_ml, only : NMET,PPBINV,PPTINV, KMAX_MID, MasterProc &
                              ,KMAX_BND,PT,ATWAIR, NPROC, DEBUG => DEBUG_SITES &
                              ,DomainName, RUNDOMAIN, IOU_INST, SOURCE_RECEPTOR
use NetCDF_ml,         only : Create_CDF_sondes,Out_CDF_sondes
use Par_ml,            only : li0,lj0,li1,lj1 &
                              ,GIMAX,GJMAX,IRUNBEG,JRUNBEG&
                              ,GI0,GI1,GJ0,GJ1,me,MAXLIMAX,MAXLJMAX
use SmallUtils_ml,     only : find_index
use Tabulations_ml,    only : tab_esat_Pa
use TimeDate_ml,       only : current_date
use TimeDate_ExtraUtil_ml,   only : date2string
use KeyValueTypes,       only : KeyVal, KeyValue, LENKEYVAL

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
integer, public, save :: nglobal_sites, nlocal_sites
integer, private, save :: nglobal_sondes, nlocal_sondes

! site_gindex stores the global index n asociated
! with each processor and local site

integer, private, save, allocatable,dimension (:,:)  :: site_gindex
integer, private, save, allocatable,dimension (:,:) :: sonde_gindex

integer, public, save, dimension (NSITES_MAX) :: &
        site_x, site_y, site_z      &! local coordinates
       , site_gn                        ! number in global
real, public, save, dimension (NSITES_MAX) :: Sites_lon= -999, Sites_lat= -999

integer, private, save, dimension (NSITES_MAX) :: &
         site_gx, site_gy, site_gz    ! global coordinates
integer, private, save, dimension (NSONDES_MAX) ::  &
         sonde_gx, sonde_gy   &        ! global coordinates
       , sonde_x, sonde_y     &        ! local coordinates
       , sonde_gn                       ! number in global
real, public, save, dimension (NSONDES_MAX) :: Sondes_lon= -999, Sondes_lat= -999, ps_sonde=0.0

! Values from My_Outputs_ml gives ... =>
integer, private, parameter :: & ! Total No., without counting levels
   NSPC_SITE  = NADV_SITE + NSHL_SITE + NXTRA_SITE_MISC + NXTRA_SITE_D2D &
  ,NSPC_SONDE = NADV_SONDE + NSHL_SONDE + NXTRA_SONDE
integer, public, parameter :: & ! Total No., levels included
   NOUT_SITE  = NSPC_SITE * 1 &
  ,NOUT_SONDE = NSPC_SONDE* NLEVELS_SONDE

character(len=50), public, save, dimension(NSITES_MAX) :: site_name
character(len=50), private, save, dimension(NSONDES_MAX):: sonde_name
character(len=20), private, save, dimension(NSPC_SITE)  :: site_species
character(len=20), private, save, dimension(NSPC_SONDE) :: sonde_species

character(len=70), private :: errmsg ! Message text
integer, private :: d                 ! processor index
integer, private :: i, n, nloc, ioerr ! general integers
integer, parameter, private :: Spec_Att_Size=20,N_Spec_Att_MAX=5,NSPECMAX=max(NSPC_SITE,NSPC_SONDE)
character(len=Spec_Att_Size)  :: Spec_AttributeNames(NSPECMAX,N_Spec_Att_MAX)
character(len=Spec_Att_Size)  :: Spec_AttributeValues(NSPECMAX,N_Spec_Att_MAX)
integer, private :: i_Att !Spec attribute index
integer :: NSpec_Att !number of Spec attributes defined


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
        site_x, site_y, site_z, site_gn, &
        site_name)

  call Init_sites("sondes",IO_SONDES,NSONDES_MAX, &
        nglobal_sondes,nlocal_sondes, &
        sonde_gindex, sonde_gx, sonde_gy, sonde_gz, &
        sonde_x, sonde_y, sonde_z, sonde_gn, &
        sonde_name)

  call set_species(SITE_ADV,SITE_SHL,SITE_XTRA,site_species)
  call set_species(SONDE_ADV,SONDE_SHL,SONDE_XTRA,sonde_species)

  if ( DEBUG ) then
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
      lon=-999.0
      lat=-999.0
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

      if(trim(fname)=="sites")then
         if(lon>-990)Sites_lon(n) = lon
         if(lat>-990)Sites_lat(n) = lat
      endif
      if(trim(fname)=="sondes")then
         if(lon>-990)Sondes_lon(n) = lon
         if(lat>-990)Sondes_lat(n) = lat
      endif

      s_name(n)  = s !!! remove comments// comment
      if (DEBUG) write(6,"(a,i3,i4,a)") "sitesdef s_name : ", me, n, trim(s_name(n))
    endif

  enddo SITELOOP

  nglobal = n

  ! NSITES/SONDES_MAX must be _greater_ than the number used, for safety

  call CheckStop(n >= NMAX, &
      "Error in Sites_ml/Init_sites: increase NGLOBAL_SITES_MAX!")

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

      if (DEBUG) &
        write(6,"(a,i3,a,2i3,3i4,a,3i4)") "sitesdef Site on me : ", me, &
         " Nos. ", n, nlocal, s_gx(n), s_gy(n) , s_gz(n), " =>  ", &
          s_x(nlocal), s_y(nlocal), s_z(nlocal)
        write(6,"(a,i3,a,2i3,4a)") "SPODsite : "// trim(fname), me, &
         " Nos. ", n, nlocal, " ", trim(s_name(n)), " => ", trim(s_name(nlocal))

     endif

  enddo ! nglobal

  ! inform me=0 of local array indices:
  if(DEBUG) write(6,*) "sitesdef ", fname, " before gc NLOCAL_SITES", &
                           me, nlocal

  if ( .not.MasterProc ) then
    call MPI_SEND(nlocal, 4*1, MPI_BYTE, 0, 333, MPI_COMM_WORLD, INFO)
    if(nlocal>0) call MPI_SEND(s_n, 4*nlocal, MPI_BYTE, 0, 334, &
                               MPI_COMM_WORLD, INFO)
  else
    if(DEBUG) write(6,*) "sitesdef for me =0 LOCAL_SITES", me, nlocal
    do n = 1, nlocal
      s_gindex(me,n) = s_n(n)
    enddo
    do d = 1, NPROC-1
      call MPI_RECV(nloc, 4*1, MPI_BYTE, d, 333, MPI_COMM_WORLD,STATUS, INFO)
      if(nloc>0) call MPI_RECV(s_n_recv, 4*nloc, MPI_BYTE, d, 334, &
                               MPI_COMM_WORLD,STATUS, INFO)
      if(DEBUG) write(6,*) "sitesdef: recv d ", fname, d,  &
                  " zzzz nloc : ", nloc, " zzzz me0 nlocal", nlocal
      do n = 1, nloc
        s_gindex(d,n) = s_n_recv(n)
        if(DEBUG) write(6,*) "sitesdef: for d =", fname, d, &
          " nloc = ", nloc, " n: ",  n,  " gives nglob ", s_gindex(d,n)
      enddo ! n
    enddo ! d
  endif ! MasterProc

  if ( DEBUG ) write(6,*) 'sitesdef on me', me, ' = ', nlocal

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

  if ( DEBUG ) then
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
  endif ! DEBUG

  ! assign local data to out

  i_Att=0
  NSpec_Att=1 !number of Spec attributes defined
  do i = 1, nlocal_sites
    ix = site_x(i)
    iy = site_y(i)
    iz = site_z(i)

    i_Att=0
    do ispec = 1, NADV_SITE
      if (iz == KMAX_MID ) then ! corrected to surface
        out(ispec,i) = xn_adv( SITE_ADV(ispec) ,ix,iy,KMAX_MID ) * &
                       cfac( SITE_ADV(ispec),ix,iy) * PPBINV
      else                      ! Mountain sites not corrected to surface
        out(ispec,i)  = xn_adv( SITE_ADV(ispec) ,ix,iy,iz ) * PPBINV
      endif
      i_Att=i_Att+1
      Spec_AttributeNames(i_Att,1)='units'
      Spec_AttributeValues(i_Att,1)='ppb'
 
    enddo


    do ispec = 1, NSHL_SITE
      out(NADV_SITE+ispec,i)  = xn_shl( SITE_SHL(ispec) ,ix,iy,iz )
      i_Att=i_Att+1
      Spec_AttributeNames(i_Att,1)='units'
      Spec_AttributeValues(i_Att,1)='ppb'
    end do

    ! XTRA parameters, usually the temmp or pressure
    if (.not. SOURCE_RECEPTOR) then
      nn = NADV_SITE + NSHL_SITE
      do ispec = 1, NXTRA_SITE_MISC
        nn=nn+1
        select case ( trim(SITE_XTRA_MISC(ispec)) )
          case ( "T2" )
            out(nn,i)   = t2_nwp(ix,iy,1) - 273.15
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='K'
          case ( "th" )
            out(nn,i)   = th(ix,iy,iz,1)
!         case ( "hmix" )
!           out(nn,i)   = pzpbl(ix,iy)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='K'
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
           if(MasterProc .and. my_first_call) write(*,*) &
                 "WARNING: SITES D2D NOT FOUND"//trim(d2code)
           !cycle
          i_Att=i_Att+1
           out(nn,i)   = -999.9
          Spec_AttributeNames(i_Att,1)='output'
          Spec_AttributeValues(i_Att,1)='undefined'!to improve?
        else
           !call CheckStop( d2index<1, "SITES D2D NOT FOUND"//trim(d2code) )
           out(nn,i)   = d_2d(d2index,ix,iy,IOU_INST)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='?'!to improve?
        end if

        !May25 call CheckStop( d2index<1, "SITES D2D NOT FOUND"//trim(d2code) )
        !May25 out(nn,i)   = d_2d(d2index,ix,iy,IOU_INST)
        if( DEBUG ) &
          write(6,"(a,3i3,a,i4,es10.3)") "DEBUG ", me, nn, i,&
            trim(d2code), d2index, out(nn,i)
        call CheckStop( abs(out(nn,i))>1.0e99, &
          "ABS(SITES OUT: '"//trim(SITE_XTRA_D2D(ispec))//"') TOO BIG" )
      enddo
    endif
  enddo

!


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
  ! IO_SONDES is set in io_ml to be 30
  ! -------------------------------------------------------------------

  real, dimension(:,:,:,:), intent(in) ::  xn_adv
  real, dimension(:,:,:,:), intent(in) ::  xn_shl

  ! Output variables - none

  ! Local variables
  integer :: n, i, k,  ix, iy, nn, ispec   ! Site and chem indices
  integer ::  KTOP_SONDE 
  integer, dimension(NLEVELS_SONDE)      :: itemp
  real, dimension(KMAX_MID)              :: pp, temp, qsat, rh, sum_PM, sum_NOy
  real, dimension(NOUT_SONDE,NSONDES_MAX):: out

  ! Consistency check
  KTOP_SONDE = KMAX_MID - NLEVELS_SONDE + 1
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
      out(nn+1:nn+NLEVELS_SONDE,i) = PPBINV *  &
          xn_adv( SONDE_ADV(ispec) , ix,iy,KMAX_MID:KTOP_SONDE:-1)
      nn = nn + NLEVELS_SONDE
      i_Att=i_Att+1
      Spec_AttributeNames(i_Att,1)='units'
      Spec_AttributeValues(i_Att,1)='ppb'
    enddo

    do ispec = 1, NSHL_SONDE    !/ xn_shl  in molecules/cm3
      out(nn+1:nn+NLEVELS_SONDE,i) = xn_shl( SONDE_SHL(ispec) , &
          ix,iy,KMAX_MID:KTOP_SONDE:-1)
      nn = nn + NLEVELS_SONDE
      i_Att=i_Att+1
      Spec_AttributeNames(i_Att,1)='units'
      Spec_AttributeValues(i_Att,1)='molecules/cm3'
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
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='ug/m3'

        case ( "PMco" ) !!  PM data converted to ug m-3
          sum_PM(:) = 0.
          do k = 1, KMAX_MID
            sum_PM(k) = dot_product(xn_adv(PMCO_GROUP-NSPEC_SHL,ix,iy,k),&
                                  to_ug_ADV(PMCO_GROUP-NSPEC_SHL)) &
                      * roa(ix,iy,k,1)
          enddo !k
          out(nn+1:nn+KMAX_MID,i) = sum_PM(KMAX_MID:1:-1)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='ug/m3'

        case ( "NOy" )
          sum_NOy(:) = 0.
          do k = 1, KMAX_MID
            sum_NOy(k) = sum(xn_adv(OXN_GROUP-NSPEC_SHL,ix,iy,k))
          enddo
          out(nn+1:nn+KMAX_MID,i) = PPBINV * sum_NOy(KMAX_MID:1:-1)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='mix_ratio'

        case ( "RH   " )
          do k = 1,KMAX_MID
            pp(k) = A_mid(k) + B_mid(k)*ps(ix,iy,1)
            temp(k) = th(ix,iy,k,1)* Tpot_2_T( pp(k) )
            itemp(k) = nint( temp(k) )
            qsat(k)  = 0.622 * tab_esat_Pa( itemp(k) ) / pp(k)
            rh(k) = min( q(ix,iy,k,1)/qsat(k) , 1.0)
          end do
          out(nn+1:nn+NLEVELS_SONDE,i) =  rh(KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='fraction'

        case ( "z_mid" )
          out(nn+1:nn+NLEVELS_SONDE,i) =  z_mid(ix,iy,KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='m'

        case ( "p_mid" )
          out(nn+1:nn+NLEVELS_SONDE,i) = A_mid(KMAX_MID:KTOP_SONDE:-1) + &
                                    B_mid(KMAX_MID:KTOP_SONDE:-1)*ps(ix,iy,1)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='Pa'

        case ( "Kz_m2s" )
          out(nn+1:nn+NLEVELS_SONDE,i) =  Kz_m2s(ix,iy,KMAX_MID:KTOP_SONDE:-1)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='?/m2/s'

        case ( "th" )
          out(nn+1:nn+NLEVELS_SONDE,i) =  th(ix,iy,KMAX_MID:KTOP_SONDE:-1,1)
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='K'

        case ( "U" )
          out(nn+1:nn+NLEVELS_SONDE,i) = 0.5 &
             *( u_xmj(ix,iy,KMAX_MID:KTOP_SONDE:-1,1) &
             +u_xmj(ix-1,iy,KMAX_MID:KTOP_SONDE:-1,1) )
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='m/s'

        case ( "V" )
          out(nn+1:nn+NLEVELS_SONDE,i) = 0.5 &
             *( v_xmi(ix,iy,KMAX_MID:KTOP_SONDE:-1,1) &
             +v_xmi(ix,iy-1,KMAX_MID:KTOP_SONDE:-1,1) )
          i_Att=i_Att+1
          Spec_AttributeNames(i_Att,1)='units'
          Spec_AttributeValues(i_Att,1)='m/s'

        case ( "D3D" )
            call StopAll("D3D Sites out not defined")

      end select

      nn = nn + NLEVELS_SONDE

    enddo ! ispec (NXTRA_SONDE)

    ps_sonde(i)=ps(ix,iy,1)!surface pressure always needed to define the vertical levels

  enddo ! i (nlocal_sondes)


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
  integer, save, dimension(NTYPES):: prev_month = (/ -99, -99 /) ! Initialise
  integer, save, dimension(NTYPES):: prev_year = (/ -99, -99 /) ! Initialise
  integer :: ii,nn

  integer, parameter :: NattributesMAX=100
  character(len=100),allocatable  :: SpecName(:)
  character(len=100),allocatable  :: AttributeNames(:,:)
  character(len=100),allocatable  :: AttributeValues(:,:),CoordNames(:,:)
  character(len=100)              :: fileName
  real,allocatable  :: CoordValues(:,:)
  integer,allocatable  :: NAttributes(:),NCoords(:)
  integer  :: Nlevels,ispec,NSPEC,NStations,Ndoublestations
  real ::Values(KMAX_MID)
  integer ::i_Att_MPI

  select case (fname)
  case ("sites" )
     type = 1
  case ("sondes" )
     type = 2
  case default
     write(6,*) "non-possible type in siteswrt_out for ", fname
     return
  end select

  write(suffix,fmt="(i4)") prev_year(type)
  fileName = fname // "_" // suffix // ".nc"!Name of the NetCDF file. Will overwrite any preexisting file

  !  if ( MasterProc .and. current_date%month /= prev_month(type)) then
  if (current_date%year /= prev_year(type) & 
                                !We do not consider midnight as a new year
       .and.(current_date%month/=1.or.current_date%hour/=0.or.current_date%seconds/=0)&
       ) then
     prev_year(type) = current_date%year
     write(suffix,fmt="(i4)") current_date%year
     fileName = fname // "_" // suffix // ".nc"!Name of the NetCDF file. Will overwrite any preexisting file

     if (MasterProc  ) then

        if ( prev_year(type) > 0 ) close(io_num)  ! Close last-year file
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
        write(io_num,'(9999a)')"site,date", (",", (trim(s_species(i)) ),i=1,size(s_species))

        !defintions of file for NetCDF output

        if(trim(fname)=="sondes")then
           NLevels = NLEVELS_SONDE !number of vertical levels (counting from surface)
           NSPEC=NSPC_SONDE!number of species defined for sondes
        else
           NLevels = 1
           NSPEC=NSPC_SITE!number of species defined for sites
        endif

        NStations = nglobal!number of sondes or sites defined


        allocate(SpecName(NSPEC),AttributeNames(0:NStations,NattributesMAX))
        allocate(AttributeValues(0:NStations,NattributesMAX),CoordNames(0:NStations,NattributesMAX))
        allocate(CoordValues(0:NStations,NattributesMAX))
        allocate(NAttributes(0:NStations),NCoords(0:NStations))

        n = 0!index for global attributes. Do not modify this line.
        NAttributes(n) = 2 !number of global string attributes
        call CheckStop(NAttributes(n)>NattributesMAX,'NattributesMAX too small')
        AttributeNames(n,1) = "File_Type"!name of the global attribute. For instance "File_Type"
        AttributeValues(n,1)= trim(fname)!string for that attribute. For instance "Sondes"

        AttributeNames(n,2) = 'meteo_source'
        AttributeValues(n,2)= trim(meteo)

        NCoords(n) = 4 !number of real-valued global attributes defined
        call CheckStop(NCoords(n)>NattributesMAX,'Coords: NattributesMAX too small')
        CoordNames(n,1) = "Number_of_hours_bewtween_outputs"!name of the global attribute. For instance "Number of hours bewtween outputs"
        CoordValues(n,1)= 1.0 !value for that attribute. For instance "1.0"
        CoordNames(n,2) = "Number_of_stations_defined"
        CoordValues(n,2)= NStations
        CoordNames(n,3) = "Model_domain_x_size"
        CoordValues(n,3)= GIMAX
        CoordNames(n,4) = "Model_domain_y_size"
        CoordValues(n,4)= GJMAX
        if(IRUNBEG>1)then
           NCoords(n) = NCoords(n)+1
           CoordNames(n,NCoords(n)) = "Model_domain_x_shift_origin"
           CoordValues(n,NCoords(n))= IRUNBEG
        endif
        if(JRUNBEG>1)then
           NCoords(n) = NCoords(n)+1
           CoordNames(n,NCoords(n)) = "Model_domain_y_shift_origin"
           CoordValues(n,NCoords(n))= JRUNBEG
        endif

        do ispec=1,NSPEC
           SpecName(ispec)=trim(site_species(ispec))!name of the  species
        enddo


        do n = 1, NStations

           !check to avoid two stations with same name, since this will give problems when writing in 
           !the netcdf output
           Ndoublestations=0
           do nn = n+1, NStations
              if(trim(s_name(n))==trim(s_name(nn)))then
                 if(Ndoublestations<1)write(*,*)'ERROR, ', trim(fname),' with same name: '
                 Ndoublestations=Ndoublestations+1
                 write(*,*) 'station ',trim(s_name(n)),', number ',n,' and',nn
              endif
           enddo
           call CheckStop(Ndoublestations>0,&
                "Two stations with identical name. Remove or rename one of them")

           NAttributes(n) = 1 !number of string attributes defined for the variable
           call CheckStop(NAttributes(n)>NattributesMAX,'NattributesMAX too small')
           !NB: attribute with index 1 MUST be the name of the station
           AttributeNames(n,1) = 'Name_of_station'!name of the attribute. For instance "Station_Type"
           AttributeValues(n,1)= trim(s_name(n))!string for that attribute. For instance "Urban"


           if(trim(fname)=="sites")then
              NCoords(n)=5 !number of real-valued attributes defined for the sondes variable
           else
              NCoords(n)=4 !number of real-valued attributes defined for the sites variable
           endif
           call CheckStop(NCoords(n)>NattributesMAX,'Coords: NattributesMAX too small')

           CoordNames(n,1) = 'model_x_coordinate'!name of the attribute. For instance "Station_longitude"
           CoordValues(n,1)= s_gx(n)!value for that attribute. For instance "20.6"
           CoordNames(n,2) = 'model_y_coordinate'
           CoordValues(n,2)= s_gy(n)
           if(trim(fname)=="sites")then
              CoordNames(n,3) = 'model_level'
              CoordValues(n,3)= s_gz(n)
              CoordNames(n,4) = 'site_longitude_coordinate'
              CoordValues(n,4)= Sites_lon(n)
              CoordNames(n,5) = 'site_latitude_coordinate'
              CoordValues(n,5)= Sites_lat(n)
           endif
           if(trim(fname)=="sondes")then
              CoordNames(n,3) = 'sonde_longitude_coordinate'
              CoordValues(n,3)= Sondes_lon(n)
              CoordNames(n,4) = 'sonde_latitude_coordinate'
              CoordValues(n,4)= Sondes_lat(n)
           endif
        enddo
        !take Spec_Attributes from any processor with at least one site/sonde
        if(i_Att>0.and.i_Att/=NSPEC)then
           write(*,*)'MISSING species attribute? ',i_Att,NSPEC
        endif
        do d = 1, NPROC-1
           call MPI_RECV(i_Att_MPI, 4*1, MPI_BYTE, d, 746, MPI_COMM_WORLD,STATUS, INFO)
           if(i_Att_MPI>0)then
              if(i_Att_MPI/=NSPEC)then
                 write(*,*)'MISSING species attribute? ',i_Att_MPI,NSPEC
              endif
              call MPI_RECV(Spec_AttributeNames, Spec_Att_Size*N_Spec_Att_MAX*NSPECMAX, &
                   MPI_BYTE, d, 747, MPI_COMM_WORLD,STATUS, INFO)
              call MPI_RECV(Spec_AttributeValues, Spec_Att_Size*N_Spec_Att_MAX*NSPECMAX, &
                   MPI_BYTE, d, 748, MPI_COMM_WORLD,STATUS, INFO)
           endif
        enddo

        call Create_CDF_sondes(fileName,SpecName,NSPEC,NStations,AttributeNames,AttributeValues,&
             NAttributes,CoordNames,CoordValues,NCoords,Spec_AttributeNames,Spec_AttributeValues,&
             NSpec_Att,NLevels,debug=.false.)

        deallocate(SpecName,AttributeNames)
        deallocate(AttributeValues,CoordNames)
        deallocate(CoordValues)
        deallocate(NAttributes,NCoords)
        write(*,*)'Created ',trim(fileName)
     else
        !not MasterProc
        i_Att_MPI=i_Att
        call MPI_SEND(i_Att_MPI, 4*1, MPI_BYTE, 0, 746, MPI_COMM_WORLD, INFO)
        if(i_Att>0)then
           call MPI_SEND(Spec_AttributeNames, Spec_Att_Size*N_Spec_Att_MAX*NSPECMAX, MPI_BYTE, 0, 747, MPI_COMM_WORLD, INFO)
           call MPI_SEND(Spec_AttributeValues, Spec_Att_Size*N_Spec_Att_MAX*NSPECMAX, MPI_BYTE, 0, 748, MPI_COMM_WORLD, INFO)
        endif
        prev_year(type) = current_date%year
     endif ! MasterProc 

  endif ! current_date%year /= prev_year(type)



  if ( .not.MasterProc ) then   ! send data to me=0 (MasterProc)

     call MPI_SEND(nlocal, 4*1, MPI_BYTE, 0, 346, MPI_COMM_WORLD, INFO)
     call MPI_SEND(out, 8*nout*nlocal, MPI_BYTE, 0, 347, MPI_COMM_WORLD, INFO)
     if(trim(fname)=="sondes")&
          call MPI_SEND(ps_sonde, 8*nlocal, MPI_BYTE, 0, 347, MPI_COMM_WORLD, INFO)

  else ! MasterProc

     ! first, assign me=0 local data to g_out
     if ( DEBUG ) print *, "ASSIGNS ME=0 NLOCAL_SITES", me, nlocal

     do n = 1, nlocal
        nglob = s_gindex(0,n)
        g_out(:,nglob) = out(:,n)
        if(trim(fname)=="sondes")g_ps(n) = ps_sonde(n)
     enddo ! n

     do d = 1, NPROC-1
        call MPI_RECV(nloc, 4*1, MPI_BYTE, d, 346, MPI_COMM_WORLD,STATUS, INFO)
        call MPI_RECV(out, 8*nout*nloc, MPI_BYTE, d, 347, MPI_COMM_WORLD, &
             STATUS, INFO)
        if(trim(fname)=="sondes")call MPI_RECV(ps_sonde, 8*nloc, MPI_BYTE, d, &
             347, MPI_COMM_WORLD, STATUS, INFO)
        do n = 1, nloc
           nglob = s_gindex(d,n)
           g_out(:,nglob) = out(:,n)
           if(trim(fname)=="sondes")g_ps(nglob) = ps_sonde(n)
        enddo ! n
     enddo ! d

     ! some computers print out e.g. "2.23-123" instead of "2.23e-123"
     ! when numbes get too small. Here we make a correction for this:
     where( abs(g_out)>0.0 .and. abs(g_out)<1.0e-99 ) g_out = 0.0

     ! Final output
     do n = 1, nglobal

        !! Massimo Vieno change the ouput style make the output csv
        !! Oct 2012 Formatting changed (DS,AMV) for gfortran compliance
        !! and compactness. 
        write (io_num,'(a,9999(:,",",es10.3))') & 
             trim(s_name(n)) // date2string(", DD/MM/YYYY hh:00",current_date),& 
             ( g_out(ii,n), ii =1, nout ) 
        ! (The ':' format control item will stop processing once the g_out
        !  is done, avoiding runtime warnings.)

     enddo

     if(trim(fname)=="sondes")then
        NLevels = NLEVELS_SONDE !number of vertical levels (counting from surface)
        NSPEC=NSPC_SONDE!number of species defined for sondes
     else
        NLevels=1
        NSPEC=NSPC_SITE!number of species defined for sites
     endif
     allocate(SpecName(NSPEC))

     do ispec=1,NSPEC
        SpecName(ispec)=trim(site_species(ispec))!name of the variable for one sites/sonde and species          
     enddo ! n

     call Out_CDF_sondes(fileName,SpecName,NSPEC,g_out,NLevels,g_ps,debug=.false.)
     deallocate(SpecName)

  endif ! MasterProc

end subroutine siteswrt_out
!==================================================================== >
end module Sites_ml
