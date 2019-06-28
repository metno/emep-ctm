! <PointSource_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
!> MODULE PointSource_mod
!!============================================================================
!! Implementation of NILU and other plume-rise methods, with Peter's
!! spread-plume
!!
!! The code will read coordinates and emissions from a file PointSources.txt.
!! Plume rise and hence effective emissions height  is calculated from the
!! stack parameters and meteorology, and the spread into the z-layers using
!! an assumed standard deviation of the emissions (=25% of the plume-rise
!! delta). Can work with any z-layers (no hard coding :-)
!!
!! Matthias Karl and Dave, Dec. 2012
!!============================================================================
!
!ToDo
!
!   Currently the meteorology assumes lowest layer has 90m depth. Should in
!   future also make wind-speed and dt/dz depend on z values.
!
!   Currently emissions make use of EmisNat arrays, since this was easier
!   for non SNAP emissions. Will modify.

module PointSource_mod

! --- prelim module to introduce point sources.
!
use CheckStop_mod,       only: CheckStop
use ChemSpecs_mod,       only: species
use Config_module,  only: KMAX_MID, PT, Pref, MasterProc, USES
use Debug_module,        only: DEBUG => DEBUG_EMISSTACKS
use Functions_mod,       only: StandardAtmos_kPa_2_km, Tpot_2_T
use GridValues_mod,      only: sigma_bnd, debug_proc &
                         , coord_in_processor, GridArea_m2
use Io_Progs_mod,        only: open_file, ios, read_line, datewrite
use LocalVariables_mod,  only: Grid ! Grid-meteorology
use MetFields_mod,       only: pzpbl, z_bnd, z_mid
use Par_mod,             only: me, LIMAX, LJMAX, limax, ljmax, &
                              IRUNBEG,JRUNBEG, & ! TMP for debug
                              gi0, gi1, gj0, gj1
use PlumeRise_mod,       only: Plume_PreggerFriedrich, Plume_ASME, Plume_NILU
use ZchemData_mod,  only: rcemis, temp, pp
use SmallUtils_mod,      only: find_index, wordsplit
use TimeDate_mod,        only:  nydays
use PhysicalConstants_mod,only: AVOG,PI
implicit none
private

public  :: readstacks       ! Read LPS data
public  :: get_pointsources ! get emissions
private :: spread_plume     ! From Peter
private :: gauss_integral   ! "     "
private :: errorfunctionc   ! "     "

integer, parameter, private :: NzMax=20   ! number of layers, tmp
integer, parameter, private :: MAX_STACKS=9999  ! number of lines in input files
integer, parameter, private :: MAX_POLLS =  20  ! number of emitted pollutants
integer, parameter, private :: NeMax=10    ! number of emitted tracers
real, private, dimension(NzMax-1), save  :: layer_z ! layer boundaries

type, public :: stacktype
  character(len=20) :: name   ! used as label
  integer :: cc      ! country-code (could be an LPS code too in future)
  integer :: snap
  real    :: lat, long
  integer :: i,j     ! coords in EMEP domain
  real    :: hs      ! Physical ht., m
  real    :: d       ! diameter, m
  real    :: Vs      ! exit velocity
  real    :: Ts      !  stack temp (K)
  real    :: flow    ! emission
  real    :: M       ! buoyancy param
  real    :: dh4     ! Calculated plume rise at 4 m/s
  real    :: flue    ! flue gas rate (m^3/s)
  real, dimension(MAX_POLLS) :: emis
  real    :: stdev   ! std. dev in plume 
end type stacktype
type(stacktype), private, dimension(0:MAX_STACKS) :: stack  ! for input data

logical, public, allocatable,dimension(:,:) :: pointsources
integer, private, save :: nstacks = 0
integer, private, save, dimension(MAX_POLLS) :: ispec_emis
integer, public, save :: nemis_found

contains
!--------------------------------------------------------------------------
subroutine readstacks(io)
  integer, intent(in) :: io

  character(len=20) :: txtcode   ! Just set as CCSNAP for now
  integer :: cc, snap, i, j, iloc, jloc
  real    :: he            ! Eff. ht.
  integer :: n, ns, k, kk, iocode
  real :: lat, long        ! coords
  real :: hs,Ts, diam, Vs  ! physical ht, Tempm, Diameter, flow-speed
  real :: dh4          !  Plume rise at 4m/s from Preggers+Friedrich
  character(len=100) :: fname
  character(len=200) :: txtinput  ! Big enough to contain one input record
  real   , dimension(MAX_POLLS) :: emis
  integer, save, dimension(MAX_POLLS) :: icol_emis  ! column numbers of used emissions
  
  real :: he_max  ! max limit for emissions, found from max(he+2*stdev)
  real :: p
  logical :: first_header=.true., debug_flag=.false.
  character(len=20), dimension(MAX_POLLS+9) :: words
  integer :: errcode, ncol, emcol, nemis_cols, ispec, nn, nwords, iemis
  debug_flag = DEBUG !! .and. MasterProc 

  ! The plume rise methodology is very approximate, so we calculate
  ! layer thickness once, for a standard atmosphere
  ! We also set k=1 for the level nearest the ground
  do kk = 2, KMAX_MID
    k = KMAX_MID - kk  +1
    p = PT + sigma_bnd(kk)*(Pref -PT)
    layer_z(k) = 1.0e3 * StandardAtmos_kPa_2_km( 0.001*p) 
    if( MasterProc) write(*,"(a,i3,2f12.3)") "StackLayer", kk, layer_z(k)
  end do

  he_max = 0.0
  ns  = 0 ! number of stack-data entry, this processor

  allocate( pointsources(LIMAX, LJMAX))
  pointsources(:,:) = .false.

  if(MasterProc)then
    fname = "PointSources.txt" !TEST
    call open_file(io,"r",fname,needed=.true.)
    call CheckStop(ios,"open_file error on " // fname )
  end if

  ! First, we read file to get number and list of country codes
  n = 0
  do 
    call read_line(io,txtinput,iocode)
    if(iocode< 0) exit               ! End of file


    ! 1) Header line
    if(txtinput(1:1)=="#")then    ! Figure out emissions indices from header
       if(.not.first_header)cycle ! Currently units line
       first_header = .false.    

       call wordsplit(txtinput,MAX_POLLS+9,words,nwords,errcode,&
         separator=",")

       nemis_found =0           ! Total number of used emissions
       nemis_cols  = nwords - 9 ! number of emission cols in input data

      do ncol = 10,  nwords
        ispec  = find_index( words(ncol), species(:)%name ) 

        if( ispec > 0 ) then
          nemis_found = nemis_found + 1
          emcol=ncol-9  
          ispec_emis(nemis_found) = ispec
          icol_emis(nemis_found)  = emcol
          if(MasterProc) print *, "STACKem ", nemis_found, iemis, species(ispec)%name
        end if
      end do
      cycle
    else

    ! 2) Data lines
      read(unit=txtinput, fmt=*, iostat=iocode ) txtcode, cc, snap, hs, &
         lat, long, Ts, diam, Vs, ( emis(nn), nn=1,nemis_cols )
      if(MasterProc) print *, "STACK n HERE  ", n, ns, trim(txtcode),  lat, long

    end if

    n  = n + 1     ! number, this emis file
    ns = ns + 1    ! number, all emis files
    if(MasterProc) &
      call CheckStop(n>MAX_STACKS,"Too many stacks? "//trim(fname))

    ! check if coords are on this processor and get local i/j index
    if(.not.coord_in_processor(long,lat,iloc=iloc,jloc=jloc)) cycle
    ! on the correct processor
    if(debug_flag)&
      print "(a,99i6)",'Stack Local coords are',ns,iloc,jloc

    stack(ns)%name   = trim(txtcode)
    stack(ns)%cc     = cc
    stack(ns)%hs     = hs 
    stack(ns)%d      = diam
    stack(ns)%Vs     = Vs
    stack(ns)%Ts     = Ts  
    stack(ns)%snap   = snap

    stack(ns)%lat    = lat 
    stack(ns)%long   = long

    stack(ns)%i      = iloc
    stack(ns)%j      = jloc
    pointsources(iloc,jloc)  = .true.

    ! Use air-temp = 283.15 for now, see Pregger+Friedrich
    ! Bouyancy flux
    stack(ns)%flow   =  0.25 * 9.81 *  diam**2 * Vs * (Ts-283.15)/Ts
      
    ! Flue gas flow (Pregger&Friedrich)
    stack(ns)%flue   =  0.25 * PI * Vs * diam**2

    do nn = 1, nemis_found
      ncol = icol_emis(nn)
      stack(ns)%emis(nn) = emis(ncol)
      if( MasterProc ) print *, "STACK SETS ", ns,  nn, ncol, emis(ncol)
    end do

    if(debug_flag ) write(*,"(a,i5,es10.3)") " StacksEmis ", &
      ispec_emis(1), ispec_emis(2)    !DS iemist, emisa1
    write(*,"(a,3i5,4x,2i4,4x,4i4)") " StacksCrds ", &
        me, i, j,iloc,jloc,  limax, ljmax
    if(debug_flag ) write(*,"(a,3i5,4x,2i4, 2f8.3,10f10.3)") " StacksProc ", &
       me, i, j,iloc,jloc, lat, long,  hs, dh4, he, stack(ns)%flow

      !if ( he > he_max) he_max = he
  end do 
  nstacks = ns
  if(debug_flag) write(*,*) "Stacks Read. Done"

  if(MasterProc) close(io)

! Find the maximum height of emissions across all  pollutants
! given as kup = k upwards from the ground
!   if ( layer_z(k) > he_max  )  then
!      Stacks%kup    = k
!      exit
!   end if
! end do
end subroutine readstacks
!----------------------------------------------------------------------------
subroutine get_pointsources(i,j, debug_flag)
  integer, intent(in) :: i,j
  logical, intent(in), optional :: debug_flag
  integer :: n, Nz=10,k,kk
  integer :: ispec, nn
  real :: he, dh, stddev, Ta,  Ta_C, dtdz ! FAKE
  real, dimension(NzMax) :: fraction_per_layer
  real :: emiss
  real    :: Mh,ObukhovL
  real    :: uconv, uconv1
  logical :: myfirstcall = .true.

  logical :: mydebug

  mydebug = .false.
  if(present(debug_flag)) mydebug = debug_flag

  Ta = Grid%theta_ref * Tpot_2_T( pp(KMAX_MID) )
  Ta_C = Ta - 273.15

  do n = 1, nstacks
    if(stack(n)%i/=i .or. stack(n)%j/=j) cycle

    select case(USES%PlumeMethod)
    case("ASME")
      dtdz =  ( temp(KMAX_MID) - Grid%t2 ) / layer_z(1)  ! THETA or T?

      dh =  Plume_ASME( stack(n)%hs, stack(n)%flow, Grid%u_ref, Ta, dtdz )
      he = stack(n)%hs + dh
      if(mydebug .and. myfirstcall ) then
        write(*,"(a,i3,f10.4,2f10.4,a,2f10.3,es10.2)") &
        "Plume_ASME: stack he dh",n,he,dh, dtdz, &
        "Tpot?",  Grid%theta_ref, Ta, pp(KMAX_MID)
      end if

    case("NILU")
      ObukhovL = 1.0e4  ! safety
      if (abs(Grid%invL) > 1.0e-4)  ObukhovL = 1.0/Grid%invL

      he = Plume_NILU( stack(n)%hs, Grid%ustar, stack(n)%Vs, stack(n)%d, &
          stack(n)%Ts, Ta, z_mid(i,j,KMAX_MID), pzpbl(i,j), Grid%u_ref,   &
          ObukhovL ) 

      ! For now we don't allow downwash. Too uncertain
      he = max( he, Stack(n)%hs)
      if(mydebug .and. myfirstcall ) then
        if (n==1) write(*,"(a,i3,f10.4,f10.4)") "Plume_NILU: stack he dh",  &
               n,he,he - Stack(n)%hs
       !if (n==1) write(*,*) i,j,GridArea_m2(i,j)
      end if        

    case("PVDI")
      ! Mh: Emitted heat flux in MW, VDI(1985)       
      Mh = 1.36e-3 * stack(n)%flue * (stack(n)%Ts-283.15)
      dh = Plume_PreggerFriedrich( Mh, Grid%u_ref, stack(n)%Vs, stack(n)%d )
      he = stack(n)%hs + dh      
      if(DEBUG .and. n==1) &
        write(*,"(a,i3,f10.4,f10.4)") "Plume_PVDI: stack he dh",n,he,dh

    case default
      call CheckStop("Unknown PlumeMethod: "//trim(USES%PlumeMethod))
    end select

    dh     = he - Stack(n)%hs

    stddev = 0.25 * dh   ! Bieser et al. assumed plume top and bottom at
                         ! 0.5 dh, so we assume this was 2 std-dev
    stddev = max( stddev, 10.0 )     ! Just to force some mixing
    
    fraction_per_layer = 0.0
    call spread_plume( he, stddev, layer_z(1:Nz-1), fraction_per_layer(1:Nz),Nz)

    if(mydebug)&
     call datewrite("STACK-"//trim(USES%PlumeMethod) //"-"//trim(stack(n)%name)&
            // " DD/MM hh:00", &  ! Format for date
      (/ he, Ta-273.15, Grid%u_ref,  Grid%invL, ObukhovL, &
           (fraction_per_layer(k), k=1,6) /), txt_pattern=.true.)

    ! CONVERT EMISSION FROM kg/year to (molec/cm2/s)/dz(cm)      
    uconv=0.1 /(nydays*24.0*3600.0)                ! Kg/yr --> 10^4 g/s
    uconv=uconv/GridArea_m2(i,j)                   ! --> g/s/cm2

    do nn = 1, nemis_found
      ispec = ispec_emis(nn)
      uconv1=uconv*AVOG/species(ispec)%molwt          ! --> molecules/s/cm2
      emiss = stack(n)%emis(nn) *uconv1

      do kk=1,Nz
        k = KMAX_MID - kk + 1
        rcemis(ispec,k)  = rcemis(ispec,k) + fraction_per_layer(kk) *emiss &
              / ( 100.0*(z_bnd(i,j,k) - z_bnd(i,j,k+1)) )  ! width in cm
       if(mydebug) write(*,"(a20,i3,f7.3,3es10.3)") &
         " StacksEmisRC "//trim(species(ispec)%name),   &
         k, fraction_per_layer(kk), emiss, rcemis(ispec,k)
      end do !kk
    end do ! iemis 

    if(mydebug ) write(*,"(a,es10.3,es10.3,a,2es10.3)") " StacksEmisConv ",   &
      uconv,GridArea_m2(i,j),trim(species(ispec)%name),species(ispec)%molwt,emiss

   end do
   myfirstcall = .false.
end subroutine get_pointsources
!----------------------------------------------------------------------------
subroutine spread_plume( h, w, layer_z, fraction_per_layer,Nz)
!gives the fraction of a gaussian distribution between layer_z levels
!gaussian has sigma=w and is centered at h
!fraction_per_layer(1) is integral until layer_z(1)
!fraction_per_layer(Nz) is integral from layer_z(Nz-1) until infinity
  implicit none
  integer, intent(in) ::Nz !number of layers
  real, intent(in) :: h, w !plume ht, width
  real, intent(in), dimension(Nz-1) :: layer_z !(z boundaries)
  real, intent(out), dimension(Nz) :: fraction_per_layer
  real :: sum,large
  integer ::i
  
  large=1.0E10

  fraction_per_layer(1)=gauss_integral(-large,layer_z(1),h,w)

  do i=2,Nz-1
    fraction_per_layer(i)=gauss_integral(layer_z(i-1),layer_z(i),h,w)
  end do
  fraction_per_layer(Nz)=gauss_integral(layer_z(Nz-1),large,h,w)

  !normalize (only to get last digits right)
  sum=0.0
  do i=1,Nz
    sum=sum+fraction_per_layer(i)
  end do

  do i=1,Nz
    fraction_per_layer(i)=fraction_per_layer(i)/sum

    if( fraction_per_layer(i) > 1.0 ) then
      print *, "FRAC WRONG ", i, h,w,Nz, fraction_per_layer(i)
      call CheckStop("FRAC WRONG in spread_plume")
    end if
  end do
end subroutine spread_plume

!--------------------------------------------------------------------------
real function gauss_integral(x1,x2,x0,sigma) result(integral)
  !integrates the area between x1 and x2 of a normalized gauss function 
  !         1/(sqrt(2PI)*sigma)*exp(-(x-x0)**2/(2sigma**2))
  
  real, intent(in)::x1,x2,x0,sigma
  real :: integ1,integ2,isigma_sqrt2

  isigma_sqrt2=1.0/(sqrt(2.0)*sigma) !scale factor for sigma
  
  call CheckStop( x2< x1, "WARNING: negative integral")
  !if(x2<x1)write(*,*)'WARNING: negative integral', x2, x1
  integ1=0.5*errorfunctionc((x1-x0)*isigma_sqrt2)
  if(x1-x0>0)integ1=1.0-integ1
  
  integ2=0.5*errorfunctionc((x2-x0)*isigma_sqrt2)
  if(x2-x0>0)integ2=1.0-integ2
  
  integral=integ2-integ1

end function gauss_integral

!--------------------------------------------------------------------------
real function errorfunctionc(x_in) result(erfc)
  !complementary error function
  !only defined for x>=0. For x<0 gives erfc(abs(x))
  !erf(0)=1 erfc(infinity)=0.0
  !gives approx 7 correct digits.

  implicit none
  real ,intent(in)::x_in
  integer ::i
  real ::d2x2,y,fac,x

  fac=1.0
  x=abs(x_in)

  if(abs(x)<3.0)then
    !formula for small x
    erfc=x  
    fac=-3*2
    y=1.0
    do i=1,40,1
      y=-y*x*x/i
      erfc=erfc+y*x/(2*i+1)
    end do
    erfc=1.0-erfc*2/sqrt(pi)
!   write(*,*)'result small x',x_in,erfc
  else
    !formula for large x
    d2x2=1.0/(2*x*x)
    erfc=1.0-d2x2
    fac=-1.0
    y=d2x2
    do i=3,30,2
      fac=-fac*i
      y=y*d2x2
      erfc=erfc+fac*y
    end do
    erfc=exp(-x*x)*erfc/(x*sqrt(pi))
!   write(*,*)'result large x',x_in,erfc
  end if
end function  errorfunctionc

!--------------------------------------------------------------------------
  

!program ReadPlume
!  use EmisDef_mod,  only :NEMIS_FILE, EMIS_FILE
!  use StackRead_mod
!  implicit none
!  integer, parameter ::Nz=20 !number of layers
!  real, dimension(Nz-1) :: layer_z !(z boundaries)
!  integer :: i
!  integer :: iemis, snap, cc
!
!  integer, parameter :: IONUM = 10
!  do i=1,Nz-1
!     layer_z(i)=50+i*50!example
!  end do
!
!  call readstacks(IONUM, Nz, layer_z)
!
!  do iemis = 1, NEMIS_FILE
!    cc = 10
!  do i = 1, Nz-1
!    write(*,"(2a8,i4,f9.2,6f8.3)") "Final ", trim(EMIS_FILE(iemis)),  i, &
!         layer_z(i), (emis_vfac(cc,iemis,snap,i), snap=1,6)
!    write(*,"(2a8,i4,f9.2,6f8.3)") "Avg:  ", trim(EMIS_FILE(iemis)),  i, &
!         layer_z(i), (emis_vfac(99,iemis,snap,i)/(1.0e-10+sum(emis_vfac(99,iemis,snap,:))), snap=1,6)
!  end do  !i z
!  end do  !iemis 
!end program ReadPlume

end module PointSource_mod
