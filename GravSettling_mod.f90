! <GravSettling_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.45>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2022 met.no
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
module GravSettling_mod

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!EXPERIMENTAL!! Not in production use !!!
! Testing for gravitational settling
! Code modified from original ash module written by Birthe Steensen Still
! EXPERIMENTAL!!! In particular, the numerics are probably too simple
! and diffuse too much. Tests for coarse sea-salt show only a limited
! impact, but probably more attention is needed to get larger particle
! treated properly..
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use CheckStop_mod,         only: CheckStop
use Chemfields_mod,        only: xn_adv
use ChemDims_mod,          only: NSPEC_SHL
use ChemGroups_mod,        only: chemgroups
use ChemSpecs_mod,         only: species_adv
use DerivedFields_mod,     only: f_3d,d_3d ! debug output
use GridValues_mod,        only: A_mid,B_mid,A_bnd,B_bnd
use MetFields_mod,         only: roa,th,ps
use Config_module,         only: KMAX_MID,KMAX_BND,dt_advec,MasterProc,&
                                IOU_INST,num_lev3d,lev3d
use Par_mod,               only: MAXLIMAX,MAXLJMAX,li0,li1,lj0,lj1
use PhysicalConstants_mod, only: GRAV
use SmallUtils_mod,        only: find_index

use GasParticleCoeffs_mod, only: nddep, DDspec, DDmapping !DS
use GridValues_mod,        only:  debug_proc, debug_li, debug_lj !DS
use ChemSpecs_mod, only: species_adv

implicit none
private

public :: gravSettling

real, parameter :: &
  slinnfac = 1.0  ,&
  density = 2.5E03,& !3.0E03
  F = 0.8, wil_hua = F**(-0.828) + 2*SQRT(1.07-F)

contains

subroutine gravSettling()
  real                        :: ztempx,vt,zsedl,tempc,knut
  real                        :: Re, vt_old, ztemp
  real,dimension(KMAX_MID)    :: zvis,p_mid,zlair,zflux,zdp1,num_sed
  real,dimension(KMAX_BND)    :: p_full
  integer                     :: i,j,k,ash,n,b
  integer,save                :: bins=0
  logical                     :: first_call = .true.
  integer :: icmp, idef, ispec, nadv !DS
  real :: DpgV, tmpVal !DS
  logical :: dbgIJ

! need to calculate
! zvis -> dynamic viscosity of air.. dependent on temperature

!DS  type :: sediment
!DS    integer :: spec,n3d_gravset
!DS    real    :: diameter
!DS  end type sediment
!DS  type(sediment),save,allocatable,dimension(:) :: grav_sed


  if (first_call) then
    if(MasterProc) &
      write(*,*) "GravSettling called!", nddep, trim(DDspec(1)%name)

!DS    ash=find_index("ASH",chemgroups(:)%name)
!DS    bins=0
!DS    if(ash>0)&
!DS      bins=size(chemgroups(ash)%specs)
!DS    select case(bins)
!DS    case(7)
!DS      allocate(grav_sed(bins))
!DS      grav_sed(:)%spec = chemgroups(ash)%specs(:)-NSPEC_SHL
!DS      grav_sed(:)%diameter = [0.1,0.3,1.0,3.0,10.0,30.0,100.0]*1e-6
!DS    case(9)
!DS      allocate(grav_sed(bins))
!DS      grav_sed(:)%spec = chemgroups(ash)%specs(:)-NSPEC_SHL
!DS      grav_sed(:)%diameter = [4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,25.0]*1e-6
!DS!   case(10)
!DS!     allocate(grav_sed(bins))
!DS!     grav_sed(:)%spec = chemgroups(ash)%specs(:)-NSPEC_SHL
!DS!     grav_sed(:)%diameter = [2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,25.0]*1e-6
!DS    case default
!DS      if(MasterProc) &
!DS        write(*,"(A,I0,A)") "Unsupported number of ASH bins ",bins,", skip gravset."
!DS    end select
!DS
!DS    if(allocated(grav_sed))then
!DS      grav_sed(:)%n3d_gravset=0
!DS      do n=1,size(f_3d)
!DS        if(f_3d(n)%class/='USET')cycle
!DS        b=find_index(f_3d(n)%txt,species_adv(grav_sed(:)%spec)%name)
!DS        if(b<1)cycle
!DS        select case(f_3d(n)%subclass)
!DS        case('gravset_3D')
!DS          grav_sed(b)%n3d_gravset=n
!DS          f_3d(n)%unit="m/s"
!DS          f_3d(n)%scale=1.0
!DS        end select
!DS      end do
!DS    end if

!DS    first_call = .false.
  end if !first_call
!DS  if(.not.allocated(grav_sed))&
!DS    return

  do j = lj0,lj1
    do i = li0,li1
      dbgIJ = ( debug_proc .and. i==debug_li .and. j==debug_lj )
      do k = 1,KMAX_MID
        ! dynamic viscosity of air after Prup.Klett in [Pa s]
        tempc = th(i,j,k,1) - 273.15
        if (tempc >= 0.0 ) then
          zvis(k) = (1.718 + 0.0049*tempc)*1.E-5
        else
          zvis(k) = (1.718 + 0.0049*tempc - 1.2E-05*(tempc**2))*1.E-5
        end if

        ! mean free path of air (Prupp. Klett) in [10^-6 m]
        p_mid(k) = A_mid(k)+B_mid(k)*ps(i,j,1)
        zlair(k) = 0.066 *(1.01325E+5/p_mid(k))*(th(i,j,k,1)/293.15)*1.E-06

        ! air mass auxiliary  variable --> zdp1 [kg/(m^2 *s)]
        p_full(k) = A_bnd(k)+B_bnd(k)*ps(i,j,1)
      end do
      p_full(KMAX_BND) = A_bnd(KMAX_BND)+B_bnd(KMAX_BND)*ps(i,j,1)

      do k = 1,KMAX_MID
        zdp1(k)=(p_full(k+1) - p_full(k))/(GRAV*dt_advec) ! do outside of k-loop????
      end do

      !DS do b = 1,bins
      do icmp = 1, nddep  ! size(DDmapping(:)%name )

        !f ( dbgIJ) write(*,*)'dbgVS',icmp, trim(DDspec(icmp)%name)
        if ( DDspec(icmp)%is_gas ) CYCLE

        DpgV=DDspec(icmp)%DpgV
        !if ( dbgIJ) write(*,*)'dbgVSV',icmp, trim(DDspec(icmp)%name), DpgV
                
        if ( DpgV < 1.0e-6 ) CYCLE ! too small 

        !idef = DDmapping(icmp)%idef
        SPECLOOP: do ispec = 1, size(DDmapping(icmp)%advspecs)  ! Real species now

          nadv = DDmapping(icmp)%advspecs(ispec)  ! Real species now

          if ( maxval( xn_adv(nadv,i,j,:) )<1.0e-20) CYCLE
          if ( dbgIJ.and.first_call) write(*,*)'dbgVSX',icmp,ispec, trim(DDspec(icmp)%name), trim(species_adv(nadv)%name)
          if ( dbgIJ .and. trim(species_adv(nadv)%name) == 'NO3_c') tmpVal = xn_adv(nadv,i,j,19)


        do k = 1,KMAX_MID-1
          knut = 2*zlair(k)/DpgV !DS grav_sed(b)%diameter
          !DS ztemp = 2.*((grav_sed(b)%diameter/2)**2)*(density-roa(i,j,k,1))*GRAV/ &! roa [kg m-3]
          ztemp = 2.*((DpgV/2)**2)*(density-roa(i,j,k,1))*GRAV/ &! roa [kg m-3]
                  (9.*zvis(k))![m/s]
          ! with Cunningham slip-flow correction
          vt = ztemp*slinnfac*     &
              (1.+ 1.257*knut+0.4*knut*EXP(-1.1/(knut))) ![m/s]
          !DS Re = grav_sed(b)%diameter*vt/(zvis(k)/roa(i,j,k,1))
          Re = DpgV*vt/(zvis(k)/roa(i,j,k,1))
          vt_old = vt
          vt = vt/wil_hua
          num_sed(k)= vt

          ! calculation of sedimentation flux zflux[kg/(m^2 s)]=zsedl*zdp1
          ! definition of  zflux=vt*ztm1(:,:,jt)*zdens
          ! compute flux in terms of mixing ratio zsedl= zflux/zdp1 -->>zsedl [kg/kg]
          ! change of tracer tendency according to loss of tracer
          ! due to sedimentation from the box
          ! unit of zdp1 kg of air m-2 s-1

          ztempx = min(1.0,vt*roa(i,j,k,1)/zdp1(k)) ! 1, loss is limited to content of box
          !DS zsedl = ztempx*xn_adv(grav_sed(b)%spec,i,j,k)        ! kg kg-1, loss in terms of mixing ratio ! blir likt uavhengig av vt
          zsedl = ztempx*xn_adv(nadv,i,j,k)        ! kg kg-1, loss in terms of mixing ratio ! blir likt uavhengig av vt
          zflux(k) = zsedl*zdp1(k)                     ! --> [kg m-2 s-1]

          ! loss of mass in layer
          !DS xn_adv(grav_sed(b)%spec,i,j,k) = xn_adv(grav_sed(b)%spec,i,j,k) - zsedl !  kg kg-1
          xn_adv(nadv,i,j,k) = xn_adv(nadv,i,j,k) - zsedl !  kg kg-1
        end do

        ! teste å gjøre det i en ny k-loop så det ikke blir så effektivt

        do  k = 1,KMAX_MID-1
          ! "arrival" of sedimented mass in box below
          if(k<KMAX_MID)then
            !DS xn_adv(grav_sed(b)%spec,i,j,k+1) =  xn_adv(grav_sed(b)%spec,i,j,k+1) +  zflux(k)/zdp1(k+1)
            xn_adv(nadv,i,j,k+1) =  xn_adv(nadv,i,j,k+1) +  zflux(k)/zdp1(k+1)
          end if
        end do
          if ( dbgIJ .and. trim(species_adv(nadv)%name) == 'NO3_c' &
              .and. xn_adv(nadv,i,j,19)>1.e-20 ) write(*,'(a,2es12.3,f8.3)') &
               'dbgVSba:', tmpVal, xn_adv(nadv,i,j,19), 100*(tmpVal-xn_adv(nadv,i,j,19))/xn_adv(nadv,i,j,19)

        !DS n=grav_sed(b)%n3d_gravset
        !DS if(n>0)&
        !DS   d_3d(n,i,j,:,IOU_INST)=num_sed(lev3d(:num_lev3d))

        end do SPECLOOP
      end do
    end do
  end do
   first_call = .false.

! Multilayer crossing is here no ralised!!!
! sedimentation velocity is in effect limited to z/delt
! sedimentation to the ground from first layer sflx --> [kg m-2 s-1]
! sflux = zflux(:,1) ! sedimenterer ikke på bakken, bare sender det til nederste laget
end subroutine gravSettling


end module GravSettling_mod
