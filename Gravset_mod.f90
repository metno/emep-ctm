! <Gravset_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Gravset_mod

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Testing for ash gravitational settling
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

implicit none
private

public :: gravset

real, parameter :: &
  slinnfac = 1.0  ,&
  density = 2.5E03,& !3.0E03
  F = 0.8, wil_hua = F**(-0.828) + 2*SQRT(1.07-F)

contains

subroutine gravset()
  real                        :: ztempx,vt,zsedl,tempc,knut
  real                        :: Re, vt_old, ztemp
  real,dimension(KMAX_MID)    :: zvis,p_mid,zlair,zflux,zdp1,num_sed
  real,dimension(KMAX_BND)    :: p_full
  integer                     :: i,j,k,ash,n,b
  integer,save                :: bins=0
  logical                     :: first_call = .true.

! need to calculate
! zvis -> dynamic viscosity of air.. dependent on temperature

  type :: sediment
    integer :: spec,n3d_gravset
    real    :: diameter
  end type sediment
  type(sediment),save,allocatable,dimension(:) :: grav_sed


  if (first_call) then
    if(MasterProc) &
      write(*,*) "Gravset called!"

    ash=find_index("ASH",chemgroups(:)%name)
    bins=0
    if(ash>0)&
      bins=size(chemgroups(ash)%specs)
    select case(bins)
    case(7)
      allocate(grav_sed(bins))
      grav_sed(:)%spec = chemgroups(ash)%specs(:)-NSPEC_SHL
      grav_sed(:)%diameter = [0.1,0.3,1.0,3.0,10.0,30.0,100.0]*1e-6
    case(9)
      allocate(grav_sed(bins))
      grav_sed(:)%spec = chemgroups(ash)%specs(:)-NSPEC_SHL
      grav_sed(:)%diameter = [4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,25.0]*1e-6
!   case(10)
!     allocate(grav_sed(bins))
!     grav_sed(:)%spec = chemgroups(ash)%specs(:)-NSPEC_SHL
!     grav_sed(:)%diameter = [2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,25.0]*1e-6
    case default
      if(MasterProc) &
        write(*,"(A,I0,A)") "Unsupported number of ASH bins ",bins,", skip gravset."
    end select

    if(allocated(grav_sed))then
      grav_sed(:)%n3d_gravset=0
      do n=1,size(f_3d)
        if(f_3d(n)%class/='USET')cycle
        b=find_index(f_3d(n)%txt,species_adv(grav_sed(:)%spec)%name)
        if(b<1)cycle
        select case(f_3d(n)%subclass)
        case('gravset_3D')
          grav_sed(b)%n3d_gravset=n
          f_3d(n)%unit="m/s"
          f_3d(n)%scale=1.0
        end select
      end do
    end if

    first_call = .false.
  end if !first_call
  if(.not.allocated(grav_sed))&
    return

  do j = lj0,lj1
    do i = li0,li1
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

      do b = 1,bins
        do k = 1,KMAX_MID-1
          knut = 2*zlair(k)/grav_sed(b)%diameter
          ztemp = 2.*((grav_sed(b)%diameter/2)**2)*(density-roa(i,j,k,1))*GRAV/ &! roa [kg m-3]
                  (9.*zvis(k))![m/s]
          ! with Cunningham slip-flow correction
          vt = ztemp*slinnfac*     &
              (1.+ 1.257*knut+0.4*knut*EXP(-1.1/(knut))) ![m/s]
          Re = grav_sed(b)%diameter*vt/(zvis(k)/roa(i,j,k,1))
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
          zsedl = ztempx*xn_adv(grav_sed(b)%spec,i,j,k)        ! kg kg-1, loss in terms of mixing ratio ! blir likt uavhengig av vt
          zflux(k) = zsedl*zdp1(k)                     ! --> [kg m-2 s-1]

          ! loss of mass in layer
          xn_adv(grav_sed(b)%spec,i,j,k) = xn_adv(grav_sed(b)%spec,i,j,k) - zsedl !  kg kg-1
        end do

        ! teste å gjøre det i en ny k-loop så det ikke blir så effektivt

        do  k = 1,KMAX_MID-1
          ! "arrival" of sedimented mass in box below
          if(k<KMAX_MID)then
            xn_adv(grav_sed(b)%spec,i,j,k+1) =  xn_adv(grav_sed(b)%spec,i,j,k+1) +  zflux(k)/zdp1(k+1)
          end if
        end do

        n=grav_sed(b)%n3d_gravset
        if(n>0)&
          d_3d(n,i,j,:,IOU_INST)=num_sed(lev3d(:num_lev3d))

      end do
    end do
  end do

! Multilayer crossing is here no ralised!!!
! sedimentation velocity is in effect limited to z/delt
! sedimentation to the ground from first layer sflx --> [kg m-2 s-1]
! sflux = zflux(:,1) ! sedimenterer ikke på bakken, bare sender det til nederste laget
end subroutine gravset


end module Gravset_mod
