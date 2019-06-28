! <Convection_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Convection_mod
!If subsidence is included, ps3d could actually be any constant (in k), and 
! will recover the same value after the subsidience.
!
!In the global report 2010 ps3d=1, and xn_adv is the mixing ratio
!
!We define the mass of pollutant in a level to be proportionnal to xn_adv*dp(k)
!Note that dp(k) is defined as fixed under convection, i.e. the amount of air 
!and dp are not consistent during convection, but rather in a transition state. 
!In principle we can imagine that the convection is called many times with a 
!small dt_conv, but on the other hand this will give instantaneous mixing 
!between convection and subsidience, which is not corrrect either
!-----------------------------------------------------------------------------------

use Chemfields_mod,        only: xn_adv
use ChemDims_mod,          only: NSPEC_ADV
use Config_module,         only: KMAX_BND,KMAX_MID,PT,Pref
use MetFields_mod ,        only: ps,SigmaKz,u_xmj,v_xmi,cnvuf,cnvdf
use GridValues_mod,        only: dA, dB, sigma_bnd
use Par_mod,               only: LIMAX,LJMAX,limax,ljmax,li0,li1,lj0,lj1
use PhysicalConstants_mod, only: GRAV

public :: convection_pstar!in sigma coordinates
public :: convection_Eta!in more general hybrid coordinates (Eta)

contains

subroutine convection_pstar(ps3d,dt_conv)

  implicit none
  real ,intent(inout):: ps3d(LIMAX,LJMAX,KMAX_MID),dt_conv
 
  real ::xn_in_core(NSPEC_ADV,KMAX_MID+1)
  real ::mass_air_grid(KMAX_MID),mass_air_grid0(KMAX_MID), mass_air_core(KMAX_MID)
  real ::mass_exchanged,mass
  real :: mass_air_grid_k_temp,xn_buff(NSPEC_ADV,KMAX_MID)
  real :: dk(KMAX_MID),dp(KMAX_MID),totdk
  integer ::k,i,j,k_fill,k1
  


  do k=1,KMAX_MID
     dk(k)=sigma_bnd(k+1)-sigma_bnd(k)
  end do
  totdk=sigma_bnd(KMAX_MID+1)-sigma_bnd(1)!=1 in sigma coordinates

!UPWARD

    do j=lj0,lj1
      do i=li0,li1
         xn_in_core = 0.0!concentration null below surface
         mass=0.0
!
!ps3d=PS-PT=dp/dksi if sigma coordinates
!
!mass(k)=dp/dksi *dksi/g
!
         do k=1,KMAX_MID
            mass=mass+ps3d(i,j,k)*dk(k)
         end do
         mass=mass/totdk

         do k=1,KMAX_MID
            dp(k)=dA(k)+dB(k)*ps(i,j,1)
         end do

         mass_air_core=0.0
         do k=KMAX_MID,1,-1
            !-- mass_air=(dp/g)*gridarea 
            
            mass_air_grid0(k) =mass ! average density (in /dksi unit). Used only if subsidience is included
            mass_air_grid(k) = ps3d(i,j,k)! start density = dp/dksi where ksi is the vertical coordinate


            k1=k+1
            k1=min(k1,KMAX_MID)
!mass moves from cell k_1 to k_2:
!density (mass_air_core and xn_in_core) changes with a factor of dp(k_2)/dp(k_1)

            if(k<KMAX_MID)then
               mass_air_core(k)=mass_air_core(k1)*dp(k+1)/dp(k)!flux from below
               mass_air_core(k+1)=0.0
               xn_in_core(:,k) =xn_in_core(:,k1)*dp(k+1)/dp(k)!flux from below
               xn_in_core(:,k+1) =0.0
            end if

!fraction of grid moved to core:
! df/(dp/g)   df=horizontal flux   dp/g= total mass (/m2) in grid

!fraction og core moved to grid:
! df/f1  f1=total mass in core df=part which is exchanged horizontally

!horizontal flux
            if(cnvuf(i,j,k+1)-cnvuf(i,j,k)<=0.0)then
               !mass from grid to core - horizontal exchange
               mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))*GRAV*dt_conv/dp(k)*mass_air_grid(k)
               if(mass_exchanged<-mass_air_grid(k))then
                  !limit fluxes
                  cnvuf(i,j,k+1)=0.99*dp(k)/(GRAV*dt_conv)+cnvuf(i,j,k+1)!0.99 to determine
                  mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))*GRAV*dt_conv/dp(k)*mass_air_grid(k)
               end if
            else
               !mass from core to grid - horizontal exchange
               mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*mass_air_core(k)

             end if

            !horizontal exchange
            if(cnvuf(i,j,k+1)-cnvuf(i,j,k)<=0.0)then

               !mass from grid to core - horizontal exchange
               !NB change xn_in_core before xn_adv

              xn_in_core(:,k) = xn_in_core(:,k)-(cnvuf(i,j,k+1)-cnvuf(i,j,k))*xn_adv(:,i,j,k)*GRAV*dt_conv/dp(k)
               mass_air_core(k)=mass_air_core(k)-mass_exchanged  
               xn_adv(:,i,j,k)=xn_adv(:,i,j,k)+(cnvuf(i,j,k+1)-cnvuf(i,j,k))*xn_adv(:,i,j,k)*GRAV*dt_conv/dp(k)
               mass_air_grid(k) = mass_air_grid(k)+mass_exchanged
            else
               !NB change xn_adv before xn_in_core 
               xn_adv(:,i,j,k)=xn_adv(:,i,j,k)+(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*xn_in_core(:,k)
               xn_in_core(:,k) = xn_in_core(:,k)-(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*xn_in_core(:,k)
               mass_air_core(k)=mass_air_core(k)-mass_exchanged  
               mass_air_grid(k) = mass_air_grid(k)+mass_exchanged
            end if

        end do

!DOWNWARD
      if(.true.)then

         mass_air_core=0.0
         xn_in_core=0.0
         do k=1,KMAX_MID
            !-- mass_air=(dp/g)*gridarea 
            
            k1=k+1
            k1=min(k1,KMAX_MID)
            xn_in_core(:,k) = 0.0

            !vertical exchange
            if(k>1)then
               mass_air_core(k)=mass_air_core(k-1)*dp(k-1)/dp(k)!flux from above
               mass_air_core(k-1)=0.0
               xn_in_core(:,k) = xn_in_core(:,k-1)*dp(k-1)/dp(k)!flux from above
               xn_in_core(:,k-1) =0.0
            end if

            if(cnvdf(i,j,k+1)-cnvdf(i,j,k)<=0.0)then
               !mass from grid to core - horizontal exchange
               mass_exchanged=(cnvdf(i,j,k+1)-cnvdf(i,j,k))*mass_air_grid(k)*GRAV*dt_conv/dp(k)
               if(mass_exchanged<-mass_air_grid(k))then
               !limit fluxes
                  cnvdf(i,j,k+1)=-0.99*dp(k)/(GRAV*dt_conv)+cnvdf(i,j,k)!0.99 to determine
                  mass_exchanged=(cnvdf(i,j,k+1)-cnvdf(i,j,k))*mass_air_grid(k)*GRAV*dt_conv/dp(k)
               end if
            else
!NB: cnvdf < 0
               mass_exchanged=-(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*mass_air_core(k)
            end if

            !horizontal exchange
!NB: cnvdf < 0
            if(cnvdf(i,j,k+1)-cnvdf(i,j,k)<=0.0)then
               !mass from grid to core - horizontal exchange
               !NB change xn_in_core before xn_adv
               xn_in_core(:,k) = xn_in_core(:,k)-(cnvdf(i,j,k+1)-cnvdf(i,j,k))*xn_adv(:,i,j,k)*GRAV*dt_conv/dp(k)
               mass_air_core(k)=mass_air_core(k)-mass_exchanged  
               xn_adv(:,i,j,k)= xn_adv(:,i,j,k)+(cnvdf(i,j,k+1)-cnvdf(i,j,k))*xn_adv(:,i,j,k)*GRAV*dt_conv/dp(k)
               mass_air_grid(k) = mass_air_grid(k)+mass_exchanged
            else
               !mass from core to grid - horizontal exchange
               !NB change xn_adv before xn_in_core 
               xn_adv(:,i,j,k) = xn_adv(:,i,j,k)-(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*xn_in_core(:,k)
               xn_in_core(:,k) = xn_in_core(:,k)+(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*xn_in_core(:,k)
               mass_air_core(k) = mass_air_core(k)-mass_exchanged  

               mass_air_grid(k) = mass_air_grid(k)+mass_exchanged

            end if

         end do

         end if

         if(.true.)then
!diffusion free method
!distribute mass among level starting from top
!Allows "CFL number" larger than 1

         k_fill=1
         do k=1,KMAX_MID
            mass_air_grid_k_temp=0.0
            !fill level k with available mass
            !put new xn_adv in xn_buff because xn_adv should not be changed while still used 
            xn_buff(:,k) = 0.0
            do while (mass_air_grid_k_temp +mass_air_grid(k_fill)*dk(k_fill) <mass_air_grid0(k)*dk(k).and.k_fill<KMAX_MID)
               xn_buff(:,k) =  xn_buff(:,k)+ xn_adv(:,i,j,k_fill)*dk(k_fill)
               xn_adv(:,i,j,k_fill) =  0.0
               mass_air_grid_k_temp=mass_air_grid_k_temp+mass_air_grid(k_fill)*dk(k_fill)
               mass_air_grid(k_fill)=mass_air_grid(k_fill)-mass_air_grid(k_fill)!ZERO
               k_fill=k_fill+1            
            end do

            xn_buff(:,k)=xn_buff(:,k)+ xn_adv(:,i,j,k_fill)*dk(k_fill)*&
                 (mass_air_grid0(k)*dk(k)-mass_air_grid_k_temp)/(mass_air_grid(k_fill)*dk(k_fill))

            xn_adv(:,i,j,k_fill) = xn_adv(:,i,j,k_fill)- xn_adv(:,i,j,k_fill)*&
                 (mass_air_grid0(k)*dk(k)-mass_air_grid_k_temp)/(mass_air_grid(k_fill)*dk(k_fill))

            mass_air_grid(k_fill)=mass_air_grid(k_fill)-&
                 (mass_air_grid0(k)*dk(k)-mass_air_grid_k_temp)/dk(k_fill)

            ps3d(i,j,k) = mass_air_grid0(k)!=(mass_air_grid_k_temp+(mass_air_grid0(k)*dk(k)-mass_air_grid_k_temp))/dk(k)
            
        end do
        do k=1,KMAX_MID
           xn_adv(:,i,j,k)=xn_buff(:,k)/dk(k)
        end do
!check that all mass is distributed
!         if(abs(mass_air_grid(k_fill))>1.0.or.k_fill/=KMAX_MID)then
!            if(ME==0)write(*,*)'ERRORMASS',ME,i,j,k_fill,mass_air_grid(k_fill),mass_air_grid_k_temp,mass_air_grid0(k)
!         end if
         else
            do k=1,KMAX_MID
               ps3d(i,j,k) = mass_air_grid(k)
            end do

         end if

      end do
   end do

 end subroutine convection_pstar
 subroutine convection_Eta(dpdeta,dt_conv)

    implicit none
    real ,intent(inout):: dpdeta(LIMAX,LJMAX,KMAX_MID),dt_conv
   
    real ::xn_in_core(NSPEC_ADV,KMAX_MID+1)
    real ::mass_air_grid(KMAX_MID),mass_air_grid0(KMAX_MID), mass_air_core(KMAX_MID)
    real ::mass_exchanged,mass
    real :: mass_air_grid_k_temp,xn_buff(NSPEC_ADV,KMAX_MID)
    real :: dk(KMAX_MID),dp(KMAX_MID),totdk
    integer ::k,i,j,k_fill,k1
   INCLUDE 'mpif.h'
  

   totdk=0.0
   do k=1,KMAX_MID
      dk(k)=dA(k)/Pref+dB(k)
      totdk=totdk+dk(k)
   end do

!UPWARD

    do j=lj0,lj1
      do i=li0,li1
         xn_in_core = 0.0!concentration null below surface
         mass=0.0
!
!ps3d=PS-PT=dp/dksi if sigma coordinates
!
!mass(k)=dp/dksi *dksi/g
!
         do k=1,KMAX_MID
            mass=mass+dpdeta(i,j,k)*dk(k)
         end do
         mass=mass/totdk

         do k=1,KMAX_MID
            dp(k)=dA(k)+dB(k)*ps(i,j,1)
         end do

         mass_air_core=0.0
         do k=KMAX_MID,1,-1
            !-- mass_air=(dp/g)*gridarea 
            
            mass_air_grid0(k) =mass ! average density (in /dksi unit). Used only if subsidience is included
            mass_air_grid(k) = dpdeta(i,j,k)! start density = dp/dksi where ksi is the vertical coordinate


            k1=k+1
            k1=min(k1,KMAX_MID)
!mass moves from cell k_1 to k_2:
!density (mass_air_core and xn_in_core) changes with a factor of dp(k_2)/dp(k_1)

            if(k<KMAX_MID)then
               mass_air_core(k)=mass_air_core(k1)*dp(k+1)/dp(k)!flux from below
               mass_air_core(k+1)=0.0
               xn_in_core(:,k) =xn_in_core(:,k1)*dp(k+1)/dp(k)!flux from below
               xn_in_core(:,k+1) =0.0
            end if

!fraction of grid moved to core:
! df/(dp/g)   df=horizontal flux   dp/g= total mass (/m2) in grid

!fraction og core moved to grid:
! df/f1  f1=total mass in core df=part which is exchanged horizontally

!horizontal flux
            if(cnvuf(i,j,k+1)-cnvuf(i,j,k)<=0.0)then
               !mass from grid to core - horizontal exchange
               mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))*GRAV*dt_conv/dp(k)*mass_air_grid(k)
               if(mass_exchanged<-mass_air_grid(k))then
                  !limit fluxes
                  cnvuf(i,j,k+1)=0.99*dp(k)/(GRAV*dt_conv)+cnvuf(i,j,k+1)!0.99 to determine
                  mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))*GRAV*dt_conv/dp(k)*mass_air_grid(k)
               end if
            else
               !mass from core to grid - horizontal exchange
               mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*mass_air_core(k)

             end if

            !horizontal exchange
            if(cnvuf(i,j,k+1)-cnvuf(i,j,k)<=0.0)then

               !mass from grid to core - horizontal exchange
               !NB change xn_in_core before xn_adv

              xn_in_core(:,k) = xn_in_core(:,k)-(cnvuf(i,j,k+1)-cnvuf(i,j,k))*xn_adv(:,i,j,k)*GRAV*dt_conv/dp(k)
               mass_air_core(k)=mass_air_core(k)-mass_exchanged  
               xn_adv(:,i,j,k)=xn_adv(:,i,j,k)+(cnvuf(i,j,k+1)-cnvuf(i,j,k))*xn_adv(:,i,j,k)*GRAV*dt_conv/dp(k)
               mass_air_grid(k) = mass_air_grid(k)+mass_exchanged
            else
               !NB change xn_adv before xn_in_core 
               xn_adv(:,i,j,k)=xn_adv(:,i,j,k)+(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*xn_in_core(:,k)
               xn_in_core(:,k) = xn_in_core(:,k)-(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*xn_in_core(:,k)
               mass_air_core(k)=mass_air_core(k)-mass_exchanged  
               mass_air_grid(k) = mass_air_grid(k)+mass_exchanged
            end if

        end do

!DOWNWARD
      if(.true.)then

         mass_air_core=0.0
         xn_in_core=0.0
         do k=1,KMAX_MID
            !-- mass_air=(dp/g)*gridarea 
            
            k1=k+1
            k1=min(k1,KMAX_MID)
            xn_in_core(:,k) = 0.0

            !vertical exchange
            if(k>1)then
               mass_air_core(k)=mass_air_core(k-1)*dp(k-1)/dp(k)!flux from above
               mass_air_core(k-1)=0.0
               xn_in_core(:,k) = xn_in_core(:,k-1)*dp(k-1)/dp(k)!flux from above
               xn_in_core(:,k-1) =0.0
            end if

            if(cnvdf(i,j,k+1)-cnvdf(i,j,k)<=0.0)then
               !mass from grid to core - horizontal exchange
               mass_exchanged=(cnvdf(i,j,k+1)-cnvdf(i,j,k))*mass_air_grid(k)*GRAV*dt_conv/dp(k)
               if(mass_exchanged<-mass_air_grid(k))then
               !limit fluxes
                  cnvdf(i,j,k+1)=-0.99*dp(k)/(GRAV*dt_conv)+cnvdf(i,j,k)!0.99 to determine
                  mass_exchanged=(cnvdf(i,j,k+1)-cnvdf(i,j,k))*mass_air_grid(k)*GRAV*dt_conv/dp(k)
               end if
            else
!NB: cnvdf < 0
               mass_exchanged=-(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*mass_air_core(k)
            end if

            !horizontal exchange
!NB: cnvdf < 0
            if(cnvdf(i,j,k+1)-cnvdf(i,j,k)<=0.0)then
               !mass from grid to core - horizontal exchange
               !NB change xn_in_core before xn_adv
               xn_in_core(:,k) = xn_in_core(:,k)-(cnvdf(i,j,k+1)-cnvdf(i,j,k))*xn_adv(:,i,j,k)*GRAV*dt_conv/dp(k)
               mass_air_core(k)=mass_air_core(k)-mass_exchanged  
               xn_adv(:,i,j,k)= xn_adv(:,i,j,k)+(cnvdf(i,j,k+1)-cnvdf(i,j,k))*xn_adv(:,i,j,k)*GRAV*dt_conv/dp(k)
               mass_air_grid(k) = mass_air_grid(k)+mass_exchanged
            else
               !mass from core to grid - horizontal exchange
               !NB change xn_adv before xn_in_core 
               xn_adv(:,i,j,k) = xn_adv(:,i,j,k)-(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*xn_in_core(:,k)
               xn_in_core(:,k) = xn_in_core(:,k)+(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*xn_in_core(:,k)
               mass_air_core(k) = mass_air_core(k)-mass_exchanged  

               mass_air_grid(k) = mass_air_grid(k)+mass_exchanged

            end if

         end do

         end if

         if(.true.)then
!diffusion free method
!distribute mass among level starting from top
!Allows "CFL number" larger than 1

         k_fill=1
         do k=1,KMAX_MID
            mass_air_grid_k_temp=0.0
            !fill level k with available mass
            !put new xn_adv in xn_buff because xn_adv should not be changed while still used 
            xn_buff(:,k) = 0.0
            do while (mass_air_grid_k_temp +mass_air_grid(k_fill)*dk(k_fill) <mass_air_grid0(k)*dk(k).and.k_fill<KMAX_MID)
               xn_buff(:,k) =  xn_buff(:,k)+ xn_adv(:,i,j,k_fill)*dk(k_fill)
               xn_adv(:,i,j,k_fill) =  0.0
               mass_air_grid_k_temp=mass_air_grid_k_temp+mass_air_grid(k_fill)*dk(k_fill)
               mass_air_grid(k_fill)=mass_air_grid(k_fill)-mass_air_grid(k_fill)!ZERO
               k_fill=k_fill+1            
            end do

            xn_buff(:,k)=xn_buff(:,k)+ xn_adv(:,i,j,k_fill)*dk(k_fill)*&
                 (mass_air_grid0(k)*dk(k)-mass_air_grid_k_temp)/(mass_air_grid(k_fill)*dk(k_fill))

            xn_adv(:,i,j,k_fill) = xn_adv(:,i,j,k_fill)- xn_adv(:,i,j,k_fill)*&
                 (mass_air_grid0(k)*dk(k)-mass_air_grid_k_temp)/(mass_air_grid(k_fill)*dk(k_fill))

            mass_air_grid(k_fill)=mass_air_grid(k_fill)-&
                 (mass_air_grid0(k)*dk(k)-mass_air_grid_k_temp)/dk(k_fill)

            dpdeta(i,j,k) = mass_air_grid0(k)!=(mass_air_grid_k_temp+(mass_air_grid0(k)*dk(k)-mass_air_grid_k_temp))/dk(k)
            
        end do
        do k=1,KMAX_MID
           xn_adv(:,i,j,k)=xn_buff(:,k)/dk(k)
        end do
!check that all mass is distributed
!         if(abs(mass_air_grid(k_fill))>1.0.or.k_fill/=KMAX_MID)then
!            if(ME==0)write(*,*)'ERRORMASS',ME,i,j,k_fill,mass_air_grid(k_fill),mass_air_grid_k_temp,mass_air_grid0(k)
!         end if
         else
            do k=1,KMAX_MID
               dpdeta(i,j,k) = mass_air_grid(k)
            end do

         end if

      end do
   end do

 end subroutine convection_Eta
end module Convection_mod
