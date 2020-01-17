! <Convection_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.34>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2020 met.no
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
use MetFields_mod ,        only: ps,cnvuf,cnvdf
use GridValues_mod,        only: dA, dB, sigma_bnd
use Par_mod,               only: LIMAX,LJMAX,limax,ljmax,li0,li1,lj0,lj1
use PhysicalConstants_mod, only: GRAV

public :: convection_Eta!in more general hybrid coordinates (Eta)

contains

!_________________________________________________________________________________________
  subroutine convection_Eta(dt_conv)
!Pollutants are conserved. We put a lock over top. Pollutants do not move out of the top.
!Assumes xn_adv in units of mixing ratio
!See also EMEP technical report MSC-W 1/2010 p. 11-14

!xn*dp(k) = mass in level k (up to a fixed constant *g*gridarea) 

    implicit none
    real ,intent(inout):: dt_conv

    real ::xn_in_core(NSPEC_ADV,KMAX_MID+1)
    real ::xn_air_grid(KMAX_MID),mass_air_grid0(KMAX_MID), xn_air_core(KMAX_MID)
    real ::mass_exchanged, total, totalmass,x
    real :: mass_air_grid_k_temp,xn_buff(NSPEC_ADV,KMAX_MID)
    real :: dp(KMAX_MID)
    integer ::k,i,j,k_fill,k1,kk,n
    logical, parameter ::masstest=.false.!normally, the mass should be conserved mathematically. Mostly to make demo.
    INCLUDE 'mpif.h'

    n=3 !IXADV of species to test for mass balance (masstest)
    !UPWARD

    do j=lj0,lj1
       do i=li0,li1

          !-- mass_air=(dp/g)*gridarea 

          totalmass = 0.0!total colomn mass of species n
          do k=1,KMAX_MID
             dp(k)=dA(k)+dB(k)*ps(i,j,1)
             xn_air_grid(k) = 1.0 !unit mixing ratio
             mass_air_grid0(k) =xn_air_grid(k)*dp(k) ! original mass of air 
             if(masstest)totalmass = totalmass + xn_adv(n,i,j,k)*dp(k)!mass of species n (without g*gridarea factor)
          end do

!air is moved horizontally between a "core" and the "grid". The air from core moves upward. 
          xn_air_core = 0.0
          xn_in_core = 0.0
          do k=KMAX_MID,1,-1

             !mass moves from cell k_1 to k_2:
             !mixing ratios (xn_air_core and xn_in_core) change with a factor of dp(k_1)/dp(k_2)

             if(k<KMAX_MID)then
                !xn_air_core(k+1) is defined at the preceding k
                xn_air_core(k)=xn_air_core(k+1)*dp(k+1)/dp(k)!flux from below
                xn_air_core(k+1)=0.0
                xn_in_core(:,k) =xn_in_core(:,k+1)*dp(k+1)/dp(k)!flux from below
                xn_in_core(:,k+1) =0.0
             end if

             !fraction of grid moved to core:
             ! df/(dp/g)   df=horizontal flux   dp/g= total mass (/m2) in grid

             !fraction og core moved to grid:
             ! df/f1  f1=total mass in core df=part which is exchanged horizontally

             !horizontal flux
             if(cnvuf(i,j,k+1)-cnvuf(i,j,k)<=0.0)then
                !mass from grid to core - horizontal exchange
                mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))*GRAV*dt_conv/dp(k)*xn_air_grid(k)!in units of fraction of cell content
                if(mass_exchanged<-xn_air_grid(k))then
                   !limit fluxes
                   cnvuf(i,j,k+1)=-0.9*dp(k)/(GRAV*dt_conv)+cnvuf(i,j,k)!0.9 to determine
                   mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))*GRAV*dt_conv/dp(k)*xn_air_grid(k)
                end if
             else
                !mass from core to grid - horizontal exchange
                mass_exchanged=(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*xn_air_core(k)               
             end if

             !horizontal exchange
             if(cnvuf(i,j,k+1)-cnvuf(i,j,k)<=0.0)then
                !mass from grid to core - horizontal exchange
                !NB change xn_in_core before xn_adv
                xn_in_core(:,k) = xn_in_core(:,k)-mass_exchanged*xn_adv(:,i,j,k)
                xn_air_core(k)=xn_air_core(k)-mass_exchanged  
                xn_adv(:,i,j,k)=xn_adv(:,i,j,k)+mass_exchanged*xn_adv(:,i,j,k)
                xn_air_grid(k) = xn_air_grid(k)+mass_exchanged
             else
                !NB change xn_adv before xn_in_core 
                xn_adv(:,i,j,k)=xn_adv(:,i,j,k)+(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*xn_in_core(:,k)
                xn_in_core(:,k) = xn_in_core(:,k)-(cnvuf(i,j,k+1)-cnvuf(i,j,k))/cnvuf(i,j,k+1)*xn_in_core(:,k)
                xn_air_core(k)=xn_air_core(k)-mass_exchanged  
                xn_air_grid(k) = xn_air_grid(k)+mass_exchanged
             end if

          end do

          !we put the air left in xn_in_core(:,k=2) into the top level:
          xn_adv(:,i,j,1)=xn_adv(:,i,j,1)+xn_in_core(:,1)
          xn_in_core(:,1)=0.0

          if(masstest)then
             !check that total mass of air is unchanged
             total = 0.0
             do k=1,KMAX_MID
                total = total + xn_adv(n,i,j,k)*dp(k)
                if(k>1 .and. abs(xn_in_core(n,k))>1.E-20)write(*,*)'ERROR CORE',i,j,k,xn_in_core(n,k)
             enddo
             if(abs(total - totalmass)>1.E-8)write(*,*)'ERROR MASS up',i,j,total, totalmass
          end if


          !DOWNWARD
          xn_air_core = 0.0
          xn_in_core = 0.0
          do k=1,KMAX_MID
             !-- mass_air=(dp/g)*gridarea             

             xn_in_core(:,k) = 0.0

             !vertical exchange
             if(k>1)then
                xn_air_core(k)=xn_air_core(k-1)*dp(k-1)/dp(k)!flux from above
                xn_air_core(k-1)=0.0
                xn_in_core(:,k) = xn_in_core(:,k-1)*dp(k-1)/dp(k)!flux from above
                xn_in_core(:,k-1) =0.0
             end if

             if(cnvdf(i,j,k+1)-cnvdf(i,j,k)<=0.0)then
                !mass from grid to core - horizontal exchange
                mass_exchanged=(cnvdf(i,j,k+1)-cnvdf(i,j,k))*xn_air_grid(k)*GRAV*dt_conv/dp(k)
                if(mass_exchanged<-xn_air_grid(k))then
                   !limit fluxes
                   cnvdf(i,j,k) = cnvdf(i,j,k+1)+0.9*dp(k)/(GRAV*dt_conv)!0.9 to determine
                   mass_exchanged=(cnvdf(i,j,k+1)-cnvdf(i,j,k))*xn_air_grid(k)*GRAV*dt_conv/dp(k)
                end if
             else
                !NB: cnvdf < 0
                mass_exchanged=-(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*xn_air_core(k)
             end if

             !horizontal exchange
             !NB: cnvdf < 0
             if(cnvdf(i,j,k+1)-cnvdf(i,j,k)<=0.0)then
                !mass from grid to core - horizontal exchange
                !NB change xn_in_core before xn_adv
                xn_in_core(:,k) = xn_in_core(:,k)-mass_exchanged*xn_adv(:,i,j,k)
                xn_air_core(k)=xn_air_core(k)-mass_exchanged  
                xn_adv(:,i,j,k)= xn_adv(:,i,j,k)+mass_exchanged*xn_adv(:,i,j,k)
                xn_air_grid(k) = xn_air_grid(k)+mass_exchanged
             else
                !mass from core to grid - horizontal exchange
                !NB change xn_adv before xn_in_core 
                xn_adv(:,i,j,k) = xn_adv(:,i,j,k)-(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*xn_in_core(:,k)
                xn_in_core(:,k) = xn_in_core(:,k)+(cnvdf(i,j,k+1)-cnvdf(i,j,k))/cnvdf(i,j,k)*xn_in_core(:,k)
                xn_air_core(k) = xn_air_core(k)-mass_exchanged  
                xn_air_grid(k) = xn_air_grid(k)+mass_exchanged
             end if

          end do
          if(masstest)then
             !check that total mass of air is unchanged
             total = 0.0
             do k=1,KMAX_MID
                total = total + xn_adv(n,i,j,k)*dp(k)
                if(k>1 .and. abs(xn_in_core(n,k))>1.E-20)write(*,*)'ERROR CORE down',i,j,k,xn_in_core(n,k)
             enddo
             if(abs(total - totalmass)>1.E-8)write(*,*)'ERROR MASS down',i,j,total, totalmass
          endif
!end DOWNWARD

!Air have been moved between vertical levels, but the net quantity of air in each level has changed.
!the air is "relaxed" -> subsidence

!SUBSIDENCE:

          !diffusion free method
          !distribute mass among levels starting from top
          !Allows "CFL number" larger than 1

          ! we go from top to bottom, and stop when the mass_air_grid_k_temp = mass_air_grid0(k)
          ! when xn_air_grid is moved vertically, a factor dp(k1)/dp(k2) must be applied to take into account difference in pressure and thickness

          !we take air from xn_air_grid and put into mass_air_grid_k_temp
          !we take pollutants from xn_adv and put into xn_buff

          k_fill = 1
          xn_buff = 0.0
          do k=1,KMAX_MID
             mass_air_grid_k_temp=0.0
             !fill level k with available mass
             !put new xn_adv in xn_buff because xn_adv should not be changed while still used 
             if(masstest)then
                total = 0.0
                do kk=1,KMAX_MID
                   total = total + xn_adv(n,i,j,kk)*dp(kk)+ xn_buff(n,kk)*dp(kk)
                enddo
                if(abs(total - totalmass)>1.E-8)then
                   write(*,*)'ERROR MASS subsidence',i,j,k,total, totalmass
                   stop
                endif
             endif
             do while (mass_air_grid_k_temp +xn_air_grid(k_fill)*dp(k_fill) <mass_air_grid0(k).and.k_fill<KMAX_MID)
                !fill with entire level k_fill
                mass_air_grid_k_temp=mass_air_grid_k_temp+xn_air_grid(k_fill)*dp(k_fill) !
                xn_air_grid(k_fill)=xn_air_grid(k_fill)-xn_air_grid(k_fill)!ZERO remove all air from level 
                xn_buff(:,k) =  xn_buff(:,k)+ xn_adv(:,i,j,k_fill)*dp(k_fill)/dp(k)
                xn_adv(:,i,j,k_fill) =  0.0 !now the pollutants have been removed from xn_adv, and stored into xn_buff
                k_fill=k_fill+1
             end do

             !fill with part of the level, just enough so that  xn_air_grid(k)=mass_air_grid0(k)
             mass_exchanged = mass_air_grid0(k)-mass_air_grid_k_temp

             mass_air_grid_k_temp = mass_air_grid_k_temp + mass_exchanged !just to show fluxes

             !fraction of mass from level k_fill that is moved to level k:
             !  (mass_air_grid0(k)-mass_air_grid_k_temp)/(xn_air_grid(k_fill)*dp(k_fill))
             !must apply a factor dp(k_fill)/dp(k) when moving mixing ratios
             ! -> (mass_air_grid0(k)-mass_air_grid_k_temp)/(xn_air_grid(k_fill)* dp(k))

             xn_buff(:,k)=xn_buff(:,k)+ xn_adv(:,i,j,k_fill)*&
                  mass_exchanged/(xn_air_grid(k_fill)*dp(k))

             xn_adv(:,i,j,k_fill) = xn_adv(:,i,j,k_fill)- xn_adv(:,i,j,k_fill)*&
                  mass_exchanged/(xn_air_grid(k_fill)*dp(k_fill))

             !(next line must be set at end, since xn_air_grid is used just above)
             xn_air_grid(k_fill)=xn_air_grid(k_fill)-mass_exchanged/dp(k_fill)

             if(masstest)then
                total = 0.0
                do kk=1,KMAX_MID
                   total = total + xn_adv(n,i,j,kk)*dp(kk)+ xn_buff(n,kk)*dp(kk)
                enddo
                if(abs(total - totalmass)>1.E-8)then
                   write(*,*)'ERROR AIRMASS subsidence',i,j,k,k_fill,total, totalmass,x
                   stop
                endif
             endif
          end do
          
          if(masstest)then
             !check that all mass is distributed into xn_buff. xn_air_grid and xn_adv should be empty
             if(abs(xn_air_grid(k_fill))>1.0E-12.or.k_fill/=KMAX_MID )then
                write(*,*)'ERROR AIR MASS',i,j,k_fill,xn_air_grid(k_fill),mass_air_grid_k_temp,mass_air_grid0(k_fill)
                stop
             end if
             do k=1,KMAX_MID
                if(xn_adv(n,i,j,k)>1.E-8)then
                   write(*,*)'ERROR MASS',i,j,k,xn_adv(n,i,j,k),xn_buff(n,k)
                   stop
                end if
             end do
          endif
          
!copy the end results in xn_adv          
          do k=1,KMAX_MID
             xn_adv(:,i,j,k)=xn_buff(:,k)
          end do
          
          if(masstest)then
             total = 0.0
             do k=1,KMAX_MID
                total = total + xn_adv(n,i,j,k)*dp(k)
             enddo
             if(abs(total - totalmass)>1.E-8)then
                write(*,*)'ERROR MASS final',i,j,total, totalmass
                stop
             endif
          end if
          
       end do
    end do

  end subroutine convection_Eta
end module Convection_mod
