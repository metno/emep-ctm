! <N2O5_hydrolysis_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________!
 module N2O5_hydrolysis_ml

! N2O5 -> nitrate calculation
!===========================================================================
! N2O5 -> nitrate calculation. Some constants for
! calculation of volume fraction of sulphate aerosol, and rate of uptake
! 
!
! The first order reaction coefficient K (corrected for gas phase diffusion, 
! Schwartz, 1986) is given by
!
! K= A* alpha* v/4
!    alpha=sticking coeff. for N2O5 =0.02
!    v=mean molecular speed for N2O5
!    A=aerosol surfac
!
! The surface area of the aerosols can be calculated as
! 
! A = V * surface/volume of aerosols
!     V=volume fraction of sulphate (cm3 aerosol/cm3 air)
!     (similar for nitrate and ammonium):
!
!     e.g.
!     V = (so4 in moleculescm-3) x atw sulphate
!         ---------------------------------------------------------
!        AVOG X specific density of aerosols (assumed 2g/cm3*rh correction)
!
!    Or, shorter, V = S x M0/(AVOG*rho)
!
!    where S is conc. e.g. sulphate (molecule/cm3), M0 is molwt. 
!
!
!    We do not want to include  concentrations  or rho yet, so:
!
!     Let VOL =  M0/AVOG
!   
! The surface/volume ratio is calculated using Whitby particle distribution
! with number mean radius 0.034  and standars deviation (Sigma)=2. 
! Then surface/volume=3/r *  exp( -5/2 *(lnSigma)^2)=26.54 
! 3* exp( -5/2 *(lnSigma)^2)=0.90236
! (monodisperse aerosols; 4*pi*r^2/(4/3 pi*r^3)= 3/r =88.2)
!
! Then 
!      A = VOL * S * 0.90236 /(0.034e-6*rho) 
! and
!      K = VOL * S * 0.90236 /(0.034e-6*rho)    * alpha* v/4
! Set
!      VOLFAC= VOL*0.90236/0.034e-6 *alpha    
! Then
!      K = VOLFAC *S *v/(4*rho)
!
! rcmisc(8,k) (in My_Chem)=v/(4*rho) 
!
!      K = VOLFAC *rcmisc(8,k) *S
! According to Riemer et al, 2003, we weight the reaction probability
! according to the composition of the aerosol
!
! alpha(N2O5)=f*alpha1 +(1-f)alpha2
!   alpha1=0.02
!   alpha2=0.002
!   f= Mso4/(Mso4+Mno3), M=aerosol mass concentration
 
 use ModelConstants_ml,     only :  KMAX_MID, KCHEMTOP
 use PhysicalConstants_ml,  only :AVOG

  implicit none
  private



   real, public, dimension(KCHEMTOP:KMAX_MID), save :: &
         f_Riemer              ! weighting factor for N2O5 hydrolysis
                                ! Mass of sulfate relative to sulfate+nitrate
                                ! according to  Riemer N, Vogel H, Vogel B, 
                                ! Schell B, Ackermann I, Kessler C, Hass H
                                ! JGR 108 (D4): FEB 27 2003 




  real, parameter, public  :: VOLFACSO4 = 96.0/(AVOG) * 0.90236 *0.02/0.034e-6 
  real, parameter, public  :: VOLFACNO3 = 62.0/(AVOG) * 0.90236 *0.02/0.034e-6 
  real, parameter, public  :: VOLFACNH4 = 18.0/(AVOG) * 0.90236 *0.02/0.034e-6 



 end module N2O5_hydrolysis_ml
