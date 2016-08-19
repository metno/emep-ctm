! <Aero_Rb_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!==============================================================================
  module Aero_DryDep_ml
!==============================================================================

  ! DESCRIPTION
  ! Calculates laminar sub-layer resistance (rb) and gravitational settling 
  ! velocity (vs) for particles
  ! Finally: vd= vs+1/(ra+rb+ra*rb*vs),   where
  ! vs - gravitational settling velocity,
  ! ra - aerodynamic resistance, rb - viscous sub-layer resistance,
  ! rc - surface resistance (assumed zero for particles)
  !---------------------------------------------------------------------------
 use My_Aerosols_ml,        only : NSIZE
 use LandDefs_ml,       only:  LandType
 use PhysicalConstants_ml , only : PI, GRAV, KARMAN, VISCO, BOLTZMANN, FREEPATH     

 implicit none
  private

  public  :: Aero_Rb

contains

 !  ===========================================================
    subroutine Aero_Rb (ustar, conv, roa, v50, lu,     &  ! IN
                        snow, wetarea, tsK,            &  ! IN
                        vs, rb, rbw)                      ! OUT
 
  !-------------------------------------------------------------------
  ! Calculates size dependent laminar layer resistance and gravitation  
  ! settling for particles (presently only for (1)PM2.5 and (2)PM10 )
  !     rb  - over dry surface (bounce-off for coarse PM)
  !     rbw - over wet surface (rain last 3 hours)
  !----------------------------------------------------------- ------                       
                               
! 	parameter (VISCO=1.46e-5, BOLTZ=1.381e-23      &
!                 ,FREEPATH=0.065e-6, DENSPART=2.2e3 )

  integer, intent(in)  :: lu, snow 
  real, intent(in)     :: ustar, conv, roa, v50, tsK
  real, intent(in)     :: wetarea
  real, intent(out)    :: vs(NSIZE), rb(NSIZE)  & ! over dry surface
                                   , rbw(NSIZE)   ! over wet surface

 !== local
  real, parameter, dimension(NSIZE) ::   &
                                       diam   = (/ 0.3e-6, 4.0e-6 /)   &
                                     , sigma  = (/ 1.8, 2.0 /)         &
                                     , PMdens = (/ 1600., 2200. /)
  real, parameter ::   AHAT = 1.e-3  !! charact. "radius" of grass blades, needlies etc.
  integer :: imod 
  real  :: stdlog,sig, dg, knut,slip, &
                  Di1,Di, vind, &
                  stoke, schmidt, &  ! Stoke's and Schmidt numbers
                  vsmo, vs1,      &  ! Settling velocity
                   coleff, reb, convfac

  real,    save :: log10
  logical, save :: my_first_call = .true.

  if ( my_first_call ) then
     log10  = log(10.0)   !ds
     my_first_call = .false.
  end if

!================================================

  MODEloop: do imod = 1, NSIZE

	stdlog = log(sigma(imod))
	sig = stdlog * stdlog     ! (log(STD))^2

!... mass median diameter -> geometric diameter 

        dg = exp (log(diam(imod)) - 3.* sig )

	knut = 2.*FREEPATH/dg   ! Knut's number
!... slip correction coefficient  
!	slipmo= 1.+ knut*       &               ! for monodisperse
!                (1.257+0.4*exp(-1.1* /knut))
	slip =  1.+ 1.246*knut                  ! for polydisperse

!== monodisperse aerosols =====
!     Dimo =BOLTZMANN*tsK*slipmo/(3*PI*dg *VISCO*roa)        ! diffusion coefficient
!     vsmo =dg*dg *PMdens(imod) *GRAV*slipmo/(18.*VISCO*roa) ! gravitational settling

!== polydisperse aerosols (log-normal size distribution) =====
     Di1 =BOLTZMANN*tsK/(3*PI*dg *VISCO *roa)              
     Di = Di1*(exp(-2.5*sig)+1.246*knut*exp(-4.*sig))  ! diffusion coefficient

     vs1=dg*dg * PMdens(imod)*GRAV/18./VISCO/roa  
     vs(imod) = vs1*(exp(8.*sig)+1.246*knut*exp(3.5*sig))    ! gravitational settling

! -------------------------------------------------------------------

!// Stokes and Schmidt numbers:

   ! == monodisperse ======
	! STmo=vsmo*ustar*ustar/VISCO/GRAV
	! SCmo=VISCO/dimo
   ! == polydisperse ======
	 schmidt = VISCO/Di              ! Schmidt number
	 stoke = vs(imod)*ustar*ustar/VISCO/GRAV ! Stoke number(based on depth 
                                              ! of laminar layer)

 !// collection efficiency  =======================
 !      coleff=1./sc**(2./3.) + 1./10.**(3/stoke)
 !=================================================

         vind = max( 0.005, v50 )   ! wind at 50m height



      if( LandType(lu)%is_water )    then !//===  WATER surface  ( Slinn & Slinn, 1980 ) =

      	   coleff= ustar / (KARMAN * vind) *      &          ! polydisperse
                  (exp(-0.5*log(schmidt)) + exp(-3./stoke*log10) )   

      elseif ( LandType(lu)%is_conif )  then !//===  CONIFEROUS ==============

           stoke = vs(imod)*ustar/(AHAT*GRAV)   ! vegetation (Slinn, 1982)
           coleff= exp(-2./3.*log(schmidt)) + stoke/(1.+ stoke*stoke)  !  Slinn 

      elseif ( LandType(lu)%is_veg ) then !//===  other VEGETATIVE surfaces ======

          if ( snow > 0 .or. tsK <= 273.)   then   !... covered with snow or frozen 

              coleff= exp(-2./3.*log(schmidt)) + exp(-3./stoke*log10)  ! polydisperse  

          else   !... snowfree  

               stoke = vs(imod)*ustar/(AHAT*GRAV)     ! Stoke for vegetation (Slinn, 1982)
               coleff= exp(-2./3.*log(schmidt)) + stoke/(1.+ stoke*stoke)  !  Slinn 
          endif

      else    !//====  urban/desert/ice always ===============================
              !..   Slinn at al(1978), Seinfeld(1997), Binkowski 

           coleff= exp(-2./3.*log(schmidt)) + exp(-3./stoke*log10)  ! polydisperse

      endif        

! .. laminar layer resistance .....................................
    !  rb= 1./ustar/colef                  !     Seinfeld
    !  rb= 0.4*vind/(ustar*ustar*colef)      !     Slinn

!// ==  bounce-off for coarse particles (Slinn, 1982) ===============
   !... for fine aerosol

       reb = 1.
 
   !... for coarse aerosol 

      if(imod == NSIZE )  then

        if( .not. (LandType(lu)%is_water .and. wetarea == 0.0))  & ! not on water/wet surface

            reb = max (1.e-7, exp(-2. * sqrt(stoke)))           
      endif

!//== enhanced dry.dep under convective conditions (Wesely at al,1985) ====

      convfac = 0. 
     
!  only for (low) vegetation 

      if (LandType(lu)%is_veg .and. snow == 0)  convfac = conv  

!// == sub-laminar layer resistance ====================

     rb (imod)  = 1. /(ustar*(1.+ 0.24 *convfac)) /(coleff*reb)
     rbw (imod) = 1. /(ustar*(1.+ 0.24 *convfac)) /coleff   ! no bounce-off on
                                                            ! wet surfaces 

!..monodisperse: rb1= 1./(ustar*(1. + 0.24 * conv))/(colef1*reb) 

       end do MODEloop

   end subroutine Aero_Rb
! =================================================================

end module Aero_DryDep_ml

