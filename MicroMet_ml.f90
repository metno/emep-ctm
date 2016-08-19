! <MicroMet_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2011 met.no
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
module Micromet_ml
!____________________________________________________________________
! Miscellaneous collection of "standard" micromet functions
! Including PsiM, PsiH, AerRes
! Some based upon code from Juha-Pekka Tuovinen, based upon Garrett
!____________________________________________________________________
!
!** includes
!
!   Depends on: none - self-contained.
!   Language: F
!____________________________________________________________________
  implicit none
  !F private

 !/-- Micromet (Aerodynamic) routines

  public :: rh2vpd

  public :: AerRes

  public :: AerResM

  public :: PsiH

  public :: PsiM

  public :: wind_at_h   !wind for given height

  !/-- define PI here rather than use PhysicalCOnstants_ml, to
  !    preserve self-sufficiency

  real, public, parameter  ::    &
       PI      = 3.14159265358979312000   ! pi, from 4.0*atan(1.) on cray


  !========================================
  contains
  !========================================


 !=======================================================================
  !--------------------------------------------------------------------
  function rh2vpd(T,rh) result (vpd_res)
  !This function is not currently in use.

    real, intent(in) ::  T    ! Temperature (K)
    real, intent(in) ::  rh   ! relative humidity (%)
    real :: vpd_res     ! vpd   = water vapour pressure deficit (Pa)

    !   Local:
    real :: vpSat       ! vpSat = saturated water vapour pressure (Pa)
    real :: arg

    arg   = 17.67 * (T-273.15)/(T-29.65)
    vpSat = 611.2 * exp(arg)
    vpd_res   = (1.0 - rh/100.0) * vpSat

  end function rh2vpd

  !--------------------------------------------------------------------
  function AerRes(z1,z2,uStar,Linv,Karman) result (Ra)
!...
!   Ref: Garratt, 1994, pp.55-58
!   In:
    real, intent(in) ::   z1     ! lower height (m), equivalent to h-d+1 or h-d+3
    real, intent(in) ::   z2     ! upper height (m), equivalent to z-d
    real, intent(in) ::   uStar  ! friction velocity (m/s)
    real, intent(in) ::   Linv   ! inverse of the Obukhov length (1/m)
    
    real, intent(in) ::   Karman ! von Karman's constant 
!   For AerRes, the above dummy argument is replaced by the actual argument 
!   KARMAN in the module GetMet_ml.

!   Out:
    real :: Ra      ! =  aerodynamic resistance to transfer of sensible heat
                    !from z2 to z1 (s/m)

!   uses functions:
!   PsiH   = integral flux-gradient stability function for heat 
!...

    if ( z1 > z2 ) then
      Ra = -999.0
    else
      Ra = log(z2/z1) - PsiH(z2*Linv) + PsiH(z1*Linv)
      Ra = Ra/(Karman*uStar)
    end if

  end function AerRes

  !--------------------------------------------------------------------
  function AerResM(z1,z2,uStar,Linv,Karman) result (Ra)
!...
!   Ref: Garratt, 1994, pp.55-58
!   In:
    real, intent(in) ::   z1     ! lower height (m), equivalent to h-d+1 or h-d+3
    real, intent(in) ::   z2     ! upper height (m), equivalent to z-d
    real, intent(in) ::   uStar  ! friction velocity (m/s)
    real, intent(in) ::   Linv   ! inverse of the Obukhov length (1/m)
    
    real, intent(in) ::   Karman ! von Karman's constant 
!   For AerRes, the above dummy argument is replaced by the actual argument 
!   KARMAN in the module GetMet_ml.

!   Out:
    real :: Ra     ! =  aerodynamic resistance to transfer of momentum
                    !from z2 to z1 (s/m)

!   uses functions:
!   PsiM   = integral flux-gradient stability function for momentum
!...

    Ra = log(z2/z1) - PsiM(z2*Linv) + PsiM(z1*Linv)
    Ra = Ra/(Karman*uStar)

  end function AerResM

  !--------------------------------------------------------------------
  function PsiH(zL) result (stab_h)
    !  PsiH = integral flux-gradient stability function for heat 
    !  Ref: Garratt, 1994, pp52-54

    ! In:
    real, intent(in) :: zL   ! surface layer stability parameter, (z-d)/L 
    
    ! Out:
    real :: stab_h         !   PsiH(zL) 
    
   ! Local
   real :: x
 
    if (zL <  0) then !unstable
        x    = sqrt(1.0 - 16.0 * zL)
        stab_h = 2.0 * log( (1.0 + x)/2.0 )
    else             !stable
        stab_h = -5.0 * zL
    end if

  end function PsiH

  !--------------------------------------------------------------------
  function PsiM(zL) result (stab_m)
   !   Out:
   !   PsiM = integral flux-gradient stability function for momentum 
   !   Ref: Garratt, 1994, pp52-54

    real, intent(in) ::  zL    ! = surface layer stability parameter, (z-d)/L 
                               ! notation must be preserved         
    real :: stab_m
    real  :: x
 
    if( zL < 0) then !unstable
       x    = sqrt(sqrt(1.0 - 16.0*zL))
       stab_m = log( 0.125*(1.0+x)*(1.0+x)*(1.0+x*x) ) +  PI/2.0 - 2.0*atan(x)
    else             !stable
       stab_m = -5.0 * zL
    end if

  end function PsiM

!--------------------------------------------------------------------
  function Wind_at_h(u_ref, z_ref, zh, d, z0, Linv) result (u_zh)
!...
!   Ref: Garratt, 1994, 
!   In:
    real, intent(in) ::   u_ref  ! windspeed at z_ref
    real, intent(in) ::   z_ref  ! centre of call, ca. 45m (m)
    real, intent(in) ::   zh     ! height required (m)
    real, intent(in) ::   d      ! displacement height (m)
    real, intent(in) ::   z0     ! roughness height (m)
    real, intent(in) ::   Linv   ! inverse of the Obukhov length (1/m)

!   Out:
    real :: u_zh    ! =   wind-speed at height h (m/s)

     u_zh = u_ref *  &
          ( log((zh-d)/z0)  -PsiM((zh-d)*Linv) + PsiM(z0*Linv)  )/ &
          ( log((z_ref-d)/z0) -PsiM((z_ref-d)*Linv) + PsiM(z0*Linv))

    !NB - COULD USE INSTEAD: Ra = log(z2/z1) - PsiM(z2*Linv) + PsiM(z1*Linv)
    ! Or could optimise with explicit PsiM, etc.

  end function Wind_at_h

end module Micromet_ml
