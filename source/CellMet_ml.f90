! <CellMet_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module CellMet_ml
!=============================================================================
!+
! Description
!  Module for setting some near-surface meteorology params
!  **  calls SubMet_ml ** 
!  for calculating sub-grid meteorology for each land-use. 
!=============================================================================


use CheckStop_ml, only : CheckStop
use GridValues_ml, only :  sigma_bnd,dA,dB
use Landuse_ml, only : LandCover, ice_landcover    ! Provides SGS, hveg, LAI ....
use Landuse_ml, only : mainly_sea
use LocalVariables_ml, only: Grid, ResetSub
use MicroMet_ml, only :  PsiH, PsiM, AerRes    !functions
use MetFields_ml, only: ps, u_ref
use MetFields_ml, only: cc3dmax, sdepth,ice_nwp, surface_precip, &
   fh,fl,z_mid, z_bnd, q, roa, rh2m, rho_surf, th, pzpbl, t2_nwp, ustar_nwp,&
    zen, coszen, Idirect, Idiffuse
use ModelConstants_ml,    only : KMAX_MID, KMAX_BND, PT, USE_ZREF
use PhysicalConstants_ml, only : PI, RGAS_KG, CP, GRAV, KARMAN, CHARNOCK, T0
use SoilWater_ml, only : fSW
use SubMet_ml, only : Get_SubMet, Sub
use TimeDate_ml, only: current_date

implicit none
private

!Subroutines

public :: Get_CellMet   ! sets Grid-average (e.g. NWP) near-surface met, and
                        ! calls Get_Submet routines

logical, private, parameter ::  MY_DEBUG = .false.

contains
!=======================================================================

  subroutine Get_CellMet(i,j,debug_flag)
    integer, intent(in) :: i,j
   logical, intent(in) :: debug_flag   ! set true for wanted grid square
    integer :: lu, ilu, nlu

!---------------------------------------------------------------


  !   We assume that the area of grid which is wet is proportional to 
  !   cloud-cover. To avoid some compiler/numerical issues when
  !   prec almost equal to zero, we allow a small build-up phase, with
  !   linear increase from wetarea=0 to wetarea = cc3dmax for values of
  !   prec between 1.0e-8 (near-zero!) to 0.01.

     if ( surface_precip(i,j) > 1.0d-2 ) then
          Grid%is_wet = .true.
          Grid%wetarea = cc3dmax(i,j,KMAX_MID) 
     else if ( surface_precip(i,j) > 1.0d-8 ) then
          Grid%is_wet = .true.
          Grid%wetarea =  surface_precip(i,j)/1.0d-2 * cc3dmax(i,j,KMAX_MID) 
     else 
          Grid%is_wet = .false.
          Grid%wetarea = 0.0
     end if


     Grid%i        = i
     Grid%j        = j
     Grid%psurf    = ps(i,j,1)    ! Surface pressure, Pa
     Grid%z_mid    = z_mid(i,j,KMAX_MID)   ! NB! Approx, updated every 3h

     ! Have option to use a different reference ht:
     if ( USE_ZREF ) then 
        Grid%z_ref    = &
           min( 0.1*pzpbl(i,j),  z_mid(i,j,KMAX_MID) )   ! within or top of SL
     else
        Grid%z_ref    = z_mid(i,j,KMAX_MID)    ! within or top of SL
     end if

! More exact for thickness of bottom layer, since used for emissions
! from  dp = g. rho . dz and d sigma = dp/pstar
! we get dz = d sigma . pstar/(g.rho)
!     Grid%DeltaZ  &!  = z_bnd(i,j,KMAX_BND-1) ! NB! Approx,updated every 3h
!                = (1.0 - sigma_bnd(20) ) * (ps(i,j,1)-PT) /(GRAV*roa(i,j,20,1))
!Eta coordinates:
!dp = dA+dB*Ps
     Grid%DeltaZ  &!  = z_bnd(i,j,KMAX_BND-1) ! NB! Approx,updated every 3h
                = (dA(KMAX_MID)+dB(KMAX_MID)*ps(i,j,1) )/(GRAV*roa(i,j,KMAX_MID,1))

     Grid%u_ref    = u_ref(i,j)
     Grid%qw_ref    =  q(i,j,KMAX_MID,1)   ! specific humidity
     Grid%rho_ref  = roa(i,j,KMAX_MID,1)
     Grid%zen = zen(i,j)
     Grid%coszen = coszen(i,j)
     Grid%izen = max( 1, int ( Grid%zen + 0.5 ) )! 1 avoids zero in indices.
     Grid%Idirect  =  Idirect(i,j)
     Grid%Idiffuse =  Idiffuse(i,j)

     !**  prefer micromet signs and terminology here:

     Grid%Hd    = -fh(i,j,1)       ! Heat flux, *away from* surface
     Grid%LE    = -fl(i,j,1)       ! Heat flux, *away from* surface
     Grid%ustar = ustar_nwp(i,j)   !  u*
     Grid%t2    = t2_nwp(i,j,1)    ! t2 , K
     Grid%t2C   = Grid%t2 - 273.15 ! deg C
     Grid%theta_ref = th(i,j,KMAX_MID,1) 
     Grid%rh2m  = rh2m(i,j,1)      ! 
     Grid%rho_s = rho_surf(i,j)    ! Should replace Met_ml calc. in future

     Grid%is_mainlysea = mainly_sea(i,j)
     Grid%is_allsea = ( mainly_sea(i,j) .and. LandCover(i,j)%ncodes == 1)
     Grid%sdepth    = sdepth(i,j,1)
     Grid%ice_nwp   = max( ice_nwp(i,j,1), ice_landcover(i,j) ) 
     Grid%snowice   = ( Grid%sdepth  > 1.0e-10 .or. Grid%ice_nwp > 1.0e-10 )

     Grid%fSW       = fSW(i,j)

    ! we limit u* to a physically plausible value
    ! to prevent numerical problems

     Grid%ustar = max( Grid%ustar, 0.1 )

     Grid%invL  = -1* KARMAN * GRAV * Grid%Hd & ! -Grid%Hd disliked by gfortran
            / (CP*Grid%rho_s * Grid%ustar*Grid%ustar*Grid%ustar * Grid%t2 )

    !.. we limit the range of 1/L to prevent numerical and printout problems
    !.. and because we don't trust HIRLAM or other NWPs enough.
    !   This range is very wide anyway.

    ! Grid%invL  = max( -1.0, Grid%invL ) !! limit very unstable
    ! Grid%invL  = min(  1.0, Grid%invL ) !! limit very stable


    ! wstar for particle deposition, based on Wesely
    if(Grid%Hd >  0.0 ) then          ! unstable stratification
        Grid%wstar = ( GRAV * pzpbl(i,j) * Grid%Hd /      &
        (Grid%rho_ref * CP * th(i,j,KMAX_MID,1))) ** (1./3.)
    else
         Grid%wstar = 0.
    end if

  nlu = LandCover(i,j)%ncodes
  ! Added for safety
  Sub(:)          = ResetSub !
  Sub(:)%coverage = 0.0
  Sub(:)%LAI      = 0.0
  Sub(:)%SAI      = 0.0
  Sub(:)%hveg     = 0.0

    LULOOP: do ilu= 1, nlu
        lu      = LandCover(i,j)%codes(ilu)

        Sub(lu)%coverage = LandCover(i,j)%fraction(ilu)
        Sub(lu)%LAI      = LandCover(i,j)%LAI(ilu)
        Sub(lu)%SAI      = LandCover(i,j)%SAI(ilu)
        Sub(lu)%hveg     = LandCover(i,j)%hveg(ilu)

       !=======================

        call Get_SubMet(lu, debug_flag )

        Sub(lu)%SWP   =  0.0  ! Not yet implemented
       !=======================
        end do LULOOP

  end subroutine Get_CellMet
  !=======================================================================

end module CellMet_ml
