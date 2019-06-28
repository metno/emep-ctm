! <CellMet_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module CellMet_mod
!=============================================================================
!+
! Description
!  Module for setting some near-surface meteorology params
!  **  calls SubMet_mod **
!  for calculating sub-grid meteorology for each land-use.
!=============================================================================

use CheckStop_mod,     only: CheckStop
use DerivedFields_mod, only: d_2d, f_2d
use GridValues_mod,    only: dA,dB, glat, glon
use Landuse_mod,       only: LandCover, ice_landcover ! Provides SGS,hveg,LAI,...
use Landuse_mod,       only: mainly_sea
use LandDefs_mod,      only: LandType
use LocalVariables_mod,only: Grid, ResetSub
use MicroMet_mod,      only: PsiH, PsiM, AerRes       ! functions
use MetFields_mod,     only: ps, u_ref, cc3dmax, sdepth, surface_precip, &
                            ice_nwp,fh, fl, z_mid, z_bnd, q, roa, rh2m, sst, &
                            rho_surf, th, pzpbl, t2_nwp, ustar_nwp, zen,&
                            coszen, Idirect, Idiffuse
use Config_module,    only: KMAX_MID, KMAX_BND, PT, USE_ZREF, IOU_INST
use PhysicalConstants_mod, only: PI, CP, GRAV, KARMAN
use SoilWater_mod,     only: fSW
use SubMet_mod,        only: Get_SubMet, Sub

implicit none
private

!Subroutines

public :: Get_CellMet   ! sets Grid-average (e.g. NWP) near-surface met, and
                        ! calls Get_Submet routines

integer, save, public :: z0_out_ix = -1, invL_out_ix = -1

contains
!=======================================================================

subroutine Get_CellMet(i,j,debug_flag)
  integer, intent(in) :: i,j
  logical, intent(in) :: debug_flag  ! DEBUG%RUNCHEM + wanted i,j
  integer :: lu, ilu, nlu
  real :: land_frac
  character(len=*), parameter :: dtxt='GetCell:'
!---------------------------------------------------------------

     if(z0_out_ix>0) d_2d(z0_out_ix,i,j,IOU_INST) = 0.0

! We assume that the area of grid which is wet is proportional to
! cloud-cover. To avoid some compiler/numerical issues when
! prec almost equal to zero, we allow a small build-up phase, with
! linear increase from wetarea=0 to wetarea = cc3dmax for values of
! prec between 1.0e-8 (near-zero!) to 0.01.

  if ( surface_precip(i,j) > 1.0d-2 ) then
    Grid%is_wet = .true.
    Grid%wetarea = cc3dmax(i,j,KMAX_MID)
  elseif ( surface_precip(i,j) > 1.0d-8 ) then
    Grid%is_wet = .true.
    Grid%wetarea =  surface_precip(i,j)/1.0d-2 * cc3dmax(i,j,KMAX_MID)
  else
    Grid%is_wet = .false.
    Grid%wetarea = 0.0
  end if


  Grid%i        = i
  Grid%j        = j
  Grid%latitude  = glat(i,j) !SPOD tests
  Grid%longitude = glon(i,j)

  Grid%psurf    = ps(i,j,1)             ! Surface pressure, Pa
  Grid%z_mid    = z_mid(i,j,KMAX_MID)   ! NB: Approx, updated every 3h
  Grid%rh2m = rh2m(i,j,1) !NWP value !now used in BiDir
  Grid%sst  = sst(i,j,1)  !NWP value !to be used in BiDir
  Grid%is_frozen = Grid%t2 < 271.15 ! BiDir extra


  ! Have option to use a different reference ht:
  if ( USE_ZREF ) then
    Grid%z_ref    = &
    min( 0.1*pzpbl(i,j),  z_mid(i,j,KMAX_MID) )   ! within or top of SL
  else
    Grid%z_ref    = z_mid(i,j,KMAX_MID)    ! within or top of SL
  end if

 ! The biggest trees in the new CLM ssytem are 35m high, giving displacement
 ! hts of 24.5m. We ensure that z_ref -d > z0

  Grid%z_ref    = max(30.0, Grid%z_ref ) ! for trees d=14m,

  ! More exact for thickness of bottom layer, since used for emissions
  ! from dp = dA+dB*Ps (eta coordinates)
  Grid%DeltaZ &!  = z_bnd(i,j,KMAX_BND-1) ! NB! Approx,updated every 3h
            = (dA(KMAX_MID)+dB(KMAX_MID)*ps(i,j,1) )/(GRAV*roa(i,j,KMAX_MID,1))
  Grid%u_ref    = u_ref(i,j)
  Grid%qw_ref   =  q(i,j,KMAX_MID,1)   ! specific humidity
  Grid%rho_ref  = roa(i,j,KMAX_MID,1)
  Grid%zen = zen(i,j)
  Grid%coszen = coszen(i,j)
  Grid%izen = max( 1, int ( Grid%zen + 0.5 ) )! 1 avoids zero in indices.
  Grid%Idirect  =  Idirect(i,j)
  Grid%Idiffuse =  Idiffuse(i,j)

  !**  prefer micromet signs and terminology here:
  Grid%Hd    = -fh(i,j,1)       ! Heat flux, *away from* surface
if( debug_flag ) write(*,"(a,3es12.3)") 'CellHd', Grid%Hd, &
   maxval(fh(:,:,1)), minval(fh(:,:,1))
  Grid%LE    = -fl(i,j,1)       ! Heat flux, *away from* surface
  Grid%ustar = ustar_nwp(i,j)   !  u*
  Grid%t2    = t2_nwp(i,j,1)    ! t2 , K
  Grid%t2C   = Grid%t2 - 273.15 ! deg C
  Grid%theta_ref = th(i,j,KMAX_MID,1)
  Grid%rh2m  = rh2m(i,j,1)      !
  Grid%rho_s = rho_surf(i,j)    ! Should replace Met_mod calc. in future

  Grid%is_mainlysea = mainly_sea(i,j)
  Grid%is_allsea = .true. !set to false below if land found
  
  Grid%sdepth    = sdepth(i,j,1)
  Grid%ice_nwp   = max( ice_nwp(i,j,1), ice_landcover(i,j) )
  Grid%snowice   = ( Grid%sdepth  > 1.0e-10 .or. Grid%ice_nwp > 1.0e-10 )

  Grid%fSW       = fSW(i,j)

  ! we limit u* to a physically plausible value
  ! to prevent numerical problems
  Grid%ustar = max( Grid%ustar, 0.1 )

  !NB: invL_nwp is already defined, with similar defintion 
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
  Sub(:)          = ResetSub
  Sub(:)%coverage = 0.0
  Sub(:)%LAI      = 0.0
  Sub(:)%SAI      = 0.0
  Sub(:)%hveg     = 0.0

  
  LULOOP: do ilu= 1, nlu
    lu = LandCover(i,j)%codes(ilu)

    if((.not. LandType(lu)%is_water) .and. LandCover(i,j)%fraction(ilu)>0.0) Grid%is_allsea = .false. 

    Sub(lu)%coverage = LandCover(i,j)%fraction(ilu)
    Sub(lu)%LAI      = LandCover(i,j)%LAI(ilu)
    Sub(lu)%SAI      = LandCover(i,j)%SAI(ilu)
    Sub(lu)%hveg     = LandCover(i,j)%hveg(ilu)

    !=======================
    call Get_SubMet(lu, debug_flag )
        
    Sub(lu)%SWP = 0.0  ! Not yet implemented
    
    !=======================
  end do LULOOP


  !if requested, make weighted average for output
  if(z0_out_ix>0 .or. invL_out_ix>0)then
     !we need a separate loop, beacause Grid%is_allsea is set in the first one
     land_frac=0.0
     do ilu= 1, nlu        
        lu = LandCover(i,j)%codes(ilu)
       if(Grid%is_allsea)then
           if(z0_out_ix>0) d_2d(z0_out_ix,i,j,IOU_INST) = log(Sub(lu)%z0)
           if(invL_out_ix>0) d_2d(invL_out_ix,i,j,IOU_INST) = Sub(lu)%invL 
        else if(.not.LandType(lu)%is_water) then
           land_frac = land_frac+LandCover(i,j)%fraction(ilu)
           if(z0_out_ix>0) d_2d(z0_out_ix,i,j,IOU_INST) = &
                d_2d(z0_out_ix,i,j,IOU_INST) + log(Sub(lu)%z0)*LandCover(i,j)%fraction(ilu)
           if(invL_out_ix>0)  d_2d(invL_out_ix,i,j,IOU_INST) = &
                d_2d(invL_out_ix,i,j,IOU_INST) + Sub(lu)%invL*LandCover(i,j)%fraction(ilu)
        endif
     enddo
     if(.not.Grid%is_allsea)then
        !renormalize to land
        if(land_frac >= 1.E-6)then
           if(z0_out_ix>0) d_2d(z0_out_ix,i,j,IOU_INST) = &
                d_2d(z0_out_ix,i,j,IOU_INST)/land_frac
           if(invL_out_ix>0)  d_2d(invL_out_ix,i,j,IOU_INST) = &
                d_2d(invL_out_ix,i,j,IOU_INST)/land_frac
        else
           write(*,*) dtxt//'WARNING: found grid with no sea and no land',&
                     i,j,nlu,land_frac
           do ilu= 1, nlu  
              lu = LandCover(i,j)%codes(ilu)
              write(*,*)dtxt, lu,LandCover(i,j)%fraction(ilu),LandType(lu)%is_water
          enddo
       endif
    endif
     
  endif
  

end subroutine Get_CellMet
!=======================================================================

end module CellMet_mod
