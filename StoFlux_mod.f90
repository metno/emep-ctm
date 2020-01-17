! <StoFlux_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.34>
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
!> StoFlux_mod.f90 - A component of the EMEP MSC-W Chemical transport Model

module StoFlux_mod
  use CheckStop_mod
  use Config_module,         only: NLANDUSEMAX, dt_advec, MasterProc
  use Debug_module,          only: DEBUG, DebugCell ! -> DEBUG%GRIDVALUES
  use DO3SE_mod,             only: do3se, nSumVPD, SumVPD_LC
  use GasParticleCoeffs_mod, only: DDspec
  use Io_Progs_mod,          only: current_date, datewrite
  use LandDefs_mod,          only: LandType, STUBBLE, iLC_grass
  use LocalVariables_mod,    only: L, Grid
  use MicroMet_mod,          only: AerRes, Wind_at_h
  use Par_mod,               only: LIMAX, LJMAX
  use PhysicalConstants_mod, only: AVOG, KARMAN
  use SmallUtils_mod,        only: find_index
  use SubMet_mod,            only: Sub
  implicit none
  private

  public :: Setup_StoFlux
  public :: Calc_StoFlux

  logical, public,  dimension(NLANDUSEMAX), save :: luflux_wanted

  real, private :: u_hveg ! Values at canopy top, for fluxes
  real, public,  dimension(NLANDUSEMAX), save :: &
      lai_flux,      & ! Fluxes to total LAI
      unit_flux        ! Fluxes per m2 of leaf area (flag-leaf)

  real,   public,save, allocatable,dimension(:,:,:) :: &
       SumVPD ,   &   ! For critical VPD calcs, reset each day  
       old_gsun       !
  integer, private, save, dimension(NLANDUSEMAX) :: indSumVPD

  real, private, save :: gext_leaf = 1.0/2500.0
  real, private :: rc_leaf, rb_leaf
  integer, private, save :: idepO3 



contains
 !----------------------------------------------------------------------------
  subroutine Setup_StoFlux(i,j)

    integer, intent(in) :: i,j !  coordinates
    logical, save :: my_first_call = .true.
    integer ::  istat, iL
    character(len=*), parameter :: dtxt='SetupFST:'

     if ( my_first_call ) then
       idepO3 = find_index('O3',DDspec(:)%name)
       allocate(SumVPD(LIMAX,LJMAX,nSumVPD))
       allocate(old_gsun(LIMAX,LJMAX,nSumVPD))
       SumVPD   = 0.0
       old_gsun = 1.0e99
       do iL = 1, NLANDUSEMAX
         if ( do3se(iL)%VPDcrit > 0.0  ) then
           indSumVPD(iL) = find_index( iL, SumVPD_LC )
         end if
       end do
       my_first_call = .false.
     end if

     Sub(:)%FstO3 = 0.0

    ! resets whole grid at local night
    if ( Grid%Idirect < 1.0e-6 .and. old_gsun(i,j,1)<1.0e6  ) then
    !if ( Grid%Idirect < 1.0e-6 ) then
        SumVPD(i,j,:)        = 0.0    ! For Critical VPD stuff, wheat
        old_gsun(i,j,:)      = 1.0e99 ! "     "
    end if
    if ( DebugCell ) call datewrite(dtxt, [Grid%Idirect, old_gsun(i,j,1) ])

  end subroutine Setup_StoFlux

 !----------------------------------------------------------------------------

  subroutine Calc_StoFlux(nLC,iL_used,debug_flag)
    integer, intent(in) :: nLC
    integer, dimension(nLC), intent(in) :: iL_used
    logical, intent(in) :: debug_flag
    logical :: dbg, dbghh
    character(len=*), parameter :: dtxt='CalcFST:'

    real :: tmp_gsun
    integer :: i,j, iiL, iL, ivpdLC
    ! Evapotranspiration (not currently used!) needs:
    real :: gv, gvcms  ! conductace for water vapour, mmole/ms/s and cm/s

    i = Grid%i
    j = Grid%j
    dbg = DEBUG%STOFLUX .and. DebugCell
    dbghh = dbg .and. current_date%seconds == 0

    LC_LOOP: do iiL = 1, nLC
      iL = iL_used(iiL) 
      L = Sub(iL)

      Sub(iL)%FstO3       = 0.0

       ! take care of  temperate crops, outside growing season
      if ( L%hveg < 1.1 * L%z0 ) then 

        Sub(iL)%cano3_ppb   = 0.0  !! Can't do better?
        Sub(iL)%EvapTransp  = 0.0   ! evapo-transpiration. Not used anyway..
        if ( dbghh ) call datewrite(dtxt//" hveg < z0 ", iL, [ L%hveg, L%z0 ] )

      else !=======================

       ! The fraction going to the stomata = g_sto/g_sur = g_sto * R_sur.
       ! Vg*nmole_o3 is the total deposition flux of ozone, but
       ! we calculate the actual flux later (once we know DepLoss(O3)).
       ! For now we just calculate the g_sto*R_sur bit:
       ! (Caution - g_sto is for O3 only)

       if( L%g_sun > 0.0 ) then ! J20

          !Could be coded faster with Ra....

         u_hveg = Wind_at_h( Grid%u_ref, Grid%z_ref, L%hveg,L%d,L%z0,L%invL )

         rc_leaf = 1.0/(L%g_sto+ gext_leaf)

         !McNaughton + van den Hurk:
 
         if ( do3se(iL)%Lw > 0 )  then
           rb_leaf = 1.3 * 150.0 * sqrt(do3se(iL)%Lw/u_hveg)
         else ! default (CAREFUL!)
           rb_leaf = 1.3 * 150.0 * sqrt(0.05/u_hveg)
         end if

        ! VPD limitation for wheat

         if ( do3se(iL)%VPDcrit > 0.0  ) then
           ivpdLC = indSumVPD(iL)
           SumVPD(i,j,ivpdLC) = SumVPD(i,j,ivpdLC) + L%vpd*dt_advec/3600.0
           tmp_gsun = L%g_sun

           if ( SumVPD(i,j,ivpdLC) > do3se(iL)%VPDcrit ) then
               L%g_sun = min( L%g_sun, old_gsun(i,j,ivpdLC) )
              if( dbghh .and. abs(tmp_gsun-L%g_sun)>1.0e-6 ) call datewrite(dtxt//"SUMVPD", iL, &
               [ do3se(iL)%VPDcrit, L%rh, L%t2C,  L%vpd, SumVPD(i,j,ivpdLC), &
                 tmp_gsun, L%g_sun, old_gsun(i,j,ivpdLC) ] ) 
            end if

           old_gsun(i,j,ivpdLC) = L%g_sun
         end if

       ! Flux in nmole/m2/s:

         Sub(iL)%FstO3 = L%cano3_nmole * rc_leaf/(rb_leaf+rc_leaf) * L%g_sun 

       end if ! g_sun > 0.0
       !if( dbg ) call datewrite(dtxt//" O3 ", &
       !  [ iiL, nLC, iL ], [ L%hveg, L%g_sun, L%cano3_nmole, L%cano3_ppb ] )

! ======   CLOVER  ===========================================================
      ! For Clover we have a very special procedure, using O3 from grassland
      ! to scale the fluxes. As grassland is entered in Inputs.Landuse before
      ! clover we can assume that Fst and O3 for grassland are available and 
      ! correct
        if( LandType(iL)%is_clover) then

           if ( L%g_sun>0) Sub(iL)%FstO3 = &
              Sub(iLC_grass)%cano3_ppb/Sub(iL)%cano3_ppb * Sub(iL)%FstO3

           if (  dbghh ) call datewrite(dtxt//"CLOVER ", iL, &
                 [ Sub(iL)%FstO3, Sub(iLC_grass)%cano3_ppb, &
                    Sub(iLC_grass)%cano3_ppb/Sub(iL)%cano3_ppb ] )
        end if ! clover
! ======   CLOVER  =========================================================

! ======   Evapo-transpiration =============================================
!
       ! E = g * (Cleaf - Cair) = g * D(kPa)/101(KPa), cf. Cambell+Norman,p82
       ! If E in mole/m2/s

       ! and use 1/1.51 instead of 1/1.6. Checl
       ! Lisa had R in s/mm, hence 1000.0
       !
       !gv = 1/(Rb+Rc) = 1(Rb+1/gsto)
       !  Rb(O3) = 2.0 * Rb_cor(WES_O3)/(KARMAN*ustar) 
       !
       ! Step 1: Get gv in m/s units:

        gv = 0.0
        gvcms = 0.0
        Sub(iL)%EvapTransp  = 0.0   ! evapo-transpiration. Not used anyway..

        if( L%g_sto > 1.0e-10 ) then
          gv = 1.0/ (2.0 * DDspec(idepO3)%Rb_cor /(KARMAN*L%ustar) &
             + 1.0/ ( L%g_sto * L%LAI ) )
          gvcms = gv

       ! Step 2: Convert to mole/m2/s and for H2O
       ! mmol2sm = 8.3144e-8*L%t2  for O3, plus factor 1.6 for H2O
       ! * 0.001 -> mole/m2/s
       ! ms2molm2s = 1.6*1.0e-3/(8.3144e-8*L%t2)
       
          gv = gv * 1.6e-3/(8.3144e-8*L%t2)
         
       ! Step 3: 
       ! Mass flux density is E x 0.018 kg/mole -> kg/m2/s
       ! 1 kg/m2 = 1mm
       ! ie ms2kgm2s  = 1.6*1.0e-3/(8.3144e-8*L%t2)
       ! ie ms2mm     = 1.6*1.0e-3/(8.3144e-8*L%t2)

          Sub(iL)%EvapTransp = 0.018 *  L%vpd/101.0  * gv
        end if

       !!L%g_sto * L%LAI  * 1.6e-3/(8.3144e-8*L%t2)   ! Evapo-transpiration
       !                 (rb_leaf/1.6 + 0.0224*(L%t2/273.0) 


        if ( dbghh ) call datewrite(dtxt//"VALS ", iL, &
            [ L%hveg,L%LAI,L%g_sto,L%g_sun,&
              u_hveg, Sub(iL)%cano3_ppb, Sub(iL)%FstO3 ] ) ! skip gvcms

      end if

    end do LC_LOOP
  end subroutine Calc_StoFlux

!..............................................................................
end module StoFlux_mod
