! <BoundaryConditions_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module BoundaryConditions_mod
! -----------------------------------------------------------------------
! This module is the main driver module for defining and setting
! initial & boundary conditions (ICs & BCs) for the chemical species.
!
! The main code calls up these routines with just:
!   call BoundaryConditions(month)  !once per month, after say newmonth
!
! On first call, this routine runs some intialisation routines in related
! modules, reads the global data, and sets full 3-D concentration fields
! of the advected and background concentration fields (xn_adv, xn_bgn).
! On subsequent calls (first_call=.false.), the routine reads new
! global input data, and sets the concentrations at the top level and
! lateral boundaries for advected species. For background species it reset
! full 3-D concentration fields.

! The main related IC/BC modules/files are:
!
!   1. GetBICData, reads or define BIC from tables and functions (replaces older GlobalBC_mod)
!
!   2. CM_BoundaryConditions.inc - assigns mappings, telling which
!      Unified EMEP model species the global model fdata are
!      assigned to (bc2xn_adv, bc2xn_bgn arrays).
!
!   3. NB: Nest_mod.f90 - In nested runs (such as forecast runs), the ICs & BCs
!      are reseted by readxn (Nest_mod), superseeding the IC/BCs in this module.
! -----------------------------------------------------------------------
! 
!    The routines make use of a "feature" of the model: that the concentration
!    (xn) values  along boundaries are not changed due to advection or chemistry.
!    Thus, bc values only need to be set once per month, firstly for the whole 3-D
!    domain, then monthly for the sides and top.
!    Background species must be reset in 3-D each month.
! -----------------------------------------------------------------------

use CheckStop_mod,      only: CheckStop
use Chemfields_mod,     only: xn_adv, xn_bgn, NSPEC_BGN  ! emep model concs.
use ChemDims_mod,       only: NSPEC_SHL,NSPEC_ADV
use ChemSpecs_mod                ! provide species names, IXADV_*
use Config_module, only: KMAX_MID  &  ! Number of levels in vertical
                     ,iyr_trend &  ! Used for e.g. future scenarios
                     ! Two options for CH4. BGND_CH4 has priority.
                     ,BGND_CH4  &  ! If positive, replaces defaults 
                     ,fileName_CH4_ibcs & ! If present, replaces uses iyr_trend
                     ,USES, MasterProc, PPB, Pref, LoganO3File, DustFile
use Debug_module,    only: DEBUG   ! -> DEBUG%BCS
use Functions_mod,   only: StandardAtmos_kPa_2_km ! for use in Hz scaling
use GridValues_mod,     only: glon, glat   & ! full domain lat, long
                            ,debug_proc, debug_li, debug_lj & ! debugging
                            ,i_fdom, j_fdom,A_mid,B_mid  !
use Io_mod,             only: open_file, ios, IO_TMP
use Io_Progs_mod,       only: datewrite, PrintLog
use Landuse_mod,        only: mainly_sea
use LocalVariables_mod, only: Grid
use MetFields_mod,      only: roa
use MPI_Groups_mod,     only: MPI_DOUBLE_PRECISION, MPI_SUM,MPI_INTEGER, &
                             MPI_COMM_CALC, IERROR
use NetCDF_mod,         only: ReadField_CDF,vertical_interpolate
use Par_mod,          only: &
   LIMAX, LJMAX, limax, ljmax, me &
  ,neighbor, NORTH, SOUTH, EAST, WEST   &  ! domain neighbours
  ,NOPROC&
  ,IRUNBEG,JRUNBEG,li1,li0,lj0,lj1
use PhysicalConstants_mod, only: PI, ATWAIR
use SmallUtils_mod, only : find_index 
use TimeDate_mod,    only: daynumber
use TimeDate_ExtraUtil_mod,only: date2string

implicit none
private

integer, public, parameter :: &
   IBC_O3       =  1   &
  ,IBC_NO       =  2   &
  ,IBC_NO2      =  3   &
  ,IBC_PAN      =  4   &
  ,IBC_HNO3     =  5   &  ! used for nitrate too
  ,IBC_SO2      =  6   &
  ,IBC_SO4      =  7   &
  ,IBC_CO       =  8   &
  ,IBC_C2H6     =  9   &
  ,IBC_C4H10    = 10   &
  ,IBC_HCHO     = 11   &
  ,IBC_CH3CHO   = 12   &
  ,IBC_H2O2     = 13   &
  ,IBC_NH4_f    = 14   &
  ,IBC_NO3_f    = 15   &
  ,IBC_NO3_c    = 16   &
  ,IBC_SeaSalt_f= 17   &
  ,IBC_SeaSalt_c= 18   &
  ,IBC_SeaSalt_g= 19   &
  ,IBC_Dust_f   = 20   &      ! Dust
  ,IBC_Dust_c   = 21   &      ! Dust
  ,NGLOB_BC     = IBC_Dust_c  ! Totan no. species setup in this module

! -- subroutines in this module:
public  :: BoundaryConditions         ! call every month
private :: GetBICData, &
           Set_bcmap,               & ! sets xn2adv_changed, etc.
           My_bcmap                   ! sets bc2xn_adv, bc2xn_bc, and  misc_bc

!/- Allow different behaviour on 1st call - full 3-D asimilation done
logical, private, save :: first_call  = .true.

 !/ - for debugging
logical, private, parameter :: DEBUG_MYBC = .false.

! Set indices
! -----------------------------------------------------------------------
! For species which have constant mixing ratios:
integer, public, parameter :: &
  NMISC_BC = 2,               &
  IBC_H2  = NGLOB_BC + 1,     &
  IBC_CH4 = NGLOB_BC + 2,     &
  NTOT_BC  = NGLOB_BC + NMISC_BC

! misc_bc specifies concentrations of these species:
real, public, allocatable,save, dimension(:,:) :: misc_bc

! Define mapping arrays
! -----------------------------------------------------------------------
! The mapping is done through the arrays bc2xn_adv and bc2xn_bgn, such that
! the emep species are given along the x-dimension and the bc species along
! the y.  e.g., the statement
!   bc2xn_adv(IBC_NOX,IXADV_NO2) = 0.55
! would assign the BC concentration of NOX  to the EMEP model concentration
! of NO2 after multiplication with a factor 0.55.
! (The CTM2 concentration of NOx used as BC after a multilication
!  with a factor 0.55.)
! -----------------------------------------------------------------------

real, public, save :: bc2xn_adv(NTOT_BC,NSPEC_ADV), & ! see above
                      bc2xn_bgn(NTOT_BC,NSPEC_BGN)

! Arrays for mapping from global bc to emep xn concentrations:
! ---------------------------------------------------------------------------
integer, private,save, dimension(NTOT_BC) :: &
  bc_used,bc_used_adv,bc_used_bgn ! set to 1 if bc used

integer, private, save  ::  &
  num_adv_changed,          &  ! Num. adv. species that have bc's
  num_used_adv,             &  ! Max times a conc. from e.g CTM2 is used as bc
  num_bgn_changed,num_used_bgn, &
  num_changed                   ! sum of adv and bgn

!In general "changed" means "bc used for this species"
logical, private, save ::     &
  xn_adv_changed(NSPEC_ADV),  & ! true if emep xn_adv changed by bcs
  xn_bgn_changed(NSPEC_BGN)     ! true if emep xn_bgn changed by bcs

integer, private, save ::     &
  spc_adv2changed(NSPEC_ADV), & ! index of advected specie is converted to index
  spc_bgn2changed(NSPEC_BGN)    ! in the row of advected species that have bc

integer, allocatable, dimension(:),save ::  &
  spc_changed2adv, &           ! index of adv. specie that have bc is converted
  spc_changed2bgn              ! to index in the row of advected species

integer, allocatable, dimension(:,:),save :: &
  spc_used_adv, &          ! 1.dimension(ibc) runs through bc adv.species
  spc_used_bgn             ! 2.dim.(i) through those who get same conc. from
                           ! bc (eg i=1,2 for ibc=HNO3 when HNO3 from CTM2
                           ! is used as bc both for HNO3 and SO4) spc_used_adv
                           ! gives the index in the row of advected species
real, allocatable,dimension(:,:,:),save   :: O3_logan,O3_logan_emep
real, allocatable,dimension(:,:,:),save   :: Dust_3D, Dust_3D_emep

contains

  subroutine BoundaryConditions(year,month)
    ! ---------------------------------------------------------------------------
    ! Read in monthly-average global mixing ratios, and if found, collect the
    ! data in  bc_adv, bc_bgn arrays for later interpolations
    ! NOTES
    ! 1.- If mixing ratio by mass the scale by molcular weight)
    ! 2.- So far no scaling is done, but this could be done
    !     in Set_bcmap with atomic weights
    ! 3.- On the first call, we also run the setup-subroutines
    ! 4.- Year is now obtained from the iyr_trend set in run.pl.
    !     This allows, e.g. runs with BCs for 2100 and met of 1990.
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: year         ! "meteorology" year
    integer, intent(in) :: month
    integer :: ibc, iem, k, i, j ,n, ntot ! loop variables
    integer :: info                     ! used in rsend
    real    :: bc_fac      ! Set to 1.0, except sea-salt over land = 0.01
    logical :: bc_seaspec  ! if sea-salt species

    integer  :: errcode
    integer, save :: idebug=0, itest=1, i_test=0, j_test=0
    real :: bc_data(LIMAX,LJMAX,KMAX_MID)

    if (first_call) then
       if (DEBUG%BCS) write(*,"(a,I3,1X,a,i5)") &
            "FIRST CALL TO BOUNDARY CONDITIONS, me: ", me,  "TREND YR ", iyr_trend
       allocate(misc_bc(NGLOB_BC+1:NTOT_BC,KMAX_MID))
       call My_bcmap(iyr_trend)      ! assigns bc2xn_adv and bc2xn_bgn mappings
       call Set_bcmap()              ! assigns xn2adv_changed, etc.

       num_changed = num_adv_changed + num_bgn_changed   !u1
       if (DEBUG%BCS) write(*, "((A,I0,1X))")           &
            "BCs: num_adv_changed: ", num_adv_changed,  &
            "BCs: num_bgn_changed: ", num_bgn_changed,  &
            "BCs: num     changed: ", num_changed

       bc_data=0.0

    end if ! first call
    if (DEBUG%BCS) write(*, "((A,I0,1X))")           &
         "CALL TO BOUNDARY CONDITIONS, me:", me, &
         "month ", month, "TREND2 YR ", iyr_trend, "me ", me

    if (num_changed==0) then
       write(*,*) "BCs: No species requested"
       return
    end if

    errcode = 0
    if (DEBUG%BCS.and.debug_proc) then
       do i = 1, limax
          do j = 1, ljmax
             if (i_fdom(i)==DEBUG%IJ(1).and.j_fdom(j)==DEBUG%IJ(2)) then
                i_test = i
                j_test = j
             end if
          end do
       end do
    end if

    if (first_call) then
       idebug = 1
       if (DEBUG%BCS) write(*,*) "RESET 3D BOUNDARY CONDITIONS", me
       do k = 1, KMAX_MID
          do j = 1, ljmax
             do i = 1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             end do
          end do
       end do
    else       
       if (DEBUG%BCS.and.MasterProc) write(*,*) "RESET LATERAL BOUNDARIES"
       do k = 2, KMAX_MID
          do j = lj0, lj1
             !left
             do i = 1, li0-1
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             end do
             !right
             do i = li1+1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             end do
          end do
          !lower
          do j = 1, lj0-1
             do i = 1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             end do
          end do
          !upper
          do j = lj1+1, ljmax
             do i = 1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             end do
          end do
       end do
       !top
       do k = 1, 1
          do j = 1, ljmax
             do i = 1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             end do
          end do
       end do
    end if
    !== BEGIN READ_IN OF GLOBAL DATA

    do ibc = 1, NGLOB_BC
       if(bc_used(ibc) == 0)cycle
 
       call GetBICData(year,month,ibc,bc_used(ibc),bc_data,errcode)

       if (first_call) then

          ! Set 3-D arrays of new BCs
          do n = 1, bc_used_adv(ibc)

             iem = spc_used_adv(ibc,n)
             ntot = iem + NSPEC_SHL 

             ! Sea-salt. 
             !  If SeaSalt isn't called from mk.GenChem, we don't have the
             !  SS_GROUP, so we search for the simple SEASALT name.
             bc_seaspec = .false.
             if ( USES%SEASALT .and. &
                  ( index( species(ntot)%name, "SeaSalt_" ) > 0 ) ) then
                bc_seaspec = .true.
             end if

             if ( debug_proc ) write (*,*) "SEAINDEX", &
                  trim(species(ntot)%name), n, ntot, bc_seaspec,&
                  index( species(ntot)%name, "SeaSalt_")

             do k = 1, KMAX_MID
                do j = 1, ljmax
                   do i = 1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USES%SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                                !                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)
                           bc_data(i,j,k)*bc2xn_adv(ibc,iem)
                   end do ! i
                end do ! j
             end do ! k
          end do !n

          do n = 1,bc_used_bgn(ibc)
             iem = spc_used_bgn(ibc,n)

             !/- Non-advected background species
             do k = 1, KMAX_MID
                do j = 1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                                !                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                           +  bc_data(i,j,k)*bc2xn_bgn(ibc,iem)
                      !                      !        bc_bgn(n,(i_fdom(i)-IRUNBEG+1),(j_fdom(j)-JRUNBEG+1),k)
                   end do ! i
                end do ! j
             end do ! k
          end do
       else

          ! Set LATERAL (edge and top) arrays of new BCs

          idebug = idebug + 1
          do n = 1, bc_used_adv(ibc)
             iem = spc_used_adv(ibc,n)
             ntot = iem + NSPEC_SHL 
             bc_seaspec = .false.
             if ( USES%SEASALT .and. ( index( species(ntot)%name, "SeaSalt_" ) > 0 ) ) then
                bc_seaspec = .true.
             end if

             do k = 2, KMAX_MID
                do j = lj0, lj1
                   !left
                   do i = 1, li0-1
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USES%SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                                !                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)
                           bc_data(i,j,k)*bc2xn_adv(ibc,iem)
                   end do
                   !right
                   do i = li1+1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USES%SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                                !                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)
                           bc_data(i,j,k)*bc2xn_adv(ibc,iem)

                   end do
                end do
                !lower
                do j = 1, lj0-1
                   do i = 1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USES%SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                                !                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)
                           bc_data(i,j,k)*bc2xn_adv(ibc,iem)

                   end do
                end do
                !upper
                do j = lj1+1, ljmax
                   do i = 1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USES%SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                                !                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)
                           bc_data(i,j,k)*bc2xn_adv(ibc,iem)

                   end do
                end do
             end do
             !top
             do k = 1, 1
                do j = 1, ljmax
                   do i = 1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USES%SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                                !                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)
                           bc_data(i,j,k)*bc2xn_adv(ibc,iem)

                   end do
                end do
             end do

          end do !n

          !/- Non-advected background species
          do n = 1,bc_used_bgn(ibc)
             iem = spc_used_bgn(ibc,n)

             do k = 2, KMAX_MID
                do j = lj0, lj1
                   !left
                   do i = 1, li0-1
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                                !                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                           +  bc_data(i,j,k)*bc2xn_bgn(ibc,iem)
                   end do
                   !right
                   do i = li1+1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                                !                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                           +  bc_data(i,j,k)*bc2xn_bgn(ibc,iem)
                   end do
                end do
                !lower
                do j = 1, lj0-1
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                                !                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                           +  bc_data(i,j,k)*bc2xn_bgn(ibc,iem)
                   end do
                end do
                !upper
                do j = lj1+1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                                !                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                           +  bc_data(i,j,k)*bc2xn_bgn(ibc,iem)
                   end do
                end do
             end do
             !top
             do k = 1, 1
                do j = 1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                                !                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                           +  bc_data(i,j,k)*bc2xn_bgn(ibc,iem)
                   end do
                end do
             end do
          end do
       end if
    end do  ! ibc
    if (first_call) then
       !3D misc
       do ibc = NGLOB_BC+1, NTOT_BC
          do n = 1,bc_used_bgn(ibc)
             iem = spc_used_bgn(ibc,n)            
             !/- Non-advected background misc species
             do k = 1, KMAX_MID
                do j = 1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   end do ! i
                end do ! j
             end do ! k
          end do
          do n = 1,bc_used_adv(ibc)
             iem = spc_used_adv(ibc,n)

             !/- Advected misc species
             do k = 1, KMAX_MID
                do j = 1, ljmax
                   do i = 1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                      !                      !        bc_bgn(n,(i_fdom(i)-IRUNBEG+1),(j_fdom(j)-JRUNBEG+1),k)
                   end do ! i
                end do ! j
             end do ! k
          end do!n
       end do!ibc
    else
       !LATERAL misc
       do ibc = NGLOB_BC+1, NTOT_BC
          do n = 1,bc_used_bgn(ibc)
             iem = spc_used_bgn(ibc,n)            
             !/- Non-advected background misc species
             do k = 2, KMAX_MID
                do j = lj0, lj1
                   !left
                   do i = 1, li0-1
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   end do
                   !right
                   do i = li1+1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   end do
                end do
                !lower
                do j = 1, lj0-1
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   end do
                end do
                !upper
                do j = lj1+1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   end do
                end do
             end do
             !top
             do k = 1, 1
                do j = 1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   end do
                end do
             end do
          end do
          !/- Advected misc species
          do n = 1,bc_used_adv(ibc)
             iem = spc_used_adv(ibc,n)
             do k = 2, KMAX_MID
                do j = lj0, lj1
                   !left
                   do i = 1, li0-1
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   end do
                   !right
                   do i = li1+1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   end do
                end do
                !lower
                do j = 1, lj0-1
                   do i = 1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   end do
                end do
                !upper
                do j = lj1+1, ljmax
                   do i = 1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   end do
                end do
             end do
             !top
             do k = 1, 1
                do j = 1, ljmax
                   do i = 1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   end do
                end do
             end do
          end do!n
       end do!ibc
    end if


    if (DEBUG%BCS.and.debug_proc.and.i_test>0) then
       i = i_test
       j = j_test
       print "(a20,3i4,2f8.2)","DEBUG BCS Rorvik", me, i,j,glon(i,j),glat(i,j)
       print "(a20,3i4)","DEBUG BCS Rorvik DIMS",num_adv_changed
       do k = 1, KMAX_MID
          print "(a20,i4,f8.2)","DEBUG O3  Debug-site ", k, &
               xn_adv(IXADV_O3,i_test,j_test,k)/PPB
       end do
    end if ! DEBUG

    if (DEBUG%BCS.and.debug_proc) then
       itest = 1
       print *,"BoundaryConditions: No CALLS TO BOUND Cs", first_call,idebug
       !*** the following uses hard-coded  IXADV_ values for testing.
       !    Remove later **
       info = 1   ! index for ozone in bcs
       print *,"BCs: bc2xn(info,itest) : ", bc2xn_adv(info,itest)


       info = 43   ! index for NO in bcs
       print *,"BCs: NSPECS: BC, ADV, BG, ", NTOT_BC, NSPEC_ADV, NSPEC_BGN
       print *,"BCs: Number  of bc_used: ", sum(bc_used)
       print *,"BCs: limax, ljmax",  limax, ljmax

       if (NSPEC_BGN>0) then
          do k = KMAX_MID, 1, -1
             print "(a23,i3,e14.4)","BCs NO :",k,xn_bgn(itest,i_test,j_test,k)/PPB
          end do
       else
          print "(a)","No SET BACKGROUND BCs"
       end if
    end if !  DEBUG

    if (first_call) first_call = .false.

  end subroutine BoundaryConditions

subroutine My_bcmap(iyr_trend)
! ---------------------------------------------------------------------------
! sets bc2xn_adv, bc2xn_bc, and  misc_bc
! ---------------------------------------------------------------------------
  integer, intent(in) :: iyr_trend !ds Year for which BCs are wanted
  real :: trend_ch4
  integer :: ii,i,k, io_ch4, yr_rcp= -999
  real :: decrease_factor(NGLOB_BC+1:NTOT_BC) ! Decrease factor for misc bc's
        ! Gives the factor for how much of the top-layer conc. that is left
        ! at bottom layer
  character(len=120) :: txt
  real     ::  ch4_rcp
  character(len=*),parameter :: dtxt='BCs_Mybcmap:'

  real :: top_misc_bc(NGLOB_BC+1:NTOT_BC) ! Conc. at top of misc bc
!    real :: ratio_length(KMAX_MID)    ! Vertical length of the actual layer
                                      ! divided by length from midpoint of
                                      ! layer 1 to layer KMAX_MID

  ! - Initialise
  misc_bc   = 0.0
  bc2xn_adv = 0.0
  bc2xn_bgn = 0.0

  ! Own (constant mixing ratio) boundary conditions
  ! NOTE: these species have to have the bc2xn_ indices set to 1.0 for either
  ! the advected or the background concentrations, in order that the
  ! concentrations specified in misc_bc are transferred correctly into the
  ! boundary conditions.

  ! CH4 IBCs. Default is ACP2012 values unless overridden by BGND_CH4 or 
  ! RCP settings
  ! Default:
  ! set values of 1625 in 1980, 1780 in 1990, 1820 in 2000, and 1970 in
  ! 2010. Interpolate
  ! between these for other years. Values from EMEP Rep 3/97, Table 6.2 for
  ! 1980, 1990, and from CDIAC (Mace Head) data for 2000.
  ! 2010 also from Mace Head

  if ( BGND_CH4 == -1 ) then
    if( iyr_trend >= 2010) then
      top_misc_bc(IBC_CH4) =  1870.0
    else if ( iyr_trend >= 2000) then
      top_misc_bc(IBC_CH4) = 1820 + (iyr_trend-2000)*0.1*(1870-1820) 
    else if ( iyr_trend >= 1990 ) then
      top_misc_bc(IBC_CH4) = 1780.0 + (iyr_trend-1990)*0.1*(1820-1780.0)
    else
      top_misc_bc(IBC_CH4) = 1780.0 * exp(-0.01*0.91*(1990-iyr_trend)) ! Zander,1975-1990
                                 !exp(-0.01*0.6633*(1975-iyr_trend)) ! Zander,1951-1975
    end if

  else if ( fileName_CH4_ibcs /= 'NOTSET'  ) then ! use RCP26, 45 or 85, set in config_emep

    call open_file(IO_TMP,'r',fileName_CH4_ibcs,needed=.true.)
    call CheckStop(ios,dtxt//"CH4_ibcs error in "//trim(fileName_CH4_ibcs) )
    do i=1, 9999  ! has 750 records while(.true.)
       read(IO_TMP,'(a80)') txt
       if ( txt(1:1) == '#' ) cycle
       read(txt, *) yr_rcp, ch4_rcp
       if ( yr_rcp == iyr_trend) exit
    end do
    close(IO_TMP)
    top_misc_bc(IBC_CH4) =  ch4_rcp
    if ( MasterProc ) write(*,*) dtxt//'CH4 SET from RCPs for CH4:', &
      trim(fileName_CH4_ibcs), iyr_trend, ch4_rcp
  end if

  ! Reset with namelist values if set
  if ( BGND_CH4 > 0 ) then
    if ( MasterProc ) write(*,*) dtxt//'CH4 OVERRIDE for CH4:', BGND_CH4
     top_misc_bc(IBC_CH4) = BGND_CH4
  end if

  trend_ch4 = top_misc_bc(IBC_CH4)/1780.0

  if (MasterProc) then
     write(txt,"(a,3i6,f8.1,f8.3)") dtxt//" CH4 settings (iyr,nml,ch4,trend): ",&
             iyr_trend, yr_rcp, nint(BGND_CH4), top_misc_bc(IBC_CH4),trend_ch4
     call PrintLog(txt)
  end if

  top_misc_bc(IBC_CH4) =  top_misc_bc(IBC_CH4) * PPB
  top_misc_bc(IBC_H2)  =  600.0 * PPB

  decrease_factor(IBC_H2)  = 0.0 ! No increase/decrease with height
  decrease_factor(IBC_CH4) = 0.0 ! No increase/decrease with height

! Vertical profile: Function of sigma_mid. Assumes that top_misc_bc is given for
! top(sigma_bnd(1)=0., and that at ground (sigma_bnd(KMAX_BND)=1.) the conc.
! is top_misc_bc -decrease_factor*top_misc_bc. Since the choice to set the
! concentration as a factor of sigma_mid, the concentration in the lowest
! grid cell will not be excactly top_misc_bc -decrease_factor*top_misc_bc, but close.
!pw March 2013 sigma_mid replaced by B_mid (equal for sigma coordinates)
  do ii=NGLOB_BC+1,NTOT_BC
    do k=1,KMAX_MID
      misc_bc(ii,k) = top_misc_bc(ii)*(1.0-decrease_factor(ii)*B_mid(k))
      if (MasterProc.and.DEBUG_MYBC) print "(a20,2es12.4,i4)",&
        "height,misc_vert,k",B_mid(k),misc_bc(ii,k),k
    end do
  end do

  bc2xn_adv(IBC_H2,  IXADV_H2)    = 1.0
  bc2xn_adv(IBC_CH4, IXADV_CH4)   = 1.0

  !/- completeness check
  if (DEBUG_MYBC ) then
    print *, "In My_bcmap, NGLOB_BC, NTOT_BC is", NGLOB_BC, NTOT_BC
    do i = NGLOB_BC+1 , NTOT_BC
      print *, "In My_bcmap, sum-adv", i, " is", sum(bc2xn_adv(i,:))
      print *, "In My_bcmap, sum-bgn", i, " is", sum(bc2xn_bgn(i,:))
    end do
  end if ! DEBUG

  do i = NGLOB_BC+1 , NTOT_BC
    call CheckStop(sum(bc2xn_adv(i,:))+sum(bc2xn_bgn(i,:))/=1.0,&
      "BCproblem - My_bcmap")
  end do

  !/- mappings for species from Logan + obs model given with IBC index.
  include 'BoundaryConditions_CM.inc'
end subroutine My_bcmap

subroutine Set_bcmap()
! ---------------------------------------------------------------------------
! Returns some 1-D arrays which say if a bc is used or if an
! emep xn_adv or xn_bgn is affected.  This information is derived from
! the mapping arrays bc2xn_adv and bc2xn_bgn, such that the emep species
! are given along the x-dimension and the bc species along the y.
! e.g., the statement
!   bc2xn_adv(IBC_NOX,IXADV_NO2) = 0.55
! would assign the BC concentration of NOX to the EMEP model concentration
! of NO2 after multiplication with a factor 0.55.
! These arrays have been set in the CM_BoundaryConditions.inc file.
! ---------------------------------------------------------------------------
  integer :: ibc, iem   ! local loop variables
  integer :: i

  ! Initialise
  bc_used      = 0
  bc_used_adv  = 0
  bc_used_bgn  = 0

  !/- bc_used set to one where a bc species is to be used.
  do ibc = 1, NTOT_BC
    if (any(bc2xn_adv(ibc,:)>0).or. &
        any(bc2xn_bgn(ibc,:)>0)) bc_used(ibc) = 1
    do iem = 1, NSPEC_ADV
      if(bc2xn_adv(ibc,iem)>0) bc_used_adv(ibc) = bc_used_adv(ibc)+1
    end do
    do iem = 1, NSPEC_BGN
      if(bc2xn_bgn(ibc,iem)>0) bc_used_bgn(ibc) = bc_used_bgn(ibc)+1
    end do
  end do ! ibc
  num_used_adv = maxval(bc_used_adv)
  num_used_bgn = maxval(bc_used_bgn)

  ! Initialise
  xn_adv_changed = .false.
  xn_bgn_changed = .false.
  num_adv_changed = 0
  num_bgn_changed = 0

  do iem = 1, NSPEC_ADV
    if (any(bc2xn_adv(:,iem)>0)) then
      xn_adv_changed(iem) = .true.
      num_adv_changed = num_adv_changed + 1
    end if
  end do ! iem

  do iem = 1, NSPEC_BGN
    if (any(bc2xn_bgn(:,iem)>0)) then
      xn_bgn_changed(iem) = .true.
      num_bgn_changed = num_bgn_changed + 1
    end if
  end do ! iem

  if (DEBUG%BCS) write(*,*) "TEST SET_BCMAP bc_used: ",&
    (bc_used(ibc),ibc=1, NTOT_BC)
  if (MasterProc.and.DEBUG%BCS) write(*,*)"Finished Set_bcmap: Nbcused is ", sum(bc_used)

  allocate(spc_changed2adv(num_adv_changed))
  allocate(spc_changed2bgn(num_bgn_changed))
  i = 0
  spc_adv2changed = 0
  do iem = 1, NSPEC_ADV
    if(xn_adv_changed(iem))then
      i = i+1
      spc_changed2adv(i) = iem
      spc_adv2changed(iem) = i
    end if
  end do
  i = 0
  spc_bgn2changed = 0
  do iem = 1, NSPEC_BGN
    if(xn_bgn_changed(iem))then
      i = i+1
      spc_changed2bgn(i) = iem
      spc_bgn2changed(iem) = i
    end if
  end do

  allocate(spc_used_adv(NTOT_BC,num_used_adv))
  allocate(spc_used_bgn(NTOT_BC,num_used_bgn))
  spc_used_adv = 0
  spc_used_bgn = 0
  do ibc = 1, NTOT_BC
    if ( bc_used(ibc) > 0 ) then
      ! - set bc_adv: advected species
      i = 0
      do iem = 1, NSPEC_ADV
        if (bc2xn_adv(ibc,iem)>0.0) then
          i = i+1
          spc_used_adv(ibc,i) = iem
        end if
      end do

      ! - set bc_bgn: background (prescribed) species
      i = 0
      do iem = 1, NSPEC_BGN
        if ( bc2xn_bgn(ibc,iem) > 0.0 ) then
          i = i+1
          spc_used_bgn(ibc,i) = iem
        end if
      end do
    end if    ! bc_used
  end do  ! ibc
end subroutine Set_bcmap

subroutine GetBICData(year,month,ibc,used,bc_data,errcode)
logical, parameter :: &
  DEBUG_Logan  = .false., &
  DEBUG_HZ     = .false.
! we define some concentrations in terms of sine curves and other simple data:
type :: sineconc
  real :: surf       ! Mean surface conc. (ppb)
  integer :: dmax    ! Day when concentrations peak
  real :: amp        ! amplitude of surface conc. (ppb)
  real :: hz         ! Scale-height (km) - height to drop 1/e concentration
  real :: vmin       ! background, minimum conc., in vertical direction
  real :: hmin       ! background, minimum conc., in horiz direction
  real :: conv_fac   ! factor to convert input data to mixing ratio
end type sineconc
type(sineconc), save, dimension(NGLOB_BC) :: SpecBC

type :: SIAfac ! trends in boundary conditions
  integer :: year
  real:: so2,nox,nh4
end type SIAfac

!temporary used by BoundaryConditions
real :: O3fix=0.0
real :: trend_o3=1.0, trend_co, trend_voc

! -----------------------------------------------------------------------
! HANDLES READ_IN OF GLOBAL DATA. We read in the raw data from the
! global model, and do the vertical interpolation to EMEP k values
! here if the species is to be used.
! -----------------------------------------------------------------------
  integer,             intent(in) :: year       ! for Mace Head correction
  integer,             intent(in) :: month
  integer,             intent(in) :: ibc        ! Index of BC
  integer,             intent(in) :: used       ! set to 1 if species wanted
  real, dimension(LIMAX,LJMAX,KMAX_MID), &
                      intent(out) :: bc_data   ! BC Data defined here
  integer,          intent(inout) :: errcode   ! i/o number

  logical, save :: first_call = .true.

  integer, allocatable,dimension(:,:), save :: lat5     ! for latfunc below
  real, dimension(NGLOB_BC,6:14), save  :: latfunc  ! lat. function
  real, save ::  twopi_yr, cosfac                   ! for time-variations
  !---------------------------------------------------------------------------
  ! Mace Head ozone concentrations for backgroudn sectors
  ! from Fig 5.,  Derwent et al., 1998, AE Vol. 32, No. 2, pp 145-157
  integer, parameter :: MH_YEAR1 = 1990, MH_YEAR2 = 2017
  real, dimension(12,MH_YEAR1:MH_YEAR2), parameter :: macehead_year=reshape(&
   [35.3,36.3,38.4,43.0,41.2,33.4,35.1,27.8,33.7,36.2,28.4,37.7,& !1990
    36.1,38.7,37.7,45.8,38.8,36.3,29.6,33.1,33.4,35.7,37.3,36.7,& !1991
    36.1,37.3,41.8,39.6,41.2,31.5,28.3,30.3,31.3,34.2,36.1,34.9,& !1992
    37.6,40.4,44.4,42.6,43.4,29.2,28.5,29.6,32.2,37.3,37.3,38.3,& !1993
    38.6,37.3,45.7,43.8,42.9,35.1,30.8,30.5,33.8,36.5,34.0,37.3,& !1994
    37.5,37.1,41.6,42.4,41.1,33.1,29.1,28.7,33.7,34.8,35.0,36.0,& !1995
    37.0,40.1,42.9,44.6,41.3,38.3,29.3,29.4,35.6,38.4,37.8,38.4,& !1996
    36.2,41.9,41.8,40.4,40.6,34.4,26.2,29.3,31.3,35.2,25.7,39.5,& !1997
    38.6,42.0,44.6,45.1,44.2,33.0,29.7,32.9,35.7,38.8,39.7,40.4,& !1998
    39.9,44.5,49.4,45.0,42.8,34.3,29.0,30.0,31.8,36.9,39.6,39.2,& !1999
    39.5,42.1,41.8,43.8,43.4,34.5,28.0,27.3,33.6,37.4,35.6,35.8,& !2000
    37.3,38.0,42.2,44.8,42.6,34.9,28.9,29.4,29.9,35.3,37.3,37.5,& !2001
  ! Preliminary BCs generated using Mace Head CFC and other greenhouse gases
  ! data to define clean air masses. Data cover all of 2002 and 9 months
  ! of 2003. What to do for Oct-Dec 2003?
  ! Could use (1) 2002 data or (2) 10-year average?
  ! Simmonds paper would support (1), simplicity (2).
  ! After seeing earlier 2003 plots, chose (2).
    42.4,44.4,45.5,45.0,45.9,39.8,32.5,28.7,37.7,39.3,40.5,42.3,& !2002
    39.8,40.1,44.7,45.4,45.7,41.7,33.3,31.0,35.7,37.9,40.9,38.1,& !2003
    40.8,42.0,48.3,46.6,39.9,31.9,32.4,32.1,33.9,36.7,40.2,39.8,& !2004
    40.9,41.4,44.1,45.6,42.7,32.9,26.7,30.0,33.2,37.7,39.5,38.0,& !2005
  ! 2006 and 2007 are calculated with using IE31 O3 data and
  ! trajectory sectors (based on PARLAM-PS and HIRLAM20 met) for resp. year
    39.8,42.4,44.2,48.3,41.3,39.0,31.9,29.5,34.8,37.4,41.9,39.9,& !2006
    40.7,38.2,46.1,46.4,40.9,34.5,31.2,28.8,33.3,36.1,40.6,41.7,& !2007
  ! 2008 Mace Head correction calculated using IE31 O3 data and
  ! trajectory sectors (based on HIRLAM20 met) for 2008
    41.0,45.1,48.0,46.3,44.2,37.1,30.8,31.3,34.3,37.5,37.9,40.0,& !2008
  ! 2009 to 2011 Mace Head correction calculated using IE31 O3 data and
  ! trajectory sectors (based on ECMWF met) for respective year
    37.7,43.3,46.5,46.2,41.6,39.1,31.0,29.0,34.5,34.4,40.5,38.4,& !2009
    36.8,38.9,43.9,46.4,41.7,35.5,31.0,31.3,35.6,36.7,33.4,33.8,& !2010
    36.5,42.4,43.3,44.5,40.2,34.6,30.1,30.8,32.0,34.7,37.7,38.1,& !2011
    35.0,40.2,41.0,46.8,43.1,34.0,29.6,33.8,34.9,33.3,37.9,38.7,& !2012
    38.8,42.8,45.1,46.7,43.3,31.8,31.0,33.3,32.8,39.0,39.5,42.7,& !2013
    41.4,42.9,43.5,46.4,42.4,35.1,28.6,32.6,33.8,37.1,38.1,41.1,& !2014
    41.0,43.3,43.8,42.5,39.4,33.6,31.5,35.3,35.8,42.1,40.4,41.0,& !2015
    40.4,42.5,43.7,43.6,42.4,29.7,27.5,28.6,32.0,37.7,40.5,42.5,& !2016
    41.1,45.2,46.1,45.5,40.2,33.2,28.7,32.6,34.1,39.4,41.2,39.5]& !2017
    ,[12,MH_YEAR2-MH_YEAR1+1])
  real, dimension(12), parameter :: macehead_default=&
  ! Defaults from 1998-2010 average
    (/39.8,41.9,45.4,46.5,43.2,36.2,30.5,30.1,34.1,37.0,39.0,38.5/)
  real, dimension(12):: macehead_O3=macehead_default
  !---------------------------------------------------------------------------
  integer :: i, j, k, i0, i1, Nlevel_logan, Nlevel_Dust, ierror
  real    :: f0, f1             ! interpolation factors
  character(len=30) :: fname    ! input filename
  character(len=99) :: txtmsg   ! error messages
  real,allocatable,save, dimension(:) :: p_kPa, h_km  !Use of standard atmosphere

  real :: scale_old, scale_new
  real, parameter :: macehead_lat = 53.3 !latitude of Macehead station
  real, parameter :: macehead_lon = -9.9 !longitude of Macehead station
  character(len = 100) ::varname
  real :: count_loc,O3fix_loc, mpi_rcv(2),mpi_snd(2)
  real :: conv_fac

!----------------------------------------------------------
!Trends 1980-2003 derived from EPA emissions of so2,nox.
! nh4 derived from 2/3so3+1/3nox
!Support for SO2 can be found in Hicks, Artz, Meyer and Hosker, 2002
! Figure 7 (Eastern US) which show 'close' correspondance between national
! emissions and concentration trend

!1920-1970 BCs derived from:
!  NH4: nh3 emissions
!  SOx: winter ice cores, Col du dome
!  NOx: winter ice cores
!1890-1920: trends from emissions for SOx,NOx,NH3, Aardenne USA
! Updated: April 2013
! - use data above to 1980, then EPA download of April 2013
! - then IIASA/ECLAIRE/ECLIPSE
 type(SIAfac), dimension(37), save :: SIAtrends = (/ &
    SIAfac(1890,0.12,0.15,0.44) &
   ,SIAfac(1900,0.18,0.20,0.48) &
   ,SIAfac(1910,0.27,0.27,0.52) &
   ,SIAfac(1920,0.32,0.33,0.59) &
   ,SIAfac(1930,0.35,0.33,0.55) &
   ,SIAfac(1940,0.46,0.25,0.59) &
   ,SIAfac(1950,0.59,0.33,0.69) &
   ,SIAfac(1960,0.76,0.50,0.76) &
   ,SIAfac(1970,0.95,0.75,0.90) &
   ,SIAfac(1980,   1.000,   1.000,   1.000)&
   ,SIAfac(1985,   0.899,   0.951,   0.989)&
   ,SIAfac(1990,   0.890,   0.943,   0.920)&
   ,SIAfac(1991,   0.863,   0.930,   0.934)&
   ,SIAfac(1992,   0.852,   0.933,   0.947)&
   ,SIAfac(1993,   0.840,   0.936,   0.963)&
   ,SIAfac(1994,   0.823,   0.936,   0.978)&
   ,SIAfac(1995,   0.718,   0.922,   0.993)&
   ,SIAfac(1996,   0.709,   0.915,   1.007)&
   ,SIAfac(1997,   0.727,   0.912,   1.027)&
   ,SIAfac(1998,   0.731,   0.899,   1.052)&
   ,SIAfac(1999,   0.677,   0.844,   1.035)&
   ,SIAfac(2000,   0.631,   0.835,   1.046)&
   ,SIAfac(2001,   0.615,   0.796,   0.786)&
   ,SIAfac(2002,   0.570,   0.781,   0.880)&
   ,SIAfac(2003,   0.568,   0.753,   0.877)&
   ,SIAfac(2004,   0.565,   0.726,   0.874)&
   ,SIAfac(2005,   0.572,   0.703,   0.870)&
   ,SIAfac(2006,   0.514,   0.681,   0.885)&
   ,SIAfac(2007,   0.456,   0.658,   0.900)&
   ,SIAfac(2008,   0.399,   0.635,   0.930)&
   ,SIAfac(2009,   0.320,   0.579,   0.928)&
   ,SIAfac(2010,   0.292,   0.543,   0.925)&
   ,SIAfac(2011,   0.265,   0.488,   0.921)&
   ,SIAfac(2012,   0.213,   0.421,   0.917)&
   ! Default here from IIASA ECLAIRE/ECLIPSE
   ! related to 2005 emissions as base
   ! (Created by mk.UStrends, April 2013)
   ,SIAfac(2030,   0.155,   0.252,   0.953)&
   ,SIAfac(2050,   0.225,   0.276,   0.977)&! Last year which works
   ,SIAfac(2200,   0.225,   0.276,   0.977)&! FAKE for interp
 /)
  type(SIAfac), save :: SIAtrend

  if (iyr_trend < SIAtrends(1)%year .or. iyr_trend >= SIAtrends(37)%year ) then
    write(unit=txtmsg,fmt=*) "Unspecified trend BCs for this year:", ibc, year
    call CheckStop(txtmsg)
  end if

!================================================================== 
! Interpolate between boundary condition years if needed
  i0 = 1
  do i = 1, size(SIAtrends(:)%year) 
    if ( iyr_trend >= SIAtrends(i)%year ) i0 = i
    !if(MasterProc) print "(a,5i6)", "USAsrch: ", i, i0, BCtrend(i)%year, BCtrend(i0)%year
  end do

  i1= i0 + 1
  f0 =     (SIAtrends(i1)%year - iyr_trend)/&
       real(SIAtrends(i1)%year - SIAtrends(i0)%year )
  f1 =     (iyr_trend        - SIAtrends(i0)%year )/&
       real(SIAtrends(i1)%year - SIAtrends(i0)%year )

  SIAtrend%so2 =f0*SIAtrends(i0)%so2 + f1*SIAtrends(i1)%so2
  SIAtrend%nox =f0*SIAtrends(i0)%nox + f1*SIAtrends(i1)%nox
  SIAtrend%nh4 =f0*SIAtrends(i0)%nh4 + f1*SIAtrends(i1)%nh4

!  if (MasterProc.and.first_call) then
!     write(unit=txtmsg,fmt="(a,2i5,3f8.3)") &
!       "BC:trends SOx,NOx,NH3 ", iyr_trend, SIAtrend
!     call PrintLog(txtmsg)
!  end if

!================================================================== 
! Trends - derived from EMEP report 3/97
! adjustment for years outside the range 1990-2000.

!June 2013 svn 2619:
! Assume O3 increases to 2000, consistent with obs.
! Keep the 1990 base-year for CO and VOC.
  select case(iyr_trend)
  case(2000:)
    trend_o3 = 1.0
    trend_co = 1.0
    trend_voc= 1.0
  case(1990:1999)
  if( USES%MACEHEADFIX ) then
       trend_o3 = 1.0
    else
       trend_o3 = exp(-0.01*1.0 *(2000-iyr_trend))
    end if
    trend_co = 1.0
    trend_voc= 1.0
  case default
    trend_o3 = exp(-0.01*1.0 *(2000-iyr_trend))
    trend_co = exp(-0.01*0.85*(1990-iyr_trend)) ! Zander:CO
    trend_voc= exp(-0.01*0.85*(1990-iyr_trend)) ! Zander,1975-1990
  end select
  if (MasterProc.and.first_call) then
    write(unit=txtmsg,fmt="(a,i5,3f8.3,13f9.4)") "BC:trends O3,CO,VOC,SOx,NOx,NH3: ", &
       iyr_trend, trend_o3, trend_co, trend_voc, SIAtrend%so2, SIAtrend%nox, SIAtrend%nh4
    call PrintLog(txtmsg)
  end if

!=========== BCs Generated from Mace Head Data ====================
!
! Here we use the meteorology year to get a reaslistic O3.
!      Later we use iyr_trend to adjust for other years, say for 2050.
! For 2020 "trend" runs  - use 13 yr average as base-O3 (macehead_default)
! then later scale by trend_o3:
  if((iyr_trend==year).and.(year>=MH_YEAR1).and.(year<=MH_YEAR2))then
    macehead_O3=macehead_year(:,year)
    write(unit=txtmsg,fmt="(a,i5)") "BC: O3 Mace Head correction for year ", year
  else
    macehead_O3=macehead_default
    write(unit=txtmsg,fmt="(a)") "BC: O3 default Mace Head correction"
  end if
  if (MasterProc.and.first_call) then
    call PrintLog(txtmsg)
  end if
!=========== Generated from Mace Head Data =======================

  errcode = 0
  txtmsg = "ok"
  if (DEBUG_Logan) print *,"DEBUG_LOgan ibc, mm", ibc, month

! ========= first call =========================================
  if ( first_call ) then
    ! Set up arrays to contain Logan's grid as lat/long
    !/ COnversions derived from emeplat2Logan etc.:
     allocate(lat5(LIMAX,LJMAX))
     allocate(p_kPa(KMAX_MID), h_km(KMAX_MID))
    twopi_yr = 2.0 * PI / 365.25

!    call GlobalPosition  !get glat for global domaib
    forall(i=1:limax,j=1:ljmax) ! Don't bother with south pole complications
      lat5(i,j) = glat(i,j)/5   ! lat/5 used in latfunc below
      lat5(i,j) = max(lat5(i,j),6)   ! Min value in latfunc
      lat5(i,j) = min(lat5(i,j),14)  ! Max value in latfunc
    endforall
    ! Define concs where a simple  specification based on lat/mm
    ! etc. will be given
    !                           surf   dmax   amp   hz    vmin  hmin conv_fac!ref
    !                            ppb          ppb   km   hmin,vmin:same units as input data=conv_fac
    SpecBC(IBC_SO2  )  = sineconc( 0.15 , 15.0, 0.05, 999.9, 0.03, 0.03,PPB)!W99, bcKz vmin
!pwds    SpecBC(IBC_SO4  )  = sineconc( 0.15 ,180.0, 0.00, 1.6  , 0.05, 0.03,PPB)!W99
    SpecBC(IBC_SO4  )  = sineconc( 0.15 ,180.0, 0.00, 999.9, 0.05, 0.03,PPB)!W99
    SpecBC(IBC_NO   )  = sineconc( 0.1  , 15.0, 0.03, 4.0  , 0.03, 0.02,PPB)
    SpecBC(IBC_NO2  )  = sineconc( 0.1  , 15.0, 0.03, 4.0  , 0.05, 0.04,PPB)
    SpecBC(IBC_PAN  )  = sineconc( 0.20 ,120.0, 0.15, 999.9, 0.20, 0.1 ,PPB)!Kz change vmin
    SpecBC(IBC_CO   )  = sineconc( 125.0, 75.0, 35.0, 25.0 , 70.0, 30.0,PPB)!JEJ-W
!st 14.05.2014    SpecBC(IBC_SeaSalt_F)=sineconc( 0.5  , 15.0,  0.3,  1.6 , 0.01, 0.01,PPB)
!st 14.05.2014    SpecBC(IBC_SeaSalt_C)=sineconc( 3.0  , 15.0,  1.0,  1.6 , 0.01, 0.01,PPB)
    SpecBC(IBC_SeaSalt_f)=sineconc( 0.2  , 15.0,  0.05,  1.6 , 0.01, 0.01,PPB)
    SpecBC(IBC_SeaSalt_c)=sineconc( 1.5  , 15.0,  0.25,  1.6 , 0.01, 0.01,PPB)
    SpecBC(IBC_SeaSalt_g)=sineconc( 1.0  , 15.0,  0.5,  1.0 , 0.01, 0.01,PPB)
    SpecBC(IBC_C2H6 )  = sineconc( 2.0  , 75.0, 1.0 , 10.0 , 0.05, 0.05,PPB)
    SpecBC(IBC_C4H10)  = sineconc( 2.0  , 45.0, 1.0 , 6.0  , 0.05, 0.05,PPB)
    SpecBC(IBC_HCHO )  = sineconc( 0.7  ,180.0, 0.3 , 6.0  , 0.05, 0.05,PPB)
    SpecBC(IBC_CH3CHO) = sineconc( 0.3  ,180.0, 0.05 , 6.0  , 0.005, 0.005,PPB) !NAMBLEX,Solberg,etc.
    SpecBC(IBC_HNO3 )  = sineconc( 0.07 ,180.0, 0.03, 999.9,0.025, 0.03,PPB)
                         !~=NO3, but with opposite seasonal var.
    SpecBC(IBC_NO3_f ) = sineconc( 0.07 , 15.0, 0.03, 1.6  ,0.025, 0.02,PPB) !ACE-2
    SpecBC(IBC_NO3_c ) = sineconc( 0.07 , 15.0, 0.00, 1.6  ,0.025, 0.02,PPB) !ACE-2
    SpecBC(IBC_NH4_f ) = sineconc( 0.15 ,180.0, 0.00, 1.6  , 0.05, 0.03,PPB) !ACE-2(SO4/NH4=1)
 ! all BCs read in are in mix. ratio, thus hmin,vmin needs to be in mix. ratio for thosetio for those
    SpecBC(IBC_O3   )  = sineconc(-99.9 ,-99.9,-99.9,-99.9 ,-99.9,10.0*PPB  ,1.)!N1
    SpecBC(IBC_H2O2 )  = sineconc(-99.9 ,-99.9,-99.9,-99.9 ,-99.9,0.01*PPB  ,1.)
  ! Dust: the factor PPB converts from PPB to mixing ratio.
    SpecBC(IBC_Dust_c)=sineconc(-99.9 ,-99.9,-99.9,-99.9 ,-99.9,1.0e-15,1.0)
    SpecBC(IBC_Dust_f)=sineconc(-99.9 ,-99.9,-99.9,-99.9 ,-99.9,1.0e-15,1.0)
    !refs:
    ! N1 - for ozone we read Logan's data, so the only paramater specified
    !      is a min value of 10 ppb. I hope this doesn't come into effect in
    !      Europe as presumably any such min values are in the S. hemisphere.
    !      Still, giving O3 such a value let's us use the same code for
    !      all species.
    ! W99: Warneck, Chemistry of the Natural Atmosphere, 2nd edition, 1999
    !    Academic Press. Fig 10-6 for SO2, SO4.
    ! JEJ - Joffen's suggestions from Mace/Head, UiO and other data..
    !      with scale height estimated large from W99, Isaksen+Hov (1987)
    ! M -Mozart-obs comparison
    ! ACE-2 Lots of conflicting measurements exist, from NH4/SO4=2 to NH4/SO4=0.5
    ! A 'mean' value of NH4/SO4=1 is therefore selected. Otherwise NH4 is assumed to
    ! act as SO4
    ! aNO3 is assumed to act as SO4, but with 1/2 concentrations and seasonal var.
    ! pNO3 is assumed to act like seasalt, with decreasing conc with height, with
    ! approx same conc. as fine nitrate.
    ! The seasonal var of HNO3 is now assumed to be opposite of aNO3.

    ! Consistency check:
    if (DEBUG%GLOBBC) print *, "SPECBC NGLB ",NGLOB_BC
    do i = 1, NGLOB_BC
      if (DEBUG%GLOBBC) print *,"SPECBC i, hmin ",i,SpecBC(i)%surf,SpecBC(i)%hmin
      if( SpecBC(i)%hmin*SpecBC(i)%conv_fac < 1.0e-17) then
        write(unit=txtmsg,fmt="(A,I0)") "PECBC: Error: No SpecBC set for species ", i
        call CheckStop(txtmsg)
      end if
    end do

    ! Latitude functions taken from Lagrangian model, see Simpson (1992)
    latfunc(:,6:14) = 1.0    ! default
    if(me==0)write(*,*)'WARNING SET LATFUNC TO CONSTANT 1'
!Dave, Peter 10th Feb 2015: simplify and set to 1!
                              !  30        40        50       60         70 degN
!    latfunc(IBC_SO2 ,6:14) = (/ 0.05,0.15,0.3 ,0.8 ,1.0 ,0.6 ,0.2 ,0.12,0.05/)
!    latfunc(IBC_HNO3,6:14) = (/ 1.00,1.00,1.00,0.85,0.7 ,0.55,0.4 ,0.3 ,0.2 /)
!    latfunc(IBC_PAN ,6:14) = (/ 0.15,0.33,0.5 ,0.8 ,1.0 ,0.75,0.5 ,0.3 ,0.1 /)
!    latfunc(IBC_CO  ,6:14) = (/ 0.6 ,0.7 ,0.8 ,0.9 ,1.0 ,1.0 ,0.95,0.85,0.8 /)
!
!    latfunc(IBC_SO4   ,:) = latfunc(IBC_SO2 ,:)
!    latfunc(IBC_NO    ,:) = latfunc(IBC_SO2 ,:)
!    latfunc(IBC_NO2   ,:) = latfunc(IBC_SO2 ,:)
!    latfunc(IBC_HCHO  ,:) = latfunc(IBC_HNO3,:)
!    latfunc(IBC_CH3CHO,:) = latfunc(IBC_HNO3,:)
!    latfunc(IBC_NH4_f ,:) = latfunc(IBC_SO2 ,:)
!    latfunc(IBC_NO3_f ,:) = latfunc(IBC_SO2 ,:)
!    latfunc(IBC_NO3_c ,:) = latfunc(IBC_SO2 ,:)

    ! Use Standard Atmosphere to get average heights of layers
    p_kPa(:) = 0.001*( A_mid(:) + B_mid(:)*Pref ) ! Pressure in kPa
    h_km = StandardAtmos_kPa_2_km(p_kPa)

    first_call = .false.
  end if ! first_call
! ========= end of first call ===================================
!+
!  Specifies concentrations for a fake set of Logan data.

  fname = "none"           ! dummy for printout
  select case (ibc)
  case (IBC_O3)
     Nlevel_logan=30
     if(.not.allocated(O3_logan))allocate(O3_logan(Nlevel_logan,LIMAX,LJMAX))
     if(.not.allocated(O3_logan_emep))allocate(O3_logan_emep(LIMAX,LJMAX,KMAX_MID))
     O3_logan=0.0
     O3_logan_emep=0.0
 
     varname='O3'
     if(me==0)write(*,*)'reading IBC for O3 from ',trim(LoganO3File)
     call  ReadField_CDF(LoganO3File,varname,O3_logan,nstart=month,kstart=1,kend=Nlevel_logan,interpol='zero_order', &
          needed=.true.,debug_flag=.false.)
     !interpolate vertically
     call vertical_interpolate(LoganO3File,O3_logan,Nlevel_logan,O3_logan_emep,debug=.false.)
     do k = 1, KMAX_MID
        do j = 1, ljmax
           do i = 1, limax
              !                   bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)=O3_logan_emep(i,j,k)
              bc_data(i,j,k)=O3_logan_emep(i,j,k)
              
           end do
        end do
     end do
    ! Mace Head adjustment: get mean ozone from Eastern sector
    O3fix_loc=0.0
    count_loc=0
    if(USES%MACEHEADFIX)then
       do j=1,ljmax
          do i=1,limax
             if(glat(i,j)<macehead_lat+20.0.and.&
                glat(i,j)>macehead_lat-25.0.and.&
                glon(i,j)<macehead_lon     .and.&
                glon(i,j)>macehead_lon-40.0)then
                O3fix_loc=O3fix_loc+bc_data(i,j,KMAX_MID)
                count_loc=count_loc+1
             end if
          end do
       end do
       mpi_snd(1)=O3fix_loc
       mpi_snd(2)=count_loc
       call MPI_ALLREDUCE(mpi_snd, mpi_rcv, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
       O3fix=0.0
       if(mpi_rcv(2)>0.5)O3fix=mpi_rcv(1)/mpi_rcv(2) - macehead_O3(month)*PPB
       if (me==0)write(*,"(a,4f8.3)")'Mace Head correction for O3, trend and Mace Head value',&
            -O3fix/PPB,trend_o3,macehead_O3(month)
       bc_data = max(15.0*PPB,bc_data-O3fix)
    end if
  case ( IBC_H2O2 )

     bc_data=1.0E-25

  case (IBC_NO  ,IBC_NO2  ,IBC_HNO3,IBC_CO, &
        IBC_C2H6,IBC_C4H10,IBC_PAN ,IBC_NO3_c,&
        IBC_SO2   , IBC_SO4  , IBC_HCHO , &
        IBC_SeaSalt_f,IBC_SeaSalt_c, IBC_SeaSalt_g, &
        IBC_CH3CHO, IBC_NH4_f, IBC_NO3_f)
    ! NB since we only call once per month we add 15 days to
    ! day-number to get a mid-month value
    cosfac = cos( twopi_yr * (daynumber+15.0-SpecBC(ibc)%dmax))
    bc_data(:,:,KMAX_MID) = SpecBC(ibc)%surf + SpecBC(ibc)%amp*cosfac

    if(SpecBC(ibc)%hz<100.0)then
    !/ - correct for other heights
    do k = 1, KMAX_MID-1
      scale_new = exp( -h_km(k)/SpecBC(ibc)%hz )
      bc_data(:,:,k) = bc_data(:,:,KMAX_MID)*scale_new
      if (DEBUG_HZ) then
        scale_old = exp( -(KMAX_MID-k)/SpecBC(ibc)%hz )
        write(*,"(a8,2i3,2f8.3,i4,f8.2,f8.3,2f8.3)") &
         "SCALE-HZ ", month, ibc, SpecBC(ibc)%surf, SpecBC(ibc)%hz, k,&
          h_km(k), p_kPa(k), scale_old, scale_new
      end if ! DEBUG_HZ
    end do

    else    
       do k = 1, KMAX_MID-1
          bc_data(:,:,k) = bc_data(:,:,KMAX_MID)
       end do
    end if

    !/ - min value after vertical factors, before latitude factor
    bc_data = max( bc_data, SpecBC(ibc)%vmin )

    !/ - correct for latitude functions
    forall(i=1:LIMAX,j=1:LJMAX)
      bc_data(i,j,:) = bc_data(i,j,:) * latfunc(ibc,lat5(i,j))
    endforall

    case (IBC_DUST_C,IBC_DUST_F)
       if(USES%DUST)then
!         bc_data(:,:,:) = 0.0

!dust are read from the results of a Global run
         Nlevel_Dust=20
         if(.not.allocated(Dust_3D))allocate(Dust_3D(Nlevel_Dust,LIMAX,LJMAX))
         if(.not.allocated(Dust_3D_emep))allocate(Dust_3D_emep(LIMAX,LJMAX,KMAX_MID))
         Dust_3D=0.0
         Dust_3D_emep=0.0
 
         if(ibc==IBC_DUST_C)then
            varname='D3_ug_DUST_WB_C'  ! QUERY : ARGH!   RECODE more flexibly!
          if(me==0)write(*,*)'coarse DUST BIC read from climatological file'
         else if(ibc==IBC_DUST_F)then
            varname='D3_ug_DUST_WB_F'
            if(me==0)write(*,*)'fine DUST BIC read from climatological file'
         else
            call CheckStop('IBC dust case error')
         end if
         call  ReadField_CDF(DustFile,varname,Dust_3D,nstart=month,kstart=1,kend=Nlevel_Dust,&
              interpol='zero_order', needed=.true.,debug_flag=.false.)

         !interpolate vertically
         call vertical_interpolate(DustFile,Dust_3D,Nlevel_Dust,Dust_3D_emep,debug=.false.)

!has to convert from ug/m3 into mixing ratio. NB: Dust in Netcdf file has molwt = 200 g/mol
         conv_fac=ATWAIR/200.*1.E-9
         do k = 1, KMAX_MID
            do j = 1, ljmax
               do i = 1, limax
                  bc_data(i,j,k)=Dust_3D_emep(i,j,k)*conv_fac/roa(i,j,k,1)           
               end do
            end do
         end do
         else
            bc_data=0.0
         end if

    case  default
      print *,"Error with specified BCs:", ibc
      txtmsg = "BC Error UNSPEC"
    end select
!================== end select ==================================

  call CheckStop(txtmsg)
  if (DEBUG_Logan) then
    print "(a15,3i4,f8.3)","DEBUG:LOGAN: ",ibc, used, month, cosfac
    print *,"LOGAN BC MAX ", maxval ( bc_data ), &
                    " MIN ", minval ( bc_data )
    do k = KMAX_MID, 1, -1        ! print out a random column
      print "(i4,f12.3)", k, bc_data(5,5,k)
    end do
  end if ! DEBUG


  !/ - min value after latitude factors , but before trends
  bc_data = max( bc_data, SpecBC(ibc)%hmin )


    !/ trend adjustments
   select case (ibc)
   case ( IBC_O3 )
      bc_data = bc_data*trend_o3
   case (IBC_C4H10 , IBC_C2H6 )
      bc_data =  bc_data*trend_voc
   case ( IBC_CO )
      bc_data =  bc_data*trend_co
   case ( IBC_SO2,IBC_SO4)
      bc_data = bc_data*SIAtrend%so2
   case( IBC_NH4_f)
      bc_data = bc_data*SIAtrend%nh4
   case ( IBC_NO3_f,IBC_NO3_c,IBC_HNO3,IBC_NO2,IBC_NO,IBC_PAN)
      bc_data = bc_data*SIAtrend%nox
   end select

   bc_data = bc_data * SpecBC(ibc)%conv_fac !Convert to mixing ratio


end subroutine GetBICData

end module BoundaryConditions_mod
