! <BoundaryConditions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
! -----------------------------------------------------------------------
module BoundaryConditions_ml
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
! global input data, and rsets the concentrations at the top level and
! lateral boundaries for advected species. For background species it reset
! full 3-D concentration fields.
!
! A bilinear interpolation routine is used to extrapolate from the coarser
! global data to the EMEP model arrays, from the module Interpolations_ml.
! Vertical interpolation is done on reading in the data.
!
! The main related IC/BC modules/files are:
!
!   1. GlobalBCs_ml.f90 - sets indices of global model data, e.g. IBC_O3,
!      as well as the number of global model fdata (NGLOB_BC).
!
!   2. CM_BoundaryConditions.inc - assigns mappings, telling which
!      Unified EMEP model species the global model fdata are
!      assigned to (bc2xn_adv, bc2xn_bgn arrays).
!
!   3. Nest_ml.f90 - In nested runs (such as FORECAST mode), the ICs & BCs
!      are reseted by readxn (Nest_ml), superseeding the IC/BCs in this module.
! -----------------------------------------------------------------------
! IMPORTANT NOTES:
! 1. The routines given here are constructed around the global model
!    fields from the University of Oslo (T21) global model. In order
!    to use other models as BCs then usually these routines will have to be
!    replaced by model-specific routines. The important thing is
!    that the inputs and outputs from the routine are independant of the
!    global module ufor one bc speciessed.
! 2. The routines make use of a "feature" of the model: that the concentration
!    (xn) values  along boundaries are not changed due to advection or chemistry.
!    Thus, bc values only need to be set once per month, firstly for the whole 3-D
!    domain, then monthly for the sides and top.
!    Background species must be reset in 3-D each month.
! 3. Time varying BCs and model specific ICs
! -----------------------------------------------------------------------

use CheckStop_ml,      only: CheckStop
use Chemfields_ml,     only: xn_adv, xn_bgn, NSPEC_BGN  ! emep model concs.
use ChemSpecs                ! provide NSPEC_ADV and IXADV_*
!CMR use ChemChemicals_ml         ! provide species names
!CMR use ChemSpecs_adv_ml         ! provide NSPEC_ADV and IXADV_*
!CMR use ChemSpecs_shl_ml         ! provide NSPEC_SHL
use GlobalBCs_ml,      only:  &
   NGLOB_BC                   &  ! Number of species from global-model
  ,GetGlobalData              &  ! Sub., reads global data+vert interp.
  ,IBC_SO2, IBC_SO4, IBC_HCHO, IBC_CH3CHO &
  ,IBC_O3,IBC_HNO3,IBC_PAN,IBC_CO,IBC_C2H6   &
  ,IBC_C4H10, IBC_NO ,IBC_NO2,IBC_NH4_f,IBC_NO3_f,IBC_NO3_c&
  ,IBC_H2O2, IBC_DUST_f, IBC_DUST_c,IBC_SEASALT_F, IBC_SEASALT_C &
  ,IBC_SEASALT_G ,setgl_actarray&
  ,O3fix,trend_o3!temporary
use GridValues_ml,     only: glon, glat   & ! full domain lat, long
                            ,sigma_mid    & !sigma layer midpoint
                            ,debug_proc, debug_li, debug_lj & ! debugging
                            ,i_fdom, j_fdom,B_mid  ! for debugging
use Io_Progs_ml,       only: datewrite, PrintLog
use Landuse_ml,        only: mainly_sea
use LocalVariables_ml, only: Grid
use MetFields_ml,      only: z_mid      ! height of half layers
use ModelConstants_ml, only: KMAX_MID  &  ! Number of levels in vertical
                            ,iyr_trend &  ! Used for e.g. future scenarios
                            ,BGND_CH4  &  ! If positive, replaces defaults
                            ,USE_SEASALT & 
                            ,USES,DEBUG  & ! %BCs
                            ,MasterProc, PPB
use NetCDF_ml,         only:ReadField_CDF,vertical_interpolate
use Par_ml,          only: &
   MAXLIMAX, MAXLJMAX, limax, ljmax, me &
  ,neighbor, NORTH, SOUTH, EAST, WEST   &  ! domain neighbours
  ,NOPROC&
  ,IRUNBEG,JRUNBEG,li1,li0,lj0,lj1
use SmallUtils_ml, only : find_index 

implicit none
private

! -- subroutines in this module:
public  :: BoundaryConditions         ! call every month
private :: Set_bcmap,               & ! sets xn2adv_changed, etc.
           MiscBoundaryConditions,  & ! misc bcs, not from global model.
           Set_BoundaryConditions,  & ! assigns concentrations (xn) from bcs
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

INCLUDE 'mpif.h'
INTEGER STATUS(MPI_STATUS_SIZE),INFO

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
    integer :: ibc, iem, k, iem1, i, j ,n, nadv,ntot ! loop variables
    integer :: info                     ! used in rsend
    integer :: io_num                 !  i/o number used for reading global data
    integer :: alloc_err
    real    :: bc_fac      ! Set to 1.0, except sea-salt over land = 0.01
    logical :: bc_seaspec  ! if sea-salt species

    !/ data arrays for boundary data (BCs) - quite large, so NOT saved
    real, allocatable,dimension(:,:,:)   :: bc_data   ! for one bc species
!    real, allocatable,dimension(:,:,:,:) :: bc_adv,bc_bgn
    ! Dimensions correspond to:
    !   bc_data(IGLOB,JGLOB,KMAX_MID)
    !   bc_adv(NSPEC_ADV,IGLOB,JGLOB,KMAX_MID)
    !   bc_bgn(NSPEC_BGN,IGLOB,JGLOB,KMAX_MID)

    integer  :: iglobact, jglobact, errcode, Nlevel_logan
    integer, save :: idebug=0, itest=1, i_test=0, j_test=0
    real, allocatable,dimension(:,:,:)   :: O3_logan,O3_logan_emep
    character(len = 100) ::fileName,varname
    logical :: NewLogan=.false.! under testing

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

    endif ! first call
    if (DEBUG%BCS) write(*, "((A,I0,1X))")           &
         "CALL TO BOUNDARY CONDITIONS, me:", me, &
         "month ", month, "TREND2 YR ", iyr_trend, "me ", me

    if (num_changed==0) then
       write(*,*) "BCs: No species requested"
       return
    endif

    !MUST CONTAIN DECIDED DIMENSION FOR READ-IN DATA
    ! iglobac and jglobac are now the actual domains (the chosen domain)
    ! given in the same coord as the data we read
    call setgl_actarray(iglobact,jglobact)

    allocate(bc_data(iglobact,jglobact,KMAX_MID),stat=alloc_err)
    call CheckStop(alloc_err, "alloc1 failed in BoundaryConditions_ml")

!    allocate(bc_adv(num_adv_changed,iglobact,jglobact,KMAX_MID),stat=alloc_err)
!    call CheckStop(alloc_err, "alloc2 failed in BoundaryConditions_ml")
!    bc_adv(:,:,:,:) = 0.0

!    allocate(bc_bgn(num_bgn_changed,iglobact,jglobact,KMAX_MID),stat=alloc_err)
!    call CheckStop(alloc_err, "alloc3 failed in BoundaryConditions_ml")
!    bc_bgn(:,:,:,:) = 0.0

    errcode = 0
    if (DEBUG%BCS.and.debug_proc) then
       do i = 1, limax
          do j = 1, ljmax
             if (i_fdom(i)==DEBUG%IJ(1).and.j_fdom(j)==DEBUG%IJ(2)) then
                i_test = i
                j_test = j
             endif
          enddo
       enddo
    endif

    if (first_call) then
       idebug = 1
       if (DEBUG%BCS) write(*,*) "RESET 3D BOUNDARY CONDITIONS", me
       do k = 1, KMAX_MID
          do j = 1, ljmax
             do i = 1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             enddo
          enddo
       enddo
    else       
       if (DEBUG%BCS.and.MasterProc) write(*,*) "RESET LATERAL BOUNDARIES"
       do k = 2, KMAX_MID
          do j = lj0, lj1
             !left
             do i = 1, li0-1
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             enddo
             !right
             do i = li1+1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             enddo
          enddo
          !lower
          do j = 1, lj0-1
             do i = 1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             enddo
          enddo
          !upper
          do j = lj1+1, ljmax
             do i = 1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             enddo
          enddo
       enddo
       !top
       do k = 1, 1
          do j = 1, ljmax
             do i = 1, limax
                xn_adv(:,i,j,k)=0.0
                xn_bgn(:,i,j,k)=0.0
             enddo
          enddo
       enddo
    endif
    !== BEGIN READ_IN OF GLOBAL DATA

    do ibc = 1, NGLOB_BC
      
       if (MasterProc) call GetGlobalData(year,month,ibc,bc_used(ibc), &
            iglobact,jglobact,bc_data,io_num,errcode)
       
       if (DEBUG%BCS.and.MasterProc) &
            write(*, *)'Calls GetGlobalData: year,iyr_trend,ibc,month,bc_used=', &
            year,iyr_trend,ibc,month,bc_used(ibc)
       
       call CheckStop(ibc==1.and.errcode/= 0,&
            "ERRORBCs: GetGlobalData, failed in BoundaryConditions_ml")
       
       !-- If the read-in bcs are required, we broadcast and use:
       if ( bc_used(ibc) > 0 ) then
          CALL MPI_BCAST(bc_data,8*iglobact*jglobact*KMAX_MID,MPI_BYTE,0,&
               MPI_COMM_WORLD,INFO)
          
          ! - set bc_adv: advected species
          !          do i = 1, bc_used_adv(ibc)
          !             iem = spc_used_adv(ibc,i)
          !             iem1 = spc_adv2changed(iem)
          !             bc_adv (iem1,:,:,:) = bc_adv(iem1,:,:,:) &
          !                  + bc_data(:,:,:)*bc2xn_adv(ibc,iem)
          !          enddo
          
          ! - set bc_bgn: background (prescribed) species
          !          do i = 1, bc_used_bgn(ibc)
          !             iem = spc_used_bgn(ibc,i)
          !             iem1 = spc_bgn2changed(iem)
          !             bc_bgn(iem1,:,:,:) = bc_bgn(iem1,:,:,:) &
          !                  +  bc_data(:,:,:)*bc2xn_bgn(ibc,iem)
          !          enddo
       endif    ! bc_used
       
       !   if (MasterProc) close(io_num)

       if(ibc==IBC_O3 .and. (NewLogan.or.KMAX_MID/=20))then !temporary fix, assumes IBC_O3=1
!This should have been in GetGlobalData, but GetGlobalData is called only by MasterPoroc.
!So we overwrite whatever O3 is read in from  GetGlobalData
          if(Masterproc)write(*,*)'OVERWRITING LOGAN'
 
          !Read Logan BC in pressure coordinates
          Nlevel_logan=30
          if(.not.allocated(O3_logan))allocate(O3_logan(Nlevel_logan,MAXLIMAX,MAXLJMAX))
          if(.not.allocated(O3_logan_emep))allocate(O3_logan_emep(MAXLIMAX,MAXLJMAX,KMAX_MID))
          filename='Logan_P.nc'!will be put in run.pl in due time
          varname='O3'
          call  ReadField_CDF(fileName,varname,O3_logan,nstart=month,kstart=1,kend=Nlevel_logan,interpol='zero_order', &
              needed=.true.,debug_flag=.true.)
          CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
           !interpolate vertically
          call vertical_interpolate(filename,O3_logan,Nlevel_logan,O3_logan_emep,Masterproc)
          do k = 1, KMAX_MID
             do j = 1, ljmax
                do i = 1, limax
                   bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)=O3_logan_emep(i,j,k)
                enddo
             enddo
          enddo
          if(USES%MACEHEADFIX)then
             !MaceHead correction
             CALL MPI_BCAST(trend_o3,8,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
             CALL MPI_BCAST(O3fix,8,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
             if(masterProc)write(*,*)'O3fix,trend_o3 ',O3fix,trend_o3
             bc_data = max(15.0*PPB,bc_data-O3fix)
             bc_data = bc_data*trend_o3
          endif
       endif

       if (first_call) then

          ! Set 3-D arrays of new BCs
          do n = 1, bc_used_adv(ibc)

             iem = spc_used_adv(ibc,n)
             ntot = iem + NSPEC_SHL 

            ! Sea-salt. 
            !  If SeaSalt isn't called from mk.GenChem, we don't have the
            !  SS_GROUP, so we search for the simple SEASALT name.
             bc_seaspec = .false.
             if ( USE_SEASALT .and. &
                  ( index( species(ntot)%name, "SEASALT_" ) > 0 ) ) then
                bc_seaspec = .true.
             end if

             if ( debug_proc ) write (*,*) "SEAINDEX", &
                  trim(species(ntot)%name), n, ntot, bc_seaspec,&
                       index( species(ntot)%name, "SEASALT_")

             do k = 1, KMAX_MID
                do j = 1, ljmax
                   do i = 1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USE_SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)
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
                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                      !                      !        bc_bgn(n,(i_fdom(i)-IRUNBEG+1),(j_fdom(j)-JRUNBEG+1),k)
                   end do ! i
                end do ! j
             end do ! k
          enddo
       else

          ! Set LATERAL (edge and top) arrays of new BCs

          !       call MiscBoundaryConditions(iglobact,jglobact,bc_adv,bc_bgn)
          !       call Set_BoundaryConditions("lateral",iglobact,jglobact,bc_adv,bc_bgn)
          idebug = idebug + 1
          do n = 1, bc_used_adv(ibc)
             iem = spc_used_adv(ibc,n)
             ntot = iem + NSPEC_SHL 
             bc_seaspec = .false.
             if ( USE_SEASALT .and. ( index( species(ntot)%name, "SEASALT_" ) > 0 ) ) then
                bc_seaspec = .true.
             end if

             do k = 2, KMAX_MID
                do j = lj0, lj1
                   !left
                   do i = 1, li0-1
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USE_SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)

                   enddo
                   !right
                   do i = li1+1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USE_SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)

                   enddo
                enddo
                !lower
                do j = 1, lj0-1
                   do i = 1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USE_SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)

                   enddo
                enddo
                !upper
                do j = lj1+1, ljmax
                   do i = 1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USE_SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)

                   enddo
                enddo
             enddo
             !top
             do k = 1, 1
                do j = 1, ljmax
                   do i = 1, limax
                      bc_fac     = 1.0

                      if ( bc_seaspec ) then
                         if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                         if ( .not. USE_SEASALT )  bc_fac = 0.0   ! not wanted!
                      end if

                      xn_adv(iem,i,j,k) =   xn_adv(iem,i,j,k) +&
                           bc_fac * &  ! used for sea-salt species 
                           bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_adv(ibc,iem)

                   enddo
                enddo
             enddo

          end do !n

          !/- Non-advected background species
          do n = 1,bc_used_bgn(ibc)
             iem = spc_used_bgn(ibc,n)

             do k = 2, KMAX_MID
                do j = lj0, lj1
                   !left
                   do i = 1, li0-1
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                   enddo
                   !right
                   do i = li1+1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                   enddo
                enddo
                !lower
                do j = 1, lj0-1
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                   enddo
                enddo
                !upper
                do j = lj1+1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                   enddo
                enddo
             enddo
             !top
             do k = 1, 1
                do j = 1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) &
                           +  bc_data(i_fdom(i)-IRUNBEG+1,j_fdom(j)-JRUNBEG+1,k)*bc2xn_bgn(ibc,iem)
                   enddo
                enddo
             enddo
          enddo
       endif
    enddo  ! ibc
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
          enddo
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
          enddo!n
       enddo!ibc
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
                   enddo
                   !right
                   do i = li1+1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   enddo
                enddo
                !lower
                do j = 1, lj0-1
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   enddo
                enddo
                !upper
                do j = lj1+1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   enddo
                enddo
             enddo
             !top
             do k = 1, 1
                do j = 1, ljmax
                   do i = 1, limax
                      xn_bgn(iem,i,j,k) = xn_bgn(iem,i,j,k) +misc_bc(ibc,k)
                   enddo
                enddo
             enddo
          enddo
          !/- Advected misc species
          do n = 1,bc_used_adv(ibc)
             iem = spc_used_adv(ibc,n)
             do k = 2, KMAX_MID
                do j = lj0, lj1
                   !left
                   do i = 1, li0-1
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   enddo
                   !right
                   do i = li1+1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   enddo
                enddo
                !lower
                do j = 1, lj0-1
                   do i = 1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   enddo
                enddo
                !upper
                do j = lj1+1, ljmax
                   do i = 1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   enddo
                enddo
             enddo
             !top
             do k = 1, 1
                do j = 1, ljmax
                   do i = 1, limax
                      xn_adv(iem,i,j,k) =  xn_adv(iem,i,j,k) + misc_bc(ibc,k)! 
                   enddo
                enddo
             enddo
          enddo!n
      enddo!ibc
    endif


    if (DEBUG%BCS.and.debug_proc.and.i_test>0) then
       i = i_test
       j = j_test
       print "(a20,3i4,2f8.2)","DEBUG BCS Rorvik", me, i,j,glon(i,j),glat(i,j)
       print "(a20,3i4)","DEBUG BCS Rorvik DIMS",num_adv_changed,iglobact,jglobact
       do k = 1, KMAX_MID
          print "(a20,i4,f8.2)","DEBUG O3  Debug-site ", k, &
               xn_adv(IXADV_O3,i_test,j_test,k)/PPB
       enddo
    endif ! DEBUG

    if (DEBUG%BCS.and.debug_proc) then
       itest = 1
       print *,"BoundaryConditions: No CALLS TO BOUND Cs", first_call,idebug
       !/** the following uses hard-coded  IXADV_ values for testing.
       !    Remove later **/
       info = 1   ! index for ozone in bcs
       print *,"BCs: bc2xn(info,itest) : ", bc2xn_adv(info,itest)


       info = 43   ! index for NO in bcs
       print *,"BCs: NSPECS: BC, ADV, BG, ", NTOT_BC, NSPEC_ADV, NSPEC_BGN
       print *,"BCs: Number  of bc_used: ", sum(bc_used)
       print *,"BCs: limax, ljmax",  limax, ljmax

       if (NSPEC_BGN>0) then
          do k = KMAX_MID, 1, -1
             print "(a23,i3,e14.4)","BCs NO :",k,xn_bgn(itest,i_test,j_test,k)/PPB
          enddo
       else
          print "(a)","No SET BACKGROUND BCs"
       endif
    endif !  DEBUG

    deallocate(bc_data,stat=alloc_err)
    call CheckStop(alloc_err,"de-alloc1 failed in BoundaryConditions_ml")
!    if (num_adv_changed>0) then
!       deallocate(bc_adv,stat=alloc_err)
!       call CheckStop(alloc_err,"de-alloc2 failed in BoundaryConditions_ml")
!    endif
!    if (num_bgn_changed>0) then
!       deallocate(bc_bgn,stat=alloc_err)
!       call CheckStop(alloc_err,"de-alloc3 failed in BoundaryConditions_ml")
!    endif

    if (first_call) first_call = .false.

  end subroutine BoundaryConditions

subroutine My_bcmap(iyr_trend)
! ---------------------------------------------------------------------------
! sets bc2xn_adv, bc2xn_bc, and  misc_bc
! ---------------------------------------------------------------------------
  integer, intent(in) :: iyr_trend !ds Year for which BCs are wanted
  real :: trend_ch4
  integer :: ii,i,k
  real :: decrease_factor(NGLOB_BC+1:NTOT_BC) ! Decrease factor for misc bc's
        ! Gives the factor for how much of the top-layer conc. that is left
        ! at bottom layer
  character(len=80) :: txt

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

  ! set values of 1625 in 1980, 1780 in 1990, 1820 in 2000, and 1970 in
  ! 2010. Interpolate
  ! between these for other years. Values from EMEP Rep 3/97, Table 6.2 for
  ! 1980, 1990, and from CDIAC (Mace Head) data for 2000.
  ! 2010 also from Mace Head

  if( iyr_trend >= 2010) then
    top_misc_bc(IBC_CH4) =  1870.0
  else if ( iyr_trend >= 2000) then
    top_misc_bc(IBC_CH4) = 1820 + (iyr_trend-2000)*0.1*(1870-1820) 
  else if ( iyr_trend >= 1990 ) then
    top_misc_bc(IBC_CH4) = 1780.0 + (iyr_trend-1990)*0.1*(1820-1780.0)
  else
    top_misc_bc(IBC_CH4) = 1780.0 * exp(-0.01*0.91*(1990-iyr_trend)) ! Zander,1975-1990
                                 !exp(-0.01*0.6633*(1975-iyr_trend)) ! Zander,1951-1975
  endif


  ! Reset with namelist values if set
  if ( BGND_CH4 > 0 ) then
     top_misc_bc(IBC_CH4) = BGND_CH4
  end if

  trend_ch4 = top_misc_bc(IBC_CH4)/1780.0

  if (MasterProc) then
     write(txt,"(a,2i6,f8.1,f8.3)") "BC: CH4 settings (iyr,nml,ch4,trend): ",&
             iyr_trend, nint(BGND_CH4), top_misc_bc(IBC_CH4),trend_ch4
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
    enddo
  enddo

  bc2xn_adv(IBC_H2,  IXADV_H2)    = 1.0
  bc2xn_adv(IBC_CH4, IXADV_CH4)   = 1.0

  !/- completeness check
  if (DEBUG_MYBC ) then
    print *, "In My_bcmap, NGLOB_BC, NTOT_BC is", NGLOB_BC, NTOT_BC
    do i = NGLOB_BC+1 , NTOT_BC
      print *, "In My_bcmap, sum-adv", i, " is", sum(bc2xn_adv(i,:))
      print *, "In My_bcmap, sum-bgn", i, " is", sum(bc2xn_bgn(i,:))
    enddo
  endif ! DEBUG

  do i = NGLOB_BC+1 , NTOT_BC
    call CheckStop(sum(bc2xn_adv(i,:))+sum(bc2xn_bgn(i,:))/=1.0,&
      "BCproblem - My_bcmap")
  enddo

  !/- mappings for species from Logan + obs model given with IBC index.
  include 'CM_BoundaryConditions.inc'
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
    enddo
    do iem = 1, NSPEC_BGN
      if(bc2xn_bgn(ibc,iem)>0) bc_used_bgn(ibc) = bc_used_bgn(ibc)+1
    enddo
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
    endif
  enddo ! iem

  do iem = 1, NSPEC_BGN
    if (any(bc2xn_bgn(:,iem)>0)) then
      xn_bgn_changed(iem) = .true.
      num_bgn_changed = num_bgn_changed + 1
    endif
  enddo ! iem

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
    endif
  enddo
  i = 0
  spc_bgn2changed = 0
  do iem = 1, NSPEC_BGN
    if(xn_bgn_changed(iem))then
      i = i+1
      spc_changed2bgn(i) = iem
      spc_bgn2changed(iem) = i
    endif
  enddo

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
        endif
      enddo

      ! - set bc_bgn: background (prescribed) species
      i = 0
      do iem = 1, NSPEC_BGN
        if ( bc2xn_bgn(ibc,iem) > 0.0 ) then
          i = i+1
          spc_used_bgn(ibc,i) = iem
        endif
      enddo
    endif    ! bc_used
  end do  ! ibc
end subroutine Set_bcmap

subroutine MiscBoundaryConditions(iglobact,jglobact,bc_adv,bc_bgn)
! ---------------------------------------------------------------------------
! Set bc_adv, bc_bgn
! Note: This subroutine is only used for species which have constant
! mixing ratios to start with - here CH4, H2.
! More complex variations (e.g. vertical gradients) could be also set here.
! ---------------------------------------------------------------------------
  integer, intent(in) :: iglobact,jglobact
  real, intent(inout) :: bc_adv(num_adv_changed,iglobact,jglobact,KMAX_MID), &
                         bc_bgn(num_bgn_changed,iglobact,jglobact,KMAX_MID)

  integer :: ibc, iem, i, iem1, k ! local loop variables
  integer :: itest                ! Used to specify species index

       do ibc = NGLOB_BC+1, NTOT_BC
        do i = 1,bc_used_adv(ibc)
          iem = spc_used_adv(ibc,i)
          iem1 = spc_adv2changed(iem)
          if(me==0)write(*,*)'bc_adv misc ',ibc,i,iem1
          enddo
          enddo
 if (NTOT_BC>NGLOB_BC) then
    do k=1,KMAX_MID
      do ibc = NGLOB_BC+1, NTOT_BC
        do i = 1,bc_used_adv(ibc)
          iem = spc_used_adv(ibc,i)
          iem1 = spc_adv2changed(iem)
!          bc_adv(iem1,:,:,k) = misc_bc(ibc,k)
        enddo
        do i = 1,bc_used_bgn(ibc)
          iem = spc_used_bgn(ibc,i)
          iem1 = spc_bgn2changed(iem)
!          bc_bgn(iem1,:,:,k) = misc_bc(ibc,k)
        enddo
     enddo
    enddo
  endif

  itest = 1
  if (DEBUG%BCS.and.debug_proc) write(*,*) "(a50,i4,/,(5es12.4))", &
    "From MiscBoundaryConditions: ITEST (ppb): ",&
    itest, ((bc_adv(spc_adv2changed(itest),1,1,k)/1.0e-9),k=1,20)
end subroutine MiscBoundaryConditions

subroutine Set_BoundaryConditions(mode,iglobact,jglobact,bc_adv,bc_bgn)
! ---------------------------------------------------------------------------
! Assign the global values to the interior model domain only first time
! (mode=3d) or only the edge and top boundaries are reset for other
! calls (mode=lateral). The emep model concentrations (xn) are only
! changed for those species where bcs are available, as given in the
! xn_adv_changed and xn_bgn_changed arrays.
!
! Programming note:  we use F95 "forall" and logical masks to say which grid
! squares are to be assigned.  This is probably slower than the more-explicit
! do-loops used previously, but is neater and shouldn't make too much
! difference while this reset is done only once per month.
! ---------------------------------------------------------------------------
  character(len=*), intent(in) :: mode            ! "3d" or "lateral"
  integer,          intent(in) :: iglobact,jglobact
  real, intent(in) :: bc_adv(num_adv_changed,iglobact,jglobact,KMAX_MID), &
                      bc_bgn(num_bgn_changed,iglobact,jglobact,KMAX_MID)
  logical, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: mask
  integer :: i, j, k, n, nadv, ntot
  real    :: bc_fac      ! Set to 1.0, except sea-salt over land = 0.01
  logical :: bc_seaspec  ! if sea-salt species
  character(len=20) :: txtout 

   if (mode=="3d") then
    mask(:,:,:) = .true.    ! We set everything

   elseif (mode=="lateral") then
    mask(:,:,:)  = .false.        ! Initial

    ! Set edges (except on the top)
    ! there may be no neighbor, but no external boundary (Poles in lat lon)
    if(neighbor(SOUTH)==NOPROC) mask(:,1:(lj0-1),2:KMAX_MID)     = .true.
    if(neighbor(NORTH)==NOPROC) mask(:,(lj1+1):ljmax,2:KMAX_MID) = .true.
    if(neighbor(EAST) ==NOPROC) mask((li1+1):limax,:,2:KMAX_MID) = .true.
    if(neighbor(WEST) ==NOPROC) mask(1:(li0-1),:,2:KMAX_MID)     = .true.

    mask(:,:,1) = .true.        !Set top layer
  else
    call CheckStop("BCs:Illegal option failed in BoundaryConditions_ml")
  endif

  !/- Set concentrations (xn) from boundary conditions (bcs)

  ! Note on domains: although the geographical stuff has been specified
  ! for the whole MAXLIMAX,MAXLJMAX grid, the interpolations take time
  ! and are only needed for the sub-domain actually used, i.e. for
  ! limax, ljmax.

  !/- Advected species. Sea-salt is special as we only want BICs over sea-areas
  !forall(i=1:limax, j=1:ljmax, k=1:KMAX_MID, n=1:num_adv_changed, mask(i,j,k))
  do n = 1, num_adv_changed
    nadv = spc_changed2adv(n)
    ntot = nadv + NSPEC_SHL 

    bc_seaspec = .false.
    if ( USE_SEASALT .and. ( index( species(ntot)%name, "SEASALT_" ) > 0 ) ) then
      bc_seaspec = .true.
    end if
    if ( debug_proc ) write (*,*) "SEAINDEX", &
           trim(species(ntot)%name), n, ntot, bc_seaspec

    do k = 1, KMAX_MID
      do j = 1, ljmax
        do i = 1, limax
          if ( mask(i,j,k) ) then

            bc_fac     = 1.0
            if ( bc_seaspec ) then
                  if ( .not. mainly_sea(i,j))  bc_fac = 0.001 ! low over land
                  if ( .not. USE_SEASALT )  bc_fac = 0.0   ! not wanted!
            end if

            xn_adv(nadv,i,j,k) =   &
               bc_fac * &  ! used for sea-salt species 
                 bc_adv(n,(i_fdom(i)-IRUNBEG+1),(j_fdom(j)-JRUNBEG+1),k)
          end if !mask
        end do ! i
      end do ! j
    end do ! k
   !endforall
    if ( DEBUG%BCS .and. debug_proc ) then
     i=debug_li
     j=debug_lj
     k=KMAX_MID
     txtout = "BCSET:" // trim(species(ntot)%name)
     call datewrite( trim(txtout), n, (/ xn_adv(nadv,i,j,k),  &
           bc_adv(n,(i_fdom(i)-IRUNBEG+1),(j_fdom(j)-JRUNBEG+1),k) /) )
    end if ! DEBUG
 end do !n
 

  !/- Non-advected background species
!  forall(i=1:limax, j=1:ljmax, k=1:KMAX_MID, n=1:num_bgn_changed)
!    xn_bgn(spc_changed2bgn(n),i,j,k) = &
!        bc_bgn(n,(i_fdom(i)-IRUNBEG+1),(j_fdom(j)-JRUNBEG+1),k)
!  endforall
end subroutine Set_BoundaryConditions    ! call every 3-hours

end module BoundaryConditions_ml
