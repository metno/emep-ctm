! <MassBudget_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! ----------------------------------------------------------------------------
module   MassBudget_ml
! ----------------------------------------------------------------------------
! DESCRIPTION
! Routine to cross check the mass balance of the model
!_____________________________________________________________________________
use CheckStop_ml,       only: CheckStop
use ChemChemicals_ml,   only: species_adv   ! species identifier (advected only)
use ChemSpecs_adv_ml,   only: NSPEC_ADV     ! No. species (long-lived)
use ChemSpecs_shl_ml,   only: NSPEC_SHL     ! No. species (shorshort-lived)
use Chemfields_ml,      only: xn_adv        ! advected species
use GridValues_ml,      only: carea,xmd, &  ! cell area, 1/xm2 where xm2 is
                                  ! the area factor in the middle of the cell
                              gridwidth_m,dA,dB,debug_proc,debug_li,debug_lj
use Io_ml,              only: IO_RES, PrintLog, datewrite
use MetFields_ml,       only: ps            ! surface pressure
use ModelConstants_ml,  only: KMAX_MID,KCHEMTOP,& ! Start and upper k for 1d fields
                              MasterProc,       & ! Master processor
                              dt_advec,         & ! time-step
                              PT,               & ! Pressure at top
                              ATWAIR,           & ! Mol. weight of air(Jones,1992)
                              DEBUG_MASS,EXTENDEDMASSBUDGET
use Par_ml,             only: &
  li0,li1,& ! First/Last local index in longitude when outer boundary is excluded
  lj0,lj1   ! First/Last local index in latitude  when outer boundary is excluded
use PhysicalConstants_ml,only: GRAV
use Setup_1dfields_ml,  only: amk, rcemis ! Air concentrations , emissions

implicit none
private
INCLUDE 'mpif.h'
INTEGER STATUS(MPI_STATUS_SIZE),INFO
real    MPIbuff(NSPEC_ADV*KMAX_MID)

! Some work arrays used in Aqueous_ml and (in future) DryDry:
! Use ADV index, as Dry/WetDep makes no seance for SHL.
real, public, save, dimension(NSPEC_ADV) ::   &
wdeploss=0.0, ddeploss=0.0

! The following parameters are used to check the global mass budget:
! Initialise here also.
real, public, save, dimension(NSPEC_ADV) ::   &
  sumint   = 0.0,  & !  initial mass
  fluxin   = 0.0,  & !  mass in  across lateral boundaries
  fluxout  = 0.0,  & !  mass out across lateral boundaries
  totddep  = 0.0,  & !  total dry dep
  totwdep  = 0.0,  & !  total wet dep
  totem    = 0.0,  & !  total emissions
  totox    = 0.0,  & !  total oxidation
  totldep  = 0.0     !  local deposition (not in use - Lagrangian)

real, public, save, dimension(NSPEC_ADV) ::  &
  amax = -2.0,  &  ! maximum concentration in field -2
  amin =  2.0      ! minimum concentration in field  2

public :: Init_massbudget
public :: massbudget
public :: emis_massbudget_1d
!public :: DryDep_Budget

contains

!----------------------------------------------------------------------------
subroutine Init_massbudget()
! Initialise mass-budget - calculate mass of concentrations fields
! within 3-D grid, after boundary conditions
!
!----------------------------------------------------------------------
  integer i, j, k, n, info    ! lon,lat,lev indexes
                              ! n - No. of species
                              ! info - printing info
  real rwork

  do k=2,KMAX_MID
    do j=lj0,lj1
      do i=li0,li1
        rwork = carea(k)* xmd(i,j)*(ps(i,j,1) - PT)
        sumint(:) = sumint(:) + xn_adv(:,i,j,k)*rwork  ! sumint in kg
      enddo
    enddo
  enddo

  MPIbuff(1:NSPEC_ADV)= sumint (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, sumint , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)

  if(MasterProc.and.EXTENDEDMASSBUDGET)then
    do n = 1,NSPEC_ADV
      if(sumint(n)<=0.) cycle
      write(IO_RES,"(a15,i4,4x,e10.3)") "Initial mass",n,sumint(n)
      write(*,"(a15,i4,4x,e10.3)") "Initial mass",n,sumint(n)
    enddo
  endif

 endsubroutine Init_massbudget
!----------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
subroutine emis_massbudget_1d(i,j)
  integer, intent(in) :: i,j    ! coordinates of column
  integer    :: k, iadv, itot   ! loop variables
  real :: scaling, scaling_k

  !Mass Budget calculations
  !   Adding up the emissions in each timestep

  !Do not include values on outer frame
  if(i<li0.or.i>li1.or.j<lj0.or.j>lj1)return

  scaling = dt_advec * xmd(i,j)* gridwidth_m*gridwidth_m / GRAV

  do k = KCHEMTOP,KMAX_MID
    scaling_k = scaling * (dA(k) + dB(k)*ps(i,j,1))/amk(k)
    if(all((/DEBUG_MASS,debug_proc,i==debug_li,j==debug_lj/)))&
      call datewrite("MASSRC ",k,(/dB(k)*ps(i,j,1),xmd(i,j),ps(i,j,1),scaling_k/))

    do iadv = 1, NSPEC_ADV
      itot = iadv + NSPEC_SHL
      totem(iadv) = totem(iadv) + rcemis( itot, k ) * scaling_k
    enddo
  enddo ! k loop

endsubroutine emis_massbudget_1d
!----------------------------------------------------------------------------
subroutine massbudget()
! sums over all sulphur and nitrogen, so is model independant.

  integer ::  i, j, k, n, nn, info  ! lon,lat,lev indexes
                                    ! n - No. of species
                                    ! nn - Total no. of short lived and advected species
                                    ! info - printing info
  integer :: ifam                                 ! family index
  real, dimension(NSPEC_ADV,KMAX_MID) ::  sumk    ! total mass in each layer
  integer, parameter :: NFAMILIES = 3             ! No. of families
  character(len=*), dimension(NFAMILIES), parameter :: &
    family_name = (/ "Sulphur ", "Nitrogen", "Carbon  " /)
  character(len=200) :: logtxt

  real, dimension(NFAMILIES) ::&
    family_init,    & ! initial total mass of species family
    family_mass,    & ! total family mass at the end of the model run
    family_inflow,  & ! total family mass flowing in
    family_outflow, & ! total family mass flowing out
    family_ddep,    & ! total family mass dry dep.
    family_wdep,    & ! total family mass wet dep.
    family_em,      & ! total family mass emitted
    family_input,   & ! total family mass input
    family_fracmass  ! mass fraction (should be 1.0)

  real, dimension(NSPEC_ADV) :: &
    xmax, xmin,         & ! min and max value for the individual species
    sum_mass,           & ! total mass of species
    frac_mass,          & ! mass budget fraction (should=1) for groups of species
    gfluxin,gfluxout,   & ! flux in  and out
    gtotem,             & ! total emission
    gtotddep, gtotwdep, & ! total dry and wet deposition
    gtotldep,           & ! local dry deposition
    gtotox,             & ! oxidation of SO2
    natoms                ! number of S, N or C atoms

  real :: totdiv,helsum

  sum_mass(:)   = 0.0
  frac_mass(:)  = 0.0
  xmax(:)       =-2.0
  xmin (:)      = 2.0
  gfluxin(:)    = fluxin(:)
  gfluxout(:)   = fluxout(:)
  gtotem(:)     = totem(:)
  gtotddep(:)   = totddep(:)
  gtotwdep(:)   = totwdep(:)
  gtotldep(:)   = totldep(:)
  gtotox(:)     = totox(:)
  sumk(:,:)     = 0.0

  do k = 1,KMAX_MID
    do j = lj0,lj1
      do i = li0,li1
        helsum  = carea(k)*xmd(i,j) * (ps(i,j,1) - PT)
        xmax(:) = amax1(xmax(:),xn_adv(:,i,j,k))
        xmin(:) = amin1(xmin(:),xn_adv(:,i,j,k))
        sumk(:,k) = sumk(:,k) + xn_adv(:,i,j,k)*helsum

        if(all((/DEBUG_MASS,debug_proc,i==debug_li,j==debug_lj/)))&
          call datewrite("MASSBUD",k,(/carea(k),ps(i,j,1),PT,xmd(i,j)/))
      enddo
    enddo
  enddo

  MPIbuff(1:NSPEC_ADV)= xmax(1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, xmax, NSPEC_ADV,&
    MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, INFO)
  MPIbuff(1:NSPEC_ADV)= xmin   (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, xmin   , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO)
  MPIbuff(1:NSPEC_ADV)= gfluxin (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, gfluxin , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
  MPIbuff(1:NSPEC_ADV)= gfluxout (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, gfluxout , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
  MPIbuff(1:NSPEC_ADV)= gtotem (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, gtotem , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
  MPIbuff(1:NSPEC_ADV)= gtotddep (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, gtotddep , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
  MPIbuff(1:NSPEC_ADV)= gtotwdep (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, gtotwdep , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
  MPIbuff(1:NSPEC_ADV)= gtotldep (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, gtotldep , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
  MPIbuff(1:NSPEC_ADV)= gtotox (1:NSPEC_ADV)
  CALL MPI_ALLREDUCE(MPIbuff, gtotox , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)
  j=0
  do k=1,KMAX_MID
    do i=1,NSPEC_ADV
      j=j+1
      MPIbuff(j)= sumk(i,k)
    enddo
  enddo
  CALL MPI_ALLREDUCE(MPIbuff, sumk , NSPEC_ADV*KMAX_MID, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO)

!     make some temporary variables used to hold the sum over all
!     domains. Remember that sumint already holds the sum over all
!     domains, see inass
!
  amax(:) = max( amax(:), xmax(:) )
  amin(:) = min( amin(:), xmin(:) )
  do k = 2,KMAX_MID
    sum_mass(:) = sum_mass(:)+sumk(:,k)
  enddo

  do n = 1,NSPEC_ADV
    totdiv = sumint(n) + gtotem(n) + gfluxin(n)
    frac_mass(n) = sum_mass(n) + (gtotddep(n)+gtotwdep(n))*ATWAIR + gfluxout(n)
    if(totdiv>0.0) frac_mass(n) = frac_mass(n)/totdiv
  enddo


  if(MasterProc) then   ! printout from node 0
    if(EXTENDEDMASSBUDGET)then
      do n=1,NSPEC_ADV
        if(gtotem(n)>0.0) write(*,*)'tot. emission of species ',n,gtotem(n)
      enddo
    endif

    call PrintLog('++++++++++++++++++++++++++++++++++++++++++++++++')
    do ifam = 1, 3
      write(logtxt,"(a,i3,a12)") 'Mass balance ', ifam, family_name(ifam)
      call PrintLog(logtxt)
      select case(ifam)
        case(1);natoms = real(species_adv(:)%sulphurs)
        case(2);natoms = real(species_adv(:)%nitrogens)
        case(3);natoms = real(species_adv(:)%carbons)
      endselect

      family_init(ifam)   = dot_product(sumint(:)  ,natoms(:))
      family_mass(ifam)   = dot_product(sum_mass(:),natoms(:))
      family_inflow(ifam) = dot_product(gfluxin(:) ,natoms(:))
      family_outflow(ifam)= dot_product(gfluxout(:),natoms(:))
      family_ddep(ifam)   = dot_product(gtotddep(:),natoms(:))
      family_wdep(ifam)   = dot_product(gtotwdep(:),natoms(:))
      family_em(ifam)     = dot_product(gtotem(:)  ,natoms(:))

      family_input(ifam) = family_init(ifam)    &
                         + family_inflow(ifam)  &
                         + family_em(ifam)

      if(family_input(ifam)>0.0) &
        family_fracmass(ifam) = (family_mass(ifam)         &
                              +  family_outflow(ifam)      &
                              +  family_ddep(ifam)*ATWAIR  &
                              +  family_wdep(ifam)*ATWAIR) &
                              / family_input(ifam)


      call PrintLog('++++++++++++++++++++++++++++++++++++++++++++++++')
      write(logtxt,"(a9,5a12)")" ","sumint","summas","fluxout","fluxin","fracmass"
      call PrintLog(logtxt)

      write(logtxt,"(a9,5es12.4)") family_name(ifam), &
        family_init(ifam), family_mass(ifam), family_outflow(ifam), &
        family_inflow(ifam), family_fracmass(ifam)
      call PrintLog(logtxt)

      write(logtxt,"(a9,3a14)")"ifam","totddep","totwdep","totem"
      call PrintLog(logtxt)
      write(logtxt,"(i9,3es14.3)") ifam, &
        family_ddep(ifam)*ATWAIR, family_wdep(ifam)*ATWAIR, family_em(ifam)
      call PrintLog(logtxt)
      call PrintLog('++++++++++++++++++++++++++++++++++++++++++++++++')
    enddo  ! ifam = 1,3
  endif

  if(MasterProc.and.EXTENDEDMASSBUDGET) then     ! printout from node 0
    !/.. now use species array which is set in My_MassBudget_ml
    do n=1,NSPEC_ADV
      write(IO_RES,*)
      write(*,*)
      do k=1,KMAX_MID
        write(IO_RES,"(' Spec ',i3,2x,a12,5x,'k= ',i2,5x,es12.5)")&
          n,species_adv(n)%name, k,sumk(n,k)
        write(*     ,"(' Spec ',i3,2x,a12,5x,'k= ',i2,5x,es12.5)")&
          n,species_adv(n)%name, k,sumk(n,k)
      enddo
    enddo
    do n=1,NSPEC_ADV
      write(*,*)
      write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*)
      write(*,"(a3,6a12)")" n ", "Spec", &
        "sumint", "summas", "fluxout", "fluxin", "fracmass"
      write(*,"(i3,1x,a11,5es12.4)") n,species_adv(n)%name, &
        sumint(n), sum_mass(n), gfluxout(n), gfluxin(n), frac_mass(n)
      write(*,*)
      write(*,"(a3,6a12)")  " n ", "species", &
        "totox", "totddep", "totwdep", "totem", "totldep"
      write(*,"(i3,1x,a11,5es12.4)") n, species_adv(n)%name, &
        gtotox(n), gtotddep(n), gtotwdep(n), gtotem(n), gtotldep(n)
      write(*,*)
      write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
    enddo
  endif  ! MasterProc
endsubroutine massbudget
!--------------------------------------------------------------------------
 end module MassBudget_ml
!--------------------------------------------------------------------------
