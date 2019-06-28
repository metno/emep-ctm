! <MassBudget_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module   MassBudget_mod
! ----------------------------------------------------------------------------
! DESCRIPTION
! Routine to cross check the mass balance of the model
!_____________________________________________________________________________
use CheckStop_mod,    only: CheckStop
use ChemDims_mod,     only: NSPEC_ADV, NSPEC_SHL
use ChemSpecs_mod,    only: species_adv
use Chemfields_mod,   only: xn_adv        ! advected species
use Config_module,    only: KMAX_MID,KCHEMTOP,& ! Start and upper k for 1d fields
                              MasterProc,       & ! Master processor
                              dt_advec,         & ! time-step
                              PT,               & ! Pressure at top
                              USE_OCEAN_NH3,USE_OCEAN_DMS,FOUND_OCEAN_DMS,&
                              EXTENDEDMASSBUDGET
use Debug_module,     only: DEBUG_MASS
use EmisDef_mod,      only: O_NH3, O_DMS
use GridValues_mod,   only: xmd, &  
                              gridwidth_m,dA,dB,debug_proc,debug_li,debug_lj
use Io_mod,           only: IO_LOG, PrintLog, datewrite
use MetFields_mod,    only: ps            ! surface pressure
! do not use "only", because MPI_IN_PLACE does not behave well on certain
! versions of gfortran(?), and MPI stuff clearly inidcated anyway
use MPI_Groups_mod !   , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_REAL8, &
                  !            MPI_INTEGER, MPI_LOGICAL, &
                  !            MPI_MIN, MPI_MAX, MPI_SUM, MPI_IN_PLACE, &
                  !            MPI_COMM_CALC, MPI_COMM_WORLD, MPISTATUS,&
                  ! IERROR, ME_MPI, NPROC_MPI
use Par_mod,          only: &
  li0,li1,& ! First/Last local index in long. when outer boundary is excluded
  lj0,lj1   ! First/Last local index in lat.  when outer boundary is excluded
use PhysicalConstants_mod,only: GRAV,ATWAIR! Mol. weight of air(Jones,1992)
use SmallUtils_mod,   only: find_index
use ZchemData_mod,    only: M, rcemis ! Air concentrations , emissions
implicit none
private

! Some work arrays used in Aqueous_mod and (in future) DryDry:
! Use ADV index, as Dry/WetDep makes no seance for SHL.
real, public, save, dimension(NSPEC_ADV) ::   &
wdeploss=0.0, ddeploss=0.0

! The following parameters are used to check the global mass budget:
! Initialise here also.
real, public, save, dimension(NSPEC_ADV) ::   &
  sumini   = 0.0,  & !  initial mass
  fluxin   = 0.0,  & !  mass in  across lateral boundaries
  fluxout  = 0.0,  & !  mass out across lateral boundaries
  fluxin_top  = 0.0,  & !  mass in  across top
  fluxout_top = 0.0,  & !  mass out across top
  totddep  = 0.0,  & !  total dry dep
  totwdep  = 0.0,  & !  total wet dep
  totem    = 0.0     !  total emissions

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
  integer i, j, k, n          ! lon,lat,lev indexes, n - No. of species
  real rwork,fac,wgt_fac

  fac = GRIDWIDTH_M*GRIDWIDTH_M/GRAV
  do k=2,KMAX_MID
    do j=lj0,lj1
      do i=li0,li1
        rwork = fac*(dA(k)+dB(k)*ps(i,j,1))* xmd(i,j)
        sumini(:) = sumini(:) + xn_adv(:,i,j,k)*rwork  ! sumini in kg
      end do
    end do
  end do

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, sumini , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)

  if(MasterProc.and.EXTENDEDMASSBUDGET)then
    do n = 1,NSPEC_ADV
      if(sumini(n)<=0.) cycle
      wgt_fac=species_adv(n)%molwt/ATWAIR
      write(IO_LOG,"(a15,i4,4x,e10.3)") "Initial mass",n,sumini(n)*wgt_fac
      write(*,"(a15,i4,4x,e10.3)") "Initial mass",n,sumini(n)*wgt_fac
    end do
  end if

 end subroutine Init_massbudget
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
    scaling_k = scaling * (dA(k) + dB(k)*ps(i,j,1))/M(k)
    if(all((/DEBUG_MASS,debug_proc,i==debug_li,j==debug_lj/)))&
      call datewrite("MASSRC ",k,(/dB(k)*ps(i,j,1),xmd(i,j),&
          ps(i,j,1),scaling_k/))

    do iadv = 1, NSPEC_ADV
      itot = iadv + NSPEC_SHL
      totem(iadv) = totem(iadv) + rcemis( itot, k ) * scaling_k
    end do
  end do ! k loop

end subroutine emis_massbudget_1d
!----------------------------------------------------------------------------
subroutine massbudget()
! sums over all sulphur and nitrogen, so is model independant.

  integer ::  i, j, k, n            ! lon,lat,lev indexes, n - No. of species
  integer :: ix_o3, ifam            ! family index
  integer :: iomb                   ! for MassBugetSummary.txt Table
  real,  dimension(NSPEC_ADV,KMAX_MID) ::  sumk   ! total mass in each layer
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
    frac_mass,          & ! mass budget frac. (should=1) for groups of species
    gfluxin,gfluxout,   & ! flux in  and out
    gtotem,             & ! total emission
    gtotddep, gtotwdep, & ! total dry and wet deposition
    natoms                ! number of S, N or C atoms

  real :: totdiv,helsum,fac,o3_fac,wgt_fac


  fac=GRIDWIDTH_M*GRIDWIDTH_M/GRAV

  sum_mass(:)   = 0.0
  frac_mass(:)  = 0.0
  xmax(:)       =-2.0
  xmin (:)      = 2.0
  gfluxin(:)    = fluxin(:)+fluxin_top(:)
  gfluxout(:)   = fluxout(:)+fluxout_top(:)
  gtotem(:)     = totem(:)
  gtotddep(:)   = totddep(:)
  gtotwdep(:)   = totwdep(:)
  sumk(:,:)     = 0.0

  do k = 1,KMAX_MID
    do j = lj0,lj1
      do i = li0,li1
        helsum  = fac*(dA(k)+dB(k)*ps(i,j,1))* xmd(i,j)
        xmax(:) = amax1(xmax(:),xn_adv(:,i,j,k))
        xmin(:) = amin1(xmin(:),xn_adv(:,i,j,k))
        sumk(:,k) = sumk(:,k) + xn_adv(:,i,j,k)*helsum

        if(all((/DEBUG_MASS,debug_proc,i==debug_li,j==debug_lj/)))&
          call datewrite("MASSBUD",k,(/(dA(k)*dB(k)*ps(i,j,1))*xmd(i,j)/&
          GRAV*GRIDWIDTH_M*GRIDWIDTH_M,ps(i,j,1),PT,xmd(i,j)/))
      end do
    end do
  end do

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, xmax, NSPEC_ADV,&
    MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, xmin   , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, gfluxin , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, gfluxout , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, fluxin_top , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, fluxout_top , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, fluxin , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, fluxout , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, gtotem , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, gtotddep , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, gtotwdep , NSPEC_ADV, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, sumk , NSPEC_ADV*KMAX_MID, &
    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)

!     make some temporary variables used to hold the sum over all
!     domains. Remember that sumini already holds the sum over all
!     domains, see inass
!
  amax(:) = max( amax(:), xmax(:) )
  amin(:) = min( amin(:), xmin(:) )
  do k = 2,KMAX_MID
    sum_mass(:) = sum_mass(:)+sumk(:,k)
  end do


!O3 flux are also printed out
  if(MasterProc)then
     ix_o3=find_index( 'O3', species_adv(:)%name )
     if(ix_o3>0)then
        o3_fac=species_adv(ix_o3)%molwt/ATWAIR
        fluxin_top(ix_o3)=fluxin_top(ix_o3)*o3_fac
        fluxout_top(ix_o3)=fluxout_top(ix_o3)*o3_fac
        fluxin(ix_o3)=fluxin(ix_o3)*o3_fac
        fluxout(ix_o3)=fluxout(ix_o3)*o3_fac
        33 FORMAT(5(A,es10.3))
        write(*,*)'Ozone fluxes (kg):'   
        write(*,33)'Net from top = ', fluxin_top(ix_o3)-fluxout_top(ix_o3),&
          '  in from top = ',fluxin_top(ix_o3),' out of top = ',fluxout_top(ix_o3)
        write(*,33)'Net lateral faces = ', fluxin(ix_o3)-fluxout(ix_o3),&
          '  in lateral faces = '  ,fluxin(ix_o3), &
          ' out lateral faces = ',  fluxout (ix_o3)
        write(*,33)'O3 in atmosphere at start of run = ', &
          sumini(ix_o3)*o3_fac, ' at end of run = ', sum_mass(ix_o3)*o3_fac
        write(*,33)'O3 dry deposited = ',&
             gtotddep(ix_o3)*species_adv(ix_o3)%molwt
     else
        write(*,*)'O3 index not found'
     end if
  end if


  do n = 1,NSPEC_ADV
    totdiv = sumini(n) + gtotem(n) + gfluxin(n)
    frac_mass(n) = sum_mass(n) + (gtotddep(n)+gtotwdep(n))*ATWAIR + gfluxout(n)
    if(totdiv>0.0) frac_mass(n) = frac_mass(n)/totdiv
  end do


  if(MasterProc) then   ! printout from node 0
    if(EXTENDEDMASSBUDGET)then
      do n=1,NSPEC_ADV
        wgt_fac=species_adv(n)%molwt/ATWAIR
        if(gtotem(n)>0.0) write(*,*)'tot. emission of '//&
               trim(species_adv(n)%name)//' ',gtotem(n)*wgt_fac
      end do
    end if

    call PrintLog('++++++++++++++++++++++++++++++++++++++++++++++++')
    do ifam = 1, 3
      write(logtxt,"(a,i3,a12)") 'Mass balance ', ifam, family_name(ifam)
      call PrintLog(logtxt)
      select case(ifam)
        case(1);natoms = real(species_adv(:)%sulphurs)
        case(2);natoms = real(species_adv(:)%nitrogens)
        case(3);natoms = real(species_adv(:)%carbons)
      end select

      family_init(ifam)   = dot_product(sumini(:)  ,natoms(:))
      family_mass(ifam)   = dot_product(sum_mass(:),natoms(:))
      family_inflow(ifam) = dot_product(gfluxin(:) ,natoms(:))
      family_outflow(ifam)= dot_product(gfluxout(:),natoms(:))
      family_ddep(ifam)   = dot_product(gtotddep(:),natoms(:))
      family_wdep(ifam)   = dot_product(gtotwdep(:),natoms(:))
      family_em(ifam)     = dot_product(gtotem(:)  ,natoms(:))

    !convert into kg
      select case(ifam)
        case(1);wgt_fac=32/ATWAIR!sulphurs
        case(2);wgt_fac=14/ATWAIR!nitrogens
        case(3);wgt_fac=12/ATWAIR!carbons
      end select
      family_init(ifam)=family_init(ifam)*wgt_fac
      family_inflow(ifam)=family_inflow(ifam)*wgt_fac
      family_em(ifam)=family_em(ifam)*wgt_fac
      family_mass(ifam)=family_mass(ifam)*wgt_fac
      family_outflow(ifam)=family_outflow(ifam)*wgt_fac
      family_ddep(ifam)=family_ddep(ifam)*wgt_fac*ATWAIR
      family_wdep(ifam)=family_wdep(ifam)*wgt_fac*ATWAIR

      family_input(ifam) = family_init(ifam)    &
                         + family_inflow(ifam)  &
                         + family_em(ifam)

      if(family_input(ifam)>0.0) &
        family_fracmass(ifam) = (family_mass(ifam)         &
                              +  family_outflow(ifam)      &
                              +  family_ddep(ifam)  &
                              +  family_wdep(ifam)) &
                              / family_input(ifam)


      call PrintLog('++++++++++++++++++++++++++++++++++++++++++++++++')
      write(logtxt,"(a9,5a12)")" ","sumini","summas","fluxout","fluxin",&
        "fracmass"
      call PrintLog(logtxt)

      write(logtxt,"(a9,5es12.4)") family_name(ifam), &
        family_init(ifam), family_mass(ifam), family_outflow(ifam), &
        family_inflow(ifam), family_fracmass(ifam)
      call PrintLog(logtxt)

      write(logtxt,"(a9,3a14)")"ifam","totddep","totwdep","totem"
      call PrintLog(logtxt)
      write(logtxt,"(i9,3es14.3)") ifam, &
        family_ddep(ifam), family_wdep(ifam), family_em(ifam)
      call PrintLog(logtxt)
      call PrintLog('++++++++++++++++++++++++++++++++++++++++++++++++')
    end do  ! ifam = 1,3
  end if

  if(MasterProc) then     ! printout from node 0

    if( EXTENDEDMASSBUDGET) then     ! printout from node 0
      !/.. now use species array which is set in My_MassBudget_mod
      do n=1,NSPEC_ADV
        wgt_fac=species_adv(n)%molwt/ATWAIR
        write(IO_LOG,*)
        write(*,*)
        do k=1,KMAX_MID
          write(IO_LOG,"(' Spec ',i3,2x,a12,5x,'k= ',i2,5x,es12.5)")&
            n,species_adv(n)%name, k,sumk(n,k)*wgt_fac
          write(*     ,"(' Spec ',i3,2x,a12,5x,'k= ',i2,5x,es12.5)")&
            n,species_adv(n)%name, k,sumk(n,k)*wgt_fac
        end do
      end do
     end if ! EXTENDED

    ! SUMMARY TABLE:
     open(newunit=iomb,file="MassBudgetSummary.txt")
     write(iomb,'(a)') ' # Mass Budget. Units are kg with MW used here. ADJUST! if needed' 
     write(iomb,'(a3,1x,a14,a7,99a12)') '#n', 'Spec       ', &
       'usedMW', 'emis', 'ddep', 'wdep', 'init', 'sum_mass', &
       'fluxout', 'fluxin', 'frac_mass'

     do n=1,NSPEC_ADV
        wgt_fac=species_adv(n)%molwt/ATWAIR
        if( EXTENDEDMASSBUDGET) then     ! printout from node 0
           !/.. now use species array which is set in My_MassBudget_mod
             write(*,*)
             write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
             write(*,*)
    
             write(*,"(a3,6a12)")" n ", "Spec", &
               "sumini", "summas", "fluxout", "fluxin", "fracmass"
             write(*,"(i3,1x,a11,5es12.4)") n,species_adv(n)%name, &
               sumini(n)*wgt_fac, sum_mass(n)*wgt_fac, gfluxout(n)*wgt_fac, &
                 gfluxin(n)*wgt_fac, frac_mass(n)
             write(*,*)
             write(*,"(a3,6a12)") "n ", "species", "totddep", "totwdep", "totem"
             write(*,"(i3,1x,a11,5es12.4)") n, species_adv(n)%name, &
                gtotddep(n)*wgt_fac*ATWAIR, gtotwdep(n)*wgt_fac*ATWAIR, &
                gtotem(n)*wgt_fac
             write(*,*)
             write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        end if ! EXTENDED

        write(iomb,'(i3,1x,a14,f7.1,99es12.4)') n,adjustl(species_adv(n)%name), &
          species_adv(n)%molwt, &  ! REMEMBER. Sometimes a dummy value, e.g. 1.0
          gtotem(n)*wgt_fac, gtotddep(n)*wgt_fac*ATWAIR, &
          gtotwdep(n)*wgt_fac*ATWAIR, sumini(n)*wgt_fac, sum_mass(n)*wgt_fac,&
          gfluxout(n)*wgt_fac, gfluxin(n)*wgt_fac, frac_mass(n)

     end do
     close(iomb)
  end if  ! MasterProc

  if(FOUND_OCEAN_DMS)then
     ! update dms budgets
     CALL MPI_ALLREDUCE(MPI_IN_PLACE, O_DMS%sum_month, 1,&
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
     CALL MPI_ALLREDUCE(MPI_IN_PLACE, O_NH3%sum_month, 1,&
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR)
     !natso2 already allreduced
     
     if(MasterProc)then
        write(*,59)'SO2 from ocean DMS cdf file last month ',O_DMS%sum_month
        
        O_DMS%sum_year=O_DMS%sum_year+O_DMS%sum_month
        O_DMS%sum_month=0.0
        
59      format(A,6F14.5)
        write(*,59)'SO2 from ocean DMS cdf file ',O_DMS%sum_year
     end if
  end if
  if(MasterProc)then
     if(USE_OCEAN_NH3)then
        write(*,59)'NH3 emisions from ocean cdf file (Gg)',O_NH3%sum_year
     end if
  end if
end subroutine massbudget
!--------------------------------------------------------------------------
 end module MassBudget_mod
!--------------------------------------------------------------------------
