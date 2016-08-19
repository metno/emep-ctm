! <Solver.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                          module Chemsolver_ml
! MOD OD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !=======================================================================!
  ! The following chemical solver uses variable chemical timesteps and
  ! is based on the scheme suggested in J.G. Verwer and D. Simpson (1995) 
  ! "Explicit methods for stiff ODEs from atmospheric chemistry", 
  ! Aplied Numerical Mathematics 18 (1995) 413.
  !
  ! Note that the exact formula used have been re-arranged for greater
  ! efficiency (Steffen Unger).
  ! Variable Dchem is used to keep track of changes from call to call.
  ! Note: decoupling of (NO3,N2O5), (PAN,CH3COO2), (MPAN,MACRO2)
  ! variable timestep (Peter Wind)
  !=======================================================================!
 

    use Aqueous_ml,        only: aqrck, ICLOHSO2, ICLRC1, ICLRC2, ICLRC3   
    use Biogenics_ml,      only: BIO_ISOP, BIO_TERP
    use CheckStop_ml,      only: CheckStop
    use DefPhotolysis_ml         ! => IDHNO3, etc.
    use Emissions_ml,      only: KEMISTOP    
    use GenSpec_tot_ml           ! => NSPEC_TOT, O3, NO2, etc.
    use GenSpec_bgn_ml           ! => IXBGN_  indices and xn_2d_bgn values
    use GenRates_rct_ml,   only: set_night_rct, ONLY_NIGHT
    use ModelConstants_ml, only: KMAX_MID, KCHEMTOP, dt_advec,dt_advec_inv
    use My_Aerosols_ml,    only: SEASALT
    use My_Emis_ml                        ! => QRCNO, etc.
    use OrganicAerosol_ml, only: Fgas
    use Par_ml,            only: me, MAXLIMAX, MAXLJMAX  ! me for TEST
    use Setup_1dfields_ml, only: rcemis,        & ! photolysis, emissions
                                 rcbio,         & ! biogenic emis
                                 rc_Rn222,      & ! Pb210
                                 rct, rcmisc,   & ! reaction rate coeffients
                                 xn_2d,         & 
                                 rh,            & 
                                 rcss,amk         ! Sea salt emission rate
    use N2O5_hydrolysis_ml, only :VOLFACSO4,VOLFACNO3,VOLFACNH4,&
                                 f_Riemer! to weight the hydrolysis of N2O5 with NO3,SO4 mass
  implicit none

  private
  public  :: chemistry        ! Runs chemical solver

  INCLUDE 'mpif.h'

  integer::  STATUS(MPI_STATUS_SIZE),INFO
  integer, parameter:: nchemMAX=12

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  subroutine chemistry(i,j)

    !.. In
    integer, intent(in) ::  i,j       ! Coordinates (needed for Dchem)

    real, dimension(NSPEC_TOT,KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX), save :: &
                    Dchem=0.0  ! Concentration increments due to chemistry

    logical, save ::  first_call = .true.

    real, parameter ::  CPINIT = 0.0 ! 1.0e-30  ! small value for init

    !  Local
    integer, dimension(KCHEMTOP:KMAX_MID) :: toiter
    integer ::  k, ichem, iter,n    ! Loop indices
    integer, save ::  nchem         ! No chem time-steps
    real    ::  dt2 
    real    ::  P, L                ! Production, loss terms
    real    :: psd_h2o2 ! Pseudo H2O2 concentration (lower when high so2)
    real    :: xextrapol, L1,L2,P1,P2,C1,C2,DIVID    !help variable

    ! Concentrations : xold=old, x=current, xnew=predicted
    ! - dimensioned to have same size as "x" 

    real, dimension(NSPEC_TOT)      :: &
                        x, xold ,xnew   ! Working array [molecules/cm3]
    real, dimension(nchemMAX), save :: &
                        dti             ! variable timestep*(c+1)/(c+2)
    real, dimension(nchemMAX), save :: &
                        coeff1,coeff2,cc ! coefficients for variable timestep
    integer :: nextraiter

!======================================================

    if ( first_call ) then
       call makedt(dti,nchem,coeff1,coeff2,cc,dt_advec)
       first_call = .false.
    endif

!======================================================


    !**  toiter gives the number of iterations used in TWOSTEP. 
    !**  Use more iterations near ground:

    toiter(KCHEMTOP:5)        = 1    ! Upper levels - slow chemistry
    toiter(6:KEMISTOP-1)      = 2    ! Medium and cloud levels 
    toiter(KEMISTOP:KMAX_MID) = 3    ! Near-ground, emis levels



    !** Comments: Only NO2+O3->H+ +NO3- at night time 
    !   and in the8 lowest layers and if rh>0.5

    if (ONLY_NIGHT) call set_night_rct(rct,rh,i,j)  ! Only for ACID version


    !** Establishment of initial conditions:
    !   Previous concentrations are estimated by the current
    !   minus Dchem because the current may be changed by
    !   processes outside the chemistry:

    do k = 2, KMAX_MID

       xnew(:) = xn_2d(:,k)

       x(:)    = xn_2d(:,k) - Dchem(:,k,i,j)*dti(1)*1.5
       x(:)    = max (x(:), 0.0)

 
       !*************************************
       !     Start of integration loop      *
       !*************************************


       do ichem = 1, nchem

          do n=1,NSPEC_TOT

             xextrapol = xnew(n) + (xnew(n)-x(n)) *cc(ichem)
             xold(n) = coeff1(ichem)*xnew(n) - coeff2(ichem)*x(n)
             x(n) = xnew(n)
             xnew(n) = xextrapol

          enddo

          dt2  =  dti(ichem) !*(1.0+cc(ichem))/(1.0+2.0*cc(ichem))

          where ( xnew(:) < CPINIT  ) 
             xnew(:) = CPINIT
          end where

          !== Here comes all chemical reactions

            include 'My_Reactions.inc' 

       end do 
 
       !*************************************
       !     End of integration loop        *
       !*************************************


       !**  Saves tendencies Dchem and returns the new concentrations:

            Dchem(:,k,i,j) = (xnew(:) - xn_2d(:,k))*dt_advec_inv
            xn_2d(:,k) = xnew(:)

    enddo ! End of vertical k-loop

  end subroutine chemistry

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine  makedt(dti,nchem,coeff1,coeff2,cc,dt_tot)

!=====================================================================
! Makes coefficients for two-step (written by Peter Wind, Febr. 2003)
! The formulas for coeff1, coeff2 and dti can be found in:
! J.G. Verwer and D. Simpson, "Explicit methods for stiff ODEs from
!  atmospheric chemistry", Aplied Numerical Mathematics 18 (1995) 413
!
! Note: It is better to take first some small steps, and then
!       larger steps, than to increase the timestep gradually.
!=====================================================================

 implicit none

 real, intent(in) :: dt_tot
 real, dimension(nchemMAX),intent(out) :: dti,coeff1,coeff2,cc
 integer,                  intent(out) :: nchem

 real    :: ttot,dt_first,dt_max,dtleft,tleft,step,dt(nchemMAX)
 integer :: i,j
!_________________________

  nchem=12 !number of chemical timesteps inside dt_advec

!/ Used only for 50km resolution and dt_advec=1200 seconds:
!.. timesteps from 6 to 12
  dt=(dt_advec-100.0)/(nchem-5)
!.. first five timesteps 
  dt(1)=20.0
  dt(2)=20.0
  dt(3)=20.0
  dt(4)=20.0
  dt(5)=20.0

!/ ** For smaller scales, but not tested
   if(dt_advec<520.0)then
      nchem=5+int((dt_advec-100.0)/60.0)
      dt=(dt_advec-100.0)/(nchem-5)
      dt(1)=20.0
      dt(2)=20.0
      dt(3)=20.0
      dt(4)=20.0
      dt(5)=20.0
   endif
   if(dt_advec<=100.)then
      nchem=int(dt_advec/20.0)+1
      dt=(dt_advec)/(nchem)
   endif
!/ **

   call CheckStop(dt_advec<20.0,"Error in Solver/makedt: dt_advec too small!")

   call CheckStop(nchem>nchemMAX,"Error in Solver/makedt: nchemMAX too small!")

   nchem=min(nchemMAX,nchem)

    if(me == 0) then

      write(*,*)'Number of timesteps in Solver: ',nchem
      27 format('timestep ',I,F13.6,' total: ',F13.6)

      ttot=0.0
      do i=1,nchem
         ttot=ttot+dt(i)
         write(*,27)i,dt(i),ttot
      enddo

      !check that we are using consistent timesteps
      call CheckStop(abs(ttot-dt_advec)>1.E-5, &
              "Error in Solver/makedt: dt_advec and dt not compatible")

    endif

!.. Help variables from Verwer & Simpson 
       cc(1)=1.0
       coeff2(1)=1.0/(cc(1)**2+2*cc(1))
       coeff1(1)=(cc(1)+1)**2*coeff2(1)
       dti(1)=((cc(1)+1)/(cc(1)+2))*dt(1)

       do i=2,nchem
         cc(i)=dt(i-1)/dt(i)
         coeff2(i)=1.0/(cc(i)**2+2.0*cc(i))
         coeff1(i)=((cc(i)+1.0)**2)*coeff2(i)
         dti(i)=((cc(i)+1.0)/(cc(i)+2.0))*dt(i)
         cc(i)=1.0/cc(i)
       enddo

end subroutine makedt
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module Chemsolver_ml
