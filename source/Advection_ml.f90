! <Advection_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
                    Module Advection_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DESCRIPTION
!
! This module contains the routines for advection and diffusion.
! The sequence of advection in x y or z direction and vertical diffusion,
! is controlled in the advecdiff routine.
!
! The horizontal advection is performed in the advx and advy routines using
! "Bott's fourth order scheme". The routine preadvx and preadvy take care of
! the transfer of information between processors before the advection step.
!
! The advvk routine performs the vertical advection. Bott's second order 
! scheme with variable grid distance is used. 
! The calculation of the coefficients used for this scheme is done in the 
! routine vgrid.
!
! Notes from Peter; 7/11/01
! About the division by the surface pressure (p*):
! In my opinion the best way is to divide by the advected p*, (corresponding
! to the option where ADVEC_TYPE==1).  This should ensure
! that, in the case of a uniform mixing ratio, we end up with a uniform mixing
! ratio, whatever the meteo.  The problem is that if the meteo is not
! consistent (air is "created", or surface pressure does not corespond to the
! quantity of air) the total weight of species may vary, creating problems for
! the mass budget. I see however no simple solution for this problem.
!
!
! Peter, January 2002: The advecdiff routine has been completely reorganised,
! in order to allow for flexible timesteps. The timestep can now be large.
! The advecdiff routine will divide dt_advec in several advection steps if the
! CFL condition is not met. For small dt_advec (600s in a 50x50 km2 grid))
! these changes should usually have no effect on the result.
!
! Peter, January 2003: The vertical diffusion and the division by p* have been
! extracted out of the advvdifvk routine. Now they are done separately.
! The number of diffusion iterations can be chosen (ndiff).
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use Chemfields_ml,     only : xn_adv
  use ChemSpecs,         only : NSPEC_ADV
!CRM  use ChemSpecs_adv_ml , only : NSPEC_ADV
  use CheckStop_ml,      only : CheckStop
  use Convection_ml,     only : convection_pstar,convection_Eta
  use GridValues_ml,     only : GRIDWIDTH_M,xm2,xmd,xm2ji,xmdji, &
                                carea,xm_i, Pole_Singular,dA,dB
  use Io_ml,             only : datewrite
  use ModelConstants_ml, only : KMAX_BND,KMAX_MID,NMET, nstep, nmax, &
                  dt_advec, dt_advec_inv,  PT,Pref, KCHEMTOP, NPROCX,NPROCY,NPROC, &
                  FORECAST,& 
                  USE_CONVECTION,DEBUG_ADV
  use MetFields_ml,      only : ps,sdot,Etadot,SigmaKz,EtaKz,u_xmj,v_xmi,cnvuf,cnvdf&
                                ,uw,ue,vs,vn
  use MassBudget_ml,     only : fluxin,fluxout
  use My_Timing_ml,      only : Code_timer, Add_2timing, tim_before,tim_after
  use Par_ml,            only : MAXLIMAX,MAXLJMAX,GJMAX,GIMAX,me,mex,mey,&
            li0,li1,lj0,lj1 ,limax,ljmax, gi0, IRUNBEG,gj0, JRUNBEG &
           ,neighbor,WEST,EAST,SOUTH,NORTH,NOPROC            &
           ,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2
  use PhysicalConstants_ml, only : GRAV ! gravity

  implicit none
  private

  INCLUDE 'mpif.h'
  INTEGER MPISTATUS(MPI_STATUS_SIZE)
  integer, private, parameter :: NADVS      =  3

  real, private, save, allocatable,dimension(:)  ::  dhs1, dhs1i, dhs2i

!  for vertical advection (nonequidistant spacing)
  real, private, save, allocatable, dimension(:,:,:)  ::  alfnew
  real, private, save, dimension(3)  ::  alfbegnew,alfendnew

!  real, private,save,allocatable, dimension(:,:,:) :: uw,ue
!  real, private,save,allocatable, dimension(:,:,:) :: vs,vn

  integer, public, parameter :: ADVEC_TYPE = 1 ! Divides by advected p*
! integer, public, parameter :: ADVEC_TYPE = 2 ! Divides by "meteorologically" 
                                               ! advected p* 

  public :: assign_dtadvec
  public :: assign_nmax
  public :: alloc_adv_arrays
  public :: vgrid
  public :: vgrid_Eta
  public :: advecdiff
  public :: advecdiff_poles
  public :: advecdiff_Eta
!  public :: adv_var
!  public :: adv_int

  private :: advvk
  private :: advx
  private :: advy
  private :: preadvx
  private :: preadvy

   ! Checks & warnings
   ! introduced after getting Nan when using "poor" meteo can give this too.
   !  ps3d can get zero values when winds are extremely divergent (empty a 
   !  gridcell for air). This seems to happen only very occasionally (one 
   !  gridcell, once every week for instance); & does not harm results 
   !  significantly, at least much less than the poor metdata does anyway.
   !  Still, we need to know about it.

    integer, private, save :: nWarnings = 0
    integer, private, parameter :: MAX_WARNINGS = 100

  contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine assign_dtadvec(GRIDWIDTH_M)
!
! dt_advec is set according to the grid resolution
! The choosed timestep should lead to a Courant number <1 for
! "normal" wind speeds, but this is not a strict limitation.
!
! The values of dt_advec must be an integer fraction of 3600
!
! The values put here are only suggestions
!

    implicit none
    real, intent(in) ::GRIDWIDTH_M

    if(GRIDWIDTH_M>76000.0) dt_advec=1800.0
    if(GRIDWIDTH_M<61000.0) dt_advec=1200.0
    if(GRIDWIDTH_M<21000.0) dt_advec= 900.0
    if(GRIDWIDTH_M<11000.0) dt_advec= 600.0
    if(GRIDWIDTH_M< 6000.0) dt_advec= 300.0

! GEMS025 domain 0.25 deg resol --> GRIDWIDTH_M~=27.8 km --> dt_advec=1200.0
! MACC02  domain 0.20 deg resol --> GRIDWIDTH_M~=22.2 km --> dt_advec=1200.0

!check that it is allowed:
    call CheckStop(mod(3600,nint(dt_advec)).ne.0, "3600/dt_advec must be an integer")

    dt_advec_inv=1.0/dt_advec

   if(me==0)write(*,fmt="(a,F8.1,a)")' advection time step (dt_advec) set to: ',dt_advec,' seconds'

   call alloc_adv_arrays!should be moved elsewhere

  end subroutine assign_dtadvec

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine assign_nmax(metstep)

    implicit none
    integer, intent(in) :: metstep

!     Assigne number of time-steps for the inner time-loop (over 3 hours)
!     from dt_advec

    call CheckStop(mod(3600*metstep,nint(dt_advec)).ne.0, "3600*metstep/dt_advec must be an integer")

   ! Use nint for safety anyway:

    nmax = nint(  (3600*metstep)/dt_advec )

    if (me .eq. 0) then
!      write(6,*)
!      write(6,*)'**********************************************'
      write(6,fmt="(I3,a,I2,a)")nmax,' advection steps within each metstep (',metstep,' hours)'
!      write(6,*)'**********************************************'
!      write(6,*)
    endif

  end subroutine assign_nmax
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine advecdiff
!___________________________________________________________________________________
!Flexible timestep. Peter Wind january-2002
!
! dt_advec : time interval between two advections calls
! (controls time splitting between advection and chemistry )
!
! dt_xys : time intervall between vertical and horizontal advection steps
!(controls time splitting between vertical and horizontal advection)
! There is one sequence (z),(x,y,y,x),(z) during each dt_xys
!
! dt_xy : time intervall for horizontal advection iterations
!(controls time splitting between x and y advection)
! There is one sequence x,y,y,x during each dt_xy
!
! dt_s : time intervall for vertical advection iterations
!
! dt_advec >= dt_xys >= max(dt_xy, dt_s)
!

    implicit none

!    local

    integer i,j,k,info
    real dth
    real xntop(NSPEC_ADV,MAXLIMAX,MAXLJMAX)
    real xnw(3*NSPEC_ADV,MAXLJMAX),xne(3*NSPEC_ADV,MAXLJMAX)
    real xnn(3*NSPEC_ADV,MAXLIMAX),xns(3*NSPEC_ADV,MAXLIMAX)
    real ps3d(MAXLIMAX,MAXLJMAX,KMAX_MID),psi
    real psw(3,MAXLJMAX),pse(3,MAXLJMAX)
    real psn(3,MAXLIMAX),pss(3,MAXLIMAX)
    real ds3(2:KMAX_MID),ds4(2:KMAX_MID)

    real xcmax(KMAX_MID),ycmax(KMAX_MID),scmax,sdcmax,c_max
    real dt_xysmax,dt_xymax(KMAX_MID),dt_smax
    real dt_xys,dt_xy(KMAX_MID),dt_s,div,minps3d
    integer niterxys,niterxy(KMAX_MID),niters,nxy,ndiff
    integer iterxys,iterxy,iters


    call Code_timer(tim_before)

    if(KCHEMTOP==2)then
      xntop(:,:,:)=xn_adv(:,:,:,1)
    endif

!   convert from mixing ratio to concentration before advection

    do k = 1,KMAX_MID
      do j = 1,ljmax
        do i = 1,limax
          xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*(ps(i,j,1)-PT)
          ps3d(i,j,k) = ps(i,j,1) - PT
        end do
      end do
    end do

    call Add_2timing(25,tim_after,tim_before,"advecdiff:ps")

!     time-splitting is used for the physical and chemical operators.
!     second-order accuracy in time is obtained by alternating the order
!     of the advx and advy operators from one time-step to another.

!
! Determine timestep for horizontal advection.
!
! Courant criterion, which takes into account the mapping factor xm2:
! left face:    xm2(i)    u_xmj(i)  dt/dx < 1  when u_xmj(i) > 0
! right face:   xm2(i) |u_xmj(i-1)| dt/dx < 1  when u_xmj(i-1) < 0
!
! In the case where the flux is streaming out of the cell i from both faces,
! then the total should be < 1:
! xm2(i) |u_xmj(i-1)| dt/dx + xm2(i) u_xmj(i) dt/dx < 1    for u_xmj(i-1)<0 and u_xmj(i)>0
!
! The three conditions can be written as:
!
! max(xm2(i)*u_xmj(i)*dt/dx , 0.0) - min(xm2(i)*u_xmj(i-1)*dt/dx , 0.0) < 1
!
! or equivalently:
! dt < dx / ( max(xm2(i)*u_xmj(i) , 0.0) - min(xm2(i)*u_xmj(i-1) , 0.0) )
!
! In the case of variable cell size, dx is defined as dx(i) in these formula.
!
! The value 1.e-30 is to ensure that we don't divide by 0 when all the velocities are 0.

    dth = dt_advec/GRIDWIDTH_M

    do k=1,KMAX_MID
      xcmax(k) = maxval(                                                     &
                 max(u_xmj(1:limax  ,1:ljmax,k,1)*xm2(1:limax,1:ljmax),1.e-30)   &
                -min(u_xmj(0:limax-1,1:ljmax,k,1)*xm2(1:limax,1:ljmax),0.0   ))
      ycmax(k) = maxval(                                                     &
                 max(v_xmi(1:limax,1:ljmax  ,k,1)*xm2(1:limax,1:ljmax),1.e-30)   &
                -min(v_xmi(1:limax,0:ljmax-1,k,1)*xm2(1:limax,1:ljmax),0.0   ))
      xcmax(k) = max(xcmax(k),ycmax(k))
    enddo

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,xcmax,KMAX_MID,MPI_DOUBLE_PRECISION, &
                       MPI_MAX,MPI_COMM_WORLD,INFO)

    do k=1,KMAX_MID
      dt_xymax(k)=GRIDWIDTH_M/xcmax(k)
    enddo

!Courant number in vertical sigma coordinates:  sigmadot*dt/deltasigma
!
!Note that dhs1(k+1) denotes thickness of layer k
!     and sdot(k+1) denotes sdot at the boundary between layer k and k+1
!
!flux through wall k+1:  sdot(k+1) *dt/dhs1(k+1)<1   for sdot(k+1)>0
!                       |sdot(k+1)|*dt/dhs1(k+2)<1   for sdot(k+1)<0
!
!layer k: sdot(k+1)*dt/dhs1(k+1) + |sdot(k)|*dt/dhs1(k+1) <1 for sdot(k+1)>0 and sdot(k)<0
!
!total out of layer k: max(sdot(1:limax,1:ljmax,k+1,1),0.0)-min(sdot(1:limax,1:ljmax,k,1),0.0)
!
    scmax = 1.e-30
    do k = 1,KMAX_MID
      sdcmax = maxval(max(sdot(1:limax,1:ljmax,k+1,1),0.0)   &
                     -min(sdot(1:limax,1:ljmax,k  ,1),0.0))
      scmax  = max(sdcmax/dhs1(k+1),scmax)
    enddo

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,scmax,1,MPI_DOUBLE_PRECISION, &
                       MPI_MAX,MPI_COMM_WORLD,INFO)
    dt_smax = 1./scmax

    dt_xysmax = max(dt_smax,maxval(dt_xymax(1:KMAX_MID)))
    niterxys = int(dt_advec/dt_xysmax)+1

    dt_xys = dt_advec/real(niterxys)
    niters = int(dt_xys/dt_smax)+1

    dt_s = dt_xys/real(niters)
!   if(me.eq.0)then
!     write(*,45)dt_xysmax,dt_xys,niterxys
!     write(*,45)dt_smax,dt_s,niters
!   endif
    nxy=0
    do k=1,KMAX_MID
      niterxy(k) = int(dt_xys/dt_xymax(k))+1
      dt_xy(k) = dt_xys/real(niterxy(k))
      nxy=nxy+niterxy(k)
!     if(me.eq.0)then
!       write(*,46)k,dt_xymax(k),dt_xy(k),niterxy(k)
!     endif
      enddo
      if(me.eq.0)then
        c_max=maxval(xcmax(1:KMAX_MID)*dt_xy(1:KMAX_MID)/GRIDWIDTH_M)
        if (DEBUG_ADV) write(*,47)niterxys-1,nxy-KMAX_MID,niters-1,       &
                  c_max,xcmax(KMAX_MID)*dt_xy(KMAX_MID)/GRIDWIDTH_M
      endif
!     if(me.eq.0)then
!       write(*,47)nxy,nxy/20.
!     endif
!45    format(2F12.4,I6)
!46    format('k = ',I6,2F12.4,I6,2F12.4,I6)
47    format('extra iterations (xyz,xy,z), C_max, C_max_surface:',3I3,2F7.3)

      call Add_2timing(20,tim_after,tim_before,"advecdiff:synchronization")

! Start xys advection loop:
      iterxys = 0
      do while (iterxys < niterxys)
        if(mod(nstep,2) /= 0 .or. iterxys /= 0)then !start a xys sequence

          iterxys = iterxys + 1
          do k = 1,KMAX_MID
            dth = dt_xy(k)/GRIDWIDTH_M
            do iterxy=1,niterxy(k)

              ! send/receive in x-direction
              call preadvx(1100+k                           &
                  ,xn_adv(1,1,1,k),ps3d(1,1,k),u_xmj(0,1,k,1)   &
                  ,xnw,xne                                  &
                  ,psw,pse)

              ! x-direction
              do j = lj0,lj1
                call advx(                                   &
                       u_xmj(0,j,k,1),uw(j,k,1),ue(j,k,1)        &
                      ,xn_adv(1,1,j,k),xnw(1,j),xne(1,j)     &
                      ,ps3d(1,j,k),psw(1,j),pse(1,j)         &
                      ,xm2(0,j),xmd(0,j)                     &
                      ,dth,carea(k))
              enddo

              call Add_2timing(21,tim_after,tim_before,"advecdiff:preadvx,advx")

              ! send/receive in y-direction
              call preadvy(1300+k                            &
                  ,xn_adv(1,1,1,k),ps3d(1,1,k),v_xmi(1,0,k,1)    &
                  ,xns, xnn                                  &
                  ,pss, psn)

              ! y-direction
              do i = li0,li1
                call advy(                                   &
                       v_xmi(i,0,k,1),vs(i,k,1),vn(i,k,1)        &
                      ,xn_adv(1,i,1,k),xns(1,i),xnn(1,i)     &
                      ,ps3d(i,1,k),pss(1,i),psn(1,i)         &
                      ,xm2ji(0,i),xmdji(0,i)                 &
                      ,dth,carea(k))
              enddo

              call Add_2timing(23,tim_after,tim_before,"advecdiff:preadvy,advy")  
            enddo !iterxy horizontal (xy) advection
          enddo !k horizontal (xy) advection

          do iters=1,niters

            ! perform only vertical advection
            do j = lj0,lj1
              do i = li0,li1
                call advvk(xn_adv(1,i,j,1),ps3d(i,j,1),sdot(i,j,1,1),dt_s)
              enddo
            enddo

          enddo ! vertical (s) advection

          call Add_2timing(24,tim_after,tim_before,"advecdiff:advvk")

        else  !start a yxs sequence

          iterxys = iterxys + 1

          do k = 1,KMAX_MID
            dth = dt_xy(k)/GRIDWIDTH_M
            do iterxy=1,niterxy(k)

              ! send/receive in y-direction
              call preadvy(1300+k                            &
                  ,xn_adv(1,1,1,k),ps3d(1,1,k),v_xmi(1,0,k,1)    &
                  ,xns, xnn                                  &
                  ,pss, psn)

              ! y-direction
              do i = li0,li1
                call advy(                                   &
                       v_xmi(i,0,k,1),vs(i,k,1),vn(i,k,1)        &
                      ,xn_adv(1,i,1,k),xns(1,i),xnn(1,i)     &
                      ,ps3d(i,1,k),pss(1,i),psn(1,i)         &
                      ,xm2ji(0,i),xmdji(0,i)                 &
                      ,dth,carea(k))
              enddo

              call Add_2timing(23,tim_after,tim_before,"advecdiff:preadvy,advy")

              ! send/receive in x-direction
              call preadvx(1100+k                            &
                  ,xn_adv(1,1,1,k),ps3d(1,1,k),u_xmj(0,1,k,1)    &
                  ,xnw,xne                                   &
                  ,psw,pse)

              ! x-direction
              do j = lj0,lj1
                call advx(                                   &
                       u_xmj(0,j,k,1),uw(j,k,1),ue(j,k,1)        &
                      ,xn_adv(1,1,j,k),xnw(1,j),xne(1,j)     &
                      ,ps3d(1,j,k),psw(1,j),pse(1,j)         &
                      ,xm2(0,j),xmd(0,j)                     &
                      ,dth,carea(k))
              enddo

              call Add_2timing(21,tim_after,tim_before,"advecdiff:preadvx,advx")

            enddo !iterxy horizontal (xy) advection
          enddo !k horizontal (xy) advection

          do iters=1,niters

            ! perform only vertical advection
            do j = lj0,lj1
              do i = li0,li1
                call advvk(xn_adv(1,i,j,1),ps3d(i,j,1),sdot(i,j,1,1),dt_s)
              enddo
            enddo

          enddo ! vertical (s) advection

          call Add_2timing(24,tim_after,tim_before,"advecdiff:advvk")

        endif ! yxs sequence
      enddo

      if(USE_CONVECTION)then
         do k=1,KMAX_MID
            do j = lj0,lj1
               do i = li0,li1
                  psi = (ps(i,j,1) - PT)/max(ps3d(i,j,k),1.0)
                  xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                  ps3d(i,j,k) = ps(i,j,1) - PT
               enddo
            enddo
            minps3d = minval( ps3d(li0:li1, lj0:lj1, k) )
            if ( nWarnings < MAX_WARNINGS  .and.minps3d<1.0) then
              call datewrite("WARNING ps3d < 1",k,  (/ minps3d /) )
              nWarnings = nWarnings + 1
            end if
         enddo

         call convection_pstar(ps3d,dt_advec)

      endif


! division by p*, to transform back into mixing ratios units
      if(ADVEC_TYPE==2)then !this option is "mass conservative" (uses metadvected p*)
        div = 1./real(nmax-(nstep-1))
        do j = lj0,lj1
          do i = li0,li1
            psi = 1./(ps(i,j,1)+(ps(i,j,2) - ps(i,j,1))*div-PT)
            do k=1,KMAX_MID
              xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
            enddo
          enddo
        enddo

      elseif(ADVEC_TYPE==1)then !this option is recommended (?).
                                !It is "mixing ratio  conservative" (uses advected p*)
        do k=1,KMAX_MID
          do j = lj0,lj1
            do i = li0,li1
              psi = 1./max(1.0,ps3d(i,j,k))
              xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
            enddo
          enddo
            minps3d = minval( ps3d(li0:li1, lj0:lj1, k) )
            if ( nWarnings < MAX_WARNINGS .and.minps3d<1.0 ) then
              call datewrite("WARNING:B ps3d < 1",k,  (/ minps3d /) )
              nWarnings = nWarnings + 1
            end if
        enddo

      else
        do k=1,KMAX_MID
          do j = lj0,lj1
            do i = li0,li1
              xn_adv(:,i,j,k) = xn_adv(:,i,j,k)/(ps(i,j,1)-PT)
            enddo
          enddo
        enddo

      endif

      call Add_2timing(25,tim_after,tim_before,"advecdiff:ps")

!vertical diffusion
      ndiff = 1 !number of vertical diffusion iterations (the larger the better)
      do k = 2,KMAX_MID
        ds3(k) = dt_advec*dhs1i(k)*dhs2i(k)
        ds4(k) = dt_advec*dhs1i(k+1)*dhs2i(k)
      enddo

! sum is conserved under vertical diffusion
!     sum = 0.
!     do k=1,KMAX_MID
!       sum = sum + xn_adv(1,4,4,k)/dhs1i(k+1)
!     enddo
!     write(*,*)'sum before diffusion ',me,sum
      do j = lj0,lj1
        do i = li0,li1
          call vertdiffn(xn_adv(1,i,j,1),SigmaKz(i,j,1,1),ds3,ds4,ndiff)
        enddo
      enddo
!     sum = 0.
!     do k=1,KMAX_MID
!        sum = sum + xn_adv(1,4,4,k)/dhs1i(k+1)
!     enddo
!     write(*,*)'sum after diffusion ',me,sum
      call Add_2timing(22,tim_after,tim_before,"advecdiff:diffusion")


    if(lj0.ne.1)then
      do k=KCHEMTOP,KMAX_MID
        do i = 1,limax
          xn_adv(:,i,1,k) = xn_adv(:,i,1,k)/(ps(i,1,1)-PT)
        enddo
      enddo
    endif
    if(li0.ne.1)then
      do k=KCHEMTOP,KMAX_MID
        do j=lj0,lj1
          xn_adv(:,1,j,k) = xn_adv(:,1,j,k)/(ps(1,j,1)-PT)
        enddo
      enddo
    endif
    if(li1.ne.limax)then
      do k=KCHEMTOP,KMAX_MID
        do j=lj0,lj1
          xn_adv(:,limax,j,k) = xn_adv(:,limax,j,k)/(ps(limax,j,1)-PT)
        enddo
      enddo
    endif
    if(lj1.ne.ljmax)then
      do k=KCHEMTOP,KMAX_MID
        do i = 1,limax
          xn_adv(:,i,ljmax,k) = xn_adv(:,i,ljmax,k)/(ps(i,ljmax,1)-PT)
        enddo
      enddo
    endif

    if(KCHEMTOP==2)then

!since the xn_adv are changed it corresponds to a flux in or out of the system:

      do i  = li0,li1
        do j = lj0,lj1
          where(xn_adv(:,i,j,1) .gt. xntop(:,i,j))
            fluxout(:) = fluxout(:) + &
              (xn_adv(:,i,j,1)-xntop(:,i,j))*(ps(i,j,1)-PT)*carea(1)*xmd(i,j)
          elsewhere
            fluxin(:) = fluxin(:) + &
              (xntop(:,i,j)-xn_adv(:,i,j,1))*(ps(i,j,1)-PT)*carea(1)*xmd(i,j)
          end where
        enddo
      enddo
      xn_adv(:,:,:,1) = xntop(:,:,:)

    endif

    return

  end subroutine advecdiff

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine advecdiff_poles
    !___________________________________________________________________________________
    !Uses more robust options:
    !1)Advect i,j directions independently and 1D with own timestep
    !2)Do not advect but only "mix" the concentrations near poles ("near"
    !  poles is determined by NITERXMAX.
    !
    !1/10/2012: divide by ps3d (p*) after each partial advection (x,y or z 
    !direction)
    !
    !Flexible timestep. Peter Wind january-2002
    !
    ! dt_advec : time interval between two advections calls
    ! (controls time splitting between advection and chemistry )
    !
    ! dt_xys : time intervall between vertical and horizontal advection steps
    !(controls time splitting between vertical and horizontal advection)
    ! There is one sequence (z),(x,y,y,x),(z) during each dt_xys
    !
    ! dt_xy : time intervall for horizontal advection iterations
    !(controls time splitting between x and y advection)
    ! There is one sequence x,y,y,x during each dt_xy
    !
    ! dt_s : time intervall for vertical advection iterations
    !
    ! dt_advec >= dt_xys >= max(dt_xy, dt_s)
    !

    implicit none

    !    local

    integer i,j,k,n,info
    real dth
    real xntop(NSPEC_ADV,MAXLIMAX,MAXLJMAX)
    real xnw(3*NSPEC_ADV),xne(3*NSPEC_ADV)
    real xnn(3*NSPEC_ADV),xns(3*NSPEC_ADV)
    real ps3d(MAXLIMAX,MAXLJMAX,KMAX_MID),psi
    real psw(3),pse(3)
    real psn(3),pss(3)
    real ds3(2:KMAX_MID),ds4(2:KMAX_MID)

    real xcmax(KMAX_MID,GJMAX),ycmax(KMAX_MID,GIMAX),scmax,sdcmax
    real dt_smax,dt_s,div
    real dt_x(MAXLJMAX,KMAX_MID),dt_y(MAXLIMAX,KMAX_MID)
    real dt_xmax(MAXLJMAX,KMAX_MID),dt_ymax(MAXLIMAX,KMAX_MID)
    integer niterx(MAXLJMAX,KMAX_MID),nitery(MAXLIMAX,KMAX_MID)

    integer niterxys,niters,nxy,ndiff
    integer iterxys,iters,iterx,itery,nxx,nxxmin,nyy
    logical,save :: firstcall = .true.

    integer ::isum,isumtot,iproc
    real :: xn_advjktot(NSPEC_ADV),xn_advjk(NSPEC_ADV),rfac,minps3d


    !NITERXMAX=max value of iterations accepted for fourth order Bott scheme.
    !If the calculated number of iterations (determined from Courant number)
    !exceeds NITERXMAX, the advection is not done, but instead all the mixing
    !ratio along that line are averaged (1D).
    !This case can arises where there is a singularity close to the
    !poles in long-lat coordinates.
    integer,parameter :: NITERXMAX=10


    call Code_timer(tim_before)

    if(firstcall)then
       if(NPROCY>2.and.me==0.and.Pole_Singular>1)then
          write(*,*)&
               'COMMENT: Advection routine will work faster if NDY = 2 (or 1)'
       elseif(NPROCY>1.and.me==0.and.Pole_Singular==1)then
          write(*,*)&
               'COMMENT: Advection routine will work faster if NDY = 1'
       endif
    endif

    if(KCHEMTOP==2)then
       xntop(:,:,:)=xn_adv(:,:,:,1)
    endif

    ! convert from mixing ratio to concentration before advection
    do k = 1,KMAX_MID
       do j = 1,ljmax
          do i = 1,limax
             xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*(ps(i,j,1)-PT)
             ps3d(i,j,k) = ps(i,j,1) - PT
          end do
       end do
    end do

    call Add_2timing(25,tim_after,tim_before,"advecdiff:ps")

    !     time-splitting is used for the physical and chemical operators.
    !     second-order accuracy in time is obtained by alternating the order
    !     of the advx and advy operators from one time-step to another.

    !
    ! Determine timestep for horizontal advection.
    !
    ! Courant criterion, which takes into account the mapping factor xm2:
    ! left face:    xm2(i)    u_xmj(i)  dt/dx < 1  when u_xmj(i) > 0
    ! right face:   xm2(i) |u_xmj(i-1)| dt/dx < 1  when u_xmj(i-1) < 0
    !
    ! In the case where the flux is streaming out of the cell i from both faces,
    ! then the total should be < 1:
    ! xm2(i) |u_xmj(i-1)| dt/dx + xm2(i) u_xmj(i) dt/dx < 1    for u_xmj(i-1)<0 and u_xmj(i)>0
    !
    ! The three conditions can be written as:
    !
    ! max(xm2(i)*u_xmj(i)*dt/dx , 0.0) - min(xm2(i)*u_xmj(i-1)*dt/dx , 0.0) < 1
    !
    ! or equivalently:
    ! dt < dx / ( max(xm2(i)*u_xmj(i) , 0.0) - min(xm2(i)*u_xmj(i-1) , 0.0) )
    !
    ! In the case of variable cell size, dx is defined as dx(i) in these formula.
    !
    ! The value 1.e-30 is to ensure that we don't divide by 0 when all the velocities are 0.

    dth = dt_advec/GRIDWIDTH_M
    xcmax=0.0
    ycmax=0.0
    do k=1,KMAX_MID
       do j=1,ljmax
          xcmax(k,j+gj0-1) = maxval(                                         &
               max(u_xmj(1:limax  ,j,k,1)*xm2(1:limax,j),1.e-30)   &
               -min(u_xmj(0:limax-1,j,k,1)*xm2(1:limax,j),0.0   ))
       enddo
       do i=1,limax
          ycmax(k,i+gi0-1) = maxval(                                         &
               max(v_xmi(i,1:ljmax  ,k,1)*xm2(i,1:ljmax),1.e-30)   &
               -min(v_xmi(i,0:ljmax-1,k,1)*xm2(i,1:ljmax),0.0   ))
       enddo
    enddo

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,xcmax,KMAX_MID*gjmax,MPI_DOUBLE_PRECISION, &
         MPI_MAX,MPI_COMM_WORLD,INFO)

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,ycmax,KMAX_MID*gimax,MPI_DOUBLE_PRECISION, &
         MPI_MAX,MPI_COMM_WORLD, INFO)

    do i=1,limax
       do k=1,KMAX_MID
          dt_ymax(i,k)=GRIDWIDTH_M/ycmax(k,i+gi0-1)
       enddo
    enddo
    do j=1,ljmax
       do k=1,KMAX_MID
          dt_xmax(j,k)=GRIDWIDTH_M/xcmax(k,j+gj0-1)
       enddo
    enddo

    niterx=1
    do k=1,KMAX_MID
       do j=1,ljmax
          niterx(j,k) = int(dt_advec/dt_xmax(j,k))+1
          dt_x(j,k) = dt_advec/real(niterx(j,k))
          !if(me==0)write(*,*)'x',me,j,k,niterx(j,k),xcmax(k,j+gj0-1)
       enddo
    enddo

    do k=1,KMAX_MID
       do i=1,limax
          nitery(i,k) = int(dt_advec/dt_ymax(i,k))+1
          dt_y(i,k) = dt_advec/real(nitery(i,k))
       enddo
    enddo

    !Courant number in vertical sigma coordinates:  sigmadot*dt/deltasigma
    !
    !Note that dhs1(k+1) denotes thickness of layer k
    !     and sdot(k+1) denotes sdot at the boundary between layer k and k+1
    !
    !flux through wall k+1:  sdot(k+1) *dt/dhs1(k+1)<1   for sdot(k+1)>0
    !                       |sdot(k+1)|*dt/dhs1(k+2)<1   for sdot(k+1)<0
    !
    !layer k: sdot(k+1)*dt/dhs1(k+1) + |sdot(k)|*dt/dhs1(k+1) <1 for sdot(k+1)>0 and sdot(k)<0
    !
    !total out of layer k: max(sdot(1:limax,1:ljmax,k+1,1),0.0)-min(sdot(1:limax,1:ljmax,k,1),0.0)
    !
    scmax = 1.e-30
    do k = 1,KMAX_MID
       sdcmax = maxval(max(sdot(1:limax,1:ljmax,k+1,1),0.0)   &
            -min(sdot(1:limax,1:ljmax,k  ,1),0.0))
       scmax  = max(sdcmax/dhs1(k+1),scmax)
    enddo

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,scmax,1,MPI_DOUBLE_PRECISION, &
         MPI_MAX,MPI_COMM_WORLD,INFO)
    dt_smax = 1./scmax
42  FORMAT(A,F10.2)
    if(me==0.and. firstcall.and.DEBUG_ADV)write(*,42)'dt_smax',dt_smax
    niters = int(dt_advec/dt_smax)+1
    dt_s = dt_advec/real(niters)

    niterxys = 1

    nxy=0
    nxx=0
    nxxmin=0
    nyy=0
    do k=1,KMAX_MID
       do j=1,ljmax
          nxy=nxy+niterx(j,k)-1
          nxx=nxx+niterx(j,k)-1
          if(niterx(j,k)>NITERXMAX)then
             nxxmin=nxxmin+niterx(j,k)
          endif
       enddo
       do i=1,limax
          nxy=nxy+nitery(i,k)-1
          nyy=nyy+nitery(i,k)-1
       enddo

    enddo
    if(me.eq.0)then
       !          write(*,43)KMAX_MID*ljmax,nxx,nxxmin,KMAX_MID*limax,nyy,niters
    endif
    !43    format('total iterations x, y, k: ',I4,' +',I4,' -',I4,', ',I5,' +',I3,',',I4)

    ! stop

    call Add_2timing(20,tim_after,tim_before,"advecdiff:synchronization")

    ! Start xys advection loop:
    iterxys = 0
    do while (iterxys < niterxys)
       if(mod(nstep,2) /= 0 .or. iterxys /= 0)then !start a xys sequence

          iterxys = iterxys + 1
          do k = 1,KMAX_MID
             do j = lj0,lj1
                if(niterx(j,k)<=NITERXMAX)then
                   dth = dt_x(j,k)/GRIDWIDTH_M
                   do iterx=1,niterx(j,k)

                      ! send/receive in x-direction
                      call preadvx2(110+k+KMAX_MID*j               &
                           ,xn_adv(1,1,j,k),ps3d(1,j,k),u_xmj(0,j,k,1)&
                           ,xnw,xne                               &
                           ,psw,pse)

                      ! x-direction
                      call advx(                                   &
                           u_xmj(0,j,k,1),uw(j,k,1),ue(j,k,1)        &
                           ,xn_adv(1,1,j,k),xnw,xne               &
                           ,ps3d(1,j,k),psw,pse                   &
                           ,xm2(0,j),xmd(0,j)                     &
                           ,dth,carea(k))
                      do i = li0,li1
                         psi =(ps(i,j,1)-PT)/max(ps3d(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                         ps3d(i,j,k) = ps(i,j,1) - PT
                      enddo

                   enddo !iter

                endif
             enddo !j
             !          enddo !k horizontal (x) advection


             call Add_2timing(21,tim_after,tim_before,"advecdiff:advx")

             ! y-direction
             !          do k = 1,KMAX_MID
             do i = li0,li1
                dth = dt_y(i,k)/GRIDWIDTH_M
                do itery=1,nitery(i,k)

                   ! send/receive in y-direction
                   call preadvy2(520+k                            &
                        ,xn_adv(1,1,1,k),ps3d(1,1,k),v_xmi(1,0,k,1)    &
                        ,xns, xnn                                  &
                        ,pss, psn,i)

                   call advy(                                     &
                        v_xmi(i,0,k,1),vs(i,k,1),vn(i,k,1)            &
                        ,xn_adv(1,i,1,k),xns,xnn                   &
                        ,ps3d(i,1,k),pss,psn                       &
                        ,xm2ji(0,i),xmdji(0,i)                     &
                        ,dth,carea(k))
                   do j = lj0,lj1
                      psi =(ps(i,j,1)-PT)/max(ps3d(i,j,k),1.0)
                      xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      ps3d(i,j,k) = ps(i,j,1) - PT
                   enddo

                enddo !iter
             enddo !i
          enddo !k horizontal (y) advection


          call Add_2timing(23,tim_after,tim_before,"advecdiff:advy")



          do iters=1,niters

             ! perform only vertical advection
             do j = lj0,lj1
                do i = li0,li1
                   call advvk(xn_adv(1,i,j,1),ps3d(i,j,1),sdot(i,j,1,1),dt_s)
                enddo
             enddo
             if(iters<niters .or. iterxys < niterxys)then
                do k=1,KMAX_MID
                   do j = lj0,lj1
                      do i = li0,li1
                         psi =(ps(i,j,1)-PT)/max(ps3d(i,j,k),1.0)
                         ps3d(i,j,k) = ps(i,j,1) - PT
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      enddo
                   enddo
                enddo
             else
                !advection finished
                do k=1,KMAX_MID
                   do j = lj0,lj1
                      do i = li0,li1
                         psi =1.0/max(ps3d(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      enddo
                   enddo
                enddo
             endif
          enddo ! vertical (s) advection


          call Add_2timing(24,tim_after,tim_before,"advecdiff:advvk")

       else  !start a yxs sequence

          iterxys = iterxys + 1
          do k = 1,KMAX_MID
             do i = li0,li1
                dth = dt_y(i,k)/GRIDWIDTH_M
                do itery=1,nitery(i,k)

                   ! send/receive in y-direction
                   call preadvy2(13000+k+KMAX_MID*itery+1000*i    &
                        ,xn_adv(1,1,1,k),ps3d(1,1,k),v_xmi(1,0,k,1)    &
                        ,xns, xnn                                  &
                        ,pss, psn,i)

                   ! y-direction

                   call advy(                                    &
                        v_xmi(i,0,k,1),vs(i,k,1),vn(i,k,1)           &
                        ,xn_adv(1,i,1,k),xns,xnn                  &
                        ,ps3d(i,1,k),pss,psn                      &
                        ,xm2ji(0,i),xmdji(0,i)                    &
                        ,dth,carea(k))

                   do j = lj0,lj1
                      psi =(ps(i,j,1)-PT)/max(ps3d(i,j,k),1.0)
                      xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      ps3d(i,j,k) = ps(i,j,1) - PT
                   enddo

                enddo !iter
             enddo !i
             !         enddo !k horizontal (y) advection

 

             call Add_2timing(23,tim_after,tim_before,"advecdiff:preadvy,advy")

             !          do k = 1,KMAX_MID
             do j = lj0,lj1
                if(niterx(j,k)<=NITERXMAX)then
                   dth = dt_x(j,k)/GRIDWIDTH_M
                   do iterx=1,niterx(j,k)

                      ! send/receive in x-direction
                      call preadvx2(21000+k+KMAX_MID*iterx+1000*j  &
                           ,xn_adv(1,1,j,k),ps3d(1,j,k),u_xmj(0,j,k,1)&
                           ,xnw,xne                               &
                           ,psw,pse)

                      ! x-direction
                      call advx(                                   &
                           u_xmj(0,j,k,1),uw(j,k,1),ue(j,k,1)        &
                           ,xn_adv(1,1,j,k),xnw,xne               &
                           ,ps3d(1,j,k),psw,pse                   &
                           ,xm2(0,j),xmd(0,j)                     &
                           ,dth,carea(k))

                      do i = li0,li1
                         psi =(ps(i,j,1)-PT)/max(ps3d(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                         ps3d(i,j,k) = ps(i,j,1) - PT
                      enddo

                   enddo !iter
                endif
             enddo !j
          enddo !k horizontal (x) advection

          call Add_2timing(21,tim_after,tim_before,"advecdiff:preadvx,advx") 


          do iters=1,niters

             ! perform only vertical advection
             do j = lj0,lj1
                do i = li0,li1
                   call advvk(xn_adv(1,i,j,1),ps3d(i,j,1),sdot(i,j,1,1),dt_s)
                enddo
             enddo
             if(iters<niters .or. iterxys < niterxys)then
                do k=1,KMAX_MID
                   do j = lj0,lj1
                      do i = li0,li1
                         psi =(ps(i,j,1)-PT)/max(ps3d(i,j,k),1.0)
                         ps3d(i,j,k) = ps(i,j,1) - PT
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      enddo
                   enddo
                enddo
             else
                !advection finished
                do k=1,KMAX_MID
                   do j = lj0,lj1
                      do i = li0,li1
                         psi =1.0/max(ps3d(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      enddo
                   enddo
                enddo
             endif
          enddo ! vertical (s) advection


          call Add_2timing(24,tim_after,tim_before,"advecdiff:advvk")

       endif ! yxs sequence
    enddo



    if(USE_CONVECTION)then

       call CheckStop(ADVEC_TYPE/=1, "ADVEC_TYPE no longer supported")

       do k=1,KMAX_MID
          do j = lj0,lj1
             do i = li0,li1
                psi = (ps(i,j,1) - PT)!/max(ps3d(i,j,k),1.0)
                ps3d(i,j,k) = ps(i,j,1) - PT
                xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
             enddo
          enddo
       enddo

       call convection_pstar(ps3d,dt_advec)
       do k=1,KMAX_MID
          do j = lj0,lj1
             do i = li0,li1
                psi = 1.0/max(ps3d(i,j,k),1.0)
                xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
             enddo
          enddo
          minps3d = minval( ps3d(li0:li1, lj0:lj1, k) )
          if ( nWarnings < MAX_WARNINGS  .and.minps3d<1.0) then
             call datewrite("WARNING:C ps3d < 1",k,  (/ minps3d /) )
             nWarnings = nWarnings + 1
          end if
       enddo

    endif


    do k=1,KMAX_MID
       do j = lj0,lj1
          if(niterx(j,k)>NITERXMAX)then
             !if(mex==0)write(*,*)'Simplified advection',k,j,niterx(j,k)
             !simplified "advection": average the mixing ratios over all x
             !average 1D
             xn_advjk=0.0
             isum=0
             do i = li0,li1
                !           psi=1.0/(ps(i,j,1)-PT)
                xn_advjk(:) = xn_advjk(:)+xn_adv(:,i,j,k)!*psi
                isum=isum+1
             enddo
             !sum over all processors along i direction. mex=0 collects the sum
             !me = mex + mey*NPROCX
             if(mex>0)then

                CALL MPI_SEND(xn_advjk,8*NSPEC_ADV,MPI_BYTE, &
                     mey*NPROCX,100*mey+j+1000,MPI_COMM_WORLD,INFO)
                !receive averages from mex=0
                CALL MPI_RECV(xn_advjk,8*NSPEC_ADV,MPI_BYTE, &
                     mey*NPROCX,100*mey+j+3000,MPI_COMM_WORLD,MPISTATUS,INFO)

             else
                xn_advjktot(:) = xn_advjk(:)
                isumtot=isum
                do iproc=1,NPROCX-1
                   CALL MPI_RECV(xn_advjk,8*NSPEC_ADV,MPI_BYTE, &
                        iproc+mey*NPROCX,100*mey+j+1000,MPI_COMM_WORLD,MPISTATUS,INFO)
                   xn_advjktot(:) = xn_advjktot(:)+xn_advjk(:)
                   !             isumtot=isumtot+isum
                enddo
                rfac=1.0/GIMAX
                xn_advjk(:) = xn_advjktot(:)*rfac
                !           write(*,*)'ISUM',mey,isumtot,isum,GIMAX
                !send result to all processors in i direction
                do iproc=1,NPROCX-1
                   CALL MPI_SEND(xn_advjk,8*NSPEC_ADV,MPI_BYTE, &
                        iproc+mey*NPROCX,100*mey+j+3000,MPI_COMM_WORLD,INFO)
                enddo
             endif

             do i = li0,li1
                xn_adv(:,i,j,k)= xn_advjk(:)
             enddo

          endif
       enddo
    enddo

    call Add_2timing(25,tim_after,tim_before,"advecdiff:ps")

    ! vertical diffusion
    ndiff = 1 !number of vertical diffusion iterations (the larger the better)
    do k = 2,KMAX_MID
       ds3(k) = dt_advec*dhs1i(k)*dhs2i(k)
       ds4(k) = dt_advec*dhs1i(k+1)*dhs2i(k)
    enddo

    ! sum is conserved under vertical diffusion
    !   sum = 0.
    !   do k=1,KMAX_MID
    !      sum = sum + xn_adv(1,4,4,k)/dhs1i(k+1)
    !   enddo
    !   write(*,*)'sum before diffusion ',me,sum
    do j = lj0,lj1
       do i = li0,li1
          call vertdiffn(xn_adv(1,i,j,1),SigmaKz(i,j,1,1),ds3,ds4,ndiff)
       enddo
    enddo
    !   sum = 0.
    !   do k=1,KMAX_MID
    !      sum = sum + xn_adv(1,4,4,k)/dhs1i(k+1)
    !   enddo
    !   write(*,*)'sum after diffusion ',me,sum
    call Add_2timing(22,tim_after,tim_before,"advecdiff:diffusion")


    if(lj0.ne.1)then
       do k=KCHEMTOP,KMAX_MID
          do i = 1,limax
             xn_adv(:,i,1,k) = xn_adv(:,i,1,k)/(ps(i,1,1)-PT)
          enddo
       enddo
    endif
    if(li0.ne.1)then
       do k=KCHEMTOP,KMAX_MID
          do j=lj0,lj1
             xn_adv(:,1,j,k) = xn_adv(:,1,j,k)/(ps(1,j,1)-PT)
          enddo
       enddo
    endif
    if(li1.ne.limax)then
       do k=KCHEMTOP,KMAX_MID
          do j=lj0,lj1
             xn_adv(:,limax,j,k) = xn_adv(:,limax,j,k)/(ps(limax,j,1)-PT)
          enddo
       enddo
    endif
    if(lj1.ne.ljmax)then
       do k=KCHEMTOP,KMAX_MID
          do i = 1,limax
             xn_adv(:,i,ljmax,k) = xn_adv(:,i,ljmax,k)/(ps(i,ljmax,1)-PT)
          enddo
       enddo
    endif

    if(KCHEMTOP==2)then

       !pw since the xn_adv are changed it corresponds to a flux in or
       !   out of the system:
       do i = li0,li1
          do j = lj0,lj1
             where(xn_adv(:,i,j,1) .gt. xntop(:,i,j))
                fluxout(:) = fluxout(:) + &
                     (xn_adv(:,i,j,1)-xntop(:,i,j))*(ps(i,j,1)-PT)*carea(1)*xmd(i,j)
             elsewhere
                fluxin(:) = fluxin(:) + &
                     (xntop(:,i,j)-xn_adv(:,i,j,1))*(ps(i,j,1)-PT)*carea(1)*xmd(i,j)
             end where
          enddo
       enddo
       xn_adv(:,:,:,1) = xntop(:,:,:)

    endif

    firstcall=.false.
    return

  end subroutine advecdiff_poles

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine advecdiff_Eta
    !___________________________________________________________________________________
    !Uses more robust options:
    !1)Advect i,j directions independently and 1D with own timestep
    !2)Do not advect but only "mix" the concentrations near poles ("near"
    !  poles is determined by NITERXMAX.
    !
    !1/10/2012: divide by ps3d (p*) after each partial advection (x,y or z 
    !direction)
    !
    !Flexible timestep. Peter Wind january-2002
    !
    ! dt_advec : time interval between two advections calls
    ! (controls time splitting between advection and chemistry )
    !
    ! dt_xys : time intervall between vertical and horizontal advection steps
    !(controls time splitting between vertical and horizontal advection)
    ! There is one sequence (z),(x,y,y,x),(z) during each dt_xys
    !
    ! dt_xy : time intervall for horizontal advection iterations
    !(controls time splitting between x and y advection)
    ! There is one sequence x,y,y,x during each dt_xy
    !
    ! dt_s : time intervall for vertical advection iterations
    !
    ! dt_advec >= dt_xys >= max(dt_xy, dt_s)
    !
    ! March 2013: Eta coordinates
    ! P* = Ps-PT is replaced by (dA+dB*Ps)/(dA/Pref+dB)
    ! Both are defined by dP/dEta, but in general Eta coordinates it is not height independent 

    implicit none

    !    local

    integer i,j,k,n,info
    real dth
    real xntop(NSPEC_ADV,MAXLIMAX,MAXLJMAX)
    real xnw(3*NSPEC_ADV),xne(3*NSPEC_ADV)
    real xnn(3*NSPEC_ADV),xns(3*NSPEC_ADV)
!    real ps3d(MAXLIMAX,MAXLJMAX,KMAX_MID)
    real dpdeta(MAXLIMAX,MAXLJMAX,KMAX_MID),psi
    real psw(3),pse(3)
    real psn(3),pss(3)
    real ds3(2:KMAX_MID),ds4(2:KMAX_MID)

    real xcmax(KMAX_MID,GJMAX),ycmax(KMAX_MID,GIMAX),scmax,sdcmax
    real dt_smax,dt_s,div
    real dt_x(MAXLJMAX,KMAX_MID),dt_y(MAXLIMAX,KMAX_MID)
    real dt_xmax(MAXLJMAX,KMAX_MID),dt_ymax(MAXLIMAX,KMAX_MID)
    integer niterx(MAXLJMAX,KMAX_MID),nitery(MAXLIMAX,KMAX_MID)

    integer niterxys,niters,nxy,ndiff
    integer iterxys,iters,iterx,itery,nxx,nxxmin,nyy
    logical,save :: firstcall = .true.

    integer ::isum,isumtot,iproc
    real :: xn_advjktot(NSPEC_ADV),xn_advjk(NSPEC_ADV),rfac
    real :: dpdeta0,mindpdeta,xxdg,fac1

    !NITERXMAX=max value of iterations accepted for fourth order Bott scheme.
    !If the calculated number of iterations (determined from Courant number)
    !exceeds NITERXMAX, the advection is not done, but instead all the mixing
    !ratio along that line are averaged (1D).
    !This case can arises where there is a singularity close to the
    !poles in long-lat coordinates.
    integer,parameter :: NITERXMAX=10

    xxdg=GRIDWIDTH_M*GRIDWIDTH_M/GRAV !constant used in loops

    call Code_timer(tim_before)

    if(firstcall)then
       if(NPROCY>2.and.me==0.and.Pole_Singular>1)then
          write(*,*)&
               'COMMENT: Advection routine will work faster if NDY = 2 (or 1)'
       elseif(NPROCY>1.and.me==0.and.Pole_Singular==1)then
          write(*,*)&
               'COMMENT: Advection routine will work faster if NDY = 1'
       endif
       !Overwrite the cooefficients for vertical advection, with Eta-adpated values
       call vgrid_Eta 
    endif

    if(KCHEMTOP==2)then
       xntop(:,:,:)=xn_adv(:,:,:,1)
    endif

    ! convert from mixing ratio to concentration before advection
    do k = 1,KMAX_MID
       do j = 1,ljmax
          do i = 1,limax
             dpdeta(i,j,k) = (dA(k)+dB(k)*ps(i,j,1))/(dA(k)/Pref+dB(k))
             xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*dpdeta(i,j,k)
          end do
       end do
    end do


    call Add_2timing(25,tim_after,tim_before,"advecdiff:ps")

    !     time-splitting is used for the physical and chemical operators.
    !     second-order accuracy in time is obtained by alternating the order
    !     of the advx and advy operators from one time-step to another.

    !
    ! Determine timestep for horizontal advection.
    !
    ! Courant criterion, which takes into account the mapping factor xm2:
    ! left face:    xm2(i)    u_xmj(i)  dt/dx < 1  when u_xmj(i) > 0
    ! right face:   xm2(i) |u_xmj(i-1)| dt/dx < 1  when u_xmj(i-1) < 0
    !
    ! In the case where the flux is streaming out of the cell i from both faces,
    ! then the total should be < 1:
    ! xm2(i) |u_xmj(i-1)| dt/dx + xm2(i) u_xmj(i) dt/dx < 1    for u_xmj(i-1)<0 and u_xmj(i)>0
    !
    ! The three conditions can be written as:
    !
    ! max(xm2(i)*u_xmj(i)*dt/dx , 0.0) - min(xm2(i)*u_xmj(i-1)*dt/dx , 0.0) < 1
    !
    ! or equivalently:
    ! dt < dx / ( max(xm2(i)*u_xmj(i) , 0.0) - min(xm2(i)*u_xmj(i-1) , 0.0) )
    !
    ! In the case of variable cell size, dx is defined as dx(i) in these formula.
    !
    ! The value 1.e-30 is to ensure that we don't divide by 0 when all the velocities are 0.

    dth = dt_advec/GRIDWIDTH_M
    xcmax=0.0
    ycmax=0.0
    do k=1,KMAX_MID
       do j=1,ljmax
          xcmax(k,j+gj0-1) = maxval(                                         &
               max(u_xmj(1:limax  ,j,k,1)*xm2(1:limax,j),1.e-30)   &
               -min(u_xmj(0:limax-1,j,k,1)*xm2(1:limax,j),0.0   ))
       enddo
       do i=1,limax
          ycmax(k,i+gi0-1) = maxval(                                         &
               max(v_xmi(i,1:ljmax  ,k,1)*xm2(i,1:ljmax),1.e-30)   &
               -min(v_xmi(i,0:ljmax-1,k,1)*xm2(i,1:ljmax),0.0   ))
       enddo
    enddo

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,xcmax,KMAX_MID*gjmax,MPI_DOUBLE_PRECISION, &
         MPI_MAX,MPI_COMM_WORLD,INFO)

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,ycmax,KMAX_MID*gimax,MPI_DOUBLE_PRECISION, &
         MPI_MAX,MPI_COMM_WORLD, INFO)

    do i=1,limax
       do k=1,KMAX_MID
          dt_ymax(i,k)=GRIDWIDTH_M/ycmax(k,i+gi0-1)
       enddo
    enddo
    do j=1,ljmax
       do k=1,KMAX_MID
          dt_xmax(j,k)=GRIDWIDTH_M/xcmax(k,j+gj0-1)
       enddo
    enddo

    niterx=1
    do k=1,KMAX_MID
       do j=1,ljmax
          niterx(j,k) = int(dt_advec/dt_xmax(j,k))+1
          dt_x(j,k) = dt_advec/real(niterx(j,k))
          !if(me==0)write(*,*)'x',me,j,k,niterx(j,k),xcmax(k,j+gj0-1)
       enddo
    enddo

    do k=1,KMAX_MID
       do i=1,limax
          nitery(i,k) = int(dt_advec/dt_ymax(i,k))+1
          dt_y(i,k) = dt_advec/real(nitery(i,k))
       enddo
    enddo

    !Courant number in vertical sigma coordinates:  sigmadot*dt/deltasigma
    !
    !Note that dhs1(k+1) denotes thickness of layer k
    !     and sdot(k+1) denotes sdot at the boundary between layer k and k+1
    !
    !flux through wall k+1:  sdot(k+1) *dt/dhs1(k+1)<1   for sdot(k+1)>0
    !                       |sdot(k+1)|*dt/dhs1(k+2)<1   for sdot(k+1)<0
    !
    !layer k: sdot(k+1)*dt/dhs1(k+1) + |sdot(k)|*dt/dhs1(k+1) <1 for sdot(k+1)>0 and sdot(k)<0
    !
    !total out of layer k: max(sdot(1:limax,1:ljmax,k+1,1),0.0)-min(sdot(1:limax,1:ljmax,k,1),0.0)
    !
    scmax = 1.e-30
    do k = 1,KMAX_MID
       sdcmax = maxval(max(Etadot(1:limax,1:ljmax,k+1,1),0.0)   &
            -min(Etadot(1:limax,1:ljmax,k  ,1),0.0))
       scmax  = max(sdcmax/dhs1(k+1),scmax)
    enddo

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,scmax,1,MPI_DOUBLE_PRECISION, &
         MPI_MAX,MPI_COMM_WORLD,INFO)
    dt_smax = 1./scmax
42  FORMAT(A,F10.2)
    if(me==0.and. firstcall.and.DEBUG_ADV)write(*,42)'dt_smax',dt_smax
    niters = int(dt_advec/dt_smax)+1
    dt_s = dt_advec/real(niters)

    niterxys = 1

    nxy=0
    nxx=0
    nxxmin=0
    nyy=0
    do k=1,KMAX_MID
       do j=1,ljmax
          nxy=nxy+niterx(j,k)-1
          nxx=nxx+niterx(j,k)-1
          if(niterx(j,k)>NITERXMAX)then
             nxxmin=nxxmin+niterx(j,k)
          endif
       enddo
       do i=1,limax
          nxy=nxy+nitery(i,k)-1
          nyy=nyy+nitery(i,k)-1
       enddo

    enddo
    if(me.eq.0)then
       !          write(*,43)KMAX_MID*ljmax,nxx,nxxmin,KMAX_MID*limax,nyy,niters
    endif
    !43    format('total iterations x, y, k: ',I4,' +',I4,' -',I4,', ',I5,' +',I3,',',I4)

    ! stop

    call Add_2timing(20,tim_after,tim_before,"advecdiff:synchronization")

    ! Start xys advection loop:
    iterxys = 0
    do while (iterxys < niterxys)
       if(mod(nstep,2) /= 0 .or. iterxys /= 0)then !start a xys sequence

          iterxys = iterxys + 1
          do k = 1,KMAX_MID
             fac1=(dA(k)/Pref+dB(k))*xxdg
             do j = lj0,lj1
                if(niterx(j,k)<=NITERXMAX)then
                   dth = dt_x(j,k)/GRIDWIDTH_M
                   do iterx=1,niterx(j,k)

                      ! send/receive in x-direction
                      call preadvx2(110+k+KMAX_MID*j               &
                           ,xn_adv(1,1,j,k),dpdeta(1,j,k),u_xmj(0,j,k,1)&
                           ,xnw,xne                               &
                           ,psw,pse)

                      ! x-direction
                      call advx(                                   &
                           u_xmj(0,j,k,1),uw(j,k,1),ue(j,k,1)        &
                           ,xn_adv(1,1,j,k),xnw,xne               &
                           ,dpdeta(1,j,k),psw,pse                   &
                           ,xm2(0,j),xmd(0,j)                     &
                           ,dth,fac1)!carea(k))
                      do i = li0,li1
                         dpdeta0=(dA(k)+dB(k)*ps(i,j,1))/(dA(k)/Pref+dB(k))
                         psi = dpdeta0/max(dpdeta(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                         dpdeta(i,j,k) = dpdeta0
                      enddo

                   enddo !iter

                endif
             enddo !j
             !          enddo !k horizontal (x) advection


             call Add_2timing(21,tim_after,tim_before,"advecdiff:advx")

             ! y-direction
             !          do k = 1,KMAX_MID
             do i = li0,li1
                dth = dt_y(i,k)/GRIDWIDTH_M
                do itery=1,nitery(i,k)

                   ! send/receive in y-direction
                   call preadvy2(520+k                            &
                        ,xn_adv(1,1,1,k),dpdeta(1,1,k),v_xmi(1,0,k,1)    &
                        ,xns, xnn                                  &
                        ,pss, psn,i)

                   call advy(                                     &
                        v_xmi(i,0,k,1),vs(i,k,1),vn(i,k,1)            &
                        ,xn_adv(1,i,1,k),xns,xnn                   &
                        ,dpdeta(i,1,k),pss,psn                       &
                        ,xm2ji(0,i),xmdji(0,i)                     &
                        ,dth,fac1)!carea(k))
                   do j = lj0,lj1
                      dpdeta0=(dA(k)+dB(k)*ps(i,j,1))/(dA(k)/Pref+dB(k))
                      psi = dpdeta0/max(dpdeta(i,j,k),1.0)
                      xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      dpdeta(i,j,k) = dpdeta0
                   enddo

                enddo !iter
             enddo !i
          enddo !k horizontal (y) advection


          call Add_2timing(23,tim_after,tim_before,"advecdiff:advy")



          do iters=1,niters

             ! perform only vertical advection
             do j = lj0,lj1
                do i = li0,li1
                   call advvk(xn_adv(1,i,j,1),dpdeta(i,j,1),Etadot(i,j,1,1),dt_s)
                enddo
             enddo
             if(iters<niters .or. iterxys < niterxys)then
                do k=1,KMAX_MID
                   do j = lj0,lj1
                      do i = li0,li1
                         dpdeta0=(dA(k)+dB(k)*ps(i,j,1))/(dA(k)/Pref+dB(k))
                         psi = dpdeta0/max(dpdeta(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                         dpdeta(i,j,k) = dpdeta0
                      enddo
                   enddo
                enddo
             else
                !advection finished
                do k=1,KMAX_MID
                   do j = lj0,lj1
                      do i = li0,li1
                         psi =1.0/max(dpdeta(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      enddo
                   enddo
                enddo
             endif
          enddo ! vertical (s) advection


          call Add_2timing(24,tim_after,tim_before,"advecdiff:advvk")

       else  !start a yxs sequence

          iterxys = iterxys + 1
          do k = 1,KMAX_MID
             fac1=(dA(k)/Pref+dB(k))*xxdg
             do i = li0,li1
                dth = dt_y(i,k)/GRIDWIDTH_M
                do itery=1,nitery(i,k)

                   ! send/receive in y-direction
                   call preadvy2(13000+k+KMAX_MID*itery+1000*i    &
                        ,xn_adv(1,1,1,k),dpdeta(1,1,k),v_xmi(1,0,k,1)    &
                        ,xns, xnn                                  &
                        ,pss, psn,i)

                   ! y-direction
                   call advy(                                    &
                        v_xmi(i,0,k,1),vs(i,k,1),vn(i,k,1)           &
                        ,xn_adv(1,i,1,k),xns,xnn                  &
                        ,dpdeta(i,1,k),pss,psn                      &
                        ,xm2ji(0,i),xmdji(0,i)                    &
                        ,dth,fac1)!carea(k))

                   do j = lj0,lj1
                      dpdeta0=(dA(k)+dB(k)*ps(i,j,1))/(dA(k)/Pref+dB(k))
                      psi = dpdeta0/max(dpdeta(i,j,k),1.0)
                      xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      dpdeta(i,j,k) = dpdeta0
                   enddo

                enddo !iter
             enddo !i
             !         enddo !k horizontal (y) advection

 

             call Add_2timing(23,tim_after,tim_before,"advecdiff:preadvy,advy")

             !          do k = 1,KMAX_MID
             do j = lj0,lj1
                if(niterx(j,k)<=NITERXMAX)then
                   dth = dt_x(j,k)/GRIDWIDTH_M
                   do iterx=1,niterx(j,k)

                      ! send/receive in x-direction
                      call preadvx2(21000+k+KMAX_MID*iterx+1000*j  &
                           ,xn_adv(1,1,j,k),dpdeta(1,j,k),u_xmj(0,j,k,1)&
                           ,xnw,xne                               &
                           ,psw,pse)

                      ! x-direction
                      call advx(                                   &
                           u_xmj(0,j,k,1),uw(j,k,1),ue(j,k,1)        &
                           ,xn_adv(1,1,j,k),xnw,xne               &
                           ,dpdeta(1,j,k),psw,pse                   &
                           ,xm2(0,j),xmd(0,j)                     &
                           ,dth,fac1)!carea(k))

                      do i = li0,li1
                         dpdeta0=(dA(k)+dB(k)*ps(i,j,1))/(dA(k)/Pref+dB(k))
                         psi = dpdeta0/max(dpdeta(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                         dpdeta(i,j,k) = dpdeta0
                      enddo

                   enddo !iter
                endif
             enddo !j
          enddo !k horizontal (x) advection

          call Add_2timing(21,tim_after,tim_before,"advecdiff:preadvx,advx") 


          do iters=1,niters

             ! perform only vertical advection
             do j = lj0,lj1
                do i = li0,li1
                   call advvk(xn_adv(1,i,j,1),dpdeta(i,j,1),Etadot(i,j,1,1),dt_s)
                enddo
             enddo
             if(iters<niters .or. iterxys < niterxys)then
                do k=1,KMAX_MID
                   do j = lj0,lj1
                      do i = li0,li1
                         dpdeta0=(dA(k)+dB(k)*ps(i,j,1))/(dA(k)/Pref+dB(k))
                         psi = dpdeta0/max(dpdeta(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                         dpdeta(i,j,k) = dpdeta0
                      enddo
                   enddo
                enddo
             else
                !advection finished
                do k=1,KMAX_MID
                   do j = lj0,lj1
                      do i = li0,li1
                         psi =1.0/max(dpdeta(i,j,k),1.0)
                         xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
                      enddo
                   enddo
                enddo
             endif
          enddo ! vertical (s) advection


          call Add_2timing(24,tim_after,tim_before,"advecdiff:advvk")

       endif ! yxs sequence
    enddo



    if(USE_CONVECTION)then

       call CheckStop(ADVEC_TYPE/=1, "ADVEC_TYPE no longer supported")

       do k=1,KMAX_MID
          do j = lj0,lj1
             do i = li0,li1
                dpdeta(i,j,k) = (dA(k)+dB(k)*ps(i,j,1))/(dA(k)/Pref+dB(k))
                xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*dpdeta(i,j,k)
             enddo
          enddo
       enddo

       call convection_Eta(dpdeta,dt_advec)
       do k=1,KMAX_MID
          do j = lj0,lj1
             do i = li0,li1
                psi = 1.0/max(dpdeta(i,j,k),1.0)
                xn_adv(:,i,j,k) = xn_adv(:,i,j,k)*psi
             enddo
          enddo
          mindpdeta = minval( dpdeta(li0:li1, lj0:lj1, k) )
          if ( nWarnings < MAX_WARNINGS  .and.mindpdeta<1.0) then
             call datewrite("WARNING:C dpdeta < 1",k,  (/ mindpdeta /) )
             nWarnings = nWarnings + 1
          end if
       enddo

    endif


    do k=1,KMAX_MID
       do j = lj0,lj1
          if(niterx(j,k)>NITERXMAX)then
             !if(mex==0)write(*,*)'Simplified advection',k,j,niterx(j,k)
             !simplified "advection": average the mixing ratios over all x
             !average 1D
             xn_advjk=0.0
             isum=0
             do i = li0,li1
                !           psi=1.0/(ps(i,j,1)-PT)
                xn_advjk(:) = xn_advjk(:)+xn_adv(:,i,j,k)!*psi
                isum=isum+1
             enddo
             !sum over all processors along i direction. mex=0 collects the sum
             !me = mex + mey*NPROCX
             if(mex>0)then

                CALL MPI_SEND(xn_advjk,8*NSPEC_ADV,MPI_BYTE, &
                     mey*NPROCX,100*mey+j+1000,MPI_COMM_WORLD,INFO)
                !receive averages from mex=0
                CALL MPI_RECV(xn_advjk,8*NSPEC_ADV,MPI_BYTE, &
                     mey*NPROCX,100*mey+j+3000,MPI_COMM_WORLD,MPISTATUS,INFO)

             else
                xn_advjktot(:) = xn_advjk(:)
                isumtot=isum
                do iproc=1,NPROCX-1
                   CALL MPI_RECV(xn_advjk,8*NSPEC_ADV,MPI_BYTE, &
                        iproc+mey*NPROCX,100*mey+j+1000,MPI_COMM_WORLD,MPISTATUS,INFO)
                   xn_advjktot(:) = xn_advjktot(:)+xn_advjk(:)
                   !             isumtot=isumtot+isum
                enddo
                rfac=1.0/GIMAX
                xn_advjk(:) = xn_advjktot(:)*rfac
                !           write(*,*)'ISUM',mey,isumtot,isum,GIMAX
                !send result to all processors in i direction
                do iproc=1,NPROCX-1
                   CALL MPI_SEND(xn_advjk,8*NSPEC_ADV,MPI_BYTE, &
                        iproc+mey*NPROCX,100*mey+j+3000,MPI_COMM_WORLD,INFO)
                enddo
             endif

             do i = li0,li1
                xn_adv(:,i,j,k)= xn_advjk(:)
             enddo

          endif
       enddo
    enddo

    call Add_2timing(25,tim_after,tim_before,"advecdiff:ps")

    ! vertical diffusion
    ndiff = 1 !number of vertical diffusion iterations (the larger the better)
    do k = 2,KMAX_MID
       ds3(k) = dt_advec*dhs1i(k)*dhs2i(k)
       ds4(k) = dt_advec*dhs1i(k+1)*dhs2i(k)
    enddo

    ! sum is conserved under vertical diffusion
    !   sum = 0.
    !   do k=1,KMAX_MID
    !      sum = sum + xn_adv(1,4,4,k)/dhs1i(k+1)
    !   enddo
    !   write(*,*)'sum before diffusion ',me,sum
    do j = lj0,lj1
       do i = li0,li1
          call vertdiffn(xn_adv(1,i,j,1),EtaKz(i,j,1,1),ds3,ds4,ndiff)
       enddo
    enddo
    !   sum = 0.
    !   do k=1,KMAX_MID
    !      sum = sum + xn_adv(1,4,4,k)/dhs1i(k+1)
    !   enddo
    !   write(*,*)'sum after diffusion ',me,sum
    call Add_2timing(22,tim_after,tim_before,"advecdiff:diffusion")


    if(lj0.ne.1)then
       do k=KCHEMTOP,KMAX_MID
          do i = 1,limax
             xn_adv(:,i,1,k) = xn_adv(:,i,1,k)/((dA(k)+dB(k)*ps(i,1,1))/(dA(k)/Pref+dB(k)))
          enddo
       enddo
    endif
    if(li0.ne.1)then
       do k=KCHEMTOP,KMAX_MID
          do j=lj0,lj1
             xn_adv(:,1,j,k) = xn_adv(:,1,j,k)/((dA(k)+dB(k)*ps(1,j,1))/(dA(k)/Pref+dB(k)))
          enddo
       enddo
    endif
    if(li1.ne.limax)then
       do k=KCHEMTOP,KMAX_MID
          do j=lj0,lj1
             xn_adv(:,limax,j,k) = xn_adv(:,limax,j,k)/((dA(k)+dB(k)*ps(limax,j,1))/(dA(k)/Pref+dB(k)))
          enddo
       enddo
    endif
    if(lj1.ne.ljmax)then
       do k=KCHEMTOP,KMAX_MID
          do i = 1,limax
             xn_adv(:,i,ljmax,k) = xn_adv(:,i,ljmax,k)/((dA(k)+dB(k)*ps(i,ljmax,1))/(dA(k)/Pref+dB(k)))
          enddo
       enddo
    endif

    if(KCHEMTOP==2)then

       !pw since the xn_adv are changed it corresponds to a flux in or
       !   out of the system:
       do i = li0,li1
          do j = lj0,lj1
             where(xn_adv(:,i,j,1) .gt. xntop(:,i,j))
                fluxout(:) = fluxout(:) + &
                     (xn_adv(:,i,j,1)-xntop(:,i,j))*(dA(1)+dB(1)*ps(i,j,1))*xxdg*xmd(i,j)
             elsewhere
                fluxin(:) = fluxin(:) + &
                     (xntop(:,i,j)-xn_adv(:,i,j,1))*(dA(1)+dB(1)*ps(i,j,1))*xxdg*xmd(i,j)
             end where
          enddo
       enddo
       xn_adv(:,:,:,1) = xntop(:,:,:)

    endif

    firstcall=.false.
    return

  end subroutine advecdiff_Eta

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine vgrid
!
!     inclusion of the variable grid spacing when interpolating the
!     polynominal is done by introducing new local coordinates, cor.
!
!
!     modified by pw january 2002: alfnew is modified such that
!     a Courant number of one corresponds exactly to "empty" a cell.
!     (small effects on results: less than 1%)
!


    use GridValues_ml, only : sigma_bnd,sigma_mid
    implicit none

    integer k
    real cor1, cor2, dcorl
    real hscor1(KMAX_BND+2),hscor2(KMAX_MID+2)
    real alfa1(NADVS), alfa2(NADVS), dei
    real corl1(KMAX_BND), corl2(KMAX_BND)

    real alf(9,2:KMAX_BND)

    do  k=1,KMAX_MID
      hscor1(k+1) = sigma_bnd(k)
      hscor2(k+1) = sigma_mid(k)
    enddo
    hscor1(KMAX_BND+1) = sigma_bnd(KMAX_BND)

    hscor1(1) = - sigma_bnd(2)
    hscor1(KMAX_BND+2) = 2.*sigma_bnd(KMAX_BND) - sigma_bnd(KMAX_BND-1)
    hscor2(1) = 2.*sigma_mid(1) - sigma_mid(2)
    hscor2(KMAX_MID+2) = 2.*sigma_mid(KMAX_MID) - sigma_mid(KMAX_MID-1)

    do  k=1,KMAX_BND
      dhs1(k) = hscor1(k+1) - hscor1(k)
      dhs1i(k) = 1./dhs1(k)
      dhs2i(k) = 1./(hscor2(k+1) - hscor2(k))
    enddo

    do k=2,KMAX_BND

      corl1(k) = (hscor1(k) - hscor2(k))*dhs1i(k)
      corl2(k) = (hscor1(k+1) - hscor2(k))*dhs1i(k)
      dcorl = corl2(k) - corl1(k)
      alfa1(NADVS-1) = (corl2(k)**2 - corl1(k)**2)/(2.*dcorl)
      alfa2(NADVS-1) = (corl2(k)**3 - corl1(k)**3)/(3.*dcorl)

      cor1 = (hscor1(k-1) - hscor2(k))*dhs1i(k)
!     cor2 = (hscor1(k) - hscor2(k))*dhs1i(k)
      dcorl = corl1(k) - cor1
      alfa1(NADVS-2) = (corl1(k)**2 - cor1**2)/(2.*dcorl)
      alfa2(NADVS-2) = (corl1(k)**3 - cor1**3)/(3.*dcorl)

!     cor1 = (hscor1(k+1) - hscor2(k))*dhs1i(k)
      cor2 = (hscor1(k+2) - hscor2(k))*dhs1i(k)
      dcorl = cor2 - corl2(k)
      alfa1(NADVS) = (cor2**2 - corl2(k)**2)/(2.*dcorl)
      alfa2(NADVS) = (cor2**3 - corl2(k)**3)/(3.*dcorl)

      dei = alfa1(NADVS-1)*alfa2(NADVS)                  &
          - alfa1(NADVS)  *alfa2(NADVS-1)                &
          - alfa1(NADVS-2)*(alfa2(NADVS)-alfa2(NADVS-1)) &
          + alfa2(NADVS-2)*(alfa1(NADVS)-alfa1(NADVS-1))
      dei = 1./dei

      alf(1,k) = dei*(alfa1(NADVS-1)*alfa2(NADVS)        &
                     -alfa1(NADVS)  *alfa2(NADVS-1))
      alf(4,k) = dei*(alfa2(NADVS-2)*alfa1(NADVS)        &
                     -alfa2(NADVS)  *alfa1(NADVS-2))
      alf(7,k) = dei*(alfa2(NADVS-1)*alfa1(NADVS-2)      &
                     -alfa2(NADVS-2)*alfa1(NADVS-1))
      alf(2,k) = dei*(alfa2(NADVS-1)-alfa2(NADVS))  /2.
      alf(5,k) = dei*(alfa2(NADVS)  -alfa2(NADVS-2))/2.
      alf(8,k) = dei*(alfa2(NADVS-2)-alfa2(NADVS-1))/2.
      alf(3,k) = dei*(alfa1(NADVS)  -alfa1(NADVS-1))/3.
      alf(6,k) = dei*(alfa1(NADVS-2)-alfa1(NADVS))  /3.
      alf(9,k) = dei*(alfa1(NADVS-1)-alfa1(NADVS-2))/3.
    enddo

    do k=2,KMAX_MID
      alfnew(1,k,0) = alf(1,k) + 2.*alf(2,k)*corl2(k)                &
                    + 3.*alf(3,k)*corl2(k)*corl2(k)
      alfnew(1,k,1) = -(alf(1,k+1) + 2.*alf(2,k+1)*corl1(k+1)        &
                    + 3.*alf(3,k+1)*corl1(k+1)*corl1(k+1))
      alfnew(4,k,0) = alf(4,k) + 2.*alf(5,k)*corl2(k)                &
                    + 3.*alf(6,k)*corl2(k)*corl2(k)
      alfnew(4,k,1) = -(alf(4,k+1) + 2.*alf(5,k+1)*corl1(k+1)        &
                    + 3.*alf(6,k+1)*corl1(k+1)*corl1(k+1))
      alfnew(7,k,0) = alf(7,k) + 2.*alf(8,k)*corl2(k)                &
                    + 3.*alf(9,k)*corl2(k)*corl2(k)
      alfnew(7,k,1) = -(alf(7,k+1) + 2.*alf(8,k+1)*corl1(k+1)        &
                    + 3.*alf(9,k+1)*corl1(k+1)*corl1(k+1))
!pw   alfnew(2,k,0) = -(alf(2,k)   + 3.*alf(3,k)  *corl2(k))  *dhs2i(k)
!     alfnew(2,k,1) =  (alf(2,k+1) + 3.*alf(3,k+1)*corl1(k+1))*dhs2i(k)
      alfnew(2,k,0) = -(alf(2,k)   + 3.*alf(3,k)  *corl2(k))  *dhs1i(k)
      alfnew(2,k,1) =  (alf(2,k+1) + 3.*alf(3,k+1)*corl1(k+1))*dhs1i(k+1)
!pw   alfnew(5,k,0) = -(alf(5,k)   + 3.*alf(6,k)  *corl2(k))  *dhs2i(k)
!     alfnew(5,k,1) =  (alf(5,k+1) + 3.*alf(6,k+1)*corl1(k+1))*dhs2i(k)
      alfnew(5,k,0) = -(alf(5,k)   + 3.*alf(6,k)  *corl2(k))  *dhs1i(k)
      alfnew(5,k,1) =  (alf(5,k+1) + 3.*alf(6,k+1)*corl1(k+1))*dhs1i(k+1)
!pw   alfnew(8,k,0) = -(alf(8,k)   + 3.*alf(9,k)  *corl2(k))  *dhs2i(k)
!     alfnew(8,k,1) = (alf(8,k+1) + 3.*alf(9,k+1)*corl1(k+1))*dhs2i(k)
      alfnew(8,k,0) = -(alf(8,k)   + 3.*alf(9,k)  *corl2(k))  *dhs1i(k)
      alfnew(8,k,1) =  (alf(8,k+1) + 3.*alf(9,k+1)*corl1(k+1))*dhs1i(k+1)
!pw   alfnew(3,k,0) =  alf(3,k)  *dhs2i(k)*dhs2i(k)
!     alfnew(3,k,1) = -alf(3,k+1)*dhs2i(k)*dhs2i(k)
!     alfnew(6,k,0) =  alf(6,k)  *dhs2i(k)*dhs2i(k)
!     alfnew(6,k,1) = -alf(6,k+1)*dhs2i(k)*dhs2i(k)
!     alfnew(9,k,0) =  alf(9,k)  *dhs2i(k)*dhs2i(k)
!     alfnew(9,k,1) = -alf(9,k+1)*dhs2i(k)*dhs2i(k)
      alfnew(3,k,0) =  alf(3,k)  *dhs1i(k)  *dhs1i(k)
      alfnew(3,k,1) = -alf(3,k+1)*dhs1i(k+1)*dhs1i(k+1)
      alfnew(6,k,0) =  alf(6,k)  *dhs1i(k)  *dhs1i(k)
      alfnew(6,k,1) = -alf(6,k+1)*dhs1i(k+1)*dhs1i(k+1)
      alfnew(9,k,0) =  alf(9,k)  *dhs1i(k)  *dhs1i(k)
      alfnew(9,k,1) = -alf(9,k+1)*dhs1i(k+1)*dhs1i(k+1)
    enddo

    alfbegnew(1) = alfnew(1,2,0)+alfnew(4,2,0)
    alfbegnew(2) = alfnew(2,2,0)+alfnew(5,2,0)
    alfbegnew(3) = alfnew(3,2,0)+alfnew(6,2,0)
    alfendnew(1) = alfnew(4,KMAX_MID,1)+alfnew(7,KMAX_MID,1)
    alfendnew(2) = alfnew(5,KMAX_MID,1)+alfnew(8,KMAX_MID,1)
    alfendnew(3) = alfnew(6,KMAX_MID,1)+alfnew(9,KMAX_MID,1)

  end subroutine vgrid

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine vgrid_Eta
!
!     inclusion of the variable grid spacing when interpolating the
!     polynominal is done by introducing new local coordinates, cor.
!
!
!     modified by pw january 2002: alfnew is modified such that
!     a Courant number of one corresponds exactly to "empty" a cell.
!     (small effects on results: less than 1%)
!
!     Adapted for pressure coordinates


    use GridValues_ml, only : Eta_bnd,Eta_mid
    implicit none

    integer k
    real cor1, cor2, dcorl
    real hscor1(KMAX_BND+2),hscor2(KMAX_MID+2)
    real alfa1(NADVS), alfa2(NADVS), dei
    real corl1(KMAX_BND), corl2(KMAX_BND)

    real alf(9,2:KMAX_BND)

    do  k=1,KMAX_MID
      hscor1(k+1) = Eta_bnd(k)
      hscor2(k+1) = Eta_mid(k)
    enddo
    hscor1(KMAX_BND+1) = Eta_bnd(KMAX_BND)

    hscor1(1) = - Eta_bnd(2)
    hscor1(KMAX_BND+2) = 2.*Eta_bnd(KMAX_BND) - Eta_bnd(KMAX_BND-1)
    hscor2(1) = 2.*Eta_mid(1) - Eta_mid(2)
    hscor2(KMAX_MID+2) = 2.*Eta_mid(KMAX_MID) - Eta_mid(KMAX_MID-1)

    do  k=1,KMAX_BND
      dhs1(k) = hscor1(k+1) - hscor1(k)
      dhs1i(k) = 1./dhs1(k)
      dhs2i(k) = 1./(hscor2(k+1) - hscor2(k))
    enddo

    do k=2,KMAX_BND

      corl1(k) = (hscor1(k) - hscor2(k))*dhs1i(k)
      corl2(k) = (hscor1(k+1) - hscor2(k))*dhs1i(k)
      dcorl = corl2(k) - corl1(k)
      alfa1(NADVS-1) = (corl2(k)**2 - corl1(k)**2)/(2.*dcorl)
      alfa2(NADVS-1) = (corl2(k)**3 - corl1(k)**3)/(3.*dcorl)

      cor1 = (hscor1(k-1) - hscor2(k))*dhs1i(k)
!     cor2 = (hscor1(k) - hscor2(k))*dhs1i(k)
      dcorl = corl1(k) - cor1
      alfa1(NADVS-2) = (corl1(k)**2 - cor1**2)/(2.*dcorl)
      alfa2(NADVS-2) = (corl1(k)**3 - cor1**3)/(3.*dcorl)

!     cor1 = (hscor1(k+1) - hscor2(k))*dhs1i(k)
      cor2 = (hscor1(k+2) - hscor2(k))*dhs1i(k)
      dcorl = cor2 - corl2(k)
      alfa1(NADVS) = (cor2**2 - corl2(k)**2)/(2.*dcorl)
      alfa2(NADVS) = (cor2**3 - corl2(k)**3)/(3.*dcorl)

      dei = alfa1(NADVS-1)*alfa2(NADVS)                  &
          - alfa1(NADVS)  *alfa2(NADVS-1)                &
          - alfa1(NADVS-2)*(alfa2(NADVS)-alfa2(NADVS-1)) &
          + alfa2(NADVS-2)*(alfa1(NADVS)-alfa1(NADVS-1))
      dei = 1./dei

      alf(1,k) = dei*(alfa1(NADVS-1)*alfa2(NADVS)        &
                     -alfa1(NADVS)  *alfa2(NADVS-1))
      alf(4,k) = dei*(alfa2(NADVS-2)*alfa1(NADVS)        &
                     -alfa2(NADVS)  *alfa1(NADVS-2))
      alf(7,k) = dei*(alfa2(NADVS-1)*alfa1(NADVS-2)      &
                     -alfa2(NADVS-2)*alfa1(NADVS-1))
      alf(2,k) = dei*(alfa2(NADVS-1)-alfa2(NADVS))  /2.
      alf(5,k) = dei*(alfa2(NADVS)  -alfa2(NADVS-2))/2.
      alf(8,k) = dei*(alfa2(NADVS-2)-alfa2(NADVS-1))/2.
      alf(3,k) = dei*(alfa1(NADVS)  -alfa1(NADVS-1))/3.
      alf(6,k) = dei*(alfa1(NADVS-2)-alfa1(NADVS))  /3.
      alf(9,k) = dei*(alfa1(NADVS-1)-alfa1(NADVS-2))/3.
    enddo

    do k=2,KMAX_MID
      alfnew(1,k,0) = alf(1,k) + 2.*alf(2,k)*corl2(k)                &
                    + 3.*alf(3,k)*corl2(k)*corl2(k)
      alfnew(1,k,1) = -(alf(1,k+1) + 2.*alf(2,k+1)*corl1(k+1)        &
                    + 3.*alf(3,k+1)*corl1(k+1)*corl1(k+1))
      alfnew(4,k,0) = alf(4,k) + 2.*alf(5,k)*corl2(k)                &
                    + 3.*alf(6,k)*corl2(k)*corl2(k)
      alfnew(4,k,1) = -(alf(4,k+1) + 2.*alf(5,k+1)*corl1(k+1)        &
                    + 3.*alf(6,k+1)*corl1(k+1)*corl1(k+1))
      alfnew(7,k,0) = alf(7,k) + 2.*alf(8,k)*corl2(k)                &
                    + 3.*alf(9,k)*corl2(k)*corl2(k)
      alfnew(7,k,1) = -(alf(7,k+1) + 2.*alf(8,k+1)*corl1(k+1)        &
                    + 3.*alf(9,k+1)*corl1(k+1)*corl1(k+1))
!pw   alfnew(2,k,0) = -(alf(2,k)   + 3.*alf(3,k)  *corl2(k))  *dhs2i(k)
!     alfnew(2,k,1) =  (alf(2,k+1) + 3.*alf(3,k+1)*corl1(k+1))*dhs2i(k)
      alfnew(2,k,0) = -(alf(2,k)   + 3.*alf(3,k)  *corl2(k))  *dhs1i(k)
      alfnew(2,k,1) =  (alf(2,k+1) + 3.*alf(3,k+1)*corl1(k+1))*dhs1i(k+1)
!pw   alfnew(5,k,0) = -(alf(5,k)   + 3.*alf(6,k)  *corl2(k))  *dhs2i(k)
!     alfnew(5,k,1) =  (alf(5,k+1) + 3.*alf(6,k+1)*corl1(k+1))*dhs2i(k)
      alfnew(5,k,0) = -(alf(5,k)   + 3.*alf(6,k)  *corl2(k))  *dhs1i(k)
      alfnew(5,k,1) =  (alf(5,k+1) + 3.*alf(6,k+1)*corl1(k+1))*dhs1i(k+1)
!pw   alfnew(8,k,0) = -(alf(8,k)   + 3.*alf(9,k)  *corl2(k))  *dhs2i(k)
!     alfnew(8,k,1) = (alf(8,k+1) + 3.*alf(9,k+1)*corl1(k+1))*dhs2i(k)
      alfnew(8,k,0) = -(alf(8,k)   + 3.*alf(9,k)  *corl2(k))  *dhs1i(k)
      alfnew(8,k,1) =  (alf(8,k+1) + 3.*alf(9,k+1)*corl1(k+1))*dhs1i(k+1)
!pw   alfnew(3,k,0) =  alf(3,k)  *dhs2i(k)*dhs2i(k)
!     alfnew(3,k,1) = -alf(3,k+1)*dhs2i(k)*dhs2i(k)
!     alfnew(6,k,0) =  alf(6,k)  *dhs2i(k)*dhs2i(k)
!     alfnew(6,k,1) = -alf(6,k+1)*dhs2i(k)*dhs2i(k)
!     alfnew(9,k,0) =  alf(9,k)  *dhs2i(k)*dhs2i(k)
!     alfnew(9,k,1) = -alf(9,k+1)*dhs2i(k)*dhs2i(k)
      alfnew(3,k,0) =  alf(3,k)  *dhs1i(k)  *dhs1i(k)
      alfnew(3,k,1) = -alf(3,k+1)*dhs1i(k+1)*dhs1i(k+1)
      alfnew(6,k,0) =  alf(6,k)  *dhs1i(k)  *dhs1i(k)
      alfnew(6,k,1) = -alf(6,k+1)*dhs1i(k+1)*dhs1i(k+1)
      alfnew(9,k,0) =  alf(9,k)  *dhs1i(k)  *dhs1i(k)
      alfnew(9,k,1) = -alf(9,k+1)*dhs1i(k+1)*dhs1i(k+1)
    enddo

    alfbegnew(1) = alfnew(1,2,0)+alfnew(4,2,0)
    alfbegnew(2) = alfnew(2,2,0)+alfnew(5,2,0)
    alfbegnew(3) = alfnew(3,2,0)+alfnew(6,2,0)
    alfendnew(1) = alfnew(4,KMAX_MID,1)+alfnew(7,KMAX_MID,1)
    alfendnew(2) = alfnew(5,KMAX_MID,1)+alfnew(8,KMAX_MID,1)
    alfendnew(3) = alfnew(6,KMAX_MID,1)+alfnew(9,KMAX_MID,1)

  end subroutine vgrid_Eta

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine advvk(xn_adv,ps3d,sdot,dt_s)

!     executes advection with a. bott's integreated flux-form
!     using 2'nd order polynomial in the vertical.

    use ModelConstants_ml   , only : EPSIL, dt_advec
    use ChemSpecs,         only : NSPEC_ADV
    implicit none

!    input
    real,intent(in)::  sdot(0:MAXLIMAX*MAXLJMAX*KMAX_BND-1),dt_s

!    input+output
    real ,intent(inout):: xn_adv(NSPEC_ADV,0:MAXLIMAX*MAXLJMAX*KMAX_MID-1)
    real ,intent(inout):: ps3d(0:MAXLIMAX*MAXLJMAX*KMAX_MID-1)

!    local
    real fluxk(NSPEC_ADV,KMAX_MID),fluxps(KMAX_MID),fc(KMAX_MID)

!    local
    integer  k, n1k,k1
    integer klimlow,klimhig
    real zzfl1,zzfl2,zzfl3,totk(NSPEC_ADV),totps
    real fc1,fc2,fc3

    do k = 1,KMAX_MID-1
      fc(k) = sdot(k*MAXLIMAX*MAXLJMAX)*dt_s
    enddo

    fc(KMAX_MID) = -1.
!-------------- calculate the advection ----------------------------

    klimlow = 1
    if(fc(1).ge.0.)klimlow=2
      klimhig = KMAX_MID-1
      if(fc(KMAX_MID-1).lt.0.)klimhig = KMAX_MID-2

        fluxk(:,1) = 0.
        fluxps(1) = 0.

        if(fc(1).ge.0.)then

          fc1 = fc(1)
          fc2 = fc1*fc1
          fc3 = fc1*fc2
          zzfl2 = alfbegnew(1)*fc1            &
                + alfbegnew(2)*fc2            &
                + alfbegnew(3)*fc3
          zzfl3 = alfnew(7,2,0)*fc1           &
                + alfnew(8,2,0)*fc2           &
                + alfnew(9,2,0)*fc3

          fluxk(:,2) = max(0.,xn_adv(:,0)*zzfl2    &
               +xn_adv(:,MAXLIMAX*MAXLJMAX)*zzfl3)
          fluxps(2)  = max(0.,ps3d(0)*zzfl2        &
               +ps3d(MAXLIMAX*MAXLJMAX)*zzfl3)

        endif
        do k = klimlow,klimhig

        fc1 = fc(k)
        fc2 = fc1*fc1
        fc3 = fc1*fc2
        n1k = 0
        if(fc1.lt.0)n1k=1
!pw bug corrected 29/8-2002 (emep1.2beta):
!       zzfl1 = alfnew(1,k,n1k)*fc1           &
!             + alfnew(2,k,n1k)*fc2           &
!             + alfnew(3,k,n1k)*fc3
!       zzfl2 = alfnew(4,k,n1k)*fc1           &
!             + alfnew(5,k,n1k)*fc2           &
!             + alfnew(6,k,n1k)*fc3
!       zzfl3 = alfnew(7,k,n1k)*fc1           &
!             + alfnew(8,k,n1k)*fc2           &
!             + alfnew(9,k,n1k)*fc3
        zzfl1 = alfnew(1,k+1,n1k)*fc1         &
              + alfnew(2,k+1,n1k)*fc2         &
              + alfnew(3,k+1,n1k)*fc3
        zzfl2 = alfnew(4,k+1,n1k)*fc1         &
              + alfnew(5,k+1,n1k)*fc2         &
              + alfnew(6,k+1,n1k)*fc3
        zzfl3 = alfnew(7,k+1,n1k)*fc1         &
              + alfnew(8,k+1,n1k)*fc2         &
              + alfnew(9,k+1,n1k)*fc3

        k1 = k-1+n1k

        fluxk(:,k+1) = max(0.,                            &
               xn_adv(:,(k1-1)*MAXLIMAX*MAXLJMAX)*zzfl1   &
              +xn_adv(:, k1   *MAXLIMAX*MAXLJMAX)*zzfl2   &
              +xn_adv(:,(k1+1)*MAXLIMAX*MAXLJMAX)*zzfl3)
        fluxps(k+1) = max(0.,                             &
               ps3d((k1-1)*MAXLIMAX*MAXLJMAX)*zzfl1       &
              +ps3d( k1   *MAXLIMAX*MAXLJMAX)*zzfl2       &
              +ps3d((k1+1)*MAXLIMAX*MAXLJMAX)*zzfl3)

    enddo
      if(fc(KMAX_MID-1).lt.0.)then

        fc1 = fc(KMAX_MID-1)
        fc2 = fc1*fc1
        fc3 = fc1*fc2
        zzfl1 = alfnew(1,KMAX_MID,1)*fc1      &
              + alfnew(2,KMAX_MID,1)*fc2      &
              + alfnew(3,KMAX_MID,1)*fc3
        zzfl2 = alfendnew(1)*fc1              &
              + alfendnew(2)*fc2              &
              + alfendnew(3)*fc3

        fluxk(:,KMAX_MID) =                                          &
            max(0.,xn_adv(:,(KMAX_MID-2)*MAXLIMAX*MAXLJMAX)*zzfl1    &
                  +xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*zzfl2)
        fluxps(KMAX_MID) =                                           &
            max(0.,ps3d((KMAX_MID-2)*MAXLIMAX*MAXLJMAX)*zzfl1        &
                  +ps3d((KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*zzfl2)

    endif

    k=1
    do while(k.lt.KMAX_MID)
      if(fc(k).lt.0.) then
        if(fc(k+1).ge.0.) then
          totk(:) = min(xn_adv(:,k*MAXLIMAX*MAXLJMAX)*dhs1(k+2)       &
                      /(fluxk(:,k+1) + fluxk(:,k+2)+ EPSIL),1.)
          fluxk(:,k+1) = -fluxk(:,k+1)*totk(:)
          fluxk(:,k+2) =  fluxk(:,k+2)*totk(:)
          xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) =                         &
                 max(0.,xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)             &
                      -(fluxk(:,k+1) - fluxk(:,k))*dhs1i(k+1))
          xn_adv(:, k   *MAXLIMAX*MAXLJMAX) =                         &
                 max(0.,xn_adv(:,k*MAXLIMAX*MAXLJMAX)                 &
                      -(fluxk(:,k+2) - fluxk(:,k+1))*dhs1i(k+2))

          totps = min(ps3d(k*MAXLIMAX*MAXLJMAX)*dhs1(k+2)             &
                    /(fluxps(k+1) + fluxps(k+2)+ EPSIL),1.)
          fluxps(k+1) = -fluxps(k+1)*totps
          fluxps(k+2) =  fluxps(k+2)*totps
          ps3d((k-1)*MAXLIMAX*MAXLJMAX) =                             &
               max(0.,ps3d((k-1)*MAXLIMAX*MAXLJMAX)                   &
                    -(fluxps(k+1) - fluxps(k))*dhs1i(k+1))
          ps3d( k   *MAXLIMAX*MAXLJMAX) =                             &
               max(0.,ps3d(k*MAXLIMAX*MAXLJMAX)                       &
                    -(fluxps(k+2) - fluxps(k+1))*dhs1i(k+2))
          k = k+2
        else
          fluxk(:,k+1) =                                                 &
              -min(xn_adv(:,k*MAXLIMAX*MAXLJMAX)*dhs1(k+2),fluxk(:,k+1))
          xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) =                            &
               max(0.,xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)                  &
                    -(fluxk(:,k+1) - fluxk(:,k))*dhs1i(k+1))
          fluxps(k+1) =                                                  &
              -min(ps3d(k*MAXLIMAX*MAXLJMAX)*dhs1(k+2),fluxps(k+1))
          ps3d((k-1)*MAXLIMAX*MAXLJMAX) =                                &
               max(0.,ps3d((k-1)*MAXLIMAX*MAXLJMAX)                      &
                    -(fluxps(k+1) - fluxps(k))*dhs1i(k+1))
          k = k+1
        endif
      else
        fluxk(:,k+1) =                                                   &
            min(xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)*dhs1(k+1),fluxk(:,k+1))
        xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) =                              &
            max(0.,xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)                     &
                 -(fluxk(:,k+1) - fluxk(:,k))*dhs1i(k+1))
        fluxps(k+1) =                                                    &
            min(ps3d((k-1)*MAXLIMAX*MAXLJMAX)*dhs1(k+1),fluxps(k+1))
        ps3d((k-1)*MAXLIMAX*MAXLJMAX) =                                  &
            max(0.,ps3d((k-1)*MAXLIMAX*MAXLJMAX)                         &
                 -(fluxps(k+1) - fluxps(k))*dhs1i(k+1))
        k = k+1
      endif
    enddo

    xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX) =               &
            max(0.,xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)  &
                  +fluxk(:,KMAX_MID)*dhs1i(KMAX_MID+1))
    ps3d((KMAX_MID-1)*MAXLIMAX*MAXLJMAX) =                   &
            max(0.,ps3d((KMAX_MID-1)*MAXLIMAX*MAXLJMAX)      &
                  +fluxps(KMAX_MID)*dhs1i(KMAX_MID+1))

  end subroutine advvk

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine vertdiff(xn_adv,SigmaKz,ds3,ds4)

!     executes vertical diffusion

    use ModelConstants_ml  , only : KCHEMTOP, EPSIL
    use ChemSpecs,         only : NSPEC_ADV

    implicit none

!    input
    real,intent(in)::  SigmaKz(0:MAXLIMAX*MAXLJMAX*KMAX_BND-1)
    real,intent(in)::  ds3(KMAX_MID-1),ds4(KMAX_MID-1)

!    output
    real ,intent(inout):: xn_adv(NSPEC_ADV,0:MAXLIMAX*MAXLJMAX*(KMAX_MID-1))

!    local

    integer  k
    real, dimension(KMAX_MID) :: adif,bdif,cdif,e1

    do k = 1,KMAX_MID-1
      adif(k) = SigmaKz(k*MAXLIMAX*MAXLJMAX)*ds3(k)
      bdif(k+1) = SigmaKz(k*MAXLIMAX*MAXLJMAX)*ds4(k)
    enddo

    cdif(KMAX_MID) = 1./(1. + bdif(KMAX_MID))
    e1(KMAX_MID) = bdif(KMAX_MID)*cdif(KMAX_MID)
    xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX) = &
      xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*cdif(KMAX_MID)

    do k = KMAX_MID-1,2,-1
      cdif(k) = 1./(1. + bdif(k) + adif(k) - adif(k)*e1(k+1))
      e1(k) = bdif(k)*cdif(k)
      xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) =                 &
          (xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)              &
         +adif(k)*xn_adv(:,(k)*MAXLIMAX*MAXLJMAX))*cdif(k)
    enddo

    cdif(1) = 1./(1. + adif(1) - adif(1)*e1(2))
    xn_adv(:,0) = (xn_adv(:,0) + adif(1)*xn_adv(:,MAXLIMAX*MAXLJMAX))*cdif(1)

    do k = 2,KMAX_MID
      xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) =                &
          e1(k)*xn_adv(:,(k-2)*MAXLIMAX*MAXLJMAX)        &
         +xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)
    enddo

  end subroutine vertdiff

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine vertdiffn(xn_adv,SigmaKz,ds3,ds4,ndiff)

!     executes vertical diffusion ndiff times

! SigmaKz(k) mixes xn_adv(k) and xn_adv(k-1)
!
! adif(k) -> mixing of layers k and k+1:
!            SigmaKz(k+1)*ds3(k+1)= SigmaKz(k+1)*dt_advec*dhs1i(k+1)*dhs2i(k+1)
!            = SigmaKz(k+1)*dt_advec/(sigma_bnd(k+1)-sigma_bnd(k))/(sigma_mid(k+1)-sigma_mid(k))
!
! bdif(k) -> mixing of layers k and k-1:
!            SigmaKz(k)*ds4(k)= SigmaKz(k)*dt_advec*dhs1i(k+1)*dhs2i(k)
!            = SigmaKz(k+1)*dt_advec/(sigma_bnd(k+1)-sigma_bnd(k))/(sigma_mid(k)-sigma_mid(k-1))

    use ModelConstants_ml  , only : KCHEMTOP, EPSIL
    use ChemSpecs,         only : NSPEC_ADV

    implicit none

!    input
    real,intent(in)::  SigmaKz(0:MAXLIMAX*MAXLJMAX*KMAX_BND-1)
    real,intent(in)::  ds3(KMAX_MID-1),ds4(KMAX_MID-1)
    integer,intent(in)::  ndiff

!    output
    real ,intent(inout):: xn_adv(NSPEC_ADV,0:MAXLIMAX*MAXLJMAX*(KMAX_MID-1))

!    local

    integer  k,n

    real, dimension(0:KMAX_MID-1) :: adif,bdif,cdif,e1  

    real ndiffi

    ndiffi=1./ndiff

    do k = 1,KMAX_MID-1
      adif(k-1) = SigmaKz(k*MAXLIMAX*MAXLJMAX)*ds3(k)*ndiffi
      bdif(k) = SigmaKz(k*MAXLIMAX*MAXLJMAX)*ds4(k)*ndiffi
    enddo

    cdif(KMAX_MID-1) = 1./(1. + bdif(KMAX_MID-1))
    e1(KMAX_MID-1) = bdif(KMAX_MID-1)*cdif(KMAX_MID-1)

    do k = KMAX_MID-2,1,-1
      cdif(k) = 1./(1. + bdif(k) + adif(k) - adif(k)*e1(k+1))
      e1(k) = bdif(k)*cdif(k)
    enddo

    cdif(0) = 1./(1. + adif(0) - adif(0)*e1(1))

    do n=1,ndiff
      xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX) = &
         xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*cdif(KMAX_MID-1)
      do k = KMAX_MID-2,0,-1
         xn_adv(:,k*MAXLIMAX*MAXLJMAX) =         &
           (xn_adv(:,k*MAXLIMAX*MAXLJMAX)        &
           +adif(k)*xn_adv(:,(k+1)*MAXLIMAX*MAXLJMAX))*cdif(k)
      enddo

      do k = 1,KMAX_MID-1
         xn_adv(:,k*MAXLIMAX*MAXLJMAX) =            &
            e1(k)*xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX) &
           +xn_adv(:,k*MAXLIMAX*MAXLJMAX)
      enddo

    enddo ! ndiff

  end subroutine vertdiffn

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine vertdiffn2(xn_adv,SigmaKz,ds3,ds4,ndiff)

!     executes vertical diffusion ndiff times

    use ModelConstants_ml  , only : KCHEMTOP, EPSIL
    use ChemSpecs,         only : NSPEC_ADV

    implicit none

!    input
    real,intent(in)::  SigmaKz(0:MAXLIMAX*MAXLJMAX*KMAX_BND-1)
    real,intent(in)::  ds3(KMAX_MID-1),ds4(KMAX_MID-1)
    integer,intent(in)::  ndiff

!    output
    real ,intent(inout):: xn_adv(NSPEC_ADV,0:MAXLIMAX*MAXLJMAX*(KMAX_MID-1))

!    local

    integer  k,n

    real, dimension(KMAX_MID) :: adif,bdif,cdif,e1 

    real ndiffi

    ndiffi=1./ndiff

    do k = 1,KMAX_MID-1
      adif(k) = SigmaKz(k*MAXLIMAX*MAXLJMAX)*ds3(k)*ndiffi
      bdif(k+1) = SigmaKz(k*MAXLIMAX*MAXLJMAX)*ds4(k)*ndiffi
    enddo

    cdif(KMAX_MID) = 1./(1. + bdif(KMAX_MID))
    e1(KMAX_MID) = bdif(KMAX_MID)*cdif(KMAX_MID)

    do k = KMAX_MID-1,2,-1
      cdif(k) = 1./(1. + bdif(k) + adif(k) - adif(k)*e1(k+1))
      e1(k) = bdif(k)*cdif(k)
    enddo

    cdif(1) = 1./(1. + adif(1) - adif(1)*e1(2))

    do n=1,ndiff

      xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX) = &
        xn_adv(:,(KMAX_MID-1)*MAXLIMAX*MAXLJMAX)*cdif(KMAX_MID)
      do k = KMAX_MID-2,0,-1
        xn_adv(:,k*MAXLIMAX*MAXLJMAX) =          &
          (xn_adv(:,k*MAXLIMAX*MAXLJMAX)         &
          +adif(k+1)*xn_adv(:,(k+1)*MAXLIMAX*MAXLJMAX))*cdif(k+1)
      enddo

      do k = 1,KMAX_MID-1
        xn_adv(:,k*MAXLIMAX*MAXLJMAX) =                 &
           e1(k+1)*xn_adv(:,(k-1)*MAXLIMAX*MAXLJMAX)    &
          +xn_adv(:,k*MAXLIMAX*MAXLJMAX)
      enddo

    enddo ! ndiff

  end subroutine vertdiffn2

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine advx(vel,velbeg,velend        &
      ,xn_adv,xnbeg,xnend                  &
      ,ps3d,psbeg,psend                    &
      ,xm2loc,xmdloc                       &
      ,dth,fac1)

!     executes advection with a.bott's integrated flux method using
!     4'th order polynomials in the y-direction.
!
!     modified by pw february 2002:  Takes into account the mapping factor
!     in such a way that a Courant number of one corresponds exactly to "empty" a cell.
!     (small effects on results: less than 1%)

    use Par_ml   , only : me,li0,li1,limax
    use ChemSpecs,         only : NSPEC_ADV
    use MassBudget_ml , only : fluxin,fluxout
    implicit none

!    parameter:
!    input
    real,intent(in) :: vel(0:MAXLIMAX),velbeg, velend
    real,intent(in),dimension(NSPEC_ADV,3) :: xnbeg,xnend 
    real,intent(in),dimension(3)           :: psbeg,psend 
    real,intent(in),dimension(0:MAXLIMAX+1):: xm2loc,xmdloc
    real,intent(in) :: dth,fac1

!    input+output
    real ,intent(inout)::xn_adv(NSPEC_ADV,MAXLIMAX)
    real ,intent(inout)::ps3d(MAXLIMAX)

!      output fluxin,fluxout

!    local

    integer ij, ijn,ijll
    integer limtlow,limthig
    integer lijb,lije
    real ijn1
    real x1, x2, hh3,hh4
    real y0,y1,y2,y3
    real zzfc(5,-1:MAXLIMAX+1)
    real fc(-1:MAXLIMAX+1)
    real flux(NSPEC_ADV,-1:MAXLIMAX+1)
    real fluxps(-1:MAXLIMAX+1)
    real hel1(NSPEC_ADV),hel2(NSPEC_ADV)
    real hel1ps,hel2ps
    integer ijpasses
    integer ijb1(MAXLIMAX),ije1(MAXLIMAX)
    integer ijb2(MAXLIMAX),ije2(MAXLIMAX),ijb3(MAXLIMAX)
    logical ijdoend

!-----------------------------------------------------------------------

    limtlow = li0-1
    if (li0.eq.1) then
      if (vel(0) .gt. 0..and.velbeg.lt.0.) then
        fc(-1) = velbeg*dth
        fc(-1) = min( 1.0, fc(-1))
        fc(-1) = max(-1.0, fc(-1))
        limtlow = -1

        y0 = fc(-1)
        x1 = 1.+2.*y0*xm2loc(0)
        x2 = x1*x1
        y3 = xmdloc(0)*(1.-x2)/3840.
        y1 = 5.*y3
        y2 = x1*y3
        hh3 = (116.-4.*x2)*y2
        hh4 = (2.*x2-66.)*y1
        zzfc(3,-1) = - y0 - (214. - 6.*x2)*y2
        zzfc(5,-1) = (y2-y1)*(x2-9.)
        zzfc(1,-1) = (y2+y1)*(x2-9.)
        zzfc(4,-1) = hh3+hh4
        zzfc(2,-1) = hh3-hh4

      endif
    endif

    do 10 ij = li0-1,li1
      fc(ij) = vel(ij)*dth
      fc(ij) = min( 1.0, fc(ij))
      fc(ij) = max(-1.0, fc(ij))

      ijn1 = sign(1.,fc(ij))
      ijn = ij + nint(0.5*(1-ijn1))

      y0 = ijn1*fc(ij)
      x1 = 1.-2.*y0*xm2loc(ijn)
      x2 = x1*x1
      y3 = xmdloc(ijn)*(1.-x2)/3840.
      y1 = 5.*ijn1*y3
      y2 = x1*y3
      hh3 = (116.-4.*x2)*y2
      hh4 = (66.-2.*x2)*y1
      zzfc(3,ij) = y0 - (214. - 6.*x2)*y2
      zzfc(5,ij) = (y2+y1)*(x2-9.)
      zzfc(1,ij) = (y2-y1)*(x2-9.)
      zzfc(4,ij) = hh3+hh4
      zzfc(2,ij) = hh3-hh4

10  continue

    limthig = li1
    if (li1.eq.limax) then
      if (vel(li1).lt.0..and.velend.gt.0.)then
        fc(li1+1) = velend*dth
        fc(li1+1) = min( 1.0, fc(li1+1))
        fc(li1+1) = max(-1.0, fc(li1+1))

        limthig = li1+1

        y0 = fc(li1+1)
        x1 = 1.-2.*y0*xm2loc(li1+1)
        x2 = x1*x1
        y3 = xmdloc(li1+1)*(1.-x2)/3840.
        y1 = 5.*y3
        y2 = x1*y3
        hh3 = (116.-4.*x2)*y2
        hh4 = (66.-2.*x2)*y1
        zzfc(3,li1+1) = y0 - (214.-6.*x2)*y2
        zzfc(5,li1+1) = (y2+y1)*(x2-9.)
        zzfc(1,li1+1) = (y2-y1)*(x2-9.)
        zzfc(4,li1+1) = hh3+hh4
        zzfc(2,li1+1) = hh3-hh4

      endif
    endif

!------- boundary treatment -----------------------------------------

!        helping values at the boundaries are found by linear
!        extrapolation in cases of outflow, and by assuming constant
!        values in inflow cases.

!        calculate the coefficients in the polynomial, the
!        normalized fluxes, and limit them for positivness

    if(limtlow.eq.-1)then

!     integrated flux form

      flux(:,-1) = max(0.,xn_adv(:,2)*zzfc(5,-1)                &
                        + xn_adv(:,1)*zzfc(4,-1)                &
                        + xnbeg(:,3) *zzfc(3,-1)                &
                        + xnbeg(:,2) *zzfc(2,-1)                &
                        + xnbeg(:,1) *zzfc(1,-1))
      flux(:,0)  = max(0.,xn_adv(:,2)*zzfc(5,0)                 &
                        + xn_adv(:,1)*zzfc(4,0)                 &
                        + xnbeg(:,3) *zzfc(3,0)                 &
                        + xnbeg(:,2) *zzfc(2,0)                 &
                        + xnbeg(:,1)*zzfc(1,0))
      fluxps(-1) = max(0.,ps3d(2) *zzfc(5,-1)                   &
                        + ps3d(1) *zzfc(4,-1)                   &
                        + psbeg(3)*zzfc(3,-1)                   &
                        + psbeg(2)*zzfc(2,-1)                   &
                        + psbeg(1)*zzfc(1,-1))
      fluxps(0)  = max(0.,ps3d(2) *zzfc(5,0)                    &
                        + ps3d(1) *zzfc(4,0)                    &
                        + psbeg(3)*zzfc(3,0)                    &
                        + psbeg(2)*zzfc(2,0)                    &
                        + psbeg(1)*zzfc(1,0))

    else

!     integrated flux form

      if(fc(li0-1).ge.0.)then
        flux(:,li0-1) = max(0.,xn_adv(:,li0+1)*zzfc(5,li0-1)    &
                             + xn_adv(:,li0)  *zzfc(4,li0-1)    &
                             + xnbeg(:,3)     *zzfc(3,li0-1)    &
                             + xnbeg(:,2)     *zzfc(2,li0-1)    &
                             + xnbeg(:,1)     *zzfc(1,li0-1))
        fluxps(li0-1) = max(0.,ps3d(li0+1)*zzfc(5,li0-1)        &
                             + ps3d(li0)  *zzfc(4,li0-1)        &
                             + psbeg(3)   *zzfc(3,li0-1)        &
                             + psbeg(2)   *zzfc(2,li0-1)        &
                             + psbeg(1)   *zzfc(1,li0-1))
      else
        flux(:,li0-1) = max(0.,xn_adv(:,li0+2)*zzfc(5,li0-1)    &
                             + xn_adv(:,li0+1)*zzfc(4,li0-1)    &
                             + xn_adv(:,li0)  *zzfc(3,li0-1)    &
                             + xnbeg(:,3)     *zzfc(2,li0-1)    &
                             + xnbeg(:,2)     *zzfc(1,li0-1))
        fluxps(li0-1) = max(0.,ps3d(li0+2)*zzfc(5,li0-1)        &
                             + ps3d(li0+1)*zzfc(4,li0-1)        &
                             + ps3d(li0)  *zzfc(3,li0-1)        &
                             + psbeg(3)   *zzfc(2,li0-1)        &
                             + psbeg(2)   *zzfc(1,li0-1))
      endif
    endif

!     integrated flux form

    if(fc(li0).ge.0.)then
      flux(:,li0) = max(0.,xn_adv(:,li0+2)*zzfc(5,li0)          &
                         + xn_adv(:,li0+1)*zzfc(4,li0)          &
                         + xn_adv(:,li0)  *zzfc(3,li0)          &
                         + xnbeg(:,3)     *zzfc(2,li0)          &
                         + xnbeg(:,2)     *zzfc(1,li0))
      fluxps(li0) = max(0.,ps3d(li0+2)*zzfc(5,li0)              &
                         + ps3d(li0+1)*zzfc(4,li0)              &
                         + ps3d(li0)  *zzfc(3,li0)              &
                         + psbeg(3)   *zzfc(2,li0)              &
                         + psbeg(2)   *zzfc(1,li0))
    else
      flux(:,li0) = max(0.,xn_adv(:,li0+3)*zzfc(5,li0)          &
                         + xn_adv(:,li0+2)*zzfc(4,li0)          &
                         + xn_adv(:,li0+1)*zzfc(3,li0)          &
                         + xn_adv(:,li0)  *zzfc(2,li0)          &
                         + xnbeg(:,3)     *zzfc(1,li0))
      fluxps(li0) = max(0.,ps3d(li0+3)*zzfc(5,li0)              &
                         + ps3d(li0+2)*zzfc(4,li0)              &
                         + ps3d(li0+1)*zzfc(3,li0)              &
                         + ps3d(li0)  *zzfc(2,li0)              &
                         + psbeg(3)   *zzfc(1,li0))
    endif

    if(fc(li0+1).ge.0.)then

!     integrated flux form

      flux(:,li0+1) = max(0.,xn_adv(:,li0+3)*zzfc(5,li0+1)      &
                           + xn_adv(:,li0+2)*zzfc(4,li0+1)      &
                           + xn_adv(:,li0+1)*zzfc(3,li0+1)      &
                           + xn_adv(:,li0)  *zzfc(2,li0+1)      &
                           + xnbeg(:,3)     *zzfc(1,li0+1))
      fluxps(li0+1) = max(0.,ps3d(li0+3)*zzfc(5,li0+1)          &
                           + ps3d(li0+2)*zzfc(4,li0+1)          &
                           + ps3d(li0+1)*zzfc(3,li0+1)          &
                           + ps3d(li0)  *zzfc(2,li0+1)          &
                           + psbeg(3)   *zzfc(1,li0+1))
    endif

    lijb = li0+2
    if(fc(li0+1).lt.0.)lijb = li0+1
    lije = li1-3
    if(fc(li1-2).ge.0.)lije = li1-2

    do ij = lijb,lije

      ijn1 = sign(1.,fc(ij))

!     integrated flux form

      ijn = ij+nint(0.5*(1.-ijn1))
      flux(:,ij) = max(0.,xn_adv(:,ijn+2)*zzfc(5,ij)            &
                        + xn_adv(:,ijn+1)*zzfc(4,ij)            &
                        + xn_adv(:,ijn)  *zzfc(3,ij)            &
                        + xn_adv(:,ijn-1)*zzfc(2,ij)            &
                        + xn_adv(:,ijn-2)*zzfc(1,ij))
      fluxps(ij) = max(0.,ps3d(ijn+2)*zzfc(5,ij)                &
                        + ps3d(ijn+1)*zzfc(4,ij)                &
                        + ps3d(ijn)  *zzfc(3,ij)                &
                        + ps3d(ijn-1)*zzfc(2,ij)                &
                        + ps3d(ijn-2)*zzfc(1,ij))

    enddo

    if(fc(li1-2).lt.0)then

!     integrated flux form

      flux(:,li1-2) = max(0.,xnend(:,1)     *zzfc(5,li1-2)      &
                           + xn_adv(:,li1)  *zzfc(4,li1-2)      &
                           + xn_adv(:,li1-1)*zzfc(3,li1-2)      &
                           + xn_adv(:,li1-2)*zzfc(2,li1-2)      &
                           + xn_adv(:,li1-3)*zzfc(1,li1-2))
      fluxps(li1-2) = max(0.,psend(1)   *zzfc(5,li1-2)          &
                           + ps3d(li1)  *zzfc(4,li1-2)          &
                           + ps3d(li1-1)*zzfc(3,li1-2)          &
                           + ps3d(li1-2)*zzfc(2,li1-2)          &
                           + ps3d(li1-3)*zzfc(1,li1-2))
    endif

!     integrated flux form

    if(fc(li1-1).ge.0)then
      flux(:,li1-1) = max(0.,xnend(:,1)     *zzfc(5,li1-1)      &
                           + xn_adv(:,li1)  *zzfc(4,li1-1)      &
                           + xn_adv(:,li1-1)*zzfc(3,li1-1)      &
                           + xn_adv(:,li1-2)*zzfc(2,li1-1)      &
                           + xn_adv(:,li1-3)*zzfc(1,li1-1))
      fluxps(li1-1) = max(0.,psend(1)   *zzfc(5,li1-1)          &
                           + ps3d(li1)  *zzfc(4,li1-1)          &
                           + ps3d(li1-1)*zzfc(3,li1-1)          &
                           + ps3d(li1-2)*zzfc(2,li1-1)          &
                           + ps3d(li1-3)*zzfc(1,li1-1))
    else
      flux(:,li1-1) = max(0.,xnend(:,2)     *zzfc(5,li1-1)     &
                           + xnend(:,1)     *zzfc(4,li1-1)     &
                           + xn_adv(:,li1)  *zzfc(3,li1-1)     &
                           + xn_adv(:,li1-1)*zzfc(2,li1-1)     &
                           + xn_adv(:,li1-2)*zzfc(1,li1-1))
      fluxps(li1-1) = max(0.,psend(2)   *zzfc(5,li1-1)         &
                           + psend(1)   *zzfc(4,li1-1)         &
                           + ps3d(li1)  *zzfc(3,li1-1)         &
                           + ps3d(li1-1)*zzfc(2,li1-1)         &
                           + ps3d(li1-2)*zzfc(1,li1-1))
    endif

!     integrated flux form

    if(limthig.eq.li1)then

      if(fc(li1).ge.0)then
        flux(:,li1) = max(0.,xnend(:,2)     *zzfc(5,li1)        &
                           + xnend(:,1)     *zzfc(4,li1)        &
                           + xn_adv(:,li1)  *zzfc(3,li1)        &
                           + xn_adv(:,li1-1)*zzfc(2,li1)        &
                           + xn_adv(:,li1-2)*zzfc(1,li1))
        fluxps(li1) = max(0.,psend(2)   *zzfc(5,li1)            &
                           + psend(1)   *zzfc(4,li1)            &
                           + ps3d(li1)  *zzfc(3,li1)            &
                           + ps3d(li1-1)*zzfc(2,li1)            &
                           + ps3d(li1-2)*zzfc(1,li1))
      else
        flux(:,li1) = max(0.,xnend(:,3)     *zzfc(5,li1)        &
                           + xnend(:,2)     *zzfc(4,li1)        &
                           + xnend(:,1)     *zzfc(3,li1)        &
                           + xn_adv(:,li1)  *zzfc(2,li1)        &
                           + xn_adv(:,li1-1)*zzfc(1,li1))
        fluxps(li1) = max(0.,psend(3)   *zzfc(5,li1)            &
                           + psend(2)   *zzfc(4,li1)            &
                           + psend(1)   *zzfc(3,li1)            &
                           + ps3d(li1)  *zzfc(2,li1)            &
                           + ps3d(li1-1)*zzfc(1,li1))
      endif

    else

!     integrated flux form

      flux(:,li1) = max(0.,xnend(:,3)     *zzfc(5,li1)          &
                         + xnend(:,2)     *zzfc(4,li1)          &
                         + xnend(:,1)     *zzfc(3,li1)          &
                         + xn_adv(:,li1)  *zzfc(2,li1)          &
                         + xn_adv(:,li1-1)*zzfc(1,li1))
      flux(:,li1+1) = max(0.,xnend(:,3)     *zzfc(5,li1+1)      &
                           + xnend(:,2)     *zzfc(4,li1+1)      &
                           + xnend(:,1)     *zzfc(3,li1+1)      &
                           + xn_adv(:,li1)  *zzfc(2,li1+1)      &
                           + xn_adv(:,li1-1)*zzfc(1,li1+1))
      fluxps(li1) = max(0.,psend(3)   *zzfc(5,li1)              &
                         + psend(2)   *zzfc(4,li1)              &
                         + psend(1)   *zzfc(3,li1)              &
                         + ps3d(li1)  *zzfc(2,li1)              &
                         + ps3d(li1-1)*zzfc(1,li1))
      fluxps(li1+1) = max(0.,psend(3) *zzfc(5,li1+1)            &
                         + psend(2)   *zzfc(4,li1+1)            &
                         + psend(1)   *zzfc(3,li1+1)            &
                         + ps3d(li1)  *zzfc(2,li1+1)            &
                         + ps3d(li1-1)*zzfc(1,li1+1))

    endif


    if(limtlow.eq.-1)then
      hel1(:) = xnbeg(:,3)*xmdloc(0)
      hel2(:) = flux(:,0) +  flux(:,-1)
      where(hel1(:).lt.hel2(:)) flux(:,0)=flux(:,0)*hel1(:)/(hel2(:)+1.0E-100)
      hel1ps = psbeg(3)*xmdloc(0)
      hel2ps = fluxps(0) +  fluxps(-1)
      if(hel1ps.lt.hel2ps) fluxps(0) = fluxps(0)*hel1ps/hel2ps
      ij = 1
    else
      if(fc(li0-1).ge.0.) then
        flux(:,li0-1) = min(xnbeg(:,3)*xmdloc(li0-1),flux(:,li0-1))
        fluxps(li0-1) = min(psbeg(3)*xmdloc(li0-1),fluxps(li0-1))
        ij = li0
      else
        if(fc(li0).lt.0.) then
          flux(:,li0-1) =-min(xn_adv(:,li0)*xmdloc(li0),flux(:,li0-1))
          fluxps(li0-1) =-min(ps3d(li0)*xmdloc(li0),fluxps(li0-1))
          ij = li0
        else
          hel1(:) = xn_adv(:,li0)*xmdloc(li0)
          hel2(:) = flux(:,li0) + flux(:,li0-1)
          where(hel1(:).lt.hel2(:))
            flux(:,li0-1) =-flux(:,li0-1)*hel1(:)/(hel2(:)+1.0E-100)
            flux(:,li0)   = flux(:,li0  )*hel1(:)/(hel2(:)+1.0E-100)
            xn_adv(:,li0) = 0.
          elsewhere
            flux(:,li0-1) =-flux(:,li0-1)
            xn_adv(:,li0) = xm2loc(li0)*(hel1(:)-hel2(:))
          end where
          hel1ps = ps3d(li0)*xmdloc(li0)
          hel2ps = fluxps(li0) +  fluxps(li0-1)
          if(hel1ps.lt.hel2ps)then
            fluxps(li0-1) =-fluxps(li0-1)*hel1ps/hel2ps
            fluxps(li0)   = fluxps(li0)  *hel1ps/hel2ps
            ps3d(li0) = 0.
          else
            fluxps(li0-1) =-fluxps(li0-1)
            ps3d(li0) =xm2loc(li0)*(hel1ps-hel2ps)
          endif
          ij = li0+1
        endif
      endif
    endif

    ijpasses = 0
    do while(.true.)
      ijpasses = ijpasses+1
      ijb1(ijpasses) = ij
      ije1(ijpasses) = -5
      do while(fc(ij).ge.0.)
        ije1(ijpasses) = ij
        ij = ij+1
        if(ij.gt.li1-1)then
          ijb2(ijpasses) = ij
          ije2(ijpasses) = -5
          ijb3(ijpasses) = -5
          goto 257
        endif
      enddo
      ijb2(ijpasses) = ij
      ije2(ijpasses) = -5
      do while(fc(ij+1).lt.0.)
        ije2(ijpasses) = ij
        ij = ij+1
        if(ij.gt.li1-1)then
          ijb3(ijpasses) = -5
          goto 257
        endif
      enddo
      ijb3(ijpasses) = ij
      ij = ij+2
      if(ij.gt.li1-1)goto 257
    enddo

257 continue
    ijdoend = .false.
    if(ij.eq.li1)ijdoend=.true.

    do ijll = 1,ijpasses
      do ij = ijb1(ijll),ije1(ijll)
          flux(:,ij)   = min(xn_adv(:,ij)*xmdloc(ij),flux(:,ij))
          xn_adv(:,ij) = max(0.,xn_adv(:,ij)                            &
                               -xm2loc(ij)*(flux(:,ij)-flux(:,ij-1)))
          fluxps(ij)   = min(ps3d(ij)*xmdloc(ij),fluxps(ij))
          ps3d(ij)     = max(0.,ps3d(ij)                                &
                               -xm2loc(ij)*(fluxps(ij)-fluxps(ij-1)))
      enddo
      do ij = ijb2(ijll),ije2(ijll)
          flux(:,ij)   =-min(xn_adv(:,ij+1)*xmdloc(ij+1),flux(:,ij))
          xn_adv(:,ij) = max(0.,xn_adv(:,ij)                            &
                               -xm2loc(ij)*(flux(:,ij)-flux(:,ij-1)))
          fluxps(ij)   =-min(ps3d(ij+1)*xmdloc(ij+1),fluxps(ij))
          ps3d(ij)     = max(0.,ps3d(ij)                                &
                               -xm2loc(ij)*(fluxps(ij)-fluxps(ij-1)))
      enddo
      ij = ijb3(ijll)
      if(ij.lt.-3) goto 357
        hel1(:) = xn_adv(:,ij+1)*xmdloc(ij+1)
        hel2(:) = flux(:,ij+1) +  flux(:,ij)

        where(hel1(:).lt.hel2(:))
!On IBM machine the division can give overflow if hel2 is too small
          flux(:,ij)   =-(flux(:,ij)  *hel1(:))/(hel2(:)+1.0E-100)
          flux(:,ij+1) = (flux(:,ij+1)*hel1(:))/(hel2(:)+1.0E-100)
          xn_adv(:,ij+1) = 0.
        elsewhere
          flux(:,ij) = -flux(:,ij)
          xn_adv(:,ij+1) = xm2loc(ij+1)*(hel1(:)-hel2(:))
        end where


        xn_adv(:,ij) = max(0.,xn_adv(:,ij)                              &
                             -xm2loc(ij)*(flux(:,ij)-flux(:,ij-1)))
        hel1ps = ps3d(ij+1)*xmdloc(ij+1)
        hel2ps = fluxps(ij+1) +  fluxps(ij)
        if(hel1ps.lt.hel2ps)then
          fluxps(ij)   =-fluxps(ij)  *hel1ps/hel2ps
          fluxps(ij+1) = fluxps(ij+1)*hel1ps/hel2ps
          ps3d(ij+1) = 0.
        else
          fluxps(ij) = -fluxps(ij)
          ps3d(ij+1) = xm2loc(ij+1)*(hel1ps-hel2ps)
        endif
        ps3d(ij) = max(0.,ps3d(ij)-xm2loc(ij)*(fluxps(ij)-fluxps(ij-1)))
    enddo

357 continue

    if(ijdoend)then
      if(limthig.eq.li1+1)then
        hel1(:) = xnend(:,1)*xmdloc(li1+1)
        hel2(:) = flux(:,li1+1) + flux(:,li1)
        where(hel1(:).lt.hel2(:))
          flux(:,li1) =-flux(:,li1)*hel1(:)/(hel2(:)+1.0E-100)
        elsewhere
          flux(:,li1) =-flux(:,li1)
        end where
        xn_adv(:,li1) = max(0.,xn_adv(:,li1)                            &
                              -xm2loc(li1)*(flux(:,li1)-flux(:,li1-1)))
        hel1ps = psend(1)*xmdloc(li1+1)
        hel2ps = fluxps(li1+1) + fluxps(li1)
        if(hel1ps.lt.hel2ps)then
          fluxps(li1) = -fluxps(li1)*hel1ps/hel2ps
        else
          fluxps(li1) = -fluxps(li1)
        endif
        ps3d(li1) =max(0.,ps3d(li1)                                     &
                         -xm2loc(li1)*(fluxps(li1)-fluxps(li1-1)))

      else

        if(fc(li1).ge.0.) then
          flux(:,li1)  = min(xn_adv(:,li1)*xmdloc(li1),flux(:,li1))
          xn_adv(:,li1)= max(0.,xn_adv(:,li1)                           &
                               -xm2loc(li1)*(flux(:,li1)-flux(:,li1-1)))
          fluxps(li1)  = min(ps3d(li1)*xmdloc(li1),fluxps(li1))
          ps3d(li1)    = max(0.,ps3d(li1)                               &
                                -xm2loc(li1)*(fluxps(li1)-fluxps(li1-1)))
        else
          flux(:,li1)  =-min(xnend(:,1)*xmdloc(li1+1),flux(:,li1))
          xn_adv(:,li1)= max(0.,xn_adv(:,li1)                           &
                               -xm2loc(li1)*(flux(:,li1)-flux(:,li1-1)))
          fluxps(li1)  =-min(psend(1)*xmdloc(li1+1),fluxps(li1))
          ps3d(li1)    = max(0.,ps3d(li1)                               &
                               -xm2loc(li1)*(fluxps(li1)-fluxps(li1-1)))
        endif
      endif
    endif

!     accumulation of the boundary fluxes

    if (li0.eq.2) then
      if(fc(1).ge.0.)then
        fluxin(:)  = fluxin(:)  + flux(:,1)*fac1
      else
        fluxout(:) = fluxout(:) - flux(:,1)*fac1
      endif
    endif

    if (li1.eq.limax-1) then
      if(fc(li1).ge.0.)then
        fluxout(:) = fluxout(:) + flux(:,li1)*fac1
      else
        fluxin(:)  = fluxin(:)  - flux(:,li1)*fac1
      endif
    endif

  end subroutine advx

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine advy(vel,velbeg,velend            &
      ,xn_adv,xnbeg,xnend                      &
      ,ps3d,psbeg,psend                        &
      ,xm2loc,xmdloc                           &
      ,dth,fac1)

!     executes advection with a.bott's integrated flux method using
!     4'th order polynomials in the y-direction.
!
!     modified by pw february 2002:  Takes into account the mapping factor
!     in such a way that a Courant number of one corresponds exactly to "empty" a cell.
!     (small effects on results: less than 1%)
!

    use Par_ml   , only : lj0,lj1,ljmax
    use ChemSpecs,         only : NSPEC_ADV
    use MassBudget_ml , only : fluxin,fluxout
    implicit none

!    parameter:
!    input
    real,intent(in) :: vel(0:MAXLIMAX*MAXLJMAX),velbeg, velend
    real,intent(in),dimension(NSPEC_ADV,3) :: xnbeg,xnend  
    real,intent(in),dimension(3)           :: psbeg,psend  
    real,intent(in),dimension(0:MAXLJMAX+1):: xm2loc,xmdloc
    real,intent(in):: dth,fac1

!    input+output
    real ,intent(inout)::xn_adv(NSPEC_ADV,MAXLIMAX:MAXLIMAX*MAXLJMAX)
    real ,intent(inout)::ps3d(MAXLIMAX:MAXLIMAX*MAXLJMAX)

!      output fluxin,fluxout

!    local

    integer ij, ijn,ijll
    integer limtlow,limthig
    integer lijb,lije
    real ijn1
    real x1, x2, hh3,hh4
    real y0,y1,y2,y3
    real zzfc(5,-1:MAXLJMAX+1)
    real fc(-1:MAXLJMAX+1)
    real flux(NSPEC_ADV,-1:MAXLJMAX+1)
    real fluxps(-1:MAXLJMAX+1)
    real hel1(NSPEC_ADV),hel2(NSPEC_ADV)
    real hel1ps,hel2ps
    integer ijpasses
    integer ijb1(MAXLJMAX),ije1(MAXLJMAX)
    integer ijb2(MAXLJMAX),ije2(MAXLJMAX),ijb3(MAXLJMAX)
    logical ijdoend

!-----------------------------------------------------------------------

    limtlow = lj0-1
    if (lj0.eq.1) then
      if (vel(0) .gt. 0..and.velbeg.lt.0.) then
        fc(-1) = velbeg*dth
        fc(-1) = min( 1.0, fc(-1))
        fc(-1) = max(-1.0, fc(-1))
        limtlow = -1

        y0 = fc(-1)
        x1 = 1.+2.*y0*xm2loc(0)
        x2 = x1*x1
        y3 = xmdloc(0)*(1.-x2)/3840.
        y1 = 5.*y3
        y2 = x1*y3
        hh3 = (116.-4.*x2)*y2
        hh4 = (2.*x2-66.)*y1
        zzfc(3,-1) = - y0 - (214.-6.*x2)*y2
        zzfc(5,-1) = (y2-y1)*(x2-9.)
        zzfc(1,-1) = (y2+y1)*(x2-9.)
        zzfc(4,-1) = hh3+hh4
        zzfc(2,-1) = hh3-hh4

      endif
    endif

    do 10 ij = lj0-1,lj1
      fc(ij) = vel(ij*MAXLIMAX)*dth
      fc(ij) = min( 1.0, fc(ij))
      fc(ij) = max(-1.0, fc(ij))

      ijn1 = sign(1.,fc(ij))
      ijn = ij + nint(0.5*(1-ijn1))

      y0 = ijn1*fc(ij)
      x1 = 1.-2.*y0*xm2loc(ijn)
      x2 = x1*x1
      y3 = xmdloc(ijn)*(1.-x2)/3840.
      y1 = 5.*ijn1*y3
      y2 = x1*y3
      hh3 = (116.-4.*x2)*y2
      hh4 = (66.-2.*x2)*y1
      zzfc(3,ij) = y0 - (214.-6.*x2)*y2
      zzfc(5,ij) = (y2+y1)*(x2-9.)
      zzfc(1,ij) = (y2-y1)*(x2-9.)
      zzfc(4,ij) = hh3+hh4
      zzfc(2,ij) = hh3-hh4

10    continue

    limthig = lj1
    if (lj1.eq.ljmax) then
      if (vel(lj1*MAXLIMAX).lt.0..and.velend.gt.0.)then
        fc(lj1+1) = velend*dth
        fc(lj1+1) = min( 1.0, fc(lj1+1))
        fc(lj1+1) = max(-1.0, fc(lj1+1))
        limthig = lj1+1

        y0 = fc(lj1+1)
        x1 = 1.-2.*y0*xm2loc(lj1+1)
        x2 = x1*x1
        y3 = xmdloc(lj1+1)*(1.-x2)/3840.
        y1 = 5.*y3
        y2 = x1*y3
        hh3 = (116.-4.*x2)*y2
        hh4 = (66.-2.*x2)*y1
        zzfc(3,lj1+1) = y0 - (214.-6.*x2)*y2
        zzfc(5,lj1+1) = (y2+y1)*(x2-9.)
        zzfc(1,lj1+1) = (y2-y1)*(x2-9.)
        zzfc(4,lj1+1) = hh3+hh4
        zzfc(2,lj1+1) = hh3-hh4

      endif
    endif

!------- boundary treatment -----------------------------------------

!        helping values at the boundaries are found by linear
!        extrapolation in cases of outflow, and by assuming constant
!        values in inflow cases.

!        calculate the coefficients in the polynomial, the
!        normalized fluxes, and limit them for positivness

    if(limtlow.eq.-1) then

!     integrated flux form

      flux(:,-1) = max(0.,xn_adv(:,2*MAXLIMAX)*zzfc(5,-1)    &
                        + xn_adv(:,MAXLIMAX)  *zzfc(4,-1)    &
                        + xnbeg(:,3)          *zzfc(3,-1)    &
                        + xnbeg(:,2)          *zzfc(2,-1)    &
                        + xnbeg(:,1)          *zzfc(1,-1))
      flux(:,0) = max(0.,xn_adv(:,(2)*MAXLIMAX)*zzfc(5,0)    &
                       + xn_adv(:,1*MAXLIMAX)  *zzfc(4,0)    &
                       + xnbeg(:,3)            *zzfc(3,0)    &
                       + xnbeg(:,2)            *zzfc(2,0)    &
                       + xnbeg(:,1)            *zzfc(1,0))
      fluxps(-1) = max(0.,ps3d(2*MAXLIMAX)*zzfc(5,-1)        &
                        + ps3d(MAXLIMAX)  *zzfc(4,-1)        &
                        + psbeg(3)        *zzfc(3,-1)        &
                        + psbeg(2)        *zzfc(2,-1)        &
                        + psbeg(1)*zzfc(1,-1))
      fluxps(0) = max(0.,ps3d((2)*MAXLIMAX)*zzfc(5,0)        &
                       + ps3d(1*MAXLIMAX)  *zzfc(4,0)        &
                       + psbeg(3)          *zzfc(3,0)        &
                       + psbeg(2)          *zzfc(2,0)        &
                       + psbeg(1)*zzfc(1,0))

    else

!     integrated flux form

      if(fc(lj0-1).ge.0.)then
        flux(:,lj0-1) = max(0.,xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(5,lj0-1) &
                             + xn_adv(:, lj0   *MAXLIMAX)*zzfc(4,lj0-1) &
                             + xnbeg(:,3)                *zzfc(3,lj0-1) &
                             + xnbeg(:,2)                *zzfc(2,lj0-1) &
                             + xnbeg(:,1)                *zzfc(1,lj0-1))
        fluxps(lj0-1) = max(0.,ps3d((lj0+1)*MAXLIMAX)*zzfc(5,lj0-1)     &
                             + ps3d( lj0   *MAXLIMAX)*zzfc(4,lj0-1)     &
                             + psbeg(3)              *zzfc(3,lj0-1)     &
                             + psbeg(2)              *zzfc(2,lj0-1)     &
                             + psbeg(1)              *zzfc(1,lj0-1))
      else
        flux(:,lj0-1) = max(0.,xn_adv(:,(lj0+2)*MAXLIMAX)*zzfc(5,lj0-1) &
                             + xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(4,lj0-1) &
                             + xn_adv(:, lj0   *MAXLIMAX)*zzfc(3,lj0-1) &
                             + xnbeg(:,3)                *zzfc(2,lj0-1) &
                             + xnbeg(:,2)                *zzfc(1,lj0-1))
        fluxps(lj0-1) = max(0.,ps3d((lj0+2)*MAXLIMAX)*zzfc(5,lj0-1)     &
                             + ps3d((lj0+1)*MAXLIMAX)*zzfc(4,lj0-1)     &
                             + ps3d( lj0   *MAXLIMAX)*zzfc(3,lj0-1)     &
                             + psbeg(3)              *zzfc(2,lj0-1)     &
                             + psbeg(2)              *zzfc(1,lj0-1))
      endif
    endif

!     integrated flux form

    if(fc(lj0).ge.0.)then
      flux(:,lj0) = max(0.,xn_adv(:,(lj0+2)*MAXLIMAX)*zzfc(5,lj0)       &
                         + xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(4,lj0)       &
                         + xn_adv(:, lj0   *MAXLIMAX)*zzfc(3,lj0)       &
                         + xnbeg(:,3)                *zzfc(2,lj0)       &
                         + xnbeg(:,2)                *zzfc(1,lj0))
      fluxps(lj0) = max(0.,ps3d((lj0+2)*MAXLIMAX)*zzfc(5,lj0)           &
                         + ps3d((lj0+1)*MAXLIMAX)*zzfc(4,lj0)           &
                         + ps3d( lj0   *MAXLIMAX)*zzfc(3,lj0)           &
                         + psbeg(3)              *zzfc(2,lj0)           &
                         + psbeg(2)              *zzfc(1,lj0))
    else
      flux(:,lj0) = max(0.,xn_adv(:,(lj0+3)*MAXLIMAX)*zzfc(5,lj0)       &
                         + xn_adv(:,(lj0+2)*MAXLIMAX)*zzfc(4,lj0)       &
                         + xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(3,lj0)       &
                         + xn_adv(:, lj0   *MAXLIMAX)*zzfc(2,lj0)       &
                         + xnbeg(:,3)                *zzfc(1,lj0))
      fluxps(lj0) = max(0.,ps3d((lj0+3)*MAXLIMAX)*zzfc(5,lj0)           &
                         + ps3d((lj0+2)*MAXLIMAX)*zzfc(4,lj0)           &
                         + ps3d((lj0+1)*MAXLIMAX)*zzfc(3,lj0)           &
                         + ps3d( lj0   *MAXLIMAX)*zzfc(2,lj0)           &
                         + psbeg(3)              *zzfc(1,lj0))
    endif

    if(fc(lj0+1).ge.0.)then

!     integrated flux form

      flux(:,lj0+1) = max(0.,xn_adv(:,(lj0+3)*MAXLIMAX)*zzfc(5,lj0+1)   &
                           + xn_adv(:,(lj0+2)*MAXLIMAX)*zzfc(4,lj0+1)   &
                           + xn_adv(:,(lj0+1)*MAXLIMAX)*zzfc(3,lj0+1)   &
                           + xn_adv(:, lj0   *MAXLIMAX)*zzfc(2,lj0+1)   &
                           + xnbeg(:,3)                *zzfc(1,lj0+1))
      fluxps(lj0+1) = max(0.,ps3d((lj0+3)*MAXLIMAX)*zzfc(5,lj0+1)      &
                           + ps3d((lj0+2)*MAXLIMAX)*zzfc(4,lj0+1)      &
                           + ps3d((lj0+1)*MAXLIMAX)*zzfc(3,lj0+1)      &
                           + ps3d( lj0   *MAXLIMAX)*zzfc(2,lj0+1)      &
                           + psbeg(3)              *zzfc(1,lj0+1))
    endif


    lijb = lj0+2
    if(fc(lj0+1).lt.0.)lijb = lj0+1
    lije = lj1-3
    if(fc(lj1-2).ge.0.)lije = lj1-2

    do ij = lijb,lije

      ijn1 = sign(1.,fc(ij))

!     integrated flux form

      ijn = ij+nint(0.5*(1.-ijn1))
      flux(:,ij) = max(0.,xn_adv(:,(ijn+2)*MAXLIMAX)*zzfc(5,ij)        &
                        + xn_adv(:,(ijn+1)*MAXLIMAX)*zzfc(4,ij)        &
                        + xn_adv(:, ijn   *MAXLIMAX)*zzfc(3,ij)        &
                        + xn_adv(:,(ijn-1)*MAXLIMAX)*zzfc(2,ij)        &
                        + xn_adv(:,(ijn-2)*MAXLIMAX)*zzfc(1,ij))
      fluxps(ij) = max(0.,ps3d((ijn+2)*MAXLIMAX)*zzfc(5,ij)            &
                        + ps3d((ijn+1)*MAXLIMAX)*zzfc(4,ij)            &
                        + ps3d( ijn   *MAXLIMAX)*zzfc(3,ij)            &
                        + ps3d((ijn-1)*MAXLIMAX)*zzfc(2,ij)            &
                        + ps3d((ijn-2)*MAXLIMAX)*zzfc(1,ij))

    enddo

    if(fc(lj1-2).lt.0.)then

!     integrated flux form


      flux(:,lj1-2) = max(0.,xnend(:,1)                *zzfc(5,lj1-2)   &
                           + xn_adv(:, lj1   *MAXLIMAX)*zzfc(4,lj1-2)   &
                           + xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(3,lj1-2)   &
                           + xn_adv(:,(lj1-2)*MAXLIMAX)*zzfc(2,lj1-2)   &
                           + xn_adv(:,(lj1-3)*MAXLIMAX)*zzfc(1,lj1-2))
      fluxps(lj1-2) = max(0.,psend(1)*zzfc(5,lj1-2)                     &
                           + ps3d( lj1   *MAXLIMAX)*zzfc(4,lj1-2)       &
                           + ps3d((lj1-1)*MAXLIMAX)*zzfc(3,lj1-2)       &
                           + ps3d((lj1-2)*MAXLIMAX)*zzfc(2,lj1-2)       &
                           + ps3d((lj1-3)*MAXLIMAX)*zzfc(1,lj1-2))

    endif

!     integrated flux form

    if(fc(lj1-1).ge.0.)then

      flux(:,lj1-1) = max(0.,xnend(:,1)                *zzfc(5,lj1-1)   &
                           + xn_adv(:, lj1   *MAXLIMAX)*zzfc(4,lj1-1)   &
                           + xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(3,lj1-1)   &
                           + xn_adv(:,(lj1-2)*MAXLIMAX)*zzfc(2,lj1-1)   &
                           + xn_adv(:,(lj1-3)*MAXLIMAX)*zzfc(1,lj1-1))
      fluxps(lj1-1) = max(0.,psend(1)    *zzfc(5,lj1-1)                 &
                           + ps3d( lj1   *MAXLIMAX)*zzfc(4,lj1-1)       &
                           + ps3d((lj1-1)*MAXLIMAX)*zzfc(3,lj1-1)       &
                           + ps3d((lj1-2)*MAXLIMAX)*zzfc(2,lj1-1)       &
                           + ps3d((lj1-3)*MAXLIMAX)*zzfc(1,lj1-1))

    else

      flux(:,lj1-1) = max(0.,xnend(:,2)      *zzfc(5,lj1-1)             &
                           + xnend(:,1)      *zzfc(4,lj1-1)             &
                           + xn_adv(:, lj1   *MAXLIMAX)*zzfc(3,lj1-1)   &
                           + xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(2,lj1-1)   &
                           + xn_adv(:,(lj1-2)*MAXLIMAX)*zzfc(1,lj1-1))
      fluxps(lj1-1) = max(0.,psend(2)*zzfc(5,lj1-1)                     &
                           + psend(1)*zzfc(4,lj1-1)                     &
                           + ps3d(lj1*MAXLIMAX)*zzfc(3,lj1-1)           &
                           + ps3d((lj1-1)*MAXLIMAX)*zzfc(2,lj1-1)       &
                           + ps3d((lj1-2)*MAXLIMAX)*zzfc(1,lj1-1))

    endif

!     integrated flux form

    if(limthig.eq.lj1)then
      if(fc(lj1).ge.0.)then

        flux(:,lj1) = max(0.,xnend(:,2)                *zzfc(5,lj1)   &
                           + xnend(:,1)                *zzfc(4,lj1)   &
                           + xn_adv(:, lj1   *MAXLIMAX)*zzfc(3,lj1)   &
                           + xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(2,lj1)   &
                           + xn_adv(:,(lj1-2)*MAXLIMAX)*zzfc(1,lj1))
        fluxps(lj1) = max(0.,psend(2)*zzfc(5,lj1)                     &
                           + psend(1)              *zzfc(4,lj1)       &
                           + ps3d( lj1   *MAXLIMAX)*zzfc(3,lj1)       &
                           + ps3d((lj1-1)*MAXLIMAX)*zzfc(2,lj1)       &
                           + ps3d((lj1-2)*MAXLIMAX)*zzfc(1,lj1))

      else

        flux(:,lj1) = max(0.,xnend(:,3)                *zzfc(5,lj1)   &
                           + xnend(:,2)                *zzfc(4,lj1)   &
                           + xnend(:,1)                *zzfc(3,lj1)   &
                           + xn_adv(:, lj1   *MAXLIMAX)*zzfc(2,lj1)   &
                           + xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(1,lj1))
        fluxps(lj1) = max(0.,psend(3)              *zzfc(5,lj1)       &
                           + psend(2)              *zzfc(4,lj1)       &
                           + psend(1)              *zzfc(3,lj1)       &
                           + ps3d( lj1   *MAXLIMAX)*zzfc(2,lj1)       &
                           + ps3d((lj1-1)*MAXLIMAX)*zzfc(1,lj1))

      endif

    else

!     integrated flux form

      flux(:,lj1) = max(0.,xnend(:,3)                *zzfc(5,lj1)     &
                         + xnend(:,2)                *zzfc(4,lj1)     &
                         + xnend(:,1)                *zzfc(3,lj1)     &
                         + xn_adv(:, lj1   *MAXLIMAX)*zzfc(2,lj1)     &
                         + xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(1,lj1))
      flux(:,lj1+1) = max(0.,xnend(:,3)                *zzfc(5,lj1+1) &
                           + xnend(:,2)                *zzfc(4,lj1+1) &
                           + xnend(:,1)                *zzfc(3,lj1+1) &
                           + xn_adv(:, lj1   *MAXLIMAX)*zzfc(2,lj1+1) &
                           + xn_adv(:,(lj1-1)*MAXLIMAX)*zzfc(1,lj1+1))
      fluxps(lj1) = max(0.,psend(3)              *zzfc(5,lj1)         &
                         + psend(2)              *zzfc(4,lj1)         &
                         + psend(1)              *zzfc(3,lj1)         &
                         + ps3d( lj1   *MAXLIMAX)*zzfc(2,lj1)         &
                         + ps3d((lj1-1)*MAXLIMAX)*zzfc(1,lj1))
      fluxps(lj1+1) = max(0.,psend(3)              *zzfc(5,lj1+1)     &
                           + psend(2)              *zzfc(4,lj1+1)     &
                           + psend(1)              *zzfc(3,lj1+1)     &
                           + ps3d( lj1   *MAXLIMAX)*zzfc(2,lj1+1)     &
                           + ps3d((lj1-1)*MAXLIMAX)*zzfc(1,lj1+1))

    endif

    if(limtlow.eq.-1)then
      hel1(:) = xnbeg(:,3)*xmdloc(0)
      hel2(:) = flux(:,0) +  flux(:,-1)
      where(hel1(:).lt.hel2(:)) flux(:,0)=flux(:,0)*hel1(:)/(hel2(:)+1.0E-100)
      hel1ps = psbeg(3)*xmdloc(0)
      hel2ps = fluxps(0) +  fluxps(-1)
      if(hel1ps.lt.hel2ps) fluxps(0)=fluxps(0)*hel1ps/hel2ps
      ij = 1
    else
      if(fc(lj0-1).ge.0.) then
        flux(:,lj0-1) = min(xnbeg(:,3)*xmdloc(lj0-1),flux(:,lj0-1))
        fluxps(lj0-1) = min(psbeg(3)*xmdloc(lj0-1),fluxps(lj0-1))
        ij = lj0
      else
        if(fc(lj0).lt.0.) then
          flux(:,lj0-1)=-min(xn_adv(:,lj0*MAXLIMAX)*xmdloc(lj0),flux(:,lj0-1))
          fluxps(lj0-1)=-min(ps3d(lj0*MAXLIMAX)*xmdloc(lj0),fluxps(lj0-1))
          ij = lj0
        else
          hel1(:) = xn_adv(:,lj0*MAXLIMAX)*xmdloc(lj0)
          hel2(:) = flux(:,lj0) +  flux(:,lj0-1)
          where(hel1(:).lt.hel2(:))
            flux(:,lj0-1) =-flux(:,lj0-1)*hel1(:)/(hel2(:)+1.0E-100)
            flux(:,lj0)   = flux(:,lj0)  *hel1(:)/(hel2(:)+1.0E-100)
            xn_adv(:,lj0*MAXLIMAX) = 0.
          elsewhere
            flux(:,lj0-1) =-flux(:,lj0-1)
            xn_adv(:,lj0*MAXLIMAX) =xm2loc(lj0)*(hel1(:)-hel2(:))
          end where
          hel1ps = ps3d(lj0*MAXLIMAX)*xmdloc(lj0)
          hel2ps = fluxps(lj0) + fluxps(lj0-1)
          if(hel1ps.lt.hel2ps)then
            fluxps(lj0-1) =-fluxps(lj0-1)*hel1ps/hel2ps
            fluxps(lj0)   = fluxps(lj0)  *hel1ps/hel2ps
            ps3d(lj0*MAXLIMAX) = 0.
          else
            fluxps(lj0-1) =-fluxps(lj0-1)
            ps3d(lj0*MAXLIMAX) =xm2loc(lj0)*(hel1ps-hel2ps)
          endif
          ij = lj0+1
        endif
      endif
    endif

    ijpasses = 0
    do while(.true.)

      ijpasses = ijpasses+1
      ijb1(ijpasses) = ij
      ije1(ijpasses) = -5
      do while(fc(ij).ge.0.)
        ije1(ijpasses) = ij
        ij = ij+1
        if(ij.gt.lj1-1)then
          ijb2(ijpasses) = ij
          ije2(ijpasses) = -5
          ijb3(ijpasses) = -5
          goto 257
        endif
      enddo
      ijb2(ijpasses) = ij
      ije2(ijpasses) = -5
      do while(fc(ij+1).lt.0.)
        ije2(ijpasses) = ij
        ij = ij+1
        if(ij.gt.lj1-1)then
          ijb3(ijpasses) = -5
          goto 257
        endif
      enddo
      ijb3(ijpasses) = ij
      ij = ij+2
      if(ij.gt.lj1-1)goto 257
    enddo

257 continue
    ijdoend = .false.
    if(ij.eq.lj1)ijdoend=.true.

    do ijll = 1,ijpasses

      do ij = ijb1(ijll),ije1(ijll)
        flux(:,ij)= min(xn_adv(:,ij*MAXLIMAX)*xmdloc(ij),flux(:,ij))
        xn_adv(:,ij*MAXLIMAX) =                                         &
                    max(0.,xn_adv(:,ij*MAXLIMAX)                        &
                          -xm2loc(ij)*(flux(:,ij)-flux(:,ij-1)))
        fluxps(ij)= min(ps3d(ij*MAXLIMAX)*xmdloc(ij),fluxps(ij))
        ps3d(ij*MAXLIMAX) =                                             &
                    max(0.,ps3d(ij*MAXLIMAX)                            &
                          -xm2loc(ij)*(fluxps(ij)-fluxps(ij-1)))
      enddo
      do ij = ijb2(ijll),ije2(ijll)
        flux(:,ij)=-min(xn_adv(:,(ij+1)*MAXLIMAX)*xmdloc(ij+1),flux(:,ij))
        xn_adv(:,ij*MAXLIMAX) =                                         &
                    max(0.,xn_adv(:,ij*MAXLIMAX)                        &
                          -xm2loc(ij)*(flux(:,ij)-flux(:,ij-1)))
        fluxps(ij)=-min(ps3d((ij+1)*MAXLIMAX)*xmdloc(ij+1),fluxps(ij))
        ps3d(ij*MAXLIMAX) =                                             &
                    max(0.,ps3d(ij*MAXLIMAX)                            &
                          -xm2loc(ij)*(fluxps(ij)-fluxps(ij-1)))
      enddo
      ij = ijb3(ijll)
      if(ij.lt.-3) goto 357
      hel1(:) = xn_adv(:,(ij+1)*MAXLIMAX)*xmdloc(ij+1)
      hel2(:) = flux(:,ij+1) +  flux(:,ij)
      where(hel1(:).lt.hel2(:))
!On IBM machine the division can give overflow if hel2 is too small
        flux(:,ij)   =-flux(:,ij)  *hel1(:)/(hel2(:)+1.0E-100)
        flux(:,ij+1) = flux(:,ij+1)*hel1(:)/(hel2(:)+1.0E-100)
        xn_adv(:,(ij+1)*MAXLIMAX) = 0.
      elsewhere
        flux(:,ij)   =-flux(:,ij)
        xn_adv(:,(ij+1)*MAXLIMAX) = xm2loc(ij+1)*(hel1(:)-hel2(:))
      end where
      xn_adv(:,ij*MAXLIMAX) =                                           &
                    max(0.,xn_adv(:,ij*MAXLIMAX)                        &
                          -xm2loc(ij)*(flux(:,ij)-flux(:,ij-1)))
      hel1ps = ps3d((ij+1)*MAXLIMAX)*xmdloc(ij+1)
      hel2ps = fluxps(ij+1) +  fluxps(ij)
      if(hel1ps.lt.hel2ps)then
        fluxps(ij)   =-fluxps(ij)  *hel1ps/hel2ps
        fluxps(ij+1) = fluxps(ij+1)*hel1ps/hel2ps
        ps3d((ij+1)*MAXLIMAX) = 0.
      else
        fluxps(ij) = -fluxps(ij)
        ps3d((ij+1)*MAXLIMAX) = xm2loc(ij+1)*(hel1ps-hel2ps)
      endif
    ps3d(ij*MAXLIMAX) =                                                 &
                    max(0.,ps3d(ij*MAXLIMAX)                            &
                          -xm2loc(ij)*(fluxps(ij)-fluxps(ij-1)))
    enddo

357 continue

    if(ijdoend)then
      if(limthig.eq.lj1+1)then

        hel1(:) = xnend(:,1)*xmdloc(lj1+1)
        hel2(:) = flux(:,lj1+1) + flux(:,lj1)
        where(hel1(:).lt.hel2(:))
          flux(:,lj1) =-flux(:,lj1)*hel1(:)/(hel2(:)+1.0E-100)
        elsewhere
          flux(:,lj1) =-flux(:,lj1)
        end where
        xn_adv(:,lj1*MAXLIMAX)=                                         &
                    max(0.,xn_adv(:,lj1*MAXLIMAX)                       &
                          -xm2loc(lj1)*(flux(:,lj1)-flux(:,lj1-1)))
        hel1ps = psend(1)*xmdloc(lj1+1)
        hel2ps = fluxps(lj1+1) + fluxps(lj1)
        if(hel1ps.lt.hel2ps)then
          fluxps(lj1) =-fluxps(lj1)*hel1ps/hel2ps
        else
          fluxps(lj1) =-fluxps(lj1)
        endif
        ps3d(lj1*MAXLIMAX) =                                            &
                    max(0.,ps3d(lj1*MAXLIMAX)                           &
                          -xm2loc(lj1)*(fluxps(lj1)-fluxps(lj1-1)))

      else

        if(fc(lj1).ge.0.) then
          flux(:,lj1) =                                                 &
                    min(xn_adv(:,lj1*MAXLIMAX)*xmdloc(lj1),flux(:,lj1))
          xn_adv(:,lj1*MAXLIMAX) =                                      &
                    max(0.,xn_adv(:,lj1*MAXLIMAX)                       &
                          -xm2loc(lj1)*(flux(:,lj1)-flux(:,lj1-1)))
          fluxps(lj1) =                                                 &
                    min(ps3d(lj1*MAXLIMAX)*xmdloc(lj1),fluxps(lj1))
          ps3d(lj1*MAXLIMAX) =                                          &
                    max(0.,ps3d(lj1*MAXLIMAX)                           &
                          -xm2loc(lj1)*(fluxps(lj1)-fluxps(lj1-1)))
        else
          flux(:,lj1)=                                                  &
                   -min(xnend(:,1)*xmdloc(lj1+1),flux(:,lj1))
          xn_adv(:,lj1*MAXLIMAX) =                                      &
                    max(0.,xn_adv(:,lj1*MAXLIMAX)                       &
                          -xm2loc(lj1)*(flux(:,lj1)-flux(:,lj1-1)))
          fluxps(lj1)=                                                  &
                   -min(psend(1)*xmdloc(lj1+1),fluxps(lj1))
          ps3d(lj1*MAXLIMAX) =                                          &
                    max(0.,ps3d(lj1*MAXLIMAX)                           &
                          -xm2loc(lj1)*(fluxps(lj1)-fluxps(lj1-1)))
        endif
      endif
    endif

!     accumulation of the boundary fluxes

    if (lj0.eq.2) then
      if(fc(1).ge.0.)then
        fluxin(:)  = fluxin(:)  + flux(:,1)*fac1
      else
        fluxout(:) = fluxout(:) - flux(:,1)*fac1
      endif
    endif

    if (lj1.eq.ljmax-1) then
      if(fc(lj1).ge.0.)then
        fluxout(:) = fluxout(:) + flux(:,lj1)*fac1
      else
        fluxin(:)  = fluxin(:)  - flux(:,lj1)*fac1
      endif
    endif

  end subroutine advy

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine preadvx(msgnr                &
                    ,xn_adv,ps3d,vel      &
                    ,xnbeg, xnend         &
                    ,psbeg, psend)

    use Par_ml , only : lj0,lj1,li1,neighbor,WEST,EAST
    use ChemSpecs,         only : NSPEC_ADV
    implicit none

!    input
    integer,intent(in):: msgnr
    real,intent(in):: xn_adv(NSPEC_ADV,MAXLIMAX:MAXLIMAX*(MAXLJMAX+1))
    real,intent(in):: ps3d(MAXLIMAX:MAXLIMAX*(MAXLJMAX+1))      &
                     ,vel(MAXLIMAX+1:(MAXLIMAX+1)*(MAXLJMAX+1))

!    output
    real,intent(out),dimension(NSPEC_ADV,3,MAXLJMAX) :: xnend,xnbeg
    real,intent(out),dimension(3,MAXLJMAX)           :: psend,psbeg

!    local
    integer  i, info

    integer request_ps_w, request_ps_e, &
            request_xn_w, request_xn_e
    real,dimension(NSPEC_ADV, 3, MAXLJMAX) :: buf_xn_w,buf_xn_e
    real,dimension(3, MAXLJMAX)            :: buf_ps_w,buf_ps_e

!     Initialize arrays holding boundary slices

!     send to WEST neighbor if any

    if (neighbor(WEST).ge.0) then
      do i = lj0,lj1
        buf_xn_w(:,1,i) = xn_adv(:,i*MAXLIMAX)
        buf_xn_w(:,2,i) = xn_adv(:,i*MAXLIMAX+1)
        buf_xn_w(:,3,i) = xn_adv(:,i*MAXLIMAX+2)

        buf_ps_w(1,i) = ps3d(i*MAXLIMAX)
        buf_ps_w(2,i) = ps3d(i*MAXLIMAX+1)
        buf_ps_w(3,i) = ps3d(i*MAXLIMAX+2)
      enddo

      CALL MPI_ISEND( buf_xn_w, 8*3*MAXLJMAX*NSPEC_ADV, MPI_BYTE, &
          neighbor(WEST), msgnr    , MPI_COMM_WORLD, request_xn_w, INFO)
      CALL MPI_ISEND( buf_ps_w, 8*3*MAXLJMAX          ,MPI_BYTE, &
          neighbor(WEST), msgnr+100, MPI_COMM_WORLD, request_ps_w, INFO)
    endif

    if (neighbor(EAST).ge.0) then
      do i = lj0,lj1
        buf_xn_e(:,1,i) = xn_adv(:,i*MAXLIMAX+li1-3)
        buf_xn_e(:,2,i) = xn_adv(:,i*MAXLIMAX+li1-2)
        buf_xn_e(:,3,i) = xn_adv(:,i*MAXLIMAX+li1-1)

        buf_ps_e(1,i) = ps3d(i*MAXLIMAX+li1-3)
        buf_ps_e(2,i) = ps3d(i*MAXLIMAX+li1-2)
        buf_ps_e(3,i) = ps3d(i*MAXLIMAX+li1-1)
      enddo

      CALL MPI_ISEND( buf_xn_e, 8*3*MAXLJMAX*NSPEC_ADV, MPI_BYTE, &
          neighbor(EAST), msgnr+200, MPI_COMM_WORLD, request_xn_e, INFO)
      CALL MPI_ISEND( buf_ps_e, 8*3*MAXLJMAX          , MPI_BYTE, &
          neighbor(EAST), msgnr+300, MPI_COMM_WORLD, request_ps_e, INFO)
    endif


    if (neighbor(WEST).lt.0) then
      do i = lj0,lj1
        if(vel(i*(MAXLIMAX+1)+1).lt.0)then
          xnbeg(:,2,i) = 3.*xn_adv(:,i*MAXLIMAX+1)    &
              -2.*xn_adv(:,i*MAXLIMAX+2)
          xnbeg(:,3,i) = 2.*xn_adv(:,i*MAXLIMAX+1)    &
              -xn_adv(:,i*MAXLIMAX+2)

          psbeg(2,i) = 3.*ps3d(i*MAXLIMAX+1)-2.*ps3d(i*MAXLIMAX+2)
          psbeg(3,i) = 2.*ps3d(i*MAXLIMAX+1)-ps3d(i*MAXLIMAX+2)
        else
          xnbeg(:,1,i) = xn_adv(:,i*MAXLIMAX)
          xnbeg(:,2,i) = xn_adv(:,i*MAXLIMAX)
          xnbeg(:,3,i) = xn_adv(:,i*MAXLIMAX)

          psbeg(1,i) = ps3d(i*MAXLIMAX)
          psbeg(2,i) = ps3d(i*MAXLIMAX)
          psbeg(3,i) = ps3d(i*MAXLIMAX)
        endif
      enddo
    else

      CALL MPI_RECV( xnbeg, 8*MAXLJMAX*3*NSPEC_ADV, MPI_BYTE, &
          neighbor(WEST), msgnr+200, MPI_COMM_WORLD, MPISTATUS, INFO)
      CALL MPI_RECV( psbeg, 8*MAXLJMAX*3          , MPI_BYTE, &
          neighbor(WEST), msgnr+300, MPI_COMM_WORLD, MPISTATUS, INFO)
    endif

    if (neighbor(EAST).lt.0) then
      do i = lj0,lj1
        if(vel(i*(MAXLIMAX+1)+li1).ge.0)then
          xnend(:,1,i) = 2.*xn_adv(:,i*MAXLIMAX+li1-1)    &
              -xn_adv(:,i*MAXLIMAX+li1-2)
          xnend(:,2,i) = 3.*xn_adv(:,i*MAXLIMAX+li1-1)    &
              -2.*xn_adv(:,i*MAXLIMAX+li1-2)

          psend(1,i) = 2.*ps3d(i*MAXLIMAX+li1-1)    &
              -ps3d(i*MAXLIMAX+li1-2)
          psend(2,i) = 3.*ps3d(i*MAXLIMAX+li1-1)    &
              -2.*ps3d(i*MAXLIMAX+li1-2)
        else
          xnend(:,1,i) = xn_adv(:,i*MAXLIMAX+li1)
          xnend(:,2,i) = xn_adv(:,i*MAXLIMAX+li1)
          xnend(:,3,i) = xn_adv(:,i*MAXLIMAX+li1)

          psend(1,i) = ps3d(i*MAXLIMAX+li1)
          psend(2,i) = ps3d(i*MAXLIMAX+li1)
          psend(3,i) = ps3d(i*MAXLIMAX+li1)
        endif
      enddo
    else

      CALL MPI_RECV( xnend, 8*MAXLJMAX*3*NSPEC_ADV, MPI_BYTE, &
          neighbor(EAST), msgnr    , MPI_COMM_WORLD, MPISTATUS, INFO)
      CALL MPI_RECV( psend, 8*MAXLJMAX*3          , MPI_BYTE, &
          neighbor(EAST), msgnr+100, MPI_COMM_WORLD, MPISTATUS, INFO)
    endif

    !  synchronizing sent buffers (must be done for all ISENDs!!!)
    if (neighbor(WEST) .ge. 0) then
      CALL MPI_WAIT(request_xn_w, MPISTATUS, INFO)
      CALL MPI_WAIT(request_ps_w, MPISTATUS, INFO)
    endif
    if (neighbor(EAST) .ge. 0) then
      CALL MPI_WAIT(request_xn_e, MPISTATUS, INFO)
      CALL MPI_WAIT(request_ps_e, MPISTATUS, INFO)
    endif
  end subroutine preadvx

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine preadvx2(msgnr                &
                     ,xn_adv,ps3d,vel      &
                     ,xnbeg, xnend         &
                     ,psbeg, psend)

!send only one row

    use Par_ml , only : lj0,lj1,li1,neighbor,WEST,EAST
    use ChemSpecs,         only : NSPEC_ADV
    implicit none

!    input
    integer,intent(in):: msgnr
    real,intent(in):: xn_adv(NSPEC_ADV,MAXLIMAX:MAXLIMAX*(MAXLJMAX+1))
    real,intent(in):: ps3d(MAXLIMAX:MAXLIMAX*(MAXLJMAX+1))        &
                     ,vel(MAXLIMAX+1:(MAXLIMAX+1)*(MAXLJMAX+1))

!    output
    real,intent(out),dimension(NSPEC_ADV,3) :: xnend,xnbeg 
    real,intent(out),dimension(3)           :: psend,psbeg 

!    local
    integer  i, info

  integer request_ps_w, request_ps_e, &
          request_xn_w, request_xn_e
  real,dimension(NSPEC_ADV,3) :: buf_xn_w,buf_xn_e  
  real,dimension(3)           :: buf_ps_w,buf_ps_e  

!     Initialize arrays holding boundary slices

!     send to WEST neighbor if any

    if (neighbor(WEST).ge.0) then
      do i = 1,1!lj0,lj1
        buf_xn_w(:,1) = xn_adv(:,i*MAXLIMAX)
        buf_xn_w(:,2) = xn_adv(:,i*MAXLIMAX+1)
        buf_xn_w(:,3) = xn_adv(:,i*MAXLIMAX+2)

        buf_ps_w(1) = ps3d(i*MAXLIMAX)
        buf_ps_w(2) = ps3d(i*MAXLIMAX+1)
        buf_ps_w(3) = ps3d(i*MAXLIMAX+2)
      enddo

      CALL MPI_ISEND( buf_xn_w, 8*3*NSPEC_ADV, MPI_BYTE,&
          neighbor(WEST), msgnr    , MPI_COMM_WORLD, request_xn_w, INFO)
      CALL MPI_ISEND( buf_ps_w, 8*3          , MPI_BYTE,&
          neighbor(WEST), msgnr+100, MPI_COMM_WORLD, request_ps_w, INFO)
    endif

    if (neighbor(EAST).ge.0) then
      do i = 1,1!lj0,lj1
        buf_xn_e(:,1) = xn_adv(:,i*MAXLIMAX+li1-3)
        buf_xn_e(:,2) = xn_adv(:,i*MAXLIMAX+li1-2)
        buf_xn_e(:,3) = xn_adv(:,i*MAXLIMAX+li1-1)

        buf_ps_e(1) = ps3d(i*MAXLIMAX+li1-3)
        buf_ps_e(2) = ps3d(i*MAXLIMAX+li1-2)
        buf_ps_e(3) = ps3d(i*MAXLIMAX+li1-1)
      enddo

      CALL MPI_ISEND( buf_xn_e, 8*3*NSPEC_ADV, MPI_BYTE,&
          neighbor(EAST), msgnr+200, MPI_COMM_WORLD, request_xn_e, INFO)
      CALL MPI_ISEND( buf_ps_e, 8*3          , MPI_BYTE,&
          neighbor(EAST), msgnr+300, MPI_COMM_WORLD, request_ps_e, INFO)
    endif


    if (neighbor(WEST).lt.0) then
      do i = 1,1!lj0,lj1
        if(vel(i*(MAXLIMAX+1)+1).lt.0)then
          xnbeg(:,2) = 3.*xn_adv(:,i*MAXLIMAX+1)    &
              -2.*xn_adv(:,i*MAXLIMAX+2)
          xnbeg(:,3) = 2.*xn_adv(:,i*MAXLIMAX+1)    &
              -xn_adv(:,i*MAXLIMAX+2)

          psbeg(2) = 3.*ps3d(i*MAXLIMAX+1)-2.*ps3d(i*MAXLIMAX+2)
          psbeg(3) = 2.*ps3d(i*MAXLIMAX+1)-ps3d(i*MAXLIMAX+2)
        else
          xnbeg(:,1) = xn_adv(:,i*MAXLIMAX)
          xnbeg(:,2) = xn_adv(:,i*MAXLIMAX)
          xnbeg(:,3) = xn_adv(:,i*MAXLIMAX)

          psbeg(1) = ps3d(i*MAXLIMAX)
          psbeg(2) = ps3d(i*MAXLIMAX)
          psbeg(3) = ps3d(i*MAXLIMAX)
        endif
      enddo
    else

      CALL MPI_RECV( xnbeg, 8*3*NSPEC_ADV, MPI_BYTE, &
          neighbor(WEST), msgnr+200, MPI_COMM_WORLD, MPISTATUS, INFO)
      CALL MPI_RECV( psbeg, 8*3          , MPI_BYTE, &
          neighbor(WEST), msgnr+300, MPI_COMM_WORLD, MPISTATUS, INFO)
    endif

    if (neighbor(EAST).lt.0) then
      do i = 1,1!lj0,lj1
        if(vel(i*(MAXLIMAX+1)+li1).ge.0)then
          xnend(:,1) = 2.*xn_adv(:,i*MAXLIMAX+li1-1)    &
              -xn_adv(:,i*MAXLIMAX+li1-2)
          xnend(:,2) = 3.*xn_adv(:,i*MAXLIMAX+li1-1)    &
              -2.*xn_adv(:,i*MAXLIMAX+li1-2)

          psend(1) = 2.*ps3d(i*MAXLIMAX+li1-1)    &
              -ps3d(i*MAXLIMAX+li1-2)
          psend(2) = 3.*ps3d(i*MAXLIMAX+li1-1)    &
              -2.*ps3d(i*MAXLIMAX+li1-2)
        else
          xnend(:,1) = xn_adv(:,i*MAXLIMAX+li1)
          xnend(:,2) = xn_adv(:,i*MAXLIMAX+li1)
          xnend(:,3) = xn_adv(:,i*MAXLIMAX+li1)

          psend(1) = ps3d(i*MAXLIMAX+li1)
          psend(2) = ps3d(i*MAXLIMAX+li1)
          psend(3) = ps3d(i*MAXLIMAX+li1)
        endif
      enddo
    else

      CALL MPI_RECV( xnend, 8*3*NSPEC_ADV, MPI_BYTE, &
          neighbor(EAST), msgnr    , MPI_COMM_WORLD, MPISTATUS, INFO)
      CALL MPI_RECV( psend, 8*3          , MPI_BYTE, &
          neighbor(EAST), msgnr+100, MPI_COMM_WORLD, MPISTATUS, INFO)
    endif

    !  synchronizing sent buffers (must be done for all ISENDs!!!)
    if (neighbor(WEST) .ge. 0) then
      CALL MPI_WAIT(request_xn_w, MPISTATUS, INFO)
      CALL MPI_WAIT(request_ps_w, MPISTATUS, INFO)
    endif
    if (neighbor(EAST) .ge. 0) then
      CALL MPI_WAIT(request_xn_e, MPISTATUS, INFO)
      CALL MPI_WAIT(request_ps_e, MPISTATUS, INFO)
    endif
  end subroutine preadvx2

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine preadvy(msgnr                &
                    ,xn_adv,ps3d,vel      &
                    ,xnbeg, xnend         &
                    ,psbeg, psend)

    use Par_ml , only : li0,li1,lj0,lj1,ljmax,neighbor,NORTH,SOUTH
    use ChemSpecs,         only : NSPEC_ADV
    implicit none

!    input
    integer,intent(in):: msgnr
    real,intent(in):: xn_adv(NSPEC_ADV,MAXLIMAX*MAXLJMAX)
    real,intent(in):: ps3d(MAXLIMAX*MAXLJMAX)                &
                     ,vel(MAXLIMAX*(MAXLJMAX+1))

!    output
    real,intent(out),dimension(NSPEC_ADV,3,MAXLIMAX) :: xnend,xnbeg 
    real,intent(out),dimension(3,MAXLIMAX)           :: psend,psbeg 

!    local
    integer  i, info

  integer request_ps_n, request_ps_s, &
          request_xn_n, request_xn_s
  real,dimension(NSPEC_ADV,3,MAXLIMAX) :: buf_xn_n,buf_xn_s
  real,dimension(3,MAXLIMAX)           :: buf_ps_n,buf_ps_s

!     Initialize arrays holding boundary slices

!     send to SOUTH neighbor if any

    if (neighbor(SOUTH) .ge. 0) then
      do i = li0,li1
        buf_xn_s(:,1,i) = xn_adv(:,i)
        buf_xn_s(:,2,i) = xn_adv(:,i+MAXLIMAX)
        buf_xn_s(:,3,i) = xn_adv(:,i+2*MAXLIMAX)

        buf_ps_s(1,i) = ps3d(i)
        buf_ps_s(2,i) = ps3d(i+MAXLIMAX)
        buf_ps_s(3,i) = ps3d(i+2*MAXLIMAX)
      enddo

      CALL MPI_ISEND( buf_xn_s, 8*3*MAXLIMAX*NSPEC_ADV, MPI_BYTE,&
            neighbor(SOUTH), msgnr    , MPI_COMM_WORLD, request_xn_s, INFO)
      CALL MPI_ISEND( buf_ps_s, 8*3*MAXLIMAX          , MPI_BYTE,&
            neighbor(SOUTH), msgnr+100, MPI_COMM_WORLD, request_ps_s, INFO)
    endif

    if (neighbor(NORTH) .ge. 0) then
      do i = li0,li1
        buf_xn_n(:,1,i) = xn_adv(:,i+(lj1-3)*MAXLIMAX)
        buf_xn_n(:,2,i) = xn_adv(:,i+(lj1-2)*MAXLIMAX)
        buf_xn_n(:,3,i) = xn_adv(:,i+(lj1-1)*MAXLIMAX)

        buf_ps_n(1,i) = ps3d(i+(lj1-3)*MAXLIMAX)
        buf_ps_n(2,i) = ps3d(i+(lj1-2)*MAXLIMAX)
        buf_ps_n(3,i) = ps3d(i+(lj1-1)*MAXLIMAX)
      enddo

      CALL MPI_ISEND( buf_xn_n, 8*3*MAXLIMAX*NSPEC_ADV, MPI_BYTE,&
            neighbor(NORTH), msgnr    , MPI_COMM_WORLD, request_xn_n, INFO)
      CALL MPI_ISEND( buf_ps_n, 8*3*MAXLIMAX          , MPI_BYTE, &
            neighbor(NORTH), msgnr+100, MPI_COMM_WORLD, request_ps_n, INFO)
    endif

!     receive from SOUTH neighbor if any

    if (neighbor(SOUTH).lt.0) then
      do i = li0,li1
        if(vel(i+MAXLIMAX).lt.0.and.lj0==2)then
          xnbeg(:,3,i) = 2.*xn_adv(:,i+MAXLIMAX)        &
                    -xn_adv(:,i+2*MAXLIMAX)
          xnbeg(:,2,i) = 3.*xn_adv(:,i+MAXLIMAX)        &
                    -2.*xn_adv(:,i+2*MAXLIMAX)
          xnbeg(:,1,i) =  xnbeg(:,2,i)

          psbeg(3,i) = 2.*ps3d(i+MAXLIMAX)-ps3d(i+2*MAXLIMAX)
          psbeg(2,i) = 3.*ps3d(i+MAXLIMAX)-2.*ps3d(i+2*MAXLIMAX)
          psbeg(1,i) = psbeg(2,i)
        else

          xnbeg(:,1,i) = xn_adv(:,i)
          xnbeg(:,2,i) = xn_adv(:,i)
          xnbeg(:,3,i) = xn_adv(:,i)

          psbeg(1,i) = ps3d(i)
          psbeg(2,i) = ps3d(i)
          psbeg(3,i) = ps3d(i)

        endif
      enddo
    else

      CALL MPI_RECV( xnbeg, 8*MAXLIMAX*3*NSPEC_ADV, MPI_BYTE,&
            neighbor(SOUTH), msgnr    , MPI_COMM_WORLD, MPISTATUS, INFO)
      CALL MPI_RECV( psbeg, 8*MAXLIMAX*3          , MPI_BYTE,&
            neighbor(SOUTH), msgnr+100, MPI_COMM_WORLD, MPISTATUS, INFO)
    endif

    if (neighbor(NORTH).lt.0) then
      do i = li0,li1
        if(vel(i+lj1*MAXLIMAX).ge.0.and.ljmax/=lj1)then
          xnend(:,1,i) = 2.*xn_adv(:,i+(lj1-1)*MAXLIMAX)        &
                    -xn_adv(:,i+(lj1-2)*MAXLIMAX)
          xnend(:,2,i) = 3.*xn_adv(:,i+(lj1-1)*MAXLIMAX)        &
                    -2.*xn_adv(:,i+(lj1-2)*MAXLIMAX)
          xnend(:,3,i) = xnend(:,2,i)

          psend(1,i) = 2.*ps3d(i+(lj1-1)*MAXLIMAX)        &
                    -ps3d(i+(lj1-2)*MAXLIMAX)
          psend(2,i) = 3.*ps3d(i+(lj1-1)*MAXLIMAX)        &
                    -2.*ps3d(i+(lj1-2)*MAXLIMAX)
          psend(3,i) = psend(2,i)
        else
          xnend(:,1,i) = xn_adv(:,i+(ljmax-1)*MAXLIMAX)
          xnend(:,2,i) = xn_adv(:,i+(ljmax-1)*MAXLIMAX)
          xnend(:,3,i) = xn_adv(:,i+(ljmax-1)*MAXLIMAX)

          psend(1,i) = ps3d(i+(ljmax-1)*MAXLIMAX)
          psend(2,i) = ps3d(i+(ljmax-1)*MAXLIMAX)
          psend(3,i) = ps3d(i+(ljmax-1)*MAXLIMAX)
        endif
      enddo
    else

      CALL MPI_RECV( xnend, 8*MAXLIMAX*3*NSPEC_ADV, MPI_BYTE,&
            neighbor(NORTH), msgnr    , MPI_COMM_WORLD, MPISTATUS, INFO)
      CALL MPI_RECV( psend, 8*MAXLIMAX*3          , MPI_BYTE,&
            neighbor(NORTH), msgnr+100, MPI_COMM_WORLD, MPISTATUS, INFO)
    endif

!  synchronizing sent buffers (must be done for all ISENDs!!!)
    if (neighbor(SOUTH) .ge. 0) then
      CALL MPI_WAIT(request_xn_s, MPISTATUS, INFO)
      CALL MPI_WAIT(request_ps_s, MPISTATUS, INFO)
    endif
    if (neighbor(NORTH) .ge. 0) then
      CALL MPI_WAIT(request_xn_n, MPISTATUS, INFO)
      CALL MPI_WAIT(request_ps_n, MPISTATUS, INFO)
    endif

  end subroutine preadvy

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine preadvy2(msgnr                &
                     ,xn_adv,ps3d,vel      &
                     ,xnbeg, xnend         &
                     ,psbeg, psend,i_send)

    use Par_ml , only : li0,li1,lj0,lj1,ljmax,neighbor,NORTH,SOUTH
    use ChemSpecs,         only : NSPEC_ADV
    implicit none

!    input
    integer,intent(in):: msgnr,i_send
    real,intent(in):: xn_adv(NSPEC_ADV,MAXLIMAX*MAXLJMAX)
    real,intent(in):: ps3d(MAXLIMAX*MAXLJMAX)                &
                     ,vel(MAXLIMAX*(MAXLJMAX+1))

!    output

    real,intent(out),dimension(NSPEC_ADV,3) :: xnend,xnbeg 
    real,intent(out),dimension(3)           :: psend,psbeg 

!    local
    integer  i, info

    integer request_ps_n, request_ps_s, &
            request_xn_n, request_xn_s
    real,dimension(NSPEC_ADV,3) :: buf_xn_n,buf_xn_s 
    real,dimension(3)           :: buf_ps_n,buf_ps_s 

!     Initialize arrays holding boundary slices

!     send to SOUTH neighbor if any

    if (neighbor(SOUTH) .ge. 0) then
      do i = i_send,i_send
        buf_xn_s(:,1) = xn_adv(:,i)
        buf_xn_s(:,2) = xn_adv(:,i+MAXLIMAX)
        buf_xn_s(:,3) = xn_adv(:,i+2*MAXLIMAX)

        buf_ps_s(1) = ps3d(i)
        buf_ps_s(2) = ps3d(i+MAXLIMAX)
        buf_ps_s(3) = ps3d(i+2*MAXLIMAX)
      enddo

      CALL MPI_ISEND( buf_xn_s, 8*3*NSPEC_ADV, MPI_BYTE,&
            neighbor(SOUTH), msgnr    , MPI_COMM_WORLD, request_xn_s, INFO)
      CALL MPI_ISEND( buf_ps_s, 8*3          , MPI_BYTE,&
            neighbor(SOUTH), msgnr+100, MPI_COMM_WORLD, request_ps_s, INFO)
    endif

    if (neighbor(NORTH) .ge. 0) then
      do i = i_send,i_send
        buf_xn_n(:,1) = xn_adv(:,i+(lj1-3)*MAXLIMAX)
        buf_xn_n(:,2) = xn_adv(:,i+(lj1-2)*MAXLIMAX)
        buf_xn_n(:,3) = xn_adv(:,i+(lj1-1)*MAXLIMAX)

        buf_ps_n(1) = ps3d(i+(lj1-3)*MAXLIMAX)
        buf_ps_n(2) = ps3d(i+(lj1-2)*MAXLIMAX)
        buf_ps_n(3) = ps3d(i+(lj1-1)*MAXLIMAX)
      enddo

      CALL MPI_ISEND( buf_xn_n, 8*3*NSPEC_ADV, MPI_BYTE,&
            neighbor(NORTH), msgnr    , MPI_COMM_WORLD, request_xn_n, INFO)
      CALL MPI_ISEND( buf_ps_n, 8*3          , MPI_BYTE,&
            neighbor(NORTH), msgnr+100, MPI_COMM_WORLD, request_ps_n, INFO)
    endif

!     receive from SOUTH neighbor if any

    if (neighbor(SOUTH).lt.0) then
      do i = i_send,i_send
        if(vel(i+MAXLIMAX).lt.0.and.lj0==2)then
          xnbeg(:,2) = 3.*xn_adv(:,i+MAXLIMAX)        &
                    -2.*xn_adv(:,i+2*MAXLIMAX)
          xnbeg(:,3) = 2.*xn_adv(:,i+MAXLIMAX)        &
                    -xn_adv(:,i+2*MAXLIMAX)
          xnbeg(:,1) =  xnbeg(:,2)

          psbeg(2) = 3.*ps3d(i+MAXLIMAX)-2.*ps3d(i+2*MAXLIMAX)
          psbeg(3) = 2.*ps3d(i+MAXLIMAX)-ps3d(i+2*MAXLIMAX)
          psbeg(1) = psbeg(2)
        else
          xnbeg(:,1) = xn_adv(:,i)
          xnbeg(:,2) = xnbeg(:,1)
          xnbeg(:,3) = xnbeg(:,1)

          psbeg(1) = ps3d(i)
          psbeg(2) =  psbeg(1)
          psbeg(3) =  psbeg(1)
        endif
      enddo
    else

      CALL MPI_RECV( xnbeg, 8*3*NSPEC_ADV, MPI_BYTE,&
            neighbor(SOUTH), msgnr    , MPI_COMM_WORLD, MPISTATUS, INFO)
      CALL MPI_RECV( psbeg, 8*3         , MPI_BYTE,&
            neighbor(SOUTH), msgnr+100, MPI_COMM_WORLD, MPISTATUS, INFO)
    endif

    if (neighbor(NORTH).lt.0) then
      do i = i_send,i_send
        if(vel(i+lj1*MAXLIMAX).ge.0.and.ljmax/=lj1)then
          xnend(:,1) = 2.*xn_adv(:,i+(lj1-1)*MAXLIMAX)        &
                    -xn_adv(:,i+(lj1-2)*MAXLIMAX)
          xnend(:,2) = 3.*xn_adv(:,i+(lj1-1)*MAXLIMAX)        &
                    -2.*xn_adv(:,i+(lj1-2)*MAXLIMAX)
          xnend(:,3) = xnend(:,2)

          psend(1) = 2.*ps3d(i+(lj1-1)*MAXLIMAX)        &
                    -ps3d(i+(lj1-2)*MAXLIMAX)
          psend(2) = 3.*ps3d(i+(lj1-1)*MAXLIMAX)        &
                    -2.*ps3d(i+(lj1-2)*MAXLIMAX)
          psend(3) = psend(2)
        else
          xnend(:,1) = xn_adv(:,i+(ljmax-1)*MAXLIMAX)
          xnend(:,2) = xnend(:,1)
          xnend(:,3) = xnend(:,1)

          psend(1) = ps3d(i+(ljmax-1)*MAXLIMAX)
          psend(2) = psend(1)
          psend(3) = psend(1)
        endif
      enddo
    else

      CALL MPI_RECV( xnend, 8*3*NSPEC_ADV, MPI_BYTE,&
            neighbor(NORTH), msgnr    , MPI_COMM_WORLD, MPISTATUS, INFO)
      CALL MPI_RECV( psend, 8*3          , MPI_BYTE,&
            neighbor(NORTH), msgnr+100, MPI_COMM_WORLD, MPISTATUS, INFO)
    endif

!  synchronizing sent buffers (must be done for all ISENDs!!!)
    if (neighbor(SOUTH) .ge. 0) then
      CALL MPI_WAIT(request_xn_s, MPISTATUS, INFO)
      CALL MPI_WAIT(request_ps_s, MPISTATUS, INFO)
    endif
    if (neighbor(NORTH) .ge. 0) then
      CALL MPI_WAIT(request_xn_n, MPISTATUS, INFO)
      CALL MPI_WAIT(request_ps_n, MPISTATUS, INFO)
    endif

  end subroutine preadvy2

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!  subroutine convection_pstar(ps3d,dt_conv)
!  moved to Convection_ml.f90


 
  subroutine alloc_adv_arrays

    !allocate the arrays once
!    allocate(uw(MAXLJMAX,KMAX_MID,NMET),ue(MAXLJMAX,KMAX_MID,NMET))
!    allocate(vs(MAXLIMAX,KMAX_MID,NMET),vn(MAXLIMAX,KMAX_MID,NMET))
    allocate(dhs1(KMAX_BND), dhs1i(KMAX_BND), dhs2i(KMAX_BND))
    allocate(alfnew(9,2:KMAX_MID,0:1))

  end subroutine alloc_adv_arrays

end module Advection_ml
