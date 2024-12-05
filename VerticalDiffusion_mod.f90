! < VerticalDiffusionm_mod.f90 - A component of the EMEP MSC-W Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2024 met.no
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
module VerticalDiffusion_mod
! AIMS to mimic from advec
! call vertdiffn(xn_adv(1,i,j,1),NSPEC_ADV,LIMAX*LJMAX,1,EtaKz(i,j,1,1),ds3,ds4,ndiff)
! but move ds3, ds4 into routine
! note that dhs1(k+1) denotes thickness of layer k
  use ChemDims_mod,       only: NSPEC_ADV
  use Config_module,only: KMAX_MID, KMAX_BND
  use GridValues_mod,     only: dhs1, dhs1i, dhs2i, &
                               debug_proc, debug_li, debug_lj  !DSKZ
  use Par_mod,           only: me,LIMAX,LJMAX
  implicit none
  private

  public :: vertdiffn

  !Older routines, kept for reminder
 public :: vertdiff_1d


contains

  
  !DSKz subroutine vertdiffn(xn_adv,NSPEC,Nij,KMIN_in,SigmaKz,ds3,ds4,ndiff)
  subroutine vertdiffn(xn_adv,NSPEC,Nij,KMIN_in,SigmaKz,dt_diff,ndiff)

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
!
! KMIN is the minimum value of k to include for diffusion (1 is top)

!    input
    real,intent(in)::  SigmaKz(0:LIMAX*LJMAX*KMAX_BND-1)
    real, intent(in) :: dt_diff ! DSKz
!DSKz    real,intent(in)::  ds3(KMAX_MID-1),ds4(KMAX_MID-1)
    integer,intent(in)::  NSPEC,ndiff,Nij
    integer,intent(in):: KMIN_in

!    output
    real ,intent(inout):: xn_adv(NSPEC,0:Nij*(KMAX_MID-1))

!    local

    integer  k,n,KMIN

    real, dimension(0:KMAX_MID-1) :: adif,bdif,cdif,e1
    !DSKz?? real::  ds3(KMAX_MID-1),ds4(KMAX_MID-1) ! DSKz
    real::  ds3(2:KMAX_MID),ds4(2:KMAX_MID) ! DSKz

    real ndiffi

    if(KMIN_in==KMAX_MID)return!no diffusion from one cell to itself
    KMIN = KMIN_in-1
    KMIN = max(1,KMIN)

    ndiffi=1./ndiff

    !DSKz add here
    ! dsh things are vericla map factors
    ! dhs1(KMAX_BND), dhs1i(KMAX_BND), dhs2i(KMAX_BND)
    do k = 2,KMAX_MID
       !ds3(k) = dt_advec*dhs1i(k)*dhs2i(k)
       !ds4(k) = dt_advec*dhs1i(k+1)*dhs2i(k)
       ds3(k) = dt_diff*dhs1i(k)*dhs2i(k)
       ds4(k) = dt_diff*dhs1i(k+1)*dhs2i(k)
    end do

    !END DSKz add here

    do k = KMIN,KMAX_MID-1
      !DSKz adif(k-1) = SigmaKz(k*LIMAX*LJMAX)*ds3(k)*ndiffi
      !DSKz bdif(k) = SigmaKz(k*LIMAX*LJMAX)*ds4(k)*ndiffi
      adif(k-1) = SigmaKz(k*LIMAX*LJMAX)*ds3(k+1)*ndiffi
      bdif(k) = SigmaKz(k*LIMAX*LJMAX)*ds4(k+1)*ndiffi
    end do

   cdif(KMAX_MID-1) = 1./(1. + bdif(KMAX_MID-1))
    e1(KMAX_MID-1) = bdif(KMAX_MID-1)*cdif(KMAX_MID-1)

    do k = KMAX_MID-2,KMIN,-1
      cdif(k) = 1./(1. + bdif(k) + adif(k) - adif(k)*e1(k+1))
      e1(k) = bdif(k)*cdif(k)
    end do

    cdif(KMIN-1) = 1./(1. + adif(KMIN-1) - adif(KMIN-1)*e1(KMIN))

    do n=1,ndiff

      xn_adv(:,Nij*(KMAX_MID-1)) = &
         xn_adv(:,Nij*(KMAX_MID-1))*cdif(KMAX_MID-1)

      do k = KMAX_MID-2,KMIN-1,-1
         xn_adv(:,Nij*k) =         &
           (xn_adv(:,Nij*k)        &
           +adif(k)*xn_adv(:,Nij*(k+1)))*cdif(k)
      end do

      do k = KMIN,KMAX_MID-1
         xn_adv(:,Nij*k) =            &
            e1(k)*xn_adv(:,Nij*(k-1)) &
           +xn_adv(:,Nij*k)
      end do

    end do ! ndiff

  end subroutine vertdiffn

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! OLDER ROUTINES!!! Not used now
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine vertdiff(xn_adv,SigmaKz,ds3,ds4)

!     executes vertical diffusion
!    input
    real,intent(in)::  SigmaKz(0:LIMAX*LJMAX*KMAX_BND-1)
    real,intent(in)::  ds3(KMAX_MID-1),ds4(KMAX_MID-1)

!    output
    real ,intent(inout):: xn_adv(NSPEC_ADV,0:LIMAX*LJMAX*(KMAX_MID-1))

!    local

    integer  k
    real, dimension(KMAX_MID) :: adif,bdif,cdif,e1

    do k = 1,KMAX_MID-1
      adif(k) = SigmaKz(k*LIMAX*LJMAX)*ds3(k)
      bdif(k+1) = SigmaKz(k*LIMAX*LJMAX)*ds4(k)
    end do

    cdif(KMAX_MID) = 1./(1. + bdif(KMAX_MID))
    e1(KMAX_MID) = bdif(KMAX_MID)*cdif(KMAX_MID)
    xn_adv(:,(KMAX_MID-1)*LIMAX*LJMAX) = &
      xn_adv(:,(KMAX_MID-1)*LIMAX*LJMAX)*cdif(KMAX_MID)

    do k = KMAX_MID-1,2,-1
      cdif(k) = 1./(1. + bdif(k) + adif(k) - adif(k)*e1(k+1))
      e1(k) = bdif(k)*cdif(k)
      xn_adv(:,(k-1)*LIMAX*LJMAX) =                 &
          (xn_adv(:,(k-1)*LIMAX*LJMAX)              &
         +adif(k)*xn_adv(:,(k)*LIMAX*LJMAX))*cdif(k)
    end do

    cdif(1) = 1./(1. + adif(1) - adif(1)*e1(2))
    xn_adv(:,0) = (xn_adv(:,0) + adif(1)*xn_adv(:,LIMAX*LJMAX))*cdif(1)

    do k = 2,KMAX_MID
      xn_adv(:,(k-1)*LIMAX*LJMAX) =                &
          e1(k)*xn_adv(:,(k-2)*LIMAX*LJMAX)        &
         +xn_adv(:,(k-1)*LIMAX*LJMAX)
    end do

  end subroutine vertdiff
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine vertdiff_1d(xn_adv,SigmaKz,ds3,ds4,ndiff)

!     executes vertical diffusion

    implicit none

!    input
    real,intent(in)::  SigmaKz(0:LIMAX*LJMAX*KMAX_BND-1)
    real,intent(in)::  ds3(KMAX_MID-1),ds4(KMAX_MID-1)

!    output
    real ,intent(inout):: xn_adv(KMAX_MID)

!    local
    integer, intent(in)::ndiff
    integer  k,n
    real, dimension(KMAX_MID) :: adif,bdif,cdif,e1
    real ndiffi

!    ndiff=1
    ndiffi=1./ndiff

    do k = 1,KMAX_MID-1
      adif(k) = SigmaKz(k*LIMAX*LJMAX)*ds3(k)*ndiffi
      bdif(k+1) = SigmaKz(k*LIMAX*LJMAX)*ds4(k)*ndiffi
    end do

    cdif(KMAX_MID) = 1./(1. + bdif(KMAX_MID))
    e1(KMAX_MID) = bdif(KMAX_MID)*cdif(KMAX_MID)

    do k = KMAX_MID-1,2,-1
      cdif(k) = 1./(1. + bdif(k) + adif(k) - adif(k)*e1(k+1))
      e1(k) = bdif(k)*cdif(k)
    end do

    cdif(1) = 1./(1. + adif(1) - adif(1)*e1(2))

    do n=1,ndiff

    xn_adv(KMAX_MID) = &
      xn_adv(KMAX_MID)*cdif(KMAX_MID)

    do k = KMAX_MID-1,2,-1
      xn_adv(k) =                 &
          (xn_adv(k)              &
         +adif(k)*xn_adv(k+1))*cdif(k)
    end do

    xn_adv(1) = (xn_adv(1) + adif(1)*xn_adv(2))*cdif(1)

    do k = 2,KMAX_MID
      xn_adv(k) =                &
          e1(k)*xn_adv(k-1)        &
         +xn_adv(k)
    end do
    end do

  end subroutine vertdiff_1d
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine vertdiffn2(xn_adv,SigmaKz,ds3,ds4,ndiff)

!     executes vertical diffusion ndiff times

!    input
    real,intent(in)::  SigmaKz(0:LIMAX*LJMAX*KMAX_BND-1)
    real,intent(in)::  ds3(KMAX_MID-1),ds4(KMAX_MID-1)
    integer,intent(in)::  ndiff

!    output
    real ,intent(inout):: xn_adv(NSPEC_ADV,0:LIMAX*LJMAX*(KMAX_MID-1))

!    local

    integer  k,n

    real, dimension(KMAX_MID) :: adif,bdif,cdif,e1

    real ndiffi

    ndiffi=1./ndiff

    do k = 1,KMAX_MID-1
      adif(k) = SigmaKz(k*LIMAX*LJMAX)*ds3(k)*ndiffi
      bdif(k+1) = SigmaKz(k*LIMAX*LJMAX)*ds4(k)*ndiffi
    end do

    cdif(KMAX_MID) = 1./(1. + bdif(KMAX_MID))
    e1(KMAX_MID) = bdif(KMAX_MID)*cdif(KMAX_MID)

    do k = KMAX_MID-1,2,-1
      cdif(k) = 1./(1. + bdif(k) + adif(k) - adif(k)*e1(k+1))
      e1(k) = bdif(k)*cdif(k)
    end do

    cdif(1) = 1./(1. + adif(1) - adif(1)*e1(2))

    do n=1,ndiff

      xn_adv(:,(KMAX_MID-1)*LIMAX*LJMAX) = &
        xn_adv(:,(KMAX_MID-1)*LIMAX*LJMAX)*cdif(KMAX_MID)
      do k = KMAX_MID-2,0,-1
        xn_adv(:,k*LIMAX*LJMAX) =          &
          (xn_adv(:,k*LIMAX*LJMAX)         &
          +adif(k+1)*xn_adv(:,(k+1)*LIMAX*LJMAX))*cdif(k+1)
      end do

      do k = 1,KMAX_MID-1
        xn_adv(:,k*LIMAX*LJMAX) =                 &
           e1(k+1)*xn_adv(:,(k-1)*LIMAX*LJMAX)    &
          +xn_adv(:,k*LIMAX*LJMAX)
      end do

    end do ! ndiff

  end subroutine vertdiffn2

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module VerticalDiffusion_mod
