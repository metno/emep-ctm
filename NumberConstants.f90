! <NumberConstants.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module NumberConstants
  implicit none
  private

! from KPP system
! Definition of different levels of accuracy
! for REAL variables using KIND parameterization
!
! KPP SP - Single precision kind
  integer, public, parameter :: sp = selected_real_kind(6,30)
! KPP DP - Double precision kind
  integer, public, parameter :: dp = selected_real_kind(14,300)
! CYGWIN can't handle quad precision, so we re-define
! KPP QP - Quadruple precision kind
  integer, public, parameter :: qp = dp ! selected_real_kind(18,400)

!DEWS working precision can be changed here
! Typically should be dp, but qp for testing

  integer, public, parameter :: wp = qp

!integer, public, parameter :: dp = kind(0.0d0)  ! Double precision real(qp) kind

! Sentinel values
!real(qp), public, parameter :: UNDEF_D = -huge(0.0_dp)
    real,     public, parameter :: UNDEF_R = -huge(0.0)
    real(wp), public, parameter :: UNDEF_WP = -huge(0.0_wp)
    integer,  public, parameter :: UNDEF_I = -huge(0)

end module NumberConstants
