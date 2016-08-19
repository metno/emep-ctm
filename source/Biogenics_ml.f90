! <Biogenics_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Biogenics_ml
  !/-- Reads in BVOC emisions factors (for "standard" conditions, 
  !    30 deg C and sunlit).
  !    The effects of temperature and light on the biogenic stuff is calculated
  !     in Setup_1d
  !---------------------------------------------------------------------------
  use My_Emis_ml       , only : NBVOC, BVOC_USED

  use CheckStop_ml,      only: CheckStop
  use GridValues_ml    , only : xm2, gb, &
          i_fdom,j_fdom,debug_proc,debug_li,debug_lj
  use Io_ml            , only : IO_FORES, open_file, ios, Read2DN
  use KeyValue_ml,       only : KeyVal,KeyValue
  use ModelConstants_ml, only : NPROC
  use Par_ml   , only : me, MAXLIMAX,MAXLJMAX,MSG_READ1,li0,li1,lj0,lj1
  implicit none
  private

  !/-- subroutines
  public ::  Init_BVOC

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
  logical, private, parameter:: DEBUG = .false.
  integer, public, save :: BIO_ISOP, BIO_TERP

  real, public, save, dimension(MAXLIMAX,MAXLJMAX,NBVOC) :: &
      emforest    & !  Gridded standard (30deg. C, full light) emissions
     ,emnat         !  Gridded std. emissions after scaling with density, etc.

  !/-- Canopy environmental correction factors-----------------------------
  !
  !    - to correct for temperature and light of the canopy 
  !    - from Guenther's papers. (Limit to 0<T<40 deg C.)

  real, public, save, dimension(NBVOC,40) :: &
        canopy_ecf  ! Canopy env. factors
                                                        
  !/-- DMS factors

  integer, public, parameter :: IQ_DMS = 35  ! code for DMS emissions

  logical, public, save :: first_dms_read

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    subroutine Init_BVOC()

!    Read natural BVOC emission potentials
!-----------------------------------------------------------------------------
!   Emissions now read from 50x50 landuse file, forests.dat, derived from 
!   landuse.mar2004. Emission rates now based upon Simpson et al., 1999, JGR,
!   Vol 104, D7, 8113-8152.


    real  :: biofac(NBVOC)     !Scaling factor
    integer i, j, n ,d, info, it
    real sumland

    real ::  bvocsum, bvocsum1
    real agts, agr, agtm, agct1, agct2, agct ,itk, fac

    ! Specify the assumed coords and units - Read2DN will check that the data
    ! conform to these.
    type(keyval), dimension(2) :: CheckValues = (/ keyval("Units","ug/m2/h"), &
                                                  keyval("Coords","ModelCoords") /)

   ! Check names, and derive factor used in Emissions_ml to get from ug/m2/s
   ! to molecules/cm2/s  (needs more documentation!)
    do i = 1, NBVOC
      if ( BVOC_USED(i) == "isoprene" ) then
            BIO_ISOP = i
            biofac(i)  = 1.e-9/(68.*3600.)
      else if ( BVOC_USED(i) == "terpene " ) then
            BIO_TERP = i
            biofac(i)  = 1.e-9/(136.*3600.)
      else
            call CheckStop( " BIO ERROR " )
      end if
    end do

    emforest = 0.0

   !========= Read in Standard (30 deg C, full sunlight emissions factors = !

    call Read2DN("Inputs.BVOC",2,emforest,CheckValues)

   !========================================================================!


      if( debug_proc .and. DEBUG) then
          write(*,"(a8,i3,2i4,4f18.4)") "BIONEW ", NBVOC, &
              i_fdom(debug_li), j_fdom(debug_lj), &
             ( emforest(debug_li,debug_lj,i), i=1,NBVOC)
      end if

     !output sums. Remember that "shadow" area not included here.
      do i = 1, NBVOC 
         bvocsum   = sum ( emforest(li0:li1,lj0:lj1,i) )
         CALL MPI_ALLREDUCE(bvocsum,bvocsum1, 1, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 

         if (me == 0 .and. DEBUG) then
            write(6,"(a20,i4,2es12.4)") 'Biogenics_ml, ibio, sum1',&
              me, bvocsum, bvocsum1
         end if
      end do

    ! And now we scale emforest with biofac:

   
     do i = 1, NBVOC
       emforest(:,:,i) = emforest(:,:,i) * biofac(i)
     end do

    !
    !/-- Tabulate canopy environmental correction factors
    !-----------------------------------------------------
    ! ag suffix  = Alex Guenther's parameter values

    agts = 303.
    agr = 8.314
    agtm = 314.
    agct1 = 95000.
    agct2 = 230000.

    do it = 1,40
      itk = it + 273.15
      agct = exp(agct1*(itk - agts)/(agr*agts*itk)) / &
                (1. + exp(agct2*(itk - agtm)/(agr*agts*itk)))

      canopy_ecf(BIO_ISOP,it) = agct

      ! Terpenes
      agct = exp( 0.09*(itk-agts) )
      ! for terpene fac = 0.5*fac(iso): as mass terpene = 2xmass isoprene

      canopy_ecf(BIO_TERP,it) = agct

      if(DEBUG .and. me == 0) &
             write(6,"(A12,i4,5g12.3)") 'Biogenic ecfs: ', &
                  it, ( canopy_ecf(i,it), i=1, NBVOC)
    end do

   end subroutine Init_BVOC
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module Biogenics_ml
