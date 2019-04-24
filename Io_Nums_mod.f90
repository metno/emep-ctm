! <Io_Nums_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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
module Io_Nums_mod

!_____________________________________________________________________________
! The idea is to keep all unit numbers used for input and output
! stored here. So that programmers can quickly see which numbers are
! in use or not. (Warning - some of the routines mentioned are outdated)
!
! Assign unit number, e.g. io_xxx, here, and use read(io_xxx,*) in
! main program.
!
! (c) - file opened AND closed in subroutine
! (o) - file remains open in subroutine
!_____________________________________________________________________________
implicit none


  integer, parameter, public  :: &
    IO_LOG      = 7   &! General output log (o)
   ,IO_SITES    = 8   &! sites module, first for input(c)
   ,IO_MYTIM    = 20  &! emepctm.f90(c)-output mytim.out 
   ,IO_RES      = 25  &! o3mod,massbud(o) 
   ,IO_DEBUG    = 26  &! keep open, use where needed
   ,IO_TMP      = 27  &! General IO number (files *must* be type (c))
   ,IO_NML      = 28   ! For namelist file (o)


  integer, parameter, public  :: &
    IO_SONDES   = 30  &! siteswrt_mod(o)  for output of sonde data
   ,IO_WRTCHEM  = 118 &! Used in Wrtchem (c) for AOT and BCs
   ,IO_HOURLY   = 119 &! hourly_out(o)
   ,IO_SPOD     = 120  ! Used in DryDep (c) for ozone-flux
      ! ****************** !
      ! CAREFUL. Code uses IO_SPOD + me, so potentally 120..very high IO num
      ! ****************** !


  !(some subroutine names are a bit outdated:)
  integer, parameter, public  :: &
    IO_FORES    = 49  &! rforest.f(c)-read land use %
   ,IO_AIRN     = 49  &! airnox.f(c) - read aircr. em.
   ,IO_LIGHT    = 49  &! lightning.f(c) - read lightning. emiss.
   ,IO_JOST     = 49  &! newjostinit(c) - read  global mixing ratios
   ,IO_GLOBBC   = 49  &! read  global mixing ratios e.g. Logan
   ,IO_GLOBBC2  = 91  &! read  global mixing ratios e.g. h2o2
   ,IO_DO3SE    = 51  &! for DO3SE inputs(c)
   ,IO_ROUGH    = 52  &! inpar.f -reads roughn. class   for landsea masl
   ,IO_VOLC     = 54  &
   ,IO_DJ       = 55  &! readdiss.f(c) - inp. solar r.
   ,IO_AIRCR    = 66  &! phyche.f(c) - write aircraft conc.
   ,IO_OUT      = 80  &! (c)write outday etc.
   ,IO_UKDEP    = 81  &! (o)write fluxes, etc.
   ,IO_STAB     = 82  &! (o)write fluxes, etc.
   ,IO_EMIS     = 84  &! Used for femis , emis_split(c)
   ,IO_TIMEFACS = 85  &! Used for monthly
   ,IO_NEST     = 88  &!   
   ,IO_DMS      = 90  &!  Emissions(c): for DMS 
   ,IO_NH3      = 92  &! hb NH3emis
   ,IO_CLAY     = 92  &! clay
   ,IO_SAND     = 93  &! sand
   ,IO_NH3_DEB  = 99   ! hb NH3emis

end module Io_Nums_mod
