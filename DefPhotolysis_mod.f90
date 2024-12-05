! <DefPhotolysis_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.5>
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
! <DefPhotolysis_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
!+ Photolysis coefficients
!-------------------------------------------------------------------------------

module DefPhotolysis_mod
  !-------------------------------------------------------------------------------
  !
  !   11/12/23 - Willem van Caspel; removed tabulated system. Now only for 
  !                    photolysis rate index declarations.
  !              
  !-------------------------------------------------------------------------------
    
  implicit none
  private

  integer, public, parameter :: NRCPHOTextended  = 172   ! currently max number, used for CRIv2R5

  ! Below photolysis indices are used to map CloudJ Jvalues to the chemistry photolysis reactions.
  ! Set to > 0 if found in CloudJ input file. If not found, and present in chemistry, the model stops.
  ! Matching/mapping of the Jvalues is described in some more detail in the genChem chemical pre-processor system.
  integer, public, save :: & 
     IDO3_O1D  = -1,IDO3_O3P = -1,IDH2O2 = -1,IDNO2     = -1,IDNO3_NO  = -1 &
    ,IDNO3_NO2 = -1,IDHONO   = -1,IDHNO3 = -1,IDHCHO_H  = -1,IDHCHO_H2 = -1 &
    ,IDCH3CHO  = -1 &
    ,MCM_J15   = -1, MCM_J17  = -1, MCM_J18 = -1, MCM_J20 = -1 &
    ,MCM_J22   = -1, IDMEK    = -1, MCM_J54 = -1, MCM_J51 = -1 &
    ,MCM_J23   = -1, MCM_J52  = -1, MCM_J56 = -1, MCM_J53 = -1 & 
    ,IDACETON  = -1, IDCH3COY = -1   & ! Based on CH3COY = BIACET (in OSLO CTM3)

    ,IDCHOCHO_2CO  = -1 &! 
    ,IDCHOCHO_HCHO = -1 &! 
    ,IDCHOCHO_2CHO = -1 &!

    ,IDGLYOXA = -1 & ! these are sometimes used in place of the above three; same reactions
    ,IDGLYOXB = -1 &
    ,IDGLYOXC = -1 &

    ,IDCHOCHO      = -1 & ! sometimes present as lumped reaction
    ,IDNO3         = -1 & ! sometimes present as lumped reaction

    ,IDCH3COCH3 = -1, IDRCOCHO  = -1 &
    ,IDCH3O2H   = -1 &
    ,IDETP      = -1, IDETHP    = -1, IDATOOH    = -1  &  
    ,IDR4P      = -1, IDRIPC    = -1, IDIDHPE    = -1  &  
    ,IDINPD     = -1, IDMAP     = -1, IDRP       = -1  &  
    ,IDPRALDP   = -1, IDITCN    = -1, IDPIP      = -1  &  

    ,IDiC3H7ONO2 = -1 & ! Used for ISON
    ,IDCH3COO2H  = -1 & 
    ,IDMVKN      = -1 & ! MCM 56 scaled by 1.6 following CRI R1A, unique in CloudJ
    ,IDR4N2      = -1 & ! following CRI R1A notes, unique in CloudJ
    ,IDHO2NO2    = -1 & ! IDHO2NO2 not in MCM -- use scaled MEK photolysis rate (assuming jHO2NO2 can be approximated by 1.8*jMEK)
    ,IDN2O5      = -1 & ! IDN2O5 not in MCM -- use scaled H2O2 photolysis rate (assuming jN2O5 can be approximated by 7*jH2O2)
    ,IDC2H5CHO   = -1 & 
    ,IDACETOL    = -1 & ! https://mcm.york.ac.uk/MCM/species/ACETOL   (J22)
    ,IDINPB      = -1 & ! IDMEK in CRI R1A, unique in CloudJ
    ,IDIHN3      = -1 & ! MCM53 based on CRI R1A, unique in CloudJ
  
    ,IDGLYALD    = -1 & ! https://mcm.york.ac.uk/MCM/species/HOCH2CHO (J15 in MCM, but unique in CloudJ)
    ,IDMCRENOL   = -1 & ! same MCM J15 as for GLYALD following RB (but unique in CloudJ)
    ,IDICN       = -1 & ! https://acp.copernicus.org/preprints/acp-2019-328/acp-2019-328-supplement.pdf (ICN lumped isoprene caronyl nitrate isomers)
                        ! photolysis major sink for isoprene-derived carbonyl groups, doi:10.5194/acp-14-2497-2014
                        ! Bergstrom 2022 scaled MEK by 0.778.

    ,IDPAN = -1  

end module DefPhotolysis_mod