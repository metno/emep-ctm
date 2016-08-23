! <CM_ChemGroups_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2015 met.no
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
!>_________________________________________________________<

module ChemGroups_ml
!-----------------------------------------------------------

use ChemSpecs_tot_ml  ! => species indices
use OwnDataTypes_ml   ! => typ_sp
implicit none
private
! Assignment of groups from GenIn.species:
public :: Init_ChemGroups

integer, public, target, save, dimension(2) :: &
  DDEP_SS_GROUP = (/ SEASALT_F,SEASALT_C /)

integer, public, target, save, dimension(4) :: &
  WDEP_OXN_GROUP = (/ HNO3,HONO,NO3_F,NO3_C /)

integer, public, target, save, dimension(11) :: &
  WDEP_PPM10_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C /)

integer, public, target, save, dimension(6) :: &
  DUST_GROUP = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(10) :: &
  WDEP_BSOA_GROUP = (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

integer, public, target, save, dimension(2) :: &
  DDEP_NOX_GROUP = (/ NO2,SHIPNOX /)

integer, public, target, save, dimension(4) :: &
  PPM_C_GROUP = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C /)

integer, public, target, save, dimension(2) :: &
  DDEP_DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(3) :: &
  PPM25_FIRE_GROUP = (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25 /)

integer, public, target, save, dimension(1) :: &
  WDEP_SVWOODOA25_GROUP = (/ WOODOA_NG10 /)

integer, public, target, save, dimension(15) :: &
  DDEP_PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,ASH_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(70) :: &
  NMVOC_GROUP = (/ PAN,MPAN,CH3COO2,MACR,GLYOX,MGLYOX,MAL,MEK,MVK,HCHO,CH3CHO,C2H6,NC4H10,C2H4,C3H6,OXYL,C5H8,APINENE,CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,MEKO2H,MALO2H,MACROOH,MACO3H,MACO2H,CH3COO2H,CH3OH,C2H5OH,ACETOL,GAS_ASOA_OC,PART_ASOA_OC,GAS_BSOA_OC,PART_BSOA_OC,PART_FFUELOA25_OC,PART_WOODOA25_OC,PART_FFIREOA25_OC,PART_OC10,PART_OC25,NONVOL_FFUELOC25,NONV_FFUELOC_COARSE,NONVOL_WOODOC25,NONVOL_BGNDOC,NONVOL_FFIREOC25,POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

integer, public, target, save, dimension(7) :: &
  WDEP_PPM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

integer, public, target, save, dimension(25) :: &
  WDEP_PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,ASH_F,ASH_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(3) :: &
  DDEP_OX_GROUP = (/ O3,NO2,SHIPNOX /)

integer, public, target, save, dimension(2) :: &
  DDEP_ECCOARSE_GROUP = (/ EC_C_WOOD,EC_C_FFUEL /)

integer, public, target, save, dimension(1) :: &
  NVWOODOC25_GROUP = (/ POM_F_WOOD /)

integer, public, target, save, dimension(2) :: &
  ASH_GROUP = (/ ASH_F,ASH_C /)

integer, public, target, save, dimension(1) :: &
  WDEP_NVFFUELOC_COARSE_GROUP = (/ POM_C_FFUEL /)

integer, public, target, save, dimension(26) :: &
  PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,ASH_F,ASH_C,PART_OM_F,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(1) :: &
  WDEP_PWOODOA25_GROUP = (/ WOODOA_NG10 /)

integer, public, target, save, dimension(10) :: &
  DDEP_ASOA_GROUP = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

integer, public, target, save, dimension(3) :: &
  OX_GROUP = (/ O3,NO2,SHIPNOX /)

integer, public, target, save, dimension(8) :: &
  DDEP_OXN_GROUP = (/ NO2,SHIPNOX,PAN,MPAN,HNO3,HONO,NO3_F,NO3_C /)

integer, public, target, save, dimension(1) :: &
  WDEP_PFFUELOA25_GROUP = (/ FFFUEL_NG10 /)

integer, public, target, save, dimension(11) :: &
  DDEP_PPM10_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C /)

integer, public, target, save, dimension(4) :: &
  WDEP_PPM_C_GROUP = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C /)

integer, public, target, save, dimension(3) :: &
  DDEP_PPM25_FIRE_GROUP = (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25 /)

integer, public, target, save, dimension(5) :: &
  WDEP_EC_F_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC /)

integer, public, target, save, dimension(25) :: &
  DDEP_PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,ASH_F,ASH_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(4) :: &
  WDEP_SIA_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F /)

integer, public, target, save, dimension(2) :: &
  BVOC_GROUP = (/ C5H8,APINENE /)

integer, public, target, save, dimension(1) :: &
  PWOODOA25_GROUP = (/ WOODOA_NG10 /)

integer, public, target, save, dimension(5) :: &
  EC_F_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC /)

integer, public, target, save, dimension(10) :: &
  DDEP_NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC /)

integer, public, target, save, dimension(1) :: &
  WDEP_NVFFIREOC25_GROUP = (/ FFIRE_OM /)

integer, public, target, save, dimension(2) :: &
  SOX_GROUP = (/ SO2,SO4 /)

integer, public, target, save, dimension(1) :: &
  DUST_ANT_F_GROUP = (/ DUST_ROAD_F /)

integer, public, target, save, dimension(16) :: &
  PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,ASH_F,PART_OM_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(3) :: &
  DDEP_NVABSOM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM /)

integer, public, target, save, dimension(1) :: &
  WDEP_SVFFIREOA25_GROUP = (/ FFIREOA_NG10 /)

integer, public, target, save, dimension(10) :: &
  NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC /)

integer, public, target, save, dimension(2) :: &
  WDEP_DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(10) :: &
  PMCO_GROUP = (/ NO3_C,ASH_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(1) :: &
  DDEP_DUST_ANT_F_GROUP = (/ DUST_ROAD_F /)

integer, public, target, save, dimension(1) :: &
  DDEP_OMCOARSE_GROUP = (/ POM_C_FFUEL /)

integer, public, target, save, dimension(1) :: &
  NVFFIREOC25_GROUP = (/ FFIRE_OM /)

integer, public, target, save, dimension(2) :: &
  WDEP_RDN_GROUP = (/ NH3,NH4_F /)

integer, public, target, save, dimension(10) :: &
  WDEP_ASOA_GROUP = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

integer, public, target, save, dimension(1) :: &
  WDEP_NVWOODOC25_GROUP = (/ POM_F_WOOD /)

integer, public, target, save, dimension(1) :: &
  DDEP_NVWOODOC25_GROUP = (/ POM_F_WOOD /)

integer, public, target, save, dimension(4) :: &
  WDEP_ECFINE_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE /)

integer, public, target, save, dimension(2) :: &
  WDEP_ECCOARSE_GROUP = (/ EC_C_WOOD,EC_C_FFUEL /)

integer, public, target, save, dimension(14) :: &
  OXN_GROUP = (/ NO,NO2,SHIPNOX,PAN,MPAN,NO3,N2O5,ISONO3,HNO3,HONO,ISNI,ISNIR,NO3_F,NO3_C /)

integer, public, target, save, dimension(1) :: &
  WDEP_DUST_ANT_F_GROUP = (/ DUST_ROAD_F /)

integer, public, target, save, dimension(1) :: &
  DDEP_NVFFUELOC25_GROUP = (/ POM_F_FFUEL /)

integer, public, target, save, dimension(4) :: &
  SIA_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F /)

integer, public, target, save, dimension(2) :: &
  DDEP_ASH_GROUP = (/ ASH_F,ASH_C /)

integer, public, target, save, dimension(1) :: &
  DDEP_NVFFIREOC25_GROUP = (/ FFIRE_OM /)

integer, public, target, save, dimension(2) :: &
  DDEP_TNO3_GROUP = (/ NO3_F,NO3_C /)

integer, public, target, save, dimension(10) :: &
  DDEP_PMCO_GROUP = (/ NO3_C,ASH_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(10) :: &
  DDEP_BSOA_GROUP = (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

integer, public, target, save, dimension(28) :: &
  DDEP_PCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

integer, public, target, save, dimension(1) :: &
  DDEP_NVFFUELOC_COARSE_GROUP = (/ POM_C_FFUEL /)

integer, public, target, save, dimension(2) :: &
  ECCOARSE_GROUP = (/ EC_C_WOOD,EC_C_FFUEL /)

integer, public, target, save, dimension(3) :: &
  WOODEC_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD /)

integer, public, target, save, dimension(2) :: &
  WDEP_TNO3_GROUP = (/ NO3_F,NO3_C /)

integer, public, target, save, dimension(1) :: &
  DDEP_DUST_ANT_C_GROUP = (/ DUST_ROAD_C /)

integer, public, target, save, dimension(10) :: &
  WDEP_PMCO_GROUP = (/ NO3_C,ASH_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(3) :: &
  DDEP_FFUELEC_GROUP = (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL /)

integer, public, target, save, dimension(1) :: &
  WDEP_NVFFUELOC25_GROUP = (/ POM_F_FFUEL /)

integer, public, target, save, dimension(3) :: &
  WDEP_PM10ANTHR_GROUP = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL /)

integer, public, target, save, dimension(1) :: &
  WDEP_SVFFUELOA25_GROUP = (/ FFFUEL_NG10 /)

integer, public, target, save, dimension(7) :: &
  DDEP_PPM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

integer, public, target, save, dimension(3) :: &
  WDEP_PPM25_FIRE_GROUP = (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25 /)

integer, public, target, save, dimension(1) :: &
  FFIREBC_GROUP = (/ FFIRE_BC /)

integer, public, target, save, dimension(1) :: &
  WDEP_FFIREBC_GROUP = (/ FFIRE_BC /)

integer, public, target, save, dimension(3) :: &
  NOX_GROUP = (/ NO,NO2,SHIPNOX /)

integer, public, target, save, dimension(2) :: &
  SS_GROUP = (/ SEASALT_F,SEASALT_C /)

integer, public, target, save, dimension(2) :: &
  DUST_NAT_F_GROUP = (/ DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(6) :: &
  DDEP_DUST_GROUP = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(14) :: &
  RO2_GROUP = (/ HO2,CH3O2,C2H5O2,SECC4H9O2,ISRO2,ETRO2,PRRO2,OXYO2,MEKO2,MALO2,MVKO2,MACRO2,MACO3,TERPPEROXY /)

integer, public, target, save, dimension(5) :: &
  DDEP_EC_F_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC /)

integer, public, target, save, dimension(16) :: &
  ROOH_GROUP = (/ CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,OXYO2H,MEKO2H,MALO2H,MVKO2H,MACROOH,MACO3H,ISRO2H,H2O2,CH3COO2H,ISONO3H,ISNIRH /)

integer, public, target, save, dimension(1) :: &
  DUST_ANT_C_GROUP = (/ DUST_ROAD_C /)

integer, public, target, save, dimension(24) :: &
  AOD_GROUP = (/ ASH_C,ASH_F,DUST_ROAD_C,DUST_ROAD_F,DUST_SAH_C,DUST_SAH_F,DUST_WB_C,DUST_WB_F,EC_F_FFUEL_AGE,EC_F_FFUEL_NEW,EC_F_WOOD_AGE,EC_F_WOOD_NEW,FFIRE_BC,FFIRE_OM,FFIRE_REMPPM25,NH4_F,NO3_C,NO3_F,PART_OM_F,REMPPM25,REMPPM_C,SEASALT_C,SEASALT_F,SO4 /)

integer, public, target, save, dimension(4) :: &
  DDEP_PPM_C_GROUP = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C /)

integer, public, target, save, dimension(3) :: &
  WDEP_NVABSOM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM /)

integer, public, target, save, dimension(1) :: &
  SVFFUELOA25_GROUP = (/ FFFUEL_NG10 /)

integer, public, target, save, dimension(2) :: &
  DDEP_SOX_GROUP = (/ SO2,SO4 /)

integer, public, target, save, dimension(2) :: &
  WDEP_DUST_NAT_F_GROUP = (/ DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(2) :: &
  WDEP_SOX_GROUP = (/ SO2,SO4 /)

integer, public, target, save, dimension(1) :: &
  PFFUELOA25_GROUP = (/ FFFUEL_NG10 /)

integer, public, target, save, dimension(3) :: &
  FFUELEC_GROUP = (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL /)

integer, public, target, save, dimension(3) :: &
  PM10ANTHR_GROUP = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL /)

integer, public, target, save, dimension(1) :: &
  WDEP_DUST_ANT_C_GROUP = (/ DUST_ROAD_C /)

integer, public, target, save, dimension(2) :: &
  DDEP_DUST_NAT_F_GROUP = (/ DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(3) :: &
  DDEP_PM10ANTHR_GROUP = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL /)

integer, public, target, save, dimension(3) :: &
  DDEP_WOODEC_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD /)

integer, public, target, save, dimension(3) :: &
  WDEP_WOODEC_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD /)

integer, public, target, save, dimension(31) :: &
  WDEP_PCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

integer, public, target, save, dimension(2) :: &
  WDEP_SS_GROUP = (/ SEASALT_F,SEASALT_C /)

integer, public, target, save, dimension(2) :: &
  DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(10) :: &
  WDEP_NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC /)

integer, public, target, save, dimension(7) :: &
  PPM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

integer, public, target, save, dimension(10) :: &
  ASOA_GROUP = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

integer, public, target, save, dimension(3) :: &
  NVABSOM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM /)

integer, public, target, save, dimension(11) :: &
  PPM10_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C /)

integer, public, target, save, dimension(10) :: &
  BSOA_GROUP = (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

integer, public, target, save, dimension(4) :: &
  ECFINE_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE /)

integer, public, target, save, dimension(1) :: &
  DDEP_FFIREBC_GROUP = (/ FFIRE_BC /)

integer, public, target, save, dimension(6) :: &
  WDEP_DUST_GROUP = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(1) :: &
  NVFFUELOC_COARSE_GROUP = (/ POM_C_FFUEL /)

integer, public, target, save, dimension(4) :: &
  DDEP_ECFINE_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE /)

integer, public, target, save, dimension(1) :: &
  WDEP_ROOH_GROUP = (/ H2O2 /)

integer, public, target, save, dimension(31) :: &
  PCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

integer, public, target, save, dimension(1) :: &
  SVWOODOA25_GROUP = (/ WOODOA_NG10 /)

integer, public, target, save, dimension(4) :: &
  DDEP_SIA_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F /)

integer, public, target, save, dimension(1) :: &
  WDEP_PFFIREOA25_GROUP = (/ FFIREOA_NG10 /)

integer, public, target, save, dimension(20) :: &
  PCM_HELP_GROUP = (/ GAS_ASOA_OC,PART_ASOA_OC,PART_ASOA_OM,GAS_BSOA_OC,PART_BSOA_OC,PART_BSOA_OM,PART_FFUELOA25_OC,PART_FFUELOA25_OM,PART_WOODOA25_OC,PART_WOODOA25_OM,PART_FFIREOA25_OC,PART_FFIREOA25_OM,PART_OC10,PART_OC25,NONVOL_FFUELOC25,NONV_FFUELOC_COARSE,NONVOL_WOODOC25,NONVOL_BGNDOC,NONVOL_FFIREOC25,PART_OM_F /)

integer, public, target, save, dimension(1) :: &
  WDEP_OMCOARSE_GROUP = (/ POM_C_FFUEL /)

integer, public, target, save, dimension(2) :: &
  TNO3_GROUP = (/ NO3_F,NO3_C /)

integer, public, target, save, dimension(1) :: &
  SVFFIREOA25_GROUP = (/ FFIREOA_NG10 /)

integer, public, target, save, dimension(3) :: &
  DDEP_DAOBS_GROUP = (/ O3,NO2,SHIPNOX /)

integer, public, target, save, dimension(1) :: &
  NVFFUELOC25_GROUP = (/ POM_F_FFUEL /)

integer, public, target, save, dimension(1) :: &
  OMCOARSE_GROUP = (/ POM_C_FFUEL /)

integer, public, target, save, dimension(3) :: &
  DDEP_ROOH_GROUP = (/ CH3O2H,C2H5OOH,H2O2 /)

integer, public, target, save, dimension(15) :: &
  WDEP_PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,ASH_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(2) :: &
  DDEP_RDN_GROUP = (/ NH3,NH4_F /)

integer, public, target, save, dimension(1) :: &
  PFFIREOA25_GROUP = (/ FFIREOA_NG10 /)

integer, public, target, save, dimension(2) :: &
  RDN_GROUP = (/ NH3,NH4_F /)

integer, public, target, save, dimension(2) :: &
  WDEP_ASH_GROUP = (/ ASH_F,ASH_C /)

integer, public, target, save, dimension(3) :: &
  DAOBS_GROUP = (/ O3,NO2,SHIPNOX /)

integer, public, target, save, dimension(3) :: &
  WDEP_FFUELEC_GROUP = (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL /)


! ------- RO2 Pool     species ------------------
integer, public, parameter :: SIZE_RO2_POOL      = 1
integer, public, parameter, dimension(1) :: &
  RO2_POOL      = (/ -99 /)
type(typ_sp), dimension(134), public, save :: chemgroups


!-----------------------------------------------------------
contains
subroutine Init_ChemGroups()

  chemgroups(1)%name="DDEP_SS"
  chemgroups(1)%ptr=>DDEP_SS_GROUP

  chemgroups(2)%name="WDEP_OXN"
  chemgroups(2)%ptr=>WDEP_OXN_GROUP

  chemgroups(3)%name="WDEP_PPM10"
  chemgroups(3)%ptr=>WDEP_PPM10_GROUP

  chemgroups(4)%name="DUST"
  chemgroups(4)%ptr=>DUST_GROUP

  chemgroups(5)%name="WDEP_BSOA"
  chemgroups(5)%ptr=>WDEP_BSOA_GROUP

  chemgroups(6)%name="DDEP_NOX"
  chemgroups(6)%ptr=>DDEP_NOX_GROUP

  chemgroups(7)%name="PPM_C"
  chemgroups(7)%ptr=>PPM_C_GROUP

  chemgroups(8)%name="DDEP_DUST_NAT_C"
  chemgroups(8)%ptr=>DDEP_DUST_NAT_C_GROUP

  chemgroups(9)%name="PPM25_FIRE"
  chemgroups(9)%ptr=>PPM25_FIRE_GROUP

  chemgroups(10)%name="WDEP_SVWOODOA25"
  chemgroups(10)%ptr=>WDEP_SVWOODOA25_GROUP

  chemgroups(11)%name="DDEP_PMFINE"
  chemgroups(11)%ptr=>DDEP_PMFINE_GROUP

  chemgroups(12)%name="NMVOC"
  chemgroups(12)%ptr=>NMVOC_GROUP

  chemgroups(13)%name="WDEP_PPM25"
  chemgroups(13)%ptr=>WDEP_PPM25_GROUP

  chemgroups(14)%name="WDEP_PM10"
  chemgroups(14)%ptr=>WDEP_PM10_GROUP

  chemgroups(15)%name="DDEP_OX"
  chemgroups(15)%ptr=>DDEP_OX_GROUP

  chemgroups(16)%name="DDEP_ECCOARSE"
  chemgroups(16)%ptr=>DDEP_ECCOARSE_GROUP

  chemgroups(17)%name="NVWOODOC25"
  chemgroups(17)%ptr=>NVWOODOC25_GROUP

  chemgroups(18)%name="ASH"
  chemgroups(18)%ptr=>ASH_GROUP

  chemgroups(19)%name="WDEP_NVFFUELOC_COARSE"
  chemgroups(19)%ptr=>WDEP_NVFFUELOC_COARSE_GROUP

  chemgroups(20)%name="PM10"
  chemgroups(20)%ptr=>PM10_GROUP

  chemgroups(21)%name="WDEP_PWOODOA25"
  chemgroups(21)%ptr=>WDEP_PWOODOA25_GROUP

  chemgroups(22)%name="DDEP_ASOA"
  chemgroups(22)%ptr=>DDEP_ASOA_GROUP

  chemgroups(23)%name="OX"
  chemgroups(23)%ptr=>OX_GROUP

  chemgroups(24)%name="DDEP_OXN"
  chemgroups(24)%ptr=>DDEP_OXN_GROUP

  chemgroups(25)%name="WDEP_PFFUELOA25"
  chemgroups(25)%ptr=>WDEP_PFFUELOA25_GROUP

  chemgroups(26)%name="DDEP_PPM10"
  chemgroups(26)%ptr=>DDEP_PPM10_GROUP

  chemgroups(27)%name="WDEP_PPM_C"
  chemgroups(27)%ptr=>WDEP_PPM_C_GROUP

  chemgroups(28)%name="DDEP_PPM25_FIRE"
  chemgroups(28)%ptr=>DDEP_PPM25_FIRE_GROUP

  chemgroups(29)%name="WDEP_EC_F"
  chemgroups(29)%ptr=>WDEP_EC_F_GROUP

  chemgroups(30)%name="DDEP_PM10"
  chemgroups(30)%ptr=>DDEP_PM10_GROUP

  chemgroups(31)%name="WDEP_SIA"
  chemgroups(31)%ptr=>WDEP_SIA_GROUP

  chemgroups(32)%name="BVOC"
  chemgroups(32)%ptr=>BVOC_GROUP

  chemgroups(33)%name="PWOODOA25"
  chemgroups(33)%ptr=>PWOODOA25_GROUP

  chemgroups(34)%name="EC_F"
  chemgroups(34)%ptr=>EC_F_GROUP

  chemgroups(35)%name="DDEP_NONVOLPCM"
  chemgroups(35)%ptr=>DDEP_NONVOLPCM_GROUP

  chemgroups(36)%name="WDEP_NVFFIREOC25"
  chemgroups(36)%ptr=>WDEP_NVFFIREOC25_GROUP

  chemgroups(37)%name="SOX"
  chemgroups(37)%ptr=>SOX_GROUP

  chemgroups(38)%name="DUST_ANT_F"
  chemgroups(38)%ptr=>DUST_ANT_F_GROUP

  chemgroups(39)%name="PMFINE"
  chemgroups(39)%ptr=>PMFINE_GROUP

  chemgroups(40)%name="DDEP_NVABSOM"
  chemgroups(40)%ptr=>DDEP_NVABSOM_GROUP

  chemgroups(41)%name="WDEP_SVFFIREOA25"
  chemgroups(41)%ptr=>WDEP_SVFFIREOA25_GROUP

  chemgroups(42)%name="NONVOLPCM"
  chemgroups(42)%ptr=>NONVOLPCM_GROUP

  chemgroups(43)%name="WDEP_DUST_NAT_C"
  chemgroups(43)%ptr=>WDEP_DUST_NAT_C_GROUP

  chemgroups(44)%name="PMCO"
  chemgroups(44)%ptr=>PMCO_GROUP

  chemgroups(45)%name="DDEP_DUST_ANT_F"
  chemgroups(45)%ptr=>DDEP_DUST_ANT_F_GROUP

  chemgroups(46)%name="DDEP_OMCOARSE"
  chemgroups(46)%ptr=>DDEP_OMCOARSE_GROUP

  chemgroups(47)%name="NVFFIREOC25"
  chemgroups(47)%ptr=>NVFFIREOC25_GROUP

  chemgroups(48)%name="WDEP_RDN"
  chemgroups(48)%ptr=>WDEP_RDN_GROUP

  chemgroups(49)%name="WDEP_ASOA"
  chemgroups(49)%ptr=>WDEP_ASOA_GROUP

  chemgroups(50)%name="WDEP_NVWOODOC25"
  chemgroups(50)%ptr=>WDEP_NVWOODOC25_GROUP

  chemgroups(51)%name="DDEP_NVWOODOC25"
  chemgroups(51)%ptr=>DDEP_NVWOODOC25_GROUP

  chemgroups(52)%name="WDEP_ECFINE"
  chemgroups(52)%ptr=>WDEP_ECFINE_GROUP

  chemgroups(53)%name="WDEP_ECCOARSE"
  chemgroups(53)%ptr=>WDEP_ECCOARSE_GROUP

  chemgroups(54)%name="OXN"
  chemgroups(54)%ptr=>OXN_GROUP

  chemgroups(55)%name="WDEP_DUST_ANT_F"
  chemgroups(55)%ptr=>WDEP_DUST_ANT_F_GROUP

  chemgroups(56)%name="DDEP_NVFFUELOC25"
  chemgroups(56)%ptr=>DDEP_NVFFUELOC25_GROUP

  chemgroups(57)%name="SIA"
  chemgroups(57)%ptr=>SIA_GROUP

  chemgroups(58)%name="DDEP_ASH"
  chemgroups(58)%ptr=>DDEP_ASH_GROUP

  chemgroups(59)%name="DDEP_NVFFIREOC25"
  chemgroups(59)%ptr=>DDEP_NVFFIREOC25_GROUP

  chemgroups(60)%name="DDEP_TNO3"
  chemgroups(60)%ptr=>DDEP_TNO3_GROUP

  chemgroups(61)%name="DDEP_PMCO"
  chemgroups(61)%ptr=>DDEP_PMCO_GROUP

  chemgroups(62)%name="DDEP_BSOA"
  chemgroups(62)%ptr=>DDEP_BSOA_GROUP

  chemgroups(63)%name="DDEP_PCM"
  chemgroups(63)%ptr=>DDEP_PCM_GROUP

  chemgroups(64)%name="DDEP_NVFFUELOC_COARSE"
  chemgroups(64)%ptr=>DDEP_NVFFUELOC_COARSE_GROUP

  chemgroups(65)%name="ECCOARSE"
  chemgroups(65)%ptr=>ECCOARSE_GROUP

  chemgroups(66)%name="WOODEC"
  chemgroups(66)%ptr=>WOODEC_GROUP

  chemgroups(67)%name="WDEP_TNO3"
  chemgroups(67)%ptr=>WDEP_TNO3_GROUP

  chemgroups(68)%name="DDEP_DUST_ANT_C"
  chemgroups(68)%ptr=>DDEP_DUST_ANT_C_GROUP

  chemgroups(69)%name="WDEP_PMCO"
  chemgroups(69)%ptr=>WDEP_PMCO_GROUP

  chemgroups(70)%name="DDEP_FFUELEC"
  chemgroups(70)%ptr=>DDEP_FFUELEC_GROUP

  chemgroups(71)%name="WDEP_NVFFUELOC25"
  chemgroups(71)%ptr=>WDEP_NVFFUELOC25_GROUP

  chemgroups(72)%name="WDEP_PM10ANTHR"
  chemgroups(72)%ptr=>WDEP_PM10ANTHR_GROUP

  chemgroups(73)%name="WDEP_SVFFUELOA25"
  chemgroups(73)%ptr=>WDEP_SVFFUELOA25_GROUP

  chemgroups(74)%name="DDEP_PPM25"
  chemgroups(74)%ptr=>DDEP_PPM25_GROUP

  chemgroups(75)%name="WDEP_PPM25_FIRE"
  chemgroups(75)%ptr=>WDEP_PPM25_FIRE_GROUP

  chemgroups(76)%name="FFIREBC"
  chemgroups(76)%ptr=>FFIREBC_GROUP

  chemgroups(77)%name="WDEP_FFIREBC"
  chemgroups(77)%ptr=>WDEP_FFIREBC_GROUP

  chemgroups(78)%name="NOX"
  chemgroups(78)%ptr=>NOX_GROUP

  chemgroups(79)%name="SS"
  chemgroups(79)%ptr=>SS_GROUP

  chemgroups(80)%name="DUST_NAT_F"
  chemgroups(80)%ptr=>DUST_NAT_F_GROUP

  chemgroups(81)%name="DDEP_DUST"
  chemgroups(81)%ptr=>DDEP_DUST_GROUP

  chemgroups(82)%name="RO2"
  chemgroups(82)%ptr=>RO2_GROUP

  chemgroups(83)%name="DDEP_EC_F"
  chemgroups(83)%ptr=>DDEP_EC_F_GROUP

  chemgroups(84)%name="ROOH"
  chemgroups(84)%ptr=>ROOH_GROUP

  chemgroups(85)%name="DUST_ANT_C"
  chemgroups(85)%ptr=>DUST_ANT_C_GROUP

  chemgroups(86)%name="AOD"
  chemgroups(86)%ptr=>AOD_GROUP

  chemgroups(87)%name="DDEP_PPM_C"
  chemgroups(87)%ptr=>DDEP_PPM_C_GROUP

  chemgroups(88)%name="WDEP_NVABSOM"
  chemgroups(88)%ptr=>WDEP_NVABSOM_GROUP

  chemgroups(89)%name="SVFFUELOA25"
  chemgroups(89)%ptr=>SVFFUELOA25_GROUP

  chemgroups(90)%name="DDEP_SOX"
  chemgroups(90)%ptr=>DDEP_SOX_GROUP

  chemgroups(91)%name="WDEP_DUST_NAT_F"
  chemgroups(91)%ptr=>WDEP_DUST_NAT_F_GROUP

  chemgroups(92)%name="WDEP_SOX"
  chemgroups(92)%ptr=>WDEP_SOX_GROUP

  chemgroups(93)%name="PFFUELOA25"
  chemgroups(93)%ptr=>PFFUELOA25_GROUP

  chemgroups(94)%name="FFUELEC"
  chemgroups(94)%ptr=>FFUELEC_GROUP

  chemgroups(95)%name="PM10ANTHR"
  chemgroups(95)%ptr=>PM10ANTHR_GROUP

  chemgroups(96)%name="WDEP_DUST_ANT_C"
  chemgroups(96)%ptr=>WDEP_DUST_ANT_C_GROUP

  chemgroups(97)%name="DDEP_DUST_NAT_F"
  chemgroups(97)%ptr=>DDEP_DUST_NAT_F_GROUP

  chemgroups(98)%name="DDEP_PM10ANTHR"
  chemgroups(98)%ptr=>DDEP_PM10ANTHR_GROUP

  chemgroups(99)%name="DDEP_WOODEC"
  chemgroups(99)%ptr=>DDEP_WOODEC_GROUP

  chemgroups(100)%name="WDEP_WOODEC"
  chemgroups(100)%ptr=>WDEP_WOODEC_GROUP

  chemgroups(101)%name="WDEP_PCM"
  chemgroups(101)%ptr=>WDEP_PCM_GROUP

  chemgroups(102)%name="WDEP_SS"
  chemgroups(102)%ptr=>WDEP_SS_GROUP

  chemgroups(103)%name="DUST_NAT_C"
  chemgroups(103)%ptr=>DUST_NAT_C_GROUP

  chemgroups(104)%name="WDEP_NONVOLPCM"
  chemgroups(104)%ptr=>WDEP_NONVOLPCM_GROUP

  chemgroups(105)%name="PPM25"
  chemgroups(105)%ptr=>PPM25_GROUP

  chemgroups(106)%name="ASOA"
  chemgroups(106)%ptr=>ASOA_GROUP

  chemgroups(107)%name="NVABSOM"
  chemgroups(107)%ptr=>NVABSOM_GROUP

  chemgroups(108)%name="PPM10"
  chemgroups(108)%ptr=>PPM10_GROUP

  chemgroups(109)%name="BSOA"
  chemgroups(109)%ptr=>BSOA_GROUP

  chemgroups(110)%name="ECFINE"
  chemgroups(110)%ptr=>ECFINE_GROUP

  chemgroups(111)%name="DDEP_FFIREBC"
  chemgroups(111)%ptr=>DDEP_FFIREBC_GROUP

  chemgroups(112)%name="WDEP_DUST"
  chemgroups(112)%ptr=>WDEP_DUST_GROUP

  chemgroups(113)%name="NVFFUELOC_COARSE"
  chemgroups(113)%ptr=>NVFFUELOC_COARSE_GROUP

  chemgroups(114)%name="DDEP_ECFINE"
  chemgroups(114)%ptr=>DDEP_ECFINE_GROUP

  chemgroups(115)%name="WDEP_ROOH"
  chemgroups(115)%ptr=>WDEP_ROOH_GROUP

  chemgroups(116)%name="PCM"
  chemgroups(116)%ptr=>PCM_GROUP

  chemgroups(117)%name="SVWOODOA25"
  chemgroups(117)%ptr=>SVWOODOA25_GROUP

  chemgroups(118)%name="DDEP_SIA"
  chemgroups(118)%ptr=>DDEP_SIA_GROUP

  chemgroups(119)%name="WDEP_PFFIREOA25"
  chemgroups(119)%ptr=>WDEP_PFFIREOA25_GROUP

  chemgroups(120)%name="PCM_HELP"
  chemgroups(120)%ptr=>PCM_HELP_GROUP

  chemgroups(121)%name="WDEP_OMCOARSE"
  chemgroups(121)%ptr=>WDEP_OMCOARSE_GROUP

  chemgroups(122)%name="TNO3"
  chemgroups(122)%ptr=>TNO3_GROUP

  chemgroups(123)%name="SVFFIREOA25"
  chemgroups(123)%ptr=>SVFFIREOA25_GROUP

  chemgroups(124)%name="DDEP_DAOBS"
  chemgroups(124)%ptr=>DDEP_DAOBS_GROUP

  chemgroups(125)%name="NVFFUELOC25"
  chemgroups(125)%ptr=>NVFFUELOC25_GROUP

  chemgroups(126)%name="OMCOARSE"
  chemgroups(126)%ptr=>OMCOARSE_GROUP

  chemgroups(127)%name="DDEP_ROOH"
  chemgroups(127)%ptr=>DDEP_ROOH_GROUP

  chemgroups(128)%name="WDEP_PMFINE"
  chemgroups(128)%ptr=>WDEP_PMFINE_GROUP

  chemgroups(129)%name="DDEP_RDN"
  chemgroups(129)%ptr=>DDEP_RDN_GROUP

  chemgroups(130)%name="PFFIREOA25"
  chemgroups(130)%ptr=>PFFIREOA25_GROUP

  chemgroups(131)%name="RDN"
  chemgroups(131)%ptr=>RDN_GROUP

  chemgroups(132)%name="WDEP_ASH"
  chemgroups(132)%ptr=>WDEP_ASH_GROUP

  chemgroups(133)%name="DAOBS"
  chemgroups(133)%ptr=>DAOBS_GROUP

  chemgroups(134)%name="WDEP_FFUELEC"
  chemgroups(134)%ptr=>WDEP_FFUELEC_GROUP

endsubroutine Init_ChemGroups
 !-----------------------------------------------------------
endmodule ChemGroups_ml
 !-----------------------------------------------------------
