! <CM_ChemGroups_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_10(3282)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2016 met.no
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

integer, public, target, save, dimension(27) :: &
  OM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,OM25_BGND,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

integer, public, target, save, dimension(2) :: &
  DDEP_DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(3) :: &
  PPM25_FIRE_GROUP = (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25 /)

integer, public, target, save, dimension(1) :: &
  WDEP_SVWOODOA25_GROUP = (/ WOODOA_NG10 /)

integer, public, target, save, dimension(15) :: &
  DDEP_PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,ASH_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(57) :: &
  NMVOC_GROUP = (/ PAN,MPAN,CH3COO2,MACR,GLYOX,MGLYOX,MAL,MEK,MVK,HCHO,CH3CHO,C2H6,NC4H10,C2H4,C3H6,OXYL,C5H8,APINENE,CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,MEKO2H,MALO2H,MACROOH,MACO3H,MACO2H,CH3COO2H,CH3OH,C2H5OH,ACETOL,POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,OM25_BGND,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

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
  PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,ASH_F,ASH_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,OM25_P,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(23) :: &
  DDEP_OM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

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

integer, public, target, save, dimension(11) :: &
  DDEP_NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC /)

integer, public, target, save, dimension(1) :: &
  WDEP_NVFFIREOC25_GROUP = (/ FFIRE_OM /)

integer, public, target, save, dimension(2) :: &
  SOX_GROUP = (/ SO2,SO4 /)

integer, public, target, save, dimension(1) :: &
  DUST_ANT_F_GROUP = (/ DUST_ROAD_F /)

integer, public, target, save, dimension(16) :: &
  PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,ASH_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,OM25_P,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(3) :: &
  DDEP_NVABSOM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM /)

integer, public, target, save, dimension(1) :: &
  WDEP_SVFFIREOA25_GROUP = (/ FFIREOA_NG10 /)

integer, public, target, save, dimension(12) :: &
  NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,OM25_BGND /)

integer, public, target, save, dimension(1) :: &
  WDEP_DAOBS_GROUP = (/ SO2 /)

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
  DUST_NAT_F_GROUP = (/ DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(2) :: &
  SS_GROUP = (/ SEASALT_F,SEASALT_C /)

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
  AOD_GROUP = (/ ASH_C,ASH_F,DUST_ROAD_C,DUST_ROAD_F,DUST_SAH_C,DUST_SAH_F,DUST_WB_C,DUST_WB_F,EC_F_FFUEL_AGE,EC_F_FFUEL_NEW,EC_F_WOOD_AGE,EC_F_WOOD_NEW,FFIRE_BC,FFIRE_OM,FFIRE_REMPPM25,NH4_F,NO3_C,NO3_F,OM25_P,REMPPM25,REMPPM_C,SEASALT_C,SEASALT_F,SO4 /)

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

integer, public, target, save, dimension(11) :: &
  WDEP_NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC /)

integer, public, target, save, dimension(7) :: &
  PPM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

integer, public, target, save, dimension(10) :: &
  ASOA_GROUP = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

integer, public, target, save, dimension(4) :: &
  NVABSOM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,OM25_BGND /)

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
  GRPS_GROUP = (/ DUMMY /)

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

integer, public, target, save, dimension(1) :: &
  WDEP_OMCOARSE_GROUP = (/ POM_C_FFUEL /)

integer, public, target, save, dimension(2) :: &
  TNO3_GROUP = (/ NO3_F,NO3_C /)

integer, public, target, save, dimension(1) :: &
  SVFFIREOA25_GROUP = (/ FFIREOA_NG10 /)

integer, public, target, save, dimension(3) :: &
  DDEP_DAOBS_GROUP = (/ O3,NO2,SO2 /)

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
  DAOBS_GROUP = (/ O3,NO2,SO2 /)

integer, public, target, save, dimension(3) :: &
  WDEP_FFUELEC_GROUP = (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL /)

integer, public, target, save, dimension(26) :: &
  WDEP_OM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)


! ------- RO2 Pool     species ------------------
integer, public, parameter :: SIZE_RO2_POOL      = 1
integer, public, parameter, dimension(1) :: &
  RO2_POOL      = (/ -99 /)
type(typ_sp), dimension(138), public, save :: chemgroups


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

  chemgroups(8)%name="OM25"
  chemgroups(8)%ptr=>OM25_GROUP

  chemgroups(9)%name="DDEP_DUST_NAT_C"
  chemgroups(9)%ptr=>DDEP_DUST_NAT_C_GROUP

  chemgroups(10)%name="PPM25_FIRE"
  chemgroups(10)%ptr=>PPM25_FIRE_GROUP

  chemgroups(11)%name="WDEP_SVWOODOA25"
  chemgroups(11)%ptr=>WDEP_SVWOODOA25_GROUP

  chemgroups(12)%name="DDEP_PMFINE"
  chemgroups(12)%ptr=>DDEP_PMFINE_GROUP

  chemgroups(13)%name="NMVOC"
  chemgroups(13)%ptr=>NMVOC_GROUP

  chemgroups(14)%name="WDEP_PPM25"
  chemgroups(14)%ptr=>WDEP_PPM25_GROUP

  chemgroups(15)%name="WDEP_PM10"
  chemgroups(15)%ptr=>WDEP_PM10_GROUP

  chemgroups(16)%name="DDEP_OX"
  chemgroups(16)%ptr=>DDEP_OX_GROUP

  chemgroups(17)%name="DDEP_ECCOARSE"
  chemgroups(17)%ptr=>DDEP_ECCOARSE_GROUP

  chemgroups(18)%name="NVWOODOC25"
  chemgroups(18)%ptr=>NVWOODOC25_GROUP

  chemgroups(19)%name="ASH"
  chemgroups(19)%ptr=>ASH_GROUP

  chemgroups(20)%name="WDEP_NVFFUELOC_COARSE"
  chemgroups(20)%ptr=>WDEP_NVFFUELOC_COARSE_GROUP

  chemgroups(21)%name="PM10"
  chemgroups(21)%ptr=>PM10_GROUP

  chemgroups(22)%name="DDEP_OM25"
  chemgroups(22)%ptr=>DDEP_OM25_GROUP

  chemgroups(23)%name="WDEP_PWOODOA25"
  chemgroups(23)%ptr=>WDEP_PWOODOA25_GROUP

  chemgroups(24)%name="DDEP_ASOA"
  chemgroups(24)%ptr=>DDEP_ASOA_GROUP

  chemgroups(25)%name="OX"
  chemgroups(25)%ptr=>OX_GROUP

  chemgroups(26)%name="DDEP_OXN"
  chemgroups(26)%ptr=>DDEP_OXN_GROUP

  chemgroups(27)%name="WDEP_PFFUELOA25"
  chemgroups(27)%ptr=>WDEP_PFFUELOA25_GROUP

  chemgroups(28)%name="DDEP_PPM10"
  chemgroups(28)%ptr=>DDEP_PPM10_GROUP

  chemgroups(29)%name="WDEP_PPM_C"
  chemgroups(29)%ptr=>WDEP_PPM_C_GROUP

  chemgroups(30)%name="DDEP_PPM25_FIRE"
  chemgroups(30)%ptr=>DDEP_PPM25_FIRE_GROUP

  chemgroups(31)%name="WDEP_EC_F"
  chemgroups(31)%ptr=>WDEP_EC_F_GROUP

  chemgroups(32)%name="DDEP_PM10"
  chemgroups(32)%ptr=>DDEP_PM10_GROUP

  chemgroups(33)%name="WDEP_SIA"
  chemgroups(33)%ptr=>WDEP_SIA_GROUP

  chemgroups(34)%name="BVOC"
  chemgroups(34)%ptr=>BVOC_GROUP

  chemgroups(35)%name="PWOODOA25"
  chemgroups(35)%ptr=>PWOODOA25_GROUP

  chemgroups(36)%name="EC_F"
  chemgroups(36)%ptr=>EC_F_GROUP

  chemgroups(37)%name="DDEP_NONVOLPCM"
  chemgroups(37)%ptr=>DDEP_NONVOLPCM_GROUP

  chemgroups(38)%name="WDEP_NVFFIREOC25"
  chemgroups(38)%ptr=>WDEP_NVFFIREOC25_GROUP

  chemgroups(39)%name="SOX"
  chemgroups(39)%ptr=>SOX_GROUP

  chemgroups(40)%name="DUST_ANT_F"
  chemgroups(40)%ptr=>DUST_ANT_F_GROUP

  chemgroups(41)%name="PMFINE"
  chemgroups(41)%ptr=>PMFINE_GROUP

  chemgroups(42)%name="DDEP_NVABSOM"
  chemgroups(42)%ptr=>DDEP_NVABSOM_GROUP

  chemgroups(43)%name="WDEP_SVFFIREOA25"
  chemgroups(43)%ptr=>WDEP_SVFFIREOA25_GROUP

  chemgroups(44)%name="NONVOLPCM"
  chemgroups(44)%ptr=>NONVOLPCM_GROUP

  chemgroups(45)%name="WDEP_DAOBS"
  chemgroups(45)%ptr=>WDEP_DAOBS_GROUP

  chemgroups(46)%name="WDEP_DUST_NAT_C"
  chemgroups(46)%ptr=>WDEP_DUST_NAT_C_GROUP

  chemgroups(47)%name="PMCO"
  chemgroups(47)%ptr=>PMCO_GROUP

  chemgroups(48)%name="DDEP_DUST_ANT_F"
  chemgroups(48)%ptr=>DDEP_DUST_ANT_F_GROUP

  chemgroups(49)%name="DDEP_OMCOARSE"
  chemgroups(49)%ptr=>DDEP_OMCOARSE_GROUP

  chemgroups(50)%name="NVFFIREOC25"
  chemgroups(50)%ptr=>NVFFIREOC25_GROUP

  chemgroups(51)%name="WDEP_RDN"
  chemgroups(51)%ptr=>WDEP_RDN_GROUP

  chemgroups(52)%name="WDEP_ASOA"
  chemgroups(52)%ptr=>WDEP_ASOA_GROUP

  chemgroups(53)%name="WDEP_NVWOODOC25"
  chemgroups(53)%ptr=>WDEP_NVWOODOC25_GROUP

  chemgroups(54)%name="DDEP_NVWOODOC25"
  chemgroups(54)%ptr=>DDEP_NVWOODOC25_GROUP

  chemgroups(55)%name="WDEP_ECFINE"
  chemgroups(55)%ptr=>WDEP_ECFINE_GROUP

  chemgroups(56)%name="WDEP_ECCOARSE"
  chemgroups(56)%ptr=>WDEP_ECCOARSE_GROUP

  chemgroups(57)%name="OXN"
  chemgroups(57)%ptr=>OXN_GROUP

  chemgroups(58)%name="WDEP_DUST_ANT_F"
  chemgroups(58)%ptr=>WDEP_DUST_ANT_F_GROUP

  chemgroups(59)%name="DDEP_NVFFUELOC25"
  chemgroups(59)%ptr=>DDEP_NVFFUELOC25_GROUP

  chemgroups(60)%name="SIA"
  chemgroups(60)%ptr=>SIA_GROUP

  chemgroups(61)%name="DDEP_ASH"
  chemgroups(61)%ptr=>DDEP_ASH_GROUP

  chemgroups(62)%name="DDEP_NVFFIREOC25"
  chemgroups(62)%ptr=>DDEP_NVFFIREOC25_GROUP

  chemgroups(63)%name="DDEP_TNO3"
  chemgroups(63)%ptr=>DDEP_TNO3_GROUP

  chemgroups(64)%name="DDEP_PMCO"
  chemgroups(64)%ptr=>DDEP_PMCO_GROUP

  chemgroups(65)%name="DDEP_BSOA"
  chemgroups(65)%ptr=>DDEP_BSOA_GROUP

  chemgroups(66)%name="DDEP_PCM"
  chemgroups(66)%ptr=>DDEP_PCM_GROUP

  chemgroups(67)%name="DDEP_NVFFUELOC_COARSE"
  chemgroups(67)%ptr=>DDEP_NVFFUELOC_COARSE_GROUP

  chemgroups(68)%name="ECCOARSE"
  chemgroups(68)%ptr=>ECCOARSE_GROUP

  chemgroups(69)%name="WOODEC"
  chemgroups(69)%ptr=>WOODEC_GROUP

  chemgroups(70)%name="WDEP_TNO3"
  chemgroups(70)%ptr=>WDEP_TNO3_GROUP

  chemgroups(71)%name="DDEP_DUST_ANT_C"
  chemgroups(71)%ptr=>DDEP_DUST_ANT_C_GROUP

  chemgroups(72)%name="WDEP_PMCO"
  chemgroups(72)%ptr=>WDEP_PMCO_GROUP

  chemgroups(73)%name="DDEP_FFUELEC"
  chemgroups(73)%ptr=>DDEP_FFUELEC_GROUP

  chemgroups(74)%name="WDEP_NVFFUELOC25"
  chemgroups(74)%ptr=>WDEP_NVFFUELOC25_GROUP

  chemgroups(75)%name="WDEP_PM10ANTHR"
  chemgroups(75)%ptr=>WDEP_PM10ANTHR_GROUP

  chemgroups(76)%name="WDEP_SVFFUELOA25"
  chemgroups(76)%ptr=>WDEP_SVFFUELOA25_GROUP

  chemgroups(77)%name="DDEP_PPM25"
  chemgroups(77)%ptr=>DDEP_PPM25_GROUP

  chemgroups(78)%name="WDEP_PPM25_FIRE"
  chemgroups(78)%ptr=>WDEP_PPM25_FIRE_GROUP

  chemgroups(79)%name="FFIREBC"
  chemgroups(79)%ptr=>FFIREBC_GROUP

  chemgroups(80)%name="WDEP_FFIREBC"
  chemgroups(80)%ptr=>WDEP_FFIREBC_GROUP

  chemgroups(81)%name="NOX"
  chemgroups(81)%ptr=>NOX_GROUP

  chemgroups(82)%name="DUST_NAT_F"
  chemgroups(82)%ptr=>DUST_NAT_F_GROUP

  chemgroups(83)%name="SS"
  chemgroups(83)%ptr=>SS_GROUP

  chemgroups(84)%name="DDEP_DUST"
  chemgroups(84)%ptr=>DDEP_DUST_GROUP

  chemgroups(85)%name="RO2"
  chemgroups(85)%ptr=>RO2_GROUP

  chemgroups(86)%name="DDEP_EC_F"
  chemgroups(86)%ptr=>DDEP_EC_F_GROUP

  chemgroups(87)%name="ROOH"
  chemgroups(87)%ptr=>ROOH_GROUP

  chemgroups(88)%name="DUST_ANT_C"
  chemgroups(88)%ptr=>DUST_ANT_C_GROUP

  chemgroups(89)%name="AOD"
  chemgroups(89)%ptr=>AOD_GROUP

  chemgroups(90)%name="DDEP_PPM_C"
  chemgroups(90)%ptr=>DDEP_PPM_C_GROUP

  chemgroups(91)%name="WDEP_NVABSOM"
  chemgroups(91)%ptr=>WDEP_NVABSOM_GROUP

  chemgroups(92)%name="SVFFUELOA25"
  chemgroups(92)%ptr=>SVFFUELOA25_GROUP

  chemgroups(93)%name="DDEP_SOX"
  chemgroups(93)%ptr=>DDEP_SOX_GROUP

  chemgroups(94)%name="WDEP_DUST_NAT_F"
  chemgroups(94)%ptr=>WDEP_DUST_NAT_F_GROUP

  chemgroups(95)%name="WDEP_SOX"
  chemgroups(95)%ptr=>WDEP_SOX_GROUP

  chemgroups(96)%name="PFFUELOA25"
  chemgroups(96)%ptr=>PFFUELOA25_GROUP

  chemgroups(97)%name="FFUELEC"
  chemgroups(97)%ptr=>FFUELEC_GROUP

  chemgroups(98)%name="PM10ANTHR"
  chemgroups(98)%ptr=>PM10ANTHR_GROUP

  chemgroups(99)%name="WDEP_DUST_ANT_C"
  chemgroups(99)%ptr=>WDEP_DUST_ANT_C_GROUP

  chemgroups(100)%name="DDEP_DUST_NAT_F"
  chemgroups(100)%ptr=>DDEP_DUST_NAT_F_GROUP

  chemgroups(101)%name="DDEP_PM10ANTHR"
  chemgroups(101)%ptr=>DDEP_PM10ANTHR_GROUP

  chemgroups(102)%name="DDEP_WOODEC"
  chemgroups(102)%ptr=>DDEP_WOODEC_GROUP

  chemgroups(103)%name="WDEP_WOODEC"
  chemgroups(103)%ptr=>WDEP_WOODEC_GROUP

  chemgroups(104)%name="WDEP_PCM"
  chemgroups(104)%ptr=>WDEP_PCM_GROUP

  chemgroups(105)%name="WDEP_SS"
  chemgroups(105)%ptr=>WDEP_SS_GROUP

  chemgroups(106)%name="DUST_NAT_C"
  chemgroups(106)%ptr=>DUST_NAT_C_GROUP

  chemgroups(107)%name="WDEP_NONVOLPCM"
  chemgroups(107)%ptr=>WDEP_NONVOLPCM_GROUP

  chemgroups(108)%name="PPM25"
  chemgroups(108)%ptr=>PPM25_GROUP

  chemgroups(109)%name="ASOA"
  chemgroups(109)%ptr=>ASOA_GROUP

  chemgroups(110)%name="NVABSOM"
  chemgroups(110)%ptr=>NVABSOM_GROUP

  chemgroups(111)%name="PPM10"
  chemgroups(111)%ptr=>PPM10_GROUP

  chemgroups(112)%name="BSOA"
  chemgroups(112)%ptr=>BSOA_GROUP

  chemgroups(113)%name="ECFINE"
  chemgroups(113)%ptr=>ECFINE_GROUP

  chemgroups(114)%name="DDEP_FFIREBC"
  chemgroups(114)%ptr=>DDEP_FFIREBC_GROUP

  chemgroups(115)%name="WDEP_DUST"
  chemgroups(115)%ptr=>WDEP_DUST_GROUP

  chemgroups(116)%name="GRPS"
  chemgroups(116)%ptr=>GRPS_GROUP

  chemgroups(117)%name="NVFFUELOC_COARSE"
  chemgroups(117)%ptr=>NVFFUELOC_COARSE_GROUP

  chemgroups(118)%name="DDEP_ECFINE"
  chemgroups(118)%ptr=>DDEP_ECFINE_GROUP

  chemgroups(119)%name="WDEP_ROOH"
  chemgroups(119)%ptr=>WDEP_ROOH_GROUP

  chemgroups(120)%name="PCM"
  chemgroups(120)%ptr=>PCM_GROUP

  chemgroups(121)%name="SVWOODOA25"
  chemgroups(121)%ptr=>SVWOODOA25_GROUP

  chemgroups(122)%name="DDEP_SIA"
  chemgroups(122)%ptr=>DDEP_SIA_GROUP

  chemgroups(123)%name="WDEP_PFFIREOA25"
  chemgroups(123)%ptr=>WDEP_PFFIREOA25_GROUP

  chemgroups(124)%name="WDEP_OMCOARSE"
  chemgroups(124)%ptr=>WDEP_OMCOARSE_GROUP

  chemgroups(125)%name="TNO3"
  chemgroups(125)%ptr=>TNO3_GROUP

  chemgroups(126)%name="SVFFIREOA25"
  chemgroups(126)%ptr=>SVFFIREOA25_GROUP

  chemgroups(127)%name="DDEP_DAOBS"
  chemgroups(127)%ptr=>DDEP_DAOBS_GROUP

  chemgroups(128)%name="NVFFUELOC25"
  chemgroups(128)%ptr=>NVFFUELOC25_GROUP

  chemgroups(129)%name="OMCOARSE"
  chemgroups(129)%ptr=>OMCOARSE_GROUP

  chemgroups(130)%name="DDEP_ROOH"
  chemgroups(130)%ptr=>DDEP_ROOH_GROUP

  chemgroups(131)%name="WDEP_PMFINE"
  chemgroups(131)%ptr=>WDEP_PMFINE_GROUP

  chemgroups(132)%name="DDEP_RDN"
  chemgroups(132)%ptr=>DDEP_RDN_GROUP

  chemgroups(133)%name="PFFIREOA25"
  chemgroups(133)%ptr=>PFFIREOA25_GROUP

  chemgroups(134)%name="RDN"
  chemgroups(134)%ptr=>RDN_GROUP

  chemgroups(135)%name="WDEP_ASH"
  chemgroups(135)%ptr=>WDEP_ASH_GROUP

  chemgroups(136)%name="DAOBS"
  chemgroups(136)%ptr=>DAOBS_GROUP

  chemgroups(137)%name="WDEP_FFUELEC"
  chemgroups(137)%ptr=>WDEP_FFUELEC_GROUP

  chemgroups(138)%name="WDEP_OM25"
  chemgroups(138)%ptr=>WDEP_OM25_GROUP

endsubroutine Init_ChemGroups
 !-----------------------------------------------------------
endmodule ChemGroups_ml
 !-----------------------------------------------------------
