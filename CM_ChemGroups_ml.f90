! <CM_ChemGroups_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2018 met.no
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

integer, public, target, save, dimension(6) :: &
  WDEP_OXN_GROUP = (/ HO2NO2,N2O5,HNO3,HONO,NO3_F,NO3_C /)

integer, public, target, save, dimension(11) :: &
  WDEP_PPM10_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C /)

integer, public, target, save, dimension(6) :: &
  DUST_GROUP = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(11) :: &
  WDEP_BSOA_GROUP = (/ SQT_SOA_NV,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

integer, public, target, save, dimension(2) :: &
  DDEP_NOX_GROUP = (/ NO2,SHIPNOX /)

integer, public, target, save, dimension(4) :: &
  PPM_C_GROUP = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C /)

integer, public, target, save, dimension(28) :: &
  OM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,OM25_BGND,SQT_SOA_NV,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

integer, public, target, save, dimension(2) :: &
  DDEP_DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(3) :: &
  PPM25_FIRE_GROUP = (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25 /)

integer, public, target, save, dimension(1) :: &
  WDEP_SVWOODOA25_GROUP = (/ WOODOA_NG10 /)

integer, public, target, save, dimension(15) :: &
  DDEP_PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,ASH_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(69) :: &
  NMVOC_GROUP = (/ PAN,CH3COO2,GLYOX,MGLYOX,MAL,MEK,HCHO,CH3CHO,C2H6,NC4H10,C2H4,C3H6,OXYL,C5H8,APINENE,BPINENE,XTERP,BIOTERP,CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,MEKO2H,MALO2H,CH3COO2H,CH3OH,C2H5OH,ACETOL,ISO2,MACRO2,MACR,MACROOH,HACET,ISOOH,ISON,HCOOH,MPAN,NALD,HPALD,PACALD,MVK,FFIRE_CO,POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,OM25_BGND,SQT_SOA_NV,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

integer, public, target, save, dimension(7) :: &
  WDEP_PPM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

integer, public, target, save, dimension(3) :: &
  DDEP_RCHO_GROUP = (/ NALD,HPALD,PACALD /)

integer, public, target, save, dimension(25) :: &
  WDEP_PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,ASH_F,ASH_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(2) :: &
  DDEP_OX_GROUP = (/ O3,NO2 /)

integer, public, target, save, dimension(2) :: &
  DDEP_ECCOARSE_GROUP = (/ EC_C_WOOD,EC_C_FFUEL /)

integer, public, target, save, dimension(1) :: &
  NVWOODOC25_GROUP = (/ POM_F_WOOD /)

integer, public, target, save, dimension(1) :: &
  WDEP_TMPOX_GROUP = (/ HO2NO2 /)

integer, public, target, save, dimension(2) :: &
  ASH_GROUP = (/ ASH_F,ASH_C /)

integer, public, target, save, dimension(1) :: &
  WDEP_NVFFUELOC_COARSE_GROUP = (/ POM_C_FFUEL /)

integer, public, target, save, dimension(26) :: &
  PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,ASH_F,ASH_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,OM25_P,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(24) :: &
  DDEP_OM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,SQT_SOA_NV,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

integer, public, target, save, dimension(1) :: &
  WDEP_PWOODOA25_GROUP = (/ WOODOA_NG10 /)

integer, public, target, save, dimension(10) :: &
  DDEP_ASOA_GROUP = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

integer, public, target, save, dimension(4) :: &
  TMPOX_GROUP = (/ O3,NO2,HO2NO2,SHIPNOX /)

integer, public, target, save, dimension(2) :: &
  OX_GROUP = (/ O3,NO2 /)

integer, public, target, save, dimension(10) :: &
  DDEP_OXN_GROUP = (/ NO2,HO2NO2,SHIPNOX,PAN,N2O5,HNO3,HONO,MPAN,NO3_F,NO3_C /)

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
  BVOC_GROUP = (/ C5H8,BIOTERP /)

integer, public, target, save, dimension(1) :: &
  PWOODOA25_GROUP = (/ WOODOA_NG10 /)

integer, public, target, save, dimension(5) :: &
  EC_F_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC /)

integer, public, target, save, dimension(12) :: &
  DDEP_NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,SQT_SOA_NV /)

integer, public, target, save, dimension(3) :: &
  RCHO_GROUP = (/ NALD,HPALD,PACALD /)

integer, public, target, save, dimension(1) :: &
  WDEP_NVFFIREOC25_GROUP = (/ FFIRE_OM /)

integer, public, target, save, dimension(2) :: &
  SOX_GROUP = (/ SO2,SO4 /)

integer, public, target, save, dimension(1) :: &
  DUST_ANT_F_GROUP = (/ DUST_ROAD_F /)

integer, public, target, save, dimension(16) :: &
  PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,ASH_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,OM25_P,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(4) :: &
  DDEP_NVABSOM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,SQT_SOA_NV /)

integer, public, target, save, dimension(1) :: &
  WDEP_SVFFIREOA25_GROUP = (/ FFIREOA_NG10 /)

integer, public, target, save, dimension(13) :: &
  NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,OM25_BGND,SQT_SOA_NV /)

integer, public, target, save, dimension(1) :: &
  WDEP_DAOBS_GROUP = (/ SO2 /)

integer, public, target, save, dimension(4) :: &
  DDEP_TMPOX_GROUP = (/ O3,NO2,HO2NO2,SHIPNOX /)

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

integer, public, target, save, dimension(12) :: &
  OXN_GROUP = (/ NO,NO2,HO2NO2,SHIPNOX,PAN,NO3,N2O5,HNO3,HONO,MPAN,NO3_F,NO3_C /)

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

integer, public, target, save, dimension(11) :: &
  DDEP_BSOA_GROUP = (/ SQT_SOA_NV,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

integer, public, target, save, dimension(29) :: &
  DDEP_PCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,SQT_SOA_NV,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

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

integer, public, target, save, dimension(15) :: &
  RO2_GROUP = (/ HO2,CH3O2,C2H5O2,SECC4H9O2,ISRO2,ETRO2,PRRO2,OXYO2,MEKO2,MALO2,MVKO2,TERPO2,XMTO3_RO2,ISO2,MACRO2 /)

integer, public, target, save, dimension(5) :: &
  DDEP_EC_F_GROUP = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC /)

integer, public, target, save, dimension(14) :: &
  ROOH_GROUP = (/ CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,OXYO2H,MEKO2H,MALO2H,H2O2,CH3COO2H,MACROOH,ISOOH,HCOOH,TERPOOH /)

integer, public, target, save, dimension(1) :: &
  DUST_ANT_C_GROUP = (/ DUST_ROAD_C /)

integer, public, target, save, dimension(24) :: &
  AOD_GROUP = (/ ASH_C,ASH_F,DUST_ROAD_C,DUST_ROAD_F,DUST_SAH_C,DUST_SAH_F,DUST_WB_C,DUST_WB_F,EC_F_FFUEL_AGE,EC_F_FFUEL_NEW,EC_F_WOOD_AGE,EC_F_WOOD_NEW,FFIRE_BC,FFIRE_OM,FFIRE_REMPPM25,NH4_F,NO3_C,NO3_F,OM25_P,REMPPM25,REMPPM_C,SEASALT_C,SEASALT_F,SO4 /)

integer, public, target, save, dimension(4) :: &
  DDEP_PPM_C_GROUP = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C /)

integer, public, target, save, dimension(4) :: &
  WDEP_NVABSOM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,SQT_SOA_NV /)

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

integer, public, target, save, dimension(32) :: &
  WDEP_PCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,SQT_SOA_NV,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

integer, public, target, save, dimension(2) :: &
  WDEP_SS_GROUP = (/ SEASALT_F,SEASALT_C /)

integer, public, target, save, dimension(2) :: &
  DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(3) :: &
  MONOTERP_GROUP = (/ APINENE,BPINENE,XTERP /)

integer, public, target, save, dimension(12) :: &
  WDEP_NONVOLPCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,SQT_SOA_NV /)

integer, public, target, save, dimension(7) :: &
  PPM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

integer, public, target, save, dimension(10) :: &
  ASOA_GROUP = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

integer, public, target, save, dimension(11) :: &
  PPM10_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C /)

integer, public, target, save, dimension(5) :: &
  NVABSOM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,OM25_BGND,SQT_SOA_NV /)

integer, public, target, save, dimension(11) :: &
  BSOA_GROUP = (/ SQT_SOA_NV,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

integer, public, target, save, dimension(2) :: &
  CHET2_GROUP = (/ IEPOX,HACET /)

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

integer, public, target, save, dimension(32) :: &
  PCM_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,SQT_SOA_NV,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

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

integer, public, target, save, dimension(4) :: &
  DDEP_ROOH_GROUP = (/ CH3O2H,C2H5OOH,H2O2,TERPOOH /)

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

integer, public, target, save, dimension(27) :: &
  WDEP_OM25_GROUP = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,SQT_SOA_NV,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)


! ------- RO2 Pool     species ------------------
integer, public, parameter :: SIZE_RO2_POOL      = 1
integer, public, parameter, dimension(1) :: &
  RO2_POOL      = (/ -99 /)
type(typ_sp), dimension(145), public, save :: chemgroups


!-----------------------------------------------------------
contains
subroutine Init_ChemGroups()

  chemgroups(1)%name="DDEP_SS"
  chemgroups(1)%specs=>DDEP_SS_GROUP

  chemgroups(2)%name="WDEP_OXN"
  chemgroups(2)%specs=>WDEP_OXN_GROUP

  chemgroups(3)%name="WDEP_PPM10"
  chemgroups(3)%specs=>WDEP_PPM10_GROUP

  chemgroups(4)%name="DUST"
  chemgroups(4)%specs=>DUST_GROUP

  chemgroups(5)%name="WDEP_BSOA"
  chemgroups(5)%specs=>WDEP_BSOA_GROUP

  chemgroups(6)%name="DDEP_NOX"
  chemgroups(6)%specs=>DDEP_NOX_GROUP

  chemgroups(7)%name="PPM_C"
  chemgroups(7)%specs=>PPM_C_GROUP

  chemgroups(8)%name="OM25"
  chemgroups(8)%specs=>OM25_GROUP

  chemgroups(9)%name="DDEP_DUST_NAT_C"
  chemgroups(9)%specs=>DDEP_DUST_NAT_C_GROUP

  chemgroups(10)%name="PPM25_FIRE"
  chemgroups(10)%specs=>PPM25_FIRE_GROUP

  chemgroups(11)%name="WDEP_SVWOODOA25"
  chemgroups(11)%specs=>WDEP_SVWOODOA25_GROUP

  chemgroups(12)%name="DDEP_PMFINE"
  chemgroups(12)%specs=>DDEP_PMFINE_GROUP

  chemgroups(13)%name="NMVOC"
  chemgroups(13)%specs=>NMVOC_GROUP

  chemgroups(14)%name="WDEP_PPM25"
  chemgroups(14)%specs=>WDEP_PPM25_GROUP

  chemgroups(15)%name="DDEP_RCHO"
  chemgroups(15)%specs=>DDEP_RCHO_GROUP

  chemgroups(16)%name="WDEP_PM10"
  chemgroups(16)%specs=>WDEP_PM10_GROUP

  chemgroups(17)%name="DDEP_OX"
  chemgroups(17)%specs=>DDEP_OX_GROUP

  chemgroups(18)%name="DDEP_ECCOARSE"
  chemgroups(18)%specs=>DDEP_ECCOARSE_GROUP

  chemgroups(19)%name="NVWOODOC25"
  chemgroups(19)%specs=>NVWOODOC25_GROUP

  chemgroups(20)%name="WDEP_TMPOX"
  chemgroups(20)%specs=>WDEP_TMPOX_GROUP

  chemgroups(21)%name="ASH"
  chemgroups(21)%specs=>ASH_GROUP

  chemgroups(22)%name="WDEP_NVFFUELOC_COARSE"
  chemgroups(22)%specs=>WDEP_NVFFUELOC_COARSE_GROUP

  chemgroups(23)%name="PM10"
  chemgroups(23)%specs=>PM10_GROUP

  chemgroups(24)%name="DDEP_OM25"
  chemgroups(24)%specs=>DDEP_OM25_GROUP

  chemgroups(25)%name="WDEP_PWOODOA25"
  chemgroups(25)%specs=>WDEP_PWOODOA25_GROUP

  chemgroups(26)%name="DDEP_ASOA"
  chemgroups(26)%specs=>DDEP_ASOA_GROUP

  chemgroups(27)%name="TMPOX"
  chemgroups(27)%specs=>TMPOX_GROUP

  chemgroups(28)%name="OX"
  chemgroups(28)%specs=>OX_GROUP

  chemgroups(29)%name="DDEP_OXN"
  chemgroups(29)%specs=>DDEP_OXN_GROUP

  chemgroups(30)%name="WDEP_PFFUELOA25"
  chemgroups(30)%specs=>WDEP_PFFUELOA25_GROUP

  chemgroups(31)%name="DDEP_PPM10"
  chemgroups(31)%specs=>DDEP_PPM10_GROUP

  chemgroups(32)%name="WDEP_PPM_C"
  chemgroups(32)%specs=>WDEP_PPM_C_GROUP

  chemgroups(33)%name="DDEP_PPM25_FIRE"
  chemgroups(33)%specs=>DDEP_PPM25_FIRE_GROUP

  chemgroups(34)%name="WDEP_EC_F"
  chemgroups(34)%specs=>WDEP_EC_F_GROUP

  chemgroups(35)%name="DDEP_PM10"
  chemgroups(35)%specs=>DDEP_PM10_GROUP

  chemgroups(36)%name="WDEP_SIA"
  chemgroups(36)%specs=>WDEP_SIA_GROUP

  chemgroups(37)%name="BVOC"
  chemgroups(37)%specs=>BVOC_GROUP

  chemgroups(38)%name="PWOODOA25"
  chemgroups(38)%specs=>PWOODOA25_GROUP

  chemgroups(39)%name="EC_F"
  chemgroups(39)%specs=>EC_F_GROUP

  chemgroups(40)%name="DDEP_NONVOLPCM"
  chemgroups(40)%specs=>DDEP_NONVOLPCM_GROUP

  chemgroups(41)%name="RCHO"
  chemgroups(41)%specs=>RCHO_GROUP

  chemgroups(42)%name="WDEP_NVFFIREOC25"
  chemgroups(42)%specs=>WDEP_NVFFIREOC25_GROUP

  chemgroups(43)%name="SOX"
  chemgroups(43)%specs=>SOX_GROUP

  chemgroups(44)%name="DUST_ANT_F"
  chemgroups(44)%specs=>DUST_ANT_F_GROUP

  chemgroups(45)%name="PMFINE"
  chemgroups(45)%specs=>PMFINE_GROUP

  chemgroups(46)%name="DDEP_NVABSOM"
  chemgroups(46)%specs=>DDEP_NVABSOM_GROUP

  chemgroups(47)%name="WDEP_SVFFIREOA25"
  chemgroups(47)%specs=>WDEP_SVFFIREOA25_GROUP

  chemgroups(48)%name="NONVOLPCM"
  chemgroups(48)%specs=>NONVOLPCM_GROUP

  chemgroups(49)%name="WDEP_DAOBS"
  chemgroups(49)%specs=>WDEP_DAOBS_GROUP

  chemgroups(50)%name="DDEP_TMPOX"
  chemgroups(50)%specs=>DDEP_TMPOX_GROUP

  chemgroups(51)%name="WDEP_DUST_NAT_C"
  chemgroups(51)%specs=>WDEP_DUST_NAT_C_GROUP

  chemgroups(52)%name="PMCO"
  chemgroups(52)%specs=>PMCO_GROUP

  chemgroups(53)%name="DDEP_DUST_ANT_F"
  chemgroups(53)%specs=>DDEP_DUST_ANT_F_GROUP

  chemgroups(54)%name="DDEP_OMCOARSE"
  chemgroups(54)%specs=>DDEP_OMCOARSE_GROUP

  chemgroups(55)%name="NVFFIREOC25"
  chemgroups(55)%specs=>NVFFIREOC25_GROUP

  chemgroups(56)%name="WDEP_RDN"
  chemgroups(56)%specs=>WDEP_RDN_GROUP

  chemgroups(57)%name="WDEP_ASOA"
  chemgroups(57)%specs=>WDEP_ASOA_GROUP

  chemgroups(58)%name="WDEP_NVWOODOC25"
  chemgroups(58)%specs=>WDEP_NVWOODOC25_GROUP

  chemgroups(59)%name="DDEP_NVWOODOC25"
  chemgroups(59)%specs=>DDEP_NVWOODOC25_GROUP

  chemgroups(60)%name="WDEP_ECFINE"
  chemgroups(60)%specs=>WDEP_ECFINE_GROUP

  chemgroups(61)%name="WDEP_ECCOARSE"
  chemgroups(61)%specs=>WDEP_ECCOARSE_GROUP

  chemgroups(62)%name="OXN"
  chemgroups(62)%specs=>OXN_GROUP

  chemgroups(63)%name="WDEP_DUST_ANT_F"
  chemgroups(63)%specs=>WDEP_DUST_ANT_F_GROUP

  chemgroups(64)%name="DDEP_NVFFUELOC25"
  chemgroups(64)%specs=>DDEP_NVFFUELOC25_GROUP

  chemgroups(65)%name="SIA"
  chemgroups(65)%specs=>SIA_GROUP

  chemgroups(66)%name="DDEP_ASH"
  chemgroups(66)%specs=>DDEP_ASH_GROUP

  chemgroups(67)%name="DDEP_NVFFIREOC25"
  chemgroups(67)%specs=>DDEP_NVFFIREOC25_GROUP

  chemgroups(68)%name="DDEP_TNO3"
  chemgroups(68)%specs=>DDEP_TNO3_GROUP

  chemgroups(69)%name="DDEP_PMCO"
  chemgroups(69)%specs=>DDEP_PMCO_GROUP

  chemgroups(70)%name="DDEP_BSOA"
  chemgroups(70)%specs=>DDEP_BSOA_GROUP

  chemgroups(71)%name="DDEP_PCM"
  chemgroups(71)%specs=>DDEP_PCM_GROUP

  chemgroups(72)%name="DDEP_NVFFUELOC_COARSE"
  chemgroups(72)%specs=>DDEP_NVFFUELOC_COARSE_GROUP

  chemgroups(73)%name="ECCOARSE"
  chemgroups(73)%specs=>ECCOARSE_GROUP

  chemgroups(74)%name="WOODEC"
  chemgroups(74)%specs=>WOODEC_GROUP

  chemgroups(75)%name="WDEP_TNO3"
  chemgroups(75)%specs=>WDEP_TNO3_GROUP

  chemgroups(76)%name="DDEP_DUST_ANT_C"
  chemgroups(76)%specs=>DDEP_DUST_ANT_C_GROUP

  chemgroups(77)%name="WDEP_PMCO"
  chemgroups(77)%specs=>WDEP_PMCO_GROUP

  chemgroups(78)%name="DDEP_FFUELEC"
  chemgroups(78)%specs=>DDEP_FFUELEC_GROUP

  chemgroups(79)%name="WDEP_NVFFUELOC25"
  chemgroups(79)%specs=>WDEP_NVFFUELOC25_GROUP

  chemgroups(80)%name="WDEP_PM10ANTHR"
  chemgroups(80)%specs=>WDEP_PM10ANTHR_GROUP

  chemgroups(81)%name="WDEP_SVFFUELOA25"
  chemgroups(81)%specs=>WDEP_SVFFUELOA25_GROUP

  chemgroups(82)%name="DDEP_PPM25"
  chemgroups(82)%specs=>DDEP_PPM25_GROUP

  chemgroups(83)%name="WDEP_PPM25_FIRE"
  chemgroups(83)%specs=>WDEP_PPM25_FIRE_GROUP

  chemgroups(84)%name="FFIREBC"
  chemgroups(84)%specs=>FFIREBC_GROUP

  chemgroups(85)%name="WDEP_FFIREBC"
  chemgroups(85)%specs=>WDEP_FFIREBC_GROUP

  chemgroups(86)%name="NOX"
  chemgroups(86)%specs=>NOX_GROUP

  chemgroups(87)%name="DUST_NAT_F"
  chemgroups(87)%specs=>DUST_NAT_F_GROUP

  chemgroups(88)%name="SS"
  chemgroups(88)%specs=>SS_GROUP

  chemgroups(89)%name="DDEP_DUST"
  chemgroups(89)%specs=>DDEP_DUST_GROUP

  chemgroups(90)%name="RO2"
  chemgroups(90)%specs=>RO2_GROUP

  chemgroups(91)%name="DDEP_EC_F"
  chemgroups(91)%specs=>DDEP_EC_F_GROUP

  chemgroups(92)%name="ROOH"
  chemgroups(92)%specs=>ROOH_GROUP

  chemgroups(93)%name="DUST_ANT_C"
  chemgroups(93)%specs=>DUST_ANT_C_GROUP

  chemgroups(94)%name="AOD"
  chemgroups(94)%specs=>AOD_GROUP

  chemgroups(95)%name="DDEP_PPM_C"
  chemgroups(95)%specs=>DDEP_PPM_C_GROUP

  chemgroups(96)%name="WDEP_NVABSOM"
  chemgroups(96)%specs=>WDEP_NVABSOM_GROUP

  chemgroups(97)%name="SVFFUELOA25"
  chemgroups(97)%specs=>SVFFUELOA25_GROUP

  chemgroups(98)%name="DDEP_SOX"
  chemgroups(98)%specs=>DDEP_SOX_GROUP

  chemgroups(99)%name="WDEP_DUST_NAT_F"
  chemgroups(99)%specs=>WDEP_DUST_NAT_F_GROUP

  chemgroups(100)%name="WDEP_SOX"
  chemgroups(100)%specs=>WDEP_SOX_GROUP

  chemgroups(101)%name="PFFUELOA25"
  chemgroups(101)%specs=>PFFUELOA25_GROUP

  chemgroups(102)%name="FFUELEC"
  chemgroups(102)%specs=>FFUELEC_GROUP

  chemgroups(103)%name="PM10ANTHR"
  chemgroups(103)%specs=>PM10ANTHR_GROUP

  chemgroups(104)%name="WDEP_DUST_ANT_C"
  chemgroups(104)%specs=>WDEP_DUST_ANT_C_GROUP

  chemgroups(105)%name="DDEP_DUST_NAT_F"
  chemgroups(105)%specs=>DDEP_DUST_NAT_F_GROUP

  chemgroups(106)%name="DDEP_PM10ANTHR"
  chemgroups(106)%specs=>DDEP_PM10ANTHR_GROUP

  chemgroups(107)%name="DDEP_WOODEC"
  chemgroups(107)%specs=>DDEP_WOODEC_GROUP

  chemgroups(108)%name="WDEP_WOODEC"
  chemgroups(108)%specs=>WDEP_WOODEC_GROUP

  chemgroups(109)%name="WDEP_PCM"
  chemgroups(109)%specs=>WDEP_PCM_GROUP

  chemgroups(110)%name="WDEP_SS"
  chemgroups(110)%specs=>WDEP_SS_GROUP

  chemgroups(111)%name="DUST_NAT_C"
  chemgroups(111)%specs=>DUST_NAT_C_GROUP

  chemgroups(112)%name="MONOTERP"
  chemgroups(112)%specs=>MONOTERP_GROUP

  chemgroups(113)%name="WDEP_NONVOLPCM"
  chemgroups(113)%specs=>WDEP_NONVOLPCM_GROUP

  chemgroups(114)%name="PPM25"
  chemgroups(114)%specs=>PPM25_GROUP

  chemgroups(115)%name="ASOA"
  chemgroups(115)%specs=>ASOA_GROUP

  chemgroups(116)%name="PPM10"
  chemgroups(116)%specs=>PPM10_GROUP

  chemgroups(117)%name="NVABSOM"
  chemgroups(117)%specs=>NVABSOM_GROUP

  chemgroups(118)%name="BSOA"
  chemgroups(118)%specs=>BSOA_GROUP

  chemgroups(119)%name="CHET2"
  chemgroups(119)%specs=>CHET2_GROUP

  chemgroups(120)%name="ECFINE"
  chemgroups(120)%specs=>ECFINE_GROUP

  chemgroups(121)%name="DDEP_FFIREBC"
  chemgroups(121)%specs=>DDEP_FFIREBC_GROUP

  chemgroups(122)%name="WDEP_DUST"
  chemgroups(122)%specs=>WDEP_DUST_GROUP

  chemgroups(123)%name="GRPS"
  chemgroups(123)%specs=>GRPS_GROUP

  chemgroups(124)%name="NVFFUELOC_COARSE"
  chemgroups(124)%specs=>NVFFUELOC_COARSE_GROUP

  chemgroups(125)%name="DDEP_ECFINE"
  chemgroups(125)%specs=>DDEP_ECFINE_GROUP

  chemgroups(126)%name="WDEP_ROOH"
  chemgroups(126)%specs=>WDEP_ROOH_GROUP

  chemgroups(127)%name="PCM"
  chemgroups(127)%specs=>PCM_GROUP

  chemgroups(128)%name="SVWOODOA25"
  chemgroups(128)%specs=>SVWOODOA25_GROUP

  chemgroups(129)%name="DDEP_SIA"
  chemgroups(129)%specs=>DDEP_SIA_GROUP

  chemgroups(130)%name="WDEP_PFFIREOA25"
  chemgroups(130)%specs=>WDEP_PFFIREOA25_GROUP

  chemgroups(131)%name="WDEP_OMCOARSE"
  chemgroups(131)%specs=>WDEP_OMCOARSE_GROUP

  chemgroups(132)%name="TNO3"
  chemgroups(132)%specs=>TNO3_GROUP

  chemgroups(133)%name="SVFFIREOA25"
  chemgroups(133)%specs=>SVFFIREOA25_GROUP

  chemgroups(134)%name="DDEP_DAOBS"
  chemgroups(134)%specs=>DDEP_DAOBS_GROUP

  chemgroups(135)%name="NVFFUELOC25"
  chemgroups(135)%specs=>NVFFUELOC25_GROUP

  chemgroups(136)%name="OMCOARSE"
  chemgroups(136)%specs=>OMCOARSE_GROUP

  chemgroups(137)%name="DDEP_ROOH"
  chemgroups(137)%specs=>DDEP_ROOH_GROUP

  chemgroups(138)%name="WDEP_PMFINE"
  chemgroups(138)%specs=>WDEP_PMFINE_GROUP

  chemgroups(139)%name="DDEP_RDN"
  chemgroups(139)%specs=>DDEP_RDN_GROUP

  chemgroups(140)%name="PFFIREOA25"
  chemgroups(140)%specs=>PFFIREOA25_GROUP

  chemgroups(141)%name="RDN"
  chemgroups(141)%specs=>RDN_GROUP

  chemgroups(142)%name="WDEP_ASH"
  chemgroups(142)%specs=>WDEP_ASH_GROUP

  chemgroups(143)%name="DAOBS"
  chemgroups(143)%specs=>DAOBS_GROUP

  chemgroups(144)%name="WDEP_FFUELEC"
  chemgroups(144)%specs=>WDEP_FFUELEC_GROUP

  chemgroups(145)%name="WDEP_OM25"
  chemgroups(145)%specs=>WDEP_OM25_GROUP

endsubroutine Init_ChemGroups
 !-----------------------------------------------------------
endmodule ChemGroups_ml
 !-----------------------------------------------------------
