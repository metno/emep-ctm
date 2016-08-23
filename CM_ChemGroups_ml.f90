!>_________________________________________________________<

  module  ChemGroups_ml
!-----------------------------------------------------------

  use ChemSpecs_tot_ml  ! => species indices
use OwnDataTypes_ml   ! => typ_sp
  implicit none
  private
! Assignment of groups from GenIn.species:
 public :: Init_ChemGroups

! ------- Gas/particle species ------------------

  integer, public, parameter ::  INDEX_DDEP_SS_GROUP = 1
  integer, public, target, save, dimension(2) :: &
                 DDEP_SS_GROUP     = (/ SEASALT_F,SEASALT_C /)

  integer, public, parameter ::  INDEX_WDEP_OXN_GROUP = 2
  integer, public, target, save, dimension(4) :: &
                 WDEP_OXN_GROUP     = (/ HNO3,HONO,NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_WDEP_PPM10_GROUP = 3
  integer, public, target, save, dimension(11) :: &
                 WDEP_PPM10_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C /)

  integer, public, parameter ::  INDEX_ASH_C_GROUP = 4
  integer, public, target, save, dimension(1) :: &
                 ASH_C_GROUP     = (/ V1702A02B_C /)

  integer, public, parameter ::  INDEX_DUST_GROUP = 5
  integer, public, target, save, dimension(6) :: &
                 DUST_GROUP     = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_WDEP_AOD_GROUP = 6
  integer, public, target, save, dimension(13) :: &
                 WDEP_AOD_GROUP     = (/ SO4,NO3_F,NH4_F,V1702A02B_F,V1702A02B_C,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_WDEP_BSOA_GROUP = 7
  integer, public, target, save, dimension(10) :: &
                 WDEP_BSOA_GROUP     = (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

  integer, public, parameter ::  INDEX_DDEP_NOX_GROUP = 8
  integer, public, target, save, dimension(1) :: &
                 DDEP_NOX_GROUP     = (/ NO2 /)

  integer, public, parameter ::  INDEX_PPM_C_GROUP = 9
  integer, public, target, save, dimension(4) :: &
                 PPM_C_GROUP     = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C /)

  integer, public, parameter ::  INDEX_DDEP_DUST_NAT_C_GROUP = 10
  integer, public, target, save, dimension(2) :: &
                 DDEP_DUST_NAT_C_GROUP     = (/ DUST_WB_C,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_PPM25_FIRE_GROUP = 11
  integer, public, target, save, dimension(3) :: &
                 PPM25_FIRE_GROUP     = (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25 /)

  integer, public, parameter ::  INDEX_WDEP_SVWOODOA25_GROUP = 12
  integer, public, target, save, dimension(1) :: &
                 WDEP_SVWOODOA25_GROUP     = (/ WOODOA_NG10 /)

  integer, public, parameter ::  INDEX_DDEP_PMFINE_GROUP = 13
  integer, public, target, save, dimension(15) :: &
                 DDEP_PMFINE_GROUP     = (/ SO4,NO3_F,NH4_F,V1702A02B_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

  integer, public, parameter ::  INDEX_WDEP_PPM25_GROUP = 14
  integer, public, target, save, dimension(7) :: &
                 WDEP_PPM25_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

  integer, public, parameter ::  INDEX_WDEP_PM10_GROUP = 15
  integer, public, target, save, dimension(25) :: &
                 WDEP_PM10_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F,V1702A02B_F,V1702A02B_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_DDEP_OX_GROUP = 16
  integer, public, target, save, dimension(2) :: &
                 DDEP_OX_GROUP     = (/ O3,NO2 /)

  integer, public, parameter ::  INDEX_DDEP_ECCOARSE_GROUP = 17
  integer, public, target, save, dimension(2) :: &
                 DDEP_ECCOARSE_GROUP     = (/ EC_C_WOOD,EC_C_FFUEL /)

  integer, public, parameter ::  INDEX_NVWOODOC25_GROUP = 18
  integer, public, target, save, dimension(1) :: &
                 NVWOODOC25_GROUP     = (/ POM_F_WOOD /)

  integer, public, parameter ::  INDEX_ASH_GROUP = 19
  integer, public, target, save, dimension(2) :: &
                 ASH_GROUP     = (/ V1702A02B_F,V1702A02B_C /)

  integer, public, parameter ::  INDEX_WDEP_NVFFUELOC_COARSE_GROUP = 20
  integer, public, target, save, dimension(1) :: &
                 WDEP_NVFFUELOC_COARSE_GROUP     = (/ POM_C_FFUEL /)

  integer, public, parameter ::  INDEX_PM10_GROUP = 21
  integer, public, target, save, dimension(26) :: &
                 PM10_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F,V1702A02B_F,V1702A02B_C,PART_OM_F,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_WDEP_PWOODOA25_GROUP = 22
  integer, public, target, save, dimension(1) :: &
                 WDEP_PWOODOA25_GROUP     = (/ WOODOA_NG10 /)

  integer, public, parameter ::  INDEX_DDEP_ASOA_GROUP = 23
  integer, public, target, save, dimension(10) :: &
                 DDEP_ASOA_GROUP     = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

  integer, public, parameter ::  INDEX_OX_GROUP = 24
  integer, public, target, save, dimension(2) :: &
                 OX_GROUP     = (/ O3,NO2 /)

  integer, public, parameter ::  INDEX_DDEP_OXN_GROUP = 25
  integer, public, target, save, dimension(7) :: &
                 DDEP_OXN_GROUP     = (/ NO2,PAN,MPAN,HNO3,HONO,NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_WDEP_PFFUELOA25_GROUP = 26
  integer, public, target, save, dimension(1) :: &
                 WDEP_PFFUELOA25_GROUP     = (/ FFFUEL_NG10 /)

  integer, public, parameter ::  INDEX_DDEP_PPM10_GROUP = 27
  integer, public, target, save, dimension(11) :: &
                 DDEP_PPM10_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C /)

  integer, public, parameter ::  INDEX_WDEP_PPM_C_GROUP = 28
  integer, public, target, save, dimension(4) :: &
                 WDEP_PPM_C_GROUP     = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C /)

  integer, public, parameter ::  INDEX_DDEP_PPM25_FIRE_GROUP = 29
  integer, public, target, save, dimension(3) :: &
                 DDEP_PPM25_FIRE_GROUP     = (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25 /)

  integer, public, parameter ::  INDEX_WDEP_EC_F_GROUP = 30
  integer, public, target, save, dimension(5) :: &
                 WDEP_EC_F_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC /)

  integer, public, parameter ::  INDEX_DDEP_PM10_GROUP = 31
  integer, public, target, save, dimension(25) :: &
                 DDEP_PM10_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F,V1702A02B_F,V1702A02B_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_WDEP_SIA_GROUP = 32
  integer, public, target, save, dimension(4) :: &
                 WDEP_SIA_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F /)

  integer, public, parameter ::  INDEX_BVOC_GROUP = 33
  integer, public, target, save, dimension(2) :: &
                 BVOC_GROUP     = (/ C5H8,APINENE /)

  integer, public, parameter ::  INDEX_PWOODOA25_GROUP = 34
  integer, public, target, save, dimension(1) :: &
                 PWOODOA25_GROUP     = (/ WOODOA_NG10 /)

  integer, public, parameter ::  INDEX_EC_F_GROUP = 35
  integer, public, target, save, dimension(5) :: &
                 EC_F_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC /)

  integer, public, parameter ::  INDEX_DDEP_NONVOLPCM_GROUP = 36
  integer, public, target, save, dimension(10) :: &
                 DDEP_NONVOLPCM_GROUP     = (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC /)

  integer, public, parameter ::  INDEX_WDEP_NVFFIREOC25_GROUP = 37
  integer, public, target, save, dimension(1) :: &
                 WDEP_NVFFIREOC25_GROUP     = (/ FFIRE_OM /)

  integer, public, parameter ::  INDEX_SOX_GROUP = 38
  integer, public, target, save, dimension(2) :: &
                 SOX_GROUP     = (/ SO2,SO4 /)

  integer, public, parameter ::  INDEX_DUST_ANT_F_GROUP = 39
  integer, public, target, save, dimension(1) :: &
                 DUST_ANT_F_GROUP     = (/ DUST_ROAD_F /)

  integer, public, parameter ::  INDEX_PMFINE_GROUP = 40
  integer, public, target, save, dimension(16) :: &
                 PMFINE_GROUP     = (/ SO4,NO3_F,NH4_F,V1702A02B_F,PART_OM_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

  integer, public, parameter ::  INDEX_DDEP_NVABSOM_GROUP = 41
  integer, public, target, save, dimension(3) :: &
                 DDEP_NVABSOM_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM /)

  integer, public, parameter ::  INDEX_WDEP_SVFFIREOA25_GROUP = 42
  integer, public, target, save, dimension(1) :: &
                 WDEP_SVFFIREOA25_GROUP     = (/ FFIREOA_NG10 /)

  integer, public, parameter ::  INDEX_DDEP_AOD_GROUP = 43
  integer, public, target, save, dimension(13) :: &
                 DDEP_AOD_GROUP     = (/ SO4,NO3_F,NH4_F,V1702A02B_F,V1702A02B_C,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_NONVOLPCM_GROUP = 44
  integer, public, target, save, dimension(10) :: &
                 NONVOLPCM_GROUP     = (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC /)

  integer, public, parameter ::  INDEX_WDEP_V1702A02B_GROUP = 45
  integer, public, target, save, dimension(2) :: &
                 WDEP_V1702A02B_GROUP     = (/ V1702A02B_F,V1702A02B_C /)

  integer, public, parameter ::  INDEX_WDEP_ASH_F_GROUP = 46
  integer, public, target, save, dimension(1) :: &
                 WDEP_ASH_F_GROUP     = (/ V1702A02B_F /)

  integer, public, parameter ::  INDEX_WDEP_DUST_NAT_C_GROUP = 47
  integer, public, target, save, dimension(2) :: &
                 WDEP_DUST_NAT_C_GROUP     = (/ DUST_WB_C,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_PMCO_GROUP = 48
  integer, public, target, save, dimension(10) :: &
                 PMCO_GROUP     = (/ NO3_C,V1702A02B_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_DDEP_DUST_ANT_F_GROUP = 49
  integer, public, target, save, dimension(1) :: &
                 DDEP_DUST_ANT_F_GROUP     = (/ DUST_ROAD_F /)

  integer, public, parameter ::  INDEX_DDEP_OMCOARSE_GROUP = 50
  integer, public, target, save, dimension(1) :: &
                 DDEP_OMCOARSE_GROUP     = (/ POM_C_FFUEL /)

  integer, public, parameter ::  INDEX_NVFFIREOC25_GROUP = 51
  integer, public, target, save, dimension(1) :: &
                 NVFFIREOC25_GROUP     = (/ FFIRE_OM /)

  integer, public, parameter ::  INDEX_WDEP_RDN_GROUP = 52
  integer, public, target, save, dimension(2) :: &
                 WDEP_RDN_GROUP     = (/ NH3,NH4_F /)

  integer, public, parameter ::  INDEX_WDEP_ASOA_GROUP = 53
  integer, public, target, save, dimension(10) :: &
                 WDEP_ASOA_GROUP     = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

  integer, public, parameter ::  INDEX_WDEP_NVWOODOC25_GROUP = 54
  integer, public, target, save, dimension(1) :: &
                 WDEP_NVWOODOC25_GROUP     = (/ POM_F_WOOD /)

  integer, public, parameter ::  INDEX_WDEP_ASH_C_GROUP = 55
  integer, public, target, save, dimension(1) :: &
                 WDEP_ASH_C_GROUP     = (/ V1702A02B_C /)

  integer, public, parameter ::  INDEX_DDEP_NVWOODOC25_GROUP = 56
  integer, public, target, save, dimension(1) :: &
                 DDEP_NVWOODOC25_GROUP     = (/ POM_F_WOOD /)

  integer, public, parameter ::  INDEX_WDEP_ECFINE_GROUP = 57
  integer, public, target, save, dimension(4) :: &
                 WDEP_ECFINE_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE /)

  integer, public, parameter ::  INDEX_WDEP_ECCOARSE_GROUP = 58
  integer, public, target, save, dimension(2) :: &
                 WDEP_ECCOARSE_GROUP     = (/ EC_C_WOOD,EC_C_FFUEL /)

  integer, public, parameter ::  INDEX_OXN_GROUP = 59
  integer, public, target, save, dimension(13) :: &
                 OXN_GROUP     = (/ NO,NO2,PAN,MPAN,NO3,N2O5,ISONO3,HNO3,HONO,ISNI,ISNIR,NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_WDEP_DUST_ANT_F_GROUP = 60
  integer, public, target, save, dimension(1) :: &
                 WDEP_DUST_ANT_F_GROUP     = (/ DUST_ROAD_F /)

  integer, public, parameter ::  INDEX_DDEP_NVFFUELOC25_GROUP = 61
  integer, public, target, save, dimension(1) :: &
                 DDEP_NVFFUELOC25_GROUP     = (/ POM_F_FFUEL /)

  integer, public, parameter ::  INDEX_SIA_GROUP = 62
  integer, public, target, save, dimension(4) :: &
                 SIA_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F /)

  integer, public, parameter ::  INDEX_DDEP_ASH_GROUP = 63
  integer, public, target, save, dimension(2) :: &
                 DDEP_ASH_GROUP     = (/ V1702A02B_F,V1702A02B_C /)

  integer, public, parameter ::  INDEX_DDEP_NVFFIREOC25_GROUP = 64
  integer, public, target, save, dimension(1) :: &
                 DDEP_NVFFIREOC25_GROUP     = (/ FFIRE_OM /)

  integer, public, parameter ::  INDEX_DDEP_TNO3_GROUP = 65
  integer, public, target, save, dimension(2) :: &
                 DDEP_TNO3_GROUP     = (/ NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_DDEP_ASH_F_GROUP = 66
  integer, public, target, save, dimension(1) :: &
                 DDEP_ASH_F_GROUP     = (/ V1702A02B_F /)

  integer, public, parameter ::  INDEX_DDEP_PMCO_GROUP = 67
  integer, public, target, save, dimension(10) :: &
                 DDEP_PMCO_GROUP     = (/ NO3_C,V1702A02B_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_DDEP_BSOA_GROUP = 68
  integer, public, target, save, dimension(10) :: &
                 DDEP_BSOA_GROUP     = (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

  integer, public, parameter ::  INDEX_DDEP_PCM_GROUP = 69
  integer, public, target, save, dimension(28) :: &
                 DDEP_PCM_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

  integer, public, parameter ::  INDEX_DDEP_NVFFUELOC_COARSE_GROUP = 70
  integer, public, target, save, dimension(1) :: &
                 DDEP_NVFFUELOC_COARSE_GROUP     = (/ POM_C_FFUEL /)

  integer, public, parameter ::  INDEX_ECCOARSE_GROUP = 71
  integer, public, target, save, dimension(2) :: &
                 ECCOARSE_GROUP     = (/ EC_C_WOOD,EC_C_FFUEL /)

  integer, public, parameter ::  INDEX_WOODEC_GROUP = 72
  integer, public, target, save, dimension(3) :: &
                 WOODEC_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD /)

  integer, public, parameter ::  INDEX_WDEP_TNO3_GROUP = 73
  integer, public, target, save, dimension(2) :: &
                 WDEP_TNO3_GROUP     = (/ NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_DDEP_DUST_ANT_C_GROUP = 74
  integer, public, target, save, dimension(1) :: &
                 DDEP_DUST_ANT_C_GROUP     = (/ DUST_ROAD_C /)

  integer, public, parameter ::  INDEX_WDEP_PMCO_GROUP = 75
  integer, public, target, save, dimension(10) :: &
                 WDEP_PMCO_GROUP     = (/ NO3_C,V1702A02B_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_DDEP_FFUELEC_GROUP = 76
  integer, public, target, save, dimension(3) :: &
                 DDEP_FFUELEC_GROUP     = (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL /)

  integer, public, parameter ::  INDEX_DDEP_ASH_C_GROUP = 77
  integer, public, target, save, dimension(1) :: &
                 DDEP_ASH_C_GROUP     = (/ V1702A02B_C /)

  integer, public, parameter ::  INDEX_WDEP_NVFFUELOC25_GROUP = 78
  integer, public, target, save, dimension(1) :: &
                 WDEP_NVFFUELOC25_GROUP     = (/ POM_F_FFUEL /)

  integer, public, parameter ::  INDEX_WDEP_PM10ANTHR_GROUP = 79
  integer, public, target, save, dimension(3) :: &
                 WDEP_PM10ANTHR_GROUP     = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL /)

  integer, public, parameter ::  INDEX_WDEP_SVFFUELOA25_GROUP = 80
  integer, public, target, save, dimension(1) :: &
                 WDEP_SVFFUELOA25_GROUP     = (/ FFFUEL_NG10 /)

  integer, public, parameter ::  INDEX_DDEP_PPM25_GROUP = 81
  integer, public, target, save, dimension(7) :: &
                 DDEP_PPM25_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

  integer, public, parameter ::  INDEX_WDEP_PPM25_FIRE_GROUP = 82
  integer, public, target, save, dimension(3) :: &
                 WDEP_PPM25_FIRE_GROUP     = (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25 /)

  integer, public, parameter ::  INDEX_FFIREBC_GROUP = 83
  integer, public, target, save, dimension(1) :: &
                 FFIREBC_GROUP     = (/ FFIRE_BC /)

  integer, public, parameter ::  INDEX_WDEP_FFIREBC_GROUP = 84
  integer, public, target, save, dimension(1) :: &
                 WDEP_FFIREBC_GROUP     = (/ FFIRE_BC /)

  integer, public, parameter ::  INDEX_NOX_GROUP = 85
  integer, public, target, save, dimension(2) :: &
                 NOX_GROUP     = (/ NO,NO2 /)

  integer, public, parameter ::  INDEX_DUST_NAT_F_GROUP = 86
  integer, public, target, save, dimension(2) :: &
                 DUST_NAT_F_GROUP     = (/ DUST_WB_F,DUST_SAH_F /)

  integer, public, parameter ::  INDEX_SS_GROUP = 87
  integer, public, target, save, dimension(2) :: &
                 SS_GROUP     = (/ SEASALT_F,SEASALT_C /)

  integer, public, parameter ::  INDEX_DDEP_DUST_GROUP = 88
  integer, public, target, save, dimension(6) :: &
                 DDEP_DUST_GROUP     = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_RO2_GROUP = 89
  integer, public, target, save, dimension(14) :: &
                 RO2_GROUP     = (/ HO2,CH3O2,C2H5O2,SECC4H9O2,ISRO2,ETRO2,PRRO2,OXYO2,MEKO2,MALO2,MVKO2,MACRO2,MACO3,TERPPEROXY /)

  integer, public, parameter ::  INDEX_DDEP_EC_F_GROUP = 90
  integer, public, target, save, dimension(5) :: &
                 DDEP_EC_F_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC /)

  integer, public, parameter ::  INDEX_ROOH_GROUP = 91
  integer, public, target, save, dimension(16) :: &
                 ROOH_GROUP     = (/ CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,OXYO2H,MEKO2H,MALO2H,MVKO2H,MACROOH,MACO3H,ISRO2H,H2O2,CH3COO2H,ISONO3H,ISNIRH /)

  integer, public, parameter ::  INDEX_DUST_ANT_C_GROUP = 92
  integer, public, target, save, dimension(1) :: &
                 DUST_ANT_C_GROUP     = (/ DUST_ROAD_C /)

  integer, public, parameter ::  INDEX_AOD_GROUP = 93
  integer, public, target, save, dimension(13) :: &
                 AOD_GROUP     = (/ SO4,NO3_F,NH4_F,V1702A02B_F,V1702A02B_C,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_DDEP_PPM_C_GROUP = 94
  integer, public, target, save, dimension(4) :: &
                 DDEP_PPM_C_GROUP     = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C /)

  integer, public, parameter ::  INDEX_WDEP_NVABSOM_GROUP = 95
  integer, public, target, save, dimension(3) :: &
                 WDEP_NVABSOM_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM /)

  integer, public, parameter ::  INDEX_SVFFUELOA25_GROUP = 96
  integer, public, target, save, dimension(1) :: &
                 SVFFUELOA25_GROUP     = (/ FFFUEL_NG10 /)

  integer, public, parameter ::  INDEX_DDEP_SOX_GROUP = 97
  integer, public, target, save, dimension(2) :: &
                 DDEP_SOX_GROUP     = (/ SO2,SO4 /)

  integer, public, parameter ::  INDEX_WDEP_DUST_NAT_F_GROUP = 98
  integer, public, target, save, dimension(2) :: &
                 WDEP_DUST_NAT_F_GROUP     = (/ DUST_WB_F,DUST_SAH_F /)

  integer, public, parameter ::  INDEX_ASH_F_GROUP = 99
  integer, public, target, save, dimension(1) :: &
                 ASH_F_GROUP     = (/ V1702A02B_F /)

  integer, public, parameter ::  INDEX_WDEP_SOX_GROUP = 100
  integer, public, target, save, dimension(2) :: &
                 WDEP_SOX_GROUP     = (/ SO2,SO4 /)

  integer, public, parameter ::  INDEX_PFFUELOA25_GROUP = 101
  integer, public, target, save, dimension(1) :: &
                 PFFUELOA25_GROUP     = (/ FFFUEL_NG10 /)

  integer, public, parameter ::  INDEX_FFUELEC_GROUP = 102
  integer, public, target, save, dimension(3) :: &
                 FFUELEC_GROUP     = (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL /)

  integer, public, parameter ::  INDEX_PM10ANTHR_GROUP = 103
  integer, public, target, save, dimension(3) :: &
                 PM10ANTHR_GROUP     = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL /)

  integer, public, parameter ::  INDEX_WDEP_DUST_ANT_C_GROUP = 104
  integer, public, target, save, dimension(1) :: &
                 WDEP_DUST_ANT_C_GROUP     = (/ DUST_ROAD_C /)

  integer, public, parameter ::  INDEX_DDEP_DUST_NAT_F_GROUP = 105
  integer, public, target, save, dimension(2) :: &
                 DDEP_DUST_NAT_F_GROUP     = (/ DUST_WB_F,DUST_SAH_F /)

  integer, public, parameter ::  INDEX_DDEP_PM10ANTHR_GROUP = 106
  integer, public, target, save, dimension(3) :: &
                 DDEP_PM10ANTHR_GROUP     = (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL /)

  integer, public, parameter ::  INDEX_DDEP_WOODEC_GROUP = 107
  integer, public, target, save, dimension(3) :: &
                 DDEP_WOODEC_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD /)

  integer, public, parameter ::  INDEX_V1702A02B_GROUP = 108
  integer, public, target, save, dimension(2) :: &
                 V1702A02B_GROUP     = (/ V1702A02B_F,V1702A02B_C /)

  integer, public, parameter ::  INDEX_WDEP_WOODEC_GROUP = 109
  integer, public, target, save, dimension(3) :: &
                 WDEP_WOODEC_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD /)

  integer, public, parameter ::  INDEX_WDEP_PCM_GROUP = 110
  integer, public, target, save, dimension(31) :: &
                 WDEP_PCM_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

  integer, public, parameter ::  INDEX_WDEP_SS_GROUP = 111
  integer, public, target, save, dimension(2) :: &
                 WDEP_SS_GROUP     = (/ SEASALT_F,SEASALT_C /)

  integer, public, parameter ::  INDEX_DUST_NAT_C_GROUP = 112
  integer, public, target, save, dimension(2) :: &
                 DUST_NAT_C_GROUP     = (/ DUST_WB_C,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_WDEP_NONVOLPCM_GROUP = 113
  integer, public, target, save, dimension(10) :: &
                 WDEP_NONVOLPCM_GROUP     = (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC /)

  integer, public, parameter ::  INDEX_PPM25_GROUP = 114
  integer, public, target, save, dimension(7) :: &
                 PPM25_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25 /)

  integer, public, parameter ::  INDEX_ASOA_GROUP = 115
  integer, public, target, save, dimension(10) :: &
                 ASOA_GROUP     = (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3 /)

  integer, public, parameter ::  INDEX_PPM10_GROUP = 116
  integer, public, target, save, dimension(11) :: &
                 PPM10_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C /)

  integer, public, parameter ::  INDEX_NVABSOM_GROUP = 117
  integer, public, target, save, dimension(3) :: &
                 NVABSOM_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM /)

  integer, public, parameter ::  INDEX_BSOA_GROUP = 118
  integer, public, target, save, dimension(10) :: &
                 BSOA_GROUP     = (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3 /)

  integer, public, parameter ::  INDEX_ECFINE_GROUP = 119
  integer, public, target, save, dimension(4) :: &
                 ECFINE_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE /)

  integer, public, parameter ::  INDEX_DDEP_FFIREBC_GROUP = 120
  integer, public, target, save, dimension(1) :: &
                 DDEP_FFIREBC_GROUP     = (/ FFIRE_BC /)

  integer, public, parameter ::  INDEX_WDEP_DUST_GROUP = 121
  integer, public, target, save, dimension(6) :: &
                 WDEP_DUST_GROUP     = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

  integer, public, parameter ::  INDEX_NVFFUELOC_COARSE_GROUP = 122
  integer, public, target, save, dimension(1) :: &
                 NVFFUELOC_COARSE_GROUP     = (/ POM_C_FFUEL /)

  integer, public, parameter ::  INDEX_DDEP_ECFINE_GROUP = 123
  integer, public, target, save, dimension(4) :: &
                 DDEP_ECFINE_GROUP     = (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE /)

  integer, public, parameter ::  INDEX_WDEP_ROOH_GROUP = 124
  integer, public, target, save, dimension(1) :: &
                 WDEP_ROOH_GROUP     = (/ H2O2 /)

  integer, public, parameter ::  INDEX_PCM_GROUP = 125
  integer, public, target, save, dimension(31) :: &
                 PCM_GROUP     = (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /)

  integer, public, parameter ::  INDEX_SVWOODOA25_GROUP = 126
  integer, public, target, save, dimension(1) :: &
                 SVWOODOA25_GROUP     = (/ WOODOA_NG10 /)

  integer, public, parameter ::  INDEX_DDEP_SIA_GROUP = 127
  integer, public, target, save, dimension(4) :: &
                 DDEP_SIA_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F /)

  integer, public, parameter ::  INDEX_WDEP_PFFIREOA25_GROUP = 128
  integer, public, target, save, dimension(1) :: &
                 WDEP_PFFIREOA25_GROUP     = (/ FFIREOA_NG10 /)

  integer, public, parameter ::  INDEX_PCM_HELP_GROUP = 129
  integer, public, target, save, dimension(20) :: &
                 PCM_HELP_GROUP     = (/ GAS_ASOA_OC,PART_ASOA_OC,PART_ASOA_OM,GAS_BSOA_OC,PART_BSOA_OC,PART_BSOA_OM,PART_FFUELOA25_OC,PART_FFUELOA25_OM,PART_WOODOA25_OC,PART_WOODOA25_OM,PART_FFIREOA25_OC,PART_FFIREOA25_OM,PART_OC10,PART_OC25,NONVOL_FFUELOC25,NONV_FFUELOC_COARSE,NONVOL_WOODOC25,NONVOL_BGNDOC,NONVOL_FFIREOC25,PART_OM_F /)

  integer, public, parameter ::  INDEX_WDEP_OMCOARSE_GROUP = 130
  integer, public, target, save, dimension(1) :: &
                 WDEP_OMCOARSE_GROUP     = (/ POM_C_FFUEL /)

  integer, public, parameter ::  INDEX_TNO3_GROUP = 131
  integer, public, target, save, dimension(2) :: &
                 TNO3_GROUP     = (/ NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_SVFFIREOA25_GROUP = 132
  integer, public, target, save, dimension(1) :: &
                 SVFFIREOA25_GROUP     = (/ FFIREOA_NG10 /)

  integer, public, parameter ::  INDEX_NVFFUELOC25_GROUP = 133
  integer, public, target, save, dimension(1) :: &
                 NVFFUELOC25_GROUP     = (/ POM_F_FFUEL /)

  integer, public, parameter ::  INDEX_DDEP_V1702A02B_GROUP = 134
  integer, public, target, save, dimension(2) :: &
                 DDEP_V1702A02B_GROUP     = (/ V1702A02B_F,V1702A02B_C /)

  integer, public, parameter ::  INDEX_OMCOARSE_GROUP = 135
  integer, public, target, save, dimension(1) :: &
                 OMCOARSE_GROUP     = (/ POM_C_FFUEL /)

  integer, public, parameter ::  INDEX_DDEP_ROOH_GROUP = 136
  integer, public, target, save, dimension(3) :: &
                 DDEP_ROOH_GROUP     = (/ CH3O2H,C2H5OOH,H2O2 /)

  integer, public, parameter ::  INDEX_WDEP_PMFINE_GROUP = 137
  integer, public, target, save, dimension(15) :: &
                 WDEP_PMFINE_GROUP     = (/ SO4,NO3_F,NH4_F,V1702A02B_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

  integer, public, parameter ::  INDEX_DDEP_RDN_GROUP = 138
  integer, public, target, save, dimension(2) :: &
                 DDEP_RDN_GROUP     = (/ NH3,NH4_F /)

  integer, public, parameter ::  INDEX_PFFIREOA25_GROUP = 139
  integer, public, target, save, dimension(1) :: &
                 PFFIREOA25_GROUP     = (/ FFIREOA_NG10 /)

  integer, public, parameter ::  INDEX_RDN_GROUP = 140
  integer, public, target, save, dimension(2) :: &
                 RDN_GROUP     = (/ NH3,NH4_F /)

  integer, public, parameter ::  INDEX_WDEP_ASH_GROUP = 141
  integer, public, target, save, dimension(2) :: &
                 WDEP_ASH_GROUP     = (/ V1702A02B_F,V1702A02B_C /)

  integer, public, parameter ::  INDEX_WDEP_FFUELEC_GROUP = 142
  integer, public, target, save, dimension(3) :: &
                 WDEP_FFUELEC_GROUP     = (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL /)


!GROUP ARRAY SIZE 142 MAXN 31 

  type, public :: gtype
       character(len=20) :: name
       integer :: Ngroup
       integer, dimension(31) :: itot   ! indices from xn_tot arrays
  end type gtype

  type(gtype), public, parameter, dimension(142) :: &
       GROUP_ARRAY = (/ &
 gtype( "DDEP_SS", 2, (/ SEASALT_F,SEASALT_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_OXN", 4, (/ HNO3,HONO,NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PPM10", 11, (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "ASH_C", 1, (/ V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DUST", 6, (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_AOD", 13, (/ SO4,NO3_F,NH4_F,V1702A02B_F,V1702A02B_C,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_BSOA", 10, (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_NOX", 1, (/ NO2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PPM_C", 4, (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_DUST_NAT_C", 2, (/ DUST_WB_C,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PPM25_FIRE", 3, (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_SVWOODOA25", 1, (/ WOODOA_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PMFINE", 15, (/ SO4,NO3_F,NH4_F,V1702A02B_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PPM25", 7, (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PM10", 25, (/ SO4,NO3_F,NO3_C,NH4_F,V1702A02B_F,V1702A02B_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_OX", 2, (/ O3,NO2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_ECCOARSE", 2, (/ EC_C_WOOD,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "NVWOODOC25", 1, (/ POM_F_WOOD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "ASH", 2, (/ V1702A02B_F,V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_NVFFUELOC_COARSE", 1, (/ POM_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PM10", 26, (/ SO4,NO3_F,NO3_C,NH4_F,V1702A02B_F,V1702A02B_C,PART_OM_F,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0 /) ) &
, gtype( "WDEP_PWOODOA25", 1, (/ WOODOA_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_ASOA", 10, (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "OX", 2, (/ O3,NO2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_OXN", 7, (/ NO2,PAN,MPAN,HNO3,HONO,NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PFFUELOA25", 1, (/ FFFUEL_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PPM10", 11, (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PPM_C", 4, (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PPM25_FIRE", 3, (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_EC_F", 5, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PM10", 25, (/ SO4,NO3_F,NO3_C,NH4_F,V1702A02B_F,V1702A02B_C,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_SIA", 4, (/ SO4,NO3_F,NO3_C,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "BVOC", 2, (/ C5H8,APINENE,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PWOODOA25", 1, (/ WOODOA_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "EC_F", 5, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_NONVOLPCM", 10, (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_NVFFIREOC25", 1, (/ FFIRE_OM,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "SOX", 2, (/ SO2,SO4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DUST_ANT_F", 1, (/ DUST_ROAD_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PMFINE", 16, (/ SO4,NO3_F,NH4_F,V1702A02B_F,PART_OM_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_NVABSOM", 3, (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_SVFFIREOA25", 1, (/ FFIREOA_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_AOD", 13, (/ SO4,NO3_F,NH4_F,V1702A02B_F,V1702A02B_C,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "NONVOLPCM", 10, (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_V1702A02B", 2, (/ V1702A02B_F,V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_ASH_F", 1, (/ V1702A02B_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_DUST_NAT_C", 2, (/ DUST_WB_C,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PMCO", 10, (/ NO3_C,V1702A02B_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_DUST_ANT_F", 1, (/ DUST_ROAD_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_OMCOARSE", 1, (/ POM_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "NVFFIREOC25", 1, (/ FFIRE_OM,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_RDN", 2, (/ NH3,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_ASOA", 10, (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_NVWOODOC25", 1, (/ POM_F_WOOD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_ASH_C", 1, (/ V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_NVWOODOC25", 1, (/ POM_F_WOOD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_ECFINE", 4, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_ECCOARSE", 2, (/ EC_C_WOOD,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "OXN", 13, (/ NO,NO2,PAN,MPAN,NO3,N2O5,ISONO3,HNO3,HONO,ISNI,ISNIR,NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_DUST_ANT_F", 1, (/ DUST_ROAD_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_NVFFUELOC25", 1, (/ POM_F_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "SIA", 4, (/ SO4,NO3_F,NO3_C,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_ASH", 2, (/ V1702A02B_F,V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_NVFFIREOC25", 1, (/ FFIRE_OM,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_TNO3", 2, (/ NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_ASH_F", 1, (/ V1702A02B_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PMCO", 10, (/ NO3_C,V1702A02B_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_BSOA", 10, (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PCM", 28, (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,0,0,0 /) ) &
, gtype( "DDEP_NVFFUELOC_COARSE", 1, (/ POM_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "ECCOARSE", 2, (/ EC_C_WOOD,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WOODEC", 3, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_TNO3", 2, (/ NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_DUST_ANT_C", 1, (/ DUST_ROAD_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PMCO", 10, (/ NO3_C,V1702A02B_C,POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_FFUELEC", 3, (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_ASH_C", 1, (/ V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_NVFFUELOC25", 1, (/ POM_F_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PM10ANTHR", 3, (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_SVFFUELOA25", 1, (/ FFFUEL_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PPM25", 7, (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PPM25_FIRE", 3, (/ FFIRE_OM,FFIRE_BC,FFIRE_REMPPM25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "FFIREBC", 1, (/ FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_FFIREBC", 1, (/ FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "NOX", 2, (/ NO,NO2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DUST_NAT_F", 2, (/ DUST_WB_F,DUST_SAH_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "SS", 2, (/ SEASALT_F,SEASALT_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_DUST", 6, (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "RO2", 14, (/ HO2,CH3O2,C2H5O2,SECC4H9O2,ISRO2,ETRO2,PRRO2,OXYO2,MEKO2,MALO2,MVKO2,MACRO2,MACO3,TERPPEROXY,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_EC_F", 5, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "ROOH", 16, (/ CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,OXYO2H,MEKO2H,MALO2H,MVKO2H,MACROOH,MACO3H,ISRO2H,H2O2,CH3COO2H,ISONO3H,ISNIRH,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DUST_ANT_C", 1, (/ DUST_ROAD_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "AOD", 13, (/ SO4,NO3_F,NH4_F,V1702A02B_F,V1702A02B_C,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PPM_C", 4, (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,REMPPM_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_NVABSOM", 3, (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "SVFFUELOA25", 1, (/ FFFUEL_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_SOX", 2, (/ SO2,SO4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_DUST_NAT_F", 2, (/ DUST_WB_F,DUST_SAH_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "ASH_F", 1, (/ V1702A02B_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_SOX", 2, (/ SO2,SO4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PFFUELOA25", 1, (/ FFFUEL_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "FFUELEC", 3, (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PM10ANTHR", 3, (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_DUST_ANT_C", 1, (/ DUST_ROAD_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_DUST_NAT_F", 2, (/ DUST_WB_F,DUST_SAH_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PM10ANTHR", 3, (/ POM_C_FFUEL,EC_C_WOOD,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_WOODEC", 3, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "V1702A02B", 2, (/ V1702A02B_F,V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_WOODEC", 3, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PCM", 31, (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /) ) &
, gtype( "WDEP_SS", 2, (/ SEASALT_F,SEASALT_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DUST_NAT_C", 2, (/ DUST_WB_C,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_NONVOLPCM", 10, (/ POM_F_WOOD,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_OM,FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PPM25", 7, (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "ASOA", 10, (/ ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PPM10", 11, (/ POM_F_WOOD,POM_F_FFUEL,POM_C_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,REMPPM25,REMPPM_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "NVABSOM", 3, (/ POM_F_WOOD,POM_F_FFUEL,FFIRE_OM,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "BSOA", 10, (/ BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "ECFINE", 4, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_FFIREBC", 1, (/ FFIRE_BC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_DUST", 6, (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "NVFFUELOC_COARSE", 1, (/ POM_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_ECFINE", 4, (/ EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_ROOH", 1, (/ H2O2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PCM", 31, (/ POM_F_WOOD,POM_F_FFUEL,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,FFIRE_OM,FFIRE_BC,ASOC_NG100,ASOC_UG1,ASOC_UG10,ASOC_UG1E2,ASOC_UG1E3,NON_C_ASOA_NG100,NON_C_ASOA_UG1,NON_C_ASOA_UG10,NON_C_ASOA_UG1E2,NON_C_ASOA_UG1E3,BSOC_NG100,BSOC_UG1,BSOC_UG10,BSOC_UG1E2,BSOC_UG1E3,NON_C_BSOA_NG100,NON_C_BSOA_UG1,NON_C_BSOA_UG10,NON_C_BSOA_UG1E2,NON_C_BSOA_UG1E3,FFFUEL_NG10,WOODOA_NG10,FFIREOA_NG10 /) ) &
, gtype( "SVWOODOA25", 1, (/ WOODOA_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_SIA", 4, (/ SO4,NO3_F,NO3_C,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PFFIREOA25", 1, (/ FFIREOA_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PCM_HELP", 20, (/ GAS_ASOA_OC,PART_ASOA_OC,PART_ASOA_OM,GAS_BSOA_OC,PART_BSOA_OC,PART_BSOA_OM,PART_FFUELOA25_OC,PART_FFUELOA25_OM,PART_WOODOA25_OC,PART_WOODOA25_OM,PART_FFIREOA25_OC,PART_FFIREOA25_OM,PART_OC10,PART_OC25,NONVOL_FFUELOC25,NONV_FFUELOC_COARSE,NONVOL_WOODOC25,NONVOL_BGNDOC,NONVOL_FFIREOC25,PART_OM_F,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_OMCOARSE", 1, (/ POM_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "TNO3", 2, (/ NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "SVFFIREOA25", 1, (/ FFIREOA_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "NVFFUELOC25", 1, (/ POM_F_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_V1702A02B", 2, (/ V1702A02B_F,V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "OMCOARSE", 1, (/ POM_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_ROOH", 3, (/ CH3O2H,C2H5OOH,H2O2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PMFINE", 15, (/ SO4,NO3_F,NH4_F,V1702A02B_F,EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,REMPPM25,FFIRE_BC,FFIRE_REMPPM25,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_RDN", 2, (/ NH3,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PFFIREOA25", 1, (/ FFIREOA_NG10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "RDN", 2, (/ NH3,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_ASH", 2, (/ V1702A02B_F,V1702A02B_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_FFUELEC", 3, (/ EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
   /)

! ------- Dry dep      species ------------------
  integer, public, parameter, dimension(7) :: &
               DDEP_OXNGROUP = (/ NO2,PAN,MPAN,HNO3,HONO,NO3_f,NO3_c /)
  integer, public, parameter, dimension(2) :: &
               DDEP_SOXGROUP = (/ SO2,SO4 /)
  integer, public, parameter, dimension(2) :: &
               DDEP_RDNGROUP = (/ NH3,NH4_f /)

  integer, public, parameter :: NMAX_DDEP = 7


! ------- Wet dep      species ------------------
  integer, public, parameter, dimension(4) :: &
               WDEP_OXNGROUP = (/ HNO3,HONO,NO3_f,NO3_c /)
  integer, public, parameter, dimension(2) :: &
               WDEP_SOXGROUP = (/ SO2,SO4 /)
  integer, public, parameter, dimension(2) :: &
               WDEP_SSALTGROUP = (/ SeaSalt_f,SeaSalt_c /)
  integer, public, parameter, dimension(2) :: &
               WDEP_RDNGROUP = (/ NH3,NH4_f /)

  integer, public, parameter :: NMAX_WDEP = 4


! ------- RO2 Pool     species ------------------
  integer, public, parameter :: SIZE_RO2_POOL      = 1
  integer, public, parameter, dimension(1) :: &
     RO2_POOL      = (/ -99 /)
   type(typ_sp), dimension(142), public, save :: chemgroups


!-----------------------------------------------------------
  contains
 subroutine Init_ChemGroups()

   integer, dimension(:), pointer :: p

 p => DDEP_SS_GROUP 
 chemgroups(1) = typ_sp("DDEP_SS", p )

 p => WDEP_OXN_GROUP 
 chemgroups(2) = typ_sp("WDEP_OXN", p )

 p => WDEP_PPM10_GROUP 
 chemgroups(3) = typ_sp("WDEP_PPM10", p )

 p => ASH_C_GROUP 
 chemgroups(4) = typ_sp("ASH_C", p )

 p => DUST_GROUP 
 chemgroups(5) = typ_sp("DUST", p )

 p => WDEP_AOD_GROUP 
 chemgroups(6) = typ_sp("WDEP_AOD", p )

 p => WDEP_BSOA_GROUP 
 chemgroups(7) = typ_sp("WDEP_BSOA", p )

 p => DDEP_NOX_GROUP 
 chemgroups(8) = typ_sp("DDEP_NOX", p )

 p => PPM_C_GROUP 
 chemgroups(9) = typ_sp("PPM_C", p )

 p => DDEP_DUST_NAT_C_GROUP 
 chemgroups(10) = typ_sp("DDEP_DUST_NAT_C", p )

 p => PPM25_FIRE_GROUP 
 chemgroups(11) = typ_sp("PPM25_FIRE", p )

 p => WDEP_SVWOODOA25_GROUP 
 chemgroups(12) = typ_sp("WDEP_SVWOODOA25", p )

 p => DDEP_PMFINE_GROUP 
 chemgroups(13) = typ_sp("DDEP_PMFINE", p )

 p => WDEP_PPM25_GROUP 
 chemgroups(14) = typ_sp("WDEP_PPM25", p )

 p => WDEP_PM10_GROUP 
 chemgroups(15) = typ_sp("WDEP_PM10", p )

 p => DDEP_OX_GROUP 
 chemgroups(16) = typ_sp("DDEP_OX", p )

 p => DDEP_ECCOARSE_GROUP 
 chemgroups(17) = typ_sp("DDEP_ECCOARSE", p )

 p => NVWOODOC25_GROUP 
 chemgroups(18) = typ_sp("NVWOODOC25", p )

 p => ASH_GROUP 
 chemgroups(19) = typ_sp("ASH", p )

 p => WDEP_NVFFUELOC_COARSE_GROUP 
 chemgroups(20) = typ_sp("WDEP_NVFFUELOC_COARSE", p )

 p => PM10_GROUP 
 chemgroups(21) = typ_sp("PM10", p )

 p => WDEP_PWOODOA25_GROUP 
 chemgroups(22) = typ_sp("WDEP_PWOODOA25", p )

 p => DDEP_ASOA_GROUP 
 chemgroups(23) = typ_sp("DDEP_ASOA", p )

 p => OX_GROUP 
 chemgroups(24) = typ_sp("OX", p )

 p => DDEP_OXN_GROUP 
 chemgroups(25) = typ_sp("DDEP_OXN", p )

 p => WDEP_PFFUELOA25_GROUP 
 chemgroups(26) = typ_sp("WDEP_PFFUELOA25", p )

 p => DDEP_PPM10_GROUP 
 chemgroups(27) = typ_sp("DDEP_PPM10", p )

 p => WDEP_PPM_C_GROUP 
 chemgroups(28) = typ_sp("WDEP_PPM_C", p )

 p => DDEP_PPM25_FIRE_GROUP 
 chemgroups(29) = typ_sp("DDEP_PPM25_FIRE", p )

 p => WDEP_EC_F_GROUP 
 chemgroups(30) = typ_sp("WDEP_EC_F", p )

 p => DDEP_PM10_GROUP 
 chemgroups(31) = typ_sp("DDEP_PM10", p )

 p => WDEP_SIA_GROUP 
 chemgroups(32) = typ_sp("WDEP_SIA", p )

 p => BVOC_GROUP 
 chemgroups(33) = typ_sp("BVOC", p )

 p => PWOODOA25_GROUP 
 chemgroups(34) = typ_sp("PWOODOA25", p )

 p => EC_F_GROUP 
 chemgroups(35) = typ_sp("EC_F", p )

 p => DDEP_NONVOLPCM_GROUP 
 chemgroups(36) = typ_sp("DDEP_NONVOLPCM", p )

 p => WDEP_NVFFIREOC25_GROUP 
 chemgroups(37) = typ_sp("WDEP_NVFFIREOC25", p )

 p => SOX_GROUP 
 chemgroups(38) = typ_sp("SOX", p )

 p => DUST_ANT_F_GROUP 
 chemgroups(39) = typ_sp("DUST_ANT_F", p )

 p => PMFINE_GROUP 
 chemgroups(40) = typ_sp("PMFINE", p )

 p => DDEP_NVABSOM_GROUP 
 chemgroups(41) = typ_sp("DDEP_NVABSOM", p )

 p => WDEP_SVFFIREOA25_GROUP 
 chemgroups(42) = typ_sp("WDEP_SVFFIREOA25", p )

 p => DDEP_AOD_GROUP 
 chemgroups(43) = typ_sp("DDEP_AOD", p )

 p => NONVOLPCM_GROUP 
 chemgroups(44) = typ_sp("NONVOLPCM", p )

 p => WDEP_V1702A02B_GROUP 
 chemgroups(45) = typ_sp("WDEP_V1702A02B", p )

 p => WDEP_ASH_F_GROUP 
 chemgroups(46) = typ_sp("WDEP_ASH_F", p )

 p => WDEP_DUST_NAT_C_GROUP 
 chemgroups(47) = typ_sp("WDEP_DUST_NAT_C", p )

 p => PMCO_GROUP 
 chemgroups(48) = typ_sp("PMCO", p )

 p => DDEP_DUST_ANT_F_GROUP 
 chemgroups(49) = typ_sp("DDEP_DUST_ANT_F", p )

 p => DDEP_OMCOARSE_GROUP 
 chemgroups(50) = typ_sp("DDEP_OMCOARSE", p )

 p => NVFFIREOC25_GROUP 
 chemgroups(51) = typ_sp("NVFFIREOC25", p )

 p => WDEP_RDN_GROUP 
 chemgroups(52) = typ_sp("WDEP_RDN", p )

 p => WDEP_ASOA_GROUP 
 chemgroups(53) = typ_sp("WDEP_ASOA", p )

 p => WDEP_NVWOODOC25_GROUP 
 chemgroups(54) = typ_sp("WDEP_NVWOODOC25", p )

 p => WDEP_ASH_C_GROUP 
 chemgroups(55) = typ_sp("WDEP_ASH_C", p )

 p => DDEP_NVWOODOC25_GROUP 
 chemgroups(56) = typ_sp("DDEP_NVWOODOC25", p )

 p => WDEP_ECFINE_GROUP 
 chemgroups(57) = typ_sp("WDEP_ECFINE", p )

 p => WDEP_ECCOARSE_GROUP 
 chemgroups(58) = typ_sp("WDEP_ECCOARSE", p )

 p => OXN_GROUP 
 chemgroups(59) = typ_sp("OXN", p )

 p => WDEP_DUST_ANT_F_GROUP 
 chemgroups(60) = typ_sp("WDEP_DUST_ANT_F", p )

 p => DDEP_NVFFUELOC25_GROUP 
 chemgroups(61) = typ_sp("DDEP_NVFFUELOC25", p )

 p => SIA_GROUP 
 chemgroups(62) = typ_sp("SIA", p )

 p => DDEP_ASH_GROUP 
 chemgroups(63) = typ_sp("DDEP_ASH", p )

 p => DDEP_NVFFIREOC25_GROUP 
 chemgroups(64) = typ_sp("DDEP_NVFFIREOC25", p )

 p => DDEP_TNO3_GROUP 
 chemgroups(65) = typ_sp("DDEP_TNO3", p )

 p => DDEP_ASH_F_GROUP 
 chemgroups(66) = typ_sp("DDEP_ASH_F", p )

 p => DDEP_PMCO_GROUP 
 chemgroups(67) = typ_sp("DDEP_PMCO", p )

 p => DDEP_BSOA_GROUP 
 chemgroups(68) = typ_sp("DDEP_BSOA", p )

 p => DDEP_PCM_GROUP 
 chemgroups(69) = typ_sp("DDEP_PCM", p )

 p => DDEP_NVFFUELOC_COARSE_GROUP 
 chemgroups(70) = typ_sp("DDEP_NVFFUELOC_COARSE", p )

 p => ECCOARSE_GROUP 
 chemgroups(71) = typ_sp("ECCOARSE", p )

 p => WOODEC_GROUP 
 chemgroups(72) = typ_sp("WOODEC", p )

 p => WDEP_TNO3_GROUP 
 chemgroups(73) = typ_sp("WDEP_TNO3", p )

 p => DDEP_DUST_ANT_C_GROUP 
 chemgroups(74) = typ_sp("DDEP_DUST_ANT_C", p )

 p => WDEP_PMCO_GROUP 
 chemgroups(75) = typ_sp("WDEP_PMCO", p )

 p => DDEP_FFUELEC_GROUP 
 chemgroups(76) = typ_sp("DDEP_FFUELEC", p )

 p => DDEP_ASH_C_GROUP 
 chemgroups(77) = typ_sp("DDEP_ASH_C", p )

 p => WDEP_NVFFUELOC25_GROUP 
 chemgroups(78) = typ_sp("WDEP_NVFFUELOC25", p )

 p => WDEP_PM10ANTHR_GROUP 
 chemgroups(79) = typ_sp("WDEP_PM10ANTHR", p )

 p => WDEP_SVFFUELOA25_GROUP 
 chemgroups(80) = typ_sp("WDEP_SVFFUELOA25", p )

 p => DDEP_PPM25_GROUP 
 chemgroups(81) = typ_sp("DDEP_PPM25", p )

 p => WDEP_PPM25_FIRE_GROUP 
 chemgroups(82) = typ_sp("WDEP_PPM25_FIRE", p )

 p => FFIREBC_GROUP 
 chemgroups(83) = typ_sp("FFIREBC", p )

 p => WDEP_FFIREBC_GROUP 
 chemgroups(84) = typ_sp("WDEP_FFIREBC", p )

 p => NOX_GROUP 
 chemgroups(85) = typ_sp("NOX", p )

 p => DUST_NAT_F_GROUP 
 chemgroups(86) = typ_sp("DUST_NAT_F", p )

 p => SS_GROUP 
 chemgroups(87) = typ_sp("SS", p )

 p => DDEP_DUST_GROUP 
 chemgroups(88) = typ_sp("DDEP_DUST", p )

 p => RO2_GROUP 
 chemgroups(89) = typ_sp("RO2", p )

 p => DDEP_EC_F_GROUP 
 chemgroups(90) = typ_sp("DDEP_EC_F", p )

 p => ROOH_GROUP 
 chemgroups(91) = typ_sp("ROOH", p )

 p => DUST_ANT_C_GROUP 
 chemgroups(92) = typ_sp("DUST_ANT_C", p )

 p => AOD_GROUP 
 chemgroups(93) = typ_sp("AOD", p )

 p => DDEP_PPM_C_GROUP 
 chemgroups(94) = typ_sp("DDEP_PPM_C", p )

 p => WDEP_NVABSOM_GROUP 
 chemgroups(95) = typ_sp("WDEP_NVABSOM", p )

 p => SVFFUELOA25_GROUP 
 chemgroups(96) = typ_sp("SVFFUELOA25", p )

 p => DDEP_SOX_GROUP 
 chemgroups(97) = typ_sp("DDEP_SOX", p )

 p => WDEP_DUST_NAT_F_GROUP 
 chemgroups(98) = typ_sp("WDEP_DUST_NAT_F", p )

 p => ASH_F_GROUP 
 chemgroups(99) = typ_sp("ASH_F", p )

 p => WDEP_SOX_GROUP 
 chemgroups(100) = typ_sp("WDEP_SOX", p )

 p => PFFUELOA25_GROUP 
 chemgroups(101) = typ_sp("PFFUELOA25", p )

 p => FFUELEC_GROUP 
 chemgroups(102) = typ_sp("FFUELEC", p )

 p => PM10ANTHR_GROUP 
 chemgroups(103) = typ_sp("PM10ANTHR", p )

 p => WDEP_DUST_ANT_C_GROUP 
 chemgroups(104) = typ_sp("WDEP_DUST_ANT_C", p )

 p => DDEP_DUST_NAT_F_GROUP 
 chemgroups(105) = typ_sp("DDEP_DUST_NAT_F", p )

 p => DDEP_PM10ANTHR_GROUP 
 chemgroups(106) = typ_sp("DDEP_PM10ANTHR", p )

 p => DDEP_WOODEC_GROUP 
 chemgroups(107) = typ_sp("DDEP_WOODEC", p )

 p => V1702A02B_GROUP 
 chemgroups(108) = typ_sp("V1702A02B", p )

 p => WDEP_WOODEC_GROUP 
 chemgroups(109) = typ_sp("WDEP_WOODEC", p )

 p => WDEP_PCM_GROUP 
 chemgroups(110) = typ_sp("WDEP_PCM", p )

 p => WDEP_SS_GROUP 
 chemgroups(111) = typ_sp("WDEP_SS", p )

 p => DUST_NAT_C_GROUP 
 chemgroups(112) = typ_sp("DUST_NAT_C", p )

 p => WDEP_NONVOLPCM_GROUP 
 chemgroups(113) = typ_sp("WDEP_NONVOLPCM", p )

 p => PPM25_GROUP 
 chemgroups(114) = typ_sp("PPM25", p )

 p => ASOA_GROUP 
 chemgroups(115) = typ_sp("ASOA", p )

 p => PPM10_GROUP 
 chemgroups(116) = typ_sp("PPM10", p )

 p => NVABSOM_GROUP 
 chemgroups(117) = typ_sp("NVABSOM", p )

 p => BSOA_GROUP 
 chemgroups(118) = typ_sp("BSOA", p )

 p => ECFINE_GROUP 
 chemgroups(119) = typ_sp("ECFINE", p )

 p => DDEP_FFIREBC_GROUP 
 chemgroups(120) = typ_sp("DDEP_FFIREBC", p )

 p => WDEP_DUST_GROUP 
 chemgroups(121) = typ_sp("WDEP_DUST", p )

 p => NVFFUELOC_COARSE_GROUP 
 chemgroups(122) = typ_sp("NVFFUELOC_COARSE", p )

 p => DDEP_ECFINE_GROUP 
 chemgroups(123) = typ_sp("DDEP_ECFINE", p )

 p => WDEP_ROOH_GROUP 
 chemgroups(124) = typ_sp("WDEP_ROOH", p )

 p => PCM_GROUP 
 chemgroups(125) = typ_sp("PCM", p )

 p => SVWOODOA25_GROUP 
 chemgroups(126) = typ_sp("SVWOODOA25", p )

 p => DDEP_SIA_GROUP 
 chemgroups(127) = typ_sp("DDEP_SIA", p )

 p => WDEP_PFFIREOA25_GROUP 
 chemgroups(128) = typ_sp("WDEP_PFFIREOA25", p )

 p => PCM_HELP_GROUP 
 chemgroups(129) = typ_sp("PCM_HELP", p )

 p => WDEP_OMCOARSE_GROUP 
 chemgroups(130) = typ_sp("WDEP_OMCOARSE", p )

 p => TNO3_GROUP 
 chemgroups(131) = typ_sp("TNO3", p )

 p => SVFFIREOA25_GROUP 
 chemgroups(132) = typ_sp("SVFFIREOA25", p )

 p => NVFFUELOC25_GROUP 
 chemgroups(133) = typ_sp("NVFFUELOC25", p )

 p => DDEP_V1702A02B_GROUP 
 chemgroups(134) = typ_sp("DDEP_V1702A02B", p )

 p => OMCOARSE_GROUP 
 chemgroups(135) = typ_sp("OMCOARSE", p )

 p => DDEP_ROOH_GROUP 
 chemgroups(136) = typ_sp("DDEP_ROOH", p )

 p => WDEP_PMFINE_GROUP 
 chemgroups(137) = typ_sp("WDEP_PMFINE", p )

 p => DDEP_RDN_GROUP 
 chemgroups(138) = typ_sp("DDEP_RDN", p )

 p => PFFIREOA25_GROUP 
 chemgroups(139) = typ_sp("PFFIREOA25", p )

 p => RDN_GROUP 
 chemgroups(140) = typ_sp("RDN", p )

 p => WDEP_ASH_GROUP 
 chemgroups(141) = typ_sp("WDEP_ASH", p )

 p => WDEP_FFUELEC_GROUP 
 chemgroups(142) = typ_sp("WDEP_FFUELEC", p )

   nullify(p)


 end subroutine Init_ChemGroups
 !-----------------------------------------------------------


 end module ChemGroups_ml
 !-----------------------------------------------------------
