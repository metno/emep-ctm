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

integer, public, target, save, dimension(1) :: &
  DDEP_DUST_ANT_C_GROUP = (/ DUST_ROAD_C /)

integer, public, target, save, dimension(3) :: &
  WDEP_OXN_GROUP = (/ HNO3,NO3_F,NO3_C /)

integer, public, target, save, dimension(10) :: &
  WDEP_PMCO_GROUP = (/ NO3_C,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(5) :: &
  DDEP_ASH_C_GROUP = (/ V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(5) :: &
  ASH_C_GROUP = (/ V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(6) :: &
  DUST_GROUP = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(2) :: &
  DDEP_NOX_GROUP = (/ NO2,SHIPNOX /)

integer, public, target, save, dimension(3) :: &
  NOX_GROUP = (/ NO,NO2,SHIPNOX /)

integer, public, target, save, dimension(2) :: &
  DDEP_DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(2) :: &
  DUST_NAT_F_GROUP = (/ DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(2) :: &
  SS_GROUP = (/ SEASALT_F,SEASALT_C /)

integer, public, target, save, dimension(6) :: &
  DDEP_DUST_GROUP = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(12) :: &
  DDEP_PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,PPM25_FIRE,V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(2) :: &
  RO2_GROUP = (/ HO2,RO2 /)

integer, public, target, save, dimension(1) :: &
  DUST_ANT_C_GROUP = (/ DUST_ROAD_C /)

integer, public, target, save, dimension(1) :: &
  AROMATIC_GROUP = (/ OXYL /)

integer, public, target, save, dimension(15) :: &
  NMVOC_GROUP = (/ VOC,CH3CHO,PAN,CH3COO2,HCHO,C2H6,NC4H10,OXYL,C5H8,APINENE,C3H6,C2H4,C2H5OH,CH3OH,UNREAC /)

integer, public, target, save, dimension(22) :: &
  AOD_GROUP = (/ DUST_ROAD_C,DUST_ROAD_F,DUST_SAH_C,DUST_SAH_F,DUST_WB_C,DUST_WB_F,NH4_F,NO3_C,NO3_F,PPM25_FIRE,SEASALT_C,SEASALT_F,SO4,V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(3) :: &
  DDEP_OX_GROUP = (/ O3,NO2,SHIPNOX /)

integer, public, target, save, dimension(2) :: &
  DDEP_SOX_GROUP = (/ SO2,SO4 /)

integer, public, target, save, dimension(22) :: &
  WDEP_PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,PPM25_FIRE,V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(2) :: &
  WDEP_DUST_NAT_F_GROUP = (/ DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(4) :: &
  ASH_F_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4 /)

integer, public, target, save, dimension(2) :: &
  WDEP_SOX_GROUP = (/ SO2,SO4 /)

integer, public, target, save, dimension(9) :: &
  ASH_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(22) :: &
  PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,PPM25_FIRE,V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(1) :: &
  WDEP_DUST_ANT_C_GROUP = (/ DUST_ROAD_C /)

integer, public, target, save, dimension(2) :: &
  DDEP_DUST_NAT_F_GROUP = (/ DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(3) :: &
  OX_GROUP = (/ O3,NO2,SHIPNOX /)

integer, public, target, save, dimension(6) :: &
  DDEP_OXN_GROUP = (/ NO2,SHIPNOX,HNO3,PAN,NO3_F,NO3_C /)

integer, public, target, save, dimension(9) :: &
  V1702A02B_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(22) :: &
  DDEP_PM10_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F,PPM25_FIRE,V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9,SEASALT_F,SEASALT_C,DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(4) :: &
  WDEP_SIA_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F /)

integer, public, target, save, dimension(2) :: &
  BVOC_GROUP = (/ C5H8,APINENE /)

integer, public, target, save, dimension(2) :: &
  WDEP_SS_GROUP = (/ SEASALT_F,SEASALT_C /)

integer, public, target, save, dimension(2) :: &
  DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(2) :: &
  SOX_GROUP = (/ SO2,SO4 /)

integer, public, target, save, dimension(12) :: &
  PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,PPM25_FIRE,V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(1) :: &
  DUST_ANT_F_GROUP = (/ DUST_ROAD_F /)

integer, public, target, save, dimension(9) :: &
  WDEP_V1702A02B_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(1) :: &
  HOX_GROUP = (/ H2 /)

integer, public, target, save, dimension(6) :: &
  WDEP_DUST_GROUP = (/ DUST_ROAD_F,DUST_ROAD_C,DUST_WB_F,DUST_WB_C,DUST_SAH_F,DUST_SAH_C /)

integer, public, target, save, dimension(4) :: &
  WDEP_ASH_F_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4 /)

integer, public, target, save, dimension(10) :: &
  PMCO_GROUP = (/ NO3_C,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(1) :: &
  ALKENE_GROUP = (/ C3H6 /)

integer, public, target, save, dimension(2) :: &
  WDEP_DUST_NAT_C_GROUP = (/ DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(1) :: &
  DDEP_DUST_ANT_F_GROUP = (/ DUST_ROAD_F /)

integer, public, target, save, dimension(2) :: &
  WDEP_RDN_GROUP = (/ NH3,NH4_F /)

integer, public, target, save, dimension(4) :: &
  DDEP_SIA_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F /)

integer, public, target, save, dimension(5) :: &
  WDEP_ASH_C_GROUP = (/ V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(7) :: &
  OXN_GROUP = (/ NO,NO2,SHIPNOX,HNO3,PAN,NO3_F,NO3_C /)

integer, public, target, save, dimension(2) :: &
  TNO3_GROUP = (/ NO3_F,NO3_C /)

integer, public, target, save, dimension(1) :: &
  DDEP_DAOBS_GROUP = (/ SHIPNOX /)

integer, public, target, save, dimension(9) :: &
  DDEP_V1702A02B_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(1) :: &
  WDEP_DUST_ANT_F_GROUP = (/ DUST_ROAD_F /)

integer, public, target, save, dimension(4) :: &
  SIA_GROUP = (/ SO4,NO3_F,NO3_C,NH4_F /)

integer, public, target, save, dimension(9) :: &
  DDEP_ASH_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(12) :: &
  WDEP_PMFINE_GROUP = (/ SO4,NO3_F,NH4_F,PPM25_FIRE,V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,SEASALT_F,DUST_ROAD_F,DUST_WB_F,DUST_SAH_F /)

integer, public, target, save, dimension(2) :: &
  DDEP_TNO3_GROUP = (/ NO3_F,NO3_C /)

integer, public, target, save, dimension(2) :: &
  DDEP_RDN_GROUP = (/ NH3,NH4_F /)

integer, public, target, save, dimension(4) :: &
  DDEP_ASH_F_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4 /)

integer, public, target, save, dimension(10) :: &
  DDEP_PMCO_GROUP = (/ NO3_C,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9,SEASALT_C,DUST_ROAD_C,DUST_WB_C,DUST_SAH_C /)

integer, public, target, save, dimension(2) :: &
  RDN_GROUP = (/ NH3,NH4_F /)

integer, public, target, save, dimension(9) :: &
  WDEP_ASH_GROUP = (/ V1702A02B_1,V1702A02B_2,V1702A02B_3,V1702A02B_4,V1702A02B_5,V1702A02B_6,V1702A02B_7,V1702A02B_8,V1702A02B_9 /)

integer, public, target, save, dimension(1) :: &
  DAOBS_GROUP = (/ SHIPNOX /)

integer, public, target, save, dimension(2) :: &
  WDEP_TNO3_GROUP = (/ NO3_F,NO3_C /)


! ------- RO2 Pool     species ------------------
integer, public, parameter :: SIZE_RO2_POOL      = 1
integer, public, parameter, dimension(1) :: &
  RO2_POOL      = (/ -99 /)
type(typ_sp), dimension(67), public, save :: chemgroups


!-----------------------------------------------------------
contains
subroutine Init_ChemGroups()

  chemgroups(1)%name="DDEP_SS"
  chemgroups(1)%ptr=>DDEP_SS_GROUP

  chemgroups(2)%name="DDEP_DUST_ANT_C"
  chemgroups(2)%ptr=>DDEP_DUST_ANT_C_GROUP

  chemgroups(3)%name="WDEP_OXN"
  chemgroups(3)%ptr=>WDEP_OXN_GROUP

  chemgroups(4)%name="WDEP_PMCO"
  chemgroups(4)%ptr=>WDEP_PMCO_GROUP

  chemgroups(5)%name="DDEP_ASH_C"
  chemgroups(5)%ptr=>DDEP_ASH_C_GROUP

  chemgroups(6)%name="ASH_C"
  chemgroups(6)%ptr=>ASH_C_GROUP

  chemgroups(7)%name="DUST"
  chemgroups(7)%ptr=>DUST_GROUP

  chemgroups(8)%name="DDEP_NOX"
  chemgroups(8)%ptr=>DDEP_NOX_GROUP

  chemgroups(9)%name="NOX"
  chemgroups(9)%ptr=>NOX_GROUP

  chemgroups(10)%name="DDEP_DUST_NAT_C"
  chemgroups(10)%ptr=>DDEP_DUST_NAT_C_GROUP

  chemgroups(11)%name="DUST_NAT_F"
  chemgroups(11)%ptr=>DUST_NAT_F_GROUP

  chemgroups(12)%name="SS"
  chemgroups(12)%ptr=>SS_GROUP

  chemgroups(13)%name="DDEP_DUST"
  chemgroups(13)%ptr=>DDEP_DUST_GROUP

  chemgroups(14)%name="DDEP_PMFINE"
  chemgroups(14)%ptr=>DDEP_PMFINE_GROUP

  chemgroups(15)%name="RO2"
  chemgroups(15)%ptr=>RO2_GROUP

  chemgroups(16)%name="DUST_ANT_C"
  chemgroups(16)%ptr=>DUST_ANT_C_GROUP

  chemgroups(17)%name="AROMATIC"
  chemgroups(17)%ptr=>AROMATIC_GROUP

  chemgroups(18)%name="NMVOC"
  chemgroups(18)%ptr=>NMVOC_GROUP

  chemgroups(19)%name="AOD"
  chemgroups(19)%ptr=>AOD_GROUP

  chemgroups(20)%name="DDEP_OX"
  chemgroups(20)%ptr=>DDEP_OX_GROUP

  chemgroups(21)%name="DDEP_SOX"
  chemgroups(21)%ptr=>DDEP_SOX_GROUP

  chemgroups(22)%name="WDEP_PM10"
  chemgroups(22)%ptr=>WDEP_PM10_GROUP

  chemgroups(23)%name="WDEP_DUST_NAT_F"
  chemgroups(23)%ptr=>WDEP_DUST_NAT_F_GROUP

  chemgroups(24)%name="ASH_F"
  chemgroups(24)%ptr=>ASH_F_GROUP

  chemgroups(25)%name="WDEP_SOX"
  chemgroups(25)%ptr=>WDEP_SOX_GROUP

  chemgroups(26)%name="ASH"
  chemgroups(26)%ptr=>ASH_GROUP

  chemgroups(27)%name="PM10"
  chemgroups(27)%ptr=>PM10_GROUP

  chemgroups(28)%name="WDEP_DUST_ANT_C"
  chemgroups(28)%ptr=>WDEP_DUST_ANT_C_GROUP

  chemgroups(29)%name="DDEP_DUST_NAT_F"
  chemgroups(29)%ptr=>DDEP_DUST_NAT_F_GROUP

  chemgroups(30)%name="OX"
  chemgroups(30)%ptr=>OX_GROUP

  chemgroups(31)%name="DDEP_OXN"
  chemgroups(31)%ptr=>DDEP_OXN_GROUP

  chemgroups(32)%name="V1702A02B"
  chemgroups(32)%ptr=>V1702A02B_GROUP

  chemgroups(33)%name="DDEP_PM10"
  chemgroups(33)%ptr=>DDEP_PM10_GROUP

  chemgroups(34)%name="WDEP_SIA"
  chemgroups(34)%ptr=>WDEP_SIA_GROUP

  chemgroups(35)%name="BVOC"
  chemgroups(35)%ptr=>BVOC_GROUP

  chemgroups(36)%name="WDEP_SS"
  chemgroups(36)%ptr=>WDEP_SS_GROUP

  chemgroups(37)%name="DUST_NAT_C"
  chemgroups(37)%ptr=>DUST_NAT_C_GROUP

  chemgroups(38)%name="SOX"
  chemgroups(38)%ptr=>SOX_GROUP

  chemgroups(39)%name="PMFINE"
  chemgroups(39)%ptr=>PMFINE_GROUP

  chemgroups(40)%name="DUST_ANT_F"
  chemgroups(40)%ptr=>DUST_ANT_F_GROUP

  chemgroups(41)%name="WDEP_V1702A02B"
  chemgroups(41)%ptr=>WDEP_V1702A02B_GROUP

  chemgroups(42)%name="HOX"
  chemgroups(42)%ptr=>HOX_GROUP

  chemgroups(43)%name="WDEP_DUST"
  chemgroups(43)%ptr=>WDEP_DUST_GROUP

  chemgroups(44)%name="WDEP_ASH_F"
  chemgroups(44)%ptr=>WDEP_ASH_F_GROUP

  chemgroups(45)%name="PMCO"
  chemgroups(45)%ptr=>PMCO_GROUP

  chemgroups(46)%name="ALKENE"
  chemgroups(46)%ptr=>ALKENE_GROUP

  chemgroups(47)%name="WDEP_DUST_NAT_C"
  chemgroups(47)%ptr=>WDEP_DUST_NAT_C_GROUP

  chemgroups(48)%name="DDEP_DUST_ANT_F"
  chemgroups(48)%ptr=>DDEP_DUST_ANT_F_GROUP

  chemgroups(49)%name="WDEP_RDN"
  chemgroups(49)%ptr=>WDEP_RDN_GROUP

  chemgroups(50)%name="DDEP_SIA"
  chemgroups(50)%ptr=>DDEP_SIA_GROUP

  chemgroups(51)%name="WDEP_ASH_C"
  chemgroups(51)%ptr=>WDEP_ASH_C_GROUP

  chemgroups(52)%name="OXN"
  chemgroups(52)%ptr=>OXN_GROUP

  chemgroups(53)%name="TNO3"
  chemgroups(53)%ptr=>TNO3_GROUP

  chemgroups(54)%name="DDEP_DAOBS"
  chemgroups(54)%ptr=>DDEP_DAOBS_GROUP

  chemgroups(55)%name="DDEP_V1702A02B"
  chemgroups(55)%ptr=>DDEP_V1702A02B_GROUP

  chemgroups(56)%name="WDEP_DUST_ANT_F"
  chemgroups(56)%ptr=>WDEP_DUST_ANT_F_GROUP

  chemgroups(57)%name="SIA"
  chemgroups(57)%ptr=>SIA_GROUP

  chemgroups(58)%name="DDEP_ASH"
  chemgroups(58)%ptr=>DDEP_ASH_GROUP

  chemgroups(59)%name="WDEP_PMFINE"
  chemgroups(59)%ptr=>WDEP_PMFINE_GROUP

  chemgroups(60)%name="DDEP_TNO3"
  chemgroups(60)%ptr=>DDEP_TNO3_GROUP

  chemgroups(61)%name="DDEP_RDN"
  chemgroups(61)%ptr=>DDEP_RDN_GROUP

  chemgroups(62)%name="DDEP_ASH_F"
  chemgroups(62)%ptr=>DDEP_ASH_F_GROUP

  chemgroups(63)%name="DDEP_PMCO"
  chemgroups(63)%ptr=>DDEP_PMCO_GROUP

  chemgroups(64)%name="RDN"
  chemgroups(64)%ptr=>RDN_GROUP

  chemgroups(65)%name="WDEP_ASH"
  chemgroups(65)%ptr=>WDEP_ASH_GROUP

  chemgroups(66)%name="DAOBS"
  chemgroups(66)%ptr=>DAOBS_GROUP

  chemgroups(67)%name="WDEP_TNO3"
  chemgroups(67)%ptr=>WDEP_TNO3_GROUP

endsubroutine Init_ChemGroups
 !-----------------------------------------------------------
endmodule ChemGroups_ml
 !-----------------------------------------------------------
