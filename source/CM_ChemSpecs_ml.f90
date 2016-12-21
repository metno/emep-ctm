!>_________________________________________________________<

module ChemSpecs_adv_ml
!-----------------------------------------------------------


implicit none
!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 52 
 


   integer, public, parameter ::   & 
     IXADV_TRACER1     =   1   &
  ,  IXADV_TRACER2     =   2   &
  ,  IXADV_O3          =   3   &
  ,  IXADV_NO          =   4   &
  ,  IXADV_NO2         =   5   &
  ,  IXADV_SHIPNOX     =   6   &
  ,  IXADV_HNO3        =   7   &
  ,  IXADV_VOC         =   8   &
  ,  IXADV_CH3CHO      =   9

   integer, public, parameter ::   & 
     IXADV_PAN         =  10   &
  ,  IXADV_CH3COO2     =  11   &
  ,  IXADV_SO2         =  12   &
  ,  IXADV_SO4         =  13   &
  ,  IXADV_NH3         =  14   &
  ,  IXADV_NO3_F       =  15   &
  ,  IXADV_NO3_C       =  16   &
  ,  IXADV_NH4_F       =  17   &
  ,  IXADV_CO          =  18   &
  ,  IXADV_HCHO        =  19

   integer, public, parameter ::   & 
     IXADV_CH4         =  20   &
  ,  IXADV_H2          =  21   &
  ,  IXADV_C2H6        =  22   &
  ,  IXADV_NC4H10      =  23   &
  ,  IXADV_OXYL        =  24   &
  ,  IXADV_C5H8        =  25   &
  ,  IXADV_APINENE     =  26   &
  ,  IXADV_C3H6        =  27   &
  ,  IXADV_C2H4        =  28   &
  ,  IXADV_C2H5OH      =  29

   integer, public, parameter ::   & 
     IXADV_CH3OH       =  30   &
  ,  IXADV_PPM25_FIRE  =  31   &
  ,  IXADV_MEK         =  32   &
  ,  IXADV_GLYOX       =  33   &
  ,  IXADV_MGLYOX      =  34   &
  ,  IXADV_UNREAC      =  35   &
  ,  IXADV_V1702A02B_1 =  36   &
  ,  IXADV_V1702A02B_2 =  37   &
  ,  IXADV_V1702A02B_3 =  38   &
  ,  IXADV_V1702A02B_4 =  39

   integer, public, parameter ::   & 
     IXADV_V1702A02B_5 =  40   &
  ,  IXADV_V1702A02B_6 =  41   &
  ,  IXADV_V1702A02B_7 =  42   &
  ,  IXADV_V1702A02B_8 =  43   &
  ,  IXADV_V1702A02B_9 =  44   &
  ,  IXADV_SEASALT_F   =  45   &
  ,  IXADV_SEASALT_C   =  46   &
  ,  IXADV_DUST_ROAD_F =  47   &
  ,  IXADV_DUST_ROAD_C =  48   &
  ,  IXADV_DUST_WB_F   =  49

   integer, public, parameter ::   & 
     IXADV_DUST_WB_C   =  50   &
  ,  IXADV_DUST_SAH_F  =  51   &
  ,  IXADV_DUST_SAH_C  =  52

 !-----------------------------------------------------------
  end module ChemSpecs_adv_ml
!>_________________________________________________________<

module ChemSpecs_shl_ml
!-----------------------------------------------------------


implicit none
!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 5 
 


   integer, public, parameter ::   & 
     IXSHL_OD          =   1   &
  ,  IXSHL_OP          =   2   &
  ,  IXSHL_OH          =   3   &
  ,  IXSHL_HO2         =   4   &
  ,  IXSHL_RO2         =   5

 !-----------------------------------------------------------
  end module ChemSpecs_shl_ml
!>_________________________________________________________<

module ChemSpecs_tot_ml
!-----------------------------------------------------------


implicit none
!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_TOT = 57 
 
  ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL=0,   &!   Number of aerosol species
                FIRST_SEMIVOL=-999, &!   First aerosol species
                LAST_SEMIVOL=-999     !   Last  aerosol species  



   integer, public, parameter ::   & 
     OD          =   1   &
  ,  OP          =   2   &
  ,  OH          =   3   &
  ,  HO2         =   4   &
  ,  RO2         =   5   &
  ,  TRACER1     =   6   &
  ,  TRACER2     =   7   &
  ,  O3          =   8   &
  ,  NO          =   9

   integer, public, parameter ::   & 
     NO2         =  10   &
  ,  SHIPNOX     =  11   &
  ,  HNO3        =  12   &
  ,  VOC         =  13   &
  ,  CH3CHO      =  14   &
  ,  PAN         =  15   &
  ,  CH3COO2     =  16   &
  ,  SO2         =  17   &
  ,  SO4         =  18   &
  ,  NH3         =  19

   integer, public, parameter ::   & 
     NO3_F       =  20   &
  ,  NO3_C       =  21   &
  ,  NH4_F       =  22   &
  ,  CO          =  23   &
  ,  HCHO        =  24   &
  ,  CH4         =  25   &
  ,  H2          =  26   &
  ,  C2H6        =  27   &
  ,  NC4H10      =  28   &
  ,  OXYL        =  29

   integer, public, parameter ::   & 
     C5H8        =  30   &
  ,  APINENE     =  31   &
  ,  C3H6        =  32   &
  ,  C2H4        =  33   &
  ,  C2H5OH      =  34   &
  ,  CH3OH       =  35   &
  ,  PPM25_FIRE  =  36   &
  ,  MEK         =  37   &
  ,  GLYOX       =  38   &
  ,  MGLYOX      =  39

   integer, public, parameter ::   & 
     UNREAC      =  40   &
  ,  V1702A02B_1 =  41   &
  ,  V1702A02B_2 =  42   &
  ,  V1702A02B_3 =  43   &
  ,  V1702A02B_4 =  44   &
  ,  V1702A02B_5 =  45   &
  ,  V1702A02B_6 =  46   &
  ,  V1702A02B_7 =  47   &
  ,  V1702A02B_8 =  48   &
  ,  V1702A02B_9 =  49

   integer, public, parameter ::   & 
     SEASALT_F   =  50   &
  ,  SEASALT_C   =  51   &
  ,  DUST_ROAD_F =  52   &
  ,  DUST_ROAD_C =  53   &
  ,  DUST_WB_F   =  54   &
  ,  DUST_WB_C   =  55   &
  ,  DUST_SAH_F  =  56   &
  ,  DUST_SAH_C  =  57

 !-----------------------------------------------------------
  end module ChemSpecs_tot_ml
!>_________________________________________________________<

module ChemChemicals_ml
!-----------------------------------------------------------


use ChemSpecs_tot_ml  ! => NSPEC_TOT, species indices
use ChemSpecs_shl_ml, only: NSPEC_SHL
use ChemSpecs_adv_ml, only: NSPEC_ADV
implicit none
private
!/--   Characteristics of species:
!/--   Number, name, molwt, carbon num, nmhc (1) or not(0)

public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc

type, public :: Chemical
     character(len=20) :: name
     real              :: molwt
     integer           :: nmhc      ! nmhc (1) or not(0)
     integer           :: carbons   ! Carbon-number
     real              :: nitrogens ! Nitrogen-number
     integer           :: sulphurs  ! Sulphur-number
     real              :: CiStar    ! VBS param
     real              :: DeltaH    ! VBS param
endtype Chemical
type(Chemical), public, dimension(NSPEC_TOT), target :: species
type(Chemical), public, dimension(:), pointer :: &
  species_shl=>null(),&             ! => species(..short lived..)
  species_adv=>null()               ! => species(..advected..)

contains
subroutine define_chemicals()
!+
! Pointers to short lived and advected portions of species
!
  species_shl=>species(1:NSPEC_SHL)
  species_adv=>species(NSPEC_SHL+1:NSPEC_SHL+NSPEC_ADV)
!+
! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
! array, using indices from total list of species (advected + short-lived).
!                                           MW  NM   C    N   S  C*  dH
    species(OD          ) = Chemical("OD          ",  16.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(OP          ) = Chemical("OP          ",  16.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(OH          ) = Chemical("OH          ",  17.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(HO2         ) = Chemical("HO2         ",  33.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(RO2         ) = Chemical("RO2         ",  47.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(TRACER1     ) = Chemical("TRACER1     ",  14.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(TRACER2     ) = Chemical("TRACER2     ",  14.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(O3          ) = Chemical("O3          ",  48.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(NO          ) = Chemical("NO          ",  30.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NO2         ) = Chemical("NO2         ",  46.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(SHIPNOX     ) = Chemical("SHIPNOX     ",  46.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(HNO3        ) = Chemical("HNO3        ",  63.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(VOC         ) = Chemical("VOC         ",  58.0000,  1,  4,   0,  0,  0.0000,    0.0 ) 
    species(CH3CHO      ) = Chemical("CH3CHO      ",  44.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(PAN         ) = Chemical("PAN         ", 121.0000,  0,  2,   1,  0,  0.0000,    0.0 ) 
    species(CH3COO2     ) = Chemical("CH3COO2     ",  75.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(SO2         ) = Chemical("SO2         ",  64.0000,  0,  0,   0,  1,  0.0000,    0.0 ) 
    species(SO4         ) = Chemical("SO4         ",  96.0000,  0,  0,   0,  1,  0.0000,    0.0 ) 
    species(NH3         ) = Chemical("NH3         ",  17.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NO3_F       ) = Chemical("NO3_F       ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NO3_C       ) = Chemical("NO3_C       ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NH4_F       ) = Chemical("NH4_F       ",  18.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(CO          ) = Chemical("CO          ",  28.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(HCHO        ) = Chemical("HCHO        ",  30.0300,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(CH4         ) = Chemical("CH4         ",  16.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(H2          ) = Chemical("H2          ",   2.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(C2H6        ) = Chemical("C2H6        ",  30.0000,  1,  2,   0,  0,  0.0000,    0.0 ) 
    species(NC4H10      ) = Chemical("NC4H10      ",  58.0000,  1,  4,   0,  0,  0.0000,    0.0 ) 
    species(OXYL        ) = Chemical("OXYL        ", 106.1700,  1,  8,   0,  0,  0.0000,    0.0 ) 
    species(C5H8        ) = Chemical("C5H8        ",  68.0000,  1,  5,   0,  0,  0.0000,    0.0 ) 
    species(APINENE     ) = Chemical("APINENE     ", 136.0000,  1, 10,   0,  0,  0.0000,    0.0 ) 
    species(C3H6        ) = Chemical("C3H6        ",  42.0800,  1,  3,   0,  0,  0.0000,    0.0 ) 
    species(C2H4        ) = Chemical("C2H4        ",  28.0000,  1,  2,   0,  0,  0.0000,    0.0 ) 
    species(C2H5OH      ) = Chemical("C2H5OH      ",  46.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(CH3OH       ) = Chemical("CH3OH       ",  32.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PPM25_FIRE  ) = Chemical("PPM25_FIRE  ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(MEK         ) = Chemical("MEK         ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(GLYOX       ) = Chemical("GLYOX       ",  16.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(MGLYOX      ) = Chemical("MGLYOX      ",  16.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(UNREAC      ) = Chemical("UNREAC      ",  26.0000,  0,  1,   1,  0,  0.0000,    0.0 ) 
    species(V1702A02B_1 ) = Chemical("V1702A02B_1 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(V1702A02B_2 ) = Chemical("V1702A02B_2 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(V1702A02B_3 ) = Chemical("V1702A02B_3 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(V1702A02B_4 ) = Chemical("V1702A02B_4 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(V1702A02B_5 ) = Chemical("V1702A02B_5 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(V1702A02B_6 ) = Chemical("V1702A02B_6 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(V1702A02B_7 ) = Chemical("V1702A02B_7 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(V1702A02B_8 ) = Chemical("V1702A02B_8 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(V1702A02B_9 ) = Chemical("V1702A02B_9 ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(SEASALT_F   ) = Chemical("SEASALT_F   ",  58.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(SEASALT_C   ) = Chemical("SEASALT_C   ",  58.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_ROAD_F ) = Chemical("DUST_ROAD_F ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_ROAD_C ) = Chemical("DUST_ROAD_C ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_WB_F   ) = Chemical("DUST_WB_F   ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_WB_C   ) = Chemical("DUST_WB_C   ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_SAH_F  ) = Chemical("DUST_SAH_F  ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_SAH_C  ) = Chemical("DUST_SAH_C  ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
  end subroutine define_chemicals
end module ChemChemicals_ml
 !-----------------------------------------------------------
