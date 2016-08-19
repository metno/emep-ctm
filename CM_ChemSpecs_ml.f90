!>_________________________________________________________<

  module  ChemSpecs_adv_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 67 
 


   integer, public, parameter ::   & 
     IXADV_O3          =   1   &
  ,  IXADV_NO          =   2   &
  ,  IXADV_NO2         =   3   &
  ,  IXADV_PAN         =   4   &
  ,  IXADV_MPAN        =   5   &
  ,  IXADV_NO3         =   6   &
  ,  IXADV_N2O5        =   7   &
  ,  IXADV_ISONO3      =   8   &
  ,  IXADV_HNO3        =   9

   integer, public, parameter ::   & 
     IXADV_HONO        =  10   &
  ,  IXADV_CH3COO2     =  11   &
  ,  IXADV_MACR        =  12   &
  ,  IXADV_ISNI        =  13   &
  ,  IXADV_ISNIR       =  14   &
  ,  IXADV_GLYOX       =  15   &
  ,  IXADV_MGLYOX      =  16   &
  ,  IXADV_MAL         =  17   &
  ,  IXADV_MEK         =  18   &
  ,  IXADV_MVK         =  19

   integer, public, parameter ::   & 
     IXADV_HCHO        =  20   &
  ,  IXADV_CH3CHO      =  21   &
  ,  IXADV_C2H6        =  22   &
  ,  IXADV_NC4H10      =  23   &
  ,  IXADV_C2H4        =  24   &
  ,  IXADV_C3H6        =  25   &
  ,  IXADV_OXYL        =  26   &
  ,  IXADV_C5H8        =  27   &
  ,  IXADV_APINENE     =  28   &
  ,  IXADV_CH3O2H      =  29

   integer, public, parameter ::   & 
     IXADV_C2H5OOH     =  30   &
  ,  IXADV_BURO2H      =  31   &
  ,  IXADV_ETRO2H      =  32   &
  ,  IXADV_PRRO2H      =  33   &
  ,  IXADV_OXYO2H      =  34   &
  ,  IXADV_MEKO2H      =  35   &
  ,  IXADV_MALO2H      =  36   &
  ,  IXADV_MVKO2H      =  37   &
  ,  IXADV_MACROOH     =  38   &
  ,  IXADV_MACO3H      =  39

   integer, public, parameter ::   & 
     IXADV_MACO2H      =  40   &
  ,  IXADV_ISRO2H      =  41   &
  ,  IXADV_H2O2        =  42   &
  ,  IXADV_CH3COO2H    =  43   &
  ,  IXADV_ISONO3H     =  44   &
  ,  IXADV_ISNIRH      =  45   &
  ,  IXADV_CH3OH       =  46   &
  ,  IXADV_C2H5OH      =  47   &
  ,  IXADV_ACETOL      =  48   &
  ,  IXADV_H2          =  49

   integer, public, parameter ::   & 
     IXADV_CO          =  50   &
  ,  IXADV_CH4         =  51   &
  ,  IXADV_SO2         =  52   &
  ,  IXADV_SO4         =  53   &
  ,  IXADV_NH3         =  54   &
  ,  IXADV_NO3_F       =  55   &
  ,  IXADV_NO3_C       =  56   &
  ,  IXADV_NH4_F       =  57   &
  ,  IXADV_PPM25       =  58   &
  ,  IXADV_PPM25_FIRE  =  59

   integer, public, parameter ::   & 
     IXADV_PPM_C       =  60   &
  ,  IXADV_SEASALT_F   =  61   &
  ,  IXADV_SEASALT_C   =  62   &
  ,  IXADV_SEASALT_G   =  63   &
  ,  IXADV_DUST_NAT_F  =  64   &
  ,  IXADV_DUST_NAT_C  =  65   &
  ,  IXADV_RN222       =  66   &
  ,  IXADV_PB210       =  67

 !-----------------------------------------------------------
  end module ChemSpecs_adv_ml
!>_________________________________________________________<

  module  ChemSpecs_shl_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 16 
 


   integer, public, parameter ::   & 
     IXSHL_OD          =   1   &
  ,  IXSHL_OP          =   2   &
  ,  IXSHL_OH          =   3   &
  ,  IXSHL_HO2         =   4   &
  ,  IXSHL_CH3O2       =   5   &
  ,  IXSHL_C2H5O2      =   6   &
  ,  IXSHL_SECC4H9O2   =   7   &
  ,  IXSHL_ISRO2       =   8   &
  ,  IXSHL_ETRO2       =   9

   integer, public, parameter ::   & 
     IXSHL_PRRO2       =  10   &
  ,  IXSHL_OXYO2       =  11   &
  ,  IXSHL_MEKO2       =  12   &
  ,  IXSHL_MALO2       =  13   &
  ,  IXSHL_MVKO2       =  14   &
  ,  IXSHL_MACRO2      =  15   &
  ,  IXSHL_MACO3       =  16

 !-----------------------------------------------------------
  end module ChemSpecs_shl_ml
!>_________________________________________________________<

  module  ChemSpecs_tot_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_TOT = 83 
 
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
  ,  CH3O2       =   5   &
  ,  C2H5O2      =   6   &
  ,  SECC4H9O2   =   7   &
  ,  ISRO2       =   8   &
  ,  ETRO2       =   9

   integer, public, parameter ::   & 
     PRRO2       =  10   &
  ,  OXYO2       =  11   &
  ,  MEKO2       =  12   &
  ,  MALO2       =  13   &
  ,  MVKO2       =  14   &
  ,  MACRO2      =  15   &
  ,  MACO3       =  16   &
  ,  O3          =  17   &
  ,  NO          =  18   &
  ,  NO2         =  19

   integer, public, parameter ::   & 
     PAN         =  20   &
  ,  MPAN        =  21   &
  ,  NO3         =  22   &
  ,  N2O5        =  23   &
  ,  ISONO3      =  24   &
  ,  HNO3        =  25   &
  ,  HONO        =  26   &
  ,  CH3COO2     =  27   &
  ,  MACR        =  28   &
  ,  ISNI        =  29

   integer, public, parameter ::   & 
     ISNIR       =  30   &
  ,  GLYOX       =  31   &
  ,  MGLYOX      =  32   &
  ,  MAL         =  33   &
  ,  MEK         =  34   &
  ,  MVK         =  35   &
  ,  HCHO        =  36   &
  ,  CH3CHO      =  37   &
  ,  C2H6        =  38   &
  ,  NC4H10      =  39

   integer, public, parameter ::   & 
     C2H4        =  40   &
  ,  C3H6        =  41   &
  ,  OXYL        =  42   &
  ,  C5H8        =  43   &
  ,  APINENE     =  44   &
  ,  CH3O2H      =  45   &
  ,  C2H5OOH     =  46   &
  ,  BURO2H      =  47   &
  ,  ETRO2H      =  48   &
  ,  PRRO2H      =  49

   integer, public, parameter ::   & 
     OXYO2H      =  50   &
  ,  MEKO2H      =  51   &
  ,  MALO2H      =  52   &
  ,  MVKO2H      =  53   &
  ,  MACROOH     =  54   &
  ,  MACO3H      =  55   &
  ,  MACO2H      =  56   &
  ,  ISRO2H      =  57   &
  ,  H2O2        =  58   &
  ,  CH3COO2H    =  59

   integer, public, parameter ::   & 
     ISONO3H     =  60   &
  ,  ISNIRH      =  61   &
  ,  CH3OH       =  62   &
  ,  C2H5OH      =  63   &
  ,  ACETOL      =  64   &
  ,  H2          =  65   &
  ,  CO          =  66   &
  ,  CH4         =  67   &
  ,  SO2         =  68   &
  ,  SO4         =  69

   integer, public, parameter ::   & 
     NH3         =  70   &
  ,  NO3_F       =  71   &
  ,  NO3_C       =  72   &
  ,  NH4_F       =  73   &
  ,  PPM25       =  74   &
  ,  PPM25_FIRE  =  75   &
  ,  PPM_C       =  76   &
  ,  SEASALT_F   =  77   &
  ,  SEASALT_C   =  78   &
  ,  SEASALT_G   =  79

   integer, public, parameter ::   & 
     DUST_NAT_F  =  80   &
  ,  DUST_NAT_C  =  81   &
  ,  RN222       =  82   &
  ,  PB210       =  83

 !-----------------------------------------------------------
  end module ChemSpecs_tot_ml
!>_________________________________________________________<

  module  ChemChemicals_ml
!-----------------------------------------------------------

   use ChemSpecs_tot_ml  ! => NSPEC_TOT, species indices
  implicit none
  private

  !/--   Characteristics of species: 
  !/--   Number, name, molwt, carbon num, nmhc (1) or not(0)
 
  public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc

  type, public :: Chemical 
       character(len=20) :: name
       integer           :: molwt
       integer           :: nmhc      ! nmhc (1) or not(0)
       integer           :: carbons   ! Carbon-number
       real              :: nitrogens ! Nitrogen-number
       integer           :: sulphurs  ! Sulphur-number
       real              :: ExtC      ! Extinction coef (aerosols)
       real              :: CiStar     ! VBS param
       real              :: DeltaH    ! VBS param
  end type Chemical
  type(Chemical), public, dimension(NSPEC_TOT) :: species

  contains
    subroutine define_chemicals()
    !+
    ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
    ! array, using indices from total list of species (advected + short-lived).
    !                                           MW  NM   C    N   S  ExtC C*  dH
     species(OD) = Chemical("OD          ",  16,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(OP) = Chemical("OP          ",  16,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(OH) = Chemical("OH          ",  17,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(HO2) = Chemical("HO2         ",  33,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(CH3O2) = Chemical("CH3O2       ",  47,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(C2H5O2) = Chemical("C2H5O2      ",  61,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(SECC4H9O2) = Chemical("SECC4H9O2   ",  89,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ISRO2) = Chemical("ISRO2       ", 101,  0,  5,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ETRO2) = Chemical("ETRO2       ",  77,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(PRRO2) = Chemical("PRRO2       ",  91,  0,  3,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(OXYO2) = Chemical("OXYO2       ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MEKO2) = Chemical("MEKO2       ", 103,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MALO2) = Chemical("MALO2       ", 147,  0,  5,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MVKO2) = Chemical("MVKO2       ", 119,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MACRO2) = Chemical("MACRO2      ", 119,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MACO3) = Chemical("MACO3       ", 101,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(O3) = Chemical("O3          ",  48,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(NO) = Chemical("NO          ",  30,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(NO2) = Chemical("NO2         ",  46,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(PAN) = Chemical("PAN         ", 121,  0,  2,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(MPAN) = Chemical("MPAN        ", 132,  0,  4,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(NO3) = Chemical("NO3         ",  62,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(N2O5) = Chemical("N2O5        ", 108,  0,  0,   2,  0,  0.0,  0.0000,    0.0 ) 
     species(ISONO3) = Chemical("ISONO3      ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(HNO3) = Chemical("HNO3        ",  63,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(HONO) = Chemical("HONO        ",  47,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(CH3COO2) = Chemical("CH3COO2     ",  75,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MACR) = Chemical("MACR        ",  70,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ISNI) = Chemical("ISNI        ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ISNIR) = Chemical("ISNIR       ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(GLYOX) = Chemical("GLYOX       ",  58,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MGLYOX) = Chemical("MGLYOX      ",  72,  0,  3,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MAL) = Chemical("MAL         ",  98,  0,  5,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MEK) = Chemical("MEK         ",  72,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MVK) = Chemical("MVK         ",  70,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(HCHO) = Chemical("HCHO        ",  30,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(CH3CHO) = Chemical("CH3CHO      ",  44,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(C2H6) = Chemical("C2H6        ",  30,  1,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(NC4H10) = Chemical("NC4H10      ",  58,  1,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(C2H4) = Chemical("C2H4        ",  28,  1,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(C3H6) = Chemical("C3H6        ",  42,  1,  3,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(OXYL) = Chemical("OXYL        ", 106,  1,  8,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(C5H8) = Chemical("C5H8        ",  68,  1,  5,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(APINENE) = Chemical("APINENE     ", 136,  1, 10,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(CH3O2H) = Chemical("CH3O2H      ",  48,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(C2H5OOH) = Chemical("C2H5OOH     ",  62,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(BURO2H) = Chemical("BURO2H      ",  90,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ETRO2H) = Chemical("ETRO2H      ",  78,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(PRRO2H) = Chemical("PRRO2H      ",  92,  0,  3,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(OXYO2H) = Chemical("OXYO2H      ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MEKO2H) = Chemical("MEKO2H      ", 104,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MALO2H) = Chemical("MALO2H      ", 147,  0,  5,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MVKO2H) = Chemical("MVKO2H      ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MACROOH) = Chemical("MACROOH     ", 120,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MACO3H) = Chemical("MACO3H      ", 102,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(MACO2H) = Chemical("MACO2H      ",  86,  0,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ISRO2H) = Chemical("ISRO2H      ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(H2O2) = Chemical("H2O2        ",  34,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(CH3COO2H) = Chemical("CH3COO2H    ",  76,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ISONO3H) = Chemical("ISONO3H     ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ISNIRH) = Chemical("ISNIRH      ",   1,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(CH3OH) = Chemical("CH3OH       ",  32,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(C2H5OH) = Chemical("C2H5OH      ",  46,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(ACETOL) = Chemical("ACETOL      ",  74,  0,  3,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(H2) = Chemical("H2          ",   2,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(CO) = Chemical("CO          ",  28,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(CH4) = Chemical("CH4         ",  16,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(SO2) = Chemical("SO2         ",  64,  0,  0,   0,  1,  0.0,  0.0000,    0.0 ) 
     species(SO4) = Chemical("SO4         ",  96,  0,  0,   0,  1,  8.5,  0.0000,    0.0 ) 
     species(NH3) = Chemical("NH3         ",  17,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(NO3_F) = Chemical("NO3_F       ",  62,  0,  0,   1,  0,  8.5,  0.0000,    0.0 ) 
     species(NO3_C) = Chemical("NO3_C       ",  62,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
     species(NH4_F) = Chemical("NH4_F       ",  18,  0,  0,   1,  0,  8.5,  0.0000,    0.0 ) 
     species(PPM25) = Chemical("PPM25       ",  12,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(PPM25_FIRE) = Chemical("PPM25_FIRE  ",  12,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(PPM_C) = Chemical("PPM_C       ",  12,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(SEASALT_F) = Chemical("SEASALT_F   ",  58,  0,  0,   0,  0,  3.0,  0.0000,    0.0 ) 
     species(SEASALT_C) = Chemical("SEASALT_C   ",  58,  0,  0,   0,  0,  0.4,  0.0000,    0.0 ) 
     species(SEASALT_G) = Chemical("SEASALT_G   ",  58,  0,  0,   0,  0,  0.4,  0.0000,    0.0 ) 
     species(DUST_NAT_F) = Chemical("DUST_NAT_F  ", 200,  0,  0,   0,  0,  1.0,  0.0000,    0.0 ) 
     species(DUST_NAT_C) = Chemical("DUST_NAT_C  ", 200,  0,  0,   0,  0,  0.3,  0.0000,    0.0 ) 
     species(RN222) = Chemical("RN222       ", 222,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
     species(PB210) = Chemical("PB210       ", 210,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
   end subroutine define_chemicals
 end module ChemChemicals_ml
 !-----------------------------------------------------------
