! <CM_ChemSpecs_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_10(3282)>
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

module ChemSpecs_adv_ml
!-----------------------------------------------------------


implicit none
!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 111 
 


   integer, public, parameter ::   & 
     IXADV_O3          =   1   &
  ,  IXADV_NO          =   2   &
  ,  IXADV_NO2         =   3   &
  ,  IXADV_SHIPNOX     =   4   &
  ,  IXADV_PAN         =   5   &
  ,  IXADV_MPAN        =   6   &
  ,  IXADV_NO3         =   7   &
  ,  IXADV_N2O5        =   8   &
  ,  IXADV_ISONO3      =   9

   integer, public, parameter ::   & 
     IXADV_HNO3        =  10   &
  ,  IXADV_HONO        =  11   &
  ,  IXADV_CH3COO2     =  12   &
  ,  IXADV_MACR        =  13   &
  ,  IXADV_ISNI        =  14   &
  ,  IXADV_ISNIR       =  15   &
  ,  IXADV_GLYOX       =  16   &
  ,  IXADV_MGLYOX      =  17   &
  ,  IXADV_MAL         =  18   &
  ,  IXADV_MEK         =  19

   integer, public, parameter ::   & 
     IXADV_MVK         =  20   &
  ,  IXADV_HCHO        =  21   &
  ,  IXADV_CH3CHO      =  22   &
  ,  IXADV_C2H6        =  23   &
  ,  IXADV_NC4H10      =  24   &
  ,  IXADV_C2H4        =  25   &
  ,  IXADV_C3H6        =  26   &
  ,  IXADV_OXYL        =  27   &
  ,  IXADV_C5H8        =  28   &
  ,  IXADV_APINENE     =  29

   integer, public, parameter ::   & 
     IXADV_CH3O2H      =  30   &
  ,  IXADV_C2H5OOH     =  31   &
  ,  IXADV_BURO2H      =  32   &
  ,  IXADV_ETRO2H      =  33   &
  ,  IXADV_PRRO2H      =  34   &
  ,  IXADV_OXYO2H      =  35   &
  ,  IXADV_MEKO2H      =  36   &
  ,  IXADV_MALO2H      =  37   &
  ,  IXADV_MVKO2H      =  38   &
  ,  IXADV_MACROOH     =  39

   integer, public, parameter ::   & 
     IXADV_MACO3H      =  40   &
  ,  IXADV_MACO2H      =  41   &
  ,  IXADV_ISRO2H      =  42   &
  ,  IXADV_H2O2        =  43   &
  ,  IXADV_CH3COO2H    =  44   &
  ,  IXADV_ISONO3H     =  45   &
  ,  IXADV_ISNIRH      =  46   &
  ,  IXADV_CH3OH       =  47   &
  ,  IXADV_C2H5OH      =  48   &
  ,  IXADV_ACETOL      =  49

   integer, public, parameter ::   & 
     IXADV_H2          =  50   &
  ,  IXADV_CO          =  51   &
  ,  IXADV_CH4         =  52   &
  ,  IXADV_SO2         =  53   &
  ,  IXADV_SO4         =  54   &
  ,  IXADV_NH3         =  55   &
  ,  IXADV_NO3_F       =  56   &
  ,  IXADV_NO3_C       =  57   &
  ,  IXADV_NH4_F       =  58   &
  ,  IXADV_DUMMY       =  59

   integer, public, parameter ::   & 
     IXADV_ASH_F       =  60   &
  ,  IXADV_ASH_C       =  61   &
  ,  IXADV_POM_F_WOOD  =  62   &
  ,  IXADV_POM_F_FFUEL =  63   &
  ,  IXADV_POM_C_FFUEL =  64   &
  ,  IXADV_EC_F_WOOD_NEW=  65   &
  ,  IXADV_EC_F_WOOD_AGE=  66   &
  ,  IXADV_EC_C_WOOD   =  67   &
  ,  IXADV_EC_F_FFUEL_NEW=  68   &
  ,  IXADV_EC_F_FFUEL_AGE=  69

   integer, public, parameter ::   & 
     IXADV_EC_C_FFUEL  =  70   &
  ,  IXADV_REMPPM25    =  71   &
  ,  IXADV_REMPPM_C    =  72   &
  ,  IXADV_FFIRE_OM    =  73   &
  ,  IXADV_FFIRE_BC    =  74   &
  ,  IXADV_FFIRE_REMPPM25=  75   &
  ,  IXADV_OM25_BGND   =  76   &
  ,  IXADV_OM25_P      =  77   &
  ,  IXADV_ASOC_NG100  =  78   &
  ,  IXADV_ASOC_UG1    =  79

   integer, public, parameter ::   & 
     IXADV_ASOC_UG10   =  80   &
  ,  IXADV_ASOC_UG1E2  =  81   &
  ,  IXADV_ASOC_UG1E3  =  82   &
  ,  IXADV_NON_C_ASOA_NG100=  83   &
  ,  IXADV_NON_C_ASOA_UG1=  84   &
  ,  IXADV_NON_C_ASOA_UG10=  85   &
  ,  IXADV_NON_C_ASOA_UG1E2=  86   &
  ,  IXADV_NON_C_ASOA_UG1E3=  87   &
  ,  IXADV_BSOC_NG100  =  88   &
  ,  IXADV_BSOC_UG1    =  89

   integer, public, parameter ::   & 
     IXADV_BSOC_UG10   =  90   &
  ,  IXADV_BSOC_UG1E2  =  91   &
  ,  IXADV_BSOC_UG1E3  =  92   &
  ,  IXADV_NON_C_BSOA_NG100=  93   &
  ,  IXADV_NON_C_BSOA_UG1=  94   &
  ,  IXADV_NON_C_BSOA_UG10=  95   &
  ,  IXADV_NON_C_BSOA_UG1E2=  96   &
  ,  IXADV_NON_C_BSOA_UG1E3=  97   &
  ,  IXADV_FFFUEL_NG10 =  98   &
  ,  IXADV_WOODOA_NG10 =  99

   integer, public, parameter ::   & 
     IXADV_FFIREOA_NG10= 100   &
  ,  IXADV_SEASALT_F   = 101   &
  ,  IXADV_SEASALT_C   = 102   &
  ,  IXADV_DUST_ROAD_F = 103   &
  ,  IXADV_DUST_ROAD_C = 104   &
  ,  IXADV_DUST_WB_F   = 105   &
  ,  IXADV_DUST_WB_C   = 106   &
  ,  IXADV_DUST_SAH_F  = 107   &
  ,  IXADV_DUST_SAH_C  = 108   &
  ,  IXADV_RN222       = 109

   integer, public, parameter ::   & 
     IXADV_RNWATER     = 110   &
  ,  IXADV_PB210       = 111

 !-----------------------------------------------------------
  end module ChemSpecs_adv_ml
!>_________________________________________________________<

module ChemSpecs_shl_ml
!-----------------------------------------------------------


implicit none
!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 17 
 


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
  ,  IXSHL_MACO3       =  16   &
  ,  IXSHL_TERPPEROXY  =  17

 !-----------------------------------------------------------
  end module ChemSpecs_shl_ml
!>_________________________________________________________<

module ChemSpecs_tot_ml
!-----------------------------------------------------------


implicit none
!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_TOT = 128 
 
  ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL=23,   &!   Number of aerosol species
                FIRST_SEMIVOL=95, &!   First aerosol species
                LAST_SEMIVOL=117     !   Last  aerosol species  



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
  ,  TERPPEROXY  =  17   &
  ,  O3          =  18   &
  ,  NO          =  19

   integer, public, parameter ::   & 
     NO2         =  20   &
  ,  SHIPNOX     =  21   &
  ,  PAN         =  22   &
  ,  MPAN        =  23   &
  ,  NO3         =  24   &
  ,  N2O5        =  25   &
  ,  ISONO3      =  26   &
  ,  HNO3        =  27   &
  ,  HONO        =  28   &
  ,  CH3COO2     =  29

   integer, public, parameter ::   & 
     MACR        =  30   &
  ,  ISNI        =  31   &
  ,  ISNIR       =  32   &
  ,  GLYOX       =  33   &
  ,  MGLYOX      =  34   &
  ,  MAL         =  35   &
  ,  MEK         =  36   &
  ,  MVK         =  37   &
  ,  HCHO        =  38   &
  ,  CH3CHO      =  39

   integer, public, parameter ::   & 
     C2H6        =  40   &
  ,  NC4H10      =  41   &
  ,  C2H4        =  42   &
  ,  C3H6        =  43   &
  ,  OXYL        =  44   &
  ,  C5H8        =  45   &
  ,  APINENE     =  46   &
  ,  CH3O2H      =  47   &
  ,  C2H5OOH     =  48   &
  ,  BURO2H      =  49

   integer, public, parameter ::   & 
     ETRO2H      =  50   &
  ,  PRRO2H      =  51   &
  ,  OXYO2H      =  52   &
  ,  MEKO2H      =  53   &
  ,  MALO2H      =  54   &
  ,  MVKO2H      =  55   &
  ,  MACROOH     =  56   &
  ,  MACO3H      =  57   &
  ,  MACO2H      =  58   &
  ,  ISRO2H      =  59

   integer, public, parameter ::   & 
     H2O2        =  60   &
  ,  CH3COO2H    =  61   &
  ,  ISONO3H     =  62   &
  ,  ISNIRH      =  63   &
  ,  CH3OH       =  64   &
  ,  C2H5OH      =  65   &
  ,  ACETOL      =  66   &
  ,  H2          =  67   &
  ,  CO          =  68   &
  ,  CH4         =  69

   integer, public, parameter ::   & 
     SO2         =  70   &
  ,  SO4         =  71   &
  ,  NH3         =  72   &
  ,  NO3_F       =  73   &
  ,  NO3_C       =  74   &
  ,  NH4_F       =  75   &
  ,  DUMMY       =  76   &
  ,  ASH_F       =  77   &
  ,  ASH_C       =  78   &
  ,  POM_F_WOOD  =  79

   integer, public, parameter ::   & 
     POM_F_FFUEL =  80   &
  ,  POM_C_FFUEL =  81   &
  ,  EC_F_WOOD_NEW=  82   &
  ,  EC_F_WOOD_AGE=  83   &
  ,  EC_C_WOOD   =  84   &
  ,  EC_F_FFUEL_NEW=  85   &
  ,  EC_F_FFUEL_AGE=  86   &
  ,  EC_C_FFUEL  =  87   &
  ,  REMPPM25    =  88   &
  ,  REMPPM_C    =  89

   integer, public, parameter ::   & 
     FFIRE_OM    =  90   &
  ,  FFIRE_BC    =  91   &
  ,  FFIRE_REMPPM25=  92   &
  ,  OM25_BGND   =  93   &
  ,  OM25_P      =  94   &
  ,  ASOC_NG100  =  95   &
  ,  ASOC_UG1    =  96   &
  ,  ASOC_UG10   =  97   &
  ,  ASOC_UG1E2  =  98   &
  ,  ASOC_UG1E3  =  99

   integer, public, parameter ::   & 
     NON_C_ASOA_NG100= 100   &
  ,  NON_C_ASOA_UG1= 101   &
  ,  NON_C_ASOA_UG10= 102   &
  ,  NON_C_ASOA_UG1E2= 103   &
  ,  NON_C_ASOA_UG1E3= 104   &
  ,  BSOC_NG100  = 105   &
  ,  BSOC_UG1    = 106   &
  ,  BSOC_UG10   = 107   &
  ,  BSOC_UG1E2  = 108   &
  ,  BSOC_UG1E3  = 109

   integer, public, parameter ::   & 
     NON_C_BSOA_NG100= 110   &
  ,  NON_C_BSOA_UG1= 111   &
  ,  NON_C_BSOA_UG10= 112   &
  ,  NON_C_BSOA_UG1E2= 113   &
  ,  NON_C_BSOA_UG1E3= 114   &
  ,  FFFUEL_NG10 = 115   &
  ,  WOODOA_NG10 = 116   &
  ,  FFIREOA_NG10= 117   &
  ,  SEASALT_F   = 118   &
  ,  SEASALT_C   = 119

   integer, public, parameter ::   & 
     DUST_ROAD_F = 120   &
  ,  DUST_ROAD_C = 121   &
  ,  DUST_WB_F   = 122   &
  ,  DUST_WB_C   = 123   &
  ,  DUST_SAH_F  = 124   &
  ,  DUST_SAH_C  = 125   &
  ,  RN222       = 126   &
  ,  RNWATER     = 127   &
  ,  PB210       = 128

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
    species(CH3O2       ) = Chemical("CH3O2       ",  47.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(C2H5O2      ) = Chemical("C2H5O2      ",  61.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(SECC4H9O2   ) = Chemical("SECC4H9O2   ",  89.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(ISRO2       ) = Chemical("ISRO2       ", 101.0000,  0,  5,   0,  0,  0.0000,    0.0 ) 
    species(ETRO2       ) = Chemical("ETRO2       ",  77.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(PRRO2       ) = Chemical("PRRO2       ",  91.0000,  0,  3,   0,  0,  0.0000,    0.0 ) 
    species(OXYO2       ) = Chemical("OXYO2       ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(MEKO2       ) = Chemical("MEKO2       ", 103.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MALO2       ) = Chemical("MALO2       ", 147.0000,  0,  5,   0,  0,  0.0000,    0.0 ) 
    species(MVKO2       ) = Chemical("MVKO2       ", 119.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MACRO2      ) = Chemical("MACRO2      ", 119.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MACO3       ) = Chemical("MACO3       ", 101.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(TERPPEROXY  ) = Chemical("TERPPEROXY  ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(O3          ) = Chemical("O3          ",  48.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(NO          ) = Chemical("NO          ",  30.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NO2         ) = Chemical("NO2         ",  46.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(SHIPNOX     ) = Chemical("SHIPNOX     ",  46.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(PAN         ) = Chemical("PAN         ", 121.0000,  0,  2,   1,  0,  0.0000,    0.0 ) 
    species(MPAN        ) = Chemical("MPAN        ", 132.0000,  0,  4,   1,  0,  0.0000,    0.0 ) 
    species(NO3         ) = Chemical("NO3         ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(N2O5        ) = Chemical("N2O5        ", 108.0000,  0,  0,   2,  0,  0.0000,    0.0 ) 
    species(ISONO3      ) = Chemical("ISONO3      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(HNO3        ) = Chemical("HNO3        ",  63.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(HONO        ) = Chemical("HONO        ",  47.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(CH3COO2     ) = Chemical("CH3COO2     ",  75.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(MACR        ) = Chemical("MACR        ",  70.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(ISNI        ) = Chemical("ISNI        ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(ISNIR       ) = Chemical("ISNIR       ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(GLYOX       ) = Chemical("GLYOX       ",  58.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(MGLYOX      ) = Chemical("MGLYOX      ",  72.0000,  0,  3,   0,  0,  0.0000,    0.0 ) 
    species(MAL         ) = Chemical("MAL         ",  98.0000,  0,  5,   0,  0,  0.0000,    0.0 ) 
    species(MEK         ) = Chemical("MEK         ",  72.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MVK         ) = Chemical("MVK         ",  70.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(HCHO        ) = Chemical("HCHO        ",  30.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(CH3CHO      ) = Chemical("CH3CHO      ",  44.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(C2H6        ) = Chemical("C2H6        ",  30.0000,  1,  2,   0,  0,  0.0000,    0.0 ) 
    species(NC4H10      ) = Chemical("NC4H10      ",  58.0000,  1,  4,   0,  0,  0.0000,    0.0 ) 
    species(C2H4        ) = Chemical("C2H4        ",  28.0000,  1,  2,   0,  0,  0.0000,    0.0 ) 
    species(C3H6        ) = Chemical("C3H6        ",  42.0000,  1,  3,   0,  0,  0.0000,    0.0 ) 
    species(OXYL        ) = Chemical("OXYL        ", 106.0000,  1,  8,   0,  0,  0.0000,    0.0 ) 
    species(C5H8        ) = Chemical("C5H8        ",  68.0000,  1,  5,   0,  0,  0.0000,    0.0 ) 
    species(APINENE     ) = Chemical("APINENE     ", 136.0000,  1, 10,   0,  0,  0.0000,    0.0 ) 
    species(CH3O2H      ) = Chemical("CH3O2H      ",  48.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(C2H5OOH     ) = Chemical("C2H5OOH     ",  62.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(BURO2H      ) = Chemical("BURO2H      ",  90.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(ETRO2H      ) = Chemical("ETRO2H      ",  78.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(PRRO2H      ) = Chemical("PRRO2H      ",  92.0000,  0,  3,   0,  0,  0.0000,    0.0 ) 
    species(OXYO2H      ) = Chemical("OXYO2H      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(MEKO2H      ) = Chemical("MEKO2H      ", 104.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MALO2H      ) = Chemical("MALO2H      ", 147.0000,  0,  5,   0,  0,  0.0000,    0.0 ) 
    species(MVKO2H      ) = Chemical("MVKO2H      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(MACROOH     ) = Chemical("MACROOH     ", 120.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MACO3H      ) = Chemical("MACO3H      ", 102.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MACO2H      ) = Chemical("MACO2H      ",  86.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(ISRO2H      ) = Chemical("ISRO2H      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(H2O2        ) = Chemical("H2O2        ",  34.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(CH3COO2H    ) = Chemical("CH3COO2H    ",  76.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(ISONO3H     ) = Chemical("ISONO3H     ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(ISNIRH      ) = Chemical("ISNIRH      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(CH3OH       ) = Chemical("CH3OH       ",  32.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(C2H5OH      ) = Chemical("C2H5OH      ",  46.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(ACETOL      ) = Chemical("ACETOL      ",  74.0000,  0,  3,   0,  0,  0.0000,    0.0 ) 
    species(H2          ) = Chemical("H2          ",   2.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(CO          ) = Chemical("CO          ",  28.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(CH4         ) = Chemical("CH4         ",  16.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(SO2         ) = Chemical("SO2         ",  64.0000,  0,  0,   0,  1,  0.0000,    0.0 ) 
    species(SO4         ) = Chemical("SO4         ",  96.0000,  0,  0,   0,  1,  0.0000,    0.0 ) 
    species(NH3         ) = Chemical("NH3         ",  17.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NO3_F       ) = Chemical("NO3_F       ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NO3_C       ) = Chemical("NO3_C       ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NH4_F       ) = Chemical("NH4_F       ",  18.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(DUMMY       ) = Chemical("DUMMY       ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(ASH_F       ) = Chemical("ASH_F       ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(ASH_C       ) = Chemical("ASH_C       ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(POM_F_WOOD  ) = Chemical("POM_F_WOOD  ",  20.4000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(POM_F_FFUEL ) = Chemical("POM_F_FFUEL ",  15.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(POM_C_FFUEL ) = Chemical("POM_C_FFUEL ",  15.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(EC_F_WOOD_NEW) = Chemical("EC_F_WOOD_NEW",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(EC_F_WOOD_AGE) = Chemical("EC_F_WOOD_AGE",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(EC_C_WOOD   ) = Chemical("EC_C_WOOD   ",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(EC_F_FFUEL_NEW) = Chemical("EC_F_FFUEL_NEW",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(EC_F_FFUEL_AGE) = Chemical("EC_F_FFUEL_AGE",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(EC_C_FFUEL  ) = Chemical("EC_C_FFUEL  ",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(REMPPM25    ) = Chemical("REMPPM25    ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(REMPPM_C    ) = Chemical("REMPPM_C    ",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(FFIRE_OM    ) = Chemical("FFIRE_OM    ",  20.4000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(FFIRE_BC    ) = Chemical("FFIRE_BC    ",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(FFIRE_REMPPM25) = Chemical("FFIRE_REMPPM25",  12.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(OM25_BGND   ) = Chemical("OM25_BGND   ",  24.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(OM25_P      ) = Chemical("OM25_P      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(ASOC_NG100  ) = Chemical("ASOC_NG100  ",  12.0000,  0,  1,   0,  0,  0.1000,   30.0 ) 
    species(ASOC_UG1    ) = Chemical("ASOC_UG1    ",  12.0000,  0,  1,   0,  0,  1.0000,   30.0 ) 
    species(ASOC_UG10   ) = Chemical("ASOC_UG10   ",  12.0000,  0,  1,   0,  0, 10.0000,   30.0 ) 
    species(ASOC_UG1E2  ) = Chemical("ASOC_UG1E2  ",  12.0000,  0,  1,   0,  0,100.0000,   30.0 ) 
    species(ASOC_UG1E3  ) = Chemical("ASOC_UG1E3  ",  12.0000,  0,  1,   0,  0,1000.0000,   30.0 ) 
    species(NON_C_ASOA_NG100) = Chemical("NON_C_ASOA_NG100",   1.0000,  0,  0,   0,  0,  0.1000,   30.0 ) 
    species(NON_C_ASOA_UG1) = Chemical("NON_C_ASOA_UG1",   1.0000,  0,  0,   0,  0,  1.0000,   30.0 ) 
    species(NON_C_ASOA_UG10) = Chemical("NON_C_ASOA_UG10",   1.0000,  0,  0,   0,  0, 10.0000,   30.0 ) 
    species(NON_C_ASOA_UG1E2) = Chemical("NON_C_ASOA_UG1E2",   1.0000,  0,  0,   0,  0,100.0000,   30.0 ) 
    species(NON_C_ASOA_UG1E3) = Chemical("NON_C_ASOA_UG1E3",   1.0000,  0,  0,   0,  0,1000.0000,   30.0 ) 
    species(BSOC_NG100  ) = Chemical("BSOC_NG100  ",  12.0000,  0,  1,   0,  0,  0.1000,   30.0 ) 
    species(BSOC_UG1    ) = Chemical("BSOC_UG1    ",  12.0000,  0,  1,   0,  0,  1.0000,   30.0 ) 
    species(BSOC_UG10   ) = Chemical("BSOC_UG10   ",  12.0000,  0,  1,   0,  0, 10.0000,   30.0 ) 
    species(BSOC_UG1E2  ) = Chemical("BSOC_UG1E2  ",  12.0000,  0,  1,   0,  0,100.0000,   30.0 ) 
    species(BSOC_UG1E3  ) = Chemical("BSOC_UG1E3  ",  12.0000,  0,  1,   0,  0,1000.0000,   30.0 ) 
    species(NON_C_BSOA_NG100) = Chemical("NON_C_BSOA_NG100",   1.0000,  0,  0,   0,  0,  0.1000,   30.0 ) 
    species(NON_C_BSOA_UG1) = Chemical("NON_C_BSOA_UG1",   1.0000,  0,  0,   0,  0,  1.0000,   30.0 ) 
    species(NON_C_BSOA_UG10) = Chemical("NON_C_BSOA_UG10",   1.0000,  0,  0,   0,  0, 10.0000,   30.0 ) 
    species(NON_C_BSOA_UG1E2) = Chemical("NON_C_BSOA_UG1E2",   1.0000,  0,  0,   0,  0,100.0000,   30.0 ) 
    species(NON_C_BSOA_UG1E3) = Chemical("NON_C_BSOA_UG1E3",   1.0000,  0,  0,   0,  0,1000.0000,   30.0 ) 
    species(FFFUEL_NG10 ) = Chemical("FFFUEL_NG10 ",  15.0000,  0,  1,   0,  0,  0.0100,  112.0 ) 
    species(WOODOA_NG10 ) = Chemical("WOODOA_NG10 ",  20.4000,  0,  1,   0,  0,  0.0100,  112.0 ) 
    species(FFIREOA_NG10) = Chemical("FFIREOA_NG10",  20.4000,  0,  1,   0,  0,  0.0100,  112.0 ) 
    species(SEASALT_F   ) = Chemical("SEASALT_F   ",  58.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(SEASALT_C   ) = Chemical("SEASALT_C   ",  58.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_ROAD_F ) = Chemical("DUST_ROAD_F ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_ROAD_C ) = Chemical("DUST_ROAD_C ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_WB_F   ) = Chemical("DUST_WB_F   ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_WB_C   ) = Chemical("DUST_WB_C   ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_SAH_F  ) = Chemical("DUST_SAH_F  ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(DUST_SAH_C  ) = Chemical("DUST_SAH_C  ", 200.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(RN222       ) = Chemical("RN222       ", 222.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(RNWATER     ) = Chemical("RNWATER     ", 222.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(PB210       ) = Chemical("PB210       ", 210.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
  end subroutine define_chemicals
end module ChemChemicals_ml
 !-----------------------------------------------------------
