! <CM_ChemSpecs_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_5(2809)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-201409 met.no
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

   integer, public, parameter ::  NSPEC_ADV = 125 
 


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
  ,  IXADV_GAS_ASOA_OC =  58   &
  ,  IXADV_PART_ASOA_OC=  59

   integer, public, parameter ::   & 
     IXADV_PART_ASOA_OM=  60   &
  ,  IXADV_GAS_BSOA_OC =  61   &
  ,  IXADV_PART_BSOA_OC=  62   &
  ,  IXADV_PART_BSOA_OM=  63   &
  ,  IXADV_PART_FFUELOA25_OC=  64   &
  ,  IXADV_PART_FFUELOA25_OM=  65   &
  ,  IXADV_PART_WOODOA25_OC=  66   &
  ,  IXADV_PART_WOODOA25_OM=  67   &
  ,  IXADV_PART_FFIREOA25_OC=  68   &
  ,  IXADV_PART_FFIREOA25_OM=  69

   integer, public, parameter ::   & 
     IXADV_PART_OC10   =  70   &
  ,  IXADV_PART_OC25   =  71   &
  ,  IXADV_NONVOL_FFUELOC25=  72   &
  ,  IXADV_NONV_FFUELOC_COARSE=  73   &
  ,  IXADV_NONVOL_WOODOC25=  74   &
  ,  IXADV_NONVOL_BGNDOC=  75   &
  ,  IXADV_NONVOL_FFIREOC25=  76   &
  ,  IXADV_PART_OM_F   =  77   &
  ,  IXADV_POM_F_WOOD  =  78   &
  ,  IXADV_POM_F_FFUEL =  79

   integer, public, parameter ::   & 
     IXADV_POM_C_FFUEL =  80   &
  ,  IXADV_EC_F_WOOD_NEW=  81   &
  ,  IXADV_EC_F_WOOD_AGE=  82   &
  ,  IXADV_EC_C_WOOD   =  83   &
  ,  IXADV_EC_F_FFUEL_NEW=  84   &
  ,  IXADV_EC_F_FFUEL_AGE=  85   &
  ,  IXADV_EC_C_FFUEL  =  86   &
  ,  IXADV_REMPPM25    =  87   &
  ,  IXADV_REMPPM_C    =  88   &
  ,  IXADV_FFIRE_OM    =  89

   integer, public, parameter ::   & 
     IXADV_FFIRE_BC    =  90   &
  ,  IXADV_FFIRE_REMPPM25=  91   &
  ,  IXADV_ASOC_NG100  =  92   &
  ,  IXADV_ASOC_UG1    =  93   &
  ,  IXADV_ASOC_UG10   =  94   &
  ,  IXADV_ASOC_UG1E2  =  95   &
  ,  IXADV_ASOC_UG1E3  =  96   &
  ,  IXADV_NON_C_ASOA_NG100=  97   &
  ,  IXADV_NON_C_ASOA_UG1=  98   &
  ,  IXADV_NON_C_ASOA_UG10=  99

   integer, public, parameter ::   & 
     IXADV_NON_C_ASOA_UG1E2= 100   &
  ,  IXADV_NON_C_ASOA_UG1E3= 101   &
  ,  IXADV_BSOC_NG100  = 102   &
  ,  IXADV_BSOC_UG1    = 103   &
  ,  IXADV_BSOC_UG10   = 104   &
  ,  IXADV_BSOC_UG1E2  = 105   &
  ,  IXADV_BSOC_UG1E3  = 106   &
  ,  IXADV_NON_C_BSOA_NG100= 107   &
  ,  IXADV_NON_C_BSOA_UG1= 108   &
  ,  IXADV_NON_C_BSOA_UG10= 109

   integer, public, parameter ::   & 
     IXADV_NON_C_BSOA_UG1E2= 110   &
  ,  IXADV_NON_C_BSOA_UG1E3= 111   &
  ,  IXADV_FFFUEL_NG10 = 112   &
  ,  IXADV_WOODOA_NG10 = 113   &
  ,  IXADV_FFIREOA_NG10= 114   &
  ,  IXADV_SEASALT_F   = 115   &
  ,  IXADV_SEASALT_C   = 116   &
  ,  IXADV_DUST_ROAD_F = 117   &
  ,  IXADV_DUST_ROAD_C = 118   &
  ,  IXADV_DUST_WB_F   = 119

   integer, public, parameter ::   & 
     IXADV_DUST_WB_C   = 120   &
  ,  IXADV_DUST_SAH_F  = 121   &
  ,  IXADV_DUST_SAH_C  = 122   &
  ,  IXADV_RN222       = 123   &
  ,  IXADV_RNWATER     = 124   &
  ,  IXADV_PB210       = 125

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

   integer, public, parameter ::  NSPEC_TOT = 142 
 
  ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL=23,   &!   Number of aerosol species
                FIRST_SEMIVOL=109, &!   First aerosol species
                LAST_SEMIVOL=131     !   Last  aerosol species  



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
  ,  PAN         =  21   &
  ,  MPAN        =  22   &
  ,  NO3         =  23   &
  ,  N2O5        =  24   &
  ,  ISONO3      =  25   &
  ,  HNO3        =  26   &
  ,  HONO        =  27   &
  ,  CH3COO2     =  28   &
  ,  MACR        =  29

   integer, public, parameter ::   & 
     ISNI        =  30   &
  ,  ISNIR       =  31   &
  ,  GLYOX       =  32   &
  ,  MGLYOX      =  33   &
  ,  MAL         =  34   &
  ,  MEK         =  35   &
  ,  MVK         =  36   &
  ,  HCHO        =  37   &
  ,  CH3CHO      =  38   &
  ,  C2H6        =  39

   integer, public, parameter ::   & 
     NC4H10      =  40   &
  ,  C2H4        =  41   &
  ,  C3H6        =  42   &
  ,  OXYL        =  43   &
  ,  C5H8        =  44   &
  ,  APINENE     =  45   &
  ,  CH3O2H      =  46   &
  ,  C2H5OOH     =  47   &
  ,  BURO2H      =  48   &
  ,  ETRO2H      =  49

   integer, public, parameter ::   & 
     PRRO2H      =  50   &
  ,  OXYO2H      =  51   &
  ,  MEKO2H      =  52   &
  ,  MALO2H      =  53   &
  ,  MVKO2H      =  54   &
  ,  MACROOH     =  55   &
  ,  MACO3H      =  56   &
  ,  MACO2H      =  57   &
  ,  ISRO2H      =  58   &
  ,  H2O2        =  59

   integer, public, parameter ::   & 
     CH3COO2H    =  60   &
  ,  ISONO3H     =  61   &
  ,  ISNIRH      =  62   &
  ,  CH3OH       =  63   &
  ,  C2H5OH      =  64   &
  ,  ACETOL      =  65   &
  ,  H2          =  66   &
  ,  CO          =  67   &
  ,  CH4         =  68   &
  ,  SO2         =  69

   integer, public, parameter ::   & 
     SO4         =  70   &
  ,  NH3         =  71   &
  ,  NO3_F       =  72   &
  ,  NO3_C       =  73   &
  ,  NH4_F       =  74   &
  ,  GAS_ASOA_OC =  75   &
  ,  PART_ASOA_OC=  76   &
  ,  PART_ASOA_OM=  77   &
  ,  GAS_BSOA_OC =  78   &
  ,  PART_BSOA_OC=  79

   integer, public, parameter ::   & 
     PART_BSOA_OM=  80   &
  ,  PART_FFUELOA25_OC=  81   &
  ,  PART_FFUELOA25_OM=  82   &
  ,  PART_WOODOA25_OC=  83   &
  ,  PART_WOODOA25_OM=  84   &
  ,  PART_FFIREOA25_OC=  85   &
  ,  PART_FFIREOA25_OM=  86   &
  ,  PART_OC10   =  87   &
  ,  PART_OC25   =  88   &
  ,  NONVOL_FFUELOC25=  89

   integer, public, parameter ::   & 
     NONV_FFUELOC_COARSE=  90   &
  ,  NONVOL_WOODOC25=  91   &
  ,  NONVOL_BGNDOC=  92   &
  ,  NONVOL_FFIREOC25=  93   &
  ,  PART_OM_F   =  94   &
  ,  POM_F_WOOD  =  95   &
  ,  POM_F_FFUEL =  96   &
  ,  POM_C_FFUEL =  97   &
  ,  EC_F_WOOD_NEW=  98   &
  ,  EC_F_WOOD_AGE=  99

   integer, public, parameter ::   & 
     EC_C_WOOD   = 100   &
  ,  EC_F_FFUEL_NEW= 101   &
  ,  EC_F_FFUEL_AGE= 102   &
  ,  EC_C_FFUEL  = 103   &
  ,  REMPPM25    = 104   &
  ,  REMPPM_C    = 105   &
  ,  FFIRE_OM    = 106   &
  ,  FFIRE_BC    = 107   &
  ,  FFIRE_REMPPM25= 108   &
  ,  ASOC_NG100  = 109

   integer, public, parameter ::   & 
     ASOC_UG1    = 110   &
  ,  ASOC_UG10   = 111   &
  ,  ASOC_UG1E2  = 112   &
  ,  ASOC_UG1E3  = 113   &
  ,  NON_C_ASOA_NG100= 114   &
  ,  NON_C_ASOA_UG1= 115   &
  ,  NON_C_ASOA_UG10= 116   &
  ,  NON_C_ASOA_UG1E2= 117   &
  ,  NON_C_ASOA_UG1E3= 118   &
  ,  BSOC_NG100  = 119

   integer, public, parameter ::   & 
     BSOC_UG1    = 120   &
  ,  BSOC_UG10   = 121   &
  ,  BSOC_UG1E2  = 122   &
  ,  BSOC_UG1E3  = 123   &
  ,  NON_C_BSOA_NG100= 124   &
  ,  NON_C_BSOA_UG1= 125   &
  ,  NON_C_BSOA_UG10= 126   &
  ,  NON_C_BSOA_UG1E2= 127   &
  ,  NON_C_BSOA_UG1E3= 128   &
  ,  FFFUEL_NG10 = 129

   integer, public, parameter ::   & 
     WOODOA_NG10 = 130   &
  ,  FFIREOA_NG10= 131   &
  ,  SEASALT_F   = 132   &
  ,  SEASALT_C   = 133   &
  ,  DUST_ROAD_F = 134   &
  ,  DUST_ROAD_C = 135   &
  ,  DUST_WB_F   = 136   &
  ,  DUST_WB_C   = 137   &
  ,  DUST_SAH_F  = 138   &
  ,  DUST_SAH_C  = 139

   integer, public, parameter ::   & 
     RN222       = 140   &
  ,  RNWATER     = 141   &
  ,  PB210       = 142

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
    species(GAS_ASOA_OC ) = Chemical("GAS_ASOA_OC ",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_ASOA_OC) = Chemical("PART_ASOA_OC",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_ASOA_OM) = Chemical("PART_ASOA_OM",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(GAS_BSOA_OC ) = Chemical("GAS_BSOA_OC ",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_BSOA_OC) = Chemical("PART_BSOA_OC",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_BSOA_OM) = Chemical("PART_BSOA_OM",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(PART_FFUELOA25_OC) = Chemical("PART_FFUELOA25_OC",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_FFUELOA25_OM) = Chemical("PART_FFUELOA25_OM",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(PART_WOODOA25_OC) = Chemical("PART_WOODOA25_OC",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_WOODOA25_OM) = Chemical("PART_WOODOA25_OM",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(PART_FFIREOA25_OC) = Chemical("PART_FFIREOA25_OC",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_FFIREOA25_OM) = Chemical("PART_FFIREOA25_OM",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(PART_OC10   ) = Chemical("PART_OC10   ",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_OC25   ) = Chemical("PART_OC25   ",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(NONVOL_FFUELOC25) = Chemical("NONVOL_FFUELOC25",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(NONV_FFUELOC_COARSE) = Chemical("NONV_FFUELOC_COARSE",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(NONVOL_WOODOC25) = Chemical("NONVOL_WOODOC25",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(NONVOL_BGNDOC) = Chemical("NONVOL_BGNDOC",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(NONVOL_FFIREOC25) = Chemical("NONVOL_FFIREOC25",  12.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(PART_OM_F   ) = Chemical("PART_OM_F   ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
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
