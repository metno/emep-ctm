! <CM_ChemSpecs_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.15>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2017 met.no
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

   integer, public, parameter ::  NSPEC_ADV = 118 
 


   integer, public, parameter ::   & 
     IXADV_O3          =   1   &
  ,  IXADV_NO          =   2   &
  ,  IXADV_NO2         =   3   &
  ,  IXADV_HO2NO2      =   4   &
  ,  IXADV_SHIPNOX     =   5   &
  ,  IXADV_PAN         =   6   &
  ,  IXADV_NO3         =   7   &
  ,  IXADV_N2O5        =   8   &
  ,  IXADV_HNO3        =   9

   integer, public, parameter ::   & 
     IXADV_HONO        =  10   &
  ,  IXADV_CH3COO2     =  11   &
  ,  IXADV_GLYOX       =  12   &
  ,  IXADV_MGLYOX      =  13   &
  ,  IXADV_MAL         =  14   &
  ,  IXADV_MEK         =  15   &
  ,  IXADV_HCHO        =  16   &
  ,  IXADV_CH3CHO      =  17   &
  ,  IXADV_C2H6        =  18   &
  ,  IXADV_NC4H10      =  19

   integer, public, parameter ::   & 
     IXADV_C2H4        =  20   &
  ,  IXADV_C3H6        =  21   &
  ,  IXADV_OXYL        =  22   &
  ,  IXADV_C5H8        =  23   &
  ,  IXADV_APINENE     =  24   &
  ,  IXADV_BPINENE     =  25   &
  ,  IXADV_XTERP       =  26   &
  ,  IXADV_BIOTERP     =  27   &
  ,  IXADV_CH3O2H      =  28   &
  ,  IXADV_C2H5OOH     =  29

   integer, public, parameter ::   & 
     IXADV_BURO2H      =  30   &
  ,  IXADV_ETRO2H      =  31   &
  ,  IXADV_PRRO2H      =  32   &
  ,  IXADV_OXYO2H      =  33   &
  ,  IXADV_MEKO2H      =  34   &
  ,  IXADV_MALO2H      =  35   &
  ,  IXADV_H2O2        =  36   &
  ,  IXADV_CH3COO2H    =  37   &
  ,  IXADV_CH3OH       =  38   &
  ,  IXADV_C2H5OH      =  39

   integer, public, parameter ::   & 
     IXADV_ACETOL      =  40   &
  ,  IXADV_H2          =  41   &
  ,  IXADV_CO          =  42   &
  ,  IXADV_CH4         =  43   &
  ,  IXADV_SO2         =  44   &
  ,  IXADV_ISO2        =  45   &
  ,  IXADV_MACRO2      =  46   &
  ,  IXADV_MACR        =  47   &
  ,  IXADV_MACROOH     =  48   &
  ,  IXADV_IEPOX       =  49

   integer, public, parameter ::   & 
     IXADV_HACET       =  50   &
  ,  IXADV_ISOOH       =  51   &
  ,  IXADV_ISON        =  52   &
  ,  IXADV_HCOOH       =  53   &
  ,  IXADV_MPAN        =  54   &
  ,  IXADV_NALD        =  55   &
  ,  IXADV_HPALD       =  56   &
  ,  IXADV_PACALD      =  57   &
  ,  IXADV_MVK         =  58   &
  ,  IXADV_TERPOOH     =  59

   integer, public, parameter ::   & 
     IXADV_SO4         =  60   &
  ,  IXADV_NH3         =  61   &
  ,  IXADV_NO3_F       =  62   &
  ,  IXADV_NO3_C       =  63   &
  ,  IXADV_NH4_F       =  64   &
  ,  IXADV_DUMMY       =  65   &
  ,  IXADV_ASH_F       =  66   &
  ,  IXADV_ASH_C       =  67   &
  ,  IXADV_POM_F_WOOD  =  68   &
  ,  IXADV_POM_F_FFUEL =  69

   integer, public, parameter ::   & 
     IXADV_POM_C_FFUEL =  70   &
  ,  IXADV_EC_F_WOOD_NEW=  71   &
  ,  IXADV_EC_F_WOOD_AGE=  72   &
  ,  IXADV_EC_C_WOOD   =  73   &
  ,  IXADV_EC_F_FFUEL_NEW=  74   &
  ,  IXADV_EC_F_FFUEL_AGE=  75   &
  ,  IXADV_EC_C_FFUEL  =  76   &
  ,  IXADV_REMPPM25    =  77   &
  ,  IXADV_REMPPM_C    =  78   &
  ,  IXADV_FFIRE_OM    =  79

   integer, public, parameter ::   & 
     IXADV_FFIRE_BC    =  80   &
  ,  IXADV_FFIRE_REMPPM25=  81   &
  ,  IXADV_OM25_BGND   =  82   &
  ,  IXADV_OM25_P      =  83   &
  ,  IXADV_SQT_SOA_NV  =  84   &
  ,  IXADV_ASOC_NG100  =  85   &
  ,  IXADV_ASOC_UG1    =  86   &
  ,  IXADV_ASOC_UG10   =  87   &
  ,  IXADV_ASOC_UG1E2  =  88   &
  ,  IXADV_ASOC_UG1E3  =  89

   integer, public, parameter ::   & 
     IXADV_NON_C_ASOA_NG100=  90   &
  ,  IXADV_NON_C_ASOA_UG1=  91   &
  ,  IXADV_NON_C_ASOA_UG10=  92   &
  ,  IXADV_NON_C_ASOA_UG1E2=  93   &
  ,  IXADV_NON_C_ASOA_UG1E3=  94   &
  ,  IXADV_BSOC_NG100  =  95   &
  ,  IXADV_BSOC_UG1    =  96   &
  ,  IXADV_BSOC_UG10   =  97   &
  ,  IXADV_BSOC_UG1E2  =  98   &
  ,  IXADV_BSOC_UG1E3  =  99

   integer, public, parameter ::   & 
     IXADV_NON_C_BSOA_NG100= 100   &
  ,  IXADV_NON_C_BSOA_UG1= 101   &
  ,  IXADV_NON_C_BSOA_UG10= 102   &
  ,  IXADV_NON_C_BSOA_UG1E2= 103   &
  ,  IXADV_NON_C_BSOA_UG1E3= 104   &
  ,  IXADV_FFFUEL_NG10 = 105   &
  ,  IXADV_WOODOA_NG10 = 106   &
  ,  IXADV_FFIREOA_NG10= 107   &
  ,  IXADV_SEASALT_F   = 108   &
  ,  IXADV_SEASALT_C   = 109

   integer, public, parameter ::   & 
     IXADV_DUST_ROAD_F = 110   &
  ,  IXADV_DUST_ROAD_C = 111   &
  ,  IXADV_DUST_WB_F   = 112   &
  ,  IXADV_DUST_WB_C   = 113   &
  ,  IXADV_DUST_SAH_F  = 114   &
  ,  IXADV_DUST_SAH_C  = 115   &
  ,  IXADV_RN222       = 116   &
  ,  IXADV_RNWATER     = 117   &
  ,  IXADV_PB210       = 118

 !-----------------------------------------------------------
  end module ChemSpecs_adv_ml
!>_________________________________________________________<

module ChemSpecs_shl_ml
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
  ,  IXSHL_TERPO2      =  15   &
  ,  IXSHL_XMTO3_RO2   =  16

 !-----------------------------------------------------------
  end module ChemSpecs_shl_ml
!>_________________________________________________________<

module ChemSpecs_tot_ml
!-----------------------------------------------------------


implicit none
!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_TOT = 134 
 
  ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL=23,   &!   Number of aerosol species
                FIRST_SEMIVOL=101, &!   First aerosol species
                LAST_SEMIVOL=123     !   Last  aerosol species  



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
  ,  TERPO2      =  15   &
  ,  XMTO3_RO2   =  16   &
  ,  O3          =  17   &
  ,  NO          =  18   &
  ,  NO2         =  19

   integer, public, parameter ::   & 
     HO2NO2      =  20   &
  ,  SHIPNOX     =  21   &
  ,  PAN         =  22   &
  ,  NO3         =  23   &
  ,  N2O5        =  24   &
  ,  HNO3        =  25   &
  ,  HONO        =  26   &
  ,  CH3COO2     =  27   &
  ,  GLYOX       =  28   &
  ,  MGLYOX      =  29

   integer, public, parameter ::   & 
     MAL         =  30   &
  ,  MEK         =  31   &
  ,  HCHO        =  32   &
  ,  CH3CHO      =  33   &
  ,  C2H6        =  34   &
  ,  NC4H10      =  35   &
  ,  C2H4        =  36   &
  ,  C3H6        =  37   &
  ,  OXYL        =  38   &
  ,  C5H8        =  39

   integer, public, parameter ::   & 
     APINENE     =  40   &
  ,  BPINENE     =  41   &
  ,  XTERP       =  42   &
  ,  BIOTERP     =  43   &
  ,  CH3O2H      =  44   &
  ,  C2H5OOH     =  45   &
  ,  BURO2H      =  46   &
  ,  ETRO2H      =  47   &
  ,  PRRO2H      =  48   &
  ,  OXYO2H      =  49

   integer, public, parameter ::   & 
     MEKO2H      =  50   &
  ,  MALO2H      =  51   &
  ,  H2O2        =  52   &
  ,  CH3COO2H    =  53   &
  ,  CH3OH       =  54   &
  ,  C2H5OH      =  55   &
  ,  ACETOL      =  56   &
  ,  H2          =  57   &
  ,  CO          =  58   &
  ,  CH4         =  59

   integer, public, parameter ::   & 
     SO2         =  60   &
  ,  ISO2        =  61   &
  ,  MACRO2      =  62   &
  ,  MACR        =  63   &
  ,  MACROOH     =  64   &
  ,  IEPOX       =  65   &
  ,  HACET       =  66   &
  ,  ISOOH       =  67   &
  ,  ISON        =  68   &
  ,  HCOOH       =  69

   integer, public, parameter ::   & 
     MPAN        =  70   &
  ,  NALD        =  71   &
  ,  HPALD       =  72   &
  ,  PACALD      =  73   &
  ,  MVK         =  74   &
  ,  TERPOOH     =  75   &
  ,  SO4         =  76   &
  ,  NH3         =  77   &
  ,  NO3_F       =  78   &
  ,  NO3_C       =  79

   integer, public, parameter ::   & 
     NH4_F       =  80   &
  ,  DUMMY       =  81   &
  ,  ASH_F       =  82   &
  ,  ASH_C       =  83   &
  ,  POM_F_WOOD  =  84   &
  ,  POM_F_FFUEL =  85   &
  ,  POM_C_FFUEL =  86   &
  ,  EC_F_WOOD_NEW=  87   &
  ,  EC_F_WOOD_AGE=  88   &
  ,  EC_C_WOOD   =  89

   integer, public, parameter ::   & 
     EC_F_FFUEL_NEW=  90   &
  ,  EC_F_FFUEL_AGE=  91   &
  ,  EC_C_FFUEL  =  92   &
  ,  REMPPM25    =  93   &
  ,  REMPPM_C    =  94   &
  ,  FFIRE_OM    =  95   &
  ,  FFIRE_BC    =  96   &
  ,  FFIRE_REMPPM25=  97   &
  ,  OM25_BGND   =  98   &
  ,  OM25_P      =  99

   integer, public, parameter ::   & 
     SQT_SOA_NV  = 100   &
  ,  ASOC_NG100  = 101   &
  ,  ASOC_UG1    = 102   &
  ,  ASOC_UG10   = 103   &
  ,  ASOC_UG1E2  = 104   &
  ,  ASOC_UG1E3  = 105   &
  ,  NON_C_ASOA_NG100= 106   &
  ,  NON_C_ASOA_UG1= 107   &
  ,  NON_C_ASOA_UG10= 108   &
  ,  NON_C_ASOA_UG1E2= 109

   integer, public, parameter ::   & 
     NON_C_ASOA_UG1E3= 110   &
  ,  BSOC_NG100  = 111   &
  ,  BSOC_UG1    = 112   &
  ,  BSOC_UG10   = 113   &
  ,  BSOC_UG1E2  = 114   &
  ,  BSOC_UG1E3  = 115   &
  ,  NON_C_BSOA_NG100= 116   &
  ,  NON_C_BSOA_UG1= 117   &
  ,  NON_C_BSOA_UG10= 118   &
  ,  NON_C_BSOA_UG1E2= 119

   integer, public, parameter ::   & 
     NON_C_BSOA_UG1E3= 120   &
  ,  FFFUEL_NG10 = 121   &
  ,  WOODOA_NG10 = 122   &
  ,  FFIREOA_NG10= 123   &
  ,  SEASALT_F   = 124   &
  ,  SEASALT_C   = 125   &
  ,  DUST_ROAD_F = 126   &
  ,  DUST_ROAD_C = 127   &
  ,  DUST_WB_F   = 128   &
  ,  DUST_WB_C   = 129

   integer, public, parameter ::   & 
     DUST_SAH_F  = 130   &
  ,  DUST_SAH_C  = 131   &
  ,  RN222       = 132   &
  ,  RNWATER     = 133   &
  ,  PB210       = 134

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
    species(TERPO2      ) = Chemical("TERPO2      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(XMTO3_RO2   ) = Chemical("XMTO3_RO2   ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(O3          ) = Chemical("O3          ",  48.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(NO          ) = Chemical("NO          ",  30.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(NO2         ) = Chemical("NO2         ",  46.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(HO2NO2      ) = Chemical("HO2NO2      ",  79.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(SHIPNOX     ) = Chemical("SHIPNOX     ",  46.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(PAN         ) = Chemical("PAN         ", 121.0000,  0,  2,   1,  0,  0.0000,    0.0 ) 
    species(NO3         ) = Chemical("NO3         ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(N2O5        ) = Chemical("N2O5        ", 108.0000,  0,  0,   2,  0,  0.0000,    0.0 ) 
    species(HNO3        ) = Chemical("HNO3        ",  63.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(HONO        ) = Chemical("HONO        ",  47.0000,  0,  0,   1,  0,  0.0000,    0.0 ) 
    species(CH3COO2     ) = Chemical("CH3COO2     ",  75.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(GLYOX       ) = Chemical("GLYOX       ",  58.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(MGLYOX      ) = Chemical("MGLYOX      ",  72.0000,  0,  3,   0,  0,  0.0000,    0.0 ) 
    species(MAL         ) = Chemical("MAL         ",  98.0000,  0,  5,   0,  0,  0.0000,    0.0 ) 
    species(MEK         ) = Chemical("MEK         ",  72.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(HCHO        ) = Chemical("HCHO        ",  30.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(CH3CHO      ) = Chemical("CH3CHO      ",  44.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(C2H6        ) = Chemical("C2H6        ",  30.0000,  1,  2,   0,  0,  0.0000,    0.0 ) 
    species(NC4H10      ) = Chemical("NC4H10      ",  58.0000,  1,  4,   0,  0,  0.0000,    0.0 ) 
    species(C2H4        ) = Chemical("C2H4        ",  28.0000,  1,  2,   0,  0,  0.0000,    0.0 ) 
    species(C3H6        ) = Chemical("C3H6        ",  42.0000,  1,  3,   0,  0,  0.0000,    0.0 ) 
    species(OXYL        ) = Chemical("OXYL        ", 106.0000,  1,  8,   0,  0,  0.0000,    0.0 ) 
    species(C5H8        ) = Chemical("C5H8        ",  68.0000,  1,  5,   0,  0,  0.0000,    0.0 ) 
    species(APINENE     ) = Chemical("APINENE     ", 136.0000,  1, 10,   0,  0,  0.0000,    0.0 ) 
    species(BPINENE     ) = Chemical("BPINENE     ", 136.0000,  1, 10,   0,  0,  0.0000,    0.0 ) 
    species(XTERP       ) = Chemical("XTERP       ", 136.0000,  1, 10,   0,  0,  0.0000,    0.0 ) 
    species(BIOTERP     ) = Chemical("BIOTERP     ", 136.0000,  1, 10,   0,  0,  0.0000,    0.0 ) 
    species(CH3O2H      ) = Chemical("CH3O2H      ",  48.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(C2H5OOH     ) = Chemical("C2H5OOH     ",  62.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(BURO2H      ) = Chemical("BURO2H      ",  90.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(ETRO2H      ) = Chemical("ETRO2H      ",  78.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(PRRO2H      ) = Chemical("PRRO2H      ",  92.0000,  0,  3,   0,  0,  0.0000,    0.0 ) 
    species(OXYO2H      ) = Chemical("OXYO2H      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(MEKO2H      ) = Chemical("MEKO2H      ", 104.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MALO2H      ) = Chemical("MALO2H      ", 147.0000,  0,  5,   0,  0,  0.0000,    0.0 ) 
    species(H2O2        ) = Chemical("H2O2        ",  34.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(CH3COO2H    ) = Chemical("CH3COO2H    ",  76.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(CH3OH       ) = Chemical("CH3OH       ",  32.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(C2H5OH      ) = Chemical("C2H5OH      ",  46.0000,  0,  2,   0,  0,  0.0000,    0.0 ) 
    species(ACETOL      ) = Chemical("ACETOL      ",  74.0000,  0,  3,   0,  0,  0.0000,    0.0 ) 
    species(H2          ) = Chemical("H2          ",   2.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(CO          ) = Chemical("CO          ",  28.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(CH4         ) = Chemical("CH4         ",  16.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(SO2         ) = Chemical("SO2         ",  64.0000,  0,  0,   0,  1,  0.0000,    0.0 ) 
    species(ISO2        ) = Chemical("ISO2        ", 101.0000,  0,  5,   0,  0,  0.0000,    0.0 ) 
    species(MACRO2      ) = Chemical("MACRO2      ", 119.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MACR        ) = Chemical("MACR        ",  70.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(MACROOH     ) = Chemical("MACROOH     ", 120.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(IEPOX       ) = Chemical("IEPOX       ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
    species(HACET       ) = Chemical("HACET       ",  29.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(ISOOH       ) = Chemical("ISOOH       ",  60.0000,  1,  5,   0,  0,  0.0000,    0.0 ) 
    species(ISON        ) = Chemical("ISON        ",  60.0000,  1,  5,   0,  0,  0.0000,    0.0 ) 
    species(HCOOH       ) = Chemical("HCOOH       ",  46.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
    species(MPAN        ) = Chemical("MPAN        ", 132.0000,  0,  4,   1,  0,  0.0000,    0.0 ) 
    species(NALD        ) = Chemical("NALD        ",  60.0000,  1,  5,   0,  0,  0.0000,    0.0 ) 
    species(HPALD       ) = Chemical("HPALD       ",  60.0000,  1,  5,   0,  0,  0.0000,    0.0 ) 
    species(PACALD      ) = Chemical("PACALD      ",  60.0000,  1,  5,   0,  0,  0.0000,    0.0 ) 
    species(MVK         ) = Chemical("MVK         ",  70.0000,  0,  4,   0,  0,  0.0000,    0.0 ) 
    species(TERPOOH     ) = Chemical("TERPOOH     ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 ) 
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
    species(SQT_SOA_NV  ) = Chemical("SQT_SOA_NV  ", 302.0000,  0, 14,   0,  0,  0.0000,    0.0 ) 
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
