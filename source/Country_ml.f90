! <Country_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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
module Country_ml

 ! Sets country index numbers (IC_xx), code, time-zones, and names
 !
 ! Regions external Atlantic (70) and external Russia (71) were outside 
 ! the original EMEP grid (132x111 cells), thus they were defined separately, 
 ! as total emissions for Russia and the Atlantic often are reported for 
 ! the old EMEP domain only and then gridded according to this total.
 !
 ! Area 71 became part of the extended EMEP grid from 2008, thus should 
 ! be included for Russia. Area 70 is not included in the extended grid either.
 !
 ! Special areas have been defined for other projects than EMEP 
 ! outside the EMEP domain (see comments in the module).  
 !
 ! timefac_index under cc as defined below assigns timefactors to 
 ! country/region/emission_type. As an example defining Bavaria as 
 ! a separate region with timefactors as in Germany.

  implicit none

  public :: Country_Init     ! sets country details

  integer, parameter, public  :: NLAND = 350  ! max number of countries 
  logical, parameter, private :: T = .true.   ! shorthand
  logical, parameter, private :: F = .false.  ! shorthand


  ! Some regions
!QUer AL, HR, CS.....

   character(len=10), public, parameter :: &
     EU15(15) = (/ "AT", "BE", "DK", "FI", "FR", "DE", "GR", "IE", "IT", &
               "NL", "PT", "ES", "SE",  "GB", "LU" /),&
     EU27(27) = (/ EU15,"HU", "PL", "CY", "CZ", "EE", "LT", "LV", "MT", &
                        "SK", "SI", "BG", "RO" /),&
     EEA(30)  = (/ EU27, "NO", "IS", "LI" /)  ,&
  ! Countries fully inside MACC2 emission area, excluding EEA.
  !                     1     2     3    4     5     6     7     8     9   10
     XMACC2(10) = (/ "CH", "MC", "TR", "MD", "GE", "AM", "AZ", "BA", "UA", "BY" /),&  ! all?
     EUMACC2(40)  = (/ EEA, XMACC2 /)


  !/ to be set in Country_Init:

  type, public :: cc
     character(len=10)  :: code          ! up to 3 letter land code
     integer           :: index         ! index number (corresponds to 
                                        ! numbering in emission files)
     logical           :: is_sea        ! 1 for sea area, 0 otherwise
     integer           :: timefac_index ! see explanation above
     integer           :: timezone      ! timezone, deviation from UTC time
     character(len=30) :: name          ! name of country/region
  end type cc

  type(cc), public, save, dimension(NLAND) :: Country

  integer, parameter, public ::  IC_AL =   1   ! Albania                       
  integer, parameter, public ::  IC_AT =   2   ! Austria                       
  integer, parameter, public ::  IC_BE =   3   ! Belgium                       
  integer, parameter, public ::  IC_BG =   4   ! Bulgaria                      
  integer, parameter, public :: IC_FCS =   5   ! Former Czechoslovakia     
  integer, parameter, public ::  IC_DK =   6   ! Denmark                       
  integer, parameter, public ::  IC_FI =   7   ! Finland                       
  integer, parameter, public ::  IC_FR =   8   ! France                        
  integer, parameter, public :: IC_GDR =   9   ! Former East Germany           
  integer, parameter, public :: IC_FRG =  10   ! Former West Germany           
  integer, parameter, public ::  IC_GR =  11   ! Greece                        
  integer, parameter, public ::  IC_HU =  12   ! Hungary                       
  integer, parameter, public ::  IC_IS =  13   ! Iceland                       
  integer, parameter, public ::  IC_IE =  14   ! Ireland                       
  integer, parameter, public ::  IC_IT =  15   ! Italy                         
  integer, parameter, public ::  IC_LU =  16   ! Luxembourg                    
  integer, parameter, public ::  IC_NL =  17   ! Netherlands                   
  integer, parameter, public ::  IC_NO =  18   ! Norway                        
  integer, parameter, public ::  IC_PL =  19   ! Poland                        
  integer, parameter, public ::  IC_PT =  20   ! Portugal                      
  integer, parameter, public ::  IC_RO =  21   ! Romania                       
  integer, parameter, public ::  IC_ES =  22   ! Spain                         
  integer, parameter, public ::  IC_SE =  23   ! Sweden                        
  integer, parameter, public ::  IC_CH =  24   ! Switzerland                   
  integer, parameter, public ::  IC_TR =  25   ! Turkey                        
  integer, parameter, public ::  IC_SU =  26   ! Former USSE                   
  integer, parameter, public ::  IC_GB =  27   ! United Kingdom                
  integer, parameter, public :: IC_VUL =  28   ! Vulcanoes                     
  integer, parameter, public :: IC_REM =  29   ! Remaining Areas              
  integer, parameter, public :: IC_BAS =  30   ! The Baltic Sea                
  integer, parameter, public :: IC_NOS =  31   ! The North Sea                 
  integer, parameter, public :: IC_ATL =  32   ! NE Atlantic Ocean (within EMEP domain)  
  integer, parameter, public :: IC_MED =  33   ! The Mediterranean Sea         
  integer, parameter, public :: IC_BLS =  34   ! The Black Sea                 
  integer, parameter, public :: IC_NAT =  35   ! Natural marine sources        
  integer, parameter, public :: IC_RUO =  36   ! Kola/Karelia                  
  integer, parameter, public :: IC_RUP =  37   ! St.Petersburg/Novgorod-Pskov  
  integer, parameter, public :: IC_RUA =  38   ! Kaliningrad                   
  integer, parameter, public ::  IC_BY =  39   ! Belarus                       
  integer, parameter, public ::  IC_UA =  40   ! Ukraine                       
  integer, parameter, public ::  IC_MD =  41   ! Moldova,                      
  integer, parameter, public :: IC_RUR =  42   ! Rest                          
  integer, parameter, public ::  IC_EE =  43   ! Estonia                       
  integer, parameter, public ::  IC_LV =  44   ! Latvia
  integer, parameter, public ::  IC_LT =  45   ! Lithuania                     
  integer, parameter, public ::  IC_CZ =  46   ! Czech                         
  integer, parameter, public ::  IC_SK =  47   ! Slovakia                      
  integer, parameter, public ::  IC_SI =  48   ! Slovenia                      
  integer, parameter, public ::  IC_HR =  49   ! Croatia                       
  integer, parameter, public ::  IC_BA =  50   ! Bosnia                        
  integer, parameter, public ::  IC_CS =  51   ! Serbia and Montenegro    
  integer, parameter, public ::  IC_MK =  52   ! Macedonia                    
  integer, parameter, public ::  IC_KZ =  53   ! Kazakstan                     
  integer, parameter, public ::  IC_GE =  54   ! Georgia                       
  integer, parameter, public ::  IC_CY =  55   ! Cyprus                        
  integer, parameter, public ::  IC_AM =  56   ! Armenia                       
  integer, parameter, public ::  IC_MT =  57   ! Malta                         
  integer, parameter, public :: IC_ASI =  58   ! Other Asian Areas             
  integer, parameter, public ::  IC_LI =  59   ! Lihtenstein                   
  integer, parameter, public ::  IC_DE =  60   ! Germany                       
  integer, parameter, public ::  IC_RU =  61   ! Russian                       
  integer, parameter, public ::  IC_MC =  62   ! Monaco                        
  integer, parameter, public :: IC_NOA =  63   ! North Africa                  
  integer, parameter, public ::  IC_EU =  64   ! European Union
  integer, parameter, public ::  IC_US =  65   ! USA
  integer, parameter, public ::  IC_CA =  66   ! Canada
  integer, parameter, public ::  IC_DUMMY1 =  67 ! Not-defined
  integer, parameter, public ::  IC_KG =  68   ! Kyrgyzstan 
  integer, parameter, public ::  IC_AZ =  69   ! Azerbaijan                 
  integer, parameter, public :: IC_ATX =  70   ! ATL outside EMEP domain
  integer, parameter, public :: IC_RUX =  71   ! RU outside old EMEP domain
  integer, parameter, public :: IC_RS  =  72   ! Serbia
  integer, parameter, public :: IC_ME  =  73   ! Montenegro
  integer, parameter, public :: IC_RUE =  93   ! Russian Federation in the extended EMEP domain (RU+RFE+RUX)

!Extra cc for rest CityZen
  integer, parameter, public :: IC_RAA  =  170   ! Rest of Africa and Asia
  integer, parameter, public :: IC_SEA  =  171   ! Ship

  ! Biomass-burnung (Wild-fires etc.) allocated to a country-number
  ! Allows easy use of emissplits to allocate speciation

  integer, parameter, public :: IC_BB  =  101 

! Extra from IIASA/ECLIPSE/ECLAIRE global
integer, parameter, public :: IC_AFGH = 201 ! Afghanistan
integer, parameter, public :: IC_ARGE = 202 ! Argentina
integer, parameter, public :: IC_AUTR = 203 ! Australia
integer, parameter, public :: IC_BANG = 204 ! Bangladesh
integer, parameter, public :: IC_BHUT = 205 ! Bhutan
integer, parameter, public :: IC_BRAZ = 206 ! Brazil
integer, parameter, public :: IC_BRUN = 207 ! Brunei
integer, parameter, public :: IC_CAMB = 208 ! Cambodia
integer, parameter, public :: IC_CHIL = 209 ! Chile
integer, parameter, public :: IC_CHIN = 210 ! China
integer, parameter, public :: IC_FSUA = 211 ! Former_USSR_(Asia)_Tajikistan_Turkmenistan_Uzbekistan
integer, parameter, public :: IC_INDI = 212 ! India
integer, parameter, public :: IC_INDO = 213 ! Indonesia
integer, parameter, public :: IC_ISRA = 214 ! Israel
integer, parameter, public :: IC_JAPA = 215 ! Japan
integer, parameter, public :: IC_LAOS = 216 ! Laos
integer, parameter, public :: IC_MALA = 217 ! Malaysia
integer, parameter, public :: IC_MEXI = 218 ! Mexico
integer, parameter, public :: IC_MIDE = 219 ! Middle_East
integer, parameter, public :: IC_MONG = 220 ! Mongolia
integer, parameter, public :: IC_MYAN = 221 ! Myanmar
integer, parameter, public :: IC_NEPA = 222 ! Nepal
integer, parameter, public :: IC_NZEL = 223 ! New_Zealand
integer, parameter, public :: IC_NAFR = 224 ! North_Africa_Libya_Tunisia_Algeria_Sudan_Morocco
integer, parameter, public :: IC_KORN = 225 ! North_Korea
integer, parameter, public :: IC_OAFR = 226 ! Other_Africa
integer, parameter, public :: IC_OLAM = 227 ! Other_Latin_America
integer, parameter, public :: IC_PAKI = 228 ! Pakistan
integer, parameter, public :: IC_PHIL = 229 ! Philippines
integer, parameter, public :: IC_SING = 230 ! Singapore
integer, parameter, public :: IC_SAFR = 231 ! South_Africa
integer, parameter, public :: IC_KORS = 232 ! South_Korea
integer, parameter, public :: IC_SRIL = 233 ! Sri_Lanka
integer, parameter, public :: IC_TAIW = 234 ! Taiwan
integer, parameter, public :: IC_THAI = 235 ! Thailand
integer, parameter, public :: IC_VIET = 236 ! Vietnam
integer, parameter, public :: IC_EGYP = 237 ! Egypt

  ! extra subdivisions of ship emissions into shipping categories:
  ! Baltic Sea  (30)
  integer, parameter, public :: IC_BA2 = 302
  integer, parameter, public :: IC_BA3 = 303
  integer, parameter, public :: IC_BA4 = 304
  integer, parameter, public :: IC_BA5 = 305
  integer, parameter, public :: IC_BA6 = 306
  integer, parameter, public :: IC_BA7 = 307
  integer, parameter, public :: IC_BA8 = 308
  integer, parameter, public :: IC_BA9 = 309

  ! North Sea  (31)
  integer, parameter, public :: IC_NS2 = 312
  integer, parameter, public :: IC_NS3 = 313
  integer, parameter, public :: IC_NS4 = 314
  integer, parameter, public :: IC_NS5 = 315
  integer, parameter, public :: IC_NS6 = 316
  integer, parameter, public :: IC_NS7 = 317
  integer, parameter, public :: IC_NS8 = 318
  integer, parameter, public :: IC_NS9 = 319

  ! NE Atlantic  (32)
  integer, parameter, public :: IC_AT2 = 322
  integer, parameter, public :: IC_AT3 = 323
  integer, parameter, public :: IC_AT4 = 324
  integer, parameter, public :: IC_AT5 = 325
  integer, parameter, public :: IC_AT6 = 326
  integer, parameter, public :: IC_AT7 = 327
  integer, parameter, public :: IC_AT8 = 328
  integer, parameter, public :: IC_AT9 = 329

  ! Mediterranean   (33)
  integer, parameter, public :: IC_ME2 = 332
  integer, parameter, public :: IC_ME3 = 333
  integer, parameter, public :: IC_ME4 = 334
  integer, parameter, public :: IC_ME5 = 335
  integer, parameter, public :: IC_ME6 = 336
  integer, parameter, public :: IC_ME7 = 337
  integer, parameter, public :: IC_ME8 = 338
  integer, parameter, public :: IC_ME9 = 339

  ! Black Sea   (34)
  integer, parameter, public :: IC_BL2 = 342
  integer, parameter, public :: IC_BL3 = 343
  integer, parameter, public :: IC_BL4 = 344
  integer, parameter, public :: IC_BL5 = 345
  integer, parameter, public :: IC_BL6 = 346
  integer, parameter, public :: IC_BL7 = 347
  integer, parameter, public :: IC_BL8 = 348
  integer, parameter, public :: IC_BL9 = 349

  ! Ship emissions when sea areas are not divided
  ! Eg. TNO emissions (added on 25th March 2009)
  integer, parameter, public :: IC_INTSHIPS = 350  

  ! New codes defined for the extended EMEP area in 2008
  integer, parameter, public :: IC_RFE = 74 ! Rest of extended Russian Federation 
                                            ! (in the extended EMEP domain)
  integer, parameter, public :: IC_KZE = 75 ! Rest of Kazakhstan 
                                            ! (in the extended EMEP domain)
  integer, parameter, public :: IC_UZ  = 76 ! Uzbekistan (in the orig. EMEP domain)
  integer, parameter, public :: IC_TM  = 77 ! Turkmenistan (in the orig. EMEP domain)
  integer, parameter, public :: IC_UZE = 78 ! Rest of  Uzbekistan 
                                            ! (in the extended EMEP domain)
  integer, parameter, public :: IC_TME = 79 ! Rest of Turkmenistan 
                                            ! (in the extended EMEP domain)
  integer, parameter, public :: IC_CAS = 80 ! Caspian Sea (in the orig. EMEP domain)
  integer, parameter, public :: IC_TJ  = 81 ! Tajikistan (in the extended EMEP domain)
  integer, parameter, public :: IC_ARL = 82 ! Aral Lake (in the orig. EMEP domain)
  integer, parameter, public :: IC_ARE = 83 ! Rest of Aral Lake 
                                            ! (in the extended EMEP domain)
  integer, parameter, public :: IC_ASM = 84 ! Modified remaining Asian areas 
                                            ! (in the original EMEP domain)
  integer, parameter, public :: IC_ASE = 85 ! Remaining extended Asian areas 
                                            ! (in the extended EMEP domain)
  integer, parameter, public :: IC_AOE = 86 ! Arctic Ocean (in the extended EMEP domain)

 ! New external areas (outside the 132x159 grid), these are normally not used
 ! a) Domains: x = 160-170 y = 1-132 and x = -16-0  y = 123-170
  integer, parameter, public :: IC_RFX = 87  ! Extended EMEP-external part of 
                                             ! Russian Federation
  integer, parameter, public :: IC_ASX = 88  ! Extended EMEP-ext. part of Asia
  integer, parameter, public :: IC_PAX = 89  ! Extended EMEP-ext. part of Pacific Ocean
  integer, parameter, public :: IC_AOX = 90  ! Extended EMEP-ext. part of Arctic Ocean

 !b) Domain x = -16-132 y = -11-0
  integer, parameter, public :: IC_NAX = 91  ! EMEP-external part of North Africa

 ! New code introduced for the NMR-NH3 project, not used in other projects
 ! NH3Emis x=44-75, y=35-66
  integer, parameter, public :: IC_NMR = 98  ! EMEP NMR-NH3 temporal emissions

  contains
  
  subroutine Country_Init()

  ! Set the country details. Note that time-zones for some areas are either
  ! difficult (e.g. Russia should be 3 to 12) or not relevant (e.g. sea areas,
  ! volcanoes). 

  integer :: iland

  ! First define all countries as undefined
  do iland=1,NLAND
    Country(iland) = cc(  "N/A" , iland ,F,  17 , 0  , "Not_defined                   " )
  enddo

  !--------------  code index sea timefac_index timezone  Name  ------------!
  
Country( IC_AL ) = cc(  "AL " ,  1 ,F,  1,  1  , "Albania                       " )
Country( IC_AT ) = cc(  "AT " ,  2 ,F,  2,  1  , "Austria                       " )
Country( IC_BE ) = cc(  "BE " ,  3 ,F,  3,  1  , "Belgium                       " )
Country( IC_BG ) = cc(  "BG " ,  4 ,F,  4,  2  , "Bulgaria                      " )
Country( IC_FCS ) = cc(  "FCS " ,  5 ,F,  5,  1 , "Former Czechoslovakia         " )
Country( IC_DK ) = cc(  "DK " ,  6 ,F,  6,  1  , "Denmark                       " )
Country( IC_FI ) = cc(  "FI " ,  7 ,F,  7,  2  , "Finland                       " )
Country( IC_FR ) = cc(  "FR " ,  8 ,F,  8,  1  , "France                        " )
Country( IC_GDR) = cc(  "GDR" ,  9 ,F,  9,  1  , "Former East Germany           " )
Country( IC_FRG) = cc(  "FRG" , 10 ,F, 10,  1  , "Former Fed. Rep. of Germany   " )
Country( IC_GR ) = cc(  "GR " , 11 ,F, 11,  2  , "Greece                        " )
Country( IC_HU ) = cc(  "HU " , 12 ,F, 12,  1  , "Hungary                       " )
Country( IC_IS ) = cc(  "IS " , 13 ,F, 13,  0  , "Iceland                       " )
Country( IC_IE ) = cc(  "IE " , 14 ,F, 14,  0  , "Ireland                       " )
Country( IC_IT ) = cc(  "IT " , 15 ,F, 15,  1  , "Italy                         " )
Country( IC_LU ) = cc(  "LU " , 16 ,F, 16,  1  , "Luxembourg                    " )
Country( IC_NL ) = cc(  "NL " , 17 ,F, 17,  1  , "Netherlands                   " )
Country( IC_NO ) = cc(  "NO " , 18 ,F, 18,  1  , "Norway                        " )
Country( IC_PL ) = cc(  "PL " , 19 ,F, 19,  1  , "Poland                        " )
Country( IC_PT ) = cc(  "PT " , 20 ,F, 20,  0  , "Portugal                      " )
Country( IC_RO ) = cc(  "RO " , 21 ,F, 21,  2  , "Romania                       " )
Country( IC_ES ) = cc(  "ES " , 22 ,F, 22,  1  , "Spain                         " )
Country( IC_SE ) = cc(  "SE " , 23 ,F, 23,  1  , "Sweden                        " )
Country( IC_CH ) = cc(  "CH " , 24 ,F, 24,  1  , "Switzerland                   " )
Country( IC_TR ) = cc(  "TR " , 25 ,F, 25,  2  , "Turkey                        " )
Country( IC_SU ) = cc(  "SU " , 26 ,F, 26,  -100  , "Former USSR                   " )
Country( IC_GB ) = cc(  "GB " , 27 ,F, 27,  0  , "United Kingdom                " )
Country( IC_VUL) = cc(  "VUL" , 28 ,F, 28,  1  , "Volcanoes                     " )
Country( IC_REM) = cc(  "REM" , 29 ,F, 29,  1  , "Remaining land areas          " )
!NB:
!Fix needed for following sea-areas (BAS,NOS,ATL,MED,BLS)in GEA runs done in Emissions_ml
!if ( DomainName == "HIRHAM" .and. IIFULLDOM == 182 ) then ! Special fix for HIRHAM/GEA
!if ( SEAFIX_GEA_NEEDED ) then ! Special fix for HIRHAM/GEA
Country( IC_BAS) = cc(  "BAS" , 30 ,T, 30,  1  , "The Baltic Sea                " )
Country( IC_NOS) = cc(  "NOS" , 31 ,T, 31,  1  , "The North Sea                 " )
Country( IC_ATL) = cc(  "ATL" , 32 ,T, 32,  -100  , "Remaining NE Atlantic Ocean   " )
Country( IC_MED) = cc(  "MED" , 33 ,T, 33,  1  , "The Mediterranean Sea         " )
Country( IC_BLS) = cc(  "BLS" , 34 ,T, 34,  2  , "The Black Sea                 " )

!end if ! HIRHAM/GEA fix

Country( IC_NAT) = cc(  "NAT" , 35 ,F, 35,  -100  , "Natural marine sources        " )
Country( IC_RUO) = cc(  "RUO" , 36 ,F, 36,  4  , "Kola/Karelia                  " )
Country( IC_RUP) = cc(  "RUP" , 37 ,F, 37,  4  , "St.Petersburg/Novgorod-Pskov  " )
Country( IC_RUA) = cc(  "RUA" , 38 ,F, 38,  3  , "Kaliningrad                   " )
Country( IC_BY ) = cc(  "BY " , 39 ,F, 39,  3  , "Belarus                       " )
Country( IC_UA ) = cc(  "UA " , 40 ,F, 40,  2  , "Ukraine                       " )
Country( IC_MD ) = cc(  "MD " , 41 ,F, 41,  2  , "Moldova, Republic of          " )
Country( IC_RUR) = cc(  "RUR" , 42 ,F, 42,  4  , "Rest of Russia                " )
Country( IC_EE ) = cc(  "EE " , 43 ,F, 43,  2  , "Estonia                       " )
Country( IC_LV ) = cc(  "LV " , 44 ,F, 44,  2  , "Latvia                        " )
Country( IC_LT ) = cc(  "LT " , 45 ,F, 45,  2  , "Lithuania                     " )
Country( IC_CZ ) = cc(  "CZ " , 46 ,F, 46,  1  , "Czech                         " )
Country( IC_SK ) = cc(  "SK " , 47 ,F, 47,  1  , "Slovakia                      " )
Country( IC_SI ) = cc(  "SI " , 48 ,F, 48,  1  , "Slovenia                      " )
Country( IC_HR ) = cc(  "HR " , 49 ,F, 49,  1  , "Croatia                       " )
Country( IC_BA ) = cc(  "BA " , 50 ,F, 50,  1  , "Bosnia and Herzegovina        " )
Country( IC_CS ) = cc(  "CS " , 51 ,F, 51,  1  , "Serbia and Montenegro    " )
Country( IC_MK ) = cc(  "MK " , 52 ,F, 52,  1  , "Macedonia, The F.Yugo.Rep. of " )
Country( IC_KZ ) = cc(  "KZ " , 53 ,F, 53,  -100  , "Kazakstan                     " )
Country( IC_GE ) = cc(  "GE " , 54 ,F, 54,  4  , "Georgia                       " )
Country( IC_CY ) = cc(  "CY " , 55 ,F, 55,  2  , "Cyprus                        " )
Country( IC_AM ) = cc(  "AM " , 56 ,F, 56,  4  , "Armenia                       " )
Country( IC_MT ) = cc(  "MT " , 57 ,F, 57,  1  , "Malta                         " )
Country( IC_ASI) = cc(  "ASI" , 58 ,F, 58,  -100  , "Other Asian areas             " )
Country( IC_LI ) = cc(  "LI " , 59 ,F, 59,  1  , "Lichtenstein                  " )
Country( IC_DE ) = cc(  "DE " , 60 ,F, 60,  1  , "Germany                       " )
Country( IC_RU ) = cc(  "RU " , 61 ,F, 61,  -100  , "Russian Federation            " )
Country( IC_MC ) = cc(  "MC " , 62 ,F, 62,  1  , "Monaco                        " )
Country( IC_NOA) = cc(  "NOA" , 63 ,F, 63,  1  , "North Africa                  " )
Country( IC_EU ) = cc(  "EU " , 64 ,F, 64,  1  , "European Community            " )
Country( IC_US ) = cc(  "US " , 65 ,F, 65,  -100  , "USA                           " )
Country( IC_CA ) = cc(  "CA " , 66 ,F, 66,  -100  , "Canada                        " )
Country( IC_DUMMY1 ) &
                 = cc(  "N/A" , 67 ,F,  67, -100  , "Not_defined                   " )
Country( IC_KG ) = cc(  "KG " , 68 ,F,  68, 6  , "Kyrgyzstan                    " )
Country( IC_AZ ) = cc(  "AZ " , 69 ,F,  69, 4  , "Azerbaijan                    " )
Country( IC_ATX) = cc(  "ATX" , 70 ,T,  32, -100  , "Atlantic outside. EMEP        " )
Country( IC_RUX) = cc(  "RUX" , 71 ,F,  42, -100  , "Russian Fed. outside emep     " )
Country( IC_RS)  = cc(  "RS " , 72 ,F,  72, 1  , "Serbia                        " )
Country( IC_ME)  = cc(  "ME " , 73 ,F,  73, 1  , "Montenegro                    " )
Country( IC_RUE) = cc(  "RUE" , 93 ,F,  93, -100 , "Russian Federeation (all)   " )   
! Biomass burning
Country( IC_BB)  = cc(  "BB ", 101,F,  101, -100  , "Biomass burning (wild)        " )

! Sea areas split according to innside/outside 12 nautical mile zone, 
! ferries/cargo ships, registred inside/outside EU
Country( IC_BA2 ) = cc(  "BA2" ,302 ,T,  30, 1  , "Baltic EU cargo outs.12      " )
Country( IC_BA3 ) = cc(  "BA3" ,303 ,T,  30, 1  , "Baltic ROW cargo outs. 12    " )
Country( IC_BA4 ) = cc(  "BA4" ,304 ,T,  30, 1  , "Baltic EU cargo ins. 12      " )
Country( IC_BA5 ) = cc(  "BA5" ,305 ,T,  30, 1  , "Baltic ROW cargo ins. 12     " )
Country( IC_BA6 ) = cc(  "BA6" ,306 ,T,  30, 1  , "Baltic EU ferries outs.12    " )
Country( IC_BA7 ) = cc(  "BA7" ,307 ,T,  30, 1  , "Baltic ROW ferries outs. 12  " )
Country( IC_BA8 ) = cc(  "BA8" ,308 ,T,  30, 1  , "Baltic EU ferries ins. 12    " )
Country( IC_BA9 ) = cc(  "BA9" ,309 ,T,  30, 1  , "Baltic ROW ferries ins. 12   " )

Country( IC_NS2 ) = cc(  "NS2" ,312 ,T,  31, 1  , "N. Sea EU cargo outs.12      " )
Country( IC_NS3 ) = cc(  "NS3" ,313 ,T,  31, 1  , "N. Sea ROW cargo outs. 12    " )
Country( IC_NS4 ) = cc(  "NS4" ,314 ,T,  31, 1  , "N. Sea EU cargo ins. 12      " )
Country( IC_NS5 ) = cc(  "NS5" ,315 ,T,  31, 1  , "N. Sea ROW cargo ins. 12     " )
Country( IC_NS6 ) = cc(  "NS6" ,316 ,T,  31, 1  , "N. Sea EU ferries outs.12    " )
Country( IC_NS7 ) = cc(  "NS7" ,317 ,T,  31, 1  , "N. Sea ROW ferries outs. 12  " )
Country( IC_NS8 ) = cc(  "NS8" ,318 ,T,  31, 1  , "N. Sea EU ferries ins. 12    " )
Country( IC_NS9 ) = cc(  "NS9" ,319 ,T,  31, 1  , "N. Sea ROW ferries ins. 12   " )

Country( IC_AT2 ) = cc(  "AT2" ,322 ,T,  32, 1  , "Atlant EU cargo outs.12      " )
Country( IC_AT3 ) = cc(  "AT3" ,323 ,T,  32, 1  , "Atlant ROW cargo outs. 12    " )
Country( IC_AT4 ) = cc(  "AT4" ,324 ,T,  32, 1  , "Atlant EU cargo ins. 12      " )
Country( IC_AT5 ) = cc(  "AT5" ,325 ,T,  32, 1  , "Atlant ROW cargo ins. 12     " )
Country( IC_AT6 ) = cc(  "AT6" ,326 ,T,  32, 1  , "Atlant EU ferries outs.12    " )
Country( IC_AT7 ) = cc(  "AT7" ,327 ,T,  32, 1  , "Atlant ROW ferries outs. 12  " )
Country( IC_AT8 ) = cc(  "AT8" ,328 ,T,  32, 1  , "Atlant EU ferries ins. 12    " )
Country( IC_AT9 ) = cc(  "AT9" ,329 ,T,  32, 1  , "Atlant ROW ferries ins. 12   " )

Country( IC_ME2 ) = cc(  "ME2" ,332 ,T,  33, 1  , "Medite EU cargo outs.12      " )
Country( IC_ME3 ) = cc(  "ME3" ,333 ,T,  33, 1  , "Medite ROW cargo outs. 12    " )
Country( IC_ME4 ) = cc(  "ME4" ,334 ,T,  33, 1  , "Medite EU cargo ins. 12      " )
Country( IC_ME5 ) = cc(  "ME5" ,335 ,T,  33, 1  , "Medite ROW cargo ins. 12     " )
Country( IC_ME6 ) = cc(  "ME6" ,336 ,T,  33, 1  , "Medite EU ferries outs.12    " )
Country( IC_ME7 ) = cc(  "ME7" ,337 ,T,  33, 1  , "Medite ROW ferries outs. 12  " )
Country( IC_ME8 ) = cc(  "ME8" ,338 ,T,  33, 1  , "Medite EU ferries ins. 12    " )
Country( IC_ME9 ) = cc(  "ME9" ,339 ,T,  33, 1  , "Medite ROW ferries ins. 12   " )

Country( IC_BL2 ) = cc(  "BL2" ,342 ,T,  34, 2  , "B. Sea EU cargo outs.12      " )
Country( IC_BL3 ) = cc(  "BL3" ,343 ,T,  34, 2  , "B. Sea ROW cargo outs. 12    " )
Country( IC_BL4 ) = cc(  "BL4" ,344 ,T,  34, 2  , "B. Sea EU cargo ins. 12      " )
Country( IC_BL5 ) = cc(  "BL5" ,345 ,T,  34, 2  , "B. Sea ROW cargo ins. 12     " )
Country( IC_BL6 ) = cc(  "BL6" ,346 ,T,  34, 2  , "B. Sea EU ferries outs.12    " )
Country( IC_BL7 ) = cc(  "BL7" ,347 ,T,  34, 2  , "B. Sea ROW ferries outs. 12  " )
Country( IC_BL8 ) = cc(  "BL8" ,348 ,T,  34, 2  , "B. Sea EU ferries ins. 12    " )
Country( IC_BL9 ) = cc(  "BL9" ,349 ,T,  34, 2  , "B. Sea ROW ferries ins. 12   " )

! Extended EMEP domain
Country( IC_RFE ) = cc(  "RFE" , 74 ,F, 74, -100  , "Rest of extended Russian Federation (in the extended EMEP domain)" )
Country( IC_KZE ) = cc(  "KZE" , 75 ,F, 75, -100  , "Rest of Kazakhstan (in the extended EMEP domain)                 " )
Country( IC_UZ  ) = cc(  "UZ"  , 76 ,F, 76, 5  , "Uzbekistan (in the original EMEP domain)                         " )
Country( IC_TM  ) = cc(  "TM"  , 77 ,F, 77, 5  , "Turkmenistan (in the original EMEP domain)                       " )
Country( IC_UZE ) = cc(  "UZE" , 78 ,F, 78, 5  , "Rest of  Uzbekistan (in the extended EMEP domain)                " )
Country( IC_TME ) = cc(  "TME" , 79 ,F, 79, 5  , "Rest of Turkmenistan (in the extended EMEP domain)               " )
Country( IC_CAS ) = cc(  "CAS" , 80 ,F, 80, 4  , "Caspian Sea (in the original EMEP domain)                        " )
Country( IC_TJ  ) = cc(  "TJ"  , 81 ,F, 81, 5   ,"Tajikistan (in the extended EMEP domain)                         " )
Country( IC_ARL ) = cc(  "ARL" , 82 ,F, 82, 5  , "Aral Lake (in the original EMEP domain)                          " )
Country( IC_ARE ) = cc(  "ARE" , 83 ,F, 83, 5  , "Rest of Aral Lake (in the extended EMEP domain)                  " )
Country( IC_ASM ) = cc(  "ASM" , 84 ,F, 84, -100  , "Modified remaining Asian areas (in the original EMEP domain)     " )
Country( IC_ASE ) = cc(  "ASE" , 85 ,F, 85, -100  , "Remaining extended Asian areas (in the extended EMEP domain)     " )
Country( IC_AOE ) = cc(  "AOE" , 86 ,F, 86, 9  , "Arctic Ocean (in the extended EMEP domain)                       " )

! NH3Emis new land code for NMR-NH3 project
Country( IC_NMR ) = cc(  "NMR" , 98 ,F, 98, 1  , "Area with temporal NMR-NH3 emissions              " )

! New external areas (outside the 132x159 grid), these are normally not used
! a) Domains: x = 160-170 y = 1-132 and x =  -16-0  y = 123-170
Country( IC_RFX ) = cc(  "RFX" , 87 ,F,  87, -100  ,"Extended EMEP-external part of Russian Federation" )
Country( IC_ASX ) = cc(  "ASX" , 88 ,F,  88, -100   ,"Extended EMEP-external part of Asia              " )
Country( IC_PAX ) = cc(  "PAX" , 89 ,F,  89, -100  ,"Extended EMEP-external part of Pacific Ocean     " )
Country( IC_AOX ) = cc(  "AOX" , 90 ,F,  90, 9  ,"Extended EMEP-external part of Arctic Ocean      " )
! b) Domain x = -16-132 y = -11-0
Country( IC_NAX ) = cc(  "NAX" , 91 ,F,  91, -100   ,"EMEP-external part of North Africa               " )

!Extra cc for rest CityZen
Country( IC_RAA ) = cc(  "RAA" , 170 ,F,  170, -100, "Rest of Africa and Asia" )
Country( IC_SEA ) = cc(  "SEA" , 171 ,F,  171, -100, "Ships" )

! Extra from IIASA/ECLIPSE/ECLAIRE global

Country(IC_AFGH) = cc( "AFGH", 201, F,201, -100, "Afghanistan") 
Country(IC_ARGE) = cc( "ARGE", 202, F,202, -100, "Argentina") 
Country(IC_AUTR) = cc( "AUTR", 203, F,203, -100, "Australia") 
Country(IC_BANG) = cc( "BANG", 204, F,204, -100, "Bangladesh") 
Country(IC_BHUT) = cc( "BHUT", 205, F,205, -100, "Bhutan") 
Country(IC_BRAZ) = cc( "BRAZ", 206, F,206, -100, "Brazil") 
Country(IC_BRUN) = cc( "BRUN", 207, F,207, -100, "Brunei") 
Country(IC_CAMB) = cc( "CAMB", 208, F,208, -100, "Cambodia") 
Country(IC_CHIL) = cc( "CHIL", 209, F,209, -100, "Chile") 
Country(IC_CHIN) = cc( "CHIN", 210, F,210, -100, "China") 
Country(IC_FSUA) = cc( "FSUA", 211, F,211, -100, "Former_USSR_(Asia)_Tajikistan_Turkmenistan_Uzbekistan") 
Country(IC_INDI) = cc( "INDI", 212, F,212, -100, "India") 
Country(IC_INDO) = cc( "INDO", 213, F,213, -100, "Indonesia") 
Country(IC_ISRA) = cc( "ISRA", 214, F,214, -100, "Israel") 
Country(IC_JAPA) = cc( "JAPA", 215, F,215, -100, "Japan") 
Country(IC_LAOS) = cc( "LAOS", 216, F,216, -100, "Laos") 
Country(IC_MALA) = cc( "MALA", 217, F,217, -100, "Malaysia") 
Country(IC_MEXI) = cc( "MEXI", 218, F,218, -100, "Mexico") 
Country(IC_MIDE) = cc( "MIDE", 219, F,219, -100, "Middle_East") 
Country(IC_MONG) = cc( "MONG", 220, F,220, -100, "Mongolia") 
Country(IC_MYAN) = cc( "MYAN", 221, F,221, -100, "Myanmar") 
Country(IC_NEPA) = cc( "NEPA", 222, F,222, -100, "Nepal") 
Country(IC_NZEL) = cc( "NZEL", 223, F,223, -100, "New_Zealand") 
Country(IC_NAFR) = cc( "NAFR", 224, F,224, -100, "North_Africa_Libya_Tunisia_Algeria_Sudan_Morocco") 
Country(IC_KORN) = cc( "KORN", 225, F,225, -100, "North_Korea") 
Country(IC_OAFR) = cc( "OAFR", 226, F,226, -100, "Other_Africa") 
Country(IC_OLAM) = cc( "OLAM", 227, F,227, -100, "Other_Latin_America") 
Country(IC_PAKI) = cc( "PAKI", 228, F,228, -100, "Pakistan") 
Country(IC_PHIL) = cc( "PHIL", 229, F,229, -100, "Philippines") 
Country(IC_SING) = cc( "SING", 230, F,230, -100, "Singapore") 
Country(IC_SAFR) = cc( "SAFR", 231, F,231, -100, "South_Africa") 
Country(IC_KORS) = cc( "KORS", 232, F,232, -100, "South_Korea") 
Country(IC_SRIL) = cc( "SRIL", 233, F,233, -100, "Sri_Lanka") 
Country(IC_TAIW) = cc( "TAIW", 234, F,234, -100, "Taiwan") 
Country(IC_THAI) = cc( "THAI", 235, F,235, -100, "Thailand") 
Country(IC_VIET) = cc( "VIET", 236, F,236, -100, "Vietnam") 
Country(IC_EGYP) = cc( "EGYP", 237, F,237, -100, "Egypt")
Country(IC_INTSHIPS ) = cc(  "INTSHIPS" ,350 ,T, 350, -100  , "International ships, RCP6" )
  end subroutine Country_Init

end module Country_ml
