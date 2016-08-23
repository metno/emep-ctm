! <Country_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2015 met.no
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

 ! Sets country index numbers (IC_xx), code, icode, time-zones, and names
 !
 ! "icode" should match the number used in the emission files
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

  integer, parameter, public  :: MAXNLAND = 350  ! max number of countries 
  integer,  public            :: NLAND  ! actaua number of countries defined
  logical, parameter, private :: T = .true.   ! shorthand
  logical, parameter, private :: F = .false.  ! shorthand


  ! Some regions
!QUer AL, HR, CS.....

   character(len=10), public, parameter :: &
     EU15(15) = (/ "AT", "BE", "DK", "FI", "FR", "DE", "GR", "IE", "IT", &
               "NL", "PT", "ES", "SE",  "GB", "LU" /),&
     EU27(27) = (/ EU15,"HU", "PL", "CY", "CZ", "EE", "LT", "LV", "MT", &
                        "SK", "SI", "BG", "RO" /),&
     EU28(28) = (/ EU27, "HR" /),&
     EEA(31)  = (/ EU28, "NO", "IS", "LI" /)  ,&
  ! Countries fully inside MACC2 emission area, excluding EEA.
  !                     1     2     3    4     5     6     7     8     9   10
     XMACC2(10) = (/ "CH", "MC", "TR", "MD", "GE", "AM", "AZ", "BA", "UA", "BY" /),&  ! all?
     EUMACC2(41)  = (/ EEA, XMACC2 /)


  !/ to be set in Country_Init:

  type, public :: cc
     character(len=10) :: code          ! up to 3 letter land code
     integer           :: icode         ! integer number for land code (corresponds to 
                                        ! country code number in emission files)
     logical           :: is_sea        ! 1 for sea area, 0 otherwise
     integer           :: timefac_index ! see explanation above
     integer           :: timezone      ! timezone, deviation from UTC time
     character(len=30) :: name          ! name of country/region
  end type cc

  type(cc), public, save, dimension(MAXNLAND) :: Country

  integer, public ::  IC_AL   ! Albania                       
  integer, public ::  IC_AT   ! Austria                       
  integer, public ::  IC_BE   ! Belgium                       
  integer, public ::  IC_BG   ! Bulgaria                      
  integer, public :: IC_FCS   ! Former Czechoslovakia     
  integer, public ::  IC_DK   ! Denmark                       
  integer, public ::  IC_FI   ! Finland                       
  integer, public ::  IC_FR   ! France                        
  integer, public :: IC_GDR   ! Former East Germany           
  integer, public :: IC_FRG   ! Former West Germany           
  integer, public ::  IC_GR   ! Greece                        
  integer, public ::  IC_HU   ! Hungary                       
  integer, public ::  IC_IS   ! Iceland                       
  integer, public ::  IC_IE   ! Ireland                       
  integer, public ::  IC_IT   ! Italy                         
  integer, public ::  IC_LU   ! Luxembourg                    
  integer, public ::  IC_NL   ! Netherlands                   
  integer, public ::  IC_NO   ! Norway                        
  integer, public ::  IC_PL   ! Poland                        
  integer, public ::  IC_PT   ! Portugal                      
  integer, public ::  IC_RO   ! Romania                       
  integer, public ::  IC_ES   ! Spain                         
  integer, public ::  IC_SE   ! Sweden                        
  integer, public ::  IC_CH   ! Switzerland                   
  integer, public ::  IC_TR   ! Turkey                        
  integer, public ::  IC_SU   ! Former USSE                   
  integer, public ::  IC_GB   ! United Kingdom                
  integer, public :: IC_VUL   ! Vulcanoes                     
  integer, public :: IC_REM   ! Remaining Areas              
  integer, public :: IC_BAS   ! The Baltic Sea                
  integer, public :: IC_NOS   ! The North Sea                 
  integer, public :: IC_ATL   ! NE Atlantic Ocean (within EMEP domain)  
  integer, public :: IC_MED   ! The Mediterranean Sea         
  integer, public :: IC_BLS   ! The Black Sea                 
  integer, public :: IC_NAT   ! Natural marine sources        
  integer, public :: IC_RUO   ! Kola/Karelia                  
  integer, public :: IC_RUP   ! St.Petersburg/Novgorod-Pskov  
  integer, public :: IC_RUA   ! Kaliningrad                   
  integer, public ::  IC_BY   ! Belarus                       
  integer, public ::  IC_UA   ! Ukraine                       
  integer, public ::  IC_MD   ! Moldova,                      
  integer, public :: IC_RUR   ! Rest                          
  integer, public ::  IC_EE   ! Estonia                       
  integer, public ::  IC_LV   ! Latvia
  integer, public ::  IC_LT   ! Lithuania                     
  integer, public ::  IC_CZ   ! Czech                         
  integer, public ::  IC_SK   ! Slovakia                      
  integer, public ::  IC_SI   ! Slovenia                      
  integer, public ::  IC_HR   ! Croatia                       
  integer, public ::  IC_BA   ! Bosnia                        
  integer, public ::  IC_CS   ! Serbia and Montenegro    
  integer, public ::  IC_MK   ! Macedonia                    
  integer, public ::  IC_KZ   ! Kazakstan                     
  integer, public ::  IC_GE   ! Georgia                       
  integer, public ::  IC_CY   ! Cyprus                        
  integer, public ::  IC_AM   ! Armenia                       
  integer, public ::  IC_MT   ! Malta                         
  integer, public :: IC_ASI   ! Other Asian Areas             
  integer, public ::  IC_LI   ! Lihtenstein                   
  integer, public ::  IC_DE   ! Germany                       
  integer, public ::  IC_RU   ! Russian                       
  integer, public ::  IC_MC   ! Monaco                        
  integer, public :: IC_NOA   ! North Africa                  
  integer, public ::  IC_EU   ! European Union
  integer, public ::  IC_US   ! USA
  integer, public ::  IC_CA   ! Canada
  integer, public ::  IC_DUMMY1  ! Not-defined
  integer, public ::  IC_KG   ! Kyrgyzstan 
  integer, public ::  IC_AZ   ! Azerbaijan                 
  integer, public :: IC_ATX   ! ATL outside EMEP domain
  integer, public :: IC_RUX   ! RU outside old EMEP domain
  integer, public :: IC_RS    ! Serbia
  integer, public :: IC_ME    ! Montenegro
  integer, public :: IC_RUE   ! Russian Federation in the extended EMEP domain (RU+RFE+RUX)

!Extra cc for rest CityZen
  integer,  public :: IC_RAA    ! Rest of Africa and Asia
  integer,  public :: IC_SEA    ! Ship

  ! Biomass-burnung (Wild-fires etc.) allocated to a country-number
  ! Allows easy use of emissplits to allocate speciation

  integer,  public :: IC_BB  

! Extra from IIASA/ECLIPSE/ECLAIRE global
integer, public :: IC_AFGH  ! Afghanistan
integer, public :: IC_ARGE  ! Argentina
integer, public :: IC_AUTR  ! Australia
integer, public :: IC_BANG  ! Bangladesh
integer, public :: IC_BHUT  ! Bhutan
integer, public :: IC_BRAZ  ! Brazil
integer, public :: IC_BRUN  ! Brunei
integer, public :: IC_CAMB  ! Cambodia
integer, public :: IC_CHIL  ! Chile
integer, public :: IC_CHIN  ! China
integer, public :: IC_FSUA  ! Former_USSR_(Asia)_Tajikistan_Turkmenistan_Uzbekistan
integer, public :: IC_INDI  ! India
integer, public :: IC_INDO  ! Indonesia
integer, public :: IC_ISRA  ! Israel
integer, public :: IC_JAPA  ! Japan
integer, public :: IC_LAOS  ! Laos
integer, public :: IC_MALA  ! Malaysia
integer, public :: IC_MEXI  ! Mexico
integer, public :: IC_MIDE  ! Middle_East
integer, public :: IC_MONG  ! Mongolia
integer, public :: IC_MYAN  ! Myanmar
integer, public :: IC_NEPA  ! Nepal
integer, public :: IC_NZEL  ! New_Zealand
integer, public :: IC_NAFR  ! North_Africa_Libya_Tunisia_Algeria_Sudan_Morocco
integer, public :: IC_KORN  ! North_Korea
integer, public :: IC_OAFR  ! Other_Africa
integer, public :: IC_OLAM  ! Other_Latin_America
integer, public :: IC_PAKI  ! Pakistan
integer, public :: IC_PHIL  ! Philippines
integer, public :: IC_SING  ! Singapore
integer, public :: IC_SAFR  ! South_Africa
integer, public :: IC_KORS  ! South_Korea
integer, public :: IC_SRIL  ! Sri_Lanka
integer, public :: IC_TAIW  ! Taiwan
integer, public :: IC_THAI  ! Thailand
integer, public :: IC_VIET  ! Vietnam
integer, public :: IC_EGYP  ! Egypt

  ! extra subdivisions of ship emissions into shipping categories:
  ! Baltic Sea  (30)
  integer,  public :: IC_BA2 
  integer,  public :: IC_BA3 
  integer,  public :: IC_BA4 
  integer,  public :: IC_BA5 
  integer,  public :: IC_BA6 
  integer,  public :: IC_BA7 
  integer,  public :: IC_BA8 
  integer,  public :: IC_BA9 

  ! North Sea  (31)
  integer,  public :: IC_NS2
  integer,  public :: IC_NS3
  integer,  public :: IC_NS4
  integer,  public :: IC_NS5
  integer,  public :: IC_NS6
  integer,  public :: IC_NS7
  integer,  public :: IC_NS8
  integer,  public :: IC_NS9

  ! NE Atlantic  (32)
  integer, public :: IC_AT2 
  integer, public :: IC_AT3 
  integer, public :: IC_AT4 
  integer, public :: IC_AT5 
  integer, public :: IC_AT6 
  integer, public :: IC_AT7 
  integer, public :: IC_AT8 
  integer, public :: IC_AT9 

  ! Mediterranean   (33)
  integer,  public :: IC_ME2
  integer,  public :: IC_ME3
  integer,  public :: IC_ME4
  integer,  public :: IC_ME5
  integer,  public :: IC_ME6
  integer,  public :: IC_ME7
  integer,  public :: IC_ME8
  integer,  public :: IC_ME9

  ! Black Sea   (34)
  integer,  public :: IC_BL2
  integer,  public :: IC_BL3
  integer,  public :: IC_BL4
  integer,  public :: IC_BL5
  integer,  public :: IC_BL6
  integer,  public :: IC_BL7
  integer,  public :: IC_BL8
  integer,  public :: IC_BL9

  ! Ship emissions when sea areas are not divided
  ! Eg. TNO emissions (added on 25th March 2009)
  integer,  public :: IC_INTSHIPS

  ! New codes defined for the extended EMEP area in 2008
  integer, public :: IC_RFE! Rest of extended Russian Federation 
                           ! (in the extended EMEP domain)
  integer, public :: IC_KZE! Rest of Kazakhstan 
                           ! (in the extended EMEP domain)
  integer, public :: IC_UZ ! Uzbekistan (in the orig. EMEP domain)
  integer, public :: IC_TM ! Turkmenistan (in the orig. EMEP domain)
  integer, public :: IC_UZE! Rest of  Uzbekistan 
                           ! (in the extended EMEP domain)
  integer, public :: IC_TME! Rest of Turkmenistan 
                           ! (in the extended EMEP domain)
  integer, public :: IC_CAS! Caspian Sea (in the orig. EMEP domain)
  integer, public :: IC_TJ ! Tajikistan (in the extended EMEP domain)
  integer, public :: IC_ARL! Aral Lake (in the orig. EMEP domain)
  integer, public :: IC_ARE! Rest of Aral Lake 
                           ! (in the extended EMEP domain)
  integer, public :: IC_ASM! Modified remaining Asian areas 
                           ! (in the original EMEP domain)
  integer, public :: IC_ASE! Remaining extended Asian areas 
                           ! (in the extended EMEP domain)
  integer, public :: IC_AOE! Arctic Ocean (in the extended EMEP domain)

 ! New external areas (outside the 132x159 grid), these are normally not used
 ! a) Domains: x = 160-170 y = 1-132 and x = -16-0  y = 123-170
  integer,  public :: IC_RFX ! Extended EMEP-external part of 
                             ! Russian Federation
  integer,  public :: IC_ASX ! Extended EMEP-ext. part of Asia
  integer,  public :: IC_PAX ! Extended EMEP-ext. part of Pacific Ocean
  integer,  public :: IC_AOX ! Extended EMEP-ext. part of Arctic Ocean

 !b) Domain x = -16-132 y = -11-0
  integer,  public :: IC_NAX  ! EMEP-external part of North Africa

 ! New code introduced for the NMR-NH3 project, not used in other projects
 ! NH3Emis x=44-75, y=35-66
  integer,  public :: IC_NMR  ! EMEP NMR-NH3 temporal emissions

! HTAP2 regions
  integer,  public :: IC_HTNATL
  integer,  public :: IC_HTUSCA
  integer,  public :: IC_HTEUTU
  integer,  public :: IC_HTSASI
  integer,  public :: IC_HTEASI
  integer,  public :: IC_HTSEAS
  integer,  public :: IC_HTAUST
  integer,  public :: IC_HTNAFR
  integer,  public :: IC_HTRAFR
  integer,  public :: IC_HTMIDE
  integer,  public :: IC_HT1012
  integer,  public :: IC_HT1013
  integer,  public :: IC_HT1014
  integer,  public :: IC_HT1015
  integer,  public :: IC_HT1016
  integer,  public :: IC_HT1017
  integer,  public :: IC_HT1018
  integer,  public :: IC_HT1019
  integer,  public :: IC_HT1020
  integer,  public :: IC_HT1000

  contains
  
  subroutine Country_Init()

  ! Set the country details. Note that time-zones for some areas are either
  ! difficult (e.g. Russia should be 3 to 12) or not relevant (e.g. sea areas,
  ! volcanoes). 

  integer :: iland,ix

  ! First define all countries as undefined
  do iland=1,NLAND
    Country(iland) = cc(  "N/A" , iland ,F,  17 , 0  , "Not_defined                   " )
  enddo

!The value of IC_XX is the index in Country array. Can in principle change between two runs or versions.
!The emission_code is the country code used in the emission file.
!The timefac_code is the code/index refering to the timefactor 

  !--------------  code emission_icode sea timefac_code timezone  Name  ------------!
  
ix=0 
ix=ix+1 
IC_AL=ix
Country( IC_AL ) = cc(  "AL " ,  1 ,F,  1,  1  , "Albania                       " )
ix=ix+1 
IC_AT=ix
Country( IC_AT ) = cc(  "AT " ,  2 ,F,  2,  1  , "Austria                       " )
ix=ix+1 
IC_BE=ix
Country( IC_BE ) = cc(  "BE " ,  3 ,F,  3,  1  , "Belgium                       " )
ix=ix+1 
IC_BG=ix
Country( IC_BG ) = cc(  "BG " ,  4 ,F,  4,  2  , "Bulgaria                      " )
ix=ix+1 
IC_FCS=ix
Country( IC_FCS ) = cc(  "FCS " ,  5 ,F,  5,  1 , "Former Czechoslovakia         " )
ix=ix+1 
IC_DK=ix
Country( IC_DK ) = cc(  "DK " ,  6 ,F,  6,  1  , "Denmark                       " )
ix=ix+1 
IC_FI=ix
Country( IC_FI ) = cc(  "FI " ,  7 ,F,  7,  2  , "Finland                       " )
ix=ix+1 
IC_FR=ix
Country( IC_FR ) = cc(  "FR " ,  8 ,F,  8,  1  , "France                        " )
ix=ix+1 
IC_GDR=ix
Country( IC_GDR) = cc(  "GDR" ,  9 ,F,  9,  1  , "Former East Germany           " )
ix=ix+1 
IC_FRG=ix
Country( IC_FRG) = cc(  "FRG" , 10 ,F, 10,  1  , "Former Fed. Rep. of Germany   " )
ix=ix+1 
IC_GR=ix
Country( IC_GR ) = cc(  "GR " , 11 ,F, 11,  2  , "Greece                        " )
ix=ix+1 
IC_HU=ix
Country( IC_HU ) = cc(  "HU " , 12 ,F, 12,  1  , "Hungary                       " )
ix=ix+1 
IC_IS=ix
Country( IC_IS ) = cc(  "IS " , 13 ,F, 13,  0  , "Iceland                       " )
ix=ix+1 
IC_IE=ix
Country( IC_IE ) = cc(  "IE " , 14 ,F, 14,  0  , "Ireland                       " )
ix=ix+1 
IC_IT=ix
Country( IC_IT ) = cc(  "IT " , 15 ,F, 15,  1  , "Italy                         " )
ix=ix+1 
IC_LU=ix
Country( IC_LU ) = cc(  "LU " , 16 ,F, 16,  1  , "Luxembourg                    " )
ix=ix+1 
IC_NL=ix
Country( IC_NL ) = cc(  "NL " , 17 ,F, 17,  1  , "Netherlands                   " )
ix=ix+1 
IC_NO=ix
Country( IC_NO ) = cc(  "NO " , 18 ,F, 18,  1  , "Norway                        " )
ix=ix+1 
IC_PL=ix
Country( IC_PL ) = cc(  "PL " , 19 ,F, 19,  1  , "Poland                        " )
ix=ix+1 
IC_PT=ix
Country( IC_PT ) = cc(  "PT " , 20 ,F, 20,  0  , "Portugal                      " )
ix=ix+1 
IC_RO=ix
Country( IC_RO ) = cc(  "RO " , 21 ,F, 21,  2  , "Romania                       " )
ix=ix+1 
IC_ES =ix
Country( IC_ES ) = cc(  "ES " , 22 ,F, 22,  1  , "Spain                         " )
ix=ix+1 
IC_SE=ix
Country( IC_SE ) = cc(  "SE " , 23 ,F, 23,  1  , "Sweden                        " )
ix=ix+1 
IC_CH=ix
Country( IC_CH ) = cc(  "CH " , 24 ,F, 24,  1  , "Switzerland                   " )
ix=ix+1 
IC_TR=ix
Country( IC_TR ) = cc(  "TR " , 25 ,F, 25,  2  , "Turkey                        " )
ix=ix+1 
IC_SU=ix
Country( IC_SU ) = cc(  "SU " , 26 ,F, 26,  -100  , "Former USSR                   " )
ix=ix+1 
IC_GB=ix
Country( IC_GB ) = cc(  "GB " , 27 ,F, 27,  0  , "United Kingdom                " )
ix=ix+1 
IC_VUL=ix
Country( IC_VUL) = cc(  "VUL" , 28 ,F, 28,  1  , "Volcanoes                     " )
ix=ix+1 
IC_REM=ix
Country( IC_REM) = cc(  "REM" , 29 ,F, 29,  1  , "Remaining land areas          " )
!NB:
!Fix needed for following sea-areas (BAS,NOS,ATL,MED,BLS)in GEA runs done in Emissions_ml
!if ( DomainName == "HIRHAM" .and. IIFULLDOM == 182 ) then ! Special fix for HIRHAM/GEA
!if ( SEAFIX_GEA_NEEDED ) then ! Special fix for HIRHAM/GEA
ix=ix+1 
IC_BAS=ix
Country( IC_BAS) = cc(  "BAS" , 30 ,T, 30,  1  , "The Baltic Sea                " )
ix=ix+1 
IC_NOS=ix
Country( IC_NOS) = cc(  "NOS" , 31 ,T, 31,  1  , "The North Sea                 " )
ix=ix+1 
IC_ATL=ix
Country( IC_ATL) = cc(  "ATL" , 32 ,T, 32,  -100  , "Remaining NE Atlantic Ocean   " )
ix=ix+1 
IC_MED=ix
Country( IC_MED) = cc(  "MED" , 33 ,T, 33,  1  , "The Mediterranean Sea         " )
ix=ix+1 
IC_BLS=ix
Country( IC_BLS) = cc(  "BLS" , 34 ,T, 34,  2  , "The Black Sea                 " )

!end if ! HIRHAM/GEA fix

ix=ix+1 
IC_NAT=ix
Country( IC_NAT) = cc(  "NAT" , 35 ,F, 35,  -100  , "Natural marine sources        " )
ix=ix+1
IC_RUO=ix
Country( IC_RUO) = cc(  "RUO" , 36 ,F, 36,  4  , "Kola/Karelia                  " )
ix=ix+1 
IC_RUP=ix
Country( IC_RUP) = cc(  "RUP" , 37 ,F, 37,  4  , "St.Petersburg/Novgorod-Pskov  " )
ix=ix+1 
IC_RUA=ix
Country( IC_RUA) = cc(  "RUA" , 38 ,F, 38,  3  , "Kaliningrad                   " )
ix=ix+1 
IC_BY=ix
Country( IC_BY ) = cc(  "BY " , 39 ,F, 39,  3  , "Belarus                       " )
ix=ix+1 
IC_UA=ix
Country( IC_UA ) = cc(  "UA " , 40 ,F, 40,  2  , "Ukraine                       " )
ix=ix+1 
IC_MD=ix
Country( IC_MD ) = cc(  "MD " , 41 ,F, 41,  2  , "Moldova, Republic of          " )
ix=ix+1 
IC_RUR=ix
Country( IC_RUR) = cc(  "RUR" , 42 ,F, 42,  4  , "Rest of Russia                " )
ix=ix+1 
IC_EE=ix
Country( IC_EE ) = cc(  "EE " , 43 ,F, 43,  2  , "Estonia                       " )
ix=ix+1 
IC_LV=ix
Country( IC_LV ) = cc(  "LV " , 44 ,F, 44,  2  , "Latvia                        " )
ix=ix+1 
IC_LT=ix
Country( IC_LT ) = cc(  "LT " , 45 ,F, 45,  2  , "Lithuania                     " )
ix=ix+1 
IC_CZ=ix
Country( IC_CZ ) = cc(  "CZ " , 46 ,F, 46,  1  , "Czech                         " )
ix=ix+1 
IC_SK=ix
Country( IC_SK ) = cc(  "SK " , 47 ,F, 47,  1  , "Slovakia                      " )
ix=ix+1 
IC_SI=ix
Country( IC_SI ) = cc(  "SI " , 48 ,F, 48,  1  , "Slovenia                      " )
ix=ix+1 
IC_HR=ix
Country( IC_HR ) = cc(  "HR " , 49 ,F, 49,  1  , "Croatia                       " )
ix=ix+1 
IC_BA=ix
Country( IC_BA ) = cc(  "BA " , 50 ,F, 50,  1  , "Bosnia and Herzegovina        " )
ix=ix+1 
IC_CS=ix
Country( IC_CS ) = cc(  "CS " , 51 ,F, 51,  1  , "Serbia and Montenegro    " )
ix=ix+1 
IC_MK=ix
Country( IC_MK ) = cc(  "MK " , 52 ,F, 52,  1  , "Macedonia, The F.Yugo.Rep. of " )
ix=ix+1 
IC_KZ=ix
Country( IC_KZ ) = cc(  "KZ " , 53 ,F, 53,  -100  , "Kazakstan                     " )
ix=ix+1
IC_GE=ix
Country( IC_GE ) = cc(  "GE " , 54 ,F, 54,  4  , "Georgia                       " )
ix=ix+1 
IC_CY=ix
Country( IC_CY ) = cc(  "CY " , 55 ,F, 55,  2  , "Cyprus                        " )
ix=ix+1 
IC_AM =ix
Country( IC_AM ) = cc(  "AM " , 56 ,F, 56,  4  , "Armenia                       " )
ix=ix+1 
IC_MT=ix
Country( IC_MT ) = cc(  "MT " , 57 ,F, 57,  1  , "Malta                         " )
ix=ix+1 
IC_ASI=ix
Country( IC_ASI) = cc(  "ASI" , 58 ,F, 58,  -100  , "Other Asian areas             " )
ix=ix+1 
IC_LI=ix
Country( IC_LI ) = cc(  "LI " , 59 ,F, 59,  1  , "Lichtenstein                  " )
ix=ix+1 
IC_DE=ix
Country( IC_DE ) = cc(  "DE " , 60 ,F, 60,  1  , "Germany                       " )
ix=ix+1 
IC_RU=ix
Country( IC_RU ) = cc(  "RU " , 61 ,F, 61,  -100  , "Russian Federation            " )
ix=ix+1 
IC_MC=ix
Country( IC_MC ) = cc(  "MC " , 62 ,F, 62,  1  , "Monaco                        " )
ix=ix+1 
IC_NOA=ix
Country( IC_NOA) = cc(  "NOA" , 63 ,F, 63,  1  , "North Africa                  " )
ix=ix+1 
IC_EU=ix
Country( IC_EU ) = cc(  "EU " , 64 ,F, 64,  1  , "European Community            " )
ix=ix+1 
IC_US=ix
Country( IC_US ) = cc(  "US " , 65 ,F, 65,  -100  , "USA                           " )
ix=ix+1 
IC_CA=ix
Country( IC_CA ) = cc(  "CA " , 66 ,F, 66,  -100  , "Canada                        " )
ix=ix+1 
IC_DUMMY1=ix
Country( IC_DUMMY1 ) &
                 = cc(  "N/A" , 67 ,F,  67, -100  , "Not_defined                   " )
ix=ix+1
IC_KG=ix
Country( IC_KG ) = cc(  "KG " , 68 ,F,  68, 6  , "Kyrgyzstan                    " )
ix=ix+1 
IC_AZ=ix
Country( IC_AZ ) = cc(  "AZ " , 69 ,F,  69, 4  , "Azerbaijan                    " )
ix=ix+1 
IC_ATX=ix
Country( IC_ATX) = cc(  "ATX" , 70 ,T,  32, -100  , "Atlantic outside. EMEP        " )
ix=ix+1
IC_RUX=ix
Country( IC_RUX) = cc(  "RUX" , 71 ,F,  42, -100  , "Russian Fed. outside emep     " )
ix=ix+1 
IC_RS=ix
Country( IC_RS)  = cc(  "RS " , 72 ,F,  72, 1  , "Serbia                        " )
ix=ix+1
IC_ME=ix
Country( IC_ME)  = cc(  "ME " , 73 ,F,  73, 1  , "Montenegro                    " )
! Extended EMEP domain
ix=ix+1
IC_RFE=ix
Country( IC_RFE ) = cc(  "RFE" , 74 ,F, 74, -100  , "Rest of extended Russian Federation (in the extended EMEP domain)" )
ix=ix+1 
IC_KZE=ix
Country( IC_KZE ) = cc(  "KZE" , 75 ,F, 75, -100  , "Rest of Kazakhstan (in the extended EMEP domain)                 " )
ix=ix+1 
IC_UZ=ix
Country( IC_UZ  ) = cc(  "UZ"  , 76 ,F, 76, 5  , "Uzbekistan (in the original EMEP domain)                         " )
ix=ix+1 
IC_TM=ix
Country( IC_TM  ) = cc(  "TM"  , 77 ,F, 77, 5  , "Turkmenistan (in the original EMEP domain)                       " )
ix=ix+1 
IC_UZE=ix
Country( IC_UZE ) = cc(  "UZE" , 78 ,F, 78, 5  , "Rest of  Uzbekistan (in the extended EMEP domain)                " )
ix=ix+1 
IC_TME=ix
Country( IC_TME ) = cc(  "TME" , 79 ,F, 79, 5  , "Rest of Turkmenistan (in the extended EMEP domain)               " )
ix=ix+1 
IC_CAS=ix
Country( IC_CAS ) = cc(  "CAS" , 80 ,F, 80, 4  , "Caspian Sea (in the original EMEP domain)                        " )
ix=ix+1 
IC_TJ=ix
Country( IC_TJ  ) = cc(  "TJ"  , 81 ,F, 81, 5   ,"Tajikistan (in the extended EMEP domain)                         " )
ix=ix+1 
IC_ARL=ix
Country( IC_ARL ) = cc(  "ARL" , 82 ,F, 82, 5  , "Aral Lake (in the original EMEP domain)                          " )
ix=ix+1 
IC_ARE=ix
Country( IC_ARE ) = cc(  "ARE" , 83 ,F, 83, 5  , "Rest of Aral Lake (in the extended EMEP domain)                  " )
ix=ix+1 
IC_ASM=ix
Country( IC_ASM ) = cc(  "ASM" , 84 ,F, 84, -100  , "Modified remaining Asian areas (in the original EMEP domain)     " )
ix=ix+1 
IC_ASE=ix
Country( IC_ASE ) = cc(  "ASE" , 85 ,F, 85, -100  , "Remaining extended Asian areas (in the extended EMEP domain)     " )
ix=ix+1 
IC_AOE=ix
Country( IC_AOE ) = cc(  "AOE" , 86 ,F, 86, 9  , "Arctic Ocean (in the extended EMEP domain)                       " )

! New external areas (outside the 132x159 grid), these are normally not used
! a) Domains: x = 160-170 y = 1-132 and x =  -16-0  y = 123-170
ix=ix+1 
IC_RFX=ix
Country( IC_RFX ) = cc(  "RFX" , 87 ,F,  87, -100  ,"Extended EMEP-external part of Russian Federation" )
ix=ix+1 
IC_ASX=ix
Country( IC_ASX ) = cc(  "ASX" , 88 ,F,  88, -100   ,"Extended EMEP-external part of Asia              " )
ix=ix+1 
IC_PAX=ix
Country( IC_PAX ) = cc(  "PAX" , 89 ,F,  89, -100  ,"Extended EMEP-external part of Pacific Ocean     " )
ix=ix+1 
IC_AOX=ix
Country( IC_AOX ) = cc(  "AOX" , 90 ,F,  90, 9  ,"Extended EMEP-external part of Arctic Ocean      " )
! b) Domain x = -16-132 y = -11-0
ix=ix+1 
IC_NAX=ix
Country( IC_NAX ) = cc(  "NAX" , 91 ,F,  91, -100   ,"EMEP-external part of North Africa               " )
ix=ix+1 
IC_RUE=ix
Country( IC_RUE) = cc(  "RUE" , 93 ,F,  93, -100 , "Russian Federeation (all)   " )   

! Biomass burning
ix=ix+1 
IC_BB=ix
Country( IC_BB)  = cc(  "BB ", 101,F,  101, -100  , "Biomass burning (wild)        " )

! Sea areas split according to innside/outside 12 nautical mile zone, 
! ferries/cargo ships, registred inside/outside EU
ix=ix+1 
IC_BA2=ix
Country( IC_BA2 ) = cc(  "BA2" ,302 ,T,  30, 1  , "Baltic EU cargo outs.12      " )
ix=ix+1 
IC_BA3=ix
Country( IC_BA3 ) = cc(  "BA3" ,303 ,T,  30, 1  , "Baltic ROW cargo outs. 12    " )
ix=ix+1 
IC_BA4=ix
Country( IC_BA4 ) = cc(  "BA4" ,304 ,T,  30, 1  , "Baltic EU cargo ins. 12      " )
ix=ix+1 
IC_BA5=ix
Country( IC_BA5 ) = cc(  "BA5" ,305 ,T,  30, 1  , "Baltic ROW cargo ins. 12     " )
ix=ix+1 
IC_BA6=ix
Country( IC_BA6 ) = cc(  "BA6" ,306 ,T,  30, 1  , "Baltic EU ferries outs.12    " )
ix=ix+1 
IC_BA7=ix
Country( IC_BA7 ) = cc(  "BA7" ,307 ,T,  30, 1  , "Baltic ROW ferries outs. 12  " )
ix=ix+1 
IC_BA8=ix
Country( IC_BA8 ) = cc(  "BA8" ,308 ,T,  30, 1  , "Baltic EU ferries ins. 12    " )
ix=ix+1 
IC_BA9=ix
Country( IC_BA9 ) = cc(  "BA9" ,309 ,T,  30, 1  , "Baltic ROW ferries ins. 12   " )

ix=ix+1 
IC_NS2=ix
Country( IC_NS2 ) = cc(  "NS2" ,312 ,T,  31, 1  , "N. Sea EU cargo outs.12      " )
ix=ix+1 
IC_NS3=ix
Country( IC_NS3 ) = cc(  "NS3" ,313 ,T,  31, 1  , "N. Sea ROW cargo outs. 12    " )
ix=ix+1 
IC_NS4=ix
Country( IC_NS4 ) = cc(  "NS4" ,314 ,T,  31, 1  , "N. Sea EU cargo ins. 12      " )
ix=ix+1 
IC_NS5=ix
Country( IC_NS5 ) = cc(  "NS5" ,315 ,T,  31, 1  , "N. Sea ROW cargo ins. 12     " )
ix=ix+1 
IC_NS6=ix
Country( IC_NS6 ) = cc(  "NS6" ,316 ,T,  31, 1  , "N. Sea EU ferries outs.12    " )
ix=ix+1 
IC_NS7=ix
Country( IC_NS7 ) = cc(  "NS7" ,317 ,T,  31, 1  , "N. Sea ROW ferries outs. 12  " )
ix=ix+1 
IC_NS8=ix
Country( IC_NS8 ) = cc(  "NS8" ,318 ,T,  31, 1  , "N. Sea EU ferries ins. 12    " )
ix=ix+1 
IC_NS9=ix
Country( IC_NS9 ) = cc(  "NS9" ,319 ,T,  31, 1  , "N. Sea ROW ferries ins. 12   " )

ix=ix+1 
IC_AT2=ix
Country( IC_AT2 ) = cc(  "AT2" ,322 ,T,  32, 1  , "Atlant EU cargo outs.12      " )
ix=ix+1 
IC_AT3=ix
Country( IC_AT3 ) = cc(  "AT3" ,323 ,T,  32, 1  , "Atlant ROW cargo outs. 12    " )
ix=ix+1 
IC_AT4=ix
Country( IC_AT4 ) = cc(  "AT4" ,324 ,T,  32, 1  , "Atlant EU cargo ins. 12      " )
ix=ix+1 
IC_AT5=ix
Country( IC_AT5 ) = cc(  "AT5" ,325 ,T,  32, 1  , "Atlant ROW cargo ins. 12     " )
ix=ix+1 
IC_AT6=ix
Country( IC_AT6 ) = cc(  "AT6" ,326 ,T,  32, 1  , "Atlant EU ferries outs.12    " )
ix=ix+1 
IC_AT7=ix
Country( IC_AT7 ) = cc(  "AT7" ,327 ,T,  32, 1  , "Atlant ROW ferries outs. 12  " )
ix=ix+1 
IC_AT8=ix
Country( IC_AT8 ) = cc(  "AT8" ,328 ,T,  32, 1  , "Atlant EU ferries ins. 12    " )
ix=ix+1 
IC_AT9=ix
Country( IC_AT9 ) = cc(  "AT9" ,329 ,T,  32, 1  , "Atlant ROW ferries ins. 12   " )

ix=ix+1 
IC_ME2=ix
Country( IC_ME2 ) = cc(  "ME2" ,332 ,T,  33, 1  , "Medite EU cargo outs.12      " )
ix=ix+1 
IC_ME3=ix
Country( IC_ME3 ) = cc(  "ME3" ,333 ,T,  33, 1  , "Medite ROW cargo outs. 12    " )
ix=ix+1 
IC_ME4=ix
Country( IC_ME4 ) = cc(  "ME4" ,334 ,T,  33, 1  , "Medite EU cargo ins. 12      " )
ix=ix+1 
IC_ME5=ix
Country( IC_ME5 ) = cc(  "ME5" ,335 ,T,  33, 1  , "Medite ROW cargo ins. 12     " )
ix=ix+1 
IC_ME6=ix
Country( IC_ME6 ) = cc(  "ME6" ,336 ,T,  33, 1  , "Medite EU ferries outs.12    " )
ix=ix+1 
IC_ME7=ix
Country( IC_ME7 ) = cc(  "ME7" ,337 ,T,  33, 1  , "Medite ROW ferries outs. 12  " )
ix=ix+1 
IC_ME8=ix
Country( IC_ME8 ) = cc(  "ME8" ,338 ,T,  33, 1  , "Medite EU ferries ins. 12    " )
ix=ix+1 
IC_ME9=ix
Country( IC_ME9 ) = cc(  "ME9" ,339 ,T,  33, 1  , "Medite ROW ferries ins. 12   " )

ix=ix+1 
IC_BL2=ix
Country( IC_BL2 ) = cc(  "BL2" ,342 ,T,  34, 2  , "B. Sea EU cargo outs.12      " )
ix=ix+1 
IC_BL3=ix
Country( IC_BL3 ) = cc(  "BL3" ,343 ,T,  34, 2  , "B. Sea ROW cargo outs. 12    " )
ix=ix+1 
IC_BL4=ix
Country( IC_BL4 ) = cc(  "BL4" ,344 ,T,  34, 2  , "B. Sea EU cargo ins. 12      " )
ix=ix+1 
IC_BL5=ix
Country( IC_BL5 ) = cc(  "BL5" ,345 ,T,  34, 2  , "B. Sea ROW cargo ins. 12     " )
ix=ix+1 
IC_BL6=ix
Country( IC_BL6 ) = cc(  "BL6" ,346 ,T,  34, 2  , "B. Sea EU ferries outs.12    " )
ix=ix+1 
IC_BL7=ix
Country( IC_BL7 ) = cc(  "BL7" ,347 ,T,  34, 2  , "B. Sea ROW ferries outs. 12  " )
ix=ix+1 
IC_BL8=ix
Country( IC_BL8 ) = cc(  "BL8" ,348 ,T,  34, 2  , "B. Sea EU ferries ins. 12    " )
ix=ix+1 
IC_BL9=ix
Country( IC_BL9 ) = cc(  "BL9" ,349 ,T,  34, 2  , "B. Sea ROW ferries ins. 12   " )


! NH3Emis new land code for NMR-NH3 project
ix=ix+1 
IC_NMR=ix
Country( IC_NMR ) = cc(  "NMR" , 98 ,F, 98, 1  , "Area with temporal NMR-NH3 emissions              " )


!Extra cc for rest CityZen
ix=ix+1 
IC_RAA=ix
Country( IC_RAA ) = cc(  "RAA" , 170 ,F,  170, -100, "Rest of Africa and Asia" )
ix=ix+1 
IC_SEA=ix
Country( IC_SEA ) = cc(  "SEA" , 171 ,F,  171, -100, "Ships" )

! Extra from IIASA/ECLIPSE/ECLAIRE global

ix=ix+1 
IC_AFGH=ix
Country(IC_AFGH) = cc( "AFGH", 201, F,201, -100, "Afghanistan") 
ix=ix+1 
IC_ARGE=ix
Country(IC_ARGE) = cc( "ARGE", 202, F,202, -100, "Argentina") 
ix=ix+1 
IC_AUTR=ix
Country(IC_AUTR) = cc( "AUTR", 203, F,203, -100, "Australia") 
ix=ix+1 
IC_BANG=ix
Country(IC_BANG) = cc( "BANG", 204, F,204, -100, "Bangladesh") 
ix=ix+1 
IC_BHUT=ix
Country(IC_BHUT) = cc( "BHUT", 205, F,205, -100, "Bhutan") 
ix=ix+1 
IC_BRAZ=ix
Country(IC_BRAZ) = cc( "BRAZ", 206, F,206, -100, "Brazil") 
ix=ix+1 
IC_BRUN=ix
Country(IC_BRUN) = cc( "BRUN", 207, F,207, -100, "Brunei") 
ix=ix+1 
IC_CAMB=ix
Country(IC_CAMB) = cc( "CAMB", 208, F,208, -100, "Cambodia") 
ix=ix+1 
IC_CHIL=ix
Country(IC_CHIL) = cc( "CHIL", 209, F,209, -100, "Chile") 
ix=ix+1 
IC_CHIN=ix
Country(IC_CHIN) = cc( "CHIN", 210, F,210, -100, "China") 
ix=ix+1 
IC_FSUA=ix
Country(IC_FSUA) = cc( "FSUA", 211, F,211, -100, "Former_USSR_(Asia)_Tajikistan_Turkmenistan_Uzbekistan") 
ix=ix+1 
IC_INDI=ix
Country(IC_INDI) = cc( "INDI", 212, F,212, -100, "India") 
ix=ix+1 
IC_INDO=ix
Country(IC_INDO) = cc( "INDO", 213, F,213, -100, "Indonesia") 
ix=ix+1 
IC_ISRA=ix
Country(IC_ISRA) = cc( "ISRA", 214, F,214, -100, "Israel") 
ix=ix+1 
IC_JAPA=ix
Country(IC_JAPA) = cc( "JAPA", 215, F,215, -100, "Japan") 
ix=ix+1 
IC_LAOS=ix
Country(IC_LAOS) = cc( "LAOS", 216, F,216, -100, "Laos") 
ix=ix+1 
IC_MALA=ix
Country(IC_MALA) = cc( "MALA", 217, F,217, -100, "Malaysia") 
ix=ix+1 
IC_MEXI=ix
Country(IC_MEXI) = cc( "MEXI", 218, F,218, -100, "Mexico") 
ix=ix+1 
IC_MIDE=ix
Country(IC_MIDE) = cc( "MIDE", 219, F,219, -100, "Middle_East") 
ix=ix+1 
IC_MONG=ix
Country(IC_MONG) = cc( "MONG", 220, F,220, -100, "Mongolia") 
ix=ix+1 
IC_MYAN=ix
Country(IC_MYAN) = cc( "MYAN", 221, F,221, -100, "Myanmar") 
ix=ix+1 
IC_NEPA=ix
Country(IC_NEPA) = cc( "NEPA", 222, F,222, -100, "Nepal") 
ix=ix+1 
IC_NZEL=ix
Country(IC_NZEL) = cc( "NZEL", 223, F,223, -100, "New_Zealand") 
ix=ix+1 
IC_NAFR=ix
Country(IC_NAFR) = cc( "NAFR", 224, F,224, -100, "North_Africa_Libya_Tunisia_Algeria_Sudan_Morocco") 
ix=ix+1 
IC_KORN=ix
Country(IC_KORN) = cc( "KORN", 225, F,225, -100, "North_Korea") 
ix=ix+1 
IC_OAFR=ix
Country(IC_OAFR) = cc( "OAFR", 226, F,226, -100, "Other_Africa") 
ix=ix+1 
IC_OLAM=ix
Country(IC_OLAM) = cc( "OLAM", 227, F,227, -100, "Other_Latin_America") 
ix=ix+1 
IC_PAKI=ix
Country(IC_PAKI) = cc( "PAKI", 228, F,228, -100, "Pakistan") 
ix=ix+1 
IC_PHIL=ix
Country(IC_PHIL) = cc( "PHIL", 229, F,229, -100, "Philippines") 
ix=ix+1 
IC_SING=ix
Country(IC_SING) = cc( "SING", 230, F,230, -100, "Singapore") 
ix=ix+1 
IC_SAFR=ix
Country(IC_SAFR) = cc( "SAFR", 231, F,231, -100, "South_Africa") 
ix=ix+1 
IC_KORS=ix
Country(IC_KORS) = cc( "KORS", 232, F,232, -100, "South_Korea") 
ix=ix+1 
IC_SRIL=ix
Country(IC_SRIL) = cc( "SRIL", 233, F,233, -100, "Sri_Lanka") 
ix=ix+1 
IC_TAIW=ix
Country(IC_TAIW) = cc( "TAIW", 234, F,234, -100, "Taiwan") 
ix=ix+1 
IC_THAI=ix
Country(IC_THAI) = cc( "THAI", 235, F,235, -100, "Thailand") 
ix=ix+1 
IC_VIET=ix
Country(IC_VIET) = cc( "VIET", 236, F,236, -100, "Vietnam") 
ix=ix+1 
IC_EGYP=ix
Country(IC_EGYP) = cc( "EGYP", 237, F,237, -100, "Egypt")
ix=ix+1 
IC_INTSHIPS=ix
Country(IC_INTSHIPS ) = cc(  "INTSHIPS" ,350 ,T, 350, -100  , "International ships, RCP6" )

!!  HTAP2 regions
ix=ix+1 
IC_HT1000 = ix
Country(IC_HT1000 ) = cc(  "HT1000" ,1000 ,T, -100, 2  , "HT 1000" )
ix=ix+1 
IC_HTNATL = ix
Country(IC_HTNATL ) = cc(  "N_ATL" ,1002 ,T, 32, -100  , "Int. ships, N. Atl." )
ix=ix+1 
IC_HTUSCA = ix
Country(IC_HTUSCA ) = cc(  "USCA" ,1003 ,T, 65, -100  , "USA and Canada" )
ix=ix+1 
IC_HTEUTU = ix
Country(IC_HTEUTU ) = cc(  "EU_TU" ,1004 ,T, 64, 1  , "EU and Turkey" )
ix=ix+1 
IC_HTSASI = ix
Country(IC_HTSASI ) = cc(  "S_ASIA_IP" ,1005 ,T, 212, -100  , "S. Asia India Pak." )
ix=ix+1 
IC_HTEASI = ix
Country(IC_HTEASI ) = cc(  "CHCORJAP" ,1006 ,T, 210, -100  , "China Korea Japan" )
ix=ix+1 
IC_HTSEAS = ix
Country(IC_HTSEAS ) = cc(  "TAINDOMA" ,1007 ,T, 210, -100  , "Thail. Indon. Malay" )
ix=ix+1 
IC_HTAUST = ix
Country(IC_HTAUST ) = cc(  "AUSTR" ,1008 ,T, 203, -100  , "Austr N. Zeal. ++" )
ix=ix+1 
IC_HTNAFR = ix
Country(IC_HTNAFR ) = cc(  "NAFRI" ,1009 ,T, 64, 1  , "N. Africa" )
ix=ix+1 
IC_HTRAFR = ix
Country(IC_HTRAFR ) = cc(  "RAFRI" ,1010 ,T, 64, 1  , "Rest Africa" )
ix=ix+1 
IC_HTMIDE = ix
Country(IC_HTMIDE ) = cc(  "MIDEAST" ,1011 ,T, 25, 2  , "Middle East" )
ix=ix+1 
IC_HT1012 = ix
Country(IC_HT1012 ) = cc(  "HT1012" ,1012 ,T, -100, 2  , "HT 1012" )
ix=ix+1 
IC_HT1013 = ix
Country(IC_HT1013 ) = cc(  "HT1013" ,1013 ,T, -100, 2  , "HT 1013" )
ix=ix+1 
IC_HT1014 = ix
Country(IC_HT1014 ) = cc(  "HT1014" ,1014 ,T, -100, 2  , "HT 1014" )
ix=ix+1 
IC_HT1015 = ix
Country(IC_HT1015 ) = cc(  "HT1015" ,1015 ,T, -100, 2  , "HT 1015" )
ix=ix+1 
IC_HT1016 = ix
Country(IC_HT1016 ) = cc(  "HT1016" ,1016 ,T, -100, 2  , "HT 1016" )
ix=ix+1 
IC_HT1017 = ix
Country(IC_HT1017 ) = cc(  "HT1017" ,1017 ,T, -100, 2  , "HT 1017" )
ix=ix+1 
IC_HT1018 = ix
Country(IC_HT1018 ) = cc(  "HT1018" ,1018 ,T, -100, 2  , "HT 1018" )
ix=ix+1 
IC_HT1019 = ix
Country(IC_HT1019 ) = cc(  "HT1019" ,1019 ,T, -100, 2  , "HT 1019" )
ix=ix+1 
IC_HT1020 = ix
Country(IC_HT1020 ) = cc(  "HT1020" ,1020 ,T, -100, 2  , "HT 1020" )



NLAND=ix !actual number of countries defined

  end subroutine Country_Init


end module Country_ml
