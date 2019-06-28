! <Country_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2019 met.no
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
! <Country_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
module Country_mod

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
 ! timefac_index_hourly can be defined if the hourly timefactors are defined separately

  implicit none

  public :: init_Country     ! sets country details
  public :: self_test        ! just to test numbering

  integer, parameter, public  :: MAXNLAND = 601  ! max number of countries 
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


  !/ to be set in init_Country:

  type, public :: cc
     character(len=10) :: code          ! up to 10 letter land code
     character(len=4)  :: gains         ! 4 letter GAINS code
     integer           :: icode         ! integer number for land code (corresponds to 
                                        ! country code number in emission files)
     logical           :: is_sea        ! T for sea area, F otherwise
     integer           :: timefac_index ! Country code to use for timefactors
     integer           :: timefac_index_hourly !  Country code to use for hourly timefactors
     integer           :: timezone      ! timezone, deviation from UTC time
     character(len=60) :: name          ! name of country/region
  end type cc

  type(cc), public, save, dimension(MAXNLAND) :: Country

  integer, public ::  IC_AL   ! Albania                       
  integer, public ::  IC_AT   ! Austria                       
  integer, public ::  IC_BE   ! Belgium                       
  integer, public ::  IC_BG   ! Bulgaria                      
  integer, public :: IC_FCS   ! Former Czechoslovakia     
  integer, public ::  IC_DK   ! Denmark                       
  integer, public ::  IC_GL   ! Greenland                       
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
  integer, public :: IC_DUMMY  ! Generic or undefined country
  integer, public ::  IC_KG   ! Kyrgyzstan 
  integer, public ::  IC_AZ   ! Azerbaijan                 
  integer, public :: IC_ATX   ! ATL outside EMEP domain
  integer, public :: IC_RUX   ! RU outside old EMEP domain
  integer, public :: IC_RS    ! Serbia
  integer, public :: IC_ME    ! Montenegro

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
integer, public :: IC_HANO  ! Hanoi
integer, public :: IC_NVIE  ! North Vietnam
integer, public :: IC_SVIE  ! South Vietnam
integer, public :: IC_BOLV  ! Bolivia
integer, public :: IC_CARB  ! Caribbean
integer, public :: IC_CEAM  ! Central America
integer, public :: IC_COLO  ! Colombia
integer, public :: IC_ECUA  ! Ecuador
integer, public :: IC_PARA  ! Paraguay
integer, public :: IC_PERU  ! Peru
integer, public :: IC_URUG  ! Uruguay
integer, public :: IC_VENE  ! Venezuela
integer, public :: IC_IRAN  ! Iran
integer, public :: IC_SAAR  ! Saudi Arabia
integer, public :: IC_KOSO  ! Kosovo
integer, public :: IC_OCEC  ! Oceania

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

  ! CAMS-TNO sea regions (added Nov 2018)
  integer,  public :: IC_GRS  !Greenland Sea
  integer,  public :: IC_BAR  !Barents Sea
  integer,  public :: IC_ENC  !English Channel
  integer,  public :: IC_NWS  !Norwegian Sea
  integer,  public :: IC_IRC  !Irish Sea
  integer,  public :: IC_PSG  !Persian Gulf
  integer,  public :: IC_KAR  !Kara Sea 

  ! Ship emissions when sea areas are not divided
  ! Eg. TNO emissions (added on 25th March 2009)
  integer,  public :: IC_INTSHIPS
  ! Aircraft, used in AirEmis
  integer,  public :: IC_AIRCRAFT

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

! Divided countries put together
  integer, public :: IC_RUE   ! Russian Federation in the extended EMEP domain (RU+RFE+RUX, 36, 37, 38, 42, 71, 74)
  integer, public :: IC_KZT   ! Kazakhstan (KZ+KZE, 53, 75)
  integer, public :: IC_UZT   ! Uzbekistan (UZ+UZE, 76,78)
  integer, public :: IC_TMT   ! Turkmenistan (TM+TME, 77,79)

 !b) Domain x = -16-132 y = -11-0
  integer,  public :: IC_NAX  ! EMEP-external part of North Africa

 ! New code introduced for the NMR-NH3 project, not used in other projects
 ! NH3Emis x=44-75, y=35-66
  integer,  public :: IC_NMR  ! EMEP NMR-NH3 temporal emissions

 ! 199 not country-specific land based emissions -found in PANHAM/MEIC
  integer, public :: IC_LANDX  ! 199 not country-specific land based emissions

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
  integer,  public :: IC_HTCEAM
  integer,  public :: IC_HTSOAM
  integer,  public :: IC_HTRUBU
  integer,  public :: IC_HTCASI
  integer,  public :: IC_HTPOLA
  integer,  public :: IC_HTSPOL
  integer,  public :: IC_HT1018
  integer,  public :: IC_HT1019
  integer,  public :: IC_HT1020
  integer,  public :: IC_HT1000

! UNEP SR and new GAINS regions
  integer,  public :: IC_BANG_DHAK
  integer,  public :: IC_BANG_REST
  integer,  public :: IC_CHIN_ANHU
  integer,  public :: IC_CHIN_BEIJ
  integer,  public :: IC_CHIN_CHON
  integer,  public :: IC_CHIN_FUJI 
  integer,  public :: IC_CHIN_GANS
  integer,  public :: IC_CHIN_GUAD
  integer,  public :: IC_CHIN_GUAX
  integer,  public :: IC_CHIN_GUIZ
  integer,  public :: IC_CHIN_HAIN
  integer,  public :: IC_CHIN_HEBE  
  integer,  public :: IC_CHIN_HEIL
  integer,  public :: IC_CHIN_HENA
  integer,  public :: IC_CHIN_HONG
  integer,  public :: IC_CHIN_HUBE
  integer,  public :: IC_CHIN_HUNA
  integer,  public :: IC_CHIN_JILI
  integer,  public :: IC_CHIN_JINU
  integer,  public :: IC_CHIN_JINX
  integer,  public :: IC_CHIN_LIAO
  integer,  public :: IC_CHIN_NEMO
  integer,  public :: IC_CHIN_NINX
  integer,  public :: IC_CHIN_QING
  integer,  public :: IC_CHIN_SHAA
  integer,  public :: IC_CHIN_SHAN
  integer,  public :: IC_CHIN_SHND
  integer,  public :: IC_CHIN_SHNX
  integer,  public :: IC_CHIN_SICH
  integer,  public :: IC_CHIN_TIAN
  integer,  public :: IC_CHIN_TIBE
  integer,  public :: IC_CHIN_XING
  integer,  public :: IC_CHIN_YUNN
  integer,  public :: IC_CHIN_ZHEJ
  integer,  public :: IC_INDI_ANPR
  integer,  public :: IC_INDI_ASSA
  integer,  public :: IC_INDI_BENG
  integer,  public :: IC_INDI_BIHA
  integer,  public :: IC_INDI_CHHA
  integer,  public :: IC_INDI_DELH
  integer,  public :: IC_INDI_EHIM 
  integer,  public :: IC_INDI_GOA
  integer,  public :: IC_INDI_GUJA
  integer,  public :: IC_INDI_HARY
  integer,  public :: IC_INDI_HIPR
  integer,  public :: IC_INDI_JHAR
  integer,  public :: IC_INDI_KARN
  integer,  public :: IC_INDI_KERA
  integer,  public :: IC_INDI_MAHA
  integer,  public :: IC_INDI_MAPR
  integer,  public :: IC_INDI_ORIS
  integer,  public :: IC_INDI_PUNJ
  integer,  public :: IC_INDI_RAJA
  integer,  public :: IC_INDI_TAMI
  integer,  public :: IC_INDI_UTAN
  integer,  public :: IC_INDI_UTPR
  integer,  public :: IC_INDI_WHIM
  integer,  public :: IC_INDO_JAKA 
  integer,  public :: IC_INDO_JAVA 
  integer,  public :: IC_INDO_REST 
  integer,  public :: IC_INDO_SUMA
  integer,  public :: IC_JAPA_CHSH
  integer,  public :: IC_JAPA_CHUB
  integer,  public :: IC_JAPA_HOTO
  integer,  public :: IC_JAPA_KANT
  integer,  public :: IC_JAPA_KINK
  integer,  public :: IC_JAPA_KYOK
  integer,  public :: IC_KORS_NORT
  integer,  public :: IC_KORS_PUSA
  integer,  public :: IC_KORS_SEOI
  integer,  public :: IC_KORS_SOUT
  integer,  public :: IC_MALA_KUAL
  integer,  public :: IC_MALA_PENM
  integer,  public :: IC_MALA_SASA
  integer,  public :: IC_PAKI_KARA
  integer,  public :: IC_PAKI_NMWP
  integer,  public :: IC_PAKI_PUNJ
  integer,  public :: IC_PAKI_SIND
  integer,  public :: IC_PHIL_BVMI
  integer,  public :: IC_PHIL_LUZO
  integer,  public :: IC_PHIL_MANI
  integer,  public :: IC_RUSS_ASIA
  integer,  public :: IC_RUSS_EURO
  integer,  public :: IC_THAI_BANG
  integer,  public :: IC_THAI_CVAL
  integer,  public :: IC_THAI_NEPL
  integer,  public :: IC_THAI_NHIG
  integer,  public :: IC_THAI_SPEN
  integer,  public :: IC_IRAN_REST
  integer,  public :: IC_IRAN_TEHR
  integer,  public :: IC_VIET_BNIH
  integer,  public :: IC_VIET_HYEN
  integer,  public :: IC_VIET_OCAR
  integer,  public :: IC_VIET_ONOR
  integer,  public :: IC_USAM_ALAS
  integer,  public :: IC_USAM_REST
  integer,  public :: IC_EGYP_CAIR
  integer,  public :: IC_EGYP_REST
  integer,  public :: IC_NAFR_WHOL
  integer,  public :: IC_NIGE_LAGO
  integer,  public :: IC_NIGE_OLAG
  integer,  public :: IC_NIGE_REST
  integer,  public :: IC_SAFR_JOBG
  integer,  public :: IC_SAFR_REST
  integer,  public :: IC_RSAF_WHOL
  integer,  public :: IC_EAFR_WHOL
  integer,  public :: IC_WAFR_WHOL
  integer,  public :: IC_KENY_WHOL
  integer,  public :: IC_TANZ_WHOL
  integer,  public :: IC_NIGE_WHOL 
 
  integer,  public :: IC_BERLIN
  integer,  public :: IC_BRUSSEL
  integer,  public :: IC_COPENHAGEN
  integer,  public :: IC_MADRID
  integer,  public :: IC_HELSINKI
  integer,  public :: IC_PARIS
  integer,  public :: IC_ATHENS
  integer,  public :: IC_BUDAPEST
  integer,  public :: IC_DUBLIN
  integer,  public :: IC_MILAN
  integer,  public :: IC_OSLO
  integer,  public :: IC_AMSTERDAM
  integer,  public :: IC_WARSAW
  integer,  public :: IC_LISBON
  integer,  public :: IC_LONDON
  integer,  public :: IC_STOCKHOLM
  integer,  public :: IC_GENEVA

 contains
  
  subroutine init_Country()

  ! Set the country details. Note that time-zones for some areas are either
  ! difficult (e.g. Russia should be 3 to 12) or not relevant (e.g. sea areas,
  ! volcanoes). 

  integer :: iland,ix

  ! First define all countries as undefined
  do iland=1,NLAND
    Country(iland) = cc(  "N/A" ,'-', iland ,F,  17 , 17 , -100  , "Not_defined                   " )
  end do

!The value of IC_XX is the index in Country array. Can in principle change between two runs or versions.
!The emission_code is the country code used in the emission file.
!The timefac_code is the code/index refering to the timefactor 

  !--------------  code emission_icode sea timefac_code timefac_code_hourly timezone  Name  ------------!
  
ix=0 
ix=ix+1 
IC_AL=ix
Country( IC_AL ) = cc(  "AL ",'ALBA', 1 ,F,  1,  1,  1  , "Albania                       " )
ix=ix+1 
IC_AT=ix
Country( IC_AT ) = cc(  "AT ",'AUST', 2 ,F,  2,  2,  1  , "Austria                       " )
ix=ix+1 
IC_BE=ix
Country( IC_BE ) = cc(  "BE ",'BELG', 3 ,F,  3,  3,  1  , "Belgium                       " )
ix=ix+1 
IC_BG=ix
Country( IC_BG ) = cc(  "BG ",'BULG', 4 ,F,  4,  4,  2  , "Bulgaria                      " )
ix=ix+1 
IC_FCS=ix
Country( IC_FCS ) = cc( "FCS ",'-', 5 ,F,  5,  5,  1 , "Former Czechoslovakia         " )
ix=ix+1 
IC_DK=ix
Country( IC_DK ) = cc(  "DK ",'DENM', 6 ,F,  6,  6,  1  , "Denmark                       " )
ix=ix+1 
IC_GL=ix
Country( IC_GL ) = cc(  "GL ",'-' , 601 ,F,  6,  6, -2  , "Greenland                     " )
ix=ix+1 
IC_FI=ix
Country( IC_FI ) = cc(  "FI ",'FINL', 7 ,F,  7,   7,  2  , "Finland                       " )
ix=ix+1 
IC_FR=ix
Country( IC_FR ) = cc(  "FR ",'FRAN',  8 ,F,  8,  8,  1  , "France                        " )
ix=ix+1 
IC_GDR=ix
Country( IC_GDR) = cc(  "GDR",'-',  9 ,F,  9,  9,  1  , "Former East Germany           " )
ix=ix+1 
IC_FRG=ix
Country( IC_FRG) = cc(  "FRG",'-', 10 ,F, 10,  10,  1  , "Former Fed. Rep. of Germany   " )
ix=ix+1 
IC_GR=ix
Country( IC_GR ) = cc(  "GR ",'GREE', 11 ,F, 11,  11,  2  , "Greece                        " )
ix=ix+1 
IC_HU=ix
Country( IC_HU ) = cc(  "HU ",'HUNG', 12 ,F, 12,  12,  1  , "Hungary                       " )
ix=ix+1 
IC_IS=ix
Country( IC_IS ) = cc(  "IS ",'ICEL', 13 ,F, 13,  13, 0  , "Iceland                       " )
ix=ix+1 
IC_IE=ix
Country( IC_IE ) = cc(  "IE ",'IREL', 14 ,F, 14,  14,  0  , "Ireland                       " )
ix=ix+1 
IC_IT=ix
Country( IC_IT ) = cc(  "IT ",'ITAL', 15 ,F, 15,   15,  1  , "Italy                         " )
ix=ix+1 
IC_LU=ix
Country( IC_LU ) = cc(  "LU ",'LUXE', 16 ,F, 16,  16,  1  , "Luxembourg                    " )
ix=ix+1 
IC_NL=ix
Country( IC_NL ) = cc(  "NL ",'NETH', 17 ,F, 17,  17,  1  , "Netherlands                   " )
ix=ix+1 
IC_NO=ix
Country( IC_NO ) = cc(  "NO ",'NORW', 18 ,F, 18,  18,  1  , "Norway                        " )
ix=ix+1 
IC_PL=ix
Country( IC_PL ) = cc(  "PL ",'POLA', 19 ,F, 19,  19,  1  , "Poland                        " )
ix=ix+1 
IC_PT=ix
Country( IC_PT ) = cc(  "PT ",'PORT', 20 ,F, 20,  20,  0  , "Portugal                      " )
ix=ix+1 
IC_RO=ix
Country( IC_RO ) = cc(  "RO ",'ROMA', 21 ,F, 21,  21,  2  , "Romania                       " )
ix=ix+1 
IC_ES =ix
Country( IC_ES ) = cc(  "ES ",'SPAI', 22 ,F, 22,  22,  1  , "Spain                         " )
ix=ix+1 
IC_SE=ix
Country( IC_SE ) = cc(  "SE ",'SWED', 23 ,F, 23,  23,  1  , "Sweden                        " )
ix=ix+1 
IC_CH=ix
Country( IC_CH ) = cc(  "CH ",'SWIT', 24 ,F, 24, 24,  1  , "Switzerland                   " )
ix=ix+1 
IC_TR=ix
Country( IC_TR ) = cc(  "TR ",'TURK', 25 ,F, 25,  25,  2  , "Turkey                        " )
ix=ix+1 
IC_SU=ix
Country( IC_SU ) = cc(  "SU ",'-', 26 ,F, 26, 26,  -100  , "Former USSR                   " )
ix=ix+1 
IC_GB=ix
Country( IC_GB ) = cc(  "GB " ,'UNKI', 27 ,F, 27,  27,  0  , "United Kingdom                " )
ix=ix+1 
IC_VUL=ix
Country( IC_VUL) = cc(  "VUL" ,'-', 28 ,F, 28,  28,  1  , "Volcanoes                     " )
ix=ix+1 
IC_REM=ix
Country( IC_REM) = cc(  "REM" ,'-', 29 ,F, 29,  29,  1  , "Remaining land areas          " )
!NB:
!Fix needed for following sea-areas (BAS,'-',NOS,ATL,MED,BLS)in GEA runs done in Emissions_mod
!if ( GRID == "HIRHAM" .and. IIFULLDOM == 182 ) then ! Special fix for HIRHAM/GEA
!if ( SEAFIX_GEA_NEEDED ) then ! Special fix for HIRHAM/GEA
ix=ix+1 
IC_BAS=ix
Country( IC_BAS) = cc(  "BAS" ,'-', 30 ,T, 30,  30,  1  , "The Baltic Sea                " )
ix=ix+1 
IC_NOS=ix
Country( IC_NOS) = cc(  "NOS" ,'-', 31 ,T, 31,  31,  1  , "The North Sea                 " )
ix=ix+1 
IC_ATL=ix
Country( IC_ATL) = cc(  "ATL" ,'-', 32 ,T, 32,  32,  -100  , "Remaining NE Atlantic Ocean   " )
ix=ix+1 
IC_MED=ix
Country( IC_MED) = cc(  "MED" ,'-', 33 ,T, 33,  33,  1  , "The Mediterranean Sea         " )
ix=ix+1 
IC_BLS=ix
Country( IC_BLS) = cc(  "BLS" ,'-', 34 ,T, 34,  34,  2  , "The Black Sea                 " )

!end if ! HIRHAM/GEA fix

ix=ix+1 
IC_NAT=ix
Country( IC_NAT) = cc(  "NAT",'-', 35 ,F, 35,  35,  -100  , "Natural marine sources        " )
ix=ix+1
IC_RUO=ix
Country( IC_RUO) = cc(  "RUO",'KOLK', 36 ,F, 36,  36,  4  , "Kola/Karelia                  " )
ix=ix+1 
IC_RUP=ix
! Not sure about GAINS code:
Country( IC_RUP) = cc(  "RUP",'RUSS', 37 ,F, 37,   37,  4  , "St.Petersburg/Novgorod-Pskov  " )
ix=ix+1 
IC_RUA=ix
Country( IC_RUA) = cc(  "RUA",'KALI', 38 ,F, 38, 38,  3  , "Kaliningrad                   " )
ix=ix+1 
IC_BY=ix
Country( IC_BY ) = cc(  "BY ",'BELA', 39 ,F, 39, 39,  3  , "Belarus                       " )
ix=ix+1 
IC_UA=ix
Country( IC_UA ) = cc(  "UA ",'UKRA', 40 ,F, 40, 40,  2  , "Ukraine                       " )
ix=ix+1 
IC_MD=ix
Country( IC_MD ) = cc(  "MD ",'MOLD', 41 ,F, 41, 41,  2  , "Moldova, Republic of          " )
ix=ix+1 
IC_RUR=ix
!Could also be REMR for GAINS
Country( IC_RUR) = cc(  "RUR",'RUSS', 42 ,F, 42, 42,  4  , "Rest of Russia                " )
ix=ix+1 
IC_EE=ix
Country( IC_EE ) = cc(  "EE ",'ESTO', 43 ,F, 43, 43,  2  , "Estonia                       " )
ix=ix+1 
IC_LV=ix
Country( IC_LV ) = cc(  "LV ",'LATV', 44 ,F, 44,  44,  2  , "Latvia                        " )
ix=ix+1 
IC_LT=ix
Country( IC_LT ) = cc(  "LT ",'LITH', 45 ,F, 45, 45,  2  , "Lithuania                     " )
ix=ix+1 
IC_CZ=ix
Country( IC_CZ ) = cc(  "CZ ",'CZRE', 46 ,F, 46, 46,  1  , "Czech                         " )
ix=ix+1 
IC_SK=ix
Country( IC_SK ) = cc(  "SK ",'SKRE', 47 ,F, 47, 47,  1  , "Slovakia                      " )
ix=ix+1 
IC_SI=ix
Country( IC_SI ) = cc(  "SI ",'SLOV', 48 ,F, 48,  48,  1  , "Slovenia                      " )
ix=ix+1 
IC_HR=ix
Country( IC_HR ) = cc(  "HR ",'CROA', 49 ,F, 49,  49,  1  , "Croatia                       " )
ix=ix+1 
IC_BA=ix
Country( IC_BA ) = cc(  "BA ",'BOHE', 50 ,F, 50,  50,  1  , "Bosnia and Herzegovina        " )
ix=ix+1 
IC_CS=ix
Country( IC_CS ) = cc(  "CS ",'SEMO', 51 ,F, 51, 51,  1  , "Serbia and Montenegro    " )
ix=ix+1 
IC_MK=ix
Country( IC_MK ) = cc(  "MK ",'MACE', 52 ,F, 52,  52,  1  , "Macedonia, The F.Yugo.Rep. of " )
ix=ix+1 
IC_KZ=ix
Country( IC_KZ ) = cc(  "KZ ",'KAZA', 53 ,F, 53, 53,  -100  , "Kazakstan                     " )
ix=ix+1
IC_GE=ix
Country( IC_GE ) = cc(  "GE ",'GEOR', 54 ,F, 54, 54,  4  , "Georgia                       " )
ix=ix+1 
IC_CY=ix
Country( IC_CY ) = cc(  "CY ",'CYPR', 55 ,F, 55, 55,  2  , "Cyprus                        " )
ix=ix+1 
IC_AM =ix
Country( IC_AM ) = cc(  "AM ",'ARME', 56 ,F, 56, 56,  4  , "Armenia                       " )
ix=ix+1 
IC_MT=ix
Country( IC_MT ) = cc(  "MT ",'MALT', 57 ,F, 57, 57,  1  , "Malta                         " )
ix=ix+1 
IC_ASI=ix
Country( IC_ASI) = cc(  "ASI" ,'-', 58 ,F, 58, 58,  -100  , "Other Asian areas             " )
ix=ix+1 
IC_LI=ix
Country( IC_LI ) = cc(  "LI " ,'-', 59 ,F, 59, 59,  1  , "Lichtenstein                  " )
ix=ix+1 
IC_DE=ix
Country( IC_DE ) = cc(  "DE ",'GERM', 60 ,F, 60, 60,  1  , "Germany                       " )
ix=ix+1 
IC_RU=ix
Country( IC_RU ) = cc(  "RU " ,'RUSS', 61 ,F, 61, 61,  -100  , "Russian Federation            " )
ix=ix+1 
IC_MC=ix
Country( IC_MC ) = cc(  "MC " ,'-', 62 ,F, 62, 62,  1  , "Monaco                        " )
ix=ix+1 
IC_NOA=ix
Country( IC_NOA) = cc(  "NOA" ,'-', 63 ,F, 63, 63,  1  , "North Africa                  " )
ix=ix+1 
IC_EU=ix
Country( IC_EU ) = cc(  "EU " ,'-', 64 ,F, 64, 64,  1  , "European Community            " )
ix=ix+1 
IC_US=ix
Country( IC_US ) = cc(  "US " ,'-', 65 ,F, 65,65,  -100  , "USA                           " )
ix=ix+1 
IC_CA=ix
Country( IC_CA ) = cc(  "CA " ,'-', 66 ,F, 66, 66,  -100  , "Canada                        " )
ix=ix+1 
IC_DUMMY=ix
Country( IC_DUMMY ) &
                 = cc(  "N/A" ,'-', 67 ,F,  67, 67, -100  , "Not_defined                   " )
ix=ix+1
IC_KG=ix
Country( IC_KG ) = cc(  "KG " ,'-', 68 ,F,  68, 68, 6  , "Kyrgyzstan                    " )
ix=ix+1 
IC_AZ=ix
Country( IC_AZ ) = cc(  "AZ " ,'-', 69 ,F,  69, 69, 4  , "Azerbaijan                    " )
ix=ix+1 
IC_ATX=ix
Country( IC_ATX) = cc(  "ATX" ,'-', 70 ,T,  32, 32, -100  , "Atlantic outside. EMEP        " )
ix=ix+1
IC_RUX=ix
Country( IC_RUX) = cc(  "RUX" ,'-', 71 ,F,  42, 42, -100  , "Russian Fed. outside emep     " )
ix=ix+1 
IC_RS=ix
Country( IC_RS)  = cc(  "RS " ,'-', 72 ,F,  72, 72, 1  , "Serbia                        " )
ix=ix+1
IC_ME=ix
Country( IC_ME)  = cc(  "ME " ,'-', 73 ,F,  73, 73, 1  , "Montenegro                    " )
! Extended EMEP domain
ix=ix+1
IC_RFE=ix
Country( IC_RFE ) = cc(  "RFE" ,'-', 74 ,F, 74, 74, -100  , "Rest of extended Russian Federation (in the extended EMEP domain)" )
ix=ix+1 
IC_KZE=ix
Country( IC_KZE ) = cc(  "KZE" ,'-', 75 ,F, 75, 75, -100  , "Rest of Kazakhstan (in the extended EMEP domain)                 " )
ix=ix+1 
IC_UZ=ix
Country( IC_UZ  ) = cc(  "UZ"  ,'-', 76 ,F, 76, 76, -100  , "Uzbekistan (in the original EMEP domain)                         " )
ix=ix+1 
IC_TM=ix
Country( IC_TM  ) = cc(  "TM"  ,'-', 77 ,F, 77, 77, -100  , "Turkmenistan (in the original EMEP domain)                       " )
ix=ix+1 
IC_UZE=ix
Country( IC_UZE ) = cc(  "UZE" ,'-', 78 ,F, 78, 78, -100  , "Rest of  Uzbekistan (in the extended EMEP domain)                " )
ix=ix+1 
IC_TME=ix
Country( IC_TME ) = cc(  "TME" ,'-', 79 ,F, 79, 79, -100  , "Rest of Turkmenistan (in the extended EMEP domain)               " )
ix=ix+1 
IC_CAS=ix
Country( IC_CAS ) = cc(  "CAS" ,'-', 80 ,F, 80, 80, -100  , "Caspian Sea (in the original EMEP domain)                        " )
ix=ix+1 
IC_TJ=ix
Country( IC_TJ  ) = cc(  "TJ"  ,'-', 81 ,F, 81, 81, -100   ,"Tajikistan (in the extended EMEP domain)                         " )
ix=ix+1 
IC_ARL=ix
Country( IC_ARL ) = cc(  "ARL" ,'-', 82 ,F, 82, 82, -100  , "Aral Lake (in the original EMEP domain)                          " )
ix=ix+1 
IC_ARE=ix
Country( IC_ARE ) = cc(  "ARE" ,'-', 83 ,F, 83, 83, -100  , "Rest of Aral Lake (in the extended EMEP domain)                  " )
ix=ix+1 
IC_ASM=ix
Country( IC_ASM ) = cc(  "ASM" ,'-', 84 ,F, 84, 84, -100  , "Modified remaining Asian areas (in the original EMEP domain)     " )
ix=ix+1 
IC_ASE=ix
Country( IC_ASE ) = cc(  "ASE" ,'-', 85 ,F, 85, 85, -100  , "Remaining extended Asian areas (in the extended EMEP domain)     " )
ix=ix+1 
IC_AOE=ix
Country( IC_AOE ) = cc(  "AOE" ,'-', 86 ,F, 86, 86, -100  , "Arctic Ocean (in the extended EMEP domain)                       " )

! New external areas (outside the 132x159 grid),'-', these are normally not used
! a) Domains: x = 160-170 y = 1-132 and x =  -16-0  y = 123-170
ix=ix+1 
IC_RFX=ix
Country( IC_RFX ) = cc(  "RFX" ,'-', 87 ,F,  87, 87, -100  ,"Extended EMEP-external part of Russian Federation" )
ix=ix+1 
IC_ASX=ix
Country( IC_ASX ) = cc(  "ASX" ,'-', 88 ,F,  88, 88, -100   ,"Extended EMEP-external part of Asia " )
ix=ix+1 
IC_PAX=ix
Country( IC_PAX ) = cc(  "PAX" ,'-', 89 ,F,  89, 89, -100  ,"Extended EMEP-external part of Pacific Ocean " )
ix=ix+1 
IC_AOX=ix
Country( IC_AOX ) = cc(  "AOX" ,'-', 90 ,F,  90, 90, 9  ,"Extended EMEP-external part of Arctic Ocean " )
! b) Domain x = -16-132 y = -11-0 (never used)
ix=ix+1 
IC_NAX=ix
Country( IC_NAX ) = cc(  "NAX" ,'-', 91 ,F, 91, 91, -100   ,"EMEP-external part of North Africa " )
ix=ix+1 
IC_KZT=ix
Country( IC_KZT ) = cc(  "KZT" ,'-', 92 ,F, 92, 92, -100  , "Kazakhstan (all)" )
ix=ix+1 
IC_RUE=ix
Country( IC_RUE ) = cc(  "RUE" ,'-', 93 ,F,  93, 93, -100 , "Russian Federeation (all)" )   
ix=ix+1 
IC_UZT=ix
Country( IC_UZT ) = cc(  "UZT" ,'-', 94 ,F, 94, 94, -100  , "Uzbekistan (all)" )
ix=ix+1 
IC_TMT=ix
Country( IC_TMT ) = cc(  "TMT" ,'-', 95 ,F, 95, 95, -100  , "Turkmenistan  (all)" )

! NH3Emis new land code for NMR-NH3 project
ix=ix+1 
IC_NMR=ix
Country( IC_NMR ) = cc(  "NMR" ,'-', 98 ,F, 98, 98, 1  , "Area with temporal NMR-NH3 emissions              " )

! Biomass burning
ix=ix+1 
IC_BB=ix
Country( IC_BB)  = cc(  "BB ",'-', 101,F,  101, 101, -100  , "Biomass burning (wild)        " )



!Extra cc for rest CityZen
ix=ix+1 
IC_RAA=ix
Country( IC_RAA ) = cc(  "RAA" ,'-', 170 ,F,  170, 170, -100, "Rest of Africa and Asia" )
ix=ix+1 
IC_SEA=ix
Country( IC_SEA ) = cc(  "SEA" ,'-', 171 ,F,  171, 171,-100, "Ships" )

! Extra from IIASA/ECLIPSE/ECLAIRE global
!
ix=ix+1
IC_LANDX=ix
Country( IC_LANDX)  = cc(  "LANDX ",'-', 199,F,  199, 199, -100  , "not country-specific land based emissions" )

ix=ix+1 
IC_AFGH=ix
Country(IC_AFGH) = cc( "AFGH",'-', 201, F,201, 201, -100, "Afghanistan") 
ix=ix+1 
IC_ARGE=ix
Country(IC_ARGE) = cc( "ARGE",'-', 202, F,202, 202, -100, "Argentina") 
ix=ix+1 
IC_AUTR=ix
Country(IC_AUTR) = cc( "AUTR",'-', 203, F,203, 203, -100, "Australia") 
ix=ix+1 
IC_BANG=ix
Country(IC_BANG) = cc( "BANG",'-', 204, F,204, 204, -100, "Bangladesh") 
ix=ix+1 
IC_BHUT=ix
Country(IC_BHUT) = cc( "BHUT",'-', 205, F,205, 205, -100, "Bhutan") 
ix=ix+1 
IC_BRAZ=ix
Country(IC_BRAZ) = cc( "BRAZ",'-', 206, F,206, 206, -100, "Brazil") 
ix=ix+1 
IC_BRUN=ix
Country(IC_BRUN) = cc( "BRUN",'-', 207, F,207, 207, -100, "Brunei") 
ix=ix+1 
IC_CAMB=ix
Country(IC_CAMB) = cc( "CAMB",'-', 208, F,208, 208,-100, "Cambodia") 
ix=ix+1 
IC_CHIL=ix
Country(IC_CHIL) = cc( "CHIL",'-', 209, F,209, 209, -100, "Chile") 
ix=ix+1 
IC_CHIN=ix
Country(IC_CHIN) = cc( "CHIN",'-', 210, F,210, 210, -100, "China") 
ix=ix+1 
IC_FSUA=ix
Country(IC_FSUA) = cc( "FSUA",'-', 211, F,211, 211, -100, "Former_USSR_(Asia)_Tajikistan_Turkmenistan_Uzbekistan") 
ix=ix+1 
IC_INDI=ix
Country(IC_INDI) = cc( "INDI",'-', 212, F,212, 212, -100, "India") 
ix=ix+1 
IC_INDO=ix
Country(IC_INDO) = cc( "INDO",'-', 213, F,213, 213, -100, "Indonesia") 
ix=ix+1 
IC_ISRA=ix
Country(IC_ISRA) = cc( "ISRA",'-', 214, F,214, 214, -100, "Israel") 
ix=ix+1 
IC_JAPA=ix
Country(IC_JAPA) = cc( "JAPA",'-', 215, F,215, 215, -100, "Japan") 
ix=ix+1 
IC_LAOS=ix
Country(IC_LAOS) = cc( "LAOS",'-', 216, F,216, 216,-100, "Laos") 
ix=ix+1 
IC_MALA=ix
Country(IC_MALA) = cc( "MALA",'-', 217, F,217, 217, -100, "Malaysia") 
ix=ix+1 
IC_MEXI=ix
Country(IC_MEXI) = cc( "MEXI",'-', 218, F,218, 218, -100, "Mexico") 
ix=ix+1 
IC_MIDE=ix
Country(IC_MIDE) = cc( "MIDE",'-', 219, F,219, 219, -100, "Middle_East") 
ix=ix+1 
IC_MONG=ix
Country(IC_MONG) = cc( "MONG",'-', 220, F,220, 220, -100, "Mongolia") 
ix=ix+1 
IC_MYAN=ix
Country(IC_MYAN) = cc( "MYAN",'-', 221, F,221, 221, -100, "Myanmar") 
ix=ix+1 
IC_NEPA=ix
Country(IC_NEPA) = cc( "NEPA",'-', 222, F,222, 222, -100, "Nepal") 
ix=ix+1 
IC_NZEL=ix
Country(IC_NZEL) = cc( "NZEL",'-', 223, F,223, 223, -100, "New_Zealand") 
ix=ix+1 
IC_NAFR=ix
Country(IC_NAFR) = cc( "NAFR",'-', 224, F,224, 224, -100, "North_Africa_Libya_Tunisia_Algeria_Sudan_Morocco") 
ix=ix+1 
IC_KORN=ix
Country(IC_KORN) = cc( "KORN",'-', 225, F,225, 225, -100, "North_Korea") 
ix=ix+1 
IC_OAFR=ix
Country(IC_OAFR) = cc( "OAFR",'-', 226, F,226,226, -100, "Other_Africa") 
ix=ix+1 
IC_OLAM=ix
Country(IC_OLAM) = cc( "OLAM",'-', 227, F,227, 227, -100, "Other_Latin_America") 
ix=ix+1 
IC_PAKI=ix
Country(IC_PAKI) = cc( "PAKI",'-', 228, F,228, 228, -100, "Pakistan") 
ix=ix+1 
IC_PHIL=ix
Country(IC_PHIL) = cc( "PHIL",'-', 229, F,229,229, -100, "Philippines") 
ix=ix+1 
IC_SING=ix
Country(IC_SING) = cc( "SING",'-', 230, F,230,230, -100, "Singapore") 
ix=ix+1 
IC_SAFR=ix
Country(IC_SAFR) = cc( "SAFR",'-', 231, F,231,231, -100, "South_Africa") 
ix=ix+1 
IC_KORS=ix
Country(IC_KORS) = cc( "KORS",'-', 232, F,232, 232, -100, "South_Korea") 
ix=ix+1 
IC_SRIL=ix
Country(IC_SRIL) = cc( "SRIL",'-', 233, F,233, 233, -100, "Sri_Lanka") 
ix=ix+1 
IC_TAIW=ix
Country(IC_TAIW) = cc( "TAIW",'-', 234, F,234,234, -100, "Taiwan") 
ix=ix+1 
IC_THAI=ix
Country(IC_THAI) = cc( "THAI",'-', 235, F,235, 235, -100, "Thailand") 
ix=ix+1 
IC_VIET=ix
Country(IC_VIET) = cc( "VIET",'-', 236, F,236, 236, -100, "Vietnam") 
ix=ix+1 
IC_EGYP=ix
Country(IC_EGYP) = cc( "EGYP",'-', 237, F,237, 237, -100, "Egypt")
ix=ix+1
IC_HANO=ix
Country(IC_HANO) = cc( "Hanoi",'-', 238, F, 238,238, -100, "Hanoi")
ix=ix+1
IC_NVIE=ix
Country(IC_NVIE) = cc( "NVIET",'-', 239, F, 239,239, -100, "North Vietnam")
ix=ix+1
IC_SVIE=ix
Country(IC_SVIE) = cc( "SVIET",'-', 240, F, 240,240, -100, "South Vietnam")
ix=ix+1
IC_BOLV=ix
Country(IC_BOLV) = cc( "BOLV",'-', 241, F, 241,241, -100, "Bolivia")
ix=ix+1
IC_CARB=ix
Country(IC_CARB) = cc( "CARB",'-', 242, F, 242,242, -100, "Caribbean")
ix=ix+1
IC_CEAM=ix
Country(IC_CEAM) = cc( "CEAM",'-', 243, F, 243,243, -100, "Central America")
ix=ix+1
IC_COLO=ix
Country(IC_COLO) = cc( "COLO",'-', 244, F, 244, 244, -100, "Colombia")
ix=ix+1
IC_ECUA=ix
Country(IC_ECUA) = cc( "ECUA",'-', 245, F, 245,245, -100, "Ecuador")
ix=ix+1
IC_PARA=ix
Country(IC_PARA) = cc( "PARA",'-', 246, F, 246,246, -100, "Paraguay")
ix=ix+1
IC_PERU=ix
Country(IC_PERU) = cc( "PERU",'-', 247, F, 247, 247, -100, "Peru")
ix=ix+1
IC_URUG=ix
Country(IC_URUG) = cc( "URUG",'-', 248, F, 248, 248, -100, "Uruguay")
ix=ix+1
IC_VENE=ix
Country(IC_VENE) = cc( "VENE",'-', 249, F, 249,249, -100, " Venezuela")
ix=ix+1
IC_IRAN=ix
Country(IC_IRAN) = cc( "IRAN",'-', 250, F, 250, 250,-100, "Iran")
ix=ix+1
IC_SAAR=ix
Country(IC_SAAR) = cc( "SAAR",'-', 251, F, 251,251, -100, "Saudi Arabia")
ix=ix+1 
IC_INTSHIPS=ix
Country(IC_INTSHIPS) = cc( "INTSHIPS" ,'-',350 ,T, 350,350, -100  , "International ships" )
ix=ix+1 
IC_AIRCRAFT=ix
Country(IC_AIRCRAFT) = cc( "AIRCRAFT" ,'-',900 ,T, 900,900, -100  , "International Flights" )
ix=ix+1
IC_KOSO=ix
Country(IC_KOSO) = cc( "KOSO",'KOSO', 373, F, 373, 373, -100, "Kosovo")
ix=ix+1
IC_OCEC=ix
Country(IC_OCEC) = cc( "OCEC",'-', 393, F, 393,393, -100, "Oceania")


! Sea areas split according to innside/outside 12 nautical mile zone,'-', 
! ferries/cargo ships,'-', registred inside/outside EU
ix=ix+1 
IC_BA2=ix
Country( IC_BA2 ) = cc(  "BA2" ,'-',302 ,T,  30, 30, 1  , "Baltic EU cargo outs.12      " )
ix=ix+1 
IC_BA3=ix
Country( IC_BA3 ) = cc(  "BA3" ,'-',303 ,T,  30, 30, 1  , "Baltic ROW cargo outs. 12    " )
ix=ix+1 
IC_BA4=ix
Country( IC_BA4 ) = cc(  "BA4" ,'-',304 ,T,  30, 30, 1  , "Baltic EU cargo ins. 12      " )
ix=ix+1 
IC_BA5=ix
Country( IC_BA5 ) = cc(  "BA5" ,'-',305 ,T,  30, 30, 1  , "Baltic ROW cargo ins. 12     " )
ix=ix+1 
IC_BA6=ix
Country( IC_BA6 ) = cc(  "BA6" ,'-',306 ,T,  30, 30, 1  , "Baltic EU ferries outs.12    " )
ix=ix+1 
IC_BA7=ix
Country( IC_BA7 ) = cc(  "BA7" ,'-',307 ,T,  30, 30, 1  , "Baltic ROW ferries outs. 12  " )
ix=ix+1 
IC_BA8=ix
Country( IC_BA8 ) = cc(  "BA8" ,'-',308 ,T,  30, 30, 1  , "Baltic EU ferries ins. 12    " )
ix=ix+1 
IC_BA9=ix
Country( IC_BA9 ) = cc(  "BA9" ,'-',309 ,T,  30, 30, 1  , "Baltic ROW ferries ins. 12   " )

ix=ix+1 
IC_NS2=ix
Country( IC_NS2 ) = cc(  "NS2" ,'-',312 ,T,  31, 31, 1  , "N. Sea EU cargo outs.12      " )
ix=ix+1 
IC_NS3=ix
Country( IC_NS3 ) = cc(  "NS3" ,'-',313 ,T,  31, 31, 1  , "N. Sea ROW cargo outs. 12    " )
ix=ix+1 
IC_NS4=ix
Country( IC_NS4 ) = cc(  "NS4" ,'-',314 ,T,  31, 31, 1  , "N. Sea EU cargo ins. 12      " )
ix=ix+1 
IC_NS5=ix
Country( IC_NS5 ) = cc(  "NS5" ,'-',315 ,T,  31, 31, 1  , "N. Sea ROW cargo ins. 12     " )
ix=ix+1 
IC_NS6=ix
Country( IC_NS6 ) = cc(  "NS6" ,'-',316 ,T,  31, 31, 1  , "N. Sea EU ferries outs.12    " )
ix=ix+1 
IC_NS7=ix
Country( IC_NS7 ) = cc(  "NS7" ,'-',317 ,T,  31, 31, 1  , "N. Sea ROW ferries outs. 12  " )
ix=ix+1 
IC_NS8=ix
Country( IC_NS8 ) = cc(  "NS8" ,'-',318 ,T,  31, 31, 1  , "N. Sea EU ferries ins. 12    " )
ix=ix+1 
IC_NS9=ix
Country( IC_NS9 ) = cc(  "NS9" ,'-',319 ,T,  31, 31, 1  , "N. Sea ROW ferries ins. 12   " )

ix=ix+1 
IC_AT2=ix
Country( IC_AT2 ) = cc(  "AT2" ,'-',322 ,T,  32, 32, 1  , "Atlant EU cargo outs.12      " )
ix=ix+1 
IC_AT3=ix
Country( IC_AT3 ) = cc(  "AT3" ,'-',323 ,T,  32, 32, 1  , "Atlant ROW cargo outs. 12    " )
ix=ix+1 
IC_AT4=ix
Country( IC_AT4 ) = cc(  "AT4" ,'-',324 ,T,  32, 32, 1  , "Atlant EU cargo ins. 12      " )
ix=ix+1 
IC_AT5=ix
Country( IC_AT5 ) = cc(  "AT5" ,'-',325 ,T,  32, 32, 1  , "Atlant ROW cargo ins. 12     " )
ix=ix+1 
IC_AT6=ix
Country( IC_AT6 ) = cc(  "AT6" ,'-',326 ,T,  32, 32, 1  , "Atlant EU ferries outs.12    " )
ix=ix+1 
IC_AT7=ix
Country( IC_AT7 ) = cc(  "AT7" ,'-',327 ,T,  32, 32, 1  , "Atlant ROW ferries outs. 12  " )
ix=ix+1 
IC_AT8=ix
Country( IC_AT8 ) = cc(  "AT8" ,'-',328 ,T,  32, 32, 1  , "Atlant EU ferries ins. 12    " )
ix=ix+1 
IC_AT9=ix
Country( IC_AT9 ) = cc(  "AT9" ,'-',329 ,T,  32, 32, 1  , "Atlant ROW ferries ins. 12   " )

ix=ix+1 
IC_ME2=ix
Country( IC_ME2 ) = cc(  "ME2" ,'-',332 ,T,  33, 33, 1  , "Medite EU cargo outs.12      " )
ix=ix+1 
IC_ME3=ix
Country( IC_ME3 ) = cc(  "ME3" ,'-',333 ,T,  33, 33, 1  , "Medite ROW cargo outs. 12    " )
ix=ix+1 
IC_ME4=ix
Country( IC_ME4 ) = cc(  "ME4" ,'-',334 ,T,  33, 33, 1  , "Medite EU cargo ins. 12      " )
ix=ix+1 
IC_ME5=ix
Country( IC_ME5 ) = cc(  "ME5" ,'-',335 ,T,  33, 33, 1  , "Medite ROW cargo ins. 12     " )
ix=ix+1 
IC_ME6=ix
Country( IC_ME6 ) = cc(  "ME6" ,'-',336 ,T,  33, 33, 1  , "Medite EU ferries outs.12    " )
ix=ix+1 
IC_ME7=ix
Country( IC_ME7 ) = cc(  "ME7" ,'-',337 ,T,  33, 33, 1  , "Medite ROW ferries outs. 12  " )
ix=ix+1 
IC_ME8=ix
Country( IC_ME8 ) = cc(  "ME8" ,'-',338 ,T,  33, 33, 1  , "Medite EU ferries ins. 12    " )
ix=ix+1 
IC_ME9=ix
Country( IC_ME9 ) = cc(  "ME9" ,'-',339 ,T,  33, 33, 1  , "Medite ROW ferries ins. 12   " )

ix=ix+1 
IC_BL2=ix
Country( IC_BL2 ) = cc(  "BL2" ,'-',342 ,T,  34, 34, 2  , "B. Sea EU cargo outs.12      " )
ix=ix+1 
IC_BL3=ix
Country( IC_BL3 ) = cc(  "BL3" ,'-',343 ,T,  34, 34, 2  , "B. Sea ROW cargo outs. 12    " )
ix=ix+1 
IC_BL4=ix
Country( IC_BL4 ) = cc(  "BL4" ,'-',344 ,T,  34, 34, 2  , "B. Sea EU cargo ins. 12      " )
ix=ix+1 
IC_BL5=ix
Country( IC_BL5 ) = cc(  "BL5" ,'-',345 ,T,  34, 34, 2  , "B. Sea ROW cargo ins. 12     " )
ix=ix+1 
IC_BL6=ix
Country( IC_BL6 ) = cc(  "BL6" ,'-',346 ,T,  34, 34, 2  , "B. Sea EU ferries outs.12    " )
ix=ix+1 
IC_BL7=ix
Country( IC_BL7 ) = cc(  "BL7" ,'-',347 ,T,  34, 34, 2  , "B. Sea ROW ferries outs. 12  " )
ix=ix+1 
IC_BL8=ix
Country( IC_BL8 ) = cc(  "BL8" ,'-',348 ,T,  34, 34, 2  , "B. Sea EU ferries ins. 12    " )
ix=ix+1 
IC_BL9=ix
Country( IC_BL9 ) = cc(  "BL9" ,'-',349 ,T,  34, 34, 2  , "B. Sea ROW ferries ins. 12   " )



!!  HTAP2 regions
ix=ix+1 
IC_HT1000 = ix
Country(IC_HT1000 ) = cc(  "HT1000" ,'-',1000 ,T, 30, 30, -100  , "Unefined" )
ix=ix+1 
IC_HTNATL = ix
Country(IC_HTNATL ) = cc(  "N_ATL" ,'-',1002 ,T, 32, 32, -100  , "Int. ships, N. Atl." )
ix=ix+1 
IC_HTUSCA = ix
Country(IC_HTUSCA ) = cc(  "USCA" ,'-',1003 ,T, 65,  65, -100  , "USA and Canada" )
ix=ix+1 
IC_HTEUTU = ix
Country(IC_HTEUTU ) = cc(  "EU_TU" ,'-',1004 ,T, 64, 64, 1  , "EU and Turkey" )
ix=ix+1 
IC_HTSASI = ix
Country(IC_HTSASI ) = cc(  "S_ASIA_IP" ,'-',1005 ,T, 212,212, -100  , "S. Asia India Pak." )
ix=ix+1 
IC_HTEASI = ix
Country(IC_HTEASI ) = cc(  "CHCORJAP" ,'-',1006 ,T, 210,  210, -100  , "China Korea Japan" )
ix=ix+1 
IC_HTSEAS = ix
Country(IC_HTSEAS ) = cc(  "TAINDOMA" ,'-',1007 ,T, 210,210, -100  , "Thail. Indon. Malay" )
ix=ix+1 
IC_HTAUST = ix
Country(IC_HTAUST ) = cc(  "AUSTR" ,'-',1008 ,T, 203,203, -100  , "Austr N. Zeal. ++" )
ix=ix+1 
IC_HTNAFR = ix
Country(IC_HTNAFR ) = cc(  "NAFRI" ,'-',1009 ,T, 64, 64,1  , "N. Africa" )
ix=ix+1 
IC_HTRAFR = ix
Country(IC_HTRAFR ) = cc(  "RAFRI" ,'-',1010 ,T, 64, 64, 1  , "Rest Africa" )
ix=ix+1 
IC_HTMIDE = ix
Country(IC_HTMIDE ) = cc(  "MIDEAST" ,'-',1011 ,T, 25, 25, 2  , "Middle East" )
ix=ix+1 
IC_HTCEAM = ix
Country(IC_HTCEAM ) = cc(  "HTCEAM" ,'-',1012 ,T, 243, 243, -100  , "C. Am. Carib" )
ix=ix+1 
IC_HTSOAM = ix
Country(IC_HTSOAM ) = cc(  "HTSOAM" ,'-',1013 ,T, 227, 227, -100  , "S. America" )
ix=ix+1 
IC_HTRUBU = ix
Country(IC_HTRUBU ) = cc(  "HTRUBU" ,'-',1014 ,T, 39, 39, -100  , "Russ Bel. Ukr" )
ix=ix+1 
IC_HTCASI = ix
Country(IC_HTCASI ) = cc(  "HTCASI" ,'-',1015 ,T, 92, 92, -100  , "C. Asia" )
ix=ix+1 
IC_HTPOLA = ix
Country(IC_HTPOLA ) = cc(  "HTPOLA" ,'-',1016 ,T, 30, 30, -100 , "Pol. N 66" )
ix=ix+1
IC_HTSPOL = ix
Country(IC_HTSPOL ) = cc(  "HTSPOL" ,'-',1017 ,T, 30, 30, -100  , "S. Pol .S 60" )
ix=ix+1 
IC_HT1018 = ix
Country(IC_HT1018 ) = cc(  "HT1018",'-' ,1018 ,T,  30, 30,-100  , "undefined" )
ix=ix+1 
IC_HT1019 = ix
Country(IC_HT1019 ) = cc(  "HT1019",'-' ,1019 ,T,  30, 30, -100,  "undefined" )
ix=ix+1 
IC_HT1020 = ix
Country(IC_HT1020 ) = cc(  "HT1020",'-' ,1020 ,T,  30, 30, -100,  "undefined" )

!UNEP SR and new GAINS regions'
ix=ix+1
IC_BANG_DHAK = ix
Country(IC_BANG_DHAK) = cc( "BANG_DHAK",'-',  252 , F,   252 ,252 ,-100, "Bangladesh:Dhaka")
ix=ix+1
IC_BANG_REST = ix
Country(IC_BANG_REST) = cc( "BANG_REST",'-',  253 , F,   253 ,253 ,-100, "Rest_of_Bangladesh")
ix=ix+1
IC_CHIN_ANHU = ix
Country(IC_CHIN_ANHU) = cc( "CHIN_ANHU",'-',  254 , F,   254 ,254 ,-100, "China:Anhui")
ix=ix+1
IC_CHIN_BEIJ = ix
Country(IC_CHIN_BEIJ) = cc( "CHIN_BEIJ",'-',  255 , F,   255 ,255 ,-100, "China:Beijing")
ix=ix+1
IC_CHIN_CHON = ix
Country(IC_CHIN_CHON) = cc( "CHIN_CHON",'-',  256 , F,   256 ,256 ,-100, "China:Chongqing")
ix=ix+1
IC_CHIN_FUJI = ix
Country(IC_CHIN_FUJI) = cc( "CHIN_FUJI",'-',  257 , F,   257 ,257 ,-100, "China:Fujian")
ix=ix+1
IC_CHIN_GANS = ix
Country(IC_CHIN_GANS) = cc( "CHIN_GANS",'-',  258 , F,   258 ,258 ,-100, "China:Gansu")
ix=ix+1
IC_CHIN_GUAD = ix
Country(IC_CHIN_GUAD) = cc( "CHIN_GUAD",'-',  259 , F,   259 ,259 ,-100, "China:Guangdong")
ix=ix+1
IC_CHIN_GUAX = ix
Country(IC_CHIN_GUAX) = cc( "CHIN_GUAX",'-',  260 , F,   260 ,260 ,-100, "China:Guangxi")
ix=ix+1
IC_CHIN_GUIZ = ix
Country(IC_CHIN_GUIZ) = cc( "CHIN_GUIZ",'-',  261 , F,   261 ,261 ,-100, "China:Guizhou")
ix=ix+1
IC_CHIN_HAIN = ix
Country(IC_CHIN_HAIN) = cc( "CHIN_HAIN",'-',  262 , F,   262 ,262 ,-100, "China:Hainan")
ix=ix+1
IC_CHIN_HEBE = ix
Country(IC_CHIN_HEBE) = cc( "CHIN_HEBE",'-',  263 , F,   263 ,263 ,-100, "China:Hebei")
ix=ix+1
IC_CHIN_HEIL = ix
Country(IC_CHIN_HEIL) = cc( "CHIN_HEIL",'-',  264 , F,   264 ,264 ,-100, "China:Heilongjiang")
ix=ix+1
IC_CHIN_HENA = ix
Country(IC_CHIN_HENA) = cc( "CHIN_HENA",'-',  265 , F,   265 ,265 ,-100, "China:Henan")
ix=ix+1
IC_CHIN_HONG = ix
Country(IC_CHIN_HONG) = cc( "CHIN_HONG",'-',  266 , F,   266 ,266 ,-100, "China:Hong_Kong_&_Macau")
ix=ix+1
IC_CHIN_HUBE = ix
Country(IC_CHIN_HUBE) = cc( "CHIN_HUBE",'-',  267 , F,   267 ,267 ,-100, "China:Hubei")
ix=ix+1
IC_CHIN_HUNA = ix
Country(IC_CHIN_HUNA) = cc( "CHIN_HUNA",'-',  268 , F,   268 ,268 ,-100, "China:Hunan")
ix=ix+1
IC_CHIN_JILI = ix
Country(IC_CHIN_JILI) = cc( "CHIN_JILI",'-',  269 , F,   269 ,269 ,-100, "China:Jilin")
ix=ix+1
IC_CHIN_JINU = ix
Country(IC_CHIN_JINU) = cc( "CHIN_JINU",'-',  270 , F,   270 ,270 ,-100, "China:Jiangsu")
ix=ix+1
IC_CHIN_JINX = ix
Country(IC_CHIN_JINX) = cc( "CHIN_JINX",'-',  271 , F,   271 ,271 ,-100, "China:Jiangxi")
ix=ix+1
IC_CHIN_LIAO = ix
Country(IC_CHIN_LIAO) = cc( "CHIN_LIAO",'-',  272 , F,   272 ,272 ,-100, "China:Liaoning")
ix=ix+1
IC_CHIN_NEMO = ix
Country(IC_CHIN_NEMO) = cc( "CHIN_NEMO",'-',  273 , F,   273 ,273 ,-100, "China:Inner_Mongolia")
ix=ix+1
IC_CHIN_NINX = ix
Country(IC_CHIN_NINX) = cc( "CHIN_NINX",'-',  274 , F,   274 ,274 ,-100, "China:Ningxia")
ix=ix+1
IC_CHIN_QING = ix
Country(IC_CHIN_QING) = cc( "CHIN_QING",'-',  275 , F,   275 ,275 ,-100, "China:Qinghai")
ix=ix+1
IC_CHIN_SHAA = ix
Country(IC_CHIN_SHAA) = cc( "CHIN_SHAA",'-',  276 , F,   276 ,276 ,-100, "China:Shaanxi")
ix=ix+1
IC_CHIN_SHAN = ix
Country(IC_CHIN_SHAN) = cc( "CHIN_SHAN",'-',  277 , F,   277 ,277 ,-100, "China:Shanghai")
ix=ix+1
IC_CHIN_SHND = ix
Country(IC_CHIN_SHND) = cc( "CHIN_SHND",'-',  278 , F,   278 ,278 ,-100, "China:Shandong")
ix=ix+1
IC_CHIN_SHNX = ix
Country(IC_CHIN_SHNX) = cc( "CHIN_SHNX",'-',  279 , F,   279 ,279 ,-100, "China:Shanxi")
ix=ix+1
IC_CHIN_SICH = ix
Country(IC_CHIN_SICH) = cc( "CHIN_SICH",'-',  280 , F,   280 ,280 ,-100, "China:Sichuan")
ix=ix+1
IC_CHIN_TIAN = ix
Country(IC_CHIN_TIAN) = cc( "CHIN_TIAN",'-',  281 , F,   281 ,281 ,-100, "China:Tianjin")
ix=ix+1
IC_CHIN_TIBE = ix
Country(IC_CHIN_TIBE) = cc( "CHIN_TIBE",'-',  282 , F,   282 ,282 ,-100, "China:Tibet_(Xizang)")
ix=ix+1
IC_CHIN_XING = ix
Country(IC_CHIN_XING) = cc( "CHIN_XING",'-',  283 , F,   283 ,283 ,-100, "China:Xinjiang")
ix=ix+1
IC_CHIN_YUNN = ix
Country(IC_CHIN_YUNN) = cc( "CHIN_YUNN",'-',  284 , F,   284 ,284 ,-100, "China:Yunnan")
ix=ix+1
IC_CHIN_ZHEJ = ix
Country(IC_CHIN_ZHEJ) = cc( "CHIN_ZHEJ",'-',  285 , F,   285 ,285 ,-100, "China:Zhejiang")
ix=ix+1
IC_INDI_ANPR = ix
Country(IC_INDI_ANPR) = cc( "INDI_ANPR",'-',  286 , F,   286 ,286 ,-100, "India:Andhra_Pradesh")
ix=ix+1
IC_INDI_ASSA = ix
Country(IC_INDI_ASSA) = cc( "INDI_ASSA",'-',  287 , F,   287 ,287 ,-100, "India:Assam")
ix=ix+1
IC_INDI_BENG = ix
Country(IC_INDI_BENG) = cc( "INDI_BENG",'-',  288 , F,   288 ,288 ,-100, "India:West_Bengal")
ix=ix+1
IC_INDI_BIHA = ix
Country(IC_INDI_BIHA) = cc( "INDI_BIHA",'-',  289 , F,   289 ,289 ,-100, "India:Bihar")
ix=ix+1
IC_INDI_CHHA = ix
Country(IC_INDI_CHHA) = cc( "INDI_CHHA",'-',  290 , F,   290 ,290 ,-100, "India:Chhattisgarh")
ix=ix+1
IC_INDI_DELH = ix
Country(IC_INDI_DELH) = cc( "INDI_DELH",'-',  291 , F,   291 ,291 ,-100, "India:Delhi")
ix=ix+1
IC_INDI_EHIM = ix
Country(IC_INDI_EHIM) = cc( "INDI_EHIM",'-',  292 , F,   292 ,292 ,-100, "India:North_East_(excl._Assam)")
ix=ix+1
IC_INDI_GOA = ix
Country(IC_INDI_GOA) = cc( "INDI_GOA",'-',  293 , F,   293 ,293 ,-100, "India:Goa")
ix=ix+1
IC_INDI_GUJA = ix
Country(IC_INDI_GUJA) = cc( "INDI_GUJA",'-',  294 , F,   294 ,294 ,-100, "India:Gujarat")
ix=ix+1
IC_INDI_HARY = ix
Country(IC_INDI_HARY) = cc( "INDI_HARY",'-',  295 , F,   295 ,295 ,-100, "India:Haryana")
ix=ix+1
IC_INDI_HIPR = ix
Country(IC_INDI_HIPR) = cc( "INDI_HIPR",'-',  296 , F,   296 ,296 ,-100, "India:Himachal_Pradesh")
ix=ix+1
IC_INDI_JHAR = ix
Country(IC_INDI_JHAR) = cc( "INDI_JHAR",'-',  297 , F,   297 ,297 ,-100, "India:Jharkhand")
ix=ix+1
IC_INDI_KARN = ix
Country(IC_INDI_KARN) = cc( "INDI_KARN",'-',  298 , F,   298 ,298 ,-100, "India:Karnataka")
ix=ix+1
IC_INDI_KERA = ix
Country(IC_INDI_KERA) = cc( "INDI_KERA",'-',  299 , F,   299 ,299 ,-100, "India:Kerala")
ix=ix+1
IC_INDI_MAHA = ix
Country(IC_INDI_MAHA) = cc( "INDI_MAHA",'-',  300 , F,   300 ,300 ,-100, "India:Maharashtra-Dadra-Nagar-Haveli-Daman-Diu")
ix=ix+1
IC_INDI_MAPR = ix
Country(IC_INDI_MAPR) = cc( "INDI_MAPR",'-',  351 , F,   351 ,351 ,-100, "India:Madhya_Pradesh")
ix=ix+1
IC_INDI_ORIS = ix
Country(IC_INDI_ORIS) = cc( "INDI_ORIS",'-',  352 , F,   352 ,352 ,-100, "India:Orissa")
ix=ix+1
IC_INDI_PUNJ = ix
Country(IC_INDI_PUNJ) = cc( "INDI_PUNJ",'-',  353 , F,   353 ,353 ,-100, "India:Punjab")
ix=ix+1
IC_INDI_RAJA = ix
Country(IC_INDI_RAJA) = cc( "INDI_RAJA",'-',  354 , F,   354 ,354 ,-100, "India:Rajasthan")
ix=ix+1
IC_INDI_TAMI = ix
Country(IC_INDI_TAMI) = cc( "INDI_TAMI",'-',  355 , F,   355 ,355 ,-100, "India:Tamil_Nadu")
ix=ix+1
IC_INDI_UTAN = ix
Country(IC_INDI_UTAN) = cc( "INDI_UTAN",'-',  356 , F,   356 , 356 ,-100, "India:Uttaranchal")
ix=ix+1
IC_INDI_UTPR = ix
Country(IC_INDI_UTPR) = cc( "INDI_UTPR",'-',  357 , F,   357 ,357 ,-100, "India:Uttar_Pradesh")
ix=ix+1
IC_INDI_WHIM = ix
Country(IC_INDI_WHIM) = cc( "INDI_WHIM",'-',  358 , F,   358 ,358 ,-100, "India:Jammu_and_Kashmir")
ix=ix+1
IC_INDO_JAKA = ix
Country(IC_INDO_JAKA) = cc( "INDO_JAKA",'-',  359 , F,   359 ,359 ,-100, "Indonesia:Jakarta")
ix=ix+1
IC_INDO_JAVA = ix
Country(IC_INDO_JAVA) = cc( "INDO_JAVA",'-',  360 , F,   360 ,360 ,-100, "Indonesia:Java")
ix=ix+1
IC_INDO_REST = ix
Country(IC_INDO_REST) = cc( "INDO_REST",'-',  361 , F,   361 ,361 ,-100, "Indonesia:Rest_of_Indonesia")
ix=ix+1
IC_INDO_SUMA = ix
Country(IC_INDO_SUMA) = cc( "INDO_SUMA",'-',  362 , F,   362 ,362 ,-100, "Indonesia:Sumatra")
ix=ix+1
IC_JAPA_CHSH = ix
Country(IC_JAPA_CHSH) = cc( "JAPA_CHSH",'-',  363 , F,   363 ,363 ,-100, "Japan:Chugoku-Shikoku")
ix=ix+1
IC_JAPA_CHUB = ix
Country(IC_JAPA_CHUB) = cc( "JAPA_CHUB",'-',  364 , F,   364 ,364 ,-100, "Japan:Chubu")
ix=ix+1
IC_JAPA_HOTO = ix
Country(IC_JAPA_HOTO) = cc( "JAPA_HOTO",'-',  365 , F,   365 ,365 ,-100, "Japan:Hokkaido-Tohoku")
ix=ix+1
IC_JAPA_KANT = ix
Country(IC_JAPA_KANT) = cc( "JAPA_KANT",'-',  366 , F,   366 ,366 ,-100, "Japan:Kanto")
ix=ix+1
IC_JAPA_KINK = ix
Country(IC_JAPA_KINK) = cc( "JAPA_KINK",'-',  367 , F,   367 ,367 ,-100, "Japan:Kinki")
ix=ix+1
IC_JAPA_KYOK = ix
Country(IC_JAPA_KYOK) = cc( "JAPA_KYOK",'-',  368 , F,   368 ,368 ,-100, "Japan:Kyushu-Okinawa")
ix=ix+1
IC_KORS_NORT = ix
Country(IC_KORS_NORT) = cc( "KORS_NORT",'-',  369 , F,   369 , 369 ,-100, "South_Korea:North")
ix=ix+1
IC_KORS_PUSA = ix
Country(IC_KORS_PUSA) = cc( "KORS_PUSA",'-',  370 , F,   370 ,370 ,-100, "South_Korea:Pusan")
ix=ix+1
IC_KORS_SEOI = ix
Country(IC_KORS_SEOI) = cc( "KORS_SEOI",'-',  371 , F,   371 ,371 ,-100, "South_Korea:Seoul-Inchon")
ix=ix+1
IC_KORS_SOUT = ix
Country(IC_KORS_SOUT) = cc( "KORS_SOUT",'-',  372 , F,   372 ,372 ,-100, "South_Korea:South")
ix=ix+1
IC_MALA_KUAL = ix
Country(IC_MALA_KUAL) = cc( "MALA_KUAL",'-',  374 , F,   374 ,374 ,-100, "Malaysia:Kuala_Lumpur")
ix=ix+1
IC_MALA_PENM = ix
Country(IC_MALA_PENM) = cc( "MALA_PENM",'-',  375 , F,   375 , 375 ,-100, "Malaysia:Peninsular_Malaysia")
ix=ix+1
IC_MALA_SASA = ix
Country(IC_MALA_SASA) = cc( "MALA_SASA",'-',  376 , F,   376 ,376 ,-100, "Malaysia:Sarawak-Sabah")
ix=ix+1
IC_PAKI_KARA = ix
Country(IC_PAKI_KARA) = cc( "PAKI_KARA",'-',  377 , F,   377 ,377 ,-100, "Pakistan:Karachi")
ix=ix+1
IC_PAKI_NMWP = ix
Country(IC_PAKI_NMWP) = cc( "PAKI_NMWP",'-',  378 , F,   378 ,378 ,-100, "Pakistan:NW_Frontier_Provinces-Baluchistan")
ix=ix+1
IC_PAKI_PUNJ = ix
Country(IC_PAKI_PUNJ) = cc( "PAKI_PUNJ",'-',  379 , F,   379 ,379 ,-100, "Pakistan:Punjab")
ix=ix+1
IC_PAKI_SIND = ix
Country(IC_PAKI_SIND) = cc( "PAKI_SIND",'-',  380 , F,   380 ,380 ,-100, "Pakistan:Sind")
ix=ix+1
IC_PHIL_BVMI = ix
Country(IC_PHIL_BVMI) = cc( "PHIL_BVMI",'-',  381 , F,   381 ,381 ,-100, "Philipinnes:Bicol-Visayas-Mindanao")
ix=ix+1
IC_PHIL_LUZO = ix
Country(IC_PHIL_LUZO) = cc( "PHIL_LUZO",'-',  382 , F,   382 ,382 ,-100, "Philipinnes:Luzon")
ix=ix+1
IC_PHIL_MANI = ix
Country(IC_PHIL_MANI) = cc( "PHIL_MANI",'-',  383 , F,   383 ,383 ,-100, "Philipinnes:Metropolitan_Manila")
ix=ix+1
IC_RUSS_ASIA = ix
Country(IC_RUSS_ASIA) = cc( "RUSS_ASIA",'-',  384 , F,   384 ,384 ,-100, "Russia:Asian_part")
ix=ix+1
IC_RUSS_EURO = ix
Country(IC_RUSS_EURO) = cc( "RUSS_EURO",'RUSS',  385 , F,   385 ,385 ,-100, "Russia:European_part")
ix=ix+1
IC_THAI_BANG = ix
Country(IC_THAI_BANG) = cc( "THAI_BANG",'-',  386 , F,   386 ,386 ,-100, "Thailand:Bangkok_Metropolitan_Region")
ix=ix+1
IC_THAI_CVAL = ix
Country(IC_THAI_CVAL) = cc( "THAI_CVAL",'-',  387 , F,   387 ,387 ,-100, "Thailand:Central_Valley")
ix=ix+1
IC_THAI_NEPL = ix
Country(IC_THAI_NEPL) = cc( "THAI_NEPL",'-',  388 , F,   388 ,388 ,-100, "Thailand:NE_Plateau")
ix=ix+1
IC_THAI_NHIG = ix
Country(IC_THAI_NHIG) = cc( "THAI_NHIG",'-',  389 , F,   389 ,389 ,-100, "Thailand:N_Highlands")
ix=ix+1
IC_THAI_SPEN = ix
Country(IC_THAI_SPEN) = cc( "THAI_SPEN",'-',  390 , F,   390 ,390 ,-100, "Thailand:S_Peninsula")
ix=ix+1
IC_IRAN_REST = ix
Country(IC_IRAN_REST) = cc( "IRAN_REST",'-',  391 , F,   391 ,391 ,-100, "Iran:Rest_of_Iran")
ix=ix+1
IC_IRAN_TEHR = ix
Country(IC_IRAN_TEHR) = cc( "IRAN_TEHR",'-',  392 , F,   392 ,392 ,-100, "Iran:Teheran")
ix=ix+1
IC_VIET_BNIH = ix
Country(IC_VIET_BNIH) = cc( "VIET_BNIH", '-', 470 , F,   239 ,239 ,-100, "Vietnam:Bac_Ninh")
ix=ix+1
IC_VIET_HYEN = ix
Country(IC_VIET_HYEN) = cc( "VIET_HYEN",'-',  471 , F,   239 ,239 ,-100, "Vietnam:Hung_Yen")
ix=ix+1
IC_VIET_OCAR = ix
Country(IC_VIET_OCAR) = cc( "VIET_OCAR",'-',  472 , F,   239 ,239 , -100, "Vietnam:Other_Ha_Noi_Capital_Region")
ix=ix+1
IC_VIET_ONOR = ix
Country(IC_VIET_ONOR) = cc( "VIET_ONOR",'-',  473 , F,   239 ,239 ,-100, "Vietnam:Other_Northern_Vietnam")
ix=ix+1
IC_USAM_ALAS = ix
Country(IC_USAM_ALAS) = cc( "USAM_ALAS",'-',  474 , F,    65 ,65, -100, "Alaska")
ix=ix+1
IC_USAM_REST = ix
Country(IC_USAM_REST) = cc( "USAM_REST", '-', 475 , F,    65 ,65, -100, "United_States_(without_Alaska)")
ix=ix+1
IC_EGYP_CAIR = ix
Country(IC_EGYP_CAIR) = cc( "EGYP_CAIR",'-',  476 , F,    63 , 63, -100, "Egypt:Cairo_area")
ix=ix+1
IC_EGYP_REST = ix
Country(IC_EGYP_REST) = cc( "EGYP_REST",'-',  477 , F,    63 , 63, -100, "Egypt_(remaining)")
ix=ix+1
IC_NAFR_WHOL = ix
Country(IC_NAFR_WHOL) = cc( "NAFR_WHOL",'-',  478 , F,   224 , 224,-100, "Algeria_Libya_Morocco_Tunisia_Western_Sahara")
ix=ix+1
IC_NIGE_LAGO = ix
Country(IC_NIGE_LAGO) = cc( "NIGE_LAGO",'-',  479 , F,   226 , 226, -100, "Lagos_State_(Nigeria)")
ix=ix+1
IC_NIGE_OLAG = ix
Country(IC_NIGE_OLAG) = cc( "NIGE_OLAG",'-',  480 , F,   226 , 226,-100, "Ogun_Osun_Oyo_states_(Nigeria)")
ix=ix+1
IC_NIGE_REST = ix
Country(IC_NIGE_REST) = cc( "NIGE_REST",'-',  481 , F,   226 , 226,-100, "Nigeria_excl_Lagos_Ogun_Osun_Oyo")
ix=ix+1
IC_SAFR_JOBG = ix
Country(IC_SAFR_JOBG) = cc( "SAFR_JOBG",'-',  482 , F,   231 , 231, -100, "South_Africa:Greater_Johannesburg")
ix=ix+1
IC_SAFR_REST = ix
Country(IC_SAFR_REST) = cc( "SAFR_REST",'-',  483 , F,   231 , 231, -100, "South_Africa_(remaining)")
ix=ix+1
IC_RSAF_WHOL = ix
Country(IC_RSAF_WHOL) = cc( "RSAF_WHOL",'-',  484 , F,   226 , 226, -100, "Rest_of_Southern_Africa")
ix=ix+1
IC_EAFR_WHOL = ix
Country(IC_EAFR_WHOL) = cc( "EAFR_WHOL",'-',  485 , F,   226 , 226, -100, "Eastern_Africa")
ix=ix+1
IC_WAFR_WHOL = ix
Country(IC_WAFR_WHOL) = cc( "WAFR_WHOL",'-',  486 , F,   226 , 226, -100, "Western_Africa")
ix=ix+1
IC_KENY_WHOL = ix
Country(IC_KENY_WHOL) = cc( "KENY_WHOL",'-',  487 , F,   226 , 226, -100, "Kenya")
ix=ix+1
IC_TANZ_WHOL = ix
Country(IC_TANZ_WHOL) = cc( "TANZ_WHOL",'-',  488 , F,   226 , 226, -100, "Tanzania")
ix=ix+1
IC_NIGE_WHOL = ix
Country(IC_NIGE_WHOL) = cc( "NIGE_WHOL",'-',  489 , F,   226 , 226, -100, "Nigeria")


!City codes
ix=ix+1
IC_BERLIN = ix
Country(IC_BERLIN) = cc( "BERLIN",  '-', 401 , F,   60, 401 , 1, "DE:BERLIN")
ix=ix+1
IC_BRUSSEL = ix
Country(IC_BRUSSEL) = cc( "BRUSSEL", '-',  402 , F,   3, 402 , 1, "BE:BRUSSEL")
ix=ix+1
IC_COPENHAGEN = ix
Country(IC_COPENHAGEN) = cc( "COPENHAGEN", '-',  403 , F,   6, 403 , 1, "DK:COPENHAGEN")
ix=ix+1
IC_MADRID = ix
Country(IC_MADRID) = cc( "MADRID", '-',  404 , F,   22, 404 , 1, "ES:MADRID")
ix=ix+1
IC_HELSINKI = ix
Country(IC_HELSINKI) = cc( "HELSINKI", '-',  405 , F,   7, 405 , 2, "FI:HELSINKI")
ix=ix+1
IC_PARIS = ix
Country(IC_PARIS) = cc( "PARIS", '-',  406 , F,   8, 406 , 1, "FR:PARIS")
ix=ix+1
IC_ATHENS = ix
Country(IC_ATHENS) = cc( "ATHENS", '-',  407 , F,   11, 407 , 2, "GR:ATHENS")
ix=ix+1
IC_BUDAPEST = ix
Country(IC_BUDAPEST) = cc( "BUDAPEST", '-',  408 , F,   12, 408 , 1, "HU:BUDAPEST")
ix=ix+1
IC_DUBLIN = ix
Country(IC_DUBLIN) = cc( "DUBLIN", '-',  409 , F,   14, 409 , 0, "IE:DUBLIN")
ix=ix+1
IC_MILAN = ix
Country(IC_MILAN) = cc( "MILAN", '-',  410 , F,   15, 410 , 1, "IT:MILAN")
ix=ix+1
IC_OSLO = ix
Country(IC_OSLO) = cc( "OSLO", '-',  411 , F,   18, 411 , 1, "NO:OSLO")
ix=ix+1
IC_AMSTERDAM = ix
Country(IC_AMSTERDAM) = cc( "AMSTERDAM", '-',  412 , F,   17, 412 , 1, "NL:AMSTERDAM")
ix=ix+1
IC_WARSAW = ix
Country(IC_WARSAW) = cc( "WARSAW", '-',  413 , F,   19, 413 , 1, "PL:WARSAW")
ix=ix+1
IC_LISBON = ix
Country(IC_LISBON) = cc( "LISBON", '-',  414 , F,   20, 414 , 0, "PT:LISBON")
ix=ix+1
IC_LONDON = ix
Country(IC_LONDON) = cc( "LONDON", '-',  415 , F,   27, 415 , 0, "GB:LONDON")
ix=ix+1
IC_STOCKHOLM = ix
Country(IC_STOCKHOLM) = cc( "STOCKHOLM", '-',  416 , F,   23, 416 , 1, "SE:STOCKHOLM")
ix=ix+1
IC_GENEVA = ix
Country(IC_GENEVA) = cc( "GENEVA", '-',  417 , F,   24, 417 , 1, "CH:GENEVA")

! CAMS-TNO sea regions
ix=ix+1 
IC_GRS=ix
Country( IC_GRS) = cc(  "GRS" , '-', 501 ,T, 32, 32,  -100  , "Greenland Sea")
ix=ix+1 
IC_BAR=ix
Country( IC_BAR) = cc(  "BAR" , '-', 502 ,T, 32, 32,  -100  , "Barents Sea")
ix=ix+1 
IC_NWS=ix
Country( IC_NWS) = cc(  "NWS" , '-', 503 ,T, 32, 32,  -100  , "Norwegian Sea")
ix=ix+1 
IC_ENC=ix
Country( IC_ENC) = cc(  "ENC" , '-', 504 ,T, 32, 32,  -100  , "English Channel")
ix=ix+1 
IC_IRC=ix
Country( IC_IRC) = cc(  "IRC" , '-', 505 ,T, 32, 32,  -100  , "Irish Sea")
ix=ix+1 
IC_KAR=ix
Country( IC_KAR) = cc(  "KAR" , '-', 506 ,T, 32, 32,  -100  , "Kara Sea")
ix=ix+1 
IC_PSG=ix
Country( IC_PSG) = cc(  "PSG" , '-', 507 ,T, 507, 507,  -100  , "Persian Gulf")



NLAND=ix !actual number of countries defined

  end subroutine init_Country

  subroutine self_test()
    integer :: ic
    print *, "COUNTRY TEST ==================================="
    call init_Country()
    print *, "COUNTRY TEST NLAND = ", NLAND
    do  ic = 1, NLAND
      print '(a,i3,2x,a,i5,2x,a)', "IC ", ic, Country(ic)%code, &
           Country(ic)%icode, Country(ic)%gains
    end do
  end subroutine self_test

end module Country_mod
!TSTEMX program testr
!TSTEMX use Country_mod
!TSTEMX  call init_Country()     ! sets country details
!TSTEMX  call self_test()     ! just to test numbering
!TSTEMX end program testr
