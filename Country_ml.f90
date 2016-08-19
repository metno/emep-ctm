! <Country_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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

 !+ Sets country index numbers (IC_xx), code, time-zones, and names
 !
 ! Language : F
 ! History  :
 !   4th version: 06/03/2007 - jej/ds, sea areas split into flags etc,
 !   3rd version: 13/11/2003 - jej/ds, new countries, corrected 29
 !   2nd version: 12/12/2001 - ds,  added is_sea to country type
 !   1st version: 24/1/2001 - ds, using timezones from Vigdis.
 !
 ! Notes: 
 ! Unfortunately doing the setting of the Country data needs the
 ! subroutine Set_Countries, as doing this with a neat "parameter"
 ! statement didn't work - there were too many continuation lines :-(
 ! And I would have used the word countries instead of cc, but that made
 ! the lines too long in Set-Countries....

  ! Regions Atlantic (70) and Russia (71) outside EMEP defined separately, 
  ! as total emissions for Russia and the Atlantic often are reported for 
  ! EMEP domain only and then gridded according to this total. 

  ! timefac_index under cc as defined below assigns timefactors to 
  ! country/region/emission_type. As an example defining Bavaria as 
  ! a separate region with timefactors as in Germany.

  implicit none

  public :: Country_Init     ! sets country details

  integer, parameter, public :: NLAND = 349
  logical, parameter, private :: T = .true.   ! shorthand
  logical, parameter, private :: F = .false.  ! shorthand

  !/ to be set in Country_Init:

  type, public :: cc
     character(len=3)  :: code          ! up to 3 letter land code
     integer           :: index         ! index number (corresponds to 
                                        ! numbering in emission files read inn
     logical           :: is_sea        ! 1 for sea area, 0 otherwise
     integer           :: timefac_index ! see explanation above
     integer           :: timezone      ! timezone, deviation from UTC time
     character(len=30) :: name          ! name of country/region/emission_type
  end type cc

  type(cc), public, save, dimension(NLAND) :: Country

  integer, parameter, public ::  IC_AL =   1   ! Albania                       
  integer, parameter, public ::  IC_AT =   2   ! Austria                       
  integer, parameter, public ::  IC_BE =   3   ! Belgium                       
  integer, parameter, public ::  IC_BG =   4   ! Bulgaria                      
  integer, parameter, public ::  IC_CS =   5   ! Former Yugoslavia             
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
  integer, parameter, public :: IC_REM =  29   ! Remaining                     
  integer, parameter, public :: IC_BAS =  30   ! The Baltic Sea                
  integer, parameter, public :: IC_NOS =  31   ! The North Sea                 
  integer, parameter, public :: IC_ATL =  32   ! Remaining NE Atlantic Ocean   
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
  integer, parameter, public ::  IC_YU =  51   ! Yugoslavia                    
  integer, parameter, public ::  IC_MK =  52   ! Macedonia,                    
  integer, parameter, public ::  IC_KZ =  53   ! Kazakstan                     
  integer, parameter, public ::  IC_GE =  54   ! Georgia                       
  integer, parameter, public ::  IC_CY =  55   ! Cyprus                        
  integer, parameter, public ::  IC_AM =  56   ! Armenia                       
  integer, parameter, public ::  IC_MT =  57   ! Malta                         
  integer, parameter, public :: IC_ASI =  58   ! Other Asian             
  integer, parameter, public ::  IC_LI =  59   ! Lihtenstein                   
  integer, parameter, public ::  IC_DE =  60   ! Germany                       
  integer, parameter, public ::  IC_RU =  61   ! Russian                       
  integer, parameter, public ::  IC_MC =  62   ! Monaco                        
  integer, parameter, public :: IC_NOA =  63   ! North Africa                  
  integer, parameter, public ::  IC_EU =  64   ! European
  integer, parameter, public ::  IC_US =  65   ! USA
  integer, parameter, public ::  IC_CA =  66   ! Canada
  integer, parameter, public ::  IC_DUMMY1 =  67 ! Not-defined
  integer, parameter, public :: IC_KG  =  68 ! Kyrgyzstan(outside dommain)
  integer, parameter, public :: IC_AZ  =  69   ! Azerbaijan                 
  integer, parameter, public :: IC_ATX  =  70   ! ATL outside emep
  integer, parameter, public :: IC_RUX  =  71   ! RU outside emep
  integer, parameter, public :: IC_RS  =  72   ! Serbia
  integer, parameter, public :: IC_ME  =  73   ! Montenegro


  ! extra subdivisions:
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

  ! Remaining NE Atlantic  (32)
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





  contains
  !
  subroutine Country_Init()

  ! Set the country details. Note that time-zones for some areas are either
  ! difficult (Russia should be 3 to 12) or not relevant (sea areas,
  ! volcanoes). This needs to be thought about in using these figures.

  integer :: iland

  ! first define all countries as undefined, just in case
  do iland=1,NLAND
    Country(iland) = cc(  "N/A" , iland ,F,  17 , 0  , "Not_defined                   " )
  enddo

  !--------------  code  index sea timefac_index timezone  Name  ------------!
  
Country( IC_AL ) = cc(  "AL " ,  1 ,F,  1,  1  , "Albania                       " )
Country( IC_AT ) = cc(  "AT " ,  2 ,F,  2,  1  , "Austria                       " )
Country( IC_BE ) = cc(  "BE " ,  3 ,F,  3,  1  , "Belgium                       " )
Country( IC_BG ) = cc(  "BG " ,  4 ,F,  4,  2  , "Bulgaria                      " )
Country( IC_CS ) = cc(  "CS " ,  5 ,F,  5,  1  , "Former Czechoslovakia         " )
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
Country( IC_PT ) = cc(  "PT " , 20 ,F, 20,  1  , "Portugal                      " )
Country( IC_RO ) = cc(  "RO " , 21 ,F, 21,  2  , "Romania                       " )
Country( IC_ES ) = cc(  "ES " , 22 ,F, 22,  1  , "Spain                         " )
Country( IC_SE ) = cc(  "SE " , 23 ,F, 23,  1  , "Sweden                        " )
Country( IC_CH ) = cc(  "CH " , 24 ,F, 24,  1  , "Switzerland                   " )
Country( IC_TR ) = cc(  "TR " , 25 ,F, 25,  2  , "Turkey                        " )
Country( IC_SU ) = cc(  "SU " , 26 ,F, 26,  3  , "Former USSR                   " )
Country( IC_GB ) = cc(  "GB " , 27 ,F, 27,  0  , "United Kingdom                " )
Country( IC_VUL) = cc(  "VUL" , 28 ,F, 28,  1  , "Volcanoes                     " )
Country( IC_REM) = cc(  "REM" , 29 ,F, 29,  1  , "Remaining land areas          " )
Country( IC_BAS) = cc(  "BAS" , 30 ,T, 30,  1  , "The Baltic Sea                " )
Country( IC_NOS) = cc(  "NOS" , 31 ,T, 31,  1  , "The North Sea                 " )
Country( IC_ATL) = cc(  "ATL" , 32 ,T, 32,  1  , "Remaining NE Atlantic Ocean   " )
Country( IC_MED) = cc(  "MED" , 33 ,T, 33,  1  , "The Mediterranean Sea         " )
Country( IC_BLS) = cc(  "BLS" , 34 ,T, 34,  1  , "The Black Sea                 " )
Country( IC_NAT) = cc(  "NAT" , 35 ,F, 35,  1  , "Natural marine sources        " )
Country( IC_RUO) = cc(  "RUO" , 36 ,F, 36,  3  , "Kola/Karelia                  " )
Country( IC_RUP) = cc(  "RUP" , 37 ,F, 37,  3  , "St.Petersburg/Novgorod-Pskov  " )
Country( IC_RUA) = cc(  "RUA" , 38 ,F, 38,  3  , "Kaliningrad                   " )
Country( IC_BY ) = cc(  "BY " , 39 ,F, 39,  2  , "Belarus                       " )
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
Country( IC_YU ) = cc(  "YU " , 51 ,F, 51,  1  , "Yugoslavia                    " )
Country( IC_MK ) = cc(  "MK " , 52 ,F, 52,  1  , "Macedonia, The F.Yugo.Rep. of " )
Country( IC_KZ ) = cc(  "KZ " , 53 ,F, 53,  6  , "Kazakstan                     " )
Country( IC_GE ) = cc(  "GE " , 54 ,F, 54,  4  , "Georgia                       " )
Country( IC_CY ) = cc(  "CY " , 55 ,F, 55,  2  , "Cyprus                        " )
Country( IC_AM ) = cc(  "AM " , 56 ,F, 56,  4  , "Armenia                       " )
Country( IC_MT ) = cc(  "MT " , 57 ,F, 57,  1  , "Malta                         " )
Country( IC_ASI) = cc(  "ASI" , 58 ,F, 58,  0  , "Other Asian areas             " )
Country( IC_LI ) = cc(  "LI " , 59 ,F, 59,  1  , "Lichtenstein                  " )
Country( IC_DE ) = cc(  "DE " , 60 ,F, 60,  1  , "Germany                       " )
Country( IC_RU ) = cc(  "RU " , 61 ,F, 61,  3  , "Russian Federation            " )
Country( IC_MC ) = cc(  "MC " , 62 ,F, 62,  1  , "Monaco                        " )
Country( IC_NOA) = cc(  "NOA" , 63 ,F, 63,  1  , "North Africa                  " )
Country( IC_EU ) = cc(  "EU " , 64 ,F, 64,  1  , "European Community            " )
Country( IC_US ) = cc(  "US " , 65 ,F, 65,  1  , "USA                           " )
Country( IC_CA ) = cc(  "CA " , 66 ,F, 66,  1  , "Canada                        " )
Country( IC_DUMMY1 ) &
                 = cc(  "N/A" , 67 ,F,  67, 0  , "Not_defined                   " )
Country( IC_KG ) = cc(  "KG " , 68 ,F,  68, 6  , "Kyrgyzstan                    " )
Country( IC_AZ ) = cc(  "AZ " , 69 ,F,  69, 3  , "Azerbaijan                    " )
Country( IC_ATX) = cc(  "ATX" , 70 ,T,  32, 1  , "Atlantic outside. EMEP        " )
Country( IC_RUX) = cc(  "RUX" , 71 ,F,  42, 4  , "Russian Fed. outside emep     " )
Country( IC_RS)  = cc(  "RS " , 72 ,F,  72, 1  , "Serbia                        " )
Country( IC_ME)  = cc(  "ME " , 73 ,F,  73, 1  , "Montenegro                    " )


! Sea areas splitt according to innside/outside 12 nautical mile zone, 
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

Country( IC_BL2 ) = cc(  "BL2" ,342 ,T,  34, 1  , "B. Sea EU cargo outs.12      " )
Country( IC_BL3 ) = cc(  "BL3" ,343 ,T,  34, 1  , "B. Sea ROW cargo outs. 12    " )
Country( IC_BL4 ) = cc(  "BL4" ,344 ,T,  34, 1  , "B. Sea EU cargo ins. 12      " )
Country( IC_BL5 ) = cc(  "BL5" ,345 ,T,  34, 1  , "B. Sea ROW cargo ins. 12     " )
Country( IC_BL6 ) = cc(  "BL6" ,346 ,T,  34, 1  , "B. Sea EU ferries outs.12    " )
Country( IC_BL7 ) = cc(  "BL7" ,347 ,T,  34, 1  , "B. Sea ROW ferries outs. 12  " )
Country( IC_BL8 ) = cc(  "BL8" ,348 ,T,  34, 1  , "B. Sea EU ferries ins. 12    " )
Country( IC_BL9 ) = cc(  "BL9" ,349 ,T,  34, 1  , "B. Sea ROW ferries ins. 12   " )


  end subroutine Country_Init
end module Country_ml
