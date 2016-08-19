! <My_Chem_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!########### OZONE model ###################################
!>_________________________________________________________<

  module  GenSpec_bgn_ml
!-----------------------------------------------------------
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP
  implicit none
  private

!+ Defines indices and NSPEC for bgn : Background species

 ! Species which can be specified simply for each column, e.g.
 ! as function of local meteorology or zenith angle
 !   o2, m,  and for MADE-like, oh, ch3coo2

   integer, public, parameter ::  NSPEC_BGN = 0 ! No. 3D bgn species
   integer, public, parameter ::  NSPEC_COL = 0 ! total no. prescribed specs

  !/ define xn_2d_bgn here.
   real, public, save, dimension(1,KCHEMTOP:KMAX_MID) :: xn_2d_bgn

!-----------------------------------------------------------
  end module GenSpec_bgn_ml


!>_________________________________________________________<

  module  GenSpec_adv_ml
!-----------------------------------------------------------
  implicit none
!/ ...... ..   ( from GenChem )

!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 60
 
 ! Aerosols:

   logical, public, parameter ::  ORG_AEROSOLS = .false. 

   integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  



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
     IXADV_CH2CCH3     =  10   &
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
  ,  IXADV_ISOP        =  27   &
  ,  IXADV_CH3O2H      =  28   &
  ,  IXADV_C2H5OOH     =  29

   integer, public, parameter ::   & 
     IXADV_BURO2H      =  30   &
  ,  IXADV_ETRO2H      =  31   &
  ,  IXADV_PRRO2H      =  32   &
  ,  IXADV_OXYO2H      =  33   &
  ,  IXADV_MEKO2H      =  34   &
  ,  IXADV_MALO2H      =  35   &
  ,  IXADV_MVKO2H      =  36   &
  ,  IXADV_MARO2H      =  37   &
  ,  IXADV_ISRO2H      =  38   &
  ,  IXADV_H2O2        =  39

   integer, public, parameter ::   & 
     IXADV_CH3COO2H    =  40   &
  ,  IXADV_CH2CO2HCH3  =  41   &
  ,  IXADV_ISONO3H     =  42   &
  ,  IXADV_ISNIRH      =  43   &
  ,  IXADV_CH3OH       =  44   &
  ,  IXADV_C2H5OH      =  45   &
  ,  IXADV_H2          =  46   &
  ,  IXADV_CO          =  47   &
  ,  IXADV_CH4         =  48   &
  ,  IXADV_SO2         =  49

   integer, public, parameter ::   & 
     IXADV_SO4         =  50   &
  ,  IXADV_pNO3        =  51   &
  ,  IXADV_NH3         =  52   &
  ,  IXADV_aNH4        =   53   & !total NH4
  ,  IXADV_aNO3        =   54   & !total fine particulate nitrate
  ,  IXADV_PM25        =   55  &
  ,  IXADV_PMco        =   56  &
  ,  IXADV_SSfi        =   57  &  !SeaSalt 
  ,  IXADV_SSco        =   58  & 
  ,  IXADV_Rn222       =   59  &  !
  ,  IXADV_Pb210       =   60     !
 !-----------------------------------------------------------
  end module GenSpec_adv_ml
!>_________________________________________________________<

  module  GenSpec_shl_ml
!-----------------------------------------------------------
  implicit none
!/ ...... ..   ( from GenChem )

!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 16 
 
 ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  



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
  ,  IXSHL_MACRO2      =  15 &
  ,  IXSHL_PHNO3       =  16 
 !-----------------------------------------------------------
  end module GenSpec_shl_ml
!>_________________________________________________________<

  module  GenSpec_tot_ml
!-----------------------------------------------------------

  
  implicit none
!/ ...... ..   ( from GenChem )

!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   logical, public, parameter ::  ORG_AEROSOLS = .false. 

   integer, public, parameter ::  NSPEC_TOT = 76 

 
 ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  



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
  ,  PHNO3       =  16   &
  ,  O3          =  17   &
  ,  NO          =  18   &
  ,  NO2         =  19   &
  ,  PAN         =  20

   integer, public, parameter ::   & 
     MPAN        =  21   &
  ,  NO3         =  22   &
  ,  N2O5        =  23   &
  ,  ISONO3      =  24   &
  ,  HNO3        =  25   &
  ,  CH2CCH3     =  26   &
  ,  CH3COO2     =  27   &
  ,  MACR        =  28   &
  ,  ISNI        =  29   &
  ,  ISNIR       =  30

   integer, public, parameter ::   & 
     GLYOX       =  31   &
  ,  MGLYOX      =  32   &
  ,  MAL         =  33   &
  ,  MEK         =  34   &
  ,  MVK         =  35   &
  ,  HCHO        =  36   &
  ,  CH3CHO      =  37   &
  ,  C2H6        =  38   &
  ,  NC4H10      =  39   &
  ,  C2H4        =  40

   integer, public, parameter ::   & 
     C3H6        =  41   &
  ,  OXYL        =  42   &
  ,  ISOP        =  43   &
  ,  CH3O2H      =  44   &
  ,  C2H5OOH     =  45   &
  ,  BURO2H      =  46   &
  ,  ETRO2H      =  47   &
  ,  PRRO2H      =  48   &
  ,  OXYO2H      =  49   &
  ,  MEKO2H      =  50

   integer, public, parameter ::   & 
     MALO2H      =  51   &
  ,  MVKO2H      =  52   &
  ,  MARO2H      =  53   &
  ,  ISRO2H      =  54   &
  ,  H2O2        =  55   &
  ,  CH3COO2H    =  56   &
  ,  CH2CO2HCH3  =  57   &
  ,  ISONO3H     =  58   &
  ,  ISNIRH      =  59   &
  ,  CH3OH       =  60

   integer, public, parameter ::   & 
     C2H5OH      =  61   &
  ,  H2          =  62   &
  ,  CO          =  63   &
  ,  CH4         =  64   &
  ,  SO2         =  65   &
  ,  SO4         =  66   &
  ,  pNO3        =  67   &
  ,  NH3         =  68   &
  ,  aNH4        =   69   &
  ,  aNO3        =   70  &
  ,  PM25        =   71  &
  ,  PMco        =   72  &
  ,  SSFI        =   73  &    !SeaS
  ,  SSco        =   74  &    !SeaS
  ,  Rn222       =   75  &
  ,  Pb210       =   76       !ds apr2005

 !-----------------------------------------------------------
  end module GenSpec_tot_ml
!>_________________________________________________________<
!>_________________________________________________________<

  module  GenChemicals_ml
!-----------------------------------------------------------

   use GenSpec_tot_ml, only: NSPEC_TOT    ! Total number of species for chemistry
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !/--   Characteristics of species: 
  !/--   Number, name, molwt, carbon num, nmhc (1) or not(0)
 
  public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc

  type, public :: Chemical 
       character(len=12) :: name
       real              :: molwt       
       integer           :: nmhc      ! nmhc (1) or not(0)
       integer           :: carbons   ! Carbon-number
       real              :: nitrogens ! Nitrogen-number
       integer           :: sulphurs  ! Sulphur-number
  end type Chemical
  type(Chemical), public, dimension(NSPEC_TOT) :: species

  contains
    subroutine define_chemicals()
    !+
    ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
    ! array, using indices from total list of species (advected + short-lived).
    !
    !                                           MW NMHC  C    N   S
       species(  1) = Chemical("OD          ",  16,  0,  0,   0,  0 ) 
       species(  2) = Chemical("OP          ",  16,  0,  0,   0,  0 ) 
       species(  3) = Chemical("OH          ",  17,  0,  0,   0,  0 ) 
       species(  4) = Chemical("HO2         ",  33,  0,  0,   0,  0 ) 
       species(  5) = Chemical("CH3O2       ",  47,  0,  1,   0,  0 ) 
       species(  6) = Chemical("C2H5O2      ",  61,  0,  2,   0,  0 ) 
       species(  7) = Chemical("SECC4H9O2   ",  89,  0,  4,   0,  0 ) 
       species(  8) = Chemical("ISRO2       ", 101,  0,  5,   0,  0 ) 
       species(  9) = Chemical("ETRO2       ",  77,  0,  2,   0,  0 ) 
       species( 10) = Chemical("PRRO2       ",  91,  0,  3,   0,  0 ) 
       species( 11) = Chemical("OXYO2       ",   1,  0,  8,   0,  0 ) !ds1
       species( 12) = Chemical("MEKO2       ", 103,  0,  4,   0,  0 ) 
       species( 13) = Chemical("MALO2       ", 147,  0,  5,   0,  0 ) 
       species( 14) = Chemical("MVKO2       ", 119,  0,  4,   0,  0 ) 
       species( 15) = Chemical("MACRO2      ", 102,  0,  4,   0,  0 ) 
       species( 16) = Chemical("PHNO3       ",   0,  0,  0,   0,  0 ) 
       species( 17) = Chemical("O3          ",  48,  0,  0,   0,  0 ) 
       species( 18) = Chemical("NO          ",  30,  0,  0,   1,  0 ) 
       species( 19) = Chemical("NO2         ",  46,  0,  0,   1,  0 ) 
       species( 20) = Chemical("PAN         ", 121,  0,  2,   1,  0 ) 
       species( 21) = Chemical("MPAN        ", 132,  0,  4,   1,  0 ) 
       species( 22) = Chemical("NO3         ",  62,  0,  0,   1,  0 ) 
       species( 23) = Chemical("N2O5        ", 108,  0,  0,   2,  0 ) 
       species( 24) = Chemical("ISONO3      ", 110,  0,  5,   1,  1 ) 
       species( 25) = Chemical("HNO3        ",  63,  0,  0,   1,  0 ) 
       species( 26) = Chemical("CH2CCH3     ",  73,  0,  3,   0,  0 ) 
       species( 27) = Chemical("CH3COO2     ",  75,  0,  2,   0,  0 ) 
       species( 28) = Chemical("MACR        ",  70,  0,  4,   0,  0 ) 
       species( 29) = Chemical("ISNI        ",  46,  0,  4,   1,  1 ) 
       species( 30) = Chemical("ISNIR       ",  46,  0,  4,   1,  1 )  
       species( 31) = Chemical("GLYOX       ",  58,  0,  2,   0,  0 ) 
       species( 32) = Chemical("MGLYOX      ",  72,  0,  3,   0,  0 ) 
       species( 33) = Chemical("MAL         ",  98,  0,  5,   0,  0 ) 
       species( 34) = Chemical("MEK         ",  72,  0,  4,   0,  0 ) 
       species( 35) = Chemical("MVK         ",  70,  0,  4,   0,  0 ) 
       species( 36) = Chemical("HCHO        ",  30,  0,  1,   0,  0 ) 
       species( 37) = Chemical("CH3CHO      ",  44,  0,  2,   0,  0 ) 
       species( 38) = Chemical("C2H6        ",  30,  1,  2,   0,  0 ) 
       species( 39) = Chemical("NC4H10      ",  58,  1,  4,   0,  0 ) 
       species( 40) = Chemical("C2H4        ",  28,  1,  2,   0,  0 ) 
       species( 41) = Chemical("C3H6        ",  42,  1,  3,   0,  0 ) 
       species( 42) = Chemical("OXYL        ", 106,  1,  8,   0,  0 ) 
       species( 43) = Chemical("ISOP        ",  68,  1,  5,   0,  0 ) 
       species( 44) = Chemical("CH3O2H      ",  48,  0,  1,   0,  0 ) 
       species( 45) = Chemical("C2H5OOH     ",  62,  0,  2,   0,  0 ) 
       species( 46) = Chemical("BURO2H      ",  90,  0,  4,   0,  0 ) 
       species( 47) = Chemical("ETRO2H      ",  78,  0,  2,   0,  0 ) 
       species( 48) = Chemical("PRRO2H      ",  92,  0,  3,   0,  0 ) 
       species( 49) = Chemical("OXYO2H      ",   1,  0,  8,   0,  0 ) 
       species( 50) = Chemical("MEKO2H      ", 104,  0,  4,   0,  0 ) 
       species( 51) = Chemical("MALO2H      ", 147,  0,  5,   0,  0 ) 
       species( 52) = Chemical("MVKO2H      ",   1,  0,  4,   0,  0 )
       species( 53) = Chemical("MARO2H      ",   1,  0,  5,   0,  0 ) 
       species( 54) = Chemical("ISRO2H      ",   1,  0,  5,   0,  0 )
       species( 55) = Chemical("H2O2        ",  34,  0,  0,   0,  0 ) 
       species( 56) = Chemical("CH3COO2H    ",  76,  0,  2,   0,  0 ) 
       species( 57) = Chemical("CH2CO2HCH3  ",  74,  0,  3,   0,  0 ) 
       species( 58) = Chemical("ISONO3H     ",   1,  0,  5,   0,  0 ) 
       species( 59) = Chemical("ISNIRH      ",   1,  0,  5,   0,  0 )
       species( 60) = Chemical("CH3OH       ",  32,  0,  1,   0,  0 ) 
       species( 61) = Chemical("C2H5OH      ",  46,  0,  2,   0,  0 ) 
       species( 62) = Chemical("H2          ",   2,  0,  0,   0,  0 ) 
       species( 63) = Chemical("CO          ",  28,  0,  1,   0,  0 ) 
       species( 64) = Chemical("CH4         ",  16,  0,  1,   0,  0 ) 
       species( 65) = Chemical("SO2         ",  64,  0,  0,   0,  1 ) 
       species( 66) = Chemical("SO4         ",  96,  0,  0,   0,  1 ) 
       species( 67) = Chemical("pNO3        ",  62,  0,  0,   1,  0 ) 
       species( 68) = Chemical("NH3         ",  17,  0,  0,   1,  0 ) 
       species( 69) = Chemical("aNH4        ", 18,   0,  0,   1,  0 ) 
       species( 70) = Chemical("aNO3        ", 62,   0,  0,   1,  0 ) 
       species( 71) = Chemical("PM25        ", 100,  0,  0,   0,  0 ) 
       species( 72) = Chemical("PMCO        ", 100,  0,  0,   0,  0 ) 
       species( 73) = Chemical("SSfi        ", 58,   0,  0,   0,  0 )  !SeaS
       species( 74) = Chemical("SSco        ", 58,   0,  0,   0,  0 )
       species( 75) = Chemical("Rn222       ", 222,   0,  0,   0,  0 )
       species( 76) = Chemical("Pb210       ", 210,   0,  0,   0,  0 )

   end subroutine define_chemicals
 end module GenChemicals_ml
 !-----------------------------------------------------------
!>_________________________________________________________<

  module  GenRates_rcmisc_ml
!-----------------------------------------------------------

  
  use PhysicalConstants_ml,  only : PI, RGAS_J
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP
!                                      VOLFAC        ! for N2O5-> NO3-

  use Functions_ml,          only : troe
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !+ Tabulates Rate-coefficients - complex dependancies 

    public :: set_rcmisc_rates

    integer, parameter, public :: NRCMISC = 18   !! No. coefficients

    real, save, public, dimension(NRCMISC) :: rcvmisc 
    real, save, public, dimension(366)  :: &
      tab_so2ox   ! Tabulated so2->so4 rate for 366 days (leap-year safe!)

  contains
  !------------------------------------
  subroutine set_rcmisc_rates(itemp,tinv,m,o2,h2o,rh,rcmisc) 
  integer, intent(in), dimension(KCHEMTOP:KMAX_MID) :: itemp
  real, intent(in), dimension(KCHEMTOP:KMAX_MID) :: tinv,m,o2,h2o,rh
  real, intent(out),dimension(NRCMISC,KCHEMTOP:KMAX_MID) :: rcmisc
  integer :: k ! local   
  real,  dimension(KCHEMTOP:KMAX_MID) ::  n2   ! nitrogen
  real :: lt3(KCHEMTOP:KMAX_MID)   !  - for Troe
  n2 = m - o2
 
       rcmisc(1,:) = 6.0e-34*m*o2*(300.0*tinv)**2.3 
       rcmisc(2,:) = 1.8e-11*n2*exp(107.0*tinv) 
       rcmisc(3,:) = 3.2e-11*o2*exp(67.0*tinv) 
       rcmisc(4,:) = 2.2e-10*h2o 
       rcmisc(5,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*2.3e-13*exp(600.0*tinv) 
       rcmisc(6,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*1.7e-33*exp(1000.0*tinv)*m 
       rcmisc(7,:) = 1.3e-13*(1+0.6*m/2.55e19) 


       do k = KCHEMTOP, KMAX_MID
          if ( rh(k) > 0.4) then
            rcmisc(8,k) = &
                sqrt(3.0 * RGAS_J * itemp(k) / 0.108) & ! mean molecular speed,m/s !
            /(4*(2.5 - rh(k)*1.25))!density, corrected for rh (moderate approx.)
                                   !VOLFAC now in My_Reactions
          else
            rcmisc(8,k) = 0.0
          endif

     if (rh(k) > 0.9 ) then
         rcmisc(10,k) = 1.0e-4
     else
          rcmisc(10,k) = 5.0e-6
     end if
       end do ! k

    ! - new SO2 -> SO4 method from old ACID code


    ! - troe stuff put here to simplify .....

  lt3(:) = log(300.0*tinv(:))


  rcmisc(11,:) = troe(1.0e-31*exp(1.6*lt3(:)),3.0e-11*exp(-0.3*lt3(:)), -0.1625,m(:))
  rcmisc(12,:) = troe(2.7e-30*exp(3.4*lt3(:)),2.0e-12*exp(-0.2*lt3(:)),  -1.109,m(:))
  rcmisc(13,:) = troe(1.0e-3*exp(3.5*lt3(:))*exp(-11000*tinv(:)),9.70e14*exp(-0.1*lt3(:))*exp(-11080*tinv(:)),  -1.109,m(:)) 
    rcmisc(14,:) = troe(2.6e-30*exp(2.9*lt3(:)),6.7e-11*exp(0.6*lt3(:)), -0.844,m(:))
    rcmisc(15,:) = troe(2.7e-28*exp(7.1*lt3(:)),1.2e-11*exp(0.1*lt3(:)), -1.204,m(:)) 
    rcmisc(16,:) = 1*troe(4.9e-3*exp(-12100*tinv(:)),5.4e16*exp(-13830*tinv(:)),  -1.204,m(:)) 
    rcmisc(17,:) = troe(7.0e-29*exp(3.1*lt3(:)),9.0e-12, -0.3567,m(:)) 
    rcmisc(18,:) = troe(8.0e-17*exp(3.5*lt3(:)),3.0e-11,-0.6931,m(:)) 

  end subroutine set_rcmisc_rates
end module  GenRates_rcmisc_ml
!>_________________________________________________________<

  module  GenRates_rct_ml
!-----------------------------------------------------------
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP &
                                    , CHEMTMIN, CHEMTMAX !u3
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !+ Tabulates Rate-coefficients - temperature dependant 

    public :: set_rct_rates, set_night_rct

    integer, parameter, public :: NRCT = 37   !! No. coefficients

    real, save, public, dimension(NRCT) :: rcvt 

!/ Output gas-phase chemical rates:   ! - from Tabulations

  real, save, public, &
         dimension(NRCT,CHEMTMIN:CHEMTMAX) :: rcit  ! rate-coefficients

!- added for ozone model also (consistency with ACID)
! Only nighttime NO2->NO3
   logical, public, parameter ::  ONLY_NIGHT = .false. 
       

  contains
  !------------------------------------
  subroutine set_rct_rates(tinv) 
  real, intent(in) :: tinv
       rcvt(1) = 1.8e-12*exp(-1370.0*tinv) 
       rcvt(2) = 1.2e-13*exp(-2450.0*tinv) 
       rcvt(3) = 1.9e-12*exp(-1000.0*tinv) 
       rcvt(4) = 1.4e-14*exp(-600.0*tinv) 
       rcvt(5) = 1.8e-11*exp(110.0*tinv) 
       rcvt(6) = 3.7e-12*exp(240.0*tinv) 
       rcvt(7) = 7.2e-14*exp(-1414.0*tinv) 
       rcvt(8) = 4.8e-11*exp(250.0*tinv) 
       rcvt(9) = 2.9e-12*exp(-160.0*tinv) 
       rcvt(10) = 7.7e-12*exp(-2100.0*tinv) 
       rcvt(11) = 1.05e-14*exp(785.0*tinv) 
       rcvt(12) = 3.9e-12*exp(-1765.0*tinv) 
       rcvt(13) = 4.2e-12*exp(180.0*tinv) 
       rcvt(14) = 5.9e-14*exp(509.0*tinv) 
       rcvt(15) = 7.04e-14*exp(365.0*tinv) 
       rcvt(16) = 3.1e-12*exp(-360.0*tinv) 
       rcvt(17) = 3.8e-13*exp(780.0*tinv) 
       rcvt(18) = 1e-12*exp(190.0*tinv) 
       rcvt(19) = 1.9e-12*exp(190.0*tinv) 
       rcvt(20) = 8.6e-12*exp(20.0*tinv) 
       rcvt(21) = 7.9e-12*exp(-1030.0*tinv) 
       rcvt(22) = 2.7e-13*exp(1000.0*tinv) 
       rcvt(23) = 5.8e-12*exp(190.0*tinv) 
       rcvt(24) = 5.6e-12*exp(310.0*tinv) 
       rcvt(25) = 2.8e-12*exp(530*tinv) 
       rcvt(26) = 1.3e-13*exp(1040.0*tinv) 
       rcvt(27) = 3e-13*exp(1040.0*tinv) 
       rcvt(28) = 3.69e-12*exp(-70*tinv) 
       rcvt(29) = 1.64e-11*exp(-559.0*tinv) 
       rcvt(30) = 1.2e-14*exp(-2630.0*tinv) 
       rcvt(31) = 6.5e-15*exp(-1880.0*tinv) 
       rcvt(32) = 1.23e-14*exp(-2013*tinv) 
       rcvt(33) = 2.54e-11*exp(410.0*tinv) 
       rcvt(34) = 4.13e-12*exp(452.0*tinv) 
       rcvt(35) = 1.86e-11*exp(175.0*tinv) 
       rcvt(36) = 1.34e+16*exp(-13330.0*tinv) 
       rcvt(37) = 4.32e-15*exp(-2016.0*tinv) 

  end subroutine set_rct_rates
  !------------------------------------------------------
  subroutine set_night_rct(rct,rh,i,j)
  implicit none
  integer,intent(in) :: i,j
  real,intent(in) :: rct(NRCT,KCHEMTOP:KMAX_MID)
  real,intent(in)    :: rh(KCHEMTOP:KMAX_MID)

   ! Dummy for OZONE 

  end subroutine set_night_rct
  !------------------------------------------------------
end module GenRates_rct_ml

!>_________________________________________________________<

  module  MyChem_ml
!-----------------------------------------------------------
! Module containijng initial setup routine calls (Init_mychem)
! and intended to allow the user to specify miscelanneaous
! bits of extra code as needed. Here we have so far included
! Set_2dBgnd in orer to get xn_2d:bgnd for MADE.
! 
! We have a new subroutine Init_mychem for all model versions
! which now does tabulations previously done in Tabulations_ml

  use Functions_ml,   only : Daily_sine  ! to specify so2ox
  use GenSpec_bgn_ml,        only : NSPEC_COL ! - nothing more needed
                                   ! for OZONE , xn_2d_bgn, IXBGN_OH, 

 use GenRates_rct_ml, only : &  !
                 NRCT, &         ! No. temperature dependant coefficients
                 rcvt, &         ! Temperature dependant coefficients
                 rcit, &         ! Rate coeffs as rc(n, temp(k) )
                 set_rct_rates   ! Gives RCT as function of temp t

 use GenRates_rcmisc_ml, only : tab_so2ox

  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP, KCLOUDTOP &
                                     ,CHEMTMIN, CHEMTMAX  !u3 temp. range
  use PhysicalConstants_ml,  only : PI, DEG2RAD
  implicit none
  private
                                    !depending on clouds

  public :: Init_mychem          ! Calls model-specific routines
  public :: Set_2dBgnd   ! Sets model-specific background concs.
                         ! (dummy for OZONE so far)


  contains
    !------------------------------------------------------------------

    subroutine  Init_mychem()

    !+1) Temperature-dependant rates (rct). Only needs to be called once
    !    at beginning of simulations to set up table

      integer :: it          !   Local loop variable
      real    ::  tinv       !   temperature in K

      do it = CHEMTMIN, CHEMTMAX
        tinv = 1.0/real(it)
        call set_rct_rates(tinv)
        rcit(:,it) = rcvt(:)
      end do

     !+2) 
     ! Tabulate SO2 oxidation rates with a safe 366 value for ndays
     ! Coefficients taken from Eliassen+Saltbones (1983) (also in
     ! Berge and Jakobsen, 1998

        tab_so2ox = Daily_sine(4.0e-6,2.5e-6,80+91,366)
      

    end subroutine  Init_mychem
    !------------------------------------------------------------------

    subroutine Set_2dBgnd(izen,cloud,m)
      integer, intent(in) :: izen
      real,dimension(KMAX_MID), intent(in) :: cloud ! cloud-cover fraction
      real, intent(in), dimension(KCHEMTOP:KMAX_MID) :: m ! air density

       ! Dummy for OZONE 
    end subroutine Set_2dBgnd

 end module MyChem_ml
 !-----------------------------------------------------------
!>_________________________________________________________<
