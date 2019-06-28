! <AOD_PM_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module AOD_PM_mod
!-----------------------------------------------------------------------! 
! Calculates Aerosol Optical Depth (AOD) for 0.5 um radiation based on 
! aerosol mass concentrations and specific extinction cross-sections 
! based on Tegen et al. (1997) and Kinne et al. (2005).
! Humidity dependent specific extinction (or mass extinction efficiency)
! according to Chin et. al (2002).
!-----------------------------------------------------------------------!
! References
! Tegen et al. (1997): J. Geophys. Res., 102, 23895-23915,
!   doi:10.1029/97JD01864
! Kinne et al. (2005): Atmos. Chem. Phys., 6, 1815-1834,
!   doi:10.5194/acp-6-1815-2006, 2006
! Chin et. al (2002): J. Atm.Sci., 59, 461-483,
!   doi:10.1175/1520-0469(2002)059%3C0461:TAOTFT%3E2.0.CO;2
!-----------------------------------------------------------------------!
use ChemSpecs_mod
use Chemfields_mod,        only: AOD, Extin_coeff
use ChemGroups_mod,        only: chemgroups_maps
use CheckStop_mod,         only: CheckStop
use GridValues_mod,        only: i_fdom, j_fdom
use MetFields_mod,         only: z_bnd
use Config_module,    only: KMAX_MID, KCHEMTOP, ANALYSIS, AOD_WANTED
use Par_mod,               only: LIMAX,LJMAX   ! => x, y dimensions
use PhysicalConstants_mod, only: AVOG
use ZchemData_mod,    only: xn_2d, rh
use SmallUtils_mod,        only: find_index
implicit none
private
!-----------------------------------------------------------------------!
!// Subroutines
public :: AOD_Ext,AOD_init
!// Functions
public :: Qm_grp ! mass extinction efficiency [m2/g] for any spc/group

character(len=*), parameter :: EXT_MODE="WET" ! use wet|dry extinction

! wavelengths for AOD/extinction calcualtions
integer, parameter, public :: &
  W340=1,W350=2,W380=3,W440=4,W500=5,W550=6,W675=7,W870=8,W1020=9

character(len=*), parameter, public :: &
  wavelength(W340:W1020)=["340nm ","350nm ","380nm ","440nm ","500nm ",&
                          "550nm ","675nm ","870nm ","1020nm"]

logical, public, save :: &
  wanted_wlen(W340:W1020)=.false., & ! calculate AOD/EXT if output requires it
  wanted_ext3d=.false.               ! calculate  3D EXT if output requires it

!  Note - These dry/wet extinction prototypes are for "master" or model species.
!  They do not need to be present in the chemical scheme.
!  However, the chemical scheme needs to define after one of these prototypes.
!  If you would like other characteristics, add them here.
integer, public, parameter :: &
  DRY_MODE=0,WET_MODE=1,& ! DRY/WET extinction calculation for a specie
  NUM_CEXT=9,           & ! Total number of extinction classes (as follows)
  CEXT_DDf=1,CEXT_DDc=2,& ! Desert dust: fine,coarse
  CEXT_SSf=3,CEXT_SSc=4,& ! Sea salt: fine,coarse
  CEXT_ECn=5,CEXT_ECa=6,& ! Elem. C: new,aged
  CEXT_EC =CEXT_ECa,    & ! Elem. C: ECn/ECa only diff Qm_abs (dry)
  CEXT_OC =7,           & ! Org. C
  CEXT_SO4=8,           & ! SO4
  CEXT_NH4f=CEXT_SO4,   & ! NH4_f as SO4
  CEXT_NO3f=CEXT_SO4,   & ! NO3_f as SO4
  CEXT_NO3c=9             ! NO3_c sitting on SSc:
                          ! assume rho_dry as SO4, Q and GF as SSc
integer, private, save :: NUM_EXT=0
integer, public, save, allocatable :: aod_grp(:)
type :: ExtEffMap
  integer :: itot,cext
end type ExtEffMap
type(ExtEffMap), private, save, allocatable :: ExtMap(:)

integer, parameter :: NumRH=7
real, parameter, dimension(NumRH) ::          &
  RelHum=[0.0 ,0.5 ,0.7 ,0.8 ,0.9 ,0.95,0.99],&
  RH_CNT=[1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ]
real, parameter, dimension(NUM_CEXT,NumRH) :: &
  Gf_ref=reshape(                             &
         [RH_CNT,                             & ! 1:DDf (constant)
          RH_CNT,                             & ! 2:DDc (constant)
          1.0 ,1.6 ,1.8 ,2.0 ,2.4 ,2.9 ,4.8 , & ! 3:SSf (SS)
          1.0 ,1.6 ,1.8 ,2.0 ,2.4 ,2.9 ,4.8 , & ! 4:SSc (SS) 
          RH_CNT,                             & ! 5:ECn (constant)
          1.0 ,1.0 ,1.0 ,1.2 ,1.4 ,1.5 ,1.9 , & ! 6:ECa
          1.0 ,1.2 ,1.4 ,1.5 ,1.6 ,1.8 ,2.2 , & ! 7:OC  
          1.0 ,1.4 ,1.5 ,1.6 ,1.8 ,1.9 ,2.2 , & ! 8:SO4 
          1.0 ,1.6 ,1.8 ,2.0 ,2.4 ,2.9 ,4.8 ],& ! 9:NO3c (SSc)
          [NUM_CEXT,NumRH],order=[2,1])

real, parameter, dimension(NUM_CEXT,NumRH,W340:W1020) :: &
  Qm_ref=reshape(&
  !(wet) mass extinction efficiency [m2/g] at 340 nm
  !RH:  0%    50%    70%    80%    90%    95%    99%
   [2.4900*RH_CNT,                                    & ! 1:DDf (constant)
    2.1350*RH_CNT,                                    & ! 2:DDc (constant)
    2.5250,2.3980,2.3460,2.3060,2.2690,2.2420,2.1640, & ! 3:SSf
    2.1080,2.0760,2.0810,2.0890,2.0620,2.0550,2.0440, & ! 4:SSc
    0.9230*RH_CNT,                                    & ! 5:ECn (constant)
    0.9230,0.9238,0.9012,0.8142,0.7405,0.7216,0.7413, & ! 6:ECa
    1.2670,1.3510,1.4050,1.4580,1.5720,1.7120,2.0420, & ! 7:OC 
    2.4510,2.7020,2.7720,2.8150,2.8890,2.9690,3.1000, & ! 8:SO4 (H2SO4)
    2.1080,2.0760,2.0810,2.0890,2.0620,2.0550,2.0440, & ! 9:NO3c (SSc)
  !(wet) mass extinction efficiency [m2/g] at 350 nm for EARLINET
  !RH:  0%    50%    70%    80%    90%    95%    99%
    2.5030*RH_CNT,                                    & ! 1:DDf (constant)
    2.0940*RH_CNT,                                    & ! 2:DDc (constant)
    2.5390,2.4100,2.3600,2.3080,2.2770,2.2510,2.1650, & ! 3:SSf
    2.1030,2.0750,2.0800,2.0930,2.0630,2.0550,2.0450, & ! 4:SSc
    0.8890*RH_CNT,                                    & ! 5:ECn (constant)
    0.8890,0.8897,0.8678,0.7825,0.7090,0.6893,0.7053, & ! 6:ECa
    1.2180,1.3030,1.3570,1.4100,1.5240,1.6650,2.0010, & ! 7:OC 
    2.3790,2.6690,2.7550,2.8060,2.8930,2.9860,3.1480, & ! 8:SO4 (H2SO4)
    2.1030,2.0750,2.0800,2.0930,2.0630,2.0550,2.0450, & ! 9:NO3c (SSc)
  !(wet) mass extinction efficiency [m2/g] at 380 nm
  !RH:  0%    50%    70%    80%    90%    95%    99%
    2.5480*RH_CNT,                                    & ! 1:DDf (constant)
    2.1460*RH_CNT,                                    & ! 2:DDc (constant)
    2.5650,2.4270,2.3710,2.3480,2.3080,2.2590,2.1680, & ! 3:SSf
    2.1140,2.0880,2.0810,2.0800,2.0720,2.0620,2.0460, & ! 4:SSc
    0.8070*RH_CNT,                                    & ! 5:ECn (constant)
    0.8070,0.8077,0.7875,0.7074,0.6363,0.6158,0.6251, & ! 6:ECa
    1.0900,1.1790,1.2340,1.2860,1.3960,1.5340,1.8700, & ! 7:OC 
    2.1380,2.5000,2.6240,2.7000,2.8140,2.9360,3.1650, & ! 8:SO4 (H2SO4)
    2.1140,2.0880,2.0810,2.0800,2.0720,2.0620,2.0460, & ! 9:NO3c (SSc)
  !(wet) mass extinction efficiency [m2/g] at 440 nm
  !RH:  0%    50%    70%    80%    90%    95%    99%
    2.6100*RH_CNT,                                    & ! 1:DDf (constant)
    2.1610*RH_CNT,                                    & ! 2:DDc (constant)
    2.6230,2.4930,2.4190,2.4090,2.3370,2.2930,2.1900, & ! 3:SSf
    2.1230,2.0940,2.0840,2.0900,2.0830,2.0690,2.0490, & ! 4:SSc
    0.6657*RH_CNT,                                    & ! 5:ECn (constant)
    0.6657,0.6665,0.6493,0.5791,0.5135,0.4922,0.4908, & ! 6:ECa
    0.8612,0.9563,1.0110,1.0590,1.1600,1.2870,1.6070, & ! 7:OC 
    1.6870,2.1380,2.3170,2.4280,2.5840,2.7430,3.0650, & ! 8:SO4 (H2SO4)
    2.1230,2.0940,2.0840,2.0900,2.0830,2.0690,2.0490, & ! 9:NO3c (SSc)
  !(wet) mass extinction efficiency [m2/g] at 500 nm
  !RH:  0%    50%    70%    80%    90%    95%    99%
    2.6820*RH_CNT,                                    & ! 1:DDf (constant)
    2.1780*RH_CNT,                                    & ! 2:DDc (constant)
    2.6960,2.5080,2.5180,2.4680,2.3890,2.3130,2.1990, & ! 3:SSf
    2.1430,2.1030,2.1060,2.1110,2.0840,2.0700,2.0640, & ! 4:SSc
    0.5574*RH_CNT,                                    & ! 5:ECn (constant)
    0.5574,0.5583,0.5435,0.4822,0.4226,0.4017,0.3941, & ! 6:ECa
    0.6797,0.7744,0.8271,0.8709,0.9615,1.0750,1.3650, & ! 7:OC 
    1.3430,1.7980,1.9990,2.1190,2.2910,2.4640,2.8200, & ! 8:SO4 (H2SO4)
    2.1430,2.1030,2.1060,2.1110,2.0840,2.0700,2.0640, & ! 9:NO3c (SSc)
  !(wet) mass extinction efficiency [m2/g] at 550 nm
  !RH:  0%    50%    70%    80%    90%    95%    99%
!rep15  2.7060*RH_CNT,                                    & ! 1:DDf (constant)
!rep15  2.1890*RH_CNT,                                    & ! 2:DDc (constant)
    2.6060*RH_CNT,                                    & ! 1:DDf (constant)
    2.1260*RH_CNT,                                    & ! 2:DDc (constant)
    2.6990,2.5470,2.5440,2.5080,2.4440,2.3620,2.2210, & ! 3:SSf
    2.1430,2.1030,2.0900,2.1060,2.0840,2.0700,2.0640, & ! 4:SSc
    0.4830*RH_CNT,                                    & ! 5:ECn (constant)
    0.4830,0.4840,0.4710,0.4170,0.3630,0.3430,0.3320, & ! 6:ECa
    0.5600,0.6520,0.7010,0.7410,0.8210,0.9210,1.1810, & ! 7:OC 
    1.1140,1.5450,1.7420,1.8620,2.0360,2.2060,2.5580, & ! 8:SO4 (H2SO4)
    2.1430,2.1030,2.0900,2.1060,2.0840,2.0700,2.0640, & ! 9:NO3c (SSc)
  !(wet) mass extinction efficiency [m2/g] at 675 nm
  !RH:  0%    50%    70%    80%    90%    95%    99%
    2.7720*RH_CNT,                                    & ! 1:DDf (constant)
    2.2100*RH_CNT,                                    & ! 2:DDc (constant)
    2.7450,2.6250,2.6200,2.5970,2.5270,2.4280,2.2730, & ! 3:SSf
    2.1790,2.1330,2.1180,2.1200,2.1140,2.0960,2.0640, & ! 4:SSc
    0.3618*RH_CNT,                                    & ! 5:ECn (constant)
    0.3618,0.3629,0.3530,0.3100,0.2645,0.2463,0.2305, & ! 6:ECa
    0.3551,0.4327,0.4723,0.5028,0.5623,0.6359,0.8295, & ! 7:OC 
    0.7095,1.0630,1.2350,1.3400,1.4900,1.6350,1.9400, & ! 8:SO4 (H2SO4)
    2.1790,2.1330,2.1180,2.1200,2.1140,2.0960,2.0640, & ! 9:NO3c (SSc)
  !(wet) mass extinction efficiency [m2/g] at 870 nm
  !RH:  0%    50%    70%    80%    90%    95%    99%
    2.7750*RH_CNT,                                    & ! 1:DDf (constant)
    2.2660*RH_CNT,                                    & ! 2:DDc (constant)
    2.6990,2.6680,2.6690,2.6590,2.6110,2.5350,2.3440, & ! 3:SSf
    2.2310,2.1630,2.1630,2.1470,2.1180,2.1100,2.0740, & ! 4:SSc
    0.2570*RH_CNT,                                    & ! 5:ECn (constant)
    0.2570,0.2581,0.2510,0.2185,0.1820,0.1659,0.1471, & ! 6:ECa
    0.1848,0.2388,0.2660,0.2858,0.3232,0.3688,0.4900, & ! 7:OC 
    0.3668,0.6075,0.7312,0.8078,0.9153,1.0180,1.2350, & ! 8:SO4 (H2SO4)
    2.2310,2.1630,2.1630,2.1470,2.1180,2.1100,2.0740, & ! 9:NO3c (SSc)
  !(wet) mass extinction efficiency [m2/g] at 1020 nm
  !RH:  0%    50%    70%    80%    90%    95%    99%
    2.7280*RH_CNT,                                    & ! 1:DDf (constant)
    2.2820*RH_CNT,                                    & ! 2:DDc (constant)
    2.5810,2.6120,2.6500,2.6580,2.6450,2.5780,2.3990, & ! 3:SSf
    2.2610,2.1850,2.1780,2.1620,2.1360,2.1230,2.0850, & ! 4:SSc
    0.2173*RH_CNT,                                    & ! 5:ECn (constant)
    0.2173,0.2185,0.2125,0.1844,0.1519,0.1371,0.1180, & ! 6:ECa
    0.1336,0.1744,0.1954,0.2103,0.2380,0.2718,0.3626, & ! 7:OC 
    0.2539,0.4393,0.5375,0.5985,0.6836,0.7644,0.9356, & ! 8:SO4 (H2SO4)
    2.2610,2.1850,2.1780,2.1620,2.1360,2.1230,2.0850],& ! 9:NO3c (SSc)
    [NUM_CEXT,NumRH,W1020-W340+1],order=[2,1,3])

real,parameter,dimension(NUM_CEXT) :: &
! Qm_Dabs= & ! (dry) mass absorption efficiency [m?/g] at 550 nm
!         [0.0  ,0.0  ,0.0  ,0.0  ,8.5  ,11.5 ,0.0  ,0.0  ,0.0  ], & 
  rho_dry=[2.6  ,2.6  ,2.2  ,2.2  ,1.0  ,1.0  ,1.8  ,1.6  ,1.6  ], &
  rad_eff=[0.80 ,4.5  ,0.80 ,5.73 ,0.039,0.039,0.087,0.156,5.73 ]
         ! 1:DDf 2:DDc 3:SSf 4:SSc 5:ECn 6:ECa 7:OC  8:SO4 9:NO3c

real, pointer,dimension(:,:,:,:,:),public,save :: SpecExtCross=>null()
contains

function Qm(mode,rh,wlen,debug) result(Qm_arr)
!-----------------------------------------------------------------------!
! Calculate specific extinction (or mass extinction efficiency)
!-----------------------------------------------------------------------!
! 'DRY' and 'WET' calulation modes
!   'DRY' mode: Rh independent values from Qm_Dext array
!   'WET' mode: Rh dependent calc. according to Chin et. al (2002)
!-----------------------------------------------------------------------!
  character(len=*), intent(in) :: mode
  real, intent(in) :: rh
  integer, intent(in) :: wlen
  logical, intent(in) :: debug  
  real, dimension(NUM_EXT)  :: Qm_arr
  integer,save              :: rh_n=NumRH
  real, dimension(0:1)      :: rh_w
  real, dimension(NUM_CEXT) :: Gf,ExtEff,Qm_cext

  select case(mode)
  case("d","D","dry","Dry","DRY")   ! Dry mass extinction efficiency [m2/g]
    Gf    =Gf_ref(:,1)              !  0% RH: Growth factors
    ExtEff=Qm_ref(:,1,wlen)         !  0% RH: Extinction efficiencies

  case("w","W","wet","Wet","WET")   ! Wet mass extinction efficiency [m2/g]
    if(rh<=RelHum(1))then
      Gf    =Gf_ref(:,1)            !  0% RH: Growth factors
      ExtEff=Qm_ref(:,1,wlen)       !  0% RH: Extinction efficiencies
    elseif(rh>=RelHum(NumRH))then
      Gf    =Gf_ref(:,NumRH)        ! 99% RH: Growth factors
      ExtEff=Qm_ref(:,NumRH,wlen)   ! 99% RH: Extinction efficiencies
    else
      !.. rh interpolation weights 
      ! rh_n is updated only if the one from last call does not work
      if(rh<RelHum(rh_n-1).or.rh>RelHum(rh_n))&
        rh_n=minloc(RelHum,DIM=1,MASK=(rh<=RelHum))
      rh_w(1)=(rh-RelHum(rh_n-1))/(RelHum(rh_n)-RelHum(rh_n-1))
      rh_w(0)=1.0-rh_w(1)
      if(debug) write(*,"(a15,f8.2,'%(',i3,'):',4f8.2,'=',f8.2,'%')") &
        '## Rh    =',rh*1e2,rh_n,RelHum(rh_n-1:rh_n)*1e2,rh_w(:),&
                  sum(RelHum(rh_n-1:rh_n)*rh_w)*1e2   ! check interpolation
      !.. Interpolate
      Gf    =MATMUL(Gf_ref(:,rh_n-1:rh_n)     ,rh_w)  ! Growth factors
      ExtEff=MATMUL(Qm_ref(:,rh_n-1:rh_n,wlen),rh_w)  ! Extinction efficiencies
      if(debug) write(*,'((a15,9f10.3))') &
        '## GFs   =',gf(:),'## ExtEff=',ExtEff(:)
    end if

  case default
    call CheckStop("Unknown extinction mode: "//trim(mode))
  end select

  !.. mass extinction efficiency [m2/g]
  !beta = 3/4 * ExtEff/rho_wet/rad_eff * Mwet/Mdry
  !     = 3/4 * ExtEff/rho_wet/rad_eff * Gf^3*rho_wet/rho_dry
  !     = 3/4 * ExtEff/rad_eff * Gf^3/rho_dry
  Qm_cext=0.75* ExtEff/rad_eff * Gf**3/rho_dry

  !** Backward compatibility: NO3c characterized by
  !  rho_dry(SO4),rad_eff(SSc),Gf_ref(SSc),Qm_ref(SSc), and
  !beta = 3/4 * ExtEff(SSc)/rho_wet(NO3c)/rad_eff(SSc) * Mwet(SSc)/Mdry(SSc)
  !     = 3/4 * ExtEff(SSc)/rho_wet(NO3c)/rad_eff(SSc) * Gf(SSc)^3*rho_wet(SSc)/rho_dry(SSc)
  !     = Qm_cext(SSc)*rho_wet(SSc)/rho_wet(NO3c)
  ! with rho_wet(NO3c) derived from rho_dry(SO4) and Gf(SSc)
  Qm_cext(CEXT_NO3c)=Qm_cext(CEXT_SSc)*rho_wet(CEXT_SSc)/rho_wet(CEXT_NO3c)

  !** Backward compatibility: Different fractions within NO3c
  Qm_cext(CEXT_NO3c)=0.3*Qm_cext(CEXT_NO3f)+0.7*Qm_cext(CEXT_NO3c)

  !.. mass extinction efficiency [m2/g]
  Qm_arr(:)=Qm_cext(ExtMap%cext)

contains
function rho_wet(nc)
  integer,intent(in) :: nc          ! index: 1..NUM_CEXT
  real               :: rho_wet     ! Density of wet aerosol
  real, parameter    :: RHO_H2O=1.0 ! Water density [g/cm3]
! Vfr_dry = 1/GF**3 (dry volume fraction)
! rho_wet = Vfr_dry*rho_dry + (1-Vfr_dry)*RHO_H2O
!         = (rho_dry-RHO_H2O)/GF**3 + RHO_H2O
  rho_wet = (rho_dry(nc)-RHO_H2O)/Gf(nc)**3 + RHO_H2O
end function rho_wet
end function Qm

function Qm_grp(gtot,rh,debug) result(Qm_arr)
!-----------------------------------------------------------------------!
! Returns mass extinction efficiencies [m2/g] for any spc array/group
!-----------------------------------------------------------------------!
! Returns extinction efficiency array, one value each spc index (itot)
! in spc array/group (gtot). 
! When rh is privided (optional variable) 'WET' extinction are returned,
! 'DRY' extinction are returned otherweise.
!-----------------------------------------------------------------------!
  integer,intent(in),dimension(:) :: gtot
  real,   intent(in),optional     :: rh
  logical,intent(in),optional     :: debug
  real,     dimension(size(gtot)) :: Qm_arr
  
  integer :: n=0,i=0
  logical :: my_debug=.false.
  real, dimension(NUM_EXT) :: Qm_aux

  if(present(debug)) my_debug=debug
  if(present(rh))then
    Qm_aux=Qm("WET",rh ,W550,my_debug)
  else
    Qm_aux=Qm("DRY",0.0,W550,my_debug)
  end if

  Qm_arr(:)=0
  do n=1,size(gtot)
    i=find_index(gtot(n),ExtMap(:)%itot,debug=my_debug)
    if(i>0)Qm_arr(n)=Qm_aux(i)
  end do
end function Qm_grp

subroutine AOD_init(msg,wlen,out3d)
  character(len=*), intent(in) :: msg
  character(len=*),intent(in),optional:: wlen
  logical,intent(in),optional  :: out3d
  integer :: igrp=0,n=0
!-----------------------------------------------------------------------!
! Process optional arguments
!-----------------------------------------------------------------------!
  if(present(wlen))then
    n=find_index(wlen,wavelength)! e.g. search "550nm" on array of wavelengths
    call CheckStop(n<1,&
      trim(msg)//" Unknown AOD/EXT wavelength "//trim(wlen))
    wanted_wlen(n)=.true.
  end if
  if(present(out3d))then
    wanted_ext3d=wanted_ext3d.or.out3d
  end if
!-----------------------------------------------------------------------!
! Unpack EXTINC mapping 
!  ExtMap%itot: species 
!  ExtMap%cext: extinction prototype
!  aod_grp: ExtMap%itot for public usage
!-----------------------------------------------------------------------!
  if(.not.allocated(aod_grp))then
    igrp=find_index('EXTINC',chemgroups_maps%name)
    call CheckStop(igrp<1,trim(msg)//" EXTINC mapping not found")
    NUM_EXT=size(chemgroups_maps(igrp)%species)
    allocate(aod_grp(NUM_EXT),ExtMap(NUM_EXT))
    do n=1,NUM_EXT
      aod_grp(n)=chemgroups_maps(igrp)%species(n)
      ExtMap(n)%itot=aod_grp(n)
      select case(chemgroups_maps(igrp)%maps(n))
        case('DDf' );ExtMap(n)%cext=CEXT_DDf
        case('DDc' );ExtMap(n)%cext=CEXT_DDc
        case('SSf' );ExtMap(n)%cext=CEXT_SSf
        case('SSc' );ExtMap(n)%cext=CEXT_SSc
        case('ECn' );ExtMap(n)%cext=CEXT_ECn
        case('ECa' );ExtMap(n)%cext=CEXT_ECa
        case('EC'  );ExtMap(n)%cext=CEXT_EC
        case('OC'  );ExtMap(n)%cext=CEXT_OC
        case('SO4' );ExtMap(n)%cext=CEXT_SO4
        case('NH4f');ExtMap(n)%cext=CEXT_NH4f
        case('NO3f');ExtMap(n)%cext=CEXT_NO3f
        case('NO3c');ExtMap(n)%cext=CEXT_NO3c
        case default
          call CheckStop(trim(msg)//" Unknown EXTINC mapping "//&
                         trim(chemgroups_maps(igrp)%maps(n)))
      end select
    end do
  end if   
!-----------------------------------------------------------------------!
! Consistency checks for Qm_ref array
!-----------------------------------------------------------------------!
  call CheckStop(any(Qm_ref(:,5,W440)/=&
      [2.6100,2.1610,2.3370,2.0830,0.6657,0.5135,1.1600,2.5840,2.0830]),&
    trim(msg)//" Failed Qm_ref(90%,440nm) check")
  call CheckStop(any(Qm_ref(:,4,W550)/=&
      [2.6060,2.1260,2.5080,2.1060,0.4830,0.4170,0.7410,1.8620,2.1060]),&
    trim(msg)//" Failed Qm_ref(80%,550nm) check")
  call CheckStop(any(Qm_ref(:,6,W870)/=&
      [2.7750,2.2660,2.5350,2.1100,0.2570,0.1659,0.3688,1.0180,2.1100]),&
    trim(msg)//" Failed Qm_ref(95%,870nm) check")
!-----------------------------------------------------------------------!
! Allocate AOD arrays
!-----------------------------------------------------------------------!
  !!wanted_wlen(W550)=.true.  ! calculate 550nm for debug output
  if(any(wanted_wlen(:)).and..not.allocated(AOD))then
    allocate(AOD(NUM_EXT,LIMAX,LJMAX,W340:W1020))
    AOD=0.0
  end if
  if(wanted_ext3d.and..not.allocated(Extin_coeff))then
    allocate(Extin_coeff(NUM_EXT,LIMAX,LJMAX,KMAX_MID,W340:W1020))
    Extin_coeff=0.0
  end if 
  if(ANALYSIS.and..not.associated(SpecExtCross))then
  !!wanted_wlen(W550)=.true.  ! calculate 550nm for AOD assimilation
    allocate(SpecExtCross(NUM_EXT,KMAX_MID,LIMAX,LJMAX,W340:W1020))
    SpecExtCross=0.0
  end if
end subroutine AOD_init

subroutine AOD_Ext(i,j,debug)
!-----------------------------------------------------------------------!
! Calculates AOD on 'DRY'/'WET' mode, as defined by EXT_MODE parameter
!-----------------------------------------------------------------------!
  implicit none
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug

  logical,save             :: first_call=.true.
  integer                  :: k, n, w
  real, dimension(NUM_EXT) :: kext

!.. Debug variables for individual extinction classes
  real, dimension(NUM_CEXT):: kext_cext,AOD_cext

  if(first_call)then
    call AOD_init("AOD_Ext")
    call CheckStop(AOD_WANTED.and..not.any(wanted_wlen(:)),&
      "USE_AOR=T, but no AOD/EXT output. Check config_*.nml")
    first_call=.false.
  end if
  if(debug)then
    write(*,*) '#### in AOD module  ###'
    AOD_cext(:)=0.0
  end if
 
  !===========================================================================
  ! Extinction coefficients: 
  !  Kext [1/m] = SpecExtCross [m2/g] * mass [g/m3] summed up for all comp.
  !  xn_2d(spec,k)*1.e15 *species(ispec)%molwt/AVOG  [molec/m3] -> [ng/m3]
  !                                             [ng/m3 ] * 1e-9 -> [g/m3]
  !  => xn_2d(ispec,k) * species(ispec)%molwt * 1.e6 / AVOG  [g/m3]
  !===========================================================================
  AOD(:,i,j,:) = 0.0
  do k = KCHEMTOP, KMAX_MID
    do w=W340,W1020
      if(.not.wanted_wlen(w))cycle
      
      !.. SpecExtCross [m2/g]. EXT_MODE: use dry/wet extinction coeficients
      kext(:)=Qm(EXT_MODE,rh(k),w,debug.and.((k==KCHEMTOP+1).or.(k==KMAX_MID)))
      if(associated(SpecExtCross))&
        SpecExtCross(:,k,i,j,w)=kext(:)

      !.. Specific extinction coefficients
      kext(:)=kext(:) &                                             ! [m2/g]
        *xn_2d(ExtMap%itot,k)*species(ExtMap%itot)%molwt*1.0e6/AVOG ! [g/m3]

      !.. Extinction coefficients at level:k and wavelength:w
      if(allocated(Extin_coeff))&
        Extin_coeff(:,i,j,k,w)=kext(:)

      !.. Aerosol optical depth: integral over all vertical layers
      AOD(:,i,j,w)=AOD(:,i,j,w)+kext(:)*(z_bnd(i,j,k)-z_bnd(i,j,k+1)) ! [1/m]*[m]

      if(debug.and.w==W550)then
        !.. Extinction coefficients for diferent optical groups/types
        do n = 1,NUM_CEXT
          kext_cext(n)=sum(kext(:),MASK=(ExtMap(:)%cext==n))
        end do
        !.. Aerosol optical depth for individual components
        AOD_cext(:)=AOD_cext(:)+kext_cext(:)*(z_bnd(i,j,k)-z_bnd(i,j,k+1))

        if((k==KCHEMTOP+1).or.(k==KMAX_MID))&
          write(*,"(a8,'(',i3,')=',es10.3,'=',9(es10.3,:,'+'))") &
            'EXTINCs', k, sum(kext(:)),kext_cext(:)
      end if
    end do
  end do

  if(debug) write(*,"(a24,2i5,es10.3,'=',9(es10.3,:,'+'))") &
    '>>>  AOD / AODs  <<<', i_fdom(i), j_fdom(j), sum(AOD(:,i,j,W550)), AOD_cext(:)
end subroutine AOD_Ext
endmodule AOD_PM_mod
