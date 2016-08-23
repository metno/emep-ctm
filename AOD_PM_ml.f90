! <AOD_PM_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module AOD_PM_ml
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
use ChemSpecs
!CRM use ChemSpecs_tot_ml
!CRM use ChemChemicals_ml,     only: species  
use Chemfields_ml,        only: AOD, Extin_coeff
use ChemGroups_ml,        only: chemgroups
use CheckStop_ml,         only: CheckStop
use GridValues_ml,        only: i_fdom, j_fdom
use MetFields_ml,         only: z_bnd
use ModelConstants_ml,    only: KMAX_MID, KCHEMTOP, ANALYSIS
use Par_ml,               only: MAXLIMAX,MAXLJMAX   ! => x, y dimensions
use PhysicalConstants_ml, only: AVOG
use Setup_1dfields_ml,    only: xn_2d, rh
use SmallUtils_ml,        only: find_index

implicit none
private
!-----------------------------------------------------------------------!
!// Subroutines
public :: AOD_Ext,&
          AOD_check
!// Functions
public :: Qm_grp ! mass extinction efficiency [m2/g] for any spc/group

real, parameter :: lambda = 0.55e-6
character(len=*), parameter :: EXT_MODE="WET" ! use wet|dry extinction

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
  CEXT_EC =CEXT_ECn,    & ! Elem. C: ECn/ECa only diff Qm_dry
  CEXT_OC =7,           & ! Org. C
  CEXT_SO4=8,           & ! SO4
  CEXT_NH4f=CEXT_SO4,   & ! NH4_f as SO4
  CEXT_NO3f=CEXT_SO4,   & ! NO3_f as SO4
  CEXT_NO3c=9             ! NO3_c sitting on SSc:
                          ! assume rho_dry as SO4, Q and GF as SSc
type :: ExtEffMap
  integer :: itot,cext,mode
endtype ExtEffMap

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
include 'CM_AerExt.inc'
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

real,parameter,dimension(NUM_CEXT) :: &
 Qm_dry = & ! (dry) mass extinction efficiency [m2/g] at 550 nm
         [1.0  ,0.3  ,3.0  ,0.4  ,7.5  ,11.  ,0.0  ,8.5  ,0.0  ], & 
 rho_dry=[2.6  ,2.6  ,2.2  ,2.2  ,1.0  ,1.0  ,1.8  ,1.6  ,1.6  ], &
 rad_eff=[0.80 ,4.5  ,0.80 ,5.73 ,0.039,0.039,0.087,0.156,5.73 ]
        ! 1:DDf 2:DDc 3:SSf 4:SSc 5:ECn 6:ECa 7:OC  8:SO4 9:NO3c

integer, parameter :: NumRH=7
real, parameter, dimension(NumRH) ::    &
  RelHum = [0.0  ,0.5  ,0.7  ,0.8  ,0.9  ,0.95 ,0.99 ], &
  RH_CNT = [1.0  ,1.0  ,1.0  ,1.0  ,1.0  ,1.0  ,1.0  ]
real, parameter, dimension(NUM_CEXT,NumRH) :: &
  Gf_ref=reshape(&
           [RH_CNT,                                     & ! 1:DDf (constant)
            RH_CNT,                                     & ! 2:DDc (constant)
            1.0  ,1.6  ,1.8  ,2.0  ,2.4  ,2.9  ,4.8  ,  & ! 3:SSf (SS)
            1.0  ,1.6  ,1.8  ,2.0  ,2.4  ,2.9  ,4.8  ,  & ! 4:SSc (SS) 
            1.0  ,1.0  ,1.0  ,1.2  ,1.4  ,1.5  ,1.9  ,  & ! 5:ECn (EC)
            1.0  ,1.0  ,1.0  ,1.2  ,1.4  ,1.5  ,1.9  ,  & ! 6:ECa (EC)
            1.0  ,1.2  ,1.4  ,1.5  ,1.6  ,1.8  ,2.2  ,  & ! 7:OC  
            1.0  ,1.4  ,1.5  ,1.6  ,1.8  ,1.9  ,2.2  ,  & ! 8:SO4 
            1.0  ,1.6  ,1.8  ,2.0  ,2.4  ,2.9  ,4.8  ], & ! 9:NO3c (SSc)
            [NUM_CEXT,NumRH],order=[2,1]),&
  Qm_ref=reshape(& ! (wet) mass extinction efficiency [m2/g] at 550 nm
           [RH_CNT*Qm_dry(CEXT_DDf),                    & ! 1:DDf (constant)
            RH_CNT*Qm_dry(CEXT_DDc),                    & ! 2:DDc (constant)
            2.699,2.547,2.544,2.508,2.444,2.362,2.221,  & ! 3:SSf
            2.143,2.103,2.090,2.106,2.084,2.070,2.064,  & ! 4:SSc
            0.483,0.484,0.471,0.417,0.363,0.343,0.332,  & ! 5:ECa (EC)
            0.483,0.484,0.471,0.417,0.363,0.343,0.332,  & ! 6:ECn (EC)
            0.560,0.652,0.701,0.741,0.821,0.921,1.181,  & ! 7:OC 
            1.114,1.545,1.742,1.862,2.036,2.206,2.558,  & ! 8:SO4 (H2SO4)
            2.143,2.103,2.090,2.106,2.084,2.070,2.064], & ! 9:NO3c (SSc)
            [NUM_CEXT,NumRH],order=[2,1])

real, allocatable,dimension(:,:,:,:), public,save :: SpecExtCross
contains
! <---------------------------------------------------------->
function Qm(mode,rh,debug) result(Qm_arr)
!-----------------------------------------------------------------------!
! Calculate specific extinction (or mass extinction efficiency)
!-----------------------------------------------------------------------!
! 'DRY' and 'WET' calulation modes
!   'DRY' mode: Rh independent values from Qm_dry array
!   'WET' mode: Rh dependent calc. according to Chin et. al (2002)
!-----------------------------------------------------------------------!
  character(len=*), intent(in) :: mode
  real, intent(in) :: rh
  logical, intent(in) :: debug  
  real, dimension(NUM_EXT) :: Qm_arr
  integer       :: n
  real, dimension(0:1)  :: rh_w
  real, dimension(NUM_CEXT) :: Gf,ExtEff,Qm_wet

  select case(mode)
  case("d","D","dry","Dry","DRY") ! (dry) mass extinction efficiency [m2/g]
    Qm_arr(:)=Qm_dry(ExtMap%cext)
  case("w","W","wet","Wet","WET") ! (wet) mass extinction efficiency [m2/g]
    if(rh<=RelHum(1))then
      Gf    =Gf_ref(:,1)          ! Growth factors
      ExtEff=Qm_ref(:,1)          ! (wet) Extinction efficiencies
    elseif(rh>=RelHum(NumRH))then
      Gf    =Gf_ref(:,NumRH)      ! Growth factors
      ExtEff=Qm_ref(:,NumRH)      ! (wet) Extinction efficiencies
    else
      !.. rh interpolation weights
      RHloop: do n = 2, NumRH
        if(rh>RelHum(n)) cycle RHloop
        rh_w(1)=(rh-RelHum(n-1))/(RelHum(n)-RelHum(n-1))
        rh_w(0)=1.0-rh_w(1)
        exit RHloop
      enddo RHloop
      if(debug) write(*,"(a15,f8.2,'%(',i3,'):',4f8.2,'=',f8.2,'%')") &
        '## Rh    =',rh*1e2,n,RelHum(n-1:n)*1e2,rh_w(:),&
                      sum(RelHum(n-1:n)*rh_w)*1e2 ! check interpolation
      !.. Interpolate
      Gf    =MATMUL(Gf_ref(:,n-1:n),rh_w) ! Growth factors
      ExtEff=MATMUL(Qm_ref(:,n-1:n),rh_w) ! (wet) Extinction efficiencies
      if(debug) write(*,'((a15,9f10.3))') &
        '## GFs   =',gf(:),'## ExtEff=',ExtEff(:)
    endif

    !.. mass extinction efficiency [m2/g]
    !beta = 3/4 * ExtEff/rho_wet/rad_eff * Mwet/Mdry
    !     = 3/4 * ExtEff/rho_wet/rad_eff * Gf^3*rho_wet/rho_dry
    !     = 3/4 * ExtEff/rad_eff * Gf^3/rho_dry
    Qm_wet= 0.75* ExtEff/rad_eff * Gf**3/rho_dry

    !** Backward compatibility: NO3c characterized by
    !  rho_dry(SO4),rad_eff(SSc),Gf_ref(SSc),Qm_ref(SSc), and
    !beta = 3/4 * ExtEff(SSc)/rho_wet(NO3c)/rad_eff(SSc) * Mwet(SSc)/Mdry(SSc)
    !     = 3/4 * ExtEff(SSc)/rho_wet(NO3c)/rad_eff(SSc) * Gf(SSc)^3*rho_wet(SSc)/rho_dry(SSc)
    !     = Qm_wet(SSc)*rho_wet(SSc)/rho_wet(NO3c)
    ! with rho_wet(NO3c) derived from rho_dry(SO4) and Gf(SSc)
    Qm_wet(CEXT_NO3c)=Qm_wet(CEXT_SSc)*rho_wet(CEXT_SSc)/rho_wet(CEXT_NO3c)

    !** Backward compatibility: Different fractions within NO3c
    Qm_wet(CEXT_NO3c)=0.3*Qm_wet(CEXT_NO3f)+0.7*Qm_wet(CEXT_NO3c)

    !.. mass extinction efficiency [m2/g]
    where(ExtMap%mode==WET_MODE) 
      Qm_arr(:)=Qm_wet(ExtMap%cext)
    elsewhere(ExtMap%mode==DRY_MODE)
    !** Backward compatibility: Road Dust on "DRY" mode
      Qm_arr(:)=Qm_dry(ExtMap%cext)
    elsewhere
      Qm_arr(:)=0.0
    endwhere
  case default
    call CheckStop("Unknown extinction mode: "//trim(mode))
  endselect
contains
function rho_wet(nc)
  integer,intent(in) :: nc          ! index: 1..NUM_CEXT
  real               :: rho_wet     ! Density of wet aerosol
  real, parameter    :: RHO_H2O=1.0 ! Water density [g/cm3]
! Vfr_dry = 1/GF**3 (dry volume fraction)
! rho_wet = Vfr_dry*rho_dry + (1-Vfr_dry)*RHO_H2O
!         = (rho_dry-RHO_H2O)/GF**3 + RHO_H2O
  rho_wet = (rho_dry(nc)-RHO_H2O)/Gf(nc)**3 + RHO_H2O
endfunction rho_wet
endfunction Qm
! <---------------------------------------------------------->
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
    Qm_aux=Qm("WET",rh,my_debug)
  else
    Qm_aux=Qm("DRY",0.0,my_debug)
  endif

  Qm_arr(:)=0
  do n=1,size(gtot)
    i=find_index(gtot(n),ExtMap(:)%itot,debug=my_debug)
    if(i>0)Qm_arr(n)=Qm_aux(i)
  enddo
endfunction Qm_grp
! <---------------------------------------------------------->
subroutine AOD_check(msg)
!-----------------------------------------------------------------------!
! Consistency checks for older model versions using AOD_GROUP
!-----------------------------------------------------------------------!
  character(len=*), intent(in) :: msg
  integer :: igrp=0
  integer, pointer,dimension(:) :: aod_group=>null()
  igrp=find_index('AOD',chemgroups%name)
  if(igrp<1) return   ! AOD group no longer used... nothing to check
  aod_group=>chemgroups(igrp)%ptr
  call CheckStop(size(aod_group),NUM_EXT,&
    trim(msg)//": Incompatibe AOD_GROUP size")
  call CheckStop(any(aod_group/=ExtMap%itot),&
    trim(msg)//": Incompatibe AOD_GROUP def.")
endsubroutine AOD_check
! <---------------------------------------------------------->
subroutine AOD_Ext(i,j,debug)
!-----------------------------------------------------------------------!
! Calculates AOD on 'DRY'/'WET' mode, as defined by EXT_MODE parameter
!-----------------------------------------------------------------------!
  implicit none
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug

  logical                  :: first_call=.true.
  integer                  :: k, n
  real, dimension(NUM_EXT) :: kext

!.. Debug variables for individual extinction classes
  real, dimension(NUM_CEXT):: kext_cext,AOD_cext

  if(first_call)then
    call AOD_check("AOD_Ext")
    first_call=.false.
  endif
  if(debug)then
    write(*,*) '#### in AOD module  ###'
    AOD_cext(:)=0.0
  endif
  
  if(ANALYSIS.and..not.allocated(SpecExtCross)) &
    allocate(SpecExtCross(NUM_EXT,MAXLIMAX,MAXLJMAX,KMAX_MID))

  AOD(i,j) = 0.0
  !===========================================================================
  ! Extinction coefficients: 
  !  Kext [1/m] = SpecExtCross [m2/g] * mass [g/m3] summed up for all comp.
  !  xn_2d(spec,k)*1.e15 *species(ispec)%molwt/AVOG  [molec/m3] -> [ng/m3]
  !                                             [ng/m3 ] * 1e-9 -> [g/m3]
  !  => xn_2d(ispec,k) * species(ispec)%molwt * 1.e6 / AVOG  [g/m3]
  !===========================================================================
  do k = KCHEMTOP, KMAX_MID
    !.. SpecExtCross [m2/g]. EXT_MODE: use dry/wet extinction coeficients
    kext(:)=Qm(EXT_MODE,rh(k),debug.and.((k==KCHEMTOP+1).or.(k==KMAX_MID)))
    if(allocated(SpecExtCross)) SpecExtCross(:,i,j,k)=kext(:)

    !.. Specific extinction
    kext(:)=kext(:) &                                             ! [m2/g]
      *xn_2d(ExtMap%itot,k)*species(ExtMap%itot)%molwt*1.0e6/AVOG ! [g/m3]

    !.. Extinction coefficient on level k
    Extin_coeff(i,j,k) = sum(kext(:))

    !.. Aerosol extinction optical depth: integral over all vertical layers
    AOD(i,j) = AOD(i,j) + sum(kext(:))*(z_bnd(i,j,k)-z_bnd(i,j,k+1)) ! [1/m]*[m]

    if(debug)then
      !.. Extinction coefficients for individual components
      do n = 1,NUM_CEXT
        kext_cext(n)=sum(kext(:),ExtMap(:)%cext==n)
      enddo
      !.. Aerosol optical depth for individual components
      AOD_cext(:) = AOD_cext(:) + kext_cext(:) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))

      if((k==KCHEMTOP+1).or.(k==KMAX_MID))&
        write(*,"(a8,'(',i3,')=',es10.3,'=',9(es10.3,:,'+'))") &
          'EXTINCs', k, Extin_coeff(i,j,k), kext_cext(:)
    endif
  enddo

  if(debug) write(*,"(a24,2i5,es10.3,'=',9(es10.3,:,'+'))") &
    '>>>  AOD / AODs  <<<', i_fdom(i), j_fdom(j), AOD(i,j), AOD_cext(:)
endsubroutine AOD_Ext
endmodule AOD_PM_ml

