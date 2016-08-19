! <My_DryDep_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module My_DryDep_ml    ! DryDep_ml
!+
! Module to define the dry deposition components and rates. We
! define the min (vd_min) and max dep. velocity (Vg) first and then derive the
! daytime addition (vd_day).
!/**************************************************************************
!  Specifies which of the possible species (from Wesely's list)
!  are required in the current air pollution model   
!/**************************************************************************

 use Derived_ml,    only : f_2d,   d_2d, IOU_INST

 use GenSpec_adv_ml        !   e.g. NSPEC_ADV,IXADV_O3,IXADV_H2O2,
 use LandDefs_ml,        only : LandDefs, LandType
 use Landuse_ml,         only : WheatGrowingSeason
 use LocalVariables_ml,  only : Grid  !=> izen  integer of zenith angle
 use ModelConstants_ml , only : atwS, atwN, AOT_HORIZON
 use PhysicalConstants_ml, only : AVOG
 use SmallUtils_ml,      only: find_index
 use StoFlux_ml,         only : unit_flux, lai_flux, leaf_flux
 use TimeDate_ml,        only : current_date
 use Wesely_ml
 implicit none
 private

  public :: Init_DepMap
  public :: Add_ddep


  !/** Variables used in deposition calculations
 
  ! DDEP_xx gives the index that will be used in the EMEP model
  ! WES_xx gives the index of the Wesely gas to which this corresponds


  integer, private, save :: &
    DDEP_SOX,   DDEP_OXN,   DDEP_RDN,  &
    DDEP_OXSSW, DDEP_OXSCF, DDEP_OXSDF, DDEP_OXSCR, DDEP_OXSSN, DDEP_OXSWE, &
    DDEP_OXNSW, DDEP_OXNCF, DDEP_OXNDF, DDEP_OXNCR, DDEP_OXNSN, DDEP_OXNWE, &
    DDEP_RDNSW, DDEP_RDNCF, DDEP_RDNDF, DDEP_RDNCR, DDEP_RDNSN, DDEP_RDNWE, &
    D2_AFSTDF0, D2_AFSTDF16, D2_AFSTBF0, D2_AFSTBF16, & 
    D2_AFSTCR0, D2_AFSTCR3, D2_AFSTCR6,&
    iam_medoak, iam_beech, iam_wheat, &   ! For Fluxes
    D2_O3DF,    D2_O3WH, &
    D2_EUAOT30WH, D2_EUAOT40WH, D2_EUAOT30DF, D2_EUAOT40DF, &
    D2_UNAOT30WH, D2_UNAOT40WH, D2_UNAOT30DF, D2_UNAOT40DF, &
    D2_MMAOT40WH, D2_MMAOT30WH


  ! Here we define the minimum set of species which has different
  ! deposition velocities. We calculate Vg for these, and then
  ! can use the rates for other similar species. (e.g. AMSU can use
  ! the Vg for SO4.  Must set NDRYDEP_CALC species

  !/** IMPORTANT: the variables below must match up in the sense that, for 
  ! example, if DDEP_NH3=4 then the 4th element of DRYDEP must be WES_NH3.

  integer, public, parameter :: NDRYDEP_CALC = 10  ! gases
  integer, public, parameter :: NDRYDEP_AER = 2    ! aerosols
  integer, public, parameter :: NDRYDEP_TOT = NDRYDEP_CALC + NDRYDEP_AER


  integer, public, parameter :: &
       CDEP_HNO3 = 1, CDEP_O3  = 2, CDEP_SO2 = 3  &
      ,CDEP_NH3  = 4, CDEP_NO2 = 5, CDEP_PAN  = 6 &
      ,CDEP_H2O2 = 7, CDEP_ALD = 8, CDEP_HCHO = 9, &
       CDEP_OP = 10,  CDEP_FIN = 11, CDEP_COA = 12 

  integer, public, parameter :: CDEP_SET = -99    



 ! WE NEED A FLUX_CDEP, FLUX_ADV FOR OZONE;
 ! (set to one for non-ozone models)

  logical, public, parameter :: STO_FLUXES = .true.
  integer, public, parameter :: FLUX_CDEP  = CDEP_O3
  integer, public, parameter :: FLUX_ADV   = IXADV_O3

 
  integer, public, parameter, dimension(NDRYDEP_CALC) :: &
    DRYDEP_CALC = (/ WES_HNO3, WES_O3,   WES_SO2, &
                     WES_NH3,  WES_NO2 , WES_PAN, &
                     WES_H2O2, WES_ALD, WES_HCHO, WES_OP    /)

  !/** Compensation pount approach from CEH used?:

  logical, public, parameter :: COMPENSATION_PT = .false. 



  ! We define also the number of species which will be deposited in
  ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_CALC
  ! The actual species used and their relation to the CDEP_ indices
  ! above will be defined in Init_DepMap

  integer, public, parameter ::  NDRYDEP_ADV  = 22

  !/-- we define a type to map indices of species to be deposited
  !   to the lesser number of species where Vg is calculated

   type, public :: depmap
      integer :: adv   ! Index of species in IXADV_ arrays
      integer :: calc  ! Index of species in  calculated dep arrays
      real    :: vg    ! if CDEP_SET, give vg in m/s
   end type depmap

   type(depmap), public, dimension(NDRYDEP_ADV):: Dep

   real, public, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost


   logical, private, parameter :: MY_DEBUG = .false.

contains
  subroutine Init_DepMap
   real :: cms = 0.01     ! Convert to m/s

 ! .... Define the mapping between the advected species and
 !      the specied for which the calculation needs to be done.

   Dep(1) =  depmap( IXADV_HNO3 , CDEP_HNO3, -1.)
   Dep(2) =  depmap( IXADV_PAN,   CDEP_PAN, -1. ) 
   Dep(3) =  depmap( IXADV_NO2,   CDEP_NO2, -1. )
   Dep(4) =  depmap( IXADV_SO2,   CDEP_SO2, -1. )
   Dep(5) =  depmap( IXADV_SO4,   CDEP_FIN,  -1)
   Dep(6) =  depmap( IXADV_NH3,   CDEP_NH3, -1. )
   Dep(7) =  depmap( IXADV_aNH4,  CDEP_FIN,  -1) 
   Dep(8) =  depmap( IXADV_aNO3,  CDEP_FIN,  -1) 
   Dep(9) =  depmap( IXADV_O3   , CDEP_O3  , -1.)
   Dep(10) =  depmap( IXADV_H2O2 , CDEP_H2O2, -1.)
   Dep(11) =  depmap( IXADV_MPAN , CDEP_PAN , -1.)
   Dep(12) =  depmap( IXADV_HCHO , CDEP_HCHO, -1.)
   Dep(13) =  depmap( IXADV_CH3CHO,CDEP_ALD , -1.)
   Dep(14) =  depmap( IXADV_MAL   ,CDEP_ALD , -1.)
   Dep(15) =  depmap( IXADV_CH3O2H,CDEP_OP  , -1.)
   Dep(16) =  depmap( IXADV_C2H5OOH,CDEP_OP  , -1.)
   Dep(17) =  depmap( IXADV_pNO3,  CDEP_COA, -1.)
   Dep(18) =  depmap( IXADV_PM25,  CDEP_FIN, -1. )
   Dep(19) =  depmap( IXADV_PMco,  CDEP_COA, -1. )
   Dep(20) =  depmap( IXADV_SSfi,  CDEP_FIN, -1. )
   Dep(21) =  depmap( IXADV_SSco,  CDEP_COA, -1. )
   Dep(22) =  depmap( IXADV_Pb210,  CDEP_FIN, -1. )

!#######################  NEW define indices here #######################

DDEP_SOX   = find_index("DDEP_SOX",f_2d(:)%name)
DDEP_OXN   = find_index("DDEP_OXN",f_2d(:)%name)
DDEP_RDN   = find_index("DDEP_RDN",f_2d(:)%name)

DDEP_OXSCF = find_index("DDEP_OXSCF",f_2d(:)%name)
DDEP_OXSDF = find_index("DDEP_OXSDF",f_2d(:)%name)
DDEP_OXSCR = find_index("DDEP_OXSCR",f_2d(:)%name)
DDEP_OXSSN = find_index("DDEP_OXSSN",f_2d(:)%name)

DDEP_OXNCF = find_index("DDEP_OXNCF",f_2d(:)%name)
DDEP_OXNDF = find_index("DDEP_OXNDF",f_2d(:)%name)
DDEP_OXNCR = find_index("DDEP_OXNCR",f_2d(:)%name)
DDEP_OXNSN = find_index("DDEP_OXNSN",f_2d(:)%name)

DDEP_RDNCF = find_index("DDEP_RDNCF",f_2d(:)%name)
DDEP_RDNDF = find_index("DDEP_RDNDF",f_2d(:)%name)
DDEP_RDNCR = find_index("DDEP_RDNCR",f_2d(:)%name)
DDEP_RDNSN = find_index("DDEP_RDNSN",f_2d(:)%name)

D2_AFSTDF0 = find_index("D2_AFSTDF0",f_2d(:)%name)
D2_AFSTDF16 = find_index("D2_AFSTDF16",f_2d(:)%name)

D2_AFSTBF0 = find_index("D2_AFSTBF0",f_2d(:)%name)
D2_AFSTBF16 = find_index("D2_AFSTBF16",f_2d(:)%name)

D2_AFSTCR0 = find_index("D2_AFSTCR0",f_2d(:)%name)
D2_AFSTCR3 = find_index("D2_AFSTCR3",f_2d(:)%name)
D2_AFSTCR6 = find_index("D2_AFSTCR6",f_2d(:)%name)

D2_O3DF    = find_index("D2_O3DF   ",f_2d(:)%name)
D2_O3WH    = find_index("D2_O3WH   ",f_2d(:)%name)

D2_EUAOT30WH    = find_index("D2_EUAOT30WH",f_2d(:)%name)
D2_EUAOT40WH    = find_index("D2_EUAOT40WH",f_2d(:)%name)
D2_EUAOT30DF    = find_index("D2_EUAOT30DF",f_2d(:)%name)
D2_EUAOT40DF    = find_index("D2_EUAOT40DF",f_2d(:)%name)

D2_UNAOT30WH    = find_index("D2_UNAOT30WH",f_2d(:)%name)
D2_UNAOT40WH    = find_index("D2_UNAOT40WH",f_2d(:)%name)
D2_UNAOT30DF    = find_index("D2_UNAOT30DF",f_2d(:)%name)
D2_UNAOT40DF    = find_index("D2_UNAOT40DF",f_2d(:)%name)

D2_MMAOT30WH    = find_index("D2_MMAOT30WH",f_2d(:)%name)
D2_MMAOT40WH    = find_index("D2_MMAOT40WH",f_2d(:)%name)

iam_wheat   = find_index("IAM_CR",LandDefs(:)%code)
iam_beech   = find_index("IAM_DF",LandDefs(:)%code)
iam_medoak  = find_index("IAM_MF",LandDefs(:)%code)
!####################### ds END of define indices #######################

  end subroutine Init_DepMap

  !<==========================================================================
  subroutine Add_ddep(debug_flag,dt,i,j,convfac,lossfrac,fluxfrac,c_hvegppb)

  !<==========================================================================
     ! Adds deposition losses to ddep arrays
     logical, intent(in) :: debug_flag
     real,    intent(in) :: dt              ! time-step
     integer, intent(in) :: i,j             ! coordinates
     real,    intent(in) ::  convfac, lossfrac
     real, dimension(:,:), intent(in) ::  fluxfrac   ! dim (NADV, NLANDUSE)
     real, dimension(:), intent(in) ::  c_hvegppb   ! dim (NLANDUSE)
     integer :: n, nadv, ihh, idd, imm
     real :: o3WH, o3DF   ! O3 over wheat, decid forest
     logical, parameter :: DEBUG_ECO = .false.

     integer, parameter :: N_OXS = 2        ! Number in ox. sulphur family
     real, parameter, dimension(N_OXS) :: OXS = &
             (/ IXADV_SO2, IXADV_SO4 /)
     integer, parameter :: N_OXN = 5        ! Number in ox. nitrogen family
     real, parameter, dimension(N_OXN) :: OXN = &
             (/ IXADV_HNO3, IXADV_PAN, IXADV_NO2, IXADV_aNO3, IXADV_pNO3 /)
     integer, parameter :: N_RDN = 2        ! Number in red. nitrogen family
     real, parameter, dimension(N_RDN) :: RDN = &
             (/ IXADV_NH3, IXADV_aNH4 /)

     real, parameter  :: NMOLE_M3 = 1.0e6*1.0e9/AVOG  ! Converts from 
                                                      ! mol/cm3 to nmole/m3

     real ::  to_nmole, timefrac, fstfrac 
     to_nmole =  NMOLE_M3
     timefrac = dt/3600.0
     fstfrac  = dt*1.0e-6     ! Converts also nmole to mmole


! waters and wetland categories removed after discussions with CCE/IIASA

!! OXIDIZED SULPHUR
!!-----------------------

     d_2d(DDEP_SOX,i,j,IOU_INST) = (  &
          DepLoss(IXADV_SO2) + DepLoss(IXADV_SO4) ) * convfac * atwS


     d_2d(DDEP_OXSCF,i,j,IOU_INST) = 0.0
     d_2d(DDEP_OXSDF,i,j,IOU_INST) = 0.0
     d_2d(DDEP_OXSCR,i,j,IOU_INST) = 0.0
     d_2d(DDEP_OXSSN,i,j,IOU_INST) = 0.0


     do n = 1,N_OXS
       nadv = OXS(n)

      ! ==    make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !


       d_2d(DDEP_OXSCF,i,j,IOU_INST) = d_2d(DDEP_OXSCF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_conif ) * DepLoss(nadv)


       d_2d(DDEP_OXSDF,i,j,IOU_INST) = d_2d(DDEP_OXSDF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_decid ) * DepLoss(nadv)

       d_2d(DDEP_OXSCR,i,j,IOU_INST) = d_2d(DDEP_OXSCR,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_crop ) * DepLoss(nadv)


       d_2d(DDEP_OXSSN,i,j,IOU_INST) = d_2d(DDEP_OXSSN,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_seminat ) * DepLoss(nadv)

     end do

     d_2d(DDEP_OXSCF,i,j,IOU_INST) = d_2d(DDEP_OXSCF,i,j,IOU_INST)*convfac*atwS
     d_2d(DDEP_OXSDF,i,j,IOU_INST) = d_2d(DDEP_OXSDF,i,j,IOU_INST)*convfac*atwS
     d_2d(DDEP_OXSCR,i,j,IOU_INST) = d_2d(DDEP_OXSCR,i,j,IOU_INST)*convfac*atwS
     d_2d(DDEP_OXSSN,i,j,IOU_INST) = d_2d(DDEP_OXSSN,i,j,IOU_INST)*convfac*atwS



!! OXIDIZED NITROGEN
!!-----------------------


     d_2d(DDEP_OXN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_HNO3) +  DepLoss(IXADV_PAN) +  DepLoss(IXADV_NO2) + &
          DepLoss(IXADV_aNO3)+  DepLoss(IXADV_pNO3)  ) * convfac * atwN


     d_2d(DDEP_OXNCF,i,j,IOU_INST) = 0.0
     d_2d(DDEP_OXNDF,i,j,IOU_INST) = 0.0
     d_2d(DDEP_OXNCR,i,j,IOU_INST) = 0.0
     d_2d(DDEP_OXNSN,i,j,IOU_INST) = 0.0

     do n = 1, N_OXN
       nadv = OXN(n)
         
      ! ==    make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !


       d_2d(DDEP_OXNCF,i,j,IOU_INST) = d_2d(DDEP_OXNCF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_conif ) * DepLoss(nadv)


       d_2d(DDEP_OXNDF,i,j,IOU_INST) = d_2d(DDEP_OXNDF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_decid ) * DepLoss(nadv)

       d_2d(DDEP_OXNCR,i,j,IOU_INST) = d_2d(DDEP_OXNCR,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_crop ) * DepLoss(nadv)


       d_2d(DDEP_OXNSN,i,j,IOU_INST) = d_2d(DDEP_OXNSN,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_seminat ) * DepLoss(nadv)

      ! ==                                                            ==== !


     end do

     d_2d(DDEP_OXNCF,i,j,IOU_INST) = d_2d(DDEP_OXNCF,i,j,IOU_INST)*convfac*atwN
     d_2d(DDEP_OXNDF,i,j,IOU_INST) = d_2d(DDEP_OXNDF,i,j,IOU_INST)*convfac*atwN
     d_2d(DDEP_OXNCR,i,j,IOU_INST) = d_2d(DDEP_OXNCR,i,j,IOU_INST)*convfac*atwN
     d_2d(DDEP_OXNSN,i,j,IOU_INST) = d_2d(DDEP_OXNSN,i,j,IOU_INST)*convfac*atwN



!! REDUCED NITROGEN
!!-----------------------

     d_2d(DDEP_RDN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_NH3) +  DepLoss(IXADV_aNH4)  ) * convfac * atwN



     d_2d(DDEP_RDNCF,i,j,IOU_INST) = 0.0
     d_2d(DDEP_RDNDF,i,j,IOU_INST) = 0.0
     d_2d(DDEP_RDNCR,i,j,IOU_INST) = 0.0
     d_2d(DDEP_RDNSN,i,j,IOU_INST) = 0.0

     do n = 1, N_RDN
       nadv = RDN(n)
         
      ! ==    make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !


       d_2d(DDEP_RDNCF,i,j,IOU_INST) = d_2d(DDEP_RDNCF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_conif ) * DepLoss(nadv)


       d_2d(DDEP_RDNDF,i,j,IOU_INST) = d_2d(DDEP_RDNDF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_decid ) * DepLoss(nadv)

       d_2d(DDEP_RDNCR,i,j,IOU_INST) = d_2d(DDEP_RDNCR,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_crop ) * DepLoss(nadv)


       d_2d(DDEP_RDNSN,i,j,IOU_INST) = d_2d(DDEP_RDNSN,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,:), LandType(:)%is_seminat ) * DepLoss(nadv)

      ! ==                                                            ==== !

     end do

     d_2d(DDEP_RDNCF,i,j,IOU_INST) = d_2d(DDEP_RDNCF,i,j,IOU_INST)*convfac*atwN
     d_2d(DDEP_RDNDF,i,j,IOU_INST) = d_2d(DDEP_RDNDF,i,j,IOU_INST)*convfac*atwN
     d_2d(DDEP_RDNCR,i,j,IOU_INST) = d_2d(DDEP_RDNCR,i,j,IOU_INST)*convfac*atwN
     d_2d(DDEP_RDNSN,i,j,IOU_INST) = d_2d(DDEP_RDNSN,i,j,IOU_INST)*convfac*atwN

!MAPPING_MANUAL CHANGES:
!  Use 1.6 for Beech and 3 for crops

!Beech:
     d_2d(D2_AFSTDF0,i,j,IOU_INST) =  fstfrac*leaf_flux(iam_beech)
     d_2d(D2_AFSTDF16,i,j,IOU_INST) = fstfrac* max(leaf_flux(iam_beech)-1.6,0.0)
!Med. Oak:
     d_2d(D2_AFSTBF0,i,j,IOU_INST) =  fstfrac*leaf_flux(iam_medoak)
     d_2d(D2_AFSTBF16,i,j,IOU_INST) = fstfrac* max(leaf_flux(iam_medoak)-1.6,0.0)
!Crops
     d_2d(D2_AFSTCR0,i,j,IOU_INST) =  fstfrac*leaf_flux(iam_wheat)
     d_2d(D2_AFSTCR3,i,j,IOU_INST) =  fstfrac*max(leaf_flux(iam_wheat)-3.0,0.0)
     d_2d(D2_AFSTCR6,i,j,IOU_INST) =  fstfrac*max(leaf_flux(iam_wheat)-6.0,0.0)

   !--- ecosystem specific concentrations..
   ! - use Conif forest for forests - safer for growing seasons

     imm      =    current_date%month            ! for debugging
     idd      =    current_date%day              ! for debugging
     ihh      =    current_date%hour             ! for debugging


     o3WH = c_hvegppb(iam_wheat)* lossfrac
     o3DF = c_hvegppb(iam_beech)* lossfrac

     d_2d(D2_O3DF,i,j,IOU_INST) =   o3DF
     d_2d(D2_O3WH,i,j,IOU_INST) =   o3WH

     if ( ihh >= 9 .and. ihh <= 21 ) then ! 8-20 CET, assuming summertime

        d_2d(D2_EUAOT30WH,i,j,IOU_INST) =  max(o3WH-30.0,0.0) * timefrac
        d_2d(D2_EUAOT40WH,i,j,IOU_INST) =  max(o3WH-40.0,0.0) * timefrac
        d_2d(D2_EUAOT30DF,i,j,IOU_INST) =  max(o3DF-30.0,0.0) * timefrac
        d_2d(D2_EUAOT40DF,i,j,IOU_INST) =  max(o3DF-40.0,0.0) * timefrac
     else
        d_2d(D2_EUAOT30WH,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT40WH,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT30DF,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT40DF,i,j,IOU_INST) =  0.0
     end if


    !/-- Calcuates AOT values for specific veg. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )


           if ( Grid%izen < AOT_HORIZON ) then

             d_2d(D2_UNAOT30WH,i,j,IOU_INST) =  max(o3WH-30.0,0.0) * timefrac
             d_2d(D2_UNAOT40WH,i,j,IOU_INST) =  max(o3WH-40.0,0.0) * timefrac
             d_2d(D2_UNAOT30DF,i,j,IOU_INST) =  max(o3DF-30.0,0.0) * timefrac
             d_2d(D2_UNAOT40DF,i,j,IOU_INST) =  max(o3DF-40.0,0.0) * timefrac

           else
             d_2d(D2_UNAOT30WH,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT40WH,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT30DF,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT40DF,i,j,IOU_INST) =  0.0
           end if

       ! MM AOT added (same as UNECE, but different growing season)
             d_2d(D2_MMAOT30WH,i,j,IOU_INST) = d_2d(D2_UNAOT30WH,i,j,IOU_INST)&
                     * WheatGrowingSeason(i,j)
             d_2d(D2_MMAOT40WH,i,j,IOU_INST) = d_2d(D2_UNAOT40WH,i,j,IOU_INST)&
                     * WheatGrowingSeason(i,j)

    if ( DEBUG_ECO .and. debug_flag ) then
          write(6,"(a12,3i5,f7.2,5es12.3,i3,es12.3)") "DEBUG_ECO ", &
          imm, idd, ihh, o3WH, &
             leaf_flux(iam_beech), d_2d(D2_AFSTDF0,i,j,IOU_INST), &
             leaf_flux(iam_wheat), d_2d(D2_AFSTCR0,i,j,IOU_INST), &
             d_2d(D2_UNAOT40WH,i,j,IOU_INST), &
             WheatGrowingSeason(i,j), d_2d(D2_MMAOT40WH,i,j,IOU_INST)
    end if ! DEBUG

   !---- end ecosystem specific ----------------------------------------------

  end subroutine  Add_ddep

end module My_DryDep_ml

