! <My_BoundConditions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!################ OZONE model ###############################################
module My_BoundConditions_ml
!____________________________________________________________________________
! This module specifies model-dependant:
!  (A) boundary-condition (bc) indices and 
!  (B) mapping arrays "bc2xn_adv" and "bc2xn_bgn" 
!       - which tell how  boundary conditions (bcs) are to be assigned to emep 
!         concentrations (xn_adv, xn_bgn). 
!
! THIS FILE WILL CHANGE FOR DIFFERENT CHEMISTRIES - MUST BE SUPPLIED BY USER. 
! AS A FIRST INDICATION OF THIS I HAVE SUPPLIED A "MY_MODEL" LABEL BELOW.
!
! The module BoundaryConditions_ml calls up this module with just:
!
!  call My_bcmap() 
!
! So far this module copes only with boundary condistion supplied either
! by the UiO global model (through the UiO_ml), or defined here as
! constant mixing ratios.
!
! Language: F
! History :
! ds - December 2000-January 2001
! hf - september 2001 Misc BC's added as a function of sigma
!____________________________________________________________________________
! IMPORTANT NOTE:
! The routines given here are constructed around the global model fields from 
! the University of Oslo (T21) global model. In order to use other models as 
! bcs  then usually these routines will have to be replaced by model-specific 
! routines. The important thing is that the inputs and outputs from the routine
! are independant of the global module used.
!_____________________________________________________________________________
  use GenSpec_bgn_ml, only: NSPEC_BGN
  use GenSpec_adv_ml, only: NSPEC_ADV & 
                ,IXADV_H2,IXADV_O3, IXADV_SO2, IXADV_SO4  &
                ,IXADV_HNO3,IXADV_PAN,IXADV_CO,IXADV_C2H4  &
                ,IXADV_C2H6,IXADV_NC4H10,IXADV_HCHO,IXADV_CH3CHO   &
                ,IXADV_H2O2,IXADV_CH3O2H,IXADV_ISOP,IXADV_NO,IXADV_NO2  &
                ,IXADV_CH4,IXADV_aNH4,IXADV_pNO3,IXADV_aNO3 &
                ,IXADV_CH3COO2
  use GenSpec_shl_ml, only: IXSHL_OH
  use GridValues_ml,  only: sigma_mid    !sigma layer midpoint
  use Met_ml         ,only : z_mid       ! height of half layers
  use ModelConstants_ml , only: KMAX_MID, NPROC ! No. levels in vertical, processors
  use Par_ml,         only: me
  use GlobalBCs_ml,  only:  NGLOB_BC  &!  indices from UiO model
                ,IBC_SO2, IBC_SO4, IBC_HCHO, IBC_CH3CHO &
                ,IBC_O3,IBC_HNO3,IBC_PAN,IBC_CO,IBC_C2H6   &
                ,IBC_C4H10, IBC_NO ,IBC_NO2,IBC_aNH4,IBC_aNO3,IBC_pNO3&
               ,IBC_H2O2,IBC_CH3COO2
  implicit none
  private

 !/-- subroutines
 public :: My_bcmap          ! sets bc2xn_adv, bc2xn_bc, and  misc_bc

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
 !/-- model-type
 !    for consistency checks, possibly to match label in My_model_ml??
  character(len=12), public :: MY_MODEL = "emepds"

   logical, public, parameter  :: BGN_2D = .false. !No 2d bgn species
   logical, private, parameter :: DEBUG_MYBC = .false.


 ! A. Set indices
 ! ===========================================================================
 ! For species which have constant mixing ratios:

 integer, public, parameter ::  NMISC_BC  =  2      ! H2, CH4   ! OC
 
 integer, public, parameter :: IBC_H2  = NGLOB_BC + 1   &
                              ,IBC_CH4 = NGLOB_BC + 2
                              !6s ,IBC_OC  = NGLOB_BC + 2

 integer, public, parameter ::  NTOT_BC  = NGLOB_BC + NMISC_BC

 ! We also need the array misc_bc to specify concentrations of these species:

 real, public, save, dimension(NGLOB_BC+1:NTOT_BC,KMAX_MID) :: misc_bc 
!real, public, save, dimension(NGLOB_BC+1:NTOT_BC) :: misc_bc 

 ! B. Define mapping arrays
 ! ===========================================================================
 ! The mapping is done through the arrays bc2xn_adv and bc2xn_bgn, such that 
 ! the emep species are given along the x-dimension and the bc species along 
 ! the y.  e.g., the statement
 !
 !     bc2xn_adv(IBC_NOX,IXADV_NO2) = 0.55 
 !
 ! would assign the BC concentration of NOX  to the EMEP model concentration
 ! of NO2 after multiplication with a factor 0.55.
 !(The CTM2 concentration of NOx used as BC after a multilication 
 ! with a factor 0.55.)
 !_______________________________________________________________________

    real, public, save, dimension(NTOT_BC,NSPEC_ADV) :: bc2xn_adv  ! see above
    real, public, save, dimension(NTOT_BC,NSPEC_BGN) :: bc2xn_bgn ! see above

 !-------
 contains
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine My_bcmap(iyr_trend)    ! sets bc2xn_adv, bc2xn_bc, and  misc_bc
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    integer, intent(in) :: iyr_trend !ds Year for which BCs are wanted 
    real    :: ppb = 1.0e-9
    real :: trend_ch4  !ds rv1.6.11
    integer :: ii,i,j,k 
    real :: decrease_factor(NGLOB_BC+1:NTOT_BC) ! Decrease factor for misc bc's
                                      ! Gives the factor for how much of 
                                      ! the top-layer conc. that is left 
                                      ! at bottom layer

    real :: top_misc_bc(NGLOB_BC+1:NTOT_BC) ! Conc. at top of misc bc
!    real :: ratio_length(KMAX_MID)    ! Vertical length of the actual layer
                                      ! divided by length from midpoint of 
                                      ! layer 1 to layer KMAX_MID

    misc_bc            = 0.0  ! Initialise
    bc2xn_adv          = 0.0  ! Initialise
    bc2xn_bgn          = 0.0  ! Initialise

    ! Own (constant mixing ratio) boundary conditions **********

    ! NOTE - these species have to have the bc2xn_ indices set to 1.0 for either
    ! the advected or the background concentrations, in order that the 
    ! concentrations specified in misc_bc are transferred correctly into the 
    ! boundary conditions.
    !
    ! 18.09.01 -hf- misc bc'c as function of sigma possible 

        !ds top_misc_bc(IBC_CH4) = 1760.0 * ppb

        !ds set values of 1625 in 1980, 1780 in 1990, and 1820 in 2000. Interpolate
        ! between these for other years. Values from EMEP Rep 3/97, Table 6.2 for
        ! 1980, 1990, and from CDIAC (Mace Head) data for 2000.
 
        if ( iyr_trend >= 1990 ) then

            top_misc_bc(IBC_CH4) = 1780.0 + &
                        (iyr_trend-1990) * 0.1*(1820-1780.0)

        else 

            top_misc_bc(IBC_CH4) = 1780.0 * &
                        exp(-0.01*0.91*(1990-iyr_trend)) ! Zander,1975-1990
                     !exp(-0.01*0.6633*(1975-iyr_trend)) !Zander,1951-1975, check
        end if
        trend_ch4 = top_misc_bc(IBC_CH4)/1780.0  ! Crude for now.
        if (me== 0) write(6,"(a20,i5,2f12.3)") "TREND CH4", iyr_trend, trend_ch4, top_misc_bc(IBC_CH4)

        top_misc_bc(IBC_CH4)  =  top_misc_bc(IBC_CH4) * ppb

        top_misc_bc(IBC_H2)  =  600.0 * ppb

        !! top_misc_bc(IBC_OC) =   0.0 !!! 1.0 * ppb

        decrease_factor(IBC_H2) =0.0 !No increase/decrease with height
        decrease_factor(IBC_CH4)=0.0 !No increase/decrease with height


!a) Function of height(not included yet!):
!ratio_length(i,j,k)=(z_mid(i,j,1)-z_mid(i,j,k))/ &
!                      (z_mid(i,j,1)-z_mid(i,j,KMAX_MID))
!Replace sigma_mid with ratio_length and make misc_bc 4 dimentional
!
!b)Function of sigma_mid (I assume that top_misc_bc is given for
!top(sigma_bnd(1)=0.), and that at ground (sigma_bnd(KMAX_BND)=1.) the conc.
!is top_misc_bc -decrease_factor*top_misc_bc. Since I choose to set the 
! concentration as a factor of sigma_mid, the concentration in the lowest 
! grid cell will not be 
! excactly  top_misc_bc -decrease_factor*top_misc_bc, but close.

        do ii=NGLOB_BC+1,NTOT_BC
           do k=1,KMAX_MID
              misc_bc(ii,k) = top_misc_bc(ii) - &
                  top_misc_bc(ii)*decrease_factor(ii)*sigma_mid(k) 
              if (me == 0) then
                 if (DEBUG_MYBC) write(*,"(a20,2es12.4,i4)")"height,misc_vert,k", &
                                               sigma_mid(k),misc_bc(ii,k),k
              endif
           enddo
        enddo

!        misc_bc(IBC_H2)     = 600.0 * ppb
!        misc_bc(IBC_OC)     =   1.0*ppb !!! 0.0

        bc2xn_adv(IBC_H2,  IXADV_H2)    = 1.0
        bc2xn_adv(IBC_CH4, IXADV_CH4)   = 1.0

        !/-- check, just in case we forgot something...!

        if ( DEBUG_MYBC  ) then
             print *, "In My_bcmap, NGLOB_BC, NTOT_BC is", NGLOB_BC, NTOT_BC
             do i = NGLOB_BC+1 , NTOT_BC
                print *, "In My_bcmap, sum-adv", i, " is", sum(bc2xn_adv(i,:))
                print *, "In My_bcmap, sum-bgn", i, " is", sum(bc2xn_bgn(i,:))
             end do
        end if ! DEBUG

        do i = NGLOB_BC+1 , NTOT_BC
           if ( sum(bc2xn_adv(i,:)) + sum(bc2xn_bgn(i,:)) /= 1.0 )then 
              WRITE(*,*) 'MPI_ABORT: ', "BCproblem - my" 
              call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
           endif
        end do


 ! mappings for species from LOgan + obs model ***********

  bc2xn_adv(IBC_O3      ,IXADV_O3      )   =   1.0
  bc2xn_adv(IBC_HNO3    ,IXADV_HNO3    )   =   1.0
  bc2xn_adv(IBC_SO2     ,IXADV_SO2     )   =   1.0
  bc2xn_adv(IBC_SO4     ,IXADV_SO4    )   =   1.0
  bc2xn_adv(IBC_PAN     ,IXADV_PAN     )   =   1.0
  bc2xn_adv(IBC_CO      ,IXADV_CO      )   =   1.0
  bc2xn_adv(IBC_C2H6    ,IXADV_C2H6    )   =   1.0
  bc2xn_adv(IBC_C4H10   ,IXADV_NC4H10  )   =   1.0
  bc2xn_adv(IBC_NO      ,IXADV_NO      )   =   1.0
  bc2xn_adv(IBC_NO2     ,IXADV_NO2     )   =   1.0
  bc2xn_adv(IBC_HCHO    ,IXADV_HCHO    )   =   1.0
  bc2xn_adv(IBC_CH3CHO  ,IXADV_CH3CHO  )   =   1.0
  bc2xn_adv(IBC_aNO3    ,IXADV_aNO3    )   =   1.0
  bc2xn_adv(IBC_pNO3    ,IXADV_pNO3    )   =   1.0
  bc2xn_adv(IBC_aNH4    ,IXADV_aNH4    )   =   1.0
!hfOH NEW: When a smaller domain than the full, there are BCs for these from OZONE  
!hf err  bc2xn_bgn(IBC_H2O2    ,IXADV_H2O2    )   =   1.0!hfOH
!hf err   bc2xn_bgn(IBC_CH3COO2 ,IXADV_CH3COO2 )   =   1.0 !hfOH
!hf err   bc2xn_bgn(IBC_OH      ,IXSHL_OH    )   =   1.0 !hfOH 
  bc2xn_adv(IBC_H2O2    ,IXADV_H2O2    )   =   1.0!hfOH
  bc2xn_adv(IBC_CH3COO2 ,IXADV_CH3COO2 )   =   1.0 !hfOH


! The following species are excluded either because they have no corresponding
! species in the emep model, or because they have lifetimes which are so
! short that initialisation is uncessary.
!-----------------------------------------------------------------------------
!u3  bc2xn_adv(IBC_C2H4    ,IXADV_C2H4  )   =   1.0
!u3  bc2xn_adv(IBC_C3H6    ,IXADV_NC4H10  )   =   0.75
!u3  bc2xn_adv(IBC_C6H14   ,IXADV_NC4H10  )   =   1.5   ! ds - scale by C6/C4
!u3  bc2xn_adv(IBC_CH2O    ,IXADV_HCHO    )   =   1.0   ! ds -rename
!u3  bc2xn_adv(IBC_CH3CHO  ,IXADV_CH3CHO  )   =   1.0
!u3  bc2xn_adv(IBC_H2O2    ,IXADV_H2O2    )   =   1.0
!u3  bc2xn_adv(IBC_CH3O2H  ,IXADV_CH3O2H  )   =   1.0
!u3  bc2xn_adv(IBC_ISOPRENE,IXADV_ISOP    )   =   1.0   ! ds-rename
!u3  bc2xn_adv(IBC_RCOHCO  ,IXADV_CH3CHO  )   =   1.0   ! ds - unknown scale
!u3  bc2xn_adv(IBC_CH4     ,IXADV_CH4     )   =   1.0   ! Re-included for DSMACH
!u3  bc2xn_adv(IBC_C3H8    ,IXADV_NC4H10 )   =   0.75   ! mini
!u3  bc2xn_adv(IBC_C3H8    ,IXADV_NC4H10  )   =   0.5   ! ds-split C2H6/NC4H10
!!bc2xn_adv(IBC_NOX ,IXADV_NOX     )   =  -1.0   ! Excluded, we have NO and NO2
!!bc2xn_adv(IBC_C6HXR   ,IXADV_C6HXR   )   =   1.0
!!bc2xn_adv(IBC_HO2NO2  ,IXADV_HO2NO2  )   =   1.0
!!bc2xn_adv(IBC_CH3COY  ,IXADV_CH3COY  )   =   1.0
!!bc2xn_adv(IBC_CH3COX  ,IXADV_CH3COX  )   =   1.0
!!bc2xn_adv(IBC_HO2     ,IXADV_HO2     )   =  -1.0   ! SHort-lived
!!bc2xn_adv(IBC_CH2O2OH ,IXADV_CH2O2OH )   =   1.0
!!bc2xn_adv(IBC_CH3COB  ,IXADV_CH3COB  )   =  -1.0   ! ???
!!bc2xn_adv(IBC_CH3XX   ,IXADV_CH3XX   )   =   1.0
!!bc2xn_adv(IBC_AR1     ,IXADV_AR1     )   =  -1.0   ! ???
!!bc2xn_adv(IBC_AR2     ,IXADV_AR2     )   =   1.0
!!bc2xn_adv(IBC_AR3     ,IXADV_AR3     )   =  -1.0   ! ???
!!bc2xn_adv(IBC_ISOR1   ,IXADV_ISOR1   )   =  -1.0   ! Re-Excluded
!!bc2xn_adv(IBC_ISOK    ,IXADV_ISOK    )   =  -1.0   ! ??
!!bc2xn_adv(IBC_ISOR2   ,IXADV_ISOR2   )   =  -1.0   ! SHort-lived
!!bc2xn_adv(IBC_HCOHCO  ,IXADV_HCOHCO  )   =  -1.0   ! Excluded
!!bc2xn_adv(IBC_CH3X    ,IXADV_CH3X    )   =  -1.0   ! ??
!!bc2xn_adv(IBC_NO3     ,IXADV_NO3     )   =  -1.0   ! SHort-lived
!!bc2xn_adv(IBC_N2O5    ,IXADV_N2O5    )   =   1.0   ! jej - 000927
!!bc2xn_adv(IBC_C3H7O2  ,IXADV_C3H7O2  )   =   1.0
!!bc2xn_adv(IBC_ACETONE ,IXADV_ACETON  )   =   1.0
!!bc2xn_adv(IBC_CH3COD  ,IXADV_CH3COD  )   =  -1.0   ! ??
!!bc2xn_adv(IBC_NOZ     ,IXADV_NOZ     )   =  -1.0   ! ???
!!bc2xn_adv(IBC_CH3O2   ,IXADV_CH3O2   )   =   Short-lived
!!bc2xn_adv(IBC_C2H5O2  ,IXADV_C2H5O2  )   =   Short-lived
!!bc2xn_adv(IBC_C4H9O2  ,IXADV_C4H9O2  )   =   Short-lived
!!bc2xn_adv(IBC_C6H13O2 ,IXADV_C6H13O2 )   =   Short-lived
!!bc2xn_adv(IBC_O3P     ,IXADV_O3P     )   =  -1.0   ! Short-lived
!!bc2xn_adv(IBC_O1D     ,IXADV_O1D     )   =  -1.0   ! Short-lived
!!bc2xn_adv(IBC_OH      ,IXADV_OH      )   =  -1.0   ! Short-lived
!!bc2xn_adv(IBC_O3NO    ,IXADV_O3NO    )   =  -1.0   ! Excluded
!!bc2xn_adv(IBC_DMS     ,IXADV_DMS     )   =  -1.0   ! Query ????
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 end subroutine My_bcmap
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module My_BoundConditions_ml
!_____________________________________________________________________________




