! <My_SOA_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
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
module OrganicAerosol_ml

  ! Calculates the amount of condensible species in the gas and aerosol phases. 
  !
  ! References:
  !   S2007: Simpson, D. et al., JGR, 2007, 
  !   B2012: Bergström, R. et al., ACP (in press Sept.), 2012, 
  !   S2012: Simpson, D. et al., ACP, 2012, 
  !
  !
  ! Usage: call OrganicAerosol from Runchem, after setup of column data
  !  (for meteorology, etc.). The subroutine initialises itself on the 
  !  first call and thereafter modifies two external variables:
  !   xn(k,SOA) : the concentrations of SOA   (For a 1-d column array(in EMEP))
  !   Fgas(k,X) : The fraction of X which is gas and not aeorosol
  !
  ! ----------------------------------------------------------------------------
  !
  ! From Gas/Particle theory, A/G = K.COA,
  ! therefore, Fgas = G/(G+A) = 1/(1+K.COA)
  !
  !-----------------------------------------------------------------------------
  ! NB- we exclude use of gamma for now, but leave commented out code
  !-----------------------------------------------------------------------------
  !
  ! Dave Simpson, August 2001 -- 2012
  ! Robert Bergström     2010 -- 2012
  ! 
  !--------------------------------------------------------------------------

  ! Functions + GridValues + PT only for BGNDOC
   use Functions_ml, only: StandardAtmos_kPa_2_km !ds for use in Hz scaling
use CheckStop_ml, only : StopAll  ! CHECKXN
   use ChemFields_ml,      only : Fgas3d   !  stores 3-d  between time-steps
   use ChemChemicals_ml,      only : species   ! for molwts
   use ChemSpecs_tot_ml,  S1 => FIRST_SEMIVOL , S2 => LAST_SEMIVOL

   use ChemGroups_ml, only :    &
    NONVOLPCM_GROUP &
    ,NVABSOM_GROUP & ! nonvolatile absorbing species for OA partitioning
    ,ASOA => ASOA_GROUP &
    ,BSOA => BSOA_GROUP &
    ,ECFINE => ECFINE_GROUP &                    
    ,SVFFUELOA25  => SVFFUELOA25_GROUP & ! semi-volatile FFUELOA25 
      ,SVWOODOA25  => SVWOODOA25_GROUP &    ! for VBS semivolatile WOOD BURNING OA (primary + aged!)
      ,SVFFIREOA25  => SVFFIREOA25_GROUP &    ! for VBS semivolatile Wildfire OA (primary + aged!)
  ,FFUELEC     => FFUELEC_GROUP &
  ,NVFFUELOC25 => NVFFUELOC25_GROUP & ! non-vol. FFUELOC emis. (in PM2.5)
  ,NVFFUELOCCO => NVFFUELOC_COARSE_GROUP & ! non-vol. " " , (PM2.5-10 frac.)
  ,NVWOODOC25  => NVWOODOC25_GROUP & ! non-vol. WOODOC emissions 
                                     !  (zero in VBS-PAX type runs)
  ,NVFFIREOC25 => NVFFIREOC25_GROUP  ! only non-vol. FFIREOC emissions
                                     ! (zero in VBS-PAX type runs)

   use GridValues_ml, only: A_mid,B_mid, debug_proc, debug_li, debug_lj
   use ModelConstants_ml,    only :  PT

   use ModelConstants_ml,    only : CHEMTMIN, CHEMTMAX, &
                                    MasterProc, DEBUG => DEBUG_SOA, &
                                    K2 => KMAX_MID, K1 => KCHEMTOP
   use Par_ml,               only : LIDIM => MAXLIMAX, LJDIM => MAXLJMAX, me
   use PhysicalConstants_ml, only : AVOG, RGAS_J 
   use Setup_1dfields_ml,    only : itemp, xn => xn_2d, Fgas, Fpart
   use TimeDate_ml,       only: current_date
   implicit none
   private


   !/-- subroutines

    public   :: Init_OrganicAerosol
    public   :: OrganicAerosol


   !/-- public

    logical, public, parameter :: ORGANIC_AEROSOLS = .true.

   ! We store some values in 3-D fields, to allow the next G/P partitioning
   ! calculation to  start off with values of COA, mw and Fgas which 
   ! are about right. Ensures that very few iterations are needed.

!  real,public, save, dimension(S1:S2,LIDIM,LJDIM,K1:K2) :: &
!            Grid_SOA_Fgas           !EXC Grid_SOA_gamma

  real,public, save, allocatable, dimension(:,:,:) :: Grid_COA

  real, private, allocatable, dimension(:), save :: &
        COA           & ! Org. aerosol, ug/m3  
                        ! (this version does not include EC as absorber)
       ,BGND_OC       & ! FAKE FOR NOW, 0.50 ugC/m3 at surface
       ,BGND_OA         ! Assumed OA/OC=2, -> 1 ug/m3


  !VBS real, private, dimension(S1:S2,K1:K2), save :: &
  ! TMP - we assign Fpart for all species for now, since
  ! it makes it easier to code for nonvol and vol
! From Setup_1dfields now
!  real, private, dimension(1:NSPEC_TOT,K1:K2), save :: &
!              Fgas    & ! Fraction in gas-phase
!             ,Fpart!   & ! Fraction in gas-phase
             !VBS ,tabRTpL   ! = 1.0e-6 * R.T(i)/pL(Ti) for all temps i:
               !EXC ,gamma   & ! activity coefficient

  real, parameter, public :: SMALLFN  = 1.0e-20 ! Minimum value of ug allowed

   !/-- private

   ! ug = array for aerosol masses (ug/m3). Includes non-volatile compounds:
   ! TMP??? Excluding NVOL for now?

    real, private,allocatable, dimension(:,:), save :: ug_semivol 
! - use new NONVOLOC grpup to define:
    integer, private, parameter :: NUM_NONVOLPCM = size(NONVOLPCM_GROUP)
    integer, private, parameter :: NUM_NVABSOM = size(NVABSOM_GROUP)
!    integer, private, parameter, dimension(NUM_NONVOL) ::  &
!      NONVOL = (/ NONVOLOC_GROUP, NONVOLEC_GROUP /) ! OC+EC in partitioning OM
    real, private,allocatable, dimension(:,:), save :: ug_nonvol 
    !real, private, dimension(K1:K2), save :: ug_ecf ! CityZen added 

    real,  private, save, dimension(S1:S2,CHEMTMIN:CHEMTMAX) :: tabCiStar

    integer, private, save :: NITER = 2              ! No. iterations for Ksoa
    real, private, save :: xn2molem3  ! Conversion from molec/cm3 to mole/m3
    real, private, save :: xn2ugC, ugC2xn


   !/-- DEBUG variables
   !    Usually, DEBUG_SOA gives extra outputs. debug_flag is used to allow
   !    some extra outputs for a gven i,j - set in CTM model.

    character(len=20), public, save     :: soa_errmsg     = "ok"
    character(len=*), public, parameter :: SOA_MODULE_FLAG="VBS"

   contains
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine Init_OrganicAerosol(i,j,debug_flag)
   integer, intent(in) :: i,j
   logical, intent(in) :: debug_flag
   integer :: is,  it, k
   real, parameter :: kJ = 1000.0  
   real,allocatable, dimension(:), save :: p_kPa, h_km ! for standard atmosphere 
   logical, save :: my_first_call = .true.


   if ( my_first_call ) then
      allocate(COA(K1:K2),BGND_OC(K1:K2),BGND_OA(K1:K2))
      allocate(ug_semivol(S1:S2,K1:K2))
      allocate(Grid_COA(LIDIM,LJDIM,K1:K2))
      allocate(ug_nonvol(NUM_NVABSOM,K1:K2))
      allocate(p_kPa(K2), h_km(K2))

    !=========================================================================
    ! Set up background OM 
      ! Need to convert aeros to ug/m3 or ugC/m3.  Since xn is in molecules/cm3
      ! we divide by AVOG to get mole/cm3, multiply by 1e6 to get mole/m3,
      ! by  mol. weight to get g/m3 and by 1e6 to get ug/m3

         xn2molem3 = 1.0/AVOG * 1.0e6
         xn2ugC   = xn2molem3 * 12.0 * 1.0e6
         ugC2xn   = 1/xn2ugC

  ! Use Standard Atmosphere to get average heights of layers

       p_kPa(:) = 0.001*( A_mid(:)+B_mid(:)*101325.0 ) ! Pressure in kPa
       h_km     = StandardAtmos_kPa_2_km(p_kPa)
       BGND_OC(:)= 0.2 * 1.005 ! ng/m3 !!! will give 0.2 ugC/m3 at z=0 m

       do k = K1, K2
            BGND_OC(k) = BGND_OC(k) * exp( -h_km(k)/9.1 )
            if(DEBUG .and. MasterProc ) write(*,"(a,i4,2f8.3)") &
               "BGND_OC ", k, h_km(k),  BGND_OC(k)
       end do
       BGND_OA(:) = 2*BGND_OC(:)   ! Assume OA/OC = 2 for bgnd

       do k = K1, K2
          Grid_COA(:,:,k) = BGND_OA(k)  ! Use OA, not OC here
       end do

    !=========================================================================
    ! Set up Tables for Fcond 
       ! Ci = 1.0e6*P0/RT 
       ! Now, pi(T) = Ai exp(-Hi/RT)
       ! And pi(T) = Pi(Tref) * exp( H/RT * (1/Tref - 1/T) )
       ! ->  Ci(T) = Ci(Tref) * Tref/T * exp(...)

       do is=S1,S2
         do it=CHEMTMIN,CHEMTMAX

         ! C*-values are given for 298K according to most(?) publications.
           tabCiStar(is,it) = species(is)%CiStar * 298./it * &
                  exp( species(is)%DeltaH * kJ/RGAS_J * (1.0/298. - 1.0/it) )
         end do
       end do


         if ( MasterProc ) then 
            do is = S1, S2
               write(6,"(a,i4,a20,f7.1,i3,8es10.2)") &
                " Tab SOA: MW, Carbons, C*:", is, trim(species(is)%name), &
                 species(is)%molwt, species(is)%carbons, & 
                 tabCiStar(is,273), tabCiStar(is,303)
            end do
         end if


         allocate( Fgas3d(S1:S2,LIDIM,LJDIM,K1:K2) )

       !+ initial guess (1st time-step only)
       ! Fgas3D is only defined for the semivol stuff, so no need for nonvol here
       ! We need to assume something on 1st time-step though:
         Fgas3d = 1.0

       ! Initial values. Should not change except for semi-volatiles

        Fpart(:,:)         = 0.0
        Fpart(NONVOLPCM_GROUP,:)  = 1.0
        Fgas(:,:)         = max(0.0, 1.0 - Fpart(:,:) )

           !VBS Grid_avg_mw    = 250.0              ! Da
           !EXC Grid_SOA_gamma = 1.0
    !=========================================================================
!     if (debug_proc .and. i == debug_li .and. j == debug_lj ) then
!          print *, " CHECKXN xn(PART_FFIREOA25_OC,k) ",me,i,j, xn(PART_FFIREOA25_OC,K2)
!     end if
       
    !=========================================================================
        my_first_call = .false.
    end if ! my_first_call

    ! We need to set Fgas at start of each Runchem i,j loop, as it is
    ! used for rcemis:

      Fgas(S1:S2,:) = Fgas3d(S1:S2,i,j,:)     ! Semivolatiles only in 3D Fgas
      Fpart(S1:S2,:)  = 1-Fgas(S1:S2,:)

!      Fgas(NONVOLPCM_GROUP,:) = 0.0             !  not needed, shouldn't change

  end subroutine Init_OrganicAerosol

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine OrganicAerosol(i_pos,j_pos,debug_flag)

   integer, intent(in) :: i_pos, j_pos
   logical, intent(in) :: debug_flag 

   integer :: i,  k, iter, ispec   ! loop variables 
!   real, save :: molcc2ngm3 = 1.0e9*1.0e6/AVOG  !molecules/cc-> ng/m3
   real, save :: molcc2ugm3 = 1.0e12/AVOG  !molecules/cc-> ng/m3
!   real, save :: ngm32molcc = AVOG/1.0e15  !molecules/cc <- ng/m3
!   real, save :: ugm32molcc = AVOG/1.0e12  !molecules/cc <- ng/m3
   real :: Ksoa
   integer :: nmonth, nday, nhour, seconds

  ! Outputs:
!   real :: surfOC25, surfASOA, surfBSOA, surfFFUELOC25, surfWOODOC25, surfBGNDOC, surfOFFUELOA25_OC, surfFFUELOA25_OC, surfOWOODOA25, surfWOODOA25, surfOFFIREOA25, surfFFIREOA25
   real :: surfOC25, surfASOA, surfBSOA, surfBGNDOC, surfFFUELOA25_OC, surfWOODOA25_OC, surfFFIREOA25_OC


   nmonth = current_date%month
   nday   = current_date%day
   nhour  = current_date%hour
   seconds = current_date%seconds

   if( DEBUG .and. debug_flag) write(unit=*,fmt=*) "Into SOA"

! Note that xn(SOA) is not strictly needed with the method we have, but it
! is a good first guess of the sum of the condensed phases, enables easy output
! and saves the need for iteration. 
!
! Remember also that Fgas is saved, so will initially have been set by the
! preceding call to OrganicAerosol for a different set of i,j.

 ! 1st guesses:

  COA(:)          =  Grid_COA(i_pos,j_pos,:)


 ! ============ Non-volatile species first ============================
 ! NVABSOM - Only include fine OM! That is no EC and no coarse OM!

  do i = 1, NUM_NVABSOM  ! OA/OC for POC about 1.333
    ispec = NVABSOM_GROUP(i)

    ug_nonvol(i,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt

  end do

  ! ============ SOA species now, iteration needed ===================

  do iter = 1, NITER


      ! Fgas = G/(G+A) = 1/(1+K.COA)
      ! K = tabRTpL/(mw*gamma)

       do ispec = S1, S2

          !VBS Fgas(ispec,:) = 1.0/(1.0+tabRTpL(ispec,:)*COA(:)/avg_mw(:) )  
                                         !EXC *gamma(:,ispec)) )

          Fpart(ispec,:) = COA(:)/( COA(:)+tabCiStar(ispec,itemp(:)) )

          ug_semivol(ispec,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt &
                         * Fpart(ispec,:)

       end do ! ispec

       !ng(:,:) = max(SMALLFN, ng(:,:))

     ! New estimate of COA  (in ug/m3) and avg_mw (g/mole):
     ! (nb. xn in molecules/cm3)

       do k = K1,K2

         COA(k) = sum( ug_semivol(:,k) ) + sum( ug_nonvol(:,k) ) + BGND_OA(k)

         !VBS Nmoles = ( sum( Fpart(VOL,k) * xn(VOL,k) ) + sum( xn(NONVOLOC,k) ) &
         !VBS             + BGND_OC_ng(k)*ngm32molcc/250.0 ) &  ! Assumed MW
         !VBS          * xn2molem3
         !VBS avg_mw(k)    = 1.0e-6 * COA(k)/Nmoles 

       end do  !k
     ! ====================================================================

      if( DEBUG  .and. debug_flag ) then

         if( iter == NITER .and. seconds == 0 ) then
           write(unit=6,fmt="(a,i2,a,3i3,f7.2)") "Iteration ", Niter, &
                " ======== ",nmonth, nday, nhour, itemp(K2)
           write(unit=6,fmt="(a3,a15,3a10,a4,4a10)") "SOA","Species", "xn", &
               "Ci* ", "Ki"," ", "Fpart", "ng"

           do i = 1, NUM_NONVOLPCM
              ispec = NONVOLPCM_GROUP(i)
!              write(unit=6,fmt="(a4,i3,a15,es10.2,2f10.3,a4,es10.3,f13.4)")&
!                "NVOL", ispec,&
!                species(ispec)%name, xn(ispec,K2),-999.999, &
!                -999.999, " => ", Fpart(ispec,K2), 1000.0*ug_nonvol(i, K2)
              write(unit=6,fmt="(a4,i3,a15,es10.2,2f10.3)")&
                "NVOL", ispec,&
                species(ispec)%name, xn(ispec,K2),-999.999, &
                -999.999
           end do

           do ispec = S1,S2
              ! K = tabRTpL*COA/(mw*gamma) QUERY COA===!!!!
              !Ksoa = tabRTpL(ispec,K2)*COA(K2)/(avg_mw(K2))
              Ksoa = 1.0/tabCiStar(ispec,itemp(K2)) !just for printout
              write(unit=6,fmt="(a4,i3,a15,3es10.2,a4,es10.3,f13.4)") "SOA ",ispec,&
                species(ispec)%name, xn(ispec,K2), &
                 tabCiStar(ispec,itemp(K2)),& !VBStabVpsoa(ispec,298),
                  Ksoa, " => ", Fpart(ispec,K2), 1000.0*ug_semivol(ispec, K2)
          end do ! ispec
         end if 

          write(unit=6,fmt="(a,i2,f12.6)")  "COA: ", iter, COA(K2)

       end if ! DEBUG

   end do ! ITER

 ! The above iteration has now given new values to: 
 !
 ! 1) COA(1:nz)
 ! 3) Fgas(FIRST_SEMIVOL:LAST_SEMIVOL,1:nz)
 ! 4) Fpart(FIRST_SEMIVOL:LAST_SEMIVOL,1:nz)
 ! 5) ng(FIRST_SEMIVOL:LAST_SEMIVOL,1:nz)    
 !
 ! Note:     ng(FIRST_NONVOLOC:LAST_NONVOLOC,1:nz)
 !       and xn(1:LAST_SEMIVOL,1:nz)  are unaffected by the
 !       iteration.
 !=========================================================================

  ! Set Fgas for later chemistry, and eset 3-D fields

   Fgas(S1:S2,:)  = 1.0 - Fpart(S1:S2,:)
   Grid_COA(i_pos,j_pos,:)              = COA(:)
   Fgas3d(S1:S2,i_pos,j_pos,:) = Fgas(S1:S2,:) 

! PCM_F is for output only. Has MW 1 to avoid confusion with OC
! do not use ugC outputs, just ug

   xn(PART_OM_F,:)  =  COA(:) * ugC2xn * 12.0

    !Grid_SOA_Fgas(S1:S2, i_pos,j_pos,:)  = Fgas(S1:S2,:)
    !VBS Grid_avg_mw(i_pos,j_pos,:)       = avg_mw(:)
    !SOA_gamma(i_pos,j_pos,:,:)     = gamma(:,:)

  ! Outputs, ugC/m3

  ! Would like to be able to store also total OM (not only OC) 
  !    at least for some components. And for a "total" OM and/or OM2.5 and OM10. 
  !    Total TC and EC (and TC2.5, TC10, EC2.5 and EC10) would also be useful.
  !    Also perhaps the names of these species should reflect 
  !    that they are in units of C?   

   do k = K1, K2
     xn(PART_ASOA_OC,k)  = sum ( Fpart(ASOA,k) *xn(ASOA,k)  *species(ASOA)%carbons )
     xn(GAS_ASOA_OC,k)  = sum ( Fgas(ASOA,k) * xn(ASOA,k)  *species(ASOA)%carbons )
     xn(PART_BSOA_OC,k)  = sum ( Fpart(BSOA,k) *xn(BSOA,k)  *species(BSOA)%carbons )
     xn(GAS_BSOA_OC,k)  = sum ( Fgas( BSOA,k) *xn(BSOA,k)  *species(BSOA)%carbons )
     xn(NONVOL_FFUELOC25,k) = sum ( xn(NVFFUELOC25,k) *species(NVFFUELOC25)%carbons)
     xn(NONV_FFUELOC_COARSE,k) = sum ( xn(NVFFUELOCCO,k) *species(NVFFUELOCCO)%carbons)
     xn(NONVOL_WOODOC25,k)  = sum ( xn(NVWOODOC25,k)  *species(NVWOODOC25)%carbons )
     xn(NONVOL_FFIREOC25,k)  = sum ( xn(NVFFIREOC25,k)  *species(NVFFIREOC25)%carbons )
! want to have PART_FFUELOA25/FFIREOA/WOODOA_OC working also with nonvolatile POA emissions, test this hard coded version first
     xn(PART_FFUELOA25_OC,k)  = sum ( Fpart(SVFFUELOA25,k) *xn(SVFFUELOA25,k) *species(SVFFUELOA25)%carbons ) + &
     xn(NONVOL_FFUELOC25,k)
     xn(PART_WOODOA25_OC,k)  = sum ( Fpart(SVWOODOA25,k) *xn(SVWOODOA25,k)  *species(SVWOODOA25)%carbons )  + &
     xn(NONVOL_WOODOC25,k) 
     xn(PART_FFIREOA25_OC,k)  = sum ( Fpart(SVFFIREOA25,k) *xn(SVFFIREOA25,k)  *species(SVFFIREOA25)%carbons )  + &
       xn(NONVOL_FFIREOC25,k)
!Test for storing in ug/m3, Use with caution! 
     xn(PART_ASOA_OM,k)  = sum ( ug_semivol(ASOA,k) ) * ugC2xn * 12.0 
     xn(PART_BSOA_OM,k)  = sum ( ug_semivol(BSOA,k) ) * ugC2xn * 12.0 
! want to have PART_FFUELOA25/FFIREOA/WOODOA_OM working also with nonvolatile POA emissions, test this hard coded version first
     xn(PART_FFUELOA25_OM,k)  = sum ( ug_semivol(SVFFUELOA25,k) ) * ugC2xn * 12.0 + xn(NONVOL_FFUELOC25,k) * 1.25 * 12.0 ! factor 12.0 from M(OC25-components)=12 and M(OM-components)=1, OM/OC=1.25 assumed for Primary FFUELOC emissions
     xn(PART_WOODOA25_OM,k)  = sum ( ug_semivol(SVWOODOA25,k) ) * ugC2xn * 12.0 + xn(NONVOL_WOODOC25,k) * 1.7 * 12.0 ! OM/OC=1.7 assumed for Primary WOODOC and FFIRE emissions
     xn(PART_FFIREOA25_OM,k)  = sum ( ug_semivol(SVFFIREOA25,k) ) * ugC2xn * 12.0 + xn(NONVOL_FFIREOC25,k) * 1.7 * 12.0 

!HARDCODE
!   xn(AER_TBSOA,k)  =  xn(AER_BSOA,k)  ! Just in case TBSOA is wanted for kam
!...............................................................................
!Research   xn(PART_TBSOA_OC,k)  = &
!Research     Fpart( TERP_ng100,k) *xn(TERP_ng100,k)  *species(TERP_ng100)%carbons + &
!Research     Fpart( TERP_ug1,k) *xn(TERP_ug1,k)  *species(TERP_ug1)%carbons + &
!Research     Fpart( TERP_ug10,k) *xn(TERP_ug10,k)  *species(TERP_ug10)%carbons + &
!Research     Fpart( TERP_ug1e2,k) *xn(TERP_ug1e2,k)  *species(TERP_ug1e2)%carbons + &
!Research     Fpart( TERP_ug1e3,k) *xn(TERP_ug1e3,k)  *species(TERP_ug1e3)%carbons
!Research   xn(PART_IBSOA_OC,k)  = &
!Research     Fpart( ISOP_ng100,k) *xn(ISOP_ng100,k)  *species(ISOP_ng100)%carbons + &
!Research     Fpart( ISOP_ug1,k) *xn(ISOP_ug1,k)  *species(ISOP_ug1)%carbons + &
!Research     Fpart( ISOP_ug10,k) *xn(ISOP_ug10,k)  *species(ISOP_ug10)%carbons + &
!Research     Fpart( ISOP_ug1e2,k) *xn(ISOP_ug1e2,k)  *species(ISOP_ug1e2)%carbons + &
!Research     Fpart( ISOP_ug1e3,k) *xn(ISOP_ug1e3,k)  *species(ISOP_ug1e3)%carbons
!   xn(PART_SBSOA,k)  = &
!     Fpart( SESQ_ug1,k) *xn(SESQ_ug1,k)  *species(SESQ_ug1)%carbons + &
!    Fpart( SESQ_ug10,k) *xn(SESQ_ug10,k)  *species(SESQ_ug10)%carbons + &
!    Fpart( SESQ_ug1e2,k) *xn(SESQ_ug1e2,k)  *species(SESQ_ug1e2)%carbons + &
!    Fpart( SESQ_ug1e3,k) *xn(SESQ_ug1e3,k)  *species(SESQ_ug1e3)%carbons 
!...............................................................................
   end do
   xn(NONVOL_BGNDOC,:)    = ugC2xn * BGND_OC(:)      ! FAKE FOR NOW, 0.5 ug/m3 at surface

 ! for convencience:
! WARNING! Test changing to WOODOC25, FFIREOC25 and NONVOLOC25 may cause problems here, especially if coarse components are inlcuded later! But this is rather hard coded anyway...
   xn(PART_OC10,:)  =  xn(PART_ASOA_OC,:)+xn(PART_BSOA_OC,:)+xn(NONV_FFUELOC_COARSE,:) + &
                    xn(NONVOL_BGNDOC,:)+ &
                    xn(PART_FFUELOA25_OC,:)+xn(PART_WOODOA25_OC,:)+ &
                    xn(PART_FFIREOA25_OC,:)
!WARNING! The below will NOT work for NPNA (nonvolatile) type runs. These include fine and coarse OC in NONVOL_FFUELOC. So for these runs the FFUELOC-contribution has to be added separately!!!
!Test shifting to NONVOL_FFUELOC25!
   xn(PART_OC25,:)  =  xn(PART_ASOA_OC,:)+ xn(PART_BSOA_OC,:)+ &
                    xn(NONVOL_BGNDOC,:)+ &
                    xn(PART_FFUELOA25_OC,:)+xn(PART_WOODOA25_OC,:)+ &
                    xn(PART_FFIREOA25_OC,:)

   surfASOA  = xn2ugC* xn(PART_ASOA_OC,K2) ! sum ( Fpart(ASOA,K2) * xn(ASOA,K2)*species(ASOA)%carbons )
   surfBSOA  = xn2ugC* xn(PART_BSOA_OC,K2) ! sum ( Fpart(BSOA,K2) * xn(BSOA,K2)*species(BSOA)%carbons )
!
   surfFFUELOA25_OC  = xn2ugC* xn(PART_FFUELOA25_OC,K2) ! 
   surfWOODOA25_OC = xn2ugC* xn(PART_WOODOA25_OC,K2) !
   surfFFIREOA25_OC = xn2ugC* xn(PART_FFIREOA25_OC,K2) !

   surfBGNDOC  = BGND_OC(K2)
   surfOC25    = surfASOA+surfBSOA+surfFFUELOA25_OC+surfWOODOA25_OC+surfBGNDOC+surfFFIREOA25_OC !

   !/ Sum of Biogenics -----------------------

   if ( debug_flag .and. seconds == 0 ) then

       k=20
       write(unit=6,fmt="(a,3i3,2f7.2,f5.2,11es9.2)")"xns ug ", &
         nmonth, nday, nhour, &
         xn(PART_OM_F,20)*xn2ugC, & 
         COA(20), surfOC25,surfBGNDOC, surfASOA, surfBSOA, surfFFUELOA25_OC, &
         surfWOODOA25_OC, &
         surfFFIREOA25_OC!, &
!         COA(20), surfOC25,surfBGNDOC, surfASOA, surfBSOA,surfFFUELOA25_OC, &
!         surfOFFUELOA25_OC,surfFFUELOC25, surfWOODOA25, surfOWOODOA25, surfWOODOC25, &
!         surfFFIREOA25, surfOFFIREOA25!, &
!vbs      xn2ugC*xn(PART_IBSOA,k), xn2ugC*xn(PART_TBSOA,k), xn2ugC*xn(PART_SBSOA,k)

   endif

 end subroutine OrganicAerosol

end module OrganicAerosol_ml
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

