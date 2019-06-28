! <SOA_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module OrganicAerosol_mod

  ! Calculates the amount of condensible species in the gas and aerosol phases. 
  !
  ! References:
  !   B2012: Bergström, R. et al., Atmos. Chem. Physics, 2012, 12, 8499-8527
  !   S2012: Simpson, D. et al., Atmos. Chem. Physics, 2012, 12, 7825-7865 
  !   S2007: Simpson, D. et al., JGR, 2007, 
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
  ! Dave Simpson, August 2001 -- 2019
  ! Robert Bergström     2010 -- 2019
  ! 
  !--------------------------------------------------------------------------

  ! Functions + GridValues + PT only for BGNDOC
   use CheckStop_mod,  only: StopAll, CheckStop
   use ChemDims_mod,   only: NSPEC_SHL
   use ChemFields_mod, only: Fgas3d   ! stores 3-d  between time-steps
   use ChemFields_mod, only: xn_adv   ! for OM25_BGND
   use ChemSpecs_mod,  only: species, &   ! for molwts
                        S1 => FIRST_SEMIVOL , S2 => LAST_SEMIVOL

   use ChemGroups_mod  !XSOA , only :    &

   use Config_module,  only: PT, Pref, CHEMTMIN, CHEMTMAX, &
                             MasterProc,  & ! DebugCell, &
                             K2 => KMAX_MID, K1 => KCHEMTOP
   use Debug_module,   only: DebugCell, DEBUG  ! -> DEBUG%SOA
   use Functions_mod,  only: StandardAtmos_kPa_2_km !ds for use in Hz scaling
   use GridValues_mod, only: A_mid,B_mid, debug_proc, debug_li, debug_lj
   use Par_mod,        only: LIDIM => LIMAX, LJDIM => LJMAX, me
   use PhysicalConstants_mod, only : AVOG, RGAS_J 
   use ZchemData_mod,  only: itemp, xn => xn_2d, Fgas, Fpart
   use ZchemData_mod,  only: M   ! "M" = air density
   use SmallUtils_mod, only: find_index
   use TimeDate_mod,   only: current_date
   implicit none
   private


   !/-- subroutines

    public :: Init_OrganicAerosol
    public :: OrganicAerosol
    public :: Reset_OrganicAerosol ! resets bgnd and COA after advection
    public :: Reset_3dOrganicAerosol ! resets bgnd and COA after advection


   !/-- public

    logical, public, save :: ORGANIC_AEROSOLS = S1 > 0

   ! We store OM values in 3-D fields, to allow the next G/P partitioning
   ! calculation to  start off with values of COA, mw and Fgas which 
   ! are about right. Ensures that very few iterations are needed.

  real,public, save, allocatable, dimension(:,:,:) :: Grid_COA ! ug/m3
  real, private, allocatable,dimension(:,:) :: work

  real, private, allocatable, dimension(:), save :: &
        COA           & ! Org. aerosol, ug/m3  
                        ! (this version does not include EC as absorber)
       ,BGND_OC       & ! Background OC, assumed 0.20 ugC/m3 at surface
       ,BGND_OA         ! Assumed OA/OC=2, -> 0.4 ug/m3

   integer, private, save :: itot_bgnd = -999
   integer, private, save :: itot_om25 = -999, iadv_om25 = -999
   integer, private, save :: igrp_om25 = -999

  real, parameter, public :: SMALLFN  = 1.0e-20 ! Minimum value of ug allowed

   !/-- private

   ! ug = array for aerosol masses (ug/m3). Includes non-volatile compounds:
    real, private,allocatable, dimension(:,:), save :: ug_semivol 

   ! - use new NONVOLOC grpup to define:
    integer, private, save :: NUM_NONVOLPCM = 0 !  size(NONVOLPCM_GROUP)
    integer, private, save :: NUM_NVABSOM   = 0 !  size(NVABSOM_GROUP)
    integer, private, save :: nonvolpcm = -999, nvabsom  = -999
    real, private,allocatable, dimension(:,:), save :: ug_nonvol 

    real,  private, save, dimension(S1:S2,CHEMTMIN:CHEMTMAX) :: tabCiStar

    integer, private, save :: NITER = 2              ! No. iterations for Ksoa

   ! Need to convert aeros to ug/m3 or ugC/m3.  Since xn is in molecules/cm3
   ! we divide by AVOG to get mole/cm3, multiply by 1e6 to get mole/m3,
   ! by  mol. weight to get g/m3 and by 1e6 to get ug/m3

    real, private, parameter :: &
         xn2molem3 = 1.0/AVOG * 1.0e6          &  ! from molec/cm3 to mole/m3
       , xn2ug1MW  = xn2molem3 * 1.0 * 1.0e6   &  ! to ug/m3 ... for MW 1
       , xn2ugC   = xn2molem3 * 12.0 * 1.0e6   &  ! to ug/m3 ... for MW 12
       , ugC2xn   = 1/xn2ugC                   &  ! & back...
       , ug1MW2xn = 1/xn2ug1MW

   real, private, parameter :: molcc2ugm3 = 1.0e12/AVOG  !molecules/cc-> ug/m3
    ! same as xn2ug1MW, will harmonise later


   !/-- DEBUG variables
   !    Usually, DEBUG_SOA gives extra outputs. debug_flag is used to allow
   !    some extra outputs for a gven i,j - set in CTM model.

    character(len=20), public, save     :: soa_errmsg     = "ok"
    character(len=*), public, parameter :: SOA_MODULE_FLAG="VBS"

   contains
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine Init_OrganicAerosol(i,j,first_tstep,debug_flag)
   integer, intent(in) :: i,j
   logical, intent(in) :: first_tstep
   logical, intent(in) :: debug_flag
   integer :: is, ispec, it, k
   real, parameter :: kJ = 1000.0  
   real,allocatable, dimension(:), save :: p_kPa, h_km ! for standard atmosphere 
   real :: CiStar, dH  ! VBS params
   logical :: dbg  ! debug flag
   character(len=*), parameter :: dtxt = 'InitOrgAer:'

  ! Indices of Cstar and DeltaH groups in CM_ChemGroups
   integer, save :: igrp_Cstar = -1, igrp_DeltaH = -1
   logical, save :: first_call = .true.


   dbg = first_call .and. DEBUG%SOA > 0
   if(dbg .and. first_call) write(*,*) dtxt//"begin?", ORGANIC_AEROSOLS, S1, S2

   if( .not. ORGANIC_AEROSOLS  ) RETURN

   if ( first_call ) then ! =========================================
      itot_bgnd = find_index( 'OM25_bgnd', species(:)%name , any_case=.true.) 
      itot_om25 = find_index( 'OM25_p',    species(:)%name )  !NOTE CASE!
      iadv_om25 = itot_om25 - NSPEC_SHL
      igrp_om25 = find_index( 'OM25',   chemgroups(:)%name ) 
   
      igrp_Cstar  = find_index( 'CSTAR',  chemgroups_factors(:)%name ) 
      igrp_DeltaH = find_index( 'DELTAH', chemgroups_factors(:)%name ) 
   
      nonvolpcm = find_index( 'NONVOLPCM', chemgroups(:)%name ) 
      nvabsom   = find_index( 'NVABSOM',   chemgroups(:)%name ) 
   
     ! We need all of the above to be found. Check:
   
      if ( any( [ itot_bgnd, itot_om25, igrp_om25, igrp_Cstar, &
                    igrp_DeltaH, nonvolpcm, nvabsom             ] < 1 ) ) then
         print *, dtxt//'SOANEG:',itot_bgnd, itot_om25, &
                    igrp_om25, igrp_Cstar, igrp_DeltaH, nonvolpcm, nvabsom
         call StopAll(dtxt//'SOANEG')
      end if
   
      if( nonvolpcm > 0 ) NUM_NONVOLPCM = size(chemgroups(nonvolpcm)%specs)
      if( nvabsom   > 0 ) NUM_NVABSOM   = size(chemgroups(nvabsom)%specs)
   
      if( dbg ) then
        write(*,*) dtxt//"itot_bgnd, om25sum:", itot_bgnd,itot_om25,igrp_om25
        write(*,*) dtxt//"igrp Cstar,DeltaH:", igrp_Cstar, igrp_DeltaH
        write(*,*) dtxt//"nonvol,nv:",nonvolpcm,nvabsom,NUM_NONVOLPCM,NUM_NVABSOM
      end if
      call CheckStop( nvabsom < 1 .or. nonvolpcm < 1, dtxt//' Indices not found')

      allocate(COA(K1:K2))
      allocate(BGND_OC(K1:K2))
      allocate(BGND_OA(K1:K2))
      allocate(ug_semivol(S1:S2,K1:K2))
      allocate(Grid_COA(LIDIM,LJDIM,K1:K2))
      allocate(ug_nonvol(NUM_NVABSOM,K1:K2))
      allocate(p_kPa(K2), h_km(K2))

    !=========================================================================
    ! Set up background OM 

    ! Use Standard Atmosphere to get average heights of layers

       p_kPa(:) = 0.001*( A_mid(:)+B_mid(:)*Pref ) ! Pressure in kPa
       h_km     = StandardAtmos_kPa_2_km(p_kPa)
       BGND_OC(:)= 0.2 * 1.005 ! ng/m3 !!! will give 0.2 ugC/m3 at z=0 m

       do k = K1, K2
         BGND_OC(k) = BGND_OC(k) * exp( -h_km(k)/9.1 )
         if(dbg ) write(*,"(a,i4,2f8.3)") dtxt//"BGND_OC",k,h_km(k),BGND_OC(k)
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
   
      ! CSTAR_GROUP and DELTAH_GROUP and S1:S2 should be in same order
      ! sp we can loop from S1 to S2
   
      is = 0
      do ispec=S1,S2
         is = is + 1    ! Order
         if ( dbg ) then ! Checks that GenChem worked ok
            write(*,"(a,2i6,2f10.3,1x,a)") dtxt//"CSTAR ",  is, ispec, &
              chemgroups_factors(igrp_Cstar)%factors(is), &
              chemgroups_factors(igrp_DeltaH)%factors(is), &
              trim(species(ispec)%name)
            call CheckStop( &
                 ispec /= chemgroups_factors(igrp_Cstar)%species(is) .or. &
                 ispec /= chemgroups_factors(igrp_DeltaH)%species(is), &
                 dtxt//' Cstar DeltaH order wrong')
         end if
    
         CiStar = chemgroups_factors(igrp_Cstar)%factors(is)
         dH     = chemgroups_factors(igrp_DeltaH)%factors(is)
   
        ! C*-values are given for 298K according to most publications.
        ! so we tabukate for other temperatures
   
         do it=CHEMTMIN,CHEMTMAX
            tabCiStar(ispec,it) = CiStar * 298./it * &
                 exp( dH * kJ/RGAS_J * (1.0/298. - 1.0/it) )
         end do
       end do
   
       if ( dbg ) then 
         do is = S1, S2
            write(6,"(a,i4,1x,a20,f7.1,i3,8es10.2)") &
             dtxt//" Tab: MW, Carbons, C*:", is, adjustl(species(is)%name), &
              species(is)%molwt, nint(species(is)%carbons), & 
              tabCiStar(is,273), tabCiStar(is,303)
         end do
       end if

       first_call = .false.
    end if ! FIRST CALL 

    ! Now we continue with stuff that needs to be done of every i,j on the first
    ! time-step

    ! Initial values. Should not change except for semi-volatiles
    ! Note. Fgas3D has range S1:S2 only, whereas Fgas has 1:NSPEC_TOT
   
    ! The EMEP/ESX models denote concentrations in molec/cm3 as xn, although
    ! for particles the MW is a  'dummy' value (this doesn't matter), with
    ! e.g. BGND_OM having MW 24 to give a 2:1 ratio to the assumed unit carbon
    ! content.  OM25_p is just a helper species: the sum of the particle-phase
    ! OM25, and currently has MW 1 for simplicity.
    ! We need to convert our initial OA into xn, and on first call for each i,j
    
   if( first_tstep ) then 

      Fpart(:,:)         = 0.0
      Fpart(chemgroups(nonvolpcm)%specs,:)  = 1.0
      Fgas(:,:)         = max(0.0, 1.0 - Fpart(:,:) )

      COA(:) = BGND_OA(:)  ! Good starting estimate
      xn(itot_bgnd,:) = COA(:)/(molcc2ugm3*species(itot_bgnd)%molwt)  
      xn(itot_om25,:) = COA(:)/(molcc2ugm3*species(itot_om25)%molwt)  
      if( dbg) write(*,"(a,3f8.3,2es10.2)") dtxt//"COA, MW, xn TESTS ", &
        COA(K2), species(itot_om25)%molwt, species(itot_bgnd)%molwt, &
        xn(itot_bgnd,K2), xn(itot_om25,K2)

      if ( DEBUG%SOA  >1 ) then ! Invent some concs!
        do is = S1, S2
          xn(is,:) = 0.1/(molcc2ugm3*species(is)%molwt)   ! FAKE 0.1 ug/m3
        end do
       end if
    else   ! not first_tstep
       !NOT needed Fgas3d(S1:S2,i,j,:)=Fgas(S1:S2,:)
       ! since on 1st call we don't have any of the eg SOA compounds
       ! where Fgas affects reaction rates
 
       ! We need to set Fgas at start of each Runchem i,j loop, as it is
       ! used for rcemis:

       Fgas(S1:S2,:) = Fgas3d(S1:S2,i,j,:)     ! Semivolatiles only in 3D Fgas
       Fpart(S1:S2,:)  = 1-Fgas(S1:S2,:)

!      Fgas(NONVOLPCM_GROUP,:) = 0.0             !  not needed, shouldn't change
    end if ! first_tstep , set externally

  end subroutine Init_OrganicAerosol

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine OrganicAerosol(i_pos,j_pos,first_tstep,debug_flag)

   integer, intent(in) :: i_pos, j_pos
   logical, intent(in) :: first_tstep 
   logical, intent(in) :: debug_flag 
   character(len=*), parameter :: dtxt = 'RunOrgAer:'
   logical  :: dbg0, dbg1 

   integer :: is,  k, iter, ispec   ! loop variables 
   real :: Ksoa
   real :: tmpSum
   integer :: nmonth, nday, nhour, seconds
   character(len=99) :: sfmt

   nmonth = current_date%month
   nday   = current_date%day
   nhour  = current_date%hour
   seconds = current_date%seconds

   if(first_tstep .and. DebugCell ) write(*,*) dtxt//'DSOA ', me, &
        debug_proc, DebugCell,debug_flag

   dbg0 = DebugCell .and. DEBUG%SOA > 0
   dbg1 = DebugCell .and. DEBUG%SOA > 1
   if( dbg0 ) write(unit=*,fmt=*) dtxt//"Into SOA", DEBUG%SOA, S1, S2
   if( .not. ORGANIC_AEROSOLS ) then
     if(MasterProc) write(*,*) dtxt // "skipped. ORGANIC_AEROSOLS=F"
     RETURN
   end if

! Note that xn(SOA) is not strictly needed with the method we have, but it
! is a good first guess of the sum of the condensed phases, enables easy output
! and saves the need for iteration. 
!
! Remember also that Fgas is saved, so will initially have been set by the
! preceding call to OrganicAerosol for a different set of i,j.

 ! 1st guesses:

   COA(:)          =  Grid_COA(i_pos,j_pos,:)

 ! 2)/ Advected species

  if ( first_tstep ) then
    xn(itot_bgnd,:) = COA(:)/(molcc2ugm3*species(itot_bgnd)%molwt)  
  end if


 ! ============ Non-volatile species first ============================
 ! NVABSOM - Only includes fine OM; That is no EC and no coarse OM!

  tmpSum = 0.0
  do is = 1, NUM_NVABSOM  ! OA/OC for POC about 1.333
    ispec = chemgroups(nvabsom)%specs(is)

    ug_nonvol(is,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt
    tmpSum = tmpSum  + ug_nonvol(is,K2)

    if( dbg0 ) write(unit=*,fmt="(2a,f7.1,20es12.3)") &
      "NVABSOM SOA",&
      species(ispec)%name, species(ispec)%molwt,COA(K2), &
        xn(ispec,K2)*molcc2ugm3*species(ispec)%molwt, ug_nonvol(is,K2),tmpSum

  end do

  ! ============ SOA species now, iteration needed ===================

  do iter = 1, NITER


      ! Fgas = G/(G+A) = 1/(1+K.COA)
      ! K = tabRTpL/(mw*gamma)

       do ispec = S1, S2

          Fpart(ispec,:) = COA(:)/( COA(:)+tabCiStar(ispec,itemp(:)) )

          ug_semivol(ispec,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt &
                         * Fpart(ispec,:)
          tmpSum = tmpSum  + ug_semivol(ispec,K2)

          if( dbg1) write(unit=*,fmt="(i1, 2a,f7.1,20es12.3)") iter, " ABSOM: ",&
      species(ispec)%name, species(ispec)%molwt,COA(K2), ug_semivol(ispec,K2), tmpSum

       end do ! ispec

     ! New estimate of COA  (in ug/m3) and avg_mw (g/mole):
     ! (nb. xn in molecules/cm3)
     ! (nb. BGND_OA is in xn(itot_bgnd))

       do k = K1,K2

         COA(k) = sum( ug_semivol(:,k) ) + sum( ug_nonvol(:,k) )

       end do  !k
     ! ====================================================================

      if( dbg1 ) then

         if( iter == NITER .and. seconds == 0 ) then
           sfmt= "(a4,i3,1x,a15,3es10.2,a4,es10.3,f13.4)"

           write(unit=6,fmt="(a,i2,a,3i3,i4)") "Iteration ", Niter, &
                " ======== ",nmonth, nday, nhour, itemp(K2)
           write(unit=6,fmt="(a3,a15,3a10,a4,4a10)") "SOA","Species", "xn", &
               "Ci* ", "Ki"," ", "Fpart", "ng"

           do is = 1, NUM_NVABSOM
              ispec = chemgroups(nvabsom)%specs(is)
              write(unit=6,fmt=sfmt) "NVOL", ispec,&
                species(ispec)%name, xn(ispec,K2),-999.999, &
                -999.999, ' ', 1.0 , 1000.0*ug_nonvol(is,K2)
           end do

           do ispec = S1,S2
              !Ksoa = tabRTpL(ispec,K2)*COA(K2)/(avg_mw(K2))
              Ksoa = 1.0/tabCiStar(ispec,itemp(K2)) !just for printout
              write(unit=6,fmt=sfmt) "SOA ",ispec,&
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
 ! 1) COA(1:K2)
 ! 3) Fgas(FIRST_SEMIVOL:LAST_SEMIVOL,1:K2)
 ! 4) Fpart(FIRST_SEMIVOL:LAST_SEMIVOL,1:K2)
 ! 5) ng(FIRST_SEMIVOL:LAST_SEMIVOL,1:K2)    
 !
 ! Note:     ng(FIRST_NONVOLOC:LAST_NONVOLOC,1:K2)
 !       and xn(1:LAST_SEMIVOL,1:K2)  are unaffected by the
 !       iteration.
 !=========================================================================

  ! Set Fgas for later chemistry, and eset 3-D fields

  ! S1 > 0 if SOA used:
   Fgas(S1:S2,:)               = 1.0 - Fpart(S1:S2,:)
   Grid_COA(i_pos,j_pos,:)     = COA(:)
   Fgas3d(S1:S2,i_pos,j_pos,:) = Fgas(S1:S2,:) 

 end subroutine OrganicAerosol

 ! Reset_OrganicAerosol is called after Dry and Wet deposition, and should
 ! ensure that OMN25.......
 ! We store OM25_BGND  as a species despite its simple setting above. This makes
 ! makes summation of OM25 and PM25 components easier. Apart from memory increases
 ! the main practical problem is that the advection routine modifies xn_adv values
 ! away from those used to get Fgas. We reset here.

 subroutine Reset_OrganicAerosol(i_pos,j_pos,debug_flag)
   integer, intent(in) :: i_pos, j_pos
   logical, intent(in) :: debug_flag
   logical, save :: first_call = .true.
   integer :: k, n, itot
   real :: oldval
   logical :: dbg
   character(len=*), parameter :: dtxt = 'ResetOrgAer:'

   dbg = ( DEBUG%SOA>0 .and. debug_flag) 

   if ( debug_flag ) write(*,*) dtxt//"Skip Reset?", itot_bgnd 
   if( itot_bgnd < 1 ) then
     if ( debug_flag ) write(*,*) dtxt//"Skips Reset" 
     RETURN
   end if

   if ( dbg ) then
      oldval = xn(itot_bgnd,K2)  ! just for printout
      if(first_call) write(*,*) "Into Reset Organic Aerosol?",&
        itot_bgnd , first_call, size(chemgroups(igrp_om25)%specs)
    end if


   xn(itot_bgnd,:) = BGND_OC(:)* ugC2xn

   if( first_call.and. debug_flag ) write(*,"(a,i4,f7.3,9es12.3)") dtxt// &
     "itot_bgnd C: ", itot_bgnd, BGND_OC(K2), ugC2xn, xn(itot_bgnd,K2), oldval

   ! With SOA modelling some compounds are semivolatile and others non-volatile. If
   ! in a group XXX which asks for ugPM the latter's mass is correct. If semivolatile,
   ! we need to calculate the PM fraction and just add this.

   xn(itot_om25,:) = 0.0

   if( dbg ) write(*,*) dtxt//'OFSOA scaling ', xn2ug1MW, molcc2ugm3, &
      size(chemgroups(igrp_om25)%specs)
     !cf xn(itot_bgnd,:) = COA(:)/(molcc2ugm3*species(itot_bgnd)%molwt)  

   do n = 1, size(chemgroups(igrp_om25)%specs)

      itot  = chemgroups(igrp_om25)%specs(n)

     ! NOTE !  Assumes molwt is 1.0 for itot_om25

      xn(itot_om25,:) = xn(itot_om25,:) + &
           Fpart(itot,:) * xn(itot,:) * species(itot)%molwt
   
      if(  dbg ) then
        do k = K2, K2  ! K1, K2
           write(*,"(a,3i4,1x,a15,9es12.3)") dtxt//"OFSOA fac ", n, itot, k, &
             adjustl(species(itot)%name), Fpart(itot,k), &
             Fpart(itot,k) * xn(itot,k) * species(itot)%molwt*molcc2ugm3,&
             xn(itot_om25,k) * molcc2ugm3,  Grid_COA(i_pos,j_pos,k)
         end do
      end if
   end do ! n
   Grid_COA(i_pos,j_pos,:) = xn(itot_om25,:) * molcc2ugm3
   first_call = .false.

 end subroutine Reset_OrganicAerosol


 subroutine Reset_3dOrganicAerosol(debug_flag)
   logical, intent(in) :: debug_flag
   logical, save :: first_call = .true.
   integer :: k, n, itot, iadv
   logical :: dbg
   character(len=*), parameter :: dtxt = 'Reset3dOrg:'

   dbg = ( DEBUG%SOA>0 .and. debug_flag) 

   if( itot_bgnd < 1 ) RETURN

   if ( dbg ) then
     if ( first_call) then
       write(*,*) dtxt//"COMP START:", size(chemgroups(igrp_om25)%specs)
       allocate(work(LIDIM,LJDIM))
     end if
     work   = xn_adv(iadv_om25,:,:,K2)
   end if


   xn_adv(iadv_om25,:,:,:) =  0.0

   ! (use explicit K1:K2 since xn_adv has 1:K2, Fgas3d K1:K2)
   do n = 1, size(chemgroups(igrp_om25)%specs)
      itot = chemgroups(igrp_om25)%specs(n)
      iadv = itot - NSPEC_SHL
      if ( itot >= S1 .and. itot <= S2 ) then
        xn_adv(iadv_om25,:,:,K1:K2) = xn_adv(iadv_om25,:,:,K1:K2) + & 
          (1-Fgas3d(itot,:,:,K1:K2)) * &
          xn_adv(iadv,:,:,K1:K2) * species(itot)%molwt
      else
        xn_adv(iadv_om25,:,:,K1:K2) = xn_adv(iadv_om25,:,:,K1:K2) + & 
          xn_adv(iadv,:,:,K1:K2) * species(itot)%molwt
      end if
   end do

   if ( dbg ) then
     associate (diffs =>  xn_adv(iadv_om25,:,:,K2) - work(:,:) )
      write(*,'(a,4es12.3)') dtxt//"COMP:", work(debug_li,debug_lj), &
       xn_adv(iadv_om25,debug_li,debug_lj,K2), minval(diffs), maxval(diffs)
     end associate
    end if
    first_call = .false.

 end subroutine Reset_3dOrganicAerosol

end module OrganicAerosol_mod
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!TSTEMX program tester
!TSTEMX use ChemDims_mod, only: NCHEMRATES
!TSTEMX use ChemSpecs_mod, only: define_chemicals
!TSTEMX use ChemGroups_mod
!TSTEMX use ZchemData,  only :  Alloc1Dchem, itemp !, xn => xChem, Fgas, Fpart
!TSTEMX use OrganicAerosol_mod
!TSTEMX logical :: first_tstep = .true.
!TSTEMX real, dimension(5) :: zmid = [ 50.0, 150.0, 300.0, 500.0, 1000.0 ]
!TSTEMX integer, parameter :: NEMIS_BioNat = 4 !FAKE for now,
!TSTEMX call define_chemicals()
!TSTEMX call Init_Chemgroups()
!TSTEMX call Alloc1Dchem(1,NCHEMRATES,NEMIS_BioNat,debug_level=1)
!TSTEMX print *, "Into InitOA: "
!TSTEMX !call Init_OrganicAerosol(zmid,first_tstep,dbg=.true.)
!TSTEMX stop 'CRASHING - itemp not allocated yet :-(
!TSTEMX !itemp = 298
!TSTEMX !call OrganicAerosol(dbg=.true.)
!TSTEMX end program tester

