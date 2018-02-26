! <SOA_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2018 met.no
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
  !   B2012: Bergström, R. et al., Atmos. Chem. Physics, 2012, 12, 8499-8527
  !   S2012: Simpson, D. et al., Atmos. Chem. Physics, 2012, 12, 7825-7865 
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
  ! Dave Simpson, August 2001 -- 2016
  ! Robert Bergström     2010 -- 2015
  ! 
  !--------------------------------------------------------------------------

  ! Functions + GridValues + PT only for BGNDOC
   use CheckStop_ml, only : StopAll, CheckStop
   use ChemFields_ml,      only : Fgas3d   !  stores 3-d  between time-steps
   use ChemFields_ml,      only : xn_adv   ! J16 for OM25_BGND
   use ChemChemicals_ml,      only : species   ! for molwts
   use ChemSpecs_tot_ml,  S1 => FIRST_SEMIVOL , S2 => LAST_SEMIVOL

   use ChemGroups_ml  !XSOA , only :    &

   use Functions_ml, only: StandardAtmos_kPa_2_km !ds for use in Hz scaling
   use GridValues_ml, only: A_mid,B_mid, debug_proc, debug_li, debug_lj
   use Config_module,    only :  PT

   use Config_module,    only : CHEMTMIN, CHEMTMAX, &
                                    MasterProc, DEBUG, &
                                    K2 => KMAX_MID, K1 => KCHEMTOP
   use Par_ml,               only : LIDIM => LIMAX, LJDIM => LJMAX, me
   use PhysicalConstants_ml, only : AVOG, RGAS_J 
   use Setup_1dfields_ml,    only : itemp, xn => xn_2d, Fgas, Fpart
   use Setup_1dfields_ml,    only : amk   ! "M" = air density
   use SmallUtils_ml,        only : find_index
   use TimeDate_ml,          only : current_date
   implicit none
   private


   !/-- subroutines

    public   :: Init_OrganicAerosol
    public   :: OrganicAerosol
    public   :: Reset_OrganicAerosol ! FSOA - resets bgnd and COA after advection


   !/-- public

    logical, public, save :: ORGANIC_AEROSOLS = S1 > 0

   ! We store some values in 3-D fields, to allow the next G/P partitioning
   ! calculation to  start off with values of COA, mw and Fgas which 
   ! are about right. Ensures that very few iterations are needed.

   !  real,public, save, dimension(S1:S2,LIDIM,LJDIM,K1:K2) :: &
   !            Grid_SOA_Fgas           !EXC Grid_SOA_gamma

  real,public, save, allocatable, dimension(:,:,:) :: Grid_COA

  real, private, allocatable, dimension(:), save :: &
        COA           & ! Org. aerosol, ug/m3  
                        ! (this version does not include EC as absorber)
       ,BGND_OC       & ! Background OC, assumed 0.20 ugC/m3 at surface
       ,BGND_OA         ! Assumed OA/OC=2, -> 0.4 ug/m3

   integer, private, save :: itot_bgnd = -999
   integer, private, save :: itot_om25 = -999
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
   subroutine Init_OrganicAerosol(i,j,debug_flag)
   integer, intent(in) :: i,j
   logical, intent(in) :: debug_flag
   integer :: is,  it, k
   real, parameter :: kJ = 1000.0  
   real,allocatable, dimension(:), save :: p_kPa, h_km ! for standard atmosphere 
   character(len=*), parameter :: dtxt = 'InitOrgAero'

   logical, save :: first_call = .true.

   if( .not. ORGANIC_AEROSOLS ) then
     if(MasterProc) write(*,*) dtxt // "skipped. ORGANIC_AEROSOLS=F"
     RETURN
   end if

   if( first_call ) then

     itot_bgnd = find_index( 'OM25_BGND', species(:)%name ) 
     itot_om25 = find_index( 'OM25_P',  species(:)%name ) 
     igrp_om25 = find_index( 'OM25',  chemgroups(:)%name ) 

      nonvolpcm = find_index( 'NONVOLPCM', chemgroups(:)%name ) 
      nvabsom   = find_index( 'NVABSOM',   chemgroups(:)%name ) 
      if( nonvolpcm > 0 ) NUM_NONVOLPCM = size(chemgroups(nonvolpcm)%specs)
      if( nvabsom   > 0 ) NUM_NVABSOM   = size(chemgroups(nvabsom)%specs)

      if( MasterProc ) then
         write(*,*) dtxt//"itot_bgnd, om25sum = ", itot_bgnd, itot_om25, igrp_om25
         write(*,*) dtxt // "nonvol,nv ", nonvolpcm, nvabsom,  NUM_NONVOLPCM, NUM_NVABSOM
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

       p_kPa(:) = 0.001*( A_mid(:)+B_mid(:)*101325.0 ) ! Pressure in kPa
       h_km     = StandardAtmos_kPa_2_km(p_kPa)
       BGND_OC(:)= 0.2 * 1.005 ! ng/m3 !!! will give 0.2 ugC/m3 at z=0 m

       do k = K1, K2
            BGND_OC(k) = BGND_OC(k) * exp( -h_km(k)/9.1 )
            if(DEBUG%SOA .and. MasterProc ) write(*,"(a,i4,2f8.3)") &
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
               write(6,"(a,i4,1x,a20,f7.1,i3,8es10.2)") &
                " Tab SOA: MW, Carbons, C*:", is, adjustl(species(is)%name), &
                 species(is)%molwt, species(is)%carbons, & 
                 tabCiStar(is,273), tabCiStar(is,303)
            end do
         end if

       !+ initial guess (1st time-step only)
       ! Fgas3D is only defined for the semivol stuff, so no need for nonvol here
       ! We need to assume something on 1st time-step though:
       ! nb DONE in Chem_ml:  Fgas3d = 1.0

       ! Initial values. Should not change except for semi-volatiles
       ! Note. Fgas3D has range S1:S2 only, whereas Fgas has 1:NSPEC_TOT
       

        Fpart(:,:)         = 0.0
        Fpart(chemgroups(nonvolpcm)%specs,:)  = 1.0
        Fgas(:,:)         = max(0.0, 1.0 - Fpart(:,:) )

        !NOT needed Fgas3d(S1:S2,i,j,:)=Fgas(S1:S2,:)  ! J29
        ! since on 1st call we don't have any of the eg SOA compounds
        ! where Fgas affects reaction rates

        first_call = .false.

    end if ! first_call

    ! We need to set Fgas at start of each Runchem i,j loop, as it is
    ! used for rcemis:

      Fgas(S1:S2,:) = Fgas3d(S1:S2,i,j,:)     ! Semivolatiles only in 3D Fgas
      Fpart(S1:S2,:)  = 1-Fgas(S1:S2,:)

!      Fgas(NONVOLPCM_GROUP,:) = 0.0             !  not needed, shouldn't change

  end subroutine Init_OrganicAerosol

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine OrganicAerosol(i_pos,j_pos,first_tstep,debug_flag)

   integer, intent(in) :: i_pos, j_pos
   logical, intent(in) :: first_tstep 
   logical, intent(in) :: debug_flag 
   character(len=*), parameter :: dtxt = 'RunOrgAero'

   integer :: i,  k, iter, ispec   ! loop variables 
   real :: Ksoa
   integer :: nmonth, nday, nhour, seconds
   logical :: dbg
   dbg = ( DEBUG%SOA .and. debug_flag) 

   nmonth = current_date%month
   nday   = current_date%day
   nhour  = current_date%hour
   seconds = current_date%seconds

   if( dbg ) write(unit=*,fmt=*) "Into SOA"
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
 ! NVABSOM - Only include fine OM! That is no EC and no coarse OM!

  do i = 1, NUM_NVABSOM  ! OA/OC for POC about 1.333
    ispec = chemgroups(nvabsom)%specs(i)

    ug_nonvol(i,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt

    if( dbg) write(unit=*,fmt="(2a,f7.1,20es12.3)") &
      "NVABSOM SOA",&
      species(ispec)%name, species(ispec)%molwt,COA(20), &
        xn(ispec,20)*molcc2ugm3*species(ispec)%molwt, ug_nonvol(i,20)

  end do

  ! ============ SOA species now, iteration needed ===================

  do iter = 1, NITER


      ! Fgas = G/(G+A) = 1/(1+K.COA)
      ! K = tabRTpL/(mw*gamma)

       do ispec = S1, S2

          Fpart(ispec,:) = COA(:)/( COA(:)+tabCiStar(ispec,itemp(:)) )

          ug_semivol(ispec,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt &
                         * Fpart(ispec,:)

       end do ! ispec

     ! New estimate of COA  (in ug/m3) and avg_mw (g/mole):
     ! (nb. xn in molecules/cm3)
     ! (nb. BGND_OA is in xn(itot_bgnd))

       do k = K1,K2

         COA(k) = sum( ug_semivol(:,k) ) + sum( ug_nonvol(:,k) )

       end do  !k
     ! ====================================================================

      if( dbg ) then

         if( iter == NITER .and. seconds == 0 ) then
           write(unit=6,fmt="(a,i2,a,3i3,i4)") "Iteration ", Niter, &
                " ======== ",nmonth, nday, nhour, itemp(K2)
           write(unit=6,fmt="(a3,a15,3a10,a4,4a10)") "SOA","Species", "xn", &
               "Ci* ", "Ki"," ", "Fpart", "ng"

           do i = 1, NUM_NONVOLPCM
              ispec = chemgroups(nonvolpcm)%specs(i)
              write(unit=6,fmt="(a4,i3,a15,es10.2,2f10.3)")&
                "NVOL", ispec,&
                species(ispec)%name, xn(ispec,K2),-999.999, &
                -999.999
           end do

           do ispec = S1,S2
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

   Fgas(S1:S2,:)               = 1.0 - Fpart(S1:S2,:)
   Grid_COA(i_pos,j_pos,:)     = COA(:)
   Fgas3d(S1:S2,i_pos,j_pos,:) = Fgas(S1:S2,:) 

 end subroutine OrganicAerosol

 ! Reset_OrganicAerosol is called after Dry and Wet deposition, and should
 ! ensure that OMN25.......
 ! We store OM25_BGND  as a species despite its simple setting above. This makes
 ! makes summation of OM25 and PM25 components easier. Apart from memory increases
 ! the main practuiacl problem is that the advection routine modifies xn_adv values
 ! away from those used to get Fgas. We reset here.

 subroutine Reset_OrganicAerosol(i_pos,j_pos,debug_flag)
   integer, intent(in) :: i_pos, j_pos
   logical, intent(in) :: debug_flag
   logical, save :: first_call = .true.
   integer :: k, n, itot
   real :: J16tmp
   logical :: dbg
   dbg = ( DEBUG%SOA .and. debug_flag) 

   if ( debug_flag ) write(*,*) "Skip Reset Organic Aerosol?", itot_bgnd 
   if( itot_bgnd < 1 ) then
     if ( debug_flag ) write(*,*) "Skips Reset Organic Aerosol" 
     RETURN
   end if


   if ( first_call .and. debug_proc ) then
      J16tmp = xn(itot_bgnd,20)  ! just for printout
      write(*,*) "Into Reset Organic Aerosol?",&
        itot_bgnd , first_call, size(chemgroups(igrp_om25)%specs)
    end if


   xn(itot_bgnd,:) = BGND_OC(:)* ugC2xn

   if( first_call.and. debug_proc ) write(*,"(a,i4,f7.3,9es12.3)") &
     "itot_bgnd C = ", itot_bgnd, BGND_OC(20), ugC2xn, xn(itot_bgnd,20), J16tmp

   ! With SOA modelling some compounds are semivolatile and others non-volatile. If
   ! in a group XXX which asks for ugPM the latter's mass is correct. If semivolatile,
   ! we need to calculate the PM fraction and just add this.

   xn(itot_om25,:) = 0.0

   if(  dbg  ) write(*,*) 'OFSOA scaling ', xn2ug1MW, molcc2ugm3, &
      size(chemgroups(igrp_om25)%specs)
     !cf xn(itot_bgnd,:) = COA(:)/(molcc2ugm3*species(itot_bgnd)%molwt)  

   do n = 1, size(chemgroups(igrp_om25)%specs)

        itot  = chemgroups(igrp_om25)%specs(n)

       ! NOTE !  Assumes molwt is 1.0 for itot_om25
        xn(itot_om25,:) = xn(itot_om25,:) + &
             Fpart(itot,:) * xn(itot,:) * species(itot)%molwt
   
   if(  dbg ) then
          do k = K1, K2
          !do k = K2, K2
            write(*,"(a,3i4,1x,a15,9es12.3)") "OFSOA fac ", n, itot, k, &
              adjustl(species(itot)%name), Fpart(itot,k), &
              Fpart(itot,k) * xn(itot,k) * species(itot)%molwt*molcc2ugm3,&
              xn(itot_om25,k) * molcc2ugm3,  Grid_COA(i_pos,j_pos,k)
          end do
        end if
   end do ! n
   Grid_COA(i_pos,j_pos,:) = xn(itot_om25,:) * molcc2ugm3
   first_call = .false.

 end subroutine Reset_OrganicAerosol

end module OrganicAerosol_ml
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

