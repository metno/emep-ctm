! <Runchem_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
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
!> Runchem_ml.f90 - A component of the EMEP MSC-W Chemical transport Model
!!---------------------------------------------------------------------
!! Calls for routines calculating chemical and physical processes: 
!! irreversible and equilibrium chemistry, dry and wet deposition,
!! sea salt production, particle water etc.
!!---------------------------------------------------------------------

module RunChem_ml

  use AerosolCalls,     only: AerosolEquilib & !-> My_MARS, My_EQSAM, &
                             ,Aero_water, Aero_water_MARS 
  use My_Timing_ml,     only: Code_timer, Add_2timing,  &
                              tim_before, tim_after
  use AOD_PM_ml,        only: AOD_Ext
  use Aqueous_ml,       only: Setup_Clouds, prclouds_present, WetDeposition
  use Biogenics_ml,     only: setup_bio
  use CellMet_ml,       only: Get_CellMet, z0_out_ix, invL_out_ix
  use CheckStop_ml,     only: CheckStop, StopAll
  use Chemfields_ml,    only: xn_adv    ! For DEBUG 
  use Chemsolver_ml,    only: chemistry
  use ChemSpecs                         ! DEBUG ONLY
  use ColumnSource_ml,  only: Winds, getWinds
  use Config_module,    only: USES, & 
                              MasterProc, & 
                              KMAX_MID, END_OF_EMEPDAY, nstep,  &
                              AERO, USES, & ! need USES%EMISSTACKS and more 
                              USE_FASTJ, &
                              dt_advec, &  ! for Emergency
                              DEBUG_EMISSTACKS, & ! MKPS
                              DebugCell, DEBUG    ! RUNCHEM
  use DefPhotolysis_ml, only: setup_phot
  use DerivedFields_ml, only: f_2d
  use DryDep_ml,        only: drydep
  use DustProd_ml,      only: WindDust  !DUST -> USES%DUST
  use FastJ_ml,         only: setup_phot_fastj,phot_fastj_interpolate
  use GridValues_ml,    only: debug_proc, debug_li, debug_lj, i_fdom, j_fdom
  use Io_Progs_ml,      only: datewrite
  use MassBudget_ml,    only: emis_massbudget_1d
  use OrganicAerosol_ml,only: ORGANIC_AEROSOLS, OrganicAerosol, &
                              Init_OrganicAerosol, & 
                              Reset_OrganicAerosol, & 
                              SOA_MODULE_FLAG   ! ="VBS" or "NotUsed"
  use Pollen_ml,        only: Pollen_flux
  use Par_ml,           only: lj0,lj1,li0,li1, limax, ljmax,  &
                              gi0, gj0, me,IRUNBEG, JRUNBEG  !! for testing
  use PointSource_ml,    only: pointsources, get_pointsources
  use SeaSalt_ml,       only: SeaSalt_flux
  use Setup_1d_ml,      only: setup_1d, setup_rcemis, reset_3d
  use Setup_1dfields_ml,only: first_call, &
                              amk, rcemis, xn_2d  ! DEBUG for testing
  use SmallUtils_ml,    only: find_index
  use TimeDate_ml,      only: current_date,daynumber,print_date
!--------------------------------
  implicit none
  private

  public :: runchem
  private :: check_negs

contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine runchem()

!  local
  integer :: i, j
  integer :: errcode
  integer :: nmonth, nday, nhour     
  logical ::  Jan_1st
  logical ::  debug_flag    ! =>   Set true for selected i,j
  logical, save :: first_tstep = .true. ! J16 
  character(len=*), parameter :: sub='RunChem:'
! =============================
  nmonth = current_date%month
  nday   = current_date%day
  nhour  = current_date%hour


  Jan_1st    = ( nmonth == 1 .and. nday == 1 )


  if(ORGANIC_AEROSOLS.and.first_call) &
    call CheckStop(SOA_MODULE_FLAG == "NotUsed", & ! Just safety
                   "Wrong My_SOA? Flag is "// trim(SOA_MODULE_FLAG) )

!TEMPORARY HERE could be put in Met_ml
  if( (.not. first_call) .and. USES%PreADV)then
     call getWinds
  endif

  if(first_call)then
     z0_out_ix = find_index("logz0", f_2d(:)%subclass)
     invL_out_ix = find_index("invL", f_2d(:)%subclass)
  endif

! Processes calls 
  errcode = 0

  do j = 1, ljmax
    do i = 1, limax
! do j = lj0, lj1 !  ljmax
!   do i = li0, li1 ! 1, limax

      call Code_Timer(tim_before)

      !****** debug cell set here *******
      debug_flag =  .false.  
      if(DEBUG%RUNCHEM.and.debug_proc) then
        debug_flag = (debug_li==i .and. debug_lj==j) 
        DebugCell = debug_flag
        DEBUG%datetxt = print_date(current_date)
        if(debug_flag) write(*,*) "RUNCHEM DEBUG START!"
      end if
     !write(*,"(a,4i4)") "RUNCHEM DEBUG IJTESTS", debug_li, debug_lj, i,j
     !write(*,*) "RUNCHEM DEBUG LLTESTS", me,debug_proc,debug_flag

      ! Prepare some near-surface grid and sub-scale meteorology for MicroMet
      call Get_CellMet(i,j,debug_flag) 
      call Add_2timing(24,tim_after,tim_before,"Runchem:Get_CellMet ")

      ! we need to get the gas fraction of semivols:
      if ( ORGANIC_AEROSOLS ) call Init_OrganicAerosol(i,j,debug_flag)
      call Add_2timing(25,tim_after,tim_before,"Runchem:OrganicAerosol")

      call setup_1d(i,j)     ! Extracting i,j column data
      call Add_2timing(26,tim_after,tim_before,"Runchem:setup_1d")
      call setup_rcemis(i,j) ! Sets initial rcemis=0.0
      call Add_2timing(27,tim_after,tim_before,"Runchem:setup_rcemis ")
 
      if(USES%SEASALT)  &
        call SeaSalt_flux(i,j,debug_flag) ! sets rcemis(SEASALT_...)

      if(USES%DUST)     &
        call WindDust(i,j,debug_flag)     ! sets rcemis(DUST...)

      if ( USES%EMISSTACKS ) then
         if ( pointsources(i,j) ) call get_pointsources(i,j,DEBUG_EMISSTACKS)
      end if
    
      if(USES%POLLEN) &
        call Pollen_flux(i,j,debug_flag)

      call Setup_Clouds(i,j,debug_flag)

      call setup_bio(i,j)   ! Adds bio/nat to rcemis

      call emis_massbudget_1d(i,j)   ! Adds bio/nat to rcemis

      if(USE_FASTJ)then
!         call setup_phot_fastj(i,j,errcode,0)! recalculate the column
        !interpolate (intelligently) from 3-hourly values
         call  phot_fastj_interpolate(i,j,errcode)
      else
         call setup_phot(i,j,errcode)
      end if

      call CheckStop(errcode,"setup_photerror in Runchem") 

      if(DEBUG%RUNCHEM.and.debug_flag) &
        call datewrite("Runchem Pre-Chem", (/ rcemis(NO,20), &
             !rcemis(SHIPNOX,KMAX_MID), !hardcoded chemical indice are not
             ! defined for all chem schemes, and should usually be avoided
          rcemis(C5H8,KMAX_MID), xn_2d(NO,20),xn_2d(C5H8,20) /) )
      if(DEBUG%RUNCHEM) call check_negs(i,j,'A')

      if(ORGANIC_AEROSOLS) &
        call OrganicAerosol(i,j,first_tstep,debug_flag)  ! J16 first_tstep added
      if(DEBUG%RUNCHEM) call check_negs(i,j,'B')

      call Add_2timing(28,tim_after,tim_before,"Runchem:other setups")
!     if(DEBUG%RUNCHEM.and.debug_flag) &
!       call datewrite("RUNCHEM PRE-CHEM",(/xn_2d(PPM25,20),xn_2d(AER_BGNDOC,20)/))
!     !-------------------------------------------------
!     !-------------------------------------------------
!     !-------------------------------------------------

      if( .not. USES%NOCHEM) then
        call chemistry(i,j,DEBUG%RUNCHEM.and.debug_flag)
      else
        xn_2d(NSPEC_SHL+1:NSPEC_TOT,:) =  xn_2d(NSPEC_SHL+1:NSPEC_TOT,:)  &
                                         +rcemis(NSPEC_SHL+1:NSPEC_TOT,:)*dt_advec
      end if

      if(DEBUG%RUNCHEM) call check_negs(i,j,'C')
!     !-------------------------------------------------
!     !-------------------------------------------------
!     !-------------------------------------------------
      if(DEBUG%RUNCHEM.and.debug_flag)&
        call datewrite("Runchem Post-Chem",(/xn_2d(NO,20),xn_2d(C5H8,20)/))
      !_________________________________________________

      call Add_2timing(29,tim_after,tim_before,"Runchem:chemistry")
                
      !  Alternating Dry Deposition and Equilibrium chemistry
      !  Check that one and only one eq is chosen
      if(mod(nstep,2)/=0) then 
        call AerosolEquilib(debug_flag)
        call Add_2timing(30,tim_after,tim_before,"Runchem:AerosolEquilib")
        if(DEBUG%RUNCHEM) call check_negs(i,j,'D')
        !if(AERO%EQUILIB=='EMEP' ) call ammonium() 
        !if(AERO%EQUILIB=='MARS' ) call My_MARS(debug_flag)
        !if(AERO%EQUILIB=='EQSAM') call My_EQSAM(debug_flag) 
        call DryDep(i,j)
        call Add_2timing(31,tim_after,tim_before,"Runchem:DryDep")
        if(DEBUG%RUNCHEM) call check_negs(i,j,'E')
      else !do drydep first, then eq
        call DryDep(i,j)
        call Add_2timing(31,tim_after,tim_before,"Runchem:DryDep")
        if(DEBUG%RUNCHEM) call check_negs(i,j,'F')
        call AerosolEquilib(debug_flag)
        call Add_2timing(30,tim_after,tim_before,"Runchem:AerosolEquilib")
        if(DEBUG%RUNCHEM) call check_negs(i,j,'G')
        !if(AERO%EQUILIB=='EMEP' ) call ammonium() 
        !if(AERO%EQUILIB=='MARS' ) call My_MARS(debug_flag)
        !if(AERO%EQUILIB=='EQSAM') call My_EQSAM(debug_flag) 
      end if
      !????????????????????????????????????????????????????

      if(prclouds_present) then
        call WetDeposition(i,j,debug_flag)
        call check_negs(i,j,'H')
        call Add_2timing(32,tim_after,tim_before,"Runchem:WetDeposition")
      end if

     !Should be no further concentration changes due to emissions or deposition

      if(ORGANIC_AEROSOLS) call Reset_OrganicAerosol(i,j,debug_flag)


      !// Calculate Aerosol Optical Depth
      if(USES%AOD)  &
        call AOD_Ext(i,j,debug_flag)

      !  Calculates PM water: 1. for ambient Rh and T (3D)
      !  and for filter equlibration conditions (2D at surface) 
      !  T=20C and Rh=50% for comparability with gravimetric PM
      call Aero_water_MARS(i,j, debug_flag)

!.. Water from EQSAM .......
!     ambient = .false.  ! For Rh=50%
!     call Aero_water(i,j, ambient, debug_flag)                     
!     ambient = .true.  !  For real conditions (3D) 
!     call Aero_water(i,j, ambient, debug_flag)
                   
      call check_negs(i,j,'END')
      if(i>=li0.and.i<=li1.and.j>=lj0.and.j<=lj1) then

        call reset_3d(i,j)  ! DO NOT UPDATE BC. BC are frozen

      end if 

      call Add_2timing(33,tim_after,tim_before,"Runchem:post stuff")
      first_call = .false.   ! end of first call 
    end do ! j
  end do ! i
  first_tstep = .false.   ! end of first call  over all i,j

end subroutine runchem
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine check_negs(i,j,txt)
  integer, intent(in) :: i,j
  character(len=*), intent(in) :: txt
  integer :: n
  if ( any( xn_2d(:,KMAX_MID) < 0.0 ) ) then
         print *, txt, me,i_fdom(i),j_fdom(j)
         do n = 1, NSPEC_TOT
           if( xn_2d(n,KMAX_MID) < 0.0 ) print *, txt, xn_2d(n,KMAX_MID) 
         end do
         call StopAll( txt // 'STOPPED by check_negs in RunChem_ml')
  end if
end subroutine check_negs
endmodule RunChem_ml
