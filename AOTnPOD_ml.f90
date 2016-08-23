! <AOTx_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2011 met.no
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
module AOTx_ml
  use CheckStop_ml, only : checkStop
  use Chemfields_ml, only : xn_adv, cfac
  use ChemSpecs_adv_ml, only : IXADV_O3
  use GridValues_ml, only : debug_li, debug_lj
  use Io_Progs_ml,   only : datewrite
  use LandDefs_ml,   only : LandType
  use Landuse_ml,    only : WheatGrowingSeason
  use LocalVariables_ml, only : L, Grid, Sub
  use MetFields_ml, only: zen
  use ModelConstants_ml, only : dt_advec, KMAX_MID &
     ,DEBUG_AOT  & 
     ,PPBINV ! 1.0e9, for conversion from mixing ratio to ppb
  use OwnDataTypes_ml, only : TXTLEN_DERIV, TXTLEN_SHORT
  use Par_ml, only : MAXLIMAX, MAXLJMAX, limax, ljmax
  use TimeDate_ml, only : current_date, jday => effectivdaynumber
  implicit none
  private

  public :: Calc_AOTx          ! called from AddMosaicOutput, My_DryDep
  public :: Calc_POD           ! called from AddMosaicOutput
  public :: Calc_GridAOTx      ! called from Derived_ml

! Limit of daylight zenith angle for AOTs
  integer, private,  parameter :: AOT_HORIZON  = 89         
  integer, public, parameter:: STARTMONTH_FOREST=4,ENDMONTH_FOREST=9&
                                ,STARTMONTH_CROPS=5,ENDMONTH_CROPS=7 ! EU only!
  logical, parameter, private :: T=.true., F=.false.

! VEGO3 definitions for PODY(formerly AFstY) and AOTX
!
! PODY vs Fst - not that the accumulated PODY  is processed here.
! the instantaneous Fst is set as for Canopy O3 in METCONCS
          ! N.B. AOTs have several definitions. We usually want
          ! the ICP-veg Mapping Manual (MM) ones. Other
          ! possibilities are EU (8-20daytime) or UN (May-July for
          ! crops)
   !================== 
    type, public:: O3cl
       character(len=TXTLEN_DERIV) :: name ! e.g. POD1_IAM_DF
       character(len=TXTLEN_SHORT) :: class ! POD or AOT
       real    :: Threshold     ! Threshold or CL, e.f. AOTx or AFstY
       character(len=TXTLEN_SHORT) :: defn !  MM or EU definitions
       character(len=TXTLEN_SHORT) :: TXTLC !  CF, DF, IAM_CF etc.
       logical :: RelSGS      ! true if accumulation period is relative to
                              ! start of growing season (SGS)
                              ! can be false it fixed, e.g. April 1st
       integer :: SAccPeriod  ! Start of accumulation period, either rel 
                              ! to SGS or day number, days
       integer :: EAccPeriod  ! End ...., days
    end type 

    type(O3cl), public, allocatable, dimension(:) :: &
     VEGO3_OUTPUTS

    type(O3cl), public, parameter, dimension(26) :: &
     VEGO3_DEFS    =  (/ &
         ! name               class X/Y   defn   txtLC relSGS Sacc Eacc
   O3cl( "POD1_IAM_DF   ",   "POD", 1.0,  "- ", "IAM_DF",F,0,999 ), & 
   O3cl( "POD0_IAM_DF   ",   "POD", 0.0,  "- ", "IAM_DF",F,0,999 ), &
   O3cl( "POD1_IAM_MF   ",   "POD", 1.0,  "- ", "IAM_MF",F,0,999 ), & 
   O3cl( "POD0_IAM_MF   ",   "POD", 0.0,  "- ", "IAM_DF",F,0,999 ), &
   O3cl( "POD1_DF       ",   "POD", 1.0,  "- ", "DF    ",F,0,999 ), &
   O3cl( "POD1_CF       ",   "POD", 1.0,  "- ", "CF    ",F,0,999 ), &
   O3cl( "POD3_TC       ",   "POD", 3.0,  "- ", "TC    ",F,0,999 ), & !
! 30days, 15 before mid, 15 after
   O3cl( "POD3_TC30d    ",   "POD", 3.0,  "- ", "TC    ",T,30,60 ), & 
! 55days, 15 before mid, 40 after
   O3cl( "POD3_TC55d    ",   "POD", 3.0,  "- ", "TC    ",T,30,85 ), & 
   O3cl( "POD6_TC       ",   "POD", 3.0,  "- ", "TC    ",F,0,999 ), & !
   O3cl( "POD6_TC55d    ",   "POD", 3.0,  "- ", "TC    ",T,30,85 ), & ! 55d
   O3cl( "POD3_IAM_CR   ",   "POD", 3.0,  "- ", "IAM_CR",F,0,999 ), &
   O3cl( "POD3_IAM_CR30d",   "POD", 3.0,  "- ", "IAM_CR",T,30,60 ), &
   O3cl( "POD3_IAM_CR55d",   "POD", 3.0,  "- ", "IAM_CR",T,30,85 ), &
! WARNING - POD6 is not recommended for European-scale IAM modelling
! Way too uncertain!
   O3cl( "POD6_IAM_CR   ",   "POD", 6.0,  "- ", "IAM_CR",F,0,999 ), & ! WARN!
   O3cl( "POD6_IAM_CR30d",   "POD", 6.0,  "- ", "IAM_CR",T,30,60 ), & ! WARN!
   O3cl( "POD6_IAM_CR55d",   "POD", 6.0,  "- ", "IAM_CR",T,30,85 ), & ! WARN!
   O3cl( "POD1_IAM_CR   ",   "POD", 1.0,  "- ", "IAM_CR",F,0,999 ), &
   O3cl( "POD0_IAM_CR   ",   "POD", 0.0,  "- ", "IAM_CR",F,0,999 ), &
! has 90 day GS
   O3cl( "MMAOT40_TC    ","AOT", 40.0, "MM", "TC    ",F,0,999 ), & 
   O3cl( "MMAOT40_IAM_DF","AOT", 40.0, "MM", "IAM_DF",F,0,999 ), & !
   O3cl( "MMAOT40_IAM_MF","AOT", 40.0, "MM", "IAM_MF",F,0,999 ), & !
   O3cl( "MMAOT40_IAM_CR","AOT", 40.0, "MM", "IAM_CR",F,0,999 ), &
! IAM_CR, we use 3m O3
   O3cl( "EUAOT40_Crops ", "AOT", 40.0, "EU", "IAM_CR",F,0,999 ), &  
! IAM_DF, we use 3m O3
   O3cl( "EUAOT40_Forests", "AOT", 40.0, "EU", "IAM_DF",F,0,999 ), &  
   O3cl( "MMAOT40_IAM_WH ","AOT", 40.0, "MM", "IAM_WH",F,0,999 ) &
    /) !NB -last not found. Could just be skipped, but kept
       !to show behaviour of code

contains
 !=========================================================================
  subroutine Calc_AOTx(iO3cl,iLC, aot, debug_flag, debug_txt )
    logical, intent(in) :: debug_flag
    character(len=*), intent(in), optional :: debug_txt
    integer, intent(in) :: iO3cl,iLC
    real, intent(out)    :: aot

    real    :: o3
    integer :: i,j
    real :: X

    i = Grid%i
    j = Grid%j
    X = vego3_outputs(iO3cl)%Threshold


    aot = 0.0

    if (   LandType(iLC)%is_forest .and. &
         (   current_date%month<STARTMONTH_FOREST&
         .or.current_date%month>ENDMONTH_FOREST)   )  then  
       ! for both EU and MM-EMEP so far. Reconsider for MM
             return
    end if

   ! If night, or outside growing season, we simply exit with aot=0
   ! EU AOT for 8:00 -- 20:00 CET, is really 8:00 -- 19:59 CET 
   ! Or:        7:00 -- 18:59 UTC
   ! (nb hour is integer value)

    if ( vego3_outputs(iO3cl)%defn == "EU" .and.  &
         (current_date%hour  <  7 .or. current_date%hour  > 18 )) then
          return

    else if ( Grid%izen >= AOT_HORIZON ) then  !UN or MM use daylight
          return
    end if

  ! the wheat growing season is based upon the latitude function
    if ( vego3_outputs(iO3cl)%defn == "MM" ) then
        if ( DEBUG_AOT .and. debug_flag .and. present( debug_txt )) then
          write(*,*) 'AOTiO3CL ', iLC, iO3cl, Sub(iLC)%SGS, jday
          write(*,*) 'AOTiO3CLL', LandType(iLC)%is_crop, &
              LandType(iLC)%is_iam, WheatGrowingSeason(i,j) 
        end if
       if ( jday < Sub(iLC)%SGS .or.  &
            jday > Sub(iLC)%EGS ) then
          
             return
        end if

      o3 = Sub(iLC)%cano3_ppb   ! canopy top O3 for MM

    else if ( vego3_outputs(iO3cl)%defn == "EU" ) then
      if (   LandType(iLC)%is_crop .and. &
            ( current_date%month<STARTMONTH_CROPS&
              .or.current_date%month>ENDMONTH_CROPS) ) then
           return
      end if

      o3 = Grid%surf_o3_ppb 

      else 
         if ( DEBUG_AOT )  call CheckStop("AOT not MM or EU!")
    end if

  ! Temporary setup - accumulation period for NS Clover for Mills et al paper
!    if ( defn == "CV" ) then
!      if ( current_date%month < 5 .or.  &
!           current_date%month == 5 .and.  current_date%day < 15) .or.
!           current_date%month == 5 .and.  current_date%day < 15) .or.
!
!           ) then
!             return
!      end if
!    end if


    !========== Calculate ========================

    if ( o3>X ) then  ! Add AOT, scaling for time-fraction

        aot = (o3-X)  ! dt_advec takes care of 3600s elsewhere

    end if
    if ( DEBUG_AOT .and. debug_flag .and. present( debug_txt )) then
       call datewrite("AOTxdebug"//trim(debug_txt) // "defn:" // &
         trim(vego3_outputs(iO3cl)%defn), iLC, &
           (/ real(Grid%izen), X,  o3, aot /) )
    end if

  end subroutine Calc_AOTx

 !=========================================================================
 !
  function Calc_GridAOTx( iX, debug_flag, debug_txt ) result (aot)

    !/-- Calcuates AOT values for input threshold. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )
    !    Only relevant in ozone models, so far....

    integer, intent(in) :: iX  ! usually 40.0
    logical, intent(in) :: debug_flag
    character(len=*), intent(in), optional :: debug_txt

    real, dimension(limax,ljmax) :: aot

    integer :: izen                    ! integer of zenith angle
    real :: o3, o3_ref                 ! Ozone (ppb) - needed if AOTs
    integer :: i, j 

    aot = 0.0

    !--------- ONLY April-Sept -------------------------
    if(   current_date%month<STARTMONTH_FOREST&
      .or.current_date%month>ENDMONTH_FOREST) return
    !--------- ONLY April-Sept -------------------------

    do i=1, limax
        do j=1, ljmax

           izen = max(1,int( zen(i,j) + 0.5))

           if ( izen < AOT_HORIZON ) then
                o3 = xn_adv(IXADV_O3,i,j,KMAX_MID) &
                     * cfac(IXADV_O3,i,j) * PPBINV

                aot(i,j) = max( o3 - iX , 0.0 )   ! Definition of AOTs

           end if
        end do
    end do
    if ( DEBUG_AOT .and. debug_flag .and. present( debug_txt )) then
       o3_ref = xn_adv(IXADV_O3,debug_li,debug_lj,KMAX_MID) * PPBINV
       o3   = o3_ref * cfac(IXADV_O3,debug_li,debug_lj)
       call datewrite("CalcGridAOT:"//debug_txt, (/ zen(debug_li,debug_lj), &
            o3_ref, o3, aot(debug_li,debug_lj) /) )
    end if

  end function Calc_GridAOTx

!=========================================================================
  subroutine Calc_POD(iO3cl,iLC, pod, debug_flag, debug_txt )
    logical, intent(in) :: debug_flag
    character(len=*), intent(in), optional :: debug_txt
    integer, intent(in) :: iO3cl,iLC
    real, intent(out)    :: pod

    real    :: o3
    integer :: i,j, spod, epod
    real :: Y

    pod = 0.0

    if ( Grid%izen >= 90 ) then  !UN or MM use daylight
        if ( DEBUG_AOT .and. debug_flag ) write(*,*) "PODxZen ",&
           trim(VEGO3_OUTPUTS(iO3cl)%name), iO3cl, jday, Grid%izen
        return
    end if

    ! Start and end of POD accumulation period:
    ! ( Some PODs are only for a specific period, different to 
    !   growing season (SGS), notably TC (wheat) which has 55 days.)

    if ( VEGO3_OUTPUTS(iO3cl)%RelSGS .eqv. .true. )  then
      spod = Sub(iLC)%SGS + VEGO3_OUTPUTS(iO3cl)%SAccPeriod
      epod = Sub(iLC)%SGS + VEGO3_OUTPUTS(iO3cl)%EAccPeriod
    else
      spod = Sub(iLC)%SGS
      epod = Sub(iLC)%EGS
    end if
    if ( jday < spod .or. jday > epod ) then
        if ( DEBUG_AOT .and. debug_flag ) write(*,*) "PODxJday ",&
           trim(VEGO3_OUTPUTS(iO3cl)%name), iO3cl, jday, spod, epod
       return
    end if

    i = Grid%i
    j = Grid%j
    Y = vego3_outputs(iO3cl)%Threshold ! nmole/m2/s

    o3 = Sub(iLC)%cano3_ppb   ! canopy top O3 for MM

   ! Add fluxes if Y exceeded:

     pod  = max(Sub(iLC)%FstO3 - Y,0.0)


    if ( DEBUG_AOT .and. debug_flag .and. present( debug_txt )) then
        write(6,"(a,L2,3i4)") "PODACC ", VEGO3_OUTPUTS(iO3cl)%RelSGS, &
           VEGO3_OUTPUTS(iO3cl)%SAccPeriod, VEGO3_OUTPUTS(iO3cl)%EAccPeriod, &
           jday
       call datewrite("YYY"//trim(debug_txt) // "defn:" // &
         trim(vego3_outputs(iO3cl)%defn), iLC, &
           (/ real(Grid%izen), Y,  o3, Sub(iLC)%FstO3, pod /) )
    end if

  end subroutine Calc_POD 

 !=========================================================================


! Will move SOMO and maxo3 here in future versions
!current definitions:
!SOMO35: daily max is found over 00:00 to 24:00. (not emepday)
!SOMO35: accumulated over one year
!D2_MAXO3 :  daily max is found over an EMEPDAY
!D2_MAXO3 : accumulated into yearly_output from April to September


end module AOTx_ml
