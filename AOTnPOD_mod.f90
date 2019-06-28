! <AOTnPOD_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
!> <AOTx_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!  **********************************************************************! 

module AOTx_mod
  use CheckStop_mod,  only : checkStop, StopAll
  use Chemfields_mod, only : xn_adv, cfac
  use ChemSpecs_mod,  only : IXADV_O3
  use Config_module, only : dt_advec, KMAX_MID & 
                            ,PPBINV& ! 1.0e9, for conversion from mixing ratio to ppb
                            ,OutputVegO3
  use Debug_module,   only:  DEBUG   ! -> DEBUG%AOT
  use DO3SE_mod
  use GridValues_mod, only : debug_li, debug_lj, debug_proc, i_fdom, j_fdom
  use Io_Progs_mod,   only : datewrite
  use LandDefs_mod,   only : LandType, LandDefs
  use LocalVariables_mod, only : Grid, L
  use MetFields_mod, only: zen
  use NumberConstants, only : UNDEF_R, UNDEF_I
  use OwnDataTypes_mod, only: TXTLEN_DERIV, TXTLEN_SHORT, TXTLEN_IND,TXTLEN_NAME,&
                              O3cl_t
  use Par_mod, only : LIMAX, LJMAX, limax, ljmax, me
  use TimeDate_mod, only : current_date, print_date, jday => effectivdaynumber
  implicit none
  private

  public :: Calc_AOTx          ! called from AddMosaicOutput, My_DryDep
  public :: Calc_POD           ! called from AddMosaicOutput
  public :: Calc_GridAOTx      ! called from Derived_mod

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

   integer, save, public :: nOutputVegO3 = 0

    type(O3cl_t), public, allocatable, dimension(:) :: &
     VEGO3_OUTPUTS


contains
 !=========================================================================
 ! Calc_AOTx called from MosaicOutputs, at end of DryDep calculation.
 ! after deploss
  subroutine Calc_AOTx(iO3cl,iLC, aot )
    integer, intent(in) :: iO3cl,iLC
    real, intent(out)    :: aot

    real    :: o3
    integer :: i,j, mm,hh
    real :: X
    logical :: dbg, is_MM, is_EU,  inGS
    logical, save :: first_call = .true.
    character(len=*),parameter :: dtxt='CalcAOTx:'
    character(len=TXTLEN_NAME) :: txt

   ! MM (Mapping Manual) means use daylight O3, EU uses 7:00 -- 18:59 UTC

    is_EU =  VEGO3_OUTPUTS(iO3cl)%defn == 'EU'
    is_MM =  VEGO3_OUTPUTS(iO3cl)%defn(1:2) == 'MM'  ! eg MM, MM:test..
    mm = current_date%month
    hh = current_date%hour

    i = Grid%i
    j = Grid%j
    X = VEGO3_OUTPUTS(iO3cl)%Threshold
    dbg =  DEBUG%AOT .and. debug_proc .and. i == debug_li .and. j == debug_lj

    txt=VEGO3_OUTPUTS(iO3cl)%defn
    if( dbg ) then
      txt = dtxt// trim(VEGO3_OUTPUTS(iO3cl)%name)// print_date()
      write(*,"(a,4i4,3L2)") txt//' isF,C,iam? ', iLC, iO3cl, L%SGS, jday,  &
         LandType(iLC)%is_forest, LandType(iLC)%is_crop, LandType(iLC)%is_iam
    end if

    aot = 0.0

   ! Quick tests to see if we are in growing season and day period
   ! HAD BUG HERE for both MM-EMEP, since Apr-Spe and SGS/LGS both used.
   ! If night, or outside growing season, we simply exit with aot=0
   ! EU AOT for 8:00 -- 20:00 CET, is really 8:00 -- 19:59 CET 
   ! Or:        7:00 -- 18:59 UTC  (nb hour is integer value)

    if (  is_EU ) then

      inGS = .true. 
      if ( hh  <  7 .or. hh > 18 ) RETURN
      if( LandType(iLC)%is_forest .and. &
           ( mm<STARTMONTH_FOREST .or. mm>ENDMONTH_FOREST)   ) inGS = .false.
      if (   LandType(iLC)%is_crop .and. &
            ( mm<STARTMONTH_CROPS .or.mm>ENDMONTH_CROPS) ) inGS = .false.

    else if ( is_MM ) then

      if ( Grid%izen >= AOT_HORIZON ) RETURN !UN or MM use daylight
      inGS = ( jday >= L%SGS .and. jday <= L%EGS )

    else if( first_call .and.  me==0) then
      call StopAll(txt//"Unknown flags")
    end if

    if( dbg) write(*,"(a,4i4,L2)") txt//trim(LandDefs(iLC)%name)//&
        " jd,mm, SGS-EGS ", jday, mm, L%SGS,L%EGS, inGS

    
    if ( .not. inGS )  RETURN

  ! the wheat growing season is based upon the latitude function
    if ( is_MM ) then

      o3 = L%cano3_ppb   ! canopy top O3 for MM

    else if ( is_EU ) then

      o3 = Grid%surf_o3_ppb  ! Use 3m O3 for EU

    else 
      if ( first_call ) call CheckStop(txt//"AOT not MM or EU!")
    end if

   !TESTS ======== Plus and Minus 5 ppb used for hard-coded sens. tests
      if( VEGO3_OUTPUTS(iO3cl)%defn == 'MM:Plus5' ) then  ! nmole/m2/s
         o3 = o3 + 5
      else if( VEGO3_OUTPUTS(iO3cl)%defn == 'MM:Minus5' ) then  ! nmole/m2/s
         o3 = max( 0.0, o3 - 5 )
      end if
    !========== Calculate ========================

    if ( o3>X ) then  ! Add AOT, scaling for time-fraction

        aot = (o3-X)  ! dt_advec takes care of 3600s elsewhere

    end if

    if (  dbg ) call datewrite(dtxt//"AOTxdebug" // "defn:" // &
     trim(VEGO3_OUTPUTS(iO3cl)%defn),(/ iLC, Grid%izen /), (/ X,  o3, aot /))

    first_call = .false.

  end subroutine Calc_AOTx

 !=========================================================================
 ! Calc_GridAOTx called from Derived, after deploss

  function Calc_GridAOTx( iX ) result (aot)

    !/-- Calcuates AOT values for input threshold. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )
    !    Only relevant in ozone models, so far....

    integer, intent(in) :: iX  ! usually 40 (ppb)

    real, dimension(limax,ljmax) :: aot

    integer :: izen                    ! integer of zenith angle
    real :: o3, o3_ref                 ! Ozone (ppb) - needed if AOTs
    integer :: i, j 
    character(len=*), parameter :: dtxt='CalcGridAOT:'


    aot = 0.0

    !--------- ONLY April-Sept -------------------------
!CAREFUL!!!
!Keep for now. Irrelevant for global studies though
    if(  current_date%month<STARTMONTH_FOREST&
     .or.current_date%month>ENDMONTH_FOREST   )  RETURN
    !--------- ONLY April-Sept -------------------------

    do i=1, limax
        do j=1, ljmax

           izen = max(1,int( zen(i,j) + 0.5))

           if ( izen < AOT_HORIZON ) then                 ! Use 3m O3 for EU

                o3 = xn_adv(IXADV_O3,i,j,KMAX_MID) &
                     * cfac(IXADV_O3,i,j) * PPBINV

                aot(i,j) = max( o3 - iX , 0.0 )   ! Definition of AOTs

           end if
        end do
    end do
    if ( DEBUG%AOT .and. debug_proc ) then
       i=debug_li; j=debug_lj
       o3_ref = xn_adv(IXADV_O3,i,j,KMAX_MID) * PPBINV
       o3     = o3_ref * cfac(IXADV_O3,i,j)
       call datewrite(dtxt, (/ zen(i,j), o3_ref, o3, aot(i,j) /) )
    end if

  end function Calc_GridAOTx


!=========================================================================
!Calc_POD called from DryDep -> Add_MosaicOutput,  after deploss

  subroutine Calc_POD(iO3cl,iLC, pod, debug_flag, debug_txt )
    logical, intent(in) :: debug_flag
    character(len=*), intent(in), optional :: debug_txt
    integer, intent(in) :: iO3cl,iLC
    real, intent(out)   :: pod
    character(len=*),parameter :: dtxt='CalcPOD:'
    character(len=10):: txt

    real    :: o3, o3_3m
    integer :: i,j, spod, epod
    logical :: dbg
    real :: Y

    pod = 0.0
    dbg = ( DEBUG%AOT .and. debug_flag )

    if ( Grid%izen >= AOT_HORIZON ) then  !UN or MM use daylight
        if ( dbg ) write(*,"(2a,5i5)") dtxt//"Zen ",&
           trim(VEGO3_OUTPUTS(iO3cl)%name), iO3cl, jday, Grid%izen, L%SGS, L%EGS
        return
    end if

    ! Start and end of POD accumulation period:
    ! ( Some PODs are only for a specific period, different to 
    !   growing season (SGS), notably TC (wheat) which has 55 days.)

    if ( VEGO3_OUTPUTS(iO3cl)%RelSGS .eqv. .true. )  then
      spod = L%SGS + VEGO3_OUTPUTS(iO3cl)%SAccPeriod
      epod = L%SGS + VEGO3_OUTPUTS(iO3cl)%EAccPeriod
    else
      spod = L%SGS
      epod = L%EGS
    end if
    if ( dbg ) write(*,'(2a,5i5)') dtxt//"uZen ",&
     trim(VEGO3_OUTPUTS(iO3cl)%name)//'def:'//trim(VEGO3_OUTPUTS(iO3cl)%defn),&
      iO3cl, jday, Grid%izen, spod, epod
    if ( jday < spod .or. jday > epod ) then
       if ( dbg ) write(*,*) dtxt//"Jday RETURN "
       return
    end if

    i = Grid%i
    j = Grid%j
    Y = VEGO3_OUTPUTS(iO3cl)%Threshold ! nmole/m2/s

    o3 = L%cano3_ppb   ! canopy top O3 for MM
    o3_3m = Grid%surf_o3_ppb 

   ! Add fluxes if Y exceeded:

     pod  = max(L%FstO3 - Y,0.0)

    if ( dbg ) then
       write(txt,"(a,L1)") "Rel", VEGO3_OUTPUTS(iO3cl)%RelSGS
       !txt= "Rel"
       call datewrite(dtxt//"YYY:" // trim(txt) //  &
         trim(VEGO3_OUTPUTS(iO3cl)%defn), (/ iLC,   &
              VEGO3_OUTPUTS(iO3cl)%SAccPeriod,&
              VEGO3_OUTPUTS(iO3cl)%EAccPeriod, jday, Grid%izen /), &
           (/ Y,  o3, L%FstO3, pod /) )
       if(current_date%hour > 10 .and. pod > 0.1 ) then
       end if
    end if
!if(debug_flag) write(*,"(2a,4i5,2g12.3)") "PODO3 ", trim(VEGO3_OUTPUTS(iO3cl)%name), iO3cl, jday, spod, epod, L%FstO3,L%cano3_ppb

  end subroutine Calc_POD 

! !=========================================================================


! Will move SOMO and maxo3 here in future versions
!current definitions:
!SOMO35: daily max is found over 00:00 to 24:00. (not emepday)
!SOMO35: accumulated over one year
!D2_MAXO3 :  daily max is found over an EMEPDAY
!D2_MAXO3 : accumulated into yearly_output from April to September


end module AOTx_mod
