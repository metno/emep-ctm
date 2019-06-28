! <CoDep_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module CoDep_mod
  use CheckStop_mod,  only : CheckStop
  use Chemfields_mod, only : so2nh3_24hr
  use GridValues_mod, only : debug_proc, debug_li, debug_lj
  use Config_module, only : MasterProc  ! for DEBUG
  use Par_mod,        only : LIMAX, LJMAX



  !---------------------------------------------------------------------------
  ! Calculates the  non-stomatal resistances Rns_SO2 and Rns_NH3. 
  !---------------------------------------------------------------------------
  ! For basic reference and methods, see
  !
  ! Simpson et al. 2012, ACP, 12, for general model documentation
  !
  ! Also,
  ! RIS: Smith, R.I., Fowler, D., Sutton, M.A., Flechard, C: and Coyle, M.
  !      Atmos. Environ., 34, 3757-3777
  ! + Errata + pers. comm. with R.I.Smith
  !
  ! and
  ! Smith,  Unpublished Note
  ! Eiko Nemitz papers
  !
  !---------------------------------------------------------------------------

  implicit none

  public   :: make_so2nh3_24hr ! make 24hr ratio
  public   :: CoDep_factors
  private  :: Tabulate        ! pre-calculate many values to save CPU
 
  !*** Some parameters for the so2nh3_24hr calculations

  integer, private :: nhour 
  real, private, save,allocatable ::   &             ! 24hr average ratio              
     so2nh3_hr(:,:,:)  !Predefined to 1.0 to make first hours
                                         !reasonable

  !*** Some parameters for the Rns calculations

  integer, private, parameter :: TMIN = -40, TMAX = 40 ! Allowed temp. range
  integer, private, parameter :: NTAB  = 100   ! No. intervals used for tabulation

   real, private, save, dimension(0:100)   :: tab_exp_rh  ! For eqn (8.16)
   real, private, save, dimension(0:NTAB) :: &
           tab_acidity_fac,                &
           tab_F2,&                            ! For emepctm eqn (8.16)
           tab_F4, &                           ! for Rns_NH3
           tab_F3                              !For Rns SO2 

  !/ Calculated values /outputs):
   real, public, save ::  &
       humidity_fac      &! to interpolate Gns  across different RH
       ,Rns_NH3           & ! Resistance for NH3 on ground (water) surface
       ,Rns_SO2             ! Resistance for SO2 on ground (water) surface

   logical, private, parameter :: MY_DEBUG = .false.
   real, private, parameter  :: MAX_SN = 3.0 !max SO2/NH3 ratio


contains
! =======================================================================

  subroutine CoDep_factors( so2nh3ratio24hr,so2nh3ratio, Ts_C, frh, &
                            forest, debug_flag)
! =======================================================================
!
! =======================================================================


! Input:
   
   real, intent(in) :: so2nh3ratio    ! so2/nh3 ratio
   real, intent(in) :: so2nh3ratio24hr    ! so2/nh3 ratio accumulated 
                                          ! over previous 24 hours

   real, intent(in) :: Ts_C           ! surface temp. (degrees C)
   real, intent (in):: frh            ! relative humidity (as a fraction)
   logical, intent (in):: forest      ! true if forest
   logical, intent (in):: debug_flag


 ! On first call we will tabulate Rns_NH3

   logical, save :: my_first_call = .true.

  !local terms governing intermediary calculations in the evaluation of NH3_st:

   real, parameter :: BETA=1.0/22.0   ! Rns factors, see emepctm eqn (8.16)
   real    :: F1, F2           ! Rns factors, NH3
   real    :: F3, F4           ! Rns factors for SO2
   real    :: a_SN             ! so2/nh3 after correction with 0.6
   real    :: a_SN_24hr        ! so2/nh3 24hr average 
   integer :: itemp            ! integer Temp in deg.C
   integer :: ia_SN            ! 10*a_SN
   integer :: ia_SN_24hr       ! 10*a_SN_24hr

   integer :: IRH              ! RH as percent
 
   if ( my_first_call ) then

     call Tabulate()
     my_first_call = .false.
     allocate(so2nh3_hr(24,LIMAX,LJMAX))
     so2nh3_hr=1.0
     if( debug_proc ) write(*,"(a,2es12.4,f7.2,f7.3,L1)") "First CoDep call, ",  &
           so2nh3ratio24hr,so2nh3ratio, Ts_C, frh, forest

   end if
       
   itemp     =  int( Ts_C + 0.5 )
   itemp     =  max(itemp, TMIN)   ! For safety
   itemp     =  min(itemp, TMAX)   ! For safety


    a_SN  = min(MAX_SN,so2nh3ratio)

    ia_SN = nint( NTAB * a_SN/MAX_SN )   ! Spread values from 0-3 to 0:100


   ! Cap a_SN_24hr at 3
   a_SN_24hr  = min(MAX_SN,so2nh3ratio24hr) 
   ia_SN_24hr = nint( NTAB * a_SN_24hr/MAX_SN )


   IRH   = max( 1,  int( 100.0 * frh ) )
   if ( MY_DEBUG ) then
      if ( IRH<1 .or. IRH>100 .or. ia_SN < 0 ) then 
       print *, "CODEP ERROR ", IRH, frh, ia_SN, a_SN
       call CheckStop ( IRH<1 .or. IRH>100  , "CoDep IRH ERROR")
      end if
   end if


  !/ 3) Rns_NH3   - see emepctm eqn (8.16)
  !     Rns_SO2   - Fagerli et al, in preperation


    if (Ts_C >0 ) then    ! Use "rh" - now in fraction 0..1.0

          !F1 = 10.0 * log10(Ts_C+2.0) * exp(100.0*(1.0-frh)/7.0)
           F1 = 10.0 * log10(Ts_C+2.0) * tab_exp_rh(IRH)
           F2 = tab_F2( ia_SN  )

           Rns_NH3 = BETA * F1 * F2
           Rns_NH3 = min( 200.0, Rns_NH3)  ! After discussion with Ron
           Rns_NH3 = max(  10.0,Rns_NH3)

        ! New Formulation (article to be submitted)
        ! Rns_SO2_dry = 11.84  * exp(1.1*so2nh3ratio24hr) * ( frh**(-1.67) )

           F3 = tab_F3 (ia_SN_24hr) !11.84  * exp(1.1*so2nh3ratio24hr)
           F4 = tab_F4(IRH)      !frh**(-1.67)
           Rns_SO2 = F3 * F4


        !Limit on So2/nh3 

          Rns_SO2 = min( 1000.0, Rns_SO2)  ! Set because rh goes almost 
                 ! to zero occationally over the Alps, Greenland and Svalbard
                 ! 1000 chosen sort of random


           Rns_SO2 = max(  10.0,Rns_SO2) !hf CoDep SHOULD WE LIMIT IT to 10??


       if(MY_DEBUG .and. debug_flag) then
         write(*,*) "CODEP PRE rh", IRH, frh
         write(*,*) "CODEP PRE NH3", ia_SN, a_SN, F1, F2, Rns_NH3
         write(*,*) "CODEP PRE SO2",ia_SN_24hr, a_SN_24hr, F3, F4, Rns_SO2
       end if


     else if (  Ts_C > -5 ) then
       ! In a future version, we might test for  Ts_C > -2 instead of
       ! Ts_C > 0  as we have now.

           Rns_NH3 = 100.0 ! will be modified bu low-T factor later
           Rns_SO2 = 100.0

     else 

           Rns_NH3= 500.0 ! will be modified bu low-T factor later
           Rns_SO2= 500.0

     end if !Ts_C

 end subroutine CoDep_factors

  !=======================================================================

   subroutine Tabulate()
    !***  Tabulates humidity factors, 

     real :: a_SN, a_SN_24hr
     integer :: IRH, ia_SN, ia_SN_24hr

    ! Acidity factor

     do ia_SN = 0, NTAB
       a_SN =  ia_SN * MAX_SN /real(NTAB)
       tab_acidity_fac( ia_SN )  = exp( -(2.0- a_SN) )
       tab_F2 (ia_SN)            = 10.0**( (-1.1099 * a_SN)+1.6769 )
       if(MY_DEBUG.and. MasterProc ) write(6,*) "TABIA ", ia_SN, a_SN, &
              tab_acidity_fac( ia_SN ), tab_F2(ia_SN)
     end do

     do ia_SN_24hr = 0, NTAB
       a_SN_24hr            = ia_SN_24hr * MAX_SN /real(NTAB)
       tab_F3 (ia_SN_24hr)  = 11.84  * exp(1.1 * a_SN_24hr)
       if(MY_DEBUG.and. MasterProc ) write(6,*) "TABIA24 ", ia_SN_24hr, &
            a_SN_24hr, tab_F3(ia_SN_24hr)
     end do

     do IRH = 0, 100
          tab_exp_rh(IRH) = exp( (100.0-IRH)/7.0)
          tab_F4(IRH)= (max(1,IRH)/100.)**(-1.67)
          if(MY_DEBUG.and. MasterProc ) write(6,*) "TABRH ", IRH, &
              tab_exp_rh(IRH), tab_F4(IRH)
     end do

   end subroutine Tabulate
  !=======================================================================

!=======================================================================
  subroutine make_so2nh3_24hr(hour,so2conc,nh3conc,cfac_so2,cfac_nh3)
! Calculates so2/nh3 ratio for surface, accumulated over the 
! previous 24 hours (called every hour in phychem)
! Could call more often, but then we need to keep more arrays in memory.

  implicit none

  integer, intent (in) :: hour
  integer :: i,j

  real, dimension(LIMAX,LJMAX), intent(in) :: so2conc
  real, dimension(LIMAX,LJMAX), intent(in) :: nh3conc
  real, dimension(LIMAX,LJMAX), intent(in) :: cfac_so2
  real, dimension(LIMAX,LJMAX), intent(in) :: cfac_nh3
  logical :: debug_flag        ! set true when i,j match DEBUG_i, DEBUG_j

  nhour=hour
  if (hour == 0) nhour=24
  so2nh3_hr(nhour,:,:)=so2conc*cfac_so2/max(1.0e-15,(nh3conc*cfac_nh3))
 !so2conc in mixing ratio 

  so2nh3_24hr(:,:)=0.0 !initialize each time step (hour)

  do i=1, LIMAX
     do j=1, LJMAX
        if( MY_DEBUG ) debug_flag =  &
           ( debug_proc .and. i == debug_li .and. j == debug_lj)
        do nhour=1,24
           so2nh3_24hr(i,j)=so2nh3_24hr(i,j)+so2nh3_hr(nhour,i,j)/24. 
           if( MY_DEBUG .and. debug_flag ) then
              write(*,*) "so2nh3_hr(nhour,i,j)",nhour,&
                so2nh3_hr(nhour,i,j),so2nh3_24hr(i,j)
              write(*,*) "so2nh3_24hr output", so2nh3_24hr(i,j),&
                so2conc(i,j),cfac_so2(i,j),nh3conc(i,j),cfac_nh3(i,j)
          end if
        end do ! nhour
     end do 
  end do

end subroutine make_so2nh3_24hr
!=======================================================================


end module CoDep_mod


