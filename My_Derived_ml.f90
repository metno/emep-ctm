! <My_Derived_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

!==============================================================================
module My_Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module specifies the "derived" fields, such as accumulated 
  ! precipitation 
  ! or sulphate, daily, monthly or yearly averages, depositions. These fields
  ! are all typically output as netCDF fields.
  !
  ! This module provides the user-defined setups which are used in Derived_ml.
  ! Derived fields are identified by a "class", such as "ADV" of "VOC", and
  ! the Derived_ml should perform any integrations for this.
  ! 
  ! Several often-used routines (e.g. for AOTs, acc. sulphate, are defined 
  ! in the Derived_ml.f90, but users can define their own here, since
  ! we do not use "use only" in Derived_ml. 
  !  
  !   Only text strings used here to define wanted data
  !   All data field characteristics should be defined in Derived_ml, e.g.
  !   in f_2d arrays. 
  !   Derived fields such as d_2d only exist in Derived_ml, so are
  !   accessed here through subroutine calls - using just the (i,j) part
  !   of the bigger d_2d arrays
  !---------------------------------------------------------------------------
 
use GenSpec_adv_ml        ! Use IXADV_ indices...
use GenSpec_shl_ml        ! Use IXSHL_ indices...
use GenSpec_tot_ml,  only : SO4, HCHO, CH3CHO  &   !  For mol. wts.
                           ,aNO3, pNO3, aNH4, PM25, PMCO &
                           ,SSfi, SSco  !SeaS
use GenChemicals_ml, only : species               !  For mol. wts.
use ModelConstants_ml, only : atwS, atwN, ATWAIR  &
                        , SOURCE_RECEPTOR  &  
                        , KMAX_MID & ! =>  z dimension
                        , PPBINV  &  !   1.0e9
                        , MFAC       ! converts roa (kg/m3 to M, molec/cm3)

use Chemfields_ml, only : xn_adv, xn_shl, cfac
use GenSpec_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use Met_ml,        only : z_bnd, roa    ! 6c REM: zeta
use Par_ml,    only: me, MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     limax, ljmax           ! => used x, y area 
use SmallUtils_ml,  only : AddArray, LenArray, NOT_SET_STRING, WriteArray
use TimeDate_ml,   only : current_date
implicit none
private

 public  :: Init_My_Deriv
 public  :: My_DerivFunc 

 private :: misc_xn   &          ! Miscelleaneous Sums and fractions of xn_adv
           ,pm_calc              ! Miscelleaneous PM's


   character(len=8),  public ,parameter :: model='ZD_OZONE'


   !/** Depositions are stored in separate arrays for now - to keep size of
   !    derived arrays smaller and to allow possible move to a Deposition
   !    module at a later stage.
   !  Factor 1.0e6 converts from kg/m2/a to mg/m2/a

  !        We normally distinguish source-receptor (SR) stuff from model
  !        evaluation.  The SR runs should use as few as possible outputs
  !        to keep CPU and disc-requirements down. We define first then the
  !        minimum list of outputs for use in SR, then define an extra list
  !        of parameters needed in model evaluation, or even for the base-case
  !        of SR runs.


  !============ parameters for source-receptor modelling: ===================!

    integer, public, parameter :: MAX_NUM_DERIV2D = 200
    integer, public, parameter :: MAX_NUM_DERIV3D =   5
    character(len=12), public, save, &
         dimension(MAX_NUM_DERIV2D) :: wanted_deriv2d = NOT_SET_STRING
    character(len=12), public, save, &
         dimension(MAX_NUM_DERIV3D) ::  wanted_deriv3d = NOT_SET_STRING

    integer, private, save :: mynum_deriv2d
    integer, private, save :: mynum_deriv3d

 

    character(len=12), public, parameter, dimension(49) :: &
  D2_SR = (/ &
!
!    Particles: components
       "D2_SO4      ","D2_aNO3     ","D2_pNO3     ","D2_aNH4     " &
      ,"D2_PPM25    ","D2_PPMco    ","D2_PM25_H2O " &
!
!    Particles: sums
      ,"D2_SIA      ","D2_PM25     ","D2_PM10     ","D2_PMco     " &
      ,"D2_SS       ","D2_tNO3     " &
!
!    Ozone and AOTs
      ,"D2_O3       ","D2_MAXO3    " &
      ,"D2_AOT30    ","D2_AOT40    ","D2_AOT60    " &
      ,"D2_AOT30f   ","D2_AOT40f   ","D2_AOT60f   ","D2_AOT40c   " &
      ,"D2_EUAOT30WH","D2_EUAOT30DF","D2_EUAOT40WH","D2_EUAOT40DF" &
      ,"D2_UNAOT30WH","D2_UNAOT30DF","D2_UNAOT40WH","D2_UNAOT40DF" &
      ,"D2_MMAOT30WH","D2_MMAOT40WH" &
      ,"D2_SOMO35   ","D2_SOMO0    " &
!
!    NOy-type sums 
      ,"D2_NO2      ","D2_OXN      ","D2_NOX      ","D2_NOZ      " &
      ,"D2_OX       "  &
!
!    Ecosystem - fluxes: 
      ,"D2_AFSTDF0  ","D2_AFSTDF16 ","D2_AFSTBF0  ","D2_AFSTBF16 " &
      ,"D2_AFSTCR0  ","D2_AFSTCR3  ","D2_AFSTCR6  " & !
       ,"D2_O3DF     ","D2_O3WH     " &
!
!    Surface  pressure (for cross section):
      ,"PS          " &
  /)

  !============ Extra parameters for model evaluation: ===================!

    character(len=12), public, parameter, dimension(7) :: &
  D2_EXTRA = (/ &
       "D2_SO2      ","D2_HNO3     ","D2_NH3      ","D2_VOC      "&
      ,"D2_REDN     ","D2_SSfi     ","D2_SSco     " &
  /)


!----------------------
! Less often needed:
 !exc  "D2_CO     ","D2T_HCHO  ","D2T_CH3CHO","D2_VOC    ",
 !exc ,"D2_O3CF   ","D2_O3TC   ","D2_O3GR   ","D2_ACCSU  ",
 !"D2_FRNIT  ","D2_MAXOH  ","D2_HMIX   ","D2_HMIX00 ","D2_HMIX12 " &
 !exc  "D2_PAN    ","D2_AOT20    " /)

   !======= MY_DERIVED SYSTEM ======================================

  ! use character arrays to specify which outputs are wanted

   character(len=9), public, parameter, dimension(4) :: &
       WDEP_WANTED = (/ "WDEP_PREC", "WDEP_SOX ", "WDEP_OXN ", &
                      "WDEP_RDN " /)   ! WDEP_PM not used

  !( waters and wetlands removed:)

   character(len=10), public, parameter, dimension(15) :: &
     DDEP_WANTED = (/  &
        "DDEP_SOX  ","DDEP_OXN  ","DDEP_RDN  "  &
       ,"DDEP_OXSCF","DDEP_OXSDF","DDEP_OXSCR","DDEP_OXSSN"  &
       ,"DDEP_OXNCF","DDEP_OXNDF","DDEP_OXNCR","DDEP_OXNSN"  &
       ,"DDEP_RDNCF","DDEP_RDNDF","DDEP_RDNCR","DDEP_RDNSN"  &
     /)

     character(len=13), public, parameter, dimension(2) :: &
       D3_WANTED = (/ "D3_O3        ","D3_TH        " /)


    integer, private :: i,j,k,n, ivoc, index    ! Local loop variables

   contains

 !=========================================================================
  subroutine Init_My_Deriv()

   ! Build up the array wanted_deriv2d with the required field names

     call AddArray(WDEP_WANTED, wanted_deriv2d, NOT_SET_STRING)
     call AddArray(DDEP_WANTED, wanted_deriv2d, NOT_SET_STRING)
     call AddArray( D2_SR,  wanted_deriv2d, NOT_SET_STRING)

     if ( .not. SOURCE_RECEPTOR ) then !may want extra?
        call AddArray( D2_EXTRA, wanted_deriv2d, NOT_SET_STRING)
     end if
     mynum_deriv2d  = LenArray( wanted_deriv2d, NOT_SET_STRING )

   ! ditto wanted_deriv3d....

     !if ( .not. SOURCE_RECEPTOR ) then
     !   call AddArray( D3_WANTED,  wanted_deriv3d, NOT_SET_STRING)
     !end if
     mynum_deriv3d  = LenArray( wanted_deriv3d, NOT_SET_STRING )


     if ( me == 0 ) then
       write(*,*) "Init_My_Deriv, mynum_deriv2d = ", mynum_deriv2d
       call WriteArray(wanted_deriv2d,mynum_deriv2d," Wanted 2d array is")
       write(*,*) "Init_My_Deriv, mynum_deriv3d = ", mynum_deriv3d
       call WriteArray(wanted_deriv3d,mynum_deriv3d," Wanted 3d array is")
     
     end if

  end subroutine Init_My_Deriv
 !=========================================================================
  subroutine My_DerivFunc( e_2d, n, class , timefrac, density )

    ! We define here here any functions which cannot easily be defined
    ! in the more general Derived_ml. For example, we need the 
    ! index for IXADV_O3 for AOTs, and this might not be available in the model
    ! we are running (a PM2.5 model for example), so it is better to define 
    ! this function here.

  real, dimension(:,:), intent(inout) :: e_2d  !  (i,j) 2-d extract of d_2d
  integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
  character(len=*), intent(in)    :: class       ! Class of data
  real, intent(in)    :: timefrac    ! Timestep as frationof hour, dt/3600

  real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density     
! density = 1 ( or = roa when unit ug)

  select case ( class )

      case ( "OX", "NOX", "NOZ", "TOXN", "TRDN", "FRNIT", "tNO3 ", "SSalt" )

           call misc_xn( e_2d, n, class, density )

      case ( "SIA", "PM10", "PM25", "PMco" )

          call pm_calc(e_2d, n, class,  density)

      case  default

            print *, "WARNING - REQUEST FOR UNDEFINED OUTPUT:", n, class
     end select
  


  end subroutine My_DerivFunc
 !=========================================================================

  subroutine pm_calc( pm_2d, n, class, density )

    !/--  calulates PM10 = SIA + PPM2.5 + PPMco

    real, dimension(:,:), intent(inout) :: pm_2d  ! i,j section of d_2d arrays
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    character(len=*)    :: class   ! Type of data
    real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density  

    select case ( class )

    case ( "SIA" ) 

      forall ( i=1:limax, j=1:ljmax )
        !ds d_2d( n, i,j,IOU_INST) = &
        pm_2d( i,j) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt *cfac(IXADV_SO4,i,j)  &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j))& 
         * density(i,j)
      end forall

    case ( "PM25" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt *cfac(IXADV_SO4,i,j)  &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j) &
         + xn_adv(IXADV_PM25,i,j,KMAX_MID)*species(PM25)%molwt*cfac(IXADV_PM25,i,j) & 
         + xn_adv(IXADV_SSfi,i,j,KMAX_MID)*species(SSfi)%molwt *cfac(IXADV_SSfi,i,j))&  !SeaS
         * density(i,j)
      end forall

    case ( "PMco" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_PMco,i,j,KMAX_MID)*species(PMCO)%molwt*cfac(IXADV_PMco,i,j) & 
         + xn_adv(IXADV_SSco,i,j,KMAX_MID) *species(SSco)%molwt *cfac(IXADV_SSco,i,j))&  !SeaS
          * density(i,j)
      end forall

    case ( "PM10" ) 

      forall ( i=1:limax, j=1:ljmax )
        pm_2d( i,j ) = &
         ( xn_adv(IXADV_SO4,i,j,KMAX_MID) *species(SO4)%molwt*cfac(IXADV_SO4,i,j)   &
         + xn_adv(IXADV_aNO3,i,j,KMAX_MID)*species(aNO3)%molwt*cfac(IXADV_aNO3,i,j) &
         + xn_adv(IXADV_pNO3,i,j,KMAX_MID)*species(pNO3)%molwt*cfac(IXADV_pNO3,i,j) &
         + xn_adv(IXADV_aNH4,i,j,KMAX_MID)*species(aNH4)%molwt*cfac(IXADV_aNH4,i,j) &
         + xn_adv(IXADV_PM25,i,j,KMAX_MID)*species(PM25)%molwt*cfac(IXADV_PM25,i,j) &
         + xn_adv(IXADV_PMco,i,j,KMAX_MID)*species(PMCO)%molwt*cfac(IXADV_PMco,i,j) & 
         + xn_adv(IXADV_SSfi,i,j,KMAX_MID)*species(SSfi)%molwt*cfac(IXADV_SSfi,i,j) & !SeaS
         + xn_adv(IXADV_SSco,i,j,KMAX_MID)*species(SSco)%molwt*cfac(IXADV_SSco,i,j))& !SeaS
         * density(i,j)
      end forall

    end select

  end subroutine  pm_calc
 !=========================================================================

!=========================================================================

  subroutine misc_xn( e_2d, n, class, density)
    real, dimension(:,:), intent(inout) :: e_2d  ! i,j section of d_2d arrays
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    character(len=*)    :: class   ! Type of data
    real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density  
! density = 1 ( or = roa when unit ug)


    !/--  adds up sulphate, nitrate, or whatever is defined

    select case ( class )

    case ( "TOXN" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_HNO3,i,j,KMAX_MID) * cfac(IXADV_HNO3,i,j) &
              + xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j)) &
              * density(i,j)
      end forall


! OX for O3 and NO2 trend studies

    case ( "OX" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
                xn_adv(IXADV_O3,i,j,KMAX_MID)  * cfac(IXADV_O3,i,j)   &
              + xn_adv(IXADV_NO2,i,j,KMAX_MID) * cfac(IXADV_NO2,i,j) 
      end forall

    case ( "NOX" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_NO,i,j,KMAX_MID) &
              + xn_adv(IXADV_NO2,i,j,KMAX_MID) * cfac(IXADV_NO2,i,j) &
              ) * density(i,j)
      end forall

    case ( "NOZ" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_HNO3,i,j,KMAX_MID) * cfac(IXADV_HNO3,i,j) &
              + xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j) &
              + xn_adv(IXADV_PAN,i,j,KMAX_MID) * cfac(IXADV_PAN,i,j) &
              + xn_adv(IXADV_MPAN,i,j,KMAX_MID) * cfac(IXADV_MPAN,i,j) &
              + xn_adv(IXADV_NO3,i,j,KMAX_MID) &
              + 2.0* xn_adv(IXADV_N2O5,i,j,KMAX_MID) &
              + xn_adv(IXADV_ISNI,i,j,KMAX_MID) &
              ) * density(i,j)
      end forall


    case ( "TRDN" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
               ( xn_adv(IXADV_NH3,i,j,KMAX_MID) * cfac(IXADV_NH3,i,j)    &
              +  xn_adv(IXADV_aNH4,i,j,KMAX_MID) * cfac(IXADV_aNH4,i,j))  &
               * density(i,j)
      end forall


    case ( "FRNIT" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
             ( xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j)  &
            +  xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j)) &
       /max(1E-80, (xn_adv(IXADV_HNO3,i,j,KMAX_MID) *  cfac(IXADV_HNO3,i,j))&
            +  xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j)    &
            +  xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j))
      end forall

    case ( "tNO3" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d(  i,j ) = &
              ( xn_adv(IXADV_aNO3,i,j,KMAX_MID) * cfac(IXADV_aNO3,i,j) &
              + xn_adv(IXADV_pNO3,i,j,KMAX_MID) * cfac(IXADV_pNO3,i,j) )&
              * density(i,j)
      end forall

    case ( "SSalt" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_SSfi,i,j,KMAX_MID) * cfac(IXADV_SSfi,i,j) &
              + xn_adv(IXADV_SSco,i,j,KMAX_MID) * cfac(IXADV_SSco,i,j) )&
              * density(i,j)
      end forall

    end select
  end subroutine misc_xn
 !=========================================================================
end module My_Derived_ml
