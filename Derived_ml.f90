! <Derived_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module performs the calculations associated with "derived" 2D and 3D,
  ! such as accumulated precipitation or sulphate, daily, monthly or yearly 
  ! averages, depositions. These fields are all typically output as netCDF 
  ! fields.
  !
  ! This routine defines many possible derived  outputs. 
  ! The names of the derived fields actualy required should have been specified
  !  in the user-defined My_Derived_ml.
  !

  ! User-defined routines and treatments are often needed here. Here there is 
  ! added stuff for VOC, AOTs, accsu. In 
  ! general such code should be added in such a way that it isn't activated if
  ! not needed. It then doesn't need to be commented out if not used.
  !---------------------------------------------------------------------------

!current definitions:
!SOMO35: daily max is found over 00:00 to 24:00. (not emepday)
!SOMO35: accumulated over one year
!D2_MAXO3 :  daily max is found over an EMEPDAY
!D2_MAXO3 : accumulated into yearly_output from April to September
!AOTXXc: accumulated into yearly_output from May to July
!AOTXXf: accumulated into yearly_output from April to September
!D2_EUAOTXXWH: accumulated into yearly_output from May to July
!D2_EUAOTXXDF: accumulated into yearly_output from April to September
!D2_UNAOTXXWH: accumulated into yearly_output from May to July
!D2_UNAOTXXDF: accumulated into yearly_output from April to September
!D2_MMAOTXXWH: accumulated into yearly_output over growing season
!D2_O3 is now yearly accumulated

use My_Derived_ml, only : &
            wanted_deriv2d, wanted_deriv3d  &! names of wanted derived fields
           ,Init_My_Deriv, My_DerivFunc

use CheckStop_ml,      only: CheckStop
use Chemfields_ml, only : xn_adv, xn_shl, cfac,xn_bgn, PM_water
use GenSpec_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use GenSpec_shl_ml
use GenSpec_tot_ml
use GenChemicals_ml, only : species
use GridValues_ml, only : debug_li, debug_lj, debug_proc
use Met_ml, only :   roa,pzpbl,xksig,ps,th,zen
use ModelConstants_ml, &
                   only: KMAX_MID &   ! =>  z dimension
                        , NPROC   &   ! No. processors
                        , atwS, atwN, ATWAIR  &
                        , PPBINV  &   !   1.0e9, for conversion of units 
                        , PPTINV  &   !   1.0e12, for conversion of units 
                        , MFAC    &   ! converts roa (kg/m3 to M, molec/cm3)
                        , AOT_HORIZON&! limit of daylight for AOT calcs
                        ,DEBUG_i, DEBUG_j & 
                        , SOURCE_RECEPTOR &
                        , NTDAY        !Number of 2D O3 to be saved each day (for SOMO)  
use Par_ml,    only: MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     me,                &   ! for print outs
                     gi0,gj0,IRUNBEG,JRUNBEG,&! for i_fdom, j_fdom
                     li0,lj0,limax, ljmax    ! => used x, y area 
use PhysicalConstants_ml,  only : PI
use SmallUtils_ml, only: find_index, LenArray, NOT_SET_STRING
use TimeDate_ml, only : day_of_year,daynumber,current_date

implicit none
private

 public  :: Init_Derived         !
 public  :: ResetDerived         ! Resets values to zero
 public  :: DerivedProds         ! Calculates any production terms
 private :: AddDef 
 private :: Define_Derived       !
 private :: Setups 

 public :: Derived              ! Calculations of sums, avgs etc.
 private :: Setup_VOC            ! Defines VOC group
 private :: voc_2dcalc           ! Calculates sum of VOC for 2d fields
 private :: voc_3dcalc           ! Calculates sum of VOC for 3d fields


   INCLUDE 'mpif.h'
   INTEGER STATUS(MPI_STATUS_SIZE),INFO

    type, public:: Deriv
       character(len=9) :: class ! Type of data, e.g. ADV or VOC
       logical  :: avg      ! True => average data (divide by nav at end),
                            !     else accumulate over run period
       integer  :: index    ! index in concentation array, or other
       real     :: scale    ! Scaling factor
       logical  :: rho      ! True when scale is ug (N or S)
       logical  :: inst     ! True when instantaneous values needed
       logical  :: year     ! True when yearly averages wanted
       logical  :: month    ! True when monthly averages wanted
       logical  :: day      ! True when daily averages wanted
       character(len=15) :: name ! Name of the variable (for netCDF output)
       character(len=10) :: unit ! Unit (writen in netCDF output)
    end type Deriv

   logical, private, parameter :: T = .true., F = .false. ! shorthands only

  ! Tip. For unix users, do "grep AddDef | grep -v Is3D | wc" or similar
  ! to help get the number of these:
   integer, private, parameter ::  &
       MAXDEF_DERIV2D =100 & ! Max. No. 2D derived fields to be defined
      ,MAXDEF_DERIV3D = 17   ! Max. No. 3D derived fields to be defined

   integer, public, save :: num_deriv2d, num_deriv3d
   integer, private, save :: Nadded2d = 0, Nadded3d=0 ! No. defined derived 

  ! We put definitions of **all** possible variables in def_2d, def_3d
  ! and copy the needed ones into f_xx. The data will go into d_2d, d_3d

    type(Deriv),private, dimension(MAXDEF_DERIV2D), save :: def_2d
    type(Deriv),private, dimension(MAXDEF_DERIV3D), save :: def_3d

    type(Deriv),public, allocatable, dimension(:), save :: f_2d
    type(Deriv),public, allocatable, dimension(:), save :: f_3d


  ! Define 4 output types corresponding to instantaneous,year,month,day

   integer, public, parameter ::  &
        IOU_INST=1, IOU_YEAR=2, IOU_MON=3, IOU_DAY=4, IOU_HOUR=5

  ! The 2-d and 3-d fields use the above as a time-dimension. We define
  ! LENOUTxD according to how fine resolution we want on output. For 2d
  ! fields we use daily outputs. For the big 3d fields, monthly output
  ! is sufficient.

   integer, public, parameter ::  LENOUT2D = 4  ! Allows INST..DAY for 2d fields
   integer, public, parameter ::  LENOUT3D = 4  ! Allows INST..MON for 3d fields

  !e.g. d_2d( num_deriv2d,MAXLIMAX, MAXLJMAX, LENOUT2D)
  ! &   d_3d( num_deriv3d,MAXLIMAX, MAXLJMAX, KMAX_MID, LENOUT3D )
   real, save, public, allocatable, dimension(:,:,:,:) :: d_2d
   real, save, public, allocatable, dimension(:,:,:,:,:) :: d_3d


   ! save O3 every hour during one day to find running max
    real, save,  public :: &     ! to be used for SOMO35
     D2_O3_DAY( MAXLIMAX, MAXLJMAX, NTDAY) = 0. 


  ! Counters to keep track of averaging
  ! Initialise to zero in Init.

    integer, public, allocatable, dimension(:,:), save :: nav_2d
    integer, public, allocatable, dimension(:,:), save :: nav_3d

   ! Note - previous versions did not have the LENOUT2D dimension
   ! for wet and dry deposition. Why not?  Are annual or daily
   ! depositions never printed? Since I prefer to keep all 2d
   ! fields as similar as posisble, I have kept this dimension
   ! for now - ds


   !-- some variables for the VOC sum done for ozone models
   !   (have no effect in non-ozone models - leave in code)

   integer, private, save :: nvoc   ! No. VOCs 
   integer, private, dimension(NSPEC_ADV), save :: &
             voc_index, &     ! Index of VOC in xn_adv
             voc_carbon       ! Number of C atoms

   logical, private, parameter :: MY_DEBUG = .false.
   logical, private, save :: debug_flag, Is3D
   character(len=100), private :: errmsg

   integer, private :: i,j,k,n, ivoc, index    ! Local loop variables
   integer, public, parameter:: startmonth_forest=4,endmonth_forest=9&
                                ,startmonth_crops=5,endmonth_crops=7

   contains

    !=========================================================================
    subroutine Init_Derived()

        integer :: alloc_err
          if(me==0 .and. MY_DEBUG) write(*,*) "INITIALISE My DERIVED STUFF"
          call Init_My_Deriv()  !-> wanted_deriv2d, wanted_deriv3d

         ! get lengths of wanted arrays (excludes notset values)
          num_deriv2d = LenArray(wanted_deriv2d,NOT_SET_STRING)
          num_deriv3d = LenArray(wanted_deriv3d,NOT_SET_STRING)
 
     if ( num_deriv2d > 0 ) then
          if(me==0 .and. MY_DEBUG) write(*,*) "Allocate arrays for 2d: ", num_deriv2d
          allocate(f_2d(num_deriv2d),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of f_2d")
          allocate(d_2d(num_deriv2d,MAXLIMAX,MAXLJMAX,LENOUT2D),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of d_2d")
          call CheckStop(alloc_err,"Allocation of d_3d")
          allocate(nav_2d(num_deriv2d,LENOUT2D),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of nav_2d")
          nav_2d = 0
     end if
     if ( num_deriv3d > 0 ) then
          if(me==0 .and. MY_DEBUG) write(*,*) "Allocate arrays for 3d: ", num_deriv3d
          allocate(f_3d(num_deriv3d),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of f_3d")
          allocate(d_3d(num_deriv3d,MAXLIMAX,MAXLJMAX,KMAX_MID,LENOUT3D),&
                 stat=alloc_err)
          allocate(nav_3d(num_deriv3d,LENOUT3D),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of nav_3d")
          nav_3d = 0
     end if

          call Define_Derived()
          call Setups()

    end subroutine Init_Derived

   !=========================================================================
    subroutine AddDef(class,avg,index,scale,rho,inst,year,month,day,&
           name,unit,Is3D)

       character(len=*), intent(in) :: class ! Type of data, e.g. ADV or VOC
       logical, intent(in)  :: avg      ! True => average data (divide by 
                                        ! nav at end), else accumulate over 
                                        ! run period
       integer, intent(in)  :: index    ! index in e.g. concentration array
       real, intent(in)     :: scale    ! Scaling factor
       logical, intent(in)  :: rho      ! True when scale is ug (N or S)
       logical, intent(in)  :: inst     ! True when instantaneous values needed
       logical, intent(in)  :: year     ! True when yearly averages wanted
       logical, intent(in)  :: month    ! True when monthly averages wanted
       logical, intent(in)  :: day      ! True when daily averages wanted
       character(len=*), intent(in):: name ! Name of the variable 
                                           ! (used in  netCDF output)
       character(len=*), intent(in) :: unit ! Unit (writen in netCDF output)
       logical, intent(in), optional :: Is3D

       if ( present(Is3D) .and. Is3D ) then
         Nadded3d = Nadded3d + 1
         N = Nadded3d
         if ( me == 0 .and. MY_DEBUG  ) write(*,*) "Define 3d deriv ", N, name
         call CheckStop(N>MAXDEF_DERIV3D,"Nadded3d too big!")
         def_3d(N) = Deriv(class,avg,index,scale,rho,inst,year,month,day,&
                              name,unit)
       else
         Nadded2d = Nadded2d + 1
         N = Nadded2d
         if ( me == 0 .and. MY_DEBUG ) write(*,*) "Define 2d deriv ", N, name
         call CheckStop(N>MAXDEF_DERIV2D,"Nadded2d too big!")
         def_2d(N) = Deriv(class,avg,index,scale,rho,inst,year,month,day,&
                              name,unit)
       end if

    end subroutine AddDef
   !=========================================================================
    subroutine Define_Derived()

   ! Set the parameters for the derived parameters, including the codes
   ! used by DNMI/xfelt and scaling factors. (The scaling factors may
   ! be changed later in Derived_ml.
   ! And, Initialise the fields to zero.
   
    real, save    :: ugS = atwS*PPBINV/ATWAIR
    real, save    :: ugN = atwN*PPBINV/ATWAIR
    real, save    :: ugSO4, ugHCHO, ugCH3CHO
    real, save    :: ugPMad, ugPMde, ugSS  !advected and derived PM's & SeaS


  ! - for debug  - now not affecting ModelConstants version
   integer, dimension(MAXLIMAX) :: i_fdom
   integer, dimension(MAXLJMAX) :: j_fdom
   integer :: ind


    !   same mol.wt assumed for PPM25 and PPMco

     ugPMad = species(PM25)%molwt * PPBINV /ATWAIR 
     ugPMde = PPBINV /ATWAIR
     ugSS  = species( SSfi )%molwt * PPBINV /ATWAIR  !SeaS

     ugSO4 = species( SO4 )%molwt * PPBINV /ATWAIR
     ugHCHO   = species ( HCHO )%molwt * PPBINV /ATWAIR
     ugCH3CHO = species ( CH3CHO )%molwt * PPBINV /ATWAIR

!-- Deposition fields. Define all possible fields and their xfelt codes here:

      !code class  avg? ind scale rho Inst Yr Mn Day   name      unit  

Is3D = .false.
call AddDef( "PREC ", F, -1, 1.0,   F  , F  ,T ,T ,T ,"WDEP_PREC","mm") 
call AddDef( "WDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"WDEP_SOX","mgS/m2")
call AddDef( "WDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"WDEP_OXN","mgN/m2")
call AddDef( "WDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"WDEP_RDN","mgN/m2")

    ! Dry dep. --includes fields for ecosystem specific--- 
    ! ecosystem codes: SW = sea/water, CF = conif forest, DF = decid forest, 
    !                SN = seminatural (grass/moorlande/tundra)

      !code class  avg? ind scale rho Inst Yr Mn Day   name      unit  

call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"DDEP_SOX","mgS/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"DDEP_OXN","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"DDEP_RDN","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSSW","mgS/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSCF","mgS/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSDF","mgS/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSCR","mgS/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSSN","mgS/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSWE","mgS/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNSW","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNCF","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNDF","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNCR","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNSN","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNWE","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNSW","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNCF","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNDF","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNCR","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNSN","mgN/m2")
call AddDef( "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNWE","mgN/m2")

!-- 2-D fields - the complex ones
! (multiplied with roa in layers?? ==>  rho "false" ) !ds - explain!

!       code class  avg? ind scale rho  Inst  Yr  Mn   Day  name      unit 

call AddDef( "AOT  ", F, 20, 1.0,   F  , F  ,  T , T ,  F,"D2_AOT20","ppb h")
call AddDef( "AOT  ", F, 30, 1.0,   F  , F  ,  T , T ,  F,"D2_AOT30","ppb h")
call AddDef( "AOT  ", F, 40, 1.0,   F  , F  ,  T , T ,  F,"D2_AOT40","ppb h")
call AddDef( "AOT  ", F, 60, 1.0,   F  , F  ,  T , T ,  F,"D2_AOT60","ppb h")
call AddDef( "AOT  ", F, 30, 1.0,   F  , F  ,  T , F ,  F,"D2_AOT30f","ppb h")
call AddDef( "AOT  ", F, 40, 1.0,   F  , F  ,  T , F ,  F,"D2_AOT40f","ppb h")
call AddDef( "AOT  ", F, 60, 1.0,   F  , F  ,  T , F ,  F,"D2_AOT60f","ppb h")
call AddDef( "AOT  ", F, 40, 1.0,   F  , F  ,  T , F ,  F,"D2_AOT40c","ppb h")
!
! -- simple advected species. Note that some indices need to be set as dummys
!    in ACID, e.g. IXADV_O3
!
call AddDef( "ADV  ", T, IXADV_SO2, ugS, T, F , T , T , T ,"D2_SO2","ugS/m3")
call AddDef( "ADV  ", T, IXADV_SO4, ugS, T, F , T , T , T ,"D2_SO4","ugS/m3")
call AddDef( "ADV  ", T, IXADV_HNO3,ugN, T, F , T , T , T ,"D2_HNO3","ugN/m3")
call AddDef( "ADV  ", T, IXADV_PAN, ugN, T, F , T , T , T ,"D2_PAN","ugN/m3")
call AddDef( "ADV  ", T, IXADV_NH3, ugN, T, F , T , T , T ,"D2_NH3","ugN/m3")
call AddDef( "ADV  ", T, IXADV_NO , ugN, T, F , T , T , T ,"D2_NO","ugN/m3")
call AddDef( "ADV  ", T, IXADV_NO2, ugN, T, F , T , T , T ,"D2_NO2","ugN/m3")
call AddDef( "ADV  ", T,IXADV_aNH4, ugN, T, F , T , T , T ,"D2_aNH4","ugN/m3")
call AddDef( "ADV  ",T,IXADV_O3 ,PPBINV, F, F , T, T , T ,"D2_O3","ppb")
call AddDef( "ADV  ",T,IXADV_CO ,PPBINV, F, F , T, T , T ,"D2_CO","ppb")
call AddDef( "ADV  ",T,IXADV_aNO3, ugN,  T, F , T, T , T ,"D2_aNO3","ugN/m3")
call AddDef( "ADV ", T,IXADV_pNO3, ugN,  T, F , T, T , T ,"D2_pNO3", "ugN/m3")
call AddDef( "NOX  ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_NOX","ugN/m3")
call AddDef( "NOZ  ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_NOZ","ugN/m3")
call AddDef( "OX   ", T,   -1  ,PPBINV , F , F,T,T,T,"D2_OX","ppb")
call AddDef( "ADV  ",T,IXADV_PM25, ugPMad, T, F , T, T, T,"D2_PPM25","ug/m3")
call AddDef( "ADV  ",T,IXADV_PMco, ugPMad, T, F , T, T, T,"D2_PPMco","ug/m3")
!Sea salt
call AddDef( "ADV  ",T,IXADV_SSfi, ugSS, T, F , T, T, T,"D2_SSfi","ug/m3")
call AddDef( "ADV  ",T,IXADV_SSco, ugSS, T, F , T, T, T,"D2_SSco","ug/m3")
call AddDef( "PS    ",T,  0 ,       1.0, F , T, T, T, T ,"PS","hPa")
call AddDef( "HMIX  ",T,  0 ,       1.0, T , F, T, T, T ,"D2_HMIX","m")
call AddDef( "HMIX00",T,  0 ,       1.0, T , F, T, T, T ,"D2_HMIX00","m")
call AddDef( "HMIX12",T,  0 ,       1.0, T , F, T, T, T ,"D2_HMIX12","m")
!
! drydep
!   set as "external" parameters - ie set outside Derived subroutine
!      code class   avg? ind scale rho  Inst Yr  Mn   Day    name      unit 
call AddDef( "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_AFSTDF0","mmol/m2")
call AddDef( "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_AFSTDF16","mmol/m2")
call AddDef( "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_AFSTBF0","mmol/m2")
call AddDef( "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_AFSTBF16","mmol/m2")
!
call AddDef( "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_AFSTCR0","mmol/m2")
call AddDef( "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_AFSTCR3","mmol/m2")
call AddDef( "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_AFSTCR6","mmol/m2")
!
!      code class   avg? ind scale rho Inst Yr Mn  Day   name      unit 
call AddDef( "EXT  ", T, -1, 1.   , F, F,T ,T ,T ,"D2_O3DF   ","ppb")
call AddDef( "EXT  ", T, -1, 1.   , F, F,T ,T ,T ,"D2_O3WH   ","ppb")
!
! AOT30 and AOT40 for Wheat and Beech. May need daily here.
! Also, use field for EU definition (8-20 CET) and Mapping Manual/UNECE
! (daylight hours). 
! All of these use O3 at crop height, in contrast to the older AOT30, AOT40
! as defined above, and all allow daily output.
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_EUAOT30WH","ppb h")
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_EUAOT40WH","ppb h")
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_EUAOT30DF","ppb h")
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_EUAOT40DF","ppb h")
! UNECE:
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_UNAOT30WH","ppb h")
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_UNAOT40WH","ppb h")
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_UNAOT30DF","ppb h")
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_UNAOT40DF","ppb h")
!Mapping-Manual
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_MMAOT30WH","ppb h")
call AddDef( "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_MMAOT40WH","ppb h")
!
! --  time-averages - here 8-16
!
call AddDef( "TADV ", T,IXADV_HCHO  ,ugHCHO,  T, F, T, T, T,"D2T_HCHO","ug/m3")
call AddDef( "TADV ", T,IXADV_CH3CHO,ugCH3CHO,T, F, T, T,T,"D2T_CH3CHO","ug/m3")
call AddDef( "VOC  ", T,  -1    ,PPBINV, F, F, T, T, T,"D2_VOC","ppb")
!
! -- miscellaneous user-defined functions
!
!      code class   avg? ind scale rho Inst Yr Mn  Day   name      unit 
!! ,Deriv( "TSO4 ", T,   -1  ,ugS ,    T , F,T,T,T,"D2_SOX","ugS/m3")
call AddDef( "TOXN ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_OXN","ugN/m3")
call AddDef( "TRDN ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_REDN","ugN/m3")
call AddDef( "FRNIT", T,   -1  ,1.0 ,    F , F,T,T,T,"D2_FRNIT","(1)")
call AddDef( "MAXADV", F,IXADV_O3,PPBINV, F, F,T,T,T,"D2_MAXO3","ppb")
call AddDef( "MAXSHL", F,IXSHL_OH,1.0e13,F , F,T,F,T,"D2_MAXOH","?")
!
call AddDef( "tNO3 ", T, -1, ugN,    T, F, T, T, T,"D2_tNO3", "ugN/m3")
call AddDef( "SIA  ", T, -1, ugPMde, T, F, T, T, T,"D2_SIA" , "ug/m3")
call AddDef( "PMco ", T, -1, ugPMde, T, F, T, T, T,"D2_PMco", "ug/m3")
call AddDef( "PM25 ", T, -1, ugPMde, T, F, T, T, T,"D2_PM25", "ug/m3")
call AddDef( "PM10 ", T, -1, ugPMde, T, F, T, T, T,"D2_PM10", "ug/m3")
call AddDef( "H2O  ", T, -1,   1.0 , T, F, T, T, T,"D2_PM25_H2O ", "ug/m3")
call AddDef( "SSalt", T, -1, ugSS,   T, F, T, T, T,"D2_SS  ", "ug/m3") 
call AddDef( "SOM", F, 35, 1.,   F, F, T, T, F,"D2_SOMO35", "ppb day") 
call AddDef( "SOM", F,  0, 1.,   F, F, T, T, F,"D2_SOMO0", "ppb day") 

!-- 3-D fields

Is3D = .true.
call AddDef( "TH  ",T,  0 ,       1.0, F , T, T, T, F ,"D3_TH","m",Is3D)
call AddDef( "ADV  ", T, IXADV_O3 , PPBINV, F, T, T, T, F ,"D3_O3","ppb",Is3D)
call AddDef( "ADV  ", T, IXADV_SO2, PPBINV, F, T, T, T, F ,"D3_SO2","ppb",Is3D)
call AddDef( "ADV  ", T, IXADV_PAN, PPBINV, F, T, T, T, F ,"D3_PAN","ppb",Is3D)
call AddDef( "ADV  ", T, IXADV_HNO3,PPBINV, F, T, T, T, F ,"D3_HNO3","ppb",Is3D)
call AddDef( "ADV  ", T, IXADV_aNO3,PPBINV, F, T, T, T, F ,"D3_aNO3","ppb",Is3D)
call AddDef( "ADV  ", T, IXADV_NO2, PPBINV, F, T, T, T, F ,"D3_NO2","ppb",Is3D)
call AddDef( "VOC  ", T,       -1 , PPBINV, F, T, T, T, F ,"D3_VOC","ppb",Is3D)
call AddDef( "ADV  ", T, IXADV_aNH4,PPBINV, F, T, T, T, F ,"D3_aNH4","ppb",Is3D)
call AddDef( "ADV  ", T, IXADV_SO4, PPBINV, F, T, T, T, F ,"D3_SO4","ppb",Is3D)
call AddDef( "ADV  ", T, IXADV_H2O2,PPBINV, F, T, T, T, F ,"D3_H2O2","ppb",Is3D)
!
! Set Year true to allow debug - change later
call AddDef( "SHL",   T, IXSHL_OH,  PPTINV, T, F, T, T, F ,"D3_OH","?",Is3D)
call AddDef( "ADV",   T, IXADV_CH3COO2, &
                                    PPTINV, F, F, T, T, F ,"D3_CH3COO2","?",Is3D)
call AddDef( "MAX3DSHL", T,IXSHL_OH,PPTINV, T, F, T, T, F ,"D3_MAXOH","?",Is3D)   ! rho true for shl
call AddDef( "MAX3DADV", T, IXADV_CH3COO2,&
                                    PPTINV, F, F, T, T, F ,"D3_MAXCH3COO2","?",Is3D)
call AddDef( "PHNO3   ", T, IXSHL_PHNO3,1.0e8, F, F, T, T, F ,"D3_PHNO3","?",Is3D)
call AddDef( "MAX3DADV", T, IXADV_O3,PPBINV,F, F, T, T, F ,"D3_MAXO3","?",Is3D)


     if ( SOURCE_RECEPTOR .and. num_deriv2d>0 ) then  ! We assume that no
                                              ! daily outputs are wanted.
        def_2d(:)%day = .false.
     end if


     ! Get indices of wanted fields in larger def_xx arrays:

      do i = 1, num_deriv2d
          ind = find_index( wanted_deriv2d(i), def_2d(:)%name )
          f_2d(i) = def_2d(ind)
          if ( me == 0 .and. MY_DEBUG) write(*,*) "Index f_2d ", i, " = def ", ind
      end do

      do i = 1, num_deriv3d
          ind = find_index( wanted_deriv3d(i), def_3d(:)%name )
          f_3d(i) = def_3d(ind)
          if ( me == 0 .and. MY_DEBUG) write(*,*) "Index f_3d ", i, " = def ", ind
      end do

   !Initialise to zero

      if ( num_deriv2d > 0  ) d_2d( :,:,:,:) = 0.0
      if ( num_deriv3d > 0  ) d_3d( :,:,:,:,:) = 0.0

      debug_flag = ( MY_DEBUG .and. debug_proc ) 

  end subroutine Define_Derived
 !=========================================================================
     subroutine Setups()

    !/** flexibility note. By making use of character-based tests such
    !    as for "VOC" below, we achieve code which can stay for both ACID and
    !    OZONE without having to define non-used indices. 
    !    Similarly, we avoid the previous "if NUM_ACCSU eq 1" type test,
    !    since the appropriate code will now only activate 

    !/ ** if voc wanted, set up voc_array. Works for all ozone chemistries
    !     (and code not called for MADE-type).

      if ( any(  f_2d(:)%class == "VOC" ) .or. &
           any(  f_3d(:)%class == "VOC" )  ) then
            call Setup_VOC()
            if (MY_DEBUG)then
               write(6,*) "Derived VOC setup returns ", nvoc, "vocs"
               write(6,"(a12,/,(20i3))")  "indices ", voc_index(1:nvoc)
               write(6,"(a12,/,(20i3))")  "carbons ", voc_carbon(1:nvoc)
            endif
      end if


    end subroutine Setups
    !=========================================================================

    subroutine Derived(dt,End_of_Day)

    !/** DESCRIPTION
    !  Integration and averaging of chemical fields. Intended to be
    !  a more flexible version of the old chemint routine.
    !  Includes AOT40, AOT60 if present

      real, intent(in)    :: dt           !  time-step used in intergrations
      logical, intent(in) :: End_of_Day   !  e.g. 6am for EMEP sites

      character(len=len(f_2d%class)) :: typ  !  See defs of f_2d
      real :: thour                          ! Time of day (GMT)
      real :: timefrac                       ! dt as fraction of hour (3600/dt)
      real :: dayfrac              ! fraction of day elapsed (in middle of dt)
      integer :: ntime                       ! 1...NTDAYS
      integer :: nhour                       ! hours of day (GMT) 
      real, dimension(MAXLIMAX,MAXLJMAX) :: density !  roa (kgair m-3 when 
                                                    ! scale in ug,  else 1

      real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: inv_air_density3D
                ! Inverse of No. air mols/cm3 = 1/M 
                ! where M =  roa (kgair m-3) * MFAC  when ! scale in ug,  else 1
      logical :: accumulate_2dyear !flag to know when to accumulate d_2d
                                   ! (case "EXT")

      timefrac = dt/3600.0
      thour = current_date%hour+current_date%seconds/3600.0

      daynumber=day_of_year(current_date%year,current_date%month,&
                             current_date%day)

     !/***** 2-D fields **************************

     do n = 1, num_deriv2d

        accumulate_2dyear=.true.
        typ = f_2d(n)%class


        if ( f_2d(n)%rho ) then
            forall ( i=1:limax, j=1:ljmax )
                density(i,j) = roa(i,j,KMAX_MID,1)
            end forall
        else
            density(:,:) = 1.0
        end if

        !/** user-defined time-averaging. Here we have defined TADV and TVOC
        !    so that 8-hour daytime averages will be calculated. 
        !    Just comment out if not wanted, or (better!) don't define any
        !    f_2d as TADV or TVOC

        if ( typ == "TADV" .or. typ == "TVOC" ) then
             if(thour <= 8.0 .or. thour > 16.0 ) cycle  ! Start next species
        end if

       ! hmix average at 00 and 12:

        if ( typ == "HMIX00" .or. typ == "XKSIG00" ) then
             if(thour /= 0.0 ) cycle  ! Start next species
        end if

        if ( typ == "HMIX12" .or. typ == "XKSIG12" ) then
             if(thour /= 12.0 ) cycle  ! Start next species
        end if

        index = f_2d(n)%index
        !if ( My_DEBUG ) then
        !   write(*,*) "DEBUG Derived 2d", n, f_2d(n)%name, index, typ
        !end if

        select case ( typ )

          case ( "PS" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = ps(i,j,1)*0.01
            end forall


          case ( "HMIX", "HMIX00", "HMIX12" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = pzpbl(i,j)  
            end forall

            if ( debug_flag ) then
             write(*,fmt="(a12,2i4,4f12.3)") "HMIX" , n , d_2d(n,debug_li,debug_lj,IOU_INST)       
            end if

         ! Simple advected species:
          case ( "ADV", "TADV" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = xn_adv(index,i,j,KMAX_MID)  &
                                     * cfac(index,i,j) * density(i,j)  
            end forall

            if ( debug_flag ) then
             write(*,fmt="(a12,2i4,4f12.3)") "JUST ADV" , n, index  &
              ,d_2d(n,debug_li,debug_lj,IOU_INST)*PPBINV &
              ,xn_adv(index,debug_li,debug_lj,KMAX_MID)*PPBINV &
              ,density(debug_li,debug_lj), cfac(index,debug_li,debug_lj)
            end if

          case ( "H2O" )      !water

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = PM_water(i,j,KMAX_MID)  
            end forall
  

          case ( "MAXADV" )


              d_2d( n, 1:limax,1:ljmax,IOU_DAY) = max( d_2d( n, 1:limax,1:ljmax,IOU_DAY),  &
             xn_adv(index,1:limax,1:ljmax,KMAX_MID)  &
                                     * cfac(index,1:limax,1:ljmax) * density(1:limax,1:ljmax))


            if ( debug_flag ) then
             write(*,fmt="(a12,2i4,4f12.3)") "ADV MAX. ", n, index  &
                      , d_2d(n,debug_li,debug_lj,IOU_DAY) * PPBINV      &
                      ,  xn_adv(index,debug_li,debug_lj,KMAX_MID)* PPBINV  &
                      ,  density(debug_li,debug_lj), cfac(index,debug_li,debug_lj)

            end if

            !Monthly and yearly ARE averaged over days
            if(End_of_Day)then
              d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)
              nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
              if(    current_date%month >= 4 &
                 .or.current_date%month <= 9 )then
              d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY) 
              nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
              endif
            endif


          case ( "MAXSHL" )        ! Daily maxima - short-lived

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_DAY) = max( d_2d( n, i,j,IOU_DAY), &
                                xn_shl(index,i,j,KMAX_MID)  &
                                    / (density(i,j)*MFAC) )
                                   !u4  / (roa(:,:,KMAX_MID,1)*MFAC) )
            end forall


            if ( debug_flag ) then
               write(*, *) "SHL:MAX.,MFAC ", n, index  , MFAC
               write(*,fmt="(a12,2i4,4es12.3)") "SHL MAX. ", n, index  &
                      , d_2d(n,debug_li,debug_lj,IOU_DAY) &
                      ,  xn_shl(index,debug_li,debug_lj,KMAX_MID)  &
                      ,  density(debug_li,debug_lj), MFAC
            end if

            !Monthly and yearly ARE averaged over days
            if(End_of_Day)then
              d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)
              nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
              if(    current_date%month >= 4 &
                 .or.current_date%month <= 9 )then
              d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY) 
              nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
              endif
            endif

          case ( "VOC", "TVOC" )

            call voc_2dcalc()

          case( "AOT" )

            call aot_calc( n, timefrac )

           if( debug_flag .and. i == debug_li .and. j == debug_lj ) then
              write(*,*) "GROWINDERIV? ", n, f_2d(n)%name
           end if

            if(     current_date%month<startmonth_forest&
                 .or.current_date%month>endmonth_forest)then
               if( f_2d(n)%name=="D2_AOT30f".or.& 
                   f_2d(n)%name=="D2_AOT40f".or.&
                   f_2d(n)%name=="D2_AOT60f")then
                   accumulate_2dyear=.false.
               endif

            endif
            if(     current_date%month<startmonth_crops&
               .or.current_date%month>endmonth_crops)then
               if( f_2d(n)%name=="D2_AOT30c".or.&
                   f_2d(n)%name=="D2_AOT40c".or.&
                   f_2d(n)%name=="D2_AOT60c")then
                   accumulate_2dyear=.false.
               endif
            endif

           case( "SOM" )


              !dt/7200: half a dt time step in hours
              !dayfrac "points" to the middle of the integration step
              dayfrac= (thour-(dt/7200.))/24. !must be < 1
              ntime=int(dayfrac*NTDAY )+1 !must be >=1 and <= NTDAY
              if(dayfrac<0)ntime=NTDAY !midnight

        !last value  (not averaged): 
          D2_O3_DAY( : , : , ntime) =&
           xn_adv(IXADV_O3,:,:,KMAX_MID)*cfac(IXADV_O3,:,:)*PPBINV

              if(dayfrac<0)then !only at midnight: write on d_2d

                 
                 call som_calc( n ) !  accumulate
                 d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY) 

                ! if(current_date%month>=4.and.current_date%month<=9)then
                 d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY) 
                !NB overwritten anyway D2_O3_DAY = 0.
              endif


          case ( "PREC", "WDEP", "DDEP" )
            if ( debug_flag ) write(*,"(a18,i4,a12,a4,es12.3)")"PR/DEP d_2d",&
                   n, f_2d(n)%name, " is ", d_2d(n,debug_li,debug_lj,IOU_INST)

          case ( "EXT" )

          ! Externally set for IOU_INST (in other routines); so no new work 
          ! needed except decision to accumalate to yearly or not.
          ! Used for e.g. AOT40s
             call setaccumulate_2dyear(n,accumulate_2dyear)
            if ( debug_flag ) write(*,"(a18,i4,a12,a4,es12.3)")"EXT d_2d",&
                   n, f_2d(n)%name, " is ", d_2d(n,debug_li,debug_lj,IOU_INST)

          case  default

            if ( debug_flag ) then
                 write(*,*) "My_Deriv Defaults called for n=", n, "Type ",typ, "Name ", f_2d(n)%name
                 write(*,*) "My_Deriv index?, avg?, nav? length?, class? ", index,&
                    f_2d(n)%avg, nav_2d(n,IOU_INST), len(f_2d%class), f_2d(n)%class
             end if 

             call My_DerivFunc( d_2d(n,:,:,IOU_INST), n, typ, timefrac, density ) 

        end select


        !/** add to daily, monthly and yearly average, and increment counters
        !  Note that the MAXADV and MAXSHL and SOM needn't be summed here, but
        !  since the INST values are zero it doesn't harm, and the code is 
        !  shorter. These d_2d ( MAXADV, MAXSHL, SOM) are set elsewhere

        d_2d(n,:,:,IOU_DAY )  = d_2d(n,:,:,IOU_DAY )  + d_2d(n,:,:,IOU_INST) 
        if ( f_2d(n)%avg ) nav_2d(n,IOU_DAY) = nav_2d(n,IOU_DAY) + 1
        d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_INST) 
        if ( f_2d(n)%avg ) nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
        if(accumulate_2dyear)then
           d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_INST) 
           if ( f_2d(n)%avg ) nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
        endif
 
     end do   ! num_deriv2d

     !/***** 3-D fields **************************

       if(debug_flag) then ! RUN through indices etc.
            write(*, "(a12,2i4,f12.3)") "3D3D TIME ",  me, num_deriv3d, &
                     (current_date%hour+current_date%seconds/3600.0)
        end if


     do n = 1, num_deriv3d

        index = f_3d(n)%index

       if ( f_3d(n)%rho ) then
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                inv_air_density3D(i,j,k) = 1.0/( roa(i,j,k,1) * MFAC )
            end forall
        else
            inv_air_density3D(:,:,:) = 1.0
        end if

        select case ( f_3d(n)%class )

         ! Simple advected species:
          case ( "ADV" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
            end forall

         case ( "BGN" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_bgn(index,i,j,k)
            end forall

         case ("XKSIG00", "XKSIG12" ) !hf hmix xksig

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xksig(i,j,k)
            end forall

         case ("TH  " ) !JEJ Pot. temp (needed for cross sections)

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = th(i,j,k,1)
            end forall

         case ( "PHNO3" )   !ds-hf  rv1_9_28
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_shl(index,i,j,k)
            end forall

             if(debug_flag) write(*,"(a12,i4,2es12.3)") "3D3D PHNO3", n, &
                  xn_shl(index,debug_li,debug_lj,KMAX_MID), &
                  d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST)

         case ( "MAX3DSHL" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )! Daily maxima - short-lived
              d_3d( n, i,j,k,IOU_INST) = max( d_3d( n, i,j,k,IOU_INST),&
                                      xn_shl(index,i,j,k) &
                                     * inv_air_density3D(i,j,k) )
            end forall

            if(debug_flag) write(*,"(a13,i4,f8.3,3es12.3)") "3D3D MAX3DSHL", n, thour, &
              xn_shl(index,debug_li,debug_lj,KMAX_MID), &
              1.0/inv_air_density3D(debug_li,debug_lj,KMAX_MID), &
              d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST)

          case ( "MAX3DADV" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) =  max( d_3d( n, i,j,k,IOU_INST),&
                                               xn_adv(index,i,j,k) )
            end forall

             if(debug_flag) write(*,"(a12,i4,f8.3,4es12.3)") "SET MAX3DADV", n, thour, &
                      xn_adv(index,debug_li,debug_lj,KMAX_MID), &
                      d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST)

          case ( "SHL" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) =   xn_shl(index,i,j,k) * inv_air_density3D(i,j,k)
            end forall


          case ( "VOC" )

            call voc_3dcalc()

          case  default

           write(unit=errmsg,fmt=*) "Derived 3D class NOT FOUND", n, index, &
                         f_3d(n)%name,f_3d(n)%class
           call CheckStop( errmsg )


        end select
     

      !/** add to monthly and yearly average, and increment counters
       !    ( no daily averaging done for 3-D fields so far).


       ! For the MAX3D possibilities, we store maximum value of the
       !   current day in the IOU_INST variables.
       !   These are then added into IOU_MON **only** at the end of each day.
       ! (NB there is an error made on 1st day used, since only 1st 6 hours
       !  are checked. Still, not much happens on 1st Jan.... ;-)

        if ( (f_3d(n)%class == "MAX3DSHL")  .or. &
             (f_3d(n)%class == "MAX3DADV") )then
           if (End_of_Day) then
              d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
                                     + d_3d(n,:,:,:,IOU_INST)
              d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
                                     + d_3d(n,:,:,:,IOU_INST)
              if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1 !only collected for end_of_day

              if( debug_flag ) then
                    write(*,fmt="(a20,a9,i4,f8.3,2es12.3)") "END_OF_DAY MAX3D", &
                      f_3d(n)%class, n, thour,  &
                      d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_MON ),&
                      d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST )
                    write(*,"(a20,i4,2x,6i6)") "END_OF_DAY NAV ", &
                      n, (nav_3d(n,i), i=1,LENOUT3D)
              end if

              d_3d(n,:,:,:,IOU_INST ) = 0.0  !! Reset d_3d

           endif ! End_of_Day
        else
           d_3d(n,:,:,:,IOU_DAY ) = d_3d(n,:,:,:,IOU_DAY ) &
                + d_3d(n,:,:,:,IOU_INST)
           d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
                + d_3d(n,:,:,:,IOU_INST)
           d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
                + d_3d(n,:,:,:,IOU_INST)
           if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1
        endif



!      !/** add to monthly and yearly average, and increment counters
!      !    ( no daily averaging done for 3-D fields so far).
!
!       d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
!                              + d_3d(n,:,:,:,IOU_INST)
!       d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
!                              + d_3d(n,:,:,:,IOU_INST)
!
!       if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1
 
      end do
    end subroutine Derived
    !=========================================================================

    subroutine DerivedProds(text,dt)

    !/** DESCRIPTION
    !  Calculates chemical changes by comparing values before and  after 
    !  chemistry subroutine. Intended to be a more flexible version of the old 
    !  PRODO3  calculation

      character(len=*), intent(in) :: text  ! "Before" or "After"
      real,             intent(in) :: dt    ! timestep (s)

      real :: timefrac                      ! dt as fraction of hour (3600/dt)



      if (.not. any( f_3d%class == "PROD" ) ) return

      timefrac = dt/3600.0

     !/***** 3-D fields **************************

     do n = 1, num_deriv3d

        if ( f_3d(n)%class  == "PROD " ) then
           index = f_3d(n)%index

           select case ( text )

               case ( "Before" )   !! Initialise to xn_adv

                 forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                   d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
                 end forall

               case ( "After" )    !! Calculate change

                 forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                   d_3d( n, i,j,k,IOU_INST) = &
                      d_3d( n, i,j,k,IOU_INST) - xn_adv(index,i,j,k)
                 end forall

           end select
        end if
      end do
     
    end subroutine DerivedProds
    !=========================================================================

    subroutine ResetDerived(period)
      integer, intent(in) :: period   ! Either IOU_DAY or IOU_MON

       if ( period <= LENOUT2D ) then
           nav_2d  (:,period) = 0.0
           d_2d(:,:,:,period) = 0.0
       end if 


       if ( period <= LENOUT3D ) then
           nav_3d    (:,period) = 0.0
           d_3d(:,:,:,:,period) = 0.0
       end if

    end subroutine ResetDerived
 !=========================================================================

  subroutine Setup_VOC()
      !--------------------------------------------------------
      ! Searches through the advected species and colects the
      ! index and carbon content of nmhc/voc species, as they were
      ! defined in GenOut_ml
      !
      !--------------------------------------------------------
       integer :: n
   
      do n = 1, NSPEC_ADV

        if ( species( NSPEC_SHL+n )%carbons > 0 .and. &
             species( NSPEC_SHL+n )%name   /= "CO"  .and. &
             species( NSPEC_SHL+n )%name   /= "CH4" ) then

             nvoc = nvoc + 1
             voc_index(nvoc) = n
             voc_carbon(nvoc) = species( NSPEC_SHL+n )%carbons
        end if
      end do
  end subroutine Setup_VOC
 !=========================================================================

   subroutine voc_2dcalc()

    !/-- Sums up voc species using the indices defined earlier in Setup_VOCs

     ! We initialise d_2d first, the use a simple loop
     ! over voc. Some CPU could be saved by initialising
     ! with the 1st voc, then looping over 2, nvoc, but who cares...

      
      d_2d( n, 1:limax,1:ljmax,IOU_INST) =  0.0

      do ivoc = 1, nvoc

         index = voc_index(ivoc)           ! Gives which IXADV_ to use.
         forall ( i=1:limax, j=1:ljmax )
             d_2d( n, i,j,IOU_INST) = d_2d( n, i,j,IOU_INST)      &
                                    + xn_adv(index,i,j,KMAX_MID)  &
                                    * voc_carbon(ivoc) * cfac(index,i,j)
                               ! multiplied by nr. of C and "reduced to surface"
         end forall
      end do ! ivoc
   end subroutine voc_2dcalc

 !=========================================================================
   subroutine voc_3dcalc()

    !/-- as for voc_2dcalc

      d_3d( n, 1:limax,1:ljmax,1:KMAX_MID,IOU_INST) =  0.0

      do ivoc = 1, nvoc

         index = voc_index(ivoc)
         forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
             d_3d( n, i,j,k,IOU_INST) = d_3d( n, i,j,k,IOU_INST) + &
                     xn_adv(index,i,j,k)*voc_carbon(ivoc)
         end forall
      end do ! ivoc

   end subroutine voc_3dcalc
 !=========================================================================

  subroutine aot_calc( n, timefrac )

    !/-- Calcuates AOT values for input threshold. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )
    !    Only relevant in ozone models, so far....

    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    real,    intent(in) :: timefrac    ! Timestep as fraction of hour

    real    :: threshold               ! Threshold, e.g. 40 or 60 (ppb)
    integer :: izen                    ! integer of zenith angle
    real :: o3                         ! Ozone (ppb) - needed if AOTs

     threshold = f_2d(n)%index

      do i=1,limax
        do j=1,ljmax

           izen = max(1,int( zen(i,j) + 0.5))

           if ( izen < AOT_HORIZON ) then
                o3 = xn_adv(IXADV_O3,i,j,KMAX_MID) &
                     * cfac(IXADV_O3,i,j) * PPBINV 

                o3 = max( o3 - threshold , 0.0 )   ! Definition of AOTs

             ! d_2d values will be accumulated in Derived_ml

              d_2d(n, i,j,IOU_INST ) = o3 * timefrac  

           else
               d_2d(n, i,j,IOU_INST ) = 0.0   
           end if
        end do
      end do
   end subroutine aot_calc

!=========================================================================

  subroutine som_calc( n )


    !/-- Calculates SOM (8hours) values for input threshold.  !pw rv2_1

    implicit none
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays

    real    :: threshold               ! Threshold, e.g. 35 (ppb)
    real :: o3                         ! Ozone (ppb) - needed if SOMs
    real :: sum8h
    integer, parameter :: N8h = (NTDAY*8)/24 !number of periods in 8 hours
    real, parameter :: N8h_inv=1./N8h
    integer :: nh


    threshold = f_2d(n)%index

      do i=1,limax
        do j=1,ljmax

           !find max running 8h sum O3
           sum8h=0.
           do nh=1,N8h
              sum8h = sum8h + D2_O3_DAY( i , j , nh)
           enddo
           o3=sum8h
           do nh=N8h+1,NTDAY
              sum8h =sum8h-D2_O3_DAY( i , j , nh-N8h)+D2_O3_DAY( i , j , nh)
              o3=max(o3,sum8h)
              if(n<0)write(*,*)o3 !pw fake for compiler!!
           enddo

           !divide by N8h to find 8h mean 
           o3=o3*N8h_inv

           o3 = max( o3 - threshold , 0.0 )   ! Definition of SOMs

             ! d_2d values will be accumulated in Derived_ml

           d_2d(n, i,j,IOU_DAY ) = o3  

        end do
      end do
   end subroutine som_calc

 !=========================================================================

   subroutine setaccumulate_2dyear(n,accumulate_2dyear)

! We don't want the yearly output to accumulate over the whole year
     integer, intent(in) :: n
      logical, intent(inout) :: accumulate_2dyear !flag to know when to 
                                                  !accumulate d_2d (case "EXT")

      if( f_2d(n)%name=="D2_EUAOT30DF".or.&
          f_2d(n)%name=="D2_EUAOT40DF".or.&
          f_2d(n)%name=="D2_UNAOT30DF".or.&
          f_2d(n)%name=="D2_UNAOT40DF"    &
          )then
         if(   current_date%month<startmonth_forest&
              .or.current_date%month>endmonth_forest)then
            accumulate_2dyear=.false.
         endif   
      endif

       if(f_2d(n)%name=="D2_EUAOT30WH".or.&
          f_2d(n)%name=="D2_EUAOT40WH".or.&
          f_2d(n)%name=="D2_UNAOT30WH".or.&
          f_2d(n)%name=="D2_UNAOT40WH"    &
          )then
         if(   current_date%month<startmonth_crops&
              .or.current_date%month>endmonth_crops)then
            accumulate_2dyear=.false.

         endif   
      endif

    end subroutine setaccumulate_2dyear

end module Derived_ml
