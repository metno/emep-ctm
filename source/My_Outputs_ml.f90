! <My_Outputs_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                          module  My_Outputs_ml

! MOD OD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! -----------------------------------------------------------------------
! Allows user to specify which species are output to various
!  ascii and binary output files.
!
! Sites  - surface sites,     to sites.out
! Sondes - vertical profiles, to sondes.out
! Hourly - ascii output of selected species, selcted domain
! -----------------------------------------------------------------------

  use CheckStop_ml,     only: CheckStop
  use Derived_ml,       only: f_2d, d_2d
  use GenSpec_adv_ml
  use GenSpec_shl_ml,    only: IXSHL_OH,IXSHL_HO2
  use GenChemicals_ml ,  only: species
  use ModelConstants_ml, only: PPBINV, PPTINV, ATWAIR, atwS, atwN, NPROC
  use Par_ml,            only: me, GIMAX,GJMAX,IRUNBEG,JRUNBEG
  use SmallUtils_ml,     only: find_index
  use TimeDate_ml,       only: date 
 
  implicit none

INCLUDE 'mpif.h'
INTEGER STATUS(MPI_STATUS_SIZE),INFO
logical, public, parameter :: out_binary = .false.
logical, public, parameter :: Ascii3D_WANTED = .false.


! Site outputs   (used in Sites_ml)
!==============================================================
! Specify the species to be output to the sites.out file
! For met params we have no simple index, so we use characters.
! These must be defined in Sites_ml.f90.

integer, private :: isite              ! To assign arrays, if needed
integer, public, parameter :: &
     NSITES_MAX =    99         & ! Max. no surface sites allowed
    ,FREQ_SITE  =    1          & ! Interval (hrs) between outputs
    ,NADV_SITE  =    NSPEC_ADV  & ! No. advected species (1 up to NSPEC_ADV)
    ,NSHL_SITE  =    1          & ! No. short-lived species
    ,NXTRA_SITE =    2            ! No. Misc. met. params  ( now T2)

   integer, public, parameter, dimension(NADV_SITE) :: &
    SITE_ADV =  (/ (isite, isite=1,NADV_SITE) /)  ! Everything

   integer, public, parameter, dimension(NSHL_SITE) :: &
    SITE_SHL =  (/ IXSHL_OH /)                    ! More limited

! Extra parameters - need to be coded in Sites_ml also. So far
! we can choose from T2, or th (pot. temp.)

   character(len=10), public, parameter, dimension(NXTRA_SITE) :: &
     SITE_XTRA=  (/ "hmix    ", "th      "/) 
!    SITE_XTRA=  (/ "hmix    ", "Vg_ref  ", "Vg_1m   " /) 

 !    - can access d_2d fields through index here, by
 !     setting "D2D" above and say D2_FSTCF0 here:

   integer,           public, parameter, dimension(NXTRA_SITE) :: &
    SITE_XTRA_INDEX=  (/     0,        0 /)     ! Height at mid-cell



   !/*** Aircraft outputs   (used in Polinat_ml)
   !==============================================================
   !   Specify the species to be output by Polinat for aircraft flight tracks

   integer, public, parameter :: &
     NFLIGHT_MAX =    10               &   ! Max. no sondes allowed
    ,FREQ_FLIGHT =    12               &   ! Interval (hrs) between outputs
    ,NADV_FLIGHT =    1                    ! No.  advected species

   integer, public, parameter, dimension(NADV_FLIGHT) :: &
    FLIGHT_ADV =  (/ IXADV_O3 /)


   !/*** Sonde outputs   (used in Sites_ml)
   !==============================================================
   !     Specify the species to be output to the sondes.out file
   !  We typically deal with fewer species for sonde output than
   !  surface sites, so we use a different method to specify.
   ! For met params we have no simple index, so we use characters.
   ! These must be defined in Sites_ml.f90.

   integer, public, parameter :: &
     NSONDES_MAX =    99               &   ! Max. no sondes allowed
    ,NLEVELS_SONDE =  20               &   ! No. k-levels (9 => 0--2500 m)
    ,FREQ_SONDE  =    12               &   ! Interval (hrs) between outputs
    ,NADV_SONDE  =     8                &   ! No.  advected species
    ,NSHL_SONDE  =    1                &   ! No. short-lived species
    ,NXTRA_SONDE =    4                &   ! No. Misc. met. params  (now th)
    ,N_NOy       =   10                    ! # of N species in NOy

   integer, public, parameter, dimension(NADV_SONDE) :: &
   SONDE_ADV =  (/ IXADV_O3, IXADV_NO2, IXADV_HNO3, IXADV_aNO3, & 
   IXADV_pNO3,  IXADV_SO4,  IXADV_aNH4, IXADV_NH3/)


   integer, public, parameter, dimension(N_NOy) :: &
     NOy_SPEC =  (/ IXADV_HNO3, IXADV_NO,  IXADV_NO2,  IXADV_PAN,    &
                    IXADV_MPAN, IXADV_NO3, IXADV_N2O5, IXADV_ISONO3, &
                    IXADV_ISNI, IXADV_ISNIR /)

   integer, public, parameter, dimension(NSHL_SONDE) :: &
    SONDE_SHL =  (/ IXSHL_OH /)
   character(len=10), public, parameter, dimension(NXTRA_SONDE) :: &
!    SONDE_XTRA=  (/ "PM25 ", "PMco ", "NOy  ", "z_mid", "p_mid", "th   " /) 
    SONDE_XTRA=  (/ "NOy  ", "z_mid", "p_mid", "th   " /) !Height at mid-cell

 !   can access d_3d fields through index here, by
 !   setting "D3D" above and say D3_XKSIG12 here:

   integer,           public, parameter, dimension(NXTRA_SONDE) :: &
                    SONDE_XTRA_INDEX=  (/      0, 0, 0, 0 /)



   !====================================================================
   !/*** Hourly outputs   (from hourly_out routine) to print out
   !     concentrations  or even met. parameters every hour 
   !     (or multiple: HOURLY_FREQ) for specified sub-grid.
   !     Note: as to met. parameters, only temp2m Th arespecified   
   !           so far- others need change in hourly_out.f also).

   !-------------------------------------------------------------------
   !  Possibility of multi-layer output. Specify NLEVELS_HOURLY here
   !  and in hr_out defs use either:
   !
   !  ADVppbv to get surface concentrations (only relevant for layer k=20 
   !  while gives meaningless  number for higher levels.
   !
   !  Or BCVppbv to get grid-centre concentrations (relevant for all layers)
   !----------------------------------------------------------------

    logical, public, parameter :: Hourly_ASCII = .false.
     ! Hourly_ASCII = .True. gives also Hourly files in ASCII format.

    integer, public, parameter :: NHOURLY_OUT =  1 ! No. outputs
    integer, public, parameter :: NLEVELS_HOURLY = 1 ! No. outputs
    integer, public, parameter :: FREQ_HOURLY = 1  ! 1 hours between outputs

    type, public:: Asc2D
         character(len=12):: name   ! Name (no spaces!)
         character(len=7) :: type   ! "ADVppbv" or "ADVugm3" or "SHLmcm3" 
         character(len=9) :: ofmt   ! Output format (e.g. es12.4)
         integer          :: spec   ! Species number in xn_adv or xn_shl array
                                    !  .. or other arrays
         integer          :: ix1    ! bottom-left x
         integer          :: ix2    ! upper-right x
         integer          :: iy1    ! bottom-left y
         integer          :: iy2    ! upper-right y
         integer          :: nk     ! number of vertical levels
         character(len=12) :: unit   ! Unit used 
         real             :: unitconv   !  conv. factor
         real             :: max    ! Max allowed value for output
    end type Asc2D
     
    type(Asc2D), public, dimension(NHOURLY_OUT) :: hr_out  ! Set below


  !/** wanted binary dates... specify days for which full binary
  !    output is wanted. Replaces the hard-coding which was
  !    in wrtchem:

     integer, public, parameter :: NBDATES = 3  
     type(date), public, save, dimension(NBDATES) :: wanted_dates_inst 

 !================================================================

   public :: set_output_defs

   contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine set_output_defs
   implicit none

   character(len=44) :: errmsg  ! Local error message
   integer           :: i       ! Loop index

   real                :: to_ug_S & ! conversion to ug of S
                         ,to_ug_N & ! conversion to ug of N
                         ,to_mgSIA& ! conversion to mg of N
                         ,to_ugSIA  ! conversion to ug of N 
   real, save :: m_s = 100.0 ! From cm/s to m/s

  ! introduce some integers to make specification of domain simpler
  ! and less error-prone. Numbers can be changed as desired.

   !integer, save :: ix1 = 36, ix2 = 167, iy1=12, iy2 =  122  !EMEP
   integer, save :: ix1 = 65, ix2 = 167, iy1=12, iy2 =  122  !restricted EMEP
!   integer, save :: ix1 = IRUNBEG, ix2 = IRUNBEG+GIMAX-1,  &
!                    iy1=JRUNBEG,   iy2 =  JRUNBEG+GJMAX-1   ! all

  ! For Deriv system:
   integer :: D2_O3WH, D2_O3DF,  &
    D2_AFSTDF0, D2_AFSTDF16, D2_AFSTBF0, D2_AFSTBF16, &    ! JUN06
    D2_AFSTCR0, D2_AFSTCR3, D2_AFSTCR6,&
    D2_AFSTCN0, D2_AFSTCN3, D2_AFSTCN6

! WARNING: If the specification of the subdomain is different for
!            different components (ix1=125 for ozone and ix1=98 for 
!            NH4 for example) , the variables i_EMEP, j_EMEP
!            latitude and longitude in NetCDF output will be
!            wrong. 

  to_ug_S = atwS*PPBINV/ATWAIR    ! for output in ug(S)/m3
  to_ug_N = atwN*PPBINV/ATWAIR    ! for output in ug(N)/m3
  to_mgSIA= PPBINV/ATWAIR*1000.0
  to_ugSIA= PPBINV/ATWAIR

 !/** Hourly outputs
 !    Note that the hourly output uses **lots** of disc space, so specify
 !    as few as you need and with as small format as possible (cf max value).

 ! ** REMEMBER : ADV species are mixing ratioes !!
 ! ** REMEMBER : SHL species are in molecules/cm3, not mixing ratio !!
 ! ** REMEMBER : No spaces in name, except at end !!

  !**               name     type   
  !**                ofmt   ispec     ix1 ix2  iy1 iy2  nk unit conv  max

 hr_out(1)= Asc2D("o3_3m", "ADVppbv", &
                  "(f9.4)",IXADV_o3, ix1,ix2,iy1,iy2,1, "ppbv",PPBINV,600.0)

! For deriv system
 
 D2_O3WH = find_index("D2_O3WH",f_2d(:)%name)
 D2_O3DF = find_index("D2_O3DF",f_2d(:)%name)
 D2_AFSTDF16 = find_index("D2_AFSTDF16",f_2d(:)%name)
 D2_AFSTCR3 = find_index("D2_AFSTCR3",f_2d(:)%name)

! hr_out(2)= Asc2D("O3_Wheat", "D2D", &
!                  "(f7.3)", D2_O3WH, ix1,ix2,iy1,iy2,1, "ppbv", 1.0  ,600.0)
! hr_out(3)= Asc2D("O3_Beech", "D2D", &
!                  "(f7.3)", D2_O3DF, ix1,ix2,iy1,iy2,1, "ppbv", 1.0  ,600.0)
! hr_out(4)= Asc2D("FST_DF00", "D2D", &
!                 "(f7.3)", D2_FSTDF00, ix1,ix2,iy1,iy2,1, "NNNN", 1.0  ,600.0)
! hr_out(5)= Asc2D("FST_WH00", "D2D", &
!                  "(f7.3)", D2_FSTWH00, ix1,ix2,iy1,iy2,1, "NNNN", 1.0  ,600.0)

!  hr_out(1)=  Asc2D("Ozone", "ADVppbv", &
!                  "(f9.5)",IXADV_O3, ix1,ix2,iy1,iy2, "ppb",PPBINV,600.0)
!  hr_out(2)=  Asc2D("aNH4-air", "ADVugm3", &
!                  "(f8.4)",IXADV_aNH4, ix1,ix2,iy1,iy2, "ug",to_ugSIA,600.0)
!  hr_out(3)=  Asc2D("aNO3-air", "ADVugm3", &
!                  "(f8.4)",IXADV_aNO3, ix1,ix2,iy1,iy2, "ug",to_ugSIA,600.0)
!  hr_out(4)=  Asc2D("SO4-air", "ADVugm3", &
!                  "(f8.4)",IXADV_aNO3, ix1,ix2,iy1,iy2, "ug",to_ugSIA,600.0)
!  hr_out(5)=  Asc2D("pNO3-air", "ADVugm3", &
!                  "(f8.4)",IXADV_pNO3, ix1,ix2,iy1,iy2, "ug",to_ugSIA,400.0)

!  hr_out(2)= &
!   Asc2D("ADVugm3", "(f8.4)",IXADV_aNH4, 45, 170, 1, 133, "ug",to_ugSIA,600.0)
!  hr_out(3)= &
!   Asc2D("ADVugm3", "(f8.4)",IXADV_aNO3, 45, 170, 1, 133, "ug",to_ugSIA,600.0)
!  hr_out(4)= &
!   Asc2D("ADVugm3", "(f8.4)",IXADV_SO4, 45, 170, 1, 133, "ug",to_ugSIA,400.0)
!  hr_out(5)= &
!   Asc2D("ADVugm3", "(f8.4)",IXADV_pNO3, 45, 170, 1, 133, "ug",to_ugSIA,400.0)
!
 ! Extra parameters - need to be coded in Sites_ml also. So far
 ! we can choose from T2, or th (pot. temp.)
 !  - or from d_2d arrays.

  !**           type   ofmt   ispec    ix1 ix2  iy1 iy2  unit conv    max
  !hr_out(3)= &
  !   Asc2D("D2D", "(f6.1)",   D2_HMIX, ix1,ix2,iy1,iy2, "m",1.0   ,10000.0)

!Flux stuff
!   hr_out(2)= Asc2D("Fst_TConif  ", "D2D", "(f8.5)",&
!                  D2_FSTCF0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!   hr_out(3)= Asc2D("Fst_TBroad  ", "D2D", "(f8.5)",&
!                  D2_FSTDF0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!   hr_out(4)= Asc2D("Fst_Grass   ", "D2D", "(f8.5)",&
!                  D2_FSTGR0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!   hr_out(5)= Asc2D("Fst_Wheat   ", "D2D", "(f8.5)",&
!                  D2_FSTWH0, ix1,ix2,iy1,iy2, "nmole/m2/s", 1.0  ,900.0)
!    hr_out(10)= Asc2D("O3__Conif   ", "D2D", "(f7.3)",&
!                   D2_O3CF,   ix1,ix2,iy1,iy2, "ppb       ", 1.0  ,900.0)

 !/** theta is in deg.K
 !hr_out(1)=  Asc2D("T2_C",   "T2_C   ", &
 !                "(f5.1)",     -99, ix1,ix2,iy1,iy2, "degC",1.0   ,100.0)
 !hr_out(2)=  Asc2D("Precip", "PRECIP ", &
 !                "(f11.7)",    -99, ix1,ix2,iy1,iy2, "mm/hr",1.0,  200.0)
 !hr_out(3)=  Asc2D("Idir",   "Idirect", &
 !                "(f5.1)",    -99, ix1,ix2,iy1,iy2, "umole/m2/s",1.0, 1000.0)
 !hr_out(4)=  Asc2D("Idif",   "Idiffus", &
 !                "(f5.1)",    -99, ix1,ix2,iy1,iy2, "umole/m2/s",1.0, 1000.0)



  !/** Consistency checks

   do i = 1, NHOURLY_OUT

       ! We use ix1 to see if the array has been set.

        if ( hr_out(i)%ix1 < 1 .or.  hr_out(i)%ix1 > 999 ) then

          write(unit=errmsg,fmt=*) "Failed consistency check in  &
               set_output_defs: Hourly is ",i, "Nhourly is ",NHOURLY_OUT

          call CheckStop(errmsg)
        end if


   end do

  !/** Wanted dates for instantaneous values output: 
  !    specify months,days,hours for which full output is wanted.

     wanted_dates_inst(1) = date(-1,1,1,0,0)
     wanted_dates_inst(2) = date(-1,1,1,3,0)
     wanted_dates_inst(3) = date(-1,1,1,6,0)

 end subroutine set_output_defs
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module My_Outputs_ml

