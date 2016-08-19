!>_________________________________________________________<

  module  ChemGroups_ml
!-----------------------------------------------------------

  use ChemSpecs_tot_ml  ! => species indices
use OwnDataTypes_ml   ! => typ_sp
  implicit none
  private
! Assignment of groups from GenIn.species:
 public :: Init_ChemGroups

! ------- Gas/particle species ------------------

  integer, public, parameter ::  INDEX_DDEP_SS_GROUP = 1
  integer, public, target, save, dimension(3) :: &
             DDEP_SS_GROUP     = (/ SEASALT_F,SEASALT_C,SEASALT_G /)

  integer, public, parameter ::  INDEX_SOX_GROUP = 2
  integer, public, target, save, dimension(2) :: &
             SOX_GROUP     = (/ SO2,SO4 /)

  integer, public, parameter ::  INDEX_PMFINE_GROUP = 3
  integer, public, target, save, dimension(7) :: &
             PMFINE_GROUP     = (/ SO4,NO3_F,NH4_F,PPM25,PPM25_FIRE,SEASALT_F,DUST_NAT_F /)

  integer, public, parameter ::  INDEX_WDEP_OXN_GROUP = 4
  integer, public, target, save, dimension(4) :: &
             WDEP_OXN_GROUP     = (/ HNO3,HONO,NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_WDEP_PMCO_GROUP = 5
  integer, public, target, save, dimension(5) :: &
             WDEP_PMCO_GROUP     = (/ NO3_C,PPM_C,SEASALT_C,SEASALT_G,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_DUST_GROUP = 6
  integer, public, target, save, dimension(2) :: &
             DUST_GROUP     = (/ DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_DDEP_AOD_GROUP = 7
  integer, public, target, save, dimension(8) :: &
             DDEP_AOD_GROUP     = (/ SO4,NO3_F,NH4_F,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_WDEP_AOD_GROUP = 8
  integer, public, target, save, dimension(8) :: &
             WDEP_AOD_GROUP     = (/ SO4,NO3_F,NH4_F,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_DDEP_NOX_GROUP = 9
  integer, public, target, save, dimension(1) :: &
             DDEP_NOX_GROUP     = (/ NO2 /)

  integer, public, parameter ::  INDEX_WDEP_DUST_GROUP = 10
  integer, public, target, save, dimension(2) :: &
             WDEP_DUST_GROUP     = (/ DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_NOX_GROUP = 11
  integer, public, target, save, dimension(2) :: &
             NOX_GROUP     = (/ NO,NO2 /)

  integer, public, parameter ::  INDEX_SS_GROUP = 12
  integer, public, target, save, dimension(3) :: &
             SS_GROUP     = (/ SEASALT_F,SEASALT_C,SEASALT_G /)

  integer, public, parameter ::  INDEX_DDEP_DUST_GROUP = 13
  integer, public, target, save, dimension(2) :: &
             DDEP_DUST_GROUP     = (/ DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_PMCO_GROUP = 14
  integer, public, target, save, dimension(5) :: &
             PMCO_GROUP     = (/ NO3_C,PPM_C,SEASALT_C,SEASALT_G,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_DDEP_PMFINE_GROUP = 15
  integer, public, target, save, dimension(7) :: &
             DDEP_PMFINE_GROUP     = (/ SO4,NO3_F,NH4_F,PPM25,PPM25_FIRE,SEASALT_F,DUST_NAT_F /)

  integer, public, parameter ::  INDEX_RO2_GROUP = 16
  integer, public, target, save, dimension(13) :: &
             RO2_GROUP     = (/ HO2,CH3O2,C2H5O2,SECC4H9O2,ISRO2,ETRO2,PRRO2,OXYO2,MEKO2,MALO2,MVKO2,MACRO2,MACO3 /)

  integer, public, parameter ::  INDEX_ROOH_GROUP = 17
  integer, public, target, save, dimension(16) :: &
             ROOH_GROUP     = (/ CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,OXYO2H,MEKO2H,MALO2H,MVKO2H,MACROOH,MACO3H,ISRO2H,H2O2,CH3COO2H,ISONO3H,ISNIRH /)

  integer, public, parameter ::  INDEX_WDEP_ROOH_GROUP = 18
  integer, public, target, save, dimension(1) :: &
             WDEP_ROOH_GROUP     = (/ H2O2 /)

  integer, public, parameter ::  INDEX_WDEP_RDN_GROUP = 19
  integer, public, target, save, dimension(2) :: &
             WDEP_RDN_GROUP     = (/ NH3,NH4_F /)

  integer, public, parameter ::  INDEX_AOD_GROUP = 20
  integer, public, target, save, dimension(8) :: &
             AOD_GROUP     = (/ SO4,NO3_F,NH4_F,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_DDEP_SIA_GROUP = 21
  integer, public, target, save, dimension(4) :: &
             DDEP_SIA_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F /)

  integer, public, parameter ::  INDEX_WDEP_PM10_GROUP = 22
  integer, public, target, save, dimension(12) :: &
             WDEP_PM10_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F,PPM25,PPM25_FIRE,PPM_C,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_DDEP_OX_GROUP = 23
  integer, public, target, save, dimension(2) :: &
             DDEP_OX_GROUP     = (/ O3,NO2 /)

  integer, public, parameter ::  INDEX_DDEP_SOX_GROUP = 24
  integer, public, target, save, dimension(2) :: &
             DDEP_SOX_GROUP     = (/ SO2,SO4 /)

  integer, public, parameter ::  INDEX_OXN_GROUP = 25
  integer, public, target, save, dimension(13) :: &
             OXN_GROUP     = (/ NO,NO2,PAN,MPAN,NO3,N2O5,ISONO3,HNO3,HONO,ISNI,ISNIR,NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_TNO3_GROUP = 26
  integer, public, target, save, dimension(2) :: &
             TNO3_GROUP     = (/ NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_WDEP_SOX_GROUP = 27
  integer, public, target, save, dimension(2) :: &
             WDEP_SOX_GROUP     = (/ SO2,SO4 /)

  integer, public, parameter ::  INDEX_PM10_GROUP = 28
  integer, public, target, save, dimension(12) :: &
             PM10_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F,PPM25,PPM25_FIRE,PPM_C,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_OX_GROUP = 29
  integer, public, target, save, dimension(2) :: &
             OX_GROUP     = (/ O3,NO2 /)

  integer, public, parameter ::  INDEX_SIA_GROUP = 30
  integer, public, target, save, dimension(4) :: &
             SIA_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F /)

  integer, public, parameter ::  INDEX_DDEP_OXN_GROUP = 31
  integer, public, target, save, dimension(7) :: &
             DDEP_OXN_GROUP     = (/ NO2,PAN,MPAN,HNO3,HONO,NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_DDEP_ROOH_GROUP = 32
  integer, public, target, save, dimension(3) :: &
             DDEP_ROOH_GROUP     = (/ CH3O2H,C2H5OOH,H2O2 /)

  integer, public, parameter ::  INDEX_DDEP_PM10_GROUP = 33
  integer, public, target, save, dimension(12) :: &
             DDEP_PM10_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F,PPM25,PPM25_FIRE,PPM_C,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_WDEP_PMFINE_GROUP = 34
  integer, public, target, save, dimension(7) :: &
             WDEP_PMFINE_GROUP     = (/ SO4,NO3_F,NH4_F,PPM25,PPM25_FIRE,SEASALT_F,DUST_NAT_F /)

  integer, public, parameter ::  INDEX_DDEP_TNO3_GROUP = 35
  integer, public, target, save, dimension(2) :: &
             DDEP_TNO3_GROUP     = (/ NO3_F,NO3_C /)

  integer, public, parameter ::  INDEX_WDEP_SIA_GROUP = 36
  integer, public, target, save, dimension(4) :: &
             WDEP_SIA_GROUP     = (/ SO4,NO3_F,NO3_C,NH4_F /)

  integer, public, parameter ::  INDEX_DDEP_RDN_GROUP = 37
  integer, public, target, save, dimension(2) :: &
             DDEP_RDN_GROUP     = (/ NH3,NH4_F /)

  integer, public, parameter ::  INDEX_DDEP_PMCO_GROUP = 38
  integer, public, target, save, dimension(5) :: &
             DDEP_PMCO_GROUP     = (/ NO3_C,PPM_C,SEASALT_C,SEASALT_G,DUST_NAT_C /)

  integer, public, parameter ::  INDEX_BVOC_GROUP = 39
  integer, public, target, save, dimension(2) :: &
             BVOC_GROUP     = (/ C5H8,APINENE /)

  integer, public, parameter ::  INDEX_RDN_GROUP = 40
  integer, public, target, save, dimension(2) :: &
             RDN_GROUP     = (/ NH3,NH4_F /)

  integer, public, parameter ::  INDEX_WDEP_SS_GROUP = 41
  integer, public, target, save, dimension(3) :: &
             WDEP_SS_GROUP     = (/ SEASALT_F,SEASALT_C,SEASALT_G /)

  integer, public, parameter ::  INDEX_WDEP_TNO3_GROUP = 42
  integer, public, target, save, dimension(2) :: &
             WDEP_TNO3_GROUP     = (/ NO3_F,NO3_C /)


!GROUP ARRAY SIZE 42 MAXN 16 

  type, public :: gtype 
       character(len=20) :: name
       integer :: Ngroup
       integer, dimension(16) :: itot   ! indices from xn_tot arrays
  end type gtype

  type(gtype), public, parameter, dimension(42) :: &
       GROUP_ARRAY = (/ &
 gtype( "DDEP_SS", 3, (/ SEASALT_F,SEASALT_C,SEASALT_G,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "SOX", 2, (/ SO2,SO4,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PMFINE", 7, (/ SO4,NO3_F,NH4_F,PPM25,PPM25_FIRE,SEASALT_F,DUST_NAT_F,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_OXN", 4, (/ HNO3,HONO,NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PMCO", 5, (/ NO3_C,PPM_C,SEASALT_C,SEASALT_G,DUST_NAT_C,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DUST", 2, (/ DUST_NAT_F,DUST_NAT_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_AOD", 8, (/ SO4,NO3_F,NH4_F,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_AOD", 8, (/ SO4,NO3_F,NH4_F,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_NOX", 1, (/ NO2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_DUST", 2, (/ DUST_NAT_F,DUST_NAT_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "NOX", 2, (/ NO,NO2,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "SS", 3, (/ SEASALT_F,SEASALT_C,SEASALT_G,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_DUST", 2, (/ DUST_NAT_F,DUST_NAT_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PMCO", 5, (/ NO3_C,PPM_C,SEASALT_C,SEASALT_G,DUST_NAT_C,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PMFINE", 7, (/ SO4,NO3_F,NH4_F,PPM25,PPM25_FIRE,SEASALT_F,DUST_NAT_F,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "RO2", 13, (/ HO2,CH3O2,C2H5O2,SECC4H9O2,ISRO2,ETRO2,PRRO2,OXYO2,MEKO2,MALO2,MVKO2,MACRO2,MACO3,0,0,0 /) ) &
, gtype( "ROOH", 16, (/ CH3O2H,C2H5OOH,BURO2H,ETRO2H,PRRO2H,OXYO2H,MEKO2H,MALO2H,MVKO2H,MACROOH,MACO3H,ISRO2H,H2O2,CH3COO2H,ISONO3H,ISNIRH /) ) &
, gtype( "WDEP_ROOH", 1, (/ H2O2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_RDN", 2, (/ NH3,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "AOD", 8, (/ SO4,NO3_F,NH4_F,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_SIA", 4, (/ SO4,NO3_F,NO3_C,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_PM10", 12, (/ SO4,NO3_F,NO3_C,NH4_F,PPM25,PPM25_FIRE,PPM_C,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C,0,0,0,0 /) ) &
, gtype( "DDEP_OX", 2, (/ O3,NO2,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_SOX", 2, (/ SO2,SO4,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "OXN", 13, (/ NO,NO2,PAN,MPAN,NO3,N2O5,ISONO3,HNO3,HONO,ISNI,ISNIR,NO3_F,NO3_C,0,0,0 /) ) &
, gtype( "TNO3", 2, (/ NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_SOX", 2, (/ SO2,SO4,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "PM10", 12, (/ SO4,NO3_F,NO3_C,NH4_F,PPM25,PPM25_FIRE,PPM_C,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C,0,0,0,0 /) ) &
, gtype( "OX", 2, (/ O3,NO2,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "SIA", 4, (/ SO4,NO3_F,NO3_C,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_OXN", 7, (/ NO2,PAN,MPAN,HNO3,HONO,NO3_F,NO3_C,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_ROOH", 3, (/ CH3O2H,C2H5OOH,H2O2,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PM10", 12, (/ SO4,NO3_F,NO3_C,NH4_F,PPM25,PPM25_FIRE,PPM_C,SEASALT_F,SEASALT_C,SEASALT_G,DUST_NAT_F,DUST_NAT_C,0,0,0,0 /) ) &
, gtype( "WDEP_PMFINE", 7, (/ SO4,NO3_F,NH4_F,PPM25,PPM25_FIRE,SEASALT_F,DUST_NAT_F,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_TNO3", 2, (/ NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_SIA", 4, (/ SO4,NO3_F,NO3_C,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_RDN", 2, (/ NH3,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "DDEP_PMCO", 5, (/ NO3_C,PPM_C,SEASALT_C,SEASALT_G,DUST_NAT_C,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "BVOC", 2, (/ C5H8,APINENE,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "RDN", 2, (/ NH3,NH4_F,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_SS", 3, (/ SEASALT_F,SEASALT_C,SEASALT_G,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
, gtype( "WDEP_TNO3", 2, (/ NO3_F,NO3_C,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) ) &
   /)

! ------- Dry dep      species ------------------
  integer, public, parameter, dimension(7) :: &
               DDEP_OXNGROUP = (/ NO2,PAN,MPAN,HNO3,HONO,NO3_f,NO3_c /)
  integer, public, parameter, dimension(2) :: &
               DDEP_SOXGROUP = (/ SO2,SO4 /)
  integer, public, parameter, dimension(2) :: &
               DDEP_RDNGROUP = (/ NH3,NH4_f /)

  integer, public, parameter :: NMAX_DDEP = 7


! ------- Wet dep      species ------------------
  integer, public, parameter, dimension(4) :: &
               WDEP_OXNGROUP = (/ HNO3,HONO,NO3_f,NO3_c /)
  integer, public, parameter, dimension(2) :: &
               WDEP_SOXGROUP = (/ SO2,SO4 /)
  integer, public, parameter, dimension(3) :: &
               WDEP_SSALTGROUP = (/ SeaSalt_f,SeaSalt_c,SeaSalt_g /)
  integer, public, parameter, dimension(2) :: &
               WDEP_RDNGROUP = (/ NH3,NH4_f /)

  integer, public, parameter :: NMAX_WDEP = 4


! ------- RO2 Pool     species ------------------
  integer, public, parameter :: SIZE_RO2_POOL      = 1
  integer, public, parameter, dimension(1) :: &
     RO2_POOL      = (/ -99 /)
   type(typ_sp), dimension(42), public, save :: chemgroups


!-----------------------------------------------------------
  contains
 subroutine Init_ChemGroups()

   integer, dimension(:), pointer :: p

 p => DDEP_SS_GROUP 
 chemgroups(1) = typ_sp("DDEP_SS", p )

 p => SOX_GROUP 
 chemgroups(2) = typ_sp("SOX", p )

 p => PMFINE_GROUP 
 chemgroups(3) = typ_sp("PMFINE", p )

 p => WDEP_OXN_GROUP 
 chemgroups(4) = typ_sp("WDEP_OXN", p )

 p => WDEP_PMCO_GROUP 
 chemgroups(5) = typ_sp("WDEP_PMCO", p )

 p => DUST_GROUP 
 chemgroups(6) = typ_sp("DUST", p )

 p => DDEP_AOD_GROUP 
 chemgroups(7) = typ_sp("DDEP_AOD", p )

 p => WDEP_AOD_GROUP 
 chemgroups(8) = typ_sp("WDEP_AOD", p )

 p => DDEP_NOX_GROUP 
 chemgroups(9) = typ_sp("DDEP_NOX", p )

 p => WDEP_DUST_GROUP 
 chemgroups(10) = typ_sp("WDEP_DUST", p )

 p => NOX_GROUP 
 chemgroups(11) = typ_sp("NOX", p )

 p => SS_GROUP 
 chemgroups(12) = typ_sp("SS", p )

 p => DDEP_DUST_GROUP 
 chemgroups(13) = typ_sp("DDEP_DUST", p )

 p => PMCO_GROUP 
 chemgroups(14) = typ_sp("PMCO", p )

 p => DDEP_PMFINE_GROUP 
 chemgroups(15) = typ_sp("DDEP_PMFINE", p )

 p => RO2_GROUP 
 chemgroups(16) = typ_sp("RO2", p )

 p => ROOH_GROUP 
 chemgroups(17) = typ_sp("ROOH", p )

 p => WDEP_ROOH_GROUP 
 chemgroups(18) = typ_sp("WDEP_ROOH", p )

 p => WDEP_RDN_GROUP 
 chemgroups(19) = typ_sp("WDEP_RDN", p )

 p => AOD_GROUP 
 chemgroups(20) = typ_sp("AOD", p )

 p => DDEP_SIA_GROUP 
 chemgroups(21) = typ_sp("DDEP_SIA", p )

 p => WDEP_PM10_GROUP 
 chemgroups(22) = typ_sp("WDEP_PM10", p )

 p => DDEP_OX_GROUP 
 chemgroups(23) = typ_sp("DDEP_OX", p )

 p => DDEP_SOX_GROUP 
 chemgroups(24) = typ_sp("DDEP_SOX", p )

 p => OXN_GROUP 
 chemgroups(25) = typ_sp("OXN", p )

 p => TNO3_GROUP 
 chemgroups(26) = typ_sp("TNO3", p )

 p => WDEP_SOX_GROUP 
 chemgroups(27) = typ_sp("WDEP_SOX", p )

 p => PM10_GROUP 
 chemgroups(28) = typ_sp("PM10", p )

 p => OX_GROUP 
 chemgroups(29) = typ_sp("OX", p )

 p => SIA_GROUP 
 chemgroups(30) = typ_sp("SIA", p )

 p => DDEP_OXN_GROUP 
 chemgroups(31) = typ_sp("DDEP_OXN", p )

 p => DDEP_ROOH_GROUP 
 chemgroups(32) = typ_sp("DDEP_ROOH", p )

 p => DDEP_PM10_GROUP 
 chemgroups(33) = typ_sp("DDEP_PM10", p )

 p => WDEP_PMFINE_GROUP 
 chemgroups(34) = typ_sp("WDEP_PMFINE", p )

 p => DDEP_TNO3_GROUP 
 chemgroups(35) = typ_sp("DDEP_TNO3", p )

 p => WDEP_SIA_GROUP 
 chemgroups(36) = typ_sp("WDEP_SIA", p )

 p => DDEP_RDN_GROUP 
 chemgroups(37) = typ_sp("DDEP_RDN", p )

 p => DDEP_PMCO_GROUP 
 chemgroups(38) = typ_sp("DDEP_PMCO", p )

 p => BVOC_GROUP 
 chemgroups(39) = typ_sp("BVOC", p )

 p => RDN_GROUP 
 chemgroups(40) = typ_sp("RDN", p )

 p => WDEP_SS_GROUP 
 chemgroups(41) = typ_sp("WDEP_SS", p )

 p => WDEP_TNO3_GROUP 
 chemgroups(42) = typ_sp("WDEP_TNO3", p )

   nullify(p)


 end subroutine Init_ChemGroups
 !-----------------------------------------------------------


 end module ChemGroups_ml
 !-----------------------------------------------------------
