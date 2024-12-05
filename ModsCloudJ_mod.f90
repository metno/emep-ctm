!------------------------------------------------------------------------------
!    'cmn_fjx_mod.f90'  for fast-JX code v 7.3+ (prather 2/15)
!         note that module and file begin with 'cmn_"
!         small changes: LREF=51 instead of hardwired, JVMAP replaces JMAP
!         MASFAC param added, also cloud params - see below
!------------------------------------------------------------------------------
!
! NB - ALL of these common variables are set paramters,
!    They are NOT to be used as variables for a local solution
!    Thus this entire set is 'in' only after it is initialized
!-----------------------------------------------------------------------
! !DESCRIPTION: FJX_CMN contains fast-JX variables
!
!
! !INTERFACE:
!
      MODULE FJX_CMN_MOD

        ! use Config_module,         only : KMAX_MID

        implicit none
        public
  
  !-----------------------------------------------------------------------
  
        ! JXL_: vertical(levels) dim for J-values computed within fast-JX
        integer, parameter ::  JXL_=100, JXL1_=JXL_+1
        ! JXL2_: 2*JXL_ + 2 = mx no.levels in basic FJX grid (mid-level)
        integer, parameter ::  JXL2_=2*JXL_+2
        ! WX_  = dim = no. of wavelengths in input file
        integer, parameter ::  WX_=18
        ! X_   = dim = max no. of X-section data sets (input data)  
        integer, parameter ::  X_=201
        ! A_   = dim = no. of Aerosol/cloud Mie sets (input data)
        integer, parameter ::  A_=40
        ! C_   = dim = no. of cld-data sets (input data)
        integer, parameter ::  C_=16
        ! W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only  
        integer, parameter ::  W_=12    ! W_= 8, 12 or 18
        ! N_  = no. of levels in Mie scattering arrays
        !     = 2*NC+1 = 4*(L_+1) + 1`+ 2*sum(JADDLV)
        integer, parameter ::  N_=601
        ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
        integer, parameter ::  M_=4
        ! M2_ = 2*M_ = 8, replaces MFIT
        integer, parameter ::  M2_=2*M_
  
        character*25, dimension(8), parameter :: TITCLD =  &
        ['clear sky - no clouds    ', &
         'avg cloud cover          ', &
         'avg cloud cover^3/2      ', &
         'ICAs - avg direct beam   ', &
         'ICAs - random N ICAs     ', &
         'QCAs - midpt of bins     ', &
         'QCAs - avg clouds in bins', &
         'ICAs - use all ICAs***   ']
  
  !-----------------------------------------------------------------------
        ! 4 Gauss pts = 8-stream
        real*8, DIMENSION(M_), parameter  ::  &
                           EMU = [.06943184420297d0, .33000947820757d0, &
                                  .66999052179243d0, .93056815579703d0]
        real*8, DIMENSION(M_), parameter  :: & 
                           WT  = [.17392742256873d0, .32607257743127d0, &
                                  .32607257743127d0, .17392742256873d0]
  !-----------------------------------------------------------------------
  
        ! MASFAC: Conversion factor for pressure to column density (# molec /cm2)
        real*8, parameter   ::  &
                           MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)
        ! ZZHT: scale height (cm) used above top of CTM ZHL(LPAR+1)
        real*8, parameter   :: ZZHT = 5.d5
        ! RAD: Radius of Earth (cm)
        real*8, parameter   :: RAD = 6375.d5
        ! ATAU: heating rate (factor increase from one layer to the next)
        real*8, parameter   :: ATAU = 1.120d0
        ! ATAU0: minimum heating rate
        real*8, parameter   :: ATAU0 = 0.010d0
        ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
        integer, parameter  :: JTAUMX = (N_ - 4*JXL_)/2
  
  !---- Variables in file 'FJX_spec.dat' (RD_XXX)
  
        ! WBIN: Boundaries of wavelength bins
        real*8  WBIN(WX_+1)
        ! WL: Centres of wavelength bins - 'effective wavelength'
        real*8  WL(WX_)
        ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
        ! FW: Solar flux in W/m2
        ! FP: PAR quantum action spectrum
        real*8  FL(WX_),FW(WX_),FP(WX_)
  
  
        real*8  QO2(WX_,3)   ! QO2: O2 cross-sections
        real*8  QO3(WX_,3)   ! QO3: O3 cross-sections
        real*8  Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield
  
        ! QQQ: Supplied cross sections in each wavelength bin (cm2)
        real*8  QQQ(WX_,3,X_)
        ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
        real*8  QRAYL(WX_+1)
        ! TQQ: Temperature for supplied cross sections
        real*8  TQQ(3,X_)
        ! LQQ = 1, 2, or 3 to determine interpolation with T or P
        integer LQQ(X_)
  
        ! TITLEJX: Title (short & long) for supplied cross sections, from 'FJX_spec.dat'
        CHARACTER*6  TITLEJX(X_)
        CHARACTER*16 TITLEJL(X_)
        ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
        CHARACTER*1  SQQ(X_)
  
  !---- Variables in file 'FJX_scat-aer.dat' (RD_MIE)
  
        ! TITLAA: Aerosol Mie Titles
        character*12  TITLAA(A_)
        ! QAA: Aerosol scattering phase functions
        real*8  QAA(5,A_)
        ! WAA: 5 Wavelengths for the supplied phase functions
        real*8  WAA(5,A_)
        ! PAA: Phase function: first 8 terms of expansion
        real*8  PAA(8,5,A_)
        ! RAA: Effective radius associated with aerosol type
        real*8  RAA(A_)
        ! SAA: Single scattering albedo
        real*8  SAA(5,A_)
        ! DAA: density (g/cm^3)
        real*8  DAA(A_)
        ! NAA: Number of categories for scattering phase functions
        integer NAA
  
  !---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)
  
        ! TITLCC: Cloud type titles
        character*12  TITLCC(C_)
        ! QCC: Cloud scattering phase functions
        real*8  QCC(5,C_)
        ! WCC: 5 Wavelengths for supplied phase functions
        real*8  WCC(5,C_)
        ! PCC: Phase function: first 8 terms of expansion
        real*8  PCC(8,5,C_)
        ! RCC: Effective radius associated with cloud type
        real*8  RCC(C_)
        ! SCC: Single scattering albedo
        real*8  SCC(5,C_)
        ! DCC: density (g/cm^3)
        real*8  DCC(C_)
        ! NCC: Number of categories for cloud scattering phase functions
        integer NCC
  
  !---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)
  
        ! WMM: U Michigan aerosol wavelengths
        real*8  WMM(6)
        ! UMAER: U Michigan aerosol data sets
        real*8  UMAER(3,6,21,33)
  
  !---- Variables in file 'atmos_std.dat' (RD_PROF)
  
        integer, parameter ::  LREF=51   ! layer dim. in reference profiles
        integer, parameter ::  JREF=18   ! latitude dim. in reference profiles
  
        ! T and O3 reference profiles
        real*8, DIMENSION(LREF,JREF,12) :: TREF, OREF
  
        integer NJX,NW1,NW2
  
  !-----------------NEW for FJX72 parameters for cloud grid now here------
  !
        integer, parameter :: &
            !   LPAR= 45, LWEPAR=20  &   !this can be set by CTM code. LWEPAR = number of lvls with cloud calcs
            !  ,L_=LPAR, L1_=L_+1 &   ! L_ = number of CTM layers
            !  ,L2_=2*L_+2 &        ! no. levels in the Fast-JX grid that
            !              ! includes both layer edges and layer mid-points
            !  ,JVL_=LPAR &  ! vertical(levels) dim for J-values sent to CTM
             JVN_=101 &   ! max no. of J-values
             ,AN_=25       ! # FJX aerosols in layer (needs NDX for each)

        integer, save :: LPAR, LWEPAR  &   !this can be set by CTM code. LWEPAR = number of lvls with cloud calcs
        ,L_, L1_ &   ! L_ = number of CTM layers
        ,L2_ &       ! no. levels in the Fast-JX grid that
                     ! includes both layer edges and layer mid-points
        ,JVL_        ! vertical(levels) dim for J-values sent to CTM
  
  !-----------------------------------------------------------------------
        ! variables used to map fast-JX J's onto CTM J's
  !-----------------------------------------------------------------------
  
        integer :: debug
        real*8  JFACTA(JVN_)  ! multiplication factor for fast-JX calculated J
        integer JIND(JVN_)    ! index arrays that map Jvalue(j) onto rates
        integer NRATJ         ! number of Photolysis reactions in CTM chemistry,
                              ! derived here NRATJ must be .le. JVN_
        character*6 JVMAP(JVN_)   ! label of J-value used to match w/FJX J's
        character*50 JLABEL(JVN_) ! label of J-value used in the chem model
  
  
  ! Cloud Cover parameters
  !-----------------------------------------------------------------------
  ! NB CBIN_ was set at 20, but with NRG=6 groups, 10 gives a more reasonable number of ICAs
        integer, parameter :: CBIN_ = 10     ! # of quantized cloud fration bins
        integer, parameter :: ICA_ = 20000   ! Max # of indep colm atmospheres
        integer, parameter :: NQD_ = 4       ! # of Quad Colm Atmos cloud fraction bins (4)
  
        real*8, parameter ::  CPI    = 3.141592653589793d0
        real*8, parameter ::  C2PI   = 2.d0*CPI
        real*8, parameter ::  CPI180 = CPI/180.d0
        real*8, parameter ::  G0     = 9.80665d0
        real*8,  parameter :: G100 = 100.d0/G0
  !-------data to set up the random number sequence for use in cloud-JX
        integer, parameter :: NRAN_ = 10007  ! dimension for random number
        real*4   RAN4(NRAN_)      ! Random number set
  
  
        END MODULE FJX_CMN_MOD

!------------------------------------------------------------------------------
!     'fjx_sub_mod.F90'  for fast-JX code v7.3+                                     !
!------------------------------------------------------------------------------
!
! !MODULE: FJX
!
! !DESCRIPTION: JX version 7.2  (12/13) consistent with 7.1 data and results
!          variables in call to PHOTO_JX are same as in 7.1,
!          but a logical(out) LDARK is added to to count the number of J calcs
!          Also works with new ver 7.3+ data sets if _init_ is updated.
!
! !INTERFACE:
!
      MODULE FJX_SUB_MOD
!
! !USES:
!
      USE FJX_CMN_MOD
      use Config_module, only: cloudjx_initf, USES, MasterProc

      IMPLICIT NONE
!
! !PUBLIC SUBROUTINES:
!
      PUBLIC  :: SOLAR_JX, ACLIM_FJX, JP_ATM0, PHOTO_JX, EXITC


      CONTAINS

!-----------------------------------------------------------------------
      subroutine SOLAR_JX(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!-----------------------------------------------------------------------
!     GMTIME = UT for when J-values are wanted
!           (for implicit solver this is at the end of the time step)
!     NDAY   = integer day of the year (used for solar lat and declin)
!     YGRDJ  = laitude (radians) for grid (I,J)
!     XGDRI  = longitude (radians) for grid (I,J)
!
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in) ::   GMTIME,YGRDJ,XGRDI
      integer, intent(in) ::  NDAY
      real*8, intent(out) ::  SZA,COSSZA,SOLFX
!
      real*8  LOCT
      real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
!
      SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*CPI180)
      SOLDEK = asin(SINDEC)
      COSDEC = cos(SOLDEK)
      SINLAT = sin(YGRDJ)
      SOLLAT = asin(SINLAT)
      COSLAT = cos(SOLLAT)
!
      LOCT   = (((GMTIME)*15.d0)-180.d0)*CPI180 + XGRDI
      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      SZA    = acos(COSSZA)/CPI180
!
      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*C2PI/365.d0))
!
      END SUBROUTINE SOLAR_JX


!-----------------------------------------------------------------------
      subroutine ACLIM_FJX (YLATD, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1U)
!-----------------------------------------------------------------------
!  Load fast-JX climatology for latitude & month given pressure grid
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in)  :: YLATD
      integer, intent(in) :: MONTH, L1U
      real*8, intent(in),  dimension(L1U+1) :: PPP
      real*8, intent(out), dimension(L1U+1) :: ZZZ
      real*8, intent(out), dimension(L1U) :: TTT,DDD,OOO

      real*8, dimension(LREF) :: OREF2,TREF2
      real*8, dimension(LREF+1) :: PSTD
      integer  K, L, M, N
      real*8   DLOGP,F0,T0,PB,PC,XC,SCALEH

!  Select appropriate month
      M = max(1,min(12,MONTH))
!  Select appropriate latitudinal profiles
      N = max(1, min(18, (int(YLATD+99)/10 )))
      do K = 1,LREF
        OREF2(K) = OREF(K,N,M)
        TREF2(K) = TREF(K,N,M)
      enddo

!  Apportion O3 and T on supplied climatology z levels onto CTM levels +1
!  with mass (pressure) weighting, assuming constant mixing ratio and
!  temperature half a layer on either side of the point supplied.
!   PPP(L=1:L1_)=edge-pressure of CTM layer, PPP(L1_+1)=0 (top-of-atmos)
!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
!     MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

!  Set up pressure levels for O3/T climatology - assume that value
!  given for each 2 km z* level applies from 1 km below to 1 km above,
!  so select pressures at these boundaries. Surface level values at
!  1000 mb are assumed to extend down to the actual PSURF (if > 1000)
           PSTD(1) = max(PPP(1),1000.d0)
           PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0)
           DLOGP   = 10.d0**(-2.d0/16.d0)
      do K = 3,LREF
        PSTD(K) = PSTD(K-1)*DLOGP
      enddo
        PSTD(52)  = 0.d0
      do L = 1,L1U
        F0 = 0.d0
        T0 = 0.d0
        do K = 1,LREF
          PC   = min(PPP(L),PSTD(K))
          PB   = max(PPP(L+1),PSTD(K+1))
          if (PC .gt. PB) then
            XC = (PC-PB)/(PPP(L)-PPP(L+1))
            F0 = F0 + OREF2(K)*XC
            T0 = T0 + TREF2(K)*XC
          endif
        enddo
        TTT(L)  = T0
        DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
        OOO(L) = F0*1.d-6*DDD(L)
      enddo

!  Calculate effective altitudes using scale height at each level
        ZZZ(1) = 0.d0
      do L = 1,L1U-1
        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
        ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH )
      enddo
        ZZZ(L1U+1) = ZZZ(L1U) + ZZHT

      END SUBROUTINE ACLIM_FJX

!<<<<<<<<<<<<<<<<<<<end CTM-fastJX special call subroutines<<<<<<<<<<<<<




!<<<<<<<<<<<<<<<<<<<<<<<<begin fastJX core code<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<all outside calls go through PHOTO_JX<<<<<<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
!  fastJX version 7.2 (f90) - Prather notes (Jan 2013)

!---calculates cloud optical depth in FJX-72 given water path & R-eff
!---assumes that clouds are 100% if in layer
!      IWP = ice water path (in layer, in cloud) in g/m**2
!      LWP = liquid water path (in layer, in cloud) in g/m**2
!      REFFI = effective radius of ice cloud (microns)
!      REFFL = effective radius of liquid cloud (microns)
!   ODL = LWP * 0.75 * 2.10 / REFFL
!   ODI = IWP * 0.75 * 2.00 / (REFFI * 0.917)

!>>>R-effective determined by main code, not FJX
!   REFF determined by user - some recommendations below (from Neu & Prather)
!          REFFI is function of ice water content IWC (g/m3, 0.0001 to 0.1)
!             IWC = IWP / delta-Z (of layer in m, approx OK)
! prefer Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)
!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR))
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)
!>>>indices for cloud scattering determined by FJX core, not main code.
!   ice clouds pick heagonal or irregular:
!          NDXI = 6  ! ice hexag (cold)
!             if (TCLD .ge. 233.15) then
!          NDXI = 7  ! ice irreg
!             endif
!   liquid clouds scaled to R-eff 3 - 6 - 12 - 20 microns:
!          NDXC = 1
!            do I=2,4
!             if (REFFL .gt. 0.5*(RCC(I-1)+RCC(I))) then
!          NDXC = I
!             endif
!            enddo

!-----------------------------------------------------------------------
      SUBROUTINE PHOTO_JX (U0,SZA,REFLB,SOLF, LPRTJ, PPP,ZZZ,TTT,DDD,  &
                           RRR,OOO, LWP,IWP,REFFL,REFFI,               &
                           AERSP,NDXAER, L1U,ANU, VALJXX,NJXU, LDARK)

!-----------------------------------------------------------------------
!
!  PHOTO_JX is the gateway to fast-JX calculations:
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-fo-year for sun distance, SZA (not lat or long)
!-----------------------------------------------------------------------
      implicit none

!---calling sequence variables
      integer, intent(in)                    :: L1U,ANU,NJXU
      real*8,  intent(in)                    :: U0,SZA,REFLB,SOLF
      logical, intent(in)                    :: LPRTJ
      real*8,  intent(in), dimension(L1U+1)  :: PPP,ZZZ
      real*8,  intent(in), dimension(L1U  )  :: TTT,DDD,RRR,OOO,  &
                                                LWP,IWP,REFFL,REFFI
      real*8,  intent(in), dimension(L1U,ANU):: AERSP
      integer, intent(in), dimension(L1U,ANU):: NDXAER

! reports out the JX J-values, upper level program converts to CTM chemistry J's
      real*8, intent(out), dimension(L1U-1,NJXU)::  VALJXX
      logical, intent(out)                   :: LDARK

!-----------------------------------------------------------------------
!--------key LOCAL atmospheric data needed to solve plane-parallel J----
!-----these are dimensioned JXL_, and must have JXL_ .ge. L_
      real*8, dimension(JXL1_+1) :: TTJ,DDJ,OOJ,ZZJ
      real*8, dimension(JXL1_+1) :: PPJ,RELH
      integer,dimension(JXL2_+1) :: JXTRA
      real*8, dimension(W_)      :: FJTOP,FJBOT,FSBOT,FLXD0,RFL
      real*8, dimension(JXL_, W_)   ::  AVGF, FJFLX
      real*8, dimension(JXL1_,W_)   ::  DTAUX, FLXD
      real*8, dimension(8,JXL1_,W_) ::  POMEGAX
      real*8, dimension(JXL1_)      ::  DTAU600
      real*8, dimension(8,JXL1_)    ::  POMG600
      real*8, dimension(W_,JXL1_)   ::  FFX
      real*8, dimension(W_,8)     ::  FFXNET
!---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
      real*8  FLXJ(JXL1_),FFX0,FXBOT,FABOT
      real*8             :: ODABS,ODRAY,ODI,ODL
      real*8             :: RFLECT,FREFS,FREFL,FREFI
      real*8             :: AMF2(2*JXL1_+1,2*JXL1_+1)
!------------key SCATTERING arrays for clouds+aerosols------------------
      real*8             :: OPTX(5),SSAX(5),SLEGX(8,5)
      real*8             :: OD(5,JXL1_),SSA(5,JXL1_),SLEG(8,5,JXL1_)
      real*8             :: OD600(JXL1_)
      real*8             :: PATH,RH,XTINCT
!------------key arrays AFTER solving for J's---------------------------
      real*8  FFF(W_,JXL_),VALJ(X_)
      real*8  FLXUP(W_),FLXDN(W_),DIRUP(W_),DIRDN(W_)
      real*8  VALJL(JXL_,X_) !2-D array of J_s returned by JRATET

      integer  LU,L2U, I,J,K,L,M,KMIE,KW,NAER,NDXI,NDXL, RATIO(W_)
      real*8   XQO3,XQO2,DTAUC,WAVE, TTTX
!-----------------------------------------------------------------------
      real*8   SWHEAT, SKPERD

      if (L1U .gt. JXL1_) then
        call EXITC(' PHOTO_JX: not enough levels in JX')
      endif

      LU = L1U - 1
      L2U = LU + LU + 2

      VALJXX(:,:) = 0.d0
      FFF(:,:) = 0.d0
      FREFI = 0.d0
      FREFL = 0.d0
      FREFS = 0.d0

!---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
!                        or         99.0                 80 km
      if (SZA .gt. 98.d0) then
        LDARK = .true.
        return
      else
        LDARK = .false.
      endif

!---load the amtospheric column data
      do L = 1,L1U
        PPJ(L) = PPP(L)
        TTJ(L) = TTT(L)
        DDJ(L) = DDD(L)
        OOJ(L) = OOO(L)
      enddo
        PPJ(L1U+1) = 0.d0

!---calculate spherical weighting functions (AMF: Air Mass Factor)
      do L = 1,L1U+1
        ZZJ(L) = ZZZ(L)
      enddo

        RFLECT = REFLB

!-----------------------------------------------------------------------
      call SPHERE2 (U0,RAD,ZZJ,ZZHT,AMF2, L1U,JXL1_)
!-----------------------------------------------------------------------

!---calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8))
!---  at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols
      do L = 1,L1U
        do K=1,5
          OD(K,L)  = 0.d0
          SSA(K,L) = 0.d0
          do I=1,8
            SLEG(I,K,L) = 0.d0
          enddo
        enddo
      enddo

      if (LPRTJ) then
      write(6,*)'fast-JX-(7.2)---PHOTO_JX internal print: Clouds--'
      write(6,*) '  L, P1,P2, WPath, Reff, OD, Ndx'
      endif

      do L = 1,L1U
!---Liquid Water Cloud
        if (LWP(L) .gt. 1.d-5 .and. REFFL(L) .gt. 0.1d0) then
!---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
          ODL = LWP(L) * 0.75d0 * 2.1d0 / REFFL(L)
          NDXL = 1
          do I = 2,4
            if (REFFL(L) .gt. 0.5*(RCC(I-1)+RCC(I))) then
              NDXL = I
            endif
          enddo
          call OPTICL (OPTX,SSAX,SLEGX,  ODL,NDXL)
          do K = 1,5
            OD(K,L)  = OD(K,L)  + OPTX(K)
            SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
            do I = 1,8
              SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
            enddo
          enddo
!>>>diagnostic print of cloud data:
          if (LPRTJ) then
            write(6,'(a,i3,2f8.2,f8.4,f8.2,f8.4,i4)') &
                  'Liq Cld',L,PPP(L),PPP(L+1),LWP(L),REFFL(L),ODL,NDXL
          endif
        endif

!---Ice Water Cloud
        if (IWP(L) .gt. 1.d-5 .and. REFFI(L) .gt. 0.1d0) then
          ODI = IWP(L) * 0.75d0 * 2.0d0 / (REFFI(L) * 0.917d0)
          if (TTT(L) .ge. 233.15d0) then
            NDXI = 7  ! ice irreg
          else
            NDXI = 6  ! ice hexag (cold)
          endif
          call OPTICL (OPTX,SSAX,SLEGX,  ODI,NDXI)
          do K=1,5
            OD(K,L)  = OD(K,L)  + OPTX(K)
            SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
            do I=1,8
              SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
            enddo
          enddo
!>>>diagnostic print of cloud data:
          if (LPRTJ) then
             write(6,'(a,i3,2f8.2,f8.4,f8.2,f8.4,i4)') &
                  'Ice Cld',L,PPP(L),PPP(L+1),IWP(L),REFFI(L),ODI,NDXI
          endif
        endif

!---aerosols in layer: check aerosol index
!---this uses data from climatology OR from current CTM (STT of aerosols)

!---FIND useful way to sum over different aerosol types!
          RH = RRR(L)
        do M = 1,ANU
          NAER = NDXAER(L,M)
          PATH = AERSP(L,M)

!---subroutines OPTICA & OPTICM return the same information:
!---  optical depth (OPTX), single-scat albedo (SSAX) and phase fn (SLEGX(8))
!---subs have slightly different inputs:
!---  PATH is the g/m2 in the layer, NAER in the cloud/aerosol index
!---  UMich aerosols use relative humidity (RH)
          if (PATH .gt. 0.d0) then
            if (NAER .gt.0) then
              call OPTICA (OPTX,SSAX,SLEGX,  PATH,RH, NAER)
            else
              call OPTICM (OPTX,SSAX,SLEGX,  PATH,RH,-NAER)
            endif
            do K = 1,5
              OD(K,L)  = OD(K,L)  + OPTX(K)
              SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
              do I = 1,8
                SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
              enddo
            enddo
          endif
        enddo

        do K = 1,5
          if (OD(K,L) .gt. 0.d0) then
            SSA(K,L) = SSA(K,L)/OD(K,L)
            do I = 1,8
              SLEG(I,K,L) = SLEG(I,K,L)/OD(K,L)
            enddo
          endif
        enddo

!---Include aerosol with cloud OD at 600 nm to determine added layers:
        OD600(L) = OD(4,L)

      enddo   ! end of 'do L = 1,L1U'

!---when combining with Rayleigh and O2-O3 abs, remember the SSA and
!---  phase fn SLEG are weighted by OD and OD*SSA, respectively.
!---Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add
!       additonal levels at top of clouds (now uses log spacing)

!-----------------------------------------------------------------------
      call EXTRAL(OD600,L1U,L2U,N_,JTAUMX,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------

!---set surface reflectance
      RFL(:) = max(0.d0,min(1.d0,RFLECT))

!-----------------------------------------------------------------------
!---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
!-----------------------------------------------------------------------

      do K = 1,W_
      if (FL(K) .gt. 1.d0) then
        WAVE = WL(K)
!---Pick nearest Mie wavelength to get scattering properites------------
                               KMIE=1  ! use 200 nm prop for <255 nm
        if( WAVE .gt. 255.d0 ) KMIE=2  ! use 300 nm prop for 255-355 nm
        if( WAVE .gt. 355.d0 ) KMIE=3  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500.d0 ) KMIE=4
        if( WAVE .gt. 800.d0 ) KMIE=5

!---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
!---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)
        do L = 1,L1U
          TTTX     = TTJ(L)
          call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1) &
                        ,QO2(K,2),TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2) &
                        ,QO3(K,2),TQQ(3,2),QO3(K,3), LQQ(2))
          ODABS = XQO3*OOJ(L) + XQO2*DDJ(L)*0.20948d0
          ODRAY = DDJ(L)*QRAYL(K)

          DTAUX(L,K) = OD(KMIE,L) + ODABS + ODRAY

          do I = 1,8
            POMEGAX(I,L,K) = SLEG(I,KMIE,L)*OD(KMIE,L)
          enddo
          POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0d0*ODRAY
          POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5d0*ODRAY
          do I = 1,8
            POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K)
          enddo
        enddo

      endif
      enddo

!-----------------------------------------------------------------------

      call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
              AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)

!-----------------------------------------------------------------------

         FLXUP(:) = 0.d0
         DIRUP(:) = 0.d0
         FLXDN(:) = 0.d0
         DIRDN(:) = 0.d0
         FLXJ(:) = 0.d0
         FFX(:,:) = 0.d0
         FFXNET(:,:) = 0.d0

      do K = 1,W_
      if (FL(K) .gt. 1.d0) then

!direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
!----     also at bottom (DN), does not include diffuse reflected flux.
        FLXUP(K) =  FJTOP(K)
        DIRUP(K) = -FLXD0(K)
        FLXDN(K) = -FJBOT(K)
        DIRDN(K) = -FSBOT(K)

        do L = 1,LU
          FFF(K,L) = SOLF*FL(K)*AVGF(L,K)
        enddo
        FREFI = FREFI + SOLF*FL(K)*FLXD0(K)/WL(K)
        FREFL = FREFL + SOLF*FL(K)*FJTOP(K)/WL(K)
        FREFS = FREFS + SOLF*FL(K)/WL(K)

!---for each wavelength calculate the flux budget/heating rates:
!  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L))]
!            but for spherical atmosphere!
!  FJFLX(L) = diffuse flux across top of layer L

!---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
!---     need special fix at top and bottom:
!---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
        FABOT = (1.d0-RFL(K))*(FJBOT(K)+FSBOT(K))
        FXBOT = -FJBOT(K) + RFL(K)*(FJBOT(K)+FSBOT(K))
        FLXJ(1) = FJFLX(1,K) - FXBOT
        do L = 2,LU
          FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K)
        enddo
        FLXJ(LU+1) = FJTOP(K) - FJFLX(LU,K)
!---calculate net flux deposited in each CTM layer (direct & diffuse):
        FFX0 = 0.d0
        do L=1,L1U
          FFX(K,L) = FLXD(L,K) - FLXJ(L)
          FFX0 = FFX0 + FFX(K,L)
        enddo

!  NB: the radiation level ABOVE the top CTM level is included in these budgets
!      these are the flux budget/heating terms for the column:
!  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
!  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)
!  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
!  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos
!  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos
!  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
!       these are surface fluxes to compare direct vs. diffuse:
!  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
!  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)

        FFXNET(K,1) = FLXD0(K)
        FFXNET(K,2) = FSBOT(K)
        FFXNET(K,3) = FLXD0(K) + FSBOT(K)
        FFXNET(K,4) = FJTOP(K)
        FFXNET(K,5) = FFX0
        FFXNET(K,6) = FABOT
        FFXNET(K,7) = FSBOT(K)
        FFXNET(K,8) = FJBOT(K)

!-----------------------------------------------------------------------
      endif
      enddo       ! end loop over wavelength K
!-----------------------------------------------------------------------
      FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
      FREFI = FREFI/FREFS

!---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)

!-----------------------------------------------------------------------
      call JRATET(PPJ,TTJ,FFF, VALJXX, LU,NJXU)
!-----------------------------------------------------------------------

!---mapping J-values from fast-JX onto CTM chemistry is done in main

!-----------------------------------------------------------------------
      if (LPRTJ) then
!---diagnostics below are NOT returned to the CTM code
       write(6,*)'fast-JX-(7.2)---PHOTO_JX internal print: Atmosphere--'
!---used last called values of DTAUX and POMEGAX, should be 600 nm
       do L=1,L1U
         DTAU600(L) = DTAUX(L,W_)
         do I=1,8
           POMG600(I,L) = POMEGAX(I,L,W_)
         enddo
       enddo

       call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)

!---PRINT SUMMARY of mean intensity, flux, heating rates:
       write(6,*)
       write(6,*)'fast-JX(7.2)---PHOTO_JX internal print: Mean Intens--'
       write(6,'(a,5f10.4)') &
       ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/', &
        RFLECT,SZA,U0,FREFI,FREFL

       write(6,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
       write(6,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
       write(6,'(a)') ' ---- 100000=Fsolar   MEAN INTENSITY per wvl bin'
         RATIO(:) = 0.d0
       do L = LU,1,-1
        do K=NW1,NW2
         if (FL(K) .gt. 1.d0) then
         RATIO(K) = (1.d5*FFF(K,L)/FL(K))
         endif
        enddo
         write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
       enddo

       write(6,*)
       write(6,*)'fast-JX(7.2)---PHOTO_JX internal print: Net Fluxes---'
       write(6,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
       write(6,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)

!  NB: the radiation level ABOVE the top CTM level is included in these budgets
!      these are the flux budget/heating terms for the column:
!  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
!  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)
!  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
!  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos
!  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos
!  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
!       these are surface fluxes to compare direct vs. diffuse:
!  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
!  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)


!      write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
!      write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
       write(6,*) ' ---NET FLUXES--- solar only to 850 nm'
       write(6,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3), K=NW2,NW1,-1)
       write(6,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4), K=NW2,NW1,-1)
       write(6,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5), K=NW2,NW1,-1)
       write(6,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6), K=NW2,NW1,-1)
       write(6,*) ' ---SRF FLUXES--- '
       write(6,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7), K=NW2,NW1,-1)
       write(6,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8), K=NW2,NW1,-1)

           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,3)*FW(K)
         enddo
        write(6,'(a11,f10.4)') ' sol TOTAL ',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,4)*FW(K)
         enddo
        write(6,'(a11,f10.4)') ' dif outtop',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,5)*FW(K)
         enddo
        write(6,'(a11,f10.4)') ' abs in atm',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,6)*FW(K)
         enddo
        write(6,'(a11,f10.4)') ' abs at srf',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,7)*FL(K)*FP(K)
         enddo
        write(6,'(a11,1p,e12.4)') ' PAR direct',SWHEAT
           SWHEAT = 0.d0
         do K=NW1,NW2
           SWHEAT = SWHEAT + FFXNET(K,8)*FL(K)*FP(K)
         enddo
        write(6,'(a11,1p,e12.4)') ' PAR diffus',SWHEAT

       write(6,'(2a)') ' absorb solar: W/m2/kyr K/day by wavel 10000=Fsolar', &
       '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
       do L = LU,1,-1
           SWHEAT = 0.d0
         do K=NW1,NW2
           RATIO(K) = 1.d5*FFX(K,L)
           SWHEAT = SWHEAT + FFX(K,L)*FW(K)
         enddo
           SKPERD = SWHEAT*86400.d0/(1005.d0*(PPP(L)-PPP(L+1))*100.d0/9.8d0)
         write(6,'(i9,2x,2f8.4,18i8)') L,SWHEAT,SKPERD,(RATIO(K),K=NW2,NW1,-1)
       enddo

       write(6,'(a)')
       write(6,'(a)') ' fast-JX (7.2)----J-values----'
       write(6,'(1x,a,72(a6,3x))') 'L=  ',(TITLEJX(K), K=1,NJX)
       do L = LU,1,-1
         write(6,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
       enddo

      endif

   99 continue

      END SUBROUTINE PHOTO_JX


!<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
      subroutine OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
              FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in) ::   DTAUX(JXL1_,W_),POMEGAX(8,JXL1_,W_)
      real*8, intent(in) ::   AMF2(2*JXL1_+1,2*JXL1_+1)
      real*8, intent(in) ::   U0,RFL(W_)
      integer, intent(in) ::  JXTRA(JXL2_+1), LU
      real*8, intent(out) ::FJACT(JXL_,W_),FJTOP(W_),FJBOT(W_),FSBOT(W_)
      real*8, intent(out) ::  FJFLX(JXL_,W_),FLXD(JXL1_,W_),FLXD0(W_)

      integer JNDLEV(JXL_),JNELEV(JXL1_)
      integer JADDLV(JXL2_+1),JADDTO(JXL2_+1),L2LEV(JXL2_+1)
      integer JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,ND
      integer L1U,L2U,   LZ0,LZ1,LZMID
      real*8   SUMT,SUMJ

      real*8  DTAU(JXL1_+1,W_),POMEGAJ(M2_,JXL2_+1,W_),TTAU(JXL2_+1,W_)
      real*8  FTAU2(JXL2_+1,W_),POMEGAB(M2_,W_)
      real*8  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,FJFLX0
      real*8, dimension(W_) :: TAUBTM,TAUTOP,FBTM,FTOP,ZFLUX
!--- variables used in mie code-----------------------------------------
      real*8, dimension(W_)         :: FJT,FJB
      real*8, dimension(N_,W_)      :: FJ,FZ,ZTAU
      real*8, dimension(M2_,N_,W_)  :: POMEGA
      real*8, dimension(2*JXL1_,W_)  :: FLXD2

!---there is a parallel correspondence:
!  dimension of JX arrays JXL_ .ge. dimension that CTM is using = L_
!  but calculation is done for L_=LU, L1_=L1U, L2_=L2U lengths of CTM
!
!  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
!
! in:
!     DTAUX(1:L1_,1:W_) = optical depth of each layer
!     POMEGAX(1:8,1:L1_,1:W_) = scattering phase fn (multiplied by s-s abledo)
!     U0  = cos (SZA)
!     RFL(1:W_) = Lambertian albedo of surface
!     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
!        AMF2 now does both edges and middle of CTM layers
!     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
! out:
!     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
!  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
!     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
!     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)
!     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)
!     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L
!        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
!     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
!        this should take into account sphericity, and is not just = mu0
!     FLXD0(1:W_) = sum of solar flux deposited in atmos
!        does NOT include flux on lower surface, does NOT mean absorbed!
!-----------------------------------------------------------------------
!
!     DTAU     Local optical depth of each CTM level
!     TTAU     Optical depth of air vertically above each point (to top of atm)
!     FTAU2     Attenuation of solar beam
!     POMEGAJ  Scattering phase function
!
!---new ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the
!   factor increase from sub-layer to sub-layer
!
!---------------------SET UP FOR MIE CODE-------------------------------
!
!-----------------wavelength independent--------------------------------
!
!  Transpose the ascending TTAU grid to a descending ZTAU grid.
!  Double the resolution - TTAU points become the odd points on the
!  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
!  Odd point added at top of grid for unattenuated beam   (Z='inf')
!
!  The following mapping holds for JADDLV=0
!        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
!        Top:       TTAU(L2_)  ==> ZTAU(3)
!        Infinity:     0.0     ==> ZTAU(1)
!        index: 2*(L2_+1-L2)+1 ==> LZ
!
!  Mie scattering code only used from surface to level L2_
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
!  Insert new levels, working downwards from the top of the atmosphere
!  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
!  to be incremented linearly, and the flux fz to be attenuated top-down
!    (avoiding problems where lower level fluxes are zero).
!------------------------------------------------------------------------
!
!  Ascend through atmosphere transposing grid and adding extra points
!  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
!  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
!    because we need to insert the intermediate layers (even LZ) for the
!    asymmetric scattering code.
!
!  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
!    order, expanded, doubled-level scatter grid.
!    Note that we need to deal with the expansion by JADD levels (L2L).
!      These JADDLV levels are skipped and need to be interpolated later.
!    Note that only odd LZ levels are filled,
!
!----------------------re-grid data---------------------------------------------
!  Calculate cumulative total and define levels we want J-values at.
!  Sum upwards for levels, and then downwards for Mie code readjustments.
!
!     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
!           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
!     JADDLV(L2)  Number of new levels actually added at each wavelength
!            where JADDLV = 0 when there is effectively no FTAU2
!     JADDTO(L2)   Total number of new levels to add to and above level (L2)
!     JNDLEV(L) = L2 index that maps on CTM mid-layer L
!
!---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
!---    JADDLV is taken from JXTRA, which is based on visible OD.
!---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
!---these should be fixed for all wavelengths to lock-in the array sizes

      if (LU .gt. JXL_) then
        call EXITC (' OPMIE:  JXL_ .lt. L_')
      endif

      L1U = LU + 1
      L2U = 2*LU + 2


      do L2 = 1,L2U,1
        JADDLV(L2) = JXTRA(L2)
      enddo
        JADDTO(L2U+1) = 0
      do L2 = L2U,1,-1
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
      enddo

!---expanded grid now included CTM edge and mid layers plus expanded
!---    grid to allow for finer delta-tau at tops of clouds.
!---    DIM of new grid = L2U + JADDTO(1) + 1

!---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
!     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2U+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
      enddo

!---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
!---JNELEV(L=1:L_) = L2-index for top of layer L
      do L = 1,LU
        JNDLEV(L) = L2LEV(2*L)
        JNELEV(L) = L2LEV(2*L+1)
      enddo
        JNELEV(LU+1) = 0  !need to set this to top-of-atmosphere

      ND = 2*L2U + 2*JADDTO(1) + 1

      if(ND .gt. N_) then
        call EXITC (' overflow of scatter arrays: ND > N_')
      endif

!----------------begin wavelength dependent set up------------------------------

!---Reinitialize arrays
      ZTAU(:,:)     = 0.d0
      FZ(:,:)       = 0.d0
      POMEGA(:,:,:) = 0.d0

      FJACT(:,:) = 0.d0
      FJTOP(:) = 0.d0
      FJBOT(:) = 0.d0
      FSBOT(:) = 0.d0
      FJFLX(:,:) = 0.d0
      FLXD(:,:) = 0.d0
      FLXD0(:) = 0.d0
      FJT(:) = 0.d0
      FJB(:) = 0.d0
      FJ(:,:) = 0.d0
      FZ(:,:) = 0.d0
      ZTAU(:,:) = 0.d0

      do K=1,W_
      if (FL(K) .gt. 1.d0) then

!---Set up optical depth DTAU(L)
       do L = 1,L1U
        DTAU(L,K) = DTAUX(L,K)
       enddo
        DTAU(L1U+1,K) = 0.d0

!---Define the total scattering phase fn for each CTM layer L=1:L_+1
!---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
!---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
       do L = 1,L1U
        do I = 1,M2_
          POMEGAJ(I,L,K) = POMEGAX(I,L,K)
        enddo
       enddo

!---Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
!---      at the middle & edges of the CTM layers L=1:2*L1_+1
!---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
!---  note that DTAU(L1_) is optical depth in the FULL CTM layer just above
        FTAU2(:,:) = 0.d0
        FTAU2(L2U+1,:) = 1.0d0
       do LL = 1,2*L1U+1
         L = (LL+1)/2
        if (AMF2(LL,LL) .gt. 0.0d0) then
           XLTAU = 0.0d0
         do II = 1,2*L1U+1
           I = (II+1)/2
           XLTAU = XLTAU + 0.5d0*DTAU(I,K)*AMF2(II,LL)
         enddo
         if (XLTAU .lt. 76.d0) then   ! zero out flux at 1e-33
          FTAU2(LL,K) = exp(-XLTAU)
         endif
        endif
       enddo

!---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
!---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
          FLXD2(:,:) = 0.d0
       do LL = 1,2*L1U
        if (AMF2(LL,LL) .gt. 0.d0) then
          FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL)
        endif
       enddo
        if (AMF2(1,1) .gt. 0.d0) then
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1)
        else
          FSBOT(K) = 0.d0
        endif

       do LL = 2,2*L1U,2
         L=LL/2
         FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K)
       enddo

!---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
!---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
        FLXD0(K) = 0.d0
       if (AMF2(2*L1U,2*L1U) .gt. 0.d0) then
        do L=1,L1U
         FLXD0(K) = FLXD0(K) + FLXD(L,K)
        enddo
       endif

!------------------------------------------------------------------------
!  Take optical properties on CTM layers and convert to a photolysis
!  level grid corresponding to layer centres and boundaries. This is
!  required so that J-values can be calculated for the centre of CTM
!  layers; the index of these layers is kept in the JNDLEV array.
!------------------------------------------------------------------------
!---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
!---    points (1:L_) plus 1 for the mid point of added top layer.
!---combine these edge- and mid-layer points into grid of size:
!---              L2_+1 = 2*L1_+1 = 2*L_+3
!---calculate column optical depths above each level, TTAU(1:L2_+1)
!---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD

        TTAU(L2U+1,K) = 0.0d0
       do L2 = L2U,1,-1
        L          = (L2+1)/2
        DTAUJ      = 0.5d0 * DTAU(L,K)
        TTAU(L2,K)   = TTAU(L2+1,K) + DTAUJ
       enddo

!----solar flux incident on lower boundary & Lambertian reflect factor:
       if (FSBOT(K) .gt. 0.d0) then
        ZFLUX(K) = FSBOT(K)*RFL(K)/(1.d0+RFL(K))
       else
        ZFLUX(K) = 0.d0
       endif

!  Calculate scattering properties, level centres then level boundaries
!>>>>>be careful of order, we are overwriting/shifting the 'POMEGAJ' upward
!     in index
       do L2 = L2U,2,-2
        L   = L2/2
        do I = 1,M2_
          POMEGAJ(I,L2,K) = POMEGAJ(I,L,K)
        enddo
       enddo
!---lower boundary value is set (POMEGAJ(I,1)), but set upper:
       do I = 1,M2_
         POMEGAJ(I,L2U+1,K) = POMEGAJ(I,L2U,K)
       enddo
!---now have POMEGAJ filled at even points from L2=3:L2_-1
!---use inverse interpolation for correct tau-weighted values at edges
       do L2 = 3,L2U-1,2
        TAUDN = TTAU(L2-1,K)-TTAU(L2,K)
        TAUUP = TTAU(L2,K)-TTAU(L2+1,K)
        do I = 1,M2_
          POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN + &
                 POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
        enddo
       enddo

!---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
!---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface

       do L2 = 1,L2U+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = ND + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K)
          FZ(LZ,K)   = FTAU2(L2,K)
        do I=1,M2_
          POMEGA(I,LZ,K) = POMEGAJ(I,L2,K)
        enddo
       enddo

!   Now go thru the pairs of L2 levels to see if we need JADD levels
       do L2 = 1,L2U             ! L2 = index of CTM edge- and mid-layers
         L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)
         LZ  = ND + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
         L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels

        if (L22 .gt. 0) then
          TAUBTM(K) = TTAU(L2,K)
          TAUTOP(K) = TTAU(L2+1,K)
          FBTM(K)   = FTAU2(L2,K)
          FTOP(K)   = FTAU2(L2+1,K)
         do I = 1,M2_
          POMEGAB(I,K) = POMEGAJ(I,L2,K)
         enddo

!---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
!---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM
         ATAUZ = exp(-log(TAUBTM(K)/max(TAUTOP(K),ATAU0))/float(L22+1))
         do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt arrays
          ZTAU(LZZ,K) = TAUBTM(K) * ATAUZ

!---fraction from TAUBTM=>TAUTOP
          ATAUA=(TAUBTM(K)-ZTAU(LZZ,K))/(TAUBTM(K)-TAUTOP(K))
!---solar flux at interp-levels: use exp(TAU/U0) if U0>0.02 (89 deg),
!---else scale by TAU
          if (U0 .gt. 0.02d0) then
            FZ(LZZ,K) = FTOP(K) * exp((TAUTOP(K)-ZTAU(LZZ,K))/U0)
          else
            if (FTOP(K).lt.1.d-30 .or. FBTM(K).lt.1.d-30) then
              FZ(LZZ,K) = 0.d0
            else
              FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA
            endif
          endif
          do I = 1,M2_
            POMEGA(I,LZZ,K) = POMEGAB(I,K) + &
                     ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))
          enddo
            TAUBTM(K)    = ZTAU(LZZ,K)
            FBTM(K)      = FZ(LZZ,K)
          do I = 1,M2_
            POMEGAB(I,K) = POMEGA(I,LZZ,K)
          enddo
         enddo
        endif
       enddo

!   Now fill in the even points with simple interpolation in scatter arrays:
       do LZ = 2,ND-1,2
         ZTAU(LZ,K) = 0.5d0*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
         FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
        do I=1,M2_
         POMEGA(I,LZ,K) = 0.5d0*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
        enddo
       enddo

      endif
      enddo  ! wavelength loop!

!-----------------------------------------------------------------------
       call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!-----------------------------------------------------------------------

!---Move mean intensity from scatter array FJ(LZ=1:ND)
!---              to CTM mid-level array FJACT(L=1:L_)

      do K=1,W_
      if (FL(K) .gt. 1.d0) then

!---mean intensity at mid-layer:  4*<I> + solar
!       do L = 1,LU
!        L2L = JNDLEV(L)
!        LZ  = ND+2 - 2*L2L
!        FJACT(L,K) = 4.d0*FJ(LZ,K) + FZ(LZ,K)
!       enddo

!---mean intensity averaged throughout layer:
       do L = 1,LU
         LZ0 = ND+2 - 2*JNELEV(L)
        if (L .gt. 1) then
         LZ1 = ND+2 - 2*JNELEV(L-1)
        else
         LZ1 = ND
        endif
         SUMJ = (4.d0*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+2,K)-ZTAU(LZ0,K)) &
               +(4.d0*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-2,K))
         SUMT = ZTAU(LZ0+2,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-2,K)

        do LZ = LZ0+2,LZ1-2,2
         SUMJ =SUMJ+(4.d0*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+2,K)-ZTAU(LZ-2,K))
         SUMT =SUMT + ZTAU(LZ+2,K)-ZTAU(LZ-2,K)
        enddo
        FJACT(L,K) = SUMJ/SUMT
       enddo

!---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
!---      average (tau-wtd) the h's just above and below the L-edge
       do L = 1,LU
        L2L = JNELEV(L)
        LZ  = ND+2 - 2*L2L
        FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
        FJFLX(L,K)=4.d0*(FJ(LZ-1,K)*FJFLX0 + FJ(LZ+1,K)*(1.d0-FJFLX0))
       enddo

!---diffuse fluxes reflected at top, incident at bottom
         FJTOP(K) = FJT(K)
         FJBOT(K) = FJB(K)

      endif
      enddo  ! wavelength loop!

      END SUBROUTINE OPMIE


!-----------------------------------------------------------------------
      subroutine MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_) &
                             ,RFL(W_),U0,ZFLUX(W_)
      real*8, intent(out) ::  FJ(N_,W_),FJT(W_),FJB(W_)

      real*8  PM(M_,M2_),PM0(M2_)
      integer I, IM  ,K

!-----------------------------------------------------------------------
!   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Solution of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
!-----------------------------------------------------------------------
      do I = 1,M_
       call LEGND0 (EMU(I),PM0,M2_)
       do IM = 1,M2_
         PM(I,IM) = PM0(IM)
       enddo
      enddo

       call LEGND0 (-U0,PM0,M2_)
       do IM=1,M2_
         PM0(IM) = 0.25d0*PM0(IM)
       enddo

!---BLKSLV now called with all the wavelength arrays (K=1:W_)

      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)

      END SUBROUTINE MIESCT


!-----------------------------------------------------------------------
      subroutine LEGND0 (X,PL,N)
!-----------------------------------------------------------------------
!---Calculates ORDINARY Legendre fns of X (real)
!---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      implicit none
      integer, intent(in) :: N
      real*8, intent(in)  :: X
      real*8, intent(out) :: PL(N)
      integer I
      real*8  DEN
!---Always does PL(2) = P[1]
        PL(1) = 1.d0
        PL(2) = X
        do I = 3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN)
        enddo

      END SUBROUTINE LEGND0


!-----------------------------------------------------------------------
      subroutine BLKSLV &
           (FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
!-----------------------------------------------------------------------
!  Sets up and solves the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!  This goes back to the old, dumb, fast version 5.3
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_) &
                             ,PM(M_,M2_),PM0(M2_) &
                             ,RFL(W_),ZFLUX(W_)
      real*8, intent(out) ::  FJ(N_,W_),FJTOP(W_),FJBOT(W_)

      real*8, dimension(M_,N_,W_)    ::  A,C,H,   RR

      real*8, dimension(M_,M_,N_,W_) ::  B,AA,CC,  DD
      real*8, dimension(M_,M_) ::  E
      real*8  SUMB,SUMBX,SUMT
      integer I, J, K, L

      do K = 1,W_
      if (FL(K) .gt. 1.0d0) then
       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K), &
           PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K), &
                   A(1,1,K),H(1,1,K),C(1,1,K), ND)
      endif
      enddo

      do K = 1,W_
      if (FL(K) .gt. 1.0d0) then
!-----------UPPER BOUNDARY L=1
       L = 1
        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,1,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,1,K) = -E(I,1)*CC(1,J,1,K)-E(I,2)*CC(2,J,1,K) &
                        -E(I,3)*CC(3,J,1,K)-E(I,4)*CC(4,J,1,K)
         enddo
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K) &
                    +E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
        enddo

!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,L,K) = - E(I,J)*C(J,L,K)
         enddo
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

       enddo

!---------FINAL DEPTH POINT: L=ND
       L = ND
        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) &
           + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K) &
           + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) &
           - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K) &
           - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

!---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
!---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
!---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
!---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    +E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

!-----------BACK SOLUTION
       do L = ND-1,1,-1
        do J = 1,M_
         RR(J,L,K) = RR(J,L,K) &
          + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K) &
          + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
        enddo
       enddo

!----------mean J & H
       do L = 1,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2) &
                + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       enddo
       do L = 2,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2) &
                + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       enddo

!---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
!---FJBOT = scaled diffuse flux onto surface:
!---ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
!---SUMBX = flux from Lambert reflected I+
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2) &
            + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2) &
            + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)
       SUMBX = 4.d0*SUMB*RFL(K)/(1.0d0 + RFL(K)) + ZFLUX(K)

       FJTOP(K) = 4.d0*SUMT
       FJBOT(K) = 4.d0*SUMB - SUMBX

      endif
      enddo

      END SUBROUTINE BLKSLV


!-----------------------------------------------------------------------
      subroutine GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0 &
                    ,B,CC,AA,A,H,C,  ND)
!-----------------------------------------------------------------------
!  Generates coefficient matrices for the block tri-diagonal system:
!               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  ND
      real*8, intent(in)  ::  POMEGA(M2_,N_),PM(M_,M2_),PM0(M2_)
      real*8, intent(in)  ::  ZFLUX,RFL
      real*8, intent(in),dimension(N_) :: FZ,ZTAU

      real*8, intent(out),dimension(M_,M_,N_) ::  B,AA,CC
      real*8, intent(out),dimension(M_,N_) ::  A,C,H

      integer I, J, K, L1,L2,LL
      real*8  SUM0, SUM1, SUM2, SUM3
      real*8  DELTAU, D1, D2, SURFAC
!
      real*8, dimension(M_,M_) :: S,T,U,V,W
!---------------------------------------------

!---------upper boundary:  2nd-order terms
       L1 = 1
       L2 = 2
       do I = 1,M_
        SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
       + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
       + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
       + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
       + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
       + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
       + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
       + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
       + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo
!-------------upper boundary, 2nd-order, C-matrix is full (CC)
         DELTAU = ZTAU(L2) - ZTAU(L1)
         D2 = 0.25d0*DELTAU
       do I = 1,M_
        do J = 1,M_
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
         CC(I,J,L1) = D2*U(I,J)
        enddo
         H(I,L1) = H(I,L1) + 2.0d0*D2*C(I,L1)
         A(I,L1) = 0.0d0
       enddo
       do I = 1,M_
        D1 = EMU(I)/DELTAU
        B(I,I,L1)  = B(I,I,L1) + D1
        CC(I,I,L1) = CC(I,I,L1) - D1
       enddo

!------------intermediate points:  can be even or odd, A & C diagonal
!---mid-layer h-points, Legendre terms 2,4,6,8
       do LL=2,ND-1,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4) &
         + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4) &
          +POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

!---odd-layer j-points, Legendre terms 1,3,5,7
       do LL=3,ND-2,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3) &
         + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3) &
          +POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

!---------lower boundary:  2nd-order terms
       L1 = ND
       L2 = ND-1
       do I = 1,M_
        SUM0 = &
         POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
       + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = &
         POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
       + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = &
         POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
       + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = &
         POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
       + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = &
          POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
        + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = &
          POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
        + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = &
          POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
        + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = &
          POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
        + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo

!------------lower boundary, 2nd-order, A-matrix is full (AA)
         DELTAU = ZTAU(L1) - ZTAU(L2)
         D2 = 0.25d0*DELTAU
         SURFAC = 4.0d0*RFL/(1.0d0 + RFL)
       do I = 1,M_
          D1 = EMU(I)/DELTAU
          SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4))
          SUM1 = SURFAC*SUM0
        do J = 1,M_
         AA(I,J,L1) = - D2*U(I,J)
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
        enddo
         H(I,L1) = H(I,L1) - 2.0d0*D2*C(I,L1) + SUM0*ZFLUX
       enddo

       do I = 1,M_
          D1 = EMU(I)/DELTAU
        AA(I,I,L1) = AA(I,I,L1) + D1
        B(I,I,L1)  = B(I,I,L1) + D1
        C(I,L1) = 0.0d0
       enddo

      END SUBROUTINE GEN_ID

!<<<<<<<<<<<<<<<<<<<<<<<<<<end fastJX core code<<<<<<<<<<<<<<<<<<<<<<<<<




!<<<<<begin fastJX subroutines called from PHOTO_JX or OPMIE<<<<<<<<<<<<

!------------------------------------------------------------------------------
      subroutine OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)
!------------------------------------------------------------------------------
!---set CLOUD fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!      ----FJX ver 7.0+ clouds separated
! 01 W_C02 (C1/Deir)GAMMA:r-m=2.0/alf=6 n=1.335   reff=3.000___G=19.55_rho=1.000
! 02 W_C04 (C1/Deir)GAMMA:r-m=4.0/alf=6 n=1.335   reff=6.000___G=78.19_rho=1.000
! 03 W_C08 (C1/Deir)GAMMA:r-m=8.0/alf=2 n=1.335   reff=12.00___G=301.1_rho=1.000
! 04 W_C13 (C1/Deir)GAMMA:r-m=13./alf=2 n=1.335   reff=20.00___G=472.9_rho=1.000
! 05 W_L06 (W/Lacis)GAMMA:r-m=5.5/alf=11/3        reff=10.00___G=183.9_rho=1.000
! 06 Ice-Hexagonal (Mishchencko)                  reff=50.00___G=999.9_rho=0.917
! 07 Ice-Irregular (Mishchencko)                  reff=50.00___G=999.9_rho=0.917
! 08 S-Bkg LOGN:r=.090 s=.600 n=1.514/.../1.435   reff=0.221___G=.0523_rho=1.630
! 09 S-Vol LOGN:r=.080 s=.800 n=1.514/.../1.435   reff=0.386___G=.0721_rho=1.630
!
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     ODCLD      ! optical depth of cloud layer @600 nm
      integer,intent(inout)::  NDCLD      ! index of cloud layer:  4:13

      integer I,J
      real*8  XTINCT, REFF,RHO

!---later versions should allow for interpolation, averaging of R-eff

!---default cloud type C1, Reff = 12 microns
      if (NDCLD .lt. 1 .or. NDCLD .gt. 9) then
         NDCLD = 3
      endif

!--rescale OD by Qext at 600 nm (J=4)
      do J=1,5
         OPTD(J) = ODCLD * QCC(J,NDCLD)/QCC(4,NDCLD)
         SSALB(J) = SCC(J,NDCLD)
        do I=1,8
         SLEG(I,J) =  PCC(I,J,NDCLD)
        enddo
      enddo

      END SUBROUTINE OPTICL


!------------------------------------------------------------------------------
      subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,K)
!------------------------------------------------------------------------------
!---for the UCI aerosol data sets, calculates optical properties at fast-JX's
!              std 5 wavelengths:200-300-400-600-999nm
!
!---UCI aersols optical data  v-7.0+
!01 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435  reff=0.221___G=.0523_rho=1.630
!02 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435  reff=0.386___G=.0721_rho=1.630
!03 UT-sulfate LOGN:r=0.05 s=.693 n=1.44          reff=0.166___G=.0205_rho=1.769
!04 UT-sulfate LOGN:r=0.05 s=.693 n=1.46          reff=0.166___G=.0205_rho=1.769
!05 UT-sulfatM LOGN:r=.050 s=.642 n=1.53          reff=0.140___G=.0179_rho=1.769
!06 UM-BC1     LOGN:r=.050 s=.642 n=1.80+0.50i    reff=0.140___G=.0179_rho=1.500
!07 UM-BC2     LOGN:r=.080 s=.501 n=1.80+0.50i    reff=0.150___G=.0332_rho=1.500
!08 UM-BB08 (%BC)LOGN:r=.080 s=.500 n=1.552+0.04i reff=0.149___G=.0331_rho=1.230
!09 UM-FF04 (%BC)LOGN:r=.050 s=.642 n=1.541+0.02i reff=0.140___G=.0179_rho=1.212
!10 UM-FF10 (%BC)LOGN:r=.050 s=.642 n=1.557+0.05i reff=0.140___G=.0179_rho=1.230
!11 MDust.15  (R.V. Martin generated phase fns)   reff=0.150___G=1.000_rho=2.600
!12 MDust.25  (R.V. Martin generated phase fns)   reff=0.250___G=1.000_rho=2.600
!13 MDust.40  (R.V. Martin generated phase fns)   reff=0.400___G=1.000_rho=2.600
!14 MDust.80  (R.V. Martin generated phase fns)  reff=0.800___G=1.000_rho=2.6001
!15 MDust1.5  (R.V. Martin generated phase fns)   reff=1.500___G=1.000_rho=2.600
!16 MDust2.5  (R.V. Martin generated phase fns)   reff=2.500___G=1.000_rho=2.600
!17 MDust4.0  (R.V. Martin generated phase fns)   reff=4.000___G=1.000_rho=2.600
!
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     RELH       ! relative humidity (0.00->1.00+)
      integer,intent(inout)::     K       ! index of cloud/aerosols

      integer I,J
      real*8  XTINCT, REFF,RHO

      if (K.gt.NAA .or. K.lt.1)  &
        call EXITC ('OPTICA: aerosol index out-of-range')

         REFF = RAA(K)
         RHO = DAA(K)
      do J=1,5
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
         XTINCT = 0.75d0*QAA(J,K)/(REFF*RHO)
         OPTD(J) = PATH*XTINCT
         SSALB(J) = SAA(J,K)
       do I=1,8
         SLEG(I,J) =  PAA(I,J,K)
       enddo
      enddo

      END SUBROUTINE OPTICA


!------------------------------------------------------------------------------
      subroutine OPTICM (OPTD,SSALB,SLEG, PATH,RELH,LL)
!------------------------------------------------------------------------------
!---for the U Michigan aerosol data sets, this generate fast-JX data formats.
!---NB Approximates the Legendre expansion(L) of the scattering phase fn
!---           as (2*L+1)*g**L
!---UMAER(I,J,K,L):
!   I=1:3 = [SSAbldeo, g, k-ext(m2/g)]
!   J=1:5 = [200, 300, 400, (550,) 600 , 1000 nm]
!   K=1:21= [0, 5, 10, 15, ..., 90, 95, 99 %RelHum]
!   L=1:33= UM aerosol types [SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, FF00(0%BC),
!                      FF02, ...FF14(14%BC), BB00, BB02, ...BB30(30%BC)]
      implicit none

      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real*8, intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     RELH       ! relative humidity (0.00->1.00)
      integer,intent(in)::     LL         ! index of cloud/aerosols

      integer KR,J,L
      real*8  R,FRH, GCOS, XTINCT

!---calculate fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!---extrapolate phase fn from first term (g)
      L = LL

      if (L.lt.1 .or. L.gt.33)  &
          call EXITC ('OPTICM: aerosol index out-of-range')

!---pick nearest Relative Humidity
      KR =  20.d0*RELH  + 1.5d0
      KR = max(1, min(21, KR))

      do J=1,5
       SSALB(J) = UMAER(1,J,KR,L)
         XTINCT = UMAER(3,J,KR,L)
       OPTD(J) = PATH*XTINCT
         GCOS   = UMAER(2,J,KR,L)
       SLEG(1,J) =  1.d0
       SLEG(2,J) =  3.d0*GCOS
       SLEG(3,J) =  5.d0*GCOS**2
       SLEG(4,J) =  7.d0*GCOS**3
       SLEG(5,J) =  9.d0*GCOS**4
       SLEG(6,J) = 11.d0*GCOS**5
       SLEG(7,J) = 13.d0*GCOS**6
       SLEG(8,J) = 15.d0*GCOS**7
      enddo

      END SUBROUTINE OPTICM


!-----------------------------------------------------------------------
      subroutine JRATET(PPJ,TTJ,FFF, VALJL,LU,NJXU)
!-----------------------------------------------------------------------
! in:
!        PPJ(L_+1) = pressure profile at edges
!        TTJ(L_+1) = = temperatures at mid-level
!        FFF(K=1:NW, L=1:L_) = mean actinic flux
! out:
!        VALJL(L_,JX_)  JX_ = no of dimensioned J-values in CTM code
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: LU,NJXU
      real*8, intent(in)  ::  PPJ(LU+1),TTJ(LU+1)
      real*8, intent(inout)  ::  FFF(W_,LU)
      real*8, intent(out), dimension(LU,NJXU) ::  VALJL

      real*8  VALJ(X_)
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV

      if (NJXU .lt. NJX) then
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(A,2I5)')  'NJXU<NJX',NJXU,NJX
        call EXITC(' JRATET:  CTM has not enough J-values dimensioned')
      endif
      do L = 1,LU
!need temperature, pressure, and density at mid-layer (for some quantum yields):
          TT   = TTJ(L)
         if (L .eq. 1) then
          PP = PPJ(1)
         else
          PP  = (PPJ(L)+PPJ(L+1))*0.5d0
         endif
          DD = 7.24e18*PP/TT

!---if W_=18/12, must zero bin-11/5 below 100 hPa:  216-222 & 287-291 nm mix
!   since O2 e-fold is too weak (by averaging the 287-291)
!   this combination of wavelengths works fine in stratosphere with O3 attenuation
        if (PP .gt. 100.d0) then
          if (W_ .eq. 18) then
            FFF(11,L) = 0.d0
          elseif (W_ .eq. 12) then
            FFF(5,L) = 0.d0
          endif
        endif

        do J = 1,NJX
          VALJ(J) = 0.d0
        enddo

        do K = 1,W_
          call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1), &
                         TQQ(2,1),QO2(K,2), TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1), &
                         TQQ(2,2),QO3(K,2), TQQ(3,2),QO3(K,3), LQQ(2))
          call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1), &
                         TQQ(2,3),Q1D(K,2), TQQ(3,3),Q1D(K,3), LQQ(3))
            QO31D  = QO31DY*QO3TOT
           VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
           VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
           VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
        enddo

        do J = 4,NJX
         do K = 1,W_
!---also need to allow for Pressure interpolation if SQQ(J) = 'p'
          if (SQQ(J) .eq.'p') then
           call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J), &
                     TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          else
           call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J), &
                     TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          endif
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
         enddo
        enddo

        do J=1,NJX
          VALJL(L,J) = VALJ(J)
        enddo

      enddo

      END SUBROUTINE JRATET


!-----------------------------------------------------------------------
      subroutine X_interp (TINT,XINT, T1,X1, T2,X2, T3,X3, L123)
!-----------------------------------------------------------------------
!  up-to-three-point linear interpolation function for X-sections
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in)::  TINT,T1,T2,T3, X1,X2,X3
      integer,intent(in)::  L123
      real*8,intent(out)::  XINT

      real*8  TFACT

      if (L123 .le. 1) then
           XINT = X1
      elseif (L123 .eq. 2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
      else
        if (TINT.le. T2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
        else
             TFACT = max(0.d0,min(1.d0,(TINT-T2)/(T3-T2) ))
           XINT = X2 + TFACT*(X3 - X2)
        endif
      endif

      END SUBROUTINE X_interp


!-----------------------------------------------------------------------
      subroutine JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU6,POMEG6,JXTRA,LU)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!---the CTM has L_ = LU layers and fast-JX adds layer LU+1
!---the pressure and altitude(Z) are on layer edge (LU+2)

      integer,intent(in)                  :: LU
      real*8, intent(in), dimension(LU+2) :: PPJ,ZZJ
      real*8, intent(in), dimension(LU+1) :: TTJ,DDJ,OOJ,DTAU6
      real*8, intent(in), dimension(8,LU+1) :: POMEG6
      integer,intent(in), dimension(LU+LU+3) :: JXTRA
!-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP,DAIR,DOZO

      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(4a)') '   L z(km)     p      T   ', &
       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
       '  g(cos) CTM lyr=>'

      L = LU+2
      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZZJ(L)*1.d-5,PPJ(L)

          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)

        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5
          DAIR = DDJ(L)/DELZ
          DOZO = OOJ(L)/DELZ

          if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZKM,PPJ(L),TTJ(L),DAIR,DOZO,COLO2,COLO3,DTAU6(L), &
            POMEG6(1,L),POMEG6(2,L)/3.d0, JXTRA(L+L),JXTRA(L+L-1)

        enddo

      END SUBROUTINE JP_ATM


!-----------------------------------------------------------------------
      subroutine JP_ATM0(PPJ,TTJ,DDJ,OOJ,ZZJ, LU)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!---the CTM has L_ = LU layers and fast-JX adds layer LU+1
!---the pressure and altitude(Z) are on layer edge (LU+2)

      integer,intent(in)                  :: LU
      real*8, intent(in), dimension(LU+2) :: PPJ,ZZJ
      real*8, intent(in), dimension(LU+1) :: TTJ,DDJ,OOJ
!-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP
      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(4a)') '   L z(km)     p      T   ', &
       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
       '  g(cos) CTM lyr=>'
      L = LU+2
      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZZJ(L)*1.d-5,PPJ(L)
          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)
        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5
      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
            L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ, &
            COLO2,COLO3
        enddo

      END SUBROUTINE JP_ATM0


!-----------------------------------------------------------------------
      subroutine SPHERE2(U0,RAD,ZHL,ZZHT,AMF2, L1U,LJX1U)
!-----------------------------------------------------------------------
!----new v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
!  This new AMF2 does each of the half-layers of the CTM separately,
!     whereas the original, based on the pratmo code did the whole layers
!     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
!  Since fast-JX is meant to calculate the intensity at the mid-layer, the
!     solar beam at low sun (interpolated between layer edges) was incorrect.
!  This new model does make some approximations of the geometry of the layers:
!     the CTM layer is split evenly in mass (good) and in height (approx).
!
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!-----------------------------------------------------------------------
! in:
!     U0      cos(solar zenith angle)
!     RAD     radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottome edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
! out:
!     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
!          these are calcualted for both layer middle and layer edge
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) ::   L1U, LJX1U
      real*8, intent(in)  ::   U0,RAD,ZHL(L1U+1),ZZHT
      real*8, intent(out) ::   AMF2(2*LJX1U+1,2*LJX1U+1)

      integer, parameter  ::  LSPH_ = 100

!     RZ      Distance from centre of Earth to each point (cm)
!     RQ      Square of radius ratios
!     SHADHT  Shadow height for the current SZA
!     XL      Slant path between points

      integer  I, J, K, II, L2
      real*8   XMU1,XMU2,XL,DIFF,SHADHT,RZ(LSPH_+1)
      real*8   RZ2(2*LSPH_+1),RQ2(2*LSPH_+1)

!--- must have top-of-atmos (NOT top-of-CTM) defined
!      ZHL(L1U+1) = ZHL(L1U) + ZZHT

      if (L1U .gt. LSPH_) then
        call EXITC(' SPHERE2: temp arrays not large enough')
      endif

        RZ(1) = RAD + ZHL(1)
      do II = 2,L1U+1
        RZ(II)   = RAD + ZHL(II)
      enddo

!---calculate heights for edges of split CTM-layers
      L2 = 2*L1U
      do II = 2,L2,2
        I = II/2
        RZ2(II-1) = RZ(I)
        RZ2(II) = 0.5d0*(RZ(I)+RZ(I+1))
      enddo
        RZ2(L2+1) = RZ(L1U+1)
      do II = 1,L2
        RQ2(II) = (RZ2(II)/RZ2(II+1))**2
      enddo

!---shadow height for SZA > 90
      if (U0 .lt. 0.0d0)  then
        SHADHT = RZ2(1)/dsqrt(1.0d0 - U0**2)
      else
        SHADHT = 0.d0
      endif

!---up from the surface calculating the slant paths between each level
!---  and the level above, and deriving the appropriate Air Mass Factor
         AMF2(:,:) = 0.d0

      do 16 J = 1,2*L1U+1

!  Air Mass Factors all zero if below the tangent height
        if (RZ2(J) .lt. SHADHT) goto 16
!  Ascend from layer J calculating AMF2s
        XMU1 = abs(U0)
        do I = J,2*L1U
          XMU2     = dsqrt(1.0d0 - RQ2(I)*(1.0d0-XMU1**2))
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I))
          XMU1     = XMU2
        enddo
!--fix above top-of-atmos (L=L1U+1), must set DTAU(L1U+1)=0
          AMF2(2*L1U+1,J) = 1.d0
!
!  Twilight case - Emergent Beam, calc air mass factors below layer
        if (U0 .ge. 0.0d0) goto 16

!  Descend from layer J
          XMU1       = abs(U0)
         do II = J-1,1,-1
          DIFF        = RZ2(II+1)*sqrt(1.0d0-XMU1**2)-RZ2(II)
          if (II.eq.1)  DIFF = max(DIFF,0.d0)   ! filter
!  Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0d0)  then
            XMU2      = sqrt(1.0d0 - (1.0d0-XMU1**2)/RQ2(II))
            XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2)
            AMF2(II,J) = 2.d0*XL/(RZ2(II+1)-RZ2(II))
            XMU1      = XMU2
!  Lowest level intersected by emergent beam
          else
            XL        = RZ2(II+1)*XMU1*2.0d0
            AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II))
            goto 16
          endif
         enddo

   16 continue

      END SUBROUTINE SPHERE2


!-----------------------------------------------------------------------
      subroutine EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------
!
!    new version 6.1, add sub-layers (JXTRA) to thick cloud/aerosol layers
!    this version sets up log-spaced sub-layers of increasing thickness ATAU
!
!     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
!        This can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer L
!        Set for log-spacing of tau levels, increasing top-down.
!
!     N.B. the TTAU, etc calculated here are NOT used elsewhere

!---The log-spacing parameters have been tested for convergence and chosen
!---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
!---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
!---  ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
!-----------------------------------------------------------------------
!
      implicit none
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real*8,  intent(in) ::  DTAUX(L1X)      !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0
      integer, intent(out)::  JXTRA(L2X+1)    !number of sub-layers to be added
!
      integer JTOTL,I,L,L2
      real*8  TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1
!
!---Reinitialize arrays
      TTAU(:)  = 0.d0
      JXTRA(:) = 0
!
!---combine these edge- and mid-layer points into grid of size:
!---              L2X+1 = 2*L1X+1 = 2*L_+3
!---calculate column optical depths above each level, TTAU(1:L2X+1)
!---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
!
!---Divide thick layers to achieve better accuracy in the scattering code
!---In the original fast-J, equal sub-layers were chosen, this is wasteful
!---and this new code (ver 5.3) uses log-scale:
!---        Each succesive layer (down) increase thickness by ATAU > 1
!---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
!---        4 sub-layers with ODs = 1 - 2 - 4 - 8
!---The key parameters are:
!---        ATAU = factor increase from one layer to the next
!---        ATAUMN = the smallest OD layer desired
!---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
!---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1.d0
      ATAULN = log(ATAU)
        TTAU(L2X+1)  = 0.0d0
      do L2 = L2X,1,-1
        L         = (L2+1)/2
        DTAUJ     = 0.5d0 * DTAUX(L)
        TTAU(L2)  = TTAU(L2+1) + DTAUJ
!---Now compute the number of log-spaced sub-layers to be added in
!---   the interval TTAU(L2) > TTAU(L2+1)
!---The objective is to have successive TAU-layers increasing by factor ATAU >1
!---the number of sub-layers + 1
        if (TTAU(L2) .lt. ATAU0) then
          JXTRA(L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(L2+1))
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5d0)))
        endif
      enddo

!---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2
      do L2 = L2X,1,-1
        JTOTL  = JTOTL + JXTRA(L2)
        if (JTOTL .gt. NX/2)  then
          write(6,'(A,2I5,F9.2)') 'EXTRAL: N_/L2_/L2-cutoff JXTRA:',  &
                                   NX,L2X,L2
          do L = L2,1,-1
            JXTRA(L) = 0
          enddo
          go to 10
        endif
      enddo
  10  continue

      END SUBROUTINE EXTRAL

!<<<<<<<end fastJX subroutines called from PHOTO_JX or OPMIE<<<<<<<<<<<<

!-----------------------------------------------------------------------
      subroutine EXITC(T_EXIT)
!-----------------------------------------------------------------------
      character(len=*), intent(in) ::  T_EXIT

      write(6,'(a)') T_EXIT
      stop

      END SUBROUTINE EXITC


      END MODULE FJX_SUB_MOD

!<<<<<<<<<<<<<<<<<<fastJX codes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<only access to external variable thru cmn_FJX.f and calls<<<<<<
!<<<<<<<<<<<<<<<<<<<<version 7.3+  (2/2015, mjp)<<<<<<<<<<<<<<<<<<<<<<<<
!   fixes problems back to version 6.8 with failure to collapse Xsects
!   with 1 or 3 data blocks when W_ .ne. 18
!<<<<<<<<<<<<<<<<<<<<version 7.2  (6/2013, mjp) w/ RANSET for clouds<<<<

!
! !MODULE: FJX_INIT
!
! !DESCRIPTION: FJX_INIT contains variables and routines to input fast-JX data
!
!
! !INTERFACE:
!
      MODULE FJX_INIT_MOD
!
! !USES:
!
      USE FJX_CMN_MOD

      use Config_module, only: cloudjx_initf, USES, MasterProc
      use ChemRates_mod, only: CM_schemes_ChemRates

      USE FJX_SUB_MOD, ONLY : EXITC

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC SUBROUTINES:
!
      PUBLIC  :: INIT_FJX

      CONTAINS

!-----------------------------------------------------------------------
      subroutine INIT_FJX (TITLEJXX,NJXU,NJXX)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  ::NJXU
      integer, intent(out) ::NJXX
      character*6, intent(out), dimension(NJXU) :: TITLEJXX

      integer  JXUNIT,J, RANSEED

      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,*) ' fast-JX ver-7.3 standalone CTM code'

      if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
        call EXITC(' INIT_JX: invalid no. wavelengths')
      endif

! Use channel 8 to read fastJX data files:
      JXUNIT  = 8

! Read in fast-J X-sections (spectral data)
      call RD_XXX(JXUNIT,trim(cloudjx_initf)//'/FJX_spec.dat')

! Read in cloud scattering data
      call RD_CLD(JXUNIT,trim(cloudjx_initf)//'/FJX_scat-cld.dat')

! Read in aerosols scattering data
      call RD_MIE(JXUNIT,trim(cloudjx_initf)//'/FJX_scat-aer.dat')

! Read in UMich aerosol scattering data
      call RD_UM (JXUNIT,trim(cloudjx_initf)//'/FJX_scat-UMa.dat')

! Read in T & O3 climatology used to fill e.g. upper layers or if O3 not calc.
      call RD_PROF(JXUNIT,trim(cloudjx_initf)//'/atmos_std.dat')

        NJXX = NJX
      do J = 1,NJX
        TITLEJXX(J) = TITLEJX(J)
      enddo

! Read in photolysis rates used in chemistry code and mapping onto FJX J's
!---CTM call:  read in J-values names and link to fast-JX names
      call RD_JS_JX(JXUNIT,trim(cloudjx_initf)//'/FJX_j2j.dat', TITLEJXX,NJXX)

!---setup the random number sequence RAN4
         RANSEED = 66
      call RANSET (NRAN_,RAN4,RANSEED)


      END SUBROUTINE INIT_FJX


!-----------------------------------------------------------------------
      subroutine RD_XXX(NUN,NAMFIL)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
!
!>>>>NEW v-7.3  expanded input, full names & notes
!>>>>NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!>>>>NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (cmn_FXJ.f).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!
!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, JJ, K, IW, NQRD, NWWW,   LQ
      real*8   T_FL

      character*6  TIT_J1S,TIT_J2S,TITLEJ3
      character*16 TIT_J1L
      character*120 TIT_J1N, TIT_SPEC
      character*1  T_XP

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data------------------
!         note that X_ = max # Xsects read in
!                   NJX = # fast-JX J-values derived from this (.le. X_)
! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects
! >>>> W_ = 8  <<<< extreme trop-only, discard WL #1-4 and #9-10, some X-sects
      if (W_.ne.18 .and. W_.ne.12 .and. W_.ne.8) then
       call EXITC(' no. wavelengths wrong: W_ .ne. 8,12,18')
      endif

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(a120)',err=4) TIT_SPEC
      read (NUN,*,err=4)
      read (NUN,'(i5)',err=4) NWWW
!print
         if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a120)') TIT_SPEC
         if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(i8)') NWWW
       NW1 = 1
       NW2 = NWWW
!----w-params:  1=w-bins, 2=solar(photons), 3=solar(W/m2), 4=Y-PAR spectrum, 5=Rayleigh Xsects
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
!print
         if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N
!not liked by LUMI compiler      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
!          T_XP,T_FL, (WL(IW),IW=1,NWWW)
         read (NUN,*) (WL(IW),IW=1,6)
         read (NUN,*) (WL(IW),IW=7,12)
         read (NUN,*) (WL(IW),IW=13,NWWW)

         
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
!print
         if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          T_XP,T_FL, (FL(IW),IW=1,NWWW)
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
!print
         if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          T_XP,T_FL, (FW(IW),IW=1,NWWW)
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
!print
         if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          T_XP,T_FL, (FP(IW),IW=1,NWWW)
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
!print
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          T_XP,T_FL, (QRAYL(IW),IW=1,NWWW)

!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
!---NB the O3 and q-O3-O1D are at different temperatures and cannot be combined

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,1), (QO2(IW,3),IW=1,NWWW)
        TITLEJX(1) = TIT_J1S
        TITLEJL(1) = TIT_J1L
        LQQ(1) = 3
!print
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,2), (QO3(IW,3),IW=1,NWWW)
        TITLEJX(2) = TIT_J1S
        TITLEJL(2) = TIT_J1L
        LQQ(2) = 3
!print
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)
        TITLEJX(3) = TIT_J1S
        TITLEJL(3) = TIT_J1L
        LQQ(3) = 3
!print
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N


!---Read remaining species:  X-sections at 1-2-3 T_s
!---read in 1 to 3 X-sects per J-value (JJ)
        JJ = 3
!-- read new Xsection block
    3 continue
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
!print
      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N
        if (TIT_J1S .eq. 'endofJ') goto 1
!---try to add a new Xsect
    2 continue
       JJ = JJ+1
       LQ = 1
         if (JJ .gt. X_) call EXITC(' RD_XXX: X_ not large enough')
       TITLEJX(JJ) = TIT_J1S
       TITLEJL(JJ) = TIT_J1L
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(JJ),TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
!try to read a 2nd Temperature or Pressure
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
!print
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N
        if (TIT_J1S .eq. 'endofJ') goto 1
      if (TIT_J1S .eq. TITLEJX(JJ)) then
        LQ = 2
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
        TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
!try to read a 3rd Temperature or Pressure
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
!print
         if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N
         if (TIT_J1S .eq. 'endofJ') goto 1
       if (TIT_J1S .eq. TITLEJX(JJ)) then
        LQ = 3
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
        TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
       else
        goto 2
       endif
      else
        goto 2
      endif
      goto 3
    4 continue
        call EXITC(' RD_XXX: error in read')
    1 continue
       NJX = JJ

!---read in complete, process Xsects for reduced wavelengths (Trop-Only)
!---    possibly also for WACCM >200nm-only version.
!---EmChem family of chemistry scheme; drop all other xsects to save CPU time
       if (CM_schemes_ChemRates(:7) .eq. " EmChem" .and. MasterProc) then 
        write(*,*)  &
        ' >>> Photolysis calculations only for EmChem cross-sections/J-values.'
      elseif (MasterProc) then
        write(*,*)  &
        ' >>> IMPORTANT: Photolysis calculations for all available cross-sections/J-values.'
      endif

      JJ = 0
      do J = 1,NJX
        ! include only 'e' and 'p' when EmChem family of chemistry schemes is used
        if (CM_schemes_ChemRates(:7) .ne. " EmChem" .or. SQQ(J) .eq. 'e' .or. SQQ(J) .eq. 'p') then
         ! if (SQQ(J) .eq. 'e' .or. SQQ(J) .eq. 'p') then
!---------collapse Xsects
          JJ = JJ+1
          if (JJ .lt. J) then
            TITLEJX(JJ) = TITLEJX(J)
            LQQ(JJ) = LQQ(J)
            SQQ(JJ) = SQQ(J)
            do LQ = 1,LQQ(J)
              TQQ(LQ,JJ) = TQQ(LQ,J)
              do IW = 1,NWWW
                QQQ(IW,LQ,JJ) = QQQ(IW,LQ,J)
              enddo
            enddo
          endif
        endif
      enddo
      NJX = JJ

!print-----
      do J = 1,NJX
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(a8,i5,2x,a6,2x,a16,2x,a1,i3,2x,3f6.1)') &
          ' X-sects',J,TITLEJX(J),TITLEJL(J),SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
      enddo

!---need to check that TQQ (= T(K) or p(hPa)) is monotonically increasing:
      do J = 1,NJX
        if ((LQQ(J).eq.3) .and. (TQQ(2,J).ge.TQQ(3,J))) then
            call EXITC ('TQQ out of order')
        endif
        if ((LQQ(J).eq.2) .and. (TQQ(1,J).ge.TQQ(2,J))) then
            call EXITC ('TQQ out of order')
        endif
      enddo

!---now collapse all the wavelengths for TROP-ONLY (W_ = 12 or 8)
      if (W_ .eq. 12) then
        NW2 = 12
        do IW = 1,4
           WL(IW) = WL(IW+4)
           FL(IW) = FL(IW+4)
           QRAYL(IW) = QRAYL(IW+4)
          do K = 1,3
           QO2(IW,K) = QO2(IW+4,K)
           QO3(IW,K) = QO3(IW+4,K)
           Q1D(IW,K) = Q1D(IW+4,K)
          enddo
          do J = 4,NJX
           do K = 1,LQQ(J)
            QQQ(IW,K,J) = QQQ(IW+4,K,J)
           enddo
          enddo
        enddo
        do IW = 5,12
           WL(IW) = WL(IW+6)
           FL(IW) = FL(IW+6)
           QRAYL(IW) = QRAYL(IW+6)
          do K = 1,3
           QO2(IW,K) = QO2(IW+6,K)
           QO3(IW,K) = QO3(IW+6,K)
           Q1D(IW,K) = Q1D(IW+6,K)
          enddo
          do J = 4,NJX
           do K = 1,LQQ(J)
            QQQ(IW,K,J) = QQQ(IW+6,K,J)
           enddo
          enddo
        enddo
      endif

      if (W_ .eq. 8) then
        NW2 = 8
        do IW = 1,1
           WL(IW) = WL(IW+4)
           FL(IW) = FL(IW+4) * 2.d0
           QRAYL(IW) = QRAYL(IW+4)
          do K = 1,3
           QO2(IW,K) = QO2(IW+4,K)
           QO3(IW,K) = QO3(IW+4,K)
           Q1D(IW,K) = Q1D(IW+4,K)
          enddo
          do J = 4,NJX
           do K = 1,LQQ(J)
            QQQ(IW,K,J) = QQQ(IW+4,K,J)
           enddo
          enddo
        enddo
        do IW = 2,8
           WL(IW) = WL(IW+10)
           FL(IW) = FL(IW+10)
           QRAYL(IW) = QRAYL(IW+10)
          do K = 1,3
           QO2(IW,K) = QO2(IW+10,K)
           QO3(IW,K) = QO3(IW+10,K)
           Q1D(IW,K) = Q1D(IW+10,K)
          enddo
          do J = 4,NJX
           do K = 1,LQQ(J)
            QQQ(IW,K,J) = QQQ(IW+10,K,J)
           enddo
          enddo
        enddo
      endif

      close(NUN)

      END SUBROUTINE RD_XXX


!-----------------------------------------------------------------------
      subroutine RD_CLD(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX ver 7.3+
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NCC      Number of categories for cloud scattering phase functions
!     QCC      Cloud scattering phase functions
!     WCC      5 Wavelengths for supplied phase functions
!     PCC      Phase function: first 8 terms of expansion
!     RCC      Effective radius associated with cloud type
!     SCC      Single scattering albedo
!     DCC      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, JCC
      character*120 TITLE0
!      character*12 TITLCC(C_)   ! TITLAA: Title for scattering data
      character*12 TITLCCJ
      real*8     RCCJ,DCCJ

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

        read (NUN,'(a120)',err=4) TITLE0
!print----
          if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(a120)') TITLE0
        read (NUN,*)
        read (NUN,*)
      do J = 1,C_
        read (NUN,'(i4,1x,a12,1x,2f6.3,1x,a120)',err=4) &
         JCC,TITLCCJ,RCCJ,DCCJ,TITLE0
       if (JCC.gt.0) then
         TITLCC(J) = TITLCCJ
         RCC(J) = RCCJ
         DCC(J) = DCCJ
        do K = 1,5
         read (NUN,'(f4.0,f7.4,f7.4,7f6.3)',err=4) &
          WCC(K,J),QCC(K,J),SCC(K,J),(PCC(I,K,J),I=2,8)
          PCC(1,K,J) = 1.d0
        enddo
         NCC = J
!print----
          if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(i5,1x,a12,1x,7f7.3,1x,a80)')   &
           J,TITLCCJ,RCCJ,DCCJ,(QCC(K,J),K=1,5),TITLE0
       else
          goto 2
       endif
      enddo
          goto 2

    4 continue
        call EXITC(' RD_MIE: error in read')

    2 continue
        close(NUN)
!print----
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      END SUBROUTINE RD_CLD


!-----------------------------------------------------------------------
      subroutine RD_MIE(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols scattering data set for fast-JX ver 7.3+
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!     DAA      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K , JAA
      character*120 TITLE0
!      character*12 TITLAA(A_)   ! TITLAA: Title for scattering data    NEEDS to be in COMMON
      Character*12 TITLAAJ
      real*8   RAAJ, DAAJ

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

      read (NUN,'(a120)',err=4) TITLE0
!print----
          if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(a120)') TITLE0
        read (NUN,*)
        read (NUN,*)
      do J = 1,A_
        read (NUN,'(i4,1x,a12,1x,2f6.3,1x,a120)',err=4) &
         JAA,TITLAAJ,RAAJ,DAAJ,TITLE0
       if (JAA.gt.0) then
         TITLAA(J) = TITLAAJ
         RAA(J) = RAAJ
         DAA(J) = DAAJ
        do K = 1,5
         read (NUN,'(f4.0,f7.4,f7.4,7f6.3)',err=4) &
          WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1.d0
        enddo
         NAA = J
!print----
          if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(i5,1x,a12,1x,7f7.3,1x,a80)')   &
           J,TITLAAJ,RAAJ,DAAJ,(QAA(K,J),K=1,5),TITLE0
       else
          goto 2
       endif
      enddo
          goto 2

    4 continue
        call EXITC(' RD_MIE: error in read')
    2 continue
        close(NUN)

      END SUBROUTINE RD_MIE


!-----------------------------------------------------------------------
      subroutine RD_UM(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------UMich aerosol optical data for fast-JX (ver 6.1+)
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, L
      character*120 TITLE0
      character*20 TITLUM(33)   ! TITLUM: Title for U Michigan aerosol data set

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

      read (NUN,'(a)',err=4) TITLE0
!print----
          if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(a120)') TITLE0
      read(NUN,'(5x,10f5.0)',err=4) WMM

!---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
!---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
      do L=1,33
          read(NUN,'(a4)',err=4) TITLUM(L)
!---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
        do K=1,21
!---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 5=600nm, 6=1000nm
!---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
          read(NUN,'(18f9.5)',err=4)  ((UMAER(I,J,K,L),I=1,3),J=1,6)
        enddo
      enddo

      close(NUN)

!  collapse UM wavelengths, drop 550 nm
          WMM(4) = WMM(5)
          WMM(5) = WMM(6)
       do L=1,33
       do K=1,21
       do I=1,3
          UMAER(I,4,K,L) = UMAER(I,5,K,L)
          UMAER(I,5,K,L) = UMAER(I,6,K,L)
       enddo
       enddo
       enddo

!print----
        if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(7(i5,1x,a4))') (L,TITLUM(L), L=1,33)
      goto 2
    4 continue
        call EXITC(' RD_UM: error in read')
    2 continue

      END SUBROUTINE RD_UM


!-----------------------------------------------------------------------
      subroutine RD_PROF(NJ2,NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input T and O3 reference profiles 'atmos_std.dat'
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real*8  OFAC, OFAK

      character*78 TITLE0
!
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
!      write(6,'(1X,A)') TITLE0
      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
      enddo
      close (NJ2)

!  Extend climatology to 100 km
      OFAC = exp(-2.d5/5.d5)
      do I = 32,LREF
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          OREF(I,L,M) = OREF(31,L,M)*OFAK
        enddo
        enddo
      enddo
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,LREF
        TREF(I,L,M) = TREF(41,L,M)
      enddo
      enddo
      enddo

 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_PROF


!-----------------------------------------------------------------------
      subroutine RD_JS_JX(NUNIT,NAMFIL,TITLEJX,NJX)
!-----------------------------------------------------------------------
!  Read 'FJX_j2j.dat' that defines mapping of fast-JX J's (TITLEJX(1:NJX))
!    onto the CTM reactions:  react# JJ, named T_REACT, uses fast-JX's JVMAP
!    including scaling factor JFACTA
!-----------------------------------------------------------------------
!---mapping variables stored in  block /jvchem/JFACTA,JIND,NRATJ,JLABEL,JVMAP
!           real*8  JFACTA(JVN_)          integer JIND(JVN_), NRATJ
!           character*50 JLABEL(JVN_)     character*6  JVMAP(JVN_)
!     JFACTA    multiplication factor for fast-JX calculated J
!     JLABEL    label(*50) of J-value used in the main chem model
!     JVMAP     label(*6) of J-value used to match with fast-JX J's
!     NRATJ     number of Photolysis reactions in CTM chemistry, derived here
!                   NRATJ must be .le. JVN_
!-----------------------------------------------------------------------
      implicit none
!
      integer, intent(in)                    ::  NUNIT, NJX
      character(*), intent(in)               ::  NAMFIL
      character*6, intent(in),dimension(NJX) :: TITLEJX
      integer   J,JJ,K
      character*120 CLINE
      character*50 T_REACT
      character*6  T_FJX
      real*8 F_FJX

! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
! The chemistry code title describes fully the reaction (a50)
! Blank (unfilled) chemistry J's are unmapped
! The number NRATJ is the last JJ readin that is .le. JVN
!   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JVMAP(:) = '------'
      JFACTA(:) = 0.d0

      open (NUNIT,file=NAMFIL,status='old',form='formatted')

       read (NUNIT,'(a)') CLINE
         if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(a)') CLINE
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ .gt. JVN_) exit
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JVMAP(JJ) = T_FJX
        NRATJ = JJ
      enddo

      close(NUNIT)

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
       do J = 1,NJX
        if (JVMAP(K) .eq. TITLEJX(J)) then
         JIND(K) = J
        endif
       enddo
      enddo

      if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(a,i4,a)')' Photochemistry Scheme with',NRATJ,' J-values'
      do K=1,NRATJ
       if (JVMAP(K) .ne. '------' ) then
        J = JIND(K)
        if (J.eq.0) then
          if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(i5,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K), &
               ' no mapping onto fast-JX',JVMAP(K)
        else
          if (USES%CLOUDJVERBOSE .and. MasterProc) write(6,'(i5,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K), &
               ' mapped to FJX:',J,TITLEJX(J)
        endif
       endif
      enddo

      END SUBROUTINE RD_JS_JX


!-----------------------------------------------------------------------
      SUBROUTINE RANSET (ND,RAN4L,ISTART)
!-----------------------------------------------------------------------
!  generates a sequence of real*4 pseudo-random numbers RAN4L(1:ND)
!     program RAN3 from Press, based on Knuth
      implicit none
      integer, parameter ::  MBIG=1000000000
      integer, parameter ::  MSEED=161803398
      integer, parameter ::  MZ=0
      real*4 , parameter ::  FAC=1.e-9
      integer,intent(in)    :: ND
      real*4, intent(out)   :: RAN4L(ND)
      integer,intent(inout) :: ISTART
      integer :: MA(55),MJ,MK,I,II,J,K,INEXT,INEXTP
!---initialization and/or fix of ISEED < 0
        MJ = MSEED - abs(ISTART)
        MJ = mod(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        do I=1,54
          II = mod(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          if (MK.lt.MZ) then
            MK=MK+MBIG
          endif
          MJ = MA(II)
        enddo
        do K=1,4
         do I=1,55
           MA(I)=MA(I)-MA(1+MOD(I+30,55))
           if (MA(I) .lt. MZ) then
             MA(I) = MA(I)+MBIG
           endif
         enddo
        enddo
        INEXT = 0
        INEXTP = 31
        ISTART = 1
!---generate next ND pseudo-random numbers
      do J=1,ND
         INEXT = mod(INEXT,55) +1
         INEXTP = mod(INEXTP,55) +1
         MJ = MA(INEXT) - MA(INEXTP)
        if (MJ .lt. MZ) then
          MJ=MJ+MBIG
        endif
         MA(INEXT) = MJ
         RAN4L(J) = MJ*FAC
      enddo

      END SUBROUTINE RANSET


      END MODULE FJX_INIT_MOD

!------------------------------------------------------------------------------
!          UCI fast-JX  cloud-JX v-7.3d (1/2016)                                        !
!------------------------------------------------------------------------------
!
! !DESCRIPTION: decides what to do with cloud fraction,
!    including generate Independent Column Atmospheres (ICAs)for a max-ran
!    cloud overlap algorithm, and Quadrature Colm Atmos (QCAs).
! v 7.3b corrected an indexing/segmentation error in IG2 affecting L_CLR2
!    no change in results, L_CLR2 = t or f does not affect top layer.
! v 7.3c added GLVL and reduced  cloud corr factor when there were gaps in levels
! v 7.3d  minor change in sub ICA_NR() with ZZZ(L) ==> ZZZ(L)-ZZZ(1) for Zbin categories
!             if (ZZZ(L)-ZZZ(1) .lt. Zbin(N)) then
!
      MODULE CLD_SUB_MOD

! USES:
      USE FJX_CMN_MOD

      USE FJX_SUB_MOD,  ONLY : EXITC, PHOTO_JX

      use Config_module, only: cloudjx_initf, USES, MasterProc

      IMPLICIT NONE
!-----------------------------------------------------------------------
!
! !PUBLIC SUBROUTINES:
!
      PUBLIC  :: CLOUD_JX

      CONTAINS


      SUBROUTINE CLOUD_JX (U0,SZA,REFLB,SOLF,FG0,LPRTJ,PPP,ZZZ,TTT, &
             DDD,RRR,OOO,   LWP,IWP,REFFL,REFFI, CLDF,CLDCOR,CLDIW, &
             AERSP,NDXAER,L1U,ANU,VALJXX,NJXU,                      &
             CLDFLAG,NRANDO,IRAN,LNRG,NICA,JCOUNT)

      implicit none


!---Newest recommended approach (v7.3) to use cloud correlation lengths
!   Problem with correlation is that each new cloud layer generates 2x combinations
!      Thus the possibilities are for 2**Lcloudtop ICAs - that is too many (2**35)
!   Using  correlation lengths (3 km for tropical and high clouds, 1.5 km for stratus)
!      Choose 6 bins (of thickness = correl length) as Max-Overlap, then have
!      these bins be randomly or correlated with bins above
!   For now just assume these are random as with the other 2 max-ran groups above.
!
! GRP1 = 0 - 1.5km,  GRP2 = 1.5 - 3.5km,  GRP3 = 3.5 - 6km
! GRP4 =  6 - 9km,   GRP5 = 9 - 13km,     GRP6 = 13km+   (GRP7 = separate cirrus shields)
!
!Key Refs
!  Kato, S., et al (2010), Relationships among cloud occurrence frequency, overlap,
!        and effective thickness derived from CALIPSO and CloudSat merged cloud vertical
!        profiles, J. Geophys. Res., 115, D00H28, doi:10.1029/2009JD012277.
! Pincus, R., et al. (2005), Overlap assumptions for assumed probability distribution
!        function cloud schemes in large-scale models, J. Geophys. Res., 110, D15S09,
!        doi:10.1029/2004JD005100.
! Oreopoulos, L., et al (2012) Radiative impacts of cloud heterogeneity and overlap
!        in an atmospheric General Circulation Model,  Atmos. Chem. Phys., 12, 90979111,
!        doi:10.5194/acp-12-9097-2012
! Naud, C., & A. D.  DelGenio (2006) Cloud Overlap Dependence on Atmospheric Dynamics,
!        16th ARM Science Team Meeting Proceedings, Albuquerque, NM, March 27 - 31, 2006.


!  CLOUD_JX is fractional cloud cover driver for PHOTO_JX  (fast-JX v7.2)
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-fo-year for sun distance, SZA (not lat or long)

!--CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
!       CLDFLAG = 1  :  Clear sky J's
!       CLDFLAG = 2  :  Averaged cloud cover
!       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
!       CLDFLAG = 4  :  Average direct solar beam over all ICAs, invert to get clouds
!       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
!       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
!       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
!       CLDFLAG = 8  :  Calcluate J's for ALL ICAs (up to 20,000 per cell!)

!--CLDIW = index for each lcoud layer:
!     = 0 = no cloud
!     = 1 = water cloud only
!     = 2 = ice cloud only
!     = 3 = liquid+ice cloud mix
!-----------------------------------------------------------------------

!---calling sequence variables
      integer, intent(in)                    :: L1U,ANU,NJXU, CLDFLAG,IRAN,NRANDO,LNRG
      real*8,  intent(in)                    :: U0,SZA,REFLB,SOLF,FG0, CLDCOR
      logical, intent(in)                    :: LPRTJ
      real*8,  intent(in), dimension(L1U+1)  :: PPP,ZZZ
      real*8,  intent(in), dimension(L1U  )  :: TTT,DDD,RRR,OOO,  &
                                                LWP,IWP,REFFL,REFFI
      real*8,  intent(in), dimension(L1U,ANU):: AERSP
      integer, intent(in), dimension(L1U,ANU):: NDXAER
      real*8,  intent(in), dimension(L1U  )  :: CLDF
      integer, intent(in), dimension(L1U  )  :: CLDIW
! reports out the JX J-values, upper level program converts to CTM chemistry J's
      real*8, intent(out), dimension(L1U-1,NJXU)::  VALJXX
      integer, intent(out) :: NICA,JCOUNT

!-----------------------------------------------------------------------

      logical LPRTJ0, LDARK
      integer I,II,J,L,N, LTOP
      real*8, dimension(L1U)  :: LWPX,IWPX,REFFLX,REFFIX
      real*8  CLDFR, XRAN, FSCALE, QCAOD, WTRAN, G0LIQ,G0ICE, SUM_W,SUM_O

      real*8, dimension(LWEPAR) ::  CLTL,CLTI, CLT,CLDX
      integer                   ::  NRG, IRANX
! max number of max-overlap cloud groups set at 9
      integer, dimension(LWEPAR) :: NCLDF
      integer, dimension(9) :: GBOT,GTOP,GLVL,GNR,GCMX
      integer, dimension(9,CBIN_+1) :: GFNR
      real*8, dimension(CBIN_)  ::  CFBIN
      real*8, dimension(ICA_) ::    WCOL,OCOL, OCDFS
!      real*8, dimension(LWEPAR,ICA_) :: TCOL
      integer, dimension(ICA_) :: ISORT
      real*8, dimension(LWEPAR+1)   :: TCLD,TTCOL,SOLT,TAUG
      real*8, dimension(L1U-1,NJXU)::  VALJXXX
      real*8,  dimension(NQD_)     :: WTQCA
      integer, dimension(NQD_)     :: NQ1,NQ2,NDXQS
      
!-----------------------------------------------------------------------
      LPRTJ0 = LPRTJ
      JCOUNT = 0
      NICA = 0
      do L = LWEPAR+1, L1U
        LWPX(L) = 0.d0
        IWPX(L) = 0.d0
        REFFLX(L) = 0.d0
        REFFIX(L) = 0.d0
      enddo

!---CLOUD_JX:   different cloud schemes
!-----------------------------------------------------------------------
      if (CLDFLAG.lt.1 .or. CLDFLAG.gt.8) &
         call EXITC ('>>>stop, incorrect cloud index') 

      if (CLDFLAG.le.3) then

!-----------------------------------------------------------------------
       if (CLDFLAG.eq.2) then
! 2 = average cloud cover
         do L = 1, LWEPAR
           CLDFR = CLDF(L)
           LWPX(L) = LWP(L) * CLDFR
           IWPX(L) = IWP(L) * CLDFR
           REFFLX(L) = REFFL(L)
           REFFIX(L) = REFFI(L)
         enddo

!-----------------------------------------------------------------------
       elseif (CLDFLAG.eq.3) then
! 3 = average cloud cover, adjust cloud fraction **3/2
         do L = 1, LWEPAR
           CLDFR = CLDF(L) * sqrt(CLDF(L))
           LWPX(L) = LWP(L) * CLDFR
           IWPX(L) = IWP(L) * CLDFR
           REFFLX(L) = REFFL(L)
           REFFIX(L) = REFFI(L)
         enddo

!-----------------------------------------------------------------------
       elseif (CLDFLAG.eq.1) then
! 1 = clear sky - no clouds
         do L = 1, LWEPAR
           LWPX(L) = 0.d0
           IWPX(L) = 0.d0
         enddo
       endif
!-----------------------------------------------------------------------

!----all above have only a single, simple call for fast_JX------------
         if(LPRTJ0) then
         write(6,'(2a)') ' cloud_JX (7.3) Internal print: clouds = ',&
                         TITCLD(CLDFLAG)
         endif
!-----------------------------------------------------------------------
      call PHOTO_JX (U0,SZA,REFLB,SOLF, LPRTJ0, PPP,ZZZ,TTT,  &
                  DDD,RRR,OOO, LWPX,IWPX,REFFLX,REFFIX,   &
                  AERSP,NDXAER, L1U,ANU, VALJXX,NJXU, LDARK)
        if (.not.LDARK) JCOUNT = JCOUNT + 1
!-----------------------------------------------------------------------

      else
!-----------------------------------------------------------------------
!  All below = CLDFLAG = 4:8  need to set up cloud overlap algorithm
!-----------------------------------------------------------------------
         do L = 1, LWEPAR
           CLDX(L) = CLDF(L)
           CLT(L) = 0.d0
           CLTI(L) = 0.d0
           CLTL(L) = 0.d0
           LWPX(L) = LWP(L)
           IWPX(L) = IWP(L)
           REFFLX(L) = REFFL(L)
           REFFIX(L) = REFFI(L)
         enddo
! do cloud fraction binning here and rescale IWP/LWP to conserve layer WP


! generate approx cloud visible optical depths for quadrature and sorting
!   true wavelength dependence will be recalculated in PHOTO_JX
         do L = 1,LWEPAR
            if (REFFIX(L) .gt. 0.d0) then
              CLTI(L) = IWPX(L)*0.75d0*2.d0/(REFFIX(L)*0.917d0)
              CLT(L) = CLT(L) + CLTI(L)
            endif
            if (REFFLX(L) .gt. 0.d0) then
              CLTL(L) = LWPX(L)*0.75d0*2.1d0/REFFLX(L)
              CLT(L) = CLT(L) + CLTL(L)
            endif
         enddo
         LTOP  = LWEPAR
!-------------------------------------------------------------------------
!---Generate max-ran cloud overlap groups used for CLDFLAG = 4:8
!---CLT(cloud ice+liq OD) & IWPX & LWPX adjusted to quantized cld fr
!-------------------------------------------------------------------------
      call ICA_NR(CLDX,CLT,IWPX,LWPX,ZZZ, CLDIW,LTOP,LNRG,CBIN_,ICA_, &
             CFBIN,CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA)


!---call ICA_ALL to generate the weight and cloud total OD of each ICA
!-------------------------------------------------------------------------
      call ICA_ALL(CLDX,CLT,LTOP,CBIN_,ICA_, CFBIN,     &
          CLDCOR,NCLDF,GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA,  WCOL,OCOL)


         if(LPRTJ0) then
         write(6,*) ' cloud-JX (7.3) internal print:  #ICAs = ',NICA
         endif

!-----------------------------------------------------------------------
! 4 = average direct beam over all ICAs, est. isotropic equiv for liq/ice
       if (CLDFLAG .eq. 4) then
        G0LIQ= min(0.96d0, FG0*PCC(2,3,3)/3.d0)
        G0ICE= min(0.96d0, FG0*PCC(2,3,7)/3.d0)
        call ICA_DIRECT(U0,CLDX,CLT,CLTL,CLTI, LTOP,CBIN_,ICA_,       &
                CLDCOR,NCLDF,GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA,   &
                          WCOL,TCLD,TTCOL,SOLT,G0LIQ,G0ICE,TAUG)
          if(LPRTJ0) then
           write(6,'(a,2f8.4)') ' Asymm factor g0 for Liq & Ice:', &
                                G0LIQ,G0ICE
           write(6,*)  &
               ' Average Direct beam: cld OD, cld FR - inferred OD'
           write(6,'(i5,3f10.4)') (L, CLT(L),CLDX(L),TCLD(L), L=1,LTOP)
          endif

! using TCLD as the eff-cloud OD in each layer, repartition bewteen ice & liq
         do L=1,LTOP
          if (TCLD(L).gt.1.d-7 .and. CLT(L).gt.1.d-7) then
           FSCALE = TCLD(L) / CLT(L)
           IWPX(L) = IWPX(L) * FSCALE
           LWPX(L) = LWPX(L) * FSCALE
          else
           IWPX(L) = 0.d0
           LWPX(L) = 0.d0
          endif
         enddo
         call PHOTO_JX (U0,SZA,REFLB,SOLF, LPRTJ0, PPP,ZZZ,TTT,  &
                       DDD,RRR,OOO, LWPX,IWPX,REFFLX,REFFIX,   &
                       AERSP,NDXAER, L1U,ANU, VALJXX,NJXU, LDARK)
          if (.not.LDARK) JCOUNT = JCOUNT + 1
       endif

!-----------------------------------------------------------------------
! 5 = random pick of NRANDO(#) ICAs (selected based on fractional area)
       if (CLDFLAG .eq. 5) then

          if(LPRTJ) then
          write(6,*) ' Average Js over random ICAs:',NRANDO
          endif

          VALJXXX(:,:) = 0.d0
          WTRAN = 1.d0/float(NRANDO)
          OCDFS(1) = WCOL(1)
         do I = 2,NICA
          OCDFS(I) = OCDFS(I-1) + WCOL(I)
         enddo
        do N=1,NRANDO
          IRANX = mod (IRAN+N-1, NRAN_) + 1
          XRAN = RAN4(IRANX)
            I = 1
         do while (XRAN .gt. OCDFS(I) .and. I .lt. NICA)
            I = I+1
         enddo
          call ICA_III(CLDX,CLT,LTOP,CBIN_,ICA_, I, &
              CLDCOR,NCLDF,GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)
!---zero out cloud water paths which are not in the selected random ICA
         do L = 1, LTOP
           LWPX(L) = LWP(L)
           IWPX(L) = IWP(L)
         enddo
         do L = 1,LTOP
          if(TTCOL(L) .lt. 1.d-8) then
           IWPX(L) = 0.d0
           LWPX(L) = 0.d0
          endif
         enddo

          if(LPRTJ) then
           write(6,'(a,2i6,f8.3)') ' pick random ICA:',N,I,OCOL(I)
           do L = 1,LTOP
            if (TTCOL(L) .ge. 1.d-8) then
             write(6,'(i5,f9.4,2f8.3)') L,TTCOL(L),LWPX(L),IWPX(L)
            endif
           enddo
          endif

         call PHOTO_JX (U0,SZA,REFLB,SOLF, LPRTJ0, PPP,ZZZ,TTT,  &
                       DDD,RRR,OOO, LWPX,IWPX,REFFLX,REFFIX,   &
                       AERSP,NDXAER, L1U,ANU, VALJXX,NJXU, LDARK)
          if (.not.LDARK) JCOUNT = JCOUNT + 1
          LPRTJ0 = .false.
         do J = 1,NJXU
          do L = 1,L1U-1
            VALJXXX(L,J) = VALJXXX(L,J) + VALJXX(L,J)*WTRAN
          enddo
         enddo
        enddo
         do J = 1,NJXU
          do L = 1,L1U-1
            VALJXX(L,J) = VALJXXX(L,J)
          enddo
         enddo
       endif


!-----------------------------------------------------------------------
! 6 = calculate qudrature QCAs, use all 4 mid points
       if (CLDFLAG .eq. 6) then
         call ICA_QUD(WCOL,OCOL,LTOP,ICA_,NQD_,NICA, &
                         WTQCA, ISORT,NQ1,NQ2,NDXQS)

          if (LPRTJ0) then
           write(6,*) ' quadrature QCAs(avg): wt/range/index/OD'
           do N=1,NQD_
            if (WTQCA(N).gt.0.d0) then
              write(6,'(i5,f8.4,3i8,2f10.3)')   &
               N,WTQCA(N),NQ1(N),NQ2(N),NDXQS(N),OCOL(ISORT(NDXQS(N)))
            endif
           enddo
          endif

           VALJXXX(:,:) = 0.d0
        do N = 1, NQD_
        if (WTQCA(N) .gt. 0.d0) then
          I = ISORT(NDXQS(N))
          call ICA_III(CLDX,CLT,LTOP,CBIN_,ICA_, I, &
              CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)
!---zero out cloud water paths which are not in the selected QCA
         do L = 1, LTOP
           LWPX(L) = LWP(L)
           IWPX(L) = IWP(L)
         enddo
         do L = 1,LTOP
          if (TTCOL(L) .lt. 1.d-8) then
           IWPX(L) = 0.d0
           LWPX(L) = 0.d0
          endif
         enddo
         call PHOTO_JX (U0,SZA,REFLB,SOLF, LPRTJ0, PPP,ZZZ,TTT,  &
                        DDD,RRR,OOO, LWPX,IWPX,REFFLX,REFFIX,   &
                        AERSP,NDXAER, L1U,ANU, VALJXX,NJXU, LDARK)
          if (.not.LDARK) JCOUNT = JCOUNT + 1
          LPRTJ0 = .false.
         do J = 1,NJXU
          do L = 1,L1U-1
            VALJXXX(L,J) = VALJXXX(L,J) + VALJXX(L,J)*WTQCA(N)
          enddo
         enddo
        endif
        enddo
!---
         do J = 1,NJXU
          do L = 1,L1U-1
            VALJXX(L,J) = VALJXXX(L,J)
          enddo
         enddo
       endif


!-----------------------------------------------------------------------
! 7 = calculate quadrature atmosphere - average cloud within each QCA bin.
       if (CLDFLAG .eq. 7) then

         call ICA_QUD(WCOL,OCOL,LTOP,ICA_,NQD_,NICA, &
                         WTQCA, ISORT,NQ1,NQ2,NDXQS)

           VALJXXX(:,:) = 0.d0
        do N = 1, NQD_
         if (WTQCA(N) .gt. 0.d0) then
         if (NQ2(N) .ge. NQ1(N)) then
            IWPX(:) = 0.d0
            LWPX(:) = 0.d0
            QCAOD = 0.d0
          do II = NQ1(N),NQ2(N)
            I = ISORT(II)
           call ICA_III(CLDX,CLT,LTOP,CBIN_,ICA_, I, &
              CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)

             if (LPRTJ) then
              write(6,'(a,3i5,2f8.4,f9.3)') ' N(QCA)/II/I WCOL,OCOL',  &
                        N,II,I,WCOL(I),WTQCA(N),OCOL(I)
             endif

           do L = 1,LTOP
            if (TTCOL(L) .gt. 1.d-8) then
             IWPX(L) = IWPX(L) + IWP(L)*WCOL(I)
             LWPX(L) = LWPX(L) + LWP(L)*WCOL(I)
             QCAOD = QCAOD + TTCOL(L)*WCOL(I)
            endif
           enddo
          enddo
          do L = 1,LTOP
            IWPX(L) = IWPX(L)/WTQCA(N)
            LWPX(L) = LWPX(L)/WTQCA(N)
          enddo
            QCAOD = QCAOD/WTQCA(N)

           if (LPRTJ) then
            write(6,'(a,i3,a,f10.5,f10.3)') &
               'Quad Atmos Avg #',N,' wt, tot-OD:',WTQCA(N),QCAOD
            write(6,'(a)') 'L / LWP / IWP'
            do L=1,LTOP
            if (LWPX(L)+IWPX(L) .gt. 1.d-8)  &
              write(6,'(i4,2f10.3)') L,LWPX(L),IWPX(L)
            enddo
           endif

           call PHOTO_JX (U0,SZA,REFLB,SOLF, LPRTJ0, PPP,ZZZ,TTT,  &
                          DDD,RRR,OOO, LWPX,IWPX,REFFLX,REFFIX,   &
                          AERSP,NDXAER, L1U,ANU, VALJXX,NJXU, LDARK)
          if (.not.LDARK) JCOUNT = JCOUNT + 1
          LPRTJ0 = .false.
          do J = 1,NJXU
           do L = 1,L1U-1
             VALJXXX(L,J) = VALJXXX(L,J) + VALJXX(L,J)*WTQCA(N)
           enddo
          enddo
         endif
         endif
        enddo
!---
         do J = 1,NJXU
          do L = 1,L1U-1
            VALJXX(L,J) = VALJXXX(L,J)
          enddo
         enddo
       endif


!-----------------------------------------------------------------------
! 8 = average J's over all ICAs
       if (CLDFLAG .eq. 8) then

          if(LPRTJ) then
           write(6,*) ' Average Js over all ICAs: I/ODcol/WTcol'
           write(6,'(i5,2f9.4)') (L,OCOL(L),WCOL(L), L=1,min(12,NICA-1))
           if (NICA.gt.12) write(6,'(a)') '. . .'
           write(6,'(i5,2f9.4)') NICA,OCOL(NICA),WCOL(NICA)
          endif

           VALJXXX(:,:) = 0.d0
        do I = 1, NICA
          call ICA_III(CLDX,CLT,LTOP,CBIN_,ICA_, I, &
              CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)
!---zero out cloud water paths which are not in the selected random ICA
         do L = 1, LTOP
           LWPX(L) = LWP(L)
           IWPX(L) = IWP(L)
         enddo
         do L = 1,LTOP
          if(TTCOL(L) .lt. 1.d-8) then
           IWPX(L) = 0.d0
           LWPX(L) = 0.d0
          endif
         enddo
         call PHOTO_JX (U0,SZA,REFLB,SOLF, LPRTJ0, PPP,ZZZ,TTT,  &
                       DDD,RRR,OOO, LWPX,IWPX,REFFLX,REFFIX,   &
                       AERSP,NDXAER, L1U,ANU, VALJXX,NJXU, LDARK)
          if (.not.LDARK) JCOUNT = JCOUNT + 1
          LPRTJ0 = .false.
         do J = 1,NJXU
          do L = 1,L1U-1
            VALJXXX(L,J) = VALJXXX(L,J) + VALJXX(L,J)*WCOL(I)
          enddo
         enddo
        enddo
         do J = 1,NJXU
          do L = 1,L1U-1
            VALJXX(L,J) = VALJXXX(L,J)
          enddo
         enddo
       endif
!-----------------------------------------------------------------------

      endif

      END SUBROUTINE CLOUD_JX


!-----------------------------------------------------------------------
      SUBROUTINE ICA_NR(CLDF,CLTAU,IWPX,LWPX,ZZZ,CLDIW,LTOP,LNRG,CBIN_, &
            ICA_,CFBIN,CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA)
!-----------------------------------------------------------------------
!---Read in the cloud fraction (CLDF), cloud OD (CLTAU), cloud index (CLDIW)

!---Derive max-ran cloud overlaps and set up all the ICAs (Independent Column
!          Atmos)
!   NCLDF(L) has value 1:CBIN_+1 = quantized cloud fraction (0:CBIN_). Value of 0 means NO
!          clouds
!   CFBIN(J) = cloud fraction assumed for bin L=0:CBIN_
!    (e.g., 1(=0-5%) = 0.025)
!   CLTAU(J) = is readjusted for quantum bins to preserve CLDF*CLTAU

!---Definition of ICAs is carefully laid out with key parameters below:
!-----------------------------------------------------------------------
!   NICA = no. of ICAs
!   NRG = no. of sub-groups that are randomly overlapped with each other,
!     but are maximally overlapped among the contiguous layers in the sub-group.
!         Technically, NRG can equal LTOP generating up to 2**LTOP ICAs.
!   GBOT(G=1:NRG) = lower CTM layer of max-overlap group G
!   GTOP(G=1:NRG) = upper CTM layer of max-overlap group G
!          All layers, cloudy or clear are placed in one NRG group.
!   GNR(G=1:NRG) = no. of unique quantized cloud fractions in group G
!                  (.le.CBIN_)
!          Defines the number of uniques fractions in a maximally overlapped
!          group.
!   GFNR(G=1:NRG,1:GNR(G)) = cloud fraction quantum no (value = 0 to NCBIN)
!          Stores the specific cloud fractions counted in GNR.
!
!--CLDIW = index for each lcoud layer:
!     = 0 = no cloud
!     = 1 = water cloud only
!     = 2 = ice cloud only
!     = 3 = liquid+ice cloud mix

!-----------------------------------------------------------------------
      implicit none

!---Cloud Cover parameters
!      integer, parameter ::  NQD_  = 4
!      integer, parameter ::  NRAN_ = 10007  ! dimension for random number
!      integer, parameter ::  CBIN_ = 20     ! # of quantized cld fration bins
!  may need to reduce quantum number for LNRG6 to be CBIN_ = 10
!      integer, parameter ::  ICA_  = 20000  ! # of Indep Colm Atmospheres
!---Local Cloud Cover parameters
!---define break between randomly overlapped groups
      integer, parameter  :: NG_BRK = 0
!---set up for correlated max-overlap groups based on observations
      integer, parameter ::  NRG6_ = 6
      real*8, dimension(NRG6_), parameter:: Zbin =                 &
          [0.d5, 1.5d5, 3.5d5, 6.0d5, 9.0d5, 13.d5]

      integer,intent(in) :: LTOP, LNRG, CBIN_, ICA_
      integer,intent(in),dimension(LTOP) :: CLDIW
      real*8, intent(in),dimension(LTOP) :: CLDF, ZZZ
      real*8, intent(in)                 :: CLDCOR
      real*8, intent(inout),dimension(LTOP) :: CLTAU,IWPX,LWPX

      integer, intent(out) ::  NRG, NICA
      integer, intent(out), dimension(LTOP) :: NCLDF
      integer, intent(out), dimension(9) :: GBOT,GTOP,GLVL,GNR,GCMX
      integer, intent(out), dimension(9,CBIN_+1) :: GFNR
      real*8,  intent(out), dimension(CBIN_) :: CFBIN

      real*8   FBIN, FSCALE, CLF_MIN, CLF_MAX, FSCALE2
      integer                   ::  NRGX, NICAX
      integer, dimension(9)  :: GMIN,GMAX
      integer, dimension(CBIN_) :: NSAME
      integer  I,K,L,LL,N,NC, L1,L2,L3,  LCLTOP,LCIRRUS
      logical  L1GRP,L2GRP,L3GRP, L6GRP
!-----------------------------------------------------------------------

!---quantize cloud fractions into bins to avoid excessive calculations for
!---  nearly identical maximally overlapping cloud fractions.
!---  CBIN_=20 => N=0=[0%], N=1=[0.001-5%],N=2=[5-10%], .
!                 N=19=[90-95%],N=20=[95-100%]
!---assume the upper end of the range and renormalize in-cloud TAU to preserve
!     CLDF*TAU
        FBIN = CBIN_
      do K = 1,CBIN_
        CFBIN(K) = float(K) / FBIN
      enddo

!---quantize cloud fractions into bins & adjust to conserve TAU*CF
! round up/down <2%=>0% and >98%=>100%  works for FBIN = 10 or 25 (not 40)
      CLF_MIN = 0.02d0
      CLF_MAX = 1.0d0 - CLF_MIN

! clear out any small fraction clouds and quantize the cloud fraction as NCLDF
      do L = 1,LTOP
       if (CLDF(L).lt.CLF_MIN .or. CLDIW(L).eq.0) then
         NCLDF(L) = 0
         CLTAU(L) = 0.d0
         IWPX(L) = 0.d0
         LWPX(L) = 0.d0
       elseif (CLDF(L) .gt. CLF_MAX) then
         NCLDF(L) = FBIN
       else
          FSCALE2 = FBIN
          FSCALE2 = CLDF(L)*FSCALE2 + 0.4999d0
         NCLDF(L) = max(FSCALE2,1.00001d0)
       endif
      enddo

!---find the cloud-top layer (water or ice) or identify clear sky
         LCLTOP = 0
      do L = 1,LTOP
        if (NCLDF(L).gt.0) then
          LCLTOP = L
        endif
      enddo
      if (LCLTOP .eq. 0) then
          NRG = 0
          NICA = 1
        goto 1
      endif

! rescale LWPX, IWPX, CLTAU
      do L = 1,LCLTOP
        if (NCLDF(L) .gt. 0) then
          FSCALE = CLDF(L) / CFBIN(NCLDF(L))
          CLTAU(L) = CLTAU(L) * FSCALE
          IWPX(L) = IWPX(L) * FSCALE
          LWPX(L) = LWPX(L) * FSCALE
        endif
      enddo

!---define maximally overlapping sub-groups by set levels (LNRG) or min
!   cloud fraction
      if (LNRG .eq. 0) then
!-----------------------------------------------------------------------------
!---Identify the maximally overlapped groups by breaking at a minimun NCLDF(L)
!-----------------------------------------------------------------------------
!---search from bottom to top, finding 1st level in group with cloud fraction
!    .ge. threshold, and then first level above that at which the cld fraction
!    is .lt. threshold. NRG = number of such groups.
        L = 1
        NRG = 0
        do while (L.lt.LCLTOP)
          if (NCLDF(L) .gt. NG_BRK) then
            NRG = NRG+1
           if (NRG.gt.9) exit
            GMIN(NRG) = L
            GMAX(NRG) = LCLTOP
            do LL = L+1,LCLTOP
!  look for first layer to drop below CLDF threshold = NGRBRK
              if (NCLDF(LL) .le. NG_BRK) then
                GMAX(NRG) = LL
                exit
              endif
            enddo
            L = GMAX(NRG)+1
          else
            L = L+1
          endif
        enddo

      elseif (LNRG .eq. 3) then
!-----------------------------------------------------------------------------
!---Alternative approach to fix a maximum of 3 random-overlap groups
!    (for L60/L57 CTM)
!---GRP=1 (if at all) is L=1:8 (surf to +1 km)
!---GRP=2 (if at all) is L=9 to last LWCloud
!---GRP=3 (fi at all) is L=last-LWCld+1 to LTOP
!-----------------------------------------------------------------------------
        L1 = 1
        L2 = 9
!----- L3-1 = uppermost liquid water cloud,  L3 = first of only ice-clouds
        L3 = L2
        do L = LCLTOP,L2,-1
          if (CLDIW(L).eq.1 .or. CLDIW(L).eq.3) then
            L3 = L+1
            exit
          endif
        enddo

           L1GRP = .false.
           L2GRP = .false.
           L3GRP = .false.
        do L = L1,L2-1
           L1GRP = L1GRP .or. (NCLDF(L).gt.0)
        enddo
        do L = L2,L3-1
           L2GRP = L2GRP .or. (NCLDF(L).gt.0)
        enddo
        do L = L3,LCLTOP
           L3GRP = L3GRP .or. (NCLDF(L).gt.0)
        enddo
          NRG = 0
        if (L1GRP) then
          NRG = NRG+1
          GMIN(NRG) = L1
          GMAX(NRG) = L2-1
        endif
        if (L2GRP) then
          NRG = NRG+1
          GMIN(NRG) = L2
          GMAX(NRG) = L3-1
        endif
        if (L3GRP) then
          NRG = NRG+1
          GMIN(NRG) = L3
          GMAX(NRG) = LCLTOP
        endif
           NRG = max(NRG,1)
          GMIN(1)   = 1
          GMAX(NRG) = LCLTOP

      else
!-----------------------------------------------------------------------------
!---Newest recommended approach (v7.3) to use cloud correlation lengths
!   Problem with correlation is that each new cloud layer generates 2x combinations
!      Thus the possibilities are for 2**Lcloudtop ICAs - that is too many (2**35)
!   Using  correlation lengths (3 km for tropical and high clouds, 1.5 km for stratus)
!      Choose 6 bins (of thickness = correl length) as Max-Overlap, then have
!      these bins be randomly or correlated with bins above
!   For now just assume these are random as with the other 2 max-ran groups above.
!
! GRP1 = 0 - 1.5km,  GRP2 = 1.5 - 3.5km,  GRP3 = 3.5 - 6km
! GRP4 =  6 - 9km,   GRP5 = 9 - 13km,     GRP6 = 13km+
!
!Key Refs
!  Kato, S., et al (2010), Relationships among cloud occurrence frequency, overlap,
!        and effective thickness derived from CALIPSO and CloudSat merged cloud vertical
!        profiles, J. Geophys. Res., 115, D00H28, doi:10.1029/2009JD012277.
! Pincus, R., et al. (2005), Overlap assumptions for assumed probability distribution
!        function cloud schemes in large-scale models, J. Geophys. Res., 110, D15S09,
!        doi:10.1029/2004JD005100.
! Oreopoulos, L., et al (2012) Radiative impacts of cloud heterogeneity and overlap
!        in an atmospheric General Circulation Model,  Atmos. Chem. Phys., 12, 90979111,
!        doi:10.5194/acp-12-9097-2012
! Naud, C., & A. D.  DelGenio (2006) Cloud Overlap Dependence on Atmospheric Dynamics,
!        16th ARM Science Team Meeting Proceedings, Albuquerque, NM, March 27 - 31, 2006.
!
!---Find the levels in each of the NRG6_  altitude-defined groups
       do L = 1,LCLTOP
        do N = 2,NRG6_
         if (ZZZ(L)-ZZZ(1) .lt. Zbin(N)) then
           GMAX(N-1) = L
         endif
        enddo
       enddo
         GMIN(1) = 1
       do N = 2,NRG6_
         GMIN(N) = GMAX(N-1) + 1
       enddo
         GMAX(NRG6_) = LCLTOP

!---find out if there are any clouds in each of the 6 correlated groups
          NRG = 0
       do N = 1,NRG6_
             L6GRP = .false.
        do L = GMIN(N),GMAX(N)
          if (NCLDF(L) .gt. 0) then
             L6GRP = .true.
          endif
        enddo
        if (L6GRP) then
           NRG = NRG + 1
           GMIN(NRG) = GMIN(N)
           GMAX(NRG) = GMAX(N)
           GLVL(NRG) = N
        endif
       enddo

!---pull off cirrus shields from top MAX-GRP as separate MAX-GRP to
!        allow better correlation between cumulus towers below them
       if (NRG .gt. 0) then
          LCIRRUS = 0
        do L = GMAX(NRG),GMIN(NRG),-1
         if (NCLDF(L).gt.FBIN/2 .and. CLDIW(L).eq.2) then
          LCIRRUS = L
         endif
        enddo
        if (LCIRRUS .gt.GMIN(NRG)) then
!---split the uppermost MAX-GRP
          GMIN(NRG+1) = LCIRRUS
          GMAX(NRG+1) = GMAX(NRG)
          GLVL(NRG+1) = 7
          GMAX(NRG) = LCIRRUS-1
          NRG = NRG+1
        endif
       endif

      endif
!---finished selection of max-overlap groups

!---simplify groups if no clouds with NCLDF > NG_BRK
      GBOT(:) = 0
      GTOP(:) = 0
      if (NRG .eq. 0) then
        NRG = 1
        GBOT(1) = 1
        GTOP(1) = LCLTOP
      else
!---assign levels between maximum overlap groups to group above.
        GBOT(1) = 1
        GTOP(1) = GMAX(1)
        do N=2,NRG
          GBOT(N) = max(GTOP(N-1)+1, GMIN(N))
!         GBOT(N) = min(GTOP(N-1)+1, GMIN(N))
          GTOP(N) = GMAX(N)
        enddo
        GTOP(NRG) = LCLTOP
      endif
!---for each max-overlap group calculate number of unique cloud fractions
      do N = 1,NRG
        NSAME(:) = 0
        GCMX(N) = 0
        do L = GBOT(N),GTOP(N)
          if (NCLDF(L) .gt. 0) then
            NSAME(NCLDF(L)) = 1
            GCMX(N) = max(GCMX(N),NCLDF(L))
          endif
        enddo
!---sort cloud fractions in deceasing order for each max-overlap group
!---  note that largest bin N=CBIN_ (eg, 95-100%) will be treated as 100%
        GFNR(N,1) = CBIN_
        NC = 1
        do I = CBIN_-1,1,-1
          if(NSAME(I) .gt. 0) then
            NC = NC+1
            GFNR(N,NC) = I
          endif
        enddo
        GNR(N) = NC
        GFNR(N,NC+1) = 0
      enddo
!---number of unique columns in group:  if too many ICAs, drop upper groups!
      NICA = 1
      do N = 1,NRG
        NICA = NICA*GNR(N)
        if (NICA .le. ICA_) then
          NICAX = NICA
          NRGX = N
        endif
      enddo
      if (NICA .gt. ICA_) then
        write(6,*) 'NICA greater than ICA_',NICA,ICA_,NICAX,NRG,NRGX
        NICA = NICAX
        NRG = NRGX
      endif

    1 continue

      END SUBROUTINE ICA_NR



!-----------------------------------------------------------------------
      SUBROUTINE ICA_ALL(CLF,CLT,LTOP,CBINU,ICAU, CFBIN,CLDCOR,NCLDF,  &
           GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA,  WCOL,OCOL)
!-----------------------------------------------------------------------
!    OCOL() = cloud optical depth (total) in each ICA
!    WCOL() = weight(fract area) of ICA,
!    TCOL(,) profile of cloud OD is not calcualted here.
!    ISORT() = index no. of ICA sorted by column OD from smallest to largest
!---Using the information on max-ran cloud overlap generated by ICA_NR,
!---  this usbroutine generates all the ICAs
!   GBOT(I=1:NRG) = lower CTM layer of max-overlap group I
!   GTOP(I=1:NRG) = upper CTM layer of max-overlap group I
!   GNR(I=1:NRG) = no. of unique quantized cloud fractions in group I
!    (.le.NCBINS)
!   GFNR(I=1:NRG,1:GNR(I)) = cloud fraction quantum no (value = 1 to NCBIN)
!---See JL Neu, MJ Prather, JE Penner (2007), Global atmospheric chemistry:
!      Integrating over fractional cloud cover,J. Geophys. Res., 112, D11306,
!       doi:10.1029/2006JD008007

      implicit none

      integer, intent(in) :: LTOP, CBINU, ICAU, NRG, NICA
      integer, intent(in), dimension(LTOP) :: NCLDF
      integer, intent(in), dimension(9)  :: GBOT,GTOP,GLVL,GNR,GCMX
      integer, intent(in), dimension(9,CBINU+1) :: GFNR
      real*8,  intent(in), dimension(LTOP)  :: CLF,CLT
      real*8,  intent(in), dimension(CBINU) :: CFBIN
      real*8,  intent(in)                 :: CLDCOR
      real*8,  intent(out),dimension(ICAU) :: WCOL,OCOL
!      real*8,  intent(out),dimension(LTOP,ICAU) :: TCOL

      real*8  ODCOL,WTCOL,CF0(51),  FWT(10,51),FWTC(10,51),FWTCC(10,51)
      real*8  FIG1,FIG2,GCORR,GCOWT,CORRFAC, FCMX(10) ,CLTOT(100)
      integer I, II, IG1,IG2, G, L,  IGNR(10),GCLDY(10),GRP1,GRP2
      logical L_CLR1,L_CLR2  ,LSKIP   ,LGR_CLR(10)
!-----------------------------------------------------------------------
        CLTOT(:) = 0.d0

        CF0(1) = 0.d0
      do L = 1,CBINU
        CF0(L+1) = CFBIN(L)
      enddo

      do G = 1,NRG
        FCMX(G) = CF0(GCMX(G)+1)           ! max cloud-fraction in MAX-GRP
       if (FCMX(G) .lt. 0.99d0) then
          LGR_CLR(G) = .true.          ! 1st member of MAX-GRP G = clear sky
          GCLDY(G) = 2
       else
          LGR_CLR(G) = .false.
          GCLDY(G) = 1
       endif
       do I = 1,GNR(G)        ! std weighting for each member of each MAX-Group
          FWT(G,I) = CF0(GFNR(G,I)+1) - CF0(GFNR(G,I+1)+1)
          FWTC(G,I) = FWT(G,I)
          FWTCC(G,I) = FWT(G,I)
       enddo
      enddo
        FCMX(NRG+1) = 0.d0

!  pre-calculate correl factors here:  no change if G = 100% cloud or G+1 = 100% cloud or clear
!   also no correlation fix for top MAX-GRP
      do G = 1,NRG-1
       LSKIP =   GCMX(G+1).eq.0 .or. GCMX(G+1).eq.CBINU   &
           .or.  GCMX(G).eq.0 .or. GCMX(G).eq.CBINU     !must have cloudy&clear in both MAX-GRPs
       if (.not.LSKIP) then
          FIG2 = FCMX(G+1)       ! cloudy fract of MAX-GRP just above (G+1)
          GRP2 = GLVL(G+1)       ! upper G6 group for NRG # G+1
          FIG1 = FCMX(G)         ! cloudy fract of current MAX-GRP (sum of cloudy fracts)
          GRP1 = GLVL(G)         ! current G6 group for NRG # G
          CORRFAC = CLDCOR**(GRP2-GRP1)  ! Cloud Correl Factor decreases with gap in G6 groups
        do I = 2,GNR(G)
! correlation factor: increase fract-area of cloudy member under cloudy section of upper layer (FIG2)
! Note that limits to increase depend on fract of cloud area above and the layer being increased

         GCORR = min(1.d0 + CORRFAC*(1.d0/FIG2 - 1.d0), 1.d0/FIG2, 1.d0/FIG1)

! enhance weighting for cloudy members below a cloud, reduce weighting below clear sky
         FWTC(G,I) = GCORR * FWT(G,I)
         FWTCC(G,I) = FWT(G,I) * (1.d0 - GCORR*FIG2)/(1.d0-FIG2)
        enddo
         FWTC(G,1) = 1.d0 - FIG1*GCORR
         FWTCC(G,1) = 1.d0 - FIG1*(1.d0-GCORR*FIG2)/(1.d0-FIG2)
       endif
      enddo

      do I = 1,NICA
          WTCOL = 1.d0
          ODCOL = 0.d0
! for each ICA locate the members of each GROUP that contributes to it
          II = I
        do G = 1,NRG
          IGNR(G) = mod(II-1, GNR(G)) + 1
          II = (II-1)/GNR(G) + 1
        enddo
        do G = 1,NRG
           IG1 = IGNR(G)          ! working on MAX-GRP = G, member IG1
           L_CLR1 = GFNR(G,  IG1) .gt. GCMX(G)      ! member IG1 is clear
         if (G .eq. NRG) then     ! fix of indexing error in Cloud-J 7.3 did not affect results
           L_CLR2 = .true.
         else
           IG2 = IGNR(G+1)        ! member of MAX-GRP = G+1 for this ICA
           L_CLR2 = GFNR(G+1,IG2) .gt. GCMX(G+1)    ! member above (IG2) is clear
         endif
! all of these combinations should preserve the total weighting for layer member IG1
          if (.not.L_CLR2) then
! upper layer GRP member IG2 is a cloud layer (maybe one of several)
            if (.not.L_CLR1) then
! immediate GRP layer member IG1 is a cloudy one
               GCOWT = FWTC(G,IG1)
            else
! immediate layer GRP member IG1 is the clear one (if it exists)
               GCOWT = FWTC(G,1)
            endif
          else
! upper layer GRP member IG2 is a clear layer
            if (.not.L_CLR1) then
! immediate GRP layer member IG1 is a cloudy one
               GCOWT = FWTCC(G,IG1)
            else
! immediate layer GRP member IG1 is the clear one (if it exists)
               GCOWT = FWTCC(G,1)
             endif
          endif
            WTCOL = WTCOL*GCOWT
          do L = GBOT(G),GTOP(G)
            if (NCLDF(L) .ge. GFNR(G,IG1)) then
              ODCOL = ODCOL + CLT(L)
! could store the full 2-D array of atmospheres if needed:  TCOL(L,I) = CLT(L)
            endif
          enddo
        enddo
          WCOL(I) = WTCOL
          OCOL(I) = ODCOL
      enddo

      END SUBROUTINE ICA_ALL



!-----------------------------------------------------------------------
      SUBROUTINE ICA_III(CLF,CLT,LTOP,CBINU,ICAU, III, &
               CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA,  TCOL)
!-----------------------------------------------------------------------
!    see ICA_ALL, this subroutine picks out the ICA atmosphere #III
!      and loads the REFF/WPs for a FAST_JX calculation.

      implicit none

      integer, intent(in) :: LTOP, CBINU, ICAU, NRG, NICA, III
      integer, intent(in), dimension(LTOP) :: NCLDF
      integer, intent(in), dimension(9)  :: GBOT,GTOP,GLVL,GNR,GCMX
      integer, intent(in), dimension(9,CBINU+1) :: GFNR
      real*8,  intent(in), dimension(LTOP)  :: CLF,CLT
      real*8,  intent(in)                 :: CLDCOR
      real*8,  intent(out),dimension(LTOP) :: TCOL

      integer II, IG, G, L

!-----------------------------------------------------------------------

         TCOL(:) = 0.d0
      II = max(1, min(NICA,III))
      do G = 1,NRG
          IG = mod(II-1, GNR(G)) + 1
          II = (II-1)/GNR(G) + 1
        do L = GBOT(G),GTOP(G)
          if (NCLDF(L) .ge. GFNR(G,IG)) then
            TCOL(L) = CLT(L)
          endif
        enddo
      enddo

      END SUBROUTINE ICA_III



!-----------------------------------------------------------------------
      SUBROUTINE ICA_QUD(WCOL,OCOL, LTOP,ICAU,NQDU,NICA, &
                         WTQCA, ISORT,NQ1,NQ2,NDXQS)
!-----------------------------------------------------------------------
!---Take the full set of ICAs and group into the NQD_ ranges of total OD
!---Create the Cumulative Prob Fn and select the mid-point ICA for each group
!---The Quad atmospheres have weights WTQCA
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)        :: LTOP,ICAU,NQDU,NICA
      real*8,  intent(in), dimension(ICAU)      :: WCOL,OCOL

      real*8, intent(out), dimension(NQDU)      :: WTQCA
      integer, intent(out), dimension(ICAU)     :: ISORT
      integer, intent(out), dimension(NQDU)     :: NQ1,NQ2,NDXQS

      real*8,  dimension(ICA_) :: OCDFS, OCOLS
      integer I, II, J, L, N, N1, N2

      real*8, parameter:: OD_QUAD(4) =[0.5d0, 4.0d0, 30.d0, 1.d9]
!-----------------------------------------------------------------------
      ISORT(:) = 0
      WTQCA(:)  = 0.d0
      NDXQS(:) = 0

!---sort all the Indep Column Atmos (ICAs) in order of increasing column OD
!--- ISORT is the key, giving the ICA number from smallest to largest column OD
!--- OCOLS is the column OD sorted = OCOL(ISORT(I))
!--- OCDFS is the Cum.Prob.Fn. of the successive, sorted ICA
      if (NICA .eq. 1)  then
        ISORT(1) = 1
        OCOLS(1) = OCOL(1)
      else
        call HEAPSORT_A (NICA,OCOL,OCOLS,ISORT,ICA_)
      endif
        OCDFS(1) = WCOL(ISORT(1))
      do I = 2,NICA
        OCDFS(I) = OCDFS(I-1) + WCOL(ISORT(I))
      enddo
!---find beginning/end of quad range, note NQ2 < NQ1 means nothing in that range
          I = 1
      do N = 1,NQDU
       do while (OCOLS(I).lt.OD_QUAD(N) .and. I.le.NICA)
          I = I+1
       enddo
        NQ2(N) = I-1
      enddo
        NQ1(1) = 1
      do N = 2,NQDU
        NQ1(N) = NQ2(N-1) + 1
      enddo
!---define QCA wts from cum prob, pick middle ICA as representative
      do N = 1,NQDU
          N1 = NQ1(N)
          N2 = NQ2(N)
       if (N2 .ge. N1) then
          NDXQS(N) = (N1+N2)/2
         if (N1 .gt. 1) then
            WTQCA(N) = OCDFS(N2)-OCDFS(N1-1)
         else
            WTQCA(N) = OCDFS(N2)
         endif
       endif
      enddo

      END SUBROUTINE ICA_QUD



!-----------------------------------------------------------------------
      SUBROUTINE ICA_DIRECT(U0,CLF,CLT,TAULC,TAUIC, LTOP,CBINU,ICAU,  &
                  CLDCOR,NCLDF,GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, &
                            WCOL,TCLD, TTCOL,SOLT,G0L,G0I,TAUG)
!-----------------------------------------------------------------------
      implicit none

      real*8, intent(in)  :: U0, G0L,G0I,   CLDCOR
      integer, intent(in) :: LTOP, ICAU, CBINU, NICA, NRG
      real*8,  intent(in), dimension(LTOP)  :: CLF,CLT, TAULC,TAUIC
      integer, intent(in), dimension(LTOP) :: NCLDF
      integer, intent(in), dimension(9)     :: GBOT,GTOP,GLVL,GNR,GCMX
      integer, intent(in), dimension(9,CBINU+1) :: GFNR
      real*8, intent(in), dimension(ICAU)   :: WCOL
      real*8, intent(out),dimension(LTOP+1)   :: TCLD,TTCOL,SOLT,TAUG
!---Local variables
      real*8  SOLZ, ATTEN, USZA, TAUALL
      integer II, L

!---if low sun angle, just use average cloud cover
      if (U0 .lt. 0.1) then
        do L = 1,LTOP
          TCLD(L) = CLF(L) * CLT(L)
        enddo
      else
        USZA = 1.d0/U0
        SOLT(:) = 0.d0
        TAUG(:) = 1.d0
       do L = 1,LTOP
          TAUALL = TAUIC(L) + TAULC(L)
        if (TAUALL .gt. 1.d-7) then
!  calc g0 to get equivlent isootropic OD below
          TAUG(L) = (G0L*TAULC(L) + G0I*TAUIC(L))/TAUALL
        endif
       enddo
       do II = 1,NICA
         call ICA_III(CLF,CLT,LTOP,CBIN_,ICA_, II, &
            CLDCOR,NCLDF, GFNR,GCMX,GNR,GBOT,GTOP,GLVL,NRG,NICA, TTCOL)
         SOLZ = 1.d0
        do L = LTOP,1,-1
          if (TTCOL(L) .gt. 1.d-9) then
            SOLZ = SOLZ * exp(-USZA*(1.d0-TAUG(L))*TTCOL(L))
          endif
            SOLT(L) = SOLT(L) + WCOL(II)*SOLZ
        enddo
       enddo
!---Derive effective cloud OD (w/asymm factor corr. 1-g) for avg over ICAs
            TCLD(:) = 0.d0
            SOLT(LTOP+1) = 1.d0
        do L = LTOP,1,-1
            ATTEN = SOLT(L+1)/SOLT(L)
          if (ATTEN .gt. 1.00000001d0) then
            TCLD(L) = log(ATTEN)/(USZA*(1.d0-TAUG(L)))
          endif
        enddo
      endif

      END SUBROUTINE ICA_DIRECT



!-----------------------------------------------------------------------
      SUBROUTINE HEAPSORT_A (N,A,AX,IX,ND)
!-----------------------------------------------------------------------
!  classic heapsort, sorts real*8 array A(N) into ASCENDING order,
!     places sorted array AX(N):   AX(1) .le. AX(N)
!     returns indexing IX(N) that records the location of A in sequence:
!           A(IX(J)) ==> AX(J), s.t. IX(1) = orig location of smallest A
!                           and IX(N) = original loc. of largest A
      implicit none
      integer, intent(in)  :: N, ND
      real*8, dimension(ND),intent(in)  :: A
      real*8, dimension(ND),intent(out) :: AX
      integer,dimension(ND),intent(out) :: IX
      integer :: I,J,L,IR,IA
      real*8 :: RA

      do I = 1,N
        IX(I) = I
        AX(I) = A(I)
      enddo
      L  = N/2+1
      IR = N
   10 continue
      if (L .gt. 1) then
        L = L-1
        RA = AX(L)
        IA = IX(L)
      else
        RA = AX(IR)
        IA = IX(IR)
        AX(IR) = AX(1)
        IX(IR) = IX(1)
        IR = IR-1
        if (IR .eq. 1) then
          AX(1) = RA
          IX(1) = IA
          return
        endif
      endif
      I = L
      J = L+L
   20 continue
      if (J .le. IR) then
        if (J .lt. IR) then
          if (AX(J) .lt. AX(J+1)) then
            J = J+1
          endif
        endif
        if (RA .lt. AX(J)) then
          AX(I) = AX(J)
          IX(I) = IX(J)
          I = J
          J = J+J
        else
          J = IR+1
        endif
        goto 20
      endif
        AX(I) = RA
        IX(I) = IA
      goto 10

      END SUBROUTINE HEAPSORT_A

      END MODULE CLD_SUB_MOD
