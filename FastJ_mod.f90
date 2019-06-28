! <FastJ_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
! ---------- JX71_notes.f90

!  same as JX70b except that test for solar flux = 0 so that can run 200 nm data set
!  the 2200 nm data set cuts all J's at exactly 200 nm - used to combien with WACCM tables
!  to allow J's above 60 km, including Lyman-alpha
!
! ----subroutines and calls:
!       main standalone
!   >>>   call INIT_FJX (TITLJXX,NJX_,NJXX)
!         call RD_JS_JX(9,'FJX_j2j.dat', TITLJXX,NJXX)
!         call SOLAR_JX(GMTAU,IDAY,YGRD,XGRD, SZA,U0,SOLF)
!         call ACLIM_FJX (YLAT, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1_)
!         call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
!   >>>   call PHOTO_JX(U0,SZA,REFLB,SOLF, LPRTJ,  PPP,ZZZ,TTT,DDD,RRR,OOO,
!                 LWP,IWP,REFFL,REFFI, AERSP,NDXAER,L1U,ANU,  VALJXX,NJXU)
!
! ----notes:   >>> = only two essential fast-JX calls are denoted above with >>>
!
!   >>> INIT_FJX (TITLJXX,NJX_,NJXX)
!          called once to read in data files, returns thru calling sequence:
!             TITLJXX(1:NJX_) the char*6 title for each fast-JX J-value
!             NJXX the actual number of J-values that fast-JX will compute.
!
!       RD_JS_JX(9,'FJX_j2j.dat', TITLJXX,NJXX)
!           called once after INIT_FJX to map the CTM J's onto the fast-JX J's
!           this is an example, see the data faile 'FJX_j2j.dat'
!
!   >>> PHOTO_JX(U0,SZA,REFLB,SOLF, LPRTJ,  PPP,ZZZ,TTT,DDD,RRR,OOO,
!                  CLDWP,AERSP,NDXCLD,NDXAER,L1_,AN_,    VALJXX,NJX_)
!          called every time for each Indep. Colm. Atmos. (ICA) to compute J's
!          all the information is passed through the calling sequence:
!     ----------------------------------------------------------------------------
!      L1_    CTM has layers 1:L_, for JX must add top-of-atmos layer, L1_ = L_+1
!      AN_    dimension of max number of aerosol types given to JX
!      U0     cosine of the solar zenith angle
!      SZA    solar zenith angle (degrees)
!      REFLB  Lambertian reflectance of lower boundary
!        (NB JX could be adapted to add lower layers for canopy/snow scatt/abs)
!      SOLF   Solar radiation factor, correcting for sun-earth distance
!        (NB JX could be modified for solar cycle, but needs link to internal data)
!      LPRTJ  logical to produce internal printout (stdout) of fast-JX J's, fluxes,
!         absorption, etc.  This information is now only internal to PHOTO_JX
!      PPP(1:L1_+1)   edge press (hPa)
!      ZZZ(1:L1_+1)   edge altitude (cm)
!      TTT(1:L1_)     mid-layer temp (K)
!      DDD(1:L1_)     layer dens-path (# molec /cm2)
!      RRR(1:L1_)     mid-layer relative humidity (0.0 to 1.0)
!      OOO(1:L1_)     layer O3 path (# O3 /cm2)
!      CLDWP(1:L1_)   layer cloud water path (kg/m2), liquid and ice
!      AERSP(1:L1_,1:AN_)  aerosol path (g/m2)
!      NDXCLD(1:L1_)  layer cloud index (type)
!          only a single cloud type is allowed for optical properties, pick
!          the dominant one in terms of optical depth,
!          see notes on cloud types allowed in fast-JX: 'FJX_scatt-cld.dat'
!      NDXAER(1:L1_,1:AN_) aerosol index (type)
!          sample aerosol types are in 'FJX_scat-aer.dat' and 'FJX_scat-UMa.dat'
!          the UMa data allows for relative humidity to be included
!          other aerosol types can be added.
!     ----------------------------------------------------------------------------
!      VALJXX(1:,NJX_,1:L) & NJX_ (first dimension of VALJXX) are returned
!          VALJXX is the array of fast-JX J's, the second dimension is not given
!             but is large enough to accommodate the CTM layers 1:L1_
!          the main code must use the information calcualted by RD_JS_JX to
!             re-map the VALJXX onto the CTM J's.  A useful example is given.
!     ----------------------------------------------------------------------------
!
!       SOLAR_JX calculates solar zenith angle & distance correction (if needed)
!
!       ACLIM_FJX fills in T & O3 from a climatology
!             may be needed for the layer above the CTM to account for O3 & O2
!
!       JP_ATM0 does a simple printout (stdout) of the atmosphere
!
!
! ---------- fjx70sub.f  fast-JX core routines ver 7.0+ (10/2012, mjp)
!
! ----subroutines and calls:  >>> only subroutines called from outside >>>
!      one include 'cmn_FJX.f' is common to several and has parameters, etc.
!      only other connection with CTM code is in call sequence and noted above.
!
!   >>> subroutine INIT_FJX (TITLEJXX,NJXU,NJXX)
!         call RD_XXX(JXUNIT,'FJX_spec.dat')
!         call RD_MIE(JXUNIT,'FJX_scat.dat')
!         call RD_UM (JXUNIT,'FJX_UMaer.dat')
!         call RD_PROF(JXUNIT,'atmos_std.dat')
!       subroutine RD_XXX(NJ1,NAMFIL)
!       subroutine RD_MIE(NJ1,NAMFIL)
!       subroutine RD_UM(NJ1,NAMFIL)
!       subroutine RD_PROF(NJ2,NAMFIL)
!       subroutine EXITC(T_EXIT)
!   >>> subroutine SOLAR_JX(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!   >>> subroutine ACLIM_FJX (YLATD, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1U)
!   >>> subroutine PHOTO_JX(U0,SZA,REFLB,SOLF,LPRTJ, PPP,ZZZ,TTT,DDD,RRR,OOO,
!                          CLDWP,AERSP,NDXCLD,NDXAER,L1U,ANU,  VALJXX,NJXU)
!         call SPHERE2 (U0,RAD,ZZJ,ZZHT,AMF2, L1U,JXL1_)
!         call OPTICL (OPTX,SSAX,SLEGX,  ODCLD,NDCLD)
!         call OPTICA (OPTX,SSAX,SLEGX,  PATH,RH, NAER)
!         call OPTICM (OPTX,SSAX,SLEGX,  PATH,RH,-NAER)
!         call EXTRAL(OD600,L1U,L2U,N_,JTAUMX,ATAU,ATAU0, JXTRA)
!         call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2),.,,)
!         call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2),,,,)
!         call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
!                          FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!         call JRATET(PPJ,TTJ,FFF, VALJXX, LU,NJXU)
!         call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)
!       subroutine OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
!                          FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!         call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!       subroutine MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!         call LEGND0 (EMU(I),PM0,M2_)
!         call LEGND0 (-U0,PM0,M2_)
!         call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)
!       subroutine LEGND0 (X,PL,N)
!       subroutine BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
!         call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K),PM,PM0,
!              B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K),A(1,1,K),H(1,1,K),C(1,1,K),ND)
!       subroutine GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,B,CC,AA,A,H,C,ND)
!       subroutine OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)
!       subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,L)
!       subroutine OPTICM (OPTD,SSALB,SLEG, PATH,RELH,L)
!       subroutine JRATET(PPJ,TTJ,FFF, VALJL,LU,NJXU)
!         call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1),...)
!         call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1),...)
!         call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1),...)
!         call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J),...)
!         call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J),...)
!       subroutine X_interp (TINT,XINT, T1,X1, T2,X2, T3,X3, L123)
!       subroutine JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU6,POMEG6,JXTRA,LU)
!   >>> subroutine JP_ATM0(PPJ,TTJ,DDJ,OOJ,ZZJ, LU)
!       subroutine SPHERE2(U0,RAD,ZHL,ZZHT,AMF2, L1U,LJX1U)
!       subroutine EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
! note - the calls to EXITC are not listed here
!       e.g.,  call EXITC(' INIT_JX: invalid no. wavelengths')
!
!
!
!  >>>>>>>>>>>>>>>>current code revised to JX ver 7.0+ (10/12)<<<<<<<<<<<<
!
!  fastJX version 7.0+ (f90) - Prather notes (Jan 2013)
!
!---calculation of cloud optical depth in FJX-70b !!!
!---    assumes that clouds are 100% if in layer
!
!   IWP = ice water path (in layer, in cloud) in g/m**2
!   LWP = liquid water path (in layer, in cloud) in g/m**2
!   REFFI = effective radius of ice cloud (microns)
!   REFFL = effective radius of liquid cloud (microns)
!
!>>>>method for calculating cloud OD (in layer) used by FJX core or outside
!>>>>FJX core needs only the _WP and the REFF_
!>>>> note that FJX can use correct Q's and density updates if need be.
!   ODL = LWP * 0.75 * 2.10 / REFFL
!   ODI = IWP * 0.75 * 2.00 / (REFFI * 0.917)
!
!>>>R-effective determined by main code, not FJX
!   REFF determined by user - some recommendations below (from Neu & Prather)
!       REFFI is a simple function of ice water content IWC (g/m3, 0.0001 to 0.1)
!          IWC = IWP / delta-Z (of layer in m, approx OK)
!          REFFI = 50. * (1. + 8.333 * IWC)
!   prefer Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)
!
!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR))
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)
!
!>>>indices for cloud scattering determined by FJX core, not main code.
!   NDX = cloud scattering index based on REFFL or TCLD(for ice cloud)
!          NDXI = 6  ! ice hexag (cold)
!             if (TCLD .ge. 233.15) then
!          NDXI = 7  ! ice irreg
!             end if
!          NDXC = 1
!             do I=2,4
!             if (REFFL .gt. 0.5*(RCC(I-1)+RCC(I))) then
!          NDXC = I
!             end if
!             end do
!
!
! version 7.0
!       New modular structure:
!           Designed for .f90 and CAM5
!           All set up routines grouped into one super-call
!           Separate interface with CTM and ICAs calls PHOT_J
!           PHOT_JX only sees a single column atmosphere (or ICA) thru calling params
!
! version 6.8
!       New layout and formatting of FJX_spec.dat (required)
!              allows 1,2, or 3 T's (or P's) for all X-sects
!       VOCs mostly use Pressure interp.
!
! version 6.7  (should output same J-values as ver 6.6 EXCEPT for VOCs)   (3/12)
!       New JPL-2010 update complete through all VOCs (see notes in FJX_spec.dat)
!           Most important change is incorporation of Stern-Volmer pressure-dependent
!           and wavelength-dependent quantum yields into the Temperature interpolation.
!           Acetone now has only one pair of X sections for each pathway.
!       Redo of mapping of (external) chemical model J's onto fastJX J's
!           see examples of splitting a J with fixed branching q-yields
!           see how chem model labels are not tied by the cross section labels in fastJX
!       Changes in FX_spec.dat make it incompatible with earlier versions, although
!           all of the Xsection data have the same format.
!       Now the number of X sections read in equals the number of J's that fast-JX calculates.
!       As before, all J's except for O2, O3, O3(1D) have only pairs of data at different T.
!
! version 6.6x
!       N.B. SPECIAL FIX of tropospheric bin 11 used in FULL(W_=18) and TROP-ONLY(W_=12)
!            Attenuation of this bin is too weak since it mixes 218 nm and 288 nm
!            Thus it allows a low level of J-O2 (and related 218 nm Xsects) thru-out trop.
!            The attenuation in strat (1e-4) is OK, but must zero out in trop (>100 hPa)
!            The TROP-QUICK (W_=8) does not use bin#11 and is unaffected.
!       Redo of the 4x4 matrix inversions from subroutine to inline code (faster)
!       FJX_spec.dat:  includes new JPL2010 and solar flux, but J-VOC unchanged (a to-do).
!           The J-O2 and O3 steady-state are improved, but the NOy:N2O is too small.
!           The J-NO appears too large based on both photocomp2008 & NOy values.
!             Possibly lack of thermopsheric NO absorption, or correlation with S-R bands.
!           ACTION: J-NO X-sections multiplied by 0.6 in new FJX_spec.dat
!
! version 6.5
!      J-values are now averaged over CTM layer (and not just mid-layer point)
!        This is important when cloud optical depth is large (~1).
!
! version 6.4
!     allows for shortened, speeded up troposphere versions<<<<<<
!     STD:           W_=18
!        identical results to v-6.2 if cloud OD is consistent
!     TROP-ONLY:     W_=12
!        collapses the wavelength bins from 18 to 12 (5-6-7-8 & 11-18)
!        drops many 'stratospheric' cross-sections (denoted by 'x' in 2nd title)
!        allows use of single standard spectral data set:  FJX_spec.dat
!        results close to W_=18, largest difference is J-O2 (<1% in 13-18 km!!)
!        This is recommended as accurate for troposphere only calculations.
!     TROP-QUICK:    W_=8
!        reverts to original fast-J 7-bins (12-18) plus 1 scaled UV (5) for J-O2
!        errors in 12-18 km range for J-O2, high sun are 10%, worse above.
!     ***Photolysis of O2 in the upper tropical troposphere is an important
!        source of O3.  It needs to be included in tropospheric runs.
!        TROP-ONLY is recommended, W_=8 is a quick fix if speed essential.
!
!     Major rewrite of code to minimize calls and allow better vector-type ops.
!     loop over wavelengths internal to Mie soln.
!     Driven by profiling of CTM code, may still need optimization.
!     Wavelengths can be contracted to W_=12 (trop only) and strat-only
!        X-sections are dropped.  With parm W_=18, the std fast-JX is retrieved.
!     Many call eliminated and summations in BLKSLV and GEN_ID are explicit
!     GEN_ID replaces GEN and calculates all matrix coeff's (1:L_) at once
!     RD_XXX changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 (parm_CTM.f) should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (parm_MIE.f).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!
! version 6.3
!     revise cloud/aerosol OD & wavelength properties for CTM link:
!         OPTICL is new sub for cloud optical properties, but it
!              now starts with cloud OD @ 600 nm and cloud NDX
!              fast-JX now uses cloud NDX to scale OD to other wavelengths
!         OPTICA & OPTICM are new subs to convert aerosol path (g/m2) to OD
!              A is std UCI scat data
!              M is U Michigan data tables for aerosols, includes Rel Hum effect
!     drop sub GAUSSP and put into Parameter statement (parm_MIE.f)
!
! version 6.2
!     corrects a long-standing problem at SZA > 89 degrees.
!     In prior versions the ray-tracing of the path (and air-mass functions)
!     back to the sun was done at the edges of the CTM layers (it was developed
!     for the grid-point J-value code at Harvard/GISS/UCI).  This left the
!     interpolation to the mid-layer (needed for J's) open.  The prior method
!     gave irregular fluctuations in the direct solar beam at mid-layer for
!     large SZA > 88.  This is now corrected with exact ray-tracing from
!     the mid-pt of each CTM layer.  For small SZA, there is no effective
!     difference, for large SZA, results could be erratic.
!   v-6.2 fix should be easy if you have migrated to v6.1, else some minor
!      caution may be needed:
!      replace sub SPHERE with SPHERE2, AMF2 report factors for mid and egdes.
!      replace sub OPMIE with new OPMIE, this uses the new AMF2 correctly.
!      replace sub PHOTOJ with new PHOTOJ, this just hands off AMF2 from
!            SPHERE2 to OPMIE.
!
! version 6.1 adds
!      6.1b simplifies calling sequences feeds solar factor, albedo, to PHOTOJ
!         and read LAT, LNG directly.  No substantive changes.
!      new read-in of scat data for clouds/aerosols to allow for UMich data
!      This has required substantial rewrite of some of the core subroutines:
!         OPMIE is now called for each wavelength and without aersol/cloud data
!              all subs below OPMIE are unchanged
!         OPTICD & OPTICM are new subs to convert path (g/m2) to OD and phase fn
!              D is std UCI scat data (re-ordered for clouds 1st)
!              M is U Michigan data tables for aerosols, includes Rel Hum effect
!         PHOTOJ now assembles the aerosol data (better for CTM implementation)
!      This version can reproduce earlier versions exactly, but the test input
!         is changed from OD and NDX to PATH (g/m2) and NDX.
!
! version 6.0 adds
!      new 200-nm scattering data so that stratospheric aerosols can be done!
! version 5.7
!>>>>>>>>>>>>>>>>current code revised to JX ver 7.1 (10/13)<<<<<<<<<<<<

!------------------------------------------------------------------------------
!    'FJX_CMN_MOD.f90'  for fast-JX code v 7.0+ (prather 9/12)                  !
!------------------------------------------------------------------------------
!
! NB - ALL of these common variables are set paramters,
!    They are NOT to be used as variables for a local solution
!    Thus this entire set is 'in' only after it is initialized
!-----------------------------------------------------------------------
!
! !DESCRIPTION: FJX_CMN contains fast-JX variables
!
!
! !INTERFACE:
!
      MODULE FJX_CMN_MOD

      IMPLICIT NONE
      PUBLIC

!-----------------------------------------------------------------------

      ! JXL_: vertical(levels) dim for J-values computed within fast-JX
      INTEGER, PARAMETER ::  JXL_=100, JXL1_=JXL_+1
      ! JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid (mid-level)
      INTEGER, PARAMETER ::  JXL2_=2*JXL_+2
      ! WX_  = dim = no. of wavelengths in input file
      INTEGER, PARAMETER ::  WX_=18
      ! X_   = dim = max no. of X-section data sets (input data)
      INTEGER, PARAMETER ::  X_=72
      ! A_   = dim = no. of Aerosol/cloud Mie sets (input data)
      INTEGER, PARAMETER ::  A_=40
      ! C_   = dim = no. of cld-data sets (input data)
      INTEGER, PARAMETER ::  C_=16
      ! W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
      INTEGER, PARAMETER ::  W_=18    ! W_=12 or 18
      ! N_  = no. of levels in Mie scattering arrays
      !     = 2*NC+1 = 4*(L_+1) + 1`+ 2*sum(JADDLV)
      INTEGER, PARAMETER ::  N_=601
      ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
      INTEGER, PARAMETER ::  M_=4
      ! M2_ = 2*M_ = 8, replaces MFIT
      INTEGER, PARAMETER ::  M2_=2*M_
!-----------------------------------------------------------------------
      ! 4 Gauss pts = 8-stream
      REAL*8, DIMENSION(M_), PARAMETER  ::                            &
                        EMU = [.06943184420297d0, .33000947820757d0,  &
                               .66999052179243d0, .93056815579703d0]
      REAL*8, DIMENSION(M_), PARAMETER  ::                            &
                        WT  = [.17392742256873d0, .32607257743127d0,  &
                               .32607257743127d0, .17392742256873d0]
!-----------------------------------------------------------------------

      ! ZZHT: scale height (cm)
      REAL*8, PARAMETER   :: ZZHT = 5.d5
      ! RAD: Radius of Earth (cm)
      REAL*8, PARAMETER   :: RAD = 6375.d5
      ! ATAU: heating rate (factor increase from one layer to the next)
      REAL*8, PARAMETER   :: ATAU = 1.120d0
      ! ATAU0: minimum heating rate
      REAL*8, PARAMETER   :: ATAU0 = 0.010d0
      ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
      INTEGER, PARAMETER  :: JTAUMX = (N_ - 4*JXL_)/2

!---- Variables in file 'FJX_spec.dat' (RD_XXX)

      ! WBIN: Boundaries of wavelength bins
      REAL*8  WBIN(WX_+1)
      ! WL: Centres of wavelength bins - 'effective wavelength'
      REAL*8  WL(WX_)
      ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      REAL*8  FL(WX_)

      REAL*8  QO2(WX_,3)   ! QO2: O2 cross-sections
      REAL*8  QO3(WX_,3)   ! QO3: O3 cross-sections
      REAL*8  Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield

      ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      REAL*8  QQQ(WX_,3,X_)
      ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      REAL*8  QRAYL(WX_+1)
      ! TQQ: Temperature for supplied cross sections
      REAL*8  TQQ(3,X_)
      ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      INTEGER LQQ(X_)

      ! TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*6  TITLEJX(X_)
      ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*1  SQQ(X_)

!---- Variables in file 'FJX_scat-aer.dat' (RD_MIE)

      ! QAA: Aerosol scattering phase functions
      REAL*8  QAA(5,A_)
      ! WAA: 5 Wavelengths for the supplied phase functions
      REAL*8  WAA(5,A_)
      ! PAA: Phase function: first 8 terms of expansion
      REAL*8  PAA(8,5,A_)
      ! RAA: Effective radius associated with aerosol type
      REAL*8  RAA(A_)
      ! SAA: Single scattering albedo
      REAL*8  SAA(5,A_)
      ! DAA: density (g/cm^3)
      REAL*8  DAA(A_)
      ! NAA: Number of categories for scattering phase functions
      INTEGER NAA

!---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)

      ! QCC: Cloud scattering phase functions
      REAL*8  QCC(5,C_)
      ! WCC: 5 Wavelengths for supplied phase functions
      REAL*8  WCC(5,C_)
      ! PCC: Phase function: first 8 terms of expansion
      REAL*8  PCC(8,5,C_)
      ! RCC: Effective radius associated with cloud type
      REAL*8  RCC(C_)
      ! SCC: Single scattering albedo
      REAL*8  SCC(5,C_)
      ! DCC: density (g/cm^3)
      REAL*8  DCC(C_)
      ! NCC: Number of categories for cloud scattering phase functions
      INTEGER NCC

!---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)

      ! WMM: U Michigan aerosol wavelengths
      REAL*8  WMM(6)
      ! UMAER: U Michigan aerosol data sets
      REAL*8  UMAER(3,6,21,33)

!---- Variables in file 'atmos_std.dat' (RD_PROF)

      ! T and O3 reference profiles
      REAL*8, DIMENSION(51,18,12) :: TREF, OREF


      INTEGER NJX,NW1,NW2

      END MODULE FJX_CMN_MOD
!------------------------------------------------------------------------------
!     'cmn_CTM.f'  for fast-JX code v7.1                                     !
!------------------------------------------------------------------------------
!
! !MODULE: FJX
!
! !DESCRIPTION: JX version 7.1  (10/13)  works with FL=0 for <200 nm
!
!
! !INTERFACE:
!
      MODULE FJX_SUB_MOD
!
! !USES:
!
      USE FJX_CMN_MOD

      IMPLICIT NONE
      PUBLIC

!-----------------------------------------------------------------------
      integer, parameter :: &
              L_CTM=57, L1_CTM=L_CTM+1 &   ! L_ = number of CTM layers
              ,L_=67, L1_=L_+1 &   ! L_ = max number of layers EMEP+stratosphere
!The results should be independent of the values of L_
!             ,L2_=2*L_+2 &        ! no. levels in the Fast-JX grid that
                                ! includes both layer edges and layer mid-points
             ,JVL_=57 &   ! vertical(levels) dim for J-values sent to CTM
             ,JVN_=101 &  ! max no. of J-values
             ,AN_=25      ! no of separate aerosols per layer (needs NDX for each)

!-----------------------------------------------------------------------
      ! variables used to map fast-JX J's onto CTM J's
!-----------------------------------------------------------------------
      integer::debug
      real*8  JFACTA(JVN_)  ! multiplication factor for fast-JX calculated J
      integer JIND(JVN_)    ! index arrays that map Jvalue(j) onto rates
      integer NRATJ         ! number of Photolysis reactions in CTM chemistry,
                            ! derived here NRATJ must be .le. JVN_
      character*50 JLABEL(JVN_) ! label of J-value used in the main chem model

!
! !PUBLIC SUBROUTINES:
!
      PUBLIC  :: EXITC, SOLAR_JX, ACLIM_FJX, JP_ATM0, PHOTO_JX


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
      real*8  PI, PI180, LOCT
      real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
!
      PI     = 3.141592653589793d0
      PI180  = PI/180.d0
      SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*PI180)
      SOLDEK = asin(SINDEC)
      COSDEC = cos(SOLDEK)
      SINLAT = sin(YGRDJ)
      SOLLAT = asin(SINLAT)
      COSLAT = cos(SOLLAT)
!
      LOCT   = (((GMTIME)*15.d0)-180.d0)*PI180 + XGRDI
      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      SZA    = acos(COSSZA)/PI180
!
!      write(6,*) ' XGRDI,YGRDJ',XGRDI,YGRDJ
!      write(6,*) ' LOCT (rad)',LOCT
!      write(6,*) ' SINDEC,COSDEC', SINDEC,COSDEC
!      write(6,*) ' SINLAT,COSLAT', SINLAT,COSLAT
!      write(6,*) ' COS, SZA',COSSZA,SZA
!
      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*2.d0*PI/365.d0))
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

      real*8, dimension(51) :: OREF2,TREF2
      real*8, dimension(52) :: PSTD
      integer  K, L, M, N
      real*8   DLOGP,F0,T0,PB,PC,XC,MASFAC,SCALEH

!  Select appropriate month
      M = max(1,min(12,MONTH))
!  Select appropriate latitudinal profiles
      N = max(1, min(18, (int(YLATD+99)/10 )))
      do K = 1,51
        OREF2(K) = OREF(K,N,M)
        TREF2(K) = TREF(K,N,M)
      end do

!  Apportion O3 and T on supplied climatology z levels onto CTM levels +1
!  with mass (pressure) weighting, assuming constant mixing ratio and
!  temperature half a layer on either side of the point supplied.
!   PPP(L=1:L1_)=edge-pressure of CTM layer, PPP(L1_+1)=0 (top-of-atmos)
!  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

!  Set up pressure levels for O3/T climatology - assume that value
!  given for each 2 km z* level applies from 1 km below to 1 km above,
!  so select pressures at these boundaries. Surface level values at
!  1000 mb are assumed to extend down to the actual PSURF (if > 1000)
           PSTD(1) = max(PPP(1),1000.d0)
           PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0)
           DLOGP   = 10.d0**(-2.d0/16.d0)
      do K = 3,51
        PSTD(K) = PSTD(K-1)*DLOGP
      end do
        PSTD(52)  = 0.d0
      do L = 1,L1U
        F0 = 0.d0
        T0 = 0.d0
        do K = 1,51
          PC   = min(PPP(L),PSTD(K))
          PB   = max(PPP(L+1),PSTD(K+1))
          if (PC .gt. PB) then
            XC = (PC-PB)/(PPP(L)-PPP(L+1))
            F0 = F0 + OREF2(K)*XC
            T0 = T0 + TREF2(K)*XC
          end if
        end do
        TTT(L)  = T0
        DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
        OOO(L) = F0*1.d-6*DDD(L)
      end do

!  Calculate effective altitudes using scale height at each level
        ZZZ(1) = 0.d0
      do L = 1,L1U-1
        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
        ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH )
      end do
        ZZZ(L1U+1) = ZZZ(L1U) + ZZHT

      END SUBROUTINE ACLIM_FJX

!<<<<<<<<<<<<<<<<<<<end CTM-fastJX special call subroutines<<<<<<<<<<<<<




!<<<<<<<<<<<<<<<<<<<<<<<<begin fastJX core code<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<all outside calls go through PHOTO_JX<<<<<<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
!  fastJX version 7.2 (f90) - Prather notes (Jan 2013)

!---calculation of cloud optical depth in FJX-72 !!!
!---    assumes that clouds are 100% if in layer

!   IWP = ice water path (in layer, in cloud) in g/m**2
!   LWP = liquid water path (in layer, in cloud) in g/m**2
!   REFFI = effective radius of ice cloud (microns)
!   REFFL = effective radius of liquid cloud (microns)

!>>>>method for calculating cloud OD (in layer) used by FJX core or outside
!>>>>FJX core needs only the _WP and the REFF_
!>>>> note that FJX can use correct Q's and density updates if need be.
!   ODL = LWP * 0.75 * 2.10 / REFFL
!   ODI = IWP * 0.75 * 2.00 / (REFFI * 0.917)

!>>>R-effective determined by main code, not FJX
!   REFF determined by user - some recommendations below (from Neu & Prather)
!       REFFI is a simplle function of ice water content IWC (g/m3, 0.0001 to 0.1)
!          IWC = IWP / delta-Z (of layer in m, approx OK)
!          REFFI = 50. * (1. + 8.333 * IWC)
!   prefer Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)

!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR))
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)

!>>>indices for cloud scattering determined by FJX core, not main code.
!   NDX = cloud scattering index based on REFFL or TCLD(for ice cloud)
!          NDXI = 6  ! ice hexag (cold)
!             if (TCLD .ge. 233.15) then
!          NDXI = 7  ! ice irreg
!             end if
!          NDXC = 1
!             do I=2,4
!             if (REFFL .gt. 0.5*(RCC(I-1)+RCC(I))) then
!          NDXC = I
!             end if
!             end do




!>>>>>revise PHOTO_JX to drop NDXCLD, and to include
!   LWP, REFFL, IWP, REFFI only, use T, etc to get NDX.

!-----------------------------------------------------------------------
      subroutine PHOTO_JX                                             &
           (U0,SZA,REFLB,SOLF, LPRTJ,PPP,ZZZ,TTT,DDD,RRR,OOO,            &
           LWP,IWP,REFFL,REFFI, AERSP,NDXAER,L1U,L1_SIZE,ANU,  VALJXX,NJXU, FJBOT,FSBOT,debug)
        !-----------------------------------------------------------------------
        !
        !  PHOTO_JX is the gateway to fast-JX calculations:
        !    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
        !    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
        !    needs day-fo-year for sun distance, SZA (not lat or long)
        !-----------------------------------------------------------------------
        implicit none

        !---calling sequence variables
        integer, intent(in)                   :: L1U,L1_SIZE,ANU,NJXU,debug
        real*8,  intent(in)                   :: U0,SZA,REFLB,SOLF
        logical, intent(in)                   :: LPRTJ
        real*8,  intent(in),dimension(L1_SIZE+1)  :: PPP,ZZZ
        real*8,  intent(in),dimension(L1_SIZE  )  :: TTT,DDD,RRR,OOO
        real*8,  intent(in),dimension(L1_SIZE  )  :: LWP,IWP,REFFL,REFFI
        real*8,  intent(in),dimension(L1_SIZE,ANU):: AERSP
        integer, intent(in),dimension(L1_SIZE,ANU):: NDXAER
        !---reports out the JX J-values, upper level program converts to CTM chemistry J's
        real*8, intent(out), dimension(L1_SIZE-1,NJXU)::  VALJXX
        real*8, intent(out), dimension(W_)     :: FJBOT,FSBOT

        !-----------------------------------------------------------------------
        !--------key LOCAL atmospheric data needed to solve plane-parallel J----
        !-----these are dimensioned JXL_, and must have JXL_ .ge. L_
        real*8, dimension(JXL1_+1) :: TTJ,DDJ,OOJ,ZZJ
        real*8, dimension(JXL1_+1) :: PPJ,RELH
        integer,dimension(JXL2_+1) :: JXTRA
        real*8, dimension(W_)     :: FJTOP,FLXD0,RFL
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


        if (L1U .gt. JXL1_) then
           call EXITC(' PHOTO_JX: not enough levels in JX')
        end if

        LU = L1U - 1
        L2U = LU + LU + 2

        FFF(:,:) = 0.d0
        FREFI = 0.d0
        FREFL = 0.d0
        FREFS = 0.d0

        !---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
        !                        or         99.0                 80 km
        if (SZA .gt. 98.d0)then
           VALJXX=0.0
           goto 99
        end if
        !---load the amtospheric column data
        do L = 1,L1U
           PPJ(L) = PPP(L)
           TTJ(L) = TTT(L)
           DDJ(L) = DDD(L)
           OOJ(L) = OOO(L)
        end do
        PPJ(L1U+1) = 0.d0

        !---calculate spherical weighting functions (AMF: Air Mass Factor)
        do L = 1,L1U+1
           ZZJ(L) = ZZZ(L)
        end do

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
              end do
           end do
        end do

        do L = 1,L1U

           !---Liquid Water Cloud
           if (LWP(L) .gt. 1.d-5 .and. REFFL(L) .gt. 0.1d0) then
              !---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
              ODL = LWP(L) * 0.75d0 * 2.1d0 / REFFL(L)
              NDXL = 1
              do I=2,4
                 if (REFFL(L) .gt. 0.5*(RCC(I-1)+RCC(I))) then
                    NDXL = I
                 end if
              end do
              call OPTICL (OPTX,SSAX,SLEGX,  ODL,NDXL)
              do K=1,5
                 OD(K,L)  = OD(K,L)  + OPTX(K)
                 SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
                 do I=1,8
                    SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
                 end do
              end do
              !>>>diagnostic print of cloud data:
              !        write(6,'(a,i3,2f8.2,f8.4,f8.2,f8.4,i4)') &
              !         'Liq Cld',L,PPP(L),PPP(L+1),LWP(L),REFFL(L),ODL,NDXL
           end if

           !---Ice Water Cloud
           if (IWP(L) .gt. 1.d-5 .and. REFFI(L) .gt. 0.1d0) then
              ODI = IWP(L) * 0.75d0 * 2.0d0 / (REFFI(L) * 0.917d0)
              if (TTT(L) .ge. 233.15d0) then
                 NDXI = 7  ! ice irreg
              else
                 NDXI = 6  ! ice hexag (cold)
              end if
              call OPTICL (OPTX,SSAX,SLEGX,  ODI,NDXI)
              do K=1,5
                 OD(K,L)  = OD(K,L)  + OPTX(K)
                 SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
                 do I=1,8
                    SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
                 end do
              end do
              !>>>diagnostic print of cloud data:
              !        write(6,'(a,i3,2f8.2,f8.4,f8.2,f8.4,i4)') &
              !         'Ice Cld',L,PPP(L),PPP(L+1),IWP(L),REFFI(L),ODI,NDXI
           end if

           !---aerosols in layer: check aerosol index
           !---this uses data from climatology OR from current CTM (STT of aerosols)

           !---FIND useful way to sum over different aerosol types!
           do M = 1,ANU
              NAER = NDXAER(L,M)
              PATH = AERSP(L,M)
              RH = RRR(L)
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
                 end if
                 do K=1,5
                    OD(K,L)  = OD(K,L)  + OPTX(K)
                    SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
                    do I=1,8
                       SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
                    end do
                 end do
              end if
           end do

           do K=1,5
              if (OD(K,L) .gt. 0.d0) then
                 SSA(K,L) = SSA(K,L)/OD(K,L)
                 do I=1,8
                    SLEG(I,K,L) = SLEG(I,K,L)/OD(K,L)
                 end do
              end if
           end do

           !---Include aerosol with cloud OD at 600 nm to determine added layers:
           OD600(L) = OD(4,L)

        end do

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
                 call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2), &
                      TQQ(3,1),QO2(K,3), LQQ(1))
                 call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2), &
                      TQQ(3,2),QO3(K,3), LQQ(2))
                 ODABS = XQO3*OOJ(L) + XQO2*DDJ(L)*0.20948d0
                 ODRAY = DDJ(L)*QRAYL(K)

                 DTAUX(L,K) = OD(KMIE,L) + ODABS + ODRAY

                 do I=1,8
                    POMEGAX(I,L,K) = SLEG(I,KMIE,L)*OD(KMIE,L)
                 end do
                 POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0d0*ODRAY
                 POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5d0*ODRAY
                 do I=1,8
                    POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K)
                 end do
              end do

           end if
        end do

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

              !----direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
              !----     also at bottom (DN), does not include diffuse reflected flux.
              FLXUP(K) =  FJTOP(K)
              DIRUP(K) = -FLXD0(K)
              FLXDN(K) = -FJBOT(K)
              DIRDN(K) = -FSBOT(K)

              do L = 1,LU
                 FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L,K)
              end do
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
              do L=2,LU
                 FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K)
              end do
              FLXJ(LU+1) = FJTOP(K) - FJFLX(LU,K)
              !---calculate net flux deposited in each CTM layer (direct & diffuse):
              FFX0 = 0.d0
              do L=1,L1U
                 FFX(K,L) = FLXD(L,K) - FLXJ(L)
                 FFX0 = FFX0 + FFX(K,L)
              end do

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
           end if
        end do       ! end loop over wavelength K
        !-----------------------------------------------------------------------
        FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
        FREFI = FREFI/FREFS

        !---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)

        !-----------------------------------------------------------------------
        call JRATET(PPJ,TTJ,FFF, VALJXX, LU,L1_SIZE-1,NJXU)
        !-----------------------------------------------------------------------

        !---mapping J-values from fast-JX onto CTM chemistry is done in main

        !-----------------------------------------------------------------------
        if (LPRTJ) then
           !---diagnostics below are NOT returned to the CTM code
           !      write(6,*)'fast-JX-(7.0)---PHOTO_JX internal print: Atmosphere---'
           !---used last called values of DTAUX and POMEGAX, should be 600 nm
           do L=1,L1U
              DTAU600(L) = DTAUX(L,W_)
              do I=1,8
                 POMG600(I,L) = POMEGAX(I,L,W_)
              end do
           end do

           call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)

           if(.false.)then
              !---PRINT SUMMARY of mean intensity, flux, heating rates:
              write(6,*)
              write(6,*)'fast-JX(7.0)---PHOTO_JX internal print: Mean Intens---'
              write(6,'(a,5f10.4)') &
                   ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/', &
                   RFLECT,SZA,U0,FREFI,FREFL

              write(6,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
              write(6,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
              write(6,'(a)') ' ----  100000=Fsolar   MEAN INTENSITY per wvl bin'
              RATIO(:) = 0.d0
              do L = LU,1,-1
                 do K=NW1,NW2
                    if (FL(K) .gt. 1.d0) then
                       RATIO(K) = (1.d5*FFF(K,L)/FL(K))
                    end if
                 end do
                 write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
              end do

              write(6,*)
              write(6,*)'fast-JX(7.0)---PHOTO_JX internal print: Net Fluxes---'
              write(6,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
              write(6,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
              !      write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
              !      write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
              write(6,*) ' ---NET FLUXES--- '
              write(6,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3), K=NW2,NW1,-1)
              write(6,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4), K=NW2,NW1,-1)
              write(6,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5), K=NW2,NW1,-1)
              write(6,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6), K=NW2,NW1,-1)
              write(6,*) ' ---SRF FLUXES--- '
              write(6,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7), K=NW2,NW1,-1)
              write(6,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8), K=NW2,NW1,-1)
              write(6,'(2a)') '  ---NET ABS per layer:       10000=Fsolar', &
                   '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
              do L = LU,1,-1
                 do K=NW1,NW2
                    RATIO(K) = 1.d5*FFX(K,L)
                 end do
                 write(6,'(i9,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
              end do
              write(6,'(a)')
              write(6,'(a)') ' fast-JX (7.0)----J-values----'
              write(6,'(1x,a,72(a6,3x))') 'L=  ',(TITLEJX(K), K=1,NJX)
              do L = LU,1,-1
                 write(6,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
              end do
           end if

        end if

99      continue



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
      end if

      L1U = LU + 1
      L2U = 2*LU + 2


      do L2 = 1,L2U,1
        JADDLV(L2) = JXTRA(L2)
      end do
        JADDTO(L2U+1) = 0
      do L2 = L2U,1,-1
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
      end do

!---expanded grid now included CTM edge and mid layers plus expanded
!---    grid to allow for finer delta-tau at tops of clouds.
!---    DIM of new grid = L2U + JADDTO(1) + 1

!---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
!     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2U+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
      end do

!---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
!---JNELEV(L=1:L_) = L2-index for top of layer L
      do L = 1,LU
        JNDLEV(L) = L2LEV(2*L)
        JNELEV(L) = L2LEV(2*L+1)
      end do
        JNELEV(LU+1) = 0  !need to set this to top-of-atmosphere

      ND = 2*L2U + 2*JADDTO(1) + 1

      if(ND .gt. N_) then
        call EXITC (' overflow of scatter arrays: ND > N_')
      end if

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
       end do
        DTAU(L1U+1,K) = 0.d0

!---Define the total scattering phase fn for each CTM layer L=1:L_+1
!---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
!---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
       do L = 1,L1U
        do I = 1,M2_
          POMEGAJ(I,L,K) = POMEGAX(I,L,K)
        end do
       end do

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
         end do
         if (XLTAU .lt. 76.d0) then   ! zero out flux at 1e-33
          FTAU2(LL,K) = exp(-XLTAU)
         end if
        end if
       end do

!---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
!---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
          FLXD2(:,:) = 0.d0
       do LL = 1,2*L1U
        if (AMF2(LL,LL) .gt. 0.d0) then
          FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL)
        end if
       end do
        if (AMF2(1,1) .gt. 0.d0) then
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1)
        else
          FSBOT(K) = 0.d0
        end if

       do LL = 2,2*L1U,2
         L=LL/2
         FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K)
       end do

!---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
!---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
        FLXD0(K) = 0.d0
       if (AMF2(2*L1U,2*L1U) .gt. 0.d0) then
        do L=1,L1U
         FLXD0(K) = FLXD0(K) + FLXD(L,K)
        end do
       end if

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
       end do

!----solar flux incident on lower boundary & Lambertian reflect factor:
       if (FSBOT(K) .gt. 0.d0) then
        ZFLUX(K) = FSBOT(K)*RFL(K)/(1.d0+RFL(K))
       else
        ZFLUX(K) = 0.d0
       end if

!  Calculate scattering properties, level centres then level boundaries
!>>>>>be careful of order, we are overwriting/shifting the 'POMEGAJ' upward in index
       do L2 = L2U,2,-2
        L   = L2/2
        do I = 1,M2_
          POMEGAJ(I,L2,K) = POMEGAJ(I,L,K)
        end do
       end do
!---lower boundary value is set (POMEGAJ(I,1)), but set upper:
       do I = 1,M2_
         POMEGAJ(I,L2U+1,K) = POMEGAJ(I,L2U,K)
       end do
!---now have POMEGAJ filled at even points from L2=3:L2_-1
!---use inverse interpolation for correct tau-weighted values at edges
       do L2 = 3,L2U-1,2
        TAUDN = TTAU(L2-1,K)-TTAU(L2,K)
        TAUUP = TTAU(L2,K)-TTAU(L2+1,K)
        do I = 1,M2_
          POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN + &
                 POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
        end do
       end do

!---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
!---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface

       do L2 = 1,L2U+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = ND + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K)
          FZ(LZ,K)   = FTAU2(L2,K)
        do I=1,M2_
          POMEGA(I,LZ,K) = POMEGAJ(I,L2,K)
        end do
       end do

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
         end do

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
            if (FBTM(K) .lt. 1.d-32) then
              FZ(LZZ,K) = 0.d0
            else
              FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA
            end if
          end if
          do I = 1,M2_
            POMEGA(I,LZZ,K) = POMEGAB(I,K) + &
                     ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))
          end do
            TAUBTM(K)    = ZTAU(LZZ,K)
            FBTM(K)      = FZ(LZZ,K)
          do I = 1,M2_
            POMEGAB(I,K) = POMEGA(I,LZZ,K)
          end do
         end do
        end if
       end do

!   Now fill in the even points with simple interpolation in scatter arrays:
       do LZ = 2,ND-1,2
         ZTAU(LZ,K) = 0.5d0*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
         FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
        do I=1,M2_
         POMEGA(I,LZ,K) = 0.5d0*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
        end do
       end do

      end if
      end do  ! wavelength loop!

!-----------------------------------------------------------------------
       call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!-----------------------------------------------------------------------

!---Move mean intensity from scatter array FJ(LZ=1:ND)
!---              to CTM mid-level array FJACT(L=1:L_)

      do K=1,W_
      if (FL(K) .gt. 1.0d0) then

!---mean intensity at mid-layer:  4*<I> + solar
!       do L = 1,LU
!        L2L = JNDLEV(L)
!        LZ  = ND+2 - 2*L2L
!        FJACT(L,K) = 4.d0*FJ(LZ,K) + FZ(LZ,K)
!       end do

!---mean intensity averaged throughout layer:
       do L = 1,LU
         LZ0 = ND+2 - 2*JNELEV(L)
        if (L .gt. 1) then
         LZ1 = ND+2 - 2*JNELEV(L-1)
        else
         LZ1 = ND
        end if
         SUMJ = (4.d0*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+2,K)-ZTAU(LZ0,K)) &
               +(4.d0*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-2,K))
         SUMT = ZTAU(LZ0+2,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-2,K)

        do LZ = LZ0+2,LZ1-2,2
         SUMJ =SUMJ+(4.d0*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+2,K)-ZTAU(LZ-2,K))
         SUMT =SUMT + ZTAU(LZ+2,K)-ZTAU(LZ-2,K)
        end do
        FJACT(L,K) = SUMJ/SUMT
       end do

!---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
!---      average (tau-wtd) the h's just above and below the L-edge
       do L = 1,LU
        L2L = JNELEV(L)
        LZ  = ND+2 - 2*L2L
        FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
        FJFLX(L,K)=4.d0*(FJ(LZ-1,K)*FJFLX0 + FJ(LZ+1,K)*(1.d0-FJFLX0))
       end do

!---diffuse fluxes reflected at top, incident at bottom
         FJTOP(K) = FJT(K)
         FJBOT(K) = FJB(K)

      end if
      end do  ! wavelength loop!

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
       end do
      end do

       call LEGND0 (-U0,PM0,M2_)
       do IM=1,M2_
         PM0(IM) = 0.25d0*PM0(IM)
       end do

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
        end do

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
      end if
      end do

      do K = 1,W_
      if (FL(K) .gt. 1.0d0) then
!-----------UPPER BOUNDARY L=1
       L = 1
        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,1,K)
         end do
        end do

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
         end do
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K) &
                    +E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
        end do

!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
         end do
          H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
        end do

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         end do
        end do

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
         end do
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        end do

       end do

!---------FINAL DEPTH POINT: L=ND
       L = ND
        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) &
           + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K) &
           + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
         end do
          H(J,L,K) = H(J,L,K) &
           - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K) &
           - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
        end do

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         end do
        end do

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
        end do

!-----------BACK SOLUTION
       do L = ND-1,1,-1
        do J = 1,M_
         RR(J,L,K) = RR(J,L,K) &
          + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K) &
          + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
        end do
       end do

!----------mean J & H
       do L = 1,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2) &
                + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       end do
       do L = 2,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2) &
                + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       end do

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

      end if
      end do

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
       end do

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
        end do
       end do

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       end do

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        end do
       end do
!-------------upper boundary, 2nd-order, C-matrix is full (CC)
         DELTAU = ZTAU(L2) - ZTAU(L1)
         D2 = 0.25d0*DELTAU
       do I = 1,M_
        do J = 1,M_
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
         CC(I,J,L1) = D2*U(I,J)
        end do
         H(I,L1) = H(I,L1) + 2.0d0*D2*C(I,L1)
         A(I,L1) = 0.0d0
       end do
       do I = 1,M_
        D1 = EMU(I)/DELTAU
        B(I,I,L1)  = B(I,I,L1) + D1
        CC(I,I,L1) = CC(I,I,L1) - D1
       end do

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
        end do
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4) &
          +POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         end do
        end do
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        end do
       end do

!---odd-layer j-points, Legendre terms 1,3,5,7
       do LL=3,ND-2,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
           POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3) &
         + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
        end do
        do I = 1,M_
         do J=1,I
          SUM0 = &
           POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3) &
          +POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         end do
        end do
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        end do
       end do

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
       end do

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
        end do
       end do

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
                + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       end do

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
                + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
                + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        end do
       end do

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
        end do
         H(I,L1) = H(I,L1) - 2.0d0*D2*C(I,L1) + SUM0*ZFLUX
       end do

       do I = 1,M_
          D1 = EMU(I)/DELTAU
        AA(I,I,L1) = AA(I,I,L1) + D1
        B(I,I,L1)  = B(I,I,L1) + D1
        C(I,L1) = 0.0d0
       end do

      END SUBROUTINE GEN_ID

!<<<<<<<<<<<<<<<<<<<<<<<<<<end fastJX core code<<<<<<<<<<<<<<<<<<<<<<<<<




!<<<<<begin fastJX subroutines called from PHOTO_JX or OPMIE<<<<<<<<<<<<

!------------------------------------------------------------------------------
      subroutine OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)
!------------------------------------------------------------------------------
!---set CLOUD fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!      ----FJX ver 7.0+ clouds separated
! 01 W_C02   (C1/Deir)GAMMA:r-m=2.0/alf=6 n=1.335   reff=3.000___G=19.55_rho=1.000
! 02 W_C04   (C1/Deir)GAMMA:r-m=4.0/alf=6 n=1.335   reff=6.000___G=78.19_rho=1.000
! 03 W_C08   (C1/Deir)GAMMA:r-m=8.0/alf=2 n=1.335   reff=12.00___G=301.1_rho=1.000
! 04 W_C13   (C1/Deir)GAMMA:r-m=13./alf=2 n=1.335   reff=20.00___G=472.9_rho=1.000
! 05 W_L06   (W/Lacis)GAMMA:r-m=5.5/alf=11/3        reff=10.00___G=183.9_rho=1.000
! 06 Ice-Hexagonal (Mishchencko)                    reff=50.00___G=999.9_rho=0.917
! 07 Ice-Irregular (Mishchencko)                    reff=50.00___G=999.9_rho=0.917
! 08 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435   reff=0.221___G=.0523_rho=1.630
! 09 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435   reff=0.386___G=.0721_rho=1.630
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
      end if

!--rescale OD by Qext at 600 nm (J=4)
      do J=1,5
         OPTD(J) = ODCLD * QCC(J,NDCLD)/QCC(4,NDCLD)
         SSALB(J) = SCC(J,NDCLD)
        do I=1,8
         SLEG(I,J) =  PCC(I,J,NDCLD)
        end do
      end do

      END SUBROUTINE OPTICL


!------------------------------------------------------------------------------
      subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,K)
!------------------------------------------------------------------------------
!---for the UCI aerosol data sets, calculates optical properties at fast-JX's
!              std 5 wavelengths:200-300-400-600-999nm
!
!---UCI aersols optical data  v-7.0+
! 01 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435   reff=0.221___G=.0523_rho=1.630
! 02 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435   reff=0.386___G=.0721_rho=1.630
! 03 UT-sulfate LOGN:r=0.05 s=.693 n=1.44           reff=0.166___G=.0205_rho=1.769
! 04 UT-sulfate LOGN:r=0.05 s=.693 n=1.46           reff=0.166___G=.0205_rho=1.769
! 05 UT-sulfatM LOGN:r=.050 s=.642 n=1.53           reff=0.140___G=.0179_rho=1.769
! 06 UM-BC1     LOGN:r=.050 s=.642 n=1.80+0.50i     reff=0.140___G=.0179_rho=1.500
! 07 UM-BC2     LOGN:r=.080 s=.501 n=1.80+0.50i     reff=0.150___G=.0332_rho=1.500
! 08 UM-BB08 (%BC)LOGN:r=.080 s=.500 n=1.552+0.04i  reff=0.149___G=.0331_rho=1.230
! 09 UM-FF04 (%BC) LOGN:r=.050 s=.642 n=1.541+0.02i reff=0.140___G=.0179_rho=1.212
! 10 UM-FF10 (%BC)LOGN:r=.050 s=.642 n=1.557+0.05i  reff=0.140___G=.0179_rho=1.230
! 11 MDust.15  (R.V. Martin generated phase fns)    reff=0.150___G=1.000_rho=2.600
! 12 MDust.25  (R.V. Martin generated phase fns)    reff=0.250___G=1.000_rho=2.600
! 13 MDust.40  (R.V. Martin generated phase fns)    reff=0.400___G=1.000_rho=2.600
! 14 MDust.80  (R.V. Martin generated phase fns)    reff=0.800___G=1.000_rho=2.6001
! 15 MDust1.5  (R.V. Martin generated phase fns)    reff=1.500___G=1.000_rho=2.600
! 16 MDust2.5  (R.V. Martin generated phase fns)    reff=2.500___G=1.000_rho=2.600
! 17 MDust4.0  (R.V. Martin generated phase fns)    reff=4.000___G=1.000_rho=2.600
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

      if (K .gt. NAA .or. K .lt. 1) then
         write(6,*) ' aerosol index out-of-range: K/NAA',K,NAA
         K = 18
      end if

         REFF = RAA(K)
         RHO = DAA(K)
      do J=1,5
!---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
         XTINCT = 0.75d0*QAA(J,K)/(REFF*RHO)
         OPTD(J) = PATH*XTINCT
         SSALB(J) = SAA(J,K)
       do I=1,8
         SLEG(I,J) =  PAA(I,J,K)
       end do
      end do

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
      if (L .lt. 1 .or. L .gt. 33) then
!ccc         write(6,*) ' UM aer index too large: L',L
         L = 1
      end if

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
      end do

      END SUBROUTINE OPTICM


!-----------------------------------------------------------------------
      subroutine JRATET(PPJ,TTJ,FFF, VALJL,LU,L_SIZE,NJXU)
!-----------------------------------------------------------------------
! in:
!        PPJ(L_+1) = pressure profile at edges
!        TTJ(L_+1) = = temperatures at mid-level
!        FFF(K=1:NW, L=1:L_) = mean actinic flux
! out:
!        VALJL(L_,JX_)  JX_ = no of dimensioned J-values in CTM code
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: LU,NJXU,L_SIZE
      real*8, intent(in)  ::  PPJ(LU+1),TTJ(LU+1)
      real*8, intent(inout)  ::  FFF(W_,LU)
      real*8, intent(out), dimension(L_SIZE,NJXU) ::  VALJL

      real*8  VALJ(X_)
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV

      if (NJXU .lt. NJX) then
        call EXITC(' JRATET:  CTM has not enough J-values dimensioned')
      end if
      do L = 1,LU
!---need temperature, pressure, and density at mid-layer (for some quantum yields):
          TT   = TTJ(L)
         if (L .eq. 1) then
          PP = PPJ(1)
         else
          PP  = (PPJ(L)+PPJ(L+1))*0.5d0
         end if
          DD = 7.24e18*PP/TT

!---if W_=18/12, must zero bin-11/5 below 100 hPa, since O2 e-fold is too weak
!        and does not represent the decay of 215.5-221.5 nm sunlight.
        if (PP .gt. 100.d0) then
          if (W_ .eq. 18) then
            FFF(11,L) = 0.d0
          elseif (W_ .eq. 12) then
            FFF(5,L) = 0.d0
          end if
        end if

        do J = 1,NJX
          VALJ(J) = 0.d0
        end do

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
         end do

        do J = 4,NJX
         do K = 1,W_
!---also need to allow for Pressure interpolation if SQQ(J) = 'p'
          if (SQQ(J) .eq.'p') then
           call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J), &
              TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          else
           call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J), &
              TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          end if
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
         end do
        end do

        do J=1,NJX
          VALJL(L,J) = VALJ(J)
        end do

      end do

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
        end if
      end if

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
      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP

!      write(6,'(4a)') '   L z(km)     p      T   ', &
!       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
!       '  g(cos) CTM lyr=>'

      L = LU+2
!      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!            L,ZZJ(L)*1.d-5,PPJ(L)

          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)

        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5

!      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!            L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ, &
!            COLO2,COLO3,DTAU6(L),POMEG6(1,L),POMEG6(2,L)/3.d0, &
!            JXTRA(L+L),JXTRA(L+L-1)

        end do

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
!      write(6,'(4a)') '   L z(km)     p      T   ', &
!       '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb', &
!       '  g(cos) CTM lyr=>'
      L = LU+2
!      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!            L,ZZJ(L)*1.d-5,PPJ(L)
          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)
        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5
!      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') &
!            L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ, &
!            COLO2,COLO3
        end do

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
      end if

        RZ(1) = RAD + ZHL(1)
      do II = 2,L1U+1
        RZ(II)   = RAD + ZHL(II)
      end do

!---calculate heights for edges of split CTM-layers
      L2 = 2*L1U
      do II = 2,L2,2
        I = II/2
        RZ2(II-1) = RZ(I)
        RZ2(II) = 0.5d0*(RZ(I)+RZ(I+1))
      end do
        RZ2(L2+1) = RZ(L1U+1)
      do II = 1,L2
        RQ2(II) = (RZ2(II)/RZ2(II+1))**2
      end do

!---shadow height for SZA > 90
      if (U0 .lt. 0.0d0)  then
        SHADHT = RZ2(1)/dsqrt(1.0d0 - U0**2)
      else
        SHADHT = 0.d0
      end if

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
        end do
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
          end if
         end do

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
        end if
      end do

!---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2
      do L2 = L2X,1,-1
        JTOTL  = JTOTL + JXTRA(L2)
        if (JTOTL .gt. NX/2)  then
!          write(6,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
          do L = L2,1,-1
            JXTRA(L) = 0
          end do
          go to 10
        end if
      end do
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
!<<<<<<<<<<<<<<<<<<<<version 7.0+  (3/2013, mjp)<<<<<<<<<<<<<<<<<<<<<<<<

!
! !MODULE: FJX_INIT
!
! !DESCRIPTION: FJX70_INIT contains routines to input fast-JX data
!
!
! !INTERFACE:
!
      MODULE FJX_INIT_MOD
!
! !USES:
!
      USE FJX_CMN_MOD

      USE FJX_SUB_MOD, ONLY : JVN_, NRATJ, JIND, JLABEL, JFACTA, EXITC

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

      integer  JXUNIT,J

!      write(6,*) ' fast-JX ver-7.0 standalone CTM code'

      if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
        call EXITC(' INIT_JX: invalid no. wavelengths')
      end if

! Use channel 8 to read fastJX data files:
      JXUNIT  = 8

! Read in fast-J X-sections (spectral data)
      call RD_XXX(JXUNIT,'FJX_spec.dat')

! Read in cloud scattering data
      call RD_CLD(JXUNIT,'FJX_scat-cld.dat')
! Read in aerosols scattering data
      call RD_MIE(JXUNIT,'FJX_scat-aer.dat')
! Read in UMich aerosol scattering data
      call RD_UM (JXUNIT,'FJX_scat-UMa.dat')

! Read in T and O3 climatology used to fill e.g. upper layers or if O3 not calc.
!pw      call RD_PROF(JXUNIT,'atmos_std.dat')

        NJXX = NJX
      do J = 1,NJX
        TITLEJXX(J) = TITLEJX(J)
      end do

! Read in photolysis rates used in chemistry code and mapping onto FJX J's
!---CTM call:  read in J-values names and link to fast-JX names
      call RD_JS_JX(9,'FJX_j2j.dat', TITLEJXX,NJXX)

      END SUBROUTINE INIT_FJX


!-----------------------------------------------------------------------
      subroutine RD_XXX(NUN,NAMFIL)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
!
!>>>>NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!>>>>NEW v-6.4  changed to collapse wavelengths and x-sections to Trop-only:
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
      real*8  QQ2(199), TQQ2

      character*78 TITLE0
      character*6  TITLEJ2,TITLEJ3
      character*1  TSTRAT

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data------------------
!         note that X_ = max # Xsects read in
!                   NJX = # fast-JX J-values derived from this (.le. X_)

! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

      open (NUN,FILE=NAMFIL,status='old',form='formatted')
      read (NUN,100) TITLE0

!----note that NQRD is not used any more, a read until 'endofJ' is performed
      read (NUN,101) NQRD,NWWW
         NW1 = 1
         NW2 = NWWW
!      write(6,'(1x,a)') TITLE0
!      write(6,'(i8)') NWWW
!----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      read (NUN,102) (WL(IW),IW=1,NWWW)
      read (NUN,102) (FL(IW),IW=1,NWWW)
      read (NUN,102) (QRAYL(IW),IW=1,NWWW)

!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
!---NB the O3 and q-O3-O1D are at different temperatures and cannot be combined
      read (NUN,103) TITLEJX(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

      read (NUN,103) TITLEJX(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

      read (NUN,103) TITLEJX(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

      SQQ(1) = ' '
      SQQ(2) = ' '
      SQQ(3) = ' '

      LQQ(1) = 3
      LQQ(2) = 3
      LQQ(3) = 3

!---Read remaining species:  X-sections at 1-2-3 T_s
        JJ = 3
      do I=1,9999

!--try to read in 3 X-sects per J-value (JJ)
        read (NUN,104) TITLEJ2,TSTRAT,TQQ2,(QQ2(IW),IW=1,NWWW)
        if (TITLEJ2 .eq. 'endofJ') goto 1
!---skip stratosphere only J's (denoted by 'x')if W_<18 => trop-only J's
        if (W_.eq.18 .or. TSTRAT.ne.'x') then
           if (TITLEJ2 .ne. TITLEJX(JJ)) then
              JJ = JJ+1

              if (JJ .gt. X_) then
                 call EXITC(' RD_XXX:  X_ not large enough for Xsects read in')
              end if

              TITLEJX(JJ) = TITLEJ2
              LQQ(JJ) = 1
              SQQ(JJ) = TSTRAT
              LQ = LQQ(JJ)
              TQQ(LQ,JJ) = TQQ2
             do IW = 1,NWWW
                QQQ(IW,LQ,JJ) = QQ2(IW)
             end do
          else
             LQQ(JJ) = LQQ(JJ)+1
             if (LQQ(JJ) .le. 3) then
                LQ = LQQ(JJ)
                TQQ(LQ,JJ) = TQQ2
                do IW = 1,NWWW
                   QQQ(IW,LQ,JJ) = QQ2(IW)
                end do
             end if
           end if
        end if
      end do
    1 continue
        NJX = JJ

      do J = 1,NJX
!        write(6,200) J,TITLEJX(J),SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
!---need to check that TQQ is monotonically increasing:
        if (LQQ(J) .eq. 3) then
          if (TQQ(2,J) .ge. TQQ(3,J)) then
            call EXITC ('TQQ out of order')
          end if
          if (TQQ(1,J) .ge. TQQ(2,J)) then
            call EXITC ('TQQ out of order')
          end if
        end if
        if (LQQ(J) .eq. 2) then
          if (TQQ(1,J) .ge. TQQ(2,J)) then
            call EXITC ('TQQ out of order')
          end if
        end if
      end do

!---check on doingpressure interp
!---check on consolidating Qo2 and others into
!---wrte a newFJX_J2J.dat for mapping on fjx Xsects


!---truncate number of wavelengths to do troposphere-only
      if (W_ .ne. WX_) then
!---TROP-ONLY
       if (W_ .eq. 12) then
!        write(6,'(a)')  &
!         ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
        NW2 = 12
        do IW = 1,4
          WL(IW) = WL(IW+4)
          FL(IW) = FL(IW+4)
          QRAYL(IW) = QRAYL(IW+4)
         do K = 1,3
          QO2(IW,K) = QO2(IW+4,K)
          QO3(IW,K) = QO3(IW+4,K)
          Q1D(IW,K) = Q1D(IW+4,K)
         end do
         do J = 4,NJX
          QQQ(IW,1,J) = QQQ(IW+4,1,J)
          QQQ(IW,2,J) = QQQ(IW+4,2,J)
         end do
        end do
        do IW = 5,12
          WL(IW) = WL(IW+6)
          FL(IW) = FL(IW+6)
          QRAYL(IW) = QRAYL(IW+6)
         do K = 1,3
          QO2(IW,K) = QO2(IW+6,K)
          QO3(IW,K) = QO3(IW+6,K)
          Q1D(IW,K) = Q1D(IW+6,K)
         end do
         do J = 4,NJX
          QQQ(IW,1,J) = QQQ(IW+6,1,J)
          QQQ(IW,2,J) = QQQ(IW+6,2,J)
         end do
        end do
!---TROP-QUICK  (must scale solar flux for W=5)
       elseif (W_ .eq. 8) then
!        write(6,'(a)')   &
!         ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
        NW2 = 8
        do IW = 1,1
          WL(IW) = WL(IW+4)
          FL(IW) = FL(IW+4)  * 2.d0
          QRAYL(IW) = QRAYL(IW+4)
         do K = 1,3
          QO2(IW,K) = QO2(IW+4,K)
          QO3(IW,K) = QO3(IW+4,K)
          Q1D(IW,K) = Q1D(IW+4,K)
         end do
         do J = 4,NJX
          QQQ(IW,1,J) = QQQ(IW+4,1,J)
          QQQ(IW,2,J) = QQQ(IW+4,2,J)
         end do
        end do
        do IW = 2,8
          WL(IW) = WL(IW+10)
          FL(IW) = FL(IW+10)
          QRAYL(IW) = QRAYL(IW+10)
         do K = 1,3
          QO2(IW,K) = QO2(IW+10,K)
          QO3(IW,K) = QO3(IW+10,K)
          Q1D(IW,K) = Q1D(IW+10,K)
         end do
         do J = 4,NJX
          QQQ(IW,1,J) = QQQ(IW+10,1,J)
          QQQ(IW,2,J) = QQQ(IW+10,2,J)
         end do
        end do

       else
         call EXITC(' no. wavelengths wrong: W_ .ne. 8,12,18')
       end if
      end if

      close(NUN)

  100 format(a)
  101 format(10x,5i5)
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
  103 format(a6,1x,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  104 format(a6,a1,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))

      END SUBROUTINE RD_XXX


!-----------------------------------------------------------------------
      subroutine RD_CLD(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
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

      integer  I, J, K
      character*78 TITLE0
      character*20 TITLAA(A_)   ! TITLAA: Title for scatering data

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(i2,a78)') NCC,TITLE0
        if (NCC .gt. C_) then
          call EXITC(' too many cld-data sets: NCC > C_')
        end if

!      write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      read (NUN,*)
      read (NUN,*)
      do J = 1,NCC
          read (NUN,'(3x,a8,1x,2f6.3)') TITLAA(J),RCC(J),DCC(J)
        do K = 1,5
          read (NUN,'(f4.0,f7.4,f7.4,7f6.3)')     &
        WCC(K,J),QCC(K,J),SCC(K,J),(PCC(I,K,J),I=2,8)
          PCC(1,K,J) = 1.d0
        end do
      end do

      close(NUN)

!      write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'  &
!                   ,(WCC(K,1),K=1,5)
!      write(6,*) TITLE0
      do J=1,NCC
!      write(6,'(i3,1x,a8,7f8.3)')   &
!                    J,TITLAA(J),RCC(J),DCC(J),(QCC(K,J),K=1,5)
      end do

      END SUBROUTINE RD_CLD


!-----------------------------------------------------------------------
      subroutine RD_MIE(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
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

      integer  I, J, K
      character*78 TITLE0
      character*20 TITLAA(A_)   ! TITLAA: Title for scatering data

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(i2,a78)') NAA,TITLE0
        if (NAA .gt. A_) then
          call EXITC(' too many aerosol-data sets: NAA > A_')
        end if

!      write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      read (NUN,*)
      read (NUN,*)
      do J = 1,NAA
          read (NUN,'(3x,a8,1x,2f6.3)') TITLAA(J),RAA(J),DAA(J)
        do K = 1,5
          read (NUN,'(f4.0,f7.4,f7.4,7f6.3)')    &
        WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1.d0
        end do
      end do

      close(NUN)

!      write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'  &
!                   ,(WAA(K,1),K=1,5)
!      write(6,*) TITLE0
      do J=1,NAA
!      write(6,'(i3,1x,a8,7f8.3)')   &
!                    J,TITLAA(J),RAA(J),DAA(J),(QAA(K,J),K=1,5)
      end do

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
      character*78 TITLE0
      character*20 TITLUM(33)   ! TITLUM: Title for U Michigan aerosol data set

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(a78)') TITLE0
!        write(6,*) 'UMichigan Aerosols', TITLE0
      read(NUN,'(5x,10f5.0)') WMM
!        write(6,'(a,10f7.1)') ' UMIchigan aerosol wavelengths:',WMM

!---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
!---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
      do L=1,33
          read(NUN,'(a4)') TITLUM(L)
!---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
        do K=1,21
!---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 5=600nm, 6=1000nm
!---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
          read(NUN,'(18f9.5)')  ((UMAER(I,J,K,L),I=1,3),J=1,6)
        end do
      end do

      close(NUN)

!        write(6,'(a)') 'collapse UM wavelengths, drop 550 nm'
          WMM(4) = WMM(5)
          WMM(5) = WMM(6)
       do L=1,33
       do K=1,21
       do I=1,3
          UMAER(I,4,K,L) = UMAER(I,5,K,L)
          UMAER(I,5,K,L) = UMAER(I,6,K,L)
       end do
       end do
       end do

!        write(6,'(7(i5,1x,a4))') (L,TITLUM(L), L=1,33)

      END SUBROUTINE RD_UM


!-----------------------------------------------------------------------
      subroutine RD_PROF(NJ2,NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input T and O3 reference profiles 'atmos_std.dat'
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL

      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real*8  OFAC, OFAK

      character*78 TITLE0

      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
!      write(6,'(1X,A)') TITLE0
!      write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
      end do
      close (NJ2)

!  Extend climatology to 100 km
      OFAC = exp(-2.d5/5.d5)
      do I = 32,51
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          OREF(I,L,M) = OREF(31,L,M)*OFAK
        end do
        end do
      end do
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,51
        TREF(I,L,M) = TREF(41,L,M)
      end do
      end do
      end do


      END SUBROUTINE RD_PROF


!-----------------------------------------------------------------------
      subroutine RD_JS_JX(NUNIT,NAMFIL,TITLEJX,NJX)
!-----------------------------------------------------------------------
!  Read 'FJX_j2j.dat' that defines mapping of fast-JX J's (TITLEJX(1:NJX))
!    onto the CTM reactions:  react# JJ, named T_REACT, uses fast-JX's T_FJX
!    including scaling factor F_FJX.
!-----------------------------------------------------------------------
!---mapping variables stored in  block /jvchem/JFACTA,JIND,NRATJ,JLABEL,JMAP
!           real*8  JFACTA(JVN_)          integer JIND(JVN_), NRATJ
!           character*50 JLABEL(JVN_)     character*6  JMAP(JVN_)
!     JFACTA    multiplication factor for fast-JX calculated J
!     JLABEL    label(*50) of J-value used in the main chem model
!     JMAP      label(*6) of J-value used to match with fast-JX J's
!     NRATJ     number of Photolysis reactions in CTM chemistry, derived here
!                   NRATJ must be .le. JVN_
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)                    ::  NUNIT, NJX
      character(*), intent(in)               ::  NAMFIL
      character*6, intent(in),dimension(NJX) :: TITLEJX
      integer   J,JJ,K
      character*120 CLINE
      character*50 T_REACT
      character*6  JMAP(JVN_), T_FJX
      real*8 F_FJX

! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
! The chemistry code title describes fully the reaction (a50)
! Blank (unfilled) chemistry J's are unmapped
! The number NRATJ is the last JJ readin that is .le. JVN
!   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JMAP(:)   = '------'
      JFACTA(:) = 0.d0

      open (NUNIT,file=NAMFIL,status='old',form='formatted')

       read (NUNIT,'(a)') CLINE
!         write(6,'(a)') CLINE
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ.gt.JVN_) goto 20
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JMAP(JJ) = T_FJX
        NRATJ = JJ
      end do

 20   close(NUNIT)

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
         do J = 1,NJX
            T_FJX = TITLEJX(J)
            if (JMAP(K) .eq. TITLEJX(J)) then
               JIND(K)=J
            end if
         end do
      end do
      !      if(.false.)then
      !      write(6,'(a,i4,a)')' Photochemistry Scheme with',NRATJ,' J-values'
      do K=1,NRATJ
         if (JMAP(K) .ne. '------' ) then
            J = JIND(K)
            if (J.eq.0) then
               !        write(6,'(i5,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K),      &
               !              ' no mapping onto fast-JX',JMAP(K)
            else
               !         write(6,'(i5,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K),   &
               !               ' mapped to FJX:',J,TITLEJX(J)
            end if
         end if
      end do
!      end if
      END SUBROUTINE RD_JS_JX


      END MODULE FJX_INIT_MOD

      module FastJ_mod

        use DefPhotolysis_mod!, only : IDAO3, IDBO3 ....IDACETON
        use GridValues_mod, only : glat,glon,A_bnd,B_bnd,A_mid,B_mid,dA,dB
        use LandDefs_mod,    only: LandType, LandDefs
        use Landuse_mod,    only: LandCover
        use MetFields_mod, only : ps ,foundcloudwater,q,th,lwc,cc3dmax
        use Config_module,    only : KMAX_BND,KMAX_MID,KCHEMTOP, METSTEP
        use NetCDF_mod, only :ReadField_CDF
        use Par_mod,           only: me,LIMAX, LJMAX
        use PhysicalConstants_mod, only :KAPPA, RGAS_KG, GRAV
        use Radiation_mod, only : ZenithAngleS,ZenithAngle
        use TimeDate_mod, only : daynumber,current_date
        use ZchemData_mod, only: rcphot

      USE FJX_CMN_MOD

      USE FJX_SUB_MOD

      USE FJX_INIT_MOD

      implicit none
      PUBLIC

      real*8, allocatable, save :: rcphot_3D(:,:,:,:,:)
      real*8, parameter::  PI180 = 3.141592653589793d0/180.d0
      real*8, parameter::  MASFAC=100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

      real*8, dimension(L1_CTM) :: ETAA_CTM,ETAB_CTM, RI_CTM,TI_CTM,CLDP_CTM,AER1_CTM,AER2_CTM
      integer, dimension(L1_CTM):: NCLD_CTM,NAA1_CTM,NAA2_CTM
      real*8, save, dimension(L1_) :: ETAA,ETAB, RI,TI,CLDP,AER1,AER2,OOO_CTM
      integer, dimension(L1_):: NCLD,NAA1,NAA2
      real*8 GMTAU,ALBEDO,  XLNG,YLAT,XGRD,YGRD
      real*8 PSURF, SCALEH,x

      integer  MONTH
      integer J,K,L, IDAY
!-----------------------------------------------------------------------

!--------------------key params sent to fast JX-------------------------
!---SZA = solar zenith angle, U0 = cos(SZA)
!---SOLF = solar flux factor for sun-earth distance
!---REFLB = Lambertian reflectivity at the Lower Boundary
      real*8                    :: U0,SZA,REFLB,SOLF
!---ZPJ(layer, JV#) = 2D array of J-values, indexed to CTM chemistry code
      real*8, dimension(JVL_,JVN_)::   ZPJ
!---turn internal PHOTO_JX print (unit=6) on/off
      logical                     :: LPRTJ
!---Independent Column Atmosphere data passed to PHOTO_JX:
!--- P = edge press (hPa), Z = edge alt (cm), T = layer temp (K)
!--- D = layer dens-path (# molec /cm2), O = layer O3 path (# O3 /cm2)
!---
!--- R = layer rel.hum.(fraction)
!--- CLDWP = cloud water path (g/m2), AERSP = aerosol path (g/m2)
!--- NDXCLD = cld index type (liq & size, ice & scatt-phase)
!--- NDXAER = aerosol index type
!--- NB. clouds have only a single type within an ICA layer, pick liquid or ice
!---     aerosols are dimensioned with up to AN_ different types in an ICA layer
      real*8,  dimension(L1_+1)  :: PPP,ZZZ
      real*8,  dimension(L1_  )  :: TTT,DDD,RRR,OOO,LWP,IWP,REFFL,REFFI
      real*8,  dimension(L1_,AN_):: AERSP
      integer, dimension(L1_  )  :: NDXCLD
      integer, dimension(L1_,AN_):: NDXAER

      real*8        PCLD,IWC,FACTOR

!---these are the key results coming out of fast-JX core
!---they need to be dimensioned here and passed to fast-JX to fill.
      integer, parameter             ::  NJX_ = 100
      character*6, dimension(NJX_)   ::  TITLJXX
      real*8, dimension(L_,NJX_)     ::  VALJXX
      integer :: NJXX
      integer, save :: L_CLIM_first,L_FastJ

      CONTAINS

      subroutine phot_fastj_interpolate(i_emep,j_emep,errcode)
        integer, intent(in) ::i_emep,j_emep
        integer, intent(inout) ::errcode
        logical, save::first_call=.true.

        real*8 :: weight1,weight2, y, y1, y2, Z, CosZ, latitude, longitude
        real*8 :: thour,thour1,thour2,YGRD,XGRD,SOLF

!part of this could be put somewhere else (in metint and store all results in nr=1?)
        ! we assume the rates followe a sqrt(cos(zenithangle)) function

        latitude = glat(i_emep,j_emep)
        longitude = glon(i_emep,j_emep)
         YGRD = latitude*PI180
         XGRD = longitude*PI180
        thour = real(current_date%hour) + current_date%seconds/3600.0
        thour1 = int((thour+0.000001)/METSTEP)*METSTEP !last meteo time
        thour2 = thour1 + METSTEP !next meteo time (nr=2)
!        call ZenithAngleS( longitude, latitude,daynumber, 365,thour, Z, CosZ )
        call SOLAR_JX(thour,daynumber,YGRD,XGRD, Z,CosZ,SOLF)
        y=sqrt(max(0.0,CosZ))
        call SOLAR_JX(thour1,daynumber,YGRD,XGRD, Z,CosZ,SOLF)
        y1=sqrt(max(0.0,CosZ))
        call SOLAR_JX(thour2,daynumber,YGRD,XGRD, Z,CosZ,SOLF)
        y2=sqrt(max(0.0,CosZ))

!IDAO3 IDNO3 follow a sqrt(cosZ) pattern
        weight1 = (abs(y-y2)+0.5e-8)/(abs(y-y2)+abs(y-y1)+1.0e-8)
        weight2 = 1.0 - weight1 != (abs(y-y1)+0.5e-8)/(abs(y-y2)+abs(y-y1)+1.0e-8)

!        rcphot(:,:) = rcphot_3D(:,:,i_emep,j_emep,1)*weight1 + rcphot_3D(:,:,i_emep,j_emep,2)*weight2
        rcphot(IDAO3,:) = rcphot_3D(IDAO3,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDAO3,:,i_emep,j_emep,2)*weight2
        rcphot(IDNO3,:) = rcphot_3D(IDNO3,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDNO3,:,i_emep,j_emep,2)*weight2
!IDNO2 IDH2O2 IDACH2O(?) IDBCH2O IDHCOHCO IDRCOHCO IDCH3O2H seems to follow a cosZ pattern
        y=y**2
        y1=y1**2
        y2=y2**2
        weight1 = (abs(y-y2)+0.5e-8)/(abs(y-y2)+abs(y-y1)+1.0e-8)
        weight2 = 1.0 - weight1 != (abs(y-y1)+0.5e-8)/(abs(y-y2)+abs(y-y1)+1.0e-8)
        rcphot(IDNO2,:) = rcphot_3D(IDNO2,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDNO2,:,i_emep,j_emep,2)*weight2
        rcphot(IDH2O2,:) = rcphot_3D(IDH2O2,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDH2O2,:,i_emep,j_emep,2)*weight2
        rcphot(IDACH2O,:) = rcphot_3D(IDACH2O,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDACH2O,:,i_emep,j_emep,2)*weight2
        rcphot(IDBCH2O,:) = rcphot_3D(IDBCH2O,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDBCH2O,:,i_emep,j_emep,2)*weight2
        rcphot(IDHCOHCO,:) = rcphot_3D(IDHCOHCO,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDHCOHCO,:,i_emep,j_emep,2)*weight2
        rcphot(IDCH3O2H,:) = rcphot_3D(IDCH3O2H,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDCH3O2H,:,i_emep,j_emep,2)*weight2
        rcphot(IDRCOHCO,:) = rcphot_3D(IDRCOHCO,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDRCOHCO,:,i_emep,j_emep,2)*weight2
!pattern not checked: IDHO2NO2 IDCH3COY IDACETON IDN2O5
        rcphot(IDHO2NO2,:) = rcphot_3D(IDHO2NO2,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDHO2NO2,:,i_emep,j_emep,2)*weight2
        rcphot(IDACETON,:) = rcphot_3D(IDACETON,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDACETON,:,i_emep,j_emep,2)*weight2
        rcphot(IDN2O5,:) = rcphot_3D(IDN2O5,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDN2O5,:,i_emep,j_emep,2)*weight2
        rcphot(IDCH3COY,:) = rcphot_3D(IDCH3COY,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDCH3COY,:,i_emep,j_emep,2)*weight2

!IDBO3 IDHNO3 IDCH3CHO IDCH3COX seems to follow a cosZ**2 pattern
        y=y**2
        y1=y1**2
        y2=y2**2
        weight1 = (abs(y-y2)+0.5e-8)/(abs(y-y2)+abs(y-y1)+1.0e-8)
        weight2 = 1.0 - weight1 != (abs(y-y1)+0.5e-8)/(abs(y-y2)+abs(y-y1)+1.0e-8)
        rcphot(IDBO3,:) = rcphot_3D(IDBO3,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDBO3,:,i_emep,j_emep,2)*weight2
        rcphot(IDHNO3,:) = rcphot_3D(IDHNO3,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDHNO3,:,i_emep,j_emep,2)*weight2
        rcphot(IDCH3CHO,:) = rcphot_3D(IDCH3CHO,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDCH3CHO,:,i_emep,j_emep,2)*weight2
        rcphot(IDCH3COX,:) = rcphot_3D(IDCH3COX,:,i_emep,j_emep,1)*weight1 + rcphot_3D(IDCH3COX,:,i_emep,j_emep,2)*weight2
!        rcphot(ID,:) = rcphot_3D(ID,:,i_emep,j_emep,1)*weight1 + rcphot_3D(ID,:,i_emep,j_emep,2)*weight2

!FOR TESTING
if(me==-6.and.i_emep==5.and.j_emep==5)then
  L=3 !NB: k=kmax_bnd-L
  write(*,22)'AO3 ',    rcphot(IDAO3,kmax_bnd-L) , &
                        rcphot_3D(IDAO3,kmax_bnd-L,i_emep,j_emep,1) ,&
                        rcphot_3D(IDAO3,kmax_bnd-L,i_emep,j_emep,2)
  write(*,22)'BO3 ',    rcphot(IDBO3,kmax_bnd-L)
  write(*,22)'NO2 ',    rcphot(IDNO2,kmax_bnd-L)
  write(*,22)'H2O2 ',   rcphot(IDH2O2,kmax_bnd-L)
  write(*,22)'HNO3 ',   rcphot(IDHNO3,kmax_bnd-L)
  write(*,22)'ACH2O ',  rcphot(IDACH2O,kmax_bnd-L)
  write(*,22)'BCH2O ',  rcphot(IDBCH2O,kmax_bnd-L)
  write(*,22)'CH3CHO ', rcphot(IDCH3CHO,kmax_bnd-L)
  write(*,22)'CH3COX ', rcphot(IDCH3COX,kmax_bnd-L)
  write(*,22)'HCOHCO ', rcphot(IDHCOHCO,kmax_bnd-L)
  write(*,22)'RCOHCO ', rcphot(IDRCOHCO,kmax_bnd-L)

22 format(a,15E12.4)
end if
end subroutine phot_fastj_interpolate

!put values for rcphot for one column at (i_emep,j_emep) in the emep model
subroutine setup_phot_fastj(i_emep,j_emep,errcode,mode)
!mode=0: compute only rcphot values and do not put into rcphot_3D; use nr=1 "now" for fields
!mode=1: put into rcphot_3D; use nr=1 "now" for met fields
!mode=2: rcphot_3D; use nr=2 "future" for met fields
use netcdf

    implicit none
  integer, intent(in) ::i_emep,j_emep
  integer, intent(inout) ::errcode
  integer,  intent(in) ::mode
  integer :: kmid,ilu,lu, nr_local
  real :: Pres_mid,temperature,swp

  logical, save::first_call=.true.
  integer, save::previous_month=-999

  real, allocatable, save :: temperature_clim(:,:,:),cloudliquidwater_clim(:,:,:),humidity_clim(:,:,:),o3_clim(:,:,:)
  real, allocatable, save :: etaa_CLIM(:),etab_CLIM(:)
  integer,save :: Nlevel_CLIM
  character(len=150)::filename,varname
  integer :: LCTM,L_CLIM
  real,parameter :: P0_Pa=101325.0!Pa
  real,parameter :: PS0=1013.250!hPa
  real ::P0
  integer :: ncFileID,varID,dimID
  real*8, dimension(W_)     :: FJBOT,FSBOT

  !-----------------------------------------------------------------------
  !---fast-JX:  INIT_JX is called only once to read in and store all fast-JX data:
  !PW should be moved to emepctm.f90
  if(first_call)then

     call INIT_FJX (TITLJXX,NJX_,NJXX)
     !-----------------------------------------------------------------------

  end if


  nr_local=1
  if(mode>0)nr_local=mode
  !PUT EMEP VALUES in fastj parameters
  IDAY = daynumber
  MONTH = current_date%month
  GMTAU = current_date%hour + current_date%seconds/3600.0
  if(nr_local==2)GMTAU = GMTAU + METSTEP!compute a point in the future



  !Need to describe the entire atmosphere, even if only the J values for the lowest levels are used

  filename='/global/work/mifapw/emep/Data/CLIM/clim2012.nc'!will be put in run.pl in due time

  if(first_call)then
     !find vertical coordinates (size and position compared to model coordinates) and allocate arrays

     !find vertical levels in clim file
     !Note: The pressure at the upper levels from clim file are (almost) identical with the CTM levels
     call check(nf90_open(path=trim(filename),mode=nf90_nowrite,ncid=ncFileID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "lev", dimID = dimID))
     call check(nf90_inquire_dimension(ncid=ncFileID,dimID=dimID,len=Nlevel_CLIM ))
     if(me==0)write(*,*)'Clim number of levels ',Nlevel_CLIM

     allocate(etaa_CLIM(Nlevel_CLIM+1),etab_CLIM(Nlevel_CLIM+1))
     allocate(temperature_clim(Nlevel_CLIM,LIMAX,LJMAX))
     allocate(cloudliquidwater_clim(Nlevel_CLIM,LIMAX,LJMAX))
     allocate(humidity_clim(Nlevel_CLIM,LIMAX,LJMAX))
     allocate(o3_clim(Nlevel_CLIM,LIMAX,LJMAX))
     allocate(rcphot_3D(NRCPHOT,KCHEMTOP:KMAX_MID,LIMAX,LJMAX,2))

     call check(nf90_inq_varid(ncid = ncFileID, name = "P0", varID = varID))
     call check(nf90_get_var(ncFileID, varID, P0 ))
     call check(nf90_inq_varid(ncid = ncFileID, name = "hyai", varID = varID))
     call check(nf90_get_var(ncFileID, varID,  etaa_CLIM ))
     !(clim file uses: P=hyai*101325.0+hybi*PS)
     etaa_CLIM=P0*etaa_CLIM!different definition in model and grid_Def
     call check(nf90_inq_varid(ncid = ncFileID, name = "hybi", varID = varID))
     call check(nf90_get_var(ncFileID, varID, etab_CLIM ))
     call check(nf90_close(ncFileID))
     !level 1 is P=0 corresponds to "L1_+1"

     !find number of levels above emep levels.
     !find first level with pressure at least 100 Pa smaller than topp emep level (i.e with lowest pressure)
     !NB: Pa here (not hPa)
     do L_CLIM=Nlevel_CLIM+1,1,-1
        if((etaa_CLIM(L_CLIM)+ETAB_CLIM(L_CLIM)*110000.0)< A_bnd(1)+B_bnd(1)*110000.0-100.0)exit
     end do
     L_CLIM_first=L_CLIM
     !Ensure that there are no possibilities for level crossing at high mountains
     if((etaa_CLIM(L_CLIM_first)+ETAB_CLIM(L_CLIM_first)*45000.0)> A_bnd(1)+B_bnd(1)*45000.0)then
        write(*,*)'LEVEL CROSSING!',L_CLIM,etaa_CLIM(L_CLIM_first)+ETAB_CLIM(L_CLIM_first)*45000.0,A_bnd(1)+B_bnd(1)*45000.0
        stop
     end if
     !number of FastJ (mid) levels. L_FastJ corresponds to L1_
     L_FastJ=kmax_mid+L_CLIM_first !remember in clim file L=1 is top level
     if(L_<L_FastJ)then
        write(*,*)'too few levels defined in FastJ! ',L_,L_FastJ
        stop
     else
        if(me==0)write(*,*)'FASTJ number of levels ',L_FastJ,' , above emep levels: ',L_CLIM_first
     end if

!define the boundaries of the vertical levels used in FastJ. FastJ counts from surface (L=1) to top (L =  L_FastJ+1)
     L=0
     !first put emep levels
     do k=kmax_bnd,1,-1
        L=L+1
        ETAA(L) = A_bnd(k)/100.0!Pa->hPa
        ETAB(L) = B_bnd(k)
     end do
     !then fill with clim levels
     do L_CLIM=L_CLIM_first,1,-1
        L=L+1
        ETAA(L) = etaa_CLIM(L_CLIM)/100.0!Pa->hPa
        ETAB(L) = ETAB_CLIM(L_CLIM)
     end do
     if(me==0)write(*,*)'FASTJ (standard) pressure levels'
     L=1
     if(me==0)write(*,*)L,ETAA(L)+ETAB(L)*1013.25
     do L=2,L_FastJ+1
        if(me==10)write(*,*)L,ETAA(L)+ETAB(L)*1013.25
        !check for monoticity
        if(ETAA(L)+ETAB(L)*1013.25>=ETAA(L-1)+ETAB(L-1)*1013.25)then
           write(*,*)'1 wrong level'
           stop
        end if
     end do

  end if


  YGRD = glat(i_emep,j_emep)*PI180
  XGRD = glon(i_emep,j_emep)*PI180
  PSURF = ps(i_emep,j_emep,nr_local)/100.0 !Pa-> hPa
  ALBEDO = 0.0
  do ilu= 1, LandCover(i_emep,j_emep)%ncodes
     lu      = LandCover(i_emep,j_emep)%codes(ilu)
     ALBEDO = ALBEDO + LandDefs(lu)%Albedo*0.01*LandCover(i_emep,j_emep)%fraction(ilu)
  end do

  !use fastj vertical direction, i.e. L largest at top, 1 at surface
  !first emep model levels
  L=0
  do k=kmax_mid,1,-1
     L=L+1
     Pres_mid=(A_mid(k) + B_mid(k)*ps(i_emep,j_emep,nr_local))
     !potential -> absolute temperature:
     temperature = th(i_emep,j_emep,k,nr_local)*exp(KAPPA*log(Pres_mid*1.e-5))
     !   write(*,*)'temperature ',L,temperature ,TI(L)
     TI(L) = temperature

     !specific -> relative humidity:
     swp=611.2*exp(17.67*(temperature-273.15)/(temperature-29.65))!saturation water pressure
     !   write(*,*)'humidity ',q(i_emep,j_emep,k,1)*(Pres_mid)/0.622 /swp,RI(L)
     RI(L) = q(i_emep,j_emep,k,nr_local)*(Pres_mid)/0.622 /swp

     !cloud water path (g/m2)
     !NB: cloudwater should be read directly from metdata!
     !can distinguish ice and liquid?
     if(foundcloudwater)then
        !   write(*,*)'CLD ',lwc(i_emep,j_emep,k)*1000*(dA(k)+dB(k)*ps(i_emep,j_emep,1))/GRAV,CLDP(L)
        CLDP(L) = lwc(i_emep,j_emep,k)*1000*(dA(k)+dB(k)*ps(i_emep,j_emep,1))/GRAV!kg/kg -> g/m2
     end if

     ! AERSP(1:L1_,1:AN_)  aerosol path (g/m2)
     ! second index are different types of aerosols?
     !AERSP=0.0 ! to revise
     AER1(L)=0.0
     NAA1(L)=0 !aerosol type?
     AER2(L)=0.0
     NAA2(L)=0 !aerosol type?

  end do


  !fill upper levels with climatological values:
  if(previous_month/=month)then
     !fill levels above emep levels
     !Read climatology
     varname='temperature'
     call  ReadField_CDF(fileName,varname,temperature_clim,nstart=month,kstart=1,kend=Nlevel_CLIM &
          ,interpol='zero_order', needed=.true.,debug_flag=.false.)
     varname='cloudwater'
     call  ReadField_CDF(fileName,varname,cloudliquidwater_clim,nstart=month,kstart=1,kend=Nlevel_CLIM &
          ,interpol='zero_order',needed=.true.,debug_flag=.false.)
     varname='specific_humidity'
     call  ReadField_CDF(fileName,varname,humidity_clim,nstart=month,kstart=1,kend=Nlevel_CLIM &
          ,interpol='zero_order',needed=.true.,debug_flag=.false.)
     varname='O3'
     call  ReadField_CDF(fileName,varname,o3_clim,nstart=month,kstart=1,kend=Nlevel_CLIM &
          ,interpol='zero_order',needed=.true.,debug_flag=.false.)

     previous_month=month
  end if

!NB: PPP changes because PSURF changes
  do L = 1,L_FastJ+1
     PPP(L) = ETAA(L) + ETAB(L)*PSURF
  end do
!  PPP(L_FastJ+1+1)=0.0! intergalactical altitude

  !start from topp
  L_CLIM = 0
  do L=L_FastJ,kmax_mid+1,-1
     L_CLIM = L_CLIM+1
     TI(L) = temperature_clim(L_CLIM,i_emep,j_emep)
     swp=611.2*exp(17.67*(TI(L)-273.15)/(TI(L)-29.65))!saturation water pressure in Pa
     !   write(*,*)'humidity ',q(i_emep,j_emep,kmid,1)*(Pres_mid)/0.622 /swp,RI(L)
     RI(L) =  humidity_clim(L_CLIM,i_emep,j_emep)*100*0.5*(PPP(L)+PPP(L-1))/0.622 /swp
     CLDP(L) = cloudliquidwater_clim(L_CLIM,i_emep,j_emep)
     OOO(L) = o3_clim(L_CLIM,i_emep,j_emep)*(PPP(L)-PPP(L+1))*MASFAC
     if(me==0.and.first_call)write(*,44)'for O3 level',L,' taking level',L_CLIM,' between',&
          etaa_CLIM(L_CLIM)/100+ETAB_CLIM(L_CLIM)*PSURF, ' and ',&
          etaa_CLIM(L_CLIM+1)/100+ETAB_CLIM(L_CLIM+1)*PSURF,' compared to intervall',PPP(L+1),' and ',PPP(L)
44   format(A,I4,A,I4,A,10(F7.2,A))
  end do

  !for O3 we fill also the lower levels with clim values
  !should fill with emep instantaneous values? climatological are more robust.
  !probably does not matter much since most of O3 is higher up
  L_CLIM = Nlevel_CLIM!surface
  do L=1,kmax_mid
     do while ((etaa_CLIM(L_CLIM)/100+ETAB_CLIM(L_CLIM)*PSURF)> PPP(L)-0.001.and.L_CLIM>1)
        L_CLIM = L_CLIM-1
     end do
     if(me==0.and.first_call)write(*,*)'taking O3 from level ',L_CLIM,L
     OOO(L) = o3_clim(L_CLIM,i_emep,j_emep)*(PPP(L)-PPP(L+1))*MASFAC
  end do

  !-----------------------------------------------------------------------
  !---fast-JX:  SOLAR_JX is called only once per grid-square to set U0, etc.
  call SOLAR_JX(GMTAU,IDAY,YGRD,XGRD, SZA,U0,SOLF)

  !-----------------------------------------------------------------------

  LWP(:)  =  0.d0
  IWP(:)  =  0.d0
  REFFL(:) = 0.d0
  REFFI(:) = 0.d0
  AERSP(:,:)  = 0.d0
  NDXAER(:,:) = 0

  do L = 1,L_FastJ
     TTT(L) = TI(L)
     RRR(L) = RI(L)
     if (TTT(L) .gt. 253.d0) then
        LWP(L) = CLDP(L)
     else
        IWP(L) = CLDP(L)
     end if
     NDXAER(L,1) = NAA1(L)
     AERSP(L,1)  = AER1(L)
     NDXAER(L,2) = NAA2(L)
     AERSP(L,2)  = AER2(L)
  end do
  ZZZ(1) = 0.d0
  do L = 1,L_FastJ-1
     DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
     SCALEH      = 1.3806d-19*MASFAC*TTT(L)
     ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH )
  end do
  DDD(L_FastJ)  = (PPP(L_FastJ)-0.0)*MASFAC
  ZZZ(L_FastJ+1) = ZZZ(L_FastJ) + 5.d5
  REFLB = ALBEDO
  LPRTJ = .true.

  !>>>R-effective of clouds determined by main code, not FJX
  !   REFF determined by user - some recommendations below
  !       REFFI is a simple function of ice water content IWC (g/m3, 0.0001 to 0.1)
  !          IWC = IWP / delta-Z (of layer in m, approx OK)
  !   Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a, p.1389
  !              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
  !          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at 1.0)
  !          REFFL is a simple function of pressure (PCLD):
  !            FACTOR = (PCLD - 610.) / 200.
  !            FACTOR = min(1.0, max(0.0, FACTOR))
  !          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)
  do L = 1,L_FastJ
     if (IWP(L) .gt. 1.d-5) then
        IWC = IWP(L) *100.d0 / (ZZZ(L+1)-ZZZ(L))
        IWC = max(0.001d0, IWC)
        REFFI(L) = 164.d0 * IWC**0.23d0
        !       write(6,'(a,i3,3f10.4)') 'ICE:',L,IWP(L),IWC,REFFI(L)
     end if
     if (LWP(L) .gt. 1.d-5) then
        PCLD = 0.5d0*(PPP(L)+PPP(L+1))
        FACTOR = min(1.d0, max(0.d0, (PCLD-610.d0)/200.d0))
        REFFL(L) = 9.60d0*FACTOR + 12.68d0*(1.-FACTOR)
     end if
  end do

  !        call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !---fast-JX:  PHOTO_JX is called once per ICA, calculates all J-values

  debug=0
  call PHOTO_JX                                                 &
       (U0,SZA,REFLB,SOLF, LPRTJ, PPP,ZZZ,TTT,DDD,RRR,OOO,         &
       LWP,IWP,REFFL,REFFI,AERSP,NDXAER,L_FastJ,L1_,AN_,       VALJXX,NJX_, FJBOT,FSBOT,debug)

  !pw: fetch FJBOT,FSBOT diffuse and direct fluxes at bottom (surface). Not 100% sure of definition
  !-----------------------------------------------------------------------

  !---map the J-values from fast-JX onto CTM (ZPJ) using JIND & JFACTA
  do L = 1,L_FastJ-1
     do J = 1,NRATJ
        if (JIND(J).gt.0) then
           ZPJ(L,J) = VALJXX(L,JIND(J))*JFACTA(J)
        else
           ZPJ(L,J) = 0.d0
        end if
     end do
  end do

!TESTING
  if(me==6.and.i_emep==5.and.j_emep==5)then
22   format(a,15E12.4)

     L=18
     if(SZA<90.and.SZA>0)then
        x=sqrt(cos(SZA*PI180))*0.00073
     else
        x=0.0
     end if
     write(*,*)'  k=',kmax_bnd-L
     write(*,22)'PARAMETERS ',         SZA,x,GMTAU,XGRD/PI180,YGRD/PI180!3 O3        PHOTON    O2        O(total)                1.000 /
     write(*,22)'AO3 ',         rcphot(IDAO3,kmax_bnd-L) , ZPJ(L,3), SZA,x,x**2/0.00073,x**4/0.00073/0.00073/0.00073!3 O3        PHOTON    O2        O(total)                1.000 /O3    /
     write(*,22)'BO3 ',          rcphot(IDBO3,kmax_bnd-L) , ZPJ(L,4), SZA,x*90*0.00073,x**2*90,x**4*213017751,x**6*491579426490669.0!4 O3        PHOTON    O2        O(1D)                   1.000 /O3(1D)/
     write(*,22) 'NO2 ',         rcphot(IDNO2,kmax_bnd-L) , ZPJ(L,9),SZA,x*0.00073*27692,x**2*27692,x**4*63905325443.0!9 NO2       PHOTON    N2        O                       1.000 /NO2   /
     write(*,22) 'H2O2 ',         rcphot(IDH2O2,kmax_bnd-L) , ZPJ(L,7),SZA,x*0.00073*22.38,x**2*22.38,x**4*51656804!7 H2O2      PHOTON    OH        OH                      1.000 /H2O2  /
     write(*,22) 'HNO3 ',         rcphot(IDHNO3,kmax_bnd-L) , ZPJ(L,15),SZA,x*0.00073*1.846,x**2*1.846,x**4*4260355!15 HNO3      PHOTON    NO2       OH                      1.000 /HNO3  /
     write(*,22) 'ACH2O ',         rcphot(IDACH2O,kmax_bnd-L) , ZPJ(L,5),SZA,x*0.00073*103.8461,x**2*103.8461,x**4*239644970!5 H2CO      PHOTON    HCO       H                       1.000 /H2COa /
     write(*,22) 'BCH2O ',         rcphot(IDBCH2O,kmax_bnd-L) , ZPJ(L,6),SZA,x*0.00073*161.5,x**2*161.5,x**4*372781065!6 H2CO      PHOTON    CO        H2                      1.000 /H2COb /
     write(*,22) 'CH3CHO ',         rcphot(IDCH3CHO,kmax_bnd-L) , ZPJ(L,54),SZA,x*0.00073*13.84,x**2*13.84,x**4*31952662!54 CH3CHO    PHOTON    CH3       HCO                     1.000 /ActAld/
     write(*,22) 'CH3COX ',         rcphot(IDCH3COX,kmax_bnd-L) , ZPJ(L,61)+ZPJ(L,62),SZA,x*0.00073*13.84,x**2*13.84,x**4*31952662!61 CH3COC2H5 PHOTON    C2H5      CH3CO                   0.850 /MEKeto/
     !62 CH3COC2H5 PHOTON    CH3       C2H5CO                  0.150 /MEKeto/  ?
     write(*,22) 'HCOHCO ',         rcphot(IDHCOHCO,kmax_bnd-L) , ZPJ(L,66),SZA,&
      x*0.00073*50.7692,x**2*50.7692,x**4*117159763, ZPJ(L,67),ZPJ(L,65)
     write(*,22) 'RCOHCO ',         rcphot(IDRCOHCO,kmax_bnd-L) , ZPJ(L,64),SZA,x*0.00073*576.92,x**2*576.92,x**4*1331360946.0!64 CH3COCHO  PHOTON    CH3CO     CO                      1.000 /MGlyxl/

     !should add 11 and 12 or only 12? 11 NO3       PHOTON    NO        O2                      0.114 /NO3   /
     write(*,22) 'IDNO3 ',         rcphot(IDNO3,kmax_bnd-L) , ZPJ(L,11)+ZPJ(L,12),SZA,&
      x*0.00073*623076,x**2*623076,x**4*1437869822485.0!12 NO3       PHOTON    NO2       O                       0.886 /NO3   /

     write(*,22) 'IDCH3O2H ',         rcphot(IDCH3O2H,kmax_bnd-L) , ZPJ(L,8),SZA,x*0.00073*16.6,x**2*16.6,x**4*38343195!8 CH3OOH    PHOTON    CH3O      OH                      1.000 /CH3OOH/
  end if

  if(.not.(allocated(rcphot)))allocate(rcphot(NRCPHOT,KCHEMTOP:KMAX_MID))
!could put directly into rcphot_3D in later version
  do L=1,KMAX_BND-KCHEMTOP
     !hardcoded for now
     !definitions of reactions and indices in FJX_j2j.dat
     !example of interpretation (by PeterW!)
     !  11 NO3       PHOTON    NO        O2                      0.114 /NO3   /
     !  12 NO3       PHOTON    NO2       O                       0.886 /NO3   /
     !NO3 can react with a photon in two ways:
     ! NO3 -> NO+O2 with weight 0.114 and
     ! NO3 -> NO2+O with weight 0.886
     !11 is the index for the first reaction in ZPJ
     !the reaction rate for NO3 reaction is in VALJXX

     !NB: all these mappings MUST be revised!!!

     rcphot(IDAO3,kmax_bnd-L) = ZPJ(L,3)!3 O3        PHOTON    O2        O(total)                1.000 /O3    /
     rcphot(IDBO3,kmax_bnd-L) = ZPJ(L,4)!4 O3        PHOTON    O2        O(1D)                   1.000 /O3(1D)/
     rcphot(IDNO2,kmax_bnd-L) = ZPJ(L,9)!9 NO2       PHOTON    N2        O                       1.000 /NO2   /
     rcphot(IDH2O2,kmax_bnd-L) = ZPJ(L,7)!7 H2O2      PHOTON    OH        OH                      1.000 /H2O2  /
     rcphot(IDHNO3,kmax_bnd-L) = ZPJ(L,15)!15 HNO3      PHOTON    NO2       OH                      1.000 /HNO3  /
     rcphot(IDACH2O,kmax_bnd-L) = ZPJ(L,5)!5 H2CO      PHOTON    HCO       H                       1.000 /H2COa /
     rcphot(IDBCH2O,kmax_bnd-L) = ZPJ(L,6)!6 H2CO      PHOTON    CO        H2                      1.000 /H2COb /

     !not same reaction:
     !emep CH3CHO -> CH3O2+HO2+CO    , fastj: CH3CHO-> CH3 + HCO
     rcphot(IDCH3CHO,kmax_bnd-L) = ZPJ(L,54)!54 CH3CHO    PHOTON    CH3       HCO                     1.000 /ActAld/
     rcphot(IDCH3COX,kmax_bnd-L) = ZPJ(L,61)+ZPJ(L,62)!61 CH3COC2H5 PHOTON    C2H5      CH3CO                   0.850 /MEKeto/
     !62 CH3COC2H5 PHOTON    CH3       C2H5CO                  0.150 /MEKeto/  ?
     rcphot(IDCH3COY,kmax_bnd-L) = ZPJ(L,69)!not used !69 CH3COCH3  PHOTON    CH3       CH3       CO            1.000 /Acet-b/ ?

     !emep: GLYOX=HCOHCO ->CO (1.9) +HO2(0.5)+HCHO(0.1)
     rcphot(IDHCOHCO,kmax_bnd-L) = ZPJ(L,66)!66 CHOCHO    PHOTON    H2        CO        CO            1.000 /Glyxlb/

     rcphot(IDRCOHCO,kmax_bnd-L) = ZPJ(L,64)!64 CH3COCHO  PHOTON    CH3CO     CO                      1.000 /MGlyxl/

     !should add 11 and 12 or only 12? 11 NO3       PHOTON    NO        O2                      0.114 /NO3   /
     rcphot(IDNO3,kmax_bnd-L) = ZPJ(L,11)+ZPJ(L,12)!12 NO3       PHOTON    NO2       O                       0.886 /NO3   /
     rcphot(IDN2O5,kmax_bnd-L) = ZPJ(L,13)!not used !13 N2O5      PHOTON    NO2       NO3                     1.000 /N2O5  /

     rcphot(IDCH3O2H,kmax_bnd-L) = ZPJ(L,8)!8 CH3OOH    PHOTON    CH3O      OH                      1.000 /CH3OOH/
     rcphot(IDHO2NO2,kmax_bnd-L) = ZPJ(L,16)!not used !16 HNO4      PHOTON    NO2       HO2                     1.000 /HNO4  /
     rcphot(IDACETON,kmax_bnd-L) = ZPJ(L,68)!not used !68 CH3COCH3  PHOTON    CH3CO     CH3                     1.000 /Acet-a/

  end do

  if(mode/=0)rcphot_3D(:,:,i_emep,j_emep,nr_local)=rcphot(:,:)

  first_call=.false.

end subroutine setup_phot_fastj

  subroutine check(status)
    use netcdf
    implicit none
    integer, intent ( in) :: status
    if(status /= nf90_noerr)then
    write(*,*)"Error in NetCDF fastJ " //  trim( nf90_strerror(status) )
    stop
    end if
  end subroutine check

    end module FastJ_mod
