! <EQSAM_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-201409 met.no
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
module EQSAM_v03d_ml

 implicit none
 private


 !/- subroutines:
 public   ::  eqsam_v03d


 contains

!subroutine eqsam_v03c(yi,yo,nca,nco,iopt,loop,imax,ipunit,in)
!
!implicit none
!___________________________________________________________________________________________________________________________________
!      Written by Swen Metzger 3/11/99. Modified October 2002, March 2003.
!
!      Department of Atmospheric Chemistry, Max-Planck-Institute for Chemistry.
!      email: metzger@mpch-mainz.mpg.de
!
!      COPYRIGHT 1999-2003
!
!      purpose
!      -------
!      EQSAM is a new and Simplified Aerosol Model, which allows to calculate the gas/aerosol (EQuilibrium) 
!      partitioning, including the aerosol water and aerosol composition suffieciently fast and accurate for 
!      global modeling. EQSAM is based on a parameterization of activcity coefficients (AC), i.e. an AC-RH 
!      relationship, which holds for atmospheric aerosols in equilibrium with the ambient relative humidity (RH).
!      Note that EQSAM should be regarded as a starting point for further development. Although not yet perfect, 
!      it compares rather well with more complex thermodynamic gas/aerosol equilibrium models (EQMs), such as
!      ISORROPIA, or SCAPE.
!      
!      interface
!      ---------
!      call  eqsam_v03b(yi,yo,nca,nco,iopt,loop,imax,ipunit,in)
!
!      yi = input  array (imax, nca)
!      yo = output array (imax, nco)
!      imax = max loop (e.g. time steps)
!      nca >= 11
!      nc0 >= 35
!      iopt = 1 metastable 
!      iopt = 2 solids 
!      iopt = 3 hysteresis (metastable/solids) for online calculations
!      iopt = 31 hysteresis lower branch 
!      iopt = 32 hysteresis upper branch 
!      ipunit = I/O unit (can be skipped)
!      in = array        (can be skipped)
!         
!      method
!      ------
!      equilibrium / internal mixture assumption / aw=rh
!      System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, HCl,Cl-/Na+, H2O 
!              (K+,Ca++,Mg++)
!      external
!      --------
!      program    eqmd.f90 (driver)
!      subroutine gribio.f90  (provides diagnostics output in grib/binary/ascii format)
!      
!      reference
!      ---------
!      Swen Metzger Ph.D Thesis, University Utrecht, 2000
!         http://www.mpch-mainz.mpg.de/~metzger
!
!      Metzger, S. M., F. J. Dentener, J. Lelieveld, and S. N. Pandis, 
!         GAS/AEROSOL PARTITIONING I: A COMPUTATIONALLY EFFICIENT MODEL, 
!         JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 107, NO. D16, 10.1029/2001JD001102, 2002
!      Metzger, S. M., F. J. Dentener, A. Jeuken, and M. Krol, J. Lelieveld, 
!         GAS/AEROSOL PARTITIONING II: GLOBAL MODELING RESULTS, 
!         JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 107, NO. D16, 10.1029/2001JD001103, 2002.
!___________________________________________________________________________________________________________________________

!>-------------------------------------------------------------------------------<
subroutine eqsam_v03d (SO4in, HNO3in,NO3in,NH3in,NH4in,NAin,CLin, relh,temp,pa,   &
                       aSO4out, aNO3out, aNH4out, aNaout, aClout,                 &
                       gSO4out, gNH3out, gNO3out, gClout, aH2Oout) 
!>-------------------------------------------------------------------------------<

  use ModelConstants_ml,  only :  KMAX_MID, KCHEMTOP

implicit none
 real, intent(in):: temp(KCHEMTOP:KMAX_MID),relh(KCHEMTOP:KMAX_MID),  &
                    pa(KCHEMTOP:KMAX_MID)

!hf real :: c(nx,ny,nz,nspec), ah2o(nx,ny,nz)
  real,intent(in)::   &
             SO4in(KCHEMTOP:KMAX_MID),  &
             NO3in(KCHEMTOP:KMAX_MID),  &
             NH4in(KCHEMTOP:KMAX_MID),  &
             NAin (KCHEMTOP:KMAX_MID),  &
             CLin (KCHEMTOP:KMAX_MID),  &
             HNO3in(KCHEMTOP:KMAX_MID), &
             NH3in(KCHEMTOP:KMAX_MID)

   real,intent(out):: &
             aSO4out(KCHEMTOP:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), &
             aNAout (KCHEMTOP:KMAX_MID), &
             aCLout (KCHEMTOP:KMAX_MID), &
             gSO4out(KCHEMTOP:KMAX_MID), &
             gNH3out(KCHEMTOP:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID), &
             gCLout (KCHEMTOP:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!.. local ....
 ! mean value for mixture of wet (2) and dry (1) gridboxes (needed for HYSTERESIS)
real,parameter         :: RH_HIST_DW=1.50                                   
real,parameter         :: T0=298.15, T1=298.0, &
                          AVO=6.03e23,R=82.0567e-6        ! in cu.m*atm/deg/mole
real,parameter         :: RHMAX=0.99, RHMIN=0.0001        ! restrict to max / min RH
real,parameter         :: MWNH4=18., MWSO4=96.,   &       ! mole mass of species considered
                          MWNO3=62., MWCl=35.5,   &           
                          MWNa=23.0, MWH20=55.51*18.01    ! MWCa=40.1,MWN=14.0, MWS=32.1
real,parameter         :: ZERO=0.0
real,parameter         :: GF1=0.25,GF2=0.50,GF3=0.40,GF4=1.00   ! exponents of AC-RH functions
!______________________________________________
integer,parameter                 :: NPAIR=10
!
integer                           :: ii,IHYST, k, iopt
!integer,intent(in)                :: nca,nco,imax,loop,ipunit
!integer,intent(inout)             :: iopt
!______________________________________________
!integer,dimension(6),intent(in)   :: in
!______________________________________________
real                              :: T0T,TT,RH,PX,RHD,KAN,KAC,ZIONIC,RH_HIST,GAMA,GG,GF,GFN
real                              :: X00,X01,X02,X03,X04,X05,X08,X09,X10,X11
real                              :: X0,X1,X2,X3,X4,X5,X6,XK10,XK6
real                              :: ZFLAG,ZKAN,ZKAC,PH,COEF,GAMAAN,HPLUS,AKW,XKW,MOLAL
real                              :: TNH4,TSO4,TNO3,TNa,TCl,TPo,TCa,TMg
real                              :: PNH4,PSO4,PNO3,PCl,PNa,GNO3,GNH3,GSO4,GHCl
real                              :: ASO4,ANO3,ANH4,ACl,ANa,SNH4,SSO4,SNO3,SCl,SNa
real                              :: WH2O,PMt,RINC,DON,RATIONS,GR,NO3P,NH4P !PM,PMs,
!_______________________________________________
!real,dimension(imax,nca),intent(in)  :: yi
!real,dimension(imax,nco),intent(out) :: yo
real,dimension(8)                     :: w1,w2
real,dimension(8)                     :: RHDA,RHDE,RHDX,RHDZ    ! RHD/MRHD for different aerosol types
real,dimension(NPAIR)                 :: M0,MW,NW,ZW            ! arrays of ion pairs
!
! salt solutes:
!   1 = NACl,  2 = (NA)2SO4, 3 = NANO3,  4 = (NH4)2SO4,  5 = NH4NO3, 6 = NH4CL,   7 = 2H-SO4
!   8 = NH4HSO4,   9 = NAHSO4, 10 = (NH4)3H(SO4)2
!
! mole mass of the salt solute
DATA MW(1:NPAIR) / 58.5, 142.0, 88.0, 132.0, 80.0, 53.5, 98.0, 115.0, 120.0, 247.0/
! square of max.  dissocation number (not consistent)
DATA NW(1:NPAIR) / 2.0,   2.5,  2.5,   2.5,  3.5,  1.0,  4.5,   2.0,   2.0,   2.5/ 
! exponents of water activity functions
DATA ZW(1:NPAIR) / 0.67,   1.0,  1.0,   1.0,  1.0,  1.0,  0.5,   1.0,   1.0,   1.0/ 
! RHD / MRHD values as of ISORROPIA / SCAPE (T=298.15K)
DATA RHDA(1:8) / 0.32840, 0.4906, 0.6183, 0.7997, 0.67500, 0.5000, 0.4000, 0.0000/
! Temp. coeff.
DATA RHDE(1:8) / -1860.0, -431.0, 852.00, 80.000, 262.000, 3951.0, 384.00, 0.0000/

logical, parameter :: HYSTERESIS_HISTORY = .false.
!_____________________________________________________________________________________
 IOPT = 1  ! METASTABLE aerosols

IHYST=2
IF(IOPT.EQ.31) THEN      ! SOLID HYSTORY
   IHYST=1
   IOPT=3
ELSEIF(IOPT.EQ.32) THEN  ! WET   HISTORY
   IHYST=2
   IOPT=3
ENDIF

w1=0.;w2=0.        ! init/reset
!______________________________________________________________________________________

 do k=KCHEMTOP,KMAX_MID
! get old relative humidity to calculate aerosol hysteresis (online only)

   RH_HIST = 2.                                        ! WET HISTORY (DEFAULT)
   IF(IHYST.EQ.1.OR.IOPT.EQ.2)  RH_HIST = 1.           ! SET TO SOLIDS

!  meteorology
   TT = temp(k)    ! yi(il,1)      ! T [K]
   RH = relh(k)    ! yi(il,2)      ! RH [0-1]
   PX = pa(k)      ! yi(il,11)     ! p [hPa]
!
! gas+aerosol:
   w1(1) = NAin(k)             !yi(il,6)        ! Na+ (ss  + xsod) (a)   [umol/m^3]
   w1(2) = SO4in(k)            !yi(il,4)        ! H2SO4    + SO4-- (p)   [umol/m^3]
   w1(3) = NH3in(k)+NH4in(k)   !yi(il,3)        ! NH3  (g) + NH4+  (p)   [umol/m^3]
   w1(4) = HNO3in(k)+NO3in(k)  !yi(il,5)        ! HNO3 (g) + NO3-  (p)   [umol/m^3]
   w1(5) = CLin(k)             !yi(il,7)        ! HCl  (g) + Cl-   (p)   [umol/m^3]
   w1(6) = 0. !yi(il, 8)                        ! K+   (p) from Dust     [umol/m^3]
   w1(7) = 0. !yi(il, 9)                        ! Ca++ (p) from Dust     [umol/m^3]
   w1(8) = 0. !yi(il,10)                        ! Mg++ (p) from Dust     [umol/m^3]
!______________________________________________

   zflag=1.

   w1=w1*1.0e-6                     ! [mol/m^3 air]

   TNa   = w1(1)                    ! total input sodium   (g+p) 
   TSO4  = w1(2)                    ! total input sulfate  (g+p) 
   TNH4  = w1(3)                    ! total input ammonium (g+p)
   TNO3  = w1(4)                    ! total input nitrate  (g+p) 
   TCl   = w1(5)                    ! total input chloride (g+p) 
   TPo   = w1(6)                    ! total input potasium (g+p) 
   TCa   = w1(7)                    ! total input calcium  (g+p)
   TMg   = w1(8)                    ! total input magnesium(g+p)

! SULFATE RICH

      if((TNa + TNH4 + TPo +2.*(TCa + TMg)) .le. (2.*TSO4)) then
          zflag=3.
      endif

! SULFATE VERY RICH CASE if (NH4+Na+K+2(Ca+Mg))/SO4 < 1

      if((TNa + TNH4 + TPo +2.*(TCa + TMg)) .le. TSO4) then
          zflag=4.
      endif

! SULFATE NEUTRAL CASE

      if((TNa + TNH4 + TPo +2.*(TCa + TMg)) .gt. (2.*TSO4)) then
          zflag=2.
      endif

! SULFATE POOR AND CATION POOR CASE

      if((TNa + TPo +2.*(TCa + TMg)) .gt. (2.*TSO4)) then       
          zflag=1.
      endif

      IF ( RH .LT. RHMIN ) RH=RHMIN
      IF ( RH .GT. RHMAX ) RH=RHMAX

! CALCULATE TEMPERATURE DEPENDENCY FOR SOME RHDs

      RHDX(:)=RHDA(:)*exp(RHDE(:)*(1./TT-1./T0))
      RHDZ(:)=RHDX(:)
      
! ACCOUNT FOR VARIOUS AMMOMIUM/SODIUM SULFATE SALTS ACCORDING TO MEAN VALUE AS OF ISORROPIA
      GG=2.0            ! (Na)2SO4/(NH4)2SO4 is PREFFERED SPECIES FOR SULFATE DEFICIENT CASES
      IF(ZFLAG.EQ.3.) THEN
         IF(RH.LE.RHDZ(7)) THEN    ! MIXTURE OF (NH4)2SO4(s) & NH4HSO4(s) & (NH4)3H(SO4)2(s) 
            GG=1.677               !  (Na)2SO4 &  NaHSO4
!           GG=1.5
         ELSEIF(RH.GT.RHDZ(7).AND.RH.LE.RHDZ(5)) THEN ! MAINLY (Na)2SO4/(NH4)2SO4(s) & (NH4)3H(SO4)2(s)
            GG=1.75
!           GG=1.5
         ELSEIF(RH.GE.RHDZ(5)) THEN   ! (NH4)2SO4(S) & NH4HSO4(S) & SO4-- & HSO4-
            GG=1.5                    !  (Na)2SO4 &  NaHSO4
         ENDIF
      ENDIF
      IF(ZFLAG.EQ.4.) GG=1.0          ! IF SO4 NEUTRALIZED, THEN ONLY AS NaHSO4/NH4HSO4(S)
                                      !OR HSO4- / H2SO4
      RHD=RH
      IF(IOPT.EQ.2.OR.RH_HIST.LT.RH_HIST_DW) THEN   ! GET RHD FOR SOLIDS / HYSTERESIS
!
! GET LOWEST DELIQUESCENCE RELATIVE HUMIDITIES ACCORDING TO THE CONCENTRATION DOMAIN  
! (APROXIMATION BASED ON RHD / MRHD ISORROPIA/SCAPE
!
      w2(:)=1.
      do ii=1,8
         if(w1(ii).le.1.e-12) w2(ii)=0.    ! skip compound in RHD calculation if 
                enddo                      ! concentration is zero or rather small

! GET LOWEST RHD ACCORDING TO THE CONCENTRATION DOMAIN

! zflag=1. (cation rich)  ...
! 1. sea salt      aerosol          : RHDX(1)=MgCl2
! 2. mineral dust  aerosol          : RHDX(2)=Ca(NO3)2
!
! zflag=2. (sulfate neutral) ...
! 3. ammonium + nitrate             : RHDX(3)= NH4NO3
! 4. ammonium + sulfate             : RHDX(4)=(NH4)2SO4        
! 5. ammonium + sulfate mixed salt  : RHDX(5)=(NH4)3H(SO4)2, (NH4)2SO4        
! 6. ammonium + nitrate  + sulfate  : RHDX(6)=(NH4)2SO4, NH4NO3, NA2SO4, NH4CL
!
! zflag=3. (sulfate poor) ...
! 7. ammonium + sulfate  (1:1,1.5)  : RHDX(7)= NH4HSO4
!
! zflag=4. (sulfate very poor) ...
! 8. sulfuric acid                  : RHDX(8)= H2SO4       

   IF(ZFLAG.EQ.1.)THEN

      RHD=W2(1)+W2(5)                     ! Na+  dependency
      IF(RHD.EQ.0.)  RHDX(1)=1. 
      RHD=W2(6)+W2(7)+W2(8)               ! K+/Ca++/Mg++ dependency (incl. ss)
      IF(RHD.EQ.0.)  RHDX(2)=1. 

      RHD=MINVAL(RHDX(1:2))

   ELSEIF(ZFLAG.EQ.2.)THEN

      RHD=W2(3)*W2(4)                     ! NH4+ & NO3-  dependency
      IF(RHD.EQ.0.)  RHDX(3)=1. 
      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(GG.NE.2.)   RHD=0.               ! account only for pure (NH4)2SO4
      IF(RHD.EQ.0.)  RHDX(4)=1. 
      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(RHD.EQ.0.)  RHDX(5)=1. 
      RHD=W2(2)+W2(3)+W2(4)+W2(5)         ! (NH4)2SO4, NH4NO3, NA2SO4, NH4CL dependency
      IF(RHD.EQ.0.)  RHDX(6)=1. 

!      RHD=MINVAL(RHDX(3:4))
     RHD=MINVAL(RHDX(3:6))

   ELSEIF(ZFLAG.EQ.3.)THEN

      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(RHD.EQ.0.)  RHDX(7)=1. 
      RHD=RHDX(7)                        

   ELSEIF(ZFLAG.EQ.4.)THEN

      RHD=W2(2)                           ! H2SO4 dependency (assume no dry aerosol)
      IF(RHD.EQ.0.)  RHDX(8)=1. 

      RHD=RHDX(8)

   ENDIF ! ZFLAG
   ENDIF ! SOLIDS

   
! GET WATER ACTIVITIES ACCORDING TO METZGER, 2000.
! FUNCTION DERIVED FROM ZSR RELATIONSHIP DATA (AS USED IN ISORROPIA)

      M0(:) = ((NW(:)*MWH20/MW(:)*(1./RH-1.)))**ZW(:)

! CALCULATE TEMPERATURE DEPENDENT EQUILIBRIUM CONSTANTS

      T0T=T0/TT
      COEF=1.0+LOG(T0T)-T0T

! EQUILIBRIUM CONSTANT NH4NO3(s) <==> NH3(g) + HNO3(g)[atm^2] (ISORROPIA)

      XK10 = 5.746e-17
      XK10= XK10 * EXP(-74.38*(T0T-1.0) + 6.120*COEF)
      KAN = XK10/(R*TT)/(R*TT)

! EQUILIBRIUM CONSTANT  NH4CL(s) <==> NH3(g) + HCL(g) [atm^2] (ISORROPIA)

      XK6 = 1.086e-16
      XK6 = XK6 * EXP(-71.00*(T0T-1.0) + 2.400*COEF)
      KAC = XK6/(R*TT)/(R*TT)

! CALCULATE AUTODISSOCIATION CONSTANT (KW) FOR WATER H2O <==> H(aq) + OH(aq) [mol^2/kg^2] (ISORROPIA)

      XKW  = 1.010e-14
      XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)

! GET MEAN MOLAL IONIC ACTIVITY COEFF ACCORDING TO METZGER, 2002.

      GAMA=0.0
      IF(RH.GE.RHD)  GAMA=(RH**ZFLAG/(1000./ZFLAG*(1.-RH)+ZFLAG))
      GAMA = GAMA**GF1            ! ONLY GAMA TYPE OF NH4NO3, NaCl, etc. NEEDED SO FAR

      GAMA=0.0
      GFN=K*K                      ! K=2, i.e. condensation of 2 water molecules per 1 mole ion pair
      GF=GFN*GF1                   ! = GFN[=Nw=4] * GF1[=(1*1^1+1*1^1)/2/Nw=1/4] = 1
                                   ! ONLY GAMA TYPE OF NH4NO3, NH4Cl, etc. needed so far

      IF(RH.GE.RHD) GAMA=RH**GF/((GFN*MWH20*(1./RH-1.)))**GF1

      GAMA = min(GAMA,1.0)        ! FOCUS ON 0-1 SCALE
      GAMA = max(GAMA,0.0)
      GAMA = (1.-GAMA)**K          ! transplate into aqueous phase equillibrium and account for 
                                   ! enhanced uptake of aerosol precursor gases with increasing RH
                                   ! (to match the results of ISORROPIA)


! CALCULATE RHD DEPENDENT EQ: IF RH <  RHD => NH4NO3(s) <==> NH3 (g) + HNO3(g) (ISORROPIA)
!                             IF RH >> RHD => HNO3  (g)   -> NO3 (aq)

      X00  = MAX(ZERO,MIN(TNa,GG*TSO4))       ! MAX SODIUM   SULFATE
      X0   = MAX(ZERO,MIN(TNH4,GG*TSO4-X00))  ! MAX AMMOMIUM SULFATE
      X01  = MAX(ZERO,MIN(TNa-X00, TNO3))     ! MAX SODIUM   NITRATE
      X1   = MAX(ZERO,MIN(TNH4-X0,TNO3-X01))  ! MAX AMMOMIUM NITRATE
!
      X02  = MAX(ZERO,MIN(TNa-X01-X00,TCl))   ! MAX SODIUM   CHLORIDE
      X03  = MAX(ZERO,MIN(TNH4-X0-X1,TCl-X02))! MAX AMMOMIUM CHLORIDE

      X2   = MAX(TNH4-X1-X0-X03,ZERO)         ! INTERIM RESIDUAL NH3
      X3   = MAX(TNO3-X1-X01,ZERO)            ! INTERIM RESIDUAL HNO3
      X04  = MAX(TSO4-(X0+X00)/GG,ZERO)       ! INTERIM RESIDUAL H2SO4
      X05  = MAX(TCl-X03-X02,ZERO)            ! INTERIM RESIDUAL HCl
!     X06  = MAX(TNa-X02-X01-X00,ZERO)        ! INTERIM RESIDUAL Na (should be zero for electro-neutrality in input data)
!

      ZKAN=2.
      IF(RH.GE.RHD) ZKAN=ZKAN*GAMA

      X4   = X2 + X3
!corrected SM      X5   = SQRT(X4*X4+KAN*ZKAN)
      X5   = SQRT(X4*X4+KAN*ZKAN*ZKAN) 
      X6   = 0.5*(-X4+X5)
      X6   = MIN(X1,X6)
      
      GHCl = X05                              ! INTERIM RESIDUAl HCl
      GNH3 = X2 + X6                          ! INTERIM RESIDUAl NH3
      GNO3 = X3 + X6                          ! RESIDUAl HNO3
      GSO4 = X04                              ! RESIDUAl H2SO4
      PNa  = X02 + X01 + X00                  ! RESIDUAl Na (neutralized)
      
      ZKAC=2.
      IF(RH.GE.RHD) ZKAC=ZKAC*GAMA

      X08   = GNH3 + GHCl
      X09   = SQRT(X08*X08+KAC*ZKAC*ZKAC)
      X10   = 0.5*(-X08+X09)
      X11   = MIN(X03,X10)

      GHCl = GHCl + X11                       ! RESIDUAL HCl
      GNH3 = GNH3 + X11                       ! RESIDUAL NH3

! GO SAVE ...

      IF(GHCl.LT.0.)   GHCl=0.
      IF(GSO4.LT.0.)   GSO4=0.
      IF(GNH3.LT.0.)   GNH3=0.
      IF(GNO3.LT.0.)   GNO3=0.
      IF(PNa.LT.0.)    PNa=0.
      IF(GSO4.GT.TSO4) GSO4=TSO4
      IF(GNH3.GT.TNH4) GNH3=TNH4
      IF(GNO3.GT.TNO3) GNO3=TNO3
      IF(GHCl.GT.TCl)  GHCl=TCl
      IF(PNa.GT.TNa)   PNa=TNa
!     IF(PNa.LT.TNa)   print*,il,' PNa.LT.TNa => no electro-neutrality in input data! ',PNa,TNa


! DEFINE AQUEOUSE PHASE (NO SOLID NH4NO3 IF NO3/SO4>1, TEN BRINK, ET AL., 1996, ATMOS ENV, 24, 4251-4261)

!     IF(TSO4.EQ.ZERO.AND.TNO3.GT.ZERO.OR.TNO3/TSO4.GE.1.) RHD=RH

!     IF(IOPT.EQ.2.AND.RH.LT.RHD.OR.IOPT.EQ.2.AND.RH_HIST.LT.RH_HIST_DW) THEN        ! SOLIDS / HYSTERESIS
      IF(RH_HIST.EQ.1.AND.RH.LT.RHD) THEN        ! SOLIDS / HYSTERESIS

       ! EVERYTHING DRY, ONLY H2SO4 (GSO4) REMAINS IN THE AQUEOUSE PHASE

         ANH4 = 0.
         ASO4 = 0.
         ANO3 = 0.
         ACl  = 0.
         ANa  = 0.

      ELSE  !  SUPERSATURATED SOLUTIONS NO SOLID FORMATION

         ASO4 = TSO4 - GSO4
         ANH4 = TNH4 - GNH3
         ANO3 = TNO3 - GNO3
         ACl  = TCl  - GHCl
         ANa  = PNa

      ENDIF ! SOLIDS / HYSTERESIS

! CALCULATE AEROSOL WATER [kg/m^3(air)]
!
! salt solutes:
!   1 = NACl,  2 = (NA)2SO4, 3 = NANO3,  4 = (NH4)2SO4,  5 = NH4NO3, 6 = NH4CL,   7 = 2H-SO4
!   8 = NH4HSO4,   9 = NAHSO4, 10 = (NH4)3H(SO4)2
!
      WH2O = 1.0e-6   !  small initial value 
      IF(ZFLAG.EQ.1.) WH2O = ASO4/M0( 2) + ANO3/M0(3) +  ACl/M0(6)
      IF(ZFLAG.EQ.2.) WH2O = ASO4/M0( 4) + ANO3/M0(5) +  ACl/M0(6)
      IF(ZFLAG.EQ.3.) WH2O = ASO4/M0( 8) + ANO3/M0(5) +  ACl/M0(6)
      IF(ZFLAG.EQ.4.) WH2O = ASO4/M0( 8) + GSO4/M0(7)


! CALCULATE AQUEOUS PHASE PROPERTIES

!     PH    = 9999.
      PH    = 7.
      MOLAL = 0.
      HPLUS = 0.
      ZIONIC= 0.

!hf&pw      IF(WH2O.GT.0.) THEN            
      IF(WH2O.GT.1.0e-6) THEN            

         ! CALCULATE AUTODISSOCIATION CONSTANT (KW) FOR WATER

         AKW=XKW*RH*WH2O*WH2O                  ! H2O <==> H+ + OH- with kw [mol^2/kg^2]
         AKW=AKW**0.5                          ! [OH-] = [H+] [mol]

         ! Calculate hydrogen molality [mol/kg], i.e. H+ of the ions:
         !           Na+, NH4+, NO3-, Cl-, SO4--, HH-SO4- [mol/kg(water)]
         !                                   with [OH-] = kw/[H+]

         HPLUS = (-ANa/WH2O-ANH4/WH2O+ANO3/WH2O+ACl/WH2O+GG*ASO4/WH2O+GG*GSO4/WH2O+ & 
            SQRT(( ANa/WH2O+ANH4/WH2O-ANO3/WH2O-ACl/WH2O-GG*ASO4/WH2O-GG*GSO4/WH2O)**2 &
                  +XKW/AKW*WH2O))/2.

         ! Calculate pH

         PH=-ALOG10(HPLUS)                             ! aerosol pH 

         ! Calculate ionic strength [mol/kg]

         ZIONIC=0.5*(ANa+ANH4+ANO3+ACl+ASO4*GG*GG+GSO4*GG*GG+XKW/AKW*WH2O*WH2O)
         ZIONIC=ZIONIC/WH2O                            ! ionic strength [mol/kg]
!        ZIONIC=min(ZIONIC,200.0)                      ! limit for output
!        ZIONIC=max(ZIONIC,0.0)

      ENDIF ! AQUEOUS PHASE 

!
!-------------------------------------------------------
! calculate diagnostic output consistent with other EQMs ...
!
      ASO4 = ASO4 + GSO4                        ! assuming H2SO4 remains aqueous

      TNa   = TNa  * 1.e6                       ! total input sodium   (g+p)  [umol/m^3]
      TSO4  = TSO4 * 1.e6                       ! total input sulfate  (g+p)  [umol/m^3]
      TNH4  = TNH4 * 1.e6                       ! total input ammonium (g+p)  [umol/m^3]
      TNO3  = TNO3 * 1.e6                       ! total input nitrate  (g+p)  [umol/m^3]
      TCl   = TCl  * 1.e6                       ! total input chloride (g+p)  [umol/m^3]
      TPo   = TPo  * 1.e6                       ! total input potasium (g+p)  [umol/m^3]
      TCa   = TCa  * 1.e6                       ! total input calcium  (g+p)  [umol/m^3]
      TMg   = TMg  * 1.e6                       ! total input magnesium(g+p)  [umol/m^3]
!
! residual gas:
      GNH3 = GNH3 * 1.e6                        ! residual NH3
      GSO4 = GSO4 * 1.e6                        ! residual H2SO4
      GNO3 = GNO3 * 1.e6                        ! residual HNO3
      GHCl = GHCl * 1.e6                        ! residual HCl

! total particulate matter (neutralized)
      PNH4=TNH4-GNH3                            ! particulate ammonium   [umol/m^3]
      PNO3=TNO3-GNO3                            ! particulate nitrate    [umol/m^3]
      PCl =TCl -GHCl                            ! particulate chloride   [umol/m^3]
      PNa =TNa                                  ! particulate sodium     [umol/m^3]
      PSO4=TSO4                                 ! particulate sulfate    [umol/m^3]

! liquid matter
      ANH4 = ANH4 * 1.e6                        ! aqueous phase ammonium [umol/m^3]
      ANO3 = ANO3 * 1.e6                        ! aqueous phase nitrate  [umol/m^3]
      ACl  = ACl  * 1.e6                        ! aqueous phase chloride [umol/m^3]
      ANa  = ANa  * 1.e6                        ! aqueous phase sodium   [umol/m^3]
      ASO4 = ASO4 * 1.e6                        ! aqueous phase sulfate  [umol/m^3] 

! solid matter
      SNH4=PNH4-ANH4                            ! solid phase ammonium   [umol/m^3]
      SSO4=PSO4-ASO4                            ! solid phase sulfate    [umol/m^3]
      SNO3=PNO3-ANO3                            ! solid phase nitrate    [umol/m^3]
      SCl =PCl -ACl                             ! solid phase chloride   [umol/m^3]
      SNa =PNa -ANa                             ! solid phase sodium     [umol/m^3]

! GO SAVE ...

      IF(SNH4.LT.0.)   SNH4=0.
      IF(SSO4.LT.0.)   SSO4=0.
      IF(SNO3.LT.0.)   SNO3=0.
      IF(SCl.LT.0.)    SCl=0.
      IF(SNa.LT.0.)    SNa=0.

 ! PM=SNH4+SSO4+SNO3+SNH4+SCl+SNa+ANH4+ASO4+ANO3+ACl+ANa  ! total PM [umol/m^3]
 ! PMs=SNH4*MWNH4+SSO4*MWSO4+SNO3*MWNO3+SCl*MWCl+SNa*MWNa ! dry PM   [ug/m^3]
 ! PMt=PMs+ANH4*MWNH4+ASO4*MWSO4+ANO3*MWNO3+ACl*MWCl+ ANa*MWNa ! dry+wet PM, excl.H20[ug/m^3]

      WH2O = WH2O * 1.e9                   ! convert aerosol water from [kg/m^3] to [ug/m^3]
      IF(WH2O.LT.1.e-3) WH2O=0.

! UPDATE HISTORY RH FOR HYSTERESIS (ONLINE CALCULATIONS ONLY) - not tested here!!!!
     if  (HYSTERESIS_HISTORY) then
      RH_HIST=2.                                     ! wet
      IF(WH2O.EQ.0.) RH_HIST=1.                      ! dry

! Approximate the pH (for test purposes only)
      PH    = 7.
      HPLUS = 0.
      IF(WH2O.GT.0.) &
       HPLUS=(2.*TSO4+ANO3+ACl-ANH4-ANa)/WH2O*1000.  ! hydrogen ion concentration [mol/l]
      IF(HPLUS.GT.0.) PH=-ALOG10(HPLUS)              ! aerosol pH 

      ZIONIC=0.
      IF(WH2O.GT.0.) ZIONIC=0.5*(ANa+ANH4+ANO3+ACl+ASO4*4.)  ! ionic strength [moles/kg]
      ZIONIC=ZIONIC*1.e3/WH2O
      ZIONIC=min(ZIONIC,200.0)                               ! limit for output
      ZIONIC=max(ZIONIC,0.0)

      GAMAAN=0.0
      IF(WH2O.GT.0.) GAMAAN = GAMA**GF1             ! activity coefficient (NH4NO3)
      GAMAAN=min(GAMAAN,1.0)                        ! focus on 0-1 scale
      GAMAAN=max(GAMAAN,0.0)

      RINC = 1.
      IF(PMt.GT.0.)   RINC = (WH2O/PMt+1)**(1./3.)  ! radius increase due to water uptake
      IF(RINC.EQ.0.)  RINC = 1.

      RATIONS = 0.
      IF(PSO4.GT.0.) RATIONS = PNO3/PSO4            ! nitrate / sulfate mol ratio

      GR = 0.
      IF(GNO3.GT.0.) GR = GNH3/GNO3                 ! gas ratio=residual NH3/residual HNO3[-]

      DON = 0.
      IF((PNO3+2.*PSO4).GT.0.) DON = 100.*PNH4/(PNO3+2.*PSO4) ! degree of neutralization
                                      ! by ammonia : ammonium / total nitrate + sulfate  [%]
      NO3P = 0.
      IF(TNO3.GT.0.) NO3P = 100.*PNO3/TNO3          ! nitrate partitioning=nitrate/total nitrate[%]

      NH4P = 0.
      IF(TNH4.GT.0.) NH4P = 100.*PNH4/TNH4          ! ammonium partitioning=ammonium/total ammonium[%]

!     KAN   = rks5/(r*temp)**2                      ! Keq of NH3(g)+HNO3(g)---> NH4NO3 (s) 
                                  ! [mol^2/kg]/(R[m^3*atm/deg/mole]*T[K])**2 = [m^3*atm/kg]
  endif
!
! store aerosol species for diagnostic output:
!______________________________________________________________

! Output values:
!//.. aerosols
      aSO4out(k) = PSO4                    ! particulate sulfate  (p=a+s)     [umol/m^3]
      aNO3out(k) = PNO3                    ! particulate nitrate  (p=a+s)     [umol/m^3]
      aNH4out(k) = PNH4                    ! particulate ammonium (p=a+s)     [umol/m^3]
      aNAout(k)  = PNa                     ! particulate sodium (p=a+s)       [umol/m^3]
      aClout(k)  = PCl                     ! particulate chloride (p=a+s)     [umol/m^3]
!//.. gases
      gSO4out(k) = GSO4                    ! residual H2SO4 (aq)              [umol/m^3]
      gNO3out(k) = GNO3                    ! residual HNO3  (g)               [umol/m^3]
      gNH3out(k) = GNH3                    ! residual NH3   (g)               [umol/m^3]
      gCLout(k)  = GHCL                    ! residual HCl   (g)               [umol/m^3]

!//.. aerosol water
      aH2Oout(k) = WH2O                    ! aerosol Water  (aq)              [ug/m^3]

 enddo
!
 end subroutine eqsam_v03d

end module EQSAM_v03d_ml
