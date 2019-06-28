! <EQSAM4clim_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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
module EQSAM4clim_ml

 implicit none
 private


 !/- subroutines:
 public   ::  EQSAM4clim


 contains

!____________________________________________________________________________________________
! SUBROUTINE EQSAM4clim(xTT,xAW,xsPM,xaPM,xPMS,xPMT,xRHO,xVOL,xPH,xHp,xGF,xWH2O,&
!                       xYPa,xYPs,xYMa,xYMs,xYG,imask,neq,nleq,jmeq,Dd)
!____________________________________________________________________________________________
!  IMPLICIT NONE
!  SAVE
!____________________________________________________________________________________________
! WRITTEN BY SWEN METZGER (s.metzger@cyi.ac.cy), 2012-2014
!          The Cyprus Institute, www.cyi.ac.cy
!               *** COPYRIGHT 2012-2018+ ***
!
! Contact, w.r.t. EQSAM4clim modeling:
! Swen Metzger <swen.metzger@eco-serve.de>
! Founder, ResearchConcepts io GmbH, HRB 717519
! http://www.eco-serve.de/en/Publications.html
! http://www.researchconcepts.io
!
! Metzger, S., B. Steil, L. Xu, J. E. Penner, and J. Lelieveld, 
! New representation of water activity based on a single solute specific constant to 
! parameterize the hygroscopic growth of aerosols in atmospheric models, 
! Atmos. Chem. Phys., 12, 5429–5446, doi:10.5194/acp-12-5429-2012, 2012; 
! http://www.atmos-chem-phys.net/12/5429/2012/acp-12-5429-2012.html
!
! Metzger, S., B. Steil, M. Abdelkader, K. Klingmüller, L. Xu, J.E. Penner, C. Fountoukis, 
! A. Nenes, and J. Lelieveld; Aerosol Water Parameterization: A single parameter framework; 
! Atmos. Chem. Phys., 16, 7213–7237, doi:10.5194/acp-16-7213-2016, 2016; 
! https://www.atmos-chem-phys.net/16/7213/2016/acp-16-7213-2016.html
!
! EUMETSAT ITT 15/210839 Validation Report, 12/2016: 
! Comparison of Metop - PMAp Version 2 AOD Products using Model Data
! http://www.eumetsat.int/website/wcm/idc/idcplg?IdcService=GET_FILE&dDocName=
! PDF_PMAP_V2_MODEL_COMP&RevisionSelectionMethod=LatestReleased&Rendition=Web
!____________________________________________________________________________________________
!  CHARACTER(LEN=*),PARAMETER :: modstr  = 'EQSAM4clim'    ! name of module
!  CHARACTER(LEN=*),PARAMETER :: modver  = 'v10'           ! module version
!  CHARACTER(LEN=*),PARAMETER :: moddat  = '15Jan2018'     ! last modification date
!____________________________________________________________________________________________
! ' Hydrogen |   H2O           H2SO4          HNO3            HCl            NH3    '
! '    H+    |    1              2             3               4              5     '
! '  Index   |  eqh2o          eqhsa          eqhna          eqhca          eqxam   '
! '----------|----------------------------------------------------------------------'
! ' Ammonium |(NH4)3H(SO4)2  (NH4)2SO4       NH4HSO4        NH4NO3          NH4Cl   '
! '   NH4+   |    6              7             8               9             10     '
! '  Index   |  eqalc          eqasu          eqahs          eqano          eqacl   '
! '----------|----------------------------------------------------------------------'
! ' Sodium   | Na3H(SO4)2     Na2SO4          NaHSO4         NaNO3          NaCl    '
! '   Na+    |   11             12            13              14             15     '
! '  Index   |  eqslc          eqssu          eqshs          eqsno          eqscl   '
! '----------|----------------------------------------------------------------------'
! ' Potassium|  K3H(SO4)2      K2SO4           KHSO4          KNO3           KCl    '
! '    K+    |   16             17            18              19             20     '
! '  Index   |  eqplc          eqpsu          eqphs          eqpno          eqpcl   '
! '----------|----------------------------------------------------------------------'
! ' Calcium  |   ---           CaSO4          ---           Ca(NO3)2        CaCl2   '
! '   Ca++   |   21             22            23              24             25     '
! '  Index   |  eqc01          eqcsu          eqc02          eqcno          eqccl   '
! '----------|----------------------------------------------------------------------'
! ' Magnesium|   ---           MgSO4          ---           Mg(NO3)2        MgCl2   '
! '   Mg++   |   26             27            28              29             30     '
! '  Index   |  eqm01          eqmsu          eqm02          eqmno          eqmcl   '
!____________________________________________________________________________________________
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Start modification for EMEP
!>-------------------------------------------------------------------------------<
subroutine EQSAM4clim (SO4in, HNO3in,NO3in,NH3in,NH4in,NAin,CLin, relh,temp,   &
                       aSO4out, aNO3out, aNH4out, aNaout, aClout,              &
                       gSO4out, gNH3out, gNO3out, gClout, aH2Oout, KCHEMTOP, KMAX_MID) 
!>-------------------------------------------------------------------------------<

!  use Config_module,  only :  KMAX_MID, KCHEMTOP

implicit none
 integer, intent(in):: KCHEMTOP, KMAX_MID
 real, intent(in)   :: temp(KCHEMTOP:KMAX_MID),relh(KCHEMTOP:KMAX_MID)

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

!... Swen's comments ....
! neq = inner loop and should go over your number of grid boxes (e.g., vector loop), 
!       e.g., certain number of lon or lat, or time. 
! nleq1 = outer loop and could go via levels, e.g.  nleq1=KCHEMTOP, nleq2=KMAX_MID 
!         (though the order doesn’t really matter)
   integer,parameter                          :: nleq1 = 1
   integer,parameter                          :: nleq2 = 1
   integer,parameter                          :: nleq = nleq2
!   integer,parameter                          :: neq1 = KCHEMTOP
!   integer,parameter                          :: neq2 = KMAX_MID  ! KMAX_MID isn't parameter
   integer                                    :: neq1, neq2
!   integer,parameter                          :: neq  = neq2
!   integer :: neq 
   integer,parameter                          :: jmeq = 1
 
!
! End modification for EMEP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!______________________________________________
INTEGER,PARAMETER                             :: dp = SELECTED_REAL_KIND(12,307)
!______________________________________________
REAL(dp),PARAMETER                            :: REALZERO=tiny(0._dp),ZERO=0._dp,ONE=1._dp
REAL(dp),PARAMETER                            :: TINYX=1.e-15_dp,eqT0=298.15_dp ![K]
REAL(dp),PARAMETER                            :: R=8.314409_dp            ![J/mol/K]
REAL(dp),PARAMETER                            :: eqR=R/101325._dp   ![atm*m^3/mol/K]
REAL(dp),PARAMETER                            :: sigma=0.0761_dp    ![J/m^2]
!______________________________________________
 LOGICAL,PARAMETER                            :: lke         =.FALSE.
 LOGICAL,PARAMETER                            :: lvola       =.TRUE.
 LOGICAL,PARAMETER                            :: lmixs       =.FALSE.
 LOGICAL,PARAMETER                            :: lrhdm       =.TRUE.
 LOGICAL,PARAMETER                            :: lHSO4       =.TRUE.
 LOGICAL,PARAMETER                            :: lH2SO4gas   =.FALSE.
 LOGICAL,PARAMETER                            :: lmetastable =.TRUE.  !.FALSE.
!______________________________________________
INTEGER,PARAMETER :: jMs=1
INTEGER,PARAMETER :: jDs=2
INTEGER,PARAMETER :: jZa=3
INTEGER,PARAMETER :: jWs=4
INTEGER,PARAMETER :: jNs=5
INTEGER,PARAMETER :: jNi=6
INTEGER,PARAMETER :: jRHD=7
INTEGER,PARAMETER :: jRHDc=8
INTEGER,PARAMETER :: jSOL=30
INTEGER,PARAMETER :: jVOL=2
INTEGER,PARAMETER :: jGAS=4
INTEGER,PARAMETER :: jKEQ=2
INTEGER,PARAMETER :: jAP=1
INTEGER,PARAMETER :: jDP=2
INTEGER,PARAMETER :: jGP=3
INTEGER,PARAMETER :: jSP=1
INTEGER,PARAMETER :: jSM=2
INTEGER,PARAMETER :: jSS=3
!______________________________________________
INTEGER,PARAMETER :: mciam =  1  ! ammonium        -  NH4+
INTEGER,PARAMETER :: mciso =  2  ! sodium          -  Na+
INTEGER,PARAMETER :: mcipo =  3  ! potassium       -  K+
INTEGER,PARAMETER :: mcica =  4  ! calcium         -  Ca++
INTEGER,PARAMETER :: mcimg =  5  ! magnesium       -  Mg++
INTEGER,PARAMETER :: mcati =  5
!______________________________________________
INTEGER,PARAMETER :: mdumm =  0  ! dummy
INTEGER,PARAMETER :: maisu =  1  ! sulfate         -  SO4--
INTEGER,PARAMETER :: maihs =  2  ! bisulfate       -  HSO4-
INTEGER,PARAMETER :: maino =  3  ! nitrate         -  NO3-
INTEGER,PARAMETER :: maicl =  4  ! chloride        -  Cl-
INTEGER,PARAMETER :: manio =  4
!______________________________________________
CHARACTER(len= 5), DIMENSION(0:jSOL)   :: ceqsolute
INTEGER,           DIMENSION(0:jSOL,5) :: ieqsolute
INTEGER,           DIMENSION(1:jVOL)   :: ieqvola
INTEGER,           DIMENSION(1:jGAS)   :: ieqgases
INTEGER,           DIMENSION(KCHEMTOP:KMAX_MID,nleq) :: imask
!_______________________________________________
LOGICAL,           DIMENSION(0:jSOL)   :: leqsolute
!_______________________________________________
INTEGER                                :: i,n,m,II,IJ,IK,Jp,Jm,Ip,Im,Is
!INTEGER                               :: neq,neq1,neq2,nleq,nleq1,nleq2,jmeq ! modified for EMEP
!______________________________________________
REAL(dp)                               :: T0T,AKW,X1,X2,X3,XZ,KEQ,Mw,Dw,COEF,Dd
!______________________________________________
REAL(dp),DIMENSION(0:mcati)  :: eqxpz,eqxpm
REAL(dp),DIMENSION(0:manio)  :: eqxmz,eqxmm
REAL(dp),DIMENSION(0:jSOL)   :: eqxyZ,eqxyA,eqxyB,eqxyC
REAL(dp),DIMENSION(8,0:jSOL) :: eqxy
!______________________________________________
LOGICAL :: leqskip_all,lhelp
!______________________________________________
LOGICAL, DIMENSION(KCHEMTOP:KMAX_MID) :: leqskip
!______________________________________________
INTEGER, DIMENSION(KCHEMTOP:KMAX_MID) :: NZeq,IWATeq,IDeq
!______________________________________________
INTEGER, DIMENSION(KCHEMTOP:KMAX_MID,0:jSOL) :: ieqsalt
!______________________________________________
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID) :: eqTT,eqAW,eqRH,eqsPM,eqPMs,eqaPM,eqPMt,eqWH2O,eqPH,eqGF,YY,YZ
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID) :: eqVOL,eqRHO,eqXPi,eqXMi,eqHPLUS,XRHDMIN,XRHDMAX,XRHDIFF,TSO4
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID) :: eqMolm,eqMs,eqNs,eqWs,eqNZm,eqRHDM
!_______________________________________________
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,mcati,2)   :: eqyp
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,manio,2)   :: eqym
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,jGAS)      :: eqyg
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,0:jSOL,3)  :: eqys
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,0:jSOL)    :: eqxyRHD
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,0:jSOL)    :: eqxyXS,eqxyMol,eqxyKe,eqxyGF
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,0:mcati)   :: eqxp
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,0:manio)   :: eqxm
!______________________________________________
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,nleq,4)    :: xYG
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,nleq,4)    :: xYMa,xYMs
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,nleq,5)    :: xYPa,xYPs
REAL(dp),DIMENSION(KCHEMTOP:KMAX_MID,nleq)      :: xTT,xAW,xsPM,xPMs,xaPM,xPMt,xWH2O,xPH,xVOL,xRHO,xGF,xHp
!____________________________________________________________________________________________
! cation name     |            mciam      mciso      mcipo      mcica       mcimg
DATA eqxpz(0:mcati) / 1._dp,  1.0_dp,    1.0_dp,    1.0_dp,     2.0_dp,     2.0_dp     /
DATA eqxpm(0:mcati) / 1._dp,0.01805_dp,0.02299_dp,0.039098_dp,0.040080_dp,0.024305_dp /
!____________________________________________________________________________________________
! anion  name     |            maisu      maihs      maino      maicl
DATA eqxmz(0:manio) / 1._dp,  2.0_dp,    1.0_dp,    1.0_dp,    1.0_dp /
DATA eqxmm(0:manio) / 1._dp,0.09607_dp,0.09708_dp,0.062010_dp,0.03545_dp/
!____________________________________________________________________________________________
INTEGER,PARAMETER ::     eqdum =  0, &
          eqh2o =  1,    eqhsa =  2,    eqhna =  3,    eqhca =  4,    eqxam =  5,&
          eqalc =  6,    eqasu =  7,    eqahs =  8,    eqano =  9,    eqacl = 10,&
          eqslc = 11,    eqssu = 12,    eqshs = 13,    eqsno = 14,    eqscl = 15,&
          eqplc = 16,    eqpsu = 17,    eqphs = 18,    eqpno = 19,    eqpcl = 20,&
          eqc01 = 21,    eqcsu = 22,    eqc02 = 23,    eqcno = 24,    eqccl = 25,&
          eqm01 = 26,    eqmsu = 27,    eqm02 = 28,    eqmno = 29,    eqmcl = 30
!____________________________________________________________________________________________
DATA ieqvola (1:jVOL)  /eqacl,eqano/
DATA ieqgases(1:jGAS)  /eqhna,eqhca,eqxam,eqhsa/
!____________________________________________________________________________________________
INTEGER                ::     ieqsu(mcati)
DATA                          ieqsu/eqasu,eqssu,eqpsu,eqcsu,eqmsu/
INTEGER                ::     ieqhs(mcati)
DATA                          ieqhs/eqahs,eqshs,eqphs,mdumm,mdumm/
INTEGER                ::     ieqno(mcati)
DATA                          ieqno/eqano,eqsno,eqpno,eqcno,eqmno/
INTEGER                ::     ieqcl(mcati)
DATA                          ieqcl/eqacl,eqscl,eqpcl,eqccl,eqmcl/
!____________________________________________________________________________________________
INTEGER                ::     ieqams(manio)
DATA                          ieqams/eqasu,eqahs,eqano,eqacl/
INTEGER                ::     ieqsos(manio)
DATA                          ieqsos/eqssu,eqshs,eqsno,eqscl/
INTEGER                ::     ieqpos(manio)
DATA                          ieqpos/eqpsu,eqphs,eqpno,eqpcl/
INTEGER                ::     ieqcas(manio)
DATA                          ieqcas/eqcsu,mdumm,eqcno,eqccl/
INTEGER                ::     ieqmgs(manio)
DATA                          ieqmgs/eqmsu,mdumm,eqmno,eqmcl/
!____________________________________________________________________________________________
DATA leqsolute(0:jSOL)/.FALSE.,&
         .FALSE.,       .FALSE.,       .FALSE.,       .FALSE.,       .FALSE.,&
         .TRUE.,        .TRUE.,        .TRUE.,        .TRUE.,        .TRUE.,&
         .FALSE.,       .TRUE.,        .TRUE.,        .TRUE.,        .TRUE.,&
         .FALSE.,       .TRUE.,        .TRUE.,        .TRUE.,        .TRUE.,&
         .FALSE.,       .TRUE.,        .FALSE.,       .TRUE.,        .TRUE.,&
         .FALSE.,       .TRUE.,        .FALSE.,       .TRUE.,        .TRUE./
!____________________________________________________________________________________________
DATA ieqsolute(0:jSOL,jSS)/mdumm,&
          eqh2o,         eqhsa,         eqhna,         eqhca,         eqxam,&
          eqalc,         eqasu,         eqahs,         eqano,         eqacl,&
          eqslc,         eqssu,         eqshs,         eqsno,         eqscl,&
          eqplc,         eqpsu,         eqphs,         eqpno,         eqpcl,&
          eqc01,         eqcsu,         eqc02,         eqcno,         eqccl,&
          eqm01,         eqmsu,         eqm02,         eqmno,         eqmcl/
!____________________________________________________________________________________________
DATA ieqsolute(0:jSOL,jSP)/mdumm,&
          mdumm,         mdumm,         mdumm,         mdumm,         mciam,&
          mciam,         mciam,         mciam,         mciam,         mciam,&
          mciso,         mciso,         mciso,         mciso,         mciso,&
          mcipo,         mcipo,         mcipo,         mcipo,         mcipo,&
          mdumm,         mcica,         mcica,         mcica,         mcica,&
          mdumm,         mcimg,         mdumm,         mcimg,         mcimg/
!____________________________________________________________________________________________
DATA ieqsolute(0:jSOL,jSM)/mdumm,&
          mdumm,         maisu,         maino,         maicl,         mdumm,&
          maisu,         maisu,         maihs,         maino,         maicl,&
          maisu,         maisu,         maihs,         maino,         maicl,&
          maisu,         maisu,         maihs,         maino,         maicl,&
          mdumm,         maisu,         mdumm,         maino,         maicl,&
          mdumm,         maisu,         mdumm,         maino,         maicl/
!____________________________________________________________________________________________
DATA ceqsolute(0:jSOL)/ 'eqdum',                                             &
         'eqh2o',       'eqhsa',       'eqhna',       'eqhca',       'eqxam',&
         'eqalc',       'eqasu',       'eqahs',       'eqano',       'eqacl',&
         'eqslc',       'eqssu',       'eqshs',       'eqsno',       'eqscl',&
         'eqplc',       'eqpsu',       'eqphs',       'eqpno',       'eqpcl',&
         'eqc01',       'eqcsu',       'eqc02',       'eqcno',       'eqccl',&
         'eqm01',       'eqmsu',       'eqm02',       'eqmno',       'eqmcl'/
!____________________________________________________________________________________________
!        Ms   [kg/mol],&
DATA eqxy(jMs,0:jSOL)/  ZERO, &
       0.018020_dp,   0.098090_dp,   0.063020_dp,   0.036460_dp,   0.017040_dp,&
       0.247300_dp,   0.132170_dp,   0.115130_dp,   0.080060_dp,   0.053500_dp,&
       0.262120_dp,   0.142050_dp,   0.120070_dp,   0.085000_dp,   0.058440_dp,&
       0.310444_dp,   0.174266_dp,   0.136178_dp,   0.101108_dp,   0.074548_dp,&
       0.000000_dp,   0.136150_dp,   0.000000_dp,   0.164100_dp,   0.110980_dp,&
       0.000000_dp,   0.120375_dp,   0.000000_dp,   0.148325_dp,   0.095205_dp/
!____________________________________________________________________________________________
!        Ds   [kg/m^3],&
DATA eqxy(jDs,0:jSOL)/  ZERO, &
     997.000000_dp,1830.000000_dp,1513.000000_dp,1490.000000_dp, 696.000000_dp,&
    1775.000000_dp,1770.000000_dp,1780.000000_dp,1720.000000_dp,1519.000000_dp,&
    2565.000000_dp,2700.000000_dp,2430.000000_dp,2260.000000_dp,2170.000000_dp,&
    2490.000000_dp,2660.000000_dp,2320.000000_dp,2110.000000_dp,1988.000000_dp,&
    1700.000000_dp,2960.000000_dp,1700.000000_dp,2500.000000_dp,2150.000000_dp,&
    1000.000000_dp,2660.000000_dp,1000.000000_dp,2300.000000_dp,2325.000000_dp/
!____________________________________________________________________________________________
!        Za        [-],&
DATA eqxy(jZa,0:jSOL)/  ZERO, &
       1.000000_dp,   2.000000_dp,   1.000000_dp,   1.000000_dp,   1.000000_dp,&
       3.000000_dp,   2.000000_dp,   1.000000_dp,   1.000000_dp,   1.000000_dp,&
       3.000000_dp,   2.000000_dp,   1.000000_dp,   1.000000_dp,   1.000000_dp,&
       3.000000_dp,   2.000000_dp,   1.000000_dp,   1.000000_dp,   1.000000_dp,&
       3.000000_dp,   2.000000_dp,   3.000000_dp,   2.000000_dp,   2.000000_dp,&
       1.000000_dp,   2.000000_dp,   1.000000_dp,   2.000000_dp,   2.000000_dp/
!____________________________________________________________________________________________
!        Ws (T_o)  [%],&
DATA eqxy(jWs,0:jSOL)/  ZERO, &
       0.000000_dp,  70.000000_dp,  25.000000_dp,  15.000000_dp,  30.000000_dp,&
      53.300000_dp,  43.310000_dp,  76.000000_dp,  68.050000_dp,  28.340000_dp,&
      44.060000_dp,  21.940000_dp,  66.180000_dp,  47.700000_dp,  26.470000_dp,&
       0.000000_dp,  10.710000_dp,  33.600000_dp,  27.690000_dp,  26.230000_dp,&
       0.000000_dp,   0.210000_dp,   0.000000_dp,  59.020000_dp,  44.840000_dp,&
       0.000000_dp,  26.310000_dp,   0.000000_dp,  41.590000_dp,  35.900000_dp/
!____________________________________________________________________________________________
!        nu_s     [-],&
DATA eqxy(jNs,0:jSOL)/  ZERO, &
       2.000000_dp,   3.000000_dp,   2.000000_dp,   2.000000_dp,   1.000000_dp,&
       5.000000_dp,   3.000000_dp,   2.000000_dp,   2.000000_dp,   2.000000_dp,&
       5.000000_dp,   3.000000_dp,   2.000000_dp,   2.000000_dp,   2.000000_dp,&
       5.000000_dp,   3.000000_dp,   2.000000_dp,   2.000000_dp,   2.000000_dp,&
       5.000000_dp,   2.000000_dp,   5.000000_dp,   3.000000_dp,   3.000000_dp,&
       1.000000_dp,   2.000000_dp,   1.000000_dp,   3.000000_dp,   3.000000_dp/
!____________________________________________________________________________________________
!        nu_i     [-],&
DATA eqxy(jNi,0:jSOL)/  ZERO, &
       1.000000_dp,   1.761309_dp,   1.793689_dp,   2.681333_dp,   0.527416_dp,&
       1.616356_dp,   1.274822_dp,   1.253573_dp,   1.051480_dp,   1.243054_dp,&
       1.000000_dp,   1.278762_dp,   1.293906_dp,   1.160345_dp,   1.358377_dp,&
       1.000000_dp,   1.286445_dp,   1.308499_dp,   1.014102_dp,   1.256989_dp,&
       1.000000_dp,   1.271828_dp,   1.000000_dp,   1.586562_dp,   2.024869_dp,&
       1.000000_dp,   1.435281_dp,   1.000000_dp,   1.878693_dp,   2.107772_dp/
!____________________________________________________________________________________________
!        RHD (T_o) [-],&
DATA eqxy(jRHD,0:jSOL)/  ZERO, &
       1.000000_dp,   1.000000_dp,   1.000000_dp,   1.000000_dp,   1.000000_dp,&
       0.690000_dp,   0.799700_dp,   0.400000_dp,   0.618300_dp,   0.771000_dp,&
       0.390000_dp,   0.930000_dp,   0.520000_dp,   0.737900_dp,   0.752800_dp,&
       1.000000_dp,   0.975000_dp,   0.860000_dp,   0.924800_dp,   0.842600_dp,&
       1.000000_dp,   0.990000_dp,   1.000000_dp,   0.490600_dp,   0.283000_dp,&
       1.000000_dp,   0.861300_dp,   1.000000_dp,   0.540000_dp,   0.328400_dp/
!____________________________________________________________________________________________
!        T-coef.   [-],&
DATA eqxy(jRHDc,0:jSOL)/  ZERO, &
       0.000000_dp,   0.000000_dp,   0.000000_dp,   0.000000_dp,   0.000000_dp,&
     186.000000_dp,  80.000000_dp, 384.000000_dp, 852.000000_dp, 239.000000_dp,&
       0.000000_dp,  80.000000_dp, -45.000000_dp, 304.000000_dp,  25.000000_dp,&
       0.000000_dp,  35.600000_dp,   0.000000_dp,   0.000000_dp, 159.000000_dp,&
       0.000000_dp,   0.000000_dp,   0.000000_dp, 509.400000_dp, 551.100000_dp,&
       0.000000_dp,-714.450000_dp,   0.000000_dp, 230.200000_dp,  42.230000_dp/
!____________________________________________________________________________________________
! GLOBAL INITIALIZATION 
!    nleq1=1                             ! modified for EMEP
!    nleq2=nleq                          ! modified for EMEP
!    neq1=1                              ! modified for EMEP
!    neq2=neq                            ! modified for EMEP

 neq1 = KCHEMTOP
 neq2 = KMAX_MID 
!   neq  = neq2
!______________________________________
   DO n=nleq1,nleq2    ! n=1,1
!______________________________________
   DO i=neq1,neq2      ! i=KCHEMTOP,KMAX_MID
      eqXPi(i)=ZERO
      eqXMi(i)=ZERO
      eqxp (i,0)=ZERO
      eqxp (i,1)=ZERO
      eqxp (i,2)=ZERO
      eqxp (i,3)=ZERO
      eqxp (i,4)=ZERO
      eqxp (i,5)=ZERO
      eqxm (i,0)=ZERO
      eqxm (i,1)=ZERO
      eqxm (i,2)=ZERO
      eqxm (i,3)=ZERO
      eqxm (i,4)=ZERO
      eqyg (i,1)=ZERO
      eqyg (i,2)=ZERO
      eqyg (i,3)=ZERO
      eqyg (i,4)=ZERO
      leqskip(i)=.TRUE.
      imask(i,n) = 1    ! comment line for 3-D implementation
      xYPa (i,n,:)=ZERO ! modified for EMEP
      xYMa (i,n,:)=ZERO ! modified for EMEP
      xWH2O(i,n)  =ZERO ! modified for EMEP
   END DO
!______________________________________
! INPUT - modified for EMEP
!    DO i=neq1,neq2
!       eqWH2O(i)=xWH2O(i,n)
!       eqTT  (i)=xTT  (i,n)
!       eqAW  (i)=xAW  (i,n)
!       eqRH  (i)=eqAW (i)
!   END DO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Start INPUT for EMEP:
   DO i=neq1,neq2     ! i = vertikal levels
      eqWH2O(i)=xWH2O(i,n) ! aerosol water (from previous time step)
      eqTT  (i)=temp (i)   ! T [K]
      eqAW  (i)=relh (i)   ! water activity [0-1]
      eqRH  (i)=eqAW (i)   ! equilibrium assumption (aw=RH)
      xYPa(i,n,mciam) = (NH3in(i)+NH4in(i))*1.E-6_dp    ! NH3  (g) + NH4+  (p)   [mol/m^3]
      xYPa(i,n,mciso) =  NAin(i)*1.E-6_dp               ! Na+  (ss + xsod) (p)   [mol/m^3]
      xYMa(i,n,maisu) =  SO4in(i)*1.E-6_dp              ! H2SO4(g) + SO4- + HSO4-- (p) [mol/m^3]
      xYMa(i,n,maino) = (HNO3in(i)+NO3in(i))*1.E-6_dp   ! HNO3 (g) + NO3-  (p)   [mol/m^3]
      xYMa(i,n,maicl) =  CLin(i)*1.E-6_dp               ! HCl  (g) + Cl-   (p)   [mol/m^3]
   END DO
!  END INPUT for EMEP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! CATIONS [MOL/M^3(air)]
   DO IJ=1,mcati
   DO i=neq1,neq2
      eqxp (i,IJ)    =xYPa(i,n,IJ)
      eqXPi(i)=eqXPi(i)+eqxp(i,IJ)*eqxpz(IJ)
   END DO
   END DO ! mcati
   ! ANIONS [MOL/M^3(air)]
   DO IJ=1,manio
   DO i=neq1,neq2
      eqxm (i,IJ)    =xYMa(i,n,IJ)
      eqXMi(i)=eqXMi(i)+eqxm(i,IJ)*eqxmz(IJ)
   END DO
   END DO ! manio
    leqskip_all=.TRUE.
   DO i=neq1,neq2
      IF(eqXPi(i)+eqXMi(i) > TINYX) leqskip(i)=.FALSE.
      IF(imask (i,n) == 0) leqskip(i) =.TRUE.
      IF(.NOT.leqskip(i))  leqskip_all=.FALSE.
   END DO
!______________________________________
   IF(leqskip_all) GOTO 1000
!______________________________________
! DOMAINS
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      IWATeq(i)=2 ! SOLID+LIQUID AEROSOL
      IF(lmetastable.OR.eqWH2O(i) > REALZERO) &
      IWATeq(i)=1 ! METASTABLE AEROSOL
      TSO4(i) = eqxm(i,maihs)+eqxm(i,maisu)
      X1 = eqxmz (maihs)*eqxm(i,maihs)
      X2 = eqxmz (maisu)
      IF(lHSO4.AND.eqXPi(i)>ZERO.AND.eqXPi(i)<X1+X2) THEN
         Is=eqasu
         IF(eqxp(i,mciam) >= eqxm(i,maisu)*0.5_dp) Is=eqahs
         IF(eqxp(i,mcipo) >= eqxm(i,maisu)*0.5_dp) Is=eqphs
         IF(eqxp(i,mciso) >= eqxm(i,maisu)*0.5_dp) Is=eqshs
         IF(eqRH(i) >= eqxy(jRHD,Is)) X2=1.95_dp
      END IF
      X2 =            X2*eqxm(i,maisu)
      X3 = eqxpz (mciam)*eqxp(i,mciam)
      IDeq(i) = 0 ! SULFATE POOR
      IF(eqXPi(i) < TINYX .AND. TSO4(i) > TINYX) THEN
         IDeq(i) = 5 ! ONLY SULFURIC ACID
         eqxm(i,maisu)=eqxm(i,maihs)+eqxm(i,maisu)
         eqxm(i,maihs)=ZERO
      ELSE IF(eqXPi(i) > TINYX .AND. eqXPi(i) < TSO4(i)) THEN
         IDeq(i) = 4 ! SULFATE VERY RICH
         eqxm(i,maihs)=eqxm(i,maihs)+eqxm(i,maisu)
         eqxm(i,maisu)=ZERO
      ELSE IF(eqXPi(i) >= TSO4(i) .AND. eqXPi(i) < X1+X2) THEN
         IDeq(i) = 3 ! SULFATE RICH
         eqxm(i,maisu)=eqxm(i,maihs)+eqxm(i,maisu)
         eqxm(i,maihs)=ZERO
         XZ = MAX(ZERO,X1+X2-eqXPi(i))
         eqxm(i,maihs)=eqxm(i,maihs)+XZ
         eqxm(i,maisu)=eqxm(i,maisu)-XZ
      ELSE IF(eqXPi(i) >= X1+X2 .AND. eqXPi(i)-X3 < TSO4(i)) THEN
         IDeq(i) = 2 ! SULFATE NEUTRAL
      ELSE IF(eqXPi(i)-X3 >= TSO4(i)) THEN
         IDeq(i) = 1 ! SULFATE POOR / MINERAL CATION RICH
      END IF
      IF(IDeq(i) <= 2) THEN
         eqxm(i,maisu)=eqxm(i,maisu)+eqxm(i,maihs)
         eqxm(i,maihs)=ZERO
      END IF
   END DO
!______________________________________
! LOCAL INITIALIZATION
   Mw=eqxy(jMs,eqh2o)
   Dw=eqxy(jDs,eqh2o)
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      YY     (i)=ZERO
      eqWH2O (i)=ZERO
      eqMolm (i)=ZERO
      eqMs   (i)=ZERO
      eqNs   (i)=ZERO
      eqWs   (i)=ZERO
      eqNZm  (i)=ZERO
      NZeq   (i)=0
      eqPH   (i)=ZERO
      eqHPLUS(i)=ZERO
      eqsPM  (i)=ZERO
      eqaPM  (i)=ZERO
      eqPMs  (i)=ZERO
      eqPMt  (i)=ZERO
      eqVOL  (i)=ZERO
      eqGF   (i)=ONE
      eqRHO  (i)=ONE
      eqRHDM (i)=ONE
      eqXPi  (i)=ZERO
      eqXMi  (i)=ZERO
      XRHDMAX(i)=ZERO
      XRHDMIN(i)=ONE
      XRHDIFF(i)=ONE
   END DO
   DO IJ=1,mcati
      DO i=neq1,neq2
        eqyp(i,IJ,1) = ZERO
        eqyp(i,IJ,2) = ZERO
      END DO
   END DO
   DO IJ=1,manio
      DO i=neq1,neq2
        eqym(i,IJ,1) = ZERO
        eqym(i,IJ,2) = ZERO
      END DO
   END DO
   DO IS=0,jSOL
      eqxyA(Is)=ONE
      eqxyB(Is)=ZERO
      eqxyC(Is)=ZERO
      eqxyZ(Is)=ZERO
      DO i=neq1,neq2
         ieqsalt(i,Is)=0
         eqxyGF (i,Is)=ONE
         eqxyRHD(i,Is)=ONE
         eqxyKe  (i,Is)=ONE
         eqxyXS  (i,Is)=ONE
         eqxyMol (i,Is)=ZERO
         eqys(i,Is,jAP)=ZERO
         eqys(i,Is,jDP)=ZERO
         eqys(i,Is,jGP)=ZERO
      END DO
   END DO
!______________________________________
! NEUTRALIZATION/REACTION ORDER
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      IF(IDeq(i) < 3) THEN
          ieqsalt(i, 1)=eqcsu ; ieqsalt(i, 2)=eqmsu ; ieqsalt(i, 3)=eqpsu
          ieqsalt(i, 4)=eqssu ; ieqsalt(i, 5)=eqasu ; ieqsalt(i, 6)=eqcno
          ieqsalt(i, 7)=eqmno ; ieqsalt(i, 8)=eqpno ; ieqsalt(i, 9)=eqsno
          ieqsalt(i,10)=eqano ; ieqsalt(i,11)=eqccl ; ieqsalt(i,12)=eqmcl
          ieqsalt(i,13)=eqpcl ; ieqsalt(i,14)=eqscl ; ieqsalt(i,15)=eqacl
      END IF
      IF(IDeq(i) == 3) THEN
          ieqsalt(i, 1)=eqcsu ; ieqsalt(i, 2)=eqmsu
          ieqsalt(i, 3)=eqpsu ; ieqsalt(i, 4)=eqphs
          ieqsalt(i, 5)=eqssu ; ieqsalt(i, 6)=eqshs
          ieqsalt(i, 7)=eqasu ; ieqsalt(i, 8)=eqahs
      END IF
      IF(IDeq(i) == 4) THEN
          ieqsalt(i, 1)=eqcsu ; ieqsalt(i, 2)=eqmsu
          ieqsalt(i, 3)=eqphs ; ieqsalt(i, 4)=eqshs
          ieqsalt(i, 5)=eqahs ; ieqsalt(i, 6)=eqhsa
          ieqsalt(i, 7)=eqalc
      END IF
      IF(IDeq(i) == 5) THEN
          ieqsalt(i, 1)=eqhsa
          ieqsalt(i, 2)=eqalc
      END IF
   END DO
!______________________________________
! SOLUTE MOLALITY
   DO IJ=1,jSOL
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         Is=ieqsalt(i,IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         eqxyB(Is)=ZERO
         eqxyZ(Is)=eqxy(jNi,Is)
         eqxyMol(i,Is)=((eqxyKe(i,Is)/eqRH(i)-eqxyA(Is))/MW/eqxyZ(Is))**(ONE/eqxyZ(Is))
         !1st
         eqxyGF(i,Is) = (eqxy(jDs,Is)/(eqxy(jMs,Is)*Dw*eqxyMol(i,Is))+1._dp)**(1._dp/3._dp)
         IF(lke)    eqxyKe(i,Is) = exp(4._dp*Mw*sigma/(R*eqTT(i)*Dw*eqxyGF(i,Is)*Dd))
         eqxyXS(i,Is) = ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(i,Is))+ONE)
         eqxyB   (Is) = eqxyXS(i,Is)**(ONE/(ONE+eqxyZ(Is)+eqxyXS(i,Is)))
         eqxyC   (Is) = ((eqxyKe(i,Is)/eqRH(i)-eqxyA(Is))/MW/eqxyZ(Is))**(ONE/eqxyZ(Is))
         IF(eqxyB(Is) < eqxyC(Is)) eqxyMol(i,Is) = eqxyC(Is)-eqxyB(Is)
         !2nd
         eqxyGF(i,Is) = (eqxy(jDs,Is)/(eqxy(jMs,Is)*Dw*eqxyMol(i,Is))+1._dp)**(1._dp/3._dp)
         IF(lke)    eqxyKe(i,Is) = exp(4._dp*Mw*sigma/(R*eqTT(i)*Dw*eqxyGF(i,Is)*Dd))
         eqxyXS(i,Is) = ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(i,Is))+ONE)
         eqxyB   (Is) = eqxyXS(i,Is)**(ONE/(ONE+eqxyZ(Is)+eqxyXS(i,Is)))
         eqxyC   (Is) = ((eqxyKe(i,Is)/eqRH(i)-eqxyA(Is))/MW/eqxyZ(Is))**(ONE/eqxyZ(Is))
         IF(eqxyB(Is) < eqxyC(Is)) eqxyMol(i,Is) = eqxyC(Is)-eqxyB(Is)
         !3rd
         eqxyGF(i,Is) = (eqxy(jDs,Is)/(eqxy(jMs,Is)*Dw*eqxyMol(i,Is))+1._dp)**(1._dp/3._dp)
         IF(lke)    eqxyKe(i,Is) = exp(4._dp*Mw*sigma/(R*eqTT(i)*Dw*eqxyGF(i,Is)*Dd))
         eqxyXS(i,Is) = ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(i,Is))+ONE)
         eqxyB   (Is) = eqxyXS(i,Is)**(ONE/(ONE+eqxyZ(Is)+eqxyXS(i,Is)))
         eqxyC   (Is) = ((eqxyKe(i,Is)/eqRH(i)-eqxyA(Is))/MW/eqxyZ(Is))**(ONE/eqxyZ(Is))
         IF(eqxyB(Is) < eqxyC(Is)) eqxyMol(i,Is) = eqxyC(Is)-eqxyB(Is)
         eqxyXS(i,Is) = ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(i,Is))+ONE)
      END DO
   END DO
!______________________________________
! RELATIVE HUMIDITY OF DELIQUESCENCE
   DO IJ=1,jSOL
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         Is=ieqsalt(i,IJ)
         IF(IWATeq(i) == 2) THEN
            IF(.NOT.leqsolute(Is)) CYCLE
            IF(eqxy(jRHD,Is) >= ONE) CYCLE
            eqxyRHD (i,Is)=eqxy(jRHD,Is)*EXP(eqxy(jRHDc,Is)*(ONE/eqTT(i)-ONE/eqT0))
            eqxyRHD (i,Is)=eqxyRHD(i,Is)*eqxyKe(i,Is)
            eqxyRHD (i,Is)=MAX(ZERO,MIN(eqxyRHD(i,Is),0.999_dp))
         ELSE
            eqxyRHD (i,Is)=eqRH(i)
         END IF
      END DO
   END DO
!______________________________________
! EQUILIBRIUM
   DO IJ=1,jSOL
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         Is=ieqsalt(i,IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         Ip=ieqsolute(Is,jSP)
         Im=ieqsolute(Is,jSM)
         IF(Ip == 0 .OR. Im == 0 .OR. IDeq(i)==5) CYCLE
         IF(eqxp(i,Ip)*eqxm(i,Im) > REALZERO) THEN
            NZeq(i) = NZeq(i) + 1
            XZ = MAX(ZERO, MIN(eqxpz(Ip)*eqxp(i,Ip),eqxmz(Im )*eqxm(i,Im)))
            eqxp(i,Ip)   = MAX(ZERO,eqxp(i,Ip)-XZ/eqxpz(Ip))
            eqxm(i,Im)   = MAX(ZERO,eqxm(i,Im)-XZ/eqxmz(Im))
            eqys(i,Is,jAP) = eqys(i,Is,jAP)+XZ/eqxy(jZa,Is)
            eqNs(i)=eqNs(i)+eqys(i,Is,jAP)
         ELSE
            eqxyRHD(i,Is) = ONE
         END IF
      END DO
   END DO
!______________________________________
   IF(lrhdm) THEN
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      DO IJ=1,jSOL
         Is=ieqsalt(i,IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         IF(eqys(i,Is,jAP) > TINYX .AND. IWATeq(i) == 2) THEN
            eqNZm (i)=eqNZm(i)+ONE
            eqMs  (i)=eqMs(i)+eqxy(jMs,Is)
            XZ       =(ONE/(100._dp/eqxy(jWs,Is)-ONE))/eqxy(jMs,Is)
            eqMolm(i)=eqMolm(i)+XZ
            eqWs  (i)=MAX(0.1_dp,MIN(ONE,ONE/(ONE/(eqMolm(i)*eqMs(i))+ONE)))
            YY    (i)=ONE/MAX(ZERO,(0.25_dp*log(eqWs(i))+ONE))
            eqRHDM(i)=ONE/(ONE+Mw*YY(i)*(eqMolm(i))**YY(i))
         END IF
      END DO
   END DO
   END IF
   YY(:)=ONE
!______________________________________
! GAS/AEROSOL PARTITIONING OF SEMI-VOLATILES
   IF(lvola) THEN
   ! semi-volatiles
   DO IJ=1,jVOL
      Is=ieqvola(IJ)
      IF(.NOT.leqsolute(Is)) CYCLE
      Ip=ieqsolute(Is,jSP)
      Im=ieqsolute(Is,jSM)
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         IF(IDeq(i) >2) CYCLE
         XZ = eqys(i,Is,jAP)
         IF(XZ > TINYX) THEN
            YY(i)=ONE
            ! TEMPERATURE DEPENDENT EQUILIBRIUM CONSTANTS
            T0T=eqT0/eqTT(i)
            COEF=ONE+LOG(T0T)-T0T
            X1=ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(i,Is))+ONE)
            X2=MAX(ZERO,MIN(2._dp*X1**2._dp,2._dp))
            XRHDMIN(i)=eqxyRHD(i,Is)
            IF(Is == eqacl) THEN
               ! NH4CL(S) <==> NH3(G) + HCL(G)   [ppb^2]
               X3   = 1.086E-16_dp ! [ppb^2]
               X3   = X3*EXP(-71.00_dp*(T0T-ONE)+2.400_dp*COEF)
               KEQ  = X3/(eqR*eqTT(i))/(eqR*eqTT(i)) ! [(mol^2/m^3(air))^2]
            END IF
            IF(Is == eqano) THEN
               ! NH4NO3(S) <==> NH3(G) + HNO3(G)
               ! ISORROPIA2
               X3   = 5.746E-17_dp ! [ppb^2]
               X3   = X3*EXP(-74.38_dp*(T0T-ONE)+6.120_dp*COEF)
               ! Mozurkewich (1993)
               !X3   = 4.199E-17_dp ! [ppb^2]
               !X3   = X3*EXP(-74.7351_dp*(T0T-ONE)+6.025_dp*COEF)
               ! SEQUILIB
               !X3   = 2.985e-17_dp ! [ppb^2]
               !X3   = X3*EXP(-75.11_dp*(T0T-ONE)+13.460_dp*COEF)
               KEQ   = X3/(eqR*eqTT(i))/(eqR*eqTT(i)) ! [(mol^2/m^3(air))^2]
            END IF
            COEF=X2
            IF(lmixs.AND.TSO4(i)>eqXPi(i)) THEN
!           IF(lmixs.AND.TSO4(i)>TINYX) THEN
               YY(i)=(XZ/(XZ+3._dp*TSO4(i)))**0.8_dp
               IF(lrhdm.AND.IWATeq(i)==2) &
               XRHDMIN(i)=eqRHDM(i)*YY(i)**0.25_dp+eqxyRHD(i,Is)   *(ONE-YY(i)**0.25_dp)
!              XRHDMIN(i)=eqRHDM(i)*YY(i)**0.25_dp+eqxyRHD(i,eqasu)*(ONE-YY(i)**0.25_dp)
               IF(eqRH(i)>=XRHDMIN(i)) THEN
                  COEF=X2*YY(i)
                  XZ=XZ*(ONE-COEF)
               END IF
            END IF
            IF(eqRH(i)<XRHDMIN(i)) THEN
               COEF=ONE
            ELSE
               IF(Is==eqacl) KEQ=KEQ*6.0_dp
            END IF
            KEQ=KEQ*COEF
            X1=eqxp(i,Ip)+eqxm(i,Im) ! [mol/m3(air)]
            X2=SQRT(X1*X1+4._dp*KEQ)
            X3=0.5_dp*(-X1+X2)
            X3=MIN(XZ,X3)
            eqxp(i,Ip)=eqxp(i,Ip)+X3
            eqxm(i,Im)=eqxm(i,Im)+X3
            eqys(i,Is,jAP)=MAX(0._dp,eqys(i,Is,jAP)-X3)
         END IF
      END DO
   END DO
   END IF
!______________________________________
! LIQUID/SOLID PARTITIONING
   DO IJ=1,jSOL
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         Is=ieqsalt(i,IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         IF(eqys(i,Is,jAP) < REALZERO) CYCLE
         IF(IWATeq(i) == 2) THEN
            XRHDIFF(i)=ONE
            XRHDMAX(i)=eqxyRHD(i,Is)
            IF(lrhdm.AND.eqNZm(i)>ONE.AND.Is/=eqcsu.AND.Is/=eqpsu) THEN
               XRHDMIN(i)=eqRHDM(i)
               YZ(i)=eqys(i,Is,jAP)/eqNs(i)
               XRHDMAX(i)=XRHDMIN(i)*YZ(i)**0.25_dp+XRHDMAX(i)*(ONE-YZ(i)**0.25_dp)
            ELSE
               XRHDMIN(i)=XRHDMAX(i)
            END IF
            IF(eqRH(i) < XRHDMAX(i)) THEN
               IF(eqRH(i)> XRHDMIN(i)   .AND.   XRHDMIN(i)<XRHDMAX(i)) &
               XRHDIFF(i)=(XRHDMAX(i)-eqRH(i))/(XRHDMAX(i)-XRHDMIN(i))
               eqys(i,Is,jDP) = MAX(0._dp,eqys(i,Is,jDP) + eqys(i,Is,jAP)*XRHDIFF(i))
               eqys(i,Is,jAP) = MAX(0._dp,eqys(i,Is,jAP) * (ONE-XRHDIFF(i)))
            END IF
         END IF
      END DO
   END DO
!______________________________________
! AEROSOL WATER [KG/M^3(AIR)]
   DO IJ=1,jSOL
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         Is=ieqsalt(i,IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         IF(eqxyMol(i,Is) > REALZERO) &
         eqWH2O(i) = eqWH2O(i) + eqys(i,Is,jAP)/eqxyMol(i,Is)
      END DO
   END DO
!______________________________________
! Remaining H+/OH- [MOL]
   DO IJ=1,mcati
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         eqXPi(i)=eqXPi(i)+eqxp(i,IJ)*eqxpz(IJ)
      END DO
   END DO
   DO IJ=1,manio
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         eqXMi(i)=eqXMi(i)+eqxm(i,IJ)*eqxmz(IJ)
      END DO
   END DO
!  DO i=neq1,neq2
!     IF(leqskip(i)) CYCLE
!     eqHPLUS(i)=eqXPi(i)-eqXMi(i)
!  END DO
!______________________________________
! RESIDUAL GASES
   DO IJ=1,jGAS
      Is=ieqgases(IJ)
      Ip=ieqsolute(Is,jSP)
      Im=ieqsolute(Is,jSM)
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         XZ=eqxm(i,Im)
         IF(Is == eqxam) XZ=eqxp(i,Ip)
         IF(Is == eqhsa) XZ=eqxm(i,maisu)+eqxm(i,maihs)
         IF(XZ < REALZERO) CYCLE
         IF(Is == eqhna) THEN
            eqyg(i,1)  = XZ
            eqxm(i,Im) = ZERO
         ELSE IF(Is == eqhca) THEN
            eqyg(i,2)  = XZ
            eqxm(i,Im) = ZERO
         ELSE IF(Is == eqxam) THEN
            eqyg(i,3)  = XZ
            eqxp(i,Ip) = ZERO
         ELSE IF(Is == eqhsa) THEN
            II=eqalc
            IF(eqxyMol(i,II) > REALZERO) &
            eqWH2O(i) = eqWH2O(i) + XZ/eqxyMol(i,II)
            IF(lH2SO4gas) THEN
               eqyg(i,4)  = XZ
               eqxm(i,maisu) = ZERO
               eqxm(i,maihs) = ZERO
            END IF
         END IF
      END DO
   END DO
!______________________________________
! OUTPUT
   DO IJ=1,mcati
      DO i=neq1,neq2
        IF(leqskip(i)) CYCLE
        eqyp (i,IJ,1)=eqyp(i,IJ,1)+eqxp(i,IJ)
        eqaPM(i)=eqaPM(i)+eqxp(i,IJ)
        eqPMt(i)=eqPMt(i)+eqxp(i,IJ)*eqxpm(IJ)
      END DO
   END DO
   DO IJ=1,manio
      DO i=neq1,neq2
        IF(leqskip(i)) CYCLE
        eqym (i,IJ,1)=eqym(i,IJ,1)+eqxm(i,IJ)
        eqaPM(i)=eqaPM(i)+eqxm(i,IJ)
        eqPMt(i)=eqPMt(i)+eqxm(i,IJ)*eqxmm(IJ)
      END DO
   END DO
   DO IJ=1,mcati
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         eqym(i,maisu,1) = eqym(i,maisu,1) + eqys(i,ieqsu(IJ),jAP)*eqxy(jZa,ieqsu(IJ))/eqxmz(ieqsolute(ieqsu(IJ),jSM))
         eqym(i,maisu,2) = eqym(i,maisu,2) + eqys(i,ieqsu(IJ),jDP)*eqxy(jZa,ieqsu(IJ))/eqxmz(ieqsolute(ieqsu(IJ),jSM))
         eqym(i,maihs,1) = eqym(i,maihs,1) + eqys(i,ieqhs(IJ),jAP)*eqxy(jZa,ieqhs(IJ))/eqxmz(ieqsolute(ieqhs(IJ),jSM))
         eqym(i,maihs,2) = eqym(i,maihs,2) + eqys(i,ieqhs(IJ),jDP)*eqxy(jZa,ieqhs(IJ))/eqxmz(ieqsolute(ieqhs(IJ),jSM))
         eqym(i,maino,1) = eqym(i,maino,1) + eqys(i,ieqno(IJ),jAP)*eqxy(jZa,ieqno(IJ))/eqxmz(ieqsolute(ieqno(IJ),jSM))
         eqym(i,maino,2) = eqym(i,maino,2) + eqys(i,ieqno(IJ),jDP)*eqxy(jZa,ieqno(IJ))/eqxmz(ieqsolute(ieqno(IJ),jSM))
         eqym(i,maicl,1) = eqym(i,maicl,1) + eqys(i,ieqcl(IJ),jAP)*eqxy(jZa,ieqcl(IJ))/eqxmz(ieqsolute(ieqcl(IJ),jSM))
         eqym(i,maicl,2) = eqym(i,maicl,2) + eqys(i,ieqcl(IJ),jDP)*eqxy(jZa,ieqcl(IJ))/eqxmz(ieqsolute(ieqcl(IJ),jSM))
      END DO
   END DO
   DO IJ=1,manio
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         eqyp(i,mciam,1) = eqyp(i,mciam,1) + eqys(i,ieqams(IJ),jAP)*eqxy(jZa,ieqams(IJ))/eqxpz(ieqsolute(ieqams(IJ),jSP))
         eqyp(i,mciam,2) = eqyp(i,mciam,2) + eqys(i,ieqams(IJ),jDP)*eqxy(jZa,ieqams(IJ))/eqxpz(ieqsolute(ieqams(IJ),jSP))
         eqyp(i,mciso,1) = eqyp(i,mciso,1) + eqys(i,ieqsos(IJ),jAP)*eqxy(jZa,ieqsos(IJ))/eqxpz(ieqsolute(ieqsos(IJ),jSP))
         eqyp(i,mciso,2) = eqyp(i,mciso,2) + eqys(i,ieqsos(IJ),jDP)*eqxy(jZa,ieqsos(IJ))/eqxpz(ieqsolute(ieqsos(IJ),jSP))
         eqyp(i,mcipo,1) = eqyp(i,mcipo,1) + eqys(i,ieqpos(IJ),jAP)*eqxy(jZa,ieqpos(IJ))/eqxpz(ieqsolute(ieqpos(IJ),jSP))
         eqyp(i,mcipo,2) = eqyp(i,mcipo,2) + eqys(i,ieqpos(IJ),jDP)*eqxy(jZa,ieqpos(IJ))/eqxpz(ieqsolute(ieqpos(IJ),jSP))
         eqyp(i,mcica,1) = eqyp(i,mcica,1) + eqys(i,ieqcas(IJ),jAP)*eqxy(jZa,ieqcas(IJ))/eqxpz(ieqsolute(ieqcas(IJ),jSP))
         eqyp(i,mcica,2) = eqyp(i,mcica,2) + eqys(i,ieqcas(IJ),jDP)*eqxy(jZa,ieqcas(IJ))/eqxpz(ieqsolute(ieqcas(IJ),jSP))
         eqyp(i,mcimg,1) = eqyp(i,mcimg,1) + eqys(i,ieqmgs(IJ),jAP)*eqxy(jZa,ieqmgs(IJ))/eqxpz(ieqsolute(ieqmgs(IJ),jSP))
         eqyp(i,mcimg,2) = eqyp(i,mcimg,2) + eqys(i,ieqmgs(IJ),jDP)*eqxy(jZa,ieqmgs(IJ))/eqxpz(ieqsolute(ieqmgs(IJ),jSP))
      END DO
   END DO
   ! PARTICULATE MATTER
   DO Is=1,jSOL
      IF(.NOT.leqsolute(Is)) CYCLE
      DO i=neq1,neq2
         IF(leqskip(i)) CYCLE
         eqsPM(i)=eqsPM(i)+eqys(i,Is,jDP)*eqxy(jNs,Is)
         eqaPM(i)=eqaPM(i)+eqys(i,Is,jAP)*eqxy(jNs,Is)
         eqPMt(i)=eqPMt(i)+eqys(i,Is,jAP)*eqxy(jMs,Is)
         eqPMs(i)=eqPMs(i)+eqys(i,Is,jDP)*eqxy(jMs,Is)
         eqVOL(i)=eqVOL(i)+eqys(i,Is,jAP)*eqxy(jMs,Is)/eqxy(jDs,Is)
         eqVOL(i)=eqVOL(i)+eqys(i,Is,jDP)*eqxy(jMs,Is)/eqxy(jDs,Is)
      END DO
   END DO
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      ! TOTAL PM      [KG/M^3(AIR)]
      eqPMt(i)=eqPMt(i)+eqPMs(i)
      ! TOTAL VOLUME  [M^3/M^3(AIR)]
      ! TOTAL DENSITY [KG/M^3]
      IF(eqVOL(i) > REALZERO) &
      eqRHO(i)=eqPMt(i)/eqVOL(i)
      ! AQUEOUS PHASE PROPERTIES
   END DO
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      eqGF(i) = ONE
      IF(eqWH2O(i) > 1.e-12_dp) THEN
         T0T=eqT0/eqTT(i)
         COEF=ONE+LOG(T0T)-T0T
         ! AUTODISSOCIATION CONSTANT (KW) OF WATER
         X1   = 1.010E-14_dp
         KEQ  = X1*EXP(-22.52_dp*(T0T-ONE) + 26.920_dp*COEF)
         ! H2O <==> H+ + OH- WITH KW [MOL^2/KG^2]
         AKW  = KEQ*eqRH(i)*eqWH2O(i)*eqWH2O(i)
         ! [OH-] = [H+] [MOL]
         AKW  = AKW**0.5_dp
         ! H+ PARAMETERIZATION
         XZ=eqXMi(i)-eqXPi(i)
         IF(IDeq(i) > 2) THEN
            eqHPLUS(i)=XZ+eqyg(i,3)-eqyg(i,1)
            XZ=(eqHPLUS(i)/eqWH2O(i)+AKW)*Dw*1.e-3_dp
         ELSE IF(IDeq(i) < 3) THEN
             eqHPLUS(i)=XZ+eqyg(i,3)
             XZ=(eqHPLUS(i)/eqWH2O(i)+AKW)*Dw*1.e-6_dp
         END IF
         ! AEROSOL PH
         IF (XZ > REALZERO) THEN
            ! HYDROGEN CONCENTRATION [MOL/L(H2O)]
            eqPH(i) = -LOG10(XZ)
         ELSE IF (XZ < REALZERO) THEN
            ! HYDROXY ION CONCENTRATION [MOL/L(H2O)]
            eqPH(i) = 14._dp + LOG10(-XZ)
            eqHPLUS(i)=ZERO
         END IF
         ! Growth Factor [-]
         eqGF(i) =(eqRHO(i)/Dw*eqWH2O(i)/eqPMt(i)+ONE)**(1._dp/3._dp)
      END IF
   END DO
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      ! AEROSOL WATER  [UG/M^3]
      eqWH2O(i) = eqWH2O(i)*1.E9_dp
      ! TOTAL PM      [UG/M^3(AIR)]
      eqPMt(i)=eqPMt(i)*1.E9_dp
      ! DRY PM        [UG/M^3(AIR)]
      eqPMs(i)=eqPMs(i)*1.E9_dp
      ! DRY PM        [UMOL/M^3(AIR)]
      eqsPM(i)=eqsPM(i)*1.E6_dp
      ! AQUEOUS PM    [UMOL/M^3(AIR)]
      eqaPM(i)=eqaPM(i)*1.E6_dp
      ! AQUEOUS H+    [MOL/M^3(AIR)]
      eqHPLUS(i)=eqHPLUS(i)!*1.E6_dp
      ! DENSITY [G/CM^3]
      eqRHO(i)=eqRHO(i)*1.e-3_dp
      eqVOL(i)=eqPMt(i)/eqRHO(i)
   END DO
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      xWH2O(i,n)=eqWH2O(i)
      xPMt (i,n)=eqPMt (i)
      xPMs (i,n)=eqPMs (i)
      xsPM (i,n)=eqsPM (i)
      xaPM (i,n)=eqaPM (i)
      xRHO (i,n)=eqRHO (i)
      xVOL (i,n)=eqVOL (i)
      xPH  (i,n)=eqPH  (i)
      xGF  (i,n)=eqGF  (i)
      xHp  (i,n)=eqHPLUS(i)
   END DO
   DO IJ=1,mcati
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
       xYPa(i,n,IJ)=eqyp(i,IJ,1)
       xYPs(i,n,IJ)=eqyp(i,IJ,2)
   END DO
   END DO
   DO IJ=1,manio
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
       xYMa(i,n,IJ)=eqym(i,IJ,1)
       xYMs(i,n,IJ)=eqym(i,IJ,2)
   END DO
   END DO
   DO IJ=1,jGAS
   DO i=neq1,neq2
      IF(leqskip(i)) CYCLE
      xYG(i,n,IJ)=eqyg(i,IJ)
   END DO
   END DO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Start OUTPUT for EMEP:
   DO i=neq1,neq2
!//.. aerosols [umol/m^3]
      aSO4out(i) = (eqym(i,maisu,jAP)+eqym(i,maisu,jDP))*1.E6_dp & ! particulate sulfate  (p=a+s)
                 + (eqym(i,maihs,jAP)+eqym(i,maihs,jDP))*1.E6_dp   ! particulate bisulfate(p=a+s)
      aNO3out(i) = (eqym(i,maino,jAP)+eqym(i,maino,jDP))*1.E6_dp   ! particulate nitrate  (p=a+s)
      aClout (i) = (eqym(i,maicl,jAP)+eqym(i,maicl,jDP))*1.E6_dp   ! particulate chloride (p=a+s)
      aNH4out(i) = (eqyp(i,mciam,jAP)+eqyp(i,mciam,jDP))*1.E6_dp   ! particulate ammonium (p=a+s)
      aNAout (i) = (eqyp(i,mciso,jAP)+eqyp(i,mciso,jDP))*1.E6_dp   ! particulate sodium   (p=a+s)
!//.. gases [umol/m^3]
      gNO3out(i) = eqyg(i,1)*1.E6_dp       ! residual HNO3  (g)
      gCLout (i) = eqyg(i,2)*1.E6_dp       ! residual HCl   (g)
      gNH3out(i) = eqyg(i,3)*1.E6_dp       ! residual NH3   (g)
      gSO4out(i) = eqyg(i,4)*1.E6_dp       ! residual H2SO4 (g)
!//.. aerosol water [ug/m^3]
      aH2Oout(i) = eqWH2O(i)               ! aerosol Water  (aq)
   END DO
!  END OUTPUT for EMEP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!______________________________________
1000 CONTINUE  ! leqskip_all = .TRUE.
!______________________________________
END DO ! nleq1,nleq2
!____________________________________________________________________________________________
END SUBROUTINE EQSAM4clim
!____________________________________________________________________________________________

end module EQSAM4clim_ml
