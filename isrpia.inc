!=======================================================================
! *** ISORROPIA CODE II
! *** INCLUDE FILE 'ISRPIA.INC'
! *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!     AND VARIABLES.
!
! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS
!
!=======================================================================
!
      IMPLICIT REAL (A-H,O-Z)
      PARAMETER (NCOMP=8,NIONS=10,NGASAQ=3,NSLDS=19,NPAIR=23,NZSR=100,&
                 NERRMX=25,NORG=3)
!
! *** INPUT VARIABLES **************************************************
!
      COMMON /INPT/ W(NCOMP), ORG(NORG), WAER(NCOMP), TEMP, RH, IPROB,&
                    NADJ
!
! *** WATER ACTIVITIES OF PURE SALT SOLUTIONS **************************
!
      COMMON /ZSR / AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR),&
                    AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR),&
                    AWSA(NZSR), AWLC(NZSR), AWCS(NZSR), AWCN(NZSR),&
                    AWCC(NZSR), AWPS(NZSR), AWPB(NZSR), AWPN(NZSR),&
                    AWPC(NZSR), AWMS(NZSR), AWMN(NZSR), AWMC(NZSR)
!
! *** DELIQUESCENCE RELATIVE HUMIDITIES ********************************
!
      INTEGER :: WFTYP
      COMMON /DRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3,&
                    DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL,&
                    DRCASO4,  DRCANO32, DRCACL2,  DRK2SO4,  DRKHSO4,&
                    DRKNO3,   DRKCL,    DRMGSO4,  DRMGNO32, DRMGCL2
!
      COMMON /MDRH/ DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2, &
                    DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2, &
                    DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3, &
                    DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8, &
                    DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13,&
                    WFTYP
!
      COMMON /MDRH2/ DRMO1,    DRMO2,    DRMO3,    DRML1,    DRML2,&
                     DRML3,    DRMM1,    DRMM2,    DRMP1,    DRMP2,&
                     DRMP3,    DRMP4,    DRMP5,    DRMV1

!
! *** VARIABLES FOR LIQUID AEROSOL PHASE *******************************
!

      REAL :: MOLAL, MOLALR, M0
      REAL :: IONIC
	INTEGER :: ACTMOD
      LOGICAL :: CALAOU, CALAIN, FRST, DRYF
      COMMON /IONS/ MOLAL(NIONS), MOLALR(NPAIR), GAMA(NPAIR), ZZ(NPAIR), &
                    Z(NIONS),     GAMOU(NPAIR),  GAMIN(NPAIR),M0(NPAIR), &
                    GASAQ(NGASAQ),  &
                    EPSACT,       COH,    CHNO3,  CHCL, &
                    WATER,        IONIC,  WATCMP(NPAIR+1), &
                    FRST,         CALAIN, CALAOU, DRYF, ACTMOD

!
! *** VARIABLES FOR SOLID AEROSOL PHASE ********************************
!
      COMMON /SALT/ CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4, &
                    CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC, CCASO4,&
                    CCANO32, CCACL2,  CK2SO4,  CKHSO4,  CKNO3, CKCL,&
                    CMGSO4,  CMGNO32, CMGCL2
!
! *** VARIABLES FOR GAS PHASE ******************************************
!
      COMMON /GAS / GNH3, GHNO3, GHCL
!
! *** EQUILIBRIUM CONSTANTS ********************************************
!
      COMMON /EQUK/ XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10,&
                    XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41,&
                    XK42, XK15, XK16, XK17, XK18, XK19, XK20, XK23,&
                    XK24, XK25
!     &              , XK26, XK27
!
! *** MOLECULAR WEIGHTS ************************************************
!
      REAL :: IMW
      COMMON /OTHR/ R, IMW(NIONS), WMW(NCOMP), SMW(NPAIR)
!
! *** SOLUTION/INFO VARIABLES ******************************************
!
      CHARACTER SCASE*15
      COMMON /CASE/ SULRATW, SULRAT, SODRAT, SO4RAT, CRNARAT, CRRAT, &
                    SCASE
!
      COMMON /SOLN/ EPS, MAXIT, NSWEEP, NDIV, ICLACT
!
! *** ERROR SYSTEM *****************************************************
!
      CHARACTER ERRMSG*40
      INTEGER :: ERRSTK, NOFER
      LOGICAL :: STKOFL
      COMMON /EROR/ STKOFL, NOFER, ERRSTK(NERRMX), ERRMSG(NERRMX)
!
! *** GENERIC VARIABLES **************************************************
!
      CHARACTER VERSION*15
      COMMON /CGEN/ GREAT, TINY, TINY2, ZERO, ONE, VERSION
!
! *** END OF INCLUDE FILE **********************************************
!
