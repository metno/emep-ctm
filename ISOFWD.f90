!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISRP1F
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF
!     AN AMMONIUM-SULFATE AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!     THE AMBIENT RELATIVE HUMIDITY.
!
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS
!
!=======================================================================
!
      SUBROUTINE ISRP1F (WI, RHI, TEMPI)
        INCLUDE 'isrpia.inc'
        DIMENSION WI(NCOMP)
  !
  ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
  !
        CALL INIT1 (WI, RHI, TEMPI)
  !
  ! *** CALCULATE SULFATE RATIO *******************************************
  !
        SULRAT = W(3)/W(2)
  !
  ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
  !
  ! *** SULFATE POOR
  !
        IF (2.0.LE.SULRAT) THEN
        DC   = W(3) - 2.001D0*W(2)  ! For numerical stability
        W(3) = W(3) + MAX(-DC, ZERO)
  !
           SCASE = 'A2'
           CALL CALCA2                 ! Only liquid (metastable)
  !
  ! *** SULFATE RICH (NO ACID)
  !
        ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN
  !
           SCASE = 'B4'
           CALL CALCB4                 ! Only liquid (metastable)

        CALL CALCNH3
  !
  ! *** SULFATE RICH (FREE ACID)
  !
        ELSEIF (SULRAT.LT.1.0) THEN
  !
           SCASE = 'C2'
           CALL CALCC2                 ! Only liquid (metastable)

        CALL CALCNH3
        ENDIF
  !
  ! *** RETURN POINT
  !
        RETURN
  !
  ! *** END OF SUBROUTINE ISRP1F *****************************************
  !
      END SUBROUTINE ISRP1F
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE ISRP2F
  ! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF
  !     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM.
  !     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
  !     THE AMBIENT RELATIVE HUMIDITY.
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE ISRP2F (WI, RHI, TEMPI)
        INCLUDE 'isrpia.inc'
        DIMENSION WI(NCOMP)
  !
  ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
  !
        CALL INIT2 (WI, RHI, TEMPI)
  !
  ! *** CALCULATE SULFATE RATIO *******************************************
  !
        SULRAT = W(3)/W(2)
  !
  ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
  !
  ! *** SULFATE POOR
  !
        IF (2.0.LE.SULRAT) THEN
  !
           SCASE = 'D3'
           CALL CALCD3                 ! Only liquid (metastable)

  !
  ! *** SULFATE RICH (NO ACID)
  !     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES,
  !     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
  !     SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
  !     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
  !
        ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN
  !
           SCASE = 'B4'
           CALL CALCB4                 ! Only liquid (metastable)
           SCASE = 'E4'
  !
        CALL CALCNA                 ! HNO3(g) DISSOLUTION
  !
  ! *** SULFATE RICH (FREE ACID)
  !     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES,
  !     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
  !     SUBROUTINE CALCC? IS CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
  !     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
  !
        ELSEIF (SULRAT.LT.1.0) THEN
  !
           SCASE = 'C2'
           CALL CALCC2                 ! Only liquid (metastable)
           SCASE = 'F2'
  !
        CALL CALCNA                 ! HNO3(g) DISSOLUTION
        ENDIF
  !
  ! *** RETURN POINT
  !
        RETURN
  !
  ! *** END OF SUBROUTINE ISRP2F *****************************************
  !
      END SUBROUTINE ISRP2F
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE ISRP3F
  ! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
  !     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
  !     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
  !     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE ISRP3F (WI, RHI, TEMPI)
        INCLUDE 'isrpia.inc'
        DIMENSION WI(NCOMP)
  !
  ! *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
  !
        WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
        WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
  !
  ! *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
  !
        IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
           WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
           WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
        ENDIF
  !
  ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
  !
        CALL ISOINIT3 (WI, RHI, TEMPI)
  !
  ! *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
  !
        REST = 2.D0*W(2) + W(4) + W(5)
        IF (W(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
           W(1) = (ONE-1D-6)*REST         ! Adjust Na amount
           CALL PUSHERR (0050, 'ISRP3F')  ! Warning error: Na adjusted
        ENDIF
  !
  ! *** CALCULATE SULFATE & SODIUM RATIOS *********************************
  !
        SULRAT = (W(1)+W(3))/W(2)
        SODRAT = W(1)/W(2)
  !
  ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

  ! *** SULFATE POOR ; SODIUM POOR
  !
        IF (2.0.LE.SULRAT .AND. SODRAT.LT.2.0) THEN
  !
           SCASE = 'G5'
           CALL CALCG5                 ! Only liquid (metastable)

  !
  ! *** SULFATE POOR ; SODIUM RICH
  !
        ELSE IF (SULRAT.GE.2.0 .AND. SODRAT.GE.2.0) THEN
  !
           SCASE = 'H6'
           CALL CALCH6                 ! Only liquid (metastable)
  !
  ! *** SULFATE RICH (NO ACID)
  !
        ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN
  !
           SCASE = 'I6'
           CALL CALCI6                 ! Only liquid (metastable)
  !
        CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3                !                NH3
  !
  ! *** SULFATE RICH (FREE ACID)
  !
        ELSEIF (SULRAT.LT.1.0) THEN
  !
           SCASE = 'J3'
           CALL CALCJ3                 ! Only liquid (metastable)

  !
        CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3                !                NH3
        ENDIF
  !
  ! *** RETURN POINT
  !
        RETURN
  !
  ! *** END OF SUBROUTINE ISRP3F *****************************************
  !
      END SUBROUTINE ISRP3F
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE ISRP4F
  ! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
  !     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM-CALCIUM-POTASSIUM-MAGNESIUM
  !     AEROSOL SYSTEM.
  !     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
  !     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE ISRP4F (WI, RHI, TEMPI)
        INCLUDE 'isrpia.inc'
        DIMENSION WI(NCOMP)
        REAL :: NAFRI, NO3FRI
  !
  ! *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
  !
  !      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
  !      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
  !
  ! *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
  !
  !      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
  !         WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
  !         WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
  !      ENDIF
  !
  ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
  !
        CALL INIT4 (WI, RHI, TEMPI)
  !
  ! *** CHECK IF TOO MUCH SODIUM+CRUSTALS ; ADJUST AND ISSUE ERROR MESSAGE
  !
        REST = 2.D0*W(2) + W(4) + W(5)
  !
        IF (W(1)+W(6)+W(7)+W(8).GT.REST) THEN
  !
        CCASO4I  = MIN (W(2),W(6))
        FRSO4I   = MAX (W(2) - CCASO4I, ZERO)
        CAFRI    = MAX (W(6) - CCASO4I, ZERO)
        CCANO32I = MIN (CAFRI, 0.5D0*W(4))
        CAFRI    = MAX (CAFRI - CCANO32I, ZERO)
        NO3FRI   = MAX (W(4) - 2.D0*CCANO32I, ZERO)
        CCACL2I  = MIN (CAFRI, 0.5D0*W(5))
        CLFRI    = MAX (W(5) - 2.D0*CCACL2I, ZERO)
        REST1    = 2.D0*FRSO4I + NO3FRI + CLFRI
  !
        CNA2SO4I = MIN (FRSO4I, 0.5D0*W(1))
        FRSO4I   = MAX (FRSO4I - CNA2SO4I, ZERO)
        NAFRI    = MAX (W(1) - 2.D0*CNA2SO4I, ZERO)
        CNACLI   = MIN (NAFRI, CLFRI)
        NAFRI    = MAX (NAFRI - CNACLI, ZERO)
        CLFRI    = MAX (CLFRI - CNACLI, ZERO)
        CNANO3I  = MIN (NAFRI, NO3FRI)
        NO3FR    = MAX (NO3FRI - CNANO3I, ZERO)
        REST2    = 2.D0*FRSO4I + NO3FRI + CLFRI
  !
        CMGSO4I  = MIN (FRSO4I, W(8))
        FRMGI    = MAX (W(8) - CMGSO4I, ZERO)
        FRSO4I   = MAX (FRSO4I - CMGSO4I, ZERO)
        CMGNO32I = MIN (FRMGI, 0.5D0*NO3FRI)
        FRMGI    = MAX (FRMGI - CMGNO32I, ZERO)
        NO3FRI   = MAX (NO3FRI - 2.D0*CMGNO32I, ZERO)
        CMGCL2I  = MIN (FRMGI, 0.5D0*CLFRI)
        CLFRI    = MAX (CLFRI - 2.D0*CMGCL2I, ZERO)
        REST3    = 2.D0*FRSO4I + NO3FRI + CLFRI
  !
           IF (W(6).GT.REST) THEN                       ! Ca > 2*SO4+CL+NO3 ?
               W(6) = (ONE-1D-6)*REST              ! Adjust Ca amount
               W(1)= ZERO                          ! Adjust Na amount
               W(7)= ZERO                          ! Adjust K amount
               W(8)= ZERO                          ! Adjust Mg amount
               CALL PUSHERR (0051, 'ISRP4F')       ! Warning error: Ca, Na, K, Mg in excess
  !
           ELSE IF (W(1).GT.REST1) THEN                 ! Na > 2*FRSO4+FRCL+FRNO3 ?
               W(1) = (ONE-1D-6)*REST1             ! Adjust Na amount
               W(7)= ZERO                          ! Adjust K amount
               W(8)= ZERO                          ! Adjust Mg amount
               CALL PUSHERR (0052, 'ISRP4F')       ! Warning error: Na, K, Mg in excess
  !
           ELSE IF (W(8).GT.REST2) THEN                 ! Mg > 2*FRSO4+FRCL+FRNO3 ?
               W(8) = (ONE-1D-6)*REST2             ! Adjust Mg amount
               W(7)= ZERO                          ! Adjust K amount
               CALL PUSHERR (0053, 'ISRP4F')       ! Warning error: K, Mg in excess
  !
           ELSE IF (W(7).GT.REST3) THEN                 ! K > 2*FRSO4+FRCL+FRNO3 ?
               W(7) = (ONE-1D-6)*REST3             ! Adjust K amount
               CALL PUSHERR (0054, 'ISRP4F')       ! Warning error: K in excess
           ENDIF
        ENDIF
  !
  ! *** CALCULATE RATIOS *************************************************
  !
        SO4RAT  = (W(1)+W(3)+W(6)+W(7)+W(8))/W(2)
        CRNARAT = (W(1)+W(6)+W(7)+W(8))/W(2)
        CRRAT   = (W(6)+W(7)+W(8))/W(2)
  !
  ! *** FIND CALCULATION REGIME FROM (SO4RAT, CRNARAT, CRRAT, RRH) ********
  !
  ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) POOR: R(Cr+Na)<2
  !
        IF (2.0.LE.SO4RAT .AND. CRNARAT.LT.2.0) THEN
  !
           SCASE = 'O7'
           CALL CALCO7                 ! Only liquid (metastable)

  !
  ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
  !
        ELSEIF (SO4RAT.GE.2.0 .AND. CRNARAT.GE.2.0) THEN
  !
         IF (CRRAT.LE.2.0) THEN
  !
           SCASE = 'M8'
           CALL CALCM8                 ! Only liquid (metastable)

  !        CALL CALCHCO3
  !
  ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
  !
          ELSEIF (CRRAT.GT.2.0) THEN
  !
           SCASE = 'P13'
           CALL CALCP13                 ! Only liquid (metastable)

  !        CALL CALCHCO3
         ENDIF
  !
  ! *** SULFATE RICH (NO ACID): 1<Rso4<2;
  !
        ELSEIF (1.0.LE.SO4RAT .AND. SO4RAT.LT.2.0) THEN
  !
           SCASE = 'L9'
           CALL CALCL9                ! Only liquid (metastable)

  !
        CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3                !                NH3
  !
  ! *** SULFATE SUPER RICH (FREE ACID): Rso4<1;
  !
        ELSEIF (SO4RAT.LT.1.0) THEN
  !
           SCASE = 'K4'
           CALL CALCK4                 ! Only liquid (metastable)

  !
        CALL CALCNHA                  ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3                  !                NH3
  !
        ENDIF
  !
        RETURN
        END
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCA2
  ! *** CASE A2
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT >= 2.0)
  !     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
  !
  !     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS X, THE
  !     AMOUNT OF HYDROGEN IONS (H+) FOUND IN THE LIQUID PHASE.
  !     FOR EACH ESTIMATION OF H+, FUNCTION FUNCB2A CALCULATES THE
  !     CONCENTRATION OF IONS FROM THE NH3(GAS) - NH4+(LIQ) EQUILIBRIUM.
  !     ELECTRONEUTRALITY IS USED AS THE OBJECTIVE FUNCTION.
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCA2
        INCLUDE 'isrpia.inc'
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU    =.TRUE.       ! Outer loop activity calculation flag
        OMELO     = TINY        ! Low  limit: SOLUTION IS VERY BASIC
        OMEHI     = 2.0D0*W(2)  ! High limit: FROM NH4+ -> NH3(g) + H+(aq)
  !
  ! *** CALCULATE WATER CONTENT *****************************************
  !
        MOLAL(5) = W(2)
        MOLAL(6) = ZERO
        CALL CALCMR
  !
  ! *** INITIAL VALUES FOR BISECTION ************************************
  !
        X1 = OMEHI
        Y1 = FUNCA2 (X1)
        IF (ABS(Y1).LE.EPS) RETURN
  !
  ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
  !
        DX = (OMEHI-OMELO)/FLOAT(NDIV)
        DO 10 I=1,NDIV
           X2 = MAX(X1-DX, OMELO)
           Y2 = FUNCA2 (X2)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
           X1 = X2
           Y1 = Y2
  10    CONTINUE
        IF (ABS(Y2).LE.EPS) THEN
           RETURN
        ELSE
           CALL PUSHERR (0001, 'CALCA2')    ! WARNING ERROR: NO SOLUTION
           RETURN
        ENDIF
  !
  ! *** PERFORM BISECTION ***********************************************
  !
  20    DO 30 I=1,MAXIT
           X3 = 0.5*(X1+X2)
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y3 = FUNCA2 (X3)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
              Y2    = Y3
              X2    = X3
           ELSE
              Y1    = Y3
              X1    = X3
           ENDIF
           IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
  30    CONTINUE
        CALL PUSHERR (0002, 'CALCA2')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CONVERGED ; RETURN **********************************************
  !
  40    X3 = 0.5*(X1+X2)
        CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
        Y3 = FUNCA2 (X3)
        RETURN
  !
  ! *** END OF SUBROUTINE CALCA2 ****************************************
  !
      END SUBROUTINE CALCA2

  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** FUNCTION FUNCA2
  ! *** CASE A2
  !     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE A2 ;
  !     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA2.
  !
  !=======================================================================
  !
        REAL FUNCTION FUNCA2 (OMEGI)
        INCLUDE 'isrpia.inc'
        REAL :: LAMDA
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        FRST   = .TRUE.
        CALAIN = .TRUE.
        PSI    = W(2)         ! INITIAL AMOUNT OF (NH4)2SO4 IN SOLUTION
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
           A1    = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
           A2    = XK2*R*TEMP/XKW*(GAMA(8)/GAMA(9))**2.
           A3    = XKW*RH*WATER*WATER
  !
           LAMDA = PSI/(A1/OMEGI+ONE)
           ZETA  = A3/OMEGI
  !
  ! *** SPECIATION & WATER CONTENT ***************************************
  !
           MOLAL (1) = OMEGI                                        ! HI
  !    ! slc.question
           MOLAL (5) = MAX(PSI-LAMDA,TINY)                          ! SO4I
           MOLAL (3) = MAX(W(3)/(ONE/A2/OMEGI + ONE), 2.*MOLAL(5))  ! NH4I
           MOLAL (6) = LAMDA                                        ! HSO4I
           GNH3      = MAX (W(3)-MOLAL(3), TINY)                    ! NH3GI
           COH       = ZETA                                         ! OHI
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
           IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
              CALL CALCACT
           ELSE
              GOTO 20
           ENDIF
  10    CONTINUE
  !
  ! *** CALCULATE OBJECTIVE FUNCTION ************************************
  !
  20    DENOM = (2.0*MOLAL(5)+MOLAL(6))
        FUNCA2= (MOLAL(3)/DENOM - ONE) + MOLAL(1)/DENOM
        RETURN
  !
  ! *** END OF FUNCTION FUNCA2 ********************************************
  !
      END FUNCTION FUNCA2
  !=======================================================================

  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCB4
  ! *** CASE B4
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
  !     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
  !
  !     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
  !     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
  !     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCB4
        INCLUDE 'isrpia.inc'
  !
  ! *** SOLVE EQUATIONS **************************************************
  !
        FRST       = .TRUE.
        CALAIN     = .TRUE.
        CALAOU     = .TRUE.
  !
  ! *** CALCULATE WATER CONTENT ******************************************
  !
        CALL CALCB1A         ! GET DRY SALT CONTENT, AND USE FOR WATER.
        MOLALR(13) = CLC
        MOLALR(9)  = CNH4HS4
        MOLALR(4)  = CNH42S4
        CLC        = ZERO
        CNH4HS4    = ZERO
        CNH42S4    = ZERO
        WATER      = MOLALR(13)/M0(13)+MOLALR(9)/M0(9)+MOLALR(4)/M0(4)
  !
        MOLAL(3)   = W(3)   ! NH4I
  !
        DO 20 I=1,NSWEEP
           AK1   = XK1*((GAMA(8)/GAMA(7))**2.)*(WATER/GAMA(7))
           BET   = W(2)
           GAM   = MOLAL(3)
  !
           BB    = BET + AK1 - GAM
           CC    =-AK1*BET
           DD    = BB*BB - 4.D0*CC
  !
  ! *** SPECIATION & WATER CONTENT ***************************************
  !
           MOLAL (5) = MAX(TINY,MIN(0.5*(-BB + SQRT(DD)), W(2))) ! SO4I
           MOLAL (6) = MAX(TINY,MIN(W(2)-MOLAL(5),W(2)))         ! HSO4I
           MOLAL (1) = MAX(TINY,MIN(AK1*MOLAL(6)/MOLAL(5),W(2))) ! HI
           CALL CALCMR                                           ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
           IF (.NOT.CALAIN) GOTO 30
           CALL CALCACT
  20    CONTINUE
  !
  30    RETURN
  !
  ! *** END OF SUBROUTINE CALCB4 ******************************************
  !
        END SUBROUTINE CALCB4
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCB1A
  ! *** CASE B1 ; SUBCASE 1
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE RICH
  !     2. THERE IS NO LIQUID PHASE
  !     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
  !                         BUT NOT BOTH)
  !
  !     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE AMOUNT OF LC
  !     IS CALCULATED FROM THE (NH4)2SO4 AND NH4HSO4 WHICH IS LEAST
  !     ABUNDANT (STOICHIMETRICALLY). THE REMAINING EXCESS OF SALT
  !     IS MIXED WITH THE LC.
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON
  ! UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
        SUBROUTINE CALCB1A
        INCLUDE 'isrpia.inc'
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        X = 2*W(2)-W(3)       ! Equivalent NH4HSO4
        Y = W(3)-W(2)         ! Equivalent (NH4)2SO4
  !
  ! *** CALCULATE COMPOSITION *******************************************
  !
        IF (X.LE.Y) THEN      ! LC is the MIN (x,y)
           CLC     = X        ! NH4HSO4 >= (NH4)2S04
           CNH4HS4 = ZERO
           CNH42S4 = Y-X
        ELSE
           CLC     = Y        ! NH4HSO4 <  (NH4)2S04
           CNH4HS4 = X-Y
           CNH42S4 = ZERO
        ENDIF
        RETURN
  !
  ! *** END OF SUBROUTINE CALCB1
  ! ******************************************
  !
      END SUBROUTINE CALCB1A

  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCC2
  ! *** CASE C2
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
  !     2. THERE IS ONLY A LIQUID PHASE
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCC2
        INCLUDE 'isrpia.inc'
        REAL :: LAMDA, KAPA
  !
        CALAOU =.TRUE.         ! Outer loop activity calculation flag
        FRST   =.TRUE.
        CALAIN =.TRUE.
  !
  ! *** SOLVE EQUATIONS **************************************************
  !
        LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
        PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION
        DO 20 I=1,NSWEEP
           PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.
           BB    = PSI+PARM
           CC    =-PARM*(LAMDA+PSI)
           KAPA  = 0.5*(-BB+SQRT(BB*BB-4.0*CC))
  !
  ! *** SPECIATION & WATER CONTENT ***************************************
  !
           MOLAL(1) = PSI+KAPA                               ! HI
           MOLAL(3) = LAMDA                                  ! NH4I
           MOLAL(5) = KAPA                                   ! SO4I
           MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
           CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3), ZERO)  ! Free H2SO4
           CALL CALCMR                                       ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
           IF (.NOT.CALAIN) GOTO 30
           CALL CALCACT
  20    CONTINUE
  !
  30    RETURN
  !
  ! *** END OF SUBROUTINE CALCC2 *****************************************
  !
        END SUBROUTINE CALCC2
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCD3
  ! *** CASE D3
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0)
  !     2. THERE IS OLNY A LIQUID PHASE
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCD3
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** FIND DRY COMPOSITION **********************************************
  !
        CALL CALCD1A
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CHI1 = CNH4NO3               ! Save from CALCD1 run
        CHI2 = CNH42S4
        CHI3 = GHNO3
        CHI4 = GNH3
  !
        PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
        PSI2 = CHI2
        PSI3 = ZERO
        PSI4 = ZERO
  !
        MOLAL(5) = PSI2              ! Include initial amount in water calc
        MOLAL(6) = ZERO
        MOLAL(3) = PSI1
        MOLAL(7) = PSI1
        CALL CALCMR                  ! Initial water
  !
        CALAOU = .TRUE.              ! Outer loop activity calculation flag
        PSI4LO = TINY                ! Low  limit
        PSI4HI = CHI4                ! High limit
  !
  ! *** INITIAL VALUES FOR BISECTION ************************************
  !
  60    X1 = PSI4LO
        CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
        Y1 = FUNCD3 (X1)
        IF (ABS(Y1).LE.EPS) RETURN
        YLO= Y1                 ! Save Y-value at HI position
  !
  ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
  !
        DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
        DO 10 I=1,NDIV
           X2 = X1+DX
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y2 = FUNCD3 (X2)
           IF (((Y1) .LT. ZERO) .AND. ((Y2) .GT. ZERO)) GOTO 20  ! (Y1*Y2.LT.ZERO) (slc.1.2012)
  !         IF (SIGN(1.0,Y1)*SIGN(1.0,Y2).LT.ZERO) GOTO 20
           X1 = X2
           Y1 = Y2
  10    CONTINUE
  !
  ! *** NO SUBDIVISION WITH SOLUTION FOUND
  !
        YHI= Y1                      ! Save Y-value at Hi position
        IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION
           RETURN
  !
  ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
  ! Physically I dont know when this might happen, but I have put this
  ! branch in for completeness. I assume there is no solution; all NO3 goes to the
  ! gas phase.
  !
        ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
           P4 = TINY ! PSI4LO ! CHI4
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           YY = FUNCD3(P4)
           GOTO 50
  !
  ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
  ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
  ! and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
  ! and proceed again with root tracking.
  !
        ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
           PSI4HI = PSI4LO
           PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) ! No solution; some NH3 evaporates
           IF (PSI4LO.LT.-(PSI1+PSI2)) THEN
              CALL PUSHERR (0001, 'CALCD3')  ! WARNING ERROR: NO SOLUTION
              RETURN
           ELSE
              MOLAL(5) = PSI2              ! Include sulfate in initial water calculation
              MOLAL(6) = ZERO
              MOLAL(3) = PSI1
              MOLAL(7) = PSI1
              CALL CALCMR                  ! Initial water
              GOTO 60                        ! Redo root tracking
           ENDIF
        ENDIF
  !
  ! *** PERFORM BISECTION ***********************************************
  !
  20    DO 30 I=1,MAXIT
           X3 = 0.5*(X1+X2)
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y3 = FUNCD3 (X3)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
              Y2    = Y3
              X2    = X3
           ELSE
              Y1    = Y3
              X1    = X3
           ENDIF
           IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
  30    CONTINUE
        CALL PUSHERR (0002, 'CALCD3')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CONVERGED ; RETURN **********************************************
  !
  40    X3 = 0.5*(X1+X2)
        CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
        Y3 = FUNCD3 (X3)
  !
  ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
  !
  50    CONTINUE
        IF (MOLAL(1).GT.TINY) THEN
           CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
           MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
           MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
           MOLAL(6) = DELTA                                ! HSO4 EFFECT
        ENDIF
        RETURN
  !
  ! *** END OF SUBROUTINE CALCD3 ******************************************
  !
      END SUBROUTINE CALCD3

  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** FUNCTION FUNCD3
  ! *** CASE D3
  !     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ;
  !     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.
  !
  !=======================================================================
  !
        REAL FUNCTION FUNCD3 (P4)
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        FRST   = .TRUE.
        CALAIN = .TRUE.
        PSI4   = P4
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
           A2   = XK7*(WATER/GAMA(4))**3.0
           A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.0
           A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
           A7   = XKW *RH*WATER*WATER
  !
           PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
           PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4)
           PSI3 = MIN(MAX(PSI3, ZERO), CHI3)
  !
           BB   = PSI4 - PSI3
  !CCOLD         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
  !CC         AHI  =2.0*A7/(BB+SQRT(BB*BB + 4.d0*A7)) ! Avoid overflow when HI->0
           DENM = BB+SQRT(BB*BB + 4.d0*A7)
           IF (DENM.LE.TINY) THEN       ! Avoid overflow when HI->0
              ABB  = ABS(BB)
              DENM = (BB+ABB) + 2.0*A7/ABB ! Taylor expansion of SQRT
           ENDIF
           AHI = 2.0*A7/DENM
  !
  ! *** SPECIATION & WATER CONTENT ***************************************
  !
           MOLAL (1) = AHI                             ! HI
           MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2         ! NH4I
           MOLAL (5) = PSI2                            ! SO4I
           MOLAL (6) = ZERO                            ! HSO4I
           MOLAL (7) = PSI3 + PSI1                     ! NO3I
           CNH42S4   = CHI2 - PSI2                     ! Solid (NH4)2SO4
           CNH4NO3   = ZERO                            ! Solid NH4NO3
           GHNO3     = CHI3 - PSI3                     ! Gas HNO3
           GNH3      = CHI4 - PSI4                     ! Gas NH3
           CALL CALCMR                                 ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
           IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
              CALL CALCACT
           ELSE
              GOTO 20
           ENDIF
  10    CONTINUE
  !
  ! *** CALCULATE OBJECTIVE FUNCTION ************************************
  !
  20    CONTINUE
  !CC      FUNCD3= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE
        FUNCD3= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE
        RETURN
  !
  ! *** END OF FUNCTION FUNCD3 ********************************************
  !
      END FUNCTION FUNCD3
  !========================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCD1A
  ! *** CASE D1 ; SUBCASE 1
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0)
  !     2. SOLID AEROSOL ONLY
  !     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
  !
  !     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
  !     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
  !     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
  !     THE SOLID PHASE.
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCD1A
        INCLUDE 'isrpia.inc'
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        PARM    = XK10/(R*TEMP)/(R*TEMP)
  !
  ! *** CALCULATE NH4NO3 THAT VOLATIZES *********************************
  !
        CNH42S4 = W(2)
        X       = MAX(ZERO, MIN(W(3)-2.0*CNH42S4, W(4)))  ! MAX NH4NO3
        PS      = MAX(W(3) - X - 2.0*CNH42S4, ZERO)
        OM      = MAX(W(4) - X, ZERO)
  !
        OMPS    = OM+PS
        DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              ! DIAKRINOUSA
        ZE      = MIN(X, 0.5*(-OMPS + DIAK))              ! THETIKI RIZA
  !
  ! *** SPECIATION *******************************************************
  !
        CNH4NO3 = X  - ZE    ! Solid NH4NO3
        GNH3    = PS + ZE    ! Gas NH3
        GHNO3   = OM + ZE    ! Gas HNO3
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCD1A *****************************************
  !
        END SUBROUTINE CALCD1A
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCG5
  ! *** CASE G5
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCG5
        INCLUDE 'isrpia.inc'
  !
        REAL :: LAMDA
        COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
                       PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,  &
                       A1,   A2,   A3,   A4,   A5,   A6,   A7
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU = .TRUE.
        CHI1   = 0.5*W(1)
        CHI2   = MAX (W(2)-CHI1, ZERO)
        CHI3   = ZERO
        CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
        CHI5   = W(4)
        CHI6   = W(5)
  !
        PSI1   = CHI1
        PSI2   = CHI2
        PSI6LO = TINY
        PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
  !
        WATER  = CHI2/M0(4) + CHI1/M0(2)
  !
  ! *** INITIAL VALUES FOR BISECTION ************************************
  !
        X1 = PSI6LO
        Y1 = FUNCG5A (X1)
        IF (CHI6.LE.TINY) GOTO 50
  !cc      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
  !cc      IF (WATER .LE. TINY) RETURN                    ! No water
  !
  ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
  !
        DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
        DO 10 I=1,NDIV
           X2 = X1+DX
           Y2 = FUNCG5A (X2)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
           X1 = X2
           Y1 = Y2
  10    CONTINUE
  !
  ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
  !
        IF (ABS(Y2) .GT. EPS) THEN
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y2 = FUNCG5A (PSI6LO)
        ENDIF
        GOTO 50
  !
  ! *** PERFORM BISECTION ***********************************************
  !
  20    DO 30 I=1,MAXIT
           X3 = 0.5*(X1+X2)
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y3 = FUNCG5A (X3)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
              Y2    = Y3
              X2    = X3
           ELSE
              Y1    = Y3
              X1    = X3
           ENDIF
           IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
  30    CONTINUE
        CALL PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CONVERGED ; RETURN **********************************************
  !
  40    X3 = 0.5*(X1+X2)
        CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
        Y3 = FUNCG5A (X3)
  !
  ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
  !
  50    CONTINUE
        IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  ! If quadrat.called
           CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
           MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
           MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
           MOLAL(6) = DELTA                               ! HSO4 EFFECT
        ENDIF
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCG5 *******************************************
  !
      END SUBROUTINE CALCG5

  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE FUNCG5A
  ! *** CASE G5
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        REAL FUNCTION FUNCG5A (X)
        INCLUDE 'isrpia.inc'
  !
        REAL :: LAMDA
        COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
                       PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,  &
                       A1,   A2,   A3,   A4,   A5,   A6,   A7
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        PSI6   = X
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A2  = XK7 *(WATER/GAMA(4))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        AKK = A4*A6
  !
  !  CALCULATE DISSOCIATION QUANTITIES
  !
        IF (CHI5.GE.TINY) THEN
           PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
        ELSE
           PSI5 = TINY
        ENDIF
  !
  !CC      IF(CHI4.GT.TINY) THEN
        IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
           BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
           CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
           DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
           PSI4 =0.5d0*(-BB - SQRT(DD))
        ELSE
           PSI4 = TINY
        ENDIF
  !
  ! *** CALCULATE SPECIATION ********************************************
  !
        MOLAL (2) = 2.0D0*PSI1                          ! NAI
        MOLAL (3) = 2.0*PSI2 + PSI4                     ! NH4I
        MOLAL (4) = PSI6                                ! CLI
        MOLAL (5) = PSI2 + PSI1                         ! SO4I
        MOLAL (6) = ZERO
        MOLAL (7) = PSI5                                ! NO3I
  !
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
  !
        GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
        GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
        GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl
  !
        CNH42S4   = ZERO                                ! Solid (NH4)2SO4
        CNH4NO3   = ZERO                                ! Solid NH4NO3
        CNH4CL    = ZERO                                ! Solid NH4Cl
  !
        CALL CALCMR                                     ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !
  ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
  !
  20    FUNCG5A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
  !CC         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
  !
        RETURN
  !
  ! *** END OF FUNCTION FUNCG5A *******************************************
  !
      END FUNCTION FUNCG5A

  !========================================================================
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCH6
  ! *** CASE H6
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCH6
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU = .TRUE.
        CHI1   = W(2)                                ! CNA2SO4
        CHI2   = ZERO                                ! CNH42S4
        CHI3   = ZERO                                ! CNH4CL
        FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
        CHI8   = MIN (FRNA, W(4))                    ! CNANO3
        CHI4   = W(3)                                ! NH3(g)
        CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
        CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
        CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
  !
        PSI6LO = TINY
        PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
  !
  ! *** INITIAL VALUES FOR BISECTION ************************************
  !
        X1 = PSI6LO
        Y1 = FUNCH6A (X1)
        IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
  !
  ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
  !
        DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
        DO 10 I=1,NDIV
           X2 = X1+DX
           Y2 = FUNCH6A (X2)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
           X1 = X2
           Y1 = Y2
  10    CONTINUE
  !
  ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
  !
        IF (ABS(Y2) .GT. EPS) THEN
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y2 = FUNCH6A (PSI6LO)
        ENDIF
        GOTO 50
  !
  ! *** PERFORM BISECTION ***********************************************
  !
  20    DO 30 I=1,MAXIT
           X3 = 0.5*(X1+X2)
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y3 = FUNCH6A (X3)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
              Y2    = Y3
              X2    = X3
           ELSE
              Y1    = Y3
              X1    = X3
           ENDIF
           IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
  30    CONTINUE
        CALL PUSHERR (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CONVERGED ; RETURN **********************************************
  !
  40    X3 = 0.5*(X1+X2)
        CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
        Y3 = FUNCH6A (X3)
  !
  ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
  !
  50    CONTINUE
        IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
           CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
           MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
           MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
           MOLAL(6) = DELTA                                ! HSO4 EFFECT
        ENDIF
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCH6 ******************************************
  !
      END SUBROUTINE CALCH6

  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE FUNCH6A
  ! *** CASE H6
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        REAL FUNCTION FUNCH6A (X)
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        PSI6   = X
        PSI1   = CHI1
        PSI2   = ZERO
        PSI3   = ZERO
        PSI7   = CHI7
        PSI8   = CHI8
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
        A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
  !
  !  CALCULATE DISSOCIATION QUANTITIES
  !
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MAX(PSI5, TINY)
  !
        IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
           BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
           CC   = CHI4*(PSI5+PSI6)
           DD   = BB*BB-4.d0*CC
           PSI4 =0.5d0*(-BB - SQRT(DD))
           PSI4 = MIN(PSI4,CHI4)
        ELSE
           PSI4 = TINY
        ENDIF
  !
  ! *** CALCULATE SPECIATION ********************************************
  !
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1                           ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
  !
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
  !
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
  !
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = MAX(CHI7 - PSI7, ZERO)
        CNANO3    = MAX(CHI8 - PSI8, ZERO)
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
  !
        CALL CALCMR                                    ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !
  ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
  !
  20    FUNCH6A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
  !
        RETURN
  !
  ! *** END OF FUNCTION FUNCH6A *******************************************
  !
      END FUNCTION FUNCH6A

  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCI6
  ! *** CASE I6
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
  !     2. SOLID & LIQUID AEROSOL POSSIBLE
  !     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCI6
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,    &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,    &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,  &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,      &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,&
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** FIND DRY COMPOSITION **********************************************
  !
        CALL CALCI1A
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CHI1 = CNH4HS4               ! Save from CALCI1 run
        CHI2 = CLC
        CHI3 = CNAHSO4
        CHI4 = CNA2SO4
        CHI5 = CNH42S4
  !
        PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
        PSI2 = CLC
        PSI3 = CNAHSO4
        PSI4 = CNA2SO4
        PSI5 = CNH42S4
  !
        CALAOU = .TRUE.              ! Outer loop activity calculation flag
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
  !
  !  CALCULATE DISSOCIATION QUANTITIES
  !
        BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
        CC   =-A6*(PSI2 + PSI3 + PSI1)
        DD   = BB*BB - 4.D0*CC
        PSI6 = 0.5D0*(-BB + SQRT(DD))
  !
  ! *** CALCULATE SPECIATION ********************************************
  !
        MOLAL (1) = PSI6                                    ! HI
        MOLAL (2) = 2.D0*PSI4 + PSI3                        ! NAI
        MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1            ! NH4I
        MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6               ! SO4I
        MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6               ! HSO4I
        CLC       = ZERO
        CNAHSO4   = ZERO
        CNA2SO4   = CHI4 - PSI4
        CNH42S4   = ZERO
        CNH4HS4   = ZERO
        CALL CALCMR                                         ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !
  20    RETURN
  !
  ! *** END OF SUBROUTINE CALCI6 *****************************************
  !
      END SUBROUTINE CALCI6
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCI1A
  ! *** CASE I1 ; SUBCASE 1
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
  !     2. SOLID AEROSOL ONLY
  !     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCI1A
        INCLUDE 'isrpia.inc'
  !
  ! *** CALCULATE NON VOLATILE SOLIDS ***********************************
  !
        CNA2SO4 = 0.5D0*W(1)
        CNH4HS4 = ZERO
        CNAHSO4 = ZERO
        CNH42S4 = ZERO
        FRSO4   = MAX(W(2)-CNA2SO4, ZERO)
  !
        CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
        FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
        FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)
  !
        IF (FRSO4.LE.TINY) THEN
           CLC     = MAX(CLC - FRNH4, ZERO)
           CNH42S4 = 2.D0*FRNH4

        ELSEIF (FRNH4.LE.TINY) THEN
           CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
           CLC     = MAX(CLC-FRSO4, ZERO)
           IF (CNA2SO4.GT.TINY) THEN
              FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
              CNAHSO4 = 2.D0*FRSO4
              CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
           ENDIF
        ENDIF
  !
  ! *** CALCULATE GAS SPECIES *********************************************
  !
        GHNO3 = W(4)
        GHCL  = W(5)
        GNH3  = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCI1A *****************************************
  !
        END SUBROUTINE CALCI1A
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCJ3
  ! *** CASE J3
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
  !     2. THERE IS ONLY A LIQUID PHASE
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** UPDATED BY CHRISTOS FOUNTOUKIS
  !
  !=======================================================================
  !
        SUBROUTINE CALCJ3
        INCLUDE 'isrpia.inc'
  !
        REAL :: LAMDA, KAPA
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU = .TRUE.              ! Outer loop activity calculation flag
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
        LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
        CHI1   = W(1)                           ! NA TOTAL as NaHSO4
        CHI2   = W(3)                           ! NH4 TOTAL as NH4HSO4
        PSI1   = CHI1
        PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
  !
  !  CALCULATE DISSOCIATION QUANTITIES
  !
        BB   = A3+LAMDA                        ! KAPA
        CC   =-A3*(LAMDA + PSI1 + PSI2)
        DD   = BB*BB-4.D0*CC
        KAPA = 0.5D0*(-BB+SQRT(DD))
  !
  ! *** CALCULATE SPECIATION ********************************************
  !
        MOLAL (1) = LAMDA + KAPA                 ! HI
        MOLAL (2) = PSI1                         ! NAI
        MOLAL (3) = PSI2                         ! NH4I
        MOLAL (4) = ZERO                         ! CLI
        MOLAL (5) = KAPA                         ! SO4I
        MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA   ! HSO4I
        MOLAL (7) = ZERO                         ! NO3I
  !
        CNAHSO4   = ZERO
        CNH4HS4   = ZERO
  !
        CALL CALCMR                              ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 50
        ENDIF
  10    CONTINUE
  !
  50    RETURN
  !
  ! *** END OF SUBROUTINE CALCJ3 ******************************************
  !
        END SUBROUTINE CALCJ3
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCO7
  ! *** CASE O7
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : CaSO4
  !     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MgSO4, NA2SO4, K2SO4
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCO7
        INCLUDE 'isrpia.inc'
  !
        REAL :: LAMDA
        COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
                       CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,      &
                       PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,     &
                       A5, A6, A7, A8, A9
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU = .TRUE.
        CHI9   = MIN (W(6), W(2))                     ! CCASO4
        SO4FR  = MAX (W(2)-CHI9, ZERO)
        CAFR   = MAX (W(6)-CHI9, ZERO)
        CHI7   = MIN (0.5D0*W(7), SO4FR)              ! CK2SO4
        FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
        SO4FR  = MAX (SO4FR - CHI7, ZERO)
        CHI1   = MIN (0.5D0*W(1), SO4FR)              ! NA2SO4
        NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
        SO4FR  = MAX (SO4FR - CHI1, ZERO)
        CHI8   = MIN (W(8), SO4FR)                    ! CMGSO4
        FRMG    = MAX(W(8) - CHI8, ZERO)
        SO4FR   = MAX(SO4FR - CHI8, ZERO)
        CHI3   = ZERO
        CHI5   = W(4)
        CHI6   = W(5)
        CHI2   = MAX (SO4FR, ZERO)
        CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
  !
        PSI1   = CHI1
        PSI2   = CHI2
        PSI3   = ZERO
        PSI4   = ZERO
        PSI5   = ZERO
        PSI6   = ZERO
        PSI7   = CHI7
        PSI8   = CHI8
        PSI6LO = TINY
        PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
  !
        WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
        WATER  = MAX (WATER , TINY)
  !
  ! *** INITIAL VALUES FOR BISECTION ************************************
  !
        X1 = PSI6LO
        Y1 = FUNCO7 (X1)
        IF (CHI6.LE.TINY) GOTO 50
  !cc      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
  !cc      IF (WATER .LE. TINY) RETURN                    ! No water
  !
  ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
  !
        DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
        DO 10 I=1,NDIV
           X2 = X1+DX
           Y2 = FUNCO7 (X2)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
           X1 = X2
           Y1 = Y2
  10    CONTINUE
  !
  ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
  !
        IF (ABS(Y2) .GT. EPS) THEN
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y2 = FUNCO7 (PSI6LO)
        ENDIF
        GOTO 50
  !
  ! *** PERFORM BISECTION ***********************************************
  !
  20    DO 30 I=1,MAXIT
           X3 = 0.5*(X1+X2)
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y3 = FUNCO7 (X3)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
              Y2    = Y3
              X2    = X3
           ELSE
              Y1    = Y3
              X1    = X3
           ENDIF
           IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
  30    CONTINUE
        CALL PUSHERR (0002, 'CALCO7')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CONVERGED ; RETURN **********************************************
  !
  40    X3 = 0.5*(X1+X2)
        CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
        Y3 = FUNCO7 (X3)
  !
  ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
  !
  50    CONTINUE
        IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  ! If quadrat.called
           CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
           MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
           MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
           MOLAL(6) = DELTA                               ! HSO4 EFFECT
        ENDIF
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCO7 *******************************************
  !
        END SUBROUTINE CALCO7
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE FUNCO7
  ! *** CASE O7
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : CaSO4
  !     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MgSO4, NA2SO4, K2SO4
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        REAL FUNCTION FUNCO7 (X)
        INCLUDE 'isrpia.inc'
  !
        REAL :: LAMDA
        COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
                       CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,      &
                       PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,     &
                       A5, A6, A7, A8, A9
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        PSI6   = X
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
  !
  !
        IF (CHI5.GE.TINY) THEN
           PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
           PSI5 = MIN (PSI5,CHI5)
        ELSE
           PSI5 = TINY
        ENDIF
  !
  !CC      IF(CHI4.GT.TINY) THEN
        IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
           BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
           CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
           DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
           PSI4 =0.5d0*(-BB - SQRT(DD))
           PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
        ELSE
           PSI4 = TINY
        ENDIF
  !
  ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
  !
        MOLAL (2) = 2.0D0*PSI1                       ! Na+
        MOLAL (3) = 2.0D0*PSI2 + PSI4                ! NH4I
        MOLAL (4) = PSI6                             ! CLI
        MOLAL (5) = PSI1+PSI2+PSI7+PSI8              ! SO4I
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI5                             ! NO3I
        MOLAL (8) = ZERO                             ! CaI
        MOLAL (9) = 2.0D0*PSI7                       ! KI
        MOLAL (10)= PSI8                             ! Mg
  !
  ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
  !
  !CC      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
         SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                     -MOLAL(9)-2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
  !
  ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
  !
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
  !
        CNA2SO4  = ZERO
        CNH42S4  = ZERO
        CNH4NO3  = ZERO
        CNH4Cl   = ZERO
        CK2SO4   = ZERO
        CMGSO4   = ZERO
        CCASO4   = CHI9
  !
  ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
  !
        CALL CALCMR
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !
  ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
  !
  20    FUNCO7 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
  !CC         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
  !
        RETURN
  !
  ! *** END OF FUNCTION FUNCO7 *******************************************
  !
        END FUNCTION FUNCO7
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCM8
  ! *** CASE M8
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : CaSO4
  !     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4, K2SO4
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCM8
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU = .TRUE.
        CHI11  = MIN (W(6), W(2))                    ! CCASO4
        SO4FR  = MAX(W(2)-CHI11, ZERO)
        CAFR   = MAX(W(6)-CHI11, ZERO)
        CHI9   = MIN (0.5D0*W(7), SO4FR)             ! CK2S04
        FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
        SO4FR  = MAX(SO4FR-CHI9, ZERO)
        CHI10  = MIN (W(8), SO4FR)                  ! CMGSO4
        FRMG   = MAX(W(8)-CHI10, ZERO)
        SO4FR  = MAX(SO4FR-CHI10, ZERO)
        CHI1   = MAX (SO4FR,ZERO)                    ! CNA2SO4
        CHI2   = ZERO                                ! CNH42S4
        CHI3   = ZERO                                ! CNH4CL
        FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
        CHI8   = MIN (FRNA, W(4))                    ! CNANO3
        CHI4   = W(3)                                ! NH3(g)
        CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
        CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
        CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
  !
        PSI6LO = TINY
        PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
  !
  ! *** INITIAL VALUES FOR BISECTION ************************************
  !
        X1 = PSI6LO
        Y1 = FUNCM8 (X1)
        IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
  !
  ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
  !
        DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
        DO 10 I=1,NDIV
           X2 = X1+DX
           Y2 = FUNCM8 (X2)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
           X1 = X2
           Y1 = Y2
  10    CONTINUE
  !
  ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
  !
        IF (ABS(Y2) .GT. EPS) THEN
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y2 = FUNCM8 (PSI6LO)
        ENDIF
        GOTO 50
  !
  ! *** PERFORM BISECTION ***********************************************
  !
  20    DO 30 I=1,MAXIT
           X3 = 0.5*(X1+X2)
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y3 = FUNCM8 (X3)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
              Y2    = Y3
              X2    = X3
           ELSE
              Y1    = Y3
              X1    = X3
           ENDIF
           IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
  30    CONTINUE
        CALL PUSHERR (0002, 'CALCM8')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CONVERGED ; RETURN **********************************************
  !
  40    X3 = 0.5*(X1+X2)
        CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
        Y3 = FUNCM8 (X3)
  !
  ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
  !
  50    CONTINUE
        IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
           CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
           MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
           MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
           MOLAL(6) = DELTA                                ! HSO4 EFFECT
        ENDIF
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCM8 ******************************************
  !
        END SUBROUTINE CALCM8

  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE FUNCM8
  ! *** CASE M8
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : CaSO4
  !     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4, K2SO4
  !
  ! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
  !
  !=======================================================================
  !
        REAL FUNCTION FUNCM8 (X)
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,    &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,    &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,  &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,      &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,&
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        PSI6   = X
        PSI1   = CHI1
        PSI2   = ZERO
        PSI3   = ZERO
        PSI7   = CHI7
        PSI8   = CHI8
        PSI9   = CHI9
        PSI10  = CHI10
        PSI11  = ZERO
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
  !      A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
  !      A7  = XK8 *(WATER/GAMA(1))**2.0
  !      A8  = XK9 *(WATER/GAMA(3))**2.0
  !      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
  !
  !  CALCULATE DISSOCIATION QUANTITIES
  !
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
  !
        IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
           BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
           CC   = CHI4*(PSI5+PSI6)
           DD   = MAX(BB*BB-4.d0*CC,ZERO)
           PSI4 =0.5d0*(-BB - SQRT(DD))
           PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
           PSI4 = TINY
        ENDIF
  !
  ! *** CALCULATE SPECIATION ********************************************
  !
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
        MOLAL (8) = PSI11                                 ! CAI
        MOLAL (9) = 2.D0*PSI9                             ! KI
        MOLAL (10)= PSI10                                 ! MGI
  !
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                    - MOLAL(9) - 2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
  !
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
  !
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CNA2SO4   = ZERO
        CK2SO4    = ZERO
        CMGSO4    = ZERO
        CCASO4    = CHI11
  !
        CALL CALCMR                                    ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !
  ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
  !
  !20    FUNCM8 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
  20    FUNCM8 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
  !
        RETURN
  !
  ! *** END OF FUNCTION FUNCM8 *******************************************
  !
        END FUNCTION FUNCM8

  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCP13
  ! *** CASE P13
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : CaSO4
  !     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
  !                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCP13
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,      &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,      &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,    &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,        &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,  &
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU  = .TRUE.
        CHI11   = MIN (W(2), W(6))                    ! CCASO4
        FRCA    = MAX (W(6) - CHI11, ZERO)
        FRSO4   = MAX (W(2) - CHI11, ZERO)
        CHI9    = MIN (FRSO4, 0.5D0*W(7))             ! CK2SO4
        FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
        FRSO4   = MAX (FRSO4 - CHI9, ZERO)
        CHI10   = FRSO4                               ! CMGSO4
        FRMG    = MAX (W(8) - CHI10, ZERO)
        CHI7    = MIN (W(1), W(5))                    ! CNACL
        FRNA    = MAX (W(1) - CHI7, ZERO)
        FRCL    = MAX (W(5) - CHI7, ZERO)
        CHI12   = MIN (FRCA, 0.5D0*W(4))              ! CCANO32
        FRCA    = MAX (FRCA - CHI12, ZERO)
        FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
        CHI17   = MIN (FRCA, 0.5D0*FRCL)              ! CCACL2
        FRCA    = MAX (FRCA - CHI17, ZERO)
        FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
        CHI15   = MIN (FRMG, 0.5D0*FRNO3)             ! CMGNO32
        FRMG    = MAX (FRMG - CHI15, ZERO)
        FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
        CHI16   = MIN (FRMG, 0.5D0*FRCL)              ! CMGCL2
        FRMG    = MAX (FRMG - CHI16, ZERO)
        FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
        CHI8    = MIN (FRNA, FRNO3)                   ! CNANO3
        FRNA    = MAX (FRNA - CHI8, ZERO)
        FRNO3   = MAX (FRNO3 - CHI8, ZERO)
        CHI14   = MIN (FRK, FRCL)                     ! CKCL
        FRK     = MAX (FRK - CHI14, ZERO)
        FRCL    = MAX (FRCL - CHI14, ZERO)
        CHI13   = MIN (FRK, FRNO3)                    ! CKNO3
        FRK     = MAX (FRK - CHI13, ZERO)
        FRNO3   = MAX (FRNO3 - CHI13, ZERO)
  !
        CHI5    = FRNO3                               ! HNO3(g)
        CHI6    = FRCL                                ! HCL(g)
        CHI4    = W(3)                                ! NH3(g)
  !
        CHI3    = ZERO                                ! CNH4CL
        CHI1    = ZERO
        CHI2    = ZERO
  !
        PSI6LO = TINY
        PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
  !
  ! *** INITIAL VALUES FOR BISECTION ************************************
  !
        X1 = PSI6LO
        Y1 = FUNCP13 (X1)
        IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
  !
  ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
  !
        DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
        DO 10 I=1,NDIV
           X2 = X1+DX
           Y2 = FUNCP13 (X2)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
           X1 = X2
           Y1 = Y2
  10    CONTINUE
  !
  ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
  !
        IF (ABS(Y2) .GT. EPS) THEN
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y2 = FUNCP13 (PSI6LO)
        ENDIF
        GOTO 50
  !
  ! *** PERFORM BISECTION ***********************************************
  !
  20    DO 30 I=1,MAXIT
           X3 = 0.5*(X1+X2)
           CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
           Y3 = FUNCP13 (X3)
           IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
              Y2    = Y3
              X2    = X3
           ELSE
              Y1    = Y3
              X1    = X3
           ENDIF
           IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
  30    CONTINUE
        CALL PUSHERR (0002, 'CALCP13')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CONVERGED ; RETURN **********************************************
  !
  40    X3 = 0.5*(X1+X2)
        CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
        Y3 = FUNCP13 (X3)
  !
  ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
  !
  50    CONTINUE
        IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
           CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
           MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
           MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
           MOLAL(6) = DELTA                                ! HSO4 EFFECT
        ENDIF
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCP13 ******************************************
  !
        END SUBROUTINE CALCP13

  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE FUNCP13
  ! *** CASE P13
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : CaSO4
  !     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
  !                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        REAL FUNCTION FUNCP13 (X)
        INCLUDE 'isrpia.inc'
  !
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        PSI6   = X
        PSI1   = ZERO
        PSI2   = ZERO
        PSI3   = ZERO
        PSI4   = ZERO
        PSI7   = CHI7
        PSI8   = CHI8
        PSI9   = CHI9
        PSI10  = CHI10
        PSI11  = ZERO
        PSI12  = CHI12
        PSI13  = CHI13
        PSI14  = CHI14
        PSI15  = CHI15
        PSI16  = CHI16
        PSI17  = CHI17
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
  !
  !  CALCULATE DISSOCIATION QUANTITIES
  !
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
               A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6) + PSI6 + PSI7 + PSI14 +        &
               2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
  !
        IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
           BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
           CC   = CHI4*(PSI5+PSI6)
           DD   = MAX(BB*BB-4.d0*CC,ZERO)
           PSI4 =0.5d0*(-BB - SQRT(DD))
           PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
           PSI4 = TINY
        ENDIF
  !
  ! *** CALCULATE SPECIATION *********************************************
  !
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
  !
  ! *** CALCULATE H+ *****************************************************
  !
  !      REST  = 2.D0*W(2) + W(4) + W(5)
  !C
  !      DELT1 = 0.0d0
  !      DELT2 = 0.0d0
  !      IF (W(1)+W(6)+W(7)+W(8).GT.REST) THEN
  !C
  !C *** CALCULATE EQUILIBRIUM CONSTANTS **********************************
  !C
  !      ALFA1 = XK26*RH*(WATER/1.0)                   ! CO2(aq) + H2O
  !      ALFA2 = XK27*(WATER/1.0)                      ! HCO3-
  !C
  !      X     = W(1)+W(6)+W(7)+W(8) - REST            ! EXCESS OF CRUSTALS EQUALS CO2(aq)
  !C
  !      DIAK  = SQRT( (ALFA1)**2.0 + 4.0D0*ALFA1*X)
  !      DELT1 = 0.5*(-ALFA1 + DIAK)
  !      DELT1 = MIN ( MAX (DELT1, ZERO), X)
  !      DELT2 = ALFA2
  !      DELT2 = MIN ( DELT2, DELT1)
  !      MOLAL(1) = DELT1 + DELT2                      ! H+
  !      ELSE
  !
  ! *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
  !
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
                    - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
  !      ENDIF
  !
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
  !
        CNH4NO3   = ZERO
        CNH4CL    = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CK2SO4    = ZERO
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = ZERO
        CKCL      = ZERO
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
  !
        CALL CALCMR                                    ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !
  ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
  !
  !20    FUNCP13 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
  20    FUNCP13 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
  !
        RETURN
  !
  ! *** END OF FUNCTION FUNCP13 *******************************************
  !
        END FUNCTION FUNCP13
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCL9
  ! *** CASE L9
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
  !     2. SOLID & LIQUID AEROSOL POSSIBLE
  !     3. SOLIDS POSSIBLE : CASO4
  !     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4, MGSO4, NA2SO4, K2SO4
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCL9
        INCLUDE 'isrpia.inc'
        REAL LAMDA
        COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,     &
                       CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,     &
                       CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,   &
                       PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,       &
                       PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
                       A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
  !
  ! *** FIND DRY COMPOSITION **********************************************
  !
        CALL CALCL1A
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CHI1 = CNH4HS4               ! Save from CALCL1 run
        CHI2 = CLC
        CHI3 = CNAHSO4
        CHI4 = CNA2SO4
        CHI5 = CNH42S4
        CHI6 = CK2SO4
        CHI7 = CMGSO4
        CHI8 = CKHSO4
  !
        PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
        PSI2 = CLC
        PSI3 = CNAHSO4
        PSI4 = CNA2SO4
        PSI5 = CNH42S4
        PSI6 = CK2SO4
        PSI7 = CMGSO4
        PSI8 = CKHSO4
  !
        CALAOU = .TRUE.              ! Outer loop activity calculation flag
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        A9 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
  !
  !  CALCULATE DISSOCIATION QUANTITIES
  !
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9              ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = MAX(BB*BB - 4.D0*CC, ZERO)
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
  !
  ! *** CALCULATE SPECIATION ********************************************
  !
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = PSI2 + PSI3 + PSI1 + PSI8 - LAMDA                ! HSO4I
        MOLAL(9) = PSI8 + 2.0D0*PSI6                                ! KI
        MOLAL(10)= PSI7                                             ! MGI
  !
        CLC      = ZERO
        CNAHSO4  = ZERO
        CNA2SO4  = ZERO
        CNH42S4  = ZERO
        CNH4HS4  = ZERO
        CK2SO4   = ZERO
        CMGSO4   = ZERO
        CKHSO4   = ZERO
  !
        CALL CALCMR                                         ! Water content

  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !
  20    RETURN
  !
  ! *** END OF SUBROUTINE CALCL9 *****************************************
  !
        END SUBROUTINE CALCL9
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCL1A
  ! *** CASE L1A
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
  !     2. SOLID AEROSOL ONLY
  !     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCL1A
        INCLUDE 'isrpia.inc'
  !
  ! *** CALCULATE NON VOLATILE SOLIDS ***********************************
  !
        CCASO4  = MIN (W(6), W(2))                    ! CCASO4
        FRSO4   = MAX(W(2) - CCASO4, ZERO)
        CAFR    = MAX(W(6) - CCASO4, ZERO)
        CK2SO4  = MIN (0.5D0*W(7), FRSO4)             ! CK2SO4
        FRK     = MAX(W(7) - 2.D0*CK2SO4, ZERO)
        FRSO4   = MAX(FRSO4 - CK2SO4, ZERO)
        CNA2SO4 = MIN (0.5D0*W(1), FRSO4)             ! CNA2SO4
        FRNA    = MAX(W(1) - 2.D0*CNA2SO4, ZERO)
        FRSO4   = MAX(FRSO4 - CNA2SO4, ZERO)
        CMGSO4  = MIN (W(8), FRSO4)                   ! CMGSO4
        FRMG    = MAX(W(8) - CMGSO4, ZERO)
        FRSO4   = MAX(FRSO4 - CMGSO4, ZERO)
  !
        CNH4HS4 = ZERO
        CNAHSO4 = ZERO
        CNH42S4 = ZERO
        CKHSO4  = ZERO
  !
        CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
        FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
        FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)
  !
        IF (FRSO4.LE.TINY) THEN
           CLC     = MAX(CLC - FRNH4, ZERO)
           CNH42S4 = 2.D0*FRNH4

        ELSEIF (FRNH4.LE.TINY) THEN
           CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
           CLC     = MAX(CLC-FRSO4, ZERO)
  !         IF (CK2SO4.GT.TINY) THEN
  !            FRSO4  = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
  !           CKHSO4 = 2.D0*FRSO4
  !            CK2SO4 = MAX(CK2SO4-FRSO4, ZERO)
  !         ENDIF
  !         IF (CNA2SO4.GT.TINY) THEN
  !            FRSO4   = MAX(FRSO4-CKHSO4/2.D0, ZERO)
  !            CNAHSO4 = 2.D0*FRSO4
  !            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
  !         ENDIF
  !
           IF (CNA2SO4.GT.TINY) THEN
              FRSO4  = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
              CNAHSO4 = 2.D0*FRSO4
              CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
           ENDIF
           IF (CK2SO4.GT.TINY) THEN
              FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
              CKHSO4 = 2.D0*FRSO4
              CK2SO4 = MAX(CK2SO4-FRSO4, ZERO)
         ENDIF
        ENDIF
  !
  ! *** CALCULATE GAS SPECIES ********************************************
  !
        GHNO3 = W(4)
        GHCL  = W(5)
        GNH3  = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCL1A *****************************************
  !
        END SUBROUTINE CALCL1A
  !
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCK4
  ! *** CASE K4
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : CASO4
  !
  ! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCK4
        INCLUDE 'isrpia.inc'
  !
        REAL :: LAMDA, KAPA
        COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
                       A1,   A2,   A3,   A4
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU =.TRUE.               ! Outer loop activity calculation flag
        FRST   = .TRUE.
        CALAIN = .TRUE.
  !
        CHI1   = W(3)                !  Total NH4 initially as NH4HSO4
        CHI2   = W(1)                !  Total NA initially as NaHSO4
        CHI3   = W(7)                !  Total K initially as KHSO4
        CHI4   = W(8)                !  Total Mg initially as MgSO4
  !
        LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  ! FREE H2SO4
        PSI1   = CHI1                            ! ALL NH4HSO4 DELIQUESCED
        PSI2   = CHI2                            ! ALL NaHSO4 DELIQUESCED
        PSI3   = CHI3                            ! ALL KHSO4 DELIQUESCED
        PSI4   = CHI4                            ! ALL MgSO4 DELIQUESCED
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
  !
        BB   = A4+LAMDA+PSI4                               ! KAPA
        CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
        DD   = MAX(BB*BB-4.D0*CC, ZERO)
        KAPA = 0.5D0*(-BB+SQRT(DD))
  !
  ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
  !
        MOLAL (1) = MAX(LAMDA + KAPA, TINY)                         ! HI
        MOLAL (2) = PSI2                                            ! NAI
        MOLAL (3) = PSI1                                            ! NH4I
        MOLAL (5) = MAX(KAPA + PSI4, ZERO)                          ! SO4I
        MOLAL (6) = MAX(LAMDA + PSI1 + PSI2 + PSI3 - KAPA, ZERO)    ! HSO4I
        MOLAL (9) = PSI3                                            ! KI
        MOLAL (10)= PSI4                                            ! MGI
  !
        CNH4HS4 = ZERO
        CNAHSO4 = ZERO
        CKHSO4  = ZERO
        CCASO4  = W(6)
        CMGSO4  = ZERO
  !
        CALL CALCMR                                      ! Water content
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !
  20    RETURN
  !
  ! *** END OF SUBROUTINE CALCK4
  !
        END SUBROUTINE CALCK4
