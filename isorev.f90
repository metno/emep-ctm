!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISRP1R
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF 
!     AN AMMONIUM-SULFATE AEROSOL SYSTEM. 
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY 
!     THE AMBIENT RELATIVE HUMIDITY.
!
! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      SUBROUTINE ISRP1R (WI, RHI, TEMPI)
        INCLUDE 'isrpia.inc' 
        DIMENSION WI(NCOMP)
  !
  ! *** INITIALIZE COMMON BLOCK VARIABLES *********************************
  !
        CALL INIT1 (WI, RHI, TEMPI)
  !
  ! *** CALCULATE SULFATE RATIO *******************************************
  !
        IF (RH.GE.DRNH42S4) THEN         ! WET AEROSOL, NEED NH4 AT SRATIO=2.0
           SULRATW = GETASR(WAER(2), RHI)     ! AEROSOL SULFATE RATIO
        ELSE
           SULRATW = 2.0D0                    ! DRY AEROSOL SULFATE RATIO
        ENDIF
        SULRAT  = WAER(3)/WAER(2)         ! SULFATE RATIO
  !
  ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
  !
  ! *** SULFATE POOR 
  !
        IF (SULRATW.LE.SULRAT) THEN
  !
  
           SCASE = 'S2'
           CALL CALCS2                 ! Only liquid (metastable)
  
  !
  ! *** SULFATE RICH (NO ACID)
  !
        ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN
        W(2) = WAER(2)
        W(3) = WAER(3)
  !
           SCASE = 'B4'
           CALL CALCB4                 ! Only liquid (metastable)
           SCASE = 'B4'
  !
        CALL CALCNH3P          ! Compute NH3(g)
  !
  ! *** SULFATE RICH (FREE ACID)
  !
        ELSEIF (SULRAT.LT.1.0) THEN             
        W(2) = WAER(2)
        W(3) = WAER(3)
  !
           SCASE = 'C2'
           CALL CALCC2                 ! Only liquid (metastable)
           SCASE = 'C2'
  
  ! 
        CALL CALCNH3P
  !
        ENDIF
        RETURN
  !
  ! *** END OF SUBROUTINE ISRP1R *****************************************
  !
        END SUBROUTINE ISRP1R
  
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE ISRP2R
  ! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF 
  !     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM. 
  !     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
  !     THE AMBIENT RELATIVE HUMIDITY.
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE ISRP2R (WI, RHI, TEMPI)
        INCLUDE 'isrpia.inc'
        DIMENSION WI(NCOMP)
        LOGICAL   TRYLIQ
  !
  ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
  !
        TRYLIQ = .TRUE.             ! Assume liquid phase, sulfate poor limit 
  !
  10    CALL INIT2 (WI, RHI, TEMPI)
  !
  ! *** CALCULATE SULFATE RATIO *******************************************
  !
        IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN ! *** WET AEROSOL
           SULRATW = GETASR(WAER(2), RHI)     ! LIMITING SULFATE RATIO
        ELSE
           SULRATW = 2.0D0                    ! *** DRY AEROSOL
        ENDIF
        SULRAT = WAER(3)/WAER(2)
  !
  ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
  !
  ! *** SULFATE POOR 
  !
        IF (SULRATW.LE.SULRAT) THEN                
  !
           SCASE = 'N3'
           CALL CALCN3                 ! Only liquid (metastable)
  
  ! *** SULFATE RICH (NO ACID)
  !
  !     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
  !     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE 
  !     AEROSOL EQUILIBRIUM.
  !
        ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
        W(2) = WAER(2)
        W(3) = WAER(3)
        W(4) = WAER(4)
  !
           SCASE = 'B4'
           CALL CALCB4                 ! Only liquid (metastable)
           SCASE = 'B4'
  !
  ! *** Add the NO3 to the solution now and calculate partitioning.
  !
        MOLAL(7) = WAER(4)             ! There is always water, so NO3(aer) is NO3-
        MOLAL(1) = MOLAL(1) + WAER(4)  ! Add H+ to balance out
        CALL CALCNAP            ! HNO3, NH3 dissolved
        CALL CALCNH3P
  !
  ! *** SULFATE RICH (FREE ACID)
  !
  !     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
  !     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE 
  !     AEROSOL EQUILIBRIUM.
  !
        ELSEIF (SULRAT.LT.1.0) THEN             
        W(2) = WAER(2)
        W(3) = WAER(3)
        W(4) = WAER(4)
  !
           SCASE = 'C2'
           CALL CALCC2                 ! Only liquid (metastable)
           SCASE = 'C2'
  
  !
  ! *** Add the NO3 to the solution now and calculate partitioning.
  !
        MOLAL(7) = WAER(4)             ! There is always water, so NO3(aer) is NO3-
        MOLAL(1) = MOLAL(1) + WAER(4)  ! Add H+ to balance out
  !
        CALL CALCNAP                   ! HNO3, NH3 dissolved
        CALL CALCNH3P
        ENDIF
  !
  ! *** IF SULRATW < SULRAT < 2.0 and WATER = 0 => SULFATE RICH CASE.
  !
        IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0  &
                                            .AND. WATER.LE.TINY) THEN
            TRYLIQ = .FALSE.
            GOTO 10
        ENDIF
  !
        RETURN
  !
  ! *** END OF SUBROUTINE ISRP2R *****************************************
  !
        END SUBROUTINE ISRP2R
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE ISRP3R
  ! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
  !     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM. 
  !     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM 
  !     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE ISRP3R (WI, RHI, TEMPI)
        INCLUDE 'isrpia.inc'
        DIMENSION WI(NCOMP)
        LOGICAL   TRYLIQ
  !cC
  !cC *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
  !cC
  !c      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
  !c      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
  !
  ! *** INITIALIZE ALL VARIABLES ******************************************
  !
        TRYLIQ = .TRUE.             ! Use liquid phase sulfate poor limit 
  !
  10    CALL ISOINIT3 (WI, RHI, TEMPI) ! COMMON block variables
  !cC
  !cC *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
  !cC
  !c      REST = 2.D0*WAER(2) + WAER(4) + WAER(5) 
  !c      IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
  !c         WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
  !c         CALL PUSHERR (0050, 'ISRP3R')     ! Warning error: Na adjusted
  !c      ENDIF
  !
  ! *** CALCULATE SULFATE & SODIUM RATIOS *********************************
  !
        IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN  ! ** WET AEROSOL
           FRSO4   = WAER(2) - WAER(1)/2.0D0     ! SULFATE UNBOUND BY SODIUM
           FRSO4   = MAX(FRSO4, TINY)
           SRI     = GETASR(FRSO4, RHI)          ! SULFATE RATIO FOR NH4+
           SULRATW = (WAER(1)+FRSO4*SRI)/WAER(2) ! LIMITING SULFATE RATIO
           SULRATW = MIN (SULRATW, 2.0D0)
        ELSE
           SULRATW = 2.0D0                     ! ** DRY AEROSOL
        ENDIF
        SULRAT = (WAER(1)+WAER(3))/WAER(2)
        SODRAT = WAER(1)/WAER(2)
  !
  ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
  !
  ! *** SULFATE POOR ; SODIUM POOR
  !
        IF (SULRATW.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                
  !
           SCASE = 'Q5'
           CALL CALCQ5                 ! Only liquid (metastable)
           SCASE = 'Q5'
  !
  ! *** SULFATE POOR ; SODIUM RICH
  !
        ELSE IF (SULRAT.GE.SULRATW .AND. SODRAT.GE.2.0) THEN                
  !
           SCASE = 'R6'
           CALL CALCR6                 ! Only liquid (metastable)
           SCASE = 'R6'
  !
  ! *** SULFATE RICH (NO ACID) 
  !
        ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
        DO 100 I=1,NCOMP
           W(I) = WAER(I)
  100   CONTINUE
  !
           SCASE = 'I6'
           CALL CALCI6                 ! Only liquid (metastable)
           SCASE = 'I6'
  !
        CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
        CALL CALCNH3P
  !
  ! *** SULFATE RICH (FREE ACID)
  !
        ELSEIF (SULRAT.LT.1.0) THEN             
        DO 200 I=1,NCOMP
           W(I) = WAER(I)
  200   CONTINUE
  !
           SCASE = 'J3'
           CALL CALCJ3                 ! Only liquid (metastable)
           SCASE = 'J3'
  !
        CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
        CALL CALCNH3P
  !
        ENDIF
  !
  ! *** IF AFTER CALCULATIONS, SULRATW < SULRAT < 2.0  
  !                            and WATER = 0          => SULFATE RICH CASE.
  !
        IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0  &
                              .AND. WATER.LE.TINY) THEN
            TRYLIQ = .FALSE.
            GOTO 10
        ENDIF
  !
        RETURN
  !
  ! *** END OF SUBROUTINE ISRP3R *****************************************
  !
        END SUBROUTINE ISRP3R
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE ISRP4R
  ! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
  !     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM-CALCIUM-POTTASIUM-MAGNESIUM AEROSOL SYSTEM.
  !     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
  !     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE ISRP4R (WI, RHI, TEMPI)
        INCLUDE 'isrpia.inc'
        DIMENSION WI(NCOMP)
        LOGICAL   TRYLIQ
  !cC
  !cC *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
  !cC
  !c      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
  !c      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
  !
  ! *** INITIALIZE ALL VARIABLES ******************************************
  !
        TRYLIQ  = .TRUE.             ! Use liquid phase sulfate poor limit
        IPROB   = 1            ! SOLVE REVERSE PROBLEM
  !
  10    CALL INIT4 (WI, RHI, TEMPI) ! COMMON block variables
  !cC
  !cC *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
  !cC
  !c      REST = 2.D0*WAER(2) + WAER(4) + WAER(5)
  !c      IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
  !c         WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
  !c         CALL PUSHERR (0050, 'ISRP3R')     ! Warning error: Na adjusted
  !c      ENDIF
  !
  ! *** CALCULATE SULFATE, CRUSTAL & SODIUM RATIOS ***********************
  !
        IF (TRYLIQ) THEN                               ! ** WET AEROSOL
           FRSO4   = WAER(2) - WAER(1)/2.0D0 &
                   - WAER(6) - WAER(7)/2.0D0 - WAER(8) ! SULFATE UNBOUND BY SODIUM,CALCIUM,POTTASIUM,MAGNESIUM
           FRSO4   = MAX(FRSO4, TINY)
           SRI     = GETASR(FRSO4, RHI)                ! SULFATE RATIO FOR NH4+
           SULRATW = (WAER(1)+FRSO4*SRI+WAER(6) &
                      +WAER(7)+WAER(8))/WAER(2)       ! LIMITING SULFATE RATIO
           SULRATW = MIN (SULRATW, 2.0D0)
        ELSE
           SULRATW = 2.0D0                     ! ** DRY AEROSOL
        ENDIF
        SO4RAT = (WAER(1)+WAER(3)+WAER(6)+WAER(7)+WAER(8))/WAER(2)
        CRNARAT = (WAER(1)+WAER(6)+WAER(7)+WAER(8))/WAER(2)
        CRRAT  = (WAER(6)+WAER(7)+WAER(8))/WAER(2)
  !
  ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
  !
  ! *** SULFATE POOR ; SODIUM+CRUSTALS POOR
  !
        IF (SULRATW.LE.SO4RAT .AND. CRNARAT.LT.2.0) THEN
  !
           SCASE = 'V7'
           CALL CALCV7                 ! Only liquid (metastable)
  !
  ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
  !
        ELSEIF (SO4RAT.GE.SULRATW .AND. CRNARAT.GE.2.0) THEN
  !
         IF (CRRAT.LE.2.0) THEN
  !
           SCASE = 'U8'
           CALL CALCU8                 ! Only liquid (metastable)
  !
  ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
  !
         ELSEIF (CRRAT.GT.2.0) THEN
  !
           SCASE = 'W13'
           CALL CALCW13                 ! Only liquid (metastable)
  
  !        CALL CALCNH3
         ENDIF
  !
  ! *** SULFATE RICH (NO ACID): 1<Rso4<2;
  !
        ELSEIF (1.0.LE.SO4RAT .AND. SO4RAT.LT.SULRATW) THEN
        DO 800 I=1,NCOMP
           W(I) = WAER(I)
   800  CONTINUE
  !
           SCASE = 'L9'
           CALL CALCL9                 ! Only liquid (metastable)
  
  !
        CALL CALCNHP                ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3P               !                NH3
  !
  ! *** SULFATE SUPER RICH (FREE ACID): Rso4<1;
  !
        ELSEIF (SO4RAT.LT.1.0) THEN
        DO 900 I=1,NCOMP
           W(I) = WAER(I)
   900  CONTINUE
  !
           SCASE = 'K4'
           CALL CALCK4                 ! Only liquid (metastable)
  
  !
        CALL CALCNHP                  ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3P                 !                NH3
  !
        ENDIF
  !
  ! *** IF AFTER CALCULATIONS, SULRATW < SO4RAT < 2.0
  !                            and WATER = 0          => SULFATE RICH CASE.
  !
        IF (SULRATW.LE.SO4RAT .AND. SO4RAT.LT.2.0 &
                              .AND. WATER.LE.TINY) THEN
            TRYLIQ = .FALSE.
            GOTO 10
        ENDIF
  !
        RETURN
  !
  ! *** END OF SUBROUTINE ISRP4R *****************************************
  !
        END SUBROUTINE ISRP4R 
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCS2
  ! *** CASE S2
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0)
  !     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCS2
        INCLUDE 'isrpia.inc'
        REAL :: NH4I, NH3GI, NH3AQ
  !
  ! *** SETUP PARAMETERS ************************************************
  !
        CALAOU   =.TRUE.     ! Outer loop activity calculation flag
        FRST     =.TRUE.
        CALAIN   =.TRUE.
  !
  ! *** CALCULATE WATER CONTENT *****************************************
  !
        MOLALR(4)= MIN(WAER(2), 0.5d0*WAER(3))
        WATER    = MOLALR(4)/M0(4)  ! ZSR correlation
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !C         A21  = XK21*WATER*R*TEMP
           A2   = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
           AKW  = XKW *RH*WATER*WATER
  !
           NH4I = WAER(3)
           SO4I = WAER(2)
           HSO4I= ZERO
  !
           CALL CALCPH (2.D0*SO4I - NH4I, HI, OHI)    ! Get pH
  !
           NH3AQ = ZERO                               ! AMMONIA EQUILIBRIUM
           IF (HI.LT.OHI) THEN
              CALL CALCAMAQ (NH4I, OHI, DEL)
              NH4I  = MAX (NH4I-DEL, ZERO) 
              OHI   = MAX (OHI -DEL, TINY)
              NH3AQ = DEL
              HI    = AKW/OHI
           ENDIF
  !
           CALL CALCHS4 (HI, SO4I, ZERO, DEL)         ! SULFATE EQUILIBRIUM
           SO4I  = SO4I - DEL
           HI    = HI   - DEL
           HSO4I = DEL
  !
           NH3GI = NH4I/HI/A2   !    NH3AQ/A21
  !
  ! *** SPECIATION & WATER CONTENT ***************************************
  !
           MOLAL(1) = HI
           MOLAL(3) = NH4I
           MOLAL(5) = SO4I
           MOLAL(6) = HSO4I
           COH      = OHI
           GASAQ(1) = NH3AQ
           GNH3     = NH3GI
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
  ! *** END OF SUBROUTINE CALCS2 ****************************************
  !
        END SUBROUTINE CALCS2
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCN3
  ! *** CASE N3
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0)
  !     2. THERE IS ONLY A LIQUID PHASE
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCN3
        INCLUDE 'isrpia.inc'
        REAL :: NH4I, NO3I, NH3AQ, NO3AQ
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
        CALAOU =.TRUE.              ! Outer loop activity calculation flag
        FRST   =.TRUE.
        CALAIN =.TRUE.
  !
  ! *** AEROSOL WATER CONTENT
  !
        MOLALR(4) = MIN(WAER(2),0.5d0*WAER(3))       ! (NH4)2SO4
        AML5      = MAX(WAER(3)-2.D0*MOLALR(4),ZERO) ! "free" NH4
        MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO)     ! NH4NO3=MIN("free",NO3)
        WATER     = MOLALR(4)/M0(4) + MOLALR(5)/M0(5)
        WATER     = MAX(WATER, TINY)
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
           A2    = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
  !C         A21   = XK21*WATER*R*TEMP
           A3    = XK4*R*TEMP*(WATER/GAMA(10))**2.0
           A4    = XK7*(WATER/GAMA(4))**3.0
           AKW   = XKW *RH*WATER*WATER
  !
  ! ION CONCENTRATIONS
  !
           NH4I  = WAER(3)
           NO3I  = WAER(4)
           SO4I  = WAER(2)
           HSO4I = ZERO
  !
           CALL CALCPH (2.D0*SO4I + NO3I - NH4I, HI, OHI)
  !
  ! AMMONIA ASSOCIATION EQUILIBRIUM
  !
           NH3AQ = ZERO
           NO3AQ = ZERO
           GG    = 2.D0*SO4I + NO3I - NH4I
           IF (HI.LT.OHI) THEN
              CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
              HI    = AKW/OHI
           ELSE
              HI    = ZERO
              CALL CALCNIAQ2 (GG, NO3I, HI, NO3AQ) ! HNO3
  !
  ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
  !
              CALL CALCHS4 (HI, SO4I, ZERO, DEL)
              SO4I  = SO4I  - DEL
              HI    = HI    - DEL
              HSO4I = DEL
              OHI   = AKW/HI
           ENDIF
  !
  ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
  !
           MOLAL (1) = HI
           MOLAL (3) = NH4I
           MOLAL (5) = SO4I
           MOLAL (6) = HSO4I
           MOLAL (7) = NO3I
           COH       = OHI
  !
           CNH42S4   = ZERO
           CNH4NO3   = ZERO
  !
           GASAQ(1)  = NH3AQ
           GASAQ(3)  = NO3AQ
  !
           GHNO3     = HI*NO3I/A3
           GNH3      = NH4I/HI/A2   !   NH3AQ/A21 
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
  !
           IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
              CALL CALCACT     
           ELSE
              GOTO 20
           ENDIF
  10    CONTINUE
  !
  ! *** RETURN ***********************************************************
  !
  20    RETURN
  !
  ! *** END OF SUBROUTINE CALCN3 *****************************************
  !
        END SUBROUTINE CALCN3
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCQ5
  ! *** CASE Q5
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
  !     2. LIQUID AND SOLID PHASES ARE POSSIBLE
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCQ5
        INCLUDE 'isrpia.inc'
  !
        REAL :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ
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
        FRST    =.TRUE.
        CALAIN  =.TRUE. 
        CALAOU  =.TRUE.
  !
  ! *** CALCULATE INITIAL SOLUTION ***************************************
  !
        CALL CALCQ1A
  !
        PSI1   = CNA2SO4      ! SALTS DISSOLVED
        PSI4   = CNH4CL
        PSI5   = CNH4NO3
        PSI6   = CNH42S4
  !
        CALL CALCMR           ! WATER
  !
        NH3AQ  = ZERO
        NO3AQ  = ZERO
        CLAQ   = ZERO
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
  !
  ! ION CONCENTRATIONS
  !
        NAI    = WAER(1)
        SO4I   = WAER(2)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
  !
  ! SOLUTION ACIDIC OR BASIC?
  !
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG.GT.TINY) THEN                        ! H+ in excess
           BB =-GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           HI = 0.5D0*(-BB + SQRT(DD))
           OHI= AKW/HI
        ELSE                                        ! OH- in excess
           BB = GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           OHI= 0.5D0*(-BB + SQRT(DD))
           HI = AKW/OHI
        ENDIF
  !
  ! UNDISSOCIATED SPECIES EQUILIBRIA
  !
        IF (HI.LT.OHI) THEN
           CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
           HI    = AKW/OHI
           HSO4I = ZERO
        ELSE
           GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
           GGCL  = MAX(GG-GGNO3, ZERO)
           IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
           IF (GGNO3.GT.TINY) THEN
              IF (GGCL.LE.TINY) HI = ZERO
              CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
           ENDIF
  !
  ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
  !
           CALL CALCHS4 (HI, SO4I, ZERO, DEL)
           SO4I  = SO4I  - DEL
           HI    = HI    - DEL
           HSO4I = DEL
           OHI   = AKW/HI
        ENDIF
  !
  ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
  !
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !cc      CALL PUSHERR (0002, 'CALCQ5')    ! WARNING ERROR: NO CONVERGENCE
  ! 
  ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
  !
  20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
  !
        GNH3    = NH4I/HI/A2
        GHNO3   = HI*NO3I/A3
        GHCL    = HI*CLI /A4
  !
        GASAQ(1)= NH3AQ
        GASAQ(2)= CLAQ
        GASAQ(3)= NO3AQ
  !
        CNH42S4 = ZERO
        CNH4NO3 = ZERO
        CNH4CL  = ZERO
        CNACL   = ZERO
        CNANO3  = ZERO
        CNA2SO4 = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCQ5 ******************************************
  !
        END SUBROUTINE CALCQ5
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCQ1A
  ! *** CASE Q1 ; SUBCASE 1
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
  !     2. SOLID AEROSOL ONLY
  !     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, (NH4)2SO4, NA2SO4
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCQ1A
        INCLUDE 'isrpia.inc'
  !
  ! *** CALCULATE SOLIDS **************************************************
  !
        CNA2SO4 = 0.5d0*WAER(1)
        FRSO4   = MAX (WAER(2)-CNA2SO4, ZERO)
  !
        CNH42S4 = MAX (MIN(FRSO4,0.5d0*WAER(3)), TINY)
        FRNH3   = MAX (WAER(3)-2.D0*CNH42S4, ZERO)
  !
        CNH4NO3 = MIN (FRNH3, WAER(4))
  !CC      FRNO3   = MAX (WAER(4)-CNH4NO3, ZERO)
        FRNH3   = MAX (FRNH3-CNH4NO3, ZERO)
  !
        CNH4CL  = MIN (FRNH3, WAER(5))
  !CC      FRCL    = MAX (WAER(5)-CNH4CL, ZERO)
        FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
  !
  ! *** OTHER PHASES ******************************************************
  !
        WATER   = ZERO
  !
        GNH3    = ZERO
        GHNO3   = ZERO
        GHCL    = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCQ1A *****************************************
  !
        END SUBROUTINE CALCQ1A
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCR6
  ! *** CASE R6
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
  !     2. THERE IS ONLY A LIQUID PHASE
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCR6
        INCLUDE 'isrpia.inc'
  !
        REAL :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ
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
        CALL CALCR1A
  !
        PSI1   = CNA2SO4
        PSI2   = CNANO3
        PSI3   = CNACL
        PSI4   = CNH4CL
        PSI5   = CNH4NO3
  !
        FRST   = .TRUE.
        CALAIN = .TRUE. 
        CALAOU = .TRUE. 
  !
  ! *** CALCULATE WATER **************************************************
  !
        CALL CALCMR
  !
  ! *** SETUP LIQUID CONCENTRATIONS **************************************
  !
        HSO4I  = ZERO
        NH3AQ  = ZERO
        NO3AQ  = ZERO
        CLAQ   = ZERO
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
        AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+      
  !
        NAI    = WAER(1)
        SO4I   = WAER(2)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
  !
  ! SOLUTION ACIDIC OR BASIC?
  !
        GG  = 2.D0*WAER(2) + NO3I + CLI - NAI - NH4I
        IF (GG.GT.TINY) THEN                        ! H+ in excess
           BB =-GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           HI = 0.5D0*(-BB + SQRT(DD))
           OHI= AKW/HI
        ELSE                                        ! OH- in excess
           BB = GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           OHI= 0.5D0*(-BB + SQRT(DD))
           HI = AKW/OHI
        ENDIF
  !
  ! UNDISSOCIATED SPECIES EQUILIBRIA
  !
        IF (HI.LT.OHI) THEN
           CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
           HI    = AKW/OHI
        ELSE
           GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
           GGCL  = MAX(GG-GGNO3, ZERO)
           IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
           IF (GGNO3.GT.TINY) THEN
              IF (GGCL.LE.TINY) HI = ZERO
              CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
           ENDIF
  !
  ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
  !
           CALL CALCHS4 (HI, SO4I, ZERO, DEL)
           SO4I  = SO4I  - DEL
           HI    = HI    - DEL
           HSO4I = DEL
           OHI   = AKW/HI
        ENDIF
  !
  ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
  !
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !cc      CALL PUSHERR (0002, 'CALCR6')    ! WARNING ERROR: NO CONVERGENCE
  ! 
  ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
  !
  20    A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
  !
        GNH3     = NH4I/HI/A2
        GHNO3    = HI*NO3I/A3
        GHCL     = HI*CLI /A4
  !
        GASAQ(1) = NH3AQ
        GASAQ(2) = CLAQ
        GASAQ(3) = NO3AQ
  !
        CNH42S4  = ZERO
        CNH4NO3  = ZERO
        CNH4CL   = ZERO
        CNACL    = ZERO
        CNANO3   = ZERO
        CNA2SO4  = ZERO 
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCR6 ******************************************
  !
        END SUBROUTINE CALCR6
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCR1A
  ! *** CASE R1 ; SUBCASE 1
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
  !     2. SOLID AEROSOL ONLY
  !     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4, NANO3, NACL
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCR1A
        INCLUDE 'isrpia.inc'
  !
  ! *** CALCULATE SOLIDS **************************************************
  !
        CNA2SO4 = WAER(2)
        FRNA    = MAX (WAER(1)-2*CNA2SO4, ZERO)
  !
        CNH42S4 = ZERO
  !
        CNANO3  = MIN (FRNA, WAER(4))
        FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
        FRNA    = MAX (FRNA-CNANO3, ZERO)
  !
        CNACL   = MIN (FRNA, WAER(5))
        FRCL    = MAX (WAER(5)-CNACL, ZERO)
        FRNA    = MAX (FRNA-CNACL, ZERO)
  !
        CNH4NO3 = MIN (FRNO3, WAER(3))
        FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
        FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)
  !
        CNH4CL  = MIN (FRCL, FRNH3)
        FRCL    = MAX (FRCL-CNH4CL, ZERO)
        FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
  !
  ! *** OTHER PHASES ******************************************************
  !
        WATER   = ZERO
  !
        GNH3    = ZERO
        GHNO3   = ZERO
        GHCL    = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCR1A *****************************************
  !
        END SUBROUTINE CALCR1A
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCV7
  ! *** CASE V7
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCV7
        INCLUDE 'isrpia.inc'
  !
        REAL :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI
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
        FRST    =.TRUE.
        CALAIN  =.TRUE.
        CALAOU  =.TRUE.
  !
  ! *** CALCULATE INITIAL SOLUTION ***************************************
  !
        CALL CALCV1A
  !
        CHI9   = CCASO4
  !
        PSI1   = CNA2SO4      ! SALTS DISSOLVED
        PSI4   = CNH4CL
        PSI5   = CNH4NO3
        PSI6   = CNH42S4
        PSI7   = CK2SO4
        PSI8   = CMGSO4
        PSI9   = CCASO4
  !
        CALL CALCMR           ! WATER
  !
        NH3AQ  = ZERO
        NO3AQ  = ZERO
        CLAQ   = ZERO
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
  !
  ! ION CONCENTRATIONS
  !
        NAI    = WAER(1)
        SO4I   = MAX (WAER(2) - WAER(6), ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        KI     = WAER(7)
        MGI    = WAER(8)
  !
  ! SOLUTION ACIDIC OR BASIC?
  !
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
               - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG.GT.TINY) THEN                    ! H+ in excess
           BB =-GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           HI = 0.5D0*(-BB + SQRT(DD))
           OHI= AKW/HI
        ELSE                                     ! OH- in excess
           BB = GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           OHI= 0.5D0*(-BB + SQRT(DD))
           HI = AKW/OHI
        ENDIF
  
  !
  ! UNDISSOCIATED SPECIES EQUILIBRIA
  !
        IF (HI.GT.OHI) THEN
  !         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
  !         HI    = AKW/OHI
  !         HSO4I = ZERO
  !      ELSE
  !         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I - 2.D0*CAI
  !     &           - KI - 2.D0*MGI, ZERO)
  !         GGCL  = MAX(GG-GGNO3, ZERO)
  !         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
  !         IF (GGNO3.GT.TINY) THEN
  !            IF (GGCL.LE.TINY) HI = ZERO
  !            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
  !         ENDIF
  !
  ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
  !
           CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
          del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
  !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
  !
        IF (HI.LE.TINY) THEN
        HI = SQRT(AKW)
        OHI   = AKW/HI
        ENDIF
  !
  ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
  !
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
        MOLAL(8) = CAI
        MOLAL(9) = KI
        MOLAL(10)= MGI
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !cc      CALL PUSHERR (0002, 'CALCV7')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
  !
  20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
  !
        GNH3    = NH4I/HI/A2
        GHNO3   = HI*NO3I/A3
        GHCL    = HI*CLI /A4
  !
        GASAQ(1)= NH3AQ
        GASAQ(2)= CLAQ
        GASAQ(3)= NO3AQ
  !
        CNH42S4 = ZERO
        CNH4NO3 = ZERO
        CNH4CL  = ZERO
        CNA2SO4 = ZERO
        CMGSO4  = ZERO
        CK2SO4  = ZERO
        CCASO4  = MIN (WAER(6), WAER(2))
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCV7 ******************************************
  !
        END SUBROUTINE CALCV7
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCV1A
  ! *** CASE V1A
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
  !     2. SOLID AEROSOL ONLY
  !     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCV1A
        INCLUDE 'isrpia.inc'
  !
  ! *** CALCULATE SOLIDS **************************************************
  !
        CCASO4  = MIN (WAER(6), WAER(2))                     ! CCASO4
        SO4FR   = MAX (WAER(2) - CCASO4, ZERO)
        CAFR    = MAX (WAER(6) - CCASO4, ZERO)
        CK2SO4  = MIN (0.5D0*WAER(7), SO4FR)                 ! CK2SO4
        FRK     = MAX (WAER(7) - 2.D0*CK2SO4, ZERO)
        SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
        CNA2SO4 = MIN (0.5D0*WAER(1), SO4FR)                 ! CNA2SO4
        NAFR    = MAX (WAER(1) - 2.D0*CNA2SO4, ZERO)
        SO4FR   = MAX (SO4FR - CNA2SO4, ZERO)
        CMGSO4  = MIN (WAER(8), SO4FR)                       ! CMGSO4
        FRMG    = MAX(WAER(8) - CMGSO4, ZERO)
        SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
        CNH42S4 = MAX (MIN (SO4FR , 0.5d0*WAER(3)) , TINY)
        FRNH3   = MAX (WAER(3) - 2.D0*CNH42S4, ZERO)
  !
        CNH4NO3 = MIN (FRNH3, WAER(4))
  !CC      FRNO3   = MAX (WAER(4) - CNH4NO3, ZERO)
        FRNH3   = MAX (FRNH3 - CNH4NO3, ZERO)
  !
        CNH4CL  = MIN (FRNH3, WAER(5))
  !CC      FRCL    = MAX (WAER(5) - CNH4CL, ZERO)
        FRNH3   = MAX (FRNH3 - CNH4CL, ZERO)
  !
  ! *** OTHER PHASES ******************************************************
  !
        WATER   = ZERO
  !
        GNH3    = ZERO
        GHNO3   = ZERO
        GHCL    = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCV1A *****************************************
  !
        END SUBROUTINE CALCV1A
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCU8
  ! *** CASE U8
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0); CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
  !     2. THERE IS ONLY A LIQUID PHASE
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCU8
        INCLUDE 'isrpia.inc'
  !
        REAL :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI
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
        CALL CALCU1A
  !
        CHI9   = CCASO4        ! SALTS
  !
        PSI1   = CNA2SO4
        PSI2   = CNANO3
        PSI3   = CNACL
        PSI4   = CNH4CL
        PSI5   = CNH4NO3
        PSI7   = CK2SO4
        PSI8   = CMGSO4
        PSI9   = CCASO4
  !
        FRST   = .TRUE.
        CALAIN = .TRUE.
        CALAOU = .TRUE.
  !
  ! *** CALCULATE WATER **************************************************
  !
        CALL CALCMR
  !
  ! *** SETUP LIQUID CONCENTRATIONS **************************************
  !
        HSO4I  = ZERO
        NH3AQ  = ZERO
        NO3AQ  = ZERO
        CLAQ   = ZERO
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  !
        AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+
  !
        NAI    = WAER(1)
        SO4I   = MAX(WAER(2) - WAER(6), ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        KI     = WAER(7)
        MGI    = WAER(8)
  
  !
  ! SOLUTION ACIDIC OR BASIC?
  !
        GG  = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
               - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG.GT.TINY) THEN                        ! H+ in excess
           BB =-GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           HI = 0.5D0*(-BB + SQRT(DD))
           OHI= AKW/HI
        ELSE                                        ! OH- in excess
           BB = GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           OHI= 0.5D0*(-BB + SQRT(DD))
           HI = AKW/OHI
        ENDIF
        IF (HI.LE.TINY) HI = SQRT(AKW)
  !      OHI   = AKW/HI
  !
  ! UNDISSOCIATED SPECIES EQUILIBRIA
  !
        IF (HI.GT.OHI) THEN
  !         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
  !         HI    = AKW/OHI
  !         HSO4I = ZERO
  !      ELSE
  !         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I - 2.D0*CAI
  !     &           - KI - 2.D0*MGI, ZERO)
  !         GGCL  = MAX(GG-GGNO3, ZERO)
  !         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
  !         IF (GGNO3.GT.TINY) THEN
  !            IF (GGCL.LE.TINY) HI = ZERO
  !            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
  !         ENDIF
  !
  ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
  !
           CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
          del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
  !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
  !
        IF (HI.LE.TINY) THEN
        HI = SQRT(AKW)
        OHI   = AKW/HI
        ENDIF
  !
  ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
  !
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
        MOLAL(8) = CAI
        MOLAL(9) = KI
        MOLAL(10)= MGI
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !cc      CALL PUSHERR (0002, 'CALCU8')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
  !
  20    A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
  !
        GNH3     = NH4I/HI/A2
        GHNO3    = HI*NO3I/A3
        GHCL     = HI*CLI /A4
  !
        GASAQ(1) = NH3AQ
        GASAQ(2) = CLAQ
        GASAQ(3) = NO3AQ
  !
        CNH42S4  = ZERO
        CNH4NO3  = ZERO
        CNH4CL   = ZERO
        CNACL    = ZERO
        CNANO3   = ZERO
        CNA2SO4  = ZERO
        CMGSO4   = ZERO
        CK2SO4   = ZERO
        CCASO4   = MIN (WAER(6), WAER(2))
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCU8 ******************************************
  !
        END SUBROUTINE CALCU8
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE
  ! *** SUBROUTINE CALCU1A
  ! *** CASE U1A
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0); CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
  !     2. THERE IS ONLY A SOLID PHASE
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCU1A
        INCLUDE 'isrpia.inc'
  !
  ! *** CALCULATE SOLIDS *************************************************
  !
        CCASO4  = MIN (WAER(6), WAER(2))                 ! CCASO4
        SO4FR   = MAX(WAER(2) - CCASO4, ZERO)
        CAFR    = MAX(WAER(6) - CCASO4, ZERO)
        CK2SO4  = MIN (0.5D0*WAER(7), SO4FR)             ! CK2SO4
        FRK     = MAX(WAER(7) - 2.D0*CK2SO4, ZERO)
        SO4FR   = MAX(SO4FR - CK2SO4, ZERO)
        CMGSO4  = MIN (WAER(8), SO4FR)                   ! CMGSO4
        FRMG    = MAX(WAER(8) - CMGSO4, ZERO)
        SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
        CNA2SO4 = MAX (SO4FR, ZERO)                      ! CNA2SO4
        FRNA    = MAX (WAER(1) - 2.D0*CNA2SO4, ZERO)
  !
        CNH42S4 = ZERO
  !
        CNANO3  = MIN (FRNA, WAER(4))
        FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
        FRNA    = MAX (FRNA-CNANO3, ZERO)
  !
        CNACL   = MIN (FRNA, WAER(5))
        FRCL    = MAX (WAER(5)-CNACL, ZERO)
        FRNA    = MAX (FRNA-CNACL, ZERO)
  !
        CNH4NO3 = MIN (FRNO3, WAER(3))
        FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
        FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)
  !
        CNH4CL  = MIN (FRCL, FRNH3)
        FRCL    = MAX (FRCL-CNH4CL, ZERO)
        FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
  !
  ! *** OTHER PHASES ******************************************************
  !
        WATER   = ZERO
  !
        GNH3    = ZERO
        GHNO3   = ZERO
        GHCL    = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCU1A *****************************************
  !
        END SUBROUTINE CALCU1A
  !
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCW13
  ! *** CASE W13
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
  !     2. THERE IS BOTH A LIQUID & SOLID PHASE
  !     3. SOLIDS POSSIBLE : CaSO4
  !     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
  !                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCW13
        INCLUDE 'isrpia.inc'
  !
        REAL :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI
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
        FRST    =.TRUE.
        CALAIN  =.TRUE.
        CALAOU  =.TRUE.
  !
  ! *** CALCULATE INITIAL SOLUTION ***************************************
  !
        CALL CALCW1A
  !
        CHI11   = CCASO4
  !
        PSI1   = CNA2SO4      ! SALTS DISSOLVED
        PSI5   = CNH4CL
        PSI6   = CNH4NO3
        PSI7   = CNACL
        PSI8   = CNANO3
        PSI9   = CK2SO4
        PSI10  = CMGSO4
        PSI11  = CCASO4
        PSI12  = CCANO32
        PSI13  = CKNO3
        PSI14  = CKCL
        PSI15  = CMGNO32
        PSI16  = CMGCL2
        PSI17  = CCACL2
  !
        CALL CALCMR           ! WATER
  !
        NH3AQ  = ZERO
        NO3AQ  = ZERO
        CLAQ   = ZERO
  !
  ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
  !
        DO 10 I=1,NSWEEP
  
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
  !
  ! ION CONCENTRATIONS
  !
        NAI    = WAER(1)
        SO4I   = MAX (WAER(2) - WAER(6), ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        KI     = WAER(7)
        MGI    = WAER(8)
  !
  ! SOLUTION ACIDIC OR BASIC?
  !
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
               - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG.GT.TINY) THEN                        ! H+ in excess
           BB =-GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           HI = 0.5D0*(-BB + SQRT(DD))
           OHI= AKW/HI
        ELSE                                        ! OH- in excess
           BB = GG
           CC =-AKW
           DD = BB*BB - 4.D0*CC
           OHI= 0.5D0*(-BB + SQRT(DD))
           HI = AKW/OHI
        ENDIF
  !
  ! UNDISSOCIATED SPECIES EQUILIBRIA
  !
        IF (HI.GT.OHI) THEN
  !         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
  !         HI    = AKW/OHI
  !         HSO4I = ZERO
  !      ELSE
  !         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I - 2.D0*CAI
  !     &           - KI - 2.D0*MGI, ZERO)
  !         GGCL  = MAX(GG-GGNO3, ZERO)
  !         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
  !         IF (GGNO3.GT.TINY) THEN
  !            IF (GGCL.LE.TINY) HI = ZERO
  !            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
  !         ENDIF
  !
  ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
  !
           CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
          del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
  !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
  !
        IF (HI.LE.TINY) THEN
        HI = SQRT(AKW)
        OHI   = AKW/HI
        ENDIF
  !
  ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
  !
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
        MOLAL(8) = CAI
        MOLAL(9) = KI
        MOLAL(10)= MGI
  !
  ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
  !
        IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
           CALL CALCACT
        ELSE
           GOTO 20
        ENDIF
  10    CONTINUE
  !cc      CALL PUSHERR (0002, 'CALCW13')    ! WARNING ERROR: NO CONVERGENCE
  !
  ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
  !
  20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
  !
        GNH3    = NH4I/HI/A2
        GHNO3   = HI*NO3I/A3
        GHCL    = HI*CLI /A4
  !
        GASAQ(1)= NH3AQ
        GASAQ(2)= CLAQ
        GASAQ(3)= NO3AQ
  !
        CNH42S4 = ZERO
        CNH4NO3 = ZERO
        CNH4CL  = ZERO
        CNACL   = ZERO
        CNANO3  = ZERO
        CMGSO4  = ZERO
        CK2SO4  = ZERO
        CCASO4  = MIN (WAER(6), WAER(2))
        CCANO32 = ZERO
        CKNO3   = ZERO
        KCL     = ZERO
        CMGNO32 = ZERO
        CMGCL2  = ZERO
        CCACL2  = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCW13 ******************************************
  !
        END  SUBROUTINE CALCW13
  !=======================================================================
  !
  ! *** ISORROPIA CODE II
  ! *** SUBROUTINE CALCW1A
  ! *** CASE W1A
  !
  !     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
  !     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
  !     2. SOLID AEROSOL ONLY
  !     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
  !                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
  !
  ! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
  ! *** GEORGIA INSTITUTE OF TECHNOLOGY
  ! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
  !
  !=======================================================================
  !
        SUBROUTINE CALCW1A
        INCLUDE 'isrpia.inc'
  !
  ! *** CALCULATE SOLIDS **************************************************
  !
        CCASO4  = MIN (WAER(2), WAER(6))              !SOLID CASO4
        CAFR    = MAX (WAER(6) - CCASO4, ZERO)
        SO4FR   = MAX (WAER(2) - CCASO4, ZERO)
        CK2SO4  = MIN (SO4FR, 0.5D0*WAER(7))          !SOLID K2SO4
        FRK     = MAX (WAER(7) - 2.D0*CK2SO4, ZERO)
        SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
        CMGSO4  = SO4FR                               !SOLID MGSO4
        FRMG    = MAX (WAER(8) - CMGSO4, ZERO)
        CNACL   = MIN (WAER(1), WAER(5))              !SOLID NACL
        FRNA    = MAX (WAER(1) - CNACL, ZERO)
        CLFR    = MAX (WAER(5) - CNACL, ZERO)
        CCACL2  = MIN (CAFR, 0.5D0*CLFR)              !SOLID CACL2
        CAFR    = MAX (CAFR - CCACL2, ZERO)
        CLFR    = MAX (WAER(5) - 2.D0*CCACL2, ZERO)
        CCANO32 = MIN (CAFR, 0.5D0*WAER(4))           !SOLID CA(NO3)2
        CAFR    = MAX (CAFR - CCANO32, ZERO)
        FRNO3   = MAX (WAER(4) - 2.D0*CCANO32, ZERO)
        CMGCL2  = MIN (FRMG, 0.5D0*CLFR)              !SOLID MGCL2
        FRMG    = MAX (FRMG - CMGCL2, ZERO)
        CLFR    = MAX (CLFR - 2.D0*CMGCL2, ZERO)
        CMGNO32 = MIN (FRMG, 0.5D0*FRNO3)             !SOLID MG(NO3)2
        FRMG    = MAX (FRMG - CMGNO32, ZERO)
        FRNO3   = MAX (FRNO3 - 2.D0*CMGNO32, ZERO)
        CNANO3  = MIN (FRNA, FRNO3)                   !SOLID NANO3
        FRNA    = MAX (FRNA - CNANO3, ZERO)
        FRNO3   = MAX (FRNO3 - CNANO3, ZERO)
        CKCL    = MIN (FRK, CLFR)                     !SOLID KCL
        FRK     = MAX (FRK - CKCL, ZERO)
        CLFR    = MAX (CLFR - CKCL, ZERO)
        CKNO3   = MIN (FRK, FRNO3)                    !SOLID KNO3
        FRK     = MAX (FRK - CKNO3, ZERO)
        FRNO3   = MAX (FRNO3 - CKNO3, ZERO)
  !
  ! *** OTHER PHASES ******************************************************
  !
        WATER   = ZERO
  !
        GNH3    = ZERO
        GHNO3   = ZERO
        GHCL    = ZERO
  !
        RETURN
  !
  ! *** END OF SUBROUTINE CALCW1A *****************************************
  !
        END SUBROUTINE CALCW1A