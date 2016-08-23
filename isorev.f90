! <isorev.f90 - A component of the EMEP MSC-W Chemical transport Model, version 3049(3049)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2015 met.no
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
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE ISRP1R
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
!     AN AMMONIUM-SULFATE AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!     THE AMBIENT RELATIVE HUMIDITY.

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE ISRP1R (WI, RHI, TEMPI)
    INCLUDE 'isrpia.inc'
    DIMENSION WI(NCOMP)

! *** INITIALIZE COMMON BLOCK VARIABLES *********************************

    CALL INIT1 (WI, RHI, TEMPI)

! *** CALCULATE SULFATE RATIO *******************************************

    IF (RH >= DRNH42S4) THEN         ! WET AEROSOL, NEED NH4 AT SRATIO=2.0
        SULRATW = GETASR(WAER(2), RHI)     ! AEROSOL SULFATE RATIO
    ELSE
        SULRATW = 2.0D0                    ! DRY AEROSOL SULFATE RATIO
    ENDIF
    SULRAT  = WAER(3)/WAER(2)         ! SULFATE RATIO

! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

! *** SULFATE POOR

    IF (SULRATW <= SULRAT) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'S2'
            CALL CALCS2                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH42S4) THEN
                SCASE = 'S1'
                CALL CALCS1              ! NH42SO4              ; case K1
            
            ELSEIF (DRNH42S4 <= RH) THEN
                SCASE = 'S2'
                CALL CALCS2              ! Only liquid          ; case K2
            ENDIF
        ENDIF
    
    ! *** SULFATE RICH (NO ACID)
    
    ELSEIF (1.0 <= SULRAT .AND. SULRAT < SULRATW) THEN
        W(2) = WAER(2)
        W(3) = WAER(3)
    
        IF(METSTBL == 1) THEN
            SCASE = 'B4'
            CALL CALCB4                 ! Only liquid (metastable)
            SCASE = 'B4'
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'B1'
                CALL CALCB1              ! NH4HSO4,LC,NH42SO4   ; case B1
                SCASE = 'B1'
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRLC) THEN
                SCASE = 'B2'
                CALL CALCB2              ! LC,NH42S4            ; case B2
                SCASE = 'B2'
            
            ELSEIF (DRLC <= RH .AND. RH < DRNH42S4) THEN
                SCASE = 'B3'
                CALL CALCB3              ! NH42S4               ; case B3
                SCASE = 'B3'
            
            ELSEIF (DRNH42S4 <= RH) THEN
                SCASE = 'B4'
                CALL CALCB4              ! Only liquid          ; case B4
                SCASE = 'B4'
            ENDIF
        ENDIF
    
        CALL CALCNH3P          ! Compute NH3(g)
    
    ! *** SULFATE RICH (FREE ACID)
    
    ELSEIF (SULRAT < 1.0) THEN
        W(2) = WAER(2)
        W(3) = WAER(3)
    
        IF(METSTBL == 1) THEN
            SCASE = 'C2'
            CALL CALCC2                 ! Only liquid (metastable)
            SCASE = 'C2'
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'C1'
                CALL CALCC1              ! NH4HSO4              ; case C1
                SCASE = 'C1'
            
            ELSEIF (DRNH4HS4 <= RH) THEN
                SCASE = 'C2'
                CALL CALCC2              ! Only liquid          ; case C2
                SCASE = 'C2'
            ENDIF
        ENDIF
    
        CALL CALCNH3P
    
    ENDIF
    RETURN

! *** END OF SUBROUTINE ISRP1R *****************************************

    END SUBROUTINE ISRP1R

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE ISRP2R
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
!     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!     THE AMBIENT RELATIVE HUMIDITY.

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE ISRP2R (WI, RHI, TEMPI)
    INCLUDE 'isrpia.inc'
    DIMENSION WI(NCOMP)
    LOGICAL ::   TRYLIQ

! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************

    TRYLIQ = .TRUE.             ! Assume liquid phase, sulfate poor limit

    10 CALL INIT2 (WI, RHI, TEMPI)

! *** CALCULATE SULFATE RATIO *******************************************

    IF (TRYLIQ .AND. RH >= DRNH4NO3) THEN ! *** WET AEROSOL
        SULRATW = GETASR(WAER(2), RHI)     ! LIMITING SULFATE RATIO
    ELSE
        SULRATW = 2.0D0                    ! *** DRY AEROSOL
    ENDIF
    SULRAT = WAER(3)/WAER(2)

! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

! *** SULFATE POOR

    IF (SULRATW <= SULRAT) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'N3'
            CALL CALCN3                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4NO3) THEN
                SCASE = 'N1'
                CALL CALCN1              ! NH42SO4,NH4NO3       ; case N1
            
            ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNH42S4) THEN
                SCASE = 'N2'
                CALL CALCN2              ! NH42S4               ; case N2
            
            ELSEIF (DRNH42S4 <= RH) THEN
                SCASE = 'N3'
                CALL CALCN3              ! Only liquid          ; case N3
            ENDIF
        ENDIF
    
    ! *** SULFATE RICH (NO ACID)
    
    !     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
    !     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE
    !     AEROSOL EQUILIBRIUM.
    
    ELSEIF (1.0 <= SULRAT .AND. SULRAT < SULRATW) THEN
        W(2) = WAER(2)
        W(3) = WAER(3)
        W(4) = WAER(4)
    
        IF(METSTBL == 1) THEN
            SCASE = 'B4'
            CALL CALCB4                 ! Only liquid (metastable)
            SCASE = 'B4'
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'B1'
                CALL CALCB1              ! NH4HSO4,LC,NH42SO4   ; case O1
                SCASE = 'B1'
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRLC) THEN
                SCASE = 'B2'
                CALL CALCB2              ! LC,NH42S4            ; case O2
                SCASE = 'B2'
            
            ELSEIF (DRLC <= RH .AND. RH < DRNH42S4) THEN
                SCASE = 'B3'
                CALL CALCB3              ! NH42S4               ; case O3
                SCASE = 'B3'
            
            ELSEIF (DRNH42S4 <= RH) THEN
                SCASE = 'B4'
                CALL CALCB4              ! Only liquid          ; case O4
                SCASE = 'B4'
            ENDIF
        ENDIF
    
    ! *** Add the NO3 to the solution now and calculate partitioning.
    
        MOLAL(7) = WAER(4)             ! There is always water, so NO3(aer) is NO3-
        MOLAL(1) = MOLAL(1) + WAER(4)  ! Add H+ to balance out
        CALL CALCNAP            ! HNO3, NH3 dissolved
        CALL CALCNH3P
    
    ! *** SULFATE RICH (FREE ACID)
    
    !     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
    !     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE
    !     AEROSOL EQUILIBRIUM.
    
    ELSEIF (SULRAT < 1.0) THEN
        W(2) = WAER(2)
        W(3) = WAER(3)
        W(4) = WAER(4)
    
        IF(METSTBL == 1) THEN
            SCASE = 'C2'
            CALL CALCC2                 ! Only liquid (metastable)
            SCASE = 'C2'
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'C1'
                CALL CALCC1              ! NH4HSO4              ; case P1
                SCASE = 'C1'
            
            ELSEIF (DRNH4HS4 <= RH) THEN
                SCASE = 'C2'
                CALL CALCC2              ! Only liquid          ; case P2
                SCASE = 'C2'
            ENDIF
        ENDIF
    
    ! *** Add the NO3 to the solution now and calculate partitioning.
    
        MOLAL(7) = WAER(4)             ! There is always water, so NO3(aer) is NO3-
        MOLAL(1) = MOLAL(1) + WAER(4)  ! Add H+ to balance out
    
        CALL CALCNAP                   ! HNO3, NH3 dissolved
        CALL CALCNH3P
    ENDIF

! *** IF SULRATW < SULRAT < 2.0 and WATER = 0 => SULFATE RICH CASE.

    IF (SULRATW <= SULRAT .AND. SULRAT < 2.0 &
     .AND. WATER <= TINY) THEN
        TRYLIQ = .FALSE. 
        GOTO 10
    ENDIF

    RETURN

! *** END OF SUBROUTINE ISRP2R *****************************************

    END SUBROUTINE ISRP2R
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE ISRP3R
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
!     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
!     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE ISRP3R (WI, RHI, TEMPI)
    INCLUDE 'isrpia.inc'
    DIMENSION WI(NCOMP)
    LOGICAL ::   TRYLIQ
! C
! C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
! C
!c      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
!c      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3

! *** INITIALIZE ALL VARIABLES ******************************************

    TRYLIQ = .TRUE.             ! Use liquid phase sulfate poor limit

    10 CALL ISOINIT3 (WI, RHI, TEMPI) ! COMMON block variables
! C
! C *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
! C
!c      REST = 2.D0*WAER(2) + WAER(4) + WAER(5)
!c      IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
!c         WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
!c         CALL PUSHERR (0050, 'ISRP3R')     ! Warning error: Na adjusted
!c      ENDIF

! *** CALCULATE SULFATE & SODIUM RATIOS *********************************

    IF (TRYLIQ .AND. RH >= DRNH4NO3) THEN  ! ** WET AEROSOL
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

! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

! *** SULFATE POOR ; SODIUM POOR

    IF (SULRATW <= SULRAT .AND. SODRAT < 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'Q5'
            CALL CALCQ5                 ! Only liquid (metastable)
            SCASE = 'Q5'
        ELSE
        
            IF (RH < DRNH4NO3) THEN
                SCASE = 'Q1'
                CALL CALCQ1              ! NH42SO4,NH4NO3,NH4CL,NA2SO4
            
            ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNH4CL) THEN
                SCASE = 'Q2'
                CALL CALCQ2              ! NH42SO4,NH4CL,NA2SO4
            
            ELSEIF (DRNH4CL <= RH  .AND. RH < DRNH42S4) THEN
                SCASE = 'Q3'
                CALL CALCQ3              ! NH42SO4,NA2SO4
            
            ELSEIF (DRNH42S4 <= RH  .AND. RH < DRNA2SO4) THEN
                SCASE = 'Q4'
                CALL CALCQ4              ! NA2SO4
                SCASE = 'Q4'
            
            ELSEIF (DRNA2SO4 <= RH) THEN
                SCASE = 'Q5'
                CALL CALCQ5              ! Only liquid
                SCASE = 'Q5'
            ENDIF
        ENDIF
    
    ! *** SULFATE POOR ; SODIUM RICH
    
    ELSE IF (SULRAT >= SULRATW .AND. SODRAT >= 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'R6'
            CALL CALCR6                 ! Only liquid (metastable)
            SCASE = 'R6'
        ELSE
        
            IF (RH < DRNH4NO3) THEN
                SCASE = 'R1'
                CALL CALCR1              ! NH4NO3,NH4CL,NA2SO4,NACL,NANO3
            
            ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNANO3) THEN
                SCASE = 'R2'
                CALL CALCR2              ! NH4CL,NA2SO4,NACL,NANO3
            
            ELSEIF (DRNANO3 <= RH  .AND. RH < DRNACL) THEN
                SCASE = 'R3'
                CALL CALCR3              ! NH4CL,NA2SO4,NACL
            
            ELSEIF (DRNACL <= RH   .AND. RH < DRNH4CL) THEN
                SCASE = 'R4'
                CALL CALCR4              ! NH4CL,NA2SO4
            
            ELSEIF (DRNH4CL <= RH .AND. RH < DRNA2SO4) THEN
                SCASE = 'R5'
                CALL CALCR5              ! NA2SO4
                SCASE = 'R5'
            
            ELSEIF (DRNA2SO4 <= RH) THEN
                SCASE = 'R6'
                CALL CALCR6              ! NO SOLID
                SCASE = 'R6'
            ENDIF
        ENDIF
    
    ! *** SULFATE RICH (NO ACID)
    
    ELSEIF (1.0 <= SULRAT .AND. SULRAT < SULRATW) THEN
        DO 100 I=1,NCOMP
            W(I) = WAER(I)
        100 END DO
    
        IF(METSTBL == 1) THEN
            SCASE = 'I6'
            CALL CALCI6                 ! Only liquid (metastable)
            SCASE = 'I6'
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'I1'
                CALL CALCI1              ! NA2SO4,(NH4)2SO4,NAHSO4,NH4HSO4,LC
                SCASE = 'I1'
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRNAHSO4) THEN
                SCASE = 'I2'
                CALL CALCI2              ! NA2SO4,(NH4)2SO4,NAHSO4,LC
                SCASE = 'I2'
            
            ELSEIF (DRNAHSO4 <= RH .AND. RH < DRLC) THEN
                SCASE = 'I3'
                CALL CALCI3              ! NA2SO4,(NH4)2SO4,LC
                SCASE = 'I3'
            
            ELSEIF (DRLC <= RH     .AND. RH < DRNH42S4) THEN
                SCASE = 'I4'
                CALL CALCI4              ! NA2SO4,(NH4)2SO4
                SCASE = 'I4'
            
            ELSEIF (DRNH42S4 <= RH .AND. RH < DRNA2SO4) THEN
                SCASE = 'I5'
                CALL CALCI5              ! NA2SO4
                SCASE = 'I5'
            
            ELSEIF (DRNA2SO4 <= RH) THEN
                SCASE = 'I6'
                CALL CALCI6              ! NO SOLIDS
                SCASE = 'I6'
            ENDIF
        ENDIF
    
        CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
        CALL CALCNH3P
    
    ! *** SULFATE RICH (FREE ACID)
    
    ELSEIF (SULRAT < 1.0) THEN
        DO 200 I=1,NCOMP
            W(I) = WAER(I)
        200 END DO
    
        IF(METSTBL == 1) THEN
            SCASE = 'J3'
            CALL CALCJ3                 ! Only liquid (metastable)
            SCASE = 'J3'
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'J1'
                CALL CALCJ1              ! NH4HSO4,NAHSO4
                SCASE = 'J1'
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRNAHSO4) THEN
                SCASE = 'J2'
                CALL CALCJ2              ! NAHSO4
                SCASE = 'J2'
            
            ELSEIF (DRNAHSO4 <= RH) THEN
                SCASE = 'J3'
                CALL CALCJ3
                SCASE = 'J3'
            ENDIF
        ENDIF
    
        CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
        CALL CALCNH3P
    
    ENDIF

! *** IF AFTER CALCULATIONS, SULRATW < SULRAT < 2.0
!                            and WATER = 0          => SULFATE RICH CASE.

    IF (SULRATW <= SULRAT .AND. SULRAT < 2.0 &
     .AND. WATER <= TINY) THEN
        TRYLIQ = .FALSE. 
        GOTO 10
    ENDIF

    RETURN

! *** END OF SUBROUTINE ISRP3R *****************************************

    END SUBROUTINE ISRP3R

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE ISRP4R
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
!     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM-CALCIUM-POTTASIUM-MAGNESIUM AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
!     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE ISRP4R (WI, RHI, TEMPI)
    INCLUDE 'isrpia.inc'
    DIMENSION WI(NCOMP)
    LOGICAL ::   TRYLIQ
! C
! C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
! C
!c      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
!c      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3

! *** INITIALIZE ALL VARIABLES ******************************************

    TRYLIQ  = .TRUE.             ! Use liquid phase sulfate poor limit
    IPROB   = 1            ! SOLVE REVERSE PROBLEM
!      METSTBL = 1

    10 CALL INIT4 (WI, RHI, TEMPI) ! COMMON block variables
! C
! C *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
! C
!c      REST = 2.D0*WAER(2) + WAER(4) + WAER(5)
!c      IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
!c         WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
!c         CALL PUSHERR (0050, 'ISRP3R')     ! Warning error: Na adjusted
!c      ENDIF

! *** CALCULATE SULFATE, CRUSTAL & SODIUM RATIOS ***********************

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

! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

! *** SULFATE POOR ; SODIUM+CRUSTALS POOR

    IF (SULRATW <= SO4RAT .AND. CRNARAT < 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'V7'
            CALL CALCV7                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4NO3) THEN
                SCASE = 'V1'
                CALL CALCV1              ! CaSO4, NH4NO3, NH4CL, (NH4)2SO4, MGSO4, NA2SO4, K2SO4
            
            ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNH4CL) THEN
                SCASE = 'V2'
                CALL CALCV2              ! CaSO4, NH4CL, (NH4)2SO4, MGSO4, NA2SO4, K2SO4
            
            ELSEIF (DRNH4CL <= RH  .AND. RH < DRNH42S4) THEN
                SCASE = 'V3'
                CALL CALCV3              ! CaSO4, (NH4)2SO4, MGSO4, NA2SO4, K2SO4
            
            ELSEIF (DRNH42S4 <= RH  .AND. RH < DRMGSO4) THEN
                SCASE = 'V4'
                CALL CALCV4              ! CaSO4, MGSO4, NA2SO4, K2SO4
            
            ELSEIF (DRMGSO4 <= RH .AND. RH < DRNA2SO4) THEN
                SCASE = 'V5'
                CALL CALCV5              ! CaSO4, NA2SO4, K2SO4
            
            ELSEIF (DRNA2SO4 <= RH .AND. RH < DRK2SO4) THEN
                SCASE = 'V6'
                CALL CALCV6              ! CaSO4, K2SO4
            
            ELSEIF (DRK2SO4 <= RH) THEN
                SCASE = 'V7'
                CALL CALCV7              ! CaSO4
            ENDIF
        ENDIF
    
    ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
    
    ELSEIF (SO4RAT >= SULRATW .AND. CRNARAT >= 2.0) THEN
    
        IF (CRRAT <= 2.0) THEN
        
            IF(METSTBL == 1) THEN
                SCASE = 'U8'
                CALL CALCU8                 ! Only liquid (metastable)
            ELSE
            
                IF (RH < DRNH4NO3) THEN
                    SCASE = 'U1'
                    CALL CALCU1             ! CaSO4, NH4NO3, NH4CL, MGSO4, NA2SO4, K2SO4, NACL, NANO3
                
                ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNANO3) THEN
                    SCASE = 'U2'
                    CALL CALCU2            ! CaSO4, NH4CL, MGSO4, NA2SO4, K2SO4, NACL, NANO3
                
                ELSEIF (DRNANO3 <= RH  .AND. RH < DRNACL) THEN
                    SCASE = 'U3'
                    CALL CALCU3            ! CaSO4, NH4CL, MGSO4, NA2SO4, K2SO4, NACL
                
                ELSEIF (DRNACL <= RH   .AND. RH < DRNH4Cl) THEN
                    SCASE = 'U4'
                    CALL CALCU4            ! CaSO4, NH4CL, MGSO4, NA2SO4, K2SO4
                
                ELSEIF (DRNH4Cl <= RH .AND. RH < DRMGSO4) THEN
                    SCASE = 'U5'
                    CALL CALCU5            ! CaSO4, MGSO4, NA2SO4, K2SO4
                
                ELSEIF (DRMGSO4 <= RH .AND. RH < DRNA2SO4) THEN
                    SCASE = 'U6'
                    CALL CALCU6            ! CaSO4, NA2SO4, K2SO4
                
                ELSEIF (DRNA2SO4 <= RH .AND. RH < DRK2SO4) THEN
                    SCASE = 'U7'
                    CALL CALCU7            ! CaSO4, K2SO4
                
                ELSEIF (DRK2SO4 <= RH) THEN
                    SCASE = 'U8'
                    CALL CALCU8            ! CaSO4
                ENDIF
            ENDIF
        
        ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
        
        ELSEIF (CRRAT > 2.0) THEN
        
            IF(METSTBL == 1) THEN
                SCASE = 'W13'
                CALL CALCW13                 ! Only liquid (metastable)
            ELSE
            
                IF (RH < DRCACL2) THEN
                    SCASE = 'W1'
                    CALL CALCW1             ! CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
                !                                    ! MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRCACL2 <= RH .AND. RH < DRMGCL2) THEN
                    SCASE = 'W2'
                    CALL CALCW2            ! CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
                !                                   ! MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRMGCL2 <= RH  .AND. RH < DRCANO32) THEN
                    SCASE = 'W3'
                    CALL CALCW3            ! CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
                !                                   ! MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRCANO32 <= RH   .AND. RH < DRMGNO32) THEN
                    SCASE = 'W4'
                    CALL CALCW4            ! CaSO4, K2SO4, KNO3, KCL, MGSO4,
                !                                   ! MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRMGNO32 <= RH .AND. RH < DRNH4NO3) THEN
                    SCASE = 'W5'
                    CALL CALCW5            ! CaSO4, K2SO4, KNO3, KCL, MGSO4,
                !                                   ! NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNANO3) THEN
                    SCASE = 'W6'
                    CALL CALCW6            ! CaSO4, K2SO4, KNO3, KCL, MGSO4, NANO3, NACL, NH4CL
                
                ELSEIF (DRNANO3 <= RH .AND. RH < DRNACL) THEN
                    SCASE = 'W7'
                    CALL CALCW7            ! CaSO4, K2SO4, KNO3, KCL, MGSO4, NACL, NH4CL
                
                ELSEIF (DRNACL <= RH .AND. RH < DRNH4CL) THEN
                    SCASE = 'W8'
                    CALL CALCW8            ! CaSO4, K2SO4, KNO3, KCL, MGSO4, NH4CL
                
                ELSEIF (DRNH4CL <= RH .AND. RH < DRKCL) THEN
                    SCASE = 'W9'
                    CALL CALCW9            ! CaSO4, K2SO4, KNO3, KCL, MGSO4
                
                ELSEIF (DRKCL <= RH .AND. RH < DRMGSO4) THEN
                    SCASE = 'W10'
                    CALL CALCW10            ! CaSO4, K2SO4, KNO3, MGSO4
                
                ELSEIF (DRMGSO4 <= RH .AND. RH < DRKNO3) THEN
                    SCASE = 'W11'
                    CALL CALCW11            ! CaSO4, K2SO4, KNO3
                
                ELSEIF (DRKNO3 <= RH .AND. RH < DRK2SO4) THEN
                    SCASE = 'W12'
                    CALL CALCW12            ! CaSO4, K2SO4
                
                ELSEIF (DRK2SO4 <= RH) THEN
                    SCASE = 'W13'
                    CALL CALCW13            ! CaSO4
                ENDIF
            ENDIF
        !        CALL CALCNH3
        ENDIF
    
    ! *** SULFATE RICH (NO ACID): 1<Rso4<2;
    
    ELSEIF (1.0 <= SO4RAT .AND. SO4RAT < SULRATW) THEN
        DO 800 I=1,NCOMP
            W(I) = WAER(I)
        800 END DO
    
        IF(METSTBL == 1) THEN
            SCASE = 'L9'
            CALL CALCL9                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'L1'
                CALL CALCL1            ! CASO4,K2SO4,MGSO4,KHSO4,NA2SO4,(NH4)2SO4,NAHSO4,NH4HSO4,LC
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRNAHSO4) THEN
                SCASE = 'L2'
                CALL CALCL2            ! CASO4,K2SO4,MGSO4,KHSO4,NA2SO4,(NH4)2SO4,NAHSO4,LC
            
            ELSEIF (DRNAHSO4 <= RH .AND. RH < DRLC) THEN
                SCASE = 'L3'
                CALL CALCL3            ! CASO4,K2SO4,MGSO4,KHSO4,NA2SO4,(NH4)2SO4,LC
            
            ELSEIF (DRLC <= RH .AND. RH < DRNH42S4) THEN
                SCASE = 'L4'
                CALL CALCL4            ! CASO4,K2SO4,MGSO4,KHSO4,NA2SO4,(NH4)2SO4
            
            ELSEIF (DRNH42S4 <= RH .AND. RH < DRKHSO4) THEN
                SCASE = 'L5'
                CALL CALCL5            ! CASO4,K2SO4,MGSO4,KHSO4,NA2SO4
            
            ELSEIF (DRKHSO4 <= RH .AND. RH < DRMGSO4) THEN
                SCASE = 'L6'
                CALL CALCL6            ! CASO4,K2SO4,MGSO4,NA2SO4
            
            ELSEIF (DRMGSO4 <= RH .AND. RH < DRNA2SO4) THEN
                SCASE = 'L7'
                CALL CALCL7            ! CASO4,K2SO4,NA2SO4
            
            ELSEIF (DRNA2SO4 <= RH .AND. RH < DRK2SO4) THEN
                SCASE = 'L8'
                CALL CALCL8            ! CASO4,K2SO4
            
            ELSEIF (DRK2SO4 <= RH) THEN
                SCASE = 'L9'
                CALL CALCL9            ! CaSO4
            ENDIF
        ENDIF
    
        CALL CALCNHP                ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3P               !                NH3
    
    ! *** SULFATE SUPER RICH (FREE ACID): Rso4<1;
    
    ELSEIF (SO4RAT < 1.0) THEN
        DO 900 I=1,NCOMP
            W(I) = WAER(I)
        900 END DO
    
        IF(METSTBL == 1) THEN
            SCASE = 'K4'
            CALL CALCK4                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4HS4) THEN                   ! RH < 0.4
                SCASE = 'K1'
                CALL CALCK1           ! NH4HSO4,NAHSO4,KHSO4,CASO4
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRNAHSO4) THEN
                SCASE = 'K2'
                CALL CALCK2           ! NAHSO4,KHSO4,CASO4
            
            ELSEIF (DRNAHSO4 <= RH .AND. RH < DRKHSO4) THEN
                SCASE = 'K3'
                CALL CALCK3           ! KHSO4,CASO4    0.52 < RH < 0.86
            
            ELSEIF (DRKHSO4 <= RH) THEN
                SCASE = 'K4'
                CALL CALCK4           ! CASO4
            ENDIF
        ENDIF
    
        CALL CALCNHP                  ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3P                 !                NH3
    
    ENDIF

! *** IF AFTER CALCULATIONS, SULRATW < SO4RAT < 2.0
!                            and WATER = 0          => SULFATE RICH CASE.

    IF (SULRATW <= SO4RAT .AND. SO4RAT < 2.0 &
     .AND. WATER <= TINY) THEN
        TRYLIQ = .FALSE. 
        GOTO 10
    ENDIF

    RETURN

! *** END OF SUBROUTINE ISRP4R *****************************************

    END SUBROUTINE ISRP4R
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCS2
! *** CASE S2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. LIQUID AEROSOL PHASE ONLY POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCS2
    INCLUDE 'isrpia.inc'
    DOUBLE PRECISION :: NH4I, NH3GI, NH3AQ

! *** SETUP PARAMETERS ************************************************

    CALAOU   = .TRUE.     ! Outer loop activity calculation flag
    FRST     = .TRUE. 
    CALAIN   = .TRUE. 

! *** CALCULATE WATER CONTENT *****************************************

    MOLALR(4)= MIN(WAER(2), 0.5d0*WAER(3))
    WATER    = MOLALR(4)/M0(4)  ! ZSR correlation

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    !C         A21  = XK21*WATER*R*TEMP
        A2   = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
        AKW  = XKW *RH*WATER*WATER
    
        NH4I = WAER(3)
        SO4I = WAER(2)
        HSO4I= ZERO
    
        CALL CALCPH (2.D0*SO4I - NH4I, HI, OHI)    ! Get pH
    
        NH3AQ = ZERO                               ! AMMONIA EQUILIBRIUM
        IF (HI < OHI) THEN
            CALL CALCAMAQ (NH4I, OHI, DEL)
            NH4I  = MAX (NH4I-DEL, ZERO)
            OHI   = MAX (OHI -DEL, TINY)
            NH3AQ = DEL
            HI    = AKW/OHI
        ENDIF
    
        CALL CALCHS4 (HI, SO4I, ZERO, DEL)         ! SULFATE EQUILIBRIUM
        SO4I  = SO4I - DEL
        HI    = HI   - DEL
        HSO4I = DEL
    
        NH3GI = NH4I/HI/A2   !    NH3AQ/A21
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
        MOLAL(1) = HI
        MOLAL(3) = NH4I
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        COH      = OHI
        GASAQ(1) = NH3AQ
        GNH3     = NH3GI
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

    20 RETURN

! *** END OF SUBROUTINE CALCS2 ****************************************

    END SUBROUTINE CALCS2
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCS1
! *** CASE S1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4

!     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE SOLID (NH4)2SO4
!     IS CALCULATED FROM THE SULFATES. THE EXCESS AMMONIA REMAINS IN
!     THE GAS PHASE.

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCS1
    INCLUDE 'isrpia.inc'

    CNH42S4 = MIN(WAER(2),0.5d0*WAER(3))  ! For bad input problems
    GNH3    = ZERO

    W(2)    = CNH42S4
    W(3)    = 2.D0*CNH42S4 + GNH3

    RETURN

! *** END OF SUBROUTINE CALCS1 ******************************************

    END SUBROUTINE CALCS1


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCN3
! *** CASE N3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. THERE IS ONLY A LIQUID PHASE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCN3
    INCLUDE 'isrpia.inc'
    DOUBLE PRECISION :: NH4I, NO3I, NH3AQ, NO3AQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** AEROSOL WATER CONTENT

    MOLALR(4) = MIN(WAER(2),0.5d0*WAER(3))       ! (NH4)2SO4
    AML5      = MAX(WAER(3)-2.D0*MOLALR(4),ZERO) ! "free" NH4
    MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO)     ! NH4NO3=MIN("free",NO3)
    WATER     = MOLALR(4)/M0(4) + MOLALR(5)/M0(5)
    WATER     = MAX(WATER, TINY)

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A2    = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
    !C         A21   = XK21*WATER*R*TEMP
        A3    = XK4*R*TEMP*(WATER/GAMA(10))**2.0
        A4    = XK7*(WATER/GAMA(4))**3.0
        AKW   = XKW *RH*WATER*WATER
    
    ! ION CONCENTRATIONS
    
        NH4I  = WAER(3)
        NO3I  = WAER(4)
        SO4I  = WAER(2)
        HSO4I = ZERO
    
        CALL CALCPH (2.D0*SO4I + NO3I - NH4I, HI, OHI)
    
    ! AMMONIA ASSOCIATION EQUILIBRIUM
    
        NH3AQ = ZERO
        NO3AQ = ZERO
        GG    = 2.D0*SO4I + NO3I - NH4I
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
        ELSE
            HI    = ZERO
            CALL CALCNIAQ2 (GG, NO3I, HI, NO3AQ) ! HNO3
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (1) = HI
        MOLAL (3) = NH4I
        MOLAL (5) = SO4I
        MOLAL (6) = HSO4I
        MOLAL (7) = NO3I
        COH       = OHI
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
    
        GASAQ(1)  = NH3AQ
        GASAQ(3)  = NO3AQ
    
        GHNO3     = HI*NO3I/A3
        GNH3      = NH4I/HI/A2   !   NH3AQ/A21
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** RETURN ***********************************************************

    20 RETURN

! *** END OF SUBROUTINE CALCN3 *****************************************

    END SUBROUTINE CALCN3
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCN2
! *** CASE N2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCN2
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    CHI1   = MIN(WAER(2),0.5d0*WAER(3))     ! (NH4)2SO4
    CHI2   = MAX(WAER(3) - 2.D0*CHI1, ZERO) ! "Free" NH4+
    CHI3   = MAX(WAER(4) - CHI2, ZERO)      ! "Free" NO3

    PSI2   = CHI2
    PSI3   = CHI3

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI1LO = TINY                ! Low  limit
    PSI1HI = CHI1                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI1HI
    Y1 = FUNCN2 (X1)
    IF (Y1 <= EPS) RETURN   ! IF (ABS(Y1) <= EPS .OR. Y1 <= ZERO) RETURN
    YHI= Y1                 ! Save Y-value at HI position

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, ZERO)
        Y2 = FUNCN2 (X2)
        IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION FOUND

    YLO= Y1                      ! Save Y-value at Hi position
    IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        RETURN
    
    ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
    
    ELSE IF (YLO < ZERO .AND. YHI < ZERO) THEN
        P4 = CHI4
        YY = FUNCN2(P4)
        GOTO 50
    
    ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
    
    ELSE IF (YLO > ZERO .AND. YHI > ZERO) THEN
        P4 = TINY
        YY = FUNCN2(P4)
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCN2')    ! WARNING ERROR: NO SOLUTION
        RETURN
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCN2 (X3)
        IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCN2')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCN2 (X3)
    50 CONTINUE
    RETURN

! *** END OF SUBROUTINE CALCN2 ******************************************

    END SUBROUTINE CALCN2



!======================================================================

! *** ISORROPIA CODE
! *** FUNCTION FUNCN2
! *** CASE D2
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D2 ;
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCN2.

!=======================================================================

    DOUBLE PRECISION FUNCTION FUNCN2 (P1)
    INCLUDE 'isrpia.inc'
    DOUBLE PRECISION :: NH4I, NO3I, NH3AQ, NO3AQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    PSI1   = P1

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A2    = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
    !C         A21   = XK21*WATER*R*TEMP
        A3    = XK4*R*TEMP*(WATER/GAMA(10))**2.0
        A4    = XK7*(WATER/GAMA(4))**3.0
        AKW   = XKW *RH*WATER*WATER
    
    ! ION CONCENTRATIONS
    
        NH4I  = 2.D0*PSI1 + PSI2
        NO3I  = PSI2 + PSI3
        SO4I  = PSI1
        HSO4I = ZERO
    
        CALL CALCPH (2.D0*SO4I + NO3I - NH4I, HI, OHI)
    
    ! AMMONIA ASSOCIATION EQUILIBRIUM
    
        NH3AQ = ZERO
        NO3AQ = ZERO
        GG    = 2.D0*SO4I + NO3I - NH4I
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
        ELSE
            HI    = ZERO
            CALL CALCNIAQ2 (GG, NO3I, HI, NO3AQ) ! HNO3
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (1) = HI
        MOLAL (3) = NH4I
        MOLAL (5) = SO4I
        MOLAL (6) = HSO4I
        MOLAL (7) = NO3I
        COH       = OHI
    
        CNH42S4   = CHI1 - PSI1
        CNH4NO3   = ZERO
    
        GASAQ(1)  = NH3AQ
        GASAQ(3)  = NO3AQ
    
        GHNO3     = HI*NO3I/A3
        GNH3      = NH4I/HI/A2   !   NH3AQ/A21
    
    ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 FUNCN2= NH4I*NH4I*SO4I/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCN2 ********************************************

  END FUNCTION FUNCN2
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCN1
! *** CASE N1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3

!     THERE ARE TWO REGIMES DEFINED BY RELATIVE HUMIDITY:
!     1. RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCN1A)
!     2. RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCN1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCN1A, CALCN2

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMASAN) THEN
        SCASE = 'N1 ; SUBCASE 1'
        CALL CALCN1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'N1 ; SUBCASE 1'
    ELSE
        SCASE = 'N1 ; SUBCASE 2'
        CALL CALCMDRP (RH, DRMASAN, DRNH4NO3, CALCN1A, CALCN2)
        SCASE = 'N1 ; SUBCASE 2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCN1 ******************************************

    END SUBROUTINE CALCN1



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCN1A
! *** CASE N1 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCN1A
    INCLUDE 'isrpia.inc'

! *** SETUP PARAMETERS *************************************************

! C      A1      = XK10/R/TEMP/R/TEMP

! *** CALCULATE AEROSOL COMPOSITION ************************************

! C      CHI1    = 2.D0*WAER(4)        ! Free parameter ; arbitrary value.
    PSI1    = WAER(4)

! *** The following statment is here to avoid negative NH4+ values in
!     CALCN? routines that call CALCN1A

    PSI2    = MAX(MIN(WAER(2),0.5d0*(WAER(3)-PSI1)),TINY)

    CNH4NO3 = PSI1
    CNH42S4 = PSI2

! C      GNH3    = CHI1 + PSI1 + 2.0*PSI2
! C      GHNO3   = A1/(CHI1-PSI1) + PSI1
    GNH3    = ZERO
    GHNO3   = ZERO

    W(2)    = PSI2
    W(3)    = GNH3  + PSI1 + 2.0*PSI2
    W(4)    = GHNO3 + PSI1

    RETURN

! *** END OF SUBROUTINE CALCN1A *****************************************

    END SUBROUTINE CALCN1A

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCQ5
! *** CASE Q5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
!     2. LIQUID AND SOLID PHASES ARE POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCQ5
    INCLUDE 'isrpia.inc'

    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCQ1A

    PSI1   = CNA2SO4      ! SALTS DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4

    CALL CALCMR           ! WATER

    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! ION CONCENTRATIONS
    
        NAI    = WAER(1)
        SO4I   = WAER(2)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
            HSO4I = ZERO
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCQ5')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = ZERO

    RETURN

! *** END OF SUBROUTINE CALCQ5 ******************************************

    END SUBROUTINE CALCQ5

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCQ4
! *** CASE Q4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
!     2. LIQUID AND SOLID PHASES ARE POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCQ4
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV1
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV1 = .TRUE. 
    PSI1O   =-GREAT
    ROOT3   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCQ1A

    CHI1   = CNA2SO4      ! SALTS

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A5  = XK5 *(WATER/GAMA(2))**3.         ! Na2SO4    <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A5) THEN
            BB =-(WAER(2) + WAER(1))
            CC = WAER(1)*WAER(2) + 0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*WAER(2) - A5)
            CALL POLY3(BB, CC, DD, ROOT3, ISLV)
            IF (ISLV /= 0) ROOT3 = TINY
            ROOT3 = MIN (ROOT3, WAER(1)/2.0, WAER(2), CHI1)
            ROOT3 = MAX (ROOT3, ZERO)
            PSI1  = CHI1-ROOT3
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        NAI = WAER(1) - 2.D0*ROOT3
        SO4I= WAER(2) - ROOT3
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
            HSO4I = ZERO
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV1) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCQ4')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = CHI1 - PSI1

    RETURN

! *** END OF SUBROUTINE CALCQ4 ******************************************

    END SUBROUTINE CALCQ4
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCQ3
! *** CASE Q3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCQ3
    INCLUDE 'isrpia.inc'
    LOGICAL :: EXNO, EXCL
    EXTERNAL CALCQ1A, CALCQ4

! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***

    EXNO = WAER(4) > TINY
    EXCL = WAER(5) > TINY

    IF (EXNO .OR. EXCL) THEN             ! *** NITRATE OR CHLORIDE EXISTS
        SCASE = 'Q3 ; SUBCASE 1'
        CALL CALCQ3A
        SCASE = 'Q3 ; SUBCASE 1'
    
    ELSE                                 ! *** NO CHLORIDE AND NITRATE
        IF (RH < DRMG3) THEN
            SCASE = 'Q3 ; SUBCASE 2'
            CALL CALCQ1A             ! SOLID
            SCASE = 'Q3 ; SUBCASE 2'
        ELSE
            SCASE = 'Q3 ; SUBCASE 3' ! MDRH (NH4)2SO4, NA2SO4
            CALL CALCMDRP (RH, DRMG3, DRNH42S4, CALCQ1A, CALCQ4)
            SCASE = 'Q3 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCQ3 ******************************************

    END SUBROUTINE CALCQ3



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCQ3A
! *** CASE Q3 ; SUBCASE A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
!     2. LIQUID AND SOLID PHASES ARE POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCQ3A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV1, PSCONV6
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV1 = .TRUE. 
    PSCONV6 = .TRUE. 

    PSI1O   =-GREAT
    PSI6O   =-GREAT

    ROOT1   = ZERO
    ROOT3   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCQ1A

    CHI1   = CNA2SO4      ! SALTS
    CHI4   = CNH4CL
    CHI6   = CNH42S4

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A5  = XK5 *(WATER/GAMA(2))**3.         ! Na2SO4    <==> Na+
        A7  = XK7 *(WATER/GAMA(4))**3.         ! (NH4)2SO4 <==> NH4+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A5) THEN
            BB =-(WAER(2) + WAER(1) - ROOT1)
            CC = WAER(1)*(WAER(2) - ROOT1) + 0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*(WAER(2) - ROOT1) - A5)
            CALL POLY3(BB, CC, DD, ROOT3, ISLV)
            IF (ISLV /= 0) ROOT3 = TINY
            ROOT3 = MIN (ROOT3, WAER(1)/2.0, WAER(2) - ROOT1, CHI1)
            ROOT3 = MAX (ROOT3, ZERO)
            PSI1  = CHI1-ROOT3
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! AMMONIUM SULFATE
    
        IF (NH4I*NH4I*SO4I > A7) THEN
            BB =-(WAER(2)+WAER(3)-ROOT3)
            CC =  WAER(3)*(WAER(2)-ROOT3+0.5D0*WAER(3))
            DD =-((WAER(2)-ROOT3)*WAER(3)**2.D0 + A7)/4.D0
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MIN(ROOT1, WAER(3), WAER(2)-ROOT3, CHI6)
            ROOT1 = MAX(ROOT1, ZERO)
            PSI6  = CHI6-ROOT1
        ENDIF
        PSCONV6 = ABS(PSI6-PSI6O) <= EPS*PSI6O
        PSI6O   = PSI6
    
    ! ION CONCENTRATIONS
    
        NAI = WAER(1) - 2.D0*ROOT3
        SO4I= WAER(2) - ROOT1 - ROOT3
        NH4I= WAER(3) - 2.D0*ROOT1
        NO3I= WAER(4)
        CLI = WAER(5)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
            HSO4I = ZERO
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV1 .AND. PSCONV6) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCQ3A')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = CHI6 - PSI6
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = CHI1 - PSI1

    RETURN

! *** END OF SUBROUTINE CALCQ3A *****************************************

    END SUBROUTINE CALCQ3A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCQ2
! *** CASE Q2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCQ2
    INCLUDE 'isrpia.inc'
    LOGICAL :: EXNO, EXCL
    EXTERNAL CALCQ1A, CALCQ3A, CALCQ4

! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***

    EXNO = WAER(4) > TINY
    EXCL = WAER(5) > TINY

    IF (EXNO) THEN                       ! *** NITRATE EXISTS
        SCASE = 'Q2 ; SUBCASE 1'
        CALL CALCQ2A
        SCASE = 'Q2 ; SUBCASE 1'
    
    ELSEIF ( .NOT. EXNO .AND. EXCL) THEN   ! *** ONLY CHLORIDE EXISTS
        IF (RH < DRMG2) THEN
            SCASE = 'Q2 ; SUBCASE 2'
            CALL CALCQ1A             ! SOLID
            SCASE = 'Q2 ; SUBCASE 2'
        ELSE
            SCASE = 'Q2 ; SUBCASE 3' ! MDRH (NH4)2SO4, NA2SO4, NH4CL
            CALL CALCMDRP (RH, DRMG2, DRNH4CL, CALCQ1A, CALCQ3A)
            SCASE = 'Q2 ; SUBCASE 3'
        ENDIF
    
    ELSE                                 ! *** NO CHLORIDE AND NITRATE
        IF (RH < DRMG3) THEN
            SCASE = 'Q2 ; SUBCASE 2'
            CALL CALCQ1A             ! SOLID
            SCASE = 'Q2 ; SUBCASE 2'
        ELSE
            SCASE = 'Q2 ; SUBCASE 4' ! MDRH (NH4)2SO4, NA2SO4
            CALL CALCMDRP (RH, DRMG3, DRNH42S4, CALCQ1A, CALCQ4)
            SCASE = 'Q2 ; SUBCASE 4'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCQ2 ******************************************

    END SUBROUTINE CALCQ2


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCQ2A
! *** CASE Q2 ; SUBCASE A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
!     2. LIQUID AND SOLID PHASES ARE POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCQ2A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV1, PSCONV4, PSCONV6
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV1 = .TRUE. 
    PSCONV4 = .TRUE. 
    PSCONV6 = .TRUE. 

    PSI1O   =-GREAT
    PSI4O   =-GREAT
    PSI6O   =-GREAT

    ROOT1   = ZERO
    ROOT2   = ZERO
    ROOT3   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCQ1A

    CHI1   = CNA2SO4      ! SALTS
    CHI4   = CNH4CL
    CHI6   = CNH42S4

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A5  = XK5 *(WATER/GAMA(2))**3.         ! Na2SO4    <==> Na+
        A14 = XK14*(WATER/GAMA(6))**2.         ! NH4Cl     <==> NH4+
        A7  = XK7 *(WATER/GAMA(4))**3.         ! (NH4)2SO4 <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! AMMONIUM CHLORIDE
    
        IF (NH4I*CLI > A14) THEN
            BB    =-(WAER(3) + WAER(5) - 2.D0*ROOT1)
            CC    = WAER(5)*(WAER(3) - 2.D0*ROOT1) - A14
            DD    = BB*BB - 4.D0*CC
            IF (DD < ZERO) THEN
                ROOT2 = ZERO
            ELSE
                DD    = SQRT(DD)
                ROOT2A= 0.5D0*(-BB+DD)
                ROOT2B= 0.5D0*(-BB-DD)
                IF (ZERO <= ROOT2A) THEN
                    ROOT2 = ROOT2A
                ELSE
                    ROOT2 = ROOT2B
                ENDIF
                ROOT2 = MIN(ROOT2, WAER(5), WAER(3) - 2.D0*ROOT1, CHI4)
                ROOT2 = MAX(ROOT2, ZERO)
                PSI4  = CHI4 - ROOT2
            ENDIF
        ENDIF
        PSCONV4 = ABS(PSI4-PSI4O) <= EPS*PSI4O
        PSI4O   = PSI4
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A5) THEN
            BB =-(WAER(2) + WAER(1) - ROOT1)
            CC = WAER(1)*(WAER(2) - ROOT1) + 0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*(WAER(2) - ROOT1) - A5)
            CALL POLY3(BB, CC, DD, ROOT3, ISLV)
            IF (ISLV /= 0) ROOT3 = TINY
            ROOT3 = MIN (ROOT3, WAER(1)/2.0, WAER(2) - ROOT1, CHI1)
            ROOT3 = MAX (ROOT3, ZERO)
            PSI1  = CHI1-ROOT3
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! AMMONIUM SULFATE
    
        IF (NH4I*NH4I*SO4I > A7) THEN
            BB =-(WAER(2)+WAER(3)-ROOT2-ROOT3)
            CC = (WAER(3)-ROOT2)*(WAER(2)-ROOT3+0.5D0*(WAER(3)-ROOT2))
            DD =-((WAER(2)-ROOT3)*(WAER(3)-ROOT2)**2.D0 + A7)/4.D0
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MIN(ROOT1, WAER(3)-ROOT2, WAER(2)-ROOT3, CHI6)
            ROOT1 = MAX(ROOT1, ZERO)
            PSI6  = CHI6-ROOT1
        ENDIF
        PSCONV6 = ABS(PSI6-PSI6O) <= EPS*PSI6O
        PSI6O   = PSI6
    
    ! ION CONCENTRATIONS
    
        NAI = WAER(1) - 2.D0*ROOT3
        SO4I= WAER(2) - ROOT1 - ROOT3
        NH4I= WAER(3) - ROOT2 - 2.D0*ROOT1
        NO3I= WAER(4)
        CLI = WAER(5) - ROOT2
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
            HSO4I = ZERO
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV1 .AND. PSCONV4 .AND. PSCONV6) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCQ2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = CHI6 - PSI6
    CNH4NO3 = ZERO
    CNH4CL  = CHI4 - PSI4
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = CHI1 - PSI1

    RETURN

! *** END OF SUBROUTINE CALCQ2A *****************************************

    END SUBROUTINE CALCQ2A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCQ1
! *** CASE Q1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, (NH4)2SO4, NA2SO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCQ1
    INCLUDE 'isrpia.inc'
    LOGICAL :: EXNO, EXCL
    EXTERNAL CALCQ1A, CALCQ2A, CALCQ3A, CALCQ4

! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***

    EXNO = WAER(4) > TINY
    EXCL = WAER(5) > TINY

    IF (EXNO .AND. EXCL) THEN           ! *** NITRATE & CHLORIDE EXIST
        IF (RH < DRMG1) THEN
            SCASE = 'Q1 ; SUBCASE 1'
            CALL CALCQ1A             ! SOLID
            SCASE = 'Q1 ; SUBCASE 1'
        ELSE
            SCASE = 'Q1 ; SUBCASE 2' ! MDRH (NH4)2SO4, NA2SO4, NH4CL, NH4NO3
            CALL CALCMDRP (RH, DRMG1, DRNH4NO3, CALCQ1A, CALCQ2A)
            SCASE = 'Q1 ; SUBCASE 2'
        ENDIF
    
    ELSE IF (EXNO .AND. .NOT. EXCL) THEN ! *** ONLY NITRATE EXISTS
        IF (RH < DRMQ1) THEN
            SCASE = 'Q1 ; SUBCASE 1'
            CALL CALCQ1A             ! SOLID
            SCASE = 'Q1 ; SUBCASE 1'
        ELSE
            SCASE = 'Q1 ; SUBCASE 3' ! MDRH (NH4)2SO4, NA2SO4, NH4NO3
            CALL CALCMDRP (RH, DRMQ1, DRNH4NO3, CALCQ1A, CALCQ2A)
            SCASE = 'Q1 ; SUBCASE 3'
        ENDIF
    
    ELSE IF ( .NOT. EXNO .AND. EXCL) THEN ! *** ONLY CHLORIDE EXISTS
        IF (RH < DRMG2) THEN
            SCASE = 'Q1 ; SUBCASE 1'
            CALL CALCQ1A             ! SOLID
            SCASE = 'Q1 ; SUBCASE 1'
        ELSE
            SCASE = 'Q1 ; SUBCASE 4' ! MDRH (NH4)2SO4, NA2SO4, NH4CL
            CALL CALCMDRP (RH, DRMG2, DRNH4CL, CALCQ1A, CALCQ3A)
            SCASE = 'Q1 ; SUBCASE 4'
        ENDIF
    
    ELSE                                ! *** NO CHLORIDE AND NITRATE
        IF (RH < DRMG3) THEN
            SCASE = 'Q1 ; SUBCASE 1'
            CALL CALCQ1A             ! SOLID
            SCASE = 'Q1 ; SUBCASE 1'
        ELSE
            SCASE = 'Q1 ; SUBCASE 5' ! MDRH (NH4)2SO4, NA2SO4
            CALL CALCMDRP (RH, DRMG3, DRNH42S4, CALCQ1A, CALCQ4)
            SCASE = 'Q1 ; SUBCASE 5'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCQ1 ******************************************

    END SUBROUTINE CALCQ1


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCQ1A
! *** CASE Q1 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, (NH4)2SO4, NA2SO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCQ1A
    INCLUDE 'isrpia.inc'

! *** CALCULATE SOLIDS **************************************************

    CNA2SO4 = 0.5d0*WAER(1)
    FRSO4   = MAX (WAER(2)-CNA2SO4, ZERO)

    CNH42S4 = MAX (MIN(FRSO4,0.5d0*WAER(3)), TINY)
    FRNH3   = MAX (WAER(3)-2.D0*CNH42S4, ZERO)

    CNH4NO3 = MIN (FRNH3, WAER(4))
! C      FRNO3   = MAX (WAER(4)-CNH4NO3, ZERO)
    FRNH3   = MAX (FRNH3-CNH4NO3, ZERO)

    CNH4CL  = MIN (FRNH3, WAER(5))
! C      FRCL    = MAX (WAER(5)-CNH4CL, ZERO)
    FRNH3   = MAX (FRNH3-CNH4CL, ZERO)

! *** OTHER PHASES ******************************************************

    WATER   = ZERO

    GNH3    = ZERO
    GHNO3   = ZERO
    GHCL    = ZERO

    RETURN

! *** END OF SUBROUTINE CALCQ1A *****************************************

    END SUBROUTINE CALCQ1A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR6
! *** CASE R6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS ONLY A LIQUID PHASE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR6
    INCLUDE 'isrpia.inc'

    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    CALL CALCR1A

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    CALAOU = .TRUE. 

! *** CALCULATE WATER **************************************************

    CALL CALCMR

! *** SETUP LIQUID CONCENTRATIONS **************************************

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+
    
        NAI    = WAER(1)
        SO4I   = WAER(2)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG  = 2.D0*WAER(2) + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCR6')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3     = NH4I/HI/A2
    GHNO3    = HI*NO3I/A3
    GHCL     = HI*CLI /A4

    GASAQ(1) = NH3AQ
    GASAQ(2) = CLAQ
    GASAQ(3) = NO3AQ

    CNH42S4  = ZERO
    CNH4NO3  = ZERO
    CNH4CL   = ZERO
    CNACL    = ZERO
    CNANO3   = ZERO
    CNA2SO4  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCR6 ******************************************

    END SUBROUTINE CALCR6
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR5
! *** CASE R5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!     2. LIQUID AND SOLID PHASES ARE POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR5
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

    LOGICAL ::  NEAN, NEAC, NESN, NESC

! *** SETUP PARAMETERS ************************************************

    CALL CALCR1A                             ! DRY SOLUTION

    NEAN = CNH4NO3 <= TINY    ! NH4NO3       ! Water exists?
    NEAC = CNH4CL <= TINY    ! NH4CL
    NESN = CNANO3 <= TINY    ! NANO3
    NESC = CNACL  <= TINY    ! NACL
    IF (NEAN .AND. NEAC .AND. NESN .AND. NESC) RETURN

    CHI1   = CNA2SO4

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3

    PSIO   =-GREAT

! *** CALCULATE WATER **************************************************

    CALL CALCMR

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    CALAOU = .TRUE. 
    PSCONV = .FALSE. 

! *** SETUP LIQUID CONCENTRATIONS **************************************

    NAI    = WAER(1)
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A5  = XK5*(WATER/GAMA(2))**3.                   ! Na2SO4 <==> Na+
        AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+
    
    ! SODIUM SULFATE
    
        ROOT = ZERO
        IF (NAI*NAI*SO4I > A5) THEN
            BB =-3.D0*CHI1
            CC = 3.D0*CHI1**2.0
            DD =-CHI1**3.0 + 0.25D0*A5
            CALL POLY3(BB, CC, DD, ROOT, ISLV)
            IF (ISLV /= 0) ROOT = TINY
            ROOT = MIN (MAX(ROOT,ZERO), CHI1)
            PSI1 = CHI1-ROOT
        ENDIF
        PSCONV = ABS(PSI1-PSIO) <= EPS*PSIO
        PSIO   = PSI1
    
    ! ION CONCENTRATIONS
    
        NAI  = WAER(1) - 2.D0*ROOT
        SO4I = WAER(2) - ROOT
        NH4I = WAER(3)
        NO3I = WAER(4)
        CLI  = WAER(5)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCR5')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
!C      A21      = XK21*WATER*R*TEMP
    A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3     = NH4I/HI/A2  ! NH4I*OHI/A2/AKW
    GHNO3    = HI*NO3I/A3
    GHCL     = HI*CLI /A4

    GASAQ(1) = NH3AQ
    GASAQ(2) = CLAQ
    GASAQ(3) = NO3AQ

    CNH42S4  = ZERO
    CNH4NO3  = ZERO
    CNH4CL   = ZERO
    CNACL    = ZERO
    CNANO3   = ZERO
    CNA2SO4  = CHI1 - PSI1

    RETURN

! *** END OF SUBROUTINE CALCR5 ******************************************

    END SUBROUTINE CALCR5
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR4
! *** CASE R4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR4
    INCLUDE 'isrpia.inc'
    LOGICAL ::  EXAN, EXAC, EXSN, EXSC
    EXTERNAL CALCR1A, CALCR5

! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************

    SCASE = 'R4 ; SUBCASE 2'
    CALL CALCR1A              ! SOLID
    SCASE = 'R4 ; SUBCASE 2'

    EXAN = CNH4NO3 > TINY    ! NH4NO3
    EXAC = CNH4CL > TINY    ! NH4CL
    EXSN = CNANO3 > TINY    ! NANO3
    EXSC = CNACL  > TINY    ! NACL

! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********

    IF (EXAN .OR. EXSN .OR. EXSC) THEN   ! *** NH4NO3,NANO3 EXIST
        IF (RH >= DRMH1) THEN
            SCASE = 'R4 ; SUBCASE 1'
            CALL CALCR4A
            SCASE = 'R4 ; SUBCASE 1'
        ENDIF
    
    ELSE IF (EXAC) THEN                  ! *** NH4CL EXISTS ONLY
        IF (RH >= DRMR5) THEN
            SCASE = 'R4 ; SUBCASE 3'
            CALL CALCMDRP (RH, DRMR5, DRNH4CL, CALCR1A, CALCR5)
            SCASE = 'R4 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCR4 ******************************************

    END SUBROUTINE CALCR4



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR4A
! *** CASE R4A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!     2. LIQUID AND SOLID PHASES ARE POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR4A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV1, PSCONV4
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 
    PSCONV1 = .FALSE. 
    PSCONV4 = .FALSE. 
    PSIO1   =-GREAT
    PSIO4   =-GREAT

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCR1A

    CHI1   = CNA2SO4      ! SALTS
    CHI4   = CNH4CL

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A5  = XK5 *(WATER/GAMA(2))**3.                  ! Na2SO4 <==> Na+
        A14 = XK14*(WATER/GAMA(6))**2.                  ! NH4Cl  <==> NH4+
        AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+
    
    ! SODIUM SULFATE
    
        ROOT = ZERO
        IF (NAI*NAI*SO4I > A5) THEN
            BB =-3.D0*CHI1
            CC = 3.D0*CHI1**2.0
            DD =-CHI1**3.0 + 0.25D0*A5
            CALL POLY3(BB, CC, DD, ROOT, ISLV)
            IF (ISLV /= 0) ROOT = TINY
            ROOT = MIN (MAX(ROOT,ZERO), CHI1)
            PSI1 = CHI1-ROOT
            NAI  = WAER(1) - 2.D0*ROOT
            SO4I = WAER(2) - ROOT
        ENDIF
        PSCONV1 = ABS(PSI1-PSIO1) <= EPS*PSIO1
        PSIO1   = PSI1
    
    ! AMMONIUM CHLORIDE
    
        ROOT = ZERO
        IF (NH4I*CLI > A14) THEN
            BB   =-(NH4I + CLI)
            CC   =-A14 + NH4I*CLI
            DD   = BB*BB - 4.D0*CC
            ROOT = 0.5D0*(-BB-SQRT(DD))
            IF (ROOT > TINY) THEN
                ROOT    = MIN(ROOT, CHI4)
                PSI4    = CHI4 - ROOT
                NH4I    = WAER(3) - ROOT
                CLI     = WAER(5) - ROOT
            ENDIF
        ENDIF
        PSCONV4 = ABS(PSI4-PSIO4) <= EPS*PSIO4
        PSIO4   = PSI4
    
        NO3I   = WAER(4)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV1 .AND. PSCONV4) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCR4A')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI4 - PSI4
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = CHI1 - PSI1

    RETURN

! *** END OF SUBROUTINE CALCR4A *****************************************

    END SUBROUTINE CALCR4A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR3
! *** CASE R3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR3
    INCLUDE 'isrpia.inc'
    LOGICAL ::  EXAN, EXAC, EXSN, EXSC
    EXTERNAL CALCR1A, CALCR4A, CALCR5

! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************

    SCASE = 'R3 ; SUBCASE 2'
    CALL CALCR1A              ! SOLID
    SCASE = 'R3 ; SUBCASE 2'

    EXAN = CNH4NO3 > TINY    ! NH4NO3
    EXAC = CNH4CL > TINY    ! NH4CL
    EXSN = CNANO3 > TINY    ! NANO3
    EXSC = CNACL  > TINY    ! NACL

! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********

    IF (EXAN .OR. EXSN) THEN                   ! *** NH4NO3,NANO3 EXIST
        IF (RH >= DRMH1) THEN
            SCASE = 'R3 ; SUBCASE 1'
            CALL CALCR3A
            SCASE = 'R3 ; SUBCASE 1'
        ENDIF
    
    ELSE IF ( .NOT. EXAN .AND. .NOT. EXSN) THEN   ! *** NH4NO3,NANO3 = 0
        IF      (     EXAC .AND.      EXSC) THEN
            IF (RH >= DRMR4) THEN
                SCASE = 'R3 ; SUBCASE 3'
                CALL CALCMDRP (RH, DRMR4, DRNACL, CALCR1A, CALCR4A)
                SCASE = 'R3 ; SUBCASE 3'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSC) THEN
            IF (RH >= DRMR2) THEN
                SCASE = 'R3 ; SUBCASE 4'
                CALL CALCMDRP (RH, DRMR2, DRNACL, CALCR1A, CALCR4A)
                SCASE = 'R3 ; SUBCASE 4'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR5) THEN
                SCASE = 'R3 ; SUBCASE 5'
                CALL CALCMDRP (RH, DRMR5, DRNACL, CALCR1A, CALCR5)
                SCASE = 'R3 ; SUBCASE 5'
            ENDIF
        ENDIF
    
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCR3 ******************************************

    END SUBROUTINE CALCR3


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR3A
! *** CASE R3A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!     2. LIQUID AND SOLID PHASES ARE POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR3A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV1, PSCONV3, PSCONV4
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 
    PSCONV1 = .TRUE. 
    PSCONV3 = .TRUE. 
    PSCONV4 = .TRUE. 
    PSI1O   =-GREAT
    PSI3O   =-GREAT
    PSI4O   =-GREAT
    ROOT1   = ZERO
    ROOT2   = ZERO
    ROOT3   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCR1A

    CHI1   = CNA2SO4      ! SALTS
    CHI4   = CNH4CL
    CHI3   = CNACL

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

    MOLAL(1) = ZERO
    MOLAL(2) = NAI
    MOLAL(3) = NH4I
    MOLAL(4) = CLI
    MOLAL(5) = SO4I
    MOLAL(6) = HSO4I
    MOLAL(7) = NO3I

    CALL CALCACT          ! CALCULATE ACTIVITY COEFFICIENTS

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A5  = XK5 *(WATER/GAMA(2))**3.                  ! Na2SO4 <==> Na+
        A8  = XK8 *(WATER/GAMA(1))**2.                  ! NaCl   <==> Na+
        A14 = XK14*(WATER/GAMA(6))**2.                  ! NH4Cl  <==> NH4+
        AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+
    
    ! AMMONIUM CHLORIDE
    
        IF (NH4I*CLI > A14) THEN
            BB    =-(WAER(3) + WAER(5) - ROOT3)
            CC    =-A14 + NH4I*(WAER(5) - ROOT3)
            DD    = MAX(BB*BB - 4.D0*CC, ZERO)
            ROOT2A= 0.5D0*(-BB+SQRT(DD))
            ROOT2B= 0.5D0*(-BB-SQRT(DD))
            IF (ZERO <= ROOT2A) THEN
                ROOT2 = ROOT2A
            ELSE
                ROOT2 = ROOT2B
            ENDIF
            ROOT2 = MIN(MAX(ZERO, ROOT2), MAX(WAER(5)-ROOT3,ZERO), &
            CHI4, WAER(3))
            PSI4  = CHI4 - ROOT2
        ENDIF
        PSCONV4 = ABS(PSI4-PSI4O) <= EPS*PSI4O
        PSI4O   = PSI4
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A5) THEN
            BB =-(CHI1 + WAER(1) - ROOT3)
            CC = 0.25D0*(WAER(1) - ROOT3)*(4.D0*CHI1+WAER(1)-ROOT3)
            DD =-0.25D0*(CHI1*(WAER(1)-ROOT3)**2.D0 - A5)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MIN (MAX(ROOT1,ZERO), MAX(WAER(1)-ROOT3,ZERO), &
            CHI1, WAER(2))
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! ION CONCENTRATIONS
    
        NAI = WAER(1) - (2.D0*ROOT1 + ROOT3)
        SO4I= WAER(2) - ROOT1
        NH4I= WAER(3) - ROOT2
        CLI = WAER(5) - (ROOT3 + ROOT2)
        NO3I= WAER(4)
    
    ! SODIUM CHLORIDE  ; To obtain new value for ROOT3
    
        IF (NAI*CLI > A8) THEN
            BB    =-((CHI1-2.D0*ROOT1) + (WAER(5) - ROOT2))
            CC    = (CHI1-2.D0*ROOT1)*(WAER(5) - ROOT2) - A8
            DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT3A= 0.5D0*(-BB-SQRT(DD))
            ROOT3B= 0.5D0*(-BB+SQRT(DD))
            IF (ZERO <= ROOT3A) THEN
                ROOT3 = ROOT3A
            ELSE
                ROOT3 = ROOT3B
            ENDIF
            ROOT3   = MIN(MAX(ROOT3, ZERO), CHI3)
            PSI3    = CHI3-ROOT3
        ENDIF
        PSCONV3 = ABS(PSI3-PSI3O) <= EPS*PSI3O
        PSI3O   = PSI3
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV1 .AND. PSCONV3 .AND. PSCONV4) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCR3A')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 IF (CLI <= TINY .AND. WAER(5) > TINY) THEN !No disslv Cl-;solid only
        DO 30 I=1,NIONS
            MOLAL(I) = ZERO
        30 END DO
        DO 40 I=1,NGASAQ
            GASAQ(I) = ZERO
        40 END DO
        CALL CALCR1A
    ELSE
        A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
    
        GNH3    = NH4I/HI/A2
        GHNO3   = HI*NO3I/A3
        GHCL    = HI*CLI /A4
    
        GASAQ(1)= NH3AQ
        GASAQ(2)= CLAQ
        GASAQ(3)= NO3AQ
    
        CNH42S4 = ZERO
        CNH4NO3 = ZERO
        CNH4CL  = CHI4 - PSI4
        CNACL   = CHI3 - PSI3
        CNANO3  = ZERO
        CNA2SO4 = CHI1 - PSI1
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCR3A *****************************************

    END SUBROUTINE CALCR3A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR2
! *** CASE R2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR2
    INCLUDE 'isrpia.inc'
    LOGICAL ::  EXAN, EXAC, EXSN, EXSC
    EXTERNAL CALCR1A, CALCR3A, CALCR4A, CALCR5

! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************

    SCASE = 'R2 ; SUBCASE 2'
    CALL CALCR1A              ! SOLID
    SCASE = 'R2 ; SUBCASE 2'

    EXAN = CNH4NO3 > TINY    ! NH4NO3
    EXAC = CNH4CL > TINY    ! NH4CL
    EXSN = CNANO3 > TINY    ! NANO3
    EXSC = CNACL  > TINY    ! NACL

! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********

    IF (EXAN) THEN                             ! *** NH4NO3 EXISTS
        IF (RH >= DRMH1) THEN
            SCASE = 'R2 ; SUBCASE 1'
            CALL CALCR2A
            SCASE = 'R2 ; SUBCASE 1'
        ENDIF
    
    ELSE IF ( .NOT. EXAN) THEN                   ! *** NH4NO3 = 0
        IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMH2) THEN
                SCASE = 'R2 ; SUBCASE 3'
                CALL CALCMDRP (RH, DRMH2, DRNANO3, CALCR1A, CALCR3A)
                SCASE = 'R2 ; SUBCASE 3'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR1) THEN
                SCASE = 'R2 ; SUBCASE 4'
                CALL CALCMDRP (RH, DRMR1, DRNANO3, CALCR1A, CALCR3A)
                SCASE = 'R2 ; SUBCASE 4'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR2) THEN
                SCASE = 'R2 ; SUBCASE 5'
                CALL CALCMDRP (RH, DRMR2, DRNACL, CALCR1A, CALCR4A)
                SCASE = 'R2 ; SUBCASE 5'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR3) THEN
                SCASE = 'R2 ; SUBCASE 6'
                CALL CALCMDRP (RH, DRMR3, DRNANO3, CALCR1A, CALCR3A)
                SCASE = 'R2 ; SUBCASE 6'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR4) THEN
                SCASE = 'R2 ; SUBCASE 7'
                CALL CALCMDRP (RH, DRMR4, DRNACL, CALCR1A, CALCR4A)
                SCASE = 'R2 ; SUBCASE 7'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR5) THEN
                SCASE = 'R2 ; SUBCASE 8'
                CALL CALCMDRP (RH, DRMR5, DRNH4CL, CALCR1A, CALCR5)
                SCASE = 'R2 ; SUBCASE 8'
            ENDIF

        ELSE IF (     EXAC .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR6) THEN
                SCASE = 'R2 ; SUBCASE 9'
                CALL CALCMDRP (RH, DRMR6, DRNANO3, CALCR1A, CALCR3A)
                SCASE = 'R2 ; SUBCASE 9'
            ENDIF
        ENDIF
    
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCR2 ******************************************

    END SUBROUTINE CALCR2


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR2A
! *** CASE R2A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!     2. LIQUID AND SOLID PHASES ARE POSSIBLE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR2A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV1, PSCONV2, PSCONV3, PSCONV4
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV1 = .TRUE. 
    PSCONV2 = .TRUE. 
    PSCONV3 = .TRUE. 
    PSCONV4 = .TRUE. 

    PSI1O   =-GREAT
    PSI2O   =-GREAT
    PSI3O   =-GREAT
    PSI4O   =-GREAT

    ROOT1   = ZERO
    ROOT2   = ZERO
    ROOT3   = ZERO
    ROOT4   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCR1A

    CHI1   = CNA2SO4      ! SALTS
    CHI2   = CNANO3
    CHI3   = CNACL
    CHI4   = CNH4CL

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

    MOLAL(1) = ZERO
    MOLAL(2) = NAI
    MOLAL(3) = NH4I
    MOLAL(4) = CLI
    MOLAL(5) = SO4I
    MOLAL(6) = HSO4I
    MOLAL(7) = NO3I

    CALL CALCACT          ! CALCULATE ACTIVITY COEFFICIENTS

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A5  = XK5 *(WATER/GAMA(2))**3.                  ! Na2SO4 <==> Na+
        A8  = XK8 *(WATER/GAMA(1))**2.                  ! NaCl   <==> Na+
        A9  = XK9 *(WATER/GAMA(3))**2.                  ! NaNO3  <==> Na+
        A14 = XK14*(WATER/GAMA(6))**2.                  ! NH4Cl  <==> NH4+
        AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+
    
    ! AMMONIUM CHLORIDE
    
        IF (NH4I*CLI > A14) THEN
            BB    =-(WAER(3) + WAER(5) - ROOT3)
            CC    = NH4I*(WAER(5) - ROOT3) - A14
            DD    = MAX(BB*BB - 4.D0*CC, ZERO)
            DD    = SQRT(DD)
            ROOT2A= 0.5D0*(-BB+DD)
            ROOT2B= 0.5D0*(-BB-DD)
            IF (ZERO <= ROOT2A) THEN
                ROOT2 = ROOT2A
            ELSE
                ROOT2 = ROOT2B
            ENDIF
            ROOT2 = MIN(MAX(ROOT2, ZERO), CHI4)
            PSI4  = CHI4 - ROOT2
        ENDIF
        PSCONV4 = ABS(PSI4-PSI4O) <= EPS*PSI4O
        PSI4O   = PSI4
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A5) THEN
            BB =-(WAER(2) + WAER(1) - ROOT3 - ROOT4)
            CC = WAER(1)*(2.D0*ROOT3 + 2.D0*ROOT4 - 4.D0*WAER(2) - ONE) &
            -(ROOT3 + ROOT4)**2.0 + 4.D0*WAER(2)*(ROOT3 + ROOT4)
            CC =-0.25*CC
            DD = WAER(1)*WAER(2)*(ONE - 2.D0*ROOT3 - 2.D0*ROOT4) + &
            WAER(2)*(ROOT3 + ROOT4)**2.0 - A5
            DD =-0.25*DD
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MIN (MAX(ROOT1,ZERO), CHI1)
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! SODIUM NITRATE
    
        IF (NAI*NO3I > A9) THEN
            BB    =-(WAER(4) + WAER(1) - 2.D0*ROOT1 - ROOT3)
            CC    = WAER(4)*(WAER(1) - 2.D0*ROOT1 - ROOT3) - A9
            DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT4A= 0.5D0*(-BB-DD)
            ROOT4B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT4A) THEN
                ROOT4 = ROOT4A
            ELSE
                ROOT4 = ROOT4B
            ENDIF
            ROOT4 = MIN(MAX(ROOT4, ZERO), CHI2)
            PSI2  = CHI2-ROOT4
        ENDIF
        PSCONV2 = ABS(PSI2-PSI2O) <= EPS*PSI2O
        PSI2O   = PSI2
    
    ! ION CONCENTRATIONS
    
        NAI = WAER(1) - (2.D0*ROOT1 + ROOT3 + ROOT4)
        SO4I= WAER(2) - ROOT1
        NH4I= WAER(3) - ROOT2
        NO3I= WAER(4) - ROOT4
        CLI = WAER(5) - (ROOT3 + ROOT2)
    
    ! SODIUM CHLORIDE  ; To obtain new value for ROOT3
    
        IF (NAI*CLI > A8) THEN
            BB    =-(WAER(1) - 2.D0*ROOT1 + WAER(5) - ROOT2 - ROOT4)
            CC    = (WAER(5) + ROOT2)*(WAER(1) - 2.D0*ROOT1 - ROOT4) - A8
            DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT3A= 0.5D0*(-BB-DD)
            ROOT3B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT3A) THEN
                ROOT3 = ROOT3A
            ELSE
                ROOT3 = ROOT3B
            ENDIF
            ROOT3   = MIN(MAX(ROOT3, ZERO), CHI3)
            PSI3    = CHI3-ROOT3
        ENDIF
        PSCONV3 = ABS(PSI3-PSI3O) <= EPS*PSI3O
        PSI3O   = PSI3
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI < OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
        ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL > TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3 > TINY) THEN
                IF (GGCL <= TINY) HI = ZERO
                CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL(1) = HI
        MOLAL(2) = NAI
        MOLAL(3) = NH4I
        MOLAL(4) = CLI
        MOLAL(5) = SO4I
        MOLAL(6) = HSO4I
        MOLAL(7) = NO3I
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV1 .AND. PSCONV2 .AND. PSCONV3 .AND. PSCONV4) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCR2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 IF (CLI <= TINY .AND. WAER(5) > TINY) THEN !No disslv Cl-;solid only
        DO 30 I=1,NIONS
            MOLAL(I) = ZERO
        30 END DO
        DO 40 I=1,NGASAQ
            GASAQ(I) = ZERO
        40 END DO
        CALL CALCR1A
    ELSE                                     ! OK, aqueous phase present
        A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
    
        GNH3    = NH4I/HI/A2
        GHNO3   = HI*NO3I/A3
        GHCL    = HI*CLI /A4
    
        GASAQ(1)= NH3AQ
        GASAQ(2)= CLAQ
        GASAQ(3)= NO3AQ
    
        CNH42S4 = ZERO
        CNH4NO3 = ZERO
        CNH4CL  = CHI4 - PSI4
        CNACL   = CHI3 - PSI3
        CNANO3  = CHI2 - PSI2
        CNA2SO4 = CHI1 - PSI1
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCR2A *****************************************

    END SUBROUTINE CALCR2A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR1
! *** CASE R1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR1
    INCLUDE 'isrpia.inc'
    LOGICAL ::  EXAN, EXAC, EXSN, EXSC
    EXTERNAL CALCR1A, CALCR2A, CALCR3A, CALCR4A, CALCR5

! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************

    SCASE = 'R1 ; SUBCASE 1'
    CALL CALCR1A              ! SOLID
    SCASE = 'R1 ; SUBCASE 1'

    EXAN = CNH4NO3 > TINY    ! NH4NO3
    EXAC = CNH4CL > TINY    ! NH4CL
    EXSN = CNANO3 > TINY    ! NANO3
    EXSC = CNACL  > TINY    ! NACL

! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********

    IF (EXAN .AND. EXAC .AND. EXSC .AND. EXSN) THEN  ! *** ALL EXIST
        IF (RH >= DRMH1) THEN
            SCASE = 'R1 ; SUBCASE 2'  ! MDRH
            CALL CALCMDRP (RH, DRMH1, DRNH4NO3, CALCR1A, CALCR2A)
            SCASE = 'R1 ; SUBCASE 2'
        ENDIF
    
    ELSE IF ( .NOT. EXAN) THEN                   ! *** NH4NO3 = 0
        IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMH2) THEN
                SCASE = 'R1 ; SUBCASE 3'
                CALL CALCMDRP (RH, DRMH2, DRNANO3, CALCR1A, CALCR3A)
                SCASE = 'R1 ; SUBCASE 3'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR1) THEN
                SCASE = 'R1 ; SUBCASE 4'
                CALL CALCMDRP (RH, DRMR1, DRNANO3, CALCR1A, CALCR3A)
                SCASE = 'R1 ; SUBCASE 4'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR2) THEN
                SCASE = 'R1 ; SUBCASE 5'
                CALL CALCMDRP (RH, DRMR2, DRNACL, CALCR1A, CALCR3A) !, CALCR4A)
                SCASE = 'R1 ; SUBCASE 5'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR3) THEN
                SCASE = 'R1 ; SUBCASE 6'
                CALL CALCMDRP (RH, DRMR3, DRNANO3, CALCR1A, CALCR3A)
                SCASE = 'R1 ; SUBCASE 6'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR4) THEN
                SCASE = 'R1 ; SUBCASE 7'
                CALL CALCMDRP (RH, DRMR4, DRNACL, CALCR1A, CALCR3A) !, CALCR4A)
                SCASE = 'R1 ; SUBCASE 7'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR5) THEN
                SCASE = 'R1 ; SUBCASE 8'
                CALL CALCMDRP (RH, DRMR5, DRNH4CL, CALCR1A, CALCR3A) !, CALCR5)
                SCASE = 'R1 ; SUBCASE 8'
            ENDIF

        ELSE IF (     EXAC .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR6) THEN
                SCASE = 'R1 ; SUBCASE 9'
                CALL CALCMDRP (RH, DRMR6, DRNANO3, CALCR1A, CALCR3A)
                SCASE = 'R1 ; SUBCASE 9'
            ENDIF
        ENDIF
    
    ELSE IF ( .NOT. EXAC) THEN                   ! *** NH4CL  = 0
        IF      (     EXAN .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR7) THEN
                SCASE = 'R1 ; SUBCASE 10'
                CALL CALCMDRP (RH, DRMR7, DRNH4NO3, CALCR1A, CALCR2A)
                SCASE = 'R1 ; SUBCASE 10'
            ENDIF

        ELSE IF (     EXAN .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR8) THEN
                SCASE = 'R1 ; SUBCASE 11'
                CALL CALCMDRP (RH, DRMR8, DRNH4NO3, CALCR1A, CALCR2A)
                SCASE = 'R1 ; SUBCASE 11'
            ENDIF

        ELSE IF (     EXAN .AND. .NOT. EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR9) THEN
                SCASE = 'R1 ; SUBCASE 12'
                CALL CALCMDRP (RH, DRMR9, DRNH4NO3, CALCR1A, CALCR2A)
                SCASE = 'R1 ; SUBCASE 12'
            ENDIF

        ELSE IF (     EXAN .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR10) THEN
                SCASE = 'R1 ; SUBCASE 13'
                CALL CALCMDRP (RH, DRMR10, DRNH4NO3, CALCR1A, CALCR2A)
                SCASE = 'R1 ; SUBCASE 13'
            ENDIF
        ENDIF
    
    ELSE IF ( .NOT. EXSN) THEN                  ! *** NANO3  = 0
        IF      (     EXAN .AND.      EXAC .AND.      EXSC) THEN
            IF (RH >= DRMR11) THEN
                SCASE = 'R1 ; SUBCASE 14'
                CALL CALCMDRP (RH, DRMR11, DRNH4NO3, CALCR1A, CALCR2A)
                SCASE = 'R1 ; SUBCASE 14'
            ENDIF

        ELSE IF (     EXAN .AND.      EXAC .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR12) THEN
                SCASE = 'R1 ; SUBCASE 15'
                CALL CALCMDRP (RH, DRMR12, DRNH4NO3, CALCR1A, CALCR2A)
                SCASE = 'R1 ; SUBCASE 15'
            ENDIF
        ENDIF
    
    ELSE IF ( .NOT. EXSC) THEN                  ! *** NACL   = 0
        IF      (     EXAN .AND.      EXAC .AND.      EXSN) THEN
            IF (RH >= DRMR13) THEN
                SCASE = 'R1 ; SUBCASE 16'
                CALL CALCMDRP (RH, DRMR13, DRNH4NO3, CALCR1A, CALCR2A)
                SCASE = 'R1 ; SUBCASE 16'
            ENDIF
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCR1 ******************************************

    END SUBROUTINE CALCR1


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCR1A
! *** CASE R1 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCR1A
    INCLUDE 'isrpia.inc'

! *** CALCULATE SOLIDS **************************************************

    CNA2SO4 = WAER(2)
    FRNA    = MAX (WAER(1)-2*CNA2SO4, ZERO)

    CNH42S4 = ZERO

    CNANO3  = MIN (FRNA, WAER(4))
    FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
    FRNA    = MAX (FRNA-CNANO3, ZERO)

    CNACL   = MIN (FRNA, WAER(5))
    FRCL    = MAX (WAER(5)-CNACL, ZERO)
    FRNA    = MAX (FRNA-CNACL, ZERO)

    CNH4NO3 = MIN (FRNO3, WAER(3))
    FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
    FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)

    CNH4CL  = MIN (FRCL, FRNH3)
    FRCL    = MAX (FRCL-CNH4CL, ZERO)
    FRNH3   = MAX (FRNH3-CNH4CL, ZERO)

! *** OTHER PHASES ******************************************************

    WATER   = ZERO

    GNH3    = ZERO
    GHNO3   = ZERO
    GHCL    = ZERO

    RETURN

! *** END OF SUBROUTINE CALCR1A *****************************************

    END SUBROUTINE CALCR1A
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV7
! *** CASE V7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV7
    INCLUDE 'isrpia.inc'

    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCV1A

    CHI9   = CCASO4

    PSI1   = CNA2SO4      ! SALTS DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! ION CONCENTRATIONS
    
        NAI    = WAER(1)
        SO4I   = MAX (WAER(2) - WAER(6), ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        KI     = WAER(7)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                    ! H+ in excess
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

    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCV7')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNA2SO4 = ZERO
    CMGSO4  = ZERO
    CK2SO4  = ZERO
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCV7 ******************************************

    END SUBROUTINE CALCV7

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV6
! *** CASE V6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MGSO4, NA2SO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV6
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .TRUE. 
    PSI70   =-GREAT                                 ! GREAT = 1.D10
    ROOT7   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCV1A

    CHI9   = CCASO4
    CHI7   = CK2SO4       ! SALTS

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7))
            CC = WAER(7)*(WAER(2)-WAER(6)) + 0.25D0*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*WAER(2) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MIN (ROOT7,WAER(7)/2.0,MAX(WAER(2)-WAER(6),ZERO),CHI7)
            ROOT7 = MAX (ROOT7, ZERO)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT7, ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        NAI    = WAER(1)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCV6')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNA2SO4 = ZERO
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCV6 ******************************************

    END SUBROUTINE CALCV6
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV5
! *** CASE V5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MGSO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV5
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .TRUE. 
    PSCONV1 = .TRUE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCV1A

    CHI9   = CCASO4
    CHI7   = CK2SO4       ! SALTS
    CHI1   = CNA2SO4

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
            CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX(WAER(2)-WAER(6) - ROOT1, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
            CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MAX (ROOT1, ZERO)
            ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
            MAX ((WAER(2)-WAER(6)) - ROOT7, ZERO), CHI1)
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
        SO4I   = MAX ((WAER(2)-WAER(6)) - ROOT7 - ROOT1, ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCV5')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNA2SO4 = CHI1 - PSI1
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCV5******************************************

    END SUBROUTINE CALCV5

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV4
! *** CASE V4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV4
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .TRUE. 
    PSCONV1 = .TRUE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCV1A

    CHI9   = CCASO4
    CHI7   = CK2SO4       ! SALTS
    CHI1   = CNA2SO4
    CHI8   = CMGSO4

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
            CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX((WAER(2)-WAER(6)) - ROOT1, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
            CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MAX (ROOT1, ZERO)
            ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
            MAX ((WAER(2)-WAER(6)) - ROOT7, ZERO), CHI1)
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
        SO4I   = MAX ((WAER(2)-WAER(6)) - ROOT7 - ROOT1, ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCV4')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNA2SO4 = CHI1 - PSI1
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCV4******************************************

    END SUBROUTINE CALCV4

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV3
! *** CASE V3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, (NH4)2SO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV3
    INCLUDE 'isrpia.inc'
    LOGICAL :: EXNO, EXCL
    EXTERNAL CALCV1A, CALCV4

! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***

    EXNO = WAER(4) > TINY
    EXCL = WAER(5) > TINY

    IF (EXNO .OR. EXCL) THEN             ! *** NITRATE OR CHLORIDE EXISTS
        SCASE = 'V3 ; SUBCASE 1'
        CALL CALCV3A
        SCASE = 'V3 ; SUBCASE 1'
    
    ELSE                                 ! *** NO CHLORIDE AND NITRATE
        IF (RH < DRMO3) THEN
            SCASE = 'V3 ; SUBCASE 2'
            CALL CALCV1A             ! SOLID
            SCASE = 'V3 ; SUBCASE 2'
        ELSE
            SCASE = 'V3 ; SUBCASE 3' ! MDRH (CaSO4, (NH4)2SO4, MGSO4, NA2SO4, K2SO4)
            CALL CALCMDRPII (RH, DRMO3, DRNH42S4, CALCV1A, CALCV4)
            SCASE = 'V3 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCV3 ******************************************

    END SUBROUTINE CALCV3

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV3A
! *** CASE V3A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, (NH4)2SO4
!     4. Completely dissolved: NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV3A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1, PSCONV6
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .TRUE. 
    PSCONV1 = .TRUE. 
    PSCONV6 = .TRUE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT
    PSI60   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO
    ROOT6   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCV1A

    CHI9   = CCASO4
    CHI7   = CK2SO4       ! SALTS
    CHI1   = CNA2SO4
    CHI8   = CMGSO4
    CHI6   = CNH42S4

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = WAER(2)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        A6  = XK7 *(WATER/GAMA(4))**3.0        !(NH4)2SO4  <==> NH4+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1 - ROOT6)
            CC = WAER(7)*((WAER(2) - WAER(6)) - ROOT1 - ROOT6) + &
            &         0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6))-ROOT1-ROOT6)-A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX (WAER(2)-WAER(6)-ROOT1-ROOT6, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7 - ROOT6)
            CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7 - ROOT6) + &
            &         0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6))-ROOT7-ROOT6)-A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MAX (ROOT1, ZERO)
            ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
            MAX (WAER(2)-WAER(6)-ROOT7-ROOT6, ZERO), CHI1)
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! AMMONIUM SULFATE
    
        IF (NH4I*NH4I*SO4I > A6) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(3) - ROOT7 - ROOT1)
            CC = WAER(3)*((WAER(2)-WAER(6)) - ROOT7 - ROOT1) + &
            &         0.25*WAER(3)*WAER(3)
            DD =-0.25*(WAER(3)*WAER(3)*((WAER(2)-WAER(6))-ROOT7-ROOT1)-A6)
            CALL POLY3(BB, CC, DD, ROOT6, ISLV)
            IF (ISLV /= 0) ROOT6 = TINY
            ROOT6 = MAX (ROOT6, ZERO)
            ROOT6 = MIN (ROOT6, WAER(3)/2.0, &
            MAX (WAER(2)-WAER(6)-ROOT7-ROOT1, ZERO), CHI6)
            PSI6  = CHI6-ROOT6
        ENDIF
        PSCONV6 = ABS(PSI6-PSI60) <= EPS*PSI60
        PSI60   = PSI6
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1 - ROOT6, ZERO)
        NH4I   = MAX (WAER(3) - 2.D0*ROOT6, ZERO)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV6) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCV3')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = CHI6 - PSI6
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNA2SO4 = CHI1 - PSI1
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCV3A******************************************

    END SUBROUTINE CALCV3A

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCV2
! *** CASE V2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, (NH4)2SO4, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV2
    INCLUDE 'isrpia.inc'
    LOGICAL :: EXNO, EXCL
    EXTERNAL CALCV1A, CALCV3A, CALCV4

! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***

    EXNO = WAER(4) > TINY
    EXCL = WAER(5) > TINY

    IF (EXNO) THEN                       ! *** NITRATE EXISTS
        SCASE = 'V2 ; SUBCASE 1'
        CALL CALCV2A
        SCASE = 'V2 ; SUBCASE 1'
    
    ELSEIF ( .NOT. EXNO .AND. EXCL) THEN   ! *** ONLY CHLORIDE EXISTS
        IF (RH < DRMO2) THEN
            SCASE = 'V2 ; SUBCASE 2'
            CALL CALCV1A             ! SOLID
            SCASE = 'V2 ; SUBCASE 2'
        ELSE
            SCASE = 'V2 ; SUBCASE 3' ! MDRH CaSO4, NH4CL, (NH4)2SO4, MGSO4, NA2SO4, K2SO4
            CALL CALCMDRPII (RH, DRMO2, DRNH4CL, CALCV1A, CALCV3A)
            SCASE = 'V2 ; SUBCASE 3'
        ENDIF
    
    ELSE                                 ! *** NO CHLORIDE AND NITRATE
        IF (RH < DRMO3) THEN
            SCASE = 'V2 ; SUBCASE 2'
            CALL CALCV1A             ! SOLID
            SCASE = 'V2 ; SUBCASE 2'
        ELSE
            SCASE = 'V2 ; SUBCASE 4' ! MDRH CaSO4, (NH4)2SO4, MGSO4, NA2SO4, K2SO4
            CALL CALCMDRPII (RH, DRMO3, DRNH42S4, CALCV1A, CALCV4)
            SCASE = 'V2 ; SUBCASE 4'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCV2 ******************************************

    END SUBROUTINE CALCV2

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV2A
! *** CASE V2A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, (NH4)2SO4, NH4CL
!     4. Completely dissolved: NH4NO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV2A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1, PSCONV6, PSCONV4
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .TRUE. 
    PSCONV1 = .TRUE. 
    PSCONV6 = .TRUE. 
    PSCONV4 = .TRUE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT
    PSI60   =-GREAT
    PSI40   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO
    ROOT6   = ZERO
    ROOT4   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCV1A

    CHI9   = CCASO4
    CHI7   = CK2SO4       ! SALTS
    CHI1   = CNA2SO4
    CHI8   = CMGSO4
    CHI6   = CNH42S4
    CHI4   = CNH4CL

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI6   = CNH42S4
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        A6  = XK7 *(WATER/GAMA(4))**3.0        ! (NH4)2SO4 <==> NH4+
        A14 = XK14*(WATER/GAMA(6))**2.         ! NH4Cl     <==> NH4+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! AMMONIUM CHLORIDE
    
        IF (NH4I*CLI > A14) THEN
            BB    =-(WAER(3) + WAER(5) - 2.D0*ROOT6)
            CC    = WAER(5)*(WAER(3) - 2.D0*ROOT6) - A14
            DD    = BB*BB - 4.D0*CC
            IF (DD < ZERO) THEN
                ROOT4 = ZERO
            ELSE
                DD    = SQRT(DD)
                ROOT4A= 0.5D0*(-BB+DD)
                ROOT4B= 0.5D0*(-BB-DD)
                IF (ZERO <= ROOT4A) THEN
                    ROOT4 = ROOT4A
                ELSE
                    ROOT4 = ROOT4B
                ENDIF
                ROOT4 = MAX(ROOT4, ZERO)
                ROOT4 = MIN(ROOT4, WAER(5), &
                MAX (WAER(3) - 2.D0*ROOT6, ZERO), CHI4)
                PSI4  = CHI4 - ROOT4
            ENDIF
        ENDIF
        PSCONV4 = ABS(PSI4-PSI40) <= EPS*PSI40
        PSI40   = PSI4
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2) - WAER(6)) + WAER(7) - ROOT1 - ROOT6)
            CC = WAER(7)*((WAER(2) - WAER(6)) - ROOT1 - ROOT6) &
            + 0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6))-ROOT1-ROOT6)-A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX (WAER(2)-WAER(6)-ROOT1-ROOT6, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-((WAER(2) - WAER(6)) + WAER(1) - ROOT7 - ROOT6)
            CC = WAER(1)*((WAER(2) - WAER(6)) - ROOT7 - ROOT6) + &
            &         0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6))-ROOT7-ROOT6)-A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MAX (ROOT1, ZERO)
            ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
            MAX (WAER(2)-WAER(6)-ROOT7-ROOT6, ZERO), CHI1)
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! AMMONIUM SULFATE
    
        IF (NH4I*NH4I*SO4I > A6) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(3) - ROOT7 - ROOT1 - ROOT4)
            CC = WAER(3)*((WAER(2)-WAER(6)) - ROOT7 - ROOT1) + 0.25* &
            (WAER(3)-ROOT4)**2.0 + ROOT4*(ROOT1+ROOT7-(WAER(2)-WAER(6)))
            DD =-0.25*((WAER(3)-ROOT4)**2.0 * &
            ((WAER(2)-WAER(6))-ROOT7-ROOT1) - A6)
            CALL POLY3(BB, CC, DD, ROOT6, ISLV)
            IF (ISLV /= 0) ROOT6 = TINY
            ROOT6 = MAX (ROOT6, ZERO)
            ROOT6 = MIN (ROOT6, WAER(3)/2.0, &
            MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1, ZERO), CHI6)
            PSI6  = CHI6-ROOT6
        ENDIF
        PSCONV6 = ABS(PSI6-PSI60) <= EPS*PSI60
        PSI60   = PSI6
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1 - ROOT6, ZERO)
        NH4I   = MAX (WAER(3) - 2.D0*ROOT6, ZERO)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV6 .AND. PSCONV4) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCV2')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = CHI6 - PSI6
    CNH4NO3 = ZERO
    CNH4CL  = CHI4 - PSI4
    CNA2SO4 = CHI1 - PSI1
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCV2A******************************************

    END SUBROUTINE CALCV2A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV1
! *** CASE V1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV1
    INCLUDE 'isrpia.inc'
    LOGICAL :: EXNO, EXCL
    EXTERNAL CALCV1A, CALCV2A, CALCV3A, CALCV4

! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***

    EXNO = WAER(4) > TINY
    EXCL = WAER(5) > TINY

    IF (EXNO .AND. EXCL) THEN           ! *** NITRATE & CHLORIDE EXIST
        IF (RH < DRMO1) THEN
            SCASE = 'V1 ; SUBCASE 1'
            CALL CALCV1A             ! SOLID
            SCASE = 'V1 ; SUBCASE 1'
        ELSE
            SCASE = 'V1 ; SUBCASE 2' ! MDRH (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4
            CALL CALCMDRPII (RH, DRMO1, DRNH4NO3, CALCV1A, CALCV2A)
            SCASE = 'V1 ; SUBCASE 2'
        ENDIF
    
    ELSE IF (EXNO .AND. .NOT. EXCL) THEN ! *** ONLY NITRATE EXISTS
        IF (RH < DRMV1) THEN
            SCASE = 'V1 ; SUBCASE 1'
            CALL CALCV1A             ! SOLID
            SCASE = 'V1 ; SUBCASE 1'
        ELSE
            SCASE = 'V1 ; SUBCASE 3' ! MDRH (NH4)2SO4, NH4NO3, NA2SO4, K2SO4, MGSO4, CASO4
            CALL CALCMDRPII (RH, DRMV1, DRNH4NO3, CALCV1A, CALCV2A)
            SCASE = 'V1 ; SUBCASE 3'
        ENDIF
    
    ELSE IF ( .NOT. EXNO .AND. EXCL) THEN ! *** ONLY CHLORIDE EXISTS
        IF (RH < DRMO2) THEN
            SCASE = 'V1 ; SUBCASE 1'
            CALL CALCV1A             ! SOLID
            SCASE = 'V1 ; SUBCASE 1'
        ELSE
            SCASE = 'V1 ; SUBCASE 4' ! MDRH (NH4)2SO4, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4
            CALL CALCMDRPII (RH, DRMO2, DRNH4CL, CALCV1A, CALCV3A)
            SCASE = 'V1 ; SUBCASE 4'
        ENDIF
    
    ELSE                                ! *** NO CHLORIDE AND NITRATE
        IF (RH < DRMO3) THEN
            SCASE = 'V1 ; SUBCASE 1'
            CALL CALCV1A             ! SOLID
            SCASE = 'V1 ; SUBCASE 1'
        ELSE
            SCASE = 'V1 ; SUBCASE 5' ! MDRH (NH4)2SO4, NA2SO4, K2SO4, MGSO4, CASO4
            CALL CALCMDRPII (RH, DRMO3, DRNH42S4, CALCV1A, CALCV4)
            SCASE = 'V1 ; SUBCASE 5'
        ENDIF
    ENDIF

    RETURN

!      IF (RH.LT.DRMO1) THEN
!         SCASE = 'V1 ; SUBCASE 1'
!         CALL CALCV1A              ! SOLID PHASE ONLY POSSIBLE
!         SCASE = 'V1 ; SUBCASE 1'
!      ELSE
!         SCASE = 'V1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
!         CALL CALCMDRPII (RH, DRMO1, DRNH4NO3, CALCV1A, CALCV2A)
!         SCASE = 'V1 ; SUBCASE 2'
!         ENDIF

!      RETURN

! *** END OF SUBROUTINE CALCV1 ******************************************

    END SUBROUTINE CALCV1

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCV1A
! *** CASE V1A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCV1A
    INCLUDE 'isrpia.inc'

! *** CALCULATE SOLIDS **************************************************

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

    CNH4NO3 = MIN (FRNH3, WAER(4))
! C      FRNO3   = MAX (WAER(4) - CNH4NO3, ZERO)
    FRNH3   = MAX (FRNH3 - CNH4NO3, ZERO)

    CNH4CL  = MIN (FRNH3, WAER(5))
! C      FRCL    = MAX (WAER(5) - CNH4CL, ZERO)
    FRNH3   = MAX (FRNH3 - CNH4CL, ZERO)

! *** OTHER PHASES ******************************************************

    WATER   = ZERO

    GNH3    = ZERO
    GHNO3   = ZERO
    GHCL    = ZERO

    RETURN

! *** END OF SUBROUTINE CALCV1A *****************************************

    END SUBROUTINE CALCV1A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCU8
! *** CASE U8

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
!     2. THERE IS ONLY A LIQUID PHASE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU8
    INCLUDE 'isrpia.inc'

    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    CALL CALCU1A

    CHI9   = CCASO4        ! SALTS

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    CALAOU = .TRUE. 

! *** CALCULATE WATER **************************************************

    CALL CALCMR

! *** SETUP LIQUID CONCENTRATIONS **************************************

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+
    
        NAI    = WAER(1)
        SO4I   = MAX(WAER(2) - WAER(6), ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        KI     = WAER(7)
        MGI    = WAER(8)

    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG  = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
        IF (HI <= TINY) HI = SQRT(AKW)
    !      OHI   = AKW/HI
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCU8')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3     = NH4I/HI/A2
    GHNO3    = HI*NO3I/A3
    GHCL     = HI*CLI /A4

    GASAQ(1) = NH3AQ
    GASAQ(2) = CLAQ
    GASAQ(3) = NO3AQ

    CNH42S4  = ZERO
    CNH4NO3  = ZERO
    CNH4CL   = ZERO
    CNACL    = ZERO
    CNANO3   = ZERO
    CNA2SO4  = ZERO
    CMGSO4   = ZERO
    CK2SO4   = ZERO
    CCASO4   = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCU8 ******************************************

    END SUBROUTINE CALCU8

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU7
! *** CASE U7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MGSO4, NA2SO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU7
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .TRUE. 
    PSI70   =-GREAT                                 ! GREAT = 1.D10
    ROOT7   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCU1A

    CHI7   = CK2SO4       ! SALTS
    CHI9   = CCASO4

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4


    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

    MOLAL(1) = ZERO
    MOLAL(2) = NAI
    MOLAL(3) = NH4I
    MOLAL(4) = CLI
    MOLAL(5) = SO4I
    MOLAL(6) = HSO4I
    MOLAL(7) = NO3I
    MOLAL(8) = CAI
    MOLAL(9) = KI
    MOLAL(10)= MGI

    CALL CALCACT          ! CALCULATE ACTIVITY COEFFICIENTS

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7))
            CC = WAER(7)*(WAER(2)-WAER(6)) + 0.25D0*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*(WAER(2)-WAER(6)) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7,WAER(7)/2.0,MAX(WAER(2)-WAER(6),ZERO),CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        SO4I   = MAX (WAER(2) - WAER(6) - ROOT7, ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        NAI    = WAER(1)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    !      IF (HI.LE.TINY) HI = SQRT(AKW)
    !      OHI   = AKW/HI
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCU7')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = ZERO
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCU7 ******************************************

    END SUBROUTINE CALCU7
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU6
! *** CASE U6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MGSO4

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU6
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .TRUE. 
    PSCONV1 = .TRUE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCU1A

    CHI1   = CNA2SO4            ! SALTS
    CHI7   = CK2SO4
    CHI9   = CCASO4

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

    MOLAL(1) = ZERO
    MOLAL(2) = NAI
    MOLAL(3) = NH4I
    MOLAL(4) = CLI
    MOLAL(5) = SO4I
    MOLAL(6) = HSO4I
    MOLAL(7) = NO3I
    MOLAL(8) = CAI
    MOLAL(9) = KI
    MOLAL(10)= MGI

    CALL CALCACT          ! CALCULATE ACTIVITY COEFFICIENTS

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
            CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX((WAER(2)-WAER(6)) - ROOT1,ZERO), CHI7)
            PSI7  = CHI7-ROOT7

        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
            CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MAX (ROOT1, ZERO)
            ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
            MAX((WAER(2)-WAER(6)) - ROOT7, ZERO) ,CHI1)
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
        SO4I   = MAX (WAER(2) - WAER(6) - ROOT7 - ROOT1, ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    !      IF (HI.LE.TINY) HI = SQRT(AKW)
    !      OHI   = AKW/HI
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCU6')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = CHI1 - PSI1
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCU6******************************************

    END SUBROUTINE CALCU6
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU5
! *** CASE U5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU5
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .TRUE. 
    PSCONV1 = .TRUE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCU1A

    CHI1   = CNA2SO4            ! SALTS
    CHI7   = CK2SO4
    CHI8   = CMGSO4
    CHI9   = CCASO4

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

    MOLAL(1) = ZERO
    MOLAL(2) = NAI
    MOLAL(3) = NH4I
    MOLAL(4) = CLI
    MOLAL(5) = SO4I
    MOLAL(6) = HSO4I
    MOLAL(7) = NO3I
    MOLAL(8) = CAI
    MOLAL(9) = KI
    MOLAL(10)= MGI

    CALL CALCACT          ! CALCULATE ACTIVITY COEFFICIENTS

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
            CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX(WAER(2)-WAER(6)-ROOT1, ZERO),CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
            CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MAX (ROOT1, ZERO)
            ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
            MAX(WAER(2)-WAER(6)-ROOT7, ZERO),CHI1)
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1, ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    !      IF (HI.LE.TINY) HI = SQRT(AKW)
    !      OHI   = AKW/HI
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCU5')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = CHI1 - PSI1
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCU5******************************************

    END SUBROUTINE CALCU5

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU4
! *** CASE U4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU4
    INCLUDE 'isrpia.inc'
    LOGICAL ::  EXAN, EXAC, EXSN, EXSC
    EXTERNAL CALCU1A, CALCU5

! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************

    SCASE = 'U4 ; SUBCASE 2'
    CALL CALCU1A              ! SOLID
    SCASE = 'U4 ; SUBCASE 2'

    EXAN = CNH4NO3 > TINY    ! NH4NO3
    EXAC = CNH4CL > TINY    ! NH4CL
    EXSN = CNANO3 > TINY    ! NANO3
    EXSC = CNACL  > TINY    ! NACL

! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********

    IF (EXAN .OR. EXSN .OR. EXSC) THEN   ! *** NH4NO3,NANO3 EXIST
        IF (RH >= DRMM1) THEN
            SCASE = 'U4 ; SUBCASE 1'
            CALL CALCU4A
            SCASE = 'U4 ; SUBCASE 1'
        ENDIF
    
    ELSE IF (EXAC) THEN                  ! *** NH4CL EXISTS ONLY
        IF (RH >= DRMR5) THEN
            SCASE = 'U4 ; SUBCASE 3'
            CALL CALCMDRPII (RH, DRMR5, DRNH4CL, CALCU1A, CALCU5)
            SCASE = 'U4 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCU4 ******************************************

    END SUBROUTINE CALCU4

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU4A
! *** CASE U4A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, NH4CL
!     4. Completely dissolved: NH4NO3, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU4A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1, PSCONV4
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .FALSE. 
    PSCONV1 = .FALSE. 
    PSCONV4 = .FALSE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT
    PSI40   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO
    ROOT4   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCU1A

    CHI1   = CNA2SO4            ! SALTS
    CHI4   = CNH4CL
    CHI7   = CK2SO4
    CHI8   = CMGSO4
    CHI9   = CCASO4

    PSI1   = CNA2SO4
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

    MOLAL(1) = ZERO
    MOLAL(2) = NAI
    MOLAL(3) = NH4I
    MOLAL(4) = CLI
    MOLAL(5) = SO4I
    MOLAL(6) = HSO4I
    MOLAL(7) = NO3I
    MOLAL(8) = CAI
    MOLAL(9) = KI
    MOLAL(10)= MGI

    CALL CALCACT          ! CALCULATE ACTIVITY COEFFICIENTS

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        A14 = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
            CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX(WAER(2)-WAER(6)-ROOT1, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(1) - ROOT7)
            CC = WAER(1)*((WAER(2)-WAER(6)) - ROOT7) + 0.25*WAER(1)*WAER(1)
            DD =-0.25*(WAER(1)*WAER(1)*((WAER(2)-WAER(6)) - ROOT7) - A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MAX (ROOT1, ZERO)
            ROOT1 = MIN (ROOT1, WAER(1)/2.0, &
            MAX (WAER(2)-WAER(6)-ROOT7, ZERO), CHI1)
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! AMMONIUM CHLORIDE
    
        IF (NH4I*CLI > A14) THEN
            BB   =-(NH4I + CLI)
            CC   =-A14 + NH4I*CLI
            DD   = BB*BB - 4.D0*CC
            ROOT4 = 0.5D0*(-BB-SQRT(DD))
            IF (ROOT4 > TINY) THEN
                ROOT4    = MIN(MAX (ROOT4, ZERO), CHI4)
                PSI4    = CHI4 - ROOT4
            ENDIF
        ENDIF
        PSCONV4 = ABS(PSI4-PSI40) <= EPS*PSI40
        PSI40   = PSI4
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.D0*ROOT1, ZERO)
        SO4I   = MAX (WAER(2) - WAER(6) - ROOT7 - ROOT1, ZERO)
        NH4I   = MAX (WAER(3) - ROOT4, ZERO)
        NO3I   = WAER(4)
        CLI    = MAX (WAER(5) - ROOT4, ZERO)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    !      IF (HI.LE.TINY) HI = SQRT(AKW)
    !      OHI   = AKW/HI
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV4) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCU4')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI4 - PSI4
    CNACL   = ZERO
    CNANO3  = ZERO
    CNA2SO4 = CHI1 - PSI1
    CMGSO4  = ZERO
    CK2SO4  = CHI7 - PSI7
    CCASO4  = MIN (WAER(6), WAER(2))

    RETURN

! *** END OF SUBROUTINE CALCU4A ****************************************

    END SUBROUTINE CALCU4A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU3
! *** CASE U3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, NH4CL, NANO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU3
    INCLUDE 'isrpia.inc'
    LOGICAL ::  EXAN, EXAC, EXSN, EXSC
    EXTERNAL CALCU1A, CALCU4A, CALCU5

! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************

    SCASE = 'U3 ; SUBCASE 2'
    CALL CALCU1A              ! SOLID
    SCASE = 'U3 ; SUBCASE 2'

    EXAN = CNH4NO3 > TINY    ! NH4NO3
    EXAC = CNH4CL > TINY    ! NH4CL
    EXSN = CNANO3 > TINY    ! NANO3
    EXSC = CNACL  > TINY    ! NACL

! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********

    IF (EXAN .OR. EXSN) THEN                   ! *** NH4NO3,NANO3 EXIST
        IF (RH >= DRMM1) THEN
            SCASE = 'U3 ; SUBCASE 1'
            CALL CALCU3A
            SCASE = 'U3 ; SUBCASE 1'
        ENDIF
    
    ELSE IF ( .NOT. EXAN .AND. .NOT. EXSN) THEN   ! *** NH4NO3,NANO3 = 0
        IF      (     EXAC .AND.      EXSC) THEN
            IF (RH >= DRMR4) THEN
                SCASE = 'U3 ; SUBCASE 3'
                CALL CALCMDRPII (RH, DRMR4, DRNACL, CALCU1A, CALCU4A)
                SCASE = 'U3 ; SUBCASE 3'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSC) THEN
            IF (RH >= DRMR2) THEN
                SCASE = 'U3 ; SUBCASE 4'
                CALL CALCMDRPII (RH, DRMR2, DRNACL, CALCU1A, CALCU4A)
                SCASE = 'U3 ; SUBCASE 4'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR5) THEN
                SCASE = 'U3 ; SUBCASE 5'
                CALL CALCMDRPII (RH, DRMR5, DRNACL, CALCU1A, CALCU5)
                SCASE = 'U3 ; SUBCASE 5'
            ENDIF
        ENDIF
    
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCU3 ******************************************

    END SUBROUTINE CALCU3

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU3A
! *** CASE U3A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, NH4CL, NACL
!     4. Completely dissolved: NH4NO3, NANO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU3A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1, PSCONV4, PSCONV3
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .FALSE. 
    PSCONV1 = .FALSE. 
    PSCONV4 = .FALSE. 
    PSCONV3 = .FALSE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT
    PSI40   =-GREAT
    PSI30   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO
    ROOT4   = ZERO
    ROOT3   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCU1A

    CHI1   = CNA2SO4            ! SALTS
    CHI3   = CNACL
    CHI4   = CNH4CL
    CHI7   = CK2SO4
    CHI8   = CMGSO4
    CHI9   = CCASO4

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

    MOLAL(1) = ZERO
    MOLAL(2) = NAI
    MOLAL(3) = NH4I
    MOLAL(4) = CLI
    MOLAL(5) = SO4I
    MOLAL(6) = HSO4I
    MOLAL(7) = NO3I
    MOLAL(8) = CAI
    MOLAL(9) = KI
    MOLAL(10)= MGI

    CALL CALCACT          ! CALCULATE ACTIVITY COEFFICIENTS

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        A14 = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
        A8  = XK8 *(WATER/GAMA(1))**2.0        ! NaCl      <==> Na+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
            CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX(WAER(2)-WAER(6)-ROOT1, ZERO),CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-(((WAER(2)-WAER(6))-ROOT7)*(WAER(1) - ROOT3))
            CC = ((WAER(2) - WAER(6)) - ROOT7)*(WAER(1) - ROOT3) + &
            &           0.25D0*(WAER(1) - ROOT3)**2.
            DD =-0.25D0*(((WAER(2) - WAER(6)) - ROOT7)* &
            (WAER(1) - ROOT3)**2.D0 - A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MIN (MAX(ROOT1, ZERO), MAX(WAER(1) - ROOT3, ZERO), &
            CHI1, MAX(WAER(2)-WAER(6), ZERO))
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! AMMONIUM CHLORIDE
    
        IF (NH4I*CLI > A14) THEN
            BB    =-(WAER(3) + WAER(5) - ROOT4)
            CC    =-A14 + NH4I*(WAER(5) - ROOT4)
            DD    = MAX(BB*BB - 4.D0*CC, ZERO)
            ROOT4A= 0.5D0*(-BB+SQRT(DD))
            ROOT4B= 0.5D0*(-BB-SQRT(DD))
            IF (ZERO <= ROOT4A) THEN
                ROOT4 = ROOT4A
            ELSE
                ROOT4 = ROOT4B
            ENDIF
            ROOT4 = MIN(MAX(ZERO, ROOT4), MAX(WAER(5)-ROOT3,ZERO), &
            CHI4, WAER(3))
            PSI4  = CHI4 - ROOT4
        ENDIF
        PSCONV4 = ABS(PSI4-PSI40) <= EPS*PSI40
        PSI40   = PSI4
    
    ! SODIUM CHLORIDE  ; To obtain new value for ROOT3
    
        IF (NAI*CLI > A8) THEN
            BB    =-((CHI1-2.D0*ROOT1) + (WAER(5) - ROOT4))
            CC    = (CHI1-2.D0*ROOT1)*(WAER(5) - ROOT4) - A8
            DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT3A= 0.5D0*(-BB-SQRT(DD))
            ROOT3B= 0.5D0*(-BB+SQRT(DD))
            IF (ZERO <= ROOT3A) THEN
                ROOT3 = ROOT3A
            ELSE
                ROOT3 = ROOT3B
            ENDIF
            ROOT3   = MIN(MAX(ROOT3, ZERO), CHI3)
            PSI3    = CHI3-ROOT3
        ENDIF
        PSCONV3 = ABS(PSI3-PSI30) <= EPS*PSI30
        PSI30   = PSI3
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.D0*ROOT1 - ROOT3, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT7 - ROOT1, ZERO)
        NH4I   = MAX (WAER(3) - ROOT4, ZERO)
        NO3I   = WAER(4)
        CLI    = MAX (WAER(5) - ROOT4 - ROOT3, ZERO)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    !      IF (HI.LE.TINY) HI = SQRT(AKW)
    !      OHI   = AKW/HI
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV4 .AND. PSCONV3) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCU3A')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 IF (CLI <= TINY .AND. WAER(5) > TINY) THEN !No disslv Cl-;solid only
        DO 30 I=1,NIONS
            MOLAL(I) = ZERO
        30 END DO
        DO 40 I=1,NGASAQ
            GASAQ(I) = ZERO
        40 END DO
        CALL CALCU1A
    ELSE
        A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
    
        GNH3    = NH4I/HI/A2
        GHNO3   = HI*NO3I/A3
        GHCL    = HI*CLI /A4
    
        GASAQ(1)= NH3AQ
        GASAQ(2)= CLAQ
        GASAQ(3)= NO3AQ
    
        CNH42S4 = ZERO
        CNH4NO3 = ZERO
        CNH4CL  = CHI4 - PSI4
        CNACL   = CHI3 - PSI3
        CNANO3  = ZERO
        CNA2SO4 = CHI1 - PSI1
        CMGSO4  = ZERO
        CK2SO4  = CHI7 - PSI7
        CCASO4  = MIN (WAER(6), WAER(2))
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCU3A*****************************************

    END SUBROUTINE CALCU3A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU2
! *** CASE U2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, NH4CL, NANO3, NACL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU2
    INCLUDE 'isrpia.inc'
    LOGICAL ::  EXAN, EXAC, EXSN, EXSC
    EXTERNAL CALCU1A, CALCU3A, CALCU4A, CALCU5

! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************

    SCASE = 'U2 ; SUBCASE 2'
    CALL CALCU1A              ! SOLID
    SCASE = 'U2 ; SUBCASE 2'

    EXAN = CNH4NO3 > TINY    ! NH4NO3
    EXAC = CNH4CL > TINY    ! NH4CL
    EXSN = CNANO3 > TINY    ! NANO3
    EXSC = CNACL  > TINY    ! NACL

! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********

    IF (EXAN) THEN                             ! *** NH4NO3 EXISTS
        IF (RH >= DRMM1) THEN
            SCASE = 'U2 ; SUBCASE 1'
            CALL CALCU2A
            SCASE = 'U2 ; SUBCASE 1'
        ENDIF
    
    ELSE IF ( .NOT. EXAN) THEN                   ! *** NH4NO3 = 0
        IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMM2) THEN
                SCASE = 'U2 ; SUBCASE 3'
                CALL CALCMDRPII (RH, DRMM2, DRNANO3, CALCU1A, CALCU3A)
                SCASE = 'U2 ; SUBCASE 3'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR1) THEN
                SCASE = 'U2 ; SUBCASE 4'
                CALL CALCMDRPII (RH, DRMR1, DRNANO3, CALCU1A, CALCU3A)
                SCASE = 'U2 ; SUBCASE 4'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR2) THEN
                SCASE = 'U2 ; SUBCASE 5'
                CALL CALCMDRPII (RH, DRMR2, DRNACL, CALCU1A, CALCU4A)
                SCASE = 'U2 ; SUBCASE 5'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR3) THEN
                SCASE = 'U2 ; SUBCASE 6'
                CALL CALCMDRPII (RH, DRMR3, DRNANO3, CALCU1A, CALCU3A)
                SCASE = 'U2 ; SUBCASE 6'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR4) THEN
                SCASE = 'U2 ; SUBCASE 7'
                CALL CALCMDRPII (RH, DRMR4, DRNACL, CALCU1A, CALCU4A)
                SCASE = 'U2 ; SUBCASE 7'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR5) THEN
                SCASE = 'U2 ; SUBCASE 8'
                CALL CALCMDRPII (RH, DRMR5, DRNH4CL, CALCU1A, CALCU5)
                SCASE = 'U2 ; SUBCASE 8'
            ENDIF

        ELSE IF (     EXAC .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR6) THEN
                SCASE = 'U2 ; SUBCASE 9'
                CALL CALCMDRPII (RH, DRMR6, DRNANO3, CALCU1A, CALCU3A)
                SCASE = 'U2 ; SUBCASE 9'
            ENDIF
        ENDIF
    
    ENDIF

    RETURN

!      IF (W(4).GT.TINY) THEN        ! NO3 EXISTS, WATER POSSIBLE
!         SCASE = 'U2 ; SUBCASE 1'
!         CALL CALCU2A
!         SCASE = 'U2 ; SUBCASE 1'
!      ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
!         SCASE = 'U2 ; SUBCASE 1'
!         CALL CALCU1A
!         SCASE = 'U2 ; SUBCASE 1'
!      ENDIF
!C
!      IF (WATER.LE.TINY .AND. RH.LT.DRMM2) THEN      ! DRY AEROSOL
!         SCASE = 'U2 ; SUBCASE 2'
!         CALL CALCU2A
!         SCASE = 'U2 ; SUBCASE 1'
!C
!      ELSEIF (WATER.LE.TINY .AND. RH.GE.DRMM2) THEN  ! MDRH OF M2
!         SCASE = 'U2 ; SUBCASE 3'
!         CALL CALCMDRPII (RH, DRMM2, DRNANO3, CALCU1A, CALCU3A)
!         SCASE = 'U2 ; SUBCASE 3'
!      ENDIF
!C
!      RETURN

! *** END OF SUBROUTINE CALCU2 ******************************************

    END SUBROUTINE CALCU2

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU2A
! *** CASE U2A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, NH4CL, NACL, NANO3
!     4. Completely dissolved: NH4NO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU2A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV7, PSCONV1, PSCONV4, PSCONV3, PSCONV5
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV7 = .FALSE. 
    PSCONV1 = .FALSE. 
    PSCONV4 = .FALSE. 
    PSCONV3 = .FALSE. 
    PSCONV5 = .FALSE. 

    PSI70   =-GREAT                                 ! GREAT = 1.D10
    PSI1O   =-GREAT
    PSI40   =-GREAT
    PSI30   =-GREAT
    PSI50   =-GREAT

    ROOT7   = ZERO
    ROOT1   = ZERO
    ROOT4   = ZERO
    ROOT3   = ZERO
    ROOT5   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCU1A

    CHI1   = CNA2SO4            ! SALTS
    CHI2   = CNANO3
    CHI3   = CNACL
    CHI4   = CNH4CL
    CHI7   = CK2SO4
    CHI8   = CMGSO4
    CHI9   = CCASO4

    PSI1   = CNA2SO4      ! AMOUNT DISSOLVED
    PSI2   = CNANO3
    PSI3   = CNACL
    PSI4   = CNH4CL
    PSI5   = CNH4NO3
    PSI7   = CK2SO4
    PSI8   = CMGSO4
    PSI9   = CCASO4

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

    MOLAL(1) = ZERO
    MOLAL(2) = NAI
    MOLAL(3) = NH4I
    MOLAL(4) = CLI
    MOLAL(5) = SO4I
    MOLAL(6) = HSO4I
    MOLAL(7) = NO3I
    MOLAL(8) = CAI
    MOLAL(9) = KI
    MOLAL(10)= MGI

    CALL CALCACT          ! CALCULATE ACTIVITY COEFFICIENTS

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A7  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A1  = XK5 *(WATER/GAMA(2))**3.0        ! NA2S04    <==> Na+
        A14 = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        A8  = XK8 *(WATER/GAMA(1))**2.0        ! NaCl      <==> Na+
        A9  = XK9 *(WATER/GAMA(3))**2.0        ! NaNO3     <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A7) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT1)
            CC = WAER(7)*((WAER(2)-WAER(6)) - ROOT1) + 0.25*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*((WAER(2)-WAER(6)) - ROOT1) - A7)
            CALL POLY3(BB, CC, DD, ROOT7, ISLV)
            IF (ISLV /= 0) ROOT7 = TINY
            ROOT7 = MAX (ROOT7, ZERO)
            ROOT7 = MIN (ROOT7, WAER(7)/2.0, &
            MAX(WAER(2)-WAER(6)-ROOT1, ZERO),CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI70) <= EPS*PSI70
        PSI70   = PSI7
    
    ! SODIUM SULFATE
    
        IF (NAI*NAI*SO4I > A1) THEN
            BB =-(((WAER(2)-WAER(6))-ROOT7)*(WAER(1) - ROOT3 - ROOT5))
            CC = ((WAER(2)-WAER(6)) - ROOT7)*(WAER(1) - ROOT3 - ROOT5) + &
            &          0.25D0*(WAER(1) - ROOT3 - ROOT5)**2.0
            DD =-0.25D0*(((WAER(2) - WAER(6)) - ROOT7)* &
            (WAER(1) - ROOT3 - ROOT5)**2.D0 - A1)
            CALL POLY3(BB, CC, DD, ROOT1, ISLV)
            IF (ISLV /= 0) ROOT1 = TINY
            ROOT1 = MIN (MAX(ROOT1,ZERO), MAX(WAER(1)-ROOT3-ROOT5,ZERO), &
            CHI1, MAX(WAER(2)-WAER(6),ZERO))
            PSI1  = CHI1-ROOT1
        ENDIF
        PSCONV1 = ABS(PSI1-PSI1O) <= EPS*PSI1O
        PSI1O   = PSI1
    
    ! AMMONIUM CHLORIDE
    
        IF (NH4I*CLI > A14) THEN
            BB    =-(WAER(3) + WAER(5) - ROOT4)
            CC    =-A14 + NH4I*(WAER(5) - ROOT4)
            DD    = MAX(BB*BB - 4.D0*CC, ZERO)
            ROOT4A= 0.5D0*(-BB+SQRT(DD))
            ROOT4B= 0.5D0*(-BB-SQRT(DD))
            IF (ZERO <= ROOT4A) THEN
                ROOT4 = ROOT4A
            ELSE
                ROOT4 = ROOT4B
            ENDIF
            ROOT4 = MIN(MAX(ZERO, ROOT4), MAX(WAER(5)-ROOT3,ZERO), &
            CHI4, WAER(3))
            PSI4  = CHI4 - ROOT4
        ENDIF
        PSCONV4 = ABS(PSI4-PSI40) <= EPS*PSI40
        PSI40   = PSI4
    
    ! SODIUM CHLORIDE  ; To obtain new value for ROOT3
    
        IF (NAI*CLI > A8) THEN
            BB    =-((CHI1-2.D0*ROOT1-ROOT5) + (WAER(5) - ROOT4))
            CC    = (CHI1-2.D0*ROOT1-ROOT5)*(WAER(5) - ROOT4) - A8
            DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT3A= 0.5D0*(-BB-SQRT(DD))
            ROOT3B= 0.5D0*(-BB+SQRT(DD))
            IF (ZERO <= ROOT3A) THEN
                ROOT3 = ROOT3A
            ELSE
                ROOT3 = ROOT3B
            ENDIF
            ROOT3   = MIN(MAX(ROOT3, ZERO), CHI3)
            PSI3    = CHI3-ROOT3
        ENDIF
        PSCONV3 = ABS(PSI3-PSI30) <= EPS*PSI30
        PSI30   = PSI3
    
    ! SODIUM NITRATE
    
        IF (NAI*NO3I > A9) THEN
            BB    =-(WAER(4) + WAER(1) - 2.D0*ROOT1 - ROOT3)
            CC    = WAER(4)*(WAER(1) - 2.D0*ROOT1 - ROOT3) - A9
            DD    = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT5A= 0.5D0*(-BB-DD)
            ROOT5B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT5A) THEN
                ROOT5 = ROOT5A
            ELSE
                ROOT5 = ROOT5B
            ENDIF
            ROOT5 = MIN(MAX(ROOT5, ZERO), CHI2)
            PSI2  = CHI2-ROOT5
        ENDIF
    
        PSCONV5 = ABS(PSI2-PSI20) <= EPS*PSI20
        PSI20   = PSI2
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.0D0*ROOT7, ZERO)
        NAI    = MAX (WAER(1) - 2.0D0*ROOT1 - ROOT3 - ROOT5, ZERO)
        SO4I   = MAX (WAER(2) - WAER(6) - ROOT7 - ROOT1, ZERO)
        NH4I   = MAX (WAER(3) - ROOT4, ZERO)
        NO3I   = MAX (WAER(4) - ROOT5, ZERO)
        CLI    = MAX (WAER(5) - ROOT4 - ROOT3, ZERO)
        CAI    = ZERO
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    !      IF (HI.LE.TINY) HI = SQRT(AKW)
    !      OHI   = AKW/HI
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV7 .AND. PSCONV1 .AND. PSCONV4 .AND. PSCONV3 &
             .AND. PSCONV5) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCU2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 IF (CLI <= TINY .AND. WAER(5) > TINY) THEN !No disslv Cl-;solid only
        DO 30 I=1,NIONS
            MOLAL(I) = ZERO
        30 END DO
        DO 40 I=1,NGASAQ
            GASAQ(I) = ZERO
        40 END DO
        CALL CALCU1A
    ELSE                                     ! OK, aqueous phase present
        A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
        A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
        A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
    
        GNH3    = NH4I/HI/A2
        GHNO3   = HI*NO3I/A3
        GHCL    = HI*CLI /A4
    
        GASAQ(1)= NH3AQ
        GASAQ(2)= CLAQ
        GASAQ(3)= NO3AQ
    
        CNH42S4 = ZERO
        CNH4NO3 = ZERO
        CNH4CL  = CHI4 - PSI4
        CNACL   = CHI3 - PSI3
        CNANO3  = CHI2 - PSI2
        CNA2SO4 = CHI1 - PSI1
        CMGSO4  = ZERO
        CK2SO4  = CHI7 - PSI7
        CCASO4  = MIN (WAER(6), WAER(2))
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCU2A*****************************************

    END SUBROUTINE CALCU2A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCU1
! *** CASE U1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4, MGSO4, NH4CL, NANO3, NACL, NH4NO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU1
    INCLUDE 'isrpia.inc'
    LOGICAL ::  EXAN, EXAC, EXSN, EXSC
    EXTERNAL CALCU1A, CALCU2A, CALCU3A, CALCU4A, CALCU5

! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************

    SCASE = 'U1 ; SUBCASE 1'
    CALL CALCU1A              ! SOLID
    SCASE = 'U1 ; SUBCASE 1'

    EXAN = CNH4NO3 > TINY    ! NH4NO3
    EXAC = CNH4CL > TINY    ! NH4CL
    EXSN = CNANO3 > TINY    ! NANO3
    EXSC = CNACL  > TINY    ! NACL

! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********

    IF (EXAN .OR. EXAC .OR. EXSC .OR. EXSN) THEN  ! *** WATER POSSIBLE
        IF (RH >= DRMM1) THEN
            SCASE = 'U1 ; SUBCASE 2'  ! MDRH
            CALL CALCMDRPII (RH, DRMM1, DRNH4NO3, CALCU1A, CALCU2A)
            SCASE = 'U1 ; SUBCASE 2'
        ENDIF
    
    ELSE IF ( .NOT. EXAN) THEN                   ! *** NH4NO3 = 0
        IF      (     EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMM2) THEN
                SCASE = 'U1 ; SUBCASE 3'
                CALL CALCMDRPII (RH, DRMM2, DRNANO3, CALCU1A, CALCU3A)
                SCASE = 'U1 ; SUBCASE 3'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR1) THEN
                SCASE = 'U1 ; SUBCASE 4'
                CALL CALCMDRPII (RH, DRMR1, DRNANO3, CALCU1A, CALCU3A)
                SCASE = 'U1 ; SUBCASE 4'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR2) THEN
                SCASE = 'U1 ; SUBCASE 5'
                CALL CALCMDRPII (RH, DRMR2, DRNACL, CALCU1A, CALCU3A) !, CALCR4A)
                SCASE = 'U1 ; SUBCASE 5'
            ENDIF

        ELSE IF ( .NOT. EXAC .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR3) THEN
                SCASE = 'U1 ; SUBCASE 6'
                CALL CALCMDRPII (RH, DRMR3, DRNANO3, CALCU1A, CALCU3A)
                SCASE = 'U1 ; SUBCASE 6'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR4) THEN
                SCASE = 'U1 ; SUBCASE 7'
                CALL CALCMDRPII (RH, DRMR4, DRNACL, CALCU1A, CALCU3A) !, CALCR4A)
                SCASE = 'U1 ; SUBCASE 7'
            ENDIF

        ELSE IF (     EXAC .AND. .NOT. EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR5) THEN
                SCASE = 'U1 ; SUBCASE 8'
                CALL CALCMDRPII (RH, DRMR5, DRNH4CL, CALCU1A, CALCU3A) !, CALCR5)
                SCASE = 'U1 ; SUBCASE 8'
            ENDIF

        ELSE IF (     EXAC .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR6) THEN
                SCASE = 'U1 ; SUBCASE 9'
                CALL CALCMDRPII (RH, DRMR6, DRNANO3, CALCU1A, CALCU3A)
                SCASE = 'U1 ; SUBCASE 9'
            ENDIF
        ENDIF
    
    ELSE IF ( .NOT. EXAC) THEN                   ! *** NH4CL  = 0
        IF      (     EXAN .AND.      EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR7) THEN
                SCASE = 'U1 ; SUBCASE 10'
                CALL CALCMDRPII (RH, DRMR7, DRNH4NO3, CALCU1A, CALCU2A)
                SCASE = 'U1 ; SUBCASE 10'
            ENDIF

        ELSE IF (     EXAN .AND. .NOT. EXSN .AND.      EXSC) THEN
            IF (RH >= DRMR8) THEN
                SCASE = 'U1 ; SUBCASE 11'
                CALL CALCMDRPII (RH, DRMR8, DRNH4NO3, CALCU1A, CALCU2A)
                SCASE = 'U1 ; SUBCASE 11'
            ENDIF

        ELSE IF (     EXAN .AND. .NOT. EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR9) THEN
                SCASE = 'U1 ; SUBCASE 12'
                CALL CALCMDRPII (RH, DRMR9, DRNH4NO3, CALCU1A, CALCU2A)
                SCASE = 'U1 ; SUBCASE 12'
            ENDIF

        ELSE IF (     EXAN .AND.      EXSN .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR10) THEN
                SCASE = 'U1 ; SUBCASE 13'
                CALL CALCMDRPII (RH, DRMR10, DRNH4NO3, CALCU1A, CALCU2A)
                SCASE = 'U1 ; SUBCASE 13'
            ENDIF
        ENDIF
    
    ELSE IF ( .NOT. EXSN) THEN                  ! *** NANO3  = 0
        IF      (     EXAN .AND.      EXAC .AND.      EXSC) THEN
            IF (RH >= DRMR11) THEN
                SCASE = 'U1 ; SUBCASE 14'
                CALL CALCMDRPII (RH, DRMR11, DRNH4NO3, CALCU1A, CALCU2A)
                SCASE = 'U1 ; SUBCASE 14'
            ENDIF

        ELSE IF (     EXAN .AND.      EXAC .AND. .NOT. EXSC) THEN
            IF (RH >= DRMR12) THEN
                SCASE = 'U1 ; SUBCASE 15'
                CALL CALCMDRPII (RH, DRMR12, DRNH4NO3, CALCU1A, CALCU2A)
                SCASE = 'U1 ; SUBCASE 15'
            ENDIF
        ENDIF
    
    ELSE IF ( .NOT. EXSC) THEN                  ! *** NACL   = 0
        IF      (     EXAN .AND.      EXAC .AND.      EXSN) THEN
            IF (RH >= DRMR13) THEN
                SCASE = 'U1 ; SUBCASE 16'
                CALL CALCMDRPII (RH, DRMR13, DRNH4NO3, CALCU1A, CALCU2A)
                SCASE = 'U1 ; SUBCASE 16'
            ENDIF
        ENDIF
    ENDIF

    RETURN


!      IF (RH.LT.DRMM1) THEN
!         SCASE = 'U1 ; SUBCASE 1'
!         CALL CALCU1A              ! SOLID PHASE ONLY POSSIBLE
!         SCASE = 'U1 ; SUBCASE 1'
!      ELSE
!         SCASE = 'U1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
!         CALL CALCMDRPII (RH, DRMM1, DRNH4NO3, CALCU1A, CALCU2A)
!         SCASE = 'U1 ; SUBCASE 2'
!         ENDIF
!C
!      RETURN
!C
! *** END OF SUBROUTINE CALCU1 ******************************************

    END SUBROUTINE CALCU1

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCU1A
! *** CASE U1A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0); CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
!     2. THERE IS ONLY A SOLID PHASE

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCU1A
    INCLUDE 'isrpia.inc'

! *** CALCULATE SOLIDS *************************************************

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

    CNH42S4 = ZERO

    CNANO3  = MIN (FRNA, WAER(4))
    FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
    FRNA    = MAX (FRNA-CNANO3, ZERO)

    CNACL   = MIN (FRNA, WAER(5))
    FRCL    = MAX (WAER(5)-CNACL, ZERO)
    FRNA    = MAX (FRNA-CNACL, ZERO)

    CNH4NO3 = MIN (FRNO3, WAER(3))
    FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
    FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)

    CNH4CL  = MIN (FRCL, FRNH3)
    FRCL    = MAX (FRCL-CNH4CL, ZERO)
    FRNH3   = MAX (FRNH3-CNH4CL, ZERO)

! *** OTHER PHASES ******************************************************

    WATER   = ZERO

    GNH3    = ZERO
    GHNO3   = ZERO
    GHCL    = ZERO

    RETURN

! *** END OF SUBROUTINE CALCU1A *****************************************

    END SUBROUTINE CALCU1A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW13
! *** CASE W13

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4
!     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW13
    INCLUDE 'isrpia.inc'

    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI11   = CCASO4

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

    CALL CALCMR           ! WATER

    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP

        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! ION CONCENTRATIONS
    
        NAI    = WAER(1)
        SO4I   = MAX (WAER(2) - WAER(6), ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        KI     = WAER(7)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW13')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

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

    RETURN

! *** END OF SUBROUTINE CALCW13 ******************************************

    END SUBROUTINE CALCW13
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW12
! *** CASE W12

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4
!     4. Completely dissolved: CA(NO3)2, CACL2, KNO3, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW12
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSI9O   =-GREAT                                 ! GREAT = 1.D10
    ROOT9   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI11   = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7))
            CC = WAER(7)*(WAER(2)-WAER(6)) + 0.25D0*WAER(7)*WAER(7)
            DD =-0.25*(WAER(7)*WAER(7)*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0, (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = WAER(3)
        NO3I   = WAER(4)
        CLI    = WAER(5)
        CAI    = ZERO
        NAI    = WAER(1)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW12')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = ZERO
    KCL     = ZERO
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW12 ******************************************

    END SUBROUTINE CALCW12
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW11
! *** CASE W11

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3
!     4. Completely dissolved: CA(NO3)2, CACL2, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW11
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT                                ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI11   = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13)
            CC = (WAER(7)-ROOT13)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13)**2.0
            DD =-0.25*((WAER(7)-ROOT13)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9,WAER(7)/2.0-ROOT13,(WAER(2)-WAER(6)),CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9)
            CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = WAER(3)
        NO3I   = MAX (WAER(4) - ROOT13, ZERO)
        CLI    = WAER(5)
        CAI    = ZERO
        NAI    = WAER(1)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW11')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = ZERO
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW11 ******************************************

    END SUBROUTINE CALCW11
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW10
! *** CASE W10

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4
!     4. Completely dissolved: CA(NO3)2, CACL2, KCL,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW10
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT                                ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A


    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI11   = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13)
            CC = (WAER(7)-ROOT13)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13)**2.0
            DD =-0.25*((WAER(7)-ROOT13)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9,WAER(7)/2.0-ROOT13,(WAER(2)-WAER(6)),CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9)
            CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = WAER(3)
        NO3I   = MAX (WAER(4) - ROOT13, ZERO)
        CLI    = WAER(5)
        CAI    = ZERO
        NAI    = WAER(1)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW10')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = ZERO
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW10 ******************************************

    END SUBROUTINE CALCW10
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW9
! *** CASE W9

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW9
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13, PSCONV14
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 
    PSCONV14= .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT
    PSI14O  =-GREAT                              ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO
    ROOT14  = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI14  = CKCL
    CHI11   = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        A14 = XK20 *(WATER/GAMA(20))**2.0      ! KCL       <==> K+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
            CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
            DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
            (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9 - ROOT14)
            CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! POTASSIUM CLORIDE
    
        IF (KI*CLI > A14) THEN
            BB     =-(WAER(5) + WAER(7) - 2.D0*ROOT9 - ROOT13)
            CC     = WAER(5)*(WAER(7) - 2.D0*ROOT9 - ROOT13) - A14
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT14A= 0.5D0*(-BB-DD)
            ROOT14B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT14A) THEN
                ROOT14 = ROOT14A
            ELSE
                ROOT14 = ROOT14B
            ENDIF
            ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
            PSI14  = CHI14-ROOT14
        ENDIF
        PSCONV14 = ABS(PSI14-PSI14O) <= EPS*PSI14O
        PSI14O   = PSI14
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = WAER(3)
        NO3I   = MAX (WAER(4) - ROOT13, ZERO)
        CLI    = MAX (WAER(5) - ROOT14, ZERO)
        CAI    = ZERO
        NAI    = WAER(1)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW9')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = ZERO
    CNACL   = ZERO
    CNANO3  = ZERO
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = CHI14 - PSI14
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW9 ******************************************

    END SUBROUTINE CALCW9
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW8
! *** CASE W8

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW8
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13, PSCONV14, PSCONV5
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 
    PSCONV14= .TRUE. 
    PSCONV5 = .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT
    PSI14O  =-GREAT
    PSI5O   =-GREAT                              ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO
    ROOT14  = ZERO
    ROOT5   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI14  = CKCL
    CHI5   = CNH4CL
    CHI11  = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        A14 = XK20 *(WATER/GAMA(20))**2.0      ! KCL       <==> K+
        A5  = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
            CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
            DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
            (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9 - ROOT14)
            CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! POTASSIUM CLORIDE
    
        IF (KI*CLI > A14) THEN
            BB     =-(WAER(5) - ROOT5 + WAER(7) - 2.D0*ROOT9 - ROOT13)
            CC     = (WAER(5)-ROOT5)*(WAER(7) - 2.D0*ROOT9 - ROOT13) - A14
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT14A= 0.5D0*(-BB-DD)
            ROOT14B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT14A) THEN
                ROOT14 = ROOT14A
            ELSE
                ROOT14 = ROOT14B
            ENDIF
            ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
            PSI14  = CHI14-ROOT14
        ENDIF
        PSCONV14 = ABS(PSI14-PSI14O) <= EPS*PSI14O
        PSI14O   = PSI14
    
    ! AMMONIUM CLORIDE
    
        IF (NH4I*CLI > A5) THEN
            BB     =-(WAER(5) + WAER(3) - ROOT14)
            CC     = (WAER(5)-ROOT14)*WAER(3) - A5
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT5A = 0.5D0*(-BB-DD)
            ROOT5B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT5A) THEN
                ROOT5 = ROOT5A
            ELSE
                ROOT5 = ROOT5B
            ENDIF
            ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
            PSI5  = CHI5-ROOT5
        ENDIF
        PSCONV5 = ABS(PSI5-PSI5O) <= EPS*PSI5O
        PSI5O   = PSI5
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = MAX (WAER(3) - ROOT5, ZERO)
        NO3I   = MAX (WAER(4) - ROOT13, ZERO)
        CLI    = MAX (WAER(5) - ROOT14 - ROOT5, ZERO)
        CAI    = ZERO
        NAI    = WAER(1)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW8')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI5 - PSI5
    CNACL   = ZERO
    CNANO3  = ZERO
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = CHI14 - PSI14
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW8 ******************************************

    END SUBROUTINE CALCW8
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW7
! *** CASE W7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NH4NO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW7
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 
    PSCONV14= .TRUE. 
    PSCONV5 = .TRUE. 
    PSCONV7 = .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT
    PSI14O  =-GREAT
    PSI5O  =-GREAT
    PSI7O  =-GREAT                            ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO
    ROOT14  = ZERO
    ROOT5   = ZERO
    ROOT7   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI14  = CKCL
    CHI5   = CNH4CL
    CHI7   = CNACL
    CHI11  = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        A14 = XK20 *(WATER/GAMA(20))**2.0      ! KCL       <==> K+
        A5  = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        A7  = XK8 *(WATER/GAMA(1))**2.0        ! NaCl      <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
            CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
            DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
            (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) + WAER(7) - 2.D0*ROOT9 - ROOT14)
            CC     = WAER(4)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! POTASSIUM CLORIDE
    
        IF (KI*CLI > A14) THEN
            BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
            CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT14A= 0.5D0*(-BB-DD)
            ROOT14B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT14A) THEN
                ROOT14 = ROOT14A
            ELSE
                ROOT14 = ROOT14B
            ENDIF
            ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
            PSI14  = CHI14-ROOT14
        ENDIF
        PSCONV14 = ABS(PSI14-PSI14O) <= EPS*PSI14O
        PSI14O   = PSI14
    
    ! AMMONIUM CLORIDE
    
        IF (NH4I*CLI > A5) THEN
            BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
            CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT5A = 0.5D0*(-BB-DD)
            ROOT5B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT5A) THEN
                ROOT5 = ROOT5A
            ELSE
                ROOT5 = ROOT5B
            ENDIF
            ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
            PSI5  = CHI5-ROOT5
        ENDIF
        PSCONV5 = ABS(PSI5-PSI5O) <= EPS*PSI5O
        PSI5O   = PSI5
    
    ! SODIUM CLORIDE
    
        IF (NAI*CLI > A7) THEN
            BB     =-(WAER(5) + WAER(1) - ROOT14 - ROOT5)
            CC     = (WAER(5) - ROOT14 - ROOT5)*WAER(1) - A7
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT7A = 0.5D0*(-BB-DD)
            ROOT7B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT7A) THEN
                ROOT7 = ROOT7A
            ELSE
                ROOT7 = ROOT7B
            ENDIF
            ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI7O) <= EPS*PSI7O
        PSI7O   = PSI7
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = MAX (WAER(3) - ROOT5, ZERO)
        NO3I   = MAX (WAER(4) - ROOT13, ZERO)
        CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
        CAI    = ZERO
        NAI    = MAX (WAER(1) - ROOT7, ZERO)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW7')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI5 - PSI5
    CNACL   = CHI7 - PSI7
    CNANO3  = ZERO
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = CHI14 - PSI14
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW7 ******************************************

    END SUBROUTINE CALCW7

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW6
! *** CASE W6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL, NANO3
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NH4NO3

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW6
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 
    PSCONV14= .TRUE. 
    PSCONV5 = .TRUE. 
    PSCONV7 = .TRUE. 
    PSCONV8 = .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT
    PSI14O  =-GREAT
    PSI5O   =-GREAT
    PSI7O   =-GREAT
    PSI8O   =-GREAT                     ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO
    ROOT14  = ZERO
    ROOT5   = ZERO
    ROOT7   = ZERO
    ROOT8   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI14  = CKCL
    CHI5   = CNH4CL
    CHI7   = CNACL
    CHI8   = CNANO3
    CHI11  = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        A14 = XK20 *(WATER/GAMA(20))**2.0      ! KCL       <==> K+
        A5  = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        A7  = XK8 *(WATER/GAMA(1))**2.0        ! NaCl      <==> Na+
        A8  = XK9 *(WATER/GAMA(3))**2.         ! NaNO3     <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
            CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
            DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
            (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
            CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! POTASSIUM CLORIDE
    
        IF (KI*CLI > A14) THEN
            BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
            CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT14A= 0.5D0*(-BB-DD)
            ROOT14B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT14A) THEN
                ROOT14 = ROOT14A
            ELSE
                ROOT14 = ROOT14B
            ENDIF
            ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
            PSI14  = CHI14-ROOT14
        ENDIF
        PSCONV14 = ABS(PSI14-PSI14O) <= EPS*PSI14O
        PSI14O   = PSI14
    
    ! AMMONIUM CLORIDE
    
        IF (NH4I*CLI > A5) THEN
            BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
            CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT5A = 0.5D0*(-BB-DD)
            ROOT5B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT5A) THEN
                ROOT5 = ROOT5A
            ELSE
                ROOT5 = ROOT5B
            ENDIF
            ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
            PSI5  = CHI5-ROOT5
        ENDIF
        PSCONV5 = ABS(PSI5-PSI5O) <= EPS*PSI5O
        PSI5O   = PSI5
    
    ! SODIUM CLORIDE
    
        IF (NAI*CLI > A7) THEN
            BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
            CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT7A = 0.5D0*(-BB-DD)
            ROOT7B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT7A) THEN
                ROOT7 = ROOT7A
            ELSE
                ROOT7 = ROOT7B
            ENDIF
            ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI7O) <= EPS*PSI7O
        PSI7O   = PSI7
    
    ! SODIUM NITRATE
    
        IF (NAI*NO3I > A8) THEN
            BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
            CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT8A = 0.5D0*(-BB-DD)
            ROOT8B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT8A) THEN
                ROOT8 = ROOT8A
            ELSE
                ROOT8 = ROOT8B
            ENDIF
            ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
            PSI8  = CHI8-ROOT8
        ENDIF
        PSCONV8 = ABS(PSI8-PSI8O) <= EPS*PSI8O
        PSI8O   = PSI8
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = MAX (WAER(3) - ROOT5, ZERO)
        NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
        CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
        CAI    = ZERO
        NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW6')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI5 - PSI5
    CNACL   = CHI7 - PSI7
    CNANO3  = CHI8 - PSI8
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = CHI14 - PSI14
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW6 ******************************************

    END SUBROUTINE CALCW6

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW5
! *** CASE W5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL, NANO3, NH4NO3
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW5
    INCLUDE 'isrpia.inc'

    EXTERNAL CALCW1A, CALCW6

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (WAER(4) > TINY)   THEN ! NO3 EXIST, WATER POSSIBLE
        SCASE = 'W5 ; SUBCASE 1'
        CALL CALCW5A
        SCASE = 'W5 ; SUBCASE 1'
    ELSE                                      ! NO3, CL NON EXISTANT
        SCASE = 'W1 ; SUBCASE 1'
        CALL CALCW1A
        SCASE = 'W1 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMP5) THEN        ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCW1A
            SCASE = 'W5 ; SUBCASE 2'
            RETURN
        ELSE
            SCASE = 'W5 ; SUBCASE 3'  ! MDRH REGION (CaSO4, K2SO4, KNO3, KCL, MGSO4,
        !                                                    NANO3, NACL, NH4NO3, NH4CL)
            CALL CALCMDRPII (RH, DRMP5, DRNH4NO3, CALCW1A, CALCW6)
            SCASE = 'W5 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCW5 ******************************************

    END SUBROUTINE CALCW5

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW5A
! *** CASE W5A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL,
!                          NANO3, NH4NO3
!     4. Completely dissolved: CA(NO3)2, CACL2, MG(NO3)2, MGCL2

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW5A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 
    PSCONV14= .TRUE. 
    PSCONV5 = .TRUE. 
    PSCONV7 = .TRUE. 
    PSCONV8 = .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT
    PSI14O  =-GREAT
    PSI5O   =-GREAT
    PSI7O   =-GREAT
    PSI8O   =-GREAT                ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO
    ROOT14  = ZERO
    ROOT5   = ZERO
    ROOT7   = ZERO
    ROOT8   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI14  = CKCL
    CHI5   = CNH4CL
    CHI7   = CNACL
    CHI8   = CNANO3
    CHI6   = CNH4NO3
    CHI11   = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        A14 = XK20 *(WATER/GAMA(20))**2.0      ! KCL       <==> K+
        A5  = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        A7  = XK8 *(WATER/GAMA(1))**2.0        ! NaCl      <==> Na+
        A8  = XK9 *(WATER/GAMA(3))**2.         ! NaNO3     <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
            CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
            DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
            (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
            CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! POTASSIUM CLORIDE
    
        IF (KI*CLI > A14) THEN
            BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
            CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT14A= 0.5D0*(-BB-DD)
            ROOT14B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT14A) THEN
                ROOT14 = ROOT14A
            ELSE
                ROOT14 = ROOT14B
            ENDIF
            ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
            PSI14  = CHI14-ROOT14
        ENDIF
        PSCONV14 = ABS(PSI14-PSI14O) <= EPS*PSI14O
        PSI14O   = PSI14
    
    ! AMMONIUM CLORIDE
    
        IF (NH4I*CLI > A5) THEN
            BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
            CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT5A = 0.5D0*(-BB-DD)
            ROOT5B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT5A) THEN
                ROOT5 = ROOT5A
            ELSE
                ROOT5 = ROOT5B
            ENDIF
            ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
            PSI5  = CHI5-ROOT5
        ENDIF
        PSCONV5 = ABS(PSI5-PSI5O) <= EPS*PSI5O
        PSI5O   = PSI5
    
    ! SODIUM CLORIDE
    
        IF (NAI*CLI > A7) THEN
            BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
            CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT7A = 0.5D0*(-BB-DD)
            ROOT7B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT7A) THEN
                ROOT7 = ROOT7A
            ELSE
                ROOT7 = ROOT7B
            ENDIF
            ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI7O) <= EPS*PSI7O
        PSI7O   = PSI7
    
    ! SODIUM NITRATE
    
        IF (NAI*NO3I > A8) THEN
            BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
            CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT8A = 0.5D0*(-BB-DD)
            ROOT8B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT8A) THEN
                ROOT8 = ROOT8A
            ELSE
                ROOT8 = ROOT8B
            ENDIF
            ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
            PSI8  = CHI8-ROOT8
        ENDIF
        PSCONV8 = ABS(PSI8-PSI8O) <= EPS*PSI8O
        PSI8O   = PSI8
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = MAX (WAER(3) - ROOT5, ZERO)
        NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
        CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
        CAI    = ZERO
        NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW5')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI5 - PSI5
    CNACL   = CHI7 - PSI7
    CNANO3  = CHI8 - PSI8
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = CHI14 - PSI14
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW5 ******************************************

    END SUBROUTINE CALCW5A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW4
! *** CASE W4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW4
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCW1A, CALCW5A

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (WAER(4) > TINY)   THEN ! NO3 EXIST, WATER POSSIBLE
        SCASE = 'W4 ; SUBCASE 1'
        CALL CALCW4A
        SCASE = 'W4 ; SUBCASE 1'
    ELSE                                      ! NO3, CL NON EXISTANT
        SCASE = 'W1 ; SUBCASE 1'
        CALL CALCW1A
        SCASE = 'W1 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMP4) THEN        ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCW1A
            SCASE = 'W4 ; SUBCASE 2'
            RETURN
        ELSE
            SCASE = 'W4 ; SUBCASE 3'  ! MDRH REGION (CaSO4, K2SO4, KNO3, KCL, MGSO4,
        !                                                    MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL)
            CALL CALCMDRPII (RH, DRMP4, DRMGNO32, CALCW1A, CALCW5A)
            SCASE = 'W4 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCW4 ******************************************

    END SUBROUTINE CALCW4

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW4A
! *** CASE W4A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL,
!                          NANO3, NH4NO3, MG(NO3)2
!     4. Completely dissolved: CA(NO3)2, CACL2, MGCL2

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW4A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 
    PSCONV14= .TRUE. 
    PSCONV5 = .TRUE. 
    PSCONV7 = .TRUE. 
    PSCONV8 = .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT
    PSI14O  =-GREAT
    PSI5O   =-GREAT
    PSI7O   =-GREAT
    PSI8O   =-GREAT                ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO
    ROOT14  = ZERO
    ROOT5   = ZERO
    ROOT7   = ZERO
    ROOT8   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI14  = CKCL
    CHI5   = CNH4CL
    CHI7   = CNACL
    CHI8   = CNANO3
    CHI6   = CNH4NO3
    CHI15  = CMGNO32
    CHI11   = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        A14 = XK20 *(WATER/GAMA(20))**2.0      ! KCL       <==> K+
        A5  = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        A7  = XK8 *(WATER/GAMA(1))**2.0        ! NaCl      <==> Na+
        A8  = XK9 *(WATER/GAMA(3))**2.         ! NaNO3     <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
            CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
            DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
            (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
            CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! POTASSIUM CLORIDE
    
        IF (KI*CLI > A14) THEN
            BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
            CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT14A= 0.5D0*(-BB-DD)
            ROOT14B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT14A) THEN
                ROOT14 = ROOT14A
            ELSE
                ROOT14 = ROOT14B
            ENDIF
            ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
            PSI14  = CHI14-ROOT14
        ENDIF
        PSCONV14 = ABS(PSI14-PSI14O) <= EPS*PSI14O
        PSI14O   = PSI14
    
    ! AMMONIUM CLORIDE
    
        IF (NH4I*CLI > A5) THEN
            BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
            CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT5A = 0.5D0*(-BB-DD)
            ROOT5B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT5A) THEN
                ROOT5 = ROOT5A
            ELSE
                ROOT5 = ROOT5B
            ENDIF
            ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
            PSI5  = CHI5-ROOT5
        ENDIF
        PSCONV5 = ABS(PSI5-PSI5O) <= EPS*PSI5O
        PSI5O   = PSI5
    
    ! SODIUM CLORIDE
    
        IF (NAI*CLI > A7) THEN
            BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
            CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT7A = 0.5D0*(-BB-DD)
            ROOT7B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT7A) THEN
                ROOT7 = ROOT7A
            ELSE
                ROOT7 = ROOT7B
            ENDIF
            ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI7O) <= EPS*PSI7O
        PSI7O   = PSI7
    
    ! SODIUM NITRATE
    
        IF (NAI*NO3I > A8) THEN
            BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
            CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT8A = 0.5D0*(-BB-DD)
            ROOT8B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT8A) THEN
                ROOT8 = ROOT8A
            ELSE
                ROOT8 = ROOT8B
            ENDIF
            ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
            PSI8  = CHI8-ROOT8
        ENDIF
        PSCONV8 = ABS(PSI8-PSI8O) <= EPS*PSI8O
        PSI8O   = PSI8
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = MAX (WAER(3) - ROOT5, ZERO)
        NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
        CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
        CAI    = ZERO
        NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW4')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI5 - PSI5
    CNACL   = CHI7 - PSI7
    CNANO3  = CHI8 - PSI8
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = CHI14 - PSI14
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW4A ******************************************

    END SUBROUTINE CALCW4A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW3
! *** CASE W3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW3
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCW1A, CALCW4A

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

!      IF (WAER(4).GT.TINY .AND. WAER(5).GT.TINY) THEN ! NO3,CL EXIST, WATER POSSIBLE
!         SCASE = 'W3 ; SUBCASE 1'
!         CALL CALCW3A
!         SCASE = 'W3 ; SUBCASE 1'
!      ELSE                                      ! NO3, CL NON EXISTANT
!         SCASE = 'W1 ; SUBCASE 1'
!         CALL CALCW1A
!         SCASE = 'W1 ; SUBCASE 1'
!      ENDIF

    CALL CALCW1A
          
    IF (WATER <= TINY) THEN
        IF (RH < DRMP3) THEN        ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCW1A
            SCASE = 'W3 ; SUBCASE 2'
            RETURN
        ELSE
            SCASE = 'W3 ; SUBCASE 3'  ! MDRH REGION (CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
        !                                                    MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL)
            CALL CALCMDRPII (RH, DRMP3, DRCANO32, CALCW1A, CALCW4A)
            SCASE = 'W3 ; SUBCASE 3'
        ENDIF
    ELSE                                      ! NO3, CL NON EXISTANT
        SCASE = 'W3 ; SUBCASE 1'
        CALL CALCW3A
        SCASE = 'W3 ; SUBCASE 1'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCW3 ******************************************

    END SUBROUTINE CALCW3

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW3A
! *** CASE W3A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL,
!                          NANO3, NH4NO3, CA(NO3)2, MG(NO3)2
!     4. Completely dissolved: CACL2, MGCL2

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW3A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 
    PSCONV14= .TRUE. 
    PSCONV5 = .TRUE. 
    PSCONV7 = .TRUE. 
    PSCONV8 = .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT
    PSI14O  =-GREAT
    PSI5O   =-GREAT
    PSI7O   =-GREAT
    PSI8O   =-GREAT                ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO
    ROOT14  = ZERO
    ROOT5   = ZERO
    ROOT7   = ZERO
    ROOT8   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI14  = CKCL
    CHI5   = CNH4CL
    CHI7   = CNACL
    CHI8   = CNANO3
    CHI6   = CNH4NO3
    CHI15  = CMGNO32
    CHI12  = CCANO32
    CHI11   = CCASO4
!C
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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        A14 = XK20 *(WATER/GAMA(20))**2.0      ! KCL       <==> K+
        A5  = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        A7  = XK8 *(WATER/GAMA(1))**2.0        ! NaCl      <==> Na+
        A8  = XK9 *(WATER/GAMA(3))**2.         ! NaNO3     <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
            CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
            DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
            (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
            CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, ZERO))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! POTASSIUM CLORIDE
    
        IF (KI*CLI > A14) THEN
            BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
            CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT14A= 0.5D0*(-BB-DD)
            ROOT14B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT14A) THEN
                ROOT14 = ROOT14A
            ELSE
                ROOT14 = ROOT14B
            ENDIF
            ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
            PSI14  = CHI14-ROOT14
        ENDIF
        PSCONV14 = ABS(PSI14-PSI14O) <= EPS*PSI14O
        PSI14O   = PSI14
    
    ! AMMONIUM CLORIDE
    
        IF (NH4I*CLI > A5) THEN
            BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
            CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT5A = 0.5D0*(-BB-DD)
            ROOT5B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT5A) THEN
                ROOT5 = ROOT5A
            ELSE
                ROOT5 = ROOT5B
            ENDIF
            ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
            PSI5  = CHI5-ROOT5
        ENDIF
        PSCONV5 = ABS(PSI5-PSI5O) <= EPS*PSI5O
        PSI5O   = PSI5
    
    ! SODIUM CLORIDE
    
        IF (NAI*CLI > A7) THEN
            BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
            CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT7A = 0.5D0*(-BB-DD)
            ROOT7B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT7A) THEN
                ROOT7 = ROOT7A
            ELSE
                ROOT7 = ROOT7B
            ENDIF
            ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI7O) <= EPS*PSI7O
        PSI7O   = PSI7
    
    ! SODIUM NITRATE
    
        IF (NAI*NO3I > A8) THEN
            BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
            CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT8A = 0.5D0*(-BB-DD)
            ROOT8B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT8A) THEN
                ROOT8 = ROOT8A
            ELSE
                ROOT8 = ROOT8B
            ENDIF
            ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
            PSI8  = CHI8-ROOT8
        ENDIF
        PSCONV8 = ABS(PSI8-PSI8O) <= EPS*PSI8O
        PSI8O   = PSI8
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = MAX (WAER(3) - ROOT5, ZERO)
        NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
        CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
        CAI    = ZERO
        NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW3')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI5 - PSI5
    CNACL   = CHI7 - PSI7
    CNANO3  = CHI8 - PSI8
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = CHI14 - PSI14
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW3A ******************************************

    END SUBROUTINE CALCW3A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW2
! *** CASE W2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

!     THERE ARE THREE REGIMES IN THIS CASE:
!     1. CACL2(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCL2A)
!     2. CACL2(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!     3. CACL2(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL

!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES W1A, W2B
!     RESPECTIVELY
! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================


    SUBROUTINE CALCW2
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCW1A, CALCW3A, CALCW4A, CALCW5A, CALCW6

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCW1A

! *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************

    IF (CCACL2 > TINY) THEN
        SCASE = 'W2 ; SUBCASE 1'
        CALL CALCW2A
        SCASE = 'W2 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMP2) THEN             ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCW1A
            SCASE = 'W2 ; SUBCASE 2'
        ELSE
            IF (CMGCL2 > TINY) THEN
                SCASE = 'W2 ; SUBCASE 3'    ! MDRH (CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4, MGCL2,
            !                                                  MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL)
                CALL CALCMDRPII (RH, DRMP2, DRMGCL2, CALCW1A, CALCW3A)
                SCASE = 'W2 ; SUBCASE 3'
            ENDIF
            IF (WATER <= TINY .AND. RH >= DRMP3 .AND. RH < DRMP4) THEN
                SCASE = 'W2 ; SUBCASE 4'    ! MDRH (CaSO4, K2SO4, KNO3, KCL, MGSO4, CANO32,
            !                                                  MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL)
                CALL CALCMDRPII (RH, DRMP3, DRCANO32, CALCW1A, CALCW4A)
                SCASE = 'W2 ; SUBCASE 4'
            ENDIF
            IF (WATER <= TINY .AND. RH >= DRMP4 .AND. RH < DRMP5) THEN
                SCASE = 'W2 ; SUBCASE 5'    ! MDRH (CaSO4, K2SO4, KNO3, KCL, MGSO4,
            !                                                  MGNO32, NANO3, NACL, NH4NO3, NH4CL)
                CALL CALCMDRPII (RH, DRMP4, DRMGNO32, CALCW1A, CALCW5A)
                SCASE = 'W2 ; SUBCASE 5'
            ENDIF
            IF (WATER <= TINY .AND. RH >= DRMP5) THEN
                SCASE = 'W2 ; SUBCASE 6'    ! MDRH (CaSO4, K2SO4, KNO3, KCL, MGSO4,
            !                                                  NANO3, NACL, NH4NO3, NH4CL)
                CALL CALCMDRPII (RH, DRMP5, DRNH4NO3, CALCW1A, CALCW6)
                SCASE = 'W2 ; SUBCASE 6'
            ELSE
                WATER = TINY
                DO 20 I=1,NIONS
                    MOLAL(I) = ZERO
                20 END DO
                CALL CALCW1A
                SCASE = 'W2 ; SUBCASE 2'
            ENDIF
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCW2 ******************************************

    END SUBROUTINE CALCW2

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW2A
! *** CASE W2A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL,
!                          NANO3, NH4NO3, CA(NO3)2, MG(NO3)2, MGCL2
!     4. Completely dissolved: CACL2

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW2A
    INCLUDE 'isrpia.inc'

    LOGICAL :: PSCONV9, PSCONV13, PSCONV14, PSCONV5, PSCONV7, PSCONV8
    DOUBLE PRECISION :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST    = .TRUE. 
    CALAIN  = .TRUE. 
    CALAOU  = .TRUE. 

    PSCONV9 = .TRUE. 
    PSCONV13= .TRUE. 
    PSCONV14= .TRUE. 
    PSCONV5 = .TRUE. 
    PSCONV7 = .TRUE. 
    PSCONV8 = .TRUE. 

    PSI9O   =-GREAT
    PSI13O  =-GREAT
    PSI14O  =-GREAT
    PSI5O   =-GREAT
    PSI7O   =-GREAT
    PSI8O   =-GREAT                ! GREAT = 1.D10

    ROOT9   = ZERO
    ROOT13  = ZERO
    ROOT14  = ZERO
    ROOT5   = ZERO
    ROOT7   = ZERO
    ROOT8   = ZERO

! *** CALCULATE INITIAL SOLUTION ***************************************

    CALL CALCW1A

    CHI9   = CK2SO4       ! SALTS
    CHI13  = CKNO3
    CHI10  = CMGSO4
    CHI14  = CKCL
    CHI5   = CNH4CL
    CHI7   = CNACL
    CHI8   = CNANO3
    CHI6   = CNH4NO3
    CHI15  = CMGNO32
    CHI12  = CCANO32
    CHI16  = CMGCL2
    CHI11   = CCASO4

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

    CALL CALCMR           ! WATER

    NAI    = WAER(1)      ! LIQUID CONCENTRATIONS
    SO4I   = MAX (WAER(2) - WAER(6), ZERO)
    NH4I   = WAER(3)
    NO3I   = WAER(4)
    CLI    = WAER(5)
    CAI    = WAER(6)
    KI     = WAER(7)
    MGI    = WAER(8)

    HSO4I  = ZERO
    NH3AQ  = ZERO
    NO3AQ  = ZERO
    CLAQ   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9  = XK17 *(WATER/GAMA(17))**3.0      ! K2SO4     <==> K+
        A13 = XK19 *(WATER/GAMA(19))**2.0      ! KNO3      <==> K+
        A14 = XK20 *(WATER/GAMA(20))**2.0      ! KCL       <==> K+
        A5  = XK14*(WATER/GAMA(6))**2.0        ! NH4Cl     <==> NH4+
        A7  = XK8 *(WATER/GAMA(1))**2.0        ! NaCl      <==> Na+
        A8  = XK9 *(WATER/GAMA(3))**2.         ! NaNO3     <==> Na+
        AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
    
    ! POTASSIUM SULFATE
    
        IF (KI*KI*SO4I > A9) THEN
            BB =-((WAER(2)-WAER(6)) + WAER(7) - ROOT13 - ROOT14)
            CC = (WAER(7)-ROOT13-ROOT14)*(WAER(2)-WAER(6)) + &
            &          0.25D0*(WAER(7)-ROOT13-ROOT14)**2.0
            DD =-0.25*((WAER(7)-ROOT13-ROOT14)**2.0*(WAER(2)-WAER(6)) - A9)
            CALL POLY3(BB, CC, DD, ROOT9, ISLV)
            IF (ISLV /= 0) ROOT9 = TINY
            ROOT9 = MIN (ROOT9, WAER(7)/2.0-ROOT13-ROOT14, &
            (WAER(2)-WAER(6)), CHI9)
            ROOT9 = MAX (ROOT9, ZERO)
            PSI9  = CHI9 - ROOT9
        ENDIF
        PSCONV9 = ABS(PSI9-PSI9O) <= EPS*PSI9O
        PSI9O   = PSI9
    
    ! POTASSIUM NITRATE
    
        IF (KI*NO3I > A13) THEN
            BB     =-(WAER(4) - ROOT8 + WAER(7) - 2.D0*ROOT9 - ROOT14)
            CC     = (WAER(4)-ROOT8)*(WAER(7) - 2.D0*ROOT9 - ROOT14) - A13
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT13A= 0.5D0*(-BB-DD)
            ROOT13B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT13A) THEN
                ROOT13 = ROOT13A
            ELSE
                ROOT13 = ROOT13B
            ENDIF
            ROOT13 = MIN(MAX(ROOT13, ZERO), CHI13)
            PSI13  = CHI13-ROOT13
        ENDIF
        PSCONV13 = ABS(PSI13-PSI13O) <= EPS*PSI13O
        PSI13O   = PSI13
    
    ! POTASSIUM CLORIDE
    
        IF (KI*CLI > A14) THEN
            BB     =-(WAER(5)-ROOT5-ROOT7 + WAER(7)-2.D0*ROOT9-ROOT13)
            CC     = (WAER(5)-ROOT5-ROOT7)*(WAER(7)-2.D0*ROOT9-ROOT13)-A14
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT14A= 0.5D0*(-BB-DD)
            ROOT14B= 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT14A) THEN
                ROOT14 = ROOT14A
            ELSE
                ROOT14 = ROOT14B
            ENDIF
            ROOT14 = MIN(MAX(ROOT14, ZERO), CHI14)
            PSI14  = CHI14-ROOT14
        ENDIF
        PSCONV14 = ABS(PSI14-PSI14O) <= EPS*PSI14O
        PSI14O   = PSI14
    
    ! AMMONIUM CLORIDE
    
        IF (NH4I*CLI > A5) THEN
            BB     =-(WAER(5) + WAER(3) - ROOT14 - ROOT7)
            CC     = (WAER(5) - ROOT14 - ROOT7)*WAER(3) - A5
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT5A = 0.5D0*(-BB-DD)
            ROOT5B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT5A) THEN
                ROOT5 = ROOT5A
            ELSE
                ROOT5 = ROOT5B
            ENDIF
            ROOT5 = MIN(MAX(ROOT5, ZERO), CHI5)
            PSI5  = CHI5-ROOT5
        ENDIF
        PSCONV5 = ABS(PSI5-PSI5O) <= EPS*PSI5O
        PSI5O   = PSI5
    
    ! SODIUM CLORIDE
    
        IF (NAI*CLI > A7) THEN
            BB     =-(WAER(5) + WAER(1) - ROOT8 - ROOT14 - ROOT5)
            CC     = (WAER(5) - ROOT14 - ROOT5)*(WAER(1)-ROOT8) - A7
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT7A = 0.5D0*(-BB-DD)
            ROOT7B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT7A) THEN
                ROOT7 = ROOT7A
            ELSE
                ROOT7 = ROOT7B
            ENDIF
            ROOT7 = MIN(MAX(ROOT7, ZERO), CHI7)
            PSI7  = CHI7-ROOT7
        ENDIF
        PSCONV7 = ABS(PSI7-PSI7O) <= EPS*PSI7O
        PSI7O   = PSI7
    
    ! SODIUM NITRATE
    
        IF (NAI*NO3I > A8) THEN
            BB     =-(WAER(4) - ROOT13 + WAER(1) - ROOT7)
            CC     = (WAER(4) - ROOT13)*(WAER(1)-ROOT7) - A8
            DD     = SQRT(MAX(BB*BB - 4.D0*CC, TINY))
            ROOT8A = 0.5D0*(-BB-DD)
            ROOT8B = 0.5D0*(-BB+DD)
            IF (ZERO <= ROOT8A) THEN
                ROOT8 = ROOT8A
            ELSE
                ROOT8 = ROOT8B
            ENDIF
            ROOT8 = MIN(MAX(ROOT8, ZERO), CHI8)
            PSI8  = CHI8-ROOT8
        ENDIF
        PSCONV8 = ABS(PSI8-PSI8O) <= EPS*PSI8O
        PSI8O   = PSI8
    
    ! ION CONCENTRATIONS ; CORRECTIONS
    
        KI     = MAX (WAER(7) - 2.D0*ROOT9 - ROOT13 - ROOT14, ZERO)
        SO4I   = MAX (WAER(2)-WAER(6) - ROOT9, ZERO)
        NH4I   = MAX (WAER(3) - ROOT5, ZERO)
        NO3I   = MAX (WAER(4) - ROOT13 - ROOT8, ZERO)
        CLI    = MAX (WAER(5) - ROOT14 - ROOT5 - ROOT7, ZERO)
        CAI    = ZERO
        NAI    = MAX (WAER(1) - ROOT7 - ROOT8, ZERO)
        MGI    = WAER(8)
    
    ! SOLUTION ACIDIC OR BASIC?
    
        GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I &
        - 2.D0*CAI - KI - 2.D0*MGI
        IF (GG > TINY) THEN                        ! H+ in excess
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
    
    ! UNDISSOCIATED SPECIES EQUILIBRIA
    
        IF (HI > OHI) THEN
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
        
        ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
        
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
        else
            del= zero
        ENDIF
        SO4I  = SO4I  - DEL
        HI    = HI    - DEL
        HSO4I = DEL
    !         IF (HI.LE.TINY) HI = SQRT(AKW)
        OHI   = AKW/HI
    
        IF (HI <= TINY) THEN
            HI = SQRT(AKW)
            OHI   = AKW/HI
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
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
    
    ! *** CALCULATE WATER **************************************************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            IF (PSCONV9 .AND. PSCONV13 .AND. PSCONV14 .AND. PSCONV5 &
             .AND. PSCONV7 .AND. PSCONV8) GOTO 20
        ENDIF
    10 END DO
! c      CALL PUSHERR (0002, 'CALCW2')    ! WARNING ERROR: NO CONVERGENCE

! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********

    20 A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
    A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
    A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

    GNH3    = NH4I/HI/A2
    GHNO3   = HI*NO3I/A3
    GHCL    = HI*CLI /A4

    GASAQ(1)= NH3AQ
    GASAQ(2)= CLAQ
    GASAQ(3)= NO3AQ

    CNH42S4 = ZERO
    CNH4NO3 = ZERO
    CNH4CL  = CHI5 - PSI5
    CNACL   = CHI7 - PSI7
    CNANO3  = CHI8 - PSI8
    CMGSO4  = ZERO
    CK2SO4  = CHI9 - PSI9
    CCASO4  = MIN (WAER(6), WAER(2))
    CCANO32 = ZERO
    CKNO3   = CHI13 - PSI13
    KCL     = CHI14 - PSI14
    CMGNO32 = ZERO
    CMGCL2  = ZERO
    CCACL2  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW2A ******************************************

    END SUBROUTINE CALCW2A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW1
! *** CASE W1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCP1A)

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCW1A, CALCW2A

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMP1) THEN
        SCASE = 'W1 ; SUBCASE 1'
        CALL CALCW1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'W1 ; SUBCASE 1'
    ELSE
        SCASE = 'W1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRPII (RH, DRMP1, DRCACL2, CALCW1A, CALCW2A)
        SCASE = 'W1 ; SUBCASE 2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCW1 ******************************************

    END SUBROUTINE CALCW1

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCW1A
! *** CASE W1A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCW1A
    INCLUDE 'isrpia.inc'

! *** CALCULATE SOLIDS **************************************************

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

! *** OTHER PHASES ******************************************************

    WATER   = ZERO

    GNH3    = ZERO
    GHNO3   = ZERO
    GHCL    = ZERO

    RETURN

! *** END OF SUBROUTINE CALCW1A *****************************************

    END SUBROUTINE CALCW1A
