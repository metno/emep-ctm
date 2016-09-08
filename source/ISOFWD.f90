! <ISOFWD.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4_10(3282)>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2016 met.no
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
! *** SUBROUTINE ISRP1F
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF
!     AN AMMONIUM-SULFATE AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!     THE AMBIENT RELATIVE HUMIDITY.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE ISRP1F (WI, RHI, TEMPI)
    INCLUDE 'isrpia.inc'
    DIMENSION WI(NCOMP)

! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************

    CALL INIT1 (WI, RHI, TEMPI)

! *** CALCULATE SULFATE RATIO *******************************************

    SULRAT = W(3)/W(2)

! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

! *** SULFATE POOR

    IF (2.0 <= SULRAT) THEN
        DC   = W(3) - 2.001D0*W(2)  ! For numerical stability
        W(3) = W(3) + MAX(-DC, ZERO)
    
        IF(METSTBL == 1) THEN
            SCASE = 'A2'
            CALL CALCA2                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH42S4) THEN
                SCASE = 'A1'
                CALL CALCA1              ! NH42SO4              ; case A1
            
            ELSEIF (DRNH42S4 <= RH) THEN
                SCASE = 'A2'
                CALL CALCA2              ! Only liquid          ; case A2
            ENDIF
        ENDIF
    
    ! *** SULFATE RICH (NO ACID)
    
    ELSEIF (1.0 <= SULRAT .AND. SULRAT < 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'B4'
            CALL CALCB4                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'B1'
                CALL CALCB1              ! NH4HSO4,LC,NH42SO4   ; case B1
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRLC) THEN
                SCASE = 'B2'
                CALL CALCB2              ! LC,NH42S4            ; case B2
            
            ELSEIF (DRLC <= RH .AND. RH < DRNH42S4) THEN
                SCASE = 'B3'
                CALL CALCB3              ! NH42S4               ; case B3
            
            ELSEIF (DRNH42S4 <= RH) THEN
                SCASE = 'B4'
                CALL CALCB4              ! Only liquid          ; case B4
            ENDIF
        ENDIF
        CALL CALCNH3
    
    ! *** SULFATE RICH (FREE ACID)
    
    ELSEIF (SULRAT < 1.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'C2'
            CALL CALCC2                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'C1'
                CALL CALCC1              ! NH4HSO4              ; case C1
            
            ELSEIF (DRNH4HS4 <= RH) THEN
                SCASE = 'C2'
                CALL CALCC2              ! Only liquid          ; case C2
            
            ENDIF
        ENDIF
        CALL CALCNH3
    ENDIF

! *** RETURN POINT

    RETURN

! *** END OF SUBROUTINE ISRP1F *****************************************

    END SUBROUTINE ISRP1F
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE ISRP2F
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF
!     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!     THE AMBIENT RELATIVE HUMIDITY.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE ISRP2F (WI, RHI, TEMPI)
    INCLUDE 'isrpia.inc'
    DIMENSION WI(NCOMP)

! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************

    CALL INIT2 (WI, RHI, TEMPI)

! *** CALCULATE SULFATE RATIO *******************************************

    SULRAT = W(3)/W(2)

! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

! *** SULFATE POOR

    IF (2.0 <= SULRAT) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'D3'
            CALL CALCD3                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4NO3) THEN
                SCASE = 'D1'
                CALL CALCD1              ! NH42SO4,NH4NO3       ; case D1
            
            ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNH42S4) THEN
                SCASE = 'D2'
                CALL CALCD2              ! NH42S4               ; case D2
            
            ELSEIF (DRNH42S4 <= RH) THEN
                SCASE = 'D3'
                CALL CALCD3              ! Only liquid          ; case D3
            ENDIF
        ENDIF
    
    ! *** SULFATE RICH (NO ACID)
    !     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES,
    !     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
    !     SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
    !     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
    
    ELSEIF (1.0 <= SULRAT .AND. SULRAT < 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'B4'
            CALL CALCB4                 ! Only liquid (metastable)
            SCASE = 'E4'
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'B1'
                CALL CALCB1              ! NH4HSO4,LC,NH42SO4   ; case E1
                SCASE = 'E1'
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRLC) THEN
                SCASE = 'B2'
                CALL CALCB2              ! LC,NH42S4            ; case E2
                SCASE = 'E2'
            
            ELSEIF (DRLC <= RH .AND. RH < DRNH42S4) THEN
                SCASE = 'B3'
                CALL CALCB3              ! NH42S4               ; case E3
                SCASE = 'E3'
            
            ELSEIF (DRNH42S4 <= RH) THEN
                SCASE = 'B4'
                CALL CALCB4              ! Only liquid          ; case E4
                SCASE = 'E4'
            ENDIF
        ENDIF
    
        CALL CALCNA                 ! HNO3(g) DISSOLUTION
    
    ! *** SULFATE RICH (FREE ACID)
    !     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES,
    !     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
    !     SUBROUTINE CALCC? IS CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
    !     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
    
    ELSEIF (SULRAT < 1.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'C2'
            CALL CALCC2                 ! Only liquid (metastable)
            SCASE = 'F2'
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'C1'
                CALL CALCC1              ! NH4HSO4              ; case F1
                SCASE = 'F1'
            
            ELSEIF (DRNH4HS4 <= RH) THEN
                SCASE = 'C2'
                CALL CALCC2              ! Only liquid          ; case F2
                SCASE = 'F2'
            ENDIF
        ENDIF
    
        CALL CALCNA                 ! HNO3(g) DISSOLUTION
    ENDIF

! *** RETURN POINT

    RETURN

! *** END OF SUBROUTINE ISRP2F *****************************************

    END SUBROUTINE ISRP2F
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE ISRP3F
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
!     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
!     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE ISRP3F (WI, RHI, TEMPI)
    INCLUDE 'isrpia.inc'
    DIMENSION WI(NCOMP)

! *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************

    WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
    WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3

! *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********

    IF (WI(1)+WI(2)+WI(4) <= 1d-10) THEN
        WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
        WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
    ENDIF

! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************

    CALL ISOINIT3 (WI, RHI, TEMPI)

! *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********

    REST = 2.D0*W(2) + W(4) + W(5)
    IF (W(1) > REST) THEN            ! NA > 2*SO4+CL+NO3 ?
        W(1) = (ONE-1D-6)*REST         ! Adjust Na amount
        CALL PUSHERR (0050, 'ISRP3F')  ! Warning error: Na adjusted
    ENDIF

! *** CALCULATE SULFATE & SODIUM RATIOS *********************************

    SULRAT = (W(1)+W(3))/W(2)
    SODRAT = W(1)/W(2)

! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

! *** SULFATE POOR ; SODIUM POOR

    IF (2.0 <= SULRAT .AND. SODRAT < 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'G5'
            CALL CALCG5                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4NO3) THEN
                SCASE = 'G1'
                CALL CALCG1              ! NH42SO4,NH4NO3,NH4CL,NA2SO4
            
            ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNH4CL) THEN
                SCASE = 'G2'
                CALL CALCG2              ! NH42SO4,NH4CL,NA2SO4
            
            ELSEIF (DRNH4CL <= RH  .AND. RH < DRNH42S4) THEN
                SCASE = 'G3'
                CALL CALCG3              ! NH42SO4,NA2SO4
            
            ELSEIF (DRNH42S4 <= RH  .AND. RH < DRNA2SO4) THEN
                SCASE = 'G4'
                CALL CALCG4              ! NA2SO4
            
            ELSEIF (DRNA2SO4 <= RH) THEN
                SCASE = 'G5'
                CALL CALCG5              ! Only liquid
            ENDIF
        ENDIF
    
    ! *** SULFATE POOR ; SODIUM RICH
    
    ELSE IF (SULRAT >= 2.0 .AND. SODRAT >= 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'H6'
            CALL CALCH6                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4NO3) THEN
                SCASE = 'H1'
                CALL CALCH1              ! NH4NO3,NH4CL,NA2SO4,NACL,NANO3
            
            ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNANO3) THEN
                SCASE = 'H2'
                CALL CALCH2              ! NH4CL,NA2SO4,NACL,NANO3
            
            ELSEIF (DRNANO3 <= RH  .AND. RH < DRNACL) THEN
                SCASE = 'H3'
                CALL CALCH3              ! NH4CL,NA2SO4,NACL
            
            ELSEIF (DRNACL <= RH   .AND. RH < DRNH4Cl) THEN
                SCASE = 'H4'
                CALL CALCH4              ! NH4CL,NA2SO4
            
            ELSEIF (DRNH4Cl <= RH .AND. RH < DRNA2SO4) THEN
                SCASE = 'H5'
                CALL CALCH5              ! NA2SO4
            
            ELSEIF (DRNA2SO4 <= RH) THEN
                SCASE = 'H6'
                CALL CALCH6              ! NO SOLID
            ENDIF
        ENDIF
    
    ! *** SULFATE RICH (NO ACID)
    
    ELSEIF (1.0 <= SULRAT .AND. SULRAT < 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'I6'
            CALL CALCI6                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'I1'
                CALL CALCI1              ! NA2SO4,(NH4)2SO4,NAHSO4,NH4HSO4,LC
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRNAHSO4) THEN
                SCASE = 'I2'
                CALL CALCI2              ! NA2SO4,(NH4)2SO4,NAHSO4,LC
            
            ELSEIF (DRNAHSO4 <= RH .AND. RH < DRLC) THEN
                SCASE = 'I3'
                CALL CALCI3              ! NA2SO4,(NH4)2SO4,LC
            
            ELSEIF (DRLC <= RH     .AND. RH < DRNH42S4) THEN
                SCASE = 'I4'
                CALL CALCI4              ! NA2SO4,(NH4)2SO4
            
            ELSEIF (DRNH42S4 <= RH .AND. RH < DRNA2SO4) THEN
                SCASE = 'I5'
                CALL CALCI5              ! NA2SO4
            
            ELSEIF (DRNA2SO4 <= RH) THEN
                SCASE = 'I6'
                CALL CALCI6              ! NO SOLIDS
            ENDIF
        ENDIF
    
        CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3                !                NH3
    
    ! *** SULFATE RICH (FREE ACID)
    
    ELSEIF (SULRAT < 1.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'J3'
            CALL CALCJ3                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4HS4) THEN
                SCASE = 'J1'
                CALL CALCJ1              ! NH4HSO4,NAHSO4
            
            ELSEIF (DRNH4HS4 <= RH .AND. RH < DRNAHSO4) THEN
                SCASE = 'J2'
                CALL CALCJ2              ! NAHSO4
            
            ELSEIF (DRNAHSO4 <= RH) THEN
                SCASE = 'J3'
                CALL CALCJ3
            ENDIF
        ENDIF
    
        CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3                !                NH3
    ENDIF

! *** RETURN POINT

    RETURN

! *** END OF SUBROUTINE ISRP3F *****************************************

    END SUBROUTINE ISRP3F

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE ISRP4F
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
!     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM-CALCIUM-POTASSIUM-MAGNESIUM
!     AEROSOL SYSTEM.
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
!     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE ISRP4F (WI, RHI, TEMPI)
    INCLUDE 'isrpia.inc'
    DIMENSION WI(NCOMP)
    real :: NAFRI, NO3FRI

! *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************

!      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
!      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3

! *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********

!      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
!         WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
!         WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
!      ENDIF

! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************

    CALL INIT4 (WI, RHI, TEMPI)

! *** CHECK IF TOO MUCH SODIUM+CRUSTALS ; ADJUST AND ISSUE ERROR MESSAGE

    REST = 2.D0*W(2) + W(4) + W(5)

    IF (W(1)+W(6)+W(7)+W(8) > REST) THEN
    
        CCASO4I  = MIN (W(2),W(6))
        FRSO4I   = MAX (W(2) - CCASO4I, ZERO)
        CAFRI    = MAX (W(6) - CCASO4I, ZERO)
        CCANO32I = MIN (CAFRI, 0.5D0*W(4))
        CAFRI    = MAX (CAFRI - CCANO32I, ZERO)
        NO3FRI   = MAX (W(4) - 2.D0*CCANO32I, ZERO)
        CCACL2I  = MIN (CAFRI, 0.5D0*W(5))
        CLFRI    = MAX (W(5) - 2.D0*CCACL2I, ZERO)
        REST1    = 2.D0*FRSO4I + NO3FRI + CLFRI
    
        CNA2SO4I = MIN (FRSO4I, 0.5D0*W(1))
        FRSO4I   = MAX (FRSO4I - CNA2SO4I, ZERO)
        NAFRI    = MAX (W(1) - 2.D0*CNA2SO4I, ZERO)
        CNACLI   = MIN (NAFRI, CLFRI)
        NAFRI    = MAX (NAFRI - CNACLI, ZERO)
        CLFRI    = MAX (CLFRI - CNACLI, ZERO)
        CNANO3I  = MIN (NAFRI, NO3FRI)
        NO3FR    = MAX (NO3FRI - CNANO3I, ZERO)
        REST2    = 2.D0*FRSO4I + NO3FRI + CLFRI
    
        CMGSO4I  = MIN (FRSO4I, W(8))
        FRMGI    = MAX (W(8) - CMGSO4I, ZERO)
        FRSO4I   = MAX (FRSO4I - CMGSO4I, ZERO)
        CMGNO32I = MIN (FRMGI, 0.5D0*NO3FRI)
        FRMGI    = MAX (FRMGI - CMGNO32I, ZERO)
        NO3FRI   = MAX (NO3FRI - 2.D0*CMGNO32I, ZERO)
        CMGCL2I  = MIN (FRMGI, 0.5D0*CLFRI)
        CLFRI    = MAX (CLFRI - 2.D0*CMGCL2I, ZERO)
        REST3    = 2.D0*FRSO4I + NO3FRI + CLFRI
    
        IF (W(6) > REST) THEN                       ! Ca > 2*SO4+CL+NO3 ?
            W(6) = (ONE-1D-6)*REST              ! Adjust Ca amount
            W(1)= ZERO                          ! Adjust Na amount
            W(7)= ZERO                          ! Adjust K amount
            W(8)= ZERO                          ! Adjust Mg amount
            CALL PUSHERR (0051, 'ISRP4F')       ! Warning error: Ca, Na, K, Mg in excess
        
        ELSE IF (W(1) > REST1) THEN                 ! Na > 2*FRSO4+FRCL+FRNO3 ?
            W(1) = (ONE-1D-6)*REST1             ! Adjust Na amount
            W(7)= ZERO                          ! Adjust K amount
            W(8)= ZERO                          ! Adjust Mg amount
            CALL PUSHERR (0052, 'ISRP4F')       ! Warning error: Na, K, Mg in excess
        
        ELSE IF (W(8) > REST2) THEN                 ! Mg > 2*FRSO4+FRCL+FRNO3 ?
            W(8) = (ONE-1D-6)*REST2             ! Adjust Mg amount
            W(7)= ZERO                          ! Adjust K amount
            CALL PUSHERR (0053, 'ISRP4F')       ! Warning error: K, Mg in excess
        
        ELSE IF (W(7) > REST3) THEN                 ! K > 2*FRSO4+FRCL+FRNO3 ?
            W(7) = (ONE-1D-6)*REST3             ! Adjust K amount
            CALL PUSHERR (0054, 'ISRP4F')       ! Warning error: K in excess
        ENDIF
    ENDIF

! *** CALCULATE RATIOS *************************************************

    SO4RAT  = (W(1)+W(3)+W(6)+W(7)+W(8))/W(2)
    CRNARAT = (W(1)+W(6)+W(7)+W(8))/W(2)
    CRRAT   = (W(6)+W(7)+W(8))/W(2)

! *** FIND CALCULATION REGIME FROM (SO4RAT, CRNARAT, CRRAT, RRH) ********

! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) POOR: R(Cr+Na)<2

    IF (2.0 <= SO4RAT .AND. CRNARAT < 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'O7'
            CALL CALCO7                 ! Only liquid (metastable)
        ELSE
        
            IF (RH < DRNH4NO3) THEN
                SCASE = 'O1'
                CALL CALCO1              ! CaSO4, NH4NO3, NH4CL, (NH4)2SO4, MGSO4, NA2SO4, K2SO4
            
            ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNH4CL) THEN
                SCASE = 'O2'
                CALL CALCO2              ! CaSO4, NH4CL, (NH4)2SO4, MGSO4, NA2SO4, K2SO4
            
            ELSEIF (DRNH4CL <= RH  .AND. RH < DRNH42S4) THEN
                SCASE = 'O3'
                CALL CALCO3              ! CaSO4, (NH4)2SO4, MGSO4, NA2SO4, K2SO4
            
            ELSEIF (DRNH42S4 <= RH .AND. RH < DRMGSO4) THEN
                SCASE = 'O4'
                CALL CALCO4              ! CaSO4, MGSO4, NA2SO4, K2SO4
            
            ELSEIF (DRMGSO4 <= RH .AND. RH < DRNA2SO4) THEN
                SCASE = 'O5'
                CALL CALCO5              ! CaSO4, NA2SO4, K2SO4
            
            ELSEIF (DRNA2SO4 <= RH .AND. RH < DRK2SO4) THEN
                SCASE = 'O6'
                CALL CALCO6              ! CaSO4, K2SO4
            
            ELSEIF (DRK2SO4 <= RH) THEN
                SCASE = 'O7'
                CALL CALCO7              ! CaSO4
            ENDIF
        ENDIF
    
    ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
    
    ELSEIF (SO4RAT >= 2.0 .AND. CRNARAT >= 2.0) THEN
    
        IF (CRRAT <= 2.0) THEN
        
            IF(METSTBL == 1) THEN
                SCASE = 'M8'
                CALL CALCM8                 ! Only liquid (metastable)
            ELSE
            
                IF (RH < DRNH4NO3) THEN
                    SCASE = 'M1'
                    CALL CALCM1            ! CaSO4, NH4NO3, NH4CL, MGSO4, NA2SO4, K2SO4, NACL, NANO3
                
                ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNANO3) THEN
                    SCASE = 'M2'
                    CALL CALCM2            ! CaSO4, NH4CL, MGSO4, NA2SO4, K2SO4, NACL, NANO3
                
                ELSEIF (DRNANO3 <= RH  .AND. RH < DRNACL) THEN
                    SCASE = 'M3'
                    CALL CALCM3            ! CaSO4, NH4CL, MGSO4, NA2SO4, K2SO4, NACL
                
                ELSEIF (DRNACL <= RH   .AND. RH < DRNH4Cl) THEN
                    SCASE = 'M4'
                    CALL CALCM4            ! CaSO4, NH4CL, MGSO4, NA2SO4, K2SO4
                
                ELSEIF (DRNH4Cl <= RH .AND. RH < DRMGSO4) THEN
                    SCASE = 'M5'
                    CALL CALCM5            ! CaSO4, MGSO4, NA2SO4, K2SO4
                
                ELSEIF (DRMGSO4 <= RH .AND. RH < DRNA2SO4) THEN
                    SCASE = 'M6'
                    CALL CALCM6            ! CaSO4, NA2SO4, K2SO4
                
                ELSEIF (DRNA2SO4 <= RH .AND. RH < DRK2SO4) THEN
                    SCASE = 'M7'
                    CALL CALCM7            ! CaSO4, K2SO4
                
                ELSEIF (DRK2SO4 <= RH) THEN
                    SCASE = 'M8'
                    CALL CALCM8            ! CaSO4
                ENDIF
            ENDIF
        !        CALL CALCHCO3
        
        ! *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
        
        ELSEIF (CRRAT > 2.0) THEN
        
            IF(METSTBL == 1) THEN
                SCASE = 'P13'
                CALL CALCP13                 ! Only liquid (metastable)
            ELSE
            
                IF (RH < DRCACL2) THEN
                    SCASE = 'P1'
                    CALL CALCP1             ! CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
                !                                    ! MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRCACL2 <= RH .AND. RH < DRMGCL2) THEN
                    SCASE = 'P2'
                    CALL CALCP2            ! CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
                !                                   ! MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRMGCL2 <= RH  .AND. RH < DRCANO32) THEN
                    SCASE = 'P3'
                    CALL CALCP3            ! CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
                !                                   ! MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRCANO32 <= RH   .AND. RH < DRMGNO32) THEN
                    SCASE = 'P4'
                    CALL CALCP4            ! CaSO4, K2SO4, KNO3, KCL, MGSO4,
                !                                   ! MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRMGNO32 <= RH .AND. RH < DRNH4NO3) THEN
                    SCASE = 'P5'
                    CALL CALCP5            ! CaSO4, K2SO4, KNO3, KCL, MGSO4,
                !                                   ! NANO3, NACL, NH4NO3, NH4CL
                
                ELSEIF (DRNH4NO3 <= RH .AND. RH < DRNANO3) THEN
                    SCASE = 'P6'
                    CALL CALCP6            ! CaSO4, K2SO4, KNO3, KCL, MGSO4, NANO3, NACL, NH4CL
                
                ELSEIF (DRNANO3 <= RH .AND. RH < DRNACL) THEN
                    SCASE = 'P7'
                    CALL CALCP7            ! CaSO4, K2SO4, KNO3, KCL, MGSO4, NACL, NH4CL
                
                ELSEIF (DRNACL <= RH .AND. RH < DRNH4CL) THEN
                    SCASE = 'P8'
                    CALL CALCP8            ! CaSO4, K2SO4, KNO3, KCL, MGSO4, NH4CL
                
                ELSEIF (DRNH4CL <= RH .AND. RH < DRKCL) THEN
                    SCASE = 'P9'
                    CALL CALCP9            ! CaSO4, K2SO4, KNO3, KCL, MGSO4
                
                ELSEIF (DRKCL <= RH .AND. RH < DRMGSO4) THEN
                    SCASE = 'P10'
                    CALL CALCP10            ! CaSO4, K2SO4, KNO3, MGSO4
                
                ELSEIF (DRMGSO4 <= RH .AND. RH < DRKNO3) THEN
                    SCASE = 'P11'
                    CALL CALCP11            ! CaSO4, K2SO4, KNO3
                
                ELSEIF (DRKNO3 <= RH .AND. RH < DRK2SO4) THEN
                    SCASE = 'P12'
                    CALL CALCP12            ! CaSO4, K2SO4
                
                ELSEIF (DRK2SO4 <= RH) THEN
                    SCASE = 'P13'
                    CALL CALCP13            ! CaSO4
                ENDIF
            ENDIF
        !        CALL CALCHCO3
        ENDIF
    
    ! *** SULFATE RICH (NO ACID): 1<Rso4<2;
    
    ELSEIF (1.0 <= SO4RAT .AND. SO4RAT < 2.0) THEN
    
        IF(METSTBL == 1) THEN
            SCASE = 'L9'
            CALL CALCL9                ! Only liquid (metastable)
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
    
        CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3                !                NH3
    
    ! *** SULFATE SUPER RICH (FREE ACID): Rso4<1;
    
    ELSEIF (SO4RAT < 1.0) THEN
    
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
    
        CALL CALCNHA                  ! MINOR SPECIES: HNO3, HCl
        CALL CALCNH3                  !                NH3
    
    ENDIF

    RETURN
    END SUBROUTINE ISRP4F

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCA2
! *** CASE A2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT >= 2.0)
!     2. LIQUID AEROSOL PHASE ONLY POSSIBLE

!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS X, THE
!     AMOUNT OF HYDROGEN IONS (H+) FOUND IN THE LIQUID PHASE.
!     FOR EACH ESTIMATION OF H+, FUNCTION FUNCB2A CALCULATES THE
!     CONCENTRATION OF IONS FROM THE NH3(GAS) - NH4+(LIQ) EQUILIBRIUM.
!     ELECTRONEUTRALITY IS USED AS THE OBJECTIVE FUNCTION.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCA2
    INCLUDE 'isrpia.inc'

! *** SETUP PARAMETERS ************************************************

    CALAOU    = .TRUE.       ! Outer loop activity calculation flag
    OMELO     = TINY        ! Low  limit: SOLUTION IS VERY BASIC
    OMEHI     = 2.0D0*W(2)  ! High limit: FROM NH4+ -> NH3(g) + H+(aq)

! *** CALCULATE WATER CONTENT *****************************************

    MOLAL(5) = W(2)
    MOLAL(6) = ZERO
    CALL CALCMR

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = OMEHI
    Y1 = FUNCA2 (X1)
    IF (ABS(Y1) <= EPS) RETURN

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (OMEHI-OMELO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, OMELO)
        Y2 = FUNCA2 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO
    IF (ABS(Y2) <= EPS) THEN
        RETURN
    ELSE
        CALL PUSHERR (0001, 'CALCA2')    ! WARNING ERROR: NO SOLUTION
        RETURN
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCA2 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCA2')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCA2 (X3)
    RETURN

! *** END OF SUBROUTINE CALCA2 ****************************************

    END SUBROUTINE CALCA2



!=======================================================================

! *** ISORROPIA CODE
! *** FUNCTION FUNCA2
! *** CASE A2
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE A2 ;
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA2.

!=======================================================================

    real  FUNCTION FUNCA2 (OMEGI)
    INCLUDE 'isrpia.inc'
    real :: LAMDA

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    PSI    = W(2)         ! INITIAL AMOUNT OF (NH4)2SO4 IN SOLUTION

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A1    = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
        A2    = XK2*R*TEMP/XKW*(GAMA(8)/GAMA(9))**2.
        A3    = XKW*RH*WATER*WATER
    
        LAMDA = PSI/(A1/OMEGI+ONE)
        ZETA  = A3/OMEGI
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
        MOLAL (1) = OMEGI                                        ! HI
        MOLAL (5) = MAX(PSI-LAMDA,TINY)                          ! SO4I
        MOLAL (3) = MAX(W(3)/(ONE/A2/OMEGI + ONE), 2.*MOLAL(5))  ! NH4I
        MOLAL (6) = LAMDA                                        ! HSO4I
        GNH3      = MAX (W(3)-MOLAL(3), TINY)                    ! NH3GI
        COH       = ZETA                                         ! OHI
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 DENOM = (2.0*MOLAL(5)+MOLAL(6))
    FUNCA2= (MOLAL(3)/DENOM - ONE) + MOLAL(1)/DENOM
    RETURN

! *** END OF FUNCTION FUNCA2 ********************************************

  END FUNCTION FUNCA2
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCA1
! *** CASE A1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4

!     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE SOLID (NH4)2SO4
!     IS CALCULATED FROM THE SULFATES. THE EXCESS AMMONIA REMAINS IN
!     THE GAS PHASE.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCA1
    INCLUDE 'isrpia.inc'

    CNH42S4 = W(2)
    GNH3    = MAX (W(3)-2.0*CNH42S4, ZERO)
    RETURN

! *** END OF SUBROUTINE CALCA1 ******************************************

    END SUBROUTINE CALCA1



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB4
! *** CASE B4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. LIQUID AEROSOL PHASE ONLY POSSIBLE

!     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
!     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
!     AND THAT CALCULATED FROM ELECTRONEUTRALITY.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB4
    INCLUDE 'isrpia.inc'

! *** SOLVE EQUATIONS **************************************************

    FRST       = .TRUE. 
    CALAIN     = .TRUE. 
    CALAOU     = .TRUE. 

! *** CALCULATE WATER CONTENT ******************************************

    CALL CALCB1A         ! GET DRY SALT CONTENT, AND USE FOR WATER.
    MOLALR(13) = CLC
    MOLALR(9)  = CNH4HS4
    MOLALR(4)  = CNH42S4
    CLC        = ZERO
    CNH4HS4    = ZERO
    CNH42S4    = ZERO
    WATER      = MOLALR(13)/M0(13)+MOLALR(9)/M0(9)+MOLALR(4)/M0(4)

    MOLAL(3)   = W(3)   ! NH4I

    DO 20 I=1,NSWEEP
        AK1   = XK1*((GAMA(8)/GAMA(7))**2.)*(WATER/GAMA(7))
        BET   = W(2)
        GAM   = MOLAL(3)
    
        BB    = BET + AK1 - GAM
        CC    =-AK1*BET
        DD    = BB*BB - 4.D0*CC
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
        MOLAL (5) = MAX(TINY,MIN(0.5*(-BB + SQRT(DD)), W(2))) ! SO4I
        MOLAL (6) = MAX(TINY,MIN(W(2)-MOLAL(5),W(2)))         ! HSO4I
        MOLAL (1) = MAX(TINY,MIN(AK1*MOLAL(6)/MOLAL(5),W(2))) ! HI
        CALL CALCMR                                           ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF ( .NOT. CALAIN) GOTO 30
        CALL CALCACT
    20 END DO

    30 RETURN

! *** END OF SUBROUTINE CALCB4 ******************************************

    END SUBROUTINE CALCB4
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3
! *** CASE B3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: (NH4)2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB3
    INCLUDE 'isrpia.inc'

! *** CALCULATE EQUIVALENT AMOUNT OF HSO4 AND SO4 ***********************

    X = MAX(2*W(2)-W(3), ZERO)   ! Equivalent NH4HSO4
    Y = MAX(W(3)  -W(2), ZERO)   ! Equivalent NH42SO4

! *** CALCULATE SPECIES ACCORDING TO RELATIVE ABUNDANCE OF HSO4 *********

    IF (X < Y) THEN             ! LC is the MIN (x,y)
        SCASE   = 'B3 ; SUBCASE 1'
        TLC     = X
        TNH42S4 = Y-X
        CALL CALCB3A (TLC,TNH42S4)      ! LC + (NH4)2SO4
    ELSE
        SCASE   = 'B3 ; SUBCASE 2'
        TLC     = Y
        TNH4HS4 = X-Y
        CALL CALCB3B (TLC,TNH4HS4)      ! LC + NH4HSO4
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCB3 ******************************************

    END SUBROUTINE CALCB3


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3A
! *** CASE B3 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: (NH4)2SO4

!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
!     AMOUNT OF SOLID (NH4)2SO4 DISSOLVED IN THE LIQUID PHASE.
!     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB3A CALCULATES THE
!     AMOUNT OF H+ PRODUCED (BASED ON THE SO4 RELEASED INTO THE
!     SOLUTION). THE SOLUBILITY PRODUCT OF (NH4)2SO4 IS USED AS THE
!     OBJECTIVE FUNCTION.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB3A (TLC, TNH42S4)
    INCLUDE 'isrpia.inc'

    CALAOU = .TRUE.         ! Outer loop activity calculation flag
    ZLO    = ZERO           ! MIN DISSOLVED (NH4)2SO4
    ZHI    = TNH42S4        ! MAX DISSOLVED (NH4)2SO4

! *** INITIAL VALUES FOR BISECTION (DISSOLVED (NH4)2SO4 ****************

    Z1 = ZLO
    Y1 = FUNCB3A (Z1, TLC, TNH42S4)
    IF (ABS(Y1) <= EPS) RETURN
    YLO= Y1

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************

    DZ = (ZHI-ZLO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        Z2 = Z1+DZ
        Y2 = FUNCB3A (Z2, TLC, TNH42S4)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        Z1 = Z2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION FOUND

    YHI= Y1                      ! Save Y-value at HI position
    IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        RETURN
    
    ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
    
    ELSE IF (YLO < ZERO .AND. YHI < ZERO) THEN
        Z1 = ZHI
        Z2 = ZHI
        GOTO 40
    
    ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
    
    ELSE IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Z1 = ZLO
        Z2 = ZLO
        GOTO 40
    ELSE
        CALL PUSHERR (0001, 'CALCB3A')    ! WARNING ERROR: NO SOLUTION
        RETURN
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        Z3 = 0.5*(Z1+Z2)
        Y3 = FUNCB3A (Z3, TLC, TNH42S4)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            Z2    = Z3
        ELSE
            Y1    = Y3
            Z1    = Z3
        ENDIF
        IF (ABS(Z2-Z1) <= EPS*Z1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCB3A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN ************************************************

    40 ZK = 0.5*(Z1+Z2)
    Y3 = FUNCB3A (ZK, TLC, TNH42S4)

    RETURN

! *** END OF SUBROUTINE CALCB3A ******************************************

    END SUBROUTINE CALCB3A



!=======================================================================

! *** ISORROPIA CODE
! *** FUNCTION FUNCB3A
! *** CASE B3 ; SUBCASE 1
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B3
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA3.

!=======================================================================

    real FUNCTION FUNCB3A (ZK, Y, X)
    INCLUDE 'isrpia.inc'
    real :: KK

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    DO 20 I=1,NSWEEP
        GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
        DD    = SQRT( (ZK+GRAT1+Y)**2. + 4.0*Y*GRAT1)
        KK    = 0.5*(-(ZK+GRAT1+Y) + DD )
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
        MOLAL (1) = KK                ! HI
        MOLAL (5) = KK+ZK+Y           ! SO4I
        MOLAL (6) = MAX (Y-KK, TINY)  ! HSO4I
        MOLAL (3) = 3.0*Y+2*ZK        ! NH4I
        CNH42S4   = X-ZK              ! Solid (NH4)2SO4
        CALL CALCMR                   ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 30
        ENDIF
    20 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

! C30    FUNCB3A= ( SO4I*NH4I**2.0 )/( XK7*(WATER/GAMA(4))**3.0 )
    30 FUNCB3A= MOLAL(5)*MOLAL(3)**2.0
    FUNCB3A= FUNCB3A/(XK7*(WATER/GAMA(4))**3.0) - ONE
    RETURN

! *** END OF FUNCTION FUNCB3A ********************************************

  END FUNCTION FUNCB3A



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3B
! *** CASE B3 ; SUBCASE 2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. LIQUID PHASE ONLY IS POSSIBLE

!     SPECIATION CALCULATIONS IS BASED ON THE HSO4 <--> SO4 EQUILIBRIUM.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB3B (Y, X)
    INCLUDE 'isrpia.inc'
    real :: KK

    CALAOU = .FALSE.        ! Outer loop activity calculation flag
    FRST   = .FALSE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 20 I=1,NSWEEP
        GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
        DD    = SQRT( (GRAT1+Y)**2. + 4.0*(X+Y)*GRAT1)
        KK    = 0.5*(-(GRAT1+Y) + DD )
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
        MOLAL (1) = KK                   ! HI
        MOLAL (5) = Y+KK                 ! SO4I
        MOLAL (6) = MAX (X+Y-KK, TINY)   ! HSO4I
        MOLAL (3) = 3.0*Y+X              ! NH4I
        CALL CALCMR                      ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF ( .NOT. CALAIN) GOTO 30
        CALL CALCACT
    20 END DO

    30 RETURN

! *** END OF SUBROUTINE CALCB3B ******************************************

    END SUBROUTINE CALCB3B
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : LC, (NH4)2SO4

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON THE SULFATE RATIO:
!     1. WHEN BOTH LC AND (NH4)2SO4 ARE POSSIBLE (SUBROUTINE CALCB2A)
!     2. WHEN ONLY LC IS POSSIBLE (SUBROUTINE CALCB2B).

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB2
    INCLUDE 'isrpia.inc'

! *** CALCULATE EQUIVALENT AMOUNT OF HSO4 AND SO4 ***********************

    X = MAX(2*W(2)-W(3), TINY)   ! Equivalent NH4HSO4
    Y = MAX(W(3)  -W(2), TINY)   ! Equivalent NH42SO4

! *** CALCULATE SPECIES ACCORDING TO RELATIVE ABUNDANCE OF HSO4 *********

    IF (X <= Y) THEN             ! LC is the MIN (x,y)
        SCASE = 'B2 ; SUBCASE 1'
        CALL CALCB2A (X,Y-X)      ! LC + (NH4)2SO4 POSSIBLE
    ELSE
        SCASE = 'B2 ; SUBCASE 2'
        CALL CALCB2B (Y,X-Y)      ! LC ONLY POSSIBLE
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCB2 ******************************************

    END SUBROUTINE CALCB2



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2 ; SUBCASE A.

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. SOLID PHASE ONLY POSSIBLE
!     3. SOLIDS POSSIBLE: LC, (NH4)2SO4

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE

!     FOR SOLID CALCULATIONS, A MATERIAL BALANCE BASED ON THE STOICHIMETRIC
!     PROPORTION OF AMMONIUM AND SULFATE IS DONE TO CALCULATE THE AMOUNT
!     OF LC AND (NH4)2SO4 IN THE SOLID PHASE.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB2A (TLC, TNH42S4)
    INCLUDE 'isrpia.inc'

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMLCAS) THEN
        SCASE   = 'B2 ; SUBCASE A1'    ! SOLIDS POSSIBLE ONLY
        CLC     = TLC
        CNH42S4 = TNH42S4
        SCASE   = 'B2 ; SUBCASE A1'
    ELSE
        SCASE = 'B2 ; SUBCASE A2'
        CALL CALCB2A2 (TLC, TNH42S4)   ! LIQUID & SOLID PHASE POSSIBLE
        SCASE = 'B2 ; SUBCASE A2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCB2A *****************************************

    END SUBROUTINE CALCB2A



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2A2
! *** CASE B2 ; SUBCASE A2.

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. SOLID PHASE ONLY POSSIBLE
!     3. SOLIDS POSSIBLE: LC, (NH4)2SO4

!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE CALCB2A1) AND THE
!     THE SOLID WITH LIQUID PHASE (SUBROUTINE CALCB3).

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB2A2 (TLC, TNH42S4)
    INCLUDE 'isrpia.inc'

! *** FIND WEIGHT FACTOR **********************************************

    IF (WFTYP == 0) THEN
        WF = ZERO
    ELSEIF (WFTYP == 1) THEN
        WF = 0.5D0
    ELSE
        WF = (DRLC-RH)/(DRLC-DRMLCAS)
    ENDIF
    ONEMWF  = ONE - WF

! *** FIND FIRST SECTION ; DRY ONE ************************************

    CLCO     = TLC                     ! FIRST (DRY) SOLUTION
    CNH42SO  = TNH42S4

! *** FIND SECOND SECTION ; DRY & LIQUID ******************************

    CLC     = ZERO
    CNH42S4 = ZERO
    CALL CALCB3                        ! SECOND (LIQUID) SOLUTION

! *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.

    MOLAL(1)= ONEMWF*MOLAL(1)                                   ! H+
    MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + 3.D0*(CLCO-CLC)) ! NH4+
    MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               ! SO4--
    MOLAL(6)= ONEMWF*(CLCO-CLC)                                 ! HSO4-

    WATER   = ONEMWF*WATER

    CLC     = WF*CLCO    + ONEMWF*CLC
    CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4

    RETURN

! *** END OF SUBROUTINE CALCB2A2 ****************************************

    END SUBROUTINE CALCB2A2



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2 ; SUBCASE B

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: LC

!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
!     AMOUNT OF SOLID LC DISSOLVED IN THE LIQUID PHASE.
!     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB2A CALCULATES THE
!     AMOUNT OF H+ PRODUCED (BASED ON THE HSO4, SO4 RELEASED INTO THE
!     SOLUTION). THE SOLUBILITY PRODUCT OF LC IS USED AS THE OBJECTIVE
!     FUNCTION.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB2B (TLC,TNH4HS4)
    INCLUDE 'isrpia.inc'

    CALAOU = .TRUE.       ! Outer loop activity calculation flag
    ZLO    = ZERO
    ZHI    = TLC          ! High limit: all of it in liquid phase

! *** INITIAL VALUES FOR BISECTION **************************************

    X1 = ZHI
    Y1 = FUNCB2B (X1,TNH4HS4,TLC)
    IF (ABS(Y1) <= EPS) RETURN
    YHI= Y1                        ! Save Y-value at Hi position

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ************************

    DX = (ZHI-ZLO)/NDIV
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCB2B (X2,TNH4HS4,TLC)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION FOUND

    YLO= Y1                      ! Save Y-value at LO position
    IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        RETURN
    
    ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
    
    ELSE IF (YLO < ZERO .AND. YHI < ZERO) THEN
        X1 = ZHI
        X2 = ZHI
        GOTO 40
    
    ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
    
    ELSE IF (YLO > ZERO .AND. YHI > ZERO) THEN
        X1 = ZLO
        X2 = ZLO
        GOTO 40
    ELSE
        CALL PUSHERR (0001, 'CALCB2B')    ! WARNING ERROR: NO SOLUTION
        RETURN
    ENDIF

! *** PERFORM BISECTION *************************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCB2B (X3,TNH4HS4,TLC)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCB2B')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN ************************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCB2B (X3,TNH4HS4,TLC)

    RETURN

! *** END OF SUBROUTINE CALCB2B *****************************************

    END SUBROUTINE CALCB2B



!=======================================================================

! *** ISORROPIA CODE
! *** FUNCTION FUNCB2B
! *** CASE B2 ;
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B2 ; SUBCASE 2
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCB2B.

!=======================================================================

    real FUNCTION FUNCB2B (X,TNH4HS4,TLC)
    INCLUDE 'isrpia.inc'

! *** SOLVE EQUATIONS **************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    DO 20 I=1,NSWEEP
        GRAT2 = XK1*WATER*(GAMA(8)/GAMA(7))**2./GAMA(7)
        PARM  = X+GRAT2
        DELTA = PARM*PARM + 4.0*(X+TNH4HS4)*GRAT2 ! Diakrinousa
        OMEGA = 0.5*(-PARM + SQRT(DELTA))         ! Thetiki riza (ie:H+>0)
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
        MOLAL (1) = OMEGA                         ! HI
        MOLAL (3) = 3.0*X+TNH4HS4                 ! NH4I
        MOLAL (5) = X+OMEGA                       ! SO4I
        MOLAL (6) = MAX (X+TNH4HS4-OMEGA, TINY)   ! HSO4I
        CLC       = MAX(TLC-X,ZERO)               ! Solid LC
        CALL CALCMR                               ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 30
        ENDIF
    20 END DO

! *** CALCULATE OBJECTIVE FUNCTION **************************************

! C30    FUNCB2B= ( NH4I**3.*SO4I*HSO4I )/( XK13*(WATER/GAMA(13))**5. )
    30 FUNCB2B= (MOLAL(3)**3.)*MOLAL(5)*MOLAL(6)
    FUNCB2B= FUNCB2B/(XK13*(WATER/GAMA(13))**5.) - ONE
    RETURN

! *** END OF FUNCTION FUNCB2B *******************************************

  END FUNCTION FUNCB2B


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1
! *** CASE B1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : LC, (NH4)2SO4, NH4HSO4

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCB1A)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB1
    INCLUDE 'isrpia.inc'

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMLCAB) THEN
        SCASE = 'B1 ; SUBCASE 1'
        CALL CALCB1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'B1 ; SUBCASE 1'
    ELSE
        SCASE = 'B1 ; SUBCASE 2'
        CALL CALCB1B              ! LIQUID & SOLID PHASE POSSIBLE
        SCASE = 'B1 ; SUBCASE 2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCB1 ******************************************

    END SUBROUTINE CALCB1



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1A
! *** CASE B1 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH
!     2. THERE IS NO LIQUID PHASE
!     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
!                         BUT NOT BOTH)

!     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE AMOUNT OF LC
!     IS CALCULATED FROM THE (NH4)2SO4 AND NH4HSO4 WHICH IS LEAST
!     ABUNDANT (STOICHIMETRICALLY). THE REMAINING EXCESS OF SALT
!     IS MIXED WITH THE LC.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB1A
    INCLUDE 'isrpia.inc'

! *** SETUP PARAMETERS ************************************************

    X = 2*W(2)-W(3)       ! Equivalent NH4HSO4
    Y = W(3)-W(2)         ! Equivalent (NH4)2SO4

! *** CALCULATE COMPOSITION *******************************************

    IF (X <= Y) THEN      ! LC is the MIN (x,y)
        CLC     = X        ! NH4HSO4 >= (NH4)2S04
        CNH4HS4 = ZERO
        CNH42S4 = Y-X
    ELSE
        CLC     = Y        ! NH4HSO4 <  (NH4)2S04
        CNH4HS4 = X-Y
        CNH42S4 = ZERO
    ENDIF
    RETURN

! *** END OF SUBROUTINE CALCB1 ******************************************

    END SUBROUTINE CALCB1A




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1B
! *** CASE B1 ; SUBCASE 2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
!                         BUT NOT BOTH)

!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE CALCB1A) AND THE
!     THE SOLID WITH LIQUID PHASE (SUBROUTINE CALCB2).

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCB1B
    INCLUDE 'isrpia.inc'

! *** FIND WEIGHT FACTOR **********************************************

    IF (WFTYP == 0) THEN
        WF = ZERO
    ELSEIF (WFTYP == 1) THEN
        WF = 0.5D0
    ELSE
        WF = (DRNH4HS4-RH)/(DRNH4HS4-DRMLCAB)
    ENDIF
    ONEMWF  = ONE - WF

! *** FIND FIRST SECTION ; DRY ONE ************************************

    CALL CALCB1A
    CLCO     = CLC               ! FIRST (DRY) SOLUTION
    CNH42SO  = CNH42S4
    CNH4HSO  = CNH4HS4

! *** FIND SECOND SECTION ; DRY & LIQUID ******************************

    CLC     = ZERO
    CNH42S4 = ZERO
    CNH4HS4 = ZERO
    CALL CALCB2                  ! SECOND (LIQUID) SOLUTION

! *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.

    MOLAL(1)= ONEMWF*MOLAL(1)                                   ! H+
    MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + (CNH4HSO-CNH4HS4) &
    + 3.D0*(CLCO-CLC))                          ! NH4+
    MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               ! SO4--
    MOLAL(6)= ONEMWF*(CNH4HSO-CNH4HS4 + CLCO-CLC)               ! HSO4-

    WATER   = ONEMWF*WATER

    CLC     = WF*CLCO    + ONEMWF*CLC
    CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
    CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4

    RETURN

! *** END OF SUBROUTINE CALCB1B *****************************************

    END SUBROUTINE CALCB1B


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCC2
! *** CASE C2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS ONLY A LIQUID PHASE

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCC2
    INCLUDE 'isrpia.inc'
    real :: LAMDA, KAPA

    CALAOU = .TRUE.         ! Outer loop activity calculation flag
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS **************************************************

    LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
    PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION
    DO 20 I=1,NSWEEP
        PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.
        BB    = PSI+PARM
        CC    =-PARM*(LAMDA+PSI)
        KAPA  = 0.5*(-BB+SQRT(BB*BB-4.0*CC))
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
        MOLAL(1) = PSI+KAPA                               ! HI
        MOLAL(3) = LAMDA                                  ! NH4I
        MOLAL(5) = KAPA                                   ! SO4I
        MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
        CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3), ZERO)  ! Free H2SO4
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF ( .NOT. CALAIN) GOTO 30
        CALL CALCACT
    20 END DO

    30 RETURN

! *** END OF SUBROUTINE CALCC2 *****************************************

    END SUBROUTINE CALCC2



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCC1
! *** CASE C1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE: NH4HSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCC1
    INCLUDE 'isrpia.inc'
    real :: KLO, KHI

    CALAOU = .TRUE.    ! Outer loop activity calculation flag
    KLO    = TINY
    KHI    = W(3)

! *** INITIAL VALUES FOR BISECTION *************************************

    X1 = KLO
    Y1 = FUNCC1 (X1)
    IF (ABS(Y1) <= EPS) GOTO 50
    YLO= Y1

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************

    DX = (KHI-KLO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCC1 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20 ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION FOUND

    YHI= Y2                 ! Save Y-value at HI position
    IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    
    ! *** { YLO, YHI } < 0.0  SOLUTION IS ALWAYS UNDERSATURATED WITH NH4HS04
    
    ELSE IF (YLO < ZERO .AND. YHI < ZERO) THEN
        GOTO 50
    
    ! *** { YLO, YHI } > 0.0 SOLUTION IS ALWAYS SUPERSATURATED WITH NH4HS04
    
    ELSE IF (YLO > ZERO .AND. YHI > ZERO) THEN
        X1 = KLO
        X2 = KLO
        GOTO 40
    ELSE
        CALL PUSHERR (0001, 'CALCC1')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION OF DISSOLVED NH4HSO4 **************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCC1 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCC1')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN ***********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCC1 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCC1 *****************************************

    END SUBROUTINE CALCC1



!=======================================================================

! *** ISORROPIA CODE
! *** FUNCTION FUNCC1
! *** CASE C1 ;
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE C1
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCC1.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCC1 (KAPA)
    INCLUDE 'isrpia.inc'
    real :: KAPA, LAMDA

! *** SOLVE EQUATIONS **************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

    PSI = W(2)-W(3)
    DO 20 I=1,NSWEEP
        PAR1  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
        PAR2  = XK12*(WATER/GAMA(9))**2.0
        BB    = PSI + PAR1
        CC    =-PAR1*(PSI+KAPA)
        LAMDA = 0.5*(-BB+SQRT(BB*BB-4*CC))
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY *******************************
    
        MOLAL(1) = PSI+LAMDA                    ! HI
        MOLAL(3) = KAPA                         ! NH4I
        MOLAL(5) = LAMDA                        ! SO4I
        MOLAL(6) = MAX (ZERO, PSI+KAPA-LAMDA)   ! HSO4I
        CNH4HS4  = MAX(W(3)-MOLAL(3), ZERO)     ! Solid NH4HSO4
        CH2SO4   = MAX(PSI, ZERO)               ! Free H2SO4
        CALL CALCMR                             ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 30
        ENDIF
    20 END DO

! *** CALCULATE ZERO FUNCTION *******************************************

! C30    FUNCC1= (NH4I*HSO4I/PAR2) - ONE
    30 FUNCC1= (MOLAL(3)*MOLAL(6)/PAR2) - ONE
    RETURN

! *** END OF FUNCTION FUNCC1 ********************************************

  END FUNCTION FUNCC1

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCD3
! *** CASE D3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. THERE IS OLNY A LIQUID PHASE

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCD3
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCD1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4NO3               ! Save from CALCD1 run
    CHI2 = CNH42S4
    CHI3 = GHNO3
    CHI4 = GNH3

    PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
    PSI2 = CHI2
    PSI3 = ZERO
    PSI4 = ZERO

    MOLAL(5) = ZERO
    MOLAL(6) = ZERO
    MOLAL(3) = PSI1
    MOLAL(7) = PSI1
    CALL CALCMR                  ! Initial water

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI4LO = TINY                ! Low  limit
    PSI4HI = CHI4                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    60 X1 = PSI4LO
    Y1 = FUNCD3 (X1)
    IF (ABS(Y1) <= EPS) RETURN
    YLO= Y1                 ! Save Y-value at HI position

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCD3 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION FOUND

    YHI= Y1                      ! Save Y-value at Hi position
    IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        RETURN
    
    ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
    ! Physically I dont know when this might happen, but I have put this
    ! branch in for completeness. I assume there is no solution; all NO3 goes to the
    ! gas phase.
    
    ELSE IF (YLO < ZERO .AND. YHI < ZERO) THEN
        P4 = TINY ! PSI4LO ! CHI4
        YY = FUNCD3(P4)
        GOTO 50
    
    ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
    ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
    ! and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
    ! and proceed again with root tracking.
    
    ELSE IF (YLO > ZERO .AND. YHI > ZERO) THEN
        PSI4HI = PSI4LO
        PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) ! No solution; some NH3 evaporates
        IF (PSI4LO < -(PSI1+PSI2)) THEN
            CALL PUSHERR (0001, 'CALCD3')  ! WARNING ERROR: NO SOLUTION
            RETURN
        ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
            GOTO 60                        ! Redo root tracking
        ENDIF
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCD3 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*ABS(X1)) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCD3')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCD3 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF
    RETURN

! *** END OF SUBROUTINE CALCD3 ******************************************

    END SUBROUTINE CALCD3



!=======================================================================

! *** ISORROPIA CODE
! *** FUNCTION FUNCD3
! *** CASE D3
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ;
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.

!=======================================================================

    real FUNCTION FUNCD3 (P4)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    PSI4   = P4

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A2   = XK7*(WATER/GAMA(4))**3.0
        A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.0
        A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A7   = XKW *RH*WATER*WATER
    
        PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
        PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4)
        PSI3 = MIN(MAX(PSI3, ZERO), CHI3)
    
        BB   = PSI4 - PSI3
    ! COLD         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
    ! C         AHI  =2.0*A7/(BB+SQRT(BB*BB + 4.d0*A7)) ! Avoid overflow when HI->0
        DENM = BB+SQRT(BB*BB + 4.d0*A7)
        IF (DENM <= TINY) THEN       ! Avoid overflow when HI->0
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.0*A7/ABB ! Taylor expansion of SQRT
        ENDIF
        AHI = 2.0*A7/DENM
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
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
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 CONTINUE
! C      FUNCD3= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE
    FUNCD3= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCD3 ********************************************

  END FUNCTION FUNCD3
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCD2
! *** CASE D2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCD2
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCD1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4NO3               ! Save from CALCD1 run
    CHI2 = CNH42S4
    CHI3 = GHNO3
    CHI4 = GNH3

    PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
    PSI2 = CNH42S4
    PSI3 = ZERO
    PSI4 = ZERO

    MOLAL(5) = ZERO
    MOLAL(6) = ZERO
    MOLAL(3) = PSI1
    MOLAL(7) = PSI1
    CALL CALCMR                  ! Initial water

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI4LO = TINY                ! Low  limit
    PSI4HI = CHI4                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    60 X1 = PSI4LO
    Y1 = FUNCD2 (X1)
    IF (ABS(Y1) <= EPS) RETURN
    YLO= Y1                 ! Save Y-value at HI position

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX   = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCD2 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) THEN
        
        ! This is done, in case if Y(PSI4LO)>0, but Y(PSI4LO+DX) < 0 (i.e.undersat)
        
            IF (Y1 <= Y2) GOTO 20  ! (Y1*Y2 < ZERO)
        ENDIF
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION FOUND

    YHI= Y1                      ! Save Y-value at Hi position
    IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        RETURN
    
    ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
    ! Physically I dont know when this might happen, but I have put this
    ! branch in for completeness. I assume there is no solution; all NO3 goes to the
    ! gas phase.
    
    ELSE IF (YLO < ZERO .AND. YHI < ZERO) THEN
        P4 = TINY ! PSI4LO ! CHI4
        YY = FUNCD2(P4)
        GOTO 50
    
    ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
    ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
    ! and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
    ! and proceed again with root tracking.
    
    ELSE IF (YLO > ZERO .AND. YHI > ZERO) THEN
        PSI4HI = PSI4LO
        PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) ! No solution; some NH3 evaporates
        IF (PSI4LO < -(PSI1+PSI2)) THEN
            CALL PUSHERR (0001, 'CALCD2')  ! WARNING ERROR: NO SOLUTION
            RETURN
        ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
            GOTO 60                        ! Redo root tracking
        ENDIF
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCD2 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*ABS(X1)) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCD2')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = MIN(X1,X2)   ! 0.5*(X1+X2)  ! Get "low" side, it's acidic soln.
    Y3 = FUNCD2 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF
    RETURN

! *** END OF SUBROUTINE CALCD2 ******************************************

    END SUBROUTINE CALCD2



!=======================================================================

! *** ISORROPIA CODE
! *** FUNCTION FUNCD2
! *** CASE D2
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D2 ;
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD2.

!=======================================================================

    real FUNCTION FUNCD2 (P4)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    CALL RSTGAM       ! Reset activity coefficients to 0.1
    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    PSI4   = P4
    PSI2   = CHI2

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
        A2  = XK7*(WATER/GAMA(4))**3.0
        A3  = XK4*R*TEMP*(WATER/GAMA(10))**2.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A7  = XKW *RH*WATER*WATER
    
        IF (CHI2 > TINY .AND. WATER > TINY) THEN
            PSI14 = PSI1+PSI4
            CALL POLY3 (PSI14,0.25*PSI14**2.,-A2/4.D0, PSI2, ISLV)  ! PSI2
            IF (ISLV == 0) THEN
                PSI2 = MIN (PSI2, CHI2)
            ELSE
                PSI2 = TINY
            ENDIF
        ENDIF
    
        PSI3  = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
        PSI3  = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4)
    ! c         PSI3  = MIN(MAX(PSI3, ZERO), CHI3)
    
        BB   = PSI4-PSI3 ! (BB > 0, acidic solution, <0 alkaline)
    
    ! Do not change computation scheme for H+, all others did not work well.
    
        DENM = BB+SQRT(BB*BB + 4.d0*A7)
        IF (DENM <= TINY) THEN       ! Avoid overflow when HI->0
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.d0*A7/ABB ! Taylor expansion of SQRT
        ENDIF
        AHI = 2.d0*A7/DENM
    
    ! *** SPECIATION & WATER CONTENT ***************************************
    
        MOLAL (1) = AHI                              ! HI
        MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2          ! NH4
        MOLAL (5) = PSI2                             ! SO4
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI3 + PSI1                      ! NO3
        CNH42S4   = CHI2 - PSI2                      ! Solid (NH4)2SO4
        CNH4NO3   = ZERO                             ! Solid NH4NO3
        GHNO3     = CHI3 - PSI3                      ! Gas HNO3
        GNH3      = CHI4 - PSI4                      ! Gas NH3
        CALL CALCMR                                  ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 CONTINUE
! C      FUNCD2= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE
    FUNCD2= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCD2 ********************************************

  END FUNCTION FUNCD2
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCD1
! *** CASE D1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3

!     THERE ARE TWO REGIMES DEFINED BY RELATIVE HUMIDITY:
!     1. RH < MDRH ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCD1A)
!     2. RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCD1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCD1A, CALCD2

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMASAN) THEN
        SCASE = 'D1 ; SUBCASE 1'   ! SOLID PHASE ONLY POSSIBLE
        CALL CALCD1A
        SCASE = 'D1 ; SUBCASE 1'
    ELSE
        SCASE = 'D1 ; SUBCASE 2'   ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRH (RH, DRMASAN, DRNH4NO3, CALCD1A, CALCD2)
        SCASE = 'D1 ; SUBCASE 2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCD1 ******************************************

    END SUBROUTINE CALCD1


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCD1A
! *** CASE D1 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3

!     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
!     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
!     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
!     THE SOLID PHASE.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCD1A
    INCLUDE 'isrpia.inc'

! *** SETUP PARAMETERS ************************************************

    PARM    = XK10/(R*TEMP)/(R*TEMP)

! *** CALCULATE NH4NO3 THAT VOLATIZES *********************************

    CNH42S4 = W(2)
    X       = MAX(ZERO, MIN(W(3)-2.0*CNH42S4, W(4)))  ! MAX NH4NO3
    PS      = MAX(W(3) - X - 2.0*CNH42S4, ZERO)
    OM      = MAX(W(4) - X, ZERO)

    OMPS    = OM+PS
    DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              ! DIAKRINOUSA
    ZE      = MIN(X, 0.5*(-OMPS + DIAK))              ! THETIKI RIZA

! *** SPECIATION *******************************************************

    CNH4NO3 = X  - ZE    ! Solid NH4NO3
    GNH3    = PS + ZE    ! Gas NH3
    GHNO3   = OM + ZE    ! Gas HNO3

    RETURN

! *** END OF SUBROUTINE CALCD1A *****************************************

    END SUBROUTINE CALCD1A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCG5
! *** CASE G5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCG5
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
    PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, &
    A1,   A2,   A3,   A4,   A5,   A6,   A7

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE. 
    CHI1   = 0.5*W(1)
    CHI2   = MAX (W(2)-CHI1, ZERO)
    CHI3   = ZERO
    CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
    CHI5   = W(4)
    CHI6   = W(5)

    PSI1   = CHI1
    PSI2   = CHI2
    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

    WATER  = CHI2/M0(4) + CHI1/M0(2)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCG5A (X1)
    IF (CHI6 <= TINY) GOTO 50
! c      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
! c      IF (WATER .LE. TINY) RETURN                    ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCG5A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCG5A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCG5A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCG5A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN  ! If quadrat.called
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
        MOLAL(6) = DELTA                               ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCG5 *******************************************

    END SUBROUTINE CALCG5




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCG5A
! *** CASE G5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCG5A (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
    PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, &
    A1,   A2,   A3,   A4,   A5,   A6,   A7

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A2  = XK7 *(WATER/GAMA(4))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        AKK = A4*A6
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        IF (CHI5 >= TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
        ELSE
            PSI5 = TINY
        ENDIF
    
    ! C      IF(CHI4.GT.TINY) THEN
        IF(W(2) > TINY) THEN       ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
            PSI4 =0.5d0*(-BB - SQRT(DD))
        ELSE
            PSI4 = TINY
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = 2.0D0*PSI1                          ! NAI
        MOLAL (3) = 2.0*PSI2 + PSI4                     ! NH4I
        MOLAL (4) = PSI6                                ! CLI
        MOLAL (5) = PSI2 + PSI1                         ! SO4I
        MOLAL (6) = ZERO
        MOLAL (7) = PSI5                                ! NO3I
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
        GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
        GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl
    
        CNH42S4   = ZERO                                ! Solid (NH4)2SO4
        CNH4NO3   = ZERO                                ! Solid NH4NO3
        CNH4CL    = ZERO                                ! Solid NH4Cl
    
        CALL CALCMR                                     ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCG5A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
! C         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCG5A *******************************************

  END FUNCTION FUNCG5A

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCG4
! *** CASE G4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCG4
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
    PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, &
    A1,   A2,   A3,   A4,   A5,   A6,   A7

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE. 
    CHI1   = 0.5*W(1)
    CHI2   = MAX (W(2)-CHI1, ZERO)
    CHI3   = ZERO
    CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
    CHI5   = W(4)
    CHI6   = W(5)

    PSI2   = CHI2
    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

    WATER  = CHI2/M0(4) + CHI1/M0(2)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCG4A (X1)
    IF (CHI6 <= TINY) GOTO 50
! C      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY .OR. WATER .LE. TINY) GOTO 50
! C      IF (WATER .LE. TINY) RETURN                    ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2  = X1+DX
        Y2  = FUNCG4A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1  = X2
        Y1  = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCG4A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCG4A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCG4')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCG4A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCG4 *******************************************

    END SUBROUTINE CALCG4




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCG4A
! *** CASE G4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCG4A (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA, NAI, NH4I, NO3I
    COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
    PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, &
    A1,   A2,   A3,   A4,   A5,   A6,   A7

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A2  = XK7 *(WATER/GAMA(4))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        IF (CHI5 >= TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
        ELSE
            PSI5 = TINY
        ENDIF
    
    ! C      IF(CHI4.GT.TINY) THEN
        IF(W(2) > TINY) THEN       ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO) ! Patch proposed by Uma shankar, 19/11/2001
            PSI4 =0.5d0*(-BB - SQRT(DD))
        ELSE
            PSI4 = TINY
        ENDIF
    
    !  CALCULATE CONCENTRATIONS
    
        NH4I = 2.0*PSI2 + PSI4
        CLI  = PSI6
        SO4I = PSI2 + PSI1
        NO3I = PSI5
        NAI  = 2.0D0*PSI1
    
        CALL CALCPH(2.d0*SO4I+NO3I+CLI-NAI-NH4I, HI, OHI)
    
    ! *** Na2SO4 DISSOLUTION
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN        ! PSI1
            CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1, CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ELSE
            PSI1 = ZERO
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (1) = HI
        MOLAL (2) = NAI
        MOLAL (3) = NH4I
        MOLAL (4) = CLI
        MOLAL (5) = SO4I
        MOLAL (6) = ZERO
        MOLAL (7) = NO3I
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNH4CL    = ZERO
        CNA2SO4   = MAX(CHI1-PSI1,ZERO)
    
    ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCG4A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
! C         FUNCG4A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCG4A *******************************************

  END FUNCTION FUNCG4A

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCG3
! *** CASE G3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. LIQUID & SOLID PHASE ARE BOTH POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCG3
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCG1A, CALCG4

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (W(4) > TINY .AND. W(5) > TINY) THEN ! NO3,CL EXIST, WATER POSSIBLE
        SCASE = 'G3 ; SUBCASE 1'
        CALL CALCG3A
        SCASE = 'G3 ; SUBCASE 1'
    ELSE                                      ! NO3, CL NON EXISTANT
        SCASE = 'G1 ; SUBCASE 1'
        CALL CALCG1A
        SCASE = 'G1 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMG3) THEN        ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCG1A
            SCASE = 'G3 ; SUBCASE 2'
            RETURN
        ELSE
            SCASE = 'G3 ; SUBCASE 3'  ! MDRH REGION (NA2SO4, NH42S4)
            CALL CALCMDRH (RH, DRMG3, DRNH42S4, CALCG1A, CALCG4)
            SCASE = 'G3 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCG3 ******************************************

    END SUBROUTINE CALCG3


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCG3A
! *** CASE G3 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCG3A
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
    PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, &
    A1,   A2,   A3,   A4,   A5,   A6,   A7

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE. 
    CHI1   = 0.5*W(1)
    CHI2   = MAX (W(2)-CHI1, ZERO)
    CHI3   = ZERO
    CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
    CHI5   = W(4)
    CHI6   = W(5)

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

    WATER  = TINY

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCG3A (X1)
    IF (CHI6 <= TINY) GOTO 50
! C      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY .OR. WATER .LE. TINY) GOTO 50
! C      IF (WATER .LE. TINY) RETURN                    ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2  = X1+DX
        Y2  = FUNCG3A (X2)
    
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1  = X2
        Y1  = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCG3A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCG3A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCG3A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCG3A (X3)

! *** FINAL CALCULATIONS *************************************************

    50 CONTINUE

! *** Na2SO4 DISSOLUTION

    IF (CHI1 > TINY .AND. WATER > TINY) THEN        ! PSI1
        CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
        IF (ISLV == 0) THEN
            PSI1 = MIN (PSI1, CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
    ELSE
        PSI1 = ZERO
    ENDIF
    MOLAL(2) = 2.0D0*PSI1               ! Na+  EFFECT
    MOLAL(5) = MOLAL(5) + PSI1          ! SO4  EFFECT
    CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   ! NA2SO4(s) depletion

! *** HSO4 equilibrium

    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCG3A ******************************************

    END SUBROUTINE CALCG3A




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCG3A
! *** CASE G3 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCG3A (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
    PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, &
    A1,   A2,   A3,   A4,   A5,   A6,   A7

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI2   = CHI2
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A2  = XK7 *(WATER/GAMA(4))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        IF (CHI5 >= TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
        ELSE
            PSI5 = TINY
        ENDIF
    
    ! C      IF(CHI4.GT.TINY) THEN
        IF(W(2) > TINY) THEN       ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO)  ! Patch proposed by Uma Shankar, 19/11/01
            PSI4 =0.5d0*(-BB - SQRT(DD))
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI2 > TINY .AND. WATER > TINY) THEN
            CALL POLY3 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
            IF (ISLV == 0) PSI2 = MIN (PSI20, CHI2)
        ENDIF
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        MOLAL (2) = ZERO                                ! Na
        MOLAL (3) = 2.0*PSI2 + PSI4                     ! NH4I
        MOLAL (4) = PSI6                                ! CLI
        MOLAL (5) = PSI2                                ! SO4I
        MOLAL (6) = ZERO                                ! HSO4
        MOLAL (7) = PSI5                                ! NO3I
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
        GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
        GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl
    
        CNH42S4   = CHI2 - PSI2                         ! Solid (NH4)2SO4
        CNH4NO3   = ZERO                                ! Solid NH4NO3
        CNH4CL    = ZERO                                ! Solid NH4Cl
    
        CALL CALCMR                                     ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCG3A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
! C         FUNCG3A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCG3A *******************************************

  END FUNCTION FUNCG3A

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCG2
! *** CASE G2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. LIQUID & SOLID PHASE ARE BOTH POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCG2
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCG1A, CALCG3A, CALCG4

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    IF (W(4) > TINY) THEN        ! NO3 EXISTS, WATER POSSIBLE
        SCASE = 'G2 ; SUBCASE 1'
        CALL CALCG2A
        SCASE = 'G2 ; SUBCASE 1'
    ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
        SCASE = 'G1 ; SUBCASE 1'
        CALL CALCG1A
        SCASE = 'G1 ; SUBCASE 1'
    ENDIF

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (WATER <= TINY) THEN
        IF (RH < DRMG2) THEN             ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCG1A
            SCASE = 'G2 ; SUBCASE 2'
        ELSE
            IF (W(5) > TINY) THEN
                SCASE = 'G2 ; SUBCASE 3'    ! MDRH (NH4CL, NA2SO4, NH42S4)
                CALL CALCMDRH (RH, DRMG2, DRNH4CL, CALCG1A, CALCG3A)
                SCASE = 'G2 ; SUBCASE 3'
            ENDIF
            IF (WATER <= TINY .AND. RH >= DRMG3) THEN
                SCASE = 'G2 ; SUBCASE 4'    ! MDRH (NA2SO4, NH42S4)
                CALL CALCMDRH (RH, DRMG3, DRNH42S4, CALCG1A, CALCG4)
                SCASE = 'G2 ; SUBCASE 4'
            ELSE
                WATER = TINY
                DO 20 I=1,NIONS
                    MOLAL(I) = ZERO
                20 END DO
                CALL CALCG1A
                SCASE = 'G2 ; SUBCASE 2'
            ENDIF
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCG2 ******************************************

    END SUBROUTINE CALCG2


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCG2A
! *** CASE G2 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCG2A
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
    PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, &
    A1,   A2,   A3,   A4,   A5,   A6,   A7

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE. 
    CHI1   = 0.5*W(1)
    CHI2   = MAX (W(2)-CHI1, ZERO)
    CHI3   = ZERO
    CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
    CHI5   = W(4)
    CHI6   = W(5)

    PSI6LO = TINY
    PSI6HI = CHI6-TINY

    WATER  = TINY

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCG2A (X1)
    IF (CHI6 <= TINY) GOTO 50
! C      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
! C      IF (WATER .LE. TINY) GOTO 50               ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCG2A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) WATER = TINY
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCG2A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCG2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    IF (X3 <= TINY2) THEN   ! PRACTICALLY NO NITRATES, SO DRY SOLUTION
        WATER = TINY
    ELSE
        Y3 = FUNCG2A (X3)
    ENDIF

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE

! *** Na2SO4 DISSOLUTION

    IF (CHI1 > TINY .AND. WATER > TINY) THEN        ! PSI1
        CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
        IF (ISLV == 0) THEN
            PSI1 = MIN (PSI1, CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
    ELSE
        PSI1 = ZERO
    ENDIF
    MOLAL(2) = 2.0D0*PSI1               ! Na+  EFFECT
    MOLAL(5) = MOLAL(5) + PSI1          ! SO4  EFFECT
    CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   ! NA2SO4(s) depletion

! *** HSO4 equilibrium

    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA     ! H+   AFFECT
        MOLAL(5) = MOLAL(5) - DELTA     ! SO4  AFFECT
        MOLAL(6) = DELTA                ! HSO4 AFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCG2A ******************************************

    END SUBROUTINE CALCG2A




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCG2A
! *** CASE G2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCG2A (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA, &
    PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, &
    A1,   A2,   A3,   A4,   A5,   A6,   A7

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI2   = CHI2
    PSI3   = ZERO
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A2  = XK7 *(WATER/GAMA(4))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
    
        DENO = MAX(CHI6-PSI6-PSI3, ZERO)
        PSI5 = CHI5/((A6/A5)*(DENO/PSI6) + ONE)
    
        PSI4 = MIN(PSI5+PSI6,CHI4)
    
        IF (CHI2 > TINY .AND. WATER > TINY) THEN
            CALL POLY3 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
            IF (ISLV == 0) PSI2 = MIN (PSI20, CHI2)
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (2) = ZERO                             ! NA
        MOLAL (3) = 2.0*PSI2 + PSI4                  ! NH4I
        MOLAL (4) = PSI6                             ! CLI
        MOLAL (5) = PSI2                             ! SO4I
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI5                             ! NO3I
    
    ! C      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = MAX(CHI2 - PSI2, ZERO)
        CNH4NO3   = ZERO
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
    ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 IF (CHI4 <= TINY) THEN
        FUNCG2A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
    ELSE
        FUNCG2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    ENDIF

    RETURN

! *** END OF FUNCTION FUNCG2A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCG1
! *** CASE G1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4CL, NA2SO4

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCG1A)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCG1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCG1A, CALCG2A

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMG1) THEN
        SCASE = 'G1 ; SUBCASE 1'
        CALL CALCG1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'G1 ; SUBCASE 1'
    ELSE
        SCASE = 'G1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRH (RH, DRMG1, DRNH4NO3, CALCG1A, CALCG2A)
        SCASE = 'G1 ; SUBCASE 2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCG1 ******************************************

    END SUBROUTINE CALCG1


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCG1A
! *** CASE G1 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3

!     SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
!     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
!     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
!     THE SOLID PHASE.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCG1A
    INCLUDE 'isrpia.inc'
    real :: LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2

! *** CALCULATE NON VOLATILE SOLIDS ***********************************

    CNA2SO4 = MIN (0.5*W(1), W(2))
    FRNA    = MAX(W(1) - 2.D0*CNA2SO4, ZERO)
    SO4FR   = MAX(W(2) - CNA2SO4, ZERO)
!      CNH42S4 = W(2) - CNA2SO4
    CNH42S4 = MAX (SO4FR , ZERO)                  ! CNH42S4

! *** CALCULATE VOLATILE SPECIES **************************************

    ALF     = W(3) - 2.0*CNH42S4
    BET     = W(5)
    GAM     = W(4)

    RTSQ    = R*TEMP*R*TEMP
    A1      = XK6/RTSQ
    A2      = XK10/RTSQ

    THETA1  = GAM - BET*(A2/A1)
    THETA2  = A2/A1

! QUADRATIC EQUATION SOLUTION

    BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
    CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
    DD      = BB*BB - 4.0D0*CC
    IF (DD < ZERO) GOTO 100   ! Solve each reaction seperately

! TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID

    SQDD    = SQRT(DD)
    KAPA1   = 0.5D0*(-BB+SQDD)
    KAPA2   = 0.5D0*(-BB-SQDD)
    LAMDA1  = THETA1 + THETA2*KAPA1
    LAMDA2  = THETA1 + THETA2*KAPA2

    IF (KAPA1 >= ZERO .AND. LAMDA1 >= ZERO) THEN
        IF (ALF-KAPA1-LAMDA1 >= ZERO .AND. &
        BET-KAPA1 >= ZERO .AND. GAM-LAMDA1 >= ZERO) THEN
            KAPA = KAPA1
            LAMDA= LAMDA1
            GOTO 200
        ENDIF
    ENDIF

    IF (KAPA2 >= ZERO .AND. LAMDA2 >= ZERO) THEN
        IF (ALF-KAPA2-LAMDA2 >= ZERO .AND. &
        BET-KAPA2 >= ZERO .AND. GAM-LAMDA2 >= ZERO) THEN
            KAPA = KAPA2
            LAMDA= LAMDA2
            GOTO 200
        ENDIF
    ENDIF

! SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA

    100 KAPA  = ZERO
    LAMDA = ZERO
    DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
    DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)

! NH4CL EQUILIBRIUM

    IF (DD1 >= ZERO) THEN
        SQDD1 = SQRT(DD1)
        KAPA1 = 0.5D0*(ALF+BET + SQDD1)
        KAPA2 = 0.5D0*(ALF+BET - SQDD1)
    
        IF (KAPA1 >= ZERO .AND. KAPA1 <= MIN(ALF,BET)) THEN
            KAPA = KAPA1
        ELSE IF (KAPA2 >= ZERO .AND. KAPA2 <= MIN(ALF,BET)) THEN
            KAPA = KAPA2
        ELSE
            KAPA = ZERO
        ENDIF
    ENDIF

! NH4NO3 EQUILIBRIUM

    IF (DD2 >= ZERO) THEN
        SQDD2 = SQRT(DD2)
        LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
        LAMDA2= 0.5D0*(ALF+GAM - SQDD2)
    
        IF (LAMDA1 >= ZERO .AND. LAMDA1 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1
        ELSE IF (LAMDA2 >= ZERO .AND. LAMDA2 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
        ELSE
            LAMDA = ZERO
        ENDIF
    ENDIF

! IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION

    IF (KAPA > ZERO .AND. LAMDA > ZERO) THEN
        IF (BET < LAMDA/THETA1) THEN
            KAPA = ZERO
        ELSE
            LAMDA= ZERO
        ENDIF
    ENDIF

! *** CALCULATE COMPOSITION OF VOLATILE SPECIES ***********************

    200 CONTINUE
    CNH4NO3 = LAMDA
    CNH4CL  = KAPA

    GNH3    = MAX(ALF - KAPA - LAMDA, ZERO)
    GHNO3   = MAX(GAM - LAMDA, ZERO)
    GHCL    = MAX(BET - KAPA, ZERO)

    RETURN

! *** END OF SUBROUTINE CALCG1A *****************************************

    END SUBROUTINE CALCG1A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCH6
! *** CASE H6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCH6
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCH6A (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCH6A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCH6A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCH6A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCH6A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCH6 ******************************************

    END SUBROUTINE CALCH6




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH6A
! *** CASE H6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCH6A (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
        A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MAX(PSI5, TINY)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1                           ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = MAX(CHI7 - PSI7, ZERO)
        CNANO3    = MAX(CHI8 - PSI8, ZERO)
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCH6A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCH6A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCH5
! *** CASE H5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCH5
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    IF (W(4) <= TINY .AND. W(5) <= TINY) THEN
        SCASE = 'H5'
        CALL CALCH1A
        SCASE = 'H5'
        RETURN
    ENDIF

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCH5A (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCH5A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCH5A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCH5A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCH5')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCH5A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCH5 ******************************************

    END SUBROUTINE CALCH5




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH5A
! *** CASE H5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NONE

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCH5A (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
        A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MAX(PSI5, TINY)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN     ! NA2SO4 DISSOLUTION
            AA = PSI7+PSI8
            BB = AA*AA
            CC =-A1/4.D0
            CALL POLY3 (AA, BB, CC, PSI1, ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1, CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                ! NAI
        MOLAL (3) = PSI4                                   ! NH4I
        MOLAL (4) = PSI6 + PSI7                            ! CLI
        MOLAL (5) = PSI2 + PSI1                            ! SO4I
        MOLAL (6) = ZERO
        MOLAL (7) = PSI5 + PSI8                            ! NO3I
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = MAX(CHI7 - PSI7, ZERO)
        CNANO3    = MAX(CHI8 - PSI8, ZERO)
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
    
        CALL CALCMR                               ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCH5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCH5A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCH4
! *** CASE H4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCH4
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    IF (W(4) <= TINY .AND. W(5) <= TINY) THEN
        SCASE = 'H4'
        CALL CALCH1A
        SCASE = 'H4'
        RETURN
    ENDIF

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCH4A (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCH4A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCH4A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCH4A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCH4')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCH4A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                      ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                      ! SO4  EFFECT
        MOLAL(6) = DELTA                                 ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCH4 ******************************************

    END SUBROUTINE CALCH4




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH4A
! *** CASE H4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCH4A (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
        A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MAX(PSI5, TINY)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN     ! NA2SO4 DISSOLUTION
            AA = PSI7+PSI8
            BB = AA*AA
            CC =-A1/4.D0
            CALL POLY3 (AA, BB, CC, PSI1, ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1, CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                ! NAI
        MOLAL (3) = PSI4                                   ! NH4I
        MOLAL (4) = PSI6 + PSI7                            ! CLI
        MOLAL (5) = PSI2 + PSI1                            ! SO4I
        MOLAL (6) = ZERO
        MOLAL (7) = PSI5 + PSI8                            ! NO3I
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = MAX(CHI7 - PSI7, ZERO)
        CNANO3    = MAX(CHI8 - PSI8, ZERO)
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        DELT = MIN(GNH3, GHCL)
        BB = -(GNH3+GHCL)
        CC = GNH3*GHCL-A3
        DD = BB*BB - 4.D0*CC
        PSI31 = 0.5D0*(-BB + SQRT(DD))
        PSI32 = 0.5D0*(-BB - SQRT(DD))
        IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
            PSI3 = PSI31
        ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
            PSI3 = PSI32
        ELSE
            PSI3 = ZERO
        ENDIF
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                           ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCH4A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCH4A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCH3
! *** CASE H3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCH3
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    IF (W(4) <= TINY) THEN        ! NO3 NOT EXIST, WATER NOT POSSIBLE
        SCASE = 'H3'
        CALL CALCH1A
        SCASE = 'H3'
        RETURN
    ENDIF

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCH3A (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCH3A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCH3A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCH3A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCH3')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCH3A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCH3 ******************************************

    END SUBROUTINE CALCH3




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH3A
! *** CASE H3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCH3A (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
        A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MAX(PSI5, TINY)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
            PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN     ! NA2SO4 DISSOLUTION
            AA = PSI7+PSI8
            BB = AA*AA
            CC =-A1/4.D0
            CALL POLY3 (AA, BB, CC, PSI1, ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1, CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1             ! NAI
        MOLAL (3) = PSI4                                ! NH4I
        MOLAL (4) = PSI6 + PSI7                         ! CLI
        MOLAL (5) = PSI2 + PSI1                         ! SO4I
        MOLAL (6) = ZERO
        MOLAL (7) = PSI5 + PSI8                         ! NO3I
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = MAX(CHI7 - PSI7, ZERO)
        CNANO3    = MAX(CHI8 - PSI8, ZERO)
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        DELT = MIN(GNH3, GHCL)
        BB = -(GNH3+GHCL)
        CC = GNH3*GHCL-A3
        DD = BB*BB - 4.D0*CC
        PSI31 = 0.5D0*(-BB + SQRT(DD))
        PSI32 = 0.5D0*(-BB - SQRT(DD))
        IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
            PSI3 = PSI31
        ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
            PSI3 = PSI32
        ELSE
            PSI3 = ZERO
        ENDIF
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                                 ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCH3A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCH3A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCH2
! *** CASE H2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : NH4Cl, NA2SO4, NANO3, NACL

!     THERE ARE THREE REGIMES IN THIS CASE:
!     1. NH4NO3(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCH2A)
!     2. NH4NO3(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!     3. NH4NO3(s) NOT POSSIBLE, AND RH >= MDRH. (MDRH REGION)

!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES H1A, H2B
!     RESPECTIVELY (BECAUSE MDRH POINTS COINCIDE).

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCH2
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCH1A, CALCH3

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    IF (W(4) > TINY) THEN        ! NO3 EXISTS, WATER POSSIBLE
        SCASE = 'H2 ; SUBCASE 1'
        CALL CALCH2A
        SCASE = 'H2 ; SUBCASE 1'
    ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
        SCASE = 'H2 ; SUBCASE 1'
        CALL CALCH1A
        SCASE = 'H2 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY .AND. RH < DRMH2) THEN      ! DRY AEROSOL
        SCASE = 'H2 ; SUBCASE 2'
    
    ELSEIF (WATER <= TINY .AND. RH >= DRMH2) THEN  ! MDRH OF H2
        SCASE = 'H2 ; SUBCASE 3'
        CALL CALCMDRH (RH, DRMH2, DRNANO3, CALCH1A, CALCH3)
        SCASE = 'H2 ; SUBCASE 3'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCH2 ******************************************

    END SUBROUTINE CALCH2




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCH2A
! *** CASE H2 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCH2A
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCH2A (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCH2A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (Y2 > EPS) Y2 = FUNCH2A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCH2A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCH2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCH2A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
        MOLAL(6) = DELTA                               ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCH2A ******************************************

    END SUBROUTINE CALCH2A




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH2A
! *** CASE H2 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCH2A (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
        A64 = (XK3*XK2/XKW)*(GAMA(10)/GAMA(5)/GAMA(11))**2.0
        A64 = A64*(R*TEMP*WATER)**2.0
        A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MAX(PSI5, TINY)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
            PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
        IF (CHI8 > TINY .AND. WATER > TINY) THEN     ! NANO3 DISSOLUTION
            DIAK = (PSI7-PSI5)**2.D0 + 4.D0*A8
            PSI8 = 0.5D0*( -(PSI7+PSI5) + SQRT(DIAK) )
            PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
        ENDIF
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN     ! NA2SO4 DISSOLUTION
            AA = PSI7+PSI8
            BB = AA*AA
            CC =-A1/4.D0
            CALL POLY3 (AA, BB, CC, PSI1, ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1, CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                 ! NAI
        MOLAL (3) = PSI4                                    ! NH4I
        MOLAL (4) = PSI6 + PSI7                             ! CLI
        MOLAL (5) = PSI2 + PSI1                             ! SO4I
        MOLAL (6) = ZERO                                    ! HSO4I
        MOLAL (7) = PSI5 + PSI8                             ! NO3I
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = MAX(CHI7 - PSI7, ZERO)
        CNANO3    = MAX(CHI8 - PSI8, ZERO)
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        DELT = MIN(GNH3, GHCL)
        BB = -(GNH3+GHCL)
        CC = GNH3*GHCL-A3
        DD = BB*BB - 4.D0*CC
        PSI31 = 0.5D0*(-BB + SQRT(DD))
        PSI32 = 0.5D0*(-BB - SQRT(DD))
        IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
            PSI3 = PSI31
        ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
            PSI3 = PSI32
        ELSE
            PSI3 = ZERO
        ENDIF
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                        ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCH2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A64 - ONE

    RETURN

! *** END OF FUNCTION FUNCH2A *******************************************

    END


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCH1
! *** CASE H1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCH1A)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCH1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCH1A, CALCH2A

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMH1) THEN
        SCASE = 'H1 ; SUBCASE 1'
        CALL CALCH1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'H1 ; SUBCASE 1'
    ELSE
        SCASE = 'H1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRH (RH, DRMH1, DRNH4NO3, CALCH1A, CALCH2A)
        SCASE = 'H1 ; SUBCASE 2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCH1 ******************************************

    END SUBROUTINE CALCH1


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCH1A
! *** CASE H1 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCH1A
    INCLUDE 'isrpia.inc'
    real :: LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR, &
    NO3FR

! *** CALCULATE NON VOLATILE SOLIDS ***********************************

    CNA2SO4 = W(2)
    CNH42S4 = ZERO
    NAFR    = MAX (W(1)-2*CNA2SO4, ZERO)
    CNANO3  = MIN (NAFR, W(4))
    NO3FR   = MAX (W(4)-CNANO3, ZERO)
    CNACL   = MIN (MAX(NAFR-CNANO3, ZERO), W(5))
    CLFR    = MAX (W(5)-CNACL, ZERO)

! *** CALCULATE VOLATILE SPECIES **************************************

    ALF     = W(3)                     ! FREE NH3
    BET     = CLFR                     ! FREE CL
    GAM     = NO3FR                    ! FREE NO3

    RTSQ    = R*TEMP*R*TEMP
    A1      = XK6/RTSQ
    A2      = XK10/RTSQ

    THETA1  = GAM - BET*(A2/A1)
    THETA2  = A2/A1

! QUADRATIC EQUATION SOLUTION

    BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
    CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
    DD      = BB*BB - 4.0D0*CC
    IF (DD < ZERO) GOTO 100   ! Solve each reaction seperately

! TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID

    SQDD    = SQRT(DD)
    KAPA1   = 0.5D0*(-BB+SQDD)
    KAPA2   = 0.5D0*(-BB-SQDD)
    LAMDA1  = THETA1 + THETA2*KAPA1
    LAMDA2  = THETA1 + THETA2*KAPA2

    IF (KAPA1 >= ZERO .AND. LAMDA1 >= ZERO) THEN
        IF (ALF-KAPA1-LAMDA1 >= ZERO .AND. &
        BET-KAPA1 >= ZERO .AND. GAM-LAMDA1 >= ZERO) THEN
            KAPA = KAPA1
            LAMDA= LAMDA1
            GOTO 200
        ENDIF
    ENDIF

    IF (KAPA2 >= ZERO .AND. LAMDA2 >= ZERO) THEN
        IF (ALF-KAPA2-LAMDA2 >= ZERO .AND. &
        BET-KAPA2 >= ZERO .AND. GAM-LAMDA2 >= ZERO) THEN
            KAPA = KAPA2
            LAMDA= LAMDA2
            GOTO 200
        ENDIF
    ENDIF

! SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA

    100 KAPA  = ZERO
    LAMDA = ZERO
    DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
    DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)

! NH4CL EQUILIBRIUM

    IF (DD1 >= ZERO) THEN
        SQDD1 = SQRT(DD1)
        KAPA1 = 0.5D0*(ALF+BET + SQDD1)
        KAPA2 = 0.5D0*(ALF+BET - SQDD1)
    
        IF (KAPA1 >= ZERO .AND. KAPA1 <= MIN(ALF,BET)) THEN
            KAPA = KAPA1
        ELSE IF (KAPA2 >= ZERO .AND. KAPA2 <= MIN(ALF,BET)) THEN
            KAPA = KAPA2
        ELSE
            KAPA = ZERO
        ENDIF
    ENDIF

! NH4NO3 EQUILIBRIUM

    IF (DD2 >= ZERO) THEN
        SQDD2 = SQRT(DD2)
        LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
        LAMDA2= 0.5D0*(ALF+GAM - SQDD2)
    
        IF (LAMDA1 >= ZERO .AND. LAMDA1 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1
        ELSE IF (LAMDA2 >= ZERO .AND. LAMDA2 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
        ELSE
            LAMDA = ZERO
        ENDIF
    ENDIF

! IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION

    IF (KAPA > ZERO .AND. LAMDA > ZERO) THEN
        IF (BET < LAMDA/THETA1) THEN
            KAPA = ZERO
        ELSE
            LAMDA= ZERO
        ENDIF
    ENDIF

! *** CALCULATE COMPOSITION OF VOLATILE SPECIES ***********************

    200 CONTINUE
    CNH4NO3 = LAMDA
    CNH4CL  = KAPA

    GNH3    = ALF - KAPA - LAMDA
    GHNO3   = GAM - LAMDA
    GHCL    = BET - KAPA

    RETURN

! *** END OF SUBROUTINE CALCH1A *****************************************

    END SUBROUTINE CALCH1A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI6
! *** CASE I6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI6
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCI1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCI1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = CNA2SO4
    PSI5 = CNH42S4

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
        CC   =-A6*(PSI2 + PSI3 + PSI1)
        DD   = BB*BB - 4.D0*CC
        PSI6 = 0.5D0*(-BB + SQRT(DD))
    
    ! *** CALCULATE SPECIATION ********************************************
    
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
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

    20 RETURN

! *** END OF SUBROUTINE CALCI6 *****************************************

    END SUBROUTINE CALCI6

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI5
! *** CASE I5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI5
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCI1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCI1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = ZERO
    PSI5 = CNH42S4

    CALAOU = .TRUE.               ! Outer loop activity calculation flag
    PSI4LO = ZERO                ! Low  limit
    PSI4HI = CHI4                ! High limit

! *** IF NA2SO4(S) =0, CALL FUNCI5B FOR Y4=0 ***************************

    IF (CHI4 <= TINY) THEN
        Y1 = FUNCI5A (ZERO)
        GOTO 50
    ENDIF

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI4HI
    Y1 = FUNCI5A (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 **

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCI5A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH4CL

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCI5A (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCI5')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCI5A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCI5')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCI5A (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCI5 *****************************************

    END SUBROUTINE CALCI5




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCI5A
! *** CASE I5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCI5A (P4)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4     ! PSI3 already assigned in FUNCI5A
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK5 *(WATER/GAMA(2))**3.0
        A5 = XK7 *(WATER/GAMA(4))**3.0
        A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
        CC   =-A6*(PSI2 + PSI3 + PSI1)
        DD   = BB*BB - 4.D0*CC
        PSI6 = 0.5D0*(-BB + SQRT(DD))
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (1) = PSI6                            ! HI
        MOLAL (2) = 2.D0*PSI4 + PSI3                ! NAI
        MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1    ! NH4I
        MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6       ! SO4I
        MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6       ! HSO4I
        CLC       = ZERO
        CNAHSO4   = ZERO
        CNA2SO4   = CHI4 - PSI4
        CNH42S4   = ZERO
        CNH4HS4   = ZERO
        CALL CALCMR                                 ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4     = XK5 *(WATER/GAMA(2))**3.0
    FUNCI5A= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCI5A ********************************************

    END
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI4
! *** CASE I4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI4
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCI1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCI1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = ZERO
    PSI5 = ZERO

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI4LO = ZERO                ! Low  limit
    PSI4HI = CHI4                ! High limit

! *** IF NA2SO4(S) =0, CALL FUNCI4B FOR Y4=0 ***************************

    IF (CHI4 <= TINY) THEN
        Y1 = FUNCI4A (ZERO)
        GOTO 50
    ENDIF

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI4HI
    Y1 = FUNCI4A (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 **

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCI4A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH4CL

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCI4A (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCI4')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCI4A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCI4')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCI4A (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCI4 *****************************************

    END SUBROUTINE CALCI4




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCI4A
! *** CASE I4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCI4A (P4)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4     ! PSI3 already assigned in FUNCI4A
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK5 *(WATER/GAMA(2))**3.0
        A5 = XK7 *(WATER/GAMA(4))**3.0
        A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
        A7 = SQRT(A4/A5)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
        CC   =-A6*(PSI2 + PSI3 + PSI1)
        DD   = BB*BB - 4.D0*CC
        PSI6 = 0.5D0*(-BB + SQRT(DD))
    
        PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7
        PSI5 = MAX (MIN (PSI5, CHI5), ZERO)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (1) = PSI6                            ! HI
        MOLAL (2) = 2.D0*PSI4 + PSI3                ! NAI
        MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1    ! NH4I
        MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6       ! SO4I
        MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6       ! HSO4I
        CLC       = ZERO
        CNAHSO4   = ZERO
        CNA2SO4   = CHI4 - PSI4
        CNH42S4   = CHI5 - PSI5
        CNH4HS4   = ZERO
        CALL CALCMR                                 ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4     = XK5 *(WATER/GAMA(2))**3.0
    FUNCI4A= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCI4A ********************************************

    END
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI3
! *** CASE I3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC

!     THERE ARE THREE REGIMES IN THIS CASE:
!     1.(NA,NH4)HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCI3A)
!     2.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!     3.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL

!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
!     RESPECTIVELY

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI3
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCI1A, CALCI4

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCI1A

! *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************

    IF (CNH4HS4 > TINY .OR. CNAHSO4 > TINY) THEN
        SCASE = 'I3 ; SUBCASE 1'
        CALL CALCI3A                     ! FULL SOLUTION
        SCASE = 'I3 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMI3) THEN         ! SOLID SOLUTION
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCI1A
            SCASE = 'I3 ; SUBCASE 2'
        
        ELSEIF (RH >= DRMI3) THEN     ! MDRH OF I3
            SCASE = 'I3 ; SUBCASE 3'
            CALL CALCMDRH (RH, DRMI3, DRLC, CALCI1A, CALCI4)
            SCASE = 'I3 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCI3 ******************************************

    END SUBROUTINE CALCI3



!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI3A
! *** CASE I3 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI3A
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCI1A         ! Needed when called from CALCMDRH

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCI1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = ZERO
    PSI3 = CNAHSO4
    PSI4 = ZERO
    PSI5 = ZERO

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI2LO = ZERO                ! Low  limit
    PSI2HI = CHI2                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI2HI
    Y1 = FUNCI3A (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC *********

    IF (YHI < EPS) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, PSI2LO)
        Y2 = FUNCI3A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC

    IF (Y2 > EPS) Y2 = FUNCI3A (ZERO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCI3A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCI3A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCI3A (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCI3A *****************************************

    END SUBROUTINE CALCI3A

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCI3A
! *** CASE I3 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCI3A (P2)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI2   = P2                  ! Save PSI2 in COMMON BLOCK
    PSI4LO = ZERO                ! Low  limit for PSI4
    PSI4HI = CHI4                ! High limit for PSI4

! *** IF NH3 =0, CALL FUNCI3B FOR Y4=0 ********************************

    IF (CHI4 <= TINY) THEN
        FUNCI3A = FUNCI3B (ZERO)
        GOTO 50
    ENDIF

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI4HI
    Y1 = FUNCI3B (X1)
    IF (ABS(Y1) <= EPS) GOTO 50
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 *****

    IF (YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, PSI4LO)
        Y2 = FUNCI3B (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NA2SO4

    IF (Y2 > EPS) Y2 = FUNCI3B (PSI4LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCI3B (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0004, 'FUNCI3A')    ! WARNING ERROR: NO CONVERGENCE

! *** INNER LOOP CONVERGED **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCI3B (X3)

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    50 A2      = XK13*(WATER/GAMA(13))**5.0
    FUNCI3A = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.D0/A2 - ONE
    RETURN

! *** END OF FUNCTION FUNCI3A *******************************************

    END



!=======================================================================

! *** ISORROPIA CODE
! *** FUNCTION FUNCI3B
! *** CASE I3 ; SUBCASE 2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC

!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/

!=======================================================================

    real FUNCTION FUNCI3B (P4)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK5*(WATER/GAMA(2))**3.0
        A5 = XK7*(WATER/GAMA(4))**3.0
        A6 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
        A7 = SQRT(A4/A5)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
        CC   =-A6*(PSI2 + PSI3 + PSI1)
        DD   = BB*BB - 4.D0*CC
        PSI6 = 0.5D0*(-BB + SQRT(DD))
    
        PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7
        PSI5 = MAX (MIN (PSI5, CHI5), ZERO)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = PSI6                                  ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                      ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1          ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6             ! SO4I
        MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 - PSI6, TINY)  ! HSO4I
        CLC      = MAX(CHI2 - PSI2, ZERO)
        CNAHSO4  = ZERO
        CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
        CNH42S4  = MAX(CHI5 - PSI5, ZERO)
        CNH4HS4  = ZERO
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4     = XK5 *(WATER/GAMA(2))**3.0
    FUNCI3B= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCI3B ********************************************

    END
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI2
! *** CASE I2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC

!     THERE ARE THREE REGIMES IN THIS CASE:
!     1. NH4HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCI2A)
!     2. NH4HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!     3. NH4HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL

!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
!     RESPECTIVELY

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI2
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCI1A, CALCI3A

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCI1A

! *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************

    IF (CNH4HS4 > TINY) THEN
        SCASE = 'I2 ; SUBCASE 1'
        CALL CALCI2A
        SCASE = 'I2 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMI2) THEN         ! SOLID SOLUTION ONLY
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCI1A
            SCASE = 'I2 ; SUBCASE 2'
        
        ELSEIF (RH >= DRMI2) THEN     ! MDRH OF I2
            SCASE = 'I2 ; SUBCASE 3'
            CALL CALCMDRH (RH, DRMI2, DRNAHSO4, CALCI1A, CALCI3A)
            SCASE = 'I2 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCI2 ******************************************

    END SUBROUTINE CALCI2


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI2A
! *** CASE I2 ; SUBCASE A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI2A
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCI1A    ! Needed when called from CALCMDRH

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCI1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = ZERO
    PSI3 = ZERO
    PSI4 = ZERO
    PSI5 = ZERO

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI2LO = ZERO                ! Low  limit
    PSI2HI = CHI2                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI2HI
    Y1 = FUNCI2A (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC *********

    IF (YHI < EPS) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, PSI2LO)
        Y2 = FUNCI2A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC

    IF (Y2 > EPS) Y2 = FUNCI2A (ZERO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCI2A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCI2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCI2A (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCI2A *****************************************

    END SUBROUTINE CALCI2A




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCI2A
! *** CASE I2 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCI2A (P2)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    PSI2   = P2                  ! Save PSI2 in COMMON BLOCK
    PSI3   = CHI3
    PSI4   = CHI4
    PSI5   = CHI5
    PSI6   = ZERO

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A3 = XK11*(WATER/GAMA(12))**2.0
        A4 = XK5 *(WATER/GAMA(2))**3.0
        A5 = XK7 *(WATER/GAMA(4))**3.0
        A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
        A7 = SQRT(A4/A5)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        IF (CHI5 > TINY .AND. WATER > TINY) THEN
            PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7
            PSI5 = MAX(MIN (PSI5, CHI5), TINY)
        ENDIF
    
        IF (CHI4 > TINY .AND. WATER > TINY) THEN
            AA   = PSI2+PSI5+PSI6+PSI3
            BB   = PSI3*AA
            CC   = 0.25D0*(PSI3*PSI3*(PSI2+PSI5+PSI6)-A4)
            CALL POLY3 (AA, BB, CC, PSI4, ISLV)
            IF (ISLV == 0) THEN
                PSI4 = MIN (PSI4, CHI4)
            ELSE
                PSI4 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI3 > TINY .AND. WATER > TINY) THEN
            AA   = 2.D0*PSI4 + PSI2 + PSI1 - PSI6
            BB   = 2.D0*PSI4*(PSI2 + PSI1 - PSI6) - A3
            CC   = ZERO
            CALL POLY3 (AA, BB, CC, PSI3, ISLV)
            IF (ISLV == 0) THEN
                PSI3 = MIN (PSI3, CHI3)
            ELSE
                PSI3 = ZERO
            ENDIF
        ENDIF
    
        BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
        CC   =-A6*(PSI2 + PSI3 + PSI1)
        DD   = BB*BB - 4.D0*CC
        PSI6 = 0.5D0*(-BB + SQRT(DD))
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (1) = PSI6                           ! HI
        MOLAL (2) = 2.D0*PSI4 + PSI3               ! NAI
        MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1   ! NH4I
        MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6      ! SO4I
        MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6      ! HSO4I
        CLC       = CHI2 - PSI2
        CNAHSO4   = CHI3 - PSI3
        CNA2SO4   = CHI4 - PSI4
        CNH42S4   = CHI5 - PSI5
        CNH4HS4   = ZERO
        CALL CALCMR                                ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 A2      = XK13*(WATER/GAMA(13))**5.0
    FUNCI2A = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.D0/A2 - ONE
    RETURN

! *** END OF FUNCTION FUNCI2A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI1
! *** CASE I1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCI1A)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCI1A, CALCI2A

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMI1) THEN
        SCASE = 'I1 ; SUBCASE 1'
        CALL CALCI1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'I1 ; SUBCASE 1'
    ELSE
        SCASE = 'I1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRH (RH, DRMI1, DRNH4HS4, CALCI1A, CALCI2A)
        SCASE = 'I1 ; SUBCASE 2'
    ENDIF

! *** AMMONIA IN GAS PHASE **********************************************

!      CALL CALCNH3

    RETURN

! *** END OF SUBROUTINE CALCI1 ******************************************

    END SUBROUTINE CALCI1


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCI1A
! *** CASE I1 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCI1A
    INCLUDE 'isrpia.inc'

! *** CALCULATE NON VOLATILE SOLIDS ***********************************

    CNA2SO4 = 0.5D0*W(1)
    CNH4HS4 = ZERO
    CNAHSO4 = ZERO
    CNH42S4 = ZERO
    FRSO4   = MAX(W(2)-CNA2SO4, ZERO)

    CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
    FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
    FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)

    IF (FRSO4 <= TINY) THEN
        CLC     = MAX(CLC - FRNH4, ZERO)
        CNH42S4 = 2.D0*FRNH4

    ELSEIF (FRNH4 <= TINY) THEN
        CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
        CLC     = MAX(CLC-FRSO4, ZERO)
        IF (CNA2SO4 > TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
        ENDIF
    ENDIF

! *** CALCULATE GAS SPECIES *********************************************

    GHNO3 = W(4)
    GHCL  = W(5)
    GNH3  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCI1A *****************************************

    END SUBROUTINE CALCI1A
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCJ3
! *** CASE J3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS ONLY A LIQUID PHASE

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCJ3
    INCLUDE 'isrpia.inc'

    real :: LAMDA, KAPA

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

    LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
    CHI1   = W(1)                           ! NA TOTAL as NaHSO4
    CHI2   = W(3)                           ! NH4 TOTAL as NH4HSO4
    PSI1   = CHI1
    PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        BB   = A3+LAMDA                        ! KAPA
        CC   =-A3*(LAMDA + PSI1 + PSI2)
        DD   = BB*BB-4.D0*CC
        KAPA = 0.5D0*(-BB+SQRT(DD))
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (1) = LAMDA + KAPA                 ! HI
        MOLAL (2) = PSI1                         ! NAI
        MOLAL (3) = PSI2                         ! NH4I
        MOLAL (4) = ZERO                         ! CLI
        MOLAL (5) = KAPA                         ! SO4I
        MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA   ! HSO4I
        MOLAL (7) = ZERO                         ! NO3I
    
        CNAHSO4   = ZERO
        CNH4HS4   = ZERO
    
        CALL CALCMR                              ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 50
        ENDIF
    10 END DO

    50 RETURN

! *** END OF SUBROUTINE CALCJ3 ******************************************

    END SUBROUTINE CALCJ3
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCJ2
! *** CASE J2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NAHSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCJ2
    INCLUDE 'isrpia.inc'

    real :: LAMDA, KAPA
    COMMON /CASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, &
    A1,   A2,   A3

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    CHI1   = W(1)                ! NA TOTAL
    CHI2   = W(3)                ! NH4 TOTAL
    PSI1LO = TINY                ! Low  limit
    PSI1HI = CHI1                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI1HI
    Y1 = FUNCJ2 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH42SO4 ****

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCJ2 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH42SO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCJ2 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCJ2')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCJ2 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCJ2')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCJ2 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCJ2 ******************************************

    END SUBROUTINE CALCJ2




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCJ2
! *** CASE J2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCJ2 (P1)
    INCLUDE 'isrpia.inc'

    real :: LAMDA, KAPA
    COMMON /CASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, &
    A1,   A2,   A3

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

    LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
    PSI1   = P1
    PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1 = XK11 *(WATER/GAMA(12))**2.0
        A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        BB   = A3+LAMDA                        ! KAPA
        CC   =-A3*(LAMDA + PSI1 + PSI2)
        DD   = BB*BB-4.D0*CC
        KAPA = 0.5D0*(-BB+SQRT(DD))
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (1) = LAMDA + KAPA                  ! HI
        MOLAL (2) = PSI1                          ! NAI
        MOLAL (3) = PSI2                          ! NH4I
        MOLAL (4) = ZERO                          ! CLI
        MOLAL (5) = KAPA                          ! SO4I
        MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA    ! HSO4I
        MOLAL (7) = ZERO                          ! NO3I
    
        CNAHSO4   = MAX(CHI1-PSI1,ZERO)
        CNH4HS4   = ZERO
    
        CALL CALCMR                               ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 FUNCJ2 = MOLAL(2)*MOLAL(6)/A1 - ONE

! *** END OF FUNCTION FUNCJ2 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCJ1
! *** CASE J1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    SUBROUTINE CALCJ1
    INCLUDE 'isrpia.inc'

    real :: LAMDA, KAPA
    COMMON /CASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, &
    A1,   A2,   A3

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE.               ! Outer loop activity calculation flag
    CHI1   = W(1)                ! Total NA initially as NaHSO4
    CHI2   = W(3)                ! Total NH4 initially as NH4HSO4

    PSI1LO = TINY                ! Low  limit
    PSI1HI = CHI1                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI1HI
    Y1 = FUNCJ1 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH42SO4 ****

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCJ1 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH42SO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCJ1 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCJ1')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCJ1 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCJ1')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCJ1 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCJ1 ******************************************

    END SUBROUTINE CALCJ1




!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCJ1
! *** CASE J1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS

!=======================================================================

    real FUNCTION FUNCJ1 (P1)
    INCLUDE 'isrpia.inc'
    real :: LAMDA, KAPA
    COMMON /CASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, &
    A1,   A2,   A3

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

    LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
    PSI1   = P1

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1 = XK11 *(WATER/GAMA(12))**2.0
        A2 = XK12 *(WATER/GAMA(09))**2.0
        A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
    
        PSI2 = 0.5*(-(LAMDA+PSI1) + SQRT((LAMDA+PSI1)**2.D0+4.D0*A2))  ! PSI2
        PSI2 = MIN (PSI2, CHI2)
    
        BB   = A3+LAMDA                        ! KAPA
        CC   =-A3*(LAMDA + PSI2 + PSI1)
        DD   = BB*BB-4.D0*CC
        KAPA = 0.5D0*(-BB+SQRT(DD))
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (1) = LAMDA + KAPA                  ! HI
        MOLAL (2) = PSI1                          ! NAI
        MOLAL (3) = PSI2                          ! NH4I
        MOLAL (4) = ZERO
        MOLAL (5) = KAPA                          ! SO4I
        MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA    ! HSO4I
        MOLAL (7) = ZERO
    
        CNAHSO4   = MAX(CHI1-PSI1,ZERO)
        CNH4HS4   = MAX(CHI2-PSI2,ZERO)
    
        CALL CALCMR                               ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 FUNCJ1 = MOLAL(2)*MOLAL(6)/A1 - ONE

! *** END OF FUNCTION FUNCJ1 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO7
! *** CASE O7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MgSO4, NA2SO4, K2SO4
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO7
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

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

    WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
    WATER  = MAX (WATER , TINY)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCO7 (X1)
    IF (CHI6 <= TINY) GOTO 50
! c      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
! c      IF (WATER .LE. TINY) RETURN                    ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCO7 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCO7 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCO7 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCO7')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCO7 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN  ! If quadrat.called
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
        MOLAL(6) = DELTA                               ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO7 *******************************************

    END SUBROUTINE CALCO7

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCO7
! *** CASE O7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MgSO4, NA2SO4, K2SO4
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCO7 (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
    
    
        IF (CHI5 >= TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
            PSI5 = MIN (PSI5,CHI5)
        ELSE
            PSI5 = TINY
        ENDIF
    
    ! C      IF(CHI4.GT.TINY) THEN
        IF(W(2) > TINY) THEN       ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
        ELSE
            PSI4 = TINY
        ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (2) = 2.0D0*PSI1                       ! Na+
        MOLAL (3) = 2.0D0*PSI2 + PSI4                ! NH4I
        MOLAL (4) = PSI6                             ! CLI
        MOLAL (5) = PSI1+PSI2+PSI7+PSI8              ! SO4I
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI5                             ! NO3I
        MOLAL (8) = ZERO                             ! CaI
        MOLAL (9) = 2.0D0*PSI7                       ! KI
        MOLAL (10)= PSI8                             ! Mg
    
    ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
    
    ! C      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        -MOLAL(9)-2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNA2SO4  = ZERO
        CNH42S4  = ZERO
        CNH4NO3  = ZERO
        CNH4Cl   = ZERO
        CK2SO4   = ZERO
        CMGSO4   = ZERO
        CCASO4   = CHI9
    
    ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCO7 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
! C         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCO7 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO6
! *** CASE O6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MGSO4, NA2SO4
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO6
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

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


    PSI1   = CHI1
    PSI2   = CHI2
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = CHI8
    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

    WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
    WATER  = MAX (WATER , TINY)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCO6 (X1)
    IF (CHI6 <= TINY) GOTO 50
! c      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
! c      IF (WATER .LE. TINY) RETURN                    ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCO6 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCO6 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCO6 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCO6')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCO6 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN  ! If quadrat.called
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
        MOLAL(6) = DELTA                               ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO6 *******************************************

    END SUBROUTINE CALCO6

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCO6
! *** CASE O6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4 , K2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MgSO4, NA2SO4
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCO6 (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK17 *(WATER/GAMA(17))**3.0
    
    
        IF (CHI5 >= TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
            PSI5 = MIN (PSI5,CHI5)
        ELSE
            PSI5 = TINY
        ENDIF
    
    ! C      IF(CHI4.GT.TINY) THEN
        IF(W(2) > TINY) THEN       ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN        ! PSI7
            CALL POLY3 (PSI1+PSI2+PSI8, ZERO, -A7/4.D0, PSI7, ISLV)
            IF (ISLV == 0) THEN
                PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
            ELSE
                PSI7 = ZERO
            ENDIF
        ELSE
            PSI7 = ZERO
        ENDIF
    
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (2) = 2.0D0*PSI1                       ! Na+
        MOLAL (3) = 2.0D0*PSI2 + PSI4                ! NH4I
        MOLAL (4) = PSI6                             ! CLI
        MOLAL (5) = PSI1+PSI2+PSI7+PSI8              ! SO4I
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI5                             ! NO3I
        MOLAL (8) = ZERO                             ! CaI
        MOLAL (9) = 2.0D0*PSI7                       ! KI
        MOLAL (10)= PSI8                             ! Mg

    
    ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    
    ! C      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        -MOLAL(9)-2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3     = MAX(CHI4 - PSI4, TINY)
        GHNO3    = MAX(CHI5 - PSI5, TINY)
        GHCL     = MAX(CHI6 - PSI6, TINY)
    
        CNA2SO4  = ZERO
        CNH42S4  = ZERO
        CNH4NO3  = ZERO
        CNH4Cl   = ZERO
        CK2SO4   = MAX(CHI7 - PSI7, TINY)
        CMGSO4   = ZERO
        CCASO4   = CHI9
    
    ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCO6 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
! C         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCO6 *******************************************

    END


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO5
! *** CASE O5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MGSO4
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO5
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

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

    PSI1   = ZERO
    PSI2   = CHI2
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = CHI8
    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

    WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
    WATER  = MAX (WATER , TINY)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCO5 (X1)
    IF (CHI6 <= TINY) GOTO 50
! c      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
! c      IF (WATER .LE. TINY) RETURN                    ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCO5 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCO5 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCO5 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCO5')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCO5 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN  ! If quadrat.called
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
        MOLAL(6) = DELTA                               ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO5 *******************************************

    END SUBROUTINE CALCO5

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCO5
! *** CASE O5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MGSO4
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCO5 (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK17 *(WATER/GAMA(17))**3.0
    
    
        IF (CHI5 >= TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
            PSI5 = MIN (PSI5,CHI5)
        ELSE
            PSI5 = TINY
        ENDIF
    
    ! C      IF(CHI4.GT.TINY) THEN
        IF(W(2) > TINY) THEN       ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN        ! PSI7
            CALL POLY3 ((PSI2+PSI8)/(SQRT(A1/A7)+1.0), ZERO, &
            -A7/4.D0/(SQRT(A1/A7)+1.0), PSI7, ISLV)
            IF (ISLV == 0) THEN
                PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
            ELSE
                PSI7 = ZERO
            ENDIF
        ELSE
            PSI7 = ZERO
        ENDIF
    
        IF (CHI1 >= TINY) THEN                              ! PSI1
            PSI1   = SQRT(A1/A7)*PSI7
            PSI1   = MIN(PSI1,CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
    
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (2) = 2.0D0*PSI1                       ! NaI
        MOLAL (3) = 2.0D0*PSI2 + PSI4                ! NH4I
        MOLAL (4) = PSI6                             ! CLI
        MOLAL (5) = PSI1+PSI2+PSI7+PSI8              ! SO4I
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI5                             ! NO3I
        MOLAL (8) = ZERO                             ! CaI
        MOLAL (9) = 2.0D0*PSI7                       ! KI
        MOLAL (10)= PSI8                             ! Mg

    
    ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    
    ! C      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        -MOLAL(9)-2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3     = MAX(CHI4 - PSI4, TINY)
        GHNO3    = MAX(CHI5 - PSI5, TINY)
        GHCL     = MAX(CHI6 - PSI6, TINY)
    
        CNA2SO4  = MAX(CHI1 - PSI1, TINY)
        CNH42S4  = ZERO
        CNH4NO3  = ZERO
        CNH4Cl   = ZERO
        CK2SO4   = MAX(CHI7 - PSI7, TINY)
        CMGSO4   = ZERO
        CCASO4   = CHI9
    
    ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCO5 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
! C         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCO5 *******************************************

    END


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO4
! *** CASE O4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NA2SO4, K2SO4, MGSO4, CASO4
!     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO4
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

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

    PSI2   = CHI2
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = CHI8
    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

    WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
    WATER  = MAX (WATER , TINY)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCO4 (X1)
    IF (CHI6 <= TINY) GOTO 50
! C      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
! C      IF (WATER .LE. TINY) GOTO 50               ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCO4 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCO4 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCO4 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCO4')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCO4 (X3)

! *** FINAL CALCULATIONS **********************************************

    50 CONTINUE

! *** Na2SO4 DISSOLUTION

    IF (CHI1 > TINY .AND. WATER > TINY) THEN        ! PSI1
        CALL POLY3 (PSI2+PSI7+PSI8, ZERO, -A1/4.D0, PSI1, ISLV)
        IF (ISLV == 0) THEN
            PSI1 = MIN (PSI1, CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
    ELSE
        PSI1 = ZERO
    ENDIF
    MOLAL(2) = 2.0D0*PSI1               ! Na+  EFFECT
    MOLAL(5) = MOLAL(5) + PSI1          ! SO4  EFFECT
    CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   ! NA2SO4(s) depletion

! *** HSO4 equilibrium


    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA     ! H+   AFFECT
        MOLAL(5) = MOLAL(5) - DELTA     ! SO4  AFFECT
        MOLAL(6) = DELTA                ! HSO4 AFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO4 ******************************************

    END SUBROUTINE CALCO4


!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCO4
! *** CASE O4 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCO4 (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI2   = CHI2
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK17 *(WATER/GAMA(17))**3.0
    !      A8  = XK23 *(WATER/GAMA(21))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        IF (CHI5 >= TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
            PSI5 = MIN (PSI5,CHI5)
        ELSE
            PSI5 = TINY
        ENDIF
    
    ! C      IF(CHI4.GT.TINY) THEN
        IF(W(2) > TINY) THEN       ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN        ! PSI7
            CALL POLY3 (PSI2+PSI8, ZERO, -A7/4.D0, PSI7, ISLV)
            IF (ISLV == 0) THEN
                PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
            ELSE
                PSI7 = ZERO
            ENDIF
        ELSE
            PSI7 = ZERO
        ENDIF
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        MOLAL (2) = ZERO                             ! NAI
        MOLAL (3) = 2.0D0*PSI2 + PSI4                ! NH4I
        MOLAL (4) = PSI6                             ! CLI
        MOLAL (5) = PSI2+PSI7+PSI8                   ! SO4I
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI5                             ! NO3I
        MOLAL (8) = ZERO                             ! CAI
        MOLAL (9) = 2.0D0*PSI7                       ! KI
        MOLAL (10)= PSI8                             ! MGI

    
    ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    
    ! C      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        -MOLAL(9)-2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3     = MAX(CHI4 - PSI4, TINY)
        GHNO3    = MAX(CHI5 - PSI5, TINY)
        GHCL     = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4  = ZERO
        CNH4NO3  = ZERO
        CNH4Cl   = ZERO
        CK2SO4   = MAX(CHI7 - PSI7, TINY)
        CMGSO4   = ZERO
        CCASO4   = CHI9
    
        CALL CALCMR                                     ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCO4 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
! C         FUNCO4 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCO4 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO3
! *** CASE O3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO3
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCO1A, CALCO4

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (W(4) > TINY .AND. W(5) > TINY) THEN ! NO3,CL EXIST, WATER POSSIBLE
        SCASE = 'O3 ; SUBCASE 1'
        CALL CALCO3A
        SCASE = 'O3 ; SUBCASE 1'
    ELSE                                      ! NO3, CL NON EXISTANT
        SCASE = 'O1 ; SUBCASE 1'
        CALL CALCO1A
        SCASE = 'O1 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMO3) THEN        ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCO1A
            SCASE = 'O3 ; SUBCASE 2'
            RETURN
        ELSE
            SCASE = 'O3 ; SUBCASE 3'  ! MDRH REGION (NA2SO4, NH42S4, K2SO4, MGSO4, CASO4)
            CALL CALCMDRH2 (RH, DRMO3, DRNH42S4, CALCO1A, CALCO4)
            SCASE = 'O3 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO3 ******************************************

    END SUBROUTINE CALCO3

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO3A
! *** CASE O3 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, K2SO4, MGSO4, CASO4
!     4. Completely dissolved: NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO3A
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9, A1,  A2,  A3,  A4, &
    A5,  A6,  A7,  A8,  A9

! *** SETUP PARAMETERS ************************************************

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

    PSI8   = CHI8
    PSI6LO = TINY
    PSI6HI = CHI6-TINY

    WATER  = TINY

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCO3A (X1)
    IF (CHI6 <= TINY) GOTO 50
! C      IF (ABS(Y1).LE.EPS .OR. CHI7.LE.TINY) GOTO 50
! C      IF (WATER .LE. TINY) GOTO 50               ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCO3A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCO3A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCO3A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCO3A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCO3A (X3)

! *** FINAL CALCULATIONS *************************************************

    50 CONTINUE

! *** Na2SO4 DISSOLUTION

    IF (CHI1 > TINY .AND. WATER > TINY) THEN        ! PSI1
        CALL POLY3 (PSI2+PSI7+PSI8, ZERO, -A1/4.D0, PSI1, ISLV)
        IF (ISLV == 0) THEN
            PSI1 = MIN (max (PSI1, zero), CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
    ELSE
        PSI1 = ZERO
    ENDIF
    MOLAL(2) = 2.0D0*PSI1               ! Na+  EFFECT
    MOLAL(5) = MOLAL(5) + PSI1          ! SO4  EFFECT
    CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   ! NA2SO4(s) depletion

! *** HSO4 equilibrium

    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA     ! H+   AFFECT
        MOLAL(5) = MOLAL(5) - DELTA     ! SO4  AFFECT
        MOLAL(6) = DELTA                ! HSO4 AFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO3A ******************************************

    END SUBROUTINE CALCO3A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCO3A
! *** CASE O3; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, K2SO4, MgSO4, CaSO4
!     4. Completely dissolved: NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCO3A (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

    PSI2   = CHI2
    PSI8   = CHI8
    PSI3   = ZERO
    PSI6   = X

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0D0
        A2  = XK7 *(WATER/GAMA(4))**3.0D0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0D0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0D0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0D0
        A7  = XK17 *(WATER/GAMA(17))**3.0D0
    !      A8  = XK23 *(WATER/GAMA(21))**2.0D0
        A65 = A6/A5
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        DENO = MAX(CHI6-PSI6-PSI3, ZERO)
        PSI5 = PSI6*CHI5/(A6/A5*DENO + PSI6)
        PSI5 = MIN(MAX(PSI5,ZERO),CHI5)
    
    ! C      IF(CHI4.GT.TINY) THEN                             ! PSI4
        IF(W(2) > TINY) THEN       ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO)  ! Patch proposed by Uma Shankar, 19/11/01
            PSI4 =0.5d0*(-BB - SQRT(DD))
        ELSE
            PSI4 = TINY
        ENDIF
        PSI4 = MIN (MAX (PSI4,ZERO), CHI4)
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN        ! PSI7
            CALL POLY3 (PSI2+PSI8, ZERO, -A7/4.D0, PSI7, ISLV)
            IF (ISLV == 0) THEN
                PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
            ELSE
                PSI7 = ZERO
            ENDIF
        ELSE
            PSI7 = ZERO
        ENDIF
    
        IF (CHI2 > TINY .AND. WATER > TINY) THEN
            CALL POLY3 (PSI7+PSI8+PSI4, PSI4*(PSI7+PSI8)+ &
            PSI4*PSI4/4.D0, (PSI4*PSI4*(PSI7+PSI8)-A2) &
            /4.D0,PSI20, ISLV)
            IF (ISLV == 0) PSI2 = MIN (MAX(PSI20,ZERO), CHI2)
        ENDIF
    !      PSI2 = 0.5D0*(2.0D0*SQRT(A2/A7)*PSI7 - PSI4)
    !      PSI2 = MIN (MAX(PSI2, ZERO), CHI2)
    !      ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (2) = ZERO                             ! NaI
        MOLAL (3) = 2.0D0*PSI2 + PSI4                ! NH4I
        MOLAL (4) = PSI6                             ! CLI
        MOLAL (5) = PSI2+PSI7+PSI8                   ! SO4I
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI5                             ! NO3I
        MOLAL (8) = ZERO                             ! CAI
        MOLAL (9) = 2.0D0*PSI7                       ! KI
        MOLAL (10)= PSI8                             ! MGI
    
    ! C      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        -MOLAL(9)-2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
    !      CNA2SO4  = MAX(CHI1 - PSI1, ZERO)
        CNH42S4  = MAX(CHI2 - PSI2, ZERO)
        CNH4NO3  = ZERO
        CNH4Cl   = ZERO
        CK2SO4   = MAX(CHI7 - PSI7, ZERO)
        CMGSO4   = ZERO
        CCASO4   = CHI9
    
    ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    20 FUNCO3A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE


    RETURN

! *** END OF FUNCTION FUNCO3A *******************************************

    END


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO2
! *** CASE O2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO2
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCO1A, CALCO3A, CALCO4

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    IF (W(4) > TINY) THEN        ! NO3 EXISTS, WATER POSSIBLE
        SCASE = 'O2 ; SUBCASE 1'
        CALL CALCO2A
        SCASE = 'O2 ; SUBCASE 1'
    ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
        SCASE = 'O1 ; SUBCASE 1'
        CALL CALCO1A
        SCASE = 'O1 ; SUBCASE 1'
    ENDIF

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (WATER <= TINY) THEN
        IF (RH < DRMO2) THEN             ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCO1A
            SCASE = 'O2 ; SUBCASE 2'
        ELSE
            IF (W(5) > TINY) THEN
                SCASE = 'O2 ; SUBCASE 3'    ! MDRH (NH4CL, NA2SO4, NH42S4, K2SO4, MGSO4, CASO4)
                CALL CALCMDRH2 (RH, DRMO2, DRNH4CL, CALCO1A, CALCO3A)
                SCASE = 'O2 ; SUBCASE 3'
            ENDIF
            IF (WATER <= TINY .AND. RH >= DRMO3) THEN
                SCASE = 'O2 ; SUBCASE 4'    ! MDRH (NA2SO4, NH42S4, K2SO4, MGSO4, CASO4)
                CALL CALCMDRH2 (RH, DRMO3, DRNH42S4, CALCO1A, CALCO4)
                SCASE = 'O2 ; SUBCASE 4'
            ELSE
                WATER = TINY
                DO 20 I=1,NIONS
                    MOLAL(I) = ZERO
                20 END DO
                CALL CALCO1A
                SCASE = 'O2 ; SUBCASE 2'
            ENDIF
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO2 ******************************************

    END SUBROUTINE CALCO2

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO2A
! *** CASE O2 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4, K2SO4, MgSO4, CaSO4
!     4. Completely dissolved: NH4NO3

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO2A
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS *************************************************

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

    PSI8   = CHI8
    PSI6LO = TINY
    PSI6HI = CHI6-TINY

    WATER  = TINY

! *** INITIAL VALUES FOR BISECTION *************************************

    X1 = PSI6LO
    Y1 = FUNCO2A (X1)
    IF (CHI6 <= TINY) GOTO 50
! C      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
! C      IF (WATER .LE. TINY) GOTO 50               ! No water

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCO2A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) WATER = TINY
    GOTO 50

! *** PERFORM BISECTION ************************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCO2A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCO2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN ***********************************************

    40 X3 = 0.5*(X1+X2)
    IF (X3 <= TINY2) THEN   ! PRACTICALLY NO NITRATES, SO DRY SOLUTION
        WATER = TINY
    ELSE
        Y3 = FUNCO2A (X3)
    ENDIF

! *** FINAL CALCULATIONS *************************************************

    50 CONTINUE

! *** Na2SO4 DISSOLUTION

    IF (CHI1 > TINY .AND. WATER > TINY) THEN        ! PSI1
        CALL POLY3 (PSI2+PSI7+PSI8, ZERO, -A1/4.D0, PSI1, ISLV)
        IF (ISLV == 0) THEN
            PSI1 = MIN (PSI1, CHI1)
        ELSE
            PSI1 = ZERO
        ENDIF
    ELSE
        PSI1 = ZERO
    ENDIF
    MOLAL(2) = 2.0D0*PSI1               ! Na+  EFFECT
    MOLAL(5) = MOLAL(5) + PSI1          ! SO4  EFFECT
    CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   ! NA2SO4(s) depletion

! *** HSO4 equilibrium

    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA     ! H+   AFFECT
        MOLAL(5) = MOLAL(5) - DELTA     ! SO4  AFFECT
        MOLAL(6) = DELTA                ! HSO4 AFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO2A ******************************************

    END SUBROUTINE CALCO2A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCO2A
! *** CASE O2; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4, K2SO4, MgSO4, CaSO4
!     4. Completely dissolved: NH4NO3
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCO2A (X)
    INCLUDE 'isrpia.inc'

    real :: LAMDA
    COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5, &
    PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4, &
    A5, A6, A7, A8, A9

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI2   = CHI2
    PSI3   = ZERO

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0D0
        A2  = XK7 *(WATER/GAMA(4))**3.0D0
        A3  = XK6 /(R*TEMP*R*TEMP)
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0D0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0D0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0D0
        A65 = A6/A5
        A7  = XK17 *(WATER/GAMA(17))**3.0D0
    !      A8  = XK23 *(WATER/GAMA(21))**2.0D0
    
        DENO = MAX(CHI6-PSI6-PSI3, ZERO)
        PSI5 = PSI6*CHI5/(A6/A5*DENO + PSI6)
        PSI5 = MIN(PSI5,CHI5)
    
        PSI4 = MIN(PSI5+PSI6,CHI4)
    
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN        ! PSI7
            CALL POLY3 (PSI2+PSI8, ZERO, -A7/4.D0, PSI7, ISLV)
            IF (ISLV == 0) THEN
                PSI7 = MAX (MIN (PSI7, CHI7), ZERO)
            ELSE
                PSI7 = ZERO
            ENDIF
        ELSE
            PSI7 = ZERO
        ENDIF
    
        IF (CHI2 > TINY .AND. WATER > TINY) THEN
            CALL POLY3 (PSI7+PSI8+PSI4, PSI4*(PSI7+PSI8)+ &
            PSI4*PSI4/4.D0, (PSI4*PSI4*(PSI7+PSI8)-A2) &
            /4.D0,PSI20, ISLV)
            IF (ISLV == 0) PSI2 = MIN (MAX(PSI20,ZERO), CHI2)
        ENDIF
    !      PSI2 = 0.5D0*(2.0D0*SQRT(A2/A7)*PSI7 - PSI4)
    !      PSI2 = MIN (PSI2, CHI2)
    !      ENDIF
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (2) = ZERO                             ! NaI
        MOLAL (3) = 2.0D0*PSI2 + PSI4                ! NH4I
        MOLAL (4) = PSI6                             ! CLI
        MOLAL (5) = PSI2+PSI7+PSI8                   ! SO4I
        MOLAL (6) = ZERO                             ! HSO4
        MOLAL (7) = PSI5                             ! NO3I
        MOLAL (8) = ZERO                             ! CAI
        MOLAL (9) = 2.0D0*PSI7                       ! KI
        MOLAL (10)= PSI8                             ! MGI
    
    ! C      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        -MOLAL(9)-2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
    !      CNA2SO4  = MAX(CHI1 - PSI1, ZERO)
        CNH42S4  = MAX(CHI2 - PSI2, ZERO)
        CNH4NO3  = ZERO
        CK2SO4   = MAX(CHI7 - PSI7, ZERO)
        CMGSO4   = ZERO
        CCASO4   = CHI9
              
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
    ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES *********************
    
        CALL CALCMR
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP **************************


! 0    IF (CHI4.LE.TINY) THEN
!         FUNCO2A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
!      ELSE
    20 FUNCO2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
!      ENDIF

    RETURN

! *** END OF FUNCTION FUNCO2A ****************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO1
! *** CASE O1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCO1A)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCO1A, CALCO2A

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMO1) THEN
        SCASE = 'O1 ; SUBCASE 1'
        CALL CALCO1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'O1 ; SUBCASE 1'
    ELSE
        SCASE = 'O1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRH2 (RH, DRMO1, DRNH4NO3, CALCO1A, CALCO2A)
        SCASE = 'O1 ; SUBCASE 2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCO1 ******************************************

    END SUBROUTINE CALCO1
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCO1A
! *** CASE O1A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4

!     SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
!     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
!     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
!     THE SOLID PHASE.

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCO1A
    INCLUDE 'isrpia.inc'
    real :: LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2

! *** CALCULATE NON VOLATILE SOLIDS ***********************************

    CCASO4  = MIN (W(6), W(2))                    ! CCASO4
    SO4FR   = MAX(W(2) - CCASO4, ZERO)
    CAFR    = MAX(W(6) - CCASO4, ZERO)
    CK2SO4  = MIN (0.5D0*W(7), SO4FR)             ! CK2S04
    FRK     = MAX(W(7) - 2.D0*CK2SO4, ZERO)
    SO4FR   = MAX(SO4FR - CK2SO4, ZERO)
    CNA2SO4 = MIN (0.5D0*W(1), SO4FR)             ! CNA2SO4
    FRNA    = MAX(W(1) - 2.D0*CNA2SO4, ZERO)
    SO4FR   = MAX(SO4FR - CNA2SO4, ZERO)
    CMGSO4  = MIN (W(8), SO4FR)                   ! CMGSO4
    FRMG    = MAX(W(8) - CMGSO4, ZERO)
    SO4FR   = MAX(SO4FR - CMGSO4, ZERO)

    CNH42S4 = MAX (SO4FR , ZERO)                  ! CNH42S4

! *** CALCULATE VOLATILE SPECIES **************************************

    ALF     = W(3) - 2.0D0*CNH42S4
    BET     = W(5)
    GAM     = W(4)

    RTSQ    = R*TEMP*R*TEMP
    A1      = XK6/RTSQ
    A2      = XK10/RTSQ
    print *, A2

    THETA1  = GAM - BET*(A2/A1)
    THETA2  = A2/A1

! QUADRATIC EQUATION SOLUTION

    BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
    CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
    DD      = BB*BB - 4.0D0*CC
    IF (DD < ZERO) GOTO 100   ! Solve each reaction seperately

! TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID

    SQDD    = SQRT(DD)
    KAPA1   = 0.5D0*(-BB+SQDD)
    KAPA2   = 0.5D0*(-BB-SQDD)
    LAMDA1  = THETA1 + THETA2*KAPA1
    LAMDA2  = THETA1 + THETA2*KAPA2

    IF (KAPA1 >= ZERO .AND. LAMDA1 >= ZERO) THEN
        IF (ALF-KAPA1-LAMDA1 >= ZERO .AND. &
        BET-KAPA1 >= ZERO .AND. GAM-LAMDA1 >= ZERO) THEN
            KAPA = KAPA1
            LAMDA= LAMDA1
            GOTO 200
        ENDIF
    ENDIF

    IF (KAPA2 >= ZERO .AND. LAMDA2 >= ZERO) THEN
        IF (ALF-KAPA2-LAMDA2 >= ZERO .AND. &
        BET-KAPA2 >= ZERO .AND. GAM-LAMDA2 >= ZERO) THEN
            KAPA = KAPA2
            LAMDA= LAMDA2
            GOTO 200
        ENDIF
    ENDIF

! SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA

    100 KAPA  = ZERO
    LAMDA = ZERO
    DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
    DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)

! NH4CL EQUILIBRIUM

    IF (DD1 >= ZERO) THEN
        SQDD1 = SQRT(DD1)
        KAPA1 = 0.5D0*(ALF+BET + SQDD1)
        KAPA2 = 0.5D0*(ALF+BET - SQDD1)
    
        IF (KAPA1 >= ZERO .AND. KAPA1 <= MIN(ALF,BET)) THEN
            KAPA = KAPA1
        ELSE IF (KAPA2 >= ZERO .AND. KAPA2 <= MIN(ALF,BET)) THEN
            KAPA = KAPA2
        ELSE
            KAPA = ZERO
        ENDIF
    ENDIF

! NH4NO3 EQUILIBRIUM

    IF (DD2 >= ZERO) THEN
        SQDD2 = SQRT(DD2)
        LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
        LAMDA2= 0.5D0*(ALF+GAM - SQDD2)
    
        IF (LAMDA1 >= ZERO .AND. LAMDA1 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1
        ELSE IF (LAMDA2 >= ZERO .AND. LAMDA2 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
        ELSE
            LAMDA = ZERO
        ENDIF
    ENDIF

! IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION

    IF (KAPA > ZERO .AND. LAMDA > ZERO) THEN
        IF (BET < LAMDA/THETA1) THEN
            KAPA = ZERO
        ELSE
            LAMDA= ZERO
        ENDIF
    ENDIF

! *** CALCULATE COMPOSITION OF VOLATILE SPECIES ************************

    200 CONTINUE
    CNH4NO3 = LAMDA
    CNH4CL  = KAPA

    GNH3    = MAX(ALF - KAPA - LAMDA, ZERO)
    GHNO3   = MAX(GAM - LAMDA, ZERO)
    GHCL    = MAX(BET - KAPA, ZERO)

    RETURN

! *** END OF SUBROUTINE CALCO1A *****************************************

    END SUBROUTINE CALCO1A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM8
! *** CASE M8

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4, K2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM8
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCM8 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCM8 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCM8 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCM8 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCM8')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCM8 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM8 ******************************************

    END SUBROUTINE CALCM8




!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCM8
! *** CASE M8

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4, K2SO4

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCM8 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
    !      A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
    !      A7  = XK8 *(WATER/GAMA(1))**2.0
    !      A8  = XK9 *(WATER/GAMA(3))**2.0
    !      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
        MOLAL (8) = PSI11                                 ! CAI
        MOLAL (9) = 2.D0*PSI9                             ! KI
        MOLAL (10)= PSI10                                 ! MGI
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CNA2SO4   = ZERO
        CK2SO4    = ZERO
        CMGSO4    = ZERO
        CCASO4    = CHI11
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCM8 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCM8 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCM8 *******************************************

    END


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM7
! *** CASE M7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM7
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCM7 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCM7 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCM7 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCM7 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCM7')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCM7 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM7 ******************************************

    END SUBROUTINE CALCM7


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCM7
! *** CASE M7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCM7 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
    !      A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
    !      A7  = XK8 *(WATER/GAMA(1))**2.0
    !      A8  = XK9 *(WATER/GAMA(3))**2.0
    !      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            CALL POLY3 (PSI1+PSI10,ZERO,-A9/4.D0, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MAX (MIN (PSI9,CHI9), ZERO)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
        MOLAL (8) = PSI11                                 ! CAI
        MOLAL (9) = 2.D0*PSI9                             ! KI
        MOLAL (10)= PSI10                                 ! MGI
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CNA2SO4   = ZERO
        CK2SO4    = MAX(CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCM7 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCM7 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCM7 *******************************************

    END
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM6
! *** CASE M6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM6
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCM6 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCM6 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCM6 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCM6 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCM6')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCM6 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM6 ******************************************

    END SUBROUTINE CALCM6

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCM6
! *** CASE M6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCM6 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
    !      A7  = XK8 *(WATER/GAMA(1))**2.0
    !      A8  = XK9 *(WATER/GAMA(3))**2.0
    !      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN   !NA2SO4
            RIZ = SQRT(A9/A1)
            AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.0+RIZ)*(PSI7+PSI8)) &
            /(1.0+RIZ)
            BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.0+RIZ))/(1.0+RIZ)
            CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
            -A1/4.D0)/(1.0+RIZ)
        !      AA  = PSI7+PSI8+PSI9+PSI10
        !      BB  = (PSI7+PSI8)*(PSI9+PSI10)+0.25D0*(PSI7+PSI8)**2.
        !      CC  = ((PSI7+PSI8)**2.*(PSI9+PSI10)-A1)/4.0D0
        
            CALL POLY3 (AA,BB,CC,PSI1,ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1,CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
    !      IF (CHI9.GE.TINY .AND. WATER.GT.TINY) THEN
    !         PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
    !         PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
    !      ELSE
    !         PSI9  = ZERO
    !      ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN   !K2SO4
            CALL POLY3 (PSI1+PSI10,ZERO,-A9/4.D0, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (PSI9,CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
        MOLAL (8) = PSI11                                 ! CAI
        MOLAL (9) = 2.D0*PSI9                             ! KI
        MOLAL (10)= PSI10                                 ! MGI
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
        CK2SO4    = MAX(CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCM6 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCM6 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCM6 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM5
! *** CASE M5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM5
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCM5 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCM5 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCM5 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCM5 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCM5')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCM5 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM5 ******************************************

    END SUBROUTINE CALCM5

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCM5
! *** CASE M5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4
!     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCM5 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
    !      A7  = XK8 *(WATER/GAMA(1))**2.0
    !      A8  = XK9 *(WATER/GAMA(3))**2.0
    !      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN   !NA2SO4
            RIZ = SQRT(A9/A1)
            AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.0+RIZ)*(PSI7+PSI8)) &
            /(1.0+RIZ)
            BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.0+RIZ))/(1.0+RIZ)
            CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
            -A1/4.D0)/(1.0+RIZ)
        !      AA  = PSI7+PSI8+PSI9+PSI10
        !      BB  = (PSI7+PSI8)*(PSI9+PSI10)+0.25D0*(PSI7+PSI8)**2.
        !      CC  = ((PSI7+PSI8)**2.*(PSI9+PSI10)-A1)/4.0D0
        
            CALL POLY3 (AA,BB,CC,PSI1,ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1,CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI9 >= TINY .AND. WATER > TINY) THEN
            PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
            PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
        ELSE
            PSI9  = ZERO
        ENDIF
    
    !      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN   !K2SO4
    !      CALL POLY3 (PSI1+PSI10,ZERO,-A9/4.D0, PSI9, ISLV)
    !        IF (ISLV.EQ.0) THEN
    !            PSI9 = MIN (PSI9,CHI9)
    !        ELSE
    !            PSI9 = ZERO
    !        ENDIF
    !      ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
        MOLAL (8) = PSI11                                 ! CAI
        MOLAL (9) = 2.D0*PSI9                             ! KI
        MOLAL (10)= PSI10                                 ! MGI
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
        CK2SO4    = MAX(CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCM5 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCM5 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCM5 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM4
! *** CASE M4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL
!     4. Completely dissolved: NH4NO3, NANO3, NACL

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM4
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    IF (W(4) <= TINY .AND. W(5) <= TINY) THEN
        SCASE = 'M4 ; SUBCASE 1'
        CALL CALCM1A
        SCASE = 'M4 ; SUBCASE 1'
        RETURN
    ENDIF

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCM4 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCM4 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCM4 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCM4 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCM4')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCM4 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM4 ******************************************

    END SUBROUTINE CALCM4

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCM4
! *** CASE M4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL
!     4. Completely dissolved: NH4NO3, NANO3, NACL

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCM4 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A3  = XK6 /(R*TEMP*R*TEMP)
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
    !      A7  = XK8 *(WATER/GAMA(1))**2.0
    !      A8  = XK9 *(WATER/GAMA(3))**2.0
    !      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,TINY),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN   !NA2SO4
            RIZ = SQRT(A9/A1)
            AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.0+RIZ)*(PSI7+PSI8)) &
            /(1.0+RIZ)
            BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.0+RIZ))/(1.0+RIZ)
            CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
            -A1/4.D0)/(1.0+RIZ)
        !      AA  = PSI7+PSI8+PSI9+PSI10
        !      BB  = (PSI7+PSI8)*(PSI9+PSI10)+0.25D0*(PSI7+PSI8)**2.
        !      CC  = ((PSI7+PSI8)**2.*(PSI9+PSI10)-A1)/4.0D0
        
            CALL POLY3 (AA,BB,CC,PSI1,ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1,CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI9 >= TINY .AND. WATER > TINY) THEN
            PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
            PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
        ELSE
            PSI9  = ZERO
        ENDIF
    
    !      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN   !K2SO4
    !      CALL POLY3 (PSI1+PSI10,ZERO,-A9/4.D0, PSI9, ISLV)
    !        IF (ISLV.EQ.0) THEN
    !            PSI9 = MIN (PSI9,CHI9)
    !        ELSE
    !            PSI9 = ZERO
    !        ENDIF
    !      ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
        MOLAL (8) = PSI11                                 ! CAI
        MOLAL (9) = 2.D0*PSI9                             ! KI
        MOLAL (10)= PSI10                                 ! MGI
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
        CK2SO4    = MAX(CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX (MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6), ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCM4 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCM4 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCM4 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM3
! *** CASE M3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL, NACL
!     4. Completely dissolved: NH4NO3, NANO3

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM3
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    IF (W(4) <= TINY) THEN        ! NO3 NOT EXIST, WATER NOT POSSIBLE
        SCASE = 'M3 ; SUBCASE 1'
        CALL CALCM1A
        SCASE = 'M3 ; SUBCASE 1'
        RETURN
    ENDIF

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCM3 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCM3 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCM3 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCM3 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCM3')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCM3 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM3 ******************************************

    END SUBROUTINE CALCM3

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCM3
! *** CASE M3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL, NACL
!     4. Completely dissolved: NH4NO3, NANO3

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCM3 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A3  = XK6 /(R*TEMP*R*TEMP)
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A10 = XK23 *(WATER/GAMA(21))**2.0
    !      A8  = XK9 *(WATER/GAMA(3))**2.0
    !      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,TINY),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
    !      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     ! NACL DISSOLUTION
    !         VITA = 2.D0*PSI1+PSI8+PSI6                 ! AN DE DOULEPSEI KALA VGALE PSI1 APO DW
    !         GKAMA= PSI6*(2.D0*PSI1+PSI8)-A7
    !         DIAK = MAX(VITA**2.0 - 4.0D0*GKAMA,ZERO)
    !         PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
    !         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
    !      ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
            PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    !C
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN   !NA2SO4
            RIZ = SQRT(A9/A1)
            AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.0+RIZ)*(PSI7+PSI8)) &
            /(1.0+RIZ)
            BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.0+RIZ))/(1.0+RIZ)
            CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
            -A1/4.D0)/(1.0+RIZ)
        !      AA  = PSI7+PSI8+PSI9+PSI10
        !      BB  = (PSI7+PSI8)*(PSI9+PSI10)+0.25D0*(PSI7+PSI8)**2.
        !      CC  = ((PSI7+PSI8)**2.*(PSI9+PSI10)-A1)/4.0D0
        
            CALL POLY3 (AA,BB,CC,PSI1,ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1,CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI9 >= TINY) THEN
            PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
            PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
        ELSE
            PSI9  = ZERO
        ENDIF
    
    !      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN   !K2SO4
    !      CALL POLY3 (PSI1+PSI10,ZERO,-A9/4.D0, PSI9, ISLV)
    !        IF (ISLV.EQ.0) THEN
    !            PSI9 = MIN (PSI9,CHI9)
    !        ELSE
    !            PSI9 = ZERO
    !        ENDIF
    !      ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
        MOLAL (8) = PSI11                                 ! CAI
        MOLAL (9) = 2.D0*PSI9                             ! KI
        MOLAL (10)= PSI10                                 ! MGI
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = MAX(CHI7 - PSI7, ZERO)
        CNANO3    = ZERO
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
        CK2SO4    = MAX(CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX (MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6), ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCM3 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCM3 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCM3 *******************************************

    END


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM2
! *** CASE M2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL, NACL, NANO3

!     THERE ARE THREE REGIMES IN THIS CASE:
!     1. NH4NO3(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCH2A)
!     2. NH4NO3(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!     3. NH4NO3(s) NOT POSSIBLE, AND RH >= MDRH. (MDRH REGION)

!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES M1A, M2B
!     RESPECTIVELY (BECAUSE MDRH POINTS COINCIDE).

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM2
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCM1A, CALCM3

! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************

    CALL CALCM1A

    IF (CNH4NO3 > TINY) THEN        ! NO3 EXISTS, WATER POSSIBLE
        SCASE = 'M2 ; SUBCASE 1'
        CALL CALCM2A
        SCASE = 'M2 ; SUBCASE 1'
    ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
        SCASE = 'M2 ; SUBCASE 1'
        CALL CALCM1A
        SCASE = 'M2 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY .AND. RH < DRMM2) THEN      ! DRY AEROSOL
        SCASE = 'M2 ; SUBCASE 2'
    
    ELSEIF (WATER <= TINY .AND. RH >= DRMM2) THEN  ! MDRH OF M2
        SCASE = 'M2 ; SUBCASE 3'
        CALL CALCMDRH2 (RH, DRMM2, DRNANO3, CALCM1A, CALCM3)
        SCASE = 'M2 ; SUBCASE 3'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM2 ******************************************

    END SUBROUTINE CALCM2

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM2A
! *** CASE M2A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL, NACL, NANO3
!     4. Completely dissolved: NH4NO3

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM2A
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCM2A (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCM2A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCM2A (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCM2A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCM2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCM2A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM2A ******************************************

    END SUBROUTINE CALCM2A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCM2A
! *** CASE M2A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL, NACL, NANO3
!     4. Completely dissolved: NH4NO3

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCM2A (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = CHI1
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1  = XK5 *(WATER/GAMA(2))**3.0
        A3  = XK6 /(R*TEMP*R*TEMP)
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A64 = (XK3*XK2/XKW)*(GAMA(10)/GAMA(5)/GAMA(11))**2.0
        A64 = A64*(R*TEMP*WATER)**2.0
    !      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,TINY),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
    !      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     ! NACL DISSOLUTION
    !         VITA = 2.D0*PSI1+PSI8+PSI6
    !         GKAMA= PSI6*(2.D0*PSI1+PSI8)-A7
    !         DIAK = MAX(VITA**2.0 - 4.0D0*GKAMA,ZERO)
    !         PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
    !         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
    !      ENDIF
    !C
    !      IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     ! NANO3 DISSOLUTION
    !         BIT  = 2.D0*PSI1+PSI7+PSI5
    !         GKAM = PSI5*(2.D0*PSI1+PSI8)-A8
    !         DIA  = BIT**2.0 - 4.0D0*GKAM
    !        PSI8 = 0.5D0*( -BIT + SQRT(DIA) )
    !        PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
    !      ENDIF
    !C
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
            PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
        IF (CHI8 > TINY .AND. WATER > TINY) THEN     ! NANO3 DISSOLUTION
            DIAK = (PSI7-PSI5)**2.D0 + 4.D0*A8
            PSI8 = 0.5D0*( -(PSI7+PSI5) + SQRT(DIAK) )
            PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
        ENDIF
    
        IF (CHI1 > TINY .AND. WATER > TINY) THEN   !NA2SO4
            RIZ = SQRT(A9/A1)
            AA  = (0.5D0*RIZ*(PSI7+PSI8)+PSI10+(1.0+RIZ)*(PSI7+PSI8)) &
            /(1.0+RIZ)
            BB  = ((PSI7+PSI8)*(0.5D0*RIZ*(PSI7+PSI8)+PSI10)+0.25D0* &
            (PSI7+PSI8)**2.0*(1.0+RIZ))/(1.0+RIZ)
            CC  = (0.25D0*(PSI7+PSI8)**2.0*(0.5D0*RIZ*(PSI7+PSI8)+PSI10) &
            -A1/4.D0)/(1.0+RIZ)
        
        !      AA  = PSI7+PSI8+PSI9+PSI10
        !      BB  = (PSI7+PSI8)*(PSI9+PSI10)+0.25D0*(PSI7+PSI8)**2.
        !      CC  = ((PSI7+PSI8)**2.*(PSI9+PSI10)-A1)/4.0D0
        !C
            CALL POLY3 (AA,BB,CC,PSI1,ISLV)
            IF (ISLV == 0) THEN
                PSI1 = MIN (PSI1,CHI1)
            ELSE
                PSI1 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI9 >= TINY .AND. WATER > TINY) THEN
        !         PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
            PSI9  = 0.5D0*SQRT(A9/A1)*(2.0D0*PSI1+PSI7+PSI8)
            PSI9  = MAX (MIN (PSI9,CHI9), ZERO)
        ELSE
            PSI9  = ZERO
        ENDIF
    
    !      IF (CHI9.GT.TINY .AND. WATER.GT.TINY) THEN   !K2SO4
    !      CALL POLY3 (PSI1+PSI10,ZERO,-A9/4.D0, PSI9, ISLV)
    !        IF (ISLV.EQ.0) THEN
    !            PSI9 = MAX (MIN (PSI9,CHI9), ZERO)
    !        ELSE
    !            PSI9 = ZERO
    !        ENDIF
    !      ENDIF
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
        MOLAL (3) = PSI4                                  ! NH4I
        MOLAL (4) = PSI6 + PSI7                           ! CLI
        MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
        MOLAL (6) = ZERO                                  ! HSO4I
        MOLAL (7) = PSI5 + PSI8                           ! NO3I
        MOLAL (8) = PSI11                                 ! CAI
        MOLAL (9) = 2.D0*PSI9                             ! KI
        MOLAL (10)= PSI10                                 ! MGI
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH42S4   = ZERO
        CNH4NO3   = ZERO
        CNACL     = MAX(CHI7 - PSI7, ZERO)
        CNANO3    = MAX(CHI8 - PSI8, ZERO)
        CNA2SO4   = MAX(CHI1 - PSI1, ZERO)
        CK2SO4    = MAX(CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6), ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCM2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A64 - ONE
    20 FUNCM2A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCM2A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM1
! *** CASE M1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL, NACL, NANO3, NH4NO3

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCH1A)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCM1A, CALCM2A

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMM1) THEN
        SCASE = 'M1 ; SUBCASE 1'
        CALL CALCM1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'M1 ; SUBCASE 1'
    ELSE
        SCASE = 'M1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRH2 (RH, DRMM1, DRNH4NO3, CALCM1A, CALCM2A)
        SCASE = 'M1 ; SUBCASE 2'
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCM1 ******************************************

    END SUBROUTINE CALCM1

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCM1A
! *** CASE M1A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, NA2SO4, MGSO4, NH4CL, NACL, NANO3, NH4NO3

! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCM1A
    INCLUDE 'isrpia.inc'
    real :: LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR, &
    NO3FR

! *** CALCULATE NON VOLATILE SOLIDS ***********************************

    CCASO4  = MIN (W(6), W(2))                    ! CCASO4
    SO4FR   = MAX(W(2) - CCASO4, ZERO)
    CAFR    = MAX(W(6) - CCASO4, ZERO)
    CK2SO4  = MIN (0.5D0*W(7), SO4FR)             ! CK2S04
    FRK     = MAX(W(7) - 2.D0*CK2SO4, ZERO)
    SO4FR   = MAX(SO4FR - CK2SO4, ZERO)
    CMGSO4  = MIN (W(8), SO4FR)                   ! CMGSO4
    FRMG    = MAX(W(8) - CMGSO4, ZERO)
    SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
    CNA2SO4 = MAX (SO4FR,ZERO)                    ! CNA2SO4
    NAFR    = MAX (W(1)-2.D0*CNA2SO4, ZERO)
    CNANO3  = MIN (NAFR, W(4))                    ! CNANO3
    NO3FR   = MAX (W(4)-CNANO3, ZERO)
    CNACL   = MIN (MAX(NAFR-CNANO3, ZERO), W(5))  ! CNACL
    CLFR    = MAX (W(5)-CNACL, ZERO)

! *** CALCULATE VOLATILE SPECIES **************************************

    ALF     = W(3)                     ! FREE NH3
    BET     = CLFR                     ! FREE CL
    GAM     = NO3FR                    ! FREE NO3

    RTSQ    = R*TEMP*R*TEMP
    A1      = XK6/RTSQ
    A2      = XK10/RTSQ

    THETA1  = GAM - BET*(A2/A1)
    THETA2  = A2/A1

! QUADRATIC EQUATION SOLUTION

    BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
    CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
    DD      = BB*BB - 4.0D0*CC
    IF (DD < ZERO) GOTO 100   ! Solve each reaction seperately

! TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID

    SQDD    = SQRT(DD)
    KAPA1   = 0.5D0*(-BB+SQDD)
    KAPA2   = 0.5D0*(-BB-SQDD)
    LAMDA1  = THETA1 + THETA2*KAPA1
    LAMDA2  = THETA1 + THETA2*KAPA2

    IF (KAPA1 >= ZERO .AND. LAMDA1 >= ZERO) THEN
        IF (ALF-KAPA1-LAMDA1 >= ZERO .AND. &
        BET-KAPA1 >= ZERO .AND. GAM-LAMDA1 >= ZERO) THEN
            KAPA = KAPA1
            LAMDA= LAMDA1
            GOTO 200
        ENDIF
    ENDIF

    IF (KAPA2 >= ZERO .AND. LAMDA2 >= ZERO) THEN
        IF (ALF-KAPA2-LAMDA2 >= ZERO .AND. &
        BET-KAPA2 >= ZERO .AND. GAM-LAMDA2 >= ZERO) THEN
            KAPA = KAPA2
            LAMDA= LAMDA2
            GOTO 200
        ENDIF
    ENDIF

! SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA

    100 KAPA  = ZERO
    LAMDA = ZERO
    DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
    DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)

! NH4CL EQUILIBRIUM

    IF (DD1 >= ZERO) THEN
        SQDD1 = SQRT(DD1)
        KAPA1 = 0.5D0*(ALF+BET + SQDD1)
        KAPA2 = 0.5D0*(ALF+BET - SQDD1)
    
        IF (KAPA1 >= ZERO .AND. KAPA1 <= MIN(ALF,BET)) THEN
            KAPA = KAPA1
        ELSE IF (KAPA2 >= ZERO .AND. KAPA2 <= MIN(ALF,BET)) THEN
            KAPA = KAPA2
        ELSE
            KAPA = ZERO
        ENDIF
    ENDIF

! NH4NO3 EQUILIBRIUM

    IF (DD2 >= ZERO) THEN
        SQDD2 = SQRT(DD2)
        LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
        LAMDA2= 0.5D0*(ALF+GAM - SQDD2)
    
        IF (LAMDA1 >= ZERO .AND. LAMDA1 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1
        ELSE IF (LAMDA2 >= ZERO .AND. LAMDA2 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
        ELSE
            LAMDA = ZERO
        ENDIF
    ENDIF

! IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION

    IF (KAPA > ZERO .AND. LAMDA > ZERO) THEN
        IF (BET < LAMDA/THETA1) THEN
            KAPA = ZERO
        ELSE
            LAMDA= ZERO
        ENDIF
    ENDIF

! *** CALCULATE COMPOSITION OF VOLATILE SPECIES ***********************

    200 CONTINUE
    CNH4NO3 = LAMDA
    CNH4CL  = KAPA

    GNH3    = ALF - KAPA - LAMDA
    GHNO3   = GAM - LAMDA
    GHCL    = BET - KAPA

    RETURN

! *** END OF SUBROUTINE CALCM1A *****************************************

    END SUBROUTINE CALCM1A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP13
! *** CASE P13

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4
!     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP13
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP13 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP13 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP13 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP13 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP13')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP13 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP13 ******************************************

    END SUBROUTINE CALCP13


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP13
! *** CASE P13

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4
!     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP13 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
        A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
    ! *** CALCULATE SPECIATION *********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    
    ! *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
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
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP13 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP13 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP13 *******************************************

    END
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP12
! *** CASE P12

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4
!     4. Completely dissolved: CA(NO3)2, CACL2, KNO3, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP12
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP12 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP12 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP12 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP12 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP12')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP12 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP12 ******************************************

    END SUBROUTINE CALCP12


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP12
! *** CASE P12

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4
!     4. Completely dissolved: CA(NO3)2, CACL2, KNO3, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP12 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI4   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
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

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
        A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN(MAX(PSI5, TINY),CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH4NO3   = ZERO
        CNH4CL    = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = ZERO
        CKCL      = ZERO
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP12 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP12 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP12 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP11
! *** CASE P11

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3
!     4. Completely dissolved: CA(NO3)2, CACL2, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP11
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP11 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP11 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP11 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP11 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP11')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP11 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP11 ******************************************

    END SUBROUTINE CALCP11


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP11
! *** CASE P11

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3
!     4. Completely dissolved: CA(NO3)2, CACL2, KCL, MGSO4,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP11 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = CHI14
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
        A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 =0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH4NO3   = ZERO
        CNH4CL    = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = ZERO
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP11 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP11 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP11 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP10
! *** CASE P10

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4
!     4. Completely dissolved: CA(NO3)2, CACL2, KCL,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP10
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP10 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP10 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP10 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP10 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP10')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP10 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP10 ******************************************

    END SUBROUTINE CALCP10


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP10
! *** CASE P10

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4
!     4. Completely dissolved: CA(NO3)2, CACL2, KCL,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP10 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = CHI14
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
        A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 =0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH4NO3   = ZERO
        CNH4CL    = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = ZERO
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP10 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP10 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP10 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP9
! *** CASE P9

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP9
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP9 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP9 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP9 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP9 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP9')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP9 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP9 ******************************************

    END SUBROUTINE CALCP9


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP9
! *** CASE P9

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP9 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = ZERO
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
        A14 = XK20 *(WATER/GAMA(20))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
        A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI14 > TINY .AND. WATER > TINY) THEN          !KCL
            PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
            PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
            PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH4NO3   = ZERO
        CNH4CL    = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = MAX (CHI14 - PSI14, ZERO)
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP9 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP9 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP9 *******************************************

    END
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP8
! *** CASE P8

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP8
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP8 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP8 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP8 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP8 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP8')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP8 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP8 ******************************************

    END SUBROUTINE CALCP8


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP8
! *** CASE P8

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP8 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = CHI7
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = ZERO
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
        A14 = XK20 *(WATER/GAMA(20))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
        A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI14 > TINY .AND. WATER > TINY) THEN          !KCL
            PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
            PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
            PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
    !      REST  = 2.D0*W(2) + W(4) + W(5)
    !C
    !      DELT1 = 0.0d0
    !      DELT2 = 0.0d0
    !      IF (W(1)+W(6)+W(7)+W(8).GT.REST) THEN
    !C
    !C *** CALCULATE EQUILIBRIUM CONSTANTS **********************************
    !C
    !     ALFA1 = XK26*RH*(WATER/1.0)                   ! CO2(aq) + H2O
    !     ALFA2 = XK27*(WATER/1.0)                      ! HCO3-
    
    !      X     = W(1)+W(6)+W(7)+W(8) - REST            ! EXCESS OF CRUSTALS EQUALS CO2(aq)
    !C
    !      DIAK  = SQRT( (ALFA1)**2.0 + 4.0D0*ALFA1*X)
    !      DELT1 = 0.5*(-ALFA1 + DIAK)
    !      DELT1 = MIN ( MAX (DELT1, ZERO), X)
    !      DELT2 = ALFA2
    !      DELT2 = MIN ( DELT2, DELT1)
    !      MOLAL(1) = DELT1 + DELT2                      ! H+
    !      ELSE
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH4NO3   = ZERO
    !      CNH4CL    = ZERO
        CNACL     = ZERO
        CNANO3    = ZERO
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = MAX (CHI14 - PSI14, ZERO)
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP8 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP8 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP8 *******************************************

    END
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP7
! *** CASE P7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NH4NO3

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP7
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP7 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP7 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP7 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP7 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP7')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP7 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP7 ******************************************

    END SUBROUTINE CALCP7


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP7
! *** CASE P7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NANO3, NH4NO3

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP7 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = CHI8
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = ZERO
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
        A14 = XK20 *(WATER/GAMA(20))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
        A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI14 > TINY .AND. WATER > TINY) THEN          !KCL
            PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
            PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
            PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
            GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
            DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
            PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH4NO3   = ZERO
    !      CNH4CL    = ZERO
        CNACL     = MAX (CHI7 - PSI7, ZERO)
        CNANO3    = ZERO
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = MAX (CHI14 - PSI14, ZERO)
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP7 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP7 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP7 *******************************************

    END
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP6
! *** CASE P6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL, NANO3
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NH4NO3

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP6
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP6 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP6 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP6 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP6 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP6')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP6 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP6 ******************************************

    END SUBROUTINE CALCP6


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP6
! *** CASE P6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL, NANO3
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2, NH4NO3

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP6 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = ZERO
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = ZERO
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
        A14 = XK20 *(WATER/GAMA(20))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) - &
        A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI14 > TINY .AND. WATER > TINY) THEN          !KCL
            PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
            PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
            PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
            GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
            DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
            PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
        IF (CHI8 > TINY .AND. WATER > TINY) THEN     ! NANO3 DISSOLUTION
        !         VIT  = PSI5+PSI13+PSI7+2.D0*PSI12+2.D0*PSI15
        !         GKAM = PSI7*(2.D0*PSI12+PSI5+PSI13+2.D0*PSI15)-A8
        !         DIA  = MAX(VIT*VIT - 4.0D0*GKAM,ZERO)
        !         PSI8 = 0.5D0*( -VIT + SQRT(DIA) )
            PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
            PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
            PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
    !      REST  = 2.D0*W(2) + W(4) + W(5)
    !C
    !      DELT1 = 0.0d0
    !      DELT2 = 0.0d0
    !      IF (W(1)+W(6)+W(7)+W(8).GT.REST) THEN
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
        CNH4NO3   = ZERO
    !      CNH4CL    = ZERO
        CNACL     = MAX (CHI7 - PSI7, ZERO)
        CNANO3    = MAX (CHI8 - PSI8, ZERO)
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = MAX (CHI14 - PSI14, ZERO)
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP6 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP6 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP6 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP5
! *** CASE P5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, KCL, MGSO4,
!                          NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP5
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCP1A, CALCP6

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (W(4) > TINY)   THEN ! NO3 EXIST, WATER POSSIBLE
        SCASE = 'P5 ; SUBCASE 1'
        CALL CALCP5A
        SCASE = 'P5 ; SUBCASE 1'
    ELSE                                      ! NO3, CL NON EXISTANT
        SCASE = 'P1 ; SUBCASE 1'
        CALL CALCP1A
        SCASE = 'P1 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMP5) THEN        ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCP1A
            SCASE = 'P5 ; SUBCASE 2'
            RETURN
        ELSE
            SCASE = 'P5 ; SUBCASE 3'  ! MDRH REGION (CaSO4, K2SO4, KNO3, KCL, MGSO4,
        !                                                    NANO3, NACL, NH4NO3, NH4CL)
            CALL CALCMDRH2 (RH, DRMP5, DRNH4NO3, CALCP1A, CALCP6)
            SCASE = 'P5 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP5 ******************************************

    END SUBROUTINE CALCP5

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP5A
! *** CASE P5A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL, NANO3, NH4NO3
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP5A
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP5 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP5 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP5 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP5 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP5')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP5 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP5A ******************************************

    END SUBROUTINE CALCP5A


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP5
! *** CASE P5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL, NANO3, NH4NO3
!     4. Completely dissolved: CA(NO3)2, CACL2,
!                              MG(NO3)2, MGCL2

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP5 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = ZERO
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = ZERO
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
        A14 = XK20 *(WATER/GAMA(20))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = (CHI5-PSI2)*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) &
        - A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI14 > TINY .AND. WATER > TINY) THEN          !KCL
            PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
            PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
            PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
            GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
            DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
            PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
        IF (CHI8 > TINY .AND. WATER > TINY) THEN     ! NANO3 DISSOLUTION
        !         VIT  = PSI5+PSI13+PSI7+2.D0*PSI12+2.D0*PSI15
        !         GKAM = PSI7*(2.D0*PSI12+PSI5+PSI13+2.D0*PSI15)-A8
        !         DIA  = MAX(VIT*VIT - 4.0D0*GKAM,ZERO)
        !         PSI8 = 0.5D0*( -VIT + SQRT(DIA) )
            PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
            PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
            PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    !C
    !C *** CALCULATE H+ *****************************************************
    !C
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
    !      CNH4NO3   = ZERO
    !      CNH4CL    = ZERO
        CNACL     = MAX (CHI7 - PSI7, ZERO)
        CNANO3    = MAX (CHI8 - PSI8, ZERO)
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = MAX (CHI14 - PSI14, ZERO)
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
    ! *** NH4NO3(s) calculations
    
        A2   = XK10 /(R*TEMP*R*TEMP)
        IF (GNH3*GHNO3 > A2) THEN
            DELT = MIN(GNH3, GHNO3)
            BB = -(GNH3+GHNO3)
            CC = GNH3*GHNO3-A2
            DD = BB*BB - 4.D0*CC
            PSI21 = 0.5D0*(-BB + SQRT(DD))
            PSI22 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI21 > ZERO .AND. PSI21 > ZERO) THEN
                PSI2 = PSI21
            ELSEIF (DELT-PSI22 > ZERO .AND. PSI22 > ZERO) THEN
                PSI2 = PSI22
            ELSE
                PSI2 = ZERO
            ENDIF
        ELSE
            PSI2 = ZERO
        ENDIF
        PSI2 = MAX(MIN(MIN(PSI2,CHI4-PSI4-PSI3),CHI5-PSI5), ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI2, TINY)
        GHCL    = MAX(GHNO3 - PSI2, TINY)
        CNH4NO3 = PSI2
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP5 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP5 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP5 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP4
! *** CASE P4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP4
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCP1A, CALCP5A

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (W(4) > TINY)   THEN ! NO3 EXIST, WATER POSSIBLE
        SCASE = 'P4 ; SUBCASE 1'
        CALL CALCP4A
        SCASE = 'P4 ; SUBCASE 1'
    ELSE                                      ! NO3, CL NON EXISTANT
        SCASE = 'P1 ; SUBCASE 1'
        CALL CALCP1A
        SCASE = 'P1 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMP4) THEN        ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCP1A
            SCASE = 'P4 ; SUBCASE 2'
            RETURN
        ELSE
            SCASE = 'P4 ; SUBCASE 3'  ! MDRH REGION (CaSO4, K2SO4, KNO3, KCL, MGSO4,
        !                                                    MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL)
            CALL CALCMDRH2 (RH, DRMP4, DRMGNO32, CALCP1A, CALCP5A)
            SCASE = 'P4 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP4 ******************************************

    END SUBROUTINE CALCP4

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP4A
! *** CASE P4A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL, NANO3, NH4NO3, MG(NO3)2
!     4. Completely dissolved: CA(NO3)2, CACL2, MGCL2

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP4A
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP4 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP4 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP4 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP4 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP4')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP4 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP4A ******************************************

    END SUBROUTINE CALCP4A


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP4
! *** CASE P4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL, NANO3, NH4NO3, MG(NO3)2
!     4. Completely dissolved: CA(NO3)2, CACL2, MGCL2

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP4 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = ZERO
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = ZERO
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
        A14 = XK20 *(WATER/GAMA(20))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = (CHI5-PSI2)*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) &
        - A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI14 > TINY .AND. WATER > TINY) THEN          !KCL
            PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
            PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
            PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
            GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
            DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
            PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
        IF (CHI8 > TINY .AND. WATER > TINY) THEN     ! NANO3 DISSOLUTION
        !         VIT  = PSI5+PSI13+PSI7+2.D0*PSI12+2.D0*PSI15
        !         GKAM = PSI7*(2.D0*PSI12+PSI5+PSI13+2.D0*PSI15)-A8
        !         DIA  = MAX(VIT*VIT - 4.0D0*GKAM,ZERO)
        !         PSI8 = 0.5D0*( -VIT + SQRT(DIA) )
            PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
            PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
            PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
    !      CNH4CL    = ZERO
    !      CNH4NO3   = ZERO
        CNACL     = MAX (CHI7 - PSI7, ZERO)
        CNANO3    = MAX (CHI8 - PSI8, ZERO)
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = MAX (CHI14 - PSI14, ZERO)
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
    ! *** NH4NO3(s) calculations
    
        A2   = XK10 /(R*TEMP*R*TEMP)
        IF (GNH3*GHNO3 > A2) THEN
            DELT = MIN(GNH3, GHNO3)
            BB = -(GNH3+GHNO3)
            CC = GNH3*GHNO3-A2
            DD = BB*BB - 4.D0*CC
            PSI21 = 0.5D0*(-BB + SQRT(DD))
            PSI22 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI21 > ZERO .AND. PSI21 > ZERO) THEN
                PSI2 = PSI21
            ELSEIF (DELT-PSI22 > ZERO .AND. PSI22 > ZERO) THEN
                PSI2 = PSI22
            ELSE
                PSI2 = ZERO
            ENDIF
        ELSE
            PSI2 = ZERO
        ENDIF
        PSI2 = MAX(MIN(MIN(PSI2,CHI4-PSI4-PSI3),CHI5-PSI5), ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI2, TINY)
        GHCL    = MAX(GHNO3 - PSI2, TINY)
        CNH4NO3 = PSI2
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP4 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP4 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP4 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP3
! *** CASE P3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP3
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCP1A, CALCP4A

! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************

    IF (W(4) > TINY .AND. W(5) > TINY) THEN ! NO3,CL EXIST, WATER POSSIBLE
        SCASE = 'P3 ; SUBCASE 1'
        CALL CALCP3A
        SCASE = 'P3 ; SUBCASE 1'
    ELSE                                      ! NO3, CL NON EXISTANT
        SCASE = 'P1 ; SUBCASE 1'
        CALL CALCP1A
        SCASE = 'P1 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMP3) THEN        ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCP1A
            SCASE = 'P3 ; SUBCASE 2'
            RETURN
        ELSE
            SCASE = 'P3 ; SUBCASE 3'  ! MDRH REGION (CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
        !                                                    MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL)
            CALL CALCMDRH2 (RH, DRMP3, DRCANO32, CALCP1A, CALCP4A)
            SCASE = 'P3 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP3 ******************************************

    END SUBROUTINE CALCP3

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP3A
! *** CASE P3A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL,
!                          NANO3, NH4NO3, MG(NO3)2, CA(NO3)2
!     4. Completely dissolved: CACL2, MGCL2

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP3A
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP3 (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP3 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP3 (PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP3 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP3')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP3 (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP3A ******************************************

    END SUBROUTINE CALCP3A


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP3
! *** CASE P3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL,
!                          NANO3, NH4NO3, MG(NO3)2, CA(NO3)2
!     4. Completely dissolved: CACL2, MGCL2

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP3 (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = ZERO
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = ZERO
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
        A14 = XK20 *(WATER/GAMA(20))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = (CHI5-PSI2)*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) &
        - A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI14 > TINY .AND. WATER > TINY) THEN          !KCL
            PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
            PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
            PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
            GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
            DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
            PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
        IF (CHI8 > TINY .AND. WATER > TINY) THEN     ! NANO3 DISSOLUTION
        !         VIT  = PSI5+PSI13+PSI7+2.D0*PSI12+2.D0*PSI15
        !         GKAM = PSI7*(2.D0*PSI12+PSI5+PSI13+2.D0*PSI15)-A8
        !         DIA  = MAX(VIT*VIT - 4.0D0*GKAM,ZERO)
        !         PSI8 = 0.5D0*( -VIT + SQRT(DIA) )
            PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
            PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
            PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
    !      CNH4CL    = ZERO
    !      CNH4NO3   = ZERO
        CNACL     = MAX (CHI7 - PSI7, ZERO)
        CNANO3    = MAX (CHI8 - PSI8, ZERO)
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = MAX (CHI14 - PSI14, ZERO)
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6), ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
    ! *** NH4NO3(s) calculations
    
        A2   = XK10 /(R*TEMP*R*TEMP)
        IF (GNH3*GHNO3 > A2) THEN
            DELT = MIN(GNH3, GHNO3)
            BB = -(GNH3+GHNO3)
            CC = GNH3*GHNO3-A2
            DD = BB*BB - 4.D0*CC
            PSI21 = 0.5D0*(-BB + SQRT(DD))
            PSI22 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI21 > ZERO .AND. PSI21 > ZERO) THEN
                PSI2 = PSI21
            ELSEIF (DELT-PSI22 > ZERO .AND. PSI22 > ZERO) THEN
                PSI2 = PSI22
            ELSE
                PSI2 = ZERO
            ENDIF
        ELSE
            PSI2 = ZERO
        ENDIF
        PSI2 = MAX(MIN(MIN(PSI2,CHI4-PSI4-PSI3),CHI5-PSI5),ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI2, TINY)
        GHCL    = MAX(GHNO3 - PSI2, TINY)
        CNH4NO3 = PSI2
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP3 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP3 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP3 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP2
! *** CASE P2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

!     THERE ARE THREE REGIMES IN THIS CASE:
!     1. CACL2(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCL2A)
!     2. CACL2(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!     3. CACL2(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL

!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES P1A, P2B
!     RESPECTIVELY
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================


    SUBROUTINE CALCP2
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCP1A, CALCP3A, CALCP4A, CALCP5A, CALCP6

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCP1A

! *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************

    IF (CCACL2 > TINY) THEN
        SCASE = 'P2 ; SUBCASE 1'
        CALL CALCP2A
        SCASE = 'P2 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRMP2) THEN             ! ONLY SOLIDS
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCP1A
            SCASE = 'P2 ; SUBCASE 2'
        ELSE
            IF (CMGCL2 > TINY) THEN
                SCASE = 'P2 ; SUBCASE 3'    ! MDRH (CaSO4, CA(NO3)2, K2SO4, KNO3, KCL, MGSO4, MGCL2,
            !                                                  MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL)
                CALL CALCMDRH2 (RH, DRMP2, DRMGCL2, CALCP1A, CALCP3A)
                SCASE = 'P2 ; SUBCASE 3'
            ENDIF
            IF (WATER <= TINY .AND. RH >= DRMP3 .AND. RH < DRMP4) THEN
                SCASE = 'P2 ; SUBCASE 4'    ! MDRH (CaSO4, K2SO4, KNO3, KCL, MGSO4, CANO32,
            !                                                  MG(NO3)2, NANO3, NACL, NH4NO3, NH4CL)
                CALL CALCMDRH2 (RH, DRMP3, DRCANO32, CALCP1A, CALCP4A)
                SCASE = 'P2 ; SUBCASE 4'
            ENDIF
            IF (WATER <= TINY .AND. RH >= DRMP4 .AND. RH < DRMP5) THEN
                SCASE = 'P2 ; SUBCASE 5'    ! MDRH (CaSO4, K2SO4, KNO3, KCL, MGSO4,
            !                                                  MGNO32, NANO3, NACL, NH4NO3, NH4CL)
                CALL CALCMDRH2 (RH, DRMP4, DRMGNO32, CALCP1A, CALCP5A)
                SCASE = 'P2 ; SUBCASE 5'
            ENDIF
            IF (WATER <= TINY .AND. RH >= DRMP5) THEN
                SCASE = 'P2 ; SUBCASE 6'    ! MDRH (CaSO4, K2SO4, KNO3, KCL, MGSO4,
            !                                                  NANO3, NACL, NH4NO3, NH4CL)
                CALL CALCMDRH2 (RH, DRMP5, DRNH4NO3, CALCP1A, CALCP6)
                SCASE = 'P2 ; SUBCASE 6'
            ELSE
                WATER = TINY
                DO 20 I=1,NIONS
                    MOLAL(I) = ZERO
                20 END DO
                CALL CALCP1A
                SCASE = 'P2 ; SUBCASE 2'
            ENDIF
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP2 ******************************************

    END SUBROUTINE CALCP2

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP2A
! *** CASE P2A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL,
!                          NANO3, NH4NO3, MG(NO3)2, CA(NO3)2
!     4. Completely dissolved: CACL2

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP2A
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

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

    CHI5    = FRNO3                               ! HNO3(g)
    CHI6    = FRCL                                ! HCL(g)
    CHI4    = W(3)                                ! NH3(g)

    CHI3    = ZERO                                ! CNH4CL
    CHI1    = ZERO
    CHI2    = ZERO

    PSI6LO = TINY
    PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI6LO
    Y1 = FUNCP2A (X1)
    IF (ABS(Y1) <= EPS .OR. CHI6 <= TINY) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1+DX
        Y2 = FUNCP2A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED

    IF (ABS(Y2) > EPS) Y2 = FUNCP2A(PSI6LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCP2A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCP2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCP2A (X3)

! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************

    50 CONTINUE
    IF (MOLAL(1) > TINY .AND. MOLAL(5) > TINY) THEN
        CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
        MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
        MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
        MOLAL(6) = DELTA                                ! HSO4 EFFECT
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCP2A ******************************************

    END SUBROUTINE CALCP2A


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCP2A
! *** CASE P2A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CaSO4, K2SO4, KNO3, MGSO4, KCL, NH4CL, NACL,
!                          NANO3, NH4NO3, MG(NO3)2, CA(NO3)2, MGCL2
!     4. Completely dissolved: CACL2

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCP2A (X)
    INCLUDE 'isrpia.inc'

    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = X
    PSI1   = ZERO
    PSI2   = ZERO
    PSI3   = ZERO
    PSI7   = ZERO
    PSI8   = ZERO
    PSI9   = ZERO
    PSI10  = CHI10
    PSI11  = ZERO
    PSI12  = CHI12
    PSI13  = ZERO
    PSI14  = ZERO
    PSI15  = CHI15
    PSI16  = CHI16
    PSI17  = CHI17
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
        A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
        A9  = XK17 *(WATER/GAMA(17))**3.0
        A13 = XK19 *(WATER/GAMA(19))**2.0
        A14 = XK20 *(WATER/GAMA(20))**2.0
        A7  = XK8 *(WATER/GAMA(1))**2.0
        A8  = XK9 *(WATER/GAMA(3))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = (CHI5-PSI2)*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) &
        - A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6-PSI3)
        PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7 + PSI14 + &
        &        2.D0*PSI16 + 2.D0*PSI17)
        PSI5 = MIN (MAX (PSI5, TINY) , CHI5)
    
        IF (W(3) > TINY .AND. WATER > TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = MAX(BB*BB-4.d0*CC,ZERO)
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
        ELSE
            PSI4 = TINY
        ENDIF
    
        IF (CHI13 > TINY .AND. WATER > TINY) THEN          !KNO3
            VHTA  = PSI5+PSI8+2.D0*PSI12+2.D0*PSI15+PSI14+2.D0*PSI9
            GKAMA = (PSI5+PSI8+2.D0*PSI12+2.D0*PSI15)*(2.D0*PSI9+PSI14)-A13
            DELTA = MAX(VHTA*VHTA-4.d0*GKAMA,ZERO)
            PSI13 = 0.5d0*(-VHTA + SQRT(DELTA))
            PSI13 = MIN(MAX(PSI13,ZERO),CHI13)
        ENDIF
    
        IF (CHI14 > TINY .AND. WATER > TINY) THEN          !KCL
            PSI14 = A14/A13*(PSI5+PSI8+2.D0*PSI12+PSI13+2.D0*PSI15) - &
            PSI6-PSI7-2.D0*PSI16-2.D0*PSI17
            PSI14 = MIN (MAX (PSI14, ZERO), CHI14)
        ENDIF
    
        IF (CHI9 > TINY .AND. WATER > TINY) THEN          !K2SO4
            BBP = PSI10+PSI13+PSI14
            CCP = (PSI13+PSI14)*(0.25D0*(PSI13+PSI14)+PSI10)
            DDP = 0.25D0*(PSI13+PSI14)**2.0*PSI10-A9/4.D0
            CALL POLY3 (BBP, CCP, DDP, PSI9, ISLV)
            IF (ISLV == 0) THEN
                PSI9 = MIN (MAX(PSI9,ZERO) , CHI9)
            ELSE
                PSI9 = ZERO
            ENDIF
        ENDIF
    
        IF (CHI7 > TINY .AND. WATER > TINY) THEN     ! NACL DISSOLUTION
            VITA = PSI6+PSI14+PSI8+2.D0*PSI16+2.D0*PSI17
            GKAMA= PSI8*(2.D0*PSI16+PSI6+PSI14+2.D0*PSI17)-A7
            DIAK = MAX(VITA*VITA - 4.0D0*GKAMA,ZERO)
            PSI7 = 0.5D0*( -VITA + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
        ENDIF
    
        IF (CHI8 > TINY .AND. WATER > TINY) THEN     ! NANO3 DISSOLUTION
        !         VIT  = PSI5+PSI13+PSI7+2.D0*PSI12+2.D0*PSI15
        !         GKAM = PSI7*(2.D0*PSI12+PSI5+PSI13+2.D0*PSI15)-A8
        !         DIA  = MAX(VIT*VIT - 4.0D0*GKAM,ZERO)
        !         PSI8 = 0.5D0*( -VIT + SQRT(DIA) )
            PSI8 = A8/A7*(PSI6+PSI7+PSI14+2.D0*PSI16+2.D0*PSI17)- &
            PSI5-2.D0*PSI12-PSI13-2.D0*PSI15
            PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
        ENDIF
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL (2) = PSI8 + PSI7                                     ! NAI
        MOLAL (3) = PSI4                                            ! NH4I
        MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
        MOLAL (5) = PSI9 + PSI10                                    ! SO4I
        MOLAL (6) = ZERO                                            ! HSO4I
        MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
        MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
        MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
        MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
    
    ! *** CALCULATE H+ *****************************************************
    
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
    !C
    !C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
    !C
        SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3) &
        - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
        CALL CALCPH (SMIN, HI, OHI)
        MOLAL (1) = HI
    !      ENDIF
    
        GNH3      = MAX(CHI4 - PSI4, TINY)
        GHNO3     = MAX(CHI5 - PSI5, TINY)
        GHCL      = MAX(CHI6 - PSI6, TINY)
    
    !      CNH4CL    = ZERO
    !      CNH4NO3   = ZERO
        CNACL     = MAX (CHI7 - PSI7, ZERO)
        CNANO3    = MAX (CHI8 - PSI8, ZERO)
        CK2SO4    = MAX (CHI9 - PSI9, ZERO)
        CMGSO4    = ZERO
        CCASO4    = CHI11
        CCANO32   = ZERO
        CKNO3     = MAX (CHI13 - PSI13, ZERO)
        CKCL      = MAX (CHI14 - PSI14, ZERO)
        CMGNO32   = ZERO
        CMGCL2    = ZERO
        CCACL2    = ZERO
    
    ! *** NH4Cl(s) calculations
    
        A3   = XK6 /(R*TEMP*R*TEMP)
        IF (GNH3*GHCL > A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31 > ZERO .AND. PSI31 > ZERO) THEN
                PSI3 = PSI31
            ELSEIF (DELT-PSI32 > ZERO .AND. PSI32 > ZERO) THEN
                PSI3 = PSI32
            ELSE
                PSI3 = ZERO
            ENDIF
        ELSE
            PSI3 = ZERO
        ENDIF
        PSI3 = MAX(MIN(MIN(PSI3,CHI4-PSI4),CHI6-PSI6),ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX(GNH3 - PSI3, TINY)
        GHCL    = MAX(GHCL - PSI3, TINY)
        CNH4CL  = PSI3
    
    ! *** NH4NO3(s) calculations
    
        A2   = XK10 /(R*TEMP*R*TEMP)
        IF (GNH3*GHNO3 > A2) THEN
            DELT = MIN(GNH3, GHNO3)
            BB = -(GNH3+GHNO3)
            CC = GNH3*GHNO3-A2
            DD = BB*BB - 4.D0*CC
            PSI21 = 0.5D0*(-BB + SQRT(DD))
            PSI22 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI21 > ZERO .AND. PSI21 > ZERO) THEN
                PSI2 = PSI21
            ELSEIF (DELT-PSI22 > ZERO .AND. PSI22 > ZERO) THEN
                PSI2 = PSI22
            ELSE
                PSI2 = ZERO
            ENDIF
        ELSE
            PSI2 = ZERO
        ENDIF
        PSI2 = MAX(MIN(MIN(PSI2,CHI4-PSI4-PSI3),CHI5-PSI5),ZERO)
    
    ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
    
        GNH3    = MAX (GNH3 - PSI2, TINY)
        GHCL    = MAX (GHNO3 - PSI2, TINY)
        CNH4NO3 = PSI2
    
        CALL CALCMR                                    ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

! 0    FUNCP2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
    20 FUNCP2A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

    RETURN

! *** END OF FUNCTION FUNCP2A *******************************************

    END


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP1
! *** CASE P1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCP1A)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCP1A, CALCP2A

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRMP1) THEN
        SCASE = 'P1 ; SUBCASE 1'
        CALL CALCP1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'P1 ; SUBCASE 1'
    ELSE
        SCASE = 'P1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRH2 (RH, DRMP1, DRCACL2, CALCP1A, CALCP2A)
        SCASE = 'P1 ; SUBCASE 2'
    ENDIF


    RETURN

! *** END OF SUBROUTINE CALCP1 ******************************************

    END SUBROUTINE CALCP1

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCP1A
! *** CASE P1A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
!                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCP1A
    INCLUDE 'isrpia.inc'
    real :: LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR, &
    NO3FR

! *** CALCULATE NON VOLATILE SOLIDS ***********************************

    CCASO4  = MIN (W(2), W(6))                    !SOLID CASO4
    CAFR    = MAX (W(6) - CCASO4, ZERO)
    SO4FR   = MAX (W(2) - CCASO4, ZERO)
    CK2SO4  = MIN (SO4FR, 0.5D0*W(7))             !SOLID K2SO4
    FRK     = MAX (W(7) - 2.D0*CK2SO4, ZERO)
    SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
    CMGSO4  = SO4FR                               !SOLID MGSO4
    FRMG    = MAX (W(8) - CMGSO4, ZERO)
    CNACL   = MIN (W(1), W(5))                    !SOLID NACL
    NAFR    = MAX (W(1) - CNACL, ZERO)
    CLFR    = MAX (W(5) - CNACL, ZERO)
    CCANO32 = MIN (CAFR, 0.5D0*W(4))              !SOLID CA(NO3)2
    CAFR    = MAX (CAFR - CCANO32, ZERO)
    NO3FR   = MAX (W(4) - 2.D0*CCANO32, ZERO)
    CCACL2  = MIN (CAFR, 0.5D0*CLFR)              !SOLID CACL2
    CAFR    = MAX (CAFR - CCACL2, ZERO)
    CLFR    = MAX (CLFR - 2.D0*CCACL2, ZERO)
    CMGNO32 = MIN (FRMG, 0.5D0*NO3FR)             !SOLID MG(NO3)2
    FRMG    = MAX (FRMG - CMGNO32, ZERO)
    NO3FR   = MAX (NO3FR - 2.D0*CMGNO32, ZERO)
    CMGCL2  = MIN (FRMG, 0.5D0*CLFR)              !SOLID MGCL2
    FRMG    = MAX (FRMG - CMGCL2, ZERO)
    CLFR    = MAX (CLFR - 2.D0*CMGCL2, ZERO)
    CNANO3  = MIN (NAFR, NO3FR)                   !SOLID NANO3
    NAFR    = MAX (NAFR - CNANO3, ZERO)
    NO3FR   = MAX (NO3FR - CNANO3, ZERO)
    CKCL    = MIN (FRK, CLFR)                     !SOLID KCL
    FRK     = MAX (FRK - CKCL, ZERO)
    CLFR    = MAX (CLFR - CKCL, ZERO)
    CKNO3   = MIN (FRK, NO3FR)                    !SOLID KNO3
    FRK     = MAX (FRK - CKNO3, ZERO)
    NO3FR   = MAX (NO3FR - CKNO3, ZERO)

! *** CALCULATE VOLATILE SPECIES **************************************

    ALF     = W(3)                     ! FREE NH3
    BET     = CLFR                     ! FREE CL
    GAM     = NO3FR                    ! FREE NO3

    RTSQ    = R*TEMP*R*TEMP
    A1      = XK6/RTSQ
    A2      = XK10/RTSQ

    THETA1  = GAM - BET*(A2/A1)
    THETA2  = A2/A1

! QUADRATIC EQUATION SOLUTION

    BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
    CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
    DD      = BB*BB - 4.0D0*CC
    IF (DD < ZERO) GOTO 100   ! Solve each reaction seperately

! TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID

    SQDD    = SQRT(DD)
    KAPA1   = 0.5D0*(-BB+SQDD)
    KAPA2   = 0.5D0*(-BB-SQDD)
    LAMDA1  = THETA1 + THETA2*KAPA1
    LAMDA2  = THETA1 + THETA2*KAPA2

    IF (KAPA1 >= ZERO .AND. LAMDA1 >= ZERO) THEN
        IF (ALF-KAPA1-LAMDA1 >= ZERO .AND. &
        BET-KAPA1 >= ZERO .AND. GAM-LAMDA1 >= ZERO) THEN
            KAPA = KAPA1
            LAMDA= LAMDA1
            GOTO 200
        ENDIF
    ENDIF

    IF (KAPA2 >= ZERO .AND. LAMDA2 >= ZERO) THEN
        IF (ALF-KAPA2-LAMDA2 >= ZERO .AND. &
        BET-KAPA2 >= ZERO .AND. GAM-LAMDA2 >= ZERO) THEN
            KAPA = KAPA2
            LAMDA= LAMDA2
            GOTO 200
        ENDIF
    ENDIF

! SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA

    100 KAPA  = ZERO
    LAMDA = ZERO
    DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
    DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)

! NH4CL EQUILIBRIUM

    IF (DD1 >= ZERO) THEN
        SQDD1 = SQRT(DD1)
        KAPA1 = 0.5D0*(ALF+BET + SQDD1)
        KAPA2 = 0.5D0*(ALF+BET - SQDD1)
    
        IF (KAPA1 >= ZERO .AND. KAPA1 <= MIN(ALF,BET)) THEN
            KAPA = KAPA1
        ELSE IF (KAPA2 >= ZERO .AND. KAPA2 <= MIN(ALF,BET)) THEN
            KAPA = KAPA2
        ELSE
            KAPA = ZERO
        ENDIF
    ENDIF

! NH4NO3 EQUILIBRIUM

    IF (DD2 >= ZERO) THEN
        SQDD2 = SQRT(DD2)
        LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
        LAMDA2= 0.5D0*(ALF+GAM - SQDD2)
    
        IF (LAMDA1 >= ZERO .AND. LAMDA1 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1
        ELSE IF (LAMDA2 >= ZERO .AND. LAMDA2 <= MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
        ELSE
            LAMDA = ZERO
        ENDIF
    ENDIF

! IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION

    IF (KAPA > ZERO .AND. LAMDA > ZERO) THEN
        IF (BET < LAMDA/THETA1) THEN
            KAPA = ZERO
        ELSE
            LAMDA= ZERO
        ENDIF
    ENDIF

! *** CALCULATE COMPOSITION OF VOLATILE SPECIES ***********************

    200 CONTINUE
    CNH4NO3 = LAMDA
    CNH4CL  = KAPA

    GNH3    = ALF - KAPA - LAMDA
    GHNO3   = GAM - LAMDA
    GHCL    = BET - KAPA

    RETURN

! *** END OF SUBROUTINE CALCP1A *****************************************

    END SUBROUTINE CALCP1A

!======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCL9
! *** CASE L9

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : CASO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4, MGSO4, NA2SO4, K2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL9
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCL1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCL1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4
    CHI6 = CK2SO4
    CHI7 = CMGSO4
    CHI8 = CKHSO4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = CNA2SO4
    PSI5 = CNH42S4
    PSI6 = CK2SO4
    PSI7 = CMGSO4
    PSI8 = CKHSO4

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9              ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = MAX(BB*BB - 4.D0*CC, ZERO)
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = PSI2 + PSI3 + PSI1 + PSI8 - LAMDA                ! HSO4I
        MOLAL(9) = PSI8 + 2.0D0*PSI6                                ! KI
        MOLAL(10)= PSI7                                             ! MGI
    
        CLC      = ZERO
        CNAHSO4  = ZERO
        CNA2SO4  = ZERO
        CNH42S4  = ZERO
        CNH4HS4  = ZERO
        CK2SO4   = ZERO
        CMGSO4   = ZERO
        CKHSO4   = ZERO
    
        CALL CALCMR                                         ! Water content

    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

    20 RETURN

! *** END OF SUBROUTINE CALCL9 *****************************************

    END SUBROUTINE CALCL9
!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE CALCL8
! *** CASE L8

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4, MGSO4, NA2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL8
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCL1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCL1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4
    CHI6 = CK2SO4
    CHI7 = CMGSO4
    CHI8 = CKHSO4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = CNA2SO4
    PSI5 = CNH42S4
    PSI6 = ZERO
    PSI7 = CMGSO4
    PSI8 = CKHSO4

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI6LO = ZERO                ! Low  limit
    PSI6HI = CHI6                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    IF (CHI6 <= TINY) THEN
        Y1 = FUNCL8 (ZERO)
        GOTO 50
    ENDIF

    X1 = PSI6HI
    Y1 = FUNCL8 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH K2SO4 *********

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCL8 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH K2SO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCL8 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCL8')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF
! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL8 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCL8')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL8 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCL8 *****************************************

    END SUBROUTINE CALCL8

!=======================================================================

! *** ISORROPIA CODE II
! *** FUNCTION FUNCL8
! *** CASE L8

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4, MGSO4, NA2SO4

!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL8 (P6)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI6   = P6

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9              ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = BB*BB - 4.D0*CC
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     ! HSO4I
        MOLAL(9) = PSI8 + 2.0*PSI6                                  ! KI
        MOLAL(10)= PSI7                                             ! MGI
    
        CLC      = ZERO
        CNAHSO4  = ZERO
        CNA2SO4  = ZERO
        CNH42S4  = ZERO
        CNH4HS4  = ZERO
        CK2SO4   = MAX(CHI6 - PSI6, ZERO)
        CMGSO4   = ZERO
        CKHSO4   = ZERO
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A6 = XK17*(WATER/GAMA(17))**3.0
    FUNCL8 = MOLAL(9)*MOLAL(9)*MOLAL(5)/A6 - ONE
    RETURN

! *** END OF FUNCTION FUNCL8 ****************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL7
! *** CASE L7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4, MGSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL7
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCL1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCL1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4
    CHI6 = CK2SO4
    CHI7 = CMGSO4
    CHI8 = CKHSO4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = ZERO
    PSI5 = CNH42S4
    PSI6 = ZERO
    PSI7 = CMGSO4
    PSI8 = CKHSO4

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI4LO = ZERO                ! Low  limit
    PSI4HI = CHI4                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    IF (CHI4 <= TINY) THEN
        Y1 = FUNCL7 (ZERO)
        GOTO 50
    ENDIF

    X1 = PSI4HI
    Y1 = FUNCL7 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH K2SO4 *********

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCL7 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH K2SO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCL7 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCL7')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF
! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL7 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCL7')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL7 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCL7 *****************************************

    END SUBROUTINE CALCL7

!=======================================================================

! *** ISORROPIA CODE II
! *** FUNCTION FUNCL7
! *** CASE L7

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, NA2SO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4, MGSO4

!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL7 (P4)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK5 *(WATER/GAMA(2))**3.0
        A6 = XK17*(WATER/GAMA(17))**3.0
        A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
    !      PSI6 = 0.5*(SQRT(A6/A4)*(2.D0*PSI4+PSI3)-PSI8)             ! PSI6
    !      PSI6 = MIN (MAX (PSI6, ZERO), CHI6)
    
        IF (CHI6 > TINY .AND. WATER > TINY) THEN
            AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
            BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
            CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
            CALL POLY3 (AA, BB, CC, PSI6, ISLV)
            IF (ISLV == 0) THEN
                PSI6 = MIN (PSI6, CHI6)
            ELSE
                PSI6 = ZERO
            ENDIF
        ENDIF
    
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9              ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = BB*BB - 4.D0*CC
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
    
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     ! HSO4I
        MOLAL(9) = PSI8 + 2.0*PSI6                                  ! KI
        MOLAL(10)= PSI7                                             ! MGI
    
        CLC      = ZERO
        CNAHSO4  = ZERO
        CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
        CNH42S4  = ZERO
        CNH4HS4  = ZERO
        CK2SO4   = MAX(CHI6 - PSI6, ZERO)
        CMGSO4   = ZERO
        CKHSO4   = ZERO
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4     = XK5 *(WATER/GAMA(2))**3.0
    FUNCL7 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCL7 ****************************************

    END


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL6
! *** CASE L6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, NA2SO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL6
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCL1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCL1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4
    CHI6 = CK2SO4
    CHI7 = CMGSO4
    CHI8 = CKHSO4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = ZERO
    PSI5 = CNH42S4
    PSI6 = ZERO
    PSI7 = ZERO
    PSI8 = CKHSO4

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI4LO = ZERO                ! Low  limit
    PSI4HI = CHI4                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    IF (CHI4 <= TINY) THEN
        Y1 = FUNCL6 (ZERO)
        GOTO 50
    ENDIF

    X1 = PSI4HI
    Y1 = FUNCL6 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH K2SO4 *********

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCL6 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH K2SO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCL6 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCL6')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL6 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCL6')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL6 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCL6 *****************************************

    END SUBROUTINE CALCL6

!=======================================================================

! *** ISORROPIA CODE II
! *** FUNCTION FUNCL6
! *** CASE L6

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, NA2SO4

!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL6 (P4)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK5*(WATER/GAMA(2))**3.0
        A6 = XK17*(WATER/GAMA(17))**3.0
        A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
    !      PSI6 = 0.5*(SQRT(A6/A4)*(2.D0*PSI4+PSI3)-PSI8)             ! PSI6
    !      PSI6 = MIN (MAX (PSI6, ZERO), CHI6)
    
        IF (CHI6 > TINY .AND. WATER > TINY) THEN
            AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
            BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
            CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
            CALL POLY3 (AA, BB, CC, PSI6, ISLV)
            IF (ISLV == 0) THEN
                PSI6 = MIN (PSI6, CHI6)
            ELSE
                PSI6 = ZERO
            ENDIF
        ENDIF
    
        PSI7 = CHI7
    
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = BB*BB - 4.D0*CC
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     ! HSO4I
        MOLAL(9) = PSI8 + 2.0*PSI6                                  ! KI
        MOLAL(10)= PSI7                                             ! MGI
    
        CLC      = ZERO
        CNAHSO4  = ZERO
        CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
        CNH42S4  = ZERO
        CNH4HS4  = ZERO
        CK2SO4   = MAX(CHI6 - PSI6, ZERO)
        CMGSO4   = ZERO
        CKHSO4   = ZERO
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4 = XK5 *(WATER/GAMA(2))**3.0
    FUNCL6 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCL6 ****************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL5
! *** CASE L5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NA2SO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL5
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCL1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCL1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4
    CHI6 = CK2SO4
    CHI7 = CMGSO4
    CHI8 = CKHSO4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = ZERO
    PSI5 = CNH42S4
    PSI6 = ZERO
    PSI7 = ZERO
    PSI8 = ZERO

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI4LO = ZERO                ! Low  limit
    PSI4HI = CHI4                ! High limit


! *** INITIAL VALUES FOR BISECTION ************************************

    IF (CHI4 <= TINY) THEN
        Y1 = FUNCL5 (ZERO)
        GOTO 50
    ENDIF

    X1 = PSI4HI
    Y1 = FUNCL5 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 *********

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************


    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, PSI4LO)
        Y2 = FUNCL5 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NA2SO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCL5 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCL5')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL5 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCL5')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL5 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCL5 *****************************************

    END SUBROUTINE CALCL5

!=======================================================================

! *** ISORROPIA CODE II
! *** FUNCTION FUNCL5
! *** CASE L5

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NA2SO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4

!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL5 (P4)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK5*(WATER/GAMA(2))**3.0
        A6 = XK17*(WATER/GAMA(17))**3.0
        A8 = XK18*(WATER/GAMA(18))**2.0
        A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
    !      PSI6 = 0.5*(SQRT(A6/A4)*(2.D0*PSI4+PSI3)-PSI8)             ! PSI6
    !      PSI6 = MIN (MAX (PSI6, ZERO), CHI6)
    
        IF (CHI6 > TINY .AND. WATER > TINY) THEN
            AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
            BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
            CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
            CALL POLY3 (AA, BB, CC, PSI6, ISLV)
            IF (ISLV == 0) THEN
                PSI6 = MIN (PSI6, CHI6)
            ELSE
                PSI6 = ZERO
            ENDIF
        ENDIF
    
        PSI7 = CHI7
    
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = MAX(BB*BB - 4.D0*CC, ZERO)
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
    
        BITA = PSI3 + PSI2 + PSI1 + 2.D0*PSI6 - LAMDA
        CAMA = 2.D0*PSI6*(PSI3 + PSI2 + PSI1 - LAMDA) - A8
        DELT  = MAX(BITA*BITA - 4.D0*CAMA, ZERO)
        PSI8 = 0.5D0*(-BITA + SQRT(DELT))
        PSI8 = MIN(MAX (PSI8, ZERO), CHI8)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     ! HSO4I
        MOLAL(9) = PSI8 + 2.0D0*PSI6                                ! KI
        MOLAL(10)= PSI7                                             ! MGI
    
        CLC      = ZERO
        CNAHSO4  = ZERO
        CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
        CNH42S4  = ZERO
        CNH4HS4  = ZERO
        CK2SO4   = MAX(CHI6 - PSI6, ZERO)
        CMGSO4   = ZERO
        CKHSO4   = MAX(CHI8 - PSI8, ZERO)
    
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    

        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4     = XK5 *(WATER/GAMA(2))**3.0
    FUNCL5 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE

    RETURN

! *** END OF FUNCTION FUNCL5 ****************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL4
! *** CASE L4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, (NH4)2SO4, NA2SO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL4
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCL1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCL1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4
    CHI6 = CK2SO4
    CHI7 = CMGSO4
    CHI8 = CKHSO4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = CLC
    PSI3 = CNAHSO4
    PSI4 = ZERO
    PSI5 = ZERO
    PSI6 = ZERO
    PSI7 = ZERO
    PSI8 = ZERO

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI4LO = ZERO                ! Low  limit
    PSI4HI = CHI4                ! High limit

    IF (CHI4 <= TINY) THEN
        Y1 = FUNCL4 (ZERO)
        GOTO 50
    ENDIF

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI4HI
    Y1 = FUNCL4 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 *********

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCL4 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NA2SO4 **

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCL4 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCL4')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL4 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCL4')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL4 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCL4 *****************************************

    END SUBROUTINE CALCL4

!=======================================================================

! *** ISORROPIA CODE II
! *** FUNCTION FUNCL4
! *** CASE L4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, (NH4)2SO4, NA2SO4
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC

!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL4 (P4)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK5*(WATER/GAMA(2))**3.0
        A5 = XK7*(WATER/GAMA(4))**3.0
        A6 = XK17*(WATER/GAMA(17))**3.0
        A8 = XK18*(WATER/GAMA(18))**2.0
        A9 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = (PSI3 + 2.D0*PSI4 - SQRT(A4/A5)*(3.D0*PSI2 + PSI1))  & ! psi5
        /2.D0/SQRT(A4/A5)
        PSI5 = MAX (MIN (PSI5, CHI5), ZERO)
    
        PSI7 = CHI7
    
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = MAX(BB*BB - 4.D0*CC, ZERO)
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
    
    !      PSI6 = 0.5*(SQRT(A6/A4)*(2.D0*PSI4+PSI3)-PSI8)             ! PSI6
    !      PSI6 = MIN (MAX (PSI6, ZERO), CHI6)
    
        IF (CHI6 > TINY .AND. WATER > TINY) THEN
            AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
            BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
            CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
            CALL POLY3 (AA, BB, CC, PSI6, ISLV)
            IF (ISLV == 0) THEN
                PSI6 = MIN (PSI6, CHI6)
            ELSE
                PSI6 = ZERO
            ENDIF
        ENDIF
    
        BITA = PSI3 + PSI2 + PSI1 + 2.D0*PSI6 - LAMDA
        CAMA = 2.D0*PSI6*(PSI3 + PSI2 + PSI1 - LAMDA) - A8
        DELT  = MAX(BITA*BITA - 4.D0*CAMA, ZERO)
        PSI8 = 0.5D0*(-BITA + SQRT(DELT))
        PSI8 = MIN(MAX (PSI8, ZERO), CHI8)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     ! HSO4I
        MOLAL(9) = PSI8 + 2.0D0*PSI6                                ! KI
        MOLAL(10)= PSI7                                             ! MGI
    
        CLC      = ZERO
        CNAHSO4  = ZERO
        CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
        CNH42S4  = MAX(CHI5 - PSI5, ZERO)
        CNH4HS4  = ZERO
        CK2SO4   = MAX(CHI6 - PSI6, ZERO)
        CMGSO4   = ZERO
        CKHSO4   = MAX(CHI8 - PSI8, ZERO)
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4     = XK5 *(WATER/GAMA(2))**3.0
    FUNCL4 = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCL4 ****************************************

    END
!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL3
! *** CASE L3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC

!     THERE ARE THREE REGIMES IN THIS CASE:
!     1.(NA,NH4)HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCI3A)
!     2.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!     3.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL

!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
!     RESPECTIVELY

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL3
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCL1A, CALCL4

! *** FIND DRY COMPOSITION *********************************************

    CALL CALCL1A

! *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH *********************

    IF (CNH4HS4 > TINY .OR. CNAHSO4 > TINY) THEN
        SCASE = 'L3 ; SUBCASE 1'
        CALL CALCL3A                     ! FULL SOLUTION
        SCASE = 'L3 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRML3) THEN         ! SOLID SOLUTION
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCL1A
            SCASE = 'L3 ; SUBCASE 2'
        
        ELSEIF (RH >= DRML3) THEN     ! MDRH OF L3
            SCASE = 'L3 ; SUBCASE 3'
            CALL CALCMDRH2 (RH, DRML3, DRLC, CALCL1A, CALCL4)
            SCASE = 'L3 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCL3 *****************************************

    END SUBROUTINE CALCL3

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL3A
! *** CASE L3 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, (NH4)2SO4, NA2SO4, LC
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL3A
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCL1A

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCL1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4
    CHI6 = CK2SO4
    CHI7 = CMGSO4
    CHI8 = CKHSO4

    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = ZERO
    PSI3 = CNAHSO4
    PSI4 = ZERO
    PSI5 = ZERO
    PSI6 = ZERO
    PSI7 = ZERO
    PSI8 = ZERO

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI2LO = ZERO                ! Low  limit
    PSI2HI = CHI2                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI2HI
    Y1 = FUNCL3A (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC *********

    IF (YHI < EPS) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, PSI2LO)
        Y2 = FUNCL3A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC

    IF (Y2 > EPS) Y2 = FUNCL3A (ZERO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL3A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCL3A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL3A (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCL3A *****************************************

    END SUBROUTINE CALCL3A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCL3A
! *** CASE L3 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, (NH4)2SO4, NA2SO4, LC
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL3A (P2)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************


    PSI2   = P2                  ! Save PSI2 in COMMON BLOCK
    PSI4LO = ZERO                ! Low  limit for PSI4
    PSI4HI = CHI4                ! High limit for PSI4

! *** IF NH3 =0, CALL FUNCL3B FOR Y4=0 ********************************

    IF (CHI4 <= TINY) THEN
        FUNCL3A = FUNCL3B (ZERO)
        GOTO 50
    ENDIF

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI4HI
    Y1 = FUNCL3B (X1)
    IF (ABS(Y1) <= EPS) GOTO 50
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 *********

    IF (YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, PSI4LO)
        Y2 = FUNCL3B (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NA2SO4

    IF (Y2 > EPS) Y2 = FUNCL3B (PSI4LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL3B (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0004, 'FUNCL3A')    ! WARNING ERROR: NO CONVERGENCE

! *** INNER LOOP CONVERGED **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL3B (X3)

! *** CALCULATE FUNCTION VALUE FOR INTERNAL LOOP ***************************

    50 A2      = XK13*(WATER/GAMA(13))**5.0
    FUNCL3A = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.0/A2 - ONE
    RETURN

! *** END OF FUNCTION FUNCL3A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** FUNCTION FUNCL3B
! *** CASE L3 ; SUBCASE 2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, (NH4)2SO4, NA2SO4, LC
!     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4

!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL3B (P4)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK5*(WATER/GAMA(2))**3.0
        A5 = XK7*(WATER/GAMA(4))**3.0
        A6 = XK17*(WATER/GAMA(17))**3.0
        A8 = XK18*(WATER/GAMA(18))**2.0
        A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = (PSI3 + 2.D0*PSI4 - SQRT(A4/A5)*(3.D0*PSI2 + PSI1))  & ! psi5
        /2.D0/SQRT(A4/A5)
        PSI5 = MAX (MIN (PSI5, CHI5), ZERO)
    
        PSI7 = CHI7
    
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = MAX(BB*BB - 4.D0*CC, ZERO)
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
    
    !      PSI6 = 0.5*(SQRT(A6/A4)*(2.D0*PSI4+PSI3)-PSI8)             ! PSI6
    !      PSI6 = MIN (MAX (PSI6, ZERO), CHI6)
    
        IF (CHI6 > TINY .AND. WATER > TINY) THEN
            AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
            BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
            CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
            CALL POLY3 (AA, BB, CC, PSI6, ISLV)
            IF (ISLV == 0) THEN
                PSI6 = MIN (PSI6, CHI6)
            ELSE
                PSI6 = ZERO
            ENDIF
        ENDIF
    
        BITA = PSI3 + PSI2 + PSI1 + 2.D0*PSI6 - LAMDA
        CAMA = 2.D0*PSI6*(PSI3 + PSI2 + PSI1 - LAMDA) - A8
        DELT  = MAX(BITA*BITA - 4.D0*CAMA, ZERO)
        PSI8 = 0.5D0*(-BITA + SQRT(DELT))
        PSI8 = MIN(MAX (PSI8, ZERO), CHI8)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     ! HSO4I
        MOLAL(9) = PSI8 + 2.0D0*PSI6                                ! KI
        MOLAL(10)= PSI7                                             ! MGI
    
        CLC      = MAX(CHI2 - PSI2, ZERO)
        CNAHSO4  = ZERO
        CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
        CNH42S4  = MAX(CHI5 - PSI5, ZERO)
        CNH4HS4  = ZERO
        CK2SO4   = MAX(CHI6 - PSI6, ZERO)
        CMGSO4   = MAX(CHI7 - PSI7, ZERO)
        CKHSO4   = MAX(CHI8 - PSI8, ZERO)
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4     = XK5 *(WATER/GAMA(2))**3.0
    FUNCL3B = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCL3B ****************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL2
! *** CASE L2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC

!     THERE ARE THREE REGIMES IN THIS CASE:
!     1. NH4HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCL2A)
!     2. NH4HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!     3. NH4HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL

!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES L1A, L2B
!     RESPECTIVELY

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL2
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCL1A, CALCL3A

! *** FIND DRY COMPOSITION **********************************************

    CALL CALCL1A

! *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************

    IF (CNH4HS4 > TINY) THEN
        SCASE = 'L2 ; SUBCASE 1'
        CALL CALCL2A
        SCASE = 'L2 ; SUBCASE 1'
    ENDIF

    IF (WATER <= TINY) THEN
        IF (RH < DRML2) THEN         ! SOLID SOLUTION ONLY
            WATER = TINY
            DO 10 I=1,NIONS
                MOLAL(I) = ZERO
            10 END DO
            CALL CALCL1A
            SCASE = 'L2 ; SUBCASE 2'
        
        ELSEIF (RH >= DRML2) THEN     ! MDRH OF L2
            SCASE = 'L2 ; SUBCASE 3'
            CALL CALCMDRH2 (RH, DRML2, DRNAHSO4, CALCL1A, CALCL3A)
            SCASE = 'L2 ; SUBCASE 3'
        ENDIF
    ENDIF

    RETURN

! *** END OF SUBROUTINE CALCL2 ******************************************

    END SUBROUTINE CALCL2

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL2A
! *** CASE L2 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
!     4. COMPLETELY DISSOLVED: NH4HSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL2A
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    CHI1 = CNH4HS4               ! Save from CALCL1 run
    CHI2 = CLC
    CHI3 = CNAHSO4
    CHI4 = CNA2SO4
    CHI5 = CNH42S4
    CHI6 = CK2SO4
    CHI7 = CMGSO4
    CHI8 = CKHSO4


    PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
    PSI2 = ZERO
    PSI3 = ZERO
    PSI4 = ZERO
    PSI5 = ZERO
    PSI6 = ZERO
    PSI7 = ZERO
    PSI8 = ZERO

    CALAOU = .TRUE.              ! Outer loop activity calculation flag
    PSI2LO = ZERO                ! Low  limit
    PSI2HI = CHI2                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI2HI
    Y1 = FUNCL2A (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 *********

    IF (YHI < EPS) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, PSI2LO)
        Y2 = FUNCL2A (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NA2SO4

    IF (Y2 > EPS) Y2 = FUNCL2A (ZERO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL2A (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCL2A')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL2A (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCL2A *****************************************

    END SUBROUTINE CALCL2A

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCL2A
! *** CASE L2 ; SUBCASE 1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
!     4. COMPLETELY DISSOLVED: NH4HSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL2A (P2)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************


    PSI2   = P2                  ! Save PSI3 in COMMON BLOCK
    PSI4LO = ZERO                ! Low  limit for PSI4
    PSI4HI = CHI4                ! High limit for PSI4

! *** IF NH3 =0, CALL FUNCL3B FOR Y4=0 ********************************


    IF (CHI4 <= TINY) THEN
        FUNCL2A = FUNCL2B (ZERO)
        GOTO 50
    ENDIF

! *** INITIAL VALUES FOR BISECTION ************************************


    X1 = PSI4HI
    Y1 = FUNCL2B (X1)

    IF (ABS(Y1) <= EPS) GOTO 50
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC *********

    IF (YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = MAX(X1-DX, PSI4LO)
        Y2 = FUNCL2B (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC

    IF (Y2 > EPS) Y2 = FUNCL2B (PSI4LO)
    GOTO 50

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCL2B (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0004, 'FUNCL2A')    ! WARNING ERROR: NO CONVERGENCE

! *** INNER LOOP CONVERGED **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCL2B (X3)

! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************

    50 A2      = XK13*(WATER/GAMA(13))**5.0
    FUNCL2A = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.0/A2 - ONE
    RETURN

! *** END OF FUNCTION FUNCL2A *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCL2B
! *** CASE L2 ; SUBCASE 2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
!     4. COMPLETELY DISSOLVED: NH4HSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCL2B (P4)
    INCLUDE 'isrpia.inc'
    real :: LAMDA
    COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8, &
    CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15, &
    CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, &
    PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13, &
    PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6, &
    A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

! *** SETUP PARAMETERS ************************************************

    PSI4   = P4                  ! Save PSI4 in COMMON BLOCK

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 
    PSI3   = CHI3
    PSI5   = CHI5
    LAMDA  = ZERO
    PSI6   = CHI6
    PSI8   = CHI8

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A3 = XK11*(WATER/GAMA(12))**2.0
        A4 = XK5*(WATER/GAMA(2))**3.0
        A5 = XK7*(WATER/GAMA(4))**3.0
        A6 = XK17*(WATER/GAMA(17))**3.0
        A8 = XK18*(WATER/GAMA(18))**2.0
        A9 = XK1*(WATER)*(GAMA(8)**2.0)/(GAMA(7)**3.0)
    
    !  CALCULATE DISSOCIATION QUANTITIES
    
        PSI5 = (PSI3 + 2.D0*PSI4 - SQRT(A4/A5)*(3.D0*PSI2 + PSI1))  & ! psi5
        /2.D0/SQRT(A4/A5)
        PSI5 = MAX (MIN (PSI5, CHI5), ZERO)
    
        IF (CHI3 > TINY .AND. WATER > TINY) THEN
            AA   = 2.D0*PSI4 + PSI2 + PSI1 + PSI8 - LAMDA
            BB   = 2.D0*PSI4*(PSI2 + PSI1 + PSI8 - LAMDA) - A3
            CC   = ZERO
            CALL POLY3 (AA, BB, CC, PSI3, ISLV)
            IF (ISLV == 0) THEN
                PSI3 = MIN (PSI3, CHI3)
            ELSE
                PSI3 = ZERO
            ENDIF
        ENDIF
    
        PSI7 = CHI7
    
        BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9               ! LAMDA
        CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
        DD   = MAX(BB*BB - 4.D0*CC, ZERO)
        LAMDA= 0.5D0*(-BB + SQRT(DD))
        LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
    
    !      PSI6 = 0.5*(SQRT(A6/A4)*(2.D0*PSI4+PSI3)-PSI8)             ! PSI6
    !      PSI6 = MIN (MAX (PSI6, ZERO), CHI6)
    
        IF (CHI6 > TINY .AND. WATER > TINY) THEN
            AA   = PSI5+PSI4+PSI2+PSI7+PSI8+LAMDA
            BB   = PSI8*(PSI5+PSI4+PSI2+PSI7+0.25D0*PSI8+LAMDA)
            CC   = 0.25D0*(PSI8*PSI8*(PSI5+PSI4+PSI2+PSI7+LAMDA)-A6)
            CALL POLY3 (AA, BB, CC, PSI6, ISLV)
            IF (ISLV == 0) THEN
                PSI6 = MIN (PSI6, CHI6)
            ELSE
                PSI6 = ZERO
            ENDIF
        ENDIF
    
        BITA = PSI3 + PSI2 + PSI1 + 2.D0*PSI6 - LAMDA              ! PSI8
        CAMA = 2.D0*PSI6*(PSI3 + PSI2 + PSI1 - LAMDA) - A8
        DELT  = MAX(BITA*BITA - 4.D0*CAMA, ZERO)
        PSI8 = 0.5D0*(-BITA + SQRT(DELT))
        PSI8 = MIN(MAX (PSI8, ZERO), CHI8)
    
    ! *** CALCULATE SPECIATION ********************************************
    
        MOLAL(1) = LAMDA                                            ! HI
        MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
        MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
        MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
        MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 + PSI8 - LAMDA, TINY)     ! HSO4I
        MOLAL(9) = PSI8 + 2.0D0*PSI6                                ! KI
        MOLAL(10)= PSI7                                             ! MGI
    
        CLC      = MAX(CHI2 - PSI2, ZERO)
        CNAHSO4  = MAX(CHI3 - PSI3, ZERO)
        CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
        CNH42S4  = MAX(CHI5 - PSI5, ZERO)
        CNH4HS4  = ZERO
        CK2SO4   = MAX(CHI6 - PSI6, ZERO)
        CMGSO4   = MAX(CHI7 - PSI7, ZERO)
        CKHSO4   = MAX(CHI8 - PSI8, ZERO)
        CALL CALCMR                                       ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 A4     = XK5 *(WATER/GAMA(2))**3.0
    FUNCL2B = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
    RETURN

! *** END OF FUNCTION FUNCL2B ****************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL1
! *** CASE L1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC

!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCI1A)

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL1
    INCLUDE 'isrpia.inc'
    EXTERNAL CALCL1A, CALCL2A

! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************

    IF (RH < DRML1) THEN
        SCASE = 'L1 ; SUBCASE 1'
        CALL CALCL1A              ! SOLID PHASE ONLY POSSIBLE
        SCASE = 'L1 ; SUBCASE 1'
    ELSE
        SCASE = 'L1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
        CALL CALCMDRH2 (RH, DRML1, DRNH4HS4, CALCL1A, CALCL2A)
        SCASE = 'L1 ; SUBCASE 2'
    ENDIF

! *** AMMONIA IN GAS PHASE **********************************************

!      CALL CALCNH3

    RETURN

! *** END OF SUBROUTINE CALCL1 ******************************************

    END SUBROUTINE CALCL1

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCL1A
! *** CASE L1A

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCL1A
    INCLUDE 'isrpia.inc'

! *** CALCULATE NON VOLATILE SOLIDS ***********************************

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

    CNH4HS4 = ZERO
    CNAHSO4 = ZERO
    CNH42S4 = ZERO
    CKHSO4  = ZERO

    CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
    FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
    FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)

    IF (FRSO4 <= TINY) THEN
        CLC     = MAX(CLC - FRNH4, ZERO)
        CNH42S4 = 2.D0*FRNH4

    ELSEIF (FRNH4 <= TINY) THEN
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
    
        IF (CNA2SO4 > TINY) THEN
            FRSO4  = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
        ENDIF
        IF (CK2SO4 > TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CKHSO4 = 2.D0*FRSO4
            CK2SO4 = MAX(CK2SO4-FRSO4, ZERO)
        ENDIF
    ENDIF

! *** CALCULATE GAS SPECIES ********************************************

    GHNO3 = W(4)
    GHCL  = W(5)
    GNH3  = ZERO

    RETURN

! *** END OF SUBROUTINE CALCL1A *****************************************

    END SUBROUTINE CALCL1A


!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCK4
! *** CASE K4

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : CASO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCK4
    INCLUDE 'isrpia.inc'

    real :: LAMDA, KAPA
    COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
    A1,   A2,   A3,   A4

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE.               ! Outer loop activity calculation flag
    FRST   = .TRUE. 
    CALAIN = .TRUE. 

    CHI1   = W(3)                !  Total NH4 initially as NH4HSO4
    CHI2   = W(1)                !  Total NA initially as NaHSO4
    CHI3   = W(7)                !  Total K initially as KHSO4
    CHI4   = W(8)                !  Total Mg initially as MgSO4

    LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  ! FREE H2SO4
    PSI1   = CHI1                            ! ALL NH4HSO4 DELIQUESCED
    PSI2   = CHI2                            ! ALL NaHSO4 DELIQUESCED
    PSI3   = CHI3                            ! ALL KHSO4 DELIQUESCED
    PSI4   = CHI4                            ! ALL MgSO4 DELIQUESCED

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
    
        BB   = A4+LAMDA+PSI4                               ! KAPA
        CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
        DD   = MAX(BB*BB-4.D0*CC, ZERO)
        KAPA = 0.5D0*(-BB+SQRT(DD))
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (1) = MAX(LAMDA + KAPA, TINY)                         ! HI
        MOLAL (2) = PSI2                                            ! NAI
        MOLAL (3) = PSI1                                            ! NH4I
        MOLAL (5) = MAX(KAPA + PSI4, ZERO)                          ! SO4I
        MOLAL (6) = MAX(LAMDA + PSI1 + PSI2 + PSI3 - KAPA, ZERO)    ! HSO4I
        MOLAL (9) = PSI3                                            ! KI
        MOLAL (10)= PSI4                                            ! MGI
    
        CNH4HS4 = ZERO
        CNAHSO4 = ZERO
        CKHSO4  = ZERO
        CCASO4  = W(6)
        CMGSO4  = ZERO
    
        CALL CALCMR                                      ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

    20 RETURN

! *** END OF SUBROUTINE CALCK4

    END SUBROUTINE CALCK4

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCK3
! *** CASE K3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : KHSO4, CASO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCK3
    INCLUDE 'isrpia.inc'

    real :: LAMDA, KAPA
    COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
    A1,   A2,   A3,   A4

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE.               ! Outer loop activity calculation flag
    CHI1   = W(3)                !  Total NH4 initially as NH4HSO4
    CHI2   = W(1)                !  Total NA initially as NaHSO4
    CHI3   = W(7)                !  Total K initially as KHSO4
    CHI4   = W(8)                !  Total Mg initially as MgSO4

    PSI3LO = TINY                ! Low  limit
    PSI3HI = CHI3                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI3HI
    Y1 = FUNCK3 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH KHSO4 ****

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI3HI-PSI3LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCK3 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH KHSO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCK3 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCK3')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCK3 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCK3')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCK3 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCK3 ******************************************

    END SUBROUTINE CALCK3

!=======================================================================

! *** ISORROPIA CODE
! *** SUBROUTINE FUNCK3
! *** CASE K3

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : KHSO4, CaSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCK3 (P1)
    INCLUDE 'isrpia.inc'
    real :: LAMDA, KAPA
    COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
    A1,   A2,   A3,   A4

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

    LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  ! FREE H2SO4
    PSI3   = P1
    PSI1   = CHI1                             ! ALL NH4HSO4 DELIQUESCED
    PSI2   = CHI2                             ! ALL NaHSO4 DELIQUESCED
    PSI4   = CHI4                             ! ALL MgSO4 DELIQUESCED


! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A3 = XK18 *(WATER/GAMA(18))**2.0
        A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
    
    
        BB   = A4+LAMDA+PSI4                             ! KAPA
        CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
        DD   = MAX(BB*BB-4.D0*CC, ZERO)
        KAPA = 0.5D0*(-BB+SQRT(DD))
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (1) = MAX(LAMDA + KAPA, ZERO)                ! HI
        MOLAL (2) = PSI2                                   ! NAI
        MOLAL (3) = PSI1                                   ! NH4I
        MOLAL (4) = ZERO
        MOLAL (5) = MAX(KAPA + PSI4, ZERO)                 ! SO4I
        MOLAL (6) = MAX(LAMDA+PSI1+PSI2+PSI3-KAPA,ZERO)    ! HSO4I
        MOLAL (7) = ZERO
        MOLAL (8) = ZERO
        MOLAL (9) = PSI3                                   ! KI
        MOLAL (10)= PSI4
    
        CNH4HS4 = ZERO
        CNAHSO4 = ZERO
        CKHSO4  = CHI3-PSI3
        CCASO4  = W(6)
        CMGSO4  = ZERO
    
        CALL CALCMR                                      ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 FUNCK3 = MOLAL(9)*MOLAL(6)/A3 - ONE

! *** END OF FUNCTION FUNCK3 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCK2
! *** CASE K2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NAHSO4, KHSO4, CaSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCK2
    INCLUDE 'isrpia.inc'

    real :: LAMDA, KAPA
    COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
    A1,   A2,   A3,   A4

! *** SETUP PARAMETERS ************************************************

    CALAOU = .TRUE.               ! Outer loop activity calculation flag
    CHI1   = W(3)                !  Total NH4 initially as NH4HSO4
    CHI2   = W(1)                !  Total NA initially as NaHSO4
    CHI3   = W(7)                !  Total K initially as KHSO4
    CHI4   = W(8)                !  Total Mg initially as MgSO4

    PSI3LO = TINY                ! Low  limit
    PSI3HI = CHI3                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI3HI
    Y1 = FUNCK2 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH KHSO4 ****

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI3HI-PSI3LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCK2 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH KHSO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCK2 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN   ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCK2')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCK2 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCK2')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCK2 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCK2 ******************************************

    END SUBROUTINE CALCK2

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCK2
! *** CASE K2

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NAHSO4, KHSO4, CaSO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCK2 (P1)
    INCLUDE 'isrpia.inc'
    real :: LAMDA, KAPA
    COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
    A1,   A2,   A3,   A4

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

    LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  ! FREE H2SO4
    PSI3   = P1
    PSI1   = CHI1                              ! ALL NH4HSO4 DELIQUESCED
    PSI4   = CHI4                              ! ALL MgSO4 DELIQUESCED

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A2 = XK11 *(WATER/GAMA(12))**2.0
        A3 = XK18 *(WATER/GAMA(18))**2.0
        A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
    
        PSI2 = A2/A3*PSI3                                   ! PSI2
        PSI2 = MIN(MAX(PSI2, ZERO),CHI2)
    
        BB   = A4+LAMDA+PSI4                                ! KAPA
        CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
        DD   = MAX(BB*BB-4.D0*CC, ZERO)
        KAPA = 0.5D0*(-BB+SQRT(DD))
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (1) = MAX(LAMDA + KAPA, ZERO)                ! HI
        MOLAL (2) = PSI2                                   ! NAI
        MOLAL (3) = PSI1                                   ! NH4I
        MOLAL (4) = ZERO
        MOLAL (5) = MAX(KAPA + PSI4, ZERO)                 ! SO4I
        MOLAL (6) = MAX(LAMDA+PSI1+PSI2+PSI3-KAPA,ZERO)    ! HSO4I
        MOLAL (7) = ZERO
        MOLAL (8) = ZERO
        MOLAL (9) = PSI3                                   ! KI
        MOLAL (10)= PSI4
    
        CNH4HS4 = ZERO
        CNAHSO4 = CHI2-PSI2
        CKHSO4  = CHI3-PSI3
        CCASO4  = W(6)
        CMGSO4  = ZERO
    
        CALL CALCMR                                      ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT
        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 FUNCK2 = MOLAL(9)*MOLAL(6)/A3 - ONE

! *** END OF FUNCTION FUNCK2 *******************************************

    END

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE CALCK1
! *** CASE K1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, KHSO4, CASO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    SUBROUTINE CALCK1
    INCLUDE 'isrpia.inc'

    real :: LAMDA, KAPA
    COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
    A1,   A2,   A3,   A4

! *** SETUP PARAMETERS ************************************************


    CALAOU = .TRUE.               ! Outer loop activity calculation flag
    CHI1   = W(3)                !  Total NH4 initially as NH4HSO4
    CHI2   = W(1)                !  Total NA initially as NaHSO4
    CHI3   = W(7)                !  Total K initially as KHSO4
    CHI4   = W(8)                !  Total Mg initially as MGSO4

    PSI3LO = TINY                ! Low  limit
    PSI3HI = CHI3                ! High limit

! *** INITIAL VALUES FOR BISECTION ************************************

    X1 = PSI3HI
    Y1 = FUNCK1 (X1)
    YHI= Y1                      ! Save Y-value at HI position

! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH KHSO4 ****

    IF (ABS(Y1) <= EPS .OR. YHI < ZERO) GOTO 50

! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************

    DX = (PSI3HI-PSI3LO)/FLOAT(NDIV)
    DO 10 I=1,NDIV
        X2 = X1-DX
        Y2 = FUNCK1 (X2)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y2) < ZERO) GOTO 20  ! (Y1*Y2 < ZERO)
        X1 = X2
        Y1 = Y2
    10 END DO

! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH KHSO4

    YLO= Y1                      ! Save Y-value at Hi position
    IF (YLO > ZERO .AND. YHI > ZERO) THEN
        Y3 = FUNCK1 (ZERO)
        GOTO 50
    ELSE IF (ABS(Y2) < EPS) THEN       ! X2 IS A SOLUTION
        GOTO 50
    ELSE
        CALL PUSHERR (0001, 'CALCK1')    ! WARNING ERROR: NO SOLUTION
        GOTO 50
    ENDIF

! *** PERFORM BISECTION ***********************************************

    20 DO 30 I=1,MAXIT
        X3 = 0.5*(X1+X2)
        Y3 = FUNCK1 (X3)
        IF (SIGN(1.0,Y1)*SIGN(1.0,Y3) <= ZERO) THEN  ! (Y1*Y3 <= ZERO)
            Y2    = Y3
            X2    = X3
        ELSE
            Y1    = Y3
            X1    = X3
        ENDIF
        IF (ABS(X2-X1) <= EPS*X1) GOTO 40
    30 END DO
    CALL PUSHERR (0002, 'CALCK1')    ! WARNING ERROR: NO CONVERGENCE

! *** CONVERGED ; RETURN **********************************************

    40 X3 = 0.5*(X1+X2)
    Y3 = FUNCK1 (X3)

    50 RETURN

! *** END OF SUBROUTINE CALCK1 ******************************************

    END SUBROUTINE CALCK1

!=======================================================================

! *** ISORROPIA CODE II
! *** SUBROUTINE FUNCK1
! *** CASE K1

!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE super RICH, FREE ACID (SO4RAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, KHSO4, CASO4

! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES

!=======================================================================

    real FUNCTION FUNCK1 (P1)
    INCLUDE 'isrpia.inc'
    real :: LAMDA, KAPA
    COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3, &
    A1,   A2,   A3,   A4

! *** SETUP PARAMETERS ************************************************

    FRST   = .TRUE. 
    CALAIN = .TRUE. 

    LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  ! FREE H2SO4
    PSI3   = P1
    PSI4   = CHI4                                    ! ALL MgSO4 DELIQUESCED

! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************

    DO 10 I=1,NSWEEP
    
        A1 = XK12 *(WATER/GAMA(09))**2.0
        A2 = XK11 *(WATER/GAMA(12))**2.0
        A3 = XK18 *(WATER/GAMA(18))**2.0
        A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
    
        PSI1 = A1/A3*PSI3                                   ! PSI1
        PSI1 = MIN(MAX(PSI1, ZERO),CHI1)
    
        PSI2 = A2/A3*PSI3                                   ! PSI2
        PSI2 = MIN(MAX(PSI2, ZERO),CHI2)
    
        BB   = A4+LAMDA+PSI4                                ! KAPA
        CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
        DD   = MAX(BB*BB-4.D0*CC, ZERO)
        KAPA = 0.5D0*(-BB+SQRT(DD))
    
    ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
    
        MOLAL (1) = MAX(LAMDA + KAPA, ZERO)              ! HI
        MOLAL (2) = PSI2                                 ! NAI
        MOLAL (3) = PSI1                                 ! NH4I
        MOLAL (4) = ZERO                                 ! CLI
        MOLAL (5) = MAX(KAPA + PSI4, ZERO)               ! SO4I
        MOLAL (6) = MAX(LAMDA+PSI1+PSI2+PSI3-KAPA,ZERO)  ! HSO4I
        MOLAL (7) = ZERO                                 ! NO3I
        MOLAL (8) = ZERO                                 ! CAI
        MOLAL (9) = PSI3                                 ! KI
        MOLAL (10)= PSI4                                 ! MGI
    
        CNH4HS4 = CHI1-PSI1
        CNAHSO4 = CHI2-PSI2
        CKHSO4  = CHI3-PSI3
        CCASO4  = W(6)
        CMGSO4  = ZERO
    
        CALL CALCMR                                      ! Water content
    
    ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
    
        IF (FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN) THEN
            CALL CALCACT

        ELSE
            GOTO 20
        ENDIF
    10 END DO

! *** CALCULATE OBJECTIVE FUNCTION ************************************

    20 FUNCK1 = MOLAL(9)*MOLAL(6)/A3 - ONE

! *** END OF FUNCTION FUNCK1 ****************************************

    END

